#!/bin/bash

set -e

project_dir=".."
bin_dir="${project_dir}/bin"
submission_executable="${bin_dir}/psub"
array_spawner_executable="${bin_dir}/spawn_job_array"
extra_args=()
restrict_nodes=''
dry_run=''
wtime='6:00:00'
expected_nlines=1502
max_njobs=100

usage () {
    echo ""
    echo "Usage:"
    echo "  $0 [ OPTIONS ] <PATH-TO-SIM-DIR-1> [ <PATH-TO-SIM-DIR-2> ... ]"
    echo ""
    echo "Required positional argument:"
    echo "  PATH-TO-SIM-DIR  Path to diretory containing simphycoeval output "
    echo "                   files to be analyzed with phycoeval."
    echo "Optional arguments:"
    echo "  -h|--help        Show help message and exit."
    echo "  -t|--walltime    Max time limit for job."
    echo "                   Default: $wtime."
    echo "  -r|--restrict    Restrict job to lab nodes."
    echo "  --nsub           Submit to queue for \"new\" cluster nodes."
    echo "                   Default is to use general queue."
    echo "  -l|--nlines      Expected number of lines in each state log file"
    echo "                   output by phycoeval. Default: $expected_nlines."
    echo "  -m|--max-njobs   Maximum number of jobs to be run at a time."
    echo "                   Default: $max_njobs."
    echo "  -d|--dry-run     Prepare job submission and echo command, but do"
    echo "                   not submit job array."
    echo ""
}

# process args

if [ "$(echo "$@" | grep -c "=")" -gt 0 ]
then
    echo "ERROR: Do not use '=' for arguments. For example, use"
    echo "'--nlines 1502' instead of '--nlines=1502'."
    exit 1
fi

while [ "$1" != "" ]
do
    case $1 in
        -h| --help)
            usage
            exit
            ;;
        -t| --walltime)
            shift
            wtime="$1"
            shift
            ;;
        -r| --restrict)
            shift
            restrict_nodes=1
            ;;
        --nsub)
            shift
            export myqueue="gen28"
            export myreservation="jro0014_s28"
            ;;
        -l| --nlines)
            shift
            expected_nlines="$1"
            shift
            ;;
        -m| --max-njobs)
            shift
            max_njobs="$1"
            shift
            ;;
        -d| --dry-run)
            shift
            dry_run=1
            ;;
        *)
            extra_args+=("$1")
            shift
            ;;
    esac
done

if [ ! -x "$submission_executable" ]
then
    echo "ERROR: No executable '$submission_executable'"
    exit 1
fi

# Make sure there are positional arguments
if [ ! "${#extra_args[*]}" -gt 0 ]
then
    echo "ERROR: At least one argument must be provided; paths to a sim directories"
    usage
    exit 1
fi

reruns=()
for batch_dir in "${extra_args[@]}"
do

    # Remove any trailing slashes from the path
    batch_dir="$(echo "$batch_dir" | sed 's:/*$::')"
    
    # Make sure the argument is a valid directory
    if [ ! -d "$batch_dir" ]
    then
        echo "ERROR: Path is not a valid directory: $batch_dir"
        usage
        exit 1
    fi
    
    # Make sure the directory is a batch directory
    if [ ! "$(echo "$batch_dir" | grep -c -E "batch-[0-9]{1,12}$")" -gt 0 ]
    then
        echo "ERROR: The path provided doesn't seem to be a sim batch directory:"
        echo "    $batch_dir"
        usage
        exit 1
    fi
    
    full_batch_dir="$(cd "$batch_dir" && pwd)"
    if [ ! -d "$full_batch_dir" ]
    then
        echo "ERROR: Failed to get full path to $batch_dir"
        usage
        exit 1
    fi
    
    psub_flags="-t ${wtime}"
    if [ -n "$restrict_nodes" ]
    then
        psub_flags="-r ${psub_flags}"
    fi
    
    echo "Beginning to vet and consolidate sim analysis files in:"
    echo "  '${batch_dir}'"
    for qsub_path in ${full_batch_dir}/*simphycoeval-sim-*-config-run-*-qsub.sh
    do
        to_run="${qsub_path/-qsub.sh/}"
        run_number="${to_run##*-}"
        qsub_file_name="$(basename "$qsub_path")"
        dir_path="$(dirname "$qsub_path")"
        base_prefix="${qsub_file_name%-run-*}"
        file_prefix="run-${run_number}-${base_prefix}"
        prefix="${dir_path}/${file_prefix}"
        out_file="${prefix}.yml.out"
        state_log="${prefix}-state-run-1.log"
        trees_log="${prefix}-trees-run-1.nex"
        op_log="${prefix}-operator-run-1.log"

        # Consolidate tree logs if run was restarted 
        extra_run_number=2
        while [ -e "${prefix}-trees-run-${extra_run_number}.nex" ]
        do
            mv "${prefix}-trees-run-${extra_run_number}.nex" "$trees_log"
            ((++extra_run_number))
        done
    
        # Consolidate state logs if run was restarted 
        extra_run_number=2
        while [ -e "${prefix}-state-run-${extra_run_number}.log" ]
        do
            mv "${prefix}-state-run-${extra_run_number}.log" "$state_log"
            ((++extra_run_number))
        done
    
        # Consolidate operator logs if run was restarted 
        extra_run_number=2
        while [ -e "${prefix}-operator-run-${extra_run_number}.log" ]
        do
            mv "${prefix}-operator-run-${extra_run_number}.log" "$op_log"
            ((++extra_run_number))
        done
    
        if [ ! -e "$out_file" ]
        then
            echo "No stdout: $qsub_path"
            reruns+=( "$qsub_path" )
            continue
        fi
    
        if [ ! -e "$state_log" ]
        then
            echo "No state log: $qsub_path"
            reruns+=( "$qsub_path" )
            continue
        fi

        if [ ! -e "$trees_log" ]
        then
            echo "No trees log: $qsub_path"
            reruns+=( "$qsub_path" )
            continue
        fi
    
        if [ "$(grep -c "Runtime:" "$out_file")" != 1 ]
        then 
            echo "Incomplete stdout: $qsub_path" 
            reruns+=( "$qsub_path" )
            continue
        fi

        nlines="$(wc -l "$trees_log" | awk '{print $1}')"
        if [ "$nlines" -lt "$expected_nlines" ] || [ "$(tail -n 1 "$trees_log" | grep -c "END;")" != 1 ]
        then
            echo "Incomplete trees log: $qsub_path"
            reruns+=( "$qsub_path" )
            continue
        fi
    
        nlines="$(wc -l "$state_log" | awk '{print $1}')"
        if [ "$nlines" != "$expected_nlines" ]
        then
            echo "Incomplete state log: $qsub_path"
            reruns+=( "$qsub_path" )
            continue
        fi
    
        seed_line="$(grep "seed" "$qsub_path")"
        after_seed="${seed_line##*--seed}"
        expected_seed="$(echo ${after_seed%%--*} | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
        seed_line="$(grep -i "seed" "$out_file")"
        seed="$(echo ${seed_line##*:} | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
        
        if [ "$expected_seed" != "$seed" ]
        then
            echo "Bad seed: $qsub_path"
            reruns+=( "$qsub_path" )
            continue
        fi
    done
done

number_of_reruns="${#reruns[*]}"
if [ "$number_of_reruns" = 0 ]
then
    echo "All analyses are complete and clean!"
    exit 0
fi

# Temp file for storing paths to scripts to run
# /tmp dir of head node was not visible to compute nodes, so putting the temp
# file in working directory
# qsub_path_list="$(mktemp "${TMPDIR:-/tmp/}$(basename "$array_spawner_executable").XXXXXXXXXXXX")"
qsub_path_list="$(mktemp "$(pwd)/$(basename "$array_spawner_executable").XXXXXXXXXXXX")"

for qsub_path in "${reruns[@]}"
do
    echo "$qsub_path" >> "$qsub_path_list"
    dir_path="$(dirname "$qsub_path")"
    file_name="$(basename "$qsub_path")"

    prefix="${file_name/-qsub\.sh/}"
    run_number="${prefix##*run-}"
    sim_base="${prefix%-run-*}"

    op_log_file="${dir_path}/run-${run_number}-${sim_base}-operator-run-1.log"
    state_log_file="${dir_path}/run-${run_number}-${sim_base}-state-run-1.log"
    tree_log_file="${dir_path}/run-${run_number}-${sim_base}-trees-run-1.nex"
    stdout_file="${dir_path}/run-${run_number}-${sim_base}.yml.out"

    if [ -e "$op_log_file" ]
    then
        rm "$op_log_file"
    fi

    if [ -e "$state_log_file" ]
    then
        rm "$state_log_file"
    fi

    if [ -e "$tree_log_file" ]
    then
        rm "$tree_log_file"
    fi

    if [ -e "$stdout_file" ]
    then
        rm "$stdout_file"
    fi
done

echo "Submitting analyses to queue..."
if [ "$max_njobs" -gt "$number_of_reruns" ]
then
    psub_flags="${psub_flags} -a 1-${number_of_reruns}"
else
    psub_flags="${psub_flags} -a 1-${number_of_reruns}%${max_njobs}"
fi
echo $submission_executable $psub_flags "$array_spawner_executable" "$qsub_path_list"
if [ -z "$dry_run" ]
then
    $submission_executable $psub_flags "$array_spawner_executable" "$qsub_path_list"
fi
