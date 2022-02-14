#!/bin/bash

set -e

if [ -n "$PBS_JOBNAME" ]
then
    if [ -f "${PBS_O_HOME}/.bashrc" ]
    then
        source "${PBS_O_HOME}/.bashrc"
    fi
    cd $PBS_O_WORKDIR
fi

project_dir="../.."
output_dir="${project_dir}/gekkonid-output"

source "${project_dir}/modules-to-load.sh" >/dev/null 2>&1 || echo "    No modules loaded"
source "./labels_and_calibrations.sh"

# Unzip state and trees logs
files_to_delete=()
for out_label in ${output_labels[@]}
do
    for pth in ${output_dir}/run-??-threads-*-${out_label}-[st][tr][ae][te][es]-run-1.[ln][oe][gx].gz
    do
        if [ -e "$pth" ]
        then
            file_name="$(basename "$pth")"
            (
                cd "$output_dir" && gzip -k -f -d "$file_name"
            )
            files_to_delete+=("${pth/\.gz/}")
        else
            echo "WARNING: didn't find gzipped file pattern \"$pth\""
        fi
    done
done

idx=0
for out_label in ${output_labels[@]}
do
    # sumchains currently does not support phycoeval state log files
    # pyco-sumchains ${output_dir}/run-??-threads-*-${out_label}-state-run-1.log > "${output_dir}/pyco-sumchains-${out_label}.txt"
    yml_path="${output_dir}/posterior-summary-${out_label}.yml"

    "${project_dir}/bin/sumphycoeval" -f -b 101 --mo "${output_dir}/map-tree-${out_label}.nex" --mco "${output_dir}/map-cladogram-${out_label}.nex" ${output_dir}/run-??-threads-*-${out_label}-trees-run-1.nex > "$yml_path"

    root_height="$(julia --project get_root_height.jl "$yml_path")"

    if [ -n "$root_height" ]
    then
        root_age="${root_calibrations[idx]}"
        multiplier="$(echo "${root_age}/${root_height}" | bc -l)"

        "${project_dir}/bin/sumphycoeval" -f -b 101 -m $multiplier --mo "${output_dir}/scaled-map-tree-${out_label}.nex" --mco "${output_dir}/scaled-map-cladogram-${out_label}.nex" ${output_dir}/run-??-threads-*-${out_label}-trees-run-1.nex > "${output_dir}/scaled-posterior-summary-${out_label}.yml"
    fi
    idx=$(expr $idx + 1)
done

rm ${files_to_delete[@]}
