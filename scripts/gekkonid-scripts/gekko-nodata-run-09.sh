#!/bin/bash

set -e

run=09
nthreads=1
config_file_prefix="gekko-nopoly-long"

if [ -n "$PBS_JOBNAME" ]
then
    if [ "$nthreads" != "$PBS_NUM_PPN" ] 
    then
        echo "ERROR: nthreads and PBS_NUM_PPN do not match!"
        exit 1
    fi
    if [ -f "${PBS_O_HOME}/.bashrc" ]
    then
        source "${PBS_O_HOME}/.bashrc"
    fi
    cd $PBS_O_WORKDIR
fi

project_dir="../.."
exe_path="${project_dir}/bin/phycoeval"

if [ ! -x "$exe_path" ]
then
    echo "ERROR: No executable '${exe_path}'."
    echo "       You probably need to run the project setup script."
    exit 1
fi

source "${project_dir}/modules-to-load.sh" >/dev/null 2>&1 || echo "    No modules loaded"

config_dir="${project_dir}/configs"
output_dir="${project_dir}/gekkonid-output"

if [ ! -e "$output_dir" ]
then
    mkdir "$output_dir"
fi

prefix="${output_dir}/nodata-run-${run}-threads-${nthreads}-"
config_path="${config_dir}/${config_file_prefix}.yml"
out_path="${prefix}${config_file_prefix}.out"

"$exe_path" --ignore-data --seed $run --nthreads "$nthreads" --prefix $prefix --relax-triallelic-sites --relax-missing-sites "$config_path" 1>"$out_path" 2>&1

