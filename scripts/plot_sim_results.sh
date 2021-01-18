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

project_dir=".."
source "${project_dir}/modules-to-load.sh" >/dev/null 2>&1 || echo "    No modules loaded"

julia --project plot_sim_results.jl 1> plot_sim_results.sh.out 2>&1
