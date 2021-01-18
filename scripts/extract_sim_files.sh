#!/bin/bash

set -e

usage () {
    echo ""
    echo "Usage:"
    echo "  $0 <PATH-TO-SIM-DIR> [ <PATH-TO-SIM-DIR> ... ]"
    echo ""
    echo "Required positional argument:"
    echo "  PATH-TO-SIM-DIR  Path to diretory containing archives of sim files."
    echo ""
}

if [ -n "$PBS_JOBNAME" ]
then
    if [ -f "${PBS_O_HOME}/.bashrc" ]
    then
        source "${PBS_O_HOME}/.bashrc"
    fi
    cd $PBS_O_WORKDIR
fi

# Make sure there is at least one positional argument
if [ ! "$#" -gt 0 ]
then
    echo "ERROR: At least one argument must be provided; a path to a sim directory"
    usage
    exit 1
fi

# Vet arguments before archiving/removing anything
for batch_dir in $@
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
done

for batch_dir in $@
do
    # Remove any trailing slashes from the path
    batch_dir="$(echo "$batch_dir" | sed 's:/*$::')"
    
    echo "Extracting contents of:"
    echo "  ${batch_dir}"

    (
        cd "$batch_dir"
        for p in *.tar.gz
        do
            tar xzf "$p"
        done
    )
done
