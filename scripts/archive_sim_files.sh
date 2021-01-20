#!/bin/bash

set -e

usage () {
    echo ""
    echo "Usage:"
    echo "  $0 <PATH-TO-SIM-DIR> [ <PATH-TO-SIM-DIR> ... ]"
    echo ""
    echo "Required positional argument:"
    echo "  PATH-TO-SIM-DIR  Path to diretory containing phycoeval output of "
    echo "                   analyses of data sets simulated with simphycoeval."
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
    
    echo "Cleaning up contents of:"
    echo "  ${batch_dir}"

    (
        cd "$batch_dir"
        for p in *-qsub.sh.o*; do if [ -e "$p" ]; then
            echo "Removing PBS output files..."
            rm *-qsub.sh.o*
            break
        fi; done
        for p in *-qsub.sh; do if [ -e "$p" ]; then
            echo "Archiving and removing qsub scripts..."
            tar czf sim-files-qsub-scripts.tar.gz *-qsub.sh && rm *-qsub.sh
            break
        fi; done
        for p in simphycoeval-sim-*-chars-*.yml; do if [ -e "$p" ]; then
            echo "Archiving and removing data files..."
            tar czf sim-files-data.tar.gz simphycoeval-sim-*-chars-*.yml && rm simphycoeval-sim-*-chars-*.yml
            break
        fi; done
        for p in *simphycoeval-*.yml; do if [ -e "$p" ]; then
            echo "Archiving and removing config files..."
            tar czf sim-files-configs.tar.gz *simphycoeval-*.yml && rm *simphycoeval-*.yml
            break
        fi; done
        for p in run-?-*operator*.log; do if [ -e "$p" ]; then
            echo "Archiving and removing phycoeval operator logs..."
            tar czf sim-files-op-logs.tar.gz run-?-*operator*.log && rm run-?-*operator*.log
            break
        fi; done
        for p in run-?-*state*.log; do if [ -e "$p" ]; then
            echo "Archiving and removing phycoeval state logs..."
            tar czf sim-files-state-logs.tar.gz run-?-*state*.log && rm run-?-*state*.log
            break
        fi; done
        for p in run-?-*trees*.nex; do if [ -e "$p" ]; then
            echo "Archiving and removing phycoeval tree logs..."
            tar czf sim-files-tree-logs.tar.gz run-?-*trees*.nex && rm run-?-*trees*.nex
            break
        fi; done
        for p in run-?-*.out; do if [ -e "$p" ]; then
            echo "Archiving and removing phycoeval stdout files..."
            tar czf sim-files-stdout.tar.gz run-?-*.out && rm run-?-*.out
            break
        fi; done
        for p in simphycoeval-sim-*-true-parameters.txt; do if [ -e "$p" ]; then
            echo "Archiving and removing files with true parameters..."
            tar czf sim-files-true-parameters.tar.gz simphycoeval-sim-*-true-parameters.txt && rm simphycoeval-sim-*-true-parameters.txt
            break
        fi; done
        for p in simphycoeval-sim-*-true-tree.phy; do if [ -e "$p" ]; then
            echo "Archiving and removing files with true trees..."
            tar czf sim-files-true-trees.tar.gz simphycoeval-sim-*-true-tree.phy && rm simphycoeval-sim-*-true-tree.phy
            break
        fi; done
        for p in simphycoeval*rejected-trees.phy; do if [ -e "$p" ]; then
            echo "Archiving and removing files with rejected trees..."
            tar czf sim-files-rejected-trees.tar.gz simphycoeval*rejected-trees.phy && rm simphycoeval*rejected-trees.phy
            break
        fi; done
    )
done
