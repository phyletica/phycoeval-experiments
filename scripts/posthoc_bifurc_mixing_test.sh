#!/bin/bash

set -e

script_dir="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
project_dir="$(dirname "$script_dir")"

sim_dirs=(\
    species-9-genomes-2-generalized-tree-provided-fixed \
)

archives_needed=(\
    sim-files-configs.tar.gz \
    sim-files-data.tar.gz \
    sim-files-true-parameters.tar.gz \
    sim-files-true-trees.tar.gz \
)

# Make sure we have all of the needed archives tracked by Git LFS pulled down
# locally, and the extract them
for sim_dir in "${sim_dirs[@]}"
do
    for batch_path in "../simulations/${sim_dir}/batch-?????????"
    do
        batch_dir="$(basename "$batch_path")"
        for archive_name in "${archives_needed[@]}"
        do
            archive_project_path="simulations/${sim_dir}/${batch_dir}/${archive_name}"
            # Pull down archive from LFS
            (
                cd "$project_dir"
                git lfs pull --include "$archive_project_path"
            )
            # Extract archive
            (
                cd "$batch_path"
                if [ "$archive_name" = 'sim-files-configs.tar.gz' ] 
                then
                    tar xzf "$archive_name" --wildcards 'simphycoeval-sim-??-species-9-genomes-2-bifurcating-tree-random-config.yml'
                    for cfg_path in simphycoeval-sim-??-species-9-genomes-2-bifurcating-tree-random-config.yml
                    do
                        # Create new config for mixing test
                        new_cfg_path="${cfg_path/tree-random/long-chain}"
                        sed -e "s/starting_tree: random/starting_tree: ..\/..\/..\/starting-trees\/tree-species-9-heights-3-bifurcating.nex/g" -e "s/chain_length: 15000/chain_length: 30000/g" -e "s/sample_frequency: 10/sample_frequency: 20/g" "$cfg_path" > "$new_cfg_path"
                        # Remove original config
                        rm "$cfg_path"
                    done
                else
                    tar xzf "$archive_name"
                fi
            )
        done
    done
done

# Clean up repo by removing local copies of archives tracked by Git LFS
# (
# cd "$project_dir"
# for lfs_file in $(git lfs ls-files -n)
# do
#     echo "Removing $lfs_file"
#     rm -- "$lfs_file"
#     echo "Checkout $lfs_file"
#     git checkout -- "$lfs_file"
# done
# )
