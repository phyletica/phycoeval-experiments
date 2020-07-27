#!/bin/bash

script_dir="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
project_dir="$(dirname "$script_dir")"

(
cd "$project_dir"
for lfs_file in $(git lfs ls-files -n)
do
    echo "Removing $lfs_file"
    rm -- "$lfs_file"
    echo "Checkout $lfs_file"
    git checkout -- "$lfs_file"
done
)
