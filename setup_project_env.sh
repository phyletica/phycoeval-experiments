#! /bin/bash

set -e

ecoevolity_commit="c8a3b592"

# Get path to directory of this script
project_dir="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# Load modules
echo "Loading modules..."
source ./modules-to-load.sh >/dev/null 2>&1 || echo "    No modules loaded"

echo "Setting up julia environment..."
julia -e 'import Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.precompile(); Pkg.status()'

echo "Cloning and building ecoevolity..."
git clone https://github.com/phyletica/ecoevolity.git
(
    cd ecoevolity
    git checkout -b testing "$ecoevolity_commit"
    ./build.sh --prefix "$project_dir"
    echo "    Commit $ecoevolity_commit of ecoevolity successfully built and installed"
)
echo ""
echo "Ecoevolity was successfully installed in './bin'"
echo "You can now remove the ecoevolity directory using the command:"
echo "    rm -rf ecoevolity"
