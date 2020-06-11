#! /bin/bash

set -e

if [ -n "$PBS_JOBNAME" ]
then
    if [ -f "${PBS_O_HOME}/.bashrc" ]
    then
        source "${PBS_O_HOME}/.bashrc"
    fi
fi

cd /gpfs01/scratch/jro0014/phycoeval-experiments/scripts/simphycoeval-scripts

project_dir="../.."
exe_path="${project_dir}/bin/simphycoeval"

if [ ! -x "$exe_path" ]
then
    echo "ERROR: No executable '${exe_path}'."
    echo "       You probably need to run the project setup script."
    exit 1
fi

source "${project_dir}/modules-to-load.sh" >/dev/null 2>&1 || echo "    No modules loaded"

rng_seed=809536974
number_of_reps=10
sim_name="species-9-genomes-2-bifurcating-tree-provided-fixed"
config_dir="../../configs"
config_path="${config_dir}/species-9-genomes-2-bifurcating-tree-provided.yml"
output_dir="../../simulations/${sim_name}/batch-${rng_seed}"
config_set_up_script_path="${project_dir}/scripts/set_up_configs_for_simulated_data.jl"
qsub_set_up_script_path="${project_dir}/scripts/set_up_phycoeval_qsubs.jl"

mkdir -p "$output_dir"

$exe_path --seed="$rng_seed" -n "$number_of_reps" -t "$number_of_topo_mcmc_gens_per_rep" -o "$output_dir" --fix-model -p "${config_dir}/species-9-genomes-2-bifurcating-tree-random.yml" -p "${config_dir}/species-9-genomes-2-generalized-tree-random.yml" "$config_path" && julia --project "$config_set_up_script_path" "$output_dir" && julia --project "$qsub_set_up_script_path" "$output_dir"
