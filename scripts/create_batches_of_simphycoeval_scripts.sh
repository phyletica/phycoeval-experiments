#!/bin/bash

batch_seeds=( \
    -b 022835694 \
    -b 066292670 \
    -b 182635256 \
    -b 210523283 \
    -b 247877269 \
    -b 463883070 \
    -b 590705757 \
    -b 645706022 \
    -b 651727163 \
    -b 809536974 \
)

prior_configs=(\
    -p ../configs/species-9-genomes-2-bifurcating-tree-random.yml \
    -p ../configs/species-9-genomes-2-generalized-tree-random.yml \
)
fixed_sim_configs=(\
    -f ../configs/species-9-genomes-2-bifurcating-tree-provided.yml \
    -f ../configs/species-9-genomes-2-generalized-tree-provided.yml \
)
sim_configs=(\
    ../configs/species-9-genomes-2-bifurcating-tree-random.yml \
    ../configs/species-9-genomes-2-generalized-tree-random.yml \
)

julia --project create_batch_of_simphycoeval_scripts.jl -n 10 -t 1000 \
    "${batch_seeds[@]}" \
    "${prior_configs[@]}" \
    "${fixed_sim_configs[@]}" \
    "${sim_configs[@]}"
