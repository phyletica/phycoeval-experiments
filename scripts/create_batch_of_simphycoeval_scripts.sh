#!/bin/bash

julia --project create_batch_of_simphycoeval_scripts.jl -n 10 -t 1000 \
    -p ../configs/species-9-genomes-2-bifurcating-tree-random.yml \
    -p ../configs/species-9-genomes-2-generalized-tree-random.yml \
    -f ../configs/species-9-genomes-2-bifurcating-tree-provided.yml \
    -f ../configs/species-9-genomes-2-generalized-tree-provided.yml \
    ../configs/species-9-genomes-2-bifurcating-tree-random.yml \
    ../configs/species-9-genomes-2-generalized-tree-random.yml
