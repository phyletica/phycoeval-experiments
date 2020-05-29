#!/bin/bash

nspecies=9
ngenomes=2
nchars=50000
prefix="sp"

data_dir="../data"
mkdir -p "$data_dir"

outfile="${data_dir}/species-${nspecies}-genomes-${ngenomes}-chars-${nchars}.yml"
julia --project create_dummy_data.jl \
    --nspecies "$nspecies" \
    --ngenomes "$ngenomes" \
    --ncharacters "$nchars" \
    --prefix "$prefix" \
    > "$outfile"
