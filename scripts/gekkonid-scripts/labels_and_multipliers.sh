#!/bin/bash

# A fossil-calibrated tree estimated from ND2 estimates the age of the root of
# our Cyrtodactylus tree to be 19.852 million years (the posterior mean). The
# posterior mean height of our Cyrtodactylus tree is 0.01173513 subs/site.
# Thus, the rate we will use to rescale our posterior sample of trees is:
#
#   0.01173513 / 19.85211054 = 0.0005911276 subs/site/my
#
# So, to rescale the trees we will multiply all branch lengths (div times) by:
# 
#   1 / 0.0005911276 = 1691.682 = Cyrt tree multiplier
#
# A calibrated tree estimated by MCMCTree using UCE data estimated the age of
# the root our Gekko tree to be 33.76 million years (the posterior mean). The
# posterior mean height of our Gekko tree is 0.01283976 subs/site.
# So, the rate we will use is:
#
#   0.01283976 / 33.76 = 0.0003802927 subs/site/my
#
# To rescale the posterior sample of Gekko trees, we will multiply all branch
# lengths by:
#
#   1 / 0.0003802927 = 2629.554 = Gekko tree multiplier

cyrt_root_age=19.85211054
cyrt_nopoly_mn_ht=0.01173513
cyrt_nopoly_raw_mult="$(echo "${cyrt_root_age}/${cyrt_nopoly_mn_ht}" | bc -l)"
cyrt_nopoly_multiplier="$(printf %.3f $cyrt_nopoly_raw_mult)"

gekko_root_age=33.76
gekko_nopoly_mn_ht=0.01283976
gekko_nopoly_raw_mult="$(echo "${gekko_root_age}/${gekko_nopoly_mn_ht}" | bc -l)"
gekko_nopoly_multiplier="$(printf %.3f $gekko_nopoly_raw_mult)"

export tree_multipliers=( \
    "$cyrt_nopoly_multiplier" \
    "$gekko_nopoly_multiplier" \
)

export output_labels=( \
    cyrt-nopoly \
    gekko-nopoly \
)
