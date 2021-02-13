#!/bin/bash

# A fossil-calibrated tree estimated from ND2 estimates the age of the root of
# our Cyrtodactylus tree to be 19.852 million years (the posterior mean). The
# posterior mean height of our Cyrtodactylus tree is 0.01173622 subs/site.
# Thus, the rate we will use to rescale our posterior sample of trees is:
#
#   0.01173622 / 19.85211054 = 0.00059 subs/site/my
#
# So, to rescale the trees we will multiply all branch lengths (div times) by:
# 
#   1 / 0.0005911276 = 1691.68 = Cyrt tree multiplier
#
# A calibrated tree estimated by MCMCTree using UCE data estimated the age of
# the root our Gekko tree to be 33.76 million years (the posterior mean). The
# posterior mean height of our Gekko tree is 0.01284022 subs/site.
# So, the rate we will use is:
#
#   0.01284022 / 33.76 = 0.00038 subs/site/my
#
# To rescale the posterior sample of Gekko trees, we will multiply all branch
# lengths by:
#
#   1 / 0.0003802927 = 2629.55 = Gekko tree multiplier
#

cyrt_root_age=19.85211054

gekko_root_age=33.76

export root_calibrations=( \
    "$cyrt_root_age" \
    "$cyrt_root_age" \
    "$gekko_root_age" \
    "$gekko_root_age" \
)

export output_labels=( \
    cyrt-nopoly \
    cyrt-nopoly-bif \
    gekko-nopoly \
    gekko-nopoly-bif \
)
