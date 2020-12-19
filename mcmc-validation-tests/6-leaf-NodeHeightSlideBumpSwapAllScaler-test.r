#! /usr/bin/env Rscript

# test_general_tree_topo_operators.cpp
# NodeHeightSlideBumpSwapAllScaler< BaseTree<Node> > op;
# op.set_operate_on_root(true);
# unsigned int niterations = 2000000;
# unsigned int sample_freq = 20;

# plot_binom_on_hist <- function(counts, p) {
#     k = seq(from = min(counts), to = max(counts), by = 1)
#     n = sum(counts)
#     binom_probs = dbinom(k, n, p)
#     hist(counts, freq = F)
#     lines(k, binom_probs, type = 'l')
# }

counts = c(125, 102, 96, 115, 99, 118, 88, 106, 117, 111, 110, 100, 98, 110, 106, 96, 115, 98, 105, 108, 110, 97, 87, 109, 96, 88, 95, 93, 102, 106, 125, 112, 124, 114, 107, 108, 93, 111, 127, 107, 99, 103, 103, 120, 93, 118, 114, 107, 116, 96, 117, 106, 102, 105, 103, 121, 110, 113, 102, 79, 87, 83, 107, 109, 139, 109, 105, 104, 96, 117, 98, 104, 130, 115, 121, 92, 109, 111, 98, 118, 87, 109, 112, 105, 119, 112, 129, 104, 102, 108, 101, 124, 98, 99, 100, 91, 98, 118, 91, 109, 101, 105, 114, 115, 98, 91, 106, 116, 100, 94, 115, 99, 96, 111, 113, 117, 105, 99, 106, 97, 116, 108, 109, 101, 101, 99, 123, 101, 97, 85, 102, 113, 88, 115, 104, 104, 108, 117, 110, 99, 110, 85, 125, 118, 112, 122, 96, 93, 114, 100, 114, 110, 98, 93, 114, 109, 119, 102, 109, 98, 116, 101, 124, 102, 111, 112, 102, 103, 82, 110, 97, 129, 107, 110, 112, 97, 98, 116, 89, 100, 94, 112, 87, 101, 99, 106, 105, 111, 110, 105, 103, 116, 100, 108, 101, 104, 92, 112, 98, 92, 112, 111, 118, 93, 116, 117, 104, 100, 107, 107, 116, 107, 129, 94, 98, 102, 107, 110, 118, 108, 102, 123, 103, 93, 103, 101, 113, 98, 95, 114, 86, 97, 99, 99, 105, 105, 90, 102, 113, 91, 106, 102, 97, 101, 127, 112, 86, 98, 120, 85, 87, 94, 129, 102, 95, 117, 103, 115, 107, 110, 110, 83, 93, 107, 114, 126, 94, 115, 102, 122, 104, 108, 116, 125, 98, 95, 96, 123, 94, 110, 104, 100, 105, 109, 83, 106, 103, 117, 94, 101, 105, 105, 117, 111, 114, 104, 103, 121, 101, 121, 93, 110, 102, 101, 104, 107, 100, 117, 106, 110, 111, 119, 86, 129, 103, 120, 122, 117, 102, 112, 109, 123, 101, 108, 96, 90, 110, 107, 101, 87, 117, 119, 111, 98, 88, 102, 119, 109, 102, 117, 98, 108, 83, 104, 116, 100, 107, 123, 111, 105, 116, 115, 96, 105, 98, 118, 106, 84, 101, 104, 120, 92, 107, 113, 105, 122, 124, 94, 105, 101, 105, 95, 104, 106, 116, 110, 102, 91, 107, 112, 116, 124, 113, 99, 123, 111, 103, 113, 105, 108, 88, 98, 100, 120, 122, 125, 121, 109, 107, 105, 112, 110, 100, 94, 95, 109, 104, 116, 93, 115, 96, 83, 106, 142, 99, 94, 114, 106, 118, 103, 111, 98, 97, 96, 113, 112, 98, 111, 111, 109, 122, 87, 106, 123, 109, 114, 124, 124, 112, 108, 111, 90, 82, 131, 92, 113, 123, 107, 105, 89, 120, 116, 85, 101, 81, 99, 104, 82, 103, 89, 107, 129, 116, 113, 98, 92, 91, 86, 97, 113, 105, 95, 86, 108, 102, 141, 85, 109, 92, 110, 87, 117, 114, 105, 105, 92, 125, 101, 111, 95, 105, 124, 104, 108, 103, 102, 120, 104, 96, 110, 115, 125, 93, 108, 95, 114, 128, 108, 113, 96, 108, 119, 103, 103, 106, 119, 107, 113, 113, 114, 89, 103, 87, 89, 100, 119, 102, 101, 121, 113, 112, 107, 106, 120, 107, 114, 109, 105, 100, 109, 96, 108, 109, 98, 106, 104, 109, 109, 101, 101, 116, 114, 109, 96, 96, 97, 93, 89, 97, 121, 132, 105, 93, 106, 122, 91, 114, 97, 108, 114, 111, 117, 100, 99, 109, 95, 101, 102, 81, 114, 109, 119, 100, 97, 100, 89, 104, 109, 102, 110, 88, 99, 113, 115, 110, 96, 128, 102, 105, 89, 107, 113, 98, 115, 103, 114, 110, 97, 120, 107, 118, 110, 100, 101, 92, 125, 97, 109, 108, 115, 116, 121, 133, 104, 117, 105, 123, 95, 111, 108, 116, 84, 93, 93, 107, 98, 94, 108, 103, 91, 122, 122, 112, 100, 112, 105, 125, 116, 114, 97, 104, 106, 120, 115, 110, 100, 110, 101, 98, 110, 119, 95, 100, 98, 97, 97, 82, 110, 106, 120, 102, 118, 98, 111, 107, 86, 108, 96, 97, 83, 103, 102, 108, 115, 104, 98, 103, 108, 112, 106, 98, 100, 93, 109, 101, 103, 79, 95, 81, 117, 90, 129, 96, 106, 111, 105, 103, 103, 108, 93, 103, 90, 105, 115, 93, 108, 108, 117, 96, 89, 117, 93, 112, 101, 103, 90, 132, 77, 109, 98, 110, 108, 116, 93, 109, 103, 82, 91, 104, 87, 95, 118, 131, 105, 111, 118, 107, 103, 106, 101, 95, 122, 116, 97, 104, 95, 118, 105, 100, 114, 97, 122, 107, 90, 128, 120, 94, 109, 93, 104, 121, 98, 127, 109, 110, 108, 93, 94, 89, 110, 84, 101, 116, 87, 109, 87, 97, 123, 90, 91, 103, 104, 121, 117, 106, 96, 116, 103, 106, 115, 92, 102, 96, 111, 108, 119, 109, 107, 107, 101, 103, 107, 100, 101, 107, 100, 98, 112, 105, 106, 95, 95, 106, 104, 98, 107, 114, 106, 107, 96, 111, 103, 91, 111, 110, 111, 99, 108, 103, 108, 117, 109, 106, 118, 88, 110, 113, 113, 82, 112, 128, 112, 119, 115, 136, 102, 104, 117, 105, 113, 126, 127, 113, 124, 104, 108, 113, 109, 104, 100, 118, 125, 106, 120, 96, 111, 111, 115, 99, 106, 116, 111, 107, 100, 106, 103, 107, 111, 106, 108, 99, 95, 126, 116, 97, 120, 109, 103, 95, 94, 103, 94, 103, 115, 80, 95, 99, 98, 96, 106, 117, 111, 87, 99, 103, 113, 107, 117, 104, 108, 102, 110, 101, 108, 93, 100, 99, 93, 100, 104, 111, 100, 112, 130, 105, 103, 117, 101, 99, 124, 115, 104, 121, 101, 112)
number_of_topologies = length(counts)
binomial_prob = 1.0 / number_of_topologies
topology = seq(from = 1, to = number_of_topologies, by = 1)
d = data.frame(topology = topology, count = counts)

# pdf("6-leaf-NodeHeightSlideBumpSwapAllScaler-test-topo-vs-count.pdf")
# plot(x = d$topology, y = d$count, ylim = c(0, max(counts)))
# dev.off()

# pdf("6-leaf-NodeHeightSlideBumpSwapAllScaler-test-count-distribution.pdf")
# plot_binom_on_hist(counts = counts, p = binomial_prob)
# dev.off()

# chisq.test(counts)
