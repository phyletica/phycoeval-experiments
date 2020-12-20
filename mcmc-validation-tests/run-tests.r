#! /usr/bin/env Rscript

plot_width = 5.4
plot_height = 4.6
p_val_digits = 3

source("plot-functions.r")

source("5-leaf-general-tree-test.r")

chi_test = chisq.test(counts)
chi_test_stat = format(chi_test$statistic)
chi_test_p = format(chi_test$p.value, digits = p_val_digits)

pdf("5-leaf-general-tree-test-topo-vs-count.pdf")
plot(x = d$topology, y = d$count, ylim = c(0, max(counts)))
dev.off()

pdf("5-leaf-general-tree-test-count-distribution.pdf",
    width = plot_width, height = plot_height)
# mgp adjusts space of axis labels
par(mgp = c(2.0,0.6,0))
# oma is space around all plots
par(oma = c(0,0,0,0), mar = c(3.1,3.1,1.5,0.5))

plot_binom_on_hist(counts = counts, p = binomial_prob)
mtext(bquote(chi^2 == .(chi_test_stat)),
      side=3, line = -1.5, adj = 0.02)
mtext(bquote(italic(p) == .(chi_test_p)),
      side=3, line = -2.5, adj = 0.02)
dev.off()


source("6-leaf-general-tree-test.r")

chi_test = chisq.test(counts)
chi_test_stat = format(chi_test$statistic)
chi_test_p = format(chi_test$p.value, digits = p_val_digits)

pdf("6-leaf-general-tree-test-topo-vs-count.pdf")
plot(x = d$topology, y = d$count, ylim = c(0, max(counts)))
dev.off()

pdf("6-leaf-general-tree-test-count-distribution.pdf",
    width = plot_width, height = plot_height)
# mgp adjusts space of axis labels
par(mgp = c(2.0,0.6,0))
# oma is space around all plots
par(oma = c(0,0,0,0), mar = c(3.1,3.1,1.5,0.5))

plot_binom_on_hist(counts = counts, p = binomial_prob)
mtext(bquote(chi^2 == .(chi_test_stat)),
      side=3, line = -1.5, adj = 0.02)
mtext(bquote(italic(p) == .(chi_test_p)),
      side=3, line = -2.5, adj = 0.02)
dev.off()


source("7-leaf-general-tree-test.r")

chi_test = chisq.test(counts)
chi_test_stat = format(chi_test$statistic)
chi_test_p = format(chi_test$p.value, digits = p_val_digits)

pdf("7-leaf-general-tree-test-topo-vs-count.pdf")
plot(x = d$topology, y = d$count, ylim = c(0, max(counts)))
dev.off()

pdf("7-leaf-general-tree-test-count-distribution.pdf",
    width = plot_width, height = plot_height)
# mgp adjusts space of axis labels
par(mgp = c(2.0,0.6,0))
# oma is space around all plots
par(oma = c(0,0,0,0), mar = c(3.1,3.1,1.5,0.5))

plot_binom_on_hist(counts = counts, p = binomial_prob)
mtext(bquote(chi^2 == .(chi_test_stat)),
      side=3, line = -1.5, adj = 0.02)
mtext(bquote(italic(p) == .(chi_test_p)),
      side=3, line = -2.5, adj = 0.02)
dev.off()


source("6-leaf-NeighborHeightNodeSwapAll-test.r")

chi_test = chisq.test(counts)
chi_test_stat = format(chi_test$statistic)
chi_test_p = format(chi_test$p.value, digits = p_val_digits)

pdf("6-leaf-NeighborHeightNodeSwapAll-test-topo-vs-count.pdf")
plot(x = d$topology, y = d$count, ylim = c(0, max(counts)))
dev.off()

pdf("6-leaf-NeighborHeightNodeSwapAll-test-count-distribution.pdf",
    width = plot_width, height = plot_height)
# mgp adjusts space of axis labels
par(mgp = c(2.0,0.6,0))
# oma is space around all plots
par(oma = c(0,0,0,0), mar = c(3.1,3.1,1.5,0.5))

plot_binom_on_hist(counts = counts, p = binomial_prob)
mtext(bquote(chi^2 == .(chi_test_stat)),
      side=3, line = -1.5, adj = 0.02)
mtext(bquote(italic(p) == .(chi_test_p)),
      side=3, line = -2.5, adj = 0.02)
dev.off()


source("6-leaf-NodeHeightSlideBumpSwapAllScaler-test.r")

chi_test = chisq.test(counts)
chi_test_stat = format(chi_test$statistic)
chi_test_p = format(chi_test$p.value, digits = p_val_digits)

pdf("6-leaf-NodeHeightSlideBumpSwapAllScaler-test-topo-vs-count.pdf")
plot(x = d$topology, y = d$count, ylim = c(0, max(counts)))
dev.off()

pdf("6-leaf-NodeHeightSlideBumpSwapAllScaler-test-count-distribution.pdf",
    width = plot_width, height = plot_height)
# mgp adjusts space of axis labels
par(mgp = c(2.0,0.6,0))
# oma is space around all plots
par(oma = c(0,0,0,0), mar = c(3.1,3.1,1.5,0.5))

plot_binom_on_hist(counts = counts, p = binomial_prob)
mtext(bquote(chi^2 == .(chi_test_stat)),
      side=3, line = -1.5, adj = 0.02)
mtext(bquote(italic(p) == .(chi_test_p)),
      side=3, line = -2.5, adj = 0.02)
dev.off()
