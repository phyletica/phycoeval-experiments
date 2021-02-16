# plot_binom_on_hist <- function(counts, p) {
#     k = seq(from = min(counts), to = max(counts), by = 1)
#     n = sum(counts)
#     binom_probs = dbinom(k, n, p)
#     hist(counts, freq = F)
#     lines(k, binom_probs, type = 'l')
# }
plot_binom_on_hist <- function(counts, p) {
    k = seq(from = min(counts), to = max(counts), by = 1)
    n = sum(counts)
    binom_probs = dbinom(k, n, p)
    # hist(counts, freq = F)
    d = density(counts)
    plot(d$x, d$y,
            type = 'l',
            lwd = 1.5,
            lty = 1,
            xlab = "Number of times each topology sampled",
            ylab = "Density",
            main = NULL)
    lines(k, binom_probs, type = 'l', lty = 2, lwd = 1.5)
    num_topologies = length(counts)
    binom_label = paste("Binom(", format(n), ", 1/", format(num_topologies, big.mark = ","), ")", sep = "")
    legend("bottom",
            legend=c("MCMC", binom_label),
            lty=c(1,2),
            box.lty = 0,
            horiz = FALSE,
            inset = 0.99,
            xpd = TRUE)
}
