plot_binom_on_hist <- function(counts, p,
        sample_line_type = 1,
        binom_line_type = 2,
        sample_line_wt = 1.5,
        binom_line_wt = 1.5,
        sample_line_color = "black",
        binom_line_color = "black") {
    k = seq(from = min(counts), to = max(counts), by = 1)
    n = sum(counts)
    binom_probs = dbinom(k, n, p)
    # hist(counts, freq = F)
    d = density(counts)
    plot(d$x, d$y,
            type = 'l',
            lwd = sample_line_wt,
            lty = sample_line_type,
            col = sample_line_color,
            xlab = "Number of times each topology sampled",
            ylab = "Density",
            main = NULL)
    lines(k, binom_probs, type = 'l',
            lwd = binom_line_wt,
            lty = binom_line_type,
            col = binom_line_color)
    num_topologies = length(counts)
    binom_label = paste("Binom(", format(n), ", 1/", format(num_topologies, big.mark = ","), ")", sep = "")
    legend("bottom",
            legend=c("MCMC", binom_label),
            lty=c(sample_line_type, binom_line_type),
            lwd=c(sample_line_wt, binom_line_wt),
            col=c(sample_line_color, binom_line_color),
            box.lty = 0,
            horiz = FALSE,
            inset = 0.99,
            xpd = TRUE)
}

plot_binom_qq <- function(counts, p) {
    probs = seq(from = 0.01, to = 0.99, by = 0.01)
    n = sum(counts)
    binom_quants = qbinom(probs, n, p)
    sample_quants = quantile(counts, probs)
    plot(binom_quants, sample_quants,
         xlab = "Binomial quantiles",
         ylab = "Sample quantiles",
         col = "gray45",
         panel.first = c(
            abline(a = 0, b = 1, col = "lightgray")
         )
    )
}
