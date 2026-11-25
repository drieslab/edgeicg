# Pooling of the leak coefficients across receiver-sender pairs (same rule as lrmap's lrmap_leak()).
# Each usable estimate (lambda in [0, 1], mu in [-1, 1]) is shrunk toward the pooled value of its receiver type with
# weight support / (support + S0), support = the receivers' contact shares of the sender type, summed (about the
# number of fully touching receivers). Pairs without a usable estimate get the pooled value. Pooled value: median
# over the well-supported pairs (support >= S0) of the same receiver type, else over all well-supported pairs;
# pool = "q1" uses the first quartile for lambda (mu: median). pool = "none": pairs with support < S0 or without a
# usable estimate are not corrected (NA).
.leak_pool <- function(recv, lam, mu, support, pool = "median", S0 = 30) {
  usable <- !is.na(lam) & !is.na(mu) & lam >= 0 & lam <= 1 & mu >= -1 & mu <= 1
  good <- usable & !is.na(support) & support >= S0
  status <- ifelse(usable, ifelse(good, "ok", "low support"), ifelse(is.na(lam), "no estimate", "out of range"))
  if (pool == "none" || !any(good)) {
    if (pool != "none") warning("no receiver-sender pair has enough contact to pool leak coefficients; low-support pairs are not corrected")
    return(list(lambda = ifelse(good, lam, NA_real_), mu = ifelse(good, mu, NA_real_), status = status))
  }
  ql <- if (pool == "q1") function(x) stats::quantile(x, 0.25, names = FALSE) else stats::median
  by_r <- function(x, f) { g <- f(x[good]); vapply(recv, function(r) { w <- good & recv == r; if (any(w)) f(x[w]) else g }, 0) }
  pl <- by_r(lam, ql); pm <- by_r(mu, stats::median)
  s <- ifelse(usable & !is.na(support), support, 0); wt <- s / (s + S0)
  l2 <- ifelse(usable, wt * lam + (1 - wt) * pl, pl); m2 <- ifelse(usable, wt * mu + (1 - wt) * pm, pm)
  status[status != "ok"] <- paste0("pooled: ", status[status != "ok"])
  list(lambda = l2, mu = m2, status = status)
}
