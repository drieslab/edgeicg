# Niche-scale statistics shared by the engine.

# Estimates and HC3 covariance for models with columns S (1 or 2) of Nt (residualised niche columns), all
# genes (columns of Yt) at once. Leverage decomposes: h = h0 (base design) + h(Nt[, S]).
.niche_ols <- function(Nt, Yt, h0, S) {
  Z <- Nt[, S, drop = FALSE]; Gi <- solve(crossprod(Z)); B <- Gi %*% crossprod(Z, Yt)
  M <- Z %*% Gi
  h <- pmin(h0 + rowSums(M * Z), 0.99)
  R2 <- ((Yt - Z %*% B) / (1 - h))^2
  out <- list(b = B, v11 = colSums(M[, 1]^2 * R2))
  if (length(S) > 1) { out$v22 <- colSums(M[, 2]^2 * R2); out$v12 <- colSums(M[, 1] * M[, 2] * R2) }
  out
}

# log approximate Bayes factor (normal approximation, multivariate) vs no effect, Gaussian prior with sd tau on the raw
# effect, averaged over the prior grid (rows of tau1 / tau2).
.labf <- function(f, tau1, tau2 = NULL) {
  one <- function(k) {
    if (is.null(f$v22)) { V <- f$v11; T <- tau1[k, ]^2; b <- f$b[1, ]
      return(0.5 * log(V / (V + T)) + 0.5 * b^2 * (1 / V - 1 / (V + T))) }
    b1 <- f$b[1, ]; b2 <- f$b[2, ]; a <- f$v11; d <- f$v22; c <- f$v12
    dV <- a * d - c^2; a2 <- a + tau1[k, ]^2; d2 <- d + tau2[k, ]^2; dW <- a2 * d2 - c^2
    qV <- (b1^2 * d - 2 * b1 * b2 * c + b2^2 * a) / dV; qW <- (b1^2 * d2 - 2 * b1 * b2 * c + b2^2 * a2) / dW
    0.5 * log(dV / dW) + 0.5 * (qV - qW)
  }
  L <- sapply(seq_len(nrow(tau1)), one); if (is.null(dim(L))) L <- matrix(L, nrow = 1)
  mx <- apply(L, 1, max); mx + log(rowMeans(exp(L - mx)))
}

# HC3 z for each column of Z (shifted exposures, weighted and residualised like the design), all genes:
# b = Z'Y / Z'Z; variance expanded so that every shift is three matrix products.
.shift_z <- function(Z, Y, Y2, h_base) {
  zz <- colSums(Z^2); b <- crossprod(Z, Y) / zz
  h <- pmin(h_base + sweep(Z^2, 2, zz, "/"), 0.99)
  W <- sweep(Z, 2, zz, "/")^2 / (1 - h)^2
  v <- crossprod(W, Y2) - 2 * b * crossprod(W * Z, Y) + b^2 * colSums(W * Z^2)
  b / sqrt(pmax(v, EPS))
}

# Huber regression across genes of a pair's contact coefficients on the touching-sender profile (lambda) and
# the receiver's own profile (mu); lambda >= 0.
.decompose <- function(beta, se, ak, at, huber = 1.345, iters = 50) {
  w0 <- 1 / (se^2 + stats::median(se^2) + EPS)
  Z <- cbind(ak, at)
  fit <- function(wh) {
    ww <- w0 * wh
    co <- as.vector(solve(crossprod(Z * ww, Z), crossprod(Z * ww, beta)))
    if (co[1] < 0) co <- c(0, sum(ww * at * beta) / (sum(ww * at^2) + EPS))
    co
  }
  wh <- rep(1, length(beta)); co <- fit(wh)
  for (it in seq_len(iters)) {
    r <- as.vector(beta - Z %*% co) / (se + EPS)
    s <- 1.4826 * stats::median(abs(r - stats::median(r))) + 1e-8
    u <- abs(r) / s
    wh <- ifelse(u <= huber, 1, huber / u)
    new <- fit(wh)
    if (max(abs(new - co)) < 1e-8) { co <- new; break }
    co <- new
  }
  list(lam = co[1], mu = co[2])
}
