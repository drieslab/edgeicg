#' Simulate a small tissue with known neighbour effects and leakage
#'
#' Three cell types in smooth domains (cells about 10 microns apart). Ground truth the model never sees:
#' directed contact effects (multiplicative, receiver 1 <- sender 2 up on `S1`, receiver 2 <- sender 3 down
#' on `S2`), optional niche effects (receiver 1 genes `S3` respond to the regional density of type 2), and
#' mass-conserving spillover between neighbours over a different graph than the model uses.
#'
#' @param n cells.
#' @param G genes.
#' @param seed seed.
#' @param kappa mean leaked fraction.
#' @param gamma contact log effect at full exposure.
#' @param niche_gamma niche log effect.
#' @param spacing cell spacing (microns).
#' @param region_gamma contact log effect of region-dependent genes (0 = none). Regions are the left and right
#'   halves of the tissue (`region`). Receiver 1 <- sender 2: genes `S4` respond only in the right half, genes
#'   `S5` respond in opposite directions (up on the left, down on the right). Needs G >= 150.
#' @return list(counts (genes x cells), coords, types, truth (data.table), lib, region (per cell); `truth` gains
#'   a `region` column ("all", "right" or "opposite")).
#' @export
simulate_edgeicg <- function(n = 4000, G = 200, seed = 1, kappa = 0.10, gamma = 1, niche_gamma = 0, spacing = 10,
                             region_gamma = 0) {
  set.seed(seed); K <- 3L
  L <- sqrt(n) * spacing
  coords <- cbind(stats::runif(n, 0, L), stats::runif(n, 0, L))
  f <- numeric(n)
  for (m in 1:6) {
    ang <- stats::runif(1, 0, pi); wl <- stats::runif(1, 25, 60) * spacing; ph <- stats::runif(1, 0, 2 * pi)
    f <- f + sin(2 * pi * (coords %*% c(cos(ang), sin(ang))) / wl + ph)
  }
  types <- findInterval(f, stats::quantile(f, c(0.45, 0.80))) + 1L
  fl <- stats::runif(n) < 0.10; types[fl] <- sample.int(K, sum(fl), replace = TRUE)
  p0 <- stats::rlnorm(G, 0, 1); p0 <- p0 / sum(p0)
  prof <- matrix(p0, K, G, byrow = TRUE); free <- sample.int(G)
  for (t in 1:K) { mk <- free[((t - 1) * 15 + 1):(t * 15)]; prof[t, mk] <- prof[t, mk] * 10 }
  prof <- prof / rowSums(prof)
  pool <- setdiff(order(-p0)[seq_len(round(G * 0.8))], free[1:(K * 15)])
  if (length(pool) < 30) stop("G is too small: need at least 30 non-marker genes for the effects (use G >= 100)")
  S1 <- pool[1:10]; S2 <- pool[11:20]; S3 <- pool[21:30]
  nn <- RANN::nn2(coords, k = 7); d <- nn$nn.dists[, -1]; j <- nn$nn.idx[, -1]
  w <- exp(-d / (0.8 * stats::median(d))); w <- w / rowSums(w)
  E1 <- sapply(seq_len(K), function(kk) rowSums(w * (types[j] == kk)))
  logfac <- matrix(0, n, G)
  logfac[types == 1, S1] <- gamma * E1[types == 1, 2]
  logfac[types == 2, S2] <- -gamma * E1[types == 2, 3]
  right <- coords[, 1] > L / 2
  if (region_gamma != 0) {                                           # no extra random draws: other effects unchanged
    if (length(pool) < 50) stop("region_gamma needs at least 50 non-marker genes (use G >= 150)")
    S4 <- pool[31:40]; S5 <- pool[41:50]; r1 <- types == 1
    logfac[r1 & right, S4] <- region_gamma * E1[r1 & right, 2]
    logfac[r1, S5] <- region_gamma * ifelse(right[r1], -1, 1) * E1[r1, 2]
  }
  if (niche_gamma != 0) {
    nb <- RANN::nn2(coords, k = 200); wv <- exp(-0.5 * nb$nn.dists^2 / (6 * spacing)^2)
    niche <- rowSums(wv * (types[nb$nn.idx] == 2)) / rowSums(wv)
    logfac[types == 1, S3] <- logfac[types == 1, S3] + niche_gamma * niche[types == 1]
  }
  sizes <- c(150, 350, 100)
  N <- sizes[types] * prof[types, ] * exp(logfac)
  ns <- RANN::nn2(coords, k = 11); dm <- ns$nn.dists[, -1]; im <- ns$nn.idx[, -1]
  ws <- exp(-dm / (0.7 * stats::median(dm))); ws <- ws / rowSums(ws)
  kap <- stats::rbeta(n, 2, 2 * (1 - kappa) / kappa)
  spill <- matrix(0, n, G); for (jj in 1:10) spill <- spill + ws[, jj] * N[im[, jj], ]
  mu <- (1 - kap) * N + kap * spill
  cnt <- matrix(stats::rpois(length(mu), as.vector(mu)), n, G)
  ok <- rowSums(cnt) > 0
  genes <- sprintf("g%03d", seq_len(G))
  counts <- Matrix::Matrix(t(cnt[ok, ]), sparse = TRUE, dimnames = list(genes, sprintf("cell%05d", which(ok))))
  truth <- data.table::data.table(receiver = c(rep("t1", 10), rep("t2", 10), rep("t1", 10)),
                                  sender = c(rep("t2", 10), rep("t3", 10), rep("t2", 10)),
                                  gene = genes[c(S1, S2, S3)], scale = rep(c("contact", "contact", "niche"), each = 10),
                                  dir = rep(c(1, -1, sign(niche_gamma)), each = 10))
  if (niche_gamma == 0) truth <- truth[scale != "niche"]
  truth[, region := "all"]
  if (region_gamma != 0) truth <- rbind(truth, data.table::data.table(receiver = "t1", sender = "t2", gene = genes[c(S4, S5)],
    scale = "contact", dir = sign(region_gamma), region = rep(c("right", "opposite"), each = 10)))
  list(counts = counts, coords = coords[ok, ], types = paste0("t", types[ok]), truth = truth,
       lib = Matrix::colSums(counts), region = ifelse(right[ok], "right", "left"))
}
