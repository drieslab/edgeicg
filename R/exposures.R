#' Neighbourhood exposures at three scales
#'
#' Every exposure is a *share*: for each cell, the fraction of its neighbourhood that belongs to each cell
#' type (rows sum to 1), so local cell density does not enter.
#'
#' * **contact** (`C`): share of the k nearest neighbours' weight, `w = exp(-d / median(d))`.
#' * **paracrine** (`P`): Gaussian (`sigma_p`) share among cells within `3 * sigma_p` that are not contact
#'   neighbours.
#' * **niche** (`N`): Gaussian (`sigma_n`) share from a per-type cell-count grid (bin `bin`) blurred with a
#'   separable filter, read back at each cell's bin. The grids are kept for the niche calibration.
#'
#' @param coords n x 2 matrix of cell centroids in microns.
#' @param types integer cell-type codes 1..K.
#' @param K number of types.
#' @param k contact neighbours.
#' @param sigma_p,sigma_n paracrine and niche kernel widths (microns).
#' @param bin niche grid bin (microns).
#' @param max_k cap on paracrine neighbours per cell.
#' @return list(contact = edge list (recv, send, w), C, P, N (n x K), fields (niche grids), settings).
#' @export
edgeicg_exposures <- function(coords, types, K = max(types), k = 8, sigma_p = 25, sigma_n = 150, bin = 10,
                              max_k = 400) {
  coords <- as.matrix(coords); types <- as.integer(types); n <- nrow(coords)
  ct <- .contact_knn(coords, k)
  P <- .paracrine(coords, types, K, ct, sigma_p, max_k)
  fl <- .niche_fields(coords, types, K, sigma_n, bin)
  list(contact = ct, C = .share(ct$recv, ct$send, ct$w, types, K, n), P = P,
       N = fl$Bc / pmax(fl$tot, EPS), fields = fl,
       settings = list(k = k, sigma_p = sigma_p, sigma_n = sigma_n, bin = bin))
}

.contact_knn <- function(coords, k = 8) {
  nn <- RANN::nn2(coords, k = k + 1)
  idx <- nn$nn.idx[, -1, drop = FALSE]; d <- as.vector(nn$nn.dists[, -1, drop = FALSE])
  list(recv = rep(seq_len(nrow(coords)), k), send = as.vector(idx), w = exp(-d / stats::median(d)), kind = "knn")
}

.share <- function(recv, send, w, types, K, n) {
  S <- matrix(0, n, K)
  s <- rowsum(w, recv + (types[send] - 1L) * n)
  S[as.integer(rownames(s))] <- s[, 1]
  S / pmax(rowSums(S), EPS)
}

.paracrine <- function(coords, types, K, contact, sigma = 25, max_k = 400) {
  n <- nrow(coords)
  nn <- RANN::nn2(coords, k = min(max_k, n), searchtype = "radius", radius = 3 * sigma)
  recv <- rep(seq_len(n), ncol(nn$nn.idx)); send <- as.vector(nn$nn.idx); dist <- as.vector(nn$nn.dists)
  ok <- send > 0 & send != recv
  recv <- recv[ok]; send <- send[ok]; dist <- dist[ok]
  key <- recv * (n + 1) + send                                     # numeric keys (exact below 2^53)
  ck <- contact$recv * (n + 1) + contact$send
  keep <- !(key %in% ck)                                           # exclude contact neighbours
  .share(recv[keep], send[keep], exp(-0.5 * dist[keep]^2 / sigma^2), types, K, n)
}

# per-type Gaussian-blurred count grids; Bc = values at each cell's bin, tot = their sum
.niche_fields <- function(coords, types, K, sigma = 150, bin = 10) {
  gx <- floor((coords[, 1] - min(coords[, 1])) / bin) + 1; gy <- floor((coords[, 2] - min(coords[, 2])) / bin) + 1
  nx <- max(gx); ny <- max(gy)
  r <- ceiling(3 * sigma / bin); ker <- exp(-0.5 * ((-r):r * bin)^2 / sigma^2); ker <- ker / sum(ker)
  blur1 <- function(M) {                                            # rows then columns, zero padded
    pad <- function(v) c(rep(0, r), v, rep(0, r))
    M <- t(apply(M, 1, function(v) stats::filter(pad(v), ker, sides = 2)[(r + 1):(r + length(v))]))
    if (ny == 1) M <- t(M)
    apply(M, 2, function(v) stats::filter(pad(v), ker, sides = 2)[(r + 1):(r + length(v))])
  }
  Fk <- lapply(seq_len(K), function(k) {
    M <- matrix(0, nx, ny)
    tab <- table(factor(gx[types == k], seq_len(nx)), factor(gy[types == k], seq_len(ny)))
    M[] <- as.vector(tab); blur1(M)
  })
  Bc <- vapply(Fk, function(F) F[cbind(gx, gy)], numeric(length(gx)))
  if (is.null(dim(Bc))) Bc <- matrix(Bc, ncol = K)
  list(gx = gx, gy = gy, nx = nx, ny = ny, F = Fk, Bc = Bc, tot = rowSums(Bc), sigma = sigma, bin = bin)
}

#' Touching-cell edges from cell outlines
#'
#' For every outline vertex, vertices of other cells within `tau` microns count as touching; the vertex
#' contributes half the length of its two adjacent outline segments. The shared length of (i, j) is the
#' mean of both directions. Used as the spillover graph (which neighbours can leak transcripts).
#'
#' @param poly data.frame with `cell` (index 1..n into the analysis cells), `x`, `y`: outline vertices in
#'   drawing order per cell.
#' @param n number of cells.
#' @param tau touching distance (microns).
#' @param max_k neighbouring vertices searched per vertex.
#' @return data.frame of directed edges `recv`, `send`, `w` (shared length, microns).
#' @export
contact_from_polygons <- function(poly, n, tau = 2, max_k = 6) {
  p <- data.table::as.data.table(poly)[, list(cell = as.integer(cell), x = as.double(x), y = as.double(y))]
  p[, o := seq_len(.N)]; data.table::setorder(p, cell, o)
  p[, `:=`(xn = data.table::shift(x, -1L), yn = data.table::shift(y, -1L),
           xp = data.table::shift(x, 1L), yp = data.table::shift(y, 1L)), by = cell]
  p[, `:=`(xn = data.table::fifelse(is.na(xn), x[1], xn), yn = data.table::fifelse(is.na(yn), y[1], yn),
           xp = data.table::fifelse(is.na(xp), x[.N], xp), yp = data.table::fifelse(is.na(yp), y[.N], yp)), by = cell]
  seglen <- 0.5 * (sqrt((p$xn - p$x)^2 + (p$yn - p$y)^2) + sqrt((p$xp - p$x)^2 + (p$yp - p$y)^2))
  nn <- RANN::nn2(cbind(p$x, p$y), k = max_k + 1, searchtype = "radius", radius = tau)
  ci <- p$cell; rows <- rep(seq_len(nrow(p)), max_k + 1); nbv <- as.vector(nn$nn.idx)
  ok <- nbv > 0 & ci[rows] != ci[pmax(nbv, 1)]
  d <- unique(data.table::data.table(i = ci[rows[ok]], j = ci[nbv[ok]], v = rows[ok]), by = c("v", "j"))
  d[, len := seglen[v]]
  L <- d[, list(len = sum(len)), by = list(a = pmin(i, j), b = pmax(i, j))]
  S <- L[, list(w = sum(len) / 2), by = list(a, b)]
  data.frame(recv = c(S$a, S$b), send = c(S$b, S$a), w = c(S$w, S$w))
}
