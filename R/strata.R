# Region and local effects. The global fit gives one effect per receiver x sender x gene x scale. Here the same
# effect is estimated within strata of the receiver cells: user regions (`region =`) or a grid of square patches
# (`local = TRUE`). The leak model and the design stay global; within a stratum the effect is the slope of the
# partial residual on the partial exposure, with a stratum intercept (Frisch-Waugh):
#
#   x~_j  = exposure column j residualised on the other columns (weighted; from the fit: x~_j = M_j / (w G_jj))
#   y~    = r + b_j x~_j          (r: stage-2 residual, already leak-corrected)
#   b_d   = sum_d xc y~ / sum_d xc^2,   xc = x~_j - mean_d(x~_j)
#
# In the sqrt(weight) scale every quantity is a sum over the stratum's cells of r, r^2 and xc_j r, plus
# gene-independent sums (n_d, mean_d(x~_j), sum_d xc_j^2), so strata add one sum per column and two per stratum
# to pass 2. The standard error is model-based within the stratum (residual variance of the stratum fit), so
# the variance may differ between strata. The niche screen uses the same algebra with its residualised columns
# (Nt) and residual (Yt), b = 0, and the pair's calibrated variance inflation.
#
# Between strata: Cochran's Q = sum_d (b_d - b_bar)^2 / se_d^2 (df = strata - 1, b_bar the inverse-variance
# mean), tau^2 (DerSimonian-Laird) and, for patches, empirical-Bayes shrunken patch effects.

#' Square patches for local effects
#'
#' Assigns each cell to a square patch of a grid laid over the tissue. Used by `edgeicg(local = TRUE)`.
#'
#' @param coords n x 2 cell centroids (microns).
#' @param size patch side (microns).
#' @return list(code = patch index per cell, labels, centers = data.frame(patch, x, y)).
#' @export
edgeicg_patches <- function(coords, size = 1000) {
  coords <- as.matrix(coords)
  gx <- floor((coords[, 1] - min(coords[, 1])) / size); gy <- floor((coords[, 2] - min(coords[, 2])) / size)
  key <- paste(gx, gy, sep = "_"); lab <- sort(unique(key)); code <- match(key, lab)
  ce <- data.frame(patch = lab, x = as.vector(tapply(coords[, 1], code, mean)), y = as.vector(tapply(coords[, 2], code, mean)))
  list(code = code, labels = lab, centers = ce)
}

# Setup per receiver type (z: the type's design, after Ginv / M / Nt are known): stratum code per cell (NA:
# stratum with fewer than min_cells of the type), centred partial exposures for the contact / paracrine columns
# (XC) and the niche columns (NC), and their gene-independent sums.
.strata_type_setup <- function(z, code, D, min_cells) {
  g <- code[z$idx]
  nd <- tabulate(g[!is.na(g)], D); g[!is.na(g) & nd[g] < min_cells] <- NA
  nd <- tabulate(g[!is.na(g)], D)
  cp <- which(z$cols$scale %in% c("contact", "paracrine")) + 1L       # columns of A
  dg <- diag(z$Ginv)
  xs <- sweep(z$M[, cp, drop = FALSE] / z$sw, 2, dg[cp], "/")          # x~ in the sqrt(weight) scale
  cen <- function(X) {
    if (!ncol(X)) return(list(XC = X, xbar = matrix(0, D, 0), sxx = matrix(0, D, 0)))
    ok <- !is.na(g); Ind <- Matrix::sparseMatrix(i = g[ok], j = which(ok), x = 1, dims = c(D, nrow(X)))
    xbar <- as.matrix(Ind %*% X) / pmax(nd, 1)
    XC <- X; XC[ok, ] <- X[ok, , drop = FALSE] - xbar[g[ok], , drop = FALSE]; XC[!ok, ] <- 0
    list(XC = XC, xbar = xbar, sxx = as.matrix(Ind %*% XC^2))
  }
  a <- cen(xs); b <- cen(if (length(z$nsend)) z$Nt else matrix(0, length(g), 0))
  list(g = g, D = D, nd = nd, cp = cp, XC = a$XC, xbar = a$xbar, sxx = a$sxx,
       NC = b$XC, nxbar = b$xbar, nsxx = b$sxx)
}

# Kernel part (pass 2, residual stage): per stratum sums of r, r^2 and xc_j r. r: cells x genes (sqrt-weight scale).
.strata_sums <- function(g, XC, r, D) {
  ok <- !is.na(g)
  out <- list(r1 = matrix(0, D, ncol(r)), r2 = matrix(0, D, ncol(r)), c = matrix(0, D * ncol(XC), ncol(r)))
  if (!any(ok)) return(out)
  Ind <- Matrix::sparseMatrix(i = g[ok], j = seq_len(sum(ok)), x = 1, dims = c(D, sum(ok)))
  ro <- r[ok, , drop = FALSE]
  out$r1 <- as.matrix(Ind %*% ro); out$r2 <- as.matrix(Ind %*% ro^2)
  if (ncol(XC)) out$c <- do.call(rbind, lapply(seq_len(ncol(XC)), function(j) as.matrix(Ind %*% (ro * XC[ok, j]))))
  out
}

# Stratum effects and standard errors of one column: b (global coefficient per gene, 0 for the niche screen),
# n (D), xbar, sxx (D), sums r1, r2 (D x G) and c (D x G). Returns D x G matrices (NA where untestable).
.strata_fit <- function(b, n, xbar, sxx, r1, r2, c) {
  D <- length(n); G <- ncol(r1)
  B <- matrix(b, D, G, byrow = TRUE)
  ok <- n > 2 & sxx > 1e-12
  bd <- B + c / sxx
  sy <- r1 + B * (n * xbar)
  syy <- r2 + 2 * B * (c + xbar * r1) + B^2 * (sxx + n * xbar^2)
  ss <- syy - sy^2 / pmax(n, 1) - bd^2 * sxx
  se <- sqrt(pmax(ss, 0) / pmax(n - 2, 1) / sxx)
  bd[!ok, ] <- NA; se[!ok, ] <- NA
  list(b = bd, se = se)
}

# Heterogeneity across strata for one column: b, se (D x G, NA = untestable). Returns per gene Q, df, tau2, the
# inverse-variance mean, p (chi-square), and the shrunken stratum effects and their posterior sd.
.strata_het <- function(b, se) {
  w <- 1 / se^2; w[!is.finite(w)] <- NA; okm <- !is.na(b) & !is.na(w)
  w[!okm] <- 0; bb <- b; bb[!okm] <- 0
  sw <- colSums(w); m <- colSums(w * bb) / pmax(sw, EPS)
  Q <- colSums(w * sweep(bb, 2, m)^2 * okm); df <- colSums(okm) - 1
  cc <- sw - colSums(w^2) / pmax(sw, EPS)
  tau2 <- pmax(0, (Q - df) / pmax(cc, EPS)); tau2[df < 1] <- NA
  p <- stats::pchisq(Q, pmax(df, 1), lower.tail = FALSE); p[df < 1] <- NA
  s2 <- se^2; T2 <- matrix(tau2, nrow(b), ncol(b), byrow = TRUE)
  k <- T2 / (T2 + s2)
  bs <- sweep(sweep(b, 2, m) * k, 2, m, "+"); ss <- sqrt(T2 * s2 / (T2 + s2))
  bs[!okm] <- NA; ss[!okm] <- NA
  list(Q = Q, df = df, tau2 = tau2, mean = m, p = p, bs = bs, ss = ss)
}

# Exposed receivers per stratum for each pair and scale (exposure > exposed_at), and the 10-90% exposure range
# within each stratum (for the effect size).
.strata_exposure <- function(code, D, types, pairs, smat, tested, exposed_at, need_dx) {
  out <- list()
  for (s in tested) {
    ne <- matrix(0L, nrow(pairs), D); dx <- if (need_dx) matrix(NA_real_, nrow(pairs), D) else NULL
    for (p in seq_len(nrow(pairs))) {
      t <- pairs$receiver[p]; kk <- pairs$sender[p]; i <- which(types == t & !is.na(code))
      e <- smat[[s]][i, kk]
      ne[p, ] <- tabulate(code[i][e > exposed_at], D)
      if (need_dx) dx[p, ] <- vapply(seq_len(D), function(d) {
        v <- e[code[i] == d]; if (length(v) < 3) NA_real_ else diff(stats::quantile(v, c(0.1, 0.9), names = FALSE)) }, 1)
    }
    out[[s]] <- list(ne = ne, dx = dx)
  }
  out
}

#' Region-specific effects
#'
#' For a fit with `region =`, the effect of every testable entry within each region, and a test of whether the
#' effect differs between regions.
#'
#' Within a region the effect is the slope of the receiver's expression (leak-corrected, adjusted for the other
#' senders and scales of the global model) on its exposure to the sender, with a region intercept. Between
#' regions, Cochran's Q (df = regions - 1) tests heterogeneity; `q_het` is its BH-adjusted p within scale.
#'
#' `pattern`, for entries with `q_het < q_thr`: `opposite` (significant in a region, q < `q_thr` and
#' |lfc| >= `min_lfc`, and another region goes the other way with |z| >= 2), `specific` (otherwise significant in
#' one region only), `variable` (heterogeneous otherwise); `shared` when `q_het >= q_thr`.
#'
#' @param x an `EdgeICG` object fitted with `region =`.
#' @param scale scales to report.
#' @param q_thr threshold on the region q and on `q_het`.
#' @param min_lfc effect-size floor for a region to count as significant.
#' @param all FALSE (default): only entries that are heterogeneous or significant in at least one region.
#' @return data.table with one row per entry and region: receiver, sender, gene, scale, region, n_cells,
#'   n_exposed, effect, lfc, se, z, q, global_z, global_q, Q, df, q_het, pattern.
#' @export
edgeicg_region_results <- function(x, scale = c("contact", "paracrine", "niche"), q_thr = 0.05, min_lfc = 0.25, all = FALSE) {
  stopifnot(inherits(x, "EdgeICG"))
  if (is.null(x$region)) stop("the fit has no region effects; run edgeicg(region = ...)")
  .strata_table(x, x$region, match.arg(scale, several.ok = TRUE), q_thr, min_lfc, all, "region")
}

#' Local effects (experimental)
#'
#' Experimental: the Q test over patches is calibrated but has little power for smooth or half-plane patterns.
#' For a fit with `local = TRUE`: for every testable entry, whether its effect varies across patches of the
#' tissue (Cochran's Q over patches, BH within scale) and, for entries with `q_het < q_thr`, the effect per patch
#' with empirical-Bayes shrinkage towards the mean (`effect_shrunk`, `sd_shrunk`), to map where it is strongest.
#'
#' @param x an `EdgeICG` object fitted with `local = TRUE`.
#' @param scale scales to report.
#' @param q_thr threshold on `q_het`.
#' @return list(entries = data.table(receiver, sender, gene, scale, global_effect, global_z, Q, df, tau2, q_het,
#'   local_frac), patches = data.table(receiver, sender, gene, scale, patch, x, y, n_exposed, effect, se,
#'   effect_shrunk, sd_shrunk)). `local_frac`: share of the exposed receivers that sit in patches whose shrunken
#'   effect is more than 2 posterior sd from zero, with the sign of the strongest patch.
#' @export
edgeicg_local_results <- function(x, scale = c("contact", "paracrine", "niche"), q_thr = 0.05) {
  stopifnot(inherits(x, "EdgeICG"))
  if (is.null(x$local)) stop("the fit has no local effects; run edgeicg(local = TRUE)")
  scale <- match.arg(scale, several.ok = TRUE); lev <- x$types; pr <- x$pairs; L <- x$local
  ent <- data.table::rbindlist(lapply(scale, function(s) {
    h <- L$het[[s]]; w <- which(!is.na(h$q) & h$q < q_thr, arr.ind = TRUE)
    if (!nrow(w)) return(NULL)
    f <- if (s == "niche") x$screen else x$scale_fits[[s]]
    data.table::data.table(pi = w[, 1], gi = w[, 2], receiver = lev[pr$receiver[w[, 1]]], sender = lev[pr$sender[w[, 1]]],
                           gene = x$genes[w[, 2]], scale = s, global_effect = f$effect[w], global_z = f$z[w],
                           Q = h$Q[w], df = h$df[w], tau2 = h$tau2[w], q_het = h$q[w])
  }))
  if (!nrow(ent)) return(list(
    entries = data.table::data.table(receiver = character(), sender = character(), gene = character(), scale = character(),
      local_frac = numeric(), global_effect = numeric(), global_z = numeric(), Q = numeric(), df = numeric(),
      tau2 = numeric(), q_het = numeric()),
    patches = data.table::data.table(receiver = character(), sender = character(), gene = character(), scale = character(),
      patch = character(), x = numeric(), y = numeric(), n_exposed = integer(), effect = numeric(), se = numeric(),
      effect_shrunk = numeric(), sd_shrunk = numeric())))
  P <- L$patches[ent[, list(pi, gi, scale)], on = c("pi", "gi", "scale"), nomatch = NULL]
  P[, `:=`(x = L$centers$x[patch], y = L$centers$y[patch])]
  P[, sig := abs(effect_shrunk) > 2 * sd_shrunk]
  fr <- P[, {
    sg <- sign(effect_shrunk[which.max(abs(effect_shrunk / sd_shrunk))])
    list(local_frac = sum(n_exposed[sig & sign(effect_shrunk) == sg]) / max(sum(n_exposed), 1))
  }, by = list(pi, gi, scale)]
  ent <- fr[ent, on = c("pi", "gi", "scale")]
  P <- P[, list(pi, gi, scale, patch = L$centers$patch[patch], x, y, n_exposed, effect, se, effect_shrunk, sd_shrunk)]
  P <- ent[, list(pi, gi, scale, receiver, sender, gene)][P, on = c("pi", "gi", "scale")]
  ent[, c("pi", "gi") := NULL]; P[, c("pi", "gi") := NULL]
  data.table::setcolorder(ent, c("receiver", "sender", "gene", "scale"))
  data.table::setorder(ent, scale, q_het)
  list(entries = ent[], patches = P[])
}

# Region table (see edgeicg_region_results)
.strata_table <- function(x, R, scale, q_thr, min_lfc, all, what) {
  lev <- x$types; pr <- x$pairs; D <- length(R$labels)
  data.table::rbindlist(lapply(scale, function(s) {
    E <- R$effect[[s]]; S <- R$se[[s]]; h <- R$het[[s]]
    f <- if (s == "niche") x$screen else x$scale_fits[[s]]
    Z <- E / S; Qv <- array(NA_real_, dim(Z)); okz <- !is.na(Z)
    Qv[okz] <- stats::p.adjust(2 * stats::pnorm(-abs(Z[okz])), "BH")
    LF <- R$lfc[[s]]
    sig <- !is.na(Qv) & Qv < q_thr & !is.na(LF) & abs(LF) >= min_lfc
    anyd <- function(A) apply(A, c(1, 2), function(v) any(v, na.rm = TRUE))
    nsig <- apply(sig, c(1, 2), sum); psig <- anyd(sig & E > 0); nsg <- anyd(sig & E < 0)
    opp <- (psig & anyd(Z <= -2)) | (nsg & anyd(Z >= 2))
    het <- !is.na(h$q) & h$q < q_thr
    keep <- if (all) apply(okz, c(1, 2), any) else (het | nsig > 0)
    w <- which(keep, arr.ind = TRUE); if (!nrow(w)) return(NULL)
    pat <- ifelse(!het[w], "shared", ifelse(opp[w], "opposite", ifelse(nsig[w] == 1, "specific", "variable")))
    data.table::rbindlist(lapply(seq_len(D), function(d) {
      wd <- cbind(w, d)
      data.table::data.table(receiver = lev[pr$receiver[w[, 1]]], sender = lev[pr$sender[w[, 1]]], gene = x$genes[w[, 2]],
                             scale = s, region = R$labels[d], n_cells = R$n_cells[pr$receiver[w[, 1]], d],
                             n_exposed = R$n_exposed[[s]][w[, 1], d], effect = E[wd], lfc = LF[wd], se = S[wd], z = Z[wd],
                             q = Qv[wd], global_z = f$z[w], global_q = f$q[w], Q = h$Q[w], df = h$df[w], q_het = h$q[w],
                             pattern = pat)
    }))[!is.na(effect)]
  }))
}

# Engine finalisation for one strata set (nm: "region" or "local"; k: its index in the set list).
.strata_collect <- function(nm, k, S, STR, blocks, dz, fitted, pairs, G, K, types, smat, tested, exposed_at, fits,
                            okn, im, tmr, min_exposed_region, local_p) {
  D <- length(S$labels); P_n <- nrow(pairs); local <- nm == "local"
  exq <- .strata_exposure(S$code, D, types, pairs, smat, tested, exposed_at, need_dx = !local)
  if (!local) {
    eff <- lapply(tested, function(s) array(NA_real_, c(P_n, G, D))); names(eff) <- tested; ses <- eff
  } else {
    het <- lapply(tested, function(s) list(Q = matrix(NA_real_, P_n, G), df = matrix(NA_real_, P_n, G),
                                         tau2 = matrix(NA_real_, P_n, G), p = matrix(NA_real_, P_n, G)))
    names(het) <- tested; prow <- list()
  }
  ncell <- matrix(0L, K, D)
  put <- function(s, pp, gi, fb) {                                    # mask, then store (region) or summarise (local)
    okg <- if (s == "niche") okn[pp, gi] else fits[[s]]$testable[pp, gi]
    oks <- exq[[s]]$ne[pp, ] >= min_exposed_region
    fb$b[!oks, ] <- NA; fb$se[!oks, ] <- NA; fb$b[, !okg] <- NA; fb$se[, !okg] <- NA
    if (!local) { eff[[s]][pp, gi, ] <<- t(fb$b); ses[[s]][pp, gi, ] <<- t(fb$se); return(invisible()) }
    h <- .strata_het(fb$b, fb$se)
    het[[s]]$Q[pp, gi] <<- h$Q; het[[s]]$df[pp, gi] <<- h$df; het[[s]]$tau2[pp, gi] <<- h$tau2; het[[s]]$p[pp, gi] <<- h$p
    keep <- which(!is.na(h$p) & h$p < local_p)
    if (length(keep)) {
      w <- which(!is.na(fb$b[, keep, drop = FALSE]), arr.ind = TRUE)
      prow[[length(prow) + 1L]] <<- data.table::data.table(pi = pp, gi = gi[keep][w[, 2]], scale = s, patch = w[, 1],
        n_exposed = exq[[s]]$ne[pp, w[, 1]], effect = fb$b[, keep, drop = FALSE][w], se = fb$se[, keep, drop = FALSE][w],
        effect_shrunk = h$bs[, keep, drop = FALSE][w], sd_shrunk = h$ss[, keep, drop = FALSE][w])
    }
  }
  for (t in fitted) ncell[t, ] <- dz[[t]]$sgi[[k]]$nd
  for (bk in seq_along(blocks)) {
    gi <- blocks[[bk]]
    for (t in fitted) {
      z <- dz[[t]]; e <- STR[[bk]][[t]]; if (is.null(e)) next
      g <- z$sgi[[k]]; sm <- e$str[[k]]
      for (jj in seq_along(g$cp)) {
        cl <- z$cols[g$cp[jj] - 1L, ]; pp <- z$pp[match(cl$k, z$others)]
        rows <- (jj - 1L) * D + seq_len(D)
        put(cl$scale, pp, gi, .strata_fit(e$coef[jj, ], g$nd, g$xbar[, jj], g$sxx[, jj], sm$cp$r1, sm$cp$r2,
                                          sm$cp$c[rows, , drop = FALSE]))
      }
      for (j in seq_along(z$nsend)) {
        if (is.null(sm$nt)) break
        rows <- (j - 1L) * D + seq_len(D)
        fb <- .strata_fit(rep(0, length(gi)), g$nd, g$nxbar[, j], g$nsxx[, j], sm$nt$r1, sm$nt$r2, sm$nt$c[rows, , drop = FALSE])
        fb$se <- fb$se * sqrt(im[z$npp[j]])                            # the pair's calibrated variance inflation
        put("niche", z$npp[j], gi, fb)
      }
    }
  }
  bh <- function(p) { q <- p; ok <- !is.na(p); q[ok] <- stats::p.adjust(p[ok], "BH"); q }
  if (local) {
    for (s in tested) het[[s]]$q <- bh(het[[s]]$p)
    pt <- data.table::rbindlist(prow)
    if (!nrow(pt)) pt <- data.table::data.table(pi = integer(), gi = integer(), scale = character(), patch = integer(),
      n_exposed = integer(), effect = numeric(), se = numeric(), effect_shrunk = numeric(), sd_shrunk = numeric())
    return(list(labels = S$labels, centers = S$centers, n_cells = ncell, het = het, patches = pt,
                n_exposed = lapply(exq, `[[`, "ne")))
  }
  hr <- lfc <- list()
  for (s in tested) {
    Hq <- list(Q = matrix(NA_real_, P_n, G), df = matrix(NA_real_, P_n, G), tau2 = matrix(NA_real_, P_n, G), p = matrix(NA_real_, P_n, G))
    for (pp in seq_len(P_n)) {
      b <- t(matrix(eff[[s]][pp, , ], G, D)); se <- t(matrix(ses[[s]][pp, , ], G, D))
      if (all(is.na(b))) next
      h <- .strata_het(b, se)
      Hq$Q[pp, ] <- h$Q; Hq$df[pp, ] <- h$df; Hq$tau2[pp, ] <- h$tau2; Hq$p[pp, ] <- h$p
    }
    Hq$q <- bh(Hq$p); hr[[s]] <- Hq
    dx <- exq[[s]]$dx
    lfc[[s]] <- array(vapply(seq_len(D), function(d) .edgeicg_lfc(eff[[s]][, , d], dx[, d], tmr), matrix(0, P_n, G)), c(P_n, G, D))
  }
  list(labels = S$labels, n_cells = ncell, effect = eff, se = ses, lfc = lfc, het = hr,
       n_exposed = lapply(exq, `[[`, "ne"))
}
