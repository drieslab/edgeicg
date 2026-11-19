# Gene-blocked engine. Model per receiver type t, gene g, cell i (composition scale y = counts / library * 100):
#
#   y_ig = a_tg + sum_k [ C_ik g^c_tkg + P_ik g^p_tkg + N_ik g^n_tkg ] + S_ig + e_ig
#   S_ig = sum_{j touching i, type(j) = k != t} c_ij lambda_tk x_jg + c^mu_i a_tg        (leaked transcripts)
#
# lambda_tk, mu_tk: Huber regression across genes of the stage-1 contact coefficients on the touching-sender
# profile and the receiver profile. S is removed per cell and the model refitted (WLS, weights library / mean
# library, HC3). Niche: marginal screen per sender (Frisch-Waugh on the contact + paracrine design), calibrated
# by a toroidal-shift null (pair-level variance inflation), then exact attribution over 1- and 2-sender models.
#
# Passes over the counts (gene blocks):
#   1. stage-1 fit per type -> a, contact b1 / se1, touching-sender profiles, detection, type means
#   2. leak removal, stage-2 fit (effects, HC3, omnibus), niche screen, shift-null statistics
#   3. attribution, re-reading only genes with a screen hit
# Gene-independent quantities (Gram inverse, leverage, projections of niche and shifted exposures) are computed
# once per receiver type and reused for every block.

.edgeicg_engine <- function(pv, coords, types, lib, ex, spill_edges, block_genes, min_cells, min_det, min_exposed,
                            exposed_at, huber, scale, calibrate, B, min_shift, seed, calib_cond, attribute, q_thr,
                            max_size, prior_sd, cs_level, workers, drop_degenerate, verbose) {
  t0 <- proc.time()[[3]]; tm <- c()
  genes <- pv$genes; G <- length(genes); n <- length(types); K <- max(types)
  tested <- c("contact", "paracrine", "niche")
  blocks <- split(seq_len(G), ceiling(seq_len(G) / block_genes))
  cw <- lib / mean(lib)

  if (calibrate) {                                                  # toroidal shifts
    set.seed(seed); fl <- ex$fields; bin <- fl$bin
    maxd <- sqrt((fl$nx %/% 2)^2 + (fl$ny %/% 2)^2) * bin
    if (maxd < min_shift) stop(sprintf("the tissue is too small for min_shift = %.0f um (largest wrap-around shift %.0f um); lower min_shift or set calibrate = FALSE", min_shift, maxd))
    sh <- matrix(0L, 0, 2); tries <- 0L
    while (nrow(sh) < B) {
      tries <- tries + 1L; if (tries > 1e5) stop("could not sample enough shifts of at least min_shift")
      d <- c(sample.int(fl$nx, 1) - 1L, sample.int(fl$ny, 1) - 1L); w <- pmin(d, c(fl$nx, fl$ny) - d) * bin
      if (sqrt(sum(w^2)) >= min_shift) sh <- rbind(sh, d)
    }
    sh0 <- rbind(c(0L, 0L), sh)
  }

  ct <- if (is.null(spill_edges)) ex$contact else spill_edges
  W <- .rsum(ct$w, ct$recv, n); wn <- ct$w / pmax(W[ct$recv], EPS)
  tr <- types[ct$recv]; ts <- types[ct$send]; het <- tr != ts
  smat <- list(contact = ex$C, paracrine = ex$P, niche = ex$N)

  # designs per receiver type; degenerate (constant) exposure columns dropped
  dz <- list(); pairs <- NULL; n_dropped <- 0L
  for (t in seq_len(K)) {
    idx <- which(types == t); if (length(idx) < min_cells) next
    others <- setdiff(seq_len(K), t)
    cols <- do.call(rbind, lapply(tested, function(s) data.frame(scale = s, k = others)))
    D <- do.call(cbind, lapply(tested, function(s) smat[[s]][idx, others, drop = FALSE]))
    if (drop_degenerate) {
      keep <- apply(D, 2, stats::var) > 1e-12
      n_dropped <- n_dropped + sum(!keep); D <- D[, keep, drop = FALSE]; cols <- cols[keep, , drop = FALSE]
    }
    dz[[t]] <- list(idx = idx, D = D, cols = cols, others = others)
    pairs <- rbind(pairs, data.frame(receiver = t, sender = others))
  }
  fitted <- which(!vapply(seq_len(K), function(t) t > length(dz) || is.null(dz[[t]]), TRUE))
  P_n <- nrow(pairs)
  for (t in fitted) {
    z <- dz[[t]]; idx <- z$idx
    A <- cbind(1, z$D); Aw <- A * cw[idx]; Pm <- diag(ncol(A)); Pm[1, 1] <- 0
    Ginv <- chol2inv(chol(crossprod(Aw, A) + 1e-6 * Pm)); M <- Aw %*% Ginv
    z$A <- A; z$Aw <- Aw; z$Ginv <- Ginv; z$M <- M; z$h <- pmin(rowSums(M * A), 0.99)
    z$pp <- vapply(z$others, function(kk) which(pairs$receiver == t & pairs$sender == kk), 1L)
    # column of each scale per sender (NA = dropped), and the sender's niche column within Nt
    z$ci <- lapply(z$others, function(kk) vapply(tested, function(s) {
      w <- which(z$cols$scale == s & z$cols$k == kk); if (length(w)) 1L + w else NA_integer_ }, 1L))
    nc <- which(z$cols$scale == "niche"); bc <- which(z$cols$scale != "niche")
    z$sw <- sqrt(cw[idx]); z$q <- qr(cbind(1, z$D[, bc, drop = FALSE]) * z$sw)
    z$Nt <- qr.resid(z$q, z$D[, nc, drop = FALSE] * z$sw); z$h0 <- rowSums(qr.Q(z$q)^2)
    z$nsend <- z$cols$k[nc]                                         # sender of each Nt column
    z$jn <- match(z$others, z$nsend)                                 # Nt column per sender (NA = none)
    z$npp <- z$pp[match(z$nsend, z$others)]                          # pair row per Nt column
    if (calibrate && length(nc)) {
      z$Zs <- z$qo <- z$ZsC <- z$hc <- vector("list", length(nc))
      for (j in seq_along(nc)) {
        kk <- z$nsend[j]
        raw <- vapply(seq_len(nrow(sh0)), function(b) {
          bk <- fl$F[[kk]][cbind((fl$gx[idx] - 1L + sh0[b, 1]) %% fl$nx + 1L, (fl$gy[idx] - 1L + sh0[b, 2]) %% fl$ny + 1L)]
          bk / pmax(fl$tot[idx] - fl$Bc[idx, kk] + bk, EPS) }, numeric(length(idx)))
        z$Zs[[j]] <- qr.resid(z$q, raw * z$sw)
        if (calib_cond) {
          z$qo[[j]] <- qr(z$Nt[, -j, drop = FALSE]); z$ZsC[[j]] <- qr.resid(z$qo[[j]], z$Zs[[j]])
          z$hc[[j]] <- z$h0 + rowSums(qr.Q(z$qo[[j]])^2)
        }
      }
    }
    dz[[t]] <- z
  }
  ep <- match(paste(tr, ts), paste(pairs$receiver, pairs$sender)); ue <- het & !is.na(ep)
  sw_p <- .rsum(wn[ue], ep[ue], P_n)
  E <- Matrix::sparseMatrix(i = ep[ue], j = ct$send[ue], x = wn[ue] / sw_p[ep[ue]], dims = c(P_n, n))
  .say(verbose, t0, "setup: %d genes in %d blocks, %d receiver types, %d pairs, %d degenerate columns dropped",
       G, length(blocks), length(fitted), P_n, n_dropped)

  st <- list(pv = pv, lib = lib, scale = scale, dz = dz, fitted = fitted, types = types, E = E, K = K, P_n = P_n,
             tested = tested, calibrate = calibrate, calib_cond = calib_cond)
  run <- .edgeicg_runner(workers, st)
  on.exit(run$close(), add = TRUE)

  # pass 1
  r1 <- run$map(blocks, .edgeicg_p1)
  a <- matrix(NA_real_, K, G); b1 <- se1 <- mbar <- matrix(NA_real_, P_n, G); detn <- tsum <- matrix(0, K, G)
  for (bi in seq_along(blocks)) { gi <- blocks[[bi]]; r <- r1[[bi]]
    a[, gi] <- r$a; b1[, gi] <- r$b1; se1[, gi] <- r$se1; mbar[, gi] <- r$mbar; detn[, gi] <- r$detn; tsum[, gi] <- r$tsum }
  rm(r1); mbar[sw_p == 0, ] <- NA
  tm["pass1"] <- proc.time()[[3]] - t0; .say(verbose, t0, "pass 1 done")
  lam <- mu <- rep(NA_real_, P_n)
  for (p in seq_len(P_n)) {
    t <- pairs$receiver[p]
    if (anyNA(mbar[p, ]) || anyNA(a[t, ]) || anyNA(b1[p, ])) next
    dc <- tryCatch(.decompose(b1[p, ], se1[p, ], mbar[p, ], a[t, ], huber), error = function(e) list(lam = 0, mu = 0))
    lam[p] <- dc$lam; mu[p] <- dc$mu
  }
  lam_mat <- mu_mat <- matrix(0, K, K)
  lam_mat[cbind(pairs$receiver, pairs$sender)] <- ifelse(is.na(lam), 0, lam)
  mu_mat[cbind(pairs$receiver, pairs$sender)] <- ifelse(is.na(mu), 0, mu)
  lam_e <- lam_mat[cbind(tr, ts)]; lam_e[!het] <- 0; mu_e <- mu_mat[cbind(tr, ts)]; mu_e[!het] <- 0
  Wl <- Matrix::sparseMatrix(i = ct$recv, j = ct$send, x = wn * lam_e, dims = c(n, n))
  cmu <- .rsum(wn * mu_e, ct$recv, n)
  a0 <- a; a0[is.na(a0)] <- 0

  # pass 2
  run$set(list(Wl = Wl, cmu = cmu, a0 = a0))
  r2 <- run$map(blocks, .edgeicg_p2)
  est <- lapply(tested, function(s) list(g = matrix(NA_real_, P_n, G), se = matrix(NA_real_, P_n, G))); names(est) <- tested
  omni <- bM <- seM <- Sm <- Sc <- zc0 <- matrix(NA_real_, P_n, G)
  for (bi in seq_along(blocks)) { gi <- blocks[[bi]]; r <- r2[[bi]]
    for (si in seq_along(tested)) { est[[si]]$g[, gi] <- r$g[[si]]; est[[si]]$se[, gi] <- r$se[[si]] }
    omni[, gi] <- r$omni; bM[, gi] <- r$bM; seM[, gi] <- r$seM; Sm[, gi] <- r$Sm; Sc[, gi] <- r$Sc; zc0[, gi] <- r$zc0 }
  rm(r2); run$close()
  tm["pass2"] <- proc.time()[[3]] - t0 - sum(tm); .say(verbose, t0, "pass 2 done")

  # masks and tests (BH within each scale)
  nt <- tabulate(types, K); det_t <- detn / pmax(nt, 1); det_t[setdiff(seq_len(K), fitted), ] <- 0
  det_ok <- det_t[pairs$receiver, , drop = FALSE] >= min_det
  fits <- list()
  for (s in tested) {
    ne <- mapply(function(t, kk) sum(types == t & smat[[s]][, kk] > exposed_at), pairs$receiver, pairs$sender)
    ok <- det_ok & ne >= min_exposed & !is.na(est[[s]]$g)
    zz <- est[[s]]$g / (est[[s]]$se + EPS); zz[!ok] <- NA
    fits[[s]] <- list(effect = est[[s]]$g, se = est[[s]]$se, z = zz, q = matrix(stats::p.adjust(2 * stats::pnorm(-abs(zz)), "BH"), P_n),
                      testable = ok)
    pairs[[paste0("n_exposed_", s)]] <- ne
  }
  okc <- fits$contact$testable | fits$paracrine$testable
  pv_o <- stats::pchisq(omni, df = 2, lower.tail = FALSE); pv_o[!okc] <- NA
  omni_q <- matrix(stats::p.adjust(pv_o, "BH"), P_n)
  pairs$lambda <- lam; pairs$mu <- mu

  # calibration (pair-level variance inflation from the shift null) and the niche screen
  okn <- fits$niche$testable & !is.na(bM); cal <- NULL
  pair_med <- function(S) { S[!okn] <- NA; cp <- suppressWarnings(apply(S, 1, stats::median, na.rm = TRUE)); cp[!is.finite(cp)] <- 1; cp }
  if (calibrate) {
    cm <- pair_med(Sm); cc <- if (calib_cond) pair_med(Sc) else rep(NA_real_, P_n)
    dd <- abs(zc0 - fits$niche$z); dd[!okn] <- NA
    cal <- list(infl_marg = pmax(cm, 1), infl_cond = if (calib_cond) pmax(cc, 1) else rep(1, P_n), raw_marg = cm, raw_cond = cc,
                shifts = sh, check_shift0 = if (calib_cond) suppressWarnings(max(dd, na.rm = TRUE)) else NA_real_)
  }
  im <- if (is.null(cal)) rep(1, P_n) else cal$infl_marg; ic <- if (is.null(cal)) rep(1, P_n) else cal$infl_cond
  seM <- seM * sqrt(im); zM <- bM / (seM + EPS); zM[!okn] <- NA
  qM <- matrix(stats::p.adjust(2 * stats::pnorm(-abs(zM)), "BH"), P_n)
  zC <- fits$niche$z / sqrt(ic)
  screen <- list(effect = bM, se = seM, z = zM, q = qM, testable = okn)
  cond <- list(z = zC, q = matrix(stats::p.adjust(2 * stats::pnorm(-abs(zC)), "BH"), P_n))

  # pass 3: attribution on screen-hit genes (main process)
  Att <- NULL
  if (attribute) {
    resid_block <- function(gi) { X <- .edgeicg_readX(st, gi); X - as.matrix(Wl %*% X) - cmu * a0[types, gi, drop = FALSE] }
    hit_any <- which(colSums(qM < q_thr, na.rm = TRUE) > 0); att <- list()
    for (gi in split(hit_any, ceiling(seq_along(hit_any) / block_genes))) {
      if (!length(gi)) next
      Yres <- resid_block(gi)
      for (t in fitted) {
        z <- dz[[t]]; if (!length(z$nsend)) next
        rows <- which(pairs$receiver == t)
        hl <- which(colSums(qM[rows, gi, drop = FALSE] < q_thr, na.rm = TRUE) > 0); if (!length(hl)) next
        hitg <- gi[hl]; Yt <- qr.resid(z$q, Yres[z$idx, hl, drop = FALSE] * z$sw)
        okk <- okn[z$npp, hitg, drop = FALSE]
        cand <- which(rowSums(okk) > 0); J <- length(cand); if (!J) next
        sy <- sqrt(colMeans(Yt^2)); sx <- sqrt(colMeans(z$Nt^2))
        tau <- function(j) outer(prior_sd, sy / sx[j])
        models <- c(lapply(cand, function(j) j), if (max_size > 1 && J > 1) utils::combn(cand, 2, simplify = FALSE))
        LB <- sapply(models, function(S) {
          f <- .niche_ols(z$Nt, Yt, z$h0, S)
          c1 <- im[z$npp[S[1]]]; f$v11 <- f$v11 * c1
          if (length(S) > 1) { c2 <- im[z$npp[S[2]]]; f$v22 <- f$v22 * c2; f$v12 <- f$v12 * sqrt(c1 * c2) }
          lb <- if (length(S) == 1) .labf(f, tau(S)) else .labf(f, tau(S[1]), tau(S[2]))
          lb[!apply(okk[S, , drop = FALSE], 2, all)] <- -Inf
          lb
        })
        if (is.null(dim(LB))) LB <- matrix(LB, nrow = length(hitg))
        sz <- lengths(models); pi1 <- 1 / J
        lp <- sweep(LB, 2, sz * log(pi1) + (J - sz) * log(1 - pi1), "+")
        post <- exp(lp - apply(lp, 1, max)); post <- post / rowSums(post)
        single <- which(sz == 1)
        al <- exp(LB[, single, drop = FALSE] - apply(LB[, single, drop = FALSE], 1, max)); al <- al / rowSums(al)
        for (ci in seq_along(cand)) {
          j <- cand[ci]; inm <- vapply(models, function(S) j %in% S, TRUE); pp <- z$npp[j]
          att[[length(att) + 1]] <- data.table::data.table(receiver = t, sender = z$nsend[j], gene = genes[hitg],
            z_marg = zM[pp, hitg], q_marg = qM[pp, hitg], z_cond = zC[pp, hitg],
            pip = rowSums(post[, inm, drop = FALSE]), alpha = al[, ci], p_two = rowSums(post[, sz == 2, drop = FALSE]))
        }
      }
    }
    Att <- data.table::rbindlist(att)
    if (nrow(Att)) {
      data.table::setorder(Att, receiver, gene, -alpha)
      Att[, `:=`(in_cs = (cumsum(alpha) - alpha) < cs_level, cs_size = sum((cumsum(alpha) - alpha) < cs_level)), by = list(receiver, gene)]
    }
    tm["pass3"] <- proc.time()[[3]] - t0 - sum(tm); .say(verbose, t0, "pass 3 done: %d genes re-read", length(hit_any))
  }
  type_mean <- tsum / pmax(nt, 1)                                     # mean composition (per 100) by type
  list(genes = genes, K = K, pairs = pairs, fitted = fitted, a = a, type_mean = type_mean, det = det_t,
       scale_fits = fits, omni = omni, omni_q = omni_q, screen = screen, cond = cond, attribution = Att,
       calibration = cal, timing = tm, n_dropped = n_dropped)
}

# ---------------------------------------------------------------- block functions (caller or daemons)
.edgeicg_readX <- function(st, gi) { X <- st$pv$read(gi); X / pmax(st$lib, 1) * st$scale }

.edgeicg_p1 <- function(gi, st) {
  X <- .edgeicg_readX(st, gi); B <- length(gi)
  a <- matrix(NA_real_, st$K, B); b1 <- se1 <- matrix(NA_real_, st$P_n, B)
  for (t in st$fitted) {
    z <- st$dz[[t]]; Y <- X[z$idx, , drop = FALSE]
    coef <- z$Ginv %*% crossprod(z$Aw, Y); R <- Y - z$A %*% coef
    se <- sqrt(crossprod(z$M^2, (R / (1 - z$h))^2))
    a[t, ] <- coef[1, ]
    for (j in seq_along(z$others)) {
      c1 <- z$ci[[j]][1]; if (is.na(c1)) next
      b1[z$pp[j], ] <- coef[c1, ]; se1[z$pp[j], ] <- se[c1, ]
    }
  }
  ut <- sort(unique(st$types))
  detn <- tsum <- matrix(0, st$K, B)
  detn[ut, ] <- rowsum((X > 0) * 1, st$types); tsum[ut, ] <- rowsum(X, st$types)
  list(a = a, b1 = b1, se1 = se1, mbar = as.matrix(st$E %*% X), detn = detn, tsum = tsum)
}

.edgeicg_p2 <- function(gi, st) {
  X <- .edgeicg_readX(st, gi); B <- length(gi); P_n <- st$P_n
  Yres <- X - as.matrix(st$Wl %*% X) - st$cmu * st$a0[st$types, gi, drop = FALSE]
  em <- function() matrix(NA_real_, P_n, B)
  g <- se <- lapply(st$tested, function(s) em()); omni <- bM <- seM <- Sm <- Sc <- zc0 <- em()
  for (t in st$fitted) {
    z <- st$dz[[t]]; Y <- Yres[z$idx, , drop = FALSE]
    coef <- z$Ginv %*% crossprod(z$Aw, Y); R2 <- ((Y - z$A %*% coef) / (1 - z$h))^2
    for (j in seq_along(z$others)) {
      pp <- z$pp[j]; cis <- z$ci[[j]]
      for (si in seq_along(st$tested)) if (!is.na(cis[si])) {
        g[[si]][pp, ] <- coef[cis[si], ]; se[[si]][pp, ] <- sqrt(colSums(z$M[, cis[si]]^2 * R2)) }
      c1 <- cis[1]; c2 <- cis[2]
      if (!is.na(c1) && !is.na(c2)) {
        v11 <- colSums(z$M[, c1]^2 * R2); v22 <- colSums(z$M[, c2]^2 * R2); v12 <- colSums(z$M[, c1] * z$M[, c2] * R2)
        bb1 <- coef[c1, ]; bb2 <- coef[c2, ]
        omni[pp, ] <- (bb1^2 * v22 - 2 * bb1 * bb2 * v12 + bb2^2 * v11) / pmax(v11 * v22 - v12^2, EPS)
      }
    }
    if (!length(z$nsend)) next
    Yt <- qr.resid(z$q, Y * z$sw); Y2 <- Yt^2
    for (j in seq_along(z$nsend)) {
      pp <- z$npp[j]
      f <- .niche_ols(z$Nt, Yt, z$h0, j); bM[pp, ] <- f$b[1, ]; seM[pp, ] <- sqrt(f$v11)
      if (st$calibrate) {
        zm <- .shift_z(z$Zs[[j]], Yt, Y2, z$h0); Sm[pp, ] <- colMeans(zm[-1, , drop = FALSE]^2)
        if (st$calib_cond) {
          Yo <- qr.resid(z$qo[[j]], Yt); zc <- .shift_z(z$ZsC[[j]], Yo, Yo^2, z$hc[[j]])
          Sc[pp, ] <- colMeans(zc[-1, , drop = FALSE]^2); zc0[pp, ] <- zc[1, ]
        }
      }
    }
  }
  list(g = g, se = se, omni = omni, bM = bM, seM = seM, Sm = Sm, Sc = Sc, zc0 = zc0)
}

# Sequential (workers = 0) or mirai daemons. The state goes to each daemon once (as a file); $set() adds to it.
# Daemons load the installed edgeicg namespace, so the mapped block functions are sent by reference.
.edgeicg_runner <- function(workers, st) {
  if (workers < 1) {
    env <- new.env(); env$st <- st
    return(list(map = function(x, f) lapply(x, function(gi) f(gi, env$st)),
                set = function(more) env$st[names(more)] <- more, close = function() invisible()))
  }
  if (!requireNamespace("mirai", quietly = TRUE)) stop("workers > 0 needs the mirai package")
  mirai::daemons(workers)
  closed <- FALSE
  push <- function(obj, name) {
    f <- tempfile(fileext = ".rds"); saveRDS(obj, f)
    mirai::everywhere({ x <- readRDS(f); if (nm == "st") .edgeicg_st <<- x else .edgeicg_st[names(x)] <<- x; TRUE },
                      f = f, nm = name)[]
    unlink(f)
  }
  mirai::everywhere({ Sys.setenv(OMP_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
    loadNamespace("edgeicg"); TRUE })[]
  push(st, "st")
  call1 <- function(gi, f) f(gi, .edgeicg_st); environment(call1) <- globalenv()
  list(map = function(x, f) {
         r <- mirai::mirai_map(x, call1, .args = list(f = f))[]
         bad <- vapply(r, mirai::is_mirai_error, TRUE)
         if (any(bad)) stop("daemon error: ", as.character(r[[which(bad)[1]]]))
         r },
       set = function(more) push(more, "more"),
       close = function() if (!closed) { mirai::daemons(0); closed <<- TRUE })
}
