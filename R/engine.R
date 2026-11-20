# Engine. Model per receiver type t, gene g, cell i (composition scale y = counts / library * 100):
#
#   y_ig = a_tg + sum_k [ C_ik g^c_tkg + P_ik g^p_tkg + N_ik g^n_tkg ] + S_ig + e_ig
#   S_ig = sum_{j touching i, type(j) = k != t} c_ij lambda_tk x_jg + c^mu_i a_tg        (leaked transcripts)
#
# lambda_tk, mu_tk: Huber regression across genes of the stage-1 contact coefficients on the touching-sender
# profile and the receiver profile. S is removed per cell and the model refitted (WLS, weights library / mean
# library, HC3). Niche: marginal screen per sender (Frisch-Waugh on the contact + paracrine design), calibrated
# by a toroidal-shift null (pair-level variance inflation), then exact attribution over 1- and 2-sender models.
#
# Passes (sweeps.R): 1. stage-1 fit -> a, contact b1 / se1, touching-sender profiles, detection, type means;
# 2. leak removal, stage-2 fit (effects, HC3, omnibus), niche screen, shift-null sums; 3. attribution sums on
# screen-hit genes. Each pass is a sweep of (gene block x band of cell windows) tasks returning sums over cells.
# Gene-independent quantities (Gram inverse, leverage, projections, shifted exposures) are computed once per
# receiver type in setup.

.edgeicg_engine <- function(pv, coords, types, lib, ex, spill_edges, block_genes, window_cells, min_cells, min_det,
                            min_exposed, exposed_at, huber, scale, calibrate, B, min_shift, seed, calib_cond, attribute,
                            q_thr, max_size, prior_sd, cs_level, workers, xp, drop_degenerate, verbose) {
  t0 <- proc.time()[[3]]; tm <- c()
  S <- .edgeicg_setup(pv, coords, types, lib, ex, spill_edges, block_genes, window_cells, min_cells, scale, calibrate,
                      B, min_shift, seed, calib_cond, drop_degenerate, verbose, t0)
  list2env(S, environment())                                        # genes, G, n, K, tested, blocks, win, st, ...
  rm(S)
  st$xp <- xp$xp; st$mp <- xp$mp
  bands <- .edgeicg_bands(nrow(win$w), workers)
  if (!st$single) {                                                 # provider in window order (store: bake)
    st$pv <- st$pv$order(win$ord); on.exit(st$pv$cleanup(), add = TRUE)
    .say(verbose, t0, "%d windows of ~%d cells in %d bands (counts %s)", nrow(win$w), ceiling(n / nrow(win$w)),
         length(bands), st$pv$kind)
  }
  run <- .edgeicg_runner(workers, st)
  on.exit(run$close(), add = TRUE)
  do_sweep <- function(name, groups, extra = NULL) .edgeicg_sweep(run, st, win, bands, name, groups, extra, workers)
  per_type <- function(R, f) lapply(R, function(r) lapply(seq_len(K), function(t) if (!is.null(r$t[[t]])) f(dz[[t]], r$t[[t]])))

  # pass 1
  if (st$single) R1 <- do_sweep("p1", blocks) else {
    A1 <- do_sweep("p1a", blocks); F1 <- per_type(A1, .fin1)
    B1 <- do_sweep("p1b", blocks, F1)
    R1 <- Map(function(a, f, b) { a$t <- lapply(seq_len(K), function(t) if (!is.null(f[[t]])) c(f[[t]], b$t[[t]])); a }, A1, F1, B1)
    rm(A1, F1, B1)
  }
  a <- matrix(NA_real_, K, G); b1 <- se1 <- mbar <- matrix(NA_real_, P_n, G); detn <- tsum <- matrix(0, K, G)
  for (k in seq_along(blocks)) {
    gi <- blocks[[k]]; r <- R1[[k]]
    mbar[, gi] <- r$mbar; detn[, gi] <- r$detn; tsum[, gi] <- r$tsum
    for (t in fitted) {
      z <- dz[[t]]; f <- r$t[[t]]; a[t, gi] <- f$coef[1, ]
      for (j in seq_along(z$others)) {
        c1 <- z$ci[[j]][1]; if (is.na(c1)) next
        b1[z$pp[j], gi] <- f$coef[c1, ]; se1[z$pp[j], gi] <- sqrt(f$se2[c1, ])
      }
    }
  }
  rm(R1); mbar[sw_p == 0, ] <- NA
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
  a0 <- a; a0[is.na(a0)] <- 0
  st$Wlt <- vector("list", K); for (t in fitted) st$Wlt[[t]] <- Wl[dz[[t]]$idx, , drop = FALSE]   # leak rows per type
  st$cmu <- .rsum(wn * mu_e, ct$recv, n); st$a0 <- a0
  run$set(if (st$single || workers < 1) st[c("Wlt", "cmu", "a0")] else st["a0"])   # daemons with windows: leak rows travel in the contexts
  rm(Wl)

  # pass 2
  if (st$single) R2 <- do_sweep("p2", blocks) else {
    A2 <- do_sweep("p2a", blocks); F2 <- per_type(A2, .fin2); rm(A2)
    B2 <- do_sweep("p2b", blocks, F2)
    R2 <- Map(function(f, b) list(t = lapply(seq_len(K), function(t)
      if (!is.null(f[[t]])) c(f[[t]][c("QbY", "NtYt", "coef", "B")], b$t[[t]]))), F2, B2)
    rm(F2, B2)
  }
  est <- lapply(tested, function(s) list(g = matrix(NA_real_, P_n, G), se = matrix(NA_real_, P_n, G))); names(est) <- tested
  omni <- bM <- seM <- Sm <- Sc <- zc0 <- matrix(NA_real_, P_n, G)
  QbY_all <- NtYt_all <- vector("list", K)
  for (t in fitted) {
    QbY_all[[t]] <- matrix(NA_real_, ncol(dz[[t]]$Qb), G)
    if (length(dz[[t]]$nsend)) NtYt_all[[t]] <- matrix(NA_real_, length(dz[[t]]$nsend), G)
  }
  for (k in seq_along(blocks)) {
    gi <- blocks[[k]]; r <- R2[[k]]
    for (t in fitted) {
      z <- dz[[t]]; f <- r$t[[t]]
      QbY_all[[t]][, gi] <- f$QbY
      for (j in seq_along(z$others)) {
        pp <- z$pp[j]; cis <- z$ci[[j]]
        for (si in seq_along(tested)) if (!is.na(cis[si])) {
          est[[si]]$g[pp, gi] <- f$coef[cis[si], ]; est[[si]]$se[pp, gi] <- sqrt(f$se2[cis[si], ]) }
        jj <- match(j, z$omj)
        if (!is.na(jj)) {
          c1 <- cis[1]; c2 <- cis[2]; v11 <- f$se2[c1, ]; v22 <- f$se2[c2, ]; v12 <- f$om12[jj, ]
          bb1 <- f$coef[c1, ]; bb2 <- f$coef[c2, ]
          omni[pp, gi] <- (bb1^2 * v22 - 2 * bb1 * bb2 * v12 + bb2^2 * v11) / pmax(v11 * v22 - v12^2, EPS)
        }
      }
      if (!length(z$nsend)) next
      NtYt_all[[t]][, gi] <- f$NtYt
      for (j in seq_along(z$nsend)) {
        pp <- z$npp[j]; bM[pp, gi] <- f$B[j, ]; seM[pp, gi] <- sqrt(f$v11n[j, ])
        if (calibrate) {
          zm <- .shift_fin(f$sh[[j]], z$zz[[j]], z$wz2[[j]]); Sm[pp, gi] <- colMeans(zm[-1, , drop = FALSE]^2)
          if (calib_cond) {
            zc <- .shift_fin(f$shc[[j]], z$zzC[[j]], z$wz2C[[j]])
            Sc[pp, gi] <- colMeans(zc[-1, , drop = FALSE]^2); zc0[pp, gi] <- zc[1, ]
          }
        }
      }
    }
  }
  rm(R2)
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

  # pass 3: attribution on screen-hit genes. Models and their gene-independent parts are planned here; the sweep
  # returns the HC3 sums; Bayes factors and posteriors are computed here.
  Att <- NULL
  if (attribute) {
    hit_any <- which(colSums(qM < q_thr, na.rm = TRUE) > 0)
    hblocks <- unname(split(hit_any, ceiling(seq_along(hit_any) / block_genes)))
    plans <- lapply(hblocks, function(gi) lapply(seq_len(K), function(t) {
      z <- dz[[t]]; if (is.null(z) || !length(z$nsend)) return(NULL)
      rows <- which(pairs$receiver == t)
      hl <- which(colSums(qM[rows, gi, drop = FALSE] < q_thr, na.rm = TRUE) > 0); if (!length(hl)) return(NULL)
      hitg <- gi[hl]; okk <- okn[z$npp, hitg, drop = FALSE]
      cand <- which(rowSums(okk) > 0); J <- length(cand); if (!J) return(NULL)
      models <- c(lapply(cand, function(j) j), if (max_size > 1 && J > 1) utils::combn(cand, 2, simplify = FALSE))
      Gi <- lapply(models, function(S) solve(z$NtNt[S, S, drop = FALSE]))
      list(hl = hl, models = models, Gi = Gi, QbY = QbY_all[[t]][, hitg, drop = FALSE],
           Bm = Map(function(S, g) g %*% NtYt_all[[t]][S, hitg, drop = FALSE], models, Gi))
    }))
    R3 <- if (length(hblocks)) do_sweep("p3", hblocks, plans) else list()
    att <- list()
    for (k in seq_along(hblocks)) {
      gi <- hblocks[[k]]
      for (t in fitted) {
        pl <- plans[[k]][[t]]; if (is.null(pl)) next
        z <- dz[[t]]; r <- R3[[k]]$t[[t]]; hitg <- gi[pl$hl]
        okk <- okn[z$npp, hitg, drop = FALSE]; cand <- which(rowSums(okk) > 0); J <- length(cand); models <- pl$models
        sy <- sqrt(r$yt2 / length(z$idx)); sx <- sqrt(colMeans(z$Nt^2))
        tau <- function(j) outer(prior_sd, sy / sx[j])
        LB <- sapply(seq_along(models), function(mi) {
          S <- models[[mi]]; m <- r$m[[mi]]
          f <- list(b = pl$Bm[[mi]], v11 = as.vector(m$v11))
          c1 <- im[z$npp[S[1]]]; f$v11 <- f$v11 * c1
          if (length(S) > 1) { c2 <- im[z$npp[S[2]]]; f$v22 <- as.vector(m$v22) * c2; f$v12 <- as.vector(m$v12) * sqrt(c1 * c2) }
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
  run$close()
  type_mean <- tsum / pmax(nt, 1)                                     # mean composition (per 100) by type
  list(genes = genes, K = K, pairs = pairs, fitted = fitted, a = a, type_mean = type_mean, det = det_t,
       scale_fits = fits, omni = omni, omni_q = omni_q, screen = screen, cond = cond, attribution = Att,
       calibration = cal, timing = tm, n_dropped = n_dropped, n_windows = nrow(win$w))
}

# One sweep: tasks, run, sums reduced per gene group. With one window each task is one gene group (parallel over
# genes); with several, each task is a band of windows x all groups (parallel over cells), and partial sums of
# the bands are added up.
.edgeicg_sweep <- function(run, st, win, bands, name, groups, extra, workers) {
  ship <- workers > 0 && is.function(st$pv$slice)                     # in-memory counts travel with the task
  halo <- .sweep_pass(name) != "p1"
  wspec <- function(w, ug, ctx) {
    wi <- c(list(w = w), as.list(win$w[w, ]))
    if (ctx) wi$ctx <- .ctx_build(st, wi, name)
    if (ship) wi$blk <- st$pv$slice(if (halo) wi$rlo else wi$lo, if (halo) wi$rhi else wi$hi, ug)
    wi
  }
  if (st$single) {
    tasks <- lapply(seq_along(groups), function(k)
      list(wins = list(wspec(1L, groups[[k]], FALSE)), groups = groups[k], extra = extra[k], gid = k))
  } else {
    ug <- sort(unique(unlist(groups)))
    tasks <- lapply(bands, function(b)
      list(wins = lapply(b, function(w) wspec(w, ug, workers > 0)), groups = groups, extra = extra, gid = seq_along(groups)))
  }
  res <- run$run(tasks, name); rm(tasks)
  out <- vector("list", length(groups))
  for (i in seq_along(res)) for (k in seq_along(res[[i]])) {
    g <- if (st$single) i else k
    out[[g]] <- .acc_add(out[[g]], res[[i]][[k]])
  }
  out
}

# Setup shared by every pass: shifts, window plan, per-receiver-type designs and their gene-independent quantities,
# the touching-sender averaging matrix E. Returns the variables the engine uses, and st (the state the tasks need).
.edgeicg_setup <- function(pv, coords, types, lib, ex, spill_edges, block_genes, window_cells, min_cells, scale,
                           calibrate, B, min_shift, seed, calib_cond, drop_degenerate, verbose, t0 = proc.time()[[3]]) {
  genes <- pv$genes; G <- length(genes); n <- length(types); K <- max(types)
  tested <- c("contact", "paracrine", "niche")
  blocks <- unname(split(seq_len(G), ceiling(seq_len(G) / block_genes)))
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
    # column of each scale per sender (NA = dropped); senders with both contact and paracrine (omnibus)
    z$ci <- lapply(z$others, function(kk) vapply(tested, function(s) {
      w <- which(z$cols$scale == s & z$cols$k == kk); if (length(w)) 1L + w else NA_integer_ }, 1L))
    z$omj <- which(vapply(z$ci, function(v) !is.na(v[1]) && !is.na(v[2]), TRUE))
    z$omc <- do.call(rbind, lapply(z$ci[z$omj], function(v) v[1:2]))
    nc <- which(z$cols$scale == "niche"); bc <- which(z$cols$scale != "niche")
    z$sw <- sqrt(cw[idx]); q <- qr(cbind(1, z$D[, bc, drop = FALSE]) * z$sw)
    Qf <- qr.Q(q); z$Qb <- Qf[, seq_len(q$rank), drop = FALSE]; z$h0 <- rowSums(Qf^2); rm(Qf)
    z$Nt <- qr.resid(q, z$D[, nc, drop = FALSE] * z$sw)
    z$NtN <- colSums(z$Nt^2); z$NtNt <- crossprod(z$Nt); z$NtQb <- crossprod(z$Nt, z$Qb)
    z$nsend <- z$cols$k[nc]                                         # sender of each Nt column
    z$jn <- match(z$others, z$nsend)                                 # Nt column per sender (NA = none)
    z$npp <- z$pp[match(z$nsend, z$others)]                          # pair row per Nt column
    if (calibrate && length(nc)) {
      z$Zs <- z$zz <- z$wz2 <- z$Qo <- z$hc <- z$QoZs <- z$QoQb <- z$zzC <- z$wz2C <- vector("list", length(nc))
      gsum <- function(Z, hb) {                                      # gene-independent shift sums
        zz <- colSums(Z^2); h <- pmin(hb + sweep(Z^2, 2, zz, "/"), 0.99)
        list(zz = zz, wz2 = colSums(sweep(Z, 2, zz, "/")^2 / (1 - h)^2 * Z^2))
      }
      for (j in seq_along(nc)) {
        kk <- z$nsend[j]
        raw <- vapply(seq_len(nrow(sh0)), function(b) {
          bk <- fl$F[[kk]][cbind((fl$gx[idx] - 1L + sh0[b, 1]) %% fl$nx + 1L, (fl$gy[idx] - 1L + sh0[b, 2]) %% fl$ny + 1L)]
          bk / pmax(fl$tot[idx] - fl$Bc[idx, kk] + bk, EPS) }, numeric(length(idx)))
        z$Zs[[j]] <- qr.resid(q, raw * z$sw)
        g0 <- gsum(z$Zs[[j]], z$h0); z$zz[[j]] <- g0$zz; z$wz2[[j]] <- g0$wz2
        if (calib_cond) {
          qo <- qr(z$Nt[, -j, drop = FALSE]); Qo <- qr.Q(qo)
          z$hc[[j]] <- z$h0 + rowSums(Qo^2); z$Qo[[j]] <- Qo[, seq_len(qo$rank), drop = FALSE]
          z$QoZs[[j]] <- crossprod(z$Qo[[j]], z$Zs[[j]]); z$QoQb[[j]] <- crossprod(z$Qo[[j]], z$Qb)
          gc_ <- gsum(z$Zs[[j]] - z$Qo[[j]] %*% z$QoZs[[j]], z$hc[[j]]); z$zzC[[j]] <- gc_$zz; z$wz2C[[j]] <- gc_$wz2
        }
      }
    }
    z$D <- NULL
    dz[[t]] <- z
  }
  ep <- match(paste(tr, ts), paste(pairs$receiver, pairs$sender)); ue <- het & !is.na(ep)
  sw_p <- .rsum(wn[ue], ep[ue], P_n)
  E <- Matrix::sparseMatrix(i = ep[ue], j = ct$send[ue], x = wn[ue] / sw_p[ep[ue]], dims = c(P_n, n))
  win <- .edgeicg_windows(coords, ct, n, window_cells)
  .say(verbose, t0, "setup: %d genes in %d blocks, %d receiver types, %d pairs, %d degenerate columns dropped",
       G, length(blocks), length(fitted), P_n, n_dropped)

  st <- list(pv = pv, lib = lib, scale = scale, dz = dz, fitted = fitted, types = types, E = E, K = K, P_n = P_n,
             tested = tested, calibrate = calibrate, calib_cond = calib_cond, single = nrow(win$w) == 1L,
             ord = win$ord, wtypes = .edgeicg_window_types(win, types, dz, fitted))
  sh <- if (calibrate) sh else NULL
  list(genes = genes, G = G, n = n, K = K, tested = tested, blocks = blocks, sh = sh, ct = ct, wn = wn, tr = tr,
       ts = ts, het = het, smat = smat, dz = dz, pairs = pairs, fitted = fitted, P_n = P_n, n_dropped = n_dropped,
       sw_p = sw_p, E = E, win = win, st = st)
}
