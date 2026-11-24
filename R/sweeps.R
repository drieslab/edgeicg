# Sweep kernels. Every statistic of the model is a sum over cells, so each pass is computed as partial sums over
# (cell window x gene block) tasks, added up, then finalised. A pass has an "a" part (sums that need no fitted
# coefficients) and a "b" part (residual sums that do). With a single window the task holds every cell of a type,
# so it finalises "a" itself and runs "b" straight away ("ab"); with several windows "a" sums are reduced and
# finalised in the main process first, and "b" runs as a second sweep.
#
# Dense products go through st$xp (t(A) %*% B) and st$mp (A %*% B): base R or the edgeicgrs kernels.
# Kernel arguments: st (state), ctx (window context: per receiver type, the row slices of the type's design
# arrays for the window's cells, plus global quantities), X (cells of the read range x genes, sparse, composition
# scale), gi (gene indices), ex (per type: finalised "a" quantities, or the attribution plan).

.row_fields <- list(
  p1  = c("A", "Aw", "M", "h"),
  p2a = c("Aw", "sw", "Qb", "Nt", "Qo"),
  p2b = c("A", "M", "h", "sw", "Qb", "h0", "Nt", "Zs", "Qo", "hc", "sg", "sXC", "sNC"),
  p3  = c("sw", "Qb", "h0", "Nt"))
.glob_fields <- c("Ginv", "NtN", "NtNt", "NtQb", "zz", "wz2", "QoZs", "zzC", "wz2C", "QoQb", "ci", "pp",
                  "others", "nsend", "npp", "omc", "omj")

.sweep_pass <- function(sweep) substr(sweep, 1L, 2L)
.sweep_part <- function(sweep) { p <- substring(sweep, 3L); if (p == "") "ab" else p }
.sweep_fields <- function(sweep) {
  ps <- .sweep_pass(sweep); pt <- .sweep_part(sweep)
  if (ps == "p1") return(.row_fields$p1)
  if (ps == "p3") return(.row_fields$p3)
  switch(pt, a = .row_fields$p2a, b = .row_fields$p2b, ab = union(.row_fields$p2a, .row_fields$p2b))
}

# Window context. With one window the type arrays are used as they are (no copies).
.ctx_build <- function(st, wi, sweep) {
  halo <- .sweep_pass(sweep) %in% c("p2", "p3")
  r0 <- if (halo) wi$rlo else wi$lo; r1 <- if (halo) wi$rhi else wi$hi
  read_cells <- if (is.null(st$ord)) r0:r1 else st$ord[r0:r1]
  core <- if (is.null(st$ord)) wi$lo:wi$hi else st$ord[wi$lo:wi$hi]
  if (.sweep_pass(sweep) == "p0") return(list(xc = (wi$lo:wi$hi) - r0 + 1L, ctype = st$types[core]))
  fields <- .sweep_fields(sweep)
  types <- lapply(seq_len(st$K), function(t) {
    z <- st$dz[[t]]; wt <- st$wtypes[[wi$w]][[t]]
    if (is.null(z) || is.null(wt) || !length(wt$loc)) return(NULL)
    out <- z[intersect(.glob_fields, names(z))]
    sl <- function(x) {
      if (st$single) return(x)
      if (is.list(x)) return(lapply(x, sl))
      if (is.matrix(x)) x[wt$loc, , drop = FALSE] else x[wt$loc]
    }
    for (f in intersect(fields, names(z))) out[[f]] <- sl(z[[f]])
    out$xr <- wt$p - r0 + 1L                                         # rows of the type's cells in X
    if (halo) {
      out$Wl <- if (st$single) st$Wlt[[t]] else st$Wlt[[t]][wt$loc, read_cells, drop = FALSE]
      out$cmu <- st$cmu[z$idx[wt$loc]]
    }
    out
  })
  ctx <- list(types = types, xc = (wi$lo:wi$hi) - r0 + 1L)
  if (.sweep_pass(sweep) == "p1" && .sweep_part(sweep) != "b") {
    ctx$E <- if (st$single) st$E else st$E[, core, drop = FALSE]
    ctx$ctype <- st$types[core]
  }
  ctx
}

# raw counts of the read range for genes ug -> cells x genes, sparse, composition scale
.edgeicg_read <- function(st, wi, ug, halo) {
  r0 <- if (halo) wi$rlo else wi$lo; r1 <- if (halo) wi$rhi else wi$hi
  X <- if (!is.null(wi$blk)) Matrix::t(wi$blk) else st$pv$rows(r0, r1, ug - (if (is.null(st$goff)) 0L else st$goff))
  X <- methods::as(methods::as(methods::as(X, "CsparseMatrix"), "generalMatrix"), "dMatrix")
  lib <- if (is.null(st$ord)) st$lib[r0:r1] else st$lib[st$ord[r0:r1]]
  X@x <- X@x * (st$scale / pmax(lib, 1))[X@i + 1L]
  X
}

# One task: windows of a band x gene groups. Each window is read once for the union of the groups' genes.
.edgeicg_task <- function(task, st, sweep) {
  ug <- sort(unique(unlist(task$groups)))
  halo <- .sweep_pass(sweep) %in% c("p2", "p3")
  out <- vector("list", length(task$groups))
  for (wi in task$wins) {
    ctx <- if (is.null(wi$ctx)) .ctx_build(st, wi, sweep) else wi$ctx
    Xr <- .edgeicg_read(st, wi, ug, halo)
    for (k in seq_along(task$groups)) {
      gi <- task$groups[[k]]
      X <- if (length(gi) == length(ug)) Xr else Xr[, match(gi, ug), drop = FALSE]
      r <- switch(.sweep_pass(sweep),
                  p0 = .k_p0(st, ctx, X),
                  p1 = .k_p1(st, ctx, X, gi, .sweep_part(sweep), task$extra[[k]]),
                  p2 = .k_p2(st, ctx, X, gi, .sweep_part(sweep), task$extra[[k]]),
                  p3 = .k_p3(st, ctx, X, gi, task$extra[[k]]))
      out[[k]] <- if (is.null(out[[k]])) r else .acc_add(out[[k]], r)
    }
    rm(Xr, ctx); .edgeicg_gc(st)
  }
  out
}

# Elementwise sum of two partial results (nested lists of numeric arrays; NULL entries skipped)
.acc_add <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  if (is.list(a)) { for (k in seq_along(b)) a[k] <- list(.acc_add(a[[k]], b[[k]])); return(a) }
  a + b
}

# ---------------------------------------------------------------- pass 0: detection per receiver type (prefilter)
.k_p0 <- function(st, ctx, X) {
  Xd <- X[ctx$xc, , drop = FALSE]; Xd@x <- as.numeric(Xd@x > 0)
  Ind <- Matrix::sparseMatrix(i = ctx$ctype, j = seq_along(ctx$ctype), x = 1, dims = c(st$K, length(ctx$ctype)))
  list(detn = as.matrix(Ind %*% Xd))
}

# ---------------------------------------------------------------- pass 1: stage-1 fit
.fin1 <- function(z, S) list(coef = z$Ginv %*% S$A)

.k_p1 <- function(st, ctx, X, gi, part, ex) {
  res <- list(t = vector("list", st$K))
  if (part != "b") {
    Xc <- X[ctx$xc, , drop = FALSE]
    Ind <- Matrix::sparseMatrix(i = ctx$ctype, j = seq_along(ctx$ctype), x = 1, dims = c(st$K, length(ctx$ctype)))
    Xd <- Xc; Xd@x <- as.numeric(Xd@x > 0)
    res$mbar <- as.matrix(ctx$E %*% Xc); res$detn <- as.matrix(Ind %*% Xd); res$tsum <- as.matrix(Ind %*% Xc)
  }
  for (t in st$fitted) {
    z <- ctx$types[[t]]; if (is.null(z)) next
    Y <- as.matrix(X[z$xr, , drop = FALSE])
    if (part == "a") { res$t[[t]] <- list(A = st$xp(z$Aw, Y)); next }
    f <- if (part == "ab") .fin1(z, list(A = st$xp(z$Aw, Y))) else ex[[t]]
    R <- (Y - st$mp(z$A, f$coef)) / (1 - z$h)
    res$t[[t]] <- c(if (part == "ab") f, list(se2 = st$xp(z$M^2, R^2)))
  }
  res
}

# ---------------------------------------------------------------- pass 2: leak removal, stage 2, niche screen, shift null
.edgeicg_yres <- function(st, z, X, t, gi) {
  as.matrix(X[z$xr, , drop = FALSE]) - as.matrix(z$Wl %*% X) - outer(z$cmu, st$a0[t, gi])
}

.fin2 <- function(z, S) {
  f <- list(coef = z$Ginv %*% S$A, QbY = S$Qb)
  if (length(z$nsend)) {
    f$NtYt <- S$Nt - z$NtQb %*% S$Qb; f$B <- f$NtYt / z$NtN
    if (!is.null(S$Qo)) f$QoYt <- Map(function(s, QQ) s - QQ %*% S$Qb, S$Qo, z$QoQb)
  }
  f
}

.k_p2 <- function(st, ctx, X, gi, part, ex) {
  res <- list(t = vector("list", st$K))
  for (t in st$fitted) {
    z <- ctx$types[[t]]; if (is.null(z)) next
    Y <- .edgeicg_yres(st, z, X, t, gi)
    if (part != "b") {
      Ysw <- Y * z$sw
      S <- list(A = st$xp(z$Aw, Y), Qb = st$xp(z$Qb, Ysw))
      if (length(z$nsend)) {
        S$Nt <- st$xp(z$Nt, Ysw)
        if (st$calibrate && st$calib_cond) S$Qo <- lapply(z$Qo, function(Q) st$xp(Q, Ysw))
      }
      if (part == "a") { res$t[[t]] <- S; next }
      f <- .fin2(z, S)
    } else f <- ex[[t]]
    res$t[[t]] <- c(if (part == "ab") f[c("QbY", "NtYt", "coef", "B")], .k_p2b(st, z, Y, f))
  }
  res
}

# HC3 sums of the shifted-exposure tests: for each shift s, b_s = Z_s'Y / zz_s and its HC3 variance
# v_s = W_s'Y^2 - 2 b_s (W_s Z_s)'Y + b_s^2 sum(W_s Z_s^2); the three products are sums over cells.
.shift_sums <- function(st, Z, hb, zz, Y, Y2) {
  h <- pmin(hb + sweep(Z^2, 2, zz, "/"), 0.99)
  W <- sweep(Z, 2, zz, "/")^2 / (1 - h)^2
  list(zy = st$xp(Z, Y), wy2 = st$xp(W, Y2), wzy = st$xp(W * Z, Y))
}
.shift_fin <- function(s, zz, wz2) {
  b <- s$zy / zz
  b / sqrt(pmax(s$wy2 - 2 * b * s$wzy + b^2 * wz2, EPS))
}

.k_p2b <- function(st, z, Y, f) {
  Rr <- Y - st$mp(z$A, f$coef)
  R2 <- (Rr / (1 - z$h))^2
  out <- list(se2 = st$xp(z$M^2, R2))
  if (length(z$omj)) out$om12 <- st$xp(z$M[, z$omc[, 1], drop = FALSE] * z$M[, z$omc[, 2], drop = FALSE], R2)
  Yt <- if (length(z$nsend)) Y * z$sw - st$mp(z$Qb, f$QbY) else NULL
  if (length(z$sg)) {                                                # region / local strata (strata.R)
    rs <- Rr * z$sw
    out$str <- lapply(seq_along(z$sg), function(k) list(
      cp = .strata_sums(z$sg[[k]], z$sXC[[k]], rs, st$strata_D[k]),
      nt = if (!is.null(Yt)) .strata_sums(z$sg[[k]], z$sNC[[k]], Yt, st$strata_D[k])))
  }
  rm(Rr)
  if (!length(z$nsend)) return(out)
  nn <- length(z$nsend); v11 <- matrix(0, nn, ncol(Y))
  for (j in seq_len(nn)) {
    Zj <- z$Nt[, j]; hj <- pmin(z$h0 + Zj^2 / z$NtN[j], 0.99)
    v11[j, ] <- st$xp((Zj / z$NtN[j])^2, ((Yt - outer(Zj, f$B[j, ])) / (1 - hj))^2)
  }
  out$v11n <- v11
  if (st$calibrate) {
    Y2 <- Yt^2; out$sh <- out$shc <- vector("list", nn)
    for (j in seq_len(nn)) {
      out$sh[[j]] <- .shift_sums(st, z$Zs[[j]], z$h0, z$zz[[j]], Yt, Y2)
      if (st$calib_cond) {
        Yo <- Yt - st$mp(z$Qo[[j]], f$QoYt[[j]])
        ZC <- z$Zs[[j]] - st$mp(z$Qo[[j]], z$QoZs[[j]])
        out$shc[[j]] <- .shift_sums(st, ZC, z$hc[[j]], z$zzC[[j]], Yo, Yo^2)
      }
    }
  }
  out
}

# ---------------------------------------------------------------- pass 3: attribution sums on screen-hit genes
# ex[[t]]: hl (hit columns within gi), QbY, models (sender sets, columns of Nt), Gi and Bm per model.
.k_p3 <- function(st, ctx, X, gi, ex) {
  res <- list(t = vector("list", st$K))
  for (t in st$fitted) {
    pl <- ex[[t]]; z <- ctx$types[[t]]
    if (is.null(pl) || is.null(z)) next
    Xh <- X[, pl$hl, drop = FALSE]
    Yt <- .edgeicg_yres(st, z, Xh, t, gi[pl$hl]) * z$sw - st$mp(z$Qb, pl$QbY)
    out <- list(yt2 = colSums(Yt^2), m = vector("list", length(pl$models)))
    for (mi in seq_along(pl$models)) {
      S <- pl$models[[mi]]; Zs <- z$Nt[, S, drop = FALSE]
      M <- st$mp(Zs, pl$Gi[[mi]]); h <- pmin(z$h0 + rowSums(M * Zs), 0.99)
      R2 <- ((Yt - st$mp(Zs, pl$Bm[[mi]])) / (1 - h))^2
      o <- list(v11 = st$xp(M[, 1]^2, R2))
      if (length(S) > 1) { o$v22 <- st$xp(M[, 2]^2, R2); o$v12 <- st$xp(M[, 1] * M[, 2], R2) }
      out$m[[mi]] <- o
    }
    res$t[[t]] <- out
  }
  res
}
