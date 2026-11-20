# Cell windows. The engine works on sums over cells, so the cells can be visited in windows and the sums added
# up. Windows are contiguous ranges of an internal cell order sorted by y (then x): the touching neighbours that the
# leak correction needs (the halo) then lie in a contiguous range just outside the window, so a window plus its halo
# is one contiguous row range. On a GiottoDisk store baked in this order (see counts_provider_store()), such a range
# prunes parquet row groups (GiottoDisk adr/0008); a scattered cell set would not.
#
# One window (the in-memory default) keeps the original cell order and the engine's single-pass shortcut.

# window plan: ord (internal position -> cell, NULL = identity), pos (cell -> position) and per window the core
# range [lo, hi] and the read range with halo [rlo, rhi]
.edgeicg_windows <- function(coords, ct, n, window_cells) {
  nw <- if (is.finite(window_cells)) max(1L, as.integer(ceiling(n / window_cells))) else 1L
  if (nw == 1L) return(list(ord = NULL, pos = NULL, w = data.frame(lo = 1L, hi = n, rlo = 1L, rhi = n)))
  ord <- order(coords[, 2], coords[, 1]); pos <- integer(n); pos[ord] <- seq_len(n)
  br <- round(seq(0, n, length.out = nw + 1L))
  lo <- as.integer(br[-length(br)] + 1L); hi <- as.integer(br[-1L])
  wr <- findInterval(pos[ct$recv], lo)                              # window of each edge's receiver
  ps <- pos[ct$send]
  rlo <- pmin(lo, vapply(split(ps, factor(wr, levels = seq_len(nw))), function(v) if (length(v)) min(v) else .Machine$integer.max, 1L))
  rhi <- pmax(hi, vapply(split(ps, factor(wr, levels = seq_len(nw))), function(v) if (length(v)) max(v) else 0L, 1L))
  list(ord = ord, pos = pos, w = data.frame(lo = lo, hi = hi, rlo = as.integer(rlo), rhi = as.integer(rhi)))
}

# Contiguous bands of windows, one task per band (at most `workers` bands; a worker keeps its partial sums for the
# whole band, so the number of partial results is the number of bands, not of windows).
.edgeicg_bands <- function(nw, workers) {
  nb <- max(1L, min(nw, workers))
  unname(split(seq_len(nw), ceiling(seq_len(nw) * nb / nw)))
}

# Per window and receiver type: the type's cells in the window core, as positions within the type's arrays (loc)
# and as internal positions (p).
.edgeicg_window_types <- function(win, types, dz, fitted) {
  n <- length(types); typeloc <- integer(n)
  for (t in fitted) typeloc[dz[[t]]$idx] <- seq_along(dz[[t]]$idx)
  lapply(seq_len(nrow(win$w)), function(w) {
    p <- win$w$lo[w]:win$w$hi[w]; cells <- if (is.null(win$ord)) p else win$ord[p]
    tt <- types[cells]
    out <- vector("list", max(types))
    for (t in fitted) { s <- which(tt == t); out[[t]] <- list(loc = typeloc[cells[s]], p = p[s]) }
    out
  })
}
