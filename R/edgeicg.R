#' Spillover-aware interaction-changed genes at three scales
#'
#' Tests, for every receiver type, sender type and gene, whether the receiver's expression changes with its
#' exposure to the sender at contact, paracrine and niche scale, after removing transcripts leaked from
#' touching neighbours. See the package description for the model.
#'
#' @param counts raw counts, genes x cells: a `dgCMatrix` or dense matrix, a GiottoDisk expression store, or a
#'   provider from [counts_provider()] / [counts_provider_store()]. Give all genes (or a broad set) and choose
#'   the genes to test with `genes`: the leak coefficients are fitted across the genes of `counts` (see
#'   `leak_genes`). The library size comes from `lib`.
#' @param genes genes to test (default: all genes of `counts`).
#' @param leak_genes genes the leak coefficients (lambda, mu per receiver-sender pair) are fitted on. `"all"`
#'   (default): all genes of `counts` that pass the prefilter, so the leak, and with it the contact results, do
#'   not change with the genes chosen for testing. This costs one extra pass over all genes (stage 1 only).
#'   `"tested"`: only the genes in `genes`; faster, but the leak then depends on that selection.
#' @param leak optional leak coefficients to reuse instead of fitting them: an `EdgeICG` fit or its `pairs`
#'   table (columns `receiver`, `sender`, `lambda`, `mu`). Pairs without a value get no leak correction.
#' @param coords n x 2 cell centroids (microns), in the provider's cell order.
#' @param types cell-type labels (character or factor), one per cell.
#' @param lib library size per cell over **all** genes (default: column sums of `counts` if it is a matrix).
#' @param spill_edges optional touching-cell edge list (`recv`, `send`, `w`), e.g. from
#'   [contact_from_polygons()]; the spillover graph. Default: the kNN contact graph.
#' @param ex optional exposures from [edgeicg_exposures()] (computed if NULL).
#' @param block_genes genes per block (memory: about cells of a type x block_genes x 8 bytes per copy).
#' @param window_cells cells per window. With more cells than this, cells are visited in spatially ordered windows
#'   (each with the halo of touching neighbours the leak correction needs) and every statistic is accumulated as a
#'   sum over windows, so memory no longer grows with the section. `NULL` (default): one window for in-memory
#'   counts; for a store, one window per worker, at most `getOption("edgeicg.window_cells", 50000)` cells (the
#'   store is then read once per sweep, in contiguous cell ranges). Results equal the one-window run up to
#'   rounding.
#' @param workers number of mirai daemons (0 = sequential). One window: tasks are gene blocks; several windows:
#'   tasks are bands of windows. Results are identical.
#' @param backend `"auto"` (Rust kernels from the edgeicgrs package when installed, else R), `"rust"` or `"r"`.
#'   The Rust kernels are the matrix cross-products of the sweeps, the hotspot of the method.
#' @param threads threads per process for the Rust kernels (default: 1 with workers, all cores without).
#' @param covariates optional numeric matrix (cells x q) of covariates added to every receiver type's design
#'   and adjusted for but not tested, e.g. a zonation score or [spatial_basis()] to absorb large-scale spatial
#'   trends (zonation, sub-state gradients) that otherwise correlate with exposure. Columns without variance
#'   within a receiver type are dropped (with `drop_degenerate`).
#' @param covariate_scope `"all"` (default): covariates enter every test's design; `"niche"`: only the niche
#'   screen's base design (and its calibration). Large-scale spatial trends confound mainly the niche scale; on
#'   the Atera kidney section `spatial_basis()` brought the niche pseudo-null from 4.7% to 0.5-0.7% with either
#'   scope, at about 1.4x the run time (df 4 or 6, either scope; the niche screen and its calibration dominate
#'   the extra cost). `"niche"` leaves the contact and paracrine results exactly as without covariates.
#' @param panel_genes with several windows, genes are processed in panels of at most this many (each a separately
#'   baked store when counts are on disk), which bounds the memory of the partial sums; every count is still read
#'   once per pass. Default 3000.
#' @param calibrate run the toroidal-shift calibration of the niche scale (recommended).
#' @param calib_cond also calibrate the conditional (all-sender) niche z (secondary output).
#' @param attribute attribute niche screen hits to senders.
#' @param B,min_shift number of shifts and minimum shift (microns) of the calibration.
#' @param q_thr FDR threshold for the niche screen hits that are attributed.
#' @param min_cells,min_det,min_exposed,exposed_at testability: receiver types with at least `min_cells`
#'   cells; genes detected in at least `min_det` of receivers; at least `min_exposed` receivers with exposure
#'   above `exposed_at`. `min_cells = NULL` (default) is twice the coefficients of a receiver type's fit,
#'   `2 * (1 + 3 * (K - 1) + q)` for K types and q covariates, at most 60 (44 for 8 types, 60 from 11 types).
#'   When no type has `min_cells` cells, the result is empty, with a message.
#' @param drop_degenerate drop exposure columns without variance within a receiver type (their tests are NA).
#'   FALSE keeps them in the design.
#' @param prefilter drop genes detected in fewer than `min_det` of the cells of every receiver type before the
#'   fit (one cheap detection pass). Such genes have no testable entry; without the prefilter they still enter
#'   the regression of the leak coefficients across genes.
#' @param seed seed of the calibration shifts.
#' @param verbose progress messages.
#' @param region optional region per cell (character or factor; NA = not in a region). The effect of every entry
#'   is then also estimated within each region, with a test of whether it differs between regions; see
#'   [edgeicg_region_results()]. The leak model and the other terms of the model stay global.
#' @param local experimental. TRUE: the same within square patches of `patch_size` microns
#'   ([edgeicg_patches()]), tested with Cochran's Q over patches; see [edgeicg_local_results()]. The Q test is
#'   calibrated but has little power for smooth or half-plane patterns; prefer `region` with labels.
#' @param patch_size patch side for `local` (microns).
#' @param min_cells_region,min_exposed_region a region or patch is used for a receiver type with at least
#'   `min_cells_region` of its cells, and for an entry with at least `min_exposed_region` exposed receivers.
#' @param ... further engine settings: `huber`, `scale`, `max_size`, `prior_sd`, `cs_level`, `k`, `sigma_p`,
#'   `sigma_n`, `bin`.
#' @return an `EdgeICG` object; see [edgeicg_results()].
#' @export
edgeicg <- function(counts, coords, types, lib = NULL, spill_edges = NULL, ex = NULL, block_genes = 500L,
                    window_cells = NULL, panel_genes = 3000L, covariates = NULL,
                    covariate_scope = c("all", "niche"), workers = 0L, backend = c("auto", "rust", "r"), threads = NULL, calibrate = TRUE, calib_cond = TRUE, attribute = TRUE,
                    B = 30L, min_shift = 450, q_thr = 0.05, min_cells = NULL, min_det = 0.01, min_exposed = 30L,
                    exposed_at = 0.1, drop_degenerate = TRUE, prefilter = TRUE, seed = 1L, verbose = TRUE,
                    region = NULL, local = FALSE, patch_size = 1000, min_cells_region = 50L, min_exposed_region = 20L,
                    genes = NULL, leak_genes = c("all", "tested"), leak = NULL, ...) {
  backend <- match.arg(backend); covariate_scope <- match.arg(covariate_scope); leak_genes <- match.arg(leak_genes)
  has_rs <- requireNamespace("edgeicgrs", quietly = TRUE)
  if (backend == "rust" && !has_rs) stop("backend = 'rust' needs the edgeicgrs package")
  backend <- if (backend == "r" || !has_rs) "r" else "rust"
  if (is.null(threads)) threads <- if (workers > 0) 1L else 0L
  dots <- list(...)
  opt <- function(nm, def) if (is.null(dots[[nm]])) def else dots[[nm]]
  pv <- counts_provider(counts)
  lev <- if (is.factor(types)) levels(droplevels(types)) else sort(unique(as.character(types)))
  tcode <- match(as.character(types), lev)
  if (anyNA(tcode)) stop("types has missing values")
  coords <- as.matrix(coords); n <- length(tcode)
  if (nrow(coords) != n) stop("coords and types differ in length")
  if (!is.null(pv$n_cells) && pv$n_cells != n) stop("counts and types differ in the number of cells")
  if (is.null(lib)) {
    if (!is.matrix(counts) && !inherits(counts, "Matrix")) stop("lib is required unless counts is a matrix")
    lib <- Matrix::colSums(counts)
  }
  K <- length(lev)
  if (!is.null(covariates)) {
    covariates <- as.matrix(covariates)
    if (nrow(covariates) != n || !is.numeric(covariates) || anyNA(covariates)) stop("covariates must be a numeric n x q matrix without NA")
  }
  if (is.null(min_cells)) {
    min_cells <- min(60L, 2L * (1L + 3L * (K - 1L) + if (is.null(covariates)) 0L else ncol(covariates)))
  }
  min_cells <- as.integer(min_cells)
  type_n <- tabulate(tcode, K)
  if (max(type_n) < min_cells) {
    message("edgeicg: no cell type has at least min_cells = ", min_cells, " cells (largest: ", max(type_n),
            "); nothing to test")
    genes <- if (is.null(genes)) pv$genes else unique(genes)
    r <- .edgeicg_empty(genes, K)
    r$types <- lev; r$n_cells <- n; r$type_n <- type_n
    r$settings <- list(min_cells = min_cells, calibrate = calibrate)
    return(structure(r, class = "EdgeICG"))
  }
  strata <- list()
  if (!is.null(region)) {
    if (length(region) != n) stop("region must have one value per cell")
    rl <- if (is.factor(region)) levels(droplevels(region)) else sort(unique(as.character(region[!is.na(region)])))
    if (length(rl) < 2) stop("region needs at least two regions")
    strata$region <- list(code = match(as.character(region), rl), labels = rl)
  }
  if (isTRUE(local)) strata$local <- edgeicg_patches(coords, patch_size)
  if (is.null(window_cells)) window_cells <- if (identical(pv$kind, "memory")) Inf else            # a store is read
    min(getOption("edgeicg.window_cells", 50000), ceiling(n / max(1L, workers)))                    # once per sweep
  if (is.null(ex)) ex <- edgeicg_exposures(coords, tcode, K, k = opt("k", 8), sigma_p = opt("sigma_p", 25),
                                           sigma_n = opt("sigma_n", 150), bin = opt("bin", 10))
  if (calibrate && is.null(ex$fields)) stop("ex has no niche fields; build it with edgeicg_exposures()")
  xp_b <- .edgeicg_xp(backend, as.integer(threads))
  huber <- opt("huber", 1.345); sc <- opt("scale", 100)
  gi <- NULL
  if (!is.null(genes)) {
    gi <- match(unique(genes), pv$genes)
    if (anyNA(gi)) stop(sum(is.na(gi)), " of genes not found in counts, e.g. ", unique(genes)[is.na(gi)][1])
    if (length(gi) == length(pv$genes)) gi <- NULL
  }
  leak_src <- if (!is.null(leak)) "supplied" else if (is.null(gi)) "all" else leak_genes
  if (!is.null(leak)) {
    lp <- if (inherits(leak, "EdgeICG")) {
      data.frame(receiver = leak$types[leak$pairs$receiver], sender = leak$types[leak$pairs$sender],
                 lambda = leak$pairs$lambda, mu = leak$pairs$mu)
    } else as.data.frame(leak)
    if (!all(c("receiver", "sender", "lambda", "mu") %in% names(lp))) stop("leak needs receiver, sender, lambda and mu")
    leak <- data.frame(receiver = match(as.character(lp$receiver), lev), sender = match(as.character(lp$sender), lev),
                       lambda = lp$lambda, mu = lp$mu)
  } else if (!is.null(gi) && leak_genes == "all") {                 # stage 1 on all genes, leak only
    if (verbose) message("leak: fitting lambda / mu on all ", length(pv$genes), " genes of counts")
    lk <- .edgeicg_engine(pv, coords, tcode, lib, ex, spill_edges, block_genes, window_cells, min_cells, min_det,
                          min_exposed, exposed_at, huber, sc, FALSE, B, min_shift, seed, FALSE, FALSE, q_thr,
                          opt("max_size", 2), opt("prior_sd", c(0.05, 0.1, 0.2)), opt("cs_level", 0.95),
                          as.integer(workers), xp_b, drop_degenerate, prefilter, panel_genes, verbose, covariates,
                          covariate_scope, list(), as.integer(min_cells_region), as.integer(min_exposed_region),
                          leak_only = TRUE)
    leak <- lk$pairs[, c("receiver", "sender", "lambda", "mu")]
    leak_src <- sprintf("all (%d genes)", lk$n_genes)
  }
  if (!is.null(gi)) pv <- pv$narrow(gi)
  r <- .edgeicg_engine(pv, coords, tcode, lib, ex, spill_edges, block_genes, window_cells, min_cells, min_det,
                       min_exposed, exposed_at, huber, sc, calibrate, B, min_shift, seed,
                       calib_cond, attribute, q_thr, opt("max_size", 2), opt("prior_sd", c(0.05, 0.1, 0.2)),
                       opt("cs_level", 0.95), as.integer(workers), xp_b,
                       drop_degenerate, prefilter, panel_genes, verbose, covariates, covariate_scope,
                       strata, as.integer(min_cells_region), as.integer(min_exposed_region), opt("local_p", 0.1),
                       leak = leak)
  if (leak_src %in% c("all", "tested")) leak_src <- sprintf("%s (%d genes)", leak_src, length(r$genes))
  r$types <- lev; r$n_cells <- n; r$type_n <- tabulate(tcode, K)
  r$settings <- list(min_cells = min_cells, block_genes = block_genes, window_cells = window_cells, workers = workers, backend = backend,
                     threads = threads, covariates = if (is.null(covariates)) 0L else ncol(covariates),
                     covariate_scope = covariate_scope, calibrate = calibrate,
                     calib_cond = calib_cond, B = B, min_shift = min_shift, q_thr = q_thr, exposures = ex$settings,
                     spill = if (is.null(spill_edges)) "knn" else "custom",
                     region = !is.null(region), local = isTRUE(local), patch_size = if (isTRUE(local)) patch_size,
                     min_cells_region = min_cells_region, min_exposed_region = min_exposed_region,
                     leak_genes = leak_src)
  structure(r, class = "EdgeICG")
}

# An EdgeICG result without tests: no receiver type reaches min_cells.
.edgeicg_empty <- function(genes, K) {
  m0 <- matrix(numeric(), 0L, length(genes))
  fit <- list(effect = m0, se = m0, z = m0, q = m0, lfc = m0)
  list(genes = genes, K = K, pairs = data.frame(receiver = integer(), sender = integer(), lambda = numeric(),
       mu = numeric()), fitted = integer(), type_mean = matrix(0, K, K),
       scale_fits = list(contact = fit, paracrine = fit), screen = fit, attribution = data.table::data.table(
         receiver = integer(), sender = integer(), gene = character(), z_marg = numeric(), q_marg = numeric(),
         z_cond = numeric(), pip = numeric(), alpha = numeric(), p_two = numeric(), in_cs = logical(),
         cs_size = integer()),
       calibration = NULL, timing = numeric(), n_dropped = 0L, n_prefiltered = 0L, n_windows = 0L)
}

#' Results table
#'
#' @param x an `EdgeICG` object.
#' @param scale one or more of `"contact"`, `"paracrine"`, `"niche"` (niche = calibrated screen).
#' @param q_thr keep entries with q below this (1 = all testable entries).
#' @param min_lfc keep entries whose effect size, `lfc` (log2 fold change between receivers at the 10th and 90th
#'   percentile of their exposure to the sender), is at least this in absolute value. Default 0.25: at
#'   whole-section scale significance alone also picks up small confounded associations (zonation, sub-states);
#'   on the Atera kidney section the floor cut null-gene calls from 27% to 6% without losing planted effects of
#'   at least 0.5. Use 0 for significance only.
#' @return data.table: receiver, sender, gene, scale, effect, lfc, se, z, q, plus for niche the top attributed
#'   sender (`top_sender`, `alpha`, `cs_size`).
#' @export
edgeicg_results <- function(x, scale = c("contact", "paracrine", "niche"), q_thr = 0.05, min_lfc = 0.25) {
  stopifnot(inherits(x, "EdgeICG"))
  scale <- match.arg(scale, several.ok = TRUE)
  lev <- x$types; pr <- x$pairs
  out <- data.table::rbindlist(lapply(scale, function(s) {
    f <- if (s == "niche") x$screen else x$scale_fits[[s]]
    keep <- !is.na(f$q) & f$q < q_thr
    if (min_lfc > 0) keep <- keep & !is.na(f$lfc) & abs(f$lfc) >= min_lfc
    w <- which(keep, arr.ind = TRUE)
    if (!nrow(w)) return(NULL)
    data.table::data.table(receiver = lev[pr$receiver[w[, 1]]], sender = lev[pr$sender[w[, 1]]], gene = x$genes[w[, 2]],
                           scale = s, effect = f$effect[w], lfc = f$lfc[w], se = f$se[w], z = f$z[w], q = f$q[w])
  }))
  if (nrow(out) && "niche" %in% scale && !is.null(x$attribution) && nrow(x$attribution)) {
    top <- x$attribution[, .SD[1], by = list(receiver, gene)][, list(receiver = lev[receiver], gene, top_sender = lev[sender],
                                                                     alpha, cs_size)]
    out <- merge(out, top[, scale := "niche"], by = c("receiver", "gene", "scale"), all.x = TRUE, sort = FALSE)
  }
  if (nrow(out)) data.table::setorder(out, scale, q)
  out[]
}

#' @export
print.EdgeICG <- function(x, ...) {
  cat(sprintf("EdgeICG: %d cells, %d types, %d genes, %d pairs (%s)\n", x$n_cells, length(x$types), length(x$genes),
              nrow(x$pairs), paste(names(x$timing), sprintf("%.0f s", x$timing), collapse = ", ")))
  if (!is.null(x$settings$leak_genes)) cat(sprintf("  leak fitted on: %s\n", x$settings$leak_genes))
  nh <- function(f) c(sum(f$q < 0.05 & abs(f$lfc) >= 0.25, na.rm = TRUE), sum(f$q < 0.05, na.rm = TRUE))
  cat("  hits: q < 0.05 and |log2 FC| >= 0.25 (q < 0.05 alone)\n")
  for (s in c("contact", "paracrine")) { h <- nh(x$scale_fits[[s]]); cat(sprintf("  %-9s %d (%d)\n", s, h[1], h[2])) }
  h <- nh(x$screen)
  cat(sprintf("  niche     %d (%d)%s\n", h[1], h[2], if (is.null(x$calibration)) ", uncalibrated" else ", screen"))
  if (!is.null(x$attribution) && nrow(x$attribution))
    cat(sprintf("            attributed (receiver, gene): %d\n", nrow(unique(x$attribution[, list(receiver, gene)]))))
  hq <- function(R) vapply(R$het, function(h) sum(h$q < 0.05, na.rm = TRUE), 1)
  if (!is.null(x$region)) cat(sprintf("  region-dependent (q_het < 0.05) over %d regions: %s\n", length(x$region$labels),
                                      paste(names(x$region$het), hq(x$region), collapse = ", ")))
  if (!is.null(x$local)) cat(sprintf("  local (q_het < 0.05) over %d patches: %s\n", length(x$local$labels),
                                     paste(names(x$local$het), hq(x$local), collapse = ", ")))
  invisible(x)
}

# Note on the response scale: effects are on the composition scale (counts / library). When a receiver strongly
# up-regulates highly expressed genes, the shares of its other genes fall slightly; these compositional
# shifts are real on this scale and show up as small negative effects, mostly in small panels.

# Dense product kernels of the sweeps: xp = t(A) %*% B, mp = A %*% B (A may be a vector). The Rust versions
# (edgeicgrs) are an optimised matrix multiply that does not depend on the BLAS R is linked to.
.edgeicg_xp <- function(backend, threads) {
  if (backend == "r") return(list(xp = function(A, B) crossprod(A, B), mp = function(A, B) A %*% B))
  force(threads)
  list(xp = function(A, B) edgeicgrs::crossprod_rs(A, B, threads), mp = function(A, B) edgeicgrs::matprod_rs(A, B, threads))
}

#' Spatial basis for large-scale trends
#'
#' A tensor-product natural-spline basis of the cell positions, to pass as `covariates` to [edgeicg()]: it absorbs
#' smooth spatial trends (tissue zonation, gradients of cell states within a type) so that they are not mistaken
#' for exposure effects. `df` per axis sets the scale of the trends removed: with a section of width W, trends
#' longer than about W / df are absorbed.
#'
#' @param coords n x 2 cell centroids.
#' @param df degrees of freedom per axis (default 4: 4 + 4 + 16 = 24 columns).
#' @return a numeric matrix (cells x (2 df + df^2)), columns centred.
#' @export
spatial_basis <- function(coords, df = 4L) {
  bx <- splines::ns(coords[, 1], df = df); by <- splines::ns(coords[, 2], df = df)
  X <- cbind(bx, by, do.call(cbind, lapply(seq_len(df), function(j) bx * by[, j])))
  X <- sweep(X, 2, colMeans(X))
  colnames(X) <- c(paste0("sx", seq_len(df)), paste0("sy", seq_len(df)), paste0("sxy", seq_len(df^2)))
  X
}

