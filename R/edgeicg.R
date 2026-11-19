#' Spillover-aware interaction-changed genes at three scales
#'
#' Tests, for every receiver type, sender type and gene, whether the receiver's expression changes with its
#' exposure to the sender at contact, paracrine and niche scale, after removing transcripts leaked from
#' touching neighbours. See the package description for the model.
#'
#' @param counts raw counts, genes x cells (`dgCMatrix` or dense), or a provider from [counts_provider()] /
#'   [counts_provider_store()]. Only the genes to test; the library size comes from `lib`.
#' @param coords n x 2 cell centroids (microns), in the provider's cell order.
#' @param types cell-type labels (character or factor), one per cell.
#' @param lib library size per cell over **all** genes (default: column sums of `counts` if it is a matrix).
#' @param spill_edges optional touching-cell edge list (`recv`, `send`, `w`), e.g. from
#'   [contact_from_polygons()]; the spillover graph. Default: the kNN contact graph.
#' @param ex optional exposures from [edgeicg_exposures()] (computed if NULL).
#' @param block_genes genes per block (memory: about cells x block_genes x 8 bytes per copy).
#' @param workers number of mirai daemons for the gene blocks (0 = sequential). Results are identical.
#' @param backend `"auto"` or `"r"`. Reserved for a compiled backend; currently always R.
#' @param calibrate run the toroidal-shift calibration of the niche scale (recommended).
#' @param calib_cond also calibrate the conditional (all-sender) niche z (secondary output).
#' @param attribute attribute niche screen hits to senders.
#' @param B,min_shift number of shifts and minimum shift (microns) of the calibration.
#' @param q_thr FDR threshold for the niche screen hits that are attributed.
#' @param min_cells,min_det,min_exposed,exposed_at testability: receiver types with at least `min_cells`
#'   cells; genes detected in at least `min_det` of receivers; at least `min_exposed` receivers with exposure
#'   above `exposed_at`.
#' @param drop_degenerate drop exposure columns without variance within a receiver type (their tests are NA).
#'   FALSE reproduces the 2026 prototype exactly.
#' @param seed seed of the calibration shifts.
#' @param verbose progress messages.
#' @param ... further engine settings: `huber`, `scale`, `max_size`, `prior_sd`, `cs_level`, `k`, `sigma_p`,
#'   `sigma_n`, `bin`.
#' @return an `EdgeICG` object; see [edgeicg_results()].
#' @export
edgeicg <- function(counts, coords, types, lib = NULL, spill_edges = NULL, ex = NULL, block_genes = 500L,
                    workers = 0L, backend = c("auto", "r"), calibrate = TRUE, calib_cond = TRUE, attribute = TRUE,
                    B = 30L, min_shift = 450, q_thr = 0.05, min_cells = 60L, min_det = 0.01, min_exposed = 30L,
                    exposed_at = 0.1, drop_degenerate = TRUE, seed = 1L, verbose = TRUE, ...) {
  backend <- match.arg(backend)
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
    if (is.list(counts)) stop("lib is required with a provider")
    lib <- Matrix::colSums(counts)
  }
  K <- length(lev)
  if (is.null(ex)) ex <- edgeicg_exposures(coords, tcode, K, k = opt("k", 8), sigma_p = opt("sigma_p", 25),
                                           sigma_n = opt("sigma_n", 150), bin = opt("bin", 10))
  if (calibrate && is.null(ex$fields)) stop("ex has no niche fields; build it with edgeicg_exposures()")
  r <- .edgeicg_engine(pv, coords, tcode, lib, ex, spill_edges, block_genes, min_cells, min_det, min_exposed,
                       exposed_at, opt("huber", 1.345), opt("scale", 100), calibrate, B, min_shift, seed,
                       calib_cond, attribute, q_thr, opt("max_size", 2), opt("prior_sd", c(0.05, 0.1, 0.2)),
                       opt("cs_level", 0.95), as.integer(workers), drop_degenerate, verbose)
  r$types <- lev; r$n_cells <- n; r$type_n <- tabulate(tcode, K)
  r$settings <- list(block_genes = block_genes, workers = workers, backend = "r", calibrate = calibrate,
                     calib_cond = calib_cond, B = B, min_shift = min_shift, q_thr = q_thr, exposures = ex$settings,
                     spill = if (is.null(spill_edges)) "knn" else "custom")
  structure(r, class = "EdgeICG")
}

#' Results table
#'
#' @param x an `EdgeICG` object.
#' @param scale one or more of `"contact"`, `"paracrine"`, `"niche"` (niche = calibrated screen).
#' @param q_thr keep entries with q below this (1 = all testable entries).
#' @return data.table: receiver, sender, gene, scale, effect, se, z, q, plus for niche the top attributed
#'   sender (`top_sender`, `alpha`, `cs_size`).
#' @export
edgeicg_results <- function(x, scale = c("contact", "paracrine", "niche"), q_thr = 0.05) {
  stopifnot(inherits(x, "EdgeICG"))
  scale <- match.arg(scale, several.ok = TRUE)
  lev <- x$types; pr <- x$pairs
  out <- data.table::rbindlist(lapply(scale, function(s) {
    f <- if (s == "niche") x$screen else x$scale_fits[[s]]
    w <- which(!is.na(f$q) & f$q < q_thr, arr.ind = TRUE)
    if (!nrow(w)) return(NULL)
    data.table::data.table(receiver = lev[pr$receiver[w[, 1]]], sender = lev[pr$sender[w[, 1]]], gene = x$genes[w[, 2]],
                           scale = s, effect = f$effect[w], se = f$se[w], z = f$z[w], q = f$q[w])
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
  for (s in c("contact", "paracrine")) cat(sprintf("  %-9s hits (q < 0.05): %d\n", s, sum(x$scale_fits[[s]]$q < 0.05, na.rm = TRUE)))
  cat(sprintf("  niche     screen hits (q < 0.05): %d%s\n", sum(x$screen$q < 0.05, na.rm = TRUE),
              if (is.null(x$calibration)) " (uncalibrated)" else ""))
  if (!is.null(x$attribution) && nrow(x$attribution))
    cat(sprintf("            attributed (receiver, gene): %d\n", nrow(unique(x$attribution[, list(receiver, gene)]))))
  invisible(x)
}

# Note on the response scale: effects are on the composition scale (counts / library). When a receiver strongly
# up-regulates highly expressed genes, the shares of its other genes fall slightly; these compositional
# shifts are real on this scale and show up as small negative effects, mostly in small panels.
