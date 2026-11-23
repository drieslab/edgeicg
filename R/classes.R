#' Classify hits by leakage risk
#'
#' Flags hits whose gene marks another cell type, so leaked or mis-segmented transcripts can be told apart
#' from candidate biology. Classes (first that applies):
#'
#' * `sender_marker`: the gene marks the tested sender (sender specificity > `spec_thr`, or in `markers`).
#' * `other_marker`: the gene marks a third type (not the receiver).
#' * `own_down`: a gene of the receiver goes down next to the sender, where the receiver expresses it far more
#'   than the sender (receiver specificity against the sender > `spec_thr`, or in `markers` of the receiver).
#'   Next to that sender, mixed or mis-segmented cells and dilution by the sender's transcripts lower such
#'   genes, so these hits need checking before they are read as regulation.
#' * `clean`: none of the above.
#'
#' Niche hits where the same sender and sign recur in at least `regional_n` receiver types are flagged
#' `regional` (shared tissue zone or a regional technical artefact).
#'
#' Specificity of gene g for type k against type r is `m_kg / (m_kg + m_rg)` with `m` the mean
#' composition per type from the fit (`x$type_mean`).
#'
#' @param x an `EdgeICG` object.
#' @param hits a table from [edgeicg_results()].
#' @param markers optional named list (type -> character vector of marker genes).
#' @param types_considered types used as "other" types (default: all fitted types).
#' @param spec_thr,regional_n thresholds.
#' @return `hits` with columns `class` and `regional`.
#' @export
edgeicg_classify <- function(x, hits, markers = NULL, types_considered = NULL, spec_thr = 0.8, regional_n = 3) {
  h <- data.table::copy(data.table::as.data.table(hits))
  if (!nrow(h)) return(h)
  lev <- x$types; mt <- x$type_mean; colnames(mt) <- x$genes; rownames(mt) <- lev
  oth_all <- types_considered %||% lev[x$fitted]
  gi <- match(h$gene, x$genes); ri <- match(h$receiver, lev); si <- match(h$sender, lev)
  spec <- function(k, r, g) mt[cbind(k, g)] / (mt[cbind(k, g)] + mt[cbind(r, g)] + 1e-9)
  inm <- function(type, gene) !is.null(markers) && gene %in% markers[[type]]
  h[, sender_mk := spec(si, ri, gi) > spec_thr | mapply(inm, sender, gene)]
  h[, other_mk := mapply(function(r, s, g, gidx) {
    o <- setdiff(oth_all, c(r, s)); if (!length(o)) return(FALSE)
    any(spec(match(o, lev), rep(match(r, lev), length(o)), rep(gidx, length(o))) > spec_thr) || any(vapply(o, inm, TRUE, gene = g))
  }, receiver, sender, gene, gi)]
  h[, own_down := (spec(ri, si, gi) > spec_thr | mapply(inm, receiver, gene)) & z < 0]
  h[, class := data.table::fifelse(sender_mk, "sender_marker", data.table::fifelse(other_mk, "other_marker",
                 data.table::fifelse(own_down, "own_down", "clean")))]
  h[, regional := scale == "niche" & .N >= regional_n, by = list(scale, sender, gene, sign(z))]
  h[, c("sender_mk", "other_mk", "own_down") := NULL]
  h[]
}
