#' Counts providers
#'
#' EdgeICG reads raw counts only in gene blocks. A provider is a list with `genes` (character) and
#' `read(gene_index)`, which returns a cells x genes numeric matrix of raw counts for those genes, with cells
#' in the analysis order. `counts_provider()` wraps an in-memory genes x cells matrix (sparse or dense);
#' `counts_provider_store()` wraps a GiottoDisk expression store, so whole sections never sit in memory.
#'
#' @param counts genes x cells matrix (`dgCMatrix` or dense), or an existing provider (returned as is).
#' @return a provider list.
#' @export
counts_provider <- function(counts) {
  if (is.list(counts) && is.function(counts$read)) return(counts)
  if (is.null(rownames(counts))) stop("counts needs gene names as row names")
  force(counts)
  list(genes = rownames(counts), n_cells = ncol(counts), kind = "memory",
       read = function(gi) as.matrix(Matrix::t(counts[gi, , drop = FALSE])))
}

#' @rdname counts_provider
#' @param store a GiottoDisk expression store (genes x cells), e.g.
#'   `GiottoClass::getExpression(g, values = "raw", output = "exprObj")[]`.
#' @param genes genes to analyse (store feature IDs).
#' @param cells cell IDs in analysis order (store cell IDs).
#' @export
counts_provider_store <- function(store, genes, cells) {
  if (!requireNamespace("GiottoDisk", quietly = TRUE)) stop("counts_provider_store() needs the GiottoDisk package")
  force(store); force(genes); force(cells)
  list(genes = genes, n_cells = length(cells), kind = "store",
       read = function(gi) {
         m <- GiottoDisk::storeRead(store[genes[gi], cells], output = "dgcmatrix", max_rows = Inf, max_cols = Inf)
         as.matrix(Matrix::t(m[genes[gi], cells, drop = FALSE]))
       })
}
