#' Counts providers
#'
#' EdgeICG reads raw counts only in blocks. A provider is a list with
#' * `genes` (character) and `n_cells`;
#' * `rows(lo, hi, gi)`: raw counts of cells `lo:hi` (in the provider's current cell order) for genes `gi`, as a
#'   cells x genes matrix (sparse `dgCMatrix` preferred; dense is accepted);
#' * `order(ord)`: a provider whose cell `i` is this provider's cell `ord[i]`. The engine visits cells in
#'   spatially ordered windows (see `window_cells` in [edgeicg()]) and asks the provider for that order once;
#' * `cleanup()`: removes anything `order()` created;
#' * `narrow(gi)`: a provider for genes `gi` only (used by the gene prefilter);
#' * `split(list_of_gi)`: one provider per gene panel (windowed runs with more than `panel_genes` genes);
#' * optionally `slice(lo, hi, gi)`: the same block as a compact genes x cells matrix, for providers whose data
#'   live in memory (the block is sent to the worker with its task);
#' * `read(gi)`: all cells for genes `gi` (cells x genes), kept for direct use.
#'
#' `counts_provider()` wraps an in-memory genes x cells matrix (sparse or dense). `counts_provider_store()` wraps a
#' GiottoDisk expression store, so whole sections never sit in memory: workers read their own windows from it.
#'
#' @param counts genes x cells matrix (`dgCMatrix` or dense), a GiottoDisk expression store (wrapped with
#'   `counts_provider_store()` for all its genes and cells), or an existing provider (returned as is).
#' @return a provider list.
#' @export
counts_provider <- function(counts) {
  if (is.list(counts) && is.function(counts$rows)) return(counts)
  if (methods::is(counts, "parquetExprBase"))                       # a GiottoDisk expression store (or view)
    return(counts_provider_store(counts, rownames(counts), colnames(counts)))
  if (is.null(rownames(counts))) stop("counts needs gene names as row names")
  .provider_memory(counts, NULL)
}

.provider_memory <- function(counts, perm) {
  force(counts); force(perm)
  cols <- function(lo, hi) if (is.null(perm)) lo:hi else perm[lo:hi]
  blk <- function(lo, hi, gi) methods::as(counts[gi, cols(lo, hi), drop = FALSE], "CsparseMatrix")
  list(genes = rownames(counts), n_cells = ncol(counts), kind = "memory",
       rows = function(lo, hi, gi) Matrix::t(blk(lo, hi, gi)),
       slice = blk,
       order = function(ord) .provider_memory(counts, if (is.null(perm)) ord else perm[ord]),
       narrow = function(gi) .provider_memory(counts[gi, , drop = FALSE], perm),
       split = function(gl) lapply(gl, function(gi) .provider_memory(counts[gi, , drop = FALSE], perm)),
       cleanup = function() invisible(),
       read = function(gi) Matrix::t(methods::as(counts[gi, , drop = FALSE], "CsparseMatrix")))
}

#' @rdname counts_provider
#' @param store a GiottoDisk expression store (genes x cells), e.g.
#'   `GiottoClass::getExpression(g, values = "raw", output = "exprObj")[]`.
#' @param genes genes to analyse (store feature IDs).
#' @param cells cell IDs in analysis order (store cell IDs).
#' @details For windowed runs, `order()` bakes the analysis genes into a temporary store in the engine's cell order
#'   (`GiottoDisk::storeWrite()` of the view, the step GiottoDisk's streaming PCA uses to narrow a store). The
#'   expression store is sorted by cell, so a window and its halo are then a contiguous row range that prunes
#'   parquet row groups, and only the analysis genes are scanned. The temporary store is removed by `cleanup()`.
#' @export
counts_provider_store <- function(store, genes, cells) {
  if (!requireNamespace("GiottoDisk", quietly = TRUE)) stop("counts_provider_store() needs the GiottoDisk package")
  force(store); force(genes); force(cells)
  get <- function(v) {
    m <- GiottoDisk::storeRead(v, output = "dgcmatrix", max_rows = Inf, max_cols = Inf)
    Matrix::t(m)
  }
  list(genes = genes, n_cells = length(cells), kind = "store",
       rows = function(lo, hi, gi) get(store[genes[gi], cells[lo:hi]])[cells[lo:hi], genes[gi], drop = FALSE],
       order = function(ord) {
         path <- tempfile("edgeicg_bake_", fileext = ".parquet")
         baked <- GiottoDisk::storeWrite(GiottoDisk::parquetExprStore(path = path), store[genes, cells[ord]])
         .provider_baked(baked, genes, path)
       },
       narrow = function(gi) counts_provider_store(store, genes[gi], cells),
       split = function(gl) lapply(gl, function(gi) counts_provider_store(store, genes[gi], cells)),
       cleanup = function() invisible(),
       read = function(gi) get(store[genes[gi], cells])[cells, genes[gi], drop = FALSE])
}

# A store baked by counts_provider_store()$order(): rows are the engine's internal cell order and its genes are
# exactly `genes`, so windows are read by position (contiguous cell ranges prune row groups; all genes are read
# when a task needs all of them, which is the case for the windowed sweeps).
.provider_baked <- function(baked, genes, path) {
  force(baked); force(genes); force(path)
  G <- length(genes)
  list(genes = genes, n_cells = ncol(baked), kind = "store (baked)",
       rows = function(lo, hi, gi) {
         v <- if (length(gi) == G && identical(as.integer(gi), seq_len(G))) baked[, lo:hi] else baked[gi, lo:hi]
         m <- Matrix::t(GiottoDisk::storeRead(v, output = "dgcmatrix", max_rows = Inf, max_cols = Inf))
         m[, genes[gi], drop = FALSE]
       },
       order = function(ord) stop("a baked store is already ordered"),
       narrow = function(gi) {                                       # re-bake with the kept genes only
         p2 <- tempfile("edgeicg_bake_", fileext = ".parquet")
         b2 <- GiottoDisk::storeWrite(GiottoDisk::parquetExprStore(path = p2), baked[gi, ])
         unlink(path, recursive = TRUE)
         .provider_baked(b2, genes[gi], p2)
       },
       split = function(gl) {                                        # one baked store per gene panel
         out <- lapply(gl, function(gi) {
           p2 <- tempfile("edgeicg_bake_", fileext = ".parquet")
           .provider_baked(GiottoDisk::storeWrite(GiottoDisk::parquetExprStore(path = p2), baked[gi, ]), genes[gi], p2)
         })
         unlink(path, recursive = TRUE)
         out
       },
       cleanup = function() unlink(path, recursive = TRUE),
       read = function(gi) stop("read() is not available on a baked store; use rows()"))
}
