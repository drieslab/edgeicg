test_that("in-memory provider: rows, order and slices agree", {
  m <- Matrix::rsparsematrix(30, 40, 0.3); m@x <- abs(round(m@x * 5)) + 1
  dimnames(m) <- list(paste0("g", 1:30), paste0("c", 1:40))
  pv <- counts_provider(m); ord <- sample(40)
  expect_equal(as.matrix(pv$rows(3, 9, 2:5)), as.matrix(Matrix::t(m[2:5, 3:9])))
  po <- pv$order(ord)
  expect_equal(as.matrix(po$rows(1, 10, 1:30)), as.matrix(Matrix::t(m[, ord[1:10]])))
  expect_equal(as.matrix(po$slice(5, 8, 7:9)), as.matrix(m[7:9, ord[5:8]]))
})

test_that("GiottoDisk store: counts_provider() wraps it, the bake keeps the order, and the fit equals in-memory", {
  skip_if_not_installed("GiottoDisk")
  sim <- simulate_edgeicg(n = 2500, G = 100, seed = 7, gamma = 1.2)
  path <- tempfile(fileext = ".parquet"); on.exit(unlink(path, recursive = TRUE))
  st <- GiottoDisk::storeWrite(GiottoDisk::parquetExprStore(path = path), sim$counts)
  pv <- counts_provider(st)
  expect_identical(pv$kind, "store")
  ord <- order(sim$coords[, 2]); po <- pv$order(ord); on.exit(po$cleanup(), add = TRUE)
  expect_equal(as.matrix(po$rows(11, 60, 1:100)), as.matrix(Matrix::t(sim$counts[, ord[11:60]])), ignore_attr = TRUE)
  lib <- Matrix::colSums(sim$counts)
  run <- function(x, wc) edgeicg(x, sim$coords, sim$types, lib = lib, min_shift = 150, B = 8, block_genes = 40,
                                 window_cells = wc, backend = "r", verbose = FALSE)
  one <- run(sim$counts, Inf); fs <- run(st, 700)
  expect_gt(fs$n_windows, 3L)
  for (s in c("contact", "paracrine", "niche")) expect_equal(fs$scale_fits[[s]]$z, one$scale_fits[[s]]$z, tolerance = 1e-8)
  expect_equal(fs$screen$z, one$screen$z, tolerance = 1e-8)
  expect_equal(fs$attribution$alpha, one$attribution$alpha, tolerance = 1e-8)
})
