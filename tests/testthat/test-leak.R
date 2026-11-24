# Leak coefficients fitted on a broader gene set than the tested genes (genes / leak_genes / leak)

siml <- simulate_edgeicg(n = 4000, G = 200, seed = 21)
fitl <- function(...) edgeicg(siml$counts, siml$coords, siml$types, lib = siml$lib, B = 6, min_shift = 150,
                              verbose = FALSE, ...)
sub <- rownames(siml$counts)[seq(1, 200, by = 5)]
full <- fitl()
cols <- function(f, g) match(g, f$genes)

test_that("leak_genes = 'all' fits lambda on all genes, so tested-gene effects equal the full fit", {
  f <- fitl(genes = sub)
  expect_equal(f$genes, intersect(full$genes, sub))
  expect_equal(f$pairs$lambda, full$pairs$lambda)
  expect_equal(f$pairs$mu, full$pairs$mu)
  for (s in c("contact", "paracrine")) for (m in c("effect", "se", "z"))
    expect_equal(f$scale_fits[[s]][[m]], full$scale_fits[[s]][[m]][, cols(full, f$genes), drop = FALSE])
  expect_match(f$settings$leak_genes, "^all \\(")
})

test_that("leak_genes = 'tested' equals a fit on the tested genes only", {
  f <- fitl(genes = sub, leak_genes = "tested")
  o <- edgeicg(siml$counts[sub, ], siml$coords, siml$types, lib = siml$lib, B = 6, min_shift = 150, verbose = FALSE)
  expect_equal(f$pairs, o$pairs)
  expect_equal(f$scale_fits, o$scale_fits)
  expect_false(isTRUE(all.equal(f$pairs$lambda, full$pairs$lambda)))
})

test_that("a supplied leak (fit or table) is used as is", {
  f1 <- fitl(genes = sub, leak = full)
  f2 <- fitl(genes = sub, leak = data.frame(receiver = full$types[full$pairs$receiver],
                                            sender = full$types[full$pairs$sender],
                                            lambda = full$pairs$lambda, mu = full$pairs$mu))
  f3 <- fitl(genes = sub)
  expect_equal(f1$scale_fits, f3$scale_fits); expect_equal(f2$scale_fits, f3$scale_fits)
  expect_equal(f1$settings$leak_genes, "supplied")
  expect_error(fitl(genes = c(sub, "nope")), "not found")
})

test_that("the leak pass gives the same lambda with windows and workers", {
  f <- fitl(genes = sub, window_cells = 1200)
  expect_equal(f$pairs$lambda, full$pairs$lambda, tolerance = 1e-8)
})
