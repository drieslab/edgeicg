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
                                            lambda = full$pairs$lambda, mu = full$pairs$mu,
                                            leak_status = full$pairs$leak_status))
  f3 <- fitl(genes = sub)
  expect_equal(f1$scale_fits, f3$scale_fits); expect_equal(f2$scale_fits, f3$scale_fits)
  expect_equal(f1$settings$leak_genes, "supplied")
  expect_error(fitl(genes = c(sub, "nope")), "not found")
})

test_that("the leak pass gives the same lambda with windows and workers", {
  f <- fitl(genes = sub, window_cells = 1200)
  expect_equal(f$pairs$lambda, full$pairs$lambda, tolerance = 1e-8)
})

test_that("leak pooling: pairs are shrunk toward the receiver median by contact support; raw values are kept", {
  P <- full$pairs
  expect_true(all(c("support", "leak_status", "lambda_raw", "mu_raw") %in% names(P)))
  ref <- edgeicg:::.leak_pool(P$receiver, P$lambda_raw, P$mu_raw, P$support, "median", 30)
  expect_equal(P$lambda, ref$lambda); expect_equal(P$mu, ref$mu); expect_equal(P$leak_status, ref$status)
  ok <- P$leak_status == "ok"; expect_true(all(P$support[ok] >= 30))
  expect_true(all(P$lambda >= 0 & P$lambda <= 1))
  n <- fitl(genes = sub, leak_pool = "none")
  expect_equal(n$pairs$lambda[n$pairs$leak_status == "ok"], n$pairs$lambda_raw[n$pairs$leak_status == "ok"])
  expect_true(all(is.na(n$pairs$lambda[n$pairs$leak_status != "ok"])))
  # a plain table (no status) is pooled; with pool = "none" and raw values it reproduces the unpooled fit
  raw <- data.frame(receiver = full$types[P$receiver], sender = full$types[P$sender], lambda = P$lambda_raw, mu = P$mu_raw)
  expect_equal(fitl(genes = sub, leak = raw)$pairs$lambda, P$lambda)
})
