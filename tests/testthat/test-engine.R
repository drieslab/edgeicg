sim <- simulate_edgeicg(n = 3000, G = 120, seed = 3, gamma = 1.2)

test_that("results do not depend on the block size", {
  a <- edgeicg(sim$counts, sim$coords, sim$types, lib = sim$lib, block_genes = 120, B = 8, min_shift = 150, verbose = FALSE)
  b <- edgeicg(sim$counts, sim$coords, sim$types, lib = sim$lib, block_genes = 25, B = 8, min_shift = 150, verbose = FALSE)
  for (s in c("contact", "paracrine", "niche")) expect_equal(a$scale_fits[[s]]$z, b$scale_fits[[s]]$z, tolerance = 1e-8)
  expect_equal(a$pairs$lambda, b$pairs$lambda, tolerance = 1e-10)
  expect_equal(a$screen$q, b$screen$q, tolerance = 1e-8)
  expect_equal(a$calibration$raw_marg, b$calibration$raw_marg, tolerance = 1e-10)
})

test_that("true contact effects are found with the right sign", {
  f <- edgeicg(sim$counts, sim$coords, sim$types, lib = sim$lib, B = 8, min_shift = 150, verbose = FALSE)
  r <- edgeicg_results(f, "contact", q_thr = 1)
  tr <- merge(sim$truth[scale == "contact"], r, by = c("receiver", "sender", "gene"))
  expect_gt(mean(tr$q < 0.05 & sign(tr$z) == tr$dir), 0.6)
  expect_s3_class(f, "EdgeICG")
  expect_output(print(f), "EdgeICG")
})

test_that("mirai daemons give identical results", {
  skip_if_not_installed("mirai")
  skip_on_cran()
  a <- edgeicg(sim$counts, sim$coords, sim$types, lib = sim$lib, block_genes = 40, B = 8, min_shift = 150, verbose = FALSE)
  b <- edgeicg(sim$counts, sim$coords, sim$types, lib = sim$lib, block_genes = 40, B = 8, min_shift = 150, workers = 2, verbose = FALSE)
  expect_equal(a$scale_fits$contact$z, b$scale_fits$contact$z, tolerance = 1e-10)
  expect_equal(a$screen$z, b$screen$z, tolerance = 1e-10)
})

test_that("without effects, contact and paracrine calls are rare (leakage removed)", {
  s0 <- simulate_edgeicg(n = 3000, G = 120, seed = 4, gamma = 0)
  f <- edgeicg(s0$counts, s0$coords, s0$types, lib = s0$lib, B = 8, min_shift = 150, attribute = FALSE, verbose = FALSE)
  for (sc in c("contact", "paracrine")) {
    z <- f$scale_fits[[sc]]$z
    expect_lt(mean(abs(z) > 3, na.rm = TRUE), 0.01)
  }
})
