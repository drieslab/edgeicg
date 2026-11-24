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

test_that("windowed runs (cell streaming) equal the one-window run", {
  sim <- simulate_edgeicg(n = 3000, G = 100, seed = 7, gamma = 1.2)
  run <- function(wc, w) edgeicg(sim$counts, sim$coords, sim$types, min_shift = 150, B = 8, block_genes = 40,
                                 window_cells = wc, workers = w, backend = "r", verbose = FALSE)
  one <- run(Inf, 0)
  expect_equal(one$n_windows, 1L)
  runs <- list(run(700, 0))
  if (requireNamespace("mirai", quietly = TRUE)) runs <- c(runs, list(run(700, 2)))
  for (f in runs) {
    expect_gt(f$n_windows, 3L)
    expect_equal(f$pairs$lambda, one$pairs$lambda, tolerance = 1e-8)
    for (s in c("contact", "paracrine", "niche")) expect_equal(f$scale_fits[[s]]$z, one$scale_fits[[s]]$z, tolerance = 1e-8)
    expect_equal(f$screen$z, one$screen$z, tolerance = 1e-8)
    expect_equal(f$cond$z, one$cond$z, tolerance = 1e-8)
    expect_equal(f$calibration$raw_marg, one$calibration$raw_marg, tolerance = 1e-8)
    expect_equal(f$attribution$alpha, one$attribution$alpha, tolerance = 1e-8)
  }
})

test_that("the prefilter drops only genes without a testable entry", {
  sim <- simulate_edgeicg(n = 2500, G = 100, seed = 7, gamma = 1.2)
  dead <- Matrix::Matrix(0, 20, ncol(sim$counts), sparse = TRUE,
                         dimnames = list(paste0("dead", 1:20), colnames(sim$counts)))
  dead[, 1:3] <- 1                                                   # detected in 3 cells only
  counts <- rbind(sim$counts, dead); lib <- Matrix::colSums(counts)
  run <- function(pf, wc = Inf) edgeicg(counts, sim$coords, sim$types, lib = lib, min_shift = 150, B = 8,
                                       block_genes = 40, window_cells = wc, prefilter = pf, backend = "r", verbose = FALSE)
  on <- run(TRUE); off <- run(FALSE)
  expect_equal(on$n_prefiltered, 20L)
  expect_identical(on$genes, rownames(sim$counts))
  expect_true(all(is.na(off$scale_fits$contact$z[, 101:120])))
  # kept genes: the same tests, up to the leak (lambda / mu) fit across genes, which the dropped all-zero genes
  # pulled towards zero
  expect_gt(cor(as.vector(on$scale_fits$contact$z), as.vector(off$scale_fits$contact$z[, 1:100]),
                use = "complete.obs"), 0.995)
  win <- run(TRUE, 700)
  expect_equal(win$scale_fits$contact$z, on$scale_fits$contact$z, tolerance = 1e-8)
})

test_that("gene panels (windowed, more genes than panel_genes) equal the one-window run", {
  sim <- simulate_edgeicg(n = 3000, G = 100, seed = 7, gamma = 1.2)
  run <- function(wc, w, pg = 3000L) edgeicg(sim$counts, sim$coords, sim$types, min_shift = 150, B = 8, block_genes = 20,
                                            window_cells = wc, panel_genes = pg, workers = w, backend = "r", verbose = FALSE)
  one <- run(Inf, 0)
  runs <- list(run(700, 0, 40))
  if (requireNamespace("mirai", quietly = TRUE)) runs <- c(runs, list(run(700, 2, 40)))
  for (f in runs) {
    expect_equal(f$pairs$lambda, one$pairs$lambda, tolerance = 1e-8)
    expect_equal(f$scale_fits$contact$z, one$scale_fits$contact$z, tolerance = 1e-8)
    expect_equal(f$screen$z, one$screen$z, tolerance = 1e-8)
    expect_equal(f$attribution$alpha, one$attribution$alpha, tolerance = 1e-8)
  }
})

test_that("effect sizes, the min_lfc floor and niche-only covariates", {
  sim <- simulate_edgeicg(n = 2500, G = 100, seed = 7, gamma = 1.2)
  f <- edgeicg(sim$counts, sim$coords, sim$types, min_shift = 150, B = 8, backend = "r", verbose = FALSE)
  expect_true(all(c("lfc") %in% names(f$scale_fits$contact)) && !is.null(f$screen$lfc))
  r0 <- edgeicg_results(f, "contact", min_lfc = 0); r1 <- edgeicg_results(f, "contact")
  expect_true(all(abs(r1$lfc) >= 0.25)); expect_lte(nrow(r1), nrow(r0))
  # lfc: fitted means at the exposure deciles, symmetric around the receiver mean, bounded
  p <- which(f$pairs$receiver == 1 & f$pairs$sender == 2)[1]; g <- 5
  h <- f$scale_fits$contact$effect[p, g] * f$pairs$dx_contact[p] / 2; tm <- f$type_mean[1, g]
  expect_equal(f$scale_fits$contact$lfc[p, g], log2((tm + h) / (tm - h)))
  expect_true(all(abs(f$scale_fits$contact$lfc) <= log2(1.9 / 0.1) + 1e-12, na.rm = TRUE))
  expect_equal(.edgeicg_lfc(c(-1e3, 0.01, 1e3), 1, 1), c(-log2(19), log2(1.005 / 0.995), log2(19)))
  fn <- edgeicg(sim$counts, sim$coords, sim$types, min_shift = 150, B = 8, covariates = spatial_basis(sim$coords, 3),
                covariate_scope = "niche", backend = "r", verbose = FALSE)
  expect_equal(fn$scale_fits$contact$z, f$scale_fits$contact$z)        # contact and paracrine untouched
  expect_false(isTRUE(all.equal(fn$screen$z, f$screen$z)))              # the niche screen is adjusted
})

test_that("a rare last cell type (below min_cells) is skipped, not an error", {
  ty <- sim$types; lv <- if (is.factor(ty)) levels(ty) else sort(unique(ty))
  last <- which(ty == lv[length(lv)])
  ty[last[-(1:5)]] <- lv[1]                                         # keep 5 cells of the last type
  f <- edgeicg(sim$counts, sim$coords, ty, lib = sim$lib, B = 8, min_shift = 150, min_cells = 20, verbose = FALSE)
  expect_s3_class(f, "EdgeICG")
  expect_false(lv[length(lv)] %in% as.character(f$pairs$receiver))
  w <- edgeicg(sim$counts, sim$coords, ty, lib = sim$lib, B = 8, min_shift = 150, min_cells = 20,
               window_cells = 700, backend = "r", verbose = FALSE)
  expect_equal(w$screen$z, f$screen$z, tolerance = 1e-8)
})
