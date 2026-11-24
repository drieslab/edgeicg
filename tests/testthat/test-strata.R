# Region and local effects (strata.R)

test_that("stratum slopes and standard errors equal a within-stratum lm with intercept", {
  set.seed(5); n <- 400; D <- 3; g <- sample.int(D, n, TRUE); g[1:10] <- NA
  xt <- rnorm(n); b <- c(0.3, -0.2); r <- matrix(rnorm(n * 2), n, 2)
  ok <- !is.na(g); xbar <- tapply(xt[ok], g[ok], mean); XC <- matrix(0, n, 1); XC[ok, 1] <- xt[ok] - xbar[g[ok]]
  nd <- tabulate(g[ok], D); sxx <- as.vector(tapply(XC[ok, 1]^2, g[ok], sum))
  sm <- edgeicg:::.strata_sums(g, XC, r, D)
  for (k in 1:2) {
    f <- edgeicg:::.strata_fit(b[k], nd, as.vector(xbar), sxx, sm$r1[, k, drop = FALSE], sm$r2[, k, drop = FALSE], sm$c[, k, drop = FALSE])
    for (d in 1:D) {
      i <- which(g == d); yy <- r[i, k] + b[k] * xt[i]
      cf <- summary(stats::lm(yy ~ xt[i]))$coefficients
      expect_equal(f$b[d, 1], cf[2, 1], tolerance = 1e-10)
      expect_equal(f$se[d, 1], cf[2, 2], tolerance = 1e-10)
    }
  }
})

test_that("heterogeneity: Q, tau2 and shrinkage match their definitions", {
  b <- matrix(c(0.1, 0.5, -0.2, 0.3, 0.3, 0.3), 3); se <- matrix(c(0.1, 0.2, 0.15, 0.1, 0.1, 0.1), 3)
  h <- edgeicg:::.strata_het(b, se)
  w <- 1 / se[, 1]^2; m <- sum(w * b[, 1]) / sum(w); Q <- sum(w * (b[, 1] - m)^2)
  expect_equal(h$Q[1], Q); expect_equal(h$df[1], 2)
  expect_equal(h$tau2[1], max(0, (Q - 2) / (sum(w) - sum(w^2) / sum(w))))
  expect_equal(h$Q[2], 0); expect_equal(h$tau2[2], 0); expect_equal(h$bs[, 2], rep(0.3, 3))
})

simr <- simulate_edgeicg(n = 5000, G = 160, seed = 11, region_gamma = 1.2)
fitr <- function(...) edgeicg(simr$counts, simr$coords, simr$types, lib = simr$lib, B = 8, min_shift = 150,
                              verbose = FALSE, region = simr$region, local = TRUE, patch_size = 300, ...)

test_that("region / local leave the global results unchanged", {
  f0 <- edgeicg(simr$counts, simr$coords, simr$types, lib = simr$lib, B = 8, min_shift = 150, verbose = FALSE)
  f1 <- fitr()
  expect_equal(f1$scale_fits, f0$scale_fits)
  expect_equal(f1$screen, f0$screen)
  expect_equal(f1$attribution, f0$attribution)
})

test_that("region and local results do not depend on windows or workers", {
  one <- fitr(backend = "r")
  runs <- list(fitr(window_cells = 900, backend = "r"))
  if (requireNamespace("mirai", quietly = TRUE) && !identical(Sys.getenv("NOT_CRAN"), "")) runs <- c(runs, list(fitr(window_cells = 900, workers = 2, backend = "r")))
  for (f in runs) {
    for (s in c("contact", "paracrine", "niche")) {
      expect_equal(f$region$effect[[s]], one$region$effect[[s]], tolerance = 1e-8)
      expect_equal(f$region$se[[s]], one$region$se[[s]], tolerance = 1e-8)
      expect_equal(f$local$het[[s]]$Q, one$local$het[[s]]$Q, tolerance = 1e-6)
    }
  }
})

test_that("a region-restricted effect is found in its region; the tables are well formed", {
  f <- fitr()
  R <- edgeicg_region_results(f)
  expect_true(all(c("region", "effect", "lfc", "se", "z", "q", "Q", "df", "q_het", "pattern") %in% names(R)))
  right <- simr$truth[region == "right", gene]
  hit <- unique(R[receiver == "t1" & sender == "t2" & scale == "contact" & gene %in% right & region == "right" & q < 0.05 & z > 0, gene])
  expect_gte(length(hit), 7)
  # a right-only effect is rarely read as opposite (the rule uses a nominal |z| >= 2 in the other region)
  expect_lte(mean(R[gene %in% hit & region == "right", pattern] == "opposite"), 0.2)
  L <- edgeicg_local_results(f)
  expect_true(all(c("local_frac", "q_het") %in% names(L$entries)))
  expect_true(all(c("patch", "x", "y", "effect_shrunk", "sd_shrunk") %in% names(L$patches)))
  expect_error(edgeicg_region_results(edgeicg(simr$counts, simr$coords, simr$types, lib = simr$lib, calibrate = FALSE,
                                              verbose = FALSE)), "no region effects")
})
