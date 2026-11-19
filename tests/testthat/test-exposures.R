test_that("exposures are shares", {
  set.seed(1); co <- matrix(runif(1000, 0, 300), 500); ty <- sample(1:3, 500, TRUE)
  ex <- edgeicg_exposures(co, ty, 3)
  for (m in list(ex$C, ex$P, ex$N)) expect_equal(unname(rowSums(m))[rowSums(m) > 0], rep(1, sum(rowSums(m) > 0)))
  expect_equal(dim(ex$N), c(500L, 3L))
})

test_that("contact_from_polygons finds shared edges of touching squares", {
  sq <- function(cell, x0) data.frame(cell = cell, x = x0 + c(0, 10, 10, 0), y = c(0, 0, 10, 10))
  dens <- function(p) do.call(rbind, lapply(split(p, p$cell), function(d) {   # densify to 1 um vertices
    k <- nrow(d); do.call(rbind, lapply(seq_len(k), function(i) { j <- i %% k + 1; s <- (0:9) / 10
      data.frame(cell = d$cell[1], x = d$x[i] + s * (d$x[j] - d$x[i]), y = d$y[i] + s * (d$y[j] - d$y[i])) })) }))
  p <- dens(rbind(sq(1, 0), sq(2, 10), sq(3, 50)))
  e <- contact_from_polygons(p, 3, tau = 0.5)
  expect_setequal(paste(e$recv, e$send), c("1 2", "2 1"))
  expect_true(all(abs(e$w - 11) < 1.5))                               # ~10 um shared side (+ corners)
})
