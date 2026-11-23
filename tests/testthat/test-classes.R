test_that("edgeicg_classify flags sender markers, other markers and the receiver's own genes going down", {
  # 3 types x 5 genes of mean composition: g1 marks A, g2 marks B, g3 marks C, g4 and g5 are shared
  tm <- rbind(A = c(10, 0.1, 0.1, 1, 1), B = c(0.1, 10, 0.1, 1, 1), C = c(0.1, 0.1, 10, 1, 1))
  x <- list(types = c("A", "B", "C"), genes = paste0("g", 1:5), type_mean = tm, fitted = 1:3)
  hits <- data.table::data.table(receiver = "A", sender = "B", scale = "contact",
                                 gene = c("g2", "g3", "g1", "g1", "g4"), z = c(5, 5, -5, 5, -5))
  h <- edgeicg_classify(x, hits)
  expect_identical(h$class, c("sender_marker", "other_marker", "own_down", "clean", "clean"))
  # user-supplied receiver markers also count as own genes
  h2 <- edgeicg_classify(x, hits, markers = list(A = "g4"))
  expect_identical(h2$class[5], "own_down")
})
