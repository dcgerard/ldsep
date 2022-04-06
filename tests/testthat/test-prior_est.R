test_that("llike_li() and llike_li_log() give same results", {
  data("glike")

  B <- glike[1, , ]
  A <- exp(B)

  pivec <- 1:ncol(A)
  pivec <- pivec / sum(pivec)
  lpivec <- log(pivec)

  expect_equal(llike_li(A = A, pivec = pivec),
               llike_li_log(B = B, lpivec = lpivec))

  piest1 <- em_li(A = A)
  piest2 <- em_li_log(B = B)

  expect_equal(piest1, exp(piest2))

  ## bench::mark(em_li(A = A), em_li_log(B = B), check = FALSE)
})
