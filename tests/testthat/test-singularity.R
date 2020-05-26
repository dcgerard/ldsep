context("singularity")

test_that("mle works on singularity", {
  load("./badvar.RData")
  ldout <- ldest(ga = ga, gb = gb, K = K, lang = "R")
  expect_equal(ldout[["D"]], 0)
})
