context("convert probs to LD")

test_that("Conversion works", {
  phat <- c(0.1, 0.2, 0.3, 0.4)
  temp <- convert_ld(phat = phat)

  phat <- c(0, 0.5, 0, 0.5)
  temp <- convert_ld(phat = phat)
  expect_equal(temp[["D"]], 0)
})
