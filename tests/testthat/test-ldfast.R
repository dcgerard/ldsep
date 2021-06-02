test_that("ldfast skips monomorphic sites", {


})


test_that("mycor works", {
  x <- runif(10)
  y <- runif(10)
  x[3:4] <- NA_real_
  y[4:5] <- NA_real_

  expect_equal(
    cor(x, y, use = "pairwise.complete.obs"),
    mycor(x, y)
  )
})

test_that("slcor works", {
  n <- 10
  p <- 100
  xmat <- matrix(rnorm(n * p), ncol = n)
  xmat[sample(n * p, size = 30)] <- NA_real_

  cout <- cor(xmat, use = "pairwise.complete.obs")
  slout <- slcor(xmat)

  which_compare <- !is.na(slout)

  expect_equal(
    cout[which_compare],
    slout[which_compare]
  )
})
