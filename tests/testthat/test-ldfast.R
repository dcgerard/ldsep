test_that("ldfast skips monomorphic sites", {
  set.seed(1)
  n <- 100
  p <- 10
  ploidy <- 4

  ## Generate random data
  gp <- array(stats::runif(n * p * (ploidy + 1)), dim = c(p, n, (ploidy + 1)))
  gp <- sweep(x = gp, MARGIN = c(1, 2), STATS = apply(gp, c(1, 2), sum), FUN = `/`)

  ## Add a couple monomorphic SNPs
  gp[4, , 1] <- 1
  gp[4, , 2] <- 0
  gp[4, , 3] <- 0
  gp[4, , 4] <- 0
  gp[4, , 5] <- 0

  gp[3, , 1] <- 0
  gp[3, , 2] <- 0
  gp[3, , 3] <- 0
  gp[3, , 4] <- 0
  gp[3, , 5] <- 1

  ## Check that everything only produces warnings (due to monomorphic SNPs).
  expect_warning(ldfast(gp = gp, type = "r", upper = 100))
  expect_warning(ldfast(gp = gp, type = "r2", upper = 100))
  expect_warning(ldfast(gp = gp, type = "z", upper = 100))
  expect_warning(ldfast(gp = gp, type = "D", upper = 100))
  expect_warning(ldfast(gp = gp, type = "Dprime", upper = 100))

  ## Check sliding window
  expect_error(ldfast(gp = gp, type = "r", win = 1, se = FALSE, thresh = FALSE, shrinkrr = FALSE), NA)
  expect_error(ldfast(gp = gp, type = "r2", win = 1, se = FALSE, thresh = FALSE, shrinkrr = FALSE), NA)
  expect_error(ldfast(gp = gp, type = "z", win = 1, se = FALSE, thresh = FALSE, shrinkrr = FALSE), NA)
  expect_error(ldfast(gp = gp, type = "D", win = 1, se = FALSE, thresh = FALSE, shrinkrr = FALSE), NA)
  expect_error(ldfast(gp = gp, type = "Dprime", win = 1, se = FALSE, thresh = FALSE, shrinkrr = FALSE), NA)

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
