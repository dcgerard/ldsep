context("real vs simplex")

test_that("converting from real_to_simplex works", {
  expect_equal(
    sum(real_to_simplex(y = c(4, 4, 1, -10))),
    1
  )
})

test_that("real_to_simplex and simplex_to_real are inverses", {
  x <- c(0.1, 0.6, 0.2, 0.1)
  expect_equal(
    x,
    c(real_to_simplex(simplex_to_real(x = x)))
  )

  y <- c(-1, 0, 2)
  expect_equal(
    y,
    c(simplex_to_real(real_to_simplex(y = y)))
  )
})

test_that("proballgeno and llike_geno give same results", {
  par <- c(-1, 2, 4)
  ga <- c(1, 4, 1, 1)
  gb <- c(2, 4, 2, 5)
  K <- 6
  alpha <- 1:4

  expect_equal(
    llike_geno(par = par,
               gA = ga,
               gB = gb,
               K = K,
               alpha = alpha) -
      lprior_par(par = par,
                 alpha = alpha),
    proballgeno(gA = ga,
                gB = gb,
                K = K,
                prob = real_to_simplex(par),
                log_p = TRUE)
  )

})
