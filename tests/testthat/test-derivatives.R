context("derivatives")

test_that("dmulti_dprob works", {
  x <- c(1, 2, 0, 5)
  prob <- c(0.1, 0.2, 0.3, 0.4)

  derivvec <- dmulti_dprob(x = x, prob = prob, log_p = FALSE)
  myenv <- new.env()
  assign(x = "x", value = x, envir = myenv)
  assign(x = "prob", value = prob, envir = myenv)
  nout <- stats::numericDeriv(quote(dmulti_double(x = x, prob = prob, log_p = FALSE)), "prob", myenv)

  expect_equal(
    c(attr(nout, "gradient")),
    derivvec
  )

  # microbenchmark::microbenchmark(
  #   stats::numericDeriv(quote(dmulti_double(x = x, prob = prob, log_p = FALSE)), "prob", myenv),
  #   dmulti_dprob(x = x, prob = prob, log_p = FALSE)
  # )

})

test_that("dprobgeno_dprob works", {

  gA <- 2
  gB <- 0
  K <- 4
  prob <- c(0.05, 0.25, 0.3, 0.4)
  derivvec <- dprobgeno_dprob(gA = gA, gB = gB, K = K, prob = prob)

  myenv <- new.env()
  assign(x = "gA", value = gA, envir = myenv)
  assign(x = "gB", value = gB, envir = myenv)
  assign(x = "K", value = K, envir = myenv)
  assign(x = "prob", value = prob, envir = myenv)
  nout <- stats::numericDeriv(quote(probgeno(gA = gA, gB = gB, K = K, prob = prob, log_p = TRUE)), "prob", myenv)

  expect_equal(
    c(attr(nout, "gradient")),
    derivvec
  )

  # microbenchmark::microbenchmark(
  #   stats::numericDeriv(quote(probgeno(gA = gA, gB = gB, K = K, prob = prob, log_p = TRUE)), "prob", myenv),
  #   dprobgeno_dprob(gA = gA, gB = gB, K = K, prob = prob)
  # )

})

test_that("dproballgeno_dprob works", {

  gA <- rep(0:4, each = 5)
  gB <- rep(0:4, 5)
  K <- 4
  prob <- c(0.1, 0.2, 0.3, 0.4)
  derivvec <- dproballgeno_dprob(gA = gA, gB = gB, K = K, prob = prob)

  myenv <- new.env()
  assign(x = "gA", value = gA, envir = myenv)
  assign(x = "gB", value = gB, envir = myenv)
  assign(x = "K", value = K, envir = myenv)
  assign(x = "prob", value = prob, envir = myenv)
  nout <- stats::numericDeriv(quote(proballgeno(gA = gA, gB = gB, K = K, prob = prob, log_p = TRUE)), "prob", myenv)

  expect_equal(
    c(attr(nout, "gradient")),
    derivvec,
    tolerance = 10^-5
  )

  # microbenchmark::microbenchmark(
  #   stats::numericDeriv(quote(proballgeno(gA = gA, gB = gB, K = K, prob = prob, log_p = TRUE)), "prob", myenv),
  #   dproballgeno_dprob(gA = gA, gB = gB, K = K, prob = prob)
  # )

})


test_that("derivative of real_to_simplex is correct", {
  y <- c(1, -2, 3)
  real_to_simplex(y = y)

  jacobmat <- dreal_to_simplex_dy(y = y)

  fpi <- function(y, i) {
    real_to_simplex(y = y)[i]
  }

  numerjacob <- matrix(NA_real_, nrow = length(y) + 1, ncol = length(y))
  myenv <- new.env()
  assign(x = "y", value = y, envir = myenv)
  assign(x = "i", value = 1, envir = myenv)
  nout <- stats::numericDeriv(quote(fpi(y =y, i = i)), "y", myenv)
  numerjacob[1, ] <- c(attr(nout, "gradient"))
  assign(x = "i", value = 2, envir = myenv)
  nout <- stats::numericDeriv(quote(fpi(y =y, i = i)), "y", myenv)
  numerjacob[2, ] <- c(attr(nout, "gradient"))
  assign(x = "i", value = 3, envir = myenv)
  nout <- stats::numericDeriv(quote(fpi(y =y, i = i)), "y", myenv)
  numerjacob[3, ] <- c(attr(nout, "gradient"))
  assign(x = "i", value = 4, envir = myenv)
  nout <- stats::numericDeriv(quote(fpi(y =y, i = i)), "y", myenv)
  numerjacob[4, ] <- c(attr(nout, "gradient"))

  expect_equal(
    numerjacob,
    jacobmat
  )
})

test_that("dllike_geno_dpar works", {

  gA <- rep(0:4, each = 5)
  gB <- rep(0:4, 5)
  K <- 4
  par <- c(-1, 2, 3)
  derivvec <- dllike_geno_dpar(par = par, gA = gA, gB = gB, K = K)

  myenv <- new.env()
  assign(x = "gA", value = gA, envir = myenv)
  assign(x = "gB", value = gB, envir = myenv)
  assign(x = "K", value = K, envir = myenv)
  assign(x = "par", value = par, envir = myenv)
  nout <- stats::numericDeriv(quote(llike_geno(gA = gA, gB = gB, K = K, par = par)), "par", myenv)

  expect_equal(
    c(attr(nout, "gradient")),
    derivvec,
    tolerance = 10^-5
  )

  # microbenchmark::microbenchmark(
  #   stats::numericDeriv(quote(llike_geno(gA = gA, gB = gB, K = K, par = par)), "par", myenv),
  #   dllike_geno_dpar(par = par, gA = gA, gB = gB, K = K)
  # )
})
