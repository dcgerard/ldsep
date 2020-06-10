context("derivatives")

test_that("dmulti_dprob works", {
  x <- c(1, 2, 0, 5)
  prob <- stats::runif(4)
  prob <- prob / sum(prob)

  derivvec <- dmulti_dprob(x = x, prob = prob, log_p = FALSE)
  myenv <- new.env()
  assign(x = "x", value = x, envir = myenv)
  assign(x = "prob", value = prob, envir = myenv)
  nout <- stats::numericDeriv(quote(dmulti_double(x = x,
                                                  prob = prob,
                                                  log_p = FALSE)),
                              "prob", myenv)

  expect_equal(
    c(attr(nout, "gradient")),
    c(derivvec),
    tolerance = 10^-5
  )

  # microbenchmark::microbenchmark(
  #   stats::numericDeriv(quote(dmulti_double(x = x, prob = prob, log_p = FALSE)), "prob", myenv),
  #   dmulti_dprob(x = x, prob = prob, log_p = FALSE)
  # )

})

test_that("dprobgeno_dprob works", {

  gA <- 3
  gB <- 1
  K <- 4
  prob <- stats::runif(4)
  prob <- prob / sum(prob)
  derivvec <- dprobgeno_dprob(gA = gA,
                              gB = gB,
                              K = K,
                              prob = prob)

  myenv <- new.env()
  assign(x = "gA", value = gA, envir = myenv)
  assign(x = "gB", value = gB, envir = myenv)
  assign(x = "K", value = K, envir = myenv)
  assign(x = "prob", value = prob, envir = myenv)
  nout <- stats::numericDeriv(quote(probgeno(gA = gA,
                                             gB = gB,
                                             K = K,
                                             prob = prob,
                                             log_p = TRUE)),
                              "prob", myenv)

  expect_equal(
    c(attr(nout, "gradient")),
    c(derivvec),
    tolerance = 10^-5
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
  prob <- stats::runif(4)
  prob <- prob / sum(prob)
  derivvec <- dproballgeno_dprob(gA = gA, gB = gB, K = K, prob = prob)

  myenv <- new.env()
  assign(x = "gA", value = gA, envir = myenv)
  assign(x = "gB", value = gB, envir = myenv)
  assign(x = "K", value = K, envir = myenv)
  assign(x = "prob", value = prob, envir = myenv)
  nout <- stats::numericDeriv(quote(proballgeno(gA = gA,
                                                gB = gB,
                                                K = K,
                                                prob = prob,
                                                log_p = TRUE)), "prob", myenv)

  expect_equal(
    c(attr(nout, "gradient")),
    c(derivvec),
    tolerance = 10^-5
  )

  # microbenchmark::microbenchmark(
  #   stats::numericDeriv(quote(proballgeno(gA = gA, gB = gB, K = K, prob = prob, log_p = TRUE)), "prob", myenv),
  #   dproballgeno_dprob(gA = gA, gB = gB, K = K, prob = prob)
  # )

})


test_that("derivative of real_to_simplex is correct", {
  y <- c(1, -2, 3)

  jacobmat <- dreal_to_simplex_dy(y = y)

  fpi <- function(y, i) {
    real_to_simplex(y = y)[i]
  }

  numerjacob <- matrix(NA_real_, nrow = length(y) + 1, ncol = length(y))
  myenv <- new.env()
  assign(x = "y", value = y, envir = myenv)
  assign(x = "i", value = 1, envir = myenv)
  nout <- stats::numericDeriv(quote(fpi(y = y, i = i)), "y", myenv)
  numerjacob[1, ] <- c(attr(nout, "gradient"))
  assign(x = "i", value = 2, envir = myenv)
  nout <- stats::numericDeriv(quote(fpi(y = y, i = i)), "y", myenv)
  numerjacob[2, ] <- c(attr(nout, "gradient"))
  assign(x = "i", value = 3, envir = myenv)
  nout <- stats::numericDeriv(quote(fpi(y = y, i = i)), "y", myenv)
  numerjacob[3, ] <- c(attr(nout, "gradient"))
  assign(x = "i", value = 4, envir = myenv)
  nout <- stats::numericDeriv(quote(fpi(y = y, i = i)), "y", myenv)
  numerjacob[4, ] <- c(attr(nout, "gradient"))

  expect_equal(
    numerjacob,
    jacobmat,
    tolerance = 10^-5
  )
})

test_that("derivative of simplex_to_real is correct", {
  x <- c(0.1, 0.2, 0.3, 0.4)

  jacobmat <- dsimplex_to_real_dx(x = x)

  fpi <- function(x, i) {
    simplex_to_real(x = x)[i]
  }

  numerjacob <- matrix(NA_real_, nrow = length(x) - 1, ncol = length(x))
  myenv <- new.env()
  assign(x = "x", value = x, envir = myenv)
  assign(x = "i", value = 1, envir = myenv)
  nout <- stats::numericDeriv(quote(fpi(x = x, i = i)), "x", myenv)
  numerjacob[1, ] <- c(attr(nout, "gradient"))
  assign(x = "i", value = 2, envir = myenv)
  nout <- stats::numericDeriv(quote(fpi(x = x, i = i)), "x", myenv)
  numerjacob[2, ] <- c(attr(nout, "gradient"))
  assign(x = "i", value = 3, envir = myenv)
  nout <- stats::numericDeriv(quote(fpi(x = x, i = i)), "x", myenv)
  numerjacob[3, ] <- c(attr(nout, "gradient"))

  expect_equal(
    numerjacob,
    jacobmat,
    tolerance = 10^-5
  )
})

test_that("dllike_geno_dpar works", {

  gA <- rep(0:4, each = 5)
  gB <- rep(0:4, 5)
  K <- 4
  par <- c(-1, 2, 3)
  alpha <- 1:4
  derivvec <- dllike_geno_dpar(par = par,
                               gA = gA,
                               gB = gB,
                               K = K,
                               alpha = alpha)

  myenv <- new.env()
  assign(x = "gA", value = gA, envir = myenv)
  assign(x = "gB", value = gB, envir = myenv)
  assign(x = "K", value = K, envir = myenv)
  assign(x = "par", value = par, envir = myenv)
  assign(x = "alpha", value = alpha, envir = myenv)
  nout <- stats::numericDeriv(quote(llike_geno(gA = gA,
                                               gB = gB,
                                               K = K,
                                               par = par,
                                               alpha = alpha)), "par", myenv)

  expect_equal(
    c(attr(nout, "gradient")),
    c(derivvec),
    tolerance = 10^-5
  )

  # microbenchmark::microbenchmark(
  #   stats::numericDeriv(quote(llike_geno(gA = gA, gB = gB, K = K, par = par)), "par", myenv),
  #   dllike_geno_dpar(par = par, gA = gA, gB = gB, K = K)
  # )
})


test_that("dD_dprob works", {
  D <- function(prob) { ## just the first three probabilities
    stopifnot(length(prob) == 3)
    p4 <- 1 - sum(prob)
    p4 - (prob[[2]] + p4) * (prob[[3]] + p4)
  }

  prob <- stats::runif(4)
  prob <- prob / sum(prob)
  deriv <- dD_dprob(prob)
  myenv <- new.env()
  assign(x = "prob", value = prob[1:3], envir = myenv)
  nout <- stats::numericDeriv(quote(D(prob)), "prob", myenv)

  expect_equal(
    c(attr(nout, "gradient")),
    c(deriv),
    tolerance = 10^-5
  )
})

test_that("dr2_dprob works", {
  r2 <- function(prob) { ## just the first three probabilities
    stopifnot(length(prob) == 3)
    p4 <- 1 - sum(prob)
    pA <- prob[2] + p4
    pB <- prob[3] + p4
    D <- p4 - pA * pB
    D^2 / (pA * (1 - pA) * pB * (1 - pB))
  }

  prob <- stats::runif(4)
  prob <- prob / sum(prob)
  deriv <- dr2_dprob(prob)
  myenv <- new.env()
  assign(x = "prob", value = prob[1:3], envir = myenv)
  nout <- stats::numericDeriv(quote(r2(prob)), "prob", myenv)

  expect_equal(
    c(attr(nout, "gradient")),
    c(deriv),
    tolerance = 10^-5
  )
})

test_that("drDprime_dprob works", {
  Dprime <- function(prob) { ## just the first three probabilities
    stopifnot(length(prob) == 3)
    p4 <- 1 - sum(prob)
    pA <- prob[2] + p4
    pB <- prob[3] + p4
    D <- p4 - pA * pB
    if (D < 0) {
      Dprime <- D / min(pA * pB, (1 - pA) * (1 - pB))
    } else {
      Dprime <- D / min(pA * (1 - pB), (1 - pA) * pB)
    }
    return(Dprime)
  }

  prob <- stats::runif(4)
  prob <- prob / sum(prob)
  deriv <- dDprime_dprob(prob)
  myenv <- new.env()
  assign(x = "prob", value = prob[1:3], envir = myenv)
  nout <- stats::numericDeriv(quote(Dprime(prob)), "prob", myenv)

  expect_equal(
    c(attr(nout, "gradient")),
    c(deriv),
    tolerance = 10^-5
  )
})

test_that("dprobgenolike_dprob works OK", {
  pgA <- log(c(0.1, 0.1, 0.3, 0.1, 0.05))
  pgB <- log(c(0.01, 0.1, 0.3, 0.4, 0.1))
  prob <- stats::runif(4)
  prob <- prob / sum(prob)
  derivvec <- dprobgenolike_dprob(pgA = pgA, pgB = pgB, prob = prob)

  myenv <- new.env()
  assign(x = "pgA", value = pgA, envir = myenv)
  assign(x = "pgB", value = pgB, envir = myenv)
  assign(x = "prob", value = prob, envir = myenv)
  nout <- stats::numericDeriv(quote(probgenolike(pgA = pgA,
                                                 pgB = pgB,
                                                 prob = prob)), "prob", myenv)
  expect_equal(
    c(attr(nout, "gradient")),
    c(derivvec),
    tolerance = 10^-5
  )
})


test_that("dproballgenolike_dprob works OK", {
  n <- 10
  K <- 4
  pgA <- matrix(log(stats::runif(n * K)), nrow = n)
  pgB <- matrix(log(stats::runif(n * K)), nrow = n)
  prob <- stats::runif(4)
  prob <- prob / sum(prob)
  derivvec <- dproballgenolike_dprob(pgA = pgA, pgB = pgB, prob = prob)

  myenv <- new.env()
  assign(x = "pgA", value = pgA, envir = myenv)
  assign(x = "pgB", value = pgB, envir = myenv)
  assign(x = "prob", value = prob, envir = myenv)
  nout <- stats::numericDeriv(quote(proballgenolike(pgA = pgA,
                                                    pgB = pgB,
                                                    prob = prob)),
                              "prob", myenv)
  expect_equal(
    c(attr(nout, "gradient")),
    c(derivvec),
    tolerance = 10^-5
  )
})

test_that("dllike_genolike_dpar works OK", {
  n <- 10
  K <- 4
  pgA <- matrix(log(stats::runif(n * K)), nrow = n)
  pgB <- matrix(log(stats::runif(n * K)), nrow = n)
  par <- stats::rnorm(3)
  alpha <- 1:4
  derivvec <- dllike_genolike_dpar(par = par,
                                   pgA = pgA,
                                   pgB = pgB,
                                   alpha = alpha)

  myenv <- new.env()
  assign(x = "pgA", value = pgA, envir = myenv)
  assign(x = "pgB", value = pgB, envir = myenv)
  assign(x = "par", value = par, envir = myenv)
  assign(x = "alpha", value = alpha, envir = myenv)
  nout <- stats::numericDeriv(quote(llike_genolike(pgA = pgA,
                                                   pgB = pgB,
                                                   par = par,
                                                   alpha = alpha)),
                              "par", myenv)
  expect_equal(
    c(attr(nout, "gradient")),
    c(derivvec),
    tolerance = 10^-4
  )
})


test_that("dlprior_par_dprob works OK", {
  alpha <- 1:4
  par <- stats::rnorm(3)
  derivvec <- dlprior_par_dprob(par = par, alpha = alpha)

  myenv <- new.env()
  assign(x = "alpha", value = alpha, envir = myenv)
  assign(x = "par", value = par, envir = myenv)
  nout <- stats::numericDeriv(quote(lprior_par(par = par, alpha = alpha)),
                              "par", myenv)
  expect_equal(
    c(attr(nout, "gradient")),
    c(derivvec),
    tolerance = 10^-4
  )
})
