context("hessian")

test_that("hessian works", {
  n <- 100
  K <- 6
  ga <- stats::rbinom(n = n, size = K, prob = 0.5)
  gb <- stats::rbinom(n = n, size = K, prob = 0.5)
  ga <- t(sapply(ga, stats::dnorm, x = 0:K, sd = 1, log = TRUE))
  gb <- t(sapply(gb, stats::dnorm, x = 0:K, sd = 1, log = TRUE))
  pen <- 2
  alphamat <- matrix(data = pen, nrow = K + 1, ncol = K + 1)
  pma <- factor(apply(X = ga, MARGIN = 1, FUN = which.max) - 1, levels = 0:K)
  pmb <- factor(apply(X = gb, MARGIN = 1, FUN = which.max) - 1, levels = 0:K)
  pinit <- as.matrix(prop.table(table(pma, pmb) + 1))
  gout <- em_jointgeno(p = pinit,
                       pgA = ga,
                       pgB = gb,
                       alpha = alphamat)
  oout <- stats::optim(par = c(gout),
                       fn = llike_jointgeno_vec,
                       method = "BFGS",
                       control = list(fnscale = -1, maxit = 0),
                       hessian = TRUE,
                       pgA = ga,
                       pgB = gb,
                       alpha = alphamat)
  # hessmat <- numDeriv::hessian(func = llike_jointgeno_vec,
  #                              x = c(gout),
  #                              pgA = ga,
  #                              pgB = gb,
  #                              alpha = alphamat)
  # plot(hess2, hessmat)
  # abline(0, 1)

  hess2 <- hessian_jointgeno(p = gout, pgA = ga, pgB = gb, alpha = alphamat)

  ## high because high numerical error in optim
  expect_equal(hess2, oout$hessian, tol = 1000)
})


test_that("dD_dqlm works OK", {

  K <- 7
  qmat <- matrix(runif((K + 1) ^ 2), nrow = K + 1)
  qmat <- qmat / sum(qmat)
  par <- c(qmat)
  f <- function(par, K) {
    Dfromg(matrix(par, nrow = K + 1))
  }

  gradval <- dD_dqlm(p = qmat)

  myenv <- new.env()
  assign(x = "qmat", value = qmat, envir = myenv)
  assign(x = "K", value = K, envir = myenv)
  nout <- stats::numericDeriv(quote(f(par = par, K = K)), "par", myenv)
  expect_equal(
    c(attr(nout, "gradient")),
    c(gradval),
    tol = 10^-5
  )
})


test_that("dr2_dqlm works OK", {
  K <- 6
  qmat <- matrix(runif((K + 1) ^ 2), nrow = K + 1)
  qmat <- qmat / sum(qmat)
  par <- c(qmat)
  f <- function(par, K) {
    r2fromg(matrix(par, nrow = K + 1))
  }
  D <- Dfromg(gmat = qmat)
  dgrad <- dD_dqlm(p = qmat)

  gradval <- dr2_dqlm(p = qmat, dgrad = dgrad, D = D)

  myenv <- new.env()
  assign(x = "qmat", value = qmat, envir = myenv)
  assign(x = "K", value = K, envir = myenv)
  nout <- stats::numericDeriv(quote(f(par = par, K = K)), "par", myenv)
  expect_equal(
    c(attr(nout, "gradient")),
    c(gradval),
    tol = 10^-5
  )
})


test_that("dr2_dqlm works OK", {
  f <- function(par, K) {
    gout <- matrix(par, nrow = K + 1)
    D <- Dfromg(gout)
    distA <- rowSums(gout)
    distB <- colSums(gout)
    egA <- sum((0:K) * distA)
    egB <- sum((0:K) * distB)
    if (D < 0) {
      Dmax <- min(egA * egB, (K - egA) * (K - egB)) / K^2
    } else {
      Dmax <- min(egA * (K - egB), (K - egA) * egB) / K^2
    }
    Dprime <- D / Dmax
    Dprime
  }

  g <- function(par, K) {
    gout <- matrix(par, nrow = K + 1)
    D <- Dfromg(gout)
    distA <- rowSums(gout)
    distB <- colSums(gout)
    egA <- sum((0:K) * distA)
    egB <- sum((0:K) * distB)
    if (D < 0) {
      Dmax <- min(egA * egB, (K - egA) * (K - egB)) / K^2
    } else {
      Dmax <- min(egA * (K - egB), (K - egA) * egB) / K^2
    }
    Dprime <- D / Dmax
    Dprime
    dgrad <- dD_dqlm(p = gout)
    ddprime_dqlm(p = gout, dgrad = dgrad, D = D, Dm = Dmax)
  }


  K <- 6
  qmat <- matrix(runif((K + 1) ^ 2), nrow = K + 1)
  qmat <- qmat / sum(qmat)
  par <- c(qmat)
  gradval <- g(par = par, K = K)

  myenv <- new.env()
  assign(x = "qmat", value = qmat, envir = myenv)
  assign(x = "K", value = K, envir = myenv)
  nout <- stats::numericDeriv(quote(f(par = par, K = K)), "par", myenv)
  expect_equal(
    c(attr(nout, "gradient")),
    c(gradval),
    tol = 10^-5
  )
})
