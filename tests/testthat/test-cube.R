context("bulk calculations")

test_that("get_dprobgeno_dprob_array works", {

  prob <- stats::runif(4)
  prob <- prob / sum(prob)
  K <- 4
  darray <- exp(get_dprobgeno_dprob_array(K = K, prob = prob))

  numdarray <- array(NA_real_, dim = c(K + 1, K + 1, 4))
  for (gA in 0:K) {
    for (gB in 0:K) {
      myenv <- new.env()
      assign(x = "gA", value = gA, envir = myenv)
      assign(x = "gB", value = gB, envir = myenv)
      assign(x = "prob", value = prob, envir = myenv)
      assign(x = "K", value = K, envir = myenv)
      nout <- stats::numericDeriv(quote(probgeno(gA = gA,
                                                 gB = gB,
                                                 K = K,
                                                 prob = prob,
                                                 log_p = FALSE)),
                                  "prob", myenv)

      numdarray[gA + 1, gB + 1, ] <- c(attr(nout, "gradient"))
    }
  }
  expect_equal(darray, numdarray, tolerance = 10^-5)
})


test_that("old and new progallgenolike are the same", {
  K <- 6
  n <- 100
  ga <- sample(0:K, n, TRUE)
  gb <- sample(0:K, n, TRUE)
  pgA <- t(sapply(ga, stats::dnorm, x = 0:K, sd = 2, log = TRUE))
  pgB <- t(sapply(gb, stats::dnorm, x = 0:K, sd = 2, log = TRUE))
  prob <- stats::runif(4)
  prob <- prob / sum(prob)

  expect_equal(
    proballgenolike(pgA = pgA, pgB = pgB, prob = prob, log_p = TRUE),
    proballgenolike_old(pgA = pgA, pgB = pgB, prob = prob, log_p = TRUE)
  )

  # microbenchmark::microbenchmark(
  #   proballgenolike(pgA = pgA, pgB = pgB, prob = prob, log_p = TRUE),
  #   proballgenolike_old(pgA = pgA, pgB = pgB, prob = prob, log_p = TRUE)
  # )
})
