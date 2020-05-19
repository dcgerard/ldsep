context("bulk calculations")

test_that("get_dprobgeno_dprob_array works", {

  prob <- stats::runif(4)
  prob <- prob / sum(prob)
  K <- 4
  darray <- exp(get_dprobgeno_dprob_array(K = K, prob = prob))

  numdarray <- array(NA_real_, dim = c(K + 1, K + 1, 4))
  for (gA in 0:K) {
    for(gB in 0:K) {
      myenv <- new.env()
      assign(x = "gA", value = gA, envir = myenv)
      assign(x = "gB", value = gB, envir = myenv)
      assign(x = "prob", value = prob, envir = myenv)
      assign(x = "K", value = K, envir = myenv)
      nout <- stats::numericDeriv(quote(probgeno(gA = gA, gB = gB, K = K, prob = prob, log_p = FALSE)), "prob", myenv)

      numdarray[gA + 1, gB + 1, ] <- c(attr(nout, "gradient"))
    }
  }
  expect_equal(darray, numdarray, tolerance = 10^-5)
})
