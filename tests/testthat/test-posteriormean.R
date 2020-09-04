context("Posterior moments")

test_that("gl_to_gp works", {
  data("glike", package = "ldsep")
  gp2 <- ldsep::gl_to_gp(glike)
  expect_equal(dim(glike), dim(gp2))
  expect_true(all(abs(apply(X = gp2, MARGIN = c(1, 2), FUN = sum) - 1) < 10^-5))
  expect_true(all(apply(X = exp(glike) / gp2, MARGIN = c(1, 2), FUN = var) < 10^-5))
})

test_that("ldfast versions are the same", {
  data("gp", package = "ldsep")

  c1 <- ldfast_justmean(gp = gp, type = "r", shrinkrr = FALSE)
  c2 <- ldfast(gp = gp, type = "r", se = TRUE)
  expect_equal(c1$ldmat, c2$ldmat, tolerance = 10^-5)
  expect_equal(c1$rr, c2$rr, tolerance = 10^-5)

  c1 <- ldfast_justmean(gp = gp, type = "D", shrinkrr = FALSE)
  c2 <- ldfast(gp = gp, type = "D", se = TRUE)
  expect_equal(c1$ldmat, c2$ldmat, tolerance = 10^-5)
  expect_equal(c1$rr, c2$rr, tolerance = 10^-5)

  c1 <- ldfast_justmean(gp = gp, type = "Dprime", shrinkrr = FALSE)
  c2 <- ldfast(gp = gp, type = "Dprime", se = TRUE)
  expect_equal(c1$ldmat, c2$ldmat, tolerance = 10^-5)
  expect_equal(c1$rr, c2$rr, tolerance = 10^-5)

  c1 <- ldfast_justmean(gp = gp, type = "z", shrinkrr = FALSE)
  c2 <- ldfast(gp = gp, type = "z", se = TRUE)
  expect_equal(c1$ldmat, c2$ldmat, tolerance = 10^-5)
  expect_equal(c1$rr, c2$rr, tolerance = 10^-5)

  c1 <- ldfast_justmean(gp = gp, type = "r2", shrinkrr = FALSE)
  c2 <- ldfast(gp = gp, type = "r2", se = TRUE)
  expect_equal(c1$ldmat, c2$ldmat, tolerance = 10^-5)
  expect_equal(c1$rr, c2$rr, tolerance = 10^-5)

  # microbenchmark::microbenchmark(
  #   c1 <- ldfast_justmean(gp = gp, type = "r"),
  #   c2 <- ldfast(gp = gp, type = "r")
  # )
})

test_that("NA's are not propogated", {
  data("gp", package = "ldsep")

  gp[3, 1:10, 1] <- NA

  c1 <- ldfast_justmean(gp = gp, type = "r")
  c2 <- ldfast(gp = gp, type = "r", se = TRUE)
  expect_true(all(!is.na(c2$ldmat[upper.tri(c2$ldmat)])))
  expect_true(all(!is.na(c1$ldmat[upper.tri(c1$ldmat)])))
})

test_that("gradient for delta based on m works", {
  delta_from_m <- function(m, k) {
    stopifnot(length(m) == 7)
    stopifnot(length(k) == 1)
    ((m[6] + m[2] - m[1]^2) / (m[2] - m[1]^2)) *
      ((m[7] + m[4] - m[3]^2) / (m[4] - m[3]^2)) *
      ((m[5] - m[1] * m[3]) / k)
  }

  set.seed(1)
  m <- runif(7)
  k <- 6

  myenv <- new.env()
  assign(x = "m", value = m, envir = myenv)
  assign(x = "k", value = k, envir = myenv)
  nout <- stats::numericDeriv(quote(delta_from_m(m = m, k = k)), "m", myenv)
  attr(x = nout, which = "gradient")

  grad <- rep(NA_real_, 7)
  grad_delta_m(M = m, grad = grad, pd = k)
  grad

  # microbenchmark::microbenchmark(
  #   nout <- stats::numericDeriv(quote(delta_from_m(m = m, k = k)), "m", myenv),
  #   grad_delta_m(M = m, grad = grad, pd = k)
  # )

  expect_equal(c(grad), c(attr(x = nout, which = "gradient")), tolerance = 10^-5)
})


test_that("gradient for delta prime based on m works", {
  deltaprime_from_m <- function(m, k) {
    stopifnot(length(m) == 7)
    stopifnot(length(k) == 1)
    delta <- ((m[6] + m[2] - m[1]^2) / (m[2] - m[1]^2)) *
      ((m[7] + m[4] - m[3]^2) / (m[4] - m[3]^2)) *
      ((m[5] - m[1] * m[3]) / k)
    if (m[5] < m[1] * m[3]) {
      delta_m <- min(m[1] * m[3], (k - m[1]) * (k - m[3])) / k^2
    } else {
      delta_m <- min(m[1] * (k - m[3]), (k - m[1]) * m[3]) / k^2
    }
    return(delta / delta_m)
  }

  set.seed(1)

  for (index in 1:50) {
    m <- runif(7)
    k <- 6

    myenv <- new.env()
    assign(x = "m", value = m, envir = myenv)
    assign(x = "k", value = k, envir = myenv)
    nout <- stats::numericDeriv(quote(deltaprime_from_m(m = m, k = k)), "m", myenv)

    grad <- rep(NA_real_, 7)
    grad_deltaprime_m(M = m, grad = grad, pd = k)

    expect_equal(c(grad), c(attr(x = nout, which = "gradient")), tolerance = 10^-5)
  }

  # microbenchmark::microbenchmark(
  #   nout <- stats::numericDeriv(quote(deltaprime_from_m(m = m, k = k)), "m", myenv),
  #   grad_deltaprime_m(M = m, grad = grad, pd = k)
  # )
})


test_that("gradient for rho based on m works", {
  rho_from_m <- function(m) {
    stopifnot(length(m) == 7)
    (sqrt(m[6] + m[2] - m[1]^2) / (m[2] - m[1]^2)) *
      (sqrt(m[7] + m[4] - m[3]^2) / (m[4] - m[3]^2)) *
      (m[5] - m[1] * m[3])
  }

  set.seed(1)
  m <- runif(7)

  myenv <- new.env()
  assign(x = "m", value = m, envir = myenv)
  nout <- stats::numericDeriv(quote(rho_from_m(m = m)), "m", myenv)

  grad <- rep(NA_real_, 7)
  grad_rho_m(M = m, grad = grad)

  # microbenchmark::microbenchmark(
  #   nout <- stats::numericDeriv(quote(rho_from_m(m = m)), "m", myenv),
  #   grad_rho_m(M = m, grad = grad)
  # )

  expect_equal(c(grad), c(attr(x = nout, which = "gradient")), tolerance = 10^-5)
})

test_that("posterior moments are calculated correctly", {
  data("gp")
  nsnp <- dim(gp)[[1]]
  nind <- dim(gp)[[2]]
  ploidy <- dim(gp)[[3]] - 1

  pm_mat <- apply(gp, c(1, 2), function(x) sum((0:ploidy) * x))
  pv_mat <- apply(gp * outer(X = pm_mat, Y = 0:ploidy, FUN = `-`) ^ 2, c(1, 2), sum)

  temp1 <- matrix(NA_real_, nrow = nsnp, ncol = nind)
  temp2 <- matrix(NA_real_, nrow = nsnp, ncol = nind)

  fill_pm(pm = temp1, gp = gp)
  fill_pv(pv = temp2, pm = temp1, gp = gp)

  expect_true(max(abs(pv_mat - temp2)) < 10^-5)
  expect_true(max(abs(pm_mat - temp1)) < 10^-5)


  # microbenchmark::microbenchmark(
  #   {
  #     pm_mat <- apply(gp, c(1, 2), function(x) sum((0:ploidy) * x))
  #     pv_mat <- apply(gp * outer(X = pm_mat, Y = 0:ploidy, FUN = `-`) ^ 2, c(1, 2), sum)
  #   },
  #   {
  #     temp1 <- matrix(NA_real_, nrow = nsnp, ncol = nind)
  #     temp2 <- matrix(NA_real_, nrow = nsnp, ncol = nind)
  #     fill_pm(pm = temp1, gp = gp)
  #     fill_pv(pv = temp2, pm = temp1, gp = gp)
  #   }
  # )

})

test_that("NA ses are returned at threshholds in ldfast", {
  data("gp")

  ldout <- ldfast(gp = gp, type = "r", se = TRUE)
  expect_true(is.na(ldout$se[1, 1]))

  ldout <- ldfast(gp = gp, type = "D", se = TRUE)
  expect_true(is.na(ldout$se[1, 1]))

  ldout <- ldfast(gp = gp, type = "Dprime", se = TRUE)
  expect_true(is.na(ldout$se[1, 1]))

  ldout <- ldfast(gp = gp, type = "z", se = TRUE)
  expect_true(is.na(ldout$se[1, 1]))

  ldout <- ldfast(gp = gp, type = "r2", se = TRUE)
  expect_true(is.na(ldout$se[1, 1]))
})
