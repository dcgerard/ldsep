context("r vs c optimization")

test_that("R and C++ optimization give same results", {
  set.seed(1)
  inity <- c(0, 0, 0)
  ga <- sample(0:6, 100, TRUE)
  gb <- sample(0:6, 100, TRUE)
  K <- 6

  oout1 <- stats::optim(par     = inity,
                        fn      = llike_geno,
                        gr      = dllike_geno_dpar,
                        method  = "BFGS",
                        control = list(fnscale = -1),
                        hessian = TRUE,
                        gA      = ga,
                        gB      = gb,
                        K       = K)

  oout2 <- optimize_genocor(par = inity, gA = ga, gB = gb, K = K)

  expect_equal(oout1$par, c(oout2$par))
  expect_equal(oout1$value, c(oout2$value))
  expect_equal(oout1$hessian, oout2$hessian)

  # microbenchmark::microbenchmark(
  #   oout1 <- stats::optim(par     = inity,
  #                         fn      = llike_geno,
  #                         gr      = dllike_geno_dpar,
  #                         method  = "BFGS",
  #                         control = list(fnscale = -1),
  #                         hessian = TRUE,
  #                         gA      = ga,
  #                         gB      = gb,
  #                         K       = K),
  #   oout2 <- optimize_genocor(par = inity, gA = ga, gB = gb, K = K)
  # )

})

test_that("R and C++ optimization using genotype likelihoods is the same", {
  K <- 6
  n <- 100
  ga <- sample(0:K, n, TRUE)
  gb <- sample(0:K, n, TRUE)
  pgA <- t(sapply(ga, stats::dnorm, x = 0:K, sd = 2, log = TRUE))
  pgB <- t(sapply(gb, stats::dnorm, x = 0:K, sd = 2, log = TRUE))
  inity <- stats::rnorm(3)

  oout1 <- stats::optim(par     = inity,
                        fn      = llike_genolike,
                        gr      = dllike_genolike_dpar,
                        method  = "BFGS",
                        control = list(fnscale = -1),
                        hessian = TRUE,
                        pgA      = pgA,
                        pgB      = pgB)


  oout2 <- optimize_genolikecor(par = inity, pgA = pgA, pgB = pgB)

  expect_equal(oout1$par, c(oout2$par))
  expect_equal(oout1$value, c(oout2$value))
  expect_equal(oout1$hessian, oout2$hessian)

  # microbenchmark::microbenchmark(
  #   oout1 <- stats::optim(par     = inity,
  #                         fn      = llike_genolike,
  #                         gr      = dllike_genolike_dpar,
  #                         method  = "BFGS",
  #                         control = list(fnscale = -1),
  #                         hessian = TRUE,
  #                         pgA      = pgA,
  #                         pgB      = pgB),
  #   oout2 <- optimize_genolikecor(par = inity, pgA = pgA, pgB = pgB)
  # )

})
