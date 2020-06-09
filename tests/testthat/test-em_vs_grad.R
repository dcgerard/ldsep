context("EM and gradient ascent")

test_that("EM and gradient ascent give same results", {
  set.seed(1)
  n <- 100
  K <- 4
  ga <- stats::rbinom(n = n, size = K, prob = 0.5)
  gb <- stats::rbinom(n = n, size = K, prob = 0.5)
  pgA <- t(sapply(ga, stats::dnorm, x = 0:K, sd = 1, log = TRUE))
  pgB <- t(sapply(gb, stats::dnorm, x = 0:K, sd = 1, log = TRUE))
  pen <- 2

  pinit <- stats::runif(4)
  pinit <- pinit / sum(pinit)
  inity <- simplex_to_real(pinit)

  oout <- stats::optim(par     = inity,
                       fn      = llike_genolike,
                       gr      = dllike_genolike_dpar,
                       method  = "BFGS",
                       control = list(fnscale = -1),
                       pgA      = pgA,
                       pgB      = pgB,
                       alpha    = rep(pen, 4))
  emout <- genolike_em(p = pinit,
                       pgA = pgA,
                       pgB = pgB,
                       alpha = rep(pen, 4),
                       verbose = FALSE,
                       square = FALSE)

  p_grad <- real_to_simplex(oout$par)

  p_grad
  emout

  l1 <- llike_genolike(par = simplex_to_real(emout),
                 pgA = pgA,
                 pgB = pgB,
                 alpha = rep(pen, 4))
  l2 <- llike_genolike(par = oout$par,
                 pgA = pgA,
                 pgB = pgB,
                 alpha = rep(pen, 4))

  expect_equal(l1, l2, tolerance = 10e-5)

})
