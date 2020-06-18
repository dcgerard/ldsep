context("pbnorm")

test_that("pbnorm_dist sums to 1", {
  A <- matrix(stats::rnorm(4), nrow = 2)
  sigma <- crossprod(A) * 20
  mu <- c(1, 4)
  ploidy <- 6

  expect_equal(
    sum(pbnorm_dist(mu = mu, sigma = sigma, K = ploidy, log = FALSE)),
    1
  )
})

test_that("log-likelihoods for pbnorm work", {
  set.seed(1)
  n <- 100
  K <- 6
  ga <- stats::rbinom(n = n, size = K, prob = 0.5)
  gb <- stats::rbinom(n = n, size = K, prob = 0.5)
  pgA <- t(sapply(ga, stats::dnorm, x = 0:K, sd = 1, log = TRUE))
  pgB <- t(sapply(gb, stats::dnorm, x = 0:K, sd = 1, log = TRUE))

  ega <- rowSums(sweep(x = exp(pgA), MARGIN = 2, STATS = 0:K, FUN = `*`))
  egb <- rowSums(sweep(x = exp(pgB), MARGIN = 2, STATS = 0:K, FUN = `*`))
  mu_init <- c(mean(ega), mean(egb))
  sigma_init <- stats::cov(cbind(ega, egb))
  L <- t(chol(x = sigma_init))
  par <- c(mu_init, L[lower.tri(L, diag = TRUE)])

  lval1 <- llike_pbnorm_genolike(pgA = pgA,
                                 pgB = pgB,
                                 mu = mu_init,
                                 sigma = sigma_init) +
    prior_mu(mu = mu_init, K = K) +
    prior_sigma(lvec = L[lower.tri(L, diag = TRUE)])
  lval2 <- obj_pbnorm_genolike(par = par, pgA = pgA, pgB = pgB)
  expect_equal(lval1, lval2)
})
