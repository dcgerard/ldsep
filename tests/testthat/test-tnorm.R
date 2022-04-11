test_that("tnorm estimates are reasonable", {

  xbar <- 2.3
  s <- 1.2
  ploidy <- 6

  expect_equal(
    tnorm_est(xbar = xbar, s = s, ploidy = ploidy),
    tnorm_solve(mu_obs = xbar, sigma_obs = s, ploidy = ploidy),
    tolerance = 0.1
  )
  expect_equal(mom_tnorm(mu = 2.139, sigma = 1.370, a = 0, b = ploidy),
               c(2.3, 1.2),
               tolerance = 0.001)

  # bench::mark(
  #   tnorm_est(xbar = xbar, s = s, ploidy = ploidy),
  #   tnorm_solve(mu_obs = xbar, sigma_obs = s, ploidy = ploidy),
  #   check = FALSE
  # )

})
