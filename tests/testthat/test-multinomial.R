context("multinomial")

test_that("dmulti", {
  x <- c(2, 4, 1, 0)
  pvec <- c(0.3, 0.3, 0.4, 0)
  expect_equal(
    stats::dmultinom(x = x, prob = pvec, log = TRUE),
    dmulti_double(x = x, prob = pvec, log_p = TRUE)
  )

  x <- c(2, 4, 1, 1)
  pvec <- c(0.3, 0.3, 0.4, 0)
  expect_equal(
    stats::dmultinom(x = x, prob = pvec, log = TRUE),
    dmulti_double(x = x, prob = pvec, log_p = TRUE)
  )

  x <- c(2, 4, 1, 0)
  pvec <- c(0.3, 0.3, 0.3, 0.1)
  expect_equal(
    stats::dmultinom(x = x, prob = pvec, log = TRUE),
    dmulti_double(x = x, prob = pvec, log_p = TRUE)
  )

  x <- c(2, 4, 1, 1)
  pvec <- c(0.3, 0.3, 0.3, 0.1)
  expect_equal(
    stats::dmultinom(x = x, prob = pvec, log = TRUE),
    dmulti_double(x = x, prob = pvec, log_p = TRUE)
  )

  # microbenchmark::microbenchmark(
  #   stats::dmultinom(x = x, prob = pvec, log = TRUE),
  #   dmulti_double(x = x, prob = pvec, log_p = TRUE)
  # )

})

test_that("probgeno works", {
  prob <- c(0.3, 0.3, 0.2, 0.2)
  expect_equal(
    probgeno(gA = 1, gB = 1, K = 4, prob = prob, log_p = TRUE),
    log_sum_exp_2(dmulti_double(x = c(2, 1, 1, 0), prob = prob, log_p = TRUE),
                  dmulti_double(x = c(3, 0, 0, 1), prob = prob, log_p = TRUE))
  )
})
