test_that("maxq_from_marg works", {
  set.seed(1)
  K <- 6
  pA <- stats::runif(K + 1)
  pA <- pA / sum(pA)
  pB <- stats::runif(K + 1)
  pB <- pB / sum(pB)

  qmat <- maxq_from_marg(pA, pB, dir = "max")
  expect_equal(colSums(qmat), pB)
  expect_equal(rowSums(qmat), pA)
  expect_equal(sum(qmat), 1)

  qmat <- maxq_from_marg(pA, pB, dir = "min")
  expect_equal(colSums(qmat), pB)
  expect_equal(rowSums(qmat), pA)
  expect_equal(sum(qmat), 1)
})

test_that("maxq_from_allelef works", {
  set.seed(1)
  K <- 6
  pA <- stats::runif(1)
  pB <- stats::runif(1)

  qmat <- maxq_from_allelef(pA = pA, pB = pB, K = K, dir = "max")
  expect_equal(sum(qmat), 1)
  expect_equal(sum(colSums(qmat) * 0:K) / K, pB)
  expect_equal(sum(rowSums(qmat) * 0:K) / K, pA)

  qmat <- maxq_from_allelef(pA = pA, pB = pB, K = K, dir = "min")
  expect_equal(sum(qmat), 1)
  expect_equal(sum(colSums(qmat) * 0:K) / K, pB)
  expect_equal(sum(rowSums(qmat) * 0:K) / K, pA)
})

test_that("Dprime using moments is same as optimization", {
  K <- 6
  qmat <- matrix(stats::runif((K+1)^2), nrow = K+1)
  qmat <- qmat / sum(qmat)
  dfun1 <- Dprime(qmat, type = "allele", constrain = FALSE)
  dfun2 <- Dprime(qmat, type = "allele", constrain = TRUE)

  ## manual way
  pA <- sum(rowSums(qmat) * 0:K) / K
  pB <- sum(colSums(qmat) * 0:K) / K
  D <- Dfromg(gmat = qmat)
  dir <- ifelse(D < 0, "min", "max")
  qmax <- maxq_from_allelef(pA = pA, pB = pB, K = K, dir = dir)
  Dmax <- Dfromg(gmat = qmax)

  expect_equal(Dmax / dfun1[["Dmax"]], K * sign(D))
  expect_equal(Dmax, dfun2[["Dmax"]] * sign(D))
})
