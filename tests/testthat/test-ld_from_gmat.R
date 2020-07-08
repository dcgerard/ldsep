test_that("ld_from_gmat works", {
  gmat <- readRDS("./gmat.RDS")

  ## use the function
  ldout <- ld_from_gmat(gmat = gmat)

  ## do it manually
  K <- ncol(gmat) - 1
  distA <- rowSums(gmat)
  distB <- colSums(gmat)
  egA <- sum((0:K) * distA)
  egB <- sum((0:K) * distB)
  if (ldout[["D"]] < 0) {
    Dmax <- min(egA * egB, (K - egA) * (K - egB)) / K^2
  } else {
    Dmax <- min(egA * (K - egB), (K - egA) * egB) / K^2
  }
  Dprime <- ldout[["D"]] / Dmax

  expect_equal(Dprime, ldout[["Dprime"]])
  expect_equal(Dmax, ldout[["Dmax"]])
})
