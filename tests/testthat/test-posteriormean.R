context("Posterior moments")

test_that("Posterior mean from gp works", {
  data("gp", package = "ldsep")
  ds <- matrix(data = 0, nrow = dim(gp)[[2]], ncol = dim(gp)[[1]])

  ds_from_gp(ds = ds, gp = gp)
  pmout <- apply(sweep(x = gp, MARGIN = 3, STATS = 0:4, FUN = `*`), c(1, 2), sum)
  dimnames(pmout) <- NULL
  expect_equal(ds, t(pmout))

  pmout2 <- apply(sweep(x = gp, MARGIN = 3, STATS = (0:4)^2, FUN = `*`), c(1, 2), sum)
  vout <- pmout2 - pmout^2
  dimnames(vout) <- NULL

  pv <- rep(0, length = dim(gp)[[2]])
  pv_from_gp(pv = pv, gp = gp, ds = ds, i = 0)
  expect_equal(pv, vout[1, ])

  pv <- rep(0, length = dim(gp)[[2]])
  pv_from_gp(pv = pv, gp = gp, ds = ds, i = 3)
  expect_equal(pv, vout[4, ])

  c0 <- cor(t(pmout))
  c1 <- cov2cor(ldfast_post(gp = gp, type = "D")$cor)
  c2 <- ldfast_post(gp = gp, type = "r")$cor
  expect_equal(c1, c2)
})

test_that("gl_to_gp works", {
  data("glike", package = "ldsep")
  gp2 <- ldsep::gl_to_gp(glike)
  expect_equal(dim(glike), dim(gp2))
  expect_true(all(abs(apply(X = gp2, MARGIN = c(1, 2), FUN = sum) - 1) < 10^-5))
  expect_true(all(apply(X = exp(glike) / gp2, MARGIN = c(1, 2), FUN = var) < 10^-5))
})

test_that("ldfast and ldfast_c are same", {
  data("gp", package = "ldsep")
  c1 <- ldfast(gp)
  c2 <- ldfast_c(gp)
  dimnames(c1$cor) <- NULL
  dimnames(c2$cor) <- NULL
  expect_equal(c1$cor, c2$cor, tolerance = 10^-5)
})
