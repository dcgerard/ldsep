test_that("problem data for norm", {
  pgA <- as.matrix(read.csv("./pgA.csv"))
  pgB <- as.matrix(read.csv("./pgB.csv"))
  K <- ncol(pgA) - 1

  postA <- exp(pgA)
  postA <- postA / rowSums(postA)
  postB <- exp(pgB)
  postB <- postB / rowSums(postB)
  ega <- rowSums(sweep(x = postA, MARGIN = 2, STATS = 0:K, FUN = `*`))
  egb <- rowSums(sweep(x = postB, MARGIN = 2, STATS = 0:K, FUN = `*`))
  mu_init <- c(mean(ega), mean(egb))
  sigma_init <- stats::cov(cbind(ega, egb))
  L <- t(chol(x = sigma_init))
  par <- c(mu_init, L[lower.tri(L, diag = TRUE)])

  oout <- stats::optim(par = par,
                       fn = obj_pbnorm_genolike,
                       method = "L-BFGS-B",
                       lower = c(-Inf, -Inf, 0.01, -Inf, 0.01),
                       upper = rep(Inf, 5),
                       control = list(fnscale = -1),
                       hessian = TRUE,
                       pgA = pgA,
                       pgB = pgB)
})
