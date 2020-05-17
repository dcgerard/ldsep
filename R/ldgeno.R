###################
## Functions to estimate LD directly from genotypes
###################

#' Estimate LD directly from the genotypes
#'
#' @param ga A vector of counts, containing the genotypes for each individual
#'     at the first locus.
#' @param gb vector of counts, containing the genotypes for each individual
#'     at the second locus.
#' @param K the ploidy of the species.
#'
#' @author David Gerard
#'
#' @examples
#' ga <- c(1, 4, 1, 2)
#' gb <- c(0, 3, 3, 0)
#' K <- 4
#'
#' @export
ldest_geno <- function(ga, gb, K) {
  stopifnot(length(ga) == length(gb))
  stopifnot(ga >= 0, ga <= K)
  stopifnot(gb >= 0, gb <= K)
  stopifnot(length(K) == 1)

  inity <- rep(0, 3)

  oout <- stats::optim(par     = inity,
                       fn      = llike_geno,
                       gr      = dllike_geno_dpar,
                       method  = "BFGS",
                       control = list(fnscale = -1),
                       gA      = ga,
                       gB      = gb,
                       K       = K)

  # env <- new.env()
  # env[["gA"]] <- ga
  # env[["gB"]] <- gb
  # env[["K"]] <- K
  # oout <- lbfgs::lbfgs(call_eval = llike_geno,
  #                      call_grad = dllike_geno_dpar,
  #                      vars = inity,
  #                      gA = ga,
  #                      gB = gb,
  #                      K = K)

  phat <- real_to_simplex(oout$par) # (ab, Ab, aB, AB)

  pA <- phat[[2]] + phat[[4]]
  pB <- phat[[3]] + phat[[4]]
  D  <- phat[[4]] - pA * pB
  r2 <- D ^ 2 / (pA * (1 - pA) * pB * (1 - pB))

  if (D < 0) {
    Dprime <- min(pA * pB, (1 - pA) * (1 - pB))
  } else {
    Dprime <- min(pA * (1 - pB), (1 - pA) * pB)
  }

  retvec <- c(D = D, Dprime = Dprime, r2 = r2)

  return(retvec)
}
