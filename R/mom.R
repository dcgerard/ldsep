########################
## Method of moments estimators
########################

#' Simple LD Estimators
#'
#' Provides quick and simple estimates of LD and their corresponding
#' standard errors. These estimates are just based on the sample
#' moments of the genotypes or the posterior mean genotypes.
#'
#' For large sequencing depth and large n, these estimates perform as well
#' as the MLE's but are much faster to calculate. For small n and, in
#' particular, small sequencing depth, these estimates tend to be very
#' biased.
#'
#' @param ga Either the genotype at locus 1 or the posterior mean genotype at
#'     locus 1.
#' @param gb Either the genotype at locus 2 or the posterior mean genotype at
#'     locus 2.
#' @param K The ploidy of hte species. Assumed to be the same for all
#'     individuals.
#'
#' @inherit ldest return
#'
#' @author David Gerard
#'
#' @export
#'
ldsimp <- function(ga, gb, K) {
  n <- length(ga)
  pA <- mean(ga) / K
  pB <- mean(gb) / K
  D <- stats::cov(ga, gb) / K
  sda <- stats::sd(ga)
  sdb <- stats::sd(gb)
  D_se <- sqrt(max(sda * sdb / K ^ 2 - D ^ 2, 0)) / sqrt(n)
  r <- stats::cor(ga, gb)
  r_se <- (1 - r ^ 2) / sqrt(n)
  z <- atanh(r)
  z_se <- 1 / sqrt(n - 3)
  r2 <- r ^ 2
  r2_se <- 2 * abs(r) * (1 - r ^ 2) / sqrt(n)
  g <- log(-log(r2))
  g_se <- 2 * (1 - r ^ 2) / abs(r * log(r ^ 2) * sqrt(n))

  phat <- rep(NA_real_, 4)
  phat[[4]] <- D + pA * pB
  phat[[3]] <- pB - phat[[4]]
  phat[[2]] <- pA - phat[[4]]
  phat[[1]] <- 1 - sum(phat[2:4])
  phat[phat < 0] <- 0
  phat <- phat / sum(phat)

  retvec <- c(D         = D,
              D_se      = D_se,
              r2        = r2,
              r2_se     = r2_se,
              r         = r,
              r_se      = r_se,
              z         = z,
              z_se      = z_se,
              g         = g,
              g_se      = g_se,
              p_ab      = phat[[1]],
              p_Ab      = phat[[2]],
              p_aB      = phat[[3]],
              p_AB      = phat[[4]])

  return(retvec)
}
