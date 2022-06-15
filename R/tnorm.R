###########################
## Expectations of truncated normal components.
###########################

#' @importFrom fastGHQuad aghQuad
#' @importFrom stats dnorm
#' @importFrom stats pnorm
NULL

myenv <- new.env(parent = emptyenv())
myenv$ghrule <- fastGHQuad::gaussHermiteData(n = 100)

i1 <- function(z, s, K, mu, sig) {
  ((dnorm(-z/s) - dnorm((K-z)/s)) / (pnorm((K-s)/s) - pnorm(-z/s))) * dnorm(z, mu, sqrt(s^2 + sig^2))
}

#' @param mu mean of genotypes
#' @param sig Standard deviation of genotypes
#' @param s Standard deviation of data
#' @param K Ploidy
#' @param rule Output of fastGHQuad::gaussHermiteData()
#'
#' @examples
#' mueq(1, 1, 0, 4)
#'
#' @noRd
mueq <- function(mu, sig, s, K) {
  mu + aghQuad(g = i1, muHat = mu, sigmaHat = sig, rule = myenv$ghrule, s = s, K = K, mu = mu, sig = sig)
}



