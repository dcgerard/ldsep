####################
## Truncated normal stuff
####################

#' First two moments of the truncated normal
#'
#' This is not a stable implementation.
#'
#' @param mu The mean.
#' @param sigma The standard deviation.
#' @param a The lower bound.
#' @param b The upper bound.
#'
#' @return A vector of length 2. The mean and standard deviation.
#'
#' @author David Gerard
#'
#' @examples
#' mom_tnorm(1, 2, 0, 4)
#'
#' @noRd
mom_tnorm <- function(mu, sigma, a, b) {
  alpha <- (a - mu) / sigma
  beta <- (b - mu) / sigma
  Z <- stats::pnorm(beta) - stats::pnorm(alpha)

  par <- c(
    mu + (stats::dnorm(alpha) - stats::dnorm(beta)) / Z * sigma,
    sigma * sqrt(1 + (alpha * stats::dnorm(alpha) - beta * stats::dnorm(beta)) / Z -
                   ((stats::dnorm(alpha) - stats::dnorm(beta)) / Z)^2)
  )

  return(par)
}

#' Objective function used in \code{tnorm_solve()}.
#'
#' @param par Vector of length two. First is the tentative mean, second
#'     is the tentative standard deviation.
#' @param mu_obs The mean of the truncated normal.
#' @param sigma_obs The standard deviation of the truncated normal.
#' @param a The lower bound.
#' @param b The upper bound.
#'
#' @author David Gerard
#'
#' @noRd
tnorm_obj <- function(par, mu_obs, sigma_obs, a, b) {
  sum((mom_tnorm(mu = par[[1]], sigma = par[[2]], a = a, b = b) - c(mu_obs, sigma_obs))^2)
}

#' Get Truncated normal parameters
#'
#' Returns truncated normal parameters given the mean and standard devaition
#' of the truncated normal. Assume upper bound is ploidy and lower bound is 0.
#'
#' @param mu_obs The mean of the truncated normal.
#' @param sigma_obs The standard deviation of the truncated normal.
#' @param ploidy The upper bound of the truncated normal.
#'
#' @author David Gerard
#'
#' @noRd
tnorm_solve <- function(mu_obs, sigma_obs, ploidy) {
  oout <- stats::optim(par = c(mu_obs, sigma_obs),
                       fn = tnorm_obj,
                       lower = c(0.01, 0.01),
                       upper = c(ploidy - 0.01, ploidy - 0.01),
                       method = "L-BFGS-B",
                       mu_obs = mu_obs,
                       sigma_obs = sigma_obs,
                       a = 0,
                       b = ploidy)
  return(c(mu = oout$par[[1]], sigma = oout$par[[2]]))
}

#' Find two numbers, above and below, closest to y.
#'
#' @param x A vector.
#' @param y A number.
#'
#' @author David Gerard
#'
#' @noRd
find_bounds <- function(y, x) {
  c(max(x[x <= y]), min(x[x >= y]))
}

#' Estimate parameters of truncated normal
#'
#' Scalable approximation to estimate the parameters of the truncated normal
#' by method-of-moments / maximum likelihood (equivalent methods in this case).
#' We do so by the look-up table method of Cohen (1957). The lower bound of
#' the truncated normal is assumed to be 0. Just pre and post shift the data
#' if not.
#'
#' @param xbar The sample mean.
#' @param s The sample standard deviation.
#' @param ploidy The upper bound of the truncated normal.
#'
#' @return A vector of length 2. The first element is the estimate of the
#'     mean parameter, the second element is the estimate of the standard
#'     deviation parameter.
#'
#' @references
#' \itemize{
#'   \item{Cohen, A. C. (1957). On the solution of estimating equations for truncated and censored samples from normal populations. \emph{Biometrika}, 44(1/2), 225-236. \doi{10.2307/2333256}}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' xbar <- 2
#' s <- 1
#' ploidy <- 6
#' tnorm_est(xbar = xbar, s = s, ploidy = ploidy)
#'
#' @export
tnorm_est <- function(xbar, s, ploidy) {
  vrat <- xbar / ploidy
  srat <- (s / ploidy)^2

  iclose <- which.min((tlook$vrat - vrat)^2 + (tlook$srat - srat)^2)

  xi1 <- tlook$xi1[[iclose]]
  xi2 <- tlook$xi2[[iclose]]

  s_new <- ploidy / (xi2 - xi1)
  xbar_new <- -s_new * xi1

  return(c(mu = xbar_new, sigma = s_new))
}
