###########################
## Expectations of truncated normal components.
###########################

#' @param s The standard deviation
#' @param K the upper bound
#' @author David Gerard
#' @noRd
afun <- function(s, K) {
  alpha <- -K / (2 * s)
  2 * alpha * stats::dnorm(alpha) / (1 - 2 * stats::pnorm(alpha))
}

#' @inheritParams afun
#' @author David Gerard
#' @return Approximate variance, not SD
#' @noRd
ev <- function(s, K) {
  s^2 * (1 + afun(s, K))
}

#' @inheritParams afun
#' @param muy The observed mean of posterior variances
#' @author David Gerard
#' @noRd
ev_obj <- function(s, muy, K) {
  (ev(s = s, K = K) - muy)^2
}

#' @param varx The variance of posterior means (a vector, one value is one locus).
#' @param muy The mean of posterior variances (a vector, on value is one locus).
#' @author David Gerard
#' @noRd
mom_adjust <- function(varx, muy, K) {
  stopifnot(length(muy) == length(varx))
  n <- length(muy)
  muy_new <- rep(NA_real_, length.out = n)
  varx_new <- rep(NA_real_, length.out = n)
  for (i in seq_along(varx)) {
    oout <- stats::optim(par = sqrt(muy[[i]]),
                         fn = ev_obj,
                         method = "Brent",
                         lower = 0,
                         upper = K/2,
                         muy = muy[[i]],
                         K = K)
    muy_new[[i]] <- oout$par^2
    varx_new[[i]] <- varx[[i]] / (1 + afun(s = sqrt(muy_new[[i]]), K = K))^2
  }

  return(list(muy = muy, varx = varx))
}
