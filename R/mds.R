###############################
## Functions for shrinkage of correlation matrix
###############################

#' Shrinks Fisher-z transformed correlation estimates and returns resulting
#' correlation estimates.
#'
#' This function is a wrapper for adaptive shrinkage (Stephens, 2017) on the
#' Fisher-z transformed estimates of the Pearson correlation. This approach
#' was proposed in Dey and Stephens (2018) but is re-implemented here for now
#' since the CorShrink package is not available on CRAN.
#'
#' @param zmat The matrix of Fisher-z transformed correlation estimates.
#' @param smat The matrix of standard errors of the Fisher-z transformed
#'     correlation estimates.
#'
#' @return A matrix of correlation estimates. These are posterior means
#'     of the correlation estimates after applying the CorShrink method
#'     (Dey and Stephens, 2018).
#'
#' @references
#' \itemize{
#' \item{Stephens, Matthew. "False discovery rates: a new deal."
#'       Biostatistics 18, no. 2 (2017): 275-294.}
#' \item{Dey, Kushal K., and Matthew Stephens. "CorShrink:
#'       Empirical Bayes shrinkage estimation of correlations,
#'       with applications." bioRxiv (2018): 368316.}
#' }
#'
#' @author David Gerard
#'
#' @export
zshrink <- function(zmat, smat) {
  stopifnot(dim(zmat) == dim(smat))
  betahat <- zmat[upper.tri(zmat)]
  sebetahat <- smat[upper.tri(smat)]
  aout <- ashr::ash(betahat     = betahat,
                    sebetahat   = sebetahat,
                    mixcompdist = "uniform")
  corvec <- tanh(ashr::get_pm(aout))
  cormat <- matrix(NA_real_, ncol = ncol(zmat), nrow = nrow(zmat))
  diag(cormat) <- 1
  cormat[upper.tri(cormat)] <- corvec
  cormat[lower.tri(cormat)] <- corvec
  return(cormat)
}

#' Obtain shrinkage estimates of correlation from output of
#' \code{\link{mldest}()} or \code{\link{sldest}()}.
#'
#' This will take the output of either \code{\link{mldest}()} or
#' \code{\link{sldest}()}, shrink the Fisher-z transformed
#' correlation estimates using \code{\link[ashr]{ash}()}, then return
#' the corresponding correlation estimates. You can obtain estimates of
#' r^2 by just squaring these estimates.
#'
#' @param obj An object of class \code{lddf}, usually created using
#'     either \code{\link{mldest}()} or \code{\link{sldest}()}.
#'
#' @return A correlation matrix.
#'
#' @author David Gerard
#'
#' @export
ldshrink <- function(obj) {
  stopifnot(is.lddf(obj))
  zmat <- format_lddf(obj = obj, element = "z")
  smat <- format_lddf(obj = obj, element = "z_se")
  cormat <- zshrink(zmat = zmat, smat = smat)
  return(cormat)
}
