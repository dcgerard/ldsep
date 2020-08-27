#' Linkage Disequilibrium Shrinkage Estimation for Polyploids
#'
#' Estimate gametic or composite pairwise linkage disequilibrium
#' (LD) in polyploids, using either genotypes or genotype likelihoods. Support is
#' provided to estimate the popular measures of LD: the LD coefficient D,
#' the standardized LD coefficient D', and the Pearson correlation
#' coefficient r. All estimates are returned with corresponding
#' standard errors. These estimates and standard errors can then be used
#' for shrinkage estimation.
#'
#' @section Functions:
#'
#' The main functions are:
#' \describe{
#'   \item{\code{\link{ldest}()}}{Estimates pairwise LD.}
#'   \item{\code{\link{mldest}()}}{Iteratively apply \code{\link{ldest}()}
#'       across many pairs of SNPs.}
#'   \item{\code{\link{sldest}()}}{Iteratively apply \code{\link{ldest}()}
#'       along a sliding window of fixed length.}
#'   \item{\code{\link{plot.lddf}()}}{Plot method for the output of
#'       \code{\link{mldest}()} and \code{\link{sldest}()}.}
#'   \item{\code{\link{format_lddf}()}}{Format the output of
#'       \code{\link{mldest}()} and \code{\link{sldest}()} into a matrix.}
#'   \item{\code{\link{ldshrink}()}}{Shrink correlation estimates
#'       using adaptive shrinkage (Stephens, 2017; Dey and Stephens, 2018).}
#' }
#'
#' @section Citation:
#' If you find the methods in this package useful, please run the following
#' in R for citation information: \code{citation("ldsep")}
#'
#'
#' @importFrom stats var
#' @importFrom foreach %dopar%
#' @useDynLib ldsep, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @docType package
#' @name ldsep-package
#' @aliases ldsep
#'
#' @author David Gerard
NULL
