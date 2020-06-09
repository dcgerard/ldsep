#' Linkage Disequilibrium Shrinkage Estimation for Polyploids
#'
#' Estimates pairwise linkage disequilibrium (LD), using either
#' genotype estimates or genotype likelihoods. Functions are then provided
#' for shrinkage estimation using multidimensional scaling and adaptive
#' shrinkage. The idea is to expand this to allow for a sliding window
#' estimation procedure for genome-wise LD estimation.
#'
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
