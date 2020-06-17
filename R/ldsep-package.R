#' Linkage Disequilibrium Shrinkage Estimation for Polyploids
#'
#' Estimate haplotypic or composite pairwise linkage disequilibrium
#' (LD), using either genotypes or genotype likelihoods. Support is
#' provided to estimate the popular measures of LD: the LD coefficient D,
#' the standardized LD coefficient D', and the Pearson correlation
#' coefficient r. All estimates are returned with corresponding
#' standard errors. These estimates and standard errors can then be used
#' for shrinkage estimation. The main functions are \code{\link{ldest}()},
#' \code{\link{mldest}()}, \code{\link{plot.lddf}()},
#' \code{\link{format_lddf}()}, and \code{\link{ldshrink}()}.
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
