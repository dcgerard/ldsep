#' Tests if an argument is a \code{lddf} object.
#'
#' @param x Anything.
#'
#' @return A logical. \code{TRUE} if \code{x} is a \code{lddf} object,
#'     and \code{FALSE} otherwise.
#'
#' @author David Gerard
#'
#' @export
#'
#' @examples
#' is.lddf("anything")
#' # FALSE
#'
is.lddf <- function(x) {
  inherits(x, "lddf")
}

#' Plot the output of \code{\link{mldest_geno}()} or
#' \code{\link{mldest_genolike}()} using \code{\link[corrplot]{corrplot}()}
#'
#' Uses the \code{\link[corrplot]{corrplot}} R package to visualize
#' LD estimates.
#'
#' @param x An object of class \code{lddf}, usually created using
#'     either \code{\link{mldest_geno}()} or \code{\link{mldest_genolike}()}.
#' @param element Which element of \code{x} should be plot?
#' @param type Character, \code{"full"},
#'     \code{"upper"} (default) or \code{"lower"}, display
#'     full matrix, lower triangular or upper
#'     triangular matrix.
#' @param diag Logical, whether display the correlation coefficients
#'     on the principal diagonal.
#' @param ... Additional arguments to pass to
#'     \code{\link[corrplot]{corrplot}()}. See the documentation of that
#'     function for options.
#'
#' @return (Invisibly) returns a matrix of the selected elements.
#'
#' @examples
#' set.seed(1)
#'
#' ## Simulate genotypes when true correlation is 0
#' nloci <- 5
#' nind  <- 100
#' K <- 6
#' nc <- 1
#' genomat <- matrix(sample(0:K, nind * nloci, TRUE), nrow = nloci)
#'
#' ## Haplotypic LD estimates
#' lddf <- mldest_geno(genomat = genomat,
#'                     K = K,
#'                     nc = nc,
#'                     type = "hap")
#'
#' ## Plot estimates of z
#' plot(lddf, element = "z")
#'
#' @author David Gerard
#'
#' @export
plot.lddf <- function(x,
                      element = c("z",
                                  "z_se",
                                  "D",
                                  "D_se",
                                  "Dprime",
                                  "Dprime_se",
                                  "r2",
                                  "r2_se",
                                  "r",
                                  "r_se",
                                  "p_ab",
                                  "p_Ab",
                                  "p_aB",
                                  "p_AB"),
                      type = c("upper", "full", "lower"),
                      diag = FALSE,
                      ...) {
  type <- match.arg(type)
  stopifnot(is.logical(diag))
  element <- match.arg(element)
  cormat <- format_lddf(obj = x, element = element)

  if (diag) {
    diag(cormat) <- 1
  }
  if (type != "upper") {
    cormat[lower.tri(cormat)] <- t(cormat)[lower.tri(cormat)]
  }
  if (element %in% c("z", "z_se")) {
    is.corr <- FALSE
  } else {
    is.corr <- TRUE
  }
  corrplot::corrplot(corr = cormat,
                     type = type,
                     diag = diag,
                     is.corr = is.corr,
                     ...)
}


#' Format an element of \code{\link{mldest_geno}()} or
#' \code{\link{mldest_genolike}()} into an
#' upper-triangular matrix.
#'
#' Formats the correlation estimates and standard errors output
#' from running \code{\link{mldest_geno}()} or \code{\link{mldest_genolike}()}
#' into a more conventional upper-triangular matrix.
#'
#' @param obj An object of class \code{lddf}, usually output from
#'     running either \code{\link{mldest_geno}()} or
#'     \code{\link{mldest_genolike}()}.
#' @param element Which element in \code{obj} should we format into an
#'     upper-triangular matrix?
#'
#' @return A matrix of the selected elements. Only the upper-triangle of the
#'     matrix is filled. The lower-triangle and the diagonal are \code{NA}'s.
#'
#' @author David Gerard
#'
#' @examples
#' set.seed(1)
#'
#' ## Simulate genotypes when true correlation is 0
#' nloci <- 5
#' nind  <- 100
#' K <- 6
#' nc <- 1
#' genomat <- matrix(sample(0:K, nind * nloci, TRUE), nrow = nloci)
#'
#' ## Haplotypic LD estimates
#' lddf <- mldest_geno(genomat = genomat,
#'                     K = K,
#'                     nc = nc,
#'                     type = "hap")
#'
#' ## Obtain the D estimates in matrix form
#' Dmat <- format_lddf(obj = lddf, element = "D")
#' Dmat
#'
#' @export
format_lddf <- function(obj,
                        element = c("z",
                                    "z_se",
                                    "D",
                                    "D_se",
                                    "Dprime",
                                    "Dprime_se",
                                    "r2",
                                    "r2_se",
                                    "r",
                                    "r_se",
                                    "p_ab",
                                    "p_Ab",
                                    "p_aB",
                                    "p_AB")) {
  stopifnot(is.lddf(obj))
  element <- match.arg(element)
  nloci <- max(max(obj$i), max(obj$j))
  cormat <- matrix(NA_real_, ncol = nloci, nrow = nloci)
  cormat[as.matrix(obj[, c("i", "j")])] <- obj[[element]]

  if (all(is.na(obj$snpi)) & all(is.na(obj$snpj))) {
    rownames(cormat) <- seq_len(nloci)
    colnames(cormat) <- seq_len(nloci)
  } else {
    agdf1 <- stats::aggregate(formula = snpi ~ i, data = obj, FUN = `unique`)
    agdf2 <- stats::aggregate(formula = snpj ~ j, data = obj, FUN = `unique`)
    names(agdf2) <- c("i", "snpi")
    agdf <- merge(x = agdf1,
                  y = agdf2,
                  by = c("i", "snpi"),
                  all = TRUE)
    rownames(cormat) <- rep(NA_character_, length = nloci)
    colnames(cormat) <- rep(NA_character_, length = nloci)
    rownames(cormat)[agdf$i] <- agdf$snpi
    colnames(cormat)[agdf$i] <- agdf$snpi
  }

  return(cormat)
}
