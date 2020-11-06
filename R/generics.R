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

#' Plot the output of \code{\link{mldest}()} or
#' \code{\link{sldest}()} using \code{\link[corrplot]{corrplot}()}
#'
#' Formats the LD estimates in the form of a matrix and creates a heatmap of
#' these estimates. This heatmap is created using the
#' \code{\link[corrplot]{corrplot}} R package. I've adjusted a lot of the defaults
#' to suit my visualization preferences.
#'
#' For greater plotting flexibility, see \code{\link[corrplot]{corrplot}()}
#' for the parameter options.
#'
#' @param x An object of class \code{lddf}, usually created using
#'     either \code{\link{mldest}()} or \code{\link{sldest}()}.
#' @param element Which element of \code{x} should we plot?
#' @param type Character, \code{"full"},
#'     \code{"upper"} (default) or \code{"lower"}, display
#'     full matrix, lower triangular or upper
#'     triangular matrix.
#' @param method See \code{\link[corrplot]{corrplot}()} for available options.
#'     Default value is \code{"color"}.
#' @param is.corr See \code{\link[corrplot]{corrplot}()}. Default behavior
#'     is \code{TRUE} if an element is constrained
#'     between -1 and 1 and \code{FALSE} otherwise.
#' @param tl.pos See \code{\link[corrplot]{corrplot}()}. Default value
#'     is \code{"n"}.
#' @param na.label See \code{\link[corrplot]{corrplot}()}. Default value
#'     is \code{"square"}.
#' @param diag Logical, whether display the correlation coefficients
#'     on the principal diagonal.
#' @param title What should the title be? Defaults to the element name.
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
#' lddf <- mldest(geno = genomat,
#'                K = K,
#'                nc = nc,
#'                type = "hap")
#'
#' ## Plot estimates of z
#' plot(lddf, element = "z")
#'
#' @author David Gerard
#'
#' @export
plot.lddf <- function(x,
                      element = "r2",
                      type = c("upper",
                               "full",
                               "lower"),
                      method = c("color",
                                 "circle",
                                 "square",
                                 "ellipse",
                                 "number",
                                 "shade",
                                 "pie"),
                      diag = FALSE,
                      is.corr = NULL,
                      tl.pos = "n",
                      title = NULL,
                      na.label = "square",
                      ...) {
  type <- match.arg(type)
  stopifnot(is.logical(diag))
  stopifnot(length(element) == 1)
  stopifnot(element %in% names(x))
  method <- match.arg(method)
  if (is.null(is.corr)) {
    is.corr <- all((x[[element]] >= -1) & (x[[element]] <= 1), na.rm = TRUE)
  }
  if (is.null(title)) {
    title <- element
  } else {
    stopifnot(is.character(title))
    stopifnot(length(title) == 1)
  }

  cormat <- format_lddf(obj = x, element = element)

  if (diag) {
    diag(cormat) <- 1
  }
  if (type != "upper") {
    cormat[lower.tri(cormat)] <- t(cormat)[lower.tri(cormat)]
  }

  corrplot::corrplot(corr = cormat,
                     type = type,
                     diag = diag,
                     is.corr = is.corr,
                     method = method,
                     tl.pos = tl.pos,
                     title = title,
                     na.label = na.label,
                     mar = c(0,0,1,0),
                     ...)
}


#' Format an element of \code{\link{mldest}()} or
#' \code{\link{sldest}()} into an
#' upper-triangular matrix.
#'
#' Formats the LD estimates and standard errors output
#' from running \code{\link{mldest}()} or \code{\link{sldest}()}
#' into a more conventional upper-triangular matrix.
#'
#' @param obj An object of class \code{lddf}, usually output from
#'     running either \code{\link{mldest}()} or
#'     \code{\link{sldest}()}.
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
#' lddf <- mldest(geno = genomat,
#'                K = K,
#'                nc = nc,
#'                type = "hap")
#'
#' ## Obtain the D estimates in matrix form
#' Dmat <- format_lddf(obj = lddf, element = "D")
#' Dmat
#'
#' @export
format_lddf <- function(obj,
                        element = "r2") {
  stopifnot(is.lddf(obj))
  stopifnot(length(element) == 1)
  stopifnot(element %in% names(obj))
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
