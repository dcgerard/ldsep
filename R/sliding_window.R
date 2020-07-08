#############################
## sliding window ld estimation
############################



#' Sliding window LD estimation
#'
#' This function is a wrapper for \code{\link{ldest}()} for estimating LD
#' along a sliding window of a fixed size. Support is provided for parallelization through the
#' foreach and doParallel packages.
#'
#' See \code{\link{ldest}()} for details on the different types of LD
#' estimators supported.
#'
#' @inheritParams mldest
#' @param win The window size. Pairwise LD will be estimated plus or minus
#'     these many positions. Larger sizes significantly increase the
#'     computational load.
#'
#' @inherit mldest return
#'
#' @examples
#' set.seed(1)
#'
#' ## Simulate genotypes when true correlation is 0
#' nloci <- 100
#' nind  <- 100
#' win <- 5
#' K <- 6
#' nc <- 1
#' genomat <- matrix(sample(0:K, nind * nloci, TRUE), nrow = nloci)
#'
#' ## Composite LD estimates
#' lddf <- sldest(geno = genomat,
#'                K = K,
#'                win = win,
#'                nc = nc,
#'                type = "comp")
#' plot(lddf, element = "z")
#'
#' @author David Gerard
#'
#' @export
sldest <- function(geno,
                   K,
                   win = 50,
                   nc = 1,
                   type = c("hap", "comp"),
                   model = c("norm", "flex"),
                   pen = ifelse(type == "hap", 2, 1),
                   se = TRUE) {
  model <- match.arg(model)
  type <- match.arg(type)
  if (length(dim(geno)) == 2) {
    outdf <- mldest_geno(genomat = geno,
                         K = K,
                         nc = nc,
                         pen = pen,
                         type = type,
                         se = se,
                         win = win)
  } else if (length(dim(geno)) == 3) {
    outdf <- mldest_genolike(genoarray = geno,
                             nc = nc,
                             pen = pen,
                             type = type,
                             model = model,
                             se = se,
                             win = win)
  } else {
    stop("sldest: geno needs to either be a matrix or a three-way array.")
  }
  return(outdf)
}
