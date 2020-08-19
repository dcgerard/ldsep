###################
## FAST LD Correction
###################

#' Fast bias-correction for LD Estimation
#'
#' Estimates the reliability ratios from posterior marginal moments and uses
#' these to correct the biases in linkage disequilibrium estimation
#' caused by genotype uncertainty.
#'
#' @section Details:
#'
#' Sometimes, the reliability ratio results in a correlation greater than 1
#' or less than -1 (due to sample variability). We truncate these occasions
#' at 1 and -1, respectively.
#'
#' This won't work for missing data in \code{gp} yet, but the plan is to
#' allow missing data in the future.
#'
#' @section Mathematical formulation:
#' Let
#' \itemize{
#'   \item{\eqn{r} be the sample correlation of posterior mean genotypes
#'       between loci 1 and 2,}
#'   \item{\eqn{a1} be the sample variance of posterior means at locus 1,}
#'   \item{\eqn{a2} be the sample variance of posterior means at locus 2,}
#'   \item{\eqn{b1} be the sample mean of posterior variances at locus 1, and}
#'   \item{\eqn{b2} be the sample mean of posterior variances at locus 2.}
#' }
#' Then the estimated Pearson correlation between the genotypes at
#' loci 1 and 2 is
#' \deqn{[(a1 + b1)/a1][(a2 + b2)/a2]r.}
#' This is the estimated LD when \code{pv = NULL}.
#'
#' Alternatively, suppose that we know the prior variances and we let
#' \itemize{
#'   \item{\eqn{s1} be the prior variance at locus 1,}
#'   \item{\eqn{s2} be the prior variance at locus 2, and}
#'   \item{\eqn{t} be the sample covariance of posterior mean genotypes
#'       between loci 1 and 2.}
#' }
#' Then the estimated covariance between genotypes at loci 1 and 2 is
#' \deqn{[s1/(s1 - b1)][s2/(s2-b2)]t.}
#' We obtain the estimate Pearson correlation by normalizing by the
#' prior standard deviations
#' \deqn{[\sqrt{s1}/(s1 - b1)][\sqrt{s2}/(s2-b2)]t.}
#' This is the estimated LD when \code{pv} is not \code{NULL}.
#'
#' @param gp A three-way array with dimensions SNPs by individuals by dosage.
#'     That is, \code{gp[i, j, k]} is the posterior probability of
#'     dosage \code{k-1} for individual \code{j} at SNP \code{i}.
#' @param pv A vector of numerics. Optionally provided prior variances.
#'     This should only be provided if these prior variances were estimated
#'     adaptively.
#' @param type What LD measure should we estimate? The Pearson correlation
#'     (\code{type = "r"}), the squared Pearson correlation
#'     (\code{type = "r2"}), or the LD coefficient (\code{type = "D"})?
#'     These are all \emph{composite} measures of LD.
#'
#' @seealso
#' \describe{
#' \item{\code{\link{gl_to_gp}()}}{Normalize genotype likelihoods to
#'     posterior probabilities using naive uniform prior.}
#' \item{\code{\link{pvcalc}()}}{Calculate prior variances given genotype
#'     priors.}
#' }
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{cor}}{The bias-corrected LD matrix.}
#'   \item{\code{rr}}{The estimated reliability ratio for the
#'       \emph{correlation}. This is how much you multiply the
#'       Pearson correlation by to obtain the bias-adjusted Pearson
#'       correlation estimate. Square this term to get the reliability
#'       ratio for the LD coefficient.}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' ## Load the data -----
#' data("uit") # mupdog() fits
#' pv <- pvcalc(uit$snpdf[, paste0("Pr_", 0:4)]) # prior probs
#' data("gp") # posterior probs
#'
#' ## Estimate squared correlation using just posterior moments -----
#' cout1 <- ldfast(gp = gp, type = "r2")
#' corrplot::corrplot(corr = cout1$cor,
#'                    diag = FALSE,
#'                    type = "upper",
#'                    method = "color")
#' hist(cout1$rr, xlab = "Reliability Ratio")
#'
#' ## Estimate squared correlation also using prior variances -----
#' cout2 <- ldfast(gp = gp, type = "r2", pv = pv)
#' corrplot::corrplot(corr = cout2$cor,
#'                    diag = FALSE,
#'                    type = "upper",
#'                    method = "color")
#' hist(cout2$rr, xlab = "Reliability Ratio")
#'
#' ## The adjustments are really close in these data -----
#' plot(x = cout1$rr,
#'      y = cout2$rr,
#'      xlab = "Without Prior Variance",
#'      ylab = "With Prior Variance")
#' abline(0, 1)
#'
#' @export
ldfast <- function(gp, type = c("r", "r2", "D"), pv = NULL) {
  stopifnot(inherits(gp, "array"))
  stopifnot(length(dim(gp)) == 3)
  type <- match.arg(type)
  if (!is.null(pv)) {
    stopifnot(inherits(pv, "numeric"))
    stopifnot(length(pv) == dim(gp)[[1]])
  }
  ploidy <- dim(gp)[[3]] - 1

  ds <- matrix(data = 0, nrow = dim(gp)[[2]], ncol = dim(gp)[[1]])
  ds_from_gp(ds = ds, gp = gp)
  if (is.null(pv)) {
    rr <- post_rr(gp = gp, ds = ds)
    cmat <- rr[, 1] * coop::pcor(x = ds, use = "pairwise.complete.obs") * rep(rr[, 1], each = ncol(ds))
    # cmat2 <- sweep(x = rr[, 1] * coop::pcor(x = ds, use = "pairwise.complete.obs"), MARGIN = 2, STATS = rr[, 1], FUN = `*`)
  } else {
    rr <- prior_rr(gp = gp, ds = ds, priorvar = pv)
    cmat <- (rr[, 1]^2) * coop::covar(x = ds, use = "pairwise.complete.obs") * rep(rr[, 1]^2, each = ncol(ds))
    # cmat2 <- sweep(x = (rr[, 1]^2) * coop::covar(x = ds, use = "pairwise.complete.obs"), MARGIN = 2, STATS = rr[, 1]^2, FUN = `*`)
    diag(cmat) <- rr[, 2] ^ 2
    cmat <- stats::cov2cor(cmat)
  }
  cmat[cmat > 1] <- 1
  cmat[cmat < -1] <- -1

  if (type == "r2") {
    cmat <- cmat ^ 2
  } else if (type == "D") {
    cmat <- sweep(x = rr[, 2] * cmat, MARGIN = 2, STATS = rr[, 2], FUN = `*`) / ploidy
  }

  return(list(cor = cmat, rr = rr[, 1]))
}

# alternate c++ based code
ldfast_c <- function(gp, type = c("r", "r2", "D"), pv = NULL) {
  stopifnot(length(dim(gp)) == 3)
  type <- match.arg(type)
  if (!is.null(pv)) {
    stopifnot(inherits(pv, "numeric"))
    stopifnot(length(pv) == dim(gp)[[1]])
  }

  if (type %in% c("r", "r2")) {
    cout <- ldfast_post(gp = gp, type = "r", priorvar_ = pv)
    if (type == "r2") {
      cout$cor <- cout$cor ^ 2
    }
  } else {
    cout <- ldfast_post(gp = gp, type = "D", priorvar_ = pv)
  }
  rownames(cout$cor) <- dimnames(gp)[[1]]
  colnames(cout$cor) <- dimnames(gp)[[1]]
  return(cout)
}

#' Normalize genotype likelihoods to posterior probabilities.
#'
#' This will take genotype log-likelihoods and normalize them to
#' sum to one. This corresponds to using a naive discrete uniform prior
#' over the genotypes, which is typically OK if we are not adaptively
#' estimating likelihood elements using this prior.
#'
#' @param gl A three dimensional array of genotype \emph{log}-likelihoods.
#'     Element \code{gl[i, j, k]} is the genotype log-likelihood of dosage
#'     \code{k} for individual \code{j} at SNP \code{i}.
#'
#' @return A three-dimensional array, of the same dimensions as \code{gl},
#'     containing the posterior probabilities of each dosage.
#'
#' @author David Gerard
#'
#' @examples
#' data("glike")
#' class(glike)
#' dim(glike)
#' gl_to_gp(glike)
#'
#' @export
gl_to_gp <- function(gl) {
  stopifnot(inherits(gl, "array"))
  stopifnot(length(dim(gl)) == 3)
  gp <- apply(X = gl, MARGIN = c(1, 2), FUN = function(x) exp(x - log_sum_exp(x)))
  gp <- aperm(gp, c(2, 3, 1))
  return(gp)
}

#' Calculate prior variances from a matrix of prior genotype probabilities.
#'
#' Given a matrix of prior probabilities for the genotypes at each SNP,
#' this function will calculate the prior variance of genotypes.
#'
#' @param priormat A matrix of prior genotype probabilities. Element
#'     \code{priormat[i, j]} is the prior probability of dosage \code{j}
#'     at SNP \code{i}.
#'
#' @author David Gerard
#'
#' @return A vector of prior variances.
#'
#' @examples
#' data("uit")
#' priormat <- uit$snpdf[, paste0("Pr_", 0:4)]
#' pvcalc(priormat)
#'
#' @export
pvcalc <- function(priormat) {
  stopifnot(length(dim(priormat)) == 2)
  ploidy <- ncol(priormat) - 1
  rowSums(sweep(x = priormat, MARGIN = 2, STATS = (0:ploidy)^2, FUN = `*`)) -
    rowSums(sweep(x = priormat, MARGIN = 2, STATS = 0:ploidy, FUN = `*`))^2
}
