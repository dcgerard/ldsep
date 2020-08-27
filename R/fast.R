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
#'
#' @param gp A three-way array with dimensions SNPs by individuals by dosage.
#'     That is, \code{gp[i, j, k]} is the posterior probability of
#'     dosage \code{k-1} for individual \code{j} at SNP \code{i}.
#' @param type What LD measure should we estimate?
#'     \describe{
#'       \item{\code{r}}{The Pearson correlation.}
#'       \item{\code{r2}}{The squared Pearson correlation.}
#'       \item{\code{z}}{The Fisher-z transformed Pearson correlation.}
#'       \item{\code{D}}{The LD coefficient.}
#'       \item{\code{Dprime}}{The standardized LD coefficient.}
#'     }
#'     Note that these are all \emph{composite} measures of LD.
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
#'   \item{\code{ldmat}}{The bias-corrected LD matrix.}
#'   \item{\code{semat}}{A matrix 0f standard errors of the corresponding
#'       estimators of LD.}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' ## Load the data -----
#' data("gp") # posterior probs
#' ldout <- ldfast(gp, "r")
#' ldout$ldmat
#' ldout$semat
#'
#' ldout <- ldfast(gp, "D")
#' ldout$ldmat
#' ldout$semat
#'
#' ldout <- ldfast(gp, "Dprime")
#' ldout$ldmat
#' ldout$semat
#'
#' @export
ldfast <- function(gp, type = c("r", "r2", "z", "D", "Dprime")) {
  stopifnot(inherits(gp, "array"))
  stopifnot(length(dim(gp)) == 3)
  type <- match.arg(type)

  nsnp <- dim(gp)[[1]]
  nind <- dim(gp)[[2]]
  ploidy <- dim(gp)[[3]] - 1

  cormat <- matrix(NA_real_, ncol = nsnp, nrow = nsnp)
  semat <- matrix(NA_real_, ncol = nsnp, nrow = nsnp)

  if (type == "D") {
    ldfast_calc(cormat = cormat, semat = semat, gp = gp, type = "a")
  } else if (type == "Dprime") {
    ldfast_calc(cormat = cormat, semat = semat, gp = gp, type = "c")
  } else {
    ldfast_calc(cormat = cormat, semat = semat, gp = gp, type = "b")
  }

  if (type == "r2") {
    semat <- semat * abs(cormat) * 2
    cormat <- cormat ^ 2
  } else if (type == "z") {
    semat <- semat / (1 - cormat ^ 2)
    cormat <- atanh(cormat)
  }

  return(list(ldmat = cormat, semat = semat))
}


#' same as ldfast, but just calculate the ld, not the se
#'
#'
#' @author David Gerard
#'
#' @noRd
ldfast_justmean <- function(gp, type = c("r", "r2", "z", "D", "Dprime")) {
  stopifnot(inherits(gp, "array"))
  stopifnot(length(dim(gp)) == 3)
  type <- match.arg(type)

  nsnp <- dim(gp)[[1]]
  nind <- dim(gp)[[2]]
  ploidy <- dim(gp)[[3]] - 1

  pm_mat <- matrix(NA_real_, nrow = nsnp, ncol = nind)
  pv_mat <- matrix(NA_real_, nrow = nsnp, ncol = nind)
  fill_pm(pm = pm_mat, gp = gp)
  fill_pv(pv = pv_mat, pm = pm_mat, gp = gp)

  varx <- apply(X = pm_mat, MARGIN = 1, FUN = `var`, na.rm = TRUE)
  muy <- rowMeans(x = pv_mat, na.rm = TRUE)

  ## calculate correlation
  rr <- sqrt((muy + varx) / varx)
  ldmat <- rr * stats::cor(t(pm_mat), use = "pairwise.complete.obs") * rep(rr, each = nsnp)
  if (type != "Dprime") {
    ldmat[ldmat > 1] <- 1
    ldmat[ldmat < -1] <- -1
  }

  if (type == "D") {
    sdvec <- sqrt(muy + varx)
    ldmat <- (sdvec * ldmat * rep(sdvec, each = nsnp)) / ploidy
  } else if (type == "Dprime") {
    mux <- rowMeans(pm_mat, na.rm = TRUE)
    sdvec <- sqrt(muy + varx)
    ldmat <- sdvec * ldmat * rep(sdvec, each = nsnp) / ploidy

    deltam_neg <- pmin(outer(X = mux, Y = mux, FUN = `*`),
                       outer(X = ploidy - mux, Y = ploidy - mux, FUN = `*`)) / ploidy ^ 2
    deltam_pos <- pmin(outer(X = mux, Y = ploidy - mux, FUN = `*`),
                       outer(X = ploidy - mux, Y = mux, FUN = `*`)) / ploidy ^ 2

    ldmat[ldmat < 0] <- ldmat[ldmat < 0] / deltam_neg[ldmat < 0]
    ldmat[ldmat > 0] <- ldmat[ldmat > 0] / deltam_pos[ldmat > 0]

    ldmat[ldmat > ploidy] <- ploidy
    ldmat[ldmat < -ploidy] <- -ploidy
  } else if (type == "r2") {
    ldmat <- ldmat ^ 2
  } else if (type == "z") {
    ldmat <- atanh(ldmat)
  } else {
    ## do nothing
  }

  return(ldmat)
}

ldfast_r <- function(gp, type = c("r", "r2", "D"), pv = NULL) {
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
