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
#' Returns consistent and bias-corrected estimates of Linkage Disequilibrium.
#' The usual measures of LD are implemented: D, D', r, r2, and z
#' (Fisher-z of r). These are all \emph{composite} measures of LD, not
#' gametic measures of LD. They are always appropriate measures of association
#' between loci, but only correspond to gametic measures of LD when
#' Hardy-Weinberg equilibrium is fulfilled in polyploids.
#'
#' Due to sampling varibility, the estimates sometimes lie outside of the
#' theoretical boundary of the parameters being estimated. In such cases,
#' we truncate the estimates at the boundary and return \code{NA} for the
#' standard errors.
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
#' @param se Should we also return a matrix of standard errors (\code{TRUE})
#'     or not (\code{FALSE})? It is faster to not return standard errors.
#'     Defaults to \code{TRUE}.
#'
#' @seealso
#' \describe{
#' \item{\code{\link{gl_to_gp}()}}{Normalize genotype likelihoods to
#'     posterior probabilities using naive uniform prior.}
#' }
#'
#' @return A list with some or all of the following elements:
#' \describe{
#'   \item{\code{ldmat}}{The bias-corrected LD matrix.}
#'   \item{\code{rr}}{The estimated reliability ratio for each SNP. This
#'       is the multiplicative factor applied to the naive LD estimate
#'       for each SNP.}
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
#' ldout$rr
#' ldout$semat
#'
#' ldout <- ldfast(gp, "D")
#' ldout$ldmat
#' ldout$rr
#' ldout$semat
#'
#' ldout <- ldfast(gp, "Dprime")
#' ldout$ldmat
#' ldout$rr
#' ldout$semat
#'
#' @export
ldfast <- function(gp, type = c("r", "r2", "z", "D", "Dprime"), se = TRUE) {
  stopifnot(inherits(gp, "array"))
  stopifnot(length(dim(gp)) == 3)
  stopifnot(is.logical(se))
  stopifnot(length(se) == 1)
  type <- match.arg(type)

  ## Just return the LD if asked ---------------------------------------------
  if (!se) {
    return(ldfast_justmean(gp = gp, type = type))
  }

  ## Otherwise, go the slow way ----------------------------------------------
  nsnp <- dim(gp)[[1]]
  nind <- dim(gp)[[2]]
  ploidy <- dim(gp)[[3]] - 1

  ldmat <- matrix(NA_real_, ncol = nsnp, nrow = nsnp)
  semat <- matrix(NA_real_, ncol = nsnp, nrow = nsnp)
  rr    <- rep(NA_real_, length.out = nsnp)

  if (type == "D") {
    ldfast_calc(cormat = ldmat, semat = semat, rr = rr, gp = gp, type = "a")
  } else if (type == "Dprime") {
    ldfast_calc(cormat = ldmat, semat = semat, rr = rr, gp = gp, type = "c")
  } else {
    ldfast_calc(cormat = ldmat, semat = semat, rr = rr, gp = gp, type = "b")
  }

  if (type == "r2") {
    semat <- semat * abs(ldmat) * 2
    ldmat <- ldmat ^ 2
    rr <- rr ^ 2
  } else if (type == "z") {
    semat <- semat / (1 - ldmat ^ 2)
    ldmat <- atanh(ldmat)
  }

  return(list(ldmat = ldmat, rr = rr, semat = semat))
}

#' Same as ldfast, but just calculate the ld, not the se
#'
#' @inheritParams ldfast
#' @param shrinkrr A logical. Should we use adaptive shrinkage to shrink
#'     the reliability ratios (\code{TRUE}) or keep the raw reliability
#'     ratios (\code{FALSE}). Defaults to \code{TRUE}.
#'
#' @return A list with two elements. The mean matrix and the reliability ratio.
#'
#' @author David Gerard
#'
#' @examples
#' data("gp") # posterior probs
#' ldout <- ldfast_justmean(gp, "r", TRUE)
#'
#' @noRd
ldfast_justmean <- function(gp,
                            type = c("r", "r2", "z", "D", "Dprime"),
                            shrinkrr = TRUE) {
  stopifnot(inherits(gp, "array"))
  stopifnot(length(dim(gp)) == 3)
  stopifnot(is.logical(shrinkrr))
  stopifnot(length(shrinkrr) == 1)
  type <- match.arg(type)

  nsnp <- dim(gp)[[1]]
  nind <- dim(gp)[[2]]
  ploidy <- dim(gp)[[3]] - 1

  pm_mat <- matrix(NA_real_, nrow = nsnp, ncol = nind)
  pv_mat <- matrix(NA_real_, nrow = nsnp, ncol = nind)
  fill_pm(pm = pm_mat, gp = gp)
  fill_pv(pv = pv_mat, pm = pm_mat, gp = gp)

  varx <- matrixStats::rowVars(x = pm_mat, na.rm = TRUE)
  muy <- rowMeans(x = pv_mat, na.rm = TRUE)

  ## Calculate reliability ratios
  rr <- (muy + varx) / varx
  if (shrinkrr) {
    amom <- abind::abind(pm_mat, pm_mat^2, pv_mat, along = 3)
    mbar <- apply(X = amom, MARGIN = 3, FUN = rowMeans, na.rm = TRUE)
    covarray <- apply(X = amom, MARGIN = 1, FUN = stats::cov, use = "pairwise.complete.obs")
    dim(covarray) <- c(3, 3, nsnp)
    gradmat <- matrix(NA_real_, nrow = nsnp, ncol = 3)
    gradmat[, 1] <-
      -2 * mbar[, 1] / (mbar[, 3] + mbar[, 2] - mbar[, 1]^2) +
      2 * mbar[, 1] / (mbar[, 2] - mbar[, 1]^2)
    gradmat[, 2] <-
      1 / (mbar[, 3] + mbar[, 2] - mbar[, 1]^2) -
      1 / (mbar[, 2] - mbar[, 1]^2)
    gradmat[, 3] <-
      1 / (mbar[, 3] + mbar[, 2] - mbar[, 1]^2)
    svec <- rep(NA_real_, nsnp) # Variances for log rr
    nvec <- rowSums(!is.na(pm_mat)) # missing values
    for (i in seq_len(nsnp)) {
      svec[[i]] <- gradmat[i, , drop = FALSE] %*% covarray[, , i, drop = TRUE] %*% t(gradmat[i, , drop = FALSE])
    }
    svec <- svec / nvec
    svec[svec < 0] <- 0
    lvec <- log(rr)
    modest <- modeest::hsm(x = lvec)
    ashout <- ashr::ash(betahat = lvec,
                        sebetahat = sqrt(svec),
                        mixcompdist = "uniform",
                        mode = modest)
    rr <- exp(ashr::get_pm(a = ashout) / 2) ## divide by 2 for square root
  } else {
    rr <- sqrt(rr)
  }

  ## rr should now be for correlation, not covariance

  ## calculate correlation
  ldmat <- rr * stats::cor(t(pm_mat), use = "pairwise.complete.obs") * rep(rr, each = nsnp)
  if (type != "Dprime") {
    ldmat[ldmat > 1] <- 1
    ldmat[ldmat < -1] <- -1
  }

  if (type == "D") {
    rr <- rr ^ 2
    sdvec <- sqrt(muy + varx)
    ldmat <- (sdvec * ldmat * rep(sdvec, each = nsnp)) / ploidy
  } else if (type == "Dprime") {
    rr <- rr ^ 2
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
    diag(ldmat) <- ploidy
  } else if (type == "r2") {
    rr <- rr ^ 2
    ldmat <- ldmat ^ 2
  } else if (type == "z") {
    ldmat <- atanh(ldmat)
  } else {
    ## do nothing
  }

  return(list(ldmat = ldmat, rr = rr))
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
