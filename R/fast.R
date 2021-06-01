###################
## FAST LD Correction
###################

#' Fast bias-correction for LD Estimation
#'
#' Estimates the reliability ratios from posterior marginal moments and uses
#' these to correct the biases in linkage disequilibrium estimation
#' caused by genotype uncertainty. These methods are described in
#' Gerard (2021).
#'
#' @section Details:
#'
#' Returns consistent and bias-corrected estimates of linkage disequilibrium.
#' The usual measures of LD are implemented: D, D', r, r2, and z
#' (Fisher-z of r). These are all \emph{composite} measures of LD, not
#' haplotypic measures of LD (see the description in \code{\link{ldest}()}).
#' They are always appropriate measures of association
#' between loci, but only correspond to haplotypic measures of LD when
#' Hardy-Weinberg equilibrium is fulfilled in autopolyploids.
#'
#' Calculating standard errors and performing hierarchical shrinkage of the
#' reliability ratios are both rather slow operations compared to just
#' raw method-of-moments based estimation for LD. If you don't need
#' standard errors, you can double your speed by setting
#' \code{se = FALSE}. It is not recommended that you disable the
#' hierarchical shrinkage.
#'
#' Due to sampling variability, the estimates sometime lie outside of the
#' theoretical boundaries of the parameters being estimated. In such cases,
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
#' \deqn{\sqrt{(a1 + b1)/a1}\sqrt{(a2 + b2)/a2}r.}
#' All other LD calculations are based on this equation. In particular,
#' the estimated genotype variances at loci 1 and 2 are
#' \eqn{a1 + b1} and \eqn{a2 + b2}, respectively, which can be
#' used to calculate D and D'.
#'
#' The reliability ratio for SNP i is defined by \eqn{(ai + bi)/ai}.
#' By default, we apply \code{\link[ashr]{ash}()} (Stephens, 2016)
#' to the log of these reliability ratios before adjusting the
#' Pearson correlation. Standard errors are required before using
#' \code{\link[ashr]{ash}()}, but these are easily obtained
#' using the central limit theorem and the delta-method.
#'
#' @param gp A three-way array with dimensions SNPs by individuals by dosage.
#'     That is, \code{gp[i, j, k]} is the posterior probability of
#'     dosage \code{k-1} for individual \code{j} at SNP \code{i}.
#' @param type What LD measure should we estimate?
#'     \describe{
#'       \item{\code{"r"}}{The Pearson correlation.}
#'       \item{\code{"r2"}}{The squared Pearson correlation.}
#'       \item{\code{"z"}}{The Fisher-z transformed Pearson correlation.}
#'       \item{\code{"D"}}{The LD coefficient.}
#'       \item{\code{"Dprime"}}{The standardized LD coefficient.}
#'     }
#'     Note that these are all \emph{composite} measures of LD (see
#'     the description in \code{\link{ldest}()}).
#' @param se Should we also return a matrix of standard errors (\code{TRUE})
#'     or not (\code{FALSE})? It is faster to not return standard errors.
#'     Defaults to \code{TRUE}.
#' @param shrinkrr A logical. Should we use adaptive shrinkage
#'     (Stephens, 2016) to shrink the reliability ratios (\code{TRUE})
#'     or keep the raw reliability ratios (\code{FALSE}). Defaults
#'     to \code{TRUE}.
#' @param thresh A logical. Should we apply an upper bound on the reliability
#'     ratios (\code{TRUE}) or not (\code{FALSE}).
#' @param upper The upper bound on the reliability ratios if
#'     \code{thresh = TRUE}. The default is a generous 10.
#' @param mode A character. Only applies if \code{shrinkrr = TRUE}. When using
#'     hierarchical shrinkage on the log of the reliability ratios, should
#'     we use zero as the mode (\code{mode = "zero"}) or estimate it using
#'     the procedure of Robertson and Cryer (1974)
#'     (\code{mode = "estimate"})?
#'
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{gl_to_gp}()}}{Normalize genotype likelihoods to
#'       posterior probabilities using naive uniform prior.}
#'   \item{\code{\link[ashr]{ash}()}}{Function used to perform hierarchical
#'       shrinkage on the log of the reliability ratios.}
#'   \item{\code{\link{ldest}()}, \code{\link{mldest}()}, \code{\link{sldest}()}}{Maximum likelihood estimation of linkage disequilibrium.}
#' }
#'
#' @return A list with some or all of the following elements:
#' \describe{
#'   \item{\code{ldmat}}{The bias-corrected LD matrix.}
#'   \item{\code{rr}}{The estimated reliability ratio for each SNP. This
#'       is the multiplicative factor applied to the naive LD estimate
#'       for each SNP.}
#'   \item{\code{rr_raw}}{The raw reliability ratios (for the covariance,
#'       not the correlation). Only returned if \code{shrinkrr = TRUE}.}
#'   \item{\code{rr_se}}{The standard errors for the \emph{log}-raw
#'       reliability ratios for each SNP. That is, we have
#'       sd(log(rr_raw)) ~ rr_se. Only returned if \code{shrinkrr = TRUE}.}
#'   \item{\code{semat}}{A matrix of standard errors of the corresponding
#'       estimators of LD.}
#' }
#'
#' @references
#' \itemize{
#'   \item{Gerard, David. Scalable Bias-corrected Linkage Disequilibrium Estimation Under Genotype Uncertainty. \emph{bioRxiv}, 2021. \doi{10.1101/2021.02.08.430270}}
#'   \item{T. Robertson and J. D. Cryer. An iterative procedure for estimating the mode. \emph{Journal of the American Statistical Association}, 69(348):1012–1016, 1974. \doi{10.1080/01621459.1974.10480246}.}
#'   \item{M. Stephens. False discovery rates: a new deal. \emph{Biostatistics}, 18(2):275–294, 10 2016. ISSN 1465-4644 \doi{10.1093/biostatistics/kxw041}.}
#' }
#'
#' @author David Gerard
#'
#' @examples
#' data("gp")
#'
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
ldfast <- function(gp,
                   type = c("r", "r2", "z", "D", "Dprime"),
                   shrinkrr = TRUE,
                   se = TRUE,
                   thresh = TRUE,
                   upper = 10,
                   mode = c("zero", "estimate")) {
  ## Check input -------------------------------------------------------------
  stopifnot(inherits(gp, "array"))
  stopifnot(length(dim(gp)) == 3)
  stopifnot(is.logical(shrinkrr))
  stopifnot(length(shrinkrr) == 1)
  stopifnot(is.logical(se))
  stopifnot(length(se) == 1)
  stopifnot(is.logical(thresh))
  stopifnot(length(thresh) == 1)
  stopifnot(is.numeric(upper))
  stopifnot(length(numeric) == 1)
  type <- match.arg(type)
  mode <- match.arg(mode)

  nsnp <- dim(gp)[[1]]
  nind <- dim(gp)[[2]]
  ploidy <- dim(gp)[[3]] - 1

  ## Calculate posterior moments ----------------------------------------------
  pm_mat <- matrix(NA_real_, nrow = nsnp, ncol = nind)
  pv_mat <- matrix(NA_real_, nrow = nsnp, ncol = nind)
  fill_pm(pm = pm_mat, gp = gp)
  fill_pv(pv = pv_mat, pm = pm_mat, gp = gp)

  varx <- matrixStats::rowVars(x = pm_mat, na.rm = TRUE)
  muy <- rowMeans(x = pv_mat, na.rm = TRUE)

  ## Calculate reliability ratios ---------------------------------------------
  ## after this step rr should be for correlation, *not* for covariance
  rr_raw <- (muy + varx) / varx
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
    svec <- rep(NA_real_, nsnp) # Variances for log rr_raw
    nvec <- rowSums(!is.na(pm_mat)) # missing values
    for (i in seq_len(nsnp)) {
      svec[[i]] <- gradmat[i, , drop = FALSE] %*% covarray[, , i, drop = TRUE] %*% t(gradmat[i, , drop = FALSE])
    }
    svec <- svec / nvec
    svec[svec < 0] <- 0
    if (thresh) {
      svec[rr_raw > upper] <- Inf
    }
    lvec <- log(rr_raw)
    if (mode == "estimate") {
      modest <- modeest::hsm(x = lvec)
    } else if (mode == "zero") {
      modest <- 0
    }
    ashout <- ashr::ash(betahat = lvec,
                        sebetahat = sqrt(svec),
                        mixcompdist = "uniform",
                        mode = modest)
    rr <- exp(ashr::get_pm(a = ashout) / 2) ## divide by 2 for square root
  } else {
    if (thresh) {
      rr <- rr_raw
      rr[rr > upper] <- stats::median(rr)
    }
    rr <- sqrt(rr)
  }

  ## Calculate Pearson correlation --------------------------------------------
  ldmat <- rr * stats::cor(t(pm_mat), use = "pairwise.complete.obs") * rep(rr, each = nsnp)
  if (type != "Dprime") {
    if (se) {
      which_truncate <- (ldmat > 1) | (ldmat < -1)
    }
    ldmat[ldmat > 1] <- 1
    ldmat[ldmat < -1] <- -1
  }

  ## LD adjustment based on user selection ------------------------------------
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


    if (se) {
      which_truncate <- (ldmat > ploidy) | (ldmat < -ploidy)
    }
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

  ## Return list --------------------------------------------------------------
  retlist <- list(ldmat = ldmat, rr = rr)

  ## Add standard errors of rr if shrinkrr = TRUE
  if (shrinkrr) {
    retlist$rr_raw <- rr_raw
    retlist$rr_se <- sqrt(svec)
  }

  ## Standard error calculations if option ------------------------------------
  if (se) {
    if (type == "D") {
      semat <- secalc(gp = gp, pm_mat = pm_mat, pv_mat = pv_mat, type = "a")
    } else if (type %in% c("r", "r2", "z")) {
      semat <- secalc(gp = gp, pm_mat = pm_mat, pv_mat = pv_mat, type = "b")
      if (type == "r2") {
        semat <- 2 * semat * sqrt(ldmat)
      } else if (type == "z") {
        semat <- semat / (1 - tanh(ldmat) ^ 2)
      }
    } else if (type == "Dprime") {
      semat <- secalc(gp = gp, pm_mat = pm_mat, pv_mat = pv_mat, type = "c")
    }

    semat[which_truncate] <- NA_real_
    retlist$semat <- semat
  }

  return(retlist)
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
