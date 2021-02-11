#' Pairwise LD estimation in polyploids.
#'
#' Estimates either haplotypic or composite measures of LD using either
#' genotypes are genotype likelihoods via maximum likelihood.
#' The usual measures of LD are estimated (D, D', and r) along with
#' the Fisher-z transformation of r (called "z"). All estimates
#' are returned with standard errors. See Gerard (2021) for details.
#'
#' @section Haplotypic LD:
#'
#' This section describes the methods used when \code{type = "hap"} is
#' selected.
#'
#' Haplotypic LD measures the association
#' between two loci on the same haplotype. When haplotypes are known, estimating
#' haplotypic LD is simple using just the haplotypic frequencies.
#'
#' When haplotypes are not known, we can still estimate haplotypic frequencies
#' using the genotypes or genotype likelihoods
#' \emph{in autopolyploids as long as Hardy-Weinberg equilibrium (HWE) is satisfied}. We do
#' this via maximum likelihood using gradient ascent. Gradient ascent is
#' performed over the unconstrained parameterization of the 3-simplex from
#' Betancourt (2012). The estimated haplotype frequencies are then used to
#' estimate haplotypic LD.
#'
#' Standard errors are provided using standard maximum likelihood theory.
#' In brief, the Hessian matrix of the log-likelihood is calculated at
#' the MLE's of the haplotype frequencies. The negative inverse of this
#' Hessian matrix is approximately the covariance matrix of the MLE's of the
#' haplotype frequencies. Since all haplotypic LD measures are functions
#' of the haplotype frequencies, we use the delta-method to obtain
#' the standard errors for each LD estimate.
#'
#' A Dirichlet(2,2,2,2) prior is placed over the frequencies of
#' haplotypes 00, 01, 10, and 11. This corresponds to the "add two" rule
#' of Agresti (1998). You can change this prior via the \code{pen} argument.
#'
#' When you either do not have autopolyploids or when HWE is \emph{not}
#' satisfied, then the estimates using \code{type = "hap"}
#' are nonsensical. However, the composite measures of LD are still
#' applicable (see below).
#'
#' @section Composite LD:
#'
#' This section describes the methods used when \code{type = "comp"} is
#' selected.
#'
#' When HWE is not satisfied, haplotype frequencies are not estimable. However,
#' measures of association between two loci are still estimable. These
#' associations may be caused by LD either on the same haplotype or between
#' different haplotypes. Cockerham and Weir (1977) thus called such measures
#' "composite" measures of LD.
#'
#' When the genotypes are known, these composite measures have simple
#' correspondences to well-known statistical measures of association.
#' D is the covariance of genotypes between loci divided by the ploidy.
#' r is the Pearson correlation of genotypes. D' is D divided by a
#' term that involves only mean genotypes.
#'
#' When genotypes are not known, we estimate the joint genotype frequencies
#' and use these to estimate the composite measures of LD using
#' genotype likelihoods. The distribution of genotypes is assumed to
#' either follow a proportional bivariate normal model (by default) or
#' a general categorical model.
#'
#' These estimates of composite measures of LD estimate the haplotypic
#' measures of LD when HWE is fulfilled, but are still applicable when HWE
#' is not fulfilled.
#'
#' When genotypes are known, standard errors are calculated using standard
#' moment-based approaches. When genotypes are not known, standard
#' errors are calculated using standard maximum likelihood theory,
#' same as for the haplotypic LD estimates (see above), or using
#' a bootstrap.
#'
#' @param ga One of two possible inputs:
#'     \enumerate{
#'         \item{A vector of counts, containing the genotypes for each
#'               individual at the first locus. When \code{type = "comp"},
#'               the vector of genotypes may be continuous (e.g. the
#'               posterior mean genotype).}
#'         \item{A matrix of genotype log-likelihoods at the first locus.
#'               The rows index the individuals and the columns index
#'               the genotypes. That is \code{ga[i, j]} is the genotype
#'               likelihood of individual \code{i} for genotype \code{j-1}.}
#'         }
#' @param gb One of two possible inputs:
#'     \enumerate{
#'         \item{A vector of counts, containing the genotypes for each
#'               individual at the second locus. When \code{type = "comp"},
#'               the vector of genotypes may be continuous (e.g. the
#'               posterior mean genotype).}
#'         \item{A matrix of genotype log-likelihoods at the second locus.
#'               The rows index the individuals and the columns index
#'               the genotypes. That is \code{gb[i, j]} is the genotype
#'               likelihood of individual \code{i} for genotype \code{j-1}.}
#'         }
#' @param K The ploidy of the species. Assumed to be the same for all
#'     individuals.
#' @param type The type of LD to calculate. The available types are
#'     haplotypic LD (\code{type = "hap"}) or composite LD
#'     (\code{type = "comp"}). Haplotypic LD is only appropriate for
#'     autopolyploids when the individuals are in Hardy-Weinberg
#'     equilibrium (HWE). The composite
#'     measures of LD are always applicable, and consistently estimate the
#'     usual measures of LD when HWE is fulfilled in autopolyploids.
#'     However, when HWE is not fulfilled, interpreting the
#'     composite measures of LD could be a little tricky.
#' @param model When \code{type = "comp"} and using genotype likelihoods,
#'     should we use the proportional
#'     bivariate normal model to estimate the genotype distribution
#'     (\code{model = "norm"}), or the general categorical distribution
#'     (\code{model = "flex"})? Defaults to \code{"norm"}.
#' @param pen The penalty to be applied to the likelihood. You can think about
#'     this as the prior sample size. Should be greater than 1. Does not
#'     apply if \code{model = "norm"}, \code{type = "comp"}, and using
#'     genotype likelihoods. Also does not apply when \code{type = "comp"}
#'     and using genotypes.
#' @param se A logical. Should we calculate standard errors (\code{TRUE}) or
#'     not (\code{FALSE}). Calculating standard errors can be really slow
#'     when \code{type = "comp"}, \code{model = "flex"}, and when using
#'     genotype likelihoods. Otherwise, standard error calculations
#'     should be pretty fast.
#'
#' @return A vector with some or all of the following elements:
#' \describe{
#'   \item{\code{D}}{The estimate of the LD coefficient.}
#'   \item{\code{D_se}}{The standard error of the estimate of
#'       the LD coefficient.}
#'   \item{\code{r2}}{The estimate of the squared Pearson correlation.}
#'   \item{\code{r2_se}}{The standard error of the estimate of the
#'       squared Pearson correlation.}
#'   \item{\code{r}}{The estimate of the Pearson correlation.}
#'   \item{\code{r_se}}{The standard error of the estimate of the
#'       Pearson correlation.}
#'   \item{\code{Dprime}}{The estimate of the standardized LD
#'       coefficient. When \code{type} = "comp", this corresponds
#'       to the standardization where we fix allele frequencies.}
#'   \item{\code{Dprime_se}}{The standard error of \code{Dprime}.}
#'   \item{\code{Dprimeg}}{The estimate of the standardized LD
#'       coefficient. This corresponds to the standardization where
#'       we fix genotype frequencies.}
#'   \item{\code{Dprimeg_se}}{The standard error of \code{Dprimeg}.}
#'   \item{\code{z}}{The Fisher-z transformation of \code{r}.}
#'   \item{\code{z_se}}{The standard error of the Fisher-z
#'       transformation of \code{r}.}
#'   \item{\code{p_ab}}{The estimated haplotype frequency of ab.
#'       Only returned if estimating the haplotypic LD.}
#'   \item{\code{p_Ab}}{The estimated haplotype frequency of Ab.
#'       Only returned if estimating the haplotypic LD.}
#'   \item{\code{p_aB}}{The estimated haplotype frequency of aB.
#'       Only returned if estimating the haplotypic LD.}
#'   \item{\code{p_AB}}{The estimated haplotype frequency of AB.
#'       Only returned if estimating the haplotypic LD.}
#'   \item{\code{q_ij}}{The estimated frequency of genotype i at locus 1
#'       and genotype j at locus 2. Only returned if estimating the
#'       composite LD.}
#'   \item{\code{n}}{The number of individuals used to estimate pairwise LD.}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 100 # sample size
#' K <- 6 # ploidy
#'
#' ## generate some fake genotypes when LD = 0.
#' ga <- stats::rbinom(n = n, size = K, prob = 0.5)
#' gb <- stats::rbinom(n = n, size = K, prob = 0.5)
#' head(ga)
#' head(gb)
#'
#' ## generate some fake genotype likelihoods when LD = 0.
#' gamat <- t(sapply(ga, stats::dnorm, x = 0:K, sd = 1, log = TRUE))
#' gbmat <- t(sapply(gb, stats::dnorm, x = 0:K, sd = 1, log = TRUE))
#' head(gamat)
#' head(gbmat)
#'
#' ## Haplotypic LD with genotypes
#' ldout1 <- ldest(ga = ga,
#'                 gb = gb,
#'                 K = K,
#'                 type = "hap")
#' head(ldout1)
#'
#' ## Haplotypic LD with genotype likelihoods
#' ldout2 <- ldest(ga = gamat,
#'                 gb = gbmat,
#'                 K = K,
#'                 type = "hap")
#' head(ldout2)
#'
#' ## Composite LD with genotypes
#' ldout3 <- ldest(ga = ga,
#'                 gb = gb,
#'                 K = K,
#'                 type = "comp")
#' head(ldout3)
#'
#' ## Composite LD with genotype likelihoods and normal model
#' ldout4 <- ldest(ga = gamat,
#'                 gb = gbmat,
#'                 K = K,
#'                 type = "comp",
#'                 model = "norm")
#' head(ldout4)
#'
#' ## Composite LD with genotype likelihoods and general categorical model
#' ldout5 <- ldest(ga = gamat,
#'                 gb = gbmat,
#'                 K = K,
#'                 type = "comp",
#'                 model = "flex",
#'                 se = FALSE)
#' head(ldout5)
#'
#' ldout1[["D"]]
#' ldout2[["D"]]
#' ldout3[["D"]]
#' ldout4[["D"]]
#' ldout5[["D"]]
#'
#' @author David Gerard
#'
#' @seealso
#' \describe{
#'   \item{\code{\link{ldfast}()}}{Fast, moment-based approach to LD estimation
#'       that also accounts for genotype uncertainty.}
#'   \item{\code{\link{mldest}()}}{For calculating pairwise LD among all
#'       pairs of a collection of SNPs.}
#'   \item{\code{\link{sldest}()}}{For calculating pairwise LD along a
#'       sliding window of SNPs.}
#'   \item{\code{\link{ldest_hap}()}}{For the function that directly estimates
#'       haplotypic LD when HWE is fulfilled.}
#'   \item{\code{\link{ldest_comp}()}}{For the function that directly
#'       estimates composite LD.}
#' }
#'
#' @references
#' \itemize{
#'   \item{Agresti, Alan, and Brent A. Coull. "Approximate is better than
#'         "exact" for interval estimation of binomial proportions."
#'         \emph{The American Statistician} 52, no. 2 (1998): 119-126.}
#'   \item{Betancourt, Michael. "Cruising the simplex: Hamiltonian Monte
#'         Carlo and the Dirichlet distribution." In
#'         \emph{AIP Conference Proceedings 31st}, vol. 1443, no. 1,
#'         pp. 157-164. American Institute of Physics, 2012.}
#'   \item{Cockerham, C. Clark, and B. S. Weir. "Digenic descent measures
#'         for finite populations." \emph{Genetics Research} 30, no. 2 (1977):
#'         121-147.}
#'   \item{Gerard, David. "Pairwise Linkage Disequilibrium Estimation
#'         for Polyploids." \emph{Molecular Ecology Resources}.
#'         Accepted Author Manuscript. (2021)
#'         \href{https://doi.org/10.1111/1755-0998.13349}{doi:10.1111/1755-0998.13349}}
#' }
#'
#' @export
ldest <- function(ga,
                  gb,
                  K,
                  se = TRUE,
                  type = c("hap", "comp"),
                  model = c("norm", "flex"),
                  pen = ifelse(type == "hap", 2, 1)) {
  type <- match.arg(type)
  model <- match.arg(model)
  stopifnot(is.logical(se))

  if (type == "hap") {
    retvec <- ldest_hap(ga  = ga,
                        gb  = gb,
                        K   = K,
                        pen = pen,
                        se  = se)
  } else {
    retvec <- ldest_comp(ga    = ga,
                         gb    = gb,
                         K     = K,
                         pen   = pen,
                         se    = se,
                         model = model)
  }

  return(retvec)
}


#' Wrapper for optim() conditional on input
#'
#' @inheritParams ldest
#'
#' @author David Gerard
#'
#' @noRd
find_mle <- function(ga,
                     gb,
                     K,
                     reltol = 10^-8,
                     using = c("genotypes", "likelihoods"),
                     pen = 2,
                     grid_init = FALSE) {
  using <- match.arg(using)
  if (grid_init) {
    inity_list <- lapply(
      list(
        c(0.25, 0.25, 0.25, 0.25),
        c(0.85, 0.05, 0.05, 0.05),
        c(0.05, 0.85, 0.05, 0.05),
        c(0.05, 0.05, 0.85, 0.05),
        c(0.05, 0.05, 0.05, 0.85),
        c(0.45, 0.45, 0.05, 0.05),
        c(0.45, 0.05, 0.45, 0.05),
        c(0.45, 0.05, 0.05, 0.45),
        c(0.05, 0.45, 0.45, 0.05),
        c(0.05, 0.45, 0.05, 0.45),
        c(0.05, 0.05, 0.45, 0.45),
        c(0.04, 0.32, 0.32, 0.32),
        c(0.32, 0.04, 0.32, 0.32),
        c(0.32, 0.32, 0.04, 0.32),
        c(0.32, 0.32, 0.32, 0.04)
      ), simplex_to_real)
  } else {
    if (using == "genotypes") {
      mafA <- mean(ga) / K
      mafB <- mean(gb) / K
    } else {
      mafA <- mean(apply(X = ga, MARGIN = 1, FUN = which.max) - 1) / K
      mafB <- mean(apply(X = gb, MARGIN = 1, FUN = which.max) - 1) / K
    }
    inity_list <- list(
      simplex_to_real(c((1 - mafA) * (1 - mafB),
                        mafA * (1 - mafB),
                        (1 - mafA) * mafB,
                        mafA * mafB))
    )
  }

  oldlike <- -Inf
  for (i in seq_along(inity_list)) {
    inity <- inity_list[[i]]
    if (using == "genotypes") {
      oout <- stats::optim(par     = inity,
                           fn      = llike_geno,
                           gr      = dllike_geno_dpar,
                           method  = "BFGS",
                           control = list(fnscale = -1, reltol = reltol),
                           hessian = TRUE,
                           gA      = ga,
                           gB      = gb,
                           K       = K,
                           alpha   = rep(pen, 4))
    } else {
      oout <- stats::optim(par     = inity,
                           fn      = llike_genolike,
                           gr      = dllike_genolike_dpar,
                           method  = "BFGS",
                           control = list(fnscale = -1, reltol = reltol),
                           hessian = TRUE,
                           pgA     = ga,
                           pgB     = gb,
                           alpha   = rep(pen, 4))
    }
    if (oout$value > oldlike) {
      oldlike <- oout$value
      oout_best <- oout
    }
  }
  return(oout_best)
}

#' Converts a vector of probabilities of LD estimates
#'
#' @param phat A vector of haplotype frequencies in the order (ab, Ab, aB, AB)
#'
#' @author David Gerard
#'
#' @noRd
convert_ld <- function(phat) {
  pA <- phat[[2]] + phat[[4]]
  pB <- phat[[3]] + phat[[4]]
  D  <- phat[[4]] - pA * pB
  TOL <- sqrt(.Machine$double.eps)
  if (abs(pA) > TOL &&
      abs(pA - 1) > TOL &&
      abs(pB) > TOL &&
      abs(pB - 1) > TOL) {
    ## everthing is perfect
    r2 <- D ^ 2 / (pA * (1 - pA) * pB * (1 - pB))
    if (D < 0) {
      Dprime <- D / min(pA * pB, (1 - pA) * (1 - pB))
    } else {
      Dprime <- D / min(pA * (1 - pB), (1 - pA) * pB)
    }
    r <- sqrt(r2) * sign(D)
    z <- atanh(r)
  } else {
    ## bad
    r2     <- 0
    Dprime <- 0
    r      <- 0
    z      <- 0
  }

  retvec <- c(pA = pA,
              pB = pB,
              D = D,
              r2 = r2,
              Dprime = Dprime,
              r = r,
              z = z)
  return(retvec)
}

#' Estimate haplotypic pair-wise LD using either genotypes or genotype
#' likelihoods.
#'
#' Given genotype (allele dosage) or genotype likelihood data
#' for each individual at a pair of loci, this function will
#' calculate the maximum likelihood estimates
#' and their corresponding asymptotic standard errors of some
#' measures of linkage disequilibrium (LD): D, D', the Pearson correlation,
#' the squared Pearson correlation, and the Fisher-z transformation of the
#' Pearson correlation. This function can be used for both
#' diploids and polyploids.
#'
#' Let A and a be the reference and alternative alleles, respectively, at
#' locus 1. Let B and b be the reference and alternative alleles,
#' respectively, at locus 2. Let paa, pAb, paB, and pAB be the
#' frequencies of haplotypes ab, Ab, aB, and AB, respectively.
#' Let pA = pAb + pAB and let pB = paB + pAB
#' The \code{ldest} returns estimates of the following measures
#' of LD.
#' \itemize{
#'   \item{D: pAB - pA pB}
#'   \item{D': D / Dmax, where Dmax = min(pA pB, (1 - pA) (1 - pB)) if
#'         D < 0 and Dmax = min(pA (1 - pB), pA (1 - pB)) if D > 0}
#'   \item{r-squared: The squared Pearson correlation,
#'         r^2 = D^2 / (pA (1 - pA) pB (1 - pB))}
#'   \item{r: The Pearson correlation,
#'         r = D / sqrt(pA (1 - pA) pB (1 - pB))}
#' }
#'
#' Estimates are obtained via maximum likelihood under the assumption
#' of Hardy-Weinberg equilibrium. The likelihood is calculated by
#' integrating over the possible haplotypes for each pair of genotypes.
#'
#' The resulting standard errors are based on the square roots of the inverse of the
#' negative Fisher-information. This is from standard maximum likelihood
#' theory. The Fisher-information is known to be biased low, so the actual
#' standard errors are probably a little bigger for small n (n < 20).
#' In some cases the Fisher-information matrix is singular, and so we
#' in these cases we return a bootstrap estimate of the standard error.
#'
#' The standard error estimate of the squared Pearson correlation is not
#' valid when r^2 = 0.
#'
#' In cases where either SNP is estimated to be monoallelic
#' (\code{pA %in% c(0, 1)} or \code{pB %in% c(0, 1)}), this function
#' will return LD estimates of \code{NA}.
#'
#' @inheritParams ldest
#' @param reltol The relative tolerance for the stopping criterion.
#' @param nboot Sometimes, the MLE standard errors don't exist. So we use
#'     the bootstrap as a backup. \code{nboot} specifies the number
#'     of bootstrap iterations.
#' @param useboot A logical. Optionally, you may always use the bootstrap
#'     to estimate the standard errors (\code{TRUE}). These will be more
#'     accurate but also much slower, so this defaults to \code{FALSE}. Only
#'     applicable if using genotype likelihoods.
#' @param grid_init A logical. Should we initialize the gradient ascent
#'     at a grid of initial values (\code{TRUE}) or just initialize
#'     at one value corresponding to the simplex point
#'     \code{rep(0.25, 4)} (\code{FALSE})?
#'
#' @author David Gerard
#'
#' @inherit ldest return
#'
#' @examples
#' set.seed(1)
#' n <- 100 # sample size
#' K <- 6 # ploidy
#'
#' ## generate some fake genotypes when LD = 0.
#' ga <- stats::rbinom(n = n, size = K, prob = 0.5)
#' gb <- stats::rbinom(n = n, size = K, prob = 0.5)
#' head(ga)
#' head(gb)
#'
#' ## generate some fake genotype likelihoods when LD = 0.
#' gamat <- t(sapply(ga, stats::dnorm, x = 0:K, sd = 1, log = TRUE))
#' gbmat <- t(sapply(gb, stats::dnorm, x = 0:K, sd = 1, log = TRUE))
#' head(gamat)
#' head(gbmat)
#'
#' ## Haplotypic LD with genotypes
#' ldout1 <- ldest_hap(ga = ga,
#'                     gb = gb,
#'                     K = K)
#' head(ldout1)
#'
#' ## Haplotypic LD with genotype likelihoods
#' ldout2 <- ldest_hap(ga = gamat,
#'                     gb = gbmat,
#'                     K = K)
#' head(ldout2)
#'
#' @export
ldest_hap <- function(ga,
                      gb,
                      K,
                      reltol  = 10^-8,
                      nboot   = 100,
                      useboot = FALSE,
                      pen     = 2,
                      grid_init = FALSE,
                      se = TRUE) {

  TOL <- sqrt(.Machine$double.eps)
  stopifnot(is.logical(se))
  stopifnot(length(K) == 1,
            length(nboot) == 1,
            length(reltol) == 1,
            length(useboot) == 1,
            length(grid_init) == 1)
  stopifnot(pen > 0)
  stopifnot(is.logical(useboot))
  stopifnot(is.logical(grid_init))
  stopifnot(reltol > 0, nboot > 0, K > 0)
  if (is.vector(ga) & is.vector(gb)) {
    which_bad <- is.na(ga) | is.na(gb)
    ga <- ga[!which_bad]
    gb <- gb[!which_bad]
    stopifnot(length(ga) == length(gb))
    stopifnot(ga >= 0, ga <= K)
    stopifnot(gb >= 0, gb <= K)
    using = "genotypes"
  } else if (is.matrix(ga) & is.matrix(gb)) {
    which_bad <- apply(ga, 1, function(x) any(is.na(x))) | apply(gb, 1, function(x) any(is.na(x)))
    ga <- ga[!which_bad, , drop = FALSE]
    gb <- gb[!which_bad, , drop = FALSE]
    stopifnot(dim(ga) == dim(gb))
    stopifnot(K + 1 == ncol(ga))
    using = "likelihoods"
  } else {
    stop("ldest: ga and gb must either both be vectors or both be matrices.")
  }

  ## check for monoallelic SNPs -----------------------------------------------
  if (using == "genotypes" & ((stats::sd(ga, na.rm = TRUE) < TOL) || (stats::sd(gb, na.rm = TRUE) < TOL))) {
    retvec <- nullvec_hap()
    return(retvec)
  }


  ## find MLE of proportions ---------
  oout <- find_mle(ga     = ga,
                   gb     = gb,
                   K      = K,
                   reltol = reltol,
                   using  = using,
                   pen    = pen,
                   grid_init = grid_init)

  ## Get estimates -------------------
  phat <- real_to_simplex(oout$par) # (ab, Ab, aB, AB)
  ldestvec <- convert_ld(phat = phat)
  pA     <- ldestvec[["pA"]]
  pB     <- ldestvec[["pB"]]
  D      <- ldestvec[["D"]]
  r2     <- ldestvec[["r2"]]
  Dprime <- ldestvec[["Dprime"]]
  r      <- ldestvec[["r"]]
  z      <- ldestvec[["z"]]

  ## Get variance estimates --------
  Hy <- oout$hessian ## Hessian of real parameters

  ## test for singularity
  eval <- eigen(Hy)
  if ((any(abs(eval$values) < sqrt(.Machine$double.eps)) || useboot) & se) {
    ## Bootstrap SE
    if (using == "genotypes") {
      nind <- length(ga)
    } else {
      nind <- nrow(ga)
    }
    for (iboot in seq_len(nboot)) {
      if (using == "genotypes") {
        ga_boot <- ga[sample(seq_len(nind), size = nind, replace = TRUE)]
        gb_boot <- gb[sample(seq_len(nind), size = nind, replace = TRUE)]
      } else {
        ga_boot <- ga[sample(seq_len(nind), size = nind, replace = TRUE), ]
        gb_boot <- gb[sample(seq_len(nind), size = nind, replace = TRUE), ]
      }
      oout_boot <- find_mle(ga     = ga_boot,
                            gb     = gb_boot,
                            K      = K,
                            reltol = reltol,
                            using  = using)
      phat_boot <- real_to_simplex(oout_boot$par)
      ldestvec_boot <- convert_ld(phat = phat_boot)
      if (iboot == 1) {
        bootmat <- matrix(NA_real_,
                          ncol = length(ldestvec_boot),
                          nrow = nboot)
        colnames(bootmat) <- names(ldestvec_boot)
      }
      bootmat[iboot, ] <- ldestvec_boot
    }
    sevec <- apply(bootmat, 2, stats::sd, na.rm = TRUE)
    D_se      <- sevec[["D"]]
    r2_se     <- sevec[["r2"]]
    Dprime_se <- sevec[["Dprime"]]
    z_se      <- sevec[["z"]]
    r_se      <- sevec[["r"]]
  } else if (se) {
    ## MLE theory

    ## jacobian converting from simplex to real. The last is just a column of
    ## zeros because we use the transform from p1, p2, p3 to y1, y2, y3
    J <- dsimplex_to_real_dx(phat)[, -4]
    Hp <- t(J) %*% Hy %*% J ## Hessian of first three simplex parameters
    nHp_inv <- -solve(Hp)

    dD <- dD_dprob(prob = phat)
    D_se <- sqrt(t(dD) %*% nHp_inv %*% dD)

    dDprime <- dDprime_dprob(prob = phat)
    Dprime_se <- sqrt(t(dDprime) %*% nHp_inv %*% dDprime)

    dr2 <- dr2_dprob(prob = phat)
    r2_se <- sqrt(t(dr2) %*% nHp_inv %*% dr2)

    r_se <- r2_se / sqrt(4 * r2)

    z_se <- r_se / (1 - r2)

    # g <- log(-log(r2))
    # g_se <- r2_se / abs(r2 * log(r2))
  } else {
    D_se <- NA_real_
    r2_se <- NA_real_
    r_se <- NA_real_
    Dprime_se <- NA_real_
    z_se <- NA_real_
  }

  if (using == "genotypes") {
    nind <- length(ga)
  } else {
    nind <- nrow(ga)
  }

  retvec <- c(D         = D,
              D_se      = D_se,
              r2        = r2,
              r2_se     = r2_se,
              r         = r,
              r_se      = r_se,
              Dprime    = Dprime,
              Dprime_se = Dprime_se,
              z         = z,
              z_se      = z_se,
              p_ab      = phat[[1]],
              p_Ab      = phat[[2]],
              p_aB      = phat[[3]],
              p_AB      = phat[[4]],
              n         = nind)

  return(retvec)
}

#' The null return value when estimating haplotypic LD
#'
#' @param K the ploidy of the species
#'
#' @author David Gerard
#'
#' @noRd
nullvec_hap <- function() {
  retvec <- c(D         = NA_real_,
              D_se      = NA_real_,
              r2        = NA_real_,
              r2_se     = NA_real_,
              r         = NA_real_,
              r_se      = NA_real_,
              Dprime    = NA_real_,
              Dprime_se = NA_real_,
              z         = NA_real_,
              z_se      = NA_real_,
              p_ab      = NA_real_,
              p_Ab      = NA_real_,
              p_aB      = NA_real_,
              p_AB      = NA_real_,
              n         = NA_real_)
  return(retvec)
}
