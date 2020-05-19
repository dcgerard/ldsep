###################
## Functions to estimate LD directly from genotypes
###################

#' Estimate pair-wise LD directly from the genotypes.
#'
#' Given genotype (allele dosage) data for each individual at a pair of
#' loci, this function will calculate the maximum likelihood estimates
#' and their corresponding asymptotic standard errors of a variety of
#' measures of linkage disequilibrium (LD): D, D', and the
#' squared correlation. This function can be used for both
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
#' The resulting standard errors are the square roots of the inverse of the
#' negative Fisher-information. This is from standard maximum likelihood
#' theory. The Fisher-information is known to be biased low, so the actual
#' standard errors are probably a little bigger for small n (n < 100).
#'
#' The standard error estimate of the squared Pearson correlation is not
#' valid when r^2 = 0.
#'
#' The standard error estimate of D' is not valid when
#' \itemize{
#'   \item{D' = 0,}
#'   \item{D' < 0 and pA pB = (1 - pA) (1 - pB), or}
#'   \item{D' > 0 and pA (1 - pB) = pA (1 - pB).}
#' }
#'
#' @param ga A vector of counts, containing the genotypes for each individual
#'     at the first locus.
#' @param gb A vector of counts, containing the genotypes for each individual
#'     at the second locus.
#' @param K The ploidy of the species. Assumed the same for all individuals.
#' @param reltol The relative tolerance for the stopping criterion.
#' @param lang Should we use the R interface for optim (\code{"R"}) or the
#'     C++ interface for optim through the roptim package (\code{"C++"})?
#'
#' @author David Gerard
#'
#' @return A vector with the following elements:
#' \describe{
#'   \item{\code{D}}{The MLE of D.}
#'   \item{\code{D_se}}{The standard error of the estimate of D.}
#'   \item{\code{Dprime}}{The MLE of D'.}
#'   \item{\code{Dprime_se}}{The standard error of the estimate of D'.}
#'   \item{\code{r2}}{The MLE of the squared Pearson correlation.}
#'   \item{\code{r2_se}}{The standard error of the estimate of the
#'       squared Pearson correlation.}
#'   \item{\code{r}}{he MLE of the Pearson correlation.}
#'   \item{\code{r_se}}{The standard error of the estimate of the
#'       Pearson correlation.}
#'   \item{\code{p_ab}}{The estimated haplotype frequency of ab.}
#'   \item{\code{p_Ab}}{The estimated haplotype frequency of Ab.}
#'   \item{\code{p_aB}}{The estimated haplotype frequency of aB.}
#'   \item{\code{p_AB}}{The estimated haplotype frequency of AB.}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 10
#' K <- 6
#'
#' ## If you give ldest vectors, it assumes you are using genotypes
#' ga <- sample(0:K, 100, TRUE)
#' gb <- sample(0:K, 100, TRUE)
#' head(ga)
#' head(gb)
#' ldout <- ldest(ga = ga, gb = gb, K = K)
#' ldout
#'
#' ## if you give ldest matrices, it assumes you are using genotype likelihoods.
#' gamat <- t(sapply(ga, stats::dnorm, x = 0:K, sd = 2, log = TRUE))
#' gbmat <- t(sapply(gb, stats::dnorm, x = 0:K, sd = 2, log = TRUE))
#' head(gamat)
#' head(gbmat)
#' ldout <- ldest(ga = gamat, gb = gbmat, K = K)
#' ldout
#'
#' @export
ldest <- function(ga, gb, K, reltol = 10^-8, lang = c("C++", "R")) {

  stopifnot(length(K) == 1)
  lang <- match.arg(lang)
  if (is.vector(ga) & is.vector(gb)) {
    stopifnot(length(ga) == length(gb))
    stopifnot(ga >= 0, ga <= K)
    stopifnot(gb >= 0, gb <= K)
    using = "genotypes"
  } else if (is.matrix(ga) & is.matrix(gb)) {
    stopifnot(dim(ga) == dim(gb))
    stopifnot(K + 1 == ncol(ga))
    using = "likelihoods"
  } else {
    stop("ldest: ga and gb must either both be vectors or both be matrices.")
  }


  ## find MLE of proportions ---------
  inity <- rep(0, 3)
  if (lang == "R" & using == "genotypes") {
    oout <- stats::optim(par     = inity,
                         fn      = llike_geno,
                         gr      = dllike_geno_dpar,
                         method  = "BFGS",
                         control = list(fnscale = -1, reltol = reltol),
                         hessian = TRUE,
                         gA      = ga,
                         gB      = gb,
                         K       = K)
  } else if (using == "genotypes") {
    oout <- optimize_genocor(par    = inity,
                             gA     = ga,
                             gB     = gb,
                             K      = K,
                             reltol = reltol)
  } else {
    oout <- stats::optim(par     = inity,
                         fn      = llike_genolike,
                         gr      = dllike_genolike_dpar,
                         method  = "BFGS",
                         control = list(fnscale = -1, reltol = reltol),
                         hessian = TRUE,
                         pgA      = ga,
                         pgB      = gb)

  }

  ## Get estimates -------------------
  phat <- real_to_simplex(oout$par) # (ab, Ab, aB, AB)
  pA <- phat[[2]] + phat[[4]]
  pB <- phat[[3]] + phat[[4]]
  D  <- phat[[4]] - pA * pB
  r2 <- D ^ 2 / (pA * (1 - pA) * pB * (1 - pB))
  if (D < 0) {
    Dprime <- D / min(pA * pB, (1 - pA) * (1 - pB))
  } else {
    Dprime <- D / min(pA * (1 - pB), (1 - pA) * pB)
  }

  ## Get variance estimates --------
  Hy <- oout$hessian ## Hessian of real parameters

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

  r <- sqrt(r2) * sign(D)
  r_se <- r2_se / sqrt(4 * r2)

  retvec <- c(D         = D,
              D_se      = D_se,
              Dprime    = Dprime,
              Dprime_se = Dprime_se,
              r2        = r2,
              r2_se     = r2_se,
              r         = r,
              r_se      = r_se,
              p_ab      = phat[[1]],
              p_Ab      = phat[[2]],
              p_aB      = phat[[3]],
              p_AB      = phat[[4]])

  return(retvec)
}


#' Estimate all pair-wise LD's in a collection of SNPs.
#'
#' This function will run \code{\link{ldest}()} iteratively over
#' all possible pairs of SNPs provided. Support is provided for parallelization
#' through the doParallel and foreach packages.
#'
#' @param genomat A matrix of genotypes (allele dosages). The rows index the
#'     SNPs and the columns index the individuals. That is, \code{genomat[i, j]}
#'     is the allele dosage for individual \code{j} in SNP \code{i}.
#' @param K The ploidy of the species. Assumed to be the same for all
#'     individuals at all SNPs
#' @param nc he number of computing cores to use. This should never be
#'     more than the number of cores available in your computing environment.
#'     You can determine the maximum number of available cores by running
#'     \code{parallel::detectCores()} in R. This is probably fine for a
#'     personal computer, but some environments are only
#'     able to use fewer. Ask your admins if you are unsure.
#' @param reltol The relative tolerance for the stopping criterion.
#'
#' @author David Gerard
#'
#' @examples
#' nloci <- 5
#' nind  <- 100
#' K <- 6
#' nc <- 2
#' genomat <- matrix(sample(0:K, nind * nloci, TRUE), nrow = nloci)
#' rdf <- multi_ldest_geno(genomat = genomat, K = K, nc = nc)
#'
#'
#' @export
multi_ldest_geno <- function(genomat, K, nc = 1, reltol = 10^-8) {
  stopifnot(is.matrix(genomat))
  nloci <- nrow(genomat)

  ## Register workers ----------------------------------------------------------
  if (nc == 1) {
    foreach::registerDoSEQ()
  } else {
    cl = parallel::makeCluster(nc)
    doParallel::registerDoParallel(cl = cl)
    if (foreach::getDoParWorkers() == 1) {
      stop(paste0("multi_ldest_geno: nc > 1 ",
                  "but only one core registered from ",
                  "foreach::getDoParWorkers()."))
    }
  }

  i <- 1
  outmat <- foreach::foreach(i = seq_len(nloci - 1),
                             .combine = rbind,
                             .export = c("ldest")) %dopar% {

                               for (j in (i + 1):nloci) {
                                 ldout <- ldest(ga = genomat[i, ],
                                                gb = genomat[j, ],
                                                K = K,
                                                reltol = reltol,
                                                lang = "C++")
                                 if (j == i + 1) {
                                   estmat <- matrix(NA_real_,
                                                    nrow = nloci - i,
                                                    ncol = length(ldout) + 2)
                                 }
                                 estmat[j - i, 1] <- i
                                 estmat[j - i, 2] <- j
                                 estmat[j - i, -c(1, 2)] <- ldout

                               }
                               colnames(estmat) <- c("i", "j", names(ldout))
                               estmat
                             }

  if (nc > 1) {
    parallel::stopCluster(cl)
  }

  outmat <- as.data.frame(outmat)
  class(outmat) <- c("lddf", "data.frame")
  return(outmat)
}

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


#' Format an element of \code{\link{multi_ldest_geno}()} into an
#' upper-triangular matrix.
#'
#' Formats the correlation estimates and standard errors output
#' from running \code{\link{multi_ldest_geno}()} into a more
#' conventional upper-triangular matrix.
#'
#' @param obj An object of class \code{lddf}, usually output from
#'     running \code{\link{multi_ldest_geno}()}.
#' @param element Which element in \code{obj} should we format into an
#'     upper-triangular matrix?
#'
#' @author David Gerard
#'
#' @export
format_lddf <- function(obj,
                        element = c("r2",
                                    "r2_se",
                                    "D",
                                    "D_se",
                                    "Dprime",
                                    "Dprime_se",
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
  return(cormat)
}


