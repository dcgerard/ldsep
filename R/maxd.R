#################
## Functions for detla' when maximizing covariance
#################

#' Find joint distribution that maximizes covariance conditional
#' on marginal distributions.
#'
#' This will run a linear program using \code{\link[lpsolve]{lp}()}
#' to maximize covariance conditional on marginals. It returns the
#' joint distribution in the form of a matrix.
#'
#' @param pA The marginal distribution of genotypes of locus 1.
#'     \code{pA[i]} is the marginal probability of genotype i-1 at
#'     locus 1.
#' @param pB The marginal distribution of genotypes of locus 2.
#'     \code{pB[j]} is the marginal probability of genotype j-1 at
#'     locus 2.
#' @param dir Should we maximize the covariance (\code{"max"}), or
#'     minimize the covariance (\code{"min"}).
#'
#' @return A matrix. The joint distribution that maximizes the
#'     covariance conditional on the marginals. Element (i, j) is
#'     the probability of genotype i-1 at locus 1 and j-1 at locus 2.
#'
#' @noRd
#'
#' @examples
#' set.seed(1)
#' K <- 6
#' pA <- stats::runif(K + 1)
#' pA <- pA / sum(pA)
#' pB <- stats::runif(K + 1)
#' pB <- pB / sum(pB)
#' maxq_from_marg(pA, pB)
#'
#' @author David Gerard
maxq_from_marg <- function(pA, pB, dir = c("max", "min")) {
  TOL <- sqrt(.Machine$double.eps)
  dir <- match.arg(dir)
  stopifnot(length(pA) == length(pB))
  stopifnot(is.numeric(pA), is.numeric(pB))
  stopifnot(abs(sum(pA) - 1) < TOL,
            abs(sum(pB) - 1) < TOL)
  K <- length(pA) - 1

  obvec <- c((0:K) %*% t(0:K))

  constmat <- matrix(NA_real_,
                     nrow = 2 * (K + 1) + 1,
                     ncol = length(obvec))
  for (i in 0:K) {
    tempmat <- matrix(0, nrow = K + 1, ncol = K + 1)
    tempmat[i + 1, ] <- 1
    constmat[i + 1, ] <- c(tempmat)
  }
  for (j in 0:K) {
    tempmat <- matrix(0, nrow = K + 1, ncol = K + 1)
    tempmat[, j + 1] <- 1
    constmat[(K + 1) + (j + 1), ] <- c(tempmat)
  }
  constmat[2 * (K + 1) + 1, ] <- 1

  dirvec <- rep("==", nrow(constmat))

  rhsvec <- c(pA, pB, 1)

  lpout <- lpSolve::lp(direction    = dir,
                       objective.in = obvec,
                       const.mat    = constmat,
                       const.dir    = dirvec,
                       const.rhs    = rhsvec)

  qmat <- matrix(lpout$solution, nrow = K + 1, ncol = K + 1)

  return(qmat)
}

#' Find joint distribution that maximizes covariance conditional
#' on allele frequencies.
#'
#' @param pA The allele frequency at locus 1.
#' @param pB The allele frequency at locus 2.
#' @param dir Should we maximize (\code{"max"}) or minimize
#'     (\code{"min"}) the covariance?
#' @param K The ploidy of the species.
#'
#' @examples
#' pA <- 0.5
#' pB <- 0.5
#' K <- 6
#' maxq_from_allelef(pA = pA, pB = pB, K = K, dir = "max")
#'
#' @noRd
#'
#' @author David Gerard
#'
maxq_from_allelef <- function(pA, pB, K, dir = c("max", "min")) {
  TOL <- sqrt(.Machine$double.eps)
  dir <- match.arg(dir)
  stopifnot(length(pA) == 1, length(pB) == 1, length(K) == 1)
  stopifnot(pA > -TOL, pA < 1 + TOL, pB > -TOL, pB < 1 + TOL)
  stopifnot(K > 0)

  obvec <- c((0:K) %*% t(0:K))

  constmat <- matrix(NA_real_, nrow = 3, ncol = length(obvec))
  constmat[1, ] <- rep((0:K) / K, times = K + 1)
  constmat[2, ] <- rep((0:K) / K, each = K + 1)
  constmat[3, ] <- 1

  dirvec <- rep("==", times = 3)
  rhsvec <- c(pA, pB, 1)

  lpout <- lpSolve::lp(direction    = dir,
                       objective.in = obvec,
                       const.mat    = constmat,
                       const.dir    = dirvec,
                       const.rhs    = rhsvec)

  qmat <- matrix(lpout$solution, nrow = K + 1, ncol = K + 1)

  return(qmat)
}

#' Get the standardized composite D'.
#'
#' This function will either standardize by the maximum covariance
#' conditional on the marginal genotype distribution, or by the
#' maximum covariance conditional on the marginal allele frequencies.
#'
#' Note that when \code{type = "allele"} and \code{constrain = FALSE},
#' the resulting D' is constrained to fall between -K and K, where
#' K is the ploidy of the species. However, under HWE, this measure is
#' equal to gametic D'. Using \code{constrain = TRUE} will result
#' in a measure that is constrained to lie between -1 and 1, but
#' it will not equal gametic D' under HWE.
#'
#' Using \code{type = "geno"} is its own thing and will not equal
#' D' generally under HWE. When \code{type = "geno"}, then the
#' the \code{constrain} parameter has no effect.
#'
#' @param qmat The observed joint genotype distribution.
#' @param type Should we condition on the marginal genotype distribution
#'     (\code{type = "geno"}), or should we condition on the allele frequency
#'     (\code{type = "allele"})?
#' @param constrain A logical. This option is only applicable when
#'     \code{type = "allele"}. Should return an value that is equal
#'     to D' under HWE (\code{FALSE}) or a value that is constrained
#'     to lie between -1 and 1 (\code{TRUE})? Defaults to \code{FALSE}.
#'
#' @return A vector of length 2. The first element is the estimated
#'     D'. The second element is the normalization used.
#'
#' @examples
#' K <- 6
#' qmat <- matrix(stats::runif((K+1)^2), nrow = K+1)
#' qmat <- qmat / sum(qmat)
#' Dprime(qmat, type = "geno")
#' Dprime(qmat, type = "allele")
#'
#' @author David Gerard
#'
#' @export
Dprime <- function(qmat, type = c("geno", "allele"), constrain = FALSE) {
  TOL <- sqrt(.Machine$double.eps)
  stopifnot(abs(sum(qmat) - 1) < TOL)
  stopifnot(qmat > -TOL)
  stopifnot(is.matrix(qmat))
  stopifnot(ncol(qmat) == nrow(qmat))
  stopifnot(is.logical(constrain))
  stopifnot(length(constrain) == 1)
  type <- match.arg(type)

  K <- ncol(qmat) - 1
  D <- Dfromg(gmat = qmat)
  dir <- ifelse(D < 0, "min", "max")

  if (type == "geno") {
    pA <- rowSums(qmat)
    pB <- colSums(qmat)
    qmax <- maxq_from_marg(pA = pA, pB = pB, dir = dir)
    Dmax <- Dfromg(gmat = qmax)
  } else if (type == "allele") {
    pA <- rowSums(qmat)
    pB <- colSums(qmat)
    egA <- sum((0:K) * pA)
    egB <- sum((0:K) * pB)
    if (dir == "min") {
      Dmax <- min(egA * egB, (K - egA) * (K - egB)) / K^2
    } else {
      Dmax <- min(egA * (K - egB), (K - egA) * egB) / K^2
    }

    if (constrain) {
      Dmax <- Dmax * K
    }
  } else { ## Dmax should be K times that in "allele" when constrain = FALSE
    pA <- sum(rowSums(qmat) * 0:K) / K
    pB <- sum(colSums(qmat) * 0:K) / K
    qmax <- maxq_from_allelef(pA = pA, pB = pB, K = K, dir = dir)
    Dmax <- Dfromg(gmat = qmax)
  }

  return(c(Dprime = D / abs(Dmax), Dmax = Dmax))
}
