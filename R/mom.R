########################
## Method of moments estimators
########################

#' Simple LD Estimators
#'
#' Provides quick and simple estimates of LD and their corresponding
#' standard errors. These estimates are just based on the sample
#' moments of the genotypes or the posterior mean genotypes.
#'
#' For large sequencing depth and large n, these estimates perform as well
#' as the MLE's but are much faster to calculate. For small n and, in
#' particular, small sequencing depth, these estimates tend to be very
#' biased.
#'
#' @param ga Either the genotype at locus 1 or the posterior mean genotype at
#'     locus 1.
#' @param gb Either the genotype at locus 2 or the posterior mean genotype at
#'     locus 2.
#' @param K The ploidy of hte species. Assumed to be the same for all
#'     individuals.
#'
#' @inherit ldest return
#'
#' @author David Gerard
#'
#' @noRd
#'
ldsimp <- function(ga, gb, K) {
  n     <- length(ga)
  pA    <- mean(ga) / K
  pB    <- mean(gb) / K
  D     <- stats::cov(ga, gb) / K
  vara  <- stats::var(ga)
  varb  <- stats::var(gb)
  D_se  <- sqrt(vara * varb / K ^ 2 + D ^ 2) / sqrt(n)
  if (D < 0) {
    Dmax <- min(pA * pB, (1 - pA) * (1 - pB))
  } else {
    Dmax <- min(pA * (1 - pB), pA * (1 - pB))
  }
  Dprime <- D / Dmax
  Dprime_se <- D_se / Dmax
  r     <- stats::cor(ga, gb)
  r_se  <- (1 - r ^ 2) / sqrt(n)
  z     <- atanh(r)
  z_se  <- 1 / sqrt(n - 3)
  r2    <- r ^ 2
  r2_se <- 2 * abs(r) * (1 - r ^ 2) / sqrt(n)

  qvec <- c(
    as.matrix(
      prop.table(
        table(factor(round(ga), levels = 0:K),
              factor(round(gb), levels = 0:K))
      )
    )
  )
  inddf <- expand.grid(i = 0:K, j = 0:K)
  names(qvec) <- paste0("q", inddf$i, inddf$j)

  retvec <- c(D         = D,
              D_se      = D_se,
              r2        = r2,
              r2_se     = r2_se,
              r         = r,
              r_se      = r_se,
              Dprime    = Dprime,
              Dprime_se = Dprime_se,
              z         = z,
              z_se      = z_se)
  retvec <- c(retvec, qvec)

  return(retvec)
}


#' Estimates of composite pairwise LD based either on genotype estimates or
#' genotype likelihoods.
#'
#' This function will estimate the composite LD between two loci, either
#' using genotype estimates or using genotype likelihoods. The resulting
#' measures of LD are generalizations of Burrow's "composite" LD measure.
#'
#' @inheritParams ldest
#'
#' @inherit ldest return
#'
#' @author David Gerard
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
#' ## Composite LD with genotypes
#' ldout1 <- ldest_comp(ga = ga,
#'                      gb = gb,
#'                      K = K)
#' head(ldout1)
#'
#' ## Composite LD with genotype likelihoods
#' ldout2 <- ldest_comp(ga = gamat,
#'                      gb = gbmat,
#'                      K = K)
#' head(ldout2)
#'
#' @export
ldest_comp <- function(ga,
                       gb,
                       K,
                       pen = 2) {
  stopifnot(pen > 1)
  if (is.vector(ga) & is.vector(gb)) {
    stopifnot(length(ga) == length(gb))
    stopifnot(ga >= 0, ga <= K)
    stopifnot(gb >= 0, gb <= K)
    using <- "genotypes"
  } else if (is.matrix(ga) & is.matrix(gb)) {
    stopifnot(dim(ga) == dim(gb))
    stopifnot(K + 1 == ncol(ga))
    using <- "likelihoods"
  } else {
    stop("ldest: ga and gb must either both be vectors or both be matrices.")
  }

  if (using == "genotypes") {
    retvec <- ldsimp(ga = ga, gb = gb, K = K)
  } else {
    alphamat <- matrix(data = pen, nrow = K + 1, ncol = K + 1)
    pma <- factor(apply(X = ga, MARGIN = 1, FUN = which.max) - 1, levels = 0:K)
    pmb <- factor(apply(X = gb, MARGIN = 1, FUN = which.max) - 1, levels = 0:K)
    pinit <- as.matrix(prop.table(table(pma, pmb) + 1))
    gout <- em_jointgeno(p = pinit,
                         pgA = ga,
                         pgB = gb,
                         alpha = alphamat)
    D <- Dfromg(gmat = gout)
    r2 <- r2fromg(gmat = gout)
    r <- sqrt(r2) * sign(D)
    z <- atanh(r)

    distA <- rowSums(gout)
    distB <- colSums(gout)
    egA <- sum((0:K) * distA)
    egB <- sum((0:K) * distB)
    if (D < 0) {
      Dmax <- min(egA * egB, (K - egA) * (K - egB)) / K^2
    } else {
      Dmax <- min(egA * (K - egB), (K - egA) * egB) / K^2
    }
    Dprime <- D / Dmax

    ## get asymptotic covaraince ----
    finfo <- -solve(
      hessian_jointgeno(p = gout, pgA = ga, pgB = gb, alpha = alphamat)
    )

    ## Standard errors ----
    grad_dq <- dD_dqlm(p = gout)
    D_se <- sqrt(c(t(grad_dq) %*% finfo %*% grad_dq))

    grad_r2q <- dr2_dqlm(p = gout, dgrad = grad_dq, D = D)
    r2_se <- sqrt(c(t(grad_r2q) %*% finfo %*% grad_r2q))

    grad_dprimeq <- ddprime_dqlm(p = gout, dgrad = grad_dq, D = D, Dm = Dmax)
    Dprime_se <- sqrt(c(t(grad_dprimeq) %*% finfo %*% grad_dprimeq))

    r_se <- r2_se / sqrt(4 * r2)

    z_se <- r_se / (1 - r2)

    ## set up retvec ----
    retvec <- c(D         = D,
                D_se      = D_se,
                r2        = r2,
                r2_se     = r2_se,
                r         = r,
                r_se      = r_se,
                Dprime    = Dprime,
                Dprime_se = Dprime_se,
                z         = z,
                z_se      = z_se)

    qvec <- c(gout)
    inddf <- expand.grid(i = 0:K, j = 0:K)
    names(qvec) <- paste0("q", inddf$i, inddf$j)
    retvec <- c(retvec, qvec)
  }

  return(retvec)
}


#' Obtain D measure of LD from joint genotype frequencies.
#'
#' @param gmat Element (i, j) is the probability of gentoype i on locus 1
#'     and genotype j on locus 2.
#'
#' @author David Gerard
#'
#' @noRd
Dfromg <- function(gmat) {
  TOL <- sqrt(.Machine$double.eps)
  stopifnot(nrow(gmat) == ncol(gmat))
  stopifnot(abs(sum(gmat) - 1) < TOL)
  K <- ncol(gmat) - 1
  sum(
    sweep(x = sweep(x = gmat,
                    MARGIN = 1,
                    STATS = 0:K,
                    FUN = `*`),
          MARGIN = 2,
          STATS = 0:K,
          FUN = `*`)) / K -
    sum(0:K * rowSums(gmat)) * sum(0:K * colSums(gmat)) / K
}

#' Obtain r^2 measure of LD from joint genotype frequencies.
#'
#' @param gmat Element (i, j) is the probability of gentoype i on locus 1
#'     and genotype j on locus 2.
#'
#' @author David Gerard
#'
#' @noRd
r2fromg <- function(gmat) {
  stopifnot(ncol(gmat) == nrow(gmat))
  K <- ncol(gmat) - 1
  D <- Dfromg(gmat)

  distA <- rowSums(gmat)
  distB <- colSums(gmat)

  egA <- sum((0:K) * distA)
  egB <- sum((0:K) * distB)
  egA2 <- sum((0:K)^2 * distA)
  egB2 <- sum((0:K)^2 * distB)

  vargA <- egA2 - egA^2
  vargB <- egB2 - egB^2

  return(K^2 * D^2 / (vargA * vargB))
}


#' Wrapper for llike_jointgeno so that takes vector of arguments
#'
#' @param par A vector of proportions. \code{matrix(par, K+1, K+1)}
#'     should recover the probability matrix in \code{\link{em_jointgeno}()}.
#' @inheritParams em_jointgeno
#'
#' @author David Gerard
#'
#' @noRd
llike_jointgeno_vec <- function(par, pgA, pgB, alpha) {
  stopifnot(dim(pgA) == dim(pgB))
  stopifnot(length(par) == ncol(pgA)^2)
  stopifnot(dim(alpha) == rep(ncol(pgA), 2))
  pmat <- matrix(par, nrow = ncol(pgA), ncol = ncol(pgB))
  llike_jointgeno(p = pmat, pgA = pgA, pgB = pgB, alpha = alpha)
}
