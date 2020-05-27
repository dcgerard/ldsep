###################
## Functions to estimate LD directly from genotypes
###################



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
    inity_list <- list(c(0, 0, 0))
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
                           pgA      = ga,
                           pgB      = gb,
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

#' Estimate pair-wise LD using either genotypes or genotype likelihoods.
#'
#' Given genotype (allele dosage) or genotype likelihood data
#' for each individual at a pair of loci, this function will
#' calculate the maximum likelihood estimates
#' and their corresponding asymptotic standard errors of some
#' measures of linkage disequilibrium (LD): D, the Pearson correlation,
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
#' The resulting standard errors are the square roots of the inverse of the
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
#' will return LD estimates of \code{0}.
#'
#' @param ga One of two possible inputs: (i) A vector of counts, containing the genotypes for each individual
#'     at the first locus; or (ii) A matrix of genotype log-likelihoods
#'     at the second locus. The rows index the individuals and the columns
#'     index the genotypes. That is \code{ga[i, j]} is the genotype
#'     likelihood of individual \code{i} for genotype \code{j}.
#' @param gb One of two possible inputs: (i) A vector of counts, containing the genotypes for each individual
#'     at the second locus; or (ii) A matrix of genotype log-likelihoods at
#'     the second locus. The rows index the individuals and the columns
#'     index the genotypes. That is \code{ga[i, j]} is the genotype
#'     likelihood of individual \code{i} for genotype \code{j}.
#' @param K The ploidy of the species. Assumed the same for all individuals.
#' @param reltol The relative tolerance for the stopping criterion.
#' @param nboot Sometimes, the MLE standard errors don't exist. So we use
#'     the bootstrap as a backup. \code{nboot} specifies the number
#'     of bootstrap iterations.
#' @param useboot A logical. Optionally, you may always use the bootstrap
#'     to estimate the standard errors (\code{TRUE}). These will be more
#'     accurate but also much slower, so this defaults to \code{FALSE}.
#' @param pen The penalty to be applied to the likelihood. You can think about
#'     this as the prior sample size.
#' @param grid_init A logical. Should we initialize the gradient ascent
#'     at a grid of initial values (\code{TRUE}) or just initialize
#'     at one value corresponding to the simplex point
#'     \code{rep(0.25, 4)} (\code{FALSE})?
#'
#' @author David Gerard
#'
#' @return A vector with some or all of the following elements:
#' \describe{
#'   \item{\code{D}}{The MLE of D.}
#'   \item{\code{D_se}}{The standard error of the estimate of D.}
#'   \item{\code{r2}}{The MLE of the squared Pearson correlation.}
#'   \item{\code{r2_se}}{The standard error of the estimate of the
#'       squared Pearson correlation.}
#'   \item{\code{r}}{The MLE of the Pearson correlation.}
#'   \item{\code{r_se}}{The standard error of the estimate of the
#'       Pearson correlation.}
#'   \item{\code{z}}{The Fisher-z transformation of \code{r}.}
#'   \item{\code{z_se}}{The standard error of the Fisher-z
#'       transformation of \code{r}.}
#'   \item{\code{p_ab}}{The estimated haplotype frequency of ab.}
#'   \item{\code{p_Ab}}{The estimated haplotype frequency of Ab.}
#'   \item{\code{p_aB}}{The estimated haplotype frequency of aB.}
#'   \item{\code{p_AB}}{The estimated haplotype frequency of AB.}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 100
#' K <- 6
#'
#' ## If you give ldest() vectors,
#' ## it assumes you are using genotypes.
#' ga <- stats::rbinom(n = n, size = K, prob = 0.5)
#' gb <- stats::rbinom(n = n, size = K, prob = 0.5)
#' head(ga)
#' head(gb)
#' ldout1 <- ldest(ga = ga,
#'                 gb = gb,
#'                 K = K,
#'                 grid_init = FALSE)
#' ldout1
#'
#' ## Use bootstap for standard errors instead:
#' ldout2 <- ldest(ga = ga,
#'                 gb = gb,
#'                 K = K,
#'                 useboot = TRUE,
#'                 grid_init = FALSE)
#' ldout2
#'
#' ## Standard error estimates are similar for D, r, and z
#' ldout1[c("D_se", "r2_se", "r_se", "z_se")]
#' ldout2[c("D_se", "r2_se", "r_se", "z_se")]
#'
#' ## If you give ldest() matrices,
#' ## it assumes you are using genotype likelihoods.
#' gamat <- t(sapply(ga, stats::dnorm, x = 0:K, sd = 1, log = TRUE))
#' gbmat <- t(sapply(gb, stats::dnorm, x = 0:K, sd = 1, log = TRUE))
#' head(gamat)
#' head(gbmat)
#' ldout3 <- ldest(ga = gamat, gb = gbmat, K = K)
#' ldout3
#'
#' @export
ldest <- function(ga,
                  gb,
                  K,
                  reltol  = 10^-8,
                  nboot   = 100,
                  useboot = FALSE,
                  pen     = 2,
                  grid_init = FALSE) {

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
  if (any(abs(eval$values) < sqrt(.Machine$double.eps)) || useboot) {
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
  } else {
    ## MLE theory

    ## jacobian converting from simplex to real. The last is just a column of
    ## zeros because we use the transform from p1, p2, p3 to y1, y2, y3
    J <- dsimplex_to_real_dx(phat)[, -4]
    Hp <- t(J) %*% Hy %*% J ## Hessian of first three simplex parameters
    nHp_inv <- -solve(Hp)

    dD <- dD_dprob(prob = phat)
    D_se <- sqrt(t(dD) %*% nHp_inv %*% dD)

    # dDprime <- dDprime_dprob(prob = phat)
    # Dprime_se <- sqrt(t(dDprime) %*% nHp_inv %*% dDprime)

    dr2 <- dr2_dprob(prob = phat)
    r2_se <- sqrt(t(dr2) %*% nHp_inv %*% dr2)

    r_se <- r2_se / sqrt(4 * r2)

    z_se <- r_se / (1 - r2)

    # g <- log(-log(r2))
    # g_se <- r2_se / abs(r2 * log(r2))
  }

  retvec <- c(D         = D,
              D_se      = D_se,
              r2        = r2,
              r2_se     = r2_se,
              r         = r,
              r_se      = r_se,
              z         = z,
              z_se      = z_se,
              p_ab      = phat[[1]],
              p_Ab      = phat[[2]],
              p_aB      = phat[[3]],
              p_AB      = phat[[4]])

  return(retvec)
}


#' Estimate all pair-wise LD's in a collection of SNPs using the genotypes.
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
#' @param nc The number of computing cores to use. This should never be
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
#' ## Simulate genotypes when true correlation is 0
#' nloci <- 5
#' nind  <- 100
#' K <- 6
#' nc <- 1
#' genomat <- matrix(sample(0:K, nind * nloci, TRUE), nrow = nloci)
#'
#' ## Estimate LD
#' rdf <- mldest_geno(genomat = genomat, K = K, nc = nc)
#'
#' @export
mldest_geno <- function(genomat, K, nc = 1, reltol = 10^-8) {
  stopifnot(is.matrix(genomat))
  nloci <- nrow(genomat)

  ## Register workers ----------------------------------------------------------
  if (nc == 1) {
    foreach::registerDoSEQ()
  } else {
    cl = parallel::makeCluster(nc)
    doParallel::registerDoParallel(cl = cl)
    if (foreach::getDoParWorkers() == 1) {
      stop(paste0("mldest_geno: nc > 1 ",
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
                                                reltol = reltol)
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


#' Estimate all pair-wise LD's in a collection of SNPs using the genotype
#' likleihoods.
#'
#' This function will run \code{\link{ldest}()} iteratively over
#' all possible pairs of SNPs provided. Support is provided for parallelization
#' through the doParallel and foreach packages.
#'
#' @param genoarray An three-way array of genotype \emph{log}-likelihoods.
#'     The first dimension indexes the SNPs, the second dimension indexes
#'     the individuals, and the third dimension indexes the genotypes.
#'     That is, \code{genolike_array[i, j, k]} is the genotype log-likelihood
#'     at SNP \code{i} for individual \code{j} and dosage \code{k - 1}.
#'     The ploidy (assumed to be the same for all individuals) is assumed to
#'     be one minus the size of the third dimension.
#' @param nc The number of computing cores to use. This should never be
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
#' ## Simulate some data with true correlation of 0
#' nloci <- 10
#' nind  <- 100
#' K <- 6
#' genovec <- sample(0:K, nind * nloci, TRUE)
#' genomat <- matrix(genovec, nrow = nloci)
#' genoarray <- array(sapply(genovec,
#'                           stats::dnorm,
#'                           x = 0:K,
#'                           sd = 1,
#'                           log = TRUE),
#'                    dim = c(K + 1, nloci, nind))
#' genoarray <- aperm(genoarray, c(2, 3, 1)) ## loci, individuals, dosages
#'
#' ## Verify simulation
#' locnum <- sample(seq_len(nloci), 1)
#' indnum <- sample(seq_len(nind), 1)
#' genoarray[locnum, indnum, ]
#' stats::dnorm(x = 0:K, mean = genomat[locnum, indnum], sd = 1, log = TRUE)
#'
#' ## Find pairwise LD between all loci
#' nc <- 1
#' rdf <- mldest_genolike(genoarray = genoarray, nc = nc)
#' rdf
#'
#' @export
mldest_genolike <- function(genoarray, nc = 1, reltol = 10^-8) {
  stopifnot(is.array(genoarray))
  stopifnot(length(dim(genoarray)) == 3)
  nloci <- dim(genoarray)[[1]]
  nind <- dim(genoarray)[[2]]
  K <- dim(genoarray)[[3]] - 1

  ## Register workers ----------------------------------------------------------
  if (nc == 1) {
    foreach::registerDoSEQ()
  } else {
    cl = parallel::makeCluster(nc)
    doParallel::registerDoParallel(cl = cl)
    if (foreach::getDoParWorkers() == 1) {
      stop(paste0("mldest_geno: nc > 1 ",
                  "but only one core registered from ",
                  "foreach::getDoParWorkers()."))
    }
  }

  i <- 1
  outmat <- foreach::foreach(i = seq_len(nloci - 1),
                             .combine = rbind,
                             .export = c("ldest")) %dopar% {

                               for (j in (i + 1):nloci) {
                                 ldout <- ldest(ga = genoarray[i, , ],
                                                gb = genoarray[j, , ],
                                                K = K,
                                                reltol = reltol)
                                 if (j == i + 1) {
                                   estmat <- matrix(NA_real_,
                                                    nrow = nloci - i,
                                                    ncol = length(ldout) + 2)
                                   colnames(estmat) <- c("i", "j", names(ldout))
                                 }
                                 estmat[j - i, 1] <- i
                                 estmat[j - i, 2] <- j
                                 estmat[j - i, -c(1, 2)] <- ldout
                               }
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


#' Format an element of \code{\link{mldest_geno}()} or
#' \code{\link{mldest_genolike}()} into an
#' upper-triangular matrix.
#'
#' Formats the correlation estimates and standard errors output
#' from running \code{\link{mldest_geno}()} or \code{\link{mldest_genolike}()}
#' into a more conventional upper-triangular matrix.
#'
#' @param obj An object of class \code{lddf}, usually output from
#'     running either \code{\link{mldest_geno}()} or
#'     \code{\link{mldest_genolike}()}.
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


