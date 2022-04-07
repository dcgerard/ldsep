#' A function to simulate genotype log-likelihoods under a given correlation under HWE for two loci
#'
#' A user inputs the correlation (\code{rho}), the number of individuals
#' (\code{nind}), the ploidy (\code{ploidy}), the allele frequencies at two
#' loci (\code{pA} and \code{pB}), and the "read-depth" (larger corresponds
#' to being more sure of the genotypes) (\code{rdepth}). Then this
#' function will  simulate some genotype log-likelihoods.
#'
#' @param rho The true correlation.
#' @param nind The number of individuals to simulate.
#' @param rdepth The read depth to simulate (larger means more certain).
#' @param ploidy The ploidy of the species.
#' @param pA The allele frequency of the first locus.
#' @param pB The allele frequency of the second locus.
#' @param seq The sequencing error rate.
#' @param bias The allele bias.
#' @param od The overdispersion parameter.
#'
#' @return A list with some or all of the following elements
#' \describe{
#'   \item{lA}{The genotype log-likelihoods at the first locus. \code{lA[i, k]}
#'       is the genotype log-likelihood for individual \code{i} at genotype \code{k - 1}.}
#'   \item{lB}{The genotype log-likelihoods at the second locus. \code{lB[i, K]}
#'       is the genotype log-likelihood for individual \code{i} at genotype \code{k - 1}.}
#'   \item{rho}{The true Pearson correlation that the user provided.}
#'   \item{delta}{The corresponding LD coefficient.}
#' }
#'
#' @examples
#' glsim_pairwise(rho = 0.1, nind = 100, ploidy = 6, pA = 0.5, pB = 0.5)
#'
#' @export
#'
#' @author David Gerard
glsim_pairwise <- function(rho, nind, ploidy, pA, pB, rdepth = 10,
                           seq = 0.001, bias = 1, od = 0.001) {
  stopifnot(rho >= 0, rho <= 1)
  stopifnot(nind > 0)
  stopifnot(rdepth > 0)
  stopifnot(ploidy > 0, ploidy %% 2 == 0)
  stopifnot(pA > 0, pA < 1)
  stopifnot(pB > 0, pB < 1)
  stopifnot(length(rho) == 1,
            length(nind) == 1,
            length(rdepth) == 1,
            length(ploidy) == 1,
            length(pA) == 1,
            length(pB) == 1)

  ## Get haplotype probabilities ----
  pAB = rho * sqrt(pA * (1 - pA) * pB * (1 - pB)) + pA * pB
  pAb = pA - pAB
  paB = pB - pAB
  pab = 1 - pAB - pAb - paB
  delta <- pAB - pA * pB

  stopifnot(pAB >= 0, pAb >= 0, paB >= 0, pab >= 0,
            pAB <= 1, pAb <= 1, paB <= 1, pab <= 1)

  ## Simulate haplotypes and genotypes ----
  hapval <- stats::rmultinom(n = nind,
                             size = ploidy,
                             prob = c(pAB, pAb, paB, pab))

  gA <- colSums(hapval[1:2, ])
  gB <- colSums(hapval[c(1, 3), ])
  sizevec <- rep(rdepth, nind)

  ## Simulate read-depths given genotypes ----
  refA <- updog::rflexdog(sizevec = sizevec,
                          geno    = gA,
                          ploidy  = ploidy,
                          seq     = seq,
                          bias    = bias,
                          od      = od)
  refB <- updog::rflexdog(sizevec = sizevec,
                          geno    = gB,
                          ploidy  = ploidy,
                          seq     = seq,
                          bias    = bias,
                          od      = od)

  ## Get genotype likelihoods give read-depths ----
  foutA <- updog::flexdog(refvec      = refA,
                          sizevec     = sizevec,
                          ploidy      = ploidy,
                          verbose     = FALSE,
                          bias_init   = 1,
                          update_bias = FALSE,
                          seq         = 0.01,
                          update_seq  = FALSE,
                          od          = 0.01,
                          update_od   = FALSE,
                          model       = "hw")
  foutB <- updog::flexdog(refvec      = refB,
                          sizevec     = sizevec,
                          ploidy      = ploidy,
                          verbose     = FALSE,
                          bias_init   = 1,
                          update_bias = FALSE,
                          seq         = 0.01,
                          update_seq  = FALSE,
                          od          = 0.01,
                          update_od   = FALSE,
                          model       = "hw")

  ## Extract genotype likelihoods from that output ----
  retlist <- list(lA = foutA$genologlike,
                  lB = foutB$genologlike,
                  rho = rho,
                  delta = delta)
  return(retlist)
}

#' Return the true correlation from \code{\link{hap_sim}()}.
#'
#' The correlation is \code{rho^a}, where \code{a} is the number of
#' positions away two SNPs are.
#'
#' @inheritParams hap_sim
#'
#' @return A matrix of size \code{nloc} by \code{nloc} that contains the
#'    true correlation of the data simulated in \code{\link{hap_sim}()}.
#'
#' @author David Gerard and Hanwei Hu
#'
#' @examples
#' hap_truecor(0.9, 10)
#'
#' @export
hap_truecor <- function(rho, nloc) {
  stopifnot(nloc >= 0, length(nloc) == 1)
  stopifnot(rho >= 0, rho <= 1, length(rho) == 1)
  return(rho^abs(outer(seq_len(nloc), seq_len(nloc), "-")))
}

#' Simulate haplotypes with a banded correlation structure.
#'
#' Simulates haplotypes at biallelic loci where each adjacent pair has
#' correlation \code{rho} and each locus has allele frequency \code{alpha}.
#'
#' @param rho The correlation between adjacent loci. If loci are \code{a}
#'     positions apart, then the correlation between them is \code{rho^a}.
#' @param alpha The allele frequency at each locus.
#' @param nloc The number of loci to simulate.
#'
#' @return A vector of length \code{nloc} with 0's and 1's, with a
#'     provided correlation structure.
#'
#' @author David Gerard and Hanwei Hu
#'
#' @examples
#' set.seed(1)
#' hap_sim(0.5, 0.5, 50)
#'
#' @export
hap_sim <- function(rho, alpha, nloc) {
  stopifnot(nloc >= 0, length(nloc) == 1)
  stopifnot(rho >= 0, rho <= 1, length(rho) == 1)
  stopifnot(alpha >= 0, alpha <= 1, length(alpha) == 1)

  hapvec <- rep(NA_real_, length.out = nloc)

  if (nloc == 0) {
    return(double())
  }

  p11 <- rho * (1 - alpha) + alpha # Pr(next = 1 | last = 1)
  p10 <- alpha * (1 - rho) # Pr(next = 1 | last = 0)

  if (p11 < 0 | p11 > 1 | p10 < 0 | p10 > 1) {
    stop("alpha row combo not allowed")
  }

  # Sample initial value
  hapvec[[1]] <- sample(x = c(0, 1), size = 1, prob = c(1 - alpha, alpha))

  for (ell in 2:nloc) {
    if (hapvec[[ell - 1]] == 1) {
      hapvec[[ell]] <- sample(c(0, 1), size = 1, prob = c(1 - p11, p11))
    } else {
      hapvec[[ell]] <- sample(c(0, 1), size = 1, prob = c(1 - p10, p10))
    }
  }

  return(hapvec)
}

#' Simulate genotypes with a banded correlation structure.
#'
#' Simulates genotypes at biallelic loci where each adjacent pair has
#' correlation \code{rho} and each locus has allele frequency \code{alpha}.
#'
#' @inheritParams hap_sim
#' @param ploidy The ploidy of the individual.
#'
#' @return A vector of length \code{nloc} with the dosage of the individual,
#'     with a provided correlation structure.
#'
#' @author David Gerard and Hanwei Hu
#'
#' @examples
#' set.seed(1)
#' hap_sim(0.5, 0.5, 50)
#'
#' @export
geno_sim <- function(rho, alpha, nloc, ploidy){
  genovec <- vector(mode = "numeric", length = nloc)
  for(j in 1 : ploidy) {
    genovec <- genovec + hap_sim(rho = rho, alpha = alpha, nloc = nloc)
  }
  return(genovec)
}

#' Simulate genotype likelihood array for individuals given a banded correlation structure.
#'
#' Simulates genotypes using \code{\link{geno_sim}()}, then simulates
#' sequencing data using the model of Gerard et al (2018), then fits
#' \code{\link[updog]{flexdog}()} on these data to obtain genotype likelihoods.
#'
#' @inheritParams geno_sim
#' @param nind The number of individuals to sample.
#' @param rdepth The read-depth of the sequencing reads used to create the
#'     genoype likelihoods. Lower means greater uncertainty.
#' @param seq The sequencing error rate.
#' @param bias The allele bias.
#' @param od The overdispersion rate.
#'
#' @references
#' \itemize{
#'   \item{Gerard, D., Ferrão, L. F. V., Garcia, A. A. F., & Stephens, M. (2018). Genotyping polyploids from messy sequencing data. Genetics, 210(3), 789-807. \doi{genetics.118.301468}}
#' }
#'
#' @author David Gerard and Hanwei Hu
#'
#' @examples
#' glsim(rho = 0.9, alpha = 0.5, nloc = 10, ploidy = 4, nind = 25)
#'
#' @export
glsim <- function(rho, alpha, nloc, ploidy, nind, rdepth = 10, seq = 0.001, bias = 1, od = 0.001){

  sizemat <- matrix(rdepth, ncol = nloc, nrow = nind)

  genomat <- matrix(NA_real_, ncol = nloc, nrow = nind)
  for (i in seq_len(nind)) {
    genomat[i, ] <- geno_sim(rho = rho, alpha = alpha, nloc = nloc, ploidy = ploidy)
  }

  refmat <- matrix(NA_real_, ncol = nloc, nrow = nind)
  for (j in seq_len(nloc)) {
    refmat[, j]  <- updog::rflexdog(sizevec = sizemat[,j],
                                    geno = genomat[, j],
                                    ploidy = ploidy,
                                    seq = seq,
                                    bias = bias,
                                    od = od)
  }

  rownames(refmat) <- 1:nrow(refmat)
  rownames(sizemat) <- 1:nrow(sizemat)
  colnames(refmat) <- 1:ncol(refmat)
  colnames(sizemat) <- 1:ncol(sizemat)

  trash <- utils::capture.output(
    uout <- updog::multidog(refmat = t(refmat),
                            sizemat = t(sizemat),
                            ploidy = ploidy,
                            model = "norm",
                            bias_init = bias,
                            seq = seq,
                            od = od,
                            update_bias = FALSE,
                            update_seq = FALSE,
                            update_od = FALSE)
  )

  llarray <- updog::format_multidog(x = uout, varname = paste0("logL_", 0:ploidy))

  return(llarray)
}
