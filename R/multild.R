

#' Estimate all pair-wise LD's in a collection of SNPs using genotypes or
#' genotype likelihoods.
#'
#' This function is a wrapper to run \code{\link{ldest}()} for many pairs of
#' SNPs. Support is provided for parallelization through the
#' foreach and doParallel packages.
#'
#' See \code{\link{ldest}()} for details on the different types of LD
#' estimators supported.
#'
#' @inheritParams ldest
#' @param nc The number of computing cores to use. This should never be
#'     more than the number of cores available in your computing environment.
#'     You can determine the maximum number of available cores by running
#'     \code{parallel::detectCores()} in R. This is probably fine for a
#'     personal computer, but some environments are only
#'     able to use fewer. Ask your admins if you are unsure.
#' @param geno One of two possible inputs:
#'     \itemize{
#'       \item{A matrix of genotypes (allele dosages). The rows index the
#'             SNPs and the columns index the individuals. That is,
#'             \code{genomat[i, j]} is the allele dosage for individual
#'             \code{j} in SNP \code{i}. When \code{type = "comp"}, the
#'             dosages are allowed to be continuous (e.g. posterior
#'             mean genotypes).}
#'       \item{A three-way array of genotype \emph{log}-likelihoods.
#'             The first dimension indexes the SNPs, the second dimension
#'             indexes the individuals, and the third dimension indexes
#'             the genotypes. That is, \code{genolike_array[i, j, k]}
#'             is the genotype log-likelihood at SNP \code{i} for
#'             individual \code{j} and dosage \code{k - 1}.}
#'     }
#'
#' @return A data frame of class \code{c("lddf", "data.frame")}
#'     with some or all of the following elements:
#' \describe{
#'   \item{\code{i}}{The index of the first SNP.}
#'   \item{\code{j}}{The index of the second SNP.}
#'   \item{\code{snpi}}{The row name corresponding to SNP \code{i}, if
#'       row names are provided.}
#'   \item{\code{snpj}}{The row name corresponding to SNP \code{j}, if
#'       row names are provided.}
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
#' @seealso
#' \describe{
#'   \item{\code{\link{ldest}()}}{For the base function that estimates
#'       pairwise LD.}
#'   \item{\code{\link{sldest}()}}{For estimating pairwise LD along a
#'       sliding window.}
#'   \item{\code{\link{format_lddf}()}}{For formatting the output of
#'       \code{mldest()} as a matrix.}
#'   \item{\code{\link{plot.lddf}()}}{For plotting the output of
#'       \code{mldest()}.}
#' }
#'
#' @examples
#' set.seed(1)
#'
#' ## Simulate genotypes when true correlation is 0
#' nloci <- 5
#' nind  <- 100
#' K <- 6
#' nc <- 1
#' genomat <- matrix(sample(0:K, nind * nloci, TRUE), nrow = nloci)
#'
#' ## Composite LD estimates
#' lddf <- mldest(geno = genomat,
#'                K = K,
#'                nc = nc,
#'                type = "comp")
#' lddf[1:6, 1:6]
#'
#'
#' @author David Gerard
#'
#' @export
mldest <- function(geno,
                   K,
                   nc = 1,
                   type = c("hap", "comp"),
                   model = c("norm", "flex"),
                   pen = ifelse(type == "hap", 2, 1),
                   se = TRUE) {
  model <- match.arg(model)
  type <- match.arg(type)
  if (length(dim(geno)) == 2) {
    outdf <- mldest_geno(genomat = geno,
                         K = K,
                         nc = nc,
                         pen = pen,
                         type = type,
                         se = se)
  } else if (length(dim(geno)) == 3) {
    outdf <- mldest_genolike(genoarray = geno,
                             nc = nc,
                             pen = pen,
                             type = type,
                             model = model,
                             se = se)
  } else {
    stop("mldest: geno needs to either be a matrix or a three-way array.")
  }
  return(outdf)
}

#' Estimate all pair-wise LD's in a collection of SNPs using the genotypes.
#'
#' This function will run \code{\link{ldest}()} iteratively over
#' all possible pairs of SNPs provided. Support is provided for parallelization
#' through the doParallel and foreach packages.
#'
#' @inheritParams mldest
#' @param genomat A matrix of genotypes (allele dosages). The rows index the
#'     SNPs and the columns index the individuals. That is, \code{genomat[i, j]}
#'     is the allele dosage for individual \code{j} in SNP \code{i}.
#' @param K The ploidy of the species. Assumed to be the same for all
#'     individuals at all SNPs
#' @param pen The penalty to be applied to the likelihood. You can think about
#'     this as the prior sample size.
#' @param win The window size. Pairwise LD will be estimated plus or minus
#'     these many positions. Larger sizes significantly increase the
#'     computational load.
#'
#' @author David Gerard
#'
#' @inherit mldest return
#'
#' @examples
#' set.seed(1)
#'
#' ## Simulate genotypes when true correlation is 0
#' nloci <- 5
#' nind  <- 100
#' K <- 6
#' nc <- 1
#' genomat <- matrix(sample(0:K, nind * nloci, TRUE), nrow = nloci)
#'
#' ## Haplotypic LD estimates
#' rdf_hap <- mldest_geno(genomat = genomat,
#'                        K = K,
#'                        nc = nc,
#'                        type = "hap")
#'
#' ## Composite LD estimates
#' rdf_comp <- mldest_geno(genomat = genomat,
#'                         K = K,
#'                         nc = nc,
#'                         type = "comp")
#'
#' ## Haplotypic and Composite LD are both close to 0
#' ## (because HWE is fulfilled)
#' ## But composite is more variable.
#' Dmat <- cbind(Haplotypic = rdf_hap$D, Composite = rdf_comp$D)
#' boxplot(Dmat,
#'         main = "Estimates of D",
#'         ylab = "D Estimate",
#'         xlab = "Type")
#' abline(h = 0, lty = 2, col = 2)
#'
#'
#' @noRd
mldest_geno <- function(genomat,
                        K,
                        nc = 1,
                        type = c("hap", "comp"),
                        pen = ifelse(type == "hap", 2, 1),
                        se = TRUE,
                        win = NULL) {
  type <- match.arg(type)
  stopifnot(is.matrix(genomat))
  stopifnot(is.logical(se))
  nloci <- nrow(genomat)
  if (is.null(win)) {
    win <- nloci
  } else {
    stopifnot(win > 0)
    stopifnot(length(win) == 1)
  }

  ## Register workers ---------------------------------------------------------
  if (nc == 1) {
    foreach::registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(nc)
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
                             .export = c("ldest",
                                         "nullvec_hap",
                                         "nullvec_comp")) %dopar% {

                               if (type == "hap") {
                                 ldnull <- nullvec_hap()
                               } else {
                                 ldnull <- nullvec_comp(K = K, model = "flex") ## stupid hack. "flex" otherwise you get slots for sigma and mu
                               }
                               endit <- min(nloci, i + win)

                               estmat <- matrix(NA_real_,
                                                nrow = endit - i,
                                                ncol = length(ldnull) + 4)
                               colnames(estmat) <- c("i", "j", "snpi", "snpj", names(ldnull))


                               for (j in (i + 1):endit) {
                                 estmat[j - i, 1] <- i
                                 estmat[j - i, 2] <- j
                                 tryCatch({
                                   ldout <- ldest(ga = genomat[i, ],
                                                  gb = genomat[j, ],
                                                  K = K,
                                                  type = type,
                                                  pen = pen,
                                                  se = se)
                                   estmat[j - i, -(1:4)] <- ldout
                                 }, error = function(e) NULL)

                               }
                               estmat
                             }

  if (nc > 1) {
    parallel::stopCluster(cl)
  }

  outmat <- as.data.frame(outmat)
  class(outmat) <- c("lddf", "data.frame")

  ## Check for snp names ------------------------------------------------------
  if (!is.null(rownames(genomat))) {
    snpnamevec <- rownames(genomat)
    outmat$snpi <- snpnamevec[outmat$i]
    outmat$snpj <- snpnamevec[outmat$j]
  }

  return(outmat)
}

#' Estimate all pair-wise LD's in a collection of SNPs using the genotype
#' likelihoods.
#'
#' This function will run \code{\link{ldest}()} iteratively over
#' all possible pairs of SNPs provided. Support is provided for parallelization
#' through the doParallel and foreach packages.
#'
#' @inheritParams mldest
#' @param genoarray A three-way array of genotype \emph{log}-likelihoods.
#'     The first dimension indexes the SNPs, the second dimension indexes
#'     the individuals, and the third dimension indexes the genotypes.
#'     That is, \code{genolike_array[i, j, k]} is the genotype log-likelihood
#'     at SNP \code{i} for individual \code{j} and dosage \code{k - 1}.
#'     The ploidy (assumed to be the same for all individuals) is assumed to
#'     be one minus the size of the third dimension.
#' @param pen The penalty to be applied to the likelihood. You can think about
#'     this as the prior sample size.
#' @param win The window size. Pairwise LD will be estimated plus or minus
#'     these many positions. Larger sizes significantly increase the
#'     computational load.
#'
#' @author David Gerard
#'
#' @inherit mldest return
#'
#' @examples
#' set.seed(1)
#'
#' ## Simulate some data with true correlation of 0
#' nloci <- 5
#' nind  <- 100
#' K <- 6
#' nc <- 1
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
#' ## Haplotypic LD estimates
#' rdf_hap <- mldest_genolike(genoarray = genoarray,
#'                            nc = nc,
#'                            type = "hap")
#'
#' ## Composite LD estimates
#' rdf_comp <- mldest_genolike(genoarray = genoarray,
#'                             nc = nc,
#'                             type = "comp")
#'
#' ## Haplotypic and Composite LD are both close to 0.
#' ## But composite is more variable.
#' Dmat <- cbind(Haplotypic = rdf_hap$D, Composite = rdf_comp$D)
#' boxplot(Dmat,
#'         main = "Estimates of D",
#'         ylab = "D Estimate",
#'         xlab = "Type")
#' abline(h = 0, lty = 2, col = 2)
#'
#' @noRd
mldest_genolike <- function(genoarray,
                            nc = 1,
                            type = c("hap", "comp"),
                            model = c("norm", "flex"),
                            pen = ifelse(type == "hap", 2, 1),
                            se = TRUE,
                            win = NULL) {
  type <- match.arg(type)
  model <- match.arg(model)
  stopifnot(is.logical(se))
  stopifnot(is.array(genoarray))
  stopifnot(length(dim(genoarray)) == 3)
  nloci <- dim(genoarray)[[1]]
  nind <- dim(genoarray)[[2]]
  K <- dim(genoarray)[[3]] - 1
  if (is.null(win)) {
    win <- nloci
  } else {
    stopifnot(win > 0)
    stopifnot(length(win) == 1)
  }

  ## Register workers ---------------------------------------------------------
  if (nc == 1) {
    foreach::registerDoSEQ()
  } else {
    cl <- parallel::makeCluster(nc)
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
                             .export = c("ldest",
                                         "nullvec_hap",
                                         "nullvec_comp")) %dopar% {

                               if (type == "hap") {
                                 ldnull <- nullvec_hap()
                               } else {
                                 ldnull <- nullvec_comp(K = K, model = model)
                               }

                               endit <- min(nloci, i + win)
                               estmat <- matrix(NA_real_,
                                                nrow = endit - i,
                                                ncol = length(ldnull) + 4)
                               colnames(estmat) <- c("i", "j", "snpi", "snpj", names(ldnull))

                               for (j in (i + 1):endit) {
                                 estmat[j - i, 1] <- i
                                 estmat[j - i, 2] <- j
                                 tryCatch({
                                   ldout <- ldest(ga = genoarray[i, , ],
                                                  gb = genoarray[j, , ],
                                                  K = K,
                                                  type = type,
                                                  model = model,
                                                  pen = pen,
                                                  se = se)
                                   estmat[j - i, -(1:4)] <- ldout
                                 }, error = function(e) NULL)
                               }
                               estmat
                             }

  if (nc > 1) {
    parallel::stopCluster(cl)
  }

  outmat <- as.data.frame(outmat)
  class(outmat) <- c("lddf", "data.frame")

  ## Check for snp names ------------------------------------------------------
  if (!is.null(dimnames(genoarray))) {
    snpnamevec <- dimnames(genoarray)[[1]]
    outmat$snpi <- snpnamevec[outmat$i]
    outmat$snpj <- snpnamevec[outmat$j]
  }

  return(outmat)
}
