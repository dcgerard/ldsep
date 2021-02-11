

#' Updog fits on the data from Uitdewilligen et. al. (2013)
#'
#' 10 SNPs from the "PGSC0003DMB000000062" super scaffold were genotyped
#' using the \code{\link[updog]{multidog}()} function from the updog R package.
#' These data are the resulting output.
#'
#' @format An object of class \code{\link[updog]{multidog}()}.
#'     See the documentation from the updog R package.
#'
#' @source \doi{10.1371/journal.pone.0062355}
#'
#' @references
#' \itemize{
#'   \item{Uitdewilligen, Jan GAML, Anne-Marie A. Wolters, B. Bjorn, Theo JA Borm, Richard GF Visser, and Herman J. Van Eck. "A next-generation sequencing method for genotyping-by-sequencing of highly heterozygous autotetraploid potato." \emph{PloS one} 8, no. 5 (2013): e62355. \doi{10.1371/journal.pone.0062355}}
#' }
#'
"uit"


#' Posterior probabilities from \code{\link{uit}}
#'
#' Contains an array of posterior probabilities of the genotypes from
#' the \code{\link{uit}} dataset. Element \code{gp[i, j, k]} is the
#' posterior probability of dosage \code{k-1} for individual \code{j}
#' at SNP \code{i}.
#'
#' @format A three-dimensional \code{array} object.
#'
#' @source \doi{10.1371/journal.pone.0062355}
#'
#' @references
#' \itemize{
#'   \item{Uitdewilligen, Jan GAML, Anne-Marie A. Wolters, B. Bjorn, Theo JA Borm, Richard GF Visser, and Herman J. Van Eck. "A next-generation sequencing method for genotyping-by-sequencing of highly heterozygous autotetraploid potato." \emph{PloS one} 8, no. 5 (2013): e62355. \doi{10.1371/journal.pone.0062355}}
#' }
#'
#' @seealso \code{\link{uit}} for the full \code{multidog()} fit.
#'
"gp"

#' Genotype log-likelihoods from \code{\link{uit}}
#'
#' Contains an array of genotype log-likelihoods from
#' the \code{\link{uit}} dataset. Element \code{gp[i, j, k]} is the
#' log-likelihood of dosage \code{k-1} for individual \code{j}
#' at SNP \code{i}.
#'
#' @format A three-dimensional \code{array} object.
#'
#' @source \doi{10.1371/journal.pone.0062355}
#'
#' @references
#' \itemize{
#'   \item{Uitdewilligen, Jan GAML, Anne-Marie A. Wolters, B. Bjorn, Theo JA Borm, Richard GF Visser, and Herman J. Van Eck. "A next-generation sequencing method for genotyping-by-sequencing of highly heterozygous autotetraploid potato." \emph{PloS one} 8, no. 5 (2013): e62355. \doi{10.1371/journal.pone.0062355}}
#' }
#'
#' @seealso \code{\link{uit}} for the full \code{multidog()} fit.
#'
"glike"
