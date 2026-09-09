# Estimate haplotypic pair-wise LD using either genotypes or genotype likelihoods.

Given genotype (allele dosage) or genotype likelihood data for each
individual at a pair of loci, this function will calculate the maximum
likelihood estimates and their corresponding asymptotic standard errors
of some measures of linkage disequilibrium (LD): D, D', the Pearson
correlation, the squared Pearson correlation, and the Fisher-z
transformation of the Pearson correlation. This function can be used for
both diploids and polyploids.

## Usage

``` r
ldest_hap(
  ga,
  gb,
  K,
  reltol = 10^-8,
  nboot = 100,
  useboot = FALSE,
  pen = 2,
  grid_init = FALSE,
  se = TRUE
)
```

## Arguments

- ga:

  One of two possible inputs:

  1.  A vector of counts, containing the genotypes for each individual
      at the first locus. When `type = "comp"`, the vector of genotypes
      may be continuous (e.g. the posterior mean genotype).

  2.  A matrix of genotype log-likelihoods at the first locus. The rows
      index the individuals and the columns index the genotypes. That is
      `ga[i, j]` is the genotype likelihood of individual `i` for
      genotype `j-1`.

- gb:

  One of two possible inputs:

  1.  A vector of counts, containing the genotypes for each individual
      at the second locus. When `type = "comp"`, the vector of genotypes
      may be continuous (e.g. the posterior mean genotype).

  2.  A matrix of genotype log-likelihoods at the second locus. The rows
      index the individuals and the columns index the genotypes. That is
      `gb[i, j]` is the genotype likelihood of individual `i` for
      genotype `j-1`.

- K:

  The ploidy of the species. Assumed to be the same for all individuals.

- reltol:

  The relative tolerance for the stopping criterion.

- nboot:

  Sometimes, the MLE standard errors don't exist. So we use the
  bootstrap as a backup. `nboot` specifies the number of bootstrap
  iterations.

- useboot:

  A logical. Optionally, you may always use the bootstrap to estimate
  the standard errors (`TRUE`). These will be more accurate but also
  much slower, so this defaults to `FALSE`. Only applicable if using
  genotype likelihoods.

- pen:

  The penalty to be applied to the likelihood. You can think about this
  as the prior sample size. Should be greater than 1. Does not apply if
  `model = "norm"`, `type = "comp"`, and using genotype likelihoods.
  Also does not apply when `type = "comp"` and using genotypes.

- grid_init:

  A logical. Should we initialize the gradient ascent at a grid of
  initial values (`TRUE`) or just initialize at one value corresponding
  to the simplex point `rep(0.25, 4)` (`FALSE`)?

- se:

  A logical. Should we calculate standard errors (`TRUE`) or not
  (`FALSE`). Calculating standard errors can be really slow when
  `type = "comp"`, `model = "flex"`, and when using genotype
  likelihoods. Otherwise, standard error calculations should be pretty
  fast.

## Value

A vector with some or all of the following elements:

- `D`:

  The estimate of the LD coefficient.

- `D_se`:

  The standard error of the estimate of the LD coefficient.

- `r2`:

  The estimate of the squared Pearson correlation.

- `r2_se`:

  The standard error of the estimate of the squared Pearson correlation.

- `r`:

  The estimate of the Pearson correlation.

- `r_se`:

  The standard error of the estimate of the Pearson correlation.

- `Dprime`:

  The estimate of the standardized LD coefficient. When `type` = "comp",
  this corresponds to the standardization where we fix allele
  frequencies.

- `Dprime_se`:

  The standard error of `Dprime`.

- `Dprimeg`:

  The estimate of the standardized LD coefficient. This corresponds to
  the standardization where we fix genotype frequencies.

- `Dprimeg_se`:

  The standard error of `Dprimeg`.

- `z`:

  The Fisher-z transformation of `r`.

- `z_se`:

  The standard error of the Fisher-z transformation of `r`.

- `p_ab`:

  The estimated haplotype frequency of ab. Only returned if estimating
  the haplotypic LD.

- `p_Ab`:

  The estimated haplotype frequency of Ab. Only returned if estimating
  the haplotypic LD.

- `p_aB`:

  The estimated haplotype frequency of aB. Only returned if estimating
  the haplotypic LD.

- `p_AB`:

  The estimated haplotype frequency of AB. Only returned if estimating
  the haplotypic LD.

- `q_ij`:

  The estimated frequency of genotype i at locus 1 and genotype j at
  locus 2. Only returned if estimating the composite LD.

- `n`:

  The number of individuals used to estimate pairwise LD.

## Details

Let A and a be the reference and alternative alleles, respectively, at
locus 1. Let B and b be the reference and alternative alleles,
respectively, at locus 2. Let paa, pAb, paB, and pAB be the frequencies
of haplotypes ab, Ab, aB, and AB, respectively. Let pA = pAb + pAB and
let pB = paB + pAB The `ldest` returns estimates of the following
measures of LD.

- D: pAB - pA pB

- D': D / Dmax, where Dmax = min(pA pB, (1 - pA) (1 - pB)) if D \< 0 and
  Dmax = min(pA (1 - pB), pA (1 - pB)) if D \> 0

- r-squared: The squared Pearson correlation, r^2 = D^2 / (pA (1 - pA)
  pB (1 - pB))

- r: The Pearson correlation, r = D / sqrt(pA (1 - pA) pB (1 - pB))

Estimates are obtained via maximum likelihood under the assumption of
Hardy-Weinberg equilibrium. The likelihood is calculated by integrating
over the possible haplotypes for each pair of genotypes.

The resulting standard errors are based on the square roots of the
inverse of the negative Fisher-information. This is from standard
maximum likelihood theory. The Fisher-information is known to be biased
low, so the actual standard errors are probably a little bigger for
small n (n \< 20). In some cases the Fisher-information matrix is
singular, and so we in these cases we return a bootstrap estimate of the
standard error.

The standard error estimate of the squared Pearson correlation is not
valid when r^2 = 0.

In cases where either SNP is estimated to be monoallelic
(`pA %in% c(0, 1)` or `pB %in% c(0, 1)`), this function will return LD
estimates of `NA`.

## Author

David Gerard

## Examples

``` r
set.seed(1)
n <- 100 # sample size
K <- 6 # ploidy

## generate some fake genotypes when LD = 0.
ga <- stats::rbinom(n = n, size = K, prob = 0.5)
gb <- stats::rbinom(n = n, size = K, prob = 0.5)
head(ga)
#> [1] 2 3 3 5 2 5
head(gb)
#> [1] 3 3 2 6 3 2

## generate some fake genotype likelihoods when LD = 0.
gamat <- t(sapply(ga, stats::dnorm, x = 0:K, sd = 1, log = TRUE))
gbmat <- t(sapply(gb, stats::dnorm, x = 0:K, sd = 1, log = TRUE))
head(gamat)
#>            [,1]      [,2]       [,3]       [,4]      [,5]       [,6]      [,7]
#> [1,]  -2.918939 -1.418939 -0.9189385 -1.4189385 -2.918939 -5.4189385 -8.918939
#> [2,]  -5.418939 -2.918939 -1.4189385 -0.9189385 -1.418939 -2.9189385 -5.418939
#> [3,]  -5.418939 -2.918939 -1.4189385 -0.9189385 -1.418939 -2.9189385 -5.418939
#> [4,] -13.418939 -8.918939 -5.4189385 -2.9189385 -1.418939 -0.9189385 -1.418939
#> [5,]  -2.918939 -1.418939 -0.9189385 -1.4189385 -2.918939 -5.4189385 -8.918939
#> [6,] -13.418939 -8.918939 -5.4189385 -2.9189385 -1.418939 -0.9189385 -1.418939
head(gbmat)
#>            [,1]       [,2]       [,3]       [,4]      [,5]      [,6]       [,7]
#> [1,]  -5.418939  -2.918939 -1.4189385 -0.9189385 -1.418939 -2.918939 -5.4189385
#> [2,]  -5.418939  -2.918939 -1.4189385 -0.9189385 -1.418939 -2.918939 -5.4189385
#> [3,]  -2.918939  -1.418939 -0.9189385 -1.4189385 -2.918939 -5.418939 -8.9189385
#> [4,] -18.918939 -13.418939 -8.9189385 -5.4189385 -2.918939 -1.418939 -0.9189385
#> [5,]  -5.418939  -2.918939 -1.4189385 -0.9189385 -1.418939 -2.918939 -5.4189385
#> [6,]  -2.918939  -1.418939 -0.9189385 -1.4189385 -2.918939 -5.418939 -8.9189385

## Haplotypic LD with genotypes
ldout1 <- ldest_hap(ga = ga,
                    gb = gb,
                    K = K)
head(ldout1)
#>            D         D_se           r2        r2_se            r         r_se 
#> 0.0058365062 0.0286113348 0.0005454794 0.0053480254 0.0233554997 0.1144917775 

## Haplotypic LD with genotype likelihoods
ldout2 <- ldest_hap(ga = gamat,
                    gb = gbmat,
                    K = K)
head(ldout2)
#>           D        D_se          r2       r2_se           r        r_se 
#> 0.020950505 0.097754998 0.007028751 0.065592168 0.083837645 0.391185656 
```
