# Pairwise LD estimation in polyploids.

Estimates either haplotypic or composite measures of LD using either
genotypes are genotype likelihoods via maximum likelihood. The usual
measures of LD are estimated (D, D', and r) along with the Fisher-z
transformation of r (called "z"). All estimates are returned with
standard errors. See Gerard (2021) for details.

## Usage

``` r
ldest(
  ga,
  gb,
  K,
  se = TRUE,
  type = c("hap", "comp"),
  model = c("norm", "flex"),
  pen = ifelse(type == "hap", 2, 1)
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

- se:

  A logical. Should we calculate standard errors (`TRUE`) or not
  (`FALSE`). Calculating standard errors can be really slow when
  `type = "comp"`, `model = "flex"`, and when using genotype
  likelihoods. Otherwise, standard error calculations should be pretty
  fast.

- type:

  The type of LD to calculate. The available types are haplotypic LD
  (`type = "hap"`) or composite LD (`type = "comp"`). Haplotypic LD is
  only appropriate for autopolyploids when the individuals are in
  Hardy-Weinberg equilibrium (HWE). The composite measures of LD are
  always applicable, and consistently estimate the usual measures of LD
  when HWE is fulfilled in autopolyploids. However, when HWE is not
  fulfilled, interpreting the composite measures of LD could be a little
  tricky.

- model:

  When `type = "comp"` and using genotype likelihoods, should we use the
  proportional bivariate normal model to estimate the genotype
  distribution (`model = "norm"`), or the general categorical
  distribution (`model = "flex"`)? Defaults to `"norm"`.

- pen:

  The penalty to be applied to the likelihood. You can think about this
  as the prior sample size. Should be greater than 1. Does not apply if
  `model = "norm"`, `type = "comp"`, and using genotype likelihoods.
  Also does not apply when `type = "comp"` and using genotypes.

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

## Haplotypic LD

This section describes the methods used when `type = "hap"` is selected.

Haplotypic LD measures the association between two loci on the same
haplotype. When haplotypes are known, estimating haplotypic LD is simple
using just the haplotypic frequencies.

When haplotypes are not known, we can still estimate haplotypic
frequencies using the genotypes or genotype likelihoods *in
autopolyploids as long as Hardy-Weinberg equilibrium (HWE) is
satisfied*. We do this via maximum likelihood using gradient ascent.
Gradient ascent is performed over the unconstrained parameterization of
the 3-simplex from Betancourt (2012). The estimated haplotype
frequencies are then used to estimate haplotypic LD.

Standard errors are provided using standard maximum likelihood theory.
In brief, the Hessian matrix of the log-likelihood is calculated at the
MLE's of the haplotype frequencies. The negative inverse of this Hessian
matrix is approximately the covariance matrix of the MLE's of the
haplotype frequencies. Since all haplotypic LD measures are functions of
the haplotype frequencies, we use the delta-method to obtain the
standard errors for each LD estimate.

A Dirichlet(2,2,2,2) prior is placed over the frequencies of haplotypes
00, 01, 10, and 11. This corresponds to the "add two" rule of Agresti
(1998). You can change this prior via the `pen` argument.

When you either do not have autopolyploids or when HWE is *not*
satisfied, then the estimates using `type = "hap"` are nonsensical.
However, the composite measures of LD are still applicable (see below).

## Composite LD

This section describes the methods used when `type = "comp"` is
selected.

When HWE is not satisfied, haplotype frequencies are not estimable.
However, measures of association between two loci are still estimable.
These associations may be caused by LD either on the same haplotype or
between different haplotypes. Cockerham and Weir (1977) thus called such
measures "composite" measures of LD.

When the genotypes are known, these composite measures have simple
correspondences to well-known statistical measures of association. D is
the covariance of genotypes between loci divided by the ploidy. r is the
Pearson correlation of genotypes. D' is D divided by a term that
involves only mean genotypes.

When genotypes are not known, we estimate the joint genotype frequencies
and use these to estimate the composite measures of LD using genotype
likelihoods. The distribution of genotypes is assumed to either follow a
proportional bivariate normal model (by default) or a general
categorical model.

These estimates of composite measures of LD estimate the haplotypic
measures of LD when HWE is fulfilled, but are still applicable when HWE
is not fulfilled.

When genotypes are known, standard errors are calculated using standard
moment-based approaches. When genotypes are not known, standard errors
are calculated using standard maximum likelihood theory, same as for the
haplotypic LD estimates (see above), or using a bootstrap.

## References

- Agresti, Alan, and Brent A. Coull. "Approximate is better than "exact"
  for interval estimation of binomial proportions." *The American
  Statistician* 52, no. 2 (1998): 119-126.
  [doi:10.1080/00031305.1998.10480550](https://doi.org/10.1080/00031305.1998.10480550)

- Betancourt, Michael. "Cruising the simplex: Hamiltonian Monte Carlo
  and the Dirichlet distribution." In *AIP Conference Proceedings 31st*,
  vol. 1443, no. 1, pp. 157-164. American Institute of Physics, 2012.
  [doi:10.1063/1.3703631](https://doi.org/10.1063/1.3703631)

- Cockerham, C. Clark, and B. S. Weir. "Digenic descent measures for
  finite populations." *Genetics Research* 30, no. 2 (1977): 121-147.
  [doi:10.1017/S0016672300017547](https://doi.org/10.1017/S0016672300017547)

- Gerard, David. "Pairwise Linkage Disequilibrium Estimation for
  Polyploids." *Molecular Ecology Resources* 21, no. 4 (2021):
  1230-1242.
  [doi:10.1111/1755-0998.13349](https://doi.org/10.1111/1755-0998.13349)

## See also

- [`ldfast()`](https://dcgerard.github.io/ldsep/reference/ldfast.md):

  Fast, moment-based approach to LD estimation that also accounts for
  genotype uncertainty.

- [`mldest()`](https://dcgerard.github.io/ldsep/reference/mldest.md):

  For calculating pairwise LD among all pairs of a collection of SNPs.

- [`sldest()`](https://dcgerard.github.io/ldsep/reference/sldest.md):

  For calculating pairwise LD along a sliding window of SNPs.

- [`ldest_hap()`](https://dcgerard.github.io/ldsep/reference/ldest_hap.md):

  For the function that directly estimates haplotypic LD when HWE is
  fulfilled.

- [`ldest_comp()`](https://dcgerard.github.io/ldsep/reference/ldest_comp.md):

  For the function that directly estimates composite LD.

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
ldout1 <- ldest(ga = ga,
                gb = gb,
                K = K,
                type = "hap")
head(ldout1)
#>            D         D_se           r2        r2_se            r         r_se 
#> 0.0058365062 0.0286113348 0.0005454794 0.0053480254 0.0233554997 0.1144917775 

## Haplotypic LD with genotype likelihoods
ldout2 <- ldest(ga = gamat,
                gb = gbmat,
                K = K,
                type = "hap")
head(ldout2)
#>           D        D_se          r2       r2_se           r        r_se 
#> 0.020950505 0.097754998 0.007028751 0.065592168 0.083837645 0.391185656 

## Composite LD with genotypes
ldout3 <- ldest(ga = ga,
                gb = gb,
                K = K,
                type = "comp")
head(ldout3)
#>            D         D_se           r2        r2_se            r         r_se 
#> 0.0044612795 0.0221319876 0.0004064944 0.0040307019 0.0201617053 0.0999593506 

## Composite LD with genotype likelihoods and normal model
ldout4 <- ldest(ga = gamat,
                gb = gbmat,
                K = K,
                type = "comp",
                model = "norm")
head(ldout4)
#>          D       D_se         r2      r2_se          r       r_se 
#> 0.02582342 0.01216744 0.67448607 0.26764449 0.82127101 0.16294530 

## Composite LD with genotype likelihoods and general categorical model
ldout5 <- ldest(ga = gamat,
                gb = gbmat,
                K = K,
                type = "comp",
                model = "flex",
                se = FALSE)
head(ldout5)
#>           D        D_se          r2       r2_se           r        r_se 
#> 0.008882068          NA 0.018772196          NA 0.137011665          NA 

ldout1[["D"]]
#> [1] 0.005836506
ldout2[["D"]]
#> [1] 0.02095051
ldout3[["D"]]
#> [1] 0.004461279
ldout4[["D"]]
#> [1] 0.02582342
ldout5[["D"]]
#> [1] 0.008882068
```
