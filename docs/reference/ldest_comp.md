# Estimates of composite pairwise LD based either on genotype estimates or genotype likelihoods.

This function will estimate the composite LD between two loci, either
using genotype estimates or using genotype likelihoods. The resulting
measures of LD are generalizations of Burrow's "composite" LD measure.

## Usage

``` r
ldest_comp(
  ga,
  gb,
  K,
  pen = 1,
  useboot = TRUE,
  nboot = 50,
  se = TRUE,
  model = c("norm", "flex")
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

- pen:

  The penalty to be applied to the likelihood. You can think about this
  as the prior sample size. Should be greater than 1. Does not apply if
  `model = "norm"`, `type = "comp"`, and using genotype likelihoods.
  Also does not apply when `type = "comp"` and using genotypes.

- useboot:

  Should we use bootstrap standard errors `TRUE` or not `FALSE`? Only
  applicable if using genotype likelihoods and `model = "flex"`

- nboot:

  The number of bootstrap iterations to use is `boot = TRUE`. Only
  applicable if using genotype likelihoods and `model = "flex"`.

- se:

  A logical. Should we calculate standard errors (`TRUE`) or not
  (`FALSE`). Calculating standard errors can be really slow when
  `type = "comp"`, `model = "flex"`, and when using genotype
  likelihoods. Otherwise, standard error calculations should be pretty
  fast.

- model:

  Should we assume the class of joint genotype distributions is from the
  proportional bivariate normal (`model = "norm"`) or from the general
  categorical distribution (`model = "flex"`). Only applicable if using
  genotype likelihoods.

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

## Composite LD with genotypes
ldout1 <- ldest_comp(ga = ga,
                     gb = gb,
                     K = K)
head(ldout1)
#>            D         D_se           r2        r2_se            r         r_se 
#> 0.0044612795 0.0221319876 0.0004064944 0.0040307019 0.0201617053 0.0999593506 

## Composite LD with genotype likelihoods
ldout2 <- ldest_comp(ga = gamat,
                     gb = gbmat,
                     K = K,
                     se = FALSE,
                     model = "flex")
head(ldout2)
#>           D        D_se          r2       r2_se           r        r_se 
#> 0.008882068          NA 0.018772196          NA 0.137011665          NA 

## Composite LD with genotype likelihoods and proportional bivariate normal
ldout3 <- ldest_comp(ga = gamat,
                     gb = gbmat,
                     K = K,
                     model = "norm")
head(ldout3)
#>          D       D_se         r2      r2_se          r       r_se 
#> 0.02582342 0.01216744 0.67448607 0.26764449 0.82127101 0.16294530 
```
