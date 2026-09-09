# Sliding window LD estimation

This function is a wrapper for
[`ldest()`](https://dcgerard.github.io/ldsep/reference/ldest.md) for
estimating LD along a sliding window of a fixed size. Support is
provided for parallelization through the foreach and doParallel
packages.

## Usage

``` r
sldest(
  geno,
  K,
  win = 50,
  nc = 1,
  type = c("hap", "comp"),
  model = c("norm", "flex"),
  pen = ifelse(type == "hap", 2, 1),
  se = TRUE
)
```

## Arguments

- geno:

  One of two possible inputs:

  - A matrix of genotypes (allele dosages). The rows index the SNPs and
    the columns index the individuals. That is, `genomat[i, j]` is the
    allele dosage for individual `j` in SNP `i`. When `type = "comp"`,
    the dosages are allowed to be continuous (e.g. posterior mean
    genotypes).

  - A three-way array of genotype *log*-likelihoods. The first dimension
    indexes the SNPs, the second dimension indexes the individuals, and
    the third dimension indexes the genotypes. That is,
    `genolike_array[i, j, k]` is the genotype log-likelihood at SNP `i`
    for individual `j` and dosage `k - 1`.

- K:

  The ploidy of the species. Assumed to be the same for all individuals.

- win:

  The window size. Pairwise LD will be estimated plus or minus these
  many positions. Larger sizes significantly increase the computational
  load.

- nc:

  The number of computing cores to use. This should never be more than
  the number of cores available in your computing environment. You can
  determine the maximum number of available cores by running
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  in R. This is probably fine for a personal computer, but some
  environments are only able to use fewer. Ask your admins if you are
  unsure.

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

- se:

  A logical. Should we calculate standard errors (`TRUE`) or not
  (`FALSE`). Calculating standard errors can be really slow when
  `type = "comp"`, `model = "flex"`, and when using genotype
  likelihoods. Otherwise, standard error calculations should be pretty
  fast.

## Value

A data frame of class `c("lddf", "data.frame")` with some or all of the
following elements:

- `i`:

  The index of the first SNP.

- `j`:

  The index of the second SNP.

- `snpi`:

  The row name corresponding to SNP `i`, if row names are provided.

- `snpj`:

  The row name corresponding to SNP `j`, if row names are provided.

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

See [`ldest()`](https://dcgerard.github.io/ldsep/reference/ldest.md) for
details on the different types of LD estimators supported.

## See also

- [`ldest()`](https://dcgerard.github.io/ldsep/reference/ldest.md):

  For the base function that estimates pairwise LD.

- [`mldest()`](https://dcgerard.github.io/ldsep/reference/mldest.md):

  For estimating pairwise LD between *all* provided SNPs.

- [`ldfast()`](https://dcgerard.github.io/ldsep/reference/ldfast.md):

  Fast, moment-based approach to LD estimation that also accounts for
  genotype uncertainty.

- [`format_lddf()`](https://dcgerard.github.io/ldsep/reference/format_lddf.md):

  For formatting the output of `sldest()` as a matrix.

- [`plot.lddf()`](https://dcgerard.github.io/ldsep/reference/plot.lddf.md):

  For plotting the output of `sldest()`.

## Author

David Gerard

## Examples

``` r
set.seed(1)

## Simulate genotypes when true correlation is 0
nloci <- 100
nind  <- 100
win <- 5
K <- 6
nc <- 1
genomat <- matrix(sample(0:K, nind * nloci, TRUE), nrow = nloci)

## Composite LD estimates
lddf <- sldest(geno = genomat,
               K = K,
               win = win,
               nc = nc,
               type = "comp")
plot(lddf, element = "z")

```
