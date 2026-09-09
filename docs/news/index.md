# Changelog

## ldsep 2.1.7

- Remove dependency on [modeest](https://github.com/paulponcet/modeest)
  by implementing my own version of
  [`hsm()`](https://dcgerard.github.io/ldsep/reference/hsm.md).

## ldsep 2.1.6

CRAN release: 2025-09-18

- Remove dependency on C++11.

## ldsep 2.1.5

CRAN release: 2022-10-18

- Bug fix: Use `&&` instead of `&` in C++ source.

## ldsep 2.1.4

CRAN release: 2022-02-11

- Use unnamed formula argument in
  [`aggregate()`](https://rdrr.io/r/stats/aggregate.html) to account for
  change in R.

## ldsep 2.1.2

CRAN release: 2021-08-10

- `DOUBLE_EPS` -\> `DOUBLE_EPSILON`.
- Modified citation to point to *Heredity*.
- Adds notes about the necessity of using adaptive priors in
  [`ldfast()`](https://dcgerard.github.io/ldsep/reference/ldfast.md).

## ldsep 2.1.0

CRAN release: 2021-06-11

- Includes a `win` argument in
  [`ldfast()`](https://dcgerard.github.io/ldsep/reference/ldfast.md),
  which implements the moment-based LD correction along a sliding
  window.
- Fixes a bug where monomorphic SNPs were causing errors when
  `type = "Dprime"` was selected in
  [`ldfast()`](https://dcgerard.github.io/ldsep/reference/ldfast.md).
  Now we just return `NA`’s for LD with monomorphic SNPs.
- Uses the complete reference of Gerard (2021)
  (<doi:10.1111/1755-0998.13349>)

## ldsep 2.0.2

CRAN release: 2021-02-17

- Removes `ldfast_old()` and `ldfast_calc()`, which were not used in any
  exported functions, because these functions had memory issues,
  detected by valgrind.

## ldsep 2.0.1

CRAN release: 2021-02-12

- Added
  [`ldfast()`](https://dcgerard.github.io/ldsep/reference/ldfast.md), a
  new LD estimation approach based on sample moments of marginal
  posterior genotype moments.
- Unlike
  [`ldest()`](https://dcgerard.github.io/ldsep/reference/ldest.md),
  [`mldest()`](https://dcgerard.github.io/ldsep/reference/mldest.md),
  and
  [`sldest()`](https://dcgerard.github.io/ldsep/reference/sldest.md),
  the new approach implemented in
  [`ldfast()`](https://dcgerard.github.io/ldsep/reference/ldfast.md) is
  scalable to genome-wide applications, as these new estimators can be
  calculated in linear time in the sample size.
- Citation of MLE approach points to MER article.

## ldsep 1.1.0

CRAN release: 2020-11-14

- I have changed the terminology from “gametic LD” to “haplotypic” LD,
  and so all instances of “gametic” have changed to “haplotypic”. A
  breaking change is that all options that were `"gam"` are now `"hap"`.
- Fixed an issue where the title in
  [`plot.lddf()`](https://dcgerard.github.io/ldsep/reference/plot.lddf.md)
  was being cut off.
- Added a reference to the preprint where the methodology is developed.
- Updated the vignette to also take a user through uploading a VCF file
  into R using the `VariantAnnotation` package. We also provided
  examples on formatting genotype likelihoods from `updog` and
  `fitpoly`.

## ldsep 1.0.0

CRAN release: 2020-08-02

- Initial release of package.
