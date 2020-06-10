
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ldsep: Linkage Disequilibrium Shrinkage Estimation for Polyploids

<!-- badges: start -->

[![R build
status](https://github.com/dcgerard/ldsep/workflows/R-CMD-check/badge.svg)](https://github.com/dcgerard/ldsep/actions)
[![Codecov test
coverage](https://codecov.io/gh/dcgerard/ldsep/branch/master/graph/badge.svg)](https://codecov.io/gh/dcgerard/ldsep?branch=master)
[![CRAN
status](https://www.r-pkg.org/badges/version/ldsep)](https://CRAN.R-project.org/package=ldsep)
<!-- badges: end -->

Estimates haplotypic or component pairwise linkage disequilibrium (LD),
using either genotypes or genotype likelihoods. Estimates are returned
with standard errors.

The main functions are:

  - `ldest()`: Estimates pairwise LD.
  - `mldest()`: Iteratively apply `ldest()` across many pairs of SNPs.
  - `ldshrink()`: For shrinking correlation estimates using adaptive
    shrinkage.

## Installation

You can install the released version of ldsep from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ldsep")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("dcgerard/ldsep")
```

## Code of Conduct

Please note that the ldsep project is released with a [Contributor Code
of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## References

  - Gerard, David. “Pairwise Linkage Disequilibrium Estimation for
    Polyploids” (unpublished manuscript).
  - Stephens, Matthew. “False discovery rates: a new deal.”
    Biostatistics 18, no. 2 (2017): 275-294.
