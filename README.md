
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ldsep: Linkage Disequilibrium Shrinkage Estimation for Polyploids <a href='https://dcgerard.github.io/ldsep/'><img src='man/figures/logo.png' align="right" height="138" alt="ldsep logo"/></a>

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![R-CMD-check](https://github.com/dcgerard/ldsep/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dcgerard/ldsep/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/dcgerard/ldsep/graph/badge.svg)](https://app.codecov.io/gh/dcgerard/ldsep)
[![CRAN
status](https://www.r-pkg.org/badges/version/ldsep)](https://CRAN.R-project.org/package=ldsep)
[![](http://cranlogs.r-pkg.org/badges/grand-total/ldsep)](https://cran.r-project.org/package=ldsep)
<!-- badges: end -->

Estimate haplotypic or composite pairwise linkage disequilibrium (LD) in
polyploids, using either genotypes or genotype likelihoods. Support is
provided to estimate the popular measures of LD: the LD coefficient $D$,
the standardized LD coefficient $D'$, and the Pearson correlation
coefficient $r$. All estimates are returned with corresponding standard
errors. These estimates and standard errors can then be used for
shrinkage estimation. The methods are described in Gerard (2021a) and
Gerard (2021b).

The main functions are:

- `ldfast()`: Fast, moment-based approach to estimate pairwise LD in the
  presence of genotype uncertainty.
- `ldest()`: Estimates pairwise LD via maximum likelihood.
- `mldest()`: Iteratively apply `ldest()` across many pairs of SNPs.
- `sldest()`: Iteratively apply `ldest()` along a sliding window of
  fixed length.
- `plot.lddf()`: Plot method for the output of `mldest()` and
  `sldest()`.
- `format_lddf()`: Format the output of `mldest()` and `sldest()` into a
  matrix.
- `ldshrink()`: Shrink correlation estimates using adaptive shrinkage
  (Stephens, 2017; Dey and Stephens, 2018).

## Installation

You can install the released version of ldsep from
[CRAN](https://cran.r-project.org/package=ldsep) with:

``` r
install.packages("ldsep")
```

And the development version from
[GitHub](https://github.com/dcgerard/ldsep) with:

``` r
# install.packages("devtools")
devtools::install_github("dcgerard/ldsep")
```

## Citation

To cite `ldsep` in publications use:

> Gerard, David (2021). “Pairwise Linkage Disequilibrium Estimation for
> Polyploids.” *Molecular Ecology Resources*, *21*(4), 1230–1242.
> [doi:10.1111/1755-0998.13349](https://doi.org/10.1111/1755-0998.13349).

A BibTeX entry for LaTeX users is

``` tex
@Article{,
  title = {Pairwise Linkage Disequilibrium Estimation for Polyploids},
  author = {David Gerard},
  journal = {Molecular Ecology Resources},
  year = {2021},
  doi = {10.1111/1755-0998.13349},
  volume = {21},
  number = {4},
  pages = {1230--1242},
}
```

If you use `ldfast()`, please cite:

> Gerard, David (2021). “Scalable Bias-corrected Linkage Disequilibrium
> Estimation Under Genotype Uncertainty.” *Heredity*, *127*(4), 357–362.
> [doi:10.1038/s41437-021-00462-5](https://doi.org/10.1038/s41437-021-00462-5).

A BibTeX entry for LaTeX users is

``` tex
@Article{,
  title = {Scalable Bias-corrected Linkage Disequilibrium Estimation Under Genotype Uncertainty},
  author = {David Gerard},
  journal = {Heredity},
  year = {2021},
  volume = {127},
  number = {4},
  pages = {357--362},
  doi = {10.1038/s41437-021-00462-5},
}
```

## Code of Conduct

Please note that the ldsep project is released with a [Contributor Code
of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## References

- Dey, Kushal K., and Matthew Stephens (2018). “CorShrink: Empirical
  Bayes shrinkage estimation of correlations, with applications.”
  *bioRxiv*. [doi:10.1101/368316](https://doi.org/10.1101/368316)

- Gerard, David (2021a). “Pairwise Linkage Disequilibrium Estimation for
  Polyploids.” *Molecular Ecology Resources*, *21*(4), 1230–1242.
  [doi:10.1111/1755-0998.13349](https://doi.org/10.1111/1755-0998.13349).

- Gerard, David (2021b). “Scalable Bias-corrected Linkage Disequilibrium
  Estimation Under Genotype Uncertainty.” *Heredity*, *127*(4), 357–362.
  [doi:10.1038/s41437-021-00462-5](https://doi.org/10.1038/s41437-021-00462-5).

- Stephens, Matthew (2017). “False discovery rates: a new deal.”
  *Biostatistics* 18(2), 275–294.
  [doi:10.1093/biostatistics/kxw041](https://doi.org/10.1093/biostatistics/kxw041)
