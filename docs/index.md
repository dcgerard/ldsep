# ldsep: Linkage Disequilibrium Shrinkage Estimation for Polyploids

Estimate haplotypic or composite pairwise linkage disequilibrium (LD) in
polyploids, using either genotypes or genotype likelihoods. Support is
provided to estimate the popular measures of LD: the LD coefficient
$`D`$, the standardized LD coefficient $`D'`$, and the Pearson
correlation coefficient $`r`$. All estimates are returned with
corresponding standard errors. These estimates and standard errors can
then be used for shrinkage estimation. The methods are described in
Gerard (2021a) and Gerard (2021b).

The main functions are:

- [`ldfast()`](https://dcgerard.github.io/ldsep/reference/ldfast.md):
  Fast, moment-based approach to estimate pairwise LD in the presence of
  genotype uncertainty.
- [`ldest()`](https://dcgerard.github.io/ldsep/reference/ldest.md):
  Estimates pairwise LD via maximum likelihood.
- [`mldest()`](https://dcgerard.github.io/ldsep/reference/mldest.md):
  Iteratively apply
  [`ldest()`](https://dcgerard.github.io/ldsep/reference/ldest.md)
  across many pairs of SNPs.
- [`sldest()`](https://dcgerard.github.io/ldsep/reference/sldest.md):
  Iteratively apply
  [`ldest()`](https://dcgerard.github.io/ldsep/reference/ldest.md) along
  a sliding window of fixed length.
- [`plot.lddf()`](https://dcgerard.github.io/ldsep/reference/plot.lddf.md):
  Plot method for the output of
  [`mldest()`](https://dcgerard.github.io/ldsep/reference/mldest.md) and
  [`sldest()`](https://dcgerard.github.io/ldsep/reference/sldest.md).
- [`format_lddf()`](https://dcgerard.github.io/ldsep/reference/format_lddf.md):
  Format the output of
  [`mldest()`](https://dcgerard.github.io/ldsep/reference/mldest.md) and
  [`sldest()`](https://dcgerard.github.io/ldsep/reference/sldest.md)
  into a matrix.
- [`ldshrink()`](https://dcgerard.github.io/ldsep/reference/ldshrink.md):
  Shrink correlation estimates using adaptive shrinkage (Stephens, 2017;
  Dey and Stephens, 2018).

## Installation

You can install the released version of ldsep from
[CRAN](https://cran.r-project.org/package=ldsep) with:

\
[`install.packages`](https://rdrr.io/r/utils/install.packages.html)`(``"ldsep"``)`

And the development version from
[GitHub](https://github.com/dcgerard/ldsep) with:

\
`# install.packages("devtools")`\
`devtools``::`[`install_github`](https://devtools.r-lib.org/reference/install-deprecated.html)`(``"dcgerard/ldsep"``)`

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

If you use
[`ldfast()`](https://dcgerard.github.io/ldsep/reference/ldfast.md),
please cite:

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
