# Fast LD Estimation

## Abstract

We demonstrate our new moment-based linkage disequilibrium (LD)
estimators, as implemented in the
[`ldfast()`](https://dcgerard.github.io/ldsep/reference/ldfast.md)
function. These methods are described in detail in Gerard (2021). The
following code chunks are only evaluated if
[updog](https://dcgerard.github.io/updog/) is installed with a version
of at least 2.0.2.

## Analysis

We will load the packages and data that we will need

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`updog`](https://dcgerard.github.io/updog/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`ldsep`](https://dcgerard.github.io/ldsep/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`corrplot`](https://github.com/taiyun/corrplot)`)`\
`#> corrplot 0.95 loaded`\
[`data`](https://rdrr.io/r/utils/data.html)`(``"uit"``, package ``=`` ``"ldsep"``)`\
[`class`](https://rdrr.io/r/base/class.html)`(``uit``)`\
`#> [1] "multidog"`

We need the posterior probabilities for each genotype for each
individual at each SNP. In [updog](https://dcgerard.github.io/updog/),
this corresponds to the `Pr_n` variables in its output, where `n` varies
from 0 to the ploidy. We can use
[`updog::format_multidog()`](https://dcgerard.github.io/updog/reference/format_multidog.html)
to extract these values into an array with the proper dimensions.

\
`ploidy`` ``<-`` ``4`\
`gp`` ``<-`` `[`format_multidog`](https://dcgerard.github.io/updog/reference/format_multidog.html)`(``x ``=`` ``uit``, varname ``=`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"Pr_"``, ``0``:``ploidy``)``)`\
[`class`](https://rdrr.io/r/base/class.html)`(``gp``)`\
`#> [1] "array"`\
[`dim`](https://rdrr.io/r/base/dim.html)`(``gp``)`\
`#> [1] 10 84  5`

The first dimension indexes the SNPs, the second dimension indexes the
individuals, and the third dimension indexes the possible dosages. So
`gp[i, j, k]` contains the posterior probability that individual `j` has
dosage `k-1` at SNP `i`.

You can use other genotyping programs to obtain these posterior
probabilities, but they have to be formatted as an array with the
dimensions in the correct order (SNPs by individuals by dosages).

**However**, in order for these estimates to perform well, you need to
use posterior genotype probabilities that have been calculated using
adaptive priors, i.e. empirical/hierarchical Bayes approaches. There are
many approaches that do this, such as
[`{updog}`](https://cran.r-project.org/package=updog),
[`{polyRAD}`](https://cran.r-project.org/package=polyRAD),
[`{fitPoly}`](https://cran.r-project.org/package=fitPoly), or
[`{SuperMASSA}`](https://github.com/guilherme-pereira/vcf2sm). Note that
GATK uses a uniform prior, so would be inappropriate for use in
[`ldfast()`](https://dcgerard.github.io/ldsep/reference/ldfast.md).

Let’s use our genotype posterior probabilities (`gp`) to estimate LD.

\
`ldout`` ``<-`` `[`ldfast`](https://dcgerard.github.io/ldsep/reference/ldfast.md)`(``gp ``=`` ``gp``, type ``=`` ``"r2"``)`

`ldout` contains:

- `ldmat`: The estimated squared correlation coefficients, adjusted for
  genotype uncertainty.

  \
  [`corrplot`](https://rdrr.io/pkg/corrplot/man/corrplot.html)`(``corr ``=`` ``ldout``$``ldmat``, `\
  `         method ``=`` ``"color"``, `\
  `         type ``=`` ``"upper"``,`\
  `         diag ``=`` ``FALSE``,`\
  `         tl.pos ``=`` ``"n"``)`

  ![Heat map of r-squared](fast_files/figure-html/unnamed-chunk-4-1.png)

- `semat`: The estimated standard errors for the squared correlation
  coefficient estimates.

  \
  [`corrplot`](https://rdrr.io/pkg/corrplot/man/corrplot.html)`(``corr ``=`` ``ldout``$``semat``, `\
  `         method ``=`` ``"color"``, `\
  `         type ``=`` ``"upper"``,`\
  `         diag ``=`` ``FALSE``,`\
  `         tl.pos ``=`` ``"n"``, `\
  `         is.corr ``=`` ``FALSE``)`

  ![Heat map of standard errors of
  r-squared](fast_files/figure-html/unnamed-chunk-5-1.png)

- `rr`: The estimated “reliability ratios”, the multiplicative effects
  applied to the naive LD estimates using just the posterior mean
  genotypes.

  \
  `graphics``::`[`hist`](https://rdrr.io/r/graphics/hist.html)`(``x ``=`` ``ldout``$``rr``, `\
  `               main ``=`` ``"Reliability Ratios"``, `\
  `               xlab ``=`` ``"Reliability Ratios"``)`

  ![Histogram of reliability
  ratios](fast_files/figure-html/unnamed-chunk-6-1.png)

The slowest part of
[`ldfast()`](https://dcgerard.github.io/ldsep/reference/ldfast.md) is
using hierarchical shrinkage on the log of the reliability ratios
(`shrinkrr = TRUE`), but I wouldn’t disable this as some reliability
ratios have high variance and are greatly improved through shrinkage.

Other LD estimators are available: `D`, `Dprime`, `r`, and `z`. See the
help file of
[`ldfast()`](https://dcgerard.github.io/ldsep/reference/ldfast.md) for
details.

## References

- Gerard, David (2021). “Scalable Bias-corrected Linkage Disequilibrium
  Estimation Under Genotype Uncertainty.” *Heredity*, *127*(4), 357–362.
  [doi:10.1038/s41437-021-00462-5](https://doi.org/10.1038/s41437-021-00462-5).
