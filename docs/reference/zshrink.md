# Shrinks Fisher-z transformed correlation estimates and returns resulting correlation estimates.

This function is a wrapper for adaptive shrinkage (Stephens, 2017) on
the Fisher-z transformed estimates of the Pearson correlation. This
approach was proposed in Dey and Stephens (2018) but is re-implemented
here for now since the CorShrink package is not available on CRAN.

## Usage

``` r
zshrink(zmat, smat, ...)
```

## Arguments

- zmat:

  The matrix of Fisher-z transformed correlation estimates.

- smat:

  The matrix of standard errors of the Fisher-z transformed correlation
  estimates.

- ...:

  Additional arguments to pass to
  [`ash()`](https://rdrr.io/pkg/ashr/man/ash.html).

## Value

A matrix of correlation estimates. These are posterior means of the
correlation estimates after applying the CorShrink method (Dey and
Stephens, 2018).

## References

- Stephens, Matthew. "False discovery rates: a new deal." Biostatistics
  18, no. 2 (2017): 275-294.

- Dey, Kushal K., and Matthew Stephens. "CorShrink: Empirical Bayes
  shrinkage estimation of correlations, with applications." bioRxiv
  (2018): 368316.

## Author

David Gerard
