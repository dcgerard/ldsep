# Obtain shrinkage estimates of correlation from output of [`mldest()`](https://dcgerard.github.io/ldsep/reference/mldest.md) or [`sldest()`](https://dcgerard.github.io/ldsep/reference/sldest.md).

This will take the output of either
[`mldest()`](https://dcgerard.github.io/ldsep/reference/mldest.md) or
[`sldest()`](https://dcgerard.github.io/ldsep/reference/sldest.md),
shrink the Fisher-z transformed correlation estimates using
[`ash()`](https://rdrr.io/pkg/ashr/man/ash.html) (Stephens, 2017; Dey
and Stephens, 2018), then return the corresponding correlation
estimates. You can obtain estimates of r^2 by just squaring these
estimates.

## Usage

``` r
ldshrink(obj, ...)
```

## Arguments

- obj:

  An object of class `lddf`, usually created using either
  [`mldest()`](https://dcgerard.github.io/ldsep/reference/mldest.md) or
  [`sldest()`](https://dcgerard.github.io/ldsep/reference/sldest.md).

- ...:

  Additional arguments to pass to
  [`ash()`](https://rdrr.io/pkg/ashr/man/ash.html).

## Value

A correlation matrix.

## References

- Stephens, Matthew. "False discovery rates: a new deal." Biostatistics
  18, no. 2 (2017): 275-294.

- Dey, Kushal K., and Matthew Stephens. "CorShrink: Empirical Bayes
  shrinkage estimation of correlations, with applications." bioRxiv
  (2018): 368316.

## Author

David Gerard
