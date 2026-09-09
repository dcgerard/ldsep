# Format an element of [`mldest()`](https://dcgerard.github.io/ldsep/reference/mldest.md) or [`sldest()`](https://dcgerard.github.io/ldsep/reference/sldest.md) into an upper-triangular matrix.

Formats the LD estimates and standard errors output from running
[`mldest()`](https://dcgerard.github.io/ldsep/reference/mldest.md) or
[`sldest()`](https://dcgerard.github.io/ldsep/reference/sldest.md) into
a more conventional upper-triangular matrix.

## Usage

``` r
format_lddf(obj, element = "r2")
```

## Arguments

- obj:

  An object of class `lddf`, usually output from running either
  [`mldest()`](https://dcgerard.github.io/ldsep/reference/mldest.md) or
  [`sldest()`](https://dcgerard.github.io/ldsep/reference/sldest.md).

- element:

  Which element in `obj` should we format into an upper-triangular
  matrix?

## Value

A matrix of the selected elements. Only the upper-triangle of the matrix
is filled. The lower-triangle and the diagonal are `NA`'s.

## Author

David Gerard

## Examples

``` r
set.seed(1)

## Simulate genotypes when true correlation is 0
nloci <- 5
nind  <- 100
K <- 6
nc <- 1
genomat <- matrix(sample(0:K, nind * nloci, TRUE), nrow = nloci)

## Haplotypic LD estimates
lddf <- mldest(geno = genomat,
               K = K,
               nc = nc,
               type = "hap")

## Obtain the D estimates in matrix form
Dmat <- format_lddf(obj = lddf, element = "D")
Dmat
#>    1         2          3            4           5
#> 1 NA -0.015268 0.01127912 -0.004729285 -0.02260301
#> 2 NA        NA 0.00394976  0.009062196 -0.02153841
#> 3 NA        NA         NA  0.023425711 -0.01808722
#> 4 NA        NA         NA           NA -0.01764326
#> 5 NA        NA         NA           NA          NA
```
