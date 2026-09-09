# Calculate prior variances from a matrix of prior genotype probabilities.

Given a matrix of prior probabilities for the genotypes at each SNP,
this function will calculate the prior variance of genotypes.

## Usage

``` r
pvcalc(priormat)
```

## Arguments

- priormat:

  A matrix of prior genotype probabilities. Element `priormat[i, j]` is
  the prior probability of dosage `j` at SNP `i`.

## Value

A vector of prior variances.

## Author

David Gerard

## Examples

``` r
data("uit")
priormat <- uit$snpdf[, paste0("Pr_", 0:4)]
pvcalc(priormat)
#>       248       249       251       252       253       254       257       258 
#> 0.9612945 0.9612945 1.0763415 0.8411998 0.9008180 0.9008180 1.0941723 1.0941723 
#>       262       266 
#> 0.2012746 0.7882709 
```
