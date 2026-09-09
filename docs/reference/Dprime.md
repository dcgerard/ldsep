# Get the standardized composite D'.

This function will either standardize by the maximum covariance
conditional on the marginal genotype distribution, or by the maximum
covariance conditional on the marginal allele frequencies.

## Usage

``` r
Dprime(qmat, type = c("allele", "geno"), constrain = FALSE)
```

## Arguments

- qmat:

  The observed joint genotype distribution.

- type:

  Should we condition on the marginal genotype distribution
  (`type = "geno"`), or should we condition on the allele frequency
  (`type = "allele"`)?

- constrain:

  A logical. This option is only applicable when `type = "allele"`.
  Should return an value that is equal to D' under HWE (`FALSE`) or a
  value that is constrained to lie between -1 and 1 (`TRUE`)? Defaults
  to `FALSE`.

## Value

A vector of length 2. The first element is the estimated D'. The second
element is the normalization used.

## Details

Note that when `type = "allele"` and `constrain = FALSE`, the resulting
D' is constrained to fall between -K and K, where K is the ploidy of the
species. However, under HWE, this measure is equal to haplotypic D'.
Using `constrain = TRUE` will result in a measure that is constrained to
lie between -1 and 1, but it will not equal haplotypic D' under HWE.

Using `type = "geno"` is its own thing and will not equal D' generally
under HWE. When `type = "geno"`, then the the `constrain` parameter has
no effect.

## Author

David Gerard

## Examples

``` r
K <- 6
qmat <- matrix(stats::runif((K+1)^2), nrow = K+1)
qmat <- qmat / sum(qmat)
Dprime(qmat, type = "geno")
#>     Dprime       Dmax 
#> 0.08697558 0.66851339 
Dprime(qmat, type = "allele")
#>    Dprime      Dmax 
#> 0.2644744 0.2198487 
```
