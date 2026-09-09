# Returns distribution of proportional bivariate normal.

Returns distribution of proportional bivariate normal.

## Usage

``` r
pbnorm_dist(mu, sigma, K, log = FALSE)
```

## Arguments

- mu:

  A vector of length 2 containing the mean.

- sigma:

  A 2-by-2 positive definite covariance matrix

- K:

  The ploidy of the individual.

- log:

  A logical. If `TRUE`, log probabilities are returned.

## Value

A matrix. Element (i,j) is the (log) probability of genotype i-1 at
locus 1 and j-1 at locus 2.

## Author

David Gerard
