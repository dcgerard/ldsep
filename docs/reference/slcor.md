# Sliding window correlation

Calculates the pairwise Pearson correlation between all columns within a
fixed window size (`win`) using the `use = "pairwise.complete.obs"`
option from [`cor()`](https://rdrr.io/r/stats/cor.html). That is, the
correlation between each pair of variables is computed using all
complete pairs of observations on those variables.

## Usage

``` r
slcor(x, win = 1L)
```

## Arguments

- x:

  A numeric matrix. The variables index the columns.

- win:

  The size of the window. Defaults to 1.

## Value

A correlation matrix with only the observations within a window
containing calculated correlations.

## Author

David Gerard

## Examples

``` r
set.seed(1)
n <- 10
p <- 100
xmat <- matrix(rnorm(n * p), ncol = n)
xmat[sample(n * p, size = 30)] <- NA_real_
slcor(xmat, win = 2)
#>              [,1]        [,2]        [,3]        [,4]        [,5]        [,6]
#>  [1,]  1.00000000 -0.01630404  0.04259209          NA          NA          NA
#>  [2,] -0.01630404  1.00000000 -0.07108527 -0.00698654          NA          NA
#>  [3,]  0.04259209 -0.07108527  1.00000000  0.06059139 -0.02902797          NA
#>  [4,]          NA -0.00698654  0.06059139  1.00000000  0.08572247 -0.13856139
#>  [5,]          NA          NA -0.02902797  0.08572247  1.00000000 -0.04491246
#>  [6,]          NA          NA          NA -0.13856139 -0.04491246  1.00000000
#>  [7,]          NA          NA          NA          NA  0.03849965  0.01209165
#>  [8,]          NA          NA          NA          NA          NA -0.12981911
#>  [9,]          NA          NA          NA          NA          NA          NA
#> [10,]          NA          NA          NA          NA          NA          NA
#>              [,7]        [,8]        [,9]      [,10]
#>  [1,]          NA          NA          NA         NA
#>  [2,]          NA          NA          NA         NA
#>  [3,]          NA          NA          NA         NA
#>  [4,]          NA          NA          NA         NA
#>  [5,]  0.03849965          NA          NA         NA
#>  [6,]  0.01209165 -0.12981911          NA         NA
#>  [7,]  1.00000000  0.12268533 -0.01350859         NA
#>  [8,]  0.12268533  1.00000000 -0.04209290 0.06735998
#>  [9,] -0.01350859 -0.04209290  1.00000000 0.21561209
#> [10,]          NA  0.06735998  0.21561209 1.00000000
```
