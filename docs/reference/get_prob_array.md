# Obtain the distribution of genotypes given haplotype frequencies under HWE

This function will calculate the (log) probabilities for all genotype
combinations at two loci given just the haplotype frequencies. This is
under the assumptions of HWE.

## Usage

``` r
get_prob_array(K, prob, log_p = TRUE)
```

## Arguments

- K:

  The ploidy of the species.

- prob:

  Haplotype frequencies in the order of (ab, Ab, aB, AB).

- log_p:

  A logical. Should we return the log-probabilities (`TRUE`) or the
  probabilities (`FALSE`). Defaults to `TRUE`.

## Value

Element (i, j) is the (log) probability of genotype i-1 at locus 1 and
genotype j-1 at locus 2.

## Author

David Gerard

## Examples

``` r
get_prob_array(K = 6, prob = c(0.1, 0.2, 0.3, 0.4), log_p = FALSE)
#>          [,1]     [,2]     [,3]    [,4]     [,5]     [,6]     [,7]
#> [1,] 0.000001 0.000018 0.000135 0.00054 0.001215 0.001458 0.000729
#> [2,] 0.000012 0.000204 0.001440 0.00540 0.011340 0.012636 0.005832
#> [3,] 0.000060 0.000960 0.006360 0.02232 0.043740 0.045360 0.019440
#> [4,] 0.000160 0.002400 0.014880 0.04880 0.089280 0.086400 0.034560
#> [5,] 0.000240 0.003360 0.019440 0.05952 0.101760 0.092160 0.034560
#> [6,] 0.000192 0.002496 0.013440 0.03840 0.061440 0.052224 0.018432
#> [7,] 0.000064 0.000768 0.003840 0.01024 0.015360 0.012288 0.004096
```
