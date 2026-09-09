# Fast bias-correction for LD Estimation

Estimates the reliability ratios from posterior marginal moments and
uses these to correct the biases in linkage disequilibrium estimation
caused by genotype uncertainty. These methods are described in Gerard
(2021).

## Usage

``` r
ldfast(
  gp,
  type = c("r", "r2", "z", "D", "Dprime"),
  shrinkrr = TRUE,
  se = TRUE,
  thresh = TRUE,
  upper = 10,
  mode = c("zero", "estimate"),
  win = NULL
)
```

## Arguments

- gp:

  A three-way array with dimensions SNPs by individuals by dosage. That
  is, `gp[i, j, k]` is the posterior probability of dosage `k-1` for
  individual `j` at SNP `i`.

- type:

  What LD measure should we estimate?

  `"r"`

  : The Pearson correlation.

  `"r2"`

  : The squared Pearson correlation.

  `"z"`

  : The Fisher-z transformed Pearson correlation.

  `"D"`

  : The LD coefficient.

  `"Dprime"`

  : The standardized LD coefficient.

  Note that these are all *composite* measures of LD (see the
  description in
  [`ldest()`](https://dcgerard.github.io/ldsep/reference/ldest.md)).

- shrinkrr:

  A logical. Should we use adaptive shrinkage (Stephens, 2016) to shrink
  the reliability ratios (`TRUE`) or keep the raw reliability ratios
  (`FALSE`). Defaults to `TRUE`.

- se:

  Should we also return a matrix of standard errors (`TRUE`) or not
  (`FALSE`)? It is faster to not return standard errors. Defaults to
  `TRUE`.

- thresh:

  A logical. Should we apply an upper bound on the reliability ratios
  (`TRUE`) or not (`FALSE`).

- upper:

  The upper bound on the reliability ratios if `thresh = TRUE`. The
  default is a generous 10.

- mode:

  A character. Only applies if `shrinkrr = TRUE`. When using
  hierarchical shrinkage on the log of the reliability ratios, should we
  use zero as the mode (`mode = "zero"`) or estimate it using the
  procedure of Robertson and Cryer (1974) (`mode = "estimate"`)?

- win:

  A positive integer. The window size. This will constrain the
  correlations calculated to those +/- the window size. This will only
  improve speed if the window size is *much* less than the number of
  SNPs.

## Value

A list with some or all of the following elements:

- `ldmat`:

  The bias-corrected LD matrix.

- `rr`:

  The estimated reliability ratio for each SNP. This is the
  multiplicative factor applied to the naive LD estimate for each SNP.

- `rr_raw`:

  The raw reliability ratios (for the covariance, not the correlation).
  Only returned if `shrinkrr = TRUE`.

- `rr_se`:

  The standard errors for the *log*-raw reliability ratios for each SNP.
  That is, we have sd(log(rr_raw)) ~ rr_se. Only returned if
  `shrinkrr = TRUE`.

- `semat`:

  A matrix of standard errors of the corresponding estimators of LD.

## Details

Returns consistent and bias-corrected estimates of linkage
disequilibrium. The usual measures of LD are implemented: D, D', r, r2,
and z (Fisher-z of r). These are all *composite* measures of LD, not
haplotypic measures of LD (see the description in
[`ldest()`](https://dcgerard.github.io/ldsep/reference/ldest.md)). They
are always appropriate measures of association between loci, but only
correspond to haplotypic measures of LD when Hardy-Weinberg equilibrium
is fulfilled in autopolyploids.

In order for these estimates to perform well, you need to use posterior
genotype probabilities that have been calculated using adaptive priors,
i.e. empirical/hierarchical Bayes approaches. There are many approaches
that do this, such as
[`updog`](https://cran.r-project.org/package=updog),
[`polyRAD`](https://cran.r-project.org/package=polyRAD),
[`fitPoly`](https://cran.r-project.org/package=fitPoly), or
[`SuperMASSA`](https://github.com/guilherme-pereira/vcf2sm). Note that
GATK uses a uniform prior, so would be inappropriate for use in
`ldfast()`.

Calculating standard errors and performing hierarchical shrinkage of the
reliability ratios are both rather slow operations compared to just raw
method-of-moments based estimation for LD. If you don't need standard
errors, you can double your speed by setting `se = FALSE`. It is not
recommended that you disable the hierarchical shrinkage.

Due to sampling variability, the estimates sometime lie outside of the
theoretical boundaries of the parameters being estimated. In such cases,
we truncate the estimates at the boundary and return `NA` for the
standard errors.

## Mathematical formulation

Let

- \\r\\ be the sample correlation of posterior mean genotypes between
  loci 1 and 2,

- \\a1\\ be the sample variance of posterior means at locus 1,

- \\a2\\ be the sample variance of posterior means at locus 2,

- \\b1\\ be the sample mean of posterior variances at locus 1, and

- \\b2\\ be the sample mean of posterior variances at locus 2.

Then the estimated Pearson correlation between the genotypes at loci 1
and 2 is \$\$\sqrt{(a1 + b1)/a1}\sqrt{(a2 + b2)/a2}r.\$\$ All other LD
calculations are based on this equation. In particular, the estimated
genotype variances at loci 1 and 2 are \\a1 + b1\\ and \\a2 + b2\\,
respectively, which can be used to calculate D and D'.

The reliability ratio for SNP i is defined by \\(ai + bi)/ai\\. By
default, we apply [`ash()`](https://rdrr.io/pkg/ashr/man/ash.html)
(Stephens, 2016) to the log of these reliability ratios before adjusting
the Pearson correlation. Standard errors are required before using
[`ash()`](https://rdrr.io/pkg/ashr/man/ash.html), but these are easily
obtained using the central limit theorem and the delta-method.

## References

- Gerard, David. Scalable Bias-corrected Linkage Disequilibrium
  Estimation Under Genotype Uncertainty. *Heredity*, 127(4),
  357–362, 2021.
  [doi:10.1038/s41437-021-00462-5](https://doi.org/10.1038/s41437-021-00462-5)
  .

- T. Robertson and J. D. Cryer. An iterative procedure for estimating
  the mode. *Journal of the American Statistical Association*,
  69(348):1012–1016, 1974.
  [doi:10.1080/01621459.1974.10480246](https://doi.org/10.1080/01621459.1974.10480246)
  .

- M. Stephens. False discovery rates: a new deal. *Biostatistics*,
  18(2):275–294, 10 2016.
  [doi:10.1093/biostatistics/kxw041](https://doi.org/10.1093/biostatistics/kxw041)
  .

## See also

- [`ash()`](https://rdrr.io/pkg/ashr/man/ash.html):

  Function used to perform hierarchical shrinkage on the log of the
  reliability ratios.

- [`ldest()`](https://dcgerard.github.io/ldsep/reference/ldest.md),
  [`mldest()`](https://dcgerard.github.io/ldsep/reference/mldest.md),
  [`sldest()`](https://dcgerard.github.io/ldsep/reference/sldest.md):

  Maximum likelihood estimation of linkage disequilibrium.

## Author

David Gerard

## Examples

``` r
data("gp")

ldout <- ldfast(gp, "r")
ldout$ldmat
#>            [,1]      [,2]       [,3]         [,4]      [,5]      [,6]
#>  [1,] 1.0000000 1.0000000 0.43269941  0.831256306 0.6028276 0.6028276
#>  [2,] 1.0000000 1.0000000 0.43269941  0.831256306 0.6028276 0.6028276
#>  [3,] 0.4326994 0.4326994 1.00000000  0.573455065 0.8927188 0.8927188
#>  [4,] 0.8312563 0.8312563 0.57345506  1.000000000 0.5484621 0.5484621
#>  [5,] 0.6028276 0.6028276 0.89271875  0.548462093 1.0000000 1.0000000
#>  [6,] 0.6028276 0.6028276 0.89271875  0.548462093 1.0000000 1.0000000
#>  [7,] 0.4281823 0.4281823 0.97896686  0.581033314 0.8885195 0.8885195
#>  [8,] 0.4281823 0.4281823 0.97896686  0.581033314 0.8885195 0.8885195
#>  [9,] 0.1036140 0.1036140 0.09310332 -0.001987959 0.1340953 0.1340953
#> [10,] 0.5571129 0.5571129 0.82400859  0.540610084 0.9605484 0.9605484
#>             [,7]       [,8]         [,9]     [,10]
#>  [1,] 0.42818226 0.42818226  0.103613990 0.5571129
#>  [2,] 0.42818226 0.42818226  0.103613990 0.5571129
#>  [3,] 0.97896686 0.97896686  0.093103317 0.8240086
#>  [4,] 0.58103331 0.58103331 -0.001987959 0.5406101
#>  [5,] 0.88851952 0.88851952  0.134095304 0.9605484
#>  [6,] 0.88851952 0.88851952  0.134095304 0.9605484
#>  [7,] 1.00000000 1.00000000  0.079131179 0.8512614
#>  [8,] 1.00000000 1.00000000  0.079131179 0.8512614
#>  [9,] 0.07913118 0.07913118  1.000000000 0.1020791
#> [10,] 0.85126139 0.85126139  0.102079061 1.0000000
ldout$rr
#>  [1] 1.005771 1.005771 1.003719 1.003953 1.001617 1.001617 1.003960 1.003960
#>  [9] 1.062503 1.032458
ldout$semat
#>             [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
#>  [1,]         NA         NA 0.09332343 0.03572633 0.06723054 0.06723054
#>  [2,]         NA         NA 0.09332343 0.03572633 0.06723054 0.06723054
#>  [3,] 0.09332343 0.09332343         NA 0.06374431 0.02960981 0.02960981
#>  [4,] 0.03572633 0.03572633 0.06374431         NA 0.06308449 0.06308449
#>  [5,] 0.06723054 0.06723054 0.02960981 0.06308449         NA         NA
#>  [6,] 0.06723054 0.06723054 0.02960981 0.06308449         NA         NA
#>  [7,] 0.09234704 0.09234704 0.01306726 0.06062951 0.02913919 0.02913919
#>  [8,] 0.09234704 0.09234704 0.01306726 0.06062951 0.02913919 0.02913919
#>  [9,] 0.11330929 0.11330929 0.11338916 0.10962925 0.10304775 0.10304775
#> [10,] 0.07796815 0.07796815 0.04525697 0.06934184 0.01956130 0.01956130
#>             [,7]       [,8]      [,9]      [,10]
#>  [1,] 0.09234704 0.09234704 0.1133093 0.07796815
#>  [2,] 0.09234704 0.09234704 0.1133093 0.07796815
#>  [3,] 0.01306726 0.01306726 0.1133892 0.04525697
#>  [4,] 0.06062951 0.06062951 0.1096292 0.06934184
#>  [5,] 0.02913919 0.02913919 0.1030478 0.01956130
#>  [6,] 0.02913919 0.02913919 0.1030478 0.01956130
#>  [7,]         NA         NA 0.1152616 0.04149810
#>  [8,]         NA         NA 0.1152616 0.04149810
#>  [9,] 0.11526155 0.11526155        NA 0.11264141
#> [10,] 0.04149810 0.04149810 0.1126414         NA

ldout <- ldfast(gp, "D")
ldout$ldmat
#>             [,1]       [,2]       [,3]          [,4]       [,5]       [,6]
#>  [1,] 0.24316990 0.24316990 0.11133968  0.1890926874 0.14190970 0.14190970
#>  [2,] 0.24316990 0.24316990 0.11133968  0.1890926874 0.14190970 0.14190970
#>  [3,] 0.11133968 0.11133968 0.27228111  0.1380362082 0.22237575 0.22237575
#>  [4,] 0.18909269 0.18909269 0.13803621  0.2127990758 0.12078013 0.12078013
#>  [5,] 0.14190970 0.14190970 0.22237575  0.1207801316 0.22789145 0.22789145
#>  [6,] 0.14190970 0.14190970 0.22237575  0.1207801316 0.22789145 0.22789145
#>  [7,] 0.11108594 0.11108594 0.26875234  0.1410137301 0.22315493 0.22315493
#>  [8,] 0.11108594 0.11108594 0.26875234  0.1410137301 0.22315493 0.22315493
#>  [9,] 0.01152076 0.01152076 0.01095422 -0.0002067759 0.01443395 0.01443395
#> [10,] 0.12264546 0.12264546 0.19195255  0.1113326112 0.20470887 0.20470887
#>              [,7]        [,8]          [,9]      [,10]
#>  [1,] 0.111085938 0.111085938  0.0115207554 0.12264546
#>  [2,] 0.111085938 0.111085938  0.0115207554 0.12264546
#>  [3,] 0.268752337 0.268752337  0.0109542225 0.19195255
#>  [4,] 0.141013730 0.141013730 -0.0002067759 0.11133261
#>  [5,] 0.223154931 0.223154931  0.0144339511 0.20470887
#>  [6,] 0.223154931 0.223154931  0.0144339511 0.20470887
#>  [7,] 0.276790383 0.276790383  0.0093870857 0.19993638
#>  [8,] 0.276790383 0.276790383  0.0093870857 0.19993638
#>  [9,] 0.009387086 0.009387086  0.0508411413 0.01027538
#> [10,] 0.199936377 0.199936377  0.0102753754 0.19929976
ldout$rr
#>  [1] 1.011575 1.011575 1.007451 1.007921 1.003237 1.003237 1.007936 1.007936
#>  [9] 1.128912 1.065970
ldout$semat
#>             [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
#>  [1,]         NA         NA 0.03145416 0.03197325 0.03082308 0.03082308
#>  [2,]         NA         NA 0.03145416 0.03197325 0.03082308 0.03082308
#>  [3,] 0.03145416 0.03145416         NA 0.02713108 0.03349861 0.03349861
#>  [4,] 0.03197325 0.03197325 0.02713108         NA 0.02596956 0.02596956
#>  [5,] 0.03082308 0.03082308 0.03349861 0.02596956         NA         NA
#>  [6,] 0.03082308 0.03082308 0.03349861 0.02596956         NA         NA
#>  [7,] 0.03101405 0.03101405 0.03348837 0.02704158 0.03328525 0.03328525
#>  [8,] 0.03101405 0.03101405 0.03348837 0.02704158 0.03328525 0.03328525
#>  [9,] 0.01277533 0.01277533 0.01328870 0.01127557 0.01137312 0.01137312
#> [10,] 0.02665070 0.02665070 0.02646479 0.02337153 0.02738681 0.02738681
#>             [,7]       [,8]       [,9]      [,10]
#>  [1,] 0.03101405 0.03101405 0.01277533 0.02665070
#>  [2,] 0.03101405 0.03101405 0.01277533 0.02665070
#>  [3,] 0.03348837 0.03348837 0.01328870 0.02646479
#>  [4,] 0.02704158 0.02704158 0.01127557 0.02337153
#>  [5,] 0.03328525 0.03328525 0.01137312 0.02738681
#>  [6,] 0.03328525 0.03328525 0.01137312 0.02738681
#>  [7,]         NA         NA 0.01360690 0.02608113
#>  [8,]         NA         NA 0.01360690 0.02608113
#>  [9,] 0.01360690 0.01360690         NA 0.01132382
#> [10,] 0.02608113 0.02608113 0.01132382         NA

ldout <- ldfast(gp, "Dprime")
ldout$ldmat
#>             [,1]       [,2]       [,3]          [,4]      [,5]      [,6]
#>  [1,] 4.00000000 0.99234186 0.68861855  0.9333391533 1.1182071 1.1182071
#>  [2,] 0.99234186 4.00000000 0.68861855  0.9333391533 1.1182071 1.1182071
#>  [3,] 0.68861855 0.68861855 4.00000000  1.0445593471 1.3592571 1.3592571
#>  [4,] 0.93333915 0.93333915 1.04455935  4.0000000000 1.1644397 1.1644397
#>  [5,] 1.11820712 1.11820712 1.35925711  1.1644397472 4.0000000 1.2815177
#>  [6,] 1.11820712 1.11820712 1.35925711  1.1644397472 1.2815177 4.0000000
#>  [7,] 0.67835286 0.67835286 1.29638385  1.0535843861 1.3714162 1.3714162
#>  [8,] 0.67835286 0.67835286 1.29638385  1.0535843861 1.3714162 1.3714162
#>  [9,] 0.04904818 0.04904818 0.06492008 -0.0008680408 0.1089843 0.1089843
#> [10,] 0.96807877 0.96807877 1.17532263  1.0752091511 1.1494195 1.1494195
#>             [,7]       [,8]          [,9]      [,10]
#>  [1,] 0.67835286 0.67835286  0.0490481819 0.96807877
#>  [2,] 0.67835286 0.67835286  0.0490481819 0.96807877
#>  [3,] 1.29638385 1.29638385  0.0649200760 1.17532263
#>  [4,] 1.05358439 1.05358439 -0.0008680408 1.07520915
#>  [5,] 1.37141620 1.37141620  0.1089843435 1.14941949
#>  [6,] 1.37141620 1.37141620  0.1089843435 1.14941949
#>  [7,] 4.00000000 1.32871897  0.0549282885 1.23084578
#>  [8,] 1.32871897 4.00000000  0.0549282885 1.23084578
#>  [9,] 0.05492829 0.05492829  4.0000000000 0.07771872
#> [10,] 1.23084578 1.23084578  0.0777187155 4.00000000
ldout$rr
#>  [1] 1.011575 1.011575 1.007451 1.007921 1.003237 1.003237 1.007936 1.007936
#>  [9] 1.128912 1.065970
ldout$semat
#>             [,1]       [,2]       [,3]      [,4]       [,5]       [,6]
#>  [1,]         NA 0.13995724 0.19150934 0.1526757 0.21646378 0.21646378
#>  [2,] 0.13995724         NA 0.19150934 0.1526757 0.21646378 0.21646378
#>  [3,] 0.19150934 0.19150934         NA 0.1889567 0.17151103 0.17151103
#>  [4,] 0.15267573 0.15267573 0.18895674        NA 0.22732911 0.22732911
#>  [5,] 0.21646378 0.21646378 0.17151103 0.2273291         NA 0.16529356
#>  [6,] 0.21646378 0.21646378 0.17151103 0.2273291 0.16529356         NA
#>  [7,] 0.18673492 0.18673492 0.14851597 0.1876919 0.17147607 0.17147607
#>  [8,] 0.18673492 0.18673492 0.14851597 0.1876919 0.17147607 0.17147607
#>  [9,] 0.05460539 0.05460539 0.07852983 0.0473413 0.08452456 0.08452456
#> [10,] 0.19999687 0.19999687 0.15592059 0.2163060 0.14446347 0.14446347
#>             [,7]       [,8]       [,9]      [,10]
#>  [1,] 0.18673492 0.18673492 0.05460539 0.19999687
#>  [2,] 0.18673492 0.18673492 0.05460539 0.19999687
#>  [3,] 0.14851597 0.14851597 0.07852983 0.15592059
#>  [4,] 0.18769185 0.18769185 0.04734130 0.21630603
#>  [5,] 0.17147607 0.17147607 0.08452456 0.14446347
#>  [6,] 0.17147607 0.17147607 0.08452456 0.14446347
#>  [7,]         NA 0.14500916 0.07949034 0.15287303
#>  [8,] 0.14500916         NA 0.07949034 0.15287303
#>  [9,] 0.07949034 0.07949034         NA 0.08538073
#> [10,] 0.15287303 0.15287303 0.08538073         NA
```
