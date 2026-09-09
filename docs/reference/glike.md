# Genotype log-likelihoods from [`uit`](https://dcgerard.github.io/ldsep/reference/uit.md)

Contains an array of genotype log-likelihoods from the
[`uit`](https://dcgerard.github.io/ldsep/reference/uit.md) dataset.
Element `gp[i, j, k]` is the log-likelihood of dosage `k-1` for
individual `j` at SNP `i`.

## Usage

``` r
glike
```

## Format

A three-dimensional `array` object.

## Source

[doi:10.1371/journal.pone.0062355](https://doi.org/10.1371/journal.pone.0062355)

## References

- Uitdewilligen, Jan GAML, Anne-Marie A. Wolters, B. Bjorn, Theo JA
  Borm, Richard GF Visser, and Herman J. Van Eck. "A next-generation
  sequencing method for genotyping-by-sequencing of highly heterozygous
  autotetraploid potato." *PloS one* 8, no. 5 (2013): e62355.
  [doi:10.1371/journal.pone.0062355](https://doi.org/10.1371/journal.pone.0062355)

## See also

[`uit`](https://dcgerard.github.io/ldsep/reference/uit.md) for the full
`multidog()` fit.
