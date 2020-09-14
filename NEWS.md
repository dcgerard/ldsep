# ldsep 2.0.0

- Added `ldfast()`, a new LD estimation approach based on sample
  moments of marginal posterior genotype moments.
- Unlike `ldest()`, `mldest()`, and `sldest()`, the new approach
  implemented in `ldfast()` is scalable to genome-wide applications,
  as these new estimators can be calculated in linear time in the 
  sample size.

# ldsep 1.0.0

- Initial release of package.
