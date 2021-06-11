# ldsep 2.1.0

- Includes a `win` argument in `ldfast()`, which implements the moment-based LD correction along a sliding window.
- Fixes a bug where monomorphic SNPs were causing errors when `type = "Dprime` was selected in `ldfast()`. Now we just return `NA`'s for LD with monomorphic SNPs.
- Uses the complete reference of Gerard (2021) [doi:10.1111/1755-0998.13349](https://www.doi.org/10.1111/1755-0998.13349).

# ldsep 2.0.2

- Removes `ldfast_old()` and `ldfast_calc()`, which were not used in any 
  exported functions, because these functions had memory issues, detected 
  by valgrind.

# ldsep 2.0.1

- Added `ldfast()`, a new LD estimation approach based on sample
  moments of marginal posterior genotype moments.
- Unlike `ldest()`, `mldest()`, and `sldest()`, the new approach
  implemented in `ldfast()` is scalable to genome-wide applications,
  as these new estimators can be calculated in linear time in the 
  sample size.
- Citation of MLE approach points to MER article.

# ldsep 1.1.0

* I have changed the terminology from "gametic LD" to "haplotypic" LD,
  and so all instances of "gametic" have changed to "haplotypic". A breaking
  change is that all options that were `"gam"` are now `"hap"`.
* Fixed an issue where the title in `plot.lddf()` was being cut off.
* Added a reference to the preprint where the methodology is developed.
* Updated the vignette to also take a user through uploading a VCF file
  into R using the `VariantAnnotation` package. We also provided
  examples on formatting genotype likelihoods from `updog` and 
  `fitpoly`.

# ldsep 1.0.0

- Initial release of package.
