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

* Initial release of package.
