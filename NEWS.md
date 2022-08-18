# snpR 1.2.3

## Features

### Major
* Added `calc_prop_poly()` to calculate proportion of polymorphic loci for a given pop.
* Renamed `plot_tree()` to `calc_tree()` and removed automatic plot generation, since `ggtree`, which that depended on, can behave a bit oddly sometimes. An example for generating a plot with `ggtree` was added to the examples section of `calc_tree()`'s documentation. `ggtree` is no longer a suggested dependency.

### Minor
* Added support for FWE (multiple comparisons p-value adjustment) to `filter_snps()` HWE filtering.
* Added singleton fitering to `filter_snps()`.
* Added iPCA `ncp` and `ncp.max` options to `run_random_forest()` and `run_genomic_prediction()`. Previously, selecting the `iPCA` option would work, but run with the default `iPCA` options and thus determine `ncp` internally, which is pretty slow.
* Added rug plotting to `plot_manhattan()` to allow for easy plotting of gene positions, etc under plots. Both ribbon-style and classic rug style supported.

## Documentation
* Fixed some outdated documentation in `calc_ne()`.
* Fixed some typos in the documentation.
* Fixed some parameter input descriptions in the documentation for Colony.
* Cleaned up some documentation for fwe correction across multiple functions.

## Bug fixes
* Fixed a bug with `calc_ne()` where pop names were not being properly handled and het/coan method results were not being returned by `get.snpR.stats()`.
* Fixed a rare internal bug where `.get.task.list()` would sometimes add a space between facet levels when pasting due to weird `t()` bug. To the user, this might have occasionally resulted in weird behavior when using multiple SNP facets with mixed numeric and character classes.
* Fixed a bug with `RefManageR` apparently introduced recently where `bibentry` objects wouldn't correctly write. Eliminated `RefManageR` dependency, now just uses `rbibutils` for reading and writing, and, when calling `citations()`, just spits out the full citation rather than the inline for the "Citation: " line. Not ideal, but doesn't need `RefManageR`.
* Fixed a bug with `calc_genetic_distances()` where the "Nei" method wouldn't work (due to a typo in the code).
* Fixed an issue where the ID column wasn't written to vcf files.
* Fixed interpolation not checking for `missMDA` installation or reporting `iPCA` as an option if an invalid method provided.
* Fixed a bug where the `subfacet` and `facet` columns in the `stats` slot of a `snpRdata` object would get flipped in order during `calc_association()`, resulting in the addition of a bunch of empty data rows when something else was merged in. This would produce a downstream error during `calc_pairwise_fst` (and potentially elsewhere) due to attempting to cbind the `stats` slot to the `facet_meta` slot. Note that this wouldn't cause any bad stats to be calculated or returned anywhere due to the way that `get.snpR.stats()` functions to fetch results!

# snpR 1.2.2

## UI

## Features

### Major
* Added `read_structure()` to read STRUCTURE formatted files. Added auto-read of .str files to `import.snpR.data()`.
* Added `calc_abba_baba()` to do ABBA/BABA tests, including block jackknifing for significance.
* Added `calc_diagnostic()` to plot basic diagnostic plots for snpRdata objects.

### Minor
* Now redistributing part of Roy Francis's `pophelper` package, since it's not on CRAN. Switched to GPL license to allow this. The package is still cited automatically when generating a citation during `plot_structure()` or `plot_structure_map`(). Added a dependency on `stats` to allow for the re-packaging.
* Cleanup option added to `calc_ne()`.
* Removed dependencies: `stringi`, `pkgcond`, `stringr`, `tidyr`, `CATT`. All had only one or a few used functions that were not time intensive, and so could be home-brewed easily to avoid the dependency.
* Adjusted `stickSNPs` to be smaller--only 100 loci and 100 samples now.
* Added a `verbose` argument to `calc_ne()`.

### Documentation
* Removed the `snpR_association` vignette, since it required the `GMMAT`, `BGLR`, and `ranger` R packages to be installed, but they are only suggested, not required. This could cause vignette building to fail. May re-tool later to just use the internally implemented association tests.
* Updated documentation for sfs related functions to more explicitly talk about ref and anc columns for folding.

### Bug fixes
* Fixed a bug where an incorrect armitage stat would be calculated with `calc_association` if the major allele differed between the case and control.
* `calc_ne()` and `run_colony()` actually cleanup now if asked on Windows.
* Fixed the `verbose` argument on `run_colony()` to actually work on Windows.
* Fixed a bug trying to merge new and old data with the same columns in different orders could cause NAs to be filled by the wrong numbers. Internal only problem that only effected the new `calc_abba_baba()` function.


# snpR 1.2.1

## UI

* Added a verbose option to many functions to avoid putting quite so much noise onto the console.

## Features

### Major

* Added `calc_he()` for traditional H<sub>E</sub> = 2pq calculation. Note that this produces results *almost* identical to `calc_pi()`.
* Added `calc_hs()` for Coltman et al (1999)'s individual heterozygosity.

### Minor

* Renamed `calc_SFS()` and `calc_pairwise_LD_heatmap()` to lowercase for consistency.
* Got rid of `sfs` arguments in `plot_sfs()` and `calc_direcitonality()`. They now just take provided sfs objects as x, consistent with other functions.
* Reworked `calc_pairwise_fst()` bootstrapping to be more memory efficient.
* Added a `chr_order` argument to `plot_manhattan()` to allow for manual resorting of chromosomes (since factors are coerced away in `snpRdata` objects).
* Added a `highlight_style` argument to `plot_manhattan()` to allow for coloring SNPs instead of labeling them if highlighted.
* Added a `verbose` option (defaulting to `TRUE`) for `filter_snps()` to suppress all of the filtering reports.

## Documentation

* Added a `NEWS.md` file to track changes to the package.
* Cleaned up the association testing vignette.

## Bug fixes

* Fixed a bug in `run_random_forest()` where formulas would be incorrectly specified the first time a line of code was run due to a weird environment scope issue.
* Fixed a bug in `get.snpR.stats()` when requesting fst values from both a facet with and without fst calculated would throw an error during fst matrix construction. Implemented a test.
* Added ggtree citation to plot_tree.
* Fixed a bug where running no formula but generating importance estimates would throw an error in `run_random_forest()`.
* Fixed a bug where filtering a `snpRdata` object with `filter_snps()` such that no individuals or SNPs remained would result in an uninformative error. Added a test to check error messages here and in `susbet_snpR_data()` for this.
* Fixed a bug in `format_snps()` to allow for BYmax and BYmin columns to get translated to BY.max and BY.min, since periods aren't allowed in snpRdata metadata columns.
* Fixed a bug in `calc_het_hom_ratio()` complex facets would throw an error.
