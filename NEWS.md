# snpR 1.2.5

## Features

### Major


### Minor


## Documentation
* Changed "sample metadata" to "SNP metadata" in the description of the `header_cols` argument to `import.snpR.data()`.
* Edited the warning message when importing structure data without setting header_cols properly to more explicitly reference that problem.
* Updated the `calc_tajimas_d()` examples to use the updated `get.snpR.stats()` syntax.

## Bug fixes
* Changed `import.snpR.data()` to return an error if the genotypic dimensions are not correct for the provided SNP and sample meta data instead of proceeding and returning a bogus result.
* Fixed a bug where not using a pop facet with `plot_structure()` with the `method = "structure"` option would cause a bug due to an incorrectly set "LOCISPOP" flag (which is checked even if not using the locprior option).


# snpR 1.2.4 -- hot-fix

## Bug fixes
* Fixed a nasty bug during imputation of sn data where the wrong maf was found. Shouldn't directly bias clustering, but good to fix asap.


# snpR 1.2.3

## Features

### Major
* Added `calc_prop_poly()` to calculate proportion of polymorphic loci for a given pop.
* Renamed `plot_tree()` to `calc_tree()` and removed automatic plot generation, since `ggtree`, which that depended on, can behave a bit oddly sometimes. An example for generating a plot with `ggtree` was added to the examples section of `calc_tree()`'s documentation. `ggtree` is no longer a suggested dependency.
* Added `summarize_facets()`, which provides information on either the possible facets or provided specific facets for a given
`snpRdata` object.

### Minor
* Added support for FWE (multiple comparisons p-value adjustment) to `filter_snps()` HWE filtering.
* Added singleton filtering to `filter_snps()`.
* Added iPCA `ncp` and `ncp.max` options to `run_random_forest()` and `run_genomic_prediction()`. Previously, selecting the `iPCA` option would work, but run with the default `iPCA` options and thus determine `ncp` internally, which is pretty slow.
* Added rug plotting to `plot_manhattan()` to allow for easy plotting of gene positions, etc under plots. Both ribbon-style and classic rug style supported.
* Moved code for `run_sequoia()` and `plot_structure_map()` to a secondary github repo, since both of these use CRAN unfriendly dependencies (sequoia and sf, respectively). These functions now source this code (after asking for permission) in order to run, allowing the dependencies to be dropped. This should not change the user experience whatsoever.
* Set `get.snpR.stats()` to return an informative warning if an empty list is returned.
* Added GenAlEx outfile support to `format_snps()`. `openxlsx` added to the `Suggests` field of the `DESCRIPTION` as needed for `.xlsx` file creation.

## Documentation
* Fixed some outdated documentation in `calc_ne()`.
* Fixed some typos in the documentation.
* Fixed some parameter input descriptions in the documentation for Colony.
* Cleaned up some documentation for FWE correction across multiple functions.
* Fixed spelling issues throughout package.
* Added a note on import format guessing, cleaned up the file format descriptions to more clearly note the extension `snpR` expects, and added the missing STRUCTRE import file description to the documentation for `import.snpR.data()`.

## Bug fixes
* Fixed a bug with `calc_ne()` where pop names were not being properly handled and het/coan method results were not being returned by `get.snpR.stats()`.
* Fixed a rare internal bug where `.get.task.list()` would sometimes add a space between facet levels when pasting due to weird `t()` bug. To the user, this might have occasionally resulted in weird behavior when using multiple SNP facets with mixed numeric and character classes.
* Fixed a bug with `RefManageR` apparently introduced recently where `bibentry` objects wouldn't correctly write. Eliminated `RefManageR` dependency, now just uses `rbibutils` for reading and writing, and, when calling `citations()`, just spits out the full citation rather than the inline for the "Citation: " line. Not ideal, but doesn't need `RefManageR`.
* Fixed a bug with `calc_genetic_distances()` where the "Nei" method wouldn't work (due to a typo in the code).
* Fixed an issue where the ID column wasn't written to vcf files.
* Fixed interpolation not checking for `missMDA` installation or reporting `iPCA` as an option if an invalid method provided.
* Fixed a bug where the `subfacet` and `facet` columns in the `stats` slot of a `snpRdata` object would get flipped in order during `calc_association()`, resulting in the addition of a bunch of empty data rows when something else was merged in. This would produce a downstream error during `calc_pairwise_fst` (and potentially elsewhere) due to attempting to `cbind` the `stats` slot to the `facet_meta` slot. Note that this wouldn't cause any bad stats to be calculated or returned anywhere due to the way that `get.snpR.stats()` functions to fetch results!
* Fixed a bug where trying to calculate windowed averages on top of previously calculated Tajima's D would throw a merge error.
* Fixed a bug where `plot_manhattan()` wouldn't correctly plot Tajima's D (didn't look for it in the right place)
* Fixed a bug where `plot_manhattan()` would display the facet name even if there was only one possible level (for example, if the sample facet was the .base level).
* Fixed a bug where importing VCF files with `read_vcf()` would fail due to a typo in a sanity check.
* Fixed a bug where reading in files with, for example, "GC" AND "CG" genotypes somewhere would cause tabulated genotypes to treat these as two different genotypes. Not a common bug except where reading in, for example, VCF files.
* Fixed a bug where asking to return fst for a facet that it hasn't been calculated would return an uniformative error instead of an empty list and a warning.
* Not so much a bug, but added a warning when reading in genepop/fstat/etc file formats that do not specify the allelic identity (ATCG), sine snpR assumes A/C SNPs in this case. This could potentially cause downstream issues with *other* packages or with stuff like `calc_sfs()` with `fold = FALSE` and `calc_abba_baba()`.

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
