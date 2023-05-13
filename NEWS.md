# snpR 1.2.7
## Features

### Minor
* Changed the behavior of `calc_fis()` to take the ratio of average variance components instead of the average of ratios following the recommendations of [Bahtia et al. 2013](https://doi.org/10.1101/gr.154831.113) and in line with the new behavior of `calc_pairwise_fst()`.
* To aid in calculating weighted mean $F_{IS}$ values for cases where the dataset is too big to run at once with snpR (as may be the case for large WGS data sets, for example), the `keep_components` argument was added to `calc_fis()` to return the "b" and "c" variance components for each locus for later processing. Brief instructions were added to the documentation for `calc_fis()` to explain this process. This brings `calc_fis()` behavior fully in line with `calc_pairwise_fst()` for bi-allelic markers, although it still needs to be fixed for poly-allelic markers (which are not yet supported on the front-end).
* Changed the internal behavior of `sample.meta()<-` (setting new sample meta) to intelligently update `snpRdata()` objects by removing only calculated statistics and summary tabulations that applied to any changed facets instead
of simply re-importing the entire dataset as before. This should *substantially* speed up this function for large
data sets.
* Added the `smart_PCA` option to `plot_clusters()`. This will use Patterson et al. (2006)'s methods (with Price et al. (2006)'s allele frequency estimation) for centering and scaling genotypic data for PCA/tSNE/umap construction. This generally doesn't change much unless there is a lot of missing data, since this approach avoids imputation.
* Incorporated update of `sequoia V 2.5.3` which adds 'Year.last' - a cutoff for an individuals reproductive window into `format_snps()`. Returned `sequoia` to 'suggests'. Incorporated `sequoia` function `GetMaybeRel` in the `run_sequoia()` wrapper.
* Changed `plot_structure_map()` to take additional `ggplot2` layers directly and plot them prior to the pie charts instead of taking `sf` objects and trying to guess what the user wanted to do with them. This makes things considerably more flexible and makes it much easier to do things like plot precipitation/etc under the pie charts, although it means the user needs to be a bit more savy. Updated documentation to reflect.
* Added `nsnps` argument to `calc_ne()` to do automatic subsetting to run with less SNPs while still merging results into the original dataset. Note that this isn't terribly quick at the moment since it passes to the still somewhat inefficient subset operator `[`. With reasonably small numbers of SNPs, like what is usually suggested for LDNe, it should be fine.
* Added a few aliases to `get.snpR.stats()` for IBD ("ibd") and Tajima's D ("tsd", "d") to make it easier to fetch the correct values. Also adjusted the requested statistics to send everything to lower case, so things like "LD", "D", or "He" will still work.

## Bug fixes
* Fixed a weird bug where metadata class conversion during statistic calculation could result in merge errors.
* Fixed an un-informative error message when nothing remains after merging in `merge_snpRdata()`.
* Fixed an issue where having no matching column names during merging in `merge_snpRdata()` would produce weird results.
* Fixed an issue where having matching but not merged by column names during merging in `merge_snpRdata()` would produce an error.
* Fixed a niche bug in small data where attempting to build `geno.table` data during `snpRdata` facet tabulation when a facet level has no non-missing genotypes would error.
* Fixed a bug where running the `.base` facet *alongside* other facets with `calc_pairwise_ld()` with the `CLD` option
would cause an error during merging.
* Fixed a bug where running the `.base` facet *alongside* other facets with different snp level facets during `calc_smoothed_averages()` would cause errors.
* Fixed a bug where running `calc_pairwise_ld()` would fail to return a proximity table if there were `NA` values in the sample metadata.
* Fixed a bug where SFS construction by `calc_sfs()` and other SFS functions would error with some but not all data sets due to issues when adding `anc` and `ref` columns without using `snp.meta(x)<-`. Existing tests didn't catch this because it didn't occur with the `stickSNPs` test data or other test data sets.
* Fixed a bug where `plot_diagnostic()` would fail if run with a new facet and the `maf` plot option but not the `fis` plot option.
* Fixed a bug where an inappropriate error would be returned if `plot_structure()` was run with `facet = ".base"`, which should be treated like `facet = NULL`.


# snpR 1.2.6.1 -- hotfix
### bug fixes
* Fixed a bug where `calc_smoothed_averages()` and `calc_tajimas_d()` would do a step 100 times larger than expected if the default was used! Thus the hotfix.



# snpR 1.2.6

## Features

### Minor
* Added he/ho plotting to `plot_diagnostic()` and added manual control of which plots to generate. SFS skipped by default. Improved documentation and testing a bit.
* Changed the behavior of `calc_pairwise_fst()` to take the ratio of average variance components instead of the average of ratios following the recommendations of [Bahtia et al. 2013](https://doi.org/10.1101/gr.154831.113).
* To aid in calculating weighted mean $F_{ST}$ values for cases where the dataset is too big to run at once with snpR (as may be the case for large WGS data sets, for example), the `keep_components` argument was added to `calc_pairwise_fst()` to return the "a", "b", and "c" variance components for each locus for later processing. Brief instructions were added to the documentation for `calc_pairwise_fst()` to explain this process.


# snpR 1.2.5

## Features

### Major
* Rework to allow some features to be used with microsatellite or other non-biallelic data types such as microhaplotypes is underway and will be enabled once basic filtering, diversity calculations, and plotting function support is finished.
* Added support for DAPC to `plot_clusters()` via interface to adegenet and some code pulled in from (the thankfully GPL-v3) ade4 package. Licence note adjusted to reflect ade4 code.
* Improved memory efficiency of `snpRdata` objects by removing some duplicated information. Old objects should still work fine, but new objects will be considerably smaller. This also improves `snpRdata` object
creation times!
* Added `merge_snpRdata()` to merge `snpRdata` objects using syntax equivalent to base R's `merge()` function. This can still be made more computationally efficient in the future by avoiding some internal summary tabulation.

### Minor
* Reworked `calc_smoothed_averages()` to be more memory efficient (but slightly slower) when working with large datasets. Added `triple_sigma` and `gaussian` arguments that determine if $\sigma$ is tripled to have windows with a full size of 6 x sigma and determine if gaussian smoothing is actually used. Added more info to `get.snpR.stats()` window returns.
* Minor rework to importing and facet creation to be more memory efficient (but slightly slower) when working with large datasets.
* `calc_tajimas_d()` will now also return the number of raw segregating sites per window (which was already internally calculated, since it is a part of Watterson's Theta, but not returned).
* Changed the defaults for smoothing and tajima's D to `2*sigma` (non-overlapping windows)
* Added the `mac` argument to `filter_snps()` to filter by minor allele count instead of minor allele frequency. Currently doesn't support faceting, which will probably be added later depending on user need. The `singletons` argument is now depreciated.
* Added the `hwe_excess_side` argument ot `filter_snps()`, enabling users to only remove SNPs out of HWE that have either het or hom excesses. They default behavior is still to do both.
* Added $ZF_{ST}$ and $F_{ST}/(1-F_{ST})$ to `calc_pairwise_fst()`.
* Added the explicit option to recode chromosome names as simple numeric values to `format_snps()` with the `plink` option. For some cases with many scaffolds/etc, this may be necessary. Before this was the default behavior, not an option. To account for this, also added checks to ensure that, if this option is not used, that no chromosome names start in numbers (leading numbers will be replaced with a character equivalent -- 0 -> A, 1 -> B, 9 -> I).
* Changed the second facets used in `plot_clusters()` to use different shapes instead of different fill/color combos as long as
there are less than 25 levels (the number of unique point shapes). This makes for substantially easier to interpret levels in plots! Which
facet gets shapes is controlled by the `shape_has_more_levels` argument.
* Added a `verbose` option to `check_duplicates()`. It's still a slow function and could use some work or parallelization.
* Added `lambda_gc_correction` arguments to `plot_qq()` and `plot_manhattan()` to generate $\lambda_{GC}$ genomic stratification measures and correct for them in plots (see [this paper](https://doi.org/10.1038/nrg2813)).
* Added an informational shout-out on launch directing users towards the issues page and here.

## Documentation
* Changed "sample metadata" to "SNP metadata" in the description of the `header_cols` argument to `import.snpR.data()`.
* Edited the warning message when importing structure data without setting header_cols properly to more explicitly reference that problem.
* Updated the `calc_tajimas_d()` examples to use the updated `get.snpR.stats()` syntax.
* Added warnings/messages to `calc_ne` when it is run without the `chr` argument or with more than 5,000 SNPs.
* Added a more informative error message to `calc_ne` that shows if no output files are generated.
* Added a note to the documentation for `calc_prop_poly` to note that `calc_tajimas_d` will also calculate the number of segragating sites (and the number of snps) in a window, which can be used to easily calculate the prop_poly per window.
* Extended the error message when importing finds a bad genotype format to suggest removing SNP meta data from your genotypes.

## Bug fixes
* Changed `import.snpR.data()` to return an error if the genotypic dimensions are not correct for the provided SNP and sample meta data instead of proceeding and returning a bogus result.
* Fixed a bug where not using a pop facet with `plot_structure()` with the `method = "structure"` option would cause a bug due to an incorrectly set "LOCISPOP" flag (which is checked even if not using the locprior option).
* Fixed a bug where PSIXct data in the sample meta data could cause issues with `calc_hs()`.
* Fixed a bug where vcf files with missing data would occasionally get read in as "."s instead of NAs by `vcfR`, which meant that they weren't being properly accounted for as missing data by snpR. Most functions actually still work fine, but a few, like `plot_structure()` were unhappy.
* Added a check to tab-delimited input for a column of NAs at the end (automatically remove if found). Common when importing ANGSD outputs.
* Added `normalizePath()` to calc_ne to normalize the neestimator path.
* Added automatic facet checking for bad characters when doing facet operations. This is a bit slower (when there are *many* SNPs or individuals), but will return an informative error instead of an uninformative one and so improves usability.
* Fixed bugs that would occur when NeEstimator is one with only one pcrit and when only one pcrit + one population (both during parsing).
* Fixed a bug in the `coan` option for NeEstimator.
* Fixed a bug where asking to filter by the `.base` facet with `filter_snps()` would cause an error and a spurious warning.
* Standardized the defaults of both `facet` arguments to `filter_snps()` to use `NULL` as a default.

# snpR 1.2.4 -- hot-fix

## Features

### Minor
* Added the `gradient_colors` argument to `plot_pairwise_ld_heatmap()` to allow for custom sets of colors for the scale.
* Updated the readme for github with a function table of contents.
* Added another sanity check to `calc_ne()` to catch full file paths given to `outfile`.

## Bug fixes
* Fixed a nasty bug during imputation of sn data where the wrong maf was found. Shouldn't directly bias clustering, but good to fix asap.
* Swapped the `geom` used in `plot_pairwise_ld_heatmap()` to `geom_bin2d()` to prevent points from disappearing if too many loci are plotted on a chr.
* Fixed `read_structure()` (and `import.snpR.data()`) to never assume sample names, just loci names (which is the standard) if noted.
* Added better error checking for invalid NeEstimator exe files to `calc_ne()`.
* Fixed a bug where vcf files with missing data would occasionally get read in as "."s instead of NAs by `vcfR`, which meant that they weren't being properly accounted for as missing data by snpR. Most functions actually still work fine, but a few, like `plot_structure()` were unhappy.


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
