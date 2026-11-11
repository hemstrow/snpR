# snpR

snpR is an R package for analyzing call Single Nucleotide Polymorphism
(SNP) genotypes containing most basic stats including pairwise LD,
gaussian sliding window analysis tools, plotting options, clustering
analysis, colony interface, Ne estimation, formatting, filtering, and
more! It is built primarily to be user-friendly and handle many levels
of SNP and sample metadata without the need for complicated file or
object management. Please see the example below!

## Installation

snpR can be installed from [GitHub](https://github.com/hemstrow/snpR):

``` r

# install.packages("remotes")
remotes::install_github("hemstrow/snpR")
```

To install the vignettes as well (recommended for new users), instead
use:

``` r

remotes::install_github("hemstrow/snpR", build_vignettes = T) # linux
remotes::install_github("hemstrow/snpR", build_vignettes = T, build_opts = c("--no-resave-data", "--no-manual")) # windows
```

If you wish to try out the latest features or bug fixes, the dev version
can be installed from [GitHub](https://github.com/hemstrow/snpR) as
well:

``` r

# install.packages("remotes")
remotes::install_github("hemstrow/snpR", ref = "dev")
```

A CRAN version should be available soon.

## Function Table of Contents

### Import:

- [`import.snpR.data()`](https://hemstrow.github.io/snpR/reference/import.snpR.data.md):
  generic read function, takes many file types by extension or R
  `data.frames()`.
- Wrappers for specific file types:
  - [`read_structure()`](https://hemstrow.github.io/snpR/reference/snpR_import_wrappers.md):
    Reads STRUCTURE “.str” files.
  - [`read_vcf()`](https://hemstrow.github.io/snpR/reference/snpR_import_wrappers.md):
    Reads VCF “.vcf” files.
  - [`read_FSTAT()`](https://hemstrow.github.io/snpR/reference/snpR_import_wrappers.md):
    Reads FSTAT “.fstat” files.
  - [`read_ms()`](https://hemstrow.github.io/snpR/reference/snpR_import_wrappers.md):
    Reads ms “.ms” files.
  - [`read_delimited_snps()`](https://hemstrow.github.io/snpR/reference/snpR_import_wrappers.md):
    Reads tab-delimited “NN” or “0000” data.
  - [`read_genepop()`](https://hemstrow.github.io/snpR/reference/snpR_import_wrappers.md):
    Reads genepop “.genepop” files.
  - [`read_plink()`](https://hemstrow.github.io/snpR/reference/snpR_import_wrappers.md):
    Reads PLINK! “.bed”, “.fam”, and “.bim” files.
  - [`convert_genlight()`](https://hemstrow.github.io/snpR/reference/snpR_import_wrappers.md):
    Converts `adegenet` `genlight` class objects.
  - [`convert_genind()`](https://hemstrow.github.io/snpR/reference/snpR_import_wrappers.md):
    Converts `adegenet` `genind` class objects.
  - [`convert_vcfR()`](https://hemstrow.github.io/snpR/reference/snpR_import_wrappers.md):
    Converts `vcfR` class objects.

### Utility:

- [`filter_snps()`](https://hemstrow.github.io/snpR/reference/filter_snps.md):
  Filter data.
- [`format_snps()`](https://hemstrow.github.io/snpR/reference/format_snps.md):
  Format data into other export formats.
- [`summarize_facets()`](https://hemstrow.github.io/snpR/reference/summarize_facets.md):
  Summarized available facets.
- [`citations()`](https://hemstrow.github.io/snpR/reference/citations.md):
  Fetch citations for all methods used in calculations for a specific
  `snpRdata` object.
- [`check_duplicates()`](https://hemstrow.github.io/snpR/reference/check_duplicates.md):
  Check data for potentially duplicated samples.
- [`gap_snps()`](https://hemstrow.github.io/snpR/reference/gap_snps.md):
  Select a SNP every *n* bases (simple physical LD filtering).
- [`vcf2beagle()`](https://hemstrow.github.io/snpR/reference/vcf_format_utils.md)
  and `vcf2PL()`: Utilities to convert vcf files to either beagle or
  pl/mgpl files. The latter can also optionally prepare inputs for the
  program
  [entropy](https://bitbucket.org/buerklelab/mixedploidy-entropy/src/master/).

### Object Access and Manipulation:

- Dimensions:
  - [`nsnps()`](https://hemstrow.github.io/snpR/reference/snpRdata_dims.md)
    and
    [`nrow()`](https://hemstrow.github.io/snpR/reference/snpRdata_dims.md):
    Get the number of SNPs in an object.
  - [`nsamps()`](https://hemstrow.github.io/snpR/reference/snpRdata_dims.md)
    and
    [`ncol()`](https://hemstrow.github.io/snpR/reference/snpRdata_dims.md):
    Get the number of samples in an object.
  - [`dim()`](https://hemstrow.github.io/snpR/reference/snpRdata_dims.md):
    Get number of SNPs and samples in an object.
- Access:
  - [`get.snpR.stats()`](https://hemstrow.github.io/snpR/reference/get.snpR.stats.md):
    Fetch any calculated statistics from an object.
  - [`genotypes()`](https://hemstrow.github.io/snpR/reference/extract_snpRdata.md):
    Fetch genotypes.
  - [`sample.meta()`](https://hemstrow.github.io/snpR/reference/extract_snpRdata.md):
    Fetch (or reasign with `<-`) sample metadata.
  - [`snp.meta()`](https://hemstrow.github.io/snpR/reference/extract_snpRdata.md):
    Fetch (or reasign with `<-`) SNP metadata.
- Subetting:
  - `[`: The usual bracket operator. Subset by SNP or sample index, or
    by facet.
  - [`subset_snpR_data()`](https://hemstrow.github.io/snpR/reference/subset_snpRdata.md):
    Wrapper for the bracket operator.

### Statistics:

- Basic statistics
  - [`calc_pi()`](https://hemstrow.github.io/snpR/reference/calc_single_stats.md):
    Nucleotide diversity.
  - [`calc_ho()`](https://hemstrow.github.io/snpR/reference/calc_single_stats.md):
    Observed heterozygosity.
  - [`calc_he()`](https://hemstrow.github.io/snpR/reference/calc_single_stats.md):
    Expected heterosygosity.
  - [`calc_hwe()`](https://hemstrow.github.io/snpR/reference/calc_single_stats.md):
    Hardy-Weinburg equilibrium (HWE).
  - [`calc_hs()`](https://hemstrow.github.io/snpR/reference/individual_heterozygosity.md):
    Standardized individual heterozygosity.
    - `calc_het_hom_ratios()`: Alternative, raw heterozygote/homozygote
      ratios within individuals.
  - [`calc_ne()`](https://hemstrow.github.io/snpR/reference/calc_ne.md):
    Effective population size.
  - [`calc_prop_poly()`](https://hemstrow.github.io/snpR/reference/calc_prop_poly.md):
    The proportion of polymorphic loci.
  - [`calc_maf()`](https://hemstrow.github.io/snpR/reference/calc_single_stats.md):
    Minor allele frequencies, calculated automatically when any facet
    operations are performed.
  - [`calc_private()`](https://hemstrow.github.io/snpR/reference/calc_single_stats.md):
    Rarefaction-corrected detection of private alleles across facet
    levels.
  - [`calc_seg_sites()`](https://hemstrow.github.io/snpR/reference/calc_seg_sites.md):
    Rarefaction-corrected estimates of the number of segregating sites
    per facet level.
  - `calc_allelic_richness():`: Rarefaction-corrected estimates of
    allele counts per locus per facet level.
  - [`calc_genetic_distances()`](https://hemstrow.github.io/snpR/reference/calc_genetic_distances.md):
    Genetic distances between individuals.
  - [`calc_fis()`](https://hemstrow.github.io/snpR/reference/calc_fis.md):
    $`F_{IS}`$ (inbreeding coefficients).
  - [`calc_roh()`](https://hemstrow.github.io/snpR/reference/calc_roh.md):
    Run-of-Homozygosity detection and $`F_{ROH}`$ calculation.
  - [`calc_pairwise_fst()`](https://hemstrow.github.io/snpR/reference/calc_fst.md):
    Pairwise $`F_{ST}`$ between facet levels.
  - [`calc_global_fst()`](https://hemstrow.github.io/snpR/reference/calc_fst.md):
    Global $`F_{ST}`$ across facet levels.
  - [`calc_pairwise_ld()`](https://hemstrow.github.io/snpR/reference/calc_pairwise_ld.md):
    Pairwise LD between SNPs.
  - [`calc_abba_baba()`](https://hemstrow.github.io/snpR/reference/calc_abba_baba.md):
    ABBA/BABA tests.
  - [`calc_tajimas_d()`](https://hemstrow.github.io/snpR/reference/calc_tajimas_d.md):
    Tajima’s D globally (or across windows). Also calculates Watterson’s
    and Tajima’s $`\theta`$.
- Association:
  - [`calc_association()`](https://hemstrow.github.io/snpR/reference/calc_association.md):
    Association testing against a phenotype.
  - [`run_random_forest()`](https://hemstrow.github.io/snpR/reference/run_random_forest.md):
    Run a random forest prediction/association test against a phenotype.
  - [`run_random_forest()`](https://hemstrow.github.io/snpR/reference/run_random_forest.md):
    Run genomic prediction against a phenotype.
    - [`cross_validate_genomic_prediction()`](https://hemstrow.github.io/snpR/reference/cross_validate_genomic_prediction.md):
      Bare-bones cross-validation for genomic predictions.
- Site-frequency Spectra:
  - [`calc_sfs()`](https://hemstrow.github.io/snpR/reference/calc_sfs.md):
    Generate a 1 or 2d site frequency spectra.
    - `make_sfs()`: Wrapper function that uses an external `dadi`
      formatted file to generate an sfs.
  - [`calc_directionality()`](https://hemstrow.github.io/snpR/reference/calc_directionality.md):
    Peter and Slatkin’s directionality index.
  - [`calc_origin_of_expansion()`](https://hemstrow.github.io/snpR/reference/calc_origin_of_expansion.md):
    Estimate the origin point of a range expansion based on
    directionality.
- Other:
  - [`calc_isolation_by_distance()`](https://hemstrow.github.io/snpR/reference/calc_isolation_by_distance.md):
    Run an IBD mantel test.
  - [`calc_tree()`](https://hemstrow.github.io/snpR/reference/calc_tree.md):
    Generate a tree based on individual or facet-level relatedness.
  - [`tabulate_allele_frequency_matrix()`](https://hemstrow.github.io/snpR/reference/tabulate_allele_frequency_matrix.md):
    Generate an allele frequency matrix.

### Windows:

- [`calc_smoothed_averages()`](https://hemstrow.github.io/snpR/reference/calc_smoothed_averages.md):
  Core function to do sliding window analysis using a gaussian smoothing
  kernal.
- [`calc_tajimas_d()`](https://hemstrow.github.io/snpR/reference/calc_tajimas_d.md):
  Tajima’s D across sliding windows (or globally). Also calculates
  Watterson’s and Tajima’s $`\theta`$.
- Bootstrapping:
  - [`do_bootstraps()`](https://hemstrow.github.io/snpR/reference/do_bootstraps.md):
    Core function to generate bootstrapped significance values for
    smoothed windows (elevation or reduction vs genomic background).
    - [`calc_p_from_bootstraps()`](https://hemstrow.github.io/snpR/reference/calc_p_from_bootstraps.md):
      Calc p-values from bootstraps. Run automatically by
      `do_boostraps()`.

### Plotting:

- [`plot_clusters()`](https://hemstrow.github.io/snpR/reference/plot_clusters.md):
  PCA, UMAP, and tSNE plots.
- [`plot_structure()`](https://hemstrow.github.io/snpR/reference/plot_structure.md):
  Run STRUCTURE or several alternatives OR read in existing “q” files
  and generate plots.
- [`plot_structure_map()`](https://hemstrow.github.io/snpR/reference/plot_structure_map.md):
  Plots
  [`plot_structure()`](https://hemstrow.github.io/snpR/reference/plot_structure.md)
  or parsed in q file results on a map given coordinates for
  populations.
- [`plot_diagnostic()`](https://hemstrow.github.io/snpR/reference/plot_diagnostic.md):
  A suite of useful diagnostic plots.
- [`plot_manhattan()`](https://hemstrow.github.io/snpR/reference/plot_manhattan.md):
  Manhattan plots from calculated statistics or a
  [`data.frame()`](https://rdrr.io/r/base/data.frame.html). Excellent
  for visualizing most statistics genome-wide (not just association
  tests!)
- [`plot_qq()`](https://hemstrow.github.io/snpR/reference/plot_qq.md):
  Quantile-quantile (qq) plots from calculated association test results.
- [`plot_pairwise_fst_heatmap()`](https://hemstrow.github.io/snpR/reference/plot_pairwise_fst_heatmap.md):
  Heatmap of FST scores between facet levels.
- [`plot_pairwise_ld_heatmap()`](https://hemstrow.github.io/snpR/reference/plot_pairwise_ld_heatmap.md):
  Heatmap of LD scores between SNPs.
- [`plot_sfs()`](https://hemstrow.github.io/snpR/reference/plot_sfs.md):
  Plots site frequency spectra either from pre-consturcted SFS or
  directly from `snpRdata` object.
- [`plot_roh()`](https://hemstrow.github.io/snpR/reference/plot_roh.md):
  Plots ROHs across the genome.

### Parentage:

- Colony:
  - [`run_colony()`](https://hemstrow.github.io/snpR/reference/colony_interface.md):
    All-in-one function to make a colony import file, run colony, and
    parse results.
  - `write_colony()`,
    [`call_colony()`](https://hemstrow.github.io/snpR/reference/colony_interface.md),
    [`parse_colony()`](https://hemstrow.github.io/snpR/reference/colony_interface.md):
    Write input files, call colony, and parse results as seperate
    functions.
- Sequoia:
  - [`run_sequoia()`](https://hemstrow.github.io/snpR/reference/run_sequoia.md):
    Run a basic parentage assessment with the `sequoia` package.

## Example

snpR is focused on ease-of-use. Primarily, it achieves this via the use
of *facets*, which describe sample or SNP metadata. snpR is built to
automatically split up analysis by facet. For example, calculating
observed heterozygosity for each population or family, or for each
population/family combination is easy!

### Statistic calculation with facets

``` r

library(snpR)
## basic example code

x <- calc_ho(stickSNPs, facets = c("pop")) # split by pop (stickSNPs is an example dataset included in snpR)
x <- calc_ho(x, facets = c("fam")) # split by family
x <- calc_ho(x, facets = c("pop.fam")) # split by combinations of family and pop
```

snpR also facilitates ease-of-use by being *overwrite safe*. As above,
new analyses are added to an existing object. Results can be fetched
using the get.snpR.stats handler.

``` r

res <- get.snpR.stats(x, facets = "pop", stats = "ho")
```

Functions in snpR are consistently named: functions that calculate
statistics are prefixed `calc_`, functions that do plots are prefixed
`plot_`, and functions that run external tools (like COLONY), are named
`run_`. Typing `snpR::calc` into the console on Rstudio will bring up a
helpful list of all of the statistical functions!

### Basic plotting

`snpR` provides a suite of plotting tools, all of which are named
`plot_...()`. For example, a PCA can be plotted using
[`plot_clusters()`](https://hemstrow.github.io/snpR/reference/plot_clusters.md):

``` r

plot_clusters(x, "pop", plot_type = "pca")
```

### Citation tools

`snpR` automatically tracks citations for all of the methods used on a
`snpRdata` object and can provide them or write a `.bib` file on request
using the
[`citations()`](https://hemstrow.github.io/snpR/reference/citations.md)
function:

``` r

citations(x,
          outbib = FALSE, # if TRUE, writes a .bib file for all methods
          return_bib = FALSE # if TRUE, returns a list of .bib entries
          )
```

## Vignette

For a full introduction, check the snpR_introduction vignette.

``` r

# remotes::install_github("hemstrow/snpR", build_vignettes = T, build_opts = c("--no-resave-data", "--no-manual"))

vignette("snpR_introduction")
```
