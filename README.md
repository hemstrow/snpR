
<!-- README.md is generated from README.Rmd. Please edit that file -->

# snpR

<!-- badges: start -->

[![packageversion](https://img.shields.io/badge/Package%20version-1.2.4-orange.svg?style=flat-square)](commits/master)
[![CRAN
status](https://www.r-pkg.org/badges/version/snpR)](https://CRAN.R-project.org/package=snpR)
[![R-CMD-check](https://github.com/hemstrow/snpR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hemstrow/snpR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

snpR is an R package for analyzing call Single Nucleotide Polymorphism
(SNP) genotypes containing most basic stats including pairwise LD,
gaussian sliding window analysis tools, plotting options, clustering
analysis, colony interface, Ne estimation, formatting, filtering, and
more!

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

-   `import.snpR.data()`: generic read function, takes many file types
    by extension or R `data.frames()`.
-   Wrappers for specific file types:
    -   `read_structure()`: Reads STRUCTURE “.str” files.
    -   `read_vcf()`: Reads VCF “.vcf” files.
    -   `read_FSTAT()`: Reads FSTAT “.fstat” files.
    -   `read_ms()`: Reads ms “.ms” files.
    -   `read_delimited_snps()`: Reads tab-delimited “NN” or “0000”
        data.
    -   `read_genepop()`: Reads genepop “.genepop” files.
    -   `read_plink()`: Reads PLINK! “.bed”, “.fam”, and “.bim” files.
    -   `convert_genlight()`: Converts `adegenet` `genlight` class
        objects.
    -   `convert_genind()`: Converts `adegenet` `genind` class objects.
    -   `convert_vcfR()`: Converts `vcfR` class objects.

### Utility:

-   `filter_snps()`: Filter data.
-   `format_snps()`: Format data into other export formats.
-   `summarize_facets()`: Summarized available facets.
-   `citations()`: Fetch citations for all methods used in calculations
    for a specific `snpRdata` object.
-   `check_duplicates()`: Check data for potentially duplicated samples.
-   `gap_snps()`: Select a SNP every *n* bases (simple physical LD
    filtering).

### Object Access and Manipulation:

-   Dimensions:
    -   `nsnps()` and `nrow()`: Get the number of SNPs in an object.
    -   `nsamps()` and `ncol()`: Get the number of samples in an object.
    -   `dim()`: Get number of SNPs and samples in an object.
-   Access:
    -   `get.snpR.stats()`: Fetch any calculated statistics from an
        object.
    -   `genotypes()`: Fetch genotypes.
    -   `sample.meta()`: Fetch (or reasign with `<-`) sample metadata.
    -   `snp.meta()`: Fetch (or reasign with `<-`) SNP metadata.
-   Subetting:
    -   `[`: The usual bracket operator. Subset by SNP or sample index,
        or by facet.
    -   `subset_snpR_data()`: Wrapper for the bracket operator.

### Statistics:

-   Basic statistics
    -   `calc_pi()`: Nucleotide diversity.
    -   `calc_ho()`: Observed heterozygosity.
    -   `calc_he()`: Expected heterosygosity.
    -   `calc_hwe()`: HWE.
    -   `calc_hs()`: Standardized individual heterozygosity.
        -   `calc_het_hom_ratios()`: Alternative, raw
            heterozygote/homozygote ratios within individuals.
    -   `calc_ne()`: Effective population size.
    -   `calc_prop_poly()`: The proportion of polymorphic loci.
    -   `calc_maf()`: Minor allele frequencies, calculated automatically
        when any facet operations are performed.
    -   `calc_private()`: Private alleles across facet levels.
    -   `calc_genetic_distances()`: Genetic distances between
        individuals
    -   `calc_fis()`: FIS
    -   `calc_pairwise_fst()`: Pairwise FST between facet levels.
    -   `calc_pairwise_ld()`: Pairwise LD between SNPs.
    -   `calc_abba_baba()`: ABBA/BABA tests.
-   Association:
    -   `calc_association()`: Association testing against a phenotype.
    -   `run_random_forest()`: Run a random forest
        prediction/association test against a phenotype.
    -   `run_random_forest()`: Run genomic prediction against a
        phenotype.
        -   `cross_validate_genomic_prediction()`: Bare-bones
            cross-validation for genomic predictions.
-   Site-frequency Spectra:
    -   `calc_sfs()`: Generate a 1 or 2d site frequency spectra.
        -   `make_tree()`: Wrapper function that uses an external `dadi`
            formatted file to generate an sfs.
    -   `calc_directionality()`: Peter and Slatkin’s directionality
        index.
-   Other:
    -   `calc_isolation_by_distance()`: Run an IBD mantel test.
    -   `calc_tree()`: Generate a tree based on individual or
        facet-level relatedness.
    -   `tabulate_allele_frequency_matrix()`: Generate an allele
        frequency matrix.

### Windows:

-   `calc_smoothed_averages()`: Core function to do sliding window
    analysis using a gaussian smoothing kernal.
-   `calc_tajimas_d()`: Tajima’s D across sliding windows.
-   Bootstrapping:
    -   `do_bootstraps()`: Core function to generate bootstrapped
        significance values for smoothed windows (elevation or reduction
        vs genomic background).
        -   `calc_p_from_bootstraps()`: Calc p-values from bootstraps.
            Run automatically by `do_boostraps()`.

### Plotting:

-   `plot_clusters()`: PCA, UMAP, and tSNE plots.
-   `plot_structure()`: Run STRUCTURE or several alternatives OR read in
    existing “q” files and generate plots.
-   `plot_structure_map()`: Plots `plot_structure()` or parsed in q file
    results on a map given coordinates for populations.
-   `plot_diagnostic()`: A suite of useful diagnostic plots.
-   `plot_manhattan()`: Manhattan plots from calculated statistics or a
    `data.frame()`. Excellent for visualizing most statistics
    genome-wide (not just association tests!)
-   `plot_qq()`: Quantile-quantile (qq) plots from calculated
    association test results.
-   `plot_pairwise_fst_heatmap()`: Heatmap of FST scores between facet
    levels.
-   `plot_pairwise_ld_heatmap()`: Heatmap of LD scores between SNPs.

### Parentage:

-   Colony:
    -   `run_colony()`: All-in-one function to make a colony import
        file, run colony, and parse results.
    -   `write_colony()`, `call_colony()`, `parse_colony()`: Write input
        files, call colony, and parse results as seperate functions.
-   Sequoia:
    -   `run_sequoia()`: Run a basic parentage assessment with the
        `sequoia` package.

## Example

snpR is focused on ease-of-use. Primarily, it achieves this via the use
of , which describe sample or SNP metadata. snpR is built to
automatically split up analysis by facet. For example, calculating
observed heterozygosity for each population or family, or for each
population/family combination is easy!

``` r
library(snpR)
#> Loading required package: data.table
#> Loading required package: foreach
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

For a full introduction, check the snpR_introduction vignette.

``` r
# remotes::install_github("hemstrow/snpR", build_vignettes = T, build_opts = c("--no-resave-data", "--no-manual"))

vignette("snpR_introduction")
```
