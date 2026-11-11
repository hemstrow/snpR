# Basic diagnostic plots

Create a suite of basic diagnostic plots (FIS density, 1D SFS, maf
density, PCA, % missing data per individual, and he vs ho per SNP) to
describe the condition of the data in a snpRdata object.

## Usage

``` r
plot_diagnostic(
  x,
  facet = NULL,
  projection = floor(nsnps(x)/1.2),
  fold_sfs = TRUE,
  plots = c("fis", "maf", "pca", "missingness", "heho")
)
```

## Arguments

- x:

  snpRdata object

- facet:

  character, default NULL. Categorical metadata variables by which to
  break up plots. Note that only one facet is allowed here. Missingness
  and the PCA will have individuals colored by the given sample facet.
  See
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md)
  for more details.

- projection:

  integer, default floor(nsnps(x)/1.2). A sample size to project the SFS
  to, in *number of gene copies*. Sizes too large will result in a SFS
  containing few or no SNPs.

- fold_sfs:

  logical, default TRUE. Determines if the SFS should be folded or left
  polarized. If FALSE, snp metadata columns named "ref" and "anc"
  containing the identity of the derived and ancestral alleles,
  respectively, should be present for polarization to be meaningful.

- plots:

  character vector, default all possible plots except for SFS. Plot
  options: \\ fis: density of FIS scores for all loci within each facet
  level.

  - sfs: Site Frequency Spectra for the entire dataset.

  - maf:density of minor allele frequencies for all loci within each
    facet.

  - pca: Principal Component Analysis results for the given facet.

  - missingness: Proportion of missing alleles across each individual
    withing each facet.

  - heho: expected vs. observed heterozygosity for each locus within
    each facet. Very high expected or observed heterozygosities for many
    loci can indicate genotyping issues.

## Value

A named list of diagnostic ggplot2 plots.

## Author

William Hemstrom

## Examples

``` r
if (FALSE) { # \dontrun{
# missingness and pca colored by pop
plot_diagnostic(stickSNPs, "pop")
} # }
```
