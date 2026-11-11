# Tajima's D from SNP data.

`Tajimas_D` calculates Tajima's theta/pi, Watterson's theta, and
Tajima's D over a sliding window.

## Usage

``` r
calc_tajimas_d(
  x,
  facets = NULL,
  sigma = NULL,
  step = 2 * sigma,
  par = FALSE,
  triple_sigma = FALSE,
  global = FALSE,
  verbose = FALSE
)
```

## Arguments

- x:

  snpRdata. Input SNP data.

- facets:

  character. Categorical metadata variables by which to break up
  analysis. See
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md)
  for more details. If no snp level facets are provided, the calculated
  Tajima's D will be for the entire genome.

- sigma:

  numeric, default NULL. Sliding window size, in kilobases. Each window
  will include all SNPs within 3\*sigma or sigma kilobases depending on
  the `triple_sigma` argument. If either sigma or step are NULL, the
  entire snp subfacet will be done at once (for example, the whole
  chromosome).

- step:

  numeric or NULL, default `2*sigma` (non-overlapping windows). Number
  of bases to move between each window, in kilobases. If either sigma or
  step are NULL, the entire snp subfacet will be done at once (for
  example, the whole chromosome).

- par:

  numeric or FALSE, default FALSE. If numeric, the number of cores to
  use for parallel processing.

- triple_sigma:

  logical, default FALSE If TRUE, sigma will be tripled to create
  windows of 6\*sigma total.

- global:

  logical, default FALSE. If TRUE, all window parameters will be ignored
  and the global Tajima's D across all sites will instead be calculated.
  In this instance, `global_x` values will be merged into the
  `weighted.means` slot instead of weighted mean values and no values
  will be merged into the `window.stats` slot.

- verbose:

  logical, default FALSE. If TRUE progress will be printed to the
  console.

## Value

snpRdata object, with Watterson's Theta, Tajima's Theta, and Tajima's D
for each window merged in to the window.stats slot.

## Details

Tajima's D compares estimates of theta based on either the number of
observed pairwise differences (Tajima's theta) and the number of
substitutions vs expected total tree length (Watterson's Theta). Since
low frequency minor variants contribute to these statistics and they
rely on the ratio of the number of variants vs the number of sequenced
non-polymorphic sites, this function should only be run on data that is
*unfiltered* aside from the removal of poorly sequenced bases, etc.

The data can be broken up categorically by either SNP and/or sample
metadata, as described in
[`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).

## References

Tajima, F. (1989). *Genetics*

## Author

William Hemstrom

## Examples

``` r
# slow, so not run
if (FALSE) { # \dontrun{
# broken by population, windows across linkage group
x <- calc_tajimas_d(stickSNPs, facets = "chr.pop", sigma = 200, step = 50)
get.snpR.stats(x, "chr.pop", "tajimas_d")

# the entire population at once, note that sigma and step are NULL and
# no chromosome/linkage group/scaffold/etc set.
# this will calculate overall tajima's D without a window for each population.
x <- calc_tajimas_d(stickSNPs, facets = "pop")
get.snpR.stats(x, "pop", "tajimas_d")

# for the overall dataset, note that sigma and step are NULL
# this will calculate overall tajima's D for each chr/pop
x <- calc_tajimas_d(stickSNPs, facets = "chr.pop")
get.snpR.stats(x, "pop.chr", "tajimas_d")
} # }
```
