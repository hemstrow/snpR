# Calculate the number of segregating sites.

Calculates the number segregating loci for each facet level with or
without rarefaction to account for differing sample size and missing
data levels across populations.

## Usage

``` r
calc_seg_sites(x, facets = NULL, rarefaction = TRUE, g = 0)
```

## Arguments

- x:

  snpRdata object

- facets:

  facets over which to calculate the number of segregating sites. See
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md)
  for details.

- rarefaction:

  logical, default TRUE. Should the number of segregating sites be
  estimated via rarefaction? See details.

- g:

  numeric, default 0. If doing rarefaction, controls the number of
  *genotypes* to rarefact to. If 0, this will rarefact to the smallest
  sample size per locus. If g \< 0, this will rarefact to to the
  smallest sample size per locus minus the absolute value of g. If
  positive, this will rarefact to g, and any loci where the smallest
  sample size is less than g will be dropped from the calculation.

## Value

A snpRdata object with seg_sites merged into the weighted.means slot.

## Details

Rarefaction is done by determining the probability of drawing at least
one copy of each allele given a standardized sample size, *g*, given the
observed genotype frequencies at each locus. Using the observed genotype
frequencies ensures that this is not biased by deviations from
Hardy-Weinburg equilibrium. This is then summed across all loci to get
the expected number of segregating sites for each population/subfacet.

Note that *g* will vary across loci due to differences in sequencing
coverage at those loci, equal to the smallest number of genotypes
sequenced in any population at that locus minus one.

Note no sample-specific facet is requested, rarefaction will not be
used.

## Author

William Hemstrom

## Examples

``` r
# base facet
x <- calc_seg_sites(stickSNPs, rarefaction = FALSE)
get.snpR.stats(x, stats = "seg_sites")$weighted.means
#>   facet subfacet snp.facet snp.subfacet seg_sites seg_sites_var
#> 1 .base    .base     .base        .base       100            NA

# multiple facets
x <- calc_seg_sites(stickSNPs, c("pop", "fam"))
get.snpR.stats(x, c("pop", "fam"), 
               stats = "seg_sites")$weighted.means
#>    facet subfacet snp.facet snp.subfacet seg_sites seg_sites_var
#> 1    fam        A     .base        .base  95.72737     0.2507511
#> 2    fam        B     .base        .base  98.50845     0.3610588
#> 3    fam        C     .base        .base  94.45007     0.4693410
#> 4    fam        D     .base        .base  97.60164     0.3590730
#> 5    pop      ASP     .base        .base  82.52969     5.1091053
#> 6    pop      CLF     .base        .base  78.00000     0.0000000
#> 7    pop      OPL     .base        .base  81.50088     7.4033292
#> 8    pop      PAL     .base        .base  47.60302     3.4709544
#> 9    pop      SMR     .base        .base  81.45156     8.3645599
#> 10   pop      UPD     .base        .base  84.86798     5.8373148
```
