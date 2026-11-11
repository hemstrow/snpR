# Calculate p-values from bootstrapped distributions.

`calc_p_from_bootstraps` finds p-values for observed *smoothed window*
statistics from bootstrapped distributions, such as produced by
[`do_bootstraps`](https://hemstrow.github.io/snpR/reference/do_bootstraps.md).

## Usage

``` r
calc_p_from_bootstraps(
  x,
  facets = "all",
  statistics = "all",
  alt = "two-sided",
  par = FALSE,
  fwe_method = "BY",
  fwe_case = c("by_facet", "overall")
)
```

## Arguments

- x:

  snpRdata object.

- facets:

  character, default "all". Facets to use.

- statistics:

  character, default "all". Vector naming the statistics to calculate
  p-values for. By default calculates p-values for all possible stats.

- alt:

  character, default "two-sided". Specifies the alternative hypothesis
  to be used. Options:

  - "less": probability that a bootstrapped value is as small or smaller
    than observed.

  - "greater": probability that a bootstrapped value is as large or
    larger than observed.

  - "two-sided": probability that a bootstrapped value is as or more
    extreme than observed.

- par:

  numeric or FALSE, default FALSE. If numeric, the number of cores to
  use for parallel processing.

- fwe_method:

  character, default c("bonferroni", "holm", "BH", "BY"). Type of
  Family-Wise Error correction (multiple testing correction) to use. For
  details and options, see
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

- fwe_case:

  character, default c("by_facet", "by_subfacet", "overall"). How should
  Family-Wise Error correction (multiple testing correction) be applied?

  - "by_facet": Each facet supplied (such as pop or pop.fam) is treated
    as a set of tests.

  - "by_subfacet": Each level of each subfacet is treated as a separate
    set of tests.

  - "overall": All tests are treated as a set.

## Value

snpRdata object, with p-values merged into the stats or pairwise.stats
sockets.

## Details

Calculates p-values for smoothed values of a statistic based upon a
bootstrapped null distribution of that statistic using an empirical
continuous distribution function.

p-values can be generated for specific snp or sample metadata categories
using the facets argument, as described in
[`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).
Only facets for which bootstrap data and raw statistical data have both
been calculated will be run. "all" and NULL follow the typical facet
rules.

Likewise, p-values can be generated for specific statistics using the
statistics argument. Only statistics for which bootstrap data and raw
statistical data have both been calculated will be run. By default, all
stats for which a bootstrap null distribution has been generated will be
run.

## See also

ecdf

## Author

William Hemstrom

## Examples

``` r
if (FALSE) { # \dontrun{
# add statistics and generate bootstraps
x <- calc_basic_snp_stats(stickSNPs, c("chr.pop"), sigma = 200, step = 150)
x <- do_bootstraps(x, facets = c("chr.pop"), boots = 1000, sigma = 200, step = 150)
x <- calc_p_from_bootstraps(x)
get.snpR.stats(x, "chr.pop", "single.window") # pi, ho, etc
get.snpR.stats(x, "chr.pop", "pairwise.window") # fst
} # }
```
