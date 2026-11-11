# Calculate Isolation by Distance

Calculates Isolation by Distance (IBD) for snpRdata objects by comparing
the genetic distance between samples or sets of samples to the
geographic distances between samples or sets of samples. IBD calculated
via a mantel test.

## Usage

``` r
calc_isolation_by_distance(
  x,
  facets = NULL,
  x_y = c("x", "y"),
  genetic_distance_method = "Edwards",
  interpolate = "bernoulli",
  ...
)
```

## Arguments

- x:

  snpRdata object

- facets:

  character, default NULL. Facets over which to calculate IBD, as
  described in
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).

- x_y:

  character, default c("x", "y"). Names of the columns containing
  geographic coordinates where samples were collected. There is no need
  to specify projection formats. The first should be longitude, the
  second should be latitude.

- genetic_distance_method:

  character, default "Edwards". The genetic distance method to use, see
  [`calc_genetic_distances`](https://hemstrow.github.io/snpR/reference/calc_genetic_distances.md).

- interpolate:

  Character or FALSE, default "bernoulli". If transforming to "sn"
  format, notes the interpolation method to be used to fill missing
  data. Options are "bernoulli", "af", "iPCA", or FALSE. See
  [`format_snps`](https://hemstrow.github.io/snpR/reference/format_snps.md)
  for details.

- ...:

  Additional arguments passed to
  [`mantel.randtest`](https://adeverse.github.io/ade4/reference/mantel.randtest.html)

## Value

A snpRdata object with both geographic distance and IBD results merged
into existing data.

## Details

Genetic distance is calculated via
[`calc_genetic_distances`](https://hemstrow.github.io/snpR/reference/calc_genetic_distances.md).
Geographic distances are taken as-is for individual-individual
comparisons and by finding the geographic mean of a group of samples
when sample level facets are requested using the methods described by
[`geomean`](https://rdrr.io/pkg/geosphere/man/geomean.html). Note that
this means that if many samples were collected from the same location,
and those samples compose a single level of a facet, the mean sampling
location will be that single location.

IBD is calculated by comparing geographic and genetic distances using a
mantel test via
[`mantel.randtest`](https://adeverse.github.io/ade4/reference/mantel.randtest.html).

Note that geographic distance objects are also included in the returned
values.

## See also

[`mantel.randtest`](https://adeverse.github.io/ade4/reference/mantel.randtest.html),
[`calc_genetic_distances`](https://hemstrow.github.io/snpR/reference/calc_genetic_distances.md)

## Author

William Hemstrom

## Examples

``` r
# calculate ibd for several different facets
if (FALSE) { # \dontrun{
y <- stickSNPs
sample.meta(y) <- cbind(sample.meta(y), x = rnorm(ncol(y)), y = rnorm(ncol(y)))
y <-calc_isolation_by_distance(y, facets = c(".base", "pop", "pop.chr","pop.chr.fam"))
res <- get.snpR.stats(y, "pop.chr", "ibd") # fetch results
res
plot(res$chr.pop$groupV$Edwards) # plot perms vs observed
} # }
```
