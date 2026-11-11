# Filter down to one SNPs every *n* bases.

Selects one SNPs every *n* bases. This can be used to filter for SNPs in
tight LD, such as when only one SNP per RAD-tag is desired.

## Usage

``` r
gap_snps(x, facet = NULL, n)
```

## Arguments

- x:

  snpRdata object

- facet:

  character, default NULL. SNP facet specifying chromosomes or
  scaffolds. SNP positions will be independently considered depending on
  the facet level.

- n:

  Integer. Specifies the minimum distance between selected SNPs.

## Value

A snpRdata object containing the selected SNPs.

## Details

Note that this approach takes the first SNP every *n* bases, and so is
non-random.
[`filter_snps`](https://hemstrow.github.io/snpR/reference/filter_snps.md)
can be used beforehand to ensure that the selected SNPs are above a
desired quality threshold to ensure that poor SNPs are not selected over
loci with more robust genotyping data.

It is strongly recommended to provide a SNP facet to the facet argument.
This will define "chromosomes" for the purpose of selecting SNPs. If not
set, snpR will assume that *all SNPs are on the same chromosome*, which
may produce undesired results.

## Author

William Hemstrom

## Examples

``` r
# put large gaps inbetween the example data
gapped <- gap_snps(stickSNPs, "chr", 100000)
#> Selected  100 out of  100 SNPs.

# filter first to grab only very good SNPs
vgood <- filter_snps(stickSNPs, min_ind = .8)
#> Initializing...
#> Filtering loci. Starting loci: 100 
#> Filtering non-biallelic loci...
#>  0 bad loci
#> Filtering non_polymorphic loci...
#>  0 bad loci
#> Filtering loci sequenced in few individuals...
#>  16 bad loci
#>  Ending loci: 84 
#> Re-filtering loci...
#> Filtering non_polymorphic loci...
#>  0 bad loci
#>  Final loci count: 84 
gapped <- gap_snps(vgood, "chr", 100000)
#> Selected  84 out of  84 SNPs.
```
