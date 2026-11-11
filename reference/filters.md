# Report on filters used on a snpRdata object.

snpR will automatically track the methods used on a snpRdata object.
Using `filters` on that object will provide details on the the filters
used.

## Usage

``` r
filters(x)
```

## Arguments

- x:

  snpRdata object

## Value

Prints easily human readable results to the console and also returns a
more machine readable string.

## Details

Printed outputs contain the filters used in the order they were applied
alongside the filtering stringency and any facets applied over, if
applicable. Note that some output formats from `format_snps`, like the
`vcf` format, will notes on the filters used as well.

## Author

William Hemstrom

## Examples

``` r
# filter the data
x <- filter_snps(stickSNPs, 
                 min_ind = 0.75, 
                 min_loci = 0.75,
                 maf = 0.1,
                 maf_facets = "pop")
#> Adding missing facets...
#> Initializing...
#> Filtering loci. Starting loci: 100 
#> Filtering non-biallelic loci...
#>  0 bad loci
#> Filtering non_polymorphic loci...
#>  0 bad loci
#> Filtering loci sequenced in few individuals...
#>  9 bad loci
#> Filtering low minor allele frequencies by facet.
#>  6 bad loci
#>  Ending loci: 85 
#> Filtering individuals. Starting individuals: 100 
#> Filtering out individuals sequenced in few kept loci...
#> Re-calculating and adding facets.
#>  Ending individuals: 95 
#> Re-filtering loci...
#> Filtering non_polymorphic loci...
#>  0 bad loci
#>  Final loci count: 85 

# fetch the filters used
filters(x)
#> Filter: bi-allelic 
#> 
#> Filter: bi-allelic 
#> 
#> Filter: non-polymorphic 
#> 
#> Filter: poorly sequenced loci--min_ind 
#>  Stringency: 0.75 
#>  Facet: .base 
#> 
#> Filter: Minor Allele Frequency--maf 
#>  Stringency: 0.1 
#>  Facet: pop 
#> 
#> Filter: poorly sequenced individuals--min_loci 
#>  Stringency: 0.75 
#> 
#> Filter: non-polymorphic 
#> 
#> [1] "bi-allelic;bi-allelic;non-polymorphic;min_ind=0.75,facet=.base;maf=0.1,facet=pop;min_loci=0.75;non-polymorphic"
```
