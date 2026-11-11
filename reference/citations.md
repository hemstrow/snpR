# Gather citations for methods used with a snpRdata object

snpR will automatically track the methods used for calculations on a
snpRdata object. Using `citations` on that object will provide details
on the the methods used, and can optionally write a .bib bibtex
formatted library containing citations for these methods.

## Usage

``` r
citations(x, outbib = FALSE, return_bib = FALSE)
```

## Arguments

- x:

  snpRdata object

- outbib:

  character, default FALSE. An optional file path to which a .bib bibtex
  library containing all of the citations for the used methods will be
  written.

- return_bib:

  logical, default FALSE. If TRUE, returns a list containing the bib
  entries and their details. The bib entries are formatted according to
  [`readBib`](https://geobosh.github.io/rbibutils/reference/readBib.html).

## Value

If return_bib is TRUE, a list containing four parts:

- keys: A vector of bibtex keys for each method.

- stats: A vector of the stats used.

- details: A vector of details for each method.

- bib: A `bibentry` for each citation, see
  [`readBib`](https://geobosh.github.io/rbibutils/reference/readBib.html).

## Details

Printed outputs contain the statistic calculated, the in-line citation
for the method used, a bibtex key which corresponds to the .bib library
written if the outbib argument is used, and a quick note giving any
details.

## Author

William Hemstrom

## Examples

``` r
# calculate pi
x <- calc_pi(stickSNPs)

# fetch citations
citations(x)
#> Citations for methods used thus far: 
#> ==============================================
#> Statistic:  snpR 
#> 
#> Citation:  Hemstrom W, Jones M (2023). “snpR: User friendly population genomicsfor SNP data sets with categorical metadata.” _Molecular EcologyResources_, *23*(4), 962-973. ISSN 1755-098X.doi:10.1111/1755-0998.13721 <https://doi.org/10.1111/1755-0998.13721>.<https://doi.org/10.1111/1755-0998.13721>. 
#> 
#> Bibtex key:  hemstromSnpRUserFriendly2023 
#> 
#> Details:  snpR package 
#> ==============================================
#> Statistic:  pi 
#> 
#> Citation:  Hohenlohe PA, Bassham S, Etter PD, Stiffler N, Johnson EA, Cresko WA(2010). “Population Genomics of Parallel Adaptation in ThreespineStickleback using Sequenced RAD Tags.” _PLOS Genetics_, *6*(2), 1-23.doi:10.1371/journal.pgen.1000862<https://doi.org/10.1371/journal.pgen.1000862>.<https://doi.org/10.1371/journal.pgen.1000862>. 
#> 
#> Bibtex key:  Hohenlohe2010 
#> 
#> Details:  pi, number of pairwise differences 
#> ==============================================
#> 

# fetch citations and write bibtex .bib to citations.bib
## not run
if (FALSE) { # \dontrun{
citations(x, outbib = "citations.bib")
} # }
```
