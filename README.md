
<!-- README.md is generated from README.Rmd. Please edit that file -->

# snpR

<!-- badges: start -->

[![license](https://img.shields.io/badge/license-MIT%20+%20file%20LICENSE-lightgrey.svg)](https://choosealicense.com/)
[![packageversion](https://img.shields.io/badge/Package%20version-1.2.0.9000-orange.svg?style=flat-square)](commits/master)
[![CRAN
status](https://www.r-pkg.org/badges/version/snpR)](https://CRAN.R-project.org/package=snpR)
<!-- badges: end -->

snpR is an R package for analyzing call Single Nucleotide Polymorphism
(SNP) genotypes containing most basic stats including pairwise LD,
gaussian sliding window analysis tools, plotting options, clustering
analysis, colony interface, Ne estimation, formatting, filtering, and
more!

## Installation

snpR can be installed from [GitHub](https://github.com/):

``` r
# install.packages("remotes")
remotes::install_github("hemstrow/snpR")
```

To install the vignettes as well (recommended for new users), instead
use:

``` r
remotes::install_github("hemstrow/snpR", build_vignettes = T) # linux
remotes::install_github("hemstrow/snpR", ref = "dev", build_vignettes = T, build_opts = c("--no-resave-data", "--no-manual")) # windows
```

The dev version can be installed from [GitHub](https://github.com/) as
well:

``` r
# install.packages("remotes")
remotes::install_github("hemstrow/snpR", ref = "dev")
```

A CRAN version should be available soon.

## Example

snpR is focused on ease-of-use. Primarily, it achieves this via the use
of , which describe sample or SNP metadata. snpR is built to
automatically split up analysis by facet. For example, calculating
genetic diversity (pi) for each population or family, or for each
population/family combination is easy!

``` r
library(snpR)
#> Loading required package: data.table
#> Loading required package: foreach
## basic example code

x <- calc_pi(stickSNPs, facets = c("pop")) # split by pop (stickSNPs is an example dataset included in snpR)
x <- calc_pi(x, facets = c("fam")) # split by family
x <- calc_pi(x, facets = c("pop.fam")) # split by combinations of family and pop
```

snpR also facilitates ease-of-use by being . As above, new analyses are
added to an existing object. Results can be fetched using the
get.snpR.stats handler.

``` r
head(get.snpR.stats(x))
#>   facet subfacet      chr position .snp.id pi major minor maj.count min.count
#> 1 .base    .base   groupV    42825       1 NA     C     T       791        23
#> 2 .base    .base groupXIX    67921       2 NA     C     G       515       213
#> 3 .base    .base  groupIX   100382       3 NA     G     A       739        53
#> 4 .base    .base   groupX   101821       4 NA     C     T       752        46
#> 5 .base    .base groupXIV   175941       5 NA     G     A       719        53
#> 6 .base    .base   groupI   182629       6 NA     C     T       642        22
#>          maf
#> 1 0.02825553
#> 2 0.29258242
#> 3 0.06691919
#> 4 0.05764411
#> 5 0.06865285
#> 6 0.03313253
```

For a full introduction, check the snpR_introduction vignette.

``` r
# remotes::install_github("hemstrow/snpR", ref = "dev", build_vignettes = T, build_opts = c("--no-resave-data", "--no-manual"))

vignette("snpR_introduction")
#> Warning: vignette 'snpR_introduction' not found
```
