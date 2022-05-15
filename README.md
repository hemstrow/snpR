
<!-- README.md is generated from README.Rmd. Please edit that file -->

# snpR

<!-- badges: start -->

[![packageversion](https://img.shields.io/badge/Package%20version-1.2.1.9000-orange.svg?style=flat-square)](commits/master)
[![CRAN
status](https://www.r-pkg.org/badges/version/snpR)](https://CRAN.R-project.org/package=snpR)
[![R-CMD-check](https://github.com/hemstrow/snpR/workflows/R-CMD-check/badge.svg)](https://github.com/hemstrow/snpR/actions)
[![Codecov test
coverage](https://codecov.io/gh/hemstrow/snpR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/hemstrow/snpR?branch=master)
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

snpR also facilitates ease-of-use by being . As above, new analyses are
added to an existing object. Results can be fetched using the
get.snpR.stats handler.

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
