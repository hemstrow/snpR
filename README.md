# snpR
An R package for analyzing call SNP genotypes containing most basic stats including pairwise LD, gaussian sliding window analysis tools, plotting options, clustering analysis, colony interface, Ne estimation, formatting, filtering, and more!

## Installation

snpR can be installed using the install_github functions from either the devtools or remotes packages:

remotes::install_github("hemstrow/snpR")

To install the vignettes, instead use:

remotes::install_github("hemstrow/snpR", build_vignettes = T)

for Linux, or

remotes::install_github("hemstrow/snpR", ref = "dev", build_vignettes = T, build_opts = c("--no-resave-data", "--no-manual"))

for Windows
