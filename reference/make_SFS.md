# Generate a 1-2d site frequency spectrum from a dadi input file.

Generates a 1 or 2 dimensional site frequency spectrum from a dadi input
file using the projection methods and folding methods of Marth et al
(2004) and Gutenkunst et al (2009). This code is essentially an R
re-implementation of the SFS construction methods implemented in the
program *dadi* (see Gutenkunst et al (2009)).

## Usage

``` r
make_SFS(x, pops, projection, fold = FALSE, update_bib = FALSE)
```

## Arguments

- x:

  character or data.frame. Either a path to a dadi formatted input file
  or a data.frame containing previously imported dadi formatted data.

- pops:

  character. A vector of population names of up to length 2 containing
  the names of populations for which the an SFS is to be created.

- projection:

  numeric. A vector of sample sizes to project the SFS to, in *number of
  gene copies*. Sizes too large will result in a SFS containing few or
  no SNPs.

- fold:

  logical, default FALSE. Determines if the SFS should be folded or left
  polarized.

- update_bib:

  character or FALSE, default FALSE. If a file path to an existing .bib
  library or to a valid path for a new one, will update or create a .bib
  file including any new citations for methods used. Useful given that
  this function does not return a snpRdata object, so a
  [`citations`](https://hemstrow.github.io/snpR/reference/citations.md)
  cannot be used to fetch references.

## Value

A matrix or vector containing the site frequency spectrum with a "pops"
attribute containing population IDs, such as c("POP1", "POP2"). For a 2d
SFS, the first pop is the matrix columns and the second is the matrix
rows.

## Details

Site frequency spectrums are constructed using the projection methods
detailed in Marth et al (2004) and the 2 dimensional expansion in
Gutenkunst et al (2009). Folding methods are also taken from Gutenkunst
et al (2009). Either 1 or 2d SFSs can be constructed by providing a
vector of population names and projection sizes.

## References

Gutenkunst et al (2009). Inferring the joint demographic history of
multiple populations from multidimensional SNP frequency data. *PLoS
genetics*, 5(10), e1000695.

Marth et al (2004). The allele frequency spectrum in genome-wide human
variation data reveals signals of differential demographic history in
three large world populations. *Genetics*, 166(1), 351-372.
