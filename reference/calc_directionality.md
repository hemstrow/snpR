# Calculate the directionality index from a 2d site frequency spectrum.

Calculates the directionality index based on a 2d SFS according to Peter
and Slatkin (2013). Input spectra can be created using the
[`calc_sfs`](https://hemstrow.github.io/snpR/reference/calc_sfs.md)
function, using a provided snpRdata object, or passed from other
programs. *Spectra must be not be folded*.

## Usage

``` r
calc_directionality(
  x,
  facet = NULL,
  pops = NULL,
  projection = NULL,
  update_bib = FALSE
)
```

## Arguments

- x:

  snpRdata object or matrix, default NULL. A snpRdata from which to
  calculate a SFS. Alternatively, a 2d site frequency spectra stored in
  a matrix, with an additional "pop" or "pops" attribute containing
  population IDs, such as c("POP1", "POP2"), where the first pop
  corresponds to matrix columns and the second to matrix rows. These
  objects can be produced from a dadi input file using
  [`make_SFS`](https://hemstrow.github.io/snpR/reference/make_SFS.md).
  Note that if x is a snpRdata object, snp metadata columns named "ref"
  and "anc" containing the identity of the derived and ancestral
  alleles, respectively, must be present.

- facet:

  character, default NULL. Passed to
  [`calc_sfs`](https://hemstrow.github.io/snpR/reference/calc_sfs.md) –
  see documentation there for details. Ignored if a sfs is provided.

- pops:

  character, default NULL. Passed to
  [`calc_sfs`](https://hemstrow.github.io/snpR/reference/calc_sfs.md) –
  see documentation there for details. Ignored if a sfs is provided.

- projection:

  numeric, default NULL. Passed to
  [`calc_sfs`](https://hemstrow.github.io/snpR/reference/calc_sfs.md) –
  see documentation there for details. Ignored if a sfs is provided.

- update_bib:

  character or FALSE, default FALSE. If a file path to an existing .bib
  library or to a valid path for a new one, will update or create a .bib
  file including any new citations for methods used. Useful given that
  this function does not return a snpRdata object, so a
  [`citations`](https://hemstrow.github.io/snpR/reference/citations.md)
  cannot be used to fetch references.

## Value

A numeric value giving the directionality with a "direction" attribute
designating the direction between the two populations.

## Details

Essentially, the directionality index measures the difference in derived
allele frequency between two populations to determine the directionality
of population spread between the two. Since the "destination" population
is sourced from but experienced more genetic drift than the "source"
population, it should have relatively more high-frequency derived
alleles *after the removal of fixed ancestral alleles*. See Peter and
Slatkin (2013) for details.

## References

Peter, B. M., & Slatkin, M. (2013). Detecting range expansions from
genetic data. *Evolution*, 67(11), 3274-3289.

## Examples

``` r
if (FALSE) { # \dontrun{
# directionality can be calculated without first calculating a SFS
calc_directionality(stickSNPs, facet ="pop", pops = c("ASP", "PAL"), projection = c(10, 10))

# an existing SFS can also be fed in. This may be handy if you get a SFS from elsewhere.
sfs <- calc_sfs(stickSNPs, "pop", c("ASP", "PAL"), c(10, 10), fold = FALSE)
calc_directionality(sfs)
} # }
```
