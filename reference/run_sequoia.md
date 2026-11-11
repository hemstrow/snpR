# Run Sequoia pedigree/parentage assignment with snpR

Runs the parentage assignment and pedigree construction tool from the
`sequoia` package. Note that this function *is not overwrite safe!*.

## Usage

``` r
run_sequoia(
  x,
  facets = NULL,
  run_dupcheck = FALSE,
  run_parents = FALSE,
  run_pedigree = FALSE,
  run_relatives = FALSE,
  min_maf = 0.3,
  min_ind = 0.5,
  ...
)
```

## Arguments

- x:

  snpRdata object.

- facets:

  character, default NULL. Sample-specific facets over which the sequoia
  is called to run. See
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).

- run_dupcheck:

  FALSE or TRUE, default FALSE. Uses sequoia to check for duplicate
  samples in the dataset. Duplicate samples should not be included for
  parentage and pedigree construction.

- run_parents:

  FALSE or TRUE, default FALSE. Runs parentage assignments for the
  samples. This runs quickly and is required before using the
  run_pedigree command.

- run_pedigree:

  FALSE or TRUE, default FALSE. Runs pedigree construction for the
  samples. This process can take a long time.

- run_relatives:

  FALSE or TRUE, default FALSE. Runs retrieval of other relatives which
  did not pass thresholds for assignment in the main pedigree
  construction model.

- min_maf:

  numeric in 0.25:0.5, default 0.3. Minimum allele frequency cutoff for
  analysis. Sequoia requires high minor allele frequencies for parentage
  and pedigree construction.

- min_ind:

  numeric in 0.5:1, default 0.5. Removes loci sequenced in less than
  this proportion of individuals. Note that *individuals* with genotypes
  for fewer than half of the loci will be automatically removed by
  sequoia.

- ...:

  Additional arguments passed to`sequoia` (during parentage and pedigree
  reconstruction).

## Value

A nested list with each facet specified containing sequoia output
summary information.

## Details

This is a limited integration of the program and package written by
Jisca Huisman. Note that there are many more Sequoia specific arguments
that can be added to change from the default settings (eg. ErrorM,
Tassign, Tfilt, etc.) See documentation for `sequoia`. These can be
passed to the pedigree and parentage reconstructions using the ...
argument in run_sequoia. The sequoia package has many features, and snpR
facilitates the use of a fraction of them. snpR users are encouraged to
use the sequoia R package.

## References

Huisman,J. (2017) Pedigree reconstruction from SNP data: parentage
assignment, sibship clustering and beyond. Mol. Ecol. Resour., 17,
1009–1024.

## Author

William Hemstrom

Melissa Jones

## Examples

``` r
# to follow an example using the stickSNPs example dataset you need 
# to add some variables that don't exist in the actual dataset.
a <- 2013:2015 #create a vector of possible birthyears
b <- c("M", "F", "U") #create a vector of possible sexes
stk <- stickSNPs
set.seed(4865)
sample.meta(stk)$BirthYear <- sample(x = a, size = nsamps(stickSNPs), 
                                     replace = TRUE) #create birth years
sample.meta(stk)$ID <- 1:nsamps(stk) #create unique sampleID
sample.meta(stk)$Sex <- sample(x= b, size = nsamps(stk), 
                               replace = TRUE) # create sexes

# slow, so not run here
if (FALSE) { # \dontrun{
dup <- run_sequoia(x = stk, run_dupcheck = TRUE, run_parents = FALSE, 
                   run_pedigree = FALSE, run_relatives = FALSE)
ped <- run_sequoia(x = stk, run_dupcheck = FALSE, run_parents = TRUE, 
                   run_pedigree = TRUE, run_relatives = FALSE)
rel <- run_sequoia(x = stk, run_dupcheck = FALSE, run_parents = FALSE, 
                   run_pedigree = FALSE, run_relatives = TRUE)
                   
} # }
```
