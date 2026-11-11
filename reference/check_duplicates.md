# Checks for duplicated samples in snpRdata.

Searches through a snpR dataset and, for every designated sample,
determines the proportion of identical genotypes in every other sample.
This function *is not overwrite safe*.

## Usage

``` r
check_duplicates(x, y = 1:ncol(x), id.col = NULL, verbose = FALSE)
```

## Arguments

- x:

  snpRdata object

- y:

  numeric or character, default 1:ncol(x). Designates the sample indices
  or IDs in x for which duplicates will be checked.

- id.col:

  character, default NULL. Designates a column in the sample metadata
  which contains sample IDs. If provided, y is assumed to contain sample
  IDs uniquely matching those in the the sample ID column.

- verbose:

  logical, default FALSE. If TRUE, prints detailed progress report.

## Value

A list containing:

- best_matches: Data.frame listing the best match for each sample noted
  in y and the percentage of genotypes identical between the two
  samples.

- data: A list containing the match proportion between each sample y and
  every sample in x, named for the samples y.

## Details

If an id column is specified, y should contain sample IDs matching those
contained in that column. If not, y should contain sample indices
instead. The proportion of identical genotypes between matching samples
and all other samples are calculated. By default, every sample will be
checked.

## Author

William Hemstrom

## Examples

``` r
if (FALSE) { # \dontrun{
# check for duplicates with sample 1
check_duplicates(stickSNPs, 1)

# check duplicates using the .samp.id column as sample IDs
check_duplicates(stickSNPs, 1, id.col = ".sample.id")
} # }
```
