# Utilities for direct vcf file conversion

Format vcf files as alternative formats without going through a
[`snpRdata`](https://hemstrow.github.io/snpR/reference/snpRdata-class.md)
intermediate. Used for working with likelihoods, probabilities, or other
non-GT format fields.

## Usage

``` r
vcf2beagle(file, outfile)

vcf2pl(file, outfile, init_admix = FALSE, k = NULL)
```

## Arguments

- file:

  character, path to "vcf" file.

- outfile:

  character, output filename.

- init_admix:

  logical, default FALSE. If TRUE, will write initial admix estimates
  for the program `entropy` to a set of files for each `k`.

- k:

  integer vector, default NULL. Numbers of ancestry clusters (*k*
  values)for which to estimate intial admixture proportions for the
  program `entropy`. Ignored if `init_admix` is FALSE.

## Details

These functions are provided as a utility and do not work in the general
`snpR` ecosystem. The different outputs may require specific FORMAT
tags:

- vcf2beagle: `GP` tag

- vcf2PL: `PL` tag If these are lacking, `snpR` may provide
  recommendations by which they can be added if possible.

## Functions

- `vcf2beagle()`: Format as a "beagle" likelihood file.

- `vcf2pl()`: Format as a "mgpl" phred-scaled likelihood file.
