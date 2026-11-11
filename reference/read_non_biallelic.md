# Read in non-biallelic genotype data.

Reads in non-biallelic data stored as a vcf or tab-delimited data.
Currently only supports these data types.

## Usage

``` r
read_non_biallelic(
  genotypes,
  snp.meta = NULL,
  sample.meta = NULL,
  header_cols = 0,
  mDat = "0000",
  verbose = FALSE,
  .pass_filters = FALSE
)
```

## Arguments

- genotypes:

  data.frame, vcfR object, or a file path to either vcf or tab-delimited
  file containing raw genotypes with SNPs are in rows and individual
  samples are in columns. If raw genotypes, each genotype must be stored
  in one cell (as in "001002" or "0102") with a matching length provided
  to the "mDat" argument.

- snp.meta:

  data.frame, default NULL. Metadata for each SNP, must have a number of
  rows equal to the number of SNPs in the dataset. If NULL, a single
  "snpID" column will be added.

- sample.meta:

  data.frame, default NULL. Metadata for each individual sample, must
  have a number of rows equal to the number of samples in the dataset.
  If NULL, a single "sampID" column will be added.

- header_cols:

  numeric, default 0. Number of header columns containing SNP metadata.
  Used if a tab delimited or STRUCTURE input file is provided.

- mDat:

  character, default "000000", matching the encoding of missing
  *genotypes* in the data provided to the genotypes argument.

- verbose:

  Logical, default FALSE. If TRUE, will print a few status updates and
  checks.

- .pass_filters:

  Internal, probably not for user use. Used to pass filtering history
  when sub-setting when this function is called internally.
