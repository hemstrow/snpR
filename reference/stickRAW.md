# Deschutes stickleback raw data.

A dataset containing genotypes at 100 SNP loci for 100 Three-spined
Stickleback collected in 2015 from several locations in the Deschutes
River, Oregon, USA. For more information, see Hemstrom et al (in prep).

## Usage

``` r
stickRAW
```

## Format

An object of class `data.frame` with 100 rows and 102 columns.

## Details

Genotypes noted as two characters representing bases sequenced at a
specific site (e.g. "CC" or "GT"). Each SNP occupies one row, each
sample occupies on column. Missing data coded as "NN".

Column names refer to the population of origin (first three characters)
and the sample ID.

The first three columns contain snp specific meta data:

- snp: SNP ID.

- chr: Linkage group, essentially chromosome.

- position: Position of the SNP along the linkage group, in base pairs.

## References

Hemstrom et al (in prep)
