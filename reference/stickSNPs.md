# Deschutes stickleback example snpRdata.

A dataset containing genotypes at 100 SNP loci for 100 Three-spined
Stickleback collected in 2015 from several locations in the Deschutes
River, Oregon, USA. For more information, see Hemstrom et al (in prep).

## Usage

``` r
stickSNPs
```

## Format

An object of class `snpRdata` with 100 rows and 100 columns.

## Details

Stored in the snpRdata object class. See
[`import.snpR.data`](https://hemstrow.github.io/snpR/reference/import.snpR.data.md)
for details.

Genotypes noted as two characters representing bases sequenced at a
specific site (e.g. "CC" or "GT"). Each SNP occupies one row, each
sample occupies on column. Missing data coded as "NN".

Sample and snp metadata can be found in the sample.meta and snp.meta
sockets, respectively.

## References

Hemstrom et al (in prep)
