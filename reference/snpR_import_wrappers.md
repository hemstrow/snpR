# snpRdata Import Wrappers

These functions wrap
[`import.snpR.data`](https://hemstrow.github.io/snpR/reference/import.snpR.data.md)
to import data into the snpRdata format from a range of file or object
sources.

## Usage

``` r
read_vcf(
  file,
  snp.meta = NULL,
  sample.meta = NULL,
  accept_non_biallelic = FALSE
)

read_ms(file, sample.meta = NULL, chr.length, fix_overlaps = TRUE)

read_delimited_snps(
  file,
  snp.meta = NULL,
  sample.meta = NULL,
  mDat = "NN",
  header_cols = 0
)

read_genepop(file, snp.meta = NULL, sample.meta = NULL, mDat = "0000")

read_FSTAT(file, snp.meta = NULL, sample.meta = NULL, mDat = "0000")

read_plink(file)

read_structure(
  file,
  snp.meta = NULL,
  sample.meta = NULL,
  rows_per_individual = 2,
  marker_names = FALSE,
  header_cols = 0,
  mDat = -9
)

convert_genlight(genlight, snp.meta = NULL, sample.meta = NULL)

convert_genind(genind, snp.meta = NULL, sample.meta = NULL)

convert_vcfR(vcfR, snp.meta = NULL, sample.meta = NULL)
```

## Arguments

- file:

  character, path to a file containing genotype data to import.

- snp.meta:

  data.frame or character, default NULL. Metadata for each SNP, must
  have a number of rows equal to the number of SNPs in the dataset. If
  NULL, a single "snpID" column will be added. If a character, the path
  to a file containing SNP metadata, one row per SNP, with named
  columns.

- sample.meta:

  data.frame, default NULL. Metadata for each individual sample, must
  have a number of rows equal to the number of samples in the data set.
  If NULL, a single "sampID" column will be added. If a character, the
  path to a file containing sample metadata, one row per sample, with
  named columns.

- accept_non_biallelic:

  Logical, default FALSE. If FALSE, non-biallelic markers will be
  removed. If TRUE, this function passes to
  [`read_non_biallelic`](https://hemstrow.github.io/snpR/reference/read_non_biallelic.md)
  and allows non-biallelic markers. Note that snpR currently allows only
  a reduced set of import formats in this case.

- chr.length:

  numeric, Specifies chromosome lengths. Note that a single value
  assumes that each chromosome is of equal length whereas a vector of
  values gives the length for each chromosome in order.

- fix_overlaps:

  Logical, default TRUE. If TRUE, overlapping positions will be checked
  and fixed during 'ms' file import.

- mDat:

  character, defaults "0000", "NN", or "-9" depending on method. Note,
  if the default is set but the data has genotypes stored in 6
  characters, mDat will be set to "000000".

- header_cols:

  numeric, default 0. The number of snp metadata columns prior to snp
  genotypes when importing delimited snps.

- rows_per_individual:

  numeric (1 or 2), default 2. Number of rows used for each individual.

- marker_names:

  logical, default FALSE. If TRUE, assumes that a header row of marker
  is present.

- genlight:

  genlight object to convert, see
  [`genlight`](https://rdrr.io/pkg/adegenet/man/genlight.html).

- genind:

  genind object to convert, see
  [`genind`](https://rdrr.io/pkg/adegenet/man/new.genind.html).

- vcfR:

  vcfR object to convert, see vcfR documentation.

## Details

These functions are all wrappers for
[`import.snpR.data`](https://hemstrow.github.io/snpR/reference/import.snpR.data.md),
and all are technically cross-compatible save read_ms: each other
function can actually be called with any of the supported formats
(read_vcf can be handed a genlight object without failure). These are
supported as separate functions for code readability and for ease of
discovery.

See
[`import.snpR.data`](https://hemstrow.github.io/snpR/reference/import.snpR.data.md)
for more detail.

## Functions

- `read_vcf()`: Import .vcf or .vcf.gz files.

- `read_ms()`: Import .ms files.

- `read_delimited_snps()`: Import tab delimited data where genotypes are
  stored as: NN, 0000, or snp_tab format.

- `read_genepop()`: Import genepop formatted data.

- `read_FSTAT()`: Import FSTAT formatted data.

- `read_plink()`: Import plink bed, bim, and fam data.

- `read_structure()`: Import STRUCTURE data files.

- `convert_genlight()`: Convert adegenet genlight objects.

- `convert_genind()`: Convert adegenet genind objects

- `convert_vcfR()`: Convert adegenet vcfR objects

## File import

Supports automatic import of several types of files. Options: \*.vcf or
.vcf.gz: Variant Call Format (vcf) files, supported via the vcfR
package. If not otherwise provided, snp metadata is taken from the fixed
fields in the VCF and sample metadata from the sample IDs. Note that
this only imports SNPs with called genotypes!

- .ms: Files in the ms format, as provided by many commonly used
  simulation tools.

- NN: SNP genotypes stored as actual base calls (e.g. "AA", "CT").

- 0000: SNP genotypes stored as four numeric characters (e.g. "0101",
  "0204").

- snp_tab: SNP genotypes stored with genotypes in each cell, but only a
  single nucleotide noted if homozygote and two nucleotides separated by
  a space if heterozygote (e.g. "T", "T G").

- sn: SNP genotypes stored with genotypes in each cell as 0 (homozygous
  allele 1), 1 (heterozygous), or 2 (homozyogus allele 2).

- genepop: genepop file format, with genotypes stored as either 4 or 6
  numeric characters. Works only with bi-allelic data. Genotypes will be
  converted (internally) to NN: the first allele (numerically) will be
  coded as A, the second as C.

- FSTAT: FSTAT file format, with genotypes stored as either 4 or 6
  numeric characters. Works only with bi-allelic data. Genotypes will be
  converted (internally) to NN: the first allele (numerically) will be
  coded as A, the second as C.

- plink: plink .bed, .fam, and .bim files, via
  [`read_plink`](https://rdrr.io/pkg/genio/man/read_plink.html). If any
  of these file types is provided, snpR via
  [`read_plink`](https://rdrr.io/pkg/genio/man/read_plink.html) will
  look for the other file types automatically. Sample metadata should be
  contained in the .fam file and SNP metadata in the .bim file, so
  sample or snp meta data can be provided here.

- structure: STRUCTURE import file, with individuals in rows and loci in
  columns. Can be coded either with one row per individual and two
  columns per loci or two rows per individual and two columns per loci
  using the rows_per_individual argument. Genotypes can be pretty much
  anything, although missing genotypes must be coded as -9. Must have a
  .str extension and be consistantly whitespace delimited

Sample and snp metadata can also be provided via file path, and will be
read in using [`fread`](https://rdrr.io/pkg/data.table/man/fread.html)
*with the default settings* using `read_delimited_snps`. If these
settings are not correct, please read in the metadata manually and
provide to import.snpR.data.

## Conversions from other S4 objects

Supports automatic conversions from some other popular S4 object types.
Options:

- genind: [`genind`](https://rdrr.io/pkg/adegenet/man/new.genind.html)
  objects from adegenet. Note, no need to import genepop objects, the
  equivalent statistics are calculated automatically when functions
  called with facets. Sample and SNP IDs as well as, when possible, pop
  IDs will be taken from the genind object. This data will be added too
  but will not replace data provided to the SNP or sample.meta
  arguments. Note that only *SNP* data is currently allowed, data with
  more than two alleles for loci will return an error.

- genlight: [`genlight`](https://rdrr.io/pkg/adegenet/man/genlight.html)
  objects from adegenet. Sample and SNP IDs, SNP positions, SNP
  chromosomes, and pop IDs will be taken from the genlight object if
  possible. This data will be added too but will not replace data
  provided to the SNP or sample.meta arguments.

- vcfR: vcfR objects from vcfR. If not provided, snp metadata is taken
  from the fixed fields in the VCF and sample metadata from the sample
  IDs. Note that this only imports SNPs with called genotypes!

## Author

William Hemstrom

Brent Gruber (genlight conversion re-distributed here)
