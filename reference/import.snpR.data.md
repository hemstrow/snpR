# Import genotype and metadata into a snpRdata object.

`import.snpR.data` converts genotype and meta data to the snpRdata
class, which stores raw genotype data, sample and locus specific
metadata, useful data summaries, repeatedly internally used tables,
calculated summary statistics, and sliding-window statistic data.

## Usage

``` r
import.snpR.data(
  genotypes,
  snp.meta = NULL,
  sample.meta = NULL,
  mDat = "NN",
  chr.length = NULL,
  ...,
  accept_non_biallelic = FALSE,
  header_cols = 0,
  rows_per_individual = 2,
  marker_names = FALSE,
  fix_overlaps = TRUE,
  verbose = FALSE,
  .pass_filters = FALSE,
  .skip_filters = FALSE
)
```

## Arguments

- genotypes:

  data.frame, unique S4 from other packages, or filename. If a
  data.frame, raw genotypes in a two-character format ("GG", "GA", "CT",
  "NN"), where SNPs are in rows and individual samples are in columns.
  Otherwise, see documentation for allowed S4 objects and files.

- snp.meta:

  data.frame, default NULL. Metadata for each SNP, must have a number of
  rows equal to the number of SNPs in the dataset. If NULL, a single
  "snpID" column will be added.

- sample.meta:

  data.frame, default NULL. Metadata for each individual sample, must
  have a number of rows equal to the number of samples in the dataset.
  If NULL, a single "sampID" column will be added.

- mDat:

  character, default "NN", matching the encoding of missing *genotypes*
  in the data provided to the genotypes argument.

- chr.length:

  numeric, default NULL. If a path to a .ms file is provided, specifies
  chromosome lengths. Note that a single value assumes that each
  chromosome is of equal length whereas a vector of values gives the
  length for each chromosome in order.

- ...:

  Additional arguments passed to
  [`fread`](https://rdrr.io/pkg/data.table/man/fread.html) if a
  *genotype* file name is passed that is not a vcf or ms file.

- accept_non_biallelic:

  Logical, default FALSE. If FALSE, non-biallelic markers will be
  removed. If TRUE, this function passes to
  [`read_non_biallelic`](https://hemstrow.github.io/snpR/reference/read_non_biallelic.md)
  and allows non-biallelic markers. Note that snpR currently allows only
  a reduced set of import formats in this case.

- header_cols:

  numeric, default 0. Number of header columns containing SNP metadata.
  Used if a tab delimited or STRUCTURE input file is provided.

- rows_per_individual:

  numeric (1 or 2), default 2. Number of rows used for each individual.
  For structure input files only.

- marker_names:

  logical, default FALSE. If TRUE, assumes that a header row of marker
  is present. For structure input files only.

- fix_overlaps:

  Logical, default TRUE. If TRUE, overlapping positions will be checked
  and fixed during 'ms' file import.

- verbose:

  Logical, default FALSE. If TRUE, will print a few status updates and
  checks.

- .pass_filters:

  Internal, probably not for user use. Used to pass filtering history
  when sub-setting when this function is called internally.

- .skip_filters:

  Internal, probably not for user use. Used to skip re-filtering during
  sub-setting when this function is called internally.

## Details

The snpRdata class is built to contain SNP genotype data for use by
functions in the snpR package. It inherits from the S3 class data.frame,
in which the genotypes are stored, and can be manipulated identically.
It also stores sample and locus specific metadata, genomic summary
information, and results from most snpR functions. Genotypes are stored
in the "character" format, as output by
[`format_snps`](https://hemstrow.github.io/snpR/reference/format_snps.md).
Missing data is noted with "NN".

Inputs can be provided either as pre-existing R objects (in several
different formats) or as paths to files. In both cases, `snpR` will
attempt to guess the data format from either the object classs, first
genotype or the file extension (as appropriate).

## File import

Supports automatic import of several types of files. Options:

- .vcf or .vcf.gz: Variant Call Format (vcf) files, supported via vcfR.
  If not otherwise provided, snp metadata is taken from the fixed fields
  in the VCF and sample metadata from the sample IDs. Note that this
  only imports SNPs with called genotypes!

- .ms: Files in the ms format, as provided by many commonly used
  simulation tools.

- .txt, NN: SNP genotypes stored as actual base calls (e.g. "AA", "CT").

- .txt, 0000: SNP genotypes stored as four numeric characters (e.g.
  "0101", "0204").

- .txt, snp_tab: SNP genotypes stored with genotypes in each cell, but
  only a single nucleotide noted if homozygote and two nucleotides
  separated by a space if heterozygote (e.g. "T", "T G").

- .txt, sn: SNP genotypes stored with genotypes in each cell as 0
  (homozygous allele 1), 1 (heterozygous), or 2 (homozyogus allele 2).

- .genepop or .gen: genepop file format, with genotypes stored as either
  4 or 6 numeric characters. Works only with bi-allelic data. Genotypes
  will be converted (internally) to NN: the first allele (numerically)
  will be coded as A, the second as C.

- .fstat: FSTAT file format, with genotypes stored as either 4 or 6
  numeric characters. Works only with bi-allelic data. Genotypes will be
  converted (internally) to NN: the first allele (numerically) will be
  coded as A, the second as C.

- .bed/.fam/.bim: PLINK .bed, .fam, and .bim files, via
  [`read_plink`](https://rdrr.io/pkg/genio/man/read_plink.html).

If any of these file types is provided, snpR (via
[`read_plink`](https://rdrr.io/pkg/genio/man/read_plink.html)) will look
for the other file types automatically. Sample metadata should be
contained in the .fam file and SNP metadata in the .bim file, so sample
or snp meta data provided here will be ignored.

- .str: STRUCTURE import files in either 1 or 2 rows per individual as
  defined by the `rows_per_individual` argument.

Additional arguments can be provided to import.snpR.data that will be
passed to [`fread`](https://rdrr.io/pkg/data.table/man/fread.html) when
reading in genotype data.

Sample and snp metadata can also be provided via file path, and will be
read in using [`fread`](https://rdrr.io/pkg/data.table/man/fread.html)
*with the default settings*. If these settings are not correct, please
read in the metadata manually and provide to import.snpR.data.

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

## Slots

Genotypes, metadata, and results are stored in slots and directly
accessable with the 'at' symbol operator. Slots are as follows:

- sample.meta: sample metadata (population, family, phenotype, etc.).

- snp.meta: SNP metadata (SNP ID, chromosome, linkage group, position,
  etc.).

- facet.meta: internal metadata used to track facets that have been
  previously applied to the dataset.

- mDat:missing data format.

- snp.form: number of characters per SNP.

- genotables: a list containing tabulated genotypes (gs), allele counts
  (as), and missing data (wm). facet.meta contains the corresponding
  metadata.

- facets: vector of the facets that have been added to the data.

- facet.type: classes of the added facets (snp, sample, complex, or
  .base).

- stats: data.frame containing all calculated non-pairwise single-snp
  statistics and metadata.

- window.stats: data.frame/table containing all non-pairwise statistics
  calculated for sliding windows.

- pairwise.stats: data.frame/table containing all pairwise (fst)
  single-snp statistics.

- pairwise.window.stats: data.frame/table containing all pairwise
  statistics calculated for sliding windows.

- sample.stats: data.frame/table containing statistics calculated for
  each individual sample.

- pairwise.LD:nested list containing linkage disequilibrium data (see
  [`calc_pairwise_ld`](https://hemstrow.github.io/snpR/reference/calc_pairwise_ld.md)
  for more information).

- window.bootstraps: data.frame/table containing all calculated
  bootstraps for sliding window statistics.

- sn: list containing "sn", sn formatted data, and "type" type of
  interpolation.

- calced_stats: Named list of named character vectors that tracks the
  calculated statistics for each facet (see
  [`calc_genetic_distances`](https://hemstrow.github.io/snpR/reference/calc_genetic_distances.md)
  for more information).

- genetic_distances: nested list containing genetic distance data.

- names: column names for genotypes.

- row.names: row names for genotypes.

- .Data: list of vectors containing raw genotype data.

- .S3Class: notes the inherited S3 object class.

Note that most of these slots are used primarily internally.

All calculated data can be accessed using the
[`get.snpR.stats`](https://hemstrow.github.io/snpR/reference/get.snpR.stats.md)
function. See documentation.

## Author

William Hemstrom

## Examples

``` r
# import example data as a snpRdata object
# produces data identical to that contained in the stickSNPs example dataset.
genos <- stickRAW[,-c(1:2)]
snp_meta <- stickRAW[,1:2]
sample_meta <- data.frame(pop = substr(colnames(stickRAW)[-c(1:2)], 1, 3), 
                          fam = rep(c("A", "B", "C", "D"), 
                                    length = ncol(stickRAW) - 2), 
                          stringsAsFactors = FALSE)
import.snpR.data(genos, snp.meta = snp_meta, sample.meta = sample_meta, 
                 mDat = "NN")
#>     ASP.r2.141 ASP.r2.144 ASP.r2.146 ASP.r2.147 ASP.r2.196 ASP.r2.216
#> 1           AG         GG         AA         AA         NN         AG
#> 2           GG         GG         GG         GG         GT         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AA         AA         AA         AA         AG         AA
#> 5           GG         NN         NN         AA         GG         AG
#> 6           AC         AC         AA         AA         AA         AC
#> 7           TT         TT         TT         TT         GT         TT
#> 8           GG         GG         AG         AG         AG         GG
#> 9           CC         AC         CC         AC         CC         CC
#> 10          AA         AA         AG         AA         AA         AA
#> 11          NN         TT         TT         CT         TT         CT
#> 12          AG         AG         AG         NN         AA         AA
#> 13          TT         TT         CT         TT         CT         NN
#> 14          GG         GT         GT         NN         GG         GG
#> 15          AA         AG         AG         AA         AA         AA
#> 16          NN         GG         GT         TT         GG         GT
#> 17          AA         AA         AA         AA         AA         NN
#> 18          AC         NN         CC         CC         AC         AC
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         AG         GG         GG         NN         GG
#> 21          NN         NN         CC         CC         AC         AC
#> 22          CT         NN         TT         CC         CT         NN
#> 23          GG         GG         GG         GG         GG         GG
#> 24          AA         AA         AA         NN         AA         AA
#> 25          NN         AG         AG         GG         GG         AG
#> 26          CT         CT         CT         CT         NN         TT
#> 27          AA         AA         AA         AA         AA         AA
#> 28          TT         NN         NN         NN         CT         CT
#> 29          AA         AC         AA         NN         AA         AA
#> 30          GG         GG         CG         GG         GG         GG
#> 31          CG         CG         CC         CC         CC         CC
#> 32          AA         AA         AA         NN         AA         AA
#> 33          CC         CT         CC         CC         CT         CT
#> 34          CG         CG         CG         CG         CC         NN
#> 35          AA         AA         AA         AC         AA         NN
#> 36          AA         AA         AA         AA         AA         AG
#> 37          CG         GG         GG         GG         GG         GG
#> 38          GG         TT         GG         GT         GT         TT
#> 39          TT         NN         TT         NN         AT         TT
#> 40          CG         CC         CG         NN         CG         GG
#> 41          AG         AA         AG         AG         AA         AA
#> 42          AG         NN         NN         NN         GG         GG
#> 43          AG         AG         AA         AA         AA         AA
#> 44          AT         TT         TT         TT         NN         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          NN         CC         AA         AA         NN         CC
#> 47          CC         CC         CC         CG         CG         CC
#> 48          CT         TT         TT         NN         TT         TT
#> 49          GG         GG         GG         GG         NN         GG
#> 50          NN         CC         CT         CC         CC         CT
#> 51          TT         NN         CT         CT         TT         CT
#> 52          AA         AG         AG         AA         AG         AG
#> 53          AA         AA         CC         NN         AA         AC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          GG         NN         AG         NN         AG         AG
#> 56          AG         AG         AG         AG         AA         NN
#> 57          GG         NN         GG         NN         GG         GG
#> 58          TT         CT         CC         NN         TT         CC
#> 59          AA         AG         GG         AG         GG         AA
#> 60          AA         AA         AC         NN         AC         AC
#> 61          TT         NN         TT         NN         TT         TT
#> 62          CG         GG         GG         NN         GG         GG
#> 63          AC         NN         NN         NN         AA         AC
#> 64          TT         AT         AT         TT         AT         AT
#> 65          CT         CC         CT         CC         CC         CT
#> 66          AC         AC         NN         NN         AC         CC
#> 67          GG         GG         GG         NN         GG         GG
#> 68          AA         AT         AT         AA         AA         NN
#> 69          AA         AA         NN         AA         AA         AA
#> 70          AC         AC         AA         AC         AC         NN
#> 71          AA         AA         AA         AT         AA         AA
#> 72          GG         GG         GG         GG         AG         NN
#> 73          TT         GT         TT         TT         GT         TT
#> 74          GT         GG         GT         GG         GG         GG
#> 75          AC         AC         NN         AC         CC         AC
#> 76          CC         CT         CT         CT         NN         NN
#> 77          AT         AT         AT         AT         AA         AA
#> 78          GG         GG         GG         GG         GG         GG
#> 79          TT         TT         CT         TT         CT         CT
#> 80          CT         CT         TT         TT         TT         CT
#> 81          AA         NN         AT         AT         AT         AT
#> 82          CT         TT         CT         TT         TT         TT
#> 83          AT         TT         AT         TT         TT         AT
#> 84          CC         CC         CT         CT         CC         NN
#> 85          AG         AG         GG         AA         AG         NN
#> 86          AG         AG         AG         AG         AA         AA
#> 87          CT         CC         CT         CT         TT         CT
#> 88          AG         NN         AG         NN         GG         AG
#> 89          NN         CG         CG         CG         CC         CG
#> 90          CT         CT         CC         CC         CC         CT
#> 91          CC         CC         NN         CT         CC         CT
#> 92          AA         NN         GG         NN         GG         AA
#> 93          TT         TT         TT         TT         CT         TT
#> 94          CC         CC         CC         CC         CC         NN
#> 95          GG         GG         GG         NN         GT         GG
#> 96          CC         NN         CC         NN         CC         CC
#> 97          GG         GG         GG         GG         GG         GG
#> 98          GG         GG         GG         NN         GG         GT
#> 99          AA         GG         AA         AA         AG         AG
#> 100         AT         NN         AT         AT         TT         AT
#>     ASP.r2.221 ASP.r2.227 ASP.r2.230 ASP.r2.232 ASP.r2.235 CLF.r2.040
#> 1           NN         GG         AA         GG         AG         NN
#> 2           GG         NN         GG         GG         GG         GG
#> 3           CT         CC         CC         CT         CC         CC
#> 4           AA         AA         AA         AG         AG         AA
#> 5           AG         AG         NN         GG         NN         GG
#> 6           AA         AC         AA         AC         NN         AA
#> 7           TT         TT         TT         TT         TT         GG
#> 8           GG         AG         GG         AG         GG         GG
#> 9           CC         CC         AC         CC         AA         CC
#> 10          AA         AA         AG         AA         AA         AA
#> 11          CC         CC         CT         CT         CT         CT
#> 12          AG         AA         AG         AA         AA         AA
#> 13          TT         TT         NN         TT         TT         CT
#> 14          GG         GG         GT         GG         GT         GG
#> 15          AA         AA         AG         AG         AG         AA
#> 16          GT         GG         GG         GG         GG         GT
#> 17          AC         AC         AA         AC         AA         AA
#> 18          CC         AC         AC         CC         AC         CC
#> 19          TT         TT         TT         TT         TT         GT
#> 20          GG         GG         GG         NN         GG         GG
#> 21          CC         CC         AC         CC         CC         NN
#> 22          CC         TT         CT         CC         CT         CC
#> 23          GG         AG         GG         GG         NN         GG
#> 24          AA         AA         AA         AA         AG         AA
#> 25          AG         GG         AG         AG         GG         GG
#> 26          NN         TT         CT         TT         CT         TT
#> 27          AA         AA         NN         AC         AA         AA
#> 28          TT         TT         TT         TT         TT         CT
#> 29          AA         NN         AA         AA         AA         AA
#> 30          GG         GG         GG         GG         CG         GG
#> 31          CG         CC         CG         CG         CG         CC
#> 32          AA         AC         AC         AC         AA         AC
#> 33          CT         CC         CC         CT         CT         CT
#> 34          CC         CG         CG         CG         NN         CC
#> 35          AC         AA         AA         AA         AA         AA
#> 36          AA         AG         AA         AG         NN         AG
#> 37          NN         GG         CC         GG         CG         GG
#> 38          GG         GT         NN         GG         TT         GG
#> 39          TT         TT         TT         AT         AT         TT
#> 40          CG         GG         CC         GG         CG         CG
#> 41          AG         AG         GG         AG         AG         AG
#> 42          AG         AG         AG         GG         GG         NN
#> 43          AA         AA         AA         AA         AA         AG
#> 44          TT         NN         TT         AA         TT         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          NN         AA         CC         CC         AC         AC
#> 47          CG         CG         CC         CG         CC         CC
#> 48          CT         CC         NN         TT         CT         CT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CT         CC         CC         CC
#> 51          CT         NN         CT         CC         TT         TT
#> 52          AA         AG         AG         AG         AG         GG
#> 53          AA         AA         AC         AA         AA         AC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          AG         AG         AG         AA         NN         GG
#> 56          AG         AA         AG         AA         NN         AA
#> 57          GG         AG         NN         GG         GG         NN
#> 58          NN         CT         CC         CT         CT         CC
#> 59          GG         AG         GG         AA         AG         AG
#> 60          AC         NN         NN         CC         NN         AC
#> 61          TT         NN         TT         TT         TT         TT
#> 62          GG         GG         CG         GG         GG         GG
#> 63          AC         AC         CC         AA         AA         CC
#> 64          AA         AT         TT         TT         TT         TT
#> 65          CT         CC         CC         CC         CC         TT
#> 66          CC         AC         AA         AA         AC         AC
#> 67          GG         GG         NN         GG         GG         GG
#> 68          AT         TT         AA         AA         AA         AT
#> 69          AA         AA         AA         AA         AA         NN
#> 70          AC         AC         AA         AA         AC         CC
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AG         AG         AA         GG         AA         AG
#> 73          GT         NN         GT         GT         GT         GT
#> 74          GT         GG         GG         GG         NN         GG
#> 75          AC         AC         AA         CC         NN         AA
#> 76          TT         NN         NN         CC         NN         CT
#> 77          AT         NN         AT         NN         AA         AA
#> 78          GG         GG         NN         GG         GG         TT
#> 79          CT         TT         TT         TT         TT         TT
#> 80          TT         TT         CC         TT         TT         CT
#> 81          AT         AA         AT         NN         AA         TT
#> 82          NN         TT         TT         TT         TT         TT
#> 83          TT         TT         TT         AT         TT         AT
#> 84          CT         CC         CT         CC         TT         TT
#> 85          AA         AA         AG         AA         AG         AG
#> 86          AG         NN         AG         NN         NN         AA
#> 87          CC         CT         TT         TT         CC         TT
#> 88          AA         AG         AG         AG         NN         AA
#> 89          CG         NN         CC         CG         CC         CC
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CC         CT         CC         NN         CC         CC
#> 92          AG         AG         NN         NN         AG         AG
#> 93          TT         TT         TT         CT         CT         CT
#> 94          CC         CC         CC         NN         NN         CC
#> 95          GT         GT         GG         GT         GG         GT
#> 96          CC         CC         AC         CC         NN         CC
#> 97          GG         GG         GG         GG         GG         GG
#> 98          GG         GT         NN         GG         TT         GG
#> 99          AG         AG         NN         GG         AG         AG
#> 100         TT         TT         TT         AA         AA         AT
#>     CLF.r2.048 CLF.r2.061 CLF.r2.064 CLF.r2.068 CLF.r2.069 OPL.r1.078
#> 1           AG         AG         AA         AG         AA         AG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         CC         CC         CT         CC
#> 4           AG         AA         NN         AA         AA         AG
#> 5           NN         NN         NN         GG         GG         AA
#> 6           AC         AC         AA         AC         AA         AA
#> 7           TT         TT         TT         TT         TT         TT
#> 8           AG         AG         AA         AG         GG         GG
#> 9           AC         CC         AC         CC         CC         CC
#> 10          AA         AA         AA         AG         AA         AA
#> 11          CT         CT         CT         TT         CC         CC
#> 12          AA         AG         AA         AA         AA         AA
#> 13          CT         TT         CT         TT         CC         NN
#> 14          GG         GG         GT         GT         GG         GG
#> 15          AA         AA         AA         AA         AA         NN
#> 16          GT         GG         NN         GG         GT         NN
#> 17          AA         AA         AA         AA         AA         AC
#> 18          AC         AC         NN         AC         AC         CC
#> 19          TT         TT         GT         TT         GT         NN
#> 20          GG         GG         GG         GG         AG         GG
#> 21          AC         AC         NN         AC         CC         NN
#> 22          CC         NN         TT         CT         CT         CC
#> 23          GG         GG         GG         GG         GG         GG
#> 24          AA         AA         AA         AA         AA         AA
#> 25          GG         AG         NN         GG         GG         AG
#> 26          TT         TT         TT         TT         TT         CT
#> 27          AC         AA         AA         AA         AA         AC
#> 28          TT         CT         TT         TT         TT         NN
#> 29          NN         NN         AA         AC         AA         AA
#> 30          GG         CG         GG         CG         NN         GG
#> 31          CG         CC         CC         CC         CC         CG
#> 32          AC         NN         AA         AA         AA         AA
#> 33          CC         CT         NN         CT         NN         NN
#> 34          CG         GG         CG         GG         CC         NN
#> 35          NN         AA         AA         AA         AA         AA
#> 36          AA         AA         AG         AA         AA         NN
#> 37          CG         CG         NN         GG         GG         GG
#> 38          GG         GG         GG         GG         GG         GT
#> 39          NN         NN         AT         NN         AT         NN
#> 40          GG         NN         CG         GG         GG         CG
#> 41          GG         GG         AG         NN         AG         AG
#> 42          GG         AG         GG         AG         AA         GG
#> 43          NN         AA         AA         AG         AA         NN
#> 44          TT         TT         TT         TT         TT         TT
#> 45          GG         GG         GG         GG         GG         NN
#> 46          CC         CC         NN         CC         AC         AC
#> 47          CC         CC         CC         CC         CC         NN
#> 48          CT         TT         CT         CT         CC         CT
#> 49          GG         CG         GG         CG         GG         GG
#> 50          CC         CC         CC         CT         NN         CT
#> 51          NN         TT         TT         TT         TT         CT
#> 52          AG         AA         AA         AA         AG         AG
#> 53          NN         AA         NN         AC         AC         NN
#> 54          GG         GG         GT         GG         GG         NN
#> 55          GG         GG         AG         AG         AG         AA
#> 56          AA         AG         AA         AA         AA         NN
#> 57          GG         GG         GG         GG         GG         AG
#> 58          CC         CC         CC         CC         CC         CT
#> 59          AG         AG         AG         GG         GG         AG
#> 60          AC         CC         CC         CC         AC         AC
#> 61          NN         NN         CT         TT         NN         TT
#> 62          GG         GG         GG         GG         GG         GG
#> 63          AC         AC         AA         NN         AA         NN
#> 64          TT         TT         TT         TT         AA         AT
#> 65          CC         CT         CT         CT         CC         NN
#> 66          AC         CC         AA         AC         AA         NN
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AA         AA         AT         AA         AT         AT
#> 69          AA         AA         AA         AG         AA         GG
#> 70          AA         AC         NN         AC         CC         NN
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AA         AG         GG         GG         AG         NN
#> 73          TT         GT         TT         GT         TT         NN
#> 74          GT         GG         GG         NN         GG         GT
#> 75          AC         AC         CC         AC         NN         AC
#> 76          NN         NN         NN         CC         CC         NN
#> 77          AA         AT         AA         AA         AA         AA
#> 78          GG         GG         GG         GT         GG         NN
#> 79          TT         TT         TT         TT         TT         TT
#> 80          CT         TT         TT         TT         TT         CT
#> 81          NN         AT         AT         AT         AT         AA
#> 82          NN         NN         TT         NN         TT         TT
#> 83          AT         TT         AT         AT         AT         TT
#> 84          CT         CC         CC         CT         CC         CC
#> 85          AA         NN         AG         AG         AG         NN
#> 86          AA         AG         AA         AA         AA         NN
#> 87          CT         TT         TT         TT         TT         CT
#> 88          AA         AG         AA         AG         AG         GG
#> 89          NN         GG         CC         NN         NN         CG
#> 90          CC         CC         CC         CC         CC         CC
#> 91          NN         TT         CC         NN         CT         CT
#> 92          NN         GG         GG         GG         AG         AG
#> 93          TT         CT         TT         CT         TT         TT
#> 94          CC         CC         CC         NN         CC         CG
#> 95          TT         GG         GG         GT         GG         GG
#> 96          CC         CC         AC         CC         CC         AC
#> 97          GG         GG         GG         AG         GG         AG
#> 98          GG         GG         NN         GT         TT         GG
#> 99          AA         GG         AA         GG         GG         NN
#> 100         AT         NN         AA         NN         AT         NN
#>     OPL.r1.080 OPL.r1.083 OPL.r1.086 OPL.r1.087 OPL.r1.088 OPL.r1.109
#> 1           AG         AG         AG         NN         NN         NN
#> 2           NN         GT         GG         GG         GG         GT
#> 3           CT         CT         CC         CC         CT         CC
#> 4           AG         NN         AG         NN         NN         NN
#> 5           NN         NN         AG         AG         NN         GG
#> 6           AA         AC         AA         AA         AA         NN
#> 7           NN         GT         TT         GT         TT         NN
#> 8           GG         GG         GG         AG         GG         AG
#> 9           AC         CC         CC         AC         NN         CC
#> 10          AA         AA         AA         AA         AA         AA
#> 11          NN         NN         CC         CT         CC         CC
#> 12          NN         NN         AA         AG         AG         AG
#> 13          TT         NN         TT         TT         TT         CT
#> 14          GG         GG         GG         GG         GG         GG
#> 15          AG         AA         AA         AA         AA         NN
#> 16          GT         NN         NN         NN         GT         NN
#> 17          AA         AC         NN         AC         NN         AA
#> 18          CC         CC         CC         CC         CC         CC
#> 19          TT         TT         GT         NN         TT         NN
#> 20          GG         NN         GG         GG         GG         GG
#> 21          AC         AA         NN         NN         AC         NN
#> 22          CT         CT         NN         CC         CC         NN
#> 23          GG         AG         GG         GG         GG         GG
#> 24          AA         AA         AG         AA         AA         AA
#> 25          NN         AG         GG         GG         GG         AG
#> 26          TT         TT         NN         NN         TT         TT
#> 27          AC         AA         AC         AA         AA         NN
#> 28          CT         CT         NN         NN         NN         CT
#> 29          AA         NN         AA         NN         AC         NN
#> 30          NN         GG         NN         GG         CG         NN
#> 31          CC         CG         CC         CG         CC         CC
#> 32          AA         AC         NN         NN         AA         AC
#> 33          NN         CC         CT         NN         NN         TT
#> 34          NN         CG         NN         CG         NN         NN
#> 35          AA         AA         NN         NN         AA         AA
#> 36          AG         AA         AA         AA         NN         AG
#> 37          CG         GG         CG         CG         CG         CG
#> 38          GG         NN         NN         NN         GT         GT
#> 39          NN         NN         TT         NN         NN         AT
#> 40          NN         GG         GG         NN         NN         CG
#> 41          NN         NN         NN         AG         NN         AG
#> 42          AG         AG         AG         AG         GG         AG
#> 43          AA         AA         NN         NN         AG         AG
#> 44          TT         TT         NN         TT         TT         TT
#> 45          GG         NN         GG         NN         GG         GG
#> 46          AC         NN         NN         NN         NN         AA
#> 47          CG         CC         CG         CC         CG         CG
#> 48          CC         NN         CT         NN         NN         NN
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CT         NN         NN         CC
#> 51          CC         NN         TT         NN         NN         CT
#> 52          AG         AA         NN         AG         AG         AA
#> 53          AA         AA         AC         AC         NN         AC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          NN         NN         NN         AG         GG         AG
#> 56          AA         NN         NN         AG         AG         AG
#> 57          NN         NN         NN         GG         GG         GG
#> 58          CT         CC         CC         CC         CT         CT
#> 59          AG         AG         AG         AG         AA         GG
#> 60          AC         AC         NN         NN         AC         AC
#> 61          NN         NN         TT         TT         NN         TT
#> 62          GG         CG         GG         CG         GG         GG
#> 63          AC         NN         AA         AC         AA         AA
#> 64          NN         AT         NN         AT         TT         AT
#> 65          NN         TT         NN         CT         TT         CT
#> 66          AC         AC         AC         AA         AC         NN
#> 67          NN         GG         GG         GG         GG         GG
#> 68          AT         AA         AT         NN         AA         NN
#> 69          AA         AA         NN         AG         AG         AA
#> 70          AA         NN         NN         AC         AC         AC
#> 71          AT         AA         AA         AA         AA         AA
#> 72          GG         AG         AG         AG         NN         AA
#> 73          TT         TT         TT         TT         TT         TT
#> 74          GT         GG         NN         GT         NN         GG
#> 75          AC         AC         AC         NN         AC         CC
#> 76          CT         NN         CC         NN         NN         NN
#> 77          AA         AA         NN         AT         AT         AT
#> 78          GG         NN         NN         GT         GT         NN
#> 79          TT         TT         TT         TT         TT         CT
#> 80          CT         CT         NN         NN         CT         TT
#> 81          NN         AA         NN         NN         NN         AT
#> 82          TT         TT         TT         CT         NN         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          CT         NN         NN         CC         CC         NN
#> 85          GG         NN         AG         NN         AG         GG
#> 86          AA         AA         NN         AG         AG         AG
#> 87          NN         TT         CT         CT         TT         TT
#> 88          NN         NN         NN         AG         AG         AA
#> 89          NN         CC         CG         NN         CG         NN
#> 90          CC         CC         NN         CC         CC         CC
#> 91          NN         CT         CT         CT         NN         NN
#> 92          NN         AG         AG         NN         AG         NN
#> 93          TT         TT         TT         CT         NN         TT
#> 94          CC         CC         NN         CC         CC         CG
#> 95          NN         NN         GG         GT         NN         NN
#> 96          CC         AC         NN         NN         NN         CC
#> 97          GG         GG         NN         NN         GG         NN
#> 98          GG         GT         GG         GG         GG         GT
#> 99          NN         AG         AA         AG         AG         AG
#> 100         NN         AT         AT         AT         AT         NN
#>     OPL.r1.110 OPL.r1.116 OPL.r1.125 OPL.r1.126 OPL.r1.128 OPL.r1.130
#> 1           AG         AG         NN         AG         AG         AA
#> 2           GG         GG         GG         GG         GG         GT
#> 3           NN         CT         NN         CC         CC         CC
#> 4           AG         NN         AG         AG         AG         AA
#> 5           GG         AG         GG         NN         GG         NN
#> 6           AA         AA         AA         AA         AA         AA
#> 7           NN         TT         TT         TT         TT         TT
#> 8           NN         GG         AG         NN         AG         GG
#> 9           CC         CC         CC         NN         CC         CC
#> 10          AA         NN         NN         AA         AA         AA
#> 11          NN         NN         CC         NN         CT         TT
#> 12          NN         NN         AA         AA         AA         NN
#> 13          NN         TT         TT         TT         TT         CC
#> 14          TT         GG         NN         GG         GG         GG
#> 15          AA         AG         NN         NN         AA         AG
#> 16          GT         GG         GT         TT         GG         GG
#> 17          AC         AC         AA         AA         AC         AA
#> 18          AC         CC         CC         AC         AC         AA
#> 19          TT         NN         TT         TT         NN         NN
#> 20          GG         NN         GG         NN         AG         GG
#> 21          AA         AC         AC         AC         AC         AC
#> 22          CC         CT         NN         NN         CC         CT
#> 23          GG         AG         GG         GG         GG         GG
#> 24          AA         AA         NN         AA         AG         AA
#> 25          AG         AG         NN         AG         AG         GG
#> 26          TT         NN         TT         NN         TT         TT
#> 27          AC         AA         NN         AC         AA         AA
#> 28          CT         CC         TT         CT         NN         CC
#> 29          NN         AA         NN         AA         AA         AC
#> 30          CG         GG         CG         NN         GG         GG
#> 31          NN         CG         NN         CG         CG         CC
#> 32          AA         AC         AC         AC         AC         NN
#> 33          TT         CC         CC         CT         CT         CT
#> 34          CG         CG         NN         CG         CG         CC
#> 35          AC         AA         AA         AA         AC         AA
#> 36          AA         NN         AA         AA         AA         AA
#> 37          GG         CG         NN         GG         CG         GG
#> 38          NN         GG         GT         GT         GG         GG
#> 39          NN         NN         NN         AT         TT         TT
#> 40          NN         CG         CG         GG         CG         NN
#> 41          NN         NN         GG         GG         AG         GG
#> 42          AA         NN         GG         AG         GG         AA
#> 43          NN         NN         AG         AA         AA         AA
#> 44          AT         TT         TT         TT         TT         TT
#> 45          GG         GG         GG         GG         GG         NN
#> 46          CC         CC         NN         AC         CC         AC
#> 47          CG         CG         CC         CG         CG         CC
#> 48          NN         CC         TT         CT         CT         CT
#> 49          GG         CG         GG         GG         GG         GG
#> 50          CC         CC         CC         NN         CC         CC
#> 51          CT         CT         CT         TT         TT         TT
#> 52          AA         AG         AA         AG         AG         AA
#> 53          AA         NN         AA         AC         AC         CC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          AA         GG         AG         AG         AG         GG
#> 56          AG         AA         AA         AA         AA         NN
#> 57          NN         GG         GG         GG         GG         AG
#> 58          CT         CC         CC         CC         CC         CT
#> 59          GG         AA         AG         NN         AG         AG
#> 60          AC         AA         AA         AC         AA         CC
#> 61          TT         TT         NN         NN         TT         TT
#> 62          GG         CG         NN         GG         CG         CG
#> 63          NN         AC         AA         AC         AC         AA
#> 64          NN         AT         AT         AT         TT         TT
#> 65          NN         NN         CT         CC         TT         CT
#> 66          AA         AA         NN         NN         AA         AA
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AT         AA         AT         AA         AT         AA
#> 69          AA         AA         AA         AA         AA         AG
#> 70          AA         AA         AA         NN         AA         AC
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AG         AA         NN         AA         AA         AG
#> 73          NN         GT         TT         TT         TT         GT
#> 74          NN         GT         NN         GT         GG         GT
#> 75          NN         NN         AC         AA         AC         AA
#> 76          TT         CT         NN         CC         TT         TT
#> 77          NN         AA         AA         AA         AA         AA
#> 78          NN         GG         GG         GG         GG         GT
#> 79          NN         CT         TT         TT         TT         TT
#> 80          CT         NN         TT         CT         TT         TT
#> 81          AT         AA         AT         AT         AT         TT
#> 82          TT         CT         TT         TT         TT         CT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          CT         NN         CT         NN         CC         CT
#> 85          AG         AG         AG         GG         AG         NN
#> 86          AG         NN         AA         NN         AA         AA
#> 87          TT         NN         TT         TT         TT         TT
#> 88          AG         NN         AA         AA         AA         NN
#> 89          NN         NN         NN         CG         CG         GG
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CT         CC         CT         CT         NN         CC
#> 92          GG         NN         NN         NN         GG         NN
#> 93          TT         CT         TT         TT         TT         TT
#> 94          CC         CC         CC         CC         CC         CC
#> 95          GT         NN         GG         GG         GG         GG
#> 96          CC         CC         CC         CC         AC         CC
#> 97          GG         NN         GG         GG         GG         GG
#> 98          GT         GG         NN         GG         GG         GG
#> 99          NN         NN         AG         AA         AA         AG
#> 100         NN         NN         AA         NN         AT         AA
#>     OPL.r1.139 OPL.r1.146 OPL.r1.147 OPL.r1.149 OPL.r2.007 OPL.r2.012
#> 1           AA         NN         AG         GG         AA         AG
#> 2           GG         NN         GG         NN         GG         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AG         AG         AG         AG         AG         AG
#> 5           NN         AG         NN         NN         GG         AG
#> 6           AA         AA         AC         AA         AA         AA
#> 7           NN         TT         TT         NN         GT         TT
#> 8           GG         AG         NN         AG         GG         AA
#> 9           CC         CC         CC         CC         AC         AC
#> 10          AA         AA         AA         AA         AA         AA
#> 11          CT         CT         CT         CC         TT         TT
#> 12          AG         AG         AA         AG         AG         AA
#> 13          TT         CT         CT         CC         TT         CT
#> 14          NN         GG         GG         GG         GG         GG
#> 15          AA         AG         AA         NN         NN         AA
#> 16          GG         GT         GT         GG         GT         GT
#> 17          AA         AA         AA         AA         AA         AA
#> 18          AC         CC         CC         CC         CC         CC
#> 19          GT         TT         TT         GT         TT         TT
#> 20          GG         GG         GG         GG         GG         GG
#> 21          NN         CC         CC         AC         AC         NN
#> 22          CC         NN         TT         CT         CC         CT
#> 23          GG         GG         GG         GG         GG         GG
#> 24          AA         AG         AA         AA         AA         AA
#> 25          GG         AG         GG         GG         GG         GG
#> 26          TT         TT         TT         CT         TT         NN
#> 27          AC         AC         AC         AA         NN         AA
#> 28          CC         CC         TT         CT         CC         NN
#> 29          AA         AA         AA         AC         AC         AC
#> 30          GG         GG         GG         GG         GG         GG
#> 31          CG         GG         CC         CC         CG         CG
#> 32          NN         AA         AA         AA         AA         AA
#> 33          CC         CC         NN         CT         NN         NN
#> 34          CC         NN         CC         CC         CG         CG
#> 35          NN         NN         AA         AC         AA         AA
#> 36          AG         AA         AA         AG         AA         AG
#> 37          GG         GG         CG         GG         GG         GG
#> 38          GG         GT         GG         GG         GG         GG
#> 39          TT         NN         TT         NN         AT         TT
#> 40          NN         GG         GG         GG         GG         GG
#> 41          AG         GG         GG         GG         GG         NN
#> 42          NN         AG         AG         AA         AA         AG
#> 43          AA         AA         AA         AA         NN         AA
#> 44          TT         TT         TT         TT         TT         TT
#> 45          GG         GG         NN         GG         GG         NN
#> 46          NN         NN         AC         CC         NN         AC
#> 47          CG         NN         GG         CG         NN         CG
#> 48          TT         CT         CT         CC         CC         CT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CC         CT         CC         CC
#> 51          NN         CT         CT         TT         CC         NN
#> 52          NN         AG         AA         AA         AG         AG
#> 53          AC         AA         AC         AC         AC         AA
#> 54          GG         GG         GG         GG         GG         GG
#> 55          NN         GG         GG         GG         AG         NN
#> 56          AG         AG         AA         AG         AG         AA
#> 57          AG         AG         GG         GG         GG         GG
#> 58          CT         CT         NN         CC         CC         CT
#> 59          NN         AG         NN         AG         AG         GG
#> 60          AA         NN         AC         AC         AC         AA
#> 61          CT         TT         TT         TT         NN         TT
#> 62          GG         CG         CG         GG         CG         CG
#> 63          NN         AA         AC         NN         AA         NN
#> 64          TT         NN         AT         TT         AT         NN
#> 65          TT         NN         CT         TT         TT         NN
#> 66          AC         AA         AC         AC         AA         AC
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AA         AT         AT         AA         AA         AA
#> 69          AG         AG         AA         AG         GG         AA
#> 70          AC         AC         AA         AC         NN         AA
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AA         AG         AA         AG         AA         AA
#> 73          TT         GT         TT         NN         TT         GT
#> 74          GG         NN         GT         GG         GT         NN
#> 75          AA         AC         AA         NN         AA         AA
#> 76          NN         CT         NN         CC         NN         NN
#> 77          AT         AT         AA         AT         AT         NN
#> 78          GT         GT         GG         GT         TT         GG
#> 79          TT         TT         TT         TT         TT         TT
#> 80          NN         NN         TT         CC         TT         TT
#> 81          AT         AA         AT         NN         AT         AT
#> 82          TT         TT         TT         TT         CT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          CC         NN         TT         CT         CC         CT
#> 85          AG         AG         AG         GG         AG         AG
#> 86          AG         AG         AA         AG         AG         AA
#> 87          TT         CT         CT         TT         CT         TT
#> 88          AG         AA         AA         AA         AG         NN
#> 89          CC         NN         CC         CG         CC         CC
#> 90          CC         CC         CC         CC         CC         CC
#> 91          TT         TT         CT         NN         TT         CT
#> 92          NN         AG         AG         AG         GG         NN
#> 93          TT         TT         TT         TT         TT         TT
#> 94          CC         NN         CC         CC         NN         CC
#> 95          GG         GG         GG         GG         TT         GG
#> 96          CC         CC         CC         CC         AC         CC
#> 97          GG         GG         AG         NN         GG         GG
#> 98          GT         GG         GG         GG         GG         GG
#> 99          GG         NN         AA         AG         AG         AG
#> 100         NN         TT         AA         AT         AA         AT
#>     OPL.r2.022 OPL.r2.023 OPL.r2.028 OPL.r2.031 OPL.r2.033 OPL.r2.079
#> 1           AA         NN         NN         AG         AG         AA
#> 2           GG         NN         GG         GG         GG         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           NN         AA         AG         AG         NN         AG
#> 5           GG         NN         NN         NN         NN         NN
#> 6           AA         AC         AA         AA         AA         AA
#> 7           GT         GT         TT         GT         TT         TT
#> 8           GG         GG         GG         GG         GG         GG
#> 9           CC         CC         CC         CC         AC         CC
#> 10          AA         AA         AG         AA         AA         AA
#> 11          CT         CT         CC         CT         CT         NN
#> 12          AA         NN         NN         AA         NN         AA
#> 13          CT         TT         CT         TT         TT         CT
#> 14          GG         GG         GG         GG         GG         GG
#> 15          AA         AG         GG         NN         AA         AG
#> 16          GG         GT         NN         GG         GT         GT
#> 17          AA         AA         AA         AC         AA         AA
#> 18          AC         CC         CC         AC         CC         AC
#> 19          TT         TT         NN         TT         TT         TT
#> 20          GG         GG         GG         GG         GG         GG
#> 21          CC         AC         AC         NN         AC         NN
#> 22          CT         TT         CC         CT         CC         CT
#> 23          GG         GG         GG         GG         GG         NN
#> 24          AA         AA         AA         AA         AA         AA
#> 25          GG         AG         NN         AG         GG         NN
#> 26          TT         CT         CT         TT         TT         CT
#> 27          AA         AC         AA         AA         AA         AC
#> 28          CT         NN         CC         CT         NN         CC
#> 29          AA         AA         CC         CC         NN         AA
#> 30          GG         GG         GG         GG         NN         NN
#> 31          CG         CC         CC         CC         CC         CG
#> 32          AA         AA         AA         AA         AA         AA
#> 33          CC         TT         CC         NN         CT         CC
#> 34          CG         CG         GG         GG         CC         NN
#> 35          AA         AA         AA         NN         AA         AA
#> 36          GG         AG         AA         AA         AA         AA
#> 37          CG         CG         NN         GG         CG         CG
#> 38          GT         GT         GG         GG         GG         GT
#> 39          TT         NN         NN         TT         TT         TT
#> 40          GG         GG         CG         GG         NN         GG
#> 41          AG         AG         GG         AG         AG         NN
#> 42          AA         GG         GG         NN         AG         AG
#> 43          AG         AA         AG         AG         AA         AG
#> 44          TT         TT         AT         TT         TT         TT
#> 45          GG         GG         NN         GG         GG         NN
#> 46          AC         AC         AA         CC         NN         AC
#> 47          GG         CG         GG         CG         GG         CG
#> 48          CT         TT         CC         CC         CT         NN
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CT         CC         CC         NN         CC
#> 51          NN         TT         TT         NN         CT         CT
#> 52          AG         AA         AA         GG         AG         NN
#> 53          AA         NN         NN         AC         NN         AC
#> 54          GT         GG         GG         GG         GG         GG
#> 55          GG         GG         AG         AG         AA         NN
#> 56          AA         AA         AA         AA         NN         AA
#> 57          GG         GG         NN         GG         NN         NN
#> 58          CT         NN         CC         CC         CT         CT
#> 59          AG         AG         AG         GG         GG         GG
#> 60          AA         AC         NN         AC         AA         AA
#> 61          NN         TT         NN         NN         TT         TT
#> 62          GG         GG         GG         CG         CG         GG
#> 63          AA         NN         AC         AA         NN         AC
#> 64          AT         AT         TT         AT         TT         NN
#> 65          CT         CT         NN         TT         CC         CT
#> 66          AC         NN         NN         AC         NN         CC
#> 67          GG         GG         GG         GG         NN         GG
#> 68          AT         NN         TT         AA         AA         AA
#> 69          AG         AG         NN         AG         AA         AA
#> 70          AC         AC         AC         AC         AA         AA
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AG         AG         AG         AA         AA         AG
#> 73          TT         TT         TT         TT         GT         NN
#> 74          NN         GG         NN         GT         GG         GT
#> 75          AC         AC         AC         CC         AC         AC
#> 76          CC         TT         NN         CC         CC         CC
#> 77          AA         AA         AA         AA         AT         NN
#> 78          GT         GT         GG         GT         GG         GG
#> 79          CT         CT         TT         CT         TT         CT
#> 80          TT         CT         TT         CT         CT         CT
#> 81          AT         AA         TT         AT         AA         AT
#> 82          CT         TT         TT         TT         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          TT         CT         CT         CC         NN         CT
#> 85          AG         GG         GG         GG         GG         NN
#> 86          NN         NN         AA         NN         AG         AA
#> 87          TT         CC         TT         TT         CC         TT
#> 88          AG         GG         AA         AG         NN         AG
#> 89          CC         CG         GG         CC         CG         NN
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CT         TT         CT         CC         CT         CT
#> 92          AG         NN         NN         NN         NN         AG
#> 93          TT         TT         CT         TT         TT         TT
#> 94          CC         NN         CC         CC         CC         NN
#> 95          GG         GG         GG         GT         GG         TT
#> 96          AC         AC         AC         CC         NN         AC
#> 97          GG         GG         NN         GG         GG         GG
#> 98          GG         GT         GG         GG         GG         GG
#> 99          AG         AG         AA         GG         NN         AG
#> 100         AA         AA         AA         TT         AA         AA
#>     OPL.r2.080 OPL.r2.094 PAL.r2.203 PAL.r2.237 PAL.r2.238 PAL.r2.239
#> 1           GG         AG         AG         GG         GG         GG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AG         AA         AA         AA         AA         AA
#> 5           AG         NN         AG         NN         GG         GG
#> 6           AA         AA         AA         AC         AC         NN
#> 7           TT         NN         TT         TT         TT         TT
#> 8           GG         GG         NN         GG         GG         GG
#> 9           CC         CC         CC         CC         CC         AC
#> 10          AA         AA         AA         AA         AA         AA
#> 11          CT         CT         TT         TT         TT         TT
#> 12          AA         NN         AG         AG         AG         AA
#> 13          TT         CT         CT         TT         TT         TT
#> 14          GG         TT         GG         GG         GG         GG
#> 15          AA         AA         GG         GG         AG         AG
#> 16          GG         NN         GT         NN         TT         NN
#> 17          AA         AA         CC         AA         AA         AA
#> 18          AC         AC         NN         AA         AC         AA
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         GG         NN         GG         GG         GG
#> 21          AC         AA         CC         CC         CC         CC
#> 22          CC         NN         CC         CC         CC         CC
#> 23          GG         NN         GG         GG         GG         GG
#> 24          GG         AA         AA         AA         AA         AA
#> 25          AG         AG         GG         GG         GG         GG
#> 26          CT         NN         TT         CT         CT         CC
#> 27          AA         AA         AA         AA         AA         AA
#> 28          CT         CC         TT         TT         TT         TT
#> 29          AC         AC         AC         AA         AC         AC
#> 30          GG         GG         GG         GG         GG         GG
#> 31          GG         CG         CC         CC         CC         CC
#> 32          AA         AA         AA         AC         AA         AC
#> 33          CT         CT         CC         CT         CC         NN
#> 34          NN         CC         CC         CC         CC         CC
#> 35          AA         AA         AA         AA         AA         AA
#> 36          AA         AG         AA         AA         AA         AA
#> 37          CG         NN         GG         GG         GG         GG
#> 38          GG         GG         GG         GT         GT         GT
#> 39          TT         TT         TT         AT         TT         NN
#> 40          GG         CG         GG         CG         GG         CG
#> 41          GG         AG         AA         AA         AA         AA
#> 42          NN         GG         AG         GG         GG         GG
#> 43          AA         AG         AA         AA         AA         AA
#> 44          TT         AT         AT         AT         AA         AT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          CC         AA         CC         CC         CC         CC
#> 47          CC         NN         CC         CC         CC         CC
#> 48          CT         NN         CT         CC         CT         CT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CT         CC         CT         CC         CC
#> 51          TT         CT         NN         TT         CT         CT
#> 52          GG         AA         AG         AG         NN         AA
#> 53          AA         AA         AA         AA         AA         AA
#> 54          GG         GG         GG         GG         GG         GG
#> 55          GG         GG         AA         AG         AG         AG
#> 56          NN         NN         AG         AG         AG         AA
#> 57          GG         GG         NN         NN         GG         GG
#> 58          CT         CC         CC         TT         TT         CT
#> 59          GG         AG         GG         GG         GG         GG
#> 60          AC         AC         AA         AA         AA         AA
#> 61          TT         NN         NN         TT         TT         NN
#> 62          GG         GG         GG         GG         GG         GG
#> 63          AA         NN         NN         NN         AC         CC
#> 64          TT         AT         AT         AA         NN         AT
#> 65          TT         CC         CC         CC         CT         CC
#> 66          AC         AC         CC         AC         AC         AA
#> 67          GG         GG         AG         AA         GG         NN
#> 68          AA         AA         TT         TT         NN         NN
#> 69          GG         AA         AA         AA         AA         AA
#> 70          NN         AA         AA         AA         AA         AA
#> 71          AA         AA         AA         AT         AA         AT
#> 72          AG         AG         GG         AG         AG         GG
#> 73          TT         GT         TT         TT         GT         GT
#> 74          GT         GT         GG         GG         GG         GG
#> 75          CC         AC         AA         AC         AC         AA
#> 76          CC         CC         NN         NN         NN         TT
#> 77          AA         AA         AT         AT         AT         AA
#> 78          TT         GG         GG         GG         GG         GG
#> 79          TT         TT         TT         TT         TT         CT
#> 80          CT         CT         TT         CC         CT         NN
#> 81          AT         NN         TT         NN         AT         AT
#> 82          TT         TT         TT         TT         TT         TT
#> 83          TT         NN         TT         TT         TT         TT
#> 84          CC         CT         NN         NN         TT         TT
#> 85          AG         NN         AA         AA         AA         AA
#> 86          AA         AA         AG         NN         AG         NN
#> 87          CT         CT         TT         TT         NN         TT
#> 88          AG         AA         GG         NN         GG         GG
#> 89          CG         NN         CC         NN         CC         NN
#> 90          CC         CC         CC         CC         CC         CC
#> 91          NN         CC         CC         CC         CC         CC
#> 92          AG         GG         AG         NN         GG         AG
#> 93          TT         TT         NN         CT         TT         TT
#> 94          CG         CC         CG         CG         CC         NN
#> 95          GT         GG         GG         GG         GG         GG
#> 96          CC         CC         CC         CC         CC         CC
#> 97          GG         GG         GG         GG         GG         GG
#> 98          GG         GG         GG         GG         GG         GT
#> 99          AG         AG         GG         GG         GG         GG
#> 100         AT         AA         NN         NN         AT         NN
#>     PAL.r2.240 PAL.r2.241 PAL.r2.242 PAL.r2.249 PAL.r2.253 PAL.r2.275
#> 1           GG         GG         AG         GG         GG         GG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AA         AA         AA         NN         AA         AA
#> 5           GG         GG         AG         NN         AG         NN
#> 6           AC         AC         AC         AC         AA         AC
#> 7           TT         TT         TT         TT         TT         TT
#> 8           GG         GG         GG         GG         GG         NN
#> 9           CC         CC         CC         CC         CC         CC
#> 10          AA         AA         AA         AA         AA         AA
#> 11          TT         TT         TT         NN         TT         TT
#> 12          AG         AG         AG         NN         AA         AG
#> 13          TT         NN         TT         TT         TT         NN
#> 14          GG         GG         GG         GG         GG         GG
#> 15          GG         GG         GG         GG         NN         GG
#> 16          GT         GT         GT         GG         TT         GG
#> 17          AA         AC         AA         AA         NN         AC
#> 18          AA         AC         AA         AC         AC         CC
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         GG         GG         GG         GG         GG
#> 21          CC         CC         CC         NN         CC         CC
#> 22          CC         CC         CC         CC         CC         CC
#> 23          GG         GG         NN         GG         GG         GG
#> 24          AA         AA         AA         AA         NN         AA
#> 25          GG         GG         GG         GG         GG         GG
#> 26          CC         CT         TT         CT         CC         TT
#> 27          AA         AA         AA         AA         AA         AA
#> 28          TT         TT         CT         TT         CT         TT
#> 29          AA         AA         AA         AA         AA         AC
#> 30          GG         GG         GG         GG         GG         GG
#> 31          CC         CC         CC         CC         CC         CC
#> 32          AC         AC         AC         AC         AA         AC
#> 33          CC         CT         CC         CT         CC         CC
#> 34          CC         CC         CC         CC         CC         CC
#> 35          AA         AA         AA         AA         AA         AA
#> 36          AA         AA         AA         AA         AA         AA
#> 37          NN         GG         GG         CG         GG         GG
#> 38          GT         GT         GT         GG         GT         GT
#> 39          NN         TT         AT         NN         NN         TT
#> 40          GG         CG         GG         NN         GG         CG
#> 41          AA         AA         NN         AA         AA         NN
#> 42          GG         GG         GG         NN         GG         GG
#> 43          AA         AA         AA         AA         AA         AA
#> 44          AT         AT         AT         TT         AT         AA
#> 45          GG         GG         GG         GG         GG         GG
#> 46          CC         CC         CC         CC         CC         CC
#> 47          CC         NN         CG         CG         NN         CG
#> 48          TT         CT         TT         NN         CC         TT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CC         CC         CC         CC
#> 51          NN         NN         TT         TT         TT         CT
#> 52          GG         AA         AG         NN         AG         NN
#> 53          AA         AA         AA         AA         AA         AA
#> 54          GG         GG         GG         GG         GG         GG
#> 55          NN         AG         GG         AG         GG         AG
#> 56          AG         AG         AG         AG         AA         AG
#> 57          GG         GG         NN         GG         GG         GG
#> 58          NN         CT         CT         CT         CT         NN
#> 59          GG         GG         GG         GG         GG         GG
#> 60          AA         AA         AA         AA         AA         AA
#> 61          NN         TT         TT         TT         NN         TT
#> 62          GG         GG         GG         NN         GG         GG
#> 63          AC         CC         NN         NN         NN         AC
#> 64          AA         AT         AA         AT         AA         TT
#> 65          CC         CC         CC         CC         CT         CT
#> 66          AA         CC         AC         CC         CC         CC
#> 67          GG         AG         GG         AG         GG         AG
#> 68          TT         AT         AT         AA         TT         TT
#> 69          AA         AA         AA         AA         AA         AA
#> 70          AA         NN         AA         AA         AA         AA
#> 71          NN         AA         AA         AT         AA         AA
#> 72          GG         AG         AA         NN         GG         GG
#> 73          TT         GT         GT         TT         NN         GT
#> 74          GG         GG         GG         GG         GG         GG
#> 75          AC         AC         AA         AA         AA         AC
#> 76          TT         TT         TT         NN         NN         NN
#> 77          AT         AT         AT         AT         AA         AT
#> 78          GG         GG         GG         GG         GG         GG
#> 79          CT         TT         TT         CT         CT         TT
#> 80          CT         CT         CT         TT         CT         CT
#> 81          AA         AA         AA         AT         AT         AT
#> 82          TT         TT         TT         TT         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          TT         TT         TT         TT         TT         TT
#> 85          AA         AA         NN         AG         AA         AA
#> 86          NN         AG         AG         AG         NN         AG
#> 87          TT         TT         TT         TT         TT         TT
#> 88          GG         GG         GG         GG         GG         GG
#> 89          GG         CC         CG         GG         NN         CC
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CC         CC         CC         CC         CT         CC
#> 92          NN         AA         AA         NN         AG         NN
#> 93          TT         CT         TT         TT         TT         TT
#> 94          CG         CC         NN         NN         CG         CG
#> 95          GG         GG         GG         GG         GG         GG
#> 96          CC         CC         CC         CC         CC         CC
#> 97          GG         GG         GG         GG         GG         GG
#> 98          GT         GT         GG         NN         GT         GT
#> 99          GG         GG         NN         GG         GG         GG
#> 100         TT         AA         NN         AA         NN         AT
#>     PAL.r2.276 PAL.r2.278 PAL.r2.281 PAL.r2.285 PAL.r2.287 PAL.r2.290
#> 1           AG         GG         AG         GG         GG         AG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         NN         CT         CC         CC
#> 4           AA         AA         AA         AA         AA         AA
#> 5           AG         GG         NN         AG         NN         GG
#> 6           AC         AC         AC         CC         AA         AA
#> 7           TT         TT         TT         TT         TT         TT
#> 8           GG         GG         GG         NN         GG         GG
#> 9           CC         CC         CC         CC         CC         NN
#> 10          AA         AA         AA         AA         AA         AA
#> 11          TT         TT         TT         TT         TT         TT
#> 12          AG         AA         AA         AA         AG         AA
#> 13          TT         TT         TT         TT         TT         TT
#> 14          GG         GG         GG         GG         GG         GG
#> 15          AG         GG         AG         GG         GG         GG
#> 16          GT         GT         GG         GT         GT         GT
#> 17          AA         AA         AA         AA         AC         AC
#> 18          NN         AA         NN         AA         AC         AA
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         NN         GG         GG         GG         GG
#> 21          CC         CC         CC         CC         CC         CC
#> 22          CC         CC         CC         CC         CC         CC
#> 23          GG         GG         GG         GG         GG         GG
#> 24          AA         AA         AA         AA         AA         AA
#> 25          NN         GG         GG         GG         GG         GG
#> 26          CC         NN         CT         CT         CT         NN
#> 27          AA         AA         AA         AA         AA         AA
#> 28          TT         TT         NN         NN         TT         CT
#> 29          AA         AA         NN         NN         AA         AC
#> 30          GG         GG         GG         GG         GG         GG
#> 31          CC         CC         CC         CC         CC         CC
#> 32          AC         AC         AC         NN         AC         NN
#> 33          CT         CT         CT         CC         CC         CT
#> 34          CC         CC         CC         CC         CC         CC
#> 35          AA         NN         AA         AA         AA         AA
#> 36          AA         AA         AA         AA         AA         NN
#> 37          GG         GG         GG         GG         GG         GG
#> 38          GG         GT         GT         GG         NN         GT
#> 39          TT         NN         TT         NN         TT         AT
#> 40          CG         GG         NN         GG         GG         GG
#> 41          NN         AA         NN         AA         AA         AA
#> 42          AG         GG         GG         GG         GG         GG
#> 43          AA         AA         AA         AA         AA         AA
#> 44          NN         TT         TT         NN         AT         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          CC         CC         CC         CC         CC         CC
#> 47          CC         CC         CC         CC         CC         CC
#> 48          CT         NN         NN         TT         TT         TT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         NN         CC         CC         CC
#> 51          TT         NN         NN         TT         NN         TT
#> 52          NN         AG         AG         AG         NN         AG
#> 53          AA         AA         AA         NN         AA         AA
#> 54          GG         GG         GG         GG         GG         NN
#> 55          GG         AG         AG         NN         AG         AG
#> 56          AG         AA         NN         AA         AG         AA
#> 57          GG         GG         GG         GG         NN         NN
#> 58          CC         CC         CC         NN         CC         NN
#> 59          GG         GG         GG         AG         GG         AG
#> 60          AA         AA         AA         AA         AA         AA
#> 61          TT         TT         TT         TT         TT         TT
#> 62          GG         GG         GG         GG         GG         GG
#> 63          CC         CC         NN         NN         CC         AC
#> 64          AA         TT         AT         TT         NN         NN
#> 65          CC         CT         CT         TT         CT         CC
#> 66          CC         AC         NN         AC         AC         CC
#> 67          GG         NN         AG         GG         GG         AG
#> 68          NN         AT         AT         TT         AT         NN
#> 69          AA         AA         AA         AA         AA         AA
#> 70          AA         AA         AA         AA         AA         AA
#> 71          AA         AA         AT         AA         AA         AT
#> 72          NN         GG         AA         AG         AG         AG
#> 73          TT         GT         GG         TT         NN         GT
#> 74          GG         GG         GG         GG         GG         GG
#> 75          AA         AC         NN         AC         AC         AC
#> 76          NN         TT         NN         TT         TT         NN
#> 77          AT         NN         AA         NN         AT         AA
#> 78          GG         GG         GG         GG         GG         GG
#> 79          TT         TT         CT         CT         CT         TT
#> 80          CT         CT         CC         CT         CT         NN
#> 81          AT         AT         TT         AT         AA         AT
#> 82          TT         TT         TT         TT         TT         TT
#> 83          TT         TT         TT         TT         NN         TT
#> 84          TT         TT         TT         TT         TT         TT
#> 85          AG         NN         AA         AA         AA         AG
#> 86          AG         NN         AA         NN         AG         AA
#> 87          TT         TT         TT         TT         TT         TT
#> 88          GG         GG         GG         GG         GG         NN
#> 89          CG         NN         CC         GG         NN         GG
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CC         CT         CC         CT         CC         CC
#> 92          AA         AG         NN         AG         AG         NN
#> 93          NN         TT         CT         TT         TT         CT
#> 94          CC         GG         CG         NN         CG         CC
#> 95          GG         GG         GG         GG         GG         GG
#> 96          CC         CC         CC         CC         CC         CC
#> 97          GG         GG         GG         GG         GG         GG
#> 98          GT         GG         GT         GT         TT         GT
#> 99          NN         GG         GG         GG         GG         GG
#> 100         AT         AT         AA         AT         AA         NN
#>     PAL.r2.291 SMR.r1.003 SMR.r1.004 SMR.r1.008 SMR.r1.012 SMR.r1.013
#> 1           AA         AA         GG         AA         AG         NN
#> 2           GG         GG         GT         GG         GG         GG
#> 3           CC         CC         CC         CT         CC         CC
#> 4           AA         AG         AA         AA         AG         NN
#> 5           GG         NN         NN         NN         AG         GG
#> 6           AC         AC         AA         AC         NN         AA
#> 7           TT         NN         TT         TT         NN         GT
#> 8           GG         AG         AG         GG         NN         GG
#> 9           CC         CC         AC         CC         NN         AC
#> 10          AA         AA         AA         AA         AA         NN
#> 11          TT         NN         CT         CC         CT         TT
#> 12          AG         NN         AA         AG         AA         AA
#> 13          TT         TT         TT         TT         TT         TT
#> 14          GG         GG         GG         GG         NN         GG
#> 15          GG         AA         AA         AG         AG         AA
#> 16          GG         GG         GG         GT         GG         NN
#> 17          AA         AC         AA         AC         AC         NN
#> 18          AC         AC         CC         CC         CC         NN
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         GG         GG         AG         AG         GG
#> 21          CC         AC         AC         CC         AC         AC
#> 22          CC         CC         CC         CT         CC         NN
#> 23          GG         GG         AG         GG         GG         GG
#> 24          AA         NN         AA         AG         AA         AA
#> 25          GG         AG         AA         AG         AA         AG
#> 26          CT         CT         NN         TT         TT         TT
#> 27          AA         AA         NN         AA         AA         AC
#> 28          TT         TT         TT         CC         CC         CT
#> 29          AA         NN         NN         AC         AA         AC
#> 30          GG         GG         GG         GG         CG         CG
#> 31          CC         NN         CC         GG         CG         NN
#> 32          AC         AA         AC         AC         AA         AA
#> 33          NN         CC         CT         CT         CC         CC
#> 34          CC         CC         CG         CC         CC         CG
#> 35          AA         NN         AA         AA         AA         NN
#> 36          AA         AA         AA         AG         AA         AG
#> 37          GG         CC         GG         CG         GG         CG
#> 38          TT         GT         GT         TT         GT         NN
#> 39          TT         NN         NN         AT         TT         TT
#> 40          GG         GG         GG         GG         GG         GG
#> 41          AA         AG         AG         AA         NN         AG
#> 42          GG         AG         AA         AG         GG         AG
#> 43          AA         AA         AA         AA         AA         AA
#> 44          TT         NN         TT         TT         NN         AA
#> 45          GG         GG         GG         GG         GG         GG
#> 46          NN         CC         NN         AA         CC         CC
#> 47          CC         CC         CC         CG         CC         CC
#> 48          CC         TT         TT         CC         CT         CT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CC         CC         CC         CT
#> 51          CT         NN         NN         CT         CT         NN
#> 52          AG         AA         AG         AA         AG         AA
#> 53          AA         AA         AA         AA         AA         AA
#> 54          GG         GG         GG         GG         GG         GG
#> 55          AA         GG         GG         AA         GG         AG
#> 56          AG         AA         NN         AG         AA         AA
#> 57          GG         NN         GG         GG         GG         GG
#> 58          TT         NN         CC         CC         CC         CC
#> 59          NN         NN         AG         NN         NN         AG
#> 60          AA         AA         NN         AC         AC         NN
#> 61          TT         TT         TT         TT         TT         TT
#> 62          GG         GG         GG         CG         GG         CG
#> 63          CC         NN         AC         AA         AA         AC
#> 64          AT         TT         AT         TT         AA         TT
#> 65          CT         CC         CT         CT         CT         NN
#> 66          AC         NN         AA         AC         AC         AC
#> 67          GG         NN         GG         GG         GG         GG
#> 68          TT         NN         AA         AT         AA         AT
#> 69          AA         AA         AA         AA         AA         AA
#> 70          AA         AC         AA         AA         NN         AC
#> 71          AT         AA         AA         AA         AA         AA
#> 72          AG         AG         AA         AG         AA         AA
#> 73          TT         TT         TT         TT         NN         TT
#> 74          GG         NN         GG         GT         GG         NN
#> 75          AC         AA         CC         AA         NN         CC
#> 76          TT         NN         NN         CC         CC         NN
#> 77          AT         AA         NN         AT         AA         AA
#> 78          GG         GG         GG         GG         GG         NN
#> 79          CT         TT         TT         CT         TT         NN
#> 80          CT         TT         TT         TT         CT         TT
#> 81          AT         TT         AA         AT         AA         AA
#> 82          TT         NN         NN         CT         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          TT         CC         TT         CC         CT         CT
#> 85          NN         GG         GG         AG         GG         GG
#> 86          AG         AA         AA         AG         AA         AA
#> 87          TT         CT         TT         TT         NN         CT
#> 88          NN         NN         AA         AG         AG         AA
#> 89          GG         NN         NN         CC         NN         CC
#> 90          CC         CC         CC         CC         CC         CC
#> 91          NN         CT         TT         CC         CT         TT
#> 92          AA         NN         NN         AG         GG         AG
#> 93          NN         CT         TT         TT         NN         TT
#> 94          CC         CC         CC         CC         CC         CC
#> 95          GG         NN         GT         GT         GG         GG
#> 96          CC         AC         CC         CC         CC         CC
#> 97          GG         NN         GG         GG         GG         GG
#> 98          TT         GT         GG         GG         GT         GG
#> 99          NN         AG         GG         AG         AG         AG
#> 100         AA         NN         AT         AT         AA         AT
#>     SMR.r1.023 SMR.r1.028 SMR.r1.034 SMR.r1.049 SMR.r1.052 SMR.r1.056
#> 1           AG         AA         GG         AG         AG         GG
#> 2           GG         GG         GG         NN         GG         GG
#> 3           CC         NN         CC         CC         CC         CC
#> 4           AG         AA         AG         GG         AG         AA
#> 5           NN         NN         AG         AG         GG         GG
#> 6           AA         NN         AA         AA         AA         AA
#> 7           GT         TT         TT         GT         NN         TT
#> 8           AG         AG         GG         AG         NN         AA
#> 9           AC         AC         AC         CC         CC         AC
#> 10          AA         AG         NN         AA         AG         AA
#> 11          CT         NN         CT         CC         CT         TT
#> 12          AA         NN         AG         AA         AA         NN
#> 13          TT         CT         CT         CT         CT         TT
#> 14          GG         GT         GG         GG         GT         GG
#> 15          AG         AA         AA         AG         AG         NN
#> 16          GT         GT         GG         TT         GT         NN
#> 17          AA         AA         AA         AA         AA         NN
#> 18          CC         CC         CC         CC         CC         CC
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         GG         GG         GG         GG         AG
#> 21          CC         CC         AC         AC         NN         CC
#> 22          NN         CT         CC         CC         CC         CC
#> 23          AG         GG         GG         NN         AG         GG
#> 24          AA         AA         AA         AA         AA         AG
#> 25          NN         GG         NN         NN         NN         AA
#> 26          CT         CT         TT         CT         CT         CT
#> 27          NN         AC         NN         AA         AA         AA
#> 28          CC         NN         CT         CT         TT         NN
#> 29          AA         NN         NN         AC         NN         NN
#> 30          GG         GG         NN         GG         GG         CG
#> 31          CG         CC         CC         GG         NN         CG
#> 32          AA         AC         NN         AC         AA         AA
#> 33          CT         CC         CT         CC         NN         CC
#> 34          NN         CC         CG         CC         CC         GG
#> 35          AA         AC         AA         AA         AA         AC
#> 36          AA         AA         AA         AG         AA         AA
#> 37          CG         CG         GG         CG         GG         NN
#> 38          NN         NN         GG         GG         GT         GT
#> 39          TT         NN         TT         TT         NN         NN
#> 40          CG         NN         GG         GG         GG         GG
#> 41          NN         NN         GG         AG         AG         GG
#> 42          GG         NN         AG         AG         AG         AG
#> 43          AA         AG         AA         AA         AA         AG
#> 44          NN         TT         TT         TT         TT         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          AC         CC         NN         AC         NN         CC
#> 47          CG         CG         CC         CC         CC         CC
#> 48          TT         CT         NN         CT         NN         TT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          NN         CC         CC         CC         CC         NN
#> 51          NN         CT         NN         NN         TT         CT
#> 52          AA         GG         AA         AG         AG         AA
#> 53          AA         AC         NN         AC         AA         CC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          GG         GG         GG         AG         GG         GG
#> 56          NN         AA         AG         AA         AA         NN
#> 57          GG         NN         GG         NN         GG         NN
#> 58          NN         CC         NN         CT         NN         CC
#> 59          AG         GG         GG         AG         AG         AG
#> 60          AA         AC         AA         AA         AC         AC
#> 61          TT         NN         NN         TT         CT         TT
#> 62          GG         GG         GG         GG         GG         GG
#> 63          AC         NN         AC         AC         NN         NN
#> 64          AT         TT         AT         TT         AT         AA
#> 65          CT         CC         TT         CC         CT         CC
#> 66          NN         NN         AA         NN         AC         AA
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AT         AA         AA         NN         AA         NN
#> 69          AG         AA         AA         AA         AA         NN
#> 70          AC         AA         AC         AA         AA         AA
#> 71          AA         AT         AA         AT         AA         AA
#> 72          NN         AA         AG         AA         NN         AG
#> 73          GT         TT         GT         TT         TT         TT
#> 74          GT         GT         NN         GT         NN         GG
#> 75          AC         AC         AC         AA         AC         NN
#> 76          CC         CC         NN         NN         NN         CC
#> 77          AA         AA         AT         AA         NN         AA
#> 78          GT         GG         GG         GG         GG         GG
#> 79          CT         TT         TT         TT         TT         TT
#> 80          CT         TT         CT         TT         CT         CT
#> 81          NN         AT         AT         AA         AT         AT
#> 82          TT         TT         TT         TT         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          NN         NN         CC         NN         NN         CC
#> 85          NN         GG         NN         GG         AG         AG
#> 86          AA         AA         AG         AA         NN         AA
#> 87          NN         CT         CT         CT         TT         TT
#> 88          AG         NN         NN         AG         AG         AG
#> 89          CG         NN         CC         CC         CG         NN
#> 90          CC         CC         CC         CC         CC         CC
#> 91          NN         CC         CT         CC         CC         CT
#> 92          NN         NN         NN         GG         GG         NN
#> 93          NN         NN         CT         TT         TT         TT
#> 94          CC         CC         CC         CC         NN         NN
#> 95          GT         NN         GG         GG         GG         NN
#> 96          AC         CC         CC         CC         CC         CC
#> 97          GG         GG         GG         GG         GG         NN
#> 98          NN         GT         GG         NN         GG         GT
#> 99          AA         NN         AG         AA         GG         AG
#> 100         NN         NN         AT         AA         AA         AT
#>     SMR.r1.072 SMR.r1.074 SMR.r2.001 SMR.r2.077 SMR.r2.086 SMR.r2.091
#> 1           AA         AG         AG         AG         AA         AA
#> 2           GG         GG         GG         NN         GG         GT
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AA         AG         AA         AA         AG         AA
#> 5           NN         NN         GG         NN         GG         GG
#> 6           AA         AC         AA         AA         AA         NN
#> 7           TT         TT         GT         TT         GT         TT
#> 8           AG         AG         AG         AG         NN         GG
#> 9           NN         AC         CC         CC         AC         NN
#> 10          AA         AA         AA         AA         AA         AA
#> 11          CT         CT         NN         NN         TT         CT
#> 12          AA         NN         NN         AA         AG         AG
#> 13          CT         TT         NN         TT         CC         NN
#> 14          GT         GG         GG         GT         GT         NN
#> 15          AA         AA         AG         AA         AA         AA
#> 16          GT         GG         GG         GT         GT         NN
#> 17          AC         AA         AA         AA         AC         AA
#> 18          AC         CC         CC         AC         AA         AC
#> 19          GT         TT         GT         NN         TT         TT
#> 20          GG         GG         NN         GG         GG         GG
#> 21          AC         AC         NN         AC         CC         AC
#> 22          NN         CT         CC         CC         NN         CC
#> 23          NN         GG         NN         AG         GG         GG
#> 24          AA         AA         NN         NN         AA         AA
#> 25          AG         GG         AG         NN         AG         AG
#> 26          CT         CT         CT         CT         TT         CT
#> 27          AC         AA         AC         AA         AA         AA
#> 28          NN         CT         CT         CT         NN         CT
#> 29          AA         AA         NN         NN         AA         AC
#> 30          GG         GG         GG         CG         GG         GG
#> 31          CC         CG         CG         CG         CG         NN
#> 32          AA         AA         AC         AA         NN         AC
#> 33          CC         TT         NN         CC         CC         TT
#> 34          CC         CC         NN         CC         CG         CG
#> 35          NN         AC         AA         AA         AA         AA
#> 36          AG         AA         AA         AA         AG         AA
#> 37          CG         GG         GG         CG         GG         CC
#> 38          GT         GT         GG         GG         GT         GG
#> 39          TT         NN         TT         TT         TT         AT
#> 40          CG         GG         NN         GG         NN         GG
#> 41          AG         AG         NN         NN         GG         GG
#> 42          AA         AG         AA         AG         AA         AG
#> 43          AA         AA         AA         AG         AA         AA
#> 44          TT         TT         TT         TT         TT         TT
#> 45          GG         GG         NN         NN         GG         NN
#> 46          AC         NN         NN         CC         CC         CC
#> 47          CG         CC         CC         NN         CG         CG
#> 48          CT         CT         NN         NN         CC         CT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          NN         NN         CC         CC         NN         CC
#> 51          NN         CT         CT         TT         CT         TT
#> 52          AA         AA         AA         AA         AA         AA
#> 53          AA         AA         AA         AA         AA         AA
#> 54          GG         GG         NN         GG         GT         GG
#> 55          GG         AG         NN         NN         GG         GG
#> 56          NN         AG         AA         AA         AG         AG
#> 57          GG         GG         AG         GG         GG         GG
#> 58          CC         CT         NN         CC         CC         CC
#> 59          AA         AG         AG         GG         AG         AG
#> 60          CC         AA         CC         AC         AA         AC
#> 61          TT         TT         NN         TT         TT         TT
#> 62          GG         GG         GG         GG         GG         GG
#> 63          AA         AA         NN         NN         AC         NN
#> 64          AT         AA         TT         NN         AT         NN
#> 65          CC         NN         CT         CC         CT         CT
#> 66          NN         AA         CC         AC         CC         AA
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AT         AT         AA         AA         AT         AT
#> 69          AA         AA         AA         AA         AG         AA
#> 70          AA         AA         AC         AA         AC         AC
#> 71          AA         AA         AA         AA         AT         AA
#> 72          AA         AG         AG         AA         AG         AA
#> 73          TT         GT         GT         TT         TT         TT
#> 74          GG         NN         GG         NN         GG         GT
#> 75          AC         AC         CC         AA         AC         AC
#> 76          CT         NN         NN         NN         CC         CC
#> 77          AA         AT         NN         NN         AT         AT
#> 78          GG         GG         GG         GG         GT         GG
#> 79          TT         NN         TT         TT         TT         TT
#> 80          CT         CT         CT         CT         CT         CT
#> 81          AA         NN         AA         AT         AT         AA
#> 82          TT         NN         NN         CT         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          CT         CT         CC         CC         CT         CT
#> 85          AG         GG         AG         GG         AG         AG
#> 86          AA         AG         NN         NN         AG         AG
#> 87          TT         TT         CT         TT         CT         TT
#> 88          AG         NN         NN         NN         GG         AA
#> 89          CG         CG         CC         CG         CG         GG
#> 90          CC         CC         CC         CC         CC         CC
#> 91          TT         CT         NN         NN         CT         CT
#> 92          NN         NN         NN         AG         NN         AG
#> 93          CT         TT         TT         CT         TT         TT
#> 94          CC         CC         CC         CC         CC         CC
#> 95          GT         GT         NN         GT         GG         GG
#> 96          CC         CC         CC         CC         AC         CC
#> 97          GG         NN         GG         GG         GG         AG
#> 98          GG         GG         GG         GT         GT         GG
#> 99          AA         NN         AA         AG         GG         AG
#> 100         AT         NN         AT         AT         AT         AA
#>     SMR.r2.097 SMR.r2.105 SMR.r2.139 SMR.r2.175 SMR.r2.256 SMR.r2.263
#> 1           NN         AA         AG         AG         AG         AA
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CT         CC         CC         CC         CT         CC
#> 4           AA         GG         AA         AG         AG         GG
#> 5           GG         GG         NN         AA         AG         GG
#> 6           AA         AA         AA         AA         AA         AA
#> 7           TT         TT         GT         TT         TT         GT
#> 8           GG         AA         NN         GG         GG         AG
#> 9           AC         NN         CC         AC         CC         CC
#> 10          AA         AA         AA         AA         AA         AA
#> 11          CT         CT         NN         CC         CC         CT
#> 12          AG         AA         NN         NN         AG         AA
#> 13          CT         TT         TT         NN         NN         TT
#> 14          GG         GT         GG         GT         GG         GT
#> 15          AA         AG         AG         NN         AA         NN
#> 16          GG         NN         GT         GT         GG         GT
#> 17          AA         AA         AA         AC         NN         AC
#> 18          AC         CC         CC         CC         CC         CC
#> 19          TT         TT         GT         TT         TT         TT
#> 20          GG         GG         GG         GG         GG         GG
#> 21          AC         AC         NN         NN         CC         CC
#> 22          CT         CC         CT         CC         CT         CC
#> 23          GG         GG         GG         GG         NN         GG
#> 24          AA         AA         NN         NN         AA         AA
#> 25          AA         AG         GG         GG         AG         AG
#> 26          TT         TT         CT         TT         TT         NN
#> 27          NN         AA         AA         AA         AC         AC
#> 28          CT         CT         TT         CT         CT         CC
#> 29          CC         AA         NN         AC         AA         AA
#> 30          GG         GG         NN         CG         GG         GG
#> 31          CC         CG         CG         NN         CC         CC
#> 32          AA         AA         AA         AA         NN         NN
#> 33          CT         TT         CC         TT         CC         CT
#> 34          CC         CG         GG         CC         CC         CG
#> 35          AA         AA         NN         AA         AA         AA
#> 36          NN         NN         AG         AA         AA         AG
#> 37          CG         CG         NN         GG         GG         CG
#> 38          NN         GT         GT         GT         GT         TT
#> 39          TT         TT         TT         NN         AT         TT
#> 40          CG         GG         NN         NN         CG         CG
#> 41          AG         GG         NN         GG         NN         GG
#> 42          AG         AG         NN         AG         AA         AG
#> 43          AA         AA         AA         AG         AA         AA
#> 44          AT         AT         TT         AT         TT         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          NN         CC         CC         NN         CC         CC
#> 47          GG         CC         CC         CC         CC         CG
#> 48          CT         TT         CT         CT         NN         CC
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CC         CT         CC         CC
#> 51          NN         CC         TT         TT         NN         NN
#> 52          AA         AG         AA         AA         AA         AA
#> 53          AC         AC         AA         AA         AC         AC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          NN         AG         GG         NN         GG         GG
#> 56          AG         AA         AA         AA         AG         NN
#> 57          GG         GG         NN         GG         GG         GG
#> 58          CT         CT         CC         CC         CT         CC
#> 59          GG         AG         AG         NN         AA         AG
#> 60          AC         AA         AC         AC         NN         AA
#> 61          TT         NN         TT         TT         CT         TT
#> 62          GG         GG         CG         GG         GG         GG
#> 63          CC         AA         AC         AA         AA         AC
#> 64          TT         AT         AT         AT         AT         TT
#> 65          CC         TT         CT         CC         CC         TT
#> 66          AA         AC         NN         NN         AA         CC
#> 67          GG         GG         GG         GG         GG         NN
#> 68          TT         AA         AT         NN         NN         NN
#> 69          AA         NN         AG         AA         AA         AG
#> 70          AA         AC         AC         NN         AA         AC
#> 71          AT         AT         AA         AA         AA         AT
#> 72          AA         AG         AA         AG         AA         AA
#> 73          NN         GT         TT         TT         TT         TT
#> 74          GT         GG         NN         GT         GT         GG
#> 75          AA         CC         AC         NN         NN         AC
#> 76          CC         NN         CC         CT         NN         CT
#> 77          AT         AA         AA         NN         AT         AA
#> 78          GG         NN         GT         GG         GG         GT
#> 79          TT         TT         TT         CC         TT         TT
#> 80          NN         TT         TT         TT         TT         CT
#> 81          AT         AA         AT         AA         AT         AT
#> 82          TT         TT         TT         NN         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          CC         CC         CC         NN         CC         CT
#> 85          GG         GG         AG         AG         GG         AA
#> 86          AG         AA         NN         AA         AG         AA
#> 87          CT         CT         CT         CT         TT         TT
#> 88          AA         NN         AG         AG         AG         AG
#> 89          CC         CC         CG         CG         NN         GG
#> 90          CC         CT         CC         CC         CC         CC
#> 91          TT         TT         TT         CT         NN         CT
#> 92          GG         GG         NN         NN         GG         NN
#> 93          TT         CT         NN         NN         CT         TT
#> 94          CC         CC         CC         CG         CC         CC
#> 95          GG         GG         NN         GG         GG         GG
#> 96          NN         CC         CC         CC         CC         AC
#> 97          NN         GG         GG         GG         GG         GG
#> 98          GG         GG         GT         GG         GG         GT
#> 99          AA         AG         GG         NN         NN         GG
#> 100         AT         AT         AA         AT         AT         AA
#>     SMR.r2.296 UPD.r2.120 UPD.r2.121 UPD.r2.122 UPD.r2.123 UPD.r2.124
#> 1           GG         AA         AA         AG         AG         AG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CT         CT         CC         CC         CC
#> 4           AA         AG         GG         AG         AG         AG
#> 5           AG         GG         NN         AA         AG         AG
#> 6           AA         AC         AC         AC         AC         CC
#> 7           GT         TT         TT         TT         TT         TT
#> 8           GG         NN         GG         AG         AG         NN
#> 9           CC         CC         CC         AC         CC         NN
#> 10          AA         AA         AA         AA         AA         AG
#> 11          TT         CT         CT         TT         TT         CT
#> 12          AA         AG         AA         AG         AG         AA
#> 13          TT         TT         CT         TT         TT         TT
#> 14          GG         GG         GT         GT         GT         GT
#> 15          AG         AG         AG         GG         NN         AG
#> 16          GT         GG         GG         GT         GT         GG
#> 17          AA         AC         AA         AA         AC         AA
#> 18          CC         AA         CC         AA         CC         CC
#> 19          TT         GT         TT         TT         TT         TT
#> 20          GG         GG         GG         NN         AG         GG
#> 21          AC         CC         AC         AC         NN         CC
#> 22          CC         CC         TT         CT         CT         CT
#> 23          GG         GG         GG         GG         GG         AG
#> 24          AA         AA         AA         AA         AA         AA
#> 25          AG         GG         AG         AG         AA         GG
#> 26          CT         CT         CT         CT         TT         CT
#> 27          AA         AA         AA         AC         AA         AA
#> 28          CC         TT         CT         CT         TT         TT
#> 29          AA         AA         AA         AA         AC         AA
#> 30          GG         GG         GG         GG         NN         GG
#> 31          CG         GG         CG         CG         CC         CC
#> 32          AC         AC         AA         AA         AC         AC
#> 33          CT         CC         TT         CC         TT         CC
#> 34          GG         CG         CC         CG         CC         CC
#> 35          AA         AC         AA         AA         AA         AA
#> 36          AG         AA         AA         AA         AG         AA
#> 37          CC         GG         GG         GG         GG         CG
#> 38          NN         GG         GG         GT         GG         GT
#> 39          TT         TT         AT         TT         TT         TT
#> 40          CG         CG         CG         GG         CG         GG
#> 41          AG         AA         NN         AA         NN         GG
#> 42          AG         AG         AG         GG         GG         AG
#> 43          AG         NN         AA         AG         AG         NN
#> 44          TT         TT         AT         TT         NN         TT
#> 45          NN         GG         GG         GG         GG         GG
#> 46          AA         CC         AC         CC         NN         CC
#> 47          GG         CC         CC         CG         CC         CC
#> 48          NN         CT         NN         TT         CT         CT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CC         CT         CC         CC
#> 51          TT         CC         TT         CC         CT         TT
#> 52          AA         GG         GG         AA         AG         AG
#> 53          AA         AA         CC         AA         AC         NN
#> 54          GG         GG         GG         GG         GG         GG
#> 55          AG         AA         AG         AG         AG         GG
#> 56          AA         AG         AA         AG         AG         AA
#> 57          NN         GG         GG         GG         GG         GG
#> 58          NN         NN         TT         CT         NN         NN
#> 59          AA         GG         AG         GG         AG         GG
#> 60          AA         AC         AA         NN         AC         AC
#> 61          TT         TT         CT         TT         TT         TT
#> 62          GG         GG         NN         GG         GG         GG
#> 63          NN         AC         AC         AC         AA         CC
#> 64          TT         AA         TT         TT         AT         AT
#> 65          CC         TT         CC         TT         CC         CT
#> 66          AC         AC         AC         AA         AA         AC
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AA         AA         TT         AT         AA         AT
#> 69          NN         AA         AA         AA         AA         AA
#> 70          AC         AA         AC         NN         AC         AA
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AG         AG         GG         AG         AG         NN
#> 73          NN         GT         TT         TT         GT         GT
#> 74          GG         GT         GG         GG         GG         GG
#> 75          AC         AC         AC         AC         AC         AC
#> 76          NN         CT         CC         CC         CT         CT
#> 77          AA         AT         AA         AT         AT         AA
#> 78          GT         GG         GG         GG         GG         GG
#> 79          TT         CT         TT         TT         TT         NN
#> 80          TT         TT         CT         NN         CT         CT
#> 81          AT         AT         AA         AT         NN         AT
#> 82          TT         TT         TT         TT         CT         TT
#> 83          AT         TT         TT         TT         TT         TT
#> 84          CT         CT         TT         CC         CT         CT
#> 85          NN         GG         GG         AG         GG         AG
#> 86          NN         AG         AA         AG         AG         NN
#> 87          CT         CT         TT         NN         TT         TT
#> 88          NN         NN         AG         AG         AG         AG
#> 89          CG         GG         CC         CC         CG         NN
#> 90          NN         CC         CC         CC         CT         NN
#> 91          CT         CT         CT         CC         TT         TT
#> 92          GG         NN         GG         GG         GG         AG
#> 93          CT         TT         TT         TT         CT         CT
#> 94          CC         CC         CC         CC         CC         CG
#> 95          GG         GG         GT         GG         GG         GG
#> 96          CC         CC         CC         CC         NN         CC
#> 97          AG         GG         GG         GG         GG         GG
#> 98          GT         GT         GG         GG         GT         GT
#> 99          AA         AG         GG         AA         NN         NN
#> 100         AA         AA         AT         AA         AT         AT
#>     UPD.r2.129 UPD.r2.130 UPD.r2.131 UPD.r2.133 UPD.r2.135 UPD.r2.160
#> 1           AA         AG         AG         GG         AG         NN
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AG         AG         NN         NN         AG         AG
#> 5           NN         NN         NN         AA         AG         AG
#> 6           AC         NN         AA         AC         AC         AA
#> 7           TT         TT         TT         TT         TT         TT
#> 8           GG         NN         NN         GG         AG         AG
#> 9           CC         AC         CC         AC         AC         CC
#> 10          NN         AG         AA         AA         AA         AA
#> 11          NN         CT         CT         CT         CT         CT
#> 12          NN         AA         AG         NN         NN         NN
#> 13          TT         TT         TT         CT         CT         CT
#> 14          GT         GG         GG         GT         GG         TT
#> 15          AG         AG         NN         AG         AG         AG
#> 16          GT         NN         GG         NN         GG         GG
#> 17          NN         AA         AA         AA         AC         NN
#> 18          CC         AA         CC         AC         AC         CC
#> 19          TT         TT         TT         GT         TT         TT
#> 20          GG         GG         GG         GG         GG         GG
#> 21          NN         AC         CC         AC         CC         AC
#> 22          CT         CC         NN         CC         CT         CT
#> 23          GG         GG         GG         GG         GG         GG
#> 24          AA         AA         AG         AG         AA         NN
#> 25          AG         AG         GG         AA         GG         AG
#> 26          TT         CT         CT         TT         CT         NN
#> 27          AC         AA         AA         AA         AA         AA
#> 28          TT         CT         CT         NN         TT         TT
#> 29          NN         AA         AA         NN         AA         AA
#> 30          GG         CG         CG         NN         GG         NN
#> 31          CC         CG         CG         GG         CC         GG
#> 32          AA         AA         AC         AC         AC         AA
#> 33          CT         CC         NN         CT         CC         CC
#> 34          NN         CC         CC         CC         CC         NN
#> 35          AA         AA         AA         AA         AA         AA
#> 36          AA         AA         AG         AA         NN         AA
#> 37          CG         GG         CG         CG         CG         CG
#> 38          GG         GG         TT         TT         GG         GT
#> 39          NN         NN         AT         NN         TT         AT
#> 40          NN         CG         CG         NN         CC         GG
#> 41          NN         AG         AA         AA         AG         NN
#> 42          GG         AG         AA         GG         NN         AG
#> 43          AG         AA         AG         AA         AG         AG
#> 44          AT         AT         NN         TT         NN         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          CC         AC         AA         AC         CC         AC
#> 47          CC         CC         CG         CC         CC         CC
#> 48          CT         TT         CT         CT         CT         TT
#> 49          GG         NN         CG         GG         GG         NN
#> 50          CC         CT         CC         CC         CC         NN
#> 51          CT         CT         CT         TT         CT         NN
#> 52          NN         AG         AG         AA         AA         AG
#> 53          NN         AA         AC         AC         AC         AC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          GG         AG         AA         GG         NN         NN
#> 56          AG         AA         AG         AA         AA         AA
#> 57          NN         GG         NN         GG         NN         GG
#> 58          CT         CT         CC         CC         CC         NN
#> 59          GG         GG         AG         AG         AG         AA
#> 60          AC         AC         AA         AA         AA         AC
#> 61          NN         TT         TT         TT         NN         TT
#> 62          GG         CG         CG         CG         GG         GG
#> 63          AA         AC         AC         AA         NN         AA
#> 64          TT         TT         AT         AT         AT         AT
#> 65          CC         CT         NN         CC         CC         CC
#> 66          NN         AA         AA         AA         CC         CC
#> 67          GG         GG         GG         NN         NN         AG
#> 68          AT         AA         AA         AT         AT         NN
#> 69          AA         AA         AA         AA         AA         AA
#> 70          AA         AC         AC         AA         AA         AC
#> 71          AA         AA         AA         AT         AA         AA
#> 72          AG         GG         NN         GG         AG         GG
#> 73          TT         GG         TT         TT         GT         TT
#> 74          GG         NN         GG         GG         GG         GT
#> 75          AC         AC         CC         CC         AC         CC
#> 76          CT         CT         CC         CC         TT         CT
#> 77          NN         AA         AT         AA         AA         AA
#> 78          GG         GG         GG         GG         GG         GG
#> 79          TT         TT         TT         NN         NN         TT
#> 80          CT         CT         TT         CT         TT         NN
#> 81          AT         AT         AA         AA         AT         AT
#> 82          TT         NN         CT         TT         CT         CT
#> 83          TT         TT         AT         TT         TT         TT
#> 84          CT         CT         CT         TT         CT         NN
#> 85          AA         AG         AG         AA         AG         AG
#> 86          AG         AA         AG         AA         NN         AA
#> 87          CT         CT         CC         TT         CT         CT
#> 88          AA         AA         AG         AA         NN         AG
#> 89          GG         CG         NN         CG         NN         CG
#> 90          CC         CC         CT         CC         CC         CC
#> 91          CT         CC         CC         CC         CC         CT
#> 92          NN         NN         NN         AG         GG         AG
#> 93          TT         TT         CT         TT         TT         CT
#> 94          CC         CG         CC         CC         NN         CG
#> 95          GT         GT         GT         GG         GG         NN
#> 96          CC         CC         CC         NN         AC         CC
#> 97          GG         GG         AG         GG         GG         GG
#> 98          NN         GG         GG         GT         NN         GT
#> 99          GG         AG         AG         GG         NN         NN
#> 100         AT         TT         TT         TT         NN         TT
#>     UPD.r2.162 UPD.r2.174 UPD.r2.189 UPD.r2.201
#> 1           GG         NN         GG         AG
#> 2           GG         GG         GG         GG
#> 3           CC         CC         CC         CC
#> 4           AG         AA         AA         AA
#> 5           NN         GG         NN         GG
#> 6           AC         AA         AA         AA
#> 7           NN         TT         TT         TT
#> 8           NN         AG         GG         GG
#> 9           NN         CC         AC         AA
#> 10          AA         AG         AA         AG
#> 11          CT         CT         CT         TT
#> 12          AG         AA         AG         NN
#> 13          TT         TT         NN         NN
#> 14          GT         GG         GT         GT
#> 15          AG         AG         AA         AA
#> 16          GG         TT         NN         GG
#> 17          AA         AC         AA         AA
#> 18          AC         AA         AC         AC
#> 19          TT         TT         TT         TT
#> 20          GG         GG         AG         GG
#> 21          NN         AC         AA         CC
#> 22          CT         CT         CC         CC
#> 23          GG         GG         GG         GG
#> 24          AG         AA         AA         AA
#> 25          AG         AA         AG         GG
#> 26          NN         NN         TT         TT
#> 27          AA         NN         AA         AA
#> 28          CT         CT         CT         CT
#> 29          AA         AA         AC         AC
#> 30          GG         CG         GG         CG
#> 31          CG         CG         CC         CG
#> 32          NN         AA         NN         AC
#> 33          CC         CT         NN         TT
#> 34          CG         NN         CC         NN
#> 35          NN         AA         AA         AC
#> 36          NN         AA         AA         AA
#> 37          NN         CG         NN         GG
#> 38          GG         GG         GG         NN
#> 39          TT         NN         TT         NN
#> 40          GG         NN         CC         NN
#> 41          NN         GG         AG         AG
#> 42          AG         GG         AG         AG
#> 43          AA         AA         AA         AA
#> 44          TT         AT         AT         TT
#> 45          GG         GT         GG         GT
#> 46          CC         AA         NN         AC
#> 47          GG         CC         CC         CG
#> 48          CT         NN         CC         TT
#> 49          GG         CG         GG         GG
#> 50          CC         CC         CC         CC
#> 51          CT         CT         TT         TT
#> 52          AG         AA         AA         AA
#> 53          NN         AA         AC         CC
#> 54          GG         GG         GG         GG
#> 55          AA         AG         AG         GG
#> 56          AG         AA         AG         AA
#> 57          GG         NN         NN         GG
#> 58          CC         CC         CT         CC
#> 59          AG         AA         NN         NN
#> 60          AC         AC         AA         AC
#> 61          TT         CT         TT         NN
#> 62          GG         GG         GG         GG
#> 63          AA         AA         NN         AA
#> 64          TT         AT         AT         AT
#> 65          CT         CT         CT         TT
#> 66          AC         CC         NN         AC
#> 67          GG         GG         GG         GG
#> 68          TT         AA         AA         AT
#> 69          AA         AA         AA         AA
#> 70          NN         NN         AC         AC
#> 71          AA         AT         AA         AA
#> 72          GG         GG         AG         AG
#> 73          TT         TT         NN         GT
#> 74          GG         GT         NN         GG
#> 75          CC         CC         AC         CC
#> 76          TT         CC         NN         TT
#> 77          AT         AA         AT         AA
#> 78          GG         GG         GG         GG
#> 79          CT         CT         TT         TT
#> 80          CT         TT         TT         TT
#> 81          TT         AA         AT         AT
#> 82          TT         TT         NN         CT
#> 83          TT         TT         NN         TT
#> 84          CC         CT         NN         CC
#> 85          AG         GG         AG         AG
#> 86          AG         AA         AG         AA
#> 87          CT         CT         CC         CT
#> 88          NN         AG         AA         AG
#> 89          GG         GG         GG         CC
#> 90          CC         CC         NN         CC
#> 91          NN         CC         CT         CT
#> 92          GG         NN         AG         GG
#> 93          NN         CT         CT         CT
#> 94          CC         CC         CC         CC
#> 95          GG         NN         NN         GG
#> 96          CC         CC         AC         CC
#> 97          GG         GG         GG         GG
#> 98          GG         GT         GT         TT
#> 99          AA         AA         AG         AG
#> 100         TT         AT         AT         AT

# from an adegenet genind object
ex.genind  <- adegenet::df2genind(t(stickRAW[,-c(1:2)]), 
                                  ncode = 1, NA.char = "N") # get genind data
# note, will add whatever metadata data is in the genind object to the 
# snpRdata object. 
# Could be run without the snp or sample metadatas.
import.snpR.data(ex.genind, snp_meta, sample_meta) 
#> Warning: Some sample metadata columns contain unacceptable special characters. Unaccepted characters: '.', '~', or any whitespace.
#> These can cause unexpected behaviour if the subect columns are used as facets.
#> Issues:
#> 
#> Facet: sampID    levels: ASP.r2.141, ASP.r2.144, ASP.r2.146, ASP.r2.147, ASP.r2.196, ASP.r2.216, ASP.r2.221, ASP.r2.227, ASP.r2.230, ASP.r2.232, ASP.r2.235, CLF.r2.040, CLF.r2.048, CLF.r2.061, CLF.r2.064, CLF.r2.068, CLF.r2.069, OPL.r1.078, OPL.r1.080, OPL.r1.083, OPL.r1.086, OPL.r1.087, OPL.r1.088, OPL.r1.109, OPL.r1.110, OPL.r1.116, OPL.r1.125, OPL.r1.126, OPL.r1.128, OPL.r1.130, OPL.r1.139, OPL.r1.146, OPL.r1.147, OPL.r1.149, OPL.r2.007, OPL.r2.012, OPL.r2.022, OPL.r2.023, OPL.r2.028, OPL.r2.031, OPL.r2.033, OPL.r2.079, OPL.r2.080, OPL.r2.094, PAL.r2.203, PAL.r2.237, PAL.r2.238, PAL.r2.239, PAL.r2.240, PAL.r2.241, PAL.r2.242, PAL.r2.249, PAL.r2.253, PAL.r2.275, PAL.r2.276, PAL.r2.278, PAL.r2.281, PAL.r2.285, PAL.r2.287, PAL.r2.290, PAL.r2.291, SMR.r1.003, SMR.r1.004, SMR.r1.008, SMR.r1.012, SMR.r1.013, SMR.r1.023, SMR.r1.028, SMR.r1.034, SMR.r1.049, SMR.r1.052, SMR.r1.056, SMR.r1.072, SMR.r1.074, SMR.r2.001, SMR.r2.077, SMR.r2.086, SMR.r2.091, SMR.r2.097, SMR.r2.105, SMR.r2.139, SMR.r2.175, SMR.r2.256, SMR.r2.263, SMR.r2.296, UPD.r2.120, UPD.r2.121, UPD.r2.122, UPD.r2.123, UPD.r2.124, UPD.r2.129, UPD.r2.130, UPD.r2.131, UPD.r2.133, UPD.r2.135, UPD.r2.160, UPD.r2.162, UPD.r2.174, UPD.r2.189, UPD.r2.201
#>     V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 V20 V21
#> 1   AG GG AA AA NN AG NN GG AA  GG  AG  NN  AG  AG  AA  AG  AA  AG  AG  AG  AG
#> 2   GG GG GG GG GT GG GG NN GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  GT  GG
#> 3   CC CC CC CC CC CC CT CC CC  CT  CC  CC  CC  CC  CC  CC  CT  CC  CT  CT  CC
#> 4   AA AA AA AA AG AA AA AA AA  AG  AG  AA  AG  AA  NN  AA  AA  AG  AG  NN  AG
#> 5   GG NN NN AA GG GA GA GA NN  GG  NN  GG  NN  NN  NN  GG  GG  AA  NN  NN  GA
#> 6   AC AC AA AA AA AC AA AC AA  AC  NN  AA  AC  AC  AA  AC  AA  AA  AA  AC  AA
#> 7   TT TT TT TT TG TT TT TT TT  TT  TT  GG  TT  TT  TT  TT  TT  TT  NN  TG  TT
#> 8   GG GG GA GA GA GG GG GA GG  GA  GG  GG  GA  GA  AA  GA  GG  GG  GG  GG  GG
#> 9   CC CA CC CA CC CC CC CC CA  CC  AA  CC  CA  CC  CA  CC  CC  CC  CA  CC  CC
#> 10  AA AA AG AA AA AA AA AA AG  AA  AA  AA  AA  AA  AA  AG  AA  AA  AA  AA  AA
#> 11  NN TT TT TC TT TC CC CC TC  TC  TC  TC  TC  TC  TC  TT  CC  CC  NN  NN  CC
#> 12  AG AG AG NN AA AA AG AA AG  AA  AA  AA  AA  AG  AA  AA  AA  AA  NN  NN  AA
#> 13  TT TT TC TT TC NN TT TT NN  TT  TT  TC  TC  TT  TC  TT  CC  NN  TT  NN  TT
#> 14  GG GT GT NN GG GG GG GG GT  GG  GT  GG  GG  GG  GT  GT  GG  GG  GG  GG  GG
#> 15  AA AG AG AA AA AA AA AA AG  AG  AG  AA  AA  AA  AA  AA  AA  NN  AG  AA  AA
#> 16  NN GG GT TT GG GT GT GG GG  GG  GG  GT  GT  GG  NN  GG  GT  NN  GT  NN  NN
#> 17  AA AA AA AA AA NN AC AC AA  AC  AA  AA  AA  AA  AA  AA  AA  AC  AA  AC  NN
#> 18  AC NN CC CC AC AC CC AC AC  CC  AC  CC  AC  AC  NN  AC  AC  CC  CC  CC  CC
#> 19  TT TT TT TT TT TT TT TT TT  TT  TT  TG  TT  TT  TG  TT  TG  NN  TT  TT  TG
#> 20  GG GA GG GG NN GG GG GG GG  NN  GG  GG  GG  GG  GG  GG  GA  GG  GG  NN  GG
#> 21  NN NN CC CC CA CA CC CC CA  CC  CC  NN  CA  CA  NN  CA  CC  NN  CA  AA  NN
#> 22  CT NN TT CC CT NN CC TT CT  CC  CT  CC  CC  NN  TT  CT  CT  CC  CT  CT  NN
#> 23  GG GG GG GG GG GG GG GA GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GA  GG
#> 24  AA AA AA NN AA AA AA AA AA  AA  AG  AA  AA  AA  AA  AA  AA  AA  AA  AA  AG
#> 25  NN AG AG GG GG AG AG GG AG  AG  GG  GG  GG  AG  NN  GG  GG  AG  NN  AG  GG
#> 26  CT CT CT CT NN TT NN TT CT  TT  CT  TT  TT  TT  TT  TT  TT  CT  TT  TT  NN
#> 27  AA AA AA AA AA AA AA AA NN  AC  AA  AA  AC  AA  AA  AA  AA  AC  AC  AA  AC
#> 28  TT NN NN NN TC TC TT TT TT  TT  TT  TC  TT  TC  TT  TT  TT  NN  TC  TC  NN
#> 29  AA AC AA NN AA AA AA NN AA  AA  AA  AA  NN  NN  AA  AC  AA  AA  AA  NN  AA
#> 30  GG GG GC GG GG GG GG GG GG  GG  GC  GG  GG  GC  GG  GC  NN  GG  NN  GG  NN
#> 31  CG CG CC CC CC CC CG CC CG  CG  CG  CC  CG  CC  CC  CC  CC  CG  CC  CG  CC
#> 32  AA AA AA NN AA AA AA AC AC  AC  AA  AC  AC  NN  AA  AA  AA  AA  AA  AC  NN
#> 33  CC CT CC CC CT CT CT CC CC  CT  CT  CT  CC  CT  NN  CT  NN  NN  NN  CC  CT
#> 34  CG CG CG CG CC NN CC CG CG  CG  NN  CC  CG  GG  CG  GG  CC  NN  NN  CG  NN
#> 35  AA AA AA AC AA NN AC AA AA  AA  AA  AA  NN  AA  AA  AA  AA  AA  AA  AA  NN
#> 36  AA AA AA AA AA AG AA AG AA  AG  NN  AG  AA  AA  AG  AA  AA  NN  AG  AA  AA
#> 37  CG GG GG GG GG GG NN GG CC  GG  CG  GG  CG  CG  NN  GG  GG  GG  CG  GG  CG
#> 38  GG TT GG GT GT TT GG GT NN  GG  TT  GG  GG  GG  GG  GG  GG  GT  GG  NN  NN
#> 39  TT NN TT NN TA TT TT TT TT  TA  TA  TT  NN  NN  TA  NN  TA  NN  NN  NN  TT
#> 40  CG CC CG NN CG GG CG GG CC  GG  CG  CG  GG  NN  CG  GG  GG  CG  NN  GG  GG
#> 41  AG AA AG AG AA AA AG AG GG  AG  AG  AG  GG  GG  AG  NN  AG  AG  NN  NN  NN
#> 42  AG NN NN NN GG GG AG AG AG  GG  GG  NN  GG  AG  GG  AG  AA  GG  AG  AG  AG
#> 43  AG AG AA AA AA AA AA AA AA  AA  AA  AG  NN  AA  AA  AG  AA  NN  AA  AA  NN
#> 44  AT TT TT TT NN TT TT NN TT  AA  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  NN
#> 45  GG GG GG GG GG GG GG GG GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  GG  NN  GG
#> 46  NN CC AA AA NN CC NN AA CC  CC  CA  CA  CC  CC  NN  CC  CA  CA  CA  NN  NN
#> 47  CC CC CC CG CG CC CG CG CC  CG  CC  CC  CC  CC  CC  CC  CC  NN  CG  CC  CG
#> 48  CT TT TT NN TT TT CT CC NN  TT  CT  CT  CT  TT  CT  CT  CC  CT  CC  NN  CT
#> 49  GG GG GG GG NN GG GG GG GG  GG  GG  GG  GG  GC  GG  GC  GG  GG  GG  GG  GG
#> 50  NN CC CT CC CC CT CC CC CT  CC  CC  CC  CC  CC  CC  CT  NN  CT  CC  CC  CT
#> 51  TT NN TC TC TT TC TC NN TC  CC  TT  TT  NN  TT  TT  TT  TT  TC  CC  NN  TT
#> 52  AA AG AG AA AG AG AA AG AG  AG  AG  GG  AG  AA  AA  AA  AG  AG  AG  AA  NN
#> 53  AA AA CC NN AA AC AA AA AC  AA  AA  AC  NN  AA  NN  AC  AC  NN  AA  AA  AC
#> 54  GG GG GG GG GG GG GG GG GG  GG  GG  GG  GG  GG  GT  GG  GG  NN  GG  GG  GG
#> 55  GG NN GA NN GA GA GA GA GA  AA  NN  GG  GG  GG  GA  GA  GA  AA  NN  NN  NN
#> 56  AG AG AG AG AA NN AG AA AG  AA  NN  AA  AA  AG  AA  AA  AA  NN  AA  NN  NN
#> 57  GG NN GG NN GG GG GG GA NN  GG  GG  NN  GG  GG  GG  GG  GG  GA  NN  NN  NN
#> 58  TT TC CC NN TT CC NN TC CC  TC  TC  CC  CC  CC  CC  CC  CC  TC  TC  CC  CC
#> 59  AA AG GG AG GG AA GG AG GG  AA  AG  AG  AG  AG  AG  GG  GG  AG  AG  AG  AG
#> 60  AA AA AC NN AC AC AC NN NN  CC  NN  AC  AC  CC  CC  CC  AC  AC  AC  AC  NN
#> 61  TT NN TT NN TT TT TT NN TT  TT  TT  TT  NN  NN  TC  TT  NN  TT  NN  NN  TT
#> 62  CG GG GG NN GG GG GG GG CG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  CG  GG
#> 63  AC NN NN NN AA AC AC AC CC  AA  AA  CC  AC  AC  AA  NN  AA  NN  AC  NN  AA
#> 64  TT TA TA TT TA TA AA TA TT  TT  TT  TT  TT  TT  TT  TT  AA  TA  NN  TA  NN
#> 65  CT CC CT CC CC CT CT CC CC  CC  CC  TT  CC  CT  CT  CT  CC  NN  NN  TT  NN
#> 66  AC AC NN NN AC CC CC AC AA  AA  AC  AC  AC  CC  AA  AC  AA  NN  AC  AC  AC
#> 67  GG GG GG NN GG GG GG GG NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  GG  GG
#> 68  AA AT AT AA AA NN AT TT AA  AA  AA  AT  AA  AA  AT  AA  AT  AT  AT  AA  AT
#> 69  AA AA NN AA AA AA AA AA AA  AA  AA  NN  AA  AA  AA  AG  AA  GG  AA  AA  NN
#> 70  AC AC AA AC AC NN AC AC AA  AA  AC  CC  AA  AC  NN  AC  CC  NN  AA  NN  NN
#> 71  AA AA AA AT AA AA AA AA AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AT  AA  AA
#> 72  GG GG GG GG GA NN GA GA AA  GG  AA  GA  AA  GA  GG  GG  GA  NN  GG  GA  GA
#> 73  TT TG TT TT TG TT TG NN TG  TG  TG  TG  TT  TG  TT  TG  TT  NN  TT  TT  TT
#> 74  GT GG GT GG GG GG GT GG GG  GG  NN  GG  GT  GG  GG  NN  GG  GT  GT  GG  NN
#> 75  AC AC NN AC CC AC AC AC AA  CC  NN  AA  AC  AC  CC  AC  NN  AC  AC  AC  AC
#> 76  CC CT CT CT NN NN TT NN NN  CC  NN  CT  NN  NN  NN  CC  CC  NN  CT  NN  CC
#> 77  AT AT AT AT AA AA AT NN AT  NN  AA  AA  AA  AT  AA  AA  AA  AA  AA  AA  NN
#> 78  GG GG GG GG GG GG GG GG NN  GG  GG  TT  GG  GG  GG  GT  GG  NN  GG  NN  NN
#> 79  TT TT TC TT TC TC TC TT TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 80  CT CT TT TT TT CT TT TT CC  TT  TT  CT  CT  TT  TT  TT  TT  CT  CT  CT  NN
#> 81  AA NN AT AT AT AT AT AA AT  NN  AA  TT  NN  AT  AT  AT  AT  AA  NN  AA  NN
#> 82  CT TT CT TT TT TT NN TT TT  TT  TT  TT  NN  NN  TT  NN  TT  TT  TT  TT  TT
#> 83  AT TT AT TT TT AT TT TT TT  AT  TT  AT  AT  TT  AT  AT  AT  TT  TT  TT  TT
#> 84  CC CC CT CT CC NN CT CC CT  CC  TT  TT  CT  CC  CC  CT  CC  CC  CT  NN  NN
#> 85  AG AG GG AA AG NN AA AA AG  AA  AG  AG  AA  NN  AG  AG  AG  NN  GG  NN  AG
#> 86  AG AG AG AG AA AA AG NN AG  NN  NN  AA  AA  AG  AA  AA  AA  NN  AA  AA  NN
#> 87  CT CC CT CT TT CT CC CT TT  TT  CC  TT  CT  TT  TT  TT  TT  CT  NN  TT  CT
#> 88  AG NN AG NN GG AG AA AG AG  AG  NN  AA  AA  AG  AA  AG  AG  GG  NN  NN  NN
#> 89  NN CG CG CG CC CG CG NN CC  CG  CC  CC  NN  GG  CC  NN  NN  CG  NN  CC  CG
#> 90  CT CT CC CC CC CT CC CC CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  NN
#> 91  CC CC NN CT CC CT CC CT CC  NN  CC  CC  NN  TT  CC  NN  CT  CT  NN  CT  CT
#> 92  AA NN GG NN GG AA AG AG NN  NN  AG  AG  NN  GG  GG  GG  AG  AG  NN  AG  AG
#> 93  TT TT TT TT TC TT TT TT TT  TC  TC  TC  TT  TC  TT  TC  TT  TT  TT  TT  TT
#> 94  CC CC CC CC CC NN CC CC CC  NN  NN  CC  CC  CC  CC  NN  CC  CG  CC  CC  NN
#> 95  GG GG GG NN GT GG GT GT GG  GT  GG  GT  TT  GG  GG  GT  GG  GG  NN  NN  GG
#> 96  CC NN CC NN CC CC CC CC CA  CC  NN  CC  CC  CC  CA  CC  CC  CA  CC  CA  NN
#> 97  GG GG GG GG GG GG GG GG GG  GG  GG  GG  GG  GG  GG  GA  GG  GA  GG  GG  NN
#> 98  GG GG GG NN GG GT GG GT NN  GG  TT  GG  GG  GG  NN  GT  TT  GG  GG  GT  GG
#> 99  AA GG AA AA AG AG AG AG NN  GG  AG  AG  AA  GG  AA  GG  GG  NN  NN  AG  AA
#> 100 AT NN AT AT TT AT TT TT TT  AA  AA  AT  AT  NN  AA  NN  AT  NN  NN  AT  AT
#>     V22 V23 V24 V25 V26 V27 V28 V29 V30 V31 V32 V33 V34 V35 V36 V37 V38 V39 V40
#> 1    NN  NN  NN  AG  AG  NN  AG  AG  AA  AA  NN  AG  GG  AA  AG  AA  NN  NN  AG
#> 2    GG  GG  GT  GG  GG  GG  GG  GG  GT  GG  NN  GG  NN  GG  GG  GG  NN  GG  GG
#> 3    CC  CT  CC  NN  CT  NN  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 4    NN  NN  NN  AG  NN  AG  AG  AG  AA  AG  AG  AG  AG  AG  AG  NN  AA  AG  AG
#> 5    GA  NN  GG  GG  GA  GG  NN  GG  NN  NN  GA  NN  NN  GG  GA  GG  NN  NN  NN
#> 6    AA  AA  NN  AA  AA  AA  AA  AA  AA  AA  AA  AC  AA  AA  AA  AA  AC  AA  AA
#> 7    TG  TT  NN  NN  TT  TT  TT  TT  TT  NN  TT  TT  NN  TG  TT  TG  TG  TT  TG
#> 8    GA  GG  GA  NN  GG  GA  NN  GA  GG  GG  GA  NN  GA  GG  AA  GG  GG  GG  GG
#> 9    CA  NN  CC  CC  CC  CC  NN  CC  CC  CC  CC  CC  CC  CA  CA  CC  CC  CC  CC
#> 10   AA  AA  AA  AA  NN  NN  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AG  AA
#> 11   TC  CC  CC  NN  NN  CC  NN  TC  TT  TC  TC  TC  CC  TT  TT  TC  TC  CC  TC
#> 12   AG  AG  AG  NN  NN  AA  AA  AA  NN  AG  AG  AA  AG  AG  AA  AA  NN  NN  AA
#> 13   TT  TT  TC  NN  TT  TT  TT  TT  CC  TT  TC  TC  CC  TT  TC  TC  TT  TC  TT
#> 14   GG  GG  GG  TT  GG  NN  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 15   AA  AA  NN  AA  AG  NN  NN  AA  AG  AA  AG  AA  NN  NN  AA  AA  AG  GG  NN
#> 16   NN  GT  NN  GT  GG  GT  TT  GG  GG  GG  GT  GT  GG  GT  GT  GG  GT  NN  GG
#> 17   AC  NN  AA  AC  AC  AA  AA  AC  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AC
#> 18   CC  CC  CC  AC  CC  CC  AC  AC  AA  AC  CC  CC  CC  CC  CC  AC  CC  CC  AC
#> 19   NN  TT  NN  TT  NN  TT  TT  NN  NN  TG  TT  TT  TG  TT  TT  TT  TT  NN  TT
#> 20   GG  GG  GG  GG  NN  GG  NN  GA  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 21   NN  CA  NN  AA  CA  CA  CA  CA  CA  NN  CC  CC  CA  CA  NN  CC  CA  CA  NN
#> 22   CC  CC  NN  CC  CT  NN  NN  CC  CT  CC  NN  TT  CT  CC  CT  CT  TT  CC  CT
#> 23   GG  GG  GG  GG  GA  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 24   AA  AA  AA  AA  AA  NN  AA  AG  AA  AA  AG  AA  AA  AA  AA  AA  AA  AA  AA
#> 25   GG  GG  AG  AG  AG  NN  AG  AG  GG  GG  AG  GG  GG  GG  GG  GG  AG  NN  AG
#> 26   NN  TT  TT  TT  NN  TT  NN  TT  TT  TT  TT  TT  CT  TT  NN  TT  CT  CT  TT
#> 27   AA  AA  NN  AC  AA  NN  AC  AA  AA  AC  AC  AC  AA  NN  AA  AA  AC  AA  AA
#> 28   NN  NN  TC  TC  CC  TT  TC  NN  CC  CC  CC  TT  TC  CC  NN  TC  NN  CC  TC
#> 29   NN  AC  NN  NN  AA  NN  AA  AA  AC  AA  AA  AA  AC  AC  AC  AA  AA  CC  CC
#> 30   GG  GC  NN  GC  GG  GC  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 31   CG  CC  CC  NN  CG  NN  CG  CG  CC  CG  GG  CC  CC  CG  CG  CG  CC  CC  CC
#> 32   NN  AA  AC  AA  AC  AC  AC  AC  NN  NN  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 33   NN  NN  TT  TT  CC  CC  CT  CT  CT  CC  CC  NN  CT  NN  NN  CC  TT  CC  NN
#> 34   CG  NN  NN  CG  CG  NN  CG  CG  CC  CC  NN  CC  CC  CG  CG  CG  CG  GG  GG
#> 35   NN  AA  AA  AC  AA  AA  AA  AC  AA  NN  NN  AA  AC  AA  AA  AA  AA  AA  NN
#> 36   AA  NN  AG  AA  NN  AA  AA  AA  AA  AG  AA  AA  AG  AA  AG  GG  AG  AA  AA
#> 37   CG  CG  CG  GG  CG  NN  GG  CG  GG  GG  GG  CG  GG  GG  GG  CG  CG  NN  GG
#> 38   NN  GT  GT  NN  GG  GT  GT  GG  GG  GG  GT  GG  GG  GG  GG  GT  GT  GG  GG
#> 39   NN  NN  TA  NN  NN  NN  TA  TT  TT  TT  NN  TT  NN  TA  TT  TT  NN  NN  TT
#> 40   NN  NN  CG  NN  CG  CG  GG  CG  NN  NN  GG  GG  GG  GG  GG  GG  GG  CG  GG
#> 41   AG  NN  AG  NN  NN  GG  GG  AG  GG  AG  GG  GG  GG  GG  NN  AG  AG  GG  AG
#> 42   AG  GG  AG  AA  NN  GG  AG  GG  AA  NN  AG  AG  AA  AA  AG  AA  GG  GG  NN
#> 43   NN  AG  AG  NN  NN  AG  AA  AA  AA  AA  AA  AA  AA  NN  AA  AG  AA  AG  AG
#> 44   TT  TT  TT  AT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  AT  TT
#> 45   NN  GG  GG  GG  GG  GG  GG  GG  NN  GG  GG  NN  GG  GG  NN  GG  GG  NN  GG
#> 46   NN  NN  AA  CC  CC  NN  CA  CC  CA  NN  NN  CA  CC  NN  CA  CA  CA  AA  CC
#> 47   CC  CG  CG  CG  CG  CC  CG  CG  CC  CG  NN  GG  CG  NN  CG  GG  CG  GG  CG
#> 48   NN  NN  NN  NN  CC  TT  CT  CT  CT  TT  CT  CT  CC  CC  CT  CT  TT  CC  CC
#> 49   GG  GG  GG  GG  GC  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 50   NN  NN  CC  CC  CC  CC  NN  CC  CC  CC  CC  CC  CT  CC  CC  CC  CT  CC  CC
#> 51   NN  NN  TC  TC  TC  TC  TT  TT  TT  NN  TC  TC  TT  CC  NN  NN  TT  TT  NN
#> 52   AG  AG  AA  AA  AG  AA  AG  AG  AA  NN  AG  AA  AA  AG  AG  AG  AA  AA  GG
#> 53   AC  NN  AC  AA  NN  AA  AC  AC  CC  AC  AA  AC  AC  AC  AA  AA  NN  NN  AC
#> 54   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GT  GG  GG  GG
#> 55   GA  GG  GA  AA  GG  GA  GA  GA  GG  NN  GG  GG  GG  GA  NN  GG  GG  GA  GA
#> 56   AG  AG  AG  AG  AA  AA  AA  AA  NN  AG  AG  AA  AG  AG  AA  AA  AA  AA  AA
#> 57   GG  GG  GG  NN  GG  GG  GG  GG  GA  GA  GA  GG  GG  GG  GG  GG  GG  NN  GG
#> 58   CC  TC  TC  TC  CC  CC  CC  CC  TC  TC  TC  NN  CC  CC  TC  TC  NN  CC  CC
#> 59   AG  AA  GG  GG  AA  AG  NN  AG  AG  NN  AG  NN  AG  AG  GG  AG  AG  AG  GG
#> 60   NN  AC  AC  AC  AA  AA  AC  AA  CC  AA  NN  AC  AC  AC  AA  AA  AC  NN  AC
#> 61   TT  NN  TT  TT  TT  NN  NN  TT  TT  TC  TT  TT  TT  NN  TT  NN  TT  NN  NN
#> 62   CG  GG  GG  GG  CG  NN  GG  CG  CG  GG  CG  CG  GG  CG  CG  GG  GG  GG  CG
#> 63   AC  AA  AA  NN  AC  AA  AC  AC  AA  NN  AA  AC  NN  AA  NN  AA  NN  AC  AA
#> 64   TA  TT  TA  NN  TA  TA  TA  TT  TT  TT  NN  TA  TT  TA  NN  TA  TA  TT  TA
#> 65   CT  TT  CT  NN  NN  CT  CC  TT  CT  TT  NN  CT  TT  TT  NN  CT  CT  NN  TT
#> 66   AA  AC  NN  AA  AA  NN  NN  AA  AA  AC  AA  AC  AC  AA  AC  AC  NN  NN  AC
#> 67   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 68   NN  AA  NN  AT  AA  AT  AA  AT  AA  AA  AT  AT  AA  AA  AA  AT  NN  TT  AA
#> 69   AG  AG  AA  AA  AA  AA  AA  AA  AG  AG  AG  AA  AG  GG  AA  AG  AG  NN  AG
#> 70   AC  AC  AC  AA  AA  AA  NN  AA  AC  AC  AC  AA  AC  NN  AA  AC  AC  AC  AC
#> 71   AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 72   GA  NN  AA  GA  AA  NN  AA  AA  GA  AA  GA  AA  GA  AA  AA  GA  GA  GA  AA
#> 73   TT  TT  TT  NN  TG  TT  TT  TT  TG  TT  TG  TT  NN  TT  TG  TT  TT  TT  TT
#> 74   GT  NN  GG  NN  GT  NN  GT  GG  GT  GG  NN  GT  GG  GT  NN  NN  GG  NN  GT
#> 75   NN  AC  CC  NN  NN  AC  AA  AC  AA  AA  AC  AA  NN  AA  AA  AC  AC  AC  CC
#> 76   NN  NN  NN  TT  CT  NN  CC  TT  TT  NN  CT  NN  CC  NN  NN  CC  TT  NN  CC
#> 77   AT  AT  AT  NN  AA  AA  AA  AA  AA  AT  AT  AA  AT  AT  NN  AA  AA  AA  AA
#> 78   GT  GT  NN  NN  GG  GG  GG  GG  GT  GT  GT  GG  GT  TT  GG  GT  GT  GG  GT
#> 79   TT  TT  TC  NN  TC  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TC  TC  TT  TC
#> 80   NN  CT  TT  CT  NN  TT  CT  TT  TT  NN  NN  TT  CC  TT  TT  TT  CT  TT  CT
#> 81   NN  NN  AT  AT  AA  AT  AT  AT  TT  AT  AA  AT  NN  AT  AT  AT  AA  TT  AT
#> 82   CT  NN  TT  TT  CT  TT  TT  TT  CT  TT  TT  TT  TT  CT  TT  CT  TT  TT  TT
#> 83   TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 84   CC  CC  NN  CT  NN  CT  NN  CC  CT  CC  NN  TT  CT  CC  CT  TT  CT  CT  CC
#> 85   NN  AG  GG  AG  AG  AG  GG  AG  NN  AG  AG  AG  GG  AG  AG  AG  GG  GG  GG
#> 86   AG  AG  AG  AG  NN  AA  NN  AA  AA  AG  AG  AA  AG  AG  AA  NN  NN  AA  NN
#> 87   CT  TT  TT  TT  NN  TT  TT  TT  TT  TT  CT  CT  TT  CT  TT  TT  CC  TT  TT
#> 88   AG  AG  AA  AG  NN  AA  AA  AA  NN  AG  AA  AA  AA  AG  NN  AG  GG  AA  AG
#> 89   NN  CG  NN  NN  NN  NN  CG  CG  GG  CC  NN  CC  CG  CC  CC  CC  CG  GG  CC
#> 90   CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 91   CT  NN  NN  CT  CC  CT  CT  NN  CC  TT  TT  CT  NN  TT  CT  CT  TT  CT  CC
#> 92   NN  AG  NN  GG  NN  NN  NN  GG  NN  NN  AG  AG  AG  GG  NN  AG  NN  NN  NN
#> 93   TC  NN  TT  TT  TC  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TC  TT
#> 94   CC  CC  CG  CC  CC  CC  CC  CC  CC  CC  NN  CC  CC  NN  CC  CC  NN  CC  CC
#> 95   GT  NN  NN  GT  NN  GG  GG  GG  GG  GG  GG  GG  GG  TT  GG  GG  GG  GG  GT
#> 96   NN  NN  CC  CC  CC  CC  CC  CA  CC  CC  CC  CC  CC  CA  CC  CA  CA  CA  CC
#> 97   NN  GG  NN  GG  NN  GG  GG  GG  GG  GG  GG  GA  NN  GG  GG  GG  GG  NN  GG
#> 98   GG  GG  GT  GT  GG  NN  GG  GG  GG  GT  GG  GG  GG  GG  GG  GG  GT  GG  GG
#> 99   AG  AG  AG  NN  NN  AG  AA  AA  AG  GG  NN  AA  AG  AG  AG  AG  AG  AA  GG
#> 100  AT  AT  NN  NN  NN  AA  NN  AT  AA  NN  TT  AA  AT  AA  AT  AA  AA  AA  TT
#>     V41 V42 V43 V44 V45 V46 V47 V48 V49 V50 V51 V52 V53 V54 V55 V56 V57 V58 V59
#> 1    AG  AA  GG  AG  AG  GG  GG  GG  GG  GG  AG  GG  GG  GG  AG  GG  AG  GG  GG
#> 2    GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 3    CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  NN  CT  CC
#> 4    NN  AG  AG  AA  AA  AA  AA  AA  AA  AA  AA  NN  AA  AA  AA  AA  AA  AA  AA
#> 5    NN  NN  GA  NN  GA  NN  GG  GG  GG  GG  GA  NN  GA  NN  GA  GG  NN  GA  NN
#> 6    AA  AA  AA  AA  AA  AC  AC  NN  AC  AC  AC  AC  AA  AC  AC  AC  AC  CC  AA
#> 7    TT  TT  TT  NN  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 8    GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  NN  GG  GG  GG  NN  GG
#> 9    CA  CC  CC  CC  CC  CC  CC  CA  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 10   AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 11   TC  NN  TC  TC  TT  TT  TT  TT  TT  TT  TT  NN  TT  TT  TT  TT  TT  TT  TT
#> 12   NN  AA  AA  NN  AG  AG  AG  AA  AG  AG  AG  NN  AA  AG  AG  AA  AA  AA  AG
#> 13   TT  TC  TT  TC  TC  TT  TT  TT  TT  NN  TT  TT  TT  NN  TT  TT  TT  TT  TT
#> 14   GG  GG  GG  TT  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 15   AA  AG  AA  AA  GG  GG  AG  AG  GG  GG  GG  GG  NN  GG  AG  GG  AG  GG  GG
#> 16   GT  GT  GG  NN  GT  NN  TT  NN  GT  GT  GT  GG  TT  GG  GT  GT  GG  GT  GT
#> 17   AA  AA  AA  AA  CC  AA  AA  AA  AA  AC  AA  AA  NN  AC  AA  AA  AA  AA  AC
#> 18   CC  AC  AC  AC  NN  AA  AC  AA  AA  AC  AA  AC  AC  CC  NN  AA  NN  AA  AC
#> 19   TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 20   GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  GG  GG  GG
#> 21   CA  NN  CA  AA  CC  CC  CC  CC  CC  CC  CC  NN  CC  CC  CC  CC  CC  CC  CC
#> 22   CC  CT  CC  NN  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 23   GG  NN  GG  NN  GG  GG  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG
#> 24   AA  AA  GG  AA  AA  AA  AA  AA  AA  AA  AA  AA  NN  AA  AA  AA  AA  AA  AA
#> 25   GG  NN  AG  AG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  GG  GG  GG  GG
#> 26   TT  CT  CT  NN  TT  CT  CT  CC  CC  CT  TT  CT  CC  TT  CC  NN  CT  CT  CT
#> 27   AA  AC  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 28   NN  CC  TC  CC  TT  TT  TT  TT  TT  TT  TC  TT  TC  TT  TT  TT  NN  NN  TT
#> 29   NN  AA  AC  AC  AC  AA  AC  AC  AA  AA  AA  AA  AA  AC  AA  AA  NN  NN  AA
#> 30   NN  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 31   CC  CG  GG  CG  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 32   AA  AA  AA  AA  AA  AC  AA  AC  AC  AC  AC  AC  AA  AC  AC  AC  AC  NN  AC
#> 33   CT  CC  CT  CT  CC  CT  CC  NN  CC  CT  CC  CT  CC  CC  CT  CT  CT  CC  CC
#> 34   CC  NN  NN  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 35   AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  NN  AA  AA  AA
#> 36   AA  AA  AA  AG  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 37   CG  CG  CG  NN  GG  GG  GG  GG  NN  GG  GG  CG  GG  GG  GG  GG  GG  GG  GG
#> 38   GG  GT  GG  GG  GG  GT  GT  GT  GT  GT  GT  GG  GT  GT  GG  GT  GT  GG  NN
#> 39   TT  TT  TT  TT  TT  TA  TT  NN  NN  TT  TA  NN  NN  TT  TT  NN  TT  NN  TT
#> 40   NN  GG  GG  CG  GG  CG  GG  CG  GG  CG  GG  NN  GG  CG  CG  GG  NN  GG  GG
#> 41   AG  NN  GG  AG  AA  AA  AA  AA  AA  AA  NN  AA  AA  NN  NN  AA  NN  AA  AA
#> 42   AG  AG  NN  GG  AG  GG  GG  GG  GG  GG  GG  NN  GG  GG  AG  GG  GG  GG  GG
#> 43   AA  AG  AA  AG  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 44   TT  TT  TT  AT  AT  AT  AA  AT  AT  AT  AT  TT  AT  AA  NN  TT  TT  NN  AT
#> 45   GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 46   NN  CA  CC  AA  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 47   GG  CG  CC  NN  CC  CC  CC  CC  CC  NN  CG  CG  NN  CG  CC  CC  CC  CC  CC
#> 48   CT  NN  CT  NN  CT  CC  CT  CT  TT  CT  TT  NN  CC  TT  CT  NN  NN  TT  TT
#> 49   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 50   NN  CC  CC  CT  CC  CT  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  NN  CC  CC
#> 51   TC  TC  TT  TC  NN  TT  TC  TC  NN  NN  TT  TT  TT  TC  TT  NN  NN  TT  NN
#> 52   AG  NN  GG  AA  AG  AG  NN  AA  GG  AA  AG  NN  AG  NN  NN  AG  AG  AG  NN
#> 53   NN  AC  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  NN  AA
#> 54   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 55   AA  NN  GG  GG  AA  GA  GA  GA  NN  GA  GG  GA  GG  GA  GG  GA  GA  NN  GA
#> 56   NN  AA  NN  NN  AG  AG  AG  AA  AG  AG  AG  AG  AA  AG  AG  AA  NN  AA  AG
#> 57   NN  NN  GG  GG  NN  NN  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  NN
#> 58   TC  TC  TC  CC  CC  TT  TT  TC  NN  TC  TC  TC  TC  NN  CC  CC  CC  NN  CC
#> 59   GG  GG  GG  AG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  AG  GG
#> 60   AA  AA  AC  AC  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 61   TT  TT  TT  NN  NN  TT  TT  NN  NN  TT  TT  TT  NN  TT  TT  TT  TT  TT  TT
#> 62   CG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG
#> 63   NN  AC  AA  NN  NN  NN  AC  CC  AC  CC  NN  NN  NN  AC  CC  CC  NN  NN  CC
#> 64   TT  NN  TT  TA  TA  AA  NN  TA  AA  TA  AA  TA  AA  TT  AA  TT  TA  TT  NN
#> 65   CC  CT  TT  CC  CC  CC  CT  CC  CC  CC  CC  CC  CT  CT  CC  CT  CT  TT  CT
#> 66   NN  CC  AC  AC  CC  AC  AC  AA  AA  CC  AC  CC  CC  CC  CC  AC  NN  AC  AC
#> 67   NN  GG  GG  GG  GA  AA  GG  NN  GG  GA  GG  GA  GG  GA  GG  NN  GA  GG  GG
#> 68   AA  AA  AA  AA  TT  TT  NN  NN  TT  AT  AT  AA  TT  TT  NN  AT  AT  TT  AT
#> 69   AA  AA  GG  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 70   AA  AA  NN  AA  AA  AA  AA  AA  AA  NN  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 71   AA  AA  AA  AA  AA  AT  AA  AT  NN  AA  AA  AT  AA  AA  AA  AA  AT  AA  AA
#> 72   AA  GA  GA  GA  GG  GA  GA  GG  GG  GA  AA  NN  GG  GG  NN  GG  AA  GA  GA
#> 73   TG  NN  TT  TG  TT  TT  TG  TG  TT  TG  TG  TT  NN  TG  TT  TG  GG  TT  NN
#> 74   GG  GT  GT  GT  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 75   AC  AC  CC  AC  AA  AC  AC  AA  AC  AC  AA  AA  AA  AC  AA  AC  NN  AC  AC
#> 76   CC  CC  CC  CC  NN  NN  NN  TT  TT  TT  TT  NN  NN  NN  NN  TT  NN  TT  TT
#> 77   AT  NN  AA  AA  AT  AT  AT  AA  AT  AT  AT  AT  AA  AT  AT  NN  AA  NN  AT
#> 78   GG  GG  TT  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 79   TT  TC  TT  TT  TT  TT  TT  TC  TC  TT  TT  TC  TC  TT  TT  TT  TC  TC  TC
#> 80   CT  CT  CT  CT  TT  CC  CT  NN  CT  CT  CT  TT  CT  CT  CT  CT  CC  CT  CT
#> 81   AA  AT  AT  NN  TT  NN  AT  AT  AA  AA  AA  AT  AT  AT  AT  AT  TT  AT  AA
#> 82   TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 83   TT  TT  TT  NN  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  NN
#> 84   NN  CT  CC  CT  NN  NN  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 85   GG  NN  AG  NN  AA  AA  AA  AA  AA  AA  NN  AG  AA  AA  AG  NN  AA  AA  AA
#> 86   AG  AA  AA  AA  AG  NN  AG  NN  NN  AG  AG  AG  NN  AG  AG  NN  AA  NN  AG
#> 87   CC  TT  CT  CT  TT  TT  NN  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 88   NN  AG  AG  AA  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 89   CG  NN  CG  NN  CC  NN  CC  NN  GG  CC  CG  GG  NN  CC  CG  NN  CC  GG  NN
#> 90   CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 91   CT  CT  NN  CC  CC  CC  CC  CC  CC  CC  CC  CC  CT  CC  CC  CT  CC  CT  CC
#> 92   NN  AG  AG  GG  AG  NN  GG  AG  NN  AA  AA  NN  AG  NN  AA  AG  NN  AG  AG
#> 93   TT  TT  TT  TT  NN  TC  TT  TT  TT  TC  TT  TT  TT  TT  NN  TT  TC  TT  TT
#> 94   CC  NN  CG  CC  CG  CG  CC  NN  CG  CC  NN  NN  CG  CG  CC  GG  CG  NN  CG
#> 95   GG  TT  GT  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 96   NN  CA  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 97   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 98   GG  GG  GG  GG  GG  GG  GG  GT  GT  GT  GG  NN  GT  GT  GT  GG  GT  GT  TT
#> 99   NN  AG  AG  AG  GG  GG  GG  GG  GG  GG  NN  GG  GG  GG  NN  GG  GG  GG  GG
#> 100  AA  AA  AT  AA  NN  NN  AT  NN  TT  AA  NN  AA  NN  AT  AT  AT  AA  AT  AA
#>     V60 V61 V62 V63 V64 V65 V66 V67 V68 V69 V70 V71 V72 V73 V74 V75 V76 V77 V78
#> 1    AG  AA  AA  GG  AA  AG  NN  AG  AA  GG  AG  AG  GG  AA  AG  AG  AG  AA  AA
#> 2    GG  GG  GG  GT  GG  GG  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  NN  GG  GT
#> 3    CC  CC  CC  CC  CT  CC  CC  CC  NN  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 4    AA  AA  AG  AA  AA  AG  NN  AG  AA  AG  GG  AG  AA  AA  AG  AA  AA  AG  AA
#> 5    GG  GG  NN  NN  NN  GA  GG  NN  NN  GA  GA  GG  GG  NN  NN  GG  NN  GG  GG
#> 6    AA  AC  AC  AA  AC  NN  AA  AA  NN  AA  AA  AA  AA  AA  AC  AA  AA  AA  NN
#> 7    TT  TT  NN  TT  TT  NN  TG  TG  TT  TT  TG  NN  TT  TT  TT  TG  TT  TG  TT
#> 8    GG  GG  GA  GA  GG  NN  GG  GA  GA  GG  GA  NN  AA  GA  GA  GA  GA  NN  GG
#> 9    NN  CC  CC  CA  CC  NN  CA  CA  CA  CA  CC  CC  CA  NN  CA  CC  CC  CA  NN
#> 10   AA  AA  AA  AA  AA  AA  NN  AA  AG  NN  AA  AG  AA  AA  AA  AA  AA  AA  AA
#> 11   TT  TT  NN  TC  CC  TC  TT  TC  NN  TC  CC  TC  TT  TC  TC  NN  NN  TT  TC
#> 12   AA  AG  NN  AA  AG  AA  AA  AA  NN  AG  AA  AA  NN  AA  NN  NN  AA  AG  AG
#> 13   TT  TT  TT  TT  TT  TT  TT  TT  TC  TC  TC  TC  TT  TC  TT  NN  TT  CC  NN
#> 14   GG  GG  GG  GG  GG  NN  GG  GG  GT  GG  GG  GT  GG  GT  GG  GG  GT  GT  NN
#> 15   GG  GG  AA  AA  AG  AG  AA  AG  AA  AA  AG  AG  NN  AA  AA  AG  AA  AA  AA
#> 16   GT  GG  GG  GG  GT  GG  NN  GT  GT  GG  TT  GT  NN  GT  GG  GG  GT  GT  NN
#> 17   AC  AA  AC  AA  AC  AC  NN  AA  AA  AA  AA  AA  NN  AC  AA  AA  AA  AC  AA
#> 18   AA  AC  AC  CC  CC  CC  NN  CC  CC  CC  CC  CC  CC  AC  CC  CC  AC  AA  AC
#> 19   TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TG  TT  TG  NN  TT  TT
#> 20   GG  GG  GG  GG  GA  GA  GG  GG  GG  GG  GG  GG  GA  GG  GG  NN  GG  GG  GG
#> 21   CC  CC  CA  CA  CC  CA  CA  CC  CC  CA  CA  NN  CC  CA  CA  NN  CA  CC  CA
#> 22   CC  CC  CC  CC  CT  CC  NN  NN  CT  CC  CC  CC  CC  NN  CT  CC  CC  NN  CC
#> 23   GG  GG  GG  GA  GG  GG  GG  GA  GG  GG  NN  GA  GG  NN  GG  NN  GA  GG  GG
#> 24   AA  AA  NN  AA  AG  AA  AA  AA  AA  AA  AA  AA  AG  AA  AA  NN  NN  AA  AA
#> 25   GG  GG  AG  AA  AG  AA  AG  NN  GG  NN  NN  NN  AA  AG  GG  AG  NN  AG  AG
#> 26   NN  CT  CT  NN  TT  TT  TT  CT  CT  TT  CT  CT  CT  CT  CT  CT  CT  TT  CT
#> 27   AA  AA  AA  NN  AA  AA  AC  NN  AC  NN  AA  AA  AA  AC  AA  AC  AA  AA  AA
#> 28   TC  TT  TT  TT  CC  CC  TC  CC  NN  TC  TC  TT  NN  NN  TC  TC  TC  NN  TC
#> 29   AC  AA  NN  NN  AC  AA  AC  AA  NN  NN  AC  NN  NN  AA  AA  NN  NN  AA  AC
#> 30   GG  GG  GG  GG  GG  GC  GC  GG  GG  NN  GG  GG  GC  GG  GG  GG  GC  GG  GG
#> 31   CC  CC  NN  CC  GG  CG  NN  CG  CC  CC  GG  NN  CG  CC  CG  CG  CG  CG  NN
#> 32   NN  AC  AA  AC  AC  AA  AA  AA  AC  NN  AC  AA  AA  AA  AA  AC  AA  NN  AC
#> 33   CT  NN  CC  CT  CT  CC  CC  CT  CC  CT  CC  NN  CC  CC  TT  NN  CC  CC  TT
#> 34   CC  CC  CC  CG  CC  CC  CG  NN  CC  CG  CC  CC  GG  CC  CC  NN  CC  CG  CG
#> 35   AA  AA  NN  AA  AA  AA  NN  AA  AC  AA  AA  AA  AC  NN  AC  AA  AA  AA  AA
#> 36   NN  AA  AA  AA  AG  AA  AG  AA  AA  AA  AG  AA  AA  AG  AA  AA  AA  AG  AA
#> 37   GG  GG  CC  GG  CG  GG  CG  CG  CG  GG  CG  GG  NN  CG  GG  GG  CG  GG  CC
#> 38   GT  TT  GT  GT  TT  GT  NN  NN  NN  GG  GG  GT  GT  GT  GT  GG  GG  GT  GG
#> 39   TA  TT  NN  NN  TA  TT  TT  TT  NN  TT  TT  NN  NN  TT  NN  TT  TT  TT  TA
#> 40   GG  GG  GG  GG  GG  GG  GG  CG  NN  GG  GG  GG  GG  CG  GG  NN  GG  NN  GG
#> 41   AA  AA  AG  AG  AA  NN  AG  NN  NN  GG  AG  AG  GG  AG  AG  NN  NN  GG  GG
#> 42   GG  GG  AG  AA  AG  GG  AG  GG  NN  AG  AG  AG  AG  AA  AG  AA  AG  AA  AG
#> 43   AA  AA  AA  AA  AA  AA  AA  AA  AG  AA  AA  AA  AG  AA  AA  AA  AG  AA  AA
#> 44   TT  TT  NN  TT  TT  NN  AA  NN  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 45   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  NN  GG  NN
#> 46   CC  NN  CC  NN  AA  CC  CC  CA  CC  NN  CA  NN  CC  CA  NN  NN  CC  CC  CC
#> 47   CC  CC  CC  CC  CG  CC  CC  CG  CG  CC  CC  CC  CC  CG  CC  CC  NN  CG  CG
#> 48   TT  CC  TT  TT  CC  CT  CT  TT  CT  NN  CT  NN  TT  CT  CT  NN  NN  CC  CT
#> 49   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 50   CC  CC  CC  CC  CC  CC  CT  NN  CC  CC  CC  CC  NN  NN  NN  CC  CC  NN  CC
#> 51   TT  TC  NN  NN  TC  TC  NN  NN  TC  NN  NN  TT  TC  NN  TC  TC  TT  TC  TT
#> 52   AG  AG  AA  AG  AA  AG  AA  AA  GG  AA  AG  AG  AA  AA  AA  AA  AA  AA  AA
#> 53   AA  AA  AA  AA  AA  AA  AA  AA  AC  NN  AC  AA  CC  AA  AA  AA  AA  AA  AA
#> 54   NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  GG  GT  GG
#> 55   GA  AA  GG  GG  AA  GG  GA  GG  GG  GG  GA  GG  GG  GG  GA  NN  NN  GG  GG
#> 56   AA  AG  AA  NN  AG  AA  AA  NN  AA  AG  AA  AA  NN  NN  AG  AA  AA  AG  AG
#> 57   NN  GG  NN  GG  GG  GG  GG  GG  NN  GG  NN  GG  NN  GG  GG  GA  GG  GG  GG
#> 58   NN  TT  NN  CC  CC  CC  CC  NN  CC  NN  TC  NN  CC  CC  TC  NN  CC  CC  CC
#> 59   AG  NN  NN  AG  NN  NN  AG  AG  GG  GG  AG  AG  AG  AA  AG  AG  GG  AG  AG
#> 60   AA  AA  AA  NN  AC  AC  NN  AA  AC  AA  AA  AC  AC  CC  AA  CC  AC  AA  AC
#> 61   TT  TT  TT  TT  TT  TT  TT  TT  NN  NN  TT  TC  TT  TT  TT  NN  TT  TT  TT
#> 62   GG  GG  GG  GG  CG  GG  CG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 63   AC  CC  NN  AC  AA  AA  AC  AC  NN  AC  AC  NN  NN  AA  AA  NN  NN  AC  NN
#> 64   NN  TA  TT  TA  TT  AA  TT  TA  TT  TA  TT  TA  AA  TA  AA  TT  NN  TA  NN
#> 65   CC  CT  CC  CT  CT  CT  NN  CT  CC  TT  CC  CT  CC  CC  NN  CT  CC  CT  CT
#> 66   CC  AC  NN  AA  AC  AC  AC  NN  NN  AA  NN  AC  AA  NN  AA  CC  AC  CC  AA
#> 67   GA  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 68   NN  TT  NN  AA  AT  AA  AT  AT  AA  AA  NN  AA  NN  AT  AT  AA  AA  AT  AT
#> 69   AA  AA  AA  AA  AA  AA  AA  AG  AA  AA  AA  AA  NN  AA  AA  AA  AA  AG  AA
#> 70   AA  AA  AC  AA  AA  NN  AC  AC  AA  AC  AA  AA  AA  AA  AA  AC  AA  AC  AC
#> 71   AT  AT  AA  AA  AA  AA  AA  AA  AT  AA  AT  AA  AA  AA  AA  AA  AA  AT  AA
#> 72   GA  GA  GA  AA  GA  AA  AA  NN  AA  GA  AA  NN  GA  AA  GA  GA  AA  GA  AA
#> 73   TG  TT  TT  TT  TT  NN  TT  TG  TT  TG  TT  TT  TT  TT  TG  TG  TT  TT  TT
#> 74   GG  GG  NN  GG  GT  GG  NN  GT  GT  NN  GT  NN  GG  GG  NN  GG  NN  GG  GT
#> 75   AC  AC  AA  CC  AA  NN  CC  AC  AC  AC  AA  AC  NN  AC  AC  CC  AA  AC  AC
#> 76   NN  TT  NN  NN  CC  CC  NN  CC  CC  NN  NN  NN  CC  CT  NN  NN  NN  CC  CC
#> 77   AA  AT  AA  NN  AT  AA  AA  AA  AA  AT  AA  NN  AA  AA  AT  NN  NN  AT  AT
#> 78   GG  GG  GG  GG  GG  GG  NN  GT  GG  GG  GG  GG  GG  GG  GG  GG  GG  GT  GG
#> 79   TT  TC  TT  TT  TC  TT  NN  TC  TT  TT  TT  TT  TT  TT  NN  TT  TT  TT  TT
#> 80   NN  CT  TT  TT  TT  CT  TT  CT  TT  CT  TT  CT  CT  CT  CT  CT  CT  CT  CT
#> 81   AT  AT  TT  AA  AT  AA  AA  NN  AT  AT  AA  AT  AT  AA  NN  AA  AT  AT  AA
#> 82   TT  TT  NN  NN  CT  TT  TT  TT  TT  TT  TT  TT  TT  TT  NN  NN  CT  TT  TT
#> 83   TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 84   TT  TT  CC  TT  CC  CT  CT  NN  NN  CC  NN  NN  CC  CT  CT  CC  CC  CT  CT
#> 85   AG  NN  GG  GG  AG  GG  GG  NN  GG  NN  GG  AG  AG  AG  GG  AG  GG  AG  AG
#> 86   AA  AG  AA  AA  AG  AA  AA  AA  AA  AG  AA  NN  AA  AA  AG  NN  NN  AG  AG
#> 87   TT  TT  CT  TT  TT  NN  CT  NN  CT  CT  CT  TT  TT  TT  TT  CT  TT  CT  TT
#> 88   NN  NN  NN  AA  AG  AG  AA  AG  NN  NN  AG  AG  AG  AG  NN  NN  NN  GG  AA
#> 89   GG  GG  NN  NN  CC  NN  CC  CG  NN  CC  CC  CG  NN  CG  CG  CC  CG  CG  GG
#> 90   CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 91   CC  NN  CT  TT  CC  CT  TT  NN  CC  CT  CC  CC  CT  TT  CT  NN  NN  CT  CT
#> 92   NN  AA  NN  NN  AG  GG  AG  NN  NN  NN  GG  GG  NN  NN  NN  NN  AG  NN  AG
#> 93   TC  NN  TC  TT  TT  NN  TT  NN  NN  TC  TT  TT  TT  TC  TT  TT  TC  TT  TT
#> 94   CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  NN  NN  CC  CC  CC  CC  CC  CC
#> 95   GG  GG  NN  GT  GT  GG  GG  GT  NN  GG  GG  GG  NN  GT  GT  NN  GT  GG  GG
#> 96   CC  CC  CA  CC  CC  CC  CC  CA  CC  CC  CC  CC  CC  CC  CC  CC  CC  CA  CC
#> 97   GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  GG  NN  GG  GG  GG  GA
#> 98   GT  TT  GT  GG  GG  GT  GG  NN  GT  GG  NN  GG  GT  GG  GG  GG  GT  GT  GG
#> 99   GG  NN  AG  GG  AG  AG  AG  AA  NN  AG  AA  GG  AG  AA  NN  AA  AG  GG  AG
#> 100  NN  AA  NN  AT  AT  AA  AT  NN  NN  AT  AA  AA  AT  AT  NN  AT  AT  AT  AA
#>     V79 V80 V81 V82 V83 V84 V85 V86 V87 V88 V89 V90 V91 V92 V93 V94 V95 V96 V97
#> 1    NN  AA  AG  AG  AG  AA  GG  AA  AA  AG  AG  AG  AA  AG  AG  GG  AG  NN  GG
#> 2    GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 3    CT  CC  CC  CC  CT  CC  CC  CT  CT  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 4    AA  GG  AA  AG  AG  GG  AA  AG  GG  AG  AG  AG  AG  AG  NN  NN  AG  AG  AG
#> 5    GG  GG  NN  AA  GA  GG  GA  GG  NN  AA  GA  GA  NN  NN  NN  AA  GA  GA  NN
#> 6    AA  AA  AA  AA  AA  AA  AA  AC  AC  AC  AC  CC  AC  NN  AA  AC  AC  AA  AC
#> 7    TT  TT  TG  TT  TT  TG  TG  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  NN
#> 8    GG  AA  NN  GG  GG  GA  GG  NN  GG  GA  GA  NN  GG  NN  NN  GG  GA  GA  NN
#> 9    CA  NN  CC  CA  CC  CC  CC  CC  CC  CA  CC  NN  CC  CA  CC  CA  CA  CC  NN
#> 10   AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AG  NN  AG  AA  AA  AA  AA  AA
#> 11   TC  TC  NN  CC  CC  TC  TT  TC  TC  TT  TT  TC  NN  TC  TC  TC  TC  TC  TC
#> 12   AG  AA  NN  NN  AG  AA  AA  AG  AA  AG  AG  AA  NN  AA  AG  NN  NN  NN  AG
#> 13   TC  TT  TT  NN  NN  TT  TT  TT  TC  TT  TT  TT  TT  TT  TT  TC  TC  TC  TT
#> 14   GG  GT  GG  GT  GG  GT  GG  GG  GT  GT  GT  GT  GT  GG  GG  GT  GG  TT  GT
#> 15   AA  AG  AG  NN  AA  NN  AG  AG  AG  GG  NN  AG  AG  AG  NN  AG  AG  AG  AG
#> 16   GG  NN  GT  GT  GG  GT  GT  GG  GG  GT  GT  GG  GT  NN  GG  NN  GG  GG  GG
#> 17   AA  AA  AA  AC  NN  AC  AA  AC  AA  AA  AC  AA  NN  AA  AA  AA  AC  NN  AA
#> 18   AC  CC  CC  CC  CC  CC  CC  AA  CC  AA  CC  CC  CC  AA  CC  AC  AC  CC  AC
#> 19   TT  TT  TG  TT  TT  TT  TT  TG  TT  TT  TT  TT  TT  TT  TT  TG  TT  TT  TT
#> 20   GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  GA  GG  GG  GG  GG  GG  GG  GG  GG
#> 21   CA  CA  NN  NN  CC  CC  CA  CC  CA  CA  NN  CC  NN  CA  CC  CA  CC  CA  NN
#> 22   CT  CC  CT  CC  CT  CC  CC  CC  TT  CT  CT  CT  CT  CC  NN  CC  CT  CT  CT
#> 23   GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GA  GG  GG  GG  GG  GG  GG  GG
#> 24   AA  AA  NN  NN  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AG  AG  AA  NN  AG
#> 25   AA  AG  GG  GG  AG  AG  AG  GG  AG  AG  AA  GG  AG  AG  GG  AA  GG  AG  AG
#> 26   TT  TT  CT  TT  TT  NN  CT  CT  CT  CT  TT  CT  TT  CT  CT  TT  CT  NN  NN
#> 27   NN  AA  AA  AA  AC  AC  AA  AA  AA  AC  AA  AA  AC  AA  AA  AA  AA  AA  AA
#> 28   TC  TC  TT  TC  TC  CC  CC  TT  TC  TC  TT  TT  TT  TC  TC  NN  TT  TT  TC
#> 29   CC  AA  NN  AC  AA  AA  AA  AA  AA  AA  AC  AA  NN  AA  AA  NN  AA  AA  AA
#> 30   GG  GG  NN  GC  GG  GG  GG  GG  GG  GG  NN  GG  GG  GC  GC  NN  GG  NN  GG
#> 31   CC  CG  CG  NN  CC  CC  CG  GG  CG  CG  CC  CC  CC  CG  CG  GG  CC  GG  CG
#> 32   AA  AA  AA  AA  NN  NN  AC  AC  AA  AA  AC  AC  AA  AA  AC  AC  AC  AA  NN
#> 33   CT  TT  CC  TT  CC  CT  CT  CC  TT  CC  TT  CC  CT  CC  NN  CT  CC  CC  CC
#> 34   CC  CG  GG  CC  CC  CG  GG  CG  CC  CG  CC  CC  NN  CC  CC  CC  CC  NN  CG
#> 35   AA  AA  NN  AA  AA  AA  AA  AC  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  NN
#> 36   NN  NN  AG  AA  AA  AG  AG  AA  AA  AA  AG  AA  AA  AA  AG  AA  NN  AA  NN
#> 37   CG  CG  NN  GG  GG  CG  CC  GG  GG  GG  GG  CG  CG  GG  CG  CG  CG  CG  NN
#> 38   NN  GT  GT  GT  GT  TT  NN  GG  GG  GT  GG  GT  GG  GG  TT  TT  GG  GT  GG
#> 39   TT  TT  TT  NN  TA  TT  TT  TT  TA  TT  TT  TT  NN  NN  TA  NN  TT  TA  TT
#> 40   CG  GG  NN  NN  CG  CG  CG  CG  CG  GG  CG  GG  NN  CG  CG  NN  CC  GG  GG
#> 41   AG  GG  NN  GG  NN  GG  AG  AA  NN  AA  NN  GG  NN  AG  AA  AA  AG  NN  NN
#> 42   AG  AG  NN  AG  AA  AG  AG  AG  AG  GG  GG  AG  GG  AG  AA  GG  NN  AG  AG
#> 43   AA  AA  AA  AG  AA  AA  AG  NN  AA  AG  AG  NN  AG  AA  AG  AA  AG  AG  AA
#> 44   AT  AT  TT  AT  TT  TT  TT  TT  AT  TT  NN  TT  AT  AT  NN  TT  NN  TT  TT
#> 45   GG  GG  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 46   NN  CC  CC  NN  CC  CC  AA  CC  CA  CC  NN  CC  CC  CA  AA  CA  CC  CA  CC
#> 47   GG  CC  CC  CC  CC  CG  GG  CC  CC  CG  CC  CC  CC  CC  CG  CC  CC  CC  GG
#> 48   CT  TT  CT  CT  NN  CC  NN  CT  NN  TT  CT  CT  CT  TT  CT  CT  CT  TT  CT
#> 49   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  GC  GG  GG  NN  GG
#> 50   CC  CC  CC  CT  CC  CC  CC  CC  CC  CT  CC  CC  CC  CT  CC  CC  CC  NN  CC
#> 51   NN  CC  TT  TT  NN  NN  TT  CC  TT  CC  TC  TT  TC  TC  TC  TT  TC  NN  TC
#> 52   AA  AG  AA  AA  AA  AA  AA  GG  GG  AA  AG  AG  NN  AG  AG  AA  AA  AG  AG
#> 53   AC  AC  AA  AA  AC  AC  AA  AA  CC  AA  AC  NN  NN  AA  AC  AC  AC  AC  NN
#> 54   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 55   NN  GA  GG  NN  GG  GG  GA  AA  GA  GA  GA  GG  GG  GA  AA  GG  NN  NN  AA
#> 56   AG  AA  AA  AA  AG  NN  AA  AG  AA  AG  AG  AA  AG  AA  AG  AA  AA  AA  AG
#> 57   GG  GG  NN  GG  GG  GG  NN  GG  GG  GG  GG  GG  NN  GG  NN  GG  NN  GG  GG
#> 58   TC  TC  CC  CC  TC  CC  NN  NN  TT  TC  NN  NN  TC  TC  CC  CC  CC  NN  CC
#> 59   GG  AG  AG  NN  AA  AG  AA  GG  AG  GG  AG  GG  GG  GG  AG  AG  AG  AA  AG
#> 60   AC  AA  AC  AC  NN  AA  AA  AC  AA  NN  AC  AC  AC  AC  AA  AA  AA  AC  AC
#> 61   TT  NN  TT  TT  TC  TT  TT  TT  TC  TT  TT  TT  NN  TT  TT  TT  NN  TT  TT
#> 62   GG  GG  CG  GG  GG  GG  GG  GG  NN  GG  GG  GG  GG  CG  CG  CG  GG  GG  GG
#> 63   CC  AA  AC  AA  AA  AC  NN  AC  AC  AC  AA  CC  AA  AC  AC  AA  NN  AA  AA
#> 64   TT  TA  TA  TA  TA  TT  TT  AA  TT  TT  TA  TA  TT  TT  TA  TA  TA  TA  TT
#> 65   CC  TT  CT  CC  CC  TT  CC  TT  CC  TT  CC  CT  CC  CT  NN  CC  CC  CC  CT
#> 66   AA  AC  NN  NN  AA  CC  AC  AC  AC  AA  AA  AC  NN  AA  AA  AA  CC  CC  AC
#> 67   GG  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  NN  GA  GG
#> 68   TT  AA  AT  NN  NN  NN  AA  AA  TT  AT  AA  AT  AT  AA  AA  AT  AT  NN  TT
#> 69   AA  NN  AG  AA  AA  AG  NN  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 70   AA  AC  AC  NN  AA  AC  AC  AA  AC  NN  AC  AA  AA  AC  AC  AA  AA  AC  NN
#> 71   AT  AT  AA  AA  AA  AT  AA  AA  AA  AA  AA  AA  AA  AA  AA  AT  AA  AA  AA
#> 72   AA  GA  AA  GA  AA  AA  GA  GA  GG  GA  GA  NN  GA  GG  NN  GG  GA  GG  GG
#> 73   NN  TG  TT  TT  TT  TT  NN  TG  TT  TT  TG  TG  TT  GG  TT  TT  TG  TT  TT
#> 74   GT  GG  NN  GT  GT  GG  GG  GT  GG  GG  GG  GG  GG  NN  GG  GG  GG  GT  GG
#> 75   AA  CC  AC  NN  NN  AC  AC  AC  AC  AC  AC  AC  AC  AC  CC  CC  AC  CC  CC
#> 76   CC  NN  CC  CT  NN  CT  NN  CT  CC  CC  CT  CT  CT  CT  CC  CC  TT  CT  TT
#> 77   AT  AA  AA  NN  AT  AA  AA  AT  AA  AT  AT  AA  NN  AA  AT  AA  AA  AA  AT
#> 78   GG  NN  GT  GG  GG  GT  GT  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 79   TT  TT  TT  CC  TT  TT  TT  TC  TT  TT  TT  NN  TT  TT  TT  NN  NN  TT  TC
#> 80   NN  TT  TT  TT  TT  CT  TT  TT  CT  NN  CT  CT  CT  CT  TT  CT  TT  NN  CT
#> 81   AT  AA  AT  AA  AT  AT  AT  AT  AA  AT  NN  AT  AT  AT  AA  AA  AT  AT  TT
#> 82   TT  TT  TT  NN  TT  TT  TT  TT  TT  TT  CT  TT  TT  NN  CT  TT  CT  CT  TT
#> 83   TT  TT  TT  TT  TT  TT  AT  TT  TT  TT  TT  TT  TT  TT  AT  TT  TT  TT  TT
#> 84   CC  CC  CC  NN  CC  CT  CT  CT  TT  CC  CT  CT  CT  CT  CT  TT  CT  NN  CC
#> 85   GG  GG  AG  AG  GG  AA  NN  GG  GG  AG  GG  AG  AA  AG  AG  AA  AG  AG  AG
#> 86   AG  AA  NN  AA  AG  AA  NN  AG  AA  AG  AG  NN  AG  AA  AG  AA  NN  AA  AG
#> 87   CT  CT  CT  CT  TT  TT  CT  CT  TT  NN  TT  TT  CT  CT  CC  TT  CT  CT  CT
#> 88   AA  NN  AG  AG  AG  AG  NN  NN  AG  AG  AG  AG  AA  AA  AG  AA  NN  AG  NN
#> 89   CC  CC  CG  CG  NN  GG  CG  GG  CC  CC  CG  NN  GG  CG  NN  CG  NN  CG  GG
#> 90   CC  CT  CC  CC  CC  CC  NN  CC  CC  CC  CT  NN  CC  CC  CT  CC  CC  CC  CC
#> 91   TT  TT  TT  CT  NN  CT  CT  CT  CT  CC  TT  TT  CT  CC  CC  CC  CC  CT  NN
#> 92   GG  GG  NN  NN  GG  NN  GG  NN  GG  GG  GG  AG  NN  NN  NN  AG  GG  AG  GG
#> 93   TT  TC  NN  NN  TC  TT  TC  TT  TT  TT  TC  TC  TT  TT  TC  TT  TT  TC  NN
#> 94   CC  CC  CC  CG  CC  CC  CC  CC  CC  CC  CC  CG  CC  CG  CC  CC  NN  CG  CC
#> 95   GG  GG  NN  GG  GG  GG  GG  GG  GT  GG  GG  GG  GT  GT  GT  GG  GG  NN  GG
#> 96   NN  CC  CC  CC  CC  CA  CC  CC  CC  CC  NN  CC  CC  CC  CC  NN  CA  CC  CC
#> 97   NN  GG  GG  GG  GG  GG  GA  GG  GG  GG  GG  GG  GG  GG  GA  GG  GG  GG  GG
#> 98   GG  GG  GT  GG  GG  GT  GT  GT  GG  GG  GT  GT  NN  GG  GG  GT  NN  GT  GG
#> 99   AA  AG  GG  NN  NN  GG  AA  AG  GG  AA  NN  NN  GG  AG  AG  GG  NN  NN  AA
#> 100  AT  AT  AA  AT  AT  AA  AA  AA  AT  AA  AT  AT  AT  TT  TT  TT  NN  TT  TT
#>     V98 V99 V100
#> 1    NN  GG   AG
#> 2    GG  GG   GG
#> 3    CC  CC   CC
#> 4    AA  AA   AA
#> 5    GG  NN   GG
#> 6    AA  AA   AA
#> 7    TT  TT   TT
#> 8    GA  GG   GG
#> 9    CC  CA   AA
#> 10   AG  AA   AG
#> 11   TC  TC   TT
#> 12   AA  AG   NN
#> 13   TT  NN   NN
#> 14   GG  GT   GT
#> 15   AG  AA   AA
#> 16   TT  NN   GG
#> 17   AC  AA   AA
#> 18   AA  AC   AC
#> 19   TT  TT   TT
#> 20   GG  GA   GG
#> 21   CA  AA   CC
#> 22   CT  CC   CC
#> 23   GG  GG   GG
#> 24   AA  AA   AA
#> 25   AA  AG   GG
#> 26   NN  TT   TT
#> 27   NN  AA   AA
#> 28   TC  TC   TC
#> 29   AA  AC   AC
#> 30   GC  GG   GC
#> 31   CG  CC   CG
#> 32   AA  NN   AC
#> 33   CT  NN   TT
#> 34   NN  CC   NN
#> 35   AA  AA   AC
#> 36   AA  AA   AA
#> 37   CG  NN   GG
#> 38   GG  GG   NN
#> 39   NN  TT   NN
#> 40   NN  CC   NN
#> 41   GG  AG   AG
#> 42   GG  AG   AG
#> 43   AA  AA   AA
#> 44   AT  AT   TT
#> 45   GT  GG   GT
#> 46   AA  NN   CA
#> 47   CC  CC   CG
#> 48   NN  CC   TT
#> 49   GC  GG   GG
#> 50   CC  CC   CC
#> 51   TC  TT   TT
#> 52   AA  AA   AA
#> 53   AA  AC   CC
#> 54   GG  GG   GG
#> 55   GA  GA   GG
#> 56   AA  AG   AA
#> 57   NN  NN   GG
#> 58   CC  TC   CC
#> 59   AA  NN   NN
#> 60   AC  AA   AC
#> 61   TC  TT   NN
#> 62   GG  GG   GG
#> 63   AA  NN   AA
#> 64   TA  TA   TA
#> 65   CT  CT   TT
#> 66   CC  NN   AC
#> 67   GG  GG   GG
#> 68   AA  AA   AT
#> 69   AA  AA   AA
#> 70   NN  AC   AC
#> 71   AT  AA   AA
#> 72   GG  GA   GA
#> 73   TT  NN   TG
#> 74   GT  NN   GG
#> 75   CC  AC   CC
#> 76   CC  NN   TT
#> 77   AA  AT   AA
#> 78   GG  GG   GG
#> 79   TC  TT   TT
#> 80   TT  TT   TT
#> 81   AA  AT   AT
#> 82   TT  NN   CT
#> 83   TT  NN   TT
#> 84   CT  NN   CC
#> 85   GG  AG   AG
#> 86   AA  AG   AA
#> 87   CT  CC   CT
#> 88   AG  AA   AG
#> 89   GG  GG   CC
#> 90   CC  NN   CC
#> 91   CC  CT   CT
#> 92   NN  AG   GG
#> 93   TC  TC   TC
#> 94   CC  CC   CC
#> 95   NN  NN   GG
#> 96   CC  CA   CC
#> 97   GG  GG   GG
#> 98   GT  GT   TT
#> 99   AA  AG   AG
#> 100  AT  AT   AT

# from an adegenet genlight object
num <- format_snps(stickSNPs, "sn", interpolate = FALSE)
genlight <- methods::as(t(num[,-c(1:2)]), "genlight")

## run the conversion, could be run without the snp or sample metadatas.
dat <- import.snpR.data(genlight)
#> 
  |                                                                            
  |=                                                                     |   1%
  |                                                                            
  |=                                                                     |   2%
  |                                                                            
  |==                                                                    |   3%
  |                                                                            
  |===                                                                   |   4%
  |                                                                            
  |====                                                                  |   5%
  |                                                                            
  |====                                                                  |   6%
  |                                                                            
  |=====                                                                 |   7%
  |                                                                            
  |======                                                                |   8%
  |                                                                            
  |======                                                                |   9%
  |                                                                            
  |=======                                                               |  10%
  |                                                                            
  |========                                                              |  11%
  |                                                                            
  |========                                                              |  12%
  |                                                                            
  |=========                                                             |  13%
  |                                                                            
  |==========                                                            |  14%
  |                                                                            
  |==========                                                            |  15%
  |                                                                            
  |===========                                                           |  16%
  |                                                                            
  |============                                                          |  17%
  |                                                                            
  |=============                                                         |  18%
  |                                                                            
  |=============                                                         |  19%
  |                                                                            
  |==============                                                        |  20%
  |                                                                            
  |===============                                                       |  21%
  |                                                                            
  |===============                                                       |  22%
  |                                                                            
  |================                                                      |  23%
  |                                                                            
  |=================                                                     |  24%
  |                                                                            
  |==================                                                    |  25%
  |                                                                            
  |==================                                                    |  26%
  |                                                                            
  |===================                                                   |  27%
  |                                                                            
  |====================                                                  |  28%
  |                                                                            
  |====================                                                  |  29%
  |                                                                            
  |=====================                                                 |  30%
  |                                                                            
  |======================                                                |  31%
  |                                                                            
  |======================                                                |  32%
  |                                                                            
  |=======================                                               |  33%
  |                                                                            
  |========================                                              |  34%
  |                                                                            
  |========================                                              |  35%
  |                                                                            
  |=========================                                             |  36%
  |                                                                            
  |==========================                                            |  37%
  |                                                                            
  |===========================                                           |  38%
  |                                                                            
  |===========================                                           |  39%
  |                                                                            
  |============================                                          |  40%
  |                                                                            
  |=============================                                         |  41%
  |                                                                            
  |=============================                                         |  42%
  |                                                                            
  |==============================                                        |  43%
  |                                                                            
  |===============================                                       |  44%
  |                                                                            
  |================================                                      |  45%
  |                                                                            
  |================================                                      |  46%
  |                                                                            
  |=================================                                     |  47%
  |                                                                            
  |==================================                                    |  48%
  |                                                                            
  |==================================                                    |  49%
  |                                                                            
  |===================================                                   |  50%
  |                                                                            
  |====================================                                  |  51%
  |                                                                            
  |====================================                                  |  52%
  |                                                                            
  |=====================================                                 |  53%
  |                                                                            
  |======================================                                |  54%
  |                                                                            
  |======================================                                |  55%
  |                                                                            
  |=======================================                               |  56%
  |                                                                            
  |========================================                              |  57%
  |                                                                            
  |=========================================                             |  58%
  |                                                                            
  |=========================================                             |  59%
  |                                                                            
  |==========================================                            |  60%
  |                                                                            
  |===========================================                           |  61%
  |                                                                            
  |===========================================                           |  62%
  |                                                                            
  |============================================                          |  63%
  |                                                                            
  |=============================================                         |  64%
  |                                                                            
  |==============================================                        |  65%
  |                                                                            
  |==============================================                        |  66%
  |                                                                            
  |===============================================                       |  67%
  |                                                                            
  |================================================                      |  68%
  |                                                                            
  |================================================                      |  69%
  |                                                                            
  |=================================================                     |  70%
  |                                                                            
  |==================================================                    |  71%
  |                                                                            
  |==================================================                    |  72%
  |                                                                            
  |===================================================                   |  73%
  |                                                                            
  |====================================================                  |  74%
  |                                                                            
  |====================================================                  |  75%
  |                                                                            
  |=====================================================                 |  76%
  |                                                                            
  |======================================================                |  77%
  |                                                                            
  |=======================================================               |  78%
  |                                                                            
  |=======================================================               |  79%
  |                                                                            
  |========================================================              |  80%
  |                                                                            
  |=========================================================             |  81%
  |                                                                            
  |=========================================================             |  82%
  |                                                                            
  |==========================================================            |  83%
  |                                                                            
  |===========================================================           |  84%
  |                                                                            
  |============================================================          |  85%
  |                                                                            
  |============================================================          |  86%
  |                                                                            
  |=============================================================         |  87%
  |                                                                            
  |==============================================================        |  88%
  |                                                                            
  |==============================================================        |  89%
  |                                                                            
  |===============================================================       |  90%
  |                                                                            
  |================================================================      |  91%
  |                                                                            
  |================================================================      |  92%
  |                                                                            
  |=================================================================     |  93%
  |                                                                            
  |==================================================================    |  94%
  |                                                                            
  |==================================================================    |  95%
  |                                                                            
  |===================================================================   |  96%
  |                                                                            
  |====================================================================  |  97%
  |                                                                            
  |===================================================================== |  98%
  |                                                                            
  |===================================================================== |  99%
  |                                                                            
  |======================================================================| 100%
#> Citations for methods used thus far: 
#> ==============================================
#> Statistic:  genlight import 
#> Citation:  Gruber B, Unmack PJ, Berry OF, Georges A (2018). “dartr: An r packageto facilitate analysis of SNP data generated from reducedrepresentation genome sequencing.” _Molecular Ecology Resources_,*18*(3), 691-699. ISSN 1755-098X. doi:10.1111/1755-0998.12745<https://doi.org/10.1111/1755-0998.12745>.<https://doi.org/10.1111/1755-0998.12745>. 
#> Bibtex key:  Gruber2018 
#> Details:  genlight import managed by dartR's gl2gi function 
#> ==============================================
#> 
#> Warning: Some sample metadata columns contain unacceptable special characters. Unaccepted characters: '.', '~', or any whitespace.
#> These can cause unexpected behaviour if the subect columns are used as facets.
#> Issues:
#> 
#> Facet: sampID    levels: ASP.r2.141, ASP.r2.144, ASP.r2.146, ASP.r2.147, ASP.r2.196, ASP.r2.216, ASP.r2.221, ASP.r2.227, ASP.r2.230, ASP.r2.232, ASP.r2.235, CLF.r2.040, CLF.r2.048, CLF.r2.061, CLF.r2.064, CLF.r2.068, CLF.r2.069, OPL.r1.078, OPL.r1.080, OPL.r1.083, OPL.r1.086, OPL.r1.087, OPL.r1.088, OPL.r1.109, OPL.r1.110, OPL.r1.116, OPL.r1.125, OPL.r1.126, OPL.r1.128, OPL.r1.130, OPL.r1.139, OPL.r1.146, OPL.r1.147, OPL.r1.149, OPL.r2.007, OPL.r2.012, OPL.r2.022, OPL.r2.023, OPL.r2.028, OPL.r2.031, OPL.r2.033, OPL.r2.079, OPL.r2.080, OPL.r2.094, PAL.r2.203, PAL.r2.237, PAL.r2.238, PAL.r2.239, PAL.r2.240, PAL.r2.241, PAL.r2.242, PAL.r2.249, PAL.r2.253, PAL.r2.275, PAL.r2.276, PAL.r2.278, PAL.r2.281, PAL.r2.285, PAL.r2.287, PAL.r2.290, PAL.r2.291, SMR.r1.003, SMR.r1.004, SMR.r1.008, SMR.r1.012, SMR.r1.013, SMR.r1.023, SMR.r1.028, SMR.r1.034, SMR.r1.049, SMR.r1.052, SMR.r1.056, SMR.r1.072, SMR.r1.074, SMR.r2.001, SMR.r2.077, SMR.r2.086, SMR.r2.091, SMR.r2.097, SMR.r2.105, SMR.r2.139, SMR.r2.175, SMR.r2.256, SMR.r2.263, SMR.r2.296, UPD.r2.120, UPD.r2.121, UPD.r2.122, UPD.r2.123, UPD.r2.124, UPD.r2.129, UPD.r2.130, UPD.r2.131, UPD.r2.133, UPD.r2.135, UPD.r2.160, UPD.r2.162, UPD.r2.174, UPD.r2.189, UPD.r2.201

if (FALSE) { # \dontrun{
## not run:
# from a file:
# note that the drop argument is passed to data.table::fread!
dat <- import.snpR.data(system.file("extdata", "stick_NN_input.txt", 
                                    package = "snpR"), drop = 1:2) 
# if wanted, snp and sample metadata could be provided as usual.

## not run:
# from plink:
# make plink data
format_snps(stickSNPs, "plink", outfile = "plink_test", chr = "chr")

# read plink
dat <- import.snpR.data("plink_test.bed")
} # }
```
