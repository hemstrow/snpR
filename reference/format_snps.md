# Re-format SNP data.

`format_snps` re-formats SNP data into a range of different possible
formats for use in snpR functions and elsewhere.

## Usage

``` r
format_snps(
  x,
  output = "snpRdata",
  facets = NULL,
  n_samp = NA,
  interpolate = "bernoulli",
  outfile = FALSE,
  ped = NULL,
  input_format = NULL,
  input_meta_columns = NULL,
  input_mDat = NULL,
  sample.meta = NULL,
  snp.meta = NULL,
  ncp = 2,
  ncp.max = 5,
  chr = "chr",
  position = "position",
  phenotype = "phenotype",
  plink_recode_numeric = FALSE,
  sn_remove_empty = TRUE,
  verbose = FALSE
)
```

## Arguments

- x:

  snpRdata object or data.frame. Input data, in any of the above listed
  input formats.

- output:

  Character, default "snpRdata". The desired output format, see
  description for details.

- facets:

  Character or NULL, default NULL. Facets over which to break up data
  for some output formats, such as within one file for `genepop` or
  across multiple files for `vcf`, following the format described in
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).

- n_samp:

  Integer or numeric vector, default NA. For structure or RAFM outputs.
  How many random loci should be selected? Can either be an integer or a
  numeric vector of loci to use.

- interpolate:

  Character or FALSE, default "bernoulli". If transforming to "sn" or
  "pa" format, notes the interpolation method to be used to fill missing
  data. Options are "bernoulli", "af", "iPCA", or FALSE. See details.

- outfile:

  character vector, default FALSE. If a file path is provided, a copy of
  the output will be saved to that location. For some output styles,
  such as genepop, additional lines will be added to the output to allow
  them to be immediately run on commonly used programs.

- ped:

  data.frame default NULL. Optional argument for the "plink" output
  format. A six column data frame containing Family ID, Individual ID,
  Paternal ID, Maternal ID, Sex, and Phenotype and one row per sample.
  If provided, outputs will contain information contained in ped. See
  plink documentation for more details.

- input_format:

  Character, default NULL. Format of x, by default a snpRdata object.
  See description for details.

- input_meta_columns:

  Numeric, default NULL. If x is not a snpRdata object, optionally
  specifies the number of metadata columns preceding genotypes in x. See
  details for more information.

- input_mDat:

  Character, default "NN". If x is not a snpRdata object, the coding for
  missing *genotypes* in x (typically "NN" or "0000").

- sample.meta:

  data.frame, default NULL. If x is not a snpRdata object, optionally
  specifies a data.frame containing meta data for each sample. See
  details for more information.

- snp.meta:

  data.frame, default NULL. If x is not a snpRdata object, optionally
  specifies a data.frame containing meta data for each SNP. See details
  for more information.

- ncp:

  numeric or NULL, default 2. Number of components to consider for iPCA
  sn format interpolations of missing data. If null, the optimum number
  will be estimated, with the maximum specified by ncp.max. This can be
  very slow.

- ncp.max:

  numeric, default 5. Maximum number of components to check for when
  determining the optimum number of components to use when interpolating
  sn data using the iPCA approach.

- chr:

  character, default "chr". Name of column containing chromosome
  information, for VCF or plink! output.

- position:

  character, default "position". Name of column containing position
  information, for VCF output.

- phenotype:

  character, default "phenotype". Optional name of column containing
  phenotype information, for plink! output.

- plink_recode_numeric:

  Logical, default FALSE. If FALSE, all chrs/scaffs will be renamed to
  numbers. This may be useful in some cases. If this is FALSE,
  chromosome names will be checked for leading numbers and replaced with
  a corresponding letter (0 becomes A, 1 becomes B, and so on).

- sn_remove_empty:

  Logical, default TRUE. If TRUE, loci will be removed upon conversion
  to sn format if they are called in no individuals.

- verbose:

  Logical, default FALSE. If TRUE, some progress updates will be
  reported.

## Value

A data.frame or snpRdata object with data in the correct format. May
also write a file to the specified path.

## Details

While this function can accept a few non-snpRdata input formats, it will
reformat to a snpRdata object internally. As such, it takes a facets
argument that works identically to elsewhere in the package, as
described in
[`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).
This argument is only used for output formats where facets are
important, such as the genepop format.

While this function can be used as an alternative to
[`import.snpR.data`](https://hemstrow.github.io/snpR/reference/import.snpR.data.md)
when the output argument is set to "snpRdata", this is more complicated
and not recommended. Instead, it is simpler
[`import.snpR.data`](https://hemstrow.github.io/snpR/reference/import.snpR.data.md)
or one of the wrappers in
[`snpR_import_wrappers`](https://hemstrow.github.io/snpR/reference/snpR_import_wrappers.md).
The option is kept for backwards compatibility and internal use.

If non-snpRdata is supplied, SNP and sample metadata may be provided.
SNP metadata may either be provided in the first few columns of x, the
number of which is designated by the input_meta_columns argument, or in
a data.frame given as via the snp.meta argument. Sample metadata may be
provided in a data.frame via the sample.meta argument.

Output format options:

- ac: allele count format, allele counts tabulated for all samples or
  within populations.

- genepop: genepop format, genotypes stored as four numeric characters
  (e.g. "0101", "0204"), transposed, and formatted for genepop. Rownames
  are individual IDs in genepop format, colnames are SNP ids, matching
  the first metadata column in input.

- structure: STRUCTURE format, two lines per individual: allele calls
  stored as single character numeric (e.g. "1", "2"). Allele calls per
  individual stored on two subsequent lines.

- 0000: numeric genotype tab format, genotypes stored as four numeric
  characters (e.g. "0101", "0204").

- hapmap: Migrate-n hapmap, allele counts tabulated within populations,
  in migrate-n hapmap format. Since this migrate-n implementation is
  iffy, this probably shouldn't be used much.

- NN: character genotype tab format, genotypes stored as actual base
  calls (e.g. "AA", "CT").

- pa: allele presence/absence format, presence or absence of each
  possible allele at each possible genotype noted. Interpolation
  possible, with missing data substituted with allele frequency in all
  samples or each population.

- rafm: RAFM format, two allele calls at each locus stored in subsequent
  columns, e.g. locus1.1 locus1.2.

- faststructure: fastSTRUCTURE format, identical to STRUCTURE format
  save with the addition of filler columns proceeding data such that
  exactly 6 columns proceed data. These columns can be filled with
  metadata if desired.

- dadi: dadi format SNP data format, requires two columns named "ref"
  and "anc" with the flanking bases around the SNP, e.g. "ACT" where the
  middle location is the A/C snp.

- plink: PLINK! binary input format, requires columns named "position"
  and one matching the name designated with the 'chr' argument, and may
  contain a column named "cM", "cm", or "morgans", containing linkage
  group/chr, snp ID, position in bp, and distance in cM in order to
  create .bim extended map file.

- sn: Single character numeric format. Each genotype will be listed as
  0, 1, or 2, corresponding to 0, 1, or 2 minor alleles. Can be
  interpolated to remove missing data with the 'interpolate' argument.

- sequoia: sequoia format. Each genotype is converted to 0/1/2/ or -9
  (for missing values). Requires columns ID, Sex, BirthYear (or instead
  of BirthYear - BYmin and BYmax), optional column Yearlast in sample
  metadata for running Sequoia. For more information see sequoia
  documentation.

- fasta: fasta sequence format.

- vcf: Variant Call Format, a standard format for SNPs and other genomic
  variants. Genotypes are coded as 0/0, 0/1, 1/1, or ./. (for missing
  values), with a healthy serving of additional metadata but very little
  sample metadata.

- genalex: GenAlEx format. If an outfile is requested, the data will be
  sorted according to any provided facets and written as an '.xlsx'
  object.

- snpRdata: a snpRdata object.

Note that for the "sn" format, the data can be interpolated to fill
missing data points, which is useful for PCA, genomic prediction, tSNE,
and other methods. To do so, specify interpolate = "af" to insert the
expected number of minor alleles given SNP allele frequency or
"bernoulli" to do binomial draws to determine the number of minor
alleles at each missing data point, where the probability of drawing a
minor allele is equal to the minor allele frequency. The expected number
of minor alleles based on the later method is equal to the interpolated
value from the former, but the later allows for multiple runs to
determine the impact of stochastic draws and is generally preferred and
required for some downstream analysis. It is therefore the default. As a
slower but more accurate alternative to "af" interpolation, "iPCA" may
be selected. This an iterative PCA approach to interpolate based on
SNP/SNP covariance via
[`imputePCA`](https://rdrr.io/pkg/missMDA/man/imputePCA.html). If the
ncp argument is not defined, the number of components used for
interpolation will be estimated using
[`estim_ncpPCA`](https://rdrr.io/pkg/missMDA/man/estim_ncpPCA.html). In
this case, this method is much slower than the other methods, especially
for large datasets. Setting an ncp of 2-5 generally results in
reasonable interpolations without the time constraint.

Note also that for the plink format, a .bed binary file can be
generated. If the "plink" option is selected and an outfile is
designated, R will generate a ".sh" shell file with the same name given
in the outfile argument. Running this file will create a plink.bed file.

Input formats:

- NULL or snpRdata: snpRdata object, the default.

- NN: SNP genotypes stored as actual base calls (e.g. "AA", "CT").

- 0000: SNP genotypes stored as four numeric characters (e.g. "0101",
  "0204").

- snp_tab: SNP genotypes stored with genotypes in each cell, but only a
  single nucleotide noted if homozygote and two nucleotides seperated by
  a space if heterozygote (e.g. "T", "T G").

- sn: SNP genotypes stored with genotypes in each cell as 0 (homozyogous
  allele 1), 1 (heterozygous), or 2 (homozyogus allele 2).

## Author

William Hemstrom

Melissa Jones

## Examples
