#'snpRdata Import Wrappers
#'
#'These functions wrap \code{\link{import.snpR.data}} to import data into the
#'snpRdata format from a range of file or object sources.
#'
#'These functions are all wrappers for \code{\link{import.snpR.data}}, and all
#'are technically cross-compatible save read_ms: each other function can
#'actually be called with any of the supported formats (read_vcf can be handed a
#'genlight object without failure). These are supported as separate functions
#'for code readability and for ease of discovery.
#'
#'See \code{\link{import.snpR.data}} for more detail.
#'
#'@section File import: Supports automatic import of several types of files.
#'  Options:
#'
#'  \itemize{\item{.vcf or .vcf.gz: } Variant Call Format (vcf) files, supported
#'  via \code{\link[vcfR]{vcfR}}. If not otherwise provided, snp metadata is
#'  taken from the fixed fields in the VCF and sample metadata from the sample
#'  IDs. Note that this only imports SNPs with called genotypes! \item{.ms: }
#'  Files in the ms format, as provided by many commonly used simulation tools.
#'  \item{NN: } SNP genotypes stored as actual base calls (e.g. "AA", "CT").
#'  \item{0000: }SNP genotypes stored as four numeric characters (e.g. "0101",
#'  "0204"). \item{snp_tab: }SNP genotypes stored with genotypes in each cell,
#'  but only a single nucleotide noted if homozygote and two nucleotides
#'  separated by a space if heterozygote (e.g. "T", "T G"). \item{sn: }SNP
#'  genotypes stored with genotypes in each cell as 0 (homozygous allele 1), 1
#'  (heterozygous), or 2 (homozyogus allele 2).\item{genepop: } genepop file
#'  format, with genotypes stored as either 4 or 6 numeric characters. Works
#'  only with bi-allelic data. Genotypes will be converted (internally) to NN:
#'  the first allele (numerically) will be coded as A, the second as C.
#'  \item{FSTAT: } FSTAT file
#'  format, with genotypes stored as either 4 or 6 numeric characters. Works
#'  only with bi-allelic data. Genotypes will be converted (internally) to NN:
#'  the first allele (numerically) will be coded as A, the second as C.
#'  \item{plink: } plink .bed, .fam, and .bim files, via
#'  \code{\link[genio]{read_plink}}. If any of these file types is provided,
#'  snpR via \code{\link[genio]{read_plink}} will look for the other file types
#'  automatically. Sample metadata should be contained in the .fam file and SNP
#'  metadata in the .bim file, so sample or snp meta data can be provided here.
#'  \item{structure: } STRUCTURE import file, with individuals in rows and loci
#'  in columns. Can be coded either with one row per individual and two columns
#'  per loci or two rows per individual and two columns per loci using the
#'  rows_per_individual argument. Genotypes can be pretty much anything,
#'  although missing genotypes must be coded as -9. Must have a .str extension 
#'  and be consistantly whitespace delimited.}
#'
#'  Sample and snp metadata can also be provided via file path, and will be read
#'  in using \code{\link[data.table]{fread}} \emph{with the default settings}
#'  using \code{\link{read_delimited_snps}}.
#'  If these settings are not correct, please read in the metadata manually and
#'  provide to import.snpR.data.
#'
#'@section Conversions from other S4 objects:
#'
#'Supports automatic conversions from some other popular S4 object types.
#'Options:
#'
#'\itemize{ \item{genind: } \code{\link[adegenet]{genind}} objects from
#'adegenet. Note, no need to import genepop objects, the equivalent statistics
#'are calculated automatically when functions called with facets. Sample and SNP
#'IDs as well as, when possible, pop IDs will be taken from the genind object.
#'This data will be added too but will not replace data provided to the SNP or
#'sample.meta arguments. Note that only \emph{SNP} data is currently allowed,
#'data with more than two alleles for loci will return an error. \item{genlight:
#'} \code{\link[adegenet]{genlight}} objects from adegenet. Sample and SNP IDs,
#'SNP positions, SNP chromosomes, and pop IDs will be taken from the genlight
#'object if possible. This data will be added too but will not replace data
#'provided to the SNP or sample.meta arguments. \item{vcfR: }
#'\code{\link[vcfR]{vcfR}} objects from vcfR. If not provided, snp metadata is
#'taken from the fixed fields in the VCF and sample metadata from the sample
#'IDs. Note that this only imports SNPs with called genotypes!}
#'
#'@param file character, path to a file containing genotype data to import.
#'@param genlight genlight object to convert, see
#'  \code{\link[adegenet]{genlight}}.
#'@param genind genind object to convert, see \code{\link[adegenet]{genind}}.
#'@param vcfR vcfR object to convert, see \code{\link[vcfR]{vcfR}}.
#'@param snp.meta data.frame or character, default NULL. Metadata for each SNP,
#'  must have a number of rows equal to the number of SNPs in the dataset. If
#'  NULL, a single "snpID" column will be added. If a character, the path to a
#'  file containing SNP metadata, one row per SNP, with named columns.
#'@param sample.meta data.frame, default NULL. Metadata for each individual
#'  sample, must have a number of rows equal to the number of samples in the
#'  data set. If NULL, a single "sampID" column will be added. If a character,
#'  the path to a file containing sample metadata, one row per sample, with
#'  named columns.
#'@param chr.length numeric, Specifies chromosome lengths. Note that a single
#'  value assumes that each chromosome is of equal length whereas a vector of
#'  values gives the length for each chromosome in order.
#'@param mDat character, defaults "0000", "NN", or "-9" depending on method. 
#'  Note, if the default is set but the
#'  data has genotypes stored in 6 characters, mDat will be set to "000000".
#'@param header_cols numeric, default 0. The number of snp metadata columns prior
#'  to snp genotypes when importing delimited snps.
#'@param rows_per_individual numeric (1 or 2), default 2. Number of rows used
#'  for each individual.
#'@param marker_names logical, default FALSE. If TRUE, assumes that a
#'  header row of marker is present.
#'@param fix_overlaps Logical, default TRUE. If TRUE, overlapping positions will
#'  be checked and fixed during 'ms' file import.
#'  
#'@author William Hemstrom
#'@author Brent Gruber (genlight conversion re-distributed here) 
#'
#'@aliases read_vcf read_ms read_delimited_snps read_genepop read_FSTAT convert_genlight convert_genind
#'@name snpR_import_wrappers
NULL

#' @export
#' @describeIn snpR_import_wrappers Import .vcf or .vcf.gz files.
read_vcf <- function(file, snp.meta = NULL, sample.meta = NULL){
  if(!grepl("\\.vcf$", file) & !grepl("\\.vcf\\.gz$", file)){
    stop("File extension is not .vcf or .vcf.gz. Please check that the correct file has been entered and rename if needed.\n")
  }
  return(import.snpR.data(file, snp.meta, sample.meta))
}

#' @export
#' @describeIn snpR_import_wrappers Import .ms files.
read_ms <- function(file, sample.meta = NULL, chr.length, fix_overlaps = TRUE){
  if(!grepl("\\.ms$", file)){
    stop("File extension is not .ms. Please check that the correct file has been entered and rename if needed.\n")
  }
  if(!is.numeric(chr.length)){
    stop("chr.length must be a numeric vector of either length 1 or equal to the number of chromosomes.\n")
  }
  return(import.snpR.data(file, sample.meta = sample.meta, chr.length = chr.length, fix_overlaps = fix_overlaps))
}

#' @export
#' @describeIn snpR_import_wrappers Import tab delimited data where genotypes
#'   are stored as: NN, 0000, or snp_tab format.
read_delimited_snps <- function(file, snp.meta = NULL, sample.meta = NULL, mDat = "NN", header_cols = 0){
  return(import.snpR.data(file, snp.meta, sample.meta, mDat = mDat, header_cols = header_cols))
}

#' @export
#' @describeIn snpR_import_wrappers Import genepop formatted data. 
read_genepop <- function(file, snp.meta = NULL, sample.meta = NULL, mDat = "0000"){
  if(!grepl("\\.genepop$", file) & !grepl("\\.gen$", file)){
    stop("File extension is not .genepop or .gen. Please check that the correct file has been entered and rename if needed.\n")
  }
  return(import.snpR.data(file, snp.meta, sample.meta, mDat))
}

#' @export
#' @describeIn snpR_import_wrappers Import FSTAT formatted data. 
read_FSTAT <- function(file, snp.meta = NULL, sample.meta = NULL, mDat = "0000"){
  if(!grepl("\\.fstat$", file)){
    stop("File extension is not .fstat. Please check that the correct file has been entered and rename if needed.\n")
  }
  return(import.snpR.data(file, snp.meta, sample.meta, mDat))
}

#' @export
#' @describeIn snpR_import_wrappers Import plink bed, bim, and fam data.
read_plink <- function(file){
  if(grepl("\\.bim$", file) | grepl("\\.fam$", file) | grepl("\\.bed$", file)){
    .check.installed("tools")
    return(.process_plink(tools::file_path_sans_ext(file)))
  }
  else{
    return(.process_plink(file))
  }
}

#' @export
#' @describeIn snpR_import_wrappers Import STRUCTURE data files.
read_structure <- function(file, snp.meta = NULL, sample.meta = NULL, 
                           rows_per_individual = 2, marker_names = FALSE,
                           header_cols = 0, mDat = -9){
  if(!grepl("\\.str$", file)){
    stop("File extension is not .str. Please check that the correct file has been entered and rename if needed.\n")
  }
  
  return(import.snpR.data(genotypes = file, rows_per_individual = rows_per_individual, 
                            marker_names = marker_names, 
                            header_cols = header_cols, snp.meta = snp.meta, sample.meta = sample.meta, mDat = mDat))
}

#' @export
#' @describeIn snpR_import_wrappers Convert adegenet genlight objects.
convert_genlight <- function(genlight, snp.meta = NULL, sample.meta = NULL){
  return(import.snpR.data(genlight, snp.meta, sample.meta))
}

#' @export
#' @describeIn snpR_import_wrappers Convert adegenet genind objects
convert_genind <- function(genind, snp.meta = NULL, sample.meta = NULL){
  return(import.snpR.data(genind, snp.meta, sample.meta))
}

#' @export
#' @describeIn snpR_import_wrappers Convert adegenet vcfR objects
convert_vcfR <- function(vcfR, snp.meta = NULL, sample.meta = NULL){
  return(import.snpR.data(vcfR, snp.meta, sample.meta))
}

