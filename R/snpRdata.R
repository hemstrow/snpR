#'@import data.table 
#'@importFrom foreach %dopar%
#'@import Matrix
NULL


.check.snpRdata <- function(object){
  errors <- character()

  # check that there are the correct number of snp and sample metadata rows
  if(nrow(object@snp.meta) != nrow(genotypes(object))){
    msg <- paste0("Number of snps in metadata (", nrow(object@snp.meta),
                  ") is not equal to the number of snps in data (", nrow(genotypes(object)), ").")
    errors <- c(errors, msg)
  }
  if(nrow(object@sample.meta) != ncol(genotypes(object))){
    msg <- paste0("Number of samples in metadata (", nrow(object@sample.meta),
                  ") is not equal to the number of samples in data (",  ncol(genotypes(object)), ").")
    errors <- c(errors, msg)
  }

  all.colnames <- c(colnames(object@sample.meta), colnames(object@snp.meta))
  .IDs <- which(all.colnames %in% c(".snp.id", ".sample.id"))
  if(length(.IDs) > 0){
    all.colnames <- all.colnames[-.IDs]
  }

  # check for duplicated or poorly named columns
  dups <- duplicated(all.colnames) | duplicated(all.colnames, fromLast = TRUE)
  bad.chars <- c("\\.", "\\~", " ")
  bad.chars <- paste(bad.chars, collapse = "|")
  p.char.matches <- grepl(bad.chars, all.colnames)

  bad.colnames <- dups | p.char.matches | all.colnames == "all"

  if(any(bad.colnames)){
    bad.colnames <- which(bad.colnames)
    msg <- paste0("Some unacceptable column names in snp or sample metadata. Column names must be unique and not contain either '.', '~', or whitespace, and cannot be 'all'. Bad column names: ",
                   paste0(all.colnames[bad.colnames], collapse = ", "), ".\n")

    errors <- c(errors, msg)
  }
  
  restricted.names <- c("facet", "subfacet", "facets", "subfacets", "snp.facet", "snp.subfacet", "snp.subfacets", "snp.facets")
  bad.colnames.facet <- which(all.colnames %in% restricted.names)
  if(length(bad.colnames.facet) > 0){
    msg <- paste0("Some restricted column names detected in metadata. Restricted names: ", paste0(restricted.names, collapse = ", "), ".\n")
    errors <- c(errors, msg)
  }

  # check that the facet meta and snp meta column names match
  if(!identical(colnames(object@snp.meta), colnames(object@facet.meta)[-c(1:3)])){
    errors <- c(errors, "Column names in the snp.meta and stored metadata for facets do not match. This is likely due to later change of the snp metadata which added or changed column names but didn't change facet.meta column names. Please re-initialize the snpRdata object using the metadata you require. This error should typically not appear.\n")
  }

  warns <- character()
  # check for bad entries in character columns (periods are fine in numeric!)
  chr.cols.samp <- unlist(lapply(object@sample.meta, class))
  chr.cols.samp <- chr.cols.samp[which(chr.cols.samp == "character")]
  chr.cols.snp <- unlist(lapply(object@snp.meta, class))
  chr.cols.snp <- chr.cols.snp[which(chr.cols.snp == "character")]

  if(length(chr.cols.samp) > 0){
    l1 <- lapply(object@sample.meta[,names(chr.cols.samp), drop = FALSE], grepl, pattern = bad.chars)
    l1s <- unlist(lapply(l1, sum))
    if(any(l1s > 0)){
      pst_msg <- paste0("Some sample metadata columns contain unacceptable special characters. Unaccepted characters: '.', '~', or any whitespace.\nThese can cause unexpected behaviour if the subect columns are used as facets.\nIssues:\n")
      for(q in 1:length(l1s)){
        if(l1s[q] > 0){
          pst_msg <- paste0(pst_msg, "\nFacet: ", names(l1s)[q], "\tlevels: ", 
                            paste0(unique(object@sample.meta[which(l1[[q]]),names(chr.cols.samp),drop = FALSE][,q]), collapse = ", "))
        }
      }
      warns <- c(warns, pst_msg)
    }
  }
  if(length(chr.cols.snp) > 0){
    l1 <- lapply(object@snp.meta[,names(chr.cols.snp), drop = FALSE], grepl, pattern = bad.chars)
    l1s <- unlist(lapply(l1, sum))
    if(any(l1s > 0)){
      pst_msg <- paste0("Some snp metadata columns contain unacceptable special characters. Unaccepted characters: '.', '~', or any whitespace.\nThese can cause unexpected behaviour if the subect columns are used as facets.\nIssues:\n")
      for(q in 1:length(l1s)){
        if(l1s[q] > 0){
          pst_msg <- paste0(pst_msg, "\nFacet: ", names(l1s)[q], "\tlevels: ", 
                            paste0(unique(object@snp.meta[which(l1[[q]]),names(chr.cols.snp), drop = FALSE][,q]), collapse = ", "))
        }
      }
      warns <- c(warns, pst_msg)
    }
  }
  
  
  # warn if anything repeated across sample level factors
  uniques <- lapply(object@sample.meta[,-which(names(object@sample.meta) == ".sample.id"), drop = FALSE], unique)
  uniques <- c(uniques, lapply(object@snp.meta[,-which(names(object@snp.meta) == ".snp.id"), drop = FALSE], unique))
  un <- rep(names(uniques), lengths(uniques))
  uniques <- unlist(uniques, use.names = FALSE)
  names(uniques) <- un
  rm(un)
  if(any(duplicated(uniques))){
    
    dups <- sort(uniques[which(duplicated(uniques) | duplicated(uniques, fromLast = TRUE))])
    msg <- unique(dups)
    pst_msg <- "Some levels are duplicated across multiple sample meta facets.\nThis will cause issues if those sample facets are run during analysis.\nIssues:\n"
    for(q in 1:length(msg)){
      pst_msg <- paste0(pst_msg, "\nLevel: ", msg[q], "\tin facets: ", paste0(names(dups)[which(dups == msg[q])], collapse = ", "))
    }
    warns <- c(warns, pst_msg)
  }
  
  
  # check .snp.id and .sample.id columns
  if(colnames(object@snp.meta)[ncol(object@snp.meta)] != ".snp.id"){
    errors <- c(errors, ".snp.id not last column in snp.meta. This is a developer error--please report on the github issues page with a reproducable example.\n")
  }
  if(colnames(object@sample.meta)[ncol(object@sample.meta)] != ".sample.id"){
    errors <- c(errors, ".sample.id not last column in sample.meta. This is a developer error--please report on the github issues page with a reproducable example.\n")
  }
  
  if(length(warns) > 0){warning(paste0(warns, collapse = paste0("\n\n", .console_hline(), "\n")))}

  if(length(errors) == 0){return(TRUE)}
  else{return(errors)}
}



#'Storage class for snpR data and calculated statistics.
#'
#'The snpRdata class stores both raw genotype data, sample and locus specific
#'metadata, useful data summaries, repeatedly internally used tables, calculated
#'summary statistics, and smoothed statistic data. Used by most snpR functions.
#'
#'The snpRdata class is built to contain SNP genotype data for use by functions
#'in the snpR package. It also stores sample and locus specific metadata,
#'genomic summary information, and any results from most snpR functions. Results
#'can be gathered using \code{\link{get.snpR.stats}}. Genotypes are stored in
#'the "character" format, as output by format_snps().
#'
#'For more information, see \code{\link{import.snpR.data}}, the constructor
#'function for this object class.
#'
#'
#'@slot sample.meta data.frame containing sample metadata
#'@slot snp.meta data.frame containing snp metadata
#'@slot facet.meta data.frame containing snp metadata for the geno.tables (in
#'  the order of snps in that data)
#'@slot mDat character, missing data key
#'@slot snp.form numeric, number of characters per genotype (not really used)
#'@slot geno.tables list containing three matrices: gs (genotype counts), as
#'  (allele counts), and wm (missing counts)
#'@slot facets character, vector of tabulated facets
#'@slot facet.type character, types of each tabulated facet.
#'@slot stats data.frame, all calculated single-snp non-pairwise stats.
#'@slot window.stats data.frame, all calculated window non-pairwise stats.
#'@slot pairwise.stats data.frame, all calculated pairwise single snp stats.
#'@slot pairwise.window.stats data.frame, all calculated window pairwise stats.
#'@slot sample.stats data.frame, all calculated single sample stats.
#'@slot pop.stats data.frame, all calculated population level summary stats.
#'@slot pairwise.LD nested list of matrices, all calculated pairwise LD stats.
#'@slot window.bootstraps data.frame, all calculated bootstraps
#'@slot sn list, contains sn formatted data if calculated
#'@slot calced_stats list, contains information on what statistics have been
#'  calculated for which facets
#'@slot allele_frequency_matrices list of matrices, allele frequency matrices if
#'  calculated
#'@slot genetic_distances list of matrices, genetic distance matrices if
#'  calculated
#'@slot other list, contains other miscellaneous calculated statistics that do
#'  not fit cleanly elsewhere.
#'
#'@author William Hemstrom
#'
#'@importFrom methods new
#'  
snpRdata <- setClass(Class = 'snpRdata', slots = c(sample.meta = "data.frame",
                                       snp.meta = "data.frame",
                                       facet.meta = "data.frame",
                                       mDat = "character",
                                       ploidy = "numeric",
                                       bi_allelic = "logical",
                                       data.type = "character",
                                       snp.form = "numeric",
                                       geno.tables = "list",
                                       facets = "character",
                                       facet.type = "character",
                                       stats = "data.table",
                                       window.stats = "data.table",
                                       pairwise.stats = "data.table",
                                       pairwise.window.stats = "data.table",
                                       sample.stats = "data.table",
                                       pop.stats = "data.table",
                                       pairwise.LD = "list",
                                       window.bootstraps = "data.table",
                                       sn = "list",
                                       calced_stats = "list",
                                       filters = "data.frame",
                                       allele_frequency_matrices = "list",
                                       genetic_distances = "list",
                                       weighted.means = "data.table",
                                       other = "list",
                                       citations = "list"),
         contains = c(data = "data.frame"),
         prototype = prototype(sample.meta = data.frame(NULL),
                          snp.meta = data.frame(NULL),
                          facet.meta = data.frame(NULL),
                          mDat = character(0),
                          ploidy = numeric(0),
                          bi_allelic = logical(0),
                          data.type = character(0),
                          snp.form = numeric(0),
                          geno.tables = vector("list"),
                          facets = character(0),
                          facet.type = character(0),
                          stats = data.table::data.table(NULL),
                          window.stats = data.table::data.table(NULL),
                          pairwise.stats = data.table::data.table(NULL),
                          pairwise.window.stats = data.table::data.table(NULL),
                          sample.stats = data.table::data.table(NULL),
                          pop.stats = data.table::data.table(NULL),
                          pairwise.LD = vector("list"),
                          window.bootstraps = data.table::data.table(NULL),
                          sn = vector("list"),
                          calced_stats = vector("list"),
                          filters = data.frame(NULL),
                          allele_frequency_matrices = vector("list"),
                          genetic_distances = vector("list"),
                          weighted.means = data.table::data.table(NULL),
                          other = vector("list"),
                          citations = vector("list")),
         validity = .check.snpRdata)


#'Import genotype and metadata into a snpRdata object.
#'
#'\code{import.snpR.data} converts genotype and meta data to the snpRdata class,
#'which stores raw genotype data, sample and locus specific metadata, useful
#'data summaries, repeatedly internally used tables, calculated summary
#'statistics, and sliding-window statistic data.
#'
#'The snpRdata class is built to contain SNP genotype data for use by functions
#'in the snpR package. It inherits from the S3 class data.frame, in which the
#'genotypes are stored, and can be manipulated identically. It also stores
#'sample and locus specific metadata, genomic summary information, and results
#'from most snpR functions. Genotypes are stored in the "character" format, as
#'output by \code{\link{format_snps}}. Missing data is noted with "NN".
#'
#'Inputs can be provided either as pre-existing R objects (in several different
#'formats) or as paths to files. In both cases, \code{snpR} will attempt to
#'guess the data format from either the object classs, first genotype or the
#'file extension (as appropriate).
#'
#'@section File import: Supports automatic import of several types of files.
#'  Options:
#'
#'  \itemize{\item{.vcf or .vcf.gz: } Variant Call Format (vcf) files, supported
#'  via \code{\link[vcfR]{vcfR}}. If not otherwise provided, snp metadata is
#'  taken from the fixed fields in the VCF and sample metadata from the sample
#'  IDs. Note that this only imports SNPs with called genotypes! \item{.ms: }
#'  Files in the ms format, as provided by many commonly used simulation tools.
#'  \item{.txt, NN: } SNP genotypes stored as actual base calls (e.g. "AA",
#'  "CT"). \item{.txt, 0000: }SNP genotypes stored as four numeric characters
#'  (e.g. "0101", "0204"). \item{.txt, snp_tab: }SNP genotypes stored with
#'  genotypes in each cell, but only a single nucleotide noted if homozygote and
#'  two nucleotides separated by a space if heterozygote (e.g. "T", "T G").
#'  \item{.txt, sn: }SNP genotypes stored with genotypes in each cell as 0
#'  (homozygous allele 1), 1 (heterozygous), or 2 (homozyogus allele 2).
#'  \item{.genepop or .gen: } genepop file format, with genotypes stored as either 4 or
#'  6 numeric characters. Works only with bi-allelic data. Genotypes will be
#'  converted (internally) to NN: the first allele (numerically) will be coded
#'  as A, the second as C. \item{.fstat: } FSTAT file format, with genotypes
#'  stored as either 4 or 6 numeric characters. Works only with bi-allelic data.
#'  Genotypes will be converted (internally) to NN: the first allele
#'  (numerically) will be coded as A, the second as C. \item{.bed/.fam/.bim: }
#'  PLINK .bed, .fam, and .bim files, via \code{\link[genio]{read_plink}}. If
#'  any of these file types is provided, snpR (via
#'  \code{\link[genio]{read_plink}}) will look for the other file types
#'  automatically. Sample metadata should be contained in the .fam file and SNP
#'  metadata in the .bim file, so sample or snp meta data provided here will be
#'  ignored. \item{.str: } STRUCTURE import files in either 1 or 2 rows per
#'  individual as defined by the \code{rows_per_individual} argument.}
#'
#'  Additional arguments can be provided to import.snpR.data that will be passed
#'  to \code{\link[data.table]{fread}} when reading in genotype data.
#'
#'  Sample and snp metadata can also be provided via file path, and will be read
#'  in using \code{\link[data.table]{fread}} \emph{with the default settings}.
#'  If these settings are not correct, please read in the metadata manually and
#'  provide to import.snpR.data.
#'
#'@section Conversions from other S4 objects:
#'
#'  Supports automatic conversions from some other popular S4 object types.
#'  Options:
#'
#'  \itemize{ \item{genind: } \code{\link[adegenet]{genind}} objects from
#'  adegenet. Note, no need to import genepop objects, the equivalent statistics
#'  are calculated automatically when functions called with facets. Sample and
#'  SNP IDs as well as, when possible, pop IDs will be taken from the genind
#'  object. This data will be added too but will not replace data provided to
#'  the SNP or sample.meta arguments. Note that only \emph{SNP} data is
#'  currently allowed, data with more than two alleles for loci will return an
#'  error. \item{genlight: } \code{\link[adegenet]{genlight}} objects from
#'  adegenet. Sample and SNP IDs, SNP positions, SNP chromosomes, and pop IDs
#'  will be taken from the genlight object if possible. This data will be added
#'  too but will not replace data provided to the SNP or sample.meta arguments.
#'  \item{vcfR: } \code{\link[vcfR]{vcfR}} objects from vcfR. If not provided,
#'  snp metadata is taken from the fixed fields in the VCF and sample metadata
#'  from the sample IDs. Note that this only imports SNPs with called
#'  genotypes!}
#'
#'@section Slots:
#'
#'  Genotypes, metadata, and results are stored in slots and directly accessable
#'  with the 'at' symbol operator. Slots are as follows:
#'
#'  \itemize{ \item{sample.meta: } sample metadata (population, family,
#'  phenotype, etc.). \item{snp.meta: } SNP metadata (SNP ID, chromosome,
#'  linkage group, position, etc.). \item{facet.meta: } internal metadata used
#'  to track facets that have been previously applied to the dataset.
#'  \item{mDat: } missing data format. \item{snp.form: } number of characters
#'  per SNP. \item{genotables: } a list containing tabulated genotypes (gs),
#'  allele counts (as), and missing data (wm). facet.meta contains the
#'  corresponding metadata. \item{facets: }
#'  vector of the facets that have been added to the data. \item{facet.type: }
#'  classes of the added facets (snp, sample, complex, or .base). \item{stats: }
#'  data.frame containing all calculated non-pairwise single-snp statistics and
#'  metadata. \item{window.stats: } data.frame/table containing all non-pairwise
#'  statistics calculated for sliding windows. \item{pairwise.stats: }
#'  data.frame/table containing all pairwise (fst) single-snp statistics.
#'  \item{pairwise.window.stats: } data.frame/table containing all pairwise
#'  statistics calculated for sliding windows. \item{sample.stats: }
#'  data.frame/table containing statistics calculated for each individual
#'  sample. \item{pairwise.LD: } nested list containing linkage disequilibrium
#'  data (see \code{\link{calc_pairwise_ld}} for more information).
#'  \item{window.bootstraps: } data.frame/table containing all calculated
#'  bootstraps for sliding window statistics. \item{sn: } list containing "sn",
#'  sn formatted data, and "type" type of interpolation. \item{calced_stats: }
#'  Named list of named character vectors that tracks the calculated statistics
#'  for each facet (see \code{\link{calc_genetic_distances}} for more
#'  information). \item{genetic_distances: } nested list containing genetic
#'  distance data. \item{names: } column names for genotypes. \item{row.names: }
#'  row names for genotypes. \item{.Data: } list of vectors containing raw
#'  genotype data. \item{.S3Class: } notes the inherited S3 object class. }
#'  
#'  Note that most of these slots are used primarily internally.
#'
#'  All calculated data can be accessed using the \code{\link{get.snpR.stats}}
#'  function. See documentation.
#'
#'
#'@param genotypes data.frame, unique S4 from other packages, or filename. If a
#'  data.frame, raw genotypes in a two-character format ("GG", "GA", "CT",
#'  "NN"), where SNPs are in rows and individual samples are in columns.
#'  Otherwise, see documentation for allowed S4 objects and files.
#'@param snp.meta data.frame, default NULL. Metadata for each SNP, must have a
#'  number of rows equal to the number of SNPs in the dataset. If NULL, a single
#'  "snpID" column will be added.
#'@param sample.meta data.frame, default NULL. Metadata for each individual
#'  sample, must have a number of rows equal to the number of samples in the
#'  dataset. If NULL, a single "sampID" column will be added.
#'@param mDat character, default "NN", matching the encoding of missing
#'  \emph{genotypes} in the data provided to the genotypes argument.
#'@param chr.length numeric, default NULL. If a path to a .ms file is provided,
#'  specifies chromosome lengths. Note that a single value assumes that each
#'  chromosome is of equal length whereas a vector of values gives the length
#'  for each chromosome in order.
#'@param ... Additional arguments passed to \code{\link[data.table]{fread}} if a
#'  \emph{genotype} file name is passed that is not a vcf or ms file.
#'@param header_cols numeric, default 0. Number of header columns containing
#'  SNP metadata. Used if a tab delimited or STRUCTURE input file is
#'  provided.
#'@param rows_per_individual numeric (1 or 2), default 2. Number of rows used
#'  for each individual. For structure input files only.
#'@param marker_names logical, default FALSE. If TRUE, assumes that a
#'  header row of marker is present. For structure input files only.
#'@param fix_overlaps Logical, default TRUE. If TRUE, overlapping positions will
#'  be checked and fixed during 'ms' file import.
#'@param verbose Logical, default FALSE. If TRUE, will print a few status
#'  updates and checks.
#'@param .pass_filters Internal, probably not for user use. Used to pass 
#'  filtering history when sub-setting when this function is called internally.
#'@param .skip_filters Internal, probably not for user use. Used to skip 
#'  re-filtering during sub-setting when this function is called internally.
#'  
#'@examples
#' # import example data as a snpRdata object
#' # produces data identical to that contained in the stickSNPs example dataset.
#' genos <- stickRAW[,-c(1:2)]
#' snp_meta <- stickRAW[,1:2]
#' sample_meta <- data.frame(pop = substr(colnames(stickRAW)[-c(1:2)], 1, 3), 
#'                           fam = rep(c("A", "B", "C", "D"), 
#'                                     length = ncol(stickRAW) - 2), 
#'                           stringsAsFactors = FALSE)
#' import.snpR.data(genos, snp.meta = snp_meta, sample.meta = sample_meta, 
#'                  mDat = "NN")
#'
#' # from an adegenet genind object
#' ex.genind  <- adegenet::df2genind(t(stickRAW[,-c(1:2)]), 
#'                                   ncode = 1, NA.char = "N") # get genind data
#' # note, will add whatever metadata data is in the genind object to the 
#' # snpRdata object. 
#' # Could be run without the snp or sample metadatas.
#' import.snpR.data(ex.genind, snp_meta, sample_meta) 
#'
#' # from an adegenet genlight object
#' num <- format_snps(stickSNPs, "sn", interpolate = FALSE)
#' genlight <- methods::as(t(num[,-c(1:2)]), "genlight")
#' 
#' ## run the conversion, could be run without the snp or sample metadatas.
#' dat <- import.snpR.data(genlight)
#'
#' \dontrun{
#' ## not run:
#' # from a file:
#' # note that the drop argument is passed to data.table::fread!
#' dat <- import.snpR.data(system.file("extdata", "stick_NN_input.txt", 
#'                                     package = "snpR"), drop = 1:2) 
#' # if wanted, snp and sample metadata could be provided as usual.
#' 
#' ## not run:
#' # from plink:
#' # make plink data
#' format_snps(stickSNPs, "plink", outfile = "plink_test", chr = "chr")
#'
#' # read plink
#' dat <- import.snpR.data("plink_test.bed")
#' }
#'
#'@export
#'
#'@author William Hemstrom
import.snpR.data <- function(genotypes, snp.meta = NULL, sample.meta = NULL, mDat = "NN", chr.length = NULL,
                             ..., header_cols = 0, rows_per_individual = 2, marker_names = FALSE, fix_overlaps = TRUE,
                             verbose = FALSE,
                             .pass_filters = FALSE, .skip_filters = FALSE){
  position <- .snp.id <- .sample.id <- NULL

  #======special cases========
  # sample and snp metadata
  if(is.character(sample.meta) & length(sample.meta) == 1){
    if(file.exists(sample.meta)){
      sample.meta <- as.data.frame(data.table::fread(sample.meta))
    }
    else{
      stop("Cannot locate sample.meta file.\n")
    }
  }
  else if(!is.null(sample.meta)){
    sample.meta <- try(as.data.frame(sample.meta), silent = TRUE)
    if(methods::is(sample.meta, "try-error")){
      stop(paste0("Could not convert sample.meta to data.frame. Error: \n", sample.meta))
    }
  }
  
  
  if(is.character(snp.meta) & length(snp.meta) == 1){
    if(file.exists(snp.meta)){
      snp.meta <- as.data.frame(data.table::fread(snp.meta))
    }
    else{
      stop("Cannot locate snp.meta file.\n")
    }
  }
  else if(!is.null(snp.meta)){
    snp.meta <- try(as.data.frame(snp.meta), silent = TRUE)
    if(methods::is(snp.meta, "try-error")){
      stop(paste0("Could not convert snp.meta to data.frame. Error: \n", snp.meta))
    }
  }
  
  # genotypes
  if(methods::is(genotypes, "genind")){
    return(.genind.tosnpRdata(genotypes, snp.meta, sample.meta))
  }
  if(methods::is(genotypes, "genlight")){
    return(.genlight.to.snpRdata(genotypes, snp.meta, sample.meta))
  }
  if(methods::is(genotypes, "vcfR")){
    return(.process_vcf(genotypes, snp.meta, sample.meta))
  }
  if(is.matrix(genotypes)){
    genotypes <- as.data.frame(genotypes)
  }
  if(is.character(genotypes) & length(genotypes) == 1){
    if(file.exists(genotypes)){
      # check for ms or vcf, etc file
      if(grepl("\\.vcf$", genotypes) | grepl("\\.vcf\\.gz$", genotypes)){
        return(.process_vcf(genotypes, snp.meta, sample.meta))
      }
      else if(grepl("\\.ms$", genotypes)){
        if(!is.null(snp.meta)){
          snp.meta <- NULL
          warning("Any provided snp.meta will be discarded during ms import.\n")
        }
        return(.process_ms(genotypes, chr.length, sample.meta, fix_overlaps))
      }
      else if(grepl("\\.genepop$", genotypes) | grepl("\\.gen$", genotypes)){
        return(.process_genepop(genotypes, snp.meta, sample.meta, mDat))
      }
      else if(grepl("\\.fstat$", genotypes)){
        return(.process_FSTAT(genotypes, snp.meta, sample.meta, mDat))
      }
      else if(grepl("\\.bim$", genotypes) | grepl("\\.fam$", genotypes) | grepl("\\.bed$", genotypes)){
        .check.installed("tools")
        return(.process_plink(tools::file_path_sans_ext(genotypes)))
      }
      else if(grepl("\\.str$", genotypes)){
        return(.process_structure(genotypes, 
                                  rows_per_individual = rows_per_individual, 
                                  marker_names = marker_names, 
                                  header_cols = header_cols, 
                                  snp.meta = snp.meta, 
                                  sample.meta = sample.meta,
                                  mDat = mDat))
      }
      else{
        genotypes <- as.data.frame(data.table::fread(genotypes, ...))
        if(sum(is.na(genotypes[,ncol(genotypes)])) == nrow(genotypes)){
          genotypes <- genotypes[,-ncol(genotypes)]
        }
      }
    }
    else{
      stop("File not found. Fix path or import manually and provide to import.snpR.data.\n")
    }
  }
  
  
  #=================check input format for non-special case=============================
  # NN, no need to do anything, just read in and proceed as normal.
  if(header_cols > 0){
    header_cols <- 1:header_cols
    snp.meta <- genotypes[,header_cols]
    genotypes <- genotypes[,-header_cols]
    
  }
  if(genotypes[1,1] %in% 
     c(apply(expand.grid(c("A", "T", "C", "G"), c("A", "T", "C", "G")), 1, paste, collapse=""), mDat)){
    if(verbose){cat("Assuming data is in NN format.\n")}
  }
  
  # sn
  else if(genotypes[1,1] %in% c(0, 1, 2, mDat)){
    if(verbose){cat("Assuming single nucleotide format.\n")}
    return(format_snps(genotypes, input_format = "sn", input_mDat = mDat, sample.meta = sample.meta, snp.meta = snp.meta))
  }
  
  # 0000
  else if(genotypes[1,1] %in% 
          c(apply(expand.grid(c("01", "02", "03", "04"), c("01", "02", "03", "04")), 1, paste, collapse=""), mDat)){
    if(verbose){cat("Assuming 0000 format.\n")}
    return(format_snps(genotypes, input_format = "0000", input_mDat = mDat, sample.meta = sample.meta, snp.meta = snp.meta))
  }
  
  # SNP_tab
  else if(genotypes[1,1] %in% 
          c(apply(expand.grid(c("A", "T", "C", "G"), c("A", "T", "C", "G")), 1, paste, collapse=" "), mDat, c("A", "T", "C", "G"))){
    if(verbose){cat("Assuming snp_tab format.\n")}
    return(format_snps(genotypes, input_format = "snp_tab", input_mDat = mDat, sample.meta = sample.meta, snp.meta = snp.meta))
  }
  
  #couldn't find a supported format
  else{
    stop("Genoytpes are not in a recognized format. First genotype:", genotypes[1,1], 
         "\n. Do you have SNP meta-data in your genotype file or object? If so, you may either remove those columns or use the header_cols argument to use them as your snp.meta.\n")
  }
  
  
  
  #============sanity checks and prep========
  if(is.null(snp.meta)){
    snp.meta <- data.frame(snpID = paste0("snp", 1:nrow(genotypes)))
  }
  if(is.null(sample.meta)){
    sample.meta <- data.frame(sampID = paste0("samp", 1:ncol(genotypes)))
  }
  
  # check mdat and genotype format -- only do the first sample to save on processing speed
  gtl <- unique(as.numeric(nchar(as.matrix(genotypes)[,1])))
  if(length(gtl) > 1){
    stop("All genotypes must be equal in length, including missing data.\n")
  }
  if(nchar(mDat) != gtl){
    stop("mDat must be equal in length (number of characters) to the genotype format, as in 'NN' or '00' with 'AC' genotypes.\n")
  }
  md1 <- substr(mDat, 0, nchar(mDat)/2)
  md2 <- substr(mDat, (nchar(mDat)/2) + 1, nchar(mDat))
  if(md1 != md2){
    stop("mDat must be symmetrical, as in 'NN' or '0000', NOT '01' or 'NGT' or 'NX'.\n")
  }
  
  # prepare things for addition to data
  if(any(is.na(genotypes))){
    stop("NA found in input genotypes. Often, this is in the last row or column.\n")
  }
  
  if(nrow(snp.meta) != nrow(genotypes)){
    stop(paste0("Number of rows in snp.meta (", nrow(snp.meta), ") not equal to number of SNPs in genotypes (", nrow(genotypes), "). Do you need to transpose your genotypes?\n"))
  }
  if(nrow(sample.meta) != ncol(genotypes)){
    stop(paste0("Number of rows in sample.meta (", nrow(sample.meta), ") not equal to number of samples in genotypes (", ncol(genotypes), "). Do you need to transpose your genotypes?\n"))
  }
  
  if(any(colnames(snp.meta) == "position")){
    snp.meta$position <- as.numeric(as.character(snp.meta$position))
    if(ncol(genotypes) == 1){
      genotypes <- genotypes[order(snp.meta$position),]
      genotypes <- as.data.frame(genotypes, stringsAsFactors = FALSE)
    }
    else{
      genotypes <- genotypes[order(snp.meta$position),]
    }
    snp.meta <- dplyr::arrange(snp.meta, position)
  }
  
  if(any(colnames(snp.meta) == ".snp.id")){
    if(any(duplicated(snp.meta$.snp.id))){stop("Duplicated .snp.id entries found in snp.meta.\n")}
    snp.meta <- dplyr::relocate(snp.meta, .snp.id, .after = dplyr::last_col())
  }
  else{
    snp.meta <- cbind(snp.meta, .snp.id = 1:nrow(snp.meta))
  }
  if(any(colnames(sample.meta) == ".sample.id")){
    if(any(duplicated(sample.meta$.sample.id))){stop("Duplicated .sample.id entries found in sample.meta.\n")}
    sample.meta <- dplyr::relocate(sample.meta, .sample.id, .after = dplyr::last_col())
    
  }
  else{
    sample.meta <- cbind(sample.meta, .sample.id = 1:nrow(sample.meta))
  }
  
  # fix factors
  sample.meta <- dplyr::mutate_if(.tbl = sample.meta, is.factor, as.character)
  snp.meta <- dplyr::mutate_if(.tbl = snp.meta, is.factor, as.character)
  genotypes <- dplyr::mutate_if(.tbl = genotypes, is.factor, as.character)
  
  #===========format and calculate some basics=========
  rownames(genotypes) <- 1:nrow(genotypes)
  rownames(snp.meta) <- 1:nrow(snp.meta)

  gs <- .tabulate_genotypes(genotypes, mDat = mDat, verbose = verbose)

  x <- methods::new("snpRdata", .Data = genotypes, sample.meta = sample.meta, snp.meta = snp.meta,
                    facet.meta = cbind(data.frame(facet = rep(".base", nrow(gs$gs)),
                                                  subfacet = rep(".base", nrow(gs$gs)),
                                                  facet.type = rep(".base", nrow(gs$gs)),
                                                  stringsAsFactors = FALSE),
                                       snp.meta),
                    geno.tables = gs,
                    mDat = mDat,
                    ploidy = 2,
                    bi_allelic = TRUE,
                    data.type = "genotypic",
                    stats = cbind(data.table::data.table(facet = rep(".base", nrow(gs$gs)),
                                                         subfacet = rep(".base", nrow(gs$gs)),
                                                         facet.type = rep(".base", nrow(gs$gs)),
                                                         stringsAsFactors = FALSE),
                                  snp.meta),
                    snp.form = nchar(genotypes[1,1]), row.names = rownames(genotypes),
                    sn = list(sn = NULL, type = NULL),
                    facets = ".base",
                    facet.type = ".base",
                    citations = list(snpR = list(key = "hemstromSnpRUserFriendly2023", details = "snpR package")))
  
  x@calced_stats$.base <- character()

  starting_snps <- nrow(x)

  if(!.skip_filters){
    # run essential filters (np, bi-al), since otherwise many of the downstream applications, including ac formatting, will be screwy.
    if(verbose){cat("Input data will be filtered to remove non bi-allelic data.\n")}
    .make_it_quiet(x <- filter_snps(x, non_poly = FALSE))
  }
  ending_snps <- nrow(x)
  if(starting_snps != ending_snps){
    warning(paste0(starting_snps - ending_snps, " SNPs were removed from the input data set because they were not bi-allelic.\n"))
  }
  
  # add basic maf
  .make_it_quiet(x <- calc_maf(x))
  
  # add ac
  # .make_it_quiet(x@ac <- format_snps(x, "ac")[,c("n_total", "n_alleles", "ni1", "ni2")])
  
  if(is.data.frame(.pass_filters)){
    x@filters <- rbind(x@filters, .pass_filters)
  }
  
  #========return=========
  return(x)
}


#' Pull calculated statistics from snpRdata objects.
#'
#' A convenience function that pulls statistics of any specified type at
#' particular facets from a snpRdata object.
#'
#' Facets are specified as described in \code{\link{Facets_in_snpR}}. If facets
#' = "all", data for all facets, including the base facet, will be returned. By
#' default, the base facet alone will be returned.
#'
#' Different statistics are returned either by named statistic or by type.
#'
#' @section Types:
#'
#'   \itemize{ \item{single: } non-pairwise, non-window statistics (pi, ho, 
#'   etc.)
#'   \item{pairwise: } pairwise, non-window statistics (Fst).
#'   \item{single.window: } non-pairwise, sliding window statistics.
#'   \item{pairwise.window: } pairwise, sliding window statistics. \item{LD: }
#'   linkage disequilibrium matrices and tables. \item{bootstraps: } bootstraps
#'   of window statistics. \item{genetic_distance: } genetic distances
#'   \item{allele_frequency_matrix: } allele frequency matrices. \item{geo_dist:
#'   } geographic distances. \item{ibd: } isolation by distance results.
#'   \item{sample: } sample stats. \item{pop: } population summary statistics.
#'   \item{weighted.means: } Weighted means of SNP statistics per facet level.}
#'
#' @section Statistics: 
#' 
#' \itemize{ \item{ho: } Observed heterozygosity, via
#'   \code{\link{calc_ho}}. \item{pi: } pi (expected number of pairwise
#'   differences, via \code{\link{calc_pi}}). \item{maf: } minor allele
#'   frequencies (and major/minor counts and identities), via
#'   \code{\link{calc_maf}}. \item{private: } private allele identities, via
#'   \code{\link{calc_private}}. \item{association: } results from phenotypic
#'   association tests, via \code{\link{calc_association}}. \item{hwe: }
#'   Hardy-Weinberg Equilibrium p-values, via \code{\link{calc_hwe}}.
#'   \item{tajimas_d: } Tajima's D, Watterson's Theta, and Tajima's Theta, via
#'   \code{\link{calc_tajimas_d}}. \item{fst: } Pairwise Fst, via
#'   \code{\link{calc_pairwise_fst}}. \item{het_hom_ratio: }
#'   Heterozygote/Homozygote ratios within individuals, via
#'   \code{\link{calc_het_hom_ratio}}. \item{ne: } Effective population size
#'   estimates, via \code{\link{calc_ne}}. \item{ld: } Pairwise linkage
#'   disequilibrium, via \code{\link{calc_pairwise_ld}}.
#'   \item{genetic_distances: } Genetic distances,
#'   \code{\link{calc_genetic_distances}}. \item{isolation_by_distance: }
#'   Genetic isolation by distance metrics, via
#'   \code{\link{calc_isolation_by_distance}}. \item{geographic_distance: }
#'   Geographic distances between samples--does not use genetic data. Calculated
#'   during \code{\link{calc_isolation_by_distance}}, but fetchable
#'   independently here.
#'   \item{random_forest} Random forest snp-specific responses, see
#'   \code{\link{run_random_forest}}.
#'
#'   }
#'
#' @param x snpRdata object.
#' @param facets character or NULL, default NULL. Facets for which to fetch
#'   data.
#' @param stats character or NULL, default NULL. Statistics for which to fetch
#'   data. Named identically as in the function used to calculate them (such as
#'   calc_pi). See description for list. Alternatively, a type of statistic can
#'   be requested--see description for details.
#' @param bootstraps logical, default FALSE. If FALSE, bootstraps will not be
#'   returned even if they have been created for the requested statistics. Since
#'   bootstrapped datasets are often quite large, returning these may use a lot
#'   of additional memory.
#'
#'   
#' @export
#' @author William Hemstrom
#'   
#' @examples
#' # generate some statistics
#' dat <- calc_pi(stickSNPs, "pop")
#' dat <- calc_pairwise_fst(dat, "pop")
#'
#' # fetch pi
#' get.snpR.stats(dat, "pop", "pi")
#' # fetch fst
#' get.snpR.stats(dat, "pop", "fst")
#' # fetch both
#' get.snpR.stats(dat, "pop", c("pi", "fst"))
#' 
#' # return a type of statistic instead of specific statistics
#' dat <- calc_ho(stickSNPs, "chr.pop")
#' get.snpR.stats(dat, "chr.pop", "weighted.means")
#' 
get.snpR.stats <- function(x, facets = NULL, stats = "single", bootstraps = FALSE){
  #===========quick and easy if referencing storage slot===============
  if(length(stats) == 1 & stats[1] %in% 
     c("single", "pairwise", "single.window", "pairwise.window", "LD", "bootstraps", "genetic_distance",
       "allele_frequency_matrix", "geo_dist", "ibd", "sample", "pop", "weighted.means")){
    return(.get.snpR.stats(x, facets, type = stats))
  }
  
  #===========aliases=========
  stats <- tolower(stats)
  aliases <- c(ibd = "isolation_by_distance",
               tsd = "tajimas_d",
               d = "tajimas_d",
               ar = "allelic_richness",
               richness = "allelic_richness",
               ROH = "roh",
               fROH = "roh",
               FROH = "roh",
               Froh = "roh")
  
  need_unaliased <- which(stats %in% names(aliases))
  
  if(length(need_unaliased) > 0){
    stats[need_unaliased] <- aliases[stats[need_unaliased]]
  }
  
  #====================sanity checks=====================
  msg <- character()
  if(!is.snpRdata(x)){
    stop("x must be a snpRdata object.\n")
  }
  
  if(is.null(facets[1])){
    facets <- ".base"
  }
  if(facets[1] == "all"){
    facets <- x@facets
  }
  
  facets <- .check.snpR.facet.request(x, facets, "none")
  
  stats <- unique(stats)
  
  if(!all(stats %in% names(.internal.data$statistic_index))){
    msg <- c(msg, paste0("Requested statistics: ", paste0(stats[which(!stats %in% names(.internal.data$statistic_index))], collapse = ", "),
                         "\nnot recognized. \nRecognized statistics: ",
                         paste0(names(.internal.data$statistic_index), collapse = ", "),
                         "."))
  }
  if(length(msg) > 0){
    stop(paste0(msg, collapse = "\n"))
  }
  #=======================fetch, using the internal version of .get.snpR.data========================
  out <- vector("list", length(stats))
  names(out) <- stats

  # determine the parts of the snpR object we need to harvest
  needed.parts <- purrr::map(.internal.data$statistic_index[stats], "types")
  unique.needed.parts <- unique(unlist(needed.parts))
  if(!bootstraps){
    if(any(unique.needed.parts == "bootstraps")){
      unique.needed.parts <- unique.needed.parts[-which(unique.needed.parts == "bootstraps")]
    }
  }
  out <- vector("list", length(unique.needed.parts))
  names(out) <- unique.needed.parts
  for(i in 1:length(unique.needed.parts)){
    need.this.part <- names(needed.parts)[grep(unique.needed.parts[i], needed.parts)]
    col_patterns <-  unlist(purrr::map(.internal.data$statistic_index[need.this.part], "col_pattern"))
    if(unique.needed.parts[i] == "bootstraps"){
      col_patterns <- c("nk", "step", "stat", "value")
    }
    if(!is.na(col_patterns[1])){
      suppressWarnings(res <- .get.snpR.stats(x, facets, type = unique.needed.parts[i],
                                              col_pattern = col_patterns))
    }
    else{
      suppressWarnings(res <- .get.snpR.stats(x, facets, type = unique.needed.parts[i]))
    }
    
    if(unique.needed.parts[i] == "bootstraps"){
      res <- res[res$stat %in% unlist(purrr::map(.internal.data$statistic_index[need.this.part], "col_pattern")),]
    }
    
    if(!is.null(res)){
      out[[i]] <- res
    }
  }
  
  out <- Filter(Negate(is.null), out)
  if(length(out) == 0){
    warning("No statistics located for requested stats/facets combination(s). Did you forget to run the requested stats or mistype the facets or statistic name (which should match the 'calc_' function used)?\n")
  }
  
  return(out)
}

# Internal version of get.snpR.stats, for backwards comparability if type is set or
# for fetching data for get.snpR.stats.
# 
# See documentation for get.snpR.stats.
.get.snpR.stats <- function(x, facets = NULL, type = "single", col_pattern = NULL){
  ..keep.cols <- subfacet <- facet <- ..grep.cols <- NULL
  
  # sanity check
  if(!is.snpRdata(x)){
    stop("x must be a snpRdata object.\n")
  }
  
  if(is.null(facets[1])){
    facets <- ".base"
  }
  if(facets[1] == "all"){
    facets <- x@facets
  }
  
  good.types <- c("single", "pairwise", "single.window", "pairwise.window", "LD", "bootstraps", "genetic_distance",
                  "allele_frequency_matrix", "geo_dist", "ibd", "sample", "pop", "weighted.means", "fst.matrix", "roh")
  if(!type %in% good.types){
    stop("Unaccepted stats type. Options: ", paste0(good.types, collapse = ", "), ".\nSee documentation for details.\n")
  }
  
  facets <- .check.snpR.facet.request(x, facets, "none")
  # bad.facets <- which(!facets %in% x@facets)
  # if(length(bad.facets) > 0){
  #   stop("No data statistics calculated for facets: ", paste0(facets[bad.facets], collapse = ", "), ".\n")
  # }
  
  
  #=======subfunctions======
  extract.basic <- function(y, facets, type = "standard", col_pattern = NULL){
    if(nrow(y) == 0){
      warning("No statistics calculated for provided type.\n")
      return(NULL)
    }
    if(type == "standard"){
      keep.rows <- which(y$facet %in% facets)
      keep.cols <- which(!colnames(y) %in% c("facet.type"))
    }
    else if(type == "window"){
      facets <- .check.snpR.facet.request(x, facets, "none", T)
      keep.rows <- numeric()
      keep.cols <- 1:ncol(y)
      
      # for each facet, figure out which rows conform
      for(i in 1:length(facets[[1]])){
        if(facets[[2]][i] == "snp"){
          keep.rows <- c(keep.rows, which(y$snp.facet == facets[[1]][i] & y$facet == ".base"))
        }
        else if(facets[[2]][i] == "sample"){
          keep.rows <- c(keep.rows, which(y$snp.facet == ".base" & y$facet == facets[[1]][i]))
        }
        else if(facets[[2]][i] == ".base"){
          keep.rows <- c(keep.rows, which(y$snp.facet == ".base" & y$facet == ".base"))
        }
        else{
          tfacet <- unlist(.split.facet(facets[[1]][i]))
          tfacet <- .check.snpR.facet.request(x, tfacet, "none", T)
          
          # need to paste together any snp or sample faces
          sample.facets <- paste0(tfacet[[1]][which(tfacet[[2]] == "sample")], collapse = ".")
          sample.facets <- .check.snpR.facet.request(x, sample.facets)
          snp.facets <- paste0(tfacet[[1]][which(tfacet[[2]] == "snp")], collapse = ".")
          snp.facets <- .check.snpR.facet.request(x, snp.facets, remove.type = "none")
          
          keep.rows <- c(keep.rows, which(y$snp.facet == snp.facets & y$facet == sample.facets))
        }
      }
    }
    else if(type == "comingled"){ # this type has facet and snp.facet info, mingles everything (as in weighted.means)
      facets <- .check.snpR.facet.request(x, facets, "none", T)
      
      keep.rows <- numeric()
      for(i in 1:length(facets[[1]])){
        split.part <- unlist(.split.facet(facets[[1]][i]))
        split.part <- .check.snpR.facet.request(x, split.part, remove.type = "none", TRUE)
        snp.part <- split.part[[1]][which(split.part[[2]] == "snp")]
        if(length(snp.part) == 0){
          snp.part <- ".base"
        }
        samp.part <- split.part[[1]][which(split.part[[2]] == "sample")]
        if(length(samp.part) == 0){
          samp.part <- ".base"
        }
        samp.part <- paste0(samp.part, collapse = ".")
        snp.part <- paste0(snp.part, collapse = ".")
        
        keep.rows <- c(keep.rows, which(y$facet == samp.part & y$snp.facet == snp.part))
      }
      
      keep.rows <- sort(unique(keep.rows))
      
      # remove empty columns
      empty.cols <- which(colSums(is.na(y[keep.rows,])) == nrow(y[keep.rows,]))
      keep.cols <- 1:ncol(y)
      if(any(empty.cols)){
        keep.cols <- keep.cols[-empty.cols]
      }
    }
    
    
    # adjust keep.cols to remove any unwanted stats
    if(!is.null(col_pattern)){
      keep.cols <- colnames(y)[keep.cols]
      
      good.cols <- keep.cols[which(keep.cols %in% c("facet", "subfacet", "snp.facet", "snp.subfacet", "comparison", "pop", "sigma", "step", "nk.status", "gaussian", "triple_sigma", "n_snps",
                                                    colnames(x@facet.meta), colnames(sample.meta(x))))]
      grep.cols <- numeric(0)
      for(i in 1:length(col_pattern)){
        grep.cols <- c(grep.cols, grep(col_pattern[i], colnames(y)))
      }
      if(length(grep.cols) == 0){
        return(NULL)
      }
      grep.cols <- colnames(y)[grep.cols]
      
      keep.cols <- which(colnames(y) %in% c(good.cols, grep.cols))
      
      # remove any rows that only contain NA values for the grep cols
      empty.rows <- which(rowSums(is.na(.fix..call(y[keep.rows, ..grep.cols]))) == length(grep.cols))
      if(length(empty.rows) > 0){
        keep.rows <- keep.rows[-empty.rows]
      }
    }
    

    
    if(length(keep.rows) == 0){
      warning("No statistics calculated for this facet and statistics type.\n")
      return(NULL)
    }
    
    if(all(colnames(y)[keep.cols] %in% c("facet", "subfacet", "snp.facet", "snp.subfacet", ".snp.id", "comparison"))){
      return(NULL)
    }
    return(as.data.frame(.fix..call(y[keep.rows, ..keep.cols]), stringsAsFactors = FALSE))
  }
  
  extract.LD <- function(y, facets){
    if(length(y) == 0){return(NULL)} # for window only
    # get sample facets
    samp.facets <- .check.snpR.facet.request(x, facets)
    if(length(samp.facets) != length(facets)){
      samp.facets <- c(samp.facets, ".base")
    }
    
    # get prox
    prox <- y$prox[y$prox$sample.facet %in% samp.facets,]
    if(nrow(prox) == 0){
      stop("No LD values calculated for these facets.\n")
    }
    
    # get matrices
    matrices <- y$LD_matrices[[which(names(y$LD_matrices) %in% facets)]]
    
    return(list(prox = prox, matrices = matrices))
  }
  
  extract.gd.afm <- function(y, facets) y[which(names(y) %in% facets)]
  
  extract.fst.matrix <- function(x, facets = NULL){
    p1 <- p2 <- NULL
    facets <- .check.snpR.facet.request(x, facets, "complex", return_base_when_empty = F, fill_with_base = F)
    if(length(facets) == 0){
      return(NULL)
    }
    fst <- data.table::as.data.table(.get.snpR.stats(x, facets, "weighted.means"), col_pattern = .internal.data$statistic_index$fst$col_pattern)
    if(nrow(fst) == 0){
      return(NULL)
    }
    fst[,c("p1", "p2") := tstrsplit(subfacet, "~")]
    
    cats <- unique(fst$facet)
    res <- vector("list", length(cats))
    names(res) <- cats
    
    not.pairwise <- which(is.na(fst$p1) | is.na(fst$p2))
    if(length(not.pairwise) > 0){
      fst <- fst[-not.pairwise,]
    }
    
    empties <- numeric(0)
    for(i in 1:length(cats)){
      tfst <- data.table::copy(fst[facet == cats[i],])
      levs <- unique(c(tfst$p1, tfst$p2))
      levs <- sort(levs)
      tfst[,p1 := factor(p1, levs)]
      tfst[,p2 := factor(p2, levs)]
      
      if(sum(fst$facet == cats[i]) == 0){
        empties <- c(empties, i)
        next
      }
      else if("weighted_mean_fst_p" %in% colnames(fst)){
        res[[i]] <- list(fst = data.table::dcast(tfst, p1~p2, value.var = "weighted_mean_fst"),
                         p = data.table::dcast(tfst, p1~p2, value.var = "weighted_mean_fst_p"))
      }
      else{
        res[[i]] <- data.table::dcast(tfst, p1~p2, value.var = "weighted_mean_fst")
      }
    }
    
    # remove any empty elements
    if(length(empties) > 0){
      res[[empties]] <- NULL
    }
    
    return(res)
  }
  
  extract.roh <- function(x, facets = NULL){
    nfacets <- .check.snpR.facet.request(x, facets, remove.type = "sample")
    if(any(sort(nfacets) != sort(.check.snpR.facet.request(x, facets, "none")))){
      warning(paste0("ROH is only calculated for snp facets; requested facets simplified to: ", paste0(nfacets, collapse = ",")))
    }
    res <- x@other$roh[which(x@other$roh$snp.facet %in% nfacets),]
    return(res)
  }
  
  #========prep=============
  if(!is.null(facets)){
    if(facets[1] == "all"){
      facets <- x@facets
    }
  }
  else {
    facets <- ".base"
  }
  
  facets <- .check.snpR.facet.request(x, facets, "none")
  
  #========extract data======
  if(type == "single"){
    facets <- .check.snpR.facet.request(x, facets)
    return(extract.basic(x@stats, facets, col_pattern = col_pattern))
  }
  else if(type == "pairwise"){
    facets <- .check.snpR.facet.request(x, facets)
    base.facets <- which(facets == ".base")
    if(length(base.facets) > 0){
      stop("Cannot find pairwise statistics without any facets (facets = NULL or '.base').\n")
    }
    return(extract.basic(x@pairwise.stats, facets, col_pattern = col_pattern))
  }
  else if(type == "single.window"){
    return(extract.basic(x@window.stats, facets, "window", col_pattern = col_pattern))
  }
  else if(type == "pairwise.window"){
    return(extract.basic(x@pairwise.window.stats, facets, "window", col_pattern = col_pattern))
  }
  else if(type == "LD"){
    return(extract.LD(x@pairwise.LD, facets))
  }
  else if(type == "bootstraps"){
    return(extract.basic(x@window.bootstraps, facets, "window", col_pattern = col_pattern))
  }
  else if(type == "genetic_distance"){
    return(extract.gd.afm(x@genetic_distances, facets))
  }
  else if(type == "allele_frequency_matrix"){
    return(extract.gd.afm(x@allele_frequency_matrices, facets))
  }
  else if(type == "geo_dist"){
    return(extract.gd.afm(x@other$geo_dist, .check.snpR.facet.request(x, facets)))
  }
  else if(type == "ibd"){
    return(extract.gd.afm(x@other$ibd, facets))
  }
  else if(type == "sample"){
    facets <- .check.snpR.facet.request(x, facets, "sample")
    return(extract.basic(x@sample.stats, facets, 
                         col_pattern = col_pattern, 
                         type = "comingled"))
  }
  else if(type == "pop"){
    return(extract.basic(x@pop.stats, facets, col_pattern = col_pattern))
  }
  else if(type == "weighted.means"){
    return(extract.basic(x@weighted.means, facets, type = "comingled", col_pattern = col_pattern))
  }
  else if(type == "fst.matrix"){
    return(extract.fst.matrix(x, facets))
  }
  else if(type == "roh"){
    return(extract.roh(x,facets))
  }
  
}
