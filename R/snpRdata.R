#'@import data.table 
#'@importFrom foreach %dopar%
NULL


check.snpRdata <- function(object){
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

  # check that the facet meta and snp meta column names match
  if(!identical(colnames(object@snp.meta), colnames(object@facet.meta)[-c(1:3)])){
    errors <- c(errors, "Column names in the snp.meta and stored metadata for facets do not match. This is likely due to later change of the snp metadata which added or changed column names but didn't change facet.meta column names.\n")
  }

  warns <- character()
  # check for bad entries in character columns (periods are fine in numeric!)
  chr.cols.samp <- unlist(lapply(object@sample.meta, class))
  chr.cols.samp <- chr.cols.samp[which(chr.cols.samp == "character")]
  chr.cols.snp <- unlist(lapply(object@snp.meta, class))
  chr.cols.snp <- chr.cols.snp[which(chr.cols.snp == "character")]

  if(length(chr.cols.samp) > 0){
    l1 <- unlist(lapply(object@sample.meta[,names(chr.cols.samp)], grepl, pattern = bad.chars))
    if(sum(l1) > 0){
      msg <- paste0("Some sample metadata columns contain unacceptable special characters. Unaccepted characters: '.', '~', or any whitespace.\nThese can cause unexpected behaviour if the subect columns are used as facets.\n")
      warns <- c(warns, msg)
    }
  }
  if(length(chr.cols.snp) > 0){
    l1 <- unlist(lapply(object@snp.meta[,names(chr.cols.snp)], grepl, pattern = bad.chars))
    if(sum(l1) > 0){
      msg <- paste0("Some snp metadata columns contain unacceptable special characters. Unaccepted characters: '.', '~', or any whitespace.\nThese can cause unexpected behaviour if the subect columns are used as facets.\n")
      warns <- c(warns, msg)
    }
  }
  if(length(warns) > 0){warning(paste0(warns, collapse = "\n"))}

  if(length(errors) == 0){return(TRUE)}
  else{return(errors)}
}



#'Storage class for snpR data and calculated statistics.
#'
#'\code{\link{import.snpR.data}} creates and defines the snpRdata class to store
#'both raw genotype data, sample and locus specific metadata, useful data
#'summaries, repeatedly internally used tables, calculated summary statistics,
#'and smoothed statistic data.
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
#'  (allele counts), and wm (genotype counts + missing counts)
#'@slot ac data.frame containing ac formatted data, see
#'  \code{\link{format_snps}}.
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
#'@slot other list, contains other miscellanious calculated statistics that do
#'  not fit cleanly elsewhere.
#'
#'@author William Hemstrom
#'  
snpRdata <- setClass(Class = 'snpRdata', slots = c(sample.meta = "data.frame",
                                       snp.meta = "data.frame",
                                       facet.meta = "data.frame",
                                       mDat = "character",
                                       snp.form = "numeric",
                                       geno.tables = "list",
                                       ac = "data.frame",
                                       facets = "character",
                                       facet.type = "character",
                                       stats = "data.frame",
                                       window.stats = "data.frame",
                                       pairwise.stats = "data.frame",
                                       pairwise.window.stats = "data.frame",
                                       sample.stats = "data.frame",
                                       pop.stats = "data.frame",
                                       pairwise.LD = "list",
                                       window.bootstraps = "data.frame",
                                       sn = "list",
                                       calced_stats = "list",
                                       allele_frequency_matrices = "list",
                                       genetic_distances = "list",
                                       weighted.means = "data.frame",
                                       other = "list"),
         contains = c(data = "data.frame"),
         validity = check.snpRdata)


#'Import genotype and metadata into a snpRdata object.
#'
#'\code{import.snpR.data} converts genotype and meta data to the snpRdata class,
#'which stores raw genotype data, sample and locus specific metadata, useful
#'data summaries, repeatedly internally used tables, calculated summary
#'statistics, and smoothed statistic data.
#'
#'The snpRdata class is built to contain SNP genotype data for use by functions
#'in the snpR package. It inherits from the S3 class data.frame, in which the
#'genotypes are stored, and can be manipulated identically. It also stores
#'sample and locus specific metadata, genomic summary information, and results
#'from most snpR functions. Genotypes are stored in the "character" format, as
#'output by \code{\link{format_snps}}. Missing data is noted with "NN".
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
#'  (heterozygous), or 2 (homozyogus allele 2). \item{.genepop: } genepop file
#'  format, with genotypes stored as either 4 or 6 numeric characters. Works
#'  only with bi-allelic data. Genotypes will be converted (internally) to NN:
#'  the first allele (numerically) will be coded as A, the second as C.
#'  \item{FSTAT: } FSTAT file format, with genotypes stored as either 4 or 6
#'  numeric characters. Works only with bi-allelic data. Genotypes will be
#'  converted (internally) to NN: the first allele (numerically) will be coded
#'  as A, the second as C. \item{plink: } plink .bed, .fam, and .bim files, via
#'  \code{\link[genio]{read_plink}}. If any of these file types is provided,
#'  snpR via \code{\link[genio]{read_plink}} will look for the other file types
#'  automatically. Sample metadata should be contained in the .fam file and
#'  SNP metadata in the .bim file, so sample or snp meta data provided here will
#'  be ignored.}
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
#'  adegenet. Note, no need to import genpop objects, the equivalent statistics
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
#'  corresponding metadata. \item{ac: } data in allele count format, used
#'  internally. facet.meta contains corresponding metadata. \item{facets: }
#'  vector of the facets that have been added to the data. \item{facet.type: }
#'  classes of the added facets (snp, sample, complex, or .base). \item{stats: }
#'  data.frame containing all calculated non-pairwise single-snp statistics and
#'  metadata. \item{window.stats: } data.frame/table containing all non-pairwise
#'  statistics calculated for sliding windows. \item{pairwise.stats: }
#'  data.frame/table containing all pairwise (fst) single-snp statistics.
#'  \item{pairwise.window.stats: } data.frame/table containing all pairwise
#'  statistics calculated for sliding windows. \item{sample.stats: }
#'  data.frame/table containing statistics caluclated for each individual
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
#'  
#'@examples
#' # import example data as a snpRdata object
#' # produces data identical to that contained in the stickSNPs example dataset.
#' genos <- stickRAW[,-c(1:3)]
#' snp_meta <- stickRAW[,1:3]
#' sample_meta <- data.frame(pop = substr(colnames(stickRAW)[-c(1:3)], 1, 3), fam = rep(c("A", "B", "C", "D"), length = ncol(stickRAW) - 3), stringsAsFactors = FALSE)
#' import.snpR.data(genos, snp.meta = snp_meta, sample.meta = sample_meta, mDat = "NN")
#'
#' # from an adegenet genind object
#' ex.genind  <- adegenet::df2genind(t(stickRAW[,-c(1:3)]), ncode = 1, NA.char = "N") # get genind data
#' import.snpR.data(ex.genind, snp_meta, sample_meta) # note, will add whatever metadata data is in the genind object to the snpRdata object. Could be run without the snp or sample metadatas.
#'
#' # from an adegenet genlight object
#' ## create a dummy dataset, add some metadata
#' dat <- lapply(1:50, function(i) sample(c(0,1,2, NA), 1000, prob=c(.25, .49, .25, .01), replace=TRUE))
#' names(dat) <- paste("indiv", 1:length(dat))
#' print(object.size(dat), unit="aut") # size of the original data
#' genlight <- new("genlight", dat) # conversion
#' newalleles <- character(adegenet::nLoc(genlight))
#' for(i in 1:length(newalleles)){
#'   newalleles[i] <- paste0(sample(c("a", "c", "g", "t"), 2, FALSE), collapse = "/")
#' }
#' adegenet::alleles(genlight) <- newalleles
#' adegenet::pop(genlight) <- sample(LETTERS[1:4], 50, TRUE)
#' adegenet::position(genlight) <- sample(100000, 1000, FALSE)
#' adegenet::chr(genlight) <- sample(10, 1000, TRUE)
#'
#' ## run the conversion
#' dat <- import.snpR.data(genlight)
#'
#' \dontrun{
#' ## not run:
#' # from a file:
#' dat <- import.snpR.data(system.file("extdata", "stick_NN_input.txt", package = "snpR"), drop = 1:3) # note that the drop argument is passed to data.table::fread!
#' # if wanted, snp and sample metadata could be provided as usual.
#' 
#' ## not run:
#' # from plink:
#' # make plink data
#' format_snps(stickSNPs, "plink", outfile = "plink_test", chr = "group")
#'
#' # read plink
#' dat <- import.snpR.data("plink_test.bed")
#' }
#'
#'@export
#'
#'@author William Hemstrom
import.snpR.data <- function(genotypes, snp.meta = NULL, sample.meta = NULL, mDat = "NN", chr.length = NULL,
                             ...){
  #======special cases========
  # sample and snp metadata
  if(is.character(sample.meta)){
    if(file.exists(sample.meta)){
      sample.meta <- as.data.frame(data.table::fread(sample.meta))
    }
    else{
      stop("Cannot locate sample.meta file.\n")
    }
  }
  if(is.character(snp.meta)){
    if(file.exists(snp.meta)){
      snp.meta <- as.data.frame(data.table::fread(snp.meta))
    }
    else{
      stop("Cannot locate snp.meta file.\n")
    }
  }
  
  # genotypes
  if("genind" %in% class(genotypes)){
    return(genind.to.snpRdata(genotypes, snp.meta, sample.meta))
  }
  if("genlight" %in% class(genotypes)){
    return(genlight.to.snpRdata(genotypes, snp.meta, sample.meta))
  }
  if("vcfR" %in% class(genotypes)){
    return(process_vcf(genotypes, snp.meta, sample.meta))
  }
  if(is.matrix(genotypes)){
    genotypes <- as.data.frame(genotypes)
  }
  if(is.character(genotypes) & length(genotypes) == 1){
    if(file.exists(genotypes)){
      # check for ms or vcf, etc file
      if(grepl("\\.vcf$", genotypes) | grepl("\\.vcf\\.gz$", genotypes)){
        return(process_vcf(genotypes, snp.meta, sample.meta))
      }
      else if(grepl("\\.ms$", genotypes)){
        return(format_snps(genotypes, input_format = "ms", snp.meta = snp.meta, sample.meta = sample.meta, chr.length = chr.length))
      }
      else if(grepl("\\.genepop$", genotypes)){
        return(process_genepop(genotypes, snp.meta, sample.meta, mDat))
      }
      else if(grepl("\\.fstat$", genotypes)){
        return(process_FSTAT(genotypes, snp.meta, sample.meta, mDat))
      }
      else if(grepl("\\.bim$", genotypes) | grepl("\\.fam$", genotypes) | grepl("\\.bed$", genotypes)){
        check.installed("tools")
        return(process_plink(tools::file_path_sans_ext(genotypes)))
      }
      else{
        genotypes <- as.data.frame(data.table::fread(genotypes, ...))
      }
    }
    else{
      stop("File not found. Fix path or import manually and provide to import.snpR.data.\n")
    }
  }
  
  
  #=================check input format for non-special case=============================
  # NN, no need to do anything, just read in and proceed as normal.
  if(genotypes[1,1] %in% 
     c(apply(expand.grid(c("A", "T", "C", "G"), c("A", "T", "C", "G")), 1, paste, collapse=""), mDat)){
    cat("Assuming data is in NN format.\n")
  }
  
  # sn
  else if(genotypes[1,1] %in% c(0, 1, 2, mDat)){
    cat("Assuming single nucleotide format.\n")
    return(format_snps(genotypes, input_format = "sn", input_mDat = mDat, sample.meta = sample.meta, snp.meta = snp.meta))
  }
  
  # 0000
  else if(genotypes[1,1] %in% 
          c(apply(expand.grid(c("01", "02", "03", "04"), c("01", "02", "03", "04")), 1, paste, collapse=""), mDat)){
    cat("Assuming 0000 format.\n")
    return(format_snps(genotypes, input_format = "0000", input_mDat = mDat, sample.meta = sample.meta, snp.meta = snp.meta))
  }
  
  # SNP_tab
  else if(genotypes[1,1] %in% 
          c(apply(expand.grid(c("A", "T", "C", "G"), c("A", "T", "C", "G")), 1, paste, collapse=" "), mDat, c("A", "T", "C", "G"))){
    cat("Assuming snp_tab format.\n")
    return(format_snps(genotypes, input_format = "snp_tab", input_mDat = mDat, sample.meta = sample.meta, snp.meta = snp.meta))
  }
  
  #couldn't find a supported format
  else{
    good.formats <- c(".vcf", ".vcf.gz", ".ms", ".bim", ".bed", ".fam", ".fstat", ".genepop", "NN", "0000", "sn", "snp_tab")
    stop(paste0(c("Unsupported file format. Supported formats: ", good.formats, "\n"), collapse = " "))
  }
  
  
  
  #============sanity checks and prep========
  if(is.null(snp.meta)){
    snp.meta <- data.frame(snpID = paste0("snp", 1:nrow(genotypes)))
  }
  if(is.null(sample.meta)){
    sample.meta <- data.frame(sampID = paste0("samp", 1:ncol(genotypes)))
  }
  
  # prepare things for addition to data
  if(any(is.na(genotypes))){
    stop("NA found in input genotypes. Often, this is in the last row or column.\n")
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
  }
  else{
    snp.meta <- cbind(snp.meta, .snp.id = 1:nrow(snp.meta))
  }
  if(any(colnames(sample.meta) == ".sample.id")){
    if(any(duplicated(sample.meta$.sample.id))){stop("Duplicated .sample.id entries found in sample.meta.\n")}
  }
  else{
    sample.meta <- cbind(sample.meta, .sample.id = 1:nrow(sample.meta))
  }
  
  # fix factors
  sample.meta <- dplyr::mutate_if(.tbl = sample.meta, is.factor, as.character)
  snp.meta <- dplyr::mutate_if(.tbl = snp.meta, is.factor, as.character)
  genotypes <- dplyr::mutate_if(.tbl = genotypes, is.factor, as.character)
  
  
  # warn if anything repeated across sample level factors
  uniques <- lapply(sample.meta, unique)
  uniques <- uniques[-which(names(uniques) == ".sample.id")]
  uniques <- unlist(uniques)
  if(any(duplicated(uniques))){
    warning(paste0("Some levels are duplicated across multiple sample meta facets.\nThis will cause issues if those sample facets are run during analysis.\nIssues:\n",
                paste0(uniques[which(duplicated(uniques))], "\n", collapse = "")))
  }
  
  
  #===========format and calculate some basics=========
  rownames(genotypes) <- 1:nrow(genotypes)
  rownames(snp.meta) <- 1:nrow(snp.meta)

  gs <- tabulate_genotypes(genotypes, mDat = mDat, verbose = TRUE)
  
  x <- new("snpRdata", .Data = genotypes, sample.meta = sample.meta, snp.meta = snp.meta,
           facet.meta = cbind(data.frame(facet = rep(".base", nrow(gs$gs)),
                                   subfacet = rep(".base", nrow(gs$gs)),
                                   facet.type = rep(".base", nrow(gs$gs)),
                                   stringsAsFactors = FALSE),
                              snp.meta),
           geno.tables = gs,
           mDat = mDat,
           stats = cbind(data.table::data.table(facet = rep(".base", nrow(gs$gs)),
                                                subfacet = rep(".base", nrow(gs$gs)),
                                                facet.type = rep(".base", nrow(gs$gs)),
                                                stringsAsFactors = FALSE),
                         snp.meta),
           snp.form = nchar(genotypes[1,1]), row.names = rownames(genotypes),
           sn <- list(sn = NULL, type = NULL),
           facets = ".base",
           facet.type = ".base",
           calced_stats = list(),
           allele_frequency_matrices = list(),
           genetic_distances = list(),
           weighted.means = data.frame(),
           other = list())
  
  x@calced_stats$.base <- character()
  
  
  # run essential filters (np, bi-al), since otherwise many of the downstream applications, including ac formatting, will be screwy.
  cat("Input data will be filtered to remove non bi-allelic data.\n")
  invisible(capture.output(x <- filter_snps(x, non_poly = FALSE)))
  
  # add basic maf
  invisible(capture.output(x <- calc_maf(x)))
  
  # add ac
  invisible(capture.output(x@ac <- format_snps(x, "ac")[,c("n_total", "n_alleles", "ni1", "ni2")]))
  
  
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
#'   \itemize{ \item{single: } non-pairwise, non-window statitics (pi, ho, ect.)
#'   \item{pairwise: } pairwise, non-window statistics (Fst).
#'   \item{single.window: } non-pairwise, sliding window statistics.
#'   \item{pairwise.window: } pairwise, sliding window statistics. \item{LD: }
#'   linkage disequilibrium matrices and tables. \item{bootstraps: } bootstraps
#'   of window statistics. \item{genetic_distance: } genetic distances
#'   \item{allele_frequency_matrix: } allele frequency matrices. \item{geo_dist:
#'   } geographic distances. \item{ibd: } isolation by distance results.
#'   \item{sample: } sample stats. \item{pop: } population summary statistics. }
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
#'   Hardy-Weinburg Equilibrium p-values, via \code{\link{calc_hwe}}.
#'   \item{tajimas_d: } Tajima's D, Waterson's Theta, and Tajima's Theta, via
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
#' dat <- calc_pi(stickSNPs, "group.pop")
#' dat <- calc_pairwise_fst(stickSNPs, "group.pop")
#'
#' # fetch pi
#' get.snpR.stats(dat, "group.pop", "pi")
#' # fetch fst
#' get.snpR.stats(dat, "group.pop", "fst")
#' 
#' # return a type of statistic instead of specific statistics
#' dat <- calc_ho(stickSNPs, "group.pop")
#' get.snpR.stats(dat, "group.pop", "weighted.means")
#' 
get.snpR.stats <- function(x, facets = NULL, stats = NULL, bootstraps = FALSE){
  if(length(stats) == 1 & stats[1] %in% 
     c("single", "pairwise", "single.window", "pairwise.window", "LD", "bootstraps", "genetic_distance",
       "allele_frequency_matrix", "geo_dist", "ibd", "sample", "pop", "weighted.means")){
    return(.get.snpR.stats(x, facets, type = stats))
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
  
  facets <- check.snpR.facet.request(x, facets, "none")
  
  if(!all(stats %in% names(statistic_index))){
    msg <- c(msg, paste0("Requested statistics: ", paste0(stats[which(!stats %in% names(statistic_index))], collapse = ", "),
                         "\nnot recognized. \nRecognized statistics: ",
                         paste0(names(statistic_index), collapse = ", "),
                         "."))
  }
  if(length(msg) > 0){
    stop(paste0(msg, collapse = "\n"))
  }
  #=======================fetch, using the internal version of .get.snpR.data========================
  out <- vector("list", length(stats))
  names(out) <- stats

  # determine the parts of the snpR object we need to harvest
  needed.parts <- purrr::map(statistic_index[stats], "types")
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
    col_patterns <-  unlist(purrr::map(statistic_index[need.this.part], "col_pattern"))
    if(!is.na(col_patterns[1])){
      suppressWarnings(res <- .get.snpR.stats(x, facets, type = unique.needed.parts[i],
                                              col_pattern = col_patterns))
    }
    else{
      suppressWarnings(res <- .get.snpR.stats(x, facets, type = unique.needed.parts[i]))
    }
    
    if(!is.null(res)){
      out[[i]] <- res
    }
  }
  
  out <- Filter(Negate(is.null), out)
  
  return(out)
}

#' Internal version of get.snpR.stats, for backwards comparability if type is set or
#' for fetching data for get.snpR.stats.
#' 
#' See documentation for get.snpR.stats.
.get.snpR.stats <- function(x, facets = NULL, type = "single", col_pattern = NULL){
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
                  "allele_frequency_matrix", "geo_dist", "ibd", "sample", "pop", "weighted.means", "fst.matrix")
  if(!type %in% good.types){
    stop("Unaccepted stats type. Options: ", paste0(good.types, collapse = ", "), ".\nSee documentation for details.\n")
  }
  
  facets <- check.snpR.facet.request(x, facets, "none")
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
      facets <- check.snpR.facet.request(x, facets, "none", T)
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
          tfacet <- unlist(strsplit(facets[[1]][i], "(?<!^)\\.", perl = TRUE))
          tfacet <- check.snpR.facet.request(x, tfacet, "none", T)
          
          # need to paste together any snp or sample faces
          sample.facets <- paste0(tfacet[[1]][which(tfacet[[2]] == "sample")], collapse = ".")
          sample.facets <- check.snpR.facet.request(x, sample.facets)
          snp.facets <- paste0(tfacet[[1]][which(tfacet[[2]] == "snp")], collapse = ".")
          snp.facets <- check.snpR.facet.request(x, snp.facets, remove.type = "none")
          
          keep.rows <- c(keep.rows, which(y$snp.facet == snp.facets & y$facet == sample.facets))
        }
      }
    }
    else if(type == "comingled"){ # this type has facet and snp.facet info, mingles everything (as in weighted.means)
      facets <- check.snpR.facet.request(x, facets, "none", T)
      
      keep.rows <- numeric()
      for(i in 1:length(facets[[1]])){
        split.part <- unlist(strsplit(facets[[1]][i], split = "(?<!^)\\.", perl = T))
        split.part <- check.snpR.facet.request(x, split.part, remove.type = "none", TRUE)
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
      good.cols <- which(colnames(y)[keep.cols] %in% c("facet", "subfacet", "snp.facet", "snp.subfacet",
                                                       colnames(x@facet.meta)))
      grep.cols <- numeric(0)
      for(i in 1:length(col_pattern)){
        grep.cols <- c(grep.cols, grep(col_pattern[i], colnames(y)))
      }
      if(length(grep.cols) == 0){
        return(NULL)
      }
      
      keep.cols <- keep.cols[which(keep.cols %in% c(good.cols, grep.cols))]
    }
    
    if(length(keep.rows) == 0){
      # return(NULL)
      stop("No statistics calculated for this facet and statistics type.\n")
    }
    return(as.data.frame(y[keep.rows, ..keep.cols], stringsAsFactors = FALSE))
  }
  
  extract.LD <- function(y, facets){
    # get sample facets
    samp.facets <- check.snpR.facet.request(x, facets)
    if(length(samp.facets) != length(facets)){
      samp.facets <- c(samp.facets, ".base")
    }
    
    # get prox
    prox <- y$prox[y$prox$sample.facet %in% samp.facets,]
    
    # get matrices
    matrices <- y$LD_matrices[[which(names(y$LD_matrices) %in% facets)]]
    
    return(list(prox = prox, matrices = matrices))
  }
  
  extract.gd.afm <- function(y, facets) y[which(names(y) %in% facets)]
  
  extract.fst.matrix <- function(x, facets = NULL){
    fst <- data.table::as.data.table(.get.snpR.stats(x, facets, "weighted.means"), col_pattern = statistic_index$fst$col_pattern)
    fst[,c("p1", "p2") := tstrsplit(subfacet, "~")]
    
    cats <- unique(fst$facet)
    res <- vector("list", length(cats))
    names(res) <- cats
    
    for(i in 1:length(cats)){
      if("weighted_mean_fst_p" %in% colnames(fst)){
        res[[i]] <- list(fst = data.table::dcast(fst[facet == cats[i],], p1~p2, value.var = "weighted_mean_fst"),
                         p = data.table::dcast(fst[facet == cats[i],], p1~p2, value.var = "weighted_mean_fst_p"))
      }
      else{
        res[[i]] <- data.table::dcast(fst[facet == cats[i],], p1~p2, value.var = "weighted_mean_fst")
      }
    }
    
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
  
  facets <- check.snpR.facet.request(x, facets, "none")
  
  #========extract data======
  if(type == "single"){
    facets <- check.snpR.facet.request(x, facets)
    return(extract.basic(x@stats, facets, col_pattern = col_pattern))
  }
  else if(type == "pairwise"){
    facets <- check.snpR.facet.request(x, facets)
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
    return(extract.gd.afm(x@other$geo_dist, check.snpR.facet.request(x, facets)))
  }
  else if(type == "ibd"){
    return(extract.gd.afm(x@other$ibd, facets))
  }
  else if(type == "sample"){
    facets <- check.snpR.facet.request(x, facets, "sample")
    return(extract.basic(x@sample.stats, facets, col_pattern = col_pattern, type = "comingled"))
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
  
}
