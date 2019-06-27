check.snpRdata <- function(object){
  errors <- character()

  # check that there are the correct number of snp and sample metadata rows
  if(nrow(object@snp.meta) != nrow(object)){
    msg <- paste0("Number of snps in metadata (", nrow(object@snp.meta),
                  ") is not equal to the number of snps in data (", nrow(object), ").")
    errors <- c(errors, msg)
  }
  if(nrow(object@sample.meta) != ncol(object)){
    msg <- paste0("Number of samples in metadata (", nrow(object@sample.meta),
                  ") is not equal to the number of samples in data (", ncol(object), ").")
    errors <- c(errors, msg)
  }

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
#'in the snpR package. It inherits from the S3 class data.frame, in which the
#'genotypes are stored, and can be manipulated identically. It also stores
#'sample and locus specific metadata, genomic summary information, and any
#'results from most snpR functions. The raw data for each of these latter
#'objects is accessable via the at operator. Genotypes are stored in the
#'"character" format, as output by format_snps(). Missing data is noted with
#'"NN".
#'
#'For more information, see \code{\link{import.snpR.data}}, the constructor
#'function for this object class.
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
                                       pairwise.LD = "list",
                                       window.bootstraps = "data.frame"),
         contains = c(data = "data.frame"),
         validity = check.snpRdata)

