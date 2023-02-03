#' @include snpRdata.R
NULL

#' Display snpRdata objects
#' 
#' @aliases show,snpRdata-method
#' 
#' @docType methods
#' 
#' @param object A \code{snpRdata} object.
#' @exportMethod show
#' @importFrom methods show
#' @export
setMethod("show", "snpRdata", function(object) {
  
  calced_stats_print <- character(0)
  if(length(object@calced_stats) > 0){
    for(i in 1:length(object@calced_stats)){
      calced_stats_print <- c(calced_stats_print, "Facet: ", names(object@calced_stats)[i], "\n", paste0(object@calced_stats[[i]], collapse = ", "), "\n\n")
    }
  }
  
  
  mafs <- object@stats[which(object@stats$facet == ".base"),]$maf
  
  bi_allelic <- .is.bi_allelic(object)
  
  cat(ifelse(bi_allelic, "bi-allelic", "non bi-allelic"), methods::is(object)[1], "with", nrow(object), "SNPs and", ncol(object), "samples.\n",
      "==============================================\n",
      "Average minor allele frequency:", mean(mafs), "\n",
      "Minimum minor allele frequency:", min(mafs), "\n",
      "Percent missing data:", mean(1 - (rowSums(object@geno.tables$gs[which(object@facet.meta$facet == ".base"),, drop = F])/nrow(object@sample.meta))), "\n",
      "==============================================\n",
      "Possible sample metadata facets:\n", paste0(colnames(object@sample.meta), collapse = "\t"), "\n\n",
      "Possible SNP metadata facets:\n", paste0(colnames(object@snp.meta), collapse = "\t"), "\n\n",
      "==============================================\n",
      "Currently tabulated facets:\n",
      paste0(object@facets, collapse = "\t"), "\n\n",
      "Currently calculated statistics:\n", calced_stats_print,
      "Calculated statistics can be accessed via get.snpR.stats()\n",
      "==============================================\n"
  )
})


#' Get the dimensions of a snpRdata object
#' 
#' @param x snpRdata object
#' 
#' @aliases nsnps nsamps nrow ncol dim
#' 
#' @name snpRdata_dims
NULL


setGeneric("nsnps", function(x) standardGeneric("nsnps"))
#' @export
#' @describeIn snpRdata_dims get the number of SNPs
setMethod("nsnps", "snpRdata", function(x) nrow(x))


setGeneric("nrow", function(x) standardGeneric("nrow"))
#' @export
#' @describeIn snpRdata_dims identical to nsnps
setMethod("nrow", "snpRdata", function(x) {
  nrow(snp.meta(x))
})


setGeneric("nsamps", function(x) standardGeneric("nsamps"))
#' @export
#' @describeIn snpRdata_dims get the number of samples
setMethod("nsamps", "snpRdata", function(x) ncol(x))


setGeneric("ncol", function(x) standardGeneric("ncol"))
#' @export
#' @describeIn snpRdata_dims identical to nsamps
setMethod("ncol", "snpRdata", function(x) {
  nrow(sample.meta(x))
})



#' @export
#' @describeIn snpRdata_dims get the number of SNPs and samples
setMethod("dim", "snpRdata", function(x) {
  c(nrow(x), ncol(x))
})





#' Get from or overwrite components of a snpRdata object
#'
#' Fetch or overwrite the major parts of a snpRdata object (genotypes, snp meta,
#' or sample meta). If overwritten, any calculated stats will be removed, since
#' their values may be dependant upon changes in metadata.
#'
#' @param x snpRdata object to get genotype data from.
#' @param value Genotypes, snp metadata, or sample metadata
#'
#' @name extract_snpRdata
#' @aliases genotypes genotypes<- snp.meta snp.meta<- sample.meta sample.meta<-
#'
#' @examples
#' # copy test data
#' test <- stickSNPs
#'
#' # show genotypes
#' genotypes(test)
#'
#' # show or overwrite snp meta
#' snp.meta(test)
#' snp.meta(test) <- data.frame(pos = sample(10000, nrow(test), replace = TRUE),
#'                              chr = sample(LETTERS[1:4], nrow(test), 
#'                                           replace = TRUE))
#'
#' #show or overwrite sample meta
#' sample.meta(test)
#' sample.meta(test) <- data.frame(fam = sample(LETTERS[1:4], ncol(test), 
#'                                              replace = TRUE), 
#'                                 pop = sample(LETTERS[5:8], ncol(test), 
#'                                              replace = TRUE))
NULL


setGeneric("genotypes", function(x) standardGeneric("genotypes"))
#' @export
#' @docType methods
#' @describeIn extract_snpRdata view genotypes
setMethod("genotypes", "snpRdata", function(x){
  genos <- as.data.frame(x@.Data, stringsAsFactors = FALSE)
  colnames(genos) <- x@names
  rownames(genos) <- x@row.names
  return(genos)
})


setGeneric("genotypes<-", function(x, value) standardGeneric("genotypes<-"))
#' @export
#' @describeIn extract_snpRdata set genotypes
setMethod("genotypes<-", "snpRdata", function(x, value) import.snpR.data(value, snp.meta(x), sample.meta(x), mDat = x@mDat))


setGeneric("snp.meta", function(x, value) standardGeneric("snp.meta"))
#' @export
#' @describeIn extract_snpRdata view snp meta
setMethod("snp.meta", "snpRdata", function(x, value) x@snp.meta)


setGeneric("snp.meta<-", function(x, value) standardGeneric("snp.meta<-"))
#' @export
#' @describeIn extract_snpRdata set snp meta
setMethod("snp.meta<-", "snpRdata", function(x, value) import.snpR.data(genotypes(x), value, sample.meta(x), mDat = x@mDat))


setGeneric("sample.meta", function(x, value) standardGeneric("sample.meta"))
#' @export
#' @describeIn extract_snpRdata view sample meta
setMethod("sample.meta", "snpRdata", function(x, value) x@sample.meta)


setGeneric("sample.meta<-", function(x, value) standardGeneric("sample.meta<-"))
#' @export
#' @describeIn extract_snpRdata set sample meta
setMethod("sample.meta<-", "snpRdata", function(x, value) import.snpR.data(genotypes(x), snp.meta(x), value, mDat = x@mDat))

#' @export
#' @describeIn subset_snpRdata extraction operator
#' @aliases [,snpRdata-method
#' @docType methods
setMethod("[", c("snpRdata", "ANY", "ANY", "ANY"), function(x, i, j, ..., drop = FALSE){
  if(rlang::is_missing(i)){
    i <- 1:nsnps(x)
  }
  if(rlang::is_missing(j)){
    j <- 1:nsamps(x)
  }
  
  # expand dots: for some reason extra calls don't get passed correctly without doing this...
  extra.args <- match.call()
  extra.args <- as.list(extra.args)
  extra.args <- extra.args[which(!names(extra.args) %in% c("x", "i", "j", "drop"))][-1]
  
  # pass to subset_snpR_data and return
  return(do.call(subset_snpR_data, args = c(list(x = x, .snps = i, .samps = j), extra.args)))
})

