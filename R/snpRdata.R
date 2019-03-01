#' Storage class for snpR calculated statistics and metadata for a single facet of interest.
#'
#' \code{snpR} creates and defines the snpRdata class to store both raw genotype data, sample and locus specific metadata, useful data summaries, repeatedly internally used tables, calculated summary statistics, and smoothed statistic data for a single facet of interest (such as a single population).
#'
#' The snpRdata class is built to contain SNP genotype data for use by functions in the snpR package. It inherits from the S3 class data.frame, in which the genotypes are stored, and can be manipulated identically. It also stores sample and locus specific metadata, genomic summary information, and any results from most snpR functions. The raw data for each of these latter objects is accessable via the at operator.
#' Genotypes are stored in the "character" format, as output by format_snps(). Missing data is noted with "NN".
#'
snpRdata.facet <- setClass(Class = 'snpRdata.facet', slots = c(facet.id = "character",
                                                         type = "character",
                                                         geno_tables = "list",
                                                         ac = "data.frame",
                                                         stats = "data.frame",
                                                         window.meta = "data.frame",
                                                         window.stats = "data.frame"))


#' Storage class for snpR data and calculated statistics.
#'
#' \code{snpR} creates and defines the snpRdata class to store both raw genotype data, sample and locus specific metadata, useful data summaries, repeatedly internally used tables, calculated summary statistics, and smoothed statistic data.
#'
#' The snpRdata class is built to contain SNP genotype data for use by functions in the snpR package. It inherits from the S3 class data.frame, in which the genotypes are stored, and can be manipulated identically. It also stores sample and locus specific metadata, genomic summary information, and any results from most snpR functions. The raw data for each of these latter objects is accessable via the at operator.
#' Genotypes are stored in the "character" format, as output by format_snps(). Missing data is noted with "NN".
#'
snpRdata <- setClass(Class = 'snpRdata', slots = c(sample.meta = "data.frame",
                                                   snp.meta = "data.frame",
                                                   mDat = "character",
                                                   snp.form = "numeric",
                                                   facets = "list"),
                     contains = c(data = "data.frame"))



import.snpR.data <- function(genotypes, snp.meta, sample.meta, mDat){
  # calculate some of the basic summary data, such as the genotype, allele, maf, min, ect. data

  gs <- tabulate_genotypes(genotypes, mDat = mDat, verbose = T)
  gs <- list(all = gs)
  facet.all <- snpRdata.facet(facet.id = "all", type = "both",
                              geno_tables = gs)

  x <- snpRdata(.Data = genotypes, sample.meta = sample.meta, snp.meta = snp.meta,
                facets = list(all = facet.all),
                mDat = mDat,
                snp.form = nchar(genotypes[1,1]), row.names = rownames(genotypes))
  return(x)
}


# function to add a new facet to snpRdata, generating gs, as, and wmat tables.
add.facet.snpR.data <- function(x, facets, set = "sample"){
  #=========================sanity check==========
  if(set == "sample"){
    if(!all(facets %in% colnames(x@sample.meta))){stop("Facets not found in sample metadata.\n")}
  }
  else if(set == "snp"){
    if(!all(facets %in% colnames(x@snp.meta))){stop("Facets not found in snp metadata.\n")}
  }
  else if(set == "both"){
    if(all(facets %in% colnames(x@snp.meta)) | all(facets %in% colnames(x@sample.meta)) |
       !all(facets %in% c(colnames(x@snp.meta), colnames(x@sample.meta)))){
      stop("All facets must be in either snp or sample meta data and must be spread between both if set == both.\n")
    }
  }
  else{stop("Unaccepted value for set. Options: sample, snp, or both.\n")}

  #=========================figure out unique levels for the facet==========
  if(set == "snp" | set == "both"){
    snp.meta <- x@snp.meta[colnames(x@snp.meta) %in% facets]
    snp.opts <- unique(snp.meta)
    if(!is.data.frame(snp.opts)){
      snp.opts <- as.data.frame(snp.opts, stringsAsFactors = F)
      colnames(snp.opts) <- facets[which(facets %in% colnames(x@snp.meta))]
    }
  }
  if(set == "sample" | set == "both"){
    sample.meta <- x@sample.meta[colnames(x@sample.meta) %in% facets]
    sample.opts <- unique(sample.meta)
    if(!is.data.frame(sample.opts)){
      sample.opts <- as.data.frame(sample.opts, stringsAsFactors = F)
      colnames(sample.opts) <- facets[which(facets %in% colnames(x@sample.meta))]
    }
  }

  #=========================get gs matrices==========
  gs.list <- list()
  if(set == "snp" | set == "both"){
    for(i in 1:nrow(snp.opts)){
      matches <- which(apply(snp.meta, 1, function(x) identical(as.character(x), as.character(snp.opts[i,]))))
      t.x <- x[matches,]
      if(set == "both"){
        for(j in 1:nrow(sample.opts)){
          matches <- which(apply(sample.meta, 1, function(x) identical(as.character(x), as.character(sample.opts[j,]))))
          t.x.2 <- t.x[,matches]
          gs.list[[length(gs.list) + 1]] <- tabulate_genotypes(t.x.2, x@mDat)
          names(gs.list)[length(gs.list)] <- paste0(paste0(snp.opts[i,], collapse = "."),
                                                    "_",
                                                    paste0(sample.opts[j,], collapse = "."))
        }
      }
      else{
        gs.list[[length(gs.list) + 1]] <- tabulate_genotypes(t.x, x@mDat)
        names(gs.list)[length(gs.list)] <- paste0(snp.opts[i,], collapse = ".")
      }
    }
  }
  else{
    for(i in 1:nrow(sample.opts)){
      matches <- which(apply(sample.meta, 1, function(x) identical(as.character(x), as.character(sample.opts[j,]))))
      t.x <- x[,matches]
      gs.list[[length(gs.list) + 1]] <- tabulate_genotypes(t.x, x@mDat)
      names(gs.list)[length(gs.list)] <- paste0(sample.opts[i,], collapse = ".")

    }
  }

  #=========================pack and return==========
  new.facet <- snpRdata.facet(facet.id = paste0(facets, collapse = "."),
                              type = set,
                              geno_tables = gs.list)
  x@facets <- c(x@facets, new.facet)
  names(x@facets)[length(x@facets)] <- paste0(facets, collapse = ".")

  return(x)
}
