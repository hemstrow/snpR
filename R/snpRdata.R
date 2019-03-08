#' Storage class for snpR calculated statistics and metadata for a single facet of interest.
#'
#' \code{snpR} creates and defines the snpRdata class to store both raw genotype data, sample and locus specific metadata, useful data summaries, repeatedly internally used tables, calculated summary statistics, and smoothed statistic data for a single facet of interest (such as a single population).
#'
#' The snpRdata class is built to contain SNP genotype data for use by functions in the snpR package. It inherits from the S3 class data.frame, in which the genotypes are stored, and can be manipulated identically. It also stores sample and locus specific metadata, genomic summary information, and any results from most snpR functions. The raw data for each of these latter objects is accessable via the at operator.
#' Genotypes are stored in the "character" format, as output by format_snps(). Missing data is noted with "NN".
#'
# snpRdata.facet <- setClass(Class = 'snpRdata.facet', slots = c(facet.id = "character",
#                                                          sub.facets = "character",
#                                                          type = "character",
#                                                          stats = "data.frame",
#                                                          window.meta = "data.frame",
#                                                          window.stats = "data.frame"))


#' Storage class for snpR data and calculated statistics.
#'
#' \code{snpR} creates and defines the snpRdata class to store both raw genotype data, sample and locus specific metadata, useful data summaries, repeatedly internally used tables, calculated summary statistics, and smoothed statistic data.
#'
#' The snpRdata class is built to contain SNP genotype data for use by functions in the snpR package. It inherits from the S3 class data.frame, in which the genotypes are stored, and can be manipulated identically. It also stores sample and locus specific metadata, genomic summary information, and any results from most snpR functions. The raw data for each of these latter objects is accessable via the at operator.
#' Genotypes are stored in the "character" format, as output by format_snps(). Missing data is noted with "NN".
#'
snpRdata <- setClass(Class = 'snpRdata', slots = c(sample.meta = "data.frame",
                                                   snp.meta = "data.frame",
                                                   facet.meta = "data.frame",
                                                   mDat = "character",
                                                   snp.form = "numeric",
                                                   geno.tables = "list",
                                                   ac = "data.frame",
                                                   facets = "character",
                                                   stats = "data.frame",
                                                   window.stats = "data.frame"),
                     contains = c(data = "data.frame"))



import.snpR.data <- function(genotypes, snp.meta, sample.meta, mDat){
  # calculate some of the basic summary data, such as the genotype, allele, maf, min, ect. data

  snp.meta <- cbind(snp.meta, .snp.id = 1:nrow(snp.meta))
  sample.meta <- cbind(sample.meta, .sample.id = 1:nrow(sample.meta))
  gs <- tabulate_genotypes(genotypes, mDat = mDat, verbose = T)
  ac <- format_snps(cbind(rep("temp", nrow(genotypes)), rep("temp", nrow(genotypes)), genotypes), # fix this when format_snps is updated!
                    ecs = 2, output = "ac") # add an in tab option later once format_snps is fixed.
  ac <- ac[,-c(1:2)] # remove this later

  fm <- data.frame(facet = rep("all", nrow(gs$gs)),
                   subfacet = rep("all", nrow(gs$gs)),
                   facet.type = rep("all", nrow(gs$gs)))
  fm <- cbind(fm, snp.meta)

  x <- snpRdata(.Data = genotypes, sample.meta = sample.meta, snp.meta = snp.meta,
                facet.meta = fm,
                geno.tables = gs,
                mDat = mDat,
                ac = ac,
                stats = fm,
                snp.form = nchar(genotypes[1,1]), row.names = rownames(genotypes))
  return(x)
}


# function to add a new facet to snpRdata, generating gs, as, and wmat tables, and ac formatted data.
# need to add the ac part.
add.facets.snpR.data <- function(x, facet.list){
  #=========================sanity check==========
  # check that the facet list is a list
  if(!is.list(facet.list)){
    if(is.character(facet.list)){
      warning("facet.list is not a list, and so will be treated as a single facet (for example, c(chromosome, population) would be treated as a single facet.\n")
      facet.list <- list(facet.list)
    }
    else{
      stop("Invalid facet.list.\n")
    }
  }

  # check that we haven't done these facets before, and remove any that we have.
  all.facets <- character(length = length(facet.list))
  for(i in 1:length(facet.list)){
    all.facets[i] <- paste0(facet.list[[i]], collapse = ".")
  }
  if(all(all.facets %in% x@facets)){stop("All facets already added to dataset.\n")}
  else if (any(all.facets %in% x@facets)){
    already.added <- all.facets[which(all.facets %in% x@facets)]
    warning(paste0("Some facets already present in data:\n\t", paste0(already.added, collapse = "\n\t"), "\n"))
    facet.list <- facet.list[-which(all.facets %in% x@facets)]
  }

  # check that all of the facet columns actually exist
  all.facets <- unlist(unlist(facet.list))
  opts <- c(colnames(x@sample.meta), colnames(x@snp.meta))
  used <- opts[which(opts %in% all.facets)]
  if(!all(all.facets %in% opts)){
    bad.facets <- which(!(all.facets %in% c(colnames(x@sample.meta), colnames(x@snp.meta))))
    stop(paste0("Facets ", paste0(bad.facets, collaps = " "), " not found in sample or snp metadata.\n"))
  }
  if(any(duplicated(used))){
    stop(paste0("Facets ", paste0(used[which(duplicated(used))], collapse = " "), " are duplicated in the sample and or snp metadata.\n"))
  }

  # process each facet.
  for(k in 1:length(facet.list)){
    facets <- facet.list[[k]] # column levels for this facet.
    #=========================figure out unique levels for the facet==========
    # figure out what kind of facets we are working with.
    if(any(facets %in% colnames(x@snp.meta))){
      if(any(facets %in% colnames(x@sample.meta))){
        set <- "both"
      }
      else{
        set <- "snp"
      }
    }
    else{
      set <- "sample"
    }

    # get the unique options for each facet.
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


    gs <- x@geno.tables
    #=========================get gs matrices==========
    if(set == "snp" | set == "both"){
      for(i in 1:nrow(snp.opts)){
        matches <- which(apply(snp.meta, 1, function(x) identical(as.character(x), as.character(snp.opts[i,]))))
        t.x <- x[matches,]
        if(set == "both"){
          for(j in 1:nrow(sample.opts)){
            matches2 <- which(apply(sample.meta, 1, function(x) identical(as.character(x), as.character(sample.opts[j,]))))
            t.x.2 <- t.x[,matches2]
            tgs <- tabulate_genotypes(t.x.2, x@mDat)
            gs$gs <- plyr::rbind.fill.matrix(gs$gs, tgs$gs)
            gs$as <- plyr::rbind.fill.matrix(gs$as, tgs$as)
            gs$wm <- plyr::rbind.fill.matrix(gs$wm, tgs$wm)
            x@facet.meta <- rbind(x@facet.meta,
                                  cbind(data.frame(facet = rep(paste0(facets, collapse = "."), nrow(tgs$gs)),
                                                   subfacet = rep(paste0(paste0(snp.opts[i,], collapse = "."),
                                                                         ".",
                                                                         paste0(sample.opts[j,], collapse = ".")),
                                                                  nrow(tgs$gs)),
                                                   facet.type = rep("both", nrow(tgs$gs))),
                                        x@snp.meta[matches,]))
          }
        }
        else{
          tgs <- tabulate_genotypes(t.x, x@mDat)
          gs$gs <- plyr::rbind.fill.matrix(gs$gs, tgs$gs)
          gs$as <- plyr::rbind.fill.matrix(gs$as, tgs$as)
          gs$wm <- plyr::rbind.fill.matrix(gs$wm, tgs$wm)
          x@facet.meta <- rbind(x@facet.meta,
                                cbind(data.frame(facet = rep(paste0(facets, collapse = "."), nrow(tgs$gs)),
                                                 subfacet = rep(paste0(snp.opts[i,], collapse = "."), nrow(tgs$gs)),
                                                 facet.type = rep("snp", nrow(tgs$gs))),
                                      x@snp.meta[matches,]))
        }
      }
    }
    else{
      for(i in 1:nrow(sample.opts)){
        matches <- which(apply(sample.meta, 1, function(x) identical(as.character(x), as.character(sample.opts[i,]))))
        t.x <- x[,matches]
        tgs <- tabulate_genotypes(t.x, x@mDat)
        gs$gs <- plyr::rbind.fill.matrix(gs$gs, tgs$gs)
        gs$as <- plyr::rbind.fill.matrix(gs$as, tgs$as)
        gs$wm <- plyr::rbind.fill.matrix(gs$wm, tgs$wm)
        x@facet.meta <- rbind(x@facet.meta,
                              cbind(data.frame(facet = rep(paste0(facets, collapse = "."), nrow(tgs$gs)),
                                               subfacet = rep(paste0(sample.opts[i,], collapse = "."), nrow(tgs$gs)),
                                               facet.type = rep("sample", nrow(tgs$gs))),
                                    x@snp.meta))

      }
    }

    #=========================pack and return==========
    x@geno.tables <- gs
    x@facets <- c(x@facets, paste0(facets, collapse = "."))
  }

  return(x)

}

# function to list existing facets.
find.snpR.facets <- function(x){
  facets <- vector("list", length(x@facets))
  names(facets) <- x@facets
  for(i in 1:length(facets)){
    facets[[i]] <- unique(x@facet.meta[x@facet.meta$facet == names(facets)[i],]$subfacet)
  }
  return(facets)
}

# function to pull stats for a given facet
get.snpR.stats <- function(x, facets = NULL){
  if(!is.null(facets)){
    if(facets[1] == "all"){
      facets <- dat@facets
    }
  }
  else {
    facets <- "all"
  }

  return(x@stats[which(x@stats$facet %in% facets), -which(colnames(x@stats) %in% c(".snp.id", "facet.type"))])
}

# function to apply a function across selected facets
# req: which part of the snpR.data object is required and should be pulled out? gs: genotype tables
# facets: if NULL, run on everything.
# cases: ps, per snp.
# fun: which function should be applied?
apply.snpR.facets <- function(x, facets = NULL, req, fun, case = "ps", ...){
  if(!is.null(facets)){
    if(facets[1] == "all"){
      facets <- dat@facets
    }
  }
  else {
    facets <- "all"
  }
  if(case == "ps"){
    if(req == "gs"){
      # bind metadata
      gs <- x@geno.tables
      gs <- plyr::llply(gs, function(y) cbind(x@facet.meta, y))

      # subset
      gs <- plyr::llply(gs, function(y) subset(y, y$facet %in% facets))

      # remove metadata
      gs <- plyr::llply(gs, function(y) dplyr::select(y, -c(facet, subfacet, facet.type, colnames(x@snp.meta))))

      # convert back to matrix
      gs <- plyr::llply(gs, function(y) as.matrix(y))

      # run the function indicated
      out <- fun(gs)

      # bind metadata
      out <- cbind(x@facet.meta[x@facet.meta$facet %in% facets,], out)

      # return
      return(out)
    }
  }
}

