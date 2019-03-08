#' Import genotype and metadata into a snpRdata object.
#'
#' \code{import.snpR.data} converts genotype and meta data to the snpRdata class, which stores raw genotype data, sample and locus specific metadata, useful data summaries, repeatedly internally used tables, calculated summary statistics, and smoothed statistic data.
#'
#' The snpRdata class is built to contain SNP genotype data for use by functions in the snpR package. It inherits from the S3 class data.frame, in which the genotypes are stored, and can be manipulated identically. It also stores sample and locus specific metadata, genomic summary information, and any results from most snpR functions. The raw data for each of these latter objects is accessable via the at operator.
#' Genotypes are stored in the "character" format, as output by format_snps(). Missing data is noted with "NN".
#'
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

  x <- new("snpRdata", .Data = genotypes, sample.meta = sample.meta, snp.meta = snp.meta,
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
  #===========================sanity check==========
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
    stop(paste0("Facets ", paste0(all.facets[bad.facets], collapse = " "), " not found in sample or snp metadata.\n"))
  }
  if(any(duplicated(used))){
    stop(paste0("Facets ", paste0(used[which(duplicated(used))], collapse = " "), " are duplicated in the sample and or snp metadata.\n"))
  }

  #===========================process each facet.===================
  for(k in 1:length(facet.list)){
    facets <- facet.list[[k]] # column levels for this facet.
    #=========================figure out unique levels for the facet==========
    # figure out what kind of facets we are working with.
    if(any(facets %in% colnames(x@snp.meta))){
      if(any(facets %in% colnames(x@sample.meta))){
        set <- "both"
        x@facet.type <- c(x@facet.type, "both")
      }
      else{
        set <- "snp"
        x@facet.type <- c(x@facet.type, "snp")
      }
    }
    else{
      set <- "sample"
      x@facet.type <- c(x@facet.type, "sample")
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

  # add this when I get format_snps up and running.
  # browser()
  # added.facets <- lapply(facet.list, function(x) paste0(x, collapse = "."))
  # added.facets <- unlist(added.facets)
  # tas <- x@geno.tables$as[x@facet.meta$facet %in% added.facets,]
  # ni1 <- matrixStats::rowMaxs(tas)
  # nt <- rowSums(tas)
  # x@ac <- rbind(x@ac, data.frame(n_total = nt,
  #                                n_alleles = rowSums(matrix(as.logical(tas), nrow = nrow(tas))),
  #                                ni1 = ni1,
  #                                ni2 = nt - ni1))
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

  # add any missing facets
  miss.facets <- facets[which(!(facets %in% x@facets))]
  if(length(miss.facets) != 0){
    # need to fix any multivariate facets (those with a .)
    comp.facets <- grep("\\.", miss.facets)
    run.facets <- as.list(miss.facets[-c(comp.facets)])
    run.facets <- c(run.facets, strsplit(miss.facets[comp.facets], split = "\\."))
    x <- add.facets.snpR.data(x, as.list(run.facets))
  }

  if(case == "ps"){
    if(req == "gs"){
      # bind metadata
      gs <- x@geno.tables
      gs <- plyr::llply(gs, function(y) cbind(x@facet.meta, y))

      # subset
      gs <- plyr::llply(gs, function(y) subset(y, y$facet %in% facets))

      # remove metadata
      gs$gs <- gs$gs[,-which(colnames(gs$gs) %in% c("facet", "subfacet", "facet.type", colnames(x@snp.meta)))]
      gs$as <- gs$as[,-which(colnames(gs$as) %in% c("facet", "subfacet", "facet.type", colnames(x@snp.meta)))]
      gs$wm <- gs$wm[,-which(colnames(gs$wm) %in% c("facet", "subfacet", "facet.type", colnames(x@snp.meta)))]

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




#' Tabulate allele and genotype counts at each locus.
#'
#' \code{tabulate_genotypes} creates matricies containing counts of observed alleles and genotypes at each locus.
#'
#' This function is pirmarily used interally in several other funcitons, but may occasionally be useful.
#'
#' @param x Input genotype data, where columns are individuals and rows are snps. No metadata.
#' @param mDat Character string. How are missing \emph{genotypes} noted?
#' @param verbose Logical. Should the function report progress?
#'
#' @return A list of matrices. gs is the genotype matrix, as is the allele matrix, and wm is the genotype matrix with missing genotypes.
#'
#' @examples
#' tabulate_genotypes(stickSNPs[,-c(1:3)], "NN")
#'
tabulate_genotypes <- function(x, mDat, verbose = F){

  # get a genotype table
  snp_form <- nchar(x[1,1])   # get information on data format
  x <- reshape2::melt(t(x)) # transpose and melt
  gmat <-table(x[,2:3]) # table to get a count of genotypes per locus
  tmat <- gmat[,-which(colnames(gmat) == mDat)] # remove missing data

  #get matrix of allele counts
  #initialize
  hs <- substr(colnames(tmat),1,snp_form/2) != substr(colnames(tmat), (snp_form/2 + 1), snp_form*2) # identify heterozygotes.
  if(verbose){cat("Getting allele table...\n")}
  as <- unique(unlist(strsplit(paste0(colnames(tmat)), "")))
  amat <- matrix(0, nrow(gmat), length(as))
  colnames(amat) <- as

  #fill in
  for(i in 1:length(as)){
    b <- grep(as[i], colnames(tmat))
    hom <- which(colnames(tmat) == paste0(as[i], as[i]))
    if(length(hom) == 0){
      het <- b
      amat[,i] <- rowSums(tmat[,het])
    }
    else{
      het <- b[b != hom]
      if(length(het) > 0){
        if(is.matrix(tmat[,het])){
          amat[,i] <- (tmat[,hom] * 2) + rowSums(tmat[,het])
        }
        else{
          amat[,i] <- (tmat[,hom] * 2) + tmat[,het]
        }
      }
      else{
        amat[,i] <- (tmat[,hom] * 2)
      }
    }
  }
  return(list(gs = unclass(tmat), as = amat, wm = unclass(gmat)))
}


#'Filter SNP data.
#'
#'\code{filter_snps} filters SNP data to remove loci which violate any of several assumptions and/or individuals which are sequenced at too few SNP loci.
#'
#'Description of x:
#'    Contains metadata in columns 1:ecs. Remainder of columns contain genotype calls for each individual. Each row is a different SNP, as given by format_snps output options 4 or 6.
#'
#'Possible filters:
#'\itemize{
#'    \item{maf, minor allele frequency: }{removes SNPs where the minor allele frequency is too low. Can look for mafs below #'provided either globally or search each population individually.}
#'    \item{hf_hets, high observed heterozygosity: }{removes SNPs where the observed heterozygosity is too high.}
#'    \item{min_ind, minimum individuals: }{removes SNPs that were genotyped in too few individuals.}
#'    \item{min_loci, minimum loci: }{removes individuals sequenced at too few loci.}
#'    \item{non_poly, non-polymorphic SNPs: }{removes SNPs that are not polymorphic (not true SNPs).}
#'    \item{bi_al, non-biallelic SNPs: }{removes SNPs that have more than two observed alleles.}
#'}
#'
#'Note that filtering out poorly sequenced individuals creates a possible conflict with the loci filters, since after individuals are removed, some loci may no longer pass filters. For example, if a portion of individuals in one population all carry the only instances of a rare minor allele that still passes the maf threshold, removing those individuals may cause the loci to no longer be polymorphic in the sample.
#'
#'To counter this, the "re_run" argument can be used to pass the data through a second filtering step after individuals are removed. By default, the "partial" re-run option is used, which re-runs only the non-polymorphic filter (if it was originally set), since these may cause downstream analysis errors. The "full" option re-runs all set filters. Note that re-running any of these filters may cause individuals to fail the individual filter after loci removal, and so subsequent tertiary re-running of the individual filters, followed by the loci filters, and so on, could be justified. This function stops after the second loci re-filtering, since that step is likely to be the most important to prevent downstream analytical errors.
#'
#'Via the "pop" argument, this function can filter by minor allele frequencies in either \emph{all} samples or \emph{in each population and the entire sample}. The latter should be used in instances where populaiton sizes are very different or there are \emph{many} populations, and thus common alleles of interest in one population might be otherwise filtered out. With very small populations, however, this may leave noise in the sample! In most cases, filtering the entire sample is sufficient.
#'
#' @param x data.frame. Input data, in the numeric or character format as given by format_snps options 4 or 6.
#' @param ecs Integer. Number of metadata columns at the start of x.
#' @param maf FALSE or numeric between 0 and 1, default FALSE. Minimum acceptable minor allele frequency
#' @param hf_hets FALSE or numeric between 0 and 1, default FALSE. Maximum acceptable heterozygote frequency.
#' @param min_ind FALSE or integer, default FALSE. Minimum number of individuals in which a loci must be sequenced.
#' @param min_loci FALSE or numeric between 0 and 1, default FALSE. Minimum proportion of SNPs an individual must be genotyped at.
#' @param re_run FALSE, "partial", or "full", default "partial". How should loci be re_filtered after individuals are filtered?
#' @param pop FALSE or table, default FALSE. A table with population information for individuals. Individuals must be sorted in input data in the population order given in this table.
#' @param non_poly boolean, default TRUE. Should non-polymorphic loci be removed?
#' @param bi_al boolean, default TRUE. Should non-biallelic SNPs be removed?
#' @param mDat character variable, default "NN". Format of missing \emph{genotypes}. Overall data format is infered from this. Can be either "NN" or "0000".
#' @param in.tab. FALSE or list. Option to provide tables of snp and genotype counts at each loci, used in many reformatting and filtering steps. Used internally.
#' @param out.tab. FALSE or list. Option to return tables of snp and genotype counts at each loci, used in many reformatting and filtering steps. Used internally.
#'
#' @return A data.frame in the same format as the input, with SNPs and individuals not passing the filters removed.
#'
#' @examples
#' #Strict filtering for missing individuals and unsequenced loci with partial re-run:
#' filter_snps(stickSNPs, 3, 0.05, 0.55, 250, .75)
#'
#' #Strict maf filtering with pops.
#' ##prep pop info
#' l <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' ##filter
#' filter_snps(stickSNPs, 3, 0.05, 0.55, 250, .75, "full", pop = l)
#'
filter_snps <- function(x, maf = FALSE, hf_hets = FALSE, min_ind = FALSE,
                        min_loci = FALSE, re_run = "partial", maf.facets = NULL,
                        non_poly = TRUE, bi_al = TRUE){

  #==============do sanity checks====================

  if(maf){
    if(!is.numeric(maf)){
      stop("maf must be a numeric value.")
    }
    if(length(maf) != 1){
      stop("maf must be a numeric vector of length 1.")
    }
  }

  if(hf_hets){
    if(!is.numeric(hf_hets)){
      stop("hf_hets must be a numeric value.")
    }
    if(length(hf_hets) != 1){
      stop("hf_hets must be a numeric vector of length 1.")
    }
  }

  if(min_ind){
    if(!is.numeric(min_ind)){
      stop("min_ind must be a numeric value.")
    }
    if(length(min_ind) != 1){
      stop("min_ind must be a numeric vector of length 1.")
    }
    if(min_ind > 1 | min_ind < 0){
      stop("min_ind is the minimum proportion of sequenced individuals, and so must be between 0 and 1.\n")
    }
  }

  if(min_loci){
    if(!is.numeric(min_loci) | (min_loci <= 0 | min_loci >= 1) | length(min_loci) != 1){
      stop("min_loci must be a numeric value between but not equal to 0 and 1.")
    }
  }

  if(re_run != FALSE){
    if(re_run != "partial" & re_run != "full"){
      cat("re_run must be set to partial or full if not FALSE.\n")
    }
  }

  if(!is.null(maf.facets[1])){
    # add any needed facets...
    miss.facets <- maf.facets[which(!(maf.facets %in% x@facets))]
    if(length(miss.facets) != 0){
      cat("Adding missing facets...\n")
      # need to fix any multivariate facets (those with a .)
      comp.facets <- grep("\\.", miss.facets)
      run.facets <- as.list(miss.facets[-c(comp.facets)])
      run.facets <- c(run.facets, strsplit(miss.facets[comp.facets], split = "\\."))
      x <- add.facets.snpR.data(x, as.list(run.facets))
    }

    # check for bad facets to remove (those that don't just consider samples)
    if(any(x@facet.type[x@facets %in% maf.facets] != "sample")){
      vio.facets <- maf.facets[match(maf.facets, x@facets)]
      vio.facets <- vio.facets[which(x@facet.type[x@facets %in% maf.facets] != "sample")]
      warning(paste0("Facets over which to maf.filter must be sample specific facets, not snp specific facets! Removing non-sample facets: \n", paste0(vio.facets, collapse = " "), ".\n"))
      maf.facets <- maf.facets[-which(maf.facets %in% vio.facets)]
    }
  }

  #==============set up, get values used later, clean up data a bit,define subfunctions==========
  cat("Initializing...\n")

  #get headers
  headers <- x@snp.meta
  snp_form <- x@snp.form
  mDat <- x@mDat


  #function to filter by loci, to be called before and after min ind filtering (if that is requested.)
  filt_by_loci <- function(){
    #==========================run filters========================
    vio.snps <- logical(nrow(x)) #vector to track status

    amat <- x@geno.tables$as[x@facet.meta$facet == "all",]
    gmat <- x@geno.tables$gs[x@facet.meta$facet == "all",]
    wmat <- x@geno.tables$wm[x@facet.meta$facet == "all",]

    # non-biallelic and non-polymorphic loci
    if(bi_al | non_poly){
      bimat <- ifelse(amat, TRUE, FALSE)

      if(bi_al){
        cat("Filtering non-biallelic loci...\n")
        bi <- ifelse(rowSums(bimat) > 2, T, F) # if true, should keep the allele
        bi.vio <- x@snp.meta$.snp.id[which(bi)] # IDs of violating snps.
        cat(paste0("\t", length(bi.vio), " bad loci\n"))
        vio.snps[which(bi)] <- T
      }

      if(non_poly){
        cat("Filtering non_polymorphic loci...\n")
        np <- ifelse(rowSums(bimat) < 2, T, F) # if true, should keep the allele
        np.vio <- x@snp.meta$.snp.id[which(np)]
        cat(paste0("\t", length(np.vio), " bad loci\n"))
        vio.snps[which(np)] <- T
      }

      #some tests require this, so subset the matrices and redefine things if true and some are multi-allelic
      # if((maf | hf_hets) &  sum(bi) != 0){
      #   x <- x[bi,]
      #   xv <-   xv <- as.vector(t(x))
      #   gs <- unique(xv)
      #   hs <- substr(gs,1,snp_form/2) != substr(gs, (snp_form/2 + 1), snp_form*2)
      #   mpos <- which(gs == mDat)
      #   as <- unique(c(substr(gs,1,snp_form/2), substr(gs, (snp_form/2 + 1), snp_form*2)))
      #   as <- as[as != substr(mDat, 1, snp_form/2)] #that aren't N?
      #   gmat <- gmat[bi,]
      #   tmat <- tmat[bi,]
      #   amat <- amat[bi,]
      #   headers <- headers[bi,]
      # }
      # else{
      #   keep <- keep + bi
      # }
    }

    #========min inds=======
    if(min_ind){
      cat("Filtering loci sequenced in few individuals...\n")
      mi <- wmat[,colnames(wmat) == mDat]
      mi <- (nrow(x@sample.meta) - mi)/nrow(x@sample.meta) < min_ind
      mi.vio <- x@snp.meta$.snp.id[which(mi)]
      cat(paste0("\t", length(mi.vio), " bad loci\n"))
    }

    #========minor allele frequency, both total and by pop. Should only run if bi_al = TRUE.=========
    if(maf){
      #if not filtering with multiple pops
      if(is.null(maf.facets)){
        cat("Filtering low minor allele frequencies, no pops...\n")

        # check to see if we need to calculate mafs:
        if(any(colnames(x@stats) == "maf")){ # check that mafs have been calculated, the all facet must exist
          if(any(is.na(x@stats$maf[x@stats$facet == "all"]))){ # check that mafs have been calculated for the all facet
            mafs <- 1 - matrixStats::rowMaxs(amat)/rowSums(amat)
          }
          else{
            mafs <- x@stats$maf[x@stats$facet == "all"]
          }
        }
        else{
          mafs <- 1 - matrixStats::rowMaxs(amat)/rowSums(amat)
        }


        mafs <- mafs < maf #less than required, set to true and reject.
        mafs[is.na(mafs)] <- TRUE
        maf.vio <- x@snp.meta$.snp.id[which(mafs)]
        vio.snps[which(mafs)] <- T
      }
      else{
        cat("Filtering low minor allele frequencies by facet.\n")
        # pmafs <- logical(nrow(x))

        # see if we need to calculate mafs
        if(any(colnames(x@stats) == "maf")){ # mafs have been caluclated
          # get mafs for any uncalculated facets

          # if mafs have been calculated, but not for of the requested any facets...
          if(!any(x@stats$facet %in% maf.facets)){
            x <- calc_maf(x, facets = maf.facets)
          }
          #if mafs have been calculated for some but not all of our facets
          else if(any(is.na(x@stats$maf[x@stats$facet %in% maf.facets]))){
            run.facets <- unique(x@stats$facet[which(is.na(x@stats$maf[x@stats$facet %in% maf.facets]))])
            x <- calc_maf(x, facets = run.facets)
          }
        }
        else{
          x <- calc_maf(x, facets = maf.facets)
        }

        # grab mafs
        mafs <- x@stats[x@stats$facet %in% maf.facets,]
        mafs$maf <- mafs$maf < maf

        # now, figure out in how many subfacets the maf is too low. If all, the loci violates the filter
        cmafs <- reshape2::dcast(mafs[,c(".snp.id", "maf")], fun.aggregate = sum,
                                 formula = ... ~ "maf", value.var = "maf")

        # add in the overall maf, since differential fixation would otherwise be removed.
        if(any(is.na(x@stats$maf[x@stats$facet == "all"]))){ # check that mafs have been calculated for the all facet
          a.mafs <- 1 - matrixStats::rowMaxs(amat)/rowSums(amat)
        }
        else{
          a.mafs <- x@stats$maf[x@stats$facet == "all"]
        }

        a.mafs <- a.mafs < maf
        cmafs$maf <- cmafs$maf + a.mafs

        # check vio and report
        maf.vio <- x@snp.meta$.snp.id[which(cmafs$maf == (1 + length(unique(mafs$subfacet))))]
        cat(paste0("\t", length(maf.vio), " bad loci\n"))
        vio.snps[which(cmafs$maf == (1 + length(unique(mafs$subfacet))))] <- T
      }
    }

    #========hf_hets. Should only run if bi_al = TRUE.==========
    # working here
    if(hf_hets){
      cat("Filtering high frequency heterozygote loci...\n")

      # get heterozygote frequency
      hs <- which(substr(colnames(gmat), 1, snp_form/2) != substr(colnames(gmat), (snp_form/2) + 1, snp_form))
      het_f <- rowSums(gmat[,hs])/rowSums(gmat)

      # check violation
      het_f <- het_f > hf_hets #if false, heterozygote frequency is lower than cut-off, keep locus
      het_f.vio <- x@snp.meta$.snp.id[which(het_f)]
      cat(paste0("\t", length(het_f.vio), " bad loci\n"))
      vio.snps[which(het_f)] <- T
    }


    #==========remove violating loci==================
    if(any(vio.snps)){
      vio.ids <- x@snp.meta$.snp.id[which(vio.snps)]
      ngs <- x@geno.tables$gs[-which(x@facet.meta$.snp.id %in% x@snp.meta$.snp.id[which(vio.snps)]),]
      nas <- x@geno.tables$as[-which(x@facet.meta$.snp.id %in% x@snp.meta$.snp.id[which(vio.snps)]),]
      nwm <- x@geno.tables$wm[-which(x@facet.meta$.snp.id %in% x@snp.meta$.snp.id[which(vio.snps)]),]
      ngs <- list(gs = ngs, as = nas, wm = nwm)
      rm(nas, nwm)

      x <- snpRdata(.Data = x[-which(vio.snps),],
                    sample.meta = x@sample.meta,
                    snp.meta = x@snp.meta[-which(vio.snps),],
                    facet.meta = x@facet.meta[-which(x@facet.meta$.snp.id %in% vio.ids),],
                    geno.tables = ngs,
                    # ac = x@ac[-which(x@facet.meta$.snp.id %in% vio.ids),], add this when I get ac up and running.
                    stats = x@stats[-which(x@stats$.snp.id %in% vio.ids)],
                    window.stats = x@window.stats,
                    facets = x@facets,
                    facet.type = x@facet.type,
                    row.names = x@row.names[-which(vio.snps)])
    }
    return(x)
  }

  #funciton to filter by individuals.
  min_loci_filt <- function(){
    cat("Filtering out individuals sequenced in few kept loci...\n")
    mcounts <- colSums(ifelse(x == mDat, 1, 0))
    rejects <- which(mcounts/nrow(x) >= (1 - min_loci))
    if(length(rejects) > 0){
      old.facets <- x@facets
      old.facets <- sapply(old.facets, function(x) strsplit(x, "\\."))
      invisible(capture.output(x <- import.snpR.data(x[,-rejects],
                                                     snp.meta = x@snp.meta,
                                                     sample.meta = x@sample.meta[-rejects,],
                                                     mDat = mDat)))
      cat("Re-calculating and adding facets.\n")
      x <- add.facets.snpR.data(x, old.facets)
      warning("Any calculated stats will be removed, since individuals were filtered out!\n")
    }
    return(list(x = x, rejects = rejects))
  }

  #==========================call the functions as requested.==================
  if(any(c(non_poly, bi_al, maf, hf_hets, min_ind) != FALSE)){
    cat("Filtering loci. Starting loci:", nrow(x), "\n")

    # run the filter
    x <- filt_by_loci()

    if(nrow(x) == 1){
      stop("No loci remain after filters.")
    }

    cat("\tEnding loci:", nrow(x), "\n")
  }

  # run the minimum sequenced loci filter
  if(min_loci){
    cat("Filtering individuals. Starting individuals:", ncol(x), "\n")
    x <- min_loci_filt()
    if(length(x$rejects) == 0){
      cat("No individuals removed.\n")
    }
    else{
      x <- x$x
      cat("\tEnding individuals:", ncol(x), "\n")
      if(re_run != FALSE){
        cat("Re-filtering loci...\n")

        if(re_run == "partial"){
          maf <- FALSE
          hf_hets <- FALSE
          min_ind <- FALSE
          bi_al <- FALSE
        }
        if(any(c(non_poly, bi_al, maf, hf_hets, min_ind) != FALSE)){
          x <- filt_by_loci() # re-filter loci to make sure that we don't have any surprise non-polys ect.
          cat("\tFinal loci count:", nrow(x), "\n")
        }
        else{
          cat("\tNo variables to re-fitler.\n")
        }
      }
    }
  }

  # return results
  return(x)
}

#'Re-format genomic data.
#'
#'\code{format_snps} reformats SNP data into a range of different possible formats for use in snpR functions and elsewhere. Supports microsatellite data for conversion to a few formats as well.
#'
#'
#'Format options:
#'\itemize{
#'    \item{ac: }{allele count format, allele counts tabulated for all samples or within populations.}
#'    \item{genepop: }{genepop format, genotypes stored as four numeric characters (e.g. "0101", "0204"), transposed, and formatted for genepop. Rownames are individual IDs in genepop format, colnames are SNP ids, matching the first metadata column in input.}
#'    \item{structure: }{STRUCTURE format, two lines per individual: allele calls stored as single character numeric (e.g. "1", "2"). Allele calls per individual stored on two subsequent lines.}
#'    \item{numeric: }{numeric genotype tab format, genotypes stored as four numeric characters (e.g. "0101", "0204").}
#'    \item{hapmap: }{Migrate-n hapmap, allele counts tabulated within populations, in migrate-n hapmap format. Since this migrate-n implementation is iffy, this probably shouldn't be used much.}
#'    \item{character: }{character genotype tab format, genotypes stored as actual base calls (e.g. "AA", "CT").}
#'    \item{pa: }{allele presence/absence format, presence or absence of each possible allele at each possible genotype noted. Interpolation possible, with missing data substituted with allele freqency in all samples or each population.}
#'    \item{rafm: }{RAFM format, two allele calls at each locus stored in subsequent columns, e.g. locus1.1 locus1.2.}
#'    \item{faststructure: }{fastSTRCTURE format, identical to STRUCTURE format save with the addition of filler columns proceeding data such that exactly 6 columns proceed data. These columns can be filled with metadata if desired.}
#'    \item{dadi: }{dadi format SNP data format, requires two columns named "ref" and "anc" with the flanking bases around the SNP, e.g. "ACT" where the middle location is the A/C snp.}
#'    \item{plink: }{PLINK! binary input format, requires columns named "group", "snp", and "position", and may contain a column named "cM", "cm", or "morgans", containing linkage group/chr, snp ID, position in bp, and distance in cM in order to create .bim extended map file.}
#'}
#'
#'Example datasets in each format are available in \code{\link{stickFORMATs}} in elements named for output options.
#'
#'Input formats: For now, all input formats require at least two metadata columns.
#'\itemize{
#'    \item{NN: }{SNP genotypes stored as actual base calls (e.g. "AA", "CT").}
#'    \item{0000: }{SNP genotypes stored as four numeric characters (e.g. "0101", "0204").}
#'    \item{msat_2: }{Microsatellite genotypes stored as four numeric characters (e.g. "0101", "2740").}
#'    \item{msat_3: }{Microsatellite genotypes stored as six numeric characters (e.g. "120123", "233235").}
#'    \item{snp_tab: }{SNP genotypes stored with genotypes in each cell, but only a single nucleotide noted if homozygote and two nucleotides seperated by a space if heterozygote (e.g. "T", "T G").}
#'    \item{0_geno: }{SNP genotypes stored with genotypes in each cell as 0 (homozyogous allele 1), 1 (heterozygous), or 2 (homozyogus allele 2).}
#'}
#'
#'Currently, msat_2 and msat_3 only support conversion to output option 7. 2, 3, and 4 are forthcoming.
#'
#'
#' @param x data.frame. Input data, in any of the above listed input formats.
#' @param ecs Integer. Number of extra metadata columns at the start of x.
#' @param output Character. Which of the output format should be used?
#' @param input_form Character, default "NN". Which of the above input formats should be used (e.g. "NN", "msat_2")?
#' @param mDat Character, default "NN". The coding for missing \emph{genotypes} in x (typically "NN" or "0000").
#' @param pop FALSE or table, default FALSE. A table with population information for individuals. Individuals must be sorted in input data in the population order given in this table.
#' @param FALSE or table, default FALSE. A table with population information for individuals. Individuals must be sorted in input data in the population order given in this table.
#' @param n_samp Integer or numeric vector, default NA. For output option 3. How many random loci should be selected? Can either be an integer or a numeric vector of loci to use.
#' @param interp_miss boolean, default TRUE. For output option 7. Should missing data be interpolated? Needed for PCA, ect.
#' @param lnames character vector, default NULL. For output option 7, optional vector of locus names by which to name output columns. If not provided, will use 1:nrow(x).
#' @param outfile character vector, default FALSE. If provided, the path to the output file to write to.
#' @param ped data.frame, default NULL. If provided, the six column header for plink .ped files.
#' @param in.tab FALSE or list. Option to provide tables of snp and genotype counts at each loci, used in many reformatting and filtering steps. Used internally.
#' @param out.tab FALSE or list. Option to return tables of snp and genotype counts at each loci, used in many reformatting and filtering steps. Used internally.
#'
#' @return A data.frame in the specified format.
#'
#' @examples
#' #allele count with pops:
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' format_snps(stickSNPs, 3, "ac", pop = pops)
#'
#' #genepop:
#' format_snps(stickSNPs, 3, "genepop")
#'
#' #STRUCTURE, subsetting out 100 random alleles:
#' format_snps(stickSNPs, 3, "structure", n_samp = 100)
#'
#' #STRUCTURE, subseting out the first 100 alleles:
#' format_snps(stickSNPs, 3, "structure", n_samp = 1:100)
#'
#' #fastSTRUCTURE
#' format_snps(stickSNPs, 3, "faststructure")
#'
#' #numeric:
#' format_snps(stickSNPs, 3, "numeric")
#'
#' #hapmap for migrate-n:
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' format_snps(stickSNPs, 3, "hapmap", pop = pops)
#'
#' #character:
#' num <- format_snps(stickSNPs, 3, 4)
#' format_snps(num, 3, "character", input_form = "0000", miss = "00")
#'
#' #presence/absence, SNP data:
#' format_snps(stickSNPs, 3, "pa")
#'
#' #presence/absence, 3 character microsat data (2 character is very similar):
#' format_snps(sthMSATS[seq(1, 13, by = 4),], 3, "pa", input_form = "msat_3", miss = "000")
#'
#' #RAFM, taking only 100 random snps.
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' format_snps(stickSNPs, 3, "rafm", pop = pops, n_samp = 100)
#'
#' #dadi
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' format_snps(cbind(ref = "ATA", anc = "ACT", stickSNPs), 5, "dadi", pop = pops)
#'
#' #PLINK! format
#' format_snps(stickSNPs, 3, "plink", outfile = "plink_out")
#' #from command line, then run plink_out.sh to generate plink_out.bed.
#'
#' #PLINK! format with provided ped
#' ped <- data.frame(fam = c(rep(1, 210), rep("FAM2", 210)), ind = 1:420, mat = 1:420, pat = 1:420, sex = sample(1:2, 420, T), pheno = sample(1:2, 420, T))
#' format_snps(stickSNPs, 3, "plink", outfile = "plink_out", ped = ped)
#' #from command line, then run plink_out.sh to generate plink_out.bed.
#'
format_snps <- function(x, ecs, output = 1, input_form = "NN",
                        mDat = "NN", pop = 1, n_samp = NA,
                        interp_miss = T, lnames = NULL, outfile = FALSE,
                        ped = NULL,
                        in.tab = FALSE, out.tab = FALSE){

  #redefine x as "data", since this was originally written a while ago and already contains a variable called "x"
  data <- x
  rm(x)
  #redefine mDat as miss, since this was written prior to standardization
  miss <- substr(mDat, 1, nchar(mDat)/2)


  #possible outputs, recode to numbers so that I don't have to rewrite the whole damn thing. Still allowing numbers for reverse compatiblity.
  if(!is.numeric(output)){
    pos_outs <- c("ac", "genepop", "structure", "numeric", "hapmap", "character", "pa",
                  "rafm", "faststructure", "dadi", "plink", "sn")
    if(!(output %in% pos_outs)){
      stop("Unaccepted output format specified. Check documentation. Remember, no capitalization!")
    }

    #reformat output to a number corresponding to format.
    output <- which(pos_outs == output)
  }


  if(ecs < 2){
    stop("Cannot have less than 2 metadata (ecs) columns.")
  }

  #####################################################################
  #do checks, print info
  if(output == 1){
    cat("Converting to per pop allele count format.\n")
    if(input_form != "0000" & input_form != "NN"){
      stop("Only 0000 and NN formats accepted.")
    }
    if(is.table(pop)){
      cat("Pops:")
      for (i in 1:length(pop)){
        cat("\n\t", names(pop)[i], ", size: ", as.numeric(pop[i]))
      }
      cat("\n")
    }
  }
  else if(output == 2){
    cat("Converting to genepop format.\n")
    if(input_form != "0000" & input_form != "NN" & input_form != "snp_tab"){
      stop("Only 0000, NN, and snp_tab formats accepted for now.")
    }
  }
  else if(output == 3 | output == 9){
    if(output == 3){cat("Converting to STRUCTURE format.\n")}
    else{cat("Converting to fastSTRUCTURE format.\n")}
    if(input_form != "0000" & input_form != "NN"){
      stop("Only 0000 and NN formats accepted for now.")
    }
    if(length(n_samp) > 1){
      if(is.integer(n_samp)){
        cat("Number of designated sub-samples to take:", length(n_samp), "\n")
      }
      else{
        stop("Number of sub-samples to take must be a positive integer vector.\n")
      }
    }
    else if (is.numeric(n_samp)){
      if(floor(n_samp) == n_samp & n_samp > 0){
        cat("Number of designated sub-samples to take:", n_samp, "\n")
      }
      else{
        stop("Number of sub-samples to take must be a positive integer vector.\n")
      }
    }
    else if (!is.na(n_samp)){
      stop("Number of sub-samples to take must be a positive integer vector.\n")
    }
  }
  else if(output == 4){
    cat("Converting to numeric 2 character format.\n")
    if(input_form != "NN" & input_form != "snp_tab"){
      stop("Only NN and snp_tab formats accepted for now.")
    }
  }
  else if(output == 5){
    cat("Converting to migrate-N hapmap format.\n")
    if(input_form != "NN" & input_form != "0000"){
      stop("Only 0000 and NN formats accepted.")
    }
    if(is.table(pop)){
      cat("Pops:")
      for (i in 1:length(pop)){
        cat("\n\t", names(pop)[i], ", size: ", as.numeric(pop[i]))
      }
      cat("\n")
    }
  }
  else if(output == 6){
    cat("Converting to NN format.\n")
    if(input_form != "0000" & input_form != "snp_tab"){
      stop("Only 0000 and snp_tab formats accepted.")
    }
  }
  else if(output == 7){
    cat("Converting to allele presence/absense format.\n")
  }
  else if(output == 8){
    if(input_form != "0000" & input_form != "snp_tab" & input_form != "NN"){
      stop("Only SNP data accepted for RAFM format.")
    }
    cat("Converting to RAFM format.\n")
  }
  else if(output == 10){
    if(input_form != "0000" & input_form != "snp_tab" & input_form != "NN"){
      stop("Only SNP data accepted for dadi format.")
    }
    if(!("ref" %in% colnames(data)) | !("anc" %in% colnames(data))){
      stop("Reference and ancestor/outgroup flanking bases required in columns named 'ref' and 'anc', respecitvely")
    }
    cat("Converting to dadi format...\n")
  }
  else if(output == 11){
    if(!all(c("position", "group", "snp") %in% colnames(data))){
      stop("Columns named position, group, and snp containing position in bp, chr/linkage group/scaffold, and snp ID must be present in x!")
    }
    if(!is.null(ped)){
      if(!is.data.frame(ped)){
        stop("ped must be a six column data frame containg Family ID, Individual ID, Paternal ID, Maternal ID, Sex, and Phenotype and one row per sample. See plink documentation.\n")
      }
      if(ncol(ped) != 6 | nrow(ped) != (ncol(data)-ecs)){
        stop("ped must be a six column data frame containg Family ID, Individual ID, Paternal ID, Maternal ID, Sex, and Phenotype and one row per sample. See plink documentation.\n")
      }
    }

    cat("Converting to PLINK! binary format.")
  }
  else if(output == 12){
    cat("Converting to single character numeric format.\n")
  }

  else{
    stop("Please specify output format.")
  }

  if(input_form == "NN"){
    cat("Input format: NN\n")
    if(nchar(miss) != 1){
      stop("Missing data format must be one character.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "0000"){
    cat("Input format: 0000\n")
    if(nchar(miss) != 2){
      stop("Missing data format must be two characters.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "msat_2"){
    cat("Input format: msat, 2 character (0000)\n")
    if(nchar(miss) != 2){
      stop("Missing data format must be two characters.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "msat_3"){
    cat("Input format: msat, 3 character (000000)\n")
    if(nchar(miss) != 3){
      stop("Missing data format must be three characters.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "snp_tab"){
    cat("Input format: snp_tab.\n")
    if(nchar(miss) != 1){
      stop("Missing data format must be one character.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "0_geno"){
    cat("Imput format: 0_geno\n")
    if(miss %in% c(0:2)){
      stop("Missing data format must be other than 0, 1, or 2.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else{
    stop("Unsupported input format.")
  }


  if(outfile != FALSE){
    if(is.character(outfile) & length(outfile) == 1){
      cat("Printing results to:", outfile, "\n")
      if(file.exists(outfile)){
        #ask for confirmation before proceeding, since shit will be overwritten.
        cat("Outfile already exits. ")
        resp <- "empty"
        while(resp != "y" & resp != "n"){
          cat("Overwrite? (y or n)\n")
          resp <- readLines(n = 1)
        }
        if(resp == "n"){
          stop("Please provide acceptable path to file for output.\n")
        }
        else{
          cat("\tProceeding with conversion...\n")
        }
      }
    }
    else{
      stop("Outfile must be a character vector of length 1.\n")
    }
  }

  #if an in tab is provided or an out.tab is requested but the output format isn't 1, 5, or 10, stop.
  if((is.list(in.tab) | out.tab == T) & !(output %in% c(1,5,10))){
    stop("Input and output allele/genotype tables are not supported in formats other than ac, hapmap, or dadi.")
  }


  #####################################################################
  #conversions

  #convert 0000 to NN unless output is 7 or 2. Return after if output is 6. (v)
  if(input_form == "0000" & output != 7 & output != 2){

    #Do conversion
    cat("Converting genotypes to NN form intermediary...")

    #vectorize and replace
    xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
    xv1 <- substr(xv, 1, 2)
    xv2 <- substr(xv, 3, 4)
    xv1[xv1 == "01"] <- "A"
    xv1[xv1 == "02"] <- "C"
    xv1[xv1 == "03"] <- "G"
    xv1[xv1 == "04"] <- "T"
    xv1[xv1 == miss] <- "N"
    xv2[xv2 == "01"] <- "A"
    xv2[xv2 == "02"] <- "C"
    xv2[xv2 == "03"] <- "G"
    xv2[xv2 == "04"] <- "T"
    xv2[xv2 == miss] <- "N"
    xv <- paste0(xv1, xv2)

    #rebind to matrix and remake data.
    xv <- matrix(xv, nrow(data), (ncol(data) - ecs), T)
    cnames <- colnames(data)
    data <- cbind(data[,1:ecs], as.data.frame(xv, stringsAsFactors = F))
    colnames(data) <- cnames

    cat("\nMoving on to conversion...", "\n")
    miss <- "N" #reset miss to the correct entry
    if(output == 6){
      rdata <- data #all done if just converting to NN
    }
  }

  else if (input_form == "NN"){
    cat("Ensure that all data columns are character vectors!\n")
  }

  #convert snp_tab to NN (v)
  if(input_form == "snp_tab"){
    header <- data[,1:ecs]
    xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
    nsamp <- ncol(data) - ecs
    snames <- colnames(data)[(ecs+1):ncol(data)]
    rm(data)
    ptf <- nchar(xv) == 1
    xv[ptf] <- paste0(xv[ptf], xv[ptf]) #double up homozygotes
    remove(ptf)
    xv <- gsub(" ", "", xv) #combine heterozygotes.
    xv[xv == paste0(miss, miss)] <- "NN" #replace with correct missing data
    xv <- as.data.frame(matrix(xv, nrow(header), nsamp, byrow = T), stringsAsFactors = F)
    colnames(xv) <- snames
    out <- cbind(header, xv)
    if(output == 6){
      rdata <- out
    }
    else{
      data <- out
    }
  }

  #convert "0_geno" format to NN
  if(input_form == "0_geno"){
    header <- data[,1:ecs]
    # save for plink if that's the destination format!
    if(output == 11){
      pl_0g_dat <- data[,(ecs + 1):ncol(data)]
    }
    xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
    xv[xv == 0] <- "AA"
    xv[xv == 1] <- "AC"
    xv[xv == 2] <- "CC"
    xv[xv == miss] <- "NN"

    xv <- matrix(xv, nrow(data), (ncol(data) - ecs), T)
    cnames <- colnames(data)
    data <- cbind(data[,1:ecs], as.data.frame(xv, stringsAsFactors = F))
    colnames(data) <- cnames

    cat("\nMoving on to conversion...", "\n")
    miss <- "N" #reset miss to the correct entry
    if(output == 6){
      rdata <- data #all done if just converting to NN
    }
  }


  ##convert to allele count, migrate-n, or dadi format. Migrate-n should ALWAYS have multiple pops (why else would you
  ##use it?) (v)
  if(output == 1 | output == 5 | output == 10){
    if(output == 5){cat("WARNING: Data does not have header or pop spacer rows.\n")}

    #if requested, prepare output allele/genotype matrices.
    if(out.tab == TRUE){
      if(!is.table(pop)){
        save.tab <- list(emats = vector(mode = "list", length = 3))
        names(save.tab$emats) <- c("amat", "gmat", "tmat")
        if("pop.emats" %in% names(in.tab)){
          save.tab$pop.emats <- in.tab$pop.emats
        }
      }
      else{
        save.tab <- list(pop.emats = vector(mode = "list", length = length(pop)))
        names(save.tab$pop.emats) <- names(pop)
        if("emats" %in% names(in.tab)){
          save.tab$emats <- in.tab$emats
        }
      }
    }

    #no pop info
    if (is.table(pop) == FALSE){
      if(output == 5){stop("Cannot convert to migrate-n format with only one pop.\n")}

      #prepare a genotype table
      ##if tables are provided, pull out and store the correct elements
      if(is.list(in.tab)){
        tabs <- list(as = in.tab$emats$amat, gs = in.tab$emats$tmat, wm = in.tab$emats$gmat)
      }
      ##otherwise make a new table.
      else{
        tabs <- tabulate_genotypes(data[,-c(1:ecs)], paste0(miss, miss), verbose = T)
      }

      #save results if requested
      if(out.tab == TRUE){
        save.tab$emats$amat <- tabs$as
        save.tab$emats$gmat <- tabs$wm
        save.tab$emats$tmat <- tabs$gs
      }

      #get the number of observed alleles per loci
      counts <- matrixStats::rowSums2(ifelse(tabs$as == 0, 0, 1))

      #check for bad loci
      if(any(counts != 2)){
        if(any(counts > 2)){
          vio <- which(counts > 2)
          warning("More than two alleles detected at some loci.\n List of violating loci...")
          return(vio)
        }
        else{
          vio <- which(counts != 2)
          warning("Some loci are not bi-allelic, see printed list.")
          if(output == 10){stop()}
          print(vio)
        }
      }

      #ac output
      if(output == 1 | output == 5){
        w_df <- data[,1:ecs] #intiallize w_df if doing option 1
        #fill in the table
        w_df$n_total <- rowSums(tabs$as)
        w_df$n_alleles <- rowSums(tabs$as != 0)
        #fill in ni1 and ni2.
        w_df$ni1 <- tabs$as[,"A"]
        w_df$ni1[w_df$ni1 == 0] <- tabs$as[which(w_df$ni1== 0), "G"]
        w_df$ni1[w_df$ni1 == 0] <- tabs$as[which(w_df$ni1== 0), "C"]
        w_df$ni2 <- tabs$as[,"T"]
        w_df$ni2[w_df$ni2 == 0] <- tabs$as[which(w_df$ni2== 0), "C"]
        w_df$ni2[w_df$ni2 == 0] <- tabs$as[which(w_df$ni2== 0), "G"]

        #fix anything that's not bi-allelic
        w_df$ni1[w_df$n_alleles == 1 & w_df$ni1 == 0] <- w_df$ni2[w_df$n_alleles == 1 & w_df$ni1 == 0]
        w_df$ni2[w_df$n_alleles == 1] <- 0
      }


      #dadi output
      else{
        tabs$as <- tabs$as[,order(colnames(tabs$as))] #order the colums
        tabs$as <- cbind(tabs$as, isort = 1:nrow(tabs$as))
        #any entries with a maf of 0.5? If so, need to fix those rows later!
        maf.f <- tabs$as
        maf.f[maf.f == 0] <- NA
        vars <- matrixStats::rowVars(maf.f[,1:4], na.rm = TRUE) #if variance is zero, maf is 0.5
        fix.i <- which(vars == 0)
        if(length(fix.i) > 0){
          bmafs <- tabs$as[fix.i,]
          tabs$as <- tabs$as[-fix.i,] #remove them for now!

          #figure out min/maj allele for these...
          if(is.matrix(bmafs)){
            mas <- which(bmafs[,1:4] != 0)
          }
          else{
            mas <- which(bmafs[1:4] != 0)
          }
          mas <- mas %% 4
          mas <- gsub("1", "A", mas) #replace with correct alleles
          mas <- gsub("2", "C", mas)
          mas <- gsub("3", "G", mas)
          mas <- gsub("0", "T", mas)
        }

        #figure out the major alleles
        maxes.r <- matrixStats::rowMaxs(tabs$as[,1:4]) #get the max in each row
        maxes <- rep(maxes.r, each = 4) #make a vector of these, with each element repeated four times
        max.i <- which(maxes == t(tabs$as[,1:4])) #figure out which indices match the as table.
        max.i <- max.i %% 4 #figure out which column this was in (0 means col 4)
        max.i <- gsub("1", "A", max.i) #replace with correct alleles
        max.i <- gsub("2", "C", max.i)
        max.i <- gsub("3", "G", max.i)
        max.i <- gsub("0", "T", max.i)

        #figure out the minor alleles
        min.i <- which(t(tabs$as[,1:4]) != maxes & t(tabs$as[,1:4]) != 0) #if they aren't the major, but aren't zero...
        min.i <- min.i %% 4 #figure out which column this was in (0 means col 4)
        min.i <- gsub("1", "A", min.i) #replace with correct alleles
        min.i <- gsub("2", "C", min.i)
        min.i <- gsub("3", "G", min.i)
        min.i <- gsub("0", "T", min.i)


        #add back in the maf calls
        if(length(fix.i) > 0){
          #fix major and minor allele calls
          ##make a data frame of the min.i, max.i and isort for everything
          repl_dat <- cbind(max.i, min.i, tabs$as[,5])
          repl_dat <- rbind(repl_dat, cbind(matrix(mas, ncol = 2, byrow = T), fix.i))
          repl_dat <- repl_dat[order(as.numeric(repl_dat[,3])),]
          max.i <- repl_dat[,1]
          min.i <- repl_dat[,2]

          #add them back in
          tabs$as <- rbind(tabs$as, bmafs)
          tabs$as <- tabs$as[order(tabs$as[,5]),]
          tabs$as <- tabs$as[,-5]
        }

        meta <- data[,1:ecs]

        #save stuff in the correct order
        maxes <- matrixStats::rowMaxs(tabs$as)
        w_df <- cbind(ref = data$ref,
                      anc = data$anc,
                      Allele1 = max.i,
                      pop = maxes,
                      Allele2 = min.i,
                      pop = matrixStats::rowSums2(tabs$as) - maxes,
                      meta[,!(colnames(meta) %in% c("ref", "anc"))]
        )
      }
    }
    #if populations are specified...
    else{
      pop_count <- length(pop)
      pop_sizes <- pop
      if(sum(pop) != (ncol(data) - ecs)){#checks to make sure that the pop sizes add up to agree
        #with data
        cat("Supplied population numbers do not equal the number of supplied loci columns. Exiting.")
        stop()
      }

      #initialize w_df for multiple pops ifac or migrate-n
      if(output == 1 | output == 5){
        j <- 2
        w_df <- data[,1:ecs]
        w_df$pop <- names(pop)[1] #create pop column and set first pop name
        wa_df <- w_df #set first temp w_df
        for (j in j:pop_count){ #loop through each pop name
          wb_df <- wa_df #create temp df copy of wa_df
          wb_df$pop <- names(pop)[j]#change pop to current name
          w_df <- rbind(w_df, wb_df) #bind this copy to w_df
        }
        remove(wa_df, wb_df) #remove temp df

        #initialize columns
        w_df$n_total <- NA
        w_df$n_alleles <- NA
        w_df$ni1 <- NA
        w_df$ni2 <- NA
      }



      #create allele count table for each pop, fill their section of data

      #build allele tables for each locus
      pop_as <- list()
      pals <- matrix(FALSE, nrow(data), 4)
      colnames(pals) <- c("A", "C", "G", "T")
      current <- ecs + 1
      cat("Generating or grabbing population allele count matrices, current pop:\n")
      for(i in 1:pop_count){
        cat("\t", names(pop)[i], "\n")

        #get the data for this population

        ne <- current + pop_sizes[i] - 1
        x <- data[,current:ne]
        ##if input tables for populations were provided,pull the info out
        if("pop.emats" %in% names(in.tab)){
          which_pop <- which(names(in.tab$pop.emats) == names(pop)[i])
          temp_tab <- in.tab$pop.emats[[which_pop]]
          pop_as[[i]] <- temp_tab$amat
        }
        ##otherwise get the info
        else{
          temp_tab <- tabulate_genotypes(x, paste0(miss, miss))
          pop_as[[i]] <- temp_tab$as
        }
        ##re-order
        pop_as[[i]] <- pop_as[[i]][, order(colnames(pop_as[[i]]))]


        #save table if requested
        if(out.tab == TRUE){
          which_pop <- which(names(save.tab$pop.emats) == names(pop)[i])
          save.tab$pop.emats[[which_pop]]$amat <- pop_as[[i]]
          save.tab$pop.emats[[which_pop]]$tmat <- temp_tab$gs
          save.tab$pop.emats[[which_pop]]$gmat <- temp_tab$wm
          names(save.tab$pop.emats)[i] <- names(pop)[i]
        }


        #if somehow not all four possible alleles found in ALL DATA for this pop (small?), fill in table
        if(!"A" %in% colnames(pop_as[[i]])){pop_as[[i]] <- cbind(A = numeric(nrow(pop_as)), pop_as[[i]])}
        if(!"C" %in% colnames(pop_as[[i]])){pop_as[[i]] <- cbind(pop_as[[i]][,1], C = numeric(nrow(pop_as[[i]])), pop_as[[i]][,-1])}
        if(!"G" %in% colnames(pop_as[[i]])){pop_as[[i]] <- cbind(pop_as[[i]][,1:2], G = numeric(nrow(pop_as[[i]])), pop_as[[i]][,-c(1,2)])}
        if(!"T" %in% colnames(pop_as[[i]])){pop_as[[i]] <- cbind(pop_as[[i]][,1:3], T = numeric(nrow(pop_as[[i]])))}

        #figure out which alleles are present at each locus to write these.
        pals <- pals + (pop_as[[i]] != 0)
        current <- current + pop_sizes[i]

      }

      #now write the correct counts
      cat("Tabulating and writing results...\n")
      pals <- pals != 0
      if(any(rowSums(pals) != 2)){
        if(any(rowSums(pals) > 2)){
          stop("Some loci have more than three alleles. Check/Filter data!")
        }
        else{
          warning("Some loci have less than two alleles.\n")
        }
      }

      #vector which says which are fixed
      fixed <- rowSums(pals) == 1

      #if any SNPs are fixed in all pops and doing dadi, kill the job
      if(any(fixed) & output == 10){
        stop(cat("Fixed snps detected:", which(fixed), "stopping."))
      }

      #Convert into vector which says which elements to keep.
      if(any(fixed)){
        palsv <- as.vector(t(pals[fixed == FALSE,]))
      }
      else{
        palsv <- as.vector(t(pals))
      }


      #if doing dadi, figure out major and minor before proceeding as before:
      if(output == 10){

        #in this case, we need to sum all of the pop matrices...
        tabs <- pop_as[[1]]
        for(i in 2:length(pop_as)){
          tabs <- tabs + pop_as[[i]]
        }

        #then proceed as before.
        tabs <- tabs[,order(colnames(tabs))] #order the colums
        tabs <- cbind(tabs, isort = 1:nrow(tabs))
        #any entries with a maf of 0.5? If so, need to fix those rows later!
        maf.f <- tabs
        maf.f[maf.f == 0] <- NA
        vars <- matrixStats::rowVars(maf.f[,1:4], na.rm = TRUE) #if variance is zero, maf is 0.5
        fix.i <- which(vars == 0)
        if(length(fix.i) > 0){
          bmafs <- tabs[fix.i,]
          tabs <- tabs[-fix.i,] #remove them for now!

          #figure out min/maj allele for these...
          if(is.matrix(bmafs)){
            mas <- which(bmafs[,1:4] != 0)
          }
          else{
            mas <- which(bmafs[1:4] != 0)
          }
          mas.index <- mas #save for later
          mas <- mas %% 4
          mas <- gsub("1", "A", mas) #replace with correct alleles
          mas <- gsub("2", "C", mas)
          mas <- gsub("3", "G", mas)
          mas <- gsub("0", "T", mas)
        }

        #figure out the major alleles
        maxes.r <- matrixStats::rowMaxs(tabs[,1:4]) #get the max in each row
        maxes <- rep(maxes.r, each = 4) #make a vector of these, with each element repeated four times
        max.i <- which(maxes == t(tabs[,1:4])) #figure out which indices match the as table.
        max.index <- max.i #save for later
        max.i <- max.i %% 4 #figure out which column this was in (0 means col 4)
        max.i <- gsub("1", "A", max.i) #replace with correct alleles
        max.i <- gsub("2", "C", max.i)
        max.i <- gsub("3", "G", max.i)
        max.i <- gsub("0", "T", max.i)

        #figure out the minor alleles
        min.i <- which(t(tabs[,1:4]) != maxes & t(tabs[,1:4]) != 0) #if they aren't the major, but aren't zero...
        min.index <- min.i #save for later
        min.i <- min.i %% 4 #figure out which column this was in (0 means col 4)
        min.i <- gsub("1", "A", min.i) #replace with correct alleles
        min.i <- gsub("2", "C", min.i)
        min.i <- gsub("3", "G", min.i)
        min.i <- gsub("0", "T", min.i)

        #add back in the maf calls
        if(length(fix.i) > 0){
          #fix major and minor allele calls
          ##make a data frame of the min.i, max.i and isort for everything
          repl_dat <- cbind(max.i, min.i, tabs[,5])
          repl_dat <- rbind(repl_dat, cbind(matrix(mas, ncol = 2, byrow = T), fix.i))
          repl_dat <- repl_dat[order(as.numeric(repl_dat[,3])),]
          max.i <- repl_dat[,1]
          min.i <- repl_dat[,2]

          #save indices for max and min counts for later!
          mas.index.mat <- matrix(mas.index + length(tabs[,-5]), ncol = 2, byrow = T) #the indices of the "major" and "minor" alleles, if they are tacked on to the end of the TRANSPOSED tabs matrix.
          mas.maj <- mas.index.mat[,1]
          mas.min <- mas.index.mat[,2]
          max.index.out <- c(max.index, mas.maj)
          min.index.out <- c(min.index, mas.min)


          #add the data back in and resort.
          tabs <- rbind(tabs, bmafs)
          tabs <- tabs[order(tabs[,5]),]
          tabs <- tabs[,-5]
        }
        else{
          max.index.out <- max.index.out
          min.index.out <- min.index
        }

        #also initialize output
        w_df <- matrix(NA, ncol = (ecs + 2 + length(pop)*2), nrow = nrow(data))
        w_df[,1] <- as.character(data$ref)
        w_df[,2] <- as.character(data$anc)
        w_df[,3] <- max.i
        w_df[,(4+length(pop))] <- min.i
        meta <- data[,1:ecs]
        meta <- meta[,!(colnames(meta) %in% c("ref", "anc"))]
        w_df <- as.data.frame(w_df, stringsAsFactors = F)
        w_df[,(ncol(w_df) - ecs + 3):ncol(w_df)] <- meta
        colnames(w_df) <- c("ref", "anc", "Allele1", names(pop), "Allele2", names(pop), colnames(meta))
      }



      #loop through and print data
      for(i in 1:pop_count){
        #have to do this slightly differently if there are any non-polymorphic bases
        if(any(fixed) > 0){
          #do the fixed ones first
          w_df[w_df$pop == names(pop)[i] & fixed,]$n_total <- rowSums(pop_as[[i]][fixed,])
          w_df[w_df$pop == names(pop)[i] & fixed,]$n_alleles <- 1 #in this case it is fixed in all pops
          w_df[w_df$pop == names(pop)[i] & fixed,]$ni1 <- rowSums(pop_as[[i]][fixed,])
          w_df[w_df$pop == names(pop)[i] & fixed,]$ni2 <- 0

          #do the others
          if(any(!fixed)){
            av <- as.vector(t(pop_as[[i]][!fixed,]))
            av <- av[palsv]
            av <- matrix(av, nrow(data[!fixed,]), 2, byrow = T)
            w_df[w_df$pop == names(pop)[i] & !fixed,]$n_total <- rowSums(av)
            w_df[w_df$pop == names(pop)[i] & !fixed,]$n_alleles <- 2 #not fixed in at least one pop!
            w_df[w_df$pop == names(pop)[i] & !fixed,]$ni1 <- av[,1]
            w_df[w_df$pop == names(pop)[i] & !fixed,]$ni2 <- av[,2]
          }
        }

        else{
          if(output == 1 | output == 5){
            #compare palsv to allele tables to keep correct alleles.
            av <- as.vector(t(pop_as[[i]]))
            av <- av[palsv]
            av <- matrix(av, nrow(data), 2, byrow = T)
            w_df[w_df$pop == names(pop)[i],]$n_total <- rowSums(av)
            w_df[w_df$pop == names(pop)[i],]$n_alleles <- 2 #not fixed in at least one pop!
            w_df[w_df$pop == names(pop)[i],]$ni1 <- av[,1]
            w_df[w_df$pop == names(pop)[i],]$ni2 <- av[,2]
          }

          #dadi
          else{
            #if we have values to fix
            if(length(fix.i) > 0){
              #figure out which values to take for min and maj allele calls
              t.tab <- pop_as[[i]] #get this table
              t.tab <- cbind(t.tab, index = 1:nrow(t.tab)) #add an index variable
              t.maf <- t.tab[fix.i,] #subset the ones to fix
              t.tab <- t.tab[-fix.i,] #remove them
              t.tab <- rbind(t.tab, t.maf) #bind the ones to fix to the bottom
              t.counts <- cbind(t(t.tab[,-5])[max.index.out], t(t.tab[,-5])[min.index.out]) #get out major and minor counts
              t.counts <- t.counts[order(t.tab[,5]),] #re-order them
            }
            #otherwise easier
            else{
              t.counts <- cbind(t(pop_as[[i]])[max.index.out], t(pop_as[[1]])[min.index.out]) #get out major and minor counts
            }

            #add data.
            w_df[,3 + i] <- t.counts[,1]
            w_df[,4 + length(pop) + i] <- t.counts[,2]
          }
        }
      }

      row.names(w_df) <- 1:nrow(w_df) #fix row names
      if(output == 5){ #if doing output style 5...
        w_df <- w_df[,c(1, (ecs + 1):ncol(w_df))] #remove extra columns except the first and pop columns
      }
    }
    rdata <- w_df
  }



  ##convert to genepop or numeric format (v)
  if (output == 2 | output == 4){
    cat("WARNING: For output option 3, data.frame output does not have population seperators or header information. Outfile, if requested and pop list provided, will.", "\n", "Converting genotypes...")

    if(input_form == "0000" & output == 2){
      w_df <- as.data.frame(t(data[,(ecs + 1):ncol(data)]), stringsAsFactors = F) #remove extra columns and transpose data
      row.names(w_df) <- paste0(row.names(w_df), " ,") #adding space and comma to row names, as required.
      rdata <- w_df
    }
    else if (input_form == "0000" & output == 4){
      stop("Same input form given as output selected.")
    }

    else{
      #vectorize and replace
      xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
      xv1 <- substr(xv, 1, 1)
      xv2 <- substr(xv, 2, 2)
      xv1[xv1 == "A"] <- "01"
      xv1[xv1 == "C"] <- "02"
      xv1[xv1 == "G"] <- "03"
      xv1[xv1 == "T"] <- "04"
      xv1[xv1 == miss] <- "00"
      xv2[xv2 == "A"] <- "01"
      xv2[xv2 == "C"] <- "02"
      xv2[xv2 == "G"] <- "03"
      xv2[xv2 == "T"] <- "04"
      xv2[xv2 == miss] <- "00"
      xv <- paste0(xv1, xv2)

      xv <- matrix(xv, nrow(data), (ncol(data) - ecs), T)
      inds <- colnames(data)
      data <- cbind(data[,1:ecs], as.data.frame(xv, stringsAsFactors = F))
      colnames(data) <- inds

      cat("\n", "Cleaning up...", "\n")
      if(output == 2){ #convert to genepop
        w_df <- as.data.frame(t(data[,(ecs + 1):ncol(data)]), stringsAsFactors = F) #remove extra columns and transpose data
        row.names(w_df) <- paste0(row.names(w_df), " ,") #adding space and comma to row names, as required.
        rdata <- w_df
      }
      else {#prepare numeric output, otherwise same format
        rdata <- data
      }
    }
  }


  ##convert to structure, fastStructure or RAFM format (v)
  if (output == 3 | output == 8 | output == 9){

    #subset if requested
    if(all(!is.na(n_samp))){
      cat("Subsetting ")
      if(length(n_samp) > 1){
        cat("designated SNPs.\n")
        data <- data[n_samp,]
      }
      else{
        cat(n_samp, " random SNPs.\n")
        data <- data[sample(nrow(data), n_samp, T),]
      }
    }

    #plop actual data into two vectors, one for each allele.
    xv <- as.vector(t(t(data[,(ecs + 1):ncol(data)])))
    xv1 <- substr(xv, 1, 1)
    xv2 <- substr(xv, 2, 2)
    remove(xv)

    #convert to numeric
    xv1[xv1 == "A"] <- 1
    xv1[xv1 == "C"] <- 2
    xv1[xv1 == "G"] <- 3
    xv1[xv1 == "T"] <- 4
    xv1[xv1 == miss] <- 0
    xv2[xv2 == "A"] <- 1
    xv2[xv2 == "C"] <- 2
    xv2[xv2 == "G"] <- 3
    xv2[xv2 == "T"] <- 4
    xv2[xv2 == miss] <- 0

    #bind back into matrices
    xv1 <- matrix(xv1, ncol(data) - ecs, nrow(data), T)
    xv2 <- matrix(xv2, ncol(data) - ecs, nrow(data), T)

    #create output matrix and bind the data to it, structure format
    if(output == 3 | output == 9){
      outm <- matrix(NA, 2*(ncol(data) - ecs), nrow(data))

      #fill
      outm[seq(1,nrow(outm),2),] <- xv1
      outm[seq(2,nrow(outm),2),] <- xv2

      #add sample names
      snames <- character(nrow(outm))
      snames[seq(1,nrow(outm),2)] <- colnames(data[,(ecs+1):ncol(data)])
      snames[seq(2,nrow(outm),2)] <- colnames(data[,(ecs+1):ncol(data)])
      if(output == 3){ #bind sample names in and return for structure
        out <- cbind(ind = snames, as.data.frame(outm, stringsAsFactors = F))
      }
      else{ #add a bunch of filler columns (4 if pop info is provided, 5 if it isn't) for faststructure and change missing data to -9
        outm[outm == 0] <- -9
        if(length(pop) > 1 & is.table(pop) == TRUE){
          out <- cbind(ind = snames,
                       matrix("filler", nrow(outm), 4),
                       as.data.frame(outm, stringsAsFactors = F))
          colnames(out)[2:5] <- paste0("filler", 1:4)
        }
        else{
          out <- cbind(ind = snames,
                       matrix("filler", nrow(outm), 5),
                       as.data.frame(outm, stringsAsFactors = F))
          colnames(out)[2:6] <- paste0("filler", 1:5)
        }
      }

      #add pop info if possible
      if(is.table(pop)){
        #add pop info if provided
        pop <- pop*2
        pl <- numeric(sum(pop))
        pl[1:pop[1]] <- 1
        tracker <- pop[1] + 1
        for(i in 2:length(pop)){
          pl[tracker:(sum(pop[1:i]))] <- i
          tracker <- sum(pop[1:i]) + 1
        }
        out <- cbind(out[,1], pop = pl, out[,2:ncol(out)])
        colnames(out)[1] <- "ind"
      }
    }

    #create output matrix and bind to it, RAFM format.
    else{
      outm <- matrix(NA, ncol(data) - ecs, nrow(data)*2)

      #fill
      outm[,seq(1,ncol(outm),2)] <- xv1
      outm[,seq(2,ncol(outm),2)] <- xv2

      #replance missings with NA
      outm[outm == 0] <- NA

      #add column names
      colnames(outm) <- paste0("locus", sort(rep(1:nrow(data),2)), rep(c(".1", ".2"), ncol(outm)/2))

      #add subpop numbers, if given
      if(is.table(pop)){
        pl <- numeric(sum(pop))
        pl[1:pop[1]] <- 1
        tracker <- pop[1] + 1
        for(i in 2:length(pop)){
          pl[tracker:(sum(pop[1:i]))] <- i
          tracker <- sum(pop[1:i]) + 1
        }
        out <- cbind(subpop = pl, as.data.frame(outm, stringsAsFactors = F))
      }
    }
    rdata <- out
  }


  #presence/absence format
  if(output == 7){
    x <- data[,(ecs+1):ncol(data)] #get just data
    if(input_form == "NN"){
      asize <- 1
    }
    else if (input_form == "0000" | input_form == "msat_2"){
      asize <- 2
    }
    else if (input_form == "msat_3"){
      asize <- 3
    }

    if(is.null(lnames)){
      lnames <- 1:nrow(x)
    }
    else{
      lnames <- gsub("_", "", lnames)
    }

    xv <- as.vector(t(x))

    #function to produce vectorized presence absence (as much as possible, not vectorized for NAs)
    pa_alleles <- function(xv, snp_form, mDat, nsamp, nloci){
      #get all possible genotypes
      gs <- unique(xv)

      #which genotype is the missing data?
      mpos <- which(gs == mDat)

      #what are the possible alleles at all loci?
      as <- unique(c(substr(gs,1,snp_form/2), substr(gs, (snp_form/2 + 1), snp_form*2)))
      as <- sort(as[as != substr(mDat, 1, snp_form/2)]) #that aren't N?

      #make the table
      cat("Creating presence/absence table...\n")

      #convert genotype vector to allele vector
      xva1 <- substr(xv, 1, snp_form/2)
      xva2 <- substr(xv, (snp_form/2+1),snp_form)

      #initialize
      amat <- matrix(0, nsamp, length(as)*nloci)
      colnames(amat) <- paste0(sort(rep(lnames,length(as))), "_", as) #initialize all of the columns, with locus number followed by allele. Will remove anything empty after assignment.

      #fill in
      for(i in 1:length(as)){
        pr1 <- grep(as[i], xva1)#unique rows which have the allele in either position one or position two.
        pr2 <- grep(as[i], xva2)
        amat[,grep(paste0("_", as[i]),colnames(amat))][pr1] <- amat[,grep(paste0("_", as[i]),colnames(amat))][pr1] + 1 #set the allele as present in the correct rows. This works because we look only at the G amat columns first, then put set only the individual IDs with a G to one. I think.
        amat[,grep(paste0("_", as[i]),colnames(amat))][pr2] <- amat[,grep(paste0("_", as[i]),colnames(amat))][pr2] + 1
      }

      #remove alleles not seen at loci.
      amat <- amat[,which(colSums(amat) != 0)]

      #get allele counts per loci:
      cat("Filling in missing data with NAs.\n")

      ###########
      #fill in missing data with NAs.
      #easy to vectorize if SNP data
      if(input_form == "NN" | input_form == "0000"){
        xmc <- which(xv == paste0(miss, miss)) #which samples had missing data?
        adj <- floor(xmc / nsamp) #how many loci over do I need to adjust xmc, since in amat each locus occupies two columns?
        adj[xmc%%nsamp == 0] <- adj[xmc%%nsamp == 0] - 1 #shift over anything that got screwed up by being in the last sample
        xmc <- xmc + (nsamp*adj) #adjust xmc for extra columns.
        if(any(amat[xmc] != 0) | any(amat[xmc + nsamp] != 0)){
          stop("Missing data values were not properly identified for replacement with NAs. This usually happens when SNP data is not completely bi-allelic. Try filtering out non-biallelic and non-polymorphic SNPs using filter_snps.\n")
        }
        amat[xmc] <- NA #make the first allele NA
        amat[xmc + nsamp] <- NA #make the second allele (another column over) NA.
      }

      #otherwise need to loop through loci and fill
      else{
        #get allele counts per locus
        labs <- gsub("_\\w+", "", colnames(amat))
        ll <- unique(labs)
        acounts <- rbind(label = ll, count = sapply(ll, function(x)sum(labs==x)))
        acounts <- as.numeric(acounts[2,])

        #intialize stuff, do first loci, then do the rest of the loci.
        spos <- 1 #starting column for loci
        epos <- acounts[1] #ending column for loci
        amat[, spos:epos][rowSums(amat[, spos:epos]) == 0,] <- NA
        for (i in 2:nloci){
          spos <- sum(acounts[1:(i - 1)]) + 1
          epos <- spos + acounts[i] - 1
          amat[, spos:epos][rowSums(amat[, spos:epos]) == 0,] <- NA
        }
      }
      return(amat)
    }

    amat <- pa_alleles(xv, asize*2, miss, ncol(x), nrow(x))

    if(interp_miss){
      #average number observed in columns
      cat("Interpolating missing data...\n")
      afs <- colMeans(amat, TRUE)
      temp <- which(is.na(amat))/nrow(amat)
      fill_element <- floor(temp) + 1 #get the column for each missing data point
      fill_element[which(temp %% 1 == 0)] <- fill_element[which(temp %% 1 == 0)] - 1 #correct for anything in the last row
      amat[which(is.na(amat))] <- afs[fill_element] #fill with the appropriate allele frequency.
    }
    else{cat("Finished. Warning: Missing data counts are also stored!\n")}
    amat <- cbind(samp = as.character(colnames(data)[(ecs+1):ncol(data)]), as.data.frame(amat, stringsAsFactors = F))
    rdata <- amat
  }


  #PLINK formatt
  if(output == 11){
    #============convert the genotypes into a binary string==========
    # pull out genotypes, save headers
    x <- data[,(ecs+1):ncol(data)]
    header <- data[,1:ecs]
    rm(data)

    #======================make a vector containing the bits for each genotype============
    # much easier with 0_geno!
    if(input_form == "0_geno"){

      # convert to vector
      xv <- as.vector(t(pl_0g_dat))

      # add a few filler genotype calls to fill in the empty bits later on.
      x.dup <- ncol(pl_0g_dat) %% 4
      if(x.dup  != 0){
        x.dup <- 4 - x.dup
      }
      x.dup <- matrix("FF", nrow(pl_0g_dat), x.dup)
      x.dup <- cbind(pl_0g_dat, x.dup)
      xv.dup <- as.vector(t(x.dup))

      # define the vector that we're going to fill with bits for each genotype
      xvt <- character(length(xv.dup))

      xvt[xv.dup == "0"] <- "00"
      xvt[xv.dup == "1"] <- "01"
      xvt[xv.dup == "2"] <- "11"
      xvt[xv.dup == miss] <- "10"
      xvt[xv.dup == "FF"] <- "00"

      # get allele names for map file down the line
      a.names <- matrix(c(0,1), nrow(header), 2, T)
    }

    # harder for the NN format
    else{
      # convert to vector
      xv <- as.vector(t(x))

      # add a few filler genotype calls to fill in the empty bits later on.
      x.dup <- ncol(x) %% 4
      if(x.dup  != 0){
        x.dup <- 4 - x.dup
      }
      x.dup <- matrix("FF", nrow(x), x.dup)
      x.dup <- cbind(x, x.dup)
      xv.dup <- as.vector(t(x.dup))

      # define the vector that we're going to fill with bits for each genotype
      xvt <- character(length(xv.dup))

      #pull out tab if provided
      if(is.list(in.tab)){
        tabs <- list(as = in.tab$emats$amat, gs = in.tab$emats$tmat, wm = in.tab$emats$gmat)
        as <- tabs$as
      }
      #otherwise make a new table.
      else{
        tabs <- tabulate_genotypes(x, paste0(miss, miss), verbose = T)
        as <- tabs$as
      }

      #save results if requested
      if(out.tab == TRUE){
        save.tab$emats$amat <- tabs$as
        save.tab$emats$gmat <- tabs$wm
        save.tab$emats$tmat <- tabs$gs
      }

      # if there are any alleles without any data, warn and remove
      if(any(rowSums(as) == 0)){
        missing.dat <- which(rowSums(as) == 0)
        warning(paste0("Removed ", length(missing.dat), "SNPs with no called genotypes.!\n"),
                "Removed SNPs:\n", missing.dat)
        x <- x[-missing.dat,]
        as <- as[-missing.dat,]
        xv <- as.vector(t(x))
      }

      # get the unique alleles for each loci, fixing for loci with one allele.
      unique.as <- cbind(as, N1 = rep(0, nrow(as))) # put in dummies for rows with missing haplotypes
      unique.as <- ifelse(unique.as != 0, TRUE, FALSE) # convert to logical
      if(any(rowSums(unique.as) != 2)){
        unique.as[rowSums(unique.as) != 2, 5] <- TRUE # add one missing allele note
      }
      unique.as <- which(t(unique.as)) %% ncol(unique.as) # get the column for each TRUE
      unique.as[unique.as == 0] <- 5 # fix for sixth column
      unique.as <- c(colnames(as), "M")[unique.as] # add column names.

      # match to observed genotypes
      ## these lines generate a vector equal in length to xv that has the correct genotype options per loci for all three genotypes.
      v00 <- paste0(unique.as[seq(1, length(unique.as), 2)],
                    unique.as[seq(1, length(unique.as), 2)])
      v00 <- rep(v00, each = ncol(x.dup))
      v11 <- paste0(unique.as[seq(2, length(unique.as), 2)],
                    unique.as[seq(2, length(unique.as), 2)])
      v11 <- rep(v11, each = ncol(x.dup))
      v01a <- paste0(unique.as[seq(1, length(unique.as), 2)],
                     unique.as[seq(2, length(unique.as), 2)])
      v01a <- rep(v01a, each = ncol(x.dup))
      v01b <- paste0(unique.as[seq(2, length(unique.as), 2)],
                     unique.as[seq(1, length(unique.as), 2)])
      v01b <- rep(v01b, each = ncol(x.dup))

      ## these lines create xvt, which has the correct binary doublets.
      xvt[xv.dup == v00] <- "00"
      xvt[xv.dup == v11] <- "11"
      xvt[xv.dup == v01a] <- "01"
      xvt[xv.dup == v01b] <- "01"
      xvt[xv.dup == "NN"] <- "10"
      xvt[xv.dup == "FF"] <- "00"

      # get allele names for map file down the line
      a.names <- matrix(unique.as, nrow(header), 2, T)
    }

    #=====================add magic numbers and reorder vector==============

    # order correctly for plink (weird reversed bytes)
    xvt <- paste0(xvt[seq(4, length(xvt), 4)],
                  xvt[seq(3, length(xvt), 4)],
                  xvt[seq(2, length(xvt), 4)],
                  xvt[seq(1, length(xvt), 4)])

    # add magic number and SNP-major identifier
    bed <- c("01101100", "00011011", "00000001", xvt)


    #===============make a ped file=================
    # make an empty set of ped header columns if not provided
    if(is.null(ped)){
      ped <- data.frame(fam = rep("NA", ncol(x)),
                             ind = colnames(x),
                             PatID = rep("NA", ncol(x)),
                             MatID = rep("NA", ncol(x)),
                             Sex = rep("NA", ncol(x)),
                             Phenotype = rep("NA", ncol(x)))
    }

    # save .fam
    fam <- ped

    # change missing data value and add a space between alleles.
    x.np <- as.vector(t(x))
    x.np[x.np == paste0(miss, miss)] <- "00"
    x.np <- gsub("(.)(.)", "\\1 \\2", x.np)

    # rebind
    ped <- cbind(ped, matrix(x.np, nrow(ped), nrow(x)), stringsAsFactors = F)

    #===============make an extended map file=================
    # with morgans
    if(any(colnames(header) %in% c("cM", "cm", "morgans"))){
      bim <- data.frame(chr = header$group,
                            rs = header$snp,
                            cM = header[,which(colnames(header) %in% c("cM", "cm", "morgans"))],
                            bp = header$position,
                            a1 = a.names[,1],
                            a2 = a.names[,2])
    }
    # without morgans
    else{
      bim <- data.frame(chr = header$group,
                            rs = header$snp,
                            bp = header$position,
                            a1 = a.names[,1],
                            a2 = a.names[,2])
    }

    # recode chr
    bim$chr <- as.numeric(as.factor(bim$chr))

    # grab normal map file
    map <- bim[,-c(ncol(bim)-1, ncol(bim))]

    # name output
    rdata <- list(ped = ped, bed = bed, map = map, bim = bim)
  }


  # single-character numeric format
  if(output == 12){
    x <- data[,-c(1:ecs)]

    # generate allele info to get maj_min info. The strategy here is to figure out the row maxes, check which indices
    # contain those, get their column position, and then get the name of that column. For the minor alleles, same, but
    # check that the index contains a number but isn't the maximum.
    gs <- tabulate_genotypes(x, "NN")
    maj <- rep(matrixStats::rowMaxs(gs$as), each = ncol(gs$as))
    maj <- gs$as == matrix(maj, ncol = ncol(gs$as), byrow = T)
    # deal with special cases where two alleles are of equal frequency:
    if(any(rowSums(maj) != 1)){
      eqf <- which(rowSums(maj) != 1)

      eq.min <- which(t(maj[eqf,])) %% ncol(maj)
      eq.maj <- eq.min[seq(1, length(eq.min), by = 2)]
      eq.min <- eq.min[seq(2, length(eq.min), by = 2)]

      eq.min[eq.min == 0] <- 4
      eq.maj[eq.maj == 0] <- 4
      eq.maj <- colnames(gs$as)[eq.maj] # these are the major alleles
      eq.min <- colnames(gs$as)[eq.min] # these are the minor alleles

      maj <- maj[-eqf,]

      # create the normal min
      min <- which(t(as.logical(gs$as[-eqf,]) + maj) == 1)
    }
    else{
      min <- which(t(as.logical(gs$as) + maj) == 1)
    }
    min <- min %% ncol(maj)
    min[min == 0] <- 4
    maj <- which(t(maj)) %% ncol(maj)
    maj[maj == 0] <- 4
    maj <- colnames(gs$as)[maj] # these are the major alleles
    min <- colnames(gs$as)[min] # these are the minor alleles
    # re-insert the exceptions
    if(exists("eqf")){
      maj.fix <- sort(c(1:length(maj), eqf))
      min.fix <- maj.fix
      dup <- duplicated(maj.fix,fromLast = T)
      maj.fix[dup] <- eq.maj
      maj.fix[!dup] <- maj
      maj <- maj.fix

      min.fix[dup] <- eq.min
      min.fix[!dup] <- min
      min <- min.fix
    }

    # check to see if each allele is the minor, assign a one if so
    a1 <- substr(unlist(t(x)), 1, 1)
    a2 <- substr(unlist(t(x)), 2, 2)

    a1 <- a1 == rep(min, each = ncol(x))
    a2 <- a2 == rep(min, each = ncol(x))

    # collapse to output
    rdata <- t(a1 + a2)
    rdata[as.matrix(x) == mDat] <- -1
    rdata <- cbind(data[,1:ecs], rdata, stringsAsFactors = F)
  }

  ############################################################################
  #return the final product, printing an outfile if requested.
  if(outfile != FALSE){
    cat("Writing output file...\n")
    #for genepop
    if(output == 2){
      cat("\tPreparing genepop file...\n")
      #get list of snps
      llist <- paste0("SNP", "_", 1:ncol(rdata), ",")
      llist[length(llist)] <- paste0("SNP_", ncol(rdata))

      #write output
      cat(paste0(unlist(strsplit(outfile, split =  "\\."))[1], "_genepop\n"), file = outfile)
      cat(llist, "\nPOP\n", file = outfile, append = T)

      #write the tables, splitting by pop if requested:
      if(is.table(pop)){
        cat("\tWriting genepop file seperated by populations...\t")
        #need to sort by pops. Fist loop does this:
        j <- 1
        rdata$pop <- NA
        for (i in 1:(length(pop) - 1)){
          rdata[j:(j+pop[i] - 1),]$pop <- names(pop)[i]
          j <- j + pop[i]
        }
        rdata[j:nrow(rdata),]$pop <- names(pop)[length(pop)]

        #sort and remove pop column
        rdata$rnames <- rownames(rdata)
        rdata <- dplyr::arrange(rdata, pop)
        rownames(rdata) <- rdata$rnames
        rdata <- rdata[,-which(colnames(rdata) %in% c("pop", "rnames"))]

        #second loop prints results.
        j <- 1
        pop <- pop[order(names(pop))]
        for (i in 1:(length(pop) - 1)){
          cat(names(pop)[i], "\t")
          data.table::fwrite(rdata[j:(j+pop[i] - 1),], outfile, quote = F, sep = "\t", col.names = F, row.names = T, append = T)
          cat("POP\n", file = outfile, append = T)
          j <- j + pop[i]
        }
        data.table::fwrite(rdata[j:nrow(rdata),], outfile, quote = F, sep = "\t", col.names = F, row.names = T, append = T)
        cat(names(pop)[length(pop)], "\t Done.\n")
      }
      else{
        data.table::fwrite(rdata, outfile, quote = F, sep = "\t", col.names = F, row.names = T, append = T)
      }
    }
    else if(output == 1){
      #write the raw output
      data.table::fwrite(rdata, outfile, quote = FALSE, col.names = T, sep = "\t", row.names = F)
      #write a bayescan object if pop list was provided.
      if(is.table(pop)){
        outfile <- paste0(outfile, ".bayes")
        #write the header
        cat("[loci]=", nrow(data), "\n\n[populations]=", length(pop), "\n\n", file = outfile, sep = "")

        #write the data for each population.
        for(i in 1:length(pop)){
          cat("[pop]=", i, "\n", file = outfile, append = T, sep = "") #write header
          wdat <- cbind(snp = 1:nrow(rdata[rdata$pop == names(pop)[i],]),
                        rdata[rdata$pop == names(pop)[i], which(colnames(rdata) %in% c("n_total", "n_alleles"))],
                        alleles = paste0(rdata[rdata$pop ==names(pop)[i],]$ni1, " ", rdata[rdata$pop == names(pop)[i],]$ni2, " "))
          data.table::fwrite(wdat,
                      outfile, col.names = F, row.names = F, quote = F, sep = "\t",
                      append = T) #write the data for this population.
          cat("\n", file = outfile, append = T) #add a line break
        }
      }
    }
    else if(output == 3 | output == 9){
      data.table::fwrite(rdata, outfile, quote = FALSE, col.names = F, sep = "\t", row.names = F)
    }
    else if(output == 11){
      t2 <- as.logical(as.numeric(unlist(strsplit(bed, "")))) # merge the data and convert it to a logical.
      t2 <- BMS::bin2hex(t2) # convet to hex codes
      t2 <- unlist(strsplit(t2,"")) # unlist
      t2 <- paste0("\\\\x",
                   t2[seq(1,length(t2),2)],
                   t2[seq(2,length(t2),2)],
                   collapse = "")
      writeLines(paste0("#!/bin/bash\n\n",
        "echo -n -e ", t2, " > ", outfile, ".bed"), paste0(outfile, ".sh"))
      warning(paste0("To get PLINK .bed file, run ", outfile, ".sh.\n"))
      data.table::fwrite(map, paste0(outfile, ".map"), quote = F, col.names = F, sep = "\t", row.names = F)
      data.table::fwrite(ped, paste0(outfile, ".ped"), quote = F, col.names = F, sep = "\t", row.names = F)
      data.table::fwrite(bim, paste0(outfile, ".bim"), quote = F, col.names = F, sep = "\t", row.names = F)
    }
    else{
      data.table::fwrite(rdata, outfile, quote = FALSE, col.names = T, sep = "\t", row.names = F)
    }
  }

  #return results
  if(exists("save.tab")){
    return(c(list(x = rdata), save.tab))
  }
  else{
    return(rdata)
  }
}

#function to run any command after spliting the data by group and by population. Group and population must be in columns named "group" and "pop", respectively.
#The first argument of the function to run must be the data provided to that function.
run_gp <- function(x, FUN, ...){
  w_df<- data.frame()
  x$pop <- as.character(x$pop)
  x$group <- as.character(x$group)
  pops <- unique(x[,"pop"])
  groups <- unique(x[,"group"])
  for (i in 1:length(groups)){
    w_data <- x[x[,"group"] == groups[i],]
    print(groups[i])
    for(j in 1:length(pops)){
      print(pops[j])
      x_data <- FUN(w_data[w_data[,"pop"] == pops[j],], ...)
      if(length(x_data) == 0){next}
      if(is.data.frame(x_data)){
        if(!any(colnames(x_data) == "pop"))
        x_data <- cbind(pop = pops[j], x_data, stringsAsFactors = F)
        if(!any(colnames(x_data) == "group")){
          x_data <- cbind(group = groups[i], x_data, stringsAsFactors = F)
        }
        w_df <- rbind(w_df, x_data)
      }
      else{
        if(is.data.frame(w_df)){
          w_df <- numeric(0)
        }
        w_df <- c(w_df,x_data)
      }
    }
  }
  if(!is.data.frame(w_df)){
    warning("Output is vector: if merging with meta info, check that info data frame is sorted identically:
            Group, Pop, then snp by snp!")
  }
  if(sum(colnames(w_df) == "pop") > 1){
    excols <- which(colnames(w_df) == "pop")[-1]
    w_df <- w_df[,-excols]
  }
  return(w_df)
}

#function to run any command after spliting the data by group. Group must be in a column named "group".
#The first argument of the function to run must be the data provided to that function.
run_g <- function(x, FUN, ..., y = NULL){
  w_df<- data.frame()
  x$group <- as.character(x$group)
  groups <- unique(x[,"group"])
  for (i in 1:length(groups)){
    w_data <- x[x[,"group"] == groups[i],]
    if(!is.null(y)){
      y_dat <- y[y[,"group"] == groups[i],]
      w_data <- FUN(w_data, y = y_dat, ...)
    }
    else{
      w_data <- FUN(w_data, ...)
    }
    print(groups[i])
    if(length(w_data) == 0){next}
    if(is.data.frame(w_data)){
      if(!any(colnames(w_data) == "group")){
        w_data <- cbind(group = groups[i], w_data, stringsAsFactors = F)
      }
      w_df <- rbind(w_df, w_data)
    }
    else{
      if(is.data.frame(w_df)){
        w_df <- numeric(0)
      }
      w_df <- c(w_df, w_data)
    }
  }
  if(!is.data.frame(w_df)){
    warning("Output is vector: if merging with meta info, check that info data frame is sorted identically:
            group then snp by snp!\n")
  }
  return(w_df)
}


#'Get random snps each window.
#'
#'\code{rgap_snps} grabs n random SNPs every window of a given length. Can do this across multiple linkage groups/chromosomes.
#'
#'Possible methods:
#'\itemize{
#'    \item{window: }{Takes n random SNPs every window. This can result in SNPs that are much nearer to each other than the gap length.}
#'    \item{gap: }{Takes the first possible SNP, then takes the next SNP at least gap bases away. Non-random.}
#'}
#'
#' @param x data.frame. Input data, must contain a column titled 'position' containing genomic positions of the SNPs.
#' @param gap Integer. How big are the windows from which to sample SNPs?
#' @param method Character, default "window". By which method should the SNPs be chosen?
#' @param n Integer, default 1. For "window" method. How many SNPs should be taken per window?
#' @param levels Character, default "group". Across which levels should the windows be partitioned (ex. "chromosome")?
#'
#' @return A data.frame containing subset SNPs.
#'
rgap_snps <- function(x, gap, method = "gap", n = 1, levels = "group"){

  if(method == "window"){
    gfun <- function(x){
      cs <- seq(gap/2, max(x$position) + gap/2, gap)
      starts <- cs - gap/2
      ends <- cs + gap/2
      pos <- x$position

      lmat <- outer(pos, starts, function(pos, starts) pos >= starts)
      lmat <- lmat + outer(pos, ends, function(pos, ends) pos < ends)
      colnames(lmat) <- cs
      rownames(lmat) <- pos
      lmat <- ifelse(lmat == 2, TRUE, FALSE)

      wins <- which(t(lmat) == TRUE) %% length(cs)
      wins[wins == 0] <- length(cs)

      x <- cbind(window = wins, x)

      rands <- x %>% group_by(window) %>% sample_n(size = n)

      rands <- rands[,-1]

      return(rands)
    }
  }
  else if (method == "gap"){
    gfun <- function(x){
      fsnp <- which.min(x$position)
      pos <- x$position
      cs <- pos[fsnp]
      #figure out windows
      while(cs[length(cs)] + gap <= max(pos)){
        cs <- c(cs,
                min(pos[pos >= cs[length(cs)] + gap]))
      }
      ret <- match(cs, pos)
      return(x[ret,])
    }
  }

  if(levels != FALSE){
    ret <- plyr::ddply(x, .variables = "group", .fun = gfun, .progress = "text")
  }
  else{
    ret <- gfun(x)
  }
  return(ret)
}

