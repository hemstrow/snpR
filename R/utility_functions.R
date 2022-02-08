
#'Subset snpRdata objects
#'
#'Subsets snpRdata objects by specific snps, samples, facets, subfacets, ect.
#'Statistics will need to be recalculated after subsetting to avoid confusion. 
#'The bracket operators can be alternatively, following the same syntax.
#'
#'Sample and snp facets to subset over can be provided. Filtering by facet
#'categories is performed by naming the facet as an argument name, then
#'providing the levels to keep or remove to that argument. See examples. Facets
#'designated as described in \code{\link{Facets_in_snpR}}.
#'
#'@param x snpRdata object.
#'@param .snps numeric, default \code{1:nrow(x)}. Row numbers corresponding to
#'  SNPs to keep. Negative subscripts allowed.
#'@param .samps numeric, default \code{1:ncol(x)}. Column numbers corresponding
#'  to samples to keep. Negative subscripts allowed.
#'@param i numeric, default \code{1:nrow(x)}. Row numbers corresponding to
#'  SNPs to keep. Negative subscripts allowed.
#'@param j numeric, default \code{1:nrow(x)}. Row numbers corresponding to
#'  SNPs to keep. Negative subscripts allowed.
#'@param ... Facet subsetting may be specified by providing the facet as an 
#'  argument and then providing the levels to keep or remove. Setting
#'  pop = "ASP", for example, would keep only samples in the "ASP" level of the
#'  pop facet, and setting pop.fam = "ASP.A" would keep only samples the ASP
#'  pop and the A fam. Negative subscripts are allowed: pop = -c("ASP", "PAL")
#'  would remove samples in the ASP or PAL pops. Subsetting by 
#'  multiple facets is supported, although negative and positive subscripts
#'  cannot be mixed across sample or SNP facets. They may be mixed between the
#'  two.
#'@param drop logical, default FALSE. Deprecated.
#'
#'@name subset_snpRdata
#'@aliases subset_snpR_data [
#'
#' @examples
#' # Keep only individuals in the ASP and PAL populations 
#' # and on the LGIX or LGIV chromosome.
#' subset_snpR_data(stickSNPs, pop = c("ASP", "PAL"), 
#'                  chr = c("groupIX", "groupIV"))
#'                  
#' # keep individuals/SNPs in the first 10 columns/rows.
#' subset_snpR_data(stickSNPs, 1:10, 1:10)
#' 
#' # negative subscripts: remove individuals in family B
#' subset_snpR_data(stickSNPs, fam = -"B")
#' 
#' # negative subscripts: remove individuals in pop PAL AND fam B or C, 
#' # and keep only SNPs on LGIV
#' subset_snpR_data(stickSNPs, pop.fam = -c("PAL.B", "PAL.C"),
#'                  chr = "groupIV")
#'                  
#' # using the bracket operator, same as example 1
#' stickSNPs[pop = c("ASP", "PAL"), chr = c("groupIX", "groupIV")]
#' 
#' # bracket operator, excluding first ten SNPs and only keeping pop = PAL,
#' # fam = B
#' stickSNPs[-c(1:10), pop.fam = c("PAL.B")]
NULL


#' @export
#' @describeIn subset_snpRdata subset_snpR_data
subset_snpR_data <- function(x, .snps = 1:nsnps(x), .samps = 1:nsamps(x), ...){
  #============extract facet info===============
  .argnames <- match.call()
  .argnames <- as.list(.argnames)
  .argnames <- .argnames[-1]
  
  .is.facet <- which(!names(.argnames) %in% c("x", ".snps", ".samps"))
  if(length(.is.facet) > 0){
    facets <- names(.argnames)[.is.facet]
    
    # eval the facets
    ## check for negatives, fix, and note
    .has.negatives <- grep("-", .argnames[.is.facet])
    .argnames[.is.facet][.has.negatives] <- lapply(.argnames[.is.facet][.has.negatives], function(x){
      x <- gsub("-", "", x)
      if(any(length(x) > 2)){stop("Negative subscripts cannot be provided inside c() calls. Please place negative subscripts outside c() calls.\n")}
      return(x[2])
    } )
    .single.part <- !grepl("^c\\(", .argnames[.is.facet][.has.negatives])
    .argnames[.is.facet][.has.negatives][!.single.part] <- lapply(.argnames[.is.facet][.has.negatives][!.single.part], str2lang)
    ## eval
    for(.candidate_facets in 1:length(.is.facet)){
      if(any(names(.argnames[.is.facet][.candidate_facets]) %in% c("facet", "subfacet", "facets", "subfacets", "snp.facet", "snp.subfacet", "snp.subfacets", "snp.facets"))){
        stop("Facets and subfacets are now desginated directly using, for example, pop = c('my.pop1', 'my.pop2'), not using the 'facet', 'subfacet', etc arguments.\n")
      }
      
      # the eval won't work if a string is given without "c()", so fix that
      if(!grepl("c\\(", .argnames[.is.facet][.candidate_facets])){
        .argnames[.is.facet][.candidate_facets] <- paste0("c(\"", .argnames[.is.facet][.candidate_facets], "\")")
      }
      .argnames[.is.facet][.candidate_facets] <- list(eval(parse(text = .argnames[.is.facet][.candidate_facets])))
    }
    
    
  }
  else{
    facets <- NULL
  }
  
  #===========sanity checks====================
  if(!is.snpRdata(x)){
    stop("x must be a snpRdata object.\n")
  }
  
  msg <- character(0)
  
  # check facet types
  if(!is.null(facets)){
    if(any(facets %in% c("facet", "subfacet", "facets", "subfacets", "snp.facet", "snp.subfacet", "snp.subfacets", "snp.facets"))){
      stop("Facets and subfacets are now desginated directly using, for example, pop = c('my.pop1', 'my.pop2'), not using the 'facet', 'subfacet', etc arguments.\n")
    }
    facets <- .check.snpR.facet.request(x, facets, "none", TRUE)
    if(any(c("sample", "complex") %in% facets[[2]]) & !identical(.samps, 1:nsamps(x))){
      msg <- c(msg, "Sample level facets cannot be provided alongside a vector of samples to retain/remove.\n")
    }
    if(any(c("snp", "complex") %in% facets[[2]]) %in% facets[[2]] & !identical(.snps, 1:nsnps(x))){
      msg <- c(msg, "SNP level facets cannot be provided alongside a vector of SNPs to retain/remove.\n")
    }
    if("complex" %in% facets[[2]]){
      msg <- c(msg, "Complex (sample + SNP) facets cannot be used for subsetting. Please split into SNP and sample components and rerun.\n")
    }
    # if(length(facets[[1]]) > 1){
    #   valid <- logical(length(facets[[1]]))
    #   for(i in 1:length(facets[[1]])){
    #     valid[i] <- ifelse(length(grep(facets[[1]][i], facets[[1]][-i])) == 0, TRUE, FALSE)
    #   }
    #   
    #   if(sum(valid) < length(valid)){
    #     msg <- c(msg, paste0("Some facets are listed more than once amongst arguments: ", paste0(facets[[1]][which(!valid)], collapse = ", "),
    #                          ". Each facet may only be listed once. Joint facets (such as pop.fam) will return only cases true in both, seperate facets will return cases true in either.\n"))
    #   }
    # }
    
  }
  
  # .snps and .samps
  if(is.logical(.snps)){
    if(length(.snps) != nsnps(x)){
      if(length(.snps) > nsnps(x)){
        warning("More .snps/j requested than present in data")
      }
      if(length(.snps) %% nsnps != 0){
        warning("Length of .snps/i is not a multiple of number of SNPs in x.\n")
      }
      .snps <- rep(.snps, length.out = nsnps(x))
    }
    .snps <- which(.snps)
  }
  if(is.logical(.samps)){
    if(length(.samps) != nsamps(x)){
      if(length(.samps) > nsamps(x)){
        warning("More .samps/j requested than present in data")
      }
      if(length(.samps) %% nsamps != 0){
        warning("Length of .samps/i is not a multiple of number of samples in x.\n")
      }
      .samps <- rep(.samps, length.out = nsamps(x))
    }
    .samps <- which(.samps)
  }
  
  if(!is.numeric(.snps) & !is.integer(.snps)){
    msg <- c(msg, ".snps/i must be a numeric vector.\n")
  }
  else{
    if(!all(.snps == as.integer(.snps))){
      msg <- c(msg, ".snps/i must only contain integers.\n")
    }
    if(max(.snps) > nrow(x)){
      msg <- c(msg, "All requested snps must be within 1:nsnps(x).\n")
    }
  }
  
  if(!is.numeric(.samps) & !is.integer(.samps)){
    msg <- c(msg, ".samps/j must be a numeric vector.\n")
  }
  else{
    if(!all(.samps == as.integer(.samps))){
      msg <- c(msg, ".samps/j must only contain integers.\n")
    }
    if(max(.samps) > ncol(x)){
      msg <- c(msg, "All requested snps must be within 1:nsnps(x).\n")
    }
  }

  
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  
  #===========subset===========================
  if(!is.null(facets)){
    
    # check to see if there are mixed negative and positive subsets
    if(any(facets[[2]] == "snp")){
      if(!all(which(facets[[2]] == "snp") %in% .has.negatives) & !all(!which(facets[[2]] == "snp") %in% .has.negatives)){
        stop("Cannot mix negative and positive subscripts in SNP subfacets.\n")
      }
    }
    if(any(facets[[2]] == "sample")){
      if(!all(which(facets[[2]] == "sample") %in% .has.negatives) & !all(!which(facets[[2]] == "sample") %in% .has.negatives)){
        stop("Cannot mix negative and positive subscripts in sample subfacets.\n")
      }
    }
    
    
    
    # check for flipped facet orders
    ord.flipped <- which(facets[[1]] != names(.argnames[.is.facet]))
    ## if any flipped, need to flip subfacets as well
    if(any(ord.flipped)){
      for(i in 1:length(ord.flipped)){
        split.ord.facet <- unlist(.split.facet(facets[[1]][i]))
        split.raw.facet <- unlist(.split.facet(names(.argnames[.is.facet])[i]))
        re_ord <- match(split.raw.facet, split.ord.facet)
        for(j in 1:length(.argnames[.is.facet][[i]])){
          split.lev <-  unlist(.split.facet(.argnames[.is.facet][[i]][j]))
          .argnames[.is.facet][[i]][j] <- paste0(split.lev[re_ord], collapse = ".")
        }
      }
    }
  }
  
  #====================snps====================
  if("snp" %in% facets[[2]]){
    if(any(which(facets[[2]] == "snp") %in% .has.negatives)){
      .snps <- rep(TRUE, nsnps(x))
    }
    else{
      .snps <- logical(nsnps(x))
    }
    
    snp.facets <- facets[[1]][which(facets[[2]] == "snp")]
    for(i in 1:length(snp.facets)){
      snp.subfacets <- .argnames[[.is.facet[which(facets[[2]] == "snp")][i]]]
      for(j in 1:length(snp.subfacets)){
        new.matches <- .fetch.snp.meta.matching.task.list(x, c(NA, NA, snp.facets[i], snp.subfacets[j]))
        if(length(new.matches) == 0){
          stop(paste0("No snp found matching: ", snp.facets[i], " -- ", snp.subfacets[j], "\n"))
          
        }
        if(which(facets[[2]] == "snp")[i] %in% .has.negatives){
          .snps[new.matches] <- FALSE
        }
        else{
          .snps[new.matches] <- TRUE
        }
      }
    }
    .snps <- which(.snps)
  }
  
  
  
  if("sample" %in% facets[[2]]){
    
    if(any(which(facets[[2]] == "sample") %in% .has.negatives)){
      .samps <- rep(TRUE, nsamps(x))
    }
    else{
      .samps <- logical(nsamps(x))
    }
    
    
    samp.facets <- facets[[1]][which(facets[[2]] == "sample")]
    for(i in 1:length(samp.facets)){
      samp.subfacets <- .argnames[[.is.facet[which(facets[[2]] == "sample")][i]]]
      for(j in 1:length(samp.subfacets)){
        new.matches <- .fetch.sample.meta.matching.task.list(x, c(samp.facets[i], samp.subfacets[j], NA, NA))
        if(length(new.matches) == 0){
          stop(paste0("No sample found matching: ", samp.facets[i], " -- ", samp.subfacets[j], "\n"))
        }
        if(which(facets[[2]] == "sample")[i] %in% .has.negatives){
          .samps[new.matches] <- FALSE
        }
        else{
          .samps[new.matches] <- TRUE
        }
      }
    }
    
    .samps <- which(.samps)
  }
  
  #===============return================
  nsampm <- sample.meta(x)[.samps,]

  nsnpm <- snp.meta(x)[.snps,]
  return(import.snpR.data(genotypes(x)[.snps, .samps, drop = FALSE], snp.meta = nsnpm,
                          sample.meta = nsampm, mDat = x@mDat))
}




# Old subset version with different facet specification. For internal use.
.subset_snpR_data <- function(x, snps = 1:nrow(x), samps = 1:ncol(x), facets = NULL, subfacets = NULL, snp.facets = NULL, snp.subfacets = NULL){
  #=========sanity checks========
  if(!is.snpRdata(x)){
    stop("x must be a snpRdata object.\n")
  }
  
  msg <- character(0)
  
  if(max(snps) > nrow(x)){
    msg <- c(msg, "All requested snps must be within 1:nrow(x).\n")
  }
  if(max(samps) > ncol(x)){
    msg <- c(msg, "All requested samps must be within 1:ncol(x).\n")
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  #=========subfunctions=========
  fix.for.one.snp <- function(x){
    if(nrow(x) == 1){
      # fix geno tables
      x@geno.tables <- lapply(x@geno.tables, FUN = function(y){
        a <- matrix(y, nrow = 1)
        colnames(a) <- names(y)
        return(a)}
      )
    }

    return(x)
  }

  #=========run subset===========

  # if subfacets or snp.subfacets were selected, figure out which samples and loci to keep
  if(!(is.null(snp.facets[1])) & !(is.null(snp.subfacets[1])) | !(is.null(facets[1])) & !(is.null(subfacets[1]))){

    # if snp.subfacets are requested
    if(!(is.null(snp.facets[1])) & !(is.null(snp.subfacets[1]))){
      if(!(any(snp.subfacets == ".base"))){
        t.snp.meta <- x@snp.meta

        # check for and get info on complex facets
        complex.snp.facets <- snp.facets[grep("(?<!^)\\.", snp.facets, perl = T)]
        if(length(complex.snp.facets) > 0){
          for(i in 1:length(complex.snp.facets)){
            tfacets <- unlist(.split.facet(complex.snp.facets[i]))
            tcols <- t.snp.meta[colnames(t.snp.meta) %in% tfacets]
            t.snp.meta <- cbind(t.snp.meta, .paste.by.facet(tcols, match(colnames(tcols), tfacets)))
          }
          colnames(t.snp.meta)[(ncol(t.snp.meta) - length(complex.snp.facets) + 1):ncol(t.snp.meta)] <- complex.snp.facets
        }

        # get the snps to keep
        t.snp.meta <- t.snp.meta[,colnames(t.snp.meta) %in% snp.facets]
        fsnps <- which(as.logical(rowSums(matrix(as.matrix(t.snp.meta) %in% snp.subfacets, nrow(x@snp.meta))))) # here's the snps to keep, those where at least one subfacet level is present in the provided snp.subfacets.
        snps <- snps[snps %in% fsnps]
      }
    }

    # if sample subfacets are requested
    if(!(is.null(facets[1])) & !(is.null(subfacets[1]))){
      if(!any(subfacets == ".base")){
        t.samp.meta <- x@sample.meta

        # check for and get info on complex facets
        complex.samp.facets <- facets[grep("(?<!^)\\.", facets, perl = T)]
        if(length(complex.samp.facets) > 0){
          for(i in 1:length(complex.samp.facets)){
            tfacets <- unlist(.split.facet(complex.samp.facets[i]))
            tcols <- t.samp.meta[colnames(t.samp.meta) %in% tfacets]
            t.samp.meta <- cbind(t.samp.meta, .paste.by.facet(tcols, match(colnames(tcols), tfacets)))
            
          }
          colnames(t.samp.meta)[(ncol(t.samp.meta) - length(complex.samp.facets) + 1):ncol(t.samp.meta)] <- complex.samp.facets
        }

        # get the samples to keep
        t.samp.meta <- t.samp.meta[,colnames(t.samp.meta) %in% facets]
        fsamps <- which(as.logical(rowSums(matrix(as.matrix(t.samp.meta) %in% subfacets, nrow(x@sample.meta))))) # here's the samples to keep, those where at least one subfacet level is present in the provided sample subfacets.
        samps <- samps[samps %in% fsamps]
      }
    }
  }

  # sort snps and samps according to sample/snp id
  snps <- snps[order(x@snp.meta$.snp.id[snps])]
  samps <- samps[order(x@sample.meta$.sample.id[samps])]

  # subset
  if(!identical(samps, 1:ncol(x))){
    dat <- genotypes(x)[snps, samps]
    if(length(samps) == 1){
      dat <- as.data.frame(dat, stringsAsFactors = F)
    }
    dat <- import.snpR.data(dat, snp.meta = snp.meta(x)[snps,], sample.meta = sample.meta(x)[samps,], mDat = x@mDat)
    if(any(x@facets != ".base")){
      dat <- .add.facets.snpR.data(dat, x@facets[-which(x@facets == ".base")])
    }

    if(length(x@sn$sn) != 0){
      sn <- x@sn$sn[,-c(1:(ncol(x@snp.meta) - 1))]
      sn <- sn[snps, samps]
      sn <- cbind(dat@snp.meta[,-ncol(dat@snp.meta)], sn)
      dat@sn <- list(type = x@sn$type, sn = sn)
    }

    warning("Since samples were subset, any stats will need to be recalculated.\n")
    return(dat)
  }
  else{
    if(length(x@sn) != 0){
      sn <- x@sn$sn[,-c(1:(ncol(x@snp.meta) - 1))]
      sn <- sn[snps,]
      sn <- cbind(snp.meta(x)[snps,-ncol(x@snp.meta)], sn)
      sn <- list(type = x@sn$type, sn = sn)
    }
    else{
      sn <- list()
    }

    # change snps and samps to snp and sample IDs
    keep.rows <- snps
    snps <- snp.meta(x)$.snp.id[snps]

    x <- snpRdata(.Data = genotypes(x)[which(snp.meta(x)$.snp.id %in% snps),],
                  sample.meta = sample.meta(x),
                  snp.meta = snp.meta(x)[which(snp.meta(x)$.snp.id %in% snps),],
                  facet.meta = x@facet.meta[x@facet.meta$.snp.id %in% snps,],
                  mDat = x@mDat,
                  snp.form = x@snp.form,
                  geno.tables = list(gs = x@geno.tables$gs[x@facet.meta$.snp.id %in% snps,],
                                     as = x@geno.tables$as[x@facet.meta$.snp.id %in% snps,],
                                     wm = x@geno.tables$wm[x@facet.meta$.snp.id %in% snps,]),
                  ac = x@ac[x@facet.meta$.snp.id %in% snps,],
                  facets = x@facets,
                  facet.type = x@facet.type,
                  stats = x@stats[x@stats$.snp.id %in% snps,],
                  pairwise.stats = x@pairwise.stats[x@pairwise.stats$.snp.id %in% snps,],
                  sn = sn,
                  names = x@names,
                  row.names = x@row.names[keep.rows])

    x <- fix.for.one.snp(x)



    warning("Any window stats will need to be recalculated.\n")
    return(x)
  }
}




#'Filter SNPs in snpRdata objects.
#'
#'\code{filter_snps} filters snpRdata objects to remove SNPs or individuals
#'which fail to pass user defined thresholds for several statistics. Since this
#'function removes all calculated statistics, etc. from the snpRdata object,
#'this should usually be the first step in an analysis. See details for filters.
#'
#'
#'
#'Possible filters: \itemize{ \item{maf, minor allele frequency: }{removes SNPs
#'where the minor allele frequency is too low. Can look for mafs below
#'#'provided either globally or search each population individually.}
#'\item{hf_hets, high observed heterozygosity: }{removes SNPs where the observed
#'heterozygosity is too high.} \item{min_ind, minimum individuals: }{removes
#'SNPs that were genotyped in too few individuals.} \item{min_loci, minimum
#'loci: }{removes individuals sequenced at too few loci.} \item{non_poly,
#'non-polymorphic SNPs: }{removes SNPs that are not polymorphic (not true
#'SNPs).} \item{bi_al, non-biallelic SNPs: }{removes SNPs that have more than
#'two observed alleles. This is mostly an internal argument, since the various
#'snpRdata import options use it automatically to prevent downstream errors in
#'other snpR functions. } }
#'
#'Note that filtering out poorly sequenced individuals creates a possible
#'conflict with the loci filters, since after individuals are removed, some loci
#'may no longer pass filters. For example, if a portion of individuals in one
#'population all carry the only instances of a rare minor allele that still
#'passes the maf threshold, removing those individuals may cause the loci to no
#'longer be polymorphic in the sample.
#'
#'To counter this, the "re_run" argument can be used to pass the data through a
#'second filtering step after individuals are removed. By default, the "partial"
#'re-run option is used, which re-runs only the non-polymorphic filter (if it
#'was originally set). The "full" option re-runs all set filters. Note that
#'re-running any of these filters may cause individuals to fail the individual
#'filter after loci removal, and so subsequent tertiary re-running of the
#'individual filters, followed by the loci filters, and so on, could be
#'justified. This is not done automatically here.
#'
#'Via the "maf_facets" argument, this function can filter by minor allele
#'frequencies in either \emph{all} samples or \emph{each level of a supplied
#'sample specific facet and the entire dataset}. In the latter case, any SNPs
#'that pass the maf filter in \emph{any} facet level are considered to pass the
#'filter. The latter should be used in instances where population sizes are very
#'different, there are \emph{many} populations, and/or allele frequencies are
#'very different between populations and thus common alleles of interest in one
#'population might be otherwise filtered out.
#'
#'The "hwe_facets" argument is the inverse of this: loci will be removed if they
#'fail the provided hwe filter in any facet level. In both cases, Facets should
#'be provided as described in \code{\link{Facets_in_snpR}}.
#'
#'
#'@param x snpRdata object.
#'@param maf numeric between 0 and 1 or FALSE, default FALSE. Minimum acceptable
#'  minor allele frequency.
#'@param hf_hets numeric between 0 and 1 or FALSE, default FALSE. Maximum
#'  acceptable heterozygote frequency.
#'@param hwe numeric between 0 and 1 or FALSE, default FALSE. Minimum acceptable
#'  HWE p-value.
#'@param min_ind numeric between 0 and 1 or FALSE, default FALSE. Minimum
#'  proportion of individuals in which a loci must be sequenced.
#'@param min_loci numeric between 0 and 1 or FALSE, default FALSE. Minimum
#'  proportion of SNPs at which an individual must be genotyped.
#'@param re_run character or FALSE, default "partial". When individuals are
#'  removed via min_ind, it is possible that some SNPs that initially passed
#'  filtering steps will now violate some filters. SNP filters can be re-run
#'  automatically via several methods: \itemize{ \item{partial: } Refilters for
#'  non-polymorphic loci (non_poly) only, if that filter was requested
#'  initially. \item{full: } Re-runs the full filtering scheme (save for
#'  min_loci).}
#'@param maf_facets character or FALSE, default FALSE. Defines a sample facet
#'  overwhich the minor allele frequency can be checked. SNPs will only fail the
#'  maf filter if they fail in every level of every provided facet.
#'@param hwe_facets character or FALSE, default FALSE. Defines a sample facet
#'  overwhich the hwe filter can be checked. SNPs will fail the hwe filter if
#'  they fail in any level of any provided facet.
#'@param non_poly logical, default TRUE. If TRUE, non-polymorphic loci will be
#'  removed.
#'@param bi_al logical, default TRUE. If TRUE, loci with more than two alleles
#'  will be removed. Note that this is mostly an internal argument and should
#'  rarely be used directly, since import.snpR.data and other snpRdata object
#'  creation functions all pass SNPs through this filter because many snpR
#'  functions will fail to work if there are more than two alleles at a locus.
#'
#'@return A data.frame in the same format as the input, with SNPs and
#'  individuals not passing the filters removed.
#'
#'@export
#'@author William Hemstrom
#'  
#' @examples
#' # Filter with a minor allele frequency of 0.05, maximum heterozygote 
#' # frequency of 0.55, 50% minimum individuals, and at least 75% of loci 
#' # sequenced per individual.
#' filter_snps(stickSNPs, maf = 0.05, hf_hets = 0.55, 
#'             min_ind = 0.5, min_loci = 0.75)
#'
#' # The same filters, but with minor allele frequency considered per-population
#' # and a full re-run of loci filters after individual removal.
#' filter_snps(stickSNPs, maf = 0.05, hf_hets = 0.55, min_ind = 0.5, 
#'             min_loci = 0.75, re_run = "full", maf_facets = "pop")
#'
filter_snps <- function(x, maf = FALSE, hf_hets = FALSE, hwe = FALSE, min_ind = FALSE,
                        min_loci = FALSE, re_run = "partial", maf_facets = NULL,
                        hwe_facets = FALSE,
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

  if(hwe){
    if(!is.numeric(hwe)){
      stop("hwe must be a numeric value.")
    }
    if(length(hwe) != 1){
      stop("hwe must be a numeric vector of length 1.")
    }
    if(hwe <= 0 | hwe >= 1){
      stop("hwe must be a value between 0 and 1.")
    }
    
    if(!isFALSE(hwe_facets)){
      hwe_facets <- .check.snpR.facet.request(x, hwe_facets)
    }
  }

  if(hf_hets){
    if(!is.numeric(hf_hets)){
      stop("hf_hets must be a numeric value.")
    }
    if(length(hf_hets) != 1){
      stop("hf_hets must be a numeric vector of length 1.")
    }
    if(hf_hets <= 0 | hf_hets  >= 1){
      stop("hf_hets must be a value between 0 and 1.")
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

  if(!is.null(maf_facets[1])){
    maf_facets <- .check.snpR.facet.request(x, maf_facets, "none")

    # add any needed facets...
    miss.facets <- maf_facets[which(!(maf_facets %in% x@facets))]
    if(length(miss.facets) != 0){
      cat("Adding missing facets...\n")
      # need to fix any multivariate facets (those with a .)
      x <- .add.facets.snpR.data(x, miss.facets)
    }

    # check for bad facets to remove (those that don't just consider samples)
    if(any(x@facet.type[x@facets %in% maf_facets] != "sample")){
      vio.facets <- x@facets[match(maf_facets, x@facets)]
      vio.facets <- vio.facets[which(x@facet.type[x@facets %in% maf_facets] != "sample")]
      warning(paste0("Facets over which to maf.filter must be sample specific facets, not snp specific facets! Removing non-sample facets: \n", paste0(vio.facets, collapse = " "), ".\n"))
      maf_facets <- maf_facets[-which(maf_facets %in% vio.facets)]
    }
  }

  #==============set up, get values used later, clean up data a bit,define subfunctions==========
  cat("Initializing...\n")

  #get headers
  headers <- x@snp.meta
  snp_form <- x@snp.form
  mDat <- x@mDat

  # fix a table if it only has one loci
  fix.one.loci <- function(x){
    if(is.null(nrow(x))){
      a <- matrix(x, nrow = 1)
      colnames(a) <- names(x)
      x <- a
    }
    return(x)
  }

  #function to filter by loci, to be called before and after min ind filtering (if that is requested.)
  filt_by_loci <- function(){
    # Store filter status in vio.snps. Those that are violating a filter will be marked TRUE, remove these.
    
    #==========================run filters========================
    vio.snps <- logical(nrow(x)) #vector to track status

    amat <- x@geno.tables$as[x@facet.meta$facet == ".base",]
    amat <- fix.one.loci(amat)
    gmat <- x@geno.tables$gs[x@facet.meta$facet == ".base",]
    gmat <- fix.one.loci(gmat)
    wmat <- x@geno.tables$wm[x@facet.meta$facet == ".base",]
    wmat <- fix.one.loci(wmat)

    # non-biallelic and non-polymorphic loci
    if(bi_al | non_poly){
      bimat <- ifelse(amat, TRUE, FALSE)

      if(bi_al){
        cat("Filtering non-biallelic loci...\n")
        bi <- ifelse(rowSums(bimat) > 2, T, F) # if false, should keep the allele
        cat(paste0("\t", sum(bi), " bad loci\n"))
        vio.snps[which(bi)] <- T
      }

      if(non_poly){
        cat("Filtering non_polymorphic loci...\n")
        np <- ifelse(rowSums(bimat) < 2, T, F) # if false, should keep the allele
        cat(paste0("\t", sum(np), " bad loci\n"))
        vio.snps[which(np)] <- T
      }
    }

    #========min inds=======
    if(min_ind){
      cat("Filtering loci sequenced in few individuals...\n")
      mi <- wmat[,colnames(wmat) == mDat]
      mi <- (nrow(x@sample.meta) - mi)/nrow(x@sample.meta) < min_ind
      vio.snps[which(mi)] <- T
      cat(paste0("\t", sum(mi), " bad loci\n"))
    }

    #========minor allele frequency, both total and by pop. Should only run if bi_al = TRUE.=========
    if(maf){
      #if not filtering with multiple pops
      if(is.null(maf_facets)){
        cat("Filtering low minor allele frequencies, no pops...\n")

        # check to see if we need to calculate mafs:
        if(any(colnames(x@stats) == "maf")){ # check that mafs have been calculated, the all facet must exist
          if(any(is.na(x@stats$maf[x@stats$facet == ".base"]))){ # check that mafs have been calculated for the all facet
            mafs <- 1 - matrixStats::rowMaxs(amat)/rowSums(amat)
          }
          else{
            mafs <- x@stats$maf[x@stats$facet == ".base"]
          }
        }
        else{
          mafs <- 1 - matrixStats::rowMaxs(amat)/rowSums(amat)
        }


        mafs <- mafs < maf #less than required, set to true and reject.
        mafs[is.na(mafs)] <- TRUE
        cat(paste0("\t", sum(mafs), " bad loci\n"))


        vio.snps[which(mafs)] <- T
      }
      else{
        cat("Filtering low minor allele frequencies by facet.\n")
        # pmafs <- logical(nrow(x))

        # see if we need to calculate mafs
        if(any(colnames(x@stats) == "maf")){ # mafs have been caluclated
          # get mafs for any uncalculated facets

          # if mafs have been calculated, but not for of the requested any facets...
          if(!any(x@stats$facet %in% maf_facets)){
            x <- calc_maf(x, facets = maf_facets)
          }
          #if mafs have been calculated for some but not all of our facets
          else if(any(is.na(x@stats$maf[x@stats$facet %in% maf_facets]))){
            run.facets <- unique(x@stats$facet[which(is.na(x@stats$maf[x@stats$facet %in% maf_facets]))])
            x <- calc_maf(x, facets = run.facets)
          }
        }
        else{
          x <- calc_maf(x, facets = maf_facets)
        }

        # grab mafs
        mafs <- x@stats[x@stats$facet %in% maf_facets,]
        mafs$maf <- mafs$maf < maf

        # now, figure out in how many subfacets the maf is too low. If all, the loci violates the filter
        cmafs <- reshape2::dcast(mafs[,c(".snp.id", "maf")], fun.aggregate = sum,
                                 formula = ... ~ "maf", value.var = "maf")

        # add in the overall maf, since differential fixation would otherwise be removed.
        if(any(is.na(x@stats$maf[x@stats$facet == ".base"]))){ # check that mafs have been calculated for the all facet
          a.mafs <- 1 - matrixStats::rowMaxs(amat)/rowSums(amat)
        }
        else{
          a.mafs <- x@stats$maf[x@stats$facet == ".base"]
        }

        a.mafs <- a.mafs < maf
        cmafs$maf <- cmafs$maf + a.mafs

        # check vio and report
        maf.vio <- which(cmafs$maf == (1 + length(unique(mafs$subfacet))))
        cat(paste0("\t", length(maf.vio), " bad loci\n"))
        vio.snps[maf.vio] <- T
      }
    }

    #========hf_hets. Should only run if bi_al = TRUE.==========
    if(hf_hets){
      cat("Filtering high frequency heterozygote loci...\n")

      # get heterozygote frequency
      hs <- which(substr(colnames(gmat), 1, snp_form/2) != substr(colnames(gmat), (snp_form/2) + 1, snp_form))
      het_f <- rowSums(gmat[,hs])/rowSums(gmat)

      # check violation
      het_f <- het_f > hf_hets #if false, heterozygote frequency is lower than cut-off, keep locus
      cat(paste0("\t", sum(het_f), " bad loci\n"))
      vio.snps[which(het_f)] <- T
    }

    #========hwe violation======================================
    if(hwe){
      cat("Filtering loci out of hwe...\n")
      
      # no facets, easy
      if(isFALSE(hwe_facets)){
        
        if(!.check_calced_stats(x, ".base", "hwe")$.base){
          invisible(utils::capture.output(x <- calc_hwe(x)))
        }
        phwe <- x@stats$pHWE[x@stats$facet == ".base"]
        phwe <- which(phwe < hwe)
        cat("\t", length(phwe), " bad loci\n")
        vio.snps[phwe] <- T
        
      }
      
      # facets, slightly more complicated
      else{
        
        run.facets <- .check_calced_stats(x, hwe_facets, "hwe")
        run.facets <- names(run.facets)[!unlist(run.facets)]
        if(length(run.facets) > 0){
          invisible(utils::capture.output(x <- calc_hwe(x, run.facets)))
        }
        
        # get the per-facet hwe stats, check against threshold, then condense by snp.
        phwe <- get.snpR.stats(x, hwe_facets)
        phwe$low.p <- ifelse(phwe$pHWE <= hwe, 1, 0)
        bad.loci <- tapply(phwe$low.p, phwe[,".snp.id"], sum, na.rm = TRUE)
        bad.loci <- reshape2::melt(bad.loci)
        bad.loci <- na.omit(bad.loci)
        bad.loci <- bad.loci[which(bad.loci$value > 0),]
        bad.loci <- which(snp.meta(x)$.snp.id %in% bad.loci$Var1)
        cat("\t", length(bad.loci), " bad loci\n")
        
        vio.snps[bad.loci] <- TRUE
      }
    }

    #==========remove violating loci==================
    if(any(vio.snps)){

      x <- x[-which(vio.snps),]
    }
    return(x)
  }

  #funciton to filter by individuals.
  min_loci_filt <- function(){
    cat("Filtering out individuals sequenced in few kept loci...\n")
    mcounts <- matrixStats::colSums2(ifelse(x != mDat, 1, 0))
    rejects <- which(mcounts/nrow(x) < min_loci)
    if(length(rejects) > 0){
      old.facets <- x@facets
      x <- x[,-rejects]
      cat("Re-calculating and adding facets.\n")
      if(any(old.facets != ".base")){
        x <- .add.facets.snpR.data(x, old.facets[-which(old.facets == ".base")])
      }

      warning("Any calculated stats will be removed, since individuals were filtered out!\n")
    }
    return(list(x = x, rejects = rejects))
  }

  #==========================call the functions as requested.==================
  if(any(c(non_poly, bi_al, maf, hf_hets, min_ind) != FALSE)){
    cat("Filtering loci. Starting loci:", nrow(x), "\n")

    # run the filter
    x <- filt_by_loci()

    if(nrow(x) == 0 | is.null(nrow(x))){
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
      x <- x$x
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
          hwe <- FALSE
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

#'Re-format SNP data.
#'
#'\code{format_snps} reformats SNP data into a range of different possible
#'formats for use in snpR functions and elsewhere.
#'
#'While this function can accept a few non-snpRdata input formats, it will
#'reformat to a snpRdata object internally. As such, it takes a facets argument
#'that works identically to elsewhere in the package, as described in
#'\code{\link{Facets_in_snpR}}. This argument is only used for output formats
#'where facets are important, such as the genepop format.
#'
#'While this function can be used as an alternative to
#'\code{\link{import.snpR.data}} when the output argument is set to "snpRdata",
#'this is more complicated and not recommended. Instead, it is simpler
#'\code{\link{import.snpR.data}} or one of the wrappers in
#'\code{\link{snpR_import_wrappers}}. The option is kept for backwards
#'compatability and internal use.
#'
#'If non-snpRdata is supplied, SNP and sample metadata may be provided. SNP
#'metadata may either be provided in the first few columns of x, the number of
#'which is designated by the input_meta_columns argument, or in a data.frame
#'given as via the snp.meta argument. Sample metadata may be provided in a
#'data.frame via the sample.meta argument.
#'
#'Output format options: \itemize{ \item{ac: }{allele count format, allele
#'counts tabulated for all samples or within populations.} \item{genepop:
#'}{genepop format, genotypes stored as four numeric characters (e.g. "0101",
#'"0204"), transposed, and formatted for genepop. Rownames are individual IDs in
#'genepop format, colnames are SNP ids, matching the first metadata column in
#'input.} \item{structure: }{STRUCTURE format, two lines per individual: allele
#'calls stored as single character numeric (e.g. "1", "2"). Allele calls per
#'individual stored on two subsequent lines.} \item{0000: }{numeric genotype tab
#'format, genotypes stored as four numeric characters (e.g. "0101", "0204").}
#'\item{hapmap: }{Migrate-n hapmap, allele counts tabulated within populations,
#'in migrate-n hapmap format. Since this migrate-n implementation is iffy, this
#'probably shouldn't be used much.} \item{NN: }{character genotype tab format,
#'genotypes stored as actual base calls (e.g. "AA", "CT").} \item{pa: }{allele
#'presence/absence format, presence or absence of each possible allele at each
#'possible genotype noted. Interpolation possible, with missing data substituted
#'with allele freqency in all samples or each population.} \item{rafm: }{RAFM
#'format, two allele calls at each locus stored in subsequent columns, e.g.
#'locus1.1 locus1.2.} \item{faststructure: }{fastSTRUCTURE format, identical to
#'STRUCTURE format save with the addition of filler columns proceeding data such
#'that exactly 6 columns proceed data. These columns can be filled with metadata
#'if desired.} \item{dadi: }{dadi format SNP data format, requires two columns
#'named "ref" and "anc" with the flanking bases around the SNP, e.g. "ACT" where
#'the middle location is the A/C snp.} \item{plink: }{PLINK! binary input
#'format, requires columns named "position" and one matching the name designated
#'with the 'chr' argument, and may contain a column named "cM", "cm", or
#'"morgans", containing linkage group/chr, snp ID, position in bp, and distance
#'in cM in order to create .bim extended map file.} \item{sn: }{Single character
#'numeric format. Each genotype will be listed as 0, 1, or 2, corresponding to
#'0, 1, or 2 minor alleles. Can be interpolated to remove missing data with the
#''interpolate' argument.} \item{sequoia: }{sequoia format. Each genotype is
#'converted to 0/1/2/ or -9 (for missing values). Requires columns ID, Sex,
#'BirthYear (or instead of BirthYear - BY.min and BY.max) in sample metadata 
#'for running Sequoia. For more information see sequoia documentation.} \item{fasta: }
#'{fasta sequence format.} \item{vcf:}{Variant Call Format, a standard format 
#'for SNPs and other genomic variants. Genotypes are coded as 0/0, 0/1, 1/1, or ./. 
#'(for missing values), with a healthy serving of additional metadata but very 
#'little sample metadata.}\item{snpRdata: }{a snpRdata object.} }
#'
#'Note that for the "sn" format, the data can be interpolated to fill missing
#'data points, which is useful for PCA, genomic prediction, tSNE, and other
#'methods. To do so, specify interpolate = "af" to insert the expected number of
#'minor alleles given SNP allele frequency or "bernoulli" to do binomial draws
#'to determine the number of minor alleles at each missing data point, where the
#'probability of drawing a minor allele is equal to the minor allele frequency.
#'The expected number of minor alleles based on the later method is equal to the
#'interpolated value from the former, but the later allows for multiple runs to
#'determine the impact of stochastic draws and is generally prefered and
#'required for some downstream analysis. It is therefore the default. As a
#'slower but more accurate alternative to "af" interpolation, "iPCA" may be
#'selected. This an iterative PCA approach to interpolate based on SNP/SNP
#'covariance via \code{\link[missMDA]{imputePCA}}. If the ncp arugment is not
#'defined, the number of components used for interpolation will be estimated
#'using \code{\link[missMDA]{estim_ncpPCA}}. In this case, this method is much
#'slower than the other methods, especially for large datasets. Setting an ncp
#'of 2-5 generally results in reasonable inpterpolations without the time
#'constraint.
#'
#'Note also that for the plink format, a .bed binary file can be generated. If
#'the "plink" option is selected and an outfile is designated, R will generate a
#'".sh" shell file with the same name given in the outfile argument. Running
#'this file will create a plink.bed file.
#'
#'Example datasets in each format are available in \code{\link{stickFORMATs}} in
#'elements named for output options.
#'
#'Input formats: \itemize{ \item{NULL or snpRdata: }{snpRdata object, the
#'default.} \item{NN: }{SNP genotypes stored as actual base calls (e.g. "AA",
#'"CT").} \item{0000: }{SNP genotypes stored as four numeric characters (e.g.
#'"0101", "0204").} \item{snp_tab: }{SNP genotypes stored with genotypes in each
#'cell, but only a single nucleotide noted if homozygote and two nucleotides
#'seperated by a space if heterozygote (e.g. "T", "T G").} \item{sn: }{SNP
#'genotypes stored with genotypes in each cell as 0 (homozyogous allele 1), 1
#'(heterozygous), or 2 (homozyogus allele 2).} \item{ms: }{.ms file, as output
#'from the simulation program ms.}}
#'
#'
#'
#'@param x snpRdata object or data.frame. Input data, in any of the above listed
#'  input formats.
#'@param output Character, default "snpRdata". The desired output format. A
#'  snpRdata object by default.
#'@param facets Character or NULL, default NULL. Facets overwhich to break up
#'  data for some output formats, following the format described in
#'  \code{\link{Facets_in_snpR}}.
#'@param n_samp Integer or numeric vector, default NA. For structure or RAFM
#'  outputs. How many random loci should be selected? Can either be an integer
#'  or a numeric vector of loci to use.
#'@param interpolate Character or FALSE, default "bernoulli". If transforming to
#'  "sn" format, notes the interpolation method to be used to fill missing data.
#'  Options are "bernoulli", "af", "iPCA", or FALSE. See details.
#'@param outfile character vector, default FALSE. If a filepath is provided, a
#'  copy of the output will be saved to that location. For some output styles,
#'  such as genepop, additional lines will be added to the output to allow them
#'  to be immediately run on commonly used programs.
#'@param ped data.frame default NULL. Optional argument for the "plink" output
#'  format. A six column data frame containing Family ID, Individual ID, Paternal
#'  ID, Maternal ID, Sex, and Phenotype and one row per sample. If provided,
#'  outputs will contain information contained in ped. See plink documentation
#'  for more details.
#'@param input_format Character, default NULL. Format of x, by default a
#'  snpRdata object. See description for details.
#'@param input_meta_columns Numeric, default NULL. If x is not a snpRdata
#'  object, optionally specifies the number of metadata columns preceeding
#'  genotypes in x. See details for more information.
#'@param input_mDat Character, default "NN". If x is not a snpRdata object, the
#'  coding for missing \emph{genotypes} in x (typically "NN" or "0000").
#'@param sample.meta data.frame, default NULL. If x is not a snpRdata object,
#'  optionally specifies a data.frame containing meta data for each sample. See
#'  details for more information.
#'@param snp.meta data.frame, default NULL. If x is not a snpRdata object,
#'  optionally specifies a data.frame containing meta data for each SNP. See
#'  details for more information.
#'@param chr.length numeric, default NULL. Chromosome lengths, for ms input
#'  files. Note that a single value assumes that each chromosome is of equal
#'  length whereas a vector of values gives the length for each chromosome in
#'  order.
#'@param ncp numeric or NULL, default 2. Number of components to consider for
#'  iPCA sn format interpolations of missing data. If null, the optimum number
#'  will be estimated, with the maximum specified by ncp.max. This can be very
#'  slow.
#'@param ncp.max numeric, default 5. Maximum number of components to check for
#'  when determining the optimum number of components to use when interpolating
#'  sn data using the iPCA approach.
#'@param chr character, default "chr". Name of column containing chromosome
#'  information, for VCF or plink! output.
#'@param position character, default "position". Name of column containing
#'  position information, for VCF output.
#'@param phenotype character, default "phenotype". Optional name of column
#'  containing phenotype information, for plink! output.
#'
#'@return A data.frame or snpRdata object with data in the correct format. May
#'  also write a file to the specified path.
#'
#'@export
#'
#'@author William Hemstrom
#'@author Melissa Jones
#'
#' @examples
#' #import data to a snpRdata object
#' ## get sample meta data
#' sample_meta <- 
#'     data.frame(pop = substr(colnames(stickFORMATs$`0000`)[-c(1:4)], 1, 3), 
#'                fam = rep(c("A", "B", "C", "D"), 
#'                          length = ncol(stickFORMATs$`0000`) - 4), 
#'                stringsAsFactors = FALSE)
#' format_snps(stickFORMATs$`0000`, input_format = "0000", 
#'             input_meta_columns = 4, 
#' input_mDat = "0000", sample.meta = sample_meta)
#'
#' #allele count, seperated by the pop facet.
#' format_snps(stickSNPs, "ac", facets = "pop")
#'
#' #genepop:
#' format_snps(stickSNPs, "genepop")
#'
#' #STRUCTURE, subsetting out 100 random alleles:
#' format_snps(stickSNPs, "structure", n_samp = 100)
#'
#' #STRUCTURE, subseting out the first 100 alleles:
#' format_snps(stickSNPs, "structure", n_samp = 1:100)
#'
#' #fastSTRUCTURE
#' format_snps(stickSNPs, "faststructure")
#'
#' #numeric:
#' format_snps(stickSNPs, "0000")
#'
#' #hapmap for migrate-n:
#' format_snps(stickSNPs, "hapmap", facets = "pop")
#'
#' #character:
#' format_snps(stickSNPs, "NN")
#'
#' #presence/absence, SNP data:
#' format_snps(stickSNPs, "pa")
#'
#' #RAFM, taking only 100 random snps and seperating by pop
#' format_snps(stickSNPs, "rafm", facets = "pop", n_samp = 100)
#'
#' #dadi
#' ## add ref and anc snp meta data columns to stickSNPs
#' dat <- as.data.frame(stickSNPs)
#' snp.meta(dat) <- cbind(ref = "ATA", anc = "ACT", snp.meta(stickSNPs))
#' format_snps(dat, "dadi", facets = "pop")
#'
#' #PLINK! format, not run to avoid file creation
#' \dontrun{
#' format_snps(stickSNPs, "plink", outfile = "plink_out", chr = "chr")
#' 
#' 
#' #PLINK! format with provided ped
#' ped <- data.frame(fam = c(rep(1, 210), rep("FAM2", 210)), ind = 1:420, 
#'                   mat = 1:420, pat = 1:420, 
#'                   sex = sample(1:2, 420, replace = TRUE), 
#'                   pheno = sample(1:2, 420, replace = TRUE))
#' format_snps(stickSNPs, "plink", outfile = "plink_out", 
#'             ped = ped, chr = "chr")
#' #note that a column in the sample metadata containing phenotypic information
#' #can be provided to the "phenotype" argument if wished.
#'
#' }
#'
#' #Sequoia format
#' b <- sample.meta(stickSNPs)
#' b$ID <- 1:nrow(b)
#' b$Sex <- rep(c("F", "M", "U", "no", "j"), length.out=nrow(b))
#' b$BirthYear <- round(runif(n = nrow(b), 1,1))
#' a <- stickSNPs
#' b$ID <- paste0(sample.meta(a)$pop, sample.meta(a)$fam, sample.meta(a)$ID)
#' sample.meta(a) <- b
#' format_snps(x=a, output = "sequoia")
#'
#'
#' # VCF format
#' test <- format_snps(stickSNPs, "vcf", chr = "chr")
format_snps <- function(x, output = "snpRdata", facets = NULL, n_samp = NA,
                        interpolate = "bernoulli", outfile = FALSE,
                        ped = NULL, input_format = NULL,
                        input_meta_columns = NULL, input_mDat = NULL,
                        sample.meta = NULL, snp.meta = NULL, chr.length = NULL,
                        ncp = 2, ncp.max = 5, chr = "chr", position = "position",
                        phenotype = "phenotype"){

  #======================sanity checks================
  if(!is.null(input_format)){
    if(tolower(input_format) == "snprdata"){input_format <- NULL}
  }

  # check that a useable output format is given. keming
  output <- tolower(output) # convert to lower case.
  if(output == "nn"){output <- "NN"}
  pos_outs <- c("ac", "genepop", "structure", "0000", "hapmap", "NN", "pa",
                "rafm", "faststructure", "dadi", "plink", "sn", "snprdata",
                "colony","adegenet", "fasta", "lea", "sequoia", "vcf")
  if(!(output %in% pos_outs)){
    stop("Unaccepted output format specified. Check documentation for options.\n")
  }

  # check that provided snpRdata objects are in the correct format
  if(is.null(input_format)){
    if(class(x) != "snpRdata"){
      stop("If x is not a snpRdata object, provide input data format.\n")
    }
  }

  # do checks, print info
  if(is.null(input_format) & !is.null(facets[1])){
    facet.types <- x@facet.type[match(facets, x@facets)]
    snp.facets <- which(facet.types == "snp")
    both.facets <- which(facet.types == "both")
    sample.facets <- which(facet.types == "sample")
  }
  else{
    both.facets <- character()
    snp.facets <- both.facets
    sample.facets <- both.facets
    facet.types <- both.facets
  }

  if(output == "ac"){
    cat("Converting to allele count format.\n")
  }
  else if(output == "genepop"){
    cat("Converting to genepop format.\n")
    if(length(c(both.facets, snp.facets)) != 0){
      warning("Removing invalid facet types (snp or snp and sample specific).\n")
      facets <- facets[-c(snp.facets, both.facets)]
    }
  }
  else if(output == "structure" | output == "faststructure"){
    if(output == "structure"){cat("Converting to STRUCTURE format.\n")}
    else{cat("Converting to fastSTRUCTURE format.\n")}
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
    if(length(c(both.facets, snp.facets)) != 0){
      warning("Removing invalid facet types (snp or snp and sample specific).\n")
      facets <- facets[-c(snp.facets, both.facets)]
    }
    if(output == "faststructure"){
      if(length(facets) > 5){
        stop("Too many facets selected for fastSTRUCTURE. Limit 5.\n")
      }
    }
    else{
      if(length(facets) > 1){
        stop("Only one facet permitted for structure format. Complex, sample only facets are allowed.\n")
      }
    }
  }
  else if(output == "0000"){
    cat("Converting to numeric 2 character format.\n")
  }
  else if(output == "hapmap"){
    cat("Converting to migrate-N hapmap format.\n")
  }
  else if(output == "NN"){
    cat("Converting to NN format.\n")
  }
  else if(output == "pa"){
    cat("Converting to allele presence/absense format.\n")
    if(length(c(both.facets, snp.facets)) != 0){
      warning("Removing invalid facet types (snp or snp and sample specific).\n")
      facets <- facets[-c(snp.facets, both.facets)]
    }
  }
  else if(output == "rafm"){
    cat("Converting to RAFM format.\n")
    if(length(c(both.facets, snp.facets)) != 0){
      warning("Removing invalid facet types (snp or snp and sample specific).\n")
      facets <- facets[-c(snp.facets, both.facets)]
    }
  }
  else if(output == "dadi"){
    if(is.null(input_format)){
      if(!("ref" %in% colnames(x@snp.meta)) | !("anc" %in% colnames(x@snp.meta))){
        stop("Reference and ancestor/outgroup flanking bases required in snp metadata columns named 'ref' and 'anc', respecitvely")
      }
    }
    else {
      if(!("ref" %in% colnames(x)) | !("anc" %in% colnames(x))){
        stop("Reference and ancestor/outgroup flanking bases required in columns named 'ref' and 'anc', respecitvely")
      }
    }
    cat("Converting to dadi format...\n")
    if(length(c(both.facets, snp.facets)) != 0){
      warning("Removing invalid facet types (snp or snp and sample specific).\n")
      facets <- facets[-c(snp.facets, both.facets)]
    }
  }
  else if(output == "plink"){
    .check.installed("genio")
    if(is.null(input_format)){
      if(!all(c("position") %in% colnames(x@snp.meta)) | !any(chr %in% colnames(x@snp.meta))){
        stop("A column named position and one matching the argument 'chr' containing position in bp and chr/linkage group/scaffold must be present in snp metadata!")
      }
    }
    else if(!all(c("position") %in% colnames(x)) | !any(chr %in% colnames(x))){
      stop("A column named position and one matching the argument 'chr' containing position in bp and chr/linkage group/scaffold must be present in x!")
    }
    if(!is.null(ped)){
      if(!is.data.frame(ped)){
        stop("ped must be a six column data.frame containg Family ID, Individual ID, Paternal ID, Maternal ID, Sex, and Phenotype and one row per sample. See plink documentation.\n")
      }
      if(ncol(ped) != 6 | nrow(ped) != ncol(x)){
        stop("ped must be a six column data.frame containg Family ID, Individual ID, Paternal ID, Maternal ID, Sex, and Phenotype and one row per sample. See plink documentation.\n")
      }
    }
    cat("Converting to PLINK! binary format.\n\tThis relies on the genio package--please cite this for publication.\n")
    if(length(c(both.facets, snp.facets)) != 0){
      warning("Removing invalid facet types (snp or snp and sample specific).\n")
      facets <- facets[-c(snp.facets, both.facets)]
    }
  }
  else if(output == "sn"){
    cat("Converting to single character numeric format.\n")
    if(length(c(both.facets, snp.facets)) != 0){
      warning("Removing invalid facet types (snp or snp and sample specific).\n")
      facets <- facets[-c(snp.facets, both.facets)]
    }
  }
  else if(output == "snprdata"){
    if(is.null(input_format)){
      stop("Data already in snpRdata object.\n")
    }
    if(input_format != "ms"){
      if(all(is.null(snp.meta), is.null(input_meta_columns)) | is.null(input_mDat)){
        stop("sample meta, snp meta, and input missing data format must be provided for conversion to snpRdata object.")
      }
      else if(is.null(snp.meta) & !is.null(input_meta_columns)){
        cat("Using input metadata columns as snp meta.\n")
      }
    }
    else if(is.null(input_mDat)){
      input_mDat <- -1
    }

    cat("Converting to snpRdata object.\n")
  }
  else if(output == "colony"){
    cat("Converting to colony format.\n")
  }
  else if(output == "adegenet"){
    cat("Converting to adegenet genind object. SNP metadata will be discarded.\n")
  }
  else if(output == "fasta"){
    cat("Converting to psuedo-fasta file. All snps will be treated as a single sequence.\nHeterozygotes will be randomly called as either the major or minor allele.\n")
  }
  else if(output == "lea"){
    cat("Converting to LEA geno format. All metadata will be discarded.\n")
  }
  else if(output == "sequoia"){
    cat("Converting to Sequoia format.\n")
  }
  else if(output == "vcf"){
    if(!position %in% colnames(snp.meta(x))){
      stop(paste0("No column named matching position argument (", position, ") in SNP metadata. This is required for VCF output.\n"))
    }
    if(!chr %in% colnames(snp.meta(x))){
      stop(paste0("No column named matching chr argument (", chr, ") in SNP metadata. This is required for VCF output.\n"))
    }
    
    cat("Converting to VCF format.\n")
  }

  else{
    stop("Please specify output format.")
  }
  
  if(interpolate == "iPCA"){
    .check.installed("missMDA")
  }

  #======================put data into snpRdata object if not in that format to start with================
  if(!is.null(input_format)){
    cat("Converting data to snpRdata, NN format.\n")

    if(input_format == "ms"){
      if(!is.null(snp.meta) | !is.null(input_meta_columns)){
        input_meta_columns <- NULL
        snp.meta <- NULL
        warning("For ms inputs, provided snp.meta will be ignored and will be pulled from
              input ms instead.\n")
      }
      if(is.null(sample.meta)){
        stop("For ms input, please provide sample metadata.\n")
      }
      if(!is.numeric(chr.length)){
        stop("For ms input, please provide either a single or a vector of chromosome lengths.\n")
      }
      if(!is.character(x)){
        stop("For ms input, please provide a file path to the 'x' argument.\n")
      }
      else if(length(x) != 1){
        stop("For ms input, please provide a file path to the 'x' argument.\n")
      }
      else if(!file.exists(x)){
        stop("File provided to 'x' not found.\n")
      }

      cat("Parsing ms file...")
      x <- process_ms(x, chr.length)
      snp.meta <- x$meta
      x <- x$x
      x <- .convert_2_to_1_column(x)
      cat(" done.\n")

      input_format <- "sn"
    }

    if(!is.null(input_meta_columns)){
      headers <- x[,c(1:input_meta_columns)]
      x <- x[,-c(1:input_meta_columns)]
    }
    else{
      if(is.null(snp.meta)){
        headers <- data.frame(snpID = 1:nrow(x))
      }
      else{
        headers <- snp.meta
      }
    }
    if(is.null(sample.meta)){
      sample.meta <- data.frame(ID = 1:ncol(x))
    }

    #====================checks===================
    if(is.null(input_mDat)){
      stop("Input missing data encoding required.\n")
    }

    if(input_format == "NN"){
      cat("Input format: NN\n")
      if(nchar(input_mDat) != 2){
        stop("Missing data format must be two character.\n")
      }
      else{
        cat("Missing data format: ", input_mDat, "\n")
      }
    }
    else if(input_format == "0000"){
      cat("Input format: 0000\n")
      if(nchar(input_mDat) != 4){
        stop("Missing data format must be four characters.\n")
      }
      else{
        cat("Missing data format: ", input_mDat, "\n")
      }
    }
    else if(input_format == "snp_tab"){
      cat("Input format: snp_tab.\n")
      if(nchar(input_mDat) != 2){
        stop("Missing data format must be two characters.\n")
      }
      else{
        cat("Missing data format: ", input_mDat, "\n")
      }
    }
    else if(input_format == "sn"){
      cat("Iput format: sn\n")
      if(input_mDat %in% c(0:2)){
        stop("Missing data format must be other than 0, 1, or 2.\n")
      }
      else{
        cat("Missing data format: ", input_mDat, "\n")
      }
    }
    else{
      stop("Unsupported input format.")
    }

    #====================convert to snpRdata intermediate=============
    # 0000 put into snpRdata
    if(input_format == "0000"){
      #vectorize and replace
      xv <- as.matrix(x)
      xv <- gsub("01", "A", xv)
      xv <- gsub("02", "C", xv)
      xv <- gsub("03", "G", xv)
      xv <- gsub("04", "T", xv)
      xv <- gsub(substr(input_mDat, 1, 2), "N", xv)
      x <- as.data.frame(xv, stringsAsFactors = F)
      rm(xv)

      input_mDat <- "NN"

      #rebind to matrix and remake data.
      if(output == "NN"){
        rdata <- cbind(sample.meta, x) #all done if just converting to NN
      }
      else{
        x <- import.snpR.data(x, sample.meta = sample.meta, snp.meta = headers, mDat = "NN")
      }
    }

    # convert snp_tab to snpRdata
    if(input_format == "snp_tab"){
      xv <- as.matrix(x)
      ptf <- nchar(xv) == 1
      xv[ptf] <- paste0(xv[ptf], xv[ptf]) #double up homozygotes
      remove(ptf)
      xv <- gsub(" ", "", xv) #combine heterozygotes.
      xv[xv == input_mDat] <- "NN" #replace with correct missing data
      x <- as.data.frame(xv, stringsAsFactors = F)
      if(output == "NN"){
        rdata <- cbind(headers, x)
      }
      else{
        x <- import.snpR.data(x, sample.meta = sample.meta, snp.meta = headers, mDat = "NN")
      }
    }

    #convert "sn" format to NN
    if(input_format == "sn"){
      xv <- as.matrix(x)
      xv[xv == 0] <- "AA"
      xv[xv == 1] <- "AC"
      xv[xv == 2] <- "CC"
      xv[xv == input_mDat] <- "NN"
      x <- as.data.frame(xv, stringsAsFactors = F)

      if(output == "NN"){
        rdata <- as.data.frame(x) #all done if just converting to NN
      }
      else{
        x <- import.snpR.data(x, sample.meta = sample.meta, snp.meta = headers, mDat = "NN")
      }
    }

    if(input_format == "NN"){
      x <- import.snpR.data(x, sample.meta = sample.meta, snp.meta = headers, mDat = "NN")
    }

    if(!is.null(facets)){
      x <- .add.facets.snpR.data(x, facets)
    }

    # if this is the desired output, we're done.
    if(output == "snprdata"){
      if(outfile != FALSE){
        saveRDS(x, outfile)
      }
      return(x)
    }
  }

  #======================prepare outfile==================================================================
  if(outfile != FALSE){
    if(is.character(outfile) & length(outfile) == 1){
      cat("Printing results to:", outfile, "\n")
    }
    else{
      stop("Outfile must be a character vector of length 1.\n")
    }
  }

  #======================do conversions===================================================================
  # add missing facets
  if(length(facets) == 0){facets <- NULL}

  if(!is.null(facets[1])){
    if(!all(facets %in% x@facets)){
      cat("Adding missing facets.\n")
      new.facets <- facets[which(!(facets %in% x@facets))]
      x <- .add.facets.snpR.data(x, new.facets)
    }
  }

  # set facets to ".base" if not requested
  if(is.null(facets[1])){
    facets <- ".base"
  }

  # convert to allele count, migrate-n, or dadi format. Migrate-n should ALWAYS have multiple pops (why else would you use it?)
  if(output == "ac" | output == "hapmap" | output == "dadi"){
    if(output == "hapmap"){cat("WARNING: Data does not have header or pop spacer rows.\n")}

    #=========basic ac constructor function, to be on each facet==========
    # x: stats slot of a snpR object, filtered to the desired facets
    # maj: major allele identities across all facets at the relevent snps
    # min: minor allele identities across all facets at the relevent snps
    # mis.al: missing allele coding
    get.ac <- function(x, maj, min, mis.al){
      # initialize:
      if(is.null(nrow(x))){
        temp.x <- matrix(x, ncol = length(x))
        colnames(temp.x) <- names(x)
        x <- temp.x
      }

      out <- data.frame(n_total = numeric(length(maj)),
                        n_alleles = numeric(length(maj)),
                        ni1 = numeric(length(maj)),
                        ni2 = numeric(length(maj)))


      # get the column from as matching the target allele.
      maj.col.match <- match(maj, colnames(x))
      out$ni1 <- t(x)[maj.col.match + seq(0, length(x) - ncol(x), by = ncol(x))]

      # ni1 is the rowsums minus this
      out$ni2 <- rowSums(x) - out$ni1
      out$n_total <- rowSums(x)
      out$n_alleles <- rowSums(ifelse(out[,3:4] != 0, 1, 0))
      out[is.na(out)] <- 0

      return(out)
    }

    #=========apply to requested facets=======
    # get missing maf info
    mafs_to_calc <- .check_calced_stats(x, unique(c(".base", facets)), "maf")
    if(any(!unlist(mafs_to_calc))){
      x <- calc_maf(x, names(mafs_to_calc[which(!unlist(mafs_to_calc))]))
    }

    # grab the major and minors for each snp in the analysis.
    .snp.id <- facet <- subfacet <- NULL
    x@stats <- dplyr::arrange(x@stats, .snp.id, facet, subfacet)
    maj <- x@stats$major[x@stats$facet == ".base"]
    min <- x@stats$minor[x@stats$facet == ".base"]
    maj <- rep(maj, each = length(unique(x@facet.meta[x@facet.meta$facet %in% facets,]$subfacet)))
    min <- rep(min, each = length(unique(x@facet.meta[x@facet.meta$facet %in% facets,]$subfacet)))


    # get the ac file
    ac.dat <- get.ac(x@geno.tables$as[x@facet.meta$facet %in% facets,],
                     maj = maj,
                     min = min,
                     substr(x@mDat, 1, nchar(x@mDat)/2))

    if(output == "dadi"){
      rdata <- cbind(x@stats[x@stats$facet %in% facets, c("facet", "subfacet", ".snp.id")], ac.dat)
      ni1 <- reshape2::dcast(rdata, facet + .snp.id ~ subfacet, value.var = c("ni1"))
      ni2 <- reshape2::dcast(rdata, facet + .snp.id ~ subfacet, value.var = c("ni2"))

      rdata <- cbind(ref = x@snp.meta$ref, # since everything is sorted by .snp.id, this will match.
                     anc = x@snp.meta$anc,
                     Allele1 = x@stats[x@stats$facet == ".base", "major"],
                     ni1[order(ni1$.snp.id), 3:ncol(ni1), drop = FALSE],
                     Allele2 = x@stats[x@stats$facet == ".base", "minor"],
                     ni2[order(ni2$.snp.id), 3:ncol(ni2), drop = FALSE],
                     x@snp.meta[,which(!(colnames(x@snp.meta) %in% c("ref", "anc", ".snp.id")))])
      rdata <- as.data.frame(rdata)
      colnames(rdata)[c(3,3 + length(3:ncol(ni1)) + 1)] <- c("Allele1", "Allele2")
      
      if(length(facets) == 1 & facets[1] == ".base"){
        colnames(rdata)[c(4,6)] <- ".base"
      }
    }
    else{
      rdata <- cbind(x@facet.meta[x@facet.meta$facet %in% facets,],
                     ac.dat)
      if(output == "hapmap"){
        rdata <- dplyr::arrange(rdata, facet, subfacet, .snp.id)
      }
    }
  }

  if(output == "NN"){
    rdata <- cbind(x@snp.meta, as.data.frame(x))
    colnames(rdata)[1:ncol(x@snp.meta)] <- colnames(x@snp.meta)
  }


  ##convert to genepop or numeric format (v)
  if (output == "genepop" | output == "0000"){
    # keming
    #vectorize and replace
    xv <- as.matrix(x)
    xv <- gsub("A", "01", xv)
    xv <- gsub("C", "02", xv)
    xv <- gsub("G", "03", xv)
    xv <- gsub("T", "04", xv)
    xv <- gsub(substr(x@mDat, 1, 2), "0000", xv)

    # keming
    if(output == "genepop"){ #convert to genepop
      rdata <- as.data.frame(t(xv), stringsAsFactors = F) #remove extra columns and transpose data
      row.names(rdata) <- paste0(row.names(rdata), " ,") #adding space and comma to row names, as required.
    }
    # else if(output == "baps"){}
    else {#prepare numeric output, otherwise same format
      rdata <- as.data.frame(xv, stringsAsFactors = F)
      rdata <- cbind(x@snp.meta, rdata)
      colnames(rdata)[1:ncol(x@snp.meta)] <- colnames(x@snp.meta)
    }
  }


  ##convert to structure, fastStructure or RAFM format (v)
  if (output == "structure" | output == "rafm" | output == "faststructure"){
    #subset if requested
    if(all(!is.na(n_samp))){
      cat("Subsetting ")
      if(length(n_samp) > 1){
        cat("designated SNPs.\n")
        x <- subset_snpR_data(x, n_samp)
      }
      else{
        cat(n_samp, " random SNPs.\n")
        x <- subset_snpR_data(x, sample(1:nrow(x)))
      }
    }

    # transpose, since these are sample based, then gsub
    xv <- t(as.matrix(x))
    xv <- gsub("A", 1, xv)
    xv <- gsub("C", 2, xv)
    xv <- gsub("G", 3, xv)
    xv <- gsub("T", 4, xv)
    xv <- gsub(substr(x@mDat, 1, nchar(x@mDat)/2), 0, xv)

    #split into two matrices
    xv1 <- substr(xv, 1, 1)
    xv2 <- substr(xv, 2, 2)

    #create output matrix and bind the data to it, structure format
    if(output == "structure" | output == "faststructure"){
      outm <- matrix(NA, nrow = 2*(ncol(x)), ncol =  nrow(x))

      #fill
      outm[seq(1,nrow(outm),2),] <- xv1
      outm[seq(2,nrow(outm),2),] <- xv2

      #add sample names
      snames <- character(nrow(outm))
      snames[seq(1,nrow(outm),2)] <- x@names
      snames[seq(2,nrow(outm),2)] <- x@names
      
      # fix missing to -9
      outm[outm == 0] <- -9

      
      if(output == "structure"){ # bind sample and metadata for structure

        
        #prep pop numbers
        facets <- .check.snpR.facet.request(x, facets)
        facets <- unlist(.split.facet(facets))
        
        struc.meta <- as.numeric(as.factor(
          .paste.by.facet(sample.meta(x)[rep(1:nsamps(x), each = 2),], colnames(sample.meta(x)) %in% facets)))
        
        
        if(length(struc.meta) > 0){
          rdata <- cbind(ind = snames,
                       struc.meta,
                       as.data.frame(outm, stringsAsFactors = F))
        }
        else{
          rdata <- cbind(ind = snames,
                         as.data.frame(outm, stringsAsFactors = F))
        }
        
        colnames(rdata)[2] <- paste(facets, collapse = ".")
      }
      else{ #add a bunch of filler columns for faststructure and change missing data to -9
        if(length(facets) == 1 & facets[1] == ".base"){
          rdata <- cbind(ind = snames,
                         matrix("filler", nrow(outm), 5),
                         as.data.frame(outm, stringsAsFactors = F))
          colnames(rdata)[2:6] <- paste0("filler", 1:5)
        }
        else{
          n.valid.facets <- sum(colnames(x@sample.meta) %in% facets)
          rdata <- cbind(ind = snames,
                         x@sample.meta[,colnames(x@sample.meta) %in% facets],
                         matrix("filler", nrow(outm), 5 - n.valid.facets),
                         as.data.frame(outm, stringsAsFactors = F))
          colnames(rdata)[(1 + n.valid.facets):6] <- paste0("filler", 1:(5 - n.valid.facets))
        }
      }
    }

    #create output matrix and bind to it, RAFM format.
    else{
      outm <- matrix(NA, ncol(x), nrow(x)*2)

      #fill
      outm[,seq(1,ncol(outm),2)] <- xv1
      outm[,seq(2,ncol(outm),2)] <- xv2

      #replance missings with NA
      outm[outm == 0] <- NA

      #add column names
      colnames(outm) <- paste0("locus", sort(rep(1:nrow(x),2)), rep(c(".1", ".2"), ncol(outm)/2))

      #add subpop numbers, if given
      rdata <- cbind(x@sample.meta[,colnames(x@sample.meta) %in% facets],
                     as.data.frame(outm, stringsAsFactors = F))
      colnames(rdata)[1:sum(colnames(x@sample.meta) %in% facets)] <- colnames(x@sample.meta)[colnames(x@sample.meta) %in% facets]
    }
  }


  #presence/absence format
  if(output == "pa"){
    xv <- as.matrix(x)

    #function to produce vectorized presence absence (as much as possible, not vectorized for NAs)
    pa_alleles <- function(xv, snp_form, mDat, nsamp, nloci){
      nsamp <- nrow(xv)
      nloci <- ncol(xv)
      #get all possible genotypes
      gs <- unique(as.vector(xv))

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
      colnames(amat) <- paste0(sort(rep(1:nrow(x),length(as))), "_", as) #initialize all of the columns, with locus number followed by allele. Will remove anything empty after assignment.

      #fill in
      for(i in 1:length(as)){
        pr1 <- grep(as[i], xva1) # unique rows which have the allele in either position one or position two.
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
      xmc <- which(xv == mDat) #which samples had missing data?
      adj <- floor(xmc / nsamp) #how many loci over do I need to adjust xmc, since in amat each locus occupies two columns?
      adj[xmc%%nsamp == 0] <- adj[xmc%%nsamp == 0] - 1 #shift over anything that got screwed up by being in the last sample
      xmc <- xmc + (nsamp*adj) #adjust xmc for extra columns.
      if(any(amat[xmc] != 0) | any(amat[xmc + nsamp] != 0)){
        stop("Missing data values were not properly identified for replacement with NAs. This usually happens when SNP data is not completely bi-allelic. Try filtering out non-biallelic and non-polymorphic SNPs using filter_snps.\n")
      }
      amat[xmc] <- NA #make the first allele NA
      amat[xmc + nsamp] <- NA #make the second allele (another column over) NA.

      return(amat)
    }

    amat <- pa_alleles(t(xv), 2, x@mDat)

    # if(interp_miss){
    #   #average number observed in columns
    #   cat("Interpolating missing data...\n")
    #   afs <- colMeans(amat, TRUE)
    #   temp <- which(is.na(amat))/nrow(amat)
    #   fill_element <- floor(temp) + 1 #get the column for each missing data point
    #   fill_element[which(temp %% 1 == 0)] <- fill_element[which(temp %% 1 == 0)] - 1 #correct for anything in the last row
    #   amat[which(is.na(amat))] <- afs[fill_element] #fill with the appropriate allele frequency.
    # }
    # else{cat("Finished. Warning: Missing data counts are also stored!\n")}
    amat <- cbind(samp = as.character(colnames(x)), as.data.frame(amat, stringsAsFactors = F))
    rdata <- amat
  }


  #PLINK format
  if(output == "plink"){
    #===============make a bed file using genio==========
    # remove any snps with no data
    sn <- as.matrix(format_snps(x, "sn", interpolate = FALSE)[,-c(1:(ncol(snp.meta(x)) - 1))])
    bads <- which(rowSums(is.na(sn)) == ncol(sn))
    if(length(bads) > 0){
      sn <- sn[-bads,]
      warning("Removed some loci without any called genotypes.\n")
    }
    
    if(isFALSE(outfile)){
      genio::write_bed("plink_out.bed", sn, verbose = FALSE)
      cat("Wrote bed file to plink_out.bed\n")
    }
    genio::write_bed(paste0(outfile, ".bed"), sn, verbose = FALSE)

    #===============make a ped file=================
    # make an empty set of ped header columns if not provided
    lower.sample.cols <- tolower(colnames(x@sample.meta))
    if(is.null(ped)){
      if(any(lower.sample.cols == "fam")){
        fam <- x@sample.meta[,which(lower.sample.cols == "fam")]
      }
      else{
        fam <- rep("NA", ncol(x))
      }
      if(any(lower.sample.cols == "patid")){
        PatID <- x@sample.meta[,which(lower.sample.cols == "patid")]
      }
      else{
        PatID <- rep("NA", ncol(x))
      }
      if(any(lower.sample.cols == "sex")){
        Sex <- x@sample.meta[,which(lower.sample.cols == "sex")]
      }
      else{
        Sex <- rep("NA", ncol(x))
      }
      if(any(colnames(x@sample.meta) == phenotype)){
        cat("Using phenotype column as phenotype.")
        Phenotype <- x@sample.meta[,which(colnames(x@sample.meta) == phenotype)]
        if(length(unique(na.omit(Phenotype))) == 2){
          uf <- factor(Phenotype)
          cat("Two phenotypes detected, options encoded as: \n\t",
              levels(uf)[1], " = 0\n\t",
              levles(uf)[2], " = 1\n\t",
              "\n\tNA as -9.\n")
          uf <- as.numeric(uf)
          uf[is.na(uf)] <- -9
          Phenotype <- uf
        }
        else{
          warning("Plink does not typically support more than two phenotypes. Phenotypes left as-is, with missing data as NA (plink requires -9).\n")
        }
      }
      else{
        cat("No phenotype column in sample metadata detected, filling with missing.\n")
        Phenotype <- rep("-9", ncol(x))
      }
      if(any(lower.sample.cols == "matid")){
        MatID <- x@sample.meta[,which(lower.sample.cols == "matid")]
      }
      else{
        MatID <- rep("NA", ncol(x))
      }
      ped <- data.frame(fam = fam,
                        ind = colnames(x),
                        PatID = PatID,
                        MatID = MatID,
                        Sex = Sex,
                        Phenotype = Phenotype)
    }
    else{
      if(any(colnames(ped) == "Phenotype")){
        Phenotype <- ped$Phenotype
        if(length(unique(na.omit(Phenotype))) == 2){
          uf <- factor(Phenotype)
          cat("Two phenotypes detected, options encoded as: \n\t",
              levels(uf)[1], " = 0\n\t",
              levles(uf)[2], " = 1\n\t",
              "\n\tNA as -9.\n")
          uf <- as.numeric(uf)
          uf[is.na(uf)] <- -9
          Phenotype <- uf
          ped$Phenotype <- Phenotype
        }
        else{
          warning("Plink does not typically support more than two phenotypes. Phenotypes left as-is, with missing data as NA (plink requires -9).\n")
        }
      }
    }
    

    # save .fam
    fam <- ped

    # change missing data value and add a space between alleles.
    x.np <- as.vector(t(x))
    x.np[x.np == x@mDat] <- "00"
    x.np <- gsub("(.)(.)", "\\1 \\2", x.np)

    # rebind
    ped <- cbind(ped, matrix(x.np, nrow(ped), nrow(x)), stringsAsFactors = F)

    #===============make an extended map file=================
    a.names <- get.snpR.stats(x, stats = "maf")$single[,c("major", "minor")]
    # with morgans
    if(any(colnames(x@snp.meta) %in% c("cM", "cm", "morgans"))){
      bim <- data.frame(chr = x@snp.meta[,chr],
                        rs = x@snp.meta$.snp.id,
                        cM = x@snp.meta[,which(colnames(x@snp.meta) %in% c("cM", "cm", "morgans"))][,1],
                        bp = x@snp.meta$position,
                        a1 = a.names[,1],
                        a2 = a.names[,2])
    }
    # without morgans
    else{
      bim <- data.frame(chr = x@snp.meta[,chr],
                        rs = x@snp.meta$.snp.id,
                        cM = 0,
                        bp = x@snp.meta$position,
                        a1 = a.names[,1],
                        a2 = a.names[,2])
    }

    # recode chr
    bim$chr <- as.numeric(as.factor(bim$chr))

    # remove bads
    if(length(bads) > 0){
      bim <- bim[-bads,]
    }
    
    # grab normal map file
    map <- bim[,1:4]

    # name output
    rdata <- list(ped = ped, map = map, bim = bim, fam = fam)
  }

  # colony format for 01234
  if(output == "colony"){
    #gsub the values
    rdata <- gsub(pattern = "N", replacement = "0", x = as.matrix(x))
    rdata <- gsub(pattern = "A", replacement = "1", x = rdata)
    rdata <- gsub(pattern = "C", replacement = "2", x = rdata)
    rdata <- gsub(pattern = "G", replacement = "3", x = rdata)
    rdata <- gsub(pattern = "T", replacement = "4", x = rdata)
    rdata <- gsub(pattern = "^([0-4])([0-4])", replacement = "\\1 \\2", rdata)
    # transpose
    rdata <- t(rdata)
    rdata <- cbind(colnames(x), as.data.frame(rdata))

  }

  # single-character numeric format
  if(output == "sn" | output == "lea" | output == "sequoia"){

    # grab major/minor info via calc_maf, unless already calculated
    if(!any(colnames(x@stats) == "major")){
      x <- calc_maf(x)
    }
    else if(any(is.na(x@stats$major[x@stats$facet == ".base"]))){
      x <- calc_maf(x)
    }

    min <- x@stats$minor[x@stats$facet == ".base"]
    maj <- x@stats$major[x@stats$facet == ".base"]

    # check to see if each allele is the minor, assign a one if so
    a1 <- substr(unlist(t(genotypes(x))), 1, 1)
    a2 <- substr(unlist(t(genotypes(x))), 2, 2)

    a1 <- a1 == rep(min, each = ncol(x))
    a2 <- a2 == rep(min, each = ncol(x))

    rdata <- t(a1 + a2)
    rdata[as.matrix(genotypes(x)) == x@mDat] <- NA


    # sn
    if(output == "sn"){

      # grab out metadata
      meta <- x@snp.meta

      # if interpolating and there are any bad loci with zero called genotypes, remove them
      if(interpolate != FALSE){
        bad.loci <- ifelse(is.na(rdata), 0, 1)
        bad.loci <- which(rowSums(bad.loci) == 0)

        if(length(bad.loci) > 0){
          rdata <- rdata[-bad.loci,]
          meta <- meta[-bad.loci,]
          warning("Some loci had no called genotypes and were removed: ", paste0(bad.loci, collapse = ", "), "\n")
        }
      }

      # interpolate?
      if(interpolate == "bernoulli"){
        rdata <- .interpolate_sn(rdata, "bernoulli")
      }
      else if(interpolate == "af"){
        rdata <- .interpolate_sn(rdata, "af")
      }
      else if(interpolate == "iPCA"){
        rdata <- .interpolate_sn(rdata, "iPCA", ncp = ncp, ncp.max = ncp.max)
      }

      # bind and save
      rdata <- cbind(meta[,-which(colnames(meta) == ".snp.id")], as.data.frame(rdata))
    }

    # sequoia
    else if(output == "sequoia"){
      rdata[is.na(rdata)] <- -9
      rdata <- t(rdata)
      rdata <- matrix(as.numeric(rdata), ncol = ncol(rdata))
      colnames(rdata) <- 1:ncol(rdata)
      
      fix.sex.sequoia <- function(y){
        #format sex for sequoia (1 = F,2 = M,3 = U,4 = H,NA)
        sexes <- c("F", "M", "U", "H", "1", "2", "3", "4")
        sex <- x@sample.meta$Sex
        sex[which(!sex %in% sexes)] <- 3
        sex[which(sex == "F")] <- 1
        sex[which(sex == "M")] <- 2
        sex[which(sex == "U")] <- 3
        sex[which(sex == "H")] <- 4
        return(sex)
      }
      
      #for lifehistory data input needs id, sex, year born
      if(all(c("Sex", "BirthYear", "ID") %in% colnames(x@sample.meta))){
        if(all(c("BY.min", "BY.max") %in% colnames(sample.meta(x)))){
          warning("BirthYear, BY.max, and BY.min columns all found in sample metadata. Defaulting to use BirthYear.\n")
        }
        ID <- x@sample.meta$ID
        sex <- fix.sex.sequoia(x)
        
        #ACTUALLY MAKE THE TABLE
        lhtable <- data.frame(ID=ID,
                              Sex = as.numeric(sex),
                              BirthYear = x@sample.meta$BirthYear,
                              stringsAsFactors = F)
        
      }
      else if(all(c("BY.min", "BY.max", "ID", "Sex") %in% colnames(sample.meta(x)))){
        ID <- x@sample.meta$ID
        sex <- fix.sex.sequoia(x)
        
        #ACTUALLY MAKE THE TABLE
        lhtable <- data.frame(ID=ID,
                              Sex = as.numeric(sex),
                              BY.min = sample.meta(x)$BY.min,
                              BY.max = sample.meta(x)$BY.max,
                              stringsAsFactors = F)
      }
      else{stop("Needs 'ID' and 'Sex' and either 'BirthYear' or 'BY.min' and 'BY.max' columns in sample metadata. \n")}
      
      rownames(rdata) <- ID
      rdata = list(dat=rdata, lh=lhtable)
    }

    # lea
    else{
      rdata <- 2 - rdata
      rdata[is.na(rdata)] <- 9
    }
  }

  # adegenet
  if(output == "adegenet"){
    pop.col <- which(colnames(x@sample.meta) == "pop")

    if(ncol(x@sample.meta) > 1){
      strata <- x@sample.meta[,-ncol(x@sample.meta)]
      if(!is.data.frame(strata)){
        strata <- as.data.frame(strata, stringsAsFactors = F)
        colnames(strata) <- colnames(x@sample.meta)[-ncol(x@sample.meta)]
        row.names(strata) <- colnames(x)
      }
      else{
        strata <- NULL
      }
    }

    if(length(pop.col) > 0){
      rdata <- adegenet::df2genind(t(as.data.frame(x, stringsAsFactors = F)), ncode = 1,
                                   NA.char = substr(x@mDat, 1, nchar(x@mDat)/2),
                                   strata = strata,
                                   loc.names = x@snp.meta$.snp.id,
                                   pop = x@sample.meta$pop)
    }
    else{
      rdata <- adegenet::df2genind(t(as.data.frame(x, stringsAsFactors = F)), ncode = 1,
                                   NA.char = substr(x@mDat, 1, nchar(x@mDat)/2),
                                   strata = strata,
                                   loc.names = x@snp.meta$.snp.id)
    }
  }

  # fasta
  if(output == "fasta"){

    # first, need to randomly draw either a major or minor allele for the heterozygotes
    sn <- format_snps(x, "sn", interpolate = FALSE)
    sn <- sn[,-(1:(ncol(x@snp.meta) - 1))]
    sn <- as.matrix(sn)
    hets <- which(sn == 1)
    nhets <- length(hets)
    draws <- stats::rbinom(nhets, 1, .5)
    draws[draws == 1] <- 2
    sn[hets] <- draws

    # assign major or minor back
    maf.check <- .check_calced_stats(x, ".base", stats = "maf")
    if(!unlist(maf.check)){
      x <- calc_maf(x)
    }
    s <- .get.snpR.stats(x)

    # majors
    sn[is.na(sn)] <- "N"
    majs <- which(sn == 0)
    rmajs <- majs %% nrow(sn)
    rmajs[rmajs == 0] <- nrow(sn)
    rmajs <- s$major[rmajs]
    sn[majs] <- rmajs

    # minors
    mins <- which(sn == 2)
    rmins <- mins %% nrow(sn)
    rmins[rmins == 0] <- nrow(sn)
    rmins <- s$minor[rmins]
    sn[mins] <- rmins

    # paste together by columns
    rdata <- do.call(paste0, as.data.frame(t(sn), stringsAsFactors = F)) # faster than the tidyr version!
    names(rdata) <- colnames(sn)
  }
  
  # VCF
  if(output == "vcf"){
    # meta
    vcf_meta <- c("##fileformat=VCFv4.0",
      paste0("##fileDate=", gsub("-", "", Sys.Date())),
     "##source=snpR",
     '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
     '##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele Count">',
     '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    )
    
    # meta columns
    data_meta <- data.frame(CHROM = snp.meta(x)[,chr],
                            POS = snp.meta(x)[,position],
                            REF = .get.snpR.stats(x)$major,
                            ALT = .get.snpR.stats(x)$minor,
                            QUAL = ".",
                            FILTER = "PASS",
                            INFO = paste0("NS=", .get.snpR.stats(x)$maj.count + .get.snpR.stats(x)$min.count,
                                          ";AC=", .get.snpR.stats(x)$min.count),
                            FORMAT = "GT")
    colnames(data_meta)[1] <- "#CHROM"
    malt <- which(data_meta$ALT == "N")
    mref <- which(data_meta$REF == "N")
    if(length(malt) > 0){
      data_meta$ALT[malt] <- "."
    }
    if(length(mref) > 0){
      data_meta$REF[mref] <- "."
    }
    
    # genotypes, use sn format intermediate
    rdata <- format_snps(x, "sn", interpolate = F)[,-c(1:(ncol(snp.meta(x)) - 1))]
    rdata <- as.matrix(rdata)
    rdata[rdata == 0] <- "0/0"
    rdata[rdata == 1] <- "0/1"
    rdata[rdata == 2] <- "1/1"
    rdata[is.na(rdata)] <- "./."
    rdata <- as.data.frame(rdata)
    rdata <- cbind(data_meta, rdata)
    rdata <- list(meta = vcf_meta, genotypes = rdata)

  }

  #======================return the final product, printing an outfile if requested.=============
  if(outfile != FALSE){
    cat("Writing output file...\n")

    if(any(facets == ".base")){
      facets <- facets[-which(facets == ".base")]
    }

    if(output == "genepop"){ #  if(output %in% c("genepop", "baps"))

      cat("\tPreparing genepop file...\n")
      # get list of snps
      llist <- paste0("SNP", "_", 1:ncol(rdata), ",")
      llist[length(llist)] <- paste0("SNP_", ncol(rdata))

      # write output
      cat(paste0(unlist(.split.facet(outfile))[1], "_genepop\n"), file = outfile)
      cat(llist, "\nPOP\n", file = outfile, append = T) 

      # write the tables, splitting by pop if requested:
      if(length(facets) > 0){
        cat("\tWriting genepop file seperated by populations. First provided facet is treated as pop.\t")

        # sort by pop
        write.facets <- sort(unlist(.split.facet(facets)))
        facet.cols <- match(write.facets, colnames(x@sample.meta))
        pop <- .paste.by.facet(sample.meta(x), facet.cols)
        pop <- gsub(" ", ".", pop)
        rdata$pop <- pop

        pop <- sort(unique(rdata$pop)) # for later


        # sort and remove pop column
        rdata$rnames <- rownames(rdata)
        rdata <- dplyr::arrange(rdata, pop)
        pop.rows <- rdata$pop
        rownames(rdata) <- rdata$rnames
        rdata <- rdata[,-which(colnames(rdata) %in% c("pop", "rnames"))]

        #second loop prints results.
        for (i in 1:(length(pop))){
          cat(pop[i], "\t")
          data.table::fwrite(rdata[pop.rows == pop[i],], outfile, quote = F, sep = "\t", col.names = F, row.names = T, append = T)
          if(i != length(pop)){
            cat("POP\n", file = outfile, append = T)
          }
        }
        cat("\t Done.\n")
      }
      else{
        data.table::fwrite(rdata, outfile, quote = F, sep = "\t", col.names = F, row.names = T, append = T)
      }
    }
    else if(output == "ac"){
      #write the raw output
      data.table::fwrite(rdata, outfile, quote = FALSE, col.names = T, sep = "\t", row.names = F)
      #write a bayescan object if pop list was provided.
      if(length(facets) > 0){
        outfile <- paste0(outfile, ".bayes")
        cat("\tWriting bayescan file seperated by populations. First provided facet is treated as pop.\t")
        pop <- x@sample.meta[,colnames(x@sample.meta) == facets[1]]
        u.pops <- unique(pop)

        trdat <- rdata[rdata$facet == facets[1], c("subfacet", ".snp.id", "n_total", "n_alleles", "ni1", "ni2")]

        #write the header
        cat("[loci]=", nrow(x), "\n\n[populations]=", length(u.pops), "\n\n", file = outfile, sep = "")

        #write the data for each population.
        for(i in 1:length(u.pops)){
          cat("[pop]=", i, "\n", file = outfile, append = T, sep = "") #write header

          tdat <- trdat[trdat$subfacet == u.pops[i],]

          wdat <- cbind(snp = tdat$.snp.id, tdat[,3:6])

          data.table::fwrite(wdat,
                             outfile, col.names = F, row.names = F, quote = F, sep = "\t",
                             append = T) # write the data for this population.

          cat("\n", file = outfile, append = T) # add a line break
        }
      }
    }
    else if(output == "structure" | output == "faststructure"){
      data.table::fwrite(rdata, outfile, quote = FALSE, col.names = F, sep = "\t", row.names = F)
    }
    else if(output == "plink"){
      data.table::fwrite(map, paste0(outfile, ".map"), quote = F, col.names = F, sep = "\t", row.names = F)
      data.table::fwrite(fam, paste0(outfile, ".fam"), quote = F, col.names = F, sep = "\t", row.names = F)
      data.table::fwrite(ped, paste0(outfile, ".ped"), quote = F, col.names = F, sep = "\t", row.names = F)
      data.table::fwrite(bim, paste0(outfile, ".bim"), quote = F, col.names = F, sep = "\t", row.names = F)
    }
    else if(output == "adegenet"){
      saveRDS(rdata, paste0(outfile, ".RDS"))
    }
    else if(output == "fasta"){
      writeobj <- c(rbind(paste0(">", names(rdata)), rdata))
      writeLines(writeobj, outfile)
    }
    else if(output == "lea"){
      utils::write.table(rdata, outfile, quote = FALSE, col.names = F, sep = "", row.names = F)
    }
    else if(output == "colony"){
      data.table::fwrite(rdata, outfile, quote = FALSE, col.names = F, sep = " ", row.names = F)
    }
    else if(output == "sequoia"){
      data.table::fwrite(rdata$dat, paste0("genos_", outfile), quote = FALSE, col.names = F, sep = "\t", row.names = F)
      data.table::fwrite(rdata$lh, paste0("lh_", outfile), quote = FALSE, col.names = T, sep = "\t", row.names = F)
    }
    else if(output == "vcf"){
      writeLines(rdata$meta, outfile)
      data.table::fwrite(rdata$genotypes, outfile, append = T, quote = F, sep = "\t", row.names = F, col.names = T)
    }
    else{
      data.table::fwrite(rdata, outfile, quote = FALSE, col.names = T, sep = "\t", row.names = F)
    }
  }

  #return results
  else{
    return(rdata)
  }
}


#' Checks for duplicated samples in snpRdata.
#'
#' Searches through a snpR dataset and, for every designated sample, determines
#' the proportion of identical genotypes in every other sample. This function
#' \emph{is not overwrite safe}.
#'
#' If an id column is specified, y should contain sample IDs matching those
#' contained in that column. If not, y should contain sample indices instead.
#' The proportion of identical genotypes between matching samples and all other
#' samples are calculated. By default, every sample will be checked.
#'
#' @param x snpRdata object
#' @param y numeric or character, default 1:ncol(x). Designates the sample
#'   indices or IDs in x for which duplicates will be checked.
#' @param id.col character, default NULL. Designates a column in the sample
#'   metadata which contains sample IDs. If provided, y is assumed to contain
#'   sample IDs uniquely matching those in the the sample ID column.
#'
#' @return A list containing: \itemize{ \item{best_matches: } Data.frame listing
#'   the best match for each sample noted in y and the percentage of genotypes
#'   identical between the two samples. \item{data: } A list containing the
#'   match proportion between each sample y and every sample in x, named for the
#'   samples y. }
#'
#' @author William Hemstrom
#' @export
#' @examples
#' # check for duplicates with sample 1
#' check_duplicates(stickSNPs, 1)
#'
#' # check duplicates using the .samp.id column as sample IDs
#' check_duplicates(stickSNPs, 1, id.col = ".sample.id")
check_duplicates <- function(x, y = 1:ncol(x), id.col = NULL){
  #============sanity checks============
  msg <- character()
  if(!is.null(id.col)){
    if(length(id.col) != 1){
      msg <- "Only one ID column may be provided."
    }
    if(!id.col %in% colnames(x@sample.meta)){
      msg <- c(msg,
               "ID column not found in sample metadata.")
    }
    else{
      if(identical(y, 1:ncol(x))){
        y <- x@sample.meta[,id.col]
      }
      if(length(unique(x@sample.meta[,id.col])) != length(x@sample.meta[,id.col])){
        msg <- c(msg,
                 "Each entry in the ID column must be unique.")
      }
      bad.y <- which(!y %in% x@sample.meta[,id.col])
      if(length(bad.y) > 0){
        msg <- c(msg,
                 paste0("Some samples in y not found in ID column: ", paste0(y[bad.y], collapse = ", "), "."))
      }
    }
  }
  else{
    if(!is.numeric(y)){
      msg <- c(msg,
               "y must be numeric if no sample ID column provided.")
    }
    else{
      if(any(y > ncol(x))){
        msg <- c(msg,
                 "All y values must be less than or equal to the number of samples in x.")
      }
    }
  }

  if(length(msg) > 0){
    stop(paste0(msg, collapse = "\n"))
  }


  #============run the duplicate check==========
  # initialize
  out <- vector("list", length(y))
  names(out) <- y
  out.best <- data.frame(sample = y, best_match = character(length(y)),
                         percentage = numeric(length(y)),
                         comparisons = numeric(length(y)), stringsAsFactors = F)

  # do each comparison
  for(i in 1:length(y)){

    # pick out this value
    if(!is.null(id.col)){
      t.samp.id <- which(x@sample.meta[,id.col] == y[i])
    }
    else{
      t.samp.id <- y[i]
    }
    if(length(t.samp.id) == 0){
      out.best$matches[i] <- "bad.ID"
      next()
    }


    #figure out which values are "NN"
    t.samp <- genotypes(x)[,t.samp.id]
    miss <- t.samp == "NN"

    # finish initializing
    out[[i]]$hits <- numeric(ncol(x))
    out[[i]]$comparisons <- numeric(ncol(x))
    if(!is.null(id.col)){
      names(out[[i]]$hits) <- x@sample.meta[,id.col]
    }
    else{
      names(out[[i]]$hits) <- 1:ncol(x)
    }
    names(out[[i]]$comparisons) <- names(out[[i]]$hits)

    # compare to every other sample. Because we don't count "NN" comparisons, we need to explicitly loop.
    for(j in (1:ncol(x))){

      # skip if a self comparison
      if(j == t.samp.id){
        out[[i]]$hits[j] <- 0
        out[[i]]$comparisons[j] <- 0
        next()
      }

      # compare only loci non "NN" in both samples
      c.samp <- genotypes(x)[,j]
      c.miss <- c.samp == "NN"
      u.miss <- which(c.miss | miss)

      # check the proportion of identical genotypes
      valid.comps <- length(c.samp[-u.miss])
      out[[i]]$hits[j] <- sum(t.samp[-u.miss] == c.samp[-u.miss])/valid.comps
      out[[i]]$comparisons[j] <- valid.comps
    }

    # figure out and save data on the "best", or most identical, hit.
    best <- which.max(out[[i]]$hits)
    out[[i]]$best <- out[[i]]$hits[best]
    names(out[[i]]$best) <- names(out[[i]]$hits)[best]

    # save to output summary.
    out.best$best_match[i] <- paste0(names(out[[i]]$best), collapse = ", ")
    out.best$percentage[i] <- out[[i]]$best
    out.best$comparisons[i] <- out[[i]]$comparisons[best]
  }

  # return
  return(list(best_matches = out.best, data = out))

}



#' Fetch the allele frequencies for all SNPs for each level of each requested
#' facet.
#'
#' Fetch allele frequencies for all SNPs for each level of all the requested
#' facets. Major and minor allele frequencies will be interleved, with the major
#' allele first for each locus. Note that this particular function is not
#' overwrite-safe.
#'
#' @param x snpRdata object
#' @param facets character, default NULL. The facets for which to calculate
#'   allele frequencies. See \code{\link{Facets_in_snpR}} for details.
#'
#' @return A named, nested list containing allele frequency matrices for each
#'   facet level for all requested facets.
#'
#' @author William Hemstrom
#' @export
tabulate_allele_frequency_matrix <- function(x, facets = NULL){
  
  ..ord <- NULL
  
  #==================prep and sanity check==================
  msg <- character()
  
  # check facets and get maf
  facets <- .check.snpR.facet.request(x, facets, "none", T)
  facet_types <- facets[[2]]
  facets <- facets[[1]]
  needed.sample.facets <- .check.snpR.facet.request(x, facets)
  if(any(facet_types == "snp")){
    needed.sample.facets <- c(needed.sample.facets, ".base")
    needed.sample.facets <- unique(needed.sample.facets)
  }
  ## add any missing maf data
  missing_mafs <- .check_calced_stats(x, needed.sample.facets, "maf")
  if(any(!unlist(missing_mafs))){
    x <- calc_maf(x, names(missing_mafs)[which(!unlist(missing_mafs))])
  }
  am <- .get.snpR.stats(x, needed.sample.facets)
  
  # grab column names
  maj_min <- .get.snpR.stats(x)
  maj_min <- paste0(rep(unique(maj_min$.snp.id), each = 2), "_",
                    c(maj_min$major, maj_min$minor)[rep(1:nrow(maj_min), each = 2) + (0:1) * nrow(maj_min)])
  
  #==================run============
  # intialize output
  ## need to do calcs for each sample level facet only
  sample_facet_freqs <- vector("list", length(needed.sample.facets))
  names(sample_facet_freqs) <- needed.sample.facets
  
  ## output
  out <- vector("list", length = length(facets))
  names(out) <- facets
  
  # run for each facet
  for(i in 1:length(sample_facet_freqs)){

    # melt the allele frequencies down into a matrix:
    tam <- data.table::dcast(as.data.table(am[which(am$facet == needed.sample.facets[i]),]), subfacet ~ .snp.id, value.var = "maf")
    pops <- tam[,1]
    tam <- tam[,-1]
    amb <- 1 - tam
    
    # interleve major and minor frequencies and save
    ord <- rep(1:ncol(tam), each = 2) + (0:1) * ncol(tam)
    amc <- .fix..call(cbind(amb, tam)[,..ord])
    colnames(amc) <- maj_min
    amc <- as.data.frame(amc)
    rownames(amc) <- unlist(pops)
    
    sample_facet_freqs[[i]] <- as.matrix(amc)
  }
  #==================break up for facets=========
  for(i in 1:length(facets)){
    # split up by snp levels if requested
    if(facet_types[i] == "complex" | facet_types[i] == "snp"){
      if(facet_types[i] == "snp"){
        t.samp.facet <- ".base"
      }
      else{
        t.samp.facet <- .check.snpR.facet.request(x, facets[i], "snp")
      }
      t.snp.facet <- .check.snpR.facet.request(x, facets[i], "sample")
      snp.parts <- .get.task.list(x, t.snp.facet)[,-c(1:2)]
      
      tout <- vector("list", length = nrow(snp.parts))
      names(tout) <- snp.parts[,2]
      for(j in 1:nrow(snp.parts)){
        # these are the snpIDs that are a part of this facet level
        include_snps <- unique(am$.snp.id[which(am[,snp.parts[j,1]] == snp.parts[j,2])])
        
        # grab just those snps
        tout[[j]] <- sample_facet_freqs[[which(needed.sample.facets == t.samp.facet)[1]]] # add everything
        tout[[j]] <- tout[[j]][,which(as.numeric(gsub("_.+", "", colnames(tout[[j]]))) %in% include_snps), drop = F] # subset the requested snps only
      }
      # assign back, nesting with the snp facet name
      out[[i]] <- tout
    }
    else{
      if(is.null(.check.snpR.facet.request(x, facets[i]))){
        facets[i] <- ".base"
      }
      out[[i]] <- list(.base = sample_facet_freqs[[which(needed.sample.facets == facets[i])]])
    }
  }
  
  #=================return============
  x <- .update_calced_stats(x, facets, "allele_frequency_matrix")
  return(.merge.snpR.stats(x, stats = out, type = "allele_frequency_matrices"))
}



#' Filter down to one SNPs every \emph{n} bases.
#'
#' Selects one SNPs every \emph{n} bases. This can be
#' used to filter for SNPs in tight LD, such as when only one SNP per RAD-tag is
#' desired.
#' 
#' Note that this approach takes the first SNP every \emph{n} bases, and so is
#' non-random. \code{\link{filter_snps}} can be used beforehand to ensure that
#' the selected SNPs are above a desired quality threshold to ensure that poor
#' SNPs are not selected over loci with more robust genotyping data.
#' 
#' It is strongly recommended to provide a SNP facet to the facet argument. This
#' will define "chromosomes" for the purpose of selecting SNPs. If not set, snpR
#' will assume that \emph{all SNPs are on the same chromosome}, which may
#' produce undesired results.
#'
#' @param x snpRdata object
#' @param facet character, default NULL. SNP facet specifying chromosomes or 
#'   scaffolds. SNP positions will be independantly considered depending on the
#'   facet level.
#' @param n Integer. Specifies the minium distance between selected SNPs.
#' 
#' @return A snpRdata object containing the selected SNPs.
#' 
#' @author William Hemstrom
#' @export
#' 
#' @examples 
#' # put large gaps inbetween the example data
#' gapped <- gap_snps(stickSNPs, "chr", 100000)
#'
#' # filter first to grab only very good SNPs
#' vgood <- filter_snps(stickSNPs, min_ind = .8)
#' gapped <- gap_snps(vgood, "chr", 100000)
gap_snps <- function(x, facet = NULL, n){
  #==========sanity checks=========
  if(!is.snpRdata(x)){
    stop("x must be a snpRdata object.\n")
  }
  msg <- character()
  
  
  ffacet <- .check.snpR.facet.request(x, facet, remove.type = "none", return.type = TRUE)
  
  if(any(ffacet[[2]] %in% c("sample", "complex"))){
    msg <- c(msg, "Non SNP facets are not permitted.\n")
  }
  
  facet <- .check.snpR.facet.request(x, facet, "sample")
  
  if(length(facet) != 1){
    msg <- c(msg, "Only one facet may be provided to gap_snps. Note that complex
             SNP facets may be used (e.g. chr.scaffold).")
  }
  
  if(!"position" %in% colnames(snp.meta(x))){
    msg <- c(msg, "There must be a 'position' column in the SNP metadata to use gap_snps!")
  }
  
  if(facet[1] == ".base"){
    warning("Gapping SNPs by the base facet is not recommended. Consider using a SNP facet instead.\n")
  }
  
  if(".tfacet.ext" %in% colnames(snp.meta(x))){
    msg <- c(msg, "A column named '.tfacet.ext' cannot be in the SNP metadata.\n")
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  #==============subfunction==============
  retained <- function(pos, n){
    fsnp <- which.min(pos)
    cs <- pos[fsnp]
    # use a while loop to grab the snps
    # should be a fairly quick loop, since it loops through each gap, not each 
    # snp. Might still be a bit slow with really large data, can fix later if
    # I get requests...
    while(cs[length(cs)] + n <= max(pos)){
      cs <- c(cs,
              min(pos[pos >= cs[length(cs)] + n]))
    }
    ret <- match(cs, pos)
    return(ret)
  }
  
  #==============find the SNPs to retain=====================
  # get a task list, then run for each task
  task.list <- .get.task.list(x, facet)
  keep.snps <- numeric()
  
  # no par option for now, should usually be quick. Can add if requested.
  for(i in 1:nrow(task.list)){
    tsm <- snp.meta(x)[.fetch.snp.meta.matching.task.list(x, task.list[i,]),]
    ret <- retained(tsm$position, n)
    keep.snps <- c(keep.snps, tsm$.snp.id[ret])
  }
  
  #==============subset and return========
  
  # grab the snp indexes to keep
  ids <- sort(match(keep.snps, snp.meta(x)$.snp.id))
  cat("Selected ", length(ids), "out of ", nrow(x), "SNPs.\n")
  
  # subset
  return(x[ids,])
}

#' Gather citations for methods used with a snpRdata object
#' 
#' snpR will automatically track the methods used when for calculations on a
#' snpRdata object. Using \code{\link{citations}} on that object will provide details 
#' on the the methods used, and can optionally write a .bib bibtex formatted
#' library containing citations for these methods.
#' 
#' Printed outputs contain the statistic calculated, the in-line citation for
#' the method used, a bibtex key which corresponds to the .bib library written 
#' if the outbib argument is used, and a quick note giving any details.
#' 
#' @param x snpRdata object
#' @param outbib character, default FALSE. An optional file path to which a .bib
#'   bibtex library containing all of the citations for the used methods will be
#'   written.
#' @param return_bib logical, default FALSE. If TRUE, returns a list containing
#'   the bib entries and their details. The bib entries are formatted according
#'   to \code{\link[RefManageR]{BibEntry}}.
#'   
#' @author William Hemstrom
#' @returns  If return_bib is TRUE, a list containing four parts: 
#' \itemize{\item{keys: } A vector of bibtex keys for each method.
#' \item{stats: } A vector of the stats used.
#' \item{details: } A vector of details for each method.
#' \item{bib: } A \code{\link[RefManageR]{BibEntry}} for each citation.}
#' 
#' @export
#' @examples 
#' # calculate pi
#' x <- calc_pi(stickSNPs)
#' 
#' # fetch citations
#' citations(x)
#' 
#' # fetch citations and write bibtex .bib to citations.bib
#' ## not run
#' \dontrun{
#' citations(x, outbib = "citations.bib")
#' }
citations <- function(x, outbib = FALSE, return_bib = FALSE){
  #==========sanity checks=======
  if(!is.snpRdata(x)){
    stop("x must be a snpRdata object.\n")
  }
  
  .check.installed("bibtex")
  .check.installed("RefManageR")
  
  #==========grab bib============
  bib.file <- system.file("extdata", "snpR_citations.bib", package = "snpR")
  bib <- RefManageR::ReadBib(bib.file)
  
  #==========filter bib==========
  keys <- as.character(unlist(purrr::map(x@citations, "key")))
  bib <- bib[keys]
  
  #==========shout at the user=====
  deets <- unlist(purrr::map(x@citations, "details"))
  stats <- names(x@citations)
  
  cat("Citations for methods used thus far: \n")
  for(i in 1:length(keys)){
    cat("==============================================\n")
    cat("Statistic: ", stats[i], "\n")
    cat("Citation: ", RefManageR::Cite(bib[keys[i]]), "\n")
    cat("Bibtex key: ", keys[i], "\n")
    cat("Details: ", deets[i], "\n")
  }
  cat("==============================================\n\n")
  
  #==========print bib=============
  if(!isFALSE(outbib)){
    
    if(file.exists(outbib)){
      current_bib <- RefManageR::ReadBib(outbib)
      not_in_current <- which(!keys %in% names(current_bib))
      
      
      if(length(not_in_current) > 0){
        #==========filter bib==========
        bib <- bib[keys[not_in_current]]
        
        current_bib <- c(current_bib, bib)
        
        RefManageR::WriteBib(current_bib, outbib, verbose = FALSE)
      }
    }
    
    
    else{
      RefManageR::WriteBib(bib, outbib, verbose = FALSE)
      cat(".bib file can be found at: ", outbib, "\n")
    }
    
  }
  
  if(!isFALSE(return_bib)){
    return(list(keys = keys, stats = stats, details = deets, bib = bib))
  }
}
