
#'Subset snpRdata objects
#'
#'Subsets snpRdata objects by specific snps, samples, facets, subfacets, etc.
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
#'@param i numeric, default \code{1:nrow(x)}. Row numbers corresponding to SNPs
#'  to keep. Negative subscripts allowed.
#'@param j numeric, default \code{1:nrow(x)}. Row numbers corresponding to SNPs
#'  to keep. Negative subscripts allowed.
#'@param ... Facet subsetting may be specified by providing the facet as an
#'  argument and then providing the levels to keep or remove. Setting pop =
#'  "ASP", for example, would keep only samples in the "ASP" level of the pop
#'  facet, and setting pop.fam = "ASP.A" would keep only samples the ASP pop and
#'  the A fam. Negative subscripts are allowed: pop = -c("ASP", "PAL") would
#'  remove samples in the ASP or PAL pops. Subsetting by multiple facets is
#'  supported, although negative and positive subscripts cannot be mixed across
#'  sample or SNP facets. They may be mixed between the two.
#'@param .facets Character, default NULL. Sample facet to select by. Alternative
#'  to direct specification as in \code{pop = "ASP"}.
#'@param .subfacets Character, default NULL. Sample facet level matching a level
#'  of the provided \code{.facets}. Alternative to direct specification as in
#'  \code{pop = "ASP"}.
#'@param .snp.facets Character, default NULL. SNP facet to select by.
#'  Alternative to direct specification as in \code{chr = "groupIX"}.
#'@param .snp.subfacets Character, default NULL. SNP facet level matching a
#'  level of the provided \code{.snp.facets}. Alternative to direct
#'  specification as in \code{chr = "groupIX"}.
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
#' 
#' # using the .facet etc arguments, useful when the facet is stored as an
#' # object
#' target_facet <- "pop"
#' target_subfacet <- "ASP"
#' subset_snpR_data(stickSNPs, 
#'                  .facets = target_facet, 
#'                  .subfacets = target_subfacet)
NULL


#' @export
#' @describeIn subset_snpRdata subset_snpR_data
subset_snpR_data <- function(x, .snps = 1:nsnps(x), .samps = 1:nsamps(x), ..., .facets = NULL, .subfacets = NULL, .snp.facets = NULL, .snp.subfacets = NULL){
  #============extract facet info===============
  .argnames <- match.call()
  .argnames <- as.list(.argnames)
  .argnames <- .argnames[-1]

  .is.facet <- which(!names(.argnames) %in% c("x", ".snps", ".samps", ".facets", ".subfacets", ".snp.facets", ".snp.subfacets"))
  if(length(.is.facet) > 0){
    
    if(any(!is.null(.facets), !is.null(.subfacets), 
           !is.null(.snp.facets), !is.null(.snp.subfacets))){
      stop("The '.facets', '.subfacets', '.snp.facets', or '.snp.subfacets' arguments cannot be used with arguments other than '.snps' or '.samps' (such as pop = 'ASP'). If you wish to subset using the 'facet = subfacet' approach, please do not use these arguments.\n")
    }
    
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
      
      # try to eval text in two ways: once for a basic string, once for code the evaluates to a string
      .res1 <- try(list(eval(parse(text = .argnames[.is.facet][.candidate_facets]))), silent = TRUE)
      
      if(methods::is(.res1, "try-error")){
        .argnames[.is.facet][.candidate_facets] <- list(eval(parse(text = paste0("c(\"", .argnames[.is.facet][.candidate_facets], "\")"))))
      }
      else{
        .argnames[.is.facet][.candidate_facets] <- .res1
      }
    }
    
    
  }
  else if(any(!is.null(.facets), !is.null(.subfacets), 
              !is.null(.snp.facets), !is.null(.snp.subfacets))){
    return(.subset_snpR_data(x, 
                             snps = .snps, 
                             samps = .samps, 
                             facets = .facets, 
                             subfacets = .subfacets, 
                             snp.facets = .snp.facets, 
                             snp.subfacets = .snp.subfacets))
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
      stop("Facets and subfacets are now desginated directly using, for example, pop = c('my.pop1', 'my.pop2'), not using the 'facet', 'subfacet', etc arguments. If you wish to subset using the facet = 'pop', subfacet = 'ASP'-style approach, please use the '.facets', '.subfacets', '.snp.facets', or '.snp.subfacets' arguments.\n")
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
    if(max(.snps) > nrow(x) | min(.snps) == 0){
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
    if(max(.samps) > ncol(x) | min(.samps) == 0){
      msg <- c(msg, "All requested samples must be within 1:nsnps(x).\n")
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
  
  # check that we still have data
  if(length(.snps) == 0){
    stop("No SNPs remain after subsetting!\n")
  }
  if(length(.samps) == 0){
    stop("No samples remain after subsetting!\n")
  }
  
  
  nsampm <- sample.meta(x)[.samps,]
  
  nsnpm <- snp.meta(x)[.snps,]
  
  r <- try(x@filters, silent = TRUE)
  if(methods::is(r, "try-error")){
    return(import.snpR.data(genotypes(x)[.snps, .samps, drop = FALSE], snp.meta = nsnpm,
                            sample.meta = nsampm, mDat = x@mDat,
                            .skip_filters = TRUE))
  }
  else{
    return(import.snpR.data(genotypes(x)[.snps, .samps, drop = FALSE], snp.meta = nsnpm,
                            sample.meta = nsampm, mDat = x@mDat, .pass_filters = x@filters, 
                            .skip_filters = TRUE))
  }
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
  
  # check that we still have data
  if(length(snps) == 0){
    stop("No SNPs remain after subsetting!\n")
  }
  if(length(samps) == 0){
    stop("No samples remain after subsetting!\n")
  }
  
  # subset
  if(!identical(samps, 1:ncol(x))){
    dat <- genotypes(x)[snps, samps]
    if(length(samps) == 1){
      dat <- as.data.frame(dat, stringsAsFactors = F)
    }
    dat <- import.snpR.data(dat, snp.meta = snp.meta(x)[snps,,drop = FALSE], sample.meta = sample.meta(x)[samps,,drop = FALSE], mDat = x@mDat)
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
    if(!is.null(x@sn$sn)){
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
    
    x <- snpRdata(.Data = genotypes(x)[which(snp.meta(x)$.snp.id %in% snps),, drop = FALSE],
                  sample.meta = sample.meta(x),
                  snp.meta = snp.meta(x)[which(snp.meta(x)$.snp.id %in% snps),, drop = FALSE],
                  facet.meta = x@facet.meta[x@facet.meta$.snp.id %in% snps,, drop = FALSE],
                  mDat = x@mDat,
                  snp.form = x@snp.form,
                  ploidy = .ploidy(x),
                  bi_allelic = .is.bi_allelic(x),
                  filters = x@filters,
                  data.type = ifelse("data.type" %in% methods::slotNames(x), x@data.type, "genotypic"),
                  geno.tables = list(gs = x@geno.tables$gs[x@facet.meta$.snp.id %in% snps,, drop = FALSE],
                                     as = x@geno.tables$as[x@facet.meta$.snp.id %in% snps,, drop = FALSE],
                                     wm = x@geno.tables$wm[x@facet.meta$.snp.id %in% snps,, drop = FALSE]),
                  # ac = x@ac[x@facet.meta$.snp.id %in% snps,],
                  facets = x@facets,
                  facet.type = x@facet.type,
                  stats = x@stats[x@stats$.snp.id %in% snps, drop = FALSE],
                  pairwise.stats = x@pairwise.stats[x@pairwise.stats$.snp.id %in% snps,, drop = FALSE],
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
#'To counter this, the \code{re_run} argument can be used to pass the data
#'through a second filtering step after individuals are removed. By default, the
#'"partial" re-run option is used, which re-runs only the non-polymorphic filter
#'(if it was originally set). The "full" option re-runs all set filters. Note
#'that re-running any of these filters may cause individuals to fail the
#'individual filter after loci removal, and so subsequent tertiary re-running of
#'the individual filters, followed by the loci filters, and so on, could be
#'justified. This is not done automatically here, however, if the
#'\code{inds_first} option is selected to filter poorly sequenced individuals
#'prior to applying loci filters, any re-run option will re-run individuals
#'filters after the loci filters have been applied.
#'
#'Via the \code{maf_facets} argument, this function can filter by minor allele
#'frequencies in either \emph{all} samples or \emph{each level of a supplied
#'sample specific facet and the entire dataset}. In the latter case, any SNPs
#'that pass the maf filter in \emph{any} facet level are considered to pass the
#'filter. The latter should be used in instances where population sizes are very
#'different, there are \emph{many} populations, and/or allele frequencies are
#'very different between populations and thus common alleles of interest in one
#'population might be otherwise filtered out.
#'
#'The \code{hwe_facets} argument is the inverse of this: loci will be removed if they
#'fail the provided hwe filter in any facet level. In both cases, Facets should
#'be provided as described in \code{\link{Facets_in_snpR}}.
#'
#'
#'@param x snpRdata object.
#'@param maf numeric between 0 and 1 or FALSE, default FALSE. Minimum acceptable
#'  minor allele frequency. Either maf or mac filtering are allowed, not both.
#'@param mac integer, between 0 and 1 minus the total number of individuals.
#'  Loci with less \emph{or equal to} this many minor alleles will be removed.
#'  Either maf or mac filtering are allowed, not both. \code{mac = 1} is
#'  therefore a singleton filter.
#'@param mgc integer, between 0 and the total number of individuals/2. Loci
#'  where the minor allele is present in less \emph{or equal to} this many
#'  individuals will be removed. Either mac or mgc filtering are allowed, not
#'  both. \code{mgc = 1} is analagous to a singleton filter, but will also
#'  remove loci with one homozygous minor individual.
#'@param hf_hets numeric between 0 and 1 or FALSE, default FALSE. Maximum
#'  acceptable heterozygote frequency.
#'@param hwe numeric between 0 and 1 or FALSE, default FALSE. Minimum acceptable
#'  HWE p-value.
#'@param hwe_excess_side character, default "both". Options:
#'  \itemize{\item{heterozygote:} only loci with a heterozygote excess are 
#'  removed.
#'  \item{homozygote:} only loci with a homozygote excess are removed.
#'  \itemize{\item{both:} loci with either a heterozygote or homozygote excess 
#'  are removed.}}.
#'@param fwe_method character, default "none". Option to use for Family-Wise
#'  Error rate correction for HWE filtering. If requested, only SNPs with
#'  p-values below the alpha provided to the \code{hwe} argument \emph{after FWE
#'  correction} will be removed. See \code{\link[stats]{p.adjust}} for
#'  information on method options.
#'@param singletons logical, default FALSE. Depricated, use \code{mac = 1} to
#'  remove singletons. If TRUE, removes singletons (loci where there is only a
#'  single minor allele). If population sizes are reasonably high, this is more
#'  or less redundant if \code{maf} is also set.
#'@param min_ind numeric between 0 and 1 or FALSE, default FALSE. Minimum
#'  proportion of individuals in which a loci must be sequenced.
#'@param min_loci numeric between 0 and 1 or FALSE, default FALSE. Minimum
#'proportion of SNPs at which an individual must be genotyped.
#'@param inds_first logical, default FALSE. If TRUE, individuals will be
#'  filtered out for missing data if that option is selected prior to loci being
#'  filtered. Otherwise loci are filtered first.
#'@param remove_garbage numeric between 0 and 1 or FALSE, default FALSE. Optionally
#'  do a filter to remove very poorly sequenced individuals and loci
#'  jointly before applying other filters. This can be used to reduce biases
#'  caused by very bad loci/individuals prior to full filtering. This number
#'  should be lower than either the \code{min_ind} or \code{min_loci} parameters
#'  if those arguments are used and should generally be quite permissive to 
#'  remove only truly bad loci or individuals.
#'@param re_run character or FALSE, default "partial". When individuals are
#'  removed via min_ind, it is possible that some SNPs that initially passed
#'  filtering steps will now violate some filters. SNP filters can be re-run
#'  automatically via several methods: \itemize{ \item{partial: } Re-filters for
#'  non-polymorphic loci (non_poly) only, if that filter was requested
#'  initially. \item{full: } Re-runs the full filtering scheme (save for
#'  min_loci).} Note: if \code{inds_first = TRUE}, all re-run options other than 
#'  FALSE will re-run the individual filter for missing data again after loci 
#'  filtering.
#'@param maf_facets character or NULL, default NULL. Defines a sample facet over
#'  which the minor allele frequency can be checked. SNPs will only fail the maf
#'  filter if they fail in every level of every provided facet.
#'@param hwe_facets character or NULL, default NULL. Defines a sample facet over
#'  which the hwe filter can be checked. SNPs will fail the hwe filter if they
#'  fail in any level of any provided facet.
#'@param non_poly logical, default TRUE. If TRUE, non-polymorphic loci will be
#'  removed.
#'@param bi_al logical, default TRUE. If TRUE, loci with more than two alleles
#'  will be removed. Note that this is mostly an internal argument and should
#'  rarely be used directly, since import.snpR.data and other snpRdata object
#'  creation functions all pass SNPs through this filter because many snpR
#'  functions will fail to work if there are more than two alleles at a locus.
#'@param verbose Logical, default TRUE. If TRUE, some progress updates and
#'  filtering notes will be printed to the console.
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
filter_snps <- function(x, maf = FALSE, 
                        mac = 0,
                        mgc = 0,
                        hf_hets = FALSE, 
                        hwe = FALSE, fwe_method = "none",
                        hwe_excess_side = "both",
                        singletons = FALSE,
                        min_ind = FALSE,
                        min_loci = FALSE,
                        inds_first = FALSE,
                        remove_garbage = FALSE,
                        re_run = "partial", 
                        maf_facets = NULL,
                        hwe_facets = NULL,
                        non_poly = TRUE, 
                        bi_al = TRUE,
                        verbose = TRUE){

  #==============do sanity checks====================
  if(singletons){
    warning("The singletons argument is depriceated. Please use mac = 1 instead!")
    if(mac == 0){
      mac <- 1
    }
  }
  
  if(maf){
    if(!is.numeric(maf)){
      stop("maf must be a numeric value.")
    }
    if(length(maf) != 1){
      stop("maf must be a numeric vector of length 1.")
    }
    if(mac != 0){
      stop("mac and maf cannot both be set.\n")
    }
    if(maf < 0 | maf >= .5){
      stop("maf must be greater than or equal to 0 and less than .5")
    }
    
    if(!is.null(maf_facets) & !isFALSE(maf_facets)){
      if(length(maf_facets) == 1 & maf_facets[1] == ".base"){
        maf_facets <- NULL
      }
      else{
        maf_facets <- .check.snpR.facet.request(x, maf_facets, "none")
        
        # add any needed facets...
        miss.facets <- maf_facets[which(!(maf_facets %in% x@facets))]
        if(length(miss.facets) != 0){
          if(verbose){cat("Adding missing facets...\n")}
          # need to fix any multivariate facets (those with a .)
          x <- .add.facets.snpR.data(x, miss.facets)
        }
        
        # check for bad facets to remove (those that don't just consider samples)
        if(any(!x@facet.type[x@facets %in% maf_facets] %in% c("sample", ".base"))){
          vio.facets <- x@facets[match(maf_facets, x@facets)]
          vio.facets <- vio.facets[which(x@facet.type[x@facets %in% maf_facets] != "sample")]
          warning(paste0("Facets over which to maf.filter must be sample specific facets, not snp specific facets! Removing non-sample facets: \n", paste0(vio.facets, collapse = " "), ".\n"))
          maf_facets <- maf_facets[-which(maf_facets %in% vio.facets)]
          if(length(maf_facets) == 0){maf_facets <- NULL}
        }
      }
    }
  }
  
  if(mac){
    if(mgc){
      stop("mac and mgc cannot both be set.")
    }
    if(!is.numeric(mac)){
      stop("mac must be a numeric value.")
    }
    if(length(mac) != 1){
      stop("mac must be a numeric vector of length 1.")
    }
    if(mac >= nsamps(x) | mac < 0){
      stop("mac must be greater than or equal to zero and less than the number of samples (a maf of 0.5 in fully sequenced loci).")
    }
    if(mac != as.integer(mac)){
      stop("mac must be an integer (but can be a numeric object).")
    }
  }
  
  if(mgc){
    if(!is.numeric(mgc)){
      stop("mgc must be a numeric value.")
    }
    if(length(mgc) != 1){
      stop("mgc must be a numeric vector of length 1.")
    }
    if(mgc >= nsamps(x)/2 | mgc < 0){
      stop("mgc must be greater than or equal to zero and less than half the number of samples (a maf of 0.5 in fully sequenced loci).")
    }
    if(mgc != as.integer(mgc)){
      stop("mgc must be an integer (but can be a numeric object).")
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
    
    good_hwe_sides <- c("both", "homozygote", "heterozygote")
    if(!hwe_excess_side %in% good_hwe_sides){
      stop("hwe_excess_sides must be one of: ", paste0(good_hwe_sides, collapse = ", "))
    }
    
    if(!is.null(hwe_facets) & !isFALSE(hwe_facets)){
      if(length(hwe_facets) == 1 & hwe_facets[1] == ".base"){
        hwe_facets <- NULL
      }
      else{
        hwe_facets <- .check.snpR.facet.request(x, hwe_facets)
      }
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
  
  if(!isFALSE(remove_garbage)){
    if(!is.numeric(remove_garbage)){
      stop("remove_garbage must be a numeric value.")
    }
    if(length(remove_garbage) != 1){
      stop("remove_garbage must be a numeric vector of length 1.")
    }
    if(remove_garbage > 1 | remove_garbage < 0){
      stop("remove_garbage is a minimum proportion of sequenced individuals/loci, and so must be between 0 and 1.\n")
    }
    
    if(is.numeric(min_loci)){
      if(min_loci <= remove_garbage){
        stop("The min_loci threshold should be higher than the garbage removal threshold if supplied.\n")
      }
    }
    if(is.numeric(min_ind)){
      if(min_ind <= remove_garbage){
        stop("The min_ind threshold should be higher than the garbage removal threshold if supplied.\n")
      }
    }
  }
  
  if(re_run != FALSE){
    if(re_run != "partial" & re_run != "full"){
      if(verbose){cat("re_run must be set to partial or full if not FALSE.\n")}
    }
  }
  
  #==============set up, get values used later, clean up data a bit, define subfunctions==========
  if(verbose){cat("Initializing...\n")}
  
  #get headers
  headers <- x@snp.meta
  snp_form <- x@snp.form
  mDat <- x@mDat
  
  
  #function to filter by loci, to be called before and after min ind filtering (if that is requested.)
  filt_by_loci <- function(){
    # Store filter status in vio.snps. Those that are violating a filter will be marked TRUE, remove these.
    
    #==========================run filters: bi_allelic/non_poly========================
    vio.snps <- logical(nrow(x)) #vector to track status
    
    amat <- x@geno.tables$as[x@facet.meta$facet == ".base",,drop = FALSE]
    gmat <- x@geno.tables$gs[x@facet.meta$facet == ".base",,drop = FALSE]
    wmat <- x@geno.tables$wm[x@facet.meta$facet == ".base",,drop = FALSE]
    
    # non-biallelic and non-polymorphic loci
    if(bi_al | non_poly){
      bimat <- ifelse(amat, TRUE, FALSE)
      
      if(bi_al){
        if(verbose){cat("Filtering non-biallelic loci...\n")}
        bi <- ifelse(rowSums(bimat) > 2, T, F) # if false, should keep the allele
        if(verbose){cat(paste0("\t", sum(bi), " bad loci\n"))}
        vio.snps[which(bi)] <- T
        x <- .update_filters(x, "bi-allelic", NA, NA)
      }
      
      if(non_poly){
        if(verbose){cat("Filtering non_polymorphic loci...\n")}
        np <- ifelse(rowSums(bimat) < 2, T, F) # if false, should keep the allele
        if(verbose){cat(paste0("\t", sum(np), " bad loci\n"))}
        vio.snps[which(np)] <- T
        x <- .update_filters(x, "non-polymorphic", NA, NA)
      }
    }
    
    #========min inds=======
    if(min_ind){
      if(verbose){cat("Filtering loci sequenced in few individuals...\n")}
      mi <- wmat[,colnames(wmat) == mDat]
      mi <- (nrow(x@sample.meta) - mi)/nrow(x@sample.meta) < min_ind
      vio.snps[which(mi)] <- T
      if(verbose){cat(paste0("\t", sum(mi), " bad loci\n"))}
      
      x <- .update_filters(x, "poorly sequenced loci--min_ind", min_ind, ".base")
    }
    
    #========minor allele frequency, both total and by pop. Should only run if bi_al = TRUE.=========
    if(maf){
      #if not filtering with multiple pops
      if(is.null(maf_facets) | isFALSE(maf_facets)){
        if(verbose){cat("Filtering low minor allele frequencies, no pops...\n")}
        
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
        if(verbose){cat(paste0("\t", sum(mafs), " bad loci\n"))}
        
        
        vio.snps[which(mafs)] <- T
      }
      else{
        if(verbose){cat("Filtering low minor allele frequencies by facet.\n")}
        # pmafs <- logical(nrow(x))
        
        # see if we need to calculate mafs
        if(any(colnames(x@stats) == "maf")){ # mafs have been calculated
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
        if(verbose){cat(paste0("\t", length(maf.vio), " bad loci\n"))}
        vio.snps[maf.vio] <- T
      }
      
      
      # update note-keeping
      x <- .update_filters(x, "Minor Allele Frequency--maf", maf, 
                           ifelse(is.null(maf_facets) | isFALSE(maf_facets), ".base", maf_facets))
      
    }
    
    #========mac/mgc========================================
    if(mac | mgc){
      if(mgc){
        # in this case, this is the count of *individuals with the minor allele*
        if(verbose){cat("Filtering low minor genotype counts.\n")}
        hs <- which(substr(colnames(gmat), 1, snp_form/2) != substr(colnames(gmat), (snp_form/2) + 1, snp_form))
        het_c <- matrixStats::rowSums2(gmat[,hs])
        min_g_c <- matrixStats::rowSums2(gmat[,-hs]) - matrixStats::rowMaxs(gmat[,-hs])
        
        tgac <- het_c + min_g_c
        
        mgc.vio <- which(tgac <= mgc)
        if(verbose){cat(paste0("\t", length(mgc.vio), " bad loci\n"))}
        
        vio.snps[mgc.vio] <- T
        
        x <- .update_filters(x, "minor genotype count--mgc", mgc, ".base")
      }
      else{
        if(verbose){cat("Filtering low minor genotype counts.\n")}
        
        # singletons exist wherever the total allele count - the maf count is 1.
        mac.vio <- which(matrixStats::rowSums2(amat) - matrixStats::rowMaxs(amat) <= mac)
        
        if(verbose){cat(paste0("\t", length(mac.vio), " bad loci\n"))}
        
        vio.snps[mac.vio] <- T
        
        x <- .update_filters(x, "minor allele count--mac", mac, ".base")
      }
    }
    #========hf_hets. Should only run if bi_al = TRUE.==========
    if(hf_hets){
      if(verbose){cat("Filtering high frequency heterozygote loci...\n")}
      
      # get heterozygote frequency
      hs <- which(substr(colnames(gmat), 1, snp_form/2) != substr(colnames(gmat), (snp_form/2) + 1, snp_form))
      het_f <- rowSums(gmat[,hs])/rowSums(gmat)
      
      # check violation
      het_f <- het_f > hf_hets #if false, heterozygote frequency is lower than cut-off, keep locus
      if(verbose){cat(paste0("\t", sum(het_f), " bad loci\n"))}
      vio.snps[which(het_f)] <- T
      
      x <- .update_filters(x, "high-frequency heterozygotes--hf_hets", hf_hets, ".base")
    }
    
    #========hwe violation======================================
    if(hwe){
      if(verbose){cat("Filtering loci out of hwe...\n")}
      
      # no facets, easy
      if(is.null(hwe_facets) | isFALSE(hwe_facets)){
        
        if(!.check_calced_stats(x, ".base", "hwe")$.base){
          .make_it_quiet(x <- calc_hwe(x))
        }
        if(hwe_excess_side != "both"){
          calced <- .check_calced_stats(x,".base", c("fis"))$.base
          if(!calced){
            .make_it_quiet(x <- calc_fis(x))
          }
        }
        
        
        phwe <- x@stats$pHWE[x@stats$facet == ".base"]
        if(fwe_method != "none"){
          phwe <- stats::p.adjust(phwe, method = fwe_method[1])
        }
        phwe <- phwe <= hwe
        
        # check sides
        # consider sides
        if(hwe_excess_side == "both"){
          side <- 1
        }
        else if(hwe_excess_side == "heterozygote"){
          side <- ifelse(x@stats$fis[x@stats$facet == ".base"] <= 0, 1, 0)
        }
        else{
          side <- ifelse(x@stats$fis[x@stats$facet == ".base"] >= 0, 1, 0)
        }
        
        if(any(is.na(side))){
          side[is.na(side)] <- 1
        }
        
        # recheck low.p
        phwe <- phwe * side
        phwe <- which(phwe == 1)
        
        if(verbose){cat("\t", length(phwe), " bad loci\n")}
        vio.snps[phwe] <- T
        
      }
      
      # facets, slightly more complicated
      else{
        # run again to ensure that the correct fwe method and case are used
        .make_it_quiet(x <- calc_hwe(x, hwe_facets, fwe_method = fwe_method[1], fwe_case = "by_facet"))
        if(hwe_excess_side != "both"){
          calced <- .check_calced_stats(x, hwe_facets, c("fis"))
          if(any(!unlist(calced))){
            .make_it_quiet(x <- calc_fis(x, names(calced)[!unlist(calced)]))
          }
        }
        
        # get the per-facet hwe stats, check against threshold, then condense by snp.
        phwe <- get.snpR.stats(x, hwe_facets)
        if(fwe_method != "none"){
          phwe$low.p <- ifelse(phwe[,paste0("pHWE_byfacet_", fwe_method[1])] <= hwe, 1, 0)
        }
        else{
          phwe$low.p <- ifelse(phwe$pHWE <= hwe, 1, 0)
        }
        
        # consider sides
        if(hwe_excess_side == "both"){
          phwe$side <- 1
        }
        else if(hwe_excess_side == "heterozygote"){
          phwe$side <- ifelse(phwe$fis <= 0, 1, 0)
        }
        else{
          phwe$side <- ifelse(phwe$fis >= 0, 1, 0)
        }
        
        if(any(is.na(phwe$side))){
          phwe$side[is.na(phwe$side)] <- 1
        }
        
        # recheck low.p
        phwe$low.p <- phwe$low.p * phwe$side
        
        bad.loci <- tapply(phwe$low.p, phwe[,".snp.id"], sum, na.rm = TRUE)
        bad.loci <- reshape2::melt(bad.loci)
        bad.loci <- stats::na.omit(bad.loci)
        bad.loci <- bad.loci[which(bad.loci$value > 0),]
        bad.loci <- which(snp.meta(x)$.snp.id %in% bad.loci$Var1)
        if(verbose){cat("\t", length(bad.loci), " bad loci\n")}
        
        vio.snps[bad.loci] <- TRUE
      }
      
      x <- .update_filters(x, "Hardy-Weinburg Proportions--hwe", 
                           paste0(hwe, ", excess side = ", hwe_excess_side),
                           ifelse(is.null(hwe_facets) | isFALSE(hwe_facets), ".base", hwe_facets))
    }
    
    #==========remove violating loci==================
    if(any(vio.snps)){
      if(sum(vio.snps) == nrow(x)){
        stop("No loci passed filters.\n")
      }
      x <- x[-which(vio.snps),]
    }
    return(x)
  }
  
  #funciton to filter by individuals.
  min_loci_filt <- function(){
    if(verbose){cat("Filtering out individuals sequenced in few kept loci...\n")}
    mcounts <- matrixStats::colSums2(ifelse(x != mDat, 1, 0))
    rejects <- which(mcounts/nrow(x) < min_loci)
    if(length(rejects) > 0){
      if(length(rejects) == ncol(x)){
        stop("No individuals passed filters.\n")
      }
      
      old.facets <- x@facets
      x <- x[,-rejects]
      if(verbose){cat("Re-calculating and adding facets.\n")}
      if(any(old.facets != ".base")){
        x <- .add.facets.snpR.data(x, old.facets[-which(old.facets == ".base")])
      }
    }
    
    x <- .update_filters(x, "poorly sequenced individuals--min_loci", min_loci, NA)
    
    return(list(x = x, rejects = rejects))
  }
  
  #==========================call the functions as requested.==================
  if(is.numeric(remove_garbage)){
    if(verbose){cat("Removing garbage individuals/loci.\n")}
    mcounts <- matrixStats::colSums2(ifelse(x != mDat, 1, 0))
    rejects <- which(mcounts/nrow(x) < remove_garbage)
    
    mi <- x@geno.tables$wm[,colnames(x@geno.tables$wm) == mDat] # for backwards compat
    mi <- which((nrow(x@sample.meta) - mi)/nrow(x@sample.meta) < remove_garbage)
    
    if(length(rejects) > 0){
      if(length(mi) > 0){
        .make_it_quiet(x <- x[-mi,-rejects])
      }
      else{
        .make_it_quiet(x <- x[,-rejects])
      }
    }
    else if(length(mi) > 0){
      .make_it_quiet(x <- x[-mi,])
    }
    
    
    if(verbose){cat("\tRemoved", length(mi), "bad loci.\n\tRemoved", length(rejects), "bad individuals.\n")}
    
    if(nrow(x) == 0){
      stop("No SNPs passed garbage removal.\n")
    }
    if(ncol(x) == 0){
      stop("No loci passed garbage removal.\n")
    }
    
    rm(mi, mcounts, rejects)
    
    x <- .update_filters(x, "poorly sequenced individuals and loci jointly--remove_garbage",
                         remove_garbage, NA)
  }
  
  if(!inds_first){
    if(any(c(non_poly, bi_al, maf, hf_hets, min_ind) != FALSE)){
      if(verbose){cat("Filtering loci. Starting loci:", nrow(x), "\n")}
      
      # run the filter
      x <- filt_by_loci()
      
      if(nrow(x) == 0 | is.null(nrow(x))){
        stop("No loci remain after filters.")
      }
      
      if(verbose){cat("\tEnding loci:", nrow(x), "\n")}
    }
  }
  
  
  # run the minimum sequenced loci filter
  if(min_loci){
    if(verbose){cat("Filtering individuals. Starting individuals:", ncol(x), "\n")}
    x <- min_loci_filt()
    if(length(x$rejects) == 0){
      if(verbose){cat("No individuals removed.\n")}
      x <- x$x
    }
    else{
      x <- x$x
      if(verbose){cat("\tEnding individuals:", ncol(x), "\n")}
    }
  }
  
  if(inds_first){
    if(any(c(non_poly, bi_al, maf, hf_hets, min_ind) != FALSE)){
      if(verbose){cat("Filtering loci. Starting loci:", nrow(x), "\n")}
      
      # run the filter
      x <- filt_by_loci()
      
      if(nrow(x) == 0 | is.null(nrow(x))){
        stop("No loci remain after filters.")
      }
      
      if(verbose){cat("\tEnding loci:", nrow(x), "\n")}
    }
  }
  
  
  if(re_run != FALSE){
    if(!inds_first){
      if(verbose){cat("Re-filtering loci...\n")}
      
      if(re_run == "partial"){
        maf <- FALSE
        hf_hets <- FALSE
        min_ind <- FALSE
        bi_al <- FALSE
        hwe <- FALSE
        mac <- 0
        singletons <- FALSE
      }
      if(any(c(non_poly, bi_al, maf, hf_hets, min_ind) != FALSE)){
        x <- filt_by_loci() # re-filter loci to make sure that we don't have any surprise non-polys etc.
        if(verbose){cat("\tFinal loci count:", nrow(x), "\n")}
      }
      else{
        if(verbose){cat("\tNo variables to re-fitler.\n")}
      }
    }
    else{
      if(verbose){cat("Re-Filtering individuals. Starting individuals:", ncol(x), "\n")}
      x <- min_loci_filt()
      if(length(x$rejects) == 0){
        if(verbose){cat("No individuals removed.\n")}
        x <- x$x
      }
      else{
        x <- x$x
        if(verbose){cat("\tEnding individuals:", ncol(x), "\n")}
      }
    }
    
  }
  
  # return results
  return(x)
}

#'Re-format SNP data.
#'
#'\code{format_snps} re-formats SNP data into a range of different possible
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
#'compatibility and internal use.
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
#'with allele frequency in all samples or each population.} \item{rafm: }{RAFM
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
#'BirthYear (or instead of BirthYear - BYmin and BYmax), optional column 
#'Yearlast in sample metadata for running Sequoia. For more information see 
#'sequoia documentation.} \item{fasta:} {fasta sequence format.} \item{vcf:}
#'{Variant Call Format, a standard format for SNPs and other genomic variants. 
#'Genotypes are coded as 0/0, 0/1, 1/1, or ./. (for missing values), with a 
#'healthy serving of additional metadata but very little sample metadata.} 
#'\item{genalex: }{GenAlEx format. If an outfile is requested, the data will be 
#'sorted according to any provided facets and written as an '.xlsx' object.} 
#'\item{snpRdata: }{a snpRdata object.} }
#'
#'Note that for the "sn" format, the data can be interpolated to fill missing
#'data points, which is useful for PCA, genomic prediction, tSNE, and other
#'methods. To do so, specify interpolate = "af" to insert the expected number of
#'minor alleles given SNP allele frequency or "bernoulli" to do binomial draws
#'to determine the number of minor alleles at each missing data point, where the
#'probability of drawing a minor allele is equal to the minor allele frequency.
#'The expected number of minor alleles based on the later method is equal to the
#'interpolated value from the former, but the later allows for multiple runs to
#'determine the impact of stochastic draws and is generally preferred and
#'required for some downstream analysis. It is therefore the default. As a
#'slower but more accurate alternative to "af" interpolation, "iPCA" may be
#'selected. This an iterative PCA approach to interpolate based on SNP/SNP
#'covariance via \code{\link[missMDA]{imputePCA}}. If the ncp argument is not
#'defined, the number of components used for interpolation will be estimated
#'using \code{\link[missMDA]{estim_ncpPCA}}. In this case, this method is much
#'slower than the other methods, especially for large datasets. Setting an ncp
#'of 2-5 generally results in reasonable interpolations without the time
#'constraint.
#'
#'Note also that for the plink format, a .bed binary file can be generated. If
#'the "plink" option is selected and an outfile is designated, R will generate a
#'".sh" shell file with the same name given in the outfile argument. Running
#'this file will create a plink.bed file.
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
#'@param output Character, default "snpRdata". The desired output format, see
#'  description for details.
#'@param facets Character or NULL, default NULL. Facets over which to break up
#'  data for some output formats, following the format described in
#'  \code{\link{Facets_in_snpR}}.
#'@param n_samp Integer or numeric vector, default NA. For structure or RAFM
#'  outputs. How many random loci should be selected? Can either be an integer
#'  or a numeric vector of loci to use.
#'@param interpolate Character or FALSE, default "bernoulli". If transforming to
#'  "sn" or "pa" format, notes the interpolation method to be used to fill
#'  missing data. Options are "bernoulli", "af", "iPCA", or FALSE. See details.
#'@param outfile character vector, default FALSE. If a file path is provided, a
#'  copy of the output will be saved to that location. For some output styles,
#'  such as genepop, additional lines will be added to the output to allow them
#'  to be immediately run on commonly used programs.
#'@param ped data.frame default NULL. Optional argument for the "plink" output
#'  format. A six column data frame containing Family ID, Individual ID,
#'  Paternal ID, Maternal ID, Sex, and Phenotype and one row per sample. If
#'  provided, outputs will contain information contained in ped. See plink
#'  documentation for more details.
#'@param input_format Character, default NULL. Format of x, by default a
#'  snpRdata object. See description for details.
#'@param input_meta_columns Numeric, default NULL. If x is not a snpRdata
#'  object, optionally specifies the number of metadata columns preceding
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
#'@param plink_recode_numeric Logical, default FALSE. If FALSE, all chrs/scaffs
#'  will be renamed to numbers. This may be useful in some cases. If this is
#'  FALSE, chromosome names will be checked for leading numbers and replaced
#'  with a corresponding letter (0 becomes A, 1 becomes B, and so on).
#'@param verbose Logical, default FALSE. If TRUE, some progress updates will be
#'  reported.

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
#' \dontrun{
#' #import data to a snpRdata object
#' ## get sample meta data
#' sample_meta <- 
#'     data.frame(pop = substr(colnames(stickRAW)[-c(1:3)], 1, 3), 
#'                fam = rep(c("A", "B", "C", "D"), 
#'                          length = ncol(stickRAW) - 3), 
#'                stringsAsFactors = FALSE)
#' format_snps(stickRAW, input_format = "0000", 
#'             input_meta_columns = 3, 
#' input_mDat = "0000", sample.meta = sample_meta)
#'
#' #allele count, separated by the pop facet.
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
#' #Sequoia format
#' b <- sample.meta(stickSNPs)
#' b$ID <- 1:nrow(b)
#' b$Sex <- rep(c("F", "M", "U"), length.out=nrow(b))
#' b$BirthYear <- round(runif(n = nrow(b), 1,1))
#' a <- stickSNPs
#' b$ID <- paste0(sample.meta(a)$pop, sample.meta(a)$fam, sample.meta(a)$.sample.id)
#' sample.meta(a) <- b
#' format_snps(x = a, output = "sequoia")
#' #note: if using the birth year windows BYmin and BYmax and or Yearlast, 
#' ensure that the column names are not stored as BY.min, BY.max, Year.last for
#' snpR. 
#'
#' # VCF format
#' test <- format_snps(stickSNPs, "vcf", chr = "chr")
#' 
#' # GenAlEx format (write to file to generate a facet-sorted xlsx file)
#' test <- format_snps(stickSNPs, "genalex")
#' }
format_snps <- function(x, output = "snpRdata", facets = NULL, n_samp = NA,
                        interpolate = "bernoulli", outfile = FALSE,
                        ped = NULL, input_format = NULL,
                        input_meta_columns = NULL, input_mDat = NULL,
                        sample.meta = NULL, snp.meta = NULL, chr.length = NULL,
                        ncp = 2, ncp.max = 5, chr = "chr", position = "position",
                        phenotype = "phenotype", plink_recode_numeric = FALSE, 
                        verbose = FALSE){
  if(!isTRUE(verbose)){
    cat <- function(...){}
  }
  
  #======================sanity checks================
  if(!is.null(input_format)){
    if(tolower(input_format) == "snprdata"){input_format <- NULL}
  }
  
  # check that a useable output format is given. keming
  output <- tolower(output) # convert to lower case.
  if(output == "nn"){output <- "NN"}
  pos_outs <- c("ac", "genepop", "structure", "0000", "hapmap", "NN", "pa",
                "rafm", "faststructure", "dadi", "plink", "sn", "snprdata",
                "colony","adegenet", "fasta", "lea", "sequoia", "vcf", "genalex")
  if(!(output %in% pos_outs)){
    stop("Unaccepted output format specified. Check documentation for options.\n")
  }
  
  # check that provided snpRdata objects are in the correct format
  if(is.null(input_format)){
    if(!is.snpRdata(x)){
      stop("If x is not a snpRdata object, provide input data format.\n")
    }
  }
  else if(is.snpRdata(x)){
    imput_format <- NULL
  }
  
  # bi-allelic check
  if(is.snpRdata(x)){
    if(!.is.bi_allelic(x)){
      pos_nbi_outs <- c("genepop", "0000", "pa")
      if(!output %in% pos_nbi_outs){
        stop(paste0("The output format you selected is not currently supported for non-biallelic data. Currently supported formats are:",
                    paste0(pos_nbi_outs, collapse = ", "), ".\n"))
      }
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
    if(is.null(chr)){
      stop("The chr argument must be provided for plink output, indicating which column of the data contains chromosome information.\n")
    }
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
  else if(output == "genalex"){
    cat("Converting to genalex format.")
    .check.installed("openxlsx")
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
      x <- .process_ms(x, chr.length)
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
    
    if(!.is.bi_allelic(x)){
      if(output == "genepop"){
        rdata <- as.data.frame(t(genotypes(x)))
        row.names(rdata) <- paste0(row.names(rdata), " ,") #adding space and comma to row names, as required.
      }
      else{
        rdata <- genotypes(x)
        rdata <- cbind(x@snp.meta, rdata)
        colnames(rdata)[1:ncol(x@snp.meta)] <- colnames(x@snp.meta)
      }
    }
    
    else{
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
  }
  
  
  ##convert to structure, fastStructure or RAFM format (v)
  if (output == "structure" | output == "rafm" | output == "faststructure" | output == "genalex"){
    if(length(facets) > 1){
      stop("Only one facet at a time permitted for structure/rafm/faststructure/genalex!\n")
    }
    facets <- .check.snpR.facet.request(x, facets)
    
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
    
    # transpose, since these are sample based, then split into two matrices
    xv <- t(as.matrix(x))
    xv1 <- substr(xv, 1, 1)
    xv2 <- substr(xv, 2, 2)
    
    # replace
    xv1[xv1 == "A"] <- 1
    xv1[xv1 == "C"] <- 2
    xv1[xv1 == "G"] <- 3
    xv1[xv1 == "T"] <- 4
    xv1[xv1 == substr(x@mDat, 1, nchar(x@mDat)/2)] <- 0
    
    
    xv2[xv2 == "A"] <- 1
    xv2[xv2 == "C"] <- 2
    xv2[xv2 == "G"] <- 3
    xv2[xv2 == "T"] <- 4
    xv2[xv2 == substr(x@mDat, 1, nchar(x@mDat)/2)] <- 0
    
    
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
    
    else if(output == "genalex"){
      rdata <- matrix(NA, nrow = nsamps(x), ncol = nsnps(x)*2)
      rdata[,seq(1, ncol(rdata), by = 2)] <- xv1
      rdata[,seq(2, ncol(rdata), by = 2)] <- xv2
      rdata <- matrix(as.numeric(rdata), nrow(rdata), ncol(rdata))
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
      if(any(colnames(x@sample.meta) %in% facets)){
        colnames(rdata)[1:sum(colnames(x@sample.meta) %in% facets)] <- colnames(x@sample.meta)[colnames(x@sample.meta) %in% facets]
      }
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
      if(any(grepl("_", as))){stop("Detected '_' in genotypes. Alleles cannot be coded with underscores for pa conversion.\n")}
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
      if(.is.bi_allelic(x)){ # fully vectorized
        xmc <- which(x == mDat) #which samples had missing data?
        adj <- floor(xmc / nsamp) #how many loci over do I need to adjust xmc, since in amat each locus occupies two columns?
        adj[xmc%%nsamp == 0] <- adj[xmc%%nsamp == 0] - 1 #shift over anything that got screwed up by being in the last sample
        xmc <- xmc + (nsamp*adj) #adjust xmc for extra columns.
        if(any(amat[xmc] != 0) | any(amat[xmc + nsamp] != 0)){
          stop("Missing data values were not properly identified for replacement with NAs. This usually happens when SNP data is not completely bi-allelic. Try filtering out non-biallelic and non-polymorphic SNPs using filter_snps.\n")
        }
        amat[xmc] <- NA #make the first allele NA
        amat[xmc + nsamp] <- NA #make the second allele (another column over) NA.
      }
      else{ # need to loop through loci :(
        loc_key <- gsub("_.+", "", colnames(amat))
        loc_key <- table(loc_key)
        loc_key <- loc_key[order(as.numeric(names(loc_key)))]
        loc_key <- cumsum(loc_key)
        loc_key <- c(0, loc_key)
        
        xmc <- which(as.matrix(x) == x@mDat, arr.ind = TRUE) #which samples had missing data?
        
        for(j in 1:(length(loc_key) - 1)){
          tl <- xmc[which(xmc[,1] == j),]
          if(sum(amat[tl[,2], (loc_key[j] + 1):loc_key[j + 1]]) != 0){stop("Failed to correctly fill NAs.\n")}
          amat[tl[,2], loc_key[j]:loc_key[j + 1]] <- NA
        }
      }
      
      return(amat)
    }
    
    amat <- pa_alleles(t(xv), x@snp.form, x@mDat)
    
    # interpolate?
    if(interpolate == "bernoulli"){
      amat <- t(.interpolate_sn(t(amat), "bernoulli"))
    }
    else if(interpolate == "af"){
      amat <- t(.interpolate_sn(t(amat), "af"))
    }
    else if(interpolate == "iPCA"){
      amat <- t(.interpolate_sn(t(amat), "iPCA", ncp = ncp, ncp.max = ncp.max))
    }
    
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
    
    # check for bad chr names
    if(!plink_recode_numeric){
      chrs <- summarize_facets(x, chr)$chr
      lead <- unlist(lapply(chrs, substr, stop = 1, start = 1))
      bad_chrs <- which(lead %in% 0:9)
      if(length(bad_chrs) > 0){
        lead[bad_chrs] <- LETTERS[1:10][match(lead[bad_chrs], 0:9)]
        n.chrs <- snp.meta(x)[,chr]
        if(length(bads) > 0){
          n.chrs <- n.chrs[-bads]
        }
        substr(n.chrs[which(n.chrs %in% chrs[bad_chrs])], 1, 1) <- lead[match(n.chrs[which(n.chrs %in% chrs[bad_chrs])], chrs)]
        warning("Found numbers at the start of chr/scaffold names, which plink does not allow. Replaced with letters (0:9 = A:J).\n")
      }
      else{
        n.chrs <- snp.meta(x)[,chr]
        if(length(bads) > 0){
          n.chrs <- n.chrs[-bads]
        }
      }
    }
    else{
      n.chrs <- snp.meta(x)[,chr]
      if(length(bads) > 0){
        n.chrs <- n.chrs[-bads]
      }
    }
    
    
    # sort by chr then position
    if(length(bads) > 0){
      n.ord <- order(n.chrs, snp.meta(x)$position[-bads])
    }
    else{
      n.ord <- order(n.chrs, snp.meta(x)$position)
    }
    n.chrs <- n.chrs[n.ord]
    sn <- sn[n.ord,]
    
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
        if(length(unique(stats::na.omit(Phenotype))) == 2){
          uf <- factor(Phenotype)
          cat("Two phenotypes detected, options encoded as: \n\t",
              levels(uf)[1], " = 0\n\t",
              levels(uf)[2], " = 1\n\t",
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
        if(length(unique(stats::na.omit(Phenotype))) == 2){
          uf <- factor(Phenotype)
          cat("Two phenotypes detected, options encoded as: \n\t",
              levels(uf)[1], " = 0\n\t",
              levels(uf)[2], " = 1\n\t",
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
    
    
    # remove bads
    if(length(bads) > 0){
      bim <- bim[-bads,]
    }
    
    # re-order
    bim <- bim[n.ord,]
    
    # re-name
    bim$chr <- n.chrs
    
    # recode chr
    if(plink_recode_numeric){
      bim$chr <- as.numeric(as.factor(bim$chr))
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
        #format sex for sequoia (1 = F,2 = M,3 = U,4 = H, NA)
        sexes <- c("F", "M", "U", "H", "1", "2", "3", "4")
        sex <- sample.meta(x)$Sex
        sex[which(!sex %in% sexes)] <- 3
        sex[which(sex == "F")] <- 1
        sex[which(sex == "M")] <- 2
        sex[which(sex == "U")] <- 3
        sex[which(is.na(sex) == TRUE)] <- 3
        sex[which(sex == "H")] <- 4
        return(sex)
      }
      
      #for lifehistory data input needs id, sex, year born
      if(!all(c("Sex", "BirthYear", "ID", "BYmin", "BYmax", "Yearlast") %in% colnames(x@sample.meta))){
          warning("Required columns Sex, ID, Yearlast, BirthYear, BYmax, and BYmin not found in sample metadata.\n")
        }
        
        ID <- x@sample.meta$ID
        sex <- fix.sex.sequoia(x)
        
        #ACTUALLY MAKE THE TABLE
        lhtable <- data.frame(ID=ID,
                              Sex = as.numeric(sex),
                              BirthYear = sample.meta(x)$BirthYear,
                              BY.min = sample.meta(x)$BYmin,
                              BY.max = sample.meta(x)$BYmax,
                              Year.last = sample.meta(x)$Yearlast,
                              stringsAsFactors = F)
      
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
    r <- try(x@filters, silent = TRUE)
    if(methods::is(r, "try-error")){
      .make_it_quiet(vcf_meta <- c("##fileformat=VCFv4.0",
                                   paste0("##fileDate=", gsub("-", "", Sys.Date())),
                                   "##source=snpR",
                                   '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
                                   '##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele Count">',
                                   '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'))
    }
    else{
      .make_it_quiet(vcf_meta <- c("##fileformat=VCFv4.0",
                                   paste0("##fileDate=", gsub("-", "", Sys.Date())),
                                   "##source=snpR",
                                   '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
                                   '##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele Count">',
                                   '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                                   paste0("##snpR filter_snps=", filters(x))))
    }
    
    # meta columns
    data_meta <- data.frame(CHROM = snp.meta(x)[,chr],
                            POS = snp.meta(x)[,position],
                            ID = paste0("SNP_", 1:nrow(x)),
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
    
    # reorder by chr then pos
    rdata <- dplyr::arrange(rdata, `#CHROM`, `POS`)
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
      base::cat(paste0(unlist(.split.facet(outfile))[1], "_genepop\n"), file = outfile)
      base::cat(llist, "\nPOP\n", file = outfile, append = T) 
      
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
            base::cat("POP\n", file = outfile, append = T)
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
        base::cat("[loci]=", nrow(x), "\n\n[populations]=", length(u.pops), "\n\n", file = outfile, sep = "")
        
        #write the data for each population.
        for(i in 1:length(u.pops)){
          base::cat("[pop]=", i, "\n", file = outfile, append = T, sep = "") #write header
          
          tdat <- trdat[trdat$subfacet == u.pops[i],]
          
          wdat <- cbind(snp = tdat$.snp.id, tdat[,3:6])
          
          data.table::fwrite(wdat,
                             outfile, col.names = F, row.names = F, quote = F, sep = "\t",
                             append = T) # write the data for this population.
          
          base::cat("\n", file = outfile, append = T) # add a line break
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
    else if(output == "genalex"){
      wb <- openxlsx::createWorkbook(creator = "snpR")
      openxlsx::addWorksheet(wb, sheetName = "Data")
      
      if(length(facets) > 0){
        
        # figure out samples in each facet
        levs <- .paste.by.facet(sample.meta(x), unlist(.split.facet(facets)))
        ulevs <- unique(levs)
        matches <- vector("list", length(ulevs))
        names(matches) <- ulevs
        for(i in 1:length(ulevs)){
          matches[[i]] <- .fetch.sample.meta.matching.task.list(x, c(facets, ulevs[i], ".base", ".base"))
        }
        
        # make the header
        header <- c(nsnps(x), nsamps(x), length(matches), unlist(lapply(matches, length)))
        names(header) <- NULL
        openxlsx::writeData(wb, "Data", matrix(header, 1), colNames = FALSE)
        openxlsx::writeData(wb, "Data", matrix(c("snpR exported data", "", "", ulevs), 1), startRow = 2, colNames = FALSE)
        
        # make the column titles
        col_titles <- c("SAMPLE", "POP", c(rbind(paste0("SNP", 1:nsnps(x)), "")))
        openxlsx::writeData(wb, "Data", matrix(col_titles, nrow = 1), startRow = 3, colNames = FALSE)
        
        # write the data
        prog <- 4
        for(i in 1:length(matches)){
          part <- (prog-3):(prog-4 + length(matches[[i]]))
          openxlsx::writeData(wb, "Data", cbind(paste0("Sample", part), levs[matches[[i]]]), colNames = FALSE, startRow = prog)
          openxlsx::writeData(wb, "Data", rdata[matches[[i]],], colNames = FALSE, startRow = prog, startCol = 3)
          prog <- prog + length(matches[[i]])
        }
      }
      else{
        # header
        header <- c(nsnps(x), nsamps(x), 1)
        # make the header
        openxlsx::writeData(wb, "Data", matrix(header, 1), colNames = FALSE)
        openxlsx::writeData(wb, "Data", matrix(c("snpR exported data", "", "", "base"), 1), startRow = 2, colNames = FALSE)
        
        # make the column titles
        col_titles <- c("SAMPLE", "POP", c(rbind(paste0("SNP", 1:nsnps(x)), "")))
        openxlsx::writeData(wb, "Data", matrix(col_titles, nrow = 1), startRow = 3, colNames = FALSE)
        
        # data
        openxlsx::writeData(wb, "Data", cbind(paste0("Sample", 1:nsamps(x)), "base"), colNames = FALSE, startRow = 4)
        openxlsx::writeData(wb, "Data", rdata, colNames = FALSE, startRow = 4, startCol = 3)
      }
      
      openxlsx::saveWorkbook(wb, file = outfile, overwrite = TRUE)
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
#' @param verbose logical, default FALSE. If TRUE, prints detailed progress 
#'   report.
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
#' \dontrun{
#' # check for duplicates with sample 1
#' check_duplicates(stickSNPs, 1)
#'
#' # check duplicates using the .samp.id column as sample IDs
#' check_duplicates(stickSNPs, 1, id.col = ".sample.id")
#' }
check_duplicates <- function(x, y = 1:ncol(x), id.col = NULL, verbose = FALSE){
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
    if(verbose){
      cat("Working on sample:", y[i], "\n")
    }
    
    # pick out this value
    if(!is.null(id.col)){
      t.samp.id <- which(x@sample.meta[,id.col] == y[i])
    }
    else{
      t.samp.id <- y[i]
    }
    if(length(t.samp.id) == 0){
      out.best$matches[i] <- "bad.ID"
      if(verbose){cat("\tBad ID, skipping.\n")}
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
    
    if(verbose){
      cat("\tBest hit(s):", out.best$best_match[i], "\n")
      cat("\tPercentage match(es):", out.best$percentage[i], "\n")
    }
  }
  
  # return
  return(list(best_matches = out.best, data = out))
  
}



#' Fetch the allele frequencies for all SNPs for each level of each requested
#' facet.
#'
#' Fetch allele frequencies for all SNPs for each level of all the requested
#' facets. Major and minor allele frequencies will be interleaved, with the major
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
    
    # interleave major and minor frequencies and save
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
#'   scaffolds. SNP positions will be independently considered depending on the
#'   facet level.
#' @param n Integer. Specifies the minimum distance between selected SNPs.
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
#' snpR will automatically track the methods used for calculations on a
#' snpRdata object. Using \code{\link{citations}} on that object will provide
#' details on the the methods used, and can optionally write a .bib bibtex
#' formatted library containing citations for these methods.
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
#'   to \code{\link[rbibutils]{readBib}}.
#'
#' @author William Hemstrom
#' @returns  If return_bib is TRUE, a list containing four parts:
#'   \itemize{\item{keys: } A vector of bibtex keys for each method.
#'   \item{stats: } A vector of the stats used. \item{details: } A vector of
#'   details for each method. \item{bib: } A \code{bibentry} for each citation,
#'   see \code{\link[rbibutils]{readBib}}.}
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
  
  .check.installed("rbibutils")
  
  #==========grab bib============
  bib.file <- system.file("extdata", "snpR_citations.bib", package = "snpR")
  bib <- rbibutils::readBib(bib.file)
  
  #==========filter bib==========
  keys <- as.character(unlist(purrr::map(x@citations, "key")))
  bib <- bib[keys]
  
  #==========shout at the user=====
  deets <- unlist(purrr::map(x@citations, "details"))
  stats <- names(x@citations)
  
  cat("Citations for methods used thus far: \n")
  for(i in 1:length(keys)){
    cat("==============================================\n")
    cat("Statistic: ", stats[i], "\n\n")
    cat("Citation: ", .quick_grab_cite(bib[keys[i]]), "\n\n")
    cat("Bibtex key: ", keys[i], "\n\n")
    cat("Details: ", deets[i], "\n")
  }
  cat("==============================================\n\n")
  
  #==========print bib=============
  if(!isFALSE(outbib)){
    
    if(file.exists(outbib)){
      current_bib <- rbibutils::readBib(outbib)
      not_in_current <- which(!keys %in% names(current_bib))
      
      
      if(length(not_in_current) > 0){
        #==========filter bib==========
        bib <- bib[keys[not_in_current]]
        
        current_bib <- c(current_bib, bib)
        
        rbibutils::writeBib(current_bib, outbib)
      }
    }
    
    
    else{
      rbibutils::writeBib(bib, outbib)
      cat(".bib file can be found at: ", outbib, "\n")
    }
    
  }
  
  if(!isFALSE(return_bib)){
    return(list(keys = keys, stats = stats, details = deets, bib = bib))
  }
}

#' Report on filters used on a snpRdata object.
#' 
#' snpR will automatically track the methods used on a
#' snpRdata object. Using \code{\link{filters}} on that object will provide
#' details on the the filters used.
#'
#' Printed outputs contain the filters used in the order they were applied
#' alongside the filtering stringency and any facets applied over, if 
#' applicable. Note that some output formats from \code{format_snps}, like the
#' \code{vcf} format, will notes on the filters used as well.
#'
#' @param x snpRdata object
#'
#' @author William Hemstrom
#' 
#' @return Prints easily human readable results to the console and also returns
#' a more machine readable string.
#' 
#' @export
#' @examples 
#' # filter the data
#' x <- filter_snps(stickSNPs, 
#'                  min_ind = 0.75, 
#'                  min_loci = 0.75,
#'                  maf = 0.1,
#'                  maf_facets = "pop")
#' 
#' # fetch the filters used
#' filters(x)
filters <- function(x){
  r <- try(x@filters, silent = TRUE)
  if(methods::is(r, "try-error")){
    cat("This is an old `snpRdata` object that did not track filters.\n")
  }
  else{
    line <- character(0)
    for(i in 1:nrow(r)){
      cat("Filter:", r[[1]][i], "\n")
      line_part <- unlist(strsplit(r[[1]][i], "--"))
      line_part <- line_part[length(line_part)]
      
      
      if(!is.na(r[[2]][i])){
        cat("\tStringency:", r[[2]][i], "\n")
        line_part <- c(line_part, paste0("=", r[[2]][i]))
      }
      
      
      if(!is.na(r[[3]][i])){
        cat("\tFacet:", r[[3]][i], "\n")
        line_part <- c(line_part, paste0(",facet=", r[[3]][i]))
      }
      cat("\n")
      
      line <- c(line, paste0(line_part, collapse = ""))
    }
  }
  
  return(paste0(line, collapse = ";"))
}

#' Summarize possible snpRdata object facet options
#' 
#' List either all of the possible SNP and sample facets (if called with no facets)
#' or all of the categories for each requested facet.
#' 
#' @param x snpRdata object
#' @param facets character. Categorical metadata variables by which to break up
#'  analysis. See \code{\link{Facets_in_snpR}} for more details. If NULL, the 
#'  possible SNP and sample facets will be listed. If facets are instead provided,
#'  the \emph{categories} for each facet will instead be listed.
#' 
#' @export
#' @author William Hemstrom
#' @return A named list containing either the possible SNP and sample facets or
#'   the categories for all of the requested facets.
#' @examples 
#' # list available facets
#' summarize_facets(stickSNPs)
#' 
#' # return details for a few facets
#' summarize_facets(stickSNPs, c("pop", "chr.pop", "fam.pop"))
summarize_facets <- function(x, facets = NULL){
  if(!is.snpRdata(x)){
    stop("x must be a snpRdata object.\n")
  }
  
  facets <- .check.snpR.facet.request(x, facets, remove.type = "none", fill_with_base = FALSE, return_base_when_empty = FALSE, return.type = TRUE)
  base_facets <- which(facets[[1]] == ".base")
  if(length(base_facets) > 0){
    facets[[1]] <- facets[[1]][-base_facets]
    facets[[2]] <- facets[[2]][-base_facets]
  }
  
  
  # If called without facets, list availabilities
  if(length(facets[[1]]) == 0){
    message("Returning list of facets. For more details, try asking for information about a specific facet (such as 'pop' or 'pop.chr')!")
    
    return(list(SNP = colnames(snp.meta(x)),
                sample = colnames(sample.meta(x))))
  }
  
  # otherwise return options for each facet
  else{
    
    
    out <- vector("list", length(facets[[1]]))
    names(out) <- facets[[1]]
    for(i in 1:length(out)){
      if(facets[[2]][i] == "snp"){
        out[[i]] <- .get.task.list(x, facets[[1]][i])[,4]
      }
      else if(facets[[2]][i] == "sample"){
        out[[i]] <- .paste.by.facet(sample.meta(x), unlist(.split.facet(facets[[1]][i])))
      }
      else if(facets[[2]][i] == "complex"){
        split_facet <- unlist(.split.facet(facets[[1]][i]))
        samp_part <- .check.snpR.facet.request(x, split_facet, fill_with_base = FALSE, return_base_when_empty = FALSE)
        snp_part <- .check.snpR.facet.request(x, split_facet, remove.type = "sample", fill_with_base = FALSE, return_base_when_empty = FALSE)
        samp_part <- sample.meta(x)[,samp_part, drop = FALSE]
        samp_part <- lapply(samp_part, unique)
        snp_part <- snp.meta(x)[,snp_part, drop = FALSE]
        snp_part <- lapply(snp_part, unique)
        
        opts <- expand.grid(c(samp_part, snp_part))
        out[[i]] <- .paste.by.facet(opts, split_facet)
      }
    }
  }
  
  out <- lapply(out, unique)
  
  return(out)
}

#' Merge two snpRdata objects
#' 
#' Merge two snpRdata objects using sample and SNP metadata. Functions
#' much like base R's  \code{\link[base]{merge}} function, but
#' the 'by' and 'all' options can be specified at the SNP and sample level.
#' 
#' While this function can be used essentially identically to how one might
#' use base R's \code{\link[base]{merge}} function, there are a few differences
#' to note. 
#' 
#' First, samples that are genotyped at identical loci 
#' in both data sets can be handled several ways, controlled by the
#' \code{resolve_conflicts} argument: \itemize{\item{warning:} Return a harsh 
#' warning and a data frame with more information on genotypes at identical 
#' samples/SNPs are different between \code{x} and \code{y}. 
#' \item{error: } The default, return an error when conflicts are detected.
#' \item{x} Use genotypes from \code{x} to resolve conflicts.
#' \item{y} Use genotypes from \code{y} to resolve conflicts.
#' \item{random} Randomly sample (non-missing) genotypes from \code{x}
#' and \code{y} to resolve conflicts.}
#' Note that called genotypes are always taken over un-called genotypes when
#' there are merge conflicts, and missing data in one but not the other data set
#' will not trigger an error or a warning if those options are selected.
#' 
#' Secondly, the \code{by} and \code{all} arugment families from 
#' \code{\link[base]{merge}} are extended to refer to either samples or SNPs,
#' such that all samples can be maintained but not all SNPs, for example.
#' 
#' Lastly, all of the \code{all} family of arguments default to \code{TRUE}
#' instead of \code{FALSE}, since purely overlapping genotypes/SNPs is unlikely
#' to be desired. \code{FALSE} values provided to any specific \code{all}
#' argument will sill override \code{all = TRUE}, as in 
#' \code{\link[base]{merge}}.
#' 
#' At present, \code{\link{merge_snpRdata}} is not maximally efficient in that
#' it will remove all tabulated statistics and re-tabulate all internal 
#' summaries. Improvements are in development.
#' 
#' @param x,y \code{snpRdata} objects to merge
#' @param by.sample,by.sample.x,by.sample.y Columns of sample metadata by which
#'   to merge across samples--function identically to the \code{by}, \code{by.x},
#'   and \code{by.y} arguments to \code{\link[base]{merge}}, see documentation
#'   there for details.
#' @param by.snp,by.snp.x,by.snp.y Columns of SNP metadata by which to merge
#'   across SNPs--function idetically to the \code{by}, \code{by.x}, and
#'   \code{by.y} arguments to \code{\link[base]{merge}}, see documentation there
#'   for details.
#' @param all logical, default TRUE. If TRUE, all samples and SNPs will be 
#'   maintained in the output \code{snpRdata} object, with missing data matching
#'   the missing data format of \code{x} added where genotypes are not in
#'   either \code{x} or \code{y}.
#' @param all.x.snps,all.y.snps logical, default \code{all}. Keep SNPs in the data
#'   even if they are only present in \code{x} or \code{y}, respectively.
#' @param all.x.samples,all.y.samples logical, default \code{all}. Keep samples in the
#'   data even if they are only present in \code{x} or \code{y}, respectively.
#' @param resolve_conflicts character, default 'error'. Controls how
#'   conflicting genotypic information in \code{x} and \code{y} is handled. See
#'   'Details' for options and explanation.
#' 
#' @author William Hemstrom
#' 
#' @return A merged \code{snpRdata} object.
#' 
#' @export
#' 
#' @examples
#' # create data to merge in
#' y <- data.frame(s1 = c("GG", "NN"),
#'                 s2 = c("GG", "TG"),
#'                 s3 = c("GG", "TT"),
#'                 s4 = c("GA", "TT"),
#'                 s5 = c("GG", "GT"),
#'                 s6 = c("NN", "GG"))
#'                 
#' snp.y <- data.frame(chr = c("groupVI", "test_chr"),
#'                     position = c(212436, 10))
#'                    
#' samp.y <- data.frame(pop = c("ASP", "ASP", "ASP", "test1", "test2", "test3"),
#'                      ID = c(1, 2, 3, "A1", "A2", "A3"),
#'                      fam = c("A", "B", "C", "T", "T", "T"))
#' y <- import.snpR.data(y, snp.y, samp.y)
#' 
#' x <- stickSNPs
#' sample.meta(x)$ID <- 1:ncol(x)
#' 
#' \dontrun{
#' # Not run, will error due to conflicts
#' z <- merge_snpRdata(x, y)
#' 
#' # Not run, will return a warning and report mismatches
#' z <- merge_snpRdata(x, y, resolve_conflicts = "warning")
#' }
#' 
#' # take a random genotype in the case of conflicts
#' z <- merge_snpRdata(x, y, resolve_conflicts = "random")
#' z
#' 
merge_snpRdata <- function(x, y, by.sample = intersect(names(sample.meta(x)), names(sample.meta(y))),
                           by.sample.x = by.sample, by.sample.y = by.sample,
                           by.snp = intersect(names(snp.meta(x)), names(snp.meta(y))),
                           by.snp.x = by.snp, by.snp.y = by.snp,
                           all = TRUE, all.x.snps = all, all.y.snps = all,
                           all.x.samples = all, all.y.samples = all,
                           resolve_conflicts = "error"){
  
  .sample.id.from.x <- .sample.id.from.y <- .snp.id.from.x <- .snp.id.from.y <- .source.ob <- NULL
  .sample.id <- .snp.id <- NULL
  #========sanity checks==============
  if(!is.snpRdata(x)){
    stop("x must be a snpRdata object.\n")
  }
  if(!is.snpRdata(y)){
    stop("y must be a snpRdata object.\n")
  }
  
  if(x@snp.form != y@snp.form){
    stop("x and y have snps in a different format.\n")
  }
  
  if(x@mDat != y@mDat){
    if(x@mDat %in% colnames(y@geno.tables$wm)){
      stop("x and y have missing data in a different format. The missing data format from x ", x@mDat, " is also present in the genotypes of y, so merging cannot procceed.\n")
    }
    
    warning("x and y have missing data in a different format. The missing data format from x ", x@mDat, " will be used.\n")
    gy <- genotypes(y)
    gy[gy == y@mDat] <- x@mDat
    y <- import.snpR.data(gy, snp.meta(y), sample.meta(y), x@mDat)
  }
  
  
  if(".sample.id" %in% by.sample){
    by.sample <- by.sample[-which(by.sample == ".sample.id")]
  }
  
  if(".snp.id" %in% by.snp){
    by.snp <- by.snp[-which(by.snp == ".snp.id")]
  }
  
  good.resolves <- c("warning", "error", "random", "x", "y")
  if(!resolve_conflicts %in% good.resolves){
    stop(paste0("'", resolve_conflicts, "' not a recognized option for the 'resolve_conflicts' argument.\n"))
  }
  
  #========match samples and SNPs using the usual merge tool, which also handles all of the 'all' and 'by' options naturally=============
  sm1 <- sample.meta(x)
  sm2 <- sample.meta(y)
  
  # handle the case where there are no matching sample metadata columns.
  if(length(by.sample.x) == 0 | length(by.sample.y) == 0){
    if(all.x.samples & all.y.samples){
      sm.m <- data.table::rbindlist(list(x = sm1, y = sm2), fill = TRUE, idcol = ".source.ob")
      sm.m[,.sample.id.from.x := ifelse(.source.ob == "x", .sample.id, NA)]
      sm.m[,.sample.id.from.y := ifelse(.source.ob == "y", .sample.id, NA)]
      sm.m$.source.ob <- NULL
      sm.m <- as.data.frame(sm.m)
    }
    else{
      stop("No matching sample metadata columns by which to merge.\n")
    }
  }
  else{
    sm.m <- merge(sm1, sm2, by.x = by.sample.x, by.y = by.sample.y,
                  all.x = all.x.samples, all.y = all.y.samples, suffixes = c(".from.x", ".from.y"), sort = FALSE)
  }
  
  sm.id.cols <- grep("^\\.sample\\.id", colnames(sm.m))
  
  
  snm1 <- snp.meta(x)
  snm2 <- snp.meta(y)
  
  # handle the case where there are no matching snp metadata columns.
  if(length(by.snp.x) == 0 | length(by.snp.y) == 0){
    if(all.x.snps & all.y.snps){
      snm.m <- data.table::rbindlist(list(x = snm1, y = snm2), fill = TRUE, idcol = ".source.ob")
      snm.m[,.snp.id.from.x := ifelse(.source.ob == "x", .snp.id, NA)]
      snm.m[,.snp.id.from.y := ifelse(.source.ob == "y", .snp.id, NA)]
      snm.m$.source.ob <- NULL
      snm.m <- as.data.frame(snm.m)
    }
    else{
      stop("No matching snp metadata columns by which to merge.\n")
    }
  }
  else{
    snm.m <- merge(snm1, snm2, by.x = by.snp.x, by.y = by.snp.y,
                   all.x = all.x.snps, all.y = all.y.snps, suffixes = c(".from.x", ".from.y"), sort = FALSE)
  }
  
  snm.id.cols <- grep("^\\.snp\\.id", colnames(snm.m))
  
  
  genotypes.m <- matrix(x@mDat, nrow = nrow(snm.m), ncol = nrow(sm.m))
  gt <- list(gs = matrix(nrow = nrow(snm.m), ncol = length(unique(c(colnames(x@geno.tables$gs), colnames(y@geno.tables$gs))))),
             as = matrix(nrow = nrow(snm.m), ncol = length(unique(c(colnames(x@geno.tables$as), colnames(y@geno.tables$as))))),
             wm = matrix(nrow = nrow(snm.m), ncol = 1))
  
  #=======merge genotypes by indices in merged metadata--genotypes present in both============
  if(nrow(snm.m) == 0){
    stop("No loci remain after merging.\n")
  }
  if(nrow(sm.m) == 0){
    stop("No samples remain after merging.\n")
  }
  idents.snp <- which(!is.na(snm.m$.snp.id.from.x) & !is.na(snm.m$.snp.id.from.y))
  idents.samp <- which(!is.na(sm.m$.sample.id.from.x) & !is.na(sm.m$.sample.id.from.y))
  
  if(length(idents.snp) > 0 & length(idents.samp) > 0){
    g.matches.x <- genotypes(x)[match(snm.m$.snp.id.from.x[idents.snp], snp.meta(x)$.snp.id), 
                                match(sm.m$.sample.id.from.x[idents.samp], sample.meta(x)$.sample.id)]
    g.matches.y <- genotypes(y)[match(snm.m$.snp.id.from.y[idents.snp], snp.meta(y)$.snp.id), 
                                match(sm.m$.sample.id.from.y[idents.samp], sample.meta(y)$.sample.id)]
    
    if(all(c(dim(g.matches.y), dim(g.matches.x))) > 0){
      
      # set matches to x or y if requested, filling missing genotypes when possible from other set
      if(resolve_conflicts == "y"){
        g.matches.m <- g.matches.y
        g.matches.m[g.matches.m == x@mDat] <- g.matches.x[g.matches.m == x@mDat]
      }
      else if(resolve_conflicts == "x"){
        g.matches.m <- g.matches.x
        g.matches.m[g.matches.m == x@mDat] <- g.matches.y[g.matches.m == x@mDat]
      }
      else{
        if(resolve_conflicts == "error"){
          stop("Some genotypes at identical loci sequenced in samples in both 'x' and 'y' are not identical. For more information on conflicts, merge_snpRdata() may instead be run with the 'resolve_conflicts' argument set to 'warning'\n")
        }
        
        # handle mismatches
        ## first sub-in any missing genotypes, since these will be seen as implicit matches in favor of the one without missing data.
        g.matches.m <- g.matches.x
        g.matches.m[g.matches.m == x@mDat] <- g.matches.y[g.matches.m == x@mDat]
        g.matches.y[g.matches.y == x@mDat] <- g.matches.m[g.matches.y == x@mDat]
        matching_idents <- g.matches.m == g.matches.y # should be no mDats here, since those will now match.
        
        
        # throw an error if requested
        if(resolve_conflicts == "warning"){
          if(any(!matching_idents)){
            warning("Error: Genotypic mismatches in identical samples/snps. Returning matrix of mismatches.")
            matching_idents <- cbind(snp.meta(x)[match(snm.m$.snp.id.from.x[idents.snp], snp.meta(x)$.snp.id),], matching_idents)
            colnames(matching_idents)[(ncol(snp.meta(x)) + 1): ncol(matching_idents)] <- paste0("Obj_x_sample_", match(sm.m$.sample.id.from.x[idents.samp], sample.meta(x)$.sample.id))
            return(matching_idents)
          }
        }
        
        # grab a random x or y genotype if requested
        else if(resolve_conflicts == "random"){
          take.y <- stats::rbinom(sum(!matching_idents), 1, .5)
          take.y <- as.logical(take.y)
          
          g.matches.m[!matching_idents][take.y] <- g.matches.y[!matching_idents][take.y]
        }
      }
      
      genotypes.m[idents.snp, idents.samp] <- as.matrix(g.matches.m)
    }
  }
  
  #=======merge genotypes by indices in merged metadata--genotypes not present in both============
  genotypes.m <- data.table::as.data.table(genotypes.m)
  
  # samples in only one:
  x.only <- which(is.na(sm.m$.sample.id.from.y) & !is.na(sm.m$.sample.id.from.x))
  y.only <- which(is.na(sm.m$.sample.id.from.x) & !is.na(sm.m$.sample.id.from.y))
  
  if(length(x.only) > 0){
    set(genotypes.m, which(!is.na(snm.m$.snp.id.from.x)), x.only, 
        data.table::as.data.table(genotypes(x)[match(snm.m$.snp.id.from.x[which(!is.na(snm.m$.snp.id.from.x))], snp.meta(x)$.snp.id),
                                               match(sm.m$.sample.id.from.x[x.only], sample.meta(x)$.sample.id)]))
  }
  if(length(y.only) > 0){
    set(genotypes.m, which(!is.na(snm.m$.snp.id.from.y)), y.only,
        data.table::as.data.table(genotypes(y)[match(snm.m$.snp.id.from.y[which(!is.na(snm.m$.snp.id.from.y))], snp.meta(y)$.snp.id),
                                               match(sm.m$.sample.id.from.y[y.only], sample.meta(y)$.sample.id)]))
  }
  
  # snps in only one
  x.only <- which(is.na(snm.m$.snp.id.from.y) & !is.na(snm.m$.snp.id.from.x))
  y.only <- which(is.na(snm.m$.snp.id.from.x) & !is.na(snm.m$.snp.id.from.y))
  
  if(length(x.only) > 0){
    set(genotypes.m, x.only, which(!is.na(sm.m$.sample.id.from.x)),
        data.table::as.data.table(genotypes(x)[match(snm.m$.snp.id.from.x[x.only], snp.meta(x)$.snp.id),
                                               match(sm.m$.sample.id.from.x[which(!is.na(sm.m$.sample.id.from.x))], sample.meta(x)$.sample.id)]))
  }
  if(length(y.only) > 0){
    set(genotypes.m, y.only, which(!is.na(sm.m$.sample.id.from.y)),
        data.table::as.data.table(genotypes(y)[match(snm.m$.snp.id.from.y[y.only], snp.meta(y)$.snp.id),
                                               match(sm.m$.sample.id.from.y[which(!is.na(sm.m$.sample.id.from.y))], sample.meta(y)$.sample.id)]))
    
    
  }
  
  # fix any lingering "."'s
  bad.snm.cols <- grep("\\.", colnames(snm.m))
  if(length(bad.snm.cols) > 0){
    colnames(snm.m)[bad.snm.cols] <- gsub("\\.", "_", colnames(snm.m)[bad.snm.cols])
  }
  bad.sm.cols <- grep("\\.", colnames(sm.m))
  if(length(bad.sm.cols) > 0){
    colnames(sm.m)[bad.sm.cols] <- gsub("\\.", "_", colnames(sm.m)[bad.sm.cols])
  }
  
  return(import.snpR.data(genotypes = genotypes.m, 
                          snp.meta = snm.m[,-grep("_snp_id", colnames(snm.m))],
                          sample.meta = sm.m[,-grep("_sample_id", colnames(sm.m))],
                          mDat = x@mDat))
}
