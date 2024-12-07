#' Check if an object is a snpRdata object.
#'
#' Checks if an object is a snpRdata object.
#'
#' @param x object to check
#' @return Logical, TRUE if the object is a snpRdata object.
#'
#' @author William Hemstrom
#'
#' @export
is.snpRdata <- function(x){
  return(methods::is(x, "snpRdata"))
}

# Add facets to snpRdata objects
# 
# Adds facets to snpRdata objects, following the rules described in
# \code{\link{Facets_in_snpR}}. Internal, called by many other functions when
# facets are requested.
# 
# @param x snpRdata object
# @param facets character. Facets to use.
# 
.add.facets.snpR.data <- function(x, facets = NULL){
  #===========================binding issues for unquoted variables============
  .snp.id <- facet <- subfacet <- NULL
  
  
  if(is.null(facets[1])){return(x)}
  facets <- .check.snpR.facet.request(x, facets)
  if(is.null(facets[1])){return(x)}
  
  #===========================smart rbind fill function========================
  .smart.rbind.fill.matrix <- function(x, y, nrf = 1){

    if(!methods::is(x, "sparseMatrix")){
      x <- Matrix::Matrix(x, sparse = TRUE)
    }
    
    if(!methods::is(y, "sparseMatrix")){
      y <- Matrix::Matrix(y, sparse = TRUE)
    }
    
    if(nrow(y) != 0 & nrow(x) != 0){
      return(.rbind_list_fill_sparse(list(x, y)))
    }
    
    if(nrow(y) == 0 & nrow(x) == 0){
      return(Matrix::Matrix(0, nrow = 0, ncol = 0))
    }
    
    if(nrow(y) == 0){
      fill <- Matrix::Matrix(0, nrf, ncol(x))
      colnames(fill) <- colnames(x)
      return(.rbind_list_fill_sparse(list(x, fill)))
    }
    
    if(nrow(x) == 0){
      fill <- Matrix::Matrix(0, nrf, ncol(y))
      colnames(fill) <- colnames(y)
      return(.rbind_list_fill_sparse(list(fill, y)))
    }
  }
  
  #===========================turn into list========
  # need to fix any multivariate facets (those with a .)
  comp.facets <- grep("(?<!^)\\.", facets, perl = T)
  if(length(comp.facets) != 0){
    run.facets <- as.list(facets[-c(comp.facets)])
    facet.list <- c(run.facets, .split.facet(facets[comp.facets]))
  }
  else{
    facet.list <- as.list(facets)
  }
  
  #===========================sanity checks=========
  # check that we haven't done these facets before, and remove any that we have.
  all.facets <- character(length = length(facet.list))
  for(i in 1:length(facet.list)){
    all.facets[i] <- paste0(facet.list[[i]], collapse = ".")
  }
  if(all(all.facets %in% x@facets)){return(x)}
  else if (any(all.facets %in% x@facets)){
    already.added <- all.facets[which(all.facets %in% x@facets)]
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
  # for sample acets, prep summary data and geno tables.
  # For snp facets, nothing to do here besides add them to the facet list and facet type list.
  # For complex facets, prep summary tables for the sample facet portion if they don't exist and add to facet list.
  added.facets <- character(0)
  # oac <- cbind(data.table::as.data.table(x@facet.meta[,c("facet", "subfacet", ".snp.id")]),
  #              data.table::as.data.table(x@ac)) # grab original ac with meta for later
  for(k in 1:length(facet.list)){
    facets <- facet.list[[k]] # column levels for this facet.
    #=========================figure out unique levels for the facet==========
    # figure out what kind of facets we are working with.
    
    # save info and get the unique options for each facet
    # save info
    if(length(facets) > 1){
      x@facets <- c(x@facets, paste0(facets, collapse = "."))
      added.facets <- c(added.facets, paste0(facets, collapse = "."))
    }
    else{
      x@facets <- c(x@facets, facets)
      added.facets <- c(added.facets, facets)
    }
    x@facet.type <- c(x@facet.type, "sample")
    
    
    # get unique options for this facet
    sample.meta <- x@sample.meta[colnames(x@sample.meta) %in% facets]
    sample.meta <- sample.meta[,sort(colnames(sample.meta))]
    if(!is.data.frame(sample.meta)){
      sample.meta <- as.data.frame(sample.meta, stringsAsFactors = F)
      colnames(sample.meta) <- colnames(x@sample.meta)[colnames(x@sample.meta) %in% facets]
    }
    sample.meta <- dplyr::mutate_all(sample.meta, as.character) # this fixes some really obscure bugs with integers in columns.
    sample.opts <- unique(sample.meta)
    if(!is.data.frame(sample.opts)){
      sample.opts <- as.data.frame(sample.opts, stringsAsFactors = F)
      colnames(sample.opts) <- facets[which(facets %in% colnames(x@sample.meta))]
    }
    sample.opts <- dplyr::arrange_all(sample.opts)
    
    
    gs <- x@geno.tables
    #=========================get gs matrices==========
    for(i in 1:nrow(sample.opts)){
      matches <- which(apply(sample.meta, 1, function(x) identical(as.character(x), as.character(sample.opts[i,]))))
      t.x <- genotypes(x)[,matches, drop = FALSE]
      tgs <- .tabulate_genotypes(t.x, x@mDat)
      
      gs$gs <- .smart.rbind.fill.matrix(gs$gs, tgs$gs, nrow(t.x))
      gs$as <- .smart.rbind.fill.matrix(gs$as, tgs$as, nrow(t.x))
      gs$wm <- .smart.rbind.fill.matrix(gs$wm, tgs$wm, nrow(t.x))
      # fix NAs that show up when there are less called genotype options in one facet level than in all levels!
      if(ncol(gs$gs) != ncol(tgs$gs)){
        gs <- lapply(gs, function(x){x[is.na(x)] <- 0;x})
      }
      
      x@facet.meta <- rbind(x@facet.meta,
                            cbind(data.frame(facet = rep(paste0(facets, collapse = "."), nrow(t.x)),
                                             subfacet = rep(paste0(sample.opts[i,], collapse = "."), nrow(t.x)),
                                             facet.type = rep("sample", nrow(t.x)), stringsAsFactors = F),
                                  x@snp.meta))
    }
    

    #=========================sort, pack, and return==========
    # sort
    x@facet.meta <- dplyr::mutate_if(.tbl = x@facet.meta, is.factor, as.character)
    x@facet.meta$.reorder <- 1:nrow(x@facet.meta)
    x@facet.meta <- dplyr::arrange(x@facet.meta, .snp.id, facet, subfacet)
    gs$gs <- gs$gs[x@facet.meta$.reorder,, drop = F]
    gs$as <- gs$as[x@facet.meta$.reorder,, drop = F]
    gs$wm <- gs$wm[x@facet.meta$.reorder,, drop = F]
    x@facet.meta <- x@facet.meta[,-ncol(x@facet.meta)]
    
    # output
    x@geno.tables <- gs
  }
  # add and sort ac formated data.
  # if(.is.bi_allelic(x)){
  #   .make_it_quiet(nac <- format_snps(x, output = "ac", facets = added.facets))
  #   nac <- data.table::as.data.table(nac)
  #   nac <- rbind(oac, nac[,c("facet", "subfacet", ".snp.id", "n_total","n_alleles", "ni1", "ni2")])
  #   nac <- dplyr::mutate_if(.tbl = nac, is.factor, as.character)
  #   nac <- dplyr::arrange(nac, .snp.id, facet, subfacet)
  #   nac <- as.data.frame(nac)
  #   x@ac <- nac[,-c(1:3)]
  # }
  

  # add in dummy rows to stats
  sm <- data.table::as.data.table(x@facet.meta[x@facet.meta$facet %in% added.facets, c("facet", "subfacet", "facet.type", colnames(x@snp.meta))])
  if(ncol(x@stats) - ncol(sm) > 0){
    sm <- cbind(sm, matrix(NA, nrow(sm), ncol(x@stats) - ncol(sm)))
  }
  colnames(sm) <- colnames(x@stats)
  os <- data.table::as.data.table(x@stats)
  if(ncol(os) - ncol(sm) > 0){
    os <- rbind(os, cbind(sm, matrix(NA, nrow(sm), ncol(os) - ncol(sm))))
  }
  else{
    os <- rbind(os, sm)
  }
  
  os <- dplyr::mutate_if(.tbl = os, is.factor, as.character)
  x@stats <- as.data.table(dplyr::arrange(os, .snp.id, facet, subfacet))
  
  
  return(x)
}

# Remove tabulated facets. Useful mostly for subsetting where a sample facet
# needs to be cleaned up and removed.
#
# @param x snpRdata object
# @param facets sample-level facets to remove. Will remove the full tabulation,
#   not by parts. For example, "pop.fam" will remove JUST pop-fam, not "pop" or
#   "fam".
.remove.facets.snpR.data <- function(x, facets = NULL){
  #==========check removal request=============
  if(is.null(facets)){
    return(x)
  }
  if(any(facets == ".base")){
    facets <- facets[-which(facets == ".base")]
    if(length(facets) == 0){
      return(x)
    }
  }

  # anything to remove?
  rm.facets <- x@facets[which(x@facets %in% facets)]
  if(length(rm.facets) == 0){
    return(x)
  }
  
  #==========remove every trace of that facet============
  # otherwise continue and find the rows to remove
  rm.rows <- which(x@facet.meta$facet %in% rm.facets)
  
  # remove everything
  x@geno.tables <- lapply(x@geno.tables, function(x) x[-rm.rows,])
  x@stats <- x@stats[-rm.rows,]
  
  # pairwise
  rm.pairwise.rows <- which(x@pairwise.stats$facet %in% rm.facets)
  if(length(rm.pairwise.rows) > 0){
    x@pairwise.stats <- x@pairwise.stats[-rm.pairwise.rows,]
    if(nrow(x@pairwise.window.stats) == 0){
      x@pairwise.window.stats <- data.frame()
    }
  }
  
  # window pairwise
  rm.window.pairwise.rows <- which(x@pairwise.window.stats$facet %in% rm.facets)
  if(length(rm.window.pairwise.rows) > 0){
    x@pairwise.window.stats <- x@pairwise.window.stats[-rm.window.pairwise.rows,]
    if(nrow(x@pairwise.window.stats) == 0){
      x@pairwise.window.stats <- data.frame()
    }
  }
  
  # standard window
  rm.window.rows <- which(x@window.stats$facet %in% rm.facets)
  if(length(rm.window.rows) > 0){
    x@window.stats <- x@window.stats[-rm.window.rows,]
    if(nrow(x@window.stats) == 0){
      x@window.stats <- data.frame()
    }
  }
  
  # window bootstraps
  if(nrow(x@window.bootstraps) > 0){
    rm.window.bootstraps.rows <- which(x@window.bootstraps$facet %in% rm.facets)
    if(length(rm.window.bootstraps.rows) > 0){
      x@window.bootstraps <- x@window.bootstraps[-rm.window.bootstraps.rows,]
      if(nrow(x@window.bootstraps) == 0){
        x@window.bootstraps <- data.frame()
      }
    }
  }
  
  # sample stats -- no need to remove anything (but may need to edit metadata if this is called due to facet removal!)
  
  # pop stats
  if(nrow(x@pop.stats) > 0){
    rm.pop.stat.rows <- which(x@pop.stats$facet %in% rm.facets)
    if(length(rm.pop.stat.rows) > 0){
      x@pop.stats <- x@pop.stats[-rm.pop.stat.rows,]
      if(nrow(x@pop.stats) == 0){
        x@pop.stats <- data.frame()
      }
    }
  }
  
  # LD
  if(length(x@pairwise.LD) != 0){
    prox.rm <- which(x@pairwise.LD$prox$sample.facet %in% rm.facets)
    if(length(prox.rm) > 0){
      x@pairwise.LD$prox <- x@pairwise.LD$prox[-prox.rm,]
      
      # if nothing remains, done
      if(nrow(x@pairwise.LD$prox) == 0){
        x@pairwise.LD <- list()
      }
      
      # otherwise check matrices and remove any with removed components
      else{
        bad.matrices <- .split.facet(names(x@pairwise.LD$LD_matrices))
        bad.matrices <- lapply(bad.matrices, function(y) .check.snpR.facet.request(x, y))
        bad.matrices <- lapply(bad.matrices, function(y){
          return(unlist(lapply(rm.facets, function(z){z == y})))
        })
        bad.matrices <- which(unlist(lapply(bad.matrices, function(y) ifelse(sum(y) == 0, FALSE, TRUE))))
        
        x@pairwise.LD$LD_matrices[bad.matrices] <- NULL
      }
    }
  }
  
  # calced stats
  bad.stats <- .split.facet(names(x@calced_stats))
  bad.stats <- lapply(bad.stats, function(y) .check.snpR.facet.request(x, y))
  bad.stats <- lapply(bad.stats, function(y){
    return(unlist(lapply(rm.facets, function(z){z == y})))
  })
  bad.stats <- which(unlist(lapply(bad.stats, function(y) ifelse(sum(y) == 0, FALSE, TRUE))))
  x@calced_stats[bad.stats] <- NULL
  
  
  # af matrices
  if(length(x@allele_frequency_matrices) > 0){
    bad.matrices <- .split.facet(names(x@allele_frequency_matrices))
    bad.matrices <- lapply(bad.matrices, function(y) .check.snpR.facet.request(x, y))
    bad.matrices <- lapply(bad.matrices, function(y){
      return(unlist(lapply(rm.facets, function(z){z == y})))
    })
    bad.matrices <- which(unlist(lapply(bad.matrices, function(y) ifelse(sum(y) == 0, FALSE, TRUE))))
    
    x@allele_frequency_matrices[bad.matrices] <- NULL
  }
  
  # genetic distances
  if(length(x@genetic_distances) > 0){
    bad.matrices <- .split.facet(names(x@genetic_distances))
    bad.matrices <- lapply(bad.matrices, function(y) .check.snpR.facet.request(x, y))
    bad.matrices <- lapply(bad.matrices, function(y){
      return(unlist(lapply(rm.facets, function(z){z == y})))
    })
    bad.matrices <- which(unlist(lapply(bad.matrices, function(y) ifelse(sum(y) == 0, FALSE, TRUE))))
    
    x@genetic_distances[bad.matrices] <- NULL
  }
  
  # weighted means
  if(nrow(x@weighted.means) > 0){
    rm.mean.rows <- which(x@weighted.means$facet %in% rm.facets)
    if(length(rm.mean.rows) > 0){
      x@weighted.means <- x@weighted.means[-rm.mean.rows,]
      if(nrow(x@weighted.means) == 0){
        x@weighted.means <- data.frame()
      }
    }
  }
  
  # IBD
  if("ibd" %in% names(x@other)){
    bad.ibd <- .split.facet(names(x@other$ibd))
    bad.ibd <- lapply(bad.ibd, function(y) .check.snpR.facet.request(x, y))
    bad.ibd <- lapply(bad.ibd, function(y){
      return(unlist(lapply(rm.facets, function(z){z == y})))
    })
    bad.ibd <- which(unlist(lapply(bad.ibd, function(y) ifelse(sum(y) == 0, FALSE, TRUE))))
    
    x@other$ibd[bad.ibd] <- NULL
    x@other$geo_dists[bad.ibd] <- NULL
    if(length(x@other$ibd) == 0){
      x@other$ibd <- NULL
      x@other$geo_dists <- NULL
    }
  }
  
  # facet meta
  x@facet.meta <- x@facet.meta[-rm.rows,]
  
  # facets and facet type
  rmfn <- which(x@facets %in% facets)
  x@facets <- x@facets[-rmfn]
  x@facet.type <- x@facet.type[-rmfn]

  #============return=========
  return(x)
}


# Apply functions across snpR facets.
#
# Internal function to apply functions across snpR data. Facet designations
# follow rules described in \code{\link{Facets_in_snpR}}. Null or "all" facet
# designations follow typical rules.
#
# This function should never be called externally. Raw statistics with metadata
# are returned and are not automatically merged into x.
#
# For examples, look at how this function is called in functions such as
# calc_pi, calc_pairwise_fst, etc.
#
# Options:
#
# req:
#
# \itemize{ \item{gs: } genotype tables \item{ac: } ac formatted data
# \item{meta.gs: } facet, .snp.id metadata cbound genotype tables.
# \item{ac.stats: } ac data cbound to stats \item{meta.as: } as data cbound to
# snp metadata. \item{snpRdata: } subset snpRdata object. }
#
# case:
#
# \itemize{ \item{ps: } stat calculated per snp. \item{ps.pf: } stat calculated
# per snp, but split per facet (such as for private alleles, comparisons only
# exist within a facet!) \item{facet.pairwise: } stat calculated pairwise
# between facets, per snp otherwise. \item{ps.pf.psf: } stat calculated per snp,
# but per both sample and snp facet. \item{psamp:} calculated individually in
# each sample, may be different across snp facets.}
#
# @param x snpRdata object
# @param facets character or NULL, default NULL. Facets to add.
# @param req character. Data type required by fun. See description.
# @param fun function. Function to apply to data.
# @param case character, default "ps". Type of statistic required by fun. See
#   description.
# @param par numeric or FALSE, default FALSE. Number of parallel computing cores
#   to use. Works for some cases/reqs.
# @param ... other arguments to be passed to fun.
# @param stats.type character, default "all". Other options "pairwise" or
#   "stats". Type of statistic under to pull under the ps.pf.psf option.
# @param response character, default NULL. Possible response variable name
#   passed to specific functions.
# @param maf numeric, defualt NULL. Possible maf, passed to some functions.
# @param interpolate character, default NULL. Possible interpolation option,
#   passed to some functions.
# @param verbose If TRUE, will cat updates to console.
#
# @author William Hemstrom
.apply.snpR.facets <- function(x, facets = NULL, req, fun, case = "ps", par = FALSE, ..., stats.type = "all", response = NULL, maf = FALSE, interpolate = NULL, verbose = FALSE){

  if(!is.null(facets)){
    if(facets[1] == "all"){
      facets <- x@facets
    }
  }
  else {
    facets <- ".base"
  }
  
  
  if(case == "ps"){
    facets <- .check.snpR.facet.request(x, facets)
    if(req == "gs"){
      # subset
      gs <- lapply(x@geno.tables, function(y) y[which(x@facet.meta$facet %in% facets),,drop = FALSE])
      
      # run the function indicated
      out <- data.table::as.data.table(fun(gs, ...))
      
      # bind metadata
      out <- cbind(data.table::as.data.table(x@facet.meta[x@facet.meta$facet %in% facets,]), out)
      #out <- cbind(x@facet.meta[x@facet.meta$facet %in% facets,], out)
      
      # return
      return(out)
    }
    else if(req == "meta.gs"){
      # gs with metadata on facets attached.
      # subset
      gs <- lapply(x@geno.tables, function(y) y[x@facet.meta$facet %in% facets,,drop=F])
      gs <- lapply(gs, function(y) cbind(x@facet.meta[x@facet.meta$facet %in% facets, c("facet", "subfacet", ".snp.id")], as.matrix(y)))
      
      # run the function indicated
      out <- data.table::as.data.table(fun(gs, ...))
      
      # bind metadata
      out <- cbind(data.table::as.data.table(x@facet.meta[x@facet.meta$facet %in% facets,]), out)
      
      # return
      return(out)
    }
    else if(req == "ac"){
      # subset
      ac <- x@ac[x@facet.meta$facet %in% facets,]
      
      # run the function indicated
      out <- data.table::as.data.table(fun(ac, ...))
      
      # bind metadata
      out <- cbind(data.table::as.data.table(x@facet.meta[x@facet.meta$facet %in% facets,]), out)
      
      # return
      return(out)
    }
    else if(req == "cast.gs" | req == "cast.ac"){

      # need to split by phenotype + facet
      pheno.facets <- paste(facets, response, sep = ".")
      
      ## remove any .base pheno facets, do these a bit differently
      base.pheno.facets <- grep("^\\.base", pheno.facets)
      if(length(base.pheno.facets) > 0){
        if(length(base.pheno.facets) != length(pheno.facets)){
          pheno.facets <- c(.check.snpR.facet.request(x, pheno.facets[-base.pheno.facets]), response)
        }
        else{
          pheno.facets <- response
        }
      }
      else{
        pheno.facets <- .check.snpR.facet.request(x, pheno.facets)
      }
      x <- .add.facets.snpR.data(x, pheno.facets)
      if(req == "cast.gs"){
        gs <- data.table::as.data.table(cbind(x@facet.meta, as.matrix(x@geno.tables$gs)))
      }
      else{
        gs <- format_snps(x, "ac", pheno.facets)
        gs <- gs[,-which(colnames(gs) %in% c("n_total", "n_alleles"))]
      }
      
      # initialize outputs
      comb <- vector("list", length = length(facets))

      # get the input gs data. Split by facet, since we need to extract metadata differently for each.
      for(i in 1:length(facets)){

        tgs <- gs[gs$facet == pheno.facets[i],]
        
        which.is.phenotype <- which(unlist(.split.facet(pheno.facets[i])) == response)
        
        # add phenotype and subfacet column
        id <- .split.facet(tgs$subfacet)
        l <- length(id[[1]])
        uid <- unlist(id)
        p.vals <- seq(from = which.is.phenotype, to = length(uid), by = l)
        subfacet <- lapply(id, FUN = function(x) x[-which.is.phenotype])
        subfacet <- unlist(lapply(subfacet, FUN = paste, collapse = "."))
        tgs$subfacet <- subfacet
        tgs$phenotype <- uid[p.vals]
        tgs$facet.type <- .check.snpR.facet.request(x, facets[i], "none", T)[[2]]
        
        # fix for .base
        tgs$subfacet[tgs$subfacet == ""] <- ".base"
        
        # cast
        if(req == "cast.gs"){
          value.vars <- colnames(x@geno.tables$gs)
        }
        else{
          value.vars <- c("ni1", "ni2")
        }
        comb[[i]] <- data.table::dcast(data.table::setDT(tgs), .snp.id + subfacet ~ phenotype, value.var = value.vars, fun.aggregate = sum)
        comb[[i]] <- merge(tgs[tgs$phenotype == unique(tgs$phenotype)[1],1:which(colnames(tgs) == ".snp.id")], comb[[i]], by = c(".snp.id", "subfacet"), all.y = T)
        
        # levels where there is no case/no control
        if(any(is.na(comb[[i]]$facet))){
          ..target_cols <- NULL
          fill <- merge(comb[[i]][is.na(comb[[i]]$facet),c(".snp.id", "subfacet")],
                        tgs[tgs$phenotype == unique(tgs$phenotype),1:which(colnames(tgs) == ".snp.id")],
                        by = c(".snp.id", "subfacet"), all.x = TRUE)
          target_cols <- c("subfacet", "facet", "facet.type", 
                           colnames(comb[[i]])[which(colnames(comb[[i]]) %in% colnames(snp.meta(x)) & colnames(comb[[i]]) != ".snp.id")])
          data.table::set(comb[[i]], which(is.na(comb[[i]]$facet)), target_cols,
                          .fix..call(fill[,..target_cols]))
        }
        
        comb[[i]]$facet <- facets[i]
      }
      
      comb <- data.table::rbindlist(comb)
      mcols <- 1:ncol(x@facet.meta)
      meta <- comb[,mcols, with = FALSE]
      comb <- comb[,-mcols, with = FALSE]
      
      # call function
      out <- fun(comb, ...)
      n.ord <- c("facet", "subfacet", "facet.type", colnames(snp.meta(x))[-which(colnames(snp.meta(x)) == ".snp.id")], ".snp.id")
      if(any(!n.ord %in% colnames(meta))){
        n.ord <- n.ord[-which(!n.ord %in% colnames(meta))]
      }
      meta <- meta[,n.ord, with = FALSE]
      out <- cbind(meta, out)
      if(req == "cast.gs"){
        colnames(out)[ncol(out)] <- paste0("p_armitage_", response)
      }
      else{
        colnames(out)[-c(1:ncol(meta))] <- paste0(colnames(out)[-c(1:ncol(meta))], "_", response)
      }
      
      return(out)
    }
    else if(req == "snpRdata"){
      #==========subfunctions========
      run.one.task <- function(opts, i){
        if(opts[i,"t.snp.facet"] == ".base" & opts[i,"t.sample.facet"] == ".base"){
          sub.x <- x
        }
        else if(opts[i,"t.snp.facet"] == ".base"){
          suppressWarnings(.make_it_quiet(sub.x <- .subset_snpR_data(x, facets = opts[i,1], subfacets = opts[i,2])))
        }
        else{
          suppressWarnings(.make_it_quiet(sub.x <- .subset_snpR_data(x, facets = opts[i,1], subfacets = opts[i,2],
                                                                              snp.facets = opts[i,3], snp.subfacets = opts[i,4])))
        }
        
        suppressWarnings(.make_it_quiet(sub.x <- filter_snps(sub.x, non_poly = FALSE, maf = maf)))
        

        out <- fun(sub.x = sub.x, ...)
        
        # return
        if(!is.data.frame(out)){
          out$.fm <- cbind(facet = opts[i,1], subfacet = opts[i,2], row.names = NULL)
          return(out)
        }
        else{
          out <- cbind(facet = opts[i,1], subfacet = opts[i,2], out, row.names = NULL)
          return(out)
        }
      }
      #==========run===============
      # check facets
      facets <- .check.snpR.facet.request(x, facets)
      
      # get options
      opts.list <- .get.task.list(x, facets)
      
      # run in serial
      if(par == FALSE | nrow(opts.list) == 1){
        
        # initialize out
        out <- vector("list", nrow(opts.list))
        
        #run loop
        for(i in 1:nrow(opts.list)){
          out[[i]] <- run.one.task(opts.list, i)
        }
      }
      # run in parallel
      else{
        cl <- parallel::makePSOCKcluster(par)
        doParallel::registerDoParallel(cl)
        
        #prepare reporting function
        ntasks <- nrow(opts.list)
        # if(verbose){
        #   progress <- function(n) cat(sprintf("Part %d out of", n), ntasks, "is complete.\n")
        #   opts <- list(progress=progress)
        # }
        # else{
        #   opts <- list()
        # }
        
        
        if(verbose){cat("Begining run.\n")}
        
        # run the LD calculations
        ## suppress warnings because you'll get wierd ... warnings. Not an issue in the non-parallel version.
        suppressWarnings(out <- foreach::foreach(q = 1:ntasks, .inorder = FALSE,
                                                 .export = c("data.table"), .packages = "snpR") %dopar% {
                                                   run.one.task(opts.list, q)
                                                 })
        
        #release cores
        parallel::stopCluster(cl)
      }
      
      # combine
      ## for rf
      if(!is.data.frame(out[[1]])){
        # return
        return(out)
      }
      ## for GMMAT
      else{
        suppressWarnings(out <- dplyr::bind_rows(out))
        id.col <- (which(colnames(out) == ".snp.id"))
        colnames(out)[(id.col+1):ncol(out)] <- tolower(colnames(out)[(id.col+1):ncol(out)])
        colnames(out)[(id.col+1):ncol(out)] <- paste0("gmmat_", colnames(out)[(id.col+1):ncol(out)], "_", response)
        n.ord <- which(colnames(out) %in% c("major", "minor"))
        n.ord <- c((1:ncol(out))[-n.ord], n.ord)
        out <- out[,n.ord]
        
        # return
        return(out)
      }
      
      
    }
  }
  else if(case == "ps.pf"){
    out <- data.frame() # initialize
    if(req == "meta.gs"){
      # gs with metadata on facets attached.
      # loop through each facet
      for(i in 1:length(facets)){
        # subset
        gs <- lapply(x@geno.tables, function(y) y[x@facet.meta$facet == facets[i],,drop = FALSE])
        gs <- lapply(gs, function(y) cbind(x@facet.meta[x@facet.meta$facet == facets[i], c("facet", "subfacet", ".snp.id")], as.matrix(y)))
        
        # run the function indicated
        this.out <- data.table::as.data.table(fun(gs, ...))
        
        # bind metadata
        this.out <- data.table::as.data.table(cbind(x@facet.meta[x@facet.meta$facet == facets[i],], this.out))
        
        # bind to full data
        out <- rbind(out, this.out)
      }
      # return
      return(out)
    }
  }
  else if(case == "facet.pairwise"){
    if(req == "snpRdata"){
      # initialize
      out1 <- data.frame()
      out2 <- data.frame()
      
      
      # in this case, the whole snpRdata object is requested, but should be passed several times with the correct facet noted.
      for(i in 1:length(facets)){
        out <- fun(x, facets[i], ...)
        # if this is a genepop fst output, we should expect a list of size two. Need to return the correct parts!
        if(length(out) == 2){
          suppressWarnings(out1 <- rbind(out1, cbind(data.table::as.data.table(x@snp.meta), facet = facets[i], data.table::as.data.table(out[[1]]))))
          out2 <- rbind(out2, data.frame(facet = facets[i], comparison = names(out[[2]]), weighted_mean_fst = out[[2]], stringsAsFactors = F))
        }
      }
      
      # return
      if(nrow(out2) > 0){
        return(list(out1, out2))
      }
      else{
        return(out1)
      }
    }
    
    
    else if(req == "ac.stats"){
      out <- data.frame()
      # need to loop through each facet, since they are all internally compared!
      for(i in 1:length(facets)){
        # subset
        ac <- x@ac[x@facet.meta$facet == facets[i],]
        stats <- x@stats[x@facet.meta$facet == facets[i],]
        
        # run the function indicated
        this.out <- data.table::as.data.table(fun(cbind(stats, ac), ...))
        
        # bind metadata and combine with full output
        suppressWarnings(out <- rbind(out, cbind(data.table::as.data.table(x@snp.meta), facet = facets[i], this.out)))
        
      }
      
      # return
      return(out)
      
    }
    
  }
  else if(case == "ps.pf.psf"){
    if(req == "pos.all.stats" | req == "meta.as"){
      x@facet.meta$facet <- as.character(x@facet.meta$facet)
      x@facet.meta$subfacet <- as.character(x@facet.meta$subfacet)
      
      # subfunctions:
      
      ## a function to run 'func' on one iteration/row of the task.list. q is the iteration. Par is needed only for progress printing.
      run.one.loop <- function(stats_to_use, meta, task.list, q, par){
        if(methods::is(stats_to_use, "sparseMatrix")){stats_to_use <- as.matrix(stats_to_use)}
        
        if(par == FALSE & verbose){cat("Sample Subfacet:\t", as.character(task.list[q,2]), "\tSNP Subfacet:\t", as.character(task.list[q,4]), "\n")}
        
        # get comps and run
        ## figure out which data rows contain matching sample facets
        sample.matches <- which(apply(meta[,1:2], 1, function(x) identical(as.character(x), as.character(task.list[q,1:2]))))
        
        snp.facets <- unlist(.split.facet(task.list[q,3]))

        if(snp.facets[1] != ".base"){
          # figure out which data rows contain matching snp facets
          snp.cols <- .paste.by.facet(meta, snp.facets, sep = " ")
          snp.matches <- which(snp.cols == task.list[q,4])
          
          # get the intersect and return.
          run.lines <- intersect(snp.matches, sample.matches)
          snp.res <- gsub(" ", ".", task.list[q,4])
        }
        else{
          run.lines <- sample.matches
          snp.res <- ".base"
        }
        
        # run the function and create a snp res metadata df to bind to the results.
        # assign("last.warning", NULL, envir = baseenv())
        res <- fun(cbind(meta[run.lines,,drop = FALSE], stats_to_use[run.lines,,drop = FALSE]), ...)
        
        # return
        if(nrow(res) == 1){
          return(cbind.data.frame(as.data.frame(t(task.list[rep(q, nrow(res)),1:2])),
                                  snp.facet = task.list[q,3],
                                  snp.subfacet = snp.res,
                                  res))
        }
        else{
          return(cbind.data.frame(task.list[rep(q, nrow(res)),1:2],
                                  snp.facet = rep(task.list[q,3], nrow(res)),
                                  snp.subfacet = rep(snp.res, nrow(res)),
                                  res))
        }
      }
      
      
      # get the statistics needed by the function and the task list. For now, expects only meta cbound to as or snp-wise statistics.
      if(req == "meta.as"){
        stats_to_use <- x@geno.tables$as
        task.list <- .get.task.list(x, facets)
        meta.to.use <- x@facet.meta
      }
      else{
        if(stats.type == "stats"){
          
          # task list for non-pairwise case
          ## get the columns containing statistics
          pos.col <- which(colnames(x@stats) == "position")
          start.col <- which(colnames(x@stats) == ".snp.id") + 1
          if(start.col > ncol(x@stats)){
            stop("No per-snp stats have been calculated for this dataset.\n")
          }
          cols.to.use <- start.col:ncol(x@stats)
          ## add nk and remove any non-numeric columns
          nk <- Matrix::rowSums(x@geno.tables$as)
          if(data.table::is.data.table(x@stats)){
            stats_to_use <- cbind(nk = nk, x@stats[,cols.to.use, with = FALSE])
          }
          else{
            stats_to_use <- cbind(nk = nk, x@stats[,cols.to.use])
          }
          numeric.cols <- which(sapply(stats_to_use, class) %in% c("numeric", "integer"))
          ## save
          stats_to_use <- stats_to_use[,numeric.cols, with = FALSE]
          task.list <- .get.task.list(x, facets)
          meta.to.use <- x@facet.meta
          
        }
        else{
          # grab the correct data columns from the pairwise stats dataset, and reorder so nk is first
          pos.col <- which(colnames(x@pairwise.stats) == "position")
          start.col <- which(colnames(x@pairwise.stats) == ".snp.id") + 1
          cols.to.use <- start.col:ncol(x@pairwise.stats)
          stats_to_use <- x@pairwise.stats[,cols.to.use, with = FALSE]
          nk.col <- which(colnames(stats_to_use) == "nk")
          n.col.ord <- c(nk.col, (1:ncol(stats_to_use))[-nk.col])
          stats_to_use <- stats_to_use[,n.col.ord, with = FALSE]
          task.list <- .get.task.list(x, facets, source = "pairwise.stats")
          
          # re-order the meta and save
          meta.to.use <- as.data.frame(x@pairwise.stats[,1:which(colnames(x@pairwise.stats) == ".snp.id")])
          facet.cols <- which(colnames(meta.to.use) == "facet")
          facet.cols <- c(facet.cols, which(colnames(meta.to.use) == "comparison"))
          n.col.ord <- c(facet.cols, (1:ncol(meta.to.use))[-facet.cols])
          if(data.table::is.data.table(meta.to.use)){
            meta.to.use <- meta.to.use[,n.col.ord, with = FALSE]
          }
          else{
            meta.to.use <- meta.to.use[,n.col.ord]
          }
        }
      }
      
      
      # run the looping function for each row of the task list
      if(par == FALSE){
        out <- vector("list", nrow(task.list))
        for(q in 1:nrow(task.list)){
          out[[q]] <- run.one.loop(stats_to_use, meta.to.use, task.list, q, FALSE)
        }
      }
      
      ## same, but in parallel
      else{
        cl <- parallel::makePSOCKcluster(par)
        doParallel::registerDoParallel(cl)
        
        #prepare reporting function
        ntasks <- nrow(task.list)
        # if(verbose){
        #   progress <- function(n) cat(sprintf("Part %d out of", n), ntasks, "is complete.\n")
        #   opts <- list(progress=progress)
        # }
        # else{
        #   opts <- list()
        # }
        
        
        if(verbose){cat("Begining run.\n")}

        # run the calculations
        ## suppress warnings because you'll get wierd ... warnings. Not an issue in the non-parallel version.
        suppressWarnings(out <- foreach::foreach(q = 1:ntasks, .inorder = TRUE,
                                                 .export = "data.table", .packages = "snpR") %dopar% {
                                                   run.one.loop(stats_to_use, meta.to.use, task.list, q, TRUE)
                                                 })
        
        #release cores
        parallel::stopCluster(cl)
      }
      
      # bind the results together.
      suppressWarnings(out <- dplyr::bind_rows(out))
      colnames(out)[1:2] <- c("facet", "subfacet")
      meta.cols <- c(1,2,3,4, which(colnames(out) %in% colnames(x@snp.meta)))
      meta <- out[,meta.cols]
      meta[is.na(meta)] <- ".base"
      out <- cbind(meta, out[,-meta.cols])
      return(out)
    }
  }
  
  
  else if(case == "psamp"){
    if(req == "genotypes"){
      split_facets <- .split.facet(facets)
      names(split_facets) <- facets
      
      # options across all facets
      opts <- lapply(split_facets, function(y){
        levs <- data.frame(unique(x@snp.meta[,which(colnames(x@snp.meta) %in% y)]))
        if(ncol(levs) == 1){colnames(levs) <- y}
        return(levs)
      })
      
      # loop function
      loop.func <- function(x, opts, fname){
        out <- vector("list", nrow(opts))
        for(i in 1:nrow(opts)){
          ident <- which(apply(x@snp.meta[,colnames(opts),drop = F], 1, function(x) identical(as.character(x), as.character(opts[i,]))))
          if(length(ident) > 0){
            genotypes <- genotypes(x)[ident,, drop = FALSE]
            suppressWarnings(out[[i]] <- cbind.data.frame(snp.facet = fname,
                                                          x@sample.meta,
                                                          opts[i,,drop = F],
                                                          stat = fun(genotypes, ...),
                                                          stringsAsFactors = F))
          }
        }
        return(dplyr::bind_rows(out))
      }
      
      # run
      out <- vector("list", length(facets))
      for(i in 1:length(facets)){
        if(facets[i] == ".base"){
          out[[i]] <- cbind.data.frame(facet = ".base", subfacet = ".base", snp.facet = ".base", snp.subfacet = ".base",
                                       x@sample.meta,
                                       stat = fun(x, ...),
                                       stringsAsFactors = F)
        }
        else{
          out[[i]] <- loop.func(x, opts[[i]], names(opts)[i])
          opt.names <- .check.snpR.facet.request(x, paste0(colnames(opts[[i]]), collapse = "."), remove.type = "sample")
          opt.names <- unlist(.split.facet(opt.names))
          out[[i]]$snp.subfacet <- do.call(paste, c(opts[[i]][,opt.names, drop = F], sep = "."))
          out[[i]]$facet <- ".base"
          out[[i]]$subfacet <- ".base"
          out[[i]] <- out[[i]][,c("facet", "subfacet", "snp.facet", "snp.subfacet", colnames(x@sample.meta), "stat")]
        }
      }
      suppressWarnings(out <- dplyr::bind_rows(out))
      
      # re-order and fix NAs from missing stats
      stat.col <- (which(colnames(out) == "stat"))
      out <- out[,c((1:ncol(out))[-stat.col], stat.col)]
      meta <- out[,-ncol(out)]
      meta[,c("facet", "subfacet", "snp.facet", "snp.subfacet")][is.na(meta[,c("facet", "subfacet", "snp.facet", "snp.subfacet")])] <- ".base"
      out[,1:(ncol(out) - 1)] <- meta
      return(out)
    }
  }
  
}

# Merge newly calculated stats into a snpRdata object.
# 
# Internal function to quickly and accurately merge newly calculated statistics
# into a snpRdata object. This should never be called externally. Takes and
# returns the same snpRdata object, with new data.
# 
# Mostly relies on the .smart.merge subfunction, which must be edited with care
# and testing. LD merging is independant.
# 
# Type options: stats, pairwise, window.stats, pairwise.window.stats, and LD,
# corresponding to slots of the snpRdata S4 object class.
# 
# For examples, see functions that use .merge.snpR.stats, such as calc_pi or
# calc_pairwise_ld.
# 
# @param x snpRdata object
# @param stats data.frame/table/list. New data to be added to the existing
#   snpRdata object, x.
# @param type character, default "stats". The type of statistic to merge, see
#   list in description.
.merge.snpR.stats <- function(x, stats, type = "stats"){
  .snp.id <- facet <- subfacet <- comparison <- NULL
  
  # note: not my code, pulled from:
  # https://stackoverflow.com/questions/51263678/r-combine-arbitrary-lists-element-by-element-with-respect-to-matched-names
  # should work at all depths?
  simple.merge.lists <- function(x, y){
    if(is.list(x) && is.list(y) && !is.null(names(x)) && !is.null(names(y))){
      ecom <- intersect(names(x), names(y))
      enew <- setdiff(names(y), names(x))
      res <- x
      if(length(enew) > 0){
        res <- c(res, y[enew])
      }
      if(length(ecom) > 0){
        for(i in ecom){
          res[i] <- list(simple.merge.lists(x[[i]], y[[i]]))
        }
      }
      return(res)
    }else{
      return(y)
    }
  }
  
  if(type == "stats"){
    # merge and return
    meta.cols <- c(colnames(stats)[1:(which(colnames(stats) == ".snp.id"))], colnames(x@snp.meta))
    starter.meta <- c("facet", "subfacet", "facet.type")
    n.s <- .smart.merge(stats, x@stats, meta.cols, starter.meta)
    x@stats <- n.s
  }
  else if(type == "pairwise"){
    # merge and return
    meta.cols <- c(colnames(stats)[1:(which(colnames(stats) == "comparison"))], colnames(x@snp.meta))
    starter.meta <- meta.cols
    n.s <- .smart.merge(stats, x@pairwise.stats, meta.cols, starter.meta)
    x@pairwise.stats <- n.s
  }
  else if(type == "window.stats"){
    # merge and return
    meta.cols <- c("facet", "subfacet", "position", "sigma", "snp.facet", "snp.subfacet", "step", "gaussian", "nk.status", "n_snps", "triple_sigma")
    starter.meta <- c("facet", "subfacet", "snp.facet", "snp.subfacet", "position")
    n.s <- .smart.merge(stats, x@window.stats, meta.cols, starter.meta)
    x@window.stats <- n.s
  }
  else if(type == "pairwise.window.stats"){
    meta.cols <- c("facet", "subfacet", "position", "sigma", "snp.facet", "snp.subfacet", "step", "gaussian", "nk.status", "n_snps", "triple_sigma")
    starter.meta <- c("facet", "subfacet", "snp.facet", "snp.subfacet", "position")
    n.s <- .smart.merge(stats, x@pairwise.window.stats, meta.cols, starter.meta)
    x@pairwise.window.stats <- n.s
  }
  else if(type == "sample.stats"){
    meta.cols <- c("facet", "subfacet", colnames(stats)[1:(which(colnames(stats) == ".sample.id"))])
    meta.cols <- unique(meta.cols)
    starter.meta <- meta.cols
    n.s <- .smart.merge(stats, x@sample.stats, meta.cols, starter.meta)
    
    # fix NAs that result from adding new facets
    meta.cols <- which(colnames(n.s) %in% colnames(x@snp.meta))
    for (j in meta.cols){
      set(n.s, which(is.na(n.s[[j]])), j , ".base")
    }
    x@sample.stats <- n.s
  }
  else if(type == "LD"){
    
    if(length(x@pairwise.LD) == 0){
      x@pairwise.LD <- stats
    }
    else{
      # deal with prox table using .smart.merge
      start.meta <- colnames(stats$prox)[1:which(colnames(stats$prox) == "proximity")]
      x@pairwise.LD$prox <- .smart.merge(stats$prox, x@pairwise.LD$prox,
                                        meta.names = c(start.meta, "sample.facet", "sample.subfacet"),
                                        starter.meta = start.meta)
      
      # Deal with matrices using the merge.lists utility in snpR.
      # a <- .merge.lists(x@pairwise.LD$LD_matrices, stats$LD_matrices)
      x@pairwise.LD$LD_matrices <- simple.merge.lists(x@pairwise.LD$LD_matrices, stats$LD_matrices)
    }
  }
  else if(type == "genetic_distances"){
    
    if(length(x@genetic_distances) == 0){
      x@genetic_distances <- stats
    }
    else{
      x@genetic_distances <- simple.merge.lists(x@genetic_distances, stats)
    }
  }
  else if(type == "allele_frequency_matrices"){
    
    if(length(x@allele_frequency_matrices) == 0){
      x@allele_frequency_matrices <- stats
    }
    else{
      x@allele_frequency_matrices <- simple.merge.lists(x@allele_frequency_matrices, stats)
    }
    
  }
  else if(type == "ibd"){
    if(length(x@other$ibd) == 0){
      x@other$ibd <- stats
    }
    else{
      x@other$ibd <- simple.merge.lists(x@other$ibd, stats)
    }
  }
  else if(type == "geo_dists"){
    if(length(x@other$geo_dists) == 0){
      x@other$geo_dists <- stats
    }
    else{
      x@other$geo_dists <- simple.merge.lists(x@other$geo_dists, stats)
    }
  }
  else if(type == "pop"){
    meta.names <- c("facet", "pop")
    starter.meta <- meta.names
    x@pop.stats <- .smart.merge(stats, x@pop.stats, meta.names, starter.meta)
  }
  else if(type == "weighted.means"){
    meta.names <- c("facet", "subfacet", "snp.facet", "snp.subfacet", colnames(x@facet.meta)[-c(1:3, ncol(x@facet.meta))])
    starter.meta <- meta.names
    x@weighted.means <- .smart.merge(stats, x@weighted.means, meta.names, starter.meta)
  }
  
  return(x)
}

# Merging function used in .merge.snpR.stats and elsewhere
# @param n.s new stats
# @param o.s old stats
# @param meta.names names of the metadata columns, usually everything up to .snp.id
# @param starter.meta any metadata columns that should specifically be put at the start of the output data (such as facet, subfacet, facet.type)
.smart.merge <- function(n.s, o.s, meta.names, starter.meta){
  ..meta.cols <- ..col.ord <- NULL
  
  # subfunction to sort by starter meta, then return the new data without respect to old. Used if old is empty or contains identical data.
  take_new <- function(n.s, starter.meta){
    smc <- which(colnames(n.s) %in% starter.meta)
    if(length(smc) > 0){
      smc <- factor(colnames(n.s)[smc], levels = starter.meta)
      smc <- sort(smc)
      smc <- as.character(smc)
      
      new.ord <- c(match(smc, colnames(n.s)), which(!colnames(n.s) %in% starter.meta))
      n.s <- .fix..call(n.s[,..new.ord])
    }
    return(n.s)
  }
  
  
  .snp.id <- facet <- subfacet <- comparison <- ..new.ord <- NULL
  
  
  n.s <- data.table::as.data.table(n.s)
  o.s <- data.table::as.data.table(o.s)
  
  if(all(colnames(n.s) %in% colnames(o.s))){
    if(isTRUE(all.equal(n.s, o.s, ignore.col.order = T, ignore.row.order = T, check.attributes = F))){
      return(take_new(n.s, starter.meta))
    }
  }
  
  if(nrow(o.s) == 0){
    return(take_new(n.s, starter.meta))
  }
  
  # figure out which columns contain metadata
  meta.cols.n <- which(colnames(n.s) %in% meta.names)
  meta.cols.n <- colnames(n.s)[meta.cols.n]
  meta.cols.o <- which(colnames(o.s) %in% meta.names)
  meta.cols.o <- colnames(o.s)[meta.cols.o]
  meta.cols <- unique(c(meta.cols.n, meta.cols.o))
  
  # add in any missing metadata columns to both o.s and n.s
  new.meta <- meta.cols.n[which(!(meta.cols.n %in% meta.cols.o))]
  if(length(new.meta) > 0){
    fill <- matrix(".base", nrow = nrow(o.s), ncol = length(new.meta))
    colnames(fill) <- new.meta
    o.s <- cbind(o.s[,starter.meta, with = FALSE], fill, o.s[,-starter.meta, with = FALSE])
  }
  old.meta <- meta.cols.o[which(!(meta.cols.o %in% meta.cols.n))]
  if(length(old.meta) > 0){
    fill <- matrix(".base", nrow = nrow(n.s), ncol = length(old.meta))
    colnames(fill) <- old.meta
    n.s <- cbind(n.s[,starter.meta, with = FALSE], fill, n.s[,-starter.meta, with = FALSE])
  }
  
  ## make sure metadata columns are sorted identically
  new.meta <- n.s[,meta.cols, with = FALSE]
  n.ord <- match(colnames(new.meta), meta.cols)
  new.meta <- new.meta[,n.ord, with = FALSE]
  old.meta <- o.s[,meta.cols, with = FALSE]
  o.ord <- match(colnames(old.meta), meta.cols)
  old.meta <- old.meta[,o.ord, with = FALSE]
  if(ncol(o.s) - length(old.meta) != 0){
    o.s <- cbind(old.meta, o.s[,-meta.cols, with = FALSE])
  }
  if(ncol(n.s) - length(new.meta) != 0){
    n.s <- cbind(new.meta, n.s[,-meta.cols, with = FALSE])
  }
  
  ## check classes in meta to make sure there are no conflicts
  o.class <- unlist(lapply(.fix..call(o.s[,..meta.cols]), class))
  n.class <- unlist(lapply(.fix..call(n.s[,..meta.cols]), class))
  miss.class <- which(o.class != n.class)
  if(length(miss.class) > 0){
    for(i in 1:length(miss.class)){
      tfix <- meta.cols[miss.class[i]]
      n.s[[tfix]] <- methods::as(n.s[[tfix]], unlist(o.class[miss.class[i]]))
    }
  }
  
  ## do the merge, then fix the .y and .x columns by replacing NAs in y with their value in x
  m.s <- merge(o.s, n.s, by = meta.cols, all = T)
  ### grab any columns that need to be fixed (end in .x or .y) and save any matching columns that are fine as is.
  match.cols <- colnames(m.s)[which(colnames(m.s) %in% c(colnames(o.s), colnames(n.s)))]
  stat.cols.to.fix <- m.s[,-match.cols, with = FALSE]
  stat.cols.y <- grep("\\.y$", colnames(stat.cols.to.fix))
  if(length(stat.cols.y) > 0){
    # replace NAs in y with values from x. No reason to do the vice versa.
    stat.cols.y <- as.matrix(stat.cols.to.fix[,stat.cols.y, with = FALSE])
    stat.cols.x <- grep("\\.x$", colnames(stat.cols.to.fix))
    stat.cols.x <- as.matrix(stat.cols.to.fix[,stat.cols.x, with = FALSE])
    stat.cols.x <- stat.cols.x[,order(colnames(stat.cols.x)), drop = FALSE]
    stat.cols.y <- stat.cols.y[,order(colnames(stat.cols.y)), drop = FALSE]
    NA.y <- is.na(stat.cols.y)
    stat.cols.y[NA.y] <- stat.cols.x[NA.y]
    colnames(stat.cols.y) <- gsub("\\.y$", "", colnames(stat.cols.y))
    
    # update m.s
    m.s <- cbind(m.s[,match.cols, with = FALSE], stat.cols.y)
    
    
    # fix column classes
    for(j in 1:ncol(m.s)){
      if(is.character(m.s[[j]]) & !(colnames(m.s)[j] %in% meta.names)){
        set(m.s, j = j, value = tryCatch(as.numeric(m.s[[j]]), warning = function(x) m.s[[j]]))
      }
    }
  }
  
  # sort by .snp.id if that's a thing here.
  if(any(colnames(m.s) == ".snp.id")){
    if("subfacet" %in% colnames(m.s)){
      m.s <- data.table::setorder(m.s, .snp.id, facet, subfacet)
    }
    else{
      m.s <- data.table::setorder(m.s, .snp.id, facet, comparison)
    }
  }
  
  # ensure that the starter meta is in the correct spot
  col.ord <- unique(c(starter.meta, meta.cols, colnames(m.s)))
  col.ord <- col.ord[which(col.ord %in% colnames(m.s))]
  m.s <- .fix..call(m.s[,..col.ord])
  
  return(m.s)
}

# Get details on and clean up facet arguments.
# 
# Given a snpRdata object, this function cleans up and determines the types of
# requested facets. Internal. If requested, can clean facet requests to purge a
# specific facet type (snp, sample, complex). Used in functions that run only on
# a specific type of facet.
# 
# Facet designation of NULL or "all" follows the typical rules.
# 
# Facets designated as described in \code{\link{Facets_in_snpR}}.
# 
# @param x snpRdata object to compare to
# @param facets character. facets to check.
# @param remove.type character. "snp", "sample", "complex", "simple" or anything
#   else, typically "none". Type of facet to remove.
# @param return.type logical, default FALSE. If true, returns both facets and
#   facet types ("snp", "sample", "complex", or ".base") as a length two list.
# @param fill_with_base logical, default TRUE. If true, fills the returned facets
#   with .base if nothing is left after facets are removed. Otherwise returns null in this case.
# @param return_base_when_empty logical, default TRUE. If not facets pass filters,
#   will pass .base if true.
# 
# @author William Hemstrom
.check.snpR.facet.request <- function(x, facets, remove.type = "snp", return.type = FALSE, 
                                      fill_with_base = TRUE, return_base_when_empty = TRUE, purge_duplicates = TRUE){

  if(any(facets == "all")){
    facets <- x@facets
  }
  if(is.null(facets) | isFALSE(facets) | length(facets) == 0){
    facets <- ".base"
    if(return.type){
      return(list(facets, ".base"))
    }
    else{
      return(".base")
    }
  }
  
  
  # remove the facet parts as requested.
  facets <- .split.facet(facets)
  to.remove <- logical(length(facets))
  missing.facets <- character(0)
  facet.types <- character(length(facets))
  for(i in 1:length(facets)){
    if(identical(facets[[i]], ".base")){next()}
    facets[[i]] <- sort(facets[[i]])
    samp.facets <- which(facets[[i]] %in% colnames(x@sample.meta))
    snp.facets <- which(facets[[i]] %in% colnames(x@snp.meta))
    missing.facets <- c(missing.facets, facets[[i]][which(!((1:length(facets[[i]])) %in% c(samp.facets, snp.facets)))])
    
    if(remove.type == "snp"){
      if(length(snp.facets) > 0){
        if(length(snp.facets) == length(facets[[i]])){
          to.remove[i] <- T
        }
        else{
          facets[[i]] <- facets[[i]][-snp.facets]
        }
      }
    }
    
    else if(remove.type == "sample"){
      if(length(samp.facets) > 0){
        if(length(samp.facets) == length(facets[[i]])){
          to.remove[i] <- T
        }
        else{
          facets[[i]] <- facets[[i]][-samp.facets]
        }
      }
    }
    
    else if(remove.type == "complex"){
      if(length(snp.facets) > 0 & length(samp.facets) > 0){
        to.remove[i] <- T
      }
    }
    
    else if(remove.type == "simple"){
      if(!(length(snp.facets) > 0 & length(samp.facets) > 0)){
        to.remove[i] <- T
      }
    }
    
    if(return.type){
      samp.facets <- which(facets[[i]] %in% colnames(x@sample.meta))
      snp.facets <- which(facets[[i]] %in% colnames(x@snp.meta))
      missing.facets <- c(missing.facets, facets[[i]][which(!((1:length(facets[[i]])) %in% c(samp.facets, snp.facets)))])
      
      if(length(samp.facets) > 0 & length(snp.facets) > 0){
        facet.types[i] <- "complex"
      }
      else if(length(samp.facets) > 0){
        facet.types[i] <- "sample"
      }
      else if(length(snp.facets) > 0){
        facet.types[i] <- "snp"
      }
    }
    
    # check for bad characters in our facet level (slows things down a small bit, but avoids user errors that are uninformative)
    bad.chars <- c("\\.", "\\~", " ")
    bad.chars.rep <- c("'.'", "'~'", "' ' (spaces)")
    for(j in 1:length(facets[[i]])){
      if(facets[[i]][j] %in% colnames(x@snp.meta)){
        bad <- lapply(bad.chars, function(y) any(grepl(y, x@snp.meta[,facets[[i]][j]])))
      }
      else if(facets[[i]][j] %in% colnames(x@sample.meta)){
        bad <- lapply(bad.chars, function(y) any(grepl(y, x@sample.meta[,facets[[i]][j]])))
      }
      else if(facets[[i]][j] == ".base"){
        next
      }
      else{
        # if the facet doesn't exit, an error will pop up later
        next
      }
      
      bad <- which(unlist(bad))
      if(length(bad) > 0){stop(paste0("Facet '", 
                                      facets[[i]][j], 
                                      "' cannot be used as a facet level since it contains: \n\t", 
                                      paste0(bad.chars.rep[bad], collapse = "\n\t ")))}
    }
    
    # check for duplicates in our multi-part facets (also slows down but helps catch errors)
    if(length(facets[[i]]) > 0){
      d <- vector("list", 0)
      for(j in 1:length(facets[[i]])){
        if(facets[[i]][j] %in% colnames(x@snp.meta)){
          d[length(d) + 1] <- x@snp.meta[,facets[[i]][j], drop = FALSE]
          names(d)[length(d)] <- facets[[i]][j]
        }
        else if(facets[[i]][j] %in% colnames(x@sample.meta)){
          d[length(d) + 1] <- x@sample.meta[,facets[[i]][j], drop = FALSE]
          names(d)[length(d)] <- facets[[i]][j]
        }
      }
      
      uniques <- lapply(d, unique)
      un <- rep(names(uniques), lengths(uniques))
      uniques <- unlist(uniques, use.names = FALSE)
      names(uniques) <- un
      rm(un)
      if(any(duplicated(uniques))){
        dups <- sort(uniques[which(duplicated(uniques) | duplicated(uniques, fromLast = TRUE))])
        msg <- unique(dups)
        pst_msg <- "Some levels are duplicated across multiple facets.\nThis will cause issues if those facets are run during analysis.\nIssues:\n"
        for(q in 1:length(msg)){
          pst_msg <- paste0(pst_msg, "\nLevel: ", msg[q], "\tin facets: ", paste0(names(dups)[which(dups == msg[q])], collapse = ", "))
        }
        stop(pst_msg)
      }
    }
    
    # update
    facets[[i]] <- paste(facets[[i]], collapse = ".")
  }
  
  if(length(missing.facets) > 0){
    dups <- which(duplicated(missing.facets))
    if(length(dups) > 0){
      stop("Facet(s) ", paste(missing.facets[-dups], collapse = ", "), " not found in x metadata.\n")
    }
    stop("Facet(s) ", paste(missing.facets, collapse = ", "), " not found in x metadata.\n")
  }
  
  if(any(to.remove)){
    if(fill_with_base){
      facets[which(to.remove)] <- ".base"
      facet.types <- facet.types[-which(to.remove)] <- ".base"
    }
    else{
      facets <- facets[-which(to.remove)]
      facet.types <- facet.types[-which(to.remove)]
    }
  }
  
  # fix if we've removed everything, return the base facet if return_base_when_empty is TRUE
  if(return_base_when_empty){
    if(length(facets) == 0){
      facets <- ".base"
      if(return.type){
        return(list(facets, ".base"))
      }
      else{
        return(".base")
      }
    }
  }
  
  # fix the facet type for .base
  if(return.type){
    base.facets <- which(facet.types == "")
    if(length(base.facets) > 0){
      facet.types[base.facets] <- ".base"
    }
  }
  
  # remove duplicates
  if(purge_duplicates){
    dups <- duplicated(facets)
    if(any(dups)){
      facets <- facets[-which(dups)]
      facet.types <- facet.types[-which(dups)]
    }
  }
  

  # return
  if(return.type){
    return(list(unlist(facets), facet.types))
  }
  else{
    return(unlist(facets))
    
  }
}

# Tabulate allele and genotype counts at each locus.
# 
# \code{.tabulate_genotypes} creates matricies containing counts of observed
# alleles and genotypes at each locus.
# 
# This function is pirmarily used interally in several other funcitons.
# 
# @param x Input raw genotype data, where columns are individuals and rows are
#   snps. No metadata.
# @param mDat Character. How are missing \emph{genotypes} noted?
# @param verbose Logical. Should the function report progress?
# 
# @return A list of matrices. gs is the genotype matrix, as is the allele
#   matrix, and wm is the genotype matrix with missing genotypes.
# 
# @author William Hemstrom
# 
.tabulate_genotypes <- function(x, mDat, verbose = F){
  ..all_idents <- ..rm_cols <- NULL
  
  .tab_func <- function(x, snp_form, mDat){
    ..mDat <- NULL
    x <- data.table::melt(data.table::transpose(x, keep.names = "samp"), id.vars = "samp") # transpose and melt
    
    gmat <- data.table::dcast(data.table::setDT(x), variable ~ value, value.var='value', length) # cast
    gmat <- gmat[,-1]
    
    # fix any cases where we have the same genotype but reversed allele order TA is the same as AT
    opts <- colnames(gmat)
    opts1 <- substr(opts, 1, snp_form/2)
    opts2 <- substr(opts, (snp_form/2 + 1), snp_form*2)
    rev_opts <- paste0(opts2, opts1)
    reps <- which(opts %in% rev_opts & opts1 != opts2)
    if(length(reps) > 0){
      done <- character()
      rm_cols <- numeric()
      for(i in 1:length(reps)){
        this_opt <- reps[i]
        if(opts[this_opt] %in% done){next}
        match <- try(which(rev_opts == opts[this_opt]), silent = TRUE)
        if(methods::is(match, "try-error")){stop("Error at match.\n")}
        all_idents <- c(this_opt, match)
        gmat[[this_opt]] <- .fix..call(rowSums(gmat[,..all_idents]))
        done <- c(done, unique(opts[all_idents]))
        rm_cols <- c(rm_cols, all_idents[-1])
      }
      
      gmat <- .fix..call(gmat[,-..rm_cols])
    }
    
    # make sure column names are properly sorted by allele
    opts <- colnames(gmat)
    opts1 <- substr(opts, 1, snp_form/2)
    opts2 <- substr(opts, (snp_form/2 + 1), snp_form*2)
    flips <- try(!opts1 <= opts2, silent = TRUE) # this'll do an alphabetic check since these are characters--weird but works
    if(is(flips, "try-error")){stop("Error at flips.\n")}
    if(any(flips)){
      colnames(gmat)[which(flips)] <- paste0(opts2[flips], opts1[flips])
    }
    
    # remove missing
    mis.cols <- -which(colnames(gmat) == mDat)
    if(length(mis.cols) > 0){
      tmat <- gmat[,mis.cols, with = FALSE] # remove missing data
    }
    else{
      tmat <- gmat
    }
    if(nrow(tmat) == 0){
      if(any(colnames(gmat) == mDat)){
        wm <- Matrix::Matrix(as.matrix(.fix..call(gmat[,..mDat])), sparse = TRUE)
      }
      else{
        wm <- Matrix::Matrix(0, nrow(gmat), 1, dimnames = list(NULL, mDat), sparse = TRUE)
      }
      
      return(list(gs = Matrix::Matrix(0, nrow = 0, ncol = 0, sparse = TRUE), 
                  as = Matrix::Matrix(0, nrow = 0, ncol = 0, sparse = TRUE), 
                  wm = wm))
    }
    
    #get matrix of allele counts
    #initialize
    hs <- substr(colnames(tmat),1,snp_form/2) != substr(colnames(tmat), (snp_form/2 + 1), snp_form*2) # identify heterozygotes.
    if(verbose){cat("Getting allele table...\n")}
    split_pattern <- paste0("(?<=", paste0(rep(".", snp_form/2), collapse = ""), ")") # split by the correct allele size
    as <- unique(unlist(strsplit(paste0(colnames(tmat)), split_pattern, perl = TRUE)))
    amat <- data.table::as.data.table(matrix(0, nrow(gmat), length(as)))
    colnames(amat) <- as

    #fill in
    for(i in 1:length(as)){
      b <- grep(paste0(as[i], "$"), colnames(tmat))
      b <- sort(unique(c(b, grep(paste0("^", as[i]), colnames(tmat)))))
      hom <- which(colnames(tmat) == paste0(as[i], as[i]))
      if(length(hom) == 0){
        het <- b
        set(amat, j = i, value = rowSums(tmat[,het, with = FALSE]))
      }
      else{
        het <- b[b != hom]
        if(length(het) > 0){
          if(data.table::is.data.table(tmat[,het, with = FALSE])){
            set(amat, j = i, value = (tmat[,hom, with = FALSE] * 2) + rowSums(tmat[,het, with = FALSE]))
          }
          else{
            amat[,i] <- (tmat[,hom] * 2) + tmat[,het]
          }
        }
        else{
          set(amat, j = i, value = (tmat[,hom, with = FALSE] * 2))
        }
      }
    }
    # convert to sparse matrices
    tmat <- Matrix::Matrix(as.matrix(tmat), sparse = TRUE)
    amat <- Matrix::Matrix(as.matrix(amat), sparse = TRUE)
    
    if(any(colnames(gmat) == mDat)){
      gmat <- .fix..call(gmat[,..mDat])
      gmat <- Matrix::Matrix(as.matrix(gmat), sparse = TRUE)
    }
    else{
      gmat <- Matrix::Matrix(0, nrow(gmat), 1, dimnames = list(NULL, mDat), sparse = TRUE)
    }
    
    
    return(list(gs = tmat, as = amat, wm = gmat))
  }

  # fix for if x is a vector (only one individual) and convert to data.table
  if(!is.data.frame(x)){
    x <- data.frame(samp = x)
  }
  x <- data.table::setDT(x)
  
  
  
  
  # get a genotype table
  snp_form <- nchar(x[1,1])   # get information on data format
  
  
  # loop to save memory if large
  max_rows <- 100000000 # about 1.5G of storage for the big melt
  max_snps <- ceiling(max_rows/ncol(x)) # max snps tolerable at once
  n_iters <- ceiling(nrow(x)/max_snps)
  titer <- 1
  
  # for each set of snps, get the weighted values and weights for all relevent windows.
  geno.tables <- vector("list", n_iters)
  for(i in 1:n_iters){
    end <- i*max_snps
    end <- ifelse(end > nrow(x), nrow(x), end)
    trows <- titer:end
    
    geno.tables[[i]] <- .tab_func(x[trows,], snp_form, mDat)

    titer <- i*max_snps+ 1
  }

  gs <- purrr::map(geno.tables, "gs")
  as <- purrr::map(geno.tables, "as")
  wm <- purrr::map(geno.tables, "wm")
  
  gs <- .rbind_list_fill_sparse(gs)
  as <- .rbind_list_fill_sparse(as)
  wm <- .rbind_list_fill_sparse(wm)
  
  return(list(gs = gs, as = as, wm = wm))
}




# Interpolate sn formatted data.
# 
# An internal function to interpolate sn formatted data using either the
# bernoulli or expected minor allele count approaches. Typically entirely
# internal, called via format_snps.
# 
# Interpolating missing data in sn formatted data is useful for PCA, genomic
# prediction, tSNE, and other methods. Specify method = "af" to insert the
# expected number of minor alleles given SNP allele frequency or "bernoulli" to
# do binomial draws to determine the number of minor alleles at each missing
# data point, where the probability of drawing a minor allele is equal to the
# minor allele frequency. The expected number of minor alleles based on the
# later method is equal to the interpolated value from the former, but the later
# allows for multiple runs to determine the impact of stochastic draws and is
# generally prefered and required for some downstream analysis. It is therefore
# the default. As a slower but more accurate alternative to "af" interpolation,
# "iPCA" may be selected. This an iterative PCA approach to interpolate based on
# SNP/SNP covariance via \code{\link[missMDA]{imputePCA}}. If the ncp argument
# is not defined, the number of components used for interpolation will be
# estimated using \code{\link[missMDA]{estim_ncpPCA}}. In this case, this method
# is much slower than the other methods, especially for large datasets. Setting
# an ncp of 2-5 generally results in reasonable interpolations without the time
# constraint.
# 
# @param sn data.frame. Input sn formatted data, as produced by
#   \code{\link{format_snps}}. Note that \code{\link{format_snps}} has an option
#   to automatically call this function during formatting.
# @param method character, default "bernoulli". Method to used for
#   interpolation, either bernoulli or af. See details.
# @param ncp numeric or NULL, default NULL. Number of components to consider for
#   iPCA sn format interpolations of missing data. If null, the optimum number
#   will be estimated, with the maximum specified by ncp.max. This can be very
#   slow.
# @param ncp.max numeric, default 5. Maximum number of components to check for
#   when determining the optimum number of components to use when interpolating
#   sn data using the iPCA approach.
# 
# @author William Hemstrom
.interpolate_sn <- function(sn, method = "bernoulli", ncp = NULL, ncp.max = 5){
  if(!method %in% c("bernoulli", "af", "iPCA")){
    stop("Unaccepted interpolation method. Accepted methods: bernoulli, af, IPCA.\n")
  }

  if(method %in% c("bernoulli", "af")){
    # find allele frequencies
    sn <- as.matrix(sn)
    sn <- t(sn)
    af <- colSums(sn, na.rm = T)
    af <- af/(colSums(!is.na(sn))*2) # this is the allele frequency of the "1" allele
    
    # identify all of the NAs and the columns that they belong to
    NAs <- which(is.na(sn)) # cells with NA
    NA.cols <- floor(NAs/nrow(sn)) + 1 # figure out the column
    adj <- which(NAs%%nrow(sn) == 0) # adjust for anything that sits in the last row, since it'll get assigned the wrong column
    NA.cols[adj] <- NA.cols[adj] - 1
    
    # do interpolation for each missing data point
    if(method == "bernoulli"){
      ndat <- stats::rbinom(length(NAs), 2, af[NA.cols])
    }
    else if(method == "af"){
      ndat <- af[NA.cols]*2
    }
    
    # replace
    sn[NAs] <- ndat
    
    # remove any columns (loci) with NO data and warn!
    no_dat_loci <- which(is.na(af))
    
    if(length(no_dat_loci) > 0){
      sn <- sn[,-no_dat_loci]
      warning("Some loci had no called genotypes and were removed: ", paste0(no_dat_loci, collapse = ", "), "\n")
    }
  }
  else if(method %in% c("iPCA")){
    .check.installed("missMDA")
    sn <- t(as.matrix(sn))
    
    # remove any columns (loci) with NO data and warn!
    no_dat_loci <- which(colSums(is.na(sn)) == nrow(sn))
    
    if(length(no_dat_loci) > 0){
      sn <- sn[,-no_dat_loci]
      warning("Some loci had no called genotypes and were removed: ", paste0(no_dat_loci, collapse = ", "), "\n")
    }
    
    
    # determine optimal np from 0 - max.ncp
    if(is.null(ncp)){
      ncp <- missMDA::estim_ncpPCA(sn, ncp.max = ncp.max, method.cv = "Kfold")
      cat("ncp scores:\n")
      print(ncp$criterion)
      ncp <- ncp$ncp
    }
    sn <- missMDA::imputePCA(sn, ncp)
  }
  return(t(sn))
}

# Do sanity checks for many window specific functions in snpR.
# 
# Internal function only. Since many window related functions have similar
# sanity checks, they have been integrated here. For examples, see use in
# functions like \code{\link{do_bootstraps}}.
# 
# @param x snpRdata object.
# @param sigma numeric. Window size, in kilobases.
# @param step numeric or NULL. Step size, in kilobases
# @param stats.type character, either "single" or "pairwise". The type of stats
#   being checked, see documentation for \code{\link{calc_smoothed_averages}}.
# @param nk Logical. Determines if nk values are to be used.
# @param facets character or NULL. Defines facets to run, following typical
#   rules as described in \code{\link{Facets_in_snpR}}.
# @param stats character or NULL, default NULL. Statistics (pi, etc.) used,
#   really only for bootstrapping.
# @param good.types character or NULL, default NULL. Good statistics, for use
#   with the stats arguement.
# 
# @author William Hemstrom
.sanity_check_window <- function(x, sigma, step, stats.type, nk, facets, stats = NULL, good.types = NULL){
  #================sanity checks=============
  msg <- character()
  warn.msg <- character()
  
  #nk
  if(!all(is.logical(nk) & length(nk) == 1)){
    msg <- "nk must be TRUE or FALSE."
  }
  
  #smoothing method
  good.methods <- c("single", "pairwise")
  if(sum(good.methods %in% stats.type) < 1){
    msg <- c(msg, paste0("No accepted stats types provided. Acceptable types:\n\t", paste0(good.methods, collapse = "\t")))
  }
  bad.stats <- which(!(stats.type %in% good.methods))
  if(length(bad.stats) > 0){
    msg <- c(msg, paste0("Unaccepted stats types provided:\n\t", paste0(stats.type[bad.stats], collapse = "\t")))
  }
  
  #sigma and ws
  if(!is.null(step)){
    if(!all(is.numeric(sigma), is.numeric(step), length(sigma) == 1, length(step) == 1)){
      msg <- c(msg, "sigma must be a numeric vector of length 1. step may be the same or NULL.")
    }
    if(sigma >= 500 | sigma >= 500){
      warn.msg <- c(warn.msg, "Sigma and/or ws are larger than typically expected. Reminder: sigma and ws are given in kilobases!")
    }
    else if(sigma <= 10){
      warn.msg <- c(warn.msg, paste0("Sigma is smaller that typically expected: ", sigma, " kb!"))
    }
  }
  else{
    if(!all(is.numeric(sigma), length(sigma) == 1)){
      msg <- c(msg, "sigma must be a numeric vector of length 1. step may be the same or NULL.")
    }
    if(sigma >= 500){
      warn.msg <- c(warn.msg, "Sigma is larger than typically expected. Reminder: sigma is given in kilobases!")
    }
    else if(sigma <= 10){
      warn.msg <- c(warn.msg, paste0("Sigma is smaller that typically expected: ", sigma, " kb!"))
    }
  }
  
  # position needs to be available
  if(!any(colnames(x@snp.meta) == "position")){
    msg <- c(msg, "No column named position found in snp metadata.")
  }
  else{
    if(any(is.na(snp.meta(x)$position))){msg <- c(msg, "NA values found in position data in snp metadata.\n")}
  }
  
  # facets
  if(is.null(facets[1]) & any(stats.type == "pairwise")){
    msg <- c(msg, "If no facets are provided, pairwise statistics cannot be smoothed. Please specify stats.type or statistics = 'single'")
  }
  facet.types <- .check.snpR.facet.request(x, facets, remove.type = "none", return.type = T)
  if(any(facet.types[[2]] == "snp") & any(stats.type == "pairwise")){
    msg <- c(msg, "If snp facets are provided, pairwise statistics cannot be smoothed. Please specify stats.type or statistics = 'single'")
  }
  
  # statistics, really only for bootstrapping
  if(!is.null(stats[1])){
    bad.stats <- which(!(stats %in% good.types))
    if(length(bad.stats) > 0){
      msg <- c(msg, paste0("Some of the requested statics are not acceptable for bootstrapping: ", paste0(stats[bad.stats], collapse = " "),
                           "Acceptable stats: ", paste0(good.types, collapse = " ")))
    }
    
    missing.stats <- which(!(stats %in% c(colnames(x@stats), colnames(x@pairwise.stats))))
    if(length(missing.stats) > 0){
      msg <- c(msg, paste0("Some of the statistics for which smoothing have been suggested have not been yet been calculated: ", paste0(stats[missing.stats], collapse = " "),
                           ". Please run these statistics for the supplied facets!"))
    }
  }
  
  # cannot smooth on .base with pairwise
  if("pairwise" %in% stats.type){
    sample.facets <- .check.snpR.facet.request(x, facets, "snp")
    if(".base" %in% sample.facets){
      msg <- c(msg, "Cannot conduct pairwise smoothing with a .base sample facet\n.")
    }
  }
  
  if(length(warn.msg) > 0){
    warning(warn.msg)
  }
  if(length(msg) > 0){
    stop(paste0(msg, collapse = "\n  "))
  }
  
  return(TRUE)
}



# List unique tasks/options for facets. Internal function to get a list of
# tasks to run (one task per unique sample/snp level facet!). The source
# arguement specifies what kind of statistics are being grabbed.
#
# @param x snpRdata
# @param facets Facets to generate tasks for
# @param source "stats" or "pairwise.stats", default "stats". Type of
#   comparison to get jobs for. Note that the latter only works if pairwise
#   stats have actually been calculated.
# @author William Hemstrom
.get.task.list <- function(x, facets, source = "stats"){
  ..use.facet <- NULL
  
  facets <- .check.snpR.facet.request(x, facets, "none", F)
  task.list <- matrix("", ncol = 4, nrow = 0) # sample facet, sample subfacet, snp facet, snp.subfacet
  
  if(source == "stats"){
    meta.to.use <- x@facet.meta
  }
  else if (source == "pairwise.stats"){
    meta.to.use <- as.data.frame(x@pairwise.stats, stringsAsFactors = F)
    meta.to.use$subfacet <- meta.to.use$comparison
  }
  
  snp.facet.list <- vector("list", length = length(facets))
  for(i in 1:length(facets)){
    t.facet <- facets[i]
    t.facet <- unlist(.split.facet(t.facet))
    t.facet.type <- .check.snpR.facet.request(x, t.facet, remove.type = "none", return.type = T)[[2]]
    
    # sample facets
    if(any(t.facet.type == "sample")){
      t.sample.facet <- .check.snpR.facet.request(x, facets[i], remove.type = "snp")
      t.sample.meta <- meta.to.use[meta.to.use$facet == t.sample.facet, c("subfacet")]
      sample.opts <- unique(t.sample.meta)
      t.sample.meta <- meta.to.use[,c("facet", "subfacet")]
      if(is.null(nrow(sample.opts))){
        sample.opts <- as.data.frame(sample.opts, stringsAsFactors = F)
        t.sample.meta <- as.data.frame(t.sample.meta, stringsAsFactors = F)
      }
      else if(nrow(sample.opts) == 0){
        stop(paste0("Facet: ", t.sample.facet, " not added. If you are a user seeing this, then the developers made a mistake. Please report a reproduceable example to the github issues page.\n"))
      }
    }
    else{
      t.sample.facet <- ".base"
      t.sample.meta <- meta.to.use[,c("facet", "subfacet")]
      sample.opts <- data.frame(.base = ".base", stringsAsFactors = F)
    }
    
    # snp facets
    if(any(t.facet.type == "snp")){
      use.facet <- t.facet[t.facet.type == "snp"]
      meta.to.use <- as.data.table(meta.to.use)
      t.snp.meta <- .fix..call(meta.to.use[,..use.facet])
      snp.opts <- unique(t.snp.meta)
      t.snp.facet <- .check.snpR.facet.request(x, facets[i], remove.type = "sample")
      if(is.null(nrow(snp.opts))){
        snp.opts <- as.data.frame(snp.opts)
        t.snp.meta <- as.data.frame(t.snp.meta)
        colnames(snp.opts) <- t.facet[which(t.facet %in% colnames(meta.to.use))]
        colnames(t.snp.meta) <- colnames(snp.opts)
      }
      else{
        t.snp.meta <- t.snp.meta[,match(colnames(t.snp.meta), colnames(snp.opts))]
      }
    }
    else{
      t.snp.facet <- ".base"
      t.snp.meta <- as.data.frame(rep(".base", nrow(meta.to.use)))
      snp.opts <- data.frame(.base = ".base")
    }
    # get all of the possible factor/subfactor options and make the task list for this facet
    all.opts.1 <- matrix(rep(as.matrix(sample.opts), each = nrow(snp.opts)), ncol = ncol(sample.opts))
    all.opts.1 <- do.call(paste, c(as.data.frame(all.opts.1), sep = "."),)
    all.opts.2 <- matrix(rep(t(dplyr::mutate_all(snp.opts, as.character)), nrow(sample.opts)), ncol = ncol(snp.opts), byrow = TRUE) # mutate before transpose because t will occasionally add a space before a number for mysterious reasons.
    all.opts.2 <- do.call(paste, c(as.data.frame(all.opts.2), sep = "."))
    t.task.list <- cbind(t.sample.facet, all.opts.1, t.snp.facet, all.opts.2)
    task.list <- rbind(task.list, t.task.list)
  }
  
  return(task.list)
}

# Internal to find which sample metadata rows match a row from a task list.
# @param x snpRdata object to check
# @param task.list.row A row of an object created by
#   \code{\link{.get.task.list}} to look up info about.
.fetch.sample.meta.matching.task.list <- function(x, task.list.row){
  task.list.row <- unlist(task.list.row) # in case it is passed with drop = F for some reason
  
  if(task.list.row[1] == ".base"){
    matches <- 1:nrow(x@sample.meta)
  }
  else{
    t.facets <- unlist(.split.facet(task.list.row[1]))
    matches <- which(do.call(paste, c(dat = as.data.frame(x@sample.meta[,t.facets]), sep = ".")) == task.list.row[2])
  }
  
  return(matches)
}

# Internal to find which snp metadata rows match a row from a task list.
# @param x snpRdata object to check
# @param task.list.row A row of an object created by
#   \code{\link{.get.task.list}} to look up info about.
.fetch.snp.meta.matching.task.list <- function(x, task.list.row){
  task.list.row <- unlist(task.list.row) # in case it is passed with drop = F for some reason
  
  if(task.list.row[3] == ".base"){
    matches <- 1:nrow(x@snp.meta)
  }
  else{
    t.facets <- unlist(.split.facet(task.list.row[3]))
    matches <- which(do.call(paste, c(dat = as.data.frame(x@snp.meta[,t.facets]), sep = ".")) == task.list.row[4])
  }
  
  return(matches)
}



# Update list of calculated stats for a vector of facet names
#
# @param x snpRdata object to update
# @param facets character. Facets to update, see \code{\link{Facets_in_snpR}}
# @param stats character. Name of the facet to update as calculated.
# @param remove.type character, default none. If provided, will grab only the
#   remove the requested facet type (omit the snp part, for example).
#
# @return A snpRdata object identical to x but with calced stats updated.
#
# @author William Hemstrom
.update_calced_stats <- function(x, facets, stats, remove.type = "none"){
  facets <- .check.snpR.facet.request(x, facets, remove.type = remove.type)
  
  for(i in 1:length(facets)){
    
    # add a storage vector for this facet if no stats have yet been added
    if(!facets[i] %in% names(x@calced_stats)){
      x@calced_stats <- c(x@calced_stats, list(stats))
      names(x@calced_stats)[length(x@calced_stats)] <- facets[i]
    }
    
    # update list of calculated stats for this facet
    else{
      x@calced_stats[[facets[i]]] <- unique(c(x@calced_stats[[facets[i]]], stats))
    }
  }
  
  return(x)
}

# Check if a given stat or vector of stats have been calculated any number of
# facets.
#
# @param x snpRdata object to check
# @param facets character. See \code{\link{Facets_in_snpR}}
# @param stats character. Names of stats to check.
#
# @return A named list with an entry for each facet containing a named logical
#   vector indicating if the provided stats have been calculated yet.
#
# @author William Hemstrom
.check_calced_stats <- function(x, facets, stats){
  # init storage
  out <- vector("list", length(facets))
  names(out) <- facets
  
  # check each facet to see if requested stats are calculated
  for(i in 1:length(facets)){
    out[[i]] <- stats %in% x@calced_stats[[facets[i]]]
    names(out[[i]]) <- stats
  }
  return(out)
}


# Utility to check if a required package is installed and install it if needed
# @param pkg character, name of the package to check.
# @param install.type character, default "basic". Where should the package be
#   sourced from? Options: "basic": CRAN, "github": github, "bioconductor":
#   bioconductor.
# @param source character, default NULL. If install.type = 'github', gives the
#   source from which to install.
.check.installed <- function(pkg, install.type = "basic", source = NULL){
  if(!pkg %in% rownames(utils::installed.packages())){
    say <- paste0("Package '", pkg, "' not found.")
    cat(say, "")
    cat("Install? (y or n)\n")
    
    if(!interactive()){
      stop(say)
    }
    
    resp <- readLines(n = 1)
    resp <- tolower(resp)
    
    
    if(resp != "y"){
      stop(say)
    }
    else{
      if(install.type == "basic"){
        utils::install.packages(pkg)
      }
      if(install.type == "bioconductor"){
        if(!"BiocManager" %in% utils::installed.packages()){
          utils::install.packages("BiocManager")
        }
        BiocManager::install(pkg)
      }
      if(install.type == "github"){
        if(!"remotes" %in% utils::installed.packages()){
          utils::install.packages("remotes")
        }
        remotes::install_github(source)
      }
      return(TRUE)
    }
  }
  else{return(TRUE)}
}



# Correct for multiple testing.
# @param p data.frame/data.table containing p values and any facets to split by
# @param pcol name of the column containing p values
# @param levs names of the columns containing facet info to split by
# @param methods methods to use
# @param case case (overall, within each subfacet, or within each facet) at
#   which to apply corrections
.fwe_correction <- function(p, pcol = NULL, levs = c("subfacet", "facet"), methods = c("bonferroni", "holm", "BH", "BY"), 
                           case = c("overall", "by_facet", "by_subfacet")){

  if(is.null(ncol(p))){
    case <- "overall"
    p <- as.data.table(p)
    colnames(p) <- "p"
    pcol <- p
  }
  if(ncol(p) == 1){
    case <- "overall"
    pcol <- colnames(p)[1]
  }
  
  # function to run on one set of p-values
  do_correction <- function(tp, methods){
    if(!is.null(ncol(tp))){
      ord <- tp[["ord"]]
      tp <- tp[[pcol]]
    }
    
    # for the methods already implemented in R, really easy.
    sig_tab <- sapply(methods, function(x) stats::p.adjust(tp, x)) # run those methods
    if(is.null(ncol(sig_tab))){
      sig_tab <- as.data.frame(matrix(sig_tab, 1, 1))
    }
    colnames(sig_tab) <- methods
    sig_tab <- as.data.frame(sig_tab)
    
    if(exists("ord")){
      sig_tab <- cbind(sig_tab, ord = ord)
    }
    return(sig_tab)
  }
  
  p <- data.table::as.data.table(p)

  # do the overall and or case by case
  if("overall" %in% case){
    p_overall <- p[[pcol]]
    p_overall <- do_correction(p_overall, methods)
    colnames(p_overall) <- paste0(pcol, "_overall_", colnames(p_overall))
    out <- cbind(p, p_overall)
  }
  else{
    out <- p
  }
  
  # split up completely
  if("by_subfacet" %in% case){
    p$ord <- 1:nrow(p)
    p_case <- p[,do_correction(.SD, methods), by = levs]
    p_case <- p_case[order(p_case$ord),]
    p_case$ord <- NULL
    colnames(p_case)[(ncol(p_case) - length(methods) + 1):ncol(p_case)] <- 
      paste0(pcol, "_bysubfacet_", colnames(p_case)[(ncol(p_case) - length(methods) + 1):ncol(p_case)])
    keep.cols <- which(!colnames(p_case) %in% levs)
    out <- cbind(out, p_case[,keep.cols, with = FALSE])
  }
  # split up, but only across facets (aka, split up the set run by pop and the set run by fam, but not within either)
  if("by_facet" %in% case){
    p$ord <- 1:nrow(p)
    p_case <- p[,do_correction(.SD, methods), by = "facet"]
    p_case <- p_case[order(p_case$ord),]
    p_case$ord <- NULL
    colnames(p_case)[(ncol(p_case) - length(methods) + 1):ncol(p_case)] <- 
      paste0(pcol, "_byfacet_", colnames(p_case)[(ncol(p_case) - length(methods) + 1):ncol(p_case)])
    keep.cols <- which(!colnames(p_case) %in% "facet")
    out <- cbind(out, p_case[,keep.cols, with = FALSE])
  }
  
  #return
  return(out)
}

# Tiny internal used when converting from a phased, transposed dataset (like
# from an MS output)
# @param x input matrix
.convert_2_to_1_column <- function(x){
  if(!is.matrix(x)){x <- as.matrix(x)}
  ind.genos <- x[,seq(1,ncol(x), by = 2)] + x[,seq(2,ncol(x), by = 2)]
  ind.genos <- matrix(ind.genos, ncol = ncol(x)/2) # rematrix and transpose!
  return(ind.genos)
}

# Perform a row-specific gsub/match given a key table
# 
# @param x input data
# @param lookup table of matches. Each column is something that can be gsub-ed
#   and each row corresponds to the rows of x
# @param lookup_match a character vector with length equal to the number of columns
#   of lookup. In order, the matches in x for the contents of each column
#   in lookup.
# @param patch A character string that will be used to join row ID to contents.
#   SHOULD NOT BE FOUND ANYWHERE IN X, LOOKUP, OR LOOKUP_MATCH.
.row_specific_gsub <- function(x, lookup, lookup_match, patch = "_"){
  ident <- NULL
  
  # make the lookup table row-specific
  lookup <- data.table::data.table(lookup = as.character(unlist(lookup)), row = 1:nrow(lookup), 
                       match = rep(lookup_match, each = nrow(lookup)))
  lookup[,key := paste0(row, patch, match)]
  
  # make x row-specific
  x_rows <- nrow(x)
  x_cols <- ncol(x)
  x <- data.table::data.table(ident = as.character(unlist(x)), row = 1:nrow(x))
  x[,key := paste0(row, patch, ident)]
  x$match <- lookup$lookup[match(x$key, lookup$key)]
  
  return(matrix(lookup$lookup[match(x$key, lookup$key)], x_rows, x_cols))
}


#  Calculate weighted averages of previously genetic statistics. Internal.
# 
#  Calculates a weighted average for a statistic, weighted by the number of
#  called genotypes at each locus. Works for single or pairwise statistics (pi,
#  ho, fst, etc.). Automatically calculates weighted statistic for every
#  previously calculated statistic.
# 
#  Weights are calculated using the equation \deqn{ M_{s} = \frac{\sum_{i =
#  1}^{n} s_{i} * w_{i}}{\sum_{i = 1}^{n} w_{i}}} Where\eqn{n} is the number of
#  SNPs, \eqn{s_{i}} is the value of the statistic in SNP \eqn{i}, and
#  \eqn{w_{i}} is the number of times that SNP was genotyped. Note that this
#  will correct for a range in sequencing depth within samples, but does not
#  really correct for variance in sequencing depth between populations or other
#  facet levels.
# 
#  @param x snpRdata object.
#  @param facets character, default NULL. Facets for which to calculate weighted
#    stats (broken down by category). See \code{\link{Facets_in_snpR}} for
#    details.
#  @param type character, default "single". Type of statistic to weight:
#    \itemize{\item{single: } Statistics calculated in a single subfacet, such
#    as pi. \item{pairwise: } Statistics calculated pairwise between subfacets,
#    such as Fst. }
#  @param stats_to_get character vector, names of the statistics to fetch 
#  weighted averages for.
# 
#  @return A snpR data object with weighted statistics merged in.
.calc_weighted_stats <- function(x, facets = NULL, type = "single", stats_to_get){
  ..drop_col <- ..new.ord <- snp.subfacet <- ..split.snp.part <- snp.facet <- subfacet <- facet <- ..good.cols <- weights_col <- NULL


  #===========sanity checks===============
  msg <- character(0)
  if(!is.snpRdata(x)){
    msg <- c(msg, "x is not a snpRdata object.\n")
  }
  
  if(!type %in%  c("pairwise", "single", "single.window", "sample")){
    msg <- c(msg, "Type unsupported. Options: pairwise, single, single.window, sample.\n")
  }
  if(length(msg) > 0){
    stop(msg)
  }
  #===========subfunc=======================
  clean <- function(x){
    bads <- is.na(x)
    bads[is.infinite(x)] <- TRUE
    bads[is.nan(x)] <- TRUE
    ret <- 1:length(x)
    if(any(bads)){
      ret <- ret[-which(bads)]
    }
    return(ret)
  }
  
  
  #===========calculate weighted stats======
  # setup
  facets <- .check.snpR.facet.request(x, facets, "none")
  x <- .add.facets.snpR.data(x, facets)
  calced <- .check_calced_stats(x, facets, "maf")
  calcedb <- unlist(calced)
  names(calcedb) <- names(calced)
  if(sum(calcedb) != length(calcedb)){
    x <- calc_maf(x, facets = names(calcedb)[which(!calcedb)])
  }
  
  stats <- .get.snpR.stats(x, facets, type)
  
  facets <- .check.snpR.facet.request(x, facets, "none", T)
  if(any(facets[[2]] == "complex") & type %in% c("single.window")){
    facets[[1]] <- c(facets[[1]], facets[[1]][which(facets[[2]] == "complex")])
    facets[[2]] <- c(facets[[2]], rep("special", sum(facets[[2]] == "complex")))
  }
  if(any(facets[[2]] == "snp") & type %in% c("single.window")){
    facets[[1]] <- c(facets[[1]], facets[[1]][which(facets[[2]] == "snp")])
    facets[[2]] <- c(facets[[2]], rep("regress.base", sum(facets[[2]] == "snp")))
  }
  
  # calc
  for(i in 1:length(facets[[1]])){

    split.part <- unlist(.split.facet(facets[[1]][i]))
    split.part <- .check.snpR.facet.request(x, split.part, remove.type = "none", TRUE)
    snp.part <- split.part[[1]][which(split.part[[2]] == "snp")]
    split.snp.part <- unlist(.split.facet(split.part[[1]][which(split.part[[2]] == "snp")]))
    samp.part <- split.part[[1]][which(split.part[[2]] == "sample")]
    split.samp.part <- unlist(.split.facet(split.part[[1]][which(split.part[[2]] == "sample")]))
    
    
    
    # get weights and stats
    #=======================simple case==================
    if(type == "single"){
      weights <- (stats$maj.count + stats$min.count)/2
      
      if(facets[[2]][i] == "complex"){
        split.part <- unlist(.split.facet(facets[[1]][i]))
        split.part <- .check.snpR.facet.request(x, split.part, remove.type = "none", TRUE)
        if(length(split.part[[1]]) > 2){
          snp.partp <- paste(snp.part, collapse = ".")
          samp.partp <- paste(samp.part, collapse = ".")
          split.part <- list(c(snp.partp, samp.partp), c("snp", "sample"))
        }
        keep.rows <- which(stats$facet %in% split.part[[1]][which(split.part[[2]] == "sample")])
        
        group_key <- c("facet", "subfacet", split.snp.part)
      }
      else if(facets[[2]][i] == "sample"){
        keep.rows <- which(stats$facet %in% facets[[1]][i])
        group_key <- c("facet", "subfacet")
      }
      else if(facets[[2]][i] == "snp"){
        keep.rows <- which(stats$facet == ".base")
        group_key <- split.snp.part
      }
      else{
        keep.rows <- which(stats$facet == ".base")
        group_key <- "facet"
      }
      
      
      selected_stats <- stats[keep.rows, stats_to_get, drop = F]
      if(ncol(selected_stats) == 0){
        stop("No calculated stats to weight.\n")
      }
      weights <- weights[keep.rows]
      
      group_key_tab <- stats[keep.rows,group_key, drop = F]
      group_key_tab$key <- do.call(paste, group_key_tab)
    }
    #=====================single window case======================
    else if(type == "single.window"){
      weights <- stats$n_snps
      
      
      if(facets[[2]][i] == "complex" | facets[[2]][i] == "special"){
        if(length(split.part[[1]]) > 2){
          snp.partp <- paste(snp.part, collapse = ".")
          samp.partp <- paste(samp.part, collapse = ".")
          split.part <- list(c(snp.partp, samp.partp), c("snp", "sample"))
        }
        keep.rows <- which(stats$facet %in% split.part[[1]][which(split.part[[2]] == "sample")] &
                             stats$snp.facet %in% split.part[[1]][which(split.part[[2]] == "snp")])
        
        if(facets[[2]][i] == "complex"){
          group_key <- c("facet", "subfacet", "snp.facet", "snp.subfacet")
        }
        else{
          group_key <- c("facet", "subfacet")
        }
      }
      else if(facets[[2]][i] == "sample"){
        keep.rows <- which(stats$facet %in% facets[[1]][i] & stats$snp.facet == ".base")
        group_key <- c("facet", "subfacet")
      }
      else if(facets[[2]][i] == "snp"){
        keep.rows <- which(stats$facet == ".base" & stats$snp.facet == facets[[1]][i])
        group_key <- c("snp.facet", "snp.subfacet")
      }
      else if(facets[[2]][i] == "regress.base"){
        keep.rows <- which(stats$facet == ".base" & stats$snp.facet == facets[[1]][i])
        group_key <- "facet"
      }
      else{
        keep.rows <- which(stats$facet == ".base" & stats$snp.facet == ".base")
        group_key <- "facet"
      }
      
      
      selected_stats <- stats[keep.rows, stats_to_get, drop = F]
      if(ncol(selected_stats) == 0){
        stop("No calculated stats to weight.\n")
      }
      weights <- weights[keep.rows]
      
      group_key_tab <- stats[keep.rows,group_key, drop = F]
      group_key_tab$key <- do.call(paste, group_key_tab)
    }
    #=======================pairwise case==========================
    else if(type == "pairwise"){
      weights <- stats$nk
      
      
      if(facets[[2]][i] == "complex"){
        if(length(split.part[[1]]) > 2){
          snp.partp <- paste(snp.part, collapse = ".")
          samp.partp <- paste(samp.part, collapse = ".")
          split.part <- list(c(snp.partp, samp.partp), c("snp", "sample"))
        }
        keep.rows <- which(stats$facet %in% split.part[[1]][which(split.part[[2]] == "sample")])
        
        group_key <- c("facet", "comparison", split.part[[1]][which(split.part[[2]] == "snp")])
      }
      else if(facets[[2]][i] == "sample"){
        keep.rows <- which(stats$facet %in% facets[[1]][i])
        group_key <- c("facet", "comparison")
      }
      else{
        stop("Cannot calculate pairwise stats for non-sample facets--check what is being passed as facets to calc_weighted_averages.")
      }
      
      
      selected_stats <- stats[keep.rows, stats_to_get, drop = F]
      if(ncol(selected_stats) == 0){
        stop("No calculated stats to weight.\n")
      }
      weights <- weights[keep.rows]
      
      group_key_tab <- stats[keep.rows,group_key, drop = F]
      group_key_tab$key <- do.call(paste, group_key_tab)
      
    }
    #=======================sample case========================
    else if(type == "sample"){
      if(i != 1){
        if(facets[[2]][i - 1] == "complex"){
          stats <- ostats
          rm(ostats); gc()
        }
      }
      # get the per sample missingness
      if(facets[[2]][i] == "sample"){
        # per sample, weights are easy
        weights <- nrow(x) - matrixStats::colSums2(genotypes(x) == x@mDat)
        keep.rows <- which(stats$snp.facet == ".base")
        
        group_key <- split.samp.part
        group_key_tab <- stats[keep.rows, group_key, drop = F]
        group_key_tab$key <- do.call(paste, c(group_key_tab, sep = "."))
        colnames(group_key_tab)[1] <- "subfacet"
        group_key_tab$subfacet <- group_key_tab$key
        group_key_tab$facet <- paste0(samp.part, collapse = ".")
      }
      if(facets[[2]][i] == "snp"){
        keep.rows <- which(stats$snp.facet == snp.part & stats$facet == ".base")
        group_key <- "snp.subfacet"
        group_key_tab <- stats[keep.rows, group_key, drop = F]
        group_key_tab$key <- group_key_tab$snp.subfacet
        
        
        tl <- .get.task.list(x, facets[[1]][i])
        weights <- numeric(nrow(group_key_tab))
        for(j in 1:nrow(tl)){
          snps_in_set <- which(unlist(do.call(paste, c(snp.meta(x)[,split.snp.part, drop = F], sep = "."))) ==
                                 tl[j,4])
          tweights <- length(snps_in_set) - matrixStats::colSums2(genotypes(x)[snps_in_set,] == x@mDat)
          weights[which(group_key_tab$key == tl[j,4])] <- tweights
        }
        group_key_tab$snp.facet <- snp.part
      }
      if(facets[[2]][i] == "complex"){
        make_group_key_tab <- function(stats, keep.rows, group_key){
          group_key_tab <- stats[keep.rows, group_key, drop = F]
          group_key_tab$key <- unlist(do.call(paste, c(group_key_tab, sep = ".")))
          group_key_tab$snp.facet <- snp.part
          group_key_tab$facet <- samp.part
          group_key_tab$subfacet <- unlist(do.call(paste, c(group_key_tab[,split.samp.part, drop = F], sep = ".")))
          return(group_key_tab)
        }
        group_key <- c("snp.subfacet", split.samp.part)
        keep.rows <- which(stats$snp.facet == snp.part & stats$facet == ".base")
        group_key_tab <- make_group_key_tab(stats, keep.rows, group_key)
        
        
        
        tl <- .get.task.list(x, facets[[1]][i])
        weights <- stats[keep.rows,]
        for(j in 1:nrow(tl)){
          
          snps_in_set <- which(unlist(do.call(paste, c(snp.meta(x)[,split.snp.part, drop = F], sep = "."))) ==
                                 tl[j,4])
          samps_in_set <- .fetch.sample.meta.matching.task.list(x, tl[j,])
          if(length(snps_in_set) > 0 & length(samps_in_set) > 0){
            tweights <- length(snps_in_set) - matrixStats::colSums2(genotypes(x)[snps_in_set, samps_in_set, drop = F] == x@mDat)
            m <- data.table(weights = tweights, facet = tl[j,1], subfacet = tl[j,2], snp.facet = tl[j,3], snp.subfacet = tl[j,4])
            m <- cbind(m, sample.meta(x)[samps_in_set,])
            mn <- c("snp.subfacet", "snp.subfacet", colnames(x@sample.meta))
            weights <- .smart.merge(m, weights, mn, mn)
          }
        }
        
        ostats <- stats
        stats <- weights
        keep.rows <- 1:nrow(stats)
        weights <- weights$weights
        stats <- as.data.frame(stats)
        group_key_tab <- make_group_key_tab(stats, keep.rows, group_key)
      }
      if(facets[[2]][i] == ".base"){
        weights <- nrow(x) - matrixStats::colSums2(genotypes(x) == x@mDat)
        keep.rows <- which(stats$snp.facet == ".base" & stats$facet == ".base")
        group_key <- c("facet", "subfacet", "snp.facet", "snp.subfacet")
        group_key_tab <- stats[keep.rows, group_key, drop = F]
        group_key_tab$key <- do.call(paste, c(group_key_tab, sep = "."))
      }
      
      selected_stats <- stats[keep.rows, stats_to_get, drop = F]
    }
    else{stop(".calc_weighted_stats type not recognized.")}
    
    #===================calculate weighted means using equation sum(w*s)/sum(w)=======================
    # override the group_key_tab if requested
    weighted <- data.table::as.data.table(selected_stats)
    weighted$weights_col <- weights
    weighted$key <- group_key_tab$key
    means <- weighted[,.SDcols = stats_to_get, lapply(.SD, function(x, w) stats::weighted.mean(x[clean(x)], w[clean(x)]), w = weights_col), by = key]
    
    # merge and return
    ## get the stats back in a format with facet and sub-facet, clean up, and return
    mstats <- merge(means, unique(group_key_tab), by = "key")
    drop_col <- which(colnames(mstats) == "key")
    mstats <- .fix..call(mstats[,-..drop_col])
    new.ord <- c((length(stats_to_get) + 1):(ncol(selected_stats) + ncol(group_key_tab) - 1), 1:ncol(selected_stats))
    mstats <- .fix..call(mstats[,..new.ord])
    colnames(mstats)[(ncol(group_key_tab)):ncol(mstats)] <- paste0("weighted_mean_", colnames(mstats)[(ncol(group_key_tab)):ncol(mstats)])
    
    
    
    #====================add on extra filler columns====================
    if(!"snp.subfacet" %in% colnames(mstats)){
      if(facets[[2]][i] %in% c(".base", "sample")){
        mstats <- mstats[,snp.subfacet := ".base"]
      }
      else{
        if(facets[[2]][i] == "special"){
          mstats$snp.subfacet <- ".OVERALL_MEAN"
        }
        else if(facets[[2]][i] == "regress.base"){
          mstats$subfacet <- ".base"
          mstats$snp.subfacet <- ".OVERALL_MEAN"
          mstats$snp.facet <- facets[[1]][i]
        }
        else if(length(split.snp.part) > 1){
          mstats$snp.subfacet <- do.call(paste, c(.fix..call(mstats[,..split.snp.part]), sep = "."))
        }
        else{
          mstats$snp.subfacet <- .fix..call(mstats[,..split.snp.part])
        }
      }
    }
    if(!"snp.facet" %in% colnames(mstats)){
      if(facets[[2]][i] %in%  "special"){
        snp.partp <- paste(snp.part, collapse = ".")
        mstats <- mstats[,snp.facet := snp.partp]
      }
      else if(facets[[2]][i] %in% c(".base", "sample")){
        mstats <- mstats[,snp.facet := ".base"]
      }
      else{
        mstats <- mstats[,snp.facet := .check.snpR.facet.request(x, paste0(split.snp.part, collapse = "."), remove.type = "sample")]
      }
    }
    if(!"subfacet" %in% colnames(mstats)){
      if("comparison" %in% colnames(mstats)){
        colnames(mstats)[which(colnames(mstats) == "comparison")] <- "subfacet"
      }
      else{
        mstats <- mstats[,subfacet := ".base"]
      }
    }
    if(!"facet" %in% colnames(mstats)){
      mstats <- mstats[,facet := ".base"]
    }
    
    good.cols <- c("facet", "subfacet", "snp.facet", "snp.subfacet", paste0("weighted_mean_", stats_to_get))
    mstats <- .fix..call(mstats[,..good.cols])
    mstats[,1:4] <- dplyr::mutate_all(mstats[,1:4], as.character)
    x <- .merge.snpR.stats(x, mstats, type = "weighted.means")
  }
  
  
  
  return(x)
}

# Split a single or multiple compound facet into parts
# 
# @param facet compound facet(s) to split
.split.facet <- function(facet) strsplit(facet, "(?<!^)\\.", perl = T)

# Paste together metadata according to a list of column names
# 
# @param df data.frame with data do paste
# @param facets facets to paste together. Often produced by \code{\link{.split.facet}}. Can also be a numeric vector of columns to use.
# @param sep character, default ".". Pasted facets will be split by this.
.paste.by.facet <- function(df, facets, sep = "."){
  ..facets <- NULL
  if(is.data.table(df)){
    return(do.call(paste, c(.fix..call(df[,..facets, drop = FALSE]), sep = sep)))
  }
  else{
    return(do.call(paste, c(df[,facets, drop = FALSE], sep = sep)))
  }
}

# Fixes calling scope warning in .. calls with data.table
# 
# Needed since ..x vars need to be defined as global variables for CRAN, which data.table will throw a warning about.
# 
# @param fun function to run
.fix..call <- function(fun){
  
  
  
  return(.suppress_specific_warning(fun, "variable in calling scope for clarity"))
}



# dist subfunction for pops. Adapted from adegenet code!
# @param x matrix, genotypes or other character states
# @param method character, default Edwards. Dist method.
.get_dist <- function(x, method = "Edwards"){
  #browser()
  if(method == "Edwards"){
    x <- x[,which(colSums(is.na(x)) == 0)] # remove anywhere where there is missing data!
    nloc <- ncol(x)
    x <- sqrt(as.matrix(x))
    am <- x%*%t(x)
    am <- 1 - (am / (nloc/2)) # can produce negative values, meaning an NaN distance measure
    diag(am) <- 0
    suppressWarnings(am <- sqrt(am))
    am <- stats::as.dist(am)
  }
  else if(method == "Nei"){
    d <- x %*% t(x)
    vec <- sqrt(diag(d))
    d <- d/vec[col(d)]
    d <- d/vec[row(d)]
    d <- -log(d)
    am <- stats::as.dist(d)
  }
  am <- list(am)
  names(am) <- method
  return(am)
}


# Determine minor and major alleles and get allele counts from @geno.tables
# 
# @param gs @geno.tables part of a snpRdata object
# @param m.al missing allele indicator
# @param ref If not the .base facet, a two column data.frame with column names
#   "major" and "minor" containing the allele idendities ("A", "C", etc.).
#
# @return A vector of the filenames for the new datasets.
.maf_func <- function(gs, m.al, ref = NULL){

  maj_count <- NULL
  
  # for the base facet, determine the major and minor then calculate maf
  if(is.null(ref)){

    # major alleles via max.col
    fmax <- colnames(gs$as)[max.col(gs$as, ties.method = "last")]
    lmax <- colnames(gs$as)[max.col(gs$as, ties.method = "first")]
    
    # minor alleles, via inversion, 0 replacement with -Inf, and max.col
    inv.as <- gs$as * -1
    inv.as[inv.as == 0] <- -Inf
    fmin <- colnames(gs$as)[max.col(inv.as, ties.method = "last")]
    lmin <- colnames(gs$as)[max.col(inv.as, ties.method = "first")]
    
    # special cases
    match.freq <- which(fmax != lmax) # maf = 0.5
    unseq <- which(Matrix::rowSums(gs$as) == 0) # unsequenced
    np <- which(Matrix::rowSums(methods::as(gs$as, "lMatrix")) == 1) # non-polymorphic
    

    # declare major and minor
    major <- fmax
    minor <- fmin
    ## maf = 0.05
    if(length(match.freq) != 0){
      minor[match.freq] <- lmax[match.freq]
    }
    ## unsequenced
    if(length(unseq) != 0){
      major[unseq] <- "N"
      minor[unseq] <- "N"
    }
    ## non-polymorphic
    if(length(np) != 0){
      minor[np] <- "N"
    }
    
    # grab the actual maf

    maf <- 1 - .rowMax_sparse(gs$as)/Matrix::rowSums(gs$as)
    maf[is.nan(maf)] <- 0
    
    # get the major and minor counts
    # round because repeating decimals will yeild things like 1.00000000001 instead of 1. Otherwise this approach is quick and easy, as long as things are bi-allelic (non-polymorphic and equal min maj frequencies are fine.)
    maj.count <- round(Matrix::rowSums(gs$as)*(1-maf))
    min.count <- round(Matrix::rowSums(gs$as)*(maf))
  }
  
  # for non-base facets, use the given major and minor to calculate maf
  else{
    # use a data.table function to get the major allele counts
    adt <- as.data.table(as.matrix(gs$as))
    rep.factor <- nrow(gs$as)/nrow(ref)
    major <- rep(ref$major, each = rep.factor) # rep for each facet level, since that's how they are sorted
    adt$major <-  major
    adt <- adt[, maj_count := .SD[[.BY[[1]]]], by=major] # get the counts of the major allele in each column
    ac.rows <- 1:(which(colnames(adt) == "major") - 1)
    total.allele.count <- rowSums(adt[,ac.rows, with = FALSE])
    maf <- 1 - adt$maj_count/total.allele.count # get the maf
    
    # other things for return
    minor <- rep(ref$minor, each = rep.factor)
    maj.count <- adt$maj_count
    min.count <- total.allele.count - maj.count
  }
  
  # return
  return(data.table(major = major, minor = minor, maj.count = maj.count, min.count = min.count, maf = maf, stringsAsFactors = F))
}


# Determine observed heterozygosity from @geno.tables
# 
# @param gs @geno.tables part of a snpRdata object
# 
# @return a vector of observed heterozygosity
.ho_func <- function(gs, snp_form){
  #identify heterozygote rows in genotype matrix
  genos <- colnames(gs$gs)
  hets <- which(substr(genos, 1, snp_form/2) != substr(genos, (snp_form/2) + 1, snp_form))
  
  # calculate ho
  ## if only one heterozygote...
  if(length(hets) == 1){
    ho <- gs$gs[,hets]/Matrix::rowSums(gs$gs)
  }
  ## if no heterozygotes
  else if(length(hets) == 0){
    ho <- rep(0, nrow(gs$gs))
  }
  ## normally
  else{
    ho <- Matrix::rowSums(gs$gs[,hets])/Matrix::rowSums(gs$gs)
  }
}

# Generate random bootstraps of a genepop input file
#
# New files will have random names ending in .genepop in the working directory
# 
# @param gp_filepath Character, path to file to permute
# @param n numeric, number of bootstrapped datasets to generate.
#
# @return A vector of the filenames for the new datasets.
.boot_genepop <- function(gp_filepath, n){
  ..shuff <- NULL
  
  #=========parse gp file=========
  
  .suppress_specific_warning(gp_file <- readLines(gp_filepath), "incomplete final line found on") # need to suppress warnings here because genepop::Fst will remove line endings from the input file for some reason......
  pop_headers <- grep("POP", gp_file)
  
  counts <- c(pop_headers[-1], length(gp_file)) - pop_headers
  header <- gp_file[1:(pop_headers[1] - 1)]
  gp_file <- gp_file[-c(1:(pop_headers[1] - 1), pop_headers)]
  gp_file <- strsplit(gp_file, "\t")
  gp_file <- data.frame(gp_file)
  gp_file <- t(gp_file)
  rns <- gp_file[,1]
  gp_file <- gp_file[,-1]
  gp_file <- as.data.frame(gp_file)
  row.names(gp_file) <- rns
  
  #=========boot=============
  rstrings <- .rand_strings(n, 10)
  rstrings <- paste0("./", rstrings)
  exists <- dir.exists(rstrings)
  abort <- 0
  while(sum(exists) > 0){
    rstrings[exists] <- .rand_strings(sum(exists), 10)
    rstrings <- paste0("./", rstrings)
    exists <- dir.exists(rstrings)
    abort <- abort + 1
    
    if(abort > 10000){
      stop("Cannot generate random directory names for bootstraps---there are (26+26+10)^10 possible names, so this shouldn't generally happen...\n")
    }
  }
  lapply(rstrings, dir.create)
  rstrings_files <- paste0(rstrings, "/genepop.txt")
  
  
  adj_pop_headers <- pop_headers - (pop_headers[1] - 1)
  adj_pop_headers <- adj_pop_headers - 1:length(adj_pop_headers)
  adj_pop_headers <- c(adj_pop_headers, nrow(gp_file))
  for(i in 1:n){
    shuff <- sample(nrow(gp_file), nrow(gp_file), replace = T)
    tgp <- gp_file[shuff,]
    rownames(tgp) <- rownames(gp_file)
    writeLines(header, rstrings_files[i], sep = "\n")
    for(j in 1:length(pop_headers)){
      write("POP", rstrings_files[i], append = TRUE)
      data.table::fwrite(tgp[(adj_pop_headers[j]+1):adj_pop_headers[j + 1],], 
            file = rstrings_files[i],
            append = TRUE, 
            quote = FALSE, 
            row.names = TRUE, 
            col.names = FALSE, 
            sep = "\t")
    }
  }
  
  return(rstrings)
}


# Generate random bootstraps in ac format
# 
# @param x snpRdata object to permute
# @param n numeric, number of bootstrapped datasets to generate.
# @param facet sample level facet to permute within
#
# @return A list containing permuted ac data.
.boot_ac <- function(x, n, facet){
  ..tm <- NULL

  # function to convert the output of .maf_func into ac
  maf.to.ac <- function(maf){
    ac <- data.table::data.table(n_total = maf$maj.count + maf$min.count,
                                 n_alleles = rowSums(maf[,c("maj.count", "min.count")] != 0),
                                 ni1 = maf$maj.count,
                                 ni2 = maf$min.count,
                                 ho = maf$ho)
    return(ac)
  }
  
  out <- vector("list", n)
  opts <- .get.task.list(x, facet)
  for(i in 1:n){
    
    # get the maf identities for the base facet
    shuff <- genotypes(x)[,sample(ncol(x), ncol(x), replace = TRUE)]
    shuff <- data.table::as.data.table(shuff)
    colnames(shuff) <- colnames(genotypes(x))
    tac <- vector("list", nrow(opts))
    glob_tab <- .tabulate_genotypes(shuff, x@mDat)
    glob <- .maf_func(glob_tab, m.al = substr(x@mDat, 1, floor(nchar(x@mDat)/2)))
    glob$ho <- .ho_func(glob_tab, x@snp.form)

    # if bootstrapping the base level, done
    if(facet == ".base"){
      out[[i]] <- maf.to.ac(glob)
      out[[i]]$facet <- ".base"
      out[[i]]$subfacet <- ".base"
      out[[i]]$.snp.id <- 1:length(shuff)
    }
    # otherwise need to do for each facet level.
    else{
      for(j in 1:nrow(opts)){
        tm <- .fetch.sample.meta.matching.task.list(x, opts[j,])
        ttab <- .tabulate_genotypes(.fix..call(shuff[,..tm]), x@mDat)
        tac[[j]] <- .maf_func(ttab, x@mDat, as.data.frame(glob[,c("major", "minor")]))
        tac[[j]]$ho <- .ho_func(ttab, x@snp.form)
        tac[[j]] <- maf.to.ac(tac[[j]])
        tac[[j]]$.snp.id <- x@snp.meta$.snp.id
        tac[[j]]$subfacet <- opts[j,2]
        tac[[j]]$facet <- facet
      }
      tac <- data.table::rbindlist(tac)
      out[[i]] <- tac
    }
  }
  
  return(out)
}

# Generate random bootstraps in as format
# 
# @param x snpRdata object to permute
# @param n numeric, number of bootstrapped datasets to generate.
# @param facet sample level facet to permute within
# @param by char, default "sample". Permute by samples or snps
#
# @return A list containing permuted ac data.
.boot_as <- function(x, n, facet, ret_gs = FALSE, by = "sample"){
  ..tm <- ..ord <-  NULL
  
  out <- vector("list", n)
  opts <- .get.task.list(x, facet)
  for(i in 1:n){
    
    # get the maf identities for the base facet
    if(by == "sample"){
      shuff <- genotypes(x)[,sample(ncol(x), ncol(x), replace = TRUE)]
    }
    else{
      shuff <- genotypes(x)[sample(nrow(x), nrow(x), replace = TRUE),]
    }
    shuff <- data.table::as.data.table(shuff)
    colnames(shuff) <- colnames(genotypes(x))
    tas <- vector("list", nrow(opts))
    
    # do for each facet level.
    for(j in 1:nrow(opts)){
      tm <- .fetch.sample.meta.matching.task.list(x, opts[j,])
      ttab <- .tabulate_genotypes(.fix..call(shuff[,..tm]), x@mDat)
      tas[[j]] <- data.table::as.data.table(as.matrix(ttab$as))
      
      # fix missing columns (no genotypes with this allele in this subfacet)
      missing <- which(!colnames(x@geno.tables$as) %in% colnames(tas[[j]]))
      if(length(missing) != 0){
        fill <- matrix(0, nrow(tas[[j]]), length(missing))
        fill <- as.data.table(fill)
        colnames(fill) <- colnames(x@geno.tables$as)[missing]
        tas[[j]] <- cbind(tas[[j]], fill)
      }
      ord <- colnames(x@geno.tables$as)
      tas[[j]] <- .fix..call(tas[[j]][,..ord])
      
      if(ret_gs){
        gst <- as.data.table(as.matrix(ttab$gs))

        # fix missing columns (no genotypes with this allele in this subfacet)
        missing <- which(!colnames(x@geno.tables$gs) %in% colnames(gst))
        if(length(missing) != 0){
          fill <- matrix(0, nrow(gst), length(missing))
          fill <- as.data.table(fill)
          colnames(fill) <- colnames(x@geno.tables$gs)[missing]
          gst <- cbind(gst, fill)
        }
        ord <- colnames(x@geno.tables$gs)
        gst<- .fix..call(gst[,..ord])
        
        tas[[j]] <- cbind(tas[[j]], gst)
      }
      
      tas[[j]]$ho <- .ho_func(ttab, x@snp.form)
      tas[[j]]$.snp.id <- x@snp.meta$.snp.id
      tas[[j]]$subfacet <- opts[j,2]
      tas[[j]]$facet <- facet
      
    }
    tas <- data.table::rbindlist(tas)
    tas[is.na(tas)] <- 0
    ord <- c((ncol(tas) - 3):ncol(tas), 1:(ncol(tas) - 4))
    out[[i]] <- .fix..call(tas[,..ord])
    out[[i]] <- .suppress_specific_warning(cbind(snp.meta(x)[,which(colnames(snp.meta(x)) != ".snp.id"),drop=FALSE], out[[i]]), "row names were found from a short variable")
  }
  
  return(out)
}


# variance components for fst, etc
# @param intot number of genotypes in pop 1
# @param jntot number of genotypes in pop 2
# @param ps1 allele frequency, allele 1
# @param ps2 allele frequency, allele 2
# @param r number of comparisons
# @param nbar average sample size across all populations
# @param nc If CV = coefficient of variation in sample size, nc = nbar*(1-(CV^2)/r)
# @param iho ho for pop 1
# @param jho ho for pop 2
.per_all_f_stat_components <- function(intot = NULL, jntot = NULL, ps1 = NULL, ps2 = NULL, r, nbar, nc, iho = NULL, jho = NULL, ntotm = NULL, psm = NULL, hom = NULL){
  if(r == 2 & !is.null(intot)){
    pbar <- ((intot*ps1) + (jntot*ps2))/(r*nbar) #average sample allele frequency
    ssq <- (((intot)*(ps1-pbar)^2) + ((jntot)*(ps2-pbar)^2))/((r-1)*nbar) #sample variance of allele frequencies
    hbar <- ((intot*iho) + (jntot*jho))/(r*nbar) #average heterozygote frequencies
  }
  else if (r == 1){
    pbar <- ps1 #average sample allele frequency
    ssq <- 0 #sample variance of allele frequencies
    hbar <- iho #average heterozygote frequencies
  }
  else if(r > 2 | (r == 2 & is.null(intot))){
    pbar <- rowSums(ntotm * psm)/(r*nbar) #average sample allele frequency
    ssq <- rowSums(ntotm*(psm-pbar)^2)/((r - 1)*nbar)
    hbar <- rowSums(ntotm*hom)/(r*nbar)
  }
  
  
  
  #equation parts used in both
  inner1 <- pbar*(1-pbar)
  inner2 <- ((r-1)/r)*ssq
  inner3 <- .25*hbar
  
  inner4 <- ((2*nbar - 1)/(4*nbar))*hbar
  a <- (nbar/nc) * (ssq - (1/(nbar - 1))*(inner1 - inner2 - inner3))
  b <- (nbar/(nbar-1))*(inner1 - inner2 - inner4)
  c <- .5*hbar
  
  # # browser()
  # if(.print & !is.null(intot)){
  #   saveRDS(data.frame(pbar = pbar, ssq = ssq, hbar = hbar, iho = iho, intot = intot, jntot = jntot, jho = jho, inner1 = inner1, inner2 = inner2, inner3 = inner3, inner4 = inner4, a = a, b = b, c = c), "test.RDS")
  # }
  # else if(.print){browser()}
  return(list(a = a, b = b, c = c))
  
  # weir--exactly the same
  # else{
  #   S1 <- ssq - (1/(nbar-1))*(inner1 - inner2 - inner3)
  #   S2i1 <- ((r*(nbar - nc))/nbar)*inner1
  #   S2i2 <- (1/nbar)*((nbar-1)+(r-1)*(nbar-nc))*ssq
  #   S2i3 <- ((nbar-nc)/(4*nc^2))*hbar
  #   S2 <- inner1 - (nbar/(r*(nbar-1)))*(S2i1 -S2i2 - S2i3)
  #   return(list(S1 = S1, S2 = S2))
  # }
  
  
}

# simple function to update citations
# @param x snpRdata object
# @param keys bibtex keys to use, corresponding to snpr_citations.bib in extdata
# @param stats vector of stat names
# @param details details on how each key was used 
.update_citations <- function(x, keys, stats, details){
  new.citations <- vector("list", length(keys))
  names(new.citations) <- stats
  for(i in 1:length(new.citations)){
    new.citations[[i]] <- list(key = keys[i], details = details[i])
  }
  x@citations <- c(x@citations, new.citations)
  x@citations <- x@citations[!duplicated(x@citations)]
  return(x)
}

# simple function to yell citation info and add to an external .bib
#
# For when a snpRdata object is not returned
# @param keys bibtex keys to use, corresponding to snpr_citations.bib in extdata
# @param stats vector of stat names
# @param details details on how each key was used
# @param bib_path FALSE or path to bib file to update
.yell_citation <- function(keys, stats, details, outbib = FALSE){
  #==========sanity checks=======
  
  .check.installed("rbibutils")

  #==========shout at the user=====
  bib.file <- system.file("extdata", "snpR_citations.bib", package = "snpR")
  bib <- rbibutils::readBib(bib.file)

  cat("Citations for methods used thus far: \n")
  for(i in 1:length(keys)){
    cat("==============================================\n")
    cat("Statistic: ", stats[i], "\n")
    cat("Citation: ", .quick_grab_cite(bib[keys[i]]), "\n")
    cat("Bibtex key: ", keys[i], "\n")
    cat("Details: ", details[i], "\n")
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
        
        cat(".bib file ", outbib, "updated with new citations.\n")
        
      }
      else{
        cat(".bib file ", outbib, "not updated, no new citations.\n")
      }
    }
    
    else {
      bib <- bib[keys]
      rbibutils::writeBib(bib, outbib)
      
      cat(".bib file ", outbib, "written.\n")
    }
  }
}

# update filters
#
# Note that filters are stored as a data.table instead of a named list by filter 
# type because the order they were applied matters!
#
# @param x snpRdata object
# @param type The type of filter (maf, mac, etc)
# @param stringency The filter stringency
# @param facet Any facets the filter was computed across
.update_filters <- function(x, type, stringency, facet){

  # update old objects
  r <- try(x@filters, silent = TRUE)
  if(methods::is(r, "try-error")){
    .suppress_specific_warning(x@filters <- data.table::data.table(type = character(0),
                                                                   stringency = character(0),
                                                                   facet = character(0)), "Not a validObject")
  }
  else{
    if(ncol(x@filters) == 0){
      r <- data.table::data.table(type = character(0),
                                  stringency = character(0),
                                  facet = character(0))
    }
    
    for(i in 1:length(type)){
      x@filters <- rbind(x@filters,
                         data.table(type = type, stringency = stringency, facet = facet))
    }
  }
  
  return(x)
}

.heterozygosity <- function(x, mDat, method){
  # make x into a logical for heterozygous
  xv <- as.matrix(x)
  logix <- substr(xv, 1, 1) != substr(xv, 2, 2) # true if het
  logix[xv == mDat] <- NA # NA when missing data!
  
  # get counts of hets and homs
  hets <- matrixStats::colSums2(logix, na.rm = T) # number heterozygous sites
  homs <- matrixStats::colSums2(!logix, na.rm = T) # number homozygous sites
  
  # if the raw ratio is desired:
  if(method == "ratio"){
    ratio <- hets/homs
    return(ratio)
  }
  
  # else if hs is desired, more of a pain since we need to do a loop because we omit different loci for each individual sample
  else if(method == "hs"){
    ind_percent <- (hets/(hets + homs))
    pop_means <- numeric(length(hets))
    for(i in 1:length(pop_means)){
      pop_means[i] <- mean(rowMeans(logix[which(!is.na(logix[,i])),, drop = FALSE], na.rm = TRUE))
    }
    hs <- ind_percent/pop_means
    return(hs)
  }
  
}

.suppress_specific_warning <- function(fun, warn_to_suppress){
  warn_handle <- function(warn){
    if(any(grepl(warn_to_suppress, warn))){
      invokeRestart("muffleWarning")
    }
  }
  
  withCallingHandlers(res <- fun,
                      warning = warn_handle)
  
  
  return(res)
}


.rand_strings <- function(n, length){
  chrs <- sample(c(LETTERS, letters, 1:9), n*length, replace = TRUE)
  chrs <- split(chrs, rep(1:n, length.out = n*length))
  chrs <- unlist(lapply(chrs, function(x) paste0(x, collapse = "")))
  chrs <- as.character(chrs)
  return(chrs)
}

.make_it_quiet <- function(fun){
  return(invisible(utils::capture.output(fun)))
}

.quick_grab_cite <- function(cite) paste0(utils::capture.output(cite), collapse = "")

.par_checker <- function(par, ret_1_on_FALSE = FALSE){
  if(!isFALSE(par)){
    if(is.numeric(par)){
      if(par != floor(par)){
        stop("par must be an integer.\n")
      }
      cor_avail <- parallel::detectCores()
      if(par > cor_avail){
        par <- cor_avail
      }
    }
    else{
      stop("par must be a numeric value or FALSE.\n")
    }
  }
  else{
    if(ret_1_on_FALSE){
      par <- 1
    }
  }
  
  return(par)
}

.is.bi_allelic <- function(x){
  r <- try(x@bi_allelic, silent = TRUE)
  if(methods::is(r, "try-error")){
    return(TRUE)
  }
  else{
    return(x@bi_allelic)
  }
}

.ploidy <- function(x){
  r <- try(x@ploidy, silent = TRUE)
  if(methods::is(r, "try-error")){
    return(2)
  }
  else{
    return(x@ploidy)
  }
}

.console_hline <- function(char = "="){
  return(paste0(rep(char, getOption("width")), collapse = ""))
}

.update.sample.stats.with.new.metadata <- function(x, nsm){
  if(nrow(x@sample.stats) == 0){
    return(x)
  }

  # identify columns in the new data also in the old
  osmc <- colnames(x@sample.stats)[which(colnames(x@sample.stats) %in% colnames(sample.meta(x)))]
  
  # identify and remove an removed columns
  rm.cols <- which(!osmc %in% colnames(nsm))
  if(length(rm.cols) > 0){
    t.cols <- osmc[rm.cols]
    x@sample.stats[,t.cols] <- NULL
    osmc <- osmc[-rm.cols]
  }
  
  # update maintained columns
  x@sample.stats[,osmc] <- 
    nsm[match(x@sample.stats$.sample.id, nsm$.sample.id),osmc]
  
  
  # add any new columns
  newcols <- colnames(nsm)[which(!colnames(nsm) %in% osmc)]
  if(length(newcols) > 0){
    x@sample.stats[,newcols] <- nsm[match(x@sample.stats$.sample.id, nsm$.sample.id),newcols]
    ord <- c("facet", "subfacet", "snp.facet", "snp.subfacet", 
             colnames(nsm)[-which(colnames(nsm) == ".sample.id")], ".sample.id")
    ord <- c(ord, colnames(x@sample.stats)[which(!colnames(x@sample.stats) %in% ord)])
  }
  
  return(x)
}

# g: number to rarefact to
.richness_parts <- function(gs, private = TRUE, alleles, g = 0){
  facet <- subfacet <- ..al_cols <- ..p_al_cols <- .snp.id <- .sum <- .g <- . <- .snp.id <- prob_seg <- NULL
  # equations from https://doi.org/10.1023/B:COGE.0000041021.91777.1a
  # Nij is the table
  # Nj is the rowsums
  # m is the row sum of != 0

  # step: need to calculate g for each row
  as <- gs$as
  rm(gs); gc();
  as <- data.table::as.data.table(as)
  al_cols <- alleles
  as[,.sum := rowSums(.SD), .SDcols = al_cols] # sums for each row
  
  if(g == 0){
    as[,.g := min(.sum), by = .(facet, .snp.id)] # min across all levels
  }
  else if(g < 0){
    as[,.g := min(.sum) + g, by = .(facet, .snp.id)] # min across all levels - g
  }
  else if(g > 0){
    as[,.g := g] # fixed g
    nans <- which(as$.sum < g)
    nans <- as$.snp.id[nans]
    nans <- which(as$.snp.id %in% nans)
  }
  
  
  Nj <- as$.sum
  g <- as$.g
  meta <- .fix..call(as[,-..al_cols])
  as <- as.matrix(.fix..call(as[,..al_cols]))
  Qijg <- choose(Nj - as, g)/choose(Nj, g) # eqn 2a
  Pijg <- 1 - Qijg # eqn 2b
  
  alpha_g <- rowSums(Pijg)
  if(private){ # eqn 4
    meta <- cbind(meta, data.table::as.data.table(Qijg))
    meta[, paste0(al_cols, "_p") := lapply(.SD, prod), .SDcols = c(al_cols), by = .(facet, .snp.id)] # product across all populations
    p_al_cols <- paste0(al_cols, "_p")
    
    
    # next, need to divide out the value for each, since it's actually the product across all save this population
    acs <- as.matrix(.fix..call(meta[,..al_cols]))
    prods <- as.matrix(.fix..call(meta[,..p_al_cols]))
    pi_g <- rowSums(Pijg*(prods/acs), na.rm = TRUE)
    pi_g[which(is.nan(alpha_g))] <- NaN # add back in the NaN values where our sample size was too small in one pop. This got dropped in the above na.rm, which was needed due to the 0/0 prob.
    return(list(richness = alpha_g, pa = pi_g, g = g))
  }
  else{
    alpha_g[which(g < 2)] <- NaN
    if(exists("nans")){
      if(length(nans) > 0){
        alpha_g[nans] <- NaN
      }
    }
    return(list(richness = alpha_g, g = g))
  }
}




.do_CLD <- function(genos, snp.meta, sample.facet, sample.subfacet){
  proximity <- s1_position <- s2_position <- S <- NULL
  
  melt_cld <- function(CLD, snp.meta, sample.facet, sample.subfacet, skip_bind = FALSE){
    # CLD <- as.data.table(CLD)
    # CLD[,.snp.id := snp.meta$.snp.id]
    # prox <- melt.data.table(CLD, id.vars = ".snp.id")
    # colnames(prox) <- c(".snp.id1", ".snp.id2", "CLD")
    # 
    # if(!skip_bind){
    #   prox <- cbind(as.data.table(snp.meta[match(prox$.snp.id2, snp.meta$.snp.id),]), prox, as.data.table(snp.meta[match(prox$.snp.id1, snp.meta$.snp.id),]))
    #   prox$.snp.id1 <- prox$.snp.id2 <- NULL
    #   colnames(prox) <- c(paste0("s1_", colnames(snp.meta)), "CLD", paste0("s2_", colnames(snp.meta)))
    #   prox[,proximity := abs(s1_position - s2_position)]
    #   prox[,sample.facet := sample.facet]
    #   prox[,sample.subfacet := sample.subfacet]
    #   setcolorder(prox, c(1:ncol(snp.meta),
    #                       (ncol(snp.meta) + 2):(ncol(prox) - 3),
    #                       ncol(prox) - 2,
    #                       ncol(snp.meta) + 1,
    #                       (ncol(prox) - 1):ncol(prox)))
    # }
    
    prox <- cbind(as.data.table(snp.meta), as.data.table(CLD))
    prox <- reshape2::melt(prox, id.vars = colnames(snp.meta))
    prox <- cbind(prox, as.data.table(snp.meta[rep(1:nrow(snp.meta), each = nrow(snp.meta)),]))
    bad.col <- which(colnames(prox) == "variable")
    prox <- prox[,-bad.col, with = FALSE]
    colnames(prox) <- c(paste0("s1_", colnames(snp.meta)), "CLD", paste0("s2_", colnames(snp.meta)))
    # prox <- prox[-which(is.na(prox$CLD)),]
    prox <- as.data.table(prox)
    if(!skip_bind){
      prox[,proximity := abs(s1_position - s2_position)]
      prox[,sample.facet := sample.facet]
      prox[,sample.subfacet := sample.subfacet]
      setcolorder(prox, c(1:ncol(snp.meta),
                          (ncol(snp.meta) + 2):(ncol(prox) - 3),
                          ncol(prox) - 2,
                          ncol(snp.meta) + 1,
                          (ncol(prox) - 1):ncol(prox)))
    }
    return(prox)
  }
  
  # do the CLD calculation, add column/row names, NA out the lower triangle and diag
  ## CLD
  suppressWarnings(CLD <- stats::cor(t(genos), use = "pairwise.complete.obs")^2)
  ## matrix of complete case sample sizes
  complete.cases.matrix <- !is.na(t(genos))
  complete.cases.matrix <- crossprod(complete.cases.matrix)
  
  # fill in NAs
  CLD[which(lower.tri(CLD))] <- NA
  diag(CLD) <- NA
  complete.cases.matrix[is.na(CLD)] <- NA
  
  # add metadata and melt
  prox <- melt_cld(CLD, snp.meta, sample.facet, sample.subfacet)
  prox_S <- melt_cld(complete.cases.matrix, snp.meta, sample.facet, sample.subfacet, skip_bind = TRUE)
  colnames(prox_S)[which(colnames(prox_S) == "CLD")] <- "S"
  
  # merge
  prox[,S := prox_S$S] # should always be in the same order since the same stuff is getting melted and cbound.
  # prox <- merge.data.table(prox, prox_S) shouldn't need the expensive merge--same order always.
  prox <- prox[-which(is.na(CLD)),]
  
  # add column/row names
  colnames(CLD) <- snp.meta$position
  rownames(CLD) <- snp.meta$position
  return(list(prox = prox, LD_matrix = CLD, S = complete.cases.matrix))
}

# function to figure out which snps we are comparing to each
# outputs a nested list. Each entry in the list is a unique sample facet. In each of these lists is an entry for each unique subfacet level.
# in this is an entry for each snp that lists the snps it is compared to.
# If multiple entries would write to the same sample facet and subfacet, it will just add any new comparisons needed.
# this function also does the composite LD calculations, since that's most efficiently done here for each subfacet level.
.determine.comparison.snps <- function(x, facets, facet.types, window_sigma, window_step, transpose_windows = TRUE){
  if(is.list(facet.types)){
    facet.types <- facet.types[[2]]
  }
  
  if(is.numeric(window_sigma)){
    
    
    # correct format (for sample facets, noting samples)
    sample.facets <- facets
    if(any(facet.types == "snp")){sample.facets[facet.types == "snp"] <- ".base"}
    if(any(facet.types == "complex")){sample.facets[facet.types == "complex"] <- unlist(lapply(sample.facets[facet.types == "complex"], function(tf) .check.snpR.facet.request(x, tf, "snp")))}
    
    snp.facets <- facets
    if(any(facet.types == "sample")){snp.facets[facet.types == "sample"] <- ".base"}
    if(any(facet.types == "complex")){snp.facets[facet.types == "complex"] <- unlist(lapply(snp.facets[facet.types == "complex"], function(tf) .check.snpR.facet.request(x, tf, "sample")))}
    
    usamp <- unique(sample.facets)
    out <- vector("list", length(usamp))
    names(out) <- usamp
    cl <- data.table(facet = character(), subfacet = character(), 
                     snp.facet = character(), snp.subfacet = character(),
                     center = numeric(),
                     start = numeric(),
                     end = numeric())
    
    # loop through and apply
    for(i in 1:length(usamp)){
      comps <- .determine.ld.comps.window(x, snp.facets[which(sample.facets == usamp[i])], facet.types, window_sigma, window_step, transpose_windows = transpose_windows)
      
      tl <- .get.task.list(x, usamp[i])
      out[[i]] <- vector("list", length = nrow(tl))
      names(out[[i]]) <- tl[,2]
      for(j in 1:nrow(tl)){
        out[[i]][[j]] <- list(snps = comps[[1]],
                              samps = .fetch.sample.meta.matching.task.list(x, tl[j,]))
        
        cl <- rbind(cl,
                    cbind(facet = usamp[i], subfacet = tl[j,2], comps$cl))
      }
    }
    
    
    return(list(out, cl = cl))
  }
  
  # sub-subfunctions to get the options for the snp and sample facets
  get.samp.opts <- function(x, t.facet){
    sample.meta <- x@sample.meta[colnames(x@sample.meta) %in% t.facet]
    sample.meta <- sample.meta[,sort(colnames(sample.meta))]
    if(!is.data.frame(sample.meta)){
      sample.meta <- as.data.frame(sample.meta)
      colnames(sample.meta) <- colnames(x@sample.meta)[colnames(x@sample.meta) %in% t.facet]
    }
    sample.opts <- unique(sample.meta)
    if(!is.data.frame(sample.opts)){
      sample.opts <- as.data.frame(sample.opts, stringsAsFactors = F)
      colnames(sample.opts) <- facets[which(facets %in% colnames(x@sample.meta))]
    }
    sample.opts <- dplyr::arrange_all(sample.opts)
    
    return(list(sample.opts, sample.meta))
  }
  
  get.snp.opts <- function(x, t.facet){
    snp.meta <- x@snp.meta[colnames(x@snp.meta) %in% t.facet]
    snp.meta <- snp.meta[,sort(colnames(snp.meta))]
    if(!is.data.frame(snp.meta)){
      snp.meta <- as.data.frame(snp.meta)
      colnames(snp.meta) <- colnames(x@snp.meta)[colnames(x@snp.meta) %in% t.facet]
    }
    snp.opts <- unique(snp.meta)
    if(!is.data.frame(snp.opts)){
      snp.opts <- as.data.frame(snp.opts, stringsAsFactors = F)
      colnames(snp.opts) <- facets[which(facets %in% colnames(x@snp.meta))]
    }
    snp.opts <- dplyr::arrange_all(snp.opts)
    
    return(list(snp.opts, snp.meta))
  }
  
  # pull out just the sample facets
  sample.facets <- .check.snpR.facet.request(x, facets) # the sample level facets that we are working with.
  
  # initialize the output list
  out <- vector(mode = "list", length(sample.facets))
  names(out) <- sample.facets
  if(any(facet.types == "snp")){
    out <- c(out[-which(names(out) == ".base")], list(.base = list(.base = list(snps = vector("list", nrow(x)), samps = 1:nrow(x@sample.meta)))))
  }
  
  # loop through each facet, do different things depending on facet type
  for(i in 1:length(facets)){
    
    # grab the facet level list we are writing to.
    t.samp.facet <- .check.snpR.facet.request(x, facets[i])
    write.facet <- which(names(out) == t.samp.facet)
    if(length(write.facet) == 0){
      t.samp.facet <- ".base"
      write.facet <- which(names(out) == ".base")
    }
    t.facet <- unlist(.split.facet(facets[i]))
    
    this.out <- out[[write.facet]]
    
    # for sample only, need to loop through each sample level subfacet, then loop through all snps
    if(facet.types[i] == c("sample")){
      
      # get options
      opts <- get.samp.opts(x, t.facet)
      sample.opts <- opts[[1]]
      sample.meta <- opts[[2]]
      
      if(t.facet == ".base"){
        sample.opts <- matrix(".base")
        sample.meta <- matrix(".base", nrow = nrow(x@sample.meta))
      }
      
      # add snp/snp comparisons. Since the facet is simple, we do all snps. Do so with a loop through all subfacets
      if(is.null(this.out)){
        this.out <- vector("list", nrow(sample.opts))
        names(this.out) <- do.call(paste, as.data.frame(sample.opts))
      }
      
      for(j in 1:nrow(sample.opts)){
        
        # grab the subfacet level we are writing to.
        write.subfacet <- which(names(this.out) == paste(sample.opts[j,], collapse = " "))
        this.subfacet <- this.out[[write.subfacet]]
        
        
        if(is.null(this.subfacet)){
          samps.in.subfacet <- which(apply(sample.meta, 1, function(x) identical(as.character(x), as.character(sample.opts[j,]))))
          this.subfacet <- list(snps = vector("list", nrow(x)), samps = samps.in.subfacet)
        }
        
        # add comparisons for each snp. Note that the last snp, with no comparisons to do, will recieve a NULL
        for(k in 1:(nrow(x) - 1)){
          c.comps <- this.subfacet$snps[[k]]
          c.comps <- c(c.comps, (k + 1):nrow(x))
          dups <- which(duplicated(c.comps))
          if(length(dups) > 0){
            this.subfacet$snps[[k]] <- c.comps[-dups]
          }
          else{
            this.subfacet$snps[[k]] <- c.comps
          }
        }
        
        # add back to this.out
        this.out[[write.subfacet]] <- this.subfacet
      }
      
    }
    
    # for snp only, need to loop through each snp level subfacet, then through all snps on that subfacet
    else if(facet.types[i] == "snp"){
      # get the subfacet options
      opts <- get.snp.opts(x, t.facet)
      snp.opts <- opts[[1]]
      snp.meta <- opts[[2]]
      
      # add snp/snp comparisons. Since the facet is simple, we do all samples, but pick the correct snps. This will be at the .base facet and .base subfacet!
      for(j in 1:nrow(snp.opts)){
        snps.in.subfacet <- which(apply(snp.meta, 1, function(x) identical(as.character(x), as.character(snp.opts[j,]))))
        
        # add comparisons for each snp. Note that the last snp, with no comparisons to do, will recieve a NULL
        for(k in 1:(length(snps.in.subfacet) - 1)){
          if(length(snps.in.subfacet) == 1){
            next
          }
          c.comps <- this.out$.base$snps[[snps.in.subfacet[k]]]
          c.comps <- c(c.comps, snps.in.subfacet[(k + 1):length(snps.in.subfacet)])
          dups <- which(duplicated(c.comps))
          if(length(dups) > 0){
            this.out$.base$snps[[snps.in.subfacet[k]]] <- c.comps[-dups]
          }
          else{
            this.out$.base$snps[[snps.in.subfacet[k]]] <- c.comps
          }
        }
      }
    }
    
    # for complex, need to loop through first each sample level subfacet, then through the snp level subfacet, then through each snp on that subfacet.
    else if(facet.types[i] == "complex"){
      
      # get the subfacet sample options, snp and sample
      sample.opts <- get.samp.opts(x, t.facet)
      snp.opts <- get.snp.opts(x, t.facet)
      sample.meta <- sample.opts[[2]]
      sample.opts <- sample.opts[[1]]
      snp.meta <- snp.opts[[2]]
      snp.opts <- snp.opts[[1]]
      
      
      if(is.null(this.out)){
        this.out <- vector("list", nrow(sample.opts))
        names(this.out) <- do.call(paste, as.data.frame(sample.opts))
      }
      
      
      # for each sample level option, we make sure that we compare only within snp facet level
      for(j in 1:nrow(sample.opts)){
        
        # grab the subfacet level we are writing to.
        write.subfacet <-which(names(this.out) == paste(sample.opts[j,], collapse = " "))
        this.subfacet <- this.out[[write.subfacet]]
        
        if(is.null(this.subfacet)){
          samps.in.subfacet <- which(apply(sample.meta, 1, function(x) identical(as.character(x), as.character(sample.opts[j,]))))
          this.subfacet <- list(snps = vector("list", nrow(x)), samps = samps.in.subfacet)
        }
        
        for(l in 1:nrow(snp.opts)){
          snps.in.subfacet <- which(apply(snp.meta, 1, function(x) identical(as.character(x), as.character(snp.opts[l,]))))
          
          # add comparisons for each snp. Note that the last snp, with no comparisons to do, will recieve a NULL
          if(length(snps.in.subfacet) == 1){next} # if only one snp here, no LD to calculate
          for(k in 1:(length(snps.in.subfacet) - 1)){
            c.comps <- this.subfacet$snps[[snps.in.subfacet[k]]]
            c.comps <- c(c.comps, snps.in.subfacet[(k + 1):length(snps.in.subfacet)])
            dups <- which(duplicated(c.comps))
            if(length(dups) > 0){
              this.subfacet$snps[[snps.in.subfacet[k]]] <- c.comps[-dups]
            }
            else{
              this.subfacet$snps[[snps.in.subfacet[k]]] <- c.comps
            }
          }
        }
        
        # add back to this.out
        this.out[[write.subfacet]] <- this.subfacet
      }
    }
    
    # save the output for this subfacet.
    out[[write.facet]] <- this.out
  }
  
  # return
  return(out)
}

.determine.ld.comps.window <- function(x, facets = NULL, par = FALSE, window_sigma, window_step, verbose = FALSE, transpose_windows = TRUE){
  lev <- NULL
  
  # one task per facet -- standardize via get.task.list
  tasks <- .get.task.list(x, .check.snpR.facet.request(x, facets, "sample"))
  tasks <- unique(tasks[,3])
  comps <- vector("list", length(tasks))
  names(comps) <- tasks
  center_list <- data.table(snp.facet = character(),
                            snp.subfacet = character(),
                            center = numeric(),
                            start = numeric(),
                            end = numeric())
  
  # get comps for each
  osp <- options("scipen") # change scipen to avoid mangling names
  options(scipen = 999)
  for(i in 1:length(tasks)){
    comps[[i]] <- .average_windows(snp.meta(x), window_sigma, window_step, chr = tasks[i], triple_sig = FALSE, stats = NULL, gaussian = FALSE, nk = FALSE)
    # convert .snp.ids to rows
    comps[[i]] <- rapply(comps[[i]], function(z) match(z, snp.meta(x)$.snp.id), how = "list")
    
    ncomps <- vector("list", nrow(x))
    # back-process into which comparisons to do for each snp
    ucomps <- unlist(comps[[i]], recursive = F)
    if(length(ucomps) == 0){next}
    # window info
    cinf <- .split.facet(names(ucomps))
    cinf <- t(as.data.frame(cinf))
    cinf <- as.data.table(cinf)
    cinf[,lev := .paste.by.facet(cinf, colnames(cinf)[-ncol(cinf)])]
    ..keep_col <- NULL
    keep_col <- ncol(cinf):(ncol(cinf) - 1)
    cinf <- .fix..call(cinf[,..keep_col])
    cinf[,2] <- as.numeric(cinf[[2]])
    center_list <- rbind(center_list, data.table(snp.facet = tasks[i],
                                                 snp.subfacet = cinf[[1]],
                                                 center = cinf[[2]],
                                                 start = cinf[[2]] - window_sigma,
                                                 end = cinf[[2]] + window_sigma))
    
    # process
    if(transpose_windows){
      for(k in 1:length(ucomps)){
        cs <- t(utils::combn(ucomps[[k]], 2))
        for(h in 1:(length(ucomps[[k]]))){
          csn <- cs[which(ucomps[[k]][[h]] == cs[,1]),2]
          if(length(csn) > 0){ncomps[[ucomps[[k]][[h]]]] <- csn}
        }
      }
      comps[[i]] <- lapply(ncomps, unique)
    }
  }
  
  # condense comparisons across facets
  if(!transpose_windows){return(list(comps, cl = center_list))}
  comps <- do.call(mapply, c(FUN = c, comps, SIMPLIFY = FALSE))
  comps <- lapply(comps, unique)
  
  options(scipen = osp)
  return(list(comps, cl = center_list))
}


# function to figure out window values
.window_LD_averages <- function(prox, facets, window_gaussian, window_triple_sigma, window_step, window_sigma, x){

  # stop if empty
  if(nrow(prox) == 0){
    win_res <- as.data.table(matrix("", 0, 11))
    cn <- c("facet", "subfacet", "snp.facet", "snp.subfacet", "position", "sigma",
            "step", "nk.status", "gaussian", "n_snps", "triple_sigma")
    colnames(win_res) <- cn
    win_res$position <- numeric()
    win_res$sigma <- numeric()
    win_res$step <- numeric()
    win_res$nk.status <- logical()
    win_res$gaussian <- logical()
    win_res$n_snps <- numeric()
    win_res$triple_sigma <- logical()
    
    return(win_res)
  }
  
  sample.facet <- sample.subfacet <- s1_snp_subfacet <- s2_snp_subfacet <- sigma <- step <- nk.status <- gaussian <- triple_sigma <- NULL
  
  facet.types <- .check.snpR.facet.request(x, facets, "none", TRUE)[[2]]
  sample.facets <- facets
  sample.facets[facet.types == "snp"] <- ".base"
  sample.facets[facet.types == "complex"] <- unlist(lapply(sample.facets[facet.types == "complex"], function(tf) .check.snpR.facet.request(x, tf, "snp")))
  
  snp.facets <- facets
  snp.facets[facet.types == "sample"] <- ".base"
  snp.facets[facet.types == "complex"] <- unlist(lapply(snp.facets[facet.types == "complex"], function(tf) .check.snpR.facet.request(x, tf, "sample")))
  
  prox <- as.data.table(prox)
  ..stats <- NULL
  stats <- c("rsq", "Dprime", "pval", "CLD")
  stats <- stats[which(stats %in% colnames(prox))]
  
  win_res <- vector("list", 0)
  
  # loop through facets
  for(i in 1:length(sample.facets)){
    tasks <- .get.task.list(x, sample.facets[i])
    
    # loop through each sample facet option, get window averages
    for(j in 1:nrow(tasks)){
      tprox <- prox[sample.facet == sample.facets[i] &
                      sample.subfacet == tasks[j,2],]
      if(nrow(tprox) == 0){next}
      
      # run for either the base or the snp subfacets
      if(snp.facets[i] != ".base"){
        tprox[, s1_snp_subfacet := .paste.by.facet(tprox, paste0(c("s1_"), .split.facet(snp.facets[i])))]
        tprox[, s2_snp_subfacet := .paste.by.facet(tprox, paste0(c("s2_"), .split.facet(snp.facets[i])))]
        tprox <- tprox[s1_snp_subfacet == s2_snp_subfacet,]
        dups <- duplicated(tprox[,c("s1_snp_subfacet", "s1_position", "s2_position", "sample.subfacet")])
        if(any(dups)){
          tprox <- tprox[-which(dups),]
        }
        
        win_res[[length(win_res) + 1]] <- .average_windows(as.data.frame(tprox), sigma = window_sigma, step = window_step, stats = .fix..call(prox[,..stats]), chr = "s1_snp_subfacet", 
                                                           gaussian = window_gaussian, nk = FALSE, pairwise_snps = TRUE, triple_sig = FALSE)
        win_res[[length(win_res)]]$snp.subfacet <- win_res[[length(win_res)]]$s1_snp_subfacet
        win_res[[length(win_res)]]$s1_snp_subfacet <- NULL
      }
      else{
        dups <- duplicated(tprox[,c("s1_position", "s2_position", "sample.subfacet")])
        
        if(any(dups)){
          tprox <- tprox[-which(dups),]
        }
        
        win_res[[length(win_res) + 1]] <- .average_windows(as.data.frame(tprox), sigma = window_sigma, step = window_step, stats = .fix..call(prox[,..stats]), 
                                                           gaussian = window_gaussian, nk = FALSE, pairwise_snps = TRUE, chr = ".base", triple_sig = FALSE)
        win_res[[length(win_res)]]$snp.subfacet <- win_res[[length(win_res)]]$chr
        win_res[[length(win_res)]]$chr <- NULL
      }
      
      # add facet metadata
      win_res[[length(win_res)]]$facet <- sample.facets[i]
      win_res[[length(win_res)]]$subfacet <- tasks[j,2]
      win_res[[length(win_res)]]$snp.facet <- snp.facets[i]
    }
  }
  
  # combine and finalize outputs
  win_res <- data.table::rbindlist(win_res)
  win_res[,sigma := window_sigma/1000]
  if(window_triple_sigma){
    win_res[,sigma := window_sigma/3]
  }
  win_res[,step := window_step/1000]
  win_res[,nk.status := FALSE]
  win_res[,gaussian := window_gaussian]
  win_res[,triple_sigma := window_triple_sigma]

  col_ord <- c("facet", "subfacet", "snp.facet", "snp.subfacet", "position", "sigma",
               "step", "nk.status", "gaussian", "n_snps", "triple_sigma",
               stats)
  ..col_ord <- NULL
  win_res <- .fix..call(win_res[,..col_ord])
  
  return(win_res)
}


.rbind_list_fill_sparse <- function(l, fill_val = 0){
  cn <- unique(unlist(lapply(l, colnames)))
  out <- Matrix::Matrix(0, 0, length(cn), dimnames = list(NULL, cn))
  for(i in 1:length(l)){
    
    # fix any missing columns
    mcn <- which(!cn %in% colnames(l[[i]]))
    if(length(mcn) > 0){
      fill <- Matrix::Matrix(fill_val, 
                             nrow = nrow(l[[i]]), 
                             ncol = length(mcn),
                             dimnames = list(rownames(l[[i]]),
                                             cn[mcn]))
      l[[i]] <- cbind(l[[i]], fill)[,cn, drop = FALSE]
    }
    
    # make sure colunms are sorted identically
    ncord <- colnames(out)
    out <- rbind(out, l[[i]][,ncord, drop = FALSE])
  }
  return(out)
}

.rowMax_sparse <- function(x) slam::rollup(x, 2, FUN = max)[,1]
.rowMin_sparse <- function(x) slam::rollup(x, 2, FUN = min)[,1]
