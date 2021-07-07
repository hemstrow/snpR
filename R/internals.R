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
  return("snpRdata" %in% class(x))
}


#'Add facets to snpRdata objects
#'
#'Adds facets to snpRdata objects, following the rules described in
#'\code{\link{Facets_in_snpR}}. Internal, called by many other functions when
#'facets are requested.
#'
#'@param x snpRdata object
#'@param facets character. Facets to use.
#'
add.facets.snpR.data <- function(x, facets = NULL){
  #===========================binding issues for unquoted variables============
  .snp.id <- facet <- subfacet <- NULL
  
  
  if(is.null(facets[1])){return(x)}
  facets <- check.snpR.facet.request(x, facets)
  if(is.null(facets[1])){return(x)}
  #===========================turn into list========
  # need to fix any multivariate facets (those with a .)
  comp.facets <- grep("(?<!^)\\.", facets, perl = T)
  if(length(comp.facets) != 0){
    run.facets <- as.list(facets[-c(comp.facets)])
    facet.list <- c(run.facets, strsplit(facets[comp.facets], split = "(?<!^)\\.", perl = T))
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
  oac <- cbind(data.table::as.data.table(x@facet.meta[,c("facet", "subfacet", ".snp.id")]),
               data.table::as.data.table(x@ac)) # grab original ac with meta for later
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
      t.x <- x[,matches]
      tgs <- tabulate_genotypes(t.x, x@mDat)
      gs$gs <- plyr::rbind.fill.matrix(gs$gs, tgs$gs)
      gs$as <- plyr::rbind.fill.matrix(gs$as, tgs$as)
      gs$wm <- plyr::rbind.fill.matrix(gs$wm, tgs$wm)
      # fix NAs that show up when there are less called genotype options in one facet level than in all levels!
      if(ncol(gs$gs) != ncol(tgs$gs)){
        gs <- lapply(gs, function(x){x[is.na(x)] <- 0;x})
      }
      
      x@facet.meta <- rbind(x@facet.meta,
                            cbind(data.frame(facet = rep(paste0(facets, collapse = "."), nrow(tgs$gs)),
                                             subfacet = rep(paste0(sample.opts[i,], collapse = "."), nrow(tgs$gs)),
                                             facet.type = rep("sample", nrow(tgs$gs)), stringsAsFactors = F),
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
  invisible(utils::capture.output(nac <- format_snps(x, output = "ac", facets = added.facets)))
  nac <- data.table::as.data.table(nac)
  nac <- rbind(oac, nac[,c("facet", "subfacet", ".snp.id", "n_total","n_alleles", "ni1", "ni2")])
  nac <- dplyr::mutate_if(.tbl = nac, is.factor, as.character)
  nac <- dplyr::arrange(nac, .snp.id, facet, subfacet)
  nac <- as.data.frame(nac)
  x@ac <- nac[,-c(1:3)]
  
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

#'Grab information of the facets present in snpRdata objects.
#'
#'Internal, grabs via reference to face.meta socket rather than the facets
#'socket.
#'
#'@param x snpRdata object.
#'  
find.snpR.facets <- function(x){
  facets <- vector("list", length(x@facets))
  names(facets) <- x@facets
  for(i in 1:length(facets)){
    facets[[i]] <- unique(x@facet.meta[x@facet.meta$facet == names(facets)[i],]$subfacet)
  }
  return(facets)
}

#'Apply functions across snpR facets.
#'
#'Internal function to apply functions across snpR data. Facet desicnations
#'follow rules described in \code{\link{Facets_in_snpR}}. Null or "all" facet
#'designations follow typical rules.
#'
#'This function should never be called externally. Raw statistics with metadata
#'are returned and are not automatically merged into x.
#'
#'For examples, look at how this function is called in functions such as
#'calc_pi, calc_pairwise_fst, ect.
#'
#'Options:
#'
#'req:
#'
#'\itemize{ \item{gs: } genotype tables \item{ac: } ac formatted data
#'\item{meta.gs: } facet, .snp.id metadata cbound genotype tables.
#'\item{ac.stats: } ac data cbound to stats \item{meta.ac: } ac data cbound to
#'snp metadata. \item{snpRdata: } subset snpRdata object. }
#'
#'case:
#'
#'\itemize{ \item{ps: } stat calculated per snp. \item{ps.pf: } stat calculated
#'per snp, but split per facet (such as for private alleles, comparisons only
#'exist within a facet!) \item{facet.pairwise: } stat calculated pairwise
#'between facets, per snp otherwise. \item{ps.pf.psf: } stat calculated per snp,
#'but per both sample and snp facet. }
#'
#'@param x snpRdata object
#'@param facets character or NULL, default NULL. Facets to add.
#'@param req character. Data type required by fun. See description.
#'@param fun function. Function to apply to data.
#'@param case character, default "ps". Type of statistic required by fun. See
#'  description.
#'@param par numeric or FALSE, default FALSE. Number of parallel computing cores
#'  to use. Works for some cases/reqs.
#'@param ... other arguments to be passed to fun.
#'@param stats.type character, default "all". Other options "pairwise" or
#'  "stats". Type of statistic under to pull under the ps.pf.psf option.
#'@param response character, default NULL. Possible response variable name
#'  passed to specific functions.
#'@param maf numeric, defualt NULL. Possible maf, passed to some functions.
#'@param interpolate character, default NULL. Possible interpolation option,
#'  passed to some functions.
#'
#'@author William Hemstrom
apply.snpR.facets <- function(x, facets = NULL, req, fun, case = "ps", par = FALSE, ..., stats.type = "all", response = NULL, maf = FALSE, interpolate = NULL){
  #=====global bindings====

  if(!is.null(facets)){
    if(facets[1] == "all"){
      facets <- x@facets
    }
  }
  else {
    facets <- ".base"
  }
  
  
  if(case == "ps"){
    facets <- check.snpR.facet.request(x, facets)
    if(req == "gs"){
      # subset
      gs <- plyr::llply(x@geno.tables, function(y) subset(y, x@facet.meta$facet %in% facets))
      
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
      gs <- plyr::llply(x@geno.tables, function(y) subset(y, x@facet.meta$facet %in% facets))
      gs <- plyr::llply(gs, function(y) cbind(x@facet.meta[x@facet.meta$facet %in% facets, c("facet", "subfacet", ".snp.id")], y))
      
      # run the function indicated
      out <- data.table::as.data.table(fun(gs, ...))
      
      # bind metadata
      out <- cbind(data.table::as.data.table(x@facet.meta[x@facet.meta$facet %in% facets,]), out)
      
      # return
      return(out)
    }
    else if(req == "ac"){ # really simple, just here for consistancy across functions. you just... run it on the thing.
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
          pheno.facets <- c(check.snpR.facet.request(x, pheno.facets[-base.pheno.facets]), response)
        }
        else{
          pheno.facets <- response
        }
      }
      else{
        pheno.facets <- check.snpR.facet.request(x, pheno.facets)
      }
      x <- add.facets.snpR.data(x, pheno.facets)
      if(req == "cast.gs"){
        gs <- data.table::as.data.table(cbind(x@facet.meta, x@geno.tables$gs))
      }
      else{
        gs <- data.table::as.data.table(cbind(x@facet.meta, x@ac[, c("ni1", "ni2")]))
      }
      
      # initialize outputs
      comb <- vector("list", length = length(facets))
      
      # get the input gs data. Split by facet, since we need to extract metadata differently for each.
      for(i in 1:length(facets)){
        tgs <- gs[gs$facet == pheno.facets[i],]
        
        which.is.phenotype <- which(unlist(strsplit(pheno.facets[i], split = "(?<!^)\\.", perl = T)) == response)
        
        # add phenotype and subfacet column
        id <- strsplit(tgs$subfacet, split = "(?<!^)\\.", perl = T)
        l <- length(id[[1]])
        uid <- unlist(id)
        p.vals <- seq(from = which.is.phenotype, to = length(uid), by = l)
        subfacet <- lapply(id, FUN = function(x) x[-which.is.phenotype])
        subfacet <- unlist(lapply(subfacet, FUN = paste, collapse = "."))
        tgs$subfacet <- subfacet
        tgs$phenotype <- uid[p.vals]
        tgs$facet.type <- check.snpR.facet.request(x, facets[i], "none", T)[[2]]
        
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
        comb[[i]]$facet <- facets[i]
      }
      
      comb <- rbindlist(comb)
      mcols <- 1:ncol(x@facet.meta)
      meta <- comb[,mcols, with = FALSE]
      comb <- comb[,-mcols, with = FALSE]
      
      # call function
      out <- fun(comb, ...)
      n.ord <- (1:ncol(meta))[-which(colnames(meta) == ".snp.id")]
      n.ord <- c(n.ord, which(colnames(meta) == ".snp.id"))
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
          suppressWarnings(invisible(utils::capture.output(sub.x <- subset_snpR_data(x, facets = opts[i,1], subfacets = opts[i,2]))))
        }
        else{
          suppressWarnings(invisible(utils::capture.output(sub.x <- subset_snpR_data(x, facets = opts[i,1], subfacets = opts[i,2],
                                                                              snp.facets = opts[i,3], snp.subfacets = opts[i,4]))))
        }
        
        suppressWarnings(invisible(utils::capture.output(sub.x <- filter_snps(sub.x, maf = maf))))
        out <- fun(sub.x = sub.x, opts.list = opts.list, ...)
        if(!is.data.frame(out)){
          out$importance <- cbind(facet = opts[i,1], subfacet = opts[i,2], out$importance)
          return(out)
        }
        else{
          out <- cbind(facet = opts[i,1], subfacet = opts[i,2], out)
          return(out)
        }
      }
      #==========run===============
      # check facets
      facets <- check.snpR.facet.request(x, facets)
      
      # get options
      opts.list <- get.task.list(x, facets)
      
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
        cl <- snow::makeSOCKcluster(par)
        doSNOW::registerDoSNOW(cl)
        
        #prepare reporting function
        ntasks <- nrow(opts.list)
        progress <- function(n) cat(sprintf("Part %d out of", n), ntasks, "is complete.\n")
        opts <- list(progress=progress)
        
        
        cat("Begining run.\n")
        
        # run the LD calculations
        ## suppress warnings because you'll get wierd ... warnings. Not an issue in the non-parallel version.
        suppressWarnings(out <- foreach::foreach(q = 1:ntasks, .inorder = FALSE,
                                                 .options.snow = opts, .export = c("data.table"), .packages = "snpR") %dopar% {
                                                   run.one.task(opts.list, q)
                                                 })
        
        #release cores
        parallel::stopCluster(cl)
        doSNOW::registerDoSNOW()
      }
      
      # run in parallel
      
      # combine
      ## for rf
      if(!is.data.frame(out[[1]])){
        
        # initialize
        bind_list <- vector("list", length(out))
        nvec <- character(length(out))
        
        # bind and grab names
        for(i in 1:length(bind_list)){
          bind_list[[i]] <- out[[i]]$importance
          out[[i]] <- out[[i]][-1]
          nvec[i] <- paste0(bind_list[[i]]$facet[1], "_", bind_list[[i]]$subfacet[1])
        }
        suppressWarnings(bind_list <- dplyr::bind_rows(bind_list))
        names(out) <- nvec
        
        # return
        return(list(stats = bind_list, models = out))
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
        gs <- plyr::llply(x@geno.tables, function(y) subset(y, x@facet.meta$facet == facets[i]))
        gs <- plyr::llply(gs, function(y) cbind(x@facet.meta[x@facet.meta$facet == facets[i], c("facet", "subfacet", ".snp.id")], y))
        
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
    if(req == "meta.ac" | req == "pos.all.stats"){
      x@facet.meta$facet <- as.character(x@facet.meta$facet)
      x@facet.meta$subfacet <- as.character(x@facet.meta$subfacet)
      
      # subfunctions:
      
      ## a function to run 'func' on one iteration/row of the task.list. q is the iteration. Par is needed only for progress printing.
      run.one.loop <- function(stats_to_use, meta, task.list, q, par){
        if(par == FALSE){cat("Sample Subfacet:\t", as.character(task.list[q,2]), "\tSNP Subfacet:\t", as.character(task.list[q,4]), "\n")}
        
        # get comps and run
        ## figure out which data rows contain matching sample facets
        sample.matches <- which(apply(meta[,1:2], 1, function(x) identical(as.character(x), as.character(task.list[q,1:2]))))
        
        snp.facets <- unlist(strsplit(task.list[q,3], "(?<!^)\\.", perl = T))
        
        if(snp.facets[1] != ".base"){
          # figure out which data rows contain matching snp facets
          snp.cols <- meta[,snp.facets]
          snp.cols <- do.call(paste, as.data.frame(snp.cols))
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
        assign("last.warning", NULL, envir = baseenv())
        res <- fun(cbind(meta[run.lines,], stats_to_use[run.lines,]), ...)
        
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
      
      
      # get the statistics needed by the function and the task list. For now, expects only meta cbound to ac or snp-wise statistics.
      if(req == "meta.ac"){
        stats_to_use <- x@ac
        task.list <- get.task.list(x, facets)
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
          nk <- matrixStats::rowSums2(x@geno.tables$as)
          if(data.table::is.data.table(x@stats)){
            stats_to_use <- cbind(nk = nk, x@stats[,cols.to.use, with = FALSE])
          }
          else{
            stats_to_use <- cbind(nk = nk, x@stats[,cols.to.use])
          }
          numeric.cols <- which(sapply(stats_to_use, class) %in% c("numeric", "integer"))
          ## save
          stats_to_use <- stats_to_use[,numeric.cols, with = FALSE]
          task.list <- get.task.list(x, facets)
          meta.to.use <- x@facet.meta
          
        }
        else{
          # grab the correct data columns from the pairwise stats dataset, and reorder so nk is first
          pos.col <- which(colnames(x@pairwise.stats) == "position")
          start.col <- which(colnames(x@pairwise.stats) == "comparison") + 1
          cols.to.use <- start.col:ncol(x@pairwise.stats)
          stats_to_use <- x@pairwise.stats[,cols.to.use, with = FALSE]
          nk.col <- which(colnames(stats_to_use) == "nk")
          n.col.ord <- c(nk.col, (1:ncol(stats_to_use))[-nk.col])
          stats_to_use <- stats_to_use[,n.col.ord, with = FALSE]
          task.list <- get.task.list(x, facets, source = "pairwise.stats")
          
          # re-order the meta and save
          meta.to.use <- as.data.frame(x@pairwise.stats[,1:which(colnames(x@pairwise.stats) == "comparison")])
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
        cl <- snow::makeSOCKcluster(par)
        doSNOW::registerDoSNOW(cl)
        
        #prepare reporting function
        ntasks <- nrow(task.list)
        progress <- function(n) cat(sprintf("Part %d out of", n), ntasks, "is complete.\n")
        opts <- list(progress=progress)
        
        
        cat("Begining run.\n")
        
        # run the LD calculations
        ## suppress warnings because you'll get wierd ... warnings. Not an issue in the non-parallel version.
        suppressWarnings(out <- foreach::foreach(q = 1:ntasks, .inorder = TRUE,
                                                 .options.snow = opts, .export = "data.table") %dopar% {
                                                   run.one.loop(stats_to_use, meta.to.use, task.list, q, TRUE)
                                                 })
        
        #release cores
        parallel::stopCluster(cl)
        doSNOW::registerDoSNOW()
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
  else if(case == "per_sample"){

    split_facets <- strsplit(facets, "(?<!^)\\.", perl = T)
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
          genotypes <- x[ident,, drop = FALSE]
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
        opt.names <- check.snpR.facet.request(x, paste0(colnames(opts[[i]]), collapse = "."), remove.type = "sample")
        opt.names <- unlist(strsplit(opt.names, "(?<!^)\\.", perl = T))
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
    meta[is.na(meta)] <- ".base"
    out[,1:(ncol(out) - 1)] <- meta
    return(out)
  }
  
}

#'Merge newly calculated stats into a snpRdata object.
#'
#'Internal function to quickly and accurately merge newly calculated statistics
#'into a snpRdata object. This should never be called externally. Takes and
#'returns the same snpRdata object, with new data.
#'
#'Mostly relies on the smart.merge subfunction, which must be edited with care
#'and testing. LD merging is independant.
#'
#'Type options: stats, pairwise, window.stats, pairwise.window.stats, and LD,
#'corresponding to slots of the snpRdata S4 object class.
#'
#'For examples, see functions that use merge.snpR.stats, such as calc_pi or
#'calc_pairwise_ld.
#'
#'@param x snpRdata object
#'@param stats data.frame/table/list. New data to be added to the existing
#'  snpRdata object, x.
#'@param type character, default "stats". The type of statistic to merge, see
#'  list in description.
merge.snpR.stats <- function(x, stats, type = "stats"){
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
    n.s <- smart.merge(stats, x@stats, meta.cols, starter.meta)
    x@stats <- n.s
  }
  else if(type == "pairwise"){
    # merge and return
    meta.cols <- c(colnames(stats)[1:(which(colnames(stats) == "comparison"))], colnames(x@snp.meta))
    starter.meta <- c("facet")
    n.s <- smart.merge(stats, x@pairwise.stats, meta.cols, starter.meta)
    x@pairwise.stats <- n.s
  }
  else if(type == "window.stats"){
    # merge and return
    meta.cols <- c("facet", "subfacet", "position", "sigma", "n_snps", "snp.facet", "snp.subfacet", "step", "nk.status", colnames(x@snp.meta))
    starter.meta <- c("facet", "subfacet", "snp.facet", "snp.subfacet", "position")
    n.s <- smart.merge(stats, x@window.stats, meta.cols, starter.meta)
    x@window.stats <- n.s
  }
  else if(type == "pairwise.window.stats"){
    meta.cols <- c("facet", "subfacet", "position", "sigma", "n_snps", "snp.facet", "snp.subfacet", "step", "nk.status", colnames(x@snp.meta))
    starter.meta <- c("facet", "subfacet", "snp.facet", "snp.subfacet", "position")
    n.s <- smart.merge(stats, x@pairwise.window.stats, meta.cols, starter.meta)
    x@pairwise.window.stats <- n.s
  }
  else if(type == "sample.stats"){
    meta.cols <- c("facet", "subfacet", colnames(stats)[1:(which(colnames(stats) == ".sample.id"))])
    starter.meta <- meta.cols
    n.s <- smart.merge(stats, x@sample.stats, meta.cols, starter.meta)
    
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
      # deal with prox table using smart.merge
      start.meta <- colnames(stats$prox)[1:which(colnames(stats$prox) == "proximity")]
      x@pairwise.LD$prox <- smart.merge(stats$prox, x@pairwise.LD$prox,
                                        meta.names = c(start.meta, "sample.facet", "sample.subfacet"),
                                        starter.meta = start.meta)
      
      # Deal with matrices using the merge.lists utility in snpR.
      x@pairwise.LD$LD_matrices <- merge.lists(x@pairwise.LD$LD_matrices, stats$LD_matrices)
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
    if(length(x@pop.stats) == 0){
      x@pop.stats <- data.table::as.data.table(stats)
    }
    else{
      meta.names <- c("facet", "subfacet")
      starter.meta <- meta.names
      x@pop.stats <- smart.merge(x@pop.stats, stats, meta.names, starter.meta)
    }
  }
  else if(type == "weighted.means"){
    if(nrow(x@weighted.means) == 0){
      x@weighted.means <- as.data.table(stats)
    }
    else{
      meta.names <- c("facet", "subfacet", "snp.facet", "snp.subfacet", colnames(x@facet.meta)[-c(1:3, ncol(x@facet.meta))])
      starter.meta <- meta.names
      x@weighted.means <- smart.merge(x@weighted.means, stats, meta.names, starter.meta)
    }
  }
  
  return(x)
}

#' Merging function used in merge.snpR.stats and elsewhere
#' @param n.s new stats
#' @param o.s old stats
#' @param meta.names names of the metadata columns, usually everything up to .snp.id
#' @param starter.meta any metadata columns that should specifically be put at the start of the output data (such as facet, subfacet, facet.type)
smart.merge <- function(n.s, o.s, meta.names, starter.meta){
  n.s <- data.table::as.data.table(n.s)
  o.s <- data.table::as.data.table(o.s)
  if(all(colnames(n.s) %in% colnames(o.s))){
    if(isTRUE(all.equal(n.s, o.s, ignore.col.order = T, ignore.row.order = T, check.attributes = F))){return(n.s)}
  }
  
  if(nrow(o.s) == 0){
    return(n.s)
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
  return(m.s)
}

#'Get details on and clean up facet arguments.
#'
#'Given a snpRdata object, this function cleans up and determines the types of
#'requested facets. Internal. If requested, can clean facet requests to purge a
#'specific facet type (snp, sample, complex). Used in functions that run only on
#'a specific type of facet.
#'
#'Facet designation of NULL or "all" follows the typical rules.
#'
#'Facets designated as described in \code{\link{Facets_in_snpR}}.
#'
#'@param x snpRdata object to compare to
#'@param facets character. facets to check.
#'@param remove.type character. "snp", "sample", "complex", "simple" or anything
#'  else, typically "none". Type of facet to remove.
#'@param return.type logical, default FALSE. If true, returns both facets and
#'  facet types ("snp", "sample", "complex", or ".base") as a length two list.
#'@param fill_with_base logical, default TRUE. If true, fills the returned facets
#'  with .base if nothing is left after facets are removed. Otherwise returns null in this case.
#'
#'@author William Hemstrom
check.snpR.facet.request <- function(x, facets, remove.type = "snp", return.type = FALSE, fill_with_base = TRUE){
  if(any(facets == "all")){
    facets <- x@facets
  }
  if(is.null(facets)){
    facets <- ".base"
    if(return.type){
      return(list(facets, ".base"))
    }
    else{
      return(".base")
    }
  }
  
  
  # remove the facet parts as requested.
  facets <- strsplit(facets, "(?<!^)\\.", perl = T)
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
    facets <- facets[-which(to.remove)]
    facet.types <- facet.types[-which(to.remove)]
    if(fill_with_base){
      if(!".base" %in% facets){
        facets <- c(facets, ".base")
        facet.types <- c(facets, ".base")
      }
    }
  }
  
  # fix if we've removed everything, return the base facet if fill_with_base is TRUE
  if(fill_with_base){
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
  
  # remove duplicates and return
  dups <- duplicated(facets)
  if(any(dups)){
    facets <- facets[-which(dups)]
    facet.types <- facet.types[-which(dups)]
  }
  if(return.type){
    return(list(unlist(facets), facet.types))
  }
  else{
    return(unlist(facets))
    
  }
}

#'Tabulate allele and genotype counts at each locus.
#'
#'\code{tabulate_genotypes} creates matricies containing counts of observed
#'alleles and genotypes at each locus.
#'
#'This function is pirmarily used interally in several other funcitons.
#'
#'@param x Input raw genotype data, where columns are individuals and rows are
#'  snps. No metadata.
#'@param mDat Character. How are missing \emph{genotypes} noted?
#'@param verbose Logical. Should the function report progress?
#'
#'@return A list of matrices. gs is the genotype matrix, as is the allele
#'  matrix, and wm is the genotype matrix with missing genotypes.
#'
#'@author William Hemstrom
#'
tabulate_genotypes <- function(x, mDat, verbose = F){
  
  # fix for if x is a vector (only one individual) and convert to data.table
  if(!is.data.frame(x)){
    x <- data.frame(samp = x)
  }
  x <- data.table::setDT(x)
  
  
  # get a genotype table
  snp_form <- nchar(x[1,1])   # get information on data format
  x <- data.table::melt(data.table::transpose(x, keep.names = "samp"), id.vars = "samp") # transpose and melt
  
  gmat <- data.table::dcast(data.table::setDT(x), variable ~ value, value.var='value', length) # cast
  gmat <- gmat[,-1]
  mis.cols <- -which(colnames(gmat) == mDat)
  if(length(mis.cols) > 0){
    tmat <- gmat[,mis.cols, with = FALSE] # remove missing data
  }
  else{
    tmat <- gmat
  }
  
  #get matrix of allele counts
  #initialize
  hs <- substr(colnames(tmat),1,snp_form/2) != substr(colnames(tmat), (snp_form/2 + 1), snp_form*2) # identify heterozygotes.
  if(verbose){cat("Getting allele table...\n")}
  as <- unique(unlist(strsplit(paste0(colnames(tmat)), "")))
  amat <- data.table::as.data.table(matrix(0, nrow(gmat), length(as)))
  colnames(amat) <- as
  
  #fill in
  for(i in 1:length(as)){
    b <- grep(as[i], colnames(tmat))
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
  return(list(gs = as.matrix(tmat), as = as.matrix(amat), wm = as.matrix(gmat)))
}




#'Interpolate sn formatted data.
#'
#'An internal function to interpolate sn formatted data using either the
#'bernoulli or expected minor allele count approaches. Typically entirely
#'internal, called via format_snps.
#'
#'Interpolating missing data in sn formatted data is useful for PCA, genomic
#'prediction, tSNE, and other methods. Specify method = "af" to insert the
#'expected number of minor alleles given SNP allele frequency or "bernoulli" to
#'do binomial draws to determine the number of minor alleles at each missing
#'data point, where the probability of drawing a minor allele is equal to the
#'minor allele frequency. The expected number of minor alleles based on the
#'later method is equal to the interpolated value from the former, but the later
#'allows for multiple runs to determine the impact of stochastic draws and is
#'generally prefered and required for some downstream analysis. It is therefore
#'the default. As a slower but more accurate alternative to "af" interpolation,
#'"iPCA" may be selected. This an iterative PCA approach to interpolate based on
#'SNP/SNP covariance via \code{\link[missMDA]{imputePCA}}. If the ncp argument
#'is not defined, the number of components used for interpolation will be
#'estimated using \code{\link[missMDA]{estim_ncpPCA}}. In this case, this method
#'is much slower than the other methods, especially for large datasets. Setting
#'an ncp of 2-5 generally results in reasonable interpolations without the time
#'constraint.
#'
#'@param sn data.frame. Input sn formatted data, as produced by
#'  \code{\link{format_snps}}. Note that \code{\link{format_snps}} has an option
#'  to automatically call this function during formatting.
#'@param method character, default "bernoulli". Method to used for
#'  interpolation, either bernoulli or af. See details.
#'@param ncp numeric or NULL, default NULL. Number of components to consider for
#'  iPCA sn format interpolations of missing data. If null, the optimum number
#'  will be estimated, with the maximum specified by ncp.max. This can be very
#'  slow.
#'@param ncp.max numeric, default 5. Maximum number of components to check for
#'  when determining the optimum number of components to use when interpolating
#'  sn data using the iPCA approach.
#'
#'@author William Hemstrom
interpolate_sn <- function(sn, method = "bernoulli", ncp = NULL, ncp.max = 5){
  if(!method %in% c("bernoulli", "af", "iPCA")){
    stop("Unaccepted interpolation method. Accepted methods: bernoulli, af.\n")
  }
  
  if(method %in% c("bernoulli", "af")){
    # find allele frequencies
    sn <- as.matrix(sn)
    sn <- t(sn)
    af <- colMeans(sn, na.rm = T) # this is the allele frequency of the "1" allele
    
    # identify all of the NAs and the columns that they belong to
    NAs <- which(is.na(sn)) # cells with NA
    NA.cols <- floor(NAs/nrow(sn)) + 1 # figure out the column
    adj <- which(NAs%%nrow(sn) == 0) # adjust for anything that sits in the last row, since it'll get assigned the wrong column
    NA.cols[adj] <- NA.cols[adj] - 1
    
    # do interpolation for each missing data point
    if(method == "bernoulli"){
      ndat <- rbinom(length(NAs), 2, af[NA.cols])
    }
    else if(method == "af"){
      ndat <- af[NA.cols]
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

#'Do sanity checks for many window specific functions in snpR.
#'
#'Internal function only. Since many window related functions have similar
#'sanity checks, they have been integrated here. For examples, see use in
#'functions like \code{\link{do_bootstraps}}.
#'
#'@param x snpRdata object.
#'@param sigma numeric. Window size, in kilobases.
#'@param step numeric or NULL. Step size, in kilobases
#'@param stats.type character, either "single" or "pairwise". The type of stats
#'  being checked, see documentation for \code{\link{calc_smoothed_averages}}.
#'@param nk Logical. Determines if nk values are to be used.
#'@param facets character or NULL. Defines facets to run, following typical
#'  rules as described in \code{\link{Facets_in_snpR}}.
#'@param stats character or NULL, default NULL. Statistics (pi, etc.) used,
#'  really only for bootstrapping.
#'@param good.types character or NULL, default NULL. Good statistics, for use
#'  with the stats arguement.
#'
#'@author William Hemstrom
sanity_check_window <- function(x, sigma, step, stats.type, nk, facets, stats = NULL, good.types = NULL){
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
      warn.msg <- c(warn.msg, "Sigma and/or ws are larger than typically expected. Reminder: sigma and ws are given in megabases!")
    }
    else if(sigma <= 50){
      warn.msg <- c(warn.msg, paste0("Provided sigma is very small: ", sigma, " mb!"))
    }
  }
  else{
    if(!all(is.numeric(sigma), length(sigma) == 1)){
      msg <- c(msg, "sigma must be a numeric vector of length 1. step may be the same or NULL.")
    }
    if(sigma >= 500){
      warn.msg <- c(warn.msg, "Sigma is larger than typically expected. Reminder: sigma is given in megabases!")
    }
    else if(sigma <= 50){
      warn.msg <- c(warn.msg, paste0("Provided sigma is very small: ", sigma, " mb!"))
    }
  }
  
  # position needs to be available
  if(!any(colnames(x@snp.meta) == "position")){
    msg <- c(msg, "No column named position found in snp metadata.")
  }
  
  # facets
  if(is.null(facets[1]) & any(stats.type == "pairwise")){
    msg <- c(msg, "If no facets are provided, pairwise statistics cannot be smoothed. Please specify stats.type or statistics = 'single'")
  }
  facet.types <- check.snpR.facet.request(x, facets, remove.type = "none", return.type = T)
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
  
  if(length(warn.msg) > 0){
    warning(warn.msg)
  }
  if(length(msg) > 0){
    stop(paste0(msg, collapse = "\n  "))
  }
  
  return(TRUE)
}



#' List unique tasks/options for facets. Internal function to get a list of
#' tasks to run (one task per unique sample/snp level facet!). The source
#' arguement specifies what kind of statistics are being grabbed.
#'
#' @param x snpRdata
#' @param facets Facets to generate tasks for
#' @param source "stats" or "pairwise.stats", default "stats". Type of
#'   comparison to get jobs for. Note that the latter only works if pairwise
#'   stats have actually been calculated.
#' @author William Hemstrom
get.task.list <- function(x, facets, source = "stats"){
  facets <- check.snpR.facet.request(x, facets, "none", F)
  task.list <- matrix("", ncol = 4, nrow = 0) # sample facet, sample subfacet, snp facet, snp.subfacet
  
  if(source == "stats"){
    meta.to.use <- x@facet.meta
  }
  else if (source == "pairwise.stats"){
    meta.to.use <- as.data.frame(x@pairwise.stats[,1:which(colnames(x@pairwise.stats) == "comparison")], stringsAsFactors = F)
    meta.to.use$subfacet <- meta.to.use$comparison
  }
  
  snp.facet.list <- vector("list", length = length(facets))
  for(i in 1:length(facets)){
    t.facet <- facets[i]
    t.facet <- unlist(strsplit(t.facet, split = "(?<!^)\\.", perl = T))
    t.facet.type <- check.snpR.facet.request(x, t.facet, remove.type = "none", return.type = T)[[2]]
    
    # sample facets
    if(any(t.facet.type == "sample")){
      t.sample.facet <- check.snpR.facet.request(x, facets[i], remove.type = "snp")
      t.sample.meta <- meta.to.use[meta.to.use$facet == t.sample.facet, c("subfacet")]
      sample.opts <- unique(t.sample.meta)
      t.sample.meta <- meta.to.use[,c("facet", "subfacet")]
      if(is.null(nrow(sample.opts))){
        sample.opts <- as.data.frame(sample.opts, stringsAsFactors = F)
        t.sample.meta <- as.data.frame(t.sample.meta, stringsAsFactors = F)
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
      t.snp.meta <- meta.to.use[,..use.facet]
      snp.opts <- unique(t.snp.meta)
      t.snp.facet <- check.snpR.facet.request(x, facets[i], remove.type = "sample")
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
    all.opts.1 <- do.call(paste, as.data.frame(all.opts.1))
    all.opts.2 <- matrix(rep(t(snp.opts), nrow(sample.opts)), ncol = ncol(snp.opts), byrow = TRUE)
    all.opts.2 <- do.call(paste, as.data.frame(all.opts.2))
    t.task.list <- cbind(t.sample.facet, all.opts.1, t.snp.facet, all.opts.2)
    task.list <- rbind(task.list, t.task.list)
  }
  
  return(task.list)
}

#' Internal to find which sample metadata rows match a row from a task list.
#' @param x snpRdata object to check
#' @param task.list.row A row of an object created by
#'   \code{\link{get.task.list}} to look up info about.
fetch.sample.meta.matching.task.list <- function(x, task.list.row){
  task.list.row <- unlist(task.list.row) # in case it is passed with drop = F for some reason
  
  if(task.list.row[1] == ".base"){
    matches <- 1:nrow(x@sample.meta)
  }
  else{
    t.facets <- unlist(strsplit(task.list.row[1], "(?<!^)\\.", perl = T))
    matches <- which(do.call(paste, c(dat = as.data.frame(x@sample.meta[,t.facets]), sep = ".")) == task.list.row[2])
  }
  
  return(matches)
}


#' Merge lists element-to-element
#'
#' Merges list 1 into list 2, adding any elements without matching equivalents
#' to list 2.
#'
#' Merges element to element. Elements in list 1 with matching elements in list
#' 2 will be replaced.
#'
#' Do this by getting the names and "pathways" of all of the deepest level
#' objects, then looping through list 2 data and adding from list 1 that isn't
#' present at the same "pathway".
#'
#' @param list1 list. The first list. Elements with identical names at all
#'   levels in list 2 will be *replaced*.
#' @param list2 list. The second list. Elements in list 2 with identical names
#'   found in list 1 will replace those elements.
#' @param possible_end_level_names vector of possible stopping point names
#'   (bottom level)
#'
#' @author William Hemstrom
merge.lists <- function(list1, list2, possible_end_level_names = c("Dprime", "rsq", "pval", "CLD", "S")){
  # prunes and prepares data on the names and nest levels of a list
  prune.names <- function(list){
    ln <- utils::capture.output(utils::str(list))
    ln <- ln[which(!grepl("attr\\(\\*,", ln))]
    ns <- character()
    pn <- vector("list")
    collapsed.names <- vector("list")
    
    # get the levels of each name
    lev <- gsub("\\$.+", "", ln[2:length(ln)])
    lev <- stringr::str_count(lev, "\\.\\.")
    
    # get the name at each level
    clean.names <- stringr::str_extract(ln, "\\$[^:]*")
    clean.names <- gsub("\\$", "", clean.names)
    clean.names <- gsub(":", "", clean.names)
    clean.names <- gsub("num \\[.+$", "", clean.names)
    clean.names <- gsub("chr \\[.+$", "", clean.names)
    clean.names <- gsub(" ", "", clean.names)
    clean.names[clean.names == ""] <- NA
    terminal <- which(clean.names %in% possible_end_level_names)
    cl <- numeric()
    for(i in 2:length(ln)){
      if(i %in% terminal){
        pn[[length(pn) + 1]] <- c(ns, clean.names[i])
        collapsed.names[[length(pn)]] <- paste(pn[[length(pn)]], collapse = "")
        ns <- ns[1:(lev[i] - 1)]
        cl <- c(cl, lev[i])
      }
      else{
        if(!is.na(clean.names[i])){
          if(lev[i] <= length(ns)){
            ns[lev[i]] <- clean.names[i]
          }
          else{
            ns <- c(ns, clean.names[i])
          }
        }
      }
    }
    
    return(list(n = pn, l = lev, cn = unlist(collapsed.names)))
  }
  
  
  all.names.1 <- prune.names(list1)
  all.names.2 <- prune.names(list2)
  ## add anything present in 1 but not 2 to 2.
  for(i in 1:length(all.names.1[[1]])){
    # if not there, add to stats
    if(!all.names.1$cn[i] %in% all.names.2$cn){
      loc <- paste(paste0("[['", all.names.1$n[[i]], "']]"), collapse = "")
      call <- paste0("list2",
                     loc,
                     " <- ",
                     "list1", loc)
      eval(parse(text=call))
    }
  }
  return(list2)
}



#' Update list of calculated stats for a vector of facet names
#'
#' @param x snpRdata object to update
#' @param facets character. Facets to update, see \code{\link{Facets_in_snpR}}
#' @param stats character. Name of the facet to update as calculated.
#' @param remove.type character, default none. If provided, will grab only the
#'   remove the requested facet type (omit the snp part, for example).
#'
#' @return A snpRdata object identical to x but with calced stats updated.
#'
#' @author William Hemstrom
update_calced_stats <- function(x, facets, stats, remove.type = "none"){
  facets <- check.snpR.facet.request(x, facets, remove.type = remove.type)
  
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

#' Check if a given stat or vector of stats have been calculated any number of
#' facets.
#'
#' @param x snpRdata object to check
#' @param facets character. See \code{\link{Facets_in_snpR}}
#' @param stats character. Names of stats to check.
#'
#' @return A named list with an entry for each facet containing a named logical
#'   vector indicating if the provided stats have been calculated yet.
#'
#' @author William Hemstrom
check_calced_stats <- function(x, facets, stats){
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


#' Utility to check if a required package is installed and install it if needed
#' @param pkg character, name of the package to check.
#' @param install.type character, default "basic". Where should the package be
#'   sourced from? Options: "basic": CRAN, "github": github, "bioconductor":
#'   bioconductor.
#' @param source character, default NULL. If install.type = 'github', gives the
#'   source from which to install.
check.installed <- function(pkg, install.type = "basic", source = NULL){
  if(!pkg %in% rownames(utils::installed.packages())){
    say <- paste0("Package '", pkg, "' not found.")
    cat(say, "")
    cat("Install? (y or n)\n")
    
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
        if(!"BiocManager" %in% installed.packages()){
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

#' Correct for multiple testing.
#' @param p data.frame/data.table containing p values and any facets to split by
#' @param pcol name of the column containing p values
#' @param levs names of the columns containing facet info to split by
#' @param methods methods to use
#' @param case case (overall, within each subfacet, or within each facet) at
#'   which to apply corrections
fwe_correction <- function(p, pcol = NULL, levs = c("subfacet", "facet"), methods = c("bonferroni", "holm", "BH", "BY"), 
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

#' Tiny internal used when converting from a phased, transposed dataset (like
#' from an MS output)
#' @param x input matrix
convert_2_to_1_column <- function(x){
  if(!is.matrix(x)){x <- as.matrix(x)}
  ind.genos <- x[,seq(1,ncol(x), by = 2)] + x[,seq(2,ncol(x), by = 2)]
  ind.genos <- matrix(ind.genos, ncol = ncol(x)/2) # rematrix and transpose!
  return(ind.genos)
}
