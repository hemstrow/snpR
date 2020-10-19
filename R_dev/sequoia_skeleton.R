run_sequoia <- function(x, facets = NULL, run_dupcheck = T, run_parents = T, run_pedigree = T, ...){
  #===================sanity checks===============
  check.installed("sequoia")
  
  msg <- character(0)
  
  if(run_pedigree & !run_parents){
    msg <- c(msg, "Parents must be run before pedigree construction!\n")
  }

  if(run_pedigree){
    if(MaxSibIter <= 0 | MaxSibIter >= 25 ){ 
    msg <- c(msg, "Must include MaxSibIter value greater than 0 and less than 25 for pedigree construction!\n")
    }
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  #===================prep=========================
  # prep facets
  facets <- check.snpR.facet.request(x, facets, "none")
  x <- add.facets.snpR.data(x, facets)
  tasks <- get.task.list(x, facets = facets)
  
  # initialize
  out <- vector("list", nrow(tasks))
  
  #==================run===========================
  for(i in 1:nrow(tasks)){
    # subset the data
    tmatches <- fetch.sample.meta.matching.task.list(x, tasks[i,])
    snpmatches <- which(x@snp.meta[tasks[i, 3]] == tasks[i,4])
    tdat <- subset_snpR_data(x, samps = tmatches, snps = snpmatches)
   # filter snps for running sequoia - requires high MAF >0.25 and 50% inds otherwise inds dropped from seq
    tdat <- filter_snps(x = tdat, maf = 0.3, min_ind = 0.5, re_run = "full")
    
    # format the snps for sequoia
    tdat <- format_snps(x = tdat, output = "sequoia")
     
    # run sequoia
    if(run_dupcheck){
      sequoia(GenoM=tdat$dat, LifeHistData = tdat$lh, MaxSibIter = -1)
    }
    
    if(run_parents){
    p <- sequoia(GenoM=tdat$dat, LifeHistData = tdat$lh, MaxSibIter = 0, ...)
    }
    
    if(run_pedigree){
      sequoia(GenoM=tdat$dat, LifeHistData = tdat$lh, SeqList = p, MaxSibIter = 10, ...) #need to add so that someone can specify sibiter 0<25
      # need to add prior seqlist from parentage (p?)
    }
    
    tout <- sequoia(tdat, ...)
    out[[i]] <- tout
  
  return(out)
}