run_sequoia <- function(x, facets = NULL, run_check = T, run_parents = T, run_pedigree = T, ...){
  #===================sanity checks===============
  check.installed("sequoia")
  
  msg <- character(0)
  
  if(run_pedigree & !run_parents){
    msg <- c(msg, "Parents must be run for a pedigree to be constructed!\n")
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
    
    # run sequoia
    if(run_check){
      
    }
    else{
      
    }
    
    
    tdat <- filter_snps()
    tdat <- format_snps()
    
    
    tout <- sequoia(tdat, ...)
    out[[i]] <- tout
  }
  
  return(out)
}