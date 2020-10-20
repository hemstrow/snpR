#' Run Sequoia with snpR 
#' 
#' \code{run_sequoia} a tool for parentage assignment and pedigree construction. snpR integrates a program and package written by Jisca Huisman \itemize{\item{paper: }{https://doi.org/10.1111/1755-0998.12665} \item{github page: }{https://github.com/JiscaH/sequoia}} 
#' 
#' Possible options: \itemize{ \item{facets: }{Specification for grouping of data for which analyses should be ran independently} \item{run_dupcheck: }{Sequoia function used to check for duplicate samples in the dataset. It is reccomended to drop duplicated samples from the analysis. } \item{run_parents: }{Sequoia function to assign parents. This runs quickly and is required before running the run_pedigree command. } \item{run_pedigree: }{Sequoia function to construct full pedigree for the sample set. Sibship clustering takes a lot of time. }}
#' 
#' Note that there are many more Sequoia specific arguments that can be added to change from the default settings (eg. ErrorM, Tassign, Tfilt, GetMaybeRel, etc.)
#' 
#'@param x snpRdata object.
#'@param facets FALSE or character, default FALSE. Sample-specific facets over which the sequoia is called to run.
#'@param run_dupcheck FALSE or TRUE, default FALSE. Should a duplicate check be run on this dataset?
#'@param run_parents FALSE or TRUE, default FALSE. Should parentage assignments be run with the dataset?
#'@param run_pedigree FALSE or TRUE, default FALSE. Should a pedigree be constructed with the dataset? Requires run_parents to have been completed first. 
#'
#'@return A dataframe for each facet specified with sequoia output summary information. 
#'
#'@export
#'@author William Hemstrom
#'@author Melissa Jones
#'
#'@examples
#'to follow an example using the stickSNPs example dataset you need to add some variables that don't exist in the actual dataset. 
#'
#' #' stk <- stickSNPs
#' a <- 2013:2015 #create a vector of possible birthyears
#' b <- c("M", "F", "U") #create a vector of possible sexes
#' stk@sample.meta$BirthYear <- sample(x = a, size = nrow(stk@sample.meta), replace = T) #create birthyears
#' stk@sample.meta$ID <- 1:nrow(stk@sample.meta) #create unique sampleID
#' stk@sample.meta$Sex <- sample(x= b, size = nrow(stk@sample.meta), replace = T) # create sexes
#' dup <- run_sequoia(x=stk, run_dupcheck =T)


run_sequoia <- function(x, facets = NULL, run_dupcheck = T, run_parents = T, run_pedigree = T, pMaxSibIter = 10, ...){
  #===================sanity checks===============
  check.installed("sequoia")
  
  msg <- character(0)
  
  if(run_pedigree & !run_parents){
    msg <- c(msg, "Parents must be run before pedigree construction!\n")
  }

  if(run_pedigree){
    if(MaxSibIter < 0 | MaxSibIter >= 25 ){ 
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
      sequoia::sequoia(GenoM=tdat$dat, LifeHistData = tdat$lh, MaxSibIter = -1)
    }
    
    if(run_parents){
    p <- sequoia::sequoia(GenoM=tdat$dat, LifeHistData = tdat$lh, MaxSibIter = 0, ...)
    }
    
    if(run_pedigree){
      sequoia::sequoia(GenoM=tdat$dat, LifeHistData = tdat$lh, SeqList = p, MaxSibIter = pMaxSibIter, ...) 
    }
    
    tout <- sequoia(tdat, ...)
    out[[i]] <- tout
  
  return(out)
  }
}