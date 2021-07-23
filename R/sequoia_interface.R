#' Run Sequoia pedigree/parentage assignment with snpR
#'
#' Runs the \code{\link[sequoia]{sequoia}} parentage assignment and pedigree
#' construction tool. Note that this function \emph{is not overwrite safe!}.
#'
#' This is an integration of the program and package written by Jisca Huisman.
#' Note that there are many more Sequoia specific arguments that can be added to
#' change from the default settings (eg. ErrorM, Tassign, Tfilt, GetMaybeRel,
#' etc.) See documentation for \code{\link[sequoia]{sequoia}}. These can be
#' passed to the pedigree and parentage reconstructions using the ... argument
#' to run_sequoia.
#'
#' @param x snpRdata object.
#' @param facets character, default NULL. Sample-specific facets over which the
#'   sequoia is called to run. See \code{\link{Facets_in_snpR}}.
#' @param run_dupcheck FALSE or TRUE, default FALSE. Determines if a Sequoia
#'   function used to check for duplicate samples in the dataset should be run
#'   prior to further analysis. This is generally recommended.
#' @param run_parents FALSE or TRUE, default FALSE. Determines if a Sequoia
#'   function to assign parents should be run. This runs quickly and is required
#'   before running the run_pedigree command.
#' @param run_pedigree FALSE or TRUE, default FALSE. Determines if a Sequoia
#'   function to construct full pedigree for the sample set should be run. This
#'   can take a while, since sibship clustering takes a lot of time for larger
#'   datasets.
#' @param pMaxSibIter numeric in 1:25, default 10. Maximum iterations to run for
#'   sibship clustering. Only specified for run_pedigree argument.
#' @param min_maf numeric in 0.25:0.5, default 0.3. Minimum allele frequency
#'   cutoff for analysis. Sequoia requires high minor allele frequencies for
#'   parentage and pedigree construction, so it is not generally recommended to
#'   reduce this number.
#' @param min_ind numeric in 0.5:1, default 0.5. Removes loci sequenced in less
#'   than this proportion of individuals. Note that \emph{individuals} with
#'   genotypes for less than half of the loci will be automatically removed by
#'   sequoia.
#' @param ... Additional arguments passed to \code{\link[sequoia]{sequoia}}
#'   (during parentage and pedigree reconstruction).
#'
#' @return A data.frame for each facet specified with sequoia output summary
#'   information.
#'
#' @author William Hemstrom
#' @author Melissa Jones
#'
#' @references Huisman,J. (2017) Pedigree reconstruction from SNP data:
#'   parentage assignment, sibship clustering and beyond. Mol. Ecol. Resour.,
#'   17, 1009–1024.
#'
#' @examples
#' # to follow an example using the stickSNPs example dataset you need to add some variables that don't exist in the actual dataset.
#' a <- 2013:2015 #create a vector of possible birthyears
#' b <- c("M", "F", "U") #create a vector of possible sexes
#' stk <- stickSNPs
#' sample.meta(stk)$BirthYear <- sample(x = a, size = nsamps(stickSNPs), replace = TRUE) #create birth years
#' sample.meta(stk)$ID <- 1:nsamps(stk) #create unique sampleID
#' sample.meta(stk)$Sex <- sample(x= b, size = nsamps(stk), replace = TRUE) # create sexes
#' dup <- run_sequoia(x = stk, run_dupcheck = TRUE, run_parents = FALSE, run_pedigree = FALSE)
run_sequoia <- function(x, facets = NULL, run_dupcheck = FALSE, run_parents = FALSE, run_pedigree = FALSE, pMaxSibIter = 10, min_maf = 0.3, min_ind = 0.5, ...){
  
  #===================sanity checks===============
  
  # check that provided snpRdata objects are in the correct format
    if(class(x) != "snpRdata"){
      stop("Not a snpRdata object.\n")
    }
  
  check.installed("sequoia")
  
  msg <- character(0)
  
  if(run_pedigree & !run_parents){
    msg <- c(msg, "Parents must be run before pedigree construction!\n")
  }

  if(run_pedigree){
    if(pMaxSibIter < 0 | pMaxSibIter >= 25 ){ 
      msg <- c(msg, "Must include MaxSibIter value greater than 0 and less than 25 for pedigree construction!\n")
    }
  }
  
  if(min_maf < 0.25){
    warning("Minumum minor allele frequencies below 0.25 are not recommended for sequoia.\n")
  }
  if(min_ind < 0.5){
    warning("Genotypes sequenced in less than 50% of individuals will automatically be removed by Sequoia.\n")
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
    
    names(out)[i] <- paste0(tasks[i,], collapse = "_")
    cat("Running facet: ", paste0(tasks[i,c(1,3)], collapse = " "), "\tsubfacet: ", paste0(tasks[i,c(2,4)], collapse = " "))
    
    
    # subset the data
    tmatches <- fetch.sample.meta.matching.task.list(x, tasks[i,])
    if(tasks[i,3] == ".base"){
      snpmatches <- 1:nsnps(x)
    }
    else{
      snpmatches <- which(snp.meta(x)[,tasks[i, 3]] == tasks[i,4])
    }
    tdat <- subset_snpR_data(x, samps = tmatches, snps = snpmatches)
    
    
    
    # filter snps for running sequoia - requires high MAF >0.25 and 50% inds otherwise inds dropped from seq
    tdat <- filter_snps(x = tdat, maf = min_maf, min_ind = min_ind, re_run = "full") # set as arguments later
    
    # format the snps for sequoia
    tdat <- format_snps(x = tdat, output = "sequoia")
    
    # initialize output
    out[[i]] <- vector("list")
    
    # run sequoia
    if(run_dupcheck){
      dups <- sequoia::sequoia(GenoM=tdat$dat, LifeHistData = tdat$lh, MaxSibIter = -1)
      out[[i]]$dups <- dups
    }
    
    if(run_parents){
      out[[i]]$parents <- sequoia::sequoia(GenoM=tdat$dat, LifeHistData = tdat$lh, MaxSibIter = 0, ...)
    }
    
    if(run_pedigree){
      out[[i]]$pedigree <- sequoia::sequoia(GenoM=tdat$dat, LifeHistData = tdat$lh, SeqList = out[[i]]$parents, MaxSibIter = pMaxSibIter, ...) 
    }
  }
  return(out)
}