#' Run Sequoia pedigree/parentage assignment with snpR
#'
#' Runs the parentage assignment and pedigree construction tool from the
#' \code{sequoia} package. Note that this function \emph{is not overwrite
#' safe!}.
#'
#' This is a limited integration of the program and package written by Jisca 
#' Huisman. Note that there are many more Sequoia specific arguments that can 
#' be added to change from the default settings (eg. ErrorM, Tassign, Tfilt, 
#' etc.) See documentation for \code{sequoia}. These can be passed to the 
#' pedigree and parentage reconstructions using the ... argument in run_sequoia.
#' The sequoia package has many features, and snpR facilitates the use of a 
#' fraction of them. snpR users are encouraged to use the sequoia R package.
#'
#' @param x snpRdata object.
#' @param facets character, default NULL. Sample-specific facets over which the
#'   sequoia is called to run. See \code{\link{Facets_in_snpR}}.
#' @param run_dupcheck FALSE or TRUE, default FALSE. Uses sequoia to check for 
#'   duplicate samples in the dataset. Duplicate samples should not be included 
#'   for parentage and pedigree construction. 
#' @param run_parents FALSE or TRUE, default FALSE. Runs parentage assignments 
#'   for the samples. This runs quickly and is required before using the 
#'   run_pedigree command.
#' @param run_pedigree FALSE or TRUE, default FALSE. Runs pedigree construction 
#'   for the samples. This process can take a long time. 
#' @param run_relatives FALSE or TRUE, default FALSE. Runs retrieval of other
#'   relatives which did not pass thresholds for assignment in the main pedigree 
#'   construction model.
#' @param min_maf numeric in 0.25:0.5, default 0.3. Minimum allele frequency
#'   cutoff for analysis. Sequoia requires high minor allele frequencies for
#'   parentage and pedigree construction.
#' @param min_ind numeric in 0.5:1, default 0.5. Removes loci sequenced in less
#'   than this proportion of individuals. Note that \emph{individuals} with
#'   genotypes for fewer than half of the loci will be automatically removed by
#'   sequoia.
#' @param ... Additional arguments passed to\code{sequoia}
#'   (during parentage and pedigree reconstruction).
#'
#' @return A nested list with each facet specified containing sequoia output 
#'   summary information.
#'
#' @export
#' @author William Hemstrom
#' @author Melissa Jones
#'
#' @references Huisman,J. (2017) Pedigree reconstruction from SNP data:
#'   parentage assignment, sibship clustering and beyond. Mol. Ecol. Resour.,
#'   17, 1009–1024.
#'
#' @examples
#' # to follow an example using the stickSNPs example dataset you need 
#' # to add some variables that don't exist in the actual dataset.
#' a <- 2013:2015 #create a vector of possible birthyears
#' b <- c("M", "F", "U") #create a vector of possible sexes
#' stk <- stickSNPs
#' set.seed(4865)
#' sample.meta(stk)$BirthYear <- sample(x = a, size = nsamps(stickSNPs), 
#'                                      replace = TRUE) #create birth years
#' sample.meta(stk)$ID <- 1:nsamps(stk) #create unique sampleID
#' sample.meta(stk)$Sex <- sample(x= b, size = nsamps(stk), 
#'                                replace = TRUE) # create sexes
#' 
#' # slow, so not run here
#' \dontrun{
#' dup <- run_sequoia(x = stk, run_dupcheck = TRUE, run_parents = FALSE, 
#'                    run_pedigree = FALSE, run_relatives = FALSE)
#' ped <- run_sequoia(x = stk, run_dupcheck = FALSE, run_parents = TRUE, 
#'                    run_pedigree = TRUE, run_relatives = FALSE)
#' rel <- run_sequoia(x = stk, run_dupcheck = FALSE, run_parents = FALSE, 
#'                    run_pedigree = FALSE, run_relatives = TRUE)
#'                    
#' }
run_sequoia <- function(x, facets = NULL, run_dupcheck = FALSE, 
                                  run_parents = FALSE, run_pedigree = FALSE, 
                                  run_relatives = FALSE, min_maf = 0.3, 
                                  min_ind = 0.5, ...
                                  ){
  
  #===================sanity checks===============
  # check that provided snpRdata objects are in the correct format
  if(!snpR::is.snpRdata(x)){
    stop("Not a snpRdata object.\n")
  }
  .check.installed("sequoia")
  
  msg <- character(0)
  warn <- character(0)
  
  if(run_pedigree & !run_parents){
    msg <- c(msg, "Parents must be run before pedigree construction!\n")
  }
  
  if(min_maf < 0.25){
    warn <- c(warn, "Minor allele frequencies below 0.25 are not recommended for sequoia.\n")
  }
  if(min_ind < 0.5){
    warn <- c(warn, "Genotypes sequenced in fewer than 50% of individuals will automatically be removed by Sequoia.\n")
  }
  
  req.columns <- c("ID", "Sex", "BirthYear", "BYmin", "BYmax", "Yearlast")
  if(any(!req.columns %in% colnames(sample.meta(x)))){
    msg <- c(msg, "Columns named 'ID', 'Sex', 'BirthYear', 'BYmin', 'BYmax', 'Yearlast', are required in the sample metadata.\n")
  }
  
  if(length(warn) > 0){
    warning(warn)
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  #===================prep=========================
  # prep facets
  facets <- .check.snpR.facet.request(x, facets, "none")
  x <- .add.facets.snpR.data(x, facets)
  tasks <- .get.task.list(x, facets = facets)
  
  # initialize
  out <- vector("list", nrow(tasks))
  
  #==================run===========================
  for(i in 1:nrow(tasks)){
    
    names(out)[i] <- paste0(tasks[i,], collapse = "_")
    cat("Running facet: ", paste0(tasks[i,c(1,3)], collapse = " "), "\tsubfacet: ", paste0(tasks[i,c(2,4)], collapse = " "))
    
    
    # subset the data
    tmatches <- .fetch.sample.meta.matching.task.list(x, tasks[i,])
    if(tasks[i,3] == ".base"){
      snpmatches <- 1:nsnps(x)
    }
    else{
      snpmatches <- which(snp.meta(x)[,tasks[i, 3]] == tasks[i,4])
    }
    tdat <- x[snpmatches, tmatches]
    
    
    
    # filter snps for running sequoia - requires high MAF >0.25 and 50% inds otherwise inds dropped from seq
    tdat <- filter_snps(x = tdat, maf = min_maf, min_ind = min_ind, re_run = "full") # set as arguments later
    
    # format the snps for sequoia
    tdat <- format_snps(x = tdat, output = "sequoia")
    
    # initialize output
    out[[i]] <- vector("list")
    
    # run sequoia
    if(run_dupcheck){
      dups <- sequoia::sequoia(GenoM=tdat$dat, LifeHistData = tdat$lh, Module = "dup")
      out[[i]]$dups <- dups
    }
    
    if(run_parents){
      out[[i]]$parents <- sequoia::sequoia(GenoM=tdat$dat, LifeHistData = tdat$lh, Module = "par", ...)
    }
    
    if(run_pedigree){
      out[[i]]$pedigree <- sequoia::sequoia(GenoM=tdat$dat, LifeHistData = tdat$lh, SeqList = out[[i]]$parents, Module = "ped", ...) 
    }
    
    if(run_relatives){
      out[[i]]$relatives <- sequoia::GetMaybeRel(GenoM = tdat$dat, LifeHistData = tdat$lh, ...)
    }
  }
  return(
    out)
}


