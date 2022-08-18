#' Run Sequoia pedigree/parentage assignment with snpR
#'
#' Runs the parentage assignment and pedigree construction tool from the
#' \code{sequoia} package. Note that this function \emph{is not overwrite
#' safe!}. This function is a wrapper which sources code from github.
#'
#' This is an integration of the program and package written by Jisca Huisman.
#' Note that there are many more Sequoia specific arguments that can be added to
#' change from the default settings (eg. ErrorM, Tassign, Tfilt, GetMaybeRel,
#' etc.) See documentation for \code{sequoia}. These can be
#' passed to the pedigree and parentage reconstructions using the ... argument
#' in run_sequoia.
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
#' @param min_maf numeric in 0.25:0.5, default 0.3. Minimum allele frequency
#'   cutoff for analysis. Sequoia requires high minor allele frequencies for
#'   parentage and pedigree construction.
#' @param min_ind numeric in 0.5:1, default 0.5. Removes loci sequenced in less
#'   than this proportion of individuals. Note that \emph{individuals} with
#'   genotypes for fewer than half of the loci will be automatically removed by
#'   sequoia.
#' @param ask logical, default TRUE. Should the function ask for confirmation
#'   before sourcing github code?
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
#'                    run_pedigree = FALSE)
#' ped <- run_sequoia(x = stk, run_dupcheck = FALSE, run_parents = TRUE, 
#'                    run_pedigree = TRUE)
#' }
run_sequoia <- function(x, facets = NULL, run_dupcheck = FALSE, run_parents = FALSE, 
                        run_pedigree = FALSE, min_maf = 0.3, min_ind = 0.5, ask = TRUE, ...){
  
  run_sequoia_extension <- NULL
  
  if(ask){
    # confirm we want to run this
    cat("run_sequoia depends on the 'sequoia' package, which is not currently on CRAN.\nThis function is a wrapper that sources R scripts from github to pull in the run_sequoia function.\nIt is tested and should function normally.\nProceed?\t")
    cat("(y or n)\n")
    
    resp <- readLines(n = 1)
    resp <- tolower(resp)
    
    while(resp != "y"){
      cat("(y or n)\n")
      if(resp != "n"){
        return(FALSE)
      }
      else{
        resp <- readLines(n = 1)
        resp <- tolower(resp)
      }
    }
  }
  

  # source scripts and pull up internals
  .check.installed("devtools")
  
  devtools::source_url("https://raw.githubusercontent.com/hemstrow/snpR_extensions/main/run_sequoia.R")
  
  internals <- list(.add.facets.snpR.data = .add.facets.snpR.data,
                    .check.installed = .check.installed,
                    .check.snpR.facet.request = .check.snpR.facet.request, 
                    .split.facet = .split.facet, 
                    .tabulate_genotypes = .tabulate_genotypes, 
                    .make_it_quiet = .make_it_quiet,
                    .get.task.list = .get.task.list,
                    .fetch.sample.meta.matching.task.list = .fetch.sample.meta.matching.task.list,
                    .fetch.snp.meta.matching.task.list = .fetch.snp.meta.matching.task.list)
  
  # run and return
  out <- run_sequoia_extension(x = x, 
                               facets = facets, 
                               run_dupcheck = run_dupcheck, 
                               run_parents = run_parents, 
                               run_pedigree = run_pedigree, 
                               min_maf = min_maf, 
                               min_ind = min_ind,
                               ..., internals = internals)
  
  return(out)
}
