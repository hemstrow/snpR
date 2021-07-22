
#' Interfaces with the COLONY pedigree assignment program.
#'
#' Interfaces with the command-line version of the COLONY pedigree program.
#' Requires a snpRdata object containing offspring genotypes and can optionally
#' take snpRdata objects containing maternal or paternal genotypes, or both. No
#' facet support.
#'
#' Requires that the COLONY program is installed locally. Input and output files
#' will be stored in a colony folder created in the current working directory.
#' The functions documented here can write input files, call them using colony,
#' and parse some parts of the results given the original snpRdata object. Note
#' that no facet support is currently available here due to the complexity and
#' number of possible input parameters that are difficult to handle across
#' multiple facets and facet levels. Facet support for the basic, default
#' operation with few additional options may be added in the future. For now,
#' facets can be handled during parentage and pedigree creation in snpR using
#' the \code{\link{run_sequoia}} function, which runs a notably simpler (with
#' respect to implementation) pedigree toolkit.
#'
#' These functions includes many commonly used options but not all possible
#' parameters for colony inputs. See Colony User Guide for using extra features.
#'
#' This is still in development. The defaults and the most commonly used options
#' have been tested and work well, but some of the more esoteric options haven't
#' been fully tested yet.
#' @param x snpRdata object containing offspring genotypes.
#' @param outfile character, default "colony_input". Output file name. A file
#'   path may be provided (e.g. "colony/colony_run_1.txt").
#' @param method character, default "FPLS". Pedigree reconstruction method. For
#'   more details see the Colony User Guide. Options: \itemize{\item{"FLPS":
#'   }{Pure pairwise likelihood method, combines the full likelihood and
#'   pairwise likelihood methods. A good compromise between speed and
#'   accuracy.}{\item{"FL":}{Full Likelihood. More accurate than PLS but more
#'   computationally intensive and slow to run, especially with large complex
#'   datasets.}}\item{"PLS": }{Pairwise likelihood score. Less accurate but less
#'   computationally intensive than FL.}}
#' @param run_length numeric in c(1,2,3,4), default 2. Length of run:
#'   short/medium/long/verylong.
#' @param sampleIDs character, default NULL. Name of a column in the sample
#'   metadata that designates sample identifications/"names". Each name must be
#'   unique!
#' @param sibship_prior numeric in c(0, 1, 2, 3, 4), default 0. Strength the
#'   sibship size prior (no prior, weak, medium, strong, OR determine from known
#'   prior values). Values other than 0 require additional parameters. Option
#'   for supplying value of 4, currently not implemented in snpR. See Colony
#'   User Guide for more details.
#' @param paternal_sib_size numeric, default NULL. Minimum value is 0. The
#'   number of offspring in the candidate pool that are known to share the same
#'   father. If this value is not zero, then you must include a file with the
#'   known paternal sibship/paternity.
#' @param maternal_sib_size numeric, default NULL. Minimum value is 0. The
#'   number of offspring in the candidate pool that are known to share the same
#'   mother. If this value is not zero, then you must include a file with the
#'   known paternal sibship/maternity.
#' @param nruns integer, default 1. A number of replicate runs for the dataset.
#' @param seed integer, default NULL. Supply a four digit integer (eg: 1234,
#'   9876) as a starting point for the algorithm.
#' @param maternal_genotypes snpRdata object containing maternal genotypes.
#' @param paternal_genotypes snpRdata object containing paternal genotypes.
#' @param maternal_inclusion_prob numeric in 0:1, default 0. Probability the
#'   mother is in the dataset ranging from 0 to 1.
#' @param paternal_inclusion_prob numeric in 0:1, default 0. Probability the
#'   father is in the dataset ranging from 0 to 1.
#' @param update_af character, default TRUE. Should Colony update the allele
#'   frequencies used in the calculations?
#' @param dioecious character, default TRUE. Is this species diploid/dioecious?
#'   FALSE = haploid/monoecious. Colony does not work with more ploidy.
#' @param inbreeding character, default TRUE. Should Colony assume inbreeding in
#'   the calculations?
#' @param male_monogamous character, default FALSE. Should Colony assume males
#'   are monogamous?
#' @param female_monogamous character, default FALSE. Should Colony assume
#'   females are monogamous?
#' @param clone_inference character, default FALSE. Should Colony infer clones
#'   in the sample set?
#' @param sibship_scaling character, default TRUE. Should Colony scale sibling
#'   groups?
#' @param known_af character, default FALSE. If TRUE snpR will calculate and
#'   supply mafs, else, supply a numeric vector containing the known maf for
#'   each locus.
#' @param precision integer in c(1,2,3,4), default 2. Low/Medium/High/Very High
#'   for calculating the maximum likelihood.
#' @param dropout numeric vector where each value is in 0:1, default 0.01.
#'   Supply a flatrate value for all markers, or a vector corresponding to the
#'   allelic droput rate for each marker.
#' @param genotyping_error numeric vector where each value is in 0:1, default
#'   0.01. Supply a flatrate value for all markers, or a vector corresponding to
#'   the genotyping error rate for each marker.
#' @param known_maternal_dyads character, default NULL. Supply matrix or
#'   dataframe with known maternal-offspring dyads. Offspring ID in column 1,
#'   Maternal ID in column 2.
#' @param known_paternal_dyads character, default NULL. Supply matrix or
#'   dataframe with known paternal-offspring dyads. Offspring ID in column 1,
#'   Paternal ID in column 2.
#' @param known_maternal_max_mismatches integer in c(0,1,2:nsample), default 0.
#' @param known_paternal_max_mismatches integer in c(0,1,2:nsample), default 0.
#' @param known_maternal_sibships character, default NULL. Data frame or matrix
#'   with sibship size followed by single column containing all of the sibling
#'   IDs separated by spaces.
#' @param known_paternal_sibships character, default NULL. Data frame or matrix
#'   with sibship size followed by single column containing all of the sibling
#'   IDs separated by spaces.
#' @param maternal_exclusions character, default NULL. Data.frame or matrix with
#'   column 1 the offspring ID, column 2 the number of excluded females, column
#'   3 the IDs of excluded females separated by spaces.
#' @param paternal_exclusions character, default NULL. Data.frame or matrix with
#'   column 1 the offspring ID, column 2 the number of excluded males, column 3
#'   the IDs of excluded males separated by spaces.
#' @param excluded_maternal_siblings character, default NULL. Data.frame or
#'   matrix with column 1 the offspring ID, column 2 the number of excluded
#'   siblings, column 3 the IDs of excluded siblings separated by spaces.
#' @param excluded_paternal_siblings character, default NULL. Data.frame or
#'   matrix with column 1 the offspring ID, column 2 the number of excluded
#'   siblings, column 3 the IDs of excluded siblings separated by spaces.
#' @param infile character. Path to the pre-written colony input file to be run.
#' @param colony_path character. Path to the colony executable.
#' @param path character. Path to the directory containing colony results.
#' @param x snpRdata object from which metadata for colony results can be found.
#' @param prefix character. The prefix for the colony files to be parsed.
#'
#' @author William Hemstrom
#' @author Melissa Jones
#'
#' @aliases write_colony_input call_colony parse_colony run_colony
#'
#' @name colony_interface
#'
#' @references Jones,O.R. and Wang,J. (2010) COLONY: a program for parentage and
#'   sibship inference from multilocus genotype data. Mol. Ecol. Resour., 10,
#'   551–555.
#'
#' @examples
#' 
#' # A simple example for running all individuals in the snpR object as siblings in
#' # colony. Not run to avoid clutter.
#' \dontrun{
#'   write_colony_input(x = stickSNPs, outfile = "stk.col")
#'  }
#'   
#' # A more complex example requires 1) creating and adding variables to the 
#' # stickSNPs example dataset and 2) creating subset snpR objects.
#' a <- 2013:2015 #create a vector of possible birthyears
#' b <- c("M", "F", "U") #create a vector of possible sexes
#' stk <- stickSNPs
#' sample.meta(stk)$BirthYear <- sample(x = a, size = nsamps(stk), replace = TRUE) #create birthyears
#' sample.meta(stk)$ID <- 1:nsamps(stk) #create unique sampleID
#' sample.meta(stk)$Sex <- sample(x= b, size = nsamps(stk), replace = TRUE) # create sexes
#' #generating snpR objects for male and female potential parents and offspring 
#' # (no U sexes in the potential parents in for this example)
#' sir <- subset_snpR_data(stk, Sex = "M", BirthYear = 2013)
#' dam <- subset_snpR_data(stk, Sex = "F", BirthYear = 2013)
#' off <- subset_snpR_data(stk, BirthYear = c("2014", "2015"))
#' # not run to avoid clutter
#' \dontrun{
#'   write_colony_input(x = off, outfile = "parents_example.col", 
#'                      maternal_genotypes = dam, paternal_genotypes = sir)
#' }
#' 
#' # running a simple model
#' \dontrun{
#'   ## intentionally shorter run, with a small subset of the samples
#'   test_dat <- subset_snpR_data(stickSNPs, pop = "ASP")
#'   run_colony(x = test_dat, colony_path = "/usr/bin/colony2s.exe", method = "PLS", run_length = 1) # intentionally a short, quick and dirty run.
#' }
NULL


#' @describeIn colony_interface Create a colony input file
#' @export
write_colony_input <- function(x, outfile = "colony_input", method = "FPLS", run_length = 2, sampleIDs = NULL,
                               sibship_prior = 0, paternal_sib_size = NULL, maternal_sib_size = NULL,
                               nruns = 1, seed = NULL, maternal_genotypes = NULL, paternal_genotypes = NULL,
                               maternal_inclusion_prob = 0, paternal_inclusion_prob = 0,
                               update_af = TRUE, dioecious = TRUE, inbreeding = TRUE, male_monogamous = FALSE,
                               female_monogamous = FALSE, clone_inference = FALSE, sibship_scaling = TRUE, known_af = FALSE,
                               precision = 2, dropout = 0.01, genotyping_error = 0.01,
                               known_maternal_dyads = NULL, known_paternal_dyads = NULL,
                               known_maternal_max_mismatches = 0, known_paternal_max_mismatches = 0,
                               known_maternal_sibships = NULL, known_paternal_sibships = NULL,
                               maternal_exclusions = NULL, paternal_exclusions = NULL,
                               excluded_maternal_siblings = NULL, excluded_paternal_siblings = NULL){


  #=====================initialize===============
  # initialize storage directory
  if(basename(getwd()) != "colony"){
    if(!dir.exists("colony")){
      dir.create("colony")
      setwd("colony")
    }
    else{
      setwd("colony")
    }
  }

  # write anything preliminary to output
  proj <- outfile
  if(nchar(proj) > 40){
    proj <- substr(proj, 1, 40)
  }
  outfile <- paste0(proj, ".dat")

  #=====================sanity checks============
  # check if file exists
  if(file.exists(outfile)){
    file.remove(outfile)
  }

  #=====================write first stuff=============

  # project name
  write(proj, outfile)

  # outfile
  write(proj, outfile, append = T)

  # other info:
  write(ncol(x), outfile, append = T) # sample size
  write(nrow(x), outfile, append = T) # nloci
  if(is.null(seed)){
    seed <- sample(1000000, 1)
  }
  write(seed, outfile, append = T) # seed
  write(as.numeric(update_af), outfile, append = T) # update allele frequency in bayesian model
  if(dioecious){
    dioecious <- 2
  }
  else{
    dioecious <- 1
  }
  write(dioecious, outfile, append = T) # dio or monoecious
  write(as.numeric(inbreeding), outfile, append = T) # inbreeding
  write(0, outfile, append = T) # diploid only for now
  write(paste(as.numeric(male_monogamous), as.numeric(female_monogamous), collapse = " "), outfile, append = T) # polygamy
  write(as.numeric(clone_inference), outfile, append = T) # infer clones
  write(as.numeric(sibship_scaling), outfile, append = T) # scale sibships

  # sibship scaling priors
  sib.prior <- sibship_prior
  if(!is.null(paternal_sib_size)){
    sib.prior <- paste(sib.prior, paternal_sib_size, maternal_sib_size)
  }
  write(sib.prior, outfile, append = T)


  #=====================reformat genotypes, write info on them==================
  colony_genotypes <- format_snps(x, "colony")
  if(!is.null(sampleIDs)){
    colony_genotypes[,1] <- x@sample.meta[,sampleIDs]
  }
  # allele frequencies
  ## if provided with minor allele frequencies
  if(length(known_af) == nrow(x) & is.numeric(known_af)){
    write(1, outfile, append = T)
    write(rep(2, nrow(x)), outfile, append = T, sep = " ", ncolumns = nrow(x)) # number of alleles per locus, should all be two.
    maj <- 1 - known_af
    afs <- cbind(maj, known_af)
    afs_tab <- matrix(0, nrow = 2*nrow(afs), ncol = 2)
    afs_tab[seq(2, nrow(afs_tab), by = 2),] <- afs
    afs_tab[seq(1, nrow(afs_tab), by = 2),] <- c(rep(1, nrow(afs)), rep(2, nrow(afs)))
    write.table(afs_tab, outfile, sep = " ", append = T, row.names = F, col.names = F)
  }
  # otherwise, check the logical and calculate minor allele frequencies if true
  else if(is.logical(known_af[1]) & length(known_af) == 1){
    write(as.numeric(known_af), outfile, append = T)
    if(known_af){
      write(rep(2, nrow(x)), outfile, append = T, sep = " ", ncolumns = nrow(x)) # number of alleles per locus, should all be two.

      # afs
      x <- calc_maf(x)
      afs <- x@stats[x@stats$facet == ".base",]$maf
      afs <- cbind(1-afs, afs)
      afs_tab <- matrix(0, nrow = 2*nrow(afs), ncol = 2)
      afs_tab[seq(2, nrow(afs_tab), by = 2),] <- afs
      afs_tab[seq(1, nrow(afs_tab), by = 2),] <- c(rep(1, nrow(afs)), rep(2, nrow(afs)))
      write.table(afs_tab, outfile, sep = " ", append = T, row.names = F, col.names = F)
    }
  }
  else{
    file.remove(outfile)
    stop("known_af must either be TRUE, FALSE, or a vector of minor allele frequencies.\n")
  }


  write(nruns, outfile, append = T) # nruns
  write(run_length, outfile, append = T) # run length
  write(0, outfile, append = T) # monitor method
  write(10000, outfile, append = T) # monitor interval
  write(0, outfile, append = T) # no GUI

  # analysis method
  analysis_methods <- c("PLS", "FL", "FPLS")
  method <- match(method, analysis_methods) - 1
  write(method, outfile, append = T)

  write(precision, outfile, append = T) # precision
  write("MK@", outfile, append = T) # marker names
  write("0@", outfile, append = T) # marker types

  # genotyping error and dropout
  ## dropout
  if(length(dropout) == 1){
    write(rep(dropout, nrow(x)), outfile, append = T, ncolumns = nrow(x), sep = " ")
  }
  else{
    write(dropout, outfile, append = T, ncolumns = nrow(x), sep = " ")
  }
  ## error rates
  if(length(genotyping_error) == 1){
    write(rep(genotyping_error, nrow(x)), outfile, append = T, ncolumns = nrow(x), sep = " ")
  }
  else{
    write(genotyping_error, outfile, append = T, ncolumns = nrow(x), sep = " ")
  }
  # offspring genotypes
  write.table(colony_genotypes, outfile, T, quote = F, sep = " ", row.names = F, col.names = F)

  #===================maternal and paternal genotype things==============
  # number of male and female candidates, maternal and paternal inclusion probabilities in OMS and PMS

  #there needs to be a newline between the offspring genotypes and the parent data #eg
      # offspring genotype end
      # "\n"
      # 0 0 # M P
      # 0 0 # M P
      # "\n"
    ## then a newline before first parental genotypes
  write("\n", outfile, append = T) # adding a newline between offspring genotypes and parents
  par_nums <- numeric(2)
  if(class(paternal_genotypes) == "snpRdata"){

    par_nums[1] <- ncol(paternal_genotypes) #ncol not nrow!
  }
  else{
    paternal_inclusion_prob <- 0
  }
  if(class(maternal_genotypes) == "snpRdata"){
    par_nums[2] <- ncol(maternal_genotypes) #ncol not nrow!
  }
  else{
    maternal_inclusion_prob <- 0
  }
  write(c(paternal_inclusion_prob, maternal_inclusion_prob), outfile, append = T, sep = " ") # write inclusion probs
  write(par_nums, outfile, append = T, sep = " ") # OPS and OMS sample sizes
  write("\n", outfile, append = T) #added newline after number parents and probs parents
  # paternal and maternal genotypes
  if(class(paternal_genotypes) == "snpRdata"){
    male_colony <- format_snps(paternal_genotypes, "colony")
    write.table(male_colony, outfile, T, quote = F, sep = " ", row.names = F, col.names = F)
    write("\n", outfile, append = T) #added newline for separating genotypes
  }
  write("\n", outfile, append = T) #added newline

  if(class(maternal_genotypes) == "snpRdata"){
    female_colony <- format_snps(maternal_genotypes, "colony")
    write.table(female_colony, outfile, T, quote = F, sep = " ", row.names = F, col.names = F)
    write("\n", outfile, append = T) #added newline for separating stuff
  }

  # paternal dyads:
  if(!is.null(known_paternal_dyads[1,1])){
    n_p_d <- nrow(known_paternal_dyads)
    write(c(n_p_d, known_paternal_max_mismatches), outfile, append = T, sep = " ") # write the number known + mismatch max
    write.table(known_paternal_dyads, outfile, T, quote = F, sep = " ", row.names = F, col.names = F) # write dyads. They should be a matrix or data frame, with the offspring ID in the first column and the paternal ID in the second
  }
  else{
    write(c(0, 0), outfile, append =  T, sep = " ") # no dyads
  }

  # maternal dyads:
  if(!is.null(known_maternal_dyads[1,1])){
    n_m_d <- nrow(known_maternal_dyads)
    write(c(n_m_d, known_maternal_max_mismatches), outfile, append = T, sep = " ") # write the number known + mismatch max
    write.table(known_maternal_dyads, outfile, T, quote = F, sep = " ", row.names = F, col.names = F) # write dyads. They should be a matrix or data frame, with the offspring ID in the first column and the maternal ID in the second
  }
  else{
    write(c(0, 0), outfile, append =  T, sep = " ") # no dyads
  }

  # known paternal sibships
  if(!is.null(known_paternal_sibships[1,1])){
    n_p_s <- nrow(known_paternal_sibships)
    write(n_p_s, outfile, append = T) # write the number known
    write.table(known_paternal_sibships, outfile, T, quote = F, sep = " ", row.names = F, col.names = F) # write sibsips. Data frame or matrix with sibship size followed by single column containing all of the sibling IDs.
  }
  else{
    write(0, outfile, append =  T) # no sibships
  }

  # known maternal sibships
  if(!is.null(known_maternal_sibships[1,1])){
    n_m_s <- nrow(known_maternal_sibships)
    write(n_m_s, outfile, append = T) # write the number known
    write.table(known_maternal_sibships, outfile, T, quote = F, sep = " ", row.names = F, col.names = F) # write sibsips. Data frame or matrix with sibship size followed by single column containing all of the sibling IDs.
  }
  else{
    write(0, outfile, append =  T) # no sibships
  }

  # paternal exclusions
  if(!is.null(paternal_exclusions[1,1])){
    npe <- nrow(paternal_exclusions)
    write(npe, outfile, append = T) # write number of exclusions
    write.table(paternal_exclusions, outfile, T, quote = F, sep = " ", row.names = F, col.names = F) # write the paternal exclusions. Data.frame or matrix with col one the offspring ID, col 2 the number of excluded males, col 3 the IDs of excluded males separated by spaces.
  }
  else{
    write(0, outfile, append = T) # no exclusions
  }

  # maternal exclusions
  if(!is.null(maternal_exclusions[1,1])){
    npe <- nrow(maternal_exclusions)
    write(npe, outfile, append = T) # write number of exclusions
    write.table(maternal_exclusions, outfile, T, quote = F, sep = " ", row.names = F, col.names = F) # write the maternal exclusions. Data.frame or matrix with col one the offspring ID, col 2 the number of excluded males, col 3 the IDs of excluded females separated by spaces.
  }
  else{
    write(0, outfile, append = T) # no exclusions
  }

  # excluded paternal sibships
  if(!is.null(excluded_paternal_siblings[1,1])){
    n_eps <- nrow(excluded_paternal_siblings)
    write(n_eps, outfile, append = T) # number of exclusions
    write.table(excluded_paternal_siblings, outfile, T, quote = F, sep = " ", row.names = F, col.names = F) # write the excluded paternal sibships. Data.frame or matrix with col one the offspring ID, col 2 the number of excluded siblings, col 3 the IDs of excluded siblings separated by spaces.
  }
  else{
    write(0, outfile, append = T) # no exclusions
  }

  # excluded maternal sibships
  if(!is.null(excluded_maternal_siblings[1,1])){
    n_eps <- nrow(excluded_maternal_siblings)
    write(n_eps, outfile, append = T) # number of exclusions
    write.table(excluded_maternal_siblings, outfile, T, quote = F, sep = " ", row.names = F, col.names = F) # write the excluded maternal sibships. Data.frame or matrix with col one the offspring ID, col 2 the number of excluded siblings, col 3 the IDs of excluded siblings separated by spaces.
  }
  else{
    write(0, outfile, append = T) # no exclusions
  }

  #=================return to original wd===========
  setwd("..")
  return(paste0("Colony input saved to: ./colony/", outfile, collapse = ""))
}


#' Call a prepared colony infile.
#' @describeIn colony_interface Call a colony executable to run a prepared colony input file.
#' @export
call_colony <- function(infile, colony_path){

  # check system type
  sys.type <- Sys.info()["sysname"]

  # save the call
  call <- paste0(colony_path," IFN:", infile)

  # call
  if(sys.type == "Windows"){
    shell(call)
  }
  else{
    system(call)
  }

  # move results
  pattern <- gsub(".+/", "", infile)
  pattern <- gsub("\\..+$", "", pattern)
  files <- list.files(pattern = pattern)

  if(!dir.exists("colony")){
    dir.create("colony")
  }
  
  file.rename(files, paste0("./colony/", files))
}

#' Parse colony data.
#' @describeIn colony_interface Parse a previously run colony analysis.
#' @export
parse_colony <- function(prefix, x, path = "./colony/", sampleIDs = NULL){
  #===============full and half sibs==============================
  # read in half and full sibs
  fsd <- readr::read_csv(paste0(path, prefix, ".FullSibDyad"))
  if(file.exists(paste0(path, prefix, ".HalfSibDyad"))){
    hsd <- readr::read_csv(paste0(path, prefix, ".HalfSibDyad"))
    ## combine
    if(nrow(fsd) > 0 & nrow(hsd) > 0){
      dyads <- rbind(cbind(fsd, type = "FullSib", stringsAsFactors = F),
                     cbind(hsd, type = "HalfSib", stringsAsFactors = F))
    }
    else if(nrow(fsd) > 0){
      dyads <- cbind(fsd, type = "FullSib", stringsAsFactors = F)
    }
    else if(nrow(hsd) > 0){
      dyads <- cbind(hsd, type = "HalfSib", stringsAsFactors = F)
    }
    
  }
  else{
    dyads <- cbind(fsd, type = "FullSib", stringsAsFactors = F)
  }

  colnames(dyads)[1:2] <- c("Sample1", "Sample2")

  ## add in non-sibs
  all_pairs <- as.data.frame(t(utils::combn(colnames(x), 2)), stringsAsFactors = F)
  colnames(all_pairs) <- c("Sample1", "Sample2")
  cdyads <- rbind(dyads, dyads[,c(2,1,3,4)])
  all_pairs <- merge(all_pairs, cdyads, by = c("Sample1", "Sample2"), all.x = T)
  all_pairs$type[which(is.na(all_pairs$type))] <- "none"

  #=============add male and female parent ID info to the sample meta===========
  clusters <- readr::read_table2(paste0(path, prefix, ".BestCluster"))
  colnames(clusters)[2] <- "ClusterProbability"
  if(is.null(sampleIDs)){
    x@sample.meta$.offspringID <- colnames(x)
    x@sample.meta <- merge(x@sample.meta, clusters, by.x = ".offspringID", by.y = "OffspringID", sort = F)
    x@sample.meta$.offspringID <- NULL
  }
  else{
    x@sample.meta <- merge(x@sample.meta, clusters, by.x = sampleIDs, by.y = "OffspringID", sort = F)
  }

  .sidcol <- which(colnames(x@sample.meta) == ".sample.id")
  x@sample.meta <- x@sample.meta[,c((1:ncol(x@sample.meta))[-.sidcol], .sidcol)]

  return(list(x = x, dyads = dyads, all_pairs = all_pairs, clusters = clusters))
}

#' Create and infile, run, and gather some basic results using the COLONY parentage analysis program.
#'
#' Create a COLONY infile using snpRdata sets containg offspring and possibly paternal genotypes given
#' specified parameters. Requires that the COLONY program is installed locally. Input and output files
#' will be stored in a colony folder created in the current working directory.
#'
#' @export
#' @describeIn colony_interface run colony on a snpRdata object, start to finish
run_colony <- function(x, colony_path, outfile = "colony_input", method = "FPLS", run_length = 2, sampleIDs = NULL,
                       sibship_prior = 0, paternal_sib_size = NULL, maternal_sib_size = NULL,
                       nruns = 1, seed = NULL, maternal_genotypes = NULL, paternal_genotypes = NULL,
                       maternal_inclusion_prob = 0, paternal_inclusion_prob = 0,
                       update_af = TRUE, dioecious = TRUE, inbreeding = TRUE, male_monogamous = FALSE,
                       female_monogamous = FALSE, clone_inference = FALSE, sibship_scaling = TRUE, known_af = FALSE,
                       precision = 2, dropout = 0.01, genotyping_error = 0.01,
                       known_maternal_dyads = NULL, known_paternal_dyads = NULL,
                       known_maternal_max_mismatches = 0, known_paternal_max_mismatches = 0,
                       known_maternal_sibships = NULL, known_paternal_sibships = NULL,
                       maternal_exclusions = NULL, paternal_exclusions = NULL,
                       excluded_maternal_siblings = NULL, excluded_paternal_siblings = NULL){
  msg <- character(0)
  if(!file.exists(colony_path)){
    msg <- c(msg, "Cannot find colony executable.\n")
  }
  if(!is.snpRdata(x)){
    msg <- c(msg, "x is not a snpRdata object.\n")
  }
  
  check.installed("readr")
  
  if(length(msg) > 0){
    stop(msg)
  }

  # prep
  write_colony_input(x = x,
                     outfile = outfile,
                     method = method,
                     run_length = run_length,
                     sampleIDs = sampleIDs,
                     sibship_prior = sibship_prior,
                     paternal_sib_size = paternal_sib_size,
                     maternal_sib_size = maternal_sib_size, nruns = nruns, seed = seed,
                     maternal_genotypes = maternal_genotypes,
                     paternal_genotypes = paternal_genotypes,
                     maternal_inclusion_prob = maternal_inclusion_prob,
                     paternal_inclusion_prob = paternal_inclusion_prob,
                     update_af = update_af,
                     dioecious = dioecious,
                     inbreeding = inbreeding,
                     male_monogamous = male_monogamous,
                     female_monogamous = female_monogamous,
                     clone_inference = clone_inference,
                     sibship_scaling = sibship_scaling,
                     known_af = known_af,
                     precision = precision,
                     dropout = dropout,
                     genotyping_error = genotyping_error,
                     known_maternal_dyads = known_maternal_dyads,
                     known_paternal_dyads = known_paternal_dyads,
                     known_maternal_max_mismatches = known_maternal_max_mismatches,
                     known_paternal_max_mismatches = known_paternal_max_mismatches,
                     known_maternal_sibships = known_maternal_sibships,
                     known_paternal_sibships = known_paternal_sibships,
                     maternal_exclusions = maternal_exclusions,
                     paternal_exclusions = paternal_exclusions,
                     excluded_maternal_siblings = excluded_maternal_siblings,
                     excluded_paternal_siblings = excluded_paternal_siblings)

  # run
  call_colony(infile = paste0("./colony/", outfile, ".dat"),
              colony_path = colony_path)

  # parse some basics
  return(parse_colony(prefix = outfile,
                      x = x,
                      path = "./colony/",
                      sampleIDs = sampleIDs))

}
