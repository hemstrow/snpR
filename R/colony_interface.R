
#' Write an input for the COLONY pedigree assignment program.
write_colony_input <- function(x, outfile = "colony_input", method = "FPLS", run_length = 2,
                               sibship_prior = 0, paternal_sib_size = NULL, maternal_sib_size = NULL,
                               nruns = 1, seed = NULL, maternal_genotypes = NULL, paternal_genotypes = NULL,
                               maternal_inclusion_prob = 0, paternal_inclusion_prob = 0,
                               update_af = TRUE, dioecious = TRUE, inbreeding = TRUE, male_monogamous = FALSE,
                               female_monogamous = FALSE, clone_inference = FALSE, sibship_scaling = TRUE, known_af = FALSE,
                               precision = 2, dropout = 0.01, genotyping_error = 0.01,
                               known_maternal_dryads = NULL, known_paternal_dryads = NULL,
                               known_maternal_max_mismatches = 0, known_paternal_max_mismatches = 0,
                               known_maternal_sibships = NULL, known_paternal_sibships = NULL,
                               maternal_exclusions = NULL, paternal_exclusions = NULL,
                               excluded_maternal_siblings = NULL, excluded_paternal_siblings = NULL){

  #=====================initialize===============
  # initialize storage directory
  if(!dir.exists("colony")){
    dir.create("colony")
    setwd("colony")
  }
  else{
    setwd("colony")
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
    cat("COLONY input file", outfile,  "already exits. ")
    resp <- "empty"
    while(resp != "y" & resp != "n"){
      cat("Overwrite? (y or n)\n")
      resp <- readLines(n = 1)
    }
    if(resp == "n"){
      setwd("..")
      stop("Please move or rename existing input file in the colony directory")
    }
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


  # allele frequencies
  write(as.numeric(known_af), outfile, append = T)
  if(known_af){
    write(rep(2, nrow(x)), outfile, append = T, sep = " ", ncolumns = nrow(x)) # number of alleles per locus, should all be two.

    # afs
    x <- calc_maf(x)
    afs <- x@stats[x@stats$facet == ".base",]$maf
    afs <- cbind(1-afs, afs)
    for(i in 1:length(afs)){
      write(1:2, outfile, sep = " ", append = T)
      write(afs[i,], outfile, append = T, sep = " ")
    }
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
  par_nums <- numeric(2)
  if(class(paternal_genotypes) == "snpRdata"){
    par_nums[1] <- nrow(paternal_genotypes)
  }
  else{
    paternal_inclusion_prob <- 0
  }
  if(class(maternal_genotypes) == "snpRdata"){
    par_nums[2] <- nrow(maternal_genotypes)
  }
  else{
    maternal_inclusion_prob <- 0
  }
  write(c(paternal_inclusion_prob, maternal_inclusion_prob), outfile, append = T, sep = " ") # write inclusion probs
  write(par_nums, outfile, append = T, sep = " ") # OPS and OMS sample sizes

  # paternal and maternal genotypes
  if(class(paternal_genotypes) == "snpRdata"){
    male_colony <- format_snps(paternal_genotypes, "colony")
    write.table(male_colony, outfile, T, quote = F, sep = " ", row.names = F, col.names = F)
  }
  if(class(maternal_genotypes) == "snpRdata"){
    female_colony <- format_snps(maternal_genotypes, "colony")
    write.table(female_colony, outfile, T, quote = F, sep = " ", row.names = F, col.names = F)
  }

  # paternal dryads:
  if(!is.null(known_paternal_dryads[1,1])){
    n_p_d <- nrow(known_paternal_dryads)
    write(c(n_p_d, known_paternal_max_mismatches), outfile, append = T, sep = " ") # write the number known + mismatch max
    write.table(known_paternal_dryads, outfile, T, quote = F, sep = " ", row.names = F, col.names = F) # write dryads. They should be a matrix or data frame, with the offspring ID in the first column and the paternal ID in the second
  }
  else{
    write(c(0, 0), outfile, append =  T, sep = " ") # no dryads
  }

  # maternal dryads:
  if(!is.null(known_maternal_dryads[1,1])){
    n_m_d <- nrow(known_maternal_dryads)
    write(c(n_m_d, known_maternal_max_mismatches), outfile, append = T, sep = " ") # write the number known + mismatch max
    write.table(known_maternal_dryads, outfile, T, quote = F, sep = " ", row.names = F, col.names = F) # write dryads. They should be a matrix or data frame, with the offspring ID in the first column and the maternal ID in the second
  }
  else{
    write(c(0, 0), outfile, append =  T, sep = " ") # no dryads
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
    write.table(paternal_exclusions, outfile, T, quote = F, sep = " ", row.names = F, col.names = F) # write the paternal exclusions. Data.frame or matrix with col one the offspring ID, col 2 the number of excluded males, col 3 the IDs of excluded males seperated by spaces.
  }
  else{
    write(0, outfile, append = T) # no exclusions
  }

  # maternal exclusions
  if(!is.null(maternal_exclusions[1,1])){
    npe <- nrow(maternal_exclusions)
    write(npe, outfile, append = T) # write number of exclusions
    write.table(maternal_exclusions, outfile, T, quote = F, sep = " ", row.names = F, col.names = F) # write the maternal exclusions. Data.frame or matrix with col one the offspring ID, col 2 the number of excluded males, col 3 the IDs of excluded females seperated by spaces.
  }
  else{
    write(0, outfile, append = T) # no exclusions
  }

  # excluded paternal sibships
  if(!is.null(excluded_paternal_siblings[1,1])){
    n_eps <- nrow(excluded_paternal_siblings)
    write(n_eps, outfile, append = T) # number of exclusions
    write.table(excluded_paternal_siblings, outfile, T, quote = F, sep = " ", row.names = F, col.names = F) # write the exluded paternal sibships. Data.frame or matrix with col one the offspring ID, col 2 the number of excluded siblings, col 3 the IDs of excluded siblings seperated by spaces.
  }
  else{
    write(0, outfile, append = T) # no exclusions
  }

  # excluded maternal sibships
  if(!is.null(excluded_maternal_siblings[1,1])){
    n_eps <- nrow(excluded_maternal_siblings)
    write(n_eps, outfile, append = T) # number of exclusions
    write.table(excluded_maternal_siblings, outfile, T, quote = F, sep = " ", row.names = F, col.names = F) # write the exluded maternal sibships. Data.frame or matrix with col one the offspring ID, col 2 the number of excluded siblings, col 3 the IDs of excluded siblings seperated by spaces.
  }
  else{
    write(0, outfile, append = T) # no exclusions
  }

  #=================return to original wd===========
  setwd("..")
}


#' Call a prepared colony infile.
call_colony <- function(infile, colony_path = "."){

  # adjust inputs to fix hanging /
  if(!grepl("/$", colony_path)){
    colony_path <- paste0(colony_path, "/")
  }

  # check system type
  sys.type <- Sys.info()["sysname"]

  # name the program
  if(sys.type == "Windows"){
    colony_program <- "Colony2P.exe"
  }
  else{
    colony_program <- "colony2s.ifort.out"
  }


  # save the call
  call <- paste0(colony_path, colony_program, " IFN:", infile)

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
parse_colony <- function(prefix, x, path = "./colony/"){
  #===============full and half sibs==============================
  # read in half and full sibs
  fsd <- readr::read_csv(paste0(path, prefix, ".FullSibDyad"))
  if(file.exists(paste0(path, prefix, ".HalfSibDyad"))){
    hsd <- readr::read_csv(paste0(path, prefix, ".HalfSibDyad"))
    ## combine
    dyads <- rbind(cbind(fsd, type = "FullSib", stringsAsFactors = F),
                   cbind(hsd, type = "HalfSib", stringsAsFactors = F))
  }
  else{
    dyads <- cbind(fsd, type = "FullSib", stringsAsFactors = F)
  }

  colnames(dyads)[1:2] <- c("Sample1", "Sample2")

  ## add in non-sibs
  all_pairs <- as.data.frame(t(combn(colnames(x), 2)), stringsAsFactors = F)
  colnames(all_pairs) <- c("Sample1", "Sample2")
  cdyads <- rbind(dyads, dyads[,c(2,1,3,4)])
  all_pairs <- merge(all_pairs, cdyads, by = c("Sample1", "Sample2"), all.x = T)
  all_pairs$type[which(is.na(all_pairs$type))] <- "none"

  return(list(dyads = dyads, all_pairs = all_pairs))
}

#' Create and infile, run, and gather some basic results using the COLONY parentage analysis program.
#'
#' Create a COLONY infile using snpRdata sets containg offspring and possibly paternal genotypes given
#' specified parameters. Requires that the COLONY program is installed locally. Input and output files
#' will be stored in a colony folder created in the current working directory.
run_colony <- function(x, colony_path = ".",  outfile = "colony_input", method = "FPLS", run_length = 2,
                       sibship_prior = 0, paternal_sib_size = NULL, maternal_sib_size = NULL,
                       nruns = 1, seed = NULL, maternal_genotypes = NULL, paternal_genotypes = NULL,
                       maternal_inclusion_prob = 0, paternal_inclusion_prob = 0,
                       update_af = TRUE, dioecious = TRUE, inbreeding = TRUE, male_monogamous = FALSE,
                       female_monogamous = FALSE, clone_inference = FALSE, sibship_scaling = TRUE, known_af = FALSE,
                       precision = 2, dropout = 0.01, genotyping_error = 0.01,
                       known_maternal_dryads = NULL, known_paternal_dryads = NULL,
                       known_maternal_max_mismatches = 0, known_paternal_max_mismatches = 0,
                       known_maternal_sibships = NULL, known_paternal_sibships = NULL,
                       maternal_exclusions = NULL, paternal_exclusions = NULL,
                       excluded_maternal_siblings = NULL, excluded_paternal_siblings = NULL){

  # prep
  write_colony_input(x = x,
                     outfile = outfile,
                     method = method,
                     run_length = run_length,
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
                     known_maternal_dryads = known_maternal_dryads,
                     known_paternal_dryads = known_paternal_dryads,
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
                      path = "./colony/"))

}
