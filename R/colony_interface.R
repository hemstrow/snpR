write_colony_input <- function(x, method = "FPLS", run_length = 2,
                               sibship_prior = 0, paternal_sib_size = NULL, maternal_sib_size = NULL,
                               nruns = 1, seed = NULL,
                               update_af = TRUE, dioecious = TRUE, inbreeding = FALSE, male_monogamous = FALSE,
                               female_monogamous = FALSE, clone_inference = FALSE, sibship_scaling = TRUE, known_af = FALSE,
                               precision = 2, dropout = 0.01, genotyping_error = 0.01){


  #=====================initialize and write first stuff=============
  browser()
  # initialize storage directory
  if(!dir.exists("colony")){
    dir.create("colony")
    setwd("colony")
  }
  else{
    setwd("colony")
  }

  # write anything preliminary to output
  proj <- deparse(substitute(x))
  if(nchar(proj) > 40){
    proj <- substr(proj, 1, 40)
  }
  outfile <- paste0(paste0("colony_input_", proj, collapse = "_"), ".txt")

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


  #=====================reformat genotypes==================


  # allele frequencies
  write(as.numeric(known_af), outfile, append = T)
  if(known_af){

  }

}
