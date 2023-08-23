run_easypop <- function(nloci = 100, two_sexes = TRUE, mating_system = "random",
                        starting_population_sizes = 100,
                        proportion_female = 0.5,
                        migration_model = 3,
                        migration_model_switch_gen = generations/2,
                        recombination = "free",
                        mutation_rate = 0.001,
                        mutation_model = "Kam",
                        mKam_proportions = 0.1,
                        mSsm_proportion = 0.1,
                        migration_proportions = 0.05,
                        generations = 1000,
                        replicates = 1,
                        prompt_file = NULL,
                        initial_variability = "fixed",
                        save_pedigress = FALSE,
                        outfiles = "easypop_from_snpR",
                        easypop_path = "/usr/bin/easypop.revised.windows.exe"){
  #browser()
  
  # write the prompt file--running interactively to generate this will error via system on windows...
  
  #=============mutation schemes==============
  mm <- match(mutation_model, c("Kam", "Ssm", "mKam", "mSsm"))
  if(length(mm) == length(mutation_rate)){
    mr <- mutation_rate
  }
  else{
    if(length(mm) == 1){
      mm <- rep(mm, length(mutation_rate))
      mr <- mutation_rate
    }
    else if(length(mutation_rate) == 1){
      mr <- rep(mutation_rate, length(mm))
    }
    else{
      stop("Equal numbers of mutation_rates and mutation_models must be supplied.\n")
    }
  }
  args <- c(mutation_rate = mr,
            mutation_model = mm,
            all_loci_same_mutation_scheme = ifelse(length(unique(mutation_rate)) == 1 & length(unique(mutation_model)) == 1, "y", "n"))
  
  if(any(mm == 3)){
    if(length(mKam_proportions) == 1){
      mKam_proportions <- rep(mKam_proportions, sum(mm == 3))
    }
    else if(length(mKam_proportions) != sum(mm == 3)){
      stop("One mKam_porportion must be supplied for each mKam loci.\n")
    }
    args <- c(args, locus_numbers_and_proportions_kam_mutation_events = paste0(paste0(which(mm == 3), ",", mKam_proportions),
                                                                               collapse = ";"))
  }
  
  if(any(mm = 4)){
    if(length(mSsm_proportions) == 1){
      mSsm_proportions <- rep(mSsm_proportions, sum(mm == 4))
    }
    else if(length(mSsm_proportions) != sum(mm == 4)){
      stop("One mSsm_porportion must be supplied for each mSsm loci.\n")
    }
    args <- c(args, locus_numbers_and_proportions_kam_mutation_events = paste0(paste0(which(mm == 4), ",", mSsm_proportions),
                                                                               collapse = ";"))
  }
  
  if(length(mm) > 1 | length(mr) > 1){
    args <- c(args, per_locus_number_possible_allelic_states = rep(2, length(mr)))
  }
  else{
    args <- c(args, number_possible_allelic_states = 2)
  }
  
  #============migration schemes=============
  if(length(starting_population_sizes) == 1){
    args <- c(args, migration_model = 3)
  }
  else{
    if(length(migration_model) == 1){
      args <- c(args, migration_model = which(c("1dss", "2dss", "island", "hss", "hi", "spatial") == mating_system))
    }
    if(length(migration_model) == 2){
      args <- c(args, migration_model = which(c("1dss", "2dss", "island", "hss", "hi", "spatial") == mating_system[1]))
      args <- c(args, migration_model_second_scheme = which(c("1dss", "2dss", "island", "hss", "hi", "spatial") == mating_system[2]))
      
      args <- c(args, number_of_generations_before_migration_scheme_change = migration_model_switch_gen)
    }
  }
  
  
  #============recombination schemes=========
  args <- c(args, free_recombination_between_loci = ifelse(recombination == "free", "y", "n"))
  if(args["free_recombination_between_loci"] == "n"){
    args <- c(args, recombination_rate_between_adjacent_loci = recombination)
  }
  
  
  #============remaining args================
  args <- c(ploidy = 2,
            two_sexes = two_sexes,
            mating_sytem = which(c("random", "polygyny", "monogyny") == mating_system),
            number_populations = length(starting_population_sizes),
            same_number_individuals_each_population = ifelse(length(unique(starting_population_sizes)) == 1, "y", "n"),
            number_females_each_population = proportion_female*starting_population_sizes,
            number_males_each_population = (1 - proportion_female)*starting_population_sizes,
            number_of_loci = nloci,
            variability_initial_population = which(c("random", "fixed") == mating_system),
            number_of_generations = generations,
            complete_dataset_in_dat_and_gen_files = "y",
            file_giving_pedigrees = ifelse(save_pedigress, "y", "n"),
            name_of_file = outfiles,
            number_of_replicates = replicates,
            migration_model = 3,
  )
  if(args["same_number_individuals_each_population"] == "n"){
    args <- c(args, per_population_individual_counts = starting_population_sizes)
  }
  
  for(i in 1:length(args)){
    write(paste0(names(args)[i], ":	", paste0(args[i], collapse = ",")))
  }
  
  
}