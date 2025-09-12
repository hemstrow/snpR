#' Basic forward genetic simulation with migration, mutation, and selection.
#' 
#' Simulates populations forward in time given potentially variable migration,
#' mutation, recombination, and selection processes. Requires a phased genotypic
#' input, although a snpRdata object can be phased and used automatically if
#' the program "beagle" is available.
#'
#' @param genotypes Either a snpRdata object, an list of unphased of genotype
#'   matrices, or NULL. If a matrices, one row per snp and one column per
#'   chromosome/gene copy, with allelic states designated with either 0 or 1. No
#'   missing data is allowed. Copies from individuals must be sequential. Each
#'   matrix is a population; a single matrix may be provided if only one
#'   population is desired. If a snpRdata object, requires \code{beagle_path}
#'   and \code{java_path} both be valid and point to the "beagle" software and
#'   "java", respectively.
#' @param meta data.frame of snp metadata or NULL. If a data.frame, first column
#'   must contain chromosome info, the second position in base pairs, and the
#'   third must be named "effects" and contain numeric effect sizes for each
#'   locus. If no selection is desired, each effect can be zero.
#' @param chr_length numeric vector of chromosome lengths \emph{in the same
#'   order as would be given by \code{unique(meta[,1])}}.
#' @param h numeric. Narrow-sense heritability, between 0 and 1.
#' @param BVs numeric vector, breeding values for starting individuals.
#'   Calculated from effects if not provided.
#' @param phenotypes numeric vector, phenotypes for starting individuals.
#'   Calculated from either provided BVs or effects if not provided.
#' @param gens numeric, maximum number of generations to simulate.
#' @param growth_function function, default \code{function(n) logistic_growth(n,
#'   500, 2)} the growth function to use during simulation, needs to take an
#'   argument \code{n}, the number of current individuals. See 
#'   \code{\link{population_growth_functions}}.
#' @param starting_surv_opt numeric, default NULL. Vector of initial optimum
#'   phenotypes for each population. If NULL, assumes that the mean BV in each
#'   population is the starting optimum.
#' @param mutation numeric, default 0. The per-base, per-generation mutation
#'   rate.
#' @param survival_function function, default \code{function(phenotypes,
#'   opt_pheno, ...) BL_survival(phenotypes, opt_pheno, omega = 1)}. Function
#'   for calculating survival probabilities of individuals. Needs two arguments:
#'   \code{phenotypes}: individual phenotypes and \code{opt_pheno}: optimum
#'   phenotype in a given generation. See \code{\link{survival_distributions}}.
#' @param selection_shift_function function, default \code{function(opt, iv)
#'   optimum_constant_slide(opt, iv, 0.3)}. Function for determining how much
#'   the fitness optimum changes each generation. Expects two arguments:
#'   \code{opt}: the current fitness optimum. \code{iv}: the initial genetic
#'   variance (the fitness optimum can slide as a function of the amount of
#'   initial genetic variance in the population). See
#'   \code{\link{optimum_slide_functions}}.
#' @param rec_dist function, default \code{function(n) rpois(n, lambda = 1)}.
#'   Function for determining how many recombination events occur per-generation,
#'   per-chromosome, per-cross. Requires the argument \code{n}, the number of
#'   chromosomes for which recombinations are needed.
#' @param thin logical, default FALSE. If TRUE, sites with no effect will be
#'   trimmed from analysis and reports unless all loci have no effects. This can
#'   massively reduce run time if few sites have effects. Set to FALSE if, for
#'   example, you wish to examine the effect of selection on linked, neutral
#'   loci or are running a simulation without selection.
#' @param thin_fixed logical, default TRUE. If TRUE, loci which become fixed in
#'   all population will be removed from the simulation to save on computation
#'   time.
#' @param stop_if_no_variance logical, default FALSE. If TRUE, will stop and
#'   return an error if no genetic variance remains in the population.
#' @param migration FALSE or a matrix of pairwise migration rates between
#'   populations. Each rate must be between 0 and 1, with rates given from row
#'   into column. Rows must sum to 1. Population order matches the genotypes
#'   provided or the order given by \code{\link{summarize_facets}} when
#'   requesting info on the provided \code{pop_facet}.
#' @param inbreeding numeric, default 0. A vector providing inbreeding targets
#'   in each population. Targets give the relative increase or decrease in the
#'   number of parental crosses each generation with a higher-than-average
#'   relatedness, so a value of 1 would double the number and a value of -1
#'   would half it. Population order matches the genotypes provided or the order
#'   given by \code{\link{summarize_facets}} when requesting info on the
#'   provided \code{pop_facet}. Only works when \code{do_sexes} is TRUE.
#' @param do_sexes logical, default TRUE. If TRUE, individuals are assigned
#'   randomly to one of two sexes when mating. The population will go extinct if
#'   all individuals are of a single sex in a generation.
#' @param mutation_effect_function function, default NULL. A function that 
#'   generates mutation effects given \code{n} mutations. If populations have
#'   different mutation effect generation, a list of a function for each can be
#'   provided instead.
#' @param pop_facet currently unused.
#' @param var_theta numeric, default 0. The variance in the optimum phenotype 
#'   each generation.
#' @param plot_during_progress logical, default TRUE.
#'
#' @export
simulate_populations <- 
  function(genotypes,
           meta,
           phenotypes = NULL,
           BVs = NULL,
           h,
           gens,
           chr_length,
           migration = FALSE,
           starting_surv_opt = NULL,
           mutation = 0,
           inbreeding = 0,
           ped = NULL,
           mutation_effect_function = NULL,
           growth_function = function(n) logistic_growth(n, 500, 2),
           survival_function = function(phenotypes, opt_pheno, ...) BL_survival(phenotypes, opt_pheno, omega = 1),
           selection_shift_function = function(opt, iv) optimum_constant_slide(opt, iv, 0.3),
           rec_dist = function(n) rpois(n, lambda = 1),
           var_theta = 0,
           plot_during_progress = FALSE,
           do_sexes = TRUE,
           effects = NULL,
           thin = FALSE,
           thin_fixed = TRUE,
           verbose = FALSE,
           print_all_freqs = F,
           print_all_thinned_freqs = FALSE,
           sampling_point = "migrants",
           stop_if_no_variance = FALSE,
           track_ped = FALSE,
           K_thin_post_surv = NULL,
           pop_facet = NULL){
    message("simulate_populations is still in development: please report any bugs!\n")
    
    if(plot_during_progress){
      .check.installed("ggh4x")
      .check.installed("scales")
    }
    if(inbreeding != 0){
      .check.installed("ggroups")
    }

    #============sub-fuctions============
    # function completed each loop (done this way to allow for ease of multiple populations)
    one_gen <- function(genotypes, phenotypes,
                        BVs,
                        effects,
                        opt,
                        survival_function,
                        K_thin_post_surv,
                        meta,
                        rec_dist,
                        chr_length,
                        do_sexes,
                        h.av,
                        model,
                        selection_shift_function,
                        mutation,
                        inbreeding,
                        pass_surv_genos = FALSE,
                        ped = NULL){
      
      
      #=========survival====
      if(length(unique(phenotypes)) == 1 & phenotypes[1] != 1 & stop_if_no_variance){stop("No genetic variance left.\n")}
      
      #survival:
      s <- rbinom(ncol(genotypes)/2, 1, #survive or not? Number of draws is the pop size in prev gen, survival probabilities are determined by the phenotypic variance and optimal phenotype in this gen.
                  survival_function(phenotypes, opt))
      
      #if the population has died out, stop.
      if(sum(s) <= 1){
        #adjust selection optima
        return(list(NA, NA))
      }
      
      # if doing carrying capacity on BREEDERS, not offspring, thin to K here.
      if(!is.null(K_thin_post_surv)){
        if(sum(s) > K_thin_post_surv){
          s[which(s == 1)][sample(sum(s), sum(s) - K_thin_post_surv, F)] <- 0
        }
      }
      
      #===============figure out next generation===============
      #what is the pop size after growth?
      offn <- round(growth_function(sum(s)))
      
      
      #make a new x with the survivors
      genotypes <- genotypes[, .SD, .SDcols = which(rep(s, each = 2) == 1)] #get the gene copies of survivors
      
      # save parents if needed
      if(pass_surv_genos){
        final_genotypes <- genotypes
        final_meta <- meta
        final_effects <- effects
      }
      
      # adjust ped
      if(any(s < 1) & !is.null(ped)){
        aped <- ped[ped$gen < max(ped$gen),]
        ped <- rbind(aped,
                     ped[ped$gen == max(ped$gen),][s == 1,])
      }
      
      #=============do random mating, adjust selection, get new phenotype scores, get ready for next gen====
      genotypes <- .rand.mating(x = genotypes, N.next = offn, meta = meta, rec_dist = rec_dist, chr_length, do_sexes,
                               mutation = mutation, inbreeding = inbreeding, ped = ped, verbose = verbose)
      if(is.null(genotypes)){return(list(NA, NA))} # return NA if empty (no individuals because everything was one sex)
      
      if(pass_surv_genos){
        return(list(genotypes = genotypes, final_genotypes = final_genotypes, final_meta = final_meta, final_effects = final_effects))
      }
      else{
        return(genotypes)
      }
    }
    
    do_mutation <- function(x.next, mutation, chr_length, meta, uf){
      #================mutation===================================
      # figure out number of mutations in each individual
      nmut <- lapply(x.next, function(z){
        if(!is.null(z)){return(rpois(ncol(z), mutation*sum(chr_length)))}
        else{return(NULL)}
      })
      pops <- unlist(lapply(nmut, length))
      nmut <- unlist(nmut)
      
      # exit if no mutations
      if(sum(nmut) == 0){
        meta$new <- FALSE
        return(list(genotypes = x.next, meta.next = meta, effects = effects))
      }
      
      # figure out positions for the mutations
      positions <- sample(sum(chr_length), sum(nmut), TRUE)
      chrs <- as.numeric(cut(positions, breaks = c(0, cumsum(chr_length))))
      positions <- positions - c(0, cumsum(chr_length)[-length(chr_length)])[chrs]
      inds <- rep(1:sum(pops), nmut) # which inds each mutation is in
      mut_info <- data.table::data.table(chrs = chrs, position = positions, ind = inds)
      # search for duplicates and re-sample them if they exist (unlikely)
      dups <- duplicated(mut_info[,1:3])
      loop_count <- 1
      while(any(dups)){
        if(loop_count > 20){
          stop("Tried to fix a double mutation at one locus in one individual 20 times without success--are your chromosomes too small for your mutation rate?\n")
        }
        npos <- sample(sum(chr_length), sum(dups), TRUE)
        nchrs <- as.numeric(cut(npos, breaks = c(0, cumsum(chr_length))))
        npos <- npos - c(0, cumsum(chr_length)[-length(chr_length)])[nchrs]
        chrs[dups] <- nchrs
        positions[dups] <- npos
        
        mut_info$chrs[dups] <- nchrs
        mut_info$position[dups] <- npos
        
        loop_count <- loop_count + 1
        dups <- duplicated(mut_info[,1:3])
      }
      
      # finish making the info df
      inds <- rep(1:sum(pops), nmut) # which inds each mutation is in
      mut_info <- dplyr::arrange(mut_info, inds, chrs, positions)
      mut_info_index <- unique(mut_info[,1:2])
      mut_info_index$row <- 1:nrow(mut_info_index)
      mut_info <- merge(mut_info, mut_info_index)
      
      # get effects and thin if doing so
      mut.eff <- matrix(NA, nrow(mut_info_index), npops)
      if(length(mutation_effect_function) > 1){
        for(i in 1:npops){
          mut.eff[,i] <- mutation_effect_function[[i]](nrow(mut_info_index))
        }
      }
      else{
        mut.eff[,1:npops] <- mutation_effect_function(nrow(mut_info_index))
      }
      
      
      # thin if requested
      if(thin){
        zeros <- which(rowSums(mut.eff) == 0)
        
        # cut and return if empty
        if(length(zeros) == nrow(mut_info_index)){
          meta$new <- FALSE
          return(list(x.next = x.next, meta.next = meta, effects = effects))
        }
        
        # otherwise adjust row info
        mut_info_index <- mut_info_index[-zeros,]
        mut_info <- mut_info[-which(mut_info$row %in% zeros),]
        mut_info_index$row <- 1:nrow(mut_info_index)
        mut_info$row <- NULL
        mut_info <- merge(mut_info, mut_info_index)
        mut.eff <- mut.eff[-zeros,,drop=FALSE]
      }
      
      # make the data
      mut.x <- matrix(0, nrow(mut_info_index), sum(pops))
      mut.x[mut_info$row + ((mut_info$ind - 1) * nrow(mut.x))] <- 1
      mut.x <- as.data.table(mut.x)
      mut_info$chrs <- uf[mut_info$chrs]
      
      # determine if any new sites overlap the existing ones
      overlap_i_in_ref <- which(paste0(meta[,1], "_", meta[,2]) %in% paste0(mut_info$chrs, "_", mut_info$position))
      if(length(overlap_i_in_ref) > 0){
        overlap_i_in_mut <- match(paste0(meta[overlap_i_in_ref,1], "_", meta[overlap_i_in_ref,2]),
                                  paste0(mut_info$chrs, "_", mut_info$position))
        
        # flip any overlaps, then remove these from the data
        for(i in 1:npops){
          this_pops <- which(mut_info[overlap_i_in_mut,]$ind > c(0, cumsum(pops))[i] &
                               mut_info[overlap_i_in_mut,]$ind <= c(0, cumsum(pops))[i + 1])
          # move to next pop if none in this one
          if(length(this_pops) == 0){next}
          
          # find the overlaps in this pop
          t.overlaps <- mut_info[overlap_i_in_mut,][this_pops]
          t.overlaps$ind <- t.overlaps$ind - c(0, cumsum(pops))[i]
          t.overlap.cols <- t.overlaps$ind
          
          # grab those rows, locate mutations, flip, and then update
          overlap_fill <- as.matrix(x.next[[i]][,..t.overlap.cols])
          overlap_indices <- t.overlaps$row + nrow(overlap_fill)*((1:nrow(t.overlaps)) - 1)
          overlap_fill[overlap_indices] <- ifelse(overlap_fill[overlap_indices] == 0, 1, 0)
          overlap_fill <- data.table::as.data.table(matrix(overlap_fill, nrow = nrow(x.next[[i]])))
          
          # fix the niche case where there are two mutation overlaps in one individual
          o_dups <- which(duplicated(t.overlap.cols) | duplicated(t.overlap.cols, fromLast = TRUE))
          if(length(o_dups) > 0){
            o_inds <- unique(t.overlap.cols[o_dups])
            rm_cols <- numeric(0)
            
            # for each dup, usually only one
            for(q in 1:length(o_inds)){
              toc <- o_inds[q]
              t_dup_fill <- as.matrix(x.next[[i]][,..toc])
              t_dup_fill[t.overlaps[ind == toc,]$row] <- ifelse(t_dup_fill[t.overlaps[ind == toc,]$row] == 0, 1, 0)
              data.table::set(overlap_fill, i = as.integer(1:nrow(overlap_fill)), j = which(t.overlap.cols == toc)[1], t_dup_fill)
              rm_cols <- c(rm_cols, which(t.overlap.cols == toc)[-1])
            }
            
            overlap_fill <- overlap_fill[,-..rm_cols]
            t.overlap.cols <- t.overlap.cols[-rm_cols]
          }
          
          # set
          data.table::set(x.next[[i]],
                          i = as.integer(1:nrow(x.next[[i]])),
                          j = as.integer(t.overlap.cols),
                          overlap_fill)
        }
        
        mut.eff <- mut.eff[-unique(mut_info$row[overlap_i_in_mut]),]
        mut_info_index <- mut_info_index[-which(mut_info_index$row %in% mut_info$row[overlap_i_in_mut]),]
        mut.x <- mut.x[-overlap_i_in_mut,, drop = FALSE]
      }
      
      # split into pops
      empties <- which(pops == 0)
      mut.x <- split(data.table::transpose(mut.x), rep(1:npops, pops))
      mut.x <- lapply(mut.x, data.table::transpose)
      if(length(empties) > 0){
        rep.mut.x <- vector("list", npops)
        rep.mut.x[-empties] <- mut.x
        mut.x <- rep.mut.x; rm(rep.mut.x)
      }
      
      for(i in 1:npops){
        if(!is.null(x.next[[i]])){
          x.next[[i]] <- rbind(x.next[[i]], mut.x[[i]])
        }
      }
      effects <- rbind(effects, mut.eff)
      mut_info_index$chrs <- uf[mut_info_index$chrs]
      mut_info_index$new <- TRUE
      meta$new <- FALSE
      mut_info_index$row <- NULL
      colnames(mut_info_index) <- colnames(meta)
      meta.next <- rbind(meta, mut_info_index)
      
      
      
      return(list(x = x.next, meta.next = meta.next, effects = effects))
    }
    
    update_phenotypes <- function(genotypes, effects, h, h.av){
      for(j in 1:length(genotypes)){
        pa <- .get.pheno.vals(genotypes[[j]], effect.sizes = effects[,j],
                             h = h,
                             hist.a.var = h.av[j],
                             phased = T)
        
        BVs[[j]] <- pa$a
        phenotypes[[j]] <- pa$p
      }
      
      return(list(BV = BVs, phenotypes = phenotypes))
    }
    
    #========================prep===========================
    if(verbose){cat("Initializing...\n")}
    if(!is.list(genotypes)){
      genotypes <- list(genotypes)
    }
    genotypes <- lapply(genotypes, data.table::as.data.table)
    if(length(unique(unlist(lapply(genotypes, nrow)))) != 1){
      stop("All genotype matrices must have the same number of loci.\n")
    }
    npops <- length(genotypes)
    
    if(is.data.table(meta)){meta <- as.data.frame(meta)}
    if(!is.data.frame(meta)){meta <- as.data.frame(meta)}
    if(colnames(meta)[2] != "position"){stop("The second column of meta must be `position`.\n")}
    
    # thin genotypes if possible
    zeros <- numeric()
    if("effect" %in% colnames(meta) & any(meta$effect != 0) & any(meta$effect == 0)){
      
      zeros <- which(meta$effect == 0)
      effects <- as.matrix(meta[,"effect", drop = FALSE])
    }
    else if(!is.null(effects)){
      if(!is.matrix(effects)){
        effects <- matrix(effects, ncol = 1)
      }
      if(nrow(effects) != nrow(genotypes[[1]])){
        stop("The number of locus effects is not equal to the number of loci.")
      }
      zeros <- which(rowSums(effects) == 0)
    }
    else{
      stop("Cannot locate effects in meta.\n")
    }
    
    if(thin & length(zeros) != 0 & length(zeros) != nrow(meta)){
      nmeta <- meta[-zeros,]
      # if chromosomes are missing, need to remove them from
      # the chr_length vector
      missing.chrs <- which(!unique(meta[,1]) %in% unique(nmeta[,1])) # which are now missing?
      if(length(missing.chrs) > 0){
        chr_length <- chr_length[-missing.chrs]
      }
      meta <- nmeta
      rm(nmeta)
      genotypes <- lapply(genotypes, function(x) x[-zeros,])
      effects <- effects[-zeros,,drop = FALSE]
    }
    
    
    
    if(is.null(phenotypes) | is.null(BVs)){
      need_phenos <- is.null(phenotypes)
      need_BVs <- is.null(BVs)
      
      for(i in 1:npops){
        
        p <- .get.pheno.vals(genotypes[[i]], effects[,i], h = h, phased = T)
        
        if(need_phenos){
          if(is.null(phenotypes)){
            phenotypes <- vector("list", npops)
          }
          phenotypes[[i]] <- p$p
        }
        if(need_BVs){
          if(is.null(BVs)){
            BVs <- vector("list", npops)
          }
          BVs[[i]] <- p$a
        }
        rm(p)
      }
    }
    
    colnames(effects) <- paste0("effects_", 1:npops)
    
    # sanity checks
    msg <- character()
    if(ncol(effects) != length(genotypes)){
      msg <- c(msg, "The number of effect columns must equal the number of populations.\n")
    }
    if(!is.null(starting_surv_opt)){
      if(length(starting_surv_opt) != length(genotypes)){
        msg <- c(msg, "The number of starting survival optima must equal the number of populations.\n")
      }
    }
    
    
    if(length(inbreeding) == 1){
      inbreeding <- rep(inbreeding, length(genotypes))
    }
    if(any(inbreeding != 0) | track_ped){
      if(!do_sexes){
        stop("Cannot do inbreeding without tracking sexes due to pedigree package limitations.\n")
      }
      if(is.null(ped)){
        ped <- data.frame(ID = 1:(sum(unlist(lapply(genotypes, ncol)))/2))
        ped$SIRE <- -1:(-1*nrow(ped))
        ped$DAM <- -1:(-1*nrow(ped))
        ped$gen <- 0
        ps <- unlist(lapply(genotypes, ncol))/2
        pop <- numeric(0)
        for(i in 1:length(ps)){
          pop <- c(pop, rep(i, ps[i]))
        }
        ped$pop <- pop
      }
      else{

        # check existing ped
        ## adjust gen
        ped$gen <- ped$gen - max(ped$gen)
        
        sizes <- unlist(lapply(genotypes, ncol))/2
        tab <- table(ped[ped$gen == max(ped$gen),]$pop)
        tab <- tab[sort(names(tab))]
        if(any(sizes != tab)){
          stop(paste0("Provided pedigree missing samples. Pedigree counts per pop:\n\t", 
               paste0(paste0(paste0("POP = ", names(tab)), "; n = ", tab), collapse = "\n\t"),
               "\nGenotype counts per pop:\n\t",
               paste0(paste0(paste0("POP = ", 1:length(genotypes)), "; n = ", sizes), collapse = "\n\t"), collapse = "\n"))
        }
      }
      
    }
    else{
      ped <- NULL
    }
    
    #================print out initial conditions, initialize final steps, and run===========
    #starting optimal phenotype, which is the starting mean additive genetic value.
    if(is.null(starting_surv_opt)){
      opt <- unlist(lapply(BVs, mean)) #optimum phenotype
    }
    else if(length(starting_surv_opt) == 1){
      opt <- rep(starting_surv_opt, length(genotypes))
    }
    else{
      opt <- starting_surv_opt
    }
    
    
    if(verbose){
      cat("\n\n===============done===============\n\nStarting parms:\n\tstarting optimum phenotype:", opt,
          "\n\tmean phenotypic value:", paste0(unlist(lapply(phenotypes, mean)), collapse = ", "),
          "\n\taddative genetic variance:", paste0(unlist(lapply(BVs, var)), collapse = ", "),
          "\n\tphenotypic variance:", paste0(unlist(lapply(phenotypes, var)), collapse = ", "), "\n\th:", h, "\n")
    }
    
    #make output matrix and get initial conditions
    out <- array(NA, c(gens + 1, 8, npops))
    colnames(out) <- c("N", "mu_phenotypes", "mu_BVs", "opt", "diff", "var_BVs", "stochastic_opt", "gen")
    N <- unlist(lapply(genotypes, ncol))/2 #initial pop size
    h.av <- unlist(lapply(BVs, var)) #get the historic addative genetic variance.
    h.pv <- unlist(lapply(phenotypes, var)) #historic phenotypic variance.
    
    out[1,,] <- t(as.matrix(data.frame(N,
                                       unlist(lapply(phenotypes, mean)),
                                       unlist(lapply(BVs, mean)),
                                       opt, 0, h.av, opt, 0))) #add this and the mean initial additive genetic variance
    if(plot_during_progress){
      pdat <- apply(out, 3, function(x) reshape2::melt(as.data.frame(x), id.vars = "gen"))
      pdat <- data.table::rbindlist(pdat, idcol = "pop")
      colnames(pdat) <- c("pop", "Generation", "var", "val")
      pdat$pop <- as.factor(pdat$pop)
      pdat <- na.omit(pdat)
      print(ggplot2::ggplot(pdat, ggplot2::aes(Generation, val)) + ggplot2::geom_point(na.rm = T) +
              ggh4x::facet_grid2(pop~var, scales = "free_y", independent = "y") +
              ggplot2::theme_bw() +
              ggplot2::scale_x_continuous(limits = c(1, gens), breaks = scales::pretty_breaks()) +
              ggplot2::theme(strip.placement = "outside", axis.title.y = ggplot2::element_blank(),
                             strip.background = ggplot2::element_blank(),
                             strip.text = ggplot2::element_text(size = 11)))
    }
    
    #initialize matrix to return allele frequencies if requested.
    if(print_all_freqs){
      num <- c(NA, 0)
      a.fqs <- array(num[1], c(nrow(meta), gens + 1, npops))
      a.fqs[,1,] <- unlist(lapply(genotypes, function(x) rowSums(x)/ncol(x)))
    }
    
    
    
    if(!is.list(survival_function)){
      survival_function <- list(survival_function)[rep(1, npops)]
    }
    if(!is.list(rec_dist)){
      rec_dist <- list(rec_dist)[rep(1, npops)]
    }
    if(!is.list(selection_shift_function)){
      selection_shift_function <- list(selection_shift_function)[rep(1, npops)]
    }
    
    if(length(K_thin_post_surv) == 1){
      K_thin_post_surv <- rep(K_thin_post_surv, npops)
    }
    
    if(length(chr_length) == 1){
      chr_length <- rep(chr_length, length(unique(meta[,1])))
    }
    
    if(length(mutation) == 1){
      mutation <- rep(mutation, length(genotypes))
    }
    
    if(is.list(mutation_effect_function)){
      if(length(mutation_effect_function) != npops){
        stop("Either a single mutation effect function or one for each population must be provided.\n")
      }
    }
    
    if(length(var_theta) == 1){
      var_theta <- rep(var_theta, npops)
    }
    
    
    
    if(verbose){
      cat("\nBeginning run...\n\n================================\n\n")
    }
    
    # init thinned tracker, which won't be updated if not thinning fixed loci
    track_thinned_afs <- thin_fixed & print_all_freqs & print_all_thinned_freqs
    if(track_thinned_afs){
      thinned_a.fqs <- array(NA, c(0, ncol(a.fqs), npops))
      thinned_effects <- matrix(NA, 0, ncol(effects))
    }
    else{
      thinned_a.fqs <- NULL
      thinned_effects <- NULL
    }
    
    #================loop through each additional gen, doing selection, survival, and fisher sampling of survivors====
    
    for(i in 2:(gens+1)){
      
      if(thin_fixed){
        fixed <- lapply(genotypes, function(z) rowSums(z) == 0)
        fixed <- matrix(unlist(fixed), ncol = length(genotypes))
        fixed <- which(rowSums(fixed) == ncol(fixed))
        
        if(length(fixed) > 0){
          genotypes <- lapply(genotypes, function(z) z[-fixed,])
          meta <- meta[-fixed,,drop=FALSE]
          
          # pull thinned loci out of a.fq
          if(print_all_freqs){
            
            if(track_thinned_afs){
              new_thinned_a.fqs <- array(NA, c(nrow(thinned_a.fqs) + length(fixed), ncol(thinned_a.fqs), npops))
              if(nrow(thinned_a.fqs) > 0){new_thinned_a.fqs[1:nrow(thinned_a.fqs),,] <- thinned_a.fqs}
              new_thinned_a.fqs[(nrow(thinned_a.fqs)+1):nrow(new_thinned_a.fqs),,] <- a.fqs[fixed,,]
              thinned_a.fqs <- new_thinned_a.fqs
              rm(new_thinned_a.fqs)
              thinned_effects <- rbind(thinned_effects, effects[fixed,,drop=FALSE])
            }
            a.fqs <- a.fqs[-fixed,,]
          }
          
          effects <- effects[-fixed,,drop=FALSE]
        }
      }
      
      gen_res <- vector("list", npops)
      
      # get the optimum phenotype(s) this gen
      t.opt <- opt
      for(j in 1:npops){
        t.opt[j] <- rnorm(1, opt[j], var_theta[j])
      }
      
      if(any(inbreeding != 0) | track_ped){
        nped <- ped
      }
      
      for(j in 1:npops){
        if(!is.null(genotypes[[j]])){

          if(is.null(ped)){tped <- NULL}
          else{
            tped <- ped[ped$gen < max(ped$gen) |
                          (ped$gen == max(ped$gen) & ped$pop == j),]
          }
          # grab the correct pedigree (only ancestors and invididuals in this population)
          genotypes[[j]] <- one_gen(genotypes = genotypes[[j]],
                                    phenotypes = phenotypes[[j]],
                                    BVs = BVs[[j]],
                                    effects = effects[,j],
                                    opt = t.opt[j],
                                    survival_function = survival_function[[j]],
                                    K_thin_post_surv = K_thin_post_surv[j],
                                    meta = meta,
                                    rec_dist = rec_dist[[j]],
                                    chr_length = chr_length,
                                    do_sexes = do_sexes,
                                    h.av = h.av[j],
                                    selection_shift_function = selection_shift_function[[j]],
                                    mutation = mutation[j],
                                    inbreeding = inbreeding[j],
                                    ped = tped,
                                    pass_surv_genos = ifelse(i == gens + 1 & sampling_point == "parents",
                                                             TRUE, FALSE)
          )
          
          # extract pedigree if needed and adjust
          if(!(i == gens + 1 & sampling_point == "parents")){
            tnped <- genotypes[[j]][[2]]
            genotypes[[j]] <- genotypes[[j]][[1]]
            
            if(!isFALSE(tnped)){
              tnped$ID <- tnped$ID + max(nped$ID)
              tnped$pop <- j
              nped <- rbind(nped, tnped)
            }
          }
          else{
            tnped <- genotypes[[j]][[1]][[2]]
            genotypes[[j]][[1]] <- genotypes[[j]][[1]][[1]]
            
            if(!isFALSE(tnped)){
              tnped$ID <- tnped$ID + max(nped$ID)
              tnped$pop <- j
              nped <- rbind(nped, tnped)
            }
          }
          
        }
      }
      
      # update ped
      if(any(inbreeding != 0) | track_ped){
        ped <- nped
      }
      
      if(i == gens + 1 & sampling_point == "parents"){
        final_genotypes <- purrr::map(genotypes, "final_genotypes")
        final_meta <- purrr::map(genotypes, "final_meta")
        final_effects <- purrr::map(genotypes, "final_effects")
        genotypes <- purrr::map(genotypes, "genotypes")
        
        pa <- update_phenotypes(genotypes = genotypes, effects = effects, h = h, h.av = h.av)
        
        final_BVs <- pa$BV; final_phenotypes <- pa$phenotypes
        
        final_genotypes <- lapply(final_genotypes, function(z){
          if(!is.data.table(z)){
            return(NULL)
          }
          else{
            return(z)
          }
        })
      }
      
      # replace NAs with NULLs, done this way to prevent list element removal
      genotypes <- lapply(genotypes, function(z){
        if(!is.data.table(z)){
          return(NULL)
        }
        else{
          return(z)
        }
      })
      
      if(all(unlist(lapply(genotypes, is.null)))){
        warning("All populations went extinct prior to designated number of generations.\n")
        res <- list(run_vars = out,
                    effects = effects, thinned_a.fqs = thinned_a.fqs, thinned_effects = thinned_effects,
                    meta = meta)
        if(print_all_freqs){
          res <- c(res, list(a.fqs = a.fqs))
        }
        
        return(res)
      }
      
      
      #====mutation======
      if(any(mutation > 0)){
        muts <- do_mutation(genotypes,
                            chr_length = chr_length, mutation = mutation,
                            meta = meta, uf = unique(meta[,1]))
        genotypes <- muts$x
        meta <- muts$meta.next
        meta$new <- NULL
        effects <- muts$effects
        
        # if tracking allele frequencies, add the new loci
        if(print_all_freqs & any(muts$meta.next$new)){
          new_a.fqs <- array(NA, c(nrow(meta), ncol(a.fqs), npops))
          new_a.fqs[1:nrow(a.fqs), 1:ncol(a.fqs),] <- a.fqs
          a.fqs <- new_a.fqs
          rm(new_a.fqs)
        }
        
        rm(muts); gc(FALSE)
      }
      
      #========migration================
      if(sampling_point == "offspring" & i == gens + 1){
        final_genotypes <- genotypes
        final_meta <- meta
        final_effects <- effects
        pa <- update_phenotypes(genotypes = genotypes, effects = effects, h = h, h.av = h.av)
        
        final_BVs <- pa$BV; final_phenotypes <- pa$phenotypes
      }
      
      if(!isFALSE(migration)){
        dest <- vector("list", npops)
        new_genotypes <- dest
        
        # if(any(unlist(lapply(genotypes, is.null)))){browser()}
        
        
        # assign destinations according to migration probabilities
        for(k in 1:npops){
          if(!is.null(genotypes[[k]])){
            dest[[k]] <- rmultinom(1, ncol(genotypes[[k]])/2, migration[k,])
            dest[[k]] <- rep(1:nrow(dest[[k]]), dest[[k]])
            dest[[k]] <- sample(dest[[k]], length(dest[[k]]), FALSE)
          }
        }
        
        # move
        nsizes <- table(unlist(dest))
        if(any(inbreeding != 0) | track_ped){
          nped <- ped[0,]
        }
        for(k in 1:npops){ # into k
          new_genotypes[[k]] <- data.table::as.data.table(matrix(0, ncol = nsizes[k]*2, nrow = max(unlist(lapply(genotypes, nrow)))))
          prog <- 0

          if(ncol(new_genotypes[[k]]) > 0){
            for(j in 1:npops){ # from j
              if(!is.null(genotypes[[j]])){
                movers <- which(rep(dest[[j]], each = 2) == k)
                
                if(length(movers) > 0){
                  data.table::set(new_genotypes[[k]], i = 1:nrow(new_genotypes[[k]]), j = (prog + 1):(prog + sum(dest[[j]] == k)*2),
                                  genotypes[[j]][,..movers])
                  
                  if(any(inbreeding != 0) | track_ped){
                    tnped <- ped[ped$gen == max(ped$gen) &
                                  ped$pop == j,][which(dest[[j]] == k),]
                    tnped$pop <- k
                    nped <- rbind(nped, tnped)
                  }
                  

                  prog <- (prog + (sum(dest[[j]] == k)*2))
                }
              }
            }
          }
          else{
            new_genotypes[[k]] <- 0
          }
        }
        if(any(inbreeding != 0) | track_ped){
          ped <- rbind(ped[ped$gen < max(ped$gen),],
                       nped)
        }
        
        genotypes <- new_genotypes
        rm(new_genotypes);gc(verbose = FALSE)
      }
      
      #========update phenotypes, etc==============
      # update genotypes, phenotypes, BVs, optima
      pa <- update_phenotypes(genotypes = genotypes, effects = effects, h = h, h.av = h.av)
      
      BVs <- pa$BV; phenotypes <- pa$phenotypes
      
      if(sampling_point == "migrants" & i != gens + 1){
        final_BVs <- BVs; final_phenotypes <- phenotypes
        final_meta <- meta
        final_genotypes <- genotypes
        final_effects <- effects
      }
      
      #========progress report==========
      out[i,1,] <- unlist(lapply(genotypes, function(x) sum(ncol(x)/2)))
      out[i,2,] <- unlist(lapply(phenotypes, mean))
      out[i,3,] <- unlist(lapply(BVs, mean))
      out[i,4,] <- opt
      out[i,5,] <- opt - out[i,3,]
      out[i,6,] <- unlist(lapply(BVs, var))
      out[i,7,] <- t.opt
      out[i,8,] <- i - 1
      
      if(verbose){
        cat("gen:", i - 1,
            "\tfixed_opt:", paste0(round(out[i,4,],3), collapse = ","),
            "\tstoch_opt", paste0(round(out[i,7,],3), collapse = ","),
            "\tmean(BVs):", paste0(round(out[i,3,],3), collapse = ","),
            "\tvar(BVs):", paste0(round(out[i,6,],3), collapse = ","),
            "\tlag:", paste0(round(out[i,4,],3) - round(out[i,3,],3), collapse = ","),
            "\tN:", paste0(out[i,1,], collapse = ","),"\n")
      }
      
      if(plot_during_progress){
        pdat <- apply(out, 3, function(x) reshape2::melt(as.data.frame(x), id.vars = "gen"))
        pdat <- data.table::rbindlist(pdat, idcol = "pop")
        colnames(pdat) <- c("pop", "Generation", "var", "val")
        pdat$pop <- as.factor(pdat$pop)
        pdat <- na.omit(pdat)
        print(ggplot2::ggplot(pdat, ggplot2::aes(Generation, val)) + ggplot2::geom_point(na.rm = T) +
                ggh4x::facet_grid2(pop~var, scales = "free_y", independent = "y") +
                ggplot2::theme_bw() +
                ggplot2::scale_x_continuous(limits = c(1, gens), breaks = scales::pretty_breaks()) +
                ggplot2::theme(strip.placement = "outside", axis.title.y = ggplot2::element_blank(),
                               strip.background = ggplot2::element_blank(),
                               strip.text = ggplot2::element_text(size = 11)))
      }
      
      if(print_all_freqs){
        a.fqs[,i,] <- unlist(lapply(genotypes, function(z) rowSums(z)/ncol(z)))
      }
      
      #==========update selection optima===========
      opt[j] <- selection_shift_function[[j]](opt[j], iv = sqrt(h.av[j]))
    }
    
    if(!print_all_freqs){
      num <- c(NA, 0)
      a.fqs <- array(num[1], c(max(unlist(lapply(genotypes, nrow))), 1, npops))
      for(i in 1:npops){
        if(!is.null(genotypes[[i]])){
          a.fqs[,1,i] <- rowSums(genotypes[[i]])/ncol(genotypes[[i]])
        }
      }
    }
    
    if(!exists("ped")){
      ped <- NULL
    }
    
    return(list(run_vars = out, genotypes = final_genotypes, phenotypes = final_phenotypes, BVs = final_BVs, a.fqs = a.fqs,
                effects = final_effects, thinned_a.fqs = thinned_a.fqs, thinned_effects = thinned_effects,
                meta = final_meta,
                ped = ped))
  }


.rand.mating <- function(x, N.next, meta, rec_dist, chr_length, do_sexes = TRUE, mutation = 0,
                        inbreeding = 0, ped = NULL, verbose = FALSE){
  
  #==========prep========
  if(length(unique(meta[,1])) != length(chr_length)){
    stop("The number of unique chromosomes is not equal to the number of chromsome lengths provided.\n")
  }
  if(!data.table::is.data.table(x)){
    x <- data.table::as.data.table(x)
  }
  facet <- colnames(meta)[1]

  if(inbreeding != 0){
    relatedness <- ggroups::buildA(ped)
    diag(relatedness) <- 0
    relatedness <- relatedness[which(ped$gen == max(ped$gen)),
                               which(ped$gen == max(ped$gen))]
  }
  
  #=========get parents and assign gcs for the next gen=======
  #make a new x with individuals in next gen
  ##find parents
  
  if(do_sexes){ # if there are two sexes
    sex <- rbinom(ncol(x)/2, 1, 0.5) #what are the parent sexes?
    if(sum(sex) == length(sex) | sum(sex) == 0){ # if every individual is the same sex, the population dies.
      return(NULL)
    }
    mates <- matrix(0, nrow = N.next, ncol = 2) # initialize, p1 and p2 are columns
    
    if(inbreeding != 0){
      # find valid pairings, select weighted by inbreeding depending on if they are less or greater than the average relatedness
      valid_sex_pairs <- outer(sex, sex, "+")
      valid_sex_pairs <- valid_sex_pairs == 1
      trelcut <- mean(relatedness[valid_sex_pairs])
      if(inbreeding > 0){
        sel_probs <- 1 + 1*(relatedness[which(valid_sex_pairs)] >= trelcut)*inbreeding # down or upranked by inbreeding request if greater than overall mean.
      }
      else{
        sel_probs <- 1 + 1*(relatedness[which(valid_sex_pairs)] <= trelcut)*-inbreeding # down or upranked by inbreeding request if greater than overall mean.
        
      }
      pairs <- sample(which(valid_sex_pairs), size = nrow(mates),  replace = TRUE, prob = sel_probs)
      
      if(verbose){
        cat("\t--mean relatedness:", trelcut, ";mean relatedness of pairs:", mean(relatedness[pairs]), "\n")
      }
      
      # fill
      mates[,1] <- row(relatedness)[pairs]
      mates[,2] <- col(relatedness)[pairs]
      
      # flip sexes
      flips <- which(sex[mates[,1]] == 0)
      mates[flips,] <- mates[flips,2:1]
      
    }
    else{
      mates[,1] <- which(sex == 1)[sample(sum(sex), nrow(mates), T)] #get parents of sex a
      mates[,2] <- which(sex == 0)[sample(length(sex) - sum(sex), nrow(mates), T)] #get parents of sex b
    }
    
  }
  else{ 
    mates <- matrix(sample(ncol(x)/2, N.next*2, T), ncol = 2) #p1 and p2 are columns
    selfings <- which(mates[,1] == mates[,2]) #any selfing?
    while(length(selfings) > 0){ #correct selfing
      mates[selfings,] <- sample(ncol(x)/2, length(selfings)*2, T) #get new parents
      selfings <- which(mates[,1] == mates[,2]) #any selfing remaining?
    }
  }
  
  # update ped
  if(!is.null(ped)){
    nped <- data.frame(ID = 1:nrow(mates),
                       SIRE = ped$ID[ped$gen == max(ped$gen)][mates[,1]],
                       DAM = ped$ID[ped$gen == max(ped$gen)][mates[,2]],
                       gen = max(ped$gen) + 1)
  }
  else{
    nped <- FALSE
  }
  
  
  # table with the showing the distribution of the number of offspring for each adult:
  # table(c(table(mates), rep(0, (ncol(x)/2) - length(table(mates)))))
  x.next <- data.table::as.data.table(matrix(0, 1, nrow(mates)*2))[rep(1, nrow(x))] #initialize x matrix for next gen
  
  #=========figure out which copy from each parent goes to offspring=====
  #randomly choose gene copies to push to individuals in the next gen.
  # for each individual, do they get copy 1 or copy 2 from the parent?
  uf <- unique(meta[,facet])
  chr.source.p1.i <- matrix(rbinom(nrow(mates)*length(uf), 1, .5), ncol = length(uf), byrow = T) + 1
  chr.source.p2.i <- matrix(rbinom(nrow(mates)*length(uf), 1, .5), ncol = length(uf), byrow = T) + 1
  
  # add the correct copies
  ##which column does the data come from?
  chr.source.p1 <- ifelse(chr.source.p1.i == 1, mates[,1] * 2 - 1,
                          mates[,1] * 2)
  
  chr.source.p2 <- ifelse(chr.source.p2.i == 1, mates[,2] * 2 - 1,
                          mates[,2] * 2)
  
  #the other copies?
  chr.nsource.p1 <- ifelse(chr.source.p1.i == 1, mates[,1] * 2,
                           mates[,1] * 2 - 1)
  chr.nsource.p2 <- ifelse(chr.source.p2.i == 1, mates[,2] * 2,
                           mates[,2] * 2 - 1)
  
  rm(chr.source.p1.i, chr.source.p2.i)
  
  #these now say which column in x to take the data from for each chr for each individual for bases that didn't recombine.
  
  #=========recombination and chromosome assignment======
  num.rec <- rec_dist(nrow(mates)*2*length(uf))
  rs <- sum(num.rec) #total number
  n.rec.per.chr <- tapply(num.rec, rep(chr_length, each = nrow(mates)*2), sum) # figure out how many recombination events per chromosome
  rec.pos <- runif(rs, min = 0, max = rep(chr_length, n.rec.per.chr)) # get positions for each of these events.
  
  prog <- 0 #progress through num.rec tracker.
  pos.prog <- 0 #progress through recombination events tracker.
  
  #fill for each facet
  for(j in 1:length(unique(uf))){
    # cat(j,"\n")
    #overall approach: for each parent:
    # line up copy 1 and copy 2 chromosomes in a matrix, each column is a seperate chr. Copy one is the copy that is getting passed! Copy two is the one that isn't.
    # for each recombination event, on each chr, flip which chr we are taking from. For portions where we are taking from chr2, paste into chr 1, which will be the output.
    this.chr.pos <- meta$position[meta[,facet] == uf[j]]
    
    for(k in 1:2){
      
      trec <- c(prog + 1, prog + nrow(mates)) #which recombination events are we working with here?
      prog <- prog + nrow(mates)
      
      #get the number of recombination events per chr in this set
      tnr <- num.rec[trec[1]:trec[2]]
      
      #initialize matrix
      c1.mat <- data.table::as.data.table(matrix(0, length(this.chr.pos), ncol = nrow(mates)))
      c2.mat <- data.table::as.data.table(matrix(0, length(this.chr.pos), ncol = nrow(mates)))
      
      #paste in the values from x. c1.mat contains the "passed" chr, c2.mat contains the "unpassed" chr.
      if(k == 1){
        data.table::set(c1.mat, 1:nrow(c1.mat), as.integer(1:ncol(c1.mat)),
                        x[i = which(meta[,facet] == uf[j]), .SD, .SDcols = chr.source.p1[,j]])
        data.table::set(c2.mat, 1:nrow(c2.mat), as.integer(1:ncol(c2.mat)),
                        x[i = which(meta[,facet] == uf[j]), .SD, .SDcols = chr.nsource.p1[,j]])
      }
      else{
        data.table::set(c1.mat, 1:nrow(c1.mat), as.integer(1:ncol(c1.mat)),
                        x[i = which(meta[,facet] == uf[j]), .SD, .SDcols = chr.source.p2[,j]])
        data.table::set(c2.mat, 1:nrow(c2.mat), as.integer(1:ncol(c2.mat)),
                        x[i = which(meta[,facet] == uf[j]), .SD, .SDcols = chr.nsource.p2[,j]])
      }
      
      #only the recombining entries.
      wnz <- which(tnr != 0)
      if(length(wnz) == 0){
        #no recombination, mostly if pop size is VERY small...
        if(k == 1){
          data.table::set(x.next, which(meta[,facet] == uf[j]), as.integer(seq(1, ncol(x.next), by = 2)), c1.mat)
        }
        else{
          data.table::set(x.next, which(meta[,facet] == uf[j]), as.integer(seq(2, ncol(x.next), by = 2)), c1.mat)
        }
        next()
      }
      tnr_nz <- tnr[wnz]
      
      
      
      #get the positions of the recombination events
      trpos <- rec.pos[(pos.prog + 1):(pos.prog + sum(tnr_nz))]
      pos.prog <- pos.prog + sum(tnr_nz)
      
      
      # Now need to make and assign the actual output vector. I can't think of a good way to do this in an actually vectorized way, mostly because I have to look up positions to reference against for each chr.
      sort.prog.trpos <- 0 #how many positions have we searched through?
      for(m in 1:length(tnr_nz)){
        # cat("\t\t", m, "\n")
        sort.pos <-
          c(sort(trpos[(sort.prog.trpos + 1):(sort.prog.trpos + tnr_nz[m])])) #sort the correct recomb positions and add the ending positions.
        sort.prog.trpos <- sort.prog.trpos + tnr_nz[m] #update.
        
        #now figure out which chr each position will be drawn from.
        chr.ident <- numeric(length(this.chr.pos))
        for(q in 1:length(sort.pos)){
          chr.ident <- chr.ident + (this.chr.pos <= sort.pos[q])
        }
        chr.ident <- chr.ident %% 2 #if the number of crossing over events was even, 0s mean copy 1 and 1s mean copy 2. Otherwise reversed.
        
        #assign.
        if(length(sort.pos) %% 2 == 0){
          data.table::set(c1.mat, which(chr.ident == 1), j = wnz[m], value = c2.mat[which(chr.ident == 1), wnz[m], with = FALSE])
          #assign entries where the chr flips in c1 mat to the respecitve entries in c2 mat.
        }
        else{
          data.table::set(c1.mat, which(chr.ident == 0), j = wnz[m], value = c2.mat[which(chr.ident == 0), wnz[m], with = FALSE])
        }
        
      }
      
      #================assign chromosomes to x.next.==========
      if(k == 1){
        data.table::set(x.next, which(meta[,facet] == uf[j]), as.integer(seq(1, ncol(x.next), by = 2)), c1.mat)
      }
      else{
        data.table::set(x.next, which(meta[,facet] == uf[j]), as.integer(seq(2, ncol(x.next), by = 2)), c1.mat)
      }
    }
  }
  
  return(list(x.next, nped))
}

#function to generate random environmental effects for a given set of BVs, addative genetic variance (either historic or current from BVs are typical), and heritability.
# note that standardization isn't perfect, but will result in data with a mean very close to 0 and var close to 1. The randomness of assigning random environmental effects everywhere will make it imperfect
# downstream standardization of the phenotypes will fix this if using estimated effect sizes!
.e.dist.func <- function(A1, hist.a.var, h, standardize = F){
  esd <- sqrt((hist.a.var/h)-hist.a.var) # re-arrangement of var(pheno) = var(G) + var(E) and h2 = var(G)/var(pheno)
  env.vals <- rnorm(length(A1), 0, esd)
  
  #if it standardization is requested, do so
  if(standardize){
    env.vals <- env.vals/sqrt(var(env.vals)/(1-h)) # set variance to 1 - h
    env.vals <- env.vals - mean(env.vals) # set mean to 0. Should be close, but not perfect because of the random draws.
  }
  return(env.vals)
}

#get phenotypic values given genotypes, effect sizes, and heritabilities. If hist.a.var is true, uses the amount of genomic variability this gen and h to figure out how big of an env effect to add. Otherwise uses the provided value (probably that in the first generation).
.get.pheno.vals <- function(x, effect.sizes, h, hist.a.var = "fgen", standardize = FALSE, phased = FALSE){
  # additive
  if(is.null(ncol(effect.sizes))){
    
    # if not loci/effects
    if(length(effect.sizes) == 0){
      return(list(p = rep(0, ncol(x)/2),
                  a = rep(0, ncol(x)/2)))
    }
    
    # remove zeros
    zeros <- which(effect.sizes == 0)
    if(length(zeros) != 0 & length(zeros) != length(effect.sizes)){
      effect.sizes <- effect.sizes[-zeros]
      x <- x[-zeros,, drop = FALSE]
    }
    
    #get effect of each individual:
    a.ind <- .weighted.colSums(x, effect.sizes) # faster than t(x)%*%effect.sizes!
    if(phased){
      a.ind <- a.ind[seq(1, length(a.ind), by = 2)] + a.ind[seq(2, length(a.ind), by = 2)] #add across both gene copies.
    }
    if(is.matrix(a.ind)){
      a.ind <- as.numeric(a.ind)
    }
  }
  
  
  # make sure h isn't 0, less than 0, or greater than 1
  if(h <= 0){
    h <- 0.00000001
  }
  else if(h > 1){
    h <- 1
  }
  
  #standardize the genetic variance if requested.
  if(standardize){
    a.ind <- a.ind/sqrt(var(a.ind)/h) # set the variance to h.
    a.ind <- a.ind - mean(a.ind) # set the mean to 0
  }
  
  #add environmental variance
  if(hist.a.var == "fgen"){
    pheno <- a.ind + .e.dist.func(a.ind, var(a.ind), h, standardize)
  }
  else{
    pheno <- a.ind + .e.dist.func(a.ind, hist.a.var, h, standardize)
  }
  
  return(list(p = pheno, a = a.ind))
}



#' Calculate survivability probabilities.
#'
#' For a given optimum phenotype and surivival variance, calculate survival
#' probabilities given a specific distribution. For use in
#' \code{\link{simulate_populations}}.
#'
#' Sevaral survival probablity distribution curves are provided. Each are coded
#' such that they expect a vector of phenotypes as the first argument and the
#' optimum phenotype as the second, with other arguments following, ending in
#' the `...` argument so that additional, unused arguments can be passed without
#' error. Other distributions can be coded in the same way if needed; call
#' \code{normal_survival} without parentheses to see an example.
#'
#' Options:
#'  * normal_survival: Probabilities are given by
#' \code{\link[stats]{pnorm}}. For values above the
#' optimum phenotype, the surivial is given by 1 - the resulting value. If a
#' maxium survival other than 0.5 is given, the surivivals are scaled such that
#' individuals with the optimum phenotype will have the provided survival
#' probability.
#'  * BL_survival: Survival according to Burger and Lynch's (1995) 
#' equation 1.

#' @param phenotypes numeric. Vector of phenotypes.
#' @param opt_pheno numeric. Optimum phenotype.
#' @param var numeric. Variance of the surivial distribution.
#' @param max.survival numeric. Maximum possible survival probability.
#' @param omega numeric. Strength of stabilizing selection. Smaller numbers mean
#'   stronger selection.
#' @param ... args passed internally. Self-made slide functions should have this
#'   arg to allow for unused args to be passed without error.
#' 
#' @name survival_distributions
#' @aliases normal_survival BL_survival
#'
#' @return A vector of survival probabilities for each phenotype.
NULL

#'@describeIn survival_distributions normally distributed survival probabilities
#'@export
normal_survival <- function(phenotypes, opt_pheno, var, max.survival = .5, ...){
  if(max.survival > 1){
    warning("Max surivial is above 1, will reset to 1.")
    max.survival <- 1
  }
  # get survival probabilities from normal curve
  x <- pnorm(phenotypes, opt_pheno, sqrt(var))
  x <- ifelse(x > .5, 1 - x, x)
  
  # scale for target max survival
  scale_factor <- max.survival/.5
  x <- x*scale_factor
  
  return(x)
}

#'@describeIn survival_distributions Burger and Lynch (1995) distributed
#'  survival probabilities
#'@export
BL_survival <- function(phenotypes, opt_pheno, omega, ...){
  exp(-((phenotypes - opt_pheno)^2)/(2*omega^2))
}


#' Increase the selection optimum each generation.
#' 
#' Increase selection optimum by either a given percentage of 
#' starting variance each generation or by a fixed amount. For use in
#' \code{\link{simulate_populations}}.
#' 
#' While two slide options are provided, others can be used instead by writing
#' functions which take the argument \code{opt}, the starting optimum in the
#' generation previous. Optionally, they
#' may also take the argument \code{iv}, which is the initial phenotypic
#' variation in the population, passed automatically by
#' \code{\link{simulate_populations}}. All functions must also take `...` to
#' allow for unused extra arguments pased by \code{\link{simulate_populations}},
#' such as \code{iv} for slides which do not use \code{iv}. Call
#' \code{optimum_scaled_slide} without parentheses to see an example.
#'
#' Options:
#'  * optimum_scaled_slide: Change the optimum phenotype 
#' each generation by an amount that depends on the initial 
#' amount of phenotypic variation in the population
#'  * optimum_fixed_slide: Change the optimum phenotype 
#' by a fixed amount each generation.
#'
#' @param opt numeric. Initial optimum phenotype.
#' @param iv numeric. Initial variance. The optimum will increase by some
#'   portion of this variance.
#' @param slide numeric. Proportion of the initial variance by which the optimum
#'   phenotype will slide.
#' @param ... args passed internally. Self-made slide functions should have this
#'   arg to allow for unused args to be passed without error.
#'
#' @name optimum_slide_functions
#' @aliases optimum_scaled_slide optimum_fixed_slide
#' @return The optimum phenotypic value in the next generation.
NULL

#'@describeIn optimum_slide_functions slide scaled by genetic variance
#'@export
optimum_scaled_slide <- function(opt, iv, slide = 0.3, ...){
  return(opt + iv*slide)
}

#'@describeIn optimum_slide_functions fixed slide each generation
#'@export
optimum_fixed_slide <- function(opt, slide = 0.3, ...){
  return(opt + slide)
}


#' Calculate population size changes each generation.
#' 
#' Increase population sizes each generation according to specific growth
#' rate equations. For use in \code{\link{simulate_populations}}.
#' 
#' Sevaral survival probablity distribution curves are provided. Each are coded
#' such that they expect \code{n}, the starting population size and end in
#' the `...` argument so that additional, unused arguments can be passed without
#' error. Other distributions can be coded in the same way if needed; call
#' \code{logistic_growth} without parentheses to see an example.
#'
#' Options:
#' * logistic_growth: Classical logistic growth model.
#' * BL_growth: Growth based on Burger and Lynch (1995). Assumes
#' each individual has \code{B} offspring (so each pair has 2*B). 
#' Density dependence must be supplied elsewhere, such as to the 
#' \code{K_thin_post_surv} argument to 
#' \code{\link{simulate_populations}}.
#'
#' @param n numeric, starting population size.
#' @param K numeric, carrying capacity.
#' @param r numeric, logistic growth rate.
#' @param B numeric, the number of offspring each individual has
#' @param ... args passed internally. Self-made slide functions should have this
#'   arg to allow for unused args to be passed without error.
#'
#' @name population_growth_functions
#' @aliases logistic_growth BL_growth
#' @return The number of individuals in the next generation.
NULL

#'@describeIn population_growth_functions logisitic growth
#'@export
logistic_growth <- function(n, K, r, ...){
  return((K*n*exp(r))/(K + n*(exp(r*1) - 1)))
}

#'@describeIn population_growth_functions Burger and Lynch (1995) growth
#' @export
BL_growth <- function(n, B, ...){
  return(B*n)
}


.basic_convert_2_to_1_column <- function(x){
  if(!is.matrix(x)){x <- as.matrix(x)}
  ind.genos <- x[,seq(1,ncol(x), by = 2)] + x[,seq(2,ncol(x), by = 2)]
  ind.genos <- matrix(ind.genos, nrow = ncol(x)/2, byrow = T) # rematrix and transpose!
  return(ind.genos)
}

.weighted.colSums <- function(data, weights) crossprod(as.matrix(data), weights)



