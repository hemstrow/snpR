# This contains functions that are currently depreciated. Some may be restored in the future.










#calculates a pop weighted stat. Needs columns named "n_total" and one with a name matching
#the "stat" argument, which should be a character string. "boots" argument is the number of bootstraps
#to do do get standard errors

#'Weighted stat averages/SEs.
#'
#'\code{calc_weighted_stat} calculates a weighted average statistic and standard errors for a variable via the bootstrapping method described by Gatz and Smith (1995), The standard error of a weighted mean concentration—I. Bootstrapping vs other methods.
#'
#'Description of x:
#'    Must contain colums containing the statistic of interest (such as pi from calc_pi), population ID, and the weights (such as the number of alleles sequenced, or the n_total column from in the allele count format) for each SNP. These columns must be named to match the stat argument, "pop", and either "n_total" or "nk", respectively.
#'
#' @param x Input statistic data.
#' @param stat Name of the statistic for which to calculate a weighted mean and se.
#' @param boots Number of bootstraps to perform to calcluate se. 30 is usually sufficient.
#'
#' @return A data frame with columns for the population ID, weighted mean, and weighted average.
#'
#' @examples
#' calc_weighted_stat(randPI, "pi")
#'
calc_weighted_stat <- function(x, stat, boots = 30){
  pops <- unique(x$pop)
  out <- matrix(0, length(pops), 3)
  out[,1] <- pops
  get_wm <- function(y){
    if(any(colnames(y) == "n_total")){
      w <- y$n_total
    }
    else if(any(colnames(y) == "nk")){
      w <- y$nk
    }
    else{
      stop("No weighting column found. Weighting column must be named either nk or n_total.")
    }
    return(sum(w*y[,stat])/sum(w))
  }
  cat("working on pop:\n")
  for(i in 1:length(pops)){
    #mean
    cat("\n\t", pops[i],"\n")
    y <- x[x$pop == pops[i],]
    y <- y[!is.na(y[,stat]),]
    out[i,2] <- get_wm(y)


    #se, via bootstrap.
    n <- nrow(y)
    bvals <- numeric(boots)
    for(j in 1:boots){
      if(j %% 10 == 0){cat("\t\tboot:", j, "\n")}
      b <- y[sample(n, n, T),]
      bvals[j] <- get_wm(b)
    }
    out[i,3] <- sd(bvals)
  }
  out <- as.data.frame(out)
  colnames(out) <- c("pop", paste0(stat, "_wm"), paste0(stat,"_se"))
  out[,2] <- as.numeric(as.character(out[,2]))
  out[,3] <- as.numeric(as.character(out[,3]))
  return(as.data.frame(out))
}


#estimate S in one generation

#'Selection coefficients.
#'
#'Estimates the selection coefficient \emph{S} across one generation given intial and final allele frequencies, using an equation adapted from Gillespie (2010) Population genetics: a concise guide, equation 3.2.
#'
#' @param p_before Intial frequency of allele \emph{p}.
#' @param p_after Frequency of allele \emph{p} after one generation.
#'
#' @return A numeric value, the estimate of \emph{s}.
#'
#' @examples
#' estimate_selection(.01, .011)
#'
estimate_selection <- function(p_before, p_after){
  top <- 1-(p_before/p_after)
  q <- 1 - p_before
  bottom <- q^2
  s <- top/bottom
  return(s)
}


#Function to test for deviation from HWE via perumtaiton (given small cell entries) for each snp.
#inputs: x: data, snps in subsequent columns
#        num_start: column index with the first column of data
#        miss: CHARACTER entry which encodes missing data
#        test: the test to be used.
#            options: exact: Exact test, from Wigginton et al 2005.
#                     permutation: permutation chi.square test.
#        n.reps: number of permutations to get p.values. Needed if test = "permutation".

#'Hardy-Weinburg Equilibrium from SNP data.
#'
#'Function to test for deviation from HWE for each SNP via either permutation or the exact test given by Wigginton et al 2005.
#'
#'Part of this code is edited from Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000. code available at http://csg.sph.umich.edu/abecasis/Exact/snp_hwe.r
#'
#'Description of x:
#'    Contains metadata in columns 1:ecs. Remainder of columns contain genotype calls for each individual. Each row is a different SNP, as given by format_snps output options 4 or 6. Note that this *ignores populations*, split data into distinct populations before running.
#'
#' @param x Input data, in either NN or 0000 format, as given by format_snps output options 4 or 6.
#' @param ecs The number of extra columns at the start of the input data.frame containing metadata.
#' @param miss Characters which code for missing *genotypes* ("NN" or "0000", for example).
#' @param test Which test to use? "exact" for the exact test, "permutation" for permutation test.
#' @param n.reps Number of reps to use for permutation test.
#'
#' @return A numeric vector containing p-values, which are the probability that a SNP in HWE would be as or more diverged from the expected equilibrium genotype frequencies. In the same order as input.
#'
#' @examples
#' HWE(stickSNPs, 3)
HWE <- function(x, ecs, mDat = "NN", test = "exact", n.reps = 20000){


  if(test == "exact"){
    #edited from Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
    #Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000
    #code available at http://csg.sph.umich.edu/abecasis/Exact/snp_hwe.r
    SNPHWE <- function(x){
      obs_homr <- x[1]
      obs_homc <- x[2]
      obs_hets <- x[3]
      if (obs_homr < 0 || obs_homc < 0 || obs_hets < 0){
        return(-1.0)
      }
      # total number of genotypes
      N <- obs_homr + obs_homc + obs_hets

      # number of rare allele copies
      rare  <- obs_homr * 2 + obs_hets

      # Initialize probability array
      probs <- rep(0, 1 + rare)

      # Find midpoint of the distribution
      mid <- floor(rare * ( 2 * N - rare) / (2 * N))
      if ( (mid %% 2) != (rare %% 2) ) mid <- mid + 1

      probs[mid + 1] <- 1.0
      mysum <- 1.0

      # Calculate probablities from midpoint down
      curr_hets <- mid
      curr_homr <- (rare - mid) / 2
      curr_homc <- N - curr_hets - curr_homr

      while ( curr_hets >=  2)
      {
        #equation 2
        probs[curr_hets - 1]  <- probs[curr_hets + 1] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0)  * (curr_homc + 1.0))
        mysum <- mysum + probs[curr_hets - 1]

        # 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
        curr_hets <- curr_hets - 2
        curr_homr <- curr_homr + 1
        curr_homc <- curr_homc + 1
      }

      # Calculate probabilities from midpoint up
      curr_hets <- mid
      curr_homr <- (rare - mid) / 2
      curr_homc <- N - curr_hets - curr_homr

      while ( curr_hets <= rare - 2)
      {
        probs[curr_hets + 3] <- probs[curr_hets + 1] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
        mysum <- mysum + probs[curr_hets + 3]

        # add 2 heterozygotes -> subtract 1 rare homozygtote, 1 common homozygote
        curr_hets <- curr_hets + 2
        curr_homr <- curr_homr - 1
        curr_homc <- curr_homc - 1
      }

      # P-value calculation
      target <- probs[obs_hets + 1]

      #plo <- min(1.0, sum(probs[1:obs_hets + 1]) / mysum)

      #phi <- min(1.0, sum(probs[obs_hets + 1: rare + 1]) / mysum)

      # This assignment is the last statement in the fuction to ensure
      # that it is used as the return value
      p <- min(1.0, sum(probs[probs <= target])/ mysum)
    }
  }

  #figure out snp format
  if(nchar(as.character(x[1,(ecs + 1)])) == 2){
    snp_form <- 2
  }
  else if (nchar(as.character(x[1,(ecs + 1)])) == 4){
    snp_form <- 4
  }
  else{
    stop("Genotypes must be in either 2 (NN) or 4 (0000) format.")
  }

  #functions to get allele and genotype frequencies minus mDating data.
  g_row <- function(row){
    t <- table(unlist(row))
    t <- t[names(t) != mDat]
    #fill table if mDating genotypes
    if(length(t) == 2){
      as <- unique(c(substr(names(t),1,snp_form/2), substr(names(t),(snp_form/2) + 1, snp_form)))
      if(!(paste0(as[1], as[1]) %in% names(t))){
        mDating.geno <- paste0(as[1], as[1])
      }
      else if (!(paste0(as[2], as[2]) %in% names(t))){
        mDating.geno <- paste0(as[2], as[2])
      }
      else if(!(paste0(as[1], as[2]) %in% names(t)) & !(paste0(as[2], as[1]) %in% names(t))){
        mDating.geno <- paste0(as[1], as[2])
      }
      t <- c(t, 0)
      t <- setNames(t, c(names(t)[1:2], mDating.geno))
    }
    else if (length(t) == 1){
      t <- c(t,0,0)
      names(t) <- c(names(t)[1], paste0(rep("X", snp_form/2), collapse = ""), paste0(rep("X", snp_form/2), collapse = ""))
    }
    #organize: homs first then het, with minor hom first
    as <- c(rbind(substr(names(t),1,snp_form/2), substr(names(t),(snp_form/2) + 1, snp_form))) #break down the table names into componenets via interleaving
    g1 <- ifelse(as[1] == as[2], "hom", "het") #determine if first entry is a hom or het
    g2 <- ifelse(as[3] == as[4], "hom", "het")
    g3 <- ifelse(as[5] == as[6], "hom", "het")
    gs <- setNames(c(g1,g2,g3), 1:3)
    gs <- rev(sort(gs)) #sorted by name. actual indices are stored as names
    t <- t[as.numeric(names(gs))] #reorder table, homs then hets
    hom.fs <- sort(c(t[1], t[2])) #sort the homs by frequency
    t <- c(hom.fs, t[3]) #rebind the table
    return(t)
  }

  #allele frequencies from g_row output converted to frequencies
  cat("Getting genotype counts...\n")
  drs <- x[,((ecs + 1):ncol(x))] #get just the data
  gfs <- t(apply(drs,1,g_row)) #get genotypes

  if(test == "permutation"){
    p.vals <- numeric(nrow(gfs))
    gfs.f <- gfs/rowSums(gfs)
    cat("Doing permutation tests.\nGetting allele frequencies...\n")
    afs <- cbind(gfs.f[,1] + (gfs.f[,3]/2), gfs.f[,2] + (gfs.f[,3]/2)) #allele frequencies from genotype frequencies

    #get the expected genotype counts
    gfs.e <- matrix(NA, nrow = nrow(gfs), ncol = 3)
    gfs.e[,1] <- rowSums(gfs)*(afs[,1]^2) #expected minor homozygote
    gfs.e[,2] <- rowSums(gfs)*(afs[,2]^2) #expected major homozygote
    gfs.e[,3] <- rowSums(gfs)*(2*afs[,1]*afs[,2]) #expected heterozygote
    gfs.test <- cbind(gfs, gfs.e)

    #print(gfs.test)
    #print(afs)
    #do chi-squared tests
    cat("Running chi-squared tests:\n")
    for(i in 1:nrow(gfs)){
      if(i %% 1000 == 0){cat("snp number: ", i, "\n")}
      if(any(afs[i,] == 0)){
        p.vals[i] <- 1
      }
      else{
        p.vals[i] <- chisq.test(cbind(gfs.test[i,1:3], gfs.test[i,4:6]), simulate.p.value = TRUE, B = n.reps)$p.value
      }
    }
  }
  else if (test == "exact"){
    cat("Running exact tests...\n")
    p.vals <- unlist(apply(gfs, 1, SNPHWE))
  }
  return(p.vals)
}



#'SNP data analysis with snpR
#'
#'\code{snpR.stats} calls many snpR functions based on the provided arguments to do the majority of the basic statistical process required for SNP analysis.
#'
#'
#'
#'Description of x:
#'    Contains metadata in columns 1:ecs. Remainder of columns contain genotype calls for each individual. Each row is a different SNP, as given by format_snps output options 4 or 6. Requires the column containing the position of the loci in base pairs be named "position". If population info is given, it must be contained in a column named "pop".
#'
#' @param x Data frame. Input data, in either numeric or character formats, as given by format_snps output options 2 or 6.
#' @param ecs Numeric value. Number of extra metadata columns at the start of x.
#' @param stats Character vector, default "basic" (see description). Which stats should be calculated?
#' @param filter Logical, default FALSE. Should the SNP data be filtered?
#' @param smooth Character vector or FALSE, default FALSE. Should variables be smoothed? If so, which?
#' @param bootstrap Character vector or FALSE, default FALSE. Should bootstraps of smoothed windows be conducted? If so, which variables?
#' @param graphs Character vector or FALSE, default FALSE. Which variables should be plotted? Plots \emph{smoothed} values of most statistics save population structure, LD, and Tajima's D. As such, requires the statistics to be set via the "smooth" argument.
#' @param sigma Numeric, default NULL. Smoothing statistic/window size, in kb.
#' @param ws Numeric, default NULL.  How far should the window slide when smoothing? If null, a window is produced centered at every SNP.
#' @param pop Table of length > 1, default NULL. Table containing the number of samples in each population. Samples should be sorted by population and ordered identically in data.
#' @param nk Logical, default TRUE. Should contributions to smoothed windows be weighted by the number of called alleles at each site?
#' @param input Character string, default "NN". How are genotypes noted? Supports "NN" and "0000". See \code{\link{format_snps}}.
#' @param mDat Character string, default "N". How are missing \emph{alleles} coded in x (typically "N" or "00")?
#' @param maf FALSE or numeric between 0 and 1, default FALSE. Minimum acceptable minor allele frequency during filtering.
#' @param hf_hets FALSE or numeric between 0 and 1, default 0.55 Maximum acceptable heterozygote frequency during filtering.
#' @param bi_al Logical, default TRUE. Should non-biallelic SNPs be removed during filtering?
#' @param non_poly Logical, default TRUE. Should non-polymorphic loci be removed during filtering?
#' @param min_ind FALSE or integer, default FALSE. Minimum number of individuals in which a loci must be sequenced during filtering.
#' @param min_loci FALSE or numeric between 0 and 1, default FALSE. Minimum proportion of SNPs an individual must be genotyped at to be kept during filtering.
#' @param pop_maf Logical, default FALSE. For filtering, should loci be kept if they are above the minor allele frequency cuttoff in \emph{any} population?
#' @param filt_re FALSE, "partial", or "full", default "partial". How should loci be re_filtered after individuals are filtered? See \code{\link{filter_snps}}.
#' @param Fst_method Character vector, default "genepop". Which FST estimator should be used? See documentation for options.
#' @param LD_level Character vector or FALSE, default c("group", "pop"). Names of metadata columns by which to split the data for LD calculations.
#' @param LD_chr Character vector, default "all". Names of linkage groups/chromsomes for which pairwise LD values should be calculated at all SNPs. Identifiers must be contained in a column named either "group" or "chr".
#' @param LD_g_par Logical, default NULL. Shoud LD be done in parallel across multiple linkage groups/chromsomes?
#' @param tsd_level Character vector or FALSE, default c("group", "pop"). Names of metadata columns by which to split the data for Tajima's D calculations.
#' @param smooth_levs Character vector or FALSE, default c("group", "pop"). Names of metadata columns by which to split the data for smoothing.
#' @param boot_levs Character vector or FALSE, default c("group", "pop"). Names of metadata columns by which to split the data for bootstraps.
#' @param p.alt Character string, default "two-sided". If bootstraps are done, which alternative hypothesis should be used? See documentation for options.
#' @param num_boots Numeric, default NULL. If bootstraps are done, how many should be performed per stat?
#' @param par_cores Numeric, default FALSE. If parallel processing is requested, how many cores should be used?
#'
#' @return A named list of objects corresponding to requested statistics/graphs/ect.
#'
#' @examples
#'
snpR.stats <- function(x, ecs, stats = "basic", filter = FALSE, smooth = FALSE, bootstrap = FALSE, graphs = FALSE, sigma = NULL, ws = NULL,
                       pop = NULL, nk = TRUE, input = "NN", mDat = "NN", maf = 0.05, hf_hets = 0.55,
                       bi_al = TRUE, non_poly = TRUE, min_ind = FALSE, min_loci = FALSE, pop_maf = FALSE,
                       filt_re = "partial", Fst_method = "genepop", LD_level = c("group", "pop"), LD_chr = "all",
                       LD_g_par = FALSE, tsd_level = c("group", "pop"),
                       smooth_levs = c("group", "pop"), internal_boot_level = "group", boot_split_levels = "pop",
                       p.alt = "two-sided", num_boots = NULL, par_cores = FALSE){

  x <- dplyr::arrange(x, group, position)
  stab <- FALSE #set that there is no saved out.tab currently

  ####################################################
  #figure out which stats we are getting

  #vector of possible stats to request:
  pstats <- c("pi", "ho", "fst", "pa", "ld", "tsd")

  #if the default (basic), set the stats.
  if(stats == "basic"){
    stats <- c("pi", "ho", "fst", "pa")
  }

  #otherwise, check that the options given are acceptable.
  else{
    stats <- tolower(stats) #incase they used caps
    fstats <- which(!(stats %in% pstats)) #any undefined stats?
    if(length(fstats) > 0){
      stop(cat("Incorrect stats defined:", stats[fstats], ". Acceptable stats:", pstats, "\n"))
    }
    #otherwise good to go.
  }

  #how about smoothing?
  if(is.character(smooth)){
    pstats <- pstats[-length(pstats)]
    smooth <- tolower(smooth)
    fsmooth <- which(!(smooth %in% stats) | smooth == "tsd") #any unaccepted stats?
    if(length(fsmooth) > 0){
      stop(cat("Smoothing requested on incorrect variables:", smooth[fsmooth], ". Acceptable stats:", stats[-which(stats == "tsd")], "\n"))
    }
  }

  #pop info, needs more work
  if(is.table(pop)){
    if(length(pop) > 1){
      pop <- list(names(pop), as.numeric(pop))
    }
    else{
      stop("Provided pop table must have more than one population!")
    }
  }
  else if(!is.null(pop)){
    stop("Pop info should be provided as a table, with entries in the same order as samples in x (typically alphabetical)!")
  }


  if(Fst_method == "all"){
    Fst_method <- c("genepop", "Wier", "WC", "Hoh")
  }

  ###################################################
  #initialize output list
  out <- list()

  ###################################################
  #format to NN, fix ecs with filler data if it is too low
  if(ecs < 2){
    nfill <- 2 - ecs
    fm <- matrix(paste0("filler", 1:nrow(x)), nrow(x), nfill)
    fm <- as.data.frame(fm, stringsAsFactors = F)
    colnames(nfill) <- paste0("filler", 1:nfill)
    x <- cbind(fm, x)
    remove(fm)
    ecs <- 2
  }
  else{
    nfill <- FALSE
  }

  if(input != "NN"){
    x <- format_snps(x, ecs = ecs, output = 6, input_form = input, mDat = mDat)
    if(nfill > 0){
      out$ch.form.x <- x[,-c(1:nfill)]
    }
    else{
      out$ch.form.x <- x
    }
  }

  ###################################################
  #filter the data as requested
  if(filter){
    if(pop_maf == TRUE){
      x <- filter_snps(x, ecs, maf, hf_hets, min_ind, min_loci, filt_re, pop, non_poly, bi_al, mDat, out.tab = TRUE)
    }
    else{
      x <- filter_snps(x, ecs, maf, hf_hets, min_ind, min_loci, filt_re, FALSE, non_poly, bi_al, mDat, out.tab = TRUE)
    }
    #pull out the out.tab
    stab <- x[-(names(x) == "x")]
    x <- x$x

    #add to return
    if(nfill > 0){
      out$x.flt <- x[,-c(1:nfill)]
    }
    else{
      out$x.flt <- x
    }
  }


  ###################################################
  #do each stat as requested and save output

  #get ac format, used in most analysis and for nk
  sink(tempfile())
  x_ac <- format_snps(x, ecs, 1, pop = pop, in.tab = stab)
  sink()
  if(is.list(pop) & length(pop) > 1){
    x_ac <- dplyr::arrange(x_ac, pop, group, position)
  }
  else{
    x_ac <- dplyr::arrange(x_ac, group, position)
  }

  #initialize output stat dataframe, also used for smoothing.
  sm.smooth <- smooth[smooth != "fst"]
  sm.in <- matrix(NA, nrow = nrow(x_ac), ncol = ecs + length(sm.smooth) + 1)
  sm.in <- as.data.frame(sm.in)
  sm.in[,1:ecs] <- x_ac[,1:ecs]
  colnames(sm.in) <- c(colnames(x_ac)[1:ecs], sort(sm.smooth), "nk")
  sm.in$nk <- x_ac$n_total
  #add a pop column if info provided...
  if(length(pop) > 1){
    sm.in$pop <- x_ac$pop
  }

  #pi
  if("pi" %in% stats | ("Hoh" %in% Fst_method & "fst" %in% stats)){
    pi <- calc_pi(x_ac)
    sm.in$pi <- pi
  }

  #ho
  if("ho" %in% stats | (("WC" %in% Fst_method | "Wier" %in% Fst_method) & "fst" %in% stats)){
    if(!is.null(pop)){
      ho <- calc_Ho(x, ecs)
    }
    else{
      ho <- calc_Ho(x, ecs, pop = pop)
    }
    sm.in$ho <- ho$ho
  }

  #Fst
  if("fst" %in% stats){
    fst.list <- character(0)
    if("genepop" %in% Fst_method){

      #make the genepop input
      if(file.exists("gp_fst.txt")){
        file.remove("gp_fst.txt")
      }
      gp <- format_snps(x, ecs, 2, pop = pop, outfile = "gp_fst.txt")
      remove(gp)

      #run Fst
      fst <- calc_pairwise_Fst("gp_fst.txt", method = "Genepop", pnames = pop[[1]])
      if(nfill != FALSE){
        fst$loci <- cbind(x[,(nfill+1):ecs], fst$loci)
        out$fst$genepop$snps <- fst$loci
      }
      else{
        fst$loci <- cbind(x[,1:ecs], fst$loci)
        out$fst$genepop$loci <- fst$loci
      }
      out$fst$genepop$overall <- fst$overall
      fst <- fst$loci

      fst.list <- "genepop"

      fst <- reshape2::melt(fst, id.vars = colnames(x)[1:ecs])
      colnames(fst)[colnames(fst) == "value"] <- "fst.genepop"
    }
    if("WC" %in% Fst_method | "Wier" %in% Fst_method){
      #get data ready
      ho <- reshape2::melt(ho, id.vars = colnames(x_ac)[1:ecs])
      ho$variable <- as.character(ho$variable)
      x_ac$Ho <- ho$value

      if("Wier" %in% Fst_method){
        fst.Wier <- calc_pairwise_Fst(x_ac, ecs, method = "Wier")
        out$fst$Wier <- fst.Wier
        fst.list <- c(fst.list, "Wier")
        fst.Wier <- reshape2::melt(fst.Wier, id.vars = colnames(x)[1:ecs])
        if(exists("fst")){
          fst <- cbind(fst, fst.Wier = fst.Wier$value)
        }
        else{
          fst <- fst.Wier
        }
      }
      if("WC" %in% Fst_method){
        fst.WC <- calc_pairwise_Fst(x_ac, ecs, method = "WC")
        out$fst$WC <- fst.WC
        fst.list <- c(fst.list, "WC")
        fst.WC <- reshape2::melt(fst.WC, id.vars = colnames(x)[1:ecs])
        if(exists("fst")){
          fst <- cbind(fst, fst.WC = fst.WC$value)
        }
        else{
          fst <- fst.WC
        }
      }
    }
    if("Hoh" %in% Fst_method){ #hohenlohe
      x_ac$pi <- pi
      fst.Hoh <- calc_pairwise_Fst(x_ac, ecs, method = "Hohenlohe")
      out$fst$Hoh <- fst.Hoh
      fst.list <-c(fst.list, "Hoh")
      fst.Hoh <- reshape2::melt(fst.Hoh, id.vars = colnames(x)[1:ecs])
      if(exists("fst")){
        fst <- cbind(fst, fst.Hoh = fst.Hoh$value)
      }
      else{
        fst <- fst.Hoh
      }
    }


    if("fst" %in% smooth){
      if(nk){
        pnk <- calc_pairwise_nk(x_ac, ecs)
      }
      pnk <- reshape2::melt(pnk, id.vars = colnames(x)[1:ecs])
      fst$nk <- pnk$value
      colnames(fst)[which(colnames(fst) == "variable")] <- "pop"
      out$fst$smooth <- s_ave_multi(fst,
                                    colnames(fst)[which(grepl("fst", colnames(fst)))],
                                    sigma, ws, nk, smooth_levs)
    }
  }

  #private alleles
  if("pa" %in% stats){
    pa <- check_private(x_ac)
    pa <- cbind(x[,1:ecs], pa)
    pa <- reshape2::melt(pa, id.vars = colnames(x_ac[,1:ecs]))
    colnames(pa)[colnames(pa) == "variable"] <- "pop"
    colnames(pa)[colnames(pa) == "value"] <- "pa"
    sm.in$pa <- pa$pa
  }

  #tsd
  if("tsd" %in% stats){
    warning("Tajima's D should be run on a full dataset, without filtering SNPs and containing non-segregating genotypes.")

    #are we splitting tsd by pop?
    if("p" %in% tsd_level){

      #are we also splitting by group?
      if("g" %in% tsd_level){
        tsd <- run_gp(x_ac, Tajimas_D, ws = ws, step = step)
      }

      else{
        tsd <- Tajimas_D(x_ac[x_ac$pop == pop[[1]][1],], ws, step)
        for(i in 2:length(pop[[1]])){
          tsd <- rbind(tsd, Tajimas_D(x_ac[x_ac$pop == pop[[1]][i],], ws, step))
        }
      }

    }

    #splitting by just group
    else if("g" %in% tsd_level){
      tsd <- run_g(x_ac, Tajimas_D, ws = ws, step = step)
    }


    #not splitting?
    else{
      tsd <- Tajimas_D(x_ac, ws, step)
    }

    out$tsd <- tsd
  }


  #LD, this one might be slow...
  if("ld" %in% stats){

    #are we running only certain groups/chr?
    if(LD_chr != "all"){
      xLD <- x[which(x$group %in% LD_chr),]
    }
    else{
      xLD <- x
    }

    #what level are we running this at?
    #are we splitting by pops?
    if("pop" %in% LD_level){

      #initialize pop output list
      pLD <- vector("list", length = length(pop))
      names(pLD) <- pop[[1]]

      #initialize column tracker
      tracker <- ecs + 1

      for(i in 1:length(pop)){
        pLD[[i]] <- LD_full_pairwise(x = cbind(LDx[,1:ecs], LDx[, tracker:(tracker + pop[[2]][i] - 1)]),
                                     ecs, levels = LD_level[-which(LD_level == "pop")], par = par_cores)
        tracker <- tracker + pop[[2]][i]
      }
    }


    #we aren't splitting by pops?
    else{
      pLD <- LD_full_pairwise(LDx, ecs, par = par_cores, levels = LD_level)
    }
    out$LD <- pLD
  }

  out$stats <- sm.in

  ###########################################
  #smooth any other requested stats

  browser()
  if(is.character(smooth)){
    sm.out <- s_ave_multi(sm.in, colnames(sm.in)[(ecs+1):(ecs + length(sm.smooth))], sigma, ws, nk, smooth_levs)
    out$stats_smoothed <- sm.out
  }

  ###########################################
  #bootstrap any requested stats
  browser()
  #if we are bootstrapping...
  if(bootstrap != FALSE){

    #If all fst. methods were requested, need to seperate these out.
    if(any(grepl("fst.all", bootstrap))){
      bootstrap <- c(bootstrap[-which(grepl("fst", bootstrap))], "fst.Hoh", "fst.genepop", "fst.WC", "Fst.Wier")
    }

    #make a matrix detailing which runs to do exactly:
    if(boot_split_levels == FALSE){
      boot_matrix <- matrix(bootstrap, ncol = 1)
    }
    else{
      boot_matrix <- matrix(NA, ncol = 1 + length(boot_split_levels))
      for(i in 1:length(bootstrap)){
        if(grepl("fst", bootstrap[i])){
          ul <- matrix(as.character(unique(fst[,colnames(fst) %in% boot_split_levels])), ncol = length(boot_split_levels))
          boot_matrix <- rbind(boot_matrix, cbind(bootstrap[i], ul))
        }
        else{
          ul <- matrix(as.character(unique(out$stats[,colnames(out$stats) %in% boot_split_levels])), ncol = length(boot_split_levels))
          boot_matrix <- rbind(boot_matrix, cbind(bootstrap[i], ul))
        }
      }
      boot_matrix <- boot_matrix[-1,]
    }

    #prepare things for parallel run if requested:
    if(par_cores != FALSE){
      cl <- snow::makeSOCKcluster(par_cores)
      doSNOW::registerDoSNOW(cl)

      ntasks <- nrow(boot_matrix)
      progress <- function(n) cat(sprintf("Part %d out of",n), ntasks, "is complete.\n")
      opts <- list(progress=progress)
    }

    #function to grab correct input data:
    gbdat <- function(tstat){
      if(grepl("fst", tstat)){

        #grab the correct data
        if(tstat == "fst.genepop"){
          tx <- out$fst$genepop$loci
        }
        else{
          tx <- out$fst[[which(grepl(substr(tstat, 5, nchar(tstat)), names(out$fst)))]]
        }

        #melt it to split pops
        tx <- reshape2::melt(tx, id.vars = colnames(tx)[1:ecs])
        colnames(tx)[-c(1:ecs)] <- c("comp", tstat)

        #get fws
        if(!is.null(ws)){
          tfws <- split(data.table::data.table(out$fst$smooth), by = boot_levs, sorted = TRUE)
        }
        else{
          tfws <- NULL
        }

        return(list(dat = tx, fws = tfws))
      }
      else{
        tdat <- out$stats
        tfws <- ifelse(is.null(ws), NULL, out$stats_smoothed)
        return(list(dat = tdat, fws = fws))
      }
    }

    #call if parallel:
    if(par_cores != FALSE){
      foreach::foreach(q = 1:ntasks, .inorder = TRUE,
                       .options.snow = opts) %dopar% {

                         tdat <- gbdat(tstat = boot_matrix[i,1])
                         tfws <- tdat$fws
                         tdat <- tdat$dat

                         if(ncol(boot_matrix) > 1){
                           #get only the correct data
                           tdat <- tdat[which(apply(tdat, 1, function(x) identical(x[which(colnames(x) %in% boot_split_levels),], boot_matrix[i,-1]))),]
                           tfws <- tfws[which(apply(tfws, 1, function(x) identical(x[which(colnames(x) %in% boot_split_levels),], boot_matrix[i,-1]))),]
                         }

                         b_output <- resample_long(tdat, boot_matrix[i,1], num_boots, sigma, nk, fws, TRUE, level = internal_boot_level)
                       }
    }
    else{
      b_output <- vector("list", length = nrow(boot_matrix))
      for(i in 1:nrow(boot_matrix)){

        tdat <- gbdat(tstat = boot_matrix[i,1])
        tfws <- tdat$fws
        tdat <- tdat$dat

        if(ncol(boot_matrix) > 1){
          #get only the correct data
          tdat <- tdat[which(apply(tdat, 1, function(x) identical(x[which(colnames(x) %in% boot_split_levels),], boot_matrix[i,-1]))),]
          tfws <- tfws[which(apply(tfws, 1, function(x) identical(x[which(colnames(x) %in% boot_split_levels),], boot_matrix[i,-1]))),]
        }

        b_output[[i]] <- resample_long(tdat, boot_matrix[i,1], num_boots, sigma, nk, fws, TRUE, level = internal_boot_level)
      }
    }

    #process output:

  }


  ###########################################
  #return the object
  return(out)

}




#function to run any command after spliting the data by group and by population. Group and population must be in columns named "group" and "pop", respectively.
#The first argument of the function to run must be the data provided to that function.
run_gp <- function(x, FUN, ...){
  w_df<- data.frame()
  x$pop <- as.character(x$pop)
  x$group <- as.character(x$group)
  pops <- unique(x[,"pop"])
  groups <- unique(x[,"group"])
  for (i in 1:length(groups)){
    w_data <- x[x[,"group"] == groups[i],]
    print(groups[i])
    for(j in 1:length(pops)){
      print(pops[j])
      x_data <- FUN(w_data[w_data[,"pop"] == pops[j],], ...)
      if(length(x_data) == 0){next}
      if(is.data.frame(x_data)){
        if(!any(colnames(x_data) == "pop"))
          x_data <- cbind(pop = pops[j], x_data, stringsAsFactors = F)
        if(!any(colnames(x_data) == "group")){
          x_data <- cbind(group = groups[i], x_data, stringsAsFactors = F)
        }
        w_df <- rbind(w_df, x_data)
      }
      else{
        if(is.data.frame(w_df)){
          w_df <- numeric(0)
        }
        w_df <- c(w_df,x_data)
      }
    }
  }
  if(!is.data.frame(w_df)){
    warning("Output is vector: if merging with meta info, check that info data frame is sorted identically:
            Group, Pop, then snp by snp!")
  }
  if(sum(colnames(w_df) == "pop") > 1){
    excols <- which(colnames(w_df) == "pop")[-1]
    w_df <- w_df[,-excols]
  }
  return(w_df)
}

#function to run any command after spliting the data by group. Group must be in a column named "group".
#The first argument of the function to run must be the data provided to that function.
run_g <- function(x, FUN, ..., y = NULL){
  w_df<- data.frame()
  x$group <- as.character(x$group)
  groups <- unique(x[,"group"])
  for (i in 1:length(groups)){
    w_data <- x[x[,"group"] == groups[i],]
    if(!is.null(y)){
      y_dat <- y[y[,"group"] == groups[i],]
      w_data <- FUN(w_data, y = y_dat, ...)
    }
    else{
      w_data <- FUN(w_data, ...)
    }
    print(groups[i])
    if(length(w_data) == 0){next}
    if(is.data.frame(w_data)){
      if(!any(colnames(w_data) == "group")){
        w_data <- cbind(group = groups[i], w_data, stringsAsFactors = F)
      }
      w_df <- rbind(w_df, w_data)
    }
    else{
      if(is.data.frame(w_df)){
        w_df <- numeric(0)
      }
      w_df <- c(w_df, w_data)
    }
  }
  if(!is.data.frame(w_df)){
    warning("Output is vector: if merging with meta info, check that info data frame is sorted identically:
            group then snp by snp!\n")
  }
  return(w_df)
}


#'Get random snps each window.
#'
#'\code{rgap_snps} grabs n random SNPs every window of a given length. Can do this across multiple linkage groups/chromosomes.
#'
#'Possible methods:
#'\itemize{
#'    \item{window: }{Takes n random SNPs every window. This can result in SNPs that are much nearer to each other than the gap length.}
#'    \item{gap: }{Takes the first possible SNP, then takes the next SNP at least gap bases away. Non-random.}
#'}
#'
#' @param x data.frame. Input data, must contain a column titled 'position' containing genomic positions of the SNPs.
#' @param gap Integer. How big are the windows from which to sample SNPs?
#' @param method Character, default "window". By which method should the SNPs be chosen?
#' @param n Integer, default 1. For "window" method. How many SNPs should be taken per window?
#' @param levels Character, default "group". Across which levels should the windows be partitioned (ex. "chromosome")?
#'
#' @return A data.frame containing subset SNPs.
#'
rgap_snps <- function(x, gap, method = "gap", n = 1, levels = "group"){

  if(method == "window"){
    gfun <- function(x){
      cs <- seq(gap/2, max(x$position) + gap/2, gap)
      starts <- cs - gap/2
      ends <- cs + gap/2
      pos <- x$position

      lmat <- outer(pos, starts, function(pos, starts) pos >= starts)
      lmat <- lmat + outer(pos, ends, function(pos, ends) pos < ends)
      colnames(lmat) <- cs
      rownames(lmat) <- pos
      lmat <- ifelse(lmat == 2, TRUE, FALSE)

      wins <- which(t(lmat) == TRUE) %% length(cs)
      wins[wins == 0] <- length(cs)

      x <- cbind(window = wins, x)

      rands <- x %>% group_by(window) %>% sample_n(size = n)

      rands <- rands[,-1]

      return(rands)
    }
  }
  else if (method == "gap"){
    gfun <- function(x){
      fsnp <- which.min(x$position)
      pos <- x$position
      cs <- pos[fsnp]
      #figure out windows
      while(cs[length(cs)] + gap <= max(pos)){
        cs <- c(cs,
                min(pos[pos >= cs[length(cs)] + gap]))
      }
      ret <- match(cs, pos)
      return(x[ret,])
    }
  }

  if(levels != FALSE){
    ret <- plyr::ddply(x, .variables = "group", .fun = gfun, .progress = "text")
  }
  else{
    ret <- gfun(x)
  }
  return(ret)
}
