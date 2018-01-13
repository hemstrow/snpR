#function to calc pi from data. should be in format with num alleles, total count of alleles,
#and subsequent alleles counts in columns named "n_aleles", "n_total", and "ni1", "ni2", and so on (allele count/bayescan format, as given by format_snps option one).
#Returns a VECTOR of pi values.

#'Calculate PI from SNP data.
#'
#'\code{calc_pi} Calculates pi (genetic diversity/average number of pairwise differences) according to Hohenlohe et al. (2010) from SNP data. Returns a vector of pi values for each SNP as sorted in input.
#'
#'Description of x:
#'    Must contain colums containing the number of *unique* alleles, total count of alleles sequenced in all individuals, and subsequent alleles counts for each observed allele in columns named "n_alleles", "n_total", "ni1", and "ni2". This matches the allele count/bayescan format as given by format_snps option one. Calculates pi for each row of data, and can therefore also contain a "pop" column and any other metadata columns such as SNP position.
#'
#' @param x Input SNP data, in allele count format as given by format_snps output option
#'
#' @return A vector of PI values, in the same order as the input SNPs.
#'
#' @examples
#' #get pop info, convert character data into the appropriate format, calculate pi, and add metadata.
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' l <- list(c(names(pops)), as.numeric(pops))
#' ac <- format_snps(stickSNPs, 3, pop = l)
#' pi <- calc_pi(ac)
#' pi <- cbind(ac[,1:4], pi)
#'
calc_pi <- function(x){
  nt <- as.numeric(x[,"n_total"])
  n1 <- as.numeric(x[,"ni1"])
  n2 <- as.numeric(x[,"ni2"])
  if(!all(as.numeric(x[,"n_alleles"]) != 2)){
    warning("Some loci do not have two alleles in some populations.\n")
  }
  binom_n1 <- choose(n1,2)
  binom_n2 <- choose(n2,2)
  binom_nt <- choose(nt,2)
  #print(n1)
  #print(n2)
  p <- 1 - (binom_n1 + binom_n2)/binom_nt
  #print(pi)
  return(p)
}


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

#Calculates Tajima's theta/pi, Waterson's theta, and Tajima's D over a sliding window. Takes the bayescan/allele count format. Pi calculated as in hohenlohe 2010.
#inputs: x: input data, first three columns are: snp ID, position, snp scaffold/lg/chr. Named required columns are "n_total", "ni1", and "ni2". Named column "pop" is optional, and will be returned in output if exists.
#        ws: window size in kbp
#        step: step size in kbp
#        report: give progress every "report" windows.
#output: A data frame where column one is group, two is position, three is Tajima's theta,
#        four is Waterson's theta, and five is D. Note that both thetas are reported as a
#        frequency (so the number per base). Full number per window is ws*theta.

#'Tajima's D from SNP data.
#'
#'Calculates Tajima's theta/pi, Waterson's theta, and Tajima's D over a sliding window. Pi calculated as in Hohenlohe et al. 2010. Tajima's D is calculated using the formula from Tajima (1989) Statistical Method for Testing the Neutral Mutation Hypothesis by DNA Polymorphism.
#'
#'Description of x:
#'    Must contain colums containing the number of *unique* alleles, total count of alleles sequenced in all individuals, and subsequent alleles counts for each observed allele in columns named "n_alleles", "n_total", "ni1", and "ni2". Also needs a column containing the position of each SNP, in bp. This matches the allele count/bayescan format as given by format_snps option one.
#'
#' @param x Input data, in allele count format as given by format_snps output option 1.
#' @param ws The size of each window, in kb.
#' @param step Lenght to slide between each window, in kb.
#' @param report When reporting progress, how many windows should be calculated between each report?
#'
#' @return Data frame containing metadata, Waterson's Theta, Tajima's Theta, and Tajima's D for each window.
#'
#' @examples
#' #convert example data into the correct format and run splitting by groups and pops.
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' l <- list(c(names(pops)), as.numeric(pops))
#' ac <- format_snps(stickSNPs, 3, pop = l)
#' run_gp(ac, Tajimas_D, 200, 50)
#'
Tajimas_D <- function(x, ws, step, report = 20){

  #windowing
  out <- data.frame() #initialize output
  tps <- sort(x$position) #get the snp positions, sort
  lsp <- tps[length(tps)] #get the position of the last snp to use as endpoint
  c <- 0 #set starting position
  i <- 1 #set starting iteration for writing output
  ws <- 1000*ws
  while (c <= lsp){
    start <- c - ws #change window start
    end <- c + ws #change window end
    if(i %% report == 0){cat("Window Postion:", c, "out of", lsp, "\n")}

    #take all the snps in the window, calculate T's theta, W's theta, and T's D
    wsnps <- x[x$position <= end & x$position >= start,] #get only the snps in the window
    wsnps <- wsnps[wsnps[,"ni1"] > 0 & wsnps[,"ni2"] > 0,] #remove any snps that are fixed in this pop/group/whatever
    if(nrow(wsnps) == 0){ #if no snps in window
      c <- c + step*1000 #step along
      next #go to the next window
    }
    b1s <- choose(wsnps[,"ni1"],2) #binomial for first allele
    b2s <- choose(wsnps[,"ni2"],2) #binomial for second allele
    bts <- choose(wsnps[,"n_total"],2) #binomial for all alleles
    ts.theta <- sum(1-((b1s+b2s)/bts)) #average number of pairwise differences (pi) per snp. Equivalent to sum((ndif/ncomp)) for all snps
    ts.thetaf <- ts.theta/ws #pi for the whole region, which includes all of the non-polymorphic sites. Reason why this shouldn't be run with a maf filter, probably.
    n_seg <- nrow(wsnps) #number of segregating sites
    K <- round(mean(wsnps[,"n_total"])) #average sample size for ws.theta, as in hohenlohe 2010. Alternative would make this into a function, then use lapply on the vector of K values
    #if(is.nan(K) == TRUE){browser()}
    a1 <- sum(1/seq(1:(K-1))) #get a1
    ws.theta <- n_seg/a1 #get ws.theta
    ws.thetaf <- ws.theta/ws #ws.theta fraction

    #get T's D by part. See original paper, Tajima 1989.

    a2 <- sum(1/(seq(1:(K-1))^2))
    b1 <- (K+1)/(3*(K-1))
    b2 <- (2*(K^2 + K + 3))/(9*K*(K-1))
    c1 <- b1 - (1/a1)
    c2 <- b2 - ((K+2)/(a1*K)) + (a2/(a1^2))
    e1 <- c1/a1
    e2 <- c2/(a1^2 + a2)
    D <- (ts.theta - ws.theta)/sqrt((e1*n_seg) + e2*n_seg*(n_seg - 1))

    #output result for this window, step to the next window

    out[i,"group"] <- x[1,3]
    if("pop" %in% colnames(x)){
      out[i,"pop"] = x[1,"pop"] #if a pop column is in the input, add a pop column here.
    }
    out[i,"position"] <- c
    out[i,"ws.theta"] <- ws.thetaf
    out[i,"ts.theta"] <- ts.thetaf
    out[i,"D"] <- D
    c <- c + step*1000
    i <- i + 1
  }
  return(out)
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
HWE <- function(x, ecs, miss = "NN", test = "exact", n.reps = 20000){


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

  #functions to get allele and genotype frequencies minus missing data.
  g_row <- function(row){
    t <- table(unlist(row))
    t <- t[names(t) != miss]
    #fill table if missing genotypes
    if(length(t) == 2){
      as <- unique(c(substr(names(t),1,snp_form/2), substr(names(t),(snp_form/2) + 1, snp_form)))
      if(!(paste0(as[1], as[1]) %in% names(t))){
        missing.geno <- paste0(as[1], as[1])
      }
      else if (!(paste0(as[2], as[2]) %in% names(t))){
        missing.geno <- paste0(as[2], as[2])
      }
      else if(!(paste0(as[1], as[2]) %in% names(t)) & !(paste0(as[2], as[1]) %in% names(t))){
        missing.geno <- paste0(as[1], as[2])
      }
      t <- c(t, 0)
      t <- setNames(t, c(names(t)[1:2], missing.geno))
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
        if(is.na(p.vals[i])){browser()}
      }
    }
  }
  else if (test == "exact"){
    cat("Running exact tests...\n")
    p.vals <- unlist(apply(gfs, 1, SNPHWE))
  }
  return(p.vals)
}


#Calculates pairwise fst for each pair of populations according to hohenlohe (2010). Input format is
#snp, position, group, pop, total allele count (total), number of alleles (num), ni1, ni2, and pi (diversity, from calc_pi).
#Automatically sorts data by group, position, and population.
calc_pairwise_Fst <- function(x, ecs, do.nk = FALSE, skip.FST = FALSE){
  if(!do.nk & skip.FST){
    stop("Must specify either pairwise FST, nk, or both")
  }
  if(!all(x$n_alleles == 2)){
    vio <- which(x$n_alleles != 2)
    cat("Some loci do not have two alleles. Violating loci:\n")
    print(vio)
    stop()
  }
  i <- 1
  colnames(x) <- c("snp", "position", "group", "pop", "ntot", "num", "ni1", "ni2", "pi")
  x[,"pop"] <- as.character(x[,"pop"])
  x[,"group"] <- as.character(x[,"group"])
  x <- arrange(x, group, position, pop)
  pops <- sort(unique(x[,"pop"]))
  out <- as.data.frame(matrix(NA, ncol = ecs+(length(pops)*(length(pops) - 1)/2), nrow = nrow(x)/length(pops)))

  #print(out)
  #initialize pop comparison columns.
  comps <- c()
  while (i < (length(pops))){
    j <- 1 + i
    for (j in j:length(pops)){
      comps <- c(comps, paste0(pops[i], "_", pops[j]))
      j <- j + 1
    }
    i <- i + 1
  }
  i <- 1
  colnames(out) <- c("snp", "position", "group", comps)
  out[,1:ecs] <- x[x$pop == pops[1],1:ecs] #add snp, position, group to output
  if(do.nk){ #prepare nk output if requested
    nout <- out
    colnames(nout)[4:ncol(nout)] <- paste0("nk_", comps)
  }

  #print(out)

  #loop through each comparison and caculate pairwise FST at each site
  c.col <- 4 #set starting column for pasting data
  for (i in 1:(length(pops) - 1)){ #i is the first pop
    idat <- x[x$pop == pops[i],] #get data for first pop
    j = i + 1 #intialize j as the next pop
    for (j in j:length(pops)){#j is pop being compared
      if(!skip.FST){
        jdat <- x[x$pop == pops[j],] #get data for second pop
        #n.ali <- ifelse(idat$ni1 != 0, 1, 0) + ifelse(idat$ni2 != 0, 1, 0) #get number of alleles in pop 1
        #n.alj <- ifelse(jdat$ni1 != 0, 1, 0) + ifelse(jdat$ni2 != 0, 1, 0) #get number of alleles in pop 2
        #com.top <- (choose(n.ali, 2) * idat$pi) + (choose(n.alj, 2) * jdat$pi) #get the numerator
        com.top <- (choose(idat$n_total, 2) * idat$pi) + (choose(jdat$n_total, 2) * jdat$pi)

        t.ni1 <- idat$ni1 + jdat$ni1 #get the total number of allele one
        t.ni2 <- idat$ni2 + jdat$ni2 #get the total number of allele two
        ptop <- choose(t.ni1, 2) + choose(t.ni2, 2) #get the pooled pi numerator
        pbot <- choose((t.ni1 + t.ni2), 2) #get the pooled pi denominator
        ppi <- 1 - ptop/pbot #get pooled pi
        #com.bot <- ppi * (choose(n.ali,2) + choose(n.alj,2)) #get the denominator
        com.bot <- ppi * (choose(idat$n_total,2) + choose(jdat$n_total,2)) #get the denominator
        Fst <- 1 - com.top/com.bot #get fst
        if(any(abs(Fst) >= 1, na.rm = T)){browser()}
        #Fst[t.ni1 == 0 | t.ni2 == 0] <- 0 #could uncomment this if want 0s instead of NaNs.
        out[,c.col] <- Fst #write fst
      }
      if(do.nk){
        nout[,c.col] <- idat$n_total + jdat$n_total
      }
      c.col <- c.col + 1 #agument c.col
    }
  }
  if(do.nk){
    if(skip.FST){
      return(nout)
    }
    else{
      return(list(FST = out, nk = nout))
    }
  }
  else if (!skip.FST){
    return(out)
  }
}

#wrapper for just doing pairwise.nk
calc_pairwise_nk <- function(x, ecs){
  out <- calc_pairwise_Fst(x, ecs, do.nk = TRUE, skip.FST = TRUE)
  return(out)
}

#Calculates observed heterozygosity at a snp.
#Inputs:  data: Input data. In numeric or NN character format, see format_snps output options 4 or 6.
#               Nucleotides coded as 01 to 04 or A,T,C,G, respectively.
#         ecs: Number of columns containing metadata.
#         m.dat: Missing data indicator. Should be in the same format as data.
#         pop: List with population information for individuals. Format is as produced by:
#              list(c("North", "South", "East", "West"), c(10,20,30,40)). First vector is names of pops,
#              second vector is the count of each pop. Input data MUST be in the same order as this list.
calc_Ho <- function(x, ecs, m.dat = "NN", pop = NULL){

  #set possible heterozygotes
  if(nchar(x[1,(ecs + 1)]) == 4 & nchar(m.dat) == 4){
    hl <- c(m.dat, "0101", "0202", "0303", "0404")
  }
  else if (nchar(x[1,(ecs + 1)]) == 2 & nchar(m.dat) == 2){
    hl <- c(m.dat, "AA", "TT", "CC", "GG")
  }
  else{
    stop("Data and missing signifier must be in either one (NN) or two character (0000) per allele format.")
  }

  #initalize output
  if(!is.list(pop)){
    pop = list("pi", ncol(x) - (ecs + 1) + 1)
  }
  pns <- pop[[1]]
  psz <- pop[[2]]
  out <- matrix(NA, nrow(x), length(pns))
  out <- cbind(x[,1:((ecs + 1) - 1)], out)
  colnames(out) <- c(colnames(x)[1:((ecs + 1) - 1)], pns) #set output column names


  #do each pop, need to loop here.

  c.col <- (ecs + 1)
  for (j in 1:length(pns)){
    wdat <- x[,c.col:(c.col+psz[j] - 1)]
    #with this data, figure out heterozygosity
    het.c <- rowSums(ifelse(wdat == hl[1] | wdat == hl[2]
                            | wdat == hl[3] | wdat == hl[4]
                            | wdat == hl[5], 0, 1))
    ho <- het.c/rowSums(ifelse(wdat == m.dat, 0, 1))
    out[,(ecs + 1) + j - 1] <- ho
    c.col <- c.col + psz[j]
  }
  return(out)
}

#Checks for private alleles.
#inputs: x: data, in allele count format such as that given by format_snps option 1. Expects
#           columns named "pop", "ni1", and "ni2", which contain pop designations, allele one counts,
#           and allele 2 counts.
check_private <- function(x){
  l <- unique(x$pop) #gets unique pops, requires column named pop

  a1m <- matrix(NA, nrow(x)/length(l), length(l)) #initialize a1 storage
  a2m <- matrix(NA, nrow(x)/length(l), length(l)) #initialize a2 storage
  nloc <- nrow(x)/length(l)

  #loop through pops load a1s and a2s
  count <- 1
  for(i in 1:length(l)){
    a1m[,i] <- x$ni1[count:(count + nloc - 1)]
    a2m[,i] <- x$ni2[count:(count + nloc - 1)]
    count <- count + nloc
  }

  #convert to presence absence
  a1m <- ifelse(a1m != 0, 1, 0)
  a2m <- ifelse(a2m != 0, 1, 0)

  #convert to private/not private
  a1m <- ifelse(rowSums(a1m) == 1 & a1m == 1, 1, 0)
  a2m <- ifelse(rowSums(a2m) == 1 & a2m == 1, 1, 0)

  #combine a1 and a2
  pa <- a1m + a2m

  #return data
  colnames(pa) <- paste0("pa_", l)
  return(as.data.frame(pa))
}

#Calculates Dprime, rsq, and a p-value for LD for each pair of snps.
#inputs: x: data, in either numeric or NN form. Must contain a column named "position", can contain a column named "group" and/or "pop"
#        ecs: number of extra columns (headers) before the data in x
#        prox_table: Should a proximity table be output?
#        matrix_out: Should LD matrices be created?
#        mDat: What is the missing data character? Expects marker for a single allele. ("N" or "01")
LD_full_pairwise <- function(x, ecs, prox_table = TRUE, matrix_out = TRUE, mDat = "N") {
  #new LD pairwise function. Loops through each loci, but vectorized for each comparison within that loci

  #get data and headers
  headers <- x[,1:ecs]
  x <- x[,(ecs + 1):ncol(x)]
  x <- as.matrix(x)

  #get missing genotype
  dmDat <- paste0(mDat, mDat)

  #data format
  sform <- nchar(x[1,1])/2

  #get unique alleles present at each locus (note: this made a much quicker a tab... should implement this elsewhere...)
  p1 <- substr(as.matrix(x), 1, sform)
  p2 <- substr(as.matrix(x), sform + 1, sform*2)
  pc <- melt(cbind(p1,p2))
  #amat <- with(pc, table(Var1, value))
  #amat <- amat[,colnames(amat) != mDat]
  as <- sort(unique(pc$value))
  as <- as.character(as[as != mDat])


  #need to loop through each loci and compare to everything else. Probably can't really vectorize the outer loop.
  #goal: make a haplotype table, where each row is a comparison and each column is a haplotype count
  #to count haplotypes: Double heterozgote (AC CG) = mark neither.
  #                     Double homozygote (AA CC): mark two of the combination (A with C)
  #                     homo/het (AA CG): mark one of each combination (A with C and A with G)
  #
  #function to do this need to: 1) paste together the observed genotypes of all observed genotype combinations.
  #                             2) convert this to a table of counts of these genotypes for each pairwise combo.
  #                             3) clean the table
  #                             4) get haplotype counts.
  #
  #make a function to generate this table given a starting locus:
  # inputs: x: row containing genotypes at starting locus
  #         y: rows containing genotypes at all comparison loci
  tabulate_haplotypes <- function(x, y){
    #1)
    #get the observed genotype combinations
    yv <- as.vector(t(y))
    gcv <- paste0(x, yv)
    if(!is.matrix(y)){
      gcv <- matrix(gcv, 1, length(x), byrow = T)
    }
    else{
      gcv <- matrix(gcv, nrow(y), length(x), byrow = T)
    }


    #2)
    #turn this into a genotype count table
    mgcv <- melt(gcv)
    ghapmat <- with(mgcv, table(Var1, value))

    #3) clean the table
    ##grab column names
    gl <- colnames(ghapmat)
    ##remove anything with missing data and double hets
    rgcs <- c(grep(paste0("^", dmDat), gl), #missing first locus
              grep(paste0(dmDat, "$"), gl), #missing second locus
              which(substr(gl, 1, sform) != substr(gl, (sform + 1), (sform *2)) &
                      substr(gl, (sform*2) + 1, sform*3) != substr(gl, (sform*3+1), sform*4))) #double het

    ##remove any double heterozygotes
    ghapmat <- ghapmat[,-rgcs]

    #add a filler row for the last pairwise comparison to make life easier.
    if(!is.matrix(y)){
      if(length(ghapmat) > 1){ #stop it from doing this if there is data for only one haplotype.
        ghapmat <- rbind(ghapmat, rep(c(10,0), 100)[1:length(ghapmat)])
      }
    }

    #if nothing remains, return nothing
    if(length(ghapmat) == 0){
      return(NA)
    }
    else if(!is.matrix(ghapmat)){ #if only one column remains...
      if(substr(gl[-rgcs], 1,sform) == substr(gl[-rgcs], sform + 1, sform*2)){ #double hom
        return(NA)
      }
      else{
        return(c(D = 0))
      }
    }

    #4) get hap table. Use the rules above to do this. Possible conditions, where alleles at locus 1 = 1a1b and locus 2 = 2a2b
    if(!is.vector(colnames(ghapmat))){browser()}
    ghapmat <- ghapmat[,order(colnames(ghapmat))] #put in order, just in case
    gl <- colnames(ghapmat) #get column names again
    hnames <- c(paste0(as, as[1]),
                paste0(as, as[2]),
                paste0(as, as[3]),
                paste0(as, as[4]))
    if(any(grepl("NA", hnames))){ #in the case of a missing allele...
      hnames <- hnames[-grep("NA",hnames)]
    }
    hnames <- sort(hnames)
    hapmat <- matrix(0, nrow(ghapmat), length(hnames))#initialize. all possible haplotypes
    colnames(hapmat) <- hnames


    ##figure out which are homozygotes and heterozygotes at either locus in the pairwise comparison.
    dhom <- substr(gl, 1, sform) == substr(gl, sform + 1, sform*2) &
      substr(gl, (sform*2) + 1, sform*3) == substr(gl, (sform*3+1), sform*4) #columns with double homozygotes
    dhom <- ghapmat[,dhom]
    het_l1 <- substr(gl, 1, sform) != substr(gl, sform + 1, sform*2) #columns where the first locus is het
    het_l1 <- ghapmat[,het_l1]
    het_l2 <- substr(gl, (sform*2) + 1, sform*3) != substr(gl, (sform*3+1), sform*4) #colunms where the second locus is het
    het_l2 <- ghapmat[,het_l2]

    #fix wierd cases where one of these isn't a matrix because only one haplotype falls into the category.
    if(any(!is.matrix(dhom), !is.matrix(het_l1), !is.matrix(het_l2))){
      if(!is.matrix(dhom)){
        dhom <- as.matrix(dhom)
        colnames(dhom) <- colnames(ghapmat)[substr(gl, 1, sform) == substr(gl, sform + 1, sform*2) &
                                              substr(gl, (sform*2) + 1, sform*3) == substr(gl, (sform*3+1), sform*4)] #columns with double homozygotes
      }
      if(!is.matrix(het_l1)){
        het_l1 <- as.matrix(het_l1)
        colnames(het_l1) <- colnames(ghapmat)[substr(gl, 1, sform) != substr(gl, sform + 1, sform*2)] #columns where the first locus is het
      }
      if(!is.matrix(het_l2)){
        het_l2 <- as.matrix(het_l2)
        colnames(het_l2) <- colnames(ghapmat)[substr(gl, (sform*2) + 1, sform*3) != substr(gl, (sform*3+1), sform*4)]
      }
    }

    #count up the haplotypes.
    #function to correct haplotype input matrix:
    GtoH <- function(x, n){
      m1 <- matrix(as.numeric(x), nrow(x), ncol(x))
      colnames(m1) <- colnames(x)
      m1 <- cbind(as.data.frame(t(m1)), n)
      m2 <- m1 %>% group_by(n) %>% summarise_all(funs(sum))
      m2 <- t(as.matrix(m2[,-1]))
      colnames(m2) <- sort(unique(n))
      return(m2)
    }
    ##homozygotes:
    hapmat[,colnames(hapmat) %in% paste0(substr(colnames(dhom), 1, sform),
                                         substr(colnames(dhom),(sform*2)+1,sform*3))] <- dhom*2
    ##heterozyogote locus 1
    n1 <- paste0(substr(colnames(het_l1), 1, sform),
                 substr(colnames(het_l1),(sform*2)+1,sform*3))
    n1 <- GtoH(het_l1, n1)
    n2 <- paste0(substr(colnames(het_l1),sform+1, sform*2),
                 substr(colnames(het_l1),(sform*3)+1, sform*4))
    n2 <- GtoH(het_l1, n2)
    hapmat[,colnames(hapmat) %in% colnames(n1)] <- n1 + hapmat[,colnames(hapmat) %in% colnames(n1)]
    hapmat[,colnames(hapmat) %in% colnames(n2)] <- n2 + hapmat[,colnames(hapmat) %in% colnames(n2)]

    ##heterozyogote locus 2
    n1 <- paste0(substr(colnames(het_l2), 1, sform),
                 substr(colnames(het_l2),(sform*2)+1,sform*3))
    n1 <- GtoH(het_l2, n1)
    n2 <- paste0(substr(colnames(het_l2),sform+1, sform*2),
                 substr(colnames(het_l2),(sform*3)+1, sform*4))
    n2 <- GtoH(het_l2, n2)
    hapmat[,colnames(hapmat) %in% colnames(n1)] <- n1 + hapmat[,colnames(hapmat) %in% colnames(n1)]
    hapmat[,colnames(hapmat) %in% colnames(n2)] <- n2 + hapmat[,colnames(hapmat) %in% colnames(n2)]


    #5)condense this hap table into the 1a2a, 1a2b, 1b2a, 1b2b format.
    # figure out how where haplotypes are missing. Note, do the case of two or three
    #missin haplotypes at the end.
    pmat <- ifelse(hapmat == 0, F, T)
    mmat <- pmat
    l1 <- substr(colnames(pmat), 1, sform)
    l2 <- substr(colnames(pmat), sform + 1, sform*2)
    mc <- 4 - rowSums(pmat)

    #function to see if haplotype is missing. x is the row index, m is a vector of the number missing haplotypes at each locus.
    cmhap <- function(x){
      out <- ifelse(rowSums(pmat[,l1 == l1[x]]) == 0 | rowSums(pmat[,l1 == l1[x]]) == 2, F,
                    ifelse(rowSums(pmat[,l2 == l2[x]]) > 0 & pmat[,x] != TRUE, T, F))
      return(out)
    }
    #fill the mmat.
    for(i in 1:ncol(pmat)){
      mmat[,i] <- cmhap(i)
    }

    #set the missing values in hapmat to NA, then replace those with zeros where there
    #are missing haplotypes.
    hapmat[hapmat == 0] <- NA
    hapmat[mmat == TRUE] <- 0






    #put in fillers when there are more than one haplotype is missing.
    pmat <- ifelse(is.na(hapmat), F, T)
    missing <- 4 - rowSums(pmat)
    m2 <- ifelse(missing >= 2, 0, NA)
    m3 <- ifelse(missing >= 3, 0, NA)
    m4 <- ifelse(missing == 4, 0, NA)
    mc <- cbind(m2,m2,m3,m4)

    #figure out which D, r values to give if two are missing...
    if(any(missing == 2)){
      m2mat <- hapmat[missing == 2,]
      m2matv <- as.vector(t(m2mat))
      if(is.matrix(m2mat)){
        m2matvcn <- rep(colnames(m2mat), nrow(m2mat))
      }
      else{
        m2matvcn <- rep(names(m2mat), 1)
      }
      m2matv[!is.na(m2matv)] <- m2matvcn[!is.na(m2matv)]
      m2matv <- na.omit(m2matv)
      if(is.matrix(m2mat)){
        m2mat <- matrix(m2matv, nrow(m2mat), 2, T)
        m2mat <- ifelse(substr(m2mat[,1], 1, sform) != substr(m2mat[,2], 1, sform) &
                          substr(m2mat[,1], sform + 1, sform*2) != substr(m2mat[,2], sform + 1, sform*2),
                        1,0)
      }
      else{
        m2mat <- ifelse(substr(m2matv[1], 1, sform) != substr(m2matv[2], 1, sform) &
                          substr(m2matv[1], sform + 1, sform*2) != substr(m2matv[2], sform + 1, sform*2),
                        1,0)
      }

    }
    else{
      m2mat <- "none"
    }

    hapmat <- cbind(hapmat, mc)
    hapmat <- as.vector(t(hapmat))
    hapmat <- na.omit(hapmat)
    #if(length(hapmat) %% nrow(y) != 0){
    #  browser()
    #}
    hapmat <- matrix(hapmat, nrow(ghapmat), 4, byrow = T)

    #now just have the haplotypes. These will calculate D in the case of 1 or 0 missing haplotypes.
    #when there are three missing haplotypes, D will be 0. When there are 2, D will be 0 or 1.
    return(list(hapmat = hapmat, missing = missing, m2 = m2mat))
  }


  #initialize output.
  if(prox_table){
    prox <- data.frame(p1 = numeric(1), p2 = numeric(1), rsq = numeric(1), Dprime = numeric(1), pval = numeric(1))
  }
  if(matrix_out){
    rmat <- matrix(NA, nrow(x) -1, nrow(x) - 1)
    colnames(rmat) <- headers$position[-1]
    rownames(rmat) <- headers$position[-length(headers$position)]
    Dpmat <- rmat
    pvmat <- rmat
  }
  if(!matrix_out & !prox_table){
    stop("Please specify output format.\n")
  }

  #run length prediction variables for progress reporting
  compfun <- function(x){
    return(((x-1)*x)/2)
  }
  totcomp <- compfun(nrow(x))
  cpercent <- 0

  #loop through and get haplotypes, calc LD for each locus.
  for(i in 1:(nrow(x) - 1)){
    cprog <- (totcomp - compfun(nrow(x) - i - 1))/totcomp
    if(cprog >= 0.05 + cpercent){
      cat("Progress:", paste0(round(cprog*100), "%."), "\n")
      cpercent <- cprog
    }
    haps <- tabulate_haplotypes(x[i,], x[(i+1):nrow(x),])

    #if we had only one haplotype or no haplotypes:
    if(is.na(haps[1])){
      if(prox_table){
        prox <- rbind(prox,
                      cbind(p1 = headers$position[i], p2 = headers$position[(i+1):nrow(x)],
                            rsq = NA, Dprime = NA, pval = NA))
      }
      if(matrix_out){
        #reminder: columns start at locus two, rows start at locus one (but end at nlocus - 1)
        rmat[i,] <- NA
        Dpmat[i,] <- NA
        pvmat[i,] <- NA
      }
      next()
    }
    if(length(haps) == 1){
      if(prox_table){
        prox <- rbind(prox,
                      cbind(p1 = headers$position[i], p2 = headers$position[(i+1):nrow(x)],
                            rsq = 0, Dprime = 0, pval = 0))
      }
      if(matrix_out){
        #reminder: columns start at locus two, rows start at locus one (but end at nlocus - 1)
        fill <- rep(NA, nrow(x) - (nrow(x[,i:nrow(x)]) - 1) - 1)
        Dprime <- c(fill, rep(0, length = (nrow(x[,i:nrow(x)]) - 1)))
        rsq <- c(fill, rep(0, length = (nrow(x[,i:nrow(x)]) - 1)))
        pval <- c(fill, rep(0, length = (nrow(x[,i:nrow(x)]) - 1)))
        rmat[i,] <- rsq
        Dpmat[i,] <- Dprime
        pvmat[i,] <- pval
      }
      next()
    }


    missing <- haps$missing
    m2 <- haps$m2
    haps <- haps$hapmat
    #A1B1 is col 1, A1B2 is col 2, A2B1 is col 3, A2B2 is col 4.

    #calc stats where >1 haps aren't missing
    A1B1f <- haps[,1]/rowSums(haps)
    A1B2f <- haps[,2]/rowSums(haps)
    A2B1f <- haps[,3]/rowSums(haps)
    A2B2f <- haps[,4]/rowSums(haps)
    A1f <- A1B1f + A1B2f
    A2f <- 1 - A1f
    B1f <- A1B1f + A2B1f
    B2f <- 1 - B1f
    D <- A1B1f - A1f*B1f
    D2 <- A2B1f - A2f*B1f
    Dprime <- ifelse(D > 0, D/rowMins(cbind(A1f*B2f, A2f*B1f)),
                     ifelse(D < 0, D/rowMaxs(cbind(-A1f*B1f, -A2f*B2f)),
                            0))
    rsq <- (D^2)/(A1f*A2f*B1f*B2f)

    #fix for when more missing haps.
    Dprime[missing == 3] <- 0
    Dprime[missing == 4] <- NA
    rsq[missing == 3 | missing == 4] <- NA
    if(length(m2) > 1){
      Dprime[missing == 2] <- m2
      rsq[missing == 2] <- ifelse(m2 == 1, 1, NA)
    }

    #get pvals
    chisqu <- rsq*(4)
    pval <- 1 - pchisq(chisqu, 1)

    #remove dummy filler if this was the final comparison.
    if(i == (nrow(x) - 1)){
      Dprime <- Dprime[-2]
      rsq <- rsq[-2]
      pval <- pval[-2]
    }

    #write output.
    if(prox_table){
      prox <- rbind(prox,
                    cbind(p1 = headers$position[i], p2 = headers$position[(i+1):nrow(x)],
                          rsq = rsq, Dprime = Dprime, pval = pval))
    }
    if(matrix_out){
      #reminder: columns start at locus two, rows start at locus one (but end at nlocus - 1)
      fill <- rep(NA, nrow(x) - length(Dprime) - 1)
      Dprime <- c(fill, Dprime)
      rsq <- c(fill, rsq)
      pval <- c(fill, pval)
      rmat[i,] <- rsq
      Dpmat[i,] <- Dprime
      pvmat[i,] <- pval
    }
  }
  if(prox_table){
    prox <- prox[-1,]
    if(any(colnames(headers) == "group")){
      prox$group <- headers[1,"group"]
    }
    if(any(colnames(headers) == "pop")){
      prox$pop <- headers[1,"pop"]
    }
    prox <- prox[,c(6:ncol(prox), 1:5)]
  }
  return(list(prox = prox, Dprime = Dpmat, rsq = rmat, pval = pvmat))
}


#Calls the LD_full_pairwise function for a file with multiple linkage groups. Returns a nested list of
#Dprime and rsq matrices. Requires the linkage group titled be named "group". Returns a concatenated prox_table.
Full_LD_w_groups <- function(x, ecs, prox_table = TRUE, matrix_out = TRUE, mDat = "N", report = 1){
  if (prox_table == FALSE & matrix_out == FALSE){
    stop("Please specify output format.")
  }
  w_list<- list()
  prox_out <- data.frame()
  groups <- unique(x$group)
  for (i in 1:length(groups)){
    w_data <- subset(x, group == groups[i])
    if(nrow(w_data) == 0 | nrow(w_data) == 1){
      next()
    }
    if(i %% report == 0){cat("Group #:", i, "of", length(groups), "Name:", groups[i], "\n")}
    w_data <- LD_full_pairwise(x = w_data, ecs = ecs, prox_table = prox_table, mDat = mDat, matrix_out = matrix_out)
    if(prox_table == TRUE){
      if(matrix_out == TRUE){
        prox_out <- rbind(prox_out, w_data$prox)
        w_list[[groups[i]]] <- w_data[1:2]
      }
      else{
        prox_out <- rbind(prox_out, w_data)
      }
    }
    else if (matrix_out == TRUE){
      w_list[[groups[i]]] <- w_data
    }
  }
  if (prox_table == TRUE){
    if(matrix_out == TRUE){
      w_list[["prox_table"]] <- prox_out
    }
    else{
      w_list <- prox_out
    }
  }
  return(w_list)
}


#runs LD_full_pairwise in parallel across many different linkage groups/scaffolds. Arguments are the same
#as LD_full_pairwise, save for num_cores (which is the number of cores to run on).
#If all outputs are desired, Returns a list where the first length(groups) elements are lists containing matrices for rsq and Dprime for
#the corresponding group. The last element is a prox table.
#currently giving an error, need to check it out in more detail.
Full_LD_g_par <- function(x, ecs, num_cores, prox_table = TRUE, matrix_out = TRUE, mDat = "N"){
  if (prox_table == FALSE & matrix_out == FALSE){
    stop("Please specify output format.")
  }
  cl <- makeSOCKcluster(num_cores)
  registerDoSNOW(cl)


  groups <- unique(x$group)
  ntasks <- length(groups)
  progress <- function(n) cat(sprintf("Group %d out of",n), length(groups), "is complete.\n")
  opts <- list(progress=progress)

  output <- foreach(i = 1:length(groups), .export = 'LD_full_pairwise', .packages = c("dplyr", "reshape2", "matrixStats"), .inorder = TRUE,
                    .options.snow = opts) %dopar% {
                      w_data <- x[x$group == groups[i],]
                      cat("Initializing :", groups[i], ".\n")
                      if(nrow(w_data) == 0 | nrow(w_data) == 1){
                        next()
                      }
                      w_data <- LD_full_pairwise(w_data, ecs = ecs, prox_table = prox_table, matrix_out = matrix_out, mDat = mDat)
                    }
  #loop through output and put in easier to manage form
  if(prox_table){p <- c()}
  if(matrix_out){m <- vector("list", length(groups))}
  for(i in 1:length(groups)){
    if(prox_table){
      if(matrix_out){
        p <- rbind(p, output[[i]]$prox)
      }
      else{
        p <- rbind(p, output[[i]])
      }
    }
    if(matrix_out){
      m[[i]] <- output[[i]][-1]
      names(m) <- groups
    }
  }

  #clear cluster and return stuff
  stopCluster(cl)
  registerDoSEQ()
  if(prox_table){
    if(matrix_out){
      m[[i+1]] <- p
      names(m) <- c(groups, "prox")
      return(m)
    }
    else{
      return(p)
    }
  }
  else{
    return(m)
  }
}
