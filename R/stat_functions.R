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
#' @param x Input SNP data, in allele count format as given by format_snps output option.
#' @param ecs Numeric, number of metadata columns at the start of x to keep for output.
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
calc_pi <- function(x, facets = NULL){
  pi_func <- function(x){
    nt <- as.numeric(x[,"n_total"])
    n1 <- as.numeric(x[,"ni1"])
    n2 <- as.numeric(x[,"ni2"])
    binom_n1 <- choose(n1,2)
    binom_n2 <- choose(n2,2)
    binom_nt <- choose(nt,2)
    #print(n1)
    #print(n2)
    p <- 1 - (binom_n1 + binom_n2)/binom_nt
  }

  out <- apply.snpR.facets(x, facets, "ac", pi_func, case = "ps")
  colnames(out)[ncol(out)] <- "pi"
  x <- merge.snpR.stats(x, out)

  return(return(x))
}


#function to calc pi from data. should be in format with num alleles, total count of alleles,
#and subsequent alleles counts in columns named "n_aleles", "n_total", and "ni1", "ni2", and so on (allele count/bayescan format, as given by format_snps option one).
#Returns a VECTOR of pi values.

#'Find minor and major allele frequencies.
#'
#'\code{calc_maf} determines the minor and major allele frequencies for snp data for each locus, splitting by any desired facets.
#'
#' @param x Input SNP data, a snpRdata object.
#' @param facets Character vector, default NULL. Name of a given facet to calculate maf for. If NULL calculates for all facets in x.
#'
#' @return Adds minor and major allele frequency information to a snpRdata object, accessable under the "\\@min.maj".
calc_maf <- function(x, facets = NULL){
  # function to run on whatever desired facets:
  func <- function(gs, m.al){
    as <- gs$as
    majl <- rep(matrixStats::rowMaxs(as), each = ncol(as))
    majl <- as == matrix(majl, ncol = ncol(as), byrow = T)

    # set aside special cases
    np <- which(rowSums(matrix(as.logical(as), nrow(as))) == 1) # non-polymorphic
    eqf <- which(rowSums(majl) == 2) # two alleles of equal frequency
    ns <- which(rowSums(as) == 0) # nothing sequenced
    vio <- sort(c(np, eqf, ns)) # all
    if(length(vio) != 0){
      majl.s <- majl[vio,]
      as.s <- as[vio,]
      majl <- majl[-vio,]
      as <- as[-vio,]
    }

    # major and minor alleles for everything without issues
    min <- which(t(as.logical(as) + majl) == 1)
    min <- min %% ncol(majl)
    min[min == 0] <- 4
    maj <- which(t(majl)) %% ncol(majl)
    maj[maj == 0] <- 4
    maj <- colnames(as)[maj] # these are the major alleles
    min <- colnames(as)[min] # these are the minor alleles

    # fix for special cases
    if(length(vio) != 0){
      if(is.null(nrow(majl.s))){
        majl.s <- matrix(majl.s, nrow = 1)
        colnames(majl.s) <- colnames(gs$as)
        as.s <- matrix(as.s, nrow = 1)
        colnames(as.s) <- colnames(gs$as)
      }
      maj.s <- character(nrow(majl.s))
      min.s <- maj.s

      # non-poly
      if(length(np) != 0){
        wnp <- which(rowSums(matrix(as.logical(as.s), nrow(as.s))) == 1)
        np.as <- as.s[wnp,]
        np.maj <- which(as.logical(t(np.as))) %% ncol(majl)
        np.maj[np.maj == 0] <- 4
        np.maj <- colnames(as)[np.maj] # these are the major alleles
        maj.s[wnp] <- np.maj
        min.s[wnp] <- "N"
      }

      # equal frequencies
      if(length(eqf) != 0){
        weqf <- which(rowSums(majl.s) == 2)

        eq.min <- which(t(majl.s[weqf,])) %% ncol(majl.s)
        eq.maj <- eq.min[seq(1, length(eq.min), by = 2)]
        eq.min <- eq.min[seq(2, length(eq.min), by = 2)]

        eq.min[eq.min == 0] <- 4
        eq.maj[eq.maj == 0] <- 4
        maj.s[weqf] <- colnames(as.s)[eq.maj] # these are the major alleles
        min.s[weqf] <- colnames(as.s)[eq.min] # these are the minor alleles
      }

      # not sequenced
      if(length(ns) != 0){
        wns <- which(rowSums(as.s) == 0)
        maj.s[wns] <- m.al
        min.s[wns] <- m.al
      }

      # add the major and minors for the special cases back in by binding and ordering
      f.order <- order(c((1:nrow(gs$as))[-vio], vio)) # this gets the indices of the non-special cases, binds them to the indices of the special cases, and then gets the order to put them in.
      maj <- c(maj, maj.s)
      maj <- maj[f.order]
      min <- c(min, min.s)
      min <- min[f.order]
    }

    # grab the actual maf
    maf <- 1 - matrixStats::rowMaxs(gs$as)/rowSums(gs$as)
    maf[is.nan(maf)] <- 0

    # get the major and minor counts
    # round because repeating decimals will yeild things like 1.00000000001 instead of 1. Otherwise this approach is quick and easy, as long as things are bi-allelic (non-polymorphic and equal min maj frequencies are fine.)
    maj.count <- round(rowSums(gs$as)*(1-maf))
    min.count <- round(rowSums(gs$as)*(maf))

    # return
    return(data.frame(major = maj, minor = min, maj.count = maj.count, min.count = min.count, maf = maf, stringsAsFactors = F))
  }

  # add any missing facets
  if(!any(facets == "all")){
    x <- add.facets.snpR.data(x, facets)
  }

  out <- apply.snpR.facets(x,
                           facets = facets,
                           req = "gs",
                           fun = func,
                           case = "ps",
                           m.al = substr(temp@mDat,1, nchar(temp@mDat)/2))
  return(merge.snpR.stats(x, out))
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
#'\code{Tajimas_D} calculates Tajima's theta/pi, Waterson's theta, and Tajima's D over a sliding window. Pi calculated as in Hohenlohe et al. 2010. Tajima's D is calculated using the formula from Tajima (1989) Statistical Method for Testing the Neutral Mutation Hypothesis by DNA Polymorphism.
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
#' ac <- format_snps(stickSNPs, 3, pop = pops)
#' run_gp(ac, Tajimas_D, 200, 50)
#'
Tajimas_D <- function(x, ws, step, report = 20){

  #windowing
  out <- data.frame() #initialize output
  tps <- sort(x$position) #get the site positions, sort
  lsp <- tps[length(tps)] #get the position of the last site to use as endpoint
  c <- 0 #set starting position
  i <- 1 #set starting iteration for writing output
  ws <- 1000*ws
  while (c <= lsp){
    start <- c - ws #change window start
    end <- c + ws #change window end
    if(i %% report == 0){cat("Window Postion:", c, "out of", lsp, "\n")}

    #take all the snps in the window, calculate T's theta, W's theta, and T's D
    wsnps <- x[x$position <= end & x$position >= start,] #get only the sites in the window
    wsnps <- wsnps[!(wsnps[,"ni1"] == 0 & wsnps[,"ni2"] == 0),] #remove any sites that are not sequenced in this pop/group/whatever
    if(nrow(wsnps) == 0){ #if no snps in window
      c <- c + step*1000 #step along
      next #go to the next window
    }
    b1s <- choose(wsnps[,"ni1"],2) #binomial for first allele
    b2s <- choose(wsnps[,"ni2"],2) #binomial for second allele
    bts <- choose(wsnps[,"n_total"],2) #binomial for all alleles
    ts.theta <- sum(1-((b1s+b2s)/bts)) #average number of pairwise differences (pi) per snp. Equivalent to sum((ndif/ncomp)) for all snps
    #ts.thetaf <- ts.theta/nrow(wsnps) #pi for the whole region, which includes all of the non-polymorphic sites. Reason why this shouldn't be run with a maf filter, probably.
    n_seg <- nrow(wsnps[wsnps$ni1 != 0 & wsnps$ni2 != 0,]) #number of segregating sites
    K <- round(mean(wsnps[,"n_total"])) #average sample size for ws.theta, as in hohenlohe 2010. Alternative would make this into a function, then use lapply on the vector of K values
    #if(is.nan(K) == TRUE){browser()}
    a1 <- sum(1/seq(1:(K-1))) #get a1
    ws.theta <- n_seg/a1 #get ws.theta
    #ws.thetaf <- ws.theta/nrow(wsnps) #ws.theta fraction

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

    out[i,"group"] <- x[1,"group"]
    if("pop" %in% colnames(x)){
      out[i,"pop"] = x[1,"pop"] #if a pop column is in the input, add a pop column here.
    }
    out[i,"position"] <- c
    out[i,"ws.theta"] <- ws.theta
    out[i,"ts.theta"] <- ts.theta
    out[i,"D"] <- D
    c <- c + step*1000
    i <- i + 1
  }
  return(out)
}


#Calculates pairwise fst for each pair of populations according to hohenlohe (2010). Input format is
#snp, position, group, pop, total allele count (total), number of alleles (num), ni1, ni2, and pi (diversity, from calc_pi).
#Automatically sorts data by group, position, and population.

#'Pairwise FST from SNP data.
#'
#'\code{calc_pairwise_Fst} calculates pairwise FST for each SNP for each possible pairwise combination of populations based on the methods in Wier and Cockerham 1984, Wier 1990, or Hohenlohe et al 2010. Can also return the total number of observed alleles at each SNP for each of these pairwise combinations.
#'
#'Description of x:
#'    Must contain colums containing the number of *unique* alleles, total count of alleles sequenced in all individuals, and subsequent alleles counts for each observed allele in columns named "n_alleles", "n_total", "ni1", and "ni2". This matches the allele count/bayescan format as given by format_snps option one. Should also contain columns titled "group", "position", and "pop", which contain the linkage group/chr, position in bp, and population ID for each SNP. Lastly, also needs a column containing pi (for Hohenlohe) or Ho (for WC or Wier) for each SNP, as given by calc_pi or calc_Ho, respectively.
#'
#'Smoothing note:
#'    To smooth the data resulting from this, simply convert to long format, name the column containing comparisons "pop", append a column with long format pairwise nk values if desired, and use run_gp. Conversion for the FST and nk data frames can be easily done with rehape2::melt, with all metadata columns given as id.vars, as in id.vars = c("snp", "group", "pop").#'
#'
#'Method Options:
#'\itemize{
#'    \item{"WC": }{Wier and Cockerham 1984.}
#'    \item{"Wier": }{Wier 1990.}
#'    \item{"Genepop": }{As used in genepop, Rousset 2008.}
#'    \item{"Hohenlohe": }{Hohenlohe 2010.}
#'}
#'
#' @param x Input data frame, in allele count format as given by format_snps option "ac", with an additional column containing pi or Ho estimates.
#' @param ecs Number of extra metadata columns at the start of x *not counting the column with population IDs* but counting position, ect. as normal.
#' @param do.nk Should pairwise nk (allele sample sizes) also be calculated?
#' @param skip.FST Should FST calcs be skipped (only set to TRUE if do.nk is as well).
#' @param method Which FST estimator should be used?
#' @param pnames Character vector, default NULL. Vector of population names to use for method "Genepop"
#' @param char.dat Data.frame, defualt NULL. Input data in character format, as given by format_snps option "character". Only needed if no Ho column provided to Wier or WC methods.
#' @param m.dat Character, default "NN". Missing *genotype* identifier in char.dat, used if no Ho column provided to Wier or WC methods.
#' @param pop FALSE or table, default FALSE. A table with population information for individuals. Individuals must be sorted in input data in the population order given in this table. Only needed if WC or Wier used but no "Ho" column provided.
#' @param c.d.ecs Numeric, default NULL. Number of extra metadata columns at the start of char.dat. Only needed if WC or Wier method chosen but no Ho column provided.
#'
#' @return If both do.nk is true and skip.FST is false, returns a list containing named data frames "FST" and "nk", which contain pairwise FST and nk values, respectively If only one output is requested, only that data frame is returned. If "Genepop" is chosen, returns a list with both a data.frame containing pairwise FST values and a vector of the overall weighted FST for each comparison.
#' @examples
#' #Wier and Cockerham
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' ac <- format_snps(stickSNPs, 3, pop = l)
#' Ho <- calc_Ho(stickSNPs, 3, pop = l)
#' m_Ho <- melt(Ho, id.vars = c("group", "position", "snp"))
#' ac$Ho <- m_Ho$value
#' calc_pairwise_Fst(ac, 3, do.nk = T)
#'
#' #Wier and Cockerham without binding the Ho column manually:
#' calc_pairwise_Fst(ac, 3, do.nk = T, char.dat = stickSNPs, pop = pops, c.d.ecs = 3)
#'
#'
#' #genepop
#' ##prepare pop list for both formatting and FST.
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' l <- list(c(names(pops)), as.numeric(pops))
#' ##write the genepop file.
#' format_snps(stickSNPs, 3, 2, outfile = "stickGENEP.txt", pop = l)
#' ##calculate FST:
#' FST <- calc_pairwise_Fst("stickGENEP.txt", method = "Genepop", pnames = l[[1]])

calc_pairwise_Fst <- function(x, facets, do.nk = FALSE, skip.FST = FALSE, method = "WC"){


  #============================sanity and facet checks========================
  browser()
  if(!do.nk & skip.FST){
    stop("Must specify either pairwise FST, nk, or both")
  }
  if(any(x@ac$n_alleles > 2)){
    vio <- which(x@ac$n_alleles > 2)
    vio <- unique(x@facet.meta$.snp.id[vio])
    stop(cat("Some loci have more than two alleles. Violating loci:\n", paste0(vio, collapse = "\n")))
  }

  add.facets <- which(!(facets %in% x@facets))
  if(length(add.facets) != 0){
    cat("Adding missing facets.\n")
    x <- add.facets.snpR.data(x, facets[add.facets])
  }
  # work here






  ###############################################################################
  #do genepop if requested
  if(method == "Genepop"){
    #sanity checks
    if(!is.character(x) | length(x) != 1){
      stop("For genepop, x is the path to the genepop file. Format_snps output option 2 can provide this.\n")
    }
    if(!file.exists(x)){
      stop("For genepop, x is the path to the genepop file. Format_snps output option 2 can provide this.\n")
    }
    if(is.null(pnames)){
      cat("No pop names provided, pulling from file.\n")
    }
    else if(!is.character(pnames)){
      stop("pnames must be a character vector.\n")
    }

    if(do.nk){
      stop("Set method to something other than Genepop for do.nk.\n")
    }
    if(skip.FST){
      stop("Nothing to do. Stopping.\n")
    }

    #run genepop
    cat("Genepop results:\n\t")
    genepop::Fst(x, pairs = TRUE)

    ############
    #read in the genepop output and parse the bugger.
    cat("Parsing genepop output...\n")

    #read the file in
    x <- paste0(x, ".ST2") #data file
    x <- readLines(x)

    #get the number of pops and the number of loci.
    np <- grep("Number of populations detected", x)
    np <- as.numeric(unlist(strsplit(x[np], " : "))[2])
    nl <- grep("Number of loci detected", x)
    nl <- as.numeric(unlist(strsplit(x[nl], " : "))[2])

    #check that the correct number of pop names were provided or grab new ones if they weren't.
    ##get pop data
    px <- grep("Indices for populations:", x)
    px <- x[(px+2):(px+1+np)]
    px <- unlist(strsplit(px, " +"))
    px <- px[seq(2, length(px), 2)]

    if(is.null(pnames)){
      cat("Pulling pop names from x...\n")
      pnames <- px
    }
    else if (length(pnames) != np){
      warning("Vector of provided pop names is not equal in length to the number of populations in x!\n")
      cat("Substituting pop names form x...\n")
      pnames <- px
    }
    else{
      if(all(pnames != sort(pnames))){
        warning("Pop names have been sorted, ensure input genepop file also has alphabetically sorted pop names.")
      }
      pnames <- sort(pnames)
    }


    #get the indices containing locus headers
    locs <- grep("  Locus:", x)
    locs <- c(locs, grep("Estimates for all loci", x))

    #get indices not containing data to parse and remove them.
    empts <- c(1:(locs[1]-2), locs, locs + 1, locs + 2, locs - 1, (length(x)-2):length(x))
    vals <- x[-empts]

    #initialize output.
    fmat <- matrix(NA, nrow = nl + 1, ncol = ((np-1)*np)/2)
    colnames(fmat) <- 1:ncol(fmat)

    #fill the matrix with a loop. Not a bad loop, since it only loops through each pop.
    prog <- 1 #column fill progress tracker
    for(i in 2:np){
      #grab the values
      tvals <- vals[grep(paste0("^", i, " +"), vals)] #get just the comparisons with this pop
      tvals <- gsub(paste0("^", i, " +"), "", tvals) #get just the values
      tvals <- unlist(strsplit(tvals, " +")) #spit and unlist the values
      tvals <- suppressWarnings(as.numeric(tvals)) #get them as numeric, NAs are fine, they should be NAs.

      #put them in a matrix to get their comparison ID.
      tmat <- matrix(tvals, nl + 1, i - 1, TRUE) #fill these values into a matrix
      colnames(tmat) <- paste0(pnames[1:(i-1)], "_", pnames[i]) #name the comparison in each column
      fmat[,prog:(prog + ncol(tmat) - 1)] <- tmat #save to the final output
      colnames(fmat)[prog:(prog + ncol(tmat) - 1)] <- colnames(tmat) #save the column names
      prog <- prog + ncol(tmat) #increment column progress.
    }

    fmat <- fmat[,order(colnames(fmat))] #re-organize output by column name

    if(np != 2){
      #grab out the overall Fst
      overall <- fmat[nrow(fmat),]

      #grab the per-locus estimates
      fmat <- fmat[-nrow(fmat),]
      out <- list(loci = as.data.frame(fmat, stringsAsFactors = F), overall = overall)
    }
    else{
      overall <- fmat[length(fmat)]
      fmat <- fmat[-length(fmat)]
      out <- list(loci = fmat, overall = overall)
    }
    #return, we're done.
    cat("Finished.\n")
    return(out)
  }

  ###############################################################################
  #otherwise do in house Fst calcs
  pops <- sort(unique(x[,"pop"]))
  out <- as.data.frame(matrix(NA, ncol = ecs+(length(pops)*(length(pops) - 1)/2), nrow = nrow(x)/length(pops)))

  if(method == "WC" | method == "Wier"){
    # if Ho not provided, need to calculate it and merge it into the dataset.
    if(!"Ho" %in% colnames(x)){
      ho <- calc_Ho(char.dat, ecs = c.d.ecs, mDat = mDat, pop = pop)
      ho <- reshape2::melt(ho, id.vars = colnames(x)[1:ecs])
      colnames(ho)[(ecs + 1):ncol(ho)] <- c("pop", "Ho")
      x <- merge(x, ho, by.x = c(colnames(x)[1:ecs], "pop"), by.y = colnames(ho)[-ncol(ho)], sort = F)
    }
  }
  if(method == "Hohenlohe"){
    ac$pi <- calc_pi(ac) # add a pi column.
  }


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
  colnames(out) <- c(colnames(x)[1:ecs], comps)
  out[,1:ecs] <- x[x$pop == pops[1],1:ecs] #add snp, position, group to output
  if(do.nk){ #prepare nk output if requested
    nout <- out
    colnames(nout)[(ecs+1):ncol(nout)] <- paste0("nk_", comps)
  }

  #print(out)
  #loop through each comparison and caculate pairwise FST at each site
  c.col <- ecs + 1 #set starting column for pasting data
  for (i in 1:(length(pops) - 1)){ #i is the first pop
    idat <- x[x$pop == pops[i],] #get data for first pop
    j = i + 1 #intialize j as the next pop
    for (j in j:length(pops)){#j is pop being compared
      jdat <- x[x$pop == pops[j],] #get data for second pop
      if(!skip.FST){
        if(method == "WC" | method == "Wier"){

          #allele frequencies in both i and jdat.
          ps1 <- idat$ni1/idat$n_total
          ps2 <- jdat$ni1/jdat$n_total

          #other stats
          r <- 2 #number of comps
          nbar <- (idat$n_total + jdat$n_total)/2 #average sample size
          comb_ntots <- cbind(idat$n_total, jdat$n_total)
          CV <- matrixStats::rowSds(comb_ntots)/rowMeans(comb_ntots) # coefficient of variation in sample size
          nc <- nbar*(1-(CV^2)/r)
          pbar <- ((idat$n_total*ps1) + (jdat$n_total*ps2))/(r*nbar) #average sample allele frequency
          ssq <- (((idat$n_total)*(ps1-pbar)^2) + ((jdat$n_total)*(ps2-pbar)^2))/((r-1)*nbar) #sample varaince of allele frequencies
          hbar <- ((idat$n_total*idat$Ho) + (jdat$n_total*jdat$Ho))/(r*nbar) #average heterozygote frequencies

          #equation parts used in both
          inner1 <- pbar*(1-pbar)
          inner2 <- ((r-1)/r)*ssq
          inner3 <- .25*hbar

          if(method == "WC"){
            inner4 <- ((2*nbar - 1)/(4*nbar))*hbar
            a <- (nbar/nc) * (ssq - (1/(nbar - 1))*(inner1 - inner2 - inner3))
            b <- (nbar/(nbar-1))*(inner1 - inner2 - inner4)
            c <- .5*hbar
            Fst <- a/(a+b+c)
            out[,c.col] <- Fst #write fst
          }
          else{
            S1 <- ssq - (1/(nbar-1))*(inner1 - inner2 - inner3)
            S2i1 <- ((r*(nbar - nc))/nbar)*inner1
            S2i2 <- (1/nbar)*((nbar-1)+(r-1)*(nbar-nc))*ssq
            S2i3 <- ((nbar-nc)/(4*nc^2))*hbar
            S2 <- inner1 - (nbar/(r*(nbar-1)))*(S2i1 -S2i2 - S2i3)
            Fst <- S1/S2
            out[,c.col] <- Fst #write fst
          }

        }

        else if(method == "Hohenlohe"){

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
          if(any(abs(Fst) > 1, na.rm = T)){cat("Fst > 1 at", which(Fst > 1), ". That's not good.");stop()}
          #Fst[t.ni1 == 0 | t.ni2 == 0] <- 0 #could uncomment this if want 0s instead of NaNs.
          out[,c.col] <- Fst #write fst
        }

        else{
          stop("Please select a method of calculating FST.\nOptions:\n\tWC: Weir and Cockerham (1984).\n\tWier: Wier (1990).\n\tHohenlohe: Hohenlohe et al. (2010).")
        }
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

#'Pairwise nk from SNP data.
#'
#'\code{calc_pairwise_nk} calculates number of alleles observed in each parwise combination of populations. Wrapper function for \code{\link{calc_pairwise_Fst}}
#'
#' Description of x:
#'     Must contain colums containing the number of *unique* alleles, total count of alleles sequenced in all individuals, and subsequent alleles counts for each observed allele in columns named "n_alleles", "n_total", "ni1", and "ni2". Also needs a column containing the position of each SNP, in bp. This matches the allele count/bayescan format as given by format_snps option one. Should also contain columns titled "group", "position", and "pop", which contain the linkage group/chr, position in bp, and population ID for each SNP.
#'
#' @param x Input data frame, in allele count format as given by format_snps option 1.
#' @param ecs Number of extra metadata columns at the start of x *not counting the column with population IDs* but counting position, ect. as normal.
#'
#' @return A data.frame containing the number of alleles secquenced in each pair of populations.
#'
#' @examples
#' calc_pairwise_nk(ac, 3)
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

#'Observed heterozygosity from SNP data.
#'
#'\code{calc_Ho} calculates observed heterozygosity at each SNP, potentially in many populations.
#'
#'Description of x:
#'    Contains metadata in columns 1:ecs. Remainder of columns contain genotype calls for each individual. Each row is a different SNP, as given by format_snps output options 4 or 6. If the data contains individuals from multiple populatoins, they must be in the same order as given to the pop agrument.
#'
#' @param x Input data, in "NN" or "0000" format, as given by format_snps output option 4 or 6.
#' @param ecs Number of extra metadata columns at the start of x.
#' @param m.dat Character variable matching the coding for missing *genotypes* in x (typically "NN" or "0000").
#' @param pop FALSE or table, default FALSE. A table with population information for individuals. Individuals must be sorted in input data in the population order given in this table.
#' @return Data frame containing metadata, and Ho for each population in named columns.
#'
#' @examples
#' #no seperate pops:
#' calc_Ho(stickSNPs, 3)
#'
#' #seperate pops
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' calc_Ho(stickSNPs, 3, pop = pops)
#'
calc_Ho <- function(x, ecs, mDat = "NN", pop = NULL){

  #set possible heterozygotes
  if(nchar(x[1,(ecs + 1)]) == 4 & nchar(mDat) == 4){
    hl <- c(mDat, "0101", "0202", "0303", "0404")
  }
  else if (nchar(x[1,(ecs + 1)]) == 2 & nchar(mDat) == 2){
    hl <- c(mDat, "AA", "TT", "CC", "GG")
  }
  else{
    stop("Data and missing signifier must be in either one (NN) or two character (0000) per allele format.")
  }

  #initalize output
  if(!is.table(pop)){
    pop <- table(rep("ho", ncol(x) - (ecs + 1) + 1))
  }
  pns <- names(pop)
  psz <- as.numeric(pop)
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
    ho <- het.c/rowSums(ifelse(wdat == mDat, 0, 1))
    out[,(ecs + 1) + j - 1] <- ho
    c.col <- c.col + psz[j]
  }
  return(out)
}

#Checks for private alleles.
#inputs: x: data, in allele count format such as that given by format_snps option 1. Expects
#           columns named "pop", "ni1", and "ni2", which contain pop designations, allele one counts,
#           and allele 2 counts.

#'Private Alleles from SNP data.
#'
#'\code{check_private} checks for the presence of private alleles at any loci across multiple populations.
#'
#'Description of x:
#'    Must contain colums containing columns containing the allele counts for each observed allele in columns named "ni1"" and "ni2". Also needs a column containing population IDs, in a column named "pop". Note that SNPs must be identically sorted in each pop, and sorted first by pop! For example, sorted by pop, group, position then position. Runs with output from format_snps output option 1.
#'
#' @param x Input data, in the format given by format_snps output option "ac". Must contain a column with pop info named "pop".
#' @param ecs Numeric, number of metadata columns at the start of x to maintain for output.
#'
#' @return Data frame with no metadata but sorted identically to input noting if each (row) locus is private in each (column) population, coded as 1 if private, 0 if not.
#'
#' @examples
#' #add a private allele for demonstration purposes.
#' pall <- data.frame(snp = rep(10000, 6), position = rep(1, 6), group = rep("Dup", 6), pop = c("ASP", "CLF", "PAL", "OPL", "SMR", "UPD"), n_total = rep(100, 6), n_alleles = rep (2, 6), ni1 = c(rep(0, 5), 20), ni2 = c(rep(100, 5), 80))
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' ac <- format_snps(stickSNPs, 3, pop = pops)
#' ac <- rbind(ac, pall)
#' ac <- dplyr::arrange(ac, pop, position, group)
#' check_private(ac)
#'
check_private <- function(x, ecs){
  l <- unique(x$pop) #gets unique pops, requires column named pop

  if(all(c("group", "position", "pop") %in% colnames(x))){
    x <- dplyr::arrange(x, pop, group, position)
  }
  else{warning("Data must be sorted by pop, then identicallly by other meta data, such as by group and position. If columns named pop, group, and position are given in x, this will be done automatically.")}

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
  colnames(pa) <- paste0(l)
  return(cbind(x[,1:ecs], as.data.frame(pa)))
}

#Calculates Dprime, rsq, and a p-value for LD for each pair of snps.
#inputs: x: data, in either numeric or NN form. Must contain a column named "position", can contain a column named "group" and/or "pop"
#        ecs: number of extra columns (headers) before the data in x
#        prox_table: Should a proximity table be output?
#        matrix_out: Should LD matrices be created?
#        mDat: What is the missing data character? Expects marker for a single allele. ("N" or "01")

#'Pairwise LD from SNP data.
#'
#'\code{LD_full_pairwise} calculates LD between each pair of SNPs. If called as is, assumes one linkage group/chromosome and one population!
#'
#'Matrix outputs are properly formated for LD_pairwise_heatmap.
#'
#'Description of x:
#'    Contains metadata in columns 1:ecs. Remainder of columns contain genotype calls for each individual. Each row is a different SNP, as given by format_snps output options 4 or 6. Requires the column containing the position of the loci in base pairs be named "position". Note that this \emph{ignores populations}, split data into distinct populations before running.
#'
#' @param x Input data, usually in either numeric or character formats for SNP data, as given by format_snps output options "numeric" or "character". ms filepath or phased haplotypes (such as from an imported ms output) also accepted.
#' @param ecs Number of extra metadata columns at the start of x.
#' @param prox_table If TRUE, a proximity table is produced.
#' @param matrix_out If TRUE, pairwise LD matrices are produced.
#' @param mDat Character variable matching the coding for missing *genotypes* in x (typically "NN" or "0000").
#' @param sr boolean, default FALSE. Should progress reports be surpessed?
#' @param input Character, default "SNP". Input data type. Options: "SNP", as given by format_snps output options "numeric" or "character". "ms": filepath to ms formatted data. "haplotype": imported ms data, where each column is a fully phased gene copy/chromosome.
#' @param chr.length Numeric, default NULL. Length of chromosomes, for ms inputs.
#' @param levels Levels to split the data by for LD calculations, such as chromosome or linkage group. Must match in input column name.
#' @param par Numeric or FALSE, default FALSE. Number of cores to parallelize the analysis by.
#' @param ss Numeric or FALSE, default FALSE. Number or proportion of SNPs to subsample for LD analysis.
#' @param level_report Numeric, default 1. Progress through levels is reported every level_report levels.
#' @param maf Numeric or FALSE, default FALSE. Minor allele frequency filter cuttoff, only for ms or haplotype data.
#'
#' @return Matrices containing rsq, Dprime, and p values for each SNP vs every other SNP. Can also produce a proximity table, which contains the rsq, Dprime, and p value for each pairwise comparison and the distance between the SNPs in those comparisons. Returns matrices and prox table as sequential elements in a named list.
#'
#' @examples
#' #returns prox table and LD matrices.
#' LD_full_pairwise(stickSNPs[stickSNPs$group == "groupI",1:53], ecs = 3)
#'
LD_full_pairwise <- function(x, ecs, prox_table = TRUE, matrix_out = TRUE,
                             mDat = "NN", sr = FALSE, input = "SNP", chr.length = NULL,
                             levels = FALSE, par = FALSE, ss = FALSE, level_report = 1, maf = FALSE){

  # correct missing data format, since I later standardized.
  mDat <- substr(mDat, 1, nchar(mDat)/2)

  #=====================sanity checks=============
  #sanity checks:
  if (prox_table == FALSE & matrix_out == FALSE){
    stop("Please specify output format.")
  }
  if(is.numeric(ss) & ss <= 0){
    stop("Number/proportion of sites to subsample must be larger than 0.")
  }
  else if(!(is.numeric(ss)) & ss != FALSE){
    stop("Unaccepted ss. Please provide a numeric value or set to FALSE.")
  }

  if(par != FALSE & !is.numeric(par)){
    stop("Par must be FALSE or an integer.")
  }
  if(round(par) != par){
    stop("Par must be an integer.")
  }
  if(is.numeric(par)){
    if(par > parallel::detectCores()){
      stop("Par must be equal to or less than the number of available cores.")
    }
  }


  if(is.numeric(maf)){
    if((maf >= 1 & maf <= 0) | !(input %in% c("ms", "haplotype"))){
      stop("Minor allele frequency filtering is only for ms or haplotype input and requires that the provided maf be a numeric value between 0 and 1.\nTo use other filters or filter other input data, see filter_snps().")
    }
  }
  else if(maf != FALSE & !is.numeric(maf)){
    stop("Minor allele frequency filtering is only for ms or haplotype input and requires that the provided maf be a numeric value between 0 and 1.\nTo use other filters or filter other input data, see filter_snps().")
  }

  if(!(input %in% c("haplotype", "ms", "SNP"))){
    stop("Only haplotype, ms or SNP inputs accepted.")
  }

  #=====================functions and final checks=============
  #functions to do stuff with MS or phased haplotypes.
  else if(input == "haplotype" | input == "ms"){
    library(data.table)

    #prepare input data
    if(input == "ms"){
      cat("Reading and formatting ms input file.\n")
      infile <- x #infile
      lines <- readLines(x)
      lines <- lines[-which(lines == "")] #remove empty entries
      lines <- lines[-c(1,2)] #remove header info
      nss <- grep("segsites", lines) #get the number of segsites per chr
      chrls <- gsub("segsites: ", "", lines[nss]) #parse this to get the lengths
      chrls <- as.numeric(chrls)
      lines <- lines[-nss] #remove the segsites lines
      pos <- lines[grep("positions:", lines)] #find the positions
      lines <- lines[-grep("positions:", lines)] #remove the position
      div <- grep("//", lines) #find the seperators
      gc <- div[2] - div[1] - 1 #find the number of gene copies per chr
      if(is.na(gc)){gc <- length(lines) - 1} #if there's only one chr
      dat <- lines[-div] #get the data only
      dat <- strsplit(dat, "") #split the lines by individual snp calls
      x <- matrix(NA, nrow = sum(chrls), ncol = gc) #prepare output
      meta <- matrix(NA, nrow = sum(chrls), 2)

      #process this into workable data
      pchrls <- c(0, chrls)
      pchrls <- cumsum(pchrls)
      cat("Processing chromsomes/regions:")
      for(i in 1:length(chrls)){
        cat("\n\tChr ", i)
        tg <- dat[(gc*(i-1) + 1):(gc*i)] #get only this data
        tg <- unlist(tg) #unlist
        tg <- matrix(tg, ncol = chrls[i], nrow = gc, byrow = T) #put into a matrix
        tg <- t(tg) #transpose. rows are now snps, columns are gene copies
        tpos <- unlist(strsplit(pos[i], " ")) #grap and process the positions
        tpos <- tpos[-1]
        meta[(pchrls[i] + 1):pchrls[i + 1],] <- cbind(paste0(rep("chr", length = nrow(tg)), i), tpos)
        x[(pchrls[i] + 1):pchrls[i + 1],] <- tg #add data to output
      }
      cat("\n")

      meta <- as.data.frame(meta, stringsAsFactors = F)
      meta[,2] <- as.numeric(meta[,2])
      meta[,2] <- meta[,2] * chr.length

      colnames(meta) <- c("group", "position")
      colnames(x) <- paste0("gc_", 1:ncol(x))
      rm(tg, dat, pchrls, pos, tpos, lines, div, gc, i, chrls, infile)
      gc()
      ecs <- 2
    }
    else{
      meta <- x[,1:ecs]
      x <- as.matrix(x[,(ecs+1):ncol(x)])
    }

    #do filtering (for minor allele frequency) if requested.
    if(is.numeric(maf)){
      cat("Filtering sites with low minor allele frequencies. Starting sites: ", nrow(x), "\n")
      mafs <- matrixStats::rowSums2(matrix(as.numeric(x), nrow(x))/ncol(x))
      mafs[mafs > 0.5] <- 1 - mafs[mafs > 0.5]
      vio <- which(mafs < maf)
      x <- x[-vio,]
      meta <- meta[-vio,]
      cat("Ending sites: ", nrow(x), "\n")
    }

    #functions
    tabulate_haplotypes <- function(x){
      thaps <- matrix(paste0(x[1,], t(x[-1,])), ncol = nrow(x) - 1) #convert each cell to haplotype vs row one
      mthaps <- reshape2::melt(thaps) #put this into long form for tabulation
      cnames <- levels(mthaps$value)
      hapmat <- bigtabulate::bigtabulate(mthaps, ccols = c(2,3))
      colnames(hapmat) <- cnames
      return(hapmat)
    }
    LD_func <- function(x, meta, prox_table = TRUE, matrix_out = TRUE, mDat = "N", sr = FALSE, chr.length = NULL, stop.row = nrow(x) - 1){
      #function to count the number of haplotypes vs every other given position
      #################
      #prep stuff
      if(prox_table){
        ncomps <- (nrow(x)*(nrow(x)-1)/2) - (nrow(x)-stop.row)*((nrow(x)-stop.row)-1)/2
        prox <- data.table::as.data.table(data.frame(p1 = numeric(ncomps),
                                                     p2 = numeric(ncomps),
                                                     rsq = numeric(ncomps),
                                                     Dprime = numeric(ncomps),
                                                     pval = numeric(ncomps)))
      }
      if(matrix_out){
        rmat <- matrix(as.numeric(NA), stop.row, nrow(x) - 1)
        rmat <- data.table::as.data.table(rmat)
        colnames(rmat) <- as.character(meta$position[-1])
        rownames(rmat) <- make.names(as.character(meta$position[1:stop.row]), unique = T)
        Dpmat <- data.table::copy(rmat)
        pvmat <- data.table::copy(rmat)
      }
      if(!matrix_out & !prox_table){
        stop("Please specify output format.\n")
      }

      #run length prediction variables for progress reporting
      compfun <- function(x){
        return(((x-1)*x)/2)
      }
      totcomp <- compfun(nrow(x))
      prog <- 0
      cpercent <- 0

      ##################
      #calculate Dprime, rsq, ect.
      for(i in 1:stop.row){
        prog_after <- prog + nrow(x) - i

        #report progress
        if(!sr){
          cprog <- (totcomp - compfun(nrow(x) - i - 1))/totcomp
          if(cprog >= 0.01 + cpercent){
            cat("Progress:", paste0(round(cprog*100), "%."), "\n")
            cpercent <- cprog
          }
        }

        #check that this site isn't fixed.
        if(length(unique(x[i,])) == 1){
          if(prox_table){
            data.table::set(prox, (prog + 1):prog_after, j = "p1", value = meta$position[i])
            data.table::set(prox, (prog + 1):prog_after, j = "p2", value = meta$position[(i+1):nrow(x)])
            data.table::set(prox, (prog + 1):prog_after, j = "rsq", value = NA)
            data.table::set(prox, (prog + 1):prog_after, j = "Dprime", value = NA)
            data.table::set(prox, (prog + 1):prog_after, j = "pval", value = NA)
          }
          prog <- prog_after
          next()
        }

        #get haplotypes
        haps <- tabulate_haplotypes(x[i:nrow(x),])

        #fix for very rare cases
        if(!is.matrix(haps)){
          haps <- t(as.matrix(haps))
        }

        #Fix only three haplotypes. While Dprime is 1, can't just set rsq, so fix the table and continue.
        if(ncol(haps) == 3){
          pos.a <- unique(unlist(unlist(strsplit(colnames(haps), ""))))
          if(!(any(colnames(haps) == paste0(pos.a[1], pos.a[1])))){
            tnames <- colnames(haps)
            if(nrow(haps) == 1){
              haps <- t(as.matrix(c(0, haps)))
            }
            else{
              haps <- cbind(numeric(nrow(haps)), haps)
            }
            colnames(haps)<- c(paste0(pos.a[1], pos.a[1]), tnames)
          }
          else if(!(any(colnames(haps) == paste0(pos.a[1], pos.a[2])))){
            tnames <- colnames(haps)
            if(nrow(haps) == 1){
              haps <- t(as.matrix(c(haps[1], 0, haps[2:3])))
            }
            else{
              haps <- cbind(haps[,1], numeric(nrow(haps)), haps[,2:3])
            }
            colnames(haps) <- c(tnames[1], paste0(pos.a[1], pos.a[2]), tnames[2:3])
          }
          else if(!(any(colnames(haps) == paste0(pos.a[2], pos.a[1])))){
            tnames <- colnames(haps)
            if(nrow(haps) == 1){
              haps <- t(as.matrix(c(haps[1:2], 0, haps[3])))
            }
            else{
              haps <- cbind(haps[,1:2], numeric(nrow(haps)), haps[,3])
            }
            colnames(haps) <-  c(tnames[1:2], paste0(pos.a[2], pos.a[1]), tnames[3])
          }
          else{
            tnames <- colnames(haps)
            if(nrow(haps) == 1){
              haps <- t(as.matrix(c(haps, 0)))
            }
            else{
              haps <- cbind(haps, numeric(nrow(haps)))
            }
            colnames(haps) <- c(tnames, paste0(pos.a[2], pos.a[2]))
          }
        }

        #fix two haplotype cases (NA if either snp is fixed, otherwise 1)
        if(ncol(haps) == 2){
          cn <- colnames(haps)
          check <- ifelse(substr(cn[1], 1, nchar(cn[1])/2) !=
                   substr(cn[2], 1, nchar(cn[1])/2) &
                   substr(cn[1], (nchar(cn[1])/2) + 1, nchar(cn[1])) !=
                   substr(cn[2], (nchar(cn[1])/2) + 1, nchar(cn[1])),
                 1,NA)

          if(is.null(nrow(cn))){
            Dprime <- rep(check, 1)
            rsq <- rep(check, 1)
          }
          else{
            Dprime <- rep(check, nrow(cn))
            rsq <- rep(check, nrow(cn))
          }
          chisqu <- rsq*(4)
          pval <- 1 - pchisq(chisqu, 1)

          #write:
          if(prox_table){
            data.table::set(prox, (prog + 1):prog_after, j = "p1", value = meta$position[i])
            data.table::set(prox, (prog + 1):prog_after, j = "p2", value = meta$position[(i+1):nrow(x)])
            data.table::set(prox, (prog + 1):prog_after, j = "rsq", value = rsq)
            data.table::set(prox, (prog + 1):prog_after, j = "Dprime", value = Dprime)
            data.table::set(prox, (prog + 1):prog_after, j = "pval", value = pval)
          }
          if(matrix_out){
            #reminder: columns start at locus two, rows start at locus one (but end at nlocus - 1)
            fill <- rep(NA, nrow(x) - length(Dprime) - 1)
            Dprime <- c(fill, Dprime)
            rsq <- c(fill, rsq)
            pval <- c(fill, pval)
            data.table::set(rmat, i, 1:ncol(rmat), as.list(rsq))
            data.table::set(Dpmat, i, 1:ncol(Dpmat), as.list(Dprime))
            data.table::set(pvmat, i, 1:ncol(pvmat), as.list(pval))
            rm(pval, rsq, Dprime)
          }
          prog <- prog_after
          next()
        }

        #fix three missing haplotypes (everything NA)
        if(ncol(haps) == 1){

          Dprime <- rep(NA, nrow(haps))
          rsq <- rep(NA, nrow(haps))
          pval <- rep(NA, nrow(haps))

          #write
          if(prox_table){
            data.table::set(prox, (prog + 1):prog_after, j = "p1", value = meta$position[i])
            data.table::set(prox, (prog + 1):prog_after, j = "p2", value = meta$position[(i+1):nrow(x)])
            data.table::set(prox, (prog + 1):prog_after, j = "rsq", value = rsq)
            data.table::set(prox, (prog + 1):prog_after, j = "Dprime", value = Dprime)
            data.table::set(prox, (prog + 1):prog_after, j = "pval", value = pval)
          }
          if(matrix_out){
            #reminder: columns start at locus two, rows start at locus one (but end at nlocus - 1)
            fill <- rep(NA, nrow(x) - length(Dprime) - 1)
            Dprime <- c(fill, Dprime)
            rsq <- c(fill, rsq)
            pval <- c(fill, pval)
            data.table::set(rmat, i, 1:ncol(rmat), as.list(rsq))
            data.table::set(Dpmat, i, 1:ncol(Dpmat), as.list(Dprime))
            data.table::set(pvmat, i, 1:ncol(pvmat), as.list(pval))
            rm(pval, rsq, Dprime)
          }
          prog <- prog_after
          next()
        }

        #calc Dprime, rsq
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
        Dprime <- ifelse(D > 0, D/matrixStats::rowMins(cbind(A1f*B2f, A2f*B1f)),
                         ifelse(D < 0, D/matrixStats::rowMaxs(cbind(-A1f*B1f, -A2f*B2f)),
                                0))
        rsq <- (D^2)/(A1f*A2f*B1f*B2f)

        #fix for when more missing haps.
        Dprime[which(rowSums(ifelse(haps == 0, T, F)) == 3)] <- 0 #if three missing haplotypes, call 0.
        Dprime[which(rowSums(ifelse(haps == 0, T, F)) == 4)] <- NA # if four, call NA.
        rsq[which(rowSums(ifelse(haps == 0, T, F)) %in% 3:4)] <- NA #replace rsq with NA when 3 or 4 missing haplotypes (the latter shouldn't ever happen without missing data.)



        #if two missing, harder:
        miss2 <- rowSums(ifelse(haps == 0, T, F)) == 2
        if(any(miss2)){
          #get the missing haplotypes in each row:
          tm_mat <- haps[which(miss2),] #grab the violating rows
          tm_mat[tm_mat != 0] <- NA
          tm_mat[tm_mat == 0] <- 1
          tm_mat <- t(t(tm_mat)*(1:ncol(haps)))
          tm_mat[!is.na(tm_mat)] <- colnames(haps)[as.numeric(tm_mat[!is.na(tm_mat)])]
          tm_mat <- matrix(t(tm_mat)[!is.na(t(tm_mat))], ncol = 2, byrow = T)

          #if both are actually polymorphic, assign a one, otherwise asign a 0.
          check <- ifelse(substr(tm_mat[,1], 1, nchar(tm_mat[1,1])/2) !=
                            substr(tm_mat[,2], 1, nchar(tm_mat[1,1])/2) &
                            substr(tm_mat[,1], (nchar(tm_mat[1,1])/2) + 1, nchar(tm_mat[1,1])) !=
                            substr(tm_mat[,2], (nchar(tm_mat[1,1])/2) + 1, nchar(tm_mat[1,1])),
                          1,NA)
          Dprime[which(miss2)] <- check
          rsq[which(miss2)] <- check
          rm(tm_mat)
        }

        #get pvals
        chisqu <- rsq*(4)
        pval <- 1 - pchisq(chisqu, 1)

        #remove stuff to clear memory
        rm(A1B1f, A1B2f, A2B1f, A2B2f, haps, miss2, D, D2, chisqu, B1f, B2f, A1f, A2f)

        #write
        if(prox_table){
          data.table::set(prox, (prog + 1):prog_after, j = "p1", value = meta$position[i])
          data.table::set(prox, (prog + 1):prog_after, j = "p2", value = meta$position[(i+1):nrow(x)])
          data.table::set(prox, (prog + 1):prog_after, j = "rsq", value = rsq)
          data.table::set(prox, (prog + 1):prog_after, j = "Dprime", value = Dprime)
          data.table::set(prox, (prog + 1):prog_after, j = "pval", value = pval)
        }
        if(matrix_out){
          #reminder: columns start at locus two, rows start at locus one (but end at nlocus - 1)
          fill <- rep(NA, nrow(x) - length(Dprime) - 1)
          Dprime <- c(fill, Dprime)
          rsq <- c(fill, rsq)
          pval <- c(fill, pval)
          data.table::set(rmat, i, 1:ncol(rmat), as.list(rsq))
          data.table::set(Dpmat, i, 1:ncol(Dpmat), as.list(Dprime))
          data.table::set(pvmat, i, 1:ncol(pvmat), as.list(pval))
        }
        prog <- prog_after
      }

      ###################################
      #finish and return output
      if(prox_table){
        prox$proximity <- abs(prox$p1 - prox$p2)
        if(any(colnames(meta) == "group")){
          prox$group <- meta[1,"group"]
        }
        if(any(colnames(meta) == "pop")){
          prox$pop <- meta[1,"pop"]
        }
        prox <- prox[,c(6:ncol(prox), 1:5)]
      }
      if(matrix_out){
        Dpmat <- cbind(position = meta$position[1:stop.row], Dpmat)
        rmat <- cbind(position = meta$position[1:stop.row], rmat)
        pvmat <- cbind(position = meta$position[1:stop.row], pvmat)
      }
      if(prox_table){
        if(matrix_out){
          return(list(prox = prox, Dprime = Dpmat, rsq = rmat, pval = pvmat))
        }
        else{
          return(prox = prox)
        }
      }
      else{
        return(list(Dprime = Dpmat, rsq = rmat, pval = pvmat))
      }
    }
  }


  #function to do LD with SNPs
  else{
    #prepare input and metadata.
    meta <- x[,1:ecs]
    x <- x[,(ecs + 1):ncol(x)]
    x <- as.matrix(x)

    library(dplyr)

    #functions:
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
    tabulate_haplotypes <- function(x, y, as, dmDat, sform){
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
      mgcv <- reshape2::melt(gcv)
      cnames <- levels(mgcv$value)
      ghapmat <- bigtabulate::bigtabulate(mgcv, ccols = which(colnames(mgcv) %in% c("Var1", "value")))
      colnames(ghapmat) <- cnames

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
      # missin haplotypes at the end.
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
    LD_func <- function(x, meta, prox_table = TRUE, matrix_out = TRUE, mDat = "N", sr = FALSE, chr.length = NULL, stop.row = nrow(x) - 1){
      #double check that the position variable is numeric!
      if(!is.numeric(meta$position)){meta$position <- as.numeric(meta$position)}

      dmDat <- paste0(mDat, mDat) #missing genotype
      #data format
      sform <- nchar(x[1,1])/2

      #get unique alleles present at each locus (note: this made a much quicker a tab... should implement this elsewhere...)
      #note, tabulate_haplotypes needs this...
      p1 <- substr(as.matrix(x), 1, sform)
      p2 <- substr(as.matrix(x), sform + 1, sform*2)
      pc <- reshape2::melt(cbind(p1,p2))
      #amat <- with(pc, table(Var1, value))
      #amat <- amat[,colnames(amat) != mDat]
      as <- sort(unique(pc$value))
      as <- as.character(as[as != mDat])


      #need to loop through each loci and compare to everything else. Probably can't really vectorize the outer loop.



      #initialize output.
      if(prox_table){
        prox <- data.frame(p1 = numeric(1), p2 = numeric(1), rsq = numeric(1), Dprime = numeric(1), pval = numeric(1))
      }
      if(matrix_out){
        rmat <- matrix(NA, stop.row, nrow(x) - 1)
        colnames(rmat) <- meta$position[-1]
        rownames(rmat) <- meta$position[1:stop.row]
        Dpmat <- rmat
        pvmat <- rmat
      }

      #run length prediction variables for progress reporting
      compfun <- function(x){
        return(((x-1)*x)/2)
      }
      totcomp <- compfun(nrow(x))
      cpercent <- 0

      #loop through and get haplotypes, calc LD for each locus.
      for(i in 1:stop.row){
        if(!sr){
          cprog <- (totcomp - compfun(nrow(x) - i - 1))/totcomp
          if(cprog >= 0.05 + cpercent){
            cat("Progress:", paste0(round(cprog*100), "%."), "\n")
            cpercent <- cprog
          }
        }
        haps <- tabulate_haplotypes(x[i,], x[(i+1):nrow(x),], as, dmDat, sform)

        #if we had only one haplotype or no haplotypes:
        if(is.na(haps[1])){
          if(prox_table){
            prox <- rbind(prox,
                          cbind(p1 = meta$position[i], p2 = meta$position[(i+1):nrow(x)],
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
                          cbind(p1 = meta$position[i], p2 = meta$position[(i+1):nrow(x)],
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
        Dprime <- ifelse(D > 0, D/matrixStats::rowMins(cbind(A1f*B2f, A2f*B1f)),
                         ifelse(D < 0, D/matrixStats::rowMaxs(cbind(-A1f*B1f, -A2f*B2f)),
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
                        cbind(p1 = meta$position[i], p2 = meta$position[(i+1):nrow(x)],
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
        prox$proximity <- abs(prox$p1 - prox$p2)
        if(any(colnames(meta) == "group")){
          prox$group <- meta[1,"group"]
        }
        if(any(colnames(meta) == "pop")){
          prox$pop <- meta[1,"pop"]
        }
        prox <- prox[,c(6:ncol(prox), 1:5)]
      }
      if(prox_table){
        if(matrix_out){
          return(list(prox = prox, Dprime = Dpmat, rsq = rmat, pval = pvmat))
        }
        else{
          return(prox = prox)
        }
      }
      else{
        return(list(Dprime = Dpmat, rsq = rmat, pval = pvmat))
      }
    }
  }


  #one more sanity check.
  if(ss > nrow(x)){
    stop("Number of sites to subsample cannot be large than the number of provided sites.")
  }

  #subset data if requested:
  if(is.numeric(ss)){
    #get sample
    if(ss <= 1){
      ss <- sample(nrow(x), round(nrow(x)*ss), F)
    }
    else{
      ss <- sample(nrow(x), ss, F)
    }

    #subset
    meta <- meta[ss,]
    x <- x[ss,]

    #re-sort
    meta$ord <- 1:nrow(meta)
    meta <- dplyr::arrange(meta, position)
    x <- x[meta$ord,]
    meta <- meta[,-ncol(meta)]
  }



  #=====================call functions=========
  #call these functions (possibly in parallel) according to supplied levels.

  cat("Beginning LD calculation...\n")

  #no facets
  if(levels == FALSE){
    cat("No levels specified.\n")

    #run in parallel if requested
    if(is.numeric(par)){
      cat("Running in parallel.\n\t")
      #this is tricky. Needs to be fed essentially the same data (although it can start farther along), but have different stopping points.
      #then need to knit the matrix outputs together. The prox outputs are easy, just pass them along.
      split <- (nrow(x)*(nrow(x)-1)/2)/par #number of comparisons to do per core
      split <- ceiling(split)
      cat("At least", split, "pairwise comparisons per processor.\n")

      #vector of comps per row
      nc.pr <- (nrow(x) - 1):0 #since the number is always this row vs every other.
      snc.pr <- cumsum(nc.pr) #cumulative sum

      #which rows does each processor do?
      rproc <- matrix(0, nrow = par, 2) #i is the proc, col 1 is the starting row, col 2 is the ending row
      cr <- 1
      for(i in 1:par){
        rproc[i,1] <- cr
        if(i == par){
          rproc[i,2] <- length(nc.pr) - 1
        }
        else{
          rproc[i,2] <- which(cumsum(nc.pr) >= split*i)[1]
        }
        cr <- rproc[i,2] + 1
      }

      #now need to start the parallel job:
      library(doParallel);library(foreach)
      cl <- snow::makeSOCKcluster(par)
      doSNOW::registerDoSNOW(cl)

      #prepare reporting function
      ntasks <- par
      progress <- function(n) cat(sprintf("Part %d out of",n), ntasks, "is complete.\n")
      opts <- list(progress=progress)


      #save the info as a bigmatrix if it can be safely converted to numeric. Usually this is true for ms but not necissarily other data types.
      if(as.numeric(x[1]) == x[1]){
        cat("Saving matrix as big.matrix object for quicker sharing.\n")
        xb <- bigmemory::as.big.matrix(x, type = "char")
        xbd <- bigmemory::describe(xb)
        remove(x)
      }

      cat("Begining run.\n")
      # for(q in 1:par){
      #   w_data <- LD_func(x = xb[rproc[q,1]:nrow(xb),],
      #                     meta = meta[rproc[q,1]:nrow(x),], prox_table = prox_table, mDat = mDat, matrix_out = matrix_out,
      #                     sr = T, chr.length = chr.length, stop.row = rproc[q,2] - (rproc[q,1]-1))
      # }

      #run the jobs
      output <- foreach::foreach(q = 1:ntasks, .packages = c("bigmemory", "dplyr"), .inorder = TRUE,
                                 .options.snow = opts) %dopar% {
                                   if(exists("xbd")){
                                     x <- bigmemory::attach.big.matrix(xbd)
                                   }
                                   #note that matrix out is FALSE, since we have to build one up from the prox table anyway...
                                   w_data <- LD_func(x = x[rproc[q,1]:nrow(x),],
                                                     meta = meta[rproc[q,1]:nrow(x),], prox_table = prox_table, mDat = mDat, matrix_out = matrix_out,
                                                     sr = T, chr.length = chr.length, stop.row = rproc[q,2] - (rproc[q,1]-1))
                                 }

      #release cores
      parallel::stopCluster(cl)
      doSNOW::registerDoSNOW()

      cat("LD computation completed. Preparing results.\n\t")

      #combine the results
      if(matrix_out){
        cat("Standardizing matrices. Batch:")
        #standardize column and row names, melt
        fnames <- meta$position
        fnames <- paste0(fnames, "_", 1:length(fnames))
        cnames <- fnames[-1]
        rnames <- fnames[-length(fnames)]

        for(i in 1:length(output)){
          cat(paste0("\n\t\t", i))

          #fixes for matrices and objects without a position column. Mostly for the snp data format.
          if(!data.table::is.data.table(output[[i]]$Dprime)){
            output[[i]]$Dprime <- data.table::as.data.table(output[[i]]$Dprime)
            output[[i]]$rsq <- data.table::as.data.table(output[[i]]$rsq)
            output[[i]]$pval <- data.table::as.data.table(output[[i]]$pval)
          }

          if(!any(colnames(output[[i]]$Dprime) == "position")){
            output[[i]]$Dprime <- cbind(position = rep(NA, nrow(output[[i]]$Dprime)), output[[i]]$Dprime)
            output[[i]]$rsq <- cbind(position = rep(NA, nrow(output[[i]]$Dprime)), output[[i]]$rsq)
            output[[i]]$pval <- cbind(position = rep(NA, nrow(output[[i]]$Dprime)), output[[i]]$pval)
          }

          #fix the names.
          colnames(output[[i]]$Dprime)[-1] <- cnames[rproc[i,1]:(nrow(x)-1)]
          colnames(output[[i]]$rsq) <- colnames(output[[i]]$Dprime)
          colnames(output[[i]]$pval) <- colnames(output[[i]]$Dprime)

          output[[i]]$Dprime$position <- rnames[rproc[i,1]:rproc[i,2]]
          output[[i]]$rsq$position <- output[[i]]$Dprime$position
          output[[i]]$pval$position <- output[[i]]$Dprime$position

          #melt
          output[[i]]$Dprime <- reshape2::melt(output[[i]]$Dprime, id.vars = "position")
          output[[i]]$rsq <- reshape2::melt(output[[i]]$rsq, id.vars = "position")
          output[[i]]$pval <- reshape2::melt(output[[i]]$pval, id.vars = "position")
        }
        cat("\n\tDone.\n")
      }

      #bind everything together
      cat("Combining results. Batch:")
      w_list <- output[[1]]
      cat("\n\t1")
      for (i in 2:length(output)){
        cat(paste0("\n\t\t", i))
        if(prox_table){
          if(matrix_out){
            w_list$prox <- rbind(w_list$prox, output[[i]]$prox)
          }
          else{
            w_list <- rbind(w_list, output[[i]])
          }
        }
        if(matrix_out){
          w_list$Dprime <- rbind(w_list$Dprime, output[[i]]$Dprime)
          w_list$rsq <- rbind(w_list$rsq, output[[i]]$rsq)
          w_list$pval <- rbind(w_list$pval, output[[i]]$pval)
        }
      }
      cat("\n\tDone.\n")

      #re-cast the output matrix
      if(matrix_out){
        cat("Re-casting output matrices.\n")
        w_list$Dprime <- reshape2::acast(w_list$Dprime, formula = position ~ variable)
        w_list$rsq <- reshape2::acast(w_list$rsq, formula = position ~ variable)
        w_list$pval <- reshape2::acast(w_list$pval, formula = position ~ variable)

        #fix columns and positions
        colnames(w_list$Dprime) <- gsub("_.+", "", colnames(w_list$Dprime))
        colnames(w_list$rsq) <- gsub("_.+", "", colnames(w_list$rsq))
        colnames(w_list$pval) <- gsub("_.+", "", colnames(w_list$pval))

        rownames(w_list$Dprime) <- gsub("_.+", "", rownames(w_list$Dprime))
        rownames(w_list$rsq) <- gsub("_.+", "", rownames(w_list$rsq))
        rownames(w_list$pval) <- gsub("_.+", "", rownames(w_list$pval))
        cat("\tDone.\n")
      }

      return(w_list)
    }

    #otherwise run normally
    else{
      return(LD_func(x, meta, prox_table, matrix_out, mDat, sr, chr.length))
    }
  }


  #prepare storage
  #find facets to loop through
  facets <- unique(meta[,colnames(meta) %in% levels])
  facets <- as.matrix(facets)

  #prepare output list
  w_list<- vector("list", length = prox_table + nrow(facets)*matrix_out)
  if(prox_table){names(w_list)[1] <- "prox"}
  if(matrix_out){
    if(prox_table){
      names(w_list)[-1] <- apply(facets, 1, paste, collapse = "_")
    }
    else{
      names(w_list) <- apply(facets, 1, paste, collapse = "_")
    }
  }

  #not in parallel
  if(par == FALSE){

    #loop through each set of facets
    for (i in 1:nrow(facets)){
      #grab the data matching these facets...
      matches <- which(apply(meta, 1, function(x) identical(as.character(x[colnames(meta) %in% levels]), facets[i,])))
      w_data <- x[matches,]
      w_meta <- meta[matches,]

      #if there are no or one row in data, skip it.
      if(nrow(w_data) == 0 | nrow(w_data) == 1){
        next()
      }

      #report progress
      if(i %% level_report == 0){cat("Facet #:", i, "of", nrow(facets), " Name:", paste0(facets[i,], collapse = " "), "\n")}

      #do the pariwise LD using the correct LD func
      out <-LD_func(x = w_data, meta = w_meta, prox_table = prox_table, mDat = mDat, matrix_out = matrix_out, sr = sr, chr.length = chr.length)

      #add the correct data
      if(prox_table == TRUE){
        w_list$prox <- rbind(w_list$prox, out$prox)
      }
      if(matrix_out == TRUE){
        w_list[[i + 1]] <- out[-1]
      }
    }

    #return
    return(w_list)
  }

  #in parallel
  else{
    library(doParallel);library(foreach)
    cl <- snow::makeSOCKcluster(par)
    doSNOW::registerDoSNOW(cl)

    #prepare reporting function
    ntasks <- nrow(facets)
    progress <- function(n) cat(sprintf("Facet %d out of",n), ntasks, "is complete.\n")
    opts <- list(progress=progress)

    #loop through each set of facets
    output <- foreach::foreach(i = 1:ntasks, .packages = c("dplyr", "reshape2", "matrixStats", "bigtabulate"), .inorder = TRUE,
                               .options.snow = opts) %dopar% {

                                 matches <- which(apply(meta, 1, function(x) identical(as.character(x[colnames(meta) %in% levels]), facets[i,])))
                                 w_data <- x[matches,]
                                 w_meta <- meta[matches,]

                                 if(nrow(w_data) > 1){
                                   w_data <- LD_func(x = w_data, meta = w_meta, prox_table = prox_table, mDat = mDat, matrix_out = matrix_out, sr = sr, chr.length = chr.length)
                                 }
                               }


    #release cores
    parallel::stopCluster(cl)
    doSNOW::registerDoSNOW()

    #make the output sensible
    for(i in 1:length(output)){
      ##prox table
      if(prox_table){
        if(!matrix_out){
          w_list$prox <- rbind(w_list$prox, output[[i]])
        }
        else{
          w_list$prox <- rbind(w_list$prox, output[[i]][1]$prox)
        }
      }
      ##matrices
      if(matrix_out){
        if(!prox_table){
          w_list[[i]] <- output[[i]]
        }
        else{
          w_list[[i + 1]] <- output[[i]][-1]
        }
      }
    }


    #return
    return(w_list)
  }
}



