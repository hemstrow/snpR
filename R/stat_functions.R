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
  func <- function(x){
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

  # add any missing facets
  facets <- check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    invisible(capture.output(x <- add.facets.snpR.data(x, facets)))
  }

  out <- apply.snpR.facets(x, facets, "ac", func, case = "ps")
  colnames(out)[ncol(out)] <- "pi"
  x <- merge.snpR.stats(x, out)

  return(x)
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
  facets <- check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    invisible(capture.output(x <- add.facets.snpR.data(x, facets)))
  }
  out <- apply.snpR.facets(x,
                           facets = facets,
                           req = "gs",
                           fun = func,
                           case = "ps",
                           m.al = substr(x@mDat,1, nchar(x@mDat)/2))
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
#' @param method Which FST estimator should be used?
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

calc_pairwise_fst <- function(x, facets, method = "WC"){
  func <- function(x, method, facets = NULL){
    pops <- sort(unique(x$subfacet))
    npops <- length(pops)
    pnk.length <- (npops*(npops-1))/2 * (nrow(x)/npops)
    pnk <- data.table::data.table(comparison = character(pnk.length), ntotal = integer(pnk.length))

    #===============genepop======================
    if(method == "genepop"){
      g.filename <- paste0("genepop.", facets, ".txt")
      if(file.exists(g.filename)){
        cat("Genepop input file", g.filename,  "already exits. ")
        resp <- "empty"
        while(resp != "y" & resp != "n"){
          cat("Overwrite? (y or n)\n")
          resp <- readLines(n = 1)
        }
        if(resp == "n"){
          genepop::Fst(g.filename, pairs = TRUE)
        }
        else{
          file.remove(g.filename)
          invisible(capture.output(format_snps(x, output = "genepop", facets = facets, outfile = g.filename)))
          genepop::Fst(g.filename, pairs = TRUE)
        }
      }
      else{
        invisible(capture.output(format_snps(x, output = "genepop", facets = facets, outfile = g.filename)))
        genepop::Fst(g.filename, pairs = TRUE)
      }



      #read in and parse genepop output
      cat("Parsing genepop output...\n")

      #read the file in
      y <- paste0(g.filename, ".ST2") #data file
      y <- readLines(y)

      #get the number of pops and the number of loci.
      np <- grep("Number of populations detected", y)
      np <- as.numeric(unlist(strsplit(y[np], " : "))[2])
      nl <- grep("Number of loci detected", y)
      nl <- as.numeric(unlist(strsplit(y[nl], " : "))[2])

      #check that the correct number of pop names were provided or grab new ones if they weren't.
      ##get pop data
      # py <- grep("Indices for populations:", y)
      # py <- y[(py+2):(py+1+np)]
      # py <- unlist(strsplit(py, " +"))
      # py <- py[seq(2, length(py), 2)]

      # population names
      pnames <- sort(unique(x@facet.meta$subfacet[x@facet.meta$facet == facets]))


      #get the indices containing locus headers
      locs <- grep("  Locus:", y)
      locs <- c(locs, grep("Estimates for all loci", y))

      #get indices not containing data to parse and remove them.
      empts <- c(1:(locs[1]-2), locs, locs + 1, locs + 2, locs - 1, (length(y)-2):length(y))
      vals <- y[-empts]

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
        colnames(tmat) <- paste0(pnames[1:(i-1)], "~", pnames[i]) #name the comparison in each column
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

      # melt
      suppressMessages(out$loci <- reshape2::melt(out$loci))
      colnames(out$loci) <- c("comparison", "fst")

      # get nk values:
      n_tots <- data.table::data.table(pop = x@facet.meta$subfacet[x@facet.meta$facet == facets],
                                       .snp.id = x@facet.meta$.snp.id[x@facet.meta$facet == facets],
                                       n_total = x@ac[x@facet.meta$facet == facets, "n_total"])
      n_tots <- data.table::dcast(n_tots, .snp.id ~ pop, value.var = "n_total")
      n_tots <- n_tots[,-1]


      prog <- 1
      for(i in 1:(ncol(n_tots) - 1)){
        for(j in (i+1):ncol(n_tots)){
            pnk <- set(pnk, prog:(prog + nrow(n_tots) - 1), 1L, paste0(colnames(n_tots)[i], "~", colnames(n_tots)[j]))
            pnk <- set(pnk, prog:(prog + nrow(n_tots) - 1), 2L, as.integer(n_tots[,i] + n_tots[,j]))
            prog <- prog + nrow(n_tots)
        }
      }

      out$loci <- pnk$n_total

      # return, we're done.
      cat("Finished.\n")
      return(out)
    }

    #===============others=====================
    browser()
    x <- data.table::as.data.table(x)
    data.table::setkey(x, subfacet)


    out <- data.table::as.data.table(matrix(NA, ncol = (length(pops)*(length(pops) - 1)/2), nrow = nrow(x)/length(pops)))

    #initialize pop comparison columns.
    comps <- c()
    i <- 1
    while (i < (length(pops))){
      j <- 1 + i
      for (j in j:length(pops)){
        comps <- c(comps, paste0(pops[i], "~", pops[j]))
        j <- j + 1
      }
      i <- i + 1
    }
    i <- 1

    colnames(out) <- comps

    #loop through each comparison and caculate pairwise FST at each site
    c.col <- 1
    prog <- 1
    for (i in 1:(length(pops) - 1)){ #i is the first pop
      system.time(temp <- x[x$subfacet == pops[i],]) #get data for first pop
      system.time(temp2 <- x[subfacet == pops[i],])


      idat <- x[x$subfacet == pops[i],] #get data for first pop
      j <- i + 1 #intialize j as the next pop
      for (j in j:length(pops)){#j is pop being compared
        jdat <- x[x$subfacet == pops[j],] #get data for second pop

        if(method == "wc" | method == "wier"){

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
          hbar <- ((idat$n_total*idat$ho) + (jdat$n_total*jdat$ho))/(r*nbar) #average heterozygote frequencies

          #equation parts used in both
          inner1 <- pbar*(1-pbar)
          inner2 <- ((r-1)/r)*ssq
          inner3 <- .25*hbar

          if(method == "wc"){
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

        else if(method == "hohenlohe"){

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

        # update pnk
        pnk <- set(pnk, prog:(prog + nrow(idat) - 1), 1L, paste0(idat$subfacet, "~", jdat$subfacet))
        pnk <- set(pnk, prog:(prog + nrow(idat) - 1), 2L, as.integer(idat$n_total + jdat$n_total))
        prog <- prog + nrow(idat)


        c.col <- c.col + 1 #agument c.col
      }
    }

    # melt, cbind pnk
    suppressMessages(out <- reshape2::melt(out))
    colnames(out) <- c("comparison", "fst")
    out$n_total <- pnk$n_total

    # return
    return(out)
  }

  #============================sanity and facet checks========================
  if(any(x@ac$n_alleles > 2)){
    vio <- which(x@ac$n_alleles[x@facet.meta$facet %in% facets] > 2)
    vio <- unique(x@facet.meta$.snp.id[x@facet.meta$facet %in% facets][vio])
    stop(cat("Some loci have more than two alleles. Violating loci:\n", paste0(vio, collapse = "\n")))
  }

  # add any missing facets
  facets <- check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    invisible(capture.output(x <- add.facets.snpR.data(x, facets)))
  }

  method <- tolower(method)


  # call apply.snpR.facets, slightly different for each method, since they require different stuff.
  if(method == "genepop"){
    out <- apply.snpR.facets(x, facets, req = "snpRdata", fun = func, case = "facet.pairwise", method = "genepop")
    ave.fst <- out[[2]]
    out <- out[[1]]
  }
  else if(method == "wc" | method == "wier"){
    x <- calc_ho(x, facets)
    out <- apply.snpR.facets(x, facets, req = c("ac.stats"), case = "facet.pairwise", fun = func, method = method)
  }
  else if(method == "hohenlohe"){
    x <- calc_pi(x, facets)
    out <- apply.snpR.facets(x, facets, req = c("ac.stats"), case = "facet.pairwise", fun = func, method = "hohenlohe")
  }

  # working here, just need to adjust merge.snpR.stats to work with pairwise stats.
  x <- merge.snpR.stats(x, out, type = "pairwise")

  if(method == "genepop"){
    message("Returning list: first element is the snpRdata object now containing pairwise fst values, second is the average pairwise fst for all comparisons.\n")
    return(list(x, ave.fst))
  }
  else{
    return(x)
  }


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
calc_ho <- function(x, facets = NULL){
  func <- function(gs){
    #identify heterozygote rows in genotype matrix
    genos <- colnames(gs$gs)
    hets <- which(substr(genos, 1, 1) != substr(genos, 2,2))

    # calculate ho
    ho <- rowSums(gs$gs[,hets])/rowSums(gs$gs)
  }

  # add any missing facets
  facets <- check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    invisible(capture.output(x <- add.facets.snpR.data(x, facets)))
  }

  out <- apply.snpR.facets(x,
                           facets = facets,
                           req = "gs",
                           fun = func,
                           case = "ps")
  colnames(out)[ncol(out)] <- "ho"
  return(merge.snpR.stats(x, out))
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
calc_private <- function(x, facets = NULL){
  func <- function(gs){
    # # things to add two kinds of private alleles for debugging purposes.
    # temp <- tail(gs$as)
    # temp$C <- c(0,0,0,0,0,43)
    # temp$.snp.id <- max(gs$as$.snp.id) + 1
    # temp$G <- c(temp$G[-6], 0)
    # temp2 <- head(gs$as)
    # temp2$A <- c(0,0,0,0,0,4)
    # gs$as <- rbind(gs$as, temp)
    # gs$as[1:6,] <- temp2

    out <- numeric(nrow(gs$as)) # initialize

    # no private alleles if only one level this facet
    if(length(unique(gs$as$subfacet)) == 1){
      return(out)
    }
    gs$as <- data.table::as.data.table(gs$as)

    # convert to logical, then melt down to long and cast back up to summarize the number of times each allele is observed across all populations in for each locus
    logi <- data.table::as.data.table(ifelse(gs$as[,4:ncol(gs$as)] == 0, F, T))
    logi <- cbind(gs$as[,1:3], logi)
    cgs <- data.table::melt(logi, id.vars = c("facet", "subfacet", ".snp.id"))
    cgs <- data.table::dcast(cgs, formula = .snp.id ~ variable, value.var = "value", fun.aggregate = sum)

    # find those with private alleles in any populations
    logi.cgs <- ifelse(cgs[,-1] == 1, T, F) # anything TRUE is a private allele in a population
    pa.loci <- which(rowSums(logi.cgs) != 0)


    if(length(pa.loci) != 0){
      # determine which population the private allele is in. Do so by first grabbing just the private allele logical
      # as and tabulated matrices. Then, make a comparison series that repeats the private allele series (T, F, F, F) for an A private allele for example)
      # once for each subfacet level. The row which matches this will be that for the populations that have the private allele!
      pa.cgs <- logi.cgs[pa.loci,]
      pa <- logi[logi$.snp.id %in% cgs$.snp.id[pa.loci],]
      comp.series <- rep(pa.cgs, each = length(unique(gs$as$subfacet)))
      has.private <- as.matrix(pa[,-c(1:3)])[comp.series] # here's where the private alleles are in the subset data.

      # mark as private in vector and return
      out[logi$.snp.id %in% cgs$.snp.id[pa.loci]][has.private] <- 1
    }
    # return
    return(out)
  }


  # add any missing facets
  facets <- check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    invisible(capture.output(x <- add.facets.snpR.data(x, facets)))
  }

  out <- apply.snpR.facets(x, facets, "meta.gs", func, case = "ps.pf")
  colnames(out)[ncol(out)] <- "pa"
  x <- merge.snpR.stats(x, out)
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
#' @param x snpRdata object
#' @param facets Facets to split by. Multi-level facets can be noted with a ".". "pop.family.chr", for example, will split by chr within each family and within each group/family level combo.
#' @param subfacets Subsets the facet levels to run. Given as a named list: list(fam = A) will run only fam A, list(fam = c("A", "B"), chr = 1) will run only fams A and B on chromosome 1. list(fam = "A", pop = "ASP") will run samples in either fam A or pop ASP, list(fam.pop = "A.ASP") will run only samples in fam A and pop ASP.
#' @param ss number of snps to subsample.
#' @param par number of parallel cores to run
#' @param sr should reports be supressed?
#'
#' @return Matrices containing rsq, Dprime, and p values for each SNP vs every other SNP. Can also produce a proximity table, which contains the rsq, Dprime, and p value for each pairwise comparison and the distance between the SNPs in those comparisons. Returns matrices and prox table as sequential elements in a named list.
#'
#' @examples
#' #returns prox table and LD matrices.
#' LD_full_pairwise(stickSNPs[stickSNPs$group == "groupI",1:53], ecs = 3)
#'
calc_pairwise_ld <- function(x, facets = NULL, subfacets = NULL, ss = FALSE,
                             par = FALSE, sr = FALSE){
  #========================sanity checks=============
  #sanity checks:
  # subsampling
  if(is.numeric(ss) & ss <= 0){
    stop("Number/proportion of sites to subsample must be larger than 0.")
  }
  else if(!(is.numeric(ss)) & ss != FALSE){
    stop("Unaccepted ss. Please provide a numeric value or set to FALSE.")
  }

  # parallelizing
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
  }  #one more sanity check.
  if(ss > nrow(x)){
    stop("Number of sites to subsample cannot be large than the number of provided sites.")
  }

  #========================sub-functions=============
  #function to do LD with SNPs
  library(dplyr)

  # sub functions:

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
    if(length(y) == length(x)){
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
    if(length(y) == length(x)){
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
    ##homozygotes:
    ### unless there are no double homozygotes:
    if(sum(substr(gl, 1, sform) == substr(gl, sform + 1, sform*2) &
           substr(gl, (sform*2) + 1, sform*3) == substr(gl, (sform*3+1), sform*4)) != 0){
      hapmat[,colnames(hapmat) %in% paste0(substr(colnames(dhom), 1, sform),
                                           substr(colnames(dhom),(sform*2)+1,sform*3))] <- dhom*2
    }

    ##heterozyogote locus 1
    ### unless locus one has no heterozygotes:
    if(sum(substr(gl, 1, sform) != substr(gl, sform + 1, sform*2)) != 0){
      n1 <- paste0(substr(colnames(het_l1), 1, sform),
                   substr(colnames(het_l1),(sform*2)+1,sform*3))
      n1 <- GtoH(het_l1, n1)
      n2 <- paste0(substr(colnames(het_l1),sform+1, sform*2),
                   substr(colnames(het_l1),(sform*3)+1, sform*4))
      n2 <- GtoH(het_l1, n2)
      hapmat[,colnames(hapmat) %in% colnames(n1)] <- n1 + hapmat[,colnames(hapmat) %in% colnames(n1)]
      hapmat[,colnames(hapmat) %in% colnames(n2)] <- n2 + hapmat[,colnames(hapmat) %in% colnames(n2)]
    }


    ##heterozyogote locus 2
    ### unless locus two has no heterozygotes
    if(sum(substr(gl, (sform*2) + 1, sform*3) != substr(gl, (sform*3+1), sform*4)) != 0){
      n1 <- paste0(substr(colnames(het_l2), 1, sform),
                   substr(colnames(het_l2),(sform*2)+1,sform*3))
      n1 <- GtoH(het_l2, n1)
      n2 <- paste0(substr(colnames(het_l2),sform+1, sform*2),
                   substr(colnames(het_l2),(sform*3)+1, sform*4))
      n2 <- GtoH(het_l2, n2)
      hapmat[,colnames(hapmat) %in% colnames(n1)] <- n1 + hapmat[,colnames(hapmat) %in% colnames(n1)]
      hapmat[,colnames(hapmat) %in% colnames(n2)] <- n2 + hapmat[,colnames(hapmat) %in% colnames(n2)]
    }


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

  # LD sub function, called in func
  LD_func <- function(x, meta, mDat, snp.list, sr = FALSE){
    smDat <- substr(mDat, 1, nchar(mDat)/2)

    # subset the requested samps
    x <- x[,snp.list$samps]

    if(!is.matrix(x)){x <- as.matrix(x)}

    #double check that the position variable is numeric!
    if(!is.numeric(meta$position)){meta$position <- as.numeric(meta$position)}

    #data format
    sform <- nchar(x[1,1])/2

    #get unique alleles present at each locus
    #note, tabulate_haplotypes needs this...
    p1 <- substr(x, 1, sform)
    p2 <- substr(x, sform + 1, sform*2)
    pc <- sort(unique(c(p1,p2)))
    as <- as.character(pc[pc != smDat])


    #need to loop through each loci and compare to everything else. Probably can't really vectorize the outer loop.



    #initialize output.
    prox <- matrix(NA, nrow = 0, ncol = 2*ncol(meta) + 4)
    colnames(prox) <- c(paste0("s1_", colnames(meta)), paste0("s2_", colnames(meta)), "proximity", "rsq", "Dprime", "pval")
    prox <- as.data.frame(prox)
    rmat <- matrix(NA, nrow(x), nrow(x))
    colnames(rmat) <- meta$position
    rownames(rmat) <- meta$position
    Dpmat <- rmat
    pvmat <- rmat


    #run length prediction variables for progress reporting
    compfun <- function(x){
      return(((x-1)*x)/2)
    }
    totcomp <- compfun(nrow(x))
    cpercent <- 0

    #loop through and get haplotypes, calc LD for each locus.
    for(i in 1:length(snp.list$snps)){
      if(!sr){
        cprog <- (totcomp - compfun(nrow(x) - i - 1))/totcomp
        if(cprog >= 0.05 + cpercent){
          cat("Progress:", paste0(round(cprog*100), "%."), "\n")
          cpercent <- cprog
        }
      }

      # get haplotypes
      if(is.null(snp.list$snps[[i]])){
        next()
      }
      haps <- tabulate_haplotypes(x[i,], x[snp.list$snps[[i]],], as, mDat, sform)

      #if we had only one haplotype or no haplotypes:
      if(is.na(haps[1])){
        tprox <- cbind(meta[i,],
                       meta[snp.list$snps[[i]],],
                       abs(meta[i,]$position - meta[snp.list$snps[[i]],]$position),
                       rsq = NA, Dprime = NA, pval = NA)

        colnames(tprox) <- colnames(prox)
        prox <- rbind(prox, tprox)

        #reminder: columns start at locus two, rows start at locus one (but end at nlocus - 1)
        rmat[i,snp.list$snps[[i]]] <- NA
        Dpmat[i,snp.list$snps[[i]]] <- NA
        pvmat[i,snp.list$snps[[i]]] <- NA
        next()
      }
      if(length(haps) == 1){
        tprox <- cbind(meta[i,],
                       meta[snp.list$snps[[i]],],
                       abs(meta[i,]$position - meta[snp.list$snps[[i]],]$position),
                       rsq = 0, Dprime = 0, pval = 0)

        colnames(tprox) <- colnames(prox)
        prox <- rbind(prox, tprox)

        #reminder: columns start at locus two, rows start at locus one (but end at nlocus - 1)
        rmat[i,snp.list$snps[[i]]] <- 0
        Dpmat[i,snp.list$snps[[i]]] <- 0
        pvmat[i,snp.list$snps[[i]]] <- 0
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
      # note that many sources give an inaccurate Dprime calculation--this should be correct.
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
      if(length(snp.list$snps[[i]]) == 1){
        Dprime <- Dprime[-2]
        rsq <- rsq[-2]
        pval <- pval[-2]
      }

      #write output.
      tprox <- cbind(meta[rep(i, length(Dprime)),],
                     meta[snp.list$snps[[i]],],
                     abs(meta[i,]$position - meta[snp.list$snps[[i]],]$position),
                     rsq, Dprime, pval)
      colnames(tprox) <- colnames(prox)
      prox <- rbind(prox, tprox)

      #reminder: columns start at locus two, rows start at locus one (but end at nlocus - 1)
      rmat[i,snp.list$snps[[i]]] <- rsq
      Dpmat[i,snp.list$snps[[i]]] <- Dprime
      pvmat[i,snp.list$snps[[i]]] <- pval
    }

    return(list(prox = prox, Dprime = Dpmat, rsq = rmat, pval = pvmat))

  }

  # function to figure out which snps we are comparing to each
  # outputs a nested list. Each entry in the list is a unique sample facet. In each of these lists is an entry for each unique subfacet level.
  # in this is an entry for each snp that lists the snps it is compared to.
  # If multiple entries would write to the same sample facet and subfacet, it will just add any new comparisons needed.
  determine.comparison.snps <- function(x, facets, facet.types){
    # progress: need to make this work with the .base facet

    # sub-subfunctions to get the options for the snp and sample facets
    get.samp.opts <- function(x, t.facet){
      sample.meta <- x@sample.meta[colnames(x@sample.meta) %in% t.facet]
      sample.meta <- sample.meta[,sort(colnames(sample.meta))]
      if(!is.data.frame(sample.meta)){
        sample.meta <- as.data.frame(sample.meta)
        colnames(sample.meta) <- colnames(x@sample.meta)[colnames(x@sample.meta) %in% t.facet]
      }
      sample.opts <- unique(sample.meta)
      if(!is.data.frame(sample.opts)){
        sample.opts <- as.data.frame(sample.opts, stringsAsFactors = F)
        colnames(sample.opts) <- facets[which(facets %in% colnames(x@sample.meta))]
      }
      sample.opts <- dplyr::arrange_all(sample.opts)

      return(list(sample.opts, sample.meta))
    }

    get.snp.opts <- function(x, t.facet){
      snp.meta <- x@snp.meta[colnames(x@snp.meta) %in% t.facet]
      snp.meta <- snp.meta[,sort(colnames(snp.meta))]
      if(!is.data.frame(snp.meta)){
        snp.meta <- as.data.frame(snp.meta)
        colnames(snp.meta) <- colnames(x@snp.meta)[colnames(x@snp.meta) %in% t.facet]
      }
      snp.opts <- unique(snp.meta)
      if(!is.data.frame(snp.opts)){
        snp.opts <- as.data.frame(snp.opts, stringsAsFactors = F)
        colnames(snp.opts) <- facets[which(facets %in% colnames(x@snp.meta))]
      }
      snp.opts <- dplyr::arrange_all(snp.opts)

      return(list(snp.opts, snp.meta))
    }

    # pull out just the sample facets
    sample.facets <- check.snpR.facet.request(x, facets) # the sample level facets that we are working with.

    # initialize the output list
    out <- vector(mode = "list", length(sample.facets))
    names(out) <- sample.facets
    if(any(facet.types == "snp")){
      out <- c(out, list(.base = list(.base = list(snps = vector("list", nrow(x)), samps = 1:nrow(x@sample.meta)))))
    }

    # loop through each facet, do different things depending on facet type
    for(i in 1:length(facets)){

      # grab the facet level list we are writing to.
      t.samp.facet <- check.snpR.facet.request(x, facets[i])
      write.facet <- which(names(out) == t.samp.facet)
      if(length(write.facet) == 0){
        t.samp.facet <- ".base"
        write.facet <- which(names(out) == ".base")
      }
      t.facet <- unlist(strsplit(facets[i], split = "(?<!^)\\.", perl = T))

      this.out <- out[[write.facet]]

      # for sample only, need to loop through each sample level subfacet, then loop through all snps
      if(facet.types[i] == c("sample")){

        # get options
        opts <- get.samp.opts(x, t.facet)
        sample.opts <- opts[[1]]
        sample.meta <- opts[[2]]

        if(t.facet == ".base"){
          sample.opts <- matrix(".base")
          sample.meta <- matrix(".base", nrow = nrow(x@sample.meta))
        }

        # add snp/snp comparisons. Since the facet is simple, we do all snps. Do so with a loop through all subfacets
        if(is.null(this.out)){
          this.out <- vector("list", nrow(sample.opts))
          names(this.out) <- do.call(paste, as.data.frame(sample.opts))
        }

        for(j in 1:nrow(sample.opts)){

          # grab the subfacet level we are writing to.
          write.subfacet <- which(names(this.out) == paste(sample.opts[j,], collapse = " "))
          this.subfacet <- this.out[[write.subfacet]]


          if(is.null(this.subfacet)){
            samps.in.subfacet <- which(apply(sample.meta, 1, function(x) identical(as.character(x), as.character(sample.opts[j,]))))
            this.subfacet <- list(snps = vector("list", nrow(x)), samps = samps.in.subfacet)
          }

          # add comparisons for each snp. Note that the last snp, with no comparisons to do, will recieve a NULL
          for(k in 1:(nrow(x) - 1)){
            c.comps <- this.subfacet$snps[[k]]
            c.comps <- c(c.comps, (k + 1):nrow(x))
            dups <- which(duplicated(c.comps))
            if(length(dups) > 0){
              this.subfacet$snps[[k]] <- c.comps[-dups]
            }
            else{
              this.subfacet$snps[[k]] <- c.comps
            }
          }

          # add back to this.out
          this.out[[write.subfacet]] <- this.subfacet
        }

      }

      # for snp only, need to loop through each snp level subfacet, then through all snps on that subfacet
      else if(facet.types[i] == "snp"){
        # get the subfacet options
        opts <- get.snp.opts(x, t.facet)
        snp.opts <- opts[[1]]
        snp.meta <- opts[[2]]

        # add snp/snp comparisons. Since the facet is simple, we do all samples, but pick the correct snps. This will be at the .base facet and .base subfacet!
        for(j in 1:nrow(snp.opts)){
          snps.in.subfacet <- which(apply(snp.meta, 1, function(x) identical(as.character(x), as.character(snp.opts[j,]))))

          # add comparisons for each snp. Note that the last snp, with no comparisons to do, will recieve a NULL
          for(k in 1:(length(snps.in.subfacet) - 1)){
            c.comps <- this.out$.base$snps[[snps.in.subfacet[k]]]
            c.comps <- c(c.comps, snps.in.subfacet[(k + 1):length(snps.in.subfacet)])
            dups <- which(duplicated(c.comps))
            if(length(dups) > 0){
              this.out$.base$snps[[snps.in.subfacet[k]]] <- c.comps[-dups]
            }
            else{
              this.out$.base$snps[[snps.in.subfacet[k]]] <- c.comps
            }
          }
        }
      }

      # for complex, need to loop through first each sample level subfacet, then through the snp level subfacet, then through each snp on that subfacet.
      else if(facet.types[i] == "complex"){

        # get the subfacet sample options, snp and sample
        sample.opts <- get.samp.opts(x, t.facet)
        snp.opts <- get.snp.opts(x, t.facet)
        sample.meta <- sample.opts[[2]]
        sample.opts <- sample.opts[[1]]
        snp.meta <- snp.opts[[2]]
        snp.opts <- snp.opts[[1]]


        if(is.null(this.out)){
          this.out <- vector("list", nrow(sample.opts))
          names(this.out) <- do.call(paste, as.data.frame(sample.opts))
        }


        # for each sample level option, we make sure that we compare only within snp facet level
        for(j in 1:nrow(sample.opts)){

          # grab the subfacet level we are writing to.
          write.subfacet <-which(names(this.out) == paste(sample.opts[j,], collapse = " "))
          this.subfacet <- this.out[[write.subfacet]]

          if(is.null(this.subfacet)){
            samps.in.subfacet <- which(apply(sample.meta, 1, function(x) identical(as.character(x), as.character(sample.opts[j,]))))
            this.subfacet <- list(snps = vector("list", nrow(x)), samps = samps.in.subfacet)
          }

          for(l in 1:nrow(snp.opts)){
            snps.in.subfacet <- which(apply(snp.meta, 1, function(x) identical(as.character(x), as.character(snp.opts[l,]))))

            # add comparisons for each snp. Note that the last snp, with no comparisons to do, will recieve a NULL
            if(length(snps.in.subfacet) == 1){next} # if only one snp here, no LD to calculate
            for(k in 1:(length(snps.in.subfacet) - 1)){
              cat(i, " ", j, " ", l, " ", k, "\n")
              c.comps <- this.subfacet$snps[[snps.in.subfacet[k]]]
              c.comps <- c(c.comps, snps.in.subfacet[(k + 1):length(snps.in.subfacet)])
              dups <- which(duplicated(c.comps))
              if(length(dups) > 0){
                this.subfacet$snps[[snps.in.subfacet[k]]] <- c.comps[-dups]
              }
              else{
                this.subfacet$snps[[snps.in.subfacet[k]]] <- c.comps
              }
            }
          }

          # add back to this.out
          this.out[[write.subfacet]] <- this.subfacet
        }
      }

      # save the output for this subfacet.
      out[[write.facet]] <- this.out
    }

    # return
    return(out)
  }


  # function to unpack a nested output list to parse out for snp level facets.
  decompose.LD.matrix <- function(x, LD_matrix, facets, facet.types){
    # for each facet type that included a snp.level facet, we need to split corrected matrices for sample or .base, just spit out everything
    out <- list()
    for(i in 1:length(facets)){
      # if a sample of base facet, just return it, no need for changes
      if(facet.types[i] == "sample" | facet.types[i] == ".base"){
        out[[i]] <- LD_matrix[[which(names(LD_matrix) == facets[i])]]
        names(out)[i] <- facets[i]
      }

      # otherwise need to split matrices
      else{

        # determine sample and snp parts
        samp.facet <- check.snpR.facet.request(x, facets[i])
        if(is.null(samp.facet)){samp.facet <- ".base"}
        snp.facet <- check.snpR.facet.request(x, facets[i], remove.type = "sample")
        split.snp.facet <- unlist(strsplit(snp.facet, "\\."))

        # grab the matrices for the corresponding sample level facet
        this.matrix <- LD_matrix[[which(names(LD_matrix) == samp.facet)]]

        # grab metadata and metadata options, ensuring correct column order
        this.meta <- x@snp.meta[,which(colnames(x@snp.meta) %in% split.snp.facet)]
        this.meta <- as.matrix(this.meta)
        if(ncol(this.meta) == 1){
          colnames(this.meta) <- split.snp.facet
        }
        else{
          this.meta <- this.meta[,order(colnames(this.meta))]
        }
        meta.opts <- as.matrix(unique(this.meta))
        colnames(meta.opts) <- colnames(this.meta)

        # intialize output
        out[[i]] <- vector(mode = "list", length = length(this.matrix))
        names(out[[i]]) <- names(this.matrix)
        names(out)[i] <- facets[i]

        for(k in 1:length(this.matrix)){
          # intialize and name
          out[[i]][[k]] <- vector(mode = "list", length = nrow(meta.opts))
          names(out[[i]][[k]]) <- do.call(paste, as.data.frame(meta.opts))
        }

        # loop through each meta option and subset parts of the matrix
        for(j in 1:nrow(meta.opts)){
          # which snps do we keep?
          keep.snps <- which(apply(this.meta, 1, function(x) identical(as.character(x), as.character(meta.opts[j,]))))

          # subset matrices
          for(k in 1:length(this.matrix)){
            Dprime <- this.matrix[[k]]$Dprime[keep.snps, keep.snps]
            rsq <- this.matrix[[k]]$rsq[keep.snps, keep.snps]
            pval <- this.matrix[[k]]$pval[keep.snps, keep.snps]

            # add to output
            out[[i]][[k]][[j]] <- list(Dprime = Dprime, rsq = rsq, pval = pval)
          }
        }
      }
    }

    return(out)
  }




  #========================primary looping function==========================
  # this will determine how to call the LD_func.
  # if just one level (".basic"), call the function simply, possibly par.
  # if multiple, take the output of determine.comparison.snps and loop through each subfacet level, doing the comps included.

  # the overall function. x is snpRdata object.
  func <- function(x, facets, snp.facets, par, sr){

    facet.types <- facets[[2]]
    facets <- facets[[1]]

    #=====================call functions=========
    # call these functions (possibly in parallel) according to supplied levels.

    cat("Beginning LD calculation...\n")

    #=====================no facets==============
    if( (length(facets) == 1 & facets[1] == ".base") | all(facet.types == "snp")){
      if(length(facets) == 1 & facets[1] == ".base"){
        comps <- determine.comparison.snps(x, facets, "sample")
      }

      else{
        comps <- determine.comparison.snps(x, facets, facet.types)
      }


      cat("No facets specified.\n")

      # grab metadata, mDat
      meta <- x@snp.meta
      mDat <- x@mDat


      # run in parallel if requested
      if(is.numeric(par)){
        cat("Running in parallel.\n\t")

        # each thread needs to be given a roughly equal number of comparisons to do
        ncomps <-  length(unlist(comps[[1]][[1]]$snps)) # number of comparisons
        split <- (ncomps)/par #number of comparisons to do per core
        split <- ceiling(split)
        cat("At least", split, "pairwise comparisons per processor.\n")

        # need to figure out which comps entries to null out for each processor.
        comps.per.snp <- unlist(lapply(comps[[1]][[1]]$snps, length))
        rproc <- ceiling(cumsum(as.numeric(comps.per.snp))/split) # which processor should each comparison be assigned to?

        #now need to start the parallel job:
        library(doParallel);library(foreach)
        cl <- snow::makeSOCKcluster(par)
        doSNOW::registerDoSNOW(cl)

        #prepare reporting function
        ntasks <- par
        progress <- function(n) cat(sprintf("Part %d out of",n), ntasks, "is complete.\n")
        opts <- list(progress=progress)

        # initialize and store things
        x_storage <- as.matrix(as.data.frame(x))
        na.test <- suppressWarnings(as.numeric(x_storage[1]))
        #save the info as a bigmatrix if it can be safelx_storage converted to numeric. Usuallx_storage this is true for ms but not necissarilx_storage other data tx_storagepes.
        if(!is.na(na.test)){
          if(as.numeric(x_storage[1]) == x_storage[1]){
            cat("Saving matrix as big.matrix object for quicker sharing.\n")
            xb <- bigmemorx_storage::as.big.matrix(x, tx_storagepe = "char")
            xbd <- bigmemorx_storage::describe(xb)
            remove(x)
          }
        }
        meta_storage <- x@snp.meta
        mDat_storage <- x@mDat
        t.comps <- comps[[1]][[1]]$snps

        cat("Begining run.\n")

        # run the LD calculations
        output <- foreach::foreach(q = 1:ntasks, .packages = c("bigmemory", "dplyr"), .inorder = TRUE,
                                   .options.snow = opts, .export = c("LD_func", "tabulate_haplotypes", "GtoH")) %dopar% {
                                     if(exists("xbd")){
                                       x_storage <- bigmemory::attach.big.matrix(xbd)
                                     }

                                     # get comps and run
                                     t.comps[which(rproc != q)] <- vector("list", sum(rproc != q)) # null out any comparisons we aren't doing
                                     t.comps <- list(snps = t.comps, samps = comps[[1]][[1]]$samps)

                                     LD_func(x = x_storage, snp.list = t.comps,
                                             meta = meta_storage, mDat = mDat_storage,
                                             sr = T)
                                     }

        #release cores
        parallel::stopCluster(cl)
        doSNOW::registerDoSNOW()


        cat("LD computation completed. Preparing results.\n\t")

        # combine results
        ## initialize
        prox <- data.frame()
        Dprime <- matrix(NA, nrow(x), nrow(x))
        colnames(Dprime) <- x@snp.meta$position
        row.names(Dprime) <- x@snp.meta$position
        rsq <- Dprime
        pval <- Dprime

        ## combine results
        for(i in 1:length(output)){
          prox <- rbind(prox, output[[i]]$prox)

          # just overwrite anything that is NA. If it was NA because of poor data at a legitimate comparison, it will just get overwritten with NA. Nothing with data should be overwritten like this.
          fill <- which(is.na(Dprime))
          Dprime[fill] <- output[[i]]$Dprime[fill]
          rsq[fill] <- output[[i]]$rsq[fill]
          pval[fill] <- output[[i]]$pval[fill]
        }

        # decompose and return (mostly for snp level facets)
        ## prep for decomposition function, done to make the format equal to something with sample level facets.
        LD_mats <- vector("list", 1)
        LD_mats[[1]] <- vector("list", 1)
        names(LD_mats) <- ".base"
        LD_mats[[1]][[1]] <- vector("list", 1)
        names(LD_mats[[1]]) <- ".base"
        LD_mats[[1]][[1]] <- list(Dprime = Dprime, rsq = rsq, pval = pval)

        ## decompose and return
        LD_mats <- decompose.LD.matrix(x, LD_mats, facets = facets, facet.types = facet.types)
        out <- list(prox = prox, LD_mats = LD_matrices)
        return(out)
      }

      #otherwise run normally
      else{
        browser()
        out <- LD_func(x, meta, snp.list = comps[[1]][[1]], mDat = mDat, sr)

        # decompose and return (mostly for snp level facets)
        ## prep for decomposition function, done to make the format equal to something with sample level facets.
        LD_mats <- vector("list", 1)
        LD_mats[[1]] <- vector("list", 1)
        names(LD_mats) <- ".base"
        LD_mats[[1]][[1]] <- vector("list", 1)
        names(LD_mats[[1]]) <- ".base"
        LD_mats[[1]][[1]] <- list(Dprime = out$Dprime, rsq = out$rsq, pval = out$pval)

        prox <- out$prox
        out <- decompose.LD.matrix(x, LD_mats, facets, facet.types)
        out <- list(prox = prox, LD_matrices = out)

        return(out)
      }
    }

    #=====================facets=================
    # approach/psuedo-code:
    # Each facet is a level to break down by. "pop" means break by pop, c("pop", "group") means break twice, once by pop, once by group,
    # c("pop.group") means to break by pop and group.
    # For each sample level facet, we must loop through all snps, since D values will be different depending on what samples we look at.
    # These must be looped through seperately!
    # If there are multiple snp level facets requested, there is no reason to do re-do snp/snp comparisons within each sample level facet. Just do the all of the relevent snp/snp comparisons.
    # If there are complex facets with repeated sample levels (c("pop.group", "pop.subgroup")), then same deal.

    # So:
    #     For each sample level facet:
    #       check if we've run the facet before
    #       For each level of those facets:
    #         For each snp:
    #           Figure out which snps we need to compare to across all snp level facets.
    #           Remove any comparisons that we've already done!
    #           Pass the genotypes and per snp comparison info to the LD_func (need to edit that function slightly to accomodate)
    #     Parse and output results.

    # as a part of this, need a function to determine the comparisons to do for each facet and subfacet.

    comps <- determine.comparison.snps(x, facets, facet.types)

    #prepare output list
    w_list<- vector("list", length = length(comps))
    names(w_list) <- names(comps)
    tot_subfacets <- 0
    task_list <- matrix(NA, 0, 2)
    for(i in 1:length(comps)){
      w_list[[i]] <- vector("list", length = length(comps[[i]]))
      names(w_list[[i]]) <- names(comps[[i]])
      tot_subfacets <- tot_subfacets + length(comps[[i]])
      for(j in 1:length(w_list[[i]])){
        w_list[[i]][[j]] <- list(Dprime = NULL, rsq = NULL, pvalue = NULL)
        task_list <- rbind(task_list, c(i, j))
      }
    }
    w_list <- list(prox = NULL, LD_mats = w_list)

    #not in parallel
    if(par == FALSE){

      #loop through each set of facets
      progress <- 1
      for (i in 1:length(comps)){
        for(j in 1:length(comps[[i]])){
          cat("Subfacet #:", progress, "of", tot_subfacets, " Name:", paste0(names(comps)[i], " " , names(comps[[i]])[j]), "\n")
          out <- LD_func(x, meta = x@snp.meta, mDat = x@mDat, snp.list = comps[[i]][[j]], sr = sr)
          #report progress
          progress <- progress + 1
          w_list$prox <- rbind(w_list$prox, cbind(out$prox, sample.facet = names(comps)[i]))
          w_list$LD_mats[[i]][[j]][[1]] <- out$Dprime
          w_list$LD_mats[[i]][[j]][[2]] <- out$rsq
          w_list$LD_mats[[i]][[j]][[3]] <- out$pval

        }
      }

      # split apart matrices and decompose
      prox <- w_list$prox
      mats <- decompose.LD.matrix(x, w_list$LD_mats, facets, facet.types)
      w_list <- list(prox = prox, LD_matrices = mats)

      return(w_list)
    }

    #in parallel
    else{
      library(doParallel);library(foreach)
      cl <- snow::makeSOCKcluster(par)
      doSNOW::registerDoSNOW(cl)

      #prepare reporting function
      ntasks <- tot_subfacets
      progress <- function(n) cat(sprintf("Facet %d out of", n), ntasks, "is complete.\n")
      opts <- list(progress=progress)

      x_storage <- as.matrix(as.data.frame(x))
      meta_storage <- x@snp.meta
      mDat_storage <- x@mDat


      #loop through each set of facets
      output <- foreach::foreach(i = 1:ntasks, .packages = c("dplyr", "reshape2", "matrixStats", "bigtabulate", "snpR"), .inorder = TRUE,
                                 .options.snow = opts, .export = c("LD_func", "tabulate_haplotypes", "GtoH")) %dopar% {
                                   t.task <- task_list[i,]
                                   t.facet <- t.task[1]
                                   t.subfacet <- t.task[2]
                                   LD_func(x_storage, meta = meta_storage, mDat = mDat_storage, snp.list = comps[[t.facet]][[t.subfacet]], sr = sr)
                                 }

      #release cores and clean up
      parallel::stopCluster(cl)
      doSNOW::registerDoSNOW()
      rm(x_storage, meta_storage, mDat_storage)
      gc();gc()


      # make the output sensible but putting it into the same format as from the other, then running the decompose function
      for(i in 1:ntasks){
        t.facet <- task_list[i,1]
        t.subfacet <- task_list[i,2]
        w_list$prox <- rbind(w_list$prox, cbind(output[[i]]$prox, sample.facet = names(comps)[t.facet]))
        w_list$LD_mats[[t.facet]][[t.subfacet]]$Dprime <- output[[i]]$Dprime
        w_list$LD_mats[[t.facet]][[t.subfacet]]$rsq <- output[[i]]$rsq
        w_list$LD_mats[[t.facet]][[t.subfacet]]$pvalue <- output[[i]]$pval
      }

      # split apart matrices and decompose
      prox <- w_list$prox
      mats <- decompose.LD.matrix(x, w_list$LD_mats, facets, facet.types)
      w_list <- list(prox = prox, LD_matrices = mats)

      #return
      return(w_list)
    }
  }

  #========================prepare and pass the primary function to apply.snpR.facets==================
  #subset data if requested:
  if(!(is.null(subfacets[1]))){
    ssfacets <- names(subfacets)

    # check complex facets
    complex.sfacets <- check.snpR.facet.request(x, ssfacets, remove.type = "simple")
    if(length(complex.sfacets) > 0){
      stop("Complex (snp and sample) level subfacets not accepted. Providing these as seperate snp and sample subfacets will run only snps/samples contained in both levels.\n")
    }

    # combine duplicates
    if(any(duplicated(ssfacets))){
      dup.sfacets <- subfacets[which(duplicated(ssfacets))]
      subfacets <- subfacets[-which(duplicated(ssfacets))]
      for(i in 1:length(dup.sfacets)){
        wmatch <- which(names(subfacets) == names(dup.sfacets[1]))
        subfacets[[wmatch]] <- c(subfacets[[wmatch]], dup.sfacets[[i]])
        ndup <- which(duplicated(subfacets[[wmatch]]))
        if(length(ndup) > 0){
          subfacets[[wmatch]] <- subfacets[[wmatch]][-ndup]
        }
      }
      ssfacets <- ssfacets[-which(duplicated(ssfacets))]
    }

    # get subfacet types
    ssfacet.types <- check.snpR.facet.request(x, ssfacets, "none", T)[[2]]

    invisible(capture.output(x <- subset.snpR.data(x,
                                                   facets = names(subfacets)[which(ssfacet.types == "sample")],
                                                   subfacets = subfacets[[which(ssfacet.types == "sample")]],
                                                   snp.facets = names(subfacets)[which(ssfacet.types == "snp")],
                                                   snp.subfacets = subfacets[[which(ssfacet.types == "snp")]])))
  }

  if(is.numeric(ss)){
    #get sample
    if(ss <= 1){
      ss <- sample(nrow(x), round(nrow(x)*ss), F)
    }
    else{
      ss <- sample(nrow(x), ss, F)
    }

    #subset
    x <- subset.snpR.data(x, ss)
  }


  # typical facet check, keeping all facet types but removing duplicates. Also returns the facet type for later use.
  facets <- check.snpR.facet.request(x, facets, remove.type = "none", return.type = T)

  # run the function
  out <- func(x, facets = facets, snp.facets = snp.facets, par = par, sr = sr)

  # add to snpRdata object and return
  out <- merge.snpR.stats(x, out, "LD")
  return(out)
}




#'Calculate HWE divergence from SNP data.
#'
#'\code{calc_hwe} Calculates p-values for Hardy-Weinburg Equilibrium divergence from snp.R.data objects.
#'
#' @param x Input SNP data, in the snpRdata format.
#' @param facets Facets to run. Non-sample level facets will be removed. Periods can be used to list complex facets, such as pop.fam to split by both pop and fam simultaniously.
#'
#' @return A snpRdata object with HWE data added.
#'
#' @examples
#'
calc_hwe <- function(x, facets = NULL, method = "exact"){
  func <- function(gs, method){
    # exact test to use if there are too few observations in a cell
    # edited from Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of
    # Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000
    # code available at http://csg.sph.umich.edu/abecasis/Exact/snp_hwe.r
    exact.hwe <- function(pp, qq, pq2){
      if(all(c(pp, qq, pq2) == 0)){
        return(-1.0)
      }
      obs_homr <- min(c(pp, qq))
      obs_homc <- max(c(pp, qq))
      obs_hets <- pq2
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
      if ((mid %% 2) != (rare %% 2)) mid <- mid + 1

      probs[mid + 1] <- 1.0
      mysum <- 1.0

      # Calculate probablities from midpoint down
      curr_hets <- mid
      curr_homr <- (rare - mid) / 2
      curr_homc <- N - curr_hets - curr_homr

      while (curr_hets >=  2){
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

      while (curr_hets <= rare - 2){
        probs[curr_hets + 3] <- probs[curr_hets + 1] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
        mysum <- mysum + probs[curr_hets + 3]

        # add 2 heterozygotes -> subtract 1 rare homozygtote, 1 common homozygote
        curr_hets <- curr_hets + 2
        curr_homr <- curr_homr - 1
        curr_homc <- curr_homc - 1
      }

      # P-value calculation
      target <- probs[obs_hets + 1]
      p <- min(1.0, sum(probs[probs <= target])/ mysum)
      return(p)
    }

    gs <- gs$gs

    # get observed genotype counts
    het.col <- which(substr(colnames(gs), 1, 1) != substr(colnames(gs), 2, 2))
    o2pq <- rowSums(gs[,het.col])
    opp <- matrixStats::rowMaxs(gs[,-het.col])
    oqq <- rowSums(gs) - o2pq - opp

    # if we are using a chisq test, easy and quick
    if(method == "chisq"){
      # get allele frequencies
      fp <- (opp*2 + o2pq)/(rowSums(gs)*2)
      fq <- 1 - fp

      # get expected genotype counts
      epp <- fp^2 * rowSums(gs)
      eqq <- fq^2 * rowSums(gs)
      e2pq <- 2*fp*fq * rowSums(gs)

      # calculate chi-squared
      calc.chi <- function(o,e){
        return(((o-e)^2)/e)
      }
      chi.pp <- calc.chi(opp, epp)
      chi.qq <- calc.chi(oqq, eqq)
      chi.2pq <- calc.chi(o2pq, e2pq)
      chi <- chi.pp + chi.qq + chi.2pq

      # calculate p-values
      out <- pchisq(chi, 2, lower.tail = FALSE)
      return(out)
    }

    # otherwise we have to use the looped version:
    else if(method == "exact"){
      cat("Using exact test from Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005).\n")
      out <- numeric(nrow(gs))
      for(i in 1:nrow(gs)){
        out[i] <- exact.hwe(oqq[i], opp[i], o2pq[i])
      }
      return(out)
    }
  }

  if(!(method %in% c("exact", "chisq"))){stop("Unrecognized HWE method, please use chisq or exact.\n")}

  # add any missing facets
  facets <- check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    invisible(capture.output(x <- add.facets.snpR.data(x, facets)))
  }

  out <- apply.snpR.facets(x, facets, "gs", func, case = "ps", method = method)
  colnames(out)[ncol(out)] <- "pHWE"
  x <- merge.snpR.stats(x, out)

  return(x)
}
