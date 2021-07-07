#'Calculate standard single-SNP genetic statistics
#'
#'These functions calculate a basic genetic statistics from SNP data contained
#'in snpRdata objects, splitting samples by any number of provided facets. Since
#'these statistics all use only a single SNP at a time, they ignore any SNP
#'specific facet levels.
#'
#'The data can be broken up categorically by sample metadata, as described in
#'\code{\link{Facets_in_snpR}}.
#'
#'@section pi:
#'
#'  Calculates pi (genetic diversity/average number of pairwise differences)
#'  according to Hohenlohe et al. (2010).
#'
#'@section maf:
#'
#'  Calculates minor allele frequencies and note identities and counts of major
#'  and minor alleles.
#'
#'@section ho:
#'
#'  Calculates observed heterozygosity.
#'
#'@section private alleles:
#'
#'  Determines if each SNP is a private allele across all levels in each sample
#'  facet. Will return an error of no sample  facets are provided.
#'
#'@section HWE:
#'
#'  Calculates a p-value for the null hypothesis that a population is in HWE at
#'  a given locus. Several methods available: \itemize{ \item{"exact"} Exact
#'  test according to Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005).
#'  Slightly slower. \item{"chisq"} Chi-squared test. May produce poor results
#'  when sample sizes for any observed or expected genotypes are small. }
#'
#'  For the exact test, code derived from
#'  \url{http://csg.sph.umich.edu/abecasis/Exact/snp_hwe.r}
#'
#'@param x snpRdata. Input SNP data.
#'@param facets character. Categorical metadata variables by which to break up
#'  analysis. See \code{\link{Facets_in_snpR}} for more details.
#'@param method character, default "exact". Defines the method to use for
#'  calculating p-values for HWE divergence. Options: \itemize{ \item{exact: }
#'  Uses the exact test as described in Wigginton et al (2005). \item{chisq: }
#'  Uses a chi-squared test. } See details
#'@param fwe_method character, default c("bonferroni", "holm", "BH", "BY"). Type
#'  of Family-Wise Error correction (mulitple testing correction) to use. For
#'  details and options, see \code{\link{p.adjust}}.
#'@param fwe_case character, default c("by_facet", "by_subfacet", "overall").
#'  How should Family-Wise Error correction (multiple testing correction) be
#'  applied? \itemize{\item{"by_facet":} Each facet supplied (such as pop or
#'  pop.fam) is treated as a set of tests. \item{"by_subfacet":} Each level of
#'  each subfacet is treated as a seperate set of tests. \item{"overall":} All
#'  tests are treated as a set.}
#'
#'@aliases calc_pi calc_hwe calc_ho calc_private calc_maf
#'
#'@return snpRdata object with requested stats merged into the stats socket
#'
#'@name calc_single_stats
#'
#'
#'@author William Hemstrom
#'
#'@references Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005). \emph{American
#'  Journal of Human Genetics}
#'@references Hohenlohe et al. (2010). \emph{PLOS Genetics}
#'
#' @examples
#' # base facet
#' x <- calc_pi(stickSNPs)
#' get.snpR.stats(x)
#'
#' # multiple facets
#' x <- calc_pi(stickSNPs, facets = c("pop", "pop.fam"))
#' get.snpR.stats(x, c("pop", "pop.fam"))
#'
#' # HWE with family-wise error correction
#' x <- calc_hwe(stickSNPs, facets = c("pop", "pop.fam"), method = "exact")
#' get.snpR.stats(x, c("pop", "pop.fam"))
NULL


#'Facets in snpR
#'
#'Facets are used to describe the ways in which data should be broken down/split
#'up for analysis.
#'
#'Facet designation follows a specific format. Facets are given as a character
#'vector, where each entry designates one unique combination of levels over
#'which to seperate the data. Levels within a facet are seperated by a '.'. Each
#'facet can contain multiple snp and/or sample levels. Multiple facets can be
#'run with a single line of code.
#'
#'For example, c("group.pop", "pop") would split the data first by both group
#'and pop and then by pop alone. This will produce the same result as running
#'the function with both "group.pop" and then again with "pop", although the
#'former is typically more computationally efficient.
#'
#'If multiple sample or snp levels are provided in a single facet, the data is
#'simultaniously broken up by \emph{both} levels. For example, the facet
#'c("fam.pop") would break up the data provided in \code{\link{stickSNPs}} by
#'both family and population, and would produce levels such as "ASP.A" for
#'individuals in the ASP population and in family A.
#'
#'All listed facet levels must match column names in the snp of sample level
#'facet data provided when constructing a snpRdata object with
#'\code{\link{import.snpR.data}}. The order of the provided levels within each
#'facet does not matter.
#'
#'In many cases, specifying "all" as a facet will calculate or return statistics
#'for all previously run facets.
#'
#'The base facet--that is, the entire data with no categorical devisions--can be
#'specified with ".base" and is typically the facet defaulted to when facets =
#'NULL.
#'
#'
#'See examples.
#'
#'@name Facets_in_snpR
#'
#' @examples
#' # base facet
#' x <- calc_pi(stickSNPs)
#' get.snpR.stats(x)
#'
#' # multiple facets
#' x <- calc_pi(stickSNPs, facets = c("pop", "pop.fam"))
#' get.snpR.stats(x, c("pop", "pop.fam"))
#'
NULL

#'@export
#'@describeIn calc_single_stats pi (average number of pairwise differences/expected heterozygosity)
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
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  

  # add any missing facets
  ofacets <- facets
  facets <- check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    invisible(utils::capture.output(x <- add.facets.snpR.data(x, facets)))
  }

  out <- apply.snpR.facets(x, facets, "ac", func, case = "ps")
  colnames(out)[ncol(out)] <- "pi"
  x <- merge.snpR.stats(x, out)
  x <- calc_weighted_stats(x, ofacets, type = "single", "pi")
  x <- update_calced_stats(x, facets, "pi", "snp")

  return(x)
}


#'@export
#'@describeIn calc_single_stats minor allele frequency
calc_maf <- function(x, facets = NULL){
  maj_count <- NULL
  
  
  # function to run on whatever desired facets:
  func <- function(gs, m.al, ref = NULL){
    
    # for the base facet, determine the major and minor then calculate maf
    if(is.null(ref)){
      # major alleles via max.col
      fmax <- colnames(gs$as)[max.col(gs$as, ties.method = "last")]
      lmax <- colnames(gs$as)[max.col(gs$as, ties.method = "first")]
      
      # minor alleles, via inversion, 0 replacement with -Inf, and max.col
      inv.as <- gs$as * -1
      inv.as[inv.as == 0] <- -Inf
      fmin <- colnames(gs$as)[max.col(inv.as, ties.method = "last")]
      lmin <- colnames(gs$as)[max.col(inv.as, ties.method = "first")]
      
      # special cases
      match.freq <- which(fmax != lmax) # maf = 0.5
      unseq <- which(matrixStats::rowSums2(gs$as) == 0) # unsequenced
      np <- which(matrixStats::rowSums2(matrix(as.logical(gs$as), nrow(gs$as))) == 1) # non-polymorphic
      
      # declair major and minor
      major <- fmax
      minor <- fmin
      ## maf = 0.05
      if(length(match.freq) != 0){
        minor[match.freq] <- lmax[match.freq]
      }
      ## unsequenced
      if(length(unseq) != 0){
        major[unseq] <- "N"
        minor[unseq] <- "N"
      }
      ## non-polymorphic
      if(length(np) != 0){
        minor[np] <- "N"
      }
      
      # grab the actual maf
      maf <- 1 - matrixStats::rowMaxs(gs$as)/matrixStats::rowSums2(gs$as)
      maf[is.nan(maf)] <- 0
      
      # get the major and minor counts
      # round because repeating decimals will yeild things like 1.00000000001 instead of 1. Otherwise this approach is quick and easy, as long as things are bi-allelic (non-polymorphic and equal min maj frequencies are fine.)
      maj.count <- round(rowSums(gs$as)*(1-maf))
      min.count <- round(rowSums(gs$as)*(maf))
    }
    
    # for non-base facets, use the given major and minor to calculate maf
    else{

      # use a data.table function to get the major allele counts
      adt <- as.data.table(gs$as)
      rep.factor <- nrow(gs$as)/nrow(ref)
      major <- rep(ref$major, each = rep.factor) # rep for each facet level, since that's how they are sorted
      adt$major <-  major
      adt <- adt[, maj_count := .SD[[.BY[[1]]]], by=major] # get the counts of the major allele in each column
      ac.rows <- 1:(which(colnames(adt) == "major") - 1)
      total.allele.count <- rowSums(adt[,ac.rows, with = FALSE])
      maf <- 1 - adt$maj_count/total.allele.count # get the maf
      
      # other things for return
      minor <- rep(ref$minor, each = rep.factor)
      maj.count <- adt$maj_count
      min.count <- total.allele.count - maj.count
    }
    
    # return
    return(data.table(major = major, minor = minor, maj.count = maj.count, min.count = min.count, maf = maf, stringsAsFactors = F))
  }
  
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  # add any missing facets
  facets <- check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    invisible(utils::capture.output(x <- add.facets.snpR.data(x, facets)))
  }
  
  if(facets[1] == ".base" & length(facets) == 1){
    out <- apply.snpR.facets(x,
                             facets = facets[1],
                             req = "gs",
                             fun = func,
                             case = "ps",
                             m.al = substr(x@mDat,1, nchar(x@mDat)/2))
  }
  else{
    # calculate the base facet if not yet added
    logi <- check_calced_stats(x, ".base", "maf")
    if(!logi[[".base"]]["maf"]){
      x <- calc_maf(x)
    }
    
    major_minor_base <- .get.snpR.stats(x)[,c("major", "minor")]
    out <- apply.snpR.facets(x,
                             facets = facets,
                             req = "gs",
                             fun = func,
                             case = "ps",
                             m.al = substr(x@mDat,1, nchar(x@mDat)/2),
                             ref = major_minor_base)
    
  }
  
  x <- update_calced_stats(x, facets, "maf", "snp")
  return(merge.snpR.stats(x, out))
}

#'Tajima's D from SNP data.
#'
#'\code{Tajimas_D} calculates Tajima's theta/pi, Waterson's theta, and Tajima's
#'D over a sliding window.
#'
#'Tajima's D compares estimates of theta based on either the number of observed
#'pairwise differences (Tajima's theta) and the number of substitutions vs
#'expected total tree length (Waterson's Theta). Since low frequency minor
#'variants contribute to these statistics and they rely on the ratio of the
#'number of variants vs the number of sequenced non-polymorphic sites, this
#'function should only be run on data that is \emph{unfiltered} aside from the
#'removal of poorly sequenced bases, ect.
#'
#'The data can be broken up categorically by either SNP and/or sample metadata,
#'as described in \code{\link{Facets_in_snpR}}.
#'
#'@param x snpRdata. Input SNP data.
#'@param facets character. Categorical metadata variables by which to break up
#'  analysis. See \code{\link{Facets_in_snpR}} for more details. If no snp level
#'  facets are provided, the calculated Tajima's D will be for the entire
#'  genome.
#'@param sigma numeric, default NULL. Sliding window size, in kilobases. Each
#'  window will include all SNPs within 3*sigma kilobases. If either sigma or
#'  step are NULL, the entire snp subfacet will be done at once (for example,
#'  the whole chromosome).
#'@param step numeric, default NULL. Number of bases to move between each
#'  window, in kilobases. If either sigma or step are NULL, the entire snp
#'  subfacet will be done at once (for example, the whole chromosome).
#'@param par numeric or FALSE, default FALSE. If numeric, the number of cores to
#'  use for parallel processing.
#'
#'@return snpRdata object, with Waterson's Theta, Tajima's Theta, and Tajima's D
#'  for each window merged in to the window.stats slot.
#'
#' @examples
#' # broken by population, windows across linkage group
#' x <- calc_tajimas_d(stickSNPs, facets = "group.pop", sigma = 200, step = 50)
#' get.snpR.stats(x, "group.pop", "single.window")
#'
#' # the entire population at once, note that sigma and step are NULL and no chromosome/linkage group/scaffold/etc set.
#' # this will calculate overall tajima's D without a window for each population.
#' x <- calc_tajimas_d(stickSNPs, facets = "pop")
#' get.snpR.stats(x, "pop", "single.window")
#'
#' # for the overall dataset, note that sigma and step are NULL
#' # this will calculate overall tajima's D for each group/pop
#' x <- calc_tajimas_d(stickSNPs, facets = "group.pop")
#' get.snpR.stats(x, "pop.group", "single.window")
#'
#'@export
#'@references Tajima, F. (1989). \emph{Genetics}
#'@author William Hemstrom
calc_tajimas_d <- function(x, facets = NULL, sigma = NULL, step = NULL, par = FALSE){
  if(!is.snpRdata(x)){
    stop("x must be a snpRdata object.")
  }
  
  if(!is.null(sigma) & !is.null(step)){
    sanity_check_window(x, sigma, step, stats.type = "single", nk = TRUE, facets = facets)
  }
  else{
    sigma <- NULL
    step <- NULL
  }
  if(is.null(facets[1]) | ".base" %in% facets | "sample" %in% check.snpR.facet.request(x, facets, "none", TRUE)){
    warning("Tajima's D has little meaning of snps on different chromosomes are considered together. Consider adding a snp level facet.")
  }


  func <- function(ac, par, sigma, step, report){
    out <- data.frame(position = numeric(0), sigma = numeric(0), ws.theta = numeric(0), ts.theta = numeric(0), D = numeric(0), n_snps = numeric(0)) #initialize output
    tps <- sort(ac$position) #get the site positions, sort
    lsp <- tps[length(tps)] #get the position of the last site to use as endpoint
    c <- 0 #set starting position
    i <- 1 #set starting iteration for writing output

    # if sigma is null, set to run everything in one window
    if(is.null(sigma)){
      sigma <- range(ac$position)
      c <- mean(sigma)
      step <- c + 1
      sigma <- sigma[2] - sigma[1]
    }
    else{
      sigma <- 1000*sigma
    }
    
    # run the loop
    while (c <= lsp){
      start <- c - sigma*3 #change window start
      end <- c + sigma*3 #change window end

      #take all the snps in the window, calculate T's theta, W's theta, and T's D
      wsnps <- ac[ac$position <= end & ac$position >= start,] #get only the sites in the window
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
      if("pop" %in% colnames(x)){
        out[i,"pop"] = x[1,"pop"] #if a pop column is in the input, add a pop column here.
      }
      out[i,"position"] <- c
      out[i,"ws.theta"] <- ws.theta
      out[i,"ts.theta"] <- ts.theta
      out[i,"D"] <- D
      out[i, "n_snps"] <- nrow(wsnps)
      c <- c + step*1000
      i <- i + 1
    }
    if(nrow(out) > 0){
      out$sigma <- sigma/1000
    }
    return(out)
  }

  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }

  # add any missing facets
  add.facets <- check.snpR.facet.request(x, facets)
  if(!all(add.facets %in% x@facets)){
    invisible(utils::capture.output(x <- add.facets.snpR.data(x, add.facets)))
  }
  facets <- check.snpR.facet.request(x, facets, remove.type = "none")

  out <- apply.snpR.facets(x,
                           facets = facets,
                           req = "meta.ac",
                           fun = func,
                           case = "ps.pf.psf",
                           par = par,
                           sigma = sigma,
                           step = step)

  x <- merge.snpR.stats(x, out, type = "window.stats")
  x <- update_calced_stats(x, facets, "tajimas_D")
  x <- calc_weighted_stats(x, facets, type = "single.window", "ts.theta")
  x <- calc_weighted_stats(x, facets, type = "single.window", "ws.theta")
  x <- calc_weighted_stats(x, facets, type = "single.window", "D")
  
  return(x)

}




#'Pairwise FST from SNP data.
#'
#'\code{calc_pairwise_fst} calculates pairwise FST for each SNP for each
#'possible pairwise combination of populations.
#'
#'Calculates FST according to either Wier and Cockerham 1984, Wier 1990,
#'Hohenlohe et al 2010, or using the \code{\link[genepop]{Fst}} function from
#'the genepop package (see references).
#'
#'If the genpop option is used, several intermediate files will be created in
#'the current working directory, which are not cleaned automatically.
#'\code{calc_pairwise_fst} will ask for permission to continue if these files
#'already exist.
#'
#'The Wier and Cockerham (1984), Wier (1990), and genepop methods tend to
#'produce very similar results. Generally, either of the two former options are
#'prefered for computational efficieny.
#'
#'The data can be broken up categorically by either SNP and/or sample metadata,
#'as described in \code{\link{Facets_in_snpR}}. Since this is a pairwise
#'statistic, at least a single sample level facet must be provided.
#'
#'Method Options: \itemize{ \item{"WC": }{Wier and Cockerham 1984.}
#'\item{"Wier": }{Wier 1990.} \item{"Genepop": }{As used in genepop, Rousset
#'2008.} \item{"Hohenlohe": }{Hohenlohe 2010.} }
#'
#'@param x snpRdata. Input SNP data.
#'@param facets character. Categorical metadata variables by which to break up
#'  analysis. See \code{\link{Facets_in_snpR}} for more details.
#'@param method character, default "WC". Defines the FST estimator to use.
#'  Options: \itemize{ \item{WC: } Wier and Cockerham (1984). \item{Wier: } Wier
#'  (1990) \item{Hohenlohe: } Hohenlohe et al (2010), identical to the STACKS
#'  package. \item{Genepop: } Rousset (2008), uses the genepop package. }
#'
#'@return A snpRdata object with pairwise FST as well as the number of total
#'  observations at each SNP in each comparison merged in to the pairwise.stats
#'  slot.
#'
#'@references Wier and Cockerham (1984). \emph{Evolution}
#'@references Wier (1990). Genetic data analysis. Sinauer,  Sunderland, MA
#'@references Hohenlohe et al. (2010). \emph{PLOS Genetics}
#'@references Rousset (2008). \emph{Molecular Ecology Resources}
#'
#'@author William Hemstrom
#'@export
#'
#' @examples
#' # Using Wier and Cockerham 1984's method
#' x <- calc_pairwise_fst(stickSNPs, "pop")
#' get.snpR.stats(x, "pop", "fst")
#'
#' # Using genepop, note that the overall value is part 2 of the returned list
#' x <- calc_pairwise_fst(stickSNPs, "pop", "genepop")
#' get.snpR.stats(x, "pop", "fst")
calc_pairwise_fst <- function(x, facets, method = "WC"){
  facet <- subfacet <- .snp.id <- NULL
  
  func <- function(x, method, facets = NULL){
    if(method != "genepop"){
      x <- data.table::as.data.table(x)
      data.table::setkey(x, subfacet, .snp.id)
      pops <- sort(unique(x$subfacet))
      npops <- length(pops)
      pnk.length <- (npops*(npops-1))/2 * (nrow(x)/npops)
    }
    else{
      fmeta <- data.table::setDT(x@facet.meta)
      pops <- sort(unique(fmeta[facet == facets, subfacet]))
      rm(fmeta)
      npops <- length(pops)
      pnk.length <- (npops*(npops-1))/2 * (nrow(x))
    }

    pnk <- data.table::data.table(comparison = character(pnk.length), ntotal = integer(pnk.length))

    #===============genepop======================
    if(method == "genepop"){
      g.filename <- paste0("genepop.", facets, ".txt")
      if(file.exists(g.filename)){
        file.remove(g.filename)
      }
      init.files <- list.files()
      invisible(utils::capture.output(format_snps(x, output = "genepop", facets = facets, outfile = g.filename)))
      genepop::Fst(g.filename, pairs = TRUE)




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
      prog <- 1L #column fill progress tracker
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
        prog <- prog + as.integer(ncol(tmat)) #increment column progress.
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
      if(ncol(out$loci) == 1){ # for cases where there are only two comparison groups
        out$loci <- data.frame(comparison = paste0(pnames[1], "~", pnames[2]), fst = out$loci[,1])
      }
      colnames(out$loci) <- c("comparison", "fst")

      # get nk values:
      n_tots <- data.table::data.table(pop = x@facet.meta$subfacet[x@facet.meta$facet == facets],
                                       .snp.id = x@facet.meta$.snp.id[x@facet.meta$facet == facets],
                                       nk = x@ac[x@facet.meta$facet == facets, "n_total"])
      n_tots <- data.table::dcast(n_tots, .snp.id ~ pop, value.var = "nk")
      n_tots <- n_tots[,-1]


      prog <- 1L
      for(i in 1:(ncol(n_tots) - 1)){
        for(j in (i+1):ncol(n_tots)){
          new.rows <- prog:(prog + nrow(n_tots) - 1)
          data.table::set(pnk, i = new.rows, 1L, paste0(colnames(n_tots)[i], "~", colnames(n_tots)[j]))
          suppressWarnings(data.table::set(pnk, i = new.rows, 2L, n_tots[,i, with = FALSE] + n_tots[,j, with = FALSE]))
          prog <- prog + as.integer(nrow(n_tots))

        }
      }

      out$loci$n_total <- pnk$ntotal
      if(length(out$overall) == 1){
        names(out$overall) <- paste0(pnames[1], "~", pnames[2])
      }

      # clean and return, we're done.
      cat("Finished.\n")
      fin.files <- list.files()
      file.remove(fin.files[which(!fin.files %in% init.files)])
      return(out)
    }

    #===============others=====================
    data.table::setkey(x, subfacet, .snp.id) # sort the data


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

    #loop through each comparison and calculate pairwise FST at each site
    c.col <- 1L
    prog <- 1L
    for (i in 1:(length(pops) - 1)){ #i is the first pop
      idat <- x[subfacet == pops[i]] # get data for first pop
      j <- i + 1 #initialize j as the next pop
      for (j in j:length(pops)){#j is pop being compared
        jdat <- x[subfacet == pops[j]] #get data for second pop

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
            data.table::set(out, j = c.col, value = Fst) # write fst
          }
          else{
            S1 <- ssq - (1/(nbar-1))*(inner1 - inner2 - inner3)
            S2i1 <- ((r*(nbar - nc))/nbar)*inner1
            S2i2 <- (1/nbar)*((nbar-1)+(r-1)*(nbar-nc))*ssq
            S2i3 <- ((nbar-nc)/(4*nc^2))*hbar
            S2 <- inner1 - (nbar/(r*(nbar-1)))*(S2i1 -S2i2 - S2i3)
            Fst <- S1/S2
            data.table::set(out, j = c.col, value = Fst) # write fst
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
          data.table::set(out, j = c.col, value = Fst) # write fst
        }

        else{
          stop("Please select a method of calculating FST.\nOptions:\n\tWC: Weir and Cockerham (1984).\n\tWier: Wier (1990).\n\tHohenlohe: Hohenlohe et al. (2010).")
        }

        # update pnk
        pnk <- data.table::set(pnk, prog:(prog + nrow(idat) - 1), 1L, paste0(idat$subfacet, "~", jdat$subfacet))
        pnk <- data.table::set(pnk, prog:(prog + nrow(idat) - 1), 2L, as.integer(idat$n_total + jdat$n_total))
        prog <- prog + as.integer(nrow(idat))


        c.col <- c.col + 1L #agument c.col
      }
    }

    # melt, cbind pnk
    suppressWarnings(out <- data.table::melt(out))
    colnames(out) <- c("comparison", "fst")
    out$nk <- pnk$ntotal

    # return
    return(out)
  }

  #============================sanity and facet checks========================
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  if(method == "genepop"){
    check.installed("genepop")
  }
  
  if(any(x@ac$n_alleles > 2)){
    vio <- which(x@ac$n_alleles[x@facet.meta$facet %in% facets] > 2)
    vio <- unique(x@facet.meta$.snp.id[x@facet.meta$facet %in% facets][vio])
    stop(cat("Some loci have more than two alleles. Violating loci:\n", paste0(vio, collapse = "\n")))
  }

  # add any missing facets
  ofacets <- facets
  facets <- check.snpR.facet.request(x, facets, return.type = T)
  if(all(facets[[2]] == ".base")){
    stop("At least one sample level facet is required for pairwise Fst esitmiation.")
  }
  facets <- facets[[1]]
  if(!all(facets %in% x@facets)){
    invisible(utils::capture.output(x <- add.facets.snpR.data(x, facets)))
  }
  

  method <- tolower(method)

  #=================================call apply.snpR.facets, slightly different for each method, since they require different stuff.===============
  if(method == "genepop"){
    out <- apply.snpR.facets(x, facets, req = "snpRdata", fun = func, case = "facet.pairwise", method = "genepop")
    ave.fst <- out[[2]]
    out <- out[[1]]
  }
  else if(method == "wc" | method == "wier"){
    need.ho <- which(!unlist(check_calced_stats(x, facets, "ho")))
    if(length(need.ho) > 0){
      x <- calc_ho(x, facets[need.ho])
    }
    out <- apply.snpR.facets(x, facets, req = c("ac.stats"), case = "facet.pairwise", fun = func, method = method)
  }
  else if(method == "hohenlohe"){
    need.pi <- which(!unlist(check_calced_stats(x, facets, "pi")))
    if(length(need.pi) > 0){
      x <- calc_pi(x, facets[need.pi])
    }
    out <- apply.snpR.facets(x, facets, req = c("ac.stats"), case = "facet.pairwise", fun = func, method = "hohenlohe")
  }

  #==================merge and return============================
  x <- merge.snpR.stats(x, out, type = "pairwise")
  if(method != "genepop"){
    ofacets <- check.snpR.facet.request(x, ofacets, "none", T)
    bad.ofacets <- which(ofacets[[2]] == "snp" | ofacets[[2]] == ".base")
    if(length(bad.ofacets) > 0){
      ofacets <- ofacets[[1]][-bad.ofacets]
    }
    else{
      ofacets <- ofacets[[1]]
    }
    x <- calc_weighted_stats(x, ofacets, type = "pairwise", stats_to_get = "fst")
  }
  else{
    ave.fst$snp.facet <- ".base"
    ave.fst$snp.subfacet <- ".base"
    ave.fst$subfacet <- ave.fst$comparison
    ave.fst$comparison <- NULL
    ave.fst <- ave.fst[,c("facet", "subfacet", "snp.facet", "snp.subfacet", "weighted_mean_fst")]
    x <- merge.snpR.stats(x, ave.fst, "weighted.means")
    # format averages and merge
  }
  x <- update_calced_stats(x, facets, "fst", "snp")

  return(x)
}



#'@export
#'@describeIn calc_single_stats observed heterozygosity
calc_ho <- function(x, facets = NULL){
  func <- function(gs){
    #identify heterozygote rows in genotype matrix
    genos <- colnames(gs$gs)
    hets <- which(substr(genos, 1, 1) != substr(genos, 2,2))

    # calculate ho
    ## if only one heterozygote...
    if(length(hets) == 1){
      ho <- gs$gs[,hets]/rowSums(gs$gs)
    }
    ## if no heterozygotes
    else if(length(hets) == 0){
      ho <- rep(0, nrow(gs$gs))
    }
    ## normally
    else{
      ho <- rowSums(gs$gs[,hets])/rowSums(gs$gs)
    }
  }
  
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }

  # add any missing facets
  ofacets <- facets
  facets <- check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    invisible(utils::capture.output(x <- add.facets.snpR.data(x, facets)))
  }

  out <- apply.snpR.facets(x,
                           facets = facets,
                           req = "gs",
                           fun = func,
                           case = "ps")
  colnames(out)[ncol(out)] <- "ho"
  
  x <- merge.snpR.stats(x, out)
  x <- calc_weighted_stats(x, ofacets, type = "single", "ho")
  x <- update_calced_stats(x, facets, "ho", "snp")
  return(x)
}

#'@export
#'@describeIn calc_single_stats find private alleles
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

  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }

  # add any missing facets
  facets <- check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    invisible(utils::capture.output(x <- add.facets.snpR.data(x, facets)))
  }

  out <- apply.snpR.facets(x, facets, "meta.gs", func, case = "ps.pf")
  colnames(out)[ncol(out)] <- "pa"
  x <- merge.snpR.stats(x, out)
  x <- update_calced_stats(x, facets, "pa", "snp")
  return(x)
}


#'Pairwise LD from SNP data.
#'
#'\code{calc_pairwise_ld} calculates LD between each pair of SNPs.
#'
#'Calculates pairwise linkage disequilibrium between pairs of SNPs using several
#'different methods. By default uses the Burrow's Composite Linkage
#'Disequilibrium method.
#'
#'If cld is not "only", haplotypes are estimated either via direct count after
#'removing all "0101" double heterozygote haplotypes (if use.ME is FALSE) or via
#'the Minimization-Expectation method  described in Excoffier, L., and Slatkin,
#'M. (1995). Note that while the latter method is likely more accurate, it can
#'be \emph{very} slow and often produces qualitatively equivalent results, and
#'so is not prefered during casual or preliminary analysis. Either method will
#'calculate D', r-squared, and the p-value for that r-squared.
#'
#'Since this process involves many pairwise comparisons, it can be very slow.
#'
#'In contrast, Burrow's Composite Linkage Disequilibrium (CLD) can be caluclated
#'very quickly via the \code{\link{cor}} function from base R.
#'\code{LD_full_pairwise} will perform this method alongside the other methods
#'if cld = TRUE and by itslef if cld = "only". For most analyses, this will be
#'sufficient and much faster than the other methods. This is the default
#'behavior.
#'
#'The data can be broken up categorically by either SNP and/or sample metadata,
#'as described in \code{\link{Facets_in_snpR}}.
#'
#'Heatmaps of the resulting data can be easily plotted using
#'\code{\link{plot_pairwise_LD_heatmap}}
#'
#'@param x snpRdata. Input SNP data. Note that a SNP column containing snp
#'  position in base pairs named 'position' is required.
#'@param facets character. Categorical metadata variables by which to break up
#'  analysis. See \code{\link{Facets_in_snpR}} for more details.
#'@param subfacets character, default NULL. Subsets the facet levels to run.
#'  Given as a named list: list(fam = A) will run only fam A, list(fam = c("A",
#'  "B"), chr = 1) will run only fams A and B on chromosome 1. list(fam = "A",
#'  pop = "ASP") will run samples in either fam A or pop ASP, list(fam.pop =
#'  "A.ASP") will run only samples in fam A and pop ASP.
#'@param ss numeric, default NULL. Number of snps to subsample.
#'@param par numeric or FALSE, default FALSE. If numeric, the number of cores to
#'  use for parallel processing.
#'@param sr logical, default FALSE. If TRUE, detailed progress reports for each
#'  subfacet will be reported.
#'@param CLD TRUE, FALSE, or "only", default "only". Specifies if the CLD method
#'  should be used either in addition to or instead of default methods. See
#'  details.
#'@param use.ME logical, default FALSE. Specifies if the
#'  Minimization-Expectation haplotype estimation should be used. See details.
#'@param sigma numeric, default 0.0001. If the ME method is used, specifies the
#'  minimum difference required between steps before haplotype frequencies are
#'  accepted.
#'
#'
#'@return a snpRdata object with linkage results stored in the pairwise.LD slot.
#'  Specifically, this slot will contain a list containing any LD matrices in a
#'  nested list broken down facet then by facet levels and a data.frame
#'  containing all pairwise comparisons, their metadata, and calculated
#'  statistics in long format for easy plotting.
#'
#'@references Dimitri Zaykin (2004). \emph{Genetic Epidemiology}
#'@references Excoffier, L., and Slatkin, M. (1995). \emph{Molecular Biology and
#'  Evolution}
#'@references Lewontin (1964). \emph{Genetics}
#'
#'@author William Hemstrom
#'@author Keming Su
#'@export
#'
#' @examples
#' \dontrun{
#' # not run, slow
#' ## CLD
#' x <- calc_pairwise_ld(stickSNPs, facets = "group.pop")
#' get.snpR.stats(x, "group.pop", "LD")
#'
#' ## standard haplotype frequency estimation
#' x <- calc_pairwise_ld(stickSNPs, facets = "group.pop", CLD = FALSE)
#' get.snpR.stats(x, "group.pop", "LD")
#' }
#'
#' # subset for specific subfacets (ASP and OPL, chromosome IX)
#' x <- calc_pairwise_ld(stickSNPs, facets = "group.pop",
#'                       subfacets = list(pop = c("ASP", "OPL"), 
#'                                        group = "groupIX"))
#' get.snpR.stats(x, "group.pop", "LD")
#' 
#' \dontrun{
#' ## not run, really slow
#' # ME haplotype estimation
#' x <- calc_pairwise_ld(stickSNPs, facets = "group.pop", 
#'                       CLD = FALSE, use.ME = TRUE,
#'                       subfacets = list(pop = c("ASP", "OPL"), 
#'                                        group = "groupIX"))
#' get.snpR.stats(x, "group.pop", "LD")
#' }
calc_pairwise_ld <- function(x, facets = NULL, subfacets = NULL, ss = FALSE,
                             par = FALSE, sr = FALSE, CLD = "only", use.ME = FALSE, sigma = 0.0001){
  #========================sanity checks=============
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }

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
  
  if(!"position" %in% colnames(x@snp.meta)){
    stop("A column named 'postion' containing SNP positions in bp is required in the SNP metadata.\n")
  }

  #========================sub-functions=============
  #function to do LD with SNPs

  # sub functions:

  #function to correct haplotype input matrix:
  GtoH <- function(x, n){
    m1 <- matrix(as.numeric(x), nrow(x), ncol(x))
    colnames(m1) <- colnames(x)
    m1 <- cbind(as.data.frame(t(m1)), n)
    m2 <- dplyr::group_by(m1, n)
    m2 <- dplyr::summarise_all(m2, dplyr::funs(sum))
    #m2 <- m1 %>% group_by(n) %>% summarise_all(funs(sum))
    m2 <- t(as.matrix(m2[,-1]))
    colnames(m2) <- sort(unique(n))
    return(m2)
  }

  # em haplotype estimation
  multi_haplotype_estimation <- function(x, haptable, sigma = 0.0001){

    # find the double het. Should be able to use an approach like this when this gets extended to work with everything.
    # cj values for each possible genotype:
    s1 <- substr(colnames(x), 1, 1)
    s2 <- substr(colnames(x), 2, 2)
    s3 <- substr(colnames(x), 3, 3)
    s4 <- substr(colnames(x), 4, 4)
    het.1 <- s1 != s2
    het.2 <- s3 != s4



    # First, make a guess at the starting haplotype frequencies. We'll do this by taking the unambigious haplotype frequencies,
    # then making a guess at the haplotype composition in the double heterozygote assuming that all possible haplotypes are equally likely
    doub.het <- which(het.1 + het.2 == 2) # identify double heterozygotes

    # if there are no double heterozygotes, we already no the haplotype frequencies, so we can just return those.
    if(length(doub.het) == 0){
      return(haptable/rowSums(haptable))
    }


    if(length(doub.het) != 1){
      doub.het.sum <- rowSums(x[,doub.het])
    }
    else{
      doub.het.sum <- x[,doub.het]
    }

    nhap.counts <- haptable # grab the haplotypes
    ehap.counts <- nhap.counts + .5*doub.het.sum # assuming that both haplopairs are equaly likely in the double het
    shap.freqs <- ehap.counts/rowSums(ehap.counts) # get the starting haplotype frequencies



    # now that we have our starting conditions, we will do the EM loops.
    # 1) First, we find out how many of each haplotype
    # we expect to get from our double heterozygotes given the initial haplotype frequencies we guessed above.
    # 2) Then, we use those expected frequencies to update our estimates of the haplotype frequencies.
    # 3) repeat 1 and 2 until the difference between the haplotype frequencies between loop iterations is less than sigma.


    # we'll use a while loop, which will run as long as the conditions are met. Note that this can freeze your computer if
    # the conditions are NEVER met! Try just hitting the stop sign, if that doesn't work you'll need to restart Rstudio.

    out <- matrix(NA, nrow(haptable), ncol = 4)
    completed <- logical(nrow(haptable))

    diff <- sigma + 1 # initialize the diff. Doesn't matter what it is as long as it's larger than sigma.
    check <- 1

    while(any(diff > sigma)){

      if(is.null(nrow(shap.freqs))){
        shap.freqs <- matrix(shap.freqs, 1)
        x <- matrix(x, 1)
        haptable <- matrix(haptable, 1)
      }

      # 1)
      # expectation, which is that we are drawing haplotypes (aka alleles) from a pool of options. Follows HWE, essentially,
      # but the "alleles" are actually haplotypes
      op1.e <- (2*shap.freqs[,1]*shap.freqs[,4])/
        ((2*shap.freqs[,1]*shap.freqs[,4])+(2*shap.freqs[,2]*shap.freqs[,3])) # percentage of AC/GG haplo pairs
      op2.e <- 1 - op1.e

      # maximization: given the expected haplotype frequencies, how many of each haplotype should we have? get new frequencies
      n1hap.freqs <- haptable # grab the known haplotype frequencies form the unambigious phenotypes again.
      if(nrow(x) == 1){
        n1hap.freqs[,c(1, 4)] <- n1hap.freqs[,c(1, 4)] + (rowSums(matrix(x[,doub.het], 1))*op1.e*.5) # we basically add the expected number of haplotypes for the double heterozygotes
        n1hap.freqs[,c(2, 3)] <- n1hap.freqs[,c(2, 3)] + (rowSums(matrix(x[,doub.het], 1))*op2.e*.5)
      }
      else{
        n1hap.freqs[,c(1, 4)] <- n1hap.freqs[,c(1, 4)] + (doub.het.sum*op1.e*.5) # we basically add the expected number of haplotypes for the double heterozygotes
        n1hap.freqs[,c(2, 3)] <- n1hap.freqs[,c(2, 3)] + (doub.het.sum*op2.e*.5)
      }
      n1hap.freqs <- n1hap.freqs/rowSums(n1hap.freqs)

      # for rows where we end up with NAs, just use the starting estimated haplotypes instead (we're done)
      na.rows <- which(is.na(op1.e))
      if(length(na.rows) > 0){
        n1hap.freqs[na.rows,] <- shap.freqs[na.rows,]
      }

      # calculate the diff and update
      diff <- rowSums(abs(n1hap.freqs - shap.freqs))
      check <- check + 1

      # save any completed results
      done <- which(diff <= sigma)
      if(length(done) > 0){
        if(sum(completed) > 0){
          if(nrow(n1hap.freqs) > 1){
            out[-which(completed),][done,] <- n1hap.freqs[done,]
          }
          else{
            out[-which(completed),] <- n1hap.freqs[done,]
          }
          completed[-which(completed)][done] <- T
        }
        else{
          out[done,] <- n1hap.freqs[done,]
          completed[done] <- T
        }
        n1hap.freqs <- n1hap.freqs[-done,]
        shap.freqs <- n1hap.freqs
        haptable <- haptable[-done,]
        x <- x[-done,]
        if(is.null(nrow(x))){
          doub.het.sum <- sum(x[doub.het])
        }
        else if(length(doub.het) != 1){
          doub.het.sum <- rowSums(x[,doub.het])
        }
        else{
          doub.het.sum <- x[,doub.het]
        }
      }
      else{
        shap.freqs <- n1hap.freqs
      }
    }



    # return the output
    return(out)
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
    ## Keming: this line can be edited--remove the last which statement, save as a new variable which identifes
    ##         columns containing missing data.
    rgcs <- c(grep(paste0("^", dmDat), gl), #missing first locus
              grep(paste0(dmDat, "$"), gl), #missing second locus
              which(substr(gl, 1, sform) != substr(gl, (sform + 1), (sform *2)) &
                      substr(gl, (sform*2) + 1, sform*3) != substr(gl, (sform*3+1), sform*4))) #double het


    # a version for the em method
    if(use.ME){
      new.var <- grep(dmDat, gl)
      if(is.null(nrow(ghapmat))){
        ghapmat2 <- as.matrix(ghapmat, nrow = 1)[,-new.var]
      }
      ghapmat2 <- ghapmat[,-new.var]
    }

    ##remove any double heterozygotes
    ghapmat <- ghapmat[,-rgcs]


    ## Keming: make a ghapmat with the new variable instead of rgcs. This will be x in our multihaplotype estimation function

    #add a filler row for the last pairwise comparison to make life easier.
    if(length(y) == length(x)){
      if(length(ghapmat) > 1){ #stop it from doing this if there is data for only one haplotype.
        ghapmat <- rbind(ghapmat, rep(c(10,0), 100)[1:length(ghapmat)])
        if(use.ME){
          ghapmat2 <- rbind(ghapmat2, rep(c(10,0), 100)[1:length(ghapmat2)])
        }
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
      m2matv <- stats::na.omit(m2matv)
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
    hapmat <- stats::na.omit(hapmat)
    #if(length(hapmat) %% nrow(y) != 0){
    #  browser()
    #}
    hapmat <- matrix(hapmat, nrow(ghapmat), 4, byrow = T)

    ## Keming: the hapmat variable is the second argument, haptable
    ##         now, what want to do is take an argument (let's name it use.ME). If use.ME == TURE,
    ##         we want to call our ME function, get a hapmat out from it, and return this instead of the
    ##         hapmat that we calculated above.

    if(use.ME){
      hapmat <- multi_haplotype_estimation(ghapmat2, haptable = hapmat, sigma = sigma)
      hapmat <- hapmat*rowSums(ghapmat2) # convert back to numbers rather than frequencies.
      #call our new haplotype function, x is ghapmat2, haptable is hapmat, sigma is sigma
      # overwrite the hapmat object.
    }

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
      pval <- 1 - stats::pchisq(chisqu, 1)

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
  # this function also does the composite LD calculations, since that's most efficiently done here for each subfacet level.
  determine.comparison.snps <- function(x, facets, facet.types){


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
        out[[".base"]][[".base"]] <- LD_matrix[[which(names(LD_matrix) == facets[i])]]
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
        cl <- snow::makeSOCKcluster(par)
        doSNOW::registerDoSNOW(cl)

        #prepare reporting function
        ntasks <- par
        progress <- function(n) cat(sprintf("Part %d out of",n), ntasks, "is complete.\n")
        opts <- list(progress=progress)

        # initialize and store things
        x_storage <- as.matrix(as.data.frame(x))
        na.test <- suppressWarnings(as.numeric(x_storage[1]))
        #save the info as a bigmatrix if it can be safely converted to numeric. Usually this is true for ms but not necessarily other data types.
        if(!is.na(na.test)){
          if(as.numeric(x_storage[1]) == x_storage[1]){
            cat("Saving matrix as big.matrix object for quicker sharing.\n")
            xb <- bigmemory::as.big.matrix(x, type = "char")
            xbd <- bigmemory::describe(xb)
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
        prox$sample.facet <- ".base"
        prox$sample.subfacet <- ".base"

        out <- list(prox = prox, LD_matrices = LD_mats)

        return(out)
      }

      #otherwise run normally
      else{
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
        prox$sample.facet <- ".base"
        prox$sample.subfacet <- ".base"
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
          w_list$prox <- rbind(w_list$prox, cbind(out$prox, sample.facet = names(comps)[i], sample.subfacet = names(comps[[i]])[j]))
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
        w_list$prox <- rbind(w_list$prox, cbind(output[[i]]$prox, sample.facet = names(comps)[t.facet], sample.subfacet = names(comps[[t.facet]])[t.subfacet]))
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
  # add any missing facets
  x <- add.facets.snpR.data(x, facets)

  #subset data if requested:
  if(!(is.null(subfacets[1]))){
    old.x <- x
    ssfacets <- names(subfacets)

    # check complex facets
    complex.sfacets <- check.snpR.facet.request(x, ssfacets, remove.type = "simple", fill_with_base = FALSE)
    if(length(complex.sfacets) > 0){
      stop(paste0("Complex (snp and sample) level SUBFACETS not accepted. Providing these as seperate snp and sample subfacets will run only snps/samples contained in both levels. Bad facets: ",
                  paste0(complex.sfacets, collapse = ", "), "\n"))
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
    filter_snp_facets <- F
    filter_samp_facets <- F
    
    sample.facets <- names(subfacets)[which(ssfacet.types == "sample")]
    if(length(sample.facets) > 0){
      sample.subfacets <- subfacets[[which(ssfacet.types == "sample")]]
      filter_samp_facets <- T
    }
    snp.facets <- names(subfacets)[which(ssfacet.types == "snp")]
    if(length(snp.facets) > 0){
      snp.subfacets <- subfacets[[which(ssfacet.types == "snp")]] 
      filter_snp_facets <- T
    }
    
    if(filter_samp_facets){
      if(filter_snp_facets){
        invisible(utils::capture.output(x <- subset_snpR_data(x,
                                                       facets = sample.facets,
                                                       subfacets = sample.subfacets,
                                                       snp.facets = snp.facets,
                                                       snp.subfacets = snp.subfacets)))
      }
      else{
        invisible(utils::capture.output(x <- subset_snpR_data(x,
                                                       facets = sample.facets,
                                                       subfacets = sample.subfacets)))
      }
    }
    else if(filter_snp_facets){
      invisible(utils::capture.output(x <- subset_snpR_data(x,
                                                     snp.facets = snp.facets,
                                                     snp.subfacets = snp.subfacets)))
    }
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
    x <- subset_snpR_data(x, ss)
  }

  # run non-CLD LD components:
  if(CLD != "only"){
    # typical facet check, keeping all facet types but removing duplicates. Also returns the facet type for later use.
    facets_trad <- check.snpR.facet.request(x, facets, remove.type = "none", return.type = T)

    # run the function
    out <- func(x, facets = facets_trad, snp.facets = snp.facets, par = par, sr = sr)

    # add to snpRdata object and return
    if(exists("old.x")){
      out <- merge.snpR.stats(old.x, out, "LD")

    }
    else{
      out <- merge.snpR.stats(x, out, "LD")
    }
  }
  # run CLD components
  if(CLD != F){
    # run the function
    out <- calc_CLD(x, facets, par)

    # add to snpRdata object and return
    if(CLD != "only"){
      out <- merge.snpR.stats(x, out, "LD")
    }
    else{
      if(exists("old.x")){
        out <- merge.snpR.stats(old.x, out, "LD")

      }
      else{
        out <- merge.snpR.stats(x, out, "LD")
      }
    }
  }
  out <- update_calced_stats(out, facets, "LD")

  return(out)
}

#' Calculate Burrow's composite LD. Internal.
#'
#' Called in \code{\link{calc_pairwise_ld}}. Complex function purely because the
#' output needs to be in the same format as that function's and to support
#' faceting without excessive recalculation.
#'
#' @param x snpRdata object
#' @param facets facets to run
#' @param par number of parallel cores
#'
#' @author William Hemstrom
calc_CLD <- function(x, facets = NULL, par = FALSE){
  proximity <- s1_position <- s2_position <- NULL
  
  #============subfunctions==============
  # calculate composite LD for a single facet of data.
  do_CLD <- function(genos, snp.meta, sample.facet, sample.subfacet){
    melt_cld <- function(CLD, snp.meta, sample.facet, sample.subfacet){
      prox <- cbind(as.data.table(snp.meta), as.data.table(CLD))
      prox <- reshape2::melt(prox, id.vars = colnames(snp.meta))
      prox <- cbind(prox, as.data.table(snp.meta[rep(1:nrow(snp.meta), each = nrow(snp.meta)),]))
      bad.col <- which(colnames(prox) == "variable")
      prox <- prox[,-bad.col, with = FALSE]
      colnames(prox) <- c(paste0("s1_", colnames(snp.meta)), "CLD", paste0("s2_", colnames(snp.meta)))
      prox <- stats::na.omit(prox)
      prox <- as.data.table(prox)
      prox[,proximity := abs(s1_position - s2_position)]
      prox[,sample.facet := sample.facet]
      prox[,sample.subfacet := sample.subfacet]
      setcolorder(prox, c(1:ncol(snp.meta),
                          (ncol(snp.meta) + 2):(ncol(prox) - 3),
                          ncol(prox) - 2,
                          ncol(snp.meta) + 1,
                          (ncol(prox) - 1):ncol(prox)))
      return(prox)
    }

    # do the CLD calculation, add column/row names, NA out the lower triangle and diag
    ## CLD
    suppressWarnings(CLD <- stats::cor(t(genos), use = "pairwise.complete.obs")^2)
    ## matrix of complete case sample sizes
    complete.cases.matrix <- !is.na(t(genos))
    complete.cases.matrix <- crossprod(complete.cases.matrix)

    # fill in NAs
    CLD[which(lower.tri(CLD))] <- NA
    diag(CLD) <- NA
    complete.cases.matrix[is.na(CLD)] <- NA

    # add metadata and melt
    prox <- melt_cld(CLD, snp.meta, sample.facet, sample.subfacet)
    prox_S <- melt_cld(complete.cases.matrix, snp.meta, sample.facet, sample.subfacet)
    colnames(prox_S)[which(colnames(prox_S) == "CLD")] <- "S"

    # merge
    prox <- merge.data.table(prox, prox_S)

    # add column/row names
    colnames(CLD) <- snp.meta$position
    rownames(CLD) <- snp.meta$position
    return(list(prox = prox, LD_matrix = CLD, S = complete.cases.matrix))
  }

  # take an output lists of matrices and prox tables and sort and name them for merging.
  decompose_outputs <- function(matrix_storage, prox_storage, tasks){
    # figure out the facet names
    facet.names <- paste(tasks[,1], tasks[,3], sep = ".")
    facet.names <- gsub("\\.base", "", facet.names)
    facet.names <- gsub("\\.\\.", "\\.", facet.names)
    facet.names <- gsub("^\\.", "", facet.names)
    for(i in 1:length(facet.names)){
      if(facet.names[1] != ""){
        facet.names[i] <- check.snpR.facet.request(x, facet.names[i], remove.type = "none")
      }
      else{
        facet.names[i] <- ".base"
      }
    }

    # initialize
    matrix_out <- vector("list", length(unique(facet.names)))
    names(matrix_out) <- unique(facet.names)

    # decompose each matrix
    ## first, initialize storage with all of the correctly names slots
    for(i in 1:length(unique(facet.names))){
      # grab the tasks with this facet
      these.tasks <- which(facet.names == unique(facet.names)[i])

      # pop options in this facet
      pops <- unique(tasks[these.tasks,2])
      matrix_out[[i]] <- vector("list", length(pops))
      names(matrix_out[[i]]) <- pops

      # snp facet options in this facet
      snp.levs <- unique(tasks[these.tasks,4])
      for(j in 1:length(matrix_out[[i]])){
        matrix_out[[i]][[j]] <- vector("list", length(snp.levs))
        names(matrix_out[[i]][[j]]) <- snp.levs
        for(k in 1:length(snp.levs)){
          matrix_out[[i]][[j]][[k]] <- vector("list", 2)
          names(matrix_out[[i]][[j]][[k]]) <- c("CLD", "S")
        }
      }
    }
    ## then fill
    for(i in 1:nrow(tasks)){
      matrix_out[[facet.names[i]]][[tasks[i,2]]][[tasks[i,4]]][["CLD"]] <- matrix_storage[[i]][["CLD"]]
      matrix_out[[facet.names[i]]][[tasks[i,2]]][[tasks[i,4]]][["S"]] <- matrix_storage[[i]][["S"]]
    }

    # rbind the prox together
    prox <- data.table::rbindlist(prox_storage)

    return(list(prox = prox, LD_matrices = matrix_out))
  }

  #============run=======================
  # get task list and do a conversion to sn
  cat("Preparing data...\n")

  suppressMessages(x <- add.facets.snpR.data(x, facets))
  tasks <- get.task.list(x, facets)
  x@sn$sn <- format_snps(x, "sn", interpolate = F)
  x@sn$type <- "FALSE"

  cat("Beginning LD calculation...\n")
  # run the loop
  if(par == F){

    # initialize storage
    matrix_storage <- vector("list", nrow(tasks))
    prox_storage <- vector("list", nrow(tasks))


    #loop through each set of facets
    for(i in 1:nrow(tasks)){
      # run
      cat("Task #:", i, "of", nrow(tasks),
          " Sample Facet:", paste0(tasks[i,1:2], collapse = "\t"),
          " SNP Facet:", paste0(tasks[i,3:4], collapse = "\t"), "\n")

      # can't integrate this part into do_CLD without screwing up the parallel due to s4 issues.
      suppressWarnings(y <- subset_snpR_data(x, facets = tasks[i,1],
                                             subfacets = tasks[i,2],
                                             snp.facets = tasks[i,3],
                                             snp.subfacets = tasks[i,4]))
      
      out <- do_CLD(y@sn$sn[,-c(1:(ncol(y@snp.meta) - 1))], y@snp.meta, tasks[i, 1], tasks[i, 2])
      
      # extract
      matrix_storage[[i]] <- list(CLD = out$LD_matrix, S = out$S)
      prox_storage[[i]] <- out$prox
    }

    # decompose
    return(decompose_outputs(matrix_storage, prox_storage, tasks))

  }
  else if(is.numeric(par)){
    cat("Running in parallel.\nSplitting up data...\n")

    geno.storage <- vector("list", nrow(tasks))


    # can't do this part in parallel due to s4 issues.
    for(i in 1:nrow(tasks)){
      cat("Task #:", i, "of", nrow(tasks),
          " Sample Facet:", paste0(tasks[i,1:2], collapse = "\t"),
          " SNP Facet:", paste0(tasks[i,3:4], collapse = "\t"), "\n")
      utils::capture.output(invisible(suppressWarnings(y <- subset_snpR_data(x, facets = tasks[i,1],
                                                                             subfacets = tasks[i,2],
                                                                             snp.facets = tasks[i,3],
                                                                             snp.subfacets = tasks[i,4]))))

      geno.storage[[i]] <- list(geno = y@sn$sn[,-c(1:(ncol(y@snp.meta) - 1))], snp.meta = y@snp.meta)
    }

    # split up
    tasks <- as.data.frame(tasks, stringsAsFactors = F)
    tasks$ord <- 1:nrow(tasks)
    if(par < nrow(tasks)){
      ptasks <- split(tasks, rep(1:par, length.out = nrow(tasks), each = ceiling(nrow(tasks)/par)))
    }
    else{
      par <- nrow(tasks)
      ptasks <- split(tasks, 1:par, drop = F)
    }

    cl <- snow::makeSOCKcluster(par)
    doSNOW::registerDoSNOW(cl)



    #prepare reporting function
    ntasks <- length(ptasks)
    progress <- function(n) cat(sprintf("Job %d out of", n), ntasks, "is complete.\n")
    opts <- list(progress=progress)

    #loop through each set of facets
    output <- foreach::foreach(q = 1:ntasks,
                               .packages = c("dplyr", "reshape2", "matrixStats", "bigtabulate", "snpR", "data.table"),
                               .inorder = TRUE,
                               .options.snow = opts) %dopar% {

                                 tasks <- ptasks[[q]]

                                 matrix_storage <- vector("list", nrow(tasks))
                                 prox_storage <- vector("list", nrow(tasks))

                                 for(i in 1:nrow(tasks)){

                                   # run
                                   out <- do_CLD(genos = geno.storage[[tasks[i,"ord"]]]$geno,
                                                 snp.meta = geno.storage[[tasks[i,"ord"]]]$snp.meta,
                                                 sample.facet = tasks[i, 1], sample.subfacet = tasks[i, 2])

                                   # extract
                                   matrix_storage[[i]] <- list(CLD = out$LD_matrix, S = out$S)
                                   prox_storage[[i]] <- out$prox
                                 }

                                 list(prox = prox_storage, matrix = matrix_storage)
                               }

    #release cores and clean up
    parallel::stopCluster(cl)
    doSNOW::registerDoSNOW()
    gc();gc()

    # split apart matrices and decompose
    matrix_out <- vector("list", length = length(output))
    prox <- vector("list", length = length(output))
    for(i in 1:length(output)){
      matrix_out[[i]] <- output[[i]][[seq(2, length(output[[i]]), by = 2)]]
      prox[[i]] <- output[[i]][[seq(1, length(output[[i]]), by = 2)]]
    }
    prox <- unlist(prox, recursive = F)
    matrix_out <- unlist(matrix_out, recursive = F)
    return(decompose_outputs(matrix_out, prox, dplyr::bind_rows(ptasks)[,-5]))
  }
}

#'@export
#'@describeIn calc_single_stats p-values for Hardy-Wienberg Equilibrium divergence
calc_hwe <- function(x, facets = NULL, method = "exact", 
                     fwe_method = "BY", 
                     fwe_case = c("by_facet", "overall")){
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
    o2pq <- rowSums(gs[,het.col, drop = F])
    opp <- matrixStats::rowMaxs(gs[,-het.col, drop = F])
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
      out <- stats::pchisq(chi, 2, lower.tail = FALSE)
    }

    # otherwise we have to use the looped version:
    else if(method == "exact"){
      cat("Using exact test from Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005).\n")
      out <- numeric(nrow(gs))
      for(i in 1:nrow(gs)){
        out[i] <- exact.hwe(oqq[i], opp[i], o2pq[i])
      }

      nas <- which(out == -1)
      if(length(nas) > 1){
        out[nas] <- NA
      }
    }
    return(out)
  }

  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  if(!(method %in% c("exact", "chisq"))){stop("Unrecognized HWE method, please use chisq or exact.\n")}

  # add any missing facets
  facets <- check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    invisible(utils::capture.output(x <- add.facets.snpR.data(x, facets)))
  }

  out <- apply.snpR.facets(x, facets, "gs", func, case = "ps", method = method)
  colnames(out)[ncol(out)] <- "pHWE"
  
  # apply corrections
  if(length(facets) == 1 & facets[1] == ".base"){
    fwe_case <- "overall"
  }
  out <- fwe_correction(out, levs = c("facet", "subfacet"), pcol = "pHWE", methods = fwe_method, case = fwe_case)
  
  
  
  x <- merge.snpR.stats(x, out)
  x <- update_calced_stats(x, facets, "hwe", "snp")

  return(x)
}

#'Caluclate basic SNP statistics
#'
#'Automatically calculate most basic statistics from snpRdata. Calculates maf,
#'pi, ho, pairwise Fst, HWE divergence, finds private alleles, and uses Gaussian
#'smoothing to produce per-window averages of all of these.
#'
#'The data can be broken up categorically by sample or SNP metadata, as
#'described in \code{\link{Facets_in_snpR}}. Note that Fst and private allele
#'calculations require a sample specific contrast (the pairwise part of pairwise
#'Fst), and thus will not be calculated unless a facet with a sample meta data
#'variable included is specified. The other stats, in contrast, are snp-specific
#'and thus ignore any snp meta data variables included in facets. Providing NULL
#'or "all" to the facets argument works as described in
#'\code{\link{Facets_in_snpR}}.
#'
#'@param x snpRdata object.
#'@param facets character. Categorical metadata variables by which to break up
#'  analysis. See \code{\link{Facets_in_snpR}} for more details.
#'@param fst.method character, default "WC". Defines the FST estimator to use.
#'  Options: \itemize{ \item{WC: } Wier and Cockerham (1984). \item{Wier: } Wier
#'  (1990) \item{Hohenlohe: } Hohenlohe et al (2010), identical to the STACKS
#'  package. \item{Genepop: } Rousset (2008), uses the genepop package. }
#'@param sigma numeric. Designates the width of windows in kilobases. Full
#'  window size is 6*sigma.
#'@param step numeric or NULL, default NULL. Designates the number of kilobases
#'  between each window centroid. If NULL, windows are centered on each SNP.
#'@param par numeric or FALSE, default FALSE. If numeric, the number of cores to
#'  use for parallel processing.
#'@param nk logical, default TRUE. If TRUE, weights SNP contribution to window
#'  averages by the number of observations at those SNPs.
#'
#'@examples
#'x <- calc_basic_snp_stats(stickSNPs, c("pop"))
#'get.snpR.stats(x, "pop") # view basic stats
#'get.snpR.stats(x, "pop", type = "pairwise") # view fst
#'
#'
#'@author William Hemstrom
#'@export
#'
#'@seealso calc_single_stats calc_pairwise_fst calc_smoothed_averages
#'@return A snpRdata object with all of the described statistics merged into the
#'  appropriate sockets.
#'  
calc_basic_snp_stats <- function(x, facets = NULL, fst.method = "WC", sigma = NULL, step = NULL, par = FALSE, nk = TRUE){
  #=========sanity checks=======
  if(!is.snpRdata(x)){
    stop("x must be a snpRdata object.\n")
  }
  
  if(!is.null(facets[1])){
    facet.types <- check.snpR.facet.request(x, facets, "none", T)
    snp.facets <- which(facet.types[[2]] == "snp")
  }

  if(!is.null(sigma)){
    if(is.null(facets[1])){
      sanity_check_window(x, sigma, step, nk = nk, stats.type = "single", facets = facets)
    }
    else if(any(facet.types[[2]] == "snp")){
      sanity_check_window(x, sigma, step, nk = nk, stats.type = "single", facets = facets)
    }
    else{
      sanity_check_window(x, sigma, step, nk = nk, stats.type =  c("pairwise", "single"), facets = facets)
    }
  }

  #=========stats===============
  # basic stats
  x <- calc_maf(x, facets)
  x <- calc_pi(x, facets)
  x <- calc_hwe(x, facets)
  x <- calc_ho(x, facets)
  if(!is.null(facets[1]) & !any(check.snpR.facet.request(x, facets, return.type = TRUE)[[2]] == ".base")){
    x <- calc_pairwise_fst(x, facets, method = fst.method)
    x <- calc_private(x, facets)
  }

  # if a snp facet level is requested, run everything at the base level too.
  if(!is.null(facets[1])){
    if(length(snp.facets) > 0){
      x <- calc_maf(x)
      x <- calc_pi(x)
      x <- calc_hwe(x)
      x <- calc_ho(x)
    }
  }
  
  #=========smoothing===========
  if(!is.null(sigma)){

    # if no facets, run only non-pairwise
    if(is.null(facets[1])){
      x <- calc_smoothed_averages(x, facets, sigma, step, nk, par = par, stats.type = "single")
    }

    # otherwise, need to run any snp facets with only single, everything else with pairwise.
    else{
      if(length(snp.facets) > 0){
        x <- calc_smoothed_averages(x, facets[snp.facets], sigma, step, nk, par = par, stats.type = "single")
        if(length(snp.facets) != length(facets)){
          x <- calc_smoothed_averages(x, facets[-snp.facets], sigma, step, nk, par = par)
        }
      }
      else{
        x <- calc_smoothed_averages(x, facets, sigma, step, nk, par = par)
      }
    }
  }

  #=========return==============
  return(x)
}

#' Calculate the ratio of heterozygous/homozygous sites per individual.
#'
#' Calculates the ratio of heterozygotes to homozygous sites across all SNPs
#' within each individual.
#'
#' @param x snpRdata object
#' @param facets facets over which to split snps within samples. Takes only SNP
#'   level facets. See \code{\link{Facets_in_snpR}} for details.
#' @param complex_averages logical, default FALSE. If TRUE, will compute weighted averages for
#'   complex (snp + sample metadata) facets. This can be quite time consuming, and so is generally
#'   not recommended.
#'
#' @return A snpRdata object with heterozygote/homozygote ratios merged into the
#'   sample.stats slot.
#'
#' @author William Hemstrom
#' @export
#'
#' @examples
#' # base facet
#' x <- calc_het_hom_ratio(stickSNPs)
#' get.snpR.stats(x, type = "sample")
#'
#' # facet by chromosome
#' x <- calc_het_hom_ratio(stickSNPs, "group")
#' get.snpR.stats(x, "group", type = "sample")
#'
#'
#' 
calc_het_hom_ratio <- function(x, facets = NULL, complex_averages = FALSE){
  func <- function(x, mDat){
    # make x into a logical for heterozygous
    xv <- as.matrix(x)
    logix <- substr(xv, 1, 1) != substr(xv, 2, 2) # true if het
    logix[xv == mDat] <- NA # NA when missing data!

    # get counts of hets and homs, then ratio
    hets <- matrixStats::colSums2(logix, na.rm = T) # number heterozygous sites
    homs <- matrixStats::colSums2(!logix, na.rm = T) # number homozygous sites
    ratio <- hets/homs
    return(ratio)
  }

  #============run for each facet================
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  # add any missing facets
  ofacets <- facets
  if(!complex_averages){
    ofacets <- check.snpR.facet.request(x, facets, "complex")
  }
  facets <- check.snpR.facet.request(x, facets, remove.type = "sample")

  out <- apply.snpR.facets(x, facets, "genotypes", func, case = "per_sample", mDat = x@mDat)
  colnames(out)[which(colnames(out) == "stat")] <- "Het/Hom"
  x <- merge.snpR.stats(x, out, "sample.stats")
  x <- calc_weighted_stats(x, ofacets, "sample", "Het/Hom")
  
  x <- update_calced_stats(x, facets, "ho_he_ratio")

  return(x)
}




#' Calculate effective population size.
#'
#' Calculates effective population size for any given sample-level facets via
#' interface with the NeEstimator v2 program by Do et al. (2013).
#'
#' Since physical linkage can cause mis-estimation of Ne, an optional snp-level
#' facet can be provided which designates chromosomes or linkage groups. Only
#' pairwise LD values between SNPs on different facet levels will be used.
#'
#' Ne can be calculated via three different methods: \itemize{ \item{"LD"}
#' Linkage Disequilibrium based estimation. \item{"Ht"} Heterozygote excess.
#' \item{"Coan"} Coancestry based.} For details, please see the documentation
#' for NeEstimator v2.
#'
#' @param x snpRdata object. The data for which Ne will be calculated.
#' @param facets character, default NULL. Categorical metadata variables by
#'   which to break up analysis. See \code{\link{Facets_in_snpR}} for more
#'   details. Only sample specific categories are allowed, all others will be
#'   removed. If NULL, Ne will be calculated for all samples.
#' @param chr character, default NULL. An optional but recommended SNP specific
#'   categorical metadata variable which designates chromosomes/linkage
#'   groups/etc. Pairwise LD scores for SNPs with the same level of this
#'   variable will be not be used to calculate Ne. Since physical linkage can
#'   bias Ne estimtes, providing a value here is recommended.
#' @param NeEstimator_path character, default "/usr/bin/Ne2-1.exe". Path to the
#'   NeEstimator executable.
#' @param mating character, default "random". The mating system to use. Options:
#'   \itemize{ \item{"random"} Random mating. \item{"monogamy"} Monogamous
#'   mating. }
#' @param pcrit numeric, default c(.05, .02, .01). Minimum minor allele
#'   frequencies for which to calculate Ne. Rare alleles can bias estimates, so
#'   a range of values should be checked.
#' @param methods character, default "LD". LD estimation methods to use.
#'   Options: \itemize{ \item{"LD"} Linkage Disequilibrium based estimation.
#'   \item{"Ht"} Heterozygote excess. \item{"Coan"} Coancestry based.}
#' @param temporal_methods character, default c("Pollak", "Nei", "Jorde").
#'   Methods to use for the temporal methods. See defaults for options. Not
#'   currently supported.
#' @param temporal_gens data.frame, default NULL. Work in progress, not
#'   currently supported.
#' @param max_ind_per_pop numeric, default NULL. Maximum number of individuals
#'   to consider per population.
#' @param outfile character, default "ne_out". Prefix for output files. Note
#'   that this function will return outputs, so there isn't a strong reason to
#'   check this.
#'
#' @return A named list containing estimated Ne values (named "ne") and the
#'   original provided data, possibly with additional LD values (named "x").
#'
#' @author William Hemstrom
#' @export
#'
#' @references Do C, Waples RS, Peel D, Macbeth GM, Tillett BJ, Ovenden JR. 2014
#'   NeEstimator v2: re-implementation of software for the estimation of
#'   contemporary effective population size (Ne) from genetic data. Mol. Ecol.
#'   Resour. 14, 209–214. (doi:10.1111/1755-0998.12157)
#'
#' @examples
#' \dontrun{
#' # not run, since the path to NeEstimator may vary
#' # calculate Ne, noting not to use LD between SNPs on the same chromosome equivalent ("group") for every population.
#' ne <- calc_ne(stickSNPs, facets = "pop", chr = "group") # need to set the path argument to the local NeEstimator installation.
#' get.snpR.stats(ne, "pop", type = "pop")}
#' 
#' @export
calc_ne <- function(x, facets = NULL, chr = NULL,
                    NeEstimator_path = "/usr/bin/Ne2-1.exe",
                    mating = "random",
                    pcrit = c(0.05, 0.02, 0.01),
                    methods = "LD",
                    temporal_methods = c("Pollak", "Nei", "Jorde"),
                    temporal_gens = NULL, max_ind_per_pop = NULL,
                    outfile = "ne_out"){
  #==============sanity checks and prep=================
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  msg <- character()
  if(!file.exists(NeEstimator_path)){
    msg <- c(msg, "NeEstimator executable not found at provided path.\n")
  }
  
  facets <- check.snpR.facet.request(x, facets, remove.type = "snp")
  if(length(msg) > 0){
    stop(msg)
  }
  
  
  #=============run====================================
  out <- vector("list", length(facets))
  for(i in 1:length(facets)){
    # write inputs
    write_neestimator_inputs(x = x,
                             facets = facets[i],
                             chr = chr,
                             methods = methods,
                             temporal_methods = temporal_methods,
                             temporal_gens = temporal_gens,
                             pcrit = pcrit,
                             mating = mating,
                             outfile = paste0(outfile, ".txt"),
                             max_ind_per_pop = max_ind_per_pop)
    
    # run
    run_neestimator(NeEstimator_path = NeEstimator_path,
                    data_path = "NeEstimator/")
    
    # parse
    out[[i]] <- parse_neestimator(path = "NeEstimator/",
                                  pattern = outfile,
                                  facets = facets[i],
                                  snpRdat = x)
  }
  
  names(out) <- facets
  out <- dplyr::bind_rows(out, .id = "facet")
  x <- merge.snpR.stats(x, out, "pop")
  x <- update_calced_stats(x, facets, "ne")

  return(x)
}

#' Calculate genetic distances between individuals or groups of samples.
#'
#' Calculates the genetic distances between either individuals or the levels of
#' any sample specific facets, broken apart by any requested snp level facets.
#' For details on methods, see details.
#'
#' Available methods: \itemize{\item{Edwards} Angular distance as described in
#' Edwards 1971.}
#'
#' All methods first calculate allele frequency matrices using
#' \code{\link{tabulate_allele_frequency_matrix}}, then use these matrices
#' (hereafter x) to calculate genetic distance.
#'
#' Method details: \itemize{\item{Edwards 1971:} x <- sqrt(x); x <- x%*%t(x); x
#' <- 1/(x/number.loci); diag(x) <- 0; x <- sqrt(x)}
#'
#' @param x snpRdata object.
#' @param facets character or NULL, default NULL. Facets for which to calculate
#'   genetic distances, as described in \code{\link{Facets_in_snpR}}. If snp or
#'   base level facets are requested, distances will be between individuals.
#'   Otherwise, distances will be between the levels of the sample facets.
#' @param method character, default "Edwards". Name of the method to use.
#'   Options: \itemize{\item{Edwards} Angular distance as described in Edwards
#'   1971.} See details.
#' @param interpolate character, default "bernoulli". Missing data interpolation
#'   method, solely for individual/individual distances. Options detailed in
#'   documentation for \code{\link{format_snps}}.
#'
#' @return An overwrite-safe snpRdata object with genetic distance information
#'   (a named, nested list containing distance measures named according to
#'   facets and facet levels) added.
#' @references Edwards, A. W. F. (1971). Distances between populations on the
#'   basis of gene frequencies. Biometrics, 873-881.
#'
#' @author William Hemstrom
#' @export
#'
#' @examples
#' # by pop:
#' y <- calc_genetic_distances(stickSNPs, facets = "pop", method = "Edwards")
#' get.snpR.stats(y, "pop", "genetic_distance")
#'
#' # by group and pop jointly
#' y <- calc_genetic_distances(stickSNPs, facets = "pop.group", method = "Edwards")
#' get.snpR.stats(y, "pop.group", "genetic_distance")
#'
#' # by pop and fam seperately
#' y <- calc_genetic_distances(stickSNPs, facets = c("pop", "fam"), method = "Edwards")
#' get.snpR.stats(y, c("pop", "group"), "genetic_distance")
#'
#' # individuals across all snps + plot
#' y <- calc_genetic_distances(stickSNPs)
#' heatmap(as.matrix(get.snpR.stats(y, type = "genetic_distance")$.base$.base$Edwards))
#' 
calc_genetic_distances <- function(x, facets = NULL, method = "Edwards", interpolate = "bernoulli"){
  #============sanity checks=========
  msg <- character()
  
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }

  good.methods <- c("Edwards")
  if(!method %in% good.methods){
    msg <- c(msg, paste0("Provided method not supported. Supported methods: ", paste0(good.methods, collapse = " ")))
  }
  
  # get the allele frequencies if not already provided
  facets <- check.snpR.facet.request(x, facets, "none", T)
  snp_lev <- which(facets[[2]] == "snp" | facets[[2]] == ".base")
  snp_facets <- facets[[1]][snp_lev]
  all_facets <- facets[[1]]
  
  # pull out snp level facets, since these'll be done later.
  if(length(snp_lev) != length(facets[[1]])){
    if(length(snp_lev) > 0){
      facets <- facets[[1]][-snp_lev]
    }
    else{
      facets <- facets[[1]]
    }
    needed.facets <- check_calced_stats(x, facets, "allele_frequency_matrix")
    needed.facets <- facets[which(!unlist(needed.facets))]
    if(length(needed.facets) > 0){
      x <- tabulate_allele_frequency_matrix(x, needed.facets)
    }
    sample_facets_detected <- T
  }
  else{
    sample_facets_detected <- F
  }

  if(length(msg) > 0){
    stop(paste0(msg, collapse = "\n"))
  }
  
  # grab the data we are working with for sample specific facets
  y <- x
  if(sample_facets_detected){
    x <- .get.snpR.stats(y, facets, "allele_frequency_matrix")
  }
  #=============subfunctions=========
  # dist subfunction
  get_dist <- function(x, method){
    if(method == "Edwards"){
      x <- x[,which(colSums(is.na(x)) == 0)] # remove anywhere where there is missing data!
      nloc <- ncol(x)
      x <- sqrt(as.matrix(x))
      am <- x%*%t(x)
      am <- 1 - (am / (nloc/2))
      diag(am) <- 0
      am <- sqrt(am)
      am <- stats::as.dist(am)
    }
    am <- list(am)
    names(am) <- method
    return(am)
  }
  
  #=============run for facets with sample aggregation===============
  if(sample_facets_detected){
    out <- vector("list", length(x))
    names(out) <- names(x)
    # enter lapply hell--first level unlists once, second level runs the function if the values are matrices, unlists if not, third level can always run the function
    # the nice thing with this is that it should keep all of the names from x natively!
    out <- lapply(x, function(y){
      lapply(y, function(z) {
        if("matrix" %in% class(z)){
          get_dist(z, method = method)
        }
        else{
          lapply(z, get_dist, method = method)
        }
      })
    })
  }
  
  #============run for facets without sample aggregation (snp or .base level)
  if(length(snp_facets) > 0){
    # intialize output and get sn
    out_snp <- vector("list", length(snp_facets))
    names(out_snp) <- snp_facets
    sn <- format_snps(y, "sn", interpolate = interpolate)
    sn <- sn[,-c(1:(ncol(y@snp.meta) - 1))]
    
    # for each snp facet...
    for(i in 1:length(snp_facets)){
      if(snp_facets[i] == ".base"){
        out_snp[[i]] <- vector("list", 1)
        names(out_snp[[i]]) <- ".base"
        out_snp[[i]][[1]] <- list(stats::dist(t(sn)))
        names(out_snp[[i]][[1]]) <- method
        
      }
      else{
        # get tasks
        tasks <- get.task.list(y, snp_facets[i])
        out_snp[[i]] <- vector("list", nrow(tasks))
        names(out_snp[[i]]) <- tasks[,4]
        snp_columns <- unlist(strsplit(tasks[,3][1], "(?<!^)\\.", perl = T))
        
        # run the tasks
        for(j in 1:nrow(tasks)){
          t_snp_cols <- y@snp.meta[,snp_columns, drop = F]
          t_snp_indices <- do.call(paste, c(t_snp_cols, sep = "    "))
          t_snp_indices <- which(t_snp_indices %in% tasks[j,4])
          out_snp[[i]][[j]] <- list(stats::dist(t(sn[t_snp_indices,])))
          names(out_snp[[i]][[j]]) <- method
        }
      }
      
    }
  }

  #============return=========
  if(exists("out") & exists("out_snp")){
    out <- c(out, out_snp)
  }
  else if(exists("out_snp")){
    out <- out_snp
  }
  
  y <- update_calced_stats(y, all_facets, paste0("genetic_distance_", method))
  return(merge.snpR.stats(y, out, "genetic_distances"))
  
}

#' Calculate Isolation by Distance
#'
#' Calculates Isolation by Distance (IBD) for snpRdata objects by comparing the
#' genetic distance between samples or sets of samples to the geographic
#' distances between samples or sets of sapmles. IBD caluclated via a mantel
#' test.
#'
#' Genetic distance is caluclated via \code{\link{calc_genetic_distances}}.
#' Geographic distances are taken as-is for individual-individual comparisons
#' and by finding the geographic mean of a group of samples when sample level
#' facets are requested using the methods described by
#' \code{\link[geosphere]{geomean}}. Note that this means that if many samples
#' were collected from the same location, and those samples compose a single
#' level of a facet, the mean sampling location will be that single location.
#'
#' IBD is caluclated by comparing geographic and genetic distances using a
#' mantel test via \code{\link[ade4]{mantel.randtest}}.
#'
#' Note that geographic distance objects are also included in the returned
#' values.
#'
#'
#' @param x snpRdata object
#' @param facets character, default NULL. Facets over which to calculate IBD, as
#'   described in \code{\link{Facets_in_snpR}}.
#' @param x_y character, default c("x", "y"). Names of the columns containing
#'   geographic coordinates where samples were collected. There is no need to
#'   specify projection formats. The first should be longitude, the second
#'   should be latitude.
#' @param genetic_distance_method character, default "Edwards". The genetic
#'   distance method to use, see \code{\link{calc_genetic_distances}}.
#' @param ... Additional arguments passed to \code{\link[ade4]{mantel.randtest}}
#'
#' @return A snpRdata object with both geographic distance and IBD results
#'   merged into existing data.
#'
#' @author William Hemstrom
#'
#' @seealso \code{\link[ade4]{mantel.randtest}},
#'   \code{\link{calc_genetic_distances}}
#'
#' @export
#'
#' @examples # calculate ibd for several different facets
#' y <- stickSNPs
#' sample.meta(y) <- cbind(sample.meta(y), x = rnorm(ncol(y)), y = rnorm(ncol(y)))
#' y <-calc_isolation_by_distance(y, facets = c(".base", "pop", "pop.group","pop.group.fam"))
#' res <- get.snpR.stats(y, "pop.group", "ibd") # fetch results
#' res
#' plot(res$group.pop$groupV$Edwards) # plot perms vs observed
#' 
calc_isolation_by_distance <- function(x, facets = NULL, x_y = c("x", "y"), genetic_distance_method = "Edwards", ...){
  #================sanity checks=============================================
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  msg <- character()
  
  bad.xy <- which(!x_y %in% colnames(x@sample.meta))
  if(length(bad.xy) > 0){
    msg <- c(msg, paste0("Some coordinates (arument x_y) not found in sample metadata: ", paste(x_y[bad.xy])))
  }
  
  if(any(abs(x@sample.meta[,x_y[2]]) > 90)){
    msg <- c(msg, "Latitude values greater 90 or less than -90 detected. Note that the first value in 'x_y' must specify longitude, the second must specify latitude!")
  }
  
  pkg.check <- check.installed("geosphere")
  if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
  
  pkg.check <- check.installed("ade4")
  if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
  
  
  if(length(msg) > 0){
    stop(paste0(msg, collapse = "\n"))
  }
  
  #================prep and get genetic dist matrices========================
  facets <- check.snpR.facet.request(x, facets, "none")
  x <- add.facets.snpR.data(x, facets)
  
  # calc gds if needed and fetch
  cs <- check_calced_stats(x, facets, paste0("genetic_distance_", genetic_distance_method))
  missing <- which(!unlist(cs))
  if(length(missing) > 0){
    x <- calc_genetic_distances(x, names(cs)[missing])
  }
  gd <- .get.snpR.stats(x, facets, "genetic_distance")
  
  #===============fetch geo dist matrices=======================
  # get the geo dist matrices
  snp.levs <- check.snpR.facet.request(x, facets)
  geo_mats <- vector("list", length(snp.levs))
  names(geo_mats) <- snp.levs
  
  # get the geo matrices for each comparison. note that the same geo matrix is re-used across snp level facets.
  for(i in 1:length(geo_mats)){
    # get the categories for this
    categories <- get.task.list(x, names(geo_mats)[i])[,1:2, drop = F]
    
    # if the categories are all base, easy peasy
    if(all(categories[,1] == ".base")){
      geo_mats$.base[[genetic_distance_method]] <- stats::dist(x@sample.meta[,x_y])
    }
    else{
      # otherwise need to break it down by sample level facet
      # find geographic mean coordiantes for each facet level
      gmeans <- t(apply(categories, 1, function(row){
        matches <- fetch.sample.meta.matching.task.list(x, row)
        matches <- x@sample.meta[matches,]
        matches <- matches[,x_y]
        colnames(matches) <- c("x", "y")
        geosphere::geomean(matches)
      }))
      rownames(gmeans) <- categories[,2]
      
      #get the distances
      geo_mats[[categories[1,1]]][[genetic_distance_method]] <- stats::dist(gmeans)
    }
  }
  
  #===============calculate IBD==================
  ibd <- vector("list", length(gd))
  names(ibd) <- names(gd)
  # do the calcs: at each level, simply copy metadata in
  for(i in 1:length(gd)){
    ibd[[i]] <- vector("list", length(gd[[i]]))
    names(ibd[[i]]) <- names(gd[[i]])
    
    
    for(j in 1:length(gd[[i]])){
      ibd[[i]][[j]] <- vector("list", length(gd[[i]][[j]]))
      names(ibd[[i]][[j]]) <- names(gd[[i]][[j]])
      
      for(k in 1:length(gd[[i]][[j]])){
        
        tsampfacet <- check.snpR.facet.request(x, names(gd)[i])
        
        # quick check for NAs
        if(any(is.na(gd[[i]][[j]][[k]]))){
          # find bad entries
          matgd <- as.matrix(gd[[i]][[j]][[k]])
          bad.entries <- which(rowSums(is.na(matgd)) == ncol(matgd) - 1)
          
          # remove from gene dist
          gd[[i]][[j]][[k]] <- stats::as.dist(matgd[-bad.entries, -bad.entries])
          
          # remove from geo dist
          matgeod <- as.matrix(geo_mats[[tsampfacet]][[names(gd[[i]][[j]])[k]]])
          tgeodist<- stats::as.dist(matgeod[-bad.entries, -bad.entries])
          
          warning("NAs detected in genetic distance data for facet: ", paste(names(gd)[i], names(gd[[i]])[j], names(gd[[i]][[j]])[k]),
                  "\nlevels\t\n\t", paste0(row.names(matgeod)[bad.entries], collapse = "\n\t"))
        }
        else{
          tgeodist <- geo_mats[[tsampfacet]][[names(gd[[i]][[j]])[k]]]
        }
        
        ibd[[i]][[j]][[k]] <- ade4::mantel.randtest(gd[[i]][[j]][[k]], tgeodist, ...)
        
      }
    }
  }
  
  #===============return================================
  x <- merge.snpR.stats(x, geo_mats, "geo_dists")
  x <- merge.snpR.stats(x, ibd, "ibd")
  x <- update_calced_stats(x, facets, c("geo_dist", "ibd"))
  
  return(x)
}



#' Calculate weighted averages of previously calculated genetic statistics.
#'
#' Calculates a weighted average for a statistic, weighted by the number of
#' called genotypes at each locus. Works for single or pairwise statistics (pi,
#' ho, fst, etc.). Automatically calculates weighted statistic for every
#' previously calculated statistic.
#'
#' Weights are calculated using the equation \deqn{ M_{s} = \frac{\sum_{i =
#' 1}^{n} s_{i} * w_{i}}{\sum_{i = 1}^{n} w_{i}}} Where\eqn{n} is the number of
#' SNPs, \eqn{s_{i}} is the value of the statistic in SNP \eqn{i}, and
#' \eqn{w_{i}} is the number of times that SNP was genotyped. Note that this
#' will correct for a range in sequencing depth within samples, but does not
#' really correct for variance in sequencing depth between populations or other
#' facet levels.
#'
#' @param x snpRdata object.
#' @param facets character, default NULL. Facets for which to calculate weighted
#'   stats (broken down by category). See \code{\link{Facets_in_snpR}} for
#'   details.
#' @param type character, default "single". Type of statistic to weight:
#'   \itemize{\item{single: } Statistics calculated in a single subfacet, such
#'   as pi. \item{pairwise: } Statistics calculated pairwise between subfacets,
#'   such as Fst. }
#'
#' @return A snpR data object with weighted statistics merged in, accessable via
#'   \code{\link{get.snpR.stats}} using type = "pop".
#'
#' @export
#' @examples
#' # single
#' x <- calc_basic_snp_stats(stickSNPs, "pop")
#' x <- calc_weighted_stats(x, "pop")
#' get.snpR.stats(x, "pop", "pop")
#'
#' # pairwise
#' x <- calc_weighted_stats(x, "pop", type = "pairwise") # fst calculated in last step
#' get.snpR.stats(x, "pop", "pairwise")
#' 
calc_weighted_stats <- function(x, facets = NULL, type = "single", stats_to_get){
  #===========sanity checks===============
  msg <- character(0)
  if(!is.snpRdata(x)){
    msg <- c(msg, "x is not a snpRdata object.\n")
  }
  
  if(!type %in%  c("pairwise", "single", "single.window", "sample")){
    msg <- c(msg, "Type unsupported. Options: pairwise, single, single.window, sample.\n")
  }
  if(length(msg) > 0){
    stop(msg)
  }
  
  #===========calculate weighted stats======
  facets <- check.snpR.facet.request(x, facets, "none")
  x <- add.facets.snpR.data(x, facets)
  calced <- check_calced_stats(x, facets, "maf")
  calcedb <- unlist(calced)
  names(calcedb) <- names(calced)
  if(sum(calcedb) != length(calcedb)){
    x <- calc_maf(x, facets = names(calcedb)[which(!calcedb)])
  }
  
  
  stats <- .get.snpR.stats(x, facets, type)
  facets <- check.snpR.facet.request(x, facets, "none", T)
  if(any(facets[[2]] == "complex") & type %in% c("single.window")){
    facets[[1]] <- c(facets[[1]], facets[[1]][which(facets[[2]] == "complex")])
    facets[[2]] <- c(facets[[2]], rep("special", sum(facets[[2]] == "complex")))
  }
  
  for(i in 1:length(facets[[1]])){
    split.part <- unlist(strsplit(facets[[1]][i], split = "(?<!^)\\.", perl = T))
    split.part <- check.snpR.facet.request(x, split.part, remove.type = "none", TRUE)
    snp.part <- split.part[[1]][which(split.part[[2]] == "snp")]
    split.snp.part <- unlist(strsplit(split.part[[1]][which(split.part[[2]] == "snp")], split = "(?<!^)\\.", perl = T))
    samp.part <- split.part[[1]][which(split.part[[2]] == "sample")]
    split.samp.part <- unlist(strsplit(split.part[[1]][which(split.part[[2]] == "sample")], split = "(?<!^)\\.", perl = T))
    
    
    
    # get weights and stats
    #=======================simple case==================
    if(type == "single"){
      weights <- (stats$maj.count + stats$min.count)/2
      
      if(facets[[2]][i] == "complex"){
        split.part <- unlist(strsplit(facets[[1]][i], split = "(?<!^)\\.", perl = T))
        split.part <- check.snpR.facet.request(x, split.part, remove.type = "none", TRUE)
        if(length(split.part[[1]]) > 2){
          snp.partp <- paste(snp.part, collapse = ".")
          samp.partp <- paste(samp.part, collapse = ".")
          split.part <- list(c(snp.partp, samp.partp), c("snp", "sample"))
        }
        keep.rows <- which(stats$facet %in% split.part[[1]][which(split.part[[2]] == "sample")])
        
        group_key <- c("facet", "subfacet", split.snp.part)
      }
      else if(facets[[2]][i] == "sample"){
        keep.rows <- which(stats$facet %in% facets[[1]][i])
        group_key <- c("facet", "subfacet")
      }
      else if(facets[[2]][i] == "snp"){
        keep.rows <- which(stats$facet == ".base")
        group_key <- split.snp.part
      }
      else{
        keep.rows <- which(stats$facet == ".base")
        group_key <- "facet"
      }
      
      
      selected_stats <- stats[keep.rows, stats_to_get, drop = F]
      if(ncol(selected_stats) == 0){
        stop("No calculated stats to weight.\n")
      }
      weights <- weights[keep.rows]
      
      group_key_tab <- stats[keep.rows,group_key, drop = F]
      group_key_tab$key <- do.call(paste, group_key_tab)
    }
    #=====================single window case======================
    else if(type == "single.window"){
      weights <- stats$n_snps
      
      
      if(facets[[2]][i] == "complex" | facets[[2]][i] == "special"){
        if(length(split.part[[1]]) > 2){
          snp.partp <- paste(snp.part, collapse = ".")
          samp.partp <- paste(samp.part, collapse = ".")
          split.part <- list(c(snp.partp, samp.partp), c("snp", "sample"))
        }
        keep.rows <- which(stats$facet %in% split.part[[1]][which(split.part[[2]] == "sample")] &
                             stats$snp.facet %in% split.part[[1]][which(split.part[[2]] == "snp")])
        
        if(facets[[2]][i] == "complex"){
          group_key <- c("facet", "subfacet", "snp.facet", "snp.subfacet")
        }
        else{
          group_key <- c("facet", "subfacet")
        }
      }
      else if(facets[[2]][i] == "sample"){
        keep.rows <- which(stats$facet %in% facets[[1]][i] & stats$snp.facet == ".base")
        group_key <- c("facet", "subfacet")
      }
      else if(facets[[2]][i] == "snp"){
        keep.rows <- which(stats$facet == ".base" & stats$snp.facet == facets[[1]][i])
        group_key <- c("snp.facet", "snp.subfacet")
      }
      else{
        keep.rows <- which(stats$facet == ".base" & stats$snp.facet == ".base")
        group_key <- "facet"
      }
      
      
      selected_stats <- stats[keep.rows, stats_to_get, drop = F]
      if(ncol(selected_stats) == 0){
        stop("No calculated stats to weight.\n")
      }
      weights <- weights[keep.rows]
      
      group_key_tab <- stats[keep.rows,group_key, drop = F]
      group_key_tab$key <- do.call(paste, group_key_tab)
    }
    #=======================pairwise case==========================
    else if(type == "pairwise"){
      weights <- stats$nk
      

      if(facets[[2]][i] == "complex"){
        if(length(split.part[[1]]) > 2){
          snp.partp <- paste(snp.part, collapse = ".")
          samp.partp <- paste(samp.part, collapse = ".")
          split.part <- list(c(snp.partp, samp.partp), c("snp", "sample"))
        }
        keep.rows <- which(stats$facet %in% split.part[[1]][which(split.part[[2]] == "sample")])
        
        group_key <- c("facet", "comparison", split.part[[1]][which(split.part[[2]] == "snp")])
      }
      else if(facets[[2]][i] == "sample"){
        keep.rows <- which(stats$facet %in% facets[[1]][i])
        group_key <- c("facet", "comparison")
      }
      else{
        stop("Cannot calculate pairwise stats for non-sample facets--check what is being passed as facets to calc_weighted_averages.")
      }
      
      
      selected_stats <- stats[keep.rows, stats_to_get, drop = F]
      if(ncol(selected_stats) == 0){
        stop("No calculated stats to weight.\n")
      }
      weights <- weights[keep.rows]
      
      group_key_tab <- stats[keep.rows,group_key, drop = F]
      group_key_tab$key <- do.call(paste, group_key_tab)

    }
    #=======================sample case========================
    else if(type == "sample"){
      if(i != 1){
        if(facets[[2]][i - 1] == "complex"){
          stats <- ostats
          rm(ostats); gc()
        }
      }
      # get the per sample missingness
      if(facets[[2]][i] == "sample"){
        # per sample, weights are easy
        weights <- nrow(x) - matrixStats::colSums2(genotypes(x) == x@mDat)
        keep.rows <- which(stats$snp.facet == ".base")
        
        group_key <- split.samp.part
        group_key_tab <- stats[keep.rows, group_key, drop = F]
        group_key_tab$key <- do.call(paste, c(group_key_tab, sep = "."))
        colnames(group_key_tab)[1] <- "subfacet"
        group_key_tab$subfacet <- group_key_tab$key
        group_key_tab$facet <- paste0(samp.part, collapse = ".")
      }
      if(facets[[2]][i] == "snp"){
        keep.rows <- which(stats$snp.facet == snp.part & stats$facet == ".base")
        group_key <- "snp.subfacet"
        group_key_tab <- stats[keep.rows, group_key, drop = F]
        group_key_tab$key <- group_key_tab$snp.subfacet
        
        
        tl <- get.task.list(x, facets[[1]][i])
        weights <- numeric(nrow(group_key_tab))
        for(j in 1:nrow(tl)){
          snps_in_set <- which(unlist(do.call(paste, c(snp.meta(x)[,split.snp.part, drop = F], sep = "."))) ==
                                        tl[j,4])
          tweights <- length(snps_in_set) - matrixStats::colSums2(genotypes(x)[snps_in_set,] == x@mDat)
          weights[which(group_key_tab$key == tl[j,4])] <- tweights
        }
        group_key_tab$snp.facet <- snp.part
      }
      if(facets[[2]][i] == "complex"){
        make_group_key_tab <- function(stats, keep.rows, group_key){
          group_key_tab <- stats[keep.rows, group_key, drop = F]
          group_key_tab$key <- unlist(do.call(paste, c(group_key_tab, sep = ".")))
          group_key_tab$snp.facet <- snp.part
          group_key_tab$facet <- samp.part
          group_key_tab$subfacet <- unlist(do.call(paste, c(group_key_tab[,split.samp.part, drop = F], sep = ".")))
          return(group_key_tab)
        }
        group_key <- c("snp.subfacet", split.samp.part)
        keep.rows <- which(stats$snp.facet == snp.part & stats$facet == ".base")
        group_key_tab <- make_group_key_tab(stats, keep.rows, group_key)
        
        
        
        tl <- get.task.list(x, facets[[1]][i])
        weights <- stats[keep.rows,]
        for(j in 1:nrow(tl)){
          
          snps_in_set <- which(unlist(do.call(paste, c(snp.meta(x)[,split.snp.part, drop = F], sep = "."))) ==
                                 tl[j,4])
          samps_in_set <- fetch.sample.meta.matching.task.list(x, tl[j,])
          if(length(snps_in_set) > 0 & length(samps_in_set) > 0){
            tweights <- length(snps_in_set) - matrixStats::colSums2(genotypes(x)[snps_in_set, samps_in_set, drop = F] == x@mDat)
            m <- data.table(weights = tweights, facet = tl[j,1], subfacet = tl[j,2], snp.facet = tl[j,3], snp.subfacet = tl[j,4])
            m <- cbind(m, sample.meta(x)[samps_in_set,])
            mn <- c("snp.subfacet", "snp.subfacet", colnames(x@sample.meta))
            weights <- smart.merge(m, weights, mn, mn)
          }
        }
        
        ostats <- stats
        stats <- weights
        keep.rows <- 1:nrow(stats)
        weights <- weights$weights
        stats <- as.data.frame(stats)
        group_key_tab <- make_group_key_tab(stats, keep.rows, group_key)
      }
      if(facets[[2]][i] == ".base"){
        weights <- nrow(x) - matrixStats::colSums2(genotypes(x) == x@mDat)
        keep.rows <- which(stats$snp.facet == ".base" & stats$facet == ".base")
        group_key <- c("facet", "subfacet", "snp.facet", "snp.subfacet")
        group_key_tab <- stats[keep.rows, group_key, drop = F]
        group_key_tab$key <- do.call(paste, c(group_key_tab, sep = "."))
      }
      
      selected_stats <- stats[keep.rows, stats_to_get, drop = F]
    }
    else{stop("calc_weighted_stats type not recognized.")}
    
    #===================calculate weighted means using equation sum(w*s)/sum(w)=======================
    weighted <- as.data.table(weights*selected_stats)
    group_mean_weights <- tapply(weights, group_key_tab$key, sum, na.rm = T)
    weighted$key <- group_key_tab$key
    means <- weighted[,lapply(.SD, sum, na.rm = T), by = key] # lapply isn't ideal, but data.table is quite fast anyway. Could consider other options
    set(means, j = 2:ncol(means), value = means[,-1]/group_mean_weights[match(means$key, names(group_mean_weights))])

    # merge and return
    ## get the stats back in a format with facet and sub-facet, clean up, and return
    mstats <- merge(means, unique(group_key_tab), by = "key")
    drop_col <- which(colnames(mstats) == "key")
    mstats <- mstats[,-..drop_col]
    new.ord <- c(2:(ncol(selected_stats) + ncol(group_key_tab) - 1), 1:ncol(selected_stats))
    mstats <- mstats[,..new.ord]
    colnames(mstats)[(ncol(group_key_tab)):ncol(mstats)] <- paste0("weighted_mean_", colnames(mstats)[(ncol(group_key_tab)):ncol(mstats)])
    
    
    
    #====================add on extra filler columns====================
    if(!"snp.subfacet" %in% colnames(mstats)){
      if(facets[[2]][i] %in% c(".base", "sample")){
        mstats <- mstats[,snp.subfacet := ".base"]
      }
      else{
        if(facets[[2]][i] == "special"){
          mstats$snp.subfacet <- ".base"
        }
        else if(length(split.snp.part) > 1){
          mstats$snp.subfacet <- do.call(paste, c(mstats[,..split.snp.part], sep = "."))
        }
        else{
          mstats$snp.subfacet <- mstats[,..split.snp.part]
        }
      }
    }
    if(!"snp.facet" %in% colnames(mstats)){
      if(facets[[2]][i] %in%  "special"){
        snp.partp <- paste(snp.part, collapse = ".")
        mstats <- mstats[,snp.facet := snp.partp]
      }
      else if(facets[[2]][i] %in% c(".base", "sample")){
        mstats <- mstats[,snp.facet := ".base"]
      }
      else{
        mstats <- mstats[,snp.facet := check.snpR.facet.request(x, paste0(split.snp.part, collapse = "."), remove.type = "sample")]
      }
    }
    if(!"subfacet" %in% colnames(mstats)){
      if("comparison" %in% colnames(mstats)){
        colnames(mstats)[which(colnames(mstats) == "comparison")] <- "subfacet"
      }
      else{
        mstats <- mstats[,subfacet := ".base"]
      }
    }
    if(!"facet" %in% colnames(mstats)){
      mstats <- mstats[,facet := ".base"]
    }
    
    good.cols <- c("facet", "subfacet", "snp.facet", "snp.subfacet", paste0("weighted_mean_", stats_to_get))
    mstats <- mstats[,..good.cols]
    mstats[,1:4] <- dplyr::mutate_all(mstats[,1:4], as.character)
    x <- merge.snpR.stats(x, mstats, type = "weighted.means")
  }
  
  
  
  return(x)
}




