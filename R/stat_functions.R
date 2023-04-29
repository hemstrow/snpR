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
#'@section \eqn{\pi}:
#'
#'  Calculates \eqn{\pi} (nucleotide diversity/average number of pairwise
#'  differences) according to Hohenlohe et al. (2010).
#'
#'@section \ifelse{html}{\out{H<sub>E</sub>}}{\eqn{H_E}}: 
#'  
#'  Calculates traditional
#'  expected heterozygosity \eqn{2pq}. Note that this will produce results
#'  almost identical to \eqn{\pi}.
#'
#'@section \ifelse{html}{\out{H<sub>O</sub>}}{\eqn{H_O}}:
#'
#'  Calculates observed heterozygosity.
#'
#'@section maf:
#'
#'  Calculates minor allele frequencies and note identities and counts of major
#'  and minor alleles.
#'
#'  
#'@section private alleles:
#'
#'  Determines if each SNP is a private allele across all levels in each sample
#'  facet. Will return an error if no sample facets are provided.
#'
#'@section hwe:
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
#'@param fwe_method character, default "BY". Type of Family-Wise Error
#'  correction (multiple testing correction) to use. For details and options,
#'  see \code{\link[stats]{p.adjust}}. If no correction is desired, set this argument
#'  to "none".
#'@param fwe_case character, default c("by_facet", "by_subfacet", "overall").
#'  How should Family-Wise Error correction (multiple testing correction) be
#'  applied? \itemize{\item{"by_facet":} Each facet supplied (such as pop or
#'  pop.fam) is treated as a set of tests. \item{"by_subfacet":} Each level of
#'  each subfacet is treated as a separate set of tests. \item{"overall":} All
#'  tests are treated as a set.}
#'
#'@aliases calc_pi calc_hwe calc_ho calc_private calc_maf calc_he
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
#' \dontrun{
#' x <- calc_hwe(stickSNPs, facets = c("pop", "pop.fam"))
#' get.snpR.stats(x, c("pop", "pop.fam"))
#' }
NULL


#'Facets in snpR
#'
#'Facets are used to describe the ways in which data should be broken down/split
#'up for analysis.
#'
#'Facet designation follows a specific format. Facets are given as a character
#'vector, where each entry designates one unique combination of levels over
#'which to separate the data. Levels within a facet are separated by a '.'. Each
#'facet can contain multiple snp and/or sample levels. Multiple facets can be
#'run with a single line of code.
#'
#'For example, c("chr.pop", "pop") would split the data first by both chr
#'and pop and then by pop alone. This will produce the same result as running
#'the function with both "chr.pop" and then again with "pop", although the
#'former is typically more computationally efficient.
#'
#'If multiple sample or snp levels are provided in a single facet, the data is
#'simultaneously broken up by \emph{both} levels. For example, the facet
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
#'The base facet--that is, the entire data with no categorical divisions--can be
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
#' \dontrun{
#' # multiple facets
#' x <- calc_pi(stickSNPs, facets = c("pop", "pop.fam"))
#' get.snpR.stats(x, c("pop", "pop.fam"))
#' }
NULL


#'@export
#'@describeIn calc_single_stats \eqn{\pi} (nucleotide diversity/average number of pairwise differences)
calc_pi <- function(x, facets = NULL){
  func <- function(x){
    bins <- choose(x$as, 2)
    tot <- choose(rowSums(x$as), 2)
    return(1 - (rowSums(bins))/tot)
  }
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  

  # add any missing facets
  ofacets <- facets
  facets <- .check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    .make_it_quiet(x <- .add.facets.snpR.data(x, facets))
  }

  out <- .apply.snpR.facets(x, facets, "gs", func, case = "ps")
  colnames(out)[ncol(out)] <- "pi"
  x <- .merge.snpR.stats(x, out)
  x <- .calc_weighted_stats(x, ofacets, type = "single", "pi")
  x <- .update_calced_stats(x, facets, "pi", "snp")
  x <- .update_citations(x, "Hohenlohe2010", "pi", "pi, number of pairwise differences")
  
  return(x)
}


#'@export
#'@describeIn calc_single_stats minor allele frequency
calc_maf <- function(x, facets = NULL){
  maj_count <- NULL
  
  
  # function to run on whatever desired facets (now an internal--.maf_func):
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  # add any missing facets
  facets <- .check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    .make_it_quiet(x <- .add.facets.snpR.data(x, facets))
  }
  
  if(facets[1] == ".base" & length(facets) == 1){
    out <- .apply.snpR.facets(x,
                             facets = facets[1],
                             req = "gs",
                             fun = .maf_func,
                             case = "ps",
                             m.al = substr(x@mDat,1, nchar(x@mDat)/2))
  }
  else{
    # calculate the base facet if not yet added
    logi <- .check_calced_stats(x, ".base", "maf")
    if(!logi[[".base"]]["maf"]){
      x <- calc_maf(x)
    }
    
    major_minor_base <- .get.snpR.stats(x)[,c("major", "minor")]
    out <- .apply.snpR.facets(x,
                             facets = facets,
                             req = "gs",
                             fun = .maf_func,
                             case = "ps",
                             m.al = substr(x@mDat,1, nchar(x@mDat)/2),
                             ref = major_minor_base)
    
  }
  if(!.is.bi_allelic(x)){
    out$minor <- "NA"
  }
  
  x <- .update_calced_stats(x, facets, "maf", "snp")
  return(.merge.snpR.stats(x, out))
}

#'Tajima's D from SNP data.
#'
#'\code{Tajimas_D} calculates Tajima's theta/pi, Watterson's theta, and Tajima's
#'D over a sliding window.
#'
#'Tajima's D compares estimates of theta based on either the number of observed
#'pairwise differences (Tajima's theta) and the number of substitutions vs
#'expected total tree length (Watterson's Theta). Since low frequency minor
#'variants contribute to these statistics and they rely on the ratio of the
#'number of variants vs the number of sequenced non-polymorphic sites, this
#'function should only be run on data that is \emph{unfiltered} aside from the
#'removal of poorly sequenced bases, etc.
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
#'  window will include all SNPs within 3*sigma or sigma kilobases depending on
#'  the \code{triple_sigma} argument. If either sigma or step are NULL, the
#'  entire snp subfacet will be done at once (for example, the whole
#'  chromosome).
#'@param step numeric or NULL, default \code{2*sigma} (non-overlapping windows). Number
#'  of bases to move between each window, in kilobases. If either sigma or step
#'  are NULL, the entire snp subfacet will be done at once (for example, the
#'  whole chromosome).
#'@param par numeric or FALSE, default FALSE. If numeric, the number of cores to
#'  use for parallel processing.
#'@param triple_sigma Logical, default TRUE. If TRUE, sigma will be tripled to
#'  create windows of 6*sigma total.
#'@param verbose Logical, default FALSE. If TRUE progress will be printed to the
#'  console.
#'
#'@return snpRdata object, with Watterson's Theta, Tajima's Theta, and Tajima's
#'  D for each window merged in to the window.stats slot.
#'
#' @examples
#' # slow, so not run
#' \dontrun{
#' # broken by population, windows across linkage group
#' x <- calc_tajimas_d(stickSNPs, facets = "chr.pop", sigma = 200, step = 50)
#' get.snpR.stats(x, "chr.pop", "tajimas_d")
#'
#' # the entire population at once, note that sigma and step are NULL and
#' # no chromosome/linkage group/scaffold/etc set.
#' # this will calculate overall tajima's D without a window for each population.
#' x <- calc_tajimas_d(stickSNPs, facets = "pop")
#' get.snpR.stats(x, "pop", "tajimas_d")
#'
#' # for the overall dataset, note that sigma and step are NULL
#' # this will calculate overall tajima's D for each chr/pop
#' x <- calc_tajimas_d(stickSNPs, facets = "chr.pop")
#' get.snpR.stats(x, "pop.chr", "tajimas_d")
#' }
#'@export
#'@references Tajima, F. (1989). \emph{Genetics}
#'@author William Hemstrom
calc_tajimas_d <- function(x, facets = NULL, sigma = NULL, step = 2*sigma, par = FALSE,
                           triple_sigma = TRUE,
                           verbose = FALSE){
  #===============sanity checks==========================
  if(!is.snpRdata(x)){
    stop("x must be a snpRdata object.")
  }
  
  step <- eval(step) # forces this to eval before we change sigma.
  
  if(!is.null(sigma) & !is.null(step)){
    .sanity_check_window(x, sigma, step, stats.type = "single", nk = TRUE, facets = facets)
  }
  else if(is.null(step) & !is.null(sigma)){
    sigma <- NULL
    step <- NULL
    warning("If no step size is provided but sigma is, the entire facet level will be considerd as a single window (sigma will be ignored).\n")
  }
  else{
    sigma <- NULL
    step <- NULL
  }
  if(is.null(facets[1]) | ".base" %in% facets | "sample" %in% .check.snpR.facet.request(x, facets, "none", TRUE)){
    warning("Tajima's D has little meaning of snps on different chromosomes are considered together. Consider adding a snp level facet.")
  }
  
  
  #=================subfunction=========
  func <- function(as, par, sigma, step, report){
    out <- data.frame(position = numeric(0), sigma = numeric(0), ws.theta = numeric(0), ts.theta = numeric(0), D = numeric(0), n_snps = numeric(0)) #initialize output
    tps <- sort(as$position) #get the site positions, sort
    lsp <- tps[length(tps)] #get the position of the last site to use as endpoint
    c <- 0 #set starting position
    i <- 1 #set starting iteration for writing output
    
    # if sigma is null, set to run everything in one window
    if(is.null(sigma)){
      sigma <- range(as$position)
      c <- mean(sigma)
      step <- c + 1
      sigma <- sigma[2] - sigma[1]
    }
    else{
      sigma <- 1000*sigma
    }
    
    
    # run the loop
    while (c <= lsp){
      
      start <- c - ifelse(triple_sigma, sigma*3, sigma) #change window start
      end <- c + ifelse(triple_sigma, sigma*3, sigma) #change window end
      
      # take all the snps in the window, calculate T's theta, W's theta, and T's D
      wsnps <- as[as$position <= end & as$position >= start,] # get only the sites in the window
      
      ac_cols <- which(colnames(wsnps) %in% colnames(x@geno.tables$as))
      
      wsnps <- wsnps[!(rowSums(wsnps[,ac_cols]) == 0),] # remove any sites that are not sequenced in this pop/group/whatever
      if(nrow(wsnps) == 0){ #if no snps in window
        c <- c + step*1000 #step along
        next #go to the next window
      }
      
      # binomials for each allele
      btop <- 0
      for(j in ac_cols){
        btop <- btop + choose(wsnps[,j],2) # binomial for this allele
      }
      
      ntotal <- rowSums(wsnps[,colnames(x@geno.tables$as)])
      
      bts <- choose(ntotal,2) #binomial for all alleles
      ts.theta <- sum(1-(btop/bts)) #average number of pairwise differences (pi) per snp. Equivalent to sum((ndif/ncomp)) for all snps
      #ts.thetaf <- ts.theta/nrow(wsnps) #pi for the whole region, which includes all of the non-polymorphic sites. Reason why this shouldn't be run with a maf filter, probably.
      n_seg <- nrow(wsnps[which(rowSums(wsnps[,ac_cols] != 0) > 1),]) #number of segregating sites
      K <- round(mean(ntotal)) #average sample size for ws.theta, as in hohenlohe 2010. Alternative would make this into a function, then use lapply on the vector of K values
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
      out[i,"num_seg"] <- n_seg
      out[i,"D"] <- D
      out[i, "n_snps"] <- nrow(wsnps)
      out[i, "start"] <- start
      out[i, "end"] <- end
      c <- c + step*1000
      i <- i + 1
    }
    if(nrow(out) > 0){
      out$sigma <- sigma/1000
      out$step <- step
      out$nk.status <- FALSE
      out$triple_sigma <- triple_sigma
      out$gaussian <- FALSE
    }
    return(out)
  }
  
  
  #=================prep==========
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  # add any missing facets
  add.facets <- .check.snpR.facet.request(x, facets)
  if(!all(add.facets %in% x@facets)){
    .make_it_quiet(x <- .add.facets.snpR.data(x, add.facets))
  }
  facets <- .check.snpR.facet.request(x, facets, remove.type = "none")
  
  #=============run=============
  out <- .apply.snpR.facets(x,
                            facets = facets,
                            req = "meta.as",
                            fun = func,
                            case = "ps.pf.psf",
                            par = par,
                            sigma = sigma,
                            step = step,
                            verbose = verbose)
  
  #===========merge and clean============
  x <- .merge.snpR.stats(x, out, type = "window.stats")
  x <- .update_calced_stats(x, facets, "tajimas_d")
  x <- .calc_weighted_stats(x, facets, type = "single.window", c("ws.theta", "ts.theta", "num_seg", "D"))
  x <- .update_citations(x, "Tajima1989", "Tajima's_d", "Tajima's D, as well as Watterson's and Tajima's Theta")
  
  # calc weights ignoring any snp levels (for stuff like overall population means)
  samp.facets <- .check.snpR.facet.request(x, facets)
  return(x)
  
}




#'Pairwise FST from SNP data.
#'
#'\code{calc_pairwise_fst} calculates pairwise FST for each SNP for each
#'possible pairwise combination of populations.
#'
#'Calculates FST according to either Weir and Cockerham 1984 or using the
#'\code{\link[genepop]{Fst}} function from the genepop package (see references).
#'
#'If the genepop option is used, several intermediate files will be created in
#'the default temporary directory (see \code{\link{tempfile}}).
#'
#'The Weir and Cockerham (1984) and genepop methods tend to produce very similar
#'results both per SNP and per population. Generally, the former option may be
#'preferred for computational efficiency.
#'
#'P-values for group level comparisons can be calculated via bootstrapping using
#'the boot option. Bootstraps are performed via randomly mixing individuals
#'amongst different levels of the supplied facet, and thus the null hypothesis
#'is that all groups are panmictic. P-values are calculated according to
#'\code{\link[ade4]{randtest}}, although that function is not directly called.
#'
#'The data can be broken up categorically by either SNP and/or sample metadata,
#'as described in \code{\link{Facets_in_snpR}}. Since this is a pairwise
#'statistic, at least a single sample level facet must be provided.
#'
#'Method Options: \itemize{ \item{"wc": }{Weir and Cockerham 1984.}
#'item{"Genepop": }{As used in genepop, Rousset
#'2008.}}
#'
#'@param x snpRdata. Input SNP data.
#'@param facets character. Categorical metadata variables by which to break up
#'  analysis. See \code{\link{Facets_in_snpR}} for more details.
#'@param method character, default "wc". Defines the FST estimator to use.
#'   Options: \itemize{ \item{wc: } Weir and Cockerham (1984).
#'   \item{genepop: } Rousset (2008), uses the genepop package. }
#'@param boot numeric or FALSE, default FALSE. The number of bootstraps to do.
#'  See details.
#'@param boot_par numeric or FALSE, default FALSE. If a number, bootstraps will
#'  be processed in parallel using the supplied number of cores.
#'@param zfst logical, default FALSE. If TRUE, z-distributed Fst scores (zFST)
#'  will be calculated, equal to (fst - mean(fst))/sd(fst) within each group.
#'  The resulting values will be in the column "zfst", accessible using the
#'  usual \code{\link{get.snpR.stats}} method.
#'@param fst_over_one_minus_fst logical, default FALSE. If TRUE, fst/(1-fst)
#'  will be calculated, and will be in the column "fst_id" accessible using the
#'  usual \code{\link{get.snpR.stats}} method.
#'@param keep_components logical, default FALSE. If TRUE, the variance
#'  components "a", "b", and "c" will be held and accessible from the
#'  \code{$pairwise} element (named "var_comp_a", "var_comp_b", and
#'  "var_comp_c", respectively) using the usual \code{\link{get.snpR.stats}}
#'  method. This may be useful if working with very large datasets that need to
#'  be run with separate objects for each chromosome, etc. for memory purposes.
#'  Weighted averages can be generated identically to those from snpR by taking
#'  the weighted mean (via the \code{\link[stats]{weighted.mean}}) of "a"
#'  divided by the sum of the weighted means of "a", "b", and "c" using the
#'  number of SNPs called in a comparison (returned in the "nk" column from
#'  \code{\link{get.snpR.stats}}) as weights within each population comparison.
#'  Note that this is different than taking the weighted mean of a/(a + b + c)!
#'@param cleanup logical, default TRUE. If TRUE, any new files created during
#'  FST calculation will be automatically removed.
#'@param verbose Logical, default FALSE. If TRUE, some progress updates will be
#'  reported.
#'
#'@return A snpRdata object with pairwise FST as well as the number of total
#'  observations at each SNP in each comparison merged in to the pairwise.stats
#'  slot.
#'
#'@references Weir and Cockerham (1984). \emph{Evolution}

#'@references Weir (1990). Genetic data analysis. Sinauer, Sunderland, MA
#'@references Rousset (2008). \emph{Molecular Ecology Resources}
#'
#'@author William Hemstrom
#'@export
#'
#' @examples
#' # Using Weir and Cockerham 1984's method
#' x <- calc_pairwise_fst(stickSNPs, "pop")
#' get.snpR.stats(x, "pop", "fst")
#'
#' \dontrun{
#' # Using genepop
#' x <- calc_pairwise_fst(stickSNPs, "pop", "genepop")
#' get.snpR.stats(x, "pop", "fst")
#' 
#' # bootstrap p-values for overall pairwise-Fst values
#' x <- calc_pairwise_fst(stickSNPs, "pop", boot = 5)
#' get.snpR.stats(x, "pop", "fst")
#' }
calc_pairwise_fst <- function(x, facets, method = "wc", boot = FALSE, 
                              boot_par = FALSE,
                              zfst = FALSE,
                              fst_over_one_minus_fst = FALSE,
                              keep_components = FALSE,
                              cleanup = TRUE,
                              verbose = FALSE){
  facet <- subfacet <- .snp.id <-  weighted.mean <- nk <- fst <- comparison <- ..meta.cols <- ..meta_colnames <- ..ac_cols <- ..col.ord <- fst_id <- . <- ..gc_cols <- ..het_cols_containing_k <- NULL
  a <- b <- ..component_cols <- NULL
  
  if(!isTRUE(verbose)){
    cat <- function(...){}
  }
  
  #============================sanity and facet checks========================
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  if(method == "genepop"){
    .check.installed("genepop")
  }
  
  method <- tolower(method)
  if(!method %in% c("genepop", "wc")){
    stop("Method not found. Acceptable methods: genepop, wc.")
  }
  
  # add any missing facets
  ofacets <- facets
  facets <- .check.snpR.facet.request(x, facets, return.type = T, fill_with_base = F)
  if(any(facets[[2]] == ".base")){
    stop("At least one sample level facet is required for pairwise Fst estimation.")
  }
  facets <- facets[[1]]
  if(!all(facets %in% x@facets)){
    .make_it_quiet(x <- .add.facets.snpR.data(x, facets))
  }
  
  
  method <- tolower(method)
  
  bi_allelic <- .is.bi_allelic(x)
  
  
  #============================subfunctions=========================
  func <- function(x, method, facets = NULL, g.filename = NULL, ac_cols = NULL, meta_colnames = NULL, gc_cols = NULL){
    
    if(method != "genepop"){
      x <- data.table::as.data.table(x)
      data.table::setkey(x, subfacet, .snp.id)
      pops <- sort(unique(x$subfacet))
      npops <- length(pops)
      pnk.length <- (npops*(npops-1))/2 * (nrow(x)/npops)
    }
    else{
      pops <- sort(unique(x@facet.meta[which(x@facet.meta$facet == facets), "subfacet"]))
      npops <- length(pops)
      pnk.length <- (npops*(npops-1))/2 * (nrow(x))
    }
    
    pnk <- data.table::data.table(comparison = character(pnk.length), ntotal = integer(pnk.length))
    
    #===============genepop======================
    if(method == "genepop"){
      if(dir.exists(g.filename)){
        o.dir <- getwd()
        setwd(g.filename)
        g.filename <- list.files(pattern = "genepop.txt")
      }
      .make_it_quiet(genepop::Fst(g.filename, pairs = TRUE))
      
      
      
      
      
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
        tvals <- unlist(strsplit(tvals, " +")) #split and unlist the values
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
                                       .snp.id = x@facet.meta$.snp.id[x@facet.meta$facet == facets])
      
      n_tots$nk <- rowSums(x@geno.tables$as[which(x@facet.meta$facet == facets),])
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
      
      out$loci <- cbind(out$loci, data.table::as.data.table(x@snp.meta))
      
      out$loci$n_total <- pnk$ntotal
      if(length(out$overall) == 1){
        names(out$overall) <- paste0(pnames[1], "~", pnames[2])
      }
      
      # clean and return, we're done.
      cat("Finished.\n")
      if(exists("o.dir")){
        setwd(o.dir)
      }
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
    
    out <- list(Fst = out, a = as.data.table(out), b = as.data.table(out), c = as.data.table(out)) # not wrapping in as.data.table results in each object occupying the same spot in memory and thus overwriting eachother, which is not what we want!
    
    #loop through each comparison and calculate pairwise FST at each site
    c.col <- 1L
    prog <- 1L
    
    for (i in 1:(length(pops) - 1)){ #i is the first pop
      idat <- x[subfacet == pops[i]] # get data for first pop
      j <- i + 1 #initialize j as the next pop
      for (j in j:length(pops)){#j is pop being compared
        jdat <- x[subfacet == pops[j]] #get data for second pop
        
        if(method == "wc"){
          
          
          nt1 <- .fix..call(rowSums(idat[,..ac_cols]))
          nt2 <- .fix..call(rowSums(jdat[,..ac_cols]))
          ps1_f <- .fix..call(idat[,..ac_cols]/nt1)
          ps2_f <- .fix..call(jdat[,..ac_cols]/nt2)
          intot <- nt1/2
          jntot <- nt2/2
          
          # compute variance parts
          r <- 2 # number of comps
          nbar <- (intot + jntot)/2 #average sample size in individuals
          comb_ntots <- cbind(intot, jntot)
          CV <- matrixStats::rowSds(comb_ntots)/rowMeans(comb_ntots) # coefficient of variation in sample size
          nc <- nbar*(1-(CV^2)/r)
          parts <- vector("list", ncol(ps1_f))
          gc_colnames <- colnames(idat)[gc_cols]
          if(!bi_allelic){
            snp_form <- nchar(gc_colnames[1])/2
            gc_colnames_1 <- substr(gc_colnames, 1, snp_form)
            gc_colnames_2 <- substr(gc_colnames, snp_form + 1, snp_form*2)
          }
          for(k in 1:ncol(ps1_f)){
            # if not bi-allelic, need to compute the heterozygosity for this specific locus using the gc table handed along for this reason!
            if(!bi_allelic){
              het_cols_containing_k <- which(xor(gc_colnames_1 == colnames(ps1_f)[k], gc_colnames_2 == colnames(ps1_f)[k])) # hets (xor) with this allele
              tiho <- .fix..call(rowSums(idat[,..gc_cols][,..het_cols_containing_k])/intot)
              tjho <- .fix..call(rowSums(jdat[,..gc_cols][,..het_cols_containing_k])/jntot)
            }
            
            # otherwise we have the ho already
            else{
              tiho <- idat$ho
              tjho <- jdat$ho
            }
            
            
            # where this allele is absent in both populations, ignore it's contribution.
            absent <- which(ps1_f[[k]] + ps2_f[[k]] == 0)
            if(length(absent) != 0){
              tiho[absent] <- 0
              tjho[absent] <- 0
            }
            
            parts[[k]] <-.per_all_f_stat_components(intot, jntot, ps1_f[[k]], ps2_f[[k]], r, nbar, nc, tiho, tjho)
          }
          a <- matrix(unlist(purrr::map(parts, "a")), ncol = length(parts))
          b <- matrix(unlist(purrr::map(parts, "b")), ncol = length(parts))
          c <- matrix(unlist(purrr::map(parts, "c")), ncol = length(parts))
          
          data.table::set(out$a, j = c.col, value = rowSums(a)) # write a
          data.table::set(out$b, j = c.col, value = rowSums(b)) # write b
          data.table::set(out$c, j = c.col, value = rowSums(c)) # write v
          # Fst <- rowSums(a)/rowSums(a + b + c)
          
          data.table::set(out$Fst, j = c.col, value = rowSums(a)/rowSums(a + b + c)) # write fst
        }
        
        else{
          stop("Please select a method of calculating FST.\nOptions:\n\tWC: Weir and Cockerham (1984).\n\tgenepop: Genepop")
        }
        
        
        # update pnk
        pnk <- data.table::set(pnk, prog:(prog + nrow(idat) - 1), 1L, paste0(idat$subfacet, "~", jdat$subfacet))
        pnk <- data.table::set(pnk, prog:(prog + nrow(idat) - 1), 2L, as.integer(intot + jntot)*2) # back to allele count
        prog <- prog + as.integer(nrow(idat))
        
        
        c.col <- c.col + 1L #agument c.col
      }
    }
    
    # melt, cbind pnk
    suppressWarnings(out <- cbind(data.table::melt(out$Fst),
                                  data.table::melt(out$a)[,2],
                                  data.table::melt(out$b)[,2],
                                  data.table::melt(out$c)[,2]))
    colnames(out) <- c("comparison", "fst", "a", "b", "c")
    # meta.cols <- (1:ncol(x))[-c(ac_cols, ncol(x))]
    out$nk <- pnk$ntotal
    meta.cols <- 2:6
    out <- cbind(out[,"comparison"], .fix..call(x[x$subfacet == pops[1],..meta_colnames]), .fix..call(out[,..meta.cols]))
    
    # return
    return(out)
  }
  
  
  one_wm <- function(x){
    out <- x[,stats::weighted.mean(a, w = nk, na.rm = T)/
               stats::weighted.mean(a + b + c, w = nk, na.rm = T), 
             by = list(comparison)] # ratio of averages
    
    # out <- x[,weighted.mean(a/(a + b + c), w = nk, na.rm = T), 
    #          by = list(comparison)] # this is the average of ratios...
    colnames(out)[which(colnames(out) == "V1")] <- "weighted_mean_fst"
    return(out)
  }
  
  
  one_run <- function(x, method, facet, g.filename, rac = NULL){
    if(method == "genepop"){
      real_out <- func(x, method, facets = facet, g.filename = g.filename)
      real_wm <- real_out[[2]]
      real_out <- real_out[[1]]
    }
    else{
      if(bi_allelic){
        if(is.null(rac)){
          rac <- cbind(x@facet.meta, x@geno.tables$as, ho = x@stats$ho)
        }
        real_out <- func(rac[which(rac$facet == facet),], method, facet, ac_cols = which(colnames(rac) %in% colnames(x@geno.tables$as)), meta_colnames = colnames(snp.meta(x)))
        real_wm <- one_wm(real_out)
      }
      else{
        if(is.null(rac)){
          rac <- cbind(x@facet.meta, x@geno.tables$as, x@geno.tables$gs)
        }
        real_out <- func(rac[which(rac$facet == facet),], method, facet, 
                         ac_cols = which(colnames(rac) %in% colnames(x@geno.tables$as)), 
                         meta_colnames = colnames(snp.meta(x)),
                         gc_cols = which(colnames(rac) %in% colnames(x@geno.tables$gs)))
        real_wm <- one_wm(real_out)
      }
    }
    
    return(list(real_out = real_out, real_wm = real_wm))
  }
  
  #============================run with real data===============
  real <- vector("list", length = length(facets))
  names(real) <- facets
  if(method == "wc" & bi_allelic){
    needed <- .check_calced_stats(x, facets, "ho")
    needed <- facets[which(!unlist(needed))]
    if(length(needed) > 0){
      x <- calc_ho(x, needed)
    }
  }
  
  
  for(i in 1:length(facets)){
    if(method == "genepop"){
      .make_it_quiet(format_snps(x, "genepop", facets[i], outfile = paste0(facets[i], "_genepop_input.txt")))
    }
    
    real[[i]] <- one_run(x,
                         method =  method, 
                         facet = facets[i], 
                         g.filename = paste0(facets[i], "_genepop_input.txt"))
    
    if(method == "genepop"){
      real[[i]]$real_wm <- data.frame(comparison = names(real[[i]]$real_wm), weighted_mean_fst = real[[i]]$real_wm)
    }
  }
  
  
  #============run bootstraps if requested=======================
  if(!isFALSE(boot)){
    for(f in 1:length(facets)){
      cat("Conducting bootstraps, facet = ", facets[f], "\n")
      cat("Preparing input data...\n")
      
      if(method == "genepop"){
        gp.filenames <- .boot_genepop(paste0(facets[f], "_genepop_input.txt"), boot)
        wc.inputs <- NULL
      }
      
      
      cat("Done.\nCalculating Fst...\n")
      
      # do the bootstrapps
      if(isFALSE(boot_par)){
        boots <- vector("list", boot)
        for(i in 1:boot){
          cat("Boot", i, "out of", boot, "\n")
          if(method == "wc"){
            gp.filenames <- NULL
            wc.inputs <- .boot_as(x, 1, facets[f], ret_gs = !bi_allelic)
          }
          boots[[i]] <- one_run(x, method = method, facet = facets[f], gp.filenames[i], wc.inputs[[1]])$real_wm
        }
      }
      
      
      else{
        
        boot <- 1:boot
        if(boot_par < length(boot)){
          pboot <- split(boot, rep(1:boot_par, length.out = length(boot), each = ceiling(length(boot)/boot_par)))
        }
        else{
          boot_par <- length(boot)
          pboot <- split(boot, 1:boot_par, drop = F)
        }
        
        cl <- parallel::makePSOCKcluster(boot_par)
        doParallel::registerDoParallel(cl)
        
        
        
        #prepare reporting function
        ntasks <- length(pboot)
        # progress <- function(n) cat(sprintf("Job %d out of", n), ntasks, "is complete.\n")
        # opts <- list(progress=progress)
        
        #loop through each set
        boots <- foreach::foreach(q = 1:ntasks,
                                  .packages = c("snpR", "data.table"), .export = ".make_it_quiet"
        ) %dopar% {
          
          boots <- vector("list", length(pboot[[q]]))
          for(i in 1:length(pboot[[q]])){
            if(method == "wc"){
              gp.filenames <- NULL
              wc.inputs <- .boot_as(x, 1, facets[f], ret_gs = !bi_allelic)
            }
            boots[[i]] <- one_run(x, method = method, facet = facets[f], 
                                  gp.filenames[pboot[[q]][i]], wc.inputs[[1]])$real_wm
          }
          boots
        }
        
        
        parallel::stopCluster(cl)
        boots <- unlist(boots, recursive = F)
        boot <- max(boot)
      }
      
      
      # bind
      if(method == "genepop"){
        boots <- data.frame(boots)
        colnames(boots) <- NULL
      }
      else{
        boots <- purrr::map(boots, "weighted_mean_fst")
        names(boots) <- 1:length(boots)
        boots <- dplyr::bind_cols(boots)
      }
      
      # calculate p-values vs real
      real[[f]]$real_wm$weighted_mean_fst_p <- (rowSums(as.matrix(boots) >= real[[f]]$real_wm$weighted_mean_fst) + 1)/(boot + 1)
    }
    cat("Bootstraps complete.\n")
  }
  
  #==================merge and return============================
  cat("Collating results...")
  
  # merge the means
  real_wm <- data.table::rbindlist(purrr::map(real, "real_wm"), idcol = "facet")
  colnames(real_wm)[which(colnames(real_wm) == "comparison")] <- "subfacet"
  real_wm$snp.facet <- ".base"
  real_wm$snp.subfacet <- ".base"
  x <- .merge.snpR.stats(x, real_wm, "weighted.means")
  
  
  # merge the per-snp
  real_out <- data.table::rbindlist(purrr::map(real, "real_out"), idcol = "facet")
  ## zfst and fst/1-fst
  if(zfst){
    real_out[, zfst := (fst - mean(fst, na.rm = TRUE))/stats::sd(fst, na.rm = TRUE), by = .(comparison, facet)]
  }
  if(fst_over_one_minus_fst){
    real_out[, fst_id := fst/(1 - fst)]
  }
  col.ord <- c(which(colnames(real_out) != "nk"), which(colnames(real_out) == "nk"))
  real_out <- .fix..call(real_out[,..col.ord])
  
  ## save components
  component_cols <- which(colnames(real_out) %in% c("a", "b", "c"))
  if(!keep_components){
    real_out <- .fix..call(real_out[,-..component_cols])
  }
  else{
    # need to give a more unique name than "a", "b", "c" to avoid issues with get.snpR.stats()
    colnames(real_out)[match(c("a", "b", "c"), colnames(real_out))] <- paste0("var_comp_", c("a", "b", "c"))
  }
  x <- .merge.snpR.stats(x, real_out, type = "pairwise")
  
  # cleanup
  if(cleanup){
    if(exists("gp.filenames")){
      unlink(gp.filenames, recursive = TRUE)
    }
    if(method == "genepop"){
      file.remove(paste0(facets, "_genepop_input.txt"))
      file.remove("fichier.in", "cmdline.txt")
      genepop::clean_workdir()
    }
  }
  
  
  
  # return
  x <- .update_calced_stats(x, facets, "fst", "snp")
  if(method == "wc"){
    x <- .update_citations(x, "Weir1984", "fst", "Pairwise FST calculation")
  }
  else if(method == "genepop"){
    x <- .update_citations(x, "Rousset2008", "fst", "Pairwise FST calculation")
  }
  
  return(x)
}

#' Calculate FIS for individual populations.
#'
#' Calculates FIS for each individual level in each provided sample facet
#' according to Weir and Cockerham (1984).
#'
#' Note that FIS is calculated by considering \emph{only data from individual
#' sample levels}! This means that individual and sub-population variances are
#' only considered within each sub-population. If snp facets are provided,
#' weighted means will be provided for each snp facet level, although raw FIS
#' values are calculated on a per-snp basis and thus ignore these levels.
#'
#' If the base facet (facets = NULL or facets = ".base") is requested, FIS will
#' compare individual to total variance across all samples instead (equivalent
#' to overall FIT).
#'
#' @param x snpRdata. Input SNP data.
#' @param facets character. Categorical metadata variables by which to break up
#'   analysis. See \code{\link{Facets_in_snpR}} for more details.
#' @param keep_components logical, default FALSE. If TRUE, the variance
#'   components "b" and "c" will be held and accessible from the
#'   \code{$single} element (named "var_comp_b" and
#'   "var_comp_c", respectively) using the usual \code{\link{get.snpR.stats}}
#'   method. This may be useful if working with very large datasets that need to
#'   be run with separate objects for each chromosome, etc. for memory purposes.
#'   Weighted averages can be generated identically to those from snpR by taking
#'   the weighted mean (via the \code{\link[stats]{weighted.mean}}) of "c"
#'   divided by the sum of the weighted mean of "b" + "c" using the
#'   number of SNPs called in a comparison (returned in the "nk" column from
#'   \code{\link{get.snpR.stats}}) as weights within each population comparison.
#'   Note that this is different than taking the weighted mean of a/(a + b + c)!
#'
#' @export
#' @author William Hemstrom
#' @references Weir and Cockerham (1984). \emph{Evolution}
#'
#' @examples
#' x <- calc_fis(stickSNPs, c("pop", "pop.chr"))
#' get.snpR.stats(x, c("pop", "pop.chr"), "fis")
#' 
calc_fis <- function(x, facets = NULL, keep_components = FALSE){
  ..ac.cols <- ..meta.cols <- ..nk.cols <- ..nk.col <- ..gc_cols <- ..het_cols_containing_k <- subfacet <- snp.subfacet <- NULL
  var_comp_c <- nk <- var_comp_b <- facet <- ..rm_cols <- NULL
  #============================sanity and facet checks========================
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  # add any missing facets
  ofacets <- facets
  facets <- .check.snpR.facet.request(x, facets, return.type = T, fill_with_base = F)
  facets <- facets[[1]]
  if(!all(facets %in% x@facets)){
    .make_it_quiet(x <- .add.facets.snpR.data(x, facets))
  }
  
  #===========subfunctions=====================================================
  
  # cracked the problem: ho_i needs to be the proportion of individuals heterozygous for an allele, not heterozygous for all alleles! Need to fix for Fst too. Revert these to there bi-allelic form for that case since it is faster if ac is already calced?
  func <- function(ac_i, ho_i){
    count_tot <- rowSums(ac_i)
    parts <- vector("list", ncol(ac_i))
    
    for(i in 1:ncol(ac_i)){
      ps1 <- ac_i[[i]]
      tho <- ho_i
      absent <- which(ps1 == 0)
      if(length(absent) > 0){
        tho[absent] <- 0
      }
      ps1 <- ps1/count_tot
      
      parts[[i]] <- .per_all_f_stat_components(intot = count_tot/2,
                                               jntot = NA,
                                               ps1 = ps1,
                                               ps2 = NA,
                                               r = 1,
                                               nbar = count_tot/2, 
                                               nc = 0, 
                                               iho = tho, 
                                               jho = NA)
    }
    
    b <- matrix(unlist(purrr::map(parts, "b")), ncol = length(parts))
    c <- matrix(unlist(purrr::map(parts, "c")), ncol = length(parts))
    
    
    if(any(is.na(b))){
      b[is.na(b)] <- 0
    }
    if(any(is.na(c))){
      c[is.na(c)] <- 0
    }
    
    fis <- 1 - rowSums(c)/rowSums(b + c)
    
    return(data.frame(fis = fis, var_comp_c = rowSums(c), var_comp_b = rowSums(b)))
  }
  
  #===========run=============================================================
  # add ho
  needed <- .check_calced_stats(x, facets, "ho")
  needed <- facets[which(!unlist(needed))]
  if(length(needed) > 0){
    x <- calc_ho(x, needed)
  }
  
  # grab inputs
  samp.facets <- .check.snpR.facet.request(x, facets, fill_with_base = TRUE, return_base_when_empty = TRUE)
  meta.match <- which(x@facet.meta$facet %in% samp.facets)
  
  tas <- cbind(as.data.table(x@geno.tables$as[meta.match,]), as.data.table(x@facet.meta[meta.match,]))
  
  stat.match <- which(x@stats$facet %in% samp.facets)
  ho <- as.data.table(x@stats[stat.match,])
  
  input <- merge(tas, ho, by = c("facet", "subfacet", "facet.type", colnames(x@snp.meta)))
  
  # run and finish
  ac.cols <- which(colnames(input) %in% colnames(x@geno.tables$as))
  out <- func(.fix..call(input[,..ac.cols]), input$ho)
  
  meta.cols <- which(colnames(input) %in% colnames(x@facet.meta))
  out <- cbind(.fix..call(input[,..meta.cols]),out)
  out$nk <- rowSums(.fix..call(input[,..ac.cols]))
  
  
  # get the weighted mean ratio of averages and merge in
  wm <- vector("list", length(ofacets))
  for(i in 1:length(ofacets)){
    snp.facet <- .check.snpR.facet.request(x, ofacets[i], remove.type = "sample")
    if(snp.facet != ".base"){
      split.snp.facet <- unlist(.split.facet(snp.facet))
      out$snp.subfacet <- .paste.by.facet(out, split.snp.facet)
      
      wm <- out[,1 - (stats::weighted.mean(var_comp_c, w = nk, na.rm = T)/
                        stats::weighted.mean(var_comp_b + var_comp_c, w = nk, na.rm = T)), 
                by = list(facet, subfacet, snp.subfacet)]
      wm$snp.facet <- snp.facet
      out$snp.subfacet <- NULL
    }
    else{
      wm <- out[,1 - (stats::weighted.mean(var_comp_c, w = nk, na.rm = T)/
                        stats::weighted.mean(var_comp_b + var_comp_c, w = nk, na.rm = T)), 
                by = list(facet, subfacet)] # ratio of averages
      wm$snp.subfacet <- ".base"
      wm$snp.facet <- ".base"
    }
    
    colnames(wm)[which(colnames(wm) == "V1")] <- "weighted_mean_fis"
    x <- .merge.snpR.stats(x, wm, "weighted.means")
  }
  
  
  
  
  #========return==============
  if(!keep_components){
    rm_cols <- which(colnames(out) %in% c("var_comp_c", "var_comp_b", "nk"))
    out <- .fix..call(out[,-..rm_cols])
  }
  
  x <- .merge.snpR.stats(x, out)
  x <- .update_calced_stats(x, ofacets, "fis")
  x <- .update_citations(x, "Weir1984", "fis", "FIS calculation")
  
  return(x)
}

#'@export
#'@describeIn calc_single_stats observed heterozygosity
calc_ho <- function(x, facets = NULL){
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  # add any missing facets
  ofacets <- facets
  facets <- .check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    .make_it_quiet(x <- .add.facets.snpR.data(x, facets))
  }
  
  out <- .apply.snpR.facets(x,
                            facets = facets,
                            req = "gs",
                            fun = .ho_func,
                            snp_form = x@snp.form,
                            case = "ps")
  colnames(out)[ncol(out)] <- "ho"
  
  x <- .merge.snpR.stats(x, out)
  x <- .calc_weighted_stats(x, ofacets, type = "single", "ho")
  x <- .update_calced_stats(x, facets, "ho", "snp")
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
  facets <- .check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    .make_it_quiet(x <- .add.facets.snpR.data(x, facets))
  }
  
  out <- .apply.snpR.facets(x, facets, "meta.gs", func, case = "ps.pf")
  colnames(out)[ncol(out)] <- "pa"
  x <- .merge.snpR.stats(x, out)
  x <- .update_calced_stats(x, facets, "pa", "snp")
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
#'so is not preferred during casual or preliminary analysis. Either method will
#'calculate D', r-squared, and the p-value for that r-squared.
#'
#'Since this process involves many pairwise comparisons, it can be very slow.
#'
#'In contrast, Burrow's Composite Linkage Disequilibrium (CLD) can be calculated
#'very quickly via the \code{\link{cor}} function from base R.
#'\code{LD_full_pairwise} will perform this method alongside the other methods
#'if cld = TRUE and by itself if cld = "only". For most analyses, this will be
#'sufficient and much faster than the other methods. This is the default
#'behavior.
#'
#'The data can be broken up categorically by either SNP and/or sample metadata,
#'as described in \code{\link{Facets_in_snpR}}.
#'
#'Heatmaps of the resulting data can be easily plotted using
#'\code{\link{plot_pairwise_ld_heatmap}}
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
#'@param CLD TRUE, FALSE, or "only", default "only". Specifies if the CLD method
#'  should be used either in addition to or instead of default methods. See
#'  details.
#'@param use.ME logical, default FALSE. Specifies if the
#'  Minimization-Expectation haplotype estimation should be used. See details.
#'@param sigma numeric, default 0.0001. If the ME method is used, specifies the
#'  minimum difference required between steps before haplotype frequencies are
#'  accepted.
#'@param verbose Logical, default FALSE. If TRUE, some progress updates will be
#'  reported.
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
#' x <- calc_pairwise_ld(stickSNPs, facets = "chr.pop")
#' get.snpR.stats(x, "chr.pop", "LD")
#'
#' ## standard haplotype frequency estimation
#' x <- calc_pairwise_ld(stickSNPs, facets = "chr.pop", CLD = FALSE)
#' get.snpR.stats(x, "chr.pop", "LD")
#' }
#'
#' # subset for specific subfacets (ASP and OPL, chromosome IX)
#' x <- calc_pairwise_ld(stickSNPs, facets = "chr.pop",
#'                       subfacets = list(pop = c("ASP", "OPL"), 
#'                                        chr = "groupIX"))
#' get.snpR.stats(x, "chr.pop", "LD")
#' 
#' \dontrun{
#' ## not run, really slow
#' # ME haplotype estimation
#' x <- calc_pairwise_ld(stickSNPs, facets = "chr.pop", 
#'                       CLD = FALSE, use.ME = TRUE,
#'                       subfacets = list(pop = c("ASP", "OPL"), 
#'                                        chr = "groupIX"))
#' get.snpR.stats(x, "chr.pop", "LD")
#' }
calc_pairwise_ld <- function(x, facets = NULL, subfacets = NULL, ss = FALSE,
                             par = FALSE, CLD = "only", use.ME = FALSE, sigma = 0.0001,
                             verbose = FALSE){
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
  
  if(par > 1){
    cinst <- FALSE
    if(is.null(facets)){
      cinst <- TRUE
    }
    else if(any(facets == ".base")){
      cinst <- TRUE
    }
    
    if(cinst){
      .check.installed("bigmemory")
      .check.installed("bigtabulate")
    }
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
    m2 <- dplyr::summarise_all(m2, list(~sum(.)))
    
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
    
    
    
    # First, make a guess at the starting haplotype frequencies. We'll do this by taking the unambiguous haplotype frequencies,
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
    ehap.counts <- nhap.counts + .5*doub.het.sum # assuming that both haplopairs are equally likely in the double het
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
      n1hap.freqs <- haptable # grab the known haplotype frequencies form the unambiguous phenotypes again.
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
    mgcv$value <- as.factor(mgcv$value)
    cnames <- levels(mgcv$value)
    ghapmat <- bigtabulate::bigtabulate(mgcv, ccols = which(colnames(mgcv) %in% c("Var1", "value")))
    colnames(ghapmat) <- cnames
    
    #3) clean the table
    ##grab column names
    gl <- colnames(ghapmat)
    
    ##remove anything with missing data and double hets
    ## Keming: this line can be edited--remove the last which statement, save as a new variable which identifies
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
    het_l2 <- substr(gl, (sform*2) + 1, sform*3) != substr(gl, (sform*3+1), sform*4) #columns where the second locus is het
    het_l2 <- ghapmat[,het_l2]
    
    #fix Weird cases where one of these isn't a matrix because only one haplotype falls into the category.
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
  LD_func <- function(x, meta, mDat, snp.list, verbose = FALSE){
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
      if(verbose){
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
                       rsq = NA, Dprime = NA, pval = NA, row.names = NULL)
        
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
    sample.facets <- .check.snpR.facet.request(x, facets) # the sample level facets that we are working with.
    
    # initialize the output list
    out <- vector(mode = "list", length(sample.facets))
    names(out) <- sample.facets
    if(any(facet.types == "snp")){
      out <- c(out, list(.base = list(.base = list(snps = vector("list", nrow(x)), samps = 1:nrow(x@sample.meta)))))
    }
    
    # loop through each facet, do different things depending on facet type
    for(i in 1:length(facets)){
      
      # grab the facet level list we are writing to.
      t.samp.facet <- .check.snpR.facet.request(x, facets[i])
      write.facet <- which(names(out) == t.samp.facet)
      if(length(write.facet) == 0){
        t.samp.facet <- ".base"
        write.facet <- which(names(out) == ".base")
      }
      t.facet <- unlist(.split.facet(facets[i]))
      
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
        samp.facet <- .check.snpR.facet.request(x, facets[i])
        if(is.null(samp.facet)){samp.facet <- ".base"}
        snp.facet <- .check.snpR.facet.request(x, facets[i], remove.type = "sample")
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
  func <- function(x, facets, snp.facets, par, verbose){
    
    facet.types <- facets[[2]]
    facets <- facets[[1]]
    
    #=====================call functions=========
    # call these functions (possibly in parallel) according to supplied levels.
    
    if(verbose){cat("Beginning LD calculation...\n")}
    
    #=====================no facets==============
    if( (length(facets) == 1 & facets[1] == ".base") | all(facet.types == "snp")){
      if(length(facets) == 1 & facets[1] == ".base"){
        comps <- determine.comparison.snps(x, facets, "sample")
      }
      
      else{
        comps <- determine.comparison.snps(x, facets, facet.types)
      }
      
      
      if(verbose){cat("No facets specified.\n")}
      
      # grab metadata, mDat
      meta <- x@snp.meta
      mDat <- x@mDat
      
      
      # run in parallel if requested
      if(is.numeric(par)){
        if(verbose){cat("Running in parallel.\n\t")}
        
        # each thread needs to be given a roughly equal number of comparisons to do
        ncomps <-  length(unlist(comps[[1]][[1]]$snps)) # number of comparisons
        split <- (ncomps)/par #number of comparisons to do per core
        split <- ceiling(split)
        if(verbose){cat("At least", split, "pairwise comparisons per processor.\n")}
        
        # need to figure out which comps entries to null out for each processor.
        comps.per.snp <- unlist(lapply(comps[[1]][[1]]$snps, length))
        rproc <- ceiling(cumsum(as.numeric(comps.per.snp))/split) # which processor should each comparison be assigned to?
        
        #now need to start the parallel job:
        cl <- parallel::makePSOCKcluster(par)
        doParallel::registerDoParallel(cl)
        
        #prepare reporting function
        ntasks <- par
        # if(verbose){
        #   progress <- function(n) cat(sprintf("Part %d out of",n), ntasks, "is complete.\n")
        #   opts <- list(progress=progress)
        # }
        # else{
        #   opts <- list()
        # }
        
        # initialize and store things
        x_storage <- as.matrix(as.data.frame(x))
        na.test <- suppressWarnings(as.numeric(x_storage[1]))
        #save the info as a bigmatrix if it can be safely converted to numeric. Usually this is true for ms but not necessarily other data types.
        if(!is.na(na.test)){
          if(as.numeric(x_storage[1]) == x_storage[1]){
            if(verbose){cat("Saving matrix as big.matrix object for quicker sharing.\n")}
            xb <- bigmemory::as.big.matrix(x, type = "char")
            xbd <- bigmemory::describe(xb)
            remove(x)
          }
        }
        meta_storage <- x@snp.meta
        mDat_storage <- x@mDat
        t.comps <- comps[[1]][[1]]$snps
        
        if(verbose){cat("Begining run.\n")}
        
        # run the LD calculations
        output <- foreach::foreach(q = 1:ntasks, .packages = c("bigmemory", "dplyr"), .inorder = TRUE,
                                   .export = c("LD_func", "tabulate_haplotypes", "GtoH")) %dopar% {
                                     if(exists("xbd")){
                                       x_storage <- bigmemory::attach.big.matrix(xbd)
                                     }
                                     
                                     # get comps and run
                                     t.comps[which(rproc != q)] <- vector("list", sum(rproc != q)) # null out any comparisons we aren't doing
                                     t.comps <- list(snps = t.comps, samps = comps[[1]][[1]]$samps)
                                     
                                     LD_func(x = x_storage, snp.list = t.comps,
                                             meta = meta_storage, mDat = mDat_storage,
                                             verbose = verbose)
                                   }
        
        #release cores
        parallel::stopCluster(cl)
        
        
        if(verbose){cat("LD computation completed. Preparing results.\n\t")}
        
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
        out <- LD_func(x, meta, snp.list = comps[[1]][[1]], mDat = mDat, verbose)
        
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
    # Each facet is a level to break down by. "pop" means break by pop, c("pop", "chr") means break twice, once by pop, once by chr,
    # c("pop.chr") means to break by pop and chr.
    # For each sample level facet, we must loop through all snps, since D values will be different depending on what samples we look at.
    # These must be looped through separately!
    # If there are multiple snp level facets requested, there is no reason to do re-do snp/snp comparisons within each sample level facet. Just do the all of the relevant snp/snp comparisons.
    # If there are complex facets with repeated sample levels (c("pop.chr", "pop.subchr")), then same deal.
    
    # So:
    #     For each sample level facet:
    #       check if we've run the facet before
    #       For each level of those facets:
    #         For each snp:
    #           Figure out which snps we need to compare to across all snp level facets.
    #           Remove any comparisons that we've already done!
    #           Pass the genotypes and per snp comparison info to the LD_func (need to edit that function slightly to accommodate)
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
          if(verbose){cat("Subfacet #:", progress, "of", tot_subfacets, " Name:", paste0(names(comps)[i], " " , names(comps[[i]])[j]), "\n")}
          out <- LD_func(x, meta = x@snp.meta, mDat = x@mDat, snp.list = comps[[i]][[j]], verbose = verbose)
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
      cl <- parallel::makePSOCKcluster(par)
      doParallel::registerDoParallel(cl)
      
      #prepare reporting function
      ntasks <- tot_subfacets
      # if(verbose){
      #   progress <- function(n) cat(sprintf("Facet %d out of", n), ntasks, "is complete.\n")
      #   opts <- list(progress=progress)
      # }
      # else{
      #   opts <- list()
      # }
      
      x_storage <- as.matrix(as.data.frame(x))
      meta_storage <- x@snp.meta
      mDat_storage <- x@mDat
      
      
      #loop through each set of facets
      output <- foreach::foreach(i = 1:ntasks, .packages = c("dplyr", "reshape2", "matrixStats", "bigtabulate", "snpR"), .inorder = TRUE,
                                 .export = c("LD_func", "tabulate_haplotypes", "GtoH")) %dopar% {
                                   t.task <- task_list[i,]
                                   t.facet <- t.task[1]
                                   t.subfacet <- t.task[2]
                                   LD_func(x_storage, meta = meta_storage, mDat = mDat_storage, snp.list = comps[[t.facet]][[t.subfacet]], verbose = verbose)
                                 }
      
      #release cores and clean up
      parallel::stopCluster(cl)
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
  
  #========================prepare and pass the primary function to .apply.snpR.facets==================
  # add any missing facets
  x <- .add.facets.snpR.data(x, facets)
  
  #subset data if requested:
  if(!(is.null(subfacets[1]))){
    old.x <- x
    ssfacets <- names(subfacets)
    
    # check complex facets
    complex.sfacets <- .check.snpR.facet.request(x, ssfacets, remove.type = "simple", fill_with_base = FALSE, return_base_when_empty = FALSE)
    if(length(complex.sfacets) > 0){
      stop(paste0("Complex (snp and sample) level SUBFACETS not accepted. Providing these as separate snp and sample subfacets will run only snps/samples contained in both levels. Bad facets: ",
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
    ssfacet.types <- .check.snpR.facet.request(x, ssfacets, "none", T)[[2]]
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
        suppressWarnings(.make_it_quiet(x <- .subset_snpR_data(x,
                                                               facets = sample.facets,
                                                               subfacets = sample.subfacets,
                                                               snp.facets = snp.facets,
                                                               snp.subfacets = snp.subfacets)))
      }
      else{
        suppressWarnings(.make_it_quiet(x <- .subset_snpR_data(x,
                                                               facets = sample.facets,
                                                               subfacets = sample.subfacets)))
      }
    }
    else if(filter_snp_facets){
      suppressWarnings(.make_it_quiet(x <- .subset_snpR_data(x,
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
    facets_trad <- .check.snpR.facet.request(x, facets, remove.type = "none", return.type = T)
    
    # run the function
    out <- func(x, facets = facets_trad, snp.facets = snp.facets, par = par, verbose = verbose)
    
    # add to snpRdata object and return
    if(exists("old.x")){
      out <- .merge.snpR.stats(old.x, out, "LD")
      
    }
    else{
      out <- .merge.snpR.stats(x, out, "LD")
    }
  }
  # run CLD components
  if(CLD != F){
    # run the function
    out <- .calc_CLD(x, facets, par, verbose = verbose)
    
    # add to snpRdata object and return
    if(CLD != "only"){
      out <- .merge.snpR.stats(x, out, "LD")
    }
    else{
      if(exists("old.x")){
        out <- .merge.snpR.stats(old.x, out, "LD")
        
      }
      else{
        out <- .merge.snpR.stats(x, out, "LD")
      }
    }
  }
  
  #======================update and return============
  out <- .update_calced_stats(out, facets, "LD")
  if(!isFALSE(CLD)){
    out <- .update_citations(out, "Cockerham1977", "LD_CLD", "Burrows' Composite Linkage Disequlibrium")
  }
  if(CLD != "only"){
    out <- .update_citations(out, "Lewontin1964", "LD_D", "D and D'")
    if(use.ME){
      out <- .update_citations(out, "Excoffier1995", "LD_MLE", "Maximization-Expectation Algorithim for caluculating haplotype frequencies")
    }
  }
  
  return(out)
}

# Calculate Burrow's composite LD. Internal.
#
# Called in \code{\link{calc_pairwise_ld}}. Complex function purely because the
# output needs to be in the same format as that function's and to support
# faceting without excessive recalculation.
#
# @param x snpRdata object
# @param facets facets to run
# @param par number of parallel cores
#
# @author William Hemstrom
.calc_CLD <- function(x, facets = NULL, par = FALSE, verbose = FALSE){
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
      prox <- prox[-which(is.na(prox$CLD)),]
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
    facet.names <- gsub("\\.$", "", facet.names)
    if(any(facet.names == "")){
      facet.names[which(facet.names == "")] <- ".base"
    }
    for(i in 1:length(facet.names)){
      if(facet.names[1] != ""){
        if(i == 23){
        }
        facet.names[i] <- .check.snpR.facet.request(x, facet.names[i], remove.type = "none")
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
  if(verbose){cat("Preparing data...\n")}
  
  suppressMessages(x <- .add.facets.snpR.data(x, facets))
  tasks <- .get.task.list(x, facets)
  x@sn$sn <- format_snps(x, "sn", interpolate = F)
  x@sn$type <- "FALSE"
  
  if(verbose){cat("Beginning LD calculation...\n")}
  # run the loop
  if(par == F){
    
    # initialize storage
    matrix_storage <- vector("list", nrow(tasks))
    prox_storage <- vector("list", nrow(tasks))
    
    
    #loop through each set of facets
    for(i in 1:nrow(tasks)){
      # run
      if(verbose){cat("Task #:", i, "of", nrow(tasks),
                      " Sample Facet:", paste0(tasks[i,1:2], collapse = "\t"),
                      " SNP Facet:", paste0(tasks[i,3:4], collapse = "\t"), "\n")}
      
      # can't integrate this part into do_CLD without screwing up the parallel due to s4 issues.
      suppressWarnings(y <- .subset_snpR_data(x, facets = tasks[i,1],
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
    if(verbose){cat("Running in parallel.\nSplitting up data...\n")}
    
    geno.storage <- vector("list", nrow(tasks))
    
    
    # can't do this part in parallel due to s4 issues.
    for(i in 1:nrow(tasks)){
      if(verbose){cat("Task #:", i, "of", nrow(tasks),
                      " Sample Facet:", paste0(tasks[i,1:2], collapse = "\t"),
                      " SNP Facet:", paste0(tasks[i,3:4], collapse = "\t"), "\n")}
      utils::capture.output(invisible(suppressWarnings(y <- .subset_snpR_data(x, facets = tasks[i,1],
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
    
    cl <- parallel::makePSOCKcluster(par)
    doParallel::registerDoParallel(cl)
    
    
    
    #prepare reporting function
    ntasks <- length(ptasks)
    # if(verbose){
    #   progress <- function(n) cat(sprintf("Job %d out of", n), ntasks, "is complete.\n")
    #   opts <- list(progress=progress)
    # }
    # else{
    #   opts <- list(opts)
    # }
    
    #loop through each set of facets
    output <- foreach::foreach(q = 1:ntasks,
                               .packages = c("dplyr", "reshape2", "matrixStats", "bigtabulate", "snpR", "data.table"),
                               .inorder = TRUE
    ) %dopar% {
      
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
    parallel::stopCluster(cl)
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
#'@describeIn calc_single_stats p-values for Hardy-Weinberg Equilibrium divergence
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
      nans <- is.nan(chi)
      chi[which(nans)] <- 0
      
      # calculate p-values
      out <- stats::pchisq(chi, ifelse(nans, 0, 1), lower.tail = FALSE)
    }
    
    # otherwise we have to use the looped version:
    else if(method == "exact"){
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
  facets <- .check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    .make_it_quiet(x <- .add.facets.snpR.data(x, facets))
  }
  
  out <- .apply.snpR.facets(x, facets, "gs", func, case = "ps", method = method)
  colnames(out)[ncol(out)] <- "pHWE"
  
  # apply corrections
  if(length(facets) == 1 & facets[1] == ".base"){
    fwe_case <- "overall"
  }
  if(fwe_method[1] != "none"){
    out <- .fwe_correction(out, levs = c("facet", "subfacet"), pcol = "pHWE", methods = fwe_method, case = fwe_case)
  }
  
  x <- .merge.snpR.stats(x, out)
  x <- .update_calced_stats(x, facets, "hwe", "snp")
  if(method == "exact"){
    x <- .update_citations(x, "Wigginton2005", "HWE", "Hardy-Weinburg Equilibrium, p-values via exact test")
  }
  
  return(x)
}

#'Calculate basic SNP statistics
#'
#'Automatically calculate most basic statistics from snpRdata. Calculates maf,
#'pi, ho, he, pairwise Fst, HWE divergence, finds private alleles, and uses Gaussian
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
#'  Options: \itemize{ \item{WC: } Weir and Cockerham (1984). \item{Weir: } Weir
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
#'x <- calc_basic_snp_stats(stickSNPs, "pop")
#'get.snpR.stats(x, "pop", stats = "single") # view basic stats
#'get.snpR.stats(x, "pop", stats = "fst") # view fst
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
    facet.types <- .check.snpR.facet.request(x, facets, "none", T)
    snp.facets <- which(facet.types[[2]] == "snp")
  }
  
  if(!is.null(sigma)){
    if(is.null(facets[1])){
      .sanity_check_window(x, sigma, step, nk = nk, stats.type = "single", facets = facets)
    }
    else if(any(facet.types[[2]] == "snp")){
      .sanity_check_window(x, sigma, step, nk = nk, stats.type = "single", facets = facets)
    }
    else{
      .sanity_check_window(x, sigma, step, nk = nk, stats.type =  c("pairwise", "single"), facets = facets)
    }
  }
  
  #=========stats===============
  # basic stats
  x <- calc_maf(x, facets)
  x <- calc_pi(x, facets)
  x <- calc_hwe(x, facets)
  x <- calc_ho(x, facets)
  x <- calc_he(x, facets)
  if(!is.null(facets[1]) & !any(.check.snpR.facet.request(x, facets, return.type = TRUE)[[2]] == ".base")){
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




# @param temporal_methods character, default c("Pollak", "Nei", "Jorde").
#   Methods to use for the temporal methods. See defaults for options. Not
#   currently supported.
# @param temporal_gens data.frame, default NULL. Work in progress, not
#   currently supported.


#' Calculate effective population size.
#'
#' Calculates effective population size for any given sample-level facets via
#' interface with the NeEstimator v2 program by Do et al. (2013).
#'
#' Since physical linkage can cause miss-estimation of Ne, an optional snp-level
#' facet can be provided which designates chromosomes or linkage groups. Only
#' pairwise LD values between SNPs on different facet levels will be used.
#'
#' Ne can be calculated via three different methods: \itemize{ \item{"LD"}
#' Linkage Disequilibrium based estimation. \item{"Het"} Heterozygote excess.
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
#'   bias Ne estimates, providing a value here is recommended.
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
#'   \item{"het"} Heterozygote excess. \item{"coan"} Coancestry based.}
#' @param max_ind_per_pop numeric, default NULL. Maximum number of individuals
#'   to consider per population.
#' @param outfile character, default "ne_out". Prefix for output files. Note
#'   that this function will return outputs, so there isn't a strong reason to
#'   check this. At the moment, this cannot be a full file path, just a file
#'   prefix ('test_ne' is OK, '~/tests/test_ne' is not).
#' @param verbose Logical, default FALSE. If TRUE, some progress updates will be
#'   reported. 
#' @param cleanup logical, default TRUE. If TRUE, the NeEstimator output
#'   directory will be removed following processing.
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
#' # calculate Ne, noting not to use LD between SNPs on the 
#' # same chromosome equivalent ("chr") for every population.
#' ne <- calc_ne(stickSNPs, facets = "pop", chr = "chr")
#' get.snpR.stats(ne, "pop", stat = "ne")}
#' 
#' @export
calc_ne <- function(x, facets = NULL, chr = NULL,
                    NeEstimator_path = "/usr/bin/Ne2-1.exe",
                    mating = "random",
                    pcrit = c(0.05, 0.02, 0.01),
                    methods = "LD",
                    # temporal_methods = c("Pollak", "Nei", "Jorde"),
                    # temporal_gens = NULL, 
                    max_ind_per_pop = NULL,
                    outfile = "ne_out", verbose = TRUE, cleanup = TRUE){
  #==============sanity checks and prep=================
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  msg <- character()
  NeEstimator_path <- normalizePath(NeEstimator_path)
  if(!file.exists(NeEstimator_path)){
    msg <- c(msg, "NeEstimator executable not found at provided path.\n")
  }
  else{
    suppressWarnings(test <- try(system(paste0(NeEstimator_path, " -check"), intern = TRUE), silent = TRUE))
    if(methods::is(test, "try-error")){
      msg <- c(msg, "NeEstimator executable not found at provided path.\n")
    }
    else if(test[1] != "Illegal argument!"){
      msg <- c(msg, "NeEstimator executable not found at provided path.\n")
    }
  }
  
  if("temporal" %in% tolower(methods)){msg <- c(msg, "Implementation for NeEstimator's temporal method is in progress but currently not supported.\n")}
  
  if(basename(outfile) != outfile){
    msg <- c(msg, "At the moment, the outfile name must be a file name prefix (such as 'test_ne'), not a full path (such as '~/tests/test_ne').\n")
  }
  
  facets <- .check.snpR.facet.request(x, facets, remove.type = "snp")
  if(length(msg) > 0){
    stop(msg)
  }
  
  
  wrn <- character(0)
  if(nsnps(x) > 5000){
    wrn <- c(wrn, "Large numbers of loci don't usually improve Ne estimates by a substantial margin and can massively increase run times and memory requirements. Loci counts of ~3,000 are usually sufficient to generate accurate Ne esitmates, but estimates are improved substantially by selecting loci with little missing data. See https://doi.org/10.1002/ece3.6016\n")
  }
  if(is.null(chr)){
    wrn <- c(wrn, "NeEstimator expects no physical linkage between loci. You may want to consider setting the 'chr' argument to only consider LD values between loci on different chromosomes/scaffolds/etc.\n")
  }
  
  if(length(wrn) > 0){
    message(wrn)
    warning(wrn)
  }
  
  
  
  #=============run====================================
  out <- vector("list", length(facets))
  for(i in 1:length(facets)){
    # write inputs
    write_neestimator_inputs(x = x,
                             facets = facets[i],
                             chr = chr,
                             methods = methods,
                             # temporal_methods = temporal_methods,
                             # temporal_gens = temporal_gens,
                             pcrit = pcrit,
                             mating = mating,
                             outfile = paste0(outfile, ".txt"),
                             max_ind_per_pop = max_ind_per_pop)
    
    # run
    run_neestimator(NeEstimator_path = NeEstimator_path,
                    data_path = "NeEstimator/", verbose = verbose)
    
    # parse
    out[[i]] <- parse_neestimator(path = "NeEstimator/",
                                  pattern = outfile,
                                  facets = facets[i],
                                  snpRdat = x)
  }
  
  names(out) <- facets
  out <- dplyr::bind_rows(out, .id = "facet")
  
  #=============return===========
  x <- .merge.snpR.stats(x, out, "pop")
  x <- .update_calced_stats(x, facets, "ne")
  x <- .update_citations(x, "Do2014", "ne", "Ne, via interface to NeEstimator")
  
  if(cleanup){
    if(Sys.info()["sysname"] == "Windows"){
      shell("rm -r ./NeEstimator")
    }
    else{
      system("rm -r ./NeEstimator")
    }
  }
  
  return(x)
}

#' Calculate genetic distances between individuals or groups of samples.
#'
#' Calculates the genetic distances between either individuals or the levels of
#' any sample specific facets, broken apart by any requested snp level facets.
#' For details on methods, see details.
#'
#'
#' If a sample facet is requested, distances are calculated via code adapted
#' from code derived from \code{\link[adegenet]{adegenet}}. Please cite them
#' alongside the tree-building and distance methods Available methods:
#' \itemize{\item{Edwards: } Angular distance as described in Edwards 1971.
#' \item{Nei: } Nei's (1978) genetic distance measure.}
#'
#' Otherwise, genetic distances are calculated via the \code{\link[stats]{dist}}
#' function using genotypes formatted numerically (the "sn" option in
#' \code{\link{format_snps}}). In all cases, missing genotypes are first
#' interpolated according to the method provided to the 'interpolate' argument.
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
#'
#' @references Edwards, A. W. F. (1971). Distances between populations on the
#'   basis of gene frequencies. Biometrics, 873-881.
#' @references Nei, M. (1978). Estimation of average heterozygosity and genetic
#'   distance from a small number of individuals. Genetics, 23, 341-369.
#' @references Jombart, T. (2008). adegenet: a R package for the multivariate
#'   analysis of genetic markers Bioinformatics 24: 1403-1405.
#'
#' @author William Hemstrom
#' @export
#' 
#' @examples
#' # by pop:
#' y <- calc_genetic_distances(stickSNPs, facets = "pop", method = "Edwards")
#' get.snpR.stats(y, "pop", "genetic_distance")
#'
#' # by chr and pop jointly
#' y <- calc_genetic_distances(stickSNPs, facets = "pop.chr", 
#'                             method = "Edwards")
#' get.snpR.stats(y, "pop.chr", "genetic_distance")
#'
#' \dontrun{
#' # by pop and fam separately
#' y <- calc_genetic_distances(stickSNPs, facets = c("pop", "fam"), 
#'                             method = "Edwards")
#' get.snpR.stats(y, c("pop", "chr"), "genetic_distance")
#'
#' # individuals across all snps + plot
#' y <- calc_genetic_distances(stickSNPs)
#' dat <- get.snpR.stats(y, stats = "genetic_distance")$.base$.base$Edwards
#' heatmap(as.matrix(dat))
#' }
calc_genetic_distances <- function(x, facets = NULL, method = "Edwards", interpolate = "bernoulli"){
  #============sanity checks=========
  msg <- character()
  
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  good.methods <- c("Edwards", "Nei")
  if(!method %in% good.methods){
    msg <- c(msg, paste0("Provided method not supported. Supported methods: ", paste0(good.methods, collapse = " ")))
  }
  
  # get the allele frequencies if not already provided
  facets <- .check.snpR.facet.request(x, facets, "none", T)
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
    needed.facets <- .check_calced_stats(x, facets, "allele_frequency_matrix")
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
  
  #=============run for facets with sample aggregation===============
  if(sample_facets_detected){
    out <- vector("list", length(x))
    names(out) <- names(x)
    # enter lapply hell--first level unlists once, second level runs the function if the values are matrices, unlists if not, third level can always run the function
    # the nice thing with this is that it should keep all of the names from x natively!
    out <- lapply(x, function(y){
      lapply(y, function(z) {
        if("matrix" %in% class(z)){
          .get_dist(z, method = method)
        }
        else{
          lapply(z, .get_dist, method = method)
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
        tasks <- .get.task.list(y, snp_facets[i])
        out_snp[[i]] <- vector("list", nrow(tasks))
        names(out_snp[[i]]) <- tasks[,4]
        snp_columns <- unlist(.split.facet(tasks[,3][1]))
        
        # run the tasks
        for(j in 1:nrow(tasks)){
          t_snp_cols <- y@snp.meta[,snp_columns, drop = F]
          t_snp_indices <- .paste.by.facet(t_snp_cols, snp_columns, sep = "    ")
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
  
  y <- .update_calced_stats(y, all_facets, paste0("genetic_distance_", method, "_", interpolate))
  if(method == "Edwards"){
    y <- .update_citations(y, "Edwards1971", "genetic_distance", "Edwards' Angular Genetic Distance")
  }
  else if(method == "Nei"){
    y <- .update_citations(y, "Nei1978", "genetic_distance", "Nei's genetic distance")
  }
  return(.merge.snpR.stats(y, out, "genetic_distances"))
  
}

#' Calculate Isolation by Distance
#'
#' Calculates Isolation by Distance (IBD) for snpRdata objects by comparing the
#' genetic distance between samples or sets of samples to the geographic
#' distances between samples or sets of samples. IBD calculated via a mantel
#' test.
#'
#' Genetic distance is calculated via \code{\link{calc_genetic_distances}}.
#' Geographic distances are taken as-is for individual-individual comparisons
#' and by finding the geographic mean of a group of samples when sample level
#' facets are requested using the methods described by
#' \code{\link[geosphere]{geomean}}. Note that this means that if many samples
#' were collected from the same location, and those samples compose a single
#' level of a facet, the mean sampling location will be that single location.
#'
#' IBD is calculated by comparing geographic and genetic distances using a
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
#' @param interpolate Character or FALSE, default "bernoulli". If transforming to
#'  "sn" format, notes the interpolation method to be used to fill missing data.
#'  Options are "bernoulli", "af", "iPCA", or FALSE. See 
#'  \code{\link{format_snps}} for details.
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
#' \dontrun{
#' y <- stickSNPs
#' sample.meta(y) <- cbind(sample.meta(y), x = rnorm(ncol(y)), y = rnorm(ncol(y)))
#' y <-calc_isolation_by_distance(y, facets = c(".base", "pop", "pop.chr","pop.chr.fam"))
#' res <- get.snpR.stats(y, "pop.chr", "ibd") # fetch results
#' res
#' plot(res$chr.pop$groupV$Edwards) # plot perms vs observed
#' }
calc_isolation_by_distance <- function(x, facets = NULL, x_y = c("x", "y"), genetic_distance_method = "Edwards", interpolate = "bernoulli", ...){
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
  
  pkg.check <- .check.installed("geosphere")
  if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
  
  pkg.check <- .check.installed("ade4")
  if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
  
  
  if(length(msg) > 0){
    stop(paste0(msg, collapse = "\n"))
  }
  
  #================prep and get genetic dist matrices========================
  facets <- .check.snpR.facet.request(x, facets, "none")
  x <- .add.facets.snpR.data(x, facets)
  
  # calc gds if needed and fetch
  cs <- .check_calced_stats(x, facets, paste0("genetic_distance_", genetic_distance_method, "_", interpolate))
  missing <- which(!unlist(cs))
  if(length(missing) > 0){
    x <- calc_genetic_distances(x, names(cs)[missing], genetic_distance_method, interpolate)
  }
  gd <- .get.snpR.stats(x, facets, "genetic_distance")
  
  #===============fetch geo dist matrices=======================
  # get the geo dist matrices
  snp.levs <- .check.snpR.facet.request(x, facets)
  geo_mats <- vector("list", length(snp.levs))
  names(geo_mats) <- snp.levs
  
  # get the geo matrices for each comparison. note that the same geo matrix is re-used across snp level facets.
  for(i in 1:length(geo_mats)){
    # get the categories for this
    categories <- .get.task.list(x, names(geo_mats)[i])[,1:2, drop = F]
    
    # if the categories are all base, easy peasy
    if(all(categories[,1] == ".base")){
      geo_mats$.base[[genetic_distance_method]] <- stats::dist(x@sample.meta[,x_y])
    }
    else{
      # otherwise need to break it down by sample level facet
      # find geographic mean coordiantes for each facet level
      gmeans <- t(apply(categories, 1, function(row){
        matches <- .fetch.sample.meta.matching.task.list(x, row)
        matches <- x@sample.meta[matches,, drop = FALSE]
        matches <- matches[,x_y, drop = FALSE]
        if(nrow(matches) == 1){
          return(matches)
        }
        colnames(matches) <- c("x", "y")
        return(geosphere::geomean(matches))
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
        
        tsampfacet <- .check.snpR.facet.request(x, names(gd)[i])
        
        # quick check for NAs
        if(any(is.na(gd[[i]][[j]][[k]]))){
          # find bad entries
          matgd <- as.matrix(gd[[i]][[j]][[k]])
          bad.entries <- which(rowSums(is.na(matgd)) > 0)
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
        
        if(attr(gd[[i]][[j]][[k]], "Size") != 0){
          ibd[[i]][[j]][[k]] <- ade4::mantel.randtest(gd[[i]][[j]][[k]], tgeodist, ...)
        }
        else{
          warning("All genetic distance measures NA for facet", paste(names(gd)[i], names(gd[[i]])[j], names(gd[[i]][[j]])[k]),"\n\t")
          ibd[[i]][[j]][[k]] <- NULL
        }
        
      }
    }
  }
  
  warning("Mantel tests for IBD may have problems with autocorrelation, and so should be interpreted carefully! See Patrick Meirmans (2012). The Trouble with Isolation by Distance, Molecular Ecology.\n")
  
  #===============return================================
  x <- .merge.snpR.stats(x, geo_mats, "geo_dists")
  x <- .merge.snpR.stats(x, ibd, "ibd")
  x <- .update_calced_stats(x, facets, c("geo_dist", "ibd"))
  x <- .update_citations(x, "Mantel209", "isolation_by_distance", "Mantel test used to generate p-value for genetic vs geographic distance.")
  
  return(x)
}



#'@export
#'@describeIn calc_single_stats expected heterozygosity
calc_he <- function(x, facets = NULL){
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  
  # add any missing facets
  ofacets <- facets
  facets <- .check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    .make_it_quiet(x <- .add.facets.snpR.data(x, facets))
  }
  
  # bi_allelic case, just need maf--very straightforward
  if(.is.bi_allelic(x)){
    he_func_bi <- function(as = NULL, maf){
      return(2 * maf$maf * (1 - maf$maf))
    }
    # add missing maf
    has_maf <- .check_calced_stats(x, facets, "maf")
    if(any(!unlist(has_maf))){
      x <- calc_maf(x,  .check.snpR.facet.request(x, facets)[which(!unlist(has_maf))])
    }
    
    # calculate he
    out <- .get.snpR.stats(x, facets = facets, type = "single")
    out$he <- he_func_bi(maf = out)
    ft <- .check.snpR.facet.request(x, facets, return.type = TRUE)
    out$facet.type <- ft[[2]][match(out$facet, ft[[1]])]
    out <- out[,c("facet", "subfacet", "facet.type", colnames(snp.meta(x)), "he")]
  }
  # non bi_allelic, need as but still easy
  else{
    he_func_nb <- function(as, maf = NULL){
      as <- as$as
      as <- as/rowSums(as)
      as <- as^2
      he <- 1 - rowSums(as)
      return(he)
    }
    
    out <- .apply.snpR.facets(x, facets, "gs", fun = he_func_nb, case = "ps")
    colnames(out)[ncol(out)] <- "he"
  }
  
  x <- .merge.snpR.stats(x, out)
  
  x <- .calc_weighted_stats(x, ofacets, type = "single", "he")
  x <- .update_calced_stats(x, facets, "he", "snp")
  return(x)
  
}

#' Calculate individual based heterozygosity.
#'
#' Calculates heterozygosity within individuals across all SNPs by either a
#' simple ratio of heterozygous to homozygous sites or standardized for
#' differences in genotyping success according to Coltman et al (1999).
#'
#' @section Het:Hom ratio:
#'
#'   Individual heterozygosity calculated as the number of heterozygous sites
#'   divided by the number of homozygous sites.
#'
#' @section \ifelse{html}{\out{H<sub>S</sub>}}{\eqn{H_S}}:
#'
#'   Calculates \ifelse{html}{\out{H<sub>S</sub>}}{\eqn{H_S}}, the mean
#'   heterozygosity of an individual standardized for unequal genotyping across
#'   individuals by dividing by the mean heterozygosity across all individuals
#'   *for the loci sequenced in that individual* (Coltman et al 1999). As a
#'   result, the global mean \ifelse{html}{\out{H<sub>S</sub>}}{\eqn{H_S}}
#'   should be roughly equal to 1, and that in `snpR` specifically the
#'   denominator is calculated across *all individuals in all facet levels* if
#'   facets are specified instead of within populations. As a result, the
#'   weighted mean \ifelse{html}{\out{H<sub>S</sub>}}{\eqn{H_S}} in a specific
#'   population can be substantially different from one if a population is much
#'   less heterozygous.
#'
#' @param x snpRdata object
#' @param facets facets over which to split snps within samples. Takes only SNP
#'   level facets. See \code{\link{Facets_in_snpR}} for details.
#' @param complex_averages logical, default FALSE. If TRUE, will compute
#'   weighted averages for complex (snp + sample metadata) facets. This can be
#'   quite time consuming, and so is generally not recommended unless needed.
#'
#' @return A snpRdata object with heterozygote/homozygote ratios merged into the
#'   sample.stats slot.
#'   
#' @author William Hemstrom
#' 
#' @aliases calc_het_hom_ratio calc_hs
#' @name individual_heterozygosity
#' 
#' @references 
#' Coltman, D. W., Pilkington, J. G., Smith, J. A., & Pemberton, J.
#' M. (1999). Parasite-mediated selection against inbred Soay sheep in a
#' free-living, island population. Evolution, 53(4), 1259–1267. doi:
#' 10.1111/j.1558-5646.1999.tb04538.x
#' 
#' @examples
#' # base facet
#' x <- calc_het_hom_ratio(stickSNPs)
#' get.snpR.stats(x, stats = "het_hom_ratio")
#'
#' # facet by chromosome
#' x <- calc_het_hom_ratio(stickSNPs, "chr")
#' get.snpR.stats(x, "chr", stats = "het_hom_ratio")
#' 
#' # Getting population means:
#' x <- calc_hs(stickSNPs, "pop")
#' get.snpR.stats(x, "pop", stats = "hs")
#'
NULL

#'@export
#'@describeIn individual_heterozygosity Ratio of heterozygous to homozygous sites.
calc_het_hom_ratio <- function(x, facets = NULL, complex_averages = FALSE){
  
  #============run for each facet================
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  # add any missing facets
  ofacets <- facets
  if(!complex_averages){
    ofacets <- .check.snpR.facet.request(x, facets, "none")
  }
  facets <- .check.snpR.facet.request(x, facets, remove.type = "sample")
  
  out <- .apply.snpR.facets(x, facets, "genotypes", .heterozygosity, case = "psamp", mDat = x@mDat, method = "ratio")
  colnames(out)[which(colnames(out) == "stat")] <- "Het/Hom"
  x <- .merge.snpR.stats(x, out, "sample.stats")
  x <- .calc_weighted_stats(x, ofacets, "sample", "Het/Hom")
  
  x <- .update_calced_stats(x, facets, "ho_he_ratio")
  
  return(x)
}

#'@export
#'@describeIn individual_heterozygosity \ifelse{html}{\out{H<sub>S</sub>}}{\eqn{H_S}}, Individual Heterozygosity
calc_hs <- function(x, facets = NULL, complex_averages = FALSE){
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  
  # add any missing facets
  ofacets <- facets
  if(!complex_averages){
    ofacets <- .check.snpR.facet.request(x, facets, "none")
  }
  facets <- .check.snpR.facet.request(x, facets, remove.type = "sample")
  
  # calculate hs
  # working here, need to make a sample.pf option
  out <- .apply.snpR.facets(x, facets, "genotypes", .heterozygosity, case = "psamp", mDat = x@mDat, method = "hs")
  
  colnames(out)[which(colnames(out) == "stat")] <- "hs"
  x <- .merge.snpR.stats(x, out, "sample.stats")
  x <- .calc_weighted_stats(x, ofacets, "sample", "hs")
  x <- .update_calced_stats(x, facets, "hs")
  
  x <- .update_citations(x, keys = "Coltman1999", stats = "Hs", details = "Individual Heterozygosity")
  
  return(x)
}

# better code for the abba baba doc, but it currently isn't working

#\loadmathjax An ABBA/BABA test according to Green et al (2010) tests the
# hypothesis that there are an equal number of loci where population 1
# (\mjeqn{p_{1}}{ascii}) and population 2 (\mjeqn{p_{2}}{ascii}) are more
# closely related to population 3 (\mjeqn{p_3}{ascii}). The ratio of these two
# scenarios is given as \mjeqn{D = ABBA/BABA}{ascii}, where: \mjdeqn{ABBA = (1
# - p_{1})p_{2}p_{3}}{ascii}\mjdeqn{BABA = p_{1}(1 - p_{2})p_{3}}{ascii}
# where \mjeqn{p_{1}}{ascii},\mjeqn{p_{2}}{ascii}, and \mjeqn{p_{3}}{ascii} are
# the derived allele frequencies in populations 1 through 3, respectively.
# \emph{D} values are provided for both the overall comparison and within any
# levels of provided snp facets.
#
#  \mjeqn{D_{J}}{ascii} is
# weighted by the number of SNPs removed that block (\mjeqn{m_{j}}{ascii}),
# such that: \mjdeqn{D_{J} \sum_{j = 1}^{g}{D - D_{-j}} + \sum_{j =
# 1}^{g}{\frac{m_{j}D_{-j}}{n}}}{ascii} where \mjeqn{g}{ascii} is the number of
# blocks and \mjeqn{n}{ascii} is the total number of SNPs, and
# \mjeqn{D_{-j}}{ascii} is the D value for one block with the SNPs for that
# block's window ommited. The squared-standard error for \mjeqn{D_{J}}{ascii}
# is then: \mjdeqn{\sigma^{2} = \frac{1/6}\sum_{j = 1}^{g}{\frac{(\tau_{j} -
# \theta_{D})^{2}}{h_{j} - 1}}}{ascii}where \mjeqn{h_{j} =
# \frac{n}{m_{j}}}{ascii} and \mjdeqn{\tau_{j} = h_{j}D - ((h_j{} -
# 1)D_{-j})}{ascii} A \emph{Z} value can then be calculated as usual (\mjeqn{Z
# = \frac{D}{\sqrt{\sigma^{2}}}}{ascii}), and a \emph{p}-value determined from
# using a two-sided \emph{Z}-test following a normal distribution with
# \mjeqn{\mu = 0}{ascii} and \mjeqn{\sigma = 1}{ascii}.



#' Conduct an ABBA/BABA test for gene flow.
#'
#' Estimate D according to Green et al (2010) via an ABBA/BABA test for a
#' specific set of three different populations.
#'
#' An ABBA/BABA test according to Green et al (2010) tests the hypothesis that
#' there are an equal number of loci where population 1 and population 2  are
#' more closely related to population 3. The ratio of these two scenarios is
#' given as ABBA/BABA, where: ABBA = (1 - p1)p2p3 and BABA = p1(1 - p2)p3. where
#' p1, p2, and p3 are the derived allele frequencies in populations 1 through 3,
#' respectively. \emph{D} values are provided for both the overall comparison
#' and within any levels of provided snp facets.
#'
#' \emph{p}-values for \emph{D} can be calculated by a block jackknifing
#' approach according to Maier et. al (2022) by removing SNPs from each of
#' \emph{n} genomic windows (blocks) and then calculating \emph{D} for all other
#' sites. Non-overlapping windows are "blocked" by reference to any provided SNP
#' metadata facets (usually chromosome), each with a length equal to
#' \emph{sigma}, so providing \emph{sigma = 100} will use 100kb windows.
#'
#' @param x snpRdata. Input SNP data.
#' @param facet character. Categorical metadata variables by which to break up
#'   analysis. Must contain a sample facet. See \code{\link{Facets_in_snpR}} for
#'   more details.
#' @param p1 character. Name of population 1, must match a category present in
#'   the provided facet.
#' @param p2 character. Name of population 2, must match a category present in
#'   the provided facet.
#' @param p3 character. Name of population 3, must match a category present in
#'   the provided facet.
#' @param jackknife logical, default FALSE. If TRUE, block-jackknifed
#'   significance for D will be calculated with window size sigma according to
#'   any SNP facet levels. See details.
#' @param jackknife_par numeric or FALSE, default FALSE. If numeric, jackknifes
#'   per SNP levels will be run with the requested number of processing threads.
#' @param sigma numeric or NULL, default NULL. If jackknifes are requested, the
#'   size of the windows to block by, in kb. Should be large enough to account
#'   for linkage!
#'
#' @author William Hemstrom
#' @references
#' Maier, R. et al (2022). On the limits of fitting complex models of population
#' history to genetic data. BioRxiv. doi: 10.1101/2022.05.08.491072
#' 
#' Green, ER, et al (2010). A Draft Sequence of the Neandertal Genome. Science,
#' 328(5979), 710–722. doi: 10.1126/science.1188021
#' 
#' @export
#' @examples 
#' 
#' # add the ref and anc columns
#' # note: here these results are meaningless since they are arbitrary.
#' x <- stickSNPs
#' maf <- get.snpR.stats(x, ".base", stats = "maf")$single
#' snp.meta(x)$ref <- maf$major
#' snp.meta(x)$anc <- maf$minor
#' 
#' # run with jackknifing, 1000kb windows!
#' # Test if ASP or UPD have more geneflow with PAL...
#' x <- calc_abba_baba(x, "pop.chr", "ASP", "UPD", "PAL", TRUE, sigma = 1000)
#' get.snpR.stats(x, "pop.chr", "abba_baba") # gets the per chr results
#' get.snpR.stats(x, "pop", "abba_baba") # gets the overall results
#' 
#' # smoothed windowed averages
#' x <- calc_smoothed_averages(x, "pop.chr", sigma = 200, step = 200, 
#'    nk = TRUE, stats.type = "pairwise")
#' get.snpR.stats(x, "pop.chr", "abba_baba")
calc_abba_baba <- function(x, facet, p1, p2, p3, jackknife = FALSE, jackknife_par = FALSE, sigma = NULL){
  ..keep.names <- ..rm.cols <- abba <- baba <- snp_facet_D <- NULL
  
  #============sanity checks==============
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  if(length(facet) != 1){
    stop("Only one facet at a time can be run through calc_abba_baba due to p1, p2, p3 specification. Please run more facets manually for now.\n")
  }
  
  facet <- .check.snpR.facet.request(x, facet, remove.type = "none")
  sample.facet <- .check.snpR.facet.request(x, facet, "snp")
  if(sample.facet == ".base"){
    stop("A sample facet MUST be provided.\n")
  }
  
  x <- .add.facets.snpR.data(x, facet)
  
  check_sample_tasks <- .get.task.list(x, facet)
  bad_levels <- which(!c(p1, p2, p3) %in% check_sample_tasks[,2])
  if(length(bad_levels) > 0){
    stop(paste0("Some sample categories not found in sample metadata (", paste0(c(p1, p2, p3)[bad_levels], collapse = ", "), ").\n"))
  }
  
  if(jackknife){
    if(!is.numeric(sigma)){
      stop("If jackknifes are requested, sigma must be provided (windows will be sigma*1000 bp wide).\n")
    }
  }
  
  if(any(!c("ref", "anc") %in% colnames(snp.meta(x)))){
    stop("ref (derived) and anc (ancestral) alleles must be noted in SNP metadata.\n")
  }
  
  
  #============sub-functions==============
  jack_fun <- function(as, p1, p2, p3){
    
    # get derived allele counts by looping through each option
    der <- numeric(nrow(as))
    nucs <- c("A", "C", "G", "T")[which(c("A", "C", "G","T") %in% colnames(as))]
    nuc_cols <- which(colnames(as) %in% nucs)
    for(i in 1:length(nucs)){
      matching <- which(as$ref == nucs[i])
      der[matching] <- as[matching, nucs[i]] 
    }
    
    # get derived allele frequencies
    der <- der/rowSums(as[,nucs])
    
    # get nk for smoothing later if requested
    nk <- rowSums(as[which(as$subfacet == p1), nuc_cols]) + 
      rowSums(as[which(as$subfacet == p2), nuc_cols]) + 
      rowSums(as[which(as$subfacet == p3), nuc_cols])
    
    # get p1, p2, p3, pO
    p1 <- der[which(as$subfacet == p1)]
    p2 <- der[which(as$subfacet == p2)]
    p3 <- der[which(as$subfacet == p3)]
    
    # get abba and baba
    abba <- (1 - p1) * p2 * p3
    baba <- p1 * (1 - p2) * p3
    
    D_overall <- (sum(abba) - sum(baba)) / (sum(abba) + sum(baba))
    D_per_loci <- (abba - baba)/(abba + baba)
    
    
    return(list(overall = D_overall, 
                per_loci = data.table::as.data.table(cbind(unique(as[,-which(colnames(as) %in% c("subfacet", "A", "C", "T", "G"))]), 
                                                           data.frame(D = D_per_loci, abba = abba, baba = baba, nk = nk)))))
  }
  
  # assume this is being passed only the correct meta, runs for the given facet level. Note that for this approach, step and sigma should be identical.
  one_jack_boot <- function(as, on_lev, p1, p2, p3, sigma){
    starts <- seq(min(as$position[on_lev]), max(as$position[on_lev]), by = sigma*1000)
    ends <- starts + sigma*1000 - 1
    
    
    out <- data.frame(D = numeric(length(ends)), m = numeric(length(ends))) #initialize output
    
    # run the loop
    for (i in 1:length(starts)){
      matches <- intersect(which(as$position >= starts[i] & as$position <= ends[i]), on_lev)
      if(length(matches) > 0){
        out[i,1] <- jack_fun(as[-matches,], p1, p2, p3)$overall
        out[i,2] <- length(matches)/3
      }
      else{
        out[i,1] <- NaN
        out[i,2] <- 0
      }
    }
    
    if(any(is.nan(out[,1]))){out <- out[-which(is.nan(out[,1])),]}
    
    return(out)
  }
  
  # selects the correct metadata for one run, given a task list without column 2 (sample subfacet)
  select_correct_meta <- function(x, adjusted_task_row, p1, p2, p3){
    select <- x@facet.meta$facet == adjusted_task_row[1] & 
      x@facet.meta$subfacet %in% c(p1, p2, p3)
    
    if(adjusted_task_row[2] != ".base"){
      paste_cols <- .paste.by.facet(x@facet.meta, unlist(.split.facet(adjusted_task_row[2])), ".")
      
      select <- select & paste_cols == adjusted_task_row[3]
    }
    
    
    return(select)
  }
  
  bind_meta <- function(x, select){
    return(cbind(x@facet.meta[select,], x@geno.tables$as[select,]))
  }
  
  #============prepare input data=======
  sample.facet <- .check.snpR.facet.request(x, facet)
  x <- .add.facets.snpR.data(x, facet)
  
  selected_meta <- select_correct_meta(x, c(sample.facet, ".base", ".base"), p1, p2, p3)
  
  overall <- jack_fun(bind_meta(x, selected_meta), p1, p2, p3)
  snp.facet <- .check.snpR.facet.request(x, facet, "sample")
  if(snp.facet != ".base"){
    cnames <- unlist(.split.facet(snp.facet))
    per_snp_facet <- overall$per_loci[, snp_facet_D := (sum(abba) - sum(baba))/(sum(abba) + sum(baba)), by = cnames]
    keep.names <- c("facet", "facet.type", cnames, "snp_facet_D")
    .fix..call(per_snp_facet <- unique(per_snp_facet[,..keep.names]))
  }
  
  
  #============jackknife if wanted=======
  if(jackknife){
    
    # grab tasks (snp facet levels)
    tasks <- .get.task.list(x, facet)
    tasks <- tasks[,-2]
    tasks <- unique(tasks)
    
    # serial
    if(isFALSE(jackknife_par)){
      jack_res <- vector("list", nrow(tasks))
      
      for(i in 1:nrow(tasks)){
        select_part <- select_correct_meta(x, tasks[i,], p1, p2, p3)
        select_whole <- select_correct_meta(x, c(tasks[i,1], ".base", ".base"), p1, p2, p3)
        on_lev <- match(which(select_part), which(select_whole))
        jack_res[[i]] <- one_jack_boot(bind_meta(x, select_whole), on_lev, p1, p2, p3, sigma)
      }
      
    }
    
    # parallel
    else{
      cl <- parallel::makePSOCKcluster(jackknife_par)
      doParallel::registerDoParallel(cl)
      
      #prepare reporting function
      ntasks <- nrow(tasks)
      # progress <- function(n) cat(sprintf("Job %d out of", n), ntasks, "is complete.\n")
      # opts <- list(progress=progress)
      
      #loop through each set
      jack_res <- foreach::foreach(q = 1:ntasks,
                                   .packages = c("snpR", "data.table")) %dopar% {
                                     
                                     select_part <- select_correct_meta(x, tasks[q,], p1, p2, p3)
                                     select_whole <- select_correct_meta(x, c(tasks[q,1], ".base", ".base"), p1, p2, p3)
                                     on_lev <- match(which(select_part), which(select_whole))
                                     one_jack_boot(bind_meta(x, select_whole), on_lev, p1, p2, p3, sigma)
                                   }
      
      parallel::stopCluster(cl)
    }
    
    #============calculate the p value and delta j (jackknifed D estimate)===========
    jack_res <- dplyr::bind_rows(jack_res)
    
    # according to Maier et al 2022, see https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/wjack.pdf for formulae
    
    n <- sum(jack_res$m)
    g <- nrow(jack_res)
    hj <- n/jack_res$m
    
    delta_j <- sum(overall$overall - jack_res$D) + sum((jack_res$m*jack_res$D)/n)
    pseudos <- hj*overall$overall - ((hj - 1)*jack_res$D)
    sq_err <- (1/g) * sum(((pseudos - delta_j)^2)/(hj - 1))
    
    # Z score from std.err, from there to p-value... technically we should use delta_j here not D, so if we jackknife, our returned value should actually be the jackknifed estimate of D
    Z <- delta_j/sqrt(sq_err)
    D.pval <- 2*stats::pnorm(abs(Z), lower.tail = FALSE)
    overall$overall_boot <- delta_j
  }
  
  #=============merge and return===========
  
  # citations
  x <- .update_citations(x, keys = "E.2010", stats = "D", details = "D based on ABBA/BABA test.")
  if(jackknife){
    x <- .update_citations(x, keys = "Maier2022.05.08.491072", stats = "D", details = "Jackknife test for significance of D.")
  }
  
  # merge
  comp_lab <- paste0(p1, "~", p2, "~", p3)
  
  ## per_loci
  overall$per_loci$comparison <- comp_lab
  rm.cols <- which(colnames(overall$per_loci) %in% c("snp_facet_D", "facet.type"))
  colnames(overall$per_loci)[which(colnames(overall$per_loci) == "D")] <- "D_abba_baba"
  .fix..call(overall$per_loci <- overall$per_loci[,-..rm.cols])
  data.table::setcolorder(overall$per_loci, c("facet", "comparison", colnames(overall$per_loci)[-which(colnames(overall$per_loci) %in% c("facet", "comparison"))]))
  x <- .merge.snpR.stats(x, overall$per_loci, type = "pairwise")
  
  ## overall
  new_mean <- data.frame(facet = sample.facet,
                         subfacet = comp_lab,
                         snp.facet = ".base",
                         snp.subfacet = ".base",
                         D_abba_baba = overall$overall)
  if(jackknife){
    new_mean$D_abba_baba_jackknife <- overall$overall_boot
    new_mean$D_abba_baba_p_value <- D.pval
  }
  x <- .merge.snpR.stats(x, 
                         stats = new_mean,
                         type = "weighted.means")
  
  ## per_snp_facet
  if(snp.facet != ".base"){
    snp.facet_levs <- unlist(.split.facet(.check.snpR.facet.request(x, snp.facet, "sample")))
    snp.facet_levs <- .paste.by.facet(as.data.frame(per_snp_facet), snp.facet_levs)
    
    new_mean <- data.frame(facet = sample.facet,
                           subfacet = comp_lab,
                           snp.facet = snp.facet,
                           snp.subfacet = snp.facet_levs)
    new_mean$D_abba_baba <- per_snp_facet$snp_facet_D
    x <- .merge.snpR.stats(x, new_mean, "weighted.means")
  }
  
  x <- .update_calced_stats(x, facet, "abba_baba")
  
  
  return(x)
}

#' Calculate the proportion of polymorphic loci.
#'
#' Calculates the proportion of polymorphic (genotyped) loci for any combination
#' of SNP and sample facets. Note: per window counts of polymorphic loci can be
#' calculated using \code{\link{calc_tajimas_d}}.
#'
#' Note that per window counts of polymorphic loci can be calculated using
#' \code{\link{calc_tajimas_d}} (contained in the \code{num_seg} column).
#' Dividing this by the \code{n_snps} column will return the proportion of
#' polymorphic loci in each window if that is desired!
#' 
#' @param x snpRdata object
#' @param facets facets over which to check for polymorphic loci. Can be any
#'   combination of SNP/sample facets.
#'   See \code{\link{Facets_in_snpR}} for details.

#'
#' @return A snpRdata object with prop_poly merged into the weighted.means 
#'   slot.
#'   
#' @author William Hemstrom
#' @export
#' 
#' @examples
#' # base facet
#' x <- calc_prop_poly(stickSNPs)
#' get.snpR.stats(x, stats = "prop_poly")$weighted.means
#'
#' # multiple facets
#' x <- calc_prop_poly(stickSNPs, c("pop", "fam", "pop.fam", 
#'                                  "pop.fam.chr", "chr"))
#' get.snpR.stats(x, c("pop", "fam", "pop.fam", "pop.fam.chr", "chr"), 
#'                stats = "prop_poly")$weighted.means
#'
calc_prop_poly <- function(x, facets = NULL){
  snp.subfacet <- NULL
  #==================sanity checks====================================
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  x <- .add.facets.snpR.data(x, facets)
  facets <- .check.snpR.facet.request(x, facets, "none")
  
  #==================subfunction: run once for each snp facet!========
  func <- function(as){ 
    meta <- as[,which(!colnames(as) %in% colnames(x@geno.tables$as) & !colnames(as) %in% c("facet.type", ".snp.id"))]
    
    as <- as[, which(colnames(as) %in% colnames(x@geno.tables$as))]
    as <- ifelse(as != 0, TRUE, FALSE)
    meta$poly <- rowSums(as) > 1
    meta$genotypes <- rowSums(as) > 0
    if(any(!meta$genotypes)){
      meta <- meta[-which(!meta$genotypes),]
    }
    
    meta <- data.table::as.data.table(meta)
    tab <- meta[, .SDcols = "poly", lapply(.SD, function(x) sum(x)/length(x)), by = c(colnames(meta)[-which(colnames(meta) %in% c("poly", "genotypes"))])]
    
    sf.cols <- which(!colnames(tab) %in% c("facet", "subfacet", "poly"))
    if(length(sf.cols) > 0){
      sf.names <- paste0(colnames(tab)[sf.cols], collapse = ".")
      sf.names <- .check.snpR.facet.request(x, sf.names, "none")
      tab$snp.facet <- sf.names
      tab[, snp.subfacet := .paste.by.facet(as.data.frame(tab), .split.facet(sf.names)[[1]])]
    }
    else{
      tab$snp.facet <- ".base"
      tab$snp.subfacet <- ".base"
    }
    
    tab <- tab[,c("facet", "subfacet", "snp.facet", "snp.subfacet", "poly")]
    colnames(tab)[5] <- "prop_poly"
    return(tab)
  }
  
  #=========apply=========
  # get unique snp facets
  snp.levs <- lapply(facets, function(f) .check.snpR.facet.request(x, f, remove.type = "sample"))
  psl <- unlist(lapply(snp.levs, function(f) paste0(f, collapse = ".")))
  sample.levs <- lapply(facets, function(f) .check.snpR.facet.request(x, f, remove.type = "snp"))
  sample.levs <- lapply(sample.levs, function(f) paste0(f, collapse = "."))
  
  usl <- unique(snp.levs)
  
  output <- vector("list", length(usl))
  for(i in 1:length(usl)){
    sample.levs.to.run <- paste0(usl[[i]], collapse = ".")
    sample.levs.to.run <- unlist(sample.levs[which(psl == sample.levs.to.run)])
    matches <- which(x@facet.meta$facet %in% sample.levs.to.run)
    col.matches <- which(colnames(x@facet.meta) %in% c("facet", "subfacet", .split.facet(usl[[i]])[[1]]))
    output[[i]] <- func(cbind(x@facet.meta[matches, col.matches], x@geno.tables$as[matches,]))
  }
  
  output <- dplyr::bind_rows(output)
  
  #==========merge========
  
  x <- .update_calced_stats(x, facets, stats = "prop_poly")
  x <- .merge.snpR.stats(x, output, type = "weighted.means")
}


#' Generate phylogenetic-like clustering trees.
#'
#' Generate dendritic trees in the style of a phylogenetic tree for
#' individuals or groups of individuals from snpR data. Note that this function
#' is not overwrite safe.
#'
#' Trees are generated via the \code{\link[ape]{nj}} or \code{\link[ape]{bionj}}
#' ape package for nj or bionj trees. Plots of the resulting trees can be
#' produced using other plotting tools, such as the \code{geom_tree} function
#' from \code{ggtree}. These are not produced automatically because
#' \code{ggtree} can have unexpected outputs, but the process is straightforward
#' and examples are provided in the \emph{examples} section of this
#' documentation. For more information, see the documentation for those
#' functions and packages.
#'
#' Bootstraps are conducted by re-sampling SNPs with replacement, according to
#' Felsenstein (1985). If no snp level facets are provides, loci are resampled
#' without restraint. If a snp level facet is provided, loci are only resampled
#' within the levels of that facet (e.g. within chromosomes).
#'
#' The genetic distances used to make the trees are calculated using
#' \code{\link{calc_genetic_distances}}. If a sample facet is used, that
#' function uses code derived from \code{\link[adegenet]{adegenet}}. Please cite
#' them and the actual method (e.g. Edwards, A. W. F. (1971)) alongside the
#' tree-building approach.
#'
#' Bootstrapping is done via the boot.phylo function in the ape package, and as
#' such does not support parallel runs on Windows machines.
#'
#' @param x snpRdata object.
#' @param facets character or NULL, default NULL. Facets for which to calculate
#'   genetic distances, as described in \code{\link{Facets_in_snpR}}. If snp or
#'   base level facets are requested, distances will be between individuals.
#'   Otherwise, distances will be between the levels of the sample facets. If
#'   bootstraps are requested, SNPs will
#' @param distance_method character, default "Edwards". Name of the method to
#'   use. Options: \itemize{\item{Edwards} Angular distance as described in
#'   Edwards 1971.} See \code{\link{calc_genetic_distances}}.
#' @param interpolate character, default "bernoulli". Missing data interpolation
#'   method, solely for individual/individual distances. Options detailed in
#'   documentation for \code{\link{format_snps}}.
#' @param tree_method character, default nj. Method by which the tree is
#'   constructed from genetic distances. Options: \itemize{\item{nj}
#'   Neighbor-joining trees, via \code{\link[ape]{nj}}. \item{bionj} BIONJ
#'   trees, according to Gascuel 1997, via \code{\link[ape]{bionj}}.
#'   \item{upgma} UPGMA trees, via \code{\link[stats]{hclust}}.}
#' @param root character or FALSE, default FALSE. A vector containing the
#'   requested roots for each facet. Roots are specified by a string matching
#'   either the individual sample or sample facet level by which to root. If
#'   FALSE for a given facet, trees will be unrooted. Note that all UPGMA trees
#'   are automatically rooted, so this argument is ignored for that tree type.
#' @param boot numeric or FALSE, default FALSE. The number of bootstraps to do
#'   for each facet. See details.
#' @param boot_par numeric or FALSE, default FALSE. If a number, bootstraps will
#'   be processed in parallel using the supplied number of cores.
#' @param update_bib character or FALSE, default FALSE. If a file path to an
#'   existing .bib library or to a valid path for a new one, will update or
#'   create a .bib file including any new citations for methods used. Useful
#'   given that this function does not return a snpRdata object, so a
#'   \code{\link{citations}} cannot be used to fetch references.
#'
#' @author William Hemstrom
#' @export
#'
#' @return A nested, named list containing plots, trees, and bootstraps for each
#'   facet and facet level.
#'
#' @references Felsenstein, J. (1985). Confidence Limits on Phylogenies: An
#'   Approach Using the Bootstrap. Evolution, 39(4), 783–791.
#'   https://doi.org/10.2307/2408678
#'
#'   Gascuel, O. (1997).BIONJ: an improved version of the NJ algorithm based on
#'   a simple model of sequence data. Molecular Biology and Evolution, 14(7),
#'   685–695. https://doi.org/10.1093/oxfordjournals.molbev.a025808
#'
#'   Paradis, E., Claude, J. and Strimmer, K. (2004). APE: analyses of
#'   phylogenetics and evolution in R language. Bioinformatics, 20, 289–290.
#'   
#' 
#' @examples
#' # Calculate nj trees for the base facet, each chromosome, 
#' # and for each population.
#' 
#' # note: some versions of ggtree/ggplot2 give an error relating to unable
#' # to use 'fortify()' when attempting to plot that is potentially due to
#' # ape masking something when loading snpR. Saving the tree object as
#' # a .RDS object with 'saveRDS()', restarting R, loading the object in,
#' # then plotting without loading snpR seems to fix the issue.
#' 
#' \dontrun{
#' 
#' # make the trees
#' tp <- plot_tree(stickSNPs, c(".base", "pop", "chr"), 
#'                 root = c(FALSE, "PAL", FALSE))
#' 
#' # plot using the ggtree package. Not done internally due to unpredictable
#' # ggtree behavior
#' library(ggtree)
#' ggplot(tp$pop$.base, aes(x, y)) + 
#'   geom_tree() +
#'   geom_tiplab() + theme_tree()
#' 
#' # Calculate bionj trees for pop with bootstrapping
#' tp <- plot_tree(stickSNPs, "pop", root = "PAL", boot = 5)
#' ## plot
#' ggplot(tp$pop$.base, aes(x, y)) + 
#'   geom_tree() +
#'   geom_tiplab() + theme_tree() +
#'   geom_text2(ggplot2::aes(subset = !isTip, label = label),
#'                           hjust = -.3)
#' }
#' 
calc_tree <- function(x, facets = NULL, distance_method = "Edwards", 
                      interpolate = "bernoulli", 
                      tree_method = "nj", root = FALSE,
                      boot = FALSE, boot_par = FALSE, update_bib = FALSE){
  #=======sanity checks==========
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  facets <- .check.snpR.facet.request(x, facets, remove.type = "none")
  x <- .add.facets.snpR.data(x, facets)
  
  msg <- character(0)
  
  if(!tree_method %in% c("nj", "bionj", "upgma")){
    msg <- c(msg, "Tree method not recognized. Acceptable methods: nj, bionj, or upgma.\n")
  }
  
  if(length(root) < length(facets)){
    if(length(root) == 1){
      root <- rep(root, length(facets))
    }
    else{
      msg <- c(msg, "Root strategy not specified for each facet.\n")
    }
  }
  
  # if(.check.installed("ggtree", "github", "YuLab-SMU/ggtree")){
  #   if(utils::packageVersion("ggtree") < numeric_version("3.1.2")){
  #     msg <- c(msg, "Package ggtree version 3.1.2+ required. The most recent development version can be installed via the remotes package with remotes::install_github('YuLab-SMU/ggtree')")
  #   }
  # }
  
  if(!isFALSE(boot_par) & Sys.info()["sysname"] == "Windows"){
    msg <- c(msg, "plot_tree uses the ape package parallel setup, which does not work on a Windows platform.\n")
  }
  else if(isFALSE(boot_par)){
    boot_par <- 1
  }
  
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  .check.installed("ape")
  
  #=======function=========
  # expects 
  fun <- function(x, tfacet, distance_method, tree_method, interpolate, root, boot, boot_par){
    if(root == "FALSE"){root <- FALSE}
    #============define tree method===========
    if(tree_method == "nj"){
      if(!isFALSE(root)){
        tree_fun <- function(x, root = NULL, ...) ape::root(ape::nj(x), outgroup = root)
      }
      else{
        tree_fun <- function(x, root = NULL, ...) ape::nj(x)
      }
    }
    else if(tree_method == "bionj"){
      if(!isFALSE(root)){
        tree_fun <- function(x, root = NULL, ...) ape::root(ape::bionj(x), outgroup = root)
      }
      else{
        tree_fun <- function(x, root = NULL, ...) ape::bionj(x)
      }
    }
    else if(tree_method == "upgma"){
      if(isFALSE(root)){
        warning("UPGMA trees are always rooted.")
        root <- TRUE
      }     
      tree_fun <- function(x, root = NULL, ...) ape::as.phylo(stats::hclust(x, "average"))
    }
    
    
    #=============fetch data==============
    #=========raw afms, for non-snp facets===============
    snp.facet <- .check.snpR.facet.request(x, tfacet, "sample", fill_with_base = FALSE, return_base_when_empty = FALSE)
    samp.facet <- .check.snpR.facet.request(x, tfacet, "snp", fill_with_base = FALSE, return_base_when_empty = FALSE)
    if(!is.null(snp.facet) & is.null(samp.facet)){
      opts <- .get.task.list(x, snp.facet)
      amfs <- vector("list", nrow(opts))
      sn <- format_snps(x, "sn", interpolate = interpolate)
      sn <- sn[,-c(1:(ncol(x@snp.meta) - 1))]
      
      for(i in 1:nrow(opts)){
        tsn <- .fetch.snp.meta.matching.task.list(x, opts[i,])
        amfs[[i]] <- t(sn[tsn,])
      }
      names(amfs) <- opts[,4]
      is.snp.only <- TRUE
    }
    else if(tfacet == ".base"){
      sn <- format_snps(x, "sn", interpolate = interpolate)
      sn <- sn[,-c(1:(ncol(x@snp.meta) - 1))]
      amfs <- vector("list", 1)
      amfs[[1]] <- t(sn)
      names(amfs) <- ".base"
      is.snp.only <- TRUE
    }
    else{
      amfs <- get.snpR.stats(x, tfacet, "allele_frequency_matrix")[[1]]
      is.snp.only <- FALSE
    }
    #=============make a tree=============
    if(is.snp.only){
      tdf <- function(part, distance_method = NULL, ...) stats::dist(part)
    }
    else{
      tdf <- function(part, distance_method = NULL, ...) .get_dist(part, distance_method)[[1]]
    }
    
    tree <- lapply(amfs, function(y) tree_fun(tdf(y, distance_method), root))
    
    #=========bootstrap, if requested=============
    if(!isFALSE(boot)){
      
      #======define function========
      apply_boots <- function(tree, boots, root){
        boot_val <- ape::prop.clades(tree, boots, rooted = !isFALSE(root))
        boot_val[is.na(boot_val)] <- 0
        boot_val <- paste0((boot_val/length(boots))*100, "%") 
        tree$node.label <- boot_val
        return(tree)
      }
      boot_fun <- function(part, ref, tree_fun, distance_method, is.snp.only, root, B){
        boots <- ape::boot.phylo(phy = ref, x = part, FUN = function(y) tree_fun(tdf(y, distance_method), root = root), 
                                 trees = TRUE, mc.cores = boot_par, B = B, block = ifelse(is.snp.only, 1, 2))$trees
        
        res <- apply_boots(ref, boots, root)
        return(res)
      }
      
      
      
      #======run bootstraps=========
      cat("Generating Bootstraps for", tfacet, "...\n")
      for(i in 1:length(tree)){
        cat("Facet part", i, "out of", length(tree), "\n")
        tree[[i]] <- boot_fun(amfs[[i]], tree[[i]], tree_fun, distance_method, is.snp.only, root, boot)
      }
    }
    
    
    #=========generate plot==================
    # not doing this anymore, ggtree is quite unreliable and behaves oddly in interactive vs not, so just return the tree and have this
    # in the example...
    return(tree)
    
    
    # tout <- vector("list", length(tree))
    # names(tout) <- names(tree)
    # 
    # 
    # for(i in 1:length(tree)){
    #   # make the plot
    #   browser()
    #   tp <- ggplot2::ggplot(label = label) + 
    #     ggtree::geom_tree(data = tree[[i]], mapping = ggplot2::aes(x, y)) +
    #     ggtree::geom_tiplab() + ggtree::theme_tree()
    #   
    #   
    #   # add parts
    #   ## boot values
    #   if(!isFALSE(boot)){
    #     tp <- tp + ggtree::geom_text2(ggplot2::aes(subset = !isTip, label = label),
    #                                   hjust = -.3)
    #   }
    #   
    #   tout[[i]] <- list(plot = tp, tree = tree[[i]])
    # }
    # return(tout)
  }
  
  #=======run==============
  needed_dists <- .check_calced_stats(x, facets, paste0("genetic_distance", "_", distance_method, "_", interpolate))
  needed_dists <- which(!unlist(needed_dists))
  if(length(needed_dists) > 0){
    .make_it_quiet(
      x <- calc_genetic_distances(x, facets[needed_dists], distance_method, interpolate = interpolate))
  }
  
  out <- vector("list", length = length(facets))
  for(i in 1:length(facets)){
    cat("Facet:", facets[i], "\n")
    out[[i]] <- fun(x, facets[i], distance_method, tree_method, interpolate, root[i], boot, boot_par)
  }
  names(out) <- facets
  
  
  # update bib
  keys <- character()
  stats <- character()
  details <- character()
  
  if(distance_method == "Edwards"){
    keys <- c(keys, "Edwards1971")
    stats <- c(stats, "genetic_distance")
    details <- c(details, "Edwards' Angular Genetic Distance")
  }
  else if(distance_method == "Nei"){
    keys <- c(keys, "Nei1978")
    stats <- c(stats, "genetic_distance")
    details <- c(details, "Nei's genetic distance")
  }
  
  if(tree_method == "nj"){
    keys <- c(keys, "Felsenstein1985")
    stats <- c(stats, "tree")
    details <- c(details, "Neighbor-joining tree construction method")
  }
  else if(tree_method == "bionj"){
    keys <- c(keys, "Gascuel1997")
    stats <- c(stats, "tree")
    details <- c(details, "BIONJ tree construction method")
  }
  else if(tree_method == "upgma"){
    keys <- c(keys, "Sokal1958")
    stats <- c(stats, "tree")
    details <- c(details, "UPGMA tree construction method")
  }
  
  keys <- c(keys, "Paradis2004")
  stats <- c(stats, "tree")
  if(!isFALSE(boot)){
    details <- c(details, "tree construction and bootstrapping")
  }
  else{
    details <- c(details, "tree construction")
  }
  
  if(length(keys) > 0){
    .yell_citation(keys, stats, details, update_bib)
  }
  
  # return
  return(out)
}


calc_relatedness <- function(x, facets = NULL, methods = c("ql", "r")){
  as <- x@geno.tables$as
  as <- as/rowSums(as)
  x <- as.matrix(genotypes(x))
  x1 <- substr(x, 1, 1)
  x2 <- substr(x, 2, 2)
  
  x1f <- (1:ncol(as))[match(x1, colnames(as))]
  x1f <- matrix(as[((x1f-1)*nrow(x)) + (1:nrow(x))], nrow(x), ncol(x1))
  x2f <- (1:ncol(as))[match(x2, colnames(as))]
  x2f <- matrix(as[((x2f-1)*nrow(x)) + (1:nrow(x))], nrow(x), ncol(x2))
  
  res <- expand.grid(1:ncol(x), 1:ncol(x))
  res <- res[-which(res$Var1 == res$Var2),]
  
  # I'm confused on the multilocus wording of Wang et al
  for(i in 1:ncol(x)){
    for(j in (i + 1):ncol(x)){
      indicator <- as.numeric(x1[,i] == x1[,j]) +
        as.numeric(x1[,i] == x2[,j]) +
        as.numeric(x2[,i] == x1[,j]) +
        as.numeric(x2[,i] == x2[,j])
      
      left_num <-  indicator - 
        2*(x1f[i,] + x2f[i,])
      left_dom <- 2 * (1 + as.numeric(x1[,i] == x2[,i]) - x1f[,i] - x2f[,i])
      left <- left_num/left_dom
      
      
      right_num <- indicator - 
        2*(x1f[j,] + x2f[j,])
      right_dom <- 2 * (1 + as.numeric(x1[,j] == x2[,j]) - x1f[,j] - x2f[,j])
      right <- right_num/right_dom
      
      if(any(is.nan(left) | !is.finite(left))){
        left[is.nan(left) | !is.finite(left)] <- 0
      }
      
      if(any(is.nan(right) | !is.finite(right))){
        right[is.nan(right) | !is.finite(right)] <- 0
      }
      
      tot <- (sum(left) + sum(right))/2
      
      
    }
  }
  
}
