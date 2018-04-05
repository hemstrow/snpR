#Filters the data according to a set of filters. Vectorized to handle large datasets.
#Arguments:
# data: input data frame
#
#FILTERS:
#non_poly: filter out snps with only one homozygous genotype observed, list return is npf
# maf: filter out snps with a minor allele frequency below that given, list return is maff
# hf_hets: filters out data with a heterozygote frequency above the given value, list return is hfhf
# min_ind: filters out snps genotyped at less than the given number of individuals, list return is mif
# min_loci: Numeric value [0-1]: filters out INDIVIDUALS genotyped at less than this percentage of remaining loci.
#This version is vectorized, and is much faster!


#'Filter SNP data.
#'
#'\code{filter_snps} filters SNP data to remove loci which violate any of several assumptions and/or individuals which are sequenced at too few SNP loci.
#'
#'Description of x:
#'    Contains metadata in columns 1:ecs. Remainder of columns contain genotype calls for each individual. Each row is a different SNP, as given by format_snps output options 4 or 6.
#'
#'Possible filters:
#'\itemize{
#'    \item{maf, minor allele frequency: }{removes SNPs where the minor allele frequency is too low. Can look for mafs below #'provided either globally or search each population individually.}
#'    \item{hf_hets, high observed heterozygosity: }{removes SNPs where the observed heterozygosity is too high.}
#'    \item{min_ind, minimum individuals: }{removes SNPs that were genotyped in too few individuals.}
#'    \item{min_loci, minimum loci: }{removes individuals sequenced at too few loci.}
#'    \item{non_poly, non-polymorphic SNPs: }{removes SNPs that are not polymorphic (not true SNPs).}
#'    \item{bi_al, non-biallelic SNPs: }{removes SNPs that have more than two observed alleles.}
#'}
#'
#'Note that filtering out poorly sequenced individuals creates a possible conflict with the loci filters, since after individuals are removed, some loci may no longer pass filters. For example, if a portion of individuals in one population all carry the only instances of a rare minor allele that still passes the maf threshold, removing those individuals may cause the loci to no longer be polymorphic in the sample.
#'
#'To counter this, the "re_run" argument can be used to pass the data through a second filtering step after individuals are removed. By default, the "partial" re-run option is used, which re-runs only the non-polymorphic filter (if it was originally set), since these may cause downstream analysis errors. The "full" option re-runs all set filters. Note that re-running any of these filters may cause individuals to fail the individual filter after loci removal, and so subsequent tertiary re-running of the individual filters, followed by the loci filters, and so on, could be justified. This function stops after the second loci re-filtering, since that step is likely to be the most important to prevent downstream analytical errors.
#'
#'Via the "pop" argument, this function can filter by minor allele frequencies in either \emph{all} samples or \emph{in each population and the entire sample}. The latter should be used in instances where populaiton sizes are very different or there are \emph{many} populations, and thus common alleles of interest in one population might be otherwise filtered out. With very small populations, however, this may leave noise in the sample! In most cases, filtering the entire sample is sufficient.
#'
#' @param x data.frame. Input data, in the numeric or character format as given by format_snps options 4 or 6.
#' @param ecs Integer. Number of metadata columns at the start of x.
#' @param maf FALSE or numeric between 0 and 1, default FALSE. Minimum acceptable minor allele frequency
#' @param hf_hets FALSE or numeric between 0 and 1, default FALSE. Maximum acceptable heterozygote frequency.
#' @param min_ind FALSE or integer, default FALSE. Minimum number of individuals in which a loci must be sequenced.
#' @param min_loci FALSE or numeric between 0 and 1, default FALSE. Minimum proportion of SNPs an individual must be genotyped at.
#' @param re_run FALSE, "partial", or "full", default "partial". How should loci be re_filtered after individuals are filtered?
#' @param pop FALSE or list, default FALSE. A list with population information for individuals. Format is as produced by: list(c("North", "South", "East", "West"), c(10,20,30,40)). First vector is names of pops, second vector is the count of each pop. Input data MUST be in the same order as this list.
#' @param non_poly boolean, default FALSE. Should non-polymorphic loci be removed?
#' @param bi_al boolean, default TRUE. Should non-biallelic SNPs be removed?
#' @param mDat character variable, default "NN". Format of missing \emph{genotypes}. Overall data format is infered from this. Can be either "NN" or "0000".
#'
#' @return A data.frame in the same format as the input, with SNPs and individuals not passing the filters removed.
#'
#' @examples
#' #Strict filtering for missing individuals and unsequenced loci with partial re-run:
#' filter_snps(stickSNPs, 3, 0.05, 0.55, 250, .75)
#'
#' #Strict maf filtering with pops.
#' ##prep pop info
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' l <- list(c(names(pops)), as.numeric(pops))
#' ##filter
#' filter_snps(stickSNPs, 3, 0.05, 0.55, 250, .75, "full", pop = l)
#'
filter_snps <- function(x, ecs, maf = FALSE, hf_hets = FALSE, min_ind = FALSE,
                        min_loci = FALSE, re_run = "partial", pop = FALSE,
                        non_poly = TRUE, bi_al = TRUE, mDat = "NN"){

  #############################################################################################
  #do sanity checks

  if(maf){
    if(!is.numeric(maf)){
      stop("maf must be a numeric value.")
    }
    if(length(maf) != 1){
      stop("maf must be a numeric vector of length 1.")
    }
  }

  if(hf_hets){
    if(!is.numeric(hf_hets)){
      stop("hf_hets must be a numeric value.")
    }
    if(length(hf_hets) != 1){
      stop("hf_hets must be a numeric vector of length 1.")
    }
  }

  if(min_ind){
    if(!is.numeric(min_ind)){
      stop("min_ind must be a numeric value.")
    }
    if(length(min_ind) != 1){
      stop("min_ind must be a numeric vector of length 1.")
    }
  }

  if(min_loci){
    if(!is.numeric(min_loci) | (min_loci <= 0 | min_loci >= 1) | length(min_loci) != 1){
      stop("min_loci must be a numeric value between but not equal to 0 and 1.")
    }
  }

  if(is.list(pop)){
    if(!is.numeric(pop[[2]])){
      stop("pop[[2]] must be a numeric vector.")
    }
    if(sum(pop[[2]]) != (ncol(x) - ecs)){
      stop("sum of pop[[2]] must be equal to number of individuals in x (check ecs as well).")
    }
    if(length(pop[[1]]) != length(pop[[2]])){
      stop("Number of provided pop names and pop sizes not equal.")
    }
  }

  if(re_run != FALSE){
    if(re_run != "partial" & re_run != "full"){
      cat("re_run must be set to partial or full if not FALSE.\n")
    }
  }

  ##############################################################################################
  ###set up, get values used later, clean up data a bit, set any functions used lower
  cat("Initializing...\n")

  #get headers
  headers <- x[,1:ecs]

  #seperate out the genotypes alone
  x <- x[,(ecs+1):ncol(x)]

  #function to filter by loci, to be called before and after min ind filtering (if that is requested.)
  filt_by_loci <- function(){
    #get a single vector of genotypes
    xv <- as.vector(t(x))

    #determine snp format
    snp_form <- nchar(xv[1])

    #get number of samples
    nsamp <- ncol(x)

    #create tables of genotype counts for each locus.
    #function to do this, takes the pattern to match and the data, which is as a SINGLE VECTOR, by rows.
    #needs the number of samples from global function environment
    count_genos <- function(pattern, x){
      XXs <- grep(pattern, x) #figure out which entries match pattern
      out <- table(ceiling(XXs/nsamp)) #devide the entries that match by the number of samps to get which samp they come from (when rounded up), then table that.
      return(out)
    }

    #get all possible genotypes
    gs <- unique(xv)

    #which gs are heterozygous?
    hs <- substr(gs,1,snp_form/2) != substr(gs, (snp_form/2 + 1), snp_form*2)

    #which genotype is the missing data?
    mpos <- which(gs == mDat)

    #what are the possible alleles at all loci?
    as <- unique(c(substr(gs,1,snp_form/2), substr(gs, (snp_form/2 + 1), snp_form*2)))
    as <- as[as != substr(mDat, 1, snp_form/2)] #that aren't N?

    #############################################################################################3
    ###get a table of genotype and allele counts at each locus.
    cat("Creating genotype table...\n")

    #for each element of gs, get the tables of genotype counts and add them to a matrix
    gmat <- matrix(0, nrow(x), length(gs)) #initialize matrix
    colnames(gmat) <- gs #set the matrix names

    #fill the matrix, one possible genotype at a time (hard enough to vectorize as it is).
    for(i in 1:length(gs)){
      tab <- count_genos(gs[i], xv)
      gmat[as.numeric(names(tab)),i] <- as.numeric(tab)
    }

    if(length(mpos) > 0){
      tmat <- gmat[,-c(mpos)] #gmat without missing data
    }
    else{
      tmat <- gmat
    }

    #get matrix of allele counts
    #initialize
    cat("Getting allele table...\n")
    amat <- matrix(0, nrow(gmat), length(as))
    colnames(amat) <- as

    #fill in
    for(i in 1:length(as)){
      b <- grep(as[i], colnames(tmat))
      hom <- which(colnames(tmat) == paste0(as[i], as[i]))
      if(length(hom) == 0){
        het <- b
        amat[,i] <- rowSums(tmat[,het])
      }
      else{
        het <- b[b != hom]
        if(length(het) > 0){
          if(is.matrix(tmat[,het])){
            amat[,i] <- (tmat[,hom] * 2) + rowSums(tmat[,het])
          }
          else{
            amat[,i] <- (tmat[,hom] * 2) + tmat[,het]
          }
        }
        else{
          amat[,i] <- (tmat[,hom] * 2)
        }
      }
    }


    ##############################################################################################
    ###do filters
    keep <- logical(nrow(x)) #vector to track status

    #===============================================
    ##non-biallelic loci
    if(bi_al){
      cat("Filtering non-biallelic loci...\n")
      bimat <- ifelse(amat, TRUE, FALSE)
      bi <- ifelse(rowSums(bimat) > 2, T, F) #if true, should keep the allele
      #some tests require this, so subset the matrices and redefine things if true and some are multi-allelic
      if((maf | hf_hets) &  sum(bi) != 0){
        x <- x[bi,]
        xv <-   xv <- as.vector(t(x))
        gs <- unique(xv)
        hs <- substr(gs,1,snp_form/2) != substr(gs, (snp_form/2 + 1), snp_form*2)
        mpos <- which(gs == mDat)
        as <- unique(c(substr(gs,1,snp_form/2), substr(gs, (snp_form/2 + 1), snp_form*2)))
        as <- as[as != substr(mDat, 1, snp_form/2)] #that aren't N?
        gmat <- gmat[bi,]
        tmat <- tmat[bi,]
        amat <- amat[bi,]
        headers <- headers[bi,]
      }
      else{
        keep <- keep + bi
      }
    }

    #===============================================
    ##non-poly
    if(non_poly){
      cat("Filtering non-polymorphic loci...\n")
      npmat <- ifelse(tmat == 0, F, T) #convert to logical
      gcounts <- rowSums(npmat) #how many genotypes observed?
      hgcounts <- as.logical(rowSums(npmat[,hs[-mpos]])) #are hets observed?
      gcounts <- ifelse(gcounts == 1, F, T) #are there more than one genotype observed?
      np <- gcounts + hgcounts #figure out if either or both of the above are true
      np <- !as.logical(np) #convert to logical (if at least one is true, this will be false and we should keep the snp)
      keep <- keep + np
    }

    #===============================================
    ##min inds
    if(min_ind){
      cat("Filtering loci sequenced in few individuals...\n")
      mi <- ncol(x) - gmat[,colnames(gmat) == mDat]
      mi <- !mi >= min_ind #if false, enough samples, so keep.
      keep <- keep + mi
    }


    #===============================================
    ##minor allele frequency, both total and by pop. Should only run if bi_al = TRUE.
    if(maf){
      #if not filtering with multiple pops
      if(!is.list(pop)){
        cat("Filtering low minor allele frequencies, no pops...\n")
        mafs <- 1 - matrixStats::rowMaxs(amat)/rowSums(amat)
        mafs <- mafs < maf #less than required, set to true and reject.
        mafs[is.na(mafs)] <- TRUE
        keep <- keep + mafs
      }
      else{
        cat("Filtering low minor allele frequencies, pop:\n")
        pmafs <- logical(nrow(x))
        for(i in 1:length(pop[[1]])){
          cat(pop[[1]][i], "\n")
          #re-establish matrices with each pop

          #set the input data
          if(i == 1){
            popx <- x[,1:pop[[2]][i]]
          }
          else{
            popx <- x[,(sum(pop[[2]][1:(i-1)]) + 1):(sum(pop[[2]][1:i]))]
          }

          #get a single vector of genotypes
          popxv <- as.vector(t(popx))

          #get number of samples
          nsamp <- ncol(popx)

          #get all possible genotypes
          popgs <- unique(popxv)

          #which gs are heterozygous?
          pophs <- substr(popgs,1,snp_form/2) != substr(popgs, (snp_form/2 + 1), snp_form*2)

          #which genotype is the missing data?
          popmpos <- which(popgs == mDat)

          #what are the possible alleles at all loci?
          popas <- unique(c(substr(popgs,1,snp_form/2), substr(popgs, (snp_form/2 + 1), snp_form*2)))
          popas <- popas[popas != substr(mDat, 1, snp_form/2)] #that aren't N?

          #for each element of gs, get the tables of genotype counts and add them to a matrix
          popgmat <- matrix(0, nrow(popx), length(popgs)) #initialize matrix
          colnames(popgmat) <- popgs #set the matrix names

          #fill the matrix, one possible genotype at a time (hard enough to vectorize as it is).
          for(j in 1:length(popgs)){
            poptab <- count_genos(popgs[j], popxv)
            popgmat[as.numeric(names(poptab)),j] <- as.numeric(poptab)
          }

          if(length(popmpos) > 0){
            poptmat <- popgmat[,-c(popmpos)] #gmat without missing data
          }
          else{
            poptmat <- popgmat
          }

          #get matrix of allele counts
          #initialize
          popamat <- matrix(0, nrow(popgmat), length(popas))
          colnames(popamat) <- popas

          #fill in
          for(j in 1:length(popas)){
            b <- grep(popas[j], colnames(poptmat))
            hom <- which(colnames(poptmat) == paste0(popas[j], popas[j]))
            if(length(hom) == 0){
              het <- b
              popamat[,j] <- rowSums(poptmat[,het])
            }
            else{
              het <- b[b != hom]
              if(length(het) > 0){
                if(is.matrix(poptmat[,het])){
                  popamat[,j] <- (poptmat[,hom] * 2) + rowSums(poptmat[,het])
                }
                else{
                  popamat[,j] <- (poptmat[,hom] * 2) + poptmat[,het]
                }
              }
              else{
                popamat[,j] <- (poptmat[,hom] * 2)
              }
            }
          }

          #-------------------
          #established new set for this, now do maf as above.
          popmafs <- 1 - matrixStats::rowMaxs( popamat)/rowSums( popamat)
          popmafs <-  popmafs >= maf #if greater than requested, passes check for this pop
          popmafs[is.na(popmafs)] <- FALSE #call false, no information to go on here
          pmafs <- pmafs + popmafs #add this logical to the vector containing the sum of all such vectors
        }
        pmafs <- !as.logical(pmafs) #if false, we keep the allele since it was at >= maf in at least one pop
        #add in the overall mafs, since differential fixation will be removed by this!
        mafs <- 1 - matrixStats::rowMaxs(amat)/rowSums(amat)
        mafs <- mafs < maf #less than required, set to true and reject.
        mafs[is.na(mafs)] <- TRUE
        omafs <- (mafs + pmafs) == 2 #only keep those that failed both within pop and overal maf filters.

        keep <- keep + omafs
      }
    }

    #===============================================
    ##hf_hets. Should only run if bi_al = TRUE.
    if(hf_hets){
      cat("Filtering high frequency heterozygote loci...\n")
      #get heterozygote counts
      if(sum(hs[-mpos]) > 1){
        hetsum <- rowSums(tmat[,hs[-mpos]])
      }
      else if (sum(hs[-mpos]) == 0){
        hetsum <- numeric(nrow(x))
      }
      else{
        hetsum <- tmat[,hs[-mpos]]
      }
      hf <- hetsum/rowSums(tmat)
      hf <- hf >= hf_hets #if false, heterozygote frequency is lower than cut-off, keep locus
      keep <- keep + hf
    }

    keep <- !as.logical(keep)
    x <- x[keep,]

    return(list(x = x, headers = headers[keep,]))
  }

  #funciton to filter by individuals.
  min_loci_filt <- function(){
    cat("Filtering out individuals sequenced in few kept loci...\n")
    mcounts <- colSums(ifelse(x == mDat, 1, 0))
    rejects <- which(mcounts/nrow(x) >= (1 - min_loci))
    if(length(rejects) > 0){
      x <- x[,-rejects]
    }
    return(list(x = x, rejects = rejects))
  }

  ##############################################################################################
  ###call the functions as requested.

  if(any(c(non_poly, bi_al, maf, hf_hets, min_ind) != FALSE)){
    cat("Filtering loci. Starting loci:", nrow(x), "\n")
    out <- filt_by_loci()
    x <- out$x
    headers <- out$headers
    cat("\tEnding loci:", nrow(x), "\n")
  }

  if(min_loci){
    cat("Filtering individuals. Starting individuals:", ncol(x), "\n")
    out <- min_loci_filt()
    x <- out$x
    cat("\tEnding individuals:", ncol(x), "\n")
    if(re_run != FALSE){
      cat("Re-filtering loci...\n")

      if(re_run == "partial"){
        maf <- FALSE
        hf_hets <- FALSE
        min_ind <- FALSE
        bi_al <- FALSE
      }
      if(any(c(non_poly, bi_al, maf, hf_hets, min_ind) != FALSE)){
        rejects <- out$rejects
        if(is.list(pop)){ #remake pop if necissary after accounting for the removed samples
          s <- 1
          pop_temp <- pop
          for(i in 1:length(pop[[1]])){
            pop_temp[[2]][i] <- pop[[2]][i] - sum(rejects >= s & rejects < s + pop[[2]][i])
            s <- s + pop[[2]][i]
          }
          pop <- pop_temp
          empties <- which(pop[[2]] == 0)
          if(length(empties) > 0){
            pop <- list(pop[[1]][-empties], pop[[2]][-empties])
          }
        }
        out <- filt_by_loci() #re-filter loci to make sure that we don't have any surprise non-polys ect.
        x <- out$x
        headers <- out$headers
        cat("\tFinal loci count:", nrow(x), "\n")
      }
      else{
        cat("\tNo variables to re-fitler.\n")
      }
    }
  }

  #return results
  x <- cbind(headers, x)
  return(x)
}



#Converts snp/microsat data into different formats. Current supported formats are:
#Arguments:
# x: input data. Rows are loci, columns are individuals.
# ecs: number of extra columns at the header of the input (metadata, ect)
# output: output format:
#   1) per pop allele count format for Bayescan, ect. (v)
#   2) genepop format (v)
#   3) structure format (v)
#   4) 2 character numeric format (v)
#   5) migrate-n hapmap format (v)
#   6) NN format (from numeric) (v)
#   7) allele presence absensence format
# input_form: Input format
#   "NN": Default, alleles as single characters (as in AT or CG),
#   "0000": numeric. Will step to NN first before converting.
#   "msat_n": microsat, in two or three character format.
#           If in two character format, use "msat_2". If in three, use "msat_3"
#           Currently only supported for conversion to 7; 2, 3, and 4 forthcoming.
#   "snp_tab": Converts to NN first.
# miss: missing data format (per allele), typically "N", or "00"
# pop: For format 1 and 5.
#      Number of samples in each pop as a list (as in list(c("POP1", "POP2"), c(48,50)))).
#      Samples must be in the order given in input data frame.
# n_samp: number of randomly selected snps, ect to take. Can also take a numeric vector containing SNP indices to keep. For option 3.
# l_names: Vector of locus names. For option 7.
# interp_miss: Should missing data be interpolated? T or F. For option 8.
#note: (v) under format options means that I've already vectorized it.



#'Re-format genomic data.
#'
#'\code{format_snps} reformats SNP data into a range of different possible formats for use in snpR functions and elsewhere. Supports microsatellite data for conversion to a few formats as well.
#'
#'
#'Format options:
#'\itemize{
#'    \item{ac: }{allele count format, allele counts tabulated for all samples or within populations.}
#'    \item{genepop: }{genepop format, genotypes stored as four numeric characters (e.g. "0101", "0204"), transposed, and formatted for genepop. Rownames are individual IDs in genepop format, colnames are SNP ids, matching the first metadata column in input.}
#'    \item{structure: }{STRUCTURE format, two lines per individual: allele calls stored as single character numeric (e.g. "1", "2"). Allele calls per individual stored on two subsequent lines.}
#'    \item{numeric: }{numeric genotype tab format, genotypes stored as four numeric characters (e.g. "0101", "0204").}
#'    \item{hapmap: }{Migrate-n hapmap, allele counts tabulated within populations, in migrate-n hapmap format. Since this migrate-n implementation is iffy, this probably shouldn't be used much.}
#'    \item{character: }{character genotype tab format, genotypes stored as actual base calls (e.g. "AA", "CT").}
#'    \item{pa: }{allele presence/absence format, presence or absence of each possible allele at each possible genotype noted. Interpolation possible, with missing data substituted with allele freqency in all samples or each population.}
#'    \item{rafm: }{RAFM format, two allele calls at each locus stored in subsequent columns, e.g. locus1.1 locus1.2.}
#'    \item{faststructure: }{fastSTRCTURE format, identical to STRUCTURE format save with the addition of filler columns proceeding data such that exactly 6 columns proceed data. These columns can be filled with metadata if desired.}
#'}
#'
#'Example datasets in each format are available in \code{\link{stickFORMATs}} in elements named for output options.
#'
#'Input formats: For now, all input formats require at least two metadata columns.
#'\itemize{
#'    \item{NN: }{SNP genotypes stored as actual base calls (e.g. "AA", "CT").}
#'    \item{0000: }{SNP genotypes stored as four numeric characters (e.g. "0101", "0204").}
#'    \item{msat_2: }{Microsatellite genotypes stored as four numeric characters (e.g. "0101", "2740").}
#'    \item{msat_3: }{Microsatellite genotypes stored as six numeric characters (e.g. "120123", "233235").}
#'    \item{snp_tab: }{SNP genotypes stored with genotypes in each cell, but only a single nucleotide noted if homozygote and two nucleotides seperated by a space if heterozygote (e.g. "T", "T G").}
#'}
#'
#'Currently, msat_2 and msat_3 only support conversion to output option 7. 2, 3, and 4 are forthcoming.
#'
#'
#' @param x data.frame. Input data, in any of the above listed input formats.
#' @param ecs Integer. Number of extra metadata columns at the start of x.
#' @param output Character. Which of the output format should be used?
#' @param input_form Character, default "NN". Which of the above input formats should be used (e.g. "NN", "msat_2")?
#' @param miss Character, default "N". The coding for missing \emph{alleles} in x (typically "N" or "00").
#' @param pop List or integer 1, default 1. Population information for individuals. Format is as produced by: list(c("North", "South", "East", "West"), c(10,20,30,40)). First vector is names of pops, second vector is the count of each pop. Input data MUST be in the same order as this list. If 1, assumes one population. For output options 1 and 5 or 2 if outfile requested.
#' @param n_samp Integer or numeric vector, default NA. For output option 3. How many random loci should be selected? Can either be an integer or a numeric vector of loci to use.
#' @param interp_miss boolean, default TRUE. For output option 7. Should missing data be interpolated? Needed for PCA, ect.
#' @param lnames character vector, default NULL. For output option 7, optional vector of locus names by which to name output columns. If not provided, will use 1:nrow(x).
#' @param outfile character vector, default FALSE. If provided, the path to the output file to write to.
#'
#' @return A data.frame in the specified format.
#'
#' @examples
#' #allele count with pops:
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' l <- list(c(names(pops)), as.numeric(pops))
#' format_snps(stickSNPs, 3, "ac", pop = l)
#'
#' #genepop:
#' format_snps(stickSNPs, 3, "genepop")
#'
#' #STRUCTURE, subsetting out 100 random alleles:
#' format_snps(stickSNPs, 3, "structure", n_samp = 100)
#'
#' #STRUCTURE, subseting out the first 100 alleles:
#' format_snps(stickSNPs, 3, "structure", n_samp = 1:100)
#'
#' #fastSTRUCTURE
#' format_snps(stickSNPs, 3, "faststructure")
#'
#' #numeric:
#' format_snps(stickSNPs, 3, "numeric")
#'
#' #hapmap for migrate-n:
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' l <- list(c(names(pops)), as.numeric(pops))
#' format_snps(stickSNPs, 3, "hapmap", pop = l)
#'
#' #character:
#' num <- format_snps(stickSNPs, 3, 4)
#' format_snps(num, 3, "character", input_form = "0000", miss = "00")
#'
#' #presence/absence, SNP data:
#' format_snps(stickSNPs, 3, "pa")
#'
#' #presence/absence, 3 character microsat data (2 character is very similar):
#' format_snps(sthMSATS[seq(1, 13, by = 4),], 3, "pa", input_form = "msat_3", miss = "000")
#'
#' #RAFM, taking only 100 random snps.
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' l <- list(c(names(pops)), as.numeric(pops))
#' format_snps(stickSNPs, 3, "rafm", pop = l, n_samp = 100)
#'
#'
#'
format_snps <- function(x, ecs, output = 1, input_form = "NN",
                        miss = "N", pop = 1, n_samp = NA,
                        interp_miss = T, lnames = NULL, outfile = FALSE){

  #redefine x as "data", since this was originally written a while ago and already contains a variable called "x"
  data <- x
  rm(x)


  #possible outputs, recode to numbers so that I don't have to rewrite the whole damn thing. Still allowing numbers for reverse compatiblity.
  if(!is.numeric(output)){
    pos_outs <- c("ac", "genepop", "structure", "numeric", "hapmap", "character", "pa",
                  "rafm", "faststructure")
    if(!(output %in% pos_outs)){
      stop("Unaccepted output format specified. Check documentation. Remember, no capitalization!")
    }

    #reformat output to a number corresponding to format.
    output <- which(pos_outs == output)
  }


  if(ecs < 2){
    stop("Cannot have less than 2 metadata (ecs) columns.")
  }

  #####################################################################
  #do checks, print info
  if(output == 1){
    cat("Converting to per pop allele count format.\n")
    if(input_form != "0000" & input_form != "NN"){
      stop("Only 0000 and NN formats accepted.")
    }
    if(is.list(pop)){
      cat("Pops:")
      for (i in 1:length(unlist(pop[1]))){
        cat("\n\t", unlist(pop[1])[i], ", size: ", unlist(pop[2])[i])
      }
      cat("\n")
    }
  }
  else if(output == 2){
    cat("Converting to genepop format.\n")
    if(input_form != "0000" & input_form != "NN" & input_form != "snp_tab"){
      stop("Only 0000, NN, and snp_tab formats accepted for now.")
    }
  }
  else if(output == 3 | output == 9){
    if(output == 3){cat("Converting to STRUCTURE format.\n")}
    else{cat("Converting to fastSTRUCTURE format.\n")}
    if(input_form != "0000" & input_form != "NN"){
      stop("Only 0000 and NN formats accepted for now.")
    }
    if(length(n_samp) > 1){
      if(is.integer(n_samp)){
        cat("Number of designated sub-samples to take:", length(n_samp), "\n")
      }
      else{
        stop("Number of sub-samples to take must be a positive integer vector.\n")
      }
    }
    else if (is.numeric(n_samp)){
      if(floor(n_samp) == n_samp & n_samp > 0){
        cat("Number of designated sub-samples to take:", n_samp, "\n")
      }
      else{
        stop("Number of sub-samples to take must be a positive integer vector.\n")
      }
    }
    else if (!is.na(n_samp)){
      stop("Number of sub-samples to take must be a positive integer vector.\n")
    }
  }
  else if(output == 4){
    cat("Converting to numeric 2 character format.\n")
    if(input_form != "NN" & input_form != "snp_tab"){
      stop("Only NN and snp_tab formats accepted for now.")
    }
  }
  else if(output == 5){
    cat("Converting to migrate-N hapmap format.\n")
    if(input_form != "NN" & input_form != "0000"){
      stop("Only 0000 and NN formats accepted.")
    }
    if(is.list(pop)){
      cat("Pops:")
      for (i in 1:length(unlist(pop[1]))){
        cat("\n\t", unlist(pop[1])[i], ", size: ", unlist(pop[2])[i])
      }
      cat("\n")
    }
  }
  else if(output == 6){
    cat("Converting to NN format.\n")
    if(input_form != "0000" & input_form != "snp_tab"){
      stop("Only 0000 and snp_tab formats accepted.")
    }
  }

  else if(output == 7){
    cat("Converting to allele presence/absense format.\n")
  }

  else if(output == 8){
    if(input_form != "0000" & input_form != "snp_tab" & input_form != "NN"){
      stop("Only SNP data accepted for RAFM format.")
    }
    cat("Converting to RAFM format.\n")
  }

  else{
    stop("Please specify output format.")
  }

  if(input_form == "NN"){
    cat("Input format: NN\n")
    if(nchar(miss) != 1){
      stop("Missing data format must be one character.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "0000"){
    cat("Input format: 0000\n")
    if(nchar(miss) != 2){
      stop("Missing data format must be two characters.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "msat_2"){
    cat("Input format: msat, 2 character (0000)\n")
    if(nchar(miss) != 2){
      stop("Missing data format must be two characters.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "msat_3"){
    cat("Input format: msat, 3 character (000000)\n")
    if(nchar(miss) != 3){
      stop("Missing data format must be three characters.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "snp_tab"){
    cat("Input format, snp_tab.\n")
    if(nchar(miss) != 1){
      stop("Missing data format must be one character.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else{
    stop("Unsupported input format.")
  }


  if(outfile != FALSE){
    if(is.character(outfile) & length(outfile) == 1){
      cat("Printing results to:", outfile, "\n")
      if(file.exists(outfile)){
        #ask for confirmation before proceeding, since shit will be overwritten.
        cat("Outfile already exits. ")
        resp <- "empty"
        while(resp != "y" & resp != "n"){
          cat("Overwrite? (y or n)\n")
          resp <- readLines(n = 1)
        }
        if(resp == "n"){
          stop("Please provide acceptable path to file for output.\n")
        }
        else{
          cat("\tProceeding with conversion...\n")
        }
      }
    }
    else{
      stop("Outfile must be a character vector of length 1.\n")
    }
  }


  #####################################################################
  #function to create table of allele and genotype counts from all data. From filter_snps
  #input is a long vector of the data xv,
  #snp_form, the length of each genotype.
  #mDat, the missing data format
  #nsap, the number of samples
  #nloci, the number of loci
  tabulate_genotypes <- function(xv, snp_form, mDat, nsamp, nloci){
    #create tables of genotype counts for each locus.
    #function to do this, takes the pattern to match and the data, which is as a SINGLE VECTOR, by rows.
    #needs the number of samples from global function environment
    count_genos <- function(pattern, x){
      XXs <- grep(pattern, x) #figure out which entries match pattern
      out <- table(ceiling(XXs/nsamp)) #devide the entries that match by the number of samps to get which samp they come from (when rounded up), then table that.
      return(out)
    }

    #get all possible genotypes
    gs <- unique(xv)

    #which gs are heterozygous?
    hs <- substr(gs,1,snp_form/2) != substr(gs, (snp_form/2 + 1), snp_form*2)

    #which genotype is the missing data?
    mpos <- which(gs == mDat)

    #what are the possible alleles at all loci?
    as <- unique(c(substr(gs,1,snp_form/2), substr(gs, (snp_form/2 + 1), snp_form*2)))
    as <- as[as != substr(mDat, 1, snp_form/2)] #that aren't N?

    ##############################################################################################
    ###get a table of genotype and allele counts at each locus.
    cat("Creating genotype table...\n")

    #for each element of gs, get the tables of genotype counts and add them to a matrix
    gmat <- matrix(0, nloci, length(gs)) #initialize matrix
    colnames(gmat) <- gs #set the matrix names
    #fill the matrix, one possible genotype at a time (hard enough to vectorize as it is).
    for(i in 1:length(gs)){
      tab <- count_genos(gs[i], xv)
      gmat[as.numeric(names(tab)),i] <- as.numeric(tab)
    }

    if(length(mpos) > 0){
      tmat <- gmat[,-c(mpos)] #gmat without missing data
    }
    else{
      tmat <- gmat
    }

    #get matrix of allele counts
    #initialize
    cat("Getting allele table...\n")
    amat <- matrix(0, nrow(gmat), length(as))
    colnames(amat) <- as

    #fill in
    for(i in 1:length(as)){
      b <- grep(as[i], colnames(tmat))
      hom <- which(colnames(tmat) == paste0(as[i], as[i]))
      if(length(hom) == 0){
        het <- b
        amat[,i] <- rowSums(tmat[,het])
      }
      else{
        het <- b[b != hom]
        if(length(het) > 0){
          if(is.matrix(tmat[,het])){
            amat[,i] <- (tmat[,hom] * 2) + rowSums(tmat[,het])
          }
          else{
            amat[,i] <- (tmat[,hom] * 2) + tmat[,het]
          }
        }
        else{
          amat[,i] <- (tmat[,hom] * 2)
        }
      }
    }
    return(list(gs = tmat, as = amat))
  }

  #####################################################################
  #conversions

  #convert 0000 to NN unless output is 7 or 2. Return after if output is 6. (v)
  if(input_form == "0000" & output != 7 & output != 2){

    #Do conversion
    cat("Converting genotypes to NN form intermediary...")

    #vectorize and replace
    xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
    xv1 <- substr(xv, 1, 2)
    xv2 <- substr(xv, 3, 4)
    xv1[xv1 == "01"] <- "A"
    xv1[xv1 == "02"] <- "C"
    xv1[xv1 == "03"] <- "G"
    xv1[xv1 == "04"] <- "T"
    xv1[xv1 == miss] <- "N"
    xv2[xv2 == "01"] <- "A"
    xv2[xv2 == "02"] <- "C"
    xv2[xv2 == "03"] <- "G"
    xv2[xv2 == "04"] <- "T"
    xv2[xv2 == miss] <- "N"
    xv <- paste0(xv1, xv2)

    #rebind to matrix and remake data.
    xv <- matrix(xv, nrow(data), (ncol(data) - ecs), T)
    cnames <- colnames(data)
    data <- cbind(data[,1:ecs], as.data.frame(xv, stringsAsFactors = F))
    colnames(data) <- cnames

    cat("\nMoving on to conversion...", "\n")
    miss <- "N" #reset miss to the correct entry
    if(output == 6){
      rdata <- data #all done if just converting to NN
    }
  }

  else if (input_form == "NN"){
    cat("Ensure that all data columns are character vectors!\n")
  }

  #convert snp_tab to NN (v)
  if(input_form == "snp_tab"){
    header <- data[,1:ecs]
    xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
    nsamp <- ncol(data) - ecs
    snames <- colnames(data)[(ecs+1):ncol(data)]
    rm(data)
    ptf <- nchar(xv) == 1
    xv[ptf] <- paste0(xv[ptf], xv[ptf]) #double up homozygotes
    remove(ptf)
    xv <- gsub(" ", "", xv) #combine heterozygotes.
    xv[xv == paste0(miss, miss)] <- "NN" #replace with correct missing data
    xv <- as.data.frame(matrix(xv, nrow(header), nsamp, byrow = T), stringsAsFactors = F)
    colnames(xv) <- snames
    out <- cbind(header, xv)
    if(output == 6){
      rdata <- out
    }
    else{
      data <- out
    }
  }


  ##convert to allele count or migrate-n format, migrate-n should ALWAYS have multiple pops (why else would you
  ##use it?) (v)
  if(output == 1 | output == 5){
    if(output == 5){cat("WARNING: Data does not have header or pop spacer rows.\n")}
    w_df <- data[,1:ecs] #intiallize w_df if doing option 1

    ##Make a function to sum each row
    sum_row <- function(row){
      l <- c(substr(row,1,1), substr(row,2,2))
      t <- table(l)
      return(t) #returns a table with counts of the two alleles and then an N
    }
    if (is.list(pop) == FALSE){
      if(output == 5){stop("Cannot convert to migrate-n format with only one pop.\n")}
      #prepare a genotype table
      xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
      tabs <- tabulate_genotypes(xv, nchar(xv[1]), paste0(miss, miss), ncol(data) - ecs, nrow(data))
      w_df$n_total <- rowSums(tabs$as)
      counts <- rowSums(ifelse(tabs$as == 0, 0, 1))
      if(any(counts != 2)){
        if(any(counts > 2)){
          vio <- which(counts > 2)
          warning("More than two alleles detected at some loci.\n List of violating loci...")
          return(vio)
        }
        else{
          vio <- which(counts != 2)
          warning("Some loci are not bi-allelic, see printed list.")
          print(vio)
        }
      }
      w_df$n_alleles <- rowSums(tabs$as != 0)
      #fill in ni1 and ni2.
      w_df$ni1 <- tabs$as[,"A"]
      w_df$ni1[w_df$ni1 == 0] <- tabs$as[which(w_df$ni1== 0), "G"]
      w_df$ni1[w_df$ni1 == 0] <- tabs$as[which(w_df$ni1== 0), "C"]
      w_df$ni2 <- tabs$as[,"T"]
      w_df$ni2[w_df$ni2 == 0] <- tabs$as[which(w_df$ni2== 0), "C"]
      w_df$ni2[w_df$ni2 == 0] <- tabs$as[which(w_df$ni2== 0), "G"]

      #fix anything that's not bi-allelic
      w_df$ni1[w_df$n_alleles == 1 & w_df$ni1 == 0] <- w_df$ni2[w_df$n_alleles == 1 & w_df$ni1 == 0]
      w_df$ni2[w_df$n_alleles == 1] <- 0
    }
    else{
      pop_count <- length(unlist(pop[1]))
      pop_sizes <- unlist(pop[2])
      if(sum(unlist(pop[2])) != (ncol(data) - ecs)){#checks to make sure that the pop sizes add up to agree
        #with data
        cat("Supplied population numbers do not equal the number of supplied loci columns. Exiting.")
        stop()
      }
      #initialize w_df for multiple pops
      #print(unlist(pop[1]))
      j <- 2
      w_df$pop <- unlist(pop[1])[1] #create pop column and set first pop name
      wa_df <- w_df #set first temp w_df
      for (j in j:pop_count){ #loop through each pop name
        wb_df <- wa_df #create temp df copy of wa_df
        #print(j)
        #print(unlist(pop[1])[j])
        wb_df$pop <- unlist(pop[1])[j]#change pop to current name
        w_df <- rbind(w_df, wb_df) #bind this copy to w_df
      }
      remove(wa_df, wb_df) #remove temp df

      ##loop through and count alleles for every row, splitting by pop

      #create allele count table for each pop, fill their section of data

      #build allele tables for each locus
      pop_as <- list()
      pals <- matrix(FALSE, nrow(data), 4)
      colnames(pals) <- c("A", "C", "G", "T")
      current <- ecs + 1
      for(i in 1:pop_count){
        #get the data for this population
        cat("Generating population allele count matrices, current pop:", pop[[1]][i], "\n")
        ne <- current + pop_sizes[i] - 1
        x <- data[,current:ne]
        xv <- as.vector(t(x))
        pop_as[[i]] <- tabulate_genotypes(xv, nchar(xv[1]), paste0(miss, miss), pop_sizes[i], nrow(data))$as
        pop_as[[i]] <- pop_as[[i]][, order(colnames(pop_as[[i]]))]

        #if somehow not all four possible alleles found in ALL DATA for this pop (small?), fill in table
        if(!"A" %in% colnames(pop_as[[i]])){pop_as[[i]] <- cbind(A = numeric(nrow(pop_as)), pop_as[[i]])}
        if(!"C" %in% colnames(pop_as[[i]])){pop_as[[i]] <- cbind(pop_as[[i]][,1], C = numeric(nrow(pop_as[[i]])), pop_as[[i]][,-1])}
        if(!"G" %in% colnames(pop_as[[i]])){pop_as[[i]] <- cbind(pop_as[[i]][,1:2], G = numeric(nrow(pop_as[[i]])), pop_as[[i]][,-c(1,2)])}
        if(!"T" %in% colnames(pop_as[[i]])){pop_as[[i]] <- cbind(pop_as[[i]][,1:3], T = numeric(nrow(pop_as[[i]])))}

        #figure out which alleles are present at each locus to write these.
        pals <- pals + (pop_as[[i]] != 0)
        current <- current + pop_sizes[i]

      }

      #now write the correct counts
      cat("Tabulating and writing results...\n")
      pals <- pals != 0
      if(any(rowSums(pals) != 2)){
        if(any(rowSums(pals) > 2)){
          stop("Some loci have more than three alleles. Check/Filter data!")
        }
        else{
          warning("Some loci have less than two alleles.\n")
        }
      }

      #vector which says which are fixed
      fixed <- rowSums(pals) == 1


      #Convert into vector which says which elements to keep.
      if(any(fixed)){
        palsv <- as.vector(t(pals[fixed == FALSE,]))
      }
      else{
        palsv <- as.vector(t(pals))
      }

      w_df$n_total <- NA
      w_df$n_alleles <- NA
      w_df$ni1 <- NA
      w_df$ni2 <- NA


      #loop through and print data
      for(i in 1:pop_count){
        #have to do this slightly differently if there are any non-polymorphic bases
        if(any(fixed) > 0){
          #do the fixed ones first
          w_df[w_df$pop == pop[[1]][i] & fixed,]$n_total <- rowSums(pop_as[[i]][fixed,])
          w_df[w_df$pop == pop[[1]][i] & fixed,]$n_alleles <- 1 #in this case it is fixed in all pops
          w_df[w_df$pop == pop[[1]][i] & fixed,]$ni1 <- rowSums(pop_as[[i]][fixed,])
          w_df[w_df$pop == pop[[1]][i] & fixed,]$ni2 <- 0

          #do the others
          if(any(!fixed)){
            av <- as.vector(t(pop_as[[i]][!fixed,]))
            av <- av[palsv]
            av <- matrix(av, nrow(data[!fixed,]), 2, byrow = T)
            w_df[w_df$pop == pop[[1]][i] & !fixed,]$n_total <- rowSums(av)
            w_df[w_df$pop == pop[[1]][i] & !fixed,]$n_alleles <- 2 #not fixed in at least one pop!
            w_df[w_df$pop == pop[[1]][i] & !fixed,]$ni1 <- av[,1]
            w_df[w_df$pop == pop[[1]][i] & !fixed,]$ni2 <- av[,2]
          }
        }

        else{
          #compare palsv to allele tables to keep correct alleles.
          av <- as.vector(t(pop_as[[i]]))
          av <- av[palsv]
          av <- matrix(av, nrow(data), 2, byrow = T)
          w_df[w_df$pop == pop[[1]][i],]$n_total <- rowSums(av)
          w_df[w_df$pop == pop[[1]][i],]$n_alleles <- 2 #not fixed in at least one pop!
          w_df[w_df$pop == pop[[1]][i],]$ni1 <- av[,1]
          w_df[w_df$pop == pop[[1]][i],]$ni2 <- av[,2]
        }
      }

      row.names(w_df) <- 1:nrow(w_df) #fix row names
      if(output == 5){ #if doing output style 5...
        w_df <- w_df[,c(1, (ecs + 1):ncol(w_df))] #remove extra columns except the first and pop columns
      }
    }
    rdata <- w_df
  }



  ##convert to genepop or numeric format (v)
  if (output == 2 | output == 4){
    cat("WARNING: For output option 3, data.frame output does not have population seperators or header information. Outfile, if requested and pop list provided, will.", "\n", "Converting genotypes...")

    if(input_form == "0000" & output == 2){
      w_df <- as.data.frame(t(data[,(ecs + 1):ncol(data)]), stringsAsFactors = F) #remove extra columns and transpose data
      row.names(w_df) <- paste0(row.names(w_df), " ,") #adding space and comma to row names, as required.
      rdata <- w_df
    }
    else if (input_form == "0000" & output == 4){
      stop("Same input form given as output selected.")
    }

    else{
      #vectorize and replace
      xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
      xv1 <- substr(xv, 1, 1)
      xv2 <- substr(xv, 2, 2)
      xv1[xv1 == "A"] <- "01"
      xv1[xv1 == "C"] <- "02"
      xv1[xv1 == "G"] <- "03"
      xv1[xv1 == "T"] <- "04"
      xv1[xv1 == miss] <- "00"
      xv2[xv2 == "A"] <- "01"
      xv2[xv2 == "C"] <- "02"
      xv2[xv2 == "G"] <- "03"
      xv2[xv2 == "T"] <- "04"
      xv2[xv2 == miss] <- "00"
      xv <- paste0(xv1, xv2)

      xv <- matrix(xv, nrow(data), (ncol(data) - ecs), T)
      inds <- colnames(data)
      data <- cbind(data[,1:ecs], as.data.frame(xv, stringsAsFactors = F))
      colnames(data) <- inds

      cat("\n", "Cleaning up...", "\n")
      if(output == 2){ #convert to genepop
        w_df <- as.data.frame(t(data[,(ecs + 1):ncol(data)]), stringsAsFactors = F) #remove extra columns and transpose data
        row.names(w_df) <- paste0(row.names(w_df), " ,") #adding space and comma to row names, as required.
        rdata <- w_df
      }
      else {#prepare numeric output, otherwise same format
        rdata <- data
      }
    }
  }


  ##convert to structure, fastStructure or RAFM format (v)
  if (output == 3 | output == 8 | output == 9){

    #subset if requested
    if(all(!is.na(n_samp))){
      cat("Subsetting ")
      if(length(n_samp) > 1){
        cat("designated SNPs.\n")
        data <- data[n_samp,]
      }
      else{
        cat(n_samp, " random SNPs.\n")
        data <- data[sample(nrow(data), n_samp, T),]
      }
    }

    #plop actual data into two vectors, one for each allele.
    xv <- as.vector(t(t(data[,(ecs + 1):ncol(data)])))
    xv1 <- substr(xv, 1, 1)
    xv2 <- substr(xv, 2, 2)
    remove(xv)

    #convert to numeric
    xv1[xv1 == "A"] <- 1
    xv1[xv1 == "C"] <- 2
    xv1[xv1 == "G"] <- 3
    xv1[xv1 == "T"] <- 4
    xv1[xv1 == miss] <- 0
    xv2[xv2 == "A"] <- 1
    xv2[xv2 == "C"] <- 2
    xv2[xv2 == "G"] <- 3
    xv2[xv2 == "T"] <- 4
    xv2[xv2 == miss] <- 0

    #bind back into matrices
    xv1 <- matrix(xv1, ncol(data) - ecs, nrow(data), T)
    xv2 <- matrix(xv2, ncol(data) - ecs, nrow(data), T)

    #create output matrix and bind the data to it, structure format
    if(output == 3 | output == 9){
      outm <- matrix(NA, 2*(ncol(data) - ecs), nrow(data))

      #fill
      outm[seq(1,nrow(outm),2),] <- xv1
      outm[seq(2,nrow(outm),2),] <- xv2

      #add sample names
      snames <- character(nrow(outm))
      snames[seq(1,nrow(outm),2)] <- colnames(data[,(ecs+1):ncol(data)])
      snames[seq(2,nrow(outm),2)] <- colnames(data[,(ecs+1):ncol(data)])
      if(output == 3){ #bind sample names in and return for structure
        out <- cbind(ind = snames, as.data.frame(outm, stringsAsFactors = F))
      }
      else{ #add a bunch of filler columns (4 if pop info is provided, 5 if it isn't) for faststructure
        if(length(pop) > 1 & is.list(pop) == TRUE){
          out <- cbind(ind = snames,
                       matrix("filler", nrow(outm), 4),
                       as.data.frame(outm, stringsAsFactors = F))
          colnames(out)[2:5] <- paste0("filler", 1:4)
        }
        else{
          out <- cbind(ind = snames,
                       matrix("filler", nrow(outm), 5),
                       as.data.frame(outm, stringsAsFactors = F))
          colnames(out)[2:6] <- paste0("filler", 1:5)
        }
      }

      #add pop info if possible
      if(is.list(pop)){
        #add pop info if provided
        pop[[2]] <- pop[[2]]*2
        pl <- numeric(sum(pop[[2]]))
        pl[1:pop[[2]][1]] <- 1
        tracker <- pop[[2]][1] + 1
        for(i in 2:length(pop[[2]])){
          pl[tracker:(sum(pop[[2]][1:i]))] <- i
          tracker <- sum(pop[[2]][1:i]) + 1
        }
        out <- cbind(out[,1], pop = pl, out[,2:ncol(out)])
        colnames(out)[1] <- "ind"
      }
    }

    #create output matrix and bind to it, RAFM format.
    else{
      outm <- matrix(NA, ncol(data) - ecs, nrow(data)*2)

      #fill
      outm[,seq(1,ncol(outm),2)] <- xv1
      outm[,seq(2,ncol(outm),2)] <- xv2

      #replance missings with NA
      outm[outm == 0] <- NA

      #add column names
      colnames(outm) <- paste0("locus", sort(rep(1:nrow(data),2)), rep(c(".1", ".2"), ncol(outm)/2))

      #add subpop numbers, if given
      if(is.list(pop)){
        pl <- numeric(sum(pop[[2]]))
        pl[1:pop[[2]][1]] <- 1
        tracker <- pop[[2]][1] + 1
        for(i in 2:length(pop[[2]])){
          pl[tracker:(sum(pop[[2]][1:i]))] <- i
          tracker <- sum(pop[[2]][1:i]) + 1
        }
        out <- cbind(subpop = pl, as.data.frame(outm, stringsAsFactors = F))
      }
    }


    rdata <- out
  }

  #presence/absence format
  if(output == 7){
    x <- data[,(ecs+1):ncol(data)] #get just data
    if(input_form == "NN"){
      asize <- 1
    }
    else if (input_form == "0000" | input_form == "msat_2"){
      asize <- 2
    }
    else if (input_form == "msat_3"){
      asize <- 3
    }

    if(is.null(lnames)){
      lnames <- 1:nrow(x)
    }
    else{
      lnames <- gsub("_", "", lnames)
    }

    xv <- as.vector(t(x))

    #function to produce vectorized presence absence (as much as possible, not vectorized for NAs)
    pa_alleles <- function(xv, snp_form, mDat, nsamp, nloci){
      #get all possible genotypes
      gs <- unique(xv)

      #which genotype is the missing data?
      mpos <- which(gs == mDat)

      #what are the possible alleles at all loci?
      as <- unique(c(substr(gs,1,snp_form/2), substr(gs, (snp_form/2 + 1), snp_form*2)))
      as <- sort(as[as != substr(mDat, 1, snp_form/2)]) #that aren't N?

      #make the table
      cat("Creating presence/absence table...\n")

      #convert genotype vector to allele vector
      xva1 <- substr(xv, 1, snp_form/2)
      xva2 <- substr(xv, (snp_form/2+1),snp_form)

      #initialize
      amat <- matrix(0, nsamp, length(as)*nloci)
      colnames(amat) <- paste0(sort(rep(lnames,length(as))), "_", as) #initialize all of the columns, with locus number followed by allele. Will remove anything empty after assignment.

      #fill in
      for(i in 1:length(as)){
        pr1 <- grep(as[i], xva1)#unique rows which have the allele in either position one or position two.
        pr2 <- grep(as[i], xva2)
        amat[,grep(paste0("_", as[i]),colnames(amat))][pr1] <- amat[,grep(paste0("_", as[i]),colnames(amat))][pr1] + 1 #set the allele as present in the correct rows. This works because we look only at the G amat columns first, then put set only the individual IDs with a G to one. I think.
        amat[,grep(paste0("_", as[i]),colnames(amat))][pr2] <- amat[,grep(paste0("_", as[i]),colnames(amat))][pr2] + 1
      }

      #remove alleles not seen at loci.
      amat <- amat[,which(colSums(amat) != 0)]

      #get allele counts per loci:
      cat("Filling in missing data with NAs.\n")

      ###########
      #fill in missing data with NAs.
      #easy to vectorize if SNP data
      if(input_form == "NN" | input_form == "0000"){
        xmc <- which(xv == paste0(miss, miss)) #which samples had missing data?
        adj <- floor(xmc / nsamp) #how many loci over do I need to adjust xmc, since in amat each locus occupies two columns?
        adj[xmc%%nsamp == 0] <- adj[xmc%%nsamp == 0] - 1 #shift over anything that got screwed up by being in the last sample
        xmc <- xmc + (nsamp*adj) #adjust xmc for extra columns.
        if(any(amat[xmc] != 0) | any(amat[xmc + nsamp] != 0)){
          stop("Missing data values were not properly identified for replacement with NAs. This usually happens when SNP data is not completely bi-allelic. Try filtering out non-biallelic and non-polymorphic SNPs using filter_snps.\n")
        }
        amat[xmc] <- NA #make the first allele NA
        amat[xmc + nsamp] <- NA #make the second allele (another column over) NA.
      }

      #otherwise need to loop through loci and fill
      else{
        #get allele counts per locus
        labs <- gsub("_\\w+", "", colnames(amat))
        ll <- unique(labs)
        acounts <- rbind(label = ll, count = sapply(ll, function(x)sum(labs==x)))
        acounts <- as.numeric(acounts[2,])

        #intialize stuff, do first loci, then do the rest of the loci.
        spos <- 1 #starting column for loci
        epos <- acounts[1] #ending column for loci
        amat[, spos:epos][rowSums(amat[, spos:epos]) == 0,] <- NA
        for (i in 2:nloci){
          spos <- sum(acounts[1:(i - 1)]) + 1
          epos <- spos + acounts[i] - 1
          amat[, spos:epos][rowSums(amat[, spos:epos]) == 0,] <- NA
        }
      }
      return(amat)
    }

    amat <- pa_alleles(xv, asize*2, miss, ncol(x), nrow(x))

    if(interp_miss){
      #average number observed in columns
      cat("Interpolating missing data...\n")
      afs <- colMeans(amat, TRUE)
      temp <- which(is.na(amat))/nrow(amat)
      fill_element <- floor(temp) + 1 #get the column for each missing data point
      fill_element[which(temp %% 1 == 0)] <- fill_element[which(temp %% 1 == 0)] - 1 #correct for anything in the last row
      amat[which(is.na(amat))] <- afs[fill_element] #fill with the appropriate allele frequency.
    }
    else{cat("Finished. Warning: Missing data counts are also stored!\n")}
    amat <- cbind(samp = as.character(colnames(data)[(ecs+1):ncol(data)]), as.data.frame(amat, stringsAsFactors = F))
    rdata <- amat
  }

  ############################################################################
  #return the final product, printing an outfile if requested.
  if(outfile != FALSE){
    cat("Writing output file...\n")
    #for genepop
    if(output == 2){
      cat("\tPreparing genepop file...\n")
      #get list of snps
      llist <- paste0("SNP", "_", 1:ncol(rdata), ",")
      llist[length(llist)] <- paste0("SNP_", ncol(rdata))

      #write output
      cat(paste0(unlist(strsplit(outfile, split =  "\\."))[1], "_genepop\n"), file = outfile)
      cat(llist, "\nPOP\n", file = outfile, append = T)

      #write the tables, splitting by pop if requested:
      if(is.list(pop)){
        cat("\tWriting genepop file seperated by populations...\t")
        #need to sort by pops. Fist loop does this:
        j <- 1
        rdata$pop <- NA
        for (i in 1:(length(pop[[1]]) - 1)){
          rdata[j:(j+pop[[2]][i] - 1),]$pop <- pop[[1]][i]
          j <- j + pop[[2]][i]
        }
        rdata[j:nrow(rdata),]$pop <- pop[[1]][length(pop[[1]])]

        #sort and remove pop column
        rdata$rnames <- rownames(rdata)
        rdata <- dplyr::arrange(rdata, pop)
        rownames(rdata) <- rdata$rnames
        rdata <- rdata[,-which(colnames(rdata) %in% c("pop", "rnames"))]

        #second loop prints results.
        j <- 1
        pop[[2]] <- pop[[2]][order(pop[[1]])]
        pop[[1]] <- sort(pop[[1]])
        for (i in 1:(length(pop[[1]]) - 1)){
          cat(pop[[1]][i], "\t")
          data.table::fwrite(rdata[j:(j+pop[[2]][i] - 1),], outfile, quote = F, sep = "\t", col.names = F, row.names = T, append = T)
          cat("POP\n", file = outfile, append = T)
          j <- j + pop[[2]][i]
        }
        data.table::fwrite(rdata[j:nrow(rdata),], outfile, quote = F, sep = "\t", col.names = F, row.names = T, append = T)
        cat(pop[[1]][length(pop[[1]])], "\t Done.\n")
      }
      else{
        data.table::fwrite(rdata, outfile, quote = F, sep = "\t", col.names = F, row.names = T, append = T)
      }
    }
    else if(output == 1){
      #write the raw output
      data.table:fwrite(rdata, outfile, quote = FALSE, col.names = T, sep = "\t", row.names = F)
      #write a bayescan object if pop list was provided.
      if(is.list(pop)){
        outfile <- paste0(outfile, ".bayes")
        #write the header
        cat("[loci]=", nrow(data), "\n\n[populations]=", length(pop[[1]]), "\n\n", file = outfile, sep = "")

        #write the data for each population.
        for(i in 1:length(pop[[1]])){
          cat("[pop]=", i, "\n", file = outfile, append = T, sep = "") #write header
          wdat <- cbind(snp = 1:nrow(rdata[rdata$pop == pop[[1]][i],]),
                        rdata[rdata$pop == pop[[1]][i], which(colnames(rdata) %in% c("n_total", "n_alleles"))],
                        alleles = paste0(rdata[rdata$pop == pop[[1]][i],]$ni1, " ", rdata[rdata$pop == pop[[1]][i],]$ni2, " "))
          data.table::fwrite(wdat,
                      outfile, col.names = F, row.names = F, quote = F, sep = "\t",
                      append = T) #write the data for this population.
          cat("\n", file = outfile, append = T) #add a line break
        }
      }
    }
    else{
      data.table::fwrite(rdata, outfile, quote = FALSE, col.names = T, sep = "\t", row.names = F)
    }
  }
  return(rdata)
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
