#' Generate a 1-2d site frequency spectrum from a snpRdata object.
#'
#' Generates a 1 or 2 dimensional site frequency spectrum from a dadi input file
#' using the projection methods and folding methods of Marth et al (2004) and
#' Gutenkunst et al (2009). This code is essentially an R re-implementation of
#' the SFS construction methods implemented in the program \emph{dadi} (see
#' Gutenkunst et al (2009)).
#'
#' Site frequency spectrums are constructed using the projection methods
#' detailed in Marth et al (2004) and the 2 dimensional expansion in Gutenkunst
#' et al (2009). Folding methods are also taken from Gutenkunst et al (2009).
#' Either 1 or 2d SFSs can be constructed by providing a vector of population
#' names and projection sizes.
#'
#' Note that ref and anc columns are suggested in the SNP metadata, containing
#' the derived and ancestral character states, respectively. These should
#' contain three characters each: two flanking bases and the SNP. For example,
#' for an A/C SNP flanked by a G and a T, "GCT" and "GAT" would be expected.
#' Note that if these character states are not known, the minor and major
#' alleles will be substituted. Unfolded spectra will be misleading in this
#' case.
#'
#' @param x snpRdata object. The SNP metadata must contain "ref" and "anc" data.
#' @param facet character, default NULL. Name of the sample metadata column
#'   which specifies the source population of individuals. For now, allows only
#'   a single simple facet (one column).If NULL, runs the entire dataset.
#' @param pops character, default NULL. A vector of population names of up to
#'   length 2 containing the names of populations for which the an SFS is to be
#'   created. If NULL, runs the entire dataset.
#' @param projection numeric. A vector of sample sizes to project the SFS to, in
#'   \emph{number of gene copies}. Sizes too large will result in a SFS
#'   containing few or no SNPs. Must match the length of the provided pops
#'   vector.
#' @param fold logical, default FALSE. Determines if the SFS should be folded or
#'   left polarized. If FALSE, snp metadata columns named "ref" and "anc"
#'   containing the identity of the derived and ancestral alleles, respectively,
#'   should be present for polarization to be meaningful.
#' @param update_bib character or FALSE, default FALSE. If a file path to an
#'   existing .bib library or to a valid path for a new one, will update or
#'   create a .bib file including any new citations for methods used. Useful
#'   given that this function does not return a snpRdata object, so a
#'   \code{\link{citations}} cannot be used to fetch references.
#'
#' @references Gutenkunst et al (2009). Inferring the joint demographic history
#'   of multiple populations from multidimensional SNP frequency data.
#'   \emph{PLoS genetics}, 5(10), e1000695.
#' @references Marth et al (2004). The allele frequency spectrum in genome-wide
#'   human variation data reveals signals of differential demographic history in
#'   three large world populations. \emph{Genetics}, 166(1), 351-372.
#'
#' @export
#' @return A matrix or vector containing the site frequency spectrum with a
#'   "pops" attribute containing population IDs, such as c("POP1", "POP2"). For
#'   a 2d SFS, the first pop is the matrix columns and the second is the matrix
#'   rows.
#'
#' @author William Hemstrom
#'
#' @examples
#' 
#' \dontrun{
#' # add the needed ref and anc columns, using the major and minor alleles (will fold later)
#' dat <- calc_maf(stickSNPs)
#' # note, setting ref and anc is done by default if these columns don't exist!
#' snp.meta(dat)$ref <- paste0("A", get.snpR.stats(dat)$minor, "A") 
#' snp.meta(dat)$anc <- paste0("A", get.snpR.stats(dat)$major, "A")
#' 
#' # run for two populations
#' ## call calc_sfs()
#' sfs <- calc_sfs(dat, "pop", c("ASP", "CLF"), c(10,10))
#' ## plot
#' plot_sfs(x = sfs)
#' 
#' 
#' # run for the overall dataset
#' sfs <- calc_sfs(dat, projection = 30)
#' ## plot
#' plot_sfs(x = sfs)
#' 
#' # note that plot_sfs() will take a snpRdata object, calling calc_sfs()
#' plot_sfs(dat, projection = 30)
#' }
#' 
calc_sfs <- function(x, facet = NULL, pops = NULL, projection, fold = TRUE, 
                     update_bib = FALSE){
  #=============sanity checks=================
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  if(!.is.bi_allelic(x)){
    stop("This function is not yet supported for non-bi-allelic markers.\n")
  }
  
  msg <- character(0)
  
  if(is.null(pops) & !is.null(facet)){
    msg <- c(msg, "Pops must be provided if a facet is designated.\n")
  }
  if(is.null(facet) & !is.null(pops)){
    msg <- c(msg, "A facet must be provided if pops are designated.\n")
  }
  if(!is.null(facet)){
    if(length(facet) > 1){
      msg <- c(msg, "For now, only a single facet is allowed at a time.\n")
    }
    check_facet <- .check.snpR.facet.request(x, facet, "none", T)
    if(any(check_facet[[2]] != "sample")){
      msg <- c(msg, "For now, only sample level facets are allowed.\n")
    }
  }
  
  if(!is.null(pops)){
    if(length(pops) != length(projection)){
      msg <- c(msg, "A projection size must be provided for every requested population of samples.\n")
    }
  }
  
  
  if(any(!c("ref", "anc") %in% colnames(x@snp.meta))){
    om <- snp.meta(x)
    om$ref <- paste0("A", .get.snpR.stats(x)$minor, "A")
    om$anc <- paste0("A", .get.snpR.stats(x)$major, "A")
    .suppress_specific_warning(snp.meta(x) <- om, "duplicated")

    if(fold == FALSE){
      warning("Without ancestral and derived character states, unfolded spectra will be misleading.\n")
    }
  }
  else{
    if(!all(sapply(c(x@snp.meta$ref, x@snp.meta$anc), nchar) == 3)){
      msg <- c(msg, "All ref and anc entries must be exactly three characters long. See documentation for details.\n")
    }
  }

  if(length(msg) > 0){
    stop(msg)
  }
  #==============run=========================
  # subset
  if(!is.null(pops)){
    suppressWarnings(x <- .subset_snpR_data(x, facets = facet, subfacets = pops))
  }

  # get dadi formatted data
  y <- format_snps(x, output = "dadi", facets = facet)

  # get sfs
  if(is.null(pops)){pops <- ".base"}
  sfs <- make_SFS(y, pops, projection, fold, update_bib)
  
  return(sfs)
}





#' Generate a 1-2d site frequency spectrum from a dadi input file.
#'
#' Generates a 1 or 2 dimensional site frequency spectrum from a dadi input file
#' using the projection methods and folding methods of Marth et al (2004) and
#' Gutenkunst et al (2009). This code is essentially an R re-implementation of
#' the SFS construction methods implemented in the program \emph{dadi} (see
#' Gutenkunst et al (2009)).
#'
#' Site frequency spectrums are constructed using the projection methods
#' detailed in Marth et al (2004) and the 2 dimensional expansion in Gutenkunst
#' et al (2009). Folding methods are also taken from Gutenkunst et al (2009).
#' Either 1 or 2d SFSs can be constructed by providing a vector of population
#' names and projection sizes.
#'
#' @param x character or data.frame. Either a path to a dadi formatted input
#'   file or a data.frame containing previously imported dadi formatted data.
#' @param pops character. A vector of population names of up to length 2
#'   containing the names of populations for which the an SFS is to be created.
#' @param projection numeric. A vector of sample sizes to project the SFS to, in
#'   \emph{number of gene copies}. Sizes too large will result in a SFS
#'   containing few or no SNPs.
#' @param fold logical, default FALSE. Determines if the SFS should be folded or
#'   left polarized.
#' @param update_bib character or FALSE, default FALSE. If a file path to an
#'   existing .bib library or to a valid path for a new one, will update or
#'   create a .bib file including any new citations for methods used. Useful
#'   given that this function does not return a snpRdata object, so a
#'   \code{\link{citations}} cannot be used to fetch references.
#'
#' @references Gutenkunst et al (2009). Inferring the joint demographic history
#'   of multiple populations from multidimensional SNP frequency data.
#'   \emph{PLoS genetics}, 5(10), e1000695.
#' @references Marth et al (2004). The allele frequency spectrum in genome-wide
#'   human variation data reveals signals of differential demographic history in
#'   three large world populations. \emph{Genetics}, 166(1), 351-372.
#'
#' @export
#' @return A matrix or vector containing the site frequency spectrum with a
#'   "pops" attribute containing population IDs, such as c("POP1", "POP2"). For
#'   a 2d SFS, the first pop is the matrix columns and the second is the matrix
#'   rows.
#'   
make_SFS <- function(x, pops, projection, fold = FALSE, update_bib = FALSE){
  #================sanity checks=========
  msg <- character()
  if(!is.data.frame(x)){
    if(!is.character(x)){
      msg <- c(msg, "x must either be either a dadi formatted data file or a data.frame derived from such.\n")
    }
    else if(!file.exists(x)){
      msg <- c(msg, "x must either be either a dadi formatted data file or a data.frame derived from such.\n")
    }
    else{
      x <- utils::read.table(x, header = T, stringsAsFactors = F)
    }
  }
  if(length(msg) > 0){
    stop(msg)
  }


  if(any(sapply(pops, function(y) sum(colnames(x) == y) != 2))){
    msg <- c(msg, "Each pop must match two column names in x.\n")
  }
  if(length(pops) > 2){
    msg <- c(msg, "Only one or two dimensional SFSs can be created.\n")
  }

  if(length(msg) > 0){
    stop(msg)
  }
  #================subfunctions==========
  # function to fold an sfs
  fold_sfs <- function(sfs){

    # fold a 2 dimensional sfs
    if(length(dim(sfs)) == 2){
      rev.matrix <- function(x){
        return(x[nrow(x):1, ncol(x):1])
      }

      # figure out the part that gets folded in
      sample_sizes <- matrix(0:(nrow(sfs) - 1), nrow(sfs), ncol(sfs))
      sample_sizes <- sample_sizes +
        t(matrix(0:(ncol(sfs) - 1), ncol(sfs), nrow(sfs)))

      total.samples <- (ncol(sfs) + nrow(sfs) - 2)
      folded.in <- ifelse(sample_sizes > floor(total.samples/2), T, F)

      # reverse the sfs
      reversed <- sfs
      reversed[folded.in == F] <- 0
      reversed <- rev.matrix(reversed)

      # get the folded sfs
      folded <- sfs + reversed

      # deal with ambiguous areas
      ambig <- ifelse(sample_sizes == total.samples/2, T, F)
      ambig[ambig == T] <- sfs[ambig == T]
      rev.ambig <- rev.matrix(ambig)
      folded <- folded + (-.5*ambig + .5*rev.ambig)

      folded[folded.in == T] <- NA
    }

    # fold a 1 dimensional sfs--life is easy
    else{
      rev <- rev(sfs)[1:floor(length(sfs)/2)]
      rev <- c(rev, rep(0, length(sfs) - length(rev)))
      folded <- sfs + rev
      folded[ceiling(length(sfs)/2):length(sfs)] <- NA
    }
    return(folded)
  }

  # make a projected sfs from a counts array
  make_proj_sfs <- function(counts, projection, fold = F){

    # initialize
    if(length(projection) == 2){
      sfs <- matrix(0, projection[1] + 1, projection[2] + 1)
    }
    else{
      sfs <- numeric(projection + 1)
      counts[,,2] <- 1
    }

    # project each snp
    for(i in 1:nrow(counts[,,1])){

      # if the snp isn't sequenced in any pop, skip
      if(any(counts[i,1,] == 0)){
        next
      }

      # if we aren't folding and this snp isn't polarizable (usually meaning we had NNN for the anc condition),
      # skip it
      if(!fold & !counts[i,3,1]){
        next
      }

      if(length(projection) == 2){
        pop_counts <- list(counts[i,,1], counts[i,,2])
      }
      else{
        pop_counts <- list(counts[i,,1])
      }
      contrib <- find_contribs(pop_counts, projection)
      sfs <- sfs + contrib
    }

    return(sfs)
  }

  # function to get a pop_counts array for projections. Result has three columns,
  # containing the three bits of info needed for a pop_counts list. The third
  # index lists the pop
  # input is a dadi formatted data.frame
  get_counts <- function(x, pops){
    # figure out which allele is the outgroup/derived in each row and figure out polarization status
    anc.als <- substr(x$anc, 2, 2)
    ref.als <- substr(x$ref, 2, 2)
    polarized <- rep(T, nrow(x))
    polarized[anc.als == "N"] <- F
    polarized[anc.als != x$Allele1 & anc.als != x$Allele2] <- F
    polarized[ref.als != x$Allele1 & ref.als != x$Allele2] <- F
    anc.als[anc.als == "N"] <- x$Allele1[anc.als == "N"]
    which.derived <- rep(2, nrow(x))
    which.derived[anc.als == x$Allele2] <- 1

    # initialize output
    out <- array(0, dim = c(nrow(x), 3, 2))

    # write to output
    for(i in 1:length(pops)){
      x <- data.table::as.data.table(x)
      dat.cols <- which(pops[i] == colnames(x))
      tdat <- x[,dat.cols, with = FALSE]
      tdat <- as.matrix(tdat)
      out[,1,i] <- rowSums(tdat) # total count
      t.index <- 2 * (1:nrow(tdat))
      t.index[which.derived == 1] <- t.index[which.derived == 1] - 1
      out[,2,i] <- t(tdat)[t.index] # derived allele count
      out[,3,i] <- polarized
    }

    return(out)
  }


  # function to get an sfs contribution for a single snp across 2 populations
  # takes a list of length 2, where the each element is a numeric vector containing:
  # 1) the number of sequenced alleles at the site
  # 2) the number of derived alleles at the site
  # 3) the polarization status
  #
  # also needs a vector of projection sizes for the populations
  find_contribs <- function(pop_counts, projection){
    pcontribs <- vector("list", length(pop_counts))

    for(i in 1:length(pop_counts)){
      # pull info
      hits <- pop_counts[[i]][2]
      n <- pop_counts[[i]][1]
      m <- projection[i]


      # if there are less snps sequenced here than our projection, it gets zero
      # for all contributions
      if(n < m){
        contrib <- rep(0, m + 1)
      }

      else{
        # make the hits vector
        proj_hits <- 0:m

        # do the binomials
        contrib <- choose(m, proj_hits)
        contrib <- contrib*choose(n - m , hits - proj_hits)
        contrib <- contrib/choose(n, hits)
      }
      pcontribs[[i]] <- contrib
    }

    return(Reduce(outer, pcontribs))
  }

  #=========================run the functions================
  counts <- get_counts(x, pops)
  sfs <- make_proj_sfs(counts, projection, fold)
  if(fold){
    sfs <- fold_sfs(sfs)
    if(length(projection) == 1){
      sfs <- sfs[1:floor(projection/2)]
    }
    attr(sfs, which = "folded") <- TRUE
  }
  else{
    attr(sfs, which = "folded") <- FALSE
  }

  # add a pops attribute
  attr(sfs, which = "pop") <- pops

  # return
  masked <- sfs
  if(is.matrix(masked)){
    masked[1,1] <- 0
    if(!fold){
      masked[nrow(masked), ncol(masked)] <- 0
    }
  }
  else{
    masked[1] <- 0
    if(!fold){
      masked[length(masked)] <- 0
    }
  }
  if(sum(masked, na.rm = T) == 0){
    stop("No segrgating sites remain after projection. Try decreasing projection sizes!\n")
  }
  
  cat("SFS completed with", sum(masked, na.rm = T), "segrgating sites.\n")
  
  .yell_citation("Gutenkunst2009", "SFS", "Used to project and possibly fold SFS data. Code is a direct R re-implementation.", outbib = update_bib)
  
  return(sfs)
}

#' Calculate the directionality index from a 2d site frequency spectrum.
#'
#' Calculates the directionality index based on a 2d SFS according to Peter and
#' Slatkin (2013). Input spectra can be created using the \code{\link{calc_sfs}}
#' function, using a provided snpRdata object, or passed from other programs.
#' \emph{Spectra must be not be folded}.
#'
#' Essentially, the directionality index measures the difference in derived
#' allele frequency between two populations to determine the directionality of
#' population spread between the two. Since the "destination" population is
#' sourced from but experienced more genetic drift than the "source" population,
#' it should have relatively more high-frequency derived alleles \emph{after the
#' removal of fixed ancestral alleles}. See Peter and Slatkin (2013) for
#' details.
#'
#' @param x snpRdata object or matrix, default NULL. A snpRdata from which to
#'   calculate a SFS. Alternatively, a 2d site frequency spectra stored in a
#'   matrix, with an additional "pop" or "pops" attribute containing population
#'   IDs, such as c("POP1", "POP2"), where the first pop corresponds to matrix
#'   columns and the second to matrix rows. These objects can be produced from a
#'   dadi input file using \code{\link{make_SFS}}. Note that if x is a snpRdata
#'   object, snp metadata columns named "ref" and "anc" containing the identity
#'   of the derived and ancestral alleles, respectively, must be present.
#' @param facet character, default NULL. Passed to \code{\link{calc_sfs}} -- see
#'   documentation there for details. Ignored if a sfs is provided.
#' @param pops character, default NULL. Passed to \code{\link{calc_sfs}} -- see
#'   documentation there for details. Ignored if a sfs is provided.
#' @param projection numeric, default NULL. Passed to \code{\link{calc_sfs}} --
#'   see documentation there for details. Ignored if a sfs is provided.
#' @param update_bib character or FALSE, default FALSE. If a file path to an
#'   existing .bib library or to a valid path for a new one, will update or
#'   create a .bib file including any new citations for methods used. Useful
#'   given that this function does not return a snpRdata object, so a
#'   \code{\link{citations}} cannot be used to fetch references.
#'
#' @export
#' @references Peter, B. M., & Slatkin, M. (2013). Detecting range expansions
#'   from genetic data. \emph{Evolution}, 67(11), 3274-3289.
#'
#' @return A numeric value giving the directionality with a "direction"
#'   attribute designating the direction between the two populations.
#'   
#' @examples
#' \dontrun{
#' # directionality can be calculated without first calculating a SFS
#' calc_directionality(stickSNPs, facet ="pop", pops = c("ASP", "PAL"), projection = c(10, 10))
#'
#' # an existing SFS can also be fed in. This may be handy if you get a SFS from elsewhere.
#' sfs <- calc_sfs(stickSNPs, "pop", c("ASP", "PAL"), c(10, 10), fold = FALSE)
#' calc_directionality(sfs)
#' }
calc_directionality <- function(x, facet = NULL, pops = NULL, projection = NULL, update_bib = FALSE){
  #==========sanity checks=============
  msg <- character(0)
  
  if(is.snpRdata(x)){
    
    if(!.is.bi_allelic(x)){
      stop("This function is not yet supported for non-bi-allelic markers.\n")
    }
    
    if(any(c(is.null(facet), is.null(pops), is.null(projection)))){
        msg <- c(msg, "facet, pops, and projection arguments must all be provided if an sfs is not.")
    }
    
    if(!is.null(pops)){
      if(length(pops) != 2){
        msg <- c(msg, "Exactly two pops must be listed.\n")
      }
    }
    
    x <- calc_sfs(x, facet, pops, projection, fold = F)
  }
  
  else{
    msg <- c(msg, .sanity_check_sfs(x, 2))
    if(any(is.na(x))){
      msg <- c(msg, "NAs found in provided sfs. This most likely means that the SFS is folded, which is not permitted for directionality calculation.\n")
    }
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  #===========calc direcitonality===========
  
  # flip everthing, since i should be pop 1 and j should be pop 2, and the
  # matrix is set up so that rows are pop 1 and columns are pop 2. Transposing fixes it.
  pops <- attr(x, "pop")
  x <- t(x)

  # normalize, exluding fixed sites
  x <- x[-1,-1] # remove fixed sites
  x <- x/sum(x)

  # get all of the allele frequencies in each cell
  freqs.i <- (1:(nrow(x)))/(nrow(x))
  freqs.j <- (1:(ncol(x)))/(ncol(x))

  # do the math, this is equ 1b (i - j times sfs)
  directionality <- sum(outer(freqs.i, freqs.j, "-") * x)
  if(directionality > 0){
    attr(directionality, "direction") <- paste0(pops[1], "->", pops[2])
  }
  else{
    attr(directionality, "direction") <- paste0(pops[1], "<-", pops[2])

  }
  
  .yell_citation("Peter2013", "Direcitonality", "Directionality index (psi)", update_bib)

  return(directionality)
}


#' Estimate the origin point of a population expansion using directionality.
#'
#' Calculates the origin of expansion from directionality indices from pairwise
#' population comparisons according to Peter and Slatkin (2013).
#'
#' Essentially, the directionality index measures the difference in derived
#' allele frequency between two populations to determine the directionality of
#' population spread between the two. Since the "destination" population is
#' sourced from but experienced more genetic drift than the "source" population,
#' it should have relatively more high-frequency derived alleles \emph{after the
#' removal of fixed ancestral alleles}. See Peter and Slatkin (2013) for
#' details.
#' 
#' This metric, calculated between multiple populations, can be used to estimate
#' the origin point of a population expansion using a Time Difference Of Arrival
#' approach by scaling allelic differences to geometric space. See Peter and 
#' Slatkin (2013) for details. Distances are calculated using the 
#' \code{\link[geosphere]{distGeo}} function assuming a WGS84 ellipsoid, and so
#' provided coordinates must be in that format. The sampling location for each
#' subfacet is derived using the \code{\link[geosphere]{geomean}} of each.
#' 
#' The TODA method requires an optimization procedure to estimate the point-of-
#' origin. This is done using the \code{\link[stats]{optim}} function with
#' the default parameters for maximization.
#'
#' @param x snpRdata object. A snpRdata from which to
#'   calculate SFS, directionality indices, and an expansion point-of-origin.
#'   SNP metadata columns named "ref" and "anc" containing the identity of the
#'   derived and ancestral alleles, respectively, must be present, as must
#'   columns named 'x' and 'y' containing WGS84 elipsoid scaled sampling
#'   locations for each sample in the sample metadata.
#' @param facet character, default NULL. The sample metadata facet by which
#'   to group populations and calculate the expansion origin. Only one facet
#'   currently allowed, although it can be complex (e.g. 'fam.pop').
#' @param boots numeric, default 1000. The number of bootstraps used to determine
#'   the variance of the directionality index for each pairwise comparison.
#' @param projection numeric, default NULL. The number of \emph{gene copies}
#'   to project to for each facet level. Should be a named numeric vector containing
#'   an entry for each facet level.
#' @param boot_par numeric or FALSE, default FALSE. If a number, bootstraps will
#'  be processed in parallel using the supplied number of cores.
#' @param verbose Logical, default FALSE. If TRUE, some progress updates will be
#'  reported.
#' @param update_bib character or FALSE, default FALSE. If a file path to an
#'   existing .bib library or to a valid path for a new one, will update or
#'   create a .bib file including any new citations for methods used. Useful
#'   given that this function does not return a snpRdata object, so a
#'   \code{\link{citations}} cannot be used to fetch references.
#'
#' @export
#' @references Peter, B. M., & Slatkin, M. (2013). Detecting range expansions
#'   from genetic data. \emph{Evolution}, 67(11), 3274-3289.
#'
#' @return A named list containing: 
#' * opt: A vector with the spatial/genetic distance linking coefficient 'v' as 
#'   well as the 'x' and 'y' coordinates of the estimated origin of the range
#'   expansion.
#' * 'pairwise_directionality': A data.frame containing the pairwise
#'   directionality estimates, coordinates, and the directionality variance for
#'   each pair of populations.
#'   
#' @author William Hemstrom
#'   
#' @examples
#' \dontrun{
#' # Bootstrapping is slow, so not run.
#' # set ref and anc--ideally use an outgroup for this
#' dat <- calc_maf(stickSNPs)
#' snp.meta(dat)$ref <- paste0("A", get.snpR.stats(dat)$minor, "A") 
#' snp.meta(dat)$anc <- paste0("A", get.snpR.stats(dat)$major, "A")
#' 
#' # setup x and y coords
#' long_lat <- data.frame(SMR = c(44.365931, -121.140420), 
#'                        CLF = c(44.267718, -121.255805), 
#'                        OPL = c(44.485958, -121.298360), 
#'                        ASP = c(43.891693, -121.448360), 
#'                        UPD = c(43.891755, -121.451600), 
#'                        PAL = c(43.714114, -121.272797))
#' long_lat <- t(long_lat)
#' long_lat <- long_lat[match(sample.meta(dat)$pop, rownames(long_lat)),]
#' colnames(long_lat) <- c("y", "x")
#' sample.meta(dat) <- cbind(sample.meta(dat), long_lat)
#' 
#' projection <- summarize_facets(dat, facet)[[facet]]
#' projection <- floor(projection*.8)
#' 
#' # run the calculation
#' calc_origin_of_expansion(dat, "pop", boots = 100, projection = projection, 
#'                          boot_par = 6, verbose = TRUE)
#' }
calc_origin_of_expansion <- function(x, facet, boots = 1000, projection = NULL,
                                     boot_par = FALSE,
                                     verbose = FALSE,
                                     update_bib = FALSE){
  y <- NULL

  #============sanity checks============
  .check.installed("geosphere")
  
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.n")
  }
  
  if(!.is.bi_allelic(x)){
    stop("This function is not yet supported for non-bi-allelic markers.\n")
  }
  
  # check facet and proj
  ofacet <- facet
  facet <- .check.snpR.facet.request(x, facet)
  if(length(facet) > 1){
    stop("calc_origin_of_expansion currently allows for only one facet at a time.\n")
  }
  
  if(facet == ".base"){
    stop("A sample facet must be provided.\n")
  }

  # resort complex if needed
  if(facet != ofacet){
    sf1 <- unlist(.split.facet(ofacet))
    sf2 <- unlist(.split.facet(facet))
    
    resort <- match(sf1, sf2)
    
    new.names <- .split.facet(names(projection))
    new.names <- unlist(lapply(new.names, function(y) paste(y[resort], collapse = ".")))
    names(projection) <- new.names
  }
  
  fl <- summarize_facets(x, facet)[[facet]]
  if(length(fl) < 3){
    stop("At least three populations/facet levels must be present in provided facet to calculate origin-of-expansion.\n")
  }
  missing_proj <- names(fl)[which(!names(fl) %in% names(projection))]
  if(length(missing_proj) > 0){
    stop("Some facet levels are missing in the projection vector. Missing levels:\n\t", paste0(missing_proj, collapse = ", "), "\n")
  }
  
  large_proj <- which(projection > 2*fl[match(names(projection), names(fl))])
  if(length(large_proj) > 0){
    stop("All projections must be smaller than or equal to 2N, where N is the sample size for the respective facet level.\nBad projections:\n\t",
         paste0(names(projection)[large_proj], collapse = ", "), "\n")
  }
  
  if(any(!c("ref", "anc") %in% colnames(x@snp.meta))){
    om <- snp.meta(x)
    om$ref <- paste0("A", .get.snpR.stats(x)$minor, "A")
    om$anc <- paste0("A", .get.snpR.stats(x)$major, "A")
    .suppress_specific_warning(snp.meta(x) <- om, "duplicated")
    
    warning("Without ancestral and derived character states, results will be misleading.\n")
  }
  else{
    if(!all(sapply(c(x@snp.meta$ref, x@snp.meta$anc), nchar) == 3)){
      stop("All ref and anc entries must be exactly three characters long. See documentation for details.\n")
    }
  }
  
  if(any(!c("x", "y") %in% colnames(sample.meta(x)))){
    stop("'x' and 'y' columns containing coordinates must be present in sample metadata. These are assumed conform to a WGS84 elipsoid.\n")
  }
  
  
  #============functions================
  come_to_dadi <- function(ac, ref, anc, major, minor){
    ni1 <- reshape2::dcast(ac, facet + .snp.id ~ subfacet, value.var = c("ni1"))
    ni2 <- reshape2::dcast(ac, facet + .snp.id ~ subfacet, value.var = c("ni2"))

    rdata <- cbind(ref = ref, # since everything is sorted by .snp.id, this will match.
                   anc = anc,
                   Allele1 = major,
                   ni1[order(ni1$.snp.id), 3:ncol(ni1), drop = FALSE],
                   Allele2 = minor,
                   ni2[order(ni2$.snp.id), 3:ncol(ni2), drop = FALSE])
    rdata <- as.data.frame(rdata)
    colnames(rdata)[c(3,3 + length(3:ncol(ni1)) + 1)] <- c("Allele1", "Allele2")

    return(rdata)
  }
  
  # boot_sfs <- function(sfs, boots, pop){
  #   nsnps <- sum(sfs, na.rm = TRUE)
  #   boot_draws <- rmultinom(nsnps*boots, 1, sfs)
  #   
  #   boot_draws <- array(boot_draws, c(length(sfs), nsnps, boots))
  #   boot_draws <- apply(boot_draws, 3, 
  #                       function(y) matrix(rowSums(y), nrow(sfs), ncol(sfs)), 
  #                       simplify = FALSE)
  #   boot_draws <- lapply(boot_draws, function(y){
  #     attr(y, "pop") <- pop
  #     return(y)
  #   })
  #   
  #   return(boot_draws)
  # }
  
  # Equation 4 from Peter and Slatkin
  opt_eq <- function(v, x, y, dirs, var_dirs, xi, yi, xj, yj){
    # cat("=================\nx = ", x, "\ny = ", y, "\nv = ", v, "\n")
    
    internal <- (geosphere::distGeo(cbind(xi, yi), c(x, y)) -
                   geosphere::distGeo(cbind(xj, yj), c(x, y)))
    # internal <- sqrt(((xi - x)^2) + (yi - y)^2) -
    #   sqrt(((xj - x)^2) + (yj - y)^2)
    internal <- (1/v) * internal
    internal <- internal - dirs
    internal <- (1/var_dirs) * internal
    
    res <- sum(internal)
    return(res)
  }
  
  #============setup====================
  if(verbose){
    cat("Beginning setup...\n")
  }
  facet <- .check.snpR.facet.request(x, facet)
  
  x <- .add.facets.snpR.data(x, facet)
  x <- calc_maf(x, facet)
  
  # generate bootstraps
  ## maybe do this by drawing n snps randomly from the SFS using the value in each cell as a scaled binom prob?
  ## this is fast but produces very different results...
  
  ## this is the most robust, but ends up being very slow. The alternative works well but is much quicker
  dadi <- .boot_ac(x, boots, facet)
  dadi <- lapply(dadi, function(y) come_to_dadi(y,
                                               snp.meta(x)$ref,
                                               snp.meta(x)$anc,
                                               x@stats[which(x@stats$facet == ".base" & x@stats$subfacet == ".base"), "major"],
                                               x@stats[which(x@stats$facet == ".base" & x@stats$subfacet == ".base"), "minor"]))

  
  # set up tasks and get positions
  tasks <- .get.task.list(x, facet)
  tasks <- t(utils::combn(tasks[,2], 2))
  split_facet <- unlist(.split.facet(facet))
  long_lat <- as.data.table(sample.meta(x)[,c(split_facet, "x", "y")])
  geomean_func <- function(x, y){
    m <- as.matrix(cbind(x, y))
    return(as.data.frame(matrix(geosphere::geomean(m), nrow = 1)))
  }
  
  xy <- long_lat[, geomean_func(x, y), by = split_facet]
  colnames(xy) <- c(split_facet, "x", "y")
  
  
  #==============directionality for each pair of populations=================
  dir_list <- data.table(dirs = numeric(nrow(tasks)),
                         dir_vars = numeric(nrow(tasks)),
                         xi = numeric(nrow(tasks)),
                         yi = numeric(nrow(tasks)),
                         xj = numeric(nrow(tasks)),
                         yj = numeric(nrow(tasks)),
                         comparison = character(nrow(tasks)))
  
  if(!isFALSE(boot_par)){
    cl <- parallel::makePSOCKcluster(boot_par)
    doParallel::registerDoParallel(cl)
  }
  
  
  if(verbose){cat("Complete.\nBeginning pairwise directionality bootstrapping and calculation...\nTotal Pairs:", nrow(tasks), "\n")}
  for(i in 1:nrow(tasks)){
    if(verbose){cat(paste0("Working on pair: ", paste0(tasks[i,], collapse = "~"), " (task ", i, ")\n"))}
    
    ## this is faster but doesn't produce the same outputs...
    # .make_it_quiet(sfs <- calc_sfs(x, facet, 
    #                                pops = tasks[i,], 
    #                                projection = projection[match(tasks[i,], names(projection))],
    #                                fold = FALSE))
    # sfs_boots <- boot_sfs(sfs, boots, pop = attr(sfs, "pop"))
    # 
    # .make_it_quiet(dirs <- lapply(sfs_boots, function(y) calc_directionality(y)))
    
   
    
    # get the directionality for each bootstrap
    if(isFALSE(boot_par)){
      # this is for the slower way with booting earlier
      .make_it_quiet(dirs <- lapply(dadi, function(y) make_SFS(y,
                                                               pops = tasks[i,],
                                                               projection = projection[match(tasks[i,], names(projection))],
                                                               fold = FALSE)))
      
      .make_it_quiet(dirs <- lapply(dirs, function(y) calc_directionality(y)))
    }
    else{

      if(boot_par < boots){
        pboot <- split(1:boots, rep(1:boot_par, length.out = boots, each = ceiling(boots/boot_par)))
      }
      else{
        boot_par <- boots
        pboot <- split(1:boots, 1:boot_par, drop = F)
      }
      
      ntasks <- length(pboot)
      dirs <- foreach::foreach(q = 1:ntasks,
                               .packages = c("snpR", "data.table"),
                               .export = c(".make_it_quiet"),
                               .combine = c,.inorder = TRUE) %dopar% {
                                 .make_it_quiet(dirs <- lapply(dadi[pboot[[q]]], function(y) make_SFS(y,
                                                                                                      pops = tasks[i,],
                                                                                                      projection = projection[match(tasks[i,], names(projection))],
                                                                                                      fold = FALSE)))
                                 .make_it_quiet(res <- lapply(dirs, function(y) calc_directionality(y)))
                                 res
                               }
    }
    
    
    # get the real directionality
    .make_it_quiet(real_dir <- calc_directionality(x, facet = facet, pops = tasks[i,], projection = projection[match(tasks[i,], names(projection))]))
    
    
    # fill in
    dir_list$dirs[i] <- real_dir
    dir_list$dir_vars[i] <- stats::var(unlist(dirs))
    dir_list$xi[i] <- xy$x[match(tasks[i,1], .paste.by.facet(xy, split_facet))]
    dir_list$yi[i] <- xy$y[match(tasks[i,1], .paste.by.facet(xy, split_facet))]
    dir_list$xj[i] <- xy$x[match(tasks[i,2], .paste.by.facet(xy, split_facet))]
    dir_list$yj[i] <- xy$y[match(tasks[i,2], .paste.by.facet(xy, split_facet))]
    dir_list$comparison[i] <- paste0(tasks[i,1], "~", tasks[i,2])
  }
  
  if(!isFALSE(boot_par)){
    parallel::stopCluster(cl)
  }
  
  #=============optimize and return===================
  if(verbose){cat("Complete.\nBeginning optimization.\n")}
  start_v <- 1
  start_xy <- geosphere::geomean(dir_list[,c("xi", "yi")])
  
  
  opt <- stats::optim(par = c(v = start_v, x = as.numeric(start_xy[1,1]), y = as.numeric(start_xy[1,2])), 
                      fn = function(y){
                        opt_eq(v = y[1], x = y[2], y = y[3], 
                               dirs = dir_list$dirs, var_dirs = dir_list$dir_vars, 
                               xi = dir_list$xi, yi = dir_list$yi, 
                               xj = dir_list$xj, yj = dir_list$yj)
                      }, 
                      control = list(fnscale = -1))
  
  if(verbose){cat("Complete.\n")}
  
  .yell_citation("Peter2013", "Direcitonality", "Directionality index (psi)", update_bib)
  .yell_citation("Peter2013", "Origin of Expansion", "Origin of Expansion", update_bib)
  
  colnames(dir_list)[1:2] <- c("Directionality", "Variance")
  
  return(list(opt = opt$par, pairwise_directionality = dir_list))
}


.sanity_check_sfs <- function(x, allowed_dimensions = c(1, 2)){
  msg <- character(0)
  if("pops" %in% names(attributes(x))){
    attr(x, "pop") <- attr(x, "pops")
    attr(x, "pops") <- NULL
  }
  
  dims <- length(dim(x))
  if(dims == 0){dims <- 1}
  if(!dims %in% allowed_dimensions){
    msg <- c(msg, paste0("SFS dimensionality ", dims, " not allowed for this function. Allowed dimensionalities: ", paste0(allowed_dimensions, collapse = ", ")))
    return(msg)
  }
  
  if(is.matrix(x)){
    if("pop" %in% names(attributes(x))){
      if(length(attr(x, "pop")) == 2){
        sfs <- x
      }
      else{
        msg <- c(msg, "2D SFS pop attribute must be length 2.\n")
      }
    }
    else{
      msg <- c(msg, "2D SFS must have a pop attribute with length 2.\n")
    }
  }
  else if(is.numeric(x)){
    if("pop" %in% names(attributes(x))){
      if(length(attr(x, "pop")) == 1){
        sfs <- x
      }
      else{
        msg <- c(msg, "1D SFS pop attribute must be length 1.\n")
      }
    }
    else{
      msg <- c(msg, "1D SFS must have a pop attribute with length 1.\n")
    }
  }
  else{
    msg <- c(msg, "x must be either a snpRdata object, a matrix containing a 2D SFS, or a numeric vector containing a 1D sfs.\n")
  }
  
  return(msg)
}
