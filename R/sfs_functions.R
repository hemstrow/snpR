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
#' sfs <- calc_sfs(dat, "pop", c("ASP", "CLF"), c(30,30))
#' ## plot
#' plot_sfs(sfs = sfs)
#' 
#' 
#' # run for the overall dataset
#' sfs <- calc_sfs(dat, projection = 100)
#' ## plot
#' plot_sfs(sfs = sfs)
#' 
#' # note that plot_sfs() will take a snpRdata object, calling calc_sfs()
#' plot_sfs(dat, projection = 100)
#' }
#' 
calc_sfs <- function(x, facet = NULL, pops = NULL, projection, fold = TRUE, 
                     update_bib = FALSE){
  #=============sanity checks=================
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
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
    if(length(.split.facet(facet)[[1]]) > 1){
      msg <- c(msg, "For now, only facets referring to only one column of metadata are allowed.\n")
    }
  }
  
  if(!is.null(pops)){
    if(length(pops) != length(projection)){
      msg <- c(msg, "A projection size must be provided for every requested population of samples.\n")
    }
  }
  
  
  if(any(!c("ref", "anc") %in% colnames(x@snp.meta))){
    warning("ref and anc columns are suggested in snp metadata. See documentation for details. The major allele will be subsituted for the ancestral state.\n")
    x@snp.meta$ref <- paste0("A", .get.snpR.stats(x)$minor, "A")
    x@snp.meta$anc <- paste0("A", .get.snpR.stats(x)$major, "A")

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


  if(any(sapply(pops, function(y) length(grep(y, colnames(x)))) != 2)){
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
      dat.cols <- grep(pops[i], colnames(x))
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
#' @param x snpRdata object, default NULL. If provided, snpRdata from which to
#'   calculate a SFS. Ignored if a sfs is provided directly.
#' @param sfs numeric matrix. A 2d site frequency spectra stored in a matrix,
#'   with an additional "pops" attribute containing population IDs, such as
#'   c("POP1", "POP2"), where the first pop corresponds to matrix columns and
#'   the second to matrix rows. These objects can be produced from a dadi input
#'   file using \code{\link{make_SFS}}.
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
#' calc_directionality(stickSNPs, facet ="pop", pops = c("ASP", "PAL"), projection = c(20, 20))
#'
#' # an existing SFS can also be fed in. This may be handy if you get a SFS from elsewhere.
#' sfs <- calc_sfs(stickSNPs, "pop", c("ASP", "PAL"), c(20, 20), fold = FALSE)
#' calc_directionality(sfs = sfs)
#' }
calc_directionality <- function(x = NULL, sfs = NULL, facet = NULL, pops = NULL, projection = NULL, update_bib = FALSE){
  #==========sanity checks=============
  msg <- character(0)
  if(!is.null(x)){
    if(!is.null(sfs)){
      x <- NULL
    }
    else{
      if(!is.snpRdata(x)){
        msg <- c(msg, "x is not a snpRdata object.\n")
      }
      if(any(c(is.null(facet), is.null(pops), is.null(projection)))){
        msg <- c(msg, "facet, pops, and projection arguments must all be provided if an sfs is not.")
      }
      if(!is.null(pops)){
        if(length(pops) != 2){
          msg <- c(msg, "Exactly two pops must be listed.\n")
        }
      }
    }
  }
  else{
    if(is.null(sfs)){
      msg <- c("A SFS must be provided if snpRdata is not.\n")
    }
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  #===========prep if a SFS not provided===========
  if(is.null(sfs)){
    sfs <- calc_sfs(x, facet, pops, projection, fold = F)
  }
  
  # flip everthing, since i should be pop 1 and j should be pop 2, and the
  # matrix is set up so that rows are pop 1 and columns are pop 2. Transposing fixes it.
  pops <- attr(sfs, "pop")
  sfs <- t(sfs)

  # normalize, exluding fixed sites
  sfs <- sfs[-1,-1] # remove fixed sites
  sfs <- sfs/sum(sfs)

  # get all of the allele frequencies in each cell
  freqs.i <- (1:(nrow(sfs)))/(nrow(sfs))
  freqs.j <- (1:(ncol(sfs)))/(ncol(sfs))

  # do the math, this is equ 1b (i - j times sfs)
  directionality <- sum(outer(freqs.i, freqs.j, "-") * sfs)
  if(directionality > 0){
    attr(directionality, "direction") <- paste0(pops[1], "->", pops[2])
  }
  else{
    attr(directionality, "direction") <- paste0(pops[1], "<-", pops[2])

  }
  
  .yell_citation("Peter2013", "Direcitonality", "Directionality index (psi)", update_bib)

  return(directionality)
}
