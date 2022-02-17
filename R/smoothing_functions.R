#function for calculating gaussian weight of a point

#' Calculate gaussian weights.
#'
#' Get the gaussian weight of a statistic given its position relative to the
#' center and scaling factor sigma.
#'
#' @param p numeric. Position of value to smooth.
#' @param c numeric. Position of center of window.
#' @param s numeric. Sigma, scaling factor.
#'
#' @return A numeric value, the relative weight of the point.
gaussian_weight <- function(p, c, s) {
  exp(-(p-c)^2/(2*s^2))
}

#'Gaussian smooth statistics across sliding windows.
#'
#'Calculates gaussian smoothed average values for statistics and the genomic
#'position at which it was observed, splitting by any requested variables.
#'
#'Averages for multiple statistics can be calculated at once. If the statistics
#'argument is set to c("pairwise", "single"), all calculated stats will be run.
#'If it is set to "single", then all non-pairwise statistics (pi, ho, maf, ect)
#'will be bootstrapped, if it is set to "pairwise", then all pairwise statistics
#'(fst) will be bootstrapped. Individual statistics currently cannot be
#'requested by name, since the computational differences to add additional types
#'is minimal.
#'
#'The data can be broken up categorically by snp or sample metadata, as
#'described in \code{\link{Facets_in_snpR}}. Windows will only be calculated
#'using only SNPs on the same level of any provided facets. NULL and "all"
#'facets work as normally described in \code{\link{Facets_in_snpR}}.
#'
#'As described in Hohelohe et al. (2010), the contribution of individual SNPs to
#'window averages can be weighted by the number of observations per SNP by
#'setting the nk argument to TRUE, as is the default. For bootstraps, nk values
#'are randomly drawn for each SNP in each window.
#'
#'Centers for windows can either every SNP (if no step size is provided), or
#'every step kilobases from the 0 position of each snp level facet category
#'(chromosome, etc.).
#'
#'The size of sliding windows are defined by the "sigma" argument. Note that
#'this value, as well as that provided to the "step" arguement, are given in
#'kilobases. Each window will include SNPs within 3*sigma kilobases from the
#'window center. Past this point, the effect of each additional SNP on the
#'window average would be very small, and so they are dropped for computational
#'efficiency (see Hohenlohe (2010)).
#'
#'@param x snpRdata object.
#'@param facets character or NULL, default NULL. Categories by which to break up
#'  windows.
#'@param sigma numeric. Designates the width of windows in kilobases. Full
#'  window size is 6*sigma.
#'@param step numeric or NULL, default NULL. Designates the number of kilobases
#'  between each window centroid. If NULL, windows are centered on each SNP.
#'@param nk logical, default TRUE. If TRUE, weights SNP contribution to window
#'  averages by the number of observations at those SNPs.
#'@param stats.type character, default c("single", "pairwise"). Designates the
#'  statistic(s) to smooth, either "single",  "pairwise", or c("single",
#'  "pairwise"). See details.
#'@param par numeric or FALSE, default FALSE. If numeric, the number of cores to
#'  use for parallel processing.
#'
#'@export
#'@author William Hemstrom
#'
#'@references Hohenlohe et al. (2010). \emph{PLOS Genetics}
#'
#'@return snpRdata object with smoothed averages for any requested statistics
#'  merged into the window.stats or pairwise.window.stats slots.
#'
#'@examples
#'\dontrun{
#'# add a few statistics
#'x <- calc_pi(stickSNPs, "chr.pop")
#'x <- calc_ho(x, "chr.pop")
#'x <- calc_pairwise_fst(x, "chr.pop")
#'# smooth with a fixed slide between window centers.
#'x <- calc_smoothed_averages(x, "chr.pop", sigma = 200, step = 50)
#'get.snpR.stats(x, "chr.pop", "single.window") # pi, ho
#'get.snpR.stats(x, "chr.pop", "pairwise.window") # fst
#'}
calc_smoothed_averages <- function(x, facets = NULL, sigma, step = NULL, nk = TRUE, stats.type = c("single", "pairwise"), par = FALSE) {
  #==============sanity checks============
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  msg <- character(0)
  if(any(!stats.type %in% c("single", "pairwise"))){
    msg <- c(msg, "Unaccepted stats.type, only 'single' or 'pairwise' accepted.\n")
  }
  if("pairwise" %in% stats.type){
    if(nrow(x@pairwise.stats) == 0){
      msg <- c(msg, "No pairwise stats calculated.\n")
    }
  }

  if(length(msg) > 0){
    stop(msg)
  }
  
  sig <- 1000*sigma
  cat("Smoothing Parameters:\n\twindow size = ", 3*1000*sigma, "\n\tWindow slide = ", step*1000, "\n")

  .sanity_check_window(x, sigma, step, stats.type = stats.type, nk, facets = facets)
  
  x <- .add.facets.snpR.data(x, facets)
  #================subfunction========
  # funciton to do a sliding window analysis on a data set.
  # x: a data frame containing one column with positions and one for each column to be smoothed
  # nk: logical, should nk wieghting be performed (usually yes)
  # if ws is not FALSE, uses a typical sliding window. Otherwise does a window centered on each SNP.
  func <- function(x, step, sig, nk){
    if(nk){
      scols <- (which(colnames(x) == "nk") + 1):ncol(x)
      un.genotyped <- which(x$nk== 0)
      if(length(un.genotyped) > 0){
        x <- x[-which(x$nk == 0),]
      }
    }
    else{
      scols <- 2:ncol(x)
    }

    if(nrow(x) <= 1){
      ret <- matrix(NA, 0, 5 + length(scols))
      ret <- as.data.frame(ret)
      colnames(ret) <- c("position", "sigma", "n_snps", "step", "nk.status", colnames(x)[scols])
      return(ret)
    }
    #get window centers, starts, and stops:
    ##possible positions
    pos <- x$position
    pos <- as.numeric(pos)

    #window starts, stops, and ends
    if(!is.null(step)){
      cs <- seq(0, max(pos), by = step)
    }
    else{
      cs <- pos
    }
    ends <- cs + sig*3
    starts <- cs - sig*3

    #matrix where the rows are the snps and columns are window centers. Are the snps in the windows?
    lmat <- outer(pos, starts, function(pos, starts) pos >= starts)
    lmat <- lmat + outer(pos, ends, function(pos, ends) pos <= ends)
    colnames(lmat) <- cs
    rownames(lmat) <- pos
    lmat <- ifelse(lmat == 2, TRUE, FALSE)

    #get number of snps per window for later...
    n_snps <- colSums(lmat)

    #Multiply this using the gaussian weight function to get the contribution of each snp in each window.
    gmat <- outer(pos, cs, gaussian_weight, s = sig)
    gmat <- gmat*lmat

    #remove any windows with no SNPs
    gmat <- gmat[,n_snps != 0] #doing it this way allows for windows with SNPs but a stat value of 0.
    cs <- cs[n_snps != 0]
    n_snps <- n_snps[n_snps != 0]


    #fix any NA values
    vals <- as.matrix(x[,scols])
    NAs <- which(is.na(vals))
    if(length(NAs) > 0){
      vals[NAs] <- 0
    }

    #multiply by value of the statistics and nk if requested
    if(nk){
      #set up nks. Have to do it this way because of the possibility of NAs.
      nkv <- matrix(rep(x$nk, ncol(vals)), ncol = ncol(vals))
      nkv <- nkv - 1
      if(length(NAs) > 0){
        nkv[NAs] <- 0 #make sure that these snps don't contribute to the weight!
      }

      # fix any -1s (for ungenotyped stuff)
      mdats <- which(nkv == -1)
      if(any(mdats)){
        nkv[mdats] <- 0
      }

      #run
      win_vals <- t(gmat) %*% (vals*(x$nk - 1))
      win_scales <- t(gmat) %*% nkv

    }
    else{
      win_vals <- t(gmat) %*% vals
      win_scales <- rowSums(t(gmat))
    }

    #get the weighted value of the window
    win_stats <- win_vals/win_scales

    #return
    if(is.numeric(step)){
      out <- cbind(position = cs, sigma = sig/1000, n_snps = n_snps, step = step/1000, nk.status = nk, win_stats)
    }
    else{
      out <- cbind(position = cs, sigma = sig/1000, n_snps = n_snps, step = step, nk.status = nk, win_stats)
    }
    colnames(out)[-c(1:5)] <- colnames(x)[scols]
    return(out)
  }

  #================smooth at proper levels========
  #which cols hold the stats of interest?

  if(!is.null(step)){step <- step*1000}

  if("single" %in% stats.type){
    cat("\nSmoothing single group stats...")
    out <- .apply.snpR.facets(x = x,
                             facets = facets,
                             req =  "pos.all.stats",
                             fun = func,
                             case =  "ps.pf.psf",
                             par = par,
                             step = step,
                             sig = sig,
                             stats.type = "stats",
                             nk = nk)

    x <- .merge.snpR.stats(x, out, "window.stats")

    if("pairwise" %in% stats.type){
      cat("\nSmoothing pairwise stats...")
      out <- .apply.snpR.facets(x = x,
                               facets = facets,
                               req =  "pos.all.stats",
                               fun = func,
                               case =  "ps.pf.psf",
                               par = par,
                               step = step,
                               sig = sig,
                               stats.type = "pairwise",
                               nk = nk)
      return(.merge.snpR.stats(x, out, "pairwise.window.stats"))
    }
    else{
      return(x)
    }
  }

  else{
    cat("Smoothing pairwise stats...")
    out <- .apply.snpR.facets(x = x,
                             facets = facets,
                             req =  "pos.all.stats",
                             fun = func,
                             case =  "ps.pf.psf",
                             par = par,
                             step = step,
                             sig = sig,
                             stats.type = "pairwise",
                             nk = nk)
    return(.merge.snpR.stats(x, out, "pairwise.window.stats"))
  }
}
