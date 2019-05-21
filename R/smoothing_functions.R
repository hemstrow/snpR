#function for calculating gaussian weight of a point

#'Gaussian weights.
#'
#'\code{gaussian_weight} Function to get the gaussian weight of a statistic given its position relative to the center and scaling factor sigma. Used internally, probably shouldn't be called.
#'
#' @param p Position of value to smooth.
#' @param c Position of center of window.
#' @param s Sigma, scaling factor.
#'
#' @return A numeric value, the relative weight of the point.
#'
#' @keywords internal
#'
#' @examples
#' gaussian_weight(10000, 10500, 1000)
gaussian_weight <- function(p, c, s) {
  exp(-(p-c)^2/(2*s^2))
}

#function to generate sliding averages based on gaussian weights of points around each point,
#cut off at 3sigma where weight becomes very low. Sigma should be given in kbp, although input should be
# in bp, in the second column of each row.

#'Gaussian smooth statistics.
#'
#'Wrapper for \code{\link{s_ave_multi}} called with no levels to split by.
#'
#' @param x Input SNP data frame.
#' @param parameter Name of the statistic to smooth.
#' @param sigma Smoothing statistic/window size, in kb.
#' @param nk_weight Should statistic contribution to window mean be additionally scaled by column "nk"?
#' @param fixed_window Should a window with a fixed slide distance be used? If so, provide window slide length in kb.
#'
#' @return If fixed_window == FALSE, returns the input data frame with an additional column containing smoothed average. Otherwise, returns a new data frame containing the position info for each window and smoothed value of the statistic at that window. Both outputs will also contain a column with the number of SNPs in each window.
#'
#' @examples
#' #fixed slide window:
#' smoothed_ave(randPI[randPI$pop == "A" & randPI$group == "chr1",], "pi", 200, TRUE, fixed_window = 50)
#'
#' #windows centered on each snp
#' smoothed_ave(randPI[randPI$pop == "A" & randPI$group == "chr1",], "pi", 200, TRUE)
#'
#' #wrapped in run_gp
#' run_gp(randPI, smoothed_ave, parameter = "pi", sigma = 200, nk_weight = TRUE, fixed_window = 50)
#'
smoothed_ave <- function(x, parameter, sigma, nk_weight = FALSE, fixed_window = NULL) {
  out <- s_ave_multi(x, parameter, sigma, fixed_window, nk_weight, levs = NA)
}



#'Gaussian smooth multiple statistics
#'
#'\code{sm_ave_multi} Calculates gaussian smoothed average values for statistics and the genomic position at which it was observed, splitting by any requested variables.
#'
#'Description of x:
#'    Requires columns titled columns titled "position", columns titled to match those in the "parms" argument, and those to match the "levs" argument if set. These contain positions in bp, the statistic to be smoothed for each SNP, and any variables to split the smoothing by (such as chromosome and population). If nk is TRUE, a column titled "nk" or "n_total" is required containing the weighting factor for smoothing, typically the number of observed alleles/sample size, "nk".
#'
#' @param x Input SNP data frame.
#' @param parms Character vector. Names of the statistics to smooth.
#' @param sigma Numeric value. Smoothing statistic/window size, in kb.
#' @param ws Numeric value or NULL, default NULL. Window slide length. If NULL, uses every SNP as a window center.
#' @param nk Boolean, default TRUE. Should statistic contribution to window mean be additionally scaled by column "nk"?
#' @param levs Character vector or NA, default c('group', 'pop')). Names of the columns to split the data by for smoothing.
#'
#' @return If fixed_window == FALSE, returns the input data frame with an additional column containing smoothed average. Otherwise, returns a new data frame containing the position info for each window and smoothed value of the statistic at that window. Both outputs will also contain a column with the number of SNPs in each window.
#'
#' @examples
#' #data prep:
#' t2 <- calc_pi(stickFORMATs$ac)
#' t2 <- cbind(stickFORMATs$ac, pi = t2)
#' t3 <- calc_Ho(stickFORMATs$character, 3, pop = l)
#' t3 <- reshape2::melt(t3, id.vars = c("snp", "position", "group"))
#' t2$ho <- t3$value
#'
#'
#'
#' #fixed slide window, splitting by pop and group, the defualt:
#' s_ave_multi(t2, c("pi", "ho"), 200, 150, TRUE)
#'
#' #splitting by only pop.
#' s_ave_multi(t2, c("pi", "ho"), 200, 150, TRUE, "pop")
#'
#' #no splitting
#' s_ave_multi(t2, c("pi", "ho"), 200, 150, TRUE, NA)
#'
calc_smoothed_averages <- function(x, facets, sigma, step = NULL, nk = TRUE, stats.type = c("stats", "pairwise"), par = FALSE) {
  sig <- 1000*sigma
  cat("Smoothing Parameters:\n\twindow size = ", 3*1000*sigma, "\n\tWindow slide = ", step*1000, "\n")

  sanity_check_window(x, sigma, step, stats.type = stats.type, nk)

  #================subfunction========
  #funciton to do a sliding window analysis on a data set. Vectorized!:
  # x: a data frame containing one column with positions and one for each column to be smoothed
  # nk: logical, should nk wieghting be performed (usually yes)
  # if ws is not FALSE, uses a typical sliding window. Otherwise does a window centered on each SNP.
  func <- function(x, step, sig, nk){
    scols <- (which(colnames(x) == "nk") + 1):ncol(x)
    un.genotyped <- which(x$nk== 0)
    if(length(un.genotyped) > 0){
      x <- x[-which(x$nk == 0),]

    }
    if(nrow(x) <= 1){
      ret <- matrix(NA, 0, 3 + length(scols))
      ret <- as.data.frame(ret)
      colnames(ret) <- c("position", "sigma", "n_snps", colnames(x)[scols])
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
    out <- cbind(position = cs, sigma = sig/1000, n_snps = n_snps, win_stats)
    colnames(out)[-c(1:3)] <- colnames(x)[scols]
    return(out)
  }

  #================smooth at proper levels========
  #which cols hold the stats of interest?

  if(!is.null(step)){step <- step*1000}

  if("stats" %in% stats.type){
    cat("\nSmoothing typical stats...")
    out <- apply.snpR.facets(x = x,
                             facets = facets,
                             req =  "pos.all.stats",
                             fun = func,
                             case =  "ps.pf.psf",
                             par = par,
                             step = step,
                             sig = sig,
                             stats.type = "stats",
                             nk = nk)

    x <- merge.snpR.stats(x, out, "window.stats")

    if("pairwise" %in% stats.type){
      cat("\nSmoothing pairwise stats...")
      out <- apply.snpR.facets(x = x,
                               facets = facets,
                               req =  "pos.all.stats",
                               fun = func,
                               case =  "ps.pf.psf",
                               par = par,
                               step = step,
                               sig = sig,
                               stats.type = "pairwise",
                               nk = nk)
      return(merge.snpR.stats(x, out, "pairwise.window.stats"))
    }
    else{
      return(x)
    }
  }

  else{
    cat("Smoothing pairwise stats...")
    out <- apply.snpR.facets(x = x,
                             facets = facets,
                             req =  "pos.all.stats",
                             fun = func,
                             case =  "ps.pf.psf",
                             par = par,
                             step = step,
                             sig = sig,
                             stats.type = "pairwise",
                             nk = nk)
    return(merge.snpR.stats(x, out, "pairwise.window.stats"))
  }
}
