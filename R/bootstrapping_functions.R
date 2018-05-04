#' Bootstrapped smoothed values.
#'
#' \code{resample_long} creates a distribution of bootstrapped smoothed values for a desired statistic.
#'
#'Bootstrapps are conducted as described by Hohenlohe et al. (2010), Population genomics of parallel adaptation in threespine stickleback using sequenced RAD tags. To run in parallel across multiple linkage groups/chromosomes, use \code{\link{resample_long_par}}.
#'
#'For each bootstrap, this function draws random window position, then draws random statistics from all provided SNPs to fill each observed position on that window and calculates a smoothed statistic for that window using gaussian-smoothing.
#'
#'This function can pick windows to bootstrap either from a provided data frame or around SNPs given a particular window size. It can also weight observations by random wieghts drawn from the provided SNP data. The function is built to run over multiple linkage groups/chromosomes.
#'
#' @param x Data.frame. Input data containing \emph{unsmoothed} parameter of interest and SNP position.
#' @param parm Character string. The name of the \emph{unsmoothed} paramter to smooth. Must match a column name in x.
#' @param boots Numeric. Number of bootstraps to run.
#' @param sigma Numeric. Size variable in kb for gaussian smoothing. Full window size is 6*sigma.
#' @param nk_weight Logical, default TRUE. If true, weights smoothing by randomly drawing variables from a column titled "nk" in x.
#' @param fws Data.frame, default NULL. If using a window with a fixed slide rather than snp centralized, provide a data frame of window positions, such as that given by s_ave_multi. Must contain a column named "position".
#' @param levels Character, default "group". Additional levels by which to segregate SNPs to windows other than position. When constructing windows, SNPs in x that will be used for their position in bootstrapped windows must also be in the same levels as the window center. Must match colnumn names in both x and fws (if provided).
#' @return A numeric vector containing bootstrapped smoothed values.
#'
#' @examples
#' #with provided windows
#' temp <- cbind(stickSTATs, stickFORMATs$ac$n_total)
#' colnames(temp)[8] <- "nk"
#' sm <- s_ave_multi(temp, "pi", 200, 150)
#' resample_long(temp[temp$pop == "ASP",], "pi", 100, 200, T, fws = sm[sm$pop == "ASP",])
#'
#' #without provided windows, centered on each SNP instead.
#' temp <- cbind(stickSTATs, stickFORMATs$ac$n_total)
#' colnames(temp)[8] <- "nk"
#' resample_long(temp[temp$pop == "ASP",], "pi", 100, 200, T)
#'
resample_long <- function(x, parm, boots, sigma, nk_weight = TRUE, fws = NULL, levels = "group"){

  #report a few things, intialize others
  sig <- 1000*sigma #correct sigma
  num_snps <- nrow(x)
  #print(count_dist)
  cat("Sigma:", sigma, "\nTotal number of input snps: ", num_snps, "\nNumber of boots:", boots, "\n")
  scol <- which(colnames(x) == parm)


  #draw random window centers
  if(is.null(fws)){
    csnps <- sample(1:nrow(x), boots, replace = T) #get random centroids
    cs <- x$position[csnps] #get centroid positions
    # if(n_snps){ #if n_snps is specified, draw all of the random snps to fill in around centroids on windows now.
    #   nrands <- sample(1:nrow(x),sum(x$n_snps[csnps]), replace = T)
    #   nrprog <- 1
    # }
  }
  else{ #same deal, but for fixed slide windows from provided window information
    cwin <- sample(1:nrow(fws), boots, replace = T)
    cs <- fws$position[cwin]
    # if(n_snps){
    #   nrands <- sample(1:nrow(x),sum(fws$n_snps[cwin]), replace = T)
    #   nrprog <- 1
    # }
  }

  #make a matrix noting if each snp is in each random window. Code is similar to that in s_ave_multi.
  ends <- cs + sig*3
  starts <- cs - sig*3
  pos <- x$position
  pos <- as.numeric(pos)

  #matrix where the rows are the snps and columns are window centers. Are the snps in the windows?
  lmat <- outer(pos, starts, function(pos, starts) pos >= starts) #are the snp pos greater than the start?
  lmat <- lmat*outer(pos, ends, function(pos, ends) pos <= ends) #are the snp pos less than the end?
  if(!is.null(levels)){
    for(i in 1:length(levels)){ #are the snps in the same data level?
      tlev <- levels[i] #this level
      lcol <- which(colnames(x) == tlev) #level index containing info on this column
      slev <- x[,lcol] #vector of that level
      #vector of window levels, from fws or not
      if(is.null(fws)){
        wlev <- x[csnps, lcol]
      }
      else{
        lcol <- which(colnames(fws) == tlev)
        wlev <- fws[cwin, lcol]
      }
      lmat <- lmat*outer(slev, wlev, function(slev, wlev) slev == wlev) #check for equality
    }
  }

  #get smoothed windows, as before.
  gmat <- outer(pos, cs, gaussian_weight, s = sig)
  gmat <- gmat*lmat

  #remove any snps that aren't in a window
  gmat <- gmat[-(which(rowSums(gmat) == 0)),]

  #assign stat values to SNPs and fix any NA values
  vals <- sample(x[,scol], size = nrow(gmat), replace = T)
  NAs <- which(is.na(vals)) #which are NA?
  if(length(NAs) > 0){
    vals[NAs] <- 0 #fix the NAs
  }

  #multiply by value of the statistics and nk if requested
  if(nk_weight){
    #get random nk values
    nkv <- sample(x$nk, nrow(gmat), replace = T)

    #set up nks. Have to do it this way because of the possibility of NAs.
    nkv <- nkv - 1
    if(length(NAs) > 0){
      nkv[NAs] <- 0 #make sure that these snps don't contribute to the weight!
    }
    #run
    win_vals <- t(gmat) %*% (vals*nkv)
    win_scales <- t(gmat) %*% nkv
  }
  else{
    win_vals <- t(gmat) %*% vals
    win_scales <- rowSums(t(gmat))
  }

  #get the weighted value of the window
  win_stats <- win_vals/win_scales

  return(as.vector(win_stats))
}


#calculates a p-value for a stat based on a null distribution caluclated via bootstraps of that stat.
#inputs:  x: vector of observed statistics
#         dist: vector of statistics bootstrapped from same stat.
#         alt: Alternative hypothesis for p-value. Can be less, greater, or two-sided (default).


#' Caculate p-values from bootstrapped distributions.
#'
#' \code{p_calc_boots} finds p-values for observed \emph{smoothed} statistics from bootstrapped distributions, such as that produced by \code{\link{resample_long}}.
#'
#' Calculates p-values for smoothed values of a statistic based upon a bootstrapped null distribution of that statistic, such as those produced by resample_long using an emperical continuous distribution function.
#'
#' @param x Numeric vector of smoothed values for which to calculate p-values.
#' @param dist: Numeric vector of bootstrapped smoothed values.
#' @param alt: Character string. Which alternative hypothesis should be used? Options: #' \itemize{
#'    \item "less": probability that a bootstrapped value is as small or smaller than observed.
#'    \item "greater": probability that a bootstrapped value is as large or larger than observed.
#'    \item "two-sided": probability that a bootstrapped value is as or more extreme than observed.
#' }
#' @return A numeric vector containing p-values for each input smoothed value.
#'
#' @examples
#' p_calc_boots(randSMOOTHed[randSMOOTHed$pop == "A",]$smoothed_pi, randPIBOOTS)
#' p_calc_boots(randSMOOTHed[randSMOOTHed$pop == "A",]$smoothed_pi, randPIBOOTS, alt = "less")
#'
p_calc_boots <- function(x, dist, alt = "two-sided"){
  if(!alt %in% c("less", "greater", "two-sided")){
    stop("Alt must be one of character strings: less, greater, or two-sided.\n")
  }

  #create a cumulative distibution function of the distribution
  edist <- ecdf(dist)

  #get the proportion of the bootstrapped distribution less than or equal to the observed point.
  xp <- edist(x)

  #get p-values
  if(alt == "two-sided"){
    #Take abs(.5 - xp), which is the proportion of the distribution between the point and the distribution mean.
    #Multiply that by 2 to get the proportion of the distribution between that point and the distirbution mean
    #as well as between the mean and a point equally extreme in the other direction.
    #take 1 minus that to get the proportion of the as or more extreme
    xp <- 1 - (abs(xp - .5)*2)
  }
  else if (alt == "greater"){
    #get the proportion of data greater than observed data points
    xp <- 1 - xp
  }
  #if alt is "less", no conversion needed. already gives the proportion less than observed.

  return(xp)
}

