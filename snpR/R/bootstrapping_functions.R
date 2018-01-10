#Bootstraps a listed statistic from the genome to create a null distribution for p-values. Provide
#raw, unsmoothed data from whole genome. Kernal smooths to create distribution.
#inputs:  data: Input data. Contains columns titled "group", "position," one titled to match the "stat"
#               argument, and perhaps one titled "nk" and another titled "n_snps".
#         boots: number of bootstraps to do
#         sigma: Size variable in mb for gaussian smoothing. Full window size is 6*sigma
#         nk_weight: If true, looks for a column containing sample sizes titled "nk" to weight smoothing by.
#         fws: If using a window with a fixed slide rather than snp centralized, this is a data frame which
#              has columns titled "group" and "position", and potentially "n_snps".
#         n_snps: If true, looks for a column titled n_snps in data or fws and draws all random numbers up front.
#                 This allows for faster processing. Otherwise, random numbers will be drawn after this value is determined
#                 for each bootstrap.
#         report: Progress report interval, in bootstraps.
#While vectorized where possible, this is still a long process when boots is large.
#resample_long_par does the same process in parallel (albiet with a different random sampling approach).
#Due to read/write limitations, this often slower.

#' Bootstrapped smoothed values.
#'
#' \code{resample_long} creates a distribution of bootstrapped smoothed values for a desired statistic, as described by Hohenlohe et al. (2010), Population genomics of parallel adaptation in threespine stickleback using sequenced RAD tags. To run in parallel across multiple linkage groups/chromosomes, use  \code{\link{resample_long_par}}.
#'
#' @param x Input SNP data.
#' @param boots Number of bootstraps to run.
#' @param sigma Size variable in kb for gaussian smoothing. Full window size is 6*sigma.
#' @param nk_weight If true, weights smoothing by randomly drawing variables from a column titled "nk" in x.
#' @param fws If using a window with a fixed slide rather than snp centralized, provide a data frame of window positions.
#' @param n_snps If true, draws all random numbers up front. A column containing the number of snps in each selected window must be provided, titled "n_snps", in either x or fws. This allows for faster processing. Otherwise, random numbers will be drawn after this value is determined for each bootstrap.
#' @param report Progress report interval, in bootstraps.
#' @return A numeric vector containing bootstrapped smoothed values.
#'
#' @examples
#' resample_long(randPI[randPI$pop == "A" & randPI$group == "chr1",], "pi", 100, 200, TRUE, randSMOOTHed[randSMOOTHed$pop == "A" & randSMOOTHed$group == "chr1",], TRUE, 10)
#'
resample_long <- function(x, statistic, boots, sigma = 150, nk_weight = FALSE, fws = NULL, n_snps = T, report = 10000){

  #report a few things, intialize others
  sig <- 1000*sigma #correct sigma
  num_snps <- nrow(x)
  #print(count_dist)
  smoothed_dist <- numeric(boots)
  cat("Sigma:", sigma, "\nTotal number of input snps: ", num_snps, "\nNumber of boots:", boots, "\n")


  #draw random centriods, their positions, lgs of those centroids, and number of snps in window if provided
  if(is.null(fws)){
    csnps <- sample(1:nrow(x), boots, replace = T) #get random centroids
    cs <- x$position[csnps] #get centroid positions
    grps <- x$group[csnps] #get centroid groups
    if(n_snps){ #if n_snps is specified, draw all of the random snps to fill in around centroids on windows now.
      nrands <- sample(1:nrow(x),sum(x$n_snps[csnps]), replace = T)
      nrprog <- 1
    }
  }
  else{ #same deal, but for fixed slide windows from provided window information
    cwin <- sample(1:nrow(fws), boots, replace = T)
    cs <- fws$position[cwin]
    grps <- fws$group[cwin]
    if(n_snps){
      nrands <- sample(1:nrow(x),sum(fws$n_snps[cwin]), replace = T)
      nrprog <- 1
    }
  }

  #do bootstraps
  for (j in 1:boots){
    if(j %% report == 0){cat("Bootstrap number: ", j, "\n")}

    #figure out positions to use. Must be on the same group as the centriod and within 3*sigma
    tpos <- x$position[x$position >= cs[j] - 3*sig &
                            x$position <= cs[j] + 3*sig &
                            x$group == grps[j]]

    #if random snps haven't already been drawn...
    if(!n_snps){
      #draw random snps to fill those positions from ALL snps, with replacement
      rdraws <- sample(1:nrow(x), length(tpos), replace = T)
    }
    else{
      #pull the correct randomly chosen snps from the vector of random snps.
      rdraws <- nrands[nrprog:(nrprog + length(tpos) - 1)]
      nrprog <- nrprog + length(tpos)
    }


    #get random stats and nks
    rs <- x[rdraws, statistic] #sample stats for the randomly selected snps
    if (nk_weight == TRUE){
      rnk <- x$nk[rdraws] #sample nks for the randomly selected snps
    }


    #calculate smoothed statistics from the random stats and their positions.
    if(nk_weight == FALSE){
      gwp <- gaussian_weight(tpos,cs[j],sig)*rs
      gws <- gaussian_weight(tpos,cs[j],sig)
    }
    else{
      gwp <- gaussian_weight(tpos,cs[j],sig)*rs*(rnk - 1)
      gws <- gaussian_weight(tpos,cs[j],sig)*(rnk - 1)
    }
    if(is.na(sum(gwp, na.rm = T)/sum(gws, na.rm = T))){browser()}
    smoothed_dist[j] <- (sum(gwp, na.rm = T)/sum(gws, na.rm = T)) #save the random stat
  }
  return(smoothed_dist)
}


#Bootstraps a listed statistic from the genome to create a null distribution for p-values. Provide
#raw, unsmoothed data from whole genome. Kernal smooths to create distribution.
#inputs:  data: Input data. Contains columns titled "group", "position," one titled to match the "stat"
#               argument, and perhaps one titled "nk" and another titled "n_snps".
#         num_cores: Number of computing cores/threads to use.
#         boots: number of bootstraps to do
#         sigma: Size variable in mb for gaussian smoothing. Full window size is 6*sigma
#         nk_weight: If true, looks for a column containing sample sizes titled "nk" to weight smoothing by.
#         fws: If using a window with a fixed slide rather than snp centralized, this is a data frame which
#              has columns titled "group" and "position", and potentially "n_snps".
#This version runs in parallel and uses a slightly different approach. To avoid read/write lag,
#this version does all random draws locally within each loop.
resample_long_par <- function(x, num_cores, statistic, boots, sigma = 150, nk_weight = FALSE, fws = NULL){
  sig <- 1000*sigma
  num_snps <- nrow(x)
  #print(count_dist)
  smoothed_dist <- numeric(boots)
  cat("Sigma:", sigma, "\nTotal number of input snps: ", num_snps, "\nNumber of boots:", boots, "\n", "\nNumber of Cores:", num_cores, "\n")
  #draw random centriods, their positions, lgs of those centroids
  cat("WARNING: This may run slower than non-parallel since random numbers can't be drawn up-front due to read/write limitations.\n")

  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  output <- foreach(j = 1:boots, .export = 'gaussian_weight', .inorder = FALSE, .combine = "c") %dopar% {
    #get random window/centroid snp
    if(is.null(fws)){
      samp <- sample(1:nrow(x), 1)
      cs <- x$position[samp]
      grps <- x$group[samp]
    }
    else{
      samp <- sample(1:nrow(fws), 1)
      cs <- fws$position[samp]
      grps <- fws$group[samp]
    }

    #figure out positions to use
    tpos <- x$position[x$position >= cs - 3*sig &
                            x$position <= cs + 3*sig &
                            x$group == grps]


    #draw random snps to fill those positions, with replacement
    rdraws <- sample(1:nrow(x), length(tpos), replace = T)

    #get random stats and nks
    rs <- x[rdraws, statistic] #sample random stats from the x
    if (nk_weight == TRUE){
      rnk <- x$nk[rdraws] #sample random nks from the x
    }

    #print(sample)
    if(nk_weight == FALSE){
      gwp <- gaussian_weight(tpos,cs,sig)*rs
      gws <- gaussian_weight(tpos,cs,sig)
    }
    else{
      gwp <- gaussian_weight(tpos,cs,sig)*rs*(rnk - 1)
      gws <- gaussian_weight(tpos,cs,sig)*(rnk - 1)
    }
    sum(gwp, na.rm = T)/sum(gws, na.rm = T)
  }
  stopCluster(cl)
  registerDoSEQ()
  return(output)
}


#calculates a p-value for a stat based on a null distribution caluclated via bootstraps of that stat.
#inputs:  x: vector of observed statistics
#         dist: vector of statistics bootstrapped from same stat.
#         alt: Alternative hypothesis for p-value. Can be less, greater, or two-sided (default).
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

