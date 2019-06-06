#' Bootstrapped smoothed values.
#'
#' \code{resample_long} creates a distribution of bootstrapped smoothed values for a desired statistic.
#'
#'Bootstrapps are conducted as described by Hohenlohe et al. (2010), Population genomics of parallel adaptation in threespine stickleback using sequenced RAD tags.
#'
#'For each bootstrap, this function draws random window position, then draws random statistics from all provided SNPs to fill each observed position on that window and calculates a smoothed statistic for that window using gaussian-smoothing.
#'
#'This function can pick windows to bootstrap either from a provided data frame or around SNPs given a particular window size. It can also weight observations by random wieghts drawn from the provided SNP data. The function is built to run over one level, such as linkage group or chromosome.
#'
#' @param x Input data containing \emph{unsmoothed} statistic of interest and SNP position. Note that this data must contain a column titled "group" containing the linkage group/chromosome/scaffold for each observed value, since positions are drawn from the same group only.
#' @param statistic A character string designating the \emph{unsmoothed} statistic to smooth.
#' @param boots Number of bootstraps to run.
#' @param sigma Size variable in kb for gaussian smoothing. Full window size is 6*sigma.
#' @param nk_weight If true, weights smoothing by randomly drawing variables from a column titled "nk" in x.
#' @param fws If using a window with a fixed slide rather than snp centralized, provide a data frame of window positions.
#' @param n_snps If true, draws all random numbers up front. A column containing the number of snps in each selected window must be provided, titled "n_snps", in either x or fws. This allows for faster processing. Otherwise, random numbers will be drawn after this value is determined for each bootstrap.
#' @param report Progress report interval, in bootstraps.
#' @param level Character or NULL, default "group." What level to bootstrap across?
#' @return A numeric vector containing bootstrapped smoothed values.
#'
#' @examples
#' resample_long(randPI[randPI$pop == "A" & randPI$group == "chr1",], "pi", 100, 200, TRUE, randSMOOTHed[randSMOOTHed$pop == "A" & randSMOOTHed$group == "chr1",], TRUE, 10)
#'
do_bootstraps <- function(x, facets = NULL, boots, sigma, statistics = "all", nk = T, step = NULL, report = 10000, par = FALSE){
  browser()

  #================sanity checks================
  facets <- check.snpR.facet.request(x, facets, "none")

  if(length(facets) > 1){
    stop("Please provide only one facet at a time.\n")
  }

  pairwise.types <- c("fst")
  single.types <- c("pi", "ho", "pa", "pHWE")
  all.types <- c(pairwise.types, single.types)
  if(statistics == "all"){
    statistics <- all.types[which(all.types %in% c(colnames(x@stats), colnames(x@pairwise.stats)))]
  }

  stats.type <- character()
  if(any(statistics %in% pairwise.types)){
    stats.type <- "pairwise"
  }
  if(any(statistics %in% single.types)){
    stats.type <- c(stats.type, "stats")
  }

  sanity_check_window(x, sigma, 200, stats.type, nk, facets, statistics, good.types = all.types)

  # get mafs if doing any normal stats
  if("stats" %in% stats.type & nk){
    x <- calc_maf(x, facets)
    x@stats$nk <- x@stats$maj.count + x@stats$min.count
  }

  #===========subfunctions=======
  get.grps <- function(snps){
    grps <- x@snp.meta[snps, snp.facets]
    grps <- as.data.frame(grps)
    colnames(grps) <- snp.facets
    grps <- do.call(paste0, grps)
    return(grps)
  }

  # gets a list containing the position of snps within each window. Note that names are window centroids
  # and NAs are still in each list. These should be removed when finished.
  get.window.info <- function(cs, sig, position, cgrps = NULL, all.grps = NULL){
    ends <- cs + sig*3
    starts <- cs - sig*3
    lmat <- outer(position, starts, function(pos, starts) pos >= starts)
    lmat <- lmat + outer(position, ends, function(pos, ends) pos <= ends)
    if(!is.null(cgrps)){
      lmat <- lmat + outer(all.grps, cgrps, function(all.grps, cgrps) all.grps == cgrps)
    }
    colnames(lmat) <- cs
    rownames(lmat) <- position
    if(!is.null(cgrps)){
      lmat <- ifelse(lmat == 3, TRUE, FALSE)
    }
    else{
      lmat <- ifelse(lmat == 2, TRUE, FALSE)
    }

    n_snps <- colSums(lmat)

    # map positions to TRUEs
    lmat[lmat == FALSE] <- NA
    pos.mat <- matrix(position, nrow = nrow(lmat), ncol = ncol(lmat))
    lmat <- lmat * pos.mat
    lmat <- as.list(as.data.frame(lmat))

    return(list(lmat = lmat, n_snps = n_snps))
  }

  # runs the bootstrapping on a set of data.
  run.bootstrap.set <- function(fws, nrands, stats = NULL, pairwise.stats = NULL, stats.type, sample.facet, sample.subfacet = NULL, n_snps, comparison = NULL){
    browser()

    # grab the stats
    if("stats" %in% stats.type){
      if(!data.table::is.data.table(stats)){
        stats <- data.table::as.data.table(stats)
      }

      # grab the relevent facet only!
      if(length(sample.facet) > 0){
        stats <- stats[stats$facet == sample.facet,]
      }
      if(!is.null(sample.subfacet[1])){
        stats <- stats[stats$subfacet == sample.subfacet,]
      }
      if(nk){
        keep.cols <- c(statistics[which(statistics %in% single.types)], "nk")
      }
      else{
        keep.cols <- statistics[which(statistics %in% single.types)]

      }
      stats <- stats[nrands, ..keep.cols]
    }
    if("pairwise" %in% stats.type){
      if(!data.table::is.data.table(pairwise.stats)){
        pairwise.stats <- data.table::as.data.table(pairwise.stats)
      }

      # grab the relevent facet only!
      if(length(sample.facet) > 0){
        pairwise.stats <- pairwise.stats[pairwise.stats$facet == sample.facet,]
      }
      if(!is.null(comparison[1])){
        pairwise.stats <- pairwise.stats[pairwise.stats$comparison == comparison,]
      }

      if(nk){
        keep.cols <- c(statistics[which(statistics %in% pairwise.types)], "nk")
      }
      else{
        keep.cols <- statistics[which(statistics %in% pairwise.types)]
      }

      pairwise.stats <- pairwise.stats[nrands, ..keep.cols]
    }

    # work point
    cpos <- as.numeric(names(fws))
    spos <- as.numeric(unlist(fws))

    for (j in 1:length(fws)){

      #calculate smoothed statistics from the random stats and their positions.
      if(nk){
        gwp <- gaussian_weight(fws[[j]],as.numeric(names(fws[j]),sig))*rs
        gws <- gaussian_weight(tpos,cs[j],sig)
      }
      else{
        gwp <- gaussian_weight(tpos,cs[j],sig)*rs*(rnk - 1)
        gws <- gaussian_weight(tpos,cs[j],sig)*(rnk - 1)
      }
      if(is.na(sum(gwp, na.rm = T)/sum(gws, na.rm = T))){stop("Got a NA value for a bootstrapped smoothed window. This usually happens if the nk_weight is TRUE but no column titled nk is provided.")}
      smoothed_dist[j] <- (sum(gwp, na.rm = T)/sum(gws, na.rm = T)) #save the random stat
    }
  }

  #===========initialize=========
  #report a few things, intialize others
  sig <- 1000*sigma #correct sigma
  if(!is.null(step)){
    step <- 1000*step
  }
  num_snps <- nrow(x)
  #print(count_dist)
  smoothed_dist <- numeric(boots)
  cat("Sigma:", sigma, "\nTotal number of input snps: ", num_snps, "\nNumber of boots:", boots, "\n")

  position <- x@snp.meta$position

  snp.facets <- check.snpR.facet.request(x, facets, remove.type = "sample")

  #draw random centriods including the positions of those centroids and of the snps in those centroids.
  if(is.null(step)){
    csnps <- sample(1:nrow(x), boots, replace = T) #get random centroids
    cs <- position[csnps] #get centroid positions

    # get centroid and all snp groups for this facet
    if(!is.null(snp.facets[1])){
      cgrps <- get.grps(csnps)
      all.grps <- get.grps(1:nrow(x))
    }
    # draw all of the random snps to fill in around centroids on windows now. First need to figure out snp count in each window.
    ## to do this, figure out which snps are in the csnp windows.
    fws <- get.window.info(cs, sig, position, cgrps, all.grps)
    n_snps <- fws$n_snps
    fws <- lapply(fws$lmat, na.omit)

    ## draw the random snps
    nrands <- sample(1:nrow(x), length(unlist(fws)), replace = T)
  }
  else{ #same deal, but for fixed slide windows from provided window information
    browser()
    # go through each group and get the windows
    grps <- get.grps(1:nrow(x))
    u.grps <- unique(grps)
    fws <- vector("list")
    n_snps <- numeric()
    for(i in 1:length(u.grps)){
      t.pos <- x@snp.meta$position[grps == u.grps[i]]
      cs <- seq(0, max(t.pos), by = step)
      tfws <- get.window.info(cs, sig, t.pos)
      n_snps <- c(n_snps, tfws$n_snps)
      fws <- c(fws, tfws$lmat)
    }

    # clean up the window info list
    fws <- lapply(fws, na.omit)


    # get whatever info we can before passing to a looping function
    cwin <- sample(1:length(fws), boots, replace = T)
    fws <- fws[cwin]
    cs <- as.numeric(names(fws))

    nrands <- sample(1:nrow(x), sum(n_snps), replace = T) # these are the snps to pass
  }



  #===========run the bootstraps==========

  #do bootstraps
  if(par == FALSE){
    out <- run.bootstrap.set(fws = fws, nrands = nrands,
                             stats = x@stats, pairwise.stats = x@pairwise.stats,
                             stats.type =  stats.type, sample.facet = "pop", sample.subfacet = "ASP", comparison = "ASP~CLF",
                             n_snps = n_snps)
  }
  return(smoothed_dist)
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

