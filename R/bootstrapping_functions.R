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
do_bootstraps <- function(x, facets = NULL, boots, sigma, step = NULL, statistics = "all", nk = T, par = FALSE){
  #note: it is possible to run all sample level facets at once, so something like c("pop.fam.group", "pop.group") can
  #      be run simultainously, with no need to loop across facets.
  #      However, SNP level facets create different windows, and so need to be run seperately. Essentially,
  #      we need to do everything once for each unique snp level facet defined in the data.

  #================sanity checks================
  facets <- check.snpR.facet.request(x, facets, "none")

  # figure out which stats we are using!
  pairwise.types <- c("fst")
  single.types <- c("pi", "ho", "pa", "pHWE")
  all.types <- c(pairwise.types, single.types)
  if(statistics == "all"){
    statistics <- all.types[which(all.types %in% c(colnames(x@stats), colnames(x@pairwise.stats)))]
  }
  else if(statistics == "single"){
    statistics <- single.types
  }
  else if(statistics == "pairwise"){
    statisitcs <- pairwise.types
  }

  stats.type <- character()
  if(any(statistics %in% pairwise.types)){
    stats.type <- "pairwise"
  }
  if(any(statistics %in% single.types)){
    stats.type <- c(stats.type, "stats")
  }

  # run basic sanity checks
  sanity_check_window(x, sigma, 200, stats.type, nk, facets, statistics, good.types = all.types)

  # get mafs if doing any normal stats. Can make this more effeicent by adding only where missing in the future.
  if("stats" %in% stats.type & nk){
    x <- calc_maf(x, facets)
    x@stats$nk <- x@stats$maj.count + x@stats$min.count
  }

  #===========subfunctions=======
  get.grps <- function(snps, facet){
    if(facet == ".base"){
      grps <- rep(".base", length(snps))
      return(grps)
    }
    grps <- x@snp.meta[snps, facet]
    grps <- as.data.frame(grps)
    colnames(grps) <- facet
    grps <- do.call(paste0, grps)
    return(grps)
  }

  # gaussian weights on data
  do.gaus <- function(fws, stats, nk, tnk, sig){

    if(nk){
      gws <- matrix(gaussian_weight(fws[[1]], as.numeric(names(fws)), sig), nrow(stats), ncol(stats))* (tnk - 1)
    }

    else{
      gws <- matrix(gaussian_weight(fws[[1]], as.numeric(names(fws)), sig), nrow(stats), ncol(stats))
    }

    gwp <- gws * stats

    return((colSums(gwp, na.rm = T)/colSums(gws, na.rm = T)))
  }

  # runs the bootstrapping on a set of data.
  # this function should be given all of the stats in wide form, so multiple sample level facets can be done at once!
  run.bootstrap.set <- function(fws, nrands, stats = NULL, stats.type, n_snps, nk, snk, pnk, part.cols){
    # zero fixer subfunction
    fix.zeros <- function(x){
      zeros <- x == 0
      if(any(zeros)){
        x[zeros] <- 1
      }
      return(x)
    }

    #==========run the bootstrap=======
    cpos <- as.numeric(names(fws))
    spos <- as.numeric(unlist(fws))

    # fix zeros and make a simple nk for when using just one stat.
    if(length(part.cols) == 1){
      if(!is.null(pnk[1,1])){
        use.nk <- pnk
      }
      else{
        use.nk <- snk
      }
      use.nk <- fix.zeros(use.nk)
    }
    else{
      snk <- fix.zeros(snk)
      pnk <- fix.zeros(pnk)
    }

    tracker <- 1
    percentage <- 0
    stats.out <- matrix(0, boots, ncol(stats))
    colnames(stats.out) <- colnames(stats)
    cat("Bootstrap Progress:\n0%\n")
    for (j in 1:boots){
      npercentage <- j/boots
      if(npercentage - percentage >= 0.05){
        cat(paste0(round(npercentage*100), "%"), "\n")
        percentage <- npercentage
      }

      trows <- tracker:(tracker + n_snps[j] - 1)

      # figure out nk
      if(!is.null(part.cols$stats) & !is.null(part.cols$pairwise)){
        tnk <- cbind(snk, pnk)[trows,]
      }
      else{
        tnk <- use.nk[trows,]
      }

      # get smoothed stats and save output
      stats.out[j,] <- do.gaus(fws[j], stats[trows,], nk, tnk, sig)

      tracker <- trows[length(trows)] + 1
    }

    # prepare and return output
    return(stats.out)
  }

  # gets a list containing the position of snps within each window. Note that names are window centroids
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

  # prepares a single snp facet to run
  prep.one.snp.facet <- function(x, facet.info, stats = NULL, pairwise.stats = NULL){
    snk <- NULL
    pnk <- NULL
    #============prepare a centroid list=========
    # elements named for all of the random centroids, contain a vector of the positions of snps in those windows
    # make sure to save a vector of the number of snps per window (n_snps)
    if(is.null(step)){
      csnps <- sample(1:nrow(x), boots, replace = T) #get random centroids
      cs <- position[csnps] #get centroid positions

      # get centroid and all snp groups for this facet
      if(!is.null(snp.facets[1])){
        cgrps <- get.grps(csnps, names(facet.info))
        all.grps <- get.grps(1:nrow(x), names(facet.info))
      }
      # draw all of the random snps to fill in around centroids on windows now. First need to figure out snp count in each window.
      ## to do this, figure out which snps are in the csnp windows.
      fws <- get.window.info(cs, sig, position, cgrps, all.grps)
      n_snps <- fws$n_snps


      # remove any empty windows
      empties <- n_snps == 0
      if(any(empties)){
        fws <- fws[-which(empties)]
        n_snps <- n_snps[-which(empties)]
      }

      fws <- lapply(fws$lmat, na.omit)

      ## draw the random snps
      nrands <- sample(1:nrow(x), length(unlist(fws)), replace = T)
    }
    else{ #same deal, but for fixed slide windows from provided window information
      # go through each group and get the windows
      grps <- get.grps(1:nrow(x), names(facet.info))
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

      # remove any empty windows
      empties <- n_snps == 0
      if(any(empties)){
        fws <- fws[-which(empties)]
        n_snps <- n_snps[-which(empties)]
      }

      # clean up the window info list
      fws <- lapply(fws, na.omit)


      # get whatever info we can before passing to a looping function
      cwin <- sample(1:length(fws), boots, replace = T)
      fws <- fws[cwin]
      n_snps <- n_snps[cwin]
      cs <- as.numeric(names(fws))

    }

    #============draw random snps to assign to each position in each window========
    nrands <- sample(1:nrow(x), sum(n_snps), replace = T)

    #============prepare the stats and/or pairwise stats=========
    # put them into long form across sample level facets
    # cbind stats and pairwise stats, creating an nk vector to pass and saving facet and subfacet/comparison data
    facet.meta <- character()
    subfacet.meta <- facet.meta
    stattype.meta <- character()

    # get nk info and cast each type
    if("stats" %in% stats.type){
      if(!data.table::is.data.table(stats)){
        stats <- data.table::as.data.table(stats)
      }

      # get only data for the correct sample facets
      stats <- stats[stats$facet %in% unlist(facet.info),]


      # melt to put different stat/facet+subfacets in one row via snp id
      stat.cols <- statistics[which(statistics %in% single.types)]
      meta.cols <- c("facet", "subfacet", ".snp.id")
      keep.cols <- c(meta.cols, stat.cols)
      stattype.meta <- rep("single", length(facet.meta))
      if(nk){ # grab a cast nk as well...
        nk.meta <- c(meta.cols, "nk")
        snk <- data.table::dcast(stats[,..nk.meta], .snp.id~facet + subfacet, value.var = "nk")
        snk <- snk[nrands,-".snp.id"]
        ## need to duplicate columns, since each nk value will be used for multiple stats!
        col.seq <- rep(1:ncol(snk), sum(statistics %in% single.types))
        snk <- snk[,..col.seq]
      }
      stats <- stats[,..keep.cols]
      stats <- data.table::dcast(stats, .snp.id~facet + subfacet, value.var = stat.cols)

      # subset according to the snps we are using
      stats <- stats[nrands, -".snp.id"]

      # fix column names if needed.
      if(sum(statistics %in% single.types) == 1){
        colnames(stats) <- paste0(statistics[which(statistics %in% single.types)], "_", colnames(stats))
      }
    }
    if("pairwise" %in% stats.type){
      if(!data.table::is.data.table(pairwise.stats)){
        pairwise.stats <- data.table::as.data.table(pairwise.stats)
      }

      # get only data for the correct sample facets
      pairwise.stats <- pairwise.stats[pairwise.stats$facet %in% unlist(facet.info),]

      # melt to put different stat/facet+subfacets in one row via snp id
      stat.cols <- statistics[which(statistics %in% pairwise.types)]
      meta.cols <- c("facet", "comparison", ".snp.id")
      keep.cols <- c(meta.cols, stat.cols)
      if(nk){ # grab a cast nk as well...
        nk.meta <- c(meta.cols, "nk")
        pnk <- data.table::dcast(pairwise.stats[,..nk.meta], .snp.id~facet + comparison, value.var = "nk")
        pnk <- pnk[nrands,-".snp.id"]
        col.seq <- rep(1:ncol(pnk), sum(statistics %in% pairwise.types))
        pnk <- pnk[,..col.seq]
      }
      pairwise.stats <- pairwise.stats[,..keep.cols]
      pairwise.stats <- data.table::dcast(pairwise.stats, .snp.id~facet + comparison, value.var = stat.cols)

      # subset according to the snps we are using
      pairwise.stats <- pairwise.stats[nrands, -".snp.id"]

      # fix column names if needed.
      if(sum(statistics %in% pairwise.types) == 1){
        colnames(pairwise.stats) <- paste0(statistics[which(statistics %in% pairwise.types)], "_", colnames(pairwise.stats))
      }
    }

    #============bind and prepare column type info================
    part.cols <- vector("list", length(stats.type))
    names(part.cols) <- stats.type
    if(any(statistics %in% single.types) & any(statistics %in% pairwise.types)){
      bound.stats <- cbind(stats, pairwise.stats)
      part.cols$stats <- ncol(stats)
      part.cols$pairwise <- ncol(pairwise.stats)
    }
    else if(any(statistics %in% single.types)){
      bound.stats <- stats
      part.cols$stats <- ncol(stats)
    }
    else{
      bound.stats <- pairwise.stats
      part.cols$pairwise <- ncol(pairwise.stats)
    }

    return(list(fws = fws, n_snps = n_snps, nrands = nrands, stats = bound.stats, snk = snk, pnk = pnk, part.cols = part.cols))
  }

  # melt a bootstrap output back into a long form and add some extra metadata
  melt.bootstrap <- function(boot.set, sigma, nk, step){
    # melt data
    boot.set <- data.table::as.data.table(boot.set)
    suppressWarnings(boot.set <- data.table::melt(boot.set))

    # get metadata
    facet <- unlist(strsplit(as.character(boot.set$variable), "_"))
    stats <- facet[seq(1, length(facet), by = 3)]
    subfacet <- facet[seq(3, length(facet), by = 3)]
    facet <- facet[seq(2, length(facet), by = 3)]

    # set metadata
    set(boot.set, j = "facet", value = facet)
    set(boot.set, j = "subfacet", value = subfacet)
    set(boot.set, j = "stat", value = stats)
    set(boot.set, j = "sigma", value = sigma)
    set(boot.set, j = "nk", value = nk)
    if(is.null(step)){
      set(boot.set, j = "step", value = NA)
    }
    else{
      set(boot.set, j = "step", value = step/1000)
    }

    # rearrange
    ord <- c(3, 4, 5, 6, 7, 8, 2)
    boot.set <- boot.set[,..ord]

    return(boot.set)
  }

  # wrapper function to run a single facet for a given number of boots
  boot.wrapper <- function(x, nk, stats.type, sigma, step){
    boot.out <- run.bootstrap.set(fws = x$fws,
                                  nrands = x$nrands,
                                  stats = x$stats,
                                  nk = nk,
                                  snk = x$snk, pnk = x$pnk,
                                  stats.type = stats.type,
                                  n_snps = x$n_snps,
                                  part.cols = x$part.cols)
    boot.out <- melt.bootstrap(boot.out, sigma, nk, step)
    return(boot.out)
  }

  #===========initialize=========
  # figure out what different snp level facets we have, and figure out which facets have which snp levels
  snp.facets <- check.snpR.facet.request(x, facets, remove.type = "sample")
  snp.facet.matches <- lapply(snp.facets, grep, x = facets)
  names(snp.facet.matches) <- snp.facets
  snp.facet.matches <- lapply(snp.facet.matches, function(y){
    lapply(facets[y], function(z) check.snpR.facet.request(x, z))
  })
  snp.facet.matches <- lapply(snp.facet.matches, unlist) # this is now a list with an entry for each snp level facet containing the pop level facets to run.

  # add in .base if needed, for sample level only facets
  samp.facets <- check.snpR.facet.request(x, facets, "none", T)
  samp.facets <- samp.facets[[1]][which(samp.facets[[2]] == "sample")]
  if(length(samp.facets) > 0){
     snp.facet.matches <- c(snp.facet.matches, .base = samp.facets)
  }

  # for snp level only facets, change the NULL to ".base"
  nulls <- which(unlist(lapply(snp.facet.matches, function(x) any(is.null(x)))))
  if(length(nulls) > 0){
    snp.facet.matches[[nulls]] <- ".base"
  }

  #report a few things, intialize others
  sig <- 1000*sigma #correct sigma
  if(!is.null(step)){
    step <- 1000*step
  }
  num_snps <- nrow(x)
  cat("Sigma:", sigma, "\nTotal number of input snps: ", num_snps, "\nNumber of boots:", boots, "\n")

  position <- x@snp.meta$position


  #===========run a loop to do bootstraps for each snp level facet==============
  out <- vector("list", length(snp.facet.matches))
  names(out) <- names(snp.facet.matches)
  cat("Unique snp level facets:", paste0(names(snp.facet.matches), collapse = ", "), "\n")

  for(i in 1:length(snp.facet.matches)){
    cat("Bootstrap facet:\n", names(snp.facet.matches[i]), "\n")

    # prepare data
    prepped.facet <- prep.one.snp.facet(x, snp.facet.matches[i], x@stats, x@pairwise.stats)

    # run in serial
    if(par == FALSE){
      out[[i]] <- boot.wrapper(x = prepped.facet,
                               nk = nk,
                               stats.type = stats.type,
                               sigma = sigma,
                               step = step)
      set(out[[i]], j = "snp.facet", value = names(snp.facet.matches[i]))
    }
    # run in parallel
    else{
      library(doParallel);library(foreach)
      cl <- snow::makeSOCKcluster(par)
      doSNOW::registerDoSNOW(cl)

      #prepare reporting function
      task_list <- split(1:boots, sort((1:boots)%%par))
      ntasks <- par
      progress <- function(n) cat(sprintf("\tPart %d out of", n), ntasks, "is complete.\n")
      opts <- list(progress=progress)

      tout <- foreach::foreach(q = 1:ntasks, .inorder = FALSE,
                               .options.snow = opts, .export = "data.table") %dopar% {
                                 # extract the correct parts of the prepped.facet
                                 tpf <- prepped.facet
                                 tpf$fws <- tpf$fws[task_list[[q]]]
                                 tpf$n_snps <- tpf$n_snps[task_list[[q]]]

                                 # for nrands, need to figure out where to start and end
                                 if(q != 1){
                                   prior.sum <- sum(prepped.facet$n_snps[1:(task_list[[q]][1] - 1)])
                                 }
                                 else{
                                   prior.sum <- 1
                                 }
                                 post.sum <- sum(prepped.facet$n_snps[1:(task_list[[q]][length(task_list[[q]])])])

                                 tpf$nrands <- tpf$nrands[prior.sum:post.sum]

                                 # do the bootstrapp
                                 boots <- length(task_list[[q]])
                                 boot.wrapper(x = tpf,
                                              nk = nk,
                                              stats.type = stats.type,
                                              sigma = sigma,
                                              step = step)
                               }

      #release cores
      parallel::stopCluster(cl)
      doSNOW::registerDoSNOW()

      # bind
      out[[i]] <- data.table::rbindlist(tout)
      set(out[[i]], j = "snp.facet", value = names(snp.facet.matches[i]))
    }
  }

  # bind and return
  out <- data.table::rbindlist(out)
  col.ord <- c(1, 2, ncol(out), 3:(ncol(out) - 1))
  out <- out[,..col.ord]
  x@window.bootstraps <- rbind(x@window.bootstraps, out)
  return(x)
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
calc_p_from_bootstraps <- function(x, facet = NULL, statistics = "all", alt = "two-sided", par = FALSE){
  #===========subfunctions=======
  # finds matching rows for given facet, ect (all of which are inherited)
  get.matches <- function(y){
    matches <- y$facet == facet &
      y$subfacet == subfacet &
      y$snp.facet == snp.facet &
      y$sigma == sigma &
      y$nk == nk &
      y$step == step
    if(any(colnames(y) == "stat")){
      matches <- (y$stat == statistic) & matches
    }
    return(which(matches))
  }

  # extracts an edf with the matching bootstrap info and calculates p-values
  get.one.pvalue <- function(x, facet, subfacet, snp.facet, statistic, nk, step, sigma, alt){

    # grab the matching raw statistical data
    matches <- get.matches(x@window.bootstraps)


    # create a cumulative distibution function of the distribution
    edist <- ecdf(x@window.bootstraps[matches,]$value)

    # get p-values
    ## get the matching statistic data
    type <- ifelse(statistic %in% colnames(x@window.stats), "single", "pairwise")
    if(type == "single"){
      scol <- which(colnames(x@window.stats) == statistic)
      matches <- x@window.stats[get.matches(x@window.stats), ..scol]
    }
    else{
      scol <- which(colnames(x@pairwise.window.stats) == statistic)
      matches <- x@pairwise.window.stats[get.matches(x@pairwise.window.stats), ..scol]
    }

    ## get the proportion of the bootstrapped distribution less than or equal to the observed point.
    xp <- edist(matches[[1]])
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

    return(xp)
  }

  #===========loop through all of the requested facets/statistics=========
  browser()

  # generate a matrix containing all of the possible tasks
  snp.

  return(xp)
}

