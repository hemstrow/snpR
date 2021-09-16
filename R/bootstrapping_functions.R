#'Generate bootstrapped windowed statistics
#'
#'\code{do_bootstraps} creates a distribution of bootstrapped smoothed values
#'for requested statistics.
#'
#'Bootstraps are conducted as described by Hohenlohe et al. (2010). For each
#'bootstrap, this function draws random window position, then draws random
#'statistics from all provided SNPs to fill each observed position on that
#'window and calculates a smoothed statistic for that window using
#'Gaussian-smoothing. Note that a "position" column \emph{must} be present in
#'the snp metadata of the snpRdata object to do any window calculations.
#'
#'Bootstraps for multiple statistics can be calculated at once. If the
#'statistics argument is set to "all", all calculated stats will be run. If it
#'is set to "single", then all non-pairwise statistics will be bootstrapped, if
#'it is set to "pairwise", then all pairwise statistics will be bootstrapped.
#'Individual statistics can also be requested by name ("pi", "ho", etc.). All
#'statistics bootstrapped at the same time will be calculated using \emph{the
#'same randomly filled windows}, and thus do not represent independent
#'observations between statistics. This is still computationally more efficient,
#'however.
#'
#'The data can be broken up categorically by snp or sample metadata, as
#'described in \code{\link{Facets_in_snpR}}.
#'
#'As described in Hohelohe et al. (2010), the contribution of individual SNPs to
#'window averages can be weighted by the number of observations per SNP by
#'setting the nk argument to TRUE, as is the default. For bootstraps, nk values
#'are randomly drawn for each SNP in each window.
#'
#'Possible centers for windows can either SNPs (if no step size is provided), or
#'every step kilobases from the 0 position of each snp level facet category
#'(chromosome, etc.).
#'
#'
#'If do.p is TRUE, calculates p-values for smoothed values of a statistic based 
#'upon the bootstrapped null distribution of that statistic using an empirical 
#'continuous distribution function. Note that in this case, the minimum possible
#'p-value for a window depends upon the number of bootstraps calculated (if only
#'a 1000 bootstraps were performed, the minimum possible p-value is about .001, 
#'or one in a thousand.)
#'
#'Note that this function will return an error if equivalent windowed statistics
#'have not first been calculated for the designated facets (if the "chromosome"
#'facet is requested with a sigma of 200 and a step of 50, do_bootstraps will
#'error unless \code{\link{calc_smoothed_averages}} has not yet been run with
#'the same facet, sigma, and step values).
#'
#'@param x snpRdata object.
#'@param facets character or NULL, default NULL. Categories by which to break up
#'  bootstraps.
#'@param boots numeric. Number of bootstraps to generate.
#'@param sigma numeric. Designates the width of windows in kilobases. Full
#'  window size is 6*sigma.
#'@param step numeric or NULL, default NULL. Designates the number of kilobases
#'  between each window centroid. If NULL, windows are centered on each SNP.
#'@param statistics character. Designates the statistic(s) to smooth, typically
#'  "all", "single", or "pairwise". See details for options.
#'@param nk logical, default TRUE. If TRUE, weights SNP contribution to window
#'  averages by the number of observations at those SNPs.
#'@param par numeric or FALSE, default FALSE. If numeric, the number of cores to
#'  use for parallel processing.
#'@param do.p logical, default TRUE. Determines if p-values should be calculated
#'  for sliding windows.
#'@param p.alt character, default "two-sided". Specifies the alternative
#'  hypothesis to be used. Options: \itemize{ \item "less": probability that a
#'  bootstrapped value is as small or smaller than observed. \item "greater":
#'  probability that a bootstrapped value is as large or larger than observed.
#'  \item "two-sided": probability that a bootstrapped value is as or more
#'  extreme than observed. }
#'@return A snpRdata object with bootstrapped windows merged in to the
#'  window.bootstraps slot. If do.p is TRUE, it will also merge p values in for
#'  bootstrapped statistics into the stats or window.stats sockets.
#'
#'@references Hohenlohe et al. (2010). \emph{PLOS Genetics}

#'@author William Hemstrom
#'@export
#'
#' @examples
#' # add statistics
#' dat <- calc_basic_snp_stats(stickSNPs, "group", sigma = 200, step = 150)
#' 
#' # do bootstraps
#' dat <- do_bootstraps(dat, facets = "group", boots = 1000, 
#'                      sigma = 200, step = 150)
#' 
#' # fetch results, bootstraps and then p-values on original stats
#' get.snpR.stats(dat, "group", "bootstraps")
#' get.snpR.stats(dat, "group", "single.window")
#' 
do_bootstraps <- function(x, facets = NULL, boots, sigma, step = NULL, statistics = "all", nk = TRUE, par = FALSE, do.p = TRUE, p.alt = "two-sided"){
  #note: it is possible to run all sample level facets at once, so something like c("pop.fam.group", "pop.group") can
  #      be run simultainously, with no need to loop across facets.
  #      However, SNP level facets create different windows, and so need to be run seperately. Essentially,
  #      we need to do everything once for each unique snp level facet defined in the data.

  #================sanity checks================
  msg <- character(0)
  o.facets <- facets
  facets <- check.snpR.facet.request(x, facets, "none")

  # figure out which stats we are using!
  pairwise.types <- c("fst")
  single.types <- c("pi", "ho", "pa", "pHWE")
  all.types <- c(pairwise.types, single.types)
  o.stats <- statistics
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
    stats.type <- c(stats.type, "single")
  }


  # get mafs if doing any normal stats. Can make this more efficient by adding only where missing in the future.
  if("single" %in% stats.type & nk){
    what.stats <- check_calced_stats(x, check.snpR.facet.request(x, facets), "maf")
    if(any(!unlist(what.stats))){
      x <- calc_maf(x,  check.snpR.facet.request(x, facets)[which(!unlist(what.stats))])
    }
    x@stats$nk <- x@stats$maj.count + x@stats$min.count
  }
  
  # check that we've actually calculated windowed stats for the facet we are doing!
  facet_details <- check.snpR.facet.request(x, facets, "none", "TRUE")
  for(i in 1:length(facets)){
    split.facet <- check.snpR.facet.request(x, unlist(.split.facet(facets[i])), "none", TRUE)
    snp.part <- which(split.facet[[2]] == "snp")
    snp.part <- split.facet[[1]][snp.part]
    if(length(snp.part) == 0){
      snp.part <- ".base"
    }
    samp.part <- which(split.facet[[2]] == "sample")
    samp.part <- split.facet[[1]][samp.part]
    if(length(samp.part) == 0){
      samp.part <- ".base"
    }
    
    
    
    # check pairwise
    if("pairwise" %in% stats.type){
      has.samp.part.pw <- which(x@pairwise.window.stats$facet == paste0(samp.part, collapse = "."))
      has.snp.part.pw <- which(x@pairwise.window.stats$snp.facet == paste0(snp.part, collapse = "."))
      has.step.matches.pw <- which(x@pairwise.window.stats$step == step)
      has.sigma.matches.pw <- which(x@pairwise.window.stats$sigma == sigma)
      has.both.parts.pw <- intersect(intersect(has.snp.part.pw, has.samp.part.pw), intersect(has.step.matches.pw, has.sigma.matches.pw))
      
    }
    if("single" %in% stats.type){
      has.samp.part.s <- which(x@window.stats$facet == paste0(samp.part, collapse = "."))
      has.snp.part.s <- which(x@window.stats$snp.facet == paste0(snp.part, collapse = "."))
      has.step.matches.s <- which(x@window.stats$step == step)
      has.sigma.matches.s <- which(x@window.stats$sigma == sigma)
      has.both.parts.s <- intersect(intersect(has.snp.part.s, has.samp.part.s), intersect(has.sigma.matches.s, has.step.matches.s))
    }
    
    # check single
    if("pairwise" %in% stats.type){
      if(length(has.both.parts.pw) == 0){
        msg <- c(msg, paste0("Windowed stats with the same step and sigma not calculated for any pairwise statstistics for facet: ", facets[i], "."))
      }
    }
    if("single" %in% stats.type){
      if(length(has.both.parts.s) == 0){
        msg <- c(msg, paste0("Windowed stats with the same step and sigma not calculated for any single statstistics for facet: ", facets[i], "."))
      }
    }
  }
  
  if(length(msg) > 0){
    stop(paste0(msg, "\n"))
  }
  
  # run basic sanity checks
  sanity_check_window(x, sigma, 200, stats.type, nk, facets, statistics, good.types = all.types)
  

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
      if(!is.null(part.cols$single) & !is.null(part.cols$pairwise)){
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

      fws <- lapply(fws$lmat, stats::na.omit)

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
      fws <- lapply(fws, stats::na.omit)


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
    if("single" %in% stats.type){
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
        snk <- data.table::dcast(stats[,nk.meta, with = FALSE], .snp.id~facet + subfacet, value.var = "nk")
        snk <- snk[nrands,-".snp.id"]
        ## need to duplicate columns, since each nk value will be used for multiple stats!
        col.seq <- rep(1:ncol(snk), sum(statistics %in% single.types))
        snk <- snk[,col.seq, with = FALSE]
      }
      stats <- stats[,keep.cols, with = FALSE]
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
        pnk <- data.table::dcast(pairwise.stats[,nk.meta, with = FALSE], .snp.id~facet + comparison, value.var = "nk")
        pnk <- pnk[nrands,-".snp.id"]
        col.seq <- rep(1:ncol(pnk), sum(statistics %in% pairwise.types))
        pnk <- pnk[,col.seq, with = FALSE]
      }
      pairwise.stats <- pairwise.stats[,keep.cols, with = FALSE]
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
      part.cols$single <- ncol(stats)
      part.cols$pairwise <- ncol(pairwise.stats)
    }
    else if(any(statistics %in% single.types)){
      bound.stats <- stats
      part.cols$single <- ncol(stats)
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
    boot.set <- boot.set[,ord, with = FALSE]

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
  
  # remove empty facets (sample level only facets)
  snp.facet.matches <- Filter(Negate(is.null), snp.facet.matches)
  
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
  out <- out[,col.ord, with = FALSE]
  x@window.bootstraps <- rbind(x@window.bootstraps, out)

  if(do.p){
    x <- calc_p_from_bootstraps(x, facets, o.stats, alt = p.alt, par = par)
  }

  return(x)
}

#calculates a p-value for a stat based on a null distribution caluclated via bootstraps of that stat.
#inputs:  x: vector of observed statistics
#         dist: vector of statistics bootstrapped from same stat.
#         alt: Alternative hypothesis for p-value. Can be less, greater, or two-sided (default).


#'Caculate p-values from bootstrapped distributions.
#'
#'\code{calc_p_from_bootstraps} finds p-values for observed \emph{smoothed
#'window} statistics from bootstrapped distributions, such as produced by
#'\code{\link{do_bootstraps}}.
#'
#'Calculates p-values for smoothed values of a statistic based upon a
#'bootstrapped null distribution of that statistic using an emperical continuous
#'distribution function.
#'
#'p-values can be generated for specific snp or sample metadata categories
#'susing the facets argument, as described in \code{\link{Facets_in_snpR}}. Only
#'facets for which bootstrap data and raws statistical data have both been
#'calculated will be run. "all" and NULL follow the typical facet rules.
#'
#'Likewise, p-values can be generated for specific statistics using the
#'statistics argument. Only statistics for which bootstrap data and raws
#'statistical data have both been calculated will be run. By default, all stats
#'for which a bootstrap null distribution has been generated will be run.
#'
#'@param x snpRdata object.
#'@param facets character, default "all". Facets to use.
#'@param statistics character, default "all". Vector naming the statistics to
#'  calculate p-values for. By default calculates p-values for all possible
#'  stats.
#'@param alt character, default "two-sided". Specifies the alternative
#'  hypothesis to be used. Options: \itemize{ \item "less": probability that a
#'  bootstrapped value is as small or smaller than observed. \item "greater":
#'  probability that a bootstrapped value is as large or larger than observed.
#'  \item "two-sided": probability that a bootstrapped value is as or more
#'  extreme than observed. }
#'@param par numeric or FALSE, default FALSE. If numeric, the number of cores to
#'  use for parallel processing.
#'@param fwe_method character, default c("bonferroni", "holm", "BH", "BY"). Type
#'  of Family-Wise Error correction (mulitple testing correction) to use. For
#'  details and options, see \code{\link{p.adjust}}.
#'@param fwe_case character, default c("by_facet", "by_subfacet", "overall").
#'  How should Family-Wise Error correction (multiple testing correction) be
#'  applied? \itemize{\item{"by_facet":} Each facet supplied (such as pop or
#'  pop.fam) is treated as a set of tests. \item{"by_subfacet":} Each level of
#'  each subfacet is treated as a seperate set of tests. \item{"overall":} All
#'  tests are treated as a set.}
#'
#'
#'@return snpRdata object, with p-values merged into the stats or pairwise.stats
#'  sockets.
#'
#'@seealso ecdf
#'
#'@export
#'@author William Hemstrom
#'
#' @examples
#' # add statistics and generate bootstraps
#' x <- calc_basic_snp_stats(stickSNPs, c("group.pop"), sigma = 200, step = 150)
#' x <- do_bootstraps(x, facets = c("group.pop"), boots = 1000, sigma = 200, step = 150)
#' x <- calc_p_from_bootstraps(x)
#' get.snpR.stats(x, "group.pop", "single.window") # pi, ho, etc
#' get.snpR.stats(x, "group.pop", "pairwise.window") # fst
#'
calc_p_from_bootstraps <- function(x, facets = "all", statistics = "all", alt = "two-sided", par = FALSE,
                                   fwe_method = "BY", 
                                   fwe_case = c("by_facet", "overall")){
  #==========sanity checks==========
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  msg <- character()

  # check statistics
  if(statistics != "all"){
    miss.stats <- which(!statistics %in% unique(x@window.bootstraps$stat))
    if(length(miss.stats) > 0){
      msg <- c(msg,
               paste0("Some statistics not found in bootstraps: ", paste0(miss.stats, collapse = ", "), "\n"))
    }
  }

  if(length(msg) > 0){
    stop(msg)
  }

  #===========subfunctions=======
  # finds matching rows for given facet, ect
  get.matches <- function(y, facet, subfacet, snp.facet, sigma, nk, step, statistic){
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
    matches <- get.matches(x@window.bootstraps, facet, subfacet, snp.facet, sigma, nk, step, statistic)


    # create a cumulative distibution function of the distribution
    edist <- stats::ecdf(x@window.bootstraps[matches,]$value)

    # get p-values
    ## get the matching statistic data
    type <- ifelse(statistic %in% colnames(x@window.stats), "single", "pairwise")
    if(type == "single"){
      scol <- which(colnames(x@window.stats) == statistic)
      meta.cols <- 1:which(colnames(x@window.stats) == "nk.status")
      matches <- get.matches(x@window.stats, facet, subfacet, snp.facet, sigma, nk, step, statistic)
      meta <- x@window.stats[matches, meta.cols, with = FALSE]
      matches <- x@window.stats[matches, scol, with = FALSE]
    }
    else{
      scol <- which(colnames(x@pairwise.window.stats) == statistic)
      meta.cols <- 1:which(colnames(x@window.stats) == "nk.status")
      matches <- get.matches(x@pairwise.window.stats, facet, subfacet, snp.facet, sigma, nk, step, statistic)
      meta <- x@pairwise.window.stats[matches, meta.cols, with = FALSE]
      matches <- x@pairwise.window.stats[matches,  scol, with = FALSE]
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

    # bind and return
    out <- cbind(meta, stat = statistic, p = xp)
    return(out)
  }

  #===========loop through all of the requested facets/statistics=========
  # generate a matrix containing all of the possible tasks
  u.rows <- unique(x@window.bootstraps[,1:(which(colnames(x@window.bootstraps) == "value") - 1)])

  # figure out which parts of the unique tasks are part of the requested facet
  if(facets[1] != "all"){

    keep.rows <- logical(nrow(u.rows))
    facets <- check.snpR.facet.request(x, facets, "none", T)

    # loop through each facet
    for(i in 1:length(facets[[1]])){

      # for sample facets
      if(facets[[2]][[i]] == "sample"){
        keep.rows[which(u.rows$facet == facets[[1]][[i]])] <- T
      }
      # for complex facets
      else if(facets[[2]][[i]] == "complex"){
        u.facet <-  unlist(.split.facet(facets[[1]][[i]]))
        split.parts <- check.snpR.facet.request(x, u.facet, remove.type = "none", return.type = T)
        samp.part <- split.parts[[1]][which(split.parts[[2]] == "sample")]
        snp.part <- split.parts[[1]][which(split.parts[[2]] == "snp")]

        keep.samp <- which(u.rows$facet == samp.part)
        keep.snp <- which(u.rows$snp.facet == snp.part)
        keep.rows[intersect(keep.samp, keep.snp)] <- T
      }
      else if(facets[[2]][[i]] == "snp"){
        keep.rows[which(u.rows$snp.facet == facets[[1]][[i]])] <- T
      }
      else{
        keep.rows[which(u.rows$facet == ".base" & u.rows$snp.facet == ".base")] <- T
      }
    }
    u.rows <- u.rows[keep.rows,]
  }

  # initialize
  out <- vector("list", nrow(u.rows))

  # run the loop
  if(par == FALSE){
    for(q in 1:nrow(u.rows)){
      out[[q]] <- get.one.pvalue(x, facet = u.rows$facet[q],
                                 subfacet = u.rows$subfacet[q],
                                 snp.facet = u.rows$snp.facet[q],
                                 statistic = u.rows$stat[q],
                                 nk = u.rows$nk[q],
                                 step = u.rows$step[q],
                                 sigma = u.rows$sigma[q],
                                 alt = alt)
    }
  }
  else{
    cl <- snow::makeSOCKcluster(par)
    doSNOW::registerDoSNOW(cl)

    #prepare reporting function
    ntasks <- nrow(u.rows)
    progress <- function(n) cat(sprintf("\tPart %d out of", n), ntasks, "is complete.\n")
    opts <- list(progress=progress)

    out <- foreach::foreach(q = 1:ntasks, .inorder = FALSE,
                            .options.snow = opts, .export = "data.table") %dopar% {
                              get.one.pvalue(x, facet = u.rows$facet[q],
                                             subfacet = u.rows$subfacet[q],
                                             snp.facet = u.rows$snp.facet[q],
                                             statistic = u.rows$stat[q],
                                             nk = u.rows$nk[q],
                                             step = u.rows$step[q],
                                             sigma = u.rows$sigma[q],
                                             alt = alt)
                            }

    #release cores
    parallel::stopCluster(cl)
    doSNOW::registerDoSNOW()
  }

  #===========bind, cast, and merge===========
  # bind the result
  out <- data.table::rbindlist(out)

  # melt
  meta.cols <- colnames(out)[-ncol(out)]
  meta.cols <- meta.cols[-which(meta.cols == "stat")]

  # merge
  if(any(out$stat == "fst")){
    cout <- data.table::dcast(out[which(out$stat == "fst"),], facet + subfacet + snp.facet + snp.subfacet + position + sigma + n_snps + step + nk.status ~ stat, value.var = "p")
    colnames(cout)[(which(colnames(cout) == "nk.status") + 1):ncol(cout)] <- paste0("p_", colnames(cout)[(which(colnames(cout) == "nk.status") + 1):ncol(cout)])
    
    # do fwe
    cout <- fwe_correction(cout, levs = c("facet", "subfacet"), pcol = "p_fst", methods = fwe_method, case = fwe_case)
    x <- merge.snpR.stats(x, cout, type = "pairwise.window.stats")
  }
  if(any(out$stat != "fst")){

    s.stats <- which(out$stat != "fst")
    cout <- data.table::dcast(out[s.stats,], facet + subfacet + snp.facet + snp.subfacet + position + sigma + n_snps + step + nk.status ~ stat, value.var = "p")
    colnames(cout)[(which(colnames(cout) == "nk.status") + 1):ncol(cout)] <- paste0("p_", colnames(cout)[(which(colnames(cout) == "nk.status") + 1):ncol(cout)])
    
    # do fwe
    p_cols <- colnames(cout)[grep("p_", colnames(cout))]
    cout_comb <- vector("list", length(p_cols))
    for(i in 1:length(p_cols)){
      cout_comb[[i]] <- fwe_correction(cout, levs = c("facet", "subfacet"), pcol = p_cols[i], methods = fwe_method, case = fwe_case)
    }
    invisible(suppressMessages(cout <- purrr::reduce(cout_comb, dplyr::left_join)))
    
    x <- merge.snpR.stats(x, cout, type = "window.stats")
  }


  return(x)
}