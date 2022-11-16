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

#'Gaussian smooth or average statistics across sliding windows.
#'
#'Calculates optionally gaussian smoothed average values for statistics and the
#'genomic position at which it was observed, splitting by any requested
#'variables.
#'
#'Averages for multiple statistics can be calculated at once. If the statistics
#'argument is set to c("pairwise", "single"), all calculated stats will be run.
#'If it is set to "single", then all non-pairwise statistics (pi, ho, maf, etc.)
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
#'are randomly drawn for each SNP in each window. SNPs can also optionally be
#'weighted by their proximity to the window centroid using a gaussian kernal.
#'
#'Centers for windows can either every SNP (if no step size is provided), or
#'every step kilobases from the 0 position of each snp level facet category
#'(chromosome, etc.).
#'
#'The size of sliding windows are defined by the "sigma" argument. Note that
#'this value, as well as that provided to the "step" argument, are given in
#'kilobases. Each window will include SNPs within 3*sigma kilobases from the
#'window center (if the \code{triple_sigma} option is selected). Past this
#'point, the effect of each additional SNP on the window average would be very
#'small, and so they are dropped for computational efficiency (see Hohenlohe
#'(2010)).
#'
#'@param x snpRdata object.
#'@param facets character or NULL, default NULL. Categories by which to break up
#'  windows.
#'@param sigma numeric. Designates the width of windows in kilobases. Full
#'  window size is 6*sigma or sigma kilobases depending on the
#'  \code{triple_sigma} argument.
#'@param step numeric or NULL, default NULL. Designates the number of kilobases
#'  between each window centroid. If NULL, windows are centered on each SNP.
#'@param nk logical, default TRUE. If TRUE, weights SNP contribution to window
#'  averages by the number of observations at those SNPs.
#'@param stats.type character, default c("single", "pairwise"). Designates the
#'  statistic(s) to smooth, either "single",  "pairwise", or c("single",
#'  "pairwise"). See details.
#'@param par numeric or FALSE, default FALSE. If numeric, the number of cores to
#'  use for parallel processing.
#'@param triple_sigma Logical, default TRUE. If TRUE, sigma will be tripled to
#'  create windows of 6*sigma total.
#'@param gaussian Logical, default TRUE. If TRUE, windows will be gaussian
#'  smoothed. If not, windows will be raw averages.
#'@param verbose Logical, default FALSE. If TRUE, some progress updates will be
#'  printed to the console.
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
calc_smoothed_averages <- function(x, facets = NULL, sigma, step = NULL, nk = TRUE, 
                                   stats.type = c("single", "pairwise"), 
                                   par = FALSE, triple_sigma = TRUE, gaussian = TRUE,
                                   verbose = FALSE) {
  .snp.id <- chr <- position <- start <- end <- ..col_ord <- NULL
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
  
  sigma <- 1000*sigma
  if(!is.null(step)){
    step <- step * 1000
  }
  if(verbose){cat("Smoothing Parameters:\n\twindow size = ", ifelse(triple_sigma, 6*sigma, 2*sigma), "\n\tWindow slide = ", ifelse(is.null(step), "per-snp", step*1000), "\n")}
  
  .sanity_check_window(x, sigma/1000, step/1000, stats.type = stats.type, nk, facets = facets)
  
  
  par <- .par_checker(par, TRUE)
  
  facets <- .check.snpR.facet.request(x, facets, "none")
  snp.facets <- .check.snpR.facet.request(x, facets, "sample")
  sample.facets <- .check.snpR.facet.request(x, facets, "snp")
  x <- .add.facets.snpR.data(x, facets)

  #============run smoothing==============

  if("single" %in% stats.type){
    stats <- get.snpR.stats(x, facets = sample.facets, "single")
    stats$nk <- stats$maj.count + stats$min.count
    stats$maj.count <- NULL
    stats$min.count <- NULL
    numeric.cols <- which(sapply(stats, class) %in% c("numeric", "integer") & colnames(stats) != "nk")
    numeric.cols <- intersect(numeric.cols, (which(colnames(stats) == ".snp.id") + 1):ncol(stats))
    
    
    task_list <- .get.task.list(x, facets)
    
    if(verbose){cat("Beginning run: single stats.\n")}
    
    cl <- parallel::makePSOCKcluster(min(c(par, nrow(task_list))))
    doParallel::registerDoParallel(cl)

    out <- foreach::foreach(q = 1:nrow(task_list), .inorder = FALSE,
                            .export = c("data.table", ".average_windows", ".fetch.snp.meta.matching.task.list", ".paste.by.facet", ".split.facets")) %dopar% {
                              
                              snp.matches <- .fetch.snp.meta.matching.task.list(x = x, task_list[q,])
                              snp.matches <- snp.meta(x)$.snp.id[snp.matches]
                              
                              sample.matches <- which(.paste.by.facet(stats, colnames(stats)[1:2]) == paste0(task_list[q,1], ".", task_list[q,2]))
                              matches <- intersect(sample.matches, which(stats$.snp.id %in% snp.matches))
                              if(length(matches) != 0){
                                out <- .average_windows(stats[matches,c(unlist(.split.facet(snp.facets)), "position")], 
                                                        chr =  unlist(.split.facet(snp.facets)),
                                                        sigma = sigma,
                                                        step = step, 
                                                        triple_sig = triple_sigma, 
                                                        gaussian = gaussian, 
                                                        stats = data.table::as.data.table(stats[matches,c(numeric.cols, which(colnames(stats) == "nk"))]),
                                                        nk = nk)
                                out$subfacet <- task_list[q,2]
                                out$facet <- task_list[q,1]
                                out$snp.facet <- task_list[q,3]
                                out$snp.subfacet <- task_list[q,4]
                              }
                              else{
                                out <- NULL
                              }
                              out
                            }
    
    
    parallel::stopCluster(cl)
    
    
    out <- dplyr::bind_rows(out)
    out <- dplyr::arrange(out, snp.facets, position, start, end)
    out$nk.status <- nk
    out$gaussian <- gaussian
    out$sigma <- sigma/1000
    out$step <- ifelse(is.null(step), NA, step/1000)
    out$triple_sigma <- triple_sigma
    
    header_cols <- c("facet", "subfacet", "snp.facet", "snp.subfacet", "position",
                     "sigma", "step", "nk.status", "gaussian", "n_snps", "triple_sigma")
    col_ord <- c(header_cols, colnames(out)[which(!colnames(out) %in% header_cols)])
    .fix..call(out <- out[,..col_ord])
    x <- .merge.snpR.stats(x, out, "window.stats")
    
  }
  
  
  if("pairwise" %in% stats.type){
    stats <- get.snpR.stats(x, facets = sample.facets, "pairwise")
    numeric.cols <- which(sapply(stats, class) %in% c("numeric", "integer") & colnames(stats) != "nk")
    numeric.cols <- intersect(numeric.cols, (which(colnames(stats) == ".snp.id") + 1):ncol(stats))
    
    
    task_list <- .get.task.list(x, snp.facets)
    task_list <- task_list[,3:4,drop = FALSE]
    task_list <- as.data.frame(task_list)
    
    sample_task_list <- unique(as.character(stats$comparison))
    sample_task_list <- cbind(sample.facets, sample_task_list)
    sample_task_list <- as.data.frame(sample_task_list)
    
    task_list <- merge(sample_task_list, task_list)
    
    cl <- parallel::makePSOCKcluster(min(c(par, nrow(task_list))))
    
    if(verbose){cat("Beginning run: pairwise stats.\n")}
    doParallel::registerDoParallel(cl)

    out <- foreach::foreach(q = 1:nrow(task_list), .inorder = FALSE,
                            .export = c("data.table", ".average_windows", ".fetch.snp.meta.matching.task.list", ".paste.by.facet", ".split.facets")) %dopar% {
                              
                              snp.matches <- .fetch.snp.meta.matching.task.list(x = x, task_list[q,])
                              snp.matches <- snp.meta(x)$.snp.id[snp.matches]
                              
                              sample.matches <- which(.paste.by.facet(stats, colnames(stats)[1:2]) == paste0(task_list[q,1], ".", task_list[q,2]))
                              matches <- intersect(sample.matches, which(stats$.snp.id %in% snp.matches))
                              
                              if(length(matches) != 0){
                                out <- .average_windows(stats[matches,c(unlist(.split.facet(snp.facets)), "position")], 
                                                        chr =  unlist(.split.facet(snp.facets)),
                                                        sigma = sigma,
                                                        step = step, 
                                                        triple_sig = triple_sigma, 
                                                        gaussian = gaussian, 
                                                        stats = data.table::as.data.table(stats[matches,c(numeric.cols, which(colnames(stats) == "nk"))]),
                                                        nk = nk)
                                out$subfacet <- task_list[q,2]
                                out$facet <- task_list[q,1]
                                out$snp.facet <- task_list[q,3]
                                out$snp.subfacet <- task_list[q,4]
                              }
                              else{
                                out <- NULL
                              }
                              out
                            }
    
    
    parallel::stopCluster(cl)
    
    
    out <- dplyr::bind_rows(out)
    out <- dplyr::arrange(out, snp.facets, position, start, end)
    out$nk.status <- nk
    out$gaussian <- gaussian
    out$sigma <- sigma/1000
    out$step <- ifelse(is.null(step), NA, step/1000)
    out$triple_sigma <- triple_sigma
    
    header_cols <- c("facet", "subfacet", "snp.facet", "snp.subfacet", "position",
                     "sigma", "step", "nk.status", "gaussian", "n_snps", "triple_sigma")
    col_ord <- c(header_cols, colnames(out)[which(!colnames(out) %in% header_cols)])
    .fix..call(out <- out[,..col_ord])
    x <- .merge.snpR.stats(x, out, "pairwise.window.stats")
  }
  
  
  return(x)
}



# Average stats for windows for snps. 
#
# More memory efficient but a bit slower if there are many SNPs, parts pulled from GeneArchEst.
#
# @param x data.frame. Must contain a "position" column and may contain chromosome info columns
# @param sigma numeric. Size of windows on each side of center, in BP.
# @param step numeric, step size, in BP.
# @param chr character, default "chr". Name of chromosome info column(s) in x. Should be passed through .check.snpR.facet.request first.
# @param triple_sig logical, default FALSE. If true, triples sigma (original protocol)
# @param stats data.table, stats to smooth. any columns not matching position or chr will be smoothed and must be numeric.
# @param gaussian logical, default FALSE. If true, does gaussian smoothing. If not, does averages.
# @param nk logical, default TRUE. If TRUE, weights by nk - 1. Expects 'nk' column in stats.
#
.average_windows <- function(x, sigma, step = NULL, chr = "chr", triple_sig = FALSE,
                          stats, gaussian = FALSE, nk = TRUE){
  ..scols <- ..chr <- NULL
  
  if(triple_sig){
    sigma <- sigma * 3
  }
  
  if(length(chr) > 1){
    nccn <- paste0(chr, collapse = ".")
    
    
    ncc <- .paste.by.facet(x[,chr], chr)
    x[,nccn] <- ncc
    x[,chr] <- NULL
    
    chr <- nccn
  }
  else if(chr == ".base" | is.null(chr)){
    x$chr <- ".base"
    stats$chr <- ".base"
    chr <- "chr"
  }
  
  # window start and end points
  unique.chr <- sort(unique(x[,chr]))
  chr.max <- tapply(x$position, x[,chr], max)
  chr.min <- tapply(x$position, x[,chr], min)
  chr.range <- matrix(c(chr.max, chr.min), ncol = 2)
  
  if(!is.null(step)){
    centers <- apply(chr.range, 1, function(y) seq(from = 0, to = y[1], by = step))
  }
  else{
    centers <- lapply(unique.chr, function(y) unique(x[x[,chr] == y,]$position))
  }
  
  if(is.matrix(centers)){
    centers <- as.list(as.data.frame(centers))
  }
  
  ends <- foreach::foreach(q = 1:length(chr.max), .inorder = T) %do% {
    e <- centers[[q]] + sigma
    e[e > chr.max[q]] <- chr.max[q]
    e
  }
  starts <- foreach::foreach(q = 1:length(chr.max), .inorder = T) %do% {
    s <- centers[[q]] - sigma
    s[s < 0] <- 0
    s
  }
  chrs <- foreach::foreach(q = 1:length(chr.max), .inorder = T) %do% {
    rep(unique.chr[q], length(centers[[q]]))
  }
  
  

  
  # for each chr, assign windows to all snps
  ## function per chr
  assign_windows <- function(y, starts, ends, centers, stats, chrs){
    comp_fun <- function(y, starts, ends, centers, stats, scols, nk_dt){
      # find snps in windows
      lmat <- outer(y$position, starts, function(pos, starts) pos > starts)
      lmat <- lmat * outer(y$position, ends, function(pos, ends) pos <= ends)
      
      tnsnps <- colSums(lmat)

      # weight (gaussian)
      if(gaussian){
        lmat <-  outer(y$pos, centers, gaussian_weight, s = sigma)*lmat
      }
      
      lmat[is.na(lmat)] <- 0
      

      # calculate weighted stats
      if(nk){
        build <- t(lmat) %*% as.matrix(.fix..call(stats[,..scols])*(nk_dt - 1))
        build_weights <- t(lmat) %*% as.matrix(nk_dt - 1)
      }
      else{
        build <- t(lmat) %*% as.matrix(.fix..call(stats[,..scols]))
        build_weights <- matrix(colSums(lmat), nrow = ncol(lmat), ncol = length(scols))
      }
      
      return(list(build = build, build_weights = build_weights, n_snps = tnsnps))
    }
    
    scols <- which(!colnames(stats) %in% c(chr, "position", "nk"))
    
    # check for NAs in either stats or weights (don't want them to contribute)
    nk_dt <- data.table::as.data.table(matrix(stats$nk, nrow = nrow(stats), ncol = length(scols)))
    if(nk){
      nk_NAs <- which(is.na(stats$nk))
      if(length(nk_NAs) > 0){
        data.table::set(stats, i = as.integer(nk_NAs), j = scols, value = 0)
        stats$nk[nk_NAs] <- 0
        data.table::set(nk_dt, nk_NAs, 1:ncol(nk_dt), value = 0)
      }
    }
    for(k in 1:length(scols)){
      if(nk){
        data.table::set(nk_dt, which(is.na(stats[[scols[k]]])), k, value = 0)
      }
      data.table::set(stats,which(is.na(stats[[scols[k]]])),scols[k],0)
    }
    
    if(nk){
      for(k in 1:ncol(nk_dt)){
        data.table::set(nk_dt, which(nk_dt[[k]] == 0), k, value = 1) # set missing data to 1 to prevent negative weights (since we do nk - 1)
      }
    }
    
    
    
    # if large (say, 50k snps), will iterate through in chunks, solve, and then combine results to minimize memory usage
    max_snps <- 50000
    if(nrow(y) > max_snps){
      n_iters <- ceiling(nrow(y)/max_snps)
      titer <- 1
      
      # initialize storage: build will hold the building weighted values for each stat for each window, build_weights will hold the total weights for those windows
      build <- matrix(0, length(starts), length(scols))
      build <- data.table::as.data.table(build)
      colnames(build) <- colnames(stats)[scols]
      build <- list(build = build, build_weights = build, n_snps = numeric(length(starts)))

      # for each set of snps, get the weighted values and weights for all relevent windows.
      for(i in 1:n_iters){
        end <- i*max_snps
        end <- ifelse(end > nrow(y), nrow(y), end)
        trows <- titer:end
        consider_windows <- starts <= max(y$position[trows]) & ends >= min(y$position[trows])
        bpart <- comp_fun(y[trows,], 
                          starts[consider_windows], 
                          ends[consider_windows], 
                          centers[consider_windows], 
                          stats[trows,], 
                          scols,
                          nk_dt[trows,])
        titer <- i*max_snps+ 1
        
        build$build[consider_windows,] <- build$build[consider_windows,] + bpart$build
        build$build_weights[consider_windows,] <- build$build_weights[consider_windows,] + bpart$build_weights
        build$n_snps <- build$n_snps[consider_windows] + bpart$n_snps
      }
    }
    else{
      build <- comp_fun(y, starts, ends, centers, stats, scols, nk_dt)
      build[1:2] <- lapply(build[1:2], data.table::as.data.table)
    }
    
    
    # calculate weighted values
    out <- build$build / build$build_weights
    
    # remove empty windows
    empties <- which(rowSums(dplyr::mutate_all(out, is.nan)) == ncol(out))
        
    out[,chr] <- chrs
    out$position <- centers
    out$n_snps <- build$n_snps
    
    
    if(length(empties) > 0){
      out <- out[-empties,]
    }
    return(out)
  }
  
  
  ## note: need to make sure that the window ids are unique!
  windows <- vector("list", length(unique.chr))
  for(i in 1:length(chr.max)){
    windows[[i]] <- assign_windows(x[x[,chr] == unique.chr[i],], starts[[i]], ends[[i]], centers[[i]], 
                                   stats = stats[x[,chr] == unique.chr[i],], chrs = chrs[[i]])
  }
  
  windows <- data.table::rbindlist(windows)
  
  windows[,chr] <- NULL

  return(windows)
}