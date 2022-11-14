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
#'are randomly drawn for each SNP in each window.
#'
#'Centers for windows can either every SNP (if no step size is provided), or
#'every step kilobases from the 0 position of each snp level facet category
#'(chromosome, etc.).
#'
#'The size of sliding windows are defined by the "sigma" argument. Note that
#'this value, as well as that provided to the "step" argument, are given in
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
#'@param verbose Logical, default FALSE. If TRUE, some progress updates will be
#'   printed to the console.
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
calc_smoothed_averages <- function(x, facets = NULL, sigma, step = NULL, nk = TRUE, stats.type = c("single", "pairwise"), 
                                   par = FALSE, triple_sigma = TRUE, gaussian = TRUE,
                                   verbose = FALSE) {
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
  if(!is.null(step)){
    step <- step * 1000
  }
  if(verbose){cat("Smoothing Parameters:\n\twindow size = ", 3*1000*sigma, "\n\tWindow slide = ", step*1000, "\n")}
  
  .sanity_check_window(x, sigma, step, stats.type = stats.type, nk, facets = facets)
  
  
  par <- .par_checker(par, TRUE)
  
  facets <- .check.snpR.facet.request(x, facets, "none")
  snp.facets <- .check.snpR.facet.request(x, facets, "sample")
  sample.facets <- .check.snpR.facet.request(x, facets, "snp")
  x <- .add.facets.snpR.data(x, facets)
  browser()
  
  #============get windows==============
  windows <- .mark_windows(snp.meta(x), sigma = sig, step = step, triple_sig = triple_sigma, chr = unlist(.split.facet(snp.facets)))
  
  cl <- parallel::makePSOCKcluster(par)
  doParallel::registerDoParallel(cl)
  
  browser()
  
  samp_task_list <- .get.task.list(x, sample.facets)
  
  task_window_ids <- sort(rep(1:par, length.out = nrow(windows$win_stats)))
  window_task_list <- split(windows$windows, task_window_ids)

  tout <- foreach::foreach(q = 1:length(window_task_list), .inorder = FALSE,
                           .export = "data.table") %dopar% {
    
  }
  
  
  window_meta <- cbind(data.table::as.data.table(meta[,c(chr, "position")]), effect = p, window = windows)
  window_meta <- window_meta[,as.list(basic(effect)), by = "window"]
  colnames(window_meta)[-1] <- paste0("window_", colnames(window_meta)[-1])
  window_meta <- data.table::melt(window_meta, id.vars = "window")
  window_meta <- window_meta[,as.list(basic(value)), by = "variable"]
  colnames(window_meta)[1] <- "window_stat"
  window_meta <- data.table::melt(window_meta, id.vars = "window_stat")
  window_meta$stat <- paste0(window_meta$window_stat, "_summary_", window_meta$variable)
  window_stats <- window_meta[["value"]]
  names(window_stats) <- window_meta[["stat"]]
}



# Mark windows for snps. More memory efficient but a bit slower, pulled from GeneArchEst. Call if there are a lot of SNPs.
#
# Determine which unique genomic window each snp belongs to.
#
# @param x data.frame. Must contain a "position" column and a chromosome info column.
# @param sigma numeric. Size of windows, in BP.
# @param step numeric, step size.
# @param chr character, default "chr". Name of chromosome info column in x.
#
.mark_windows <- function(x, sigma, step = NULL, chr = "chr", triple_sig = FALSE){
  if(triple_sig){
    sigma <- sigma * 3
  }
  
  if(length(chr) > 1){
    ncc <- .paste.by.facet(x[,chr], chr)
    nccn <- paste0(chr, collapse = ".")
    x[,nccn] <- ncc
    x[,chr] <- NULL
    chr <- nccn
  }
  else if(chr == ".base" | is.null(chr)){
    x$chr <- ".base"
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
  assign_windows <- function(y, starts, ends){
    comp_fun <- function(y, starts, ends){
      lmat <- outer(y, starts, function(pos, starts) pos > starts)
      lmat <- lmat * outer(y, ends, function(pos, ends) pos <= ends)
      lmat[lmat == 1] <- rep(1:nrow(lmat), ncol(lmat))[lmat == 1]
      # if(any(colSums(lmat) == 0)){
      #   warning("Colsums are 0.")
      # }
      
      lmat <- as.data.frame(lmat)
      lmat <- as.list(lmat)
      names(lmat) <- NULL
      lmat <- lapply(lmat, function(z) {z <- z[z != 0]; return(z)})
      return(lmat)
    }
    
    # if large (say, 50k snps), will iterate through in chunks, solve, and then combine results to minimize memory usage
    if(nrow(y) > 50000){
      n_iters <- ceiling(nrow(y)/50000)
      titer <- 1
      lmat <- vector("list", length(starts))
      for(i in 1:n_iters){
        end <- i*50000
        end <- ifelse(end > nrow(y), nrow(y), end)
        trows <- titer:end
        consider_windows <- starts <= max(y$position[trows]) & ends >= min(y$position[trows])
        tpart <- comp_fun(y$position[trows], starts[consider_windows], ends[consider_windows])
        tpart <- lapply(tpart, function(z) z + (titer - 1))
        lmat[consider_windows] <- foreach::foreach(q = 1:length(tpart), .inorder = T) %do% {
          c(lmat[consider_windows][[q]], tpart[[q]])
        }
        titer <- i*50000 + 1
      }
    }
    else{
      lmat <- comp_fun(y$position, starts, ends)
    }
    return(lmat)
  }
  
  
  ## note: need to make sure that the window ids are unique!
  windows <- vector("list", length(unique.chr))
  for(i in 1:length(chr.max)){
    windows[[i]] <- assign_windows(x[x[,chr] == unique.chr[i],], starts[[i]], ends[[i]])
  }
  windows <- unlist(windows, recursive = FALSE)
  empties <- lapply(windows, length) == 0
  windows <- windows[!empties]

  return(.suppress_specific_warning(list(windows = windows, 
                                         win_stats = data.frame(starts = unlist(starts)[!empties], 
                                                                ends = unlist(ends)[!empties], 
                                                                chr = unlist(chrs)[!empties])), "short variable"))
}