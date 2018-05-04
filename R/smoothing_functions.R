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
s_ave_multi <- function(x, parms, sigma, ws = NULL, nk = TRUE, levs = c("group", "pop")) {
  sig <- 1000*sigma
  cat("Smoothing Parameters:\n\tsigma = ", sig, "\n\tWindow slide = ", ws*1000, "\n\t")

  #================sanity checks=============
  #nk
  if(!all(is.logical(nk) & length(nk) == 1)){
    stop("nk must be TRUE or FALSE.")
  }
  if(nk){
    if(!any(colnames(x) == "n_total" | colnames(x) == "nk")){
      stop("If nk = TRUE, columns named 'nk' or 'n_total' must be in x.")
    }
    if(any(colnames(x) == "n_total")){
      colnames(x)[which(colnames(x) == "n_total")] <- "nk"
    }
    x[,"nk"] <- as.numeric(x[,"nk"])
    cat("nk weighting = TRUE\n\t")
  }
  else{
    cat("nk weighting = FALSE\n\t")
  }
  x$position <- as.numeric(x$position)


  #parms
  if(!all(parms %in% colnames(x))){
    stop("All entries in parms must be present in x. Column names in x must match entries in parms exactly.")
  }
  cat("Paramters:\n\t")
  for(i in 1:length(parms)){
    cat("\t",parms[i], "\n\t")
  }

  #sigma and ws
  if(!all(is.numeric(sigma), is.numeric(ws), length(sigma) == 1, length(ws) == 1)){
    stop("Sigma and ws must be numeric vectors of length 1.")
  }
  if(sigma >= 500 | ws >= 500){
    warning("Sigma and ws are the number of bases in megabases!")
  }
  else if(sig <= 100){
    warning("Provided sigma is very small:", sig, "bp!")
  }


  #levs
  if(!all(levs %in% colnames(x)) & !all(is.na(levs))){
    stop("All entries in levs must be present in x. Column names in x must match entries in levs exactly.")
  }
  else if(length(levs) > 1 & any(is.na(levs)) | is.numeric(levs) | length(levs) > 1 & any(is.null(levs))){
    stop("All entries in levs must be characters! Use NA to subset by no levels.")
  }
  cat("Smoothing levels:")
  for(i in 1:length(levs)){
    cat("\t",levs[i], "\n\t")
  }

  #================subfunction========
  #funciton to do a sliding window analysis on a data set. Vectorized!:
  sw  <- function(x, scols, ws, sig, nk){
    if(nrow(x) <= 1){
      ret <- rbind(rep(NA, length = 2 + length(parms)))
      ret <- as.data.frame(ret)
      colnames(ret)[1] <- "position"
      colnames(ret)[ncol(ret)] <- "n_snps"
      colnames(ret)[2:(ncol(ret) - 1)] <- paste0("V", 1:(ncol(ret) - 2))
      return(ret)
    }
    #get window centers, starts, and stops:
    ##possible positions
    pos <- x$position
    pos <- as.numeric(pos)

    #window starts, stops, and ends
    if(!is.null(ws)){
      cs <- seq(0, max(pos), by = ws)
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
    n_snps <- n_snps[n_snps != 0]

    #fix any NA values
    vals <- as.matrix(x[,scols])
    NAs <- which(is.na(vals))
    if(length(NAs) > 0){
      vals[NAs] <- 0 #fix the NAs
    }

    #multiply by value of the statistics and nk if requested
    if(nk){
      #set up nks. Have to do it this way because of the possibility of NAs.
      nkv <- matrix(rep(x$nk, ncol(vals)), ncol = ncol(vals))
      nkv <- nkv - 1
      if(length(NAs) > 0){
        nkv[NAs] <- 0 #make sure that these snps don't contribute to the weight!
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

    #prepare metadata
    if(!any(is.na(levs))){
      if(nrow(win_stats) == 1){
        ret <- as.data.frame(cbind(position = colnames(lmat), as.data.frame(win_stats), n_snps = n_snps))
      }
      else{
        levs <- sort(levs)
        meta <- as.character(x[1,levs])
        meta <- rep(meta, length = nrow(win_vals)*length(meta))
        meta <- matrix(meta, nrow = nrow(win_vals), byrow = TRUE)
        colnames(meta) <- levs
        meta <- as.data.frame(meta, stringsAsFactors = F)
        ret <- cbind(meta, position = as.numeric(row.names(win_vals)), as.data.frame(win_stats), n_snps = n_snps)
      }
    }
    else{
      if(nrow(win_stats) == 1){
        ret <- as.data.frame(cbind(position = colnames(lmat), as.data.frame(win_stats), n_snps = n_snps))
      }
      else{
        ret <- cbind(position = as.numeric(row.names(win_vals)), as.data.frame(win_stats), n_snps = n_snps)
      }
    }

    #return
    return(ret)
  }

  #================smooth at proper levels========
  #which cols hold the stats of interest?
  scols <- which(colnames(x) %in% parms)

  if(!is.null(ws)){ws <- ws*1000}

  cat("Smoothing...")

  #do the smoothing, could add a parallel option
  if(!any(is.na(levs))){
    out <- plyr::ddply(
      .data = x,
      .variables = levs,
      .fun = sw,
      .progress = "text",
      scols = scols, ws = ws, sig = sig, nk = nk
    )
  }
  else{
    out <- sw(x, scols, ws, sig, nk)
  }

  #remove any empty facets.
  nas <- ifelse(is.na(out), 1, 0)
  nas <- rowSums(nas)
  nas <- which(nas > 0)
  if(length(nas) > 0){
    out <- out[-nas,]
  }

  #put the correct names on the smoothed variables
  if(length((ncol(out) - length(parms)):(ncol(out) - 1)) == 0){browser()}
  colnames(out)[(ncol(out) - length(parms)):(ncol(out) - 1)] <- paste0("smoothed_", parms)

  #shouldn't ever happen save when this is called with run_gp for backwards compatibility, but may as well include it.
  if(nrow(out) == 0){
    out[1,] <- rep(NA, ncol(out))
  }

  cat("Done!\n")
  return(out)
}
