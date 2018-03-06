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
#'\code{smoothed_average} Calculates gaussian smoothed average values for a given statistic and the genomic position at which it was observed.
#'
#'Oblivious to group/chromosome and population, so should probably be typically wrapped in run_g or run_gp. Can weight each statistic by an additional factor (typically sample size, "nk", per SNP). Smoothed values are calculated for sliding windows which can have either arbitrary centers (fixed_windows) or centers defined by each observed SNP.
#'
#'Description of x:
#'    Requires columns titled columns titled "position" and one titled to match the "parameter" argument, containing positions in bp and statistic to be smoothed for each SNP. If nk_weight is set, a column titled "nk" is required containing the weighting factor for smoothing, typically the number of observed alleles/sample size, "nk".
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
  sig <- 1000*sigma
  new_col <- paste0("smoothed", "_", parameter)
  if(nk_weight){
    if(any(colnames(x) == "n_total")){
      colnames(x)[which(colnames(x) == "n_total")] <- "nk"
    }
    x[,"nk"] <- as.numeric(x[,"nk"])
  }
  x$position <- as.numeric(x$position)
  if (is.null(fixed_window)){
    out <- data.frame(position = x$position, fill = numeric(nrow(x)), n_snps = numeric(nrow(x)))
    colnames(out)[2] <- new_col
    for (i in 1:nrow(x)) {
      if(i %% 1000 == 0){cat("snp ", i, "\n")}
      c <- x[i,"position"]
      start <- c - 3*sig
      end <- c + 3*sig
      #cat("c:", c, "sig:", sig, "\n")
      ws <- x[x$position <= end & x$position >= start,]
      ws <- ws[!is.na(ws[,parameter]),] #remove na values

      if(!is.null(nk_weight)){
        gwp <- gaussian_weight(ws$position,c,sig)*ws[,parameter]*(ws[,"nk"] - 1)
        gws <- gaussian_weight(ws$position,c,sig)*(ws[,"nk"] - 1)
      }
      else{
        gwp <- gaussian_weight(ws$position,c,sig)*ws[,parameter]
        gws <- gaussian_weight(ws$position,c,sig)
      }

      #print(sum(gws))
      out[i,new_col] <- (sum(gwp)/sum(gws))
      out[i,"n_snps"] <- nrow(ws)
    }
    return(out)
  }
  else{
    #cat("Using fixed window, slide =", fixed_window, "kb.\n")
    out_const <- data.frame()
    num_snps <- nrow(x)
    these_snps <- sort(x$position)
    last_snp_pos <- these_snps[length(these_snps)]
    first_snp_pos <- these_snps[1]
    size <- last_snp_pos - first_snp_pos
    c <- 0
    i <- 1
    u <- 1
    while (c <= size){
      if(u %% 50 == 0){cat("Window Postion:", c, "out of", size, "\n")}
      #cat(c, start, end, "\n")
      start <- c - 3*sig
      end <- c + 3*sig
      ws <- x[x$position <= end & x$position >= start,]
      ws <- ws[!is.na(ws[,parameter]),] #remove na values

      if(nrow(ws) != 0){
        if(!nk_weight){
          gwp <- gaussian_weight(ws$position,c,sig)*ws[,parameter]*(ws[,"nk"] - 1)
          gws <- gaussian_weight(ws$position,c,sig)*(ws[,"nk"] - 1)
        }
        else{
          gwp <-  gaussian_weight(ws$position,c,sig)*ws[,parameter]
          gws <-  gaussian_weight(ws$position,c,sig)
        }
        #print(gwp)
        #print(gws)
        out_const[i,"position"] <- c
        out_const[i,new_col] <- (sum(gwp)/sum(gws))
        out_const[i,"n_snps"] <- nrow(ws)
      }
      c <- c + fixed_window*1000
      i <- i + 1
      u <- u + 1
    }
    return(out_const)
  }
}

