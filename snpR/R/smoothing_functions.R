#function for calculating gaussian weight of a point
gaussian_weight <- function(p, c, s) {
  exp(-(p-c)^2/(2*s^2))
}

#function to generate sliding averages based on gaussian weights of points around each point,
#cut off at 3sigma where weight becomes very low. Sigma should be given in kbp, although input should be
# in bp, in the second column of each row.
smoothed_ave <- function(x, parameter, sigma, nk_weight = FALSE, fixed_window = NULL) {
  sig <- 1000*sigma
  new_col <- paste0("smoothed", "_", parameter)
  if(!nk_weight){
    x[,"nk"] <- as.numeric(x[,"nk"])
  }
  x$position <- as.numeric(x$position)
  if (is.null(fixed_window)){
    out <- data.frame()
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
        out_const[i,"group"] <- x[1,"group"]
        out_const[i,"position"] <- c
        out_const[i,new_col] <- (sum(gwp)/sum(gws))
        out_const[i,"n_snps"] <- nrow(ws)
      }
      c <- c + fixed_window*1000
      i <- i + 1
      u <- u + 1
    }
    out_const <- out_const[!is.na(out_const$group),]
    return(out_const)
  }
}

