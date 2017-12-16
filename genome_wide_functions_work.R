library(plyr)
library(reshape2)

#function to calc pi from data. should be in format with num alleles in col 6, total count of alleles
#in col 5, and subsequent alleles counts in rows 6+ (allele count/bayescan format, as given by format_snps option one with pops). 
#If any allele
calc_pi <- function(data){
  data[,5] <- as.numeric(data[,5])
  data[,6] <- as.numeric(data[,6])
  data[,7] <- as.numeric(data[,7])
  data[,8] <- as.numeric(data[,8])
  num <- data[,6]
  if(!all(num != 2)){
    warning("Some loci do not have two alleles in some populations.\n")
  }
  
  n1 <- data[,7]
  n2 <- data[,8]
  nt <- data [,5]
  binom_n1 <- choose(n1,2)
  binom_n2 <- choose(n2,2)
  binom_nt <- choose(nt,2)
  #print(n1)
  #print(n2)
  p <- 1 - (binom_n1 + binom_n2)/binom_nt
  #print(pi)
  
  
  return(p)
}

#calculates a pop weighted stat. Needs columns named "pi", "n_total", and one with a name matching
#the "stat" argument, which should be a character string. "boots" argument is the number of bootstraps
#to do do get standard errors
calc_weighted_stat <- function(x, stat, boots = 30){
  pops <- unique(x$pop)
  out <- matrix(0, length(pops), 3)
  out[,1] <- pops
  get_wm <- function(y){
    w <- y$n_total
    return(sum(w*y[,stat])/sum(w))
  }
  cat("working on pop:\n")
  for(i in 1:length(pops)){
    #mean
    cat("\n\t", pops[i],"\n")
    y <- x[x$pop == pops[i],]
    y <- y[!is.na(y[,stat]),]
    out[i,2] <- get_wm(y)
    
    
    #se, via bootstrap.
    n <- nrow(y)
    bvals <- numeric(boots)
    for(j in 1:boots){
      if(j %% 10 == 0){cat("\t\tboot:", j, "\n")}
      b <- y[sample(n, n, T),]
      bvals[j] <- get_wm(b)
    }
    out[i,3] <- sd(bvals)
  }
  out <- as.data.frame(out)
  colnames(out) <- c("pop", paste0(stat, "_wm"), paste0(stat,"_se"))
  out[,2] <- as.numeric(as.character(out[,2]))
  out[,3] <- as.numeric(as.character(out[,3]))
  return(as.data.frame(out))
}

#function for calculating gaussian weight of a point
gaussian_weight <- function(p, c, s) {
  exp(-(p-c)^2/(2*s^2))
}

#function to generate sliding averages based on gaussian weights of points around each point,
#cut off at 3sigma where weight becomes very low. Sigma should be given in kbp, although input should be
# in bp, in the second column of each row. 
smoothed_ave <- function(data, parameter, sigma, nk_weight = NULL, fixed_window = NULL) {
  sig <- 1000*sigma
  new_col <- paste0("smoothed", "_", parameter)
  if(!is.null(nk_weight)){
    data[,nk_weight] <- as.numeric(nk_weight)
  }
  data$position <- as.numeric(data$position)
  if (is.null(fixed_window)){
    out <- data.frame()
    for (i in 1:nrow(data)) {
      if(i %% 1000 == 0){cat("snp ", i, "\n")}
      c <- data[i,"position"]
      start <- c - 3*sig
      end <- c + 3*sig
      #cat("c:", c, "sig:", sig, "\n")
      ws <- data[data$position <= end & data$position >= start,]
      ws <- ws[!is.na(ws[,parameter]),] #remove na values

      if(!is.null(nk_weight)){
        gwp <- gaussian_weight(ws$position,c,sig)*ws[,parameter]*(ws[,nk_weight] - 1)
        gws <- gaussian_weight(ws$position,c,sig)*(ws[,nk_weight] - 1)
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
    num_snps <- nrow(data)
    these_snps <- sort(data$position)
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
      ws <- data[data$position <= end & data$position >= start,]
      ws <- ws[!is.na(ws[,parameter]),] #remove na values
      
      if(nrow(ws) != 0){
        if(!is.null(nk_weight)){
          gwp <- gaussian_weight(ws$position,c,sig)*ws[,parameter]*(ws[,nk_weight] - 1)
          gws <- gaussian_weight(ws$position,c,sig)*(ws[,nk_weight] - 1)
        }
        else{  
          gwp <-  gaussian_weight(ws$position,c,sig)*ws[,parameter]
          gws <-  gaussian_weight(ws$position,c,sig)
        }
        #print(gwp)
        #print(gws)
        out_const[i,"group"] <- data[1,"group"]
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

#generate smoothed data for each linkage group.
smooth_w_groups<- function(data, parameter, sigma, nk_weight = FALSE, fixed_window = FALSE, slide = "NA"){
  par <- as.character(parameter)
  w_df<- data.frame()
  groups <- unique(data$group)
  if(fixed_window == FALSE){
    for (i in 1:length(groups)){
      w_data <- subset(data, group == groups[i])
      print(groups[i])
      w_data <- smoothed_ave(w_data, par, sigma, nk_weight)
      w_df <- rbind(w_df, w_data)
    }
  }
  else{
    for (i in 1:length(groups)){
      w_data <- subset(data, group == groups[i])
      print(groups[i])
      w_out <- smoothed_ave(data = w_data, parameter = par, sigma = sigma, nk_weight = nk_weight, fixed_window = fixed_window, slide = slide)
      w_df <- rbind(w_df, w_out)
    }
  }
 return(w_df) 
}

#function to use smoth_ave on a file with multiple groups and populations that all need to be seperately
#smoothed.
smooth_w_groups_and_pops <- function(data, parameter, sigma, nk_weight = FALSE, fixed_window = FALSE, slide = "NA"){
  par <- as.character(parameter)
  w_df<- data.frame()
  pops <- unique(data$pop)
  groups <- unique(data$group) 
  for (i in 1:length(groups)){
    w_data <- subset(data, group == groups[i])
    print(groups[i])
    for(j in 1:length(pops)){
      print(pops[j])
      x_data <- smoothed_ave(subset(w_data, pop == pops[j]), par, sigma, nk_weight, fixed_window, slide)
      w_df <- rbind(w_df, x_data)
    }
  }
  return(w_df)
}


#Returns the average number of snps per window of size 6 sigma (in kb).
window_average_snps <- function(data, sigma) {
  sig <- 1000*sigma
  for (i in 1:nrow(data)) {
    c <- data[i,"position"]
    start <- c - (3*sig)
    end <- c + (3*sig)
    count <- 0
    for (j in 1:nrow(data)) {
      p <-  data[j,"position"]
      if(start < p & p < end) {
          count = count + 1
      }
    }
    #print(sum(gws))
    data[i,"count"] <- (count)
    #print(sum(gwp)/sum(gws))
    #print(sum(gwp)/length(gwp))
    #data[i,new_col] <- mean(gwp)
  }
  remove(i)
  return(data)
} 

#Runs the window_average_snps with multiple groups of samples, denoted in the column "group"
window_w_groups<- function(data, sigma){
  w_df<- data.frame()
  groups <- unique(data$group) 
  for (i in 1:length(groups)){
    w_data <- subset(data, group == groups[i])
    print(groups[i])
    w_data <- window_average_snps(w_data, sigma)
    w_df <- rbind(w_df, w_data)
  }
  return(w_df)
}

#Runs the window_average_snps with multiple groups and pop, denoted in the columns "group" and "pop".
window_w_groups_and_pops <- function(data, sigma){
  w_df<- data.frame()
  pops <- unique(data$pop)
  groups <- unique(data$group) 
  for (i in 1:length(groups)){
    w_data <- subset(data, group == groups[i])
    print(groups[i])
    for(j in 1:length(pops)){
      print(pops[j])
      x_data <- window_average_snps(subset(w_data, pop == pops[j]), sigma = sigma)
      w_df <- rbind(w_df, x_data)
    }
  }
  return(w_df)
}

#Calculates the difference between stat values  (pi?) in identical snps between two different populations.
#Identical snps should be sequential in input:
#row one: SNP(a)position \t POP(1) \t STAT \t linkage_group
#row two: SNP(a)position \t POP(2) \t STAT \t linkage_group
#Positive values indicate that stat is greater in pop 1, negative values indicate that stat is greater in pop 2.
#Expects row names group (linkage group), pop (population), stat (statistic, supplied as argument).
#position(position of snp), and snps (name of snp)
Calc_stat_dif<- function(data, stat){
  i <- 1
  j <- 1
  pops <- unique(data$pop)
  df <- data.frame()
  new_col1 <- paste0(stat, "_", pops[1])
  new_col2 <- paste0(stat, "_", pops[2])
  while (i < nrow(data)){
    stat1 <- data[i, stat]
    stat2 <- data[i+1, stat]
    diff <- data[i, stat] - data[i+1, stat]
    df[j,"position"] <- data[i, "position"]
    df[j,"snp"] <- data[i, "snps"]
    df[j,"group"] <- data[i,"group"]
    df[j,new_col1] <- stat1
    df[j,new_col2] <- stat2
    df[j,paste0(stat, "_", "diff")] <- diff
    j <- j + 1
    i <- i + 2
  }
  return(df)
}


#Calculates average LD (both D'and r^2) between a snp and all snps in a 6*sigma window. 
#num_extra is number of extra data containing columns at the end. DO NOT COUNT ave_rsq and
#ave_Dprime on a re-run! Note, will give an error about taking the mean of an non-numeric or
#logical vector if a snp has no other snps within its 6sigma window. Often happens at the ends
#of contigs. 
LD_ave <- function(data, sigma, num_start, num_extra){
  sig <- 1000*sigma
  i <- 1
  data[,"ave_Dprime"] <- NA
  data[,"ave_rsq"] <- NA
  for(i in i:nrow(data)){
    if(i %% 100 == 0){cat("snp number: ", i, "\n")}
    g_list <- sort(unique(unlist(data[i,num_start:(ncol(data) - num_extra)])))
    if(g_list[1] == "0000"){g_list <- g_list[2:length(g_list)]}
    num <- length(g_list)
    if(num == 1 | num > 3){
      #cat(num, "genotypes detected at snp:", data[i,1], ", Position:", data[i, "position"], ". This loci will be skipped.\n")
      next
    }
    else if (num == 3 & substr(g_list[1],1,2) == substr(g_list[1],3,4) & substr(g_list[2],1,2) == substr(g_list[2],3,4) & substr(g_list[3],1,2) == substr(g_list[3],3,4)){
      #cat("Three homozygous genotypes detected at snp:", data[i,1], ", Position:", data[i, "position"], ". This loci will be skipped, and NA listed for D' and r^2.\n")
      data[i,"ave_Dprime"]<- NA
      data[i,"ave_rsq"]<- NA
      next
    }
    key1 <- as.character(data[i, "position"])
    if(matrix_out == TRUE){
      out_Dprime[i,1] <- key1
      out_rsq[i,1] <- key1
    }
    c <- data[i,"position"]
    
    As <- character(2*length(g_list))
    h <- 1
    ch <- 1
    while(h <= length(g_list)){ #get all possible allele ones and allele twos
      As[ch] <- substr(g_list[h],1,2)
      As[ch+1] <- substr(g_list[h],3,4)
      ch <- ch + 2
      h <- h + 1
    }
    As <- sort(unique(As)) #get As.
    A1 <- As[1]
    A2 <- As[2]
    j <- 1
    #print(g_list)
    for(j in j:nrow(data)){
      #cat("\nj = ", j, "\n")
      p <- data[j,"position"]
      #cat(data[i,1], data[j,1], "\n")
      #print(out_Dprime)
      #print(out_rsq)
      key2 <- as.character(data[j, "position"])
      #cat("\np", p, "c", c, abs(p-c), "\n")
      if (p<c){
        next
      }
      else if (p == c){
        Dprime <- NA
        rsq <- NA
      }
      else{
        #cat("\nPosition:", c, "vs. Position", p, "\n")
        g_list2 <- sort(unique(unlist(data[j,num_start:(ncol(data) - num_extra)])))
        #print(g_list2)
        if(g_list2[1] == "0000"){g_list2 <- g_list2[2:length(g_list2)]}
        num <- length(g_list2)
        if(num == 1 | num > 3){
          #cat("\n Too few or too many allels!\n")
          next
        }
        else if (num == 3 & substr(g_list2[1],1,2) == substr(g_list2[1],3,4) & substr(g_list2[2],1,2) == substr(g_list2[2],3,4) & substr(g_list2[3],1,2) == substr(g_list2[3],3,4)){
          stop("Three homozygous genotypes detected at snp.\n")
        }
        k <- num_start
        Bs <- character(2*length(g_list))
        h <- 1
        ch <- 1
        while(h <= length(g_list)){ #get all possible allele ones and allele twos
          Bs[ch] <- substr(g_list2[h],1,2)
          Bs[ch+1] <- substr(g_list2[h],3,4)
          ch <- ch + 2
          h <- h + 1
        }
        Bs <- sort(unique(Bs)) #get Bs.
        B1 <- Bs[1]
        B2 <- Bs[2]
        #cat("\nSnp", data[i,1], "alleles: ", A1, A2, "Snp", data[j,1], "alleles: ", B1, B2, "\n")
        A1B1 <- paste0(A1,B1)
        A2B1 <- paste0(A2,B1)
        A1B2 <- paste0(A1,B2)
        A2B2 <- paste0(A2,B2)
        w_df <- c()
        w_df[A1B1] <- 0
        w_df[A2B1] <- 0
        w_df[A1B2] <- 0
        w_df[A2B2] <- 0
        #print(g_list2)
        for (k in k:(ncol(data) - num_extra)){
          #cat(data[i,k], data[j,k], "\n")
          if(data[i,k] == "0000" | data[j,k] == "0000"){
            next()
          }
          i1 <- substr(data[i,k],1,2)
          i2 <- substr(data[i,k],3,4)
          j1 <- substr(data[j,k],1,2)
          j2 <- substr(data[j,k],3,4)
          if(i1 != i2 & j1 != j2){ #if double hets...
            next() #skip it, can't determine haplotypes
          }
          if(i1 == i2 & j1 == j2){#if double homo
            i1j1 <- paste0(i1,j1) #get haplotype
            w_df[i1j1] <- w_df[i1j1] + 2 #add two to haplotype count
            next() #move to next individual
          }
          #otherwise, have one het and one homo, add one count to each haplotype
          #list possible haplotypes, should be two unique entries
          i1j1 <- paste0(i1,j1)
          i2j1 <- paste0(i2,j1)
          i1j2 <- paste0(i1,j2)
          i2j2 <- paste0(i2,j2)
          #get unique haplotypes
          haps <- unique(c(i1j1, i2j1, i1j2, i2j2))
          if (length(haps) != 2){
            cat("Error, more than two haploypes for hom/het at loci pair:", data[i,1], data[j,1], "ind:",
                colnames(data)[k], "\n")
            stop()
          }
          #add to appropriate counts, one to heach
          w_df[haps[1]] <- w_df[haps[1]] + 1
          w_df[haps[2]] <- w_df[haps[2]] + 1
        }
        #cat(i, j, k, "\n")
        if(is.na(sum(w_df) == TRUE)){browser()}
        if(sum(w_df > 0) == 1 | sum(w_df > 0) == 0){
          next() #skip the comparison if only one observed haplotype.
        }
        A1B1f <- w_df[A1B1]/sum(w_df)
        A2B1f <- w_df[A2B1]/sum(w_df)
        A1B2f <- w_df[A1B2]/sum(w_df)
        A2B2f <- w_df[A2B2]/sum(w_df)
        A1f <- A1B1f + A1B2f
        A2f <- 1 - A1f
        B1f <- A1B1f + A2B1f
        B2f <- 1 - B1f
        D <- A1B1f - A1f*B1f
        D2 <- A2B1f - A2f*B1f
        #print(w_df)
        #cat("\n")
        #print(w_df)
        #cat("\n")
        #cat("D:", D, "\n")
        #cat("temp p1q2 and p2q1", temp, "\n")
        if(D > 0){
          #cat("D is positive.\n")
          Dmax <- min(A1f*B2f, A2f*B1f)
          Dprime <- D/Dmax
        }
        else if (D < 0 ){
          #cat("D is negative.\n")
          Dmax <- max(-A1f*B1f, -A2f*B2f)
          Dprime <- D/Dmax
        }
        else if (D == 0){
          #cat("D is zero.\n")
          Dprime <- 0
        }
        rsq <- (D^2)/(A1f*A2f*B1f*B2f)
        if(is.na(rsq) == TRUE){
          next()
        }
        chisqu <- rsq*(4)
        pval <- 1 - pchisq(chisqu, 1)
        #cat("\nDprime = ", Dprime, "rsq = ", rsq, "\n")
        Dprimes <- c(Dprimes, Dprime)
        rsqs <- c(rsqs, rsq)
        pvals <- c(pvals, pval)
      }
    }
    #if(is.na(mean(rsqs) == TRUE)){cat("\nSnp:", data[i,1], "rsqs:", rsqs, "Mean:", mean(rsqs), "\n")}
    data[i,"ave_Dprime"] <- mean(Dprimes)
    data[i,"ave_rsq"] <- mean(rsqs)
    data[i,"ave_pval"] <- mean(pval)
  }
  return(data)
}

#Calls the LD_ave function for a data file with multiple linkage groups, noted in the "group" column.
LD_w_groups<- function(data, sigma, num_extra){
  w_df<- data.frame()
  groups <- unique(data$group) 
  for (i in 1:length(groups)){
    w_data <- subset(data, group == groups[i])
    #print(ncol(w_data))
    print(groups[i])
    w_data <- LD_ave(data = w_data, sigma = sigma, num_extra = num_extra)
    w_df <- rbind(w_df, w_data)
  }
  return(w_df)
}



#Takes an input file of genes with start/end position noted and maps back the average
#statitstic of choice for that gene from another provided data set with positions. Requires the "zoo"
#package
gene_ave_stat <- function(gene_data, stat_data, stat){
  require(zoo)
  out_df <- data.frame()
  w_df <- cbind.data.frame(stat_data[,"position"], stat_data[,stat])
  colnames(w_df) <- c("position", stat)
  #print(w_df)
  mids <- c()
  count <- 1
  i <- 1
  count2 <- 1
  store_df <- data.frame()
  cat("Reading and storing data from input files...\n")
  for (i in i:nrow(gene_data)){
    if(i %% 100 == 0){cat("Reading gene number: ", i, "\n")}
    m <- mean(gene_data[i, "start"], gene_data[i,"end"])
    if(m %in% w_df$position == TRUE){
      out_df[count,"probeID"] <- gene_data[i,"probeID"]
      out_df[count,stat] <- w_df[match(m, w_df[,"position"]),stat]
      out_df[count,"position"] <- m
      out_df[count,"group"] <- gene_data[i,"group"]
      count <- count + 1
    }
    else{
      mids <- c(mids, m)
      store_df[count2,"probeID"] <- gene_data[i,"probeID"]
      store_df[count2,"mid"] <- m
      store_df[count2,"group"] <- gene_data[i,"group"]
      count2 <- count2 + 1
    }
  }
  #print(store_df)
  cat("Done.\n")
  new_stats <- cbind(position = mids, stat = "NA")
  colnames(new_stats) <- c("position", stat)
  w_df <- rbind(w_df, new_stats)
  #print(w_df)
  #print(str(w_df))
  w_df[,stat] <- as.numeric(w_df[,stat])
  w_df[,"position"] <- as.numeric(w_df[,"position"])
  w_df <- w_df[order(w_df$"position"),]
  #print(w_df)
  #print(str(w_df))
  cat("Interpolating", stat, "for gene midpoints...\n")
  zw_df <- zoo(w_df)
  index(zw_df) <- zw_df[,"position"]
  zw_df <- na.approx(zw_df)
  index(zw_df) <- 1:nrow(w_df)
  #print(zw_df)
  w_df <- as.data.frame(zw_df)
  colnames(w_df) <- c("position", stat)
  cat("Done.\nPrinting to output...\n")
  i <- 1
  #print(head(w_df))
  for(i in i:length(mids)){
    if(i %% 100 == 0){cat("Printing gene number: ", i, "\n")}
    #print(w_df[match(mids[i], w_df[,"position"]),stat])
    out_df[count,"probeID"] <- store_df[match(mids[i], store_df[,"mid"]),"probeID"]
    out_df[count, stat] <- w_df[match(mids[i], w_df[,"position"]),stat]
    out_df[count, "position"] <- mids[i]
    out_df[count, "group"] <- store_df[match(mids[i], store_df[,"mid"]),"group"]
    count <- count + 1
  }
  return(out_df)
}

#calls gene_ave_stat for data with multiple groups
gene_ave_stat_w_groups <- function(gene_data, stat_data, stat){
  groups <- unique(as.character(stat_data$group))
  i <- 1
  w_df <- data.frame()
  for(i in i:length(groups)){
    cat("\n----------------------\nGroup:", groups[i], "\n----------------------\n")
    w_data <- gene_ave_stat(gene_data = subset(gene_data, group == groups[i]), stat_data = subset(stat_data, group == groups[i]), stat = stat)
    w_df <- rbind(w_df, w_data)
  }
  return(w_df)
}


#calls smooth_w_groups for multiple FST comparisons, stored in sequential columns. For nk weighting,
#ensure that the nk_weight data set is sorted the same as the input data snp by snp and has comp pops
#sorted alphabetically, as given by pairwise_nk. Output only keeps data relevent to comparison.
FST_smooth_multicomp <- function(data, sigma, start_col, nk_weight = FALSE, nk_df = NA){
  list <- colnames(data)
  #cat(i, "next data\n")
  outlist <- list()
  i <- start_col
  data$position <- as.numeric(data$position)
  if(nk_weight == FALSE){
    for(i in i:length(list)){
      #cat(i, "next data\n")
      cat("\n\nWorking on:", list[i], "\n\n")
      name <- paste0(list[i], "_fst")
      w_df <- subset(data, data[,i] != "-")
      #cat ("before\n")
      #print(test)
      w_df[,list[i]] <- as.numeric(as.character(w_df[,list[i]]))
      x_df <- w_df[,1:3]
      #print(x_df)
      x_df[,list[i]] <- w_df[,list[i]]
      #cat("after\n")
      #print(test)
      x_df <- smooth_w_groups(data = x_df, parameter = list[i], sigma = sigma)
      #print(name)
      outlist[[name]] <- x_df
      i <- i + 1
    }
  }
  else{
    for(i in i:length(list)){
      #cat(i, "next data\n")
      cat("\n\nWorking on:", list[i], "\n\n")
      name <- paste0(list[i], "_fst")
      comp <- sort(unlist(strsplit(list[i], split = "_")))
      comp <- paste0(comp[1], "_", comp[2], "_nk")
      data$nk <- as.numeric(nk_df[,comp])
      w_df <- subset(data, data[,i] != "-")
      #cat ("before\n")
      #print(test)
      w_df[,list[i]] <- as.numeric(as.character(w_df[,list[i]]))
      x_df <- w_df[,1:3]
      #print(x_df)
      x_df[,list[i]] <- w_df[,list[i]]
      x_df[,"nk"] <- w_df[,"nk"]
      #cat("after\n")
      #print(head(x_df))
      x_df <- smooth_w_groups(data = x_df, parameter = list[i], sigma = sigma, nk_weight = TRUE)
      #print(name)
      outlist[[name]] <- x_df
      i <- i + 1
    }
  }
  return(outlist)
}



#gets the position and group of a set of data where the first column is: group_start_end and the second
#is position in tag. Accepts one or two substring group names (groupXXX_start_end_F/R or group_XXX_start_end_F/R).
get_position_group <- function(data){
  data[,1] <- as.character(data[,1])
  cat("Subseting out identifier rows...")
  w_df <- data[,1:2]
  cat("done.\n")
  #print(w_df)
  i <- 1
  for(i in i:nrow(w_df)){
    if(i %% 10000 == 0){cat("locus number: ", i, "\n")}
    w_v <- unlist(strsplit(w_df[i,1], "_"))
    if(length(w_v) == 4){
      w_df[i,"group"] <- w_v[1]
      start <- w_v[2]
    }
    else if(length(w_v) == 5){
      w_df[i,"group"] <- paste0(w_v[1], "_", w_v[2])
      start <- w_v[3]
    }
    else{
      warning("Incorrect location identifier format.")
      stop
    }
    w_df[i,"position"] <- as.numeric(start) + (w_df[i,2] - 1)
  }
  out_df <- w_df[,3:4]
  cat("done.\n")
  return(out_df)
}

#runs get_position_group over multiple processors. For large datasets, this speeds up things considerably.
get_position_group_par <- function(data, num_cores){
  data[,1] <- as.character(data[,1])
  cat("Subseting out identifier rows...")
  w_df <- data[,1:2]
  cat("done.\nRunning on ", num_cores, "cores.")
  #print(w_df)
  require(doParallel)
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  output <- foreach(i=1:nrow(w_df), .inorder = TRUE, .combine = rbind) %dopar%{
    #if(i %% 10 == 0){cat("locus number: ", i, "\n")}
    w_v <- unlist(strsplit(w_df[i,1], "_"))
    if(length(w_v) == 4){
      group <- w_v[1]
      start <- w_v[2]
    }
    else if(length(w_v) == 5){
      group <- paste0(w_v[1], "_", w_v[2])
      start <- w_v[3]
    }
    else{
      warning("Incorrect location identifier format.")
      stop
    }
    c(group, as.numeric(start) + (w_df[i,2] - 1))
  }
  #print(output)
  cat("Pasting together output...")
  name_df <- as.data.frame(output)
  colnames(name_df) <- colnames(w_df)
  out_df <- cbind(name_df,data[,3:ncol(data)])
  rownames(out_df) <- c(1:nrow(out_df))
  cat("done.\n")
  stopCluster(cl)
  registerDoSEQ()
  return(out_df)
}

#estimate S in one generation
estimate_selection <- function(p_before, p_after){
  top <- 1-(p_before/p_after)
  q <- 1 - p_before
  bottom <- q^2
  s <- top/bottom
  return(s)
}



#Genotype individuals for genomic inversion spanning multiple SNPs as single psuedo SNP. Requires a list of
#SNPs on the inversion, with a column named "snp" that identifies snps, as well as reference data
#from which the major and minor alleles (fixed differences between the inversion and the 
#ancestral genome) will be called. These can be the same file, and can also be the same as the data to be
#genotyped, but at least the minor reference must contain at least one instance of an individual with
#rarer genotype (either a copy of the inversion or ancestral genome). Out = 1 returns the inversion (or rare
#ancestral) genotype for each individual as a psuedo SNP, out = 2 calls each snp genotype for each individual as
#either Maj or Minor (anc or inv), out = 3 returns a list of counts of each MinMaj genotype across all SNPs for each individual.
#Header and Footer col counts refer to the number of extra columns at either the start or end of the data,
#including the required SNP names. The calling cuttoff is the percent (in decimals) of SNPs which must be 
#in agreement for the inversion to be genotyped. Min nonzeroes denotes the number of ungenotyped snps allowed
#for the genotype of the inversion to be called (note that if the SNPs provided in the SNP list are truely
#fixed and there is no genotyping error, one SNP is all that is neccisary. Multiple are likely more accurate,
#however.).
genotype_inversion <- function(data, SNP_list, Min_ref, Maj_ref, header_col_count, calling_cuttoff, min_non_zeros, footer_col_count, out = 1){
  i <- 1
  data <- data[,1:ncol(data)-footer_col_count]
  Min_ref <- Min_ref[,1:ncol(Min_ref)-footer_col_count]
  Maj_ref <- Maj_ref[,1:ncol(Maj_ref)-footer_col_count]
  require(plyr)
  div_snps <- subset(data, data$snp %in% SNP_list)
  div_snps <- arrange(div_snps, snp)
  Min_ref <- subset(Min_ref, Min_ref$snp %in% SNP_list)
  Maj_ref <- subset(Maj_ref, Maj_ref$snp %in% SNP_list)
  Min_ref <- arrange(Min_ref, snp)
  Maj_ref <- arrange(Maj_ref, snp)
  #print(head(Min_ref))
  #print(Maj_ref)
  div_snps$group <- as.character(div_snps$group)
  ind_state <- data.frame(matrix(ncol = ncol(div_snps), nrow = nrow(div_snps)))
  #print(ind_state)
  colnames(ind_state) <- colnames(div_snps)
  #print(names(div_snps))
  
  #Create list of major and minor alleles from the list of fixed snps and two (or one) sets of reference samples.
  FixInvMin <- apply(Min_ref[,(header_col_count + 1):length(Min_ref)],1,function(y)
    names(which.min(table(subset(c(substr(y, 1,2), substr(y,3,4)), c(substr(y, 1,2), substr(y,3,4))[] != "00")))))
  #print(FixInvMin)
  FixInvMaj <- apply(Maj_ref[,(header_col_count + 1):length(Maj_ref)],1,function(x)
    names(which.max(table(subset(c(substr(x, 1,2), substr(x,3,4)), c(substr(x, 1,2), substr(x,3,4))[] != "00")))))
  #print(FixInvMaj)
  InvMinMaj <- data.frame(Maj = FixInvMaj, Min = FixInvMin)
  remove(FixInvMin)
  remove(FixInvMaj)
  #print(InvMinMaj)
  
  #Call each genotype according to Major/Minor allele identities (01 and 02, respectively), write to ind_state.
  for(i in i:nrow(div_snps)){
    #print(div_snps[i,1:header_col_count])
    j <- header_col_count + 1
    ind_state[i,1:header_col_count] <- div_snps[i,1:header_col_count]
    for(j in j:ncol(div_snps)){
      a1 <- substr(div_snps[i,j],1,2)
      a2 <- substr(div_snps[i,j],3,4)
      if(a1 == a2){
        if(a1 == InvMinMaj[i,"Maj"]){
          state <- "0101"
        }
        else if (a1 == InvMinMaj[i,"Min"]){
          state <- "0202"
        }
        else if (a1 == "00"){
          state <- "0000"
        }
        else{
          cat("Undocumented allele at snp:", div_snps[i,1], ",individual:", j - header_col_count, "!\n")
          state <- "other"
        }
      }
      else{
        if ((a1 == InvMinMaj[i,"Maj"] | a1 == InvMinMaj[i,"Min"]) & (a2 == InvMinMaj[i,"Maj"] | a2 == InvMinMaj[i, "Min"])){
          state <- "0102"
        }
        else{
          cat("Undocumented allele at snp:", div_snps[i,1], ",individual:", j - header_col_count, "!\n")
          state <- "other"
        }
      }
      ind_state[i,j] <- state
    } 
  }
  if(out == 2){
    return(ind_state)
    stop
  }
  
  #Call individual genotype for whole inversion
  if (out == 1){
    inv_state <- data.frame(individual = numeric(0), genotype = numeric(0), percent_snp_cons = numeric(0))
  }
  i <- 1
  j <- header_col_count + 1
  if(out == 3){
    c_df <- data.frame(individual = numeric(0), g0000 = numeric(0), g0101 = numeric(0), g0102 = numeric(0), g0202 = numeric(0), other = numeric(0))
  }
  for (j in j:ncol(ind_state)){
    i <- 1
    if(out == 1){
      inv_state[j - header_col_count,"individual"] <- colnames(ind_state)[j]
    }
    else if (out == 3){
      c_df[j - header_col_count,"individual"] <- colnames(ind_state)[j]
    }
    cons0000 <- 0
    cons0101 <- 0
    cons0102 <- 0
    cons0202 <- 0
    consother <- 0
    for(i in i:nrow(ind_state)){
      if(ind_state[i,j] == "0000"){
        cons0000<- cons0000 + 1
      }
      else if(ind_state[i,j] == "0101"){
        cons0101 <- cons0101 + 1
      }
      else if(ind_state[i,j] == "0102"){
        cons0102 <- cons0102 + 1
      }
      else if(ind_state[i,j] == "0202"){
        cons0202 <- cons0202 + 1
      }
      else if(ind_state[i,j] == "other"){
        consother <- consother + 1
      }
      else{
        cat("Error in inv_state calling, ind_state col/row:", i, "/", j, ".")
        stop
      }
    }
    tot <- cons0202 + cons0101 + cons0102 + consother
    #cat("Ind: ", colnames(ind_state)[j], "\t0000:", cons0000, "\t0101: ", cons0101, "\t0102: ", cons0102, "\t0202: ", cons0202, "\tother:", consother, "\n\n")
    if(out == 3){
      c_df[j - header_col_count, "g0000"] <- cons0000
      c_df[j - header_col_count, "g0101"] <- cons0101
      c_df[j - header_col_count, "g0102"] <- cons0102
      c_df[j - header_col_count, "g0202"] <- cons0202
      c_df[j - header_col_count, "other"] <- consother
    }
    else if (out == 1){
      if (tot < min_non_zeros){
        inv_state[j - header_col_count, "genotype"] <- "0000"
        inv_state[j - header_col_count, "percent_snp_cons"] <- cons0000/(tot + cons0000)
      }
      else if (cons0101/tot >= calling_cuttoff){
        inv_state[j - header_col_count, "genotype"] <- "0101"
        inv_state[j - header_col_count, "percent_snp_cons"] <- cons0101/tot
      }
      else if (cons0102/tot >= calling_cuttoff){
        inv_state[j - header_col_count, "genotype"] <- "0102"
        inv_state[j - header_col_count, "percent_snp_cons"] <- cons0102/tot
      }
      else if (cons0202/tot >= calling_cuttoff){
        inv_state[j - header_col_count, "genotype"] <- "0202"
        inv_state[j - header_col_count, "percent_snp_cons"] <- cons0202/tot
      }
      else if (consother/tot >= calling_cuttoff){
        inv_state[j - header_col_count, "genotype"] <- "other"
        inv_state[j - header_col_count, "percent_snp_cons"] <- consother/tot
      }
      else {
        inv_state[j - header_col_count, "genotype"] <- "0000"
        inv_state[j - header_col_count, "percent_snp_cons"] <- max(c(consother, cons0202, cons0102, cons0102))/tot
      }
    }
  }
  if(out == 1){
    return(inv_state)
  }
  else if (out == 3){
    return(c_df)
  }
}



#Calculates Tajima's theta/pi, Waterson's theta, and Tajima's D over a sliding window. Takes the bayescan/allele count format. Pi calculated as in hohenlohe 2010.
#inputs: x: input data, first three columns are: snp ID, position, snp scaffold/lg/chr. Named required columns are "n_total", "ni1", and "ni2". Named column "pop" is optional, and will be returned in output if exists.
#        ws: window size in kbp
#        step: step size in kbp
#        report: give progress every "report" windows.
#output: A data frame where column one is group, two is position, three is Tajima's theta,
#        four is Waterson's theta, and five is D. Note that both thetas are reported as a
#        frequency (so the number per base). Full number per window is ws*theta.
Tajimas_D <- function(x, ws, step, report = 20){
  
  #windowing
  out <- data.frame() #initialize output
  tps <- sort(x[,2]) #get the snp positions, sort
  lsp <- tps[length(tps)] #get the position of the last snp to use as endpoint
  c <- 0 #set starting position
  i <- 1 #set starting iteration for writing output
  ws <- 1000*ws
  while (c <= lsp){
    start <- c - ws #change window start
    end <- c + ws #change window end
    if(i %% report == 0){cat("Window Postion:", c, "out of", lsp, "\n")}
    
    #take all the snps in the window, calculate T's theta, W's theta, and T's D
    wsnps <- x[x[,2] <= end & x[,2] >= start,] #get only the snps in the window
    wsnps <- wsnps[wsnps[,"ni1"] > 0 & wsnps[,"ni2"] > 0,] #remove any snps that are fixed in this pop/group/whatever
    if(nrow(wsnps) == 0){ #if no snps in window
      c <- c + step*1000 #step along
      next #go to the next window
    } 
    b1s <- choose(wsnps[,"ni1"],2) #binomial for first allele
    b2s <- choose(wsnps[,"ni2"],2) #binomial for second allele
    bts <- choose(wsnps[,"n_total"],2) #binomial for all alleles
    ts.theta <- sum(1-((b1s+b2s)/bts)) #average number of pairwise differences (pi) per snp. Equivalent to sum((ndif/ncomp)) for all snps
    ts.thetaf <- ts.theta/ws #pi for the whole region, which includes all of the non-polymorphic sites. Reason why this shouldn't be run with a maf filter, probably.
    n_seg <- nrow(wsnps) #number of segregating sites
    K <- round(mean(wsnps[,"n_total"])) #average sample size for ws.theta, as in hohenlohe 2010. Alternative would make this into a function, then use lapply on the vector of K values
    #if(is.nan(K) == TRUE){browser()}
    a1 <- sum(1/seq(1:(K-1))) #get a1
    ws.theta <- n_seg/a1 #get ws.theta
    ws.thetaf <- ws.theta/ws #ws.theta fraction

    #get T's D by part. See original paper, Tajima 1989.
    
    a2 <- sum(1/(seq(1:(K-1))^2))
    b1 <- (K+1)/(3*(K-1))
    b2 <- (2*(K^2 + K + 3))/(9*K*(K-1))
    c1 <- b1 - (1/a1)
    c2 <- b2 - ((K+2)/(a1*K)) + (a2/(a1^2))
    e1 <- c1/a1
    e2 <- c2/(a1^2 + a2)
    D <- (ts.theta - ws.theta)/sqrt((e1*n_seg) + e2*n_seg*(n_seg - 1))
    
    #output result for this window, step to the next window
    
    out[i,"group"] <- x[1,3]
    if("pop" %in% colnames(x)){
      out[i,"pop"] = x[1,"pop"] #if a pop column is in the input, add a pop column here.
    }
    out[i,"position"] <- c
    out[i,"ws.theta"] <- ws.thetaf
    out[i,"ts.theta"] <- ts.thetaf
    out[i,"D"] <- D
    c <- c + step*1000
    i <- i + 1
  }
  return(out)
}


#function to run any command after spliting the data by group and by population. Group and population must be in columns 3 and 4, respectively.
#The first argument of the function to run must be the data provided to that function.
run_gp <- function(x, FUN, ...){
  w_df<- data.frame()
  pops <- unique(x[,4])
  groups <- unique(x[,3]) 
  for (i in 1:length(groups)){
    w_data <- x[x[,3] == groups[i],]
    print(groups[i])
    for(j in 1:length(pops)){
      print(pops[j])
      x_data <- FUN(w_data[w_data[,4] == pops[j],], ...)
      if(length(x_data) == 0){next}
      if(is.data.frame(x_data)){
        x_data <- cbind(pop = pops[j], x_data)
        w_df <- rbind(w_df, x_data)
      }
      else{
        if(is.data.frame(w_df)){
          w_df <- numeric(0)
        }
        w_df <- c(w_df,x_data)
      }
    }
  }
  if(!is.data.frame(w_df)){
    warning("Output is vector: if merging with meta info, check that info data frame is sorted identically:
            Group, Pop, then snp by snp!")
  }
  return(w_df)
}

#Function to test for deviation from HWE via perumtaiton (given small cell entries) for each snp.
#inputs: x: data, snps in subsequent rows
#        num_start: column index with the first column of data
#        miss: CHARACTER entry which encodes missing data
#        test: the test to be used.
#            options: exact: Exact test, from Wigginton et al 2005.
#                     permutation: permutation chi.square test.
#        n.reps: number of permutations to get p.values. Needed if test = "permutation".
HWE <- function(x, num_start, miss, test = "exact", n.reps = 20000){
  
  
  if(test == "exact"){
    #edited from Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
    #Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76: 000 - 000
    #code available at http://csg.sph.umich.edu/abecasis/Exact/snp_hwe.r
    SNPHWE <- function(x){
      obs_homr <- x[1]
      obs_homc <- x[2]
      obs_hets <- x[3]
      if (obs_homr < 0 || obs_homc < 0 || obs_hets < 0){
        return(-1.0)
      }
      # total number of genotypes
      N <- obs_homr + obs_homc + obs_hets

      # number of rare allele copies
      rare  <- obs_homr * 2 + obs_hets
      
      # Initialize probability array
      probs <- rep(0, 1 + rare)
      
      # Find midpoint of the distribution
      mid <- floor(rare * ( 2 * N - rare) / (2 * N))
      if ( (mid %% 2) != (rare %% 2) ) mid <- mid + 1
      
      probs[mid + 1] <- 1.0
      mysum <- 1.0
      
      # Calculate probablities from midpoint down 
      curr_hets <- mid
      curr_homr <- (rare - mid) / 2
      curr_homc <- N - curr_hets - curr_homr
      
      while ( curr_hets >=  2)
      {
        #equation 2
        probs[curr_hets - 1]  <- probs[curr_hets + 1] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0)  * (curr_homc + 1.0))
        mysum <- mysum + probs[curr_hets - 1]
        
        # 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
        curr_hets <- curr_hets - 2
        curr_homr <- curr_homr + 1
        curr_homc <- curr_homc + 1
      }    
      
      # Calculate probabilities from midpoint up
      curr_hets <- mid
      curr_homr <- (rare - mid) / 2
      curr_homc <- N - curr_hets - curr_homr
      
      while ( curr_hets <= rare - 2)
      {
        probs[curr_hets + 3] <- probs[curr_hets + 1] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
        mysum <- mysum + probs[curr_hets + 3]
        
        # add 2 heterozygotes -> subtract 1 rare homozygtote, 1 common homozygote
        curr_hets <- curr_hets + 2
        curr_homr <- curr_homr - 1
        curr_homc <- curr_homc - 1
      }    
      
      # P-value calculation
      target <- probs[obs_hets + 1]
      
      #plo <- min(1.0, sum(probs[1:obs_hets + 1]) / mysum)
      
      #phi <- min(1.0, sum(probs[obs_hets + 1: rare + 1]) / mysum)
      
      # This assignment is the last statement in the fuction to ensure 
      # that it is used as the return value
      p <- min(1.0, sum(probs[probs <= target])/ mysum)
    }
  }
  
  #figure out snp format
  if(nchar(as.character(x[1,(num_start)])) == 2){
    snp_form <- 2
  }
  else if (nchar(as.character(x[1,(num_start)])) == 4){
    snp_form <- 4
  }
  else{
    stop("Genotypes must be in either 2 (NN) or 4 (0000) format.")
  }
  
  #functions to get allele and genotype frequencies minus missing data.
  g_row <- function(row){
    t <- table(unlist(row))
    t <- t[names(t) != miss]
    #fill table if missing genotypes
    if(length(t) == 2){
      as <- unique(c(substr(names(t),1,snp_form/2), substr(names(t),(snp_form/2) + 1, snp_form)))
      if(!(paste0(as[1], as[1]) %in% names(t))){
        missing.geno <- paste0(as[1], as[1])
      }
      else if (!(paste0(as[2], as[2]) %in% names(t))){
        missing.geno <- paste0(as[2], as[2])
      }
      else if(!(paste0(as[1], as[2]) %in% names(t)) & !(paste0(as[2], as[1]) %in% names(t))){
        missing.geno <- paste0(as[1], as[2])
      }
      t <- c(t, 0)
      t <- setNames(t, c(names(t)[1:2], missing.geno))
    }
    else if (length(t) == 1){
      t <- c(t,0,0)
      names(t) <- c(names(t)[1], paste0(rep("X", snp_form/2), collapse = ""), paste0(rep("X", snp_form/2), collapse = ""))
    }
    #organize: homs first then het, with minor hom first
    as <- c(rbind(substr(names(t),1,snp_form/2), substr(names(t),(snp_form/2) + 1, snp_form))) #break down the table names into componenets via interleaving
    g1 <- ifelse(as[1] == as[2], "hom", "het") #determine if first entry is a hom or het
    g2 <- ifelse(as[3] == as[4], "hom", "het")
    g3 <- ifelse(as[5] == as[6], "hom", "het")
    gs <- setNames(c(g1,g2,g3), 1:3)
    gs <- rev(sort(gs)) #sorted by name. actual indices are stored as names
    t <- t[as.numeric(names(gs))] #reorder table, homs then hets
    hom.fs <- sort(c(t[1], t[2])) #sort the homs by frequency
    t <- c(hom.fs, t[3]) #rebind the table
    return(t)
  }
  
  #allele frequencies from g_row output converted to frequencies
  cat("Getting genotype counts...\n")
  drs <- x[,(num_start:ncol(x))] #get just the data
  gfs <- t(apply(drs,1,g_row)) #get genotypes
  
  if(test == "permutation"){
    p.vals <- numeric(nrow(gfs))
    gfs.f <- gfs/rowSums(gfs)
    cat("Doing permutation tests.\nGetting allele frequencies...\n")
    afs <- cbind(gfs.f[,1] + (gfs.f[,3]/2), gfs.f[,2] + (gfs.f[,3]/2)) #allele frequencies from genotype frequencies
  
    #get the expected genotype counts
    gfs.e <- matrix(NA, nrow = nrow(gfs), ncol = 3)
    gfs.e[,1] <- rowSums(gfs)*(afs[,1]^2) #expected minor homozygote
    gfs.e[,2] <- rowSums(gfs)*(afs[,2]^2) #expected major homozygote
    gfs.e[,3] <- rowSums(gfs)*(2*afs[,1]*afs[,2]) #expected heterozygote
    gfs.test <- cbind(gfs, gfs.e)
    
    #print(gfs.test)
    #print(afs)
    #do chi-squared tests
    cat("Running chi-squared tests:\n")
    for(i in 1:nrow(gfs)){
      if(i %% 1000 == 0){cat("snp number: ", i, "\n")}
      if(any(afs[i,] == 0)){
        p.vals[i] <- 1
      }
      else{
        p.vals[i] <- chisq.test(cbind(gfs.test[i,1:3], gfs.test[i,4:6]), simulate.p.value = TRUE, B = n.reps)$p.value
        if(is.na(p.vals[i])){browser()}
      }
    }
  }
  else if (test == "exact"){
    cat("Running exact tests...\n")
    p.vals <- unlist(apply(gfs, 1, SNPHWE))
  }
  return(p.vals)
}



#Calculates pairwise fst for each pair of populations according to hohenlohe (2010). Input format is
#snp, position, group, pop, total allele count (total), number of alleles (num), ni1, ni2, and pi (diversity, from calc_pi).
#Automatically sorts data by group, position, and population.
calc_pairwise_Fst <- function(data, do.nk = FALSE, skip.FST = FALSE){
  if(!do.nk & skip.FST){
    stop("Must specify either pairwise FST, nk, or both")
  }
  require(plyr)
  if(!all(data[,6] == 2)){
    vio <- which(data[,6] != 2)
    cat("Some loci do not have two alleles. Violating loci:\n")
    print(vio)
    stop()
  }
  i <- 1
  colnames(data) <- c("snp", "position", "group", "pop", "ntot", "num", "ni1", "ni2", "pi")
  data[,"pop"] <- as.character(data[,"pop"])
  data[,"group"] <- as.character(data[,"group"])
  data <- arrange(data, group, position, pop)
  pops <- sort(unique(data[,"pop"]))
  out <- as.data.frame(matrix(NA, ncol = 3+(length(pops)*(length(pops) - 1)/2), nrow = nrow(data)/length(pops)))
  
  #print(out)
  #initialize pop comparison columns.
  comps <- c()
  while (i < (length(pops))){
    j <- 1 + i
    for (j in j:length(pops)){
      comps <- c(comps, paste0(pops[i], "_", pops[j]))
      j <- j + 1
    }
    i <- i + 1
  }
  i <- 1
  colnames(out) <- c("snp", "position", "group", comps)
  out[,1:3] <- data[data$pop == pops[1],1:3] #add snp, position, group to output
  if(do.nk){ #prepare nk output if requested
    nout <- out
    colnames(nout)[4:ncol(nout)] <- paste0("nk_", comps)
  }
  
  #print(out)
  
  #loop through each comparison and caculate pairwise FST at each site
  c.col <- 4 #set starting column for pasting data
  for (i in 1:(length(pops) - 1)){ #i is the first pop
    idat <- data[data$pop == pops[i],] #get data for first pop
    j = i + 1 #intialize j as the next pop
    for (j in j:length(pops)){#j is pop being compared
      if(!skip.FST){
        jdat <- data[data$pop == pops[j],] #get data for second pop
        #n.ali <- ifelse(idat[,7] != 0, 1, 0) + ifelse(idat[,8] != 0, 1, 0) #get number of alleles in pop 1
        #n.alj <- ifelse(jdat[,7] != 0, 1, 0) + ifelse(jdat[,8] != 0, 1, 0) #get number of alleles in pop 2
        #com.top <- (choose(n.ali, 2) * idat[,9]) + (choose(n.alj, 2) * jdat[,9]) #get the numerator
        com.top <- (choose(idat[,5], 2) * idat[,9]) + (choose(jdat[,5], 2) * jdat[,9])
        
        t.ni1 <- idat[,7] + jdat[,7] #get the total number of allele one
        t.ni2 <- idat[,8] + jdat[,8] #get the total number of allele two
        ptop <- choose(t.ni1, 2) + choose(t.ni2, 2) #get the pooled pi numerator
        pbot <- choose((t.ni1 + t.ni2), 2) #get the pooled pi denominator
        ppi <- 1 - ptop/pbot #get pooled pi
        #com.bot <- ppi * (choose(n.ali,2) + choose(n.alj,2)) #get the denominator
        com.bot <- ppi * (choose(idat[,5],2) + choose(jdat[,5],2)) #get the denominator
        Fst <- 1 - com.top/com.bot #get fst
        if(any(abs(Fst) >= 1, na.rm = T)){browser()}
        #Fst[t.ni1 == 0 | t.ni2 == 0] <- 0 #could uncomment this if want 0s instead of NaNs.
        out[,c.col] <- Fst #write fst
      }
      if(do.nk){
        nout[,c.col] <- idat[,5] + jdat[,5]
      }
      c.col <- c.col + 1 #agument c.col
    }
  }
  if(do.nk){
    if(skip.FST){
      return(nout)
    }
    else{
      return(list(FST = out, nk = nout))
    }
  }
  else if (!skip.FST){
    return(out)
  }
}

#wrapper for just doing pairwise.nk
calc_pairwise_nk <- function(data){
  out <- calc_pairwise_Fst(data, do.nk = TRUE, skip.FST = TRUE)
  return(out)
}

#Calculates observed heterozygosity at a snp.
#Inputs:  data: Input data. In numeric or NN character format, see format_snps output options 4 or 6.
#               Nucleotides coded as 01 to 04 or A,T,C,G, respectively.
#         num_start: Index of first column with genotypes. Columns before this will be added to output.
#         m.dat: Missing data indicator. Should be in the same format as data.
#         pop: List with population information for individuals. Format is as produced by:
#              list(c("North", "South", "East", "West"), c(10,20,30,40)). First vector is names of pops,
#              second vector is the count of each pop. Input data MUST be in the same order as this list.
calc_Ho <- function(data, num_start = 4, m.dat = "0000", pop = NULL){
  
  #set possible heterozygotes
  if(nchar(data[1,num_start]) == 4 & nchar(m.dat) == 4){
    hl <- c(m.dat, "0101", "0202", "0303", "0404")
  }
  else if (nchar(data[1,num_start]) == 2 & nchar(m.dat) == 2){
    hl <- c(m.dat, "AA", "TT", "CC", "GG")
  }
  else{
    stop("Data and missing signifier must be in either one (NN) or two character (0000) per allele format.")
  }
  
  #initalize output
  if(!is.list(pop)){
    pop = list("pi", ncol(data) - num_start + 1)
  }
  pns <- pop[[1]]
  psz <- pop[[2]]
  out <- matrix(NA, nrow(data), length(pns))
  out <- cbind(data[,1:(num_start - 1)], out)
  colnames(out) <- c(colnames(data)[1:(num_start - 1)], pns) #set output column names
  
  
  #do each pop, need to loop here.
  
  c.col <- num_start
  for (j in 1:length(pns)){
    wdat <- data[,c.col:(c.col+psz[j] - 1)]
    #with this data, figure out heterozygosity
    het.c <- rowSums(ifelse(wdat == hl[1] | wdat == hl[2] 
                            | wdat == hl[3] | wdat == hl[4] 
                            | wdat == hl[5], 0, 1))
    ho <- het.c/rowSums(ifelse(wdat == m.dat, 0, 1))
    out[,num_start + j - 1] <- ho
    c.col <- c.col + psz[j]
  }
  return(out)
}


#Prepares a heatmap from pariwse LD data. 
#Inputs:  x: Data, format like that given by LD_full_pairwise rsq and Dprime outputs.
#            However, first column can contain snp names rather than first column of data. Otherwise,
#            Data will be reformated in this way.
#         r: Region of the chromosome to subset and plot.
#            Given in mb in the format numeric vector c(lower, upper).
#         l.text: Legend title.
#         colors: Colors to use, character vector format c("lower", "upper").
#         title: Plot title
#         t.sizes: Text sizes, numeric vector format c(title, legend.title, legend.ticks, axis, axis.ticks)
LD_pairwise_heatmap <- function(x, r = NULL, 
                                l.text = "rsq", colors = c("white", "black"), 
                                title = "Pairwise LD", t.sizes = c(16, 13, 10, 12, 10), 
                                background = "white"){
  
  #Created in part by Nick Sard
  ### loading packages 
  require(ggplot2)
  require(reshape2)

  #remove columns and rows with no data
  x <- x[!apply(x, 1, function(y)all(is.na(y))), !apply(x, 2, function(y)all(is.na(y)))]
  
  #if first column doesn't contain positions, take the row names and set as first column.
  if(any(is.na(x[,1]))){
    x <- cbind(V1 = row.names(x), x)
  }
  
  
  #melting the df so its computer reable and fixing the names of the columns
  heatmap_x <- melt(x, id.vars = "V1")
  names(heatmap_x) <- c("SNPa", "SNPb", "value")

  #getting rid of all the zeros from snps being compared to themselves
  heatmap_x$SNPa <- as.numeric(as.character(heatmap_x$SNPa))
  heatmap_x$SNPb <- as.numeric(as.character(heatmap_x$SNPb))
  heatmap_x <- heatmap_x[!(heatmap_x$SNPa == heatmap_x$SNPb),]

  #removing NA values, since no Dups.  Note, also grabbing only section of interest
  heatmap_x <- heatmap_x[!is.na(heatmap_x$value),]
  
  #remove any other NAs
  heatmap_x <- heatmap_x[!is.na(heatmap_x$SNPa) & !is.na(heatmap_x$SNPb),]
  
  #get site names
  ms <- unique(c(colnames(x), x[,1]))
  ms <- c(x[1,1], colnames(x)[-1])
  ms <- as.numeric(as.character(ms))
  
  #subset down to the desired r if requested
  if(!is.null(r)){
    r <- r*1000000
    heatmap_x <- heatmap_x[heatmap_x$SNPa >= r[1] & heatmap_x$SNPa <= r[2] &
                           heatmap_x$SNPb >= r[1] & heatmap_x$SNPb <= r[2],]
    ms <- ms[ms >= r[1] & ms <= r[2]]
  }
  
  #finish messing with site names
  ms <- ms/1000000
  ms <- as.factor(ms)
  ms <- sort(ms)
  
  #set site positions to mb, convert to factor
  heatmap_x$SNPa <- heatmap_x$SNPa/1000000
  heatmap_x$SNPb <- heatmap_x$SNPb/1000000
  heatmap_x$SNPa <- as.factor(heatmap_x$SNPa)
  heatmap_x$SNPb <- as.factor(heatmap_x$SNPb)
  
  #reordering based on factors
  heatmap_x[["SNPa"]]<-factor(heatmap_x[["SNPa"]],levels= ms,ordered=T)
  heatmap_x[["SNPb"]]<-factor(heatmap_x[["SNPb"]],levels=rev(ms),ordered=T)
  
  #the plot
  out <- ggplot(heatmap_x, aes(x = SNPa, y=SNPb, fill=value))+
    geom_tile(color = "white")+
    scale_fill_gradient(low = colors[1], high = colors[2]) +
    theme_bw()+
    labs(x = "",y="", fill=l.text)+
    theme(legend.title= element_text(size = t.sizes[2]),
          axis.text = element_text(size = t.sizes[5]), 
          panel.grid.major = element_line(color = background),
          plot.title = element_text(size = t.sizes[1], hjust = 0.5), 
          axis.title = element_text(size = t.sizes[4]),
          legend.text = element_text(size = t.sizes[3]),
          panel.background = element_rect(fill = background, colour = background)) +
    scale_x_discrete(breaks = levels(heatmap_x$SNPa)[c(T, rep(F, 20))], label = abbreviate) +
    scale_y_discrete(breaks = levels(heatmap_x$SNPb)[c(T, rep(F, 20))], label = abbreviate) +
    ggtitle(title) + ylab("Position (Mb)") + xlab("Position (Mb)")
  return(out)
}


#Checks for private alleles.
#inputs: x: data, in allele count format such as that given by format_snps option 1. Expects
#           columns named "pop", "ni1", and "ni2", which contain pop designations, allele one counts,
#           and allele 2 counts.
check_private <- function(x){
  l <- unique(x$pop) #gets unique pops, requires column named pop
  
  a1m <- matrix(NA, nrow(x)/length(l), length(l)) #initialize a1 storage
  a2m <- matrix(NA, nrow(x)/length(l), length(l)) #initialize a2 storage
  nloc <- nrow(x)/length(l)
  
  #loop through pops load a1s and a2s
  count <- 1
  for(i in 1:length(l)){
    a1m[,i] <- x$ni1[count:(count + nloc - 1)]
    a2m[,i] <- x$ni2[count:(count + nloc - 1)]
    count <- count + nloc
  }
  
  #convert to presence absence
  a1m <- ifelse(a1m != 0, 1, 0)
  a2m <- ifelse(a2m != 0, 1, 0)
  
  #convert to private/not private
  a1m <- ifelse(rowSums(a1m) == 1 & a1m == 1, 1, 0)
  a2m <- ifelse(rowSums(a2m) == 1 & a2m == 1, 1, 0)
  
  #combine a1 and a2
  pa <- a1m + a2m
  
  #return data
  colnames(pa) <- paste0("pa_", l)
  return(as.data.frame(pa))
}



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
resample_long <- function(data, statistic, boots, sigma, nk_weight = FALSE, fws = NULL, n_snps = T, report = 10000){
  
  #report a few things, intialize others
  sig <- 1000*sigma #correct sigma
  num_snps <- nrow(data)
  #print(count_dist)
  smoothed_dist <- numeric(boots)
  cat("Sigma:", sigma, "\nTotal number of input snps: ", num_snps, "\nNumber of boots:", boots, "\n")
  
  
  #draw random centriods, their positions, lgs of those centroids, and number of snps in window if provided
  if(is.null(fws)){
    csnps <- sample(1:nrow(data), boots, replace = T) #get random centroids
    cs <- data$position[csnps] #get centroid positions
    grps <- data$group[csnps] #get centroid groups
    if(n_snps){ #if n_snps is specified, draw all of the random snps to fill in around centroids on windows now.
      nrands <- sample(1:nrow(data),sum(data$n_snps[csnps]), replace = T) 
      nrprog <- 1
    }
  }
  else{ #same deal, but for fixed slide windows from provided window information
    cwin <- sample(1:nrow(fws), boots, replace = T)
    cs <- fws$position[cwin]
    grps <- fws$group[cwin]
    if(n_snps){
      nrands <- sample(1:nrow(data),sum(fws$n_snps[cwin]), replace = T)
      nrprog <- 1
    }
  }
  
  #do bootstraps
  for (j in 1:boots){
    if(j %% report == 0){cat("Bootstrap number: ", j, "\n")}
    
    #figure out positions to use. Must be on the same group as the centriod and within 3*sigma
    tpos <- data$position[data$position >= cs[j] - 3*sig & 
                 data$position <= cs[j] + 3*sig & 
                 data$group == grps[j]]
    
    #if random snps haven't already been drawn...
    if(!n_snps){
      #draw random snps to fill those positions from ALL snps, with replacement  
      rdraws <- sample(1:nrow(data), length(tpos), replace = T)
    }
    else{
      #pull the correct randomly chosen snps from the vector of random snps.
      rdraws <- nrands[nrprog:(nrprog + length(tpos) - 1)]
      nrprog <- nrprog + length(tpos)
    }
    
    
    #get random stats and nks
    rs <- data[rdraws, statistic] #sample stats for the randomly selected snps
    if (nk_weight == TRUE){
      rnk <- data$nk[rdraws] #sample nks for the randomly selected snps
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

resample_long_par <- function(data, num_cores, statistic, boots, sigma, nk_weight = FALSE, fws = NULL){
  require(doParallel)
  sig <- 1000*sigma
  num_snps <- nrow(data)
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
      samp <- sample(1:nrow(data), 1)
      cs <- data$position[samp]
      grps <- data$group[samp]
    }
    else{
      samp <- sample(1:nrow(fws), 1)
      cs <- fws$position[samp]
      grps <- fws$group[samp]
    }
    
    #figure out positions to use
    tpos <- data$position[data$position >= cs - 3*sig & 
                            data$position <= cs + 3*sig & 
                            data$group == grps]
    
    
    #draw random snps to fill those positions, with replacement  
    rdraws <- sample(1:nrow(data), length(tpos), replace = T)
    
    #get random stats and nks
    rs <- data[rdraws, statistic] #sample random stats from the data
    if (nk_weight == TRUE){
      rnk <- data$nk[rdraws] #sample random nks from the data
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





#Filters the data according to a set of filters. Vectorized to handle large datasets.
#Arguments:
# data: input data frame
#
#FILTERS:
# non_poly: filter out snps with only one homozygous genotype observed, list return is npf
# maf: filter out snps with a minor allele frequency below that given, list return is maff
# hf_hets: filters out data with a heterozygote frequency above the given value, list return is hfhf
# min_ind: filters out snps genotyped at less than the given number of individuals, list return is mif
#This version is vectorized, and is much faster!
filter_snps <- function(data, ecs, non_poly = FALSE, bi_al = TRUE, maf = FALSE, pop = FALSE,
                        hf_hets = FALSE, min_ind = FALSE, mDat = "NN"){
  require(matrixStats)
  
  ##############################################################################################
  ###set up, get values used later, clean up data a bit, set any functions used lower
  cat("Initializing...\n")
  
  #seperate out the genotypes alone
  x <- data[,(ecs+1):ncol(data)]
  
  #get headers
  headers <- data[,1:ecs]
  remove(data)
  
  #get a single vector of genotypes
  xv <- as.vector(t(x))
  
  #determine snp format
  snp_form <- nchar(xv[1])
  
  #get number of samples
  nsamp <- ncol(x)
  
  #create tables of genotype counts for each locus.
  #function to do this, takes the pattern to match and the data, which is as a SINGLE VECTOR, by rows.
  #needs the number of samples from global function environment
  count_genos <- function(pattern, x){
    XXs <- grep(pattern, x) #figure out which entries match pattern
    out <- table(ceiling(XXs/nsamp)) #devide the entries that match by the number of samps to get which samp they come from (when rounded up), then table that.
    return(out)
  }
  
  #get all possible genotypes
  gs <- unique(xv)
  
  #which gs are heterozygous?
  hs <- substr(gs,1,snp_form/2) != substr(gs, (snp_form/2 + 1), snp_form*2)
  
  #which genotype is the missing data?
  mpos <- which(gs == mDat)
  
  #what are the possible alleles at all loci?
  as <- unique(c(substr(gs,1,snp_form/2), substr(gs, (snp_form/2 + 1), snp_form*2)))
  as <- as[as != substr(mDat, 1, snp_form/2)] #that aren't N?
  
  #############################################################################################3
  ###get a table of genotype and allele counts at each locus.
  cat("Creating genotype table...\n")
  
  #for each element of gs, get the tables of genotype counts and add them to a matrix
  gmat <- matrix(0, nrow(x), length(gs)) #initialize matrix
  colnames(gmat) <- gs #set the matrix names
  
  #fill the matrix, one possible genotype at a time (hard enough to vectorize as it is).
  for(i in 1:length(gs)){
    tab <- count_genos(gs[i], xv)
    gmat[as.numeric(names(tab)),i] <- as.numeric(tab)
  }
  
  if(length(mpos) > 0){
    tmat <- gmat[,-c(mpos)] #gmat without missing data
  }
  else{
    tmat <- gmat
  }
  
  #get matrix of allele counts
  #initialize
  cat("Getting allele table...\n")
  amat <- matrix(0, nrow(gmat), length(as))
  colnames(amat) <- as
  
  #fill in
  for(i in 1:length(as)){
    b <- grep(as[i], colnames(tmat))
    hom <- which(colnames(tmat) == paste0(as[i], as[i]))
    if(length(hom) == 0){
      het <- b
      amat[,i] <- rowSums(tmat[,het])
    }
    else{
      het <- b[b != hom]
      if(length(het) > 0){
        if(is.matrix(tmat[,het])){
          amat[,i] <- (tmat[,hom] * 2) + rowSums(tmat[,het])
        }
        else{
          amat[,i] <- (tmat[,hom] * 2) + tmat[,het]
        }
      }
      else{
        amat[,i] <- (tmat[,hom] * 2)
      }
    }
  }

  
  ##############################################################################################
  ###do filters
  keep <- logical(nrow(x)) #vector to track status
  
  #===============================================
  ##non-biallelic loci
  if(bi_al){
    cat("Filtering non-biallelic loci...\n")
    bimat <- ifelse(amat, TRUE, FALSE)
    bi <- ifelse(rowSums(bimat) > 2, T, F) #if true, should keep the allele
    #some tests require this, so subset the matrices and redefine things if true and some are multi-allelic
    if((maf | hf_hets) &  sum(bi) != 0){
      x <- x[bi,]
      xv <-   xv <- as.vector(t(x))
      gs <- unique(xv)
      hs <- substr(gs,1,snp_form/2) != substr(gs, (snp_form/2 + 1), snp_form*2)
      mpos <- which(gs == mDat)
      as <- unique(c(substr(gs,1,snp_form/2), substr(gs, (snp_form/2 + 1), snp_form*2)))
      as <- as[as != substr(mDat, 1, snp_form/2)] #that aren't N?
      gmat <- gmat[bi,]
      tmat <- tmat[bi,]
      amat <- amat[bi,]
      headers <- headers[bi,]
    }
    else{
      keep <- keep + bi
    }
  }
  
  #===============================================
  ##non-poly
  if(non_poly){
    cat("Filtering non-polymorphic loci...\n")
    npmat <- ifelse(tmat == 0, F, T) #convert to logical
    gcounts <- rowSums(npmat) #how many genotypes observed?
    hgcounts <- as.logical(rowSums(npmat[,hs[-mpos]])) #are hets observed?
    gcounts <- ifelse(gcounts == 1, F, T) #are there more than one genotype observed?
    np <- gcounts + hgcounts #figure out if either or both of the above are true
    np <- !as.logical(np) #convert to logical (if at least one is true, this will be false and we should keep the snp)
    keep <- keep + np
  }
  
  #===============================================
  ##min inds
  if(min_ind){
    cat("Filtering loci sequenced in few individuals...\n")
    mi <- ncol(x) - gmat[,colnames(gmat) == mDat]
    mi <- !mi >= min_ind #if false, enough samples, so keep.
    keep <- keep + mi
  }
  
  
  #===============================================
  ##minor allele frequency, both total and by pop. Should only run if bi_al = TRUE.
  if(maf){
    #if not filtering with multiple pops
    if(!is.list(pop)){
      cat("Filtering low minor allele frequencies, no pops...\n")
      mafs <- 1 - rowMaxs(amat)/rowSums(amat)
      mafs <- mafs < maf #less than required, set to true and reject.
      mafs[is.na(mafs)] <- TRUE
      keep <- keep + mafs
    }
    else{
      cat("Filtering low minor allele frequencies, pop:\n")
      pmafs <- logical(nrow(x))
      for(i in 1:length(pop[[1]])){
        cat(pop[[1]][i], "\n")
        #re-establish matrices with each pop
        
        #set the input data
        if(i == 1){
          popx <- x[,1:pop[[2]][i]]
        }
        else{
          popx <- x[,(sum(pop[[2]][1:(i-1)]) + 1):(sum(pop[[2]][1:i]))]
        }
        
        #get a single vector of genotypes
        popxv <- as.vector(t(popx))
        
        #get number of samples
        nsamp <- ncol(popx)
        
        #get all possible genotypes
        popgs <- unique(popxv)
        
        #which gs are heterozygous?
        pophs <- substr(popgs,1,snp_form/2) != substr(popgs, (snp_form/2 + 1), snp_form*2)
        
        #which genotype is the missing data?
        popmpos <- which(popgs == mDat)
        
        #what are the possible alleles at all loci?
        popas <- unique(c(substr(popgs,1,snp_form/2), substr(popgs, (snp_form/2 + 1), snp_form*2)))
        popas <- popas[popas != substr(mDat, 1, snp_form/2)] #that aren't N?
        
        #for each element of gs, get the tables of genotype counts and add them to a matrix
        popgmat <- matrix(0, nrow(popx), length(popgs)) #initialize matrix
        colnames(popgmat) <- popgs #set the matrix names
        
        #fill the matrix, one possible genotype at a time (hard enough to vectorize as it is).
        for(j in 1:length(popgs)){
          poptab <- count_genos(popgs[j], popxv)
          popgmat[as.numeric(names(poptab)),j] <- as.numeric(poptab)
        }
        
        if(length(popmpos) > 0){
          poptmat <- popgmat[,-c(popmpos)] #gmat without missing data
        }
        else{
          poptmat <- popgmat
        }
        
        #get matrix of allele counts
        #initialize
        popamat <- matrix(0, nrow(popgmat), length(popas))
        colnames(popamat) <- popas
        
        #fill in
        for(j in 1:length(popas)){
          b <- grep(popas[j], colnames(poptmat))
          hom <- which(colnames(poptmat) == paste0(popas[j], popas[j]))
          if(length(hom) == 0){
            het <- b
            popamat[,j] <- rowSums(poptmat[,het])
          }
          else{
            het <- b[b != hom]
            if(length(het) > 0){
              if(is.matrix(poptmat[,het])){
                popamat[,j] <- (poptmat[,hom] * 2) + rowSums(poptmat[,het])
              }
              else{
                popamat[,j] <- (poptmat[,hom] * 2) + poptmat[,het]
              }
            }
            else{
              popamat[,j] <- (poptmat[,hom] * 2)
            }
          }
        }
        
        #-------------------
        #established new set for this, now do maf as above.
        popmafs <- 1 - rowMaxs( popamat)/rowSums( popamat)
        popmafs <-  popmafs >= maf #if greater than requested, passes check for this pop
        popmafs[is.na(popmafs)] <- FALSE #call false, no information to go on here
        pmafs <- pmafs + popmafs #add this logical to the vector containing the sum of all such vectors
      }
      pmafs <- !as.logical(pmafs) #if false, we keep the allele since it was at >= maf in at least one pop
      keep <- keep + pmafs    
    }
  }
  
  #===============================================
  ##hf_hets. Should only run if bi_al = TRUE.
  if(hf_hets){
    cat("Filtering high frequency heterozygote loci...\n")
    #get heterozygote counts
    if(sum(hs[-mpos]) > 1){
      hetsum <- rowSums(tmat[,hs[-mpos]])
    }
    else if (sum(hs[-mpos]) == 0){
      hetsum <- numeric(nrow(x))
    }
    else{
      hetsum <- tmat[,hs[-mpos]]
    }
    hf <- hetsum/rowSums(tmat)
    hf <- hf >= hf_hets #if false, heterozygote frequency is lower than cut-off, keep locus
    keep <- keep + hf
  }
  
  ##############################################################################################
  ###apply results
  cat("Applying results...\n")
  keep <- !as.logical(keep)
  x <- x[keep,]
  x <- cbind(headers[keep,], x)
  return(x)
}





#Converts snp/microsat data into different formats. Current supported formats are:
#Arguments:
# data: input data. Rows are loci, columns are individuals.
# ecs: number of extra columns at the header of the input (metadata, ect)
# output: output format:
#   1) per pop allele count format for Bayescan, ect. (v)
#   2) genepop format (v)
#   3) structure format (v)
#   4) 2 character numeric format (v)
#   5) migrate-n hapmap format (v)
#   6) NN format (from numeric) (v)
#   7) allele presence absensence format
# input_form: Input format
#   "NN": Default, alleles as single characters (as in AT or CG),
#   "0000": numeric. Will step to NN first before converting.
#   "msat_n": microsat, in two or three character format. 
#           If in two character format, use "msat_2". If in three, use "msat_3" 
#           Currently only supported for conversion to 8; 2, 3, and 4 forthcoming.
#   "snp_tab": Converts to NN first.
# miss: missing data format (per allele), typically "N", or "00"
# pop: For format 1 and 5. 
#      Number of samples in each pop as a list (as in list(c("POP1", "POP2"), c(48,50)))).
#      Samples must be in the order given in input data frame.
# n_samp: number of randomly selected snps, ect to take. Can also take a numeric vector containing SNP indices to keep. For option 3.
# l_names: Vector of locus names. For option 8.
# interp_miss: Should missing data be interpolated? T or F. For option 8.
#note: (v) under format options means that I've already vectorized it.
format_snps <- function(data, ecs, output = 1, input_form = "NN", 
                        miss = "N", pop = 1, n_samp = NA, l_names = NULL,
                        interp_miss = T){
  
  require(dplyr)
  #do checks, print info
  if(output == 1){
    cat("Output 1 selected, converting to per pop allele count format.\n")
    if(input_form != "0000" & input_form != "NN"){
      stop("Only 0000 and NN formats accepted.")
    }
    if(is.list(pop)){
      cat("Pops:")
      for (i in 1:length(unlist(pop[1]))){
        cat("\n\t", unlist(pop[1])[i], ", size: ", unlist(pop[2])[i])
      }
      cat("\n")
    }
  }
  else if(output == 2){
    cat("Output 2 selected, converting to genepop format.\n")
    if(input_form != "0000" & input_form != "NN" & input_form != "snp_tab"){
      stop("Only 0000, NN, and snp_tab formats accepted for now.")
    }
  }
  else if(output == 3){
    cat("Output 3 selected, converting to STRUCTURE format.\n")
    if(input_form != "0000" & input_form != "NN"){
      stop("Only 0000 and NN formats accepted for now.")
    }
    if(!is.integer(n_samp)){
      cat("Number of sub-samples to take:", "n_samp", "\n")
    }
    else if (!is.na(n_samp)){
      stop("Number of sub-samples to take must be an integer.\n")
    }
  }
  else if(output == 4){
    cat("Output 4 selected, converting to numeric 2 character format.\n")
    if(input_form != "NN" & input_form != "snp_tab"){
      stop("Only NN and snp_tab formats accepted for now.")
    }
  }
  else if(output == 5){
    cat("Output 5 selected, converting to migrate-N hapmap format.\n")
    if(input_form != "NN" & input_form != "0000"){
      stop("Only 0000 and NN formats accepted.")
    }
    if(is.list(pop)){
      cat("Pops:")
      for (i in 1:length(unlist(pop[1]))){
        cat("\n\t", unlist(pop[1])[i], ", size: ", unlist(pop[2])[i])
      }
      cat("\n")
    }
  }
  else if(output == 6){
    cat("Output 6 selected, converting to NN format.\n")
    if(input_form != "0000" & input_form != "snp_tab"){
      stop("Only 0000 and snp_tab formats accepted.")
    }
  }

  else if(output == 7){
    cat("Output 7 selected, converting to allele presence/absense format.\n")
  }
  else{
    stop("Please specify output format.")
  }
  
  if(input_form == "NN"){
    cat("Input format: NN\n")
    if(nchar(miss) != 1){
      stop("Missing data format must be one character.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "0000"){
    cat("Input format: 0000\n")
    if(nchar(miss) != 2){
      stop("Missing data format must be two characters.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "msat_2"){
    cat("Input format: msat, 2 character (0000)\n")
    if(nchar(miss) != 2){
      stop("Missing data format must be two characters.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "msat_3"){
    cat("Input format: msat, 3 character (000000)\n")
    if(nchar(miss) != 3){
      stop("Missing data format must be three characters.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "snp_tab"){
    cat("Input format, snp_tab.\n")
    if(nchar(miss) != 1){
      stop("Missing data format must be one character.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else{
    stop("Unsupported input format.")
  }
  
  
  #####################################################
  #function to create table of allele and genotype counts from all data. From filter_snps
  #input is a long vector of the data xv,
  #snp_form, the length of each genotype.
  #mDat, the missing data format
  #nsap, the number of samples
  tabulate_genotypes <- function(xv, snp_form, mDat, nsamp){
    #create tables of genotype counts for each locus.
    #function to do this, takes the pattern to match and the data, which is as a SINGLE VECTOR, by rows.
    #needs the number of samples from global function environment
    count_genos <- function(pattern, x){
      XXs <- grep(pattern, x) #figure out which entries match pattern
      out <- table(ceiling(XXs/nsamp)) #devide the entries that match by the number of samps to get which samp they come from (when rounded up), then table that.
      return(out)
    }
    
    #get all possible genotypes
    gs <- unique(xv)
    
    #which gs are heterozygous?
    hs <- substr(gs,1,snp_form/2) != substr(gs, (snp_form/2 + 1), snp_form*2)
    
    #which genotype is the missing data?
    mpos <- which(gs == mDat)
    
    #what are the possible alleles at all loci?
    as <- unique(c(substr(gs,1,snp_form/2), substr(gs, (snp_form/2 + 1), snp_form*2)))
    as <- as[as != substr(mDat, 1, snp_form/2)] #that aren't N?
    
    #############################################################################################3
    ###get a table of genotype and allele counts at each locus.
    cat("Creating genotype table...\n")
    
    #for each element of gs, get the tables of genotype counts and add them to a matrix
    gmat <- matrix(0, nrow(data), length(gs)) #initialize matrix
    colnames(gmat) <- gs #set the matrix names
    #fill the matrix, one possible genotype at a time (hard enough to vectorize as it is).
    for(i in 1:length(gs)){
      tab <- count_genos(gs[i], xv)
      gmat[as.numeric(names(tab)),i] <- as.numeric(tab)
    }
    
    if(length(mpos) > 0){
      tmat <- gmat[,-c(mpos)] #gmat without missing data
    }
    else{
      tmat <- gmat
    }
    
    #get matrix of allele counts
    #initialize
    cat("Getting allele table...\n")
    amat <- matrix(0, nrow(gmat), length(as))
    colnames(amat) <- as
    
    #fill in
    for(i in 1:length(as)){
      b <- grep(as[i], colnames(tmat))
      hom <- which(colnames(tmat) == paste0(as[i], as[i]))
      if(length(hom) == 0){
        het <- b
        amat[,i] <- rowSums(tmat[,het])
      }
      else{
        het <- b[b != hom]
        if(length(het) > 0){
          if(is.matrix(tmat[,het])){
            amat[,i] <- (tmat[,hom] * 2) + rowSums(tmat[,het])
          }
          else{
            amat[,i] <- (tmat[,hom] * 2) + tmat[,het]
          }
        }
        else{
          amat[,i] <- (tmat[,hom] * 2)
        }
      }
    }
    return(list(gs = tmat, as = amat))
  }
  
  #convert 0000 to NN unless output is 7 or 2. Return after if output is 6. (v)
  if(input_form == "0000" & output != 7 & output != 2){

    #Do conversion
    cat("Converting genotypes to NN form intermediary...")

    #vectorize and replace
    xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
    xv1 <- substr(xv, 1, 2)
    xv2 <- substr(xv, 3, 4)
    xv1[xv1 == "01"] <- "A"
    xv1[xv1 == "02"] <- "C"
    xv1[xv1 == "03"] <- "G"
    xv1[xv1 == "04"] <- "T"
    xv1[xv1 == miss] <- "N"
    xv2[xv2 == "01"] <- "A"
    xv2[xv2 == "02"] <- "C"
    xv2[xv2 == "03"] <- "G"
    xv2[xv2 == "04"] <- "T"
    xv2[xv2 == miss] <- "N"
    xv <- paste0(xv1, xv2)
    
    #rebind to matrix and remake data.
    xv <- matrix(xv, nrow(data), (ncol(data) - ecs), T)
    data <- cbind(data[,1:ecs], as.data.frame(xv))
    
    cat("\nMoving on to conversion...", "\n")
    miss <- "N" #reset miss to the correct entry
    if(output == 6){
      return(data) #all done if just converting to NN
    }
  }
  
  else if (input_form == "NN"){
    cat("Ensure that all data columns are character vectors!\n")
  }
  
  #convert snp_tab to NN (v)
  if(input_form == "snp_tab"){
    header <- data[,1:ecs]
    xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
    nsamp <- ncol(data) - ecs
    snames <- colnames(data)[(ecs+1):ncol(data)]
    rm(data)
    ptf <- nchar(xv) == 1
    xv[ptf] <- paste0(xv[ptf], xv[ptf]) #double up homozygotes
    remove(ptf)
    xv <- gsub(" ", "", xv) #combine heterozygotes.
    xv[xv == paste0(miss, miss)] <- "NN" #replace with correct missing data
    xv <- as.data.frame(matrix(xv, nrow(header), nsamp, byrow = T), stringsAsFactors = F)
    colnames(xv) <- snames
    out <- cbind(header, xv)
    if(output == 6){
      return(out)
    }
    else{
      data <- out
    }
  }
  
  
  ##convert to allele count or migrate-n format, migrate-n should ALWAYS have multiple pops (why else would you
  ##use it?) (v)
  if(output == 1 | output == 5){
    if(output == 5){cat("WARNING: Data does not have header or pop spacer rows.\n")}
    w_df <- data[,1:ecs] #intiallize w_df if doing option 1
    
    ##Make a function to sum each row
    sum_row <- function(row){
      l <- c(substr(row,1,1), substr(row,2,2))
      t <- table(l)
      return(t) #returns a table with counts of the two alleles and then an N
    }
    if (is.list(pop) == FALSE){
      if(output == 5){stop("Cannot convert to migrate-n format with only one pop.\n")}
      
      #prepare a genotype table
      xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
      tabs <- tabulate_genotypes(xv, nchar(xv[1]), paste0(miss, miss), ncol(data) - ecs)
      w_df$n_total <- rowSums(tabs$as)
      if(any(w_df$n_total != 2)){stop("More or less than two alleles detected at some loci.\n")}
      w_df$n_alleles <- rowSums(tabs$as != 0)
      #fill in ni1 and ni2.
      w_df$ni1 <- tabs$as[,"A"]
      w_df$ni1[w_df$ni1 == 0] <- tabs$as[which(w_df$ni1== 0), "G"]
      w_df$ni1[w_df$ni1 == 0] <- tabs$as[which(w_df$ni1== 0), "C"]
      w_df$ni2 <- tabs$as[,"T"]
      w_df$ni2[w_df$ni2 == 0] <- tabs$as[which(w_df$ni2== 0), "C"]
      w_df$ni2[w_df$ni2 == 0] <- tabs$as[which(w_df$ni2== 0), "G"]
    }
    else{
      require(plyr)
      pop_count <- length(unlist(pop[1]))
      pop_sizes <- unlist(pop[2])
      if(sum(unlist(pop[2])) != (ncol(data) - ecs)){#checks to make sure that the pop sizes add up to agree 
        #with data
        cat("Supplied population numbers do not equal the number of supplied loci columns. Exiting.")
        stop()
      }
      #initialize w_df for multiple pops
      #print(unlist(pop[1]))
      j <- 2
      w_df$pop <- unlist(pop[1])[1] #create pop column and set first pop name
      wa_df <- w_df #set first temp w_df
      for (j in j:pop_count){ #loop through each pop name
        wb_df <- wa_df #create temp df copy of wa_df
        #print(j)
        #print(unlist(pop[1])[j])
        wb_df$pop <- unlist(pop[1])[j]#change pop to current name
        w_df <- rbind(w_df, wb_df) #bind this copy to w_df
      }
      remove(wa_df, wb_df) #remove temp df
      
      ##loop through and count alleles for every row, splitting by pop
      
      #create allele count table for each pop, fill their section of data
      
      #build allele tables for each locus
      pop_as <- list()
      pals <- matrix(FALSE, nrow(data), 4)
      colnames(pals) <- c("A", "C", "G", "T")
      current <- ecs + 1
      for(i in 1:pop_count){
        #get the data for this population
        cat("Generating population allele count matrices, current pop:", pop[[1]][i], "\n")
        ne <- current + pop_sizes[i] - 1
        x <- data[,current:ne]
        xv <- as.vector(t(x))
        pop_as[[i]] <- tabulate_genotypes(xv, nchar(xv[1]), paste0(miss, miss), pop_sizes[i])$as
        pop_as[[i]] <- pop_as[[i]][, order(colnames(pop_as[[i]]))]
        #figure out which alleles are present at each locus to write these.
        pals <- pals + ifelse(pop_as[[i]] != 0, TRUE, FALSE)
        current <- current + pop_sizes[i]
        
      }
      
      #now write the correct counts
      cat("Tabulating and writing results...\n")
      pals <- ifelse(pals != 0, TRUE, FALSE)
      if(any(rowSums(pals) != 2)){
        stop("More or less than two alleles detected at some loci! Check/filter input data.\n")
      }
      #Convert into vector which says which elements to keep.
      palsv <- as.vector(t(pals))
      
      w_df$n_total <- NA
      w_df$n_alleles <- NA
      w_df$ni1 <- NA
      w_df$ni2 <- NA
      
      #loop through and print data
      for(i in 1:pop_count){
        #compare palsv to allele tables to keep correct alleles.
        av <- as.vector(t(pop_as[[i]]))
        av <- av[palsv]
        av <- matrix(av, nrow(data), 2, byrow = T)
        w_df[w_df$pop == pop[[1]][i],]$n_total <- rowSums(av)
        w_df[w_df$pop == pop[[1]][i],]$n_alleles <- rowSums(ifelse(av != 0, TRUE, FALSE))
        w_df[w_df$pop == pop[[1]][i],]$ni1 <- av[,1]
        w_df[w_df$pop == pop[[1]][i],]$ni2 <- av[,2]
      }
      
     
      row.names(w_df) <- 1:nrow(w_df) #fix row names
      if(output == 5){ #if doing output style 5...
        w_df <- w_df[,c(1, (ecs + 1):ncol(w_df))] #remove extra columns except the first and pop columns
      }
    }
    return(w_df)
  }
  
  
  
  ##convert to genepop or numeric format (v)
  if (output == 2 | output == 4){
    cat("WARNING: For output option 3, output does not have population seperators or header information.", "\n", "Converting genotypes...")
    
    if(input_form == "0000"){
      if(output == 2){
        w_df <- as.data.frame(t(data[,(ecs + 1):ncol(data)])) #remove extra columns and transpose data
        #print(names(w_df))
        row.names(w_df) <- paste0(row.names(w_df), " ,") #adding space and comma to row names, as required.
        return(w_df)
      }
      else{
        stop("Please select different input and output formats.")
      }
    }
    
    #vectorize and replace
    xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
    xv1 <- substr(xv, 1, 1)
    xv2 <- substr(xv, 2, 2)
    xv1[xv1 == "A"] <- "01"
    xv1[xv1 == "C"] <- "02"
    xv1[xv1 == "G"] <- "03"
    xv1[xv1 == "T"] <- "04"
    xv1[xv1 == miss] <- "00"
    xv2[xv2 == "A"] <- "01"
    xv2[xv2 == "C"] <- "02"
    xv2[xv2 == "G"] <- "03"
    xv2[xv2 == "T"] <- "04"
    xv2[xv2 == miss] <- "00"
    xv <- paste0(xv1, xv2)
    
    xv <- matrix(xv, nrow(data), (ncol(data) - ecs), T)
    data <- cbind(data[,1:ecs], as.data.frame(xv))
    
    cat("\n", "Cleaning up...", "\n")
    if(output == 2){ #convert to genepop
      w_df <- as.data.frame(t(data[,(ecs + 1):ncol(data)])) #remove extra columns and transpose data
      #print(names(w_df))
      row.names(w_df) <- paste0(row.names(w_df), " ,") #adding space and comma to row names, as required.
      return(w_df)
    }
    else {#prepare numeric output, otherwise same format
      return(data)
    }
  }
  
  
  ##convert to structure format (v)
  if (output == 3){
    
    #subset if requested
    if(!is.na(n_samp)){
      cat("Subsetting ")
      if(length(n_samp) > 1){
        cat("designated SNPs.\n")
        data <- data[,n_samp]
      }
      else{
        cat(n_samp, " random SNPs.\n")
        data <- data[sample(nrow(data), n_samp, T),]
      }
    }
    
    #plop actual data into two vectors, one for each allele.
    xv <- as.vector(t(t(data[,(ecs + 1):ncol(data)])))
    xv1 <- substr(xv, 1, 1)
    xv2 <- substr(xv, 2, 2)
    remove(xv)
    
    #convert to numeric
    xv1[xv1 == "A"] <- 1
    xv1[xv1 == "C"] <- 2
    xv1[xv1 == "G"] <- 3
    xv1[xv1 == "T"] <- 4
    xv1[xv1 == miss] <- 0
    xv2[xv2 == "A"] <- 1
    xv2[xv2 == "C"] <- 2
    xv2[xv2 == "G"] <- 3
    xv2[xv2 == "T"] <- 4
    xv2[xv2 == miss] <- 0
    
    #bind back into matrices
    xv1 <- matrix(xv1, ncol(data) - ecs, nrow(data), T)
    xv2 <- matrix(xv2, ncol(data) - ecs, nrow(data), T)
    
    #create output matrix
    outm <- matrix(NA, 2*(ncol(data) - ecs), nrow(data))
    
    #fill
    outm[seq(1,nrow(outm),2),] <- xv1
    outm[seq(2,nrow(outm),2),] <- xv2
    
    #add sample names
    snames <- character(nrow(outm))
    snames[seq(1,nrow(outm),2)] <- colnames(data[,(ecs+1):ncol(data)])
    snames[seq(2,nrow(outm),2)] <- colnames(data[,(ecs+1):ncol(data)])
    out <- cbind(ind = snames, as.data.frame(outm))
    
    return(out)
  }
  
  #presence/absence format
  if(output == 7){
    x <- data[,(ecs+1):ncol(data)] #get just data
    if((is.null(l_names)) | length(l_names) != nrow(data) | !(is.vector(l_names))){
      stop("Locus names required for output format 8. Please give a vector of names (l_names) equal to number of loci.")
    }
    if(input_form == "NN"){
      asize <- 1
    }
    else if (input_form == "0000" | input_form == "msat_2"){
      asize <- 2
    }
    else if (input_form == "msat_3"){
      asize <- 3
    }
    
    #write a function to quickly get allele list per locus
    acount <- function(x){
      a <- unique(c(substr(x, 1, asize), substr(x, (asize + 1), (2*asize))))
      return(sort(a))
    }
    
    #function to get average allele frequencies if needed.
    afreqs <- function(x){
      tab <- table(c(substr(x, 1, asize), substr(x, (asize + 1), (2*asize))))
      tab <- tab[names(tab) != miss]
      tab <- tab/sum(tab)
      tab <- tab[sort(names(tab))]
      return(tab)
    }
    
    #loop through each row, figure out number of alleles, create and fill matrix with presence/absence, bind matrix to existing.
    omat <- matrix(NA, ncol(x), 1)
    for (i in 1:nrow(x)){
      if(i %% 1000 == 0){cat("snp number: ", i, "\n")}
      alist <- acount(x[i,])
      wmat <- matrix(0, ncol(x), length(alist)) #output matrix. Individuals in rows, alleles in columns.
      colnames(wmat) <- paste0(l_names[i], "_", alist)
      a1s <- substr(x[i,], 1, asize) #get first allele in each ind
      a2s <- substr(x[i,], asize + 1, 2*asize) #get second allele in each ind.
      a1s <- match(a1s, alist) #get column indexes to add a presence to
      a2s <- match(a2s, alist) #same
      for(j in 1:nrow(wmat)){ #for each individual, add alleles to the appropriate bins 
        wmat[j, a1s[j]] <- wmat[j, a1s[j]] + 1
        wmat[j, a2s[j]] <- wmat[j, a2s[j]] + 1
      }
      if(interp_miss){ #interpolate missing genotypes if requested
        afs <- afreqs(x[i,]) #get overall allele frequencies
        ms <- which(wmat[,grepl(paste0("_", miss), colnames(wmat))] > 0) #figure out which had missing data
        wmat <- wmat[,-grepl(paste0("_", miss), colnames(wmat))] #remove the missing data column
        wmat[ms,] <- matrix(afs, length(ms), ncol(wmat), byrow = T) #add the interpolated allele frequencies.
      }
      
      omat <- cbind(omat, wmat) #bind to the output
    }
    w_df <- as.data.frame(omat[,-1]) #remove int column, convert to data frame.
    w_df <- cbind(samp = colnames(x), w_df) #add sample names.
    if(!interp_miss){cat("Finished. Warning: Missing data counts are also stored!\n")}
    return(w_df)
  }
}

#Calculates Dprime, rsq, and a p-value for LD for each pair of snps.
#inputs: x: data, in either numeric or NN form. Must contain a column named "position", can contain a column named "group" and/or "pop"
#        ecs: number of extra columns (headers) before the data in x
#        prox_table: Should a proximity table be output?
#        matrix_out: Should LD matrices be created?
#        mDat: What is the missing data character? Expects marker for a single allele. ("N" or "01")
LD_full_pairwise <- function(x, ecs, prox_table = TRUE, matrix_out = TRUE, mDat = "N") {
  #new LD pairwise function. Loops through each loci, but vectorized for each comparison within that loci
  require("dplyr")
  require("reshape2")
  require("matrixStats")
  
  #get data and headers
  headers <- x[,1:ecs]
  x <- x[,(ecs + 1):ncol(x)]
  x <- as.matrix(x)
  
  #get missing genotype
  dmDat <- paste0(mDat, mDat)
  
  #data format
  sform <- nchar(x[1,1])/2
  
  #get unique alleles present at each locus (note: this made a much quicker a tab... should implement this elsewhere...)
  p1 <- substr(as.matrix(x), 1, sform)
  p2 <- substr(as.matrix(x), sform + 1, sform*2)
  pc <- melt(cbind(p1,p2))
  #amat <- with(pc, table(Var1, value))
  #amat <- amat[,colnames(amat) != mDat]
  as <- sort(unique(pc$value))
  as <- as.character(as[as != mDat])
  
  
  #need to loop through each loci and compare to everything else. Probably can't really vectorize the outer loop.
  #goal: make a haplotype table, where each row is a comparison and each column is a haplotype count
  #to count haplotypes: Double heterozgote (AC CG) = mark neither. 
  #                     Double homozygote (AA CC): mark two of the combination (A with C)
  #                     homo/het (AA CG): mark one of each combination (A with C and A with G)
  #
  #function to do this need to: 1) paste together the observed genotypes of all observed genotype combinations.
  #                             2) convert this to a table of counts of these genotypes for each pairwise combo.
  #                             3) clean the table
  #                             4) get haplotype counts.
  #
  #make a function to generate this table given a starting locus:
  # inputs: x: row containing genotypes at starting locus
  #         y: rows containing genotypes at all comparison loci
  tabulate_haplotypes <- function(x, y){
    library(reshape2)
    #1)
    #get the observed genotype combinations
    yv <- as.vector(t(y))
    gcv <- paste0(x, yv)
    if(!is.matrix(y)){
      gcv <- matrix(gcv, 1, length(x), byrow = T)
    }
    else{
      gcv <- matrix(gcv, nrow(y), length(x), byrow = T)
    }
    
    
    #2)
    #turn this into a genotype count table
    mgcv <- melt(gcv)
    ghapmat <- with(mgcv, table(Var1, value))
    
    #3) clean the table
    ##grab column names
    gl <- colnames(ghapmat)
    ##remove anything with missing data and double hets
    rgcs <- c(grep(paste0("^", dmDat), gl), #missing first locus
              grep(paste0(dmDat, "$"), gl), #missing second locus
              which(substr(gl, 1, sform) != substr(gl, (sform + 1), (sform *2)) &
                    substr(gl, (sform*2) + 1, sform*3) != substr(gl, (sform*3+1), sform*4))) #double het
    
    ##remove any double heterozygotes
    ghapmat <- ghapmat[,-rgcs]
    
    #add a filler row for the last pairwise comparison to make life easier.
    if(!is.matrix(y)){
      if(length(ghapmat) > 1){ #stop it from doing this if there is data for only one haplotype.
        ghapmat <- rbind(ghapmat, rep(c(10,0), 100)[1:length(ghapmat)])
      }
    }
    
    #if nothing remains, return nothing
    if(length(ghapmat) == 0){
      return(NA)
    }
    else if(!is.matrix(ghapmat)){ #if only one column remains...
      if(substr(gl[-rgcs], 1,sform) == substr(gl[-rgcs], sform + 1, sform*2)){ #double hom
        return(NA)
      }
      else{
        return(c(D = 0))
      }
    }
    
    #4) get hap table. Use the rules above to do this. Possible conditions, where alleles at locus 1 = 1a1b and locus 2 = 2a2b
    if(!is.vector(colnames(ghapmat))){browser()}
    ghapmat <- ghapmat[,order(colnames(ghapmat))] #put in order, just in case
    gl <- colnames(ghapmat) #get column names again
    hnames <- c(paste0(as, as[1]),
                paste0(as, as[2]),
                paste0(as, as[3]),
                paste0(as, as[4]))
    if(any(grepl("NA", hnames))){ #in the case of a missing allele...
      hnames <- hnames[-grep("NA",hnames)]
    }
    hnames <- sort(hnames)
    hapmat <- matrix(0, nrow(ghapmat), length(hnames))#initialize. all possible haplotypes
    colnames(hapmat) <- hnames
    
    
    ##figure out which are homozygotes and heterozygotes at either locus in the pairwise comparison.
    dhom <- substr(gl, 1, sform) == substr(gl, sform + 1, sform*2) &
            substr(gl, (sform*2) + 1, sform*3) == substr(gl, (sform*3+1), sform*4) #columns with double homozygotes
    dhom <- ghapmat[,dhom]
    het_l1 <- substr(gl, 1, sform) != substr(gl, sform + 1, sform*2) #columns where the first locus is het
    het_l1 <- ghapmat[,het_l1]
    het_l2 <- substr(gl, (sform*2) + 1, sform*3) != substr(gl, (sform*3+1), sform*4) #colunms where the second locus is het
    het_l2 <- ghapmat[,het_l2]
    
    #fix wierd cases where one of these isn't a matrix because only one haplotype falls into the category.
    if(any(!is.matrix(dhom), !is.matrix(het_l1), !is.matrix(het_l2))){
      if(!is.matrix(dhom)){
        dhom <- as.matrix(dhom)
        colnames(dhom) <- colnames(ghapmat)[substr(gl, 1, sform) == substr(gl, sform + 1, sform*2) &
          substr(gl, (sform*2) + 1, sform*3) == substr(gl, (sform*3+1), sform*4)] #columns with double homozygotes
      }
      if(!is.matrix(het_l1)){
        het_l1 <- as.matrix(het_l1)
        colnames(het_l1) <- colnames(ghapmat)[substr(gl, 1, sform) != substr(gl, sform + 1, sform*2)] #columns where the first locus is het
      }
      if(!is.matrix(het_l2)){
        het_l2 <- as.matrix(het_l2)
        colnames(het_l2) <- colnames(ghapmat)[substr(gl, (sform*2) + 1, sform*3) != substr(gl, (sform*3+1), sform*4)]
      }
    }
    
    #count up the haplotypes.
    #function to correct haplotype input matrix:
    GtoH <- function(x, n){
      m1 <- matrix(as.numeric(x), nrow(x), ncol(x))
      colnames(m1) <- colnames(x)
      m1 <- cbind(as.data.frame(t(m1)), n)
      m2 <- m1 %>% group_by(n) %>% summarise_all(funs(sum))
      m2 <- t(as.matrix(m2[,-1]))
      colnames(m2) <- sort(unique(n))
      return(m2)
    }
    ##homozygotes:
    hapmat[,colnames(hapmat) %in% paste0(substr(colnames(dhom), 1, sform),
                                       substr(colnames(dhom),(sform*2)+1,sform*3))] <- dhom*2
    ##heterozyogote locus 1
    n1 <- paste0(substr(colnames(het_l1), 1, sform),
                 substr(colnames(het_l1),(sform*2)+1,sform*3))
    n1 <- GtoH(het_l1, n1)
    n2 <- paste0(substr(colnames(het_l1),sform+1, sform*2),
                 substr(colnames(het_l1),(sform*3)+1, sform*4))
    n2 <- GtoH(het_l1, n2)
    hapmat[,colnames(hapmat) %in% colnames(n1)] <- n1 + hapmat[,colnames(hapmat) %in% colnames(n1)]
    hapmat[,colnames(hapmat) %in% colnames(n2)] <- n2 + hapmat[,colnames(hapmat) %in% colnames(n2)]
    
    ##heterozyogote locus 2
    n1 <- paste0(substr(colnames(het_l2), 1, sform),
                 substr(colnames(het_l2),(sform*2)+1,sform*3))
    n1 <- GtoH(het_l2, n1)
    n2 <- paste0(substr(colnames(het_l2),sform+1, sform*2),
                 substr(colnames(het_l2),(sform*3)+1, sform*4))
    n2 <- GtoH(het_l2, n2)
    hapmat[,colnames(hapmat) %in% colnames(n1)] <- n1 + hapmat[,colnames(hapmat) %in% colnames(n1)]
    hapmat[,colnames(hapmat) %in% colnames(n2)] <- n2 + hapmat[,colnames(hapmat) %in% colnames(n2)]
    
  
    #5)condense this hap table into the 1a2a, 1a2b, 1b2a, 1b2b format.
    # figure out how where haplotypes are missing. Note, do the case of two or three
    #missin haplotypes at the end.
    pmat <- ifelse(hapmat == 0, F, T)
    mmat <- pmat
    l1 <- substr(colnames(pmat), 1, sform)
    l2 <- substr(colnames(pmat), sform + 1, sform*2)
    mc <- 4 - rowSums(pmat)
    
    #function to see if haplotype is missing. x is the row index, m is a vector of the number missing haplotypes at each locus.
    cmhap <- function(x){
      out <- ifelse(rowSums(pmat[,l1 == l1[x]]) == 0 | rowSums(pmat[,l1 == l1[x]]) == 2, F,
                    ifelse(rowSums(pmat[,l2 == l2[x]]) > 0 & pmat[,x] != TRUE, T, F))
      return(out)
    }
    #fill the mmat.
    for(i in 1:ncol(pmat)){
      mmat[,i] <- cmhap(i)
    }
    
    #set the missing values in hapmat to NA, then replace those with zeros where there
    #are missing haplotypes.
    hapmat[hapmat == 0] <- NA
    hapmat[mmat == TRUE] <- 0
    
    
    
    
    
    
    #put in fillers when there are more than one haplotype is missing.
    pmat <- ifelse(is.na(hapmat), F, T)
    missing <- 4 - rowSums(pmat)
    m2 <- ifelse(missing >= 2, 0, NA)
    m3 <- ifelse(missing >= 3, 0, NA)
    m4 <- ifelse(missing == 4, 0, NA)
    mc <- cbind(m2,m2,m3,m4)
    
    #figure out which D, r values to give if two are missing...
    if(any(missing == 2)){
      m2mat <- hapmat[missing == 2,]
      m2matv <- as.vector(t(m2mat))
      if(is.matrix(m2mat)){
        m2matvcn <- rep(colnames(m2mat), nrow(m2mat))
      }
      else{
        m2matvcn <- rep(names(m2mat), 1)
      }
      m2matv[!is.na(m2matv)] <- m2matvcn[!is.na(m2matv)]
      m2matv <- na.omit(m2matv)
      if(is.matrix(m2mat)){
        m2mat <- matrix(m2matv, nrow(m2mat), 2, T)
        m2mat <- ifelse(substr(m2mat[,1], 1, sform) != substr(m2mat[,2], 1, sform) &
                          substr(m2mat[,1], sform + 1, sform*2) != substr(m2mat[,2], sform + 1, sform*2),
                        1,0)
      }
      else{
        m2mat <- ifelse(substr(m2matv[1], 1, sform) != substr(m2matv[2], 1, sform) &
                          substr(m2matv[1], sform + 1, sform*2) != substr(m2matv[2], sform + 1, sform*2),
                        1,0)
      }
      
    }
    else{
      m2mat <- "none"
    }
    
    hapmat <- cbind(hapmat, mc)
    hapmat <- as.vector(t(hapmat))
    hapmat <- na.omit(hapmat)
    #if(length(hapmat) %% nrow(y) != 0){
    #  browser()
    #}
    hapmat <- matrix(hapmat, nrow(ghapmat), 4, byrow = T)
    
    #now just have the haplotypes. These will calculate D in the case of 1 or 0 missing haplotypes.
    #when there are three missing haplotypes, D will be 0. When there are 2, D will be 0 or 1.
    return(list(hapmat = hapmat, missing = missing, m2 = m2mat))
  }
  
  
  #initialize output.
  if(prox_table){
    prox <- data.frame(p1 = numeric(1), p2 = numeric(1), rsq = numeric(1), Dprime = numeric(1), pval = numeric(1))
  }
  if(matrix_out){
    rmat <- matrix(NA, nrow(x) -1, nrow(x) - 1)
    colnames(rmat) <- headers$position[-1]
    rownames(rmat) <- headers$position[-length(headers$position)]
    Dpmat <- rmat
    pvmat <- rmat
  }
  if(!matrix_out & !prox_table){
    stop("Please specify output format.\n")
  }
  
  #run length prediction variables for progress reporting
  compfun <- function(x){
    return(((x-1)*x)/2)
  }
  totcomp <- compfun(nrow(x))
  cpercent <- 0
  
  #loop through and get haplotypes, calc LD for each locus.
  for(i in 1:(nrow(x) - 1)){
    cprog <- (totcomp - compfun(nrow(x) - i - 1))/totcomp
    if(cprog >= 0.05 + cpercent){
      cat("Progress:", paste0(round(cprog*100), "%."), "\n")
      cpercent <- cprog
    }
    haps <- tabulate_haplotypes(x[i,], x[(i+1):nrow(x),])
    
    #if we had only one haplotype or no haplotypes:
    if(is.na(haps[1])){
      if(prox_table){
        prox <- rbind(prox,
                      cbind(p1 = headers$position[i], p2 = headers$position[(i+1):nrow(x)],
                            rsq = NA, Dprime = NA, pval = NA))
      }
      if(matrix_out){
        #reminder: columns start at locus two, rows start at locus one (but end at nlocus - 1)
        rmat[i,] <- NA
        Dpmat[i,] <- NA
        pvmat[i,] <- NA
      }
      next()
    }
    if(length(haps) == 1){
      if(prox_table){
        prox <- rbind(prox,
                      cbind(p1 = headers$position[i], p2 = headers$position[(i+1):nrow(x)],
                            rsq = 0, Dprime = 0, pval = 0))
      }
      if(matrix_out){
        #reminder: columns start at locus two, rows start at locus one (but end at nlocus - 1)
        fill <- rep(NA, nrow(x) - (nrow(x[,i:nrow(x)]) - 1) - 1)
        Dprime <- c(fill, rep(0, length = (nrow(x[,i:nrow(x)]) - 1)))
        rsq <- c(fill, rep(0, length = (nrow(x[,i:nrow(x)]) - 1)))
        pval <- c(fill, rep(0, length = (nrow(x[,i:nrow(x)]) - 1)))
        rmat[i,] <- rsq
        Dpmat[i,] <- Dprime
        pvmat[i,] <- pval
      }
      next()
    }
    
    
    missing <- haps$missing
    m2 <- haps$m2
    haps <- haps$hapmat
    #A1B1 is col 1, A1B2 is col 2, A2B1 is col 3, A2B2 is col 4.
    
    #calc stats where >1 haps aren't missing
    A1B1f <- haps[,1]/rowSums(haps)
    A1B2f <- haps[,2]/rowSums(haps)
    A2B1f <- haps[,3]/rowSums(haps)
    A2B2f <- haps[,4]/rowSums(haps)
    A1f <- A1B1f + A1B2f
    A2f <- 1 - A1f
    B1f <- A1B1f + A2B1f
    B2f <- 1 - B1f
    D <- A1B1f - A1f*B1f
    D2 <- A2B1f - A2f*B1f
    Dprime <- ifelse(D > 0, D/rowMins(cbind(A1f*B2f, A2f*B1f)),
                   ifelse(D < 0, D/rowMaxs(cbind(-A1f*B1f, -A2f*B2f)),
                   0))
    rsq <- (D^2)/(A1f*A2f*B1f*B2f)
    
    #fix for when more missing haps.
    Dprime[missing == 3] <- 0
    Dprime[missing == 4] <- NA
    rsq[missing == 3 | missing == 4] <- NA
    if(length(m2) > 1){
      Dprime[missing == 2] <- m2
      rsq[missing == 2] <- ifelse(m2 == 1, 1, NA)
    }
    
    #get pvals
    chisqu <- rsq*(4)
    pval <- 1 - pchisq(chisqu, 1)
    
    #remove dummy filler if this was the final comparison.
    if(i == (nrow(x) - 1)){
      Dprime <- Dprime[-2]
      rsq <- rsq[-2]
      pval <- pval[-2]
    }
    
    #write output.
    if(prox_table){
      prox <- rbind(prox,
                    cbind(p1 = headers$position[i], p2 = headers$position[(i+1):nrow(x)],
                          rsq = rsq, Dprime = Dprime, pval = pval))
    }
    if(matrix_out){
      #reminder: columns start at locus two, rows start at locus one (but end at nlocus - 1)
      fill <- rep(NA, nrow(x) - length(Dprime) - 1)
      Dprime <- c(fill, Dprime)
      rsq <- c(fill, rsq)
      pval <- c(fill, pval)
      rmat[i,] <- rsq
      Dpmat[i,] <- Dprime
      pvmat[i,] <- pval
    }
  }
  if(prox_table){
    prox <- prox[-1,]
    if(any(colnames(headers) == "group")){
      prox$group <- headers[1,"group"]
    }
    if(any(colnames(headers) == "pop")){
      prox$pop <- headers[1,"pop"]
    }
    prox <- prox[,c(6:ncol(prox), 1:5)]
  }
  return(list(prox = prox, Dprime = Dpmat, rsq = rmat, pval = pvmat))
}


#Calls the LD_full_pairwise function for a file with multiple linkage groups. Returns a nested list of
#Dprime and rsq matrices. Requires the linkage group titled be named "group". Returns a concatenated prox_table.
Full_LD_w_groups <- function(x, ecs, prox_table = TRUE, matrix_out = TRUE, mDat = "N", report = 1){
  if (prox_table == FALSE & matrix_out == FALSE){
    stop("Please specify output format.")
  }
  w_list<- list()
  prox_out <- data.frame()
  groups <- unique(x$group)
  for (i in 1:length(groups)){
    w_data <- subset(x, group == groups[i])
    if(nrow(w_data) == 0 | nrow(w_data) == 1){
      next()
    }
    if(i %% report == 0){cat("Group #:", i, "of", length(groups), "Name:", groups[i], "\n")}
    w_data <- LD_full_pairwise(x = w_data, ecs = ecs, prox_table = prox_table, mDat = mDat, matrix_out = matrix_out)
    if(prox_table == TRUE){
      if(matrix_out == TRUE){
        prox_out <- rbind(prox_out, w_data$prox)
        w_list[[groups[i]]] <- w_data[1:2]
      }
      else{
        prox_out <- rbind(prox_out, w_data)
      }
    }
    else if (matrix_out == TRUE){
      w_list[[groups[i]]] <- w_data
    }
  }
  if (prox_table == TRUE){
    if(matrix_out == TRUE){
      w_list[["prox_table"]] <- prox_out
    }
    else{
      w_list <- prox_out
    }
  }
  return(w_list)
}


#runs LD_full_pairwise in parallel across many different linkage groups/scaffolds. Arguments are the same
#as LD_full_pairwise, save for num_cores (which is the number of cores to run on).
#If all outputs are desired, Returns a list where the first length(groups) elements are lists containing matrices for rsq and Dprime for
#the corresponding group. The last element is a prox table.
#currently giving an error, need to check it out in more detail.
Full_LD_g_par <- function(x, ecs, num_cores, prox_table = TRUE, matrix_out = TRUE, mDat = "N"){
  if (prox_table == FALSE & matrix_out == FALSE){
    stop("Please specify output format.")
  }
  require(doParallel)
  require(doSNOW)
  require(matrixStats)
  require(dplyr)
  require(reshape2)
  cl <- makeSOCKcluster(num_cores)
  registerDoSNOW(cl)
  
  
  groups <- unique(x$group)
  ntasks <- length(groups)
  progress <- function(n) cat(sprintf("Group %d out of",n), length(groups), "is complete.\n")
  opts <- list(progress=progress)
  
  output <- foreach(i = 1:length(groups), .export = 'LD_full_pairwise', .packages = c("dplyr", "reshape2", "matrixStats"), .inorder = TRUE,
                    .options.snow = opts) %dopar% {
                      w_data <- x[x$group == groups[i],]
                      cat("Initializing :", groups[i], ".\n")
                      if(nrow(w_data) == 0 | nrow(w_data) == 1){
                        next()
                      }
                      w_data <- LD_full_pairwise(w_data, ecs = ecs, prox_table = prox_table, matrix_out = matrix_out, mDat = mDat)
                    }
  #loop through output and put in easier to manage form
  if(prox_table){p <- c()}
  if(matrix_out){m <- vector("list", length(groups))}
  for(i in 1:length(groups)){
    if(prox_table){
      if(matrix_out){
        p <- rbind(p, output[[i]]$prox)
      }
      else{
        p <- rbind(p, output[[i]])
      }
    }
    if(matrix_out){
      m[[i]] <- output[[i]][-1]
      names(m) <- groups
    }
  }
  
  #clear cluster and return stuff
  stopCluster(cl)
  registerDoSEQ()
  if(prox_table){
    if(matrix_out){
      m[[i+1]] <- p
      names(m) <- c(groups, "prox")
      return(m)
    }
    else{
      return(p)
    }
  }
  else{
    return(m)
  }
}
