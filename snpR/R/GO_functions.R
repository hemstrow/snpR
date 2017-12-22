#Takes an input file of genes with start/end position noted and maps back the average
#statitstic of choice for that gene from another provided data set with positions. Requires the "zoo"
#package
gene_ave_stat <- function(gene_data, stat_data, stat){
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


