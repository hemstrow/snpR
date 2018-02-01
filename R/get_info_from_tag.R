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
