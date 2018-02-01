#Returns the average number of snps per window of size 6 sigma (in kb).
win_ave_SNPs <- function(data, sigma) {
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
  return(data)
}
