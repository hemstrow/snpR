x <- steelRAW

steelMSATs <- read_non_biallelic(t(x[,-1]), sample.meta = data.frame(pop = x[,1],
                                              fam = rep(c("A", "B", "C", "D"), length = nrow(x)),
                                              stringsAsFactors = F), mDat = "0000")

usethis::use_data(steelMSATs, overwrite = TRUE)
