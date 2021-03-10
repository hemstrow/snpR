x <- stickRAW

stickSNPs <- import.snpR.data(x[,-c(1:3)],
                              x[,2:3],
                              data.frame(pop = substr(colnames(x)[-c(1:3)], 1, 3),
                                         fam = rep(c("A", "B", "C", "D"), length = ncol(x) - 3),
                                         stringsAsFactors = F),
                              "NN")

usethis::use_data(stickSNPs, overwrite = TRUE)
