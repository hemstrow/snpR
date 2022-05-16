x <- stickRAW

stickSNPs <- import.snpR.data(x[,-c(1:2)],
                              x[,1:2],
                              data.frame(pop = substr(colnames(x)[-c(1:2)], 1, 3),
                                         fam = rep(c("A", "B", "C", "D"), length = ncol(x) - 2),
                                         stringsAsFactors = F),
                              "NN")

usethis::use_data(stickSNPs, overwrite = TRUE)
