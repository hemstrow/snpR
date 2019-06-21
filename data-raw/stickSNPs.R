x <- read.table("data/stick_NN_input.txt", stringsAsFactors = F, header = T, colClasses = "character")

stickSNPs <- import.snpR.data(x[,-c(1:3)],
                              x[,1:3],
                              data.frame(pop = substr(colnames(x)[-c(1:3)], 1, 3),
                                         fam = rep(c("A", "B", "C", "D"), length = ncol(x) - 3),
                                         stringsAsFactors = F),
                              "NN")

usethis::use_data(stickSNPs, overwrite = TRUE)
