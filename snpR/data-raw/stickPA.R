stickPA <- format_snps(stickSNPs, 3, 7)
stickPA$pop <- substr(stickPA$samp,1,3)
stickPA <- stickPA[,c(1, ncol(stickPA), (2:(ncol(stickPA) - 1)))]
devtools::use_data(stickPA, overwrite = T)
