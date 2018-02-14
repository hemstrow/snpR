#convert stickSNPs into each output format for an example!

 pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
 l <- list(c(names(pops)), as.numeric(pops))
 f1 <- format_snps(stickSNPs, 3, 1, pop = l)

 #option 2:
 f2 <- format_snps(stickSNPs, 3, 2)

 #option 3, subsetting out 100 random alleles:
 f3 <- format_snps(stickSNPs, 3, 3, n_samp = 1:100)

 #option 4:
 f4 <- format_snps(stickSNPs, 3, 4)

 #option 5:
 f5 <- format_snps(stickSNPs, 3, 5, pop = l)

 #option 6:
 num <- format_snps(stickSNPs, 3, 4)
 f6 <- format_snps(num, 3, 6, input_form = "0000", miss = "00")

 #option 7, SNP data:
 f7 <- format_snps(stickSNPs, 3, 7)

stickFORMATs <- list(o1 = f1 , o2 = f2, o3 = f3, o4 = f4, o5 = f5, o6 = f6, o7 = f7)

devtools::use_data(stickFORMATs, overwrite = T)
