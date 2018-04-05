#convert stickSNPs into each output format for an example!

 pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
 l <- list(c(names(pops)), as.numeric(pops))
 f1 <- format_snps(stickSNPs, 3, "ac", pop = l)

 #option 2:
 f2 <- format_snps(stickSNPs, 3, "genepop")

 #option 3, subsetting out 100 random alleles:
 f3 <- format_snps(stickSNPs, 3, "structure", n_samp = 1:100, pop = l)

 #option 4:
 f4 <- format_snps(stickSNPs, 3, "numeric")

 #option 5:
 f5 <- format_snps(stickSNPs, 3, "hapmap", pop = l)

 #option 6:
 num <- format_snps(stickSNPs, 3, "numeric")
 f6 <- format_snps(num, 3, "character", input_form = "0000", miss = "00")

 #option 7, SNP data:
 f7 <- format_snps(stickSNPs, 3, "pa")

 #option 8:
 f8 <- format_snps(stickSNPs, 3, "rafm", pop = l)

 #option 9:
 f9 <- format_snps(stickSNPs, 3, "faststructure", pop = l)

stickFORMATs <- list(ac = f1 , genepop = f2, structure = f3, numeric = f4,
                     hapmap = f5, character = f6, pa = f7, rafm = f8, faststructure = f9)

devtools::use_data(stickFORMATs, overwrite = T)
