
# generate example format data
stickFORMATs <- list()

#ac
stickFORMATs[[1]] <- format_snps(stickSNPs, "ac", facets = "pop")

#genepop:
stickFORMATs[[length(stickFORMATs) + 1]] <- format_snps(stickSNPs, "genepop")

#STRUCTURE
stickFORMATs[[length(stickFORMATs) + 1]] <- format_snps(stickSNPs, "structure")

#fastSTRUCTURE
stickFORMATs[[length(stickFORMATs) + 1]] <- format_snps(stickSNPs, "faststructure")

#numeric:
stickFORMATs[[length(stickFORMATs) + 1]] <- format_snps(stickSNPs, "0000")

#hapmap for migrate-n:
stickFORMATs[[length(stickFORMATs) + 1]] <- format_snps(stickSNPs, "hapmap", facets = "pop")

#character:
stickFORMATs[[length(stickFORMATs) + 1]] <- format_snps(stickSNPs, "NN")

#presence/absence, SNP data:
stickFORMATs[[length(stickFORMATs) + 1]] <- format_snps(stickSNPs, "pa")

#RAFM, taking only 100 random snps and seperating by pop
stickFORMATs[[length(stickFORMATs) + 1]] <- format_snps(stickSNPs, "rafm", facets = "pop")

#dadi
## add ref and anc snp meta data columns to stickSNPs
dat <- as.data.frame(stickSNPs)
dat <- import.snpR.data(dat, snp.meta = cbind(ref = "ATA", anc = "ACT", stickSNPs@snp.meta), sample.meta = stickSNPs@sample.meta, mDat = stickSNPs@mDat)
stickFORMATs[[length(stickFORMATs) + 1]] <- format_snps(dat, "dadi", facets = "pop")

#PLINK! format
stickFORMATs[[length(stickFORMATs) + 1]] <- format_snps(stickSNPs, "plink")
#from command line, then run the snpR generated plink_out.sh to generate plink_out.bed.

#sn format, bernoulli interpolation
stickFORMATs[[length(stickFORMATs) + 1]] <- format_snps(stickSNPs, "sn")

#sequoia format
b <- stickSNPs@sample.meta
b$ID <- 1:nrow(b)
b$Sex <- rep(c("F", "M", "U", "no", "j"), length.out=nrow(b))
b$BirthYear <- round(runif(n = nrow(b), 1,1))

a <- stickSNPs
a@sample.meta <- b
a@sample.meta$newID <- paste0(a@sample.meta$pop, a@sample.meta$fam, a@sample.meta$ID)

stickFORMATs[[length(stickFORMATs) + 1]] <- format_snps(x=a, output = "sequoia", sample_id = "newID")



# name them
names(stickFORMATs) <- c("ac", "genepop", "structure", "faststructure", "0000", "hapmap", "NN", "pa", "rafm", "dadi", "plink", "sn", "sequoia")

# save
usethis::use_data(stickFORMATs, overwrite = TRUE)
