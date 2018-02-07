#data source: Hemstrom, W., van de Wetering, S., & Banks, M. A. (2017). Fish ladder installation across a historic barrier asymmetrically increased conspecific introgressive hybridization between wild winter and summer run steelhead salmon in the Siletz River, Oregon. Canadian Journal of Fisheries and Aquatic Sciences. doi:10.1139/cjfas-2016-0411

genos <- read.table("~/GitHub/tSNE_data/steelhead/msats/steelhead_genotypes.txt", sep='', header=FALSE, colClasses = "character")
genos <- t(genos)
pops <- genos[1,]
genos <- genos[-1,]
colnames(genos) <- pops
lnam <- c("Oki23", "Ssa407", "mSsa408", "Ots209", "OtsG249b", "OtsG85", "Omy27", "Omy1001", "Ots243", "Ots409", "OtsG3", "Ots212", "Omm1087")

sthMSATS <- as.data.frame(genos, stringsAsFactors = F)
sthMSATS <- cbind(lnum = 1:nrow(sthMSATS), locus = lnam, sthMSATS)

devtools::use_data(sthMSATS, overwrite = T)
