stickSNPs <- read.table("../../../../Stickleback/Full Data/2017_reruns/snps_numeric_filt.txt", stringsAsFactors = F, header = T, colClasses = "character")
stickSNPs <- stickSNPs[grepl("group", stickSNPs$group),]

cols <- colnames(stickSNPs)

stickSNPs <- format_snps(stickSNPs, 3, 6, "0000", "00")
colnames(stickSNPs) <- cols

metadata <- stickSNPs[,1:3]
stickSNPs <- stickSNPs[,4:ncol(stickSNPs)]

stickSNPs <- cbind(metadata, stickSNPs[,order(colnames(stickSNPs))])

stickSNPs <- stickSNPs[seq(1, nrow(stickSNPs), by = 10),]

devtools::use_data(stickSNPs, overwrite = TRUE)
