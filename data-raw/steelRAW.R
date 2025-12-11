steelRAW <- read.table("inst/extdata/steelhead_msats.txt", header = FALSE, stringsAsFactors = F, colClasses = "character")
colnames(steelRAW) <- c("pop", paste0("msat_", 1:13))

usethis::use_data(steelRAW, overwrite = TRUE)

