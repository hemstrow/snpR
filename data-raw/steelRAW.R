steelRAW <- read.table("inst/extdata/steelhead_msats.txt", header = T, stringsAsFactors = F, colClasses = "character")

usethis::use_data(steelRAW, overwrite = TRUE)

