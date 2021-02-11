stickRAW <- read.table("inst/extdata/stick_NN_input.txt", header = T, stringsAsFactors = F)

usethis::use_data(stickRAW, overwrite = TRUE)

