stickPED <- read.csv("inst/extdata/test_ped.csv", header = T, stringsAsFactors = F)

usethis::use_data(stickPED, overwrite = TRUE)
