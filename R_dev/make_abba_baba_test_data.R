dat <- stickSNPs
dat <- add.facets.snpR.data(dat, "pop")


dat <- calc_maf(dat, "pop")
mafs <- get.snpR.stats(dat, "pop")

dat@snp.meta$ref <- get.snpR.stats(dat)$minor # note, this is done by default if these columns don't exist!
dat@snp.meta$anc <- get.snpR.stats(dat)$major # note, this is done by default if these columns don't exist!

test_data <- cbind(dat@facet.meta, dat@geno.tables$as)

test_data <- merge(test_data, dat@snp.meta, by = c("snp", "group", "position", ".snp.id"))
test_data <- test_data[-which(test_data$facet == ".base"),]

saveRDS(test_data, "R_dev/abba_baba_test_data.RDS")
