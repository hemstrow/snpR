test_snps <- subset_snpR_data(stickSNPs, 1:12, sample(nsamps(stickSNPs), 10, F))
test_snps <- filter_snps(test_snps)

sample.meta(test_snps)$pop <- rep(c("ASP", "PAL"), 5)
sample.meta(test_snps)$fam <- rep(c("A", "B"), each = 5)


usethis::use_data(test_snps, internal = T, overwrite = T)
