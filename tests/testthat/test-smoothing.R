test_that("basic",{
  # smooth
  check <- calc_ho(stickSNPs, "pop")
  check <- calc_smoothed_averages(check, facets = c("chr.pop"), stats.type = "single", sigma = 250, step = 250)
  check <- get.snpR.stats(check, "chr.pop", "maf")$single.window
  
  # tests
  expect_true("maf" %in% colnames(check))
  expect_true(all(c("facet", "subfacet", "snp.facet", "snp.subfacet", "position") %in% colnames(check)))
  expect_equal(nrow(check), 2838)
  expect_true(all(check[check$subfacet == "ASP" & check$snp.subfacet == "groupVI",]$position %in% 
                     seq(0, 14000000, by = 250*1000)))
})


test_that("pairwise",{
  # smooth
  check <- calc_pairwise_fst(stickSNPs, "pop")
  check <- calc_smoothed_averages(check, facets = c("chr.pop"), stats.type = "pairwise", sigma = 250, step = 250)
  check <- get.snpR.stats(check, "chr.pop", "fst")$pairwise.window
  
  # tests
  expect_true("fst" %in% colnames(check))
  expect_true(all(c("facet", "subfacet", "snp.facet", "snp.subfacet", "position") %in% colnames(check)))
  expect_equal(nrow(check), 6904)
  expect_true(all(check[check$subfacet == "ASP" & check$snp.subfacet == "groupVI",]$position %in% 
                    seq(0, 14000000, by = 250*1000)))
})


test_that("no_chr",{
  # smooth
  check <- calc_pairwise_fst(stickSNPs, "pop")
  check <- calc_pi(check, "pop")
  check <- calc_smoothed_averages(check, facets = "pop", sigma = 250, step = 250)
  checkp <- get.snpR.stats(check, "pop", "fst")$pairwise.window
  checks <- get.snpR.stats(check, "pop", "pi")$single.window
  
  
  # tests
  expect_true(all(checkp$snp.facet == ".base" & checkp$snp.subfacet == ".base"))
  expect_true(all(checks$snp.facet == ".base" & checks$snp.subfacet == ".base"))
})