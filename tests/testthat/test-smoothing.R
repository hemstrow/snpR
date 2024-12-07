test_that("basic",{
  # smooth
  check <- calc_ho(stickSNPs, "pop")
  check <- calc_smoothed_averages(check, facets = c("chr.pop"), stats.type = "single", sigma = 250, step = 250, triple_sigma = TRUE)
  check <- get.snpR.stats(check, "chr.pop", "maf")$single.window
  
  # tests
  expect_true("maf" %in% colnames(check))
  expect_true(all(c("facet", "subfacet", "snp.facet", "snp.subfacet", "position") %in% colnames(check)))
  expect_equal(nrow(check), 2964)
  expect_true(all(check[check$subfacet == "ASP" & check$snp.subfacet == "groupVI",]$position %in% 
                     seq(0, 14000000 + 250*1000, by = 250*1000)))
})


test_that("pairwise",{
  # smooth
  check <- calc_pairwise_fst(stickSNPs, "pop")
  check <- calc_smoothed_averages(check, facets = c("chr.pop"), stats.type = "pairwise", sigma = 250, step = 250, triple_sigma = TRUE)
  check <- get.snpR.stats(check, "chr.pop", "fst")$pairwise.window
  
  # tests
  expect_true("fst" %in% colnames(check))
  expect_true(all(c("facet", "subfacet", "snp.facet", "snp.subfacet", "position") %in% colnames(check)))
  expect_equal(nrow(check), 7208)
  expect_true(all(check[check$subfacet == "ASP" & check$snp.subfacet == "groupVI",]$position %in% 
                    seq(0, 14000000 + 250*1000, by = 250*1000)))
})


test_that("no_chr",{
  # smooth
  check <- calc_pairwise_fst(stickSNPs, "pop")
  check <- calc_pi(check, "pop")
  check <- calc_smoothed_averages(check, facets = "pop", sigma = 250, step = 250, triple_sigma = TRUE)
  checkp <- get.snpR.stats(check, "pop", "fst")$pairwise.window
  checks <- get.snpR.stats(check, "pop", "pi")$single.window
  
  
  # tests
  expect_true(all(checkp$snp.facet == ".base" & checkp$snp.subfacet == ".base"))
  expect_true(all(checks$snp.facet == ".base" & checks$snp.subfacet == ".base"))
})

test_that("default bug",{
  check <- calc_ho(stickSNPs, "pop")
  check <- calc_smoothed_averages(check, facets = c("chr.pop"), stats.type = "single", sigma = 250, triple_sigma = TRUE)
  check <- get.snpR.stats(check, "chr.pop", "maf")$single.window
  expect_true(all(check$step == 500))
})

test_that("stats.type reassignment",{
  # single
  check <- calc_ho(stickSNPs)
  expect_error(check <- calc_smoothed_averages(check,  stats.type = "pairwise", sigma = 250, triple_sigma = TRUE),
               "No pairwise stats calculated")
  check <- calc_smoothed_averages(check, sigma = 250, triple_sigma = TRUE)
  expect_true("single.window" %in% names(get.snpR.stats(check, stats = "ho")))
})