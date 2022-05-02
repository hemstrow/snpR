context("Smoothing")

test_that("basic",{
  skip_on_cran()
  # smooth
  check <- calc_ho(stickSNPs, "pop")
  check <- calc_smoothed_averages(check, facets = c("chr.pop"), stats.type = "single", sigma = 250, step = 250)
  check <- get.snpR.stats(check, "chr.pop", "maf")$single.window
  
  # tests
  expect_true("maf" %in% colnames(check))
  expect_true(all(c("facet", "subfacet", "snp.facet", "snp.subfacet", "position") %in% colnames(check)))
  expect_equal(nrow(check), 9462)
  expect_equal(check[1,]$position + 250*1000, check[2,]$position)
})


test_that("pairwise",{
  skip_on_cran()
  # smooth
  check <- calc_pairwise_fst(stickSNPs, "pop")
  check <- calc_smoothed_averages(check, facets = c("chr.pop"), stats.type = "pairwise", sigma = 250, step = 250)
  check <- get.snpR.stats(check, "chr.pop", "fst")$pairwise.window
  
  # tests
  expect_true("fst" %in% colnames(check))
  expect_true(all(c("facet", "subfacet", "snp.facet", "snp.subfacet", "position") %in% colnames(check)))
  expect_equal(nrow(check), 23655)
  expect_equal(check[1,]$position + 250*1000, check[2,]$position)
})
