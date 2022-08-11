#=========maf=======
test_that("maf", {
  check <- filter_snps(.internal.data$test_snps, maf = 0.15, verbose = FALSE)
  stats <- get.snpR.stats(check, stats = "maf")
  
  expect_true(all(stats$single$maf >= 0.15))
  
  # comp to unfiltered
  bads <- get.snpR.stats(.internal.data$test_snps, stats = "maf")$single
  expect_true(all(bads$maf[!bads$position %in% stats$single$position] < 0.15))
})

test_that("maf_facets", {
  check <- filter_snps(.internal.data$test_snps, maf = 0.15, maf_facets = "pop", verbose = FALSE)
  check <- calc_maf(check, "pop")
  stats <- get.snpR.stats(check, "pop")
  stats$bad <- stats$maf < 0.15
  st <- tapply(stats$bad, stats[,"position"], sum)
  expect_true(all(st < 2))
  
  # comp to unfiltered
  bads <- calc_maf(.internal.data$test_snps, "pop")
  bads <- get.snpR.stats(bads, "pop")
  bads$bad <- bads$maf < 0.15
  bt <- tapply(bads$bad, bads[,"position"], sum)
  expect_true(!all(names(bt)[bt == 2] %in% names(st)))
})

#===========hf_hets===================
test_that("hf_hets", {
  # correct removed
  alld <- format_snps(.internal.data$test_snps, "sn", interpolate = FALSE)
  check <- filter_snps(.internal.data$test_snps, hf_hets = 0.4, verbose = FALSE)
  goods <- alld$position[which(rowSums(alld[, -c(1:2)] == 1, na.rm = T)/rowSums(!is.na(alld[,-c(1:2)])) <= 0.4)] # gets the position of loci with less than 40% hets
  expect_equal(goods, snp.meta(check)$position)
})

#===========hwe=======================
test_that("hwe", {
  # correct removed
  check <- filter_snps(.internal.data$test_snps, hwe = 0.5, verbose = FALSE)
  goods <- calc_hwe(.internal.data$test_snps)
  goods <- get.snpR.stats(goods)
  expect_true(all(get.snpR.stats(check)$position %in% goods$position[which(goods$pHWE > .5)]))
})


test_that("hwe_facets", {
  # correct removed
  check <- filter_snps(.internal.data$test_snps, hwe = 0.5, hwe_facets = "pop", verbose = FALSE)
  check <- calc_maf(check, "pop")
  goods <- calc_hwe(.internal.data$test_snps, facets = "pop")
  goods <- get.snpR.stats(goods, "pop")
  goods$bad <- goods$pHWE <= 0.5
  gt <- tapply(goods$bad, goods[,"position"], sum)
  expect_true(all(!names(gt)[gt != 0] %in% snp.meta(check)$position))
})

test_that("hwe_with_fwe", {
  # mostly just checking that everything runs without errors...
  
  # correct removed
  check <- filter_snps(.internal.data$test_snps, hwe = 0.5, fwe_method = "BH", verbose = FALSE)
  goods <- calc_hwe(.internal.data$test_snps, fwe_method = "BH")
  goods <- get.snpR.stats(goods)
  expect_true(all(get.snpR.stats(check)$position %in% goods$position[which(goods$pHWE_overall_BH > .5)]))
  
  # correct removed with facets
  check <- filter_snps(.internal.data$test_snps, hwe = 0.5, hwe_facets = "pop", fwe_method = "BH", verbose = FALSE)
  check <- calc_maf(check, "pop")
  goods <- calc_hwe(.internal.data$test_snps, facets = "pop", fwe_method = "BH")
  goods <- get.snpR.stats(goods, "pop")
  goods$bad <- goods$pHWE_byfacet_BH <= 0.5
  gt <- tapply(goods$bad, goods[,"position"], sum)
  expect_true(all(!names(gt)[gt != 0] %in% snp.meta(check)$position))
})

#===========min_ind======================
test_that("min_ind", {
  # correct removed
  check <- filter_snps(.internal.data$test_snps, min_ind = .9, verbose = FALSE)
  check <- snp.meta(check)$position
  comp <- rowSums(.internal.data$test_snps != "NN")/ncol(.internal.data$test_snps)
  expect_true(all(check %in% snp.meta(.internal.data$test_snps)$position[comp >= .9]))
})

#===========min_loci=====================
test_that("min_loci", {
  # correct removed
  expect_warning(check <- filter_snps(.internal.data$test_snps, min_loci = .9, re_run = FALSE, verbose = FALSE), "individuals were filtered out")
  check <- row.names(sample.meta(check))
  
  comp <- colSums(.internal.data$test_snps != "NN")/nrow(.internal.data$test_snps)
  comp <- row.names(sample.meta(.internal.data$test_snps)[comp >= .9,])
  
  expect_equal(comp, check)
})

#===========non_poly=====================
test_that("min_loci", {
  # correct removed
  td <- .internal.data$test_snps
  genotypes(td)[c(1, 5, 8),] <- rep("CC", ncol(td)) # add non-poly loci
  check <- filter_snps(td, verbose = FALSE)
  
  expect_true(all(!snp.meta(td)$position[c(1, 5, 8)] %in% snp.meta(check)$position))
})

#==========errors=========================
test_that("errors",{
  td <- .internal.data$test_snps[-c(1:2, 5, 8, 9:10)]
  expect_error(filter_snps(td, min_ind = .99, verbose = FALSE), "No loci passed filters.")
  td <- .internal.data$test_snps[,-c(7:8)]
  expect_error(filter_snps(td, min_loci = .99, non_poly = FALSE, verbose = FALSE), "No individuals passed filters.")
})

