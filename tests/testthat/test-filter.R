#=========maf=======
test_that("maf", {
  check <- filter_snps(.internal.data$test_snps, maf = 0.15, verbose = FALSE)
  stats <- get.snpR.stats(check, stats = "maf")
  
  expect_true(all(stats$single$maf >= 0.15))
  
  # comp to unfiltered
  bads <- get.snpR.stats(.internal.data$test_snps, stats = "maf")$single
  expect_true(all(bads$maf[!bads$position %in% stats$single$position] < 0.15))
  
  # works with .base set
  check2 <- filter_snps(.internal.data$test_snps, maf = 0.15, maf_facets = ".base", verbose = FALSE)
  stats <- get.snpR.stats(check, stats = "maf")
  expect_identical(check, check2)
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
  expect_true(all(sort(get.snpR.stats(check)$position) == sort(goods$position[which(goods$pHWE > .5)])))
  
  # correct direction - he excess
  check_he <- filter_snps(.internal.data$test_snps, hwe = 0.5, verbose = FALSE, hwe_excess_side = "heterozygote")
  goods <- calc_hwe(.internal.data$test_snps)
  goods <- calc_fis(goods)
  goods <- get.snpR.stats(goods)
  expect_true(all(sort(get.snpR.stats(check_he)$position) == sort(goods$position[-which(goods$pHWE <= .5 & goods$fis <= 0)])))
  
  # correct direction - ho excess
  check_ho <- filter_snps(.internal.data$test_snps, hwe = 0.5, verbose = FALSE, hwe_excess_side = "homozygote")
  goods <- calc_hwe(.internal.data$test_snps)
  goods <- calc_fis(goods)
  goods <- get.snpR.stats(goods)
  expect_true(all(sort(get.snpR.stats(check_ho)$position) == sort(goods$position[-which(goods$pHWE <= .5 & goods$fis >= 0)])))
  
  # works with .base set
  check2 <- filter_snps(.internal.data$test_snps, hwe = 0.5, hwe_facets = ".base", verbose = FALSE)
  stats <- get.snpR.stats(check, stats = "maf")
  expect_identical(check, check2)
  
})


test_that("hwe_facets", {
  # correct removed
  check <- filter_snps(.internal.data$test_snps, hwe = 0.5, hwe_facets = "pop", verbose = FALSE)
  check <- calc_maf(check, "pop")
  goods <- calc_hwe(.internal.data$test_snps, facets = "pop")
  goods <- get.snpR.stats(goods, "pop")
  goods$bad <- goods$pHWE <= 0.5
  gt <- tapply(goods$bad, goods[,"position"], sum)
  expect_true(all(sort(as.numeric(names(gt)[gt == 0])) == sort(snp.meta(check)$position)))
  
  # correct direction - he excess
  check_he <- filter_snps(.internal.data$test_snps, hwe = 0.5, hwe_facets = "pop", verbose = FALSE, hwe_excess_side = "heterozygote")
  check_he <- calc_maf(check_he, "pop")
  goods <- calc_hwe(.internal.data$test_snps, facets = "pop")
  goods <- calc_fis(goods, "pop")
  goods <- get.snpR.stats(goods, "pop")
  goods$bad <- goods$pHWE <= 0.5
  goods$dir <- ifelse(goods$fis <= 0, 1, 0)
  if(any(is.na(goods$fis))){
    goods$dir[is.na(goods$fis)] <- 0
  }
  goods$bad <- goods$bad * goods$dir
  gt <- tapply(goods$bad, goods[,"position"], sum)
  expect_true(all(sort(as.numeric(names(gt)[gt == 0])) == sort(snp.meta(check_he)$position)))
  
  
  # correct direction - he excess
  check_ho <- filter_snps(.internal.data$test_snps, hwe = 0.5, hwe_facets = "pop", verbose = FALSE, hwe_excess_side = "homozygote")
  check_ho <- calc_maf(check_ho, "pop")
  goods <- calc_hwe(.internal.data$test_snps, facets = "pop")
  goods <- calc_fis(goods, "pop")
  goods <- get.snpR.stats(goods, "pop")
  goods$bad <- goods$pHWE <= 0.5
  goods$dir <- ifelse(goods$fis >= 0, 1, 0)
  if(any(is.na(goods$fis))){
    goods$dir[is.na(goods$fis)] <- 0
  }
  goods$bad <- goods$bad * goods$dir
  gt <- tapply(goods$bad, goods[,"position"], sum)
  expect_true(all(sort(as.numeric(names(gt)[gt == 0])) == sort(snp.meta(check_ho)$position)))
  
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
  check <- filter_snps(.internal.data$test_snps, min_loci = .9, re_run = FALSE, verbose = FALSE)
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

#==========singletons======================
test_that("singletons",{
  expect_warning(.make_it_quiet(td <- filter_snps(.internal.data$test_snps, singletons = TRUE)), "depriceated")
})

#==========mac=============================
test_that("mac",{
  expect_error(filter_snps(.internal.data$test_snps, maf = 0.05, mac = 1), "mac and maf cannot both be set")
  expect_error(filter_snps(.internal.data$test_snps, mac = 10), "mac must be greater than or equal to zero and less than the number of samples")
  expect_error(filter_snps(.internal.data$test_snps, mac = 5.5), "mac must be an integer")
  
  td <- filter_snps(.internal.data$test_snps, mac = 3, verbose = FALSE)
  
  bad.snps <- which(Matrix::rowSums(.internal.data$test_snps@geno.tables$as) - .rowMax_sparse(.internal.data$test_snps@geno.tables$as) <= 3)
  expect_equivalent(snp.meta(td), snp.meta(.internal.data$test_snps)[-bad.snps,])
  
  # with no non-poly removal
  test <- genotypes(stickSNPs)
  test <- rbind(test, rep("AA", ncol(test)))
  test <- import.snpR.data(test, sample.meta = sample.meta(stickSNPs))
  expect_true(nrow(filter_snps(test, non_poly = FALSE, mac = 1, verbose = FALSE)) == 101)
})


#==========mgc=============================
test_that("mgc",{
  expect_error(filter_snps(.internal.data$test_snps, mac = 1, mgc = 1), "mac and mgc cannot both be set")
  expect_error(filter_snps(.internal.data$test_snps, mgc = 5), "mgc must be greater than or")
  expect_error(filter_snps(.internal.data$test_snps, mgc = 3.5), "mgc must be an integer")
  
  td <- filter_snps(.internal.data$test_snps, mgc = 3, verbose = FALSE)
  
  hs <- colnames(.internal.data$test_snps@geno.tables$gs) %in% c("AC", "AG", "CT", "GT")
  hs_c <- Matrix::rowSums(.internal.data$test_snps@geno.tables$gs[, which(hs)])
  mg <- Matrix::rowSums(.internal.data$test_snps@geno.tables$gs[,-which(hs)]) - .rowMax_sparse(.internal.data$test_snps@geno.tables$gs[,-which(hs)])
  bad.snps <- which(hs_c + mg <= 3)
  expect_equivalent(snp.meta(td), snp.meta(.internal.data$test_snps)[-bad.snps,])
  
  # error that occured with less than 2 heterozygote options
  d <- format_snps(.internal.data$test_snps, output = "sn", interpolate = FALSE)
  d[is.na(d)] <- -9
  d <- import.snpR.data(d[,-c(1:2)], snp.meta = d[,1:2], mDat = -9)
  d <- filter_snps(d, mgc = 2, verbose = FALSE)
  expect_equal(nrow(d), 6)
  
  # with no non-poly removal
  test <- genotypes(stickSNPs)
  test <- rbind(test, rep("AA", ncol(test)))
  test <- import.snpR.data(test, sample.meta = sample.meta(stickSNPs))
  expect_true(nrow(filter_snps(test, non_poly = FALSE, mgc = 1, verbose = FALSE)) == 101)
})

#==========LD=============================
test_that("LD_pruning",{
  expect_error(filter_snps(stickSNPs, LD_prune_r = 1.1, LD_prune_sigma = 100), "0 and 1")
  expect_error(filter_snps(stickSNPs, LD_prune_r = .8, LD_prune_sigma = 100, LD_prune_method = "test"), "not recognized")
  expect_error(filter_snps(stickSNPs, LD_prune_facet = c("pop", "fam"), LD_prune_r = .8, LD_prune_sigma = 100), "Only one LD_prune_facet")
  expect_warning(filter_snps(stickSNPs, LD_prune_facet = "pop.chr", LD_prune_r = .8, LD_prune_sigma = 100, verbose = FALSE), "No valid LD loci pairs")
  
  td <- filter_snps(stickSNPs, LD_prune_facet = "pop.chr", LD_prune_r = .8, LD_prune_sigma = 1600, verbose = FALSE)
  expect_true(nrow(td) == 95) # should remove 5 SNPs
  
  td <- filter_snps(stickSNPs, LD_prune_facet = "pop.chr", LD_prune_r = .8, LD_prune_sigma = 1600, LD_prune_method = "Dprime", verbose = FALSE)
  expect_true(nrow(td) == 70) # should remove 30 SNPs
  
  # td <- filter_snps(stickSNPs, LD_prune_facet = "pop.chr", LD_prune_r = .8, LD_prune_sigma = 1600, LD_prune_method = "Dprime", LD_prune_use_ME = TRUE, verbose = FALSE)
  # expect_true(nrow(td) == 70) # should remove 29 SNPs
})

#==========garbage=========================
test_that("garbage",{
  expect_error(filter_snps(stickSNPs, min_ind = .5, remove_garbage = .8), "min_ind threshold should be higher")
  expect_error(filter_snps(stickSNPs, min_loci = .5, remove_garbage = .8), "min_loci threshold should be higher")
  expect_error(filter_snps(stickSNPs, remove_garbage = 1.2), "between 0 and 1")
  expect_error(filter_snps(stickSNPs, min_ind = .5, remove_garbage = ".8"), "must be a numeric")
  

  # for our test data this ends up not changing anything, so just check for the reports corresponding to correct removal.
  td <- capture.output(y <- filter_snps(stickSNPs, min_ind = .8, min_loci = .8, remove_garbage = .7))
  expect_true(any(grepl("Removing garbage individuals/loci", td)))
  expect_true(any(grepl("Removed 2 bad individuals", td)))
  expect_true(any(grepl("Removed 4 bad loci", td)))
  expect_true(any(grepl("Starting individuals: 98", td)))
  expect_true(any(grepl("Starting loci: 96", td)))
  
})

#==========filter order=========================
test_that("ind_first",{
  # check for the reports corresponding to correct removal order
  td <- capture.output(y <- filter_snps(stickSNPs, maf = .07, min_ind = .8, min_loci = .8))
  ind.call.line <- grep("^Filtering individuals. Starting individuals",td)
  loc.call.line <- grep("^Filtering loci. Starting loci",td)
  expect_true(ind.call.line > loc.call.line)
  expect_true(any(grepl("^Re-filtering loci", td)))
  
  # check for the reports corresponding to correct removal order
  td <- capture.output(y <- filter_snps(stickSNPs, maf = .07, min_ind = .8, min_loci = .8, inds_first = TRUE))
  ind.call.line <- grep("^Filtering individuals. Starting individuals",td)
  loc.call.line <- grep("^Filtering loci. Starting loci",td)
  expect_true(ind.call.line < loc.call.line)
  expect_true(any(grepl("^Re-Filtering individuals", td)))
})

#==========errors==========================
test_that("errors",{
  td <- .internal.data$test_snps[-c(1:2, 5, 8, 9:10)]
  expect_error(filter_snps(td, min_ind = .99, verbose = FALSE), "No loci passed filters.")
  td <- .internal.data$test_snps[,-c(7:8)]
  expect_error(filter_snps(td, min_loci = .99, non_poly = FALSE, verbose = FALSE), "No individuals passed filters.")
})

#==========filter reporting================
test_that("reporting",{
  .make_it_quiet(x <- filter_snps(stickSNPs, 
                                  min_ind = .75, 
                                  min_loc = .75, maf = .1, 
                                  maf_facets = "pop",
                                  hwe = .01, 
                                  hwe_excess_side = "heterozygote",
                                  hwe_facets = "pop",
                                  LD_prune_sigma = 1600,
                                  LD_prune_r = .8))
  .make_it_quiet(res <- filters(x))
  
  expect_equal(res, "bi-allelic;bi-allelic;non-polymorphic;min_ind=0.75,facet=.base;maf=0.1,facet=pop;hwe=0.01, excess side = heterozygote,facet=pop;LD pruning=r=0.8 window_sigma=1600 window_step=1600 LD_method=CLD,facet=.base;min_loci=0.75;non-polymorphic")
})



