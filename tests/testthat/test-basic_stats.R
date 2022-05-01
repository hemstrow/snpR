context("Maf, pi, ho, private, HWE (with a facet)")

tdm <- calc_basic_snp_stats(.internal.data$test_snps, c("pop", ".base"))

test_that("maf",{
  # run function
  expect_true(.check_calced_stats(tdm, "pop", "maf")$pop["maf"])
  maf <- head(get.snpR.stats(tdm, "pop", "maf")$single, 2)
  
  # test
  expect_equal(round(maf$maf, 5), c(.5, 0) # hand calced
               ) # correct maf?
  expect_equal(maf$major, c("C", "C"))
  expect_equal(maf$minor, c("G", "G")) # note, the 2nd to last A will be N if it fails to account for the overall minor, since that facet has none of the minor alleles.
  expect_equal(maf$maj.count, c(5, 6))
  expect_equal(maf$min.count, c(5, 0))
})


test_that("pi",{
  # run function
  expect_true(all(unlist(.check_calced_stats(tdm, c(".base", "pop"), "pi"))))
  pi <- get.snpR.stats(tdm, c(".base", "pop"), "pi")$single

  # test
  expect_equal(pi$pi, 1 - (choose(tdm@ac$ni1, 2) + choose(tdm@ac$ni2, 2) )/choose(tdm@ac$n_total, 2)) # equ from hohenlohe
})

test_that("he",{
  # run function
  expect_true(all(unlist(.check_calced_stats(tdm, c(".base", "pop"), "he"))))
  he <- get.snpR.stats(tdm, c(".base", "pop"), "he")$single
  
  # test
  expect_equal(he$he, 2 * (tdm@ac$ni1/tdm@ac$n_total) * (tdm@ac$ni2/tdm@ac$n_total)) # check against 2pq from another source
})

test_that("ho", {
  # run function
  expect_true(all(unlist(.check_calced_stats(tdm, c(".base", "pop"), "ho"))))
  ho <- get.snpR.stats(tdm, c("pop", ".base"), "ho")$single
  
  rs <- rowSums(tdm@geno.tables$gs)
  hets <- rowSums(tdm@geno.tables$gs[,which(substr(colnames(tdm@geno.tables$g), 1, 1) != substr(colnames(tdm@geno.tables$g), 2, 2))]) # hand calced
  expect_equal(ho$ho, hets/rs)
})


test_that("private", {
  tdm <- calc_private(tdm, "pop")
  expect_true(all(unlist(.check_calced_stats(tdm, c("pop"), "pa"))))
  pa <- get.snpR.stats(tdm, "pop", "private")$single
  
  expect_equal(pa$pa, c(1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0)) # hand calced
})

test_that("hwe", {
  # bigger sample size needed
  suppressWarnings(tdhwe <- subset_snpR_data(stickSNPs, 1:10, 1:10))
  tdhwe <- calc_hwe(tdhwe, ".base")
  expect_true(all(unlist(.check_calced_stats(tdhwe, ".base", "hwe"))))
  hwe.exact <- get.snpR.stats(tdhwe, ".base", "hwe")$single
  
  # exact
  expect_equal(round(hwe.exact$pHWE,3), c(1, .341, 1, 1, 1, 1, 1, 1, .167, 1)) # hand calced
  
  # chisq
  tdhwe <- calc_hwe(tdhwe, ".base", method = "chisq")
  hwe.chsq <- get.snpR.stats(tdhwe, ".base", "hwe")$single
  expect_equal(round(hwe.chsq$pHWE, 3), c(1, 0.284, 1, 1, .839, .659, .577, .860, .134, 1)) # from pegas
})

test_that("hs", {
  hs <- calc_hs(.internal.data$test_snps)
  hs <- get.snpR.stats(hs, stats = "hs")
  expect_equal(round(hs$sample$hs, 3), 
               round(c(1.2761109, 0.7468880, 1.9240313, 0.8304498, 1.2761109, 
                       1.2761109, 0.3190277, 1.1464968, 0.6380555, 0.6933553), 3)) # hand calced
  expect_equal(nrow(hs$weighted.means), 1)
  
  hs <- calc_hs(.internal.data$test_snps, "pop")
  hs <- get.snpR.stats(hs, "pop", "hs")
  expect_equal(nrow(hs$weighted.means), 2)
  
  hs <- calc_hs(.internal.data$test_snps, "pop.chr", complex_averages = TRUE)
  hs <- get.snpR.stats(hs, "pop.chr", "hs")
  expect_equal(sort(.paste.by.facet(hs$weighted.means, c("subfacet", "snp.subfacet"), "_")),
               sort(.paste.by.facet(expand.grid(unique(sample.meta(.internal.data$test_snps)$pop),
                                                unique(snp.meta(.internal.data$test_snps)$chr)),
                                    c("Var1", "Var2"), "_"))) # every level accounted for?
})

test_that("het_hom", {
  het_hom_ratio <- calc_het_hom_ratio(.internal.data$test_snps)
  het_hom_ratio <- get.snpR.stats(het_hom_ratio, stats = "het_hom_ratio")
  expect_equal(round(het_hom_ratio$sample$`Het/Hom`, 3),
               c(0.667 ,0.333 ,1.333 ,0.400 ,0.667 ,0.667 ,0.111 ,0.600 ,0.250 ,0.286))# hand calced
  
  expect_equal(nrow(het_hom_ratio$weighted.means), 1)
  
  het_hom_ratio <- calc_het_hom_ratio(.internal.data$test_snps, "pop")
  het_hom_ratio <- get.snpR.stats(het_hom_ratio, "pop", "het_hom_ratio")
  expect_equal(nrow(het_hom_ratio$weighted.means), 2)
  
  het_hom_ratio <- calc_het_hom_ratio(.internal.data$test_snps, "pop.chr", complex_averages = TRUE)
  het_hom_ratio <- get.snpR.stats(het_hom_ratio, "pop.chr", "het_hom_ratio")
  expect_equal(sort(.paste.by.facet(het_hom_ratio$weighted.means, c("subfacet", "snp.subfacet"), "_")),
               sort(.paste.by.facet(expand.grid(unique(sample.meta(.internal.data$test_snps)$pop),
                                                unique(snp.meta(.internal.data$test_snps)$chr)),
                                    c("Var1", "Var2"), "_"))) # every level accounted for?
})

  