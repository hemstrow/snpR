context("Maf, pi, ho, private, HWE (with a facet)")

tdm <- calc_basic_snp_stats(.internal.data$test_snps, c("pop", ".base"))

test_that("maf",{
  # run function
  expect_true(check_calced_stats(tdm, "pop", "maf")$pop["maf"])
  maf <- head(get.snpR.stats(tdm, "pop", "maf")$single, 2)
  
  # test
  expect_equal(round(maf$maf, 5), c(.4, .3) # hand calced
               ) # correct maf?
  expect_equal(maf$major, c("C", "C"))
  expect_equal(maf$minor, c("G", "G")) # note, the 2nd to last A will be N if it fails to account for the overall minor, since that facet has none of the minor alleles.
  expect_equal(maf$maj.count, c(6, 7))
  expect_equal(maf$min.count, c(4, 3))
})


test_that("pi",{
  # run function
  expect_true(all(unlist(check_calced_stats(tdm, c(".base", "pop"), "pi"))))
  pi <- get.snpR.stats(tdm, c(".base", "pop"), "pi")$single

  # test
  expect_equal(pi$pi, 1 - (choose(tdm@ac$ni1, 2) + choose(tdm@ac$ni2, 2) )/choose(tdm@ac$n_total, 2)) # equ from hohenlohe
})

test_that("ho", {
  # run function
  expect_true(all(unlist(check_calced_stats(tdm, c(".base", "pop"), "ho"))))
  ho <- get.snpR.stats(tdm, c("pop", ".base"), "ho")$single
  
  rs <- rowSums(tdm@geno.tables$gs)
  hets <- rowSums(tdm@geno.tables$gs[,which(substr(colnames(tdm@geno.tables$g), 1, 1) != substr(colnames(tdm@geno.tables$g), 2, 2))]) # hand calced
  expect_equal(ho$ho, hets/rs)
})


test_that("private", {
  tdm <- calc_private(tdm, "pop")
  expect_true(all(unlist(check_calced_stats(tdm, c("pop"), "pa"))))
  pa <- get.snpR.stats(tdm, "pop", "private")$single
  
  expect_equal(pa$pa, c(0, 0, 0, 1, rep(0,9), 1, 0, 1, 0, 0)) # hand calced
})

test_that("hwe", {
  # bigger sample size needed
  suppressWarnings(tdhwe <- subset_snpR_data(stickSNPs, 1:10, 1:10))
  tdhwe <- calc_hwe(tdhwe, ".base")
  expect_true(all(unlist(check_calced_stats(tdhwe, ".base", "hwe"))))
  hwe.exact <- get.snpR.stats(tdhwe, ".base", "hwe")$single
  
  # exact
  expect_equal(round(hwe.exact$pHWE,3), c(1, .341, 1, 1, 1, 1, 1, 1, .167, 1)) # hand calced
  
  # chisq
  tdhwe <- calc_hwe(tdhwe, ".base", method = "chisq")
  hwe.chsq <- get.snpR.stats(tdhwe, ".base", "hwe")$single
  expect_equal(round(hwe.chsq$pHWE, 3), c(1, 0.284, 1, 1, .839, .659, .577, .860, .134, 1)) # from pegas
})





  