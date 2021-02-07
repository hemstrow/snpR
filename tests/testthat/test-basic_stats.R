context("Maf, pi, ho, private, HWE (with a facet)")

tdat <- stickRAW[1:4, 1:23]
tdat <- import.snpR.data(tdat[,-c(1:3)], tdat[,1:3], data.frame(pop = rep(c("A", "B"), each = 10)))
tdat <- add.facets.snpR.data(tdat, "pop")

test_that("maf",{
  # run function
  tdm <- calc_maf(tdat, "pop")
  expect_true(check_calced_stats(tdm, "pop", "maf")$pop["maf"])
  maf <- tdm@stats
  
  # test
  expect_equal(round(maf$maf, 5),
               round(1  - (matrixStats::rowMaxs(tdm@geno.tables$as)/rowSums(tdm@geno.tables$as)), 5)) # correct maf?
  expect_equal(tdm@stats$major, c("C", "C", "C", "T", "T", "T", "G", "G", "G", "T", "T", "T"))
  expect_equal(tdm@stats$minor, c("T", "T", "T", "C", "C", "C", "N", "N", "N", "A", "A", "A")) # note, the 2nd to last A will be N if it fails to account for the overall minor, since that facet has none of the minor alleles.
  expect_equal(tdm@stats$maj.count, c(31, 12, 19, 17, 9, 8, 26, 16, 10, 35, 18, 17))
  expect_equal(tdm@stats$min.count, c(3, 2, 1, 15, 7, 8, 0, 0, 0, 1, 0, 1))
  expect_equal(tdm@stats$maf, tdm@stats$min.count/(tdm@stats$maj.count + tdm@stats$min.count))
  rm(tdm)
})


test_that("pi",{
  # run function
  tdp <- calc_pi(tdat)
  tdp <- calc_pi(tdp, "pop")
  expect_true(all(unlist(check_calced_stats(tdp, c(".base", "pop"), "pi"))))
  pi <- tdp@stats
  

  # test
  expect_equal(pi$pi, 1 - (choose(tdp@ac$ni1, 2) + choose(tdp@ac$ni2, 2) )/choose(tdp@ac$n_total, 2))
})

test_that("ho", {
  # run function
  tdh <- calc_ho(tdat)
  tdh <- calc_ho(tdh, "pop")
  expect_true(all(unlist(check_calced_stats(tdh, c(".base", "pop"), "ho"))))
  ho <- tdh@stats
  
  rs <- rowSums(tdh@geno.tables$gs)
  hets <- rowSums(tdh@geno.tables$gs[,which(substr(colnames(tdh@geno.tables$g), 1, 1) != substr(colnames(tdh@geno.tables$g), 2, 2))])
  expect_equal(ho$ho, hets/rs)
})


test_that("private", {
  tdpa <- calc_private(tdat)
  tdpa <- calc_private(tdpa, "pop")
  expect_true(all(unlist(check_calced_stats(tdpa, c(".base", "pop"), "pa"))))
  pa <- tdpa@stats
  
  expect_equal(pa$pa, c(rep(0, 11), 1)) # the last pop/locus should be private
})

test_that("hwe", {
  # bigger sample size needed
  tdhwe <- stickRAW[1:4, 1:203]
  tdhwe <- import.snpR.data(tdhwe[,-c(1:3)], tdhwe[,1:3], data.frame(pop = rep(c("A", "B"), each = 100)))
  tdhwe <- calc_hwe(tdhwe)
  tdhwe <- calc_hwe(tdhwe, "pop")
  expect_true(all(unlist(check_calced_stats(tdhwe, c(".base", "pop"), "hwe"))))
  hwe.exact <- tdhwe@stats
  
  # exact
  expect_equal(hwe.exact$pHWE, c(1, 1, 1, 0.4916573, 0.3667300, rep(1, 7)))
  
  # chisq
  hwe.chsq <- calc_hwe(tdhwe, c(".base", "pop"), method = "chisq")
  hwe.chsq <- hwe.chsq@stats
  expect_equal(round(hwe.chsq$pHWE, 3), 
               round(c(0.8272968, 0.8756337, 0.9395070, 0.7877324, 0.5869697, 
                       0.9926668, 0.9762565, 0.9981080,
                       0.9692610, 0.7445803, 0.9507463, 0.7376636)), 3)
  
  
})

  