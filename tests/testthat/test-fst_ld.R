context("fst and ld")

test_that("correct genepop", {
  local_edition(3)
  tdfst <- calc_pairwise_fst(.internal.data$test_snps, "pop", "genepop")
  tdfst <- get.snpR.stats(tdfst, "pop", "fst")
  expect_snapshot(tdfst) # note, run off of genepop, not internally calced. Thus checked, but should not change.
})

test_that("correct wc", {
  tdfst <- calc_pairwise_fst(.internal.data$test_snps, "pop", "wc")
  tdfst <- get.snpR.stats(tdfst, "pop", "fst")
  expect_equal(round(tdfst$pairwise$fst, 4), 
               round(c(-8.695652e-02, -1.040834e-16,  4.411765e-02, -1.111111e-01,
                       -1.333741e-01, -3.978780e-02,  0.000000e+00,  1.195278e-01,
                       -5.769231e-02 ), 4)) # values from pegas, also double checked by hand. Note that heifstat slightly disagrees on a few of these (where the allele is fixed in one loci)
})


test_that("correct cld ld",{
  tdld <- calc_pairwise_ld(.internal.data$test_snps, "pop")
  tdld <- get.snpR.stats(tdld, "pop", "ld")
  tdld <- tdld$LD$matrices$ASP$.base$CLD
  tdld <- as.numeric(tdld)
  
  suppressWarnings(ASPcor <- cor(t(format_snps(.internal.data$test_snps, "sn", interpolate = FALSE)[,-c(1:2)][,which(sample.meta(.internal.data$test_snps)$pop == "ASP")]), 
                use = "pairwise.complete.obs")^2)
  ASPcor[which(lower.tri(ASPcor))] <- NA
  diag(ASPcor) <- NA
  ASPcor <- as.numeric(ASPcor)
  expect_equal(ASPcor, tdld) # checked by hand vs definition of CLD LD.
})

test_that("correct traditional ld",{
  local_edition(3)
  tdld <- calc_pairwise_ld(.internal.data$test_snps, CLD = FALSE)
  tdld <- get.snpR.stats(tdld, stats = "ld")
  expect_snapshot(tdld) # hand checked
})

test_that("correct ME ld",{
  local_edition(3)
  tdldme <- calc_pairwise_ld(.internal.data$test_snps, CLD = FALSE, use.ME = TRUE)
  tdldme <- get.snpR.stats(tdldme, stats = "ld")
  expect_snapshot(tdldme) # hand checked
})

