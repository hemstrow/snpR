context("fst and ld")

test_that("correct genepop", {
  local_edition(3)
  skip_on_cran()
  tdfst <- calc_pairwise_fst(.internal.data$test_snps, "pop", "genepop")
  tdfst <- get.snpR.stats(tdfst, "pop", "fst")
  expect_snapshot_output(tdfst) # note, run off of genepop, not internally calced. Thus checked, but should not change.
})

test_that("correct wc", {
  tdfst <- calc_pairwise_fst(.internal.data$test_snps, "pop", "wc")
  tdfst <- get.snpR.stats(tdfst, "pop", "fst")
  expect_equal(round(tdfst$pairwise$fst, 4), 
               round(c(0.36724566,
                       0.21886514, 
                       -0.11111111,
                       0.27038627,
                       -0.17977528,
                       0.16666667,
                       -0.12957696,
                       -0.03529412,
                       0.16666667,
                       -0.08172128), 4)) # values from pegas, also double checked by hand. Note that heifstat slightly disagrees on a few of these (where the allele is fixed in one loci)
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
  # skip_on_cran()
  tdld <- calc_pairwise_ld(.internal.data$test_snps, CLD = FALSE)
  tdld <- get.snpR.stats(tdld, stats = "ld")
  prox <- tdld$LD$prox
  expect_snapshot_output(tdld$LD$matrices) # hand checked
  expect_snapshot_output(prox)
})

test_that("correct ME ld",{
  local_edition(3)
  # skip_on_cran()
  set.seed(1212)
  tdldme <- calc_pairwise_ld(.internal.data$test_snps, CLD = FALSE, use.ME = TRUE)
  tdldme <- get.snpR.stats(tdldme, stats = "ld")
  prox <- tdldme$LD$prox
  expect_snapshot_output(tdldme$LD$matrices) # hand checked
  expect_snapshot_output(prox)
})

