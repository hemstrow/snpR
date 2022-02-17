context("fst and ld")

test_that("correct genepop", {
  local_edition(3)
  skip_on_cran()
  tdfst <- calc_pairwise_fst(.internal.data$test_snps, "pop", "genepop")
  tdfst <- get.snpR.stats(tdfst, "pop", "fst")
  expect_snapshot_value(tdfst, style = "serialize") # note, run off of genepop, not internally calced. Thus checked, but should not change.
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
                       -0.08172128), 4)) # values from pegas, also double checked by hand. Note that heirfstat slightly disagrees on a few of these (where the allele is fixed in one loci)
})

test_that("correct fis",{
  # check
  x <- calc_fis(stickSNPs[1:10, pop = c("ASP")], c("pop"))
  fis <- get.snpR.stats(x, facets = "pop", stats =  "fis")
  
  # peg <- snpR::format_snps(x, facets = "pop", output = "adegenet")
  # peg <- pegas::genind2loci(peg)
  # peg_fis_ASP <- pegas::Fst(peg[peg$population == "ASP",])
  # peg_fis_ASP <- pegas::Fst(peg[peg$population == "ASP",])[2:10,]
  
  
  # checked vs pegas with above code
  expect_equivalent(round(fis$single$fis, 3), round(c(-0.05179856,
                                                      -0.02272727,
                                                      -0.03703704,
                                                      0.02222222,
                                                      -0.08641975,
                                                      -0.11904762,
                                                      -0.06976744,
                                                      -0.23888183,
                                                      -0.01492537), 3))
  
  # check that the means are the same if we do additional facets at the same time
  fis_mf <- calc_fis(stickSNPs[1:10, pop = c("ASP", "PAL")], c("pop", "pop.fam"))
  fis_mf <- get.snpR.stats(fis_mf, facets = c("pop", "pop.fam"), stats = "fis")
  expect_equal(fis_mf$weighted.means[which(fis_mf$weighted.means$subfacet == "ASP"),]$weighted_mean_fis,
               fis$weighted.means$weighted_mean_fis)
  
  # check means are working correctly
  fis <- calc_fis(.internal.data$test_snps, c("pop", "pop.chr", ".base"))
  fis_p <- get.snpR.stats(fis, "pop", "fis")$weighted.means
  expect_equal(fis_p$subfacet, c("ASP", "PAL"))
  
  fis_pc <- get.snpR.stats(fis, "pop.chr", "fis")$weighted.means
  expect_equal(unique(fis_pc$subfacet), c("ASP", "PAL"))
  expect_equal(sort(unique(fis_pc$snp.subfacet)), sort(unique(snp.meta(fis)$chr)))
  
  # check that base worked
  fis_b <- get.snpR.stats(fis, stats = "fis")
  expect_equal(nrow(fis_b$single), 10)
  expect_equal(unique(unlist(fis_b$weighted.means[1,1:4, drop = TRUE])),
               ".base")
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
  expect_snapshot_value(tdld$LD$matrices, style = "serialize") # hand checked
  expect_snapshot_value(prox, style = "serialize")
})

test_that("correct ME ld",{
  local_edition(3)
  # skip_on_cran()
  set.seed(1212)
  tdldme <- calc_pairwise_ld(.internal.data$test_snps, CLD = FALSE, use.ME = TRUE)
  tdldme <- get.snpR.stats(tdldme, stats = "ld")
  prox <- tdldme$LD$prox
  expect_snapshot_value(tdldme$LD$matrices, style = "serialize") # hand checked
  expect_snapshot_value(prox, style = "serialize")
})

