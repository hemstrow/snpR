test_that("correct genepop", {
  local_edition(3)
  skip_if_not_installed("genepop");
  tdfst <- calc_pairwise_fst(.internal.data$test_snps, "pop", "genepop")
  tdfst <- get.snpR.stats(tdfst, "pop", "fst")
  
  expect_equal(tdfst$pairwise$fst, 
               c(0.0341, 0, -0.1229, 0, -0.0937, -0.026, 0.1667, 0.1071, -0.0577, -0.087, 0.0165))
  expect_equal(tdfst$weighted.means$weighted_mean_fst, -0.0034)
})

test_that("correct wc", {
  tdfst <- calc_pairwise_fst(.internal.data$test_snps, "pop", "wc")
  tdfst <- get.snpR.stats(tdfst, "pop", "fst")
  expect_equal(round(tdfst$pairwise$fst, 4), 
               round(c(0.034091, 0, -0.122941, 0, -0.09375, -0.026012, 
                       0.166667, 0.107143, -0.057692, -0.086957, 0.016548), 
                     4)) # values from pegas, also double checked by hand. Note that heirfstat slightly disagrees on a few of these (where the allele is fixed in one loci)
  expect_equal(round(tdfst$weighted.means$weighted_mean_fst, 5), -0.00343) # hand calced due to weighting and ratio of averages approach
  
  # var comps
  expect_false(any(paste0("var_comp_", c("a", "b", "c")) %in% colnames(tdfst$pairwise))) # didn't ask for var comps
  tdfst <- calc_pairwise_fst(.internal.data$test_snps, "pop", "wc", keep_components = TRUE)
  tdfst <- get.snpR.stats(tdfst, "pop", "fst")
  expect_true(all(paste0("var_comp_", c("a", "b", "c")) %in% colnames(tdfst$pairwise))) # did ask for var comps
  
  # both weighted and unweighted fsts
  expect_true(all(c("weighted_mean_fst", "mean_fst") %in% colnames(tdfst$weighted.means)))
})



test_that("correct wc, non-bi", {
  keeps <- which(sample.meta(steelMSATs)$pop %in% c("sumhat", "winhat"))
  tdfst <- read_non_biallelic(genotypes(steelMSATs)[,keeps], sample.meta = sample.meta(steelMSATs)[keeps,])
  tdfst <- calc_pairwise_fst(tdfst, "pop", "wc")
  tdfstr <- get.snpR.stats(tdfst, "pop", "fst")$pairwise
  expect_equal(round(tdfstr$fst[which(tdfstr$comparison == "sumhat~winhat")], 4), 
               c(0.0255,0.0326,0.0276,0.0864,0.0614,0.0204,0.0284,0.0165,
                 0.1458,0.0036,0.0542,0.0445,0.0551)) # values from pegas, also double checked by hand.
  
  # check that bootstrapping works
  set.seed(1232)
  tdfst <- calc_pairwise_fst(tdfst, "pop", boot = 10)
  expect_equal(round(get.snpR.stats(tdfst, "pop", "fst")$weighted.means$weighted_mean_fst_p, 3), 
               0.091)
})

test_that("zfst and 1/1-fst",{
  tdfst <- calc_pairwise_fst(.internal.data$test_snps, "pop", "wc", zfst = TRUE, fst_over_one_minus_fst = TRUE)
  tdfst <- get.snpR.stats(tdfst, "pop", "fst")
  expect_true(all(c("zfst", "fst_id") %in% colnames(tdfst$pairwise)))
  expect_true(all(round(tapply(tdfst$pairwise$zfst, tdfst$pairwise$comparison, mean), 3) == 0))
})


test_that("fst errors",{
  x <- .internal.data$test_snps
  sample.meta(x)$test <- "fill"
  expect_error(calc_pairwise_fst(x, c("pop", "test")), "Some facets do not have more than one level")
  
  expect_error(calc_pairwise_fst(x, c("pop", "not_a_facet")), "not found in x metadata")
  
})

test_that("correct fis",{
  # check
  x <- calc_fis(stickSNPs[1:10, pop = c("ASP")], c("pop"))
  fis <- get.snpR.stats(x, facets = "pop", stats =  "fis")
  
  # peg <- snpR::format_snps(x, facets = "pop", output = "adegenet")
  # peg <- pegas::genind2loci(peg)
  # peg_fis_ASP <- pegas::Fst(peg[peg$population == "ASP",])

  
  # checked vs pegas with above code
  expect_equivalent(round(fis$single$fis, 3), round(c(0.384615, 0, -0.052632, -0.111111, 0.142857, 
                                                      -0.285714, 0, -0.25, 0.268293, -0.052632), 3))
  
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
  expect_equal(nrow(fis_b$single), 11)
  expect_equal(unique(unlist(fis_b$weighted.means[1,1:4, drop = TRUE])),
               ".base")
  
  # check that merging works fine
  fis <- calc_pi(.internal.data$test_snps)
  fis <- calc_fis(fis)
  fis <- calc_fis(fis, "pop")
  fisr <- get.snpR.stats(fis, "pop", c("fis", "ho"))$single
  a <- unique(cbind.data.frame(subfacet = fisr$subfacet, lev = .paste.by.facet(fisr, c("chr", "position"))))
  a <- dplyr::mutate_all(a, as.character)
  a <- dplyr::arrange(a, subfacet, lev)
  b <- unique(expand.grid(subfacet = unique(sample.meta(fis)$pop),
                          lev = .paste.by.facet(unique(snp.meta(fis)[,1:2]), c("chr", "position"))))
  b <- dplyr::mutate_all(b, as.character)
  b <- dplyr::arrange(b, subfacet, lev)
  expect_equivalent(a,b)
  expect_true(all(c("ho", "fis") %in% colnames(fisr)))
  
  # check that we have components if requested and not if not
  expect_false(any(c("var_comp_b", "var_comp_c", "nk") %in% colnames(fisr)))
  fis <- calc_fis(fis, "pop", keep_components = TRUE)
  fisr <- get.snpR.stats(fis, "pop", "fis")$single
  expect_true(all(c("var_comp_b", "var_comp_c", "nk") %in% colnames(fisr)))
  
  # check that we are producing the correct averages
  fis <- calc_fis(.internal.data$test_snps, c("pop", "pop.chr"), keep_components = TRUE)
  fisr <- get.snpR.stats(fis, c("pop", "pop.chr"), "fis")$weighted.means
  expect_true(nrow(fisr) == 18)
  
})

test_that("fst bootstrapping",{
  bs1 <- calc_pairwise_fst(.internal.data$test_snps, "pop", boot = 10)
  bs1_res <- get.snpR.stats(bs1, "pop", "fst")
  
  expect_true(is.numeric(unlist(bs1_res$fst.matrix$pop$p[1,2])))
  expect_true(is.numeric(bs1_res$weighted.means$weighted_mean_fst_p))
  
  skip_on_cran();
  bs1_par <- calc_pairwise_fst(.internal.data$test_snps, "pop", boot = 10)
  bs1_res <- get.snpR.stats(bs1_par, "pop", "fst")
  
  expect_true(is.numeric(unlist(bs1_res$fst.matrix$pop$p[1,2])))
  expect_true(is.numeric(bs1_res$weighted.means$weighted_mean_fst_p))
  
  bs2 <- calc_pairwise_fst(.internal.data$test_snps, "pop", method = "genepop", boot = 10)
  bs2_res <- get.snpR.stats(bs2, "pop", "fst")
  
  expect_true(is.numeric(unlist(bs2_res$fst.matrix$pop$p[1,2])))
  expect_true(is.numeric(bs2_res$weighted.means$weighted_mean_fst_p))
  
  bs2_par <- calc_pairwise_fst(.internal.data$test_snps, "pop", method = "genepop", boot = 10)
  bs2_res <- get.snpR.stats(bs2_par, "pop", "fst")
  
  expect_true(is.numeric(unlist(bs2_res$fst.matrix$pop$p[1,2])))
  expect_true(is.numeric(bs2_res$weighted.means$weighted_mean_fst_p))
  
  # complex facets
  bs1 <- calc_pairwise_fst(.internal.data$test_snps, c("pop", "fam", "pop.chr"), boot = 10)
  bs1_res <- get.snpR.stats(bs1, c("pop", "fam", "pop.chr"), "fst")
  expect_true("chr" %in% bs1_res$weighted.means$snp.facet)
  expect_true("fam" %in% bs1_res$weighted.means$facet)
  expect_true(any(bs1_res$weighted.means$snp.facet == "chr" & bs1_res$weighted.means$facet == "pop"))
  expect_false(any(bs1_res$weighted.means$snp.facet == "chr" & bs1_res$weighted.means$facet == "fam"))
  
  # par
  bs1 <- calc_pairwise_fst(.internal.data$test_snps, c("pop", "fam", "pop.chr"), boot = 10, boot_par = 2)
  bs1_res <- get.snpR.stats(bs1, c("pop", "fam", "pop.chr"), "fst")
  expect_true("chr" %in% bs1_res$weighted.means$snp.facet)
  expect_true("fam" %in% bs1_res$weighted.means$facet)
  expect_true(any(bs1_res$weighted.means$snp.facet == "chr" & bs1_res$weighted.means$facet == "pop"))
  expect_false(any(bs1_res$weighted.means$snp.facet == "chr" & bs1_res$weighted.means$facet == "fam"))
})

test_that("correct cld ld",{
  tdld <- calc_pairwise_ld(.internal.data$test_snps, "pop")
  expect_error(tdld <- get.snpR.stats(tdld, "chr", "ld"), "No LD values calculated for these facets")
  
  tdld <- calc_pairwise_ld(tdld, "chr")
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
  skip_on_cran();
  tdld <- calc_pairwise_ld(.internal.data$test_snps, CLD = FALSE)
  tdld <- get.snpR.stats(tdld, stats = "ld")
  prox <- tdld$LD$prox
  expect_snapshot_value(tdld$LD$matrices, style = "serialize") # hand checked
  expect_snapshot_value(prox, style = "serialize")
  
  n <- nsnps(.internal.data$test_snps)
  expect_equal(nrow(prox), n*(n-1)/2)
})

test_that("correct ME ld",{
  local_edition(3)
  skip_on_cran();
  set.seed(1212)
  tdldme <- calc_pairwise_ld(.internal.data$test_snps, CLD = FALSE, use.ME = TRUE)
  tdldme <- get.snpR.stats(tdldme, stats = "ld")
  prox <- tdldme$LD$prox
  expect_snapshot_value(tdldme$LD$matrices, style = "serialize") # hand checked
  expect_snapshot_value(prox, style = "serialize")
  
  n <- nsnps(.internal.data$test_snps)
  expect_equal(nrow(prox), n*(n-1)/2)
})

test_that("correct window ld",{
  skip_on_cran();
  set.seed(1212)
  tdld <- calc_pairwise_ld(stickSNPs, c("pop", "pop.chr", ".base", "chr"), window_gaussian = FALSE, window_sigma = 1600, window_triple_sigma = FALSE)
  tdld <- get.snpR.stats(tdld, c("pop", "pop.chr", ".base", "chr"), stats = "ld")$single.window
  
  expect_true(any(tdld$facet == ".base" & tdld$snp.facet == ".base"))
  expect_true(any(tdld$facet == ".base" & tdld$snp.facet == "chr"))
  expect_true(any(tdld$facet == "pop" & tdld$snp.facet == ".base"))
  expect_true(any(tdld$facet == "pop" & tdld$snp.facet == "chr"))
  expect_true(all(c("sigma", "step", "gaussian", "n_snps", "triple_sigma", "CLD") %in% colnames(tdld)))
  
  # correct facet-level comps
  prox <- calc_pairwise_ld(stickSNPs, c("pop.chr","chr"), window_gaussian = FALSE, window_sigma = 1600, window_triple_sigma = FALSE, .prox_only = TRUE)
  expect_true(all(prox$s1_chr == prox$s2_chr))
  prox <- calc_pairwise_ld(stickSNPs, c("pop.chr","chr"), window_gaussian = FALSE, window_sigma = 1600, window_triple_sigma = FALSE, .prox_only = TRUE, CLD = FALSE)
  expect_true(all(prox$s1_chr == prox$s2_chr))
  
  tdld <- calc_pairwise_ld(stickSNPs, c("pop", "pop.chr", ".base", "chr"), window_gaussian = FALSE, window_sigma = 1600, window_triple_sigma = FALSE, CLD = TRUE)
  tdld <- get.snpR.stats(tdld, c("pop", "pop.chr", ".base", "chr"), stats = "ld")$single.window
  expect_true(all(c("sigma", "step", "gaussian", "n_snps", "triple_sigma", "CLD", "rsq", "Dprime", "pval") %in% colnames(tdld)))
  
  # par the same as serial
  tdldpar <- calc_pairwise_ld(stickSNPs, c("pop", "pop.chr", ".base", "chr"), window_gaussian = FALSE, window_sigma = 1600, 
                              window_triple_sigma = FALSE, CLD = TRUE, par = 2)
  tdldpar <- get.snpR.stats(tdldpar, c("pop", "pop.chr", ".base", "chr"), stats = "ld")$single.window
  expect_equal(dplyr::arrange(tdld), dplyr::arrange(tdldpar), ignore_attr = TRUE)
})

test_that("fis bootstrapping",{
  bs1 <- calc_fis(.internal.data$test_snps, "pop", boot = 10)
  bs1_res <- get.snpR.stats(bs1, "pop", "fis")
  
  # basic
  expect_true("weighted_mean_fis_p" %in% colnames(bs1_res$weighted.means))
  expect_true("weighted_mean_fis_uCI" %in% colnames(bs1_res$weighted.means))
  expect_true("weighted_mean_fis_lCI" %in% colnames(bs1_res$weighted.means))
  skip_if_not("weighted_mean_fis_p" %in% colnames(bs1_res$weighted.means))
  
  
  
  expect_true(is.numeric(unlist(bs1_res$weighted.means$weighted_mean_fis_p)))

  # complex facets
  bs1_par <- calc_fis(.internal.data$test_snps, c("pop", "fam", "pop.chr", "chr"), boot = 10)
  bs1_res <- get.snpR.stats(bs1_par, c("pop", "fam", "pop.chr", "chr"), "fis")
  expect_true("weighted_mean_fis_p" %in% colnames(bs1_res$weighted.means))
  expect_true(all(c("fam", "pop", ".base") %in% bs1_res$weighted.means$facet))
  expect_true(all(bs1_res$weighted.means[bs1_res$weighted.means$facet == ".base",]$subfacet == ".base"))
  expect_true(all(bs1_res$weighted.means[bs1_res$weighted.means$facet == "fam",]$snp.subfacet == ".base"))
  expect_true(all(c(".base", "chr") %in% bs1_res$weighted.means[bs1_res$weighted.means$facet == "pop",]$snp.facet))
  
  # base facet
  bs1 <- calc_fis(.internal.data$test_snps, boot = 10)
  bs1_res <- get.snpR.stats(bs1, stats =  "fis")
  expect_true("weighted_mean_fis_p" %in% colnames(bs1_res$weighted.means))
  expect_true(all(unlist(bs1_res$weighted.mean[,c("facet", "subfacet", "snp.facet", "snp.subfacet")]) == ".base"))
  
  # base facet mixed in
  bs1 <- calc_fis(.internal.data$test_snps, c("pop", ".base"), boot = 10)
  bs1_res <- get.snpR.stats(bs1, c(".base", "pop"), stats =  "fis")
  expect_true("weighted_mean_fis_p" %in% colnames(bs1_res$weighted.means))
  expect_true(all(bs1_res$weighted.means[bs1_res$weighted.means$facet == ".base",]$subfacet == ".base"))
  expect_true(all(bs1_res$weighted.means[bs1_res$weighted.means$facet == "pop",]$snp.subfacet == ".base"))
  
  # parallel
  skip_on_cran();
  bs1_par <- calc_fis(.internal.data$test_snps, "pop", boot = 10, boot_par = 2)
  bs1_res <- get.snpR.stats(bs1_par, "pop", "fis")
  expect_true("weighted_mean_fis_p" %in% colnames(bs1_res$weighted.means))
})

test_that("global fst",{
  
  # should yield the same as global = FALSE for two comps
  tdfst <- calc_global_fst(.internal.data$test_snps, "pop")
  tdfst <- get.snpR.stats(tdfst, "pop", "fst")
  expect_equal(round(tdfst$pairwise$fst, 4),
               round(c(0.034091, 0, -0.122941, 0, -0.09375, -0.026012,
                       0.166667, 0.107143, -0.057692, -0.086957, 0.016548),
                     4)) # values from pegas, also double checked by hand. Note that heirfstat slightly disagrees on a few of these (where the allele is fixed in one loci)
  expect_equal(round(tdfst$weighted.means$weighted_mean_fst, 5), -0.00343) # hand calced due to weighting and ratio of averages approach
  
  # actual multipop
  tdfst <- calc_global_fst(stickSNPs, "pop")
  tdfst <- get.snpR.stats(tdfst, "pop", "fst")
  
  # peg <- snpR::format_snps(stickSNPs, facets = "pop", output = "adegenet")
  # peg <- pegas::genind2loci(peg)
  # peg_fst <- pegas::Fst(peg)
  # round(peg_fst[,2], 3)
  
  expect_equivalent(round(tdfst$pairwise$fst, 3)[1:10],
                    c(0.078, -0.002, -0.025, 0.143, 0.011, 0.136, 0.059, 0.08, 0.04, 0.025)) # from pegas, see above
  
  expect_equal(round(tdfst$weighted.means$weighted_mean_fst, 3), 0.094)
  
  # both weighted and unweighted fsts
  expect_true(all(c("weighted_mean_fst", "mean_fst") %in% colnames(tdfst$weighted.means)))
})
