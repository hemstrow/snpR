tdm <- calc_basic_snp_stats(.internal.data$test_snps, c("pop", ".base"))

test_that("maf",{
  # run function
  expect_true(.check_calced_stats(tdm, "pop", "maf")$pop["maf"])
  maf <- head(get.snpR.stats(tdm, "pop", "maf")$single, 2)
  
  # test
  expect_equal(round(maf$maf, 5), c(.2, .4) # hand calced
               ) # correct maf?
  expect_equal(maf$major, c("A", "A"))
  expect_equal(maf$minor, c("G", "G")) # note, the 2nd to last A will be N if it fails to account for the overall minor, since that facet has none of the minor alleles.
  expect_equal(maf$maj.count, c(8, 6))
  expect_equal(maf$min.count, c(2, 4))
  
  # test with facets that MAF/MIN is correctly ID'd when flipped
  check <- genotypes(tdm)
  check <- rep(c("AA", "TT"), length.out = ncol(tdm))
  check[2] <- "AT"
  check <- rbind(genotypes(tdm), check)
  check <- import.snpR.data(check, rbind(snp.meta(tdm), data.frame(chr = "test", position = 10, .snp.id = 13)),
                            sample.meta(tdm))
  check <- calc_maf(check, "pop")
  check <- get.snpR.stats(check, "pop", "maf")
  check <- check$single[which(check$single$.snp.id == 13),]
  expect_true(all(check$major == "A"))
  expect_true(all(check$minor == "T"))
  expect_true(all(check$maj.count == c(10, 1)))
  expect_true(all(check$min.count == c(0, 9)))
  
  # test with facets that MAF/MIN is correctly ID'd when flipped, but with fully fixed
  check <- genotypes(tdm)
  check <- rep(c("AA", "TT"), length.out = ncol(tdm))
  check <- rbind(genotypes(tdm), check)
  check <- import.snpR.data(check, rbind(snp.meta(tdm), data.frame(chr = "test", position = 10, .snp.id = 13)),
                            sample.meta(tdm))
  check <- calc_maf(check, "pop")
  check <- get.snpR.stats(check, "pop", "maf")
  check <- check$single[which(check$single$.snp.id == 13),]
  expect_true(all(check$major == "T"))
  expect_true(all(check$minor == "A"))
  expect_true(all(check$maj.count == c(0, 10)))
  expect_true(all(check$min.count == c(10, 0)))
  
  # cases where one locus is completely unsequenced at one pop
  test <- genotypes(.internal.data$test_snps)
  test[1,sample.meta(.internal.data$test_snps)$pop == "ASP"] <- "NN"
  testd <- import.snpR.data(test, 
                            snp.meta(.internal.data$test_snps),
                            sample.meta(.internal.data$test_snps))
  testd <- calc_maf(testd, "pop")
  res <- get.snpR.stats(testd, "pop", "maf")$single
  expect_true(is.na(res[1,]$maf))
  expect_true(res[1,]$major == "A")
  expect_true(res[1,]$minor == "G")
  expect_true(res[1,]$maj.count == 0)
  expect_true(res[1,]$min.count == 0)
  
  # cases where one locus is completely unsequenced at all pops
  test <- genotypes(.internal.data$test_snps)
  test[1,] <- "NN"
  testd <- import.snpR.data(test, 
                            snp.meta(.internal.data$test_snps),
                            sample.meta(.internal.data$test_snps))
  testd <- calc_maf(testd, "pop")
  res <- get.snpR.stats(testd, "pop", "maf")$single
  expect_true(all(is.na(res[1:2,]$maf)))
  expect_true(all(res[1:2,]$major == "N"))
  expect_true(all(res[1:2,]$minor == "N"))
  expect_true(all(is.na(res[1,]$maj.count)))
  expect_true(all(is.na(res[1,]$min.count)))
  
})


test_that("pi",{
  # run function
  expect_true(all(unlist(.check_calced_stats(tdm, c(".base", "pop"), "pi"))))
  pi <- get.snpR.stats(tdm, c(".base", "pop"), "pi")$single

  # test
  expect_equal(pi$pi, 1 - rowSums(apply(as.matrix(tdm@geno.tables$as), MARGIN = 2, function(x) choose(x, 2)))/choose(Matrix::rowSums(tdm@geno.tables$as), 2)) # equ from hohenlohe
})

test_that("he",{
  # run function
  expect_true(all(unlist(.check_calced_stats(tdm, c(".base", "pop"), "he"))))
  he <- get.snpR.stats(tdm, c(".base", "pop"), "he")$single
  
  # test
  as <- tdm@geno.tables$as
  as <- as/Matrix::rowSums(as)
  as <- as^2
  che <- 1 - Matrix::rowSums(as) # 1 minus hom freqs
  expect_equal(he$he, che) # check against 2pq from another source
  
  # non-bi-allelic
  check <- stickSNPs
  check <- calc_he(check, "pop")
  check <- get.snpR.stats(check, "pop", "he")
  
  check2 <- stickSNPs
  check2@bi_allelic <- FALSE
  check2 <- calc_he(check2, "pop")
  check2 <- get.snpR.stats(check2, "pop", "he")
  expect_equivalent(check, check2)
})

test_that("ho", {
  # run function
  expect_true(all(unlist(.check_calced_stats(tdm, c(".base", "pop"), "ho"))))
  ho <- get.snpR.stats(tdm, c("pop", ".base"), "ho")$single
  
  rs <- Matrix::rowSums(tdm@geno.tables$gs)
  hets <- Matrix::rowSums(tdm@geno.tables$gs[,which(substr(colnames(tdm@geno.tables$g), 1, 1) != substr(colnames(tdm@geno.tables$g), 2, 2))]) # hand calced
  expect_equal(ho$ho, hets/rs)
})


test_that("private", {
  x <- calc_private(stickSNPs[pop = c("ASP", "OPL")], "pop", rarefaction = FALSE)
  expect_true(all(unlist(.check_calced_stats(x, c("pop"), "pa"))))
  pa <- get.snpR.stats(x, "pop", "private")
  
  expect_equal(which(pa$single$pa_uncorrected == 1), c(37, 98, 108, 122, 138, 156, 166, 180, 188, 193)) # hand calced
  expect_equal(pa$weighted.means$total_pa_uncorrected, c(2, 8)) # hand calced
  
  
  
  x <- calc_private(x, "pop", rarefaction = TRUE, g = -1)
  pa2 <- get.snpR.stats(x, "pop", "private")
  expect_true(cor(pa2$single$pa_uncorrected, pa2$single$pa_corrected) > .5)
  expect_equal(round(pa2$weighted.means$total_pa_corrected, 3), c(1.996, 5.998))
})

test_that("hwe", {
  # bigger sample size needed
  suppressWarnings(tdhwe <- subset_snpR_data(stickSNPs, 1:10, 1:10))
  tdhwe <- calc_hwe(tdhwe, ".base")
  expect_true(all(unlist(.check_calced_stats(tdhwe, ".base", "hwe"))))
  hwe.exact <- get.snpR.stats(tdhwe, ".base", "hwe")$single
  
  # exact
  expect_equal(round(hwe.exact$pHWE,3), c(0.199, 1, 1, 1, 1, 1, 1, 1, 1, 1)) # hand calced
  
  # chisq
  tdhwe <- calc_hwe(tdhwe, ".base", method = "chisq")
  hwe.chsq <- get.snpR.stats(tdhwe, ".base", "hwe")$single
  expect_equal(round(hwe.chsq$pHWE, 3), c(0.157, 0.860, 0.725, 0.725, 0.860, 0.292, 0.868, 0.292, 0.577, 0.725)) # from pegas
})

test_that("hs", {
  hs <- calc_hs(.internal.data$test_snps)
  hs <- get.snpR.stats(hs, stats = "hs")
  expect_equal(round(hs$sample$hs, 3), 
               round(c(0.824, 2.156, 1.099, 0.947, 0.609, 0.748, 
                       0.966, 0.966, 1.028, 0.824), 3)) # hand calced
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
               c(0.429, 2.000, 0.667, 0.600, 0.250, 0.333, 0.571, 0.571, 0.667, 0.429))# hand calced
  
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

test_that("prop_poly",{
  tf <- c(".base", "chr", "chr.pop", "chr.fam", "fam", "pop", "fam.pop.chr")
  poly <- calc_prop_poly(.internal.data$test_snps, tf)
  polyc <- get.snpR.stats(poly, tf, "prop_poly")$weighted.means
  
  # quick, dirty percent poly function using genotypes to check, not efficient for large data!
  conf_one <- function(dat){
    id <- genotypes(dat)
    id <- t(id)
    id[id == "NN"] <- NA
    id1 <- substr(id, 1, 1)
    id2 <- substr(id, 2, 2)
    id <- rbind(id1, id2)
    id <- apply(id, MARGIN = 2, function(x) length(unique(na.omit(x))))
    return(sum(id != 1)/length(id))
  }
  
  expect_equal(polyc$prop_poly[1:9], c(conf_one(.internal.data$test_snps),
                                 conf_one(.internal.data$test_snps[fam = "A"]),
                                 conf_one(.internal.data$test_snps[fam = "B"]),
                                 conf_one(.internal.data$test_snps[pop = "ASP"]),
                                 conf_one(.internal.data$test_snps[pop = "PAL"]),
                                 conf_one(.internal.data$test_snps[chr = "groupVI"]),
                                 conf_one(.internal.data$test_snps[chr = "groupVI", fam = "A"]),
                                 conf_one(.internal.data$test_snps[chr = "groupVI", fam = "B"]),
                                 conf_one(.internal.data$test_snps[chr = "groupVI", fam.pop = "A.ASP"])))
  # note: hand checked for accuracy
  
  # merging
  poly <- calc_ho(poly, tf)
  expect_equal(get.snpR.stats(poly, tf, c("ho", "prop_poly"))$weighted.means$weighted_mean_ho,
               get.snpR.stats(calc_ho(.internal.data$test_snps, tf), tf, "ho")$weighted.means$weighted_mean_ho)
  
  # double snp level
  r1 <- get.snpR.stats(calc_prop_poly(.internal.data$test_snps, c("chr.position.fam")), "chr.position.fam", "prop_poly")
  r2 <- get.snpR.stats(calc_prop_poly(.internal.data$test_snps, c("position.chr.fam")), "position.chr.fam", "prop_poly")
  expect_equal(r1, r2)
  
  ## everything accounted for?
  check <- .add.facets.snpR.data(.internal.data$test_snps, "fam")
  expect_equal(sort(.paste.by.facet(r1$weighted.means, c("facet", "subfacet", "snp.subfacet"))),
               sort(.paste.by.facet(as.data.frame(.get.task.list(check, "chr.position.fam")), 
                                    c("t.sample.facet", "all.opts.1", "all.opts.2"))))
})

# tajima's D
test_that("tajimas_d",{
  tf <- c(".base", "chr", "chr.pop", "chr.fam", "fam", "pop", "fam.pop.chr")
  expect_warning(tsd <- calc_tajimas_d(.internal.data$test_snps, tf, step = 200, triple_sigma = FALSE, sigma = 400), "Consider adding a snp level facet")
  tsdc <- get.snpR.stats(tsd, tf, "tajimas_d")
  
  # check that all levels are there in both windows and weighted
  ## windows
  levs <- unique(tsdc$single.window[,1:4])
  levs_pasted <- .paste.by.facet(levs, colnames(levs), sep = ".")
  sf <- summarize_facets(.internal.data$test_snps, tf)
  expect_true(".base..base..base..base" %in% levs_pasted)
  expect_true(all(paste0(".base..base.chr.", sf$chr) %in% levs_pasted))
  expect_true(all(paste0("pop.", lapply(strsplit(sf$chr.pop, "\\."), function(x){paste0(x[2], ".chr.", x[1])})) %in% levs_pasted))
  expect_true(all(paste0("fam.", lapply(strsplit(sf$chr.fam, "\\."), function(x){paste0(x[2], ".chr.", x[1])})) %in% levs_pasted))
  expect_true(all(paste0("fam.", names(sf$fam), "..base..base") %in% levs_pasted))
  expect_true(all(paste0("pop.", names(sf$pop), "..base..base") %in% levs_pasted))
  expect_true(all(paste0("fam.pop.", lapply(strsplit(sf$chr.fam.pop, "\\."), function(x){paste0(x[2], ".", x[3], ".chr.", x[1])})) %in% levs_pasted))
  
  ## weighted
  levs <- unique(tsdc$weighted.means[,1:4])
  levs_pasted <- .paste.by.facet(levs, colnames(levs), sep = ".")
  sf <- summarize_facets(.internal.data$test_snps, tf)
  expect_true(".base..base..base..base" %in% levs_pasted)
  expect_true(all(paste0(".base..base.chr.", c(sf$chr, ".OVERALL_MEAN")) %in% levs_pasted))
  expect_true(all(paste0("pop.", lapply(strsplit(sf$chr.pop, "\\."), function(x){paste0(x[2], ".chr.", x[1])})) %in% levs_pasted))
  expect_true(all(paste0("pop.", names(sf$pop), ".chr..OVERALL_MEAN") %in% levs_pasted))
  expect_true(all(paste0("fam.", lapply(strsplit(sf$chr.fam, "\\."), function(x){paste0(x[2], ".chr.", x[1])})) %in% levs_pasted))
  expect_true(all(paste0("fam.", names(sf$fam), ".chr..OVERALL_MEAN") %in% levs_pasted))
  expect_true(all(paste0("fam.", names(sf$fam), "..base..base") %in% levs_pasted))
  expect_true(all(paste0("pop.", names(sf$pop), "..base..base") %in% levs_pasted))
  expect_true(all(paste0("fam.pop.", lapply(strsplit(sf$chr.fam.pop, "\\."), function(x){paste0(x[2], ".", x[3], ".chr.", x[1])})) %in% levs_pasted))
  sfc <- summarize_facets(.internal.data$test_snps, "pop.fam")
  expect_true(all(paste0("fam.pop.", names(sfc$fam.pop), ".chr..OVERALL_MEAN") %in% levs_pasted))
  
  # check correct window notation,  step size
  expect_true(all(tsdc$single.window$position %% 200*100 == 0))
  expect_true(all(!tsdc$single.window$triple_sigma))
  expect_true(all(!tsdc$single.window$nk.status))
  expect_identical(tsdc$single.window[tsdc$single.window$facet == ".base" & tsdc$single.window$snp.facet == ".base",]$position, seq(0, 3000000, 200*1000))
  expect_warning(tsd2 <- calc_tajimas_d(.internal.data$test_snps, ".base", step = 400, triple_sigma = TRUE, sigma = 200), "Consider adding a snp level facet")
  tsdc2 <- get.snpR.stats(tsd2, ".base", "tajimas_d")
  expect_true(all(tsdc2$single.window$triple_sigma))
  
  # default bug
  expect_warning(tsd <- calc_tajimas_d(.internal.data$test_snps, "pop", triple_sigma = FALSE, sigma = 400), "Consider adding a snp level facet")
  tsd <- get.snpR.stats(tsd, "pop", "tajimas_d")
  expect_true(all(tsd$single.window$step == 800))
  
  # global
  tsd <- calc_tajimas_d(.internal.data$test_snps, tf, step = 200, triple_sigma = FALSE, sigma = 400, global = TRUE)
  tsdc <- get.snpR.stats(tsd, tf, "tajimas_d")
  levs <- unique(tsdc$weighted.means[,1:4])
  levs_pasted <- .paste.by.facet(levs, colnames(levs), sep = ".")
  sf <- summarize_facets(.internal.data$test_snps, tf)
  expect_true(".base..base..base..base" %in% levs_pasted)
  expect_true(all(paste0("pop.", lapply(strsplit(sf$chr.pop, "\\."), function(x){paste0(x[2], ".chr.", x[1])})) %in% levs_pasted))
  expect_true(all(paste0("fam.", lapply(strsplit(sf$chr.fam, "\\."), function(x){paste0(x[2], ".chr.", x[1])})) %in% levs_pasted))
  expect_true(all(paste0("fam.", names(sf$fam), "..base..base") %in% levs_pasted))
  expect_true(all(paste0("pop.", names(sf$pop), "..base..base") %in% levs_pasted))
  expect_true(all(paste0("fam.pop.", lapply(strsplit(sf$chr.fam.pop, "\\."), function(x){paste0(x[2], ".", x[3], ".chr.", x[1])})) %in% levs_pasted))
  expect_true(all(c("global_ws.theta", "global_ts.theta", "global_D", "global_num_seg") %in% colnames(tsdc$weighted.means)))
})

test_that("richness", {
  x <- calc_allelic_richness(stickSNPs[pop = c("ASP", "OPL")], "pop", g = -1)
  expect_true(all(unlist(.check_calced_stats(x, c("pop"), "richness"))))
  ar <- get.snpR.stats(x, "pop", "allelic_richness")
 
  expect_true("richness" %in% colnames(ar$single))
  expect_equal(round(ar$weighted.means$weighted_mean_richness, 3), 
               c(1.893, 1.885))
})

test_that("seg_sites", {
  # rarefaction
  x <- calc_seg_sites(stickSNPs[pop = c("ASP", "OPL")], "pop", g = -1)
  expect_true(all(unlist(.check_calced_stats(x, c("pop"), "seg_sites"))))
  s <- get.snpR.stats(x, "pop", "seg_sites")
  expect_equal(round(s$weighted.means$seg_sites, 3), 
               c(88.946, 88.321))
  
  ## other returns
  expect_true(all(c("g_prob_seg", "prob_seg", "prob_seg_var") %in% colnames(s$single)))
  expect_true("seg_sites_var" %in% colnames(s$weighted.means))
  
  # no rarefaction
  x <- calc_seg_sites(stickSNPs[pop = c("ASP", "OPL")], "pop", FALSE, g = -1)
  s <- get.snpR.stats(x, "pop", "seg_sites")
  expect_equal(s$weighted.means$seg_sites, 
               c(90, 96))
  
  # base facet, try to rarefact
  expect_warning(x <- calc_seg_sites(stickSNPs[pop = c("ASP", "OPL")]), "There is no reason to conduct rarefaction")
  expect_true(all(unlist(.check_calced_stats(x, ".base", "seg_sites"))))
  s <- get.snpR.stats(x, stats =  "seg_sites")
  expect_equal(s$weighted.means$seg_sites, 
               98)
  
  # base facet, don't try to rarefact
  x <- calc_seg_sites(stickSNPs[pop = c("ASP", "OPL")], rarefaction = FALSE)
  expect_true(all(unlist(.check_calced_stats(x, ".base", "seg_sites"))))
  s <- get.snpR.stats(x, stats =  "seg_sites")
  expect_equal(s$weighted.means$seg_sites, 
               98)
  
  # no rarefaction, with a locus with complete missing data
  d <- genotypes(stickSNPs)
  d[nrow(d), sample.meta(stickSNPs)$pop == "ASP"] <- "NN" # add a completely missing locus for ASP
  d <- import.snpR.data(d, sample.meta = sample.meta(stickSNPs))
  d <- calc_seg_sites(d[pop = c("ASP", "OPL")], "pop", FALSE)
  s <- get.snpR.stats(d, "pop", "seg_sites")
  s2 <- get.snpR.stats(calc_seg_sites(stickSNPs[pop = c("ASP", "OPL")], "pop", FALSE), "pop", "seg_sites")
  expect_equal(c(s$weighted.means$seg_sites[1], s2$weighted.means$seg_sites[1]), 
               c(89, 90))
})

