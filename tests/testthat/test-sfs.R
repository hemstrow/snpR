#========basic sfs==========
test_that("sfs",{
  # 1D
  .make_it_quiet(sfs <- calc_sfs(stickSNPs, projection = 90))
  expect_true(sum(sfs, na.rm = T) <= nsnps(stickSNPs))
  expect_equal(attr(sfs, "pop"), ".base")
  expect_true(max(sfs, na.rm = T) <= 90/2) # folded
  expect_true(attr(sfs, "folded"))
  expect_equal(length(sfs), 45)
  
  # 2D
  .make_it_quiet(sfs2 <- calc_sfs(stickSNPs, facet = "pop", pops =  c("ASP", "UPD"), projection = c(10, 10)))

  expect_true(sum(round(sfs2, 5), na.rm = T) <= nsnps(stickSNPs)) # round because of rounding issues on MAC
  expect_equal(attr(sfs2, "pop"), c("ASP", "UPD"))
  expect_true(attr(sfs2, "folded"))
  expect_equal(dim(sfs2), c(11, 11))
  
  # unfolded
  expect_warning(.make_it_quiet(sfs <- calc_sfs(stickSNPs, projection = 90, fold = FALSE)), "ithout ancestral and derived character states, unfolded spectra will be misleading")
  expect_false(attr(sfs, "folded"))
  expect_equal(length(sfs), 91)
  
  expect_warning(.make_it_quiet(sfs2 <- calc_sfs(stickSNPs, facet = "pop", pops =  c("ASP", "UPD"), projection = c(10, 10), fold = FALSE)), "ithout ancestral and derived character states, unfolded spectra will be misleading")
  expect_false(attr(sfs2, "folded"))
  expect_equal(dim(sfs2), c(11, 11))
})

#========directionality=========
test_that("directionality",{
  # no provided sfs
  expect_error(.make_it_quiet(res <- calc_directionality(stickSNPs, facet = "pop", pops = c("ASP", "CLF"), projection = c(50, 50)), "No segrgating sites remain after projection."))
  expect_warning(.make_it_quiet(res <- calc_directionality(stickSNPs, facet = "pop", pops = c("ASP", "CLF"), projection = c(10, 10))), "Without ancestral and derived character states, unfolded spectra will be misleading.")
  
  
  
  expect_equal(as.numeric(round(res, 4)), -0.0203) # not meaningful here because of lack of polarization
  expect_equal(attr(res, "direction"), "ASP<-CLF")
  
  # provided sfs
  expect_warning(.make_it_quiet(sfs <- calc_sfs(stickSNPs, facet = "pop", pops = c("ASP", "CLF"), projection = c(10, 10), fold = FALSE)), "Without ancestral and derived character states, unfolded spectra will be misleading.")
  .make_it_quiet(res2 <- calc_directionality(sfs))
  expect_identical(res, res2)
  
  # sanity checks
  .make_it_quiet(sfs2 <- calc_sfs(stickSNPs, facet = "pop", pops =  c("ASP", "UPD"), projection = c(10, 10)))
  expect_error(calc_directionality(sfs2), "This most likely means that the SFS is folded")
  expect_warning(.make_it_quiet(sfs2 <- calc_sfs(stickSNPs, projection = c(10), fold = FALSE)), "Without ancestral and derived character states, unfolded spectra will be misleading.")
  expect_error(calc_directionality(sfs2), "SFS dimensionality 1 not allowed for this function. Allowed dimensionalities: 2")
})

#========origin of expansion=========
test_that("origin-of-expansion",{
  # NOTE: this test is slow and may need to be suppressed on CRAN.
  
  # setup
  # set ref and anc--ideally use an outgroup for this
  dat <- calc_maf(stickSNPs)
  snp.meta(dat)$ref <- paste0("A", get.snpR.stats(dat)$minor, "A")
  expect_warning(snp.meta(dat)$anc <- paste0("A", get.snpR.stats(dat)$major, "A"))
  expect_warning(dat <- dat[pop = c("ASP", "PAL", "CLF")])
  
  # setup x and y coords
  long_lat <- data.frame(CLF = c(44.267718, -121.255805),
                         ASP = c(43.891693, -121.448360),
                         PAL = c(43.714114, -121.272797))
  long_lat <- t(long_lat)
  long_lat <- long_lat[match(sample.meta(dat)$pop, rownames(long_lat)),]
  colnames(long_lat) <- c("y", "x")
  expect_warning(sample.meta(dat) <- cbind(sample.meta(dat), long_lat))
  
  

  projection <- summarize_facets(dat, "pop")[["pop"]]
  projection <- floor(projection*.8)
  
  # basic test with a few boots
  out <- calc_origin_of_expansion(dat, "pop", boots = 4, projection = projection, boot_par = 2)
  expect_true(all(names(out) == c("opt", "pairwise_directionality")))
  expect_true(all(names(out$opt) == c("v", "x", "y")))
  expect_true(all(colnames(out$pairwise_directionality) == c("Directionality", "Variance", "xi", "yi", "xj", "yj", "comparison")))
  expect_true(all(out$pairwise_directionality$comparison == c("ASP~CLF", "ASP~PAL", "CLF~PAL")))
  
  # test that serial works as well
  out_s <- calc_origin_of_expansion(dat, "pop", projection = projection, boots = 4)
  expect_true(all(names(out_s) == c("opt", "pairwise_directionality")))
  expect_true(all(names(out_s$opt) == c("v", "x", "y")))
  expect_true(all(colnames(out_s$pairwise_directionality) == c("Directionality", "Variance", "xi", "yi", "xj", "yj", "comparison")))
  expect_true(all(out_s$pairwise_directionality$comparison == c("ASP~CLF", "ASP~PAL", "CLF~PAL")))
  expect_true(all(out$pairwise_directionality$Directionality == out_s$pairwise_directionality$Directionality))
  
  # errors
  expect_error(calc_origin_of_expansion(dat, c("pop", "fam"), boots = 2, projection = projection),
               "one facet at a time")
  expect_error(calc_origin_of_expansion(dat, ".base", boots = 2, projection = projection),
               "sample facet must be provided")
  expect_error(calc_origin_of_expansion(dat, "pop", boots = 2, projection = c(ASP = 40, PAL = 13, CLF = 40)),
               "projections must be smaller than or.+ASP.+CLF")
  expect_error(calc_origin_of_expansion(dat, "pop", boots = 2, projection = c(ASP = 10, CLF = 8)),
               "facet levels are missing in the projection vector.+PAL\n$")
  expect_error(.suppress_specific_warning(calc_origin_of_expansion(dat[pop = c("ASP", "CLF")], "pop", boots = 2, projection = projection), "Some levels are duplicated"),
               "At least three populations")
  
  et1 <- dat
  suppressWarnings((sample.meta(et1)$x <- NULL))
  expect_error(calc_origin_of_expansion(et1, "pop", boots = 2, projection = projection),
               "'x' and 'y' columns.+must be present")
  
  et2 <- et1
  suppressWarnings((snp.meta(et2)$anc <- NULL))
  expect_warning(expect_error(calc_origin_of_expansion(et2, "pop", boots = 2, projection = projection),
               "'x' and 'y' columns.+must be present"), "Without ancestral")
  
  et3 <- dat
  snp.meta(et3)$ref <- substr(snp.meta(et3)$ref, 1, 2)
  expect_error(calc_origin_of_expansion(et3, c("pop"), boots = 2, projection = projection),
               "ref and anc entries must be exactly three characters long")
})

