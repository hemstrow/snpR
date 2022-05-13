#========basic sfs==========
test_that("sfs",{
  # 1D
  .make_it_quiet(sfs <- calc_sfs(stickSNPs, projection = 100))
  expect_true(sum(sfs, na.rm = T) <= nsnps(stickSNPs))
  expect_equal(attr(sfs, "pop"), ".base")
  expect_true(max(sfs, na.rm = T) <= 100/2) # folded
  
  # 2D
  .make_it_quiet(sfs2 <- calc_sfs(stickSNPs, facet = "pop", pops =  c("ASP", "UPD"), projection = c(10, 10)))
  expect_true(sum(sfs2, na.rm = T) <= nsnps(stickSNPs))
  expect_equal(attr(sfs2, "pop"), c("ASP", "UPD"))
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

