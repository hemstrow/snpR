context("sfs, directionality")


#========basic sfs==========
test_that("sfs",{
  # 1D
  sfs <- calc_sfs(stickSNPs, projection = 100)
  expect_true(sum(sfs, na.rm = T) <= nsnps(stickSNPs))
  expect_equal(attr(sfs, "pop"), ".base")
  expect_true(max(sfs, na.rm = T) <= 100/2) # folded
  
  # 2D
  sfs2 <- calc_sfs(stickSNPs, facet = "pop", pops =  c("ASP", "UPD"), projection = c(50, 50))
  expect_true(sum(sfs2, na.rm = T) <= nsnps(stickSNPs))
  expect_equal(attr(sfs2, "pop"), c("ASP", "UPD"))
})

#========directionality=========
test_that("directionality",{
  # no provided sfs
  expect_warning(res <- calc_directionality(stickSNPs, facet = "pop", pops = c("ASP", "CLF"), projection = c(50, 50)), "Without ancestral and derived character states, unfolded spectra will be misleading.")
  
  expect_equal(as.numeric(round(res, 4)), -0.0253) # not meaniful here because of lack of polarization
  expect_equal(attr(res, "direction"), "ASP<-CLF")
  
  # provided sfs
  sfs <-  expect_warning(calc_sfs(stickSNPs, facet = "pop", pops = c("ASP", "CLF"), projection = c(50, 50), fold = FALSE), "Without ancestral and derived character states, unfolded spectra will be misleading.")
  res2 <- calc_directionality(sfs)
  expect_identical(res, res2)
  
  # sanity checks
  sfs2 <- calc_sfs(stickSNPs, facet = "pop", pops =  c("ASP", "UPD"), projection = c(50, 50))
  expect_error(calc_directionality(sfs2), "This most likely means that the SFS is folded")
  expect_warning(sfs2 <- calc_sfs(stickSNPs, projection = c(50), fold = FALSE), "Without ancestral and derived character states, unfolded spectra will be misleading.")
  expect_error(calc_directionality(sfs2), "SFS dimensionality 1 not allowed for this function. Allowed dimensionalities: 2")
})

