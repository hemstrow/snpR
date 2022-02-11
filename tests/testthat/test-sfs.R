context("sfs, directionality")


#========basic sfs==========
test_that("sfs",{
  # 1D
  expect_warning(sfs <- calc_sfs(stickSNPs, projection = 100), "ref and anc columns are suggested in snp metadata. See documentation for details. The major allele will be subsituted for the ancestral state.")
  expect_true(sum(sfs, na.rm = T) <= nsnps(stickSNPs))
  expect_equal(attr(sfs, "pop"), ".base")
  expect_true(max(sfs, na.rm = T) <= 100/2) # folded
  
  # 2D
  expect_warning(sfs2 <- calc_sfs(stickSNPs, facet = "pop", pops =  c("ASP", "UPD"), projection = c(50, 50)), "ref and anc columns are suggested in snp metadata. See documentation for details. The major allele will be subsituted for the ancestral state.")
  expect_true(sum(sfs2, na.rm = T) <= nsnps(stickSNPs))
  expect_equal(attr(sfs2, "pop"), c("ASP", "UPD"))
})

#========directionality=========
test_that("directionality",{
  # no provided sfs
  expect_warning(res <- calc_directionality(stickSNPs, facet = "pop", pops = c("ASP", "CLF"), projection = c(50, 50)), "Without ancestral and derived character states, unfolded spectra will be misleading.")
  
  expect_equal(as.numeric(round(res, 4)), -0.0253) # not meaniful here because of lack of polarization
  expect_equal(attr(res, "direction"), "ASP<-CLF")
})

