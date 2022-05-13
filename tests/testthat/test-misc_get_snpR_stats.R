test_that("fst_matrix return with facet containing no fst calcs", {
  dat <- calc_fis(.internal.data$test_snps, facets = c("pop", "pop.fam"))
  dat <- calc_pairwise_fst(dat, "pop") # runs, comparing Fst scores between pops
  res <- get.snpR.stats(dat, c("pop", "pop.fam"), c("fis", "fst"))
  
  expect_equal(names(res$fst.matrix), "pop")
  expect_equal(unique(res$pairwise$facet), "pop")
  expect_equal(unique(res$single$facet), c("fam.pop", "pop"))
})
