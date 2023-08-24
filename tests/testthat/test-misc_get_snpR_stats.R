test_that("fst_matrix return with facet containing no fst calcs", {
  dat <- calc_fis(.internal.data$test_snps, facets = c("pop", "pop.fam"))
  dat <- calc_pairwise_fst(dat, "pop") # runs, comparing Fst scores between pops
  res <- get.snpR.stats(dat, c("pop", "pop.fam"), c("fis", "fst"))
  
  expect_equal(names(res$fst.matrix), "pop")
  expect_equal(unique(res$pairwise$facet), "pop")
  expect_equal(unique(res$single$facet), c("fam.pop", "pop"))
})

test_that("empty return warning", {
  expect_warning(get.snpR.stats(stickSNPs, "pop", "fst"), "statistics located for requested stats/facets")
})

test_that("summarize_facets",{
  expect_message(summarize_facets(stickSNPs), "Returning list of facets")
  
  # basic reporting
  expect_identical(suppressMessages(summarize_facets(stickSNPs)),
                   list(SNP = c("chr", "position", ".snp.id"),
                        sample = c("pop", "fam", ".sample.id")))
  
  # facet reporting
  check_sum <- summarize_facets(stickSNPs, c("pop", "pop.chr", "chr", "pop.chr.fam"))
  expect_equal(check_sum$pop, table(sample.meta(stickSNPs)$pop))
  expect_equal(check_sum$chr, unique(snp.meta(stickSNPs)$chr))
  chr_pop_check <- unique(.paste.by.facet(expand.grid(list(unique(snp.meta(stickSNPs)$chr), 
                                                           unique(sample.meta(stickSNPs)$pop))), 
                                          c("Var1", "Var2")))
  expect_equal(sort(check_sum$chr.pop), sort(chr_pop_check))
  chr_pop_fam_check <- unique(.paste.by.facet(expand.grid(list(unique(snp.meta(stickSNPs)$chr),
                                                               unique(sample.meta(stickSNPs)$fam),
                                                               unique(sample.meta(stickSNPs)$pop))),
                                          c("Var1", "Var2", "Var3")))
  expect_equal(sort(check_sum$chr.fam.pop), sort(chr_pop_fam_check))
  
})

test_that("bad facet request with duplicates",{
  
  meta <- sample.meta(stickSNPs)
  meta$dup_test <- "ASP"
  
  snpm <- snp.meta(stickSNPs)
  snpm$dup_test2 <- "ASP"
  expect_warning(test <- import.snpR.data(genotypes(stickSNPs), sample.meta = meta, snp.meta = snpm))
  
  
  expect_error(calc_pi(test, "pop.dup_test"), "This will cause issues if those facets are run during analysis.+Level: ASP\tin facets: dup_test, pop1")
  expect_error(calc_pi(test, "pop.dup_test2"), "This will cause issues if those facets are run during analysis.+Level: ASP\tin facets: dup_test2, pop1")
  
})
