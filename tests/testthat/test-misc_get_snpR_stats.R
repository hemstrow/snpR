test_that("fst_matrix return with facet containing no fst calcs", {
  dat <- calc_fis(.internal.data$test_snps, facets = c("pop", "pop.fam"))
  dat <- calc_pairwise_fst(dat, "pop") # runs, comparing Fst scores between pops
  res <- get.snpR.stats(dat, c("pop", "pop.fam"), c("fis", "fst"))
  
  expect_equal(names(res$fst.matrix), "pop")
  expect_equal(unique(res$pairwise$facet), "pop")
  expect_equal(unique(res$single$facet), c("fam.pop", "pop"))
})

test_that("fst_matrix always in upper triangle",{
  # add a dataset that was known to bork a bit
  d <- genotypes(stickSNPs)
  ptab <- data.frame(new = c("San_Luis_Mesa", "Rafael_Swell", "Seep_Ridge", "Soap_Creek", "SOS_CP292", "SOS_CP298"),
                     old = unique(sample.meta(stickSNPs)$pop))
  nmet <- sample.meta(stickSNPs)
  nmet$pop <- ptab$new[match(nmet$pop, ptab$old)]
  d <- import.snpR.data(d, snp.meta(stickSNPs), nmet)
  
  d <- calc_pairwise_fst(d, c("pop", "fam", "pop.fam"), boot = 5)
  
  # check that data is in the correct triangle of the matrix
  res <- get.snpR.stats(d, c("pop", "fam", "pop.fam"), "fst")
  upper_only_OK <- function(res, try_up = TRUE){
    c1 <- all(is.na(res[lower.tri(res)]))
    if(try_up){
      c2 <- all(!is.na(res[upper.tri(res)]))
      return(c(c1, c2))
    }
    return(c1)
  }
  
  ## fst
  expect_true(all(upper_only_OK(as.matrix(res$fst.matrix$pop$fst[,-1]))))
  expect_true(all(upper_only_OK(as.matrix(res$fst.matrix$fam$fst[,-1]))))
  expect_true(all(upper_only_OK(as.matrix(res$fst.matrix$fam.pop$fst[,-1]), 
                                try_up = FALSE))) # upper should correctly have a few NAs
  
  ## p-values
  expect_true(all(upper_only_OK(as.matrix(res$fst.matrix$pop$p[,-1]))))
  expect_true(all(upper_only_OK(as.matrix(res$fst.matrix$fam$p[,-1]))))
  expect_true(all(upper_only_OK(as.matrix(res$fst.matrix$fam.pop$p[,-1]), 
                                try_up = FALSE))) # upper should correctly have a few NAs
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
  
  
  expect_error(calc_pi(test, "pop.dup_test"), "This will cause issues if those facets are run during analysis.+Level: ASP\tin facets: dup_test, pop")
  expect_error(calc_pi(test, "pop.dup_test2"), "This will cause issues if those facets are run during analysis.+Level: ASP\tin facets: dup_test2, pop")
  
})
