test_that("sample meta fetching and updating",{
  
  # fetching
  expect_identical(stickSNPs@sample.meta, sample.meta(stickSNPs))
  
  # replacement
  x <- .internal.data$test_snps
  sample.meta(x) <- cbind.data.frame(sample.meta(x), 
                                     x = rnorm(ncol(x)),
                                     y = rnorm(ncol(x)))
  tf <- c("pop", "pop.chr", "pop.chr", "fam", ".base", "pop.fam")
  x <- calc_pi(x, tf)
  x <- calc_pairwise_fst(x, tf[-5])
  x <- calc_pairwise_ld(x, tf)
  x <- .suppress_specific_warning(calc_isolation_by_distance(x, tf), "Mantel tests for IBD")
  x <- calc_hs(x, tf)
  x <- calc_smoothed_averages(x, tf[-5], sigma =  100)
  x <- do_bootstraps(x, tf[-5], 5, sigma = 100)
  smx <- sample.meta(x)
  smx$fam <- gsub("A", "TEST", smx$fam)
  smx$new <- "test"
  y <- x
  
  # adding a new col doesn't change anything
  sample.meta(y)$new <- "test"
  expect_identical(x@stats, y@stats)
  
  # altering a non-facet col doesn't change anything
  sample.meta(y)$new <- "hi"
  expect_identical(x@stats, y@stats)
  
  # altering an old col causes stats to be lost
  y <- x
  sample.meta(y)$pop <- NULL
  eypect_false("fam" %in% c(y@facet.meta$facet, y@facets, 
                            y@stats$facet, y@pairwise.stats$facet, 
                            y@window.stats$facet, y@pairwise.window.stats$facet,
                            colnames(y@sample.stats), y@pop.stats$facet, y@pairwise.LD$prox$sample.facet,
                            y@window.bootstraps$facet,
                            unlist(.split.facet(names(y@calced_stats)))
                            ))
  
  
})