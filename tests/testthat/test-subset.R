#=======basic subsetting via index=========
test_that("index snps", {
  set.seed(1212)
  # snps
  sub.x <- .internal.data$test_snps[sample(10, 5, FALSE),]
  ## matching IDs
  str1 <- .paste.by.facet(snp.meta(sub.x), c("chr", "position"))
  str2 <- .paste.by.facet(snp.meta(.internal.data$test_snps), c("chr", "position"))
  expect_equal(snp.meta(sub.x)$.snp.id, snp.meta(.internal.data$test_snps)[match(str1, str2),]$.snp.id)
})


test_that("index samps", {
  set.seed(1212)
  id <- .internal.data$test_snps
  sample.meta(id)$id <- 1:nsamps(id)
  # snps
  sub.x <- id[,sample(10, 5, FALSE)]
  ## matching IDs
  str1 <- .paste.by.facet(sample.meta(sub.x), c("pop", "fam", "id"))
  str2 <- .paste.by.facet(sample.meta(id), c("pop", "fam", "id"))
  expect_equal(sample.meta(sub.x)$.sample.id, sample.meta(.internal.data$test_snps)[match(str1, str2),]$.sample.id)
  
})

#========by reference=========
test_that("sample facet",{
  id <- .internal.data$test_snps
  
  # error if using old syntax
  expect_error(id[facet = "pop", subfacet = "ASP"], "Facets and subfacets are now desginated directly using")
  
  # correct parts, very simple
  ids <- id[pop = "ASP"]
  check <- sample.meta(ids)
  expect_equal(unique(check$pop), "ASP")
  
  # correct parts, simple
  ids <- id[pop = c("ASP", "PAL")]
  check <- sample.meta(ids)
  expect_equal(unique(check$pop), c("ASP", "PAL"))
  
  
  # correct parts, complex
  ids <- id[pop.fam = c("ASP.A", "PAL.B")]
  check <- sample.meta(ids)
  expect_equivalent(unique(check[,c(1:2)]), data.frame(pop = c("ASP", "PAL"),
                                                      fam = c("A", "B")))
  
  
  # correct parts, via reference to environmental variables
  skip_on_cran();
  # need to assign up to the global environment, otherwise it won't find pops when testing. 
  # Very not ideal, but with the way testthat works I can't think of a better way. Skipped on cran for this reason.
  # This should always work with a normal script or interactive use, where objects are set to the global environment automatically.
  pops <<- data.frame(pop = c("ASP", "PAL"), fam = c("A", "B")) 
  
  tdat <- id[pop.fam = paste0(pops[, 1], ".", pops[, 2])]
  rm(pops, envir = globalenv()) # clean from global afterwards, don't actually want to pass this up out of the test
  
  expect_identical(ids, tdat)
})


test_that("snp facet",{
  id <- .internal.data$test_snps
  
  # error if using old syntax
  expect_error(id[snp.facet = "chr", snp.subfacet = "groupIX"], "Facets and subfacets are now desginated directly using")
  
  # correct parts, simple
  ids <- id[chr = c("groupIX", "groupXIX")]
  check <- snp.meta(ids)
  expect_equal(unique(check$chr), c("groupIX", "groupXIX"))
})


test_that("complex facet",{
  id <- .internal.data$test_snps
  
  # correct parts, simple
  ids <- id[pop = c("ASP", "PAL"), chr = c("groupXIX")]
  check <- sample.meta(ids)
  expect_equal(unique(check$pop), c("ASP", "PAL"))
  check <- snp.meta(ids)
  expect_equal(unique(check$chr), c("groupXIX"))
  
  
  # correct parts, complex
  ids <- stickSNPs[pop.fam = c("ASP.A", "PAL.B"), chr = c("groupIX", "groupXIX", "groupXII")]
  check <- sample.meta(ids)
  expect_equivalent(unique(check[,c(1:2)]), data.frame(pop = c("ASP", "PAL"),
                                                       fam = c("A", "B")))
  check <- snp.meta(ids)
  expect_equal(sort(unique(check$chr)), sort(c("groupIX", "groupXII", "groupXIX")))
})


#========errors=========
test_that("errors",{
  
  id <- .internal.data$test_snps
  expect_error(id[0,1:10], regexp = "All requested snps must be within 1:nsnps")
  expect_error(id[1:10,0], regexp = "All requested samples must be within 1:nsnps")
  expect_error(id[pop = "MAF"], regexp = "No sample found matching: pop -- MAF")
  expect_error(id[popl = "ASP"], regexp = "popl not found in x ")
  expect_error(id[popl = "ASP"], regexp = "popl not found in x ")
  
})


#========updating=======
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
  y <- x
  
  # adding a new col doesn't change anything, but does add a col to the sample.stats
  sample.meta(y)$new <- "test"
  expect_identical(x@stats, y@stats)
  expect_true("new" %in% colnames(y@sample.stats))
  
  # altering a non-facet col doesn't change anything except for sample.stats
  sample.meta(y)$new <- "hi"
  expect_identical(x@stats, y@stats)
  if("new" %in% colnames(y@sample.stats)){
    expect_true(all(y@sample.stats$new == "hi"))
  }
  
  
  # removing an old col causes that col to be purged from the data
  y <- x
  sample.meta(y)$pop <- NULL
  expect_false("pop" %in% 
                 unlist(.split.facet(c(y@facet.meta$facet, y@facets, 
                                       y@stats$facet, y@pairwise.stats$facet, 
                                       y@window.stats$facet, y@pairwise.window.stats$facet,
                                       y@pop.stats$facet, y@pairwise.LD$prox$sample.facet,
                                       y@window.bootstraps$facet,
                                       unlist(.split.facet(names(y@calced_stats))),
                                       names(y@allele_frequency_matrices),
                                       names(y@genetic_distances),
                                       y@weighted.means$facet,
                                       names(y@other$geo_dists),
                                       names(y@other$ibd)))))
  
  expect_false("pop" %in% colnames(y@sample.stats))
  
  # similar result if we change that col, just with an update in sample.stats instead of a removal
  y <- x
  sample.meta(y)$pop <- rep(c("test", "PAL"), length.out = ncol(x))
  expect_false("pop" %in% 
                 unlist(.split.facet(c(y@facet.meta$facet, y@facets, 
                                       y@stats$facet, y@pairwise.stats$facet, 
                                       y@window.stats$facet, y@pairwise.window.stats$facet,
                                       y@pop.stats$facet, y@pairwise.LD$prox$sample.facet,
                                       y@window.bootstraps$facet,
                                       unlist(.split.facet(names(y@calced_stats))),
                                       names(y@allele_frequency_matrices),
                                       names(y@genetic_distances),
                                       y@weighted.means$facet,
                                       names(y@other$geo_dists),
                                       names(y@other$ibd)))))
  expect_true("pop" %in% colnames(y@sample.stats))
  if("pop" %in% colnames(y@sample.stats)){
    expect_true(all(y@sample.stats$pop == rep(c("test", "PAL"), length.out = ncol(x))))
  }
  
  
  # error with identical
  nmeta <- sample.meta(x)
  colnames(nmeta)[1] <- "chr"
  expect_error(sample.meta(x) <- nmeta, "Some unacceptable column names")
  
})


