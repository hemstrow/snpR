context("subsetting")

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
  skip_on_cran(); skip_on_ci()
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
  ids <- id[chr = c("groupIX", "groupIV")]
  check <- snp.meta(ids)
  expect_equal(unique(check$chr), c("groupIX", "groupIV"))
})


test_that("complex facet",{
  id <- .internal.data$test_snps
  
  # correct parts, simple
  ids <- id[pop = c("ASP", "PAL"), chr = c("groupIV")]
  check <- sample.meta(ids)
  expect_equal(unique(check$pop), c("ASP", "PAL"))
  check <- snp.meta(ids)
  expect_equal(unique(check$chr), c("groupIV"))
  
  
  # correct parts, complex
  ids <- stickSNPs[pop.fam = c("ASP.A", "PAL.B"), chr = c("groupIX", "groupIV", "groupXIX")]
  check <- sample.meta(ids)
  expect_equivalent(unique(check[,c(1:2)]), data.frame(pop = c("ASP", "PAL"),
                                                       fam = c("A", "B")))
  check <- snp.meta(ids)
  expect_equal(sort(unique(check$chr)), sort(c("groupIX", "groupIV", "groupXIX")))
})


#========errors=========
test_that("errors",{
  
  id <- .internal.data$test_snps
  expect_error(id[0,1:10], regexp = "All requested snps must be within 1:nsnps")
  expect_error(id[1:10,0], regexp = "All requested snps must be within 1:nsnps(x)")
  expect_error(id[pop = "MAF"], regexp = "No sample found matching: pop -- MAF")
  expect_error(id[popl = "ASP"], regexp = "popl not found in x ")
  expect_error(id[popl = "ASP"], regexp = "popl not found in x ")
  
})