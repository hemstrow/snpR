context("subsetting")

#=======basic subsetting via index=========
test_that("index snps", {
  set.seed(1212)
  # snps
  sub.x <- .internal.data$test_snps[sample(10, 5, FALSE),]
  ## matching IDs
  str1 <- .paste.by.facet(snp.meta(sub.x), c("group", "position"))
  str2 <- .paste.by.facet(snp.meta(.internal.data$test_snps), c("group", "position"))
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
  
  # correct parts, simple
  ids <- id[pop = c("ASP", "PAL")]
  check <- sample.meta(ids)
  expect_equal(unique(check$pop), c("ASP", "PAL"))
  
  
  # correct parts, complex
  ids <- id[pop.fam = c("ASP.A", "PAL.B")]
  check <- sample.meta(ids)
  expect_equivalent(unique(check[,c(1:2)]), data.frame(pop = c("ASP", "PAL"),
                                                      fam = c("A", "B")))
})


test_that("snp facet",{
  id <- .internal.data$test_snps
  
  # error if using old syntax
  expect_error(id[snp.facet = "group", snp.subfacet = "groupIX"], "Facets and subfacets are now desginated directly using")
  
  # correct parts, simple
  ids <- id[group = c("groupIX", "groupIV")]
  check <- snp.meta(ids)
  expect_equal(unique(check$group), c("groupIX", "groupIV"))
})


test_that("complex facet",{
  id <- .internal.data$test_snps
  
  # correct parts, simple
  ids <- id[pop = c("ASP", "PAL"), group = c("groupIV")]
  check <- sample.meta(ids)
  expect_equal(unique(check$pop), c("ASP", "PAL"))
  check <- snp.meta(ids)
  expect_equal(unique(check$group), c("groupIV"))
  
  
  # correct parts, complex
  ids <- stickSNPs[pop.fam = c("ASP.A", "PAL.B"), group = c("groupIX", "groupIV", "groupXIX")]
  check <- sample.meta(ids)
  expect_equivalent(unique(check[,c(1:2)]), data.frame(pop = c("ASP", "PAL"),
                                                       fam = c("A", "B")))
  check <- snp.meta(ids)
  expect_equal(sort(unique(check$group)), sort(c("groupIX", "groupIV", "groupXIX")))
})