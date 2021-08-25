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