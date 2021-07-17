context("association")

set.seed(1212)
asdat <- .internal.data$test_snps
sample.meta(asdat)$phenotype <- rnorm(ncol(asdat))
sample.meta(asdat)$cat_phenotype <- sample(c("A", "B"), ncol(asdat), replace = TRUE)

test_that("correct gmmat", {
  local_edition(3)
  suppressWarnings(asgmmat <- calc_association(asdat, response = "phenotype"))
  asgmmat <- get.snpR.stats(asgmmat, stats = "association")
  expect_snapshot(asgmmat) # note, run off of gmmat, not internally calced. Thus checked, but should not change.
})

test_that("correct armitage", {
  local_edition(3)
  asarm <- calc_association(asdat, response = "cat_phenotype", method = "armitage")
  asarm <- get.snpR.stats(asarm, stats = "association")
  expect_snapshot(asarm) # Hand checked, should not change.
})

test_that("correct odds", {
  local_edition(3)
  asodds <- calc_association(asdat, response = "cat_phenotype", method = "odds_ratio")
  asodds <- get.snpR.stats(asodds, stats = "association")
  expect_snapshot(asodds) # Hand checked, should not change.
})


test_that("correct chisq", {
  local_edition(3)
  aschi <- calc_association(asdat, response = "cat_phenotype", method = "chisq")
  aschi <- get.snpR.stats(aschi, stats = "association")
  expect_snapshot(aschi) # Hand checked, should not change.
})