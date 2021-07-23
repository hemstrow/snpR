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
  
  # code to generate test values, not run to avoid CATT dependency
  # sn <- format_snps(asdat, "sn", interpolate = FALSE)[,-c(1:2)]
  # cc <- sample.meta(asdat)$cat_phenotype
  # cc <- as.numeric(as.factor(cc)) - 1
  # comp <- numeric(nrow(sn))
  # for(i in 1:nrow(sn)){
  #   nas <- which(is.na(sn[i,]))
  #   if(length(nas) > 0){
  #     comp[i] <- CATT::CATT(cc[-nas], sn[i,-nas])$p.value
  #   }
  #   else{
  #     comp[i] <- CATT::CATT(cc, sn[i,])$p.value
  #   }
  # }

  asarm <- get.snpR.stats(asarm, stats = "association")
  expect_equal(round(asarm$single$p_armitage_cat_phenotype, 4), 
               c(0.8402, 0.3894, 0.0651, 0.7469, 0.4652, 0.4386, 0.4076, 0.2888, 0.7782))# from CATT
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