context("association")

set.seed(1212)
asdat <- .internal.data$test_snps
sample.meta(asdat)$phenotype <- rnorm(ncol(asdat))
sample.meta(asdat)$cat_phenotype <- sample(c("A", "B"), ncol(asdat), replace = TRUE)

#=========association=========
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


#=========random forest========
test_that("random forest",{
  skip_if_not_installed("ranger")
  # basic
  rf <- run_random_forest(asdat, response = "phenotype", pvals = FALSE)
  expect_s3_class(rf$models$.base_.base$model, "ranger")
  expect_equal(dim(rf$models$.base_.base$predictions), c(nsamps(asdat), 2))
  expect_equal(length(rf$models), 1)

  # several facets
  rf <- run_random_forest(asdat, facets = "pop", response = "phenotype", pvals = FALSE)
  expect_equal(length(rf$models), 2)
  expect_s3_class(rf$models$pop_ASP$model, "ranger")
  expect_s3_class(rf$models$pop_PAL$model, "ranger")
  expect_equal(dim(rf$models$pop_ASP$predictions), c(5, 2))
  expect_equal(dim(rf$models$pop_PAL$predictions), c(5, 2))
  
  # formula specification
  rf <- run_random_forest(asdat, response = "cat_phenotype", formula = cat_phenotype ~ phenotype, pvals = FALSE)
  expect_equal(rf$models$.base_.base$covariate_importance$variable, "phenotype")
  
  # pvals
  rf <- run_random_forest(asdat, response = "cat_phenotype", formula = cat_phenotype ~ phenotype, pvals = TRUE)
  expect_true("p_val" %in% colnames(rf$models$.base_.base$covariate_importance))
  expect_true("cat_phenotype_RF_importance_pvals" %in% colnames(get.snpR.stats(rf$data, stats = "random_forest")$single))
  
  # note: all rf construction and calculation external, no need to check the numbers themselves.
})


#========genomic prediction=========
test_that("genomic prediction",{
  
  
  
  
})
