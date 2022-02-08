context("association")

set.seed(1212)
asdat <- .internal.data$test_snps
sample.meta(asdat)$phenotype <- rnorm(ncol(asdat))
sample.meta(asdat)$cat_phenotype <- sample(c("A", "B"), ncol(asdat), replace = TRUE)

#=========association=========
test_that("correct gmmat", {
  local_edition(3)
  skip_on_cran()
  suppressWarnings(asgmmat <- calc_association(asdat, response = "phenotype"))
  asgmmat <- get.snpR.stats(asgmmat, stats = "association")
  expect_snapshot(asgmmat) # note, run off of gmmat, not internally calced. Thus checked, but should not change.
})

test_that("correct armitage", {
  local_edition(3)
  skip_if_not_installed("CATT")
  asarm <- calc_association(asdat, response = "cat_phenotype", method = "armitage")
  
  #code to generate test values
  sn <- format_snps(asdat, "sn", interpolate = FALSE)[,-c(1:2)]
  cc <- sample.meta(asdat)$cat_phenotype
  cc <- as.numeric(as.factor(cc)) - 1
  comp <- numeric(nrow(sn))
  for(i in 1:nrow(sn)){
    nas <- which(is.na(sn[i,]))
    if(length(nas) > 0){
      comp[i] <- CATT::CATT(cc[-nas], sn[i,-nas])$p.value
    }
    else{
      comp[i] <- CATT::CATT(cc, sn[i,])$p.value
    }
  }

  asarm <- get.snpR.stats(asarm, stats = "association")
  expect_equal(round(asarm$single$p_armitage_cat_phenotype, 4), 
               comp) # from CATT
})

test_that("correct odds", {
  local_edition(3)
  skip_on_cran()
  asodds <- calc_association(asdat, response = "cat_phenotype", method = "odds_ratio")
  asodds <- get.snpR.stats(asodds, stats = "association")
  expect_snapshot(asodds) # Hand checked, should not change.
})


test_that("correct chisq", {
  local_edition(3)
  skip_on_cran()
  aschi <- calc_association(asdat, response = "cat_phenotype", method = "chisq")
  aschi <- get.snpR.stats(aschi, stats = "association")
  expect_snapshot(aschi) # Hand checked, should not change.
})


#=========random forest========
test_that("random forest",{
  skip_if_not_installed("ranger")
  # basic
  rf <- run_random_forest(asdat, response = "phenotype", pvals = FALSE)
  rfstats <- get.snpR.stats(rf$x, stats = "random_forest")
  expect_s3_class(rf$models$.base_.base$model, "ranger")
  expect_equal(dim(rf$models$.base_.base$predictions), c(nsamps(asdat), 2))
  expect_equal(length(rf$models), 1)
  expect_equal(unique(rfstats$single$subfacet), c(".base"))
  expect_equal(unique(rfstats$single$facet), c(".base"))
  expect_equal(colnames(rfstats$single), c("facet", "subfacet", "chr", "position", "phenotype_RF_importance"))

  # several facets
  rf <- run_random_forest(asdat, facets = "pop", response = "phenotype", pvals = FALSE)
  rfstats <- get.snpR.stats(rf$x, "pop", "random_forest")
  expect_equal(length(rf$models), 2)
  expect_s3_class(rf$models$pop_ASP$model, "ranger")
  expect_s3_class(rf$models$pop_PAL$model, "ranger")
  expect_equal(dim(rf$models$pop_ASP$predictions), c(5, 2))
  expect_equal(dim(rf$models$pop_PAL$predictions), c(5, 2))
  expect_equal(unique(rfstats$single$subfacet), c("ASP", "PAL"))
  expect_equal(unique(rfstats$single$facet), c("pop"))
  expect_true("phenotype_RF_importance" %in% colnames(rfstats$single))
  expect_equal(colnames(rfstats$single), c("facet", "subfacet", "chr", "position", "phenotype_RF_importance"))
  str <- .paste.by.facet(rfstats$single, c("chr", "position"))
  expect_equal(as.numeric(table(str)), rep(2, nrow(asdat))) # each snp has calcs for each pop
  
  
  
  # formula specification
  rf <- run_random_forest(asdat, response = "cat_phenotype", formula = cat_phenotype ~ phenotype, pvals = FALSE)
  expect_equal(rf$models$.base_.base$covariate_importance$variable, "phenotype")
  
  # pvals
  expect_warning(rf <- run_random_forest(asdat, response = "cat_phenotype", formula = cat_phenotype ~ phenotype, pvals = TRUE), "Consider the 'altmann' approach.")
  expect_true("p_val" %in% colnames(rf$models$.base_.base$covariate_importance))
  expect_true("cat_phenotype_RF_importance_pvals" %in% colnames(get.snpR.stats(rf$x, stats = "random_forest")$single))
  rfstats <- get.snpR.stats(rf$x, stats = "random_forest")
  expect_true("cat_phenotype_RF_importance_pvals" %in% colnames(rfstats$single))

  # note: all rf construction and calculation external, no need to check the numbers themselves.
})


#========genomic prediction=========
test_that("genomic prediction",{
  skip_if_not_installed("BGLR")
  
  # single
  gp <- run_genomic_prediction(asdat, response = "phenotype", iterations = 200, burn_in = 100, thin = 10)
  gpstats <- get.snpR.stats(gp$x, stats = "genomic_prediction")
  expect_s3_class(gp$models$.base_.base$model, "BGLR")
  expect_equal(dim(gp$models$.base_.base$predictions), c(nsamps(asdat), 2))
  expect_equal(length(gp$models), 1)
  expect_equal(unique(gpstats$single$subfacet), c(".base"))
  expect_equal(unique(gpstats$single$facet), c(".base"))
  expect_equal(colnames(gpstats$single), c("facet", "subfacet", "chr", "position", "phenotype_gp_effect"))

  # with facets
  gp <- run_genomic_prediction(asdat, facets = "pop", response = "phenotype", iterations = 200, burn_in = 100, thin = 10)
  gpstats <- get.snpR.stats(gp$x, "pop", "genomic_prediction")
  expect_equal(length(gp$models), 2)
  expect_s3_class(gp$models$pop_ASP$model, "BGLR")
  expect_s3_class(gp$models$pop_PAL$model, "BGLR")
  expect_equal(dim(gp$models$pop_ASP$predictions), c(5, 2))
  expect_equal(dim(gp$models$pop_PAL$predictions), c(5, 2))
  expect_equal(unique(gpstats$single$subfacet), c("ASP", "PAL"))
  expect_equal(unique(gpstats$single$facet), c("pop"))
  expect_equal(colnames(gpstats$single), c("facet", "subfacet", "chr", "position", "phenotype_gp_effect"))
  str <- .paste.by.facet(gpstats$single, c("chr", "position"))
  expect_equal(as.numeric(table(str)), rep(2, nrow(asdat))) # each snp has calcs for each pop
  
  
  # runing on a previously run set
  gp <- run_genomic_prediction(asdat, response = "phenotype", iterations = 200, burn_in = 100, thin = 10)
  gp <- run_genomic_prediction(gp$x, facets = "pop", response = "phenotype", iterations = 200, burn_in = 100, thin = 10)
  gpstats <- get.snpR.stats(gp$x, facets = "pop", stats = "genomic_prediction")
  str <- .paste.by.facet(gpstats$single, c("chr", "position"))
  expect_equal(as.numeric(table(str)), rep(2, nrow(asdat))) # each snp has calcs for each pop
  
  gpstats <- get.snpR.stats(gp$x, stats = "genomic_prediction")
  str <- .paste.by.facet(gpstats$single, c("chr", "position"))
  expect_equal(as.numeric(table(str)), rep(1, nrow(asdat))) # each snp has calcs for each pop
  
  
  # errors
  expect_error(run_genomic_prediction(gp$x, response = "weight", iterations = 200, burn_in = 100, thin = 10), regexp =  "column matching response 'weight' found in sample metadata")
  # note: all gp construction and calculation external, no need to check the numbers themselves.
})

test_that("genomic prediction CV",{
  skip_if_not_installed("BGLR")
  skip_on_cran()
  
  set.seed(1212)
  dat <- stickSNPs
  sample.meta(dat)$phenotype <- rnorm(ncol(dat))
  sample.meta(dat)$cat_phenotype <- sample(c("A", "B"), ncol(dat), replace = TRUE)
  ## run cross_validation
  res <-cross_validate_genomic_prediction(dat, response = "phenotype", iterations = 200, burn_in = 100, thin = 10)
  
  expect_s3_class(res$model$models$.base_.base$model, "BGLR")
  expect_equal(length(res$model.samples), floor(nsamps(stickSNPs)*.9))
  expect_false(any(res$model.samples %in% res$cross.samples))
  expect_equal(colnames(res$comparison), c("observed", "predicted"))
  expect_equal(nrow(res$comparison), length(res$cross.samples))
  expect_equal(cor(res$comparison$observed, res$comparison$predicted)^2, res$rsq)
  
})
