set.seed(1212)
asdat <- stickSNPs
sample.meta(asdat)$phenotype <- rnorm(ncol(asdat))
sample.meta(asdat)$cat_phenotype <- sample(c("case", "control"), ncol(asdat), replace = TRUE)

#=========association=========
test_that("correct gmmat", {
  skip_if_not_installed("GMMAT"); skip_on_cran()
  suppressWarnings(asgmmat <- calc_association(asdat, response = "phenotype"))
  asgmmat <- get.snpR.stats(asgmmat, stats = "association")
  hits <- asgmmat$single[which(asgmmat$single$gmmat_pval_phenotype < 0.05),]
  expect_equal(nrow(hits), 6)
})

test_that("correct armitage", {
  asarm <- calc_association(asdat, response = "cat_phenotype", method = "armitage")
  
  #code to generate test values
  # sn <- format_snps(asdat, "sn", interpolate = FALSE)[,-c(1:2)]
  # cc <- sample.meta(asdat)$cat_phenotype
  # cc <- as.numeric(as.factor(cc)) - 1
  # comp <- numeric(nrow(sn))
  # for(i in 1:nrow(sn)){
  #   nas <- which(is.na(sn[i,]))
  #   if(length(nas) > 0){
  #     suppressWarnings(comp[i] <- CATT::CATT(cc[-nas], sn[i,-nas])$p.value)
  #   }
  #   else{
  #     suppressWarnings(comp[i] <- CAT::CATT(cc, sn[i,])$p.value)
  #   }
  # }
  
  # result from code above
  comp <- c(0.87, 0.9349, 0.588, 0.0726, 0.0472, 0.7809, 0.6378, 0.5988, 0.0663, 
            0.1699, 0.4793, 0.2265, 0.4427, 0.9179, 0.9837, 0.7344, 0.1564, 
            0.8965, 0.528, 0.3657, 0.6657, 0.8753, 0.8975, 0.4399, 0.451, 0.4917,
            0.462, 0.2441, 0.7089, 0.8012, 0.9099, 0.8276, 0.8006, 0.1908, 0.1235, 0.1959, 
            0.3356, 0.102, 0.7218, 0.9849, 0.5368, 0.875, 0.9574, 0.0861, 0.1171, 
            0.793, 0.5591, 0.5038, 0.7329, 0.6578, 0.9435, 0.5509, 0.474, 0.4744, 0.319, 
            0.4582, 0.3031, 0.9932, 0.4605, 0.0763, 0.2394, 0.6418, 0.5003, 0.005, 
            0.9793, 0.5433, 0.3782, 0.2965, 0.9152, 0.8601, 0.4323, 0.2172, 0.1314, 
            0.3454, 0.2489, 0.4176, 0.4191, 0.2704, 0.647, 0.2577, 0.4959, 0.9261, 
            0.5158, 0.2599, 0.4756, 0.4141, 0.749, 0.8423, 0.7095, 0.5256, 1, 0.7817, 
            0.4607, 0.6198, 0.8414, 0.516, 0.5338, 0.3384, 0.7428, 0.5957)

  asarm <- get.snpR.stats(asarm, stats = "association")
  expect_equal(round(asarm$single$p_armitage_cat_phenotype, 2), 
               round(comp,2)) # from CATT
})

test_that("correct odds", {
  asodds <- calc_association(asdat, response = "cat_phenotype", method = "odds_ratio")
  asodds <- get.snpR.stats(asodds, stats = "association")
  check <- c(0.022, -0.029, -0.13, 0.253, 0.363, 0.042, 0.107, 0.085, 0.324, -0.397, -0.095, 
             -0.17, -0.136, -0.019, -0.003, 0.047, -0.247, -0.019, -0.165, -0.282, 0.061, -0.025, -0.04, 
             0.228, -0.108, -0.093, 0.144, 0.171, -0.069, 0.055, -0.016, 0.03, 0.038, -0.21, -0.442, 0.244, 0.145, 0.237, 
             -0.074, 0.003, -0.093, -0.021, 0.01, -0.337, -Inf, -0.048, 0.09, -0.091, 0.135, 0.103, 
             0.011, -0.08, -0.118, -0.369, -0.141, -0.097, -0.376, -0.001, 0.097, -0.244, -0.426, 
             0.093, -0.105, -0.376, 0.004, -0.085, -0.296, -0.155, 0.024, -0.024, 0.174, -0.164,
             0.24, -0.155, -0.135, 0.159, -0.107, 0.253, -0.083, 0.141, 0.077, -0.021, -0.17, 0.171, 
             -0.095, -0.107, 0.048, -0.028, -0.059, 0.236, 0, -0.045, 0.123, -0.115, 0.039, -0.141, 0.231,
             -0.145, 0.048, 0.073) # hand checked
  expect_equal(round(asodds$single$log_odds_ratio_cat_phenotype, 3), check)
})


test_that("correct chisq", {
  aschi <- calc_association(asdat, response = "cat_phenotype", method = "chisq")
  aschi <- get.snpR.stats(aschi, stats = "association")
  check <- c(0.867, 0.936, 0.602, 0.094, 0.044, 0.793, 0.633, 0.599, 0.072, 0.182, 
             0.482, 0.31, 0.425, 0.915, 0.983, 0.751, 0.176, 0.889, 0.541, 0.377, 
             0.692, 0.871, 0.9, 0.412, 0.456, 0.529, 0.491, 0.225, 0.703, 0.811, 
             0.909, 0.853, 0.794, 0.172, 0.136, 0.209, 0.353, 0.105, 0.743, 0.985, 
             0.513, 0.876, 0.961, 0.065, 0.119, 0.763, 0.545, 0.512, 0.736, 0.672, 
             0.945, 0.568, 0.462, 0.478, 0.321, 0.537, 0.314, 0.993, 0.476, 0.087, 
             0.25, 0.661, 0.491, 0.006, 0.978, 0.538, 0.338, 0.271, 0.907, 0.875, 
             0.453, 0.209, 0.156, 0.399, 0.304, 0.331, 0.497, 0.229, 0.657, 0.321, 0.57, 
             0.929, 0.529, 0.213, 0.476, 0.507, 0.743, 0.84, 0.69, 0.532, 1, 0.792, 
             0.503, 0.614, 0.834, 0.537, 0.541, 0.34, 0.724, 0.605)
  expect_equal(round(aschi$single$chi_p_cat_phenotype, 3), check) # Hand checked, should not change.
})

#=========random forest========
test_that("random forest",{
  skip_if_not_installed("ranger")
  # basic
  expect_warning(rf <- run_random_forest(asdat, response = "phenotype"), 
                 regexp = "No p-values calcuated. When a quantitative")
  rfstats <- get.snpR.stats(rf$x, stats = "random_forest")
  expect_s3_class(rf$models$.base_.base$model, "ranger")
  expect_equal(dim(rf$models$.base_.base$predictions), c(nsamps(asdat), 2))
  expect_equal(length(rf$models), 1)
  expect_equal(unique(rfstats$single$subfacet), c(".base"))
  expect_equal(unique(rfstats$single$facet), c(".base"))
  expect_equal(colnames(rfstats$single), c("facet", "subfacet", "chr", "position", "phenotype_RF_importance"))
  
  # check importance
  suppressWarnings(rf2 <- run_random_forest(asdat, response = "cat_phenotype")) # will sometimes throw an inaccurate p-values warning from ranger
  imp <- get.snpR.stats(rf2$x, stats = "random_forest")$single
  expect_true(all(imp$cat_phenotype_RF_importance_pvals >= 0))
  expect_true(is.numeric(imp$cat_phenotype_RF_importance))

  # several facets
  rf <- run_random_forest(asdat[pop = c("ASP", "PAL")], facets = "pop", response = "phenotype", pvals = FALSE)
  rfstats <- get.snpR.stats(rf$x, "pop", "random_forest")
  expect_equal(length(rf$models), 2)
  expect_s3_class(rf$models$pop_ASP$model, "ranger")
  expect_s3_class(rf$models$pop_PAL$model, "ranger")
  expect_equal(dim(rf$models$pop_ASP$predictions), c(sum(sample.meta(asdat)$pop == "ASP"), 2))
  expect_equal(dim(rf$models$pop_PAL$predictions), c(sum(sample.meta(asdat)$pop == "PAL"), 2))
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
  skip_on_cran()
  
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
  gp <- run_genomic_prediction(asdat[pop = c("ASP", "PAL")], facets = "pop", response = "phenotype", iterations = 200, burn_in = 100, thin = 10)
  gpstats <- get.snpR.stats(gp$x, "pop", "genomic_prediction")
  expect_equal(length(gp$models), 2)
  expect_s3_class(gp$models$pop_ASP$model, "BGLR")
  expect_s3_class(gp$models$pop_PAL$model, "BGLR")
  expect_equal(dim(gp$models$pop_ASP$predictions), c(sum(sample.meta(asdat)$pop == "ASP"), 2))
  expect_equal(dim(gp$models$pop_PAL$predictions), c(sum(sample.meta(asdat)$pop == "PAL"), 2))
  expect_equal(unique(gpstats$single$subfacet), c("ASP", "PAL"))
  expect_equal(unique(gpstats$single$facet), c("pop"))
  expect_equal(colnames(gpstats$single), c("facet", "subfacet", "chr", "position", "phenotype_gp_effect"))
  str <- .paste.by.facet(gpstats$single, c("chr", "position"))
  expect_equal(as.numeric(table(str)), rep(2, nrow(asdat))) # each snp has calcs for each pop
  
  
  # runing on a previously run set
  gp <- run_genomic_prediction(asdat[pop = c("ASP", "PAL")], response = "phenotype", iterations = 200, burn_in = 100, thin = 10)
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

  set.seed(1212)
  dat <- stickSNPs
  sample.meta(dat)$phenotype <- rnorm(ncol(dat))
  sample.meta(dat)$cat_phenotype <- sample(c("case", "control"), ncol(dat), replace = TRUE)
  ## run cross_validation
  res <-cross_validate_genomic_prediction(dat, response = "phenotype", iterations = 200, burn_in = 100, thin = 10)
  
  expect_s3_class(res$model$models$.base_.base$model, "BGLR")
  expect_equal(length(res$model.samples), floor(nsamps(stickSNPs)*.9))
  expect_false(any(res$model.samples %in% res$cross.samples))
  expect_equal(colnames(res$comparison), c("observed", "predicted"))
  expect_equal(nrow(res$comparison), length(res$cross.samples))
  expect_equal(cor(res$comparison$observed, res$comparison$predicted)^2, res$rsq)
  
})
