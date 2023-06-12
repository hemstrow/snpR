test_that("NeEstimator",{
  local_edition(3)
  skip_on_cran(); skip_on_ci()
  
  ne_path <- "C://usr/bin/Ne2-1.exe"
  skip_if(!file.exists(ne_path))
  
  
  ne <- calc_ne(stickSNPs[pop = "ASP"], NeEstimator_path = ne_path, chr = "chr", facets = "pop")
  ne <- get.snpR.stats(ne, "pop", "ne")
  
  
  expect_equal(as.numeric(unlist(ne$pop[1,-c(1:2), drop = TRUE])), c(108.3, 108.3, 83.5, 27.1, 27.1, 24.4, Inf, Inf, Inf, 20.9, 20.9, 18.3, Inf, Inf, Inf)) # not internally calced, just a check for proper prep and parsing
  
  ne <- calc_ne(stickSNPs[pop = "ASP"], pcrit = 0, NeEstimator_path = ne_path, chr = "chr", facets = "pop")
  ne <- get.snpR.stats(ne, "pop", "ne")
  expect_equal(as.numeric(unlist(ne$pop[1,-c(1:2), drop = TRUE])), c(108.3, 27.1, Inf, 20.9, Inf)) # tests two bugs, one that occurs when you have one pcrit, one that occurs when you have one pcrit AND one pop!
  
  
  ne <- calc_ne(stickSNPs[pop = "ASP"], pcrit = 0, NeEstimator_path = ne_path, chr = "chr", facets = "pop", methods = c("het", "coan", "LD"))
  ne <- get.snpR.stats(ne, "pop", "ne")
  expect_equal(as.numeric(unlist(ne$pop[1,-c(1:2), drop = TRUE])), c(10.8, 3.2, 22.8, 9.2, 5, 127.9, 108.3, 27.1, Inf, 20.9, Inf)) # tests multiple methods at once
  
  set.seed(12342)
  # subsampling -- note that this isn't a direct test, but the numbers below are what should be created if this ran correctly with 50 SNPs.
  ne2 <- calc_ne(stickSNPs[pop = "ASP"], pcrit = 0, NeEstimator_path = ne_path, chr = "chr", facets = "pop", methods = c("het", "coan", "LD"), nsnps = 50)
  ne2 <- get.snpR.stats(ne2, "pop", "ne")
  expect_equal(as.numeric(unlist(ne2$pop[1,-c(1:2), drop = TRUE])),
               c(4.6, 2.2, 8.0, 16.1, 4.9, Inf, 39.8, 10.6, Inf, 7.0, Inf))
  
  
  # temporal
  set.seed(12342)
  # errors
  expect_error(n3 <- calc_ne(stickSNPs, "pop", "chr", methods = c("LD", "Het", "Coan", "temporal"), pcrit = c(0, .02, 0.07),
                temporal_details = data.frame(t1 = c("ASP", "PAL", "OPL"),
                                              t2 = c("UPD", "CLF", "SMR"),
                                              gens = c(30, 30, 30))),
               "The temporal method cannot currently be run alongside other methods.")
  
  expect_error(n3 <- calc_ne(stickSNPs, c("pop", "fam"), "chr", methods = "temporal", pcrit = c(0, .02, 0.07),
                             temporal_details = data.frame(t1 = c("ASP", "PAL", "OPL"),
                                                           t2 = c("UPD", "CLF", "SMR"),
                                                           gens = c(30, 30, 30))),
               "Only one facet.")
  expect_error(n3 <- calc_ne(stickSNPs, "pop.fam", "chr", methods = "temporal", pcrit = c(0, .02, 0.07),
                temporal_details = data.frame(t1 = c("ASP.A"),
                                              t2 = c("UPD.A"),
                                              gens = c(30))),
               "alphabetical order")
  
  # basic
  n3 <- calc_ne(stickSNPs, "fam.pop", "chr", methods = "temporal", pcrit = c(0, .02, 0.07),
                temporal_details = data.frame(t1 = c("A.ASP"),
                                              t2 = c("A.UPD"),
                                              gens = c(30)))
  n3 <- get.snpR.stats(n3, "pop.fam", "ne")
  expect_true(ncol(n3$pop) == 47)
  expect_true(n3$pop$pop == "A.ASP~A.UPD")
  
  # check that plan I works without error (the test is actually in what is printed...)
  n3 <- calc_ne(stickSNPs, "fam.pop", "chr", methods = "temporal", pcrit = c(0, .02, 0.07),
                temporal_details = data.frame(t1 = c("A.ASP"),
                                              t2 = c("A.UPD"),
                                              gens = c(30),
                                              N = 10))
})



test_that("colony",{
  local_edition(3)
  skip_on_cran(); skip_on_ci()
  
  
  col_path <- "C://usr/bin/Colony/colony2s.exe"
  skip_if(!file.exists(col_path))
  
  .make_it_quiet(col <- run_colony(stickSNPs[pop = "ASP"], colony_path = col_path, run_length = 1, method = "PLS", cleanup = TRUE, verbose = FALSE))
  
  expect_identical(colnames(col$clusters), c("ClusterIndex", "ClusterProbability", "OffspringID", "FatherID", "MotherID"))
  snap_check <- col$dyads[col$dyads$Probability > .05, -3]
  rownames(snap_check) <- 1:nrow(snap_check)
  expect_snapshot_value(snap_check, style = "serialize")
  expect_true(is.snpRdata(col$x))
  expect_snapshot_value(col$clusters[col$clusters$ClusterProbability > .5, 3], style = "serialize") # not internally calced, just a check for proper prep and parsing
})


