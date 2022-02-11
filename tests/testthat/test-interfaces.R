context("Interfaces")


test_that("NeEstimator",{
  local_edition(3)
  skip_on_cran()
  
  ne_path <- "C://usr/bin/Ne2-1.exe"
  skip_if(!file.exists(ne_path))
  
  ne <- calc_ne(stickSNPs[pop = "ASP"], NeEstimator_path = ne_path, chr = "chr", facets = "pop")
  ne <- get.snpR.stats(ne, "pop", "ne")
  unlink("NeEstimator", recursive = T)
  
  expect_snapshot_value(ne, style = "serialize") # not internally calced, just a check for proper prep and parsing
})



test_that("colony",{
  local_edition(3)
  skip_on_cran()
  
  
  col_path <- "C://usr/bin/Colony/colony2s.exe"
  skip_if(!file.exists(col_path))
  
  col <- run_colony(stickSNPs[pop = "ASP"], colony_path = col_path, run_length = 1, method = "PLS", cleanup = TRUE)
  
  expect_identical(colnames(col$clusters), c("ClusterIndex", "ClusterProbability", "OffspringID", "FatherID", "MotherID"))
  snap_check <- col$dyads[col$dyads$Probability > .5, -3]
  rownames(snap_check) <- 1:nrow(snap_check)
  expect_snapshot_value(snap_check, style = "serialize")
  expect_true(is.snpRdata(col$x))
  expect_snapshot_value(col$clusters[col$clusters$ClusterProbability > .8, 3], style = "serialize") # not internally calced, just a check for proper prep and parsing
})


