context("IBD")

test_that("correct ibd", {
  skip_if_not_installed("ade4")
  skip_if_not_installed("geosphere")
  set.seed(1212)
  
  # run
  y <- stickSNPs[sample(nrow(stickSNPs), 100,  F), sample(ncol(stickSNPs), 100, F)]
  sample.meta(y) <- cbind(sample.meta(y), x = rnorm(ncol(y)), y = rnorm(ncol(y)))
  y <- calc_isolation_by_distance(y, facets = c(".base", "pop", "pop.group","pop.group.fam"))
  res <- get.snpR.stats(y, facets = c(".base", "pop", "pop.group","pop.group.fam"), stats = "ibd") # fetch result
  
  # check
  expect_equal(names(res), c("pop", "group.pop","fam.group.pop", ".base"))
  expect_true(all(unlist(lapply(unlist(unlist(res, recursive = F), recursive = F), function(x) "mantelrtest" %in% class(x))))) # checks that each result is a randtest
  
  # geo dists
  res <- get.snpR.stats(y, c(".base", "pop", "pop.group","pop.group.fam"), "geo_dist")
  expect_equal(names(res), c(".base", "pop", "fam.pop"))
  expect_equal(unlist(lapply(purrr::map(res, "Edwards"), class)), c(.base = "dist", pop = "dist", fam.pop = "dist"))
  expect_equivalent(res$.base$Edwards, stats::dist(sample.meta(y)[,c("x", "y")])) # check base
  geoloc <- matrix(NA, 6, 2)
  upops <- unique(sample.meta(y)$pop)
  for(i in 1:nrow(geoloc)){
    geoloc[i,] <- geosphere::geomean(sample.meta(y)[which(sample.meta(y)$pop == upops[i]), c("x", "y")])
  }
  rownames(geoloc) <- upops
  geoloc <- geoloc[sort(rownames(geoloc)),]
  expect_equivalent(res$pop$Edwards, stats::dist(geoloc)) # check pop
})

test_that("correct genetic distances", {
  set.seed(1212)
  
  # run
  y <- stickSNPs[sample(nrow(stickSNPs), 100,  F), sample(ncol(stickSNPs), 100, F)]
  res <- calc_genetic_distances(y, facets = c(".base", "pop", "pop.group","pop.group.fam"))
  res <- get.snpR.stats(res, facets = c(".base", "pop", "pop.group","pop.group.fam"), stat = "genetic_distance")
  
  check <- unlist(lapply(unlist(unlist(res, recursive = F), recursive = F), class))
  expect_true(all(check == "dist")) # checks that each result is a distance matrix
  expect_equal(length(check), 44) # each part there
  
})