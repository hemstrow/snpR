test_that("bootstrapping",{
  skip_on_cran()
  
  eval_test <- function(bdr, fst = TRUE){
    expect_true(all(bdr$stat$bootstraps %in% c("pi", "fst")))
    expect_true(all(c("pi", "p_pi_overall_BY", "p_pi_byfacet_BY") %in% colnames(bdr$single.window)))
    if(fst){
      expect_true(all(c("fst", "p_fst_overall_BY", "p_fst_byfacet_BY") %in% colnames(bdr$pairwise.window)))
    }
    expect_true(all(!bdr$bootstraps$triple_sigma))
    expect_true(all(bdr$bootstraps$gaussain))
    expect_true(all(bdr$bootstraps$step == 200))
    expect_true(all(bdr$bootstraps$sigma == 100))
    expect_true(all(bdr$bootstraps$nk))
  }

  
  d <- calc_pairwise_fst(stickSNPs, c("pop", "pop.fam"))
  d <- calc_pi(d, c("pop", "pop.fam"))
  d <- calc_ho(d, c("pop", "pop.fam"))
  d <- calc_smoothed_averages(d, c("pop.chr","pop","pop.fam"), sigma = 100, 
                              triple_sigma = FALSE, 
                              gaussian = FALSE)
  d <- calc_smoothed_averages(d, c("pop.chr"), sigma = 100, stats.type = "single",
                              triple_sigma = FALSE, 
                              gaussian = FALSE)
  
  # check the basic approach
  bd <- do_bootstraps(d,
                      facets = c("pop.chr","pop","pop.fam"),
                      boots = 100,
                      sigma = 100,
                      triple_sigma = FALSE, 
                      gaussian = FALSE)
  bdr <- get.snpR.stats(bd, c("pop.chr","pop","pop.fam"), stats = c("fst", "pi"), bootstraps = TRUE)
  eval_test(bdr)
  
  
  
  
  # in par
  bd <- do_bootstraps(d,
                      facets = c("pop.chr","pop","pop.fam"),
                      boots = 100,
                      sigma = 100,
                      triple_sigma = FALSE, 
                      gaussian = FALSE,
                      par = 2)
  bdr <- get.snpR.stats(bd, c("pop.chr","pop","pop.fam"), stats = c("fst", "pi"), bootstraps = TRUE)
  eval_test(bdr)
  
  
  # single
  bd <- do_bootstraps(d,
                      facets = "chr.pop", 
                      statistics = c("pi", "ho"),
                      boots = 100,
                      sigma = 100,
                      triple_sigma = FALSE, 
                      gaussian = FALSE)
  bdr <- get.snpR.stats(bd, "pop.chr", stats = "pi", bootstraps = TRUE)
  eval_test(bdr, fst = FALSE)
  
  # in par
  bd <- do_bootstraps(d,
                      facets = "chr.pop", 
                      statistics = c("pi", "ho"),
                      boots = 100,
                      sigma = 100,
                      triple_sigma = FALSE, 
                      gaussian = FALSE,
                      par = 2)
  bdr <- get.snpR.stats(bd, "pop.chr", stats = "pi", bootstraps = TRUE)
  eval_test(bdr, fst = FALSE)
  
  bdr <- get.snpR.stats(bd, "pop.chr", stats = c("pi", "ho"), bootstraps = TRUE)
  expect_true("ho" %in% bdr$bootstraps$stat)
  expect_true(all(c("p_ho", "p_pi", "p_ho_overall_BY", "p_ho_byfacet_BY", "p_pi_overall_BY", "p_pi_byfacet_BY")
              %in% colnames(bdr$single.window)))
  
  
  # specific bug checking
  ## "_" in facet level
  test <- stickSNPs
  sample.meta(test)$test <- paste0("A_", sample.meta(test)$fam)
  test <- calc_pairwise_fst(test, "test")
  test <- calc_pi(test, "test")
  test <- calc_smoothed_averages(test, c("test", "test.chr"), sigma = 100)
  test <- do_bootstraps(test, c("test", "test.chr"), 100, 100, statistics = c("pi", "fst"))
  bdr <- get.snpR.stats(test, c("test", "test.chr"), c("pi", "fst"), bootstraps = TRUE)
  expect_true(all(c("A_A", "A_B", "A_C", "A_D") %in% bdr$bootstraps$subfacet)) # just basic tests to make sure the names come out OK
  expect_true(all(c("A_A~A_B", "A_A~A_C", "A_A~A_C", "A_A~A_D") %in% bdr$bootstraps$subfacet))
})
  