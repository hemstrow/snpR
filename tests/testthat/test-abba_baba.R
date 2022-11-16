test_that("correct abba_baba", {
  x <- stickSNPs
  maf <- get.snpR.stats(x, ".base", stats = "maf")$single
  snp.meta(x)$ref <- maf$major
  snp.meta(x)$anc <- maf$minor

  x <- calc_abba_baba(x, "pop.chr", "ASP", "UPD", "PAL", TRUE, sigma = 1000)
  r1 <- get.snpR.stats(x, "pop.chr", "abba_baba") # gets the per chr results
  expect_equal(unique(r1$pairwise$facet), "pop")
  expect_equal(unique(r1$pairwise$comparison), "ASP~UPD~PAL")
  expect_equal(round(r1$pairwise$D_abba_baba, 3),
               c(0, 1, 0.167, -0.646, -0.18, -0.25, 1, -0.062, 
                 -0.112, -0.25, 0.117, -0.125, -0.185, -0.397, 
                 -0.489, 0.118, -0.015, -0.282, -1, -0.133, -0.429, 
                 0.18, 0.208, -0.39, -0.282, 0.097, -0.187, -0.474, 
                 -0.378, -0.333, -0.342, -0.352, -0.171, 0.602, 0.182, 
                 0.358, -0.191, 0.421, 0.062, 0.083, NaN, -0.25, -0.5, 
                 -0.136, -1, 0.19, 0.19, -0.238, -1, 0.393, -0.008, 0.014,
                 -0.481, NaN, 0.077, 0.243, 1, 0.362, 0.143, 0.226, -1,
                 -0.038, 0.239, -0.107, -0.489, 0.154, -1, -0.25, NaN, 
                 0.133, -0.2, -0.188, 0.125, 0.15, -0.302, NaN, 0.294, 
                 NaN, 0.217, -0.112, -0.128, -0.364, 0.714, NaN, 0.417, 
                 0.239, 0.143, NaN, -0.474, 0.309, -0.471, 0.6, -0.357, 
                 -1, 0.111, -0.111, -1, -0.299, 0.1, -0.015))
  expect_true(all(is.na(c(r1$weighted.means$D_abba_baba_jackknife & 
                            r1$weighted.means$D_abba_baba_p_value))))
  expect_equal(round(r1$weighted.means$D_abba_baba, 3),
               c(0.11, -0.212, -0.049, 0.136, -0.009, -0.25, 0.036, -0.157, 
                 -0.043, -0.553, 0.368, -0.302, -0.287, 0.239, 0.169, -0.012, 
                 -0.213, -0.01, 0.362, -0.227, 0.062))
  r2 <- get.snpR.stats(x, "pop", "abba_baba") # gets the overall results
  expect_equal(round(as.numeric(unlist(r2$weighted.means[1, 5:6]),3)),
               round(c(-0.0640, 0.060)))
  expect_equal(nrow(r2$weighted.means), 1)
               

  # smoothed windowed averages
  x <- calc_smoothed_averages(x, "pop.chr", sigma = 200, step = 200,
     nk = TRUE, stats.type = "pairwise")
  expect_equal(nrow(get.snpR.stats(x, "pop.chr", "abba_baba")$pairwise.window),
               487)
})
  