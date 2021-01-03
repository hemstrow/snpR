context("multiple testing correction")

p <- c(0.0001, 0.0010, 0.0062, 0.0101, 0.0214, 0.0227, 0.0273, 0.0292, 0.0311, 0.0323, 0.0441, 0.0490, 0.0573, 0.1262, 0.5794)

check <- fwe_correction(p)

test_that("bonferroni",{
  expect_equal(check$p_overall_bonferroni <= .05, c(T, T, rep(F, length(p) - 2)))
})

test_that("holm",{
  expect_equal(check$p_overall_holm <= .05, c(T, T, rep(F, length(p) - 2)))
})

test_that("BH",{
  expect_equal(check$p_overall_BH <= .05, c(rep(T, length = 10), rep(F, length(p) - 10)))
})

test_that("BY",{
  expect_equal(check$p_overall_BY <= .05, c(rep(T, length = 2), rep(F, length(p) - 2)))
})