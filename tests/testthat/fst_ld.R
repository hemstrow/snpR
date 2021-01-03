context("fst and ld")

tdat <- stickRAW[1:100, 1:103]
tdat <- import.snpR.data(tdat[,-c(1:3)], tdat[,1:3], data.frame(pop = rep(c("A", "B"), each = 50)))
tdat <- add.facets.snpR.data(tdat, "pop")

test_that("correct genepop", {
  tdfst <- calc_pairwise_fst(tdat, "pop", "genepop")
  expect_known_output(list(tdfst[[1]]@pairwise.stats, tdfst[[2]]), 
                      "tests/testthat/fst_genepop.test", print = TRUE)
})

test_that("correct wc", {
  tdfst <- calc_pairwise_fst(tdat, "pop", "wc")
  expect_known_output(tdfst@pairwise.stats, 
                      "tests/testthat/fst_wc.test", print = TRUE)
})

test_that("correct hoh", {
  tdfst <- calc_pairwise_fst(tdat, "pop", "hohenlohe")
  expect_known_output(tdfst@pairwise.stats, 
                      "tests/testthat/fst_hoh.test", print = TRUE)
})

test_that("correct cld ld",{
  tdld <- calc_pairwise_ld(tdat, "pop")
  
})