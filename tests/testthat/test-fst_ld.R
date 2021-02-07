context("fst and ld")
local_edition(3)

tdat <- stickRAW[1:100, 1:103]
tdat <- import.snpR.data(tdat[,-c(1:3)], tdat[,1:3], data.frame(pop = rep(c("A", "B"), each = 50)))
tdat <- add.facets.snpR.data(tdat, "pop")

test_that("correct genepop", {
  tdfst <- calc_pairwise_fst(tdat, "pop", "genepop")
  expect_snapshot(list(tdfst[[1]]@pairwise.stats, tdfst[[2]]))
})

test_that("correct wc", {
  tdfst <- calc_pairwise_fst(tdat, "pop", "wc")
  expect_snapshot(tdfst@pairwise.stats)
})

test_that("correct hoh", {
  tdfst <- calc_pairwise_fst(tdat, "pop", "hohenlohe")
  expect_snapshot(tdfst@pairwise.stats)
})

test_that("correct cld ld",{
  tdld <- calc_pairwise_ld(tdat, "pop")
  
})