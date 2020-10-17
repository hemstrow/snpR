context("snpRdata object creation/tabulation")

tdat <- stickRAW[1:4, 1:7]
tdat <- import.snpR.data(tdat[,-c(1:3)], tdat[,1:3], data.frame(pop = rep("ASP", 4)))

correct_geno_tables <- list(gs = as.matrix(data.frame(CC = as.integer(c(3, 0, 0, 0)),
                                                      CT =  as.integer(c(0, 1, 0, 0)),
                                                      GG =  as.integer(c(0, 0, 4, 0)),
                                                      TT =  as.integer(c(0, 2, 0, 4)))),
                            as = as.matrix(data.frame(C =  as.integer(c(6, 1, 0, 0)),
                                                      T =  as.integer(c(0, 5, 0, 8)),
                                                      G =  as.integer(c(0, 0, 8, 0)))),
                            wm = as.matrix(data.frame(CC =  as.integer(c(3, 0, 0, 0)),
                                                      CT =  as.integer(c(0, 1, 0, 0)),
                                                      GG =  as.integer(c(0, 0, 4, 0)),
                                                      NN =  as.integer(c(1, 1, 0, 0)),
                                                      TT =  as.integer(c(0, 2, 0, 4))))
)

test_that("Overall snpRdata object", {
  expect_s4_class(tdat, "snpRdata")
  expect_is(tdat, "data.frame")
})

test_that("Geno table creation", {
  # vals and classes
  expect_equal(tdat@geno.tables, correct_geno_tables)
  expect_is(tdat@geno.tables, "list")
  expect_is(tdat@geno.tables$gs, "matrix")
  expect_is(tdat@geno.tables$as, "matrix")
  expect_is(tdat@geno.tables$wm, "matrix")
  
  # missing data check
  expect_true(tdat@mDat %in% colnames(tdat@geno.tables$wm))
  expect_true(!tdat@mDat %in% colnames(tdat@geno.tables$gs))
})

test_that("AC table creation",{
  expect_is(tdat@ac, "data.frame")
  expect_equal(tdat@ac$ni1, matrixStats::rowMaxs(tdat@geno.tables$as))
  expect_equal(tdat@ac$n_total, rowSums(tdat@geno.tables$as)) # check minor by proxy
})

test_that("Stats table creation,"{
  expect_is(tdat@stats, "data.table")
  expect_true(all(c("maf", "major", "minor", "maj.count", "min.count") %in% colnames(tdat@stats))) # was maf added?
})




