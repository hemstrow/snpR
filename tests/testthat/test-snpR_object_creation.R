tdat <- stickRAW[1:4, 1:6]
tdat <- import.snpR.data(tdat[,-c(1:2)], tdat[,1:2], data.frame(pop = rep("ASP", 4)))

correct_geno_tables <- list(gs = as.matrix(data.frame(AA = as.integer(c(0, 0, 1, 0)),
                                                      CT =  as.integer(c(0, 1, 0, 0)),
                                                      GG =  as.integer(c(4, 0, 1, 3)),
                                                      TT =  as.integer(c(0, 3, 0, 0)))),
                            as = as.matrix(data.frame(A = as.integer(c(0, 0, 2, 0)),
                                                      C =  as.integer(c(0, 1, 0, 0)),
                                                      T =  as.integer(c(0, 7, 0, 0)),
                                                      G =  as.integer(c(8, 0, 2, 6)))),
                            wm = as.matrix(data.frame(AA =  as.integer(c(0, 0, 1, 0)),
                                                      CT =  as.integer(c(0, 1, 0, 0)),
                                                      GG =  as.integer(c(4, 0, 1, 3)),
                                                      NN =  as.integer(c(0, 0, 2, 1)),
                                                      TT =  as.integer(c(0, 3, 0, 0))))
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

test_that("Stats table creation", {
  expect_is(tdat@stats, "data.table")
  expect_true(all(c("maf", "major", "minor", "maj.count", "min.count") %in% colnames(tdat@stats))) # was maf added?
})




