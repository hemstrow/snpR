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
                            wm = as.matrix(data.frame(NN =  as.integer(c(0, 0, 2, 1))))
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
  
  # relevant test for some other data types, to implement fully later
  expect_true(all(rowSums(tdat@geno.tables$gs)*2 == rowSums(tdat@geno.tables$as)))
})


test_that("Stats table creation", {
  expect_is(tdat@stats, "data.table")
  expect_true(all(c("maf", "major", "minor", "maj.count", "min.count") %in% colnames(tdat@stats))) # was maf added?
})


test_that("facet error reporting",{
  test <- .internal.data$test_snps
  
  # sample meta
  expect_warning(sample.meta(test) <- cbind(sample.meta(test),
                                            test1 = sample(c("A", "B", ".C"), nsamps(test), replace = TRUE),
                                            test2 = sample(c("A", "B", "~C"), nsamps(test), replace = TRUE),
                                            test3 = sample(c("A", "B", " C"), nsamps(test), replace = TRUE),
                                            test4 = sample(c(".A", "~B", " C"), nsamps(test), replace = TRUE)),
                 "Unaccepted characters")
  
  expect_error(calc_ho(test, "test1"), "'\\.'")
  expect_error(calc_ho(test, "test2"), "'~'")
  expect_error(calc_ho(test, "test3"), "\\(spaces\\)")
  expect_error(calc_ho(test, "test4"), "\n\t'\\.'\n\t '~'\n\t ' ' \\(spaces\\)")
  
  
  # snp.meta
  test <- .internal.data$test_snps
  expect_warning(snp.meta(test) <- cbind(snp.meta(test),
                                         test1 = sample(c("A", "B", ".C"), nsnps(test), replace = TRUE),
                                         test2 = sample(c("A", "B", "~C"), nsnps(test), replace = TRUE),
                                         test3 = sample(c("A", "B", " C"), nsnps(test), replace = TRUE),
                                         test4 = sample(c(".A", "~B", " C"), nsnps(test), replace = TRUE)),
                 "Unaccepted characters")
  
  expect_error(calc_ho(test, "test1"), "'\\.'")
  expect_error(calc_ho(test, "test2"), "'~'")
  expect_error(calc_ho(test, "test3"), "\\(spaces\\)")
  expect_error(calc_ho(test, "test4"), "\n\t'\\.'\n\t '~'\n\t ' ' \\(spaces\\)")
  
  # dup facets
  expect_error(.suppress_specific_warning(snp.meta(test) <- cbind(snp.meta(test), pop = "test"), "unexpected"), "Bad column names: pop, pop\\.")
  
  # bad facet names
  expect_error(.suppress_specific_warning(snp.meta(test) <- cbind(snp.meta(test), `~pop` = "test"), "unexpected"), "Bad column names: ~pop\\.")
  expect_error(.suppress_specific_warning(snp.meta(test) <- cbind(snp.meta(test), `.pop` = "test"), "unexpected"), "Bad column names: \\.pop\\.")
  expect_error(.suppress_specific_warning(snp.meta(test) <- cbind(snp.meta(test), ` pop` = "test"), "unexpected"), "Bad column names: \\ pop\\.")
  expect_error(.suppress_specific_warning(snp.meta(test) <- cbind(snp.meta(test), all = "test"), "unexpected"), "Bad column names: all\\.")
  
  
  
  # facet contents detailed reporting
  meta <- sample.meta(stickSNPs)
  meta$test <- "."
  
  snpm <- snp.meta(stickSNPs)
  snpm$test2 <- "."
  
  ## bad characters
  expect_warning(import.snpR.data(genotypes(stickSNPs), sample.meta = sample.meta(stickSNPs), snp.meta = snpm), "ome snp metadata columns contain.+Facet: test2\tlevels: \\.")
  expect_warning(import.snpR.data(genotypes(stickSNPs), sample.meta = meta, snp.meta = snp.meta(stickSNPs)), "ome sample metadata columns contain.+Facet: test\tlevels: \\.")
  
  
  meta$dup_test <- "ASP"
  meta$test <- NULL
  snpm$dup_test2 <- "ASP"
  snpm$test2 <- NULL
  
  # dupliates
  expect_warning(import.snpR.data(genotypes(stickSNPs), sample.meta = meta, snp.meta = snpm), "ome levels are duplicated.+Level: ASP\tin facets: pop, dup_test, dup_test2")
  
})



