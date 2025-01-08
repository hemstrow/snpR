tdat <- stickRAW[1:4, 1:6]
tdat <- import.snpR.data(tdat[,-c(1:2)], tdat[,1:2], data.frame(pop = rep("ASP", 4)))

correct_geno_tables <- list(gs = Matrix::Matrix(as.matrix(data.frame(AA = as.integer(c(0, 0, 1, 0)),
                                                                     CT =  as.integer(c(0, 1, 0, 0)),
                                                                     GG =  as.integer(c(4, 0, 1, 3)),
                                                                     TT =  as.integer(c(0, 3, 0, 0))))),
                            as = Matrix::Matrix(as.matrix(data.frame(A = as.integer(c(0, 0, 2, 0)),
                                                                     C =  as.integer(c(0, 1, 0, 0)),
                                                                     T =  as.integer(c(0, 7, 0, 0)),
                                                                     G =  as.integer(c(8, 0, 2, 6))))),
                            wm = Matrix::Matrix(as.matrix(data.frame(NN =  as.integer(c(0, 0, 2, 1)))), sparse = TRUE)
)

test_that("Overall snpRdata object", {
  expect_s4_class(tdat, "snpRdata")
  expect_is(tdat, "data.frame")
})

test_that("Geno table creation", {
  # vals and classes
  expect_equal(tdat@geno.tables, correct_geno_tables)
  expect_is(tdat@geno.tables, "list")
  expect_is(tdat@geno.tables$gs, "sparseMatrix")
  expect_is(tdat@geno.tables$as, "sparseMatrix")
  expect_is(tdat@geno.tables$wm, "sparseMatrix")
  
  # missing data check
  expect_true(tdat@mDat %in% colnames(tdat@geno.tables$wm))
  expect_true(!tdat@mDat %in% colnames(tdat@geno.tables$gs))
  
  # relevant test for some other data types, to implement fully later
  expect_true(all(Matrix::rowSums(tdat@geno.tables$gs)*2 == Matrix::rowSums(tdat@geno.tables$as)))
  
  # cases where genotypes are flipped
  test <- genotypes(.internal.data$test_snps)
  test <- as.matrix(test)
  test_a1 <- substr(test, 1, 1)
  test_a2 <- substr(test, 2, 2)
  test[,sample.meta(.internal.data$test_snps)$pop == "ASP"] <- 
    paste0(test_a2[,sample.meta(.internal.data$test_snps)$pop == "ASP"], test_a1[,sample.meta(.internal.data$test_snps)$pop == "ASP"])
  
  testd <- import.snpR.data(test, 
                            snp.meta(.internal.data$test_snps),
                            sample.meta(.internal.data$test_snps))
  testd <- .add.facets.snpR.data(testd, "pop")
  test2d <- .add.facets.snpR.data(.internal.data$test_snps, "pop")
  expect_identical(test2d@geno.tables, testd@geno.tables)
  
  # cases where alleles appear in different orders in specific facet levels
  test[,c(1,5)] <- test[,c(5, 1)]
  test[1,1] <- "GG"
  testd <- import.snpR.data(test, 
                            snp.meta(.internal.data$test_snps),
                            sample.meta(.internal.data$test_snps))
  testd <- .add.facets.snpR.data(testd, "pop")
  test2 <- genotypes(.internal.data$test_snps)
  test2[1,5] <- "GG"
  test2d <- import.snpR.data(test2, 
                             snp.meta(.internal.data$test_snps),
                             sample.meta(.internal.data$test_snps))
  test2d <- .add.facets.snpR.data(test2d, "pop")
  
  
  expect_identical(test2d@geno.tables, testd@geno.tables)
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
  
  
  # tibbles
  spn <- dplyr::as_tibble(sample.meta(stickSNPs))
  d <- import.snpR.data(genotypes(stickSNPs), snp.meta(stickSNPs), spn)
  expect_true(!methods::is(d@sample.meta, "tbl"))
})



