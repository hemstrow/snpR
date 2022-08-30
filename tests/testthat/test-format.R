test_that("vcf",{
  vcf <- format_snps(.internal.data$test_snps, "vcf")
  expect_equal(colnames(vcf$genotypes), 
               c("#CHROM",
                 "POS",
                 "ID",
                 "REF",
                 "ALT",
                 "QUAL",
                 "FILTER",
                 "INFO",
                 "FORMAT",
                 colnames(genotypes(.internal.data$test_snps))))
  
  expect_equal(which(genotypes(.internal.data$test_snps) == "NN"), which(vcf$genotypes[,-c(1:9)] == "./."))
})


test_that("genalex",{
  # format
  format_snps(.internal.data$test_snps, "genalex", facets = "pop", outfile = "test.xlsx")
  
  # read and re-convert
  check <- openxlsx::read.xlsx("test.xlsx", colNames = FALSE)
  pops <- check$X2[-c(1:3)]
  expect_equivalent(as.numeric(check[1,1:5]), c(nsnps(.internal.data$test_snps), nsamps(.internal.data$test_snps), 2, 5, 5))
  check <- check[4:nrow(check),3:ncol(check)]
  check <- t(check)
  check[check == 1] <- "A"
  check[check == 2] <- "C"
  check[check == 3] <- "G"
  check[check == 4] <- "T"
  check[check == 0] <- "N"
  check <- paste0(check[seq(1, nrow(check), 2),], check[seq(2, nrow(check), 2),])
  check <- matrix(check, 11, 10)
  check <- import.snpR.data(as.data.frame(check), snp.meta(.internal.data$test_snps), sample.meta = data.frame(pop = pops))
  
  
  # ho, then eval identical
  check <- calc_ho(check, "pop")
  check2 <- calc_ho(.internal.data$test_snps, "pop")
  expect_identical(get.snpR.stats(check, "pop", "ho")$weighted.means,
                   get.snpR.stats(check2, "pop", "ho")$weighted.means)

  file.remove("test.xlsx")

})
