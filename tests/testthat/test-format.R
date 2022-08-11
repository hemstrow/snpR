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



