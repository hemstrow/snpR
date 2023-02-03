test_that("objects, NN", {
  check <- import.snpR.data(genotypes(stickSNPs), snp.meta(stickSNPs), sample.meta(stickSNPs))
  expect_true(is.snpRdata(check))
  expect_identical(genotypes(stickSNPs), genotypes(check))
  expect_identical(snp.meta(stickSNPs), snp.meta(check))
  expect_identical(sample.meta(stickSNPs), sample.meta(check))

  
  check <- import.snpR.data(genotypes(stickSNPs), snp.meta(stickSNPs))
  expect_true(is.snpRdata(check))
  expect_identical(genotypes(stickSNPs), genotypes(check))
  expect_identical(snp.meta(stickSNPs), snp.meta(check))
  expect_true(colnames(sample.meta(check))[1] == "sampID")
  
  check <- import.snpR.data(genotypes(stickSNPs), sample.meta = sample.meta(stickSNPs))
  expect_true(is.snpRdata(check))
  expect_identical(genotypes(stickSNPs), genotypes(check))
  expect_identical(sample.meta(stickSNPs), sample.meta(check))
  expect_true(colnames(snp.meta(check))[1] == "snpID")
  
})

test_that("files",{
  format_snps(stickSNPs, "NN", outfile = "test.txt")
  write.table(sample.meta(stickSNPs), "test_samp.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)
  write.table(snp.meta(stickSNPs), "test_snp.txt", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  check <- import.snpR.data("test.txt", "test_snp.txt", "test_samp.txt", header_cols = 3)
  expect_equivalent(stickSNPs, check)
  
  file.remove("test.txt", "test_snp.txt", "test_samp.txt")
})

test_that("vcf", {
  skip_if_not_installed("vcfR")
  format_snps(stickSNPs, output = "vcf", outfile = "test.vcf")
  
  expect_warning(.make_it_quiet(check <- import.snpR.data("test.vcf")), "sample metadata columns contain unacceptable special characters")
  expect_identical(check@geno.tables$as, stickSNPs@geno.tables$as)
  expect_equivalent(check@geno.tables$gs, stickSNPs@geno.tables$gs)
  expect_identical(colnames(check@geno.tables$gs), colnames(stickSNPs@geno.tables$gs))
  
  expect_warning(.make_it_quiet(check_b <- read_vcf("test.vcf")), "sample metadata columns contain unacceptable special characters")
  expect_identical(check, check_b)
  
  .make_it_quiet(check <- import.snpR.data("test.vcf", snp.meta(stickSNPs), sample.meta(stickSNPs)))
  expect_identical(snp.meta(stickSNPs), snp.meta(check))
  expect_identical(sample.meta(stickSNPs), sample.meta(check))
  
  file.remove("test.vcf")
})

test_that("structure",{
  format_snps(stickSNPs, output = "structure", outfile = "test.str")
  
  expect_warning(check <- read_structure("test.str", rows_per_individual = 2, header_cols = 1), "sample metadata columns contain")
  check <- calc_ho(check)
  y <- calc_ho(stickSNPs)
  expect_identical(get.snpR.stats(check, stats = "maf")$single$maf, get.snpR.stats(stickSNPs, stats = "maf")$single$maf) # tests for correct alleles
  expect_identical(get.snpR.stats(check, stats = "ho")$single$ho, get.snpR.stats(y, stats = "ho")$single$ho) # checks for correct genotypes
  rm(y)
  
  expect_warning(check <- read_structure("test.str", snp.meta = snp.meta(stickSNPs), sample.meta = sample.meta(stickSNPs),
                          rows_per_individual = 2, header_cols = 1), "Since allelic identities are not clear, alleles at each locus will be saved as A and C.")
  expect_identical(snp.meta(stickSNPs), snp.meta(check))
  expect_identical(sample.meta(stickSNPs), sample.meta(check))
  
  expect_identical(suppressWarnings(read_structure("test.str", snp.meta = snp.meta(stickSNPs), sample.meta = sample.meta(stickSNPs),
                                                   rows_per_individual = 2, header_cols = 1)),
                   suppressWarnings(import.snpR.data("test.str", snp.meta = snp.meta(stickSNPs), sample.meta = sample.meta(stickSNPs),
                                                     rows_per_individual = 2, header_cols = 1)))
  
  file.remove("test.str")
})

test_that("non-biallelic",{
  td <- read_non_biallelic(genotypes(steelMSATs), sample.meta = sample.meta(steelMSATs))
  expect_true(is.snpRdata(td))
  expect_false(.is.bi_allelic(td))
  expect_identical(td, steelMSATs)
})
