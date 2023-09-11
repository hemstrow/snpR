test_that("vcf",{
  vcf <- format_snps(.internal.data$test_snps, "vcf")
  expect_equal(c(colnames(vcf$data_meta), colnames(vcf$genotypes)), 
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

test_that("plink",{
  skip_if_not_installed("genio")
  
  dat <- stickSNPs
  
  matches <- dat@snp.meta$chr %in% paste0("group", c("IV", "XV", "VII"))
  snp.meta(dat)$chr[matches] <- paste0(c(50, 50, 30), dat@snp.meta$chr[matches])
  
  expect_warning(format_snps(dat, "plink", chr = "chr", outfile = "plink_test", plink_recode_numeric = FALSE), "numbers at the start of chr/scaffold names")
  
  
  expect_warning(check <- read_plink("plink_test"), " metadata columns contain unacceptable special")
  
  dat <- calc_pi(dat, "chr.fam")
  check <- calc_pi(check, "chr.fam")
  c1 <- get.snpR.stats(dat, "chr.fam", "pi")$weighted.means
  c2 <- get.snpR.stats(check, "chr.fam", "pi")$weighted.means
  expect_true(all(c1$weighted_mean_pi == c2$weighted_mean_pi)) # everything is good if this is true
  cleanups <- list.files(pattern = "plink_test")
  file.remove(cleanups)
  
  
  # check bad loci removal
  ## make bad data
  dat <- stickSNPs
  genos <- genotypes(dat)
  genos <- rbind(genos, rep("NN", ncol(genos)))
  sm <- snp.meta(dat)
  sm <- rbind(sm[,-ncol(sm)], c("test", 10))
  
  ## test
  dat <- import.snpR.data(genos, sm, sample.meta(dat))
  expect_warning(check <- format_snps(dat, "plink", "pop", outfile = "plink_test"), "loci without any called genotypes.")
  expect_warning(check <- read_plink("plink_test"), " metadata columns contain unacceptable special")
  dat <- calc_pi(stickSNPs, "chr.fam")
  check <- calc_pi(check, "chr.fam")
  c1 <- get.snpR.stats(dat, "chr.fam", "pi")$weighted.means
  c2 <- get.snpR.stats(check, "chr.fam", "pi")$weighted.means
  expect_true(all(c1$weighted_mean_pi == c2$weighted_mean_pi)) # everything is good if this is true
  cleanups <- list.files(pattern = "plink_test")
  file.remove(cleanups)
})
