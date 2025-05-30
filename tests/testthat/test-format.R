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
  
  # genotypes
  format_snps(stickSNPs, "plink", outfile = "plink_test")
  expect_warning(check <- read_plink("plink_test"), " metadata columns contain unacceptable special")
  expect_true(all(dplyr::arrange(stickSNPs@stats, chr, position)$major == dplyr::arrange(check@stats, chr, position)$major))
  
  # ped, fam, bim
  ## ped
  ped <- data.table::fread("plink_test.ped")
  A <- as.matrix(ped[,7:ncol(ped)])
  A1 <- substr(A, 1, 1)
  A2 <- substr(A, 3, 3)
  A <- paste0(A1, A2)
  A <- matrix(A, nsamps(stickSNPs), nsnps(stickSNPs))
  A <- t(A)
  A[A == "00"] <- "NN"
  G <- cbind(snp.meta(stickSNPs), genotypes(stickSNPs))
  G <- dplyr::arrange(G, chr, position)
  G <- G[,-c(1:ncol(snp.meta(stickSNPs)))]
  expect_true(all(A == as.matrix(G)))
  
  ## fam
  fam <- data.table::fread("plink_test.fam")
  expect_true(all(fam$V2 == colnames(stickSNPs)))
  
  ## bim
  bim <- data.table::fread("plink_test.bim")
  minors <- stickSNPs@stats
  minors <- dplyr::arrange(minors, chr, position)
  expect_true(all(bim$V5 == minors$minor))
  expect_true(all(bim$V6 == minors$major))
  
  # with odd names
  dat <- stickSNPs
  matches <- dat@snp.meta$chr %in% paste0("group", c("IV", "XV", "VII"))
  snp.meta(dat)$chr[matches] <- paste0(c(50, 50, 30), dat@snp.meta$chr[matches])
  
  expect_warning(format_snps(dat, "plink", chr = "chr", outfile = "plink_test", plink_recode_numeric = FALSE), "numbers at the start of chr/scaffold names")
  
  
  expect_warning(check <- read_plink("plink_test"), " metadata columns contain unacceptable special")
  
  
  
  dat <- calc_pi(dat, "chr.fam")
  check <- calc_pi(check, "chr.fam")
  c1 <- get.snpR.stats(dat, "chr.fam", "pi")$weighted.means
  c2 <- get.snpR.stats(check, "chr.fam", "pi")$weighted.means
  
  
  expect_true(all(round(c1$weighted_mean_pi, 5) == round(c2$weighted_mean_pi, 5))) # everything is good if this is true
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
  expect_true(all(round(c1$weighted_mean_pi, 5) == round(c2$weighted_mean_pi, 5))) # everything is good if this is true
  cleanups <- list.files(pattern = "plink_test")
  file.remove(cleanups)
})

test_that("sn",{
  # cases where one locus is completely unsequenced at one pop
  test <- genotypes(.internal.data$test_snps)
  test[1,sample.meta(.internal.data$test_snps)$pop == "ASP"] <- "NN"
  testd <- import.snpR.data(test, 
                            snp.meta(.internal.data$test_snps),
                            sample.meta(.internal.data$test_snps))
  res <- format_snps(testd, "ac", "pop")
  expect_true(all(res[1,c("n_total", "n_alleles", "ni1", "ni2")] == 0))
  
  
  # cases where one locus is completely unsequenced at all pops
  test <- genotypes(.internal.data$test_snps)
  test[1,] <- "NN"
  testd <- import.snpR.data(test, 
                            snp.meta(.internal.data$test_snps),
                            sample.meta(.internal.data$test_snps))
  res <- format_snps(testd, "ac", "pop")
  expect_true(all(res[1:2,c("n_total", "n_alleles", "ni1", "ni2")] == 0))
})
