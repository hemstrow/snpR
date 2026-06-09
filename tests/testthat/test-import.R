test_that("objects, NN", {
  # basic
  check <- import.snpR.data(genotypes(stickSNPs), snp.meta(stickSNPs), sample.meta(stickSNPs))
  expect_true(is.snpRdata(check))
  expect_identical(genotypes(stickSNPs), genotypes(check))
  expect_identical(snp.meta(stickSNPs), snp.meta(check))
  expect_identical(sample.meta(stickSNPs), sample.meta(check))

  # samp ID added automatically
  check <- import.snpR.data(genotypes(stickSNPs), snp.meta(stickSNPs))
  expect_true(is.snpRdata(check))
  expect_identical(genotypes(stickSNPs), genotypes(check))
  expect_identical(snp.meta(stickSNPs), snp.meta(check))
  expect_true(colnames(sample.meta(check))[1] == "sampID")
  
  # SNP ID added automatically
  check <- import.snpR.data(genotypes(stickSNPs), sample.meta = sample.meta(stickSNPs))
  expect_true(is.snpRdata(check))
  expect_identical(genotypes(stickSNPs), genotypes(check))
  expect_identical(sample.meta(stickSNPs), sample.meta(check))
  expect_true(colnames(snp.meta(check))[1] == "snpID")
  
  # non ACGT first genos
  genos <- genotypes(stickSNPs)
  a1 <- sample(c("I", "D"), ncol(genos), TRUE)
  a2 <- sample(c("I", "D"), ncol(genos), TRUE)
  genos <- rbind(paste0(a1, a2), genos)
  sm <- snp.meta(stickSNPs)
  sm <- rbind(c("test", 1), sm[,-3])
  test <- import.snpR.data(genos, sm, sample.meta(stickSNPs), mDat = "NN")
  expect_true(is.snpRdata(test))
  ## back conversion to and from 0000
  .make_it_quiet(test0000 <- format_snps(test, "0000"))
  .make_it_quiet(testNN <- import.snpR.data(test0000, snp.meta(test), sample.meta(test), 
                                            header_cols = 3, mDat = "0000"))
  expect_true(all((genotypes(test) == genotypes(testNN))[-1,])) # ATCG should convert cleanly back and forth
  expect_true(all(genotypes(testNN)[1,] %in% c("BB", "BD", "DD", "DB"))) # the insert will not; given first two blank slots in LETTERS
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
  

  
  # non bi-allelic
  ## prep test data
  nb <- fread("test.vcf")
  set.seed(1231)
  tobenb <- sort(sample(nrow(nb), 20))
  nb$ALT[tobenb] <- paste0(nb$ALT[tobenb], ",NB")
  AS <- as.character(unlist(nb[tobenb,10:ncol(nb)]))
  nal <- sum(AS != "./.")
  reps <- which(rbinom(nal, 1, .2) == 1)
  AS[AS != "./."][reps] <- sample(c("0/2", "2/2", "1/2"), length(reps), TRUE)
  ocn <- colnames(nb)[10:ncol(nb)]
  ngenos <- matrix("./.", nrow(nb), ncol = ncol(nb) - 9)
  ngenos[tobenb,] <- matrix(AS, ncol = ncol(nb) - 9, nrow = length(tobenb))
  ngenos[-tobenb,] <- as.matrix(nb[-tobenb,10:ncol(nb)])
  
  nb <- cbind(nb[,1:9], ngenos)
  colnames(nb)[10:ncol(nb)] <- ocn
  header <- readLines("test.vcf", n = 28)
  writeLines(header, "test_nb.vcf")
  fwrite(nb, "test_nb.vcf", append = TRUE, quote = FALSE, col.names = TRUE, row.names = FALSE,
         sep = "\t")
  
  # read test data
  expect_warning(.make_it_quiet(t2 <- read_vcf("test_nb.vcf", snp.meta(stickSNPs), sample.meta(stickSNPs))), "Detected 20 non-biallelic loci")
  expect_equal(nrow(t2), nrow(nb) - 20)
  expect_true(.is.bi_allelic(t2))
  
  .make_it_quiet(t2 <- read_vcf("test_nb.vcf", snp.meta(stickSNPs), sample.meta(stickSNPs), accept_non_biallelic = TRUE))
  expect_true(!.is.bi_allelic(t2))
  expect_equal(colnames(t2@geno.tables$as),c("A01", "A02", "A03"))
  expect_equal(nrow(t2), nrow(nb))
  
  
  
  file.remove("test.vcf");file.remove("test_nb.vcf")
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
                                                     rows_per_individual = 2, header_cols = 1, mDat = -9)))
  
  file.remove("test.str")
})

test_that("non-biallelic",{
  td <- read_non_biallelic(genotypes(steelMSATs), sample.meta = sample.meta(steelMSATs))
  expect_true(is.snpRdata(td))
  expect_false(.is.bi_allelic(td))
  expect_equivalent(td, steelMSATs)
})


test_that("plink",{
  format_snps(stickSNPs, output = "plink", outfile = "test")
  expect_warning(test <- read_plink("test"))
  
  expect_true("position" %in% colnames(snp.meta(test)))
  expect_true("chr" %in% colnames(snp.meta(test)))
  
  test <- calc_ho(test)
  t2 <- calc_ho(stickSNPs)
  test <- get.snpR.stats(test, stats = "ho")$single
  t2 <- get.snpR.stats(t2, stats = "ho")$single
  tm <- merge(test, t2, by = c("chr", "position"))
  expect_identical(tm$ho.x, tm$ho.y)
  
  file.remove("test.bed", "test.map", "test.fam", "test.ped", "test.bim")
})


test_that("genlight and genind",{
  # gi
  gi <- format_snps(stickSNPs, output = "adegenet")
  expect_warning(test <- convert_genind(gi))
  
  test <- calc_ho(test)
  t2 <- calc_ho(stickSNPs)
  testm <- get.snpR.stats(test, stats = "ho")$weighted.means
  t2m <- get.snpR.stats(t2, stats = "ho")$weighted.means
  expect_identical(testm, t2m)
  

  # gl--note this wraps into gi first
  num <- format_snps(stickSNPs, "sn", interpolate = FALSE)
  gl <- methods::as(t(num[,-c(1:2)]), "genlight")
  expect_warning(.make_it_quiet(t3 <- convert_genlight(gl)))
  t3 <- calc_ho(t3)
  t3m <- get.snpR.stats(t3, stats = "ho")$weighted.means
  expect_identical(testm, t3m)
})

test_that("indels",{
  format_snps(stickSNPs, "vcf", outfile = "test.vcf")
  write("groupXXI\t9067879\tINDEL1\tT\tTAAGGACA\t.\tPASS\tNS=130;AC=17\tGT\t0/0\t./.\t0/0\t./.\t0/1\t0/0\t0/0\t0/0\t0/0\t0/1\t0/1\t0/0\t./.\t./.\t0/1\t./.\t0/1\t./.\t./.\t./.\t0/0\t./.\t./.\t0/1\t./.\t./.\t./.\t0/1\t0/0\t0/0\t0/0\t./.\t0/0\t./.\t0/1\t0/0\t0/0\t./.\t./.\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/1\t0/0\t./.\t./.\t0/0\t0/1\t./.\t./.\t0/0\t0/0\t./.\t0/0\t./.\t0/0\t0/1\t0/0\t./.\t./.\t0/1\t0/0\t0/0\t0/0\t./.\t0/0\t0/0\t./.\t./.\t0/0\t./.\t0/0\t0/0\t0/0\t0/1\t0/0\t0/0\t0/0\t./.\t0/1\t0/0\t0/0\t0/0\t0/1\t0/0\t0/0\t0/0\t./.\t./.\t0/1\t./.\t0/0\t0/1\t0/0\t./.\t0/0\t./.\ngroupXXI\t9067879\tBADL1\tTASDSA\tT\t.\tPASS\tNS=130;AC=17\tGT\t0/0\t./.\t0/0\t./.\t0/1\t0/0\t0/0\t0/0\t0/0\t0/1\t0/1\t0/0\t./.\t./.\t0/1\t./.\t0/1\t./.\t./.\t./.\t0/0\t./.\t./.\t0/1\t./.\t./.\t./.\t0/1\t0/0\t0/0\t0/0\t./.\t0/0\t./.\t0/1\t0/0\t0/0\t./.\t./.\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/1\t0/0\t./.\t./.\t0/0\t0/1\t./.\t./.\t0/0\t0/0\t./.\t0/0\t./.\t0/0\t0/1\t0/0\t./.\t./.\t0/1\t0/0\t0/0\t0/0\t./.\t0/0\t0/0\t./.\t./.\t0/0\t./.\t0/0\t0/0\t0/0\t0/1\t0/0\t0/0\t0/0\t./.\t0/1\t0/0\t0/0\t0/0\t0/1\t0/0\t0/0\t0/0\t./.\t./.\t0/1\t./.\t0/0\t0/1\t0/0\t./.\t0/0\t./.\ngroupXXI\t9067879\tINDEL2\tTAGCA\tT\t.\tPASS\tNS=130;AC=17\tGT\t0/0\t./.\t0/0\t./.\t0/1\t0/0\t0/0\t0/0\t0/0\t0/1\t0/1\t0/0\t./.\t./.\t0/1\t./.\t0/1\t./.\t./.\t./.\t0/0\t./.\t./.\t0/1\t./.\t./.\t./.\t0/1\t0/0\t0/0\t0/0\t./.\t0/0\t./.\t0/1\t0/0\t0/0\t./.\t./.\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/1\t0/0\t./.\t./.\t0/0\t0/1\t./.\t./.\t0/0\t0/0\t./.\t0/0\t./.\t0/0\t0/1\t0/0\t./.\t./.\t0/1\t0/0\t0/0\t0/0\t./.\t0/0\t0/0\t./.\t./.\t0/0\t./.\t0/0\t0/0\t0/0\t0/1\t0/0\t0/0\t0/0\t./.\t0/1\t0/0\t0/0\t0/0\t0/1\t0/0\t0/0\t0/0\t./.\t./.\t0/1\t./.\t0/0\t0/1\t0/0\t./.\t0/0\t./.\ngroupXXI\t9067879\tLONGAL\tTAGCA\tTGGCC\t.\tPASS\tNS=130;AC=17\tGT\t0/0\t./.\t0/0\t./.\t0/1\t0/0\t0/0\t0/0\t0/0\t0/1\t0/1\t0/0\t./.\t./.\t0/1\t./.\t0/1\t./.\t./.\t./.\t0/0\t./.\t./.\t0/1\t./.\t./.\t./.\t0/1\t0/0\t0/0\t0/0\t./.\t0/0\t./.\t0/1\t0/0\t0/0\t./.\t./.\t0/0\t0/0\t0/0\t0/0\t0/0\t0/0\t0/1\t0/0\t./.\t./.\t0/0\t0/1\t./.\t./.\t0/0\t0/0\t./.\t0/0\t./.\t0/0\t0/1\t0/0\t./.\t./.\t0/1\t0/0\t0/0\t0/0\t./.\t0/0\t0/0\t./.\t./.\t0/0\t./.\t0/0\t0/0\t0/0\t0/1\t0/0\t0/0\t0/0\t./.\t0/1\t0/0\t0/0\t0/0\t0/1\t0/0\t0/0\t0/0\t./.\t./.\t0/1\t./.\t0/0\t0/1\t0/0\t./.\t0/0\t./.", 
        "test.vcf", 
        append = TRUE)
  expect_warning(d2 <- read_vcf("test.vcf", sample.meta = sample.meta(stickSNPs)))
  expect_true(all(colnames(d2@geno.tables$gs) %in% c('AA','AC','AG','AT','CC','CG','CT','DD','DI','GG','GT','II','LR','RR','TT')))
  expect_true(all(colnames(d2@geno.tables$as) %in% c('A', 'G', 'C', 'T', 'I', 'D', 'L', 'R')))
  
  d2 <- calc_pi(d2, "pop")
  d2 <- calc_ho(d2, "pop")
  d2 <- calc_pairwise_fst(d2, "pop")
  expect_false("BADL1" %in% snp.meta(d2)$ID)
  res <- get.snpR.stats(d2, "pop", c("pi", "ho", "fst"))
  single <- res$single[res$single$ID %in% c("SNP_39", "INDEL1", "INDEL2","LONGAL") & res$single$subfacet == "ASP",]
  expect_equal(length(unique(single$pi)), 1)
  expect_equal(length(unique(single$ho)), 1)
  pairwise <- res$pairwise[res$pairwise$ID %in% c("SNP_39", "INDEL1", "INDEL2","LONGAL") & res$pairwise$comparison == "ASP~CLF",]
  expect_equal(length(unique(pairwise$fst)), 1)
  file.remove("test.vcf")
})

test_that("sanity checks", {
  
  # bad mDat, genotype format
  ## inconsistant genotype format
  test <- genotypes(stickSNPs)
  test[6,1] <- "AACC"
  expect_error(import.snpR.data(test, snp.meta(stickSNPs), sample.meta(stickSNPs)), "All genotypes must be equal in length, including missing data.")

  ## short mDat
  expect_error(import.snpR.data(genotypes(stickSNPs), snp.meta(stickSNPs), sample.meta(stickSNPs), mDat = "N"),
               "argument correctly given as the GENOTYPE rather than ALLELE entry for missing data")
  
  ## mDat with different "alleles"
  expect_error(import.snpR.data(genotypes(stickSNPs), snp.meta(stickSNPs), sample.meta(stickSNPs), mDat = "NX"),
               "mDat must be symmetrical, as in 'NN' or '0000', NOT '01' or 'NGT' or 'NX'")
  
  # odd snp/sample.meta
  expect_error(import.snpR.data(genotypes(stickSNPs), c("HI", "BYE"), sample.meta(stickSNPs)),
               "Number of rows in snp.meta")
  
  expect_error(import.snpR.data(genotypes(stickSNPs), c("HI"), sample.meta(stickSNPs)),
               "Cannot locate snp.meta file")
  
  expect_error(import.snpR.data(genotypes(stickSNPs), sample.meta = c("HI", "BYE")),
               "Number of rows in sample.meta")
  
  expect_error(import.snpR.data(genotypes(stickSNPs), sample.meta = c("HI")),
               "Cannot locate sample.meta file")
})
