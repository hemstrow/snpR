test_that("merging",{
  gy <- data.frame(s1 = c("GG", "NN"),
                  s2 = c("GG", "TG"),
                  s3 = c("NN", "TT"),
                  s4 = c("GA", "TT"),
                  s5 = c("GG", "GT"),
                  s6 = c("NN", "GG"))
  
  snp.y <- data.frame(chr = c("groupVI", "test_chr"),
                      position = c(212436, 10))
  
  samp.y <- data.frame(pop = c("ASP", "ASP", "ASP", "test1", "test2", "test3"),
                       ID = c(1, 2, 3, "A1", "A2", "A3"),
                       fam = c("A", "B", "C", "T", "T", "T"))
  
  y <- import.snpR.data(gy, snp.y, samp.y)
  
  x <- stickSNPs
  sample.meta(x)$ID <- 1:ncol(x)
  
  # error calling
  expect_error(z <- merge_snpRdata(x, y), "Some genotypes at identical loci sequenced in samples in both 'x' and 'y' are not identical") 
  expect_warning(z <- merge_snpRdata(x, y, resolve_conflicts = "warning"), "enotypic mismatches in identical samples/snps. Returning matrix of mismatches.")
  expect_true(nrow(z) == 1)
  
  comp <- genotypes(x)[which(snp.meta(x)$chr == "groupVI" & snp.meta(x)$position == 212436),1:2] == c(gy[1,1:2])
  expect_true(all(z[,4:6] == c(comp, TRUE))) # also checks NN replacement
  
  # correct merging
  ## merge favoring x
  z <- merge_snpRdata(x, y, resolve_conflicts = "x")
  expect_equal(nrow(z), nrow(x) + 1)
  expect_equal(ncol(z), ncol(x) + 3)
  
  comp2 <- genotypes(z)[which(snp.meta(z)$chr == "groupVI" & snp.meta(z)$position == 212436),]
  expect_equivalent(comp2[101:103], c("GA", "GG", "NN"))
  expect_equivalent(comp2[1:3], genotypes(x)[1,1:3])
  
  ## merge favoring y
  z <- merge_snpRdata(x, y, resolve_conflicts = "y")
  
  comp2 <- genotypes(z)[which(snp.meta(z)$chr == "groupVI" & snp.meta(z)$position == 212436),]
  expect_equivalent(comp2[101:103], c("GA", "GG", "NN"))
  expect_equivalent(comp2[1:3], c(genotypes(y)[2,1:2], genotypes(x)[1,3])) # also check smart NN replacement!
  
  ## random favoring
  ## merge favoring y
  set.seed(1232)
  z <- merge_snpRdata(x, y, resolve_conflicts = "random")
  
  comp2 <- genotypes(z)[which(snp.meta(z)$chr == "groupVI" & snp.meta(z)$position == 212436),]
  comp2[1:3]
  
  expect_equivalent(comp2[101:103], c("GA", "GG", "NN"))
  expect_equivalent(comp2[1:3], c("GG", "GG", "AA")) # also check smart NN replacement!
  
  
  set.seed(1234)
  z <- merge_snpRdata(x, y, resolve_conflicts = "random")
  
  comp2 <- genotypes(z)[which(snp.meta(z)$chr == "groupVI" & snp.meta(z)$position == 212436),]
  comp2[1:3]
  
  expect_equivalent(comp2[101:103], c("GA", "GG", "NN"))
  expect_equivalent(comp2[1:3], c("AG", "GG", "AA")) # also check smart NN replacement!
  
  
  # all options
  z <- merge_snpRdata(x, y, all.x.snps = FALSE, resolve_conflicts = "x")
  expect_equal(dim(z), c(2, 103))
  z <- merge_snpRdata(x, y, all.y.snps = FALSE, resolve_conflicts = "x")
  expect_equal(dim(z), c(100, 103))
  z <- merge_snpRdata(x, y, all.x.samples = FALSE, resolve_conflicts = "x")
  expect_equal(dim(z), c(101, 6))
  z <- merge_snpRdata(x, y, all.y.samples = FALSE, resolve_conflicts = "x")
  expect_equal(dim(z), c(101, 100))
  z <- merge_snpRdata(x, y, all = FALSE, resolve_conflicts = "x")
  expect_equal(dim(z), c(1, 3))
  z <- merge_snpRdata(x, y, all.x.samples = FALSE, all.y.snps = FALSE, resolve_conflicts = "x")
  expect_equal(dim(z), c(100, 6))
  
})

test_that("merging errors",{
  # bad matches in SNP metadata
  gy <- data.frame(s1 = c("GG", "NN"),
                   s2 = c("GG", "TG"),
                   s3 = c("NN", "TT"),
                   s4 = c("GA", "TT"),
                   s5 = c("GG", "GT"),
                   s6 = c("NN", "GG"))
  
  snp.y <- data.frame(CHROM = c("groupVI", "test_chr"),
                      pos = c(212436, 10))
  
  samp.y <- data.frame(pop = c("ASP", "ASP", "ASP", "test1", "test2", "test3"),
                       ID = c(1, 2, 3, "A1", "A2", "A3"),
                       fam = c("A", "B", "C", "T", "T", "T"))
  
  y <- import.snpR.data(gy, snp.y, samp.y)
  
  x <- stickSNPs
  sample.meta(x)$ID <- 1:ncol(x)
  
  z <- merge_snpRdata(x, y)
  expect_equal(nrow(z), 102) # treats snps as no overlap
  expect_error(merge_snpRdata(x, y, all.x.snps = FALSE, all.y.snps = FALSE), "No matching snp metadata")
  
  # bad matches in sample metadata
  gy <- data.frame(s1 = c("GG", "NN"),
                   s2 = c("GG", "TG"),
                   s3 = c("NN", "TT"),
                   s4 = c("GA", "TT"),
                   s5 = c("GG", "GT"),
                   s6 = c("NN", "GG"))
  
  snp.y <- data.frame(chr = c("groupVI", "test_chr"),
                      position = c(212436, 10))
  
  samp.y <- data.frame(pops = c("ASP", "ASP", "ASP", "test1", "test2", "test3"),
                       IDs = c(1, 2, 3, "A1", "A2", "A3"),
                       fams = c("A", "B", "C", "T", "T", "T"))
  
  y <- import.snpR.data(gy, snp.y, samp.y)
  expect_warning(z <- merge_snpRdata(x, y), "Some levels are duplicated")
  expect_equal(ncol(z), 106)
  expect_error(merge_snpRdata(x, y, all.x.samples = FALSE, all.y.samples = FALSE), "No matching sample metadata")
  
  
  # nothing remains -- snps
  gy <- data.frame(s1 = c("GG", "NN"),
                   s2 = c("GG", "TG"),
                   s3 = c("NN", "TT"),
                   s4 = c("GA", "TT"),
                   s5 = c("GG", "GT"),
                   s6 = c("NN", "GG"))
  
  snp.y <- data.frame(chr = c("groupVI", "test_chr"),
                      position = c(212437, 11))
  
  samp.y <- data.frame(pop = c("ASP", "ASP", "ASP", "test1", "test2", "test3"),
                       ID = c(1, 2, 3, "A1", "A2", "A3"),
                       fam = c("A", "B", "C", "T", "T", "T"))
  
  y <- import.snpR.data(gy, snp.y, samp.y)
  expect_error(merge_snpRdata(x, y, all.x.snps = FALSE, all.y.snps = FALSE), "No loci remain after merging")
  
  # nothing remains -- samps
  gy <- data.frame(s1 = c("GG", "NN"),
                   s2 = c("GG", "TG"),
                   s3 = c("NN", "TT"),
                   s4 = c("GA", "TT"),
                   s5 = c("GG", "GT"),
                   s6 = c("NN", "GG"))
  
  snp.y <- data.frame(chr = c("groupVI", "test_chr"),
                      position = c(212436, 10))
  
  samp.y <- data.frame(pop = c("hi", "hi", "hi", "hello", "hello", "hello"),
                       ID = c(1, 2, 3, "A1", "A2", "A3"),
                       fam = c("A", "B", "C", "T", "T", "T"))
  
  y <- import.snpR.data(gy, snp.y, samp.y)
  expect_error(merge_snpRdata(x, y, all.x.samples = FALSE, all.y.samples = FALSE), "No samples remain after merging")
})
  
