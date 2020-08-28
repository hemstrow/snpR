#===============prep some example data================
stickSNPs@sample.meta
stickSNPs@snp.meta
dat <- stickSNPs
s <- calc_maf(dat)

which.anc <- rbinom(nrow(dat), 1, .5)
dat@snp.meta$anc <- ""
dat@snp.meta$anc[which(which.anc == 0)] <- get.snpR.stats(s)$major[which(which.anc == 0)]
dat@snp.meta$anc[which(which.anc != 0)] <- get.snpR.stats(s)$minor[which(which.anc != 0)]
# the @snp.meta slot will now contain a column with the ancestral allele identity at each locus

dat <- import.snpR.data(dat[1:nrow(dat), 1:ncol(dat)], dat@snp.meta, dat@sample.meta)

#=================make the function=================
# here's an example of the structure of the function
abba_baba(dat, facet, levels = c(p1, p2, p3, pO)) # data, then the sample metadata column that defines populations/species, then the identity of the p1, p2, p3, pO populations
abba_baba(dat, "pop", levels = c("ASP", "CLF", "SMR", "OPL")) # might look like this!


# start by making a table of the ancestral and derived allele counts in each population

dat <- add.facets.snpR.data(dat, "pop")
temp <- cbind(dat@facet.meta, dat@geno.tables$as)
temp$anc_count <- 0
temp$der_count <- 0 # if you can fill these two columns, we should be good to calculate the abba/baba stats!


# this might help, this code does something really similar, but for major and minor alleles
get.ac <- function(x, maj, min, mis.al){
  # initialize:
  if(is.null(nrow(x))){
    temp.x <- matrix(x, ncol = length(x))
    colnames(temp.x) <- names(x)
    x <- temp.x
  }
  
  out <- data.frame(n_total = numeric(length(maj)),
                    n_alleles = numeric(length(maj)),
                    ni1 = numeric(length(maj)),
                    ni2 = numeric(length(maj)))
  
  
  # get the column from as matching the target allele.
  maj.col.match <- match(maj, colnames(x))
  out$ni1 <- t(x)[maj.col.match + seq(0, length(x) - ncol(x), by = ncol(x))]
  
  # ni1 is the rowsums minus this
  out$ni2 <- rowSums(x) - out$ni1
  out$n_total <- rowSums(x)
  out$n_alleles <- rowSums(ifelse(out[,3:4] != 0, 1, 0))
  out[is.na(out)] <- 0
  
  return(out)
}