#' Convert genind object to snpRdata.
#'
#' Convert genind object to snpRdata. Internal, called by import.snpR.data when
#' provided a adegenet genind object.
#'
#' @param genind genind object.
#' @param snp.meta data.frame, default NULL. Metadata for each SNP, IDs from
#'   genind may be added.
#' @param sample.meta data.frame, default NULL. Metadata for each individual
#'   sample, IDs and pops from genind may be added.
#'
#' @author William Hemstrom
genind.to.snpRdata <- function(genind, snp.meta = NULL, sample.meta = NULL){
  #========sanity checks=============
  msg <- character(0)
  
  nal <- lapply(adegenet::alleles(genind), length) > 2
  if(any(nal)){
    msg <- c(msg, paste0("Some non-SNP alleles (>2 loci) detected in genind object: ",
                         paste0(adegenet::locNames(genind)[which(nal)], collapse = ", "),
                         "\nsnpR is designed to only work with SNP data.\n")
    )
  }
  
  unique.alleles <- unique(unlist(adegenet::alleles(genind)))
  if(any(!unique.alleles %in% c("A", "T", "C", "G"))){
    msg <- c(msg, "Some invalid allele names in the provided genind object. For now, alleles must be named A, C, G, or T.\n Names can be updated with adegenet::alleles(genotypes) <- .\n")
  }
  
  if(!is.null(snp.meta)){
    
    if(nrow(snp.meta) != adegenet::nLoc(genind)){
      snp.meta <- NULL
      warning("Provided SNP metadata does not have the same number of rows as the number of loci in provided genind object, will be discarded.\n")
    }
  }
  
  if(!is.null(sample.meta)){
    if(nrow(sample.meta) != adegenet::nInd(genind)){
      sample.meta <- NULL
      warning("Provided sample metadata does not have the same number of rows as the number of samples in provided genind object, will be discarded.\n")
    }
  }
  if(any(adegenet::ploidy(genind) != 2)){
    msg <- msg <- c(msg, "For now, snpR only converts diploid genind objects to snpRdata automatically.\n")
  }
  
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  #==========pull out stuff for snpR=====
  # genotypes
  ## approach: figure out where we have homozygotes and fill new genotype calls with correct identity, then do hets and missing.
  genotypes <- adegenet::tab(genind)
  alleles <- adegenet::alleles(genind)
  allele1ident <- rep(unlist(purrr::map(alleles, 1)), each = nrow(genotypes))
  allele2ident <- rep(unlist(purrr::map(alleles, 2)), each = nrow(genotypes))
  a1c <- genotypes[,seq(1, ncol(genotypes), by = 2)]
  a2c <- genotypes[,seq(2, ncol(genotypes), by = 2)]
  new.genotypes <- character(length(a1c))
  
  ## where we have a homozygote allele 1, fill
  a1homs <- which(a1c == 2)
  a2homs <- which(a2c == 2)
  nas <- which(is.na(a1c))
  not_hets <- union(union(a1homs, a2homs), nas)
  new.genotypes[a1homs] <- paste0(allele1ident[a1homs], allele1ident[a1homs])
  new.genotypes[a2homs] <- paste0(allele2ident[a2homs], allele2ident[a2homs])
  new.genotypes[nas] <- "NN"
  new.genotypes[-not_hets] <- paste0(allele1ident[-not_hets], allele2ident[-not_hets])
  genotypes <- as.data.frame(t(matrix(new.genotypes, nrow(genotypes), ncol(genotypes)/2)), stringsAsFactors = F)
  
  
  # snp metadata, needs to be provided, albeit we can add on snp names
  if(is.null(snp.meta)){
    snp.meta <- data.frame(snpID = adegenet::locNames(genind))
  }
  else if(!"snpID" %in% colnames(snp.meta)){
    snp.meta$ID <- adegenet::locNames(genind)
  }
  
  # sample metadata
  if(is.null(sample.meta)){
    sample.meta <- data.frame(sampID = adegenet::indNames(genind))
    sample.meta$pop <- adegenet::pop(genind)
  }
  else{
    if(!"sampID" %in% colnames(sample.meta)){
      sample.meta$sampID <- adegenet::indNames(genind)
    }
    if(!"pop" %in% colnames(sample.meta)){
      sample.meta$pop <- adegenet::pop(genind)
    }
  }
  
  #=======format and return=========
  return(import.snpR.data(genotypes, snp.meta, sample.meta))
}

#' Convert genlight object to snpRdata.
#'
#' Convert genlight object to snpRdata. Internal, called by import.snpR.data
#' when provided a adegenet genlight object.
#'
#' @param genlight genlight object.
#' @param snp.meta data.frame, default NULL. Metadata for each SNP, IDs from
#'   genind may be added.
#' @param sample.meta data.frame, default NULL. Metadata for each individual
#'   sample, IDs and pops from genind may be added.
#'
#' @author William Hemstrom
genlight.to.snpRdata <- function(genlight, snp.meta = NULL, sample.meta = NULL){
  #========sanity checks=============
  msg <- character(0)
  
  unique.alleles <- unique(adegenet::alleles(genlight))
  unique.alleles <- unique(unlist(strsplit(unique.alleles, "/")))
  unique.alleles <- toupper(unique.alleles)
  if(any(!unique.alleles %in% c("A", "T", "C", "G"))){
    msg <- c(msg, "Some invalid allele names in the provided genind object. For now, alleles must be named A, C, G, or T.\n Names can be updated with adegenet::alleles(genotypes) <- .\n")
  }
  
  if(!is.null(snp.meta)){
    if(nrow(snp.meta) != adegenet::nLoc(genlight)){
      snp.meta <- NULL
      warning("Provided SNP metadata does not have the same number of rows as the number of loci in provided genind object, will be discarded.\n")
    }
  }
  
  if(!is.null(sample.meta)){
    if(nrow(sample.meta) != adegenet::nInd(genlight)){
      sample.meta <- NULL
      warning("Provided sample metadata does not have the same number of rows as the number of samples in provided genind object, will be discarded.\n")
    }
  }
  
  if(any(adegenet::ploidy(genlight) != 2)){
    msg <- msg <- c(msg, "For now, snpR only converts diploid genind objects to snpRdata automatically.\n")
  }
  
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  #=========pull out parts for snpRdata========
  # genotypes
  ## get alleles at each locus
  genotypes <- adegenet::tab(genlight, NA.method = "asis")
  als <- adegenet::alleles(genlight)
  a1s <- gsub("/.", "", als)
  a1s <- rep(a1s, each = nrow(genotypes))
  a2s <- gsub("./", "", als)
  a2s <- rep(a2s, each = nrow(genotypes))
  
  ## fill in homs and hets with correct alleles.
  hom_a1s <- which(genotypes == 0)
  hets <- which(genotypes == 1)
  hom_a2s <- which(genotypes == 2)
  new_genos <- character(length(genotypes))
  new_genos[hom_a1s] <- paste0(a1s[hom_a1s], a1s[hom_a1s])
  new_genos[hets] <- paste0(a1s[hets], a2s[hets])
  new_genos[hom_a2s] <- paste0(a2s[hom_a2s], a2s[hom_a2s])
  new_genos[which(is.na(genotypes))] <- "NN"
  new_genos <- toupper(new_genos)
  genotypes <- data.frame(t(matrix(new_genos, nrow(genotypes), ncol(genotypes))), stringsAsFactors = F)
  
  
  # snp metadata, whatever possible
  if(is.null(snp.meta)){
    snp.meta <- data.frame(snpID = adegenet::locNames(genlight))
  }
  have.meta.cols <- colnames(snp.meta)
  
  if(!"snpID" %in% have.meta.cols){
    snp.meta$snpID <- adegenet::locNames(genlight)
  }
  
  if(!"chr" %in% have.meta.cols){
    snp.meta$chr <- adegenet::chr(genlight)
  }
  if(!"position" %in% have.meta.cols){
    snp.meta$position <- adegenet::position(genlight)
  }
  
  # sample metadata, whatever possible
  if(is.null(sample.meta)){
    sample.meta <- data.frame(sampID = adegenet::indNames(genlight))
  }
  have.meta.cols <- colnames(sample.meta)
  
  if(!"sampID" %in% have.meta.cols){
    sample.meta$sampID <- adegenet::indNames(genlight)
  }
  if(!"pop" %in% have.meta.cols){
    sample.meta$pop <- adegenet::pop(genlight)
  }
  
  #=========format and return===========
  return(import.snpR.data(genotypes, snp.meta, sample.meta, mDat = "NN"))
}


#' Internal to process a ms file
#' @param x filepath to ms file
#' @param chr.length length of the chromosome. If a single value, assumes all
#'   the same length. If a vector of the same length as number of chr, assumes
#'   those are the chr lengths in order of apperance in ms file.
#' @author William Hemstrom
process_ms <- function(x, chr.length){
  infile <- x #infile
  lines <- readLines(x)
  lines <- lines[-which(lines == "")] #remove empty entries
  lines <- lines[-c(1,2)] #remove header info
  nss <- grep("segsites", lines) #get the number of segsites per chr
  chrls <- gsub("segsites: ", "", lines[nss]) #parse this to get the lengths
  chrls <- as.numeric(chrls)
  lines <- lines[-nss] #remove the segsites lines
  pos <- lines[grep("positions:", lines)] #find the positions
  lines <- lines[-grep("positions:", lines)] #remove the position
  div <- grep("//", lines) #find the seperators
  gc <- div[2] - div[1] - 1 #find the number of gene copies per chr
  if(is.na(gc)){gc <- length(lines) - 1} #if there's only one chr
  dat <- lines[-div] #get the data only
  dat <- strsplit(dat, "") #split the lines by individual snp calls
  x <- matrix(NA, nrow = sum(chrls), ncol = gc) #prepare output
  meta <- matrix(NA, nrow = sum(chrls), 2)
  
  #process this into workable data
  pchrls <- c(0, chrls)
  pchrls <- cumsum(pchrls)
  
  # check lengths input
  if(length(chr.length) != 1){
    if(length(chr.length) != length(chrls)){
      stop("Provided vector of chromosome lengths is not equal to the number of chromosomes in ms file.\n")
    }
  }
  
  for(i in 1:length(chrls)){
    cat("\n\tChr ", i)
    tg <- dat[(gc*(i-1) + 1):(gc*i)] #get only this data
    tg <- unlist(tg) #unlist
    tg <- matrix(as.numeric(tg), ncol = chrls[i], nrow = gc, byrow = T) #put into a matrix
    tg <- t(tg) #transpose. rows are now snps, columns are gene copies
    tpos <- unlist(strsplit(pos[i], " ")) #grap and process the positions
    tpos <- tpos[-1]
    meta[(pchrls[i] + 1):pchrls[i + 1],] <- cbind(paste0(rep("chr", length = nrow(tg)), i), tpos)
    x[(pchrls[i] + 1):pchrls[i + 1],] <- tg #add data to output
  }
  
  meta <- as.data.frame(meta, stringsAsFactors = F)
  meta[,2] <- as.numeric(meta[,2])
  if(length(chr.length) == 1){
    meta[,2] <- meta[,2] * chr.length
  }
  else{
    meta[,2] <- meta[,2] * chr.length[as.numeric(substr(meta[,1], 4, 4))] # multiply by the correct chr length.
  }
  
  colnames(meta) <- c("group", "position")
  colnames(x) <- paste0("gc_", 1:ncol(x))
  
  return(list(x = x, meta = meta))
}

#' Convert a vcf file/vcfR object into a snpRdata object
#'
#' @param vcf_file character or vcfR object. Either a path to a vcf file or a
#'   vcfR object.
#' @param snp.meta data.frame or null, default null. snp metadata, will
#'   overwrite data in vcf if provided.
#' @param sample.meta data.frame or null, default null. sample metadata, will
#'   overwrite sample names in vcf if provided.
#'
#' @author William Hemstrom
process_vcf <- function(vcf_file, snp.meta = NULL, sample.meta = NULL){
  #========sanity checks, part 1=============
  check.installed("vcfR")
  
  #========import data=======================
  if(!"vcfR" %in% class(vcf_file)){
    vcf <- vcfR::read.vcfR(vcf_file)
  }
  else{
    vcf <- vcf_file
    rm(vcf_file)
  }
  
  #========sanity checks, part2==============
  msg <- character(0)
  warn <- character(0)
  
  # initialize bad loci storage
  good.loci <- logical(vcfR::nrow(vcf))
  good.loci <- !good.loci
  
  # check for lack of called genotypes
  formats <- vcf@gt[,"FORMAT"]
  no.genotypes <- which(!grepl("GT", formats))
  if(length(no.genotypes) > 0){
    warn <- c(warn, paste0(length(no.genotypes), " loci removed due to missing called genotypes.\n"))
    good.loci[no.genotypes] <- F
  }
  rm(formats)
  
  # check for not snps
  ref <- vcfR::getREF(vcf)
  alt <- vcfR::getALT(vcf)
  ok.alleles <- c(".", "A", "C", "G", "T")
  bad.ref <- !ref %in% ok.alleles
  bad.alt <- !alt %in% ok.alleles
  bad.either <- which(bad.ref | bad.alt)
  if(length(bad.either) > 1){
    warn <- c(warn, paste0(length(bad.either), " loci removed due to improper alleles (not ., A, C, T, or G).\n"))
    good.loci[bad.either] <- F
  }
  rm(ref, alt, ok.alleles, bad.ref, bad.alt, bad.either)
  
  # anything left?
  if(sum(good.loci) == 0){
    msg <- c(msg, "No loci remain after filtering.\n")
  }
  
  # check metadata
  if(!is.null(snp.meta)){
    if(nrow(snp.meta) != vcfR::nrow(vcf)){
      msg <- c(msg, "Number of rows in provided SNP meta not equal to number of SNPs in vcf.\n")
    }
  }
  if(!is.null(sample.meta)){
    if(nrow(sample.meta) != (ncol(vcf@gt)) - 1){
      msg <- c(msg, "Number of rows in provided sample meta not equal to number of samples in vcf.\n")
    }
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  #==========convert================
  # prep genotypes
  genos <- vcfR::extract.gt(vcf, element = "GT", return.alleles = T,  IDtoRowNames = F)
  genos <- genos[good.loci,]
  genos <- gsub("\\|", "", genos)
  genos <- gsub("\\/", "", genos)
  genos[is.na(genos)] <- "NN"
  
  # prep metadata
  if(is.null(snp.meta)){
    snp.meta <- vcfR::getFIX(vcf)
    snp.meta.col <- which(colnames(snp.meta) == "POS")
    colnames(snp.meta)[snp.meta.col] <- "position" # fix this specifically, since many functions want this name
    snp.meta.col <- which(colnames(snp.meta) == "POS")
  }
  snp.meta <- snp.meta[good.loci,]
  
  if(is.null(sample.meta)){
    sample.meta <- data.frame(sampID = colnames(genos))
  }
  
  # make the snpRdata object
  return(import.snpR.data(as.data.frame(genos), 
                          as.data.frame(snp.meta), 
                          as.data.frame(sample.meta), "NN"))
}

#' Convert a genepop into a snpRdata object.
#'
#' @param genepop_file character. A path to a genepop file.
#' @param snp.meta data.frame or null, default null. snp metadata.
#' @param sample.meta data.frame or null, default null. sample metadata
#' @param mDat character, default "0000". Missing data format. If the genotype
#'   is 6 instead of 4 characters, will be automatically changed to 000000 IF
#'   the default value is provided.
#'
#' @author William Hemstrom
process_genepop <- function(genepop_file, snp.meta = NULL, sample.meta = NULL, mDat = "0000"){

  #===============import data==============
  dat <- readLines(genepop_file)
  
  # find the pop headers
  poplines <- grep("POP", dat, ignore.case = TRUE)
  poplines <- poplines[which(nchar(dat[poplines]) == 3)] # ONLY contain the word: POP
  
  
  #===============prep genotypes===========
  gt <- dat[-c(1:(poplines[1] - 1), poplines)]
  gt <- strsplit(gt, split = "\\s")
  gt <- data.frame(gt)
  colnames(gt) <- gt[1,]
  gt <- gt[-c(1:2),]
  
  #===============prep sample metadata===========
  ps <- poplines - ((1:length(poplines)) + (poplines[1] - 2))
  samples_per_pop <- c(ps[-1], ncol(gt) + 1) - ps
  gp_pops <- colnames(gt)[ps] # pop headers
  
  if(!is.null(sample.meta)){
    if("genepop_pops" %in% colnames(sample.meta)){
      warning("genepop_pops column detected in provided sample.meta. This column will be overwritten with genepop pop labels.\n")
    }
    if(nrow(sample.meta) != ncol(gt)){
      stop("Number of samples in provided sample metadata does not match the number of samples in the genepop file.\n")
    }
    sample.meta$genepop_pops <- rep(gp_pops, samples_per_pop)
  }
  else{
    sample.meta <- data.frame(sampleID = colnames(gt), genepop_pop = rep(gp_pops, samples_per_pop))
  }
  
  #===============prep SNP metadata===========
  if(poplines[1] == 3){
    snp.ids <- unlist(strsplit(dat[2], ",\\s"))
  }
  else{
    snp.ids <- dat[2:(poplines[1] - 1)]
  }
  
  if(length(snp.ids) != nrow(gt)){
    stop("Number of SNP IDs in the provided genepop file does not match the number of SNPs in the dataset.\n")
  }
  
  if(!is.null(snp.meta)){
    if("genepop_snpIDs" %in% colnames(snp.meta)){
      warning("genepop_snpIDs column detected in provided sample.meta. This column will be overwritten with genepop SNP labels.\n")
    }
    if(nrow(snp.meta) != nrow(gt)){
      stop("Number of samples in provided sample metadata does not match the number of samples in the genepop file.\n")
    }
    snp.meta$genepop_snpIDs <- snp.ids
  }
  else{
    snp.meta <- data.frame(genepop_snpIDs = snp.ids)
  }
  
  #=============convert to snpRdata object=====
  form <- nchar(gt[1,1])
  if(!form %in% c(4, 6)){
    stop("Genepop inputs must have genotypes formatted as either 2 or 3 numeric values per allele (0101 or 123123, for example).\n")
  }
  if(mDat == "0000" & form == 6){
    mDat <- "000000"
    warning("6 character genotypes detected, switching missing data entry to 000000. To stop this, set the mDat argument to the correct value.\n")
  }
  
  # process allele 1 and allele 2
  a1 <- substr(unlist(gt), 1, form/2)
  a1 <- matrix(a1, nrow(gt), ncol(gt))
  a1[gt == mDat] <- NA
  
  a2 <- substr(unlist(gt), (form/2) + 1, form)
  a2 <- matrix(a2, nrow(gt), ncol(gt))
  a2[gt == mDat] <- NA
  
  
  # combine to check for the correct number of alleles per loci (using factor levels for speed)
  ac <- cbind(a1, a2)
  ac[ac == substr(mDat, 1, form/2)] <- NA
  ac <- as.data.frame(t(ac), stringsAsFactors = T)
  ac <- dplyr::mutate_all(ac, as.numeric)
  ac <- as.matrix(ac)
  
  if(any(matrixStats::colMaxs(ac, na.rm = T) > 2)){
    stop("Some loci with more than two alleles detected. snpR only accepts SNP (or SNP-like) data with two alleles per loci. This can happen in biallelic data if the missing data value is not properly set.\n")
  }
  
  # sub into NN and finish up
  ac[ac == 1] <- "A"
  ac[ac == 2] <- "C"
  ac[is.na(ac)] <- "N"
  ac <- paste0(ac[1:ncol(gt),], ac[(ncol(gt) + 1 ):nrow(ac),])
  ac <- matrix(ac, nrow(gt), ncol(gt), byrow = T)
  
  return(import.snpR.data(ac, snp.meta, sample.meta))
}