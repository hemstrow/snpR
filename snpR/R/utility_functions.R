#Filters the data according to a set of filters. Vectorized to handle large datasets.
#Arguments:
# data: input data frame
#
#FILTERS:
# non_poly: filter out snps with only one homozygous genotype observed, list return is npf
# maf: filter out snps with a minor allele frequency below that given, list return is maff
# hf_hets: filters out data with a heterozygote frequency above the given value, list return is hfhf
# min_ind: filters out snps genotyped at less than the given number of individuals, list return is mif
#This version is vectorized, and is much faster!
filter_snps <- function(x, ecs, non_poly = FALSE, bi_al = TRUE, maf = FALSE, pop = FALSE,
                        hf_hets = FALSE, min_ind = FALSE, mDat = "NN"){

  #############################################################################################
  #do sanity checks

  if(!maf){
    if(!is.numeric(maf)){
      stop("maf must be a numeric value.")
    }
    if(length(maf) != 1){
      stop("maf must be a numeric vector of length 1.")
    }
  }

  if(!hf_hets){
    if(!is.numeric(hf_hets)){
      stop("hf_hets must be a numeric value.")
    }
    if(length(hf_hets) != 1){
      stop("hf_hets must be a numeric vector of length 1.")
    }
  }

  if(!min_ind){
    if(!is.numeric(min_ind)){
      stop("min_ind must be a numeric value.")
    }
    if(length(min_ind) != 1){
      stop("min_ind must be a numeric vector of length 1.")
    }
  }

  if(!pop){
    if(!is.numeric(pop[[2]])){
      stop("pop[[2]] must be a numeric vector.")
    }
    if(sum(pop[[2]]) != (ncol(x) - ecs)){
      stop("sum of pop[[2]] must be equal to number of individuals in x (check ecs as well).")
    }
    if(length(pop[[1]]) != length(pop[[2]])){
      stop("Number of provided pop names and pop sizes not equal.")
    }
  }

  ##############################################################################################
  ###set up, get values used later, clean up data a bit, set any functions used lower
  cat("Initializing...\n")

  #get headers
  headers <- x[,1:ecs]

  #seperate out the genotypes alone
  x <- x[,(ecs+1):ncol(x)]

  #get a single vector of genotypes
  xv <- as.vector(t(x))

  #determine snp format
  snp_form <- nchar(xv[1])

  #get number of samples
  nsamp <- ncol(x)

  #create tables of genotype counts for each locus.
  #function to do this, takes the pattern to match and the data, which is as a SINGLE VECTOR, by rows.
  #needs the number of samples from global function environment
  count_genos <- function(pattern, x){
    XXs <- grep(pattern, x) #figure out which entries match pattern
    out <- table(ceiling(XXs/nsamp)) #devide the entries that match by the number of samps to get which samp they come from (when rounded up), then table that.
    return(out)
  }

  #get all possible genotypes
  gs <- unique(xv)

  #which gs are heterozygous?
  hs <- substr(gs,1,snp_form/2) != substr(gs, (snp_form/2 + 1), snp_form*2)

  #which genotype is the missing data?
  mpos <- which(gs == mDat)

  #what are the possible alleles at all loci?
  as <- unique(c(substr(gs,1,snp_form/2), substr(gs, (snp_form/2 + 1), snp_form*2)))
  as <- as[as != substr(mDat, 1, snp_form/2)] #that aren't N?

  #############################################################################################3
  ###get a table of genotype and allele counts at each locus.
  cat("Creating genotype table...\n")

  #for each element of gs, get the tables of genotype counts and add them to a matrix
  gmat <- matrix(0, nrow(x), length(gs)) #initialize matrix
  colnames(gmat) <- gs #set the matrix names

  #fill the matrix, one possible genotype at a time (hard enough to vectorize as it is).
  for(i in 1:length(gs)){
    tab <- count_genos(gs[i], xv)
    gmat[as.numeric(names(tab)),i] <- as.numeric(tab)
  }

  if(length(mpos) > 0){
    tmat <- gmat[,-c(mpos)] #gmat without missing data
  }
  else{
    tmat <- gmat
  }

  #get matrix of allele counts
  #initialize
  cat("Getting allele table...\n")
  amat <- matrix(0, nrow(gmat), length(as))
  colnames(amat) <- as

  #fill in
  for(i in 1:length(as)){
    b <- grep(as[i], colnames(tmat))
    hom <- which(colnames(tmat) == paste0(as[i], as[i]))
    if(length(hom) == 0){
      het <- b
      amat[,i] <- rowSums(tmat[,het])
    }
    else{
      het <- b[b != hom]
      if(length(het) > 0){
        if(is.matrix(tmat[,het])){
          amat[,i] <- (tmat[,hom] * 2) + rowSums(tmat[,het])
        }
        else{
          amat[,i] <- (tmat[,hom] * 2) + tmat[,het]
        }
      }
      else{
        amat[,i] <- (tmat[,hom] * 2)
      }
    }
  }


  ##############################################################################################
  ###do filters
  keep <- logical(nrow(x)) #vector to track status

  #===============================================
  ##non-biallelic loci
  if(bi_al){
    cat("Filtering non-biallelic loci...\n")
    bimat <- ifelse(amat, TRUE, FALSE)
    bi <- ifelse(rowSums(bimat) > 2, T, F) #if true, should keep the allele
    #some tests require this, so subset the matrices and redefine things if true and some are multi-allelic
    if((maf | hf_hets) &  sum(bi) != 0){
      x <- x[bi,]
      xv <-   xv <- as.vector(t(x))
      gs <- unique(xv)
      hs <- substr(gs,1,snp_form/2) != substr(gs, (snp_form/2 + 1), snp_form*2)
      mpos <- which(gs == mDat)
      as <- unique(c(substr(gs,1,snp_form/2), substr(gs, (snp_form/2 + 1), snp_form*2)))
      as <- as[as != substr(mDat, 1, snp_form/2)] #that aren't N?
      gmat <- gmat[bi,]
      tmat <- tmat[bi,]
      amat <- amat[bi,]
      headers <- headers[bi,]
    }
    else{
      keep <- keep + bi
    }
  }

  #===============================================
  ##non-poly
  if(non_poly){
    cat("Filtering non-polymorphic loci...\n")
    npmat <- ifelse(tmat == 0, F, T) #convert to logical
    gcounts <- rowSums(npmat) #how many genotypes observed?
    hgcounts <- as.logical(rowSums(npmat[,hs[-mpos]])) #are hets observed?
    gcounts <- ifelse(gcounts == 1, F, T) #are there more than one genotype observed?
    np <- gcounts + hgcounts #figure out if either or both of the above are true
    np <- !as.logical(np) #convert to logical (if at least one is true, this will be false and we should keep the snp)
    keep <- keep + np
  }

  #===============================================
  ##min inds
  if(min_ind){
    cat("Filtering loci sequenced in few individuals...\n")
    mi <- ncol(x) - gmat[,colnames(gmat) == mDat]
    mi <- !mi >= min_ind #if false, enough samples, so keep.
    keep <- keep + mi
  }


  #===============================================
  ##minor allele frequency, both total and by pop. Should only run if bi_al = TRUE.
  if(maf){
    #if not filtering with multiple pops
    if(!is.list(pop)){
      cat("Filtering low minor allele frequencies, no pops...\n")
      mafs <- 1 - rowMaxs(amat)/rowSums(amat)
      mafs <- mafs < maf #less than required, set to true and reject.
      mafs[is.na(mafs)] <- TRUE
      keep <- keep + mafs
    }
    else{
      cat("Filtering low minor allele frequencies, pop:\n")
      pmafs <- logical(nrow(x))
      for(i in 1:length(pop[[1]])){
        cat(pop[[1]][i], "\n")
        #re-establish matrices with each pop

        #set the input data
        if(i == 1){
          popx <- x[,1:pop[[2]][i]]
        }
        else{
          popx <- x[,(sum(pop[[2]][1:(i-1)]) + 1):(sum(pop[[2]][1:i]))]
        }

        #get a single vector of genotypes
        popxv <- as.vector(t(popx))

        #get number of samples
        nsamp <- ncol(popx)

        #get all possible genotypes
        popgs <- unique(popxv)

        #which gs are heterozygous?
        pophs <- substr(popgs,1,snp_form/2) != substr(popgs, (snp_form/2 + 1), snp_form*2)

        #which genotype is the missing data?
        popmpos <- which(popgs == mDat)

        #what are the possible alleles at all loci?
        popas <- unique(c(substr(popgs,1,snp_form/2), substr(popgs, (snp_form/2 + 1), snp_form*2)))
        popas <- popas[popas != substr(mDat, 1, snp_form/2)] #that aren't N?

        #for each element of gs, get the tables of genotype counts and add them to a matrix
        popgmat <- matrix(0, nrow(popx), length(popgs)) #initialize matrix
        colnames(popgmat) <- popgs #set the matrix names

        #fill the matrix, one possible genotype at a time (hard enough to vectorize as it is).
        for(j in 1:length(popgs)){
          poptab <- count_genos(popgs[j], popxv)
          popgmat[as.numeric(names(poptab)),j] <- as.numeric(poptab)
        }

        if(length(popmpos) > 0){
          poptmat <- popgmat[,-c(popmpos)] #gmat without missing data
        }
        else{
          poptmat <- popgmat
        }

        #get matrix of allele counts
        #initialize
        popamat <- matrix(0, nrow(popgmat), length(popas))
        colnames(popamat) <- popas

        #fill in
        for(j in 1:length(popas)){
          b <- grep(popas[j], colnames(poptmat))
          hom <- which(colnames(poptmat) == paste0(popas[j], popas[j]))
          if(length(hom) == 0){
            het <- b
            popamat[,j] <- rowSums(poptmat[,het])
          }
          else{
            het <- b[b != hom]
            if(length(het) > 0){
              if(is.matrix(poptmat[,het])){
                popamat[,j] <- (poptmat[,hom] * 2) + rowSums(poptmat[,het])
              }
              else{
                popamat[,j] <- (poptmat[,hom] * 2) + poptmat[,het]
              }
            }
            else{
              popamat[,j] <- (poptmat[,hom] * 2)
            }
          }
        }

        #-------------------
        #established new set for this, now do maf as above.
        popmafs <- 1 - rowMaxs( popamat)/rowSums( popamat)
        popmafs <-  popmafs >= maf #if greater than requested, passes check for this pop
        popmafs[is.na(popmafs)] <- FALSE #call false, no information to go on here
        pmafs <- pmafs + popmafs #add this logical to the vector containing the sum of all such vectors
      }
      pmafs <- !as.logical(pmafs) #if false, we keep the allele since it was at >= maf in at least one pop
      keep <- keep + pmafs
    }
  }

  #===============================================
  ##hf_hets. Should only run if bi_al = TRUE.
  if(hf_hets){
    cat("Filtering high frequency heterozygote loci...\n")
    #get heterozygote counts
    if(sum(hs[-mpos]) > 1){
      hetsum <- rowSums(tmat[,hs[-mpos]])
    }
    else if (sum(hs[-mpos]) == 0){
      hetsum <- numeric(nrow(x))
    }
    else{
      hetsum <- tmat[,hs[-mpos]]
    }
    hf <- hetsum/rowSums(tmat)
    hf <- hf >= hf_hets #if false, heterozygote frequency is lower than cut-off, keep locus
    keep <- keep + hf
  }

  ##############################################################################################
  ###apply results
  cat("Applying results...\n")
  keep <- !as.logical(keep)
  x <- x[keep,]
  x <- cbind(headers[keep,], x)
  return(x)
}





#Converts snp/microsat data into different formats. Current supported formats are:
#Arguments:
# data: input data. Rows are loci, columns are individuals.
# ecs: number of extra columns at the header of the input (metadata, ect)
# output: output format:
#   1) per pop allele count format for Bayescan, ect. (v)
#   2) genepop format (v)
#   3) structure format (v)
#   4) 2 character numeric format (v)
#   5) migrate-n hapmap format (v)
#   6) NN format (from numeric) (v)
#   7) allele presence absensence format
# input_form: Input format
#   "NN": Default, alleles as single characters (as in AT or CG),
#   "0000": numeric. Will step to NN first before converting.
#   "msat_n": microsat, in two or three character format.
#           If in two character format, use "msat_2". If in three, use "msat_3"
#           Currently only supported for conversion to 8; 2, 3, and 4 forthcoming.
#   "snp_tab": Converts to NN first.
# miss: missing data format (per allele), typically "N", or "00"
# pop: For format 1 and 5.
#      Number of samples in each pop as a list (as in list(c("POP1", "POP2"), c(48,50)))).
#      Samples must be in the order given in input data frame.
# n_samp: number of randomly selected snps, ect to take. Can also take a numeric vector containing SNP indices to keep. For option 3.
# l_names: Vector of locus names. For option 8.
# interp_miss: Should missing data be interpolated? T or F. For option 8.
#note: (v) under format options means that I've already vectorized it.
format_snps <- function(data, ecs, output = 1, input_form = "NN",
                        miss = "N", pop = 1, n_samp = NA, l_names = NULL,
                        interp_miss = T){

  #do checks, print info
  if(output == 1){
    cat("Output 1 selected, converting to per pop allele count format.\n")
    if(input_form != "0000" & input_form != "NN"){
      stop("Only 0000 and NN formats accepted.")
    }
    if(is.list(pop)){
      cat("Pops:")
      for (i in 1:length(unlist(pop[1]))){
        cat("\n\t", unlist(pop[1])[i], ", size: ", unlist(pop[2])[i])
      }
      cat("\n")
    }
  }
  else if(output == 2){
    cat("Output 2 selected, converting to genepop format.\n")
    if(input_form != "0000" & input_form != "NN" & input_form != "snp_tab"){
      stop("Only 0000, NN, and snp_tab formats accepted for now.")
    }
  }
  else if(output == 3){
    cat("Output 3 selected, converting to STRUCTURE format.\n")
    if(input_form != "0000" & input_form != "NN"){
      stop("Only 0000 and NN formats accepted for now.")
    }
    if(!is.integer(n_samp)){
      cat("Number of sub-samples to take:", "n_samp", "\n")
    }
    else if (!is.na(n_samp)){
      stop("Number of sub-samples to take must be an integer.\n")
    }
  }
  else if(output == 4){
    cat("Output 4 selected, converting to numeric 2 character format.\n")
    if(input_form != "NN" & input_form != "snp_tab"){
      stop("Only NN and snp_tab formats accepted for now.")
    }
  }
  else if(output == 5){
    cat("Output 5 selected, converting to migrate-N hapmap format.\n")
    if(input_form != "NN" & input_form != "0000"){
      stop("Only 0000 and NN formats accepted.")
    }
    if(is.list(pop)){
      cat("Pops:")
      for (i in 1:length(unlist(pop[1]))){
        cat("\n\t", unlist(pop[1])[i], ", size: ", unlist(pop[2])[i])
      }
      cat("\n")
    }
  }
  else if(output == 6){
    cat("Output 6 selected, converting to NN format.\n")
    if(input_form != "0000" & input_form != "snp_tab"){
      stop("Only 0000 and snp_tab formats accepted.")
    }
  }

  else if(output == 7){
    cat("Output 7 selected, converting to allele presence/absense format.\n")
  }
  else{
    stop("Please specify output format.")
  }

  if(input_form == "NN"){
    cat("Input format: NN\n")
    if(nchar(miss) != 1){
      stop("Missing data format must be one character.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "0000"){
    cat("Input format: 0000\n")
    if(nchar(miss) != 2){
      stop("Missing data format must be two characters.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "msat_2"){
    cat("Input format: msat, 2 character (0000)\n")
    if(nchar(miss) != 2){
      stop("Missing data format must be two characters.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "msat_3"){
    cat("Input format: msat, 3 character (000000)\n")
    if(nchar(miss) != 3){
      stop("Missing data format must be three characters.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else if(input_form == "snp_tab"){
    cat("Input format, snp_tab.\n")
    if(nchar(miss) != 1){
      stop("Missing data format must be one character.\n")
    }
    else{
      cat("Missing data format: ", miss, "\n")
    }
  }
  else{
    stop("Unsupported input format.")
  }


  #####################################################
  #function to create table of allele and genotype counts from all data. From filter_snps
  #input is a long vector of the data xv,
  #snp_form, the length of each genotype.
  #mDat, the missing data format
  #nsap, the number of samples
  tabulate_genotypes <- function(xv, snp_form, mDat, nsamp){
    #create tables of genotype counts for each locus.
    #function to do this, takes the pattern to match and the data, which is as a SINGLE VECTOR, by rows.
    #needs the number of samples from global function environment
    count_genos <- function(pattern, x){
      XXs <- grep(pattern, x) #figure out which entries match pattern
      out <- table(ceiling(XXs/nsamp)) #devide the entries that match by the number of samps to get which samp they come from (when rounded up), then table that.
      return(out)
    }

    #get all possible genotypes
    gs <- unique(xv)

    #which gs are heterozygous?
    hs <- substr(gs,1,snp_form/2) != substr(gs, (snp_form/2 + 1), snp_form*2)

    #which genotype is the missing data?
    mpos <- which(gs == mDat)

    #what are the possible alleles at all loci?
    as <- unique(c(substr(gs,1,snp_form/2), substr(gs, (snp_form/2 + 1), snp_form*2)))
    as <- as[as != substr(mDat, 1, snp_form/2)] #that aren't N?

    #############################################################################################3
    ###get a table of genotype and allele counts at each locus.
    cat("Creating genotype table...\n")

    #for each element of gs, get the tables of genotype counts and add them to a matrix
    gmat <- matrix(0, nrow(data), length(gs)) #initialize matrix
    colnames(gmat) <- gs #set the matrix names
    #fill the matrix, one possible genotype at a time (hard enough to vectorize as it is).
    for(i in 1:length(gs)){
      tab <- count_genos(gs[i], xv)
      gmat[as.numeric(names(tab)),i] <- as.numeric(tab)
    }

    if(length(mpos) > 0){
      tmat <- gmat[,-c(mpos)] #gmat without missing data
    }
    else{
      tmat <- gmat
    }

    #get matrix of allele counts
    #initialize
    cat("Getting allele table...\n")
    amat <- matrix(0, nrow(gmat), length(as))
    colnames(amat) <- as

    #fill in
    for(i in 1:length(as)){
      b <- grep(as[i], colnames(tmat))
      hom <- which(colnames(tmat) == paste0(as[i], as[i]))
      if(length(hom) == 0){
        het <- b
        amat[,i] <- rowSums(tmat[,het])
      }
      else{
        het <- b[b != hom]
        if(length(het) > 0){
          if(is.matrix(tmat[,het])){
            amat[,i] <- (tmat[,hom] * 2) + rowSums(tmat[,het])
          }
          else{
            amat[,i] <- (tmat[,hom] * 2) + tmat[,het]
          }
        }
        else{
          amat[,i] <- (tmat[,hom] * 2)
        }
      }
    }
    return(list(gs = tmat, as = amat))
  }

  #convert 0000 to NN unless output is 7 or 2. Return after if output is 6. (v)
  if(input_form == "0000" & output != 7 & output != 2){

    #Do conversion
    cat("Converting genotypes to NN form intermediary...")

    #vectorize and replace
    xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
    xv1 <- substr(xv, 1, 2)
    xv2 <- substr(xv, 3, 4)
    xv1[xv1 == "01"] <- "A"
    xv1[xv1 == "02"] <- "C"
    xv1[xv1 == "03"] <- "G"
    xv1[xv1 == "04"] <- "T"
    xv1[xv1 == miss] <- "N"
    xv2[xv2 == "01"] <- "A"
    xv2[xv2 == "02"] <- "C"
    xv2[xv2 == "03"] <- "G"
    xv2[xv2 == "04"] <- "T"
    xv2[xv2 == miss] <- "N"
    xv <- paste0(xv1, xv2)

    #rebind to matrix and remake data.
    xv <- matrix(xv, nrow(data), (ncol(data) - ecs), T)
    data <- cbind(data[,1:ecs], as.data.frame(xv))

    cat("\nMoving on to conversion...", "\n")
    miss <- "N" #reset miss to the correct entry
    if(output == 6){
      return(data) #all done if just converting to NN
    }
  }

  else if (input_form == "NN"){
    cat("Ensure that all data columns are character vectors!\n")
  }

  #convert snp_tab to NN (v)
  if(input_form == "snp_tab"){
    header <- data[,1:ecs]
    xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
    nsamp <- ncol(data) - ecs
    snames <- colnames(data)[(ecs+1):ncol(data)]
    rm(data)
    ptf <- nchar(xv) == 1
    xv[ptf] <- paste0(xv[ptf], xv[ptf]) #double up homozygotes
    remove(ptf)
    xv <- gsub(" ", "", xv) #combine heterozygotes.
    xv[xv == paste0(miss, miss)] <- "NN" #replace with correct missing data
    xv <- as.data.frame(matrix(xv, nrow(header), nsamp, byrow = T), stringsAsFactors = F)
    colnames(xv) <- snames
    out <- cbind(header, xv)
    if(output == 6){
      return(out)
    }
    else{
      data <- out
    }
  }


  ##convert to allele count or migrate-n format, migrate-n should ALWAYS have multiple pops (why else would you
  ##use it?) (v)
  if(output == 1 | output == 5){
    if(output == 5){cat("WARNING: Data does not have header or pop spacer rows.\n")}
    w_df <- data[,1:ecs] #intiallize w_df if doing option 1

    ##Make a function to sum each row
    sum_row <- function(row){
      l <- c(substr(row,1,1), substr(row,2,2))
      t <- table(l)
      return(t) #returns a table with counts of the two alleles and then an N
    }
    if (is.list(pop) == FALSE){
      if(output == 5){stop("Cannot convert to migrate-n format with only one pop.\n")}

      #prepare a genotype table
      xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
      tabs <- tabulate_genotypes(xv, nchar(xv[1]), paste0(miss, miss), ncol(data) - ecs)
      w_df$n_total <- rowSums(tabs$as)
      if(any(w_df$n_total != 2)){stop("More or less than two alleles detected at some loci.\n")}
      w_df$n_alleles <- rowSums(tabs$as != 0)
      #fill in ni1 and ni2.
      w_df$ni1 <- tabs$as[,"A"]
      w_df$ni1[w_df$ni1 == 0] <- tabs$as[which(w_df$ni1== 0), "G"]
      w_df$ni1[w_df$ni1 == 0] <- tabs$as[which(w_df$ni1== 0), "C"]
      w_df$ni2 <- tabs$as[,"T"]
      w_df$ni2[w_df$ni2 == 0] <- tabs$as[which(w_df$ni2== 0), "C"]
      w_df$ni2[w_df$ni2 == 0] <- tabs$as[which(w_df$ni2== 0), "G"]
    }
    else{
      pop_count <- length(unlist(pop[1]))
      pop_sizes <- unlist(pop[2])
      if(sum(unlist(pop[2])) != (ncol(data) - ecs)){#checks to make sure that the pop sizes add up to agree
        #with data
        cat("Supplied population numbers do not equal the number of supplied loci columns. Exiting.")
        stop()
      }
      #initialize w_df for multiple pops
      #print(unlist(pop[1]))
      j <- 2
      w_df$pop <- unlist(pop[1])[1] #create pop column and set first pop name
      wa_df <- w_df #set first temp w_df
      for (j in j:pop_count){ #loop through each pop name
        wb_df <- wa_df #create temp df copy of wa_df
        #print(j)
        #print(unlist(pop[1])[j])
        wb_df$pop <- unlist(pop[1])[j]#change pop to current name
        w_df <- rbind(w_df, wb_df) #bind this copy to w_df
      }
      remove(wa_df, wb_df) #remove temp df

      ##loop through and count alleles for every row, splitting by pop

      #create allele count table for each pop, fill their section of data

      #build allele tables for each locus
      pop_as <- list()
      pals <- matrix(FALSE, nrow(data), 4)
      colnames(pals) <- c("A", "C", "G", "T")
      current <- ecs + 1
      for(i in 1:pop_count){
        #get the data for this population
        cat("Generating population allele count matrices, current pop:", pop[[1]][i], "\n")
        ne <- current + pop_sizes[i] - 1
        x <- data[,current:ne]
        xv <- as.vector(t(x))
        pop_as[[i]] <- tabulate_genotypes(xv, nchar(xv[1]), paste0(miss, miss), pop_sizes[i])$as
        pop_as[[i]] <- pop_as[[i]][, order(colnames(pop_as[[i]]))]
        #figure out which alleles are present at each locus to write these.
        pals <- pals + ifelse(pop_as[[i]] != 0, TRUE, FALSE)
        current <- current + pop_sizes[i]

      }

      #now write the correct counts
      cat("Tabulating and writing results...\n")
      pals <- ifelse(pals != 0, TRUE, FALSE)
      if(any(rowSums(pals) != 2)){
        stop("More or less than two alleles detected at some loci! Check/filter input data.\n")
      }
      #Convert into vector which says which elements to keep.
      palsv <- as.vector(t(pals))

      w_df$n_total <- NA
      w_df$n_alleles <- NA
      w_df$ni1 <- NA
      w_df$ni2 <- NA

      #loop through and print data
      for(i in 1:pop_count){
        #compare palsv to allele tables to keep correct alleles.
        av <- as.vector(t(pop_as[[i]]))
        av <- av[palsv]
        av <- matrix(av, nrow(data), 2, byrow = T)
        w_df[w_df$pop == pop[[1]][i],]$n_total <- rowSums(av)
        w_df[w_df$pop == pop[[1]][i],]$n_alleles <- rowSums(ifelse(av != 0, TRUE, FALSE))
        w_df[w_df$pop == pop[[1]][i],]$ni1 <- av[,1]
        w_df[w_df$pop == pop[[1]][i],]$ni2 <- av[,2]
      }


      row.names(w_df) <- 1:nrow(w_df) #fix row names
      if(output == 5){ #if doing output style 5...
        w_df <- w_df[,c(1, (ecs + 1):ncol(w_df))] #remove extra columns except the first and pop columns
      }
    }
    return(w_df)
  }



  ##convert to genepop or numeric format (v)
  if (output == 2 | output == 4){
    cat("WARNING: For output option 3, output does not have population seperators or header information.", "\n", "Converting genotypes...")

    if(input_form == "0000"){
      if(output == 2){
        w_df <- as.data.frame(t(data[,(ecs + 1):ncol(data)])) #remove extra columns and transpose data
        #print(names(w_df))
        row.names(w_df) <- paste0(row.names(w_df), " ,") #adding space and comma to row names, as required.
        return(w_df)
      }
      else{
        stop("Please select different input and output formats.")
      }
    }

    #vectorize and replace
    xv <- as.vector(t(data[,(ecs + 1):ncol(data)]))
    xv1 <- substr(xv, 1, 1)
    xv2 <- substr(xv, 2, 2)
    xv1[xv1 == "A"] <- "01"
    xv1[xv1 == "C"] <- "02"
    xv1[xv1 == "G"] <- "03"
    xv1[xv1 == "T"] <- "04"
    xv1[xv1 == miss] <- "00"
    xv2[xv2 == "A"] <- "01"
    xv2[xv2 == "C"] <- "02"
    xv2[xv2 == "G"] <- "03"
    xv2[xv2 == "T"] <- "04"
    xv2[xv2 == miss] <- "00"
    xv <- paste0(xv1, xv2)

    xv <- matrix(xv, nrow(data), (ncol(data) - ecs), T)
    data <- cbind(data[,1:ecs], as.data.frame(xv))

    cat("\n", "Cleaning up...", "\n")
    if(output == 2){ #convert to genepop
      w_df <- as.data.frame(t(data[,(ecs + 1):ncol(data)])) #remove extra columns and transpose data
      #print(names(w_df))
      row.names(w_df) <- paste0(row.names(w_df), " ,") #adding space and comma to row names, as required.
      return(w_df)
    }
    else {#prepare numeric output, otherwise same format
      return(data)
    }
  }


  ##convert to structure format (v)
  if (output == 3){

    #subset if requested
    if(!is.na(n_samp)){
      cat("Subsetting ")
      if(length(n_samp) > 1){
        cat("designated SNPs.\n")
        data <- data[,n_samp]
      }
      else{
        cat(n_samp, " random SNPs.\n")
        data <- data[sample(nrow(data), n_samp, T),]
      }
    }

    #plop actual data into two vectors, one for each allele.
    xv <- as.vector(t(t(data[,(ecs + 1):ncol(data)])))
    xv1 <- substr(xv, 1, 1)
    xv2 <- substr(xv, 2, 2)
    remove(xv)

    #convert to numeric
    xv1[xv1 == "A"] <- 1
    xv1[xv1 == "C"] <- 2
    xv1[xv1 == "G"] <- 3
    xv1[xv1 == "T"] <- 4
    xv1[xv1 == miss] <- 0
    xv2[xv2 == "A"] <- 1
    xv2[xv2 == "C"] <- 2
    xv2[xv2 == "G"] <- 3
    xv2[xv2 == "T"] <- 4
    xv2[xv2 == miss] <- 0

    #bind back into matrices
    xv1 <- matrix(xv1, ncol(data) - ecs, nrow(data), T)
    xv2 <- matrix(xv2, ncol(data) - ecs, nrow(data), T)

    #create output matrix
    outm <- matrix(NA, 2*(ncol(data) - ecs), nrow(data))

    #fill
    outm[seq(1,nrow(outm),2),] <- xv1
    outm[seq(2,nrow(outm),2),] <- xv2

    #add sample names
    snames <- character(nrow(outm))
    snames[seq(1,nrow(outm),2)] <- colnames(data[,(ecs+1):ncol(data)])
    snames[seq(2,nrow(outm),2)] <- colnames(data[,(ecs+1):ncol(data)])
    out <- cbind(ind = snames, as.data.frame(outm))

    warning("Remove header line before inputing data into structure!")
    return(out)
  }

  #presence/absence format
  if(output == 7){
    x <- data[,(ecs+1):ncol(data)] #get just data
    if((is.null(l_names)) | length(l_names) != nrow(data) | !(is.vector(l_names))){
      stop("Locus names required for output format 8. Please give a vector of names (l_names) equal to number of loci.")
    }
    if(input_form == "NN"){
      asize <- 1
    }
    else if (input_form == "0000" | input_form == "msat_2"){
      asize <- 2
    }
    else if (input_form == "msat_3"){
      asize <- 3
    }

    #write a function to quickly get allele list per locus
    acount <- function(x){
      a <- unique(c(substr(x, 1, asize), substr(x, (asize + 1), (2*asize))))
      return(sort(a))
    }

    #function to get average allele frequencies if needed.
    afreqs <- function(x){
      tab <- table(c(substr(x, 1, asize), substr(x, (asize + 1), (2*asize))))
      tab <- tab[names(tab) != miss]
      tab <- tab/sum(tab)
      tab <- tab[sort(names(tab))]
      return(tab)
    }

    #loop through each row, figure out number of alleles, create and fill matrix with presence/absence, bind matrix to existing.
    omat <- matrix(NA, ncol(x), 1)
    for (i in 1:nrow(x)){
      if(i %% 1000 == 0){cat("snp number: ", i, "\n")}
      alist <- acount(x[i,])
      wmat <- matrix(0, ncol(x), length(alist)) #output matrix. Individuals in rows, alleles in columns.
      colnames(wmat) <- paste0(l_names[i], "_", alist)
      a1s <- substr(x[i,], 1, asize) #get first allele in each ind
      a2s <- substr(x[i,], asize + 1, 2*asize) #get second allele in each ind.
      a1s <- match(a1s, alist) #get column indexes to add a presence to
      a2s <- match(a2s, alist) #same
      for(j in 1:nrow(wmat)){ #for each individual, add alleles to the appropriate bins
        wmat[j, a1s[j]] <- wmat[j, a1s[j]] + 1
        wmat[j, a2s[j]] <- wmat[j, a2s[j]] + 1
      }
      if(interp_miss){ #interpolate missing genotypes if requested
        afs <- afreqs(x[i,]) #get overall allele frequencies
        ms <- which(wmat[,grepl(paste0("_", miss), colnames(wmat))] > 0) #figure out which had missing data
        wmat <- wmat[,-grepl(paste0("_", miss), colnames(wmat))] #remove the missing data column
        wmat[ms,] <- matrix(afs, length(ms), ncol(wmat), byrow = T) #add the interpolated allele frequencies.
      }

      omat <- cbind(omat, wmat) #bind to the output
    }
    w_df <- as.data.frame(omat[,-1]) #remove int column, convert to data frame.
    w_df <- cbind(samp = colnames(x), w_df) #add sample names.
    if(!interp_miss){cat("Finished. Warning: Missing data counts are also stored!\n")}
    return(w_df)
  }
}

#function to run any command after spliting the data by group and by population. Group and population must be in columns named "group" and "pop", respectively.
#The first argument of the function to run must be the data provided to that function.
run_gp <- function(x, FUN, ...){
  w_df<- data.frame()
  pops <- unique(x[,"pop"])
  groups <- unique(x[,"group"])
  for (i in 1:length(groups)){
    w_data <- x[x[,"group"] == groups[i],]
    print(groups[i])
    for(j in 1:length(pops)){
      print(pops[j])
      x_data <- FUN(w_data[w_data[,"pop"] == pops[j],], ...)
      if(length(x_data) == 0){next}
      if(is.data.frame(x_data)){
        x_data <- cbind(pop = pops[j], x_data)
        w_df <- rbind(w_df, x_data)
      }
      else{
        if(is.data.frame(w_df)){
          w_df <- numeric(0)
        }
        w_df <- c(w_df,x_data)
      }
    }
  }
  if(!is.data.frame(w_df)){
    warning("Output is vector: if merging with meta info, check that info data frame is sorted identically:
            Group, Pop, then snp by snp!")
  }
  return(w_df)
}

#function to run any command after spliting the data by group. Group must be in a column named "pop".
#The first argument of the function to run must be the data provided to that function.
run_g <- function(x, FUN, ...){
  w_df<- data.frame()
  groups <- unique(x[,"group"])
  for (i in 1:length(groups)){
    w_data <- x[x[,"group"] == groups[i],]
    print(groups[i])
    x_data <- FUN(w_data[w_data[,"pop"] == pops[j],], ...)
    if(length(x_data) == 0){next}
    if(is.data.frame(x_data)){
      w_df <- rbind(w_df, x_data)
    }
    else{
      if(is.data.frame(w_df)){
        w_df <- numeric(0)
      }
      w_df <- c(w_df, x_data)
    }
  }
  if(!is.data.frame(w_df)){
    warning("Output is vector: if merging with meta info, check that info data frame is sorted identically:
            Group, Pop, then snp by snp!")
  }
  return(w_df)
}
