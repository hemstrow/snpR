#Genotype individuals for genomic inversion spanning multiple SNPs as single psuedo SNP. Requires a list of
#SNPs on the inversion, with a column named "snp" that identifies snps, as well as reference data
#from which the major and minor alleles (fixed differences between the inversion and the
#ancestral genome) will be called. These can be the same file, and can also be the same as the data to be
#genotyped, but at least the minor reference must contain at least one instance of an individual with
#rarer genotype (either a copy of the inversion or ancestral genome). Out = 1 returns the inversion (or rare
#ancestral) genotype for each individual as a psuedo SNP, out = 2 calls each snp genotype for each individual as
#either Maj or Minor (anc or inv), out = 3 returns a list of counts of each MinMaj genotype across all SNPs for each individual.
#Header and Footer col counts refer to the number of extra columns at either the start or end of the data,
#including the required SNP names. The calling cuttoff is the percent (in decimals) of SNPs which must be
#in agreement for the inversion to be genotyped. Min nonzeroes denotes the number of ungenotyped snps allowed
#for the genotype of the inversion to be called (note that if the SNPs provided in the SNP list are truely
#fixed and there is no genotyping error, one SNP is all that is neccisary. Multiple are likely more accurate,
#however.).
genotype_inversion <- function(data, SNP_list, Min_ref, Maj_ref, header_col_count, calling_cuttoff, min_non_zeros, footer_col_count, out = 1){
  i <- 1
  data <- data[,1:ncol(data)-footer_col_count]
  Min_ref <- Min_ref[,1:ncol(Min_ref)-footer_col_count]
  Maj_ref <- Maj_ref[,1:ncol(Maj_ref)-footer_col_count]
  div_snps <- subset(data, data$snp %in% SNP_list)
  div_snps <- arrange(div_snps, snp)
  Min_ref <- subset(Min_ref, Min_ref$snp %in% SNP_list)
  Maj_ref <- subset(Maj_ref, Maj_ref$snp %in% SNP_list)
  Min_ref <- arrange(Min_ref, snp)
  Maj_ref <- arrange(Maj_ref, snp)
  #print(head(Min_ref))
  #print(Maj_ref)
  div_snps$group <- as.character(div_snps$group)
  ind_state <- data.frame(matrix(ncol = ncol(div_snps), nrow = nrow(div_snps)))
  #print(ind_state)
  colnames(ind_state) <- colnames(div_snps)
  #print(names(div_snps))

  #Create list of major and minor alleles from the list of fixed snps and two (or one) sets of reference samples.
  FixInvMin <- apply(Min_ref[,(header_col_count + 1):length(Min_ref)],1,function(y)
    names(which.min(table(subset(c(substr(y, 1,2), substr(y,3,4)), c(substr(y, 1,2), substr(y,3,4))[] != "00")))))
  #print(FixInvMin)
  FixInvMaj <- apply(Maj_ref[,(header_col_count + 1):length(Maj_ref)],1,function(x)
    names(which.max(table(subset(c(substr(x, 1,2), substr(x,3,4)), c(substr(x, 1,2), substr(x,3,4))[] != "00")))))
  #print(FixInvMaj)
  InvMinMaj <- data.frame(Maj = FixInvMaj, Min = FixInvMin)
  remove(FixInvMin)
  remove(FixInvMaj)
  #print(InvMinMaj)

  #Call each genotype according to Major/Minor allele identities (01 and 02, respectively), write to ind_state.
  for(i in i:nrow(div_snps)){
    #print(div_snps[i,1:header_col_count])
    j <- header_col_count + 1
    ind_state[i,1:header_col_count] <- div_snps[i,1:header_col_count]
    for(j in j:ncol(div_snps)){
      a1 <- substr(div_snps[i,j],1,2)
      a2 <- substr(div_snps[i,j],3,4)
      if(a1 == a2){
        if(a1 == InvMinMaj[i,"Maj"]){
          state <- "0101"
        }
        else if (a1 == InvMinMaj[i,"Min"]){
          state <- "0202"
        }
        else if (a1 == "00"){
          state <- "0000"
        }
        else{
          cat("Undocumented allele at snp:", div_snps[i,1], ",individual:", j - header_col_count, "!\n")
          state <- "other"
        }
      }
      else{
        if ((a1 == InvMinMaj[i,"Maj"] | a1 == InvMinMaj[i,"Min"]) & (a2 == InvMinMaj[i,"Maj"] | a2 == InvMinMaj[i, "Min"])){
          state <- "0102"
        }
        else{
          cat("Undocumented allele at snp:", div_snps[i,1], ",individual:", j - header_col_count, "!\n")
          state <- "other"
        }
      }
      ind_state[i,j] <- state
    }
  }
  if(out == 2){
    return(ind_state)
    stop
  }

  #Call individual genotype for whole inversion
  if (out == 1){
    inv_state <- data.frame(individual = numeric(0), genotype = numeric(0), percent_snp_cons = numeric(0))
  }
  i <- 1
  j <- header_col_count + 1
  if(out == 3){
    c_df <- data.frame(individual = numeric(0), g0000 = numeric(0), g0101 = numeric(0), g0102 = numeric(0), g0202 = numeric(0), other = numeric(0))
  }
  for (j in j:ncol(ind_state)){
    i <- 1
    if(out == 1){
      inv_state[j - header_col_count,"individual"] <- colnames(ind_state)[j]
    }
    else if (out == 3){
      c_df[j - header_col_count,"individual"] <- colnames(ind_state)[j]
    }
    cons0000 <- 0
    cons0101 <- 0
    cons0102 <- 0
    cons0202 <- 0
    consother <- 0
    for(i in i:nrow(ind_state)){
      if(ind_state[i,j] == "0000"){
        cons0000<- cons0000 + 1
      }
      else if(ind_state[i,j] == "0101"){
        cons0101 <- cons0101 + 1
      }
      else if(ind_state[i,j] == "0102"){
        cons0102 <- cons0102 + 1
      }
      else if(ind_state[i,j] == "0202"){
        cons0202 <- cons0202 + 1
      }
      else if(ind_state[i,j] == "other"){
        consother <- consother + 1
      }
      else{
        cat("Error in inv_state calling, ind_state col/row:", i, "/", j, ".")
        stop
      }
    }
    tot <- cons0202 + cons0101 + cons0102 + consother
    #cat("Ind: ", colnames(ind_state)[j], "\t0000:", cons0000, "\t0101: ", cons0101, "\t0102: ", cons0102, "\t0202: ", cons0202, "\tother:", consother, "\n\n")
    if(out == 3){
      c_df[j - header_col_count, "g0000"] <- cons0000
      c_df[j - header_col_count, "g0101"] <- cons0101
      c_df[j - header_col_count, "g0102"] <- cons0102
      c_df[j - header_col_count, "g0202"] <- cons0202
      c_df[j - header_col_count, "other"] <- consother
    }
    else if (out == 1){
      if (tot < min_non_zeros){
        inv_state[j - header_col_count, "genotype"] <- "0000"
        inv_state[j - header_col_count, "percent_snp_cons"] <- cons0000/(tot + cons0000)
      }
      else if (cons0101/tot >= calling_cuttoff){
        inv_state[j - header_col_count, "genotype"] <- "0101"
        inv_state[j - header_col_count, "percent_snp_cons"] <- cons0101/tot
      }
      else if (cons0102/tot >= calling_cuttoff){
        inv_state[j - header_col_count, "genotype"] <- "0102"
        inv_state[j - header_col_count, "percent_snp_cons"] <- cons0102/tot
      }
      else if (cons0202/tot >= calling_cuttoff){
        inv_state[j - header_col_count, "genotype"] <- "0202"
        inv_state[j - header_col_count, "percent_snp_cons"] <- cons0202/tot
      }
      else if (consother/tot >= calling_cuttoff){
        inv_state[j - header_col_count, "genotype"] <- "other"
        inv_state[j - header_col_count, "percent_snp_cons"] <- consother/tot
      }
      else {
        inv_state[j - header_col_count, "genotype"] <- "0000"
        inv_state[j - header_col_count, "percent_snp_cons"] <- max(c(consother, cons0202, cons0102, cons0102))/tot
      }
    }
  }
  if(out == 1){
    return(inv_state)
  }
  else if (out == 3){
    return(c_df)
  }
}
