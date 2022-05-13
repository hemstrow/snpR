# This contains functions that are currently depreciated. Some may be restored in the future.

# LD from hapmap and ms. Only a snippet, no supporting. Check the old versions of LD full pairwise.
LD_hap_ms <- function(x, mDat, input){
  if(input == "haplotype" | input == "ms"){
    library(data.table)
    #functions
    tabulate_haplotypes <- function(x){
      thaps <- matrix(paste0(x[1,], t(x[-1,])), ncol = nrow(x) - 1) #convert each cell to haplotype vs row one
      mthaps <- reshape2::melt(thaps) #put this into long form for tabulation
      cnames <- levels(mthaps$value)
      hapmat <- bigtabulate::bigtabulate(mthaps, ccols = c(2,3))
      colnames(hapmat) <- cnames
      return(hapmat)
    }
    LD_func <- function(x, meta, prox_table = TRUE, matrix_out = TRUE, mDat = "N", sr = FALSE, chr.length = NULL, stop.row = nrow(x) - 1){
      #function to count the number of haplotypes vs every other given position
      #################
      #prep stuff
      if(prox_table){
        ncomps <- (nrow(x)*(nrow(x)-1)/2) - (nrow(x)-stop.row)*((nrow(x)-stop.row)-1)/2
        prox <- data.table::as.data.table(data.frame(p1 = numeric(ncomps),
                                                     p2 = numeric(ncomps),
                                                     rsq = numeric(ncomps),
                                                     Dprime = numeric(ncomps),
                                                     pval = numeric(ncomps)))
      }
      if(matrix_out){
        rmat <- matrix(as.numeric(NA), stop.row, nrow(x) - 1)
        rmat <- data.table::as.data.table(rmat)
        colnames(rmat) <- as.character(meta$position[-1])
        rownames(rmat) <- make.names(as.character(meta$position[1:stop.row]), unique = T)
        Dpmat <- data.table::copy(rmat)
        pvmat <- data.table::copy(rmat)
      }
      if(!matrix_out & !prox_table){
        stop("Please specify output format.\n")
      }

      #run length prediction variables for progress reporting
      compfun <- function(x){
        return(((x-1)*x)/2)
      }
      totcomp <- compfun(nrow(x))
      prog <- 0
      cpercent <- 0

      ##################
      #calculate Dprime, rsq, ect.
      for(i in 1:stop.row){
        prog_after <- prog + nrow(x) - i

        #report progress
        if(!sr){
          cprog <- (totcomp - compfun(nrow(x) - i - 1))/totcomp
          if(cprog >= 0.01 + cpercent){
            cat("Progress:", paste0(round(cprog*100), "%."), "\n")
            cpercent <- cprog
          }
        }

        #check that this site isn't fixed.
        if(length(unique(x[i,])) == 1){
          if(prox_table){
            data.table::set(prox, (prog + 1):prog_after, j = "p1", value = meta$position[i])
            data.table::set(prox, (prog + 1):prog_after, j = "p2", value = meta$position[(i+1):nrow(x)])
            data.table::set(prox, (prog + 1):prog_after, j = "rsq", value = NA)
            data.table::set(prox, (prog + 1):prog_after, j = "Dprime", value = NA)
            data.table::set(prox, (prog + 1):prog_after, j = "pval", value = NA)
          }
          prog <- prog_after
          next()
        }

        #get haplotypes
        haps <- tabulate_haplotypes(x[i:nrow(x),])

        #fix for very rare cases
        if(!is.matrix(haps)){
          haps <- t(as.matrix(haps))
        }

        #Fix only three haplotypes. While Dprime is 1, can't just set rsq, so fix the table and continue.
        if(ncol(haps) == 3){
          pos.a <- unique(unlist(unlist(strsplit(colnames(haps), ""))))
          if(!(any(colnames(haps) == paste0(pos.a[1], pos.a[1])))){
            tnames <- colnames(haps)
            if(nrow(haps) == 1){
              haps <- t(as.matrix(c(0, haps)))
            }
            else{
              haps <- cbind(numeric(nrow(haps)), haps)
            }
            colnames(haps)<- c(paste0(pos.a[1], pos.a[1]), tnames)
          }
          else if(!(any(colnames(haps) == paste0(pos.a[1], pos.a[2])))){
            tnames <- colnames(haps)
            if(nrow(haps) == 1){
              haps <- t(as.matrix(c(haps[1], 0, haps[2:3])))
            }
            else{
              haps <- cbind(haps[,1], numeric(nrow(haps)), haps[,2:3])
            }
            colnames(haps) <- c(tnames[1], paste0(pos.a[1], pos.a[2]), tnames[2:3])
          }
          else if(!(any(colnames(haps) == paste0(pos.a[2], pos.a[1])))){
            tnames <- colnames(haps)
            if(nrow(haps) == 1){
              haps <- t(as.matrix(c(haps[1:2], 0, haps[3])))
            }
            else{
              haps <- cbind(haps[,1:2], numeric(nrow(haps)), haps[,3])
            }
            colnames(haps) <-  c(tnames[1:2], paste0(pos.a[2], pos.a[1]), tnames[3])
          }
          else{
            tnames <- colnames(haps)
            if(nrow(haps) == 1){
              haps <- t(as.matrix(c(haps, 0)))
            }
            else{
              haps <- cbind(haps, numeric(nrow(haps)))
            }
            colnames(haps) <- c(tnames, paste0(pos.a[2], pos.a[2]))
          }
        }

        #fix two haplotype cases (NA if either snp is fixed, otherwise 1)
        if(ncol(haps) == 2){
          cn <- colnames(haps)
          check <- ifelse(substr(cn[1], 1, nchar(cn[1])/2) !=
                            substr(cn[2], 1, nchar(cn[1])/2) &
                            substr(cn[1], (nchar(cn[1])/2) + 1, nchar(cn[1])) !=
                            substr(cn[2], (nchar(cn[1])/2) + 1, nchar(cn[1])),
                          1,NA)

          if(is.null(nrow(cn))){
            Dprime <- rep(check, 1)
            rsq <- rep(check, 1)
          }
          else{
            Dprime <- rep(check, nrow(cn))
            rsq <- rep(check, nrow(cn))
          }
          chisqu <- rsq*(4)
          pval <- 1 - pchisq(chisqu, 1)

          #write:
          if(prox_table){
            data.table::set(prox, (prog + 1):prog_after, j = "p1", value = meta$position[i])
            data.table::set(prox, (prog + 1):prog_after, j = "p2", value = meta$position[(i+1):nrow(x)])
            data.table::set(prox, (prog + 1):prog_after, j = "rsq", value = rsq)
            data.table::set(prox, (prog + 1):prog_after, j = "Dprime", value = Dprime)
            data.table::set(prox, (prog + 1):prog_after, j = "pval", value = pval)
          }
          if(matrix_out){
            #reminder: columns start at locus two, rows start at locus one (but end at nlocus - 1)
            fill <- rep(NA, nrow(x) - length(Dprime) - 1)
            Dprime <- c(fill, Dprime)
            rsq <- c(fill, rsq)
            pval <- c(fill, pval)
            data.table::set(rmat, i, 1:ncol(rmat), as.list(rsq))
            data.table::set(Dpmat, i, 1:ncol(Dpmat), as.list(Dprime))
            data.table::set(pvmat, i, 1:ncol(pvmat), as.list(pval))
            rm(pval, rsq, Dprime)
          }
          prog <- prog_after
          next()
        }

        #fix three missing haplotypes (everything NA)
        if(ncol(haps) == 1){

          Dprime <- rep(NA, nrow(haps))
          rsq <- rep(NA, nrow(haps))
          pval <- rep(NA, nrow(haps))

          #write
          if(prox_table){
            data.table::set(prox, (prog + 1):prog_after, j = "p1", value = meta$position[i])
            data.table::set(prox, (prog + 1):prog_after, j = "p2", value = meta$position[(i+1):nrow(x)])
            data.table::set(prox, (prog + 1):prog_after, j = "rsq", value = rsq)
            data.table::set(prox, (prog + 1):prog_after, j = "Dprime", value = Dprime)
            data.table::set(prox, (prog + 1):prog_after, j = "pval", value = pval)
          }
          if(matrix_out){
            #reminder: columns start at locus two, rows start at locus one (but end at nlocus - 1)
            fill <- rep(NA, nrow(x) - length(Dprime) - 1)
            Dprime <- c(fill, Dprime)
            rsq <- c(fill, rsq)
            pval <- c(fill, pval)
            data.table::set(rmat, i, 1:ncol(rmat), as.list(rsq))
            data.table::set(Dpmat, i, 1:ncol(Dpmat), as.list(Dprime))
            data.table::set(pvmat, i, 1:ncol(pvmat), as.list(pval))
            rm(pval, rsq, Dprime)
          }
          prog <- prog_after
          next()
        }

        #calc Dprime, rsq
        A1B1f <- haps[,1]/rowSums(haps)
        A1B2f <- haps[,2]/rowSums(haps)
        A2B1f <- haps[,3]/rowSums(haps)
        A2B2f <- haps[,4]/rowSums(haps)
        A1f <- A1B1f + A1B2f
        A2f <- 1 - A1f
        B1f <- A1B1f + A2B1f
        B2f <- 1 - B1f
        D <- A1B1f - A1f*B1f
        D2 <- A2B1f - A2f*B1f
        Dprime <- ifelse(D > 0, D/matrixStats::rowMins(cbind(A1f*B2f, A2f*B1f)),
                         ifelse(D < 0, D/matrixStats::rowMaxs(cbind(-A1f*B1f, -A2f*B2f)),
                                0))
        rsq <- (D^2)/(A1f*A2f*B1f*B2f)

        #fix for when more missing haps.
        Dprime[which(rowSums(ifelse(haps == 0, T, F)) == 3)] <- 0 #if three missing haplotypes, call 0.
        Dprime[which(rowSums(ifelse(haps == 0, T, F)) == 4)] <- NA # if four, call NA.
        rsq[which(rowSums(ifelse(haps == 0, T, F)) %in% 3:4)] <- NA #replace rsq with NA when 3 or 4 missing haplotypes (the latter shouldn't ever happen without missing data.)



        #if two missing, harder:
        miss2 <- rowSums(ifelse(haps == 0, T, F)) == 2
        if(any(miss2)){
          #get the missing haplotypes in each row:
          tm_mat <- haps[which(miss2),] #grab the violating rows
          tm_mat[tm_mat != 0] <- NA
          tm_mat[tm_mat == 0] <- 1
          tm_mat <- t(t(tm_mat)*(1:ncol(haps)))
          tm_mat[!is.na(tm_mat)] <- colnames(haps)[as.numeric(tm_mat[!is.na(tm_mat)])]
          tm_mat <- matrix(t(tm_mat)[!is.na(t(tm_mat))], ncol = 2, byrow = T)

          #if both are actually polymorphic, assign a one, otherwise asign a 0.
          check <- ifelse(substr(tm_mat[,1], 1, nchar(tm_mat[1,1])/2) !=
                            substr(tm_mat[,2], 1, nchar(tm_mat[1,1])/2) &
                            substr(tm_mat[,1], (nchar(tm_mat[1,1])/2) + 1, nchar(tm_mat[1,1])) !=
                            substr(tm_mat[,2], (nchar(tm_mat[1,1])/2) + 1, nchar(tm_mat[1,1])),
                          1,NA)
          Dprime[which(miss2)] <- check
          rsq[which(miss2)] <- check
          rm(tm_mat)
        }

        #get pvals
        chisqu <- rsq*(4)
        pval <- 1 - pchisq(chisqu, 1)

        #remove stuff to clear memory
        rm(A1B1f, A1B2f, A2B1f, A2B2f, haps, miss2, D, D2, chisqu, B1f, B2f, A1f, A2f)

        #write
        if(prox_table){
          data.table::set(prox, (prog + 1):prog_after, j = "p1", value = meta$position[i])
          data.table::set(prox, (prog + 1):prog_after, j = "p2", value = meta$position[(i+1):nrow(x)])
          data.table::set(prox, (prog + 1):prog_after, j = "rsq", value = rsq)
          data.table::set(prox, (prog + 1):prog_after, j = "Dprime", value = Dprime)
          data.table::set(prox, (prog + 1):prog_after, j = "pval", value = pval)
        }
        if(matrix_out){
          #reminder: columns start at locus two, rows start at locus one (but end at nlocus - 1)
          fill <- rep(NA, nrow(x) - length(Dprime) - 1)
          Dprime <- c(fill, Dprime)
          rsq <- c(fill, rsq)
          pval <- c(fill, pval)
          data.table::set(rmat, i, 1:ncol(rmat), as.list(rsq))
          data.table::set(Dpmat, i, 1:ncol(Dpmat), as.list(Dprime))
          data.table::set(pvmat, i, 1:ncol(pvmat), as.list(pval))
        }
        prog <- prog_after
      }

      ###################################
      #finish and return output
      if(prox_table){
        prox$proximity <- abs(prox$p1 - prox$p2)
        if(any(colnames(meta) == "group")){
          prox$group <- meta[1,"group"]
        }
        if(any(colnames(meta) == "pop")){
          prox$pop <- meta[1,"pop"]
        }
        prox <- prox[,c(6:ncol(prox), 1:5)]
      }
      if(matrix_out){
        Dpmat <- cbind(position = meta$position[1:stop.row], Dpmat)
        rmat <- cbind(position = meta$position[1:stop.row], rmat)
        pvmat <- cbind(position = meta$position[1:stop.row], pvmat)
      }
      if(prox_table){
        if(matrix_out){
          return(list(prox = prox, Dprime = Dpmat, rsq = rmat, pval = pvmat))
        }
        else{
          return(prox = prox)
        }
      }
      else{
        return(list(Dprime = Dpmat, rsq = rmat, pval = pvmat))
      }
    }
  }
}


#gets the position and group of a set of data where the first column is: group_start_end and the second
#is position in tag. Accepts one or two substring group names (groupXXX_start_end_F/R or group_XXX_start_end_F/R).
get_position_group <- function(data){
  data[,1] <- as.character(data[,1])
  cat("Subseting out identifier rows...")
  w_df <- data[,1:2]
  cat("done.\n")
  #print(w_df)
  i <- 1
  for(i in i:nrow(w_df)){
    if(i %% 10000 == 0){cat("locus number: ", i, "\n")}
    w_v <- unlist(strsplit(w_df[i,1], "_"))
    if(length(w_v) == 4){
      w_df[i,"group"] <- w_v[1]
      start <- w_v[2]
    }
    else if(length(w_v) == 5){
      w_df[i,"group"] <- paste0(w_v[1], "_", w_v[2])
      start <- w_v[3]
    }
    else{
      warning("Incorrect location identifier format.")
      stop
    }
    w_df[i,"position"] <- as.numeric(start) + (w_df[i,2] - 1)
  }
  out_df <- w_df[,3:4]
  cat("done.\n")
  return(out_df)
}

#runs get_position_group over multiple processors. For large datasets, this speeds up things considerably.
get_position_group_par <- function(data, num_cores){
  data[,1] <- as.character(data[,1])
  cat("Subseting out identifier rows...")
  w_df <- data[,1:2]
  cat("done.\nRunning on ", num_cores, "cores.")
  #print(w_df)
  cl <- makeCluster(num_cores)
  registerdoParallel::(cl)
  output <- foreach(i=1:nrow(w_df), .inorder = TRUE, .combine = rbind) %dopar%{
    #if(i %% 10 == 0){cat("locus number: ", i, "\n")}
    w_v <- unlist(strsplit(w_df[i,1], "_"))
    if(length(w_v) == 4){
      group <- w_v[1]
      start <- w_v[2]
    }
    else if(length(w_v) == 5){
      group <- paste0(w_v[1], "_", w_v[2])
      start <- w_v[3]
    }
    else{
      warning("Incorrect location identifier format.")
      stop
    }
    c(group, as.numeric(start) + (w_df[i,2] - 1))
  }
  #print(output)
  cat("Pasting together output...")
  name_df <- as.data.frame(output)
  colnames(name_df) <- colnames(w_df)
  out_df <- cbind(name_df,data[,3:ncol(data)])
  rownames(out_df) <- c(1:nrow(out_df))
  cat("done.\n")
  stopCluster(cl)
  registerDoSEQ()
  return(out_df)
}


#' Interpolates stats at genes.
#'
#' \code{gene_ave_stat} takes an input file of genes with start/end position noted and interpolates back the average statitstic of choice for that gene. Requires the "zoo" package.
#'
#'Description of gene_data:
#'    Requires columns titled "group", "start", "end", and "probeID", containing gene linkage group/chromosome, start and end positions, and name, respectively.
#'
#'Description of stat_data:
#'    Requires columns titled "group", "position", and one matching the stat argument, containing the linkage group/chromosome, position, and the value of the observed statistic at that position. Typically produced via windowed gaussian smoothing (from smoothed_ave).
#'
#'Uses spline interpolation from the "zoo" package.
#'
#' @param x Input gene data, with columns named "start", "end" and "probeID".
#' @param y Input stat data, typically from smoothed windows or (not recommended) raw statistics from SNPs.
#' @param stat Name of the statistic to interpolate.
#' @examples
#' gene_ave_stat(stickleGO, randSMOOTHed[randSMOOTHed$pop == "A",], "smoothed_pi")
#'
gene_ave_stat <- function(x, y, stat){

  gdat <- x
  sdat <- y

  if(nrow(y) == 0){
    warning("No stat data provided.\n")
    return()
  }

  #prepare data
  gdat$mid <- rowMeans(gdat[,c(2:3)])
  gdat$probeID <- as.character(gdat$probeID)
  cdf <- data.frame(mid = c(gdat$mid, sdat$position),
                    stat = c(rep(NA, nrow(gdat)), sdat[,stat]),
                    probeID = c(gdat$probeID, rep("snp", nrow(sdat))))

  #use zoo to interplate
  cdf <- dplyr::arrange(cdf, mid)
  cdf <- zoo::zoo(cdf)
  zoo::index(cdf) <- cdf$mid
  cdf$stat <- zoo::na.spline(cdf$stat)

  #clean up
  cdf <- as.data.frame(cdf, stringsAsFactors = F)
  cdf <- cdf[cdf$probeID != "snp",]
  cdf$stat <- as.numeric(cdf$stat)
  cdf$probeID <- as.character(cdf$probeID)
  cdf$mid <- as.numeric(as.character(cdf$mid))

  #remerge
  gdat <- merge(gdat, cdf, by = c("probeID", "mid"))
  gdat <- gdat[,colnames(gdat) != "mid"]
  colnames(gdat)[which(colnames(gdat) == "stat")] <- stat

  return(gdat)
}

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




#calculates a pop weighted stat. Needs columns named "n_total" and one with a name matching
#the "stat" argument, which should be a character string. "boots" argument is the number of bootstraps
#to do do get standard errors




#estimate S in one generation

#'Selection coefficients.
#'
#'Estimates the selection coefficient \emph{S} across one generation given intial and final allele frequencies, using an equation adapted from Gillespie (2010) Population genetics: a concise guide, equation 3.2.
#'
#' @param p_before Intial frequency of allele \emph{p}.
#' @param p_after Frequency of allele \emph{p} after one generation.
#'
#' @return A numeric value, the estimate of \emph{s}.
#'
#' @examples
#' estimate_selection(.01, .011)
#'
estimate_selection <- function(p_before, p_after){
  top <- 1-(p_before/p_after)
  q <- 1 - p_before
  bottom <- q^2
  s <- top/bottom
  return(s)
}


#Function to test for deviation from HWE via perumtaiton (given small cell entries) for each snp.
#inputs: x: data, snps in subsequent columns
#        num_start: column index with the first column of data
#        miss: CHARACTER entry which encodes missing data
#        test: the test to be used.
#            options: exact: Exact test, from Wigginton et al 2005.
#                     permutation: permutation chi.square test.
#        n.reps: number of permutations to get p.values. Needed if test = "permutation".

#'SNP data analysis with snpR
#'
#'\code{snpR.stats} calls many snpR functions based on the provided arguments to do the majority of the basic statistical process required for SNP analysis.
#'
#'
#'
#'Description of x:
#'    Contains metadata in columns 1:ecs. Remainder of columns contain genotype calls for each individual. Each row is a different SNP, as given by format_snps output options 4 or 6. Requires the column containing the position of the loci in base pairs be named "position". If population info is given, it must be contained in a column named "pop".
#'
#' @param x Data frame. Input data, in either numeric or character formats, as given by format_snps output options 2 or 6.
#' @param ecs Numeric value. Number of extra metadata columns at the start of x.
#' @param stats Character vector, default "basic" (see description). Which stats should be calculated?
#' @param filter Logical, default FALSE. Should the SNP data be filtered?
#' @param smooth Character vector or FALSE, default FALSE. Should variables be smoothed? If so, which?
#' @param bootstrap Character vector or FALSE, default FALSE. Should bootstraps of smoothed windows be conducted? If so, which variables?
#' @param graphs Character vector or FALSE, default FALSE. Which variables should be plotted? Plots \emph{smoothed} values of most statistics save population structure, LD, and Tajima's D. As such, requires the statistics to be set via the "smooth" argument.
#' @param sigma Numeric, default NULL. Smoothing statistic/window size, in kb.
#' @param ws Numeric, default NULL.  How far should the window slide when smoothing? If null, a window is produced centered at every SNP.
#' @param pop Table of length > 1, default NULL. Table containing the number of samples in each population. Samples should be sorted by population and ordered identically in data.
#' @param nk Logical, default TRUE. Should contributions to smoothed windows be weighted by the number of called alleles at each site?
#' @param input Character string, default "NN". How are genotypes noted? Supports "NN" and "0000". See \code{\link{format_snps}}.
#' @param mDat Character string, default "N". How are missing \emph{alleles} coded in x (typically "N" or "00")?
#' @param maf FALSE or numeric between 0 and 1, default FALSE. Minimum acceptable minor allele frequency during filtering.
#' @param hf_hets FALSE or numeric between 0 and 1, default 0.55 Maximum acceptable heterozygote frequency during filtering.
#' @param bi_al Logical, default TRUE. Should non-biallelic SNPs be removed during filtering?
#' @param non_poly Logical, default TRUE. Should non-polymorphic loci be removed during filtering?
#' @param min_ind FALSE or integer, default FALSE. Minimum number of individuals in which a loci must be sequenced during filtering.
#' @param min_loci FALSE or numeric between 0 and 1, default FALSE. Minimum proportion of SNPs an individual must be genotyped at to be kept during filtering.
#' @param pop_maf Logical, default FALSE. For filtering, should loci be kept if they are above the minor allele frequency cuttoff in \emph{any} population?
#' @param filt_re FALSE, "partial", or "full", default "partial". How should loci be re_filtered after individuals are filtered? See \code{\link{filter_snps}}.
#' @param Fst_method Character vector, default "genepop". Which FST estimator should be used? See documentation for options.
#' @param LD_level Character vector or FALSE, default c("group", "pop"). Names of metadata columns by which to split the data for LD calculations.
#' @param LD_chr Character vector, default "all". Names of linkage groups/chromsomes for which pairwise LD values should be calculated at all SNPs. Identifiers must be contained in a column named either "group" or "chr".
#' @param LD_g_par Logical, default NULL. Shoud LD be done in parallel across multiple linkage groups/chromsomes?
#' @param tsd_level Character vector or FALSE, default c("group", "pop"). Names of metadata columns by which to split the data for Tajima's D calculations.
#' @param smooth_levs Character vector or FALSE, default c("group", "pop"). Names of metadata columns by which to split the data for smoothing.
#' @param boot_levs Character vector or FALSE, default c("group", "pop"). Names of metadata columns by which to split the data for bootstraps.
#' @param p.alt Character string, default "two-sided". If bootstraps are done, which alternative hypothesis should be used? See documentation for options.
#' @param num_boots Numeric, default NULL. If bootstraps are done, how many should be performed per stat?
#' @param par_cores Numeric, default FALSE. If parallel processing is requested, how many cores should be used?
#'
#' @return A named list of objects corresponding to requested statistics/graphs/ect.
#'
#' @examples
#'
snpR.stats <- function(x, ecs, stats = "basic", filter = FALSE, smooth = FALSE, bootstrap = FALSE, graphs = FALSE, sigma = NULL, ws = NULL,
                       pop = NULL, nk = TRUE, input = "NN", mDat = "NN", maf = 0.05, hf_hets = 0.55,
                       bi_al = TRUE, non_poly = TRUE, min_ind = FALSE, min_loci = FALSE, pop_maf = FALSE,
                       filt_re = "partial", Fst_method = "genepop", LD_level = c("group", "pop"), LD_chr = "all",
                       LD_g_par = FALSE, tsd_level = c("group", "pop"),
                       smooth_levs = c("group", "pop"), internal_boot_level = "group", boot_split_levels = "pop",
                       p.alt = "two-sided", num_boots = NULL, par_cores = FALSE){

  x <- dplyr::arrange(x, group, position)
  stab <- FALSE #set that there is no saved out.tab currently

  ####################################################
  #figure out which stats we are getting

  #vector of possible stats to request:
  pstats <- c("pi", "ho", "fst", "pa", "ld", "tsd")

  #if the default (basic), set the stats.
  if(stats == "basic"){
    stats <- c("pi", "ho", "fst", "pa")
  }

  #otherwise, check that the options given are acceptable.
  else{
    stats <- tolower(stats) #incase they used caps
    fstats <- which(!(stats %in% pstats)) #any undefined stats?
    if(length(fstats) > 0){
      stop(cat("Incorrect stats defined:", stats[fstats], ". Acceptable stats:", pstats, "\n"))
    }
    #otherwise good to go.
  }

  #how about smoothing?
  if(is.character(smooth)){
    pstats <- pstats[-length(pstats)]
    smooth <- tolower(smooth)
    fsmooth <- which(!(smooth %in% stats) | smooth == "tsd") #any unaccepted stats?
    if(length(fsmooth) > 0){
      stop(cat("Smoothing requested on incorrect variables:", smooth[fsmooth], ". Acceptable stats:", stats[-which(stats == "tsd")], "\n"))
    }
  }

  #pop info, needs more work
  if(is.table(pop)){
    if(length(pop) > 1){
      pop <- list(names(pop), as.numeric(pop))
    }
    else{
      stop("Provided pop table must have more than one population!")
    }
  }
  else if(!is.null(pop)){
    stop("Pop info should be provided as a table, with entries in the same order as samples in x (typically alphabetical)!")
  }


  if(Fst_method == "all"){
    Fst_method <- c("genepop", "Wier", "WC", "Hoh")
  }

  ###################################################
  #initialize output list
  out <- list()

  ###################################################
  #format to NN, fix ecs with filler data if it is too low
  if(ecs < 2){
    nfill <- 2 - ecs
    fm <- matrix(paste0("filler", 1:nrow(x)), nrow(x), nfill)
    fm <- as.data.frame(fm, stringsAsFactors = F)
    colnames(nfill) <- paste0("filler", 1:nfill)
    x <- cbind(fm, x)
    remove(fm)
    ecs <- 2
  }
  else{
    nfill <- FALSE
  }

  if(input != "NN"){
    x <- format_snps(x, ecs = ecs, output = 6, input_form = input, mDat = mDat)
    if(nfill > 0){
      out$ch.form.x <- x[,-c(1:nfill)]
    }
    else{
      out$ch.form.x <- x
    }
  }

  ###################################################
  #filter the data as requested
  if(filter){
    if(pop_maf == TRUE){
      x <- filter_snps(x, ecs, maf, hf_hets, min_ind, min_loci, filt_re, pop, non_poly, bi_al, mDat, out.tab = TRUE)
    }
    else{
      x <- filter_snps(x, ecs, maf, hf_hets, min_ind, min_loci, filt_re, FALSE, non_poly, bi_al, mDat, out.tab = TRUE)
    }
    #pull out the out.tab
    stab <- x[-(names(x) == "x")]
    x <- x$x

    #add to return
    if(nfill > 0){
      out$x.flt <- x[,-c(1:nfill)]
    }
    else{
      out$x.flt <- x
    }
  }


  ###################################################
  #do each stat as requested and save output

  #get ac format, used in most analysis and for nk
  sink(tempfile())
  x_ac <- format_snps(x, ecs, 1, pop = pop, in.tab = stab)
  sink()
  if(is.list(pop) & length(pop) > 1){
    x_ac <- dplyr::arrange(x_ac, pop, group, position)
  }
  else{
    x_ac <- dplyr::arrange(x_ac, group, position)
  }

  #initialize output stat dataframe, also used for smoothing.
  sm.smooth <- smooth[smooth != "fst"]
  sm.in <- matrix(NA, nrow = nrow(x_ac), ncol = ecs + length(sm.smooth) + 1)
  sm.in <- as.data.frame(sm.in)
  sm.in[,1:ecs] <- x_ac[,1:ecs]
  colnames(sm.in) <- c(colnames(x_ac)[1:ecs], sort(sm.smooth), "nk")
  sm.in$nk <- x_ac$n_total
  #add a pop column if info provided...
  if(length(pop) > 1){
    sm.in$pop <- x_ac$pop
  }

  #pi
  if("pi" %in% stats | ("Hoh" %in% Fst_method & "fst" %in% stats)){
    pi <- calc_pi(x_ac)
    sm.in$pi <- pi
  }

  #ho
  if("ho" %in% stats | (("WC" %in% Fst_method | "Wier" %in% Fst_method) & "fst" %in% stats)){
    if(!is.null(pop)){
      ho <- calc_Ho(x, ecs)
    }
    else{
      ho <- calc_Ho(x, ecs, pop = pop)
    }
    sm.in$ho <- ho$ho
  }

  #Fst
  if("fst" %in% stats){
    fst.list <- character(0)
    if("genepop" %in% Fst_method){

      #make the genepop input
      if(file.exists("gp_fst.txt")){
        file.remove("gp_fst.txt")
      }
      gp <- format_snps(x, ecs, 2, pop = pop, outfile = "gp_fst.txt")
      remove(gp)

      #run Fst
      fst <- calc_pairwise_Fst("gp_fst.txt", method = "Genepop", pnames = pop[[1]])
      if(nfill != FALSE){
        fst$loci <- cbind(x[,(nfill+1):ecs], fst$loci)
        out$fst$genepop$snps <- fst$loci
      }
      else{
        fst$loci <- cbind(x[,1:ecs], fst$loci)
        out$fst$genepop$loci <- fst$loci
      }
      out$fst$genepop$overall <- fst$overall
      fst <- fst$loci

      fst.list <- "genepop"

      fst <- reshape2::melt(fst, id.vars = colnames(x)[1:ecs])
      colnames(fst)[colnames(fst) == "value"] <- "fst.genepop"
    }
    if("WC" %in% Fst_method | "Wier" %in% Fst_method){
      #get data ready
      ho <- reshape2::melt(ho, id.vars = colnames(x_ac)[1:ecs])
      ho$variable <- as.character(ho$variable)
      x_ac$Ho <- ho$value

      if("Wier" %in% Fst_method){
        fst.Wier <- calc_pairwise_Fst(x_ac, ecs, method = "Wier")
        out$fst$Wier <- fst.Wier
        fst.list <- c(fst.list, "Wier")
        fst.Wier <- reshape2::melt(fst.Wier, id.vars = colnames(x)[1:ecs])
        if(exists("fst")){
          fst <- cbind(fst, fst.Wier = fst.Wier$value)
        }
        else{
          fst <- fst.Wier
        }
      }
      if("WC" %in% Fst_method){
        fst.WC <- calc_pairwise_Fst(x_ac, ecs, method = "WC")
        out$fst$WC <- fst.WC
        fst.list <- c(fst.list, "WC")
        fst.WC <- reshape2::melt(fst.WC, id.vars = colnames(x)[1:ecs])
        if(exists("fst")){
          fst <- cbind(fst, fst.WC = fst.WC$value)
        }
        else{
          fst <- fst.WC
        }
      }
    }
    if("Hoh" %in% Fst_method){ #hohenlohe
      x_ac$pi <- pi
      fst.Hoh <- calc_pairwise_Fst(x_ac, ecs, method = "Hohenlohe")
      out$fst$Hoh <- fst.Hoh
      fst.list <-c(fst.list, "Hoh")
      fst.Hoh <- reshape2::melt(fst.Hoh, id.vars = colnames(x)[1:ecs])
      if(exists("fst")){
        fst <- cbind(fst, fst.Hoh = fst.Hoh$value)
      }
      else{
        fst <- fst.Hoh
      }
    }


    if("fst" %in% smooth){
      if(nk){
        pnk <- calc_pairwise_nk(x_ac, ecs)
      }
      pnk <- reshape2::melt(pnk, id.vars = colnames(x)[1:ecs])
      fst$nk <- pnk$value
      colnames(fst)[which(colnames(fst) == "variable")] <- "pop"
      out$fst$smooth <- s_ave_multi(fst,
                                    colnames(fst)[which(grepl("fst", colnames(fst)))],
                                    sigma, ws, nk, smooth_levs)
    }
  }

  #private alleles
  if("pa" %in% stats){
    pa <- check_private(x_ac)
    pa <- cbind(x[,1:ecs], pa)
    pa <- reshape2::melt(pa, id.vars = colnames(x_ac[,1:ecs]))
    colnames(pa)[colnames(pa) == "variable"] <- "pop"
    colnames(pa)[colnames(pa) == "value"] <- "pa"
    sm.in$pa <- pa$pa
  }

  #tsd
  if("tsd" %in% stats){
    warning("Tajima's D should be run on a full dataset, without filtering SNPs and containing non-segregating genotypes.")

    #are we splitting tsd by pop?
    if("p" %in% tsd_level){

      #are we also splitting by group?
      if("g" %in% tsd_level){
        tsd <- run_gp(x_ac, Tajimas_D, ws = ws, step = step)
      }

      else{
        tsd <- Tajimas_D(x_ac[x_ac$pop == pop[[1]][1],], ws, step)
        for(i in 2:length(pop[[1]])){
          tsd <- rbind(tsd, Tajimas_D(x_ac[x_ac$pop == pop[[1]][i],], ws, step))
        }
      }

    }

    #splitting by just group
    else if("g" %in% tsd_level){
      tsd <- run_g(x_ac, Tajimas_D, ws = ws, step = step)
    }


    #not splitting?
    else{
      tsd <- Tajimas_D(x_ac, ws, step)
    }

    out$tsd <- tsd
  }


  #LD, this one might be slow...
  if("ld" %in% stats){

    #are we running only certain groups/chr?
    if(LD_chr != "all"){
      xLD <- x[which(x$group %in% LD_chr),]
    }
    else{
      xLD <- x
    }

    #what level are we running this at?
    #are we splitting by pops?
    if("pop" %in% LD_level){

      #initialize pop output list
      pLD <- vector("list", length = length(pop))
      names(pLD) <- pop[[1]]

      #initialize column tracker
      tracker <- ecs + 1

      for(i in 1:length(pop)){
        pLD[[i]] <- LD_full_pairwise(x = cbind(LDx[,1:ecs], LDx[, tracker:(tracker + pop[[2]][i] - 1)]),
                                     ecs, levels = LD_level[-which(LD_level == "pop")], par = par_cores)
        tracker <- tracker + pop[[2]][i]
      }
    }


    #we aren't splitting by pops?
    else{
      pLD <- LD_full_pairwise(LDx, ecs, par = par_cores, levels = LD_level)
    }
    out$LD <- pLD
  }

  out$stats <- sm.in

  ###########################################
  #smooth any other requested stats

  browser()
  if(is.character(smooth)){
    sm.out <- s_ave_multi(sm.in, colnames(sm.in)[(ecs+1):(ecs + length(sm.smooth))], sigma, ws, nk, smooth_levs)
    out$stats_smoothed <- sm.out
  }

  ###########################################
  #bootstrap any requested stats
  browser()
  #if we are bootstrapping...
  if(bootstrap != FALSE){

    #If all fst. methods were requested, need to seperate these out.
    if(any(grepl("fst.all", bootstrap))){
      bootstrap <- c(bootstrap[-which(grepl("fst", bootstrap))], "fst.Hoh", "fst.genepop", "fst.WC", "Fst.Wier")
    }

    #make a matrix detailing which runs to do exactly:
    if(boot_split_levels == FALSE){
      boot_matrix <- matrix(bootstrap, ncol = 1)
    }
    else{
      boot_matrix <- matrix(NA, ncol = 1 + length(boot_split_levels))
      for(i in 1:length(bootstrap)){
        if(grepl("fst", bootstrap[i])){
          ul <- matrix(as.character(unique(fst[,colnames(fst) %in% boot_split_levels])), ncol = length(boot_split_levels))
          boot_matrix <- rbind(boot_matrix, cbind(bootstrap[i], ul))
        }
        else{
          ul <- matrix(as.character(unique(out$stats[,colnames(out$stats) %in% boot_split_levels])), ncol = length(boot_split_levels))
          boot_matrix <- rbind(boot_matrix, cbind(bootstrap[i], ul))
        }
      }
      boot_matrix <- boot_matrix[-1,]
    }

    #prepare things for parallel run if requested:
    if(par_cores != FALSE){
      cl <- parallel::makePSOCKcluster(par_cores)
      doParallel::registerDoParallel(cl)

      ntasks <- nrow(boot_matrix)
      # progress <- function(n) cat(sprintf("Part %d out of",n), ntasks, "is complete.\n")
      # opts <- list(progress=progress)
    }

    #function to grab correct input data:
    gbdat <- function(tstat){
      if(grepl("fst", tstat)){

        #grab the correct data
        if(tstat == "fst.genepop"){
          tx <- out$fst$genepop$loci
        }
        else{
          tx <- out$fst[[which(grepl(substr(tstat, 5, nchar(tstat)), names(out$fst)))]]
        }

        #melt it to split pops
        tx <- reshape2::melt(tx, id.vars = colnames(tx)[1:ecs])
        colnames(tx)[-c(1:ecs)] <- c("comp", tstat)

        #get fws
        if(!is.null(ws)){
          tfws <- split(data.table::data.table(out$fst$smooth), by = boot_levs, sorted = TRUE)
        }
        else{
          tfws <- NULL
        }

        return(list(dat = tx, fws = tfws))
      }
      else{
        tdat <- out$stats
        tfws <- ifelse(is.null(ws), NULL, out$stats_smoothed)
        return(list(dat = tdat, fws = fws))
      }
    }

    #call if parallel:
    if(par_cores != FALSE){
      foreach::foreach(q = 1:ntasks, .inorder = TRUE
                       ) %dopar% {

                         tdat <- gbdat(tstat = boot_matrix[i,1])
                         tfws <- tdat$fws
                         tdat <- tdat$dat

                         if(ncol(boot_matrix) > 1){
                           #get only the correct data
                           tdat <- tdat[which(apply(tdat, 1, function(x) identical(x[which(colnames(x) %in% boot_split_levels),], boot_matrix[i,-1]))),]
                           tfws <- tfws[which(apply(tfws, 1, function(x) identical(x[which(colnames(x) %in% boot_split_levels),], boot_matrix[i,-1]))),]
                         }

                         b_output <- resample_long(tdat, boot_matrix[i,1], num_boots, sigma, nk, fws, TRUE, level = internal_boot_level)
                       }
    }
    else{
      b_output <- vector("list", length = nrow(boot_matrix))
      for(i in 1:nrow(boot_matrix)){

        tdat <- gbdat(tstat = boot_matrix[i,1])
        tfws <- tdat$fws
        tdat <- tdat$dat

        if(ncol(boot_matrix) > 1){
          #get only the correct data
          tdat <- tdat[which(apply(tdat, 1, function(x) identical(x[which(colnames(x) %in% boot_split_levels),], boot_matrix[i,-1]))),]
          tfws <- tfws[which(apply(tfws, 1, function(x) identical(x[which(colnames(x) %in% boot_split_levels),], boot_matrix[i,-1]))),]
        }

        b_output[[i]] <- resample_long(tdat, boot_matrix[i,1], num_boots, sigma, nk, fws, TRUE, level = internal_boot_level)
      }
    }

    #process output:

  }


  ###########################################
  #return the object
  return(out)

}




#function to run any command after spliting the data by group and by population. Group and population must be in columns named "group" and "pop", respectively.
#The first argument of the function to run must be the data provided to that function.
run_gp <- function(x, FUN, ...){
  w_df<- data.frame()
  x$pop <- as.character(x$pop)
  x$group <- as.character(x$group)
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
        if(!any(colnames(x_data) == "pop"))
          x_data <- cbind(pop = pops[j], x_data, stringsAsFactors = F)
        if(!any(colnames(x_data) == "group")){
          x_data <- cbind(group = groups[i], x_data, stringsAsFactors = F)
        }
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
  if(sum(colnames(w_df) == "pop") > 1){
    excols <- which(colnames(w_df) == "pop")[-1]
    w_df <- w_df[,-excols]
  }
  return(w_df)
}

#function to run any command after spliting the data by group. Group must be in a column named "group".
#The first argument of the function to run must be the data provided to that function.
run_g <- function(x, FUN, ..., y = NULL){
  w_df<- data.frame()
  x$group <- as.character(x$group)
  groups <- unique(x[,"group"])
  for (i in 1:length(groups)){
    w_data <- x[x[,"group"] == groups[i],]
    if(!is.null(y)){
      y_dat <- y[y[,"group"] == groups[i],]
      w_data <- FUN(w_data, y = y_dat, ...)
    }
    else{
      w_data <- FUN(w_data, ...)
    }
    print(groups[i])
    if(length(w_data) == 0){next}
    if(is.data.frame(w_data)){
      if(!any(colnames(w_data) == "group")){
        w_data <- cbind(group = groups[i], w_data, stringsAsFactors = F)
      }
      w_df <- rbind(w_df, w_data)
    }
    else{
      if(is.data.frame(w_df)){
        w_df <- numeric(0)
      }
      w_df <- c(w_df, w_data)
    }
  }
  if(!is.data.frame(w_df)){
    warning("Output is vector: if merging with meta info, check that info data frame is sorted identically:
            group then snp by snp!\n")
  }
  return(w_df)
}




#' Refine a random forest model via sequential removal of uninformative SNPs.
#'
#' Improves the prediction accuracy of a random forest model via iterative
#' removal of uninformative SNPs. In each step, the SNPs with the lowest
#' absolute value importance are removed from the model. Depending on the
#' provided arguments, multiple trim percentages can be provided for different
#' SNP number cuttoffs.
#'
#' Random Forest models can fail to predict well with "noisy" data, where most
#' explanitory variables are uniformative. Since most whole-genome sequence data
#' is like this, it can be useful to "trim" a data to remove uniformative snps.
#' Since even models constructed on "noisy" data tend to pick out the most
#' important SNPs with a decent degree of accuracy, an initial random forest
#' model can be used to select SNPs to remove. Sequential removal of unimportant
#' SNPs in this way can radically improve prediction accuracy, although SNP
#' p-values stop being informative.
#'
#' Multiple trim levels can be specified, which will determine the percentage of
#' SNPs removed according to different cuttoff levels of remaining SNPs. For
#' example, providing trim levels of 0.9, 0.5, and 0.1 with cuttoffs of 1000 and
#' 100 will trim 90% of SNPs untill 1000 remain, then trim 50% untill 100
#' remain, then trim 10% thereafter.
#'
#' If less trim levels are provided than needed for each cuttoff, a single SNP
#' will be removed each step below the final cuttoff. For example, providing
#' trim levels of 0.9 and 0.5 and cuttoffs of 1000 and 100 will trim 90% of SNPs
#' untill 1000 remain, then 50% untill 100 remain, then one at a time.
#'
#' If the data contains informative SNPs, prediction accuracy should improve on
#' average untill informative SNPs begin to be removed, at which point accuracy
#' will decrease.
#'
#' mtry will be set to the number of SNPs in each run.
#'
#' As usual, facets can be requested. However, in this case, only a single facet
#' and facet level (subfacet) may be provided at once, and must match a facet
#' and level in the provided input model.
#'
#' @param rf The random forest model to be refined. List containg snpRdata
#'   object (named $data) and a \code{\link[ranger]{ranger}} model, named
#'   $models$x$model, where x is the facet/subfacet. Identical to ojects created
#'   by \code{\link{run_random_forest}}.
#'@param response character. Name of the column containing the response
#'   variable of interest. Must match a column name in sample metadata. Response
#'   must be categorical, with only two categories.
#' @param facets Character, default NULL. Facet to run. Only a single facet and
#'   facet level (subfacet) may be provided at once, and must match a facet and
#'   level in the provided input model. If NULL, runs the base level facet.
#' @param subfacet Character, default NULL. Facet level (subfacet) to run. Only
#'   a single facet and facet level (subfacet) may be provided at once, and must
#'   match a facet and level in the provided input model. If NULL, runs the base
#'   level facet.
#' @param formula charcter, default NULL. Model for the response variable, as
#'   described in \code{\link[stats]{formula}}. If NULL, the model will be
#'   equivalent to response ~ 1.
#' @param num.trees numeric, default 10000. Number of trees to grow. Higher
#'   numbers will increase model accuracy, but increase calculation time. See
#'   \code{\link[ranger]{ranger}} for details.
#' @param trim numeric, default 0.5. Percentages of SNPs to be trimmed between
#'   model iterations. Multiple trim levels can be provided corresponding to
#'   different trim cuttoffs. If less trim levels are provided than needed to
#'   describe every trim_cuttoff interval, will trim a single SNP below the
#'   final cuttoff. See details.
#' @param trim_cuttoffs numeric, default NULL. Specifies the number of SNPs
#'   below which to change trim percentages. If NULL, the default, trims at the
#'   given level untill 1 SNP remails. See details.
#' @param importance character, default "impurity_corrected". The method by
#'   which SNP importance is determined. Options: \itemize{\item{impurity}
#'   \item{impurity_corrected} \item{permutation}}. See
#'   \code{\link[ranger]{ranger}} for details.
#' @param interpolate character, default "bernoulli". Interpolation method for
#'   missing data. Options: \itemize{\item{bernoulli: }binomial draws for the
#'   minor allele. \item{af: } insertion of the average allele frequency}.
#' @param par numeric, default FALSE. Number of parallel computing cores to use
#'   for computing RFs across multiple facet levels or within a single facet if
#'   only a single category is run (either a one-category facet or no facet).
#' @param ... Additional arguments passed to \code{\link[ranger]{ranger}}.
#'
#' @return A list containing: \itemize{\item{error_delta: } A data.frame noting
#'   the number of SNPs and corresponding prediction_error in each model
#'   iteration. \item{confusion_matrices: } An array containing confusion
#'   matrices for categorical responses. The third subscript denotes model
#'   iteration ('[,,1]' would reference model 1.) \item{best_model: } The model
#'   with the lowest prediction error from the provided dataset, in the format
#'   provided by \code{\link{run_random_forest}}}
#'
#' @references Wright, Marvin N and Ziegler, Andreas. (2017). ranger: A Fast
#'   Implementation of Random Forests for High Dimensional Data in C++ and R.
#'   \emph{Journal of Statistical Software}.
#' @references Goldstein et al. (2011). Random forests for genetic association
#'   studies. \emph{Statistical Applications in Genetics and Molecular Biology}.
#'
#' @seealso \code{\link[ranger]{ranger}} \code{\link[ranger]{predict.ranger}}.
#' @seealso \code{\link{run_random_forest}}
#'
refine_rf <- function(rf, response, facets = NULL, subfacet = NULL,
                      formula = NULL, num.trees = 10000,
                      trim = 0.5, trim_cuttoffs = NULL,
                      importance = "impurity_corrected",
                      interpolate = "bernoulli",
                      par = FALSE,
                      ...){
  #==============sanity checks================
  msg <- character()
  if(length(facets) > 1){
    msg <- c(msg, "No more than one facet may be refined at a time.\n")
    if(is.null(subfacet)[1]){
      msg <- c(msg, "No more that one subfacet may be refined at a time. Please list a subfacet.\n")
    }
  }
  if(!is.null(subfacet[1])){
    if(length(subfacet) > 1){
      msg <- c(msg, "No more that one subfacet may be refined at a time.\n")
    }
  }
  if(!is.null(facets)){
    if(length(rf$models) != 1 | names(rf$models[[1]] != ".base_.base")){
      msg <- c(msg, "If the base facet is requested, please provide an input rf result with only the base facet calculated.\n")
    }
  }
  
  # trim checks:
  ## if trim is null, start at a single snp per iteration.
  if(is.null(trim[1])){
    trim <- 0.5
    trim_cuttoffs <- nrow(rf$data)
  }
  ## if trim_cuttoffs is null, just go until one snp left.
  if(is.null(trim_cuttoffs)){
    if(length(trim) > 1){
      msg <- c(msg, "Only one trim level allowed when no trim_cuttoffs provided.\n")
    }
    trim_cuttoffs <- 1
  }
  
  if(length(trim) + 1 > length(trim_cuttoffs)){
    msg <- c(msg, "Too many trim levels for the given cuttoffs. Only one more trim level may be given than cuttoffs.\n")
  }
  
  
  
  #==============the actual refinement function==================
  run_refine <- function(rf, response, formula, num.trees, trim_cuttoffs, trim, search_cuttoff, par,
                         ...){
    
    #=============subfunctions:=========
    # trim snps below a given importance quantile
    remove_trim <- function(imps, trim){
      best.imp <- quantile(imps, trim)
      best.imp <- which(imps >= best.imp[1])
      return(best.imp)
    }
    # trim to a specific number of snps
    remove_set_to_n_snps <- function(imps, n_snps){
      imps <- data.frame(imps, ord = 1:length(imps))
      imps <- dplyr::arrange(imps, desc(imps))
      imps <- imps[1:n_snps,]
      imps <- dplyr::arrange(imps, ord)
      best.imp <- imps$ord
      return(best.imp)
    }
    # remove a single snp
    remove_single <- function(imps){
      best.imp <- which.min(imps)
      best.imp <- (1:nrow(rf$data@stats))[-best.imp]
      return(best.imp)
    }
    # set to the next level, then trim a given percentage
    remove_set_and_trim <- function(imps, n_snps, trim){
      # set
      imps <- data.frame(imps, ord = 1:length(imps))
      imps <- dplyr::arrange(imps, desc(imps))
      imps <- imps[1:n_snps,]
      imps <- dplyr::arrange(imps, ord)
      
      # trim
      best.imp <- quantile(imps[,1], trim)
      best.imp <- imps[which(imps[,1] >= best.imp[1]),]
      
      # return
      best.imp$ord
    }
    
    # remove according to cuttoffs and trim levels
    remove_snps <- function(imps, trim_cuttoffs, trim){
      
      # figure out which "bin" we are in
      n_snps <- length(imps)
      bin.logi <- which(n_snps > trim_cuttoffs)
      if(length(bin.logi) > 0){
        
        # grab the current bin and trim
        bin <- min(bin.logi)
        best.imp <- remove_trim(imps, trim[bin])
        
        # check if we changed bins!
        bin.t.trim <- which(length(best.imp) > trim_cuttoffs)
        if(!identical(bin.t.trim, bin.logi)){
          
          # if we haven't hit the last bin
          if(length(bin.t.trim) != 0){
            # we either want to set this to - the current trim or to the cuttoff - the next trim, whichever has more snps.
            
            best.imp.current.trim <- remove_trim(imps, trim[bin])
            best.imp.set.plus.next.trim <- remove_set_and_trim(imps, trim_cuttoffs[min(bin.t.trim) - 1], trim[min(bin.t.trim)])
            if(length(best.imp.current.trim) > length(best.imp.set.plus.next.trim)){
              best.imp <- best.imp.current.trim
            }
            else{
              best.imp <- best.imp.set.plus.next.trim
            }
          }
          
          
          ## if we've hit the last bin and are using a single cutoff, set to single cuttoff
          else if(length(bin.t.trim) == 0 & trim_single){
            best.imp <- remove_set_to_n_snps(imps, trim_cuttoffs[length(trim_cuttoffs)])
            if(length(best.imp) == n_snps){
              best.imp <- remove_single(imps)
            }
          }
          
          ## if we've hit the last bin but aren't using a single cuttoff, keep trimming as usual
          else if(length(bin.t.trim) == 0){
            best.imp.current.trim <- remove_trim(imps, trim[bin])
            best.imp.set.plus.next.trim <- remove_set_and_trim(imps, trim_cuttoffs[bin],  trim[length(trim)])
            if(length(best.imp.current.trim) > length(best.imp.set.plus.next.trim)){
              best.imp <- best.imp.current.trim
            }
            else{
              best.imp <- best.imp.set.plus.next.trim
            }
            
          }
        }
      }
      # if we've hit the last bin...
      else{
        if(trim_single){
          best.imp <- remove_single(imps)
        }
        else{
          best.imp <- remove_trim(imps, trim[length(trim)])
        }
      }
      
      
      
      return(best.imp)
    }
    
    #==========initialize============
    
    
    
    # intialize output
    out <- matrix(NA, 1000, 2)
    colnames(out) <- c("n_snps", "prediction_error")
    out[1,] <- c(nrow(rf$data), rf$models[[1]]$model$prediction.error)
    
    # initialize output confusion array
    if(any(names(rf$models[[1]]$model) == "confusion.matrix")){
      conf.out <- array(NA, dim = c(dim(rf$models[[1]]$model$confusion.matrix), 1000))
      conf.out[,,1] <- rf$models[[1]]$model$confusion.matrix
    }
    
    
    # initialize best model output
    best.mod <- rf
    best.error <- out[1,2]
    
    # intialize difference
    diff <- 1
    
    # find the target column name
    tar.col <- paste0(response, "_RF_importance")
    
    # if we are provided with less trim_cuttoffs than trim levels, we assume no single trimming
    if(length(trim_cuttoffs) < length(trim)){
      trim_single <- FALSE
    }
    # if the same, we trim single SNPs beneath the lowest cutoff
    else if (length(trim_cuttoffs) == length(trim)){
      trim_single <- TRUE
    }
    else{
      stop("Not enough trim levels provided for given trim cuttoffs.\n")
    }
    
    #==========while loop================
    i <- 2
    continue <- T
    while(continue == T){
      
      # if we somehow reach the end of the output storage, make it bigger!
      if(i == nrow(out)){
        out <- rbind(out, matrix(NA, 100000000, 2))
        conf.out <- c(conf.out, array(NA, dim = c(dim(rf$models[[1]]$model$confusion.matrix), 1000)))
        conf.out <- array(conf.out, c(dim(rf$models[[1]]$model$confusion.matrix), 2000))
      }
      
      imps <- abs(rf$data@stats[[tar.col]])
      
      # get the snps to keep
      best.imp <- remove_snps(imps, trim_cuttoffs, trim)
      
      # subset
      suppressWarnings(input <- subset_snpR_data(rf$data, snps = best.imp))
      if(nrow(input) == 1){
        continue <- FALSE
      }
      
      # report some things
      out[i,1] <- nrow(input)
      cat("Refinement: ", i - 1, "\n\tStarting prediction error: ", out[i-1, 2],
          "\n\tNumber of remaining SNPs: ", out[i,1], "\n\tBeginning rf...\n")
      
      # run rf
      suppressWarnings(rf <- run_random_forest(x = input, response = response,
                                               formula = formula, num.trees = num.trees,
                                               mtry = nrow(input), par = par,  importance = importance,
                                               pvals = FALSE, ...))
      
      # save outputs
      out[i,2] <- rf$models[[1]]$model$prediction.error
      if(exists("conf.out")){
        conf.out[,,i] <- rf$models[[1]]$model$confusion.matrix
      }
      if(out[i,2] < best.error){
        best.mod <- rf
        best.error <- out[i,2]
      }
      
      i <- i + 1
    }
    
    #==============return==============
    # return final model and outputs
    empties <- which(is.na(out[,1]))
    out <- out[-empties,]
    if(exists("conf.out")){
      conf.out <- conf.out[,,-empties]
      return(list(error_delta = out, confusion_matrices = conf.out, best_model = best.mod))
    }
    else{
      return(list(error_delta = out, best_model = best.mod))
      
    }
  }
  
  #==============prep============================================
  # strip down to the correct facet, doing  more sanity checks
  facets <- check.snpR.facet.request(rf$data, facets, "snp")
  o.facet <- facets
  
  # subset if need be
  if(o.facet != ".base"){
    facets <- .split.facet(facets)
    facets <- unlist(facets)
    
    # figure out which of the run rf models (facets) we are using.
    use.model <- sapply(facets, function(x) grepl(pattern = x, x = names(rf$models)))
    if(length(names(rf$models)) == 1){
      if(sum(use.model) != length(use.model)){
        msg <- c(msg, "Requested facet not found in provided rf.\n")
      }
      else{
        use.model <- 1
      }
    }
    else{
      use.model <- which(rowSums(use.model) == ncol(use.model))
      if(length(use.model) == 0){
        msg <- c(msg, "Requested facet not found in provided rf.\n")
      }
    }
    
    
    # figure out which subfacet we are using if requested.
    if(!is.null(subfacet)){
      good.mods <- names(rf$models)[use.model]
      use.model <- which(grepl(subfacet, good.mods))
      
      if(length(use.model) > 1){
        msg <- c(msg, "Subfacet + facet request matches more than one model in the provided rf.\n")
      }
      else if(length(use.model)  == 0){
        msg <- c(msg, "Subfacet + facet request matches no models in the provided rf.\n")
      }
      
      else{
        use.model <- good.mods[use.model]
      }
    }
    
    # stop if errors
    if(length(msg) > 0){
      stop(msg)
    }
    
    # subset
    o.mod <- rf$models[[use.model]]
    dat <- .subset_snpR_data(rf$data, facets = o.facet, subfacets = subfacet)
    dat <- import.snpR.data(as.data.frame(dat, stringsAsFactors = F), dat@snp.meta, dat@sample.meta, dat@mDat)
    matches <- intersect(which(rf$data@stats$facet == o.facet),
                         which(rf$data@stats$subfacet == subfacet))
    imps <- rf$data@stats[matches,]
    imps$facet <- ".base"
    imps$subfacet <- ".base"
    dat <- merge.snpR.stats(dat, imps)
    
    rf <- list(data = dat, models = list(.base_.base = o.mod))
  }
  
  
  # add sn formatted data if not present
  if(length(rf$data@sn) != 0){
    if(rf$data@sn$type != interpolate){
      sn <- format_snps(rf$data, "sn", interpolate = interpolate)
      rf$data@sn <- list(type = interpolate, sn = sn)
      
    }
  }
  else{
    sn <- format_snps(rf$data, "sn", interpolate = interpolate)
    rf$data@sn <- list(type = interpolate, sn = sn)
  }
  
  #==============run the refinement============
  
  out <- run_refine(rf, response = response,
                    formula = formula,
                    num.trees = num.trees,
                    trim_cuttoffs = trim_cuttoffs,
                    trim = trim,
                    search_cuttoff = search_cuttoff,
                    par = par,
                    ...)
  
  return(out)
}
