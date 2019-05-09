# odds ratio association testing wrapper. Needs internal and apply.snpR.dat support.
calc_odds_ratio <- function(x, facet){
  func <- function(cast_ac){
    #============odds ratio==========
    # odds ratio
    a <- cast_ac[,1]
    b <- cast_ac[,2]
    c <- cast_ac[,3]
    d <- cast_ac[,4]
    s.e <- sqrt(1/a + 1/b + 1/c + 1/d)
    odds <- log((a/b)/(c/d))

    #============chisq===============
    # pearson's chi.sq
    prob.n1 <- (a + b)/rowSums(cast_ac)
    prob.n2 <- 1 - prob.n1
    prob.case <- (a + c)/rowSums(cast_ac)
    prob.control <- 1 - prob.case

    # expected
    n <- rowSums(cast_ac)
    ea <- prob.n1*prob.case*n
    eb <- prob.n1*prob.control*n
    ec <- prob.n2*prob.case*n
    ed <- prob.n2*prob.control*n

    calc.chi <- function(e, o) ((o-e)^2)/e

    chi.a <- calc.chi(ea, a)
    chi.b <- calc.chi(eb, b)
    chi.c <- calc.chi(ec, c)
    chi.d <- calc.chi(ed, d)

    chi.stat <- chi.a + chi.b + chi.c + chi.d
    chi.p <- pchisq(chi.stat, 1, lower.tail = F)

    #============fisher's===========
    #============Armitage===========


    return(data.frame(odds.ratio = log(odds), se = s.e))
  }


  facets <- check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    invisible(capture.output(x <- add.facets.snpR.data(x, facets)))
  }

  # stuff to be added to the apply function with edits!
  tac <- cbind(x@facet.meta, x@ac)
  tac <- tac[tac$facet == "phenotype",]
  ctac <- data.table::dcast(data.table::setDT(tac), .snp.id ~ subfacet, value.var = c("ni1", "ni2"))
  ctac <- as.matrix(ctac)[,-1]

  # or for gs tables!
  gs <- data.table::as.data.table(cbind(dat2@facet.meta, dat2@geno.tables$gs))
  gs <- gs[gs$facet == "pheontype",]
  cast_gs <- data.table::dcast(data.table::setDT(gs), .snp.id ~ subfacet, value.var = colnames(dat2@geno.tables$gs))
  cast_gs <- as.matrix(cast_gs)[,-1]

}


do.drift <- function(as, gen, N){
  browser()

  # a function to do random mating in one generation
  random.mating <- function(as, N){
    browser()
    # initialize, note that since R assigns column indices by column, and we need to deal with each locus and assign within that first, we need to transpose things
    new.matrix <- matrix(NA, nrow = nrow(as), ncol = ncol(as))
    colnames(new.matrix) <- colnames(as)
    new.matrix <- t(new.matrix)
    tas <- t(as)

    # calculate p counts/frequencies
    p.count <- matrixStats::rowMaxs(as)
    p <- p.count/rowSums(as)

    # sample new p and thus q alleles for each individual
    newalleles <- rbinom(nrow(as)*N, size = 2, prob = p)
    newalleles <- matrix(newalleles, ncol = N)
    new.ps <- rowSums(newalleles)
    new.qs <- N*2 - new.ps

    # fill in the major alleles (p) for each column in the output. To do so, first make a dummy matrix where each row has the major count repeated four times,
    # then compare this to the input matrix to determine which cells hold major allele counts.
    max.matrix <- matrix(rep(p.count, each = 4), ncol = ncol(tas))
    logi.matrix <- ifelse(tas == max.matrix, TRUE, FALSE)
    vio.cols <- which(colSums(logi.matrix) == 2) # need to note any rows where frequency of p = .5, since that may mess things up.
    clean.logi <- logi.matrix[,-vio.cols] # remove vio cols
    major.slots <- which(clean.logi) # these cells are the major alleles
    new.matrix[,-vio.rows][major.slots] <- new.ps[-vio.rows] # assign to the output


    # do the same thing for the minor alleles. Note that there may be columns with no minor alleles, so we need to adjust for that as well.
    min.matrix <- matrix(rep(rowSums(as) - p.count, each = 4), ncol = ncol(tas))
    fixed.cols <- which(colSums(min.matrix) == 0) # cols with no minor allele
    min.matrix[,fixed.cols] <- -1 # set these to -1, so that can't possibly match
    logi.matrix.min <- ifelse(tas == min.matrix, T, F)
    clean.logi.min <- logi.matrix.min[,-vio.cols]
    minor.slots <- which(clean.logi.min) # cells with minor alleles
    new.matrix[,-vio.rows][minor.slots] <- new.qs[-sort(c(vio.rows, fixed.rows))] # assign, note that we need to account for both fixed and violating rows now.

    # now need to fix the violating rows.
    ## make the replacement matrix
    vio.matrix <- matrix(0, ncol = length(vio.rows), nrow = 2)
    vio.matrix[1,] <- new.ps[vio.rows]
    vio.matrix[2,] <- new.qs[vio.rows]

    ## replace
    replace.slots <- which(ifelse(tas[,vio.rows] != 0, T, F)) # the slots to replace, where the as data actually has allele counts.
    new.matrix[,vio.rows][which(ifelse(tas[,vio.rows] != 0, T, F))] <- vio.matrix

    # transpose and return
    return(t(new.matrix))
  }

  # do random mating for gen generations and save the frequencies at each step
  out.matrix <- matrix(0, nrow = nrow(as), ncol = gen)
  for(i in 1:gen){ # this loop needs to be edited to work with the new output from the random mating function.
    as.n <- random.mating(as, N)
    p <- matrixStats::rowMaxs(as.n)/rowSums(as.n)
    out.matrix[,i] <- p
    as <- as.n
  }
  return(out.matrix)
}


calc_armitage<- function(a, w){#where a is the matrix you want to run the test on, and w is the weights
  #separate controls
  control.cols <- grep("control", colnames(a))
  control<- a[,control.cols]
  #separate cases
  case <- a[,-control.cols]


  #identify ps qs and pqs
  homs <- substr(colnames(control), 1, 1) == substr(colnames(control), 2, 2)
  hom_controls <- control[,which(homs)]
  het_controls <- control[,which(!homs)]

  homs <- substr(colnames(case), 1, 1) == substr(colnames(case), 2, 2)
  hom_cases <- case[,which(homs)]
  het_cases <- case[,which(!homs)]

  pp_cont <- matrixStats::rowMaxs(hom_controls)
  pp_case <- matrixStats::rowMaxs(hom_cases)
  qq_cont <- rowSums(hom_controls) - pp_cont
  qq_case <- rowSums(hom_cases) - pp_case
  pq_cont <- rowSums(het_controls)
  pq_case <- rowSums(het_cases)


  # define values for the three sums
  n <- rbind(pp_case, pq_case, qq_case)
  s1 <- n*w
  s1 <- colSums(s1)

  N <- n + rbind(pp_cont, pq_cont, qq_cont)
  s2 <- N*w
  s2 <- colSums(s2)
  s3 <- N*w^2
  s3 <- colSums(s3)

  bT <- colSums(N)
  t <- colSums(n)

  # equation 5
  b <- (bT*s1 - t*s2)/(bT*s3 - (s2^2))

  # equation 6
  Vb <- (t*(bT - t))/(bT*(bT*s3 - s2^2))

  # equation 7
  chi <- (b^2)/Vb

  # use the pchisq function with 1 degree of freedom to get a p-value and return.
  return(pchisq(chi, 1, lower.tail = F))
}


# process:
# Take the current guesses at haplotypes to estimate the true frequencies
# use these updated frequencies to make new guesses at haplotypes (doesn't have to be 1 or 0, can be .75 of one haplotype)
# repeat
estimate_haplotype_frequencies <- function(x){
  #function to correct haplotype input matrix:
  GtoH <- function(x, n){
    m1 <- matrix(as.numeric(x), nrow(x), ncol(x))
    colnames(m1) <- colnames(x)
    m1 <- cbind(as.data.frame(t(m1)), n)
    m2 <- m1 %>% group_by(n) %>% summarise_all(funs(sum))
    m2 <- t(as.matrix(m2[,-1]))
    colnames(m2) <- sort(unique(n))
    return(m2)
  }

  browser()
  #1)
  #get the observed genotype combinations
  yv <- as.vector(t(y))
  gcv <- paste0(x, yv)
  if(length(y) == length(x)){
    gcv <- matrix(gcv, 1, length(x), byrow = T)
  }
  else{
    gcv <- matrix(gcv, nrow(y), length(x), byrow = T)
  }


  #2)
  #turn this into a genotype count table
  mgcv <- reshape2::melt(gcv)
  cnames <- levels(mgcv$value)
  ghapmat <- bigtabulate::bigtabulate(mgcv, ccols = which(colnames(mgcv) %in% c("Var1", "value")))
  colnames(ghapmat) <- cnames

  #3) clean the table
  ##grab column names
  gl <- colnames(ghapmat)

  ##remove anything with missing data and double hets
  rgcs <- c(grep(paste0("^", dmDat), gl), #missing first locus
            grep(paste0(dmDat, "$"), gl), #missing second locus
            which(substr(gl, 1, sform) != substr(gl, (sform + 1), (sform *2)) &
                    substr(gl, (sform*2) + 1, sform*3) != substr(gl, (sform*3+1), sform*4))) #double het

  ##remove any double heterozygotes
  ghapmat <- ghapmat[,-rgcs]

  #add a filler row for the last pairwise comparison to make life easier.
  if(length(y) == length(x)){
    if(length(ghapmat) > 1){ #stop it from doing this if there is data for only one haplotype.
      ghapmat <- rbind(ghapmat, rep(c(10,0), 100)[1:length(ghapmat)])
    }
  }

  #if nothing remains, return nothing
  if(length(ghapmat) == 0){
    return(NA)
  }
  else if(!is.matrix(ghapmat)){ #if only one column remains...
    if(substr(gl[-rgcs], 1,sform) == substr(gl[-rgcs], sform + 1, sform*2)){ #double hom
      return(NA)
    }
    else{
      return(c(D = 0))
    }
  }

  #4) get hap table. Use the rules above to do this. Possible conditions, where alleles at locus 1 = 1a1b and locus 2 = 2a2b
  ghapmat <- ghapmat[,order(colnames(ghapmat))] #put in order, just in case
  gl <- colnames(ghapmat) #get column names again
  hnames <- c(paste0(as, as[1]),
              paste0(as, as[2]),
              paste0(as, as[3]),
              paste0(as, as[4]))
  if(any(grepl("NA", hnames))){ #in the case of a missing allele...
    hnames <- hnames[-grep("NA",hnames)]
  }
  hnames <- sort(hnames)
  hapmat <- matrix(0, nrow(ghapmat), length(hnames))#initialize. all possible haplotypes
  colnames(hapmat) <- hnames


  ##figure out which are homozygotes and heterozygotes at either locus in the pairwise comparison.
  dhom <- substr(gl, 1, sform) == substr(gl, sform + 1, sform*2) &
    substr(gl, (sform*2) + 1, sform*3) == substr(gl, (sform*3+1), sform*4) #columns with double homozygotes
  dhom <- ghapmat[,dhom]
  het_l1 <- substr(gl, 1, sform) != substr(gl, sform + 1, sform*2) #columns where the first locus is het
  het_l1 <- ghapmat[,het_l1]
  het_l2 <- substr(gl, (sform*2) + 1, sform*3) != substr(gl, (sform*3+1), sform*4) #colunms where the second locus is het
  het_l2 <- ghapmat[,het_l2]

  #fix wierd cases where one of these isn't a matrix because only one haplotype falls into the category.
  if(any(!is.matrix(dhom), !is.matrix(het_l1), !is.matrix(het_l2))){
    if(!is.matrix(dhom)){
      dhom <- as.matrix(dhom)
      colnames(dhom) <- colnames(ghapmat)[substr(gl, 1, sform) == substr(gl, sform + 1, sform*2) &
                                            substr(gl, (sform*2) + 1, sform*3) == substr(gl, (sform*3+1), sform*4)] #columns with double homozygotes
    }
    if(!is.matrix(het_l1)){
      het_l1 <- as.matrix(het_l1)
      colnames(het_l1) <- colnames(ghapmat)[substr(gl, 1, sform) != substr(gl, sform + 1, sform*2)] #columns where the first locus is het
    }
    if(!is.matrix(het_l2)){
      het_l2 <- as.matrix(het_l2)
      colnames(het_l2) <- colnames(ghapmat)[substr(gl, (sform*2) + 1, sform*3) != substr(gl, (sform*3+1), sform*4)]
    }
  }

  #count up the haplotypes.
  ##homozygotes:
  ### unless there are no double homozygotes:
  if(sum(substr(gl, 1, sform) == substr(gl, sform + 1, sform*2) &
         substr(gl, (sform*2) + 1, sform*3) == substr(gl, (sform*3+1), sform*4)) != 0){
    hapmat[,colnames(hapmat) %in% paste0(substr(colnames(dhom), 1, sform),
                                         substr(colnames(dhom),(sform*2)+1,sform*3))] <- dhom*2
  }

  ##heterozyogote locus 1
  ### unless locus one has no heterozygotes:
  if(sum(substr(gl, 1, sform) != substr(gl, sform + 1, sform*2)) != 0){
    n1 <- paste0(substr(colnames(het_l1), 1, sform),
                 substr(colnames(het_l1),(sform*2)+1,sform*3))
    n1 <- GtoH(het_l1, n1)
    n2 <- paste0(substr(colnames(het_l1),sform+1, sform*2),
                 substr(colnames(het_l1),(sform*3)+1, sform*4))
    n2 <- GtoH(het_l1, n2)
    hapmat[,colnames(hapmat) %in% colnames(n1)] <- n1 + hapmat[,colnames(hapmat) %in% colnames(n1)]
    hapmat[,colnames(hapmat) %in% colnames(n2)] <- n2 + hapmat[,colnames(hapmat) %in% colnames(n2)]
  }


  ##heterozyogote locus 2
  ### unless locus two has no heterozygotes
  if(sum(substr(gl, (sform*2) + 1, sform*3) != substr(gl, (sform*3+1), sform*4)) != 0){
    n1 <- paste0(substr(colnames(het_l2), 1, sform),
                 substr(colnames(het_l2),(sform*2)+1,sform*3))
    n1 <- GtoH(het_l2, n1)
    n2 <- paste0(substr(colnames(het_l2),sform+1, sform*2),
                 substr(colnames(het_l2),(sform*3)+1, sform*4))
    n2 <- GtoH(het_l2, n2)
    hapmat[,colnames(hapmat) %in% colnames(n1)] <- n1 + hapmat[,colnames(hapmat) %in% colnames(n1)]
    hapmat[,colnames(hapmat) %in% colnames(n2)] <- n2 + hapmat[,colnames(hapmat) %in% colnames(n2)]
  }


  #5)condense this hap table into the 1a2a, 1a2b, 1b2a, 1b2b format.
  # figure out how where haplotypes are missing. Note, do the case of two or three
  # missin haplotypes at the end.
  pmat <- ifelse(hapmat == 0, F, T)
  mmat <- pmat
  l1 <- substr(colnames(pmat), 1, sform)
  l2 <- substr(colnames(pmat), sform + 1, sform*2)
  mc <- 4 - rowSums(pmat)

  #function to see if haplotype is missing. x is the row index, m is a vector of the number missing haplotypes at each locus.
  cmhap <- function(x){
    out <- ifelse(rowSums(pmat[,l1 == l1[x]]) == 0 | rowSums(pmat[,l1 == l1[x]]) == 2, F,
                  ifelse(rowSums(pmat[,l2 == l2[x]]) > 0 & pmat[,x] != TRUE, T, F))
    return(out)
  }
  #fill the mmat.
  for(i in 1:ncol(pmat)){
    mmat[,i] <- cmhap(i)
  }

  #set the missing values in hapmat to NA, then replace those with zeros where there
  #are missing haplotypes.
  hapmat[hapmat == 0] <- NA
  hapmat[mmat == TRUE] <- 0






  #put in fillers when there are more than one haplotype is missing.
  pmat <- ifelse(is.na(hapmat), F, T)
  missing <- 4 - rowSums(pmat)
  m2 <- ifelse(missing >= 2, 0, NA)
  m3 <- ifelse(missing >= 3, 0, NA)
  m4 <- ifelse(missing == 4, 0, NA)
  mc <- cbind(m2,m2,m3,m4)

  #figure out which D, r values to give if two are missing...
  if(any(missing == 2)){
    m2mat <- hapmat[missing == 2,]
    m2matv <- as.vector(t(m2mat))
    if(is.matrix(m2mat)){
      m2matvcn <- rep(colnames(m2mat), nrow(m2mat))
    }
    else{
      m2matvcn <- rep(names(m2mat), 1)
    }
    m2matv[!is.na(m2matv)] <- m2matvcn[!is.na(m2matv)]
    m2matv <- na.omit(m2matv)
    if(is.matrix(m2mat)){
      m2mat <- matrix(m2matv, nrow(m2mat), 2, T)
      m2mat <- ifelse(substr(m2mat[,1], 1, sform) != substr(m2mat[,2], 1, sform) &
                        substr(m2mat[,1], sform + 1, sform*2) != substr(m2mat[,2], sform + 1, sform*2),
                      1,0)
    }
    else{
      m2mat <- ifelse(substr(m2matv[1], 1, sform) != substr(m2matv[2], 1, sform) &
                        substr(m2matv[1], sform + 1, sform*2) != substr(m2matv[2], sform + 1, sform*2),
                      1,0)
    }

  }
  else{
    m2mat <- "none"
  }

  hapmat <- cbind(hapmat, mc)
  hapmat <- as.vector(t(hapmat))
  hapmat <- na.omit(hapmat)
  #if(length(hapmat) %% nrow(y) != 0){
  #  browser()
  #}
  hapmat <- matrix(hapmat, nrow(ghapmat), 4, byrow = T)

  #now just have the haplotypes. These will calculate D in the case of 1 or 0 missing haplotypes.
  #when there are three missing haplotypes, D will be 0. When there are 2, D will be 0 or 1.
  return(list(hapmat = hapmat, missing = missing, m2 = m2mat))




































  # browser()
  # phenos <- which(x != 0)
  # x <- x[phenos]
  # # cj values for each possible genotype:
  # s1 <- substr(names(x), 1, 1)
  # s2 <- substr(names(x), 2, 2)
  # s3 <- substr(names(x), 3, 3)
  # s4 <- substr(names(x), 4, 4)
  #
  #
  #
  # het.1 <- s1 != s2
  # het.2 <- s3 != s4
  # sj <- het.1 + het.2
  # pair1 <- paste0(s1, s3)
  # pair2 <- paste0(s2, s4)
  # pair3 <- paste0(s1, s4)
  # pair4 <- paste0(s2, s3)
  # haplopair1 <- cbind(pair1, pair2)
  # haplopair2 <- cbind(pair3, pair4)
  #
  # haplopairs <- array(dim = c(nrow(haplopair1), 2, 2))
  # haplopairs[,,1] <- haplopair1
  # haplopairs[,,2] <- haplopair2
  #
  # # equation 2
  # cj <- 2^(sj - 1)
  # cj[sj == 0] <- 1
  #
  # # equation 3, starting prob of each possible genotype contributing to each observed phenotype
  # P0 <- 1/cj
  # P0 <- P0

  # only unkown are double hets. Everything else has a known state.

  # get.pj <- function(t, j, i, P0, hap.freqs, rep){
  #   tpair <- haplopairs[j,,i] # haplopair for phenotype j, option i (there will only be two DIFFERENT options for double hets)
  #   delta <- sum(tpair == haps[t])
  #   if(delta == 0){return(c(NA, delta))}
  #
  #   # expectation
  #   # nj <- x[j]
  #   # n <- sum(x)
  #   if(rep == 0){
  #     pj <- P0[j]
  #   }
  #   else{
  #     hap.matches <- match(tpair, haps)
  #     if(tpair[1] == tpair[2]){
  #       pj <- hap.freqs[hap.matches[1]]^2
  #     }
  #     else{
  #       pj <- 2*prod(hap.freqs[hap.matches])
  #     }
  #   }
  #   return(c(pj, delta))
  # }
  #
  # rep <- 0
  #
  # # expectation
  # # Pj is the frequency of phenotype j
  # # Pj(hkhl) is the probability of drawing the phenotype pair (p^2 or 2pq)
  # # nj is the count of the jth phenotype
  # # n is the total count
  #
  # # maximization
  # hap.freqs <- numeric(length(haps))
  # for(t in 1:nhap){
  #   browser()
  #   tcount <- 0
  #   for(j in 1:m){
  #     gpsj <- numeric(0)
  #     deltas <- numeric(0)
  #     for(i in 1:cj[m]){
  #       pjd <- get.pj(t, j, i, P0, hap.freqs, rep)
  #       gpsj <- c(gpsj, pjd[1])
  #       deltas <- c(deltas, pjd[2])
  #       # tcount <- tcount + delta*tp
  #     }
  #     pjs <- gpsj/sum(gpsj, na.rm = T)
  #     tcount <- tcount + sum(pjs*deltas, na.rm = T)
  #
  #
  #   }
  #   tcount <- .5*tcount
  #   n.freqs[t] <- tcount
  # }




}
