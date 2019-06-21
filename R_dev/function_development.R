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

  # a function to do random mating in one generation
  random.mating <- function(as, N){
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
    vio.cols <- which(colSums(logi.matrix) != 1) # need to note any cols where frequency of p = .5, since that may mess things up.
    clean.logi <- logi.matrix[,-vio.cols] # remove vio cols
    major.slots <- which(clean.logi) # these cells are the major alleles
    new.matrix[,-vio.cols][major.slots] <- new.ps[-vio.cols] # assign to the output


    # do the same thing for the minor alleles. Note that there may be columns with no minor alleles, so we need to adjust for that as well.
    min.matrix <- matrix(rep(rowSums(as) - p.count, each = 4), ncol = ncol(tas))
    fixed.cols <- which(colSums(min.matrix) == 0) # cols with no minor allele
    min.matrix[,fixed.cols] <- -1 # set these to -1, so that can't possibly match
    logi.matrix.min <- ifelse(tas == min.matrix, T, F)
    clean.logi.min <- logi.matrix.min[,-vio.cols]
    minor.slots <- which(clean.logi.min) # cells with minor alleles
    new.matrix[,-vio.cols][minor.slots] <- new.qs[-sort(c(vio.cols, fixed.cols))] # assign, note that we need to account for both fixed and violating rows now.

    # now need to fix the violating cols.
    ## make the replacement matrix
    vio.matrix <- matrix(0, ncol = length(vio.cols), nrow = 2)
    vio.matrix[1,] <- new.ps[vio.cols]
    vio.matrix[2,] <- new.qs[vio.cols]

    ## replace
    no.dat.cols <- which(colSums(tas[,vio.cols]) == 0)
    replace.slots <- which(ifelse(tas[,vio.cols][,-no.dat.cols] != 0, T, F)) # the slots to replace, where the as data actually has allele counts.
    new.matrix[,vio.cols][,-no.dat.cols][replace.slots] <- vio.matrix[,-no.dat.cols]
    new.matrix[is.na(new.matrix)] <- 0

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


# caluclate a one dimensional sfs. Not fully vectorized, loops through bins.
do.afs <- function(as, bins){
  # get p and q
  p <- matrixStats::rowMaxs(as)/rowSums(as)
  q <- 1 - p

  # make bins
  lbin <- seq(0, .5, length.out = bins + 1)

  # do a for loop
  tracker <- 1
  my.mat <- matrix(NA, 1, bins)
  for(i in lbin[-length(lbin)]){
    out <- q > i & q < lbin[tracker + 1]

    num <- sum(out, na.rm = T)

    my.mat[1, tracker] <- num

    tracker <- tracker + 1
  }

  binnames <- lbin[-length(lbin)]
  colnames(my.mat) <- binnames

  plot.dat <- as.data.frame(t(my.mat))
  plot.dat$bin <- colnames(my.mat)
  colnames(plot.dat)[1] <- "snps"

  plot <- ggplot(plot.dat, aes(x = bin, y = snps)) + geom_col() + theme_bw()

  print(plot)

  return(my.mat)
}







# x is a haplotype table, like in 'R_dev/haplotype_frequency_estimation_example.RDS'
# haptable is a haplotype identity table (00, 01, 10, 11), like in "R_dev/hapmat_example.RDS"
# sigma is the ML difference cuttoff
multi_haplotype_estimation <- function(x, haptable, sigma = 0.0001){

  # find the double het. Should be able to use an approach like this when this gets extended to work with everything.
  # cj values for each possible genotype:
  s1 <- substr(colnames(x), 1, 1)
  s2 <- substr(colnames(x), 2, 2)
  s3 <- substr(colnames(x), 3, 3)
  s4 <- substr(colnames(x), 4, 4)
  het.1 <- s1 != s2
  het.2 <- s3 != s4



  # First, make a guess at the starting haplotype frequencies. We'll do this by taking the unambigious haplotype frequencies,
  # then making a guess at the haplotype composition in the double heterozygote assuming that all possible haplotypes are equally likely
  doub.het <- which(het.1 + het.2 == 2) # identify double heterozygotes

  # if there are no double heterozygotes, we already no the haplotype frequencies, so we can just return those.
  if(length(doub.het) == 0){
    return(haptable/rowSums(haptable))
  }

  nhap.counts <- haptable # grab the haplotypes
  ehap.counts <- nhap.counts + .5*rowSums(x[,doub.het]) # assuming that both haplopairs are equaly likely in the double het
  shap.freqs <- ehap.counts/rowSums(ehap.counts) # get the starting haplotype frequencies



  # now that we have our starting conditions, we will do the EM loops.
  # 1) First, we find out how many of each haplotype
  # we expect to get from our double heterozygotes given the initial haplotype frequencies we guessed above.
  # 2) Then, we use those expected frequencies to update our estimates of the haplotype frequencies.
  # 3) repeat 1 and 2 until the difference between the haplotype frequencies between loop iterations is less than sigma.


  # we'll use a while loop, which will run as long as the conditions are met. Note that this can freeze your computer if
  # the conditions are NEVER met! Try just hitting the stop sign, if that doesn't work you'll need to restart Rstudio.

  out <- matrix(NA, nrow(haptable), ncol = 4)
  completed <- logical(nrow(hapmat))

  diff <- sigma + 1 # initialize the diff. Doesn't matter what it is as long as it's larger than sigma.
  while(any(diff > sigma)){

    if(is.null(nrow(shap.freqs))){
      shap.freqs <- matrix(shap.freqs, 1)
      x <- matrix(x, 1)
      haptable <- matrix(haptable, 1)
    }

    # 1)
    # expectation, which is that we are drawing haplotypes (aka alleles) from a pool of options. Follows HWE, essentially,
    # but the "alleles" are actually haplotypes
    op1.e <- (2*shap.freqs[,1]*shap.freqs[,4])/
      ((2*shap.freqs[,1]*shap.freqs[,4])+(2*shap.freqs[,2]*shap.freqs[,3])) # percentage of AC/GG haplo pairs
    op2.e <- 1 - op1.e

    # maximization: given the expected haplotype frequencies, how many of each haplotype should we have? get new frequencies
    n1hap.freqs <- haptable # grab the known haplotype frequencies form the unambigious phenotypes again.
    if(nrow(x) == 1){
      n1hap.freqs[,c(1, 4)] <- n1hap.freqs[,c(1, 4)] + (rowSums(matrix(x[,doub.het], 1))*op1.e*.5) # we basically add the expected number of haplotypes for the double heterozygotes
      n1hap.freqs[,c(2, 3)] <- n1hap.freqs[,c(2, 3)] + (rowSums(matrix(x[,doub.het], 1))*op2.e*.5)
    }
    else{
      n1hap.freqs[,c(1, 4)] <- n1hap.freqs[,c(1, 4)] + (rowSums(x[,doub.het])*op1.e*.5) # we basically add the expected number of haplotypes for the double heterozygotes
      n1hap.freqs[,c(2, 3)] <- n1hap.freqs[,c(2, 3)] + (rowSums(x[,doub.het])*op2.e*.5)
    }
    n1hap.freqs <- n1hap.freqs/rowSums(n1hap.freqs)

    # calculate the diff and update
    diff <- rowSums(abs(n1hap.freqs - shap.freqs))

    # save any completed results
    done <- which(diff <= sigma)
    if(length(done) > 0){
      if(sum(completed) > 0){
        if(nrow(n1hap.freqs) > 1){
          out[-which(completed),][done,] <- n1hap.freqs[done,]
        }
        else{
          out[-which(completed),] <- n1hap.freqs[done,]
        }
        completed[-which(completed)][done] <- T
      }
      else{
        out[done,] <- n1hap.freqs[done,]
        completed[done] <- T
      }
      n1hap.freqs <- n1hap.freqs[-done,]
      shap.freqs <- n1hap.freqs
      haptable <- haptable[-done,]
      x <- x[-done,]
    }
    else{
      shap.freqs <- n1hap.freqs
    }
  }


  # return the output
  return(out)
}

