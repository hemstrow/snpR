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
