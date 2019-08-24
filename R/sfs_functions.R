calcSFS <- function(x, facet = "pop", pops, projection, fold = T){

  # get dadi formatted data
  y <- format_snps(x, output = "dadi", facets = facet)

  # get sfs
  sfs <- make_sfs(y, pops, projection, fold)

  return(sfs)
}

# function to get an sfs contribution for a single snp across 2 populations
# takes a list of length 2, where the each element is a numeric vector containing:
# 1) the number of sequenced alleles at the site
# 2) the number of derived alleles at the site
# 3) the polarization status
#
# also needs a vector of projection sizes for the populations
find_contribs <- function(pop_counts, projection){
  pcontribs <- vector("list", length(pop_counts))

  for(i in 1:length(pop_counts)){
    # pull info
    hits <- pop_counts[[i]][2]
    n <- pop_counts[[i]][1]
    m <- projection[i]


    # if there are less snps sequenced here than our projection, it gets zero
    # for all contributions
    if(n < m){
      contrib <- rep(0, m + 1)
    }

    else{
      # make the hits vector
      proj_hits <- 0:m

      # do the binomials
      contrib <- choose(m, proj_hits)
      contrib <- contrib*choose(n - m , hits - proj_hits)
      contrib <- contrib/choose(n, hits)
    }
    pcontribs[[i]] <- contrib
  }

  return(Reduce(outer, pcontribs))
}

# function to get a pop_counts array for projections. Result has three columns,
# containing the three bits of info needed for a pop_counts list. The third
# index lists the pop
# input is a dadi formatted data.frame
get_counts <- function(x, pops){
  # figure out which allele is the outgroup/derived in each row and figure out polarization status
  anc.als <- substr(x$anc, 2, 2)
  polarized <- rep(T, nrow(x))
  polarized[anc.als == "N"] <- F
  anc.als[anc.als == "N"] <- x$Allele1[anc.als == "N"]
  which.derived <- rep(2, nrow(x))
  which.derived[anc.als == x$Allele2] <- 1

  # initialize output
  out <- array(0, dim = c(nrow(x), 3, 2))

  # write to output
  for(i in 1:length(pops)){
    dat.cols <- grep(pops[i], colnames(x))
    tdat <- x[,dat.cols]
    tdat <- as.matrix(tdat)
    out[,1,i] <- rowSums(tdat) # total count
    t.index <- 2 * (1:nrow(tdat))
    t.index[which.derived == 1] <- t.index[which.derived == 1] - 1
    out[,2,i] <- t(tdat)[t.index] # derived allele count
    out[,3,i] <- polarized
  }

  return(out)
}

# make a projected sfs from a counts array
make_proj_sfs <- function(counts, projection, fold = F){

  # initialize
  if(length(projection) == 2){
    sfs <- matrix(0, projection[1] + 1, projection[2] + 1)
  }
  else{
    sfs <- numeric(projection + 1)
    counts[,,2] <- 1
  }

  # project each snp
  for(i in 1:nrow(counts[,,1])){

    # if the snp isn't sequenced in any pop, skip
    if(any(counts[i,1,] == 0)){
      next
    }

    # if we aren't folding and this snp isn't polarizable (usually meaning we had NNN for the anc condition),
    # skip it
    if(!fold & !counts[i,3,1]){
      next
    }

    if(length(projection) == 2){
      pop_counts <- list(counts[i,,1], counts[i,,2])
    }
    else{
      pop_counts <- list(counts[i,,1])
    }
    contrib <- find_contribs(pop_counts, projection)
    sfs <- sfs + contrib
  }

  return(sfs)
}

# function to fold an sfs
fold_sfs <- function(sfs){

  # fold a 2 dimensional sfs
  if(length(dim(sfs)) == 2){
    rev.matrix <- function(x){
      return(x[nrow(x):1, ncol(x):1])
    }

    # figure out the part that gets folded in
    sample_sizes <- matrix(0:(nrow(sfs) - 1), nrow(sfs), ncol(sfs))
    sample_sizes <- sample_sizes +
      t(matrix(0:(ncol(sfs) - 1), ncol(sfs), nrow(sfs)))

    total.samples <- (ncol(sfs) + nrow(sfs) - 2)
    folded.in <- ifelse(sample_sizes > floor(total.samples/2), T, F)

    # reverse the sfs
    reversed <- sfs
    reversed[folded.in == F] <- 0
    reversed <- rev.matrix(reversed)

    # get the folded sfs
    folded <- sfs + reversed

    # deal with ambiguous areas
    ambig <- ifelse(sample_sizes == total.samples/2, T, F)
    ambig[ambig == T] <- sfs[ambig == T]
    rev.ambig <- rev.matrix(ambig)
    folded <- folded + (-.5*ambig + .5*rev.ambig)

    folded[folded.in == T] <- NA
  }

  # fold a 1 dimensional sfs--life is easy
  else{
    rev <- rev(sfs)[1:floor(length(sfs)/2)]
    rev <- c(rev, rep(0, length(sfs) - length(rev)))
    folded <- sfs + rev
    folded[ceiling(length(sfs)/2):length(sfs)] <- NA
  }
  return(folded)
}

# wrapper
make_sfs <- function(x, pops, projection, fold = F){
  counts <- get_counts(x, pops)
  sfs <- make_proj_sfs(counts, projection, fold)
  if(fold){
    sfs <- fold_sfs(sfs)
  }

  # add a pops attribute
  attr(sfs, which = "pop") <- pops


  cat("SFS completed with", sum(sfs, na.rm = T), "segrgating sites.\n")
  return(sfs)
}

# plot
plot_sfs <- function(sfs, viridis.option = "inferno", log = TRUE){

  # add colun names, row names, and melt
  pops <- attr(sfs, "pop")
  sfs <- as.data.frame(sfs)
  if(length(pops) == 1){
    sfs <- as.data.frame(t(sfs))
  }
  colnames(sfs) <- 0:(ncol(sfs) - 1)
  sfs$count <- 0:(nrow(sfs) - 1)
  sfs[1,1] <- NA # mask the first entry
  msfs <- reshape2::melt(sfs, id.vars = "count")
  colnames(msfs) <- c("p1", "p2", "N")
  msfs$p1 <- as.integer(msfs$p1)
  msfs$p2 <- as.integer(msfs$p2)



  # plot
  ## 2D
  if(nrow(sfs) > 1){
    if(log){
      p <- ggplot2::ggplot(msfs, ggplot2::aes(x = p1, y = p2, fill = log10(N), color = log10(N)))
    }
    else{
      p <- ggplot2::ggplot(msfs, ggplot2::aes(x = p1, y = p2, fill = N, color = N))
    }

    p <- p +
      ggplot2::geom_tile() + ggplot2::theme_bw() +
      ggplot2::scale_color_viridis_c(na.value = "white", option = viridis.option) +
      ggplot2::scale_fill_viridis_c(na.value = "white", option = viridis.option) +
      ggplot2::xlab(pops[1]) + ggplot2::ylab(pops[2]) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0))
  }

  ## 1D
  else{
    if(log){
      p <- ggplot2::ggplot(msfs, ggplot2::aes(x = p2, y = log(N)))
    }
    else{
      p <- ggplot2::ggplot(msfs, ggplot2::aes(x = p2, y= N))
    }
    p <- p +
      ggplot2::geom_line() + ggplot2::theme_bw() +
      ggplot2::xlab(pops[1])
      ggplot2::scale_x_continuous(expand = c(0, 0))
  }


  return(p)
}

calc_directionality <- function(sfs){
  # normalize, exluding fixed sites
  sfs <- sfs[-1,-1] # remove fixed sites
  sfs <- sfs/sum(sfs)

  # get all of the allele frequencies in each cell
  freqs.i <- (1:(nrow(sfs)))/(nrow(sfs))
  freqs.j <- (1:(ncol(sfs)))/(ncol(sfs))

  # do the math, this is equ 1b (i - j times sfs)
  directionality <- sum(outer(freqs.i, freqs.j, "-") * sfs)

  return(directionality)
}
