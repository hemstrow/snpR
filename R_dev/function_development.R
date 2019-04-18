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
  tac <- tac[-which(colnames(tac) %in% c("facet", "snp", "position", "group", "facet.type"))]


  ctac <- data.table::dcast(data.table::setDT(tac), .snp.id ~ subfacet, value.var = c("ni1", "ni2"))
  ctac <- as.matrix(ctac)[,-1]

}


do.drift <- function(as, gen, N){

  # a function to do random mating in one generation
  random.mating <- function(as, N){
    # figure out how many As, Cs, Gs, and Ts each individual draws
    new.matrix <- matrix(0, nrow = nrow(as), ncol = ncol(as))
    colnames(new.matrix) <- colnames(as)
    for(i in 1:nrow(as)){
      old.alleles <- as[i,]
      new.alleles <- rmultinom(N, 2, prob = old.alleles)
      new.frequencies <- rowSums(new.alleles)
      new.matrix[i,] <- new.frequencies
    }
    return(new.matrix)
  }

  # do random mating for gen generations and save the frequencies at each step
  out.matrix <- matrix(0, nrow = nrow(as), ncol = gen)
  for(i in 1:gen){
    as.n <- random.mating(as, N)
    p <- matrixStats::rowMaxs(as.n)/rowSums(as.n)
    out.matrix[,i] <- p
    as <- as.n
  }
  return(out.matrix)
}
