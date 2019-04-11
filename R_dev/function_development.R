# HWE chi-squared p-values per snp function part. Works on genotype matrix.
HWE_part <- function(gs){
  browser()

  # get observed genotype frequencies
  het.col <- which(substr(colnames(gs), 1, 1) != substr(colnames(gs), 2, 2))
  o2pq <- rowSums(genotypes[,het.col])/rowSums(genotypes)
  opp <- matrixStats::rowMaxs(hom.cols)/rowSums(genotypes)
  oqq <- 1 - o2pq - opp

  # get allele frequencies
  fp <- opp + o2pq/2
  fq <- 1 - fp

  # get expected genotype frequencies
  epp <- fp^2
  eqq <- fq^2
  e2pq <- 2*fp*fq

  # calculate chi-squared
  calc.chi <- function(o,e){
    return(((o-e)^2)/e)
  }
  chi.pp <- calc.chi(opp, epp)
  chi.qq <- calc.chi(oqq, eqq)
  chi.2pq <- calc.chi(o2pq, e2pq)
  chi <- chi.pp + chi.qq + chi.2pq

  # calculate p-values
  out <- pchisq(chi, 2, lower.tail = FALSE)
  return(out)
}

# odds ratio association testing wrapper. Needs internal and apply.snpR.dat support.
calc_odds_ratio <- function(x, facet){
  func <- function(cast_ac){}


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
