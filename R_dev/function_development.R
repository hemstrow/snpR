# odds ratio association testing wrapper. Needs internal and apply.snpR.dat support.
calc_odds_ratio <- function(x, facet){
  func <- function(cast_ac){
    a<- cast_ac[,1]
    b<- cast_ac[,2]
    c<- cast_ac[,3]
    d<- cast_ac[,4]
    odds<- (a/b)/(c/d)
    return(odds)
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
