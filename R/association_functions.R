run_genomic_prediction <- function(x, responce, iterations,
                                   burn_in, thin, formula = NULL,
                                   model = "BayesB"){
  # get a properly formatted input file
  sn <- format_snps(x, "sn")
  sn <- sn[,-c(1:(ncol(x@snp.meta) - 1))]

  # run the interpolation
  cat("Interpolating NAs.\n")
  sn <- interpolate_sn(sn)

  # prepare the BGLR input
  sn <- t(sn)
  phenotypes <- x@sample.meta[,responce]
  colnames(sn) <- paste0("m", 1:ncol(sn)) # marker names
  rownames(sn) <- paste0("s", 1:nrow(sn)) # ind IDS
  ETA <- list(list(X = sn, model = "BayesB", saveEffects = T)) # need to adjust this when I get around to allowing for more complicated models


  BGLR_mod <- BGLR::BGLR(y = phenotypes, ETA = ETA, nIter = iterations, burnIn = burn_in, thin = thin)

  # grab h2 estimate
  B <- BGLR::readBinMat('ETA_1_b.bin')
  h2 <- rep(NA,nrow(B))
  varU <- h2
  varE <- h2
  for(i in 1:length(h2)){
    u <- sn%*%B[i,]
    varU[i] <- var(u)
    varE[i] <- var(phenotypes-u)
    h2[i] <- varU[i]/(varU[i] + varE[i])
  }
  h2 <- mean(h2)
  return(list(model = BGLR_mod, h2 = h2))
}
