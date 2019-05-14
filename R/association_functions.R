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
  return(list(model = BGLR_mod, h2 = h2, transposed_interpolated_genotypes = t(sn)))
}

cross_validate_genomic_prediction <- function(x, responce, iterations,
                                              burn_in, thin, cross_percentage = 0.9, formula = NULL,
                                              model = "BayesB", cross_samples = NULL, plot = TRUE){

  # pick samples to make the model with
  if(is.numeric(cross_samples)){
    msg <- "Provided cross_samples must be sample indices.\n"
    if(as.integer(cross_samples) != cross_samples){
      stop(msg)
    }
    if(any(cross_samples > ncol(x)) | any(cross_samples < 0)){
      stop(msg)
    }
    model.samples <- (1:ncol(x))[-cross_samples]
  }
  else{
    model.samples <- sample(ncol(x), floor(ncol(x)*(cross_percentage)), replace = F)
  }


  # run the model
  capture.output(suppressWarnings(suppressMessages(sub.x <- subset.snpR.data(x, samps = model.samples))))
  model <- run_genomic_prediction(sub.x, responce, iterations,
                                  burn_in, thin, formula, model)

  # check accuracy
  cross.samples <- (1:ncol(x))[-sort(model.samples)]
  whole.sn <- format_snps(x, "sn")
  whole.sn <- whole.sn[,-(1:(ncol(x@snp.meta) - 1))]
  whole.sn <- interpolate_sn(whole.sn)
  new.sn <- whole.sn[,cross.samples]
  new.sn <- t(new.sn)
  pred.phenos <- new.sn%*%model$model$ETA[[1]]$b
  pdat <- as.data.frame(cbind(observed = x@sample.meta[cross.samples, responce], predicted = pred.phenos), stringsAsFactors = F)
  colnames(pdat) <- c("observed", "predicted")

  # plot
  if(plot == TRUE){
    tplot <- ggplot2::ggplot(as.data.frame(pdat), ggplot2::aes(observed, predicted)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "lm") +
      ggplot2::theme_bw()

    print(tplot)
  }

  return(list(model = model, model.samples = model.samples, cross.samples = cross.samples, comparison = pdat,
              rsq = summary(lm(observed~predicted, pdat))$r.squared))

}
