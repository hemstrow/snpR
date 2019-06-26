#' Interface with BGLR to run genomic prediction with snpRdata objects.
#'
#' Run genomic prediction given a single responce variable (usually a phenotype)
#' using the \code{\link[BGLR]{BGLR}} function. Unlike other snpR functions,
#' this returns the resulting model direcly, so overwrite with caution.
#'
#' This function is provided as a wrapper to plug snpRdata objects into the
#' \code{\link[BGLR]{BGLR}} function in order to easily run genomic prediction
#' on a simple model where a single, sample specific meta data variable is
#' provided as the response variable. To do so, this function formats the data
#' into a transposed "sn" format, as described in \code{\link{format_snps}}
#' using the bernoulli method to interpolate missing genotypes. Several
#' different prediction models are available, see the documentation the ETA
#' argument in \code{\link[BGLR]{BGLR}} for details. Defaults to the "BayesB"
#' model, which assumes a "spike-slab" prior for allele effects on phenotype
#' where most markers have a very small effect size and a few can have a much
#' larger effect.
#'
#' Unlike most snpR functions, this function does not support facets, since each
#' run can be very slow. Instead, an individual facet and facet level of
#' interest should be selected with \code{\link{subset.snpR.data}}. See
#' examples.
#'
#' See documentation for \code{\link[BGLR]{BGLR}} for more details and for a
#' full list of references.
#'
#' @param x snpRdata object
#' @param responce character. Name of the column containing the responce
#'   variable of interest. Must match a column name in sample metadata.
#' @param iterations numeric. Number of iterations to run the MCMC chain for.
#' @param burn_in numeric. Number of burn in iterations to run prior to the MCMC
#'   chain.
#' @param thin numeric. Number of iterations to discard between each recorded
#'   data point.
#' @param model character, default "BayesB". Prediction model to use, see
#'   description for the ETA argument in \code{\link[BGLR]{BGLR}}.
#'
#' @export
#' @author William Hemstrom
#' @references Pérez, P., and de los Campos, G. (2014). \emph{Genetics}.
#' @return A list containing: \itemize{ \item{model: } The model output from
#'   BGLR. See \code{\link[BGLR]{BGLR}}. \item{h2: } Estimated heritability of
#'   the responce variable. \item{transposed_interpolated_genotypes: }
#'   Transposed genotypes, the input for BGLR. \item{predictions: } A data.frame
#'   containing the provided phenotypes and the predicted Breeding Values (BVs)
#'   for those phenotypes. }
#'
#' @examples
#' # run and plot a basic prediction
#' ## add some dummy phenotypic data.
#' sample.meta <- cbind(weight = rnorm(ncol(dat)), stickSNPs@sample.meta)
#' dat <- import.snpR.data(as.data.frame(stickSNPs), stickSNPs@snp.meta, sample.meta)
#' ## run prediction
#' gp <- run_genomic_prediction(dat, responce = "weight", iterations = 1000, burn_in = 100, thin = 10)
#' ## dummy phenotypes vs. predicted Breeding Values for those dummy predictions.
#' with(gp$predictions, plot(phenotype, predicted_BV))
#'
run_genomic_prediction <- function(x, responce, iterations,
                                   burn_in, thin,
                                   model = "BayesB"){

  # get a properly formatted input file
  sn <- format_snps(x, "sn")
  sn <- sn[,-c(1:(ncol(x@snp.meta) - 1))]

  # prepare the BGLR input
  sn <- t(sn)
  phenotypes <- x@sample.meta[,responce]
  colnames(sn) <- paste0("m", 1:ncol(sn)) # marker names
  rownames(sn) <- paste0("s", 1:nrow(sn)) # ind IDS
  ETA <- list(list(X = sn, model = "BayesB", saveEffects = T)) # need to adjust this when I get around to allowing for more complicated models

  BGLR_mod <- BGLR::BGLR(y = phenotypes, ETA = ETA, nIter = iterations + burn_in, burnIn = burn_in, thin = thin)

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

  pdat <- data.frame(phenotype = phenotypes, predicted_BV = sn%*%BGLR_mod$ETA[[1]]$b)

  return(list(model = BGLR_mod, h2 = h2, transposed_interpolated_genotypes = sn, predictions = pdat))
}


#' Run a single cross-validation with BGLR.
#'
#' Run genomic prediction given a single responce variable (usually a phenotype)
#' using the \code{\link[BGLR]{BGLR}} function. Unlike other snpR functions,
#' this returns the resulting model direcly, so overwrite with caution. This
#' will leave a specified portion of the samples out when running the model,
#' then perform cross-validation.
#'
#' This function is provided as a wrapper to plug snpRdata objects into the
#' \code{\link[BGLR]{BGLR}} function in order to easily run genomic prediction
#' on a simple model where a single, sample specific meta data variable is
#' provided as the response variable. To do so, this function formats the data
#' into a transposed "sn" format, as described in \code{\link{format_snps}}
#' using the bernoulli method to interpolate missing genotypes. Several
#' different prediction models are available, see the documentation the ETA
#' argument in \code{\link[BGLR]{BGLR}} for details. Defaults to the "BayesB"
#' model, which assumes a "spike-slab" prior for allele effects on phenotype
#' where most markers have a very small effect size and a few can have a much
#' larger effect.
#'
#' Unlike most snpR functions, this function does not support facets, since each
#' run can be very slow. Instead, an individual facet and facet level of
#' interest should be selected with \code{\link{subset.snpR.data}}. See
#' examples.
#'
#' The portion of samples left out for cross-validation is defined by the
#' cross_percentage argument. Specifically, the proportion given will be used to
#' create the model, the remaining portion will be left out. For example, if
#' cross_percentage = 0.9, 90% of the samples will be used to run the model and
#' the remaining 10% will be used for cross-validation. Alternatively, if
#' cross_samples are provided, the specifed samples (by column index) will be
#' used used for cross validation instead. This may be useful for a systematic
#' leave-one-out cross-validation.
#'
#'
#' See documentation for \code{\link[BGLR]{BGLR}} for more details and for a
#' full list of references.
#'
#' @param x snpRdata object
#' @param responce character. Name of the column containing the responce
#'   variable of interest. Must match a column name in sample metadata.
#' @param iterations numeric. Number of iterations to run the MCMC chain for.
#' @param burn_in numeric. Number of burn in iterations to run prior to the MCMC
#'   chain.
#' @param thin numeric. Number of iterations to discard between each recorded
#'   data point.
#' @param cross_percentage numeric, default 0.9. The proportion of sample to use
#'   to create the model. Must be greater than 0 and less than 1.
#' @param model character, default "BayesB". Prediction model to use, see
#'   description for the ETA argument in \code{\link[BGLR]{BGLR}}.
#' @param cross_samples numeric, default NULL. Optional vector of sample indices
#'   to use for cross-validation. If provided, the cross_percentage argument
#'   will be ignored.
#' @param plot logical, default TRUE. If TRUE, will generate a ggplot of the
#'   cross-validation results.
#'
#' @export
#' @author William Hemstrom
#' @references Pérez, P., and de los Campos, G. (2014). \emph{Genetics}.
#' @return A list containing: \itemize{ \item{model: } The results from
#'   \code{\link{run_genomic_prediction}}. \item{model.samples: } Indices of the
#'   samples used to construct the model. \item{cross.samples: } Indices of the
#'   samples used to cross-validate the model. \item{comparison: } A data.frame
#'   containing the observed and predicted phenotypes/Breeding Values for the
#'   cross-validation samples. \item{rsq: } The r^2 value for the observed and
#'   predicted phenotypes/Breeding Values for the cross-validation samples. }
#'
#' @examples
#' # run and plot a basic prediction
#' ## add some dummy phenotypic data.
#' sample.meta <- cbind(weight = rnorm(ncol(dat)), stickSNPs@sample.meta)
#' dat <- import.snpR.data(as.data.frame(stickSNPs), stickSNPs@snp.meta, sample.meta)
#' ## run cross_validation
#' cross_validate_genomic_prediction(dat, responce = "weight", iterations = 1000, burn_in = 100, thin = 10)
#'
cross_validate_genomic_prediction <- function(x, responce, iterations,
                                              burn_in, thin, cross_percentage = 0.9,
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
  model <- run_genomic_prediction(sub.x, responce = responce, iterations =  iterations,
                                  burn_in = burn_in, thin = thin, model = model)

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
