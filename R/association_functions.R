#' Interface with BGLR to run genomic prediction with snpRdata objects.
#'
#' Run genomic prediction given a single response variable (usually a phenotype)
#' using the \code{\link[BGLR]{BGLR}} function. Unlike other snpR functions,
#' this returns the resulting model directly, so overwrite with caution.
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
#' interest should be selected with \code{\link{subset_snpR_data}}. See
#' examples.
#'
#' See documentation for \code{\link[BGLR]{BGLR}} for more details and for a
#' full list of references.
#'
#' @param x snpRdata object
#' @param facets character, default NULL. Categorical metadata variables by
#'   which to break up analysis. See \code{\link{Facets_in_snpR}} for more
#'   details.
#' @param response character. Name of the column containing the response
#'   variable of interest. Must match a column name in sample metadata.
#' @param iterations numeric. Number of iterations to run the MCMC chain for.
#' @param burn_in numeric. Number of burn in iterations to run prior to the MCMC
#'   chain.
#' @param thin numeric. Number of iterations to discard between each recorded
#'   data point.
#' @param model character, default "BayesB". Prediction model to use, see
#'   description for the ETA argument in \code{\link[BGLR]{BGLR}}.
#' @param interpolate character, default "bernoulli". Interpolation method for
#'   missing data. Options: 
#'   * bernoulli: binomial draws for the minor allele.
#'   * af: insertion of the average allele frequency
#'   * iPCA: As a slower but more accurate alternative to "af"
#'   interpolation, "iPCA" may be selected. This an iterative PCA approach to
#'   interpolate based on SNP/SNP covariance via
#'   \code{\link[missMDA]{imputePCA}}. If the ncp argument is not defined, the
#'   number of components used for interpolation will be estimated using
#'   \code{\link[missMDA]{estim_ncpPCA}}. In this case, this method is much
#'   slower than the other methods, especially for large datasets. Setting an
#'   ncp of 2-5 generally results in reasonable interpolations without the time
#'   constraint.
#' @param ncp numeric or NULL, default NULL. Used only if \code{iPCA}
#'   interpolation is selected. Number of components to consider for iPCA sn
#'   format interpolations of missing data. If null, the optimum number will be
#'   estimated, with the maximum specified by ncp.max. This can be very slow.
#' @param ncp.max numeric, default 5. Used only if \code{iPCA}
#'   interpolation is selected. Maximum number of components to check for
#'   when determining the optimum number of components to use when interpolating
#'   sn data using the iPCA approach.
#' @param par numeric or FALSE, default FALSE. If a number specifies the number
#'   of processing cores to use \emph{across facet levels}. Not used if only one
#'   facet level.
#' @param verbose Logical, default FALSE. If TRUE, some progress updates will be
#'   printed to the console.
#' @param ... additional arguments passed to \code{\link[BGLR]{BGLR}}
#'
#' @export
#' @author William Hemstrom
#' @references Pérez, P., and de los Campos, G. (2014). \emph{Genetics}.
#'
#' @return A list containing: two parts:
#'   * x: The provided snpRdata object with effect sizes merged in.
#'   * models: Other model results, a list containing:
#'     + model: The model output from  BGLR. See \code{\link[BGLR]{BGLR}}.
#'     + h2: Estimated heritability of the response variable.
#'     + predictions: A data.frame containing the provided phenotypes and the 
#'       predicted Breeding Values (BVs) for those phenotypes.
#'
#' @examples
#' # run and plot a basic prediction
#' ## add some dummy phenotypic data.
#' dat <- stickSNPs
#' sample.meta(dat) <- cbind(weight = rnorm(ncol(stickSNPs)), 
#'                           sample.meta(stickSNPs))
#' ## run prediction
#' gp <- run_genomic_prediction(dat, response = "weight", iterations = 1000, 
#'                              burn_in = 100, thin = 10)
#' ## dummy phenotypes vs. predicted Breeding Values for dummy predictions.
#' # given that weight was randomly assigned, definitely overfit!
#' with(gp$models$.base_.base$predictions, plot(phenotype, predicted_BV)) 
#' ## fetch estimated loci effects
#' get.snpR.stats(gp$x, stats = "genomic_prediction")
#' 
#' \dontrun{
#' # with facets, not run
#' gp <- run_genomic_prediction(gp$x, facets = "pop", response = "weight", 
#'                              iterations = 1000, burn_in = 100, thin = 10)
#' get.snpR.stats(gp$x, facets = "pop", stats = "genomic_prediction")
#' }
run_genomic_prediction <- function(x, facets = NULL, response, iterations,
                                   burn_in, thin,
                                   model = "BayesB", interpolate = "bernoulli",
                                   ncp = NULL, ncp.max = 5,
                                   par = FALSE, verbose = FALSE, ...){
  .snp.id <- facet <- subfacet <- NULL
  
  #===============sanity checks============
  if(!is.snpRdata(x)){
    stop("x must be a snpRdata object.\n")
  }
  if(!.is.bi_allelic(x)){
    stop("This function is not yet supported for non-bi-allelic markers.\n")
  }
  msg <- character(0)
  
  if(length(.check.snpR.facet.request(x, facets)) == 1 & !isFALSE(par)){
    par <- FALSE
  }
  
  if(isFALSE(interpolate)){
    msg <- c(msg, "Interpolation must be performed for genomic prediction.\n")
  }
  
  if(!response %in% colnames(sample.meta(x))){
    msg <- c(msg, paste0("No column matching response '", response, "' found in sample metadata.\n"))
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  .check.installed("BGLR"); .check.installed("stringi")
  
  #===============function====================
  run_BGLR <- function(sub.x, ...){

    # generate a random file handle for multiple jobs
    handle <- .rand_strings(1, length = 10)
    while(length(grep("handle", list.files())) > 0){
      handle <- .rand_strings(1, length = 10)
    }
    
    
    if(length(sub.x@sn$sn) != 0){
      if(sub.x@sn$type != interpolate){
        suppressWarnings(sub.x@sn$sn <- format_snps(sub.x, "sn", interpolate = interpolate, ncp = ncp, ncp.max = ncp.max))
        sub.x@sn$type <- interpolate
      }
    }
    else{
      suppressWarnings(sub.x@sn$sn <- format_snps(sub.x, "sn", interpolate = interpolate, ncp = ncp, ncp.max = ncp.max))
      sub.x@sn$type <- interpolate
    }
    
    sn <- sub.x@sn$sn
    sn <- sn[,-c(1:(ncol(sub.x@snp.meta) - 1))]
    
    # prepare the BGLR input
    sn <- t(sn)
    phenotypes <- sub.x@sample.meta[,response]
    if(!is.numeric(phenotypes)){
      phenotypes <- as.numeric(as.factor(phenotypes)) - 1
      if(length(unique(phenotypes)) != 2){
        stop("Only two unique phenotypes allowed if not quantitative.\n")
      }
    }
    colnames(sn) <- paste0("m", 1:ncol(sn)) # marker names
    rownames(sn) <- paste0("s", 1:nrow(sn)) # ind IDS
    ETA <- list(list(X = sn, model = "BayesB", saveEffects = T)) # need to adjust this when I get around to allowing for more complicated models
    
    BGLR_mod <- BGLR::BGLR(y = phenotypes, ETA = ETA, nIter = iterations + burn_in, 
                           burnIn = burn_in, thin = thin, saveAt = handle, verbose = verbose, ...)
    
    # grab h2 estimate
    B <- BGLR::readBinMat(paste0(handle,'ETA_1_b.bin'))
    h2 <- rep(NA,nrow(B))
    varU <- h2
    varE <- h2
    for(i in 1:length(h2)){
      u <- sn%*%B[i,]
      varU[i] <- stats::var(u)
      varE[i] <- stats::var(phenotypes-u)
      h2[i] <- varU[i]/(varU[i] + varE[i])
    }
    h2 <- mean(h2)
    
    pdat <- data.frame(phenotype = phenotypes, predicted_BV = sn%*%BGLR_mod$ETA[[1]]$b)
    
    # clean
    file.remove(list.files(pattern = handle))
    
    effects <- cbind.data.frame(snp.meta(sub.x), BGLR_mod$ETA[[1]]$b, row.names = NULL)
    
    
    # return(pdat)
    return(list(model = BGLR_mod, h2 = h2, predictions = pdat, effects = effects))
  }
  
  
  #===============apply===========
  # prep and run
  facets <- .check.snpR.facet.request(x, facets)
  x <- .add.facets.snpR.data(x, facets)

  
  # apply
  out <- .apply.snpR.facets(x, facets, req = "snpRdata", case = "ps", fun = run_BGLR, response = response,
                           interpolate = interpolate,  par = par, verbose = verbose, ...)
  
  
  
  #===========return================
  # process output
  stats <- vector("list", length = length(out))
  models <- stats
  ecol_name <- paste0(response, "_gp_effect")
  for(i in 1:length(out)){
    type <- ifelse(out[[i]]$.fm[1] == ".base", ".base", "sample")
    stats[[i]] <- cbind(out[[i]]$.fm, facet.type = type, out[[i]]$effects, row.names = NULL)
    colnames(stats[[i]])[ncol(stats[[i]])] <- ecol_name
    
    # correct for missing snps
    missing <- which(!snp.meta(x)$.snp.id %in% stats[[i]]$.snp.id)
    if(length(missing) > 0){
      missing <-  cbind.data.frame(out[[i]]$.fm, facet.type = type, snp.meta(x)[missing,], effect = NA, row.names = NULL)
      colnames(missing)[ncol(missing)] <- ecol_name
      stats[[i]] <- rbind(stats[[i]], missing)
    }
    
    
    models[[i]] <- out[[i]][-which(names(out[[i]]) %in% c("effects", ".fm"))]
    names(models)[i] <- paste0(out[[i]]$.fm[1,], collapse = "_")
  }
  
  # merge and return
  stats <- data.table::rbindlist(stats)
  dplyr::arrange(stats, .snp.id, facet, subfacet)
  
  x <- .merge.snpR.stats(x, stats)
  x <- .update_citations(x, "Perez2014", "genomic_prediction", paste0("Genomic prediction against ", response, " via BGLR with model ", model, ". You may also want to cite the model you used (BayesB, etc.). Try ?BGLR::BGLR for details."))
  

  return(list(x = x, models = models))
}


#' Run a single cross-validation with BGLR.
#'
#' Run genomic prediction given a single response variable (usually a phenotype)
#' using the \code{\link[BGLR]{BGLR}} function. Unlike other snpR functions,
#' this returns the resulting model directly, so overwrite with caution. This
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
#' interest should be selected with \code{\link{subset_snpR_data}}. See
#' examples.
#'
#' The portion of samples left out for cross-validation is defined by the
#' cross_percentage argument. Specifically, the proportion given will be used to
#' create the model, the remaining portion will be left out. For example, if
#' cross_percentage = 0.9, 90% of the samples will be used to run the model and
#' the remaining 10% will be used for cross-validation. Alternatively, if
#' cross_samples are provided, the specified samples (by column index) will be
#' used used for cross validation instead. This may be useful for a systematic
#' leave-one-out cross-validation.
#'
#'
#' See documentation for \code{\link[BGLR]{BGLR}} for more details and for a
#' full list of references.
#'
#' @param x snpRdata object
#' @param response character. Name of the column containing the response
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
#' @param interpolate character, default "bernoulli". Interpolation method for
#'   missing data. Options: 
#'   * bernoulli: binomial draws for the minor allele
#'   * af: insertion of the average allele frequency
#'
#' @export
#' @author William Hemstrom
#' @references Pérez, P., and de los Campos, G. (2014). \emph{Genetics}.
#' @return A list containing:
#'   * model: The results from \code{\link{run_genomic_prediction}}
#'   * model.samples: Indices of the samples used to construct the model. 
#'   * cross.samples: Indices of the samples used to cross-validate the model.
#'   * comparison: A data.frame containing the observed and predicted 
#'     phenotypes/Breeding Values for the cross-validation samples.
#'   * rsq: The r^2 value for the observed and predicted phenotypes/Breeding 
#'     Values for the cross-validation samples.
#'
#' @examples
#' # run and plot a basic prediction
#' ## add some dummy phenotypic data.
#' dat <- stickSNPs
#' sample.meta(dat) <- cbind(weight = rnorm(ncol(stickSNPs)), 
#'                           sample.meta(stickSNPs))
#' ## run cross_validation
#' cross_validate_genomic_prediction(dat, response = "weight", 
#'                                   iterations = 1000, burn_in = 100, 
#'                                   thin = 10)
#' 
cross_validate_genomic_prediction <- function(x, response, iterations = 10000,
                                              burn_in = 1000, thin = 100, cross_percentage = 0.9,
                                              model = "BayesB", cross_samples = NULL, plot = TRUE, 
                                              interpolate = "bernoulli"){

  #===============sanity checks============
  if(!is.snpRdata(x)){
    stop("x must be a snpRdata object.\n")
  }
  if(!.is.bi_allelic(x)){
    stop("This function is not yet supported for non-bi-allelic markers.\n")
  }
  
  .check.installed("BGLR")
  #===============run======================

  # check that the response is numeric
  phenotypes <- x@sample.meta[,response]
  if(!is.numeric(phenotypes)){
    phenotypes <- as.numeric(as.factor(phenotypes)) - 1
    if(length(unique(phenotypes)) != 2){
      stop("Only two unique phenotypes allowed if not quantitative.\n")
    }
  }
  x@sample.meta[,response] <- phenotypes


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
  utils::capture.output(suppressWarnings(suppressMessages(sub.x <- subset_snpR_data(x, .samps = model.samples))))
  model <- run_genomic_prediction(sub.x, response = response, iterations =  iterations,
                                  burn_in = burn_in, thin = thin, model = model)

  # check accuracy
  cross.samples <- (1:ncol(x))[-sort(model.samples)]
  if(length(x@sn$sn) != 0){
    if(x@sn$type != interpolate){
      suppressWarnings(x@sn$sn <- format_snps(x, "sn", interpolate = interpolate))
      x@sn$type <- interpolate
    }
  }
  else{
    suppressWarnings(x@sn$sn <- format_snps(x, "sn", interpolate = interpolate))
    x@sn$type <- interpolate
  }
  whole.sn <- x@sn$sn
  whole.sn <- whole.sn[,-(1:(ncol(x@snp.meta) - 1))]
  whole.sn <- .interpolate_sn(whole.sn)
  new.sn <- whole.sn[,cross.samples]
  new.sn <- t(new.sn)
  pred.phenos <- new.sn%*%model$models$.base_.base$model$ETA[[1]]$b
  pdat <- as.data.frame(cbind(observed = x@sample.meta[cross.samples, response], predicted = pred.phenos), stringsAsFactors = F)
  colnames(pdat) <- c("observed", "predicted")

  # plot
  if(plot == TRUE){
    observed <- predicted <- NULL
    tplot <- ggplot2::ggplot(as.data.frame(pdat), ggplot2::aes(observed, predicted)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "lm") +
      ggplot2::theme_bw()
    
    return(list(model = model, model.samples = model.samples, cross.samples = cross.samples, comparison = pdat,
                rsq = summary(stats::lm(observed~predicted, pdat))$r.squared,
                plot = tplot))
  }

  return(list(model = model, model.samples = model.samples, cross.samples = cross.samples, comparison = pdat,
              rsq = summary(stats::lm(observed~predicted, pdat))$r.squared))

}

#' Run case/control or quantitative association tests on SNP data.
#'
#' Runs several different association tests on SNP data. The response variable
#' must have only two different categories (as in case/control) for most test
#' types, although the "gmmat.score" method supports quantitative traits. Tests
#' may be broken up by sample-specific facets.
#'
#' Several methods can be used: Armitage, chi-squared, and odds ratio. For The
#' Armitage approach weights should be provided to the "w" argument, which
#' specifies the weight for each possible genotype (homozygote 1, heterozygote,
#' homozygote 2). The default, c(0,1,2), specifies an additive model. The
#' "gmmat.score" method uses the approach described in Chen et al. (2016) and
#' implemented in the \code{\link[GMMAT]{glmmkin}} and
#' \code{\link[GMMAT]{glmm.score}} functions. For this method, a 'G' genetic
#' relatedness matrix is first created using the
#' \code{\link[AGHmatrix]{Gmatrix}} function according to Yang et al 2010.
#'
#' Multi-category data is currently not supported, but is under development.
#'
#' Facets are specified as described in \code{\link{Facets_in_snpR}}. NULL and
#' "all" facet specifications function as described.
#'
#' Note that if all individuals in one facet level have one of the two possible
#' phenotypes, each association test will return NA or NaN. As a result, if the
#' response variable overlaps perfectly with the facet variable, all results
#' will be NA or NaN.
#'
#' If the chisq method is used, a column will also be returned that specifies
#' which allele (major or minor) is associated with which phenotype.
#'
#'
#'
#' @param x snpRdata object
#' @param facets character, default NULL. Categorical metadata variables by
#'   which to break up analysis. See \code{\link{Facets_in_snpR}} for more
#'   details.
#' @param response character. Name of the column containing the response
#'   variable of interest. Must match a column name in sample metadata. Response
#'   must be categorical, with only two categories.
#' @param method character, default "gmmat.score". Specifies association method.
#'   Options: 
#'   * gmmat.score: Population/family structure corrected mlm approach, based on
#'     Chen et al (2016).
#'   * armitage: Armitage association test, based on Armitage (1955). 
#'   * odds_ratio: Log odds ratio test.
#'   * chisq: Chi-squared test.
#'   See description for more details.
#' @param w numeric, default c(0, 1, 2). Weight variable for each genotype for
#'   the Armitage association method. See description for details.
#' @param formula character, default set to response ~ 1. Null formula for the
#'   response variable, as described in \code{\link[stats]{formula}}.
#' @param family.override character, default NULL. Provides an alternative model
#'   family object to use for GMMAT GWAS regression. By default, uses
#'   \code{\link[stats]{gaussian}}, link = "identity" for a quantitative
#'   phenotype and \code{\link[stats]{binomial}}, link = "logit" for a
#'   categorical phenotype.
#' @param maxiter numeric, default 500. Maximum iterations to use when fitting
#'   the glmm when using the gmmat.score option.
#' @param sampleID character, default NULL. Optional, the name of a column in
#'   the sample metadata to use as a sampleID when using the gmmat.score option.
#' @param Gmaf numeric, default NULL. If using the "GMMAT" option, can provide
#'   and optional minor allele frequency filter used when constructing the
#'   relatedness matrix.
#' @param par numeric or FALSE, default FALSE. Number of parallel cores to use
#'   for computation.
#' @param cleanup logical, default TRUE. If TRUE, files produced by and for the
#'   "gmmat" option will be removed after completion.
#' @param verbose Logical, default FALSE. If TRUE, some progress updates will be
#'   printed to the console.
#'
#' @author William Hemstrom
#' @author Keming Su
#' @author Avani Chitre
#' @export
#'
#' @return A snpRdata object with the resulting association test results merged
#'   into the stats socket.
#'
#' @references Armitage (1955). Tests for Linear Trends in Proportions and
#'   Frequencies. \emph{Biometrics}.
#' @references Chen et al. (2016). Control for Population Structure and
#'   Relatedness for Binary Traits in Genetic Association Studies via Logistic
#'   Mixed Models. \emph{American Journal of Human Genetics}.
#' @references Yang et al. (2010). Common SNPs explain a large proportion of the
#'   heritability for human height. \emph{Nature Genetics}.
#'
#' @examples
#'   # add a dummy phenotype and run an association test.
#'   x <- stickSNPs
#'   sample.meta(x)$phenotype <- sample(c("A", "B"), nsamps(stickSNPs), TRUE)
#'   x <- calc_association(x, facets = c("pop"), response = "phenotype", method = "armitage")
#'   get.snpR.stats(x, "pop", "association")
#' 
calc_association <- function(x, facets = NULL, response, method = "gmmat.score", w = c(0,1,2),
                             formula = NULL, family.override = FALSE, maxiter = 500, 
                             sampleID = NULL, Gmaf = 0, par = FALSE, 
                             cleanup = TRUE, verbose = FALSE){
  #==============sanity checks===========
  if(!is.snpRdata(x)){
    stop("x must be a snpRdata object.\n")
  }
  
  if(!.is.bi_allelic(x)){
    stop("This function is not yet supported for non-bi-allelic markers.\n")
  }
  
  # response
  msg <- character()
  if(length(response) != 1){
    msg <- c(msg,
             paste0("Only one response variable permitted."))
  }
  if(grepl("(?<!^)\\.", response, perl = TRUE)[1]){
    msg <- c(msg,
             paste0("Only one sample-specific category allowed (e.g. pop but not fam.pop)."))
  }
  if(!response[1] %in% colnames(x@sample.meta)){
    msg <- c(msg,
             paste0("Response variable must be present in sample metadata."))
  }
  else{
    if(length(unique(x@sample.meta[,response])) != 2){
      if(method %in% c("armitage", "odds_ratio", "chisq")){
        msg <- c(msg,
                 paste0("Only two categories allowed for response variable for method: ", method, "."))
      }
    }
    if(any(is.na(x@sample.meta[,response]))){
      msg <- c(msg, "Some NAs found in response. This can cause issues!")
    }
  }

  # method
  good.methods <- c("armitage", "odds_ratio", "chisq", "gmmat.score")
  if(!method %in% good.methods){
    msg <- c(msg,
             paste0("Method not accepted. Accepted methods: ", paste0(good.methods, collapse = ", "), "."))
  }

  # w
  if(method == "armitage"){
    if(!is.numeric(w)){
      msg <- c(msg,
               "w must be numeric.")
    }
    if(length(w) != 3){
      msg <- c(msg,
               "w must be length 3.")
    }
  }

  # vars for gmmat
  if(method == "gmmat.score"){
    if(!is.null(sampleID)){
      if(!sampleID %in% x@sample.meta){
        msg <- c(msg,
                 paste0("Sample ID column: ", sampleID, " not found in sample metadata."))
      }
    }

    if(is.null(formula)){
      formula <- paste0(response, " ~ 1")
    }
    else{
      res <- try(formula(formula), silent = T)
      if(methods::is(res, "try-error")){
        msg <- c(msg,
                 "formula must be a valid formula. Type ?formula for help.\n")
      }
      else{
        phen <- strsplit(as.character(formula), "~")
        phen <- phen[[1]][1]
        phen <- gsub(" ", "", phen)
        if(phen != response){
          msg <- c(msg,
                   "The response variable in the provided formula must be the same as that provided to the response argument.\n")
        }

        # see which covariates we need
        cvars <- stats::terms(formula(formula))
        cvars <- attr(cvars, "term.labels")
        bad.cvars <- which(!(cvars %in% colnames(x@sample.meta)))
        if(length(bad.cvars) > 0){
          msg <- c(msg,
                   "Some covariates specified in formula not found in sample metadata: ", paste0(cvars[bad.cvars], collapse = ", "), ".\n")
        }


      }
    }

    if(!"GMMAT" %in% utils::installed.packages()){
      .check.installed("SeqVarTools", "bioconductor")
      .check.installed("GMMAT")
      
      # msg <- c(msg,
      #          paste0("The gmmat.score method requires the GMMAT package. This can be installed via
      #                 install.packages('GMMAT'), but requires the SeqVar and SeqVarTools bioconductor packages.
      #                 These can be installed via BiocManager::install(c('SeqVar', 'SeqVarTools')).
      #                 If BiocManager is not installed, it can be installed via install.packages('BiocManager')."))
    }
    .check.installed("AGHmatrix")

  }
  
  if(length(msg) > 0){
    stop(msg)
  }

  #==============functions=============

  calc_armitage <- function(a, w, ...){#where a is the matrix you want to run the test on, and w is the weights
    a <- as.matrix(a)

    # figure out "control" and "case" names
    poss.id <- unique(gsub("^[ATCG]{2}_", "", colnames(a), perl = T))
    control.id <- poss.id[1]
    case.id <- poss.id[2]

    #separate controls
    control.cols <- grep(paste0("_", control.id, "$"), colnames(a))
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
    
    hom_totals <- hom_controls + hom_cases
    
    
    # loop through each genotype to figure out p and q--should only be four loops max!
    pp_cont <- pp_case <- 0
    maxs <- matrixStats::rowMaxs(hom_totals)
    for(i in 1:ncol(hom_totals)){
      is_p <- which(hom_totals[,i] == maxs)
      pp_cont[is_p] <- hom_controls[is_p,i]
      pp_case[is_p] <- hom_cases[is_p,i]
    }
    
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
    return(stats::pchisq(chi, 1, lower.tail = F))
  }
  odds.ratio.chisq <- function(cast_ac, method, ...){
   
    # odds ratio
    a <- cast_ac[[1]]
    b <- cast_ac[[2]]
    c <- cast_ac[[3]]
    d <- cast_ac[[4]]
    #============odds ratio==========
    if(method == "odds_ratio"){
      s.e <- sqrt(1/a + 1/b + 1/c + 1/d)
      odds <- log10((a/b)/(c/d))
      return(data.frame(log_odds_ratio = odds, se = log10(s.e)))
    }

    #============chisq===============
    else if(method == "chisq"){
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

      # figure out which allele is more associated with the case
      oms.a <- a - ea
      oms.a <- sign(oms.a) # if positive, more ni1 with case than expected, so ni1 is associated with case.
      oms.a[oms.a == 1] <- "major"
      oms.a[oms.a == -1] <- "minor"
      oms.a[oms.a == 0] <- NA
      asso.allele <- paste0(oms.a, "_", gsub("^ni1_", "", colnames(cast_ac)[1]))


      chi.stat <- chi.a + chi.b + chi.c + chi.d
      chi.p <- stats::pchisq(chi.stat, 1, lower.tail = F)
      return(data.frame(chi_stat = chi.stat, chi_p = chi.p, associated_allele = asso.allele))
    }
  }
  run_gmmat <- function(sub.x, form, iter, sampleID, family.override, Gmaf, ...){
    # sn format
    sn <- format_snps(sub.x, "sn", interpolate = FALSE)
    sn <- sn[,-c(which(colnames(sn) %in% colnames(sub.x@snp.meta)))]

    ## G matrix
    .make_it_quiet(G <- AGHmatrix::Gmatrix(t(sn), missingValue = NA, method = "Yang", maf = Gmaf))
    if(is.null(sampleID)){
      sampleID <- ".sample.id"
    }
    colnames(G) <- sub.x@sample.meta[,sampleID]
    rownames(G) <- sub.x@sample.meta[,sampleID]

    ## input data
    sub.x <- calc_maf(sub.x)
    stats <- sub.x@stats[sub.x@stats$facet == ".base",]
    asso.in <- cbind(stats$major, stats$minor, sub.x@snp.meta, sn)
    colnames(asso.in)[1:2] <- c("major", "minor")
    colnames(asso.in)[-c(1:(2 + ncol(sub.x@snp.meta)))] <- sub.x@sample.meta[,sampleID]
    phenos <- sub.x@sample.meta

    # set parameters
    ## response
    if(is.character(phenos[,response])){
      phenos[,response] <- as.numeric(as.factor(phenos[,response])) - 1 # set the phenotype to 0,1,etc
    }
    ## family
    if(family.override != FALSE){
      family <- family.override
    }
    else{
      if(length(unique(phenos[,response])) > 2){
        family <- stats::gaussian(link = "identity")
      }
      else if(length(unique(phenos[,response])) == 2){
        family <- stats::binomial(link = "logit")
      }
      else{
        stop("Less than two unique phenotypes.\n")
      }
    }


    # run null model
    .make_it_quiet(mod <- GMMAT::glmmkin(fixed = formula,
                                                          data = phenos,
                                                          kins = G,
                                                          id = sampleID,
                                                          family = family,
                                                          maxiter = iter, verbose = verbose))

    # run the test
    nmeta.col <- 2 + ncol(sub.x@snp.meta)
    utils::write.table(asso.in, "asso_in.txt", sep = "\t", quote = F, col.names = F, row.names = F)
    
    .suppress_specific_warning(score.out <- GMMAT::glmm.score(mod, "asso_in.txt",
                                                              "asso_out_score.txt",
                                                              infile.ncol.skip = nmeta.col,
                                                              infile.ncol.print = 1:nmeta.col,
                                                              infile.header.print = colnames(asso.in)[1:nmeta.col], 
                                                              verbose = verbose), 
                               warn_to_suppress = "Assuming the order of individuals in infile")
    
    score.out <- utils::read.table("asso_out_score.txt", header = T, stringsAsFactors = F)

    if(cleanup){
      file.remove(c("asso_in.txt", "asso_out_score.txt"))
    }
    return(score.out)
  }

  #==============run the function=========
  # check facets
  facets <- .check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    .make_it_quiet(x <- .add.facets.snpR.data(x, facets))
  }

  if(method == "armitage"){
    out <- .apply.snpR.facets(x, facets = facets, req = "cast.gs", case = "ps", fun = calc_armitage, response = response, w = w, verbose = verbose)
    x <- .update_citations(x, "Armitage1955", "association", paste0("Association test against ", response, "."))
  }
  else if(method == "odds_ratio" | method == "chisq"){
    out <- .apply.snpR.facets(x, facets = facets, req = "cast.ac", case = "ps", fun = odds.ratio.chisq, response = response, method = method, verbose = verbose)
  }
  else if(method == "gmmat.score"){
    out <- .apply.snpR.facets(x, facets = facets, req = "snpRdata", case = "ps", Gmaf = Gmaf, fun = run_gmmat, response = response, form = formula, iter = maxiter, sampleID = sampleID, family.override = family.override, verbose = verbose)
    x <- .update_citations(x, "Chen2016", "association", paste0("Association test against ", response, "."))
    x <- .update_citations(x, "Yang2010", "association", paste0("G-matrix creation for association test against ", response, "."))
  }
  
  x <- .update_calced_stats(x, facets, paste0("association_", method))

  x <- .merge.snpR.stats(x, out)
  return(x)
}


#' Run a RANGER random forest using snpRdata for a given phenotype or model.
#'
#' Creates forest machine learning models using snpRdata objects via interface
#' with \code{\link[ranger]{ranger}}. Models can be created either for a
#' specific phenotype with no-covariates or using a formula which follows the
#' basic format specified in \code{\link{formula}}.
#'
#' Random forest models can be created across multiple facets of the data at
#' once following the typical snpR framework explained in
#' \code{\link{Facets_in_snpR}}. Since RF models are calculated without
#' allowing for any SNP-specific categories (e.g. independent of chromosome
#' etc.), any sample level facets provided will be ignored. As usual, if facets
#' is set to NULL, an RF will be calculated for all samples without splitting
#' across any sample metadata categories.
#'
#' Since the \code{\link[ranger]{ranger}} RF implementation can behave
#' unexpectedly when given incomplete data, missing genotypes will be imputed.
#' Imputation can occur either via the insertion of the average allele frequency
#' or via binomial draws for the minor allele using the "af" or "bernoulli"
#' options for the "interpolate" argument.
#'
#' Extra arguments can be passed to \code{\link[ranger]{ranger}}.
#'
#' In general, random forest parameters should be tuned so as to reduce the
#' out-of-bag error rates (OOB-ER). This value is visible in the returned object
#' under the model lists. Simply calling a specific model will output the
#' OOB-ER, and they are also stored under the 'prediction.error' name in the
#' model. For details on tuning RF models, we recommend Goldstein et al. (2011).
#'
#' For more detail on the random forest model and ranger arguments, see
#' \code{\link[ranger]{ranger}}.
#'
#' @param x snpRdata object.
#' @param facets character, default NULL. Categorical metadata variables by
#'   which to break up analysis. See \code{\link{Facets_in_snpR}} for more
#'   details.
#' @param response character. Name of the column containing the response
#'   variable of interest. Must match a column name in sample metadata.
#' @param formula character, default NULL. Model for the response variable, as
#'   described in \code{\link[stats]{formula}}. If NULL, the model will be
#'   equivalent to response ~ 1.
#' @param num.trees numeric, default 10000. Number of trees to grow. Higher
#'   numbers will increase model accuracy, but increase calculation time. See
#'   \code{\link[ranger]{ranger}} for details.
#' @param mtry numeric, default is the square root of the number of SNPs. Number
#'   of variables (SNPs) by which to split each node. See
#'   \code{\link[ranger]{ranger}} for details.
#' @param importance character, default "impurity_corrected". The method by
#'   which SNP importance is determined. Options: 
#'   * impurity
#'   * impurity_corrected
#'   * permutation. 
#'   See \code{\link[ranger]{ranger}} for details.
#' @param interpolate character, default "bernoulli". Interpolation method for
#'   missing data. Options: 
#'   * bernoulli: binomial draws for the minor allele.
#'   * af: insertion of the average allele frequency
#'   * iPCA: As a slower but more accurate alternative to "af"
#'   interpolation, "iPCA" may be selected. This an iterative PCA approach to
#'   interpolate based on SNP/SNP covariance via
#'   \code{\link[missMDA]{imputePCA}}. If the ncp argument is not defined, the
#'   number of components used for interpolation will be estimated using
#'   \code{\link[missMDA]{estim_ncpPCA}}. In this case, this method is much
#'   slower than the other methods, especially for large datasets. Setting an
#'   ncp of 2-5 generally results in reasonable interpolations without the time
#'   constraint.
#' @param ncp numeric or NULL, default NULL. Used only if \code{iPCA}
#'   interpolation is selected. Number of components to consider for iPCA sn
#'   format interpolations of missing data. If null, the optimum number will be
#'   estimated, with the maximum specified by ncp.max. This can be very slow.
#' @param ncp.max numeric, default 5. Used only if \code{iPCA}
#'   interpolation is selected. Maximum number of components to check for
#'   when determining the optimum number of components to use when interpolating
#'   sn data using the iPCA approach.
#' @param pvals logical, default TRUE. Determines if p-values should be
#'   calculated for importance values. If the response variable is quantitative,
#'   no p-values will be returned, since they must be calculated via permutation
#'   and this is very slow. For details, see
#'   \code{\link[ranger]{importance_pvalues}}.
#' @param par numeric, default FALSE. Number of parallel computing cores to use
#'   for computing RFs across multiple facet levels or within a single facet if
#'   only a single category is run (either a one-category facet or no facet).
#' @param ... Additional arguments passed to \code{\link[ranger]{ranger}}.
#'
#' @return A list containing: 
#'   * data: A snpRdata object with RF importance values merged in to the stats 
#'     slot. 
#'   * models: A named list containing both the models and data.frames 
#'     containing the predictions vs observed phenotypes.
#'
#' @references Wright, Marvin N and Ziegler, Andreas. (2017). ranger: A Fast
#'   Implementation of Random Forests for High Dimensional Data in C++ and R.
#'   \emph{Journal of Statistical Software}.
#' @references Goldstein et al. (2011). Random forests for genetic association
#'   studies. \emph{Statistical Applications in Genetics and Molecular Biology}.
#'
#' @seealso \code{\link[ranger]{ranger}} \code{\link[ranger]{predict.ranger}}.
#'
#' @author William Hemstrom
#' @export
#'
#' @examples
#' \dontrun{
#' # run and plot a basic rf
#' ## add some dummy phenotypic data.
#' dat <- stickSNPs
#' sample.meta(dat) <- cbind(weight = rnorm(ncol(stickSNPs)), sample.meta(stickSNPs))
#' ## run rf
#' rf <- run_random_forest(dat, response = "weight", pvals = FALSE)
#' rf$models
#' ## dummy phenotypes vs. predicted
#' with(rf$models$.base_.base$predictions, plot(pheno, predicted)) # not overfit
#' }
#'
run_random_forest <- function(x, facets = NULL, response, formula = NULL,
                              num.trees = 10000, mtry = NULL,
                              importance = "impurity_corrected",
                              interpolate = "bernoulli", ncp = NULL, 
                              ncp.max = 5, pvals = TRUE, par = FALSE, ...){
  .formula <- formula
  rm(formula)
  #=========sanity checks=======================
  if(!is.snpRdata(x)){
    stop("x must be a snpRdata object.\n")
  }
  if(!.is.bi_allelic(x)){
    stop("This function is not yet supported for non-bi-allelic markers.\n")
  }
  
  x <- filter_snps(x, non_poly = TRUE, verbose = FALSE)
  
  .check.installed("ranger")
  #==========run============
  
  
  run_ranger <- function(sub.x, .formula = NULL, par = FALSE,...){
    #================grab data=====================
    ## sn format
    if(length(sub.x@sn$sn) != 0){
      if(sub.x@sn$type != interpolate){
        sn <- format_snps(sub.x, "sn", interpolate = interpolate, ncp = ncp, ncp.max = ncp.max)
      }
      else{
        sn <- sub.x@sn$sn
      }
    }
    else{
      sn <- format_snps(sub.x, "sn", interpolate = interpolate, ncp = ncp, ncp.max = ncp.max)
    }
    sn <- sn[,-c(which(colnames(sn) %in% colnames(sub.x@snp.meta)))]
    sn <- t(sn)


    ## attach phenotype, anything else in the formula if given.
    if(!is.null(.formula)){
      msg <- character()
      res <- try(formula(.formula), silent = T)
      if(methods::is(res, "try-error")){
        msg <- c(msg,
                 "formula must be a valid formula. Type ?formula for help.\n")
      }
      else{
        if(all.vars(.formula)[1] != response){
          msg <- c(msg,
                   "The response variable in the provided formula must be the same as that provided to the response argument.\n")
        }

        # see which covariates we need
        cvars <- stats::terms(formula(.formula))
        cvars <- attr(cvars, "term.labels")
        bad.cvars <- which(!(cvars %in% colnames(sub.x@sample.meta)))
        if(length(bad.cvars) > 0){
          msg <- c(msg,
                   "Some covariates specified in formula not found in sample metadata: ", paste0(cvars[bad.cvars], collapse = ", "), ".\n")
        }
        
        if(length(msg) > 0){
          stop(msg)
        }

        colnames(sn) <- paste0("loc_", 1:ncol(sn))
        ocn <- colnames(sn)
        sn <- as.data.frame(sn)
        sn <- cbind(sub.x@sample.meta[,response], sub.x@sample.meta[,cvars], sn)
        colnames(sn)[1:(length(cvars) + 1)] <- c(response, cvars)
        
        myterms <- all.vars(res)
        .formula <- stats::as.formula(paste0(myterms[1], " ~ ", paste0(cvars, collapse = " + "), " + ", paste0(ocn, collapse = " + ")))
        
        # cat("Model:\n")
        # print(paste0(as.character(res), " + ", paste0(ocn, collapse = " + ")))
        # cat("\nFitting...\n")
        # if(length(paste0(as.character(res), " + ", paste0(ocn, collapse = " + "))) > 1){
        #   browser()
        # }
        # 
        # # reset formula:
        # .formula <- formula(paste0(as.character(res), " + ", paste0(ocn, collapse = " + ")))
      }
    }
    else{
      cvars <- character()
      sn <- cbind.data.frame(sub.x@sample.meta[,response], sn, stringsAsFactors = F)
    }
    colnames(sn)[1] <- response
    # convert to a factor if 2 categories or character.
    if(length(unique(as.character(sn[,1]))) == 2 | is.character(sn[,1])){
      sn[,1] <- as.factor(sn[,1])
    }


    #================run the model:=====================
    # set par correctly, NULL if par is FALSE or we are looping through facet levels. Otherwise (.base and par provided), tpar = par.
    if(par == FALSE){
      par <- NULL
    }

    # run with or without formula
    if(is.null(.formula)){
      rout <- ranger::ranger(dependent.variable.name = response,
                             data = sn,
                             num.trees = num.trees,
                             mtry = mtry,
                             importance = importance, num.threads = par, verbose = T, ...)
    }
    else{
      rout <- ranger::ranger(formula = .formula,
                             data = sn,
                             num.trees = num.trees,
                             mtry = mtry,
                             importance = importance, num.threads = par, verbose = T, ...)
    }

    #===============make sense of output=================
    # get p-values:
    if(pvals){
      if(length(unique(as.character(x@sample.meta[,response]))) == 2){
        op <- ranger::importance_pvalues(rout, "janitza")[,2]
      }
      else{
        warning("No p-values calcuated. When a quantitative response or more than 2 categorical response options are present, p-values must be done via permutation. We suggest running this with ranger directly. See ?ranger::importance_pvalues for help.\n")
      }
    }


    # prepare objects for return
    if(!is.null(.formula)){
      imp.out <- cbind(sub.x@snp.meta, rout$variable.importance[which(!names(rout$variable.importance) %in% cvars)])
      colnames(imp.out)[ncol(imp.out)] <- paste0(response, "_", "RF_importance")
      covariate_importance <- data.frame(variable = cvars, importance = rout$variable.importance[which(names(rout$variable.importance) %in% cvars)])
      if(exists("op")){
        covariate_importance$p_val <- op[which(names(op) %in% cvars)]
        imp.out <- cbind(imp.out, op[which(!names(op) %in% cvars)])
        colnames(imp.out)[ncol(imp.out)] <- paste0(response, "_", "RF_importance_pvals")
        
      }
      
    }
    else{
      imp.out <- cbind(sub.x@snp.meta, rout$variable.importance)
      colnames(imp.out)[ncol(imp.out)] <- paste0(response, "_", "RF_importance")
      if(exists("op")){
        imp.out <- cbind(imp.out, op[which(!names(op) %in% cvars)])
        colnames(imp.out)[ncol(imp.out)] <- paste0(response, "_", "RF_importance_pvals")
      }
    }

    pred.out <- data.frame(predicted = rout$predictions, pheno = sub.x@sample.meta[,response],
                           stringsAsFactors = F)

    if(!is.null(.formula)){
      return(list(importance = imp.out, predictions = pred.out, model = rout, covariate_importance = covariate_importance))
    }
    else{
      return(list(importance = imp.out, predictions = pred.out, model = rout))
    }
  }


  if(!all(facets %in% x@facets)){
    .make_it_quiet(x <- .add.facets.snpR.data(x, facets))
  }

  facets <- .check.snpR.facet.request(x, facets)
  # run for each task
  tl <- .get.task.list(x, facets)
  out <- vector("list", nrow(tl))
  for(i in 1:nrow(tl)){
    out[[i]] <- run_ranger(suppressWarnings(.subset_snpR_data(x, facets = tl[i,1], subfacets = tl[i,2], snp.facets = tl[i,3], snp.subfacets = tl[i,4])),
                           .formula = .formula, par = par, ...)
    out[[i]]$.fm <- cbind(facet = tl[i,1], subfacet = tl[i,2], row.names = NULL)
  }
  
  # process
  models <- vector("list", length(out))
  for(i in 1:length(out)){
    type <- ifelse(out[[i]]$.fm[1] == ".base", ".base", "sample")
    out[[i]]$importance <- cbind(out[[i]]$.fm, facet.type = type, out[[i]]$importance, row.names = NULL)
    models[[i]] <- out[[i]][-which(names(out[[i]]) %in% c("importance", ".fm"))]
    names(models)[i] <- paste0(out[[i]]$.fm[1,], collapse = "_")
  }
  
  stats <- data.table::rbindlist(purrr::map(out, "importance"))
  x <- .merge.snpR.stats(x, stats)
  x <- .update_citations(x, "Wright2017", "Random Forest", paste0("Random Forest test against ", response, "."))
  

  return(list(x = x, models = models))
}





