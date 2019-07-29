#' Interface with BGLR to run genomic prediction with snpRdata objects.
#'
#' Run genomic prediction given a single response variable (usually a phenotype)
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
#' interest should be selected with \code{\link{subset_snpR_data}}. See
#' examples.
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
#' @param model character, default "BayesB". Prediction model to use, see
#'   description for the ETA argument in \code{\link[BGLR]{BGLR}}.
#'
#' @export
#' @author William Hemstrom
#' @references Pérez, P., and de los Campos, G. (2014). \emph{Genetics}.
#' @return A list containing: \itemize{ \item{model: } The model output from
#'   BGLR. See \code{\link[BGLR]{BGLR}}. \item{h2: } Estimated heritability of
#'   the response variable. \item{transposed_interpolated_genotypes: }
#'   Transposed genotypes, the input for BGLR. \item{predictions: } A data.frame
#'   containing the provided phenotypes and the predicted Breeding Values (BVs)
#'   for those phenotypes. }
#'
#' @examples
#' # run and plot a basic prediction
#' ## add some dummy phenotypic data.
#' sample.meta <- cbind(weight = rnorm(ncol(stickSNPs)), stickSNPs@sample.meta)
#' dat <- import.snpR.data(as.data.frame(stickSNPs), stickSNPs@snp.meta, sample.meta)
#' ## run prediction
#' gp <- run_genomic_prediction(dat, response = "weight", iterations = 1000, burn_in = 100, thin = 10)
#' ## dummy phenotypes vs. predicted Breeding Values for those dummy predictions.
#' with(gp$predictions, plot(phenotype, predicted_BV))
#'
run_genomic_prediction <- function(x, response, iterations,
                                   burn_in, thin,
                                   model = "BayesB"){

  # get a properly formatted input file
  if(length(x@sn) != 0){
    if(x@sn$type != "bernoulli"){
      suppressWarnings(x@sn$sn <- format_snps(x, "sn", interpolate = "bernoulli"))
      x@sn$type <- "bernoulli"
    }
  }
  else{
    suppressWarnings(x@sn$sn <- format_snps(x, "sn", interpolate = "bernoulli"))
    x@sn$type <- "bernoulli"
  }

  sn <- x@sn$sn
  sn <- sn[,-c(1:(ncol(x@snp.meta) - 1))]

  # prepare the BGLR input
  sn <- t(sn)
  phenotypes <- x@sample.meta[,response]
  if(!is.numeric(phenotypes)){
    phenotypes <- as.numeric(as.factor(phenotypes)) - 1
    if(length(unique(phenotypes)) != 2){
      stop("Only two unique phenotypes allowed if not quantitative.\n")
    }
  }
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
#' Run genomic prediction given a single response variable (usually a phenotype)
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
#' interest should be selected with \code{\link{subset_snpR_data}}. See
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
#' cross_validate_genomic_prediction(dat, response = "weight", iterations = 1000, burn_in = 100, thin = 10)
#'
cross_validate_genomic_prediction <- function(x, response, iterations = 10000,
                                              burn_in = 1000, thin = 100, cross_percentage = 0.9,
                                              model = "BayesB", cross_samples = NULL, plot = TRUE){

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
  capture.output(suppressWarnings(suppressMessages(sub.x <- subset_snpR_data(x, samps = model.samples))))
  model <- run_genomic_prediction(sub.x, response = response, iterations =  iterations,
                                  burn_in = burn_in, thin = thin, model = model)

  # check accuracy
  cross.samples <- (1:ncol(x))[-sort(model.samples)]
  if(length(x@sn) != 0){
    if(x@sn$type != "bernoulli"){
      suppressWarnings(x@sn$sn <- format_snps(x, "sn", interpolate = "bernoulli"))
      x@sn$type <- "bernoulli"
    }
  }
  else{
    suppressWarnings(x@sn$sn <- format_snps(x, "sn", interpolate = "bernoulli"))
    x@sn$type <- "bernoulli"
  }
  whole.sn <- x@sn$sn
  whole.sn <- whole.sn[,-(1:(ncol(x@snp.meta) - 1))]
  whole.sn <- interpolate_sn(whole.sn)
  new.sn <- whole.sn[,cross.samples]
  new.sn <- t(new.sn)
  pred.phenos <- new.sn%*%model$model$ETA[[1]]$b
  pdat <- as.data.frame(cbind(observed = x@sample.meta[cross.samples, response], predicted = pred.phenos), stringsAsFactors = F)
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

#' Run case/control or quantitative association tests on SNP data.
#'
#' Runs several different association tests on SNP data. The response variable
#' must have only two different categories (as in case/control) for most test types,
#' although the "gmmat.score" method supports quantitative traits. Tests may be
#' broken up by sample-specific facets.
#'
#' Several methods can be used: Armitage, chi-squared, and odds ratio. For The
#' Armitage approach weights should be provided to the "w" argument, which
#' specifies the weight for each possible genotype (homozygote 1, heterozygote,
#' homozygote 2). The default, c(0,1,2), specifies an addative  The "gmmat.score"
#' method uses the approach described in Chen et al. (2016) and implemented in
#' the \code{\link[GMMAT]{glmmkin}} and \code{\link[GMMAT]{glmm.score}} functions.
#' For this method, a 'G' genetic relatedness matrix is first created using the
#' \code{\link[AGHmatrix]{Gmatrix}} function according to Yang et al 2010.
#'
#' Multi-category data is currently not supported, but is under
#' development.
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
#'   Options: \itemize{ \item{gmmat.score: }  Population/family structure corrected mlm approach, based on
#'   Chen et al (2016). \item{armitage: } Armitage association test, based on
#'   Armitage (1955). \item{odds_ratio: } Log odds ratio test. \item{chisq: }
#'   Chi-squared test. } See description for more details.
#' @param w numeric, default c(0, 1, 2). Weight variable for each genotype for
#'   the Armitage association method. See description for details.
#' @param formula charcter, default set to response ~ 1. Null formula for the response variable, as described in \code{\link[stats]{formula}}.
#' @param family.override character, default NULL.
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
#' @references Chen et al. (2016). Control for Population Structure and Relatedness for
#'   Binary Traits in Genetic Association Studies via Logistic Mixed Models. \emph{American Journal of Human Genetics}.
#' @references Yang et al. (2010). Common SNPs explain a large proportion of the heritability for human height.
#'   \emph{Nature Genetics}.
#'
#' @examples
#'   # add a dummy phenotype
#'   sample.meta <- cbind(stickSNPs@sample.meta, phenotype = sample(c("A", "B"), nrow(stickSNPs@sample.meta), T))
#'   x <- import.snpR.data(as.data.frame(stickSNPs), stickSNPs@snp.meta, sample.meta)
#'   calc_association(x, facets = c("pop", "pop.fam"), response = "phenotype", method = "armitage")
#'
calc_association <- function(x, facets = NULL, response, method = "gmmat.score", w = c(0,1,2),
                             formula = NULL, family.override = FALSE, maxiter = 500, sampleID = NULL, maf = 0.05, par = FALSE, ...){
  #==============sanity checks===========
  # response
  msg <- character()
  if(length(response) != 1){
    msg <- c(msg,
             paste0("Only one response variable permitted."))
  }
  if(grepl("(?<!^)\\.", response)[1]){
    msg <- c(msg,
             paste0("Only one sample-specific category allowed (e.g. pop but not fam.pop)."))
  }
  if(!response[1] %in% x@sample.meta){
    msg <- c(msg,
             paste0("Response variable must be present in sample metadata."))
  }
  else(
    if(length(unique(x@sample.meta[,response])) != 2){
      if(method %in% c("armitage", "odds_ratio", "chisq")){
        msg <- c(msg,
                 paste0("Only two categories allowed for response variable for method: ", method, "."))
      }
    }
  )

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
      if(class(res) == "try-error"){
        msg <- c(msg,
                 "formula must be a valid formula. Type ?formula for help.\n")
      }
      else{
        phen <- strsplit(formula, "~")
        phen <- phen[[1]][1]
        phen <- gsub(" ", "", phen)
        if(phen != response){
          msg <- c(msg,
                   "The response variable in the provided formula must be the same as that provided to the response argument.\n")
        }

        # see which covariates we need
        cvars <- terms(formula(formula))
        cvars <- attr(cvars, "term.labels")
        bad.cvars <- which(!(cvars %in% colnames(x@sample.meta)))
        if(length(bad.cvars) > 0){
          msg <- c(msg,
                   "Some covariates specified in formula not found in sample metadata: ", paste0(cvars[bad.cvars], collapse = ", "), ".\n")
        }


      }
    }
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
      chi.p <- pchisq(chi.stat, 1, lower.tail = F)
      return(data.frame(chi_stat = chi.stat, chi_p = chi.p, associated_allele = asso.allele))
    }
  }
  run_gmmat <- function(sub.x, form, iter, sampleID, family.override, ...){
    # sn format
    sn <- format_snps(sub.x, "sn", interpolate = FALSE)
    sn <- sn[,-c(which(colnames(sn) %in% colnames(sub.x@snp.meta)))]

    ## G matrix
    G <- AGHmatrix::Gmatrix(t(sn), missingValue = NA, method = "Yang")
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
      phenos[,response] <- as.numeric(as.factor(phenos[,response])) - 1 # set the phenotype to 0,1,ect
    }
    ## family
    if(family.override != FALSE){
      family <- family.override
    }
    else{
      if(length(unique(phenos[,response])) > 2){
        family <- gaussian(link = "identity")
      }
      else if(length(unique(phenos[,response])) == 2){
        family <- binomial(link = "logit")
      }
      else{
        stop("Less than two unique phenotypes.\n")
      }
    }


    # run null model
    invisible(capture.output(mod <- GMMAT::glmmkin(fixed = formula,
                                                   data = phenos,
                                                   kins = G,
                                                   id = sampleID,
                                                   family = family)))

    # run the test
    nmeta.col <- 2 + ncol(sub.x@snp.meta)
    write.table(asso.in, "asso_in.txt", sep = "\t", quote = F, col.names = F, row.names = F)
    score.out <- GMMAT::glmm.score(mod, "asso_in.txt",
                                   "asso_out_score.txt",
                                   infile.ncol.skip = nmeta.col,
                                   infile.ncol.print = 1:nmeta.col,
                                   infile.header.print = colnames(asso.in)[1:nmeta.col])
    score.out <- read.table("asso_out_score.txt", header = T, stringsAsFactors = F)

    file.remove(c("asso_in.txt", "asso_out_score.txt"))
    return(score.out)
  }

  #==============run the function=========
  # check facets
  facets <- check.snpR.facet.request(x, facets)
  if(!all(facets %in% x@facets)){
    invisible(capture.output(x <- add.facets.snpR.data(x, facets)))
  }

  if(method == "armitage"){
    out <- apply.snpR.facets(x, facets = facets, req = "cast.gs", case = "ps", fun = calc_armitage, response = response, w = w)
  }
  else if(method == "odds_ratio" | method == "chisq"){
    out <- apply.snpR.facets(x, facets = facets, req = "cast.ac", case = "ps", fun = odds.ratio.chisq, response = response, method = method)
  }
  else if(method == "gmmat.score"){
    out <- apply.snpR.facets(x, facets = facets, req = "snpRdata", case = "ps", maf = maf, fun = run_gmmat, response = response, form = formula, iter = iter, sampleID = sampleID, family.override = family.override, ...)
  }

  x <- merge.snpR.stats(x, out)
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
#' \code{\link[snpR]{Facets_in_snpR}}. Since RF models are calculated without
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
#' In general, random forest parameters should be tuned so as to reduce the out-of-bag error rates (OOB-ER).
#' This value is visible in the returned object under the model lists. Simply calling a specific model will
#' output the OOB-ER, and they are also stored under the 'prediction.error' name in the model. For details
#' on tuning RF models, we recommend Goldstein et al. (2011).
#'
#' For more detail on the random forest model and ranger arguments, see
#' \code{\link[ranger]{ranger}}.
#'
#' @param x snpRdata object.
#' @param facets character, default NULL. Categorical metadata variables by
#'   which to break up analysis. See \code{\link{Facets_in_snpR}} for more
#'   details.
#' @param response character. Name of the column containing the response
#'   variable of interest. Must match a column name in sample metadata. Response
#' must be categorical, with only two categories.
#' @param formula charcter, default NULL. Model for the response variable, as
#'   described in \code{\link[stats]{formula}}. If NULL, the model will be
#'   equivalent to response ~ 1.
#' @param num.trees numeric, default 10000. Number of trees to grow. Higher
#'   numbers will increase model accuracy, but increase calculation time. See
#'   \code{\link[ranger]{ranger}} for details.
#' @param mtry numeric, default is the square root of the number of SNPs. Number
#'   of variables (SNPs) by which to split each node. See
#'   \code{\link[ranger]{ranger}} for details.
#' @param importance character, default "impurity_corrected". The method by
#'   which SNP importance is determined. Options: \itemize{\item{impurity}
#' \item{impurity_corrected} \item{permutation}}. See
#' \code{\link[ranger]{ranger}} for details.
#' @param interpolate character, default "bernoulli". Interpolation method for
#'   missing data. Options: \itemize{\item{bernoulli: }binomial draws for the
#'   minor allele. \item{af: } insertion of the average allele frequency}.
#' @param pvals logical, default TRUE. Determines if pvalues should be
#'   calculated for importance values. If the response variable is quantitative,
#'   no pvalues will be returned, since they must be calculated via permutation
#'   and this is very slow. For details, see
#'   \code{\link[ranger]{importance_pvalues}}.
#' @param par numeric, default FALSE. Number of parallel computing cores to use
#'   for computing RFs across multiple facet levels or within a single facet if
#'   only a single category is run (either a one-category facet or no facet).
#' @param ... Additional arguments passed to \code{\link[ranger]{ranger}}.
#'
#' @return A list containing: \itemize{\item{data: } A snpRdata object with RF
#'   importance values merged in to the stats slot. \item{models: } A named list
#' containing both the models and data.frames containing the predictions vs
#' observed phenotypes.}
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
run_random_forest <- function(x, facets = NULL, response, formula = NULL,
                              num.trees = 10000, mtry = NULL,
                              importance = "impurity_corrected",
                              interpolate = "bernoulli", pvals = TRUE, par = FALSE, ...){
  run_ranger <- function(sub.x, opts.list, ...){

    #================grab data=====================
    ## sn format
    if(length(sub.x@sn) != 0){
      if(sub.x@sn$type != interpolate){
        sn <- format_snps(sub.x, "sn", interpolate = interpolate)
      }
      else{
        sn <- sub.x@sn$sn
      }
    }
    else{
      sn <- format_snps(sub.x, "sn", interpolate = interpolate)
    }
    sn <- sn[,-c(which(colnames(sn) %in% colnames(sub.x@snp.meta)))]
    sn <- t(sn)


    ## attach phenotype, anything else in the formula if given.
    if(!is.null(formula)){
      res <- try(formula(formula), silent = T)
      if(class(res) == "try-error"){
        msg <- c(msg,
                 "formula must be a valid formula. Type ?formula for help.\n")
      }
      else{
        phen <- strsplit(formula, "~")
        phen <- phen[[1]][1]
        phen <- gsub(" ", "", phen)
        if(phen != response){
          msg <- c(msg,
                   "The response variable in the provided formula must be the same as that provided to the response argument.\n")
        }

        # see which covariates we need
        cvars <- terms(formula(formula))
        cvars <- attr(cvars, "term.labels")
        bad.cvars <- which(!(cvars %in% colnames(sub.x@sample.meta)))
        if(length(bad.cvars) > 0){
          msg <- c(msg,
                   "Some covariates specified in formula not found in sample metadata: ", paste0(cvars[bad.cvars], collapse = ", "), ".\n")
        }

        ocn <- colnames(sn)
        sn <- cbind(sub.x@sample.meta[,response], sub.x@sample.meta[,cvars], sn)

        # reset formula:
        formula <- paste0(formula, paste0(ocn, collapse = " + "))

      }
    }
    else{
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
      tpar <- NULL
    }
    else{
      if(nrow(opts.list) == 1){
        tpar <- par
      }
      else{
        tpar <- NULL
      }
    }

    # run with or without formula
    if(is.null(formula)){
      rout <- ranger::ranger(dependent.variable.name = response,
                             data = sn,
                             num.trees = num.trees,
                             mtry = mtry,
                             importance = importance, num.threads = tpar, verbose = T, ...)
    }
    else{
      rout <- ranger::ranger(formula = formula,
                             data = sn,
                             num.trees = num.trees,
                             mtry = mtry,
                             importance = importance, num.threads = tpar, verbose = T, ...)
    }

    #===============make sense of output=================
    # get p-values:
    if(pvals){
      if(length(unique(as.character(x@sample.meta[,response]))) == 2){
        op <- ranger::importance_pvalues(rout, "janitza")[,2]
      }
      else{
        warning("No p-values calcuated. P-value must be done via permutation when a quantitative response is used. See ?ranger::importance_pvalues for help. Imput data for these p-value caluculations can be generated with format_snps using the 'sn' format option.")
      }
    }


    # prepare objects for return
    imp.out <- cbind(sub.x@snp.meta, rout$variable.importance)
    colnames(imp.out)[ncol(imp.out)] <- paste0(response, "_", "RF_importance")
    if(exists("op")){
      imp.out <- cbind(imp.out, op)
      colnames(imp.out)[ncol(imp.out)] <- paste0(response, "_", "RF_importance_pvals")
    }

    pred.out <- data.frame(predicted = rout$predictions, pheno = sub.x@sample.meta[,response],
                           stringsAsFactors = F)

    return(list(importance = imp.out, predictions = pred.out, model = rout))
  }


  if(!all(facets %in% x@facets)){
    invisible(capture.output(x <- add.facets.snpR.data(x, facets)))
  }

  facets <- check.snpR.facet.request(x, facets)

  out <- apply.snpR.facets(x, facets, req = "snpRdata", case = "ps", fun = run_ranger, response = response, formula = formula,
                           interpolate = interpolate,  par = par, ...)

  x <- merge.snpR.stats(x, out$stats)

  return(list(data = x, models = out$models))
}


# function to run:
refine_rf <- function(rf, facets = NULL, response, subfacet = NULL, formula = NULL, num.trees = 10000, trim_cuttoffs = NULL,
                      trim = 0.5, par = FALSE,
                      interpolate = "bernoulli", ...){
  #==============sanity checks================
  msg <- character()
  if(length(facets) > 1){
    msg <- c(msg, "No more than one facet may be refined at a time.\n")
    if(is.null(subfacet)[1]){
      msg <- c(msg, "No more that one subfacet may be refined at a time. Please list a subfacet.\n")
    }
  }
  if(!is.null(subfacet[1])){
    if(length(subfacet) > 1){
      msg <- c(msg, "No more that one subfacet may be refined at a time.\n")
    }
  }
  if(!is.null(facets)){
    if(length(rf$models) != 1 | names(rf$models[[1]] != ".base_.base")){
      msg <- c(msg, "If the base facet is requested, please provide an input rf result with only the base facet calculated.\n")
    }
  }

  #==============the actual refinement function==================
  run_refine <- function(rf, response, formula, num.trees, trim_cuttoffs, trim, search_cuttoff, par,
                         ...){

    #=============subfunctions:=========
    # trim snps below a given importance quantile
    remove_trim <- function(imps, trim){
      best.imp <- quantile(imps, trim)
      best.imp <- which(imps >= best.imp[1])
      return(best.imp)
    }
    # trim to a specific number of snps
    remove_set_to_n_snps <- function(imps, n_snps){
      imps <- data.frame(imps, ord = 1:length(imps))
      imps <- dplyr::arrange(imps, desc(imps))
      imps <- imps[1:n_snps,]
      imps <- dplyr::arrange(imps, ord)
      best.imp <- imps$ord
      return(best.imp)
    }
    # remove a single snp
    remove_single <- function(imps){
      best.imp <- which.min(imps)
      best.imp <- (1:nrow(rf$data@stats))[-best.imp]
      return(best.imp)
    }
    # set to the next level, then trim a given percentage
    remove_set_and_trim <- function(imps, n_snps, trim){
      # set
      imps <- data.frame(imps, ord = 1:length(imps))
      imps <- dplyr::arrange(imps, desc(imps))
      imps <- imps[1:n_snps,]
      imps <- dplyr::arrange(imps, ord)

      # trim
      best.imp <- quantile(imps[,1], trim)
      best.imp <- imps[which(imps[,1] >= best.imp[1]),]

      # return
      best.imp$ord
    }

    # remove according to cuttoffs and trim levels
    remove_snps <- function(imps, trim_cuttoffs, trim){

      # figure out which "bin" we are in
      n_snps <- length(imps)
      bin.logi <- which(n_snps > trim_cuttoffs)
      if(length(bin.logi) > 0){

        # grab the current bin and trim
        bin <- min(bin.logi)
        best.imp <- remove_trim(imps, trim[bin])

        # check if we changed bins!
        bin.t.trim <- which(length(best.imp) > trim_cuttoffs)
        if(!identical(bin.t.trim, bin.logi)){

          # if we haven't hit the last bin
          if(length(bin.t.trim) != 0){
            # we either want to set this to - the current trim or to the cuttoff - the next trim, whichever has more snps.

            best.imp.current.trim <- remove_trim(imps, trim[bin])
            best.imp.set.plus.next.trim <- remove_set_and_trim(imps, trim_cuttoffs[min(bin.t.trim) - 1], trim[min(bin.t.trim)])
            if(length(best.imp.current.trim) > length(best.imp.set.plus.next.trim)){
              best.imp <- best.imp.current.trim
            }
            else{
              best.imp <- best.imp.set.plus.next.trim
            }
          }


          ## if we've hit the last bin and are using a single cutoff, set to single cuttoff
          else if(length(bin.t.trim) == 0 & trim_single){
            best.imp <- remove_set_to_n_snps(imps, trim_cuttoffs[length(trim_cuttoffs)])
            if(length(best.imp) == n_snps){
              best.imp <- remove_single(imps)
            }
          }

          ## if we've hit the last bin but aren't using a single cuttoff, keep trimming as usual
          else if(length(bin.t.trim) == 0){
            best.imp.current.trim <- remove_trim(imps, trim[bin])
            best.imp.set.plus.next.trim <- remove_set_and_trim(imps, trim_cuttoffs[bin],  trim[length(trim)])
            if(length(best.imp.current.trim) > length(best.imp.set.plus.next.trim)){
              best.imp <- best.imp.current.trim
            }
            else{
              best.imp <- best.imp.set.plus.next.trim
            }

          }
        }
      }
      # if we've hit the last bin...
      else{
        if(trim_single){
          best.imp <- remove_single(imps)
        }
        else{
          best.imp <- remove_trim(imps, trim[length(trim)])
        }
      }



      return(best.imp)
    }

    #==========initialize============



    # intialize output
    out <- matrix(NA, 1000, 2)
    colnames(out) <- c("n_snps", "prediction_error")
    out[1,] <- c(nrow(rf$data), rf$models[[1]]$model$prediction.error)

    # initialize output confusion array
    if(any(names(rf$models[[1]]$model) == "confusion.matrix")){
      conf.out <- array(NA, dim = c(dim(rf$models[[1]]$model$confusion.matrix), 1000))
      conf.out[,,1] <- rf$models[[1]]$model$confusion.matrix
    }


    # initialize best model output
    best.mod <- rf
    best.error <- out[1,2]

    # intialize difference
    diff <- 1

    # find the target column name
    tar.col <- paste0(response, "_RF_importance")

    # fix the single run_cuttoff if NULL
    if(is.null(trim_cuttoffs)){
      trim_cuttoffs <- 1
    }

    # if we are provided with less trim_cuttoffs than trim levels, we assume no single trimming
    if(length(trim_cuttoffs) < length(trim)){
      trim_single <- FALSE
    }
    # if the same, we trim single SNPs beneath the lowest cutoff
    else if (length(trim_cuttoffs) == length(trim)){
      trim_single <- TRUE
    }
    else{
      stop("Not enough trim levels provided for given trim cuttoffs.\n")
    }

    #==========while loop================
    i <- 2
    continue <- T
    while(continue == T){

      # if we somehow reach the end of the output storage, make it bigger!
      if(i == nrow(out)){
        out <- rbind(out, matrix(NA, 100000000, 2))
        conf.out <- c(conf.out, array(NA, dim = c(dim(rf$models[[1]]$model$confusion.matrix), 1000)))
        conf.out <- array(conf.out, c(dim(rf$models[[1]]$model$confusion.matrix), 2000))
      }

      imps <- abs(rf$data@stats[[tar.col]])

      # get the snps to keep
      best.imp <- remove_snps(imps, trim_cuttoffs, trim)

      # subset
      suppressWarnings(input <- subset_snpR_data(rf$data, snps = best.imp))
      if(nrow(input) == 1){
        continue <- FALSE
      }

      # report some things
      out[i,1] <- nrow(input)
      cat("Refinement: ", i - 1, "\n\tStarting prediction error: ", out[i-1, 2],
          "\n\tNumber of remaining SNPs: ", out[i,1], "\n\tBeginning rf...\n")

      # run rf
      suppressWarnings(rf <- run_random_forest(x = input, response = response, formula = formula, num.trees = num.trees, mtry = nrow(input), par = par, pvals = FALSE, ...))

      # save outputs
      out[i,2] <- rf$models[[1]]$model$prediction.error
      if(exists("conf.out")){
        conf.out[,,i] <- rf$models[[1]]$model$confusion.matrix
      }
      if(out[i,2] < best.error){
        best.mod <- rf
        best.error <- out[i,2]
      }

      i <- i + 1
    }

    #==============return==============
    # return final model and outputs
    empties <- which(is.na(out[,1]))
    out <- out[-empties,]
    if(exists("conf.out")){
      conf.out <- conf.out[,,-empties]
      return(list(refined_model = rf, error_delta = out, confusion_matrices = conf.out, best_model = best.mod))
    }
    else{
      return(list(refined_model = rf, error_delta = out, best_model = best.mod))

    }
  }

  #==============prep============================================
  # strip down to the correct facet, doing  more sanity checks
  facets <- check.snpR.facet.request(rf$data, facets, "snp")
  o.facet <- facets

  # subset if need be
  if(o.facet != ".base"){
    facets <- strsplit(facets, "(?<!^)\\.", perl = T)
    facets <- unlist(facets)

    # figure out which of the run rf models (facets) we are using.
    use.model <- sapply(facets, function(x) grepl(pattern = x, x = names(rf$models)))
    if(length(names(rf$models)) == 1){
      if(sum(use.model) != length(use.model)){
        msg <- c(msg, "Requested facet not found in provided rf.\n")
      }
      else{
        use.model <- 1
      }
    }
    else{
      use.model <- which(rowSums(use.model) == ncol(use.model))
      if(length(use.model) == 0){
        msg <- c(msg, "Requested facet not found in provided rf.\n")
      }
    }


    # figure out which subfacet we are using if requested.
    if(!is.null(subfacet)){
      good.mods <- names(rf$models)[use.model]
      use.model <- which(grepl(subfacet, good.mods))

      if(length(use.model) > 1){
        msg <- c(msg, "Subfacet + facet request matches more than one model in the provided rf.\n")
      }
      else if(length(use.model)  == 0){
        msg <- c(msg, "Subfacet + facet request matches no models in the provided rf.\n")
      }

      else{
        use.model <- good.mods[use.model]
      }
    }

    # stop if errors
    if(length(msg) > 0){
      stop(msg)
    }

    # subset
    o.mod <- rf$models[[use.model]]
    dat <- subset_snpR_data(rf$data, facets = o.facet, subfacets = subfacet)
    dat <- import.snpR.data(as.data.frame(dat, stringsAsFactors = F), dat@snp.meta, dat@sample.meta, dat@mDat)
    matches <- intersect(which(rf$data@stats$facet == o.facet),
                         which(rf$data@stats$subfacet == subfacet))
    imps <- rf$data@stats[matches,]
    imps$facet <- ".base"
    imps$subfacet <- ".base"
    dat <- merge.snpR.stats(dat, imps)

    rf <- list(data = dat, models = list(.base_.base = o.mod))
  }


  # add sn formatted data if not present
  if(length(rf$data@sn) != 0){
    if(rf$data@sn$type != interpolate){
      sn <- format_snps(rf$data, "sn", interpolate = interpolate)
      rf$data@sn <- list(type = interpolate, sn = sn)

    }
  }
  else{
    sn <- format_snps(rf$data, "sn", interpolate = interpolate)
    rf$data@sn <- list(type = interpolate, sn = sn)
  }

  #==============run the refinement============

  out <- run_refine(rf, response = response,
                    formula = formula,
                    num.trees = num.trees,
                    trim_cuttoffs = trim_cuttoffs,
                    trim = trim,
                    search_cuttoff = search_cuttoff,
                    par = par,
                    ...)

  return(out)
}


