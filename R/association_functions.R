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
#' @param charcter, default set to response ~ 1. Null formula for the response variable, as described in \code{\link[stats]{formula}}.
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
  if(grepl("\\.", response)[1]){
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
        bad.cvars <- which(!(cvars %in% colnames(sub.x@sample.meta)))
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
      odds <- log((a/b)/(c/d))
      return(data.frame(log_odds_ratio = odds, se = log(s.e)))
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



run_random_forest <- function(x, facets = NULL, response, formula = NULL,
                              num.trees = 10000, mtry = NULL,
                              importance = "impurity",
                              interpolate = "bernoulli",
                              cpar = FALSE, par = FALSE, ...){
  run_ranger <- function(sub.x, ...){

    #================grab data=====================
    ## sn format
    sn <- format_snps(sub.x, "sn", interpolate = interpolate)
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

        sn <- cbind(sub.x@sample.meta[,response], sub.x@sample.meta[,cvars], sn)

      }
    }
    else{
      sn <- cbind(sub.x@sample.meta[,response], sn)
    }
    colnames(sn)[1] <- response


    #================run the model:=====================
    if(is.null(formula)){
      rout <- ranger::ranger(dependent.variable.name = response,
                             data = sn,
                             num.trees = num.trees,
                             mtry = mtry,
                             importance = importance)
    }
    else{
      rout <- ranger::ranger(formula = formula,
                             data = sn,
                             num.trees = num.trees,
                             mtry = mtry,
                             importance = importance)
    }

    #===============make sense of output=================
    imp.out <- cbind(sub.x@snp.meta, rout$variable.importance)
    colnames(imp.out)[ncol(imp.out)] <- paste0(response, "_", "RF_importance")

    pred.out <- data.frame(predicted = rout$predictions, pheno = sub.x@sample.meta[,response],
                           stringsAsFactors = F)

    return(list(effects = effects, predictions = pred.out, model = rout))
  }


  if(!all(facets %in% x@facets)){
    invisible(capture.output(x <- add.facets.snpR.data(x, facets)))
  }

  facets <- check.snpR.facet.request(x, facets)

  out <- apply.snpR.facets(x, facets, req = "snpRdata", case = "ps", fun = run_ranger, response = response, formula = formula,
                           num.trees = num.trees, mtry = mtry, importance = importance, interpolate = interpolate,  par = par, ...)

  return(merge.snpR.stats(x, out))
}

