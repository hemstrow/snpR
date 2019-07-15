#'Create a heatmap from pairwise linkage data.
#'
#'Prepares a ggplot2 heatmap from pairwise linkage disequilibrium data stored in
#'a snpRdata object.
#'
#'Since the output is a ggplot object, options can be added or changed by adding
#'"+ function()" to the end of \code{LD_pairwise_heatmap}. Some common options
#'are also built into this function as agruments, but can be overwritten freely.
#'
#'Specific facets and facet levels can be requested as long as LD data for those
#'facets has been previously calculated. Facets should be provided as described
#'in \code{\link{Facets_in_snpR}}. Only one facet can be plotted at once. Plots
#'for individual SNP or sample-specific facet levels for the specified facet can
#'also be requested. See examples. For facets with many categories, this is
#'strongly recommended, since plots with many facet levels become
#'computationally challenging to produce and difficult to interpret. A facets
#'argument of NULL will plot the base level facet, a facet argument of "all" is
#'not accepted.
#'
#'Note that NA LD values will be colored white in the resulting plot.
#'
#'@param x snpRdata object.
#'@param facets character, default NULL. Categorical metadata variables by which
#'  to break up plots. Must match facets for which LD data has been previously
#'  calculated, and only a single facet can be plotted at once. See
#'  \code{\link{Facets_in_snpR}} for more details.
#'@param snp.subfacet character, default NULL. Specific snp-specific levels of
#'  the provided facet to plot. See examples.
#'@param sample.subfacet character, default NULL. Specific sample-specific
#'  levels of the provided facet to plot. See examples.
#'@param LD_measure character, default rsq. LD metric to plot. Must be present
#'  in the calculated LD data.
#'@param r Numeric. Region of the chromosome to subset and plot. Given in kb in
#'  the format numeric vector c(lower, upper).
#'@param l.text character, default "rsq". Legend title.
#'@param viridis.option character, default "inferno". Viridis color scale option to use.
#'  Other color scales may be subsituted by appending the scale_color_continuous
#'  and scale_fill_continuous ggplot functions to the produced plot using the
#'  '+' operator. See \code{\link[ggplot2]{scale_colour_gradient}} for details.
#'@param title character. Plot title.
#'@param t.sizes numeric, default c(16, 13, 10, 12, 10). Text sizes, given as
#'  c(title, legend.title, legend.ticks, axis, axis.ticks).
#'
#'@return A list containing: \itemize{ \item{plot: } A pairwise LD heatmap as a
#'  ggplot object. \item{dat: } Data used to generate the ggplot object. }
#'
#'@author William Hemstrom
#'@author Nicholas Sard
#'@export
#'
#' @examples
#' # get LD data
#' dat <- calc_pairwise_ld(stickSNPs, c("pop.group"))
#'
#' # produce plots for linkage group IX in the ASP and CLF populations.
#' plot_pairwise_LD_heatmap(dat3, c("pop.group"), "groupIX", c("ASP", "CLF"))
#'
#' # produce plots for every population for linkage group IV
#' plot_pairwise_LD_heatmap(dat3, c("pop.group"), "groupIV")
#'
#'
plot_pairwise_LD_heatmap <- function(x, facets = NULL, snp.subfacet = NULL, sample.subfacet = NULL, LD_measure = "rsq", r = NULL,
                                     l.text = "rsq", viridis.option = "inferno",
                                     title = NULL, t.sizes = c(16, 13, 10, 12, 10),
                                     background = "white"){
  #==============sanity checks===========
  msg <- character()

  if(length(x@pairwise.LD) == 0){
    stop("No LD data found.\n")
  }

  if(is.null(facets)){facets <- ".base"}

  # check facets
  snp.facet <- check.snpR.facet.request(x, facets, remove.type = "sample")[[1]]
  sample.facet <- check.snpR.facet.request(x, facets, remove.type = "snp")[[1]]
  facets <- check.snpR.facet.request(x, facets, remove.type = "none", return.type = T)
  facet.type <- facets[[2]]
  facets <- facets[[1]]
  bad.facets <- facets[!(facets %in% names(x@pairwise.LD$LD_matrices))]

  if(length(bad.facets) > 0){
    msg <- c(msg, paste0("LD data not found for some facets: ", paste0(bad.facets, collapse = ", "), "."))
  }

  if(length(facets) > 1){
    msg <- c(msg, "Only one facet may be plotted at once.")
  }


  # check LD measure
  good.ld.measures <- c("rsq", "Dprime", "pval")
  if(length(LD_measure) != 1){
    msg <- c(msg, "Only one LD measure may be plotted at once.")
  }
  if(!(LD_measure %in% good.ld.measures)){
    msg <- c(msg, paste0("Unaccepted LD measure. Accepted measures:", paste0(good.ld.measures, collapse = ", "), "."))
  }


  # check snp.subfacet
  if(!is.null(snp.subfacet)){
    if(length(snp.facet) == 0 & !is.null(snp.subfacet)){
      msg <- c(msg, "SNP subfacet graph requested, but no snp or complex facets listed in facets.")
    }
    else if(!(all(snp.subfacet %in% names(x@pairwise.LD$LD_matrices[[facets]][[1]])))){
      msg <- c(msg, "Requested SNP subfacet not located in possible subfacets. SNP subfacet may not be in the provided SNP metadata.")
    }
  }

  # check sample.subfacet
  if(!is.null(sample.subfacet)){
    if(length(sample.facet) == 0 & !is.null(sample.subfacet)){
      msg <- c(msg, "Sample subfacet graph requested, but no sample facets listed in facets.")
    }
    else if(!(all(sample.subfacet %in% names(x@pairwise.LD$LD_matrices[[facets]])))){
      msg <- c(msg, "Requested sample subfacet not all located in possible subfacets. Sample subfacet may not be in the provided SNP metadata.")
    }
  }

  if(length(msg) > 0){
    stop(paste0(msg, collapse = "\n  "))
  }
  #==============subfunctions================
  #function to prepare data.
  prep_hm_dat <- function(x, r = NULL){

    # remove columns and rows with no data
    x <- x[!apply(x, 1, function(y)all(is.na(y))), !apply(x, 2, function(y)all(is.na(y)))]

    # set rownames as first column.
    x <- cbind(position = rownames(x), x)
    x <- data.table::as.data.table(x)


    #melting the df and fixing the names of the columns
    heatmap_x <- data.table::melt(x, id.vars = "position")
    names(heatmap_x) <- c("SNPa", "SNPb", "value")

    #getting rid of all the zeros from snps being compared to themselves
    heatmap_x <- na.omit(heatmap_x)
    heatmap_x$SNPa <- as.numeric(as.character(heatmap_x$SNPa))/1000000
    heatmap_x$SNPb <- as.numeric(as.character(heatmap_x$SNPb))/1000000

    #make sure that for all comparisons, SNPa is less than SNPb.
    vio <- which(heatmap_x$SNPa > heatmap_x$SNPb)
    if(length(vio) > 0){
      viofix <- cbind(SNPa = heatmap_x[vio,]$SNPb, SNPb = heatmap_x[vio,]$SNPa, value = heatmap_x[vio,]$value)
      heatmap_x <- rbind(heatmap_x[-vio,], viofix)
    }

    #subset down to the desired r if requested
    if(!is.null(r)){
      heatmap_x <- heatmap_x[heatmap_x$SNPa >= r[1] & heatmap_x$SNPa <= r[2] &
                               heatmap_x$SNPb >= r[1] & heatmap_x$SNPb <= r[2],]
    }

    heatmap_x$value <- as.numeric(heatmap_x$value)
    return(heatmap_x)
  }

  # convert positions to ordered factors, get 10 unique levels from each snp.subfacet to plot!
  order_levels <- function(x){
    # convert positions to factors:
    ms <- unique(c(x$SNPa, x$SNPb))
    ms <- sort(as.numeric(as.character(ms)))

    # finish messing with site names
    ms <- ms
    ms <- sort(ms)
    ms <- as.factor(ms)

    # convert to factor
    x$SNPa <- as.factor(x$SNPa)
    x$SNPb <- as.factor(x$SNPb)

    # reordering based on factors
    x[["SNPa"]]<-factor(x[["SNPa"]],levels=ms, ordered = T)
    x[["SNPb"]]<-factor(x[["SNPb"]],levels=ms, ordered = T)

    # get 10 unique levels per snp.subfacet
    u.subfacets <- unique(x$snp.subfacet)
    comb.levs <- data.frame(pos = numeric(), source = character())
    if(length(u.subfacets) == 0){
      t.levs <- unique(as.numeric(c(levels(x$SNPa), levels(x$SNPb))))

      if(length(t.levs) >= 10){
        t.levs <- t.levs[seq(1, length(t.levs), length.out = 10)]
      }
      comb.levs <- rbind.data.frame(comb.levs, cbind.data.frame(t.levs, ".base", stringsAsFactors = F), stringsAsFactors = F)
    }
    for(i in 1:length(u.subfacets)){
      t.subfacet <- x[x$snp.subfacet == u.subfacets[i],]
      t.subfacet <- droplevels.data.frame(t.subfacet)
      t.levs <- unique(as.numeric(c(levels(t.subfacet$SNPa), levels(t.subfacet$SNPb))))
      if(length(t.levs) >= 10){
        t.levs <- t.levs[seq(1, length(t.levs), length.out = 10)]
      }

      comb.levs <- rbind.data.frame(comb.levs, cbind.data.frame(t.levs, u.subfacets[i], stringsAsFactors = F), stringsAsFactors = F)
    }

    return(list(dat = x, levs = comb.levs))
  }


  #==============collapse all of the LD data into a list, then bind that list together into a single data frame===========
  # for sample level facets:
  if(facet.type == "sample"){
    if(!is.null(sample.subfacet[1])){
      LD_mat_list <- vector("list", length = length(sample.subfacet))
      names(LD_mat_list) <- names(x@pairwise.LD$LD_matrices[[facets]][sample.subfacet])
    }
    else{
      LD_mat_list <- vector("list", length = length(x@pairwise.LD$LD_matrices[[facets]]))
      names(LD_mat_list) <- names(x@pairwise.LD$LD_matrices[[facets]])
    }

    # bind to a list
    tracker <- 1
    for(i in 1:length(x@pairwise.LD$LD_matrices[[facets]])){
      # if all missing data...
      if(all(is.null(x@pairwise.LD$LD_matrices[[facets]][[i]][[LD_measure]]))){
        next()
      }
      if(!is.null(sample.subfacet[1])){
        if(!(names(x@pairwise.LD$LD_matrices[[facets]][i]) %in% sample.subfacet)){
          next
        }
      }
      LD_mat_list[[tracker]] <- cbind(var = names(x@pairwise.LD$LD_matrices[[facets]][i]), prep_hm_dat(x@pairwise.LD$LD_matrices[[facets]][[i]][[LD_measure]], r))
      tracker <- tracker + 1
    }

    # bind the data.tables together
    LD_mats <- order_levels(data.table::rbindlist(LD_mat_list))
    rm(LD_mat_list)
  }
  else if(facet.type == "complex"){

    # intialize either with or without snp.subfacet filtering, and fix if doing sample subfacets
    if(!is.null(snp.subfacet[1])){
      LD_mat_list <- vector("list", length = length(x@pairwise.LD$LD_matrices[[facets]]))
      names(LD_mat_list) <- names(x@pairwise.LD$LD_matrices[[facets]])
      if(!is.null(sample.subfacet))?{
        LD_mat_list <- LD_mat_list[[which(names(LD_mat_list) %in% sample.subfacet)]]
      }
    }
    else{
      if(!is.null(sample.subfacet[1])){
        LD_mat_list <- vector("list",
                              length = length(x@pairwise.LD$LD_matrices[[facets]][sample.subfacet])*length(x@pairwise.LD$LD_matrices[[facets]][[1]]))
      }
      else{
        LD_mat_list <- vector("list",
                              length = length(x@pairwise.LD$LD_matrices[[facets]])*length(x@pairwise.LD$LD_matrices[[facets]][[1]]))
      }
    }


    tracker <- 1
    # add matrices to list
    for(i in 1:length(x@pairwise.LD$LD_matrices[[facets]])){

      # skip if this sample subfacet isn't supposed to be run.
      if(!is.null(sample.subfacet[1])){
        if(!(names(x@pairwise.LD$LD_matrices[[facets]][i])) %in% sample.subfacet){
          next
        }
      }

      rtd <- x@pairwise.LD$LD_matrices[[facets]][[i]] # data for this iteration

      # filter by snp.subfacet if needed
      if(!is.null(snp.subfacet[1])){
        rtd <- rtd[which(names(rtd) %in% snp.subfacet)]
      }

      # add data for each remaining snp.subfacet
      for(j in 1:length(rtd)){
        # if all data is missing...
        if(all(is.na(rtd[[j]][[LD_measure]]))){
          next()
        }
        LD_mat_list[[tracker]] <- cbind(var = names(x@pairwise.LD$LD_matrices[[facets]])[i], snp.subfacet = names(rtd)[j], prep_hm_dat(rtd[[j]][[LD_measure]], r))
        tracker <- tracker + 1
      }
    }
    # bind together
    LD_mats <- order_levels(data.table::rbindlist(LD_mat_list))
  }
  else if(facet.type == "snp"){
    rtd <- x@pairwise.LD$LD_matrices[[facets]][[".base"]]

    # filter by snp.subfacet if needed
    if(!is.null(snp.subfacet[1])){
      rtd <- rtd[which(names(rtd) %in% snp.subfacet)]
    }

    # initialize
    LD_mat_list <- vector("list", length = length(rtd))

    # add data for each remaining snp.subfacet
    for(i in 1:length(rtd)){
      # if all data is missing...
      if(all(is.na(rtd[[i]][[LD_measure]]))){
        next()
      }
      LD_mat_list[[i]] <- cbind(var = names(x@pairwise.LD$LD_matrices[[facets]]), snp.subfacet = names(rtd)[i], prep_hm_dat(rtd[[i]][[LD_measure]], r))
    }

    # bind
    LD_mats <- order_levels(data.table::rbindlist(LD_mat_list))
  }
  else{
    LD_mats <- cbind(var = ".base", snp.subfacet = ".base", prep_hm_dat(x@pairwise.LD$LD_matrices[[".base"]][[".base"]][[LD_measure]], r))

    LD_mats <- order_levels(LD_mats)
  }

  if(exists("LD_mat_list")){rm(LD_mat_list)}
  if(exists("rtd")){rm(rtd)}

  levs <- LD_mats$levs
  LD_mats <- as.data.frame(LD_mats$dat)

  out <- ggplot2::ggplot(LD_mats, ggplot2::aes(x = SNPa, y=SNPb, fill=value, color = value)) +
    ggplot2::geom_tile(color = "white")

  if(length(unique(LD_mats$snp.subfacet)) > 1 & length(unique(LD_mats$var)) > 1){
    out <- out + ggplot2::facet_wrap(var~snp.subfacet, scales = "free")
  }
  else if(length(unique(LD_mats$snp.subfacet)) > 1){
    out <- out + ggplot2::facet_wrap(~snp.subfacet, scales = "free")
  }
  else if(length(unique(LD_mats$var)) > 1){
    out <- out + ggplot2::facet_wrap(~var)
  }

  comb.plot.levels <- levs$t.levs

  out <- out +
    ggplot2::scale_fill_viridis_c(direction = -1, option = viridis.option) +
    ggplot2::scale_color_viridis_c(direction = -1, option = viridis.option) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "",y="", fill=l.text) +
    ggplot2::theme(legend.title= ggplot2::element_text(size = t.sizes[2]),
                   axis.text = ggplot2::element_text(size = t.sizes[5]),
                   panel.grid.major = ggplot2::element_line(color = background),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(hjust = 0.01, size = t.sizes[1]),
                   axis.title = ggplot2::element_text(size = t.sizes[4]),
                   legend.text = ggplot2::element_text(size = t.sizes[3]),
                   panel.background = ggplot2::element_rect(fill = background, colour = background)) +
    ggplot2::scale_x_discrete(breaks = comb.plot.levels, label = abbreviate) +
    ggplot2::scale_y_discrete(breaks = comb.plot.levels, label = abbreviate) +
    ggplot2::ylab("Position (Mb)") + ggplot2::xlab("Position (Mb)")


  if(!is.null(title)){
    out <- out + ggplot2::ggtitle(title)
  }

  out <- list(plot = out, dat = LD_mats)
  return(out)
}

#'PCA and tSNE plots from snpRdata.
#'
#'Generate a ggplot cluster plot based on PCA or the he Barnes-Hut simulation at
#'theta>0 implemented in \code{\link[Rtsne]{Rtsne}}. Works by conversion to the
#'"sn" format described in \code{\link{format_snps}} with interpolated missing
#'genotypes.
#'
#'Cluster plots can be produced via either PCA or tSNE. The PCA point
#'coordinates are calculated using \code{\link{prcomp}}. By default, the first
#'two principal coordinates are plotted. A PC matrix will also be returned for
#'easy plotting of other PCs. tSNE coordinates are calculated via
#'\code{\link[Rtsne]{Rtsne}}, which should be consulted to for more details
#'about this method. Stated simply, tSNE attempts to compress a
#'multi-dimensional PCA (PCs 1:n) into fewer dimensions while retaining as much
#'information as possible. As such, a tSNE plot can be seen as a representation
#'of many different PC axis compressed into a single two-dimensional plot. This
#'compression process is stochastic, and so plots will vary somewhat between
#'runs, and multiple runs are recommended.
#'
#'For more details on tSNE aruments, \code{\link[Rtsne]{Rtsne}} should be
#'consulted.
#'
#'Data points for individuals can be automatically colored by any sample-level
#'facet categories. Facets should be provided as described in
#'\code{\link{Facets_in_snpR}}. Up to two different sample-level facets can be
#'automatically plotted simultaniously.
#'
#'@param x snpRdata object.
#'@param facets character, default NULL. Categorical sample-level metadata
#'  variables by which to color points. Up to two different sample-specific
#'  facets may be provided. See \code{\link{Facets_in_snpR}} for more details.
#'@param plot.type character, default c("PCA", "tSNE"). Types of plots to be
#'  produced, see description.
#'@param check.duplicates logical, default FALSE. Checks for any duplicated
#'  individuals, which will cause errors. Since these rarely exist and
#'  drastically slow down function run-time, this defaults to FALSE.
#'@param minimum_percent_coverage numeric, default FALSE. Proportion of samples
#'  a SNP must be sequenced in to be used in generating plots.
#'@param interpolation_method character, default "bernoulli". Interpolation
#'  method to use for missing data. Options: \itemize{\item{bernoulli:
#'  }{Interpolated via binomial draw for each allele against minor allele
#'  frequency.} \item{af: }{Interpolated by inserting the expected number of
#'  minor alleles at missing data points given loci minor allele frequencies.}}
#'@param dims numeric, default 2. Output dimensionality, default 2.
#'@param initial_dims numeric, default 50. The number of dimensions retained in
#'  the initial PCA step.
#'@param perplexity numeric, default FALSE. Perplexity parameter, by default
#'  found by \code{\link[mmtsne]{hbeta}}, with beta = 1.
#'@param theta numeric, default 0. Theta parameter from
#'  \code{\link[Rtsne]{Rtsne}}. Default an exhaustive search.
#'@param iter numeric, default 5000. Number of tSNE iterations to perform.
#'@param ... Other arguments, passed to \code{\link[Rtsne]{Rtsne}}.
#'
#'@return A list containing: \itemize{ \item{data: } Raw PCA and/or tSNE plot
#'  data. \item{plots: } ggplot PCA and/or tSNE plots. } Each of these two list
#'  may contain one or two objects, one for each PCA or tSNE plot requested,
#'  named "pca" and "tsne", respectively.
#'
#'@author William Hemstrom
#'@author Matt Thorstensen
#'
#'@references Jesse H. Krijthe (2015). Rtsne: T-Distributed Stochastic Neighbor Embedding using a Barnes-Hut Implementation, URL: \url{https://github.com/jkrijthe/Rtsne}
#'@seealso \code{\link[mmtsne]{mmtsne}}
#'
#'@export
#'
#' @examples
#' # plot colored by population
#' plot_clusters(stickSNPs, "pop")
plot_clusters <- function(x, facets = FALSE, plot_type = c("PCA", "tSNE"), check_duplicates = FALSE,
                          minimum_percent_coverage = FALSE, interpolation_method = "bernoulli",
                          dims = 2, initial_dims = 50, perplexity = FALSE, theta = 0, iter = 1000, ...){

  #=============sanity checks============
  msg <- character(0)

  # mpc
  if(minimum_percent_coverage != FALSE){
    problem <- F
    if(!is.numeric(minimum_percent_coverage)){
      problem <- T
    }
    if(length(minimum_percent_coverage) > 1){
      problem <- T
    }
    if(minimum_percent_coverage > 1 | minimum_percent_coverage <= 0){
      problem <- T
    }
    if(problem){
      msg <- c(msg, "minimum_percent_coverage must be a numeric value between 0 and 1.")
    }
  }

  # facets
  facets <- check.snpR.facet.request(x, facets, remove.type = "none", return.type = T)
  if(length(facets[[1]]) > 1){
    msg <- c(msg, "Only one facet may be specified at a time. This facet may be complex, with up to two sample levels (e.g. pop.family).")
  }
  if(any(facets[[2]] != "sample")){
    bf <- facets[[1]][which(facets[[2]] != "sample")]
    msg <- c(msg,  paste0("Only sample level facets can be plotted. Facets ", paste0(bf, collapse = ", "), " refer to snp metadata."))
  }
  facets <- facets[[1]]
  facets <- unlist(strsplit(facets, "(?<!^)\\.", perl = T))
  facets <- check.snpR.facet.request(x, facets)

  plot_type <- tolower(plot_type)
  good_plot_types <- c("pca", "tsne")
  if(any(!(plot_type %in% good_plot_types))){
    msg <- c(msg, paste0("Unaccepted plot_type. Accepted types:", paste0(good_plot_types, collapse = ", "), "."))
  }
  if(!(length(plot_type) %in% length(good_plot_types))){
    msg <- c(msg, paste0("No more than", length(good_plot_types), "plots may be made at once."))
  }


  if(length(msg) > 0){
    stop(paste0(msg, collapse = "  \t"))
  }

  #=============prepare dataset===============
  cat("Formatting data...\n")

  # check for matching sn plot:
  if(length(x@sn) != 0){
    if(x@sn$type != interpolation_method){
      suppressWarnings(x@sn$sn <- format_snps(x, "sn", interpolate = interpolation_method))
      x@sn$type <- interpolation_method
    }
  }
  else{
    suppressWarnings(x@sn$sn <- format_snps(x, "sn", interpolate = interpolation_method))
    x@sn$type <- interpolation_method
  }

  sn <- x@sn$sn
  sn <- sn[,-c(1:(ncol(x@snp.meta) - 1))]
  meta <- x@sample.meta

  # figure out which snps should be dropped due to low coverage
  if(minimum_percent_coverage != FALSE){
    counts <- rowSums(x@geno.tables$gs[x@facet.meta$subfacet == ".base",])
    counts <- counts/ncol(x)
    bad.loci <- which(counts <= minimum_percent_coverage)
    if(length(bad.loci) > 0){
      sn <- sn[-bad.loci,]
    }
  }

  if(check_duplicates){
    cat("Checking for duplicates...\n")
    dups <- which(duplicated(sn) | duplicated(sn, fromLast=TRUE))
    if(length(dups) > 0){
      cat("Duplicates detected, indices:", dups, "\nRemoving all of these!\n")
      sn <- sn[,-dups]
    }
  }


  rm.snps <- nrow(sn)
  if(rm.snps == 0){
    stop("No remaining SNPs after filtering.\n")
  }
  if(rm.snps < 20){
    warning("Few remaining SNPs after filtering! Remaining SNPs:", rm.snps, ".\n")
  }

  cat("Plotting using", rm.snps, "loci.\n")
  sn <- t(sn)

  #=============prepare plot data=============
  plot_dats <- vector("list", length(plot_type))
  names(plot_dats) <- plot_type

  if("pca" %in% plot_type){
    cat("Preparing pca...\n")
    pca_r <- prcomp(as.matrix(sn))
    pca <- as.data.frame(pca_r$x) #grab the PCA vectors.
    plot_dats$pca <- cbind(meta, pca)  #add metadata that is present in the input.
  }
  if("tsne" %in% plot_type){
    #get perplexity if not provided
    if(perplexity == FALSE){
      cat("Estimating perplexity...\n")
      perp <- mmtsne::hbeta(sn, beta = 1)
      perplexity <- perp$H
    }

    #run the tSNE
    cat("Running tSNE...\n")
    tsne.out <- Rtsne::Rtsne(X = sn, dims = dims, initial_dims = initial_dims, perplexity = perplexity,
                             theta = theta, max_iter = iter, check_duplicates = FALSE,
                             verbose=TRUE, ...)
    colnames(tsne.out$Y) <- paste0("PC", 1:ncol(tsne.out$Y))
    plot_dats$tsne <- cbind(meta, as.data.frame(tsne.out$Y))
  }
  #=============make ggplots=====================
  plots <- vector("list", length(plot_dats))
  names(plots) <- names(plot_dats)

  if(length(facets) > 2){
    warning("Only up to two simultanious sample metadata columns can be plotted at once.\n")
  }

  for(i in 1:length(plot_dats)){
    tpd <- plot_dats[[i]]

    #Categories (pops, fathers, mothers, ect.) are given in plot.vars argument. Supports up to two!
    #make the base plot, then add categories as color and fill.
    out <- ggplot2::ggplot(tpd, ggplot2::aes(PC1, PC2)) + ggplot2::theme_bw() #initialize plot


    if(facets[1] == FALSE){
      plots[[i]] <- plot = out + ggplot2::geom_point()
      next
    }


    #add variables.
    v1 <- tpd[,which(colnames(tpd) == facets[1])] #get the factors

    # add geoms to plot
    if(length(facets) == 1){
      out <- out + ggplot2::geom_point(ggplot2::aes(color = v1))#add the factor
    }
    else{
      # if two plotting variables, prepare the second and add it as well.
      v2 <- tpd[,which(colnames(tpd) == facets[2])]
      out <- out + ggplot2::geom_point(ggplot2::aes(color = v1, fill = v2), pch = 21, size = 2.5, stroke = 1.25)
    }


    # change the color scales for the variables
    ## for the first variable
    # if the data is likely continuous:
    if(is.numeric(v1)){
      out <- out + ggplot2::scale_color_viridis_c(name = facets[1])
    }
    else{
      out <- out + ggplot2::scale_color_viridis_d(name = facets[1])
    }

    ## for the second variable if defined
    if(length(facets) > 1){
      if(is.numeric(v2)){
        out <- out + ggplot2::scale_fill_viridis_c(name = facets[2])
      }
      else{
        out <- out + ggplot2::scale_fill_viridis_d(name = facets[2])
      }
    }

    if(plot_type[i] == "pca"){
      loadings <- (pca_r$sdev^2)/sum(pca_r$sdev^2)
      loadings <- round(loadings * 100, 2)
      out <- out + ggplot2::xlab(paste0("PC1 (", loadings[1], "%)")) + ggplot2::ylab(paste0("PC2 (", loadings[2], "%)"))
    }
    else{
      out <- out + ggplot2::xlab("Dim 1") + ggplot2::ylab("Dim 2")
    }
    print(out)
    plots[[i]] <- out
  }
  return(list(data = plot_dats, plots = plots))
}

#' Generate a manhattan plot from snpRdata or a data.frame.
#'
#' Creates a ggplot-based manhattan plot, where chromosomes/scaffolds/ect are
#' concatenated along the x-axis. Can optionally highlight requested SNPs or
#' those that pass an arbitrary significance threshold and facet plots by
#' defined sample-specific variables such as population.
#'
#' Unlike most snpR functions, this function works with either a snpRdata object
#' or a data.frame. For snpRdata objects snp-specific or sliding window
#' statistics can be plotted. In both cases, the facet argument can be used to
#' define facets to plot, as described in \code{\link{Facets_in_snpR}}. For
#' typical stats, name of the snp meta-data column containing
#' chromosome/scaffold information must be supplied to the "chr" argument. For
#' windowed stats, chr is instead inferred from the snp-specific facet used to
#' create the smoothed windows. In both cases, the requested facets must exactly
#' match those used to calculate statistics! If x is a data frame, the "chr"
#' argument must also be given, and the "facets" argument will be ignored.
#'
#' A column defining the position of the SNP within the chromsome must be
#' provided, and is "position" by default.
#'
#' Specific snp and chr levels can also be requested using the chr.subfacet and
#' sample.subfacet arguments. See examples. For data.frames, sample.subfacets
#' levels must refer to a column in x titled "subfacet".
#'
#' Specific snps can be highlighted and annotated. If a significance level is
#' requested, SNPs above this level will be highlighted by default. SNPs above
#' the suggestive line can also be highlighted by providing "suggestive" to the
#' highlight argument. Alternatively, individual SNPs can be highlighted by
#' providing a numeric vector. For snpR data, this will correspond to the SNP's
#' row in the snpRdata object. For data.frames, it will correspond to a
#' ".snp.id" column if it exists, and the row number if not. The label for
#' highlighted SNPs will be either chr_bp by default or given in the column
#' named by the "snp" argument.
#'
#' @param x snpRdata or data.frame object containing the data to be plotted.
#' @param plot_var character. A character string naming the statistic to be
#'   plotted. For snpRdata, these names correspond to any previously calculated
#'   statistics.
#' @param window logical, default FALSE. If TRUE, sliding window averages will
#'   instead be plotted. These averages must have first been calculated with
#'   calc_smoothed_averags. Ignored if x is a data.frame.
#' @param facets character or NULL, default NULL. Facets by which to break
#'   plots, as described in \code{\link{Facets_in_snpR}}. For non-window stats,
#'   the any snp.specific facets will be ignored. Ignored if x is a data.frame.
#' @param chr character, default "chr". Column in either snp metadata or x (for
#'   snpRdata or data.frame objects, respectively) which defines the
#'   "chromosome" by which SNP positions will be concatenated along the x-axis.
#'   If window = TRUE and a snpRdata object, this will be ignored in favor of
#'   the SNP specific facet provided to the facets argument.
#' @param bp character, default "bp". Column in either snp metadata or x (for
#'   snpRdata or data.frame objects, respectively) which defines the position in
#'   bp of each SNP.
#' @param snp character, default NULL. Column in either snp metadata or x (for
#'   snpRdata or data.frame objects, respectively) containing snpIDs to use for
#'   highlighting. Ignored if no highlighting is requested.
#' @param chr.subfacet character, default NULL. Specific chromosomes to plot.
#'   See examples.
#' @param sample.subfacet character, default NULL. Specific sample-specific
#'   levels of the provided facet to plot. If x is a data.frame, this can refer
#'   to levels of a column titled "subfacet". See examples.
#' @param significant numeric, default NULL. Value at which a line will be drawn
#'   designating significant SNPs. If highlight = "significant", SNPs above this
#'   level will also be labeled.
#' @param suggestive numeric, default NULL. Value at which a line will be drawn
#'   designating suggestive SNPs. If highlight = "suggestive", SNPs above this
#'   level will also be labeled.
#' @param pval logical, default FALSE. If TRUE, treats values lower than the
#'   significance threshold as significant.
#' @param log.p logical, default FALSE. If TRUE, plot variables and thresholds
#'   will be transformed to -log.
#' @param highlight character, numeric, or FALSE, default "significant".
#'   Controls SNP highlighting. If either "significant" or "suggestive", SNPs
#'   above those respetive values will be highlighted. If a numeric vector, SNPs
#'   corresponding to vector entries will be highlighted. See details.
#' @param viridis.option character, default "plasma". Viridis color scale option
#'   to use for significance lines and SNP labels. See
#'   \code{\link[ggplot2]{scale_colour_gradient}} for details.
#' @param viridis.hue numeric, default c(0.2, 0.5). Two values between 0 and 1
#'   listing the hues at which to start and stop on the viridis palette defined
#'   by the viridis.option argument. Lower numbers are darker.
#' @param t.sizes numeric, default c(16, 12, 10). Text sizes, given as
#'   c(strip.title, axis, axis.ticks).
#' @param colors character, default c("black", "slategray3"). Colors to
#'   alternate across chromosomes.
#'
#' @author William Hemstrom
#' @export
#'
#' @return A list containing \itemize{\item{plot: } A ggplot manhattan plot.
#'   \item{data: } Raw plot data.}
#'
#'
#' @examples
#' # make some data
#' x <- calc_basic_snp_stats(x, "pop.group", sigma = 200, step = 50)
#'
#' # plot pi, breaking apart by population, keeping only the groupIX and
#' # groupIV chromosomes and the ASP, PAL, and SMR populations, with
#' # significant and suggestive lines plotted and SNPs
#' # with pi below the significance level labeled.
#' plot_manhattan(x, "pi", facets = "pop",
#' chr = "group", chr.subfacet = c("groupIX", "groupIV"),
#' sample.subfacet = c("ASP", "OPL", "SMR"),
#' significant = 0.05, suggestive = 0.15, pval = T)
#'
#' # plot FST for the ASP/PAL comparison across all chromosomes,
#' # labeling the first 10 SNPs in x (by row) with their ID
#' plot_manhattan(x, "fst", facets = "pop.group",
#' sample.subfacet = "ASP~PAL", highlight = 1:20,
#' chr = "group", snp = ".snp.id")
#'
#' # plot sliding-window FST between ASP and CLF
#' # and between OPL and SMR
#' plot_manhattan(x, "fst", window = T, facets = c("pop.group"),
#' chr = "group", sample.subfacet = c("ASP~CLF", "OPL~SMR"),
#' significant = .29, suggestive = .2)
#'
#' # plot using a data.frame,
#' # using log-transformed p-values
#' ## grab data
#' y <- get.snpR.stats(x, "pop")
#' ## plot
#' plot_manhattan(y, "pHWE", facets = "pop", chr = "group",
#' significant = 0.0001, suggestive = 0.001,
#' log.p = T, highlight = F)
#'
plot_manhattan <- function(x, plot_var, window = FALSE, facets = NULL,
                           chr = "chr", bp = "position", snp = NULL,
                           chr.subfacet = NULL, sample.subfacet = NULL,
                           significant = NULL, suggestive = NULL,
                           highlight = "significant",
                           pval = FALSE, log.p = FALSE,
                           viridis.option = "plasma", viridis.hue = c(.2, 0.5), t.sizes = c(16, 12, 10),
                           colors = c("black", "slategray3")){

  #=============grab the desired stats=====================
  #====if a snpRdata object========
  if(class(x) == "snpRdata"){
    if(window == FALSE){
      facets <- check.snpR.facet.request(x, facets)
    }
    else{
      facets <- check.snpR.facet.request(x, facets, "none")
    }
    if(plot_var %in% colnames(x@stats)){
      if(window){
        if(chr != "chr"){
          warning("chr variable will be set to the snp level facet provided to the facets argument for sliding windows.\n")
        }
        stats <- get.snpR.stats(x, facets = facets, type = "single.window")
        chr <- "snp.subfacet"

      }
      else{
        stats <- get.snpR.stats(x, facets = facets)
      }
    }
    else if(plot_var %in% colnames(x@pairwise.stats)){
      if(window){
        if(chr != "chr"){
          warning("chr variable will be set to the snp level facet provided to the facets argument for sliding windows.\n")
        }
        stats <- get.snpR.stats(x, facets, "pairwise.window")
        chr <- "snp.subfacet"
      }
      else{
        stats <- get.snpR.stats(x, facets, "pairwise")
      }
    }

    if(nrow(stats) == 0){
      stop("No matching statistics.\n")
    }
  }

  #====otherwise=====
  else if(is.data.frame(x)){
    stats <- x
  }
  else{
    stop("x must be a data.frame or snpRdata object.\n")
  }

  #====clean up=====
  # fix comparison column name
  if(any(colnames(stats) == "comparison")){
    colnames(stats)[which(colnames(stats) == "comparison")] <- "subfacet"
  }


  # remove unwanted snp or sample subfacets
  if(!is.null(chr.subfacet)){
    stats <- stats[which(stats[,chr] %in% chr.subfacet),]
  }
  if(!is.null(sample.subfacet)){
    stats <- stats[which(stats$subfacet %in% sample.subfacet),]
  }
  if(nrow(stats) == 0){
    stop("No matching statistics.\n")
  }

  #=============figure out adjusted positions on x axis================
  # fetch chromosome info
  if(is.factor(stats[,chr])){
    stats[,chr] <- as.character(stats[,chr])
  }

  # get adjusted x axis positions and figure out where chromsome tick marks should be placed.
  chr.info <- tapply(stats[,bp], stats[,chr], max) # chromosome lengths
  chr.centers <- (chr.info - tapply(stats[,bp], stats[,chr], min))/2 # chromosome centers
  chr.info <- cumsum(chr.info) # cumulative lengths
  cum.bp <- c(0, chr.info[-length(chr.info)]) # correct cummulative lengths
  names(cum.bp) <- names(chr.info) # rename
  cum.chr.centers <- cum.bp + chr.centers # cumulative chromosome centers
  stats$start <- cum.bp[match(stats[,chr], names(cum.bp))]
  stats$cum.bp <- stats$start + stats[,bp]

  #=============clean up==================
  colnames(stats)[which(colnames(stats) == plot_var)] <- "pvar"
  colnames(stats)[which(colnames(stats) == chr)] <- "chr"
  if(log.p){
    stats$pvar <- -log(stats$pvar)
    if(!is.null(significant)){
      significant <- -log(significant)
    }
    if(!is.null(suggestive)){
      suggestive <- -log(suggestive)
    }
  }

  #=============snps to highlight=====================
  if(highlight[1] != FALSE & !is.null(highlight[1]) & !is.na(highlight[1])){
    # unless the defaults are set...
    if(!((highlight == "significant") & is.null(significant))){
      do.highlight <- T

      # using significant
      if(highlight == "significant"){
        if(pval){
          stats$highlight <- ifelse(stats$pvar <= significant, 1, 0)
        }
        else{
          stats$highlight <- ifelse(stats$pvar >= significant, 1, 0)
        }
      }
      #using suggestive
      else if(highlight == "suggestive"){
        if(pval){
          stats$highlight <- ifelse(stats$pvar <= suggestive, 1, 0)
        }
        else{
          stats$highlight <- ifelse(stats$pvar >= suggestive, 1, 0)
        }
      }

      # using a list of snps to highlight
      else{
        stats$highlight <- 0
        if(".snp.id" %in% colnames(stats)){
          stats$highlight[which(stats$.snp.id %in% highlight)] <- 1
        }
        stats$highlight[highlight] <- 1
      }

      # labels
      if(is.null(snp)){
        stats$highlight.label <- paste0(stats$chr, "_", stats[,bp])
      }
      else{
        stats$highlight.label <- stats[,snp]
      }
    }
    else{
      do.highlight <- F
    }
  }
  else{
    do.highlight <- F
  }

  #============produce the plot========
  p <- ggplot2::ggplot(stats, ggplot2::aes(x = cum.bp, y = pvar, color = chr)) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 90, size = t.sizes[3], vjust = 0.5),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = t.sizes[3]),
                   strip.text = ggplot2::element_text(hjust = 0.01, size = t.sizes[1]),
                   axis.title = ggplot2::element_text(size = t.sizes[2])) +
    ggplot2::scale_x_continuous(label = names(cum.chr.centers), breaks = cum.chr.centers, minor_breaks = NULL) +
    ggplot2::scale_color_manual(values = rep(c(colors), length(cum.chr.centers))) +
    ggplot2::xlab(chr) + ggplot2::ylab(plot_var)

  #=============adjust the plot========
  if(length(unique(stats$subfacet) > 1)){
    p <- p + ggplot2::facet_wrap(~subfacet)
  }

  add.palette <- viridis::viridis(3, option = viridis.option, begin = viridis.hue[1], end = viridis.hue[2])

  # significant/suggestive lines
  if(!is.null(significant)){
    p <- p + ggplot2::geom_hline(yintercept = significant, color = add.palette[1])
  }
  if(!is.null(suggestive)){
    p <- p + ggplot2::geom_hline(yintercept = suggestive, color = add.palette[2])
  }

  # highlight
  if(do.highlight){
    p <- p + ggrepel::geom_label_repel(data = stats[which(stats$highlight == 1),],
                                       mapping = ggplot2::aes(label = highlight.label), color = add.palette[3],
                                       force = 1.3)
  }

  print(p)

  return(list(plot = p, data = stats))
}
