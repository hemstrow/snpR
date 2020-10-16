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
plot_pairwise_LD_heatmap <- function(x, facets = NULL, snp.subfacet = NULL, sample.subfacet = NULL, LD_measure = "CLD", r = NULL,
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
  good.ld.measures <- c("rsq", "Dprime", "pval", "CLD")
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
    ms <- as.factor(sort(as.numeric(as.character(ms))))

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
    LD_mats <- cbind(var = ".base", snp.subfacet = ".base", prep_hm_dat(x@pairwise.LD$LD_matrices[[".base"]][[".base"]][[".base"]][[LD_measure]], r))

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

#' PCA, tSNE, and umap plots from snpRdata.
#'
#' Generate a ggplot cluster plot based on PCA, the Barnes-Hut simulation at
#' theta>0 implemented in \code{\link[Rtsne]{Rtsne}}, or the Uniform Manifold
#' Approximation and Projection approach implemented in \code{\link[umap]{umap}}.
#' Works by conversion to the "sn" format described in \code{\link{format_snps}}
#' with interpolated missing genotypes.
#'
#' Cluster plots can be produced via, PCA, tSNE, or umap. The PCA point
#' coordinates are calculated using \code{\link{prcomp}}. By default, the first
#' two principal coordinates are plotted. A PC matrix will also be returned for
#' easy plotting of other PCs. tSNE coordinates are calculated via
#' \code{\link[Rtsne]{Rtsne}}, which should be consulted to for more details
#' about this method. Stated simply, tSNE attempts to compress a
#' multi-dimensional PCA (PCs 1:n) into fewer dimensions while retaining as much
#' information as possible. As such, a tSNE plot can be seen as a representation
#' of many different PC axis compressed into a single two-dimensional plot. This
#' compression process is stochastic, and so plots will vary somewhat between
#'  runs, and multiple runs are recommended. Uniform Manifold Approximation and Projection (UMAP)
#' coordinates are calculated via \code{\link[umap]{umap}}. UMAP similarly attempts to reduce multi-dimensional
#' results to a two dimensional visualization.
#'
#' For more details on tSNE aruments, \code{\link[Rtsne]{Rtsne}} should be
#' consulted.
#'
#' Additional arguments to the UMAP can be also be provided. Additional information on these
#' arguments can be found in \code{\link[umap]{umap.defaults}}.
#'
#' Data points for individuals can be automatically colored by any sample-level
#' facet categories. Facets should be provided as described in
#' \code{\link{Facets_in_snpR}}. Up to two different sample-level facets can be
#' automatically plotted simultaniously.
#'
#' @param x snpRdata object.
#' @param facets character, default NULL. Categorical sample-level metadata
#'   variables by which to color points. Up to two different sample-specific
#'   facets may be provided. See \code{\link{Facets_in_snpR}} for more details.
#' @param plot.type character, default c("PCA", "tSNE", "umap"). Types of plots
#'   to be produced, see description.
#' @param check.duplicates logical, default FALSE. Checks for any duplicated
#'   individuals, which will cause errors. Since these rarely exist and
#'   drastically slow down function run-time, this defaults to FALSE.
#' @param minimum_percent_coverage numeric, default FALSE. Proportion of samples
#'   a SNP must be sequenced in to be used in generating plots.
#' @param minimum_genotype_percentage numeric, default FALSE. Proportion of SNPs
#'   a sample must be sequenced at in order to be used in plots.
#' @param interpolation_method character, default "bernoulli". Interpolation
#'   method to use for missing data. Options: \itemize{\item{bernoulli:
#'   }{Interpolated via binomial draw for each allele against minor allele
#'   frequency.} \item{af: }{Interpolated by inserting the expected number of
#'   minor alleles at missing data points given loci minor allele frequencies.}}
#' @param dims numeric, default 2. Output dimensionality, default 2.
#' @param initial_dims numeric, default 50. The number of dimensions retained in
#'   the initial PCA step during tSNE.
#' @param perplexity numeric, default FALSE. Perplexity parameter, by default
#'   found by \code{\link[mmtsne]{hbeta}}, with beta = 1.
#' @param theta numeric, default 0. Theta parameter from
#'   \code{\link[Rtsne]{Rtsne}}. Default an exhaustive search.
#' @param iter numeric, default 1000. Number of tSNE iterations/umap epochs to
#'   perform.
#' @param viridis.option character, default "viridis". Viridis color scale option
#'   to use for significance lines and SNP labels. See
#'   \code{\link[ggplot2]{scale_colour_gradient}} for details.
#' @param alt.palette charcter or NULL, default NULL. Optional palette of colors
#'   to use instead of the viridis palette.
#' @param ncp numeric or NULL, default NULL. Number of components to consider for iPCA sn format
#'   interpolations of missing data. If null, the optimum number will be estimated, with the
#'   maximum specified by ncp.max. This can be very slow.
#' @param ncp.max numeric, default 5. Maximum number of components to check for when determining
#'   the optimum number of components to use when interpolating sn data using the iPCA approach.
#' @param ... Other arguments, passed to \code{\link[Rtsne]{Rtsne}} or \code{\link[umap]{umap}}.
#'
#' @return A list containing: \itemize{ \item{data: } Raw PCA, tSNE, and umap plot
#'  data. \item{plots: } ggplot PCA, tSNE, and/or umap plots.} Each of these two lists
#'  may contain one, two, or three objects, one for each PCA, tSNE, or umap plot requested,
#'  named "pca" and "tsne", and "umap", respectively.
#'
#' @author William Hemstrom
#' @author Matt Thorstensen
#'
#' @references Jesse H. Krijthe (2015). Rtsne: T-Distributed Stochastic Neighbor Embedding using a Barnes-Hut Implementation, URL: \url{https://github.com/jkrijthe/Rtsne}.
#' @references Van Der Maaten, L. & Hinton, G. (2008) Visualizing high-dimensional data using t-SNE. journal of machine learning research. \emph{Journal of Machine Learning Research}.
#' @references McInnes, L. & Healy (2018). UMAP: uniform manifold approximation and projection. Preprint at URL: \url{https://arxiv.org/abs/1802.03426}.
#'
#' @seealso \code{\link[mmtsne]{mmtsne}}
#'
#' @export
#'
#' @examples
#' # plot colored by population
#' plot_clusters(stickSNPs, "pop")
plot_clusters <- function(x, facets = FALSE, plot_type = c("PCA", "tSNE", "umap"), check_duplicates = FALSE,
                          minimum_percent_coverage = FALSE, minimum_genotype_percentage = FALSE, interpolation_method = "bernoulli",
                          dims = 2, initial_dims = 50, perplexity = FALSE, theta = 0, iter = 1000,
                          viridis.option = "viridis", alt.palette = NULL, ncp = NULL, ncp.max = 5, ...){

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
  good_plot_types <- c("pca", "tsne", "umap")
  if(any(!(plot_type %in% good_plot_types))){
    msg <- c(msg, paste0("Unaccepted plot_type. Accepted types:", paste0(good_plot_types, collapse = ", "), "."))
  }

  if("umap" %in% plot_type){
    pkg.check <- check.installed("umap")
    if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
  }
  
  if("tsne" %in% plot_type){
    pkg.check <- check.installed("Rtsne")
    if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
    
    pkg.check <- check.installed("mmtsne")
    if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
  }
  

  if(length(msg) > 0){
    stop(paste0(msg, collapse = "  \t"))
  }

  #=============prepare dataset===============
  cat("Formatting data...\n")

  # check for matching sn plot:
  if(length(x@sn) != 0){
    if(x@sn$type != interpolation_method){
      suppressWarnings(x@sn$sn <- format_snps(x, "sn", interpolate = interpolation_method, ncp = ncp, ncp.max = ncp.max))
      x@sn$type <- interpolation_method
    }
  }
  else{
    suppressWarnings(x@sn$sn <- format_snps(x, "sn", interpolate = interpolation_method, ncp = ncp, ncp.max = ncp.max))
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

  # figure out which samples should be dropped due to poor genotyping
  if(minimum_genotype_percentage != FALSE){
    counts <- colSums(ifelse(x == x@mDat, F, T))
    counts <- counts/nrow(x)
    bad.samples <- which(counts <= minimum_genotype_percentage)
    if(length(bad.samples) > 0){
      sn <- sn[,-bad.samples]
      meta <- meta[-bad.samples,]
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
  rm.ind <- ncol(sn)
  if(rm.snps == 0){
    stop("No remaining SNPs after filtering.\n")
  }
  if(rm.ind == 0){
    stop("No remaining samples after filtering.\n")
  }
  if(rm.snps < 20){
    warning("Few remaining SNPs after filtering! Remaining SNPs:", rm.snps, ".\n")
  }

  cat("Plotting using", rm.snps, "loci and", rm.ind, "samples.\n")
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
  if("umap" %in% plot_type){
    cat("Preparing umap...\n")
    umap_r <- umap::umap(as.matrix(sn), n_epochs = iter, verbose = T, ...)
    colnames(umap_r$layout) <- paste0("PC", 1:ncol(umap_r$layout))
    plot_dats$umap <- cbind(meta, as.data.frame(umap_r$layout))
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
      if(is.null(alt.palette)){
        out <- out + ggplot2::scale_color_viridis_c(name = facets[1], option = viridis.option)
      }
      else{
        out <- out + ggplot2::scale_color_gradient(low = alt.palette[1], high = alt.palette[2], name = facets[1])
      }
    }
    else{
      if(is.null(alt.palette)){
        out <- out + ggplot2::scale_color_viridis_d(name = facets[1], option = viridis.option)
      }
      else{
        out <- out + ggplot2::scale_color_manual(values = alt.palette, name = facets[1])
      }
    }

    ## for the second variable if defined
    if(length(facets) > 1){
      if(is.numeric(v2)){
        if(is.null(alt.palette)){
          out <- out + ggplot2::scale_fill_viridis_c(name = facets[2], option = viridis.option)
        }
        else{
          out <- out + ggplot2::scale_fill_gradient(low = alt.palette[1], high = alt.palette[2], name = facets[2])
        }
      }
      else{
        if(is.null(alt.palette)){
          out <- out + ggplot2::scale_fill_viridis_d(name = facets[2], option = viridis.option)
        }
        else{
          out <- out + ggplot2::scale_fill_manual(values = alt.palette, name = facets[2])
        }
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
#' @param sig_below logical, default FALSE. If TRUE, treats values lower than the
#'   significance threshold as significant.
#' @param log.p logical, default FALSE. If TRUE, plot variables and thresholds
#'   will be transformed to -log.
#' @param abs logical, default FALSE. If TRUE, converts the plot variable to it's absolute value.
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
#' significant = 0.05, suggestive = 0.15, sig_below = T)
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
                           sig_below = FALSE, log.p = FALSE, abs = FALSE,
                           viridis.option = "plasma", viridis.hue = c(.2, 0.5), t.sizes = c(16, 12, 10),
                           colors = c("black", "slategray3")){

  #=============sanity checks==============================
  msg <- character()
  pkg.check <- check.installed("ggrepel")
  if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
  
  if(length(msg) > 0){
    stop(msg, collapse = "\n")
  }
  #=============grab the desired stats=====================
  #====if a snpRdata object========
  if(class(x) == "snpRdata"){
    if(window == FALSE){
      facets <- check.snpR.facet.request(x, facets)
    }
    else{
      if(chr != "chr"){
        warning("chr variable will be set to the snp level facet provided to the facets argument for sliding windows.\n")
      }
      if(!is.null(facets)){
        pop.facets <- check.snpR.facet.request(x, facets, "snp")
        facets <- paste0(pop.facets, ".", chr)
      }
      else{
        facets <- chr
      }

      facets <- check.snpR.facet.request(x, facets, "none")
    }
    if(plot_var %in% colnames(x@stats)){
      if(window){

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
    stats$pvar <- -log10(stats$pvar)
    if(!is.null(significant)){
      significant <- -log10(significant)
    }
    if(!is.null(suggestive)){
      suggestive <- -log10(suggestive)
    }
  }

  if(abs){
    stats$pvar <- abs(stats$pvar)
  }

  #=============snps to highlight=====================
  if(highlight[1] != FALSE & !is.null(highlight[1]) & !is.na(highlight[1])){
    # unless the defaults are set...
    if(!((highlight[1] == "significant") & is.null(significant))){
      do.highlight <- T

      # using significant
      if(highlight[1] == "significant"){
        if(sig_below){
          stats$highlight <- ifelse(stats$pvar <= significant, 1, 0)
        }
        else{
          stats$highlight <- ifelse(stats$pvar >= significant, 1, 0)
        }
      }
      #using suggestive
      else if(highlight[1] == "suggestive"){
        if(sig_below){
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
  if(length(unique(as.character(stats$subfacet)) > 1)){
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

  return(list(plot = p, data = stats))
}

#' Create STRUCTURE-like cluster plots
#'
#' Creates ggplot-based stacked barcharts of assignment probabilities (Q) into
#' an arbitrary 'k' number of clusters like those produced by the program
#' STRUCTURE. Runs for each value of k between 2 and the given number.
#'
#' Individual cluster assignment probabilities can be calculated using several
#' different methods: \itemize{\item{snmf: } sNMF (sparse Non-Negative Matrix
#' Factorization). \item{snapclust: } Maximum-likelihood genetic clustering.}
#' These methods are not re-implemented in R, instead, this function calls the
#' \code{\link[LEA]{snmf}} and \code{\link[adegenet]{snapclust.choose.k}}
#' functions instead. Please cite the references noted in those functions if
#' using this function. For snapclust, the "ward" method is used to initialize clusters
#' if one rep is requested, otherwise the clusters are started randomly each rep. Other
#' methods can be used by providing pop.ini as an additional agument as long as only one
#' rep is requested.
#'
#' Multiple different runs can be conducted using the 'reps' argument, and the
#' results can be combined for plotting across all of these reps using the
#' clumpp option. This option calls the CLUMPP software package in order to
#' combine proportion population membership across multiple runs via
#' \code{\link[pophelper]{clumppExport}}. Again, please cite both CLUMPP and
#' pophelper if using this option.
#'
#' Since CLUMPP is run independantly for each value of K, cluster identites
#' often "flip" between K values. For example individuals that are grouped into
#' cluster 1 and K = 3 may be grouped into cluster 2 at K = 4. To adjust this,
#' cluster IDs are iteratively adjusted across K values by flipping IDs such
#' that the euclidian distances between clusters at K and K - 1 are minimized.
#' This tends to produce consistant cluster IDs across multiple runs of K.
#'
#' Individuals can be sorted into by membership proportion into different
#' clusters within populations using the qsort option.
#'
#' Since the clustering and CLUMPP processes can be time consuming and outside
#' tools (such as NGSadmix or fastSTRUCTURE) may be prefered, a nested list of Q
#' matrices, sorted by K and then rep or a character string giving a pattern
#' matching saved Q matrix files in the current working directory may provided
#' directly instead of a snpRdata object. Note that the output for this
#' funciton, if run on a snpRdata object, will return a properly formatted list
#' of Q files (named 'data') in addition to the plot and plot data. This allows
#' for the plot to be quickly re-constructed using different sorting parameters
#' or facets. In these cases, the facet argument should instead be a vector of
#' group identifications per individuals.
#'
#' Note that several files will be created in the working directory when using
#' this function that are not automatically cleaned after use.
#'
#' @param x snpRdata object, list of Q matrices (sorted by K in the first level
#'   and run in the second), or a character string designating a pattern
#'   matching Q matrix files in the current working directories.
#' @param facet character, default NULL. If provided, individuals will not be
#'   noted on the x axis. Instead, the levels of the facet will be noted. Only a
#'   single, simple, sample specific facet may be provided. Individuals must be
#'   sorted by this facet in x. If Q matrices are provided (either directly or
#'   via file path), this should instead be a vector of group identities for
#'   each individual (populations, etc.).
#' @param facet.order character, default NULL. Optional order in which the
#'   levels of the provided facet should appear on the plot, left to right.
#' @param k numeric, default 2. The maximum of k (number of clusters) for which
#' to run the clustering/assignment algorithm. The values 2:k will be run.
#' @param method character, default "snmf". The clustering/assignment method to
#'   run. Options: \itemize{\item{snmf: } sNMF (sparse Non-Negative Matrix
#'   Factorization). \item{snapclust: } Maximum-likelihood genetic clustering.}
#'   See \code{\link[LEA]{snmf}} or \code{\link[adegenet]{snapclust.choose.k}}
#'   for details, respectively.
#' @param reps numeric, default 1. The number of independent clustering
#'   repititions to run.
#' @param iterations numeric or Inf, default 1000. For snapclust, the maximum
#'   number of iterations to run.
#' @param I numeric or NULL, default NULL. For snmf, how many SNPs should be
#'   used to initialize the search? Initializing with a subset of the total SNPs
#'   can radically speed up computation time for large datasets.
#' @param alpha numeric, default 10. For sNMF, determines the regularization
#'   parameter. For small datasets, this can have a large effect, and should
#'   probably be larger than the default. See documentation for
#'   \code{\link[LEA]{snmf}}.
#' @param qsort character, numeric, or FALSE, default "last". Determines if
#'   individuals should be sorted (possibly within facet levels) by cluster
#'   assignment proportion. If not FALSE, determines which cluster to use for
#'   sorting (1:k). If "last" or "first" sorts by those clusters.
#' @param qsort_K numeric or character, default "last". If qsorting is
#'   performed, determines the reference k value by which individuals are
#'   sorted. If "first" or "last", sorts by k = 2 or k = k, respectively.
#' @param clumpp logical, default T. Specifies if CUMPP should be run to
#'   collapse results across multiple reps. If FALSE, will use only the first
#'   rep for plotting.
#' @param clumpp.opt character, default "greedy". Designates the CLUMPP method
#'   to use. Options: \itemize{ \item{fullsearch: } Search all possible
#'   configurations. Slow. \item{greedy: } The standard approach. Slow for large
#'   datasets at high k values. \item{large.k.greedy: } A fast but less accurate
#'   approach. } See CLUMPP documentation for details.
#' @param ID character or NULL, default NULL. Designates a column in the sample
#'   metadata containing sample IDs.
#' @param viridis.option character, default "viridis". Viridis color scale
#'   option. See \code{\link[ggplot2]{scale_colour_gradient}} for details.
#' @param alt.palette charcter or NULL, default NULL. Optional palette of colors
#'   to use instead of the viridis palette.
#' @param t.sizes numeric, default c(12, 12, 12). Text sizes, given as
#'   c(strip.title, axis, axis.ticks).
#' @param ... additional arguments passed to either \code{\link[LEA]{snmf}} or
#'   \code{\link[adegenet]{snapclust.choose.k}}.
#'
#' @export
#' @author William Hemstrom
#' @references Frichot E, Mathieu F, Trouillon T, Bouchard G, Francois O.
#'   (2014). Fast and Efficient Estimation of Individual Ancestry Coefficients.
#'   \emph{Genetics}, 194(4): 973–983.
#' @references Frichot, Eric, and Olivier François (2015). LEA: an R package for
#'   landscape and ecological association studies. \emph{Methods in Ecology and
#'   Evolution}, 6(8): 925-929.
#' @references Beugin, M. P., Gayet, T., Pontier, D., Devillard, S., & Jombart,
#'   T. (2018). A fast likelihood solution to the genetic clustering problem.
#'   \emph{Methods in ecology and evolution}, 9(4), 1006-1016.
#' @references Francis, R. M. (2017). pophelper: an R package and web app to
#'   analyse and visualize population structure. \emph{Molecular ecology
#'   resources}, 17(1), 27-32.
#' @references Jakobsson, M., & Rosenberg, N. A. (2007). CLUMPP: a cluster
#'   matching and permutation program for dealing with label switching and
#'   multimodality in analysis of population structure. \emph{Bioinformatics},
#'   23(14), 1801-1806.
#'
#' @return A list containing: \itemize{\item{plot: } A ggplot object.
#'   \item{data: } A nested list of the raw Q matrices, organized by K and then
#'   by run. \item{plot_data: } The raw data used in constructing the ggplot.
#'   \item{K_plot: } A data.frame containing the value suggested for use in K
#'   selection vs K value for the selected method.}
#'
plot_structure <- function(x, facet = NULL, facet.order = NULL, k = 2, method = "snmf", reps = 1, iterations = 1000,
                           I = NULL, alpha = 10, qsort = "last", qsort_K = "last", clumpp = T,
                           clumpp.opt = "greedy", ID = NULL, viridis.option = "viridis",
                           alt.palette = NULL, t.sizes = c(12, 12, 12), ...){

  #===========sanity checks===================
  msg <- character()
  provided_qlist <- FALSE

  # check if this is with a snpRdata object or a qlist and do some other checks
  if(!class(x) == "snpRdata"){
    # file pattern
    if(!is.list(x)){
      if(is.character(x) & length(x) > 1){
        msg <- c(msg, "Unaccepted input format. x must be a snpRdata object, a list of q matrices, or a string containing a pattern that matches qfiles in the current working directory.\n")
      }
      else if(!is.character(x)){
        msg <- c(msg, "Unaccepted input format. x must be a snpRdata object, a list of q matrices, or a string containing a pattern that matches qfiles in the current working directory.\n")
      }
      else{
        provided_qlist <- "parse"
      }
    }
    # qlist, typical k then r format.
    else{
      provided_qlist <- TRUE
      bad.list <- F

      # check that list depth is 2
      l1 <- lapply(x, is.list) # depth two?
      if(!all(unlist(l1))){
        bad.list <- T
      }
      else{
        l2 <- logical(length(x)) # depth three?
        for(i in 1:length(x)){
          if(any(unlist(lapply(x[[i]], function(y) is.matrix(y) | is.data.frame(y))) == F)){
            l2[i] <- T
          }

          # also check that this is k and then r (all qlists should have the same number of columns at the second level.)
          else if(length(unique(lapply(x[[i]], ncol))) != 1){
            bad.list <- T
          }
        }
        if(any(l2)){
          bad.list <- T
        }
      }
      if(bad.list){
        msg <- c(msg, "Provided qlist must be a two level list, sorted by k values then by run.\n")
      }
    }

    if(!is.null(facet[1])){
      if(!is.character(facet)){
        msg <- c(msg, "For a provided qlist, facet must be a character vector containing sample metadata.\n")
      }
      if(provided_qlist == TRUE){
        if(bad.list == F){
          if(length(facet) != nrow(qlist[[1]])){
            msg <- c(msg, "The number of samples in the qlist does not match the length of the provided sample metadata.\n")
          }
        }
      }
      sample_meta <- data.frame(d = facet, stringsAsFactors = F)
      facet <- deparse(substitute(facet))
      colnames(sample_meta) <- facet
    }
  }


  if(method == "snmf"){
    if(!("LEA" %in% rownames(installed.packages()))){
      msg <- c(msg, "The package LEA must be installed to run sNMF assignment. LEA can be installed via bioconductor with BiocManager::install('LEA')\n")
    }
  }

  if(provided_qlist == FALSE & reps == "all"){
    msg <- c(msg, "reps = 'all' uninterpretable if a qlist is not provided.\n")
  }
  if(clumpp & reps == 1){
    clumpp <- FALSE
    warning("Since only one rep is requested, clumpp will not be run.\n")
  }

  if(clumpp){
    if(!("pophelper" %in% rownames(installed.packages()))){
      msg <- c(msg, "The package pophelper must be installed to run clumpp. pophelper can be installed via the devtools or remotes packages with devtools::install_github('royfrancis/pophelper') or remotes::install_github('royfrancis/pophelper')\n")
    }
    good.clumpp.opts <- c("fullsearch", "greedy", "large.k.greedy")
    clumpp.opt <- tolower(clumpp.opt)
    if(!clumpp.opt %in% good.clumpp.opts){
      msg <- c(msg, paste0("Unaccepted clumpp option. Accepted options: ", paste0(good.clumpp.opts, collapse = ", "), "\n"))
    }
  }
  else if(provided_qlist == "parse"){
    if(!("pophelper" %in% rownames(installed.packages()))){
      msg <- c(msg, "The package pophelper must be installed to run clumpp. pophelper can be installed via the devtools or remotes packages with devtools::install_github('royfrancis/pophelper') or remotes::install_github('royfrancis/pophelper')\n")
    }
  }

  # checks for snpRdata objects only
  if(provided_qlist == FALSE){
    good.methods <- c("snapclust", "snmf")
    if(!method %in% good.methods){
      msg <- c(msg, paste0("Unaccepted clustering method. Accepted options: ", paste0(good.methods, collapse = ", "), "\n"))
    }

    if(length(facet) > 1){
      msg <- c(msg, "Only one facet may be plotted at once.\n")
    }
    if(!is.null(facet[[1]])){
      fcheck <- check.snpR.facet.request(x, facet, remove.type = "none", return.type = T)
      if(any(fcheck[[2]] != "sample")){
        msg <- c(msg, paste0("Only simple, sample level facets allowed.\n"))
      }
    }
    sample_meta <- x@sample.meta
  }

  # palette checks
  if(!is.null(alt.palette[1])){

    # is it long enough?
    if(length(alt.palette) < k){
      msg <- c(msg, "Provided alternative palette must contain at least as many colors
               as the maximum k value.\n")
    }

    # is everything a valid color (can ggplot work with it)?
    else{
      alt.palette <- alt.palette[1:k]
      tpd <- matrix(rnorm(k*2, 0, 1), k, 2)
      tpd <- cbind(as.data.frame(tpd), col = 1:k)

      color.check <- ggplot2::ggplot(tpd, ggplot2::aes(V1, V2, color = as.factor(col))) +
        ggplot2::geom_point() + ggplot2::scale_color_manual(values = alt.palette)
      tempfile <- tempfile()
      tempfile <- paste0(tempfile, ".pdf")

      suppressMessages(res <- try(ggplot2::ggsave(tempfile, color.check), silent = T))
      invisible(capture.output(file.remove(tempfile)))
      if("try-error" %in% class(res)){
        msg <- c(msg, "Alt.palette contains non-color entries or otherwise fails to properly plot.\n")
      }

      rm(res, tpd, color.check)
    }
  }

  if(is.null(facet[1])){qsort <- F}

  if(length(msg) != 0){
    stop(msg)
  }

  #===========sub-functions===================
  # fix things so that the cluster ID is the same in all sets
  fix_clust <- function(x){

    #loop through each q object
    for (i in 2:length(x)){
      #see which columns in the previous run are the most similar to each column

      #initialize mapping df
      mdf <- data.frame(tcol = 1:ncol(x[[i]]), pcol = numeric(ncol(x[[i]])),
                        ed = numeric(ncol(x[[i]])))

      #loop through each column and find where to map it.
      for (j in 1:ncol(x[[i]])){

        #intialize euc distance vector
        elist <- numeric(ncol(x[[i - 1]]))

        #compare to each other col.
        for(k in 1:ncol(x[[i-1]])){
          #save euclidian dist
          elist[k] <- sum((x[[i]][,j] - x[[i-1]][,k])^2, na.rm = T)
        }

        #save results
        mdf[j,2] <- which.min(elist)
        mdf[j,3] <- min(elist)
      }

      #reassign clusters in this qdf
      ##which is the new cluster? Probably that with the most distance to any original clusters.
      dups <- duplicated(mdf[,2]) | duplicated(mdf[,2], fromLast = T)
      nc <- which.max(mdf[dups,3])
      mdf[dups,2][nc] <- nrow(mdf)
      mdf <- mdf[order(mdf[,2]),]

      ##reasign clusters
      tdf <- x[[i]]
      tdf <- tdf[,mdf[,1]]

      ##replace object in x with the re-arranged qfile.
      colnames(tdf) <- colnames(x[[i]])
      x[[i]] <- tdf
    }

    return(x)
  }

  # sort the individuals within each population based on the qvals
  Q_sort <- function(x, pop, cluster = "first", q = "last"){


    #get which pop to use
    if(q == "last"){
      q <- length(x)
    }

    #get order to stick individual in
    lx <- x[[q]]
    upops <- unique(pop)
    lx$s <- 1:nrow(lx)

    #get the sorting cluster priority:
    if(cluster == "first"){
      cseq <- (ncol(lx)-1):1
    }
    else if (cluster == "last"){
      cseq <- 1:(ncol(lx)-1)
    }
    else if (is.numeric(cluster)){
      if(length(cluster) == ncol(lx)){
        cseq <- cluster
      }
      else{
        if(length(cluster) < ncol(lx)){
          cseq <- c((1:ncol(lx))[-which(1:ncol(lx) %in% cluster)], rev(cluster))
        }
        else{
          stop("Cluster length is longer than number of clusters in x element q.\n")
        }
      }
    }

    for(i in 1:nrow(upops)){
      tx <- lx[which(pop == upops[i,]),]
      for(j in cseq){
        tx <- tx[order(tx[,j]),]
      }
      lx[which(pop == upops[i,]),] <- tx
    }

    # return the order
    return(lx$s)
  }

  # process 1 q result by adjusting column names and melting. Takes a single qlist dataframe/matrix
  process_1_q <- function(tq, x, sample_meta = NULL){
    colnames(tq) <- 1:ncol(tq)
    tk <- ncol(tq)
    if(is.null(sample_meta)){
      sample_meta <- paste0("s", 1:nrow(tq))
    }
    tq <- cbind(tq, sample_meta)

    if(is.null(ID)){
      if(is.null(colnames(x))){
        tq$ID <- paste0("samp", 1:nrow(tq))
      }
      else{
        tq$ID <- colnames(x)
      }
    }
    else{
      tq$ID <- sample_meta[,ID]
    }

    tq <- reshape2::melt(tq, id.vars = colnames(tq)[-c(1:tk)])
    colnames(tq)[((ncol(tq) - 1):ncol(tq))] <- c("Cluster", "Percentage")
    tq$Cluster <- as.numeric(tq$Cluster)
    tq$K <- tk

    return(tq)
  }

  # write a list of q files to a directory with informative names for k and run
  prep_clumpp <- function(x){
    # make a directory and write the k files to it
    for(i in 1:length(x)){
      for(j in 1:reps){
        write.table(x[[i]][[j]], paste0("K", ncol(x[[i]][[j]]), "r", j, "qopt"), sep = " ", quote = F, col.names = F, row.names = F)
      }
    }
  }

  # run clumpp on a directory of q files using pophelper. Import the results into a list of  processed q tables
  run_clumpp <- function(pattern = "qopt"){
    # prepare files and run clumpp
    qfiles <- list.files(full.names = T, pattern = pattern)
    qlist <- pophelper::readQ(qfiles)
    if(clumpp.opt == "large.k.greedy"){
      clumpp.opt <- 3
    }
    else if(clumpp.opt == "greedy"){
      clumpp.opt <- 2
    }
    else{
      clumpp.opt <- 1
    }
    ## grab only the correct k values:
    ks <- unlist(lapply(qlist, function(x) attr(x, "k")))
    qlist <- qlist[which(ks <= k)]
    ## run
    pophelper::clumppExport(qlist, useexe = T, parammode = clumpp.opt)
    pophelper::collectClumppOutput(filetype = "both")

    # import results
    mq <- pophelper::readQ(list.files("pop-both/", full.names = T, pattern = "merged"))

    # get only the correct k values (in case this is re-running in a previous directory)
    ks <- unlist(lapply(mq, function(x) attr(x, "k")))
    mq <- mq[which(ks <= k)]

    # fix cluster IDs to match
    mq <- fix_clust(mq)
    save.q <- mq



    # sort by facet
    if(!is.null(facet)){

      if(!is.null(facet.order)){
        mq <- sort_by_pop_qfiles(mq, sample_meta, facet.order)
      }
      else{
        mq <- sort_by_pop_qfiles(mq, sample_meta, unique(sample_meta[,facet]))
      }
      sample_meta <- mq$meta
      mq <- mq$qlist
    }

    # sort by Q values if requested
    if(qsort != FALSE){
      if(!is.null(facet)){
        pop <- as.data.frame(sample_meta[,facet], stringsAsFactors = F)
      }
      else{
        pop <- as.data.frame(rep("pop1", nrow(sample_meta)), stringsAsFactors = F)
      }
      ind_ord <- Q_sort(mq, pop, qsort, qsort_K)
    }
    else{
      ind_ord <- 1:nrow(mq[[1]])
    }

    # process into plottable form
    for(i in 1:length(mq)){
      if(exists("sample_meta")){
        mq[[i]] <- process_1_q(mq[[i]], x, sample_meta)
      }
      else{
        mq[[i]] <- process_1_q(mq[[i]], x)
      }
    }

    return(list(q = mq, ord = ind_ord, qlist = save.q))
  }

  # function to parse in q files, used only if clumpp not run and a file pattern is provided
  parse_qfiles <- function(pattern){
    # read in the files
    qfiles <- list.files(full.names = T, pattern = pattern)
    qlist <- pophelper::readQ(qfiles)

    # get k values
    att <- lapply(qlist, attributes)
    att <- unlist(lapply(att, function(x) x$k))

    # inititalize
    out_q <- vector("list", length = length(unique(att)))
    names(out_q) <- paste0("K_", sort(unique(att)))

    # for each k, load up each run
    for(i in 1:length(unique(att))){

      # find matching runs and initialize
      tk <- sort(unique(att))[i]
      hits <- which(att == tk)
      out_q[[i]] <- vector("list", length(hits))
      names(out_q[[i]]) <- paste0("r_", 1:length(hits))

      # load
      for(j in 1:length(hits)){
        out_q[[i]][[j]] <- qlist[[hits[j]]]
      }
    }

    return(out_q)
  }

  # function to sort each element of a qlist by a facet (population, ect). Needed because they must
  # be ordered properly for qsorting to work correctly.
  sort_by_pop_qfiles <- function(qlist, meta, pop.ord){
    ord <- factor(meta[,facet], levels = pop.ord)
    ord <- order(ord)

    for(i in 1:length(qlist)){
      if(is.data.frame(qlist[[i]])){
        qlist[[i]] <- qlist[[i]][ord,]
      }
      else{
        for(j in 1:length(qlist[[i]])){
          qlist[[i]][[j]] <- qlist[[i]][[j]][ord,]
        }
      }
    }

    nmeta <- as.data.frame(meta[ord,], stringsAsFactors = F)
    colnames(nmeta) <- colnames(meta)

    return(list(qlist = qlist, meta = nmeta))
  }

  #===========run the assignment/clustering method===============
  # each method should return a list of q tables, unprocessed, and possibly work on a K_plot. The list should be nested, k then r.
  if(provided_qlist == FALSE){
    if(method == "snapclust"){
      # format and run snapclust
      invisible(capture.output(x.g <- format_snps(x, "adegenet")))

      # initialize
      qlist <- vector("list", length = k - 1)
      names(qlist) <- paste0("K_", 2:k)
      for(i in 1:length(qlist)){
        qlist[[i]] <- vector("list", length = reps)
        names(qlist[[i]]) <- paste0("r_", 1:reps)
      }
      K_plot <- vector("list", length(rep))

      # run once per rep
      for(i in 1:reps){
        if(reps == 1){
          res <- adegenet::snapclust.choose.k(x = x.g, max = k, IC.only = FALSE, ...)
        }
        else{
          res <- adegenet::snapclust.choose.k(x = x.g, max = k, IC.only = FALSE, pop.ini = NULL, ...)
        }

        for(j in 1:length(res$objects)){
          qlist[[j]][[i]] <- res$objects[[j]]$proba

          # fix instances where NaNs are likely actually 1s, based on https://github.com/thibautjombart/adegenet/issues/221
          nas <- which(is.na(rowSums(qlist[[j]][[i]])))
          if(length(nas) > 0){
            fixable <- qlist[[j]][[i]][nas,]
            fixable <- is.nan(fixable)
            not.fixable1 <- which(rowSums(fixable) != 1) # rows with more than one na
            not.fixable2 <- which(rowSums(qlist[[j]][[i]][nas,] == 0, na.rm = T) != ncol(qlist[[j]][[i]]) - 1) # rows where the other values aren't all zero
            not.fixable <- union(not.fixable1, not.fixable2)
            if(length(not.fixable) > 0){
              if(length(not.fixable) != length(nas)){
                qlist[[j]][[i]][nas,][-not.fixable,][is.nan(qlist[[j]][[i]][nas,][-not.fixable,])] <- 1
              }
            }
            else{
              qlist[[j]][[i]][nas,][fixable] <- 1
            }
          }
        }
        K_plot[[i]] <- data.frame(val = res[[1]], K = 1:k, rep = i, stringsAsFactors = F)
      }

      # fix K plot
      K_plot <- dplyr::bind_rows(K_plot)
      colnames(K_plot)[1] <- names(res)[1]
    }
    else if(method == "snmf"){
      # format data
      if(file.exists("lea_input.geno")){
        file.remove("lea_input.geno")
      }
      invisible(capture.output(format_snps(x, "lea", outfile = "lea_input.geno")))

      # run the analysis
      if(!is.null(I)){
        snmf.res <- LEA::snmf("lea_input.geno", K = 2:k, repetitions = reps, iterations = iterations, entropy = T, project = "new", alpha = alpha, I = I, ...)
      }
      else{
        snmf.res <- LEA::snmf("lea_input.geno", K = 2:k, repetitions = reps, iterations = iterations, entropy = T, project = "new", alpha = alpha, ...)
      }

      # read in
      qlist <- vector("list", k - 1)
      names(qlist) <-  paste0("K_", 2:k)
      K_plot <- vector("list", k - 1)
      for(i in 2:k){

        # process each rep
        qlist[[i - 1]] <- vector("list", reps)
        names(qlist[[i - 1]]) <-  paste0("r_", 1:reps)
        for(j in 1:reps){
          qlist[[i - 1]][[j]] <- as.data.frame(LEA::Q(snmf.res, i, j))
        }
        K_plot[[i - 1]] <- data.frame(Cross.Entropy = LEA::cross.entropy(snmf.res, k), K = i)
      }

      # fix K plot
      K_plot <- dplyr::bind_rows(K_plot)
      colnames(K_plot)[1] <- "Cross.Entropy"
    }
  }
  else if(provided_qlist == TRUE){
    qlist <- x
  }



  #===========run clumpp=========================================
  if(clumpp){

    if(provided_qlist != "parse"){
      if(!dir.exists("qfiles")){
        dir.create("qfiles")
      }
      setwd("qfiles")

      # prepare files
      prep_clumpp(qlist)

      # run clumpp
      pdat <- run_clumpp()
      setwd("..")
    }
    else{
      pdat <- run_clumpp(pattern = x)
    }





    # concatenate results
    ind_ord <- pdat$ord
    cq <- pdat$qlist
    pdat <- dplyr::bind_rows(pdat$q)

    ## add the clumpp qlist
    if(!exists("qlist")){
      qlist <- vector("list", length(cq))
    }
    for(i in 1:length(qlist)){
      qlist[[i]][["clumpp"]] <- cq[[i]]
    }

  }
  else{
    # parse in if provided with a pattern
    if(provided_qlist == "parse"){
      qlist <- parse_qfiles(x)
    }

    if(reps != 1){
      warning("Plotting only the first run.\n")
    }

    # sort by facet
    if(!is.null(facet)){
      if(!is.null(facet.order)){
        qlist <- sort_by_pop_qfiles(qlist, sample_meta, facet.order)
      }
      else{
        qlist <- sort_by_pop_qfiles(qlist, sample_meta, unique(sample_meta[,facet]))
      }
      sample_meta <- qlist$meta
      qlist <- qlist$qlist
    }

    # grab just the first run
    tq <- vector("list", length(qlist))
    for(i in 1:length(tq)){
      tq[[i]] <- qlist[[i]][[1]]
    }

    # process
    ## fix cluster IDs to match
    tq <- fix_clust(tq)

    ## qsort
    if(qsort != FALSE){
      if(!is.null(facet)){
        pop <- as.data.frame(sample_meta[,facet], stringsAsFactors = F)
      }
      else{
        pop <- as.data.frame(rep("pop1", nrow(sample_meta)), stringsAsFactors = F)
      }
      ind_ord <- Q_sort(tq, pop, qsort, qsort_K)
    }
    else{
      ind_ord <- 1:nrow(tq[[1]])
    }

    ## process into plottable form
    for(i in 1:length(tq)){
      if(exists("sample_meta")){
        tq[[i]] <- process_1_q(tq[[i]], x, sample_meta)
        
      }
      else{
        tq[[i]] <- process_1_q(tq[[i]], x)
      }
    }

    ## concatenate
    pdat <- dplyr::bind_rows(tq)
  }

  #===========final fixes to plot data=========
  pdat$ID <- factor(pdat$ID, levels = pdat$ID[ind_ord])
  pdat$K <- as.factor(pdat$K)
  levels(pdat$K) <- paste0("K = ", levels(pdat$K))
  pdat$Cluster <- as.factor(pdat$Cluster)

  #===========make plot==============
  p <- ggplot2::ggplot(pdat, ggplot2::aes(ID, Percentage, color = Cluster, fill = Cluster)) +
    ggplot2::facet_wrap(~K, ncol = 1, strip.position = "right") +
    ggplot2::theme_bw() +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(0.25, 0.5,0.75)) +
    ggplot2::ylab("Cluster Membership Proportion") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = t.sizes[3]),
                   strip.text = ggplot2::element_text(size = t.sizes[1]),
                   axis.title = ggplot2::element_text(size = t.sizes[2]),
                   strip.background = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   panel.spacing = ggplot2::unit(0.1, "lines"))

  # add palette
  if(!is.null(alt.palette[1])){
    p <- p + ggplot2::scale_color_manual(values = alt.palette) +
      ggplot2::scale_fill_manual(values = alt.palette)
  }
  else{
    p <- p + ggplot2::scale_color_viridis_d(option = viridis.option) +
      ggplot2::scale_fill_viridis_d(option = viridis.option)
  }

  # add facets
  if(!is.null(facet[1])){
    pops <- unique(sample_meta[,facet])
    fmt <- sapply(pops, function(x) sum(sample_meta[,facet] == x, na.rm = T))
    names(fmt) <- pops
    fmc <- cumsum(fmt)
    fm <- floor(c(0, fmc[-length(fmc)]) + fmt/2)
    breaks <- levels(pdat$ID)[fm]

    seps <- c(0, fmc) + 0.5
    seps[1] <- -.5
    p <- p +
      ggplot2::scale_x_discrete(labels = unique(pdat[,facet]), breaks = breaks, expand = c(0,0)) +
      ggplot2::geom_vline(xintercept = c(fmc[-length(fmc)]) + 0.5, color = "white", size = 1) +
      ggplot2::xlab(label = facet[1])
  }
  else{
    p <- p + ggplot2::scale_x_discrete(expand = c(0,0))
  }


  if(exists("K_plot")){
    return(list(plot = p, data = qlist, plot_data = pdat, K_plot = K_plot))
  }
  else{
    return(list(plot = p, data = qlist, plot_data = pdat))
  }


}



#' Plot 1 or 2d site frequency spectra.
#'
#' Plot 1 or 2d site frequency spectra such as those created by
#' \code{\link{make_SFS}}. Plots are made using ggplot2, and can be freely
#' modified as is usual for ggplot objects.
#'
#' The input SFS is either a 2d site frequency spectra stored in a matrix, with
#' an additional "pops" attribute containing population IDs, such as c("POP1",
#' "POP2"), where the first pop is the matrix columns and the second is the
#' matrix rows, or a 1d site frequency spectra stored as a numeric vector with a
#' similar pops attribute giving the population name. These objects can be
#' produced from a dadi input file using \code{\link{make_SFS}}.
#'
#' @param sfs matrix or numeric. Either a 2d site frequency spectra stored in a
#'   matrix, with an additional "pops" attribute containing population IDs, such
#'   as c("POP1", "POP2"), where the first pop is the matrix columns and the
#'   second is the matrix rows, or a 1d site frequency spectra stored as a
#'   numeric vector with a similar pops attribute giving the population name.
#'   These objects can be produced from a dadi input file using
#'   \code{\link{make_SFS}}.
#' @param viridis.option character, default "inferno". Viridis color scale
#'   option. See \code{\link[ggplot2]{scale_colour_gradient}} for details.
#' @param log logical, default TRUE. If TRUE, the number of SNPs in each SFS
#'   cell is log transformed.
#'
#' @return A ggplot2 plot object of the provided SFS.
#'
#' @export
plot_sfs <- function(sfs, viridis.option = "inferno", log = TRUE){

  # add colun names, row names, and melt
  pops <- attr(sfs, "pop")
  sfs <- as.data.frame(sfs)
  if(length(pops) == 1){
    sfs <- as.data.frame(t(sfs))
  }
  colnames(sfs) <- 0:(ncol(sfs) - 1)
  sfs$count <- 0:(nrow(sfs) - 1)
  sfs[1,1] <- NA # mask the first entry
  msfs <- reshape2::melt(sfs, id.vars = "count")
  colnames(msfs) <- c("p1", "p2", "N")
  msfs$p1 <- as.integer(as.character(msfs$p1))
  msfs$p2 <- as.integer(as.character(msfs$p2))



  # plot
  ## 2D
  if(nrow(sfs) > 1){
    if(log){
      p <- ggplot2::ggplot(msfs, ggplot2::aes(x = p1, y = p2, fill = log10(N), color = log10(N)))
    }
    else{
      p <- ggplot2::ggplot(msfs, ggplot2::aes(x = p1, y = p2, fill = N, color = N))
    }

    p <- p +
      ggplot2::geom_tile() + ggplot2::theme_bw() +
      ggplot2::scale_color_viridis_c(na.value = "white", option = viridis.option) +
      ggplot2::scale_fill_viridis_c(na.value = "white", option = viridis.option) +
      ggplot2::xlab(pops[1]) + ggplot2::ylab(pops[2]) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0))
  }

  ## 1D
  else{
    if(log){
      p <- ggplot2::ggplot(msfs, ggplot2::aes(x = p2, y = log10(N)))
    }
    else{
      p <- ggplot2::ggplot(msfs, ggplot2::aes(x = p2, y= N))
    }
    p <- p +
      ggplot2::geom_line() + ggplot2::theme_bw() +
      ggplot2::xlab(pops[1])
    ggplot2::scale_x_continuous(expand = c(0, 0))
  }


  return(p)
}

#' @export
plot_structure_map <- function(assignments, K, facet, pop_coordinates, map = ggplot2::map_data("world2"),
                                pop_names = T, viridis.option = "viridis", alt.palette = NULL, radius_scale = 0.05, label_scale = .75, point_padding_scale = .25){
  #===================sanity checks=================
  msg <- character()
  pkg.check <- check.installed("scatterpie")
  if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
  
  if(length(msg) > 0){stop(msg, "\n")}
  
  #==================plot===========================
  # generate plotting data.frame
  pie_dat <- as.data.frame(matrix(0, nrow = length(unique(assignments$plot_data[,1])), ncol = 3 + K))
  colnames(pie_dat) <- c("pop", "lat", "long", paste0("Cluster ", 1:K))
  tpd <- assignments$plot_data[assignments$plot_data$K == paste0("K = ", K),]
  tpd <- tpd[,c(facet, "Cluster", "Percentage")]
  tpd$Cluster <- as.numeric(tpd$Cluster)
  anc <- tapply(tpd$Percentage, tpd[,c(facet, "Cluster")], mean)

  ## check if we need to flip negative longitude scores
  if(all(map$long >= 0)){
    flip <- T
  }
  else{
    flip <- F
  }

  for(i in 1:ncol(pop_coordinates)){
    pie_dat[i,1] <- colnames(pop_coordinates)[i]
    pie_dat[i,-1] <- c(pop_coordinates[1,i], pop_coordinates[2,i], anc[colnames(pop_coordinates)[i],])
    if(pie_dat$long[i] < 0 & flip){
      pie_dat$long[i] <- 360 + pie_dat$long[i]
    }
  }

  # figure out the radius to use
  lat_range <- range(pie_dat$lat)
  long_range <- range(pie_dat$long)
  r <- min(lat_range[2] - lat_range[1], long_range[2] - long_range[1])
  r <- r*radius_scale

  # make the plot
  mp <- ggplot2::ggplot(map, ggplot2::aes(x = long, y = lat)) +
    ggplot2::geom_map(map = map, ggplot2::aes(map_id = region), fill = "grey", color = "white") +
    ggplot2::theme_bw() +
    ggplot2::xlim(c(min(pie_dat$long - r), max(pie_dat$long) + r)) +
    ggplot2::ylim(c(min(pie_dat$lat - r), max(pie_dat$lat) + r)) +
    scatterpie::geom_scatterpie(data = pie_dat, mapping = ggplot2::aes(x = long, y = lat, r = r), cols = colnames(pie_dat)[4:ncol(pie_dat)]) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")

  if(!is.null(alt.palette)){
    mp <- mp + ggplot2::scale_fill_manual(values = alt.palette)
  }
  else{
    mp <- mp + ggplot2::scale_fill_viridis_d(option = viridis.option)
  }
  if(pop_names){
    # add labels
    #mp + ggplot2::geom_label(data = pie_dat, mapping = ggplot2::aes(x = long, y = lat, label = pop), size = r*label_scale)
    mp <- mp + ggrepel::geom_label_repel(data = pie_dat, mapping = ggplot2::aes(x = long, y = lat, label = pop), size = r*label_scale, point.padding = r*point_padding_scale, max.iter = 10000)
  }

  return(mp)
}
