#'Create a heatmap from pairwise linkage data.
#'
#'Prepares a ggplot2 heatmap from pairwise linkage disequilibrium data stored in
#'a snpRdata object.
#'
#'Since the output is a ggplot object, options can be added or changed by adding
#'"+ function()" to the end of \code{LD_pairwise_heatmap}. Some common options
#'are also built into this function as arguments, but can be overwritten freely.
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
#'@param l.text character, default "CLD". Legend title.
#'@param viridis.option character, default "inferno". Viridis color scale option
#'  to use. Other color scales may be substituted by appending the
#'  scale_color_continuous and scale_fill_continuous ggplot functions to the
#'  produced plot using the '+' operator. See
#'  \code{\link[ggplot2]{scale_gradient}} for details.
#'@param gradient_colors character vector of length 2 or NULL, default NULL. If
#'  provided, the low and high LD colors for the color scale. Provided as an
#'  alternative to the viridis scales. Take care to set a \code{background}
#'  color not on this scale. For example, \code{c("white", "black")} is an
#'  excellent choice, but will make missing data inaccurately appear  appear as
#'  if it has a very low LD unless a different \code{background} color is
#'  selected.
#'@param title character. Plot title.
#'@param t.sizes numeric, default c(16, 13, 10, 12, 10). Text sizes, given as
#'  c(title, legend.title, legend.ticks, axis, axis.ticks).
#'@param background character, default "white". Background color for plot.
#'
#'@return A pairwise LD heatmap as a ggplot object.
#'
#'@author William Hemstrom
#'@author Nicholas Sard
#'@export
#'
#' @examples
#' \dontrun{
#' # get LD data
#' dat <- calc_pairwise_ld(stickSNPs, c("pop.chr"))
#'
#' # produce plots for linkage group IX in the ASP and CLF populations.
#' plot_pairwise_ld_heatmap(dat, c("pop.chr"), "groupIX", c("ASP", "CLF"))
#'
#' # produce plots for every population for linkage group IV
#' plot_pairwise_ld_heatmap(dat, c("pop.chr"), "groupIV")
#' }
#'
plot_pairwise_ld_heatmap <- function(x, facets = NULL, snp.subfacet = NULL, sample.subfacet = NULL, LD_measure = "CLD", r = NULL,
                                     l.text = "CLD", viridis.option = "inferno",
                                     gradient_colors = NULL,
                                     title = NULL, t.sizes = c(16, 13, 10, 12, 10),
                                     background = "white"){
  #==============sanity checks===========
  msg <- character()

  if(length(x@pairwise.LD) == 0){
    stop("No LD data found.\n")
  }

  if(is.null(facets)){facets <- ".base"}

  # check facets
  snp.facet <- .check.snpR.facet.request(x, facets, remove.type = "sample")[[1]]
  sample.facet <- .check.snpR.facet.request(x, facets, remove.type = "snp")[[1]]
  facets <- .check.snpR.facet.request(x, facets, remove.type = "none", return.type = T)
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
    x <- x[!apply(x, 1, function(y)all(is.na(y))), !apply(x, 2, function(y)all(is.na(y))), drop = FALSE]

    # set rownames as first column.
    x <- cbind(position = rownames(x), x)
    x <- data.table::as.data.table(x)


    #melting the df and fixing the names of the columns
    heatmap_x <- data.table::melt(x, id.vars = "position")
    names(heatmap_x) <- c("SNPa", "SNPb", "value")

    #getting rid of all the zeros from snps being compared to themselves
    heatmap_x <- stats::na.omit(heatmap_x)
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
      if(!is.null(sample.subfacet)){
        LD_mat_list <- LD_mat_list[which(names(LD_mat_list) %in% sample.subfacet)]
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

  SNPa <- SNPb <- value <- NULL
  
  # base plot
  out <- ggplot2::ggplot(LD_mats, ggplot2::aes(x = SNPa, y=SNPb, fill=value, color = value)) +
    ggplot2::geom_bin2d(stat = "identity")

  
  # wrapping
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

  # colors
  if(!is.null(gradient_colors)){
    out <- out +
      ggplot2::scale_color_gradient2(low = gradient_colors[1], high = gradient_colors[2]) +
      ggplot2::scale_fill_gradient2(low = gradient_colors[1], high = gradient_colors[2])
  }
  else{
    out <- out + 
      ggplot2::scale_fill_viridis_c(direction = -1, option = viridis.option) +
      ggplot2::scale_color_viridis_c(direction = -1, option = viridis.option)
  }
  
  # rest of the theme
  out <- out +
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
    ggplot2::ylab("Position (Mb)") + ggplot2::xlab("Position (Mb)") +
    ggplot2::guides(color = "none")

  if(!is.null(title)){
    out <- out + ggplot2::ggtitle(title)
  }

  return(out)
}

#'PCA, tSNE, and umap plots from snpRdata.
#'
#'Generate a ggplot cluster plot based on PCA, the Barnes-Hut simulation at
#'theta>0 implemented in \code{\link[Rtsne]{Rtsne}}, the Uniform Manifold
#'Approximation and Projection approach implemented in \code{\link[umap]{umap}},
#'or the Discriminant Analysis of Principal Components implemented in
#'\code{\link[adegenet]{dapc}}.
#'
#'
#'Works by conversion to the "sn" format described in \code{\link{format_snps}}
#'with interpolated missing genotypes for all methods other than DAPC.
#'
#'Cluster plots can be produced via, PCA, tSNE, umap, or DAPC. The PCA point
#'coordinates are calculated using \code{\link{prcomp}}. By default, the first
#'two principal coordinates are plotted. A PC matrix will also be returned for
#'easy plotting of other PCs. tSNE coordinates are calculated via
#'\code{\link[Rtsne]{Rtsne}}, which should be consulted to for more details
#'about this method. Stated simply, tSNE attempts to compress a
#'multi-dimensional PCA (PCs 1:n) into fewer dimensions while retaining as much
#'information as possible. As such, a tSNE plot can be seen as a representation
#'of many different PC axis compressed into a single two-dimensional plot. This
#'compression process is stochastic, and so plots will vary somewhat between
#'runs, and multiple runs are recommended. Uniform Manifold Approximation and
#'Projection (UMAP) coordinates are calculated via \code{\link[umap]{umap}}.
#'UMAP similarly attempts to reduce multi-dimensional results to a two
#'dimensional visualization. DAPC instead clusters individuals in \emph{n}
#'groups, a number which by default is interactively chosen (again using a
#'PCA framework).
#'
#'Note that clusters and relative positions of samples from both tSNE and UMAP
#'may not reliably represent the relationships present in the higher PCA
#'dimensions from which they are created. As such, it is probably not wise to
#'use these methods to draw conclusions about relationships. They are useful
#'exploratory tools, however, and so are kept available here.
#'
#'For more details on tSNE arguments, \code{\link[Rtsne]{Rtsne}} should be
#'consulted.
#'
#'Additional arguments to the UMAP can be also be provided. Additional
#'information on these arguments can be found in
#'\code{\link[umap]{umap.defaults}}.
#'
#'Data points for individuals can be automatically colored by any sample-level
#'facet categories. Facets should be provided as described in
#'\code{\link{Facets_in_snpR}}. Up to two different sample-level facets can be
#'automatically plotted simultaneously. If two facets are supplied, one level
#'will be noted by point shape and the other by color (by default the facet with
#'more options will be given shapes, behavior that can be controlled using the
#'\code{shape_has_more_levels} argument), as long as one has less than 6 total
#'levels. If both have more than 6 levels, one will be noted by point fill and
#'the other by point outline.
#'
#'@param x snpRdata object.
#'@param facets character, default NULL. Categorical sample-level metadata
#'  variables by which to color points. Up to two different sample-specific
#'  facets may be provided. See \code{\link{Facets_in_snpR}} for more details.
#'@param plot_type character, default "pca". c("pca", "tSNE", "umap", "dapc"). 
#'Types of
#'  plots to be produced. Options \itemize{\item{pca: } Principal Component
#'  Analysis, first two dimensions of variance. \item{tSNE: } t-Stochastic
#'  Neighbor Embedding, which collapses dims (see argument) dimensions of
#'  variance into two. \item{umap: } Uniform Manifold Approximation and
#'  Projection, which collapses multiple dimensions of variance into two.
#'  \item{dapc: } Discriminant analysis of principal components, 
#'  clusters individuals into groups for plotting via PCA.} See
#'  description for details.
#'@param check_duplicates logical, default FALSE. Checks for any duplicated
#'  individuals, which will cause errors. Since these rarely exist and
#'  drastically slow down function run-time, this defaults to FALSE.
#'@param minimum_percent_coverage numeric, default FALSE. Proportion of samples
#'  a SNP must be sequenced in to be used in generating plots.
#'@param minimum_genotype_percentage numeric, default FALSE. Proportion of SNPs
#'  a sample must be sequenced at in order to be used in plots.
#'@param smart_PCA logical, default TRUE. If TRUE, uses Patterson et. al 
#'  (2006)'s centering approach prior to plot construction. Note that this also
#'  avoids the need for interpolation, so interpolation is set to FALSE in this
#'  case.
#'@param interpolation_method character, default "bernoulli". Interpolation
#'  method to use for missing data. Options: \itemize{\item{bernoulli:
#'  }{Interpolated via binomial draw for each allele against minor allele
#'  frequency.} \item{af: }{Interpolated by inserting the expected number of
#'  minor alleles at missing data points given loci minor allele frequencies.}
#'  \item{iPCA: }{This an iterative PCA approach to interpolate based on SNP/SNP
#'  covariance via \code{\link[missMDA]{imputePCA}}. If the ncp argument is not
#'  defined, the number of components used for interpolation will be estimated
#'  using \code{\link[missMDA]{estim_ncpPCA}}. In this case, this method is much
#'  slower than the other methods, especially for large datasets. Setting an ncp
#'  of 2-5 generally results in reasonable interpolations without the time
#'  constraint.}} Ignored if \code{smart_PCA} is TRUE.
#'@param dims numeric, default 2. Output dimensionality, default 2.
#'@param initial_dims numeric, default 50. The number of dimensions retained in
#'  the initial PCA step during tSNE.
#'@param perplexity numeric, default FALSE. Perplexity parameter, by default
#'  found by \code{\link[mmtsne]{hbeta}}, with beta = 1.
#'@param theta numeric, default 0. Theta parameter from
#'  \code{\link[Rtsne]{Rtsne}}. Default an exhaustive search.
#'@param iter numeric, default 1000. Number of tSNE iterations/umap epochs to
#'  perform.
#'@param viridis.option character, default "viridis". Viridis color scale option
#'  to use for significance lines and SNP labels. See
#'  \code{\link[ggplot2]{scale_gradient}} for details.
#'@param alt.palette character or NULL, default NULL. Optional palette of colors
#'  to use instead of the viridis palette.
#'@param ncp numeric or NULL, default NULL. Number of components to consider for
#'  iPCA sn format interpolations of missing data. If null, the optimum number
#'  will be estimated, with the maximum specified by ncp.max. This can be very
#'  slow.
#'@param ncp.max numeric, default 5. Maximum number of components to check for
#'  when determining the optimum number of components to use when interpolating
#'  sn data using the iPCA approach.
#'@param dapc_clustering_max_n_clust numeric or NULL, default 20. If not NULL,
#'  the clustering parameters for DAPC calculation will be selected
#'  interactively, with \code{dapc_clustering_max_n_clust} max clusters
#'  considered. If NULL, the parameters \code{dapc_clustering_npca},
#'  \code{dapc_clustering_nclust}, \code{dapc_ndisc}, and \code{dapc_npca} must
#'  instead be set.
#'@param dapc_clustering_nclust numeric or NULL, default NULL. The number of
#'  clusters to use for DAPC. Interactive decision is recommended using
#'  \code{dapc_clustering_max_n_clust}.
#'@param dapc_clustering_npca numeric or NULL, default NULL. The number of PCS
#'  to use for assigning individuals to clusters with DAPC. Interactive decision
#'  is recommended using \code{dapc_clustering_max_n_clust}.
#'@param dapc_npca numeric or NULL, default NULL. The number of PCS to use for
#'  conducting the DAPC itself after assigning individuals to clusters.
#'  Interactive decision is recommended using
#'  \code{dapc_clustering_max_n_clust}.
#'@param dapc_ndisc numeric or NULL, default NULL. The number of discriminants
#'  to use for conducting the DAPC itself after assigning individuals to
#'  clusters. Interactive decision is recommended using
#'  \code{dapc_clustering_max_n_clust}.
#'@param ellipse_size numeric or NULL, default 1.5. The scaled-size of the
#'  ellipse to use for DAPC. If NULL, no ellipses will be calculated or drawn.
#'@param seg_lines logical, default TRUE. If TRUE, lines will be drawn between
#'  points and cluster centers when plotting with DAPC.
#'@param shape_has_more_levels logical, default TRUE. If TRUE and two facets are
#'  requested, the facet with more levels will plotted as shapes. If FALSE,
#'  the facet with less levels will be plotted with shapes. Ignored if the facet
#'  that would get shapes has more than 6 levels.
#'@param update_bib character or FALSE, default FALSE. If a file path to an
#'  existing .bib library or to a valid path for a new one, will update or
#'  create a .bib file including any new citations for methods used. Useful
#'  given that this function does not return a snpRdata object, so a
#'  \code{\link{citations}} cannot be used to fetch references.
#'@param verbose Logical, default FALSE. If TRUE, some progress updates may be
#'  reported.
#'@param ... Other arguments, passed to \code{\link[Rtsne]{Rtsne}} or
#'  \code{\link[umap]{umap}}.
#'
#'@return A list containing: ggplot PCA, tSNE, umap, and/or DAPC plots. May
#'  contain one to four objects, one for each PCA, tSNE, umap, or DAPC plot
#'  requested, named "pca" "tsne", "umap", and "dapc" respectively.
#'
#'@author William Hemstrom
#'@author Matt Thorstensen
#'
#'@references Jesse H. Krijthe (2015). Rtsne: T-Distributed Stochastic Neighbor
#'  Embedding using a Barnes-Hut Implementation, URL:
#'  \url{https://github.com/jkrijthe/Rtsne}.
#'@references Van Der Maaten, L. & Hinton, G. (2008) Visualizing
#'  high-dimensional data using t-SNE. journal of machine learning research.
#'  \emph{Journal of Machine Learning Research}.
#'@references McInnes, L. & Healy (2018). UMAP: uniform manifold approximation
#'  and projection. Preprint at URL: \url{https://arxiv.org/abs/1802.03426}.
#'@references Jombart, T., Devillard, S. & Balloux, F. Discriminant analysis of
#'  principal components: a new method for the analysis of genetically
#'  structured populations. BMC Genet 11, 94 (2010).
#'  \url{https://doi.org/10.1186}
#'
#'@seealso \code{\link[mmtsne]{mmtsne}}
#'@seealso \code{\link[umap]{umap}}
#'@seealso \code{\link[stats]{prcomp}}
#'@seealso \code{\link[adegenet]{dapc}}
#'
#'@export
#'
#' @examples
#' # plot colored by population
#' plot_clusters(stickSNPs, "pop")
#'
#' # plot colored by population and family
#' plot_clusters(stickSNPs, "pop.fam")
plot_clusters <- function(x, facets = NULL, plot_type = "pca", check_duplicates = FALSE,
                          minimum_percent_coverage = FALSE, minimum_genotype_percentage = FALSE,
                          smart_PCA = TRUE,
                          interpolation_method = "bernoulli",
                          dims = 2, initial_dims = 50, perplexity = FALSE, theta = 0, iter = 1000,
                          viridis.option = "viridis", alt.palette = NULL, ncp = NULL, ncp.max = 5,
                          dapc_clustering_max_n_clust = 20,
                          dapc_clustering_npca = NULL, dapc_clustering_nclust = NULL,
                          dapc_npca = NULL, dapc_ndisc = NULL, ellipse_size = 1.5,
                          seg_lines = TRUE, shape_has_more_levels = TRUE,
                          update_bib = FALSE, verbose = FALSE,...){

  #=============sanity checks============
  y <- cluster <- .cluster <- NULL
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
  facets <- .check.snpR.facet.request(x, facets, remove.type = "none", return.type = T)
  if(length(facets[[1]]) > 1){
    msg <- c(msg, "Only one facet may be specified at a time. This facet may be complex, with up to two sample levels (e.g. pop.family).")
  }
  if(any(facets[[2]] != "sample")){
    if(!(length(facets[[2]]) == 1 & facets[[2]][1] == ".base")){
      bf <- facets[[1]][which(facets[[2]] != "sample")]
      msg <- c(msg,  paste0("Only sample level facets can be plotted. Facet(s): ", paste0(bf, collapse = ", "), " refer to snp metadata."))
    }
  }
  facets <- facets[[1]]
  facets <- unlist(.split.facet(facets))
  facets <- .check.snpR.facet.request(x, facets)

  plot_type <- tolower(plot_type)
  good_plot_types <- c("pca", "tsne", "umap", "dapc")
  if(any(!(plot_type %in% good_plot_types))){
    msg <- c(msg, paste0("Unaccepted plot_type. Accepted types:", paste0(good_plot_types, collapse = ", "), "."))
  }

  if("umap" %in% plot_type){
    pkg.check <- .check.installed("umap")
    if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
  }
  
  if("tsne" %in% plot_type){
    pkg.check <- .check.installed("Rtsne")
    if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
    
    pkg.check <- .check.installed("mmtsne")
    if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
  }
  
  if("dapc" %in% plot_type){
    pkg.check <- .check.installed("adegenet")
    if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
    
    if(is.null(dapc_clustering_max_n_clust)){
      if(any(is.null(c(dapc_clustering_npca, dapc_clustering_nclust, dapc_npca, dapc_ndisc)))){
        msg <- c(msg, "If plot_clusters() is not run interactively please supply all dapc parameters.\n")
      }
    }
    else{
      if(!interactive()){
          msg <- c(msg, "If plot_clusters() is not run interactively, dapc interactive parameter picking cannot be used. Please set 'dapc_clustering_max_n_clust' to NULL and supply all dapc parameters.\n")
      }
    }
    
    if(any(c(is.null(dapc_npca), is.null(dapc_ndisc)))){
      if(!interactive()){
        msg <- c(msg, "If plot_clusters() is not run interactively please supply all dapc parameters.\n")
      }
      
      if(sum(c(is.null(dapc_npca), is.null(dapc_ndisc))) == 1){
        msg <- c(msg, "Please supply both dapc_npca and dapc_ndisc arguments or choose interactively instead.\n")
      }
    }
    
    if(sum(c(is.null(dapc_clustering_npca), is.null(dapc_clustering_nclust))) == 1){
      msg <- c(msg, "Please supply both dapc_clustering_npca and dapc_clustering_ndisc arguments or choose interactively instead.\n")
    }
  }
  
  if(!smart_PCA){
    if(isFALSE(interpolation_method) & length(plot_type) != 1 & plot_type[1] != "dapc"){
      stop("All methods require no missing data. Please enable interpolation.\n")
    }
  }
  else{
    interpolation_method <- FALSE
  }
  if(interpolation_method == "iPCA"){
    .check.installed("missMDA")
  }
  
  
  if(any(c("umap", "tsne") %in% plot_type)){
    cat("Note that clusters and relative placements of samples in UMAP and tSNE dimension reductions may not properly represent higher dimensionality clustering.\n",
        "While these methods are useful exploratory tools, they should probably not be used to draw conclusions.\n",
        "See https://doi.org/10.1101/2021.08.25.457696 for details (although currently still a pre-print).\n")
  }

  if(length(msg) > 0){
    stop(paste0(msg, collapse = "  \t"))
  }

  #=============prepare dataset===============
  cat("Formatting data...\n")

  bi_allelic <- .is.bi_allelic(x)
  if((length(plot_type) == 1 & plot_type[1] != "dapc") | length(plot_type) > 1){
    if(bi_allelic){
      # check for matching sn plot:
      if(length(x@sn$sn) != 0){
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
    }
    else{
      if(length(x@sn$sn) != 0){
        if(x@sn$type != interpolation_method){
          suppressWarnings(x@sn$sn <- format_snps(x, "pa", interpolate = interpolation_method, ncp = ncp, ncp.max = ncp.max))
          x@sn$type <- interpolation_method
        }
      }
      else{
        suppressWarnings(x@sn$sn <- format_snps(x, "pa", interpolate = interpolation_method, ncp = ncp, ncp.max = ncp.max))
        x@sn$type <- interpolation_method
      }
      
      sn <- x@sn$sn
      sn <- t(sn[,-c(1:(ncol(x@sample.meta) - 1))])
      meta <- x@sample.meta
    }
    
    
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
      if(verbose){cat("Checking for duplicates...\n")}
      sn <- as.data.table(t(sn))
      dups <- which(duplicated(sn) | duplicated(sn, fromLast=TRUE))
      if(length(dups) > 0){
        if(verbose){cat("Duplicates detected, indices:", dups, "\nRemoving all of these!\n")}
        sn <- sn[-dups,]
        meta <- meta[-dups,]
      }
      sn <- as.matrix(t(sn))
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
    
    if(verbose){cat("Plotting using", rm.snps, "loci and", rm.ind, "samples.\n")}
    sn <- t(sn)
    
    # do the centering if using smart_PCA approach
    if(smart_PCA){
      ms <- colMeans(sn, na.rm = TRUE)
      ms <- ms[col(sn)]
      afs <- (1+ms)/(2 + 2*nrow(sn)) # adjusted eqn, from Price et al 2006
      # afs <- ms/2 # according to patterson
      sn <- (sn - ms)/sqrt(afs*(1-afs)) # eqn 3, patterson et al 2006
      sn[is.na(sn)] <- 0
    }
  }
  else{
    meta <- x@sample.meta
  }
  

  #=============prepare plot data=============
  plot_dats <- vector("list", length(plot_type))
  names(plot_dats) <- plot_type

  if("pca" %in% plot_type){
    if(verbose){cat("Preparing pca...\n")}

    # did we do smart centering? if so, do svd instead of pca
    if(smart_PCA){
      pca_r <- svd(sn)
      colnames(pca_r$u) <- paste0("PC", 1:ncol(pca_r$u))
      plot_dats$pca <- cbind(meta, pca_r$u)
    }
    
    # or just pca directly?
    else{
      pca_r <- stats::prcomp(as.matrix(sn))
      pca <- as.data.frame(pca_r$x) #grab the PCA vectors.
      plot_dats$pca <- cbind(meta, pca)  #add metadata that is present in the input.
    }
    
  }
  if("tsne" %in% plot_type){
    #get perplexity if not provided
    if(perplexity == FALSE){
      if(verbose){cat("Estimating perplexity...\n")}
      perp <- mmtsne::hbeta(sn, beta = 1)
      perplexity <- perp$H
    }

    #run the tSNE
    if(verbose){cat("Running tSNE...\n")}
    tsne.out <- Rtsne::Rtsne(X = sn, dims = dims, initial_dims = initial_dims, perplexity = perplexity,
                             theta = theta, max_iter = iter, check_duplicates = FALSE,
                             verbose=verbose, ...)
    colnames(tsne.out$Y) <- paste0("PC", 1:ncol(tsne.out$Y))
    plot_dats$tsne <- cbind(meta, as.data.frame(tsne.out$Y))
  }
  if("umap" %in% plot_type){
    if(verbose){cat("Preparing umap...\n")}
    umap_r <- umap::umap(as.matrix(sn), n_epochs = iter, verbose = verbose, ...)
    colnames(umap_r$layout) <- paste0("PC", 1:ncol(umap_r$layout))
    plot_dats$umap <- cbind(meta, as.data.frame(umap_r$layout))
  }
  if("dapc" %in% plot_type){
    gi <- format_snps(x, "adegenet", "pop")
    if(!is.null(dapc_clustering_max_n_clust)){
      k <- adegenet::find.clusters.genind(x = gi, max.n.clust = dapc_clustering_max_n_clust)
    }
    else{
      k <- adegenet::find.clusters.genind(x = gi, choose.n.clust = FALSE,
                                          n.pca = dapc_clustering_npca, n.clust = dapc_clustering_nclust)
    }
    if(any(is.null(c(dapc_npca, dapc_ndisc)))){
      plot_dats$dapc <- adegenet::dapc(gi, k$grp)
    }
    else{
      plot_dats$dapc <- adegenet::dapc.genind(gi, k$grp, n.pca = dapc_npca, n.da = dapc_ndisc)
      plot_dats$dapc
    }
  }
  #=============make ggplots=====================
  plots <- vector("list", length(plot_dats))
  names(plots) <- names(plot_dats)

  if(length(facets) > 2){
    warning("Only up to two simultanious sample metadata columns can be plotted at once.\n")
  }

  for(i in 1:length(plot_dats)){
    
    if(names(plot_dats)[i] == "dapc"){
      
      
      # from ade4, two internals I need here to get ellipses, the second with edits to return values instead of plotting them
      fac2disj <- function(fac, drop = FALSE) {
        ## Returns the disjunctive table corrseponding to a factor
        n <- length(fac)
        fac <- as.factor(fac)
        if(drop)
          fac <- factor(fac)
        x <- matrix(0, n, nlevels(fac))
        x[(1:n) + n * (unclass(fac) - 1)] <- 1
        dimnames(x) <- list(names(fac), as.character(levels(fac)))
        return(data.frame(x, check.names = FALSE))
      }
      scatterutil.ellipse <- function (x, y, z, cellipse) {
        if (any(is.na(z))) 
          return(invisible())
        if (sum(z * z) == 0) 
          return(invisible())
        util.ellipse <- function(mx, my, vx, cxy, vy, coeff) {
          lig <- 100
          epsi <- 1e-10
          x <- 0
          y <- 0
          if (vx < 0) 
            vx <- 0
          if (vy < 0) 
            vy <- 0
          if (vx == 0 && vy == 0) 
            return(NULL)
          delta <- (vx - vy) * (vx - vy) + 4 * cxy * cxy
          delta <- sqrt(delta)
          l1 <- (vx + vy + delta)/2
          l2 <- vx + vy - l1
          if (l1 < 0) 
            l1 <- 0
          if (l2 < 0) 
            l2 <- 0
          l1 <- sqrt(l1)
          l2 <- sqrt(l2)
          test <- 0
          if (vx == 0) {
            a0 <- 0
            b0 <- 1
            test <- 1
          }
          if ((vy == 0) && (test == 0)) {
            a0 <- 1
            b0 <- 0
            test <- 1
          }
          if (((abs(cxy)) < epsi) && (test == 0)) {
            if(vx > vy){
              a0 <- 1
              b0 <- 0
            } else {
              a0 <- 0
              b0 <- 1
            }
            test <- 1
          }
          if (test == 0) {
            a0 <- 1
            b0 <- (l1 * l1 - vx)/cxy
            norm <- sqrt(a0 * a0 + b0 * b0)
            a0 <- a0/norm
            b0 <- b0/norm
          }
          a1 <- 2 * pi/lig
          c11 <- coeff * a0 * l1
          c12 <- (-coeff) * b0 * l2
          c21 <- coeff * b0 * l1
          c22 <- coeff * a0 * l2
          angle <- 0
          for (i in 1:lig) {
            cosinus <- cos(angle)
            sinus <- sin(angle)
            x[i] <- mx + c11 * cosinus + c12 * sinus
            y[i] <- my + c21 * cosinus + c22 * sinus
            angle <- angle + a1
          }
          return(list(x = x, y = y, seg1 = c(mx + c11, my + c21, 
                                             mx - c11, my - c21), seg2 = c(mx + c12, my + c22, 
                                                                           mx - c12, my - c22)))
        }
        z <- z/sum(z)
        m1 <- sum(x * z)
        m2 <- sum(y * z)
        v1 <- sum((x - m1) * (x - m1) * z)
        v2 <- sum((y - m2) * (y - m2) * z)
        cxy <- sum((x - m1) * (y - m2) * z)
        ell <- util.ellipse(m1, m2, v1, cxy, v2, cellipse)
        if (is.null(ell)) 
          return(invisible())
        else{
          return(ell)
        }
      }
      
      tpd <- as.data.frame(plot_dats[[i]]$ind.coord)
      colnames(tpd) <- gsub("LD", ".LD", colnames(tpd))
      tpd$.cluster <- plot_dats[[i]]$grp
      tpd <- data.table::as.data.table(tpd)
      
      # xs <- pd[ , cbind(density(LD1)[1:2]), by = cluster]
      # colnames(xs)[2:3] <- c("x", "y")
      # densities <- tapply(pd$LD1, pd$cluster, density)
      
      tpd <- as.data.table(cbind(meta, tpd))
      
      .LD1 <- .LD2 <- .GRP_LD1 <- .GRP_LD2 <- NULL
      
      out <- ggplot2::ggplot(tpd, ggplot2::aes(.LD1, .LD2)) + ggplot2::theme_bw()
      
      if(is.numeric(ellipse_size)){
        disj <- fac2disj(as.factor(tpd$.cluster))
        ellipses <- vector("list", length(k$size))
        for(j in 1:length(ellipses)){
          ellipses[[j]] <- scatterutil.ellipse(tpd$.LD1, tpd$.LD2, disj[,j], ellipse_size)
          ellipses[[j]] <- data.frame(x = ellipses[[j]]$x, y = ellipses[[j]]$y)
          ellipses[[j]]$cluster <- j
        }
        ellipses <- data.table::rbindlist(ellipses)
        
        out <- out + 
          ggplot2::geom_polygon(data = ellipses, mapping = ggplot2::aes(x, y, group = as.factor(cluster)), fill = NA, color = "black")
      }
      
      if(seg_lines){
        tpd[,c(".GRP_LD1", ".GRP_LD2") := as.data.frame(plot_dats[[i]]$grp.coord)[as.numeric(as.character(.cluster)),1:2]]
        tpd <- as.data.frame(tpd)
        out <- out + 
          ggplot2::geom_segment(ggplot2::aes(xend = .GRP_LD1, yend = .GRP_LD2))
      }
      
    }
    
    else{
      tpd <- plot_dats[[i]]
      
      #Categories (pops, fathers, mothers, etc.) are given in plot.vars argument. Supports up to two!
      #make the base plot, then add categories as color and fill.
      PC1 <- PC2 <- NULL
      out <- ggplot2::ggplot(tpd, ggplot2::aes(PC1, PC2)) + ggplot2::theme_bw() #initialize plot
      
    }
    

    if(facets[1] == ".base"){
      out <- out + ggplot2::geom_point()
    }
    else{
      #add variables.
      v1 <- tpd[,which(colnames(tpd) == facets[1])] #get the factors
      
      
      v2_has_color <- FALSE
      
      
      # add geoms to plot
      if(length(facets) == 1){
        out <- out + ggplot2::geom_point(ggplot2::aes(color = v1))#add the factor
      }
      else{
        # if two plotting variables, prepare the second and add it as well.
        v2 <- tpd[,which(colnames(tpd) == facets[2])]
        
        # Figure out which has less levels, make that v2 IF neither have more than 6.
        lv1 <- length(unique(v1))
        lv2 <- length(unique(v2))
        
        # adjust levels
        if(lv1 <= 6 | lv2 <= 6){
          # if need to be flipped
          if(lv1 > lv2){
            # only flip if v2 has less or equal to 6 levels and we are flipping
            if(!lv2 > 6 & shape_has_more_levels){
              temp <- v2
              v2 <- v1
              v1 <- temp
              facets <- rev(facets)
              
              temp <- lv2
              lv2 <- lv1
              lv1 <- temp
              rm(temp)
            }
          }
          # also flip if lv1 is less than lv1 and less than 7, but lv2 has too many levels
          else if(lv2 > 6 & lv1 <=6){
            temp <- v2
            v2 <- v1
            v1 <- temp
            facets <- rev(facets)
            
            temp <- lv2
            lv2 <- lv1
            lv1 <- temp
            rm(temp)
          }
        }
        
        if(lv2 <= 6){
          out <- out + ggplot2::geom_point(ggplot2::aes(color = v1, shape = v2), size = 2.5, stroke = 1.25) +
            ggplot2::scale_shape_discrete(name = facets[2])
        }
        else{
          out <- out + ggplot2::geom_point(ggplot2::aes(color = v1, fill = v2), pch = 21, size = 2.5, stroke = 1.25)
          v2_has_color <- TRUE
        }
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
        if(v2_has_color){
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
      }
    }


    

    if(plot_type[i] == "pca"){
      if(smart_PCA){
        loadings <- (pca_r$d^2)/sum(pca_r$d^2)
        loadings <- round(loadings * 100, 2)
      }
      else{
        loadings <- (pca_r$sdev^2)/sum(pca_r$sdev^2)
        loadings <- round(loadings * 100, 2)
      }
      out <- out + ggplot2::xlab(paste0("PC1 (", loadings[1], "%)")) + ggplot2::ylab(paste0("PC2 (", loadings[2], "%)"))
    }
    else{
      out <- out + ggplot2::xlab("Dim 1") + ggplot2::ylab("Dim 2")
    }
    plots[[i]] <- out
  }
  
  # cite
  keys <- character(0)
  stats <- character(0)
  details <- character(0)
  if(smart_PCA & !("dapc" %in% plot_type & length(plot_type) == 1)){
    keys <- c(keys, "pricePrincipalComponentsAnalysis2006", "pattersonPopulationStructureEigenanalysis2006")
    stats <- c(stats, "smart PCA", "smart PCA")
    details <- c(details, "afs estimation during scaling and centering", "smart PCA centering and scaling adjustment")
  }
  
  if("tsne" %in% plot_type){
    keys <- c(keys, "Krijthe2015", "Maatan2008")
    stats <- c(stats, "RtSNE", "tSNE")
    details <- c(details, "R package used to conduct tSNE", "t-stochastic Neighbor Embedding (tSNE) citation")
  }
  if("umap" %in% plot_type){
    keys <- c(keys, "McInnes2018")
    stats <- c(stats, "UMAP")
    details <- c(details, "Uniform Manifold Approximation and Projection (UMAP)")
  }
  
  if("dapc" %in% plot_type){
    keys <- c(keys, "Jombart2010")
    stats <- c(stats, "DAPC")
    details <- c(details, "Discriminant Analysis of Principal Components core method.")
    
    keys <- c(keys, "Jombart2008")
    stats <- c(stats, "adegenet")
    details <- c(details, "DAPC run via interface to adegenet.")
  }

  if(length(keys) > 0){
    .yell_citation(keys, stats, details, update_bib)
  }
  
  return(plots)
}

#' Generate a manhattan plot from snpRdata or a data.frame.
#'
#' Creates a ggplot-based manhattan plot, where chromosomes/scaffolds/etc are
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
#' A column defining the position of the SNP within the chromosome must be
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
#'   calc_smoothed_averages. Ignored if x is a data.frame.
#' @param facets character or NULL, default NULL. Facets by which to break
#'   plots, as described in \code{\link{Facets_in_snpR}}. For non-window stats,
#'   the any snp metadata facets will be ignored. Ignored if x is a data.frame.
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
#' @param sig_below logical, default FALSE. If TRUE, treats values lower than
#'   the significance threshold as significant.
#' @param log.p logical, default FALSE. If TRUE, plot variables and thresholds
#'   will be transformed to -log.
#' @param abs logical, default FALSE. If TRUE, converts the plot variable to
#'   it's absolute value.
#' @param highlight character, numeric, or FALSE, default "significant".
#'   Controls SNP highlighting. If either "significant" or "suggestive", SNPs
#'   above those respective values will be highlighted. If a numeric vector, SNPs
#'   corresponding to vector entries will be highlighted. See details.
#' @param highlight_style character, default "label". Highlighting options:
#'   \itemize{\item{label: }{labels with chr and position.}\item{color: }{Color
#'   (word or hex) to color points by.}}
#' @param viridis.option character, default "plasma". Viridis color scale option
#'   to use for significance lines and SNP labels. See
#'   \code{\link[ggplot2]{scale_gradient}} for details.
#' @param viridis.hue numeric, default c(0.2, 0.5). Two values between 0 and 1
#'   listing the hues at which to start and stop on the viridis palette defined
#'   by the viridis.option argument. Lower numbers are darker.
#' @param t.sizes numeric, default c(16, 12, 10). Text sizes, given as
#'   c(strip.title, axis, axis.ticks).
#' @param colors character, default c("black", "slategray3"). Colors to
#'   alternate across chromosomes.
#' @param rug_data data.frame or tbl, default NULL. Data to plot as a rug below
#'   the manhattan plot containing columns named to match the \code{chr}
#'   argument \emph{and} either the \code{bp} argument OR columns named
#'   \code{start} and \code{end} as well as, optionally, a column named to match
#'   the \code{rug_label} column. Useful for labeling the locations of candidate
#'   genes, for example.
#' @param rug_style character, default "point". Options for the style of the
#'   rug, ignored if \code{rug_data} is not provided. Options:
#'   \itemize{\item{point: } standard rug plot with vertical dashes below the
#'   plot at the indicated locations. If start and end points are supplied in
#'   \code{rug_data}, the midpoint will be plotted. \item{ribbon: } Ribbons for
#'   each point drawn below the plot, from the \code{start} to \code{end}
#'   columns. If the plotted range of \code{x} is very large (as in whole-genome
#'   or reduced representation sequencing), these may not be visible. A warning
#'   will be provided if this may be the case. Sub-setting the input and rug
#'   data to a range of interest may help in this case.}
#' @param rug_label character, default NULL. Names of additional labeling
#'   columns in \code{rug_data}, ignored if \code{rug_data} is not provided.
#'   \emph{These will not be directly plotted} (since the result is often very
#'   messy), but are available as aesthetics in the resulting plot, which can
#'   then be examined if something like the \code{ggplotly} function from
#'   \code{plotly} is used. This may change in the future if a clean plotting
#'   technique is suggested.
#' @param rug_alpha numeric between 0 and 1, default 0.3. Alpha (transparency)
#'   applied to a ribbon-style rug. Ignored if \code{rug_data} is not provided
#'   or the \code{rug_style} is not \code{ribbon}.
#' @param rug_thickness numeric or \code{grid}-style \code{unit}, default
#'   \code{ggplot2::unit(ifelse(rug_style == "point", 0.03, 6), "npc")}. The
#'   height of the rug lines (if \code{rug_style = "point"}) or ribbon (if
#'   \code{rug_style = "ribbon"}). Ignored if \code{rug_data} is not provided.
#'   Use of the \code{\link[ggplot2]{unit}} style of size choice recommended to
#'   avoid over-plotting.
#' @param chr_order character, default NULL. If provided, an ordered vector of
#'   chromosome/scaffold/etc names by which to sort output.
#' @param abbreviate_labels numeric or FALSE, default FALSE. If a numeric value,
#'   x-axis chromosome names will be abbreviated using
#'   \code{\link[base]{abbreviate}}, with each abbreviated label having the
#'   minimum length specified. Helpful when chromosome/scaffold/etc names are
#'   very long.
#' @param lambda_gc_correction Correct for inflated significance due to
#'   population and/or family structure using the \eqn{\gamma_{GC}} approach
#'   described in Price et al 2010.
#'
#' @references Price, A., Zaitlen, N., Reich, D. et al. New approaches to
#'   population stratification in genome-wide association studies. Nat Rev Genet
#'   11, 459–463 (2010). https://doi.org/10.1038/nrg2813
#'
#' @author William Hemstrom
#' @export
#'
#' @return A ggplot manhattan plot.
#'
#'
#' @examples
#' # add a dummy phenotype and run an association test.
#' x <- stickSNPs[pop = c("ASP", "SMR"), chr = c("groupIX", "groupIV")]
#' sample.meta(x)$phenotype <- sample(c("A", "B"), nsamps(x), TRUE)
#' x <- calc_association(x, response = "phenotype", method = "armitage")
#' plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
#'                log.p = TRUE)$plot
#' 
#' 
#' # other types of stats:
#' # make some data
#' x <- calc_basic_snp_stats(x, "pop.chr", sigma = 200, step = 50)
#' 
#' # plot pi, breaking apart by population, keeping only the groupIX
#' # and the ASP population, with
#' # significant and suggestive lines plotted and SNPs
#' # with pi below the significance level labeled.
#' plot_manhattan(x, "pi", facets = "pop",
#'                chr = "chr", chr.subfacet = "groupIX",
#'                sample.subfacet = "ASP",
#'                significant = 0.05, suggestive = 0.15, sig_below = TRUE)$plot
#' 
#' # plot FST for the ASP/SMR comparison across all chromosomes,
#' # labeling the first 10 SNPs in x (by row) with their ID
#' # Note that since this is thie ony comparison, we don't actually need to
#' # specify it.
#' plot_manhattan(x, "fst", facets = "pop.chr",
#'                sample.subfacet = "ASP~SMR", highlight = 1:10,
#'                chr = "chr", snp = ".snp.id")$plot
#' 
#' # plot sliding-window FST between ASP and SMR
#' # and between OPL and SMR
#' plot_manhattan(x, "fst", window = TRUE, facets = c("pop.chr"),
#'                chr = "chr", sample.subfacet = "ASP~SMR",
#'                significant = .29, suggestive = .2)$plot
#' 
#' # plot using a data.frame,
#' # using log-transformed p-values
#' ## grab data
#' y <- get.snpR.stats(x, "pop", stats = "hwe")$single
#' ## plot
#' plot_manhattan(y, "pHWE", facets = "subfacet", chr = "chr",
#'                significant = 0.0001, suggestive = 0.001,
#'                log.p = TRUE, highlight = FALSE)$plot
#' 
#' 
#' 
#' # plot with a rug
#' rug_data <- data.frame(chr = c("groupIX", "groupIV"), start = c(0, 1000000),
#'                        end = c(5000000, 6000000), gene = c("A", "B"))
#' 
#' # point style, midpoints plotted
#' plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
#'                log.p = TRUE, rug_data = rug_data)
#' 
#' # ribbon style
#' plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
#'                log.p = TRUE, rug_data = rug_data, rug_style = "ribbon")
#'                
#' # with plotly to mouse over information
#' \dontrun{
#' plotly::ggplotly(plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
#'                                 log.p = TRUE, rug_data = rug_data, 
#'                                 rug_style = "ribbon", 
#'                                 rug_label = "gene")$plot)
#' 
#' }
plot_manhattan <- function(x, plot_var, window = FALSE, facets = NULL,
                           chr = "chr", bp = "position", snp = NULL,
                           chr.subfacet = NULL, sample.subfacet = NULL,
                           significant = NULL, suggestive = NULL,
                           highlight = "significant",
                           highlight_style = "label",
                           sig_below = FALSE, log.p = FALSE, abs = FALSE,
                           viridis.option = "plasma", viridis.hue = c(.2, 0.5), t.sizes = c(16, 12, 10),
                           colors = c("black", "slategray3"),
                           rug_data = NULL,
                           rug_style = "point",
                           rug_label = NULL,
                           rug_alpha = 0.3,
                           rug_thickness = ggplot2::unit(ifelse(rug_style == "point", 0.03, 6), "npc"),
                           lambda_gc_correction = FALSE,
                           chr_order = NULL,
                           abbreviate_labels = FALSE){

  cum.bp <- cum.start <- cum.end <- y <- start <- end <-  NULL

  #=============sanity checks==============================
  msg <- character()
  if(highlight_style == "label" & !isFALSE(highlight)){
    pkg.check <- .check.installed("ggrepel")
  }
  
  .check.installed("viridis")
  
  # rug sanity checks, can do those now
  if(!is.null(rug_data)){
    if(!rug_style %in% c("point", "ribbon")){
      msg <- c(msg, ("Unrecognized rug_style argument. Recognized options: 'point', 'ribbon'.\n"))
    }
    if(!is.data.frame(rug_data) | methods::is(rug_data, "tbl")){
      msg <- c(msg, "rug_data must be a data.frame or tbl.\n")
    }
    else{
      
      
      if(!chr %in% colnames(rug_data)){
        msg <- c(msg, paste0(chr, " column not found in rug_data.\n"))
      }
      
      
      
      if(!all(c("start", "end") %in% colnames(rug_data))){
        if(rug_style == "ribbon"){
          msg <- c(msg, "Could not locate start and/or end columns in rug_data. Both are needed for ribbon plotting.\n")
        }
        else if(!bp %in% colnames(rug_data)){
          msg <- c(msg, paste0("Neither ", bp, " or start+end columns not found in rug_data.\n"))
        }
      }
      else if(rug_style == "point" & !bp %in% colnames(rug_data)){
        rug_data[,bp] <- stats::ave(rug_data$start, rug_data$end) # make the bp column if needed and possible!
      }
      
      
      
      if(rug_style == "ribbon"){
        if(rug_alpha < 0 | rug_alpha > 1){
          msg <- c(msg, "rug_alpha must be between 0 and 1.\n")
        }
      }
      
      
      if(!is.null(rug_label)){
        if(!rug_label %in% colnames(rug_data)){
          msg <- c(msg, paste0(rug_label, " not found in rug_data.\n"))
        }
      }
    }
  }
  
  if(length(msg) > 0){
    stop(msg, collapse = "\n")
  }
  
  #=============grab the desired stats=====================
  #====if a snpRdata object========
  if(is.snpRdata(x)){
    if(window == FALSE){
      facets <- .check.snpR.facet.request(x, facets)
    }
    else{
      if(!is.null(facets)){
        pop.facets <- .check.snpR.facet.request(x, facets, "snp")
        if(any(pop.facets == ".base")){
          pop.facets <- pop.facets[-which(pop.facets == ".base")]
        }
        if(length(pop.facets) > 0){
          facets <- paste0(pop.facets, ".", chr)
        }
        else{
          facets <- chr
        }
        
      }
      else{
        facets <- chr
      }

      facets <- .check.snpR.facet.request(x, facets, "none")
    }
    if(plot_var %in% colnames(x@stats)){
      if(window){

        stats <- .get.snpR.stats(x, facets = facets, type = "single.window")
        colnames(stats)[which(colnames(stats) == "snp.subfacet")] <- chr
      }
      else{
        stats <- .get.snpR.stats(x, facets = facets)
      }
    }
    else if(plot_var %in% colnames(x@pairwise.stats)){
      if(window){
        stats <- .get.snpR.stats(x, facets, "pairwise.window")
        if(!chr %in% colnames(stats)){
          colnames(stats)[which(colnames(stats) == "snp.subfacet")] <- chr
        }
      }
      else{
        stats <- .get.snpR.stats(x, facets, "pairwise")
      }
    }
    else if(window & plot_var %in% colnames(x@window.stats)){
      stats <- .get.snpR.stats(x, facets, "single.window")
      colnames(stats)[which(colnames(stats) == "snp.subfacet")] <- chr
    }
    else{
      stop("Unable to locate stat: ", plot_var, " in the provided data. Did you remember to run this statistic?\n")
    }
    if(is.null(stats)){
      stop("No matching statistics. Did you remember to smooth by your chromosome/scaffold/etc?\n")
    }

    if(nrow(stats) == 0){
      stop("No matching statistics. Did you remember to smooth by your chromosome/scaffold/etc?\n")
    }
  }
  #====otherwise=====
  else if(is.data.frame(x)){
    if(data.table::is.data.table(x)){x <- as.data.frame(x)}
    stats <- x[,which(colnames(x) %in% c(plot_var, bp, chr))]
    
    # deal with facets provided--may be column names or snpR style facets
    if(!is.null(facets)){
      
      # if length = 1, check if it's in snpR style and confirm everything is OK
      if(length(facets) == 1){
        
        if(facets %in% colnames(x)){
          stats$subfacet <- x[,facets]
        }
        
        else if(grepl("\\.", facets)){
          check <- .split.facet(facets)[[1]]
          
          if(all(check %in% colnames(x))){
            stats$subfacet <- .paste.by.facet(x, check)
          }
          else{
            stop(paste0("Some facets not found in column names of x. Bad facets: ", 
                        paste0(check[which(! check %in% colnames(x))], collapse = ", "),
                        "\n"))
          }
        }
        else{
          stop("Facet not found in column names of x.\n")
        }
        
      }
      
      # other wise check OK
      else if(any(!facets %in% colnames(x))){
        stop(paste0("Some facets not found in column names of x. Bad facets: ", 
                    paste0(facets[which(! facets %in% colnames(x))], collapse = ", "),
                    "\n"))
      }
      else{
        stats$subfacet <- .paste.by.facet(x, facets)
      }
    }
  }
  else{
    stop("x must be a data.frame or snpRdata object.\n")
  }
  
  #================sanity checks==============
  msg <- character(0)
  
  if(nrow(stats) != 0){
    nas <- which(is.na(stats[,plot_var]))
    if(length(nas) != 0){
      stats <- stats[-which(is.na(stats[,plot_var])),]
    }
  }
  
  
  if(nrow(stats) == 0){
    stop("No matching statistics. Did you remember to smooth by your chromosome/scaffold/etc?\n")
  }
  
  
  
  if(!bp %in% colnames(stats)){
    msg <- c(msg, paste0("Position column: ", bp, " not found in data/snp.meta. Define with argument bp = \n"))
  }
  
  if(!chr %in% colnames(stats)){
    msg <- c(msg, paste0("Chromosome column: ", chr, " not found in data/snp.meta. Define with argument chr = \n"))
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  
  
  #================ask user to pick option if window and multiple schemes======
  if(window){
    opts <- unique(stats[,c("sigma", "step", "nk.status", "gaussian", "triple_sigma")])
    if(nrow(opts) > 1){
      message("Multiple window schemes detected.\n")
      
      if(interactive()){
        message("Which would you like to use?")
        rownames(opts) <- 1:nrow(opts)
        print(opts)
        resp <- readline(prompt = "Select row (by number):")
        while(!resp %in% 1:nrow(opts)){
          resp <- readline(prompt = "Select row (by number):")
        }
      }
      else{
        resp <- 1
        message("Using:\n", paste0(colnames(opts), " = ", opts[1,], " "))
      }
      
      resp <- as.numeric(resp)
      match_opts <- .paste.by.facet(stats, colnames(opts), sep = "_")
      opts <- .paste.by.facet(opts, colnames(opts), sep = "_")
      stats <- stats[which(match_opts %in% opts[resp]),]
      rm(match_opts, opts, resp)
    }
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

  if(!is.null(chr_order)){
    stats[,chr] <- factor(stats[,chr], levels = chr_order)
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
  
  # apply the same adjustment to the rug_data if provided
  if(!is.null(rug_data)){
    rug_data$start.cum.bp <- cum.bp[match(rug_data[,chr], names(cum.bp))]
    if(rug_style == "point"){
      rug_data$cum.bp <- rug_data[,bp] + rug_data$start.cum.bp
    }
    else{
      rug_data$cum.start <- rug_data$start + rug_data$start.cum.bp
      rug_data$cum.end <- rug_data$end + rug_data$start.cum.bp
    }
  }

  #=============clean up==================
  colnames(stats)[which(colnames(stats) == plot_var)] <- "pvar"
  colnames(stats)[which(colnames(stats) == chr)] <- "chr"
  
  # lambda gc correction for pop structure, etc.
  if(lambda_gc_correction){
    .q <- .lam <- pvar <- NULL
    stats <- data.table::as.data.table(stats)
    stats[,.q := stats::qchisq(1-pvar,1)]
    
    # lam will differ depending on grouping by pops.
    if(length(unique(as.character(stats$subfacet))) > 1){
      subfacet <- NULL
      stats[,.lam := stats::median(.q)/stats::qchisq(0.5,1), by = subfacet]
    }
    else{
      stats[,.lam := stats::median(.q)/stats::qchisq(0.5,1)]
    }
    
    stats[,pvar := stats::pchisq(.q/.lam, 1, lower.tail = FALSE)]
    
  }
  
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
        if(!isFALSE(abbreviate_labels)){
          stats$highlight.label <- paste0(abbreviate(stats$chr, minlength = abbreviate_labels) , "_", stats[,bp])
        }
        else{
          stats$highlight.label <- paste0(stats$chr, "_", stats[,bp])
        }
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
  pvar <- NULL
  p <- ggplot2::ggplot(stats, ggplot2::aes(x = cum.bp, y = pvar, color = as.factor(chr))) +
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
    ggplot2::scale_color_manual(values = rep(c(colors), length(cum.chr.centers))) +
    ggplot2::xlab(chr) + ggplot2::ylab(plot_var)

  #=============adjust the plot========
  if(length(unique(as.character(stats$subfacet))) > 1){
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
    highlight.label <- NULL
    if(highlight_style == "label"){
      p <- p + ggrepel::geom_label_repel(data = stats[which(stats$highlight == 1),],
                                         mapping = ggplot2::aes(label = highlight.label), color = add.palette[3],
                                         force = 1.3)
    }
    else{
      p <- p + ggplot2::geom_point(data = stats[which(stats$highlight == 1),],
                                         color = highlight_style)
    }
    
  }

  # abbreviate
  if(!isFALSE(abbreviate_labels)){
    p <- p + ggplot2::scale_x_continuous(label = abbreviate(names(cum.chr.centers), 
                                                            minlength = abbreviate_labels), 
                                         breaks = cum.chr.centers, minor_breaks = NULL)
  }
  else{
    p <- p + ggplot2::scale_x_continuous(label = names(cum.chr.centers), 
                                         breaks = cum.chr.centers, minor_breaks = NULL)
  }
  
  # rug
  if(!is.null(rug_data)){

    # standard rug style
    if(rug_style == "point"){
      # note: labels aren't plotted because they tend to be very messy, so they are returned as an unplotted aesthetic for plotly use!
      if(!is.null(rug_label)){
        ochr <- chr
        orl <- rug_label
        obp <- bp
        
        chr <- ggplot2::sym(chr)
        rug_label <- ggplot2::sym(rug_label)
        bp <- ggplot2::sym(bp)
        
        .suppress_specific_warning(p <- p + ggplot2::geom_rug(data = rug_data, 
                                                              mapping = ggplot2::aes(label = rug_label, position = bp, color = chr),
                                                              length = rug_thickness),
                                   "Ignoring unknown aesthetics")
        
        chr <- ochr
        rug_label <- orl
        bp <- obp
      }
      else{
        ochr <- chr
        
        chr <- ggplot2::sym(chr)
        
        p <- p + ggplot2::geom_rug(data = rug_data, ggplot2::aes(color = chr))
        
        chr <- ochr
      }
    }
    
    # ribbon style -- note that this really doesn't make sense if you are plotting genome-wide data, so warn if a large x range
    else{
      p_min <- min(stats$pvar, na.rm = TRUE)
      p_range <- range(stats$pvar, na.rm = TRUE)
      p_range <- abs(p_range[2] - p_range[1])
      rug_ymin <- ifelse(p_min < 0, p_min + p_range*.05, p_min - p_range*.07)
      rug_ymax <- ifelse(p_min < 0, p_min + p_range*.07, p_min - p_range*.05)
      rug_med <- mean(rug_ymin, rug_ymax)
      
      rug_data$y <- rug_med
      
      if(min(abs(rug_data$cum.start - rug_data$cum.end), na.rm = TRUE) < max(stats$cum.bp)*0.005){
        warning("Some ribbon segments are very small and may not be visible. Consider using the 'point' rug_style or subset the input data down to a smaller positional range.")
      }
      

      if(!is.null(rug_label)){
        orl <- rug_label
        ochr <- chr
        
        rug_label <- ggplot2::sym(rug_label)
        chr <- ggplot2::sym(chr)
        
        
        .suppress_specific_warning(
          p <- p + ggplot2::geom_segment(data = rug_data, 
                                         mapping = ggplot2::aes(x = cum.start, 
                                                                xend = cum.end, 
                                                                y = y,
                                                                yend = y, 
                                                                label = rug_label,
                                                                start_position = start,
                                                                end_position = end,
                                                                color = chr),
                                         linewidth = rug_thickness), 
          "Ignoring unknown aesthetics")
        
        chr <- ochr
        rug_label <- orl
        
      }
      else{
        ochr <- chr
        
        chr <- ggplot2::sym(chr)
        
        p <- p + ggplot2::geom_segment(data = rug_data, 
                                       mapping = ggplot2::aes(x = cum.start, xend = cum.end, y = y, yend = y, color = chr),
                                       linewidth = rug_thickness)
        
        chr <- ochr
      }
    }
  }
  
  return(p)
}


#' Generate qq plots from p values
#'
#' Generate qq (quantile-quantile) plots for p-values, splitting by provided
#' facets, using either a snpRdata object or a data.frame. Produces a
#' \code{\link[ggplot2]{ggplot}} object, which is modifiable as usual.
#'
#' @param x snpRdata object or data.frame containing p-values to plot.
#' @param plot_var character, name of the p-values to plot. For a snpRdata
#'   object, the name should refer to the name of the column produced when
#'   fetching data via \code{\link{get.snpR.stats}}. For a data.frame, it should
#'   refer to the column name directly.
#' @param facets character, default NULL. Facets by which to split the plot. See
#'   \code{\link{Facets_in_snpR}}.
#' @param lambda_gc_correction Correct for inflated significance due to
#'   population and/or family structure using the \eqn{\gamma_{GC}} approach
#'   described in Price et al 2010.
#'
#' @references Price, A., Zaitlen, N., Reich, D. et al. New approaches to
#'   population stratification in genome-wide association studies. Nat Rev Genet
#'   11, 459–463 (2010). https://doi.org/10.1038/nrg2813
#'
#' @export
#' 
#' @author William Hemstrom
#' 
#' @examples
#' \dontrun{
#' # from a snpRdata object directly
#' x <- stickSNPs
#' sample.meta(x)$phenotype <- sample(c("case", "control"), nsamps(stickSNPs), TRUE)
#' x <- calc_association(x, c("pop.fam", "pop", ".base"), "phenotype",
#'                       method = "armitage")
#' p <- plot_qq(x, "p_armitage_phenotype", c("pop.fam", "pop", ".base"))
#'                       
#'                       
#' # from a data.frame
#' y <- get.snpR.stats(x, c("pop.fam", "pop", ".base"), "association")
#' y <- y$single
#' 
#' 
#' # with facet/subfacet columns:
#' p <- plot_qq(y, "p_armitage_phenotype", c("fam.pop", "pop", ".base"))
#' 
#' 
#' # from raw data with only one faceting column
#' z <- y[y$facet == "pop",]
#' z <- z[,c(1, 5)]
#' colnames(z)[1] <- "pop"
#' p <- plot_qq(z, "p_armitage_phenotype", "pop")
#' }
plot_qq <- function(x, plot_var, facets = NULL, lambda_gc_correction = FALSE){
  #=============sanity checks and prep============
  msg <- character(0)
  
  .o <- .p <- .e <- NULL
  
  # snpRdata
  if(is.snpRdata(x)){
    facets <- .check.snpR.facet.request(x, facets)
    
    if(plot_var %in% colnames(x@stats)){
      stats <- data.table::as.data.table(.get.snpR.stats(x, facets = facets))
    }
    else{
      msg <- c( msg, "No matching statistics in provided snpRdata object.\n")
    }
  }
  
  # data frame
  else if(is.data.frame(x)){
    if(!data.table::is.data.table(x)){x <- data.table::as.data.table(x)}
    stats <- x
    if(!plot_var %in% colnames(stats)){
      msg <- c(msg, "No matching statistics in provided data.frame.\n")
    }
    
    if(is.null(facets)){
      facets <- ".base"
      stats$facet <- ".base"
      stats$subfacet <- ".base"
    }
    else{
      if(all(c("facet", "subfacet") %in% colnames(stats))){
        bad <- which(!facets %in% stats$facet)
        
        
        if(length(bad) > 0){
          msg <- c(msg, paste0("No  data for facets ", paste0(facets[bad], collapse = ", "), " found in provided data.\n"))
        }
      }
      
      else{
        check <- .split.facet(facets)
        bad <- character(0)
        for(i in 1:length(check)){
          bad <- c(bad, check[[i]][which(!check[[i]] %in% colnames(stats))])
        }
        if(length(bad) > 0){
          msg <- c(msg, paste0("No columns matching facet part(s) ", paste0(check[[i]][bad], collapse = ", "), " found in provided data.\n"))
        }
      }
    }
  }
  
  # otherwise error
  else{
    msg <- c(msg, "x must be a data.frame or snpRdata object.\n")
  }
  
  
  if(length(msg) > 0){
    stop(msg)
  }
  #==============plot===========================
  # prep p values
  rm_rows <- which(is.na(stats[[plot_var]]) | is.nan(stats[[plot_var]]) |
                     stats[[plot_var]] < 0 | stats[[plot_var]] > 1 |
                     is.nan(stats[[plot_var]]))
  if(length(rm_rows) > 0){
    stats <- stats[-rm_rows,]
  }
  
  colnames(stats)[which(colnames(stats) == plot_var)] <- ".p"
  
  out <- vector("list", length(facets))
  names(out) <- facets
  
  for(i in 1:length(facets)){
    # grab and id facets
    split.facet.part <- .split.facet(facets[i])[[1]]
    if(all(c("facet", "subfacet") %in% colnames(stats))){
      tstats <- stats[which(stats$facet == facets[i]),]
      ncols <- unlist(.split.facet(tstats$subfacet))
      ncols <- data.table::as.data.table(matrix(ncols, ncol = length(split.facet.part), byrow = TRUE))
      colnames(ncols) <- split.facet.part
      tstats <- cbind(ncols, tstats)
    }
    else{
      tstats <- stats
    }
    
    # lambda correction
    if(lambda_gc_correction){
      .q <- .lam <- NULL
      
      tstats[,.q := stats::qchisq(1-.p,1)]
      tstats[,.lam := stats::median(.q)/stats::qchisq(0.5,1), by = c(split.facet.part)]
      tstats[,.p := stats::pchisq(.q/.lam, 1, lower.tail = FALSE)]
      
    }
    
    # get the x axis
    tstats[,.o := -log10(sort(.p)), by = c(split.facet.part)]
    tstats[,.e := -log10(stats::ppoints(length(.p))), by = c(split.facet.part)]
    
    # fix facet names to avoid quasiquoting stuff
    colnames(tstats)[colnames(tstats) %in% split.facet.part] <- paste0("facet_", 1:length(split.facet.part))
    
    # fix for CRAN
    # for(j in 1:length(split.facet.part)){
    #   assign(paste0("facet_", j), NULL)
    # }
    
    
    # plot
    p <- ggplot2::ggplot(tstats, ggplot2::aes(x = .e, y = .o)) + 
      ggplot2::geom_point() + 
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "blue") +
      ggplot2::theme_bw() +
      ggplot2::theme(strip.background = ggplot2::element_blank()) +
      ggplot2::ylab("Observed log10(p)") +
      ggplot2::xlab("Expected log10(p)")
    
    if(length(split.facet.part) == 1){
      if(split.facet.part != ".base"){
        p <- p + ggplot2::facet_wrap(~facet_1)
      }
    }
    else if(length(split.facet.part) == 2){
      p <- p + ggplot2::facet_grid(facet_1~facet_2)
    }
    else{
      warning("Cannot plot more than two facets simultaniously.\n")
      next
    }
    
    out[[i]] <- p
  }
  
  return(out)
}


#' Create STRUCTURE-like cluster plots
#'
#' Creates ggplot-based stacked bar charts of assignment probabilities (Q) into
#' an arbitrary 'k' number of clusters like those produced by the program
#' STRUCTURE. Files containing prior results can be parsed and plotted, or new
#' results can be generated via a suite of different methods. Many k values and
#' reps can be run at once, and results can be collapsed across reps using
#' CLUMPP.
#'
#' Individual cluster assignment probabilities can be calculated using several
#' different methods: \itemize{\item{snmf: } sNMF (sparse Non-Negative Matrix
#' Factorization). \item{snapclust: } Maximum-likelihood genetic clustering.
#' \item{admixture: } The ADMIXTURE program. Requires a local admixture
#' executable, and thus cannot run on a Windows platform. \item{ structure: }
#' The STRUCTURE program. Requires a local STRUCTURE executable. many additional
#' options are available for STRUCTURE via other arguments.} These methods are
#' not re-implemented in R, instead, this function calls the
#' \code{\link[LEA]{main_sNMF}}, \code{\link[adegenet]{snapclust.choose.k}}, or
#' a local executable for the ADMIXTURE or STRUCTURE program instead. Please
#' cite the references noted in those functions. For snapclust, the "ward"
#' method is used to initialize clusters if one rep is requested, otherwise the
#' clusters are started randomly each rep. Other methods can be used by
#' providing pop.ini as an additional argument as long as only one rep is
#' requested. Note that at this moment, the snapclust method /emph{is not
#' recommended for use} by the adegenet package maintainers.
#'
#' Multiple different runs can be conducted using the 'reps' argument, and the
#' results can be combined for plotting across all of these reps using the
#' clumpp option. This option calls the CLUMPP software package in order to
#' combine proportion population membership across multiple runs via
#' \code{clumppExport} from the \code{pophelper} package. Again, please cite
#' both CLUMPP and pophelper if using this option.
#'
#' Since CLUMPP is run independently for each value of K, cluster identities
#' often "flip" between K values. For example individuals that are grouped into
#' cluster 1 and K = 3 may be grouped into cluster 2 at K = 4. To adjust this,
#' cluster IDs are iteratively adjusted across K values by flipping IDs such
#' that the euclidean distances between clusters at K and K - 1 are minimized.
#' This tends to produce consistent cluster IDs across multiple runs of K.
#'
#' Individuals can be sorted into by membership proportion into different
#' clusters within populations using the qsort option.
#'
#' Since the clustering and CLUMPP processes can be time consuming and outside
#' tools (such as NGSadmix or fastSTRUCTURE) may be preferred, a nested list of
#' Q matrices, sorted by K and then rep or a character string giving a pattern
#' matching saved Q matrix files in the current working directory may provided
#' directly instead of a snpRdata object. Note that the output for this
#' function, if run on a snpRdata object, will return a properly formatted list
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
#' @param k numeric vector, default 2. The k value (number of clusters) for
#'   which to run the clustering/assignment algorithm. Each provided k value
#'   will be run.
#' @param method character, default "snmf". The clustering/assignment method to
#'   run. Options: \itemize{\item{snmf: } sNMF (sparse Non-Negative Matrix
#'   Factorization). See \code{\link[LEA]{main_sNMF}}. \item{snapclust: }
#'   Maximum-likelihood genetic clustering. See
#'   \code{\link[adegenet]{snapclust.choose.k}}. \item{admixture: } The
#'   ADMIXTURE program. Requires a local admixture executable, and thus cannot
#'   run on a Windows platform. \item{ structure: } The STRUCTURE program.
#'   Requires a local STRUCTURE executable. many additional options are
#'   available for STRUCTURE via other arguments.}
#' @param reps numeric, default 1. The number of independent clustering
#'   repetitions to run.
#' @param update_bib character or FALSE, default FALSE. If a file path to an
#'   existing .bib library or to a valid path for a new one, will update or
#'   create a .bib file including any new citations for methods used. Useful
#'   given that this function does not return a snpRdata object, so a
#'   \code{\link{citations}} cannot be used to fetch references.
#' @param iterations numeric or Inf, default 1000. For snapclust, the maximum
#'   number of iterations to run. For STRUCTURE the number of MCMC steps, should
#'   be in the 10,000+ range.
#' @param burnin numeric, default 100. For STRUCTURE, the number of burnin
#'   iterations to do prior to the main run. This should usually be in the
#'   10,000+ range.
#' @param I numeric or NULL, default NULL. For snmf, how many SNPs should be
#'   used to initialize the search? Initializing with a subset of the total SNPs
#'   can radically speed up computation time for large datasets.
#' @param alpha numeric, default 10. If method = "sNMF", determines the
#'   regularization parameter. For small datasets, this can have a large effect,
#'   and should probably be larger than the default. See documentation for
#'   \code{\link[LEA]{main_sNMF}}. If method = "structure", changes the ALPHA
#'   flag, which determines the degree of admixture. If infer_alpha = TRUE,
#'   instead sets the starting point for the alpha inference.
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
#' @param clumpp_path character, default "/usr/bin/CLUMPP.exe". Path to the
#'   clumpp executable, required if clumpp = T.
#' @param clumpp.opt character, default "greedy". Designates the CLUMPP method
#'   to use. Options: \itemize{ \item{fullsearch: } Search all possible
#'   configurations. Slow. \item{greedy: } The standard approach. Slow for large
#'   datasets at high k values. \item{large.k.greedy: } A fast but less accurate
#'   approach. } See CLUMPP documentation for details.
#' @param structure_path character, default "/usr/bin/structure". Path to the
#'   STRUCTURE executable, required if method = "structure".
#' @param admixture_path character, default "/usr/bin/admixture". Path to the
#'   admixture executable, required if method = "admixture".
#' @param admixture_cv numeric, default 5. Fold to use for cross-validation for
#'   admixture, used to determine the optimum k.
#' @param ID character or NULL, default NULL. Designates a column in the sample
#'   metadata containing sample IDs.
#' @param viridis.option character, default "viridis". Viridis color scale
#'   option. See \code{\link[ggplot2]{scale_gradient}} for details.
#' @param alt.palette character or NULL, default NULL. Optional palette of
#'   colors to use instead of the viridis palette.
#' @param t.sizes numeric, default c(12, 12, 12). Text sizes, given as
#'   c(strip.title, axis, axis.ticks).
#' @param separator_thickness numeric, default 1. Thickness of facet level
#'   separator lines. If 0, no separators drawn. Since separators currently
#'   overlap with samples somewhat, this may be desirable.
#' @param separator_color character, default "white". Color of facet level
#'   separator lines.
#' @param no_admix logical, default FALSE. Used if method = "structure". If
#'   TRUE, the NOADMIX flag in STRUCTURE will be set to 1, meaning that no
#'   admixture will be assumed between clusters.
#' @param use_pop_info logical, default FALSE. Used if method = "structure". If
#'   TRUE, the USEPOPINFO flag in STRUCTURE will be set to 1, meaning that
#'   individuals are assumed to come from the populations that they have been
#'   assigned in the facet provided, and the migrant status of individuals and
#'   their parents, grandparents, etc (going back n generations according to the
#'   gens_back argument) will be returned instead of ancestry proportions. The
#'   resulting plot will not be a typical structure plot.
#' @param loc_prior logical, default FALSE. Used if method = "structure". If
#'   TRUE, the LOCPRIOR flag in STRUCTURE will be set to 1. This will place a
#'   strong population prior on samples according to the facet provided. Useful
#'   with weak data when the population info is known to be robust.
#' @param correlated_frequencies logical, default TRUE. Used if method =
#'   "structure". If TRUE, the FREQSCORR flag in STRUCTURE will be set to 1.
#'   This assumes allele frequencies are correlated between populations. Usually
#'   true when populations have some degree of admixture.
#' @param infer_alpha logical, default TRUE. Used if method = "structure". If
#'   TRUE, the INFERALPHA flag in STRUCTURE will be set to 1, allowing the
#'   optimum alpha to be inferred. Large alpha values imply that most
#'   individuals are admixed.
#' @param separate_pop_alphas logical, default FALSE.  Used if method =
#'   "structure". If TRUE, the POPALPHAS flag in STRUCTURE will be set to 1,
#'   allowing populations to have different alpha values. Usually not
#'   recommended.
#' @param infer_lambda logical, default FALSE.  Used if method = "structure". If
#'   TRUE, the INFERLAMBDA flag in STRUCTURE will be set to 1, allowing the
#'   optimum lambda to be inferred. Smaller values imply that most alleles have
#'   either very low or very high frequencies. Not usually recommended.
#' @param infer_pop_specific_lambda  logical, default FALSE.  Used if method =
#'   "structure". If TRUE, the POPSPECIFICLAMBDA flag in STRUCTURE will be set
#'   to 1, allowing for different pops to have different lambda values.
#' @param lambda numeric, default 1. Used if method = "structure". Changes the
#'   LAMBDA flag. Smaller values imply that most alleles have either very low or
#'   very high frequencies. Used if method = "structure". The default works well
#'   in most cases.
#' @param f_prior_mean numeric, default 0.01. Used if method = "structure".
#'   Changes the FPRIORMEAN flag. F values for each k cluster are drawn from a
#'   gamma prior with this mean. The default value tends to perform well for
#'   detecting small amounts of structure, although it can slightly overestimate
#'   k.
#' @param f_prior_sd numeric, default 0.05. Used if method = "structure".
#'   Changes the FPRIORSD flag. F values for each k cluster are drawn from a
#'   gamma prior with this sd. The default value tends to perform well for
#'   detecting small amounts of structure, although it can slightly overestimate
#'   k.
#' @param uniform_alpha_prior logical, default TRUE. Used if method =
#'   "structure". If TRUE, the UNIFPRIORALPHA flag in STRUCTURE will be set to
#'   1, thus using a uniform prior for alpha between 0 and alpha_max. Usually
#'   works well. If FALSE, uses a gamma prior with mean alpha_prior_a *
#'   alpha_prior_b and variance alpha_prior_a*alpha_prior_b^2.
#' @param alpha_max numeric, default 10. Used if method = "structure". Changes
#'   the ALPHAMAX flag. If uniform_alpha_prior is TRUE, alpha will be drawn from
#'   a uniform prior for alpha between 0 and alpha_max.
#' @param alpha_prior_a numeric, default 1. Used if method = "structure".
#'   Changes the ALPHAPRIORA flag. If uniform_alpha_prior is FALSE, alpha will
#'   be drawn from a gamma prior with mean alpha_prior_a * alpha_prior_b and
#'   variance alpha_prior_a*alpha_prior_b^2.
#' @param alpha_prior_b numeric, default 2. Used if method = "structure".
#'   Changes the ALPHAPRIORA flag. If uniform_alpha_prior is FALSE, alpha will
#'   be drawn from a gamma prior with mean alpha_prior_a * alpha_prior_b and
#'   variance alpha_prior_a*alpha_prior_b^2.
#' @param gens_back numeric, default 2. Used if method = "structure". Changes
#'   the GENSBACK flag. If use_pop_info is TRUE, migration probabilities for
#'   individuals will be determined for the individual themselves plus gens_back
#'   generations prior (parents, grandparents, etc).
#' @param mig_prior numeric, default 0.01. Used if method = "structure". Changes
#'   the MIGRPRIOR flag. Changes the prior value for the migration rate
#'   hyperparameter. Values between 0.001 and 0.01 are usually reasonable.
#' @param locprior_init_r numeric, default 1. Used if method = "structure".
#'   Changes the LOCPRIORINIT flag. Sets the initial strength of the location
#'   prior (r). The default is often reasonable.
#' @param locprior_max_r numeric, default 20. Used if method = "structure".
#'   Changes the MAXLOCPRIOR flag. Sets the maximum value of the location prior
#'   (r). The minimum is always 0. The default is usually reasonable.
#' @param alpha_prop_sd numeric, default 0.025. Used if method = "structure".
#'   Changes the ALPHAPROPSD flag. Changes the sd of the Metropolis-Hastings
#'   alpha drawn during update steps.
#' @param start_at_pop_info logical, default FALSE. If TRUE, the STARTATPOPINFO
#'   flag in STRUCTURE will be set to 1, in which case the clusters will start
#'   at populations during the Metropolis-Hastings search rather than randomly.
#'   Requires a defined facet.
#' @param metro_update_freq numeric, default 10. Used if method = "structure".
#'   Changes the METROFREQ flag. Sets the rate at which Metropolis-Hastings
#'   updates are used. If 0, updates are never used.
#' @param seed integer, default sample(100000, 1). Used if method = "structure".
#'   Starting seed for analysis runs. Each additional run (k value or rep) will
#'   use a successive seed.
#' @param strip_col_names string, default NULL. An optional regular expression
#'   indicating a way to process the column names prior to plotting. Parts of
#'   names matching the strings provided will be cut. Useful for when the facet
#'   argument is something like "meta$pop". A vector of strings will strip
#'   multiple patterns.
#' @param cleanup logical, default TRUE. If TRUE, extra files created during
#'   assignment, clumpp, and plot construction will be removed. If FALSE, they
#'   will be left in the working directory.
#' @param ... additional arguments passed to either \code{\link[LEA]{main_sNMF}}
#'   or \code{\link[adegenet]{snapclust.choose.k}}.
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
#'   analyze and visualize population structure. \emph{Molecular ecology
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
#' @examples
#' \dontrun{
#' # basic sNMF, k = 2 and 3
#' plot_structure(stickSNPs, "pop", k = 2:3, clumpp = FALSE)
#'
#' # basic snapclust
#' plot_structure(stickSNPs, "pop", k = 2:3, clumpp = FALSE, method = "snapclust")
#' }
plot_structure <- function(x, facet = NULL, facet.order = NULL, k = 2, method = "snmf", reps = 1, update_bib = FALSE,
                           iterations = 1000, burnin = 100,
                           I = NULL, alpha = 5, qsort = "last", qsort_K = "last", clumpp = TRUE, clumpp_path = "/usr/bin/CLUMPP.exe",
                           clumpp.opt = "greedy", structure_path = "/usr/bin/structure", admixture_path = "/usr/bin/admixture", 
                           admixture_cv = 5, ID = NULL, viridis.option = "viridis",
                           alt.palette = NULL, t.sizes = c(12, 12, 12), separator_thickness = 1, separator_color = "white", 
                           no_admix = FALSE, use_pop_info = FALSE, loc_prior = FALSE, correlated_frequencies = TRUE,
                           infer_alpha = TRUE, separate_pop_alphas = FALSE, infer_lambda = FALSE, 
                           infer_pop_specific_lambda = FALSE, lambda = 1, f_prior_mean = 0.01, f_prior_sd = 0.05,
                           uniform_alpha_prior = TRUE, alpha_max = 10, alpha_prior_a = 1, alpha_prior_b = 2, 
                           gens_back = 2, mig_prior = 0.01, locprior_init_r = 1, locprior_max_r = 20,
                           alpha_prop_sd = 0.025, start_at_pop_info = FALSE, metro_update_freq = 10, seed = sample(100000, 1), 
                           strip_col_names = NULL, cleanup = TRUE, ...){
  
  rnorm <- V1 <- V2 <- popid <- genback <- ancestry_pop <- ancestry_probability <- V4 <- ..bpc <- ..upc <- NULL
  clean_popid <- prob_correctly_assigned <- K <- est_ln_prob <- Percentage <- Cluster <- NULL
  
  
  
  
  #===========sanity checks===================
  kmax <- max(k)
  msg <- character()
  provided_qlist <- FALSE

  if(!is.null(facet.order)){
    if(is.null(facet)){
      msg <- c(msg, "A faceting variable must be provided if a facet.order is.\n")
    }
    if(length(facet.order) != length(unique(facet.order))){
      msg <- c(msg, "facet.order must contain only unique entries, one per unique category in the provided facet.\n")
    }
  }
  
  # check if this is with a snpRdata object or a qlist and do some other checks
  if(!is.snpRdata(x)){
    if(!is.null(facet.order)){
      if(!is.null(facet)){
        cats <- unique(facet)
        
        if(length(facet.order) != length(cats)){
          msg <- c(msg, "The length of facet.order must equal the number of unique categories in the provided facet.\n")
        }
        
        cats <- unique(facet)
        if(!all(cats %in% facet.order)){
          msg <- c(msg, paste0("All categories in the requested facet must be present in facet.order. Missing categories: ",
                               paste0(cats[which(!cats %in% facet.order)], collapse = ", "), "\n"))
        }
      }
    }
    
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
      else if(provided_qlist == "parse"){
        all.files <- list.files(pattern = x)
        
        concerning_extinsions <- which(!grepl(".qopt$", all.files))
        if(length(concerning_extinsions) > 0){
          warning(paste0("Some files do not end in .qopt, and may not be the expected format:\n",
                  paste0(all.files[concerning_extinsions], collapse = "\n\t")))
        }
        
        lev1 <- all.files[1]
        if(!length(lev1) > 0 | is.na(lev1)){
          msg <- c(msg, paste0("No q files matching '", x, "' located.\n"))
        }
        else{
          if(method == "structure"){
            if(!use_pop_info){
              lev1 <- .readQStructure(lev1)[[1]]
              good <- nrow(lev1) != length(facet)
            }
            else{
              lev1 <- readLines(lev1, n = 16)
              good <- as.numeric(gsub(" ", "", gsub("individuals", "", lev1[16]))) == length(facet)
            }
          }
          else{
            lev1 <- .readQ(lev1)[[1]]
            good <- nrow(lev1) == length(facet)
          }
          if(!good){
            msg <- c(msg, "If a pattern for q files is provided alongside facet information, facet must be a character vector containing population identifiers for each sample in the q files.\n")
          }
          rm(lev1)
        }
      }
      sample_meta <- data.frame(d = facet, stringsAsFactors = F)
      facet <- deparse(substitute(facet))
      colnames(sample_meta) <- facet
    }
  }


  if(provided_qlist == FALSE & reps == "all"){
    msg <- c(msg, "reps = 'all' uninterpretable if a qlist is not provided.\n")
  }
  if(clumpp & reps == 1){
    clumpp <- FALSE
    warning("Since only one rep is requested, clumpp will not be run.\n")
  }
  if(method == "structure" & use_pop_info){
    warning("CLUMPP cannot be used with the use_pop_info option.\n")
    clumpp <- FALSE
  }

  if(clumpp){
    good.clumpp.opts <- c("fullsearch", "greedy", "large.k.greedy")
    clumpp.opt <- tolower(clumpp.opt)
    if(!clumpp.opt %in% good.clumpp.opts){
      msg <- c(msg, paste0("Unaccepted clumpp option. Accepted options: ", paste0(good.clumpp.opts, collapse = ", "), "\n"))
    }
    
    if(!file.exists(clumpp_path)){
      msg <- c(msg, paste0("Could not find clumpp executable at given clumpp_path.\n"))
    }
  }

  # checks for snpRdata objects only
  if(isFALSE(provided_qlist)){
    x <- .add.facets.snpR.data(x, facet)
    good.methods <- c("snapclust", "snmf", "admixture", "structure")
    if(!method %in% good.methods){
      msg <- c(msg, paste0("Unaccepted clustering method. Accepted options: ", paste0(good.methods, collapse = ", "), "\n"))
    }
    
    if(method == "snmf"){
      .check.installed("LEA", "bioconductor")
    }
    if(method == "snapclust"){
      .check.installed("adegenet")
    }
    if(method == "admixture"){
      if(!file.exists(admixture_path)){
        msg <- c(msg, "No file found at provided admixture path.\n")
      }
      if(Sys.info()[1] == "Windows"){
        msg <- c(msg, "Unfortunately, ADMIXTURE is not available for a Windows environment. Please use a unix based environment or pick another assignment approach.\n")
      }
    }
    if(method == "structure"){
      if(!file.exists(structure_path)){
        msg <- c(msg, "No file found at provided structure path.\n")
      }
      
      if(use_pop_info & is.null(facet)){
        stop("Cannot use population info if a facet is not provided.\n")
      }
      
      if(iterations <= 1){
        msg <- c(msg, "Cannot have one or fewer iterations.\n")
      }
    }

    if(length(facet) > 1){
      msg <- c(msg, "Only one facet may be plotted at once.\n")
    }
    if(!is.null(facet[[1]])){
      if(facet[[1]] == ".base"){
        facet <- NULL
      }
      else{
        fcheck <- .check.snpR.facet.request(x, facet, remove.type = "none", return.type = T)
        if(any(fcheck[[2]] != "sample")){
          stop("Only simple, sample level facets allowed.\n")
        }
        facets <- .check.snpR.facet.request(x, facet, remove.type = "snp")
      }
    }
    if(!is.null(facet.order)){
      cats <- .get.task.list(x, facet)
      num.cats <- length(unique(cats[,2]))
      if(num.cats != length(unique(facet.order))){
        msg <- c(msg, "The number of categories provided in facet.order must equal the number of categories in the provided facet.\n")
        if(!all(cats %in% facet.order)){
          msg <- c(msg, paste0("All categories in the requested facet must be present in facet.order. Missing categories: ",
                               paste0(cats[which(!cats %in% facet.order)], collapse = ", "), "\n"))
        }
      }
    }
    sample_meta <- x@sample.meta
    if(!is.null(facet[[1]])){
      sf <- unlist(.split.facet(facet))
      if(length(sf) > 1){
        sample_meta <- cbind(sample_meta, .paste.by.facet(sample_meta, sf))
        colnames(sample_meta)[ncol(sample_meta)] <- facet
      }
    }
  }

  # palette checks
  if(!is.null(alt.palette[1])){

    # is it long enough?
    if(length(alt.palette) < max(k)){
      msg <- c(msg, "Provided alternative palette must contain at least as many colors
               as the requested k values.\n")
    }

    # is everything a valid color (can ggplot work with it)?
    else{
      alt.palette <- alt.palette[1:max(k)]
      tpd <- matrix(stats::rnorm(length(alt.palette)*2, 0, 1), length(alt.palette), 2)
      tpd <- cbind(as.data.frame(tpd), col = 1:max(k))

      color.check <- ggplot2::ggplot(tpd, ggplot2::aes(V1, V2, color = as.factor(col))) +
        ggplot2::geom_point() + ggplot2::scale_color_manual(values = alt.palette)
      tempfile <- tempfile()
      tempfile <- paste0(tempfile, ".pdf")

      suppressMessages(res <- try(ggplot2::ggsave(tempfile, color.check), silent = T))
      .make_it_quiet(file.remove(tempfile))
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

    if(length(x) == 1){
      return(x)
    }
    
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
        for(tk in 1:ncol(x[[i-1]])){
          #save euclidean dist
          elist[tk] <- sum((x[[i]][,j] - x[[i-1]][,tk])^2, na.rm = T)
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
    lx <- as.data.frame(lx)
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
    rmfiles <- list.files(pattern = "qopt")
    file.remove(rmfiles)
    # make a directory and write the k files to it
    for(i in 1:length(x)){
      for(j in 1:reps){
        utils::write.table(x[[i]][[j]], paste0("K", ncol(x[[i]][[j]]), "r", j, "qopt"), sep = " ", quote = F, col.names = F, row.names = F)
      }
    }
  }

  # run clumpp on a directory of q files using pophelper. Import the results into a list of  processed q tables
  run_clumpp <- function(pattern = "qopt"){
    # prepare files and run clumpp
    qfiles <- list.files(full.names = T, pattern = pattern)
    qlist <- .readQ(qfiles)
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
    qlist <- qlist[which(ks %in% k)]
    ## run
    .clumppExport(qlist, parammode = clumpp.opt, exportpath = getwd())
    dirs <- list.files(".", "pop")
    run_dir <- logical(length(dirs))
    for(i in 1:length(k)){
      run_dir[grepl(paste0("_K", k[i]), dirs)] <- TRUE
    }
    dirs <- dirs[run_dir]
    
    for(i in 1:length(dirs)){
      file.copy(clumpp_path, paste0("./", dirs[i], "/"))
      setwd(dirs[i])
      system(basename(clumpp_path))
      setwd("..")
    }
    .collectClumppOutput(filetype = "both", runsdir = getwd(), newdir = "pop-both")


    # import results
    mq <- .readQ(list.files("pop-both/", full.names = T, pattern = "merged"))

    # get only the correct k values (in case this is re-running in a previous directory)
    ks <- unlist(lapply(mq, function(x) attr(x, "k")))
    mq <- mq[which(ks %in% k)]
    if(length(mq) == 0){
      stop("No provided q files have k values in requested range.\n")
    }

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
    if(!isFALSE(qsort)){
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
  
  # parse qfiles with usepopinfo. Terrible to parse in a vectorized way.......
  parse_qfiles_usepopinfo <- function(pattern){
    proc_one_qfile <- function(qf){
      # read and convert to data.table
      tqf <- readLines(qf)
      tqf <- tqf[(grep("Probability of being from assumed population", tqf) + 2)
                 :(grep("Estimated Allele Frequencies in each cluster", tqf) - 3)]
      tqf <- data.table::fread(text = tqf, sep = " ", fill = TRUE, stringsAsFactors = F)
      
      # prepare output
      pop_cols <- grep("Pop", unlist(tqf[1,,]))
      if(length(pop_cols) + 1 != kmax){
        warning(paste0("File ", qf, "matches pattern, but is for k = ", length(pop_cols) + 1, " not k = ", kmax, " as expected, and has been skipped.\n"))
        return(FALSE)
      }
      gens <- pop_cols[2] - pop_cols[1] - 3
      tqf_clean <- data.table(ID = tqf[[2]], missing = as.numeric(gsub("\\(", "", gsub("\\)", "", tqf[[3]])))/100, 
                              popid = tqf[[4]], prob_correctly_assigned = tqf[[6]], sig_mismatch = tqf[[ncol(tqf)]],
                              ancestry_probability = NA, clean_popid = sample_meta[,facet],
                              ancestry_pop = rep(unique(tqf[[4]]), each = nrow(tqf)))
      tqf_clean <- data.table(tqf_clean, genback = rep(1:gens, each = nrow(tqf_clean)))
      tqf_clean$ancestry_probability <- as.numeric(tqf_clean$ancestry_probability)
      
      # assign ancestry probs for inds in each pop
      for(i in 1:(length(pop_cols) + 1)){
        for(j in 1:gens){
          # pops before i
          if(i != 1){
            bpc <- pop_cols[i - 1] + j + 1
            suppressWarnings(tqf_clean[popid < i & genback == j & ancestry_pop == i, ancestry_probability := unlist(tqf[V4 < i, ..bpc])])
          }

          # this pop is left NA
          
          # pops after i
          if(i != length(pop_cols) + 1){
            bpc <- pop_cols[i] + j + 1
            suppressWarnings(tqf_clean[popid > i & genback == j & ancestry_pop == i, ancestry_probability := unlist(tqf[V4 > i, ..bpc])])
          }
        }
      }
      
      return(tqf_clean)
    }
    
    qfiles <- list.files(full.names = T, pattern = pattern)
    cats <- expand.grid(pop = k, rep = 1:reps, generation_back = 0:gens_back)
    cats <- .paste.by.facet(cats, 1:3, "_")
    qlist <- vector("list", length(cats))
    names(qlist) <- cats
    for(i in 1:length(qfiles)){
      qlist[[i]] <- proc_one_qfile(qfiles[i])
      if(is.data.table(qlist[[i]])){
        qlist[[i]][, rep := i]
      }
    }
    
    
    #=========short-circut whole process straight to plot, since all of the lumping, clumpping, and sorting aren't relevant.======
    # clean results
    empties <- unlist(lapply(qlist, isFALSE))
    if(any(empties)){
      qlist <- qlist[-which(empties)]
    }
    qlist <- purrr::compact(qlist)
    found_reps <- length(qlist)
    if(found_reps != reps){
      warning(paste0("Found data for ", found_reps, " reps, but parsing requested for only ", reps, "reps. Check source files matching provided pattern -- Returning data for all found reps.\n"))
    }
    
    qlist <- data.table::rbindlist(qlist)

    # plot
    qlist[, genback := paste0("generation ", qlist$genback - 1)]
    if(!is.null(facet)){
      upc <- c(3, 7)
      upc <- .fix..call(unique(qlist[,..upc]))
      qlist[, ancestry_pop := upc$clean_popid[match(qlist$ancestry_pop, upc$popid)]]
      qlist[, popid := clean_popid]
    }
    
    if(found_reps > 1){
      p <- ggplot2::ggplot(qlist, ggplot2::aes(x = ancestry_pop, y = ancestry_probability, color = prob_correctly_assigned, shape = as.factor(rep)))
      
    }
    else{
      p <- ggplot2::ggplot(qlist, ggplot2::aes(x = ancestry_pop, y = ancestry_probability, color = prob_correctly_assigned))
    }
    
    p <- p + ggplot2::geom_point() +
      ggplot2::facet_grid(genback ~ popid) +
      ggplot2::theme_bw() +
      ggplot2::theme(strip.background = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = t.sizes[3]),
                     strip.text = ggplot2::element_text(size = t.sizes[1]),
                     axis.title = ggplot2::element_text(size = t.sizes[2])) +
      ggplot2::scale_color_viridis_c(option = viridis.option, end = .75) +
      ggplot2::ylab("Migrant Probability") +
      ggplot2::xlab("Source Population")
    
    return(list(plot = p, qdata = qlist))

  }

  # function to parse in q files, used only if clumpp not run and a file pattern is provided
  parse_qfiles <- function(pattern){
    # read in the files
    qfiles <- list.files(full.names = T, pattern = pattern)
    if(method == "structure"){
      if(!use_pop_info){
        qlist <- .readQStructure(qfiles)
      }
      else{
        return(parse_qfiles_usepopinfo(pattern))
      }
    }
    else{
      qlist <- .readQ(qfiles)
    }

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

  # function to sort each element of a qlist by a facet (population, etc). Needed because they must
  # be ordered properly for qsorting to work correctly.
  sort_by_pop_qfiles <- function(qlist, meta, pop.ord){
    ord <- factor(meta[,facet], levels = pop.ord)
    ord <- order(ord)

    for(i in 1:length(qlist)){
      if(is.data.frame(qlist[[i]])){
        qlist[[i]] <- qlist[[i]][ord,, drop = F]
      }
      else{
        for(j in 1:length(qlist[[i]])){
          qlist[[i]][[j]] <- qlist[[i]][[j]][ord,, drop = F]
        }
      }
    }

    nmeta <- as.data.frame(meta[ord,], stringsAsFactors = F)
    colnames(nmeta) <- colnames(meta)

    return(list(qlist = qlist, meta = nmeta))
  }

  #===========run the assignment/clustering method===============
  # each method should return a list of q tables, unprocessed, and possibly work on a K_plot. The list should be nested, k then r.
  if(isFALSE(provided_qlist)){
    if(method == "snapclust"){
      warning("The adegenet maintainers do not currently recommend using snapclust. Consider picking a different assignment method.\n")
      if(!identical(k, 2:kmax)){
        k <- 2:kmax
        warning("Snapclust will always run all values of k between two and the maximum provided k.\n")
      }
      
      if(any(k == 1)){
        k <- k[-which(k == 1)]
        if(length(k) == 0){
          stop("Please supply k values other than k = 1 for snapclust.\n")
        }
      }
      # format and run snapclust
      .make_it_quiet(x.g <- format_snps(x, "adegenet"))

      # initialize
      qlist <- vector("list", length = length(k))
      names(qlist) <- paste0("K_", k)
      for(i in 1:length(qlist)){
        qlist[[i]] <- vector("list", length = reps)
        names(qlist[[i]]) <- paste0("r_", 1:reps)
        attr(qlist[[i]], "k") <- k[i]
      }
      K_plot <- vector("list", length(rep))

      # run once per rep
      for(i in 1:reps){
        if(reps == 1){
          res <- adegenet::snapclust.choose.k(x = x.g, max = kmax, IC.only = FALSE, ...)
        }
        else{
          res <- adegenet::snapclust.choose.k(x = x.g, max = kmax, IC.only = FALSE, pop.ini = NULL, ...)
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
        K_plot[[i]] <- data.frame(val = res[[1]], K = 1:kmax, rep = i, stringsAsFactors = F)
      }

      # fix K plot
      K_plot <- dplyr::bind_rows(K_plot)
      colnames(K_plot)[1] <- names(res)[1]
    }
    else if(method == "snmf"){

      if(any(k == 1)){
        k <- k[-which(k == 1)]
        if(length(k) == 0){
          stop("Please supply k values other than k = 1 for sNMF.\n")
        }
      }
      
      
      # format data
      if(file.exists("lea_input.geno")){
        file.remove("lea_input.geno")
      }
      .make_it_quiet(format_snps(x, "lea", outfile = "lea_input.geno"))

      # run the analysis
      if(!is.null(I)){
        snmf.res <- LEA::snmf("lea_input.geno", K = k, repetitions = reps, iterations = iterations, entropy = T, project = "new", alpha = alpha, I = I, ...)
      }
      else{
        snmf.res <- LEA::snmf("lea_input.geno", K = k, repetitions = reps, iterations = iterations, entropy = T, project = "new", alpha = alpha, ...)
      }

      # read in
      qlist <- vector("list", length(k))
      names(qlist) <-  paste0("K_", k)
      K_plot <- expand.grid(k = k, rep = 1:reps)
      K_plot$Cross.Entropy <- NA
      prog <- 1
      for(i in k){

        # process each rep
        qlist[[prog]] <- vector("list", reps)
        attr(qlist[[prog]], "k") <- i
        names(qlist[[prog]]) <-  paste0("r_", 1:reps)
        for(j in 1:reps){
          qlist[[prog]][[j]] <- as.data.frame(LEA::Q(snmf.res, i, j))
        }
        K_plot[which(K_plot$k == i),]$Cross.Entropy <- data.frame(Cross.Entropy = LEA::cross.entropy(snmf.res, i), K = i)[,1]
        prog <- prog + 1
      }
    }
    else if(method == "admixture"){
      #=========prep===========
      # write plink files
      old.snp.meta <- snp.meta(x)
      snp.meta(x)$chr <- 0
      format_snps(x, "plink", outfile = "plink_files")
      cv_storage <- expand.grid(k, 1:reps)
      colnames(cv_storage) <- c("K", "rep")
      cv_storage <- dplyr::arrange(cv_storage, K, rep)
      cv_storage$cv_error <- NA
      for(i in k){
        for(j in 1:reps){
          cmd <- paste0(admixture_path, " --cv=", admixture_cv, " plink_files.bed ", i, " | tee plink_files_log", i, "_", j, ".out")
          system(cmd)
          file.rename(paste0("plink_files.", i, ".P"), paste0("plink_files.", i, "_", j, ".P"))
          file.rename(paste0("plink_files.", i, ".Q"), paste0("plink_files.", i, "_", j, ".Q"))
          cv_err <- readLines(paste0("plink_files_log", i, "_", j, ".out"))
          cv_storage[which(cv_storage$K == i & cv_storage$rep == j), 3] <- 
            as.numeric(gsub("^CV.+: ", "", cv_err[grep("CV error ", cv_err)]))
        }
      }
      qlist <- parse_qfiles(".Q")
      K_plot <- cv_storage
    }
    else if(method == "structure"){
      if(!identical(k, min(k):kmax) & length(k) != 1){
        k <- min(k):kmax
      }
      
      if(length(k) != 1 & use_pop_info){
        stop("Cannot use pop info across multiple values of k. Please set k equal to the number of levels of the provided facet.\n")
      }
      else{
        cats <- .get.task.list(x, facet)
        num.cats <- length(unique(cats[,2]))
        if(num.cats != kmax & use_pop_info){
          stop(paste0("k must be set to the same value as the number of levels of the provided facet (", num.cats, ").\n"))
        }
      }
      
      tag <- .rand_strings(1, 10)
      
      osp <- options()$scipen
      options(scipen = 999)

      # write the mainparams file
      mainparams <- c(paste0("#define BURNIN ", burnin),
                      paste0("#define NUMREPS ", iterations),
                      "#define INFILE structure_infile",
                      "#define OUTFILE structure_outfile",
                      paste0("#define NUMINDS ", nsamps(x)),
                      paste0("#define NUMLOCI ", nsnps(x)),
                      "#define PLOIDY 2",
                      "#define MISSING -9",
                      "#define ONEROWPERIND 0",
                      "#define LABEL 1",
                      paste0("#define POPDATA ", ifelse(is.null(facet), 0, 1)),
                      "#define POPFLAG 0",
                      "#define LOCDATA 0",
                      "#define PHENOTYPE 0",
                      "#define EXTRACOLS 0",
                      "#define MARKERNAMES 0",
                      "#define RECESSIVEALLELES 0",
                      "#define MAPDISTANCES 0",
                      "#define PHASED 0",
                      "#define PHASEINFO 0",
                      "#define MARKOVPHASE 0",
                      "#define NOTAMBIGUOUS -999")
      
      # write the extraparams file
      extraparams <- c(paste0("#define NOADMIX ", as.numeric(no_admix)),
                       "#define LINKAGE 0",
                       paste0("#define USEPOPINFO ", as.numeric(use_pop_info)),
                       paste0("#define LOCPRIOR ", as.numeric(loc_prior)),
                       paste0("#define FREQSCORR ", as.numeric(correlated_frequencies)),
                       "#define ONEFST 0",
                       paste0("#define INFERALPHA ", as.numeric(infer_alpha)),
                       paste0("#define POPALPHAS ", as.numeric(separate_pop_alphas)),
                       paste0("#define ALPHA ", alpha),
                       paste0("#define INFERLAMBDA ", as.numeric(infer_lambda)),
                       paste0("#define POPSPECIFICLAMBDA ", as.numeric(infer_pop_specific_lambda)),
                       paste0("#define LAMBDA ", lambda),
                       paste0("#define FPRIORMEAN ", f_prior_mean),
                       paste0("#define FPRIORSD ", f_prior_sd),
                       paste0("#define UNIFPRIORALPHA ", as.numeric(uniform_alpha_prior)),
                       paste0("#define ALPHAMAX ", alpha_max),
                       paste0("#define ALPHAPRIORA ", alpha_prior_a),
                       paste0("#define ALPHAPRIORB ", alpha_prior_b),
                       "#define LOG10RMIN -4.0",
                       "#define LOG10RMAX 1.0",
                       "#define LOG10RPROPSD 0.1",
                       "#define LOG10RSTART -2.0",
                       paste0("#define GENSBACK ", gens_back),
                       paste0("#define MIGRPRIOR ", mig_prior),
                       "#define PFROMPOPFLAGONLY 0",
                       paste0("#define LOCISPOP ", ifelse(is.null(facet), 0, 1)),
                       paste0("#define LOCPRIORINIT ", locprior_init_r),
                       paste0("#define MAXLOCPRIOR ", locprior_max_r),
                       "#define PRINTNET     1",
                       "#define PRINTLAMBDA  1",
                       "#define PRINTQSUM    1",
                       "#define SITEBYSITE   0",
                       "#define PRINTQHAT    0",
                       "#define UPDATEFREQ   0",
                       "#define PRINTLIKES   0",
                       "#define INTERMEDSAVE 0",
                       "#define ECHODATA     0",
                       "#define ANCESTDIST   0",
                       "#define NUMBOXES   1000",
                       "#define ANCESTPINT 0.90",
                       "#define COMPUTEPROB 1",
                       "#define ADMBURNIN  0",
                       paste0("#define ALPHAPROPSD ", alpha_prop_sd),
                       paste0("#define STARTATPOPINFO ", as.logical(start_at_pop_info)),
                       "#define RANDOMIZE 0",
                       paste0("#define SEED ", seed),
                       paste0("#define METROFREQ ", metro_update_freq),
                       paste0("#define REPORTHITRATE 0")
                       )
      
      writeLines(mainparams, "mainparams")
      writeLines(extraparams, "extraparams")
      
      # write the input file
      format_snps(x, "structure", facet, outfile = "structure_infile")

      # prep cv storage
      cv_storage <- expand.grid(K = k, rep = 1:reps)
      cv_storage$est_ln_prob <- 0
      cv_storage$mean_ln_prob <- 0
      cv_storage$var <- 0
      
      # run for each k and rep
      for(j in 1:reps){
        for(i in k){
          # run structure
          outfile <- paste0("structure_outfile_k", i, "_r", j, "_", tag)
          cmd <- paste0(structure_path, " -K ", i, " -o ", outfile, " -D ", seed)
          
          system(cmd)
          
          # grab ln summary data
          suppressWarnings(lndat <- readLines(paste0(outfile, "_f")))
          sline <- grep("Estimated Ln Prob of Data", lndat)
          lndat <- lndat[sline:(sline + 2)]
          lndat <- as.numeric(gsub("^.+= ", "", lndat))
          cv_storage[which(cv_storage$K == i & cv_storage$rep == j),3:5] <- lndat

          
          
          seed <- seed + 1
        }
      }
      
      options(scipen = osp)
      
      
      
      # prep summary data for K plot/evanno
      cv_storage <- as.data.table(cv_storage)
      evanno <-cv_storage[, mean(est_ln_prob), by = K]
      colnames(evanno)[2] <- "mean_est_ln_prob"
      if(kmax >= 3 & reps > 1){
        evanno$lnpK <- NA
        evanno$lnppK <- NA
        evanno$deltaK <- NA
        evanno$sd_est_ln_prob <- cv_storage[, sqrt(stats::var(est_ln_prob)), by = K][[2]]
        evanno$lnpK[-1] <- evanno$mean_est_ln_prob[-1] - evanno$mean_est_ln_prob[-nrow(evanno)]
        evanno$lnppK[-nrow(evanno)] <- abs(evanno$lnpK[-nrow(evanno)] - evanno$lnpK[-1])
        # evanno$deltaK[-c(1, nrow(evanno))] <- abs((evanno$mean_est_ln_prob[-1][-1] - 
        #                                             2*evanno$mean_est_ln_prob[-c(1, nrow(evanno))] +
        #                                             evanno$mean_est_ln_prob[-nrow(evanno)][-(nrow(evanno) - 1)])/
        #                                           evanno$sd_est_ln_prob[-c(1, nrow(evanno))])
        evanno$deltaK[-c(1, nrow(evanno))] <- abs(evanno$lnppK)[-c(1, nrow(evanno))]/evanno$sd_est_ln_prob[-c(1, nrow(evanno))] # no reason to resolve for ln''(K)
        
        infs <- which(is.infinite(evanno$deltaK))
        if(length(infs) > 0){
          evanno$deltaK[infs] <- NA
          warning(paste0("For some values of K (", paste0((k)[infs], collapse = ", "), "), all reps had the same estimated ln(likelihood). Since calculating deltaK involves dividing by the standard deviation of the ln(likelihood) estimates across reps, this will return 'Inf', and have thus been assigned a deltaK of NA."))
        }
        
        K_plot <- list(raw = cv_storage, evanno = evanno)
      }
      else{
        K_plot <- list(raw = cv_storage)
      }
      
      
      
      # read in qlist
      qlist <- parse_qfiles(paste0(tag, "_f"))
      
      # end if doing use_pop_info
      if(use_pop_info){
        return(qlist)
      }
      

     
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
    else{
      # need to remove any k = 1, since that'll mess things up down the line
      k1 <- which(unlist(lapply(qlist, attr, which = "k")) == 1)
      if(length(k1 > 0)){
        qlist[[k1]] <- NULL
      }
    }
    
    for(i in 1:length(cq)){
      qlist[[i]][["clumpp"]] <- cq[[i]]
    }

  }
  else{
    # parse in if provided with a pattern
    if(provided_qlist == "parse"){
      qlist <- parse_qfiles(x)
      if(use_pop_info){
        return(qlist)
      }
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
    if(length(tq) > 1){
      tq <- fix_clust(tq)
    }

    ## qsort
    if(!isFALSE(qsort)){
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
  pdat <- pdat[which(pdat$K %in% k),] # double check that we didn't import extra K values, which is possible with externally run data.
  pdat$K <- as.factor(pdat$K)
  levels(pdat$K) <- paste0("K = ", levels(pdat$K))
  pdat$Cluster <- as.factor(pdat$Cluster)
  
  if(any(pdat$K == "K = 1")){
    pdat <- pdat[-which(pdat$K == "K = 1"),]
  }

  ofacet <- facet
  # strip column names if requested
  if(!is.null(strip_col_names)){
    for(pat in 1:length(strip_col_names)){
      colnames(pdat) <- gsub(strip_col_names[pat], "", colnames(pdat))
      facet <- gsub(strip_col_names[pat], "", facet)
    }
  }
  
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
    if(!is.null(facet.order)){
      pops <- facet.order
    }
    else{
      pops <- unique(sample_meta[,ofacet])
    }
    fmt <- sapply(pops, function(x) sum(sample_meta[,ofacet] == x, na.rm = T))
    names(fmt) <- pops
    fmc <- cumsum(fmt)
    fm <- floor(c(0, fmc[-length(fmc)]) + fmt/2)
    breaks <- levels(pdat$ID)[fm]

    seps <- c(0, fmc) + 0.5
    seps[1] <- -.5
    p <- p +
      ggplot2::scale_x_discrete(labels = unique(pdat[,facet]), breaks = breaks, expand = c(0,0)) +
      ggplot2::geom_vline(xintercept = c(fmc[-length(fmc)]) + 0.5, color = separator_color, linewidth = separator_thickness) +
      ggplot2::xlab(label = facet[1])
  }
  else{
    p <- p + ggplot2::scale_x_discrete(expand = c(0,0))
  }


  
  #=============clean-up and return============
  if(cleanup & isFALSE(provided_qlist)){
    if(method == "snmf"){
      file.remove("lea_input.geno", "lea_input.snmfProject")
      unlink("lea_input.snmf", recursive = T)
      unlink("qfiles", recursive = T, force = T)
    }
    else if(method == "snapclust"){
      unlink("qfiles", recursive = T, force = T)
    }
    else if(method == "structure"){
      file.remove("extraparams", "mainparams", "seed.txt", "structure_infile")
      file.remove(c(list.files(pattern = "structure_outfile_k")))
      unlink("qfiles", recursive = T, force = T)
    }
  }
  

  # cite
  keys <- character(0)
  stats <- character(0)
  details <- character(0)
  if(clumpp == TRUE & reps > 1){
    keys <- c(keys, "Jakobsson2007")
    stats <- c(stats, "CLUMPP")
    details <- c(details, paste0("CLUMPP on ", method, " results."))
  }
  if(method == "snmf"){
    keys <- c(keys, "Frichot2015")
    stats <- c(stats, "sNMF")
    details <- c(details, "sparse Non-Negative Matrix Factorization (sNMF)")
  }
  else if(method == "structure"){
    keys <- c(keys, "Pritchard945")
    stats <- c(stats, "STRUCTURE")
    details <- c(details, "STRUCTURE assignment clustering")
  }
  else if(method == "admixture"){
    keys <- c(keys, "Alexander2011")
    stats <- c(stats, "ADMIXTURE")
    details <- c(details, "ADMIXTURE assignment clustering")
  }
  else if(method == "snapclust"){
    keys <- c(keys, "Beugin2018")
    stats <- c(stats, "snapclust")
    details <- c(details, "snapclust assignment clustering")
  }
  keys <- c(keys, "Francis2017")
  stats <- c(stats, "pophelper")
  details <- c(details, "pophelper R package, used to assist file handling during structure-like plot generation")
  
  
  .yell_citation(keys, stats, details, update_bib)

  
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
#' Generates a 1 or 2 dimensional site frequency spectrum 
#' using the projection methods and folding methods of Marth et al (2004) and
#' Gutenkunst et al (2009). This code is essentially an R re-implementation of
#' the SFS construction methods implemented in the program \emph{dadi} (see
#' Gutenkunst et al (2009)).
#' 
#' @param x snpRdata object, matrix, or numeric vector. If a snpRdata object, 
#'   The SNP metadata should contain "ref" and "anc" data. 
#'   If it does not, the major allele will be assumed to be the ancestral. 
#'   Alternatively, either a 2d site frequency spectra stored in a
#'   matrix, with an additional "pops" attribute containing population IDs, such
#'   as c("POP1", "POP2"), where the first pop is the matrix columns and the
#'   second is the matrix rows, or a 1d site frequency spectra stored as a
#'   numeric vector with a similar pops attribute giving the population name.
#'   These objects can be produced from a dadi input file using
#'   \code{\link{make_SFS}}.
#' @param viridis.option character, default "inferno". Viridis color scale
#'   option. See \code{\link[ggplot2]{scale_gradient}} for details.
#' @param log logical, default TRUE. If TRUE, the number of SNPs in each SFS
#'   cell is log transformed.
#' @param facet character, default NULL. Name of the sample metadata column
#'   which specifies the source population of individuals. For now, allows only
#'   a single simple facet (one column).If NULL, runs the entire dataset.
#' @param pops character, default NULL. A vector of population names of up to
#'   length 2 containing the names of populations for which the an SFS is to be
#'   created. If NULL, runs the entire dataset.
#' @param projection numeric. A vector of sample sizes to project the SFS to, in
#'   \emph{number of gene copies}. Sizes too large will result in a SFS
#'   containing few or no SNPs. Must match the length of the provided pops
#'   vector.
#' @param fold logical, default FALSE. Determines if the SFS should be folded or
#'   left polarized. If FALSE, snp metadata columns named "ref" and "anc"
#'   containing the identity of the derived and ancestral alleles, respectively,
#'   should be present for polarization to be meaningful.
#' @param update_bib character or FALSE, default FALSE. If a file path to an
#'   existing .bib library or to a valid path for a new one, will update or
#'   create a .bib file including any new citations for methods used. Useful
#'   given that this function does not return a snpRdata object, so a
#'   \code{\link{citations}} cannot be used to fetch references. Ignored if a
#'   SFS is provided.
#'
#' @return A ggplot2 plot object of the provided SFS.
#'
#' @export
#' 
#' @examples
#' \dontrun{
#' # folded, 1D
#' plot_sfs(stickSNPs, projection = 20)
#' 
#' # unfolded, 1D, one specific population
#' plot_sfs(stickSNPs, facet = "pop", pops = "ASP", projection = 10, fold = FALSE)
#' 
#' # unfolded, two poplations
#' plot_sfs(stickSNPs, facet = "pop", pops = c("ASP", "CLF"), projection = c(10, 10))
#' 
#' # via a sfs matrix, useful for pulling in spectra from elsewhere
#' sfs <- calc_sfs(stickSNPs, facet = "pop", pops = c("ASP", "CLF"), projection = c(10, 10))
#' plot_sfs(sfs)
#' }
plot_sfs <- function(x = NULL, facet = NULL, viridis.option = "inferno", log = TRUE,
                     pops = NULL, projection = NULL, fold = TRUE, update_bib = FALSE){
  p1 <- p2 <- N <- NULL

  #==================sanity checks=============
  msg <- character()
  if(is.snpRdata(x)){
    x <- calc_sfs(x, facet, pops = pops, projection = projection, fold = fold, update_bib = update_bib)
  }
  else{
    msg <- .sanity_check_sfs(x, 1:2)
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  #=======================================================

  # add column names, row names, and melt
  pops <- attr(x, "pop")
  x <- as.data.frame(x)
  if(length(pops) == 1){
    x <- as.data.frame(t(x))
  }
  colnames(x) <- 0:(ncol(x) - 1)
  x$count <- 0:(nrow(x) - 1)
  x[1,1] <- NA # mask the first entry
  msfs <- reshape2::melt(x, id.vars = "count")
  colnames(msfs) <- c("p1", "p2", "N")
  msfs$p1 <- as.integer(as.character(msfs$p1))
  msfs$p2 <- as.integer(as.character(msfs$p2))



  # plot
  ## 2D
  if(nrow(x) > 1){
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
      ggplot2::xlab(pops[1]) +
      ggplot2::scale_x_continuous(expand = c(0, 0))
  }

  return(p)
}

#' Plot STRUCTURE like results on a map.
#'
#' Plots the mean cluster assignment for each population on a map using the
#' \code{scatterpie} package alongside any additional simple feature objects
#' (see \code{sf} from the \code{sf} package). Assignments must be given in the
#' format provided by \code{\link{plot_structure}}. This function is a wrapper
#' which sources code from github.
#'
#' Currently, this only works for simple, sample specific facets. Coordinates
#' for pie charts should be provided as an \code{sf} object, where
#' one column, named for the facet being plotted, provides the subfacet level
#' names matching those in the assignments. Additional sf objects can be
#' provided, which will also be plotted. Note that there is no need to
#' standardize the CRS across the objects, since each will be transformed to
#' match the sample coordinates.
#'
#' @param assignments Structure like results, parsed in or generated via
#'   \code{\link{plot_structure}}, which generates the needed plot data.
#' @param k numeric. Value of K (number of clusters) to plot.
#' @param facet character. The facet by which data is broken down in the passed
#'   assignments.
#' @param pop_coordinates sf object, see the documentation for \code{sf}
#'   function from the \code{sf} package. sf object containing
#'   points/coordinates for each facet level. Must contain a column of data with
#'   population labels named identically to the provided facet (for example,
#'   named "pop" if "pop" is the provided facet.)
#' @param layers list of \code{ggplot2} layer objects, default NULL. Additional
#'   ggplot layers to be plotted in order, below the pie charts, such as maps
#'   with borders, temperatures, forest cover, etc. As a special note,
#'   \code{link[ggnewscale]{new_scale_fill}} can be used to add additional
#'   \code{fill} aesthetic layers without conflict from the resulting pie
#'   charts. See examples.
#' @param pop_names logical, default T. If true, facet level names will be
#'   displayed on the map.
#' @param viridis.option character, default "viridis". Viridis color scale
#'   option. See \code{\link[ggplot2]{scale_gradient}} for details.
#' @param alt.palette character or NULL, default NULL. Optional palette of colors
#'   to use instead of viridis palette  the pie charts.
#' @param radius_scale numeric 0-1, default 0.05. Scale for pie chart radii as a
#'   proportion of the total map space.
#' @param label_args list, default NULL. Named list of arguments passed to
#'   \code{\link[ggrepel]{geom_label_repel}}. For example, passing
#'   list(max.overlaps = 14) will add the max.overlaps argument to the function
#'   call.
#' @param crop logical, default F. If TRUE, will will crop the plot around the
#'   sample points. If false will show the full extent of the data, often set by
#'   any additional sf objects being plotted.
#' @param scale_bar list or NULL, default \code{list()}. Arguments passed to the
#'   \code{\link[ggspatial]{annotation_scale}} function from \code{ggspatial} to
#'   add a scale to the plot. If NULL, no scale added.
#' @param compass list or NULL, default \code{list(style =
#'   ggspatial::north_arrow_fancy_orienteering, location = "br")}. Arguments
#'   passed to \code{\link[ggspatial]{annotation_north_arrow}} function from
#'   \code{ggspatial} to add a compass to the plot. If NULL, no compass added.
#'   Note that calls to alternative styles, like
#'   \code{\link[ggspatial]{north_arrow_fancy_orienteering}} cannot have
#'   \code{()} after them, like they would if called directly.
#' @param ask logical, default TRUE. Should the function ask for confirmation
#'   before sourcing github code?
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # get an sf of the sampling locations
#' lat_long <- data.frame(SMR = c(44.365931, -121.140420), 
#'                        CLF = c(44.267718, -121.255805), 
#'                        OPL = c(44.485958, -121.298360), 
#'                        ASP = c(43.891693, -121.448360), 
#'                        UPD = c(43.891755, -121.451600), 
#'                        PAL = c(43.714114, -121.272797)) # coords for points
#' lat_long <- t(lat_long)
#' colnames(lat_long) <- c("lat", "long")
#' lat_long <- as.data.frame(lat_long)
#' lat_long$pop <- rownames(lat_long)
#' psf <- sf::st_as_sf(as.data.frame(lat_long), coords = c("long", "lat"))
#' psf <- sf::`st_crs<-`(psf, "EPSG:4326")
#'
#' # get the assignments (STRUCTURE-like results)
#' assignments <- plot_structure(stickSNPs, "pop", alpha = 1) 
#'
#' # get a map of oregon as a background from the maps package. 
#' background <- rnaturalearth::ne_states(iso_a2 = "US", returnclass = "sp")
#' background <- sf::st_as_sf(background)
#' background <- background[background$name %in% "Oregon",]
#'
#' # make the plot
#' plot_structure_map(assignments, k = 2, facet = "pop", pop_coordinates = psf,
#'                    layers = list(ggplot2::geom_sf(data = background, 
#'                                      ggplot2::aes(fill = "example")),
#'                    ggnewscale::new_scale_fill()), radius_scale = .2)
#' }
plot_structure_map <- function(assignments, k, facet, pop_coordinates, layers = NULL,
                               pop_names = T, viridis.option = "viridis", alt.palette = NULL,
                               radius_scale = 0.05, label_args = NULL, crop = FALSE,
                               scale_bar = list(), 
                               compass = list(style = ggspatial::north_arrow_fancy_orienteering, location = "br"), 
                               ask = TRUE){
  
  plot_structure_map_extension <- NULL

  if(ask){
    # confirm we want to run this
    cat("plot_structure_map depends on the 'sf' package, which currently fails to compile on a Mac unless gdal is installed.\nThis causes packages with dependecies on it to fail CRAN checks on Mac.\nThis function is a wrapper that sources R scripts from github to pull in the plot_structure_map function.\nIt is tested and should function normally.\nProceed?\t")
    cat("(y or n)\n")
    
    resp <- readLines(n = 1)
    resp <- tolower(resp)
    
    while(resp != "y"){
      if(resp == "n"){
        return(FALSE)
      }
      cat("(y or n)\n")
      resp <- readLines(n = 1)
      resp <- tolower(resp)
    }
  }
  
  .check.installed("sf")
  .check.installed("scatterpie")
  .check.installed("viridis")
  if(!is.null(compass) | !is.null(scale_bar)){.check.installed("ggspatial")}
  
  
  
  # source scripts and pull up internals
  source_file <-  tempfile()
  utils::download.file("https://raw.githubusercontent.com/hemstrow/snpR_extensions/main/plot_structure_map.R", 
                       destfile = source_file)
  source(source_file)
  file.remove(source_file)
  
  
  
  internals <- list(.add.facets.snpR.data = .add.facets.snpR.data,
                 .check.snpR.facet.request = .check.snpR.facet.request, 
                 .split.facet = .split.facet, 
                 .tabulate_genotypes = .tabulate_genotypes, 
                 .make_it_quiet = .make_it_quiet,
                 .fetch.sample.meta.matching.task.list = .fetch.sample.meta.matching.task.list,
                 .fetch.snp.meta.matching.task.list = .fetch.snp.meta.matching.task.list)
  
  mp <- plot_structure_map_extension(assignments = assignments, 
                                     k = k, 
                                     facet = facet, 
                                     pop_coordinates = pop_coordinates,
                                     layers = layers,
                                     pop_names = pop_names, 
                                     viridis.option = viridis.option, 
                                     alt.palette = alt.palette,
                                     radius_scale = radius_scale, 
                                     label_args = label_args, 
                                     crop = crop,
                                     scale_bar = scale_bar,
                                     compass = compass,
                                     internals = internals)
  
  # return
  return(mp)
}



#' Make an heatmap of pairwise Fst values.
#' 
#' Creates a \code{\link[ggplot2]{ggplot}} heatmap of pairwise Fst values
#' previously calculated via \code{\link{calc_pairwise_fst}}. Optionally prints
#' Fst values and marks those with significant p-values.
#' 
#' 
#'@param x snpRdata object.
#'@param facets character, default NULL. Categorical metadata variables by which
#'  to break up plots. Must match facets for which Fst data has been previously
#'  calculated via \code{\link{calc_pairwise_fst}}. 
#'  See \code{\link{Facets_in_snpR}} for more details.
#'@param facet.order character, default NULL. Optional order in which the
#'  levels of the provided facet should appear on the plot, bottom to top/left
#'  to right. If multiple facets are plotted, this must be a named list, named
#'  by facet, otherwise a character vector. See examples.
#'@param viridis.option character, default "inferno". Viridis color scale option
#'  to use. Other color scales may be substituted by appending the
#'  scale_color_continuous and scale_fill_continuous ggplot functions to the
#'  produced plot using the '+' operator.
#'@param print_fst logical, default TRUE. If true, will print Fst values on
#'  plot.
#'@param mark_sig numeric or FALSE, default FALSE. If numeric, Fst values with
#'  bootstrapped p-values below the given number will be marked with an '*'.
#'  Requires that bootstrapped p-values were calculated for the requested 
#'  facets.
#'@param lab_lower logical, default FALSE. If TRUE, labels will be placed in 
#'  the lower triangle instead of on top of the heatmap.
#'  
#'@author William Hemstrom
#'
#'@return A ggplot or list of named ggplots for each facet.
#'@export
#'@examples
#'# Calculate pairwise fst
#'x <- calc_pairwise_fst(stickSNPs, c("pop", "fam"))
#'
#'# plot
#'plot_pairwise_fst_heatmap(x, "pop")
#'
#'# plot without FST values
#'plot_pairwise_fst_heatmap(x, "pop", print_fst = FALSE)
#'
#'# bootstrap some p-values and plot
#'x <- calc_pairwise_fst(x, "pop", boot = 5)
#'## using a high alpha due low number of boots
#'plot_pairwise_fst_heatmap(x, "pop", mark_sig = .2) 
#'
#'# put labels in lower triangle
#'plot_pairwise_fst_heatmap(x, "pop", mark_sig = .2, lab_lower = TRUE) 
#'
#'# provide facet orders
#'plot_pairwise_fst_heatmap(x, c("pop", "fam"), 
#'                          list(pop = c("PAL", "ASP", "UPD", 
#'                                       "CLF", "SMR", "OPL"),
#'                              fam = c("A", "B", "C", "D")),
#'                          print_fst = TRUE, lab_lower = TRUE)
plot_pairwise_fst_heatmap <- function(x, facets = NULL, facet.order = NULL,
                                      viridis.option = "inferno", 
                                      print_fst = TRUE, mark_sig = FALSE,
                                      lab_lower = FALSE){
  subfacet <- sig <- weighted_mean_fst_p <- p1 <- p2 <- weighted_mean_fst <- . <- NULL
  
  #=============sanity checks============
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  msg <- character(0)
  
  facets <- .check.snpR.facet.request(x, facets, return.type = T)
  facets <- facets[[1]][which(facets[[2]] == "sample")]
  if(length(facets) == 0){
    msg <- c(msg, "No sample-level facets requested.\n")
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  #============function=============
  make_one_plot <- function(mean_fst, facet.order = NULL){
    . <- NULL
    mean_fst <- as.data.table(mean_fst)
    mean_fst[,c("p1", "p2") := tstrsplit(subfacet, "~")]
    
    if(!is.null(facet.order)){
      levs <- facet.order
      if(!all(sort(unique(c(mean_fst$p1, mean_fst$p2))) == sort(facet.order))){
        stop(paste0("Subfacets in provided facet.order do not exactly match all of those in the provided data for facet: ", 
                    facets[i], ".\n")) # abuses lexical context, but easy.
      }
      
      i1 <- match(mean_fst$p1, levs)
      i2 <- match(mean_fst$p2, levs)
      flip <- which(i1 > i2)
      if(length(flip) > 0){
        mean_fst[flip, c("p1", "p2") := .(p2, p1)]
      }
      
      mean_fst$p1 <- factor(mean_fst$p1, levs)
      mean_fst$p2 <- factor(mean_fst$p2, levs)
    }
    else{
      levs <- unique(c(mean_fst$p1, mean_fst$p2))
      mean_fst$p1 <- factor(mean_fst$p1, levs)
      mean_fst$p2 <- factor(mean_fst$p2, levs)
    }
    

    if(!isFALSE(mark_sig)){
      if(!any(colnames(mean_fst) == "weighted_mean_fst_p")){
        mark_sig <- FALSE
        warning("Bootstraped significance values not calculated for facet:", mean_fst$facet[1])
      }
      else{
        mean_fst[,sig := ifelse(weighted_mean_fst_p <= mark_sig, "*", "")]
      }
    }
    
    p <- ggplot2::ggplot(mean_fst, ggplot2::aes(x = p1, y = p2, fill = weighted_mean_fst)) +
      ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(option = viridis.option) +
      ggplot2::theme_bw() + ggplot2::theme(axis.title = ggplot2::element_blank()) +
      ggplot2::labs(fill = "Fst")

    if(print_fst){
      if(mark_sig){
        if(lab_lower){
          p <- p +
            ggplot2::scale_x_discrete(drop = FALSE) +
            ggplot2::scale_y_discrete(drop = FALSE)
          p <- p + ggplot2::geom_label(ggplot2::aes(label = paste0(round(weighted_mean_fst, 4), sig), x = p2, y = p1), fill = "white", alpha = .5)
        }
        else{
          p <- p + ggplot2::geom_label(ggplot2::aes(label = paste0(round(weighted_mean_fst, 4), sig)), fill = "white", alpha = .5)
        }
      }
      else{
        if(lab_lower){
          p <- p +
            ggplot2::scale_x_discrete(drop = FALSE) +
            ggplot2::scale_y_discrete(drop = FALSE)
          p <- p + ggplot2::geom_label(ggplot2::aes(label = round(weighted_mean_fst, 4), x = p2, y = p1), fill = "white", alpha = .5)
        }
        else{
          p <- p + ggplot2::geom_label(ggplot2::aes(label = round(weighted_mean_fst, 4)), fill = "white", alpha = .5)
        }
      }
    }
    else if(mark_sig){
      if(lab_lower){
        p <- p +
          ggplot2::scale_x_discrete(drop = FALSE) +
          ggplot2::scale_y_discrete(drop = FALSE)
        p <- p + ggplot2::geom_label(data = mean_fst[which(mean_fst$sig == "*"),], 
                                     ggplot2::aes(label = sig, x = p2, y = p1), fill = "white", alpha = .5)
      }
      else{
        p <- p + ggplot2::geom_label(data = mean_fst[which(mean_fst$sig == "*"),], 
                                   ggplot2::aes(label = sig), fill = "white", alpha = .5)
      }
    }
    
    return(p)
  }

  #============make plots===========
  plots <- vector("list", length(facets))
  names(plots) <- facets
  if(!is.list(facet.order)){
    if(length(facets) > 1){
      stop("If more than one facet is requested and a facet.order is provided, an order for each facet must be included using a named list, see documentation.\n")
    }
    facet.order <- list(facet.order); 
    names(facet.order) <- facets
  }
  else if(length(facet.order) != length(facets)){
    stop("If more than one facet is requested and a facet.order is provided, an order for each facet must be included using a named list with no extra elements, see documentation.\n")
  }
  
  
  for(i in 1:length(facets)){
    plots[[i]] <- make_one_plot(get.snpR.stats(x, facets[i], "fst")$weighted.means, facet.order[[facets[i]]])
  }
  if(length(plots) == 1){plots <- plots[[1]]}
  
  return(plots)
}

#' Basic diagnostic plots
#'
#' Create a suite of basic diagnostic plots (FIS density, 1D SFS, maf density,
#' PCA, % missing data per individual, and he vs ho per SNP)
#' to describe the condition of the data in a snpRdata object.
#'
#' @param x snpRdata object
#' @param facet character, default NULL. Categorical metadata variables by which
#'   to break up plots. Note that only one facet is allowed here. Missingness
#'   and the PCA will have individuals colored by the given sample facet. See
#'   \code{\link{Facets_in_snpR}} for more details.
#' @param plots character vector, default all possible plots except for SFS.
#'   Plot options:
#'   \itemize{\item{fis: } density of FIS scores for all loci within each facet
#'   level.
#'   \item{sfs: } Site Frequency Spectra for the entire dataset.
#'   \item{maf: } density of minor allele frequencies for all loci within each
#'   facet.
#'   \item{pca: } Principal Component Analysis results for the given facet.
#'   \item{missingness: } Proportion of missing alleles across each individual
#'   withing each facet.
#'   \item{heho: } expected vs. observed heterozygosity for each locus within
#'   each facet. Very high expected or observed heterozygosities for many loci
#'   can indicate genotyping issues.}
#' @param projection integer, default floor(nsnps(x)/1.2). A sample size to
#'   project the SFS to, in \emph{number of gene copies}. Sizes too large will
#'   result in a SFS containing few or no SNPs.
#' @param fold_sfs logical, default TRUE. Determines if the SFS should be folded
#'   or left polarized. If FALSE, snp metadata columns named "ref" and "anc"
#'   containing the identity of the derived and ancestral alleles, respectively,
#'   should be present for polarization to be meaningful.
#'
#' @export
#' @author William Hemstrom
#'
#' @return A named list of diagnostic ggplot2 plots.
#' 
#' @examples 
#' \dontrun{
#' # missingness and pca colored by pop
#' plot_diagnostic(stickSNPs, "pop")
#' }
plot_diagnostic <- function(x, facet = NULL, projection = floor(nsnps(x)/1.2), fold_sfs = TRUE,
                            plots = c("fis", "maf", "pca", "missingness", "heho")){
  Individual <- subfacet <- ho <- he <- count <- NULL
  #================checks and init===========
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  facet <- .check.snpR.facet.request(x, facet)
  if(length(facet) > 1){
    facet <- facet[1]
    warning("Only the one facet can be plotted at a time for diagnostic plots. Plotting first provided facet only.\n")
  }
  
  if(facet != ".base"){
    split.facet <- .split.facet(facet)
    facet.vec <- .paste.by.facet(sample.meta(x), facets = unlist(split.facet))
  }
  
  good.plots <- c("fis", "sfs", "maf", "pca", "missingness", "heho")
  plots <- tolower(plots)
  bad.plots <- which(!plots %in% good.plots)
  if(length(bad.plots) > 0){
    stop("Some bad plot types noted:\n\t", paste0(plots[bad.plots], collapse = "\n\t"), "\n\n",
         "Acceptable plot types:", paste0(good.plots, collapse = "\n\t"))
  }
  
  if("heho" %in% plots){
    .check.installed("hexbin")
  }
  
  out <- list()
  
  #=================FIS======================
  if("fis" %in% plots){
    # calc
    calced <- .check_calced_stats(x, facet, "fis")
    if(!unlist(calced)){
      x <- calc_fis(x, facet)
    }
    
    # plot
    fis <- get.snpR.stats(x, facet, "fis")$single
    fis <- ggplot2::ggplot(fis, ggplot2::aes(x = fis, color = subfacet)) + ggplot2::geom_density() +
      ggplot2::theme_bw() +
      ggplot2::scale_color_viridis_d()
    out$fis <- fis
  }
  
  
  #=================plot sfs=================
  if("sfs" %in% plots){
    .make_it_quiet(sfs <- plot_sfs(x, projection = projection, fold = fold_sfs))
    if(fold_sfs){
      sfs <- sfs + ggplot2::xlab("Minor Allele Count")
    }
    else{
      sfs <- sfs + ggplot2::xlab("Derived Allele Count")
    }
    
    out$sfs <- sfs
  }
  
  
  #=================plot maf density=========
  if("maf" %in% plots){
    calced <- .check_calced_stats(x, facet, "maf")
    if(!unlist(calced)){
      x <- calc_maf(x, facet)
    }
    maf <- get.snpR.stats(x, facet, "maf")$single
    maf <- ggplot2::ggplot(maf, ggplot2::aes(x = maf, color = subfacet)) + ggplot2::geom_density() +
      ggplot2::theme_bw() + ggplot2::xlab("Minor Allele Frequency") +
      ggplot2::scale_color_viridis_d()
    
    out$maf <- maf
  }
 
  
  #=================plot pca=================
  if("pca" %in% plots){
    .make_it_quiet(pca <- plot_clusters(x, facets = facet))
    pca <- pca$pca
    
    out$pca <- pca
  }
 
  
  #=================missingness==============
  if("missingness" %in% plots){
    
    miss <- matrixStats::colSums2(ifelse(genotypes(x) == x@mDat, 1, 0))/nsnps(x)
    if(any(miss > 0)){
      if(exists("facet.vec")){
        miss <- ggplot2::ggplot(data.frame(Individual = 1:nsamps(x), miss = miss, facet = facet.vec),
                                ggplot2::aes(x = Individual, y = miss, color = facet.vec)) +
          ggplot2::scale_color_viridis_d() + ggplot2::labs(color = facet) +
          ggplot2::geom_boxplot() +
          ggplot2::geom_point(alpha = .5)
      }
      else{
        miss <- ggplot2::ggplot(data.frame(Individual = 1:nsamps(x), miss = miss, facet = ".base"),
                                ggplot2::aes(y = miss, x = Individual)) +
          ggplot2::geom_boxplot() +
          ggplot2::geom_point() +
          ggplot2::theme()
      }
      
      miss <- miss +
        ggplot2::theme_bw() +
        ggplot2::ylab("Proportion of loci with missing data")
    }
    
    else{
      cat("No missing data.\n")
      miss <- NA
    }
    
    
    out$missingness <- miss
  }
  
  
  #================heho=====================
  if("heho" %in% plots){
    calced <- .check_calced_stats(x, facet, "he")
    if(!unlist(calced)){
      x <- calc_he(x, facet)
    }
    
    calced <- .check_calced_stats(x, facet, "ho")
    if(!unlist(calced)){
      x <- calc_ho(x, facet)
    }
    
    heho <- get.snpR.stats(x, facet, c("he", "ho"))$single
    
    heho <- ggplot2::ggplot(heho, ggplot2::aes(x = ho, y = he)) + 
      ggplot2::geom_hex(ggplot2::aes(fill = ggplot2::after_stat(log(count)))) +
      ggplot2::theme_bw() + ggplot2::xlab("Observed Heterozygosity") +
      ggplot2::ylab("Expected Heterozygosity") +
      ggplot2::geom_abline(slope = 1, intercept = 0) +
      ggplot2::facet_wrap(~subfacet) +
      ggplot2::scale_fill_viridis_c()
    
    out$heho <- heho
  }
  
  return(out)
}




