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
#'@param title character. Plot title.
#'@param t.sizes numeric, default c(16, 13, 10, 12, 10). Text sizes, given as
#'  c(title, legend.title, legend.ticks, axis, axis.ticks).
#'@param background character, default "white". Background color for plot.
#'
#'@return A list containing: \itemize{ \item{plot: } A pairwise LD heatmap as a
#'  ggplot object. \item{dat: } Data used to generate the ggplot object. }
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
#' plot_pairwise_LD_heatmap(dat, c("pop.chr"), "groupIX", c("ASP", "CLF"))
#'
#' # produce plots for every population for linkage group IV
#' plot_pairwise_LD_heatmap(dat, c("pop.chr"), "groupIV")
#' }
#'
plot_pairwise_LD_heatmap <- function(x, facets = NULL, snp.subfacet = NULL, sample.subfacet = NULL, LD_measure = "CLD", r = NULL,
                                     l.text = "CLD", viridis.option = "inferno",
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
    x <- x[!apply(x, 1, function(y)all(is.na(y))), !apply(x, 2, function(y)all(is.na(y)))]

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
#' Approximation and Projection approach implemented in
#' \code{\link[umap]{umap}}. Works by conversion to the "sn" format described in
#' \code{\link{format_snps}} with interpolated missing genotypes.
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
#' runs, and multiple runs are recommended. Uniform Manifold Approximation and
#' Projection (UMAP) coordinates are calculated via \code{\link[umap]{umap}}.
#' UMAP similarly attempts to reduce multi-dimensional results to a two
#' dimensional visualization.
#' 
#' Note that clusters and relative positions of samples from both tSNE and UMAP 
#' may not reliably represent the relationships present in the higher PCA
#' dimensions from which they are created. As such, it is probably not wise to
#' use these methods to draw conclusions about relationships. They are useful
#' exploratory tools, however, and so are kept available here.
#' 
#'
#' For more details on tSNE arguments, \code{\link[Rtsne]{Rtsne}} should be
#' consulted.
#'
#' Additional arguments to the UMAP can be also be provided. Additional
#' information on these arguments can be found in
#' \code{\link[umap]{umap.defaults}}.
#'
#' Data points for individuals can be automatically colored by any sample-level
#' facet categories. Facets should be provided as described in
#' \code{\link{Facets_in_snpR}}. Up to two different sample-level facets can be
#' automatically plotted simultaneously.
#'
#' @param x snpRdata object.
#' @param facets character, default NULL. Categorical sample-level metadata
#'   variables by which to color points. Up to two different sample-specific
#'   facets may be provided. See \code{\link{Facets_in_snpR}} for more details.
#' @param plot_type character, default "pca". c("pca", "tSNE", "umap"). Types of
#'   plots to be produced. Options \itemize{\item{pca: } Principal Component
#'   Analysis, first two dimensions of variance. \item{tSNE: } t-Stochastic
#'   Neighbor Embedding, which collapses dims (see argument) dimensions of
#'   variance into two. \item{umap: } Uniform Manifold Approximation and
#'   Projection, which collapses multiple dimensions of variance into two. } See
#'   description for details.
#' @param check_duplicates logical, default FALSE. Checks for any duplicated
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
#'   minor alleles at missing data points given loci minor allele frequencies.}
#'   \item{iPCA: }{This an iterative PCA approach to interpolate based on
#'   SNP/SNP covariance via \code{\link[missMDA]{imputePCA}}. If the ncp
#'   argument is not defined, the number of components used for interpolation
#'   will be estimated using \code{\link[missMDA]{estim_ncpPCA}}. In this case,
#'   this method is much slower than the other methods, especially for large
#'   datasets. Setting an ncp of 2-5 generally results in reasonable
#'   interpolations without the time constraint.}}
#' @param dims numeric, default 2. Output dimensionality, default 2.
#' @param initial_dims numeric, default 50. The number of dimensions retained in
#'   the initial PCA step during tSNE.
#' @param perplexity numeric, default FALSE. Perplexity parameter, by default
#'   found by \code{\link[mmtsne]{hbeta}}, with beta = 1.
#' @param theta numeric, default 0. Theta parameter from
#'   \code{\link[Rtsne]{Rtsne}}. Default an exhaustive search.
#' @param iter numeric, default 1000. Number of tSNE iterations/umap epochs to
#'   perform.
#' @param viridis.option character, default "viridis". Viridis color scale
#'   option to use for significance lines and SNP labels. See
#'   \code{\link[ggplot2]{scale_gradient}} for details.
#' @param alt.palette charcter or NULL, default NULL. Optional palette of colors
#'   to use instead of the viridis palette.
#' @param ncp numeric or NULL, default NULL. Number of components to consider
#'   for iPCA sn format interpolations of missing data. If null, the optimum
#'   number will be estimated, with the maximum specified by ncp.max. This can
#'   be very slow.
#' @param ncp.max numeric, default 5. Maximum number of components to check for
#'   when determining the optimum number of components to use when interpolating
#'   sn data using the iPCA approach.
#'@param update_bib character or FALSE, default FALSE. If a file path to an
#'   existing .bib library or to a valid path for a new one, will update or
#'   create a .bib file including any new citations for methods used. Useful
#'   given that this function does not return a snpRdata object, so a
#'   \code{\link{citations}} cannot be used to fetch references.
#' @param ... Other arguments, passed to \code{\link[Rtsne]{Rtsne}} or
#'   \code{\link[umap]{umap}}.
#'
#' @return A list containing: \itemize{ \item{data: } Raw PCA, tSNE, and umap
#'   plot data. \item{plots: } ggplot PCA, tSNE, and/or umap plots.} Each of
#'   these two lists may contain one, two, or three objects, one for each PCA,
#'   tSNE, or umap plot requested, named "pca" and "tsne", and "umap",
#'   respectively.
#'
#' @author William Hemstrom
#' @author Matt Thorstensen
#'
#' @references Jesse H. Krijthe (2015). Rtsne: T-Distributed Stochastic Neighbor
#'   Embedding using a Barnes-Hut Implementation, URL:
#'   \url{https://github.com/jkrijthe/Rtsne}.
#' @references Van Der Maaten, L. & Hinton, G. (2008) Visualizing
#'   high-dimensional data using t-SNE. journal of machine learning research.
#'   \emph{Journal of Machine Learning Research}.
#' @references McInnes, L. & Healy (2018). UMAP: uniform manifold approximation
#'   and projection. Preprint at URL: \url{https://arxiv.org/abs/1802.03426}.
#'
#' @seealso \code{\link[mmtsne]{mmtsne}}
#'
#' @export
#'
#' @examples
#' # plot colored by population
#' plot_clusters(stickSNPs, "pop")
#'
#' # plot colored by population and family
#' plot_clusters(stickSNPs, "pop.fam")
plot_clusters <- function(x, facets = NULL, plot_type = "pca", check_duplicates = FALSE,
                          minimum_percent_coverage = FALSE, minimum_genotype_percentage = FALSE, interpolation_method = "bernoulli",
                          dims = 2, initial_dims = 50, perplexity = FALSE, theta = 0, iter = 1000,
                          viridis.option = "viridis", alt.palette = NULL, ncp = NULL, ncp.max = 5, update_bib = FALSE, ...){

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
  good_plot_types <- c("pca", "tsne", "umap")
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
  
  if(interpolation_method == "iPCA"){
    .check.installed("missMDA")
  }
  if(isFALSE(interpolation_method)){
    stop("All methods require no missing data. Please enable interpolation.\n")
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
    sn <- as.data.table(t(sn))
    dups <- which(duplicated(sn) | duplicated(sn, fromLast=TRUE))
    if(length(dups) > 0){
      cat("Duplicates detected, indices:", dups, "\nRemoving all of these!\n")
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

  cat("Plotting using", rm.snps, "loci and", rm.ind, "samples.\n")
  sn <- t(sn)

  #=============prepare plot data=============
  plot_dats <- vector("list", length(plot_type))
  names(plot_dats) <- plot_type

  if("pca" %in% plot_type){
    cat("Preparing pca...\n")
    pca_r <- stats::prcomp(as.matrix(sn))
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
    PC1 <- PC2 <- NULL
    out <- ggplot2::ggplot(tpd, ggplot2::aes(PC1, PC2)) + ggplot2::theme_bw() #initialize plot


    if(facets[1] == ".base"){
      out <- out + ggplot2::geom_point()
    }
    else{
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
    }


    

    if(plot_type[i] == "pca"){
      loadings <- (pca_r$sdev^2)/sum(pca_r$sdev^2)
      loadings <- round(loadings * 100, 2)
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
  
  if(length(keys) > 0){
    .yell_citation(keys, stats, details, update_bib)
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
#' @param sig_below logical, default FALSE. If TRUE, treats values lower than
#'   the significance threshold as significant.
#' @param log.p logical, default FALSE. If TRUE, plot variables and thresholds
#'   will be transformed to -log.
#' @param abs logical, default FALSE. If TRUE, converts the plot variable to
#'   it's absolute value.
#' @param highlight character, numeric, or FALSE, default "significant".
#'   Controls SNP highlighting. If either "significant" or "suggestive", SNPs
#'   above those respetive values will be highlighted. If a numeric vector, SNPs
#'   corresponding to vector entries will be highlighted. See details.
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
#' @param abbreviate_labels numeric or FALSE, default FALSE. If a numeric value,
#'   x-axis chromosome names will be abbreviated using 
#'   \code{\link[base]{abbreviate}}, with each abbreviated label having the
#'   minimum length specified. Helpful when chromosome/scaffold/etc names are
#'   very long.
#'
#' @author William Hemstrom
#' @export
#'
#' @return A list containing \itemize{\item{plot: } A ggplot manhattan plot.
#'   \item{data: } Raw plot data.}
#'
#'
#' @examples
#' \dontrun{
#' # association testing:
#' # add a dummy phenotype and run an association test.
#' x <- stickSNPs
#' sample.meta(x)$phenotype <- sample(c("A", "B"), nsamps(stickSNPs), TRUE)
#' x <- calc_association(x, response = "phenotype", method = "armitage")
#' plot_manhattan(x, "p_armitage_phenotype", chr = "chr", 
#'                log.p = TRUE)
#' 
#' 
#' # other types of stats:
#' # make some data
#' x <- calc_basic_snp_stats(stickSNPs, "pop.chr", sigma = 200, step = 50)
#'
#' # plot pi, breaking apart by population, keeping only the groupIX and
#' # groupIV chromosomes and the ASP, PAL, and SMR populations, with
#' # significant and suggestive lines plotted and SNPs
#' # with pi below the significance level labeled.
#' plot_manhattan(x, "pi", facets = "pop",
#' chr = "chr", chr.subfacet = c("groupIX", "groupIV"),
#' sample.subfacet = c("ASP", "OPL", "SMR"),
#' significant = 0.05, suggestive = 0.15, sig_below = TRUE)
#'
#' # plot FST for the ASP/PAL comparison across all chromosomes,
#' # labeling the first 10 SNPs in x (by row) with their ID
#' plot_manhattan(x, "fst", facets = "pop.chr",
#' sample.subfacet = "ASP~PAL", highlight = 1:20,
#' chr = "chr", snp = ".snp.id")
#'
#' # plot sliding-window FST between ASP and CLF
#' # and between OPL and SMR
#' plot_manhattan(x, "fst", window = TRUE, facets = c("pop.chr"),
#' chr = "chr", sample.subfacet = c("ASP~CLF", "OPL~SMR"),
#' significant = .29, suggestive = .2)
#'
#' # plot using a data.frame,
#' # using log-transformed p-values
#' ## grab data
#' y <- get.snpR.stats(x, "pop", stats = "hwe")$single
#' ## plot
#' plot_manhattan(y, "pHWE", facets = "pop", chr = "chr",
#' significant = 0.0001, suggestive = 0.001,
#' log.p = TRUE, highlight = FALSE)
#' }
plot_manhattan <- function(x, plot_var, window = FALSE, facets = NULL,
                           chr = "chr", bp = "position", snp = NULL,
                           chr.subfacet = NULL, sample.subfacet = NULL,
                           significant = NULL, suggestive = NULL,
                           highlight = "significant",
                           sig_below = FALSE, log.p = FALSE, abs = FALSE,
                           viridis.option = "plasma", viridis.hue = c(.2, 0.5), t.sizes = c(16, 12, 10),
                           colors = c("black", "slategray3"),
                           abbreviate_labels = FALSE){

  #=============sanity checks==============================
  msg <- character()
  pkg.check <- .check.installed("ggrepel")
  if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
  pkg.check <- .check.installed("viridis")
  if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
  
  if(length(msg) > 0){
    stop(msg, collapse = "\n")
  }
  
  #=============grab the desired stats=====================
  #====if a snpRdata object========
  if(class(x) == "snpRdata"){
    if(window == FALSE){
      facets <- .check.snpR.facet.request(x, facets)
    }
    else{
      if(chr != "chr"){
        warning("chr variable will be set to the snp level facet provided to the facets argument for sliding windows.\n")
      }
      if(!is.null(facets)){
        pop.facets <- .check.snpR.facet.request(x, facets, "snp")
        facets <- paste0(pop.facets, ".", chr)
      }
      else{
        facets <- chr
      }

      facets <- .check.snpR.facet.request(x, facets, "none")
    }
    if(plot_var %in% colnames(x@stats)){
      if(window){

        stats <- .get.snpR.stats(x, facets = facets, type = "single.window")
        chr <- "snp.subfacet"

      }
      else{
        stats <- .get.snpR.stats(x, facets = facets)
      }
    }
    else if(plot_var %in% colnames(x@pairwise.stats)){
      if(window){
        if(chr != "chr"){
          warning("chr variable will be set to the snp level facet provided to the facets argument for sliding windows.\n")
        }
        stats <- .get.snpR.stats(x, facets, "pairwise.window")
        chr <- "snp.subfacet"
      }
      else{
        stats <- .get.snpR.stats(x, facets, "pairwise")
      }
    }

    if(nrow(stats) == 0){
      stop("No matching statistics.\n")
    }
  }

  #====otherwise=====
  else if(is.data.frame(x)){
    if(data.table::is.data.table(x)){x <- as.data.frame(x)}
    stats <- x
  }
  else{
    stop("x must be a data.frame or snpRdata object.\n")
  }
  
  #================sanity checks==============
  msg <- character(0)
  
  if(!bp %in% colnames(stats)){
    msg <- c(msg, paste0("Position column: ", bp, " not found in data/snp.meta. Define with argument bp = \n"))
  }
  
  if(!chr %in% colnames(stats)){
    msg <- c(msg, paste0("Chromosome column: ", chr, " not found in data/snp.meta. Define with argument chr = \n"))
  }
  
  if(length(msg) > 0){
    stop(msg)
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
    highlight.label <- NULL
    p <- p + ggrepel::geom_label_repel(data = stats[which(stats$highlight == 1),],
                                       mapping = ggplot2::aes(label = highlight.label), color = add.palette[3],
                                       force = 1.3)
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
  
  return(list(plot = p, data = stats))
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
#' @param facets Character, default NULL. Facets by which to split the plot. See
#'   \code{\link{Facets_in_snpR}}.
#'
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
plot_qq <- function(x, plot_var, facets = NULL){
  #=============sanity checks and prep============
  msg <- character(0)
  
  .o <- .p <- .e <- NULL
  
  # snpRdata
  if(class(x) == "snpRdata"){
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
#' Creates ggplot-based stacked barcharts of assignment probabilities (Q) into
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
#' @param k numeric vector, default 2. The k value (number of clusters) for which
#'   to run the clustering/assignment algorithm. Each provided k value will be run.
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
#' @param alt.palette charcter or NULL, default NULL. Optional palette of colors
#'   to use instead of the viridis palette.
#' @param t.sizes numeric, default c(12, 12, 12). Text sizes, given as
#'   c(strip.title, axis, axis.ticks).
#' @param separator_thickness numeric, default 1. Thickness of facet level
#'   separator lines. If 0, no separators drawn. Since separators currently
#'   overlap with samples somewhat, this may be desirable.
#' @param separator_color character, default "white". Color of facet level
#'   separator lines.
#' @param no_admix logical, default FALSE. Used if method = "structure". If TRUE,
#'   the NOADMIX flag in STRUCTURE will be set to 1, meaning that no admixture
#'   will be assumed between clusters.
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
#'   argument is something like "meta$pop": "^.+\\$" would strip the "meta$pop"
#'   part off. A vector of strings will strip multiple patterns.
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

  .check.installed("pophelper", "github", "royfrancis/pophelper")
  
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
              lev1 <- pophelper::readQStructure(lev1)[[1]]
              good <- nrow(lev1) != length(facet)
            }
            else{
              lev1 <- readLines(lev1, n = 16)
              good <- as.numeric(gsub(" ", "", gsub("individuals", "", lev1[16]))) == length(facet)
            }
          }
          else{
            lev1 <- pophelper::readQ(lev1)[[1]]
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
      fcheck <- .check.snpR.facet.request(x, facet, remove.type = "none", return.type = T)
      if(any(fcheck[[2]] != "sample")){
        stop("Only simple, sample level facets allowed.\n")
      }
      facets <- .check.snpR.facet.request(x, facet, remove.type = "snp")
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
      invisible(utils::capture.output(file.remove(tempfile)))
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
          #save euclidian dist
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
    qlist <- qlist[which(ks %in% k)]
    ## run
    pophelper::clumppExport(qlist, parammode = clumpp.opt, exportpath = getwd())
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
    pophelper::collectClumppOutput(filetype = "both", runsdir = getwd(), newdir = "pop-both")


    # import results
    mq <- pophelper::readQ(list.files("pop-both/", full.names = T, pattern = "merged"))

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
        qlist <- pophelper::readQStructure(qfiles)
      }
      else{
        return(parse_qfiles_usepopinfo(pattern))
      }
    }
    else{
      qlist <- pophelper::readQ(qfiles)
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

  # function to sort each element of a qlist by a facet (population, ect). Needed because they must
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
      invisible(utils::capture.output(x.g <- format_snps(x, "adegenet")))

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
      invisible(utils::capture.output(format_snps(x, "lea", outfile = "lea_input.geno")))

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
      
      tag <- stringi::stri_rand_strings(1, "10")
      
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
                       "#define LOCISPOP 1",
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
      colnames(pdat) <- stringr::str_replace(colnames(pdat), strip_col_names[pat], "")
      facet <- stringr::str_replace(facet, strip_col_names[pat], "")
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
      ggplot2::geom_vline(xintercept = c(fmc[-length(fmc)]) + 0.5, color = separator_color, size = separator_thickness) +
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
  
  
  if(exists("K_plot")){
    return(list(plot = p, data = qlist, plot_data = pdat, K_plot = K_plot))
  }
  else{
    return(list(plot = p, data = qlist, plot_data = pdat))
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
    keys <- c(keys, "Frichot973")
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
#'
#' @param sfs matrix or numeric. Either a 2d site frequency spectra stored in a
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
#' @param x snpRdata object. The SNP metadata must contain "ref" and "anc" data.
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
#'   left polarized.
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
#' # folded, 1D
#' plot_sfs(stickSNPs, projection = 100)
#' 
#' # unfolded, 1D, one specific population
#' plot_sfs(stickSNPs, facet = "pop", pops = "ASP", projection = 40, fold = FALSE)
#' 
#' \dontrun{
#' # unfolded, two poplations
#' plot_sfs(stickSNPs, facet = "pop", pops = c("ASP", "CLF"), projection = c(40, 40))
#' 
#' # via a sfs matrix, useful for pulling in spectra from elsewhere
#' sfs <- calc_sfs(stickSNPs, facet = "pop", pops = c("ASP", "CLF"), projection = c(40, 40))
#' plot_sfs(sfs = sfs)
#' }
plot_sfs <- function(x = NULL, facet = NULL, sfs = NULL, viridis.option = "inferno", log = TRUE,
                     pops = NULL, projection = NULL, fold = TRUE, update_bib = FALSE){
  p1 <- p2 <- N <- NULL
  
  
  if(is.snpRdata(x)){
    sfs <- calc_sfs(x, facet, pops = pops, projection = projection, fold = fold, update_bib = update_bib)
  }
  
  # add column names, row names, and melt
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
      ggplot2::xlab(pops[1]) +
      ggplot2::scale_x_continuous(expand = c(0, 0))
  }

  return(p)
}

#' Plot STRUCTURE like results on a map.
#'
#' Plots the mean cluster assignment for each population on a map using the
#' scatterpie package alongside any additional simple feature objects
#' (\code{\link[sf]{sf}}). Assignments must be given in the format provided by
#' \code{\link{plot_structure}}.
#'
#' Currently, this only works for simple, sample specific facets. Coordinates
#' for pie charts should be provided as an \code{\link[sf]{sf}} object, where
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
#' @param pop_coordinates sf object, see \code{\link[sf]{sf}}. sf object
#'   containing points/coordinates for each facet level. Must contain a column
#'   of data with population labels named identically to the provided facet (for
#'   example, named "pop" if "pop" is the provided facet.)
#' @param sf list of sf objects, default NULL. Additional features to be plotted
#'   alongside points, such as rivers or county lines.
#' @param sf_fill_colors character vector, default "viridis". A vector of colors
#'   to use to fill each polygon sf object. By default, uses the viridis palette
#'   with an alpha of 0.2.
#' @param sf_line_colors character vector, default "viridis". A vector of colors
#'   to use to color lines in in each sf object. By default, uses the viridis
#'   palette with an alpha of 0.2.
#' @param pop_names logical, default T. If true, facet level names will be
#'   displayed on the map.
#' @param viridis.option character, default "viridis". Viridis color scale
#'   option. See \code{\link[ggplot2]{scale_gradient}} for details.
#' @param alt.palette charcter or NULL, default NULL. Optional palette of colors
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
#' @param scale_bar list or NULL, default list(dist = 4, dist_unit = "km",
#'   transform = T). Arguments passed to \code{\link[ggsn]{scalebar}} to add a
#'   scale to the plot. If NULL, no scale added.
#' @param compass list or NULL, list(symbol = 16). Arguments passed to
#'   \code{\link[ggsn]{north}} to add a compass to the plot. If NULL, no compass
#'   added.
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
#' # Note that this map is a bit odd as an sf, but works as an example.
#' background <- maps::map("state", "oregon")
#' background <- sf::st_as_sf(background)
#'
#' # make the plot
#' plot_structure_map(assignments, k = 2, facet = "pop", pop_coordinates = psf, 
#'                    sf = list(background), radius_scale = .2, 
#'                    scale_bar = list(dist = 40, dist_unit = "km", 
#'                                     transform = T), 
#'                    compass = list(symbol = 16, scale = 0.2))
#' }
plot_structure_map <- function(assignments, k, facet, pop_coordinates, sf = NULL, sf_fill_colors = "viridis", sf_line_colors = "viridis",
                               pop_names = T, viridis.option = "viridis", alt.palette = NULL,
                               radius_scale = 0.05, label_args = NULL, crop = FALSE,
                               scale_bar = list(dist = 4, dist_unit = "km", transform = T), compass = list(symbol = 16)){
  
  long <- lat <- pop <- NULL

  #===================sanity checks=================
  msg <- character()
  pkg.check <- .check.installed("scatterpie")
  pkg.check <- c(pkg.check, .check.installed("sf"))
  if(!is.null(compass) | !is.null(scale_bar)){pkg.check <- c(pkg.check, .check.installed("ggsn"))}
  pkg.check <- c(pkg.check, .check.installed("viridis"))
  
  if(is.character(pkg.check)){msg <- c(msg, pkg.check)}
  
  
  if(!is.null(sf)){
    use_crs <- sf::st_crs(pop_coordinates)
    
    # polygon palette
    is.poly <- unlist(lapply(sf, function(x) grepl("POLYGON", sf::st_geometry_type(x)[1])))
    poly.sum <- sum(is.poly)
    if(poly.sum > 0){
      if(sf_fill_colors[1] == "viridis"){
        poly_pal <- viridis::viridis(poly.sum, alpha = .2, option = viridis.option)
      }
      else{
        if(length(sf_fill_colors) != poly.sum){
          warning("The length of the provided fill colors is not the same as the number of polygon sf objects to be plotted, defaulting to viridis.\n")
          poly_pal <- viridis::viridis(poly.sum, alpha = .2, option = viridis.option)
        }
        else{
          poly_pal <- sf_fill_colors
        }
      }
      used_poly_pall <- 0
    }
    
    if(sf_line_colors[1] == "viridis"){
      sf_line_colors <- viridis::viridis(length(sf), option = viridis.option)
    }
    else{
      if(length(sf_line_colors) != length(sf)){
        warning("The length of the provided line colors is not the same as the number of sf objects to be plotted, defaulting to viridis.\n")
        sf_line_colors <- viridis::viridis(length(sf), alpha = .2, option = viridis.option)
      }
    }
  }
  
  K_opts <- unique(assignments$plot_data$K)
  K_opts <- as.numeric(gsub("K = ", "", K_opts))
  if(!k %in% K_opts){
    msg <- c(msg, "Requested value of k not found in provided assignments.\n")
  }
  
  if(length(msg) > 0){stop(msg, "\n")}
  
  #==================prep for plot ===========================
  # generate plotting data.frame
  pie_dat <- as.data.frame(matrix(0, nrow = length(unique(assignments$plot_data[,which(colnames(assignments$plot_data) == facet)])), ncol = 3 + k))
  colnames(pie_dat) <- c("pop", "lat", "long", paste0("Cluster ", 1:k))
  tpd <- assignments$plot_data[assignments$plot_data$K == paste0("K = ", k),]
  tpd <- tpd[,c(facet, "Cluster", "Percentage")]
  tpd$Cluster <- as.numeric(tpd$Cluster)
  anc <- tapply(tpd$Percentage, tpd[,c(facet, "Cluster")], mean)

  ## get the pie coordinates
  if(nrow(pop_coordinates) != nrow(pie_dat)){
    stop(paste0("The number of unique options of ", facet, " (", nrow(pie_dat), ") is not equal to the number of provided coordinates (", nrow(pop_coordinates), ").\n"))
  }
  pie_dat[,1] <- as.data.frame(pop_coordinates)[,facet]
  pie_dat[,3:2] <- sf::st_coordinates(pop_coordinates)
  pie_dat[,4:ncol(pie_dat)] <- anc[match(pie_dat[,1], row.names(anc)),]

  # figure out the radius to use
  lat_range <- range(pie_dat$lat)
  long_range <- range(pie_dat$long)
  r <- min(lat_range[2] - lat_range[1], long_range[2] - long_range[1])
  r <- r*radius_scale


  #============make the plot====================
  mp <- ggplot2::ggplot()

  # add sf overlay if requested.
  if(!is.null(sf)){
    for(i in 1:length(sf)){
      sf[[i]] <- sf::st_transform(sf[[i]], use_crs)
      
      if(is.poly[i]){
        mp <- mp + ggplot2::geom_sf(data = sf[[i]], fill = poly_pal[used_poly_pall + 1], color = sf_line_colors[i])
        used_poly_pall <- used_poly_pall + 1
      }
      else{
        mp <- mp + ggplot2::geom_sf(data = sf[[i]], color = sf_line_colors[i])
      }
    }
  }
  
  # add the scatterpie
  mp <- mp + ggplot2::theme_bw() +
    scatterpie::geom_scatterpie(data = pie_dat, mapping = ggplot2::aes(x = long, y = lat, r = r), cols = colnames(pie_dat)[4:ncol(pie_dat)]) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude")
  
  if(crop){
    xr <- c(min(pie_dat$long - r), max(pie_dat$long) + r)
    yr <- c(min(pie_dat$lat - r), max(pie_dat$lat) + r)
    mp <- mp +
      ggplot2::xlim(xr) +
      ggplot2::ylim(yr)
  }

  if(!is.null(alt.palette)){
    mp <- mp + ggplot2::scale_fill_manual(values = alt.palette)
  }
  else{
    mp <- mp + ggplot2::scale_fill_viridis_d(option = viridis.option)
  }
  if(pop_names){
    # add labels
    #mp + ggplot2::geom_label(data = pie_dat, mapping = ggplot2::aes(x = long, y = lat, label = pop), size = r*label_scale)
    label_call <- c(list(data = pie_dat, mapping = ggplot2::aes(x = long, y = lat, label = pop)), label_args)
    mp <- mp + do.call(ggrepel::geom_label_repel, label_call)
  }
  
  # scale/compass
  if((!is.null(scale_bar) | !is.null(compass))){
    if(!is.null(sf)){ 
      # make up a null sf object to set extent if needed
      dummy <- sf[[1]]
    }
    else{
      dummy <- pop_coordinates
    }
   
    # grab limit info for cropped data...
    b <- ggplot2::ggplot_build(mp)
    lims <- list(x = b$layout$panel_scales_x[[1]]$limits,
                 y = b$layout$panel_scales_y[[1]]$limits)
    
    
    
    if(!is.null(scale_bar)){
      scale_bar$data <- dummy
      if(crop){
        if(!"anchor" %in% names(scale_bar)){
          scale_bar$anchor <- c(x = lims$x[2], y = lims$y[1] + .1*abs(lims$y[1]))
        }
      }
      mp <- mp + do.call(ggsn::scalebar, args = scale_bar)
    }
    
    
    if(!is.null(compass)){
      compass$data <- dummy
      if(crop){
        if(!"anchor" %in% names(compass)){
          compass$anchor <- c(x = lims$x[2] + abs(lims$x[2]*.05), y = lims$y[2] + abs(lims$y[2]*.05))
        }
      }
      mp <- mp + do.call(ggsn::north, args = compass)
    }
  }
  
  # return
  return(mp)
}

#' Plot a phylogenetic-like clustering tree.
#'
#' Generate and plot dendritic trees in the style of a phylogenetic tree for
#' individuals or groups of individuals from snpR data. Note that this function
#' is not overwrite safe.
#'
#' Plots are generated via the \code{\link[ape]{nj}} or \code{\link[ape]{bionj}}
#' ape package for nj or bionj trees. The plots produced are ggplot objects,
#' produced using \code{\link[ggtree]{ggtree}} function, as well as a handful of
#' others from the ggtree package. For more information, see the documentation
#' for those functions and packages.
#'
#' Bootstraps are conducted by re-sampling SNPs with replacement, according to
#' Felsenstein (1985). If no snp level facets are provides, loci are resampled
#' without restraint. If a snp level facet is provided, loci are only resampled
#' within the levels of that facet (e.g. within chromosomes).
#'
#' The genetic distances used to make the trees are calculated using
#' \code{\link{calc_genetic_distances}}. If a sample facet is used, that
#' function uses code derived from \code{\link[adegenet]{adegenet}}. Please cite
#' them and the actual method (e.g. Edwards, A. W. F. (1971)) alongside the
#' tree-building approach.
#'
#' Bootstrapping is done via the boot.phylo function in the ape package, and as such does not
#' support parallel runs on Windows machines.
#'
#' @param x snpRdata object.
#' @param facets character or NULL, default NULL. Facets for which to calculate
#'   genetic distances, as described in \code{\link{Facets_in_snpR}}. If snp or
#'   base level facets are requested, distances will be between individuals.
#'   Otherwise, distances will be between the levels of the sample facets. If
#'   bootstraps are requested, SNPs will
#' @param distance_method character, default "Edwards". Name of the method to
#'   use. Options: \itemize{\item{Edwards} Angular distance as described in
#'   Edwards 1971.} See \code{\link{calc_genetic_distances}}.
#' @param interpolate character, default "bernoulli". Missing data interpolation
#'   method, solely for individual/individual distances. Options detailed in
#'   documentation for \code{\link{format_snps}}.
#' @param tree_method character, default nj. Method by which the tree is
#'   constructed from genetic distances. Options: \itemize{\item{nj}
#'   Neighbor-joining trees, via \code{\link[ape]{nj}}. \item{bionj} BIONJ
#'   trees, according to Gascuel 1997, via \code{\link[ape]{bionj}}.
#'   \item{upgma} UPGMA trees, via \code{\link[stats]{hclust}}.}
#' @param root character or FALSE, default FALSE. A vector containing the
#'   requested roots for each facet. Roots are specified by a string matching
#'   either the individual sample or sample facet level by which to root. If
#'   FALSE for a given facet, trees will be unrooted. Note that all UPGMA trees
#'   are automatically rooted, so this argument is ignored for that tree type.
#' @param boot numeric or FALSE, default FALSE. The number of bootstraps to do
#'   for each facet. See details.
#' @param boot_par numeric or FALSE, default FALSE. If a number, bootstraps will
#'   be processed in parallel using the supplied number of cores.
#' @param update_bib character or FALSE, default FALSE. If a file path to an
#'   existing .bib library or to a valid path for a new one, will update or
#'   create a .bib file including any new citations for methods used. Useful
#'   given that this function does not return a snpRdata object, so a
#'   \code{\link{citations}} cannot be used to fetch references.
#'
#' @author William Hemstrom
#' @export
#'
#' @return A nested, named list containing plots, trees, and bootstraps for each
#'   facet and facet level.
#'
#' @references Felsenstein, J. (1985). Confidence Limits on Phylogenies: An
#'   Approach Using the Bootstrap. Evolution, 39(4), 783–791.
#'   https://doi.org/10.2307/2408678 
#'   
#'   Gascuel, O. (1997).BIONJ: an improved version of the NJ algorithm based on a simple model of
#'   sequence data. Molecular Biology and Evolution, 14(7), 685–695.
#'   https://doi.org/10.1093/oxfordjournals.molbev.a025808
#'
#'   Paradis, E., Claude, J. and Strimmer, K. (2004). APE: analyses of
#'   phylogenetics and evolution in R language. Bioinformatics, 20, 289–290.
#'   
#' 
#' @examples 
#' # Calculate nj trees for the base facet, each chromosome, and for each population.
#' # Note: Examples not run due to odd ape interaction. Work interactively.
#' \dontrun{
#' tp <- plot_tree(stickSNPs, c(".base", "pop", "chr"), 
#'                 root = c(FALSE, "PAL", FALSE))
#' tp$pop$.base$plot # View the resulting plot
#' 
#' # Calculate bionj trees for pop with bootstrapping
#' tp <- plot_tree(stickSNPs, "pop", root = "PAL", boot = 5)
#' tp$pop$.base$plot # nodes now have a bootstrap support indicator
#' }
#' 

plot_tree <- function(x, facets = NULL, distance_method = "Edwards", interpolate = "bernoulli", 
                      tree_method = "nj", root = FALSE,
                      boot = FALSE, boot_par = FALSE, update_bib = FALSE){
  y <- label <- isTip <- NULL
  #=======sanity checks==========
  if(!is.snpRdata(x)){
    stop("x is not a snpRdata object.\n")
  }
  
  facets <- .check.snpR.facet.request(x, facets, remove.type = "none")
  x <- .add.facets.snpR.data(x, facets)
  
  msg <- character(0)
  
  if(!tree_method %in% c("nj", "bionj", "upgma")){
    msg <- c(msg, "Tree method not recognized. Acceptable methods: nj, bionj, or upgma.\n")
  }
  
  if(length(root) < length(facets)){
    if(length(root) == 1){
      root <- rep(root, length(facets))
    }
    else{
      msg <- c(msg, "Root strategy not specified for each facet.\n")
    }
  }
  
  if(.check.installed("ggtree", "github", "YuLab-SMU/ggtree")){
    if(utils::packageVersion("ggtree") < numeric_version("3.1.2")){
      msg <- c(msg, "Package ggtree version 3.1.2+ required. The most recent development version can be installed via remotes::install_github('YuLab-SMU/ggtree')")
    }
  }
  
  if(!isFALSE(boot_par) & Sys.info()["sysname"] == "Windows"){
    msg <- c(msg, "plot_tree uses the ape package parallel setup, which does not work on a Windows platform.\n")
  }
  else if(isFALSE(boot_par)){
    boot_par <- 1
  }
  
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  .check.installed("ape")
  
  #=======function=========
  # expects 
  fun <- function(x, tfacet, distance_method, tree_method, interpolate, root, boot, boot_par){
    if(root == "FALSE"){root <- FALSE}
    #============define tree method===========
    if(tree_method == "nj"){
      if(!isFALSE(root)){
        tree_fun <- function(x, root = NULL, ...) ape::root(ape::nj(x), outgroup = root)
      }
      else{
        tree_fun <- function(x, root = NULL, ...) ape::nj(x)
      }
    }
    else if(tree_method == "bionj"){
      if(!isFALSE(root)){
        tree_fun <- function(x, root = NULL, ...) ape::root(ape::bionj(x), outgroup = root)
      }
      else{
        tree_fun <- function(x, root = NULL, ...) ape::bionj(x)
      }
    }
    else if(tree_method == "upgma"){
      if(isFALSE(root)){
        warning("UPGMA trees are always rooted.")
        root <- TRUE
      }     
      tree_fun <- function(x, root = NULL, ...) ape::as.phylo(stats::hclust(x, "average"))
    }
    
    
    #=============fetch data==============
    #=========raw afms, for non-snp facets===============
    snp.facet <- .check.snpR.facet.request(x, tfacet, "sample", fill_with_base = FALSE, return_base_when_empty = FALSE)
    samp.facet <- .check.snpR.facet.request(x, tfacet, "snp", fill_with_base = FALSE, return_base_when_empty = FALSE)
    if(!is.null(snp.facet) & is.null(samp.facet)){
      opts <- .get.task.list(x, snp.facet)
      amfs <- vector("list", nrow(opts))
      sn <- format_snps(x, "sn", interpolate = interpolate)
      sn <- sn[,-c(1:(ncol(x@snp.meta) - 1))]
      
      for(i in 1:nrow(opts)){
        tsn <- .fetch.snp.meta.matching.task.list(x, opts[i,])
        amfs[[i]] <- t(sn[tsn,])
      }
      names(amfs) <- opts[,4]
      is.snp.only <- TRUE
    }
    else if(tfacet == ".base"){
      sn <- format_snps(x, "sn", interpolate = interpolate)
      sn <- sn[,-c(1:(ncol(x@snp.meta) - 1))]
      amfs <- vector("list", 1)
      amfs[[1]] <- t(sn)
      names(amfs) <- ".base"
      is.snp.only <- TRUE
    }
    else{
      amfs <- get.snpR.stats(x, tfacet, "allele_frequency_matrix")[[1]]
      is.snp.only <- FALSE
    }
    #=============make a tree=============
    if(is.snp.only){
      tdf <- function(part, distance_method = NULL, ...) stats::dist(part)
    }
    else{
      tdf <- function(part, distance_method = NULL, ...) .get_dist(part, distance_method)[[1]]
    }
    
    tree <- lapply(amfs, function(y) tree_fun(tdf(y, distance_method), root))
    
    #=========bootstrap, if requested=============
    if(!isFALSE(boot)){

      #======define function========
      apply_boots <- function(tree, boots, root){
        boot_val <- ape::prop.clades(tree, boots, rooted = !isFALSE(root))
        boot_val[is.na(boot_val)] <- 0
        boot_val <- paste0((boot_val/length(boots))*100, "%") 
        tree$node.label <- boot_val
        return(tree)
      }
      boot_fun <- function(part, ref, tree_fun, distance_method, is.snp.only, root, B){
        boots <- ape::boot.phylo(phy = ref, x = part, FUN = function(y) tree_fun(tdf(y, distance_method), outgroup = root), 
                                 trees = TRUE, mc.cores = boot_par, B = B, block = ifelse(is.snp.only, 1, 2))$trees
        
        res <- apply_boots(ref, boots, root)
        return(res)
      }
      
      
      
      #======run bootstraps=========
      cat("Generating Bootstraps for", tfacet, "...\n")
      for(i in 1:length(tree)){
        cat("Facet part", i, "out of", length(tree), "\n")
        tree[[i]] <- boot_fun(amfs[[i]], tree[[i]], tree_fun, distance_method, is.snp.only, root, boot)
      }
    }
    
    
    #=========generate plot==================
    tout <- vector("list", length(tree))
    names(tout) <- names(tree)
    
    for(i in 1:length(tree)){
      # make the plot
      tp <- ggplot2::ggplot(tree[[i]], ggplot2::aes(x, y), label = label) + ggtree::geom_tree() +
        ggtree::geom_tiplab() + ggtree::theme_tree()
      
      
      # add parts
      ## boot values
      if(!isFALSE(boot)){
        tp <- tp + ggtree::geom_text2(ggplot2::aes(subset = !isTip, label = label),
                                      hjust = -.3)
      }
      
      tout[[i]] <- list(plot = tp, tree = tree[[i]])
    }
    return(tout)
  }
  
  #=======run==============
  needed_dists <- .check_calced_stats(x, facets, paste0("genetic_distance", "_", distance_method, "_", interpolate))
  needed_dists <- which(!unlist(needed_dists))
  if(length(needed_dists) > 0){
    invisible(utils::capture.output(
      x <- calc_genetic_distances(x, facets[needed_dists], distance_method, interpolate = interpolate)))
  }
  
  out <- vector("list", length = length(facets))
  for(i in 1:length(facets)){
    cat("Facet:", facets[i], "\n")
    out[[i]] <- fun(x, facets[i], distance_method, tree_method, interpolate, root[i], boot, boot_par)
  }
  names(out) <- facets
  
  
  # update bib
  keys <- character()
  stats <- character()
  details <- character()
  
  if(distance_method == "Edwards"){
    keys <- c(keys, "Edwards1971")
    stats <- c(stats, "genetic_distance")
    details <- c(details, "Edwards' Angular Genetic Distance")
  }
  else if(distance_method == "Nei"){
    keys <- c(keys, "Nei1978")
    stats <- c(stats, "genetic_distance")
    details <- c(details, "Nei's genetic distance")
  }
  
  if(tree_method == "nj"){
    keys <- c(keys, "Felsenstein1985")
    stats <- c(stats, "tree")
    details <- c(details, "Neighbor-joining tree construction method")
  }
  else if(tree_method == "bionj"){
    keys <- c(keys, "Gascuel1997")
    stats <- c(stats, "tree")
    details <- c(details, "BIONJ tree construction method")
  }
  else if(tree_method == "upgma"){
    keys <- c(keys, "Sokal1958")
    stats <- c(stats, "tree")
    details <- c(details, "UPGMA tree construction method")
  }
  
  keys <- c(keys, "Paradis2004")
  stats <- c(stats, "tree")
  if(!isFALSE(boot)){
    details <- c(details, "tree construction and bootstrapping")
  }
  else{
    details <- c(details, "tree construction")
  }
  
  if(length(keys) > 0){
    .yell_citation(keys, stats, details, update_bib)
  }
  
  # return
  return(out)
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
#'x <- calc_pairwise_fst(stickSNPs, "pop")
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
plot_pairwise_fst_heatmap <- function(x, facets = NULL, 
                                      viridis.option = "inferno", 
                                      print_fst = TRUE, mark_sig = FALSE,
                                      lab_lower = FALSE){
  subfacet <- sig <- weighted_mean_fst_p <- p1 <- p2 <- weighted_mean_fst <- NULL
  
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
  make_one_plot <- function(mean_fst){
    mean_fst <- as.data.table(mean_fst)
    mean_fst[,c("p1", "p2") := tstrsplit(subfacet, "~")]
    levs <- unique(c(mean_fst$p1, mean_fst$p2))
    mean_fst$p1 <- factor(mean_fst$p1, levs)
    mean_fst$p2 <- factor(mean_fst$p2, levs)
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
          p <- p + ggplot2::geom_label(ggplot2::aes(label = paste0(round(weighted_mean_fst, 4), sig), x = p2, y = p1), fill = "white", alpha = .5)
        }
        else{
          p <- p + ggplot2::geom_label(ggplot2::aes(label = paste0(round(weighted_mean_fst, 4), sig)), fill = "white", alpha = .5)
        }
      }
      else{
        if(lab_lower){
          p <- p + ggplot2::geom_label(ggplot2::aes(label = round(weighted_mean_fst, 4), x = p2, y = p1), fill = "white", alpha = .5)
        }
        else{
          p <- p + ggplot2::geom_label(ggplot2::aes(label = round(weighted_mean_fst, 4)), fill = "white", alpha = .5)
        }
      }
    }
    else if(mark_sig){
      if(lab_lower){
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
  for(i in 1:length(facets)){
    plots[[i]] <- make_one_plot(get.snpR.stats(x, facets[i], "fst")$weighted.means)
  }
  if(length(plots) == 1){plots <- plots[[1]]}
  
  return(plots)
}
