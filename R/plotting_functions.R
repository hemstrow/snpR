#Prepares a heatmap from pariwse LD data.
#Inputs:  x: Data, format like that given by LD_full_pairwise rsq and Dprime outputs.
#            However, first column can contain snp names rather than first column of data. Otherwise,
#            Data will be reformated in this way.
#         r: Region of the chromosome to subset and plot.
#            Given in mb in the format numeric vector c(lower, upper).
#         l.text: Legend title.
#         colors: Colors to use, character vector format c("lower", "upper").
#         title: Plot title
#         t.sizes: Text sizes, numeric vector format c(title, legend.title, legend.ticks, axis, axis.ticks)

#'Create a heatmap from pairwise linkage data.
#'
#'\code{LD_pairwise_heatmap} prepares a ggplot2 heatmap from pairwise LD SNP data, in the format of that given by \code{\link{LD_full_pairwise}}. Should also work on pairwise FST data. Darker values are higher. Special thanks to Nicholas Sard for part of this code!
#'
#'Since the output is a ggplot object, options can be added or changed by adding "+ function()" to the end of \code{LD_pairwise_heatmap}. Some common options are also built into this function as agruments, but can be overwritten freely.
#'
#'Output from \code{\link{LD_full_pairwise}} runs directly or after saving and re-import without issue. Position data for the rows can be saved either as rownames or as the first column.
#'
#'This function was partially written by Nicholas Sard.
#'
#' @param x Data, format like that given by LD_full_pairwise rsq and Dprime outputs. Column names must be snp positions, as must either the first column or row names. Alternatively, a named list of the several such objects to plot and their titles can be provided.
#' @param r Region of the chromosome to subset and plot. Given in kb in the format numeric vector c(lower, upper).
#' @param l.text Legend title.
#' @param colors Colors to use in tiles, character vector format c("lower", "upper").
#' @param title Plot title.
#' @param t.sizes Text sizes, numeric vector format c(title, legend.title, legend.ticks, axis, axis.ticks)
#' @return A pairwise LD heatmap as a ggplot object. If multiple objects were provided to plot, these will be included as objects.
#'
#' @examples
#' #base plot
#' LD_pairwise_heatmap(stickLD, title = "ASP group IX Pairwise LD")
#'
#' #Change colors with agruments, add lines at SNP 10 using ggplot functions.
#' LD_pairwise_heatmap(stickLD, colors = c("green", "red"), title = "ASP group IX Pairwise LD") + geom_vline(xintercept = 10, color = "blue") + geom_hline(yintercept = 10, color = "blue")
#'
plot_pairwise_LD_heatmap <- function(x, facets = NULL, snp.subfacet = NULL, sample.subfacet = NULL, LD_measure = "rsq", r = NULL,
                                     l.text = "rsq", viridis.option = "B",
                                     title = NULL, t.sizes = c(16, 13, 10, 12, 10),
                                     background = "white"){
  #Created in part by Nick Sard
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

#'tSNE from genetic data.
#'
#'\code{tSNEfromPA} creates a ggplot object tSNE from presence/absense allelic data using the Barnes-Hut simulation at theta>0 implemented in \code{\link[Rtsne]{Rtsne}}. If individuals which were sequenced at too few loci are to be filtered out, both the mc and counts argument must be provided.
#'
#'See the documentaion for \code{\link[Rtsne]{Rtsne}} for details. Defaults match those of \code{\link[Rtsne]{Rtsne}}. Argument documentation taken from that function.
#'
#'Description of x:
#'    SNP or other allelic data in presence/absence format, as given by \code{\link{format_snps}} option 7. An additional column of population IDs titled "pop" must also be provided.
#'
#'This function results from collaboration with Matt Thorstensen.
#'
#' @param x Input presence/absence data, as described in details.
#' @param ecs Numeric. The number of extra metadata columns at the start of the input. Must be more that two to avoid errors. I really should fix that at some point. Includes the "pop" collumn.
#' @param plot.vars FALSE or character vector, default "pop". If FALSE, a basic plot is produced. Up to two variables to be plotted with names corresponding to column names in x be given as a character vector.
#' @param dims Integer, output dimensionality, default 2.
#' @param initial_dims Integer, default 50. The number of dimensions retained in the initial PCA step.
#' @param perplexity Perplexity parameter, by default found by \code{\link[mmtsne]{hbeta}}, with beta = 1.
#' @param theta Theta parameter from \code{\link[Rtsne]{Rtsne}}. Default 0, an exhaustive search.
#' @param iter Integer, default 5000. Number of tSNE iterations to perform.
#' @param c.dup boolean, default FALSE. Should duplicate individuals be searched for and removed? This is very slow if the data set is large!
#' @param mc Numeric or FALSE, default FALSE. Should poorly sequenced individuals be removed? If so, what is the minimum acceptable count of sequenced loci (as specificed in the vector provided to counts)?
#' @param counts Numeric vector, default FALSE, containing the number of loci sequenced per individual.
#' @param do.plot boolean, default TRUE. Should a plot be produced?
#' @param ... Other arguments, passed to \code{\link[Rtsne]{Rtsne}}.
#'
#' @return A list containing the raw tSNE output and a tSNE plot in the form of a ggplot graphical object. The plot can be changed as usual with ggplot objects.
#'
#' @examples
#' tSNEfromPA(stickPA, 2, c.dup = TRUE)
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
  suppressWarnings(sn <- format_snps(x, "sn", interpolate = interpolation_method))
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


  rm.snps <- ncol(sn)
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

  browser()

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
  for(i in 1:length(plot_dats)){
    browser()
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
    if(length(facets) == 2){
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


