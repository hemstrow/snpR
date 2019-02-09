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
LD_pairwise_heatmap <- function(x, r = NULL,
                                l.text = "rsq", colors = c("white", "black"),
                                title = NULL, t.sizes = c(16, 13, 10, 12, 10),
                                background = "white"){

  #Created in part by Nick Sard

  if(is.matrix(x)){
    x <- as.data.frame(x, stringsAsFactors = F)
  }

  #function to prepare data.
  prep_hm_dat <- function(x, r = NULL){
    #remove columns and rows with no data
    x <- x[!apply(x, 1, function(y)all(is.na(y))), !apply(x, 2, function(y)all(is.na(y)))]

    #if first column doesn't contain positions, take the row names and set as first column.
    if(any(is.na(x[,1]))){
      x <- cbind(V1 = row.names(x), x)
    }


    #melting the df and fixing the names of the columns
    heatmap_x <- reshape2::melt(x, id.vars = "V1")
    names(heatmap_x) <- c("SNPa", "SNPb", "value")

    #getting rid of all the zeros from snps being compared to themselves
    heatmap_x$SNPa <- as.numeric(as.character(heatmap_x$SNPa))
    heatmap_x$SNPb <- as.numeric(as.character(heatmap_x$SNPb))
    heatmap_x <- heatmap_x[!(heatmap_x$SNPa == heatmap_x$SNPb),]

    #removing NA values, since no Dups.  Note, also grabbing only section of interest
    heatmap_x <- heatmap_x[!is.na(heatmap_x$value),]

    #remove any other NAs
    heatmap_x <- heatmap_x[!is.na(heatmap_x$SNPa) & !is.na(heatmap_x$SNPb),]

    #make sure that for all comparisons, SNPa is less than SNPb.
    vio <- which(heatmap_x$SNPa > heatmap_x$SNPb)
    if(length(vio) > 0){
      viofix <- cbind(SNPa = heatmap_x[vio,]$SNPb, SNPb = heatmap_x[vio,]$SNPa, value = heatmap_x[vio,]$value)
      heatmap_x <- rbind(heatmap_x[-vio,], viofix)
    }

    #get site names
    ms <- unique(c(heatmap_x$SNPa, heatmap_x$SNPb))
    ms <- sort(as.numeric(as.character(ms)))

    #subset down to the desired r if requested
    if(!is.null(r)){
      r <- r*1000000
      heatmap_x <- heatmap_x[heatmap_x$SNPa >= r[1] & heatmap_x$SNPa <= r[2] &
                               heatmap_x$SNPb >= r[1] & heatmap_x$SNPb <= r[2],]
      ms <- ms[ms >= r[1] & ms <= r[2]]
    }

    #finish messing with site names
    ms <- ms/1000000
    ms <- sort(ms)
    ms <- as.factor(ms)

    #set site positions to mb, convert to factor
    heatmap_x$SNPa <- heatmap_x$SNPa/1000000
    heatmap_x$SNPb <- heatmap_x$SNPb/1000000
    heatmap_x$SNPa <- as.factor(heatmap_x$SNPa)
    heatmap_x$SNPb <- as.factor(heatmap_x$SNPb)

    #reordering based on factors
    heatmap_x[["SNPa"]]<-factor(heatmap_x[["SNPa"]],levels=ms, ordered = T)
    heatmap_x[["SNPb"]]<-factor(heatmap_x[["SNPb"]],levels=rev(ms), ordered = T)

    return(heatmap_x)
  }


  if(is.data.frame(x)){
    #the plot
    heatmap_x <- prep_hm_dat(x, r)
    out <- ggplot2::ggplot(heatmap_x, ggplot2::aes(x = SNPa, y=SNPb, fill=value))+
      ggplot2::geom_tile(color = "white")+
      ggplot2::scale_fill_gradient(low = colors[1], high = colors[2]) +
      ggplot2::theme_bw()+
      ggplot2::labs(x = "",y="", fill=l.text)+
      ggplot2::theme(legend.title= ggplot2::element_text(size = t.sizes[2]),
            axis.text = ggplot2::element_text(size = t.sizes[5]),
            panel.grid.major = ggplot2::element_line(color = background),
            plot.title = ggplot2::element_text(size = t.sizes[1], hjust = 0.5),
            axis.title = ggplot2::element_text(size = t.sizes[4]),
            legend.text = ggplot2::element_text(size = t.sizes[3]),
            panel.background = ggplot2::element_rect(fill = background, colour = background)) +
      ggplot2::scale_x_discrete(breaks = levels(heatmap_x$SNPa)[c(T, rep(F, 20))], label = abbreviate) +
      ggplot2::scale_y_discrete(breaks = levels(heatmap_x$SNPb)[c(T, rep(F, 20))], label = abbreviate) +
      ggplot2::ylab("Position (Mb)") + ggplot2::xlab("Position (Mb)")

    if(!is.null(title)){
      out <- out + ggplot2::ggtitle(title)
    }
    out <- list(plot = out, dat = heatmap_x)
  }

  else{
    #multple plots
    heatmap_x <- cbind(prep_hm_dat(as.data.frame(x[[1]], stringsAsFactors = F), r), var = names(x)[1])
    for(i in 2:length(x)){
      heatmap_x <- rbind(heatmap_x, cbind(prep_hm_dat(x[[i]], r), var = names(x)[i]))
    }

    #refactor
    heatmap_x$SNPa <- as.factor(as.numeric(as.character(heatmap_x$SNPa)))
    heatmap_x$SNPb <- as.factor(as.numeric(as.character(heatmap_x$SNPb)))

    #plot, with facets
    out <- ggplot2::ggplot(heatmap_x, ggplot2::aes(x = SNPa, y=SNPb, fill=value))+
      ggplot2::geom_tile(color = "white")+
      ggplot2::facet_wrap(~var) +
      ggplot2::scale_fill_gradient(low = colors[1], high = colors[2]) +
      ggplot2::theme_bw()+
      ggplot2::labs(x = "",y="", fill=l.text)+
      ggplot2::theme(legend.title= ggplot2::element_text(size = t.sizes[2]),
                     axis.text = ggplot2::element_text(size = t.sizes[5]),
                     panel.grid.major = ggplot2::element_line(color = background),
                     strip.background = ggplot2::element_blank(),
                     strip.text = ggplot2::element_text(hjust = 0.01, size = t.sizes[1]),
                     axis.title = ggplot2::element_text(size = t.sizes[4]),
                     legend.text = ggplot2::element_text(size = t.sizes[3]),
                     panel.background = ggplot2::element_rect(fill = background, colour = background)) +
      ggplot2::scale_x_discrete(breaks = levels(heatmap_x$SNPa)[c(T, rep(F, 20))], label = abbreviate) +
      ggplot2::scale_y_discrete(breaks = rev(levels(heatmap_x$SNPb)[c(T, rep(F, 20))]), label = abbreviate,
                                limits = rev(levels(heatmap_x$SNPb))) +
      ggplot2::ylab("Position (Mb)") + ggplot2::xlab("Position (Mb)")

    if(!is.null(title)){
      out <- out + ggplot2::ggtitle(title)
    }

    out <- list(plot = out, dat = heatmap_x)
  }
  return(out)
}


#'PCAs from genetic data.
#'
#'\code{PCAfromPA} creates a ggplot object PCA from presence/absense allelic data. If individuals which were sequenced at too few loci are to be filtered out, both the mc and counts argument must be provided.
#'
#'Description of x:
#'    SNP or other allelic data in presence/absence format, as given by \code{\link{format_snps}} option 7. An additional column of population IDs titled "pop" must also be provided.
#'
#' @param x Input presence/absence data, as described in details.
#' @param ecs Numeric. The number of extra metadata columns at the start of the input. Must be more that two to avoid errors. I really should fix that at some point. Includes the "pop" collumn.
#' @param do.plot FALSE or character vector, default "pop". If FALSE, no plot is produced. Up to two variables to be plotted with names corresponding to column names in x be given as a character vector.
#' @param c.dup boolean, default FALSE. Should duplicate individuals be searched for and removed? This is very slow if the data set is large!
#' @param mc Numeric or FALSE, default FALSE. Should poorly sequenced individuals be removed? If so, what is the minimum acceptable count of sequenced loci (as specificed in the vector provided to counts)?
#' @param counts Numeric vector, default FALSE, containing the number of loci sequenced per individual.
#'
#' @return A list containing the raw PCA output and a PCA plot in the form of a ggplot graphical object. The plot can be changed as usual with ggplot objects.
#'
#' @examples
#' PCAfromPA(stickPA, 2)
PCAfromPA <- function(x, ecs, do.plot = "pop", c.dup = FALSE, mc = FALSE, counts = FALSE){

  #grab metadata and data
  meta <- x[,1:ecs]
  x <- x[,(ecs+1):ncol(x)]


  ##############################
  #sanity checks...
  if((mc != FALSE & !length(counts) > 1) | (mc == FALSE & length(counts) > 1)){
    stop("Counts and mc must either both be defined or neither must be.")
  }

  if(mc != FALSE & length(counts) > 1){
    if(!is.numeric(mc) | !is.numeric(counts)){
      stop("Counts and mc must both be numeric vectors.\n")
    }
    if(length(mc) > 1){
      stop("mc must a single numeric value equal to the desired cuttoff number of sequenced loci.\n")
    }
    if(mc > (ncol(x) - ecs)/2){
      stop("mc must be less than the total number of sequenced loci.\n")
    }
    if(length(counts) != nrow(x)){
      stop("counts must be equal in length to the number of samples in x.")
    }
  }

  if(is.character(do.plot)){
    if(length(do.plot) > 2){
      stop("Only two plotting variables supported. For more, set do.plot to FALSE and plot manually.\n")
    }
  }
  else if (do.plot != FALSE){
    stop("do.plot must be either FALSE or between 1 and 2 variables to plot by.")
  }
  else if (!all(do.plot) %in% colnames(meta)){
    stop("Plotting variables specified in do.plot must match column names present in the metadata of x.\n")
  }

  ############################################
  #filter if requested
  if(is.numeric(mc) & is.numeric(counts)){
    keeps <- which(counts >= mc)
    x <- x[keeps,]
    meta <- meta[keeps,]
  }

  #check for any duplicates, which need to be removed!
  if(c.dup){
    cat("Checking for duplicates...\n")
    dups <- which(duplicated(x) | duplicated(x, fromLast=TRUE))
    if(length(dups) > 0){
      cat("Duplicates detected, indices:", dups, "\nRemoving all of these!\n")
      x <- x[-dups,]
      meta <- meta[-dups,]
    }
  }

  cat("Preparing pca...\n")
  pca_r <- prcomp(as.matrix(x))
  pca <- as.data.frame(pca_r$x) #grab the PCA vectors.
  pca <- cbind(meta, pca)  #add metadata that is present in the input.


  #return just the pca data if the plot isn't requested, which is useful for complex plots.
  if(!is.character(do.plot)){
    return(list(raw = pca_r, pca = pca))
  }

  ################################################
  #construct plot.
  cat("Preparing plot...\n")

  #Categories (pops, fathers, mothers, ect.) are given in do.plot argument. Supports up to two!
  #make the base plot, then add categories as color and fill.
  out <- ggplot2::ggplot(pca, ggplot2::aes(PC1, PC2)) + ggplot2::theme_bw() #initialize plot

  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


  #add variables.
  if(do.plot[1] != FALSE){

    v1 <- pca[,which(colnames(pca) == do.plot[1])] #get the factors
    v1u <- length(unique(v1)) #number of categories

    # add geoms to plot
    if(length(do.plot) == 1){
      out <- out + ggplot2::geom_point(ggplot2::aes(color = v1))#add the factor
    }
    else{
      # if two plotting variables, prepare the second and add it as well.
      v2 <- pca[,which(colnames(pca) == do.plot[2])]
      v2u <- length(unique(v2))
      out <- out + ggplot2::geom_point(ggplot2::aes(color = v1, fill = v2), pch = 21, size = 2.5, stroke = 1.25)
    }


    # change the color scales for the variables
    ## for the first variable
    if(v1u <= 8){#are there too many categories to color with the cbb palette?
      out <- out + ggplot2::scale_color_manual(values = cbbPalette, name = do.plot[1])
    }
    else{
      # if the data is likely continuous:
      if(is.numeric(v1)){
        warning("Plotting ", do.plot[1], " as a continuous variable.\n")
        out <- out + ggplot2::scale_color_viridis_c(name = do.plot[1])
      }
      out <- out + ggplot2::scale_color_viridis_d(name = do.plot[1])
    }

    ## for the second variable if defined
    if(length(do.plot) == 2){
      if(v2u <= 8){#are there too many categories to color with the cbb palette?
        out <- out + ggplot2::scale_fill_manual(values = cbbPalette, name = do.plot[2])
      }
      else{
        if(is.numeric(v2)){
          warning("Plotting ", do.plot[2], " as a continuous variable.\n")
          out <- out + ggplot2::scale_fill_viridis_c(name = do.plot[2])
        }
        else{
          out <- out + ggplot2::scale_fill_viridis_d(name = do.plot[2])
        }
      }
    }
  }

  loadings <- (pca_r$sdev^2)/sum(pca_r$sdev^2)
  loadings <- round(loadings * 100, 2)

  out <- out + ggplot2::xlab(paste0("PC1 (", loadings[1], "%)")) + ggplot2::ylab(paste0("PC2 (", loadings[2], "%)"))

  return(list(raw = pca_r, plot = out))
}


# tsne fxn using Rtsne, which uses Barnes-Hut simulation at theta>0 and returns more data than tsne()

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
#' @param do.plot FALSE or character vector, default "pop". If FALSE, no plot is produced. Up to two variables to be plotted with names corresponding to column names in x be given as a character vector.
#' @param dims Integer, output dimensionality, default 2.
#' @param initial_dims Integer, default 50. The number of dimensions retained in the initial PCA step.
#' @param perplex Perplexity parameter, by default found by \code{\link[mmtsne]{hbeta}}, with beta = 1.
#' @param gravity Theta parameter from \code{\link[Rtsne]{Rtsne}}. Default 0, an exhaustive search.
#' @param iter Integer, default 5000. Number of tSNE iterations to perform.
#' @param c.dup boolean, default FALSE. Should duplicate individuals be searched for and removed? This is very slow if the data set is large!
#' @param mc Numeric or FALSE, default FALSE. Should poorly sequenced individuals be removed? If so, what is the minimum acceptable count of sequenced loci (as specificed in the vector provided to counts)?
#' @param counts Numeric vector, default FALSE, containing the number of loci sequenced per individual.
#' @param ... Other arguments, passed to \code{\link[Rtsne]{Rtsne}}.
#'
#' @return A list containing the raw tSNE output and a tSNE plot in the form of a ggplot graphical object. The plot can be changed as usual with ggplot objects.
#'
#' @examples
#' tSNEfromPA(stickPA, 2, c.dup = TRUE)
tSNEfromPA <- function(x, ecs, do.plot = "pop", dims = 2, initial_dims = 50,
                       perplex = FALSE, gravity = 0, iter = 5000,
                       c.dup = FALSE, mc = FALSE, counts = FALSE, ...){
  library(ggplot2)
  #grab metadata and data
  meta <- x[,1:ecs]
  x <- x[,(ecs+1):ncol(x)]
  x <- as.matrix(x)

  ##############################
  #sanity checks...
  if((mc != FALSE & !length(counts) > 1) | (mc == FALSE & length(counts) > 1)){
    stop("Counts and mc must either both be defined or neither must be.")
  }

  if(mc != FALSE & length(counts) > 1){
    if(!is.numeric(mc) | !is.numeric(counts)){
      stop("Counts and mc must both be numeric vectors.\n")
    }
    if(length(mc) > 1){
      stop("mc must a single numeric value equal to the desired cuttoff number of sequenced loci.\n")
    }
    if(mc > (ncol(x) - ecs)/2){
      stop("mc must be less than the total number of sequenced loci.\n")
    }
    if(length(counts) != nrow(x)){
      stop("counts must be equal in length to the number of samples in x.")
    }
  }

  if(c.dup == FALSE){
    warning("If there are duplicates in x, expect wierd results! Set c.dup to TRUE to check.\n")
  }

  if(is.character(do.plot)){
    if(length(do.plot) > 2){
      stop("Only two plotting variables supported. For more, set do.plot to FALSE and plot manually.\n")
    }
  }
  else if (do.plot != FALSE){
    stop("do.plot must be either FALSE or between 1 and 2 variables to plot by.")
  }
  else if (!all(do.plot) %in% colnames(meta)){
    stop("Plotting variables specified in do.plot must match column names present in the metadata of x.\n")
  }

  ##############################
  #filter if requested
  if(is.numeric(mc) & is.numeric(counts)){
    keeps <- which(counts >= mc)
    x <- x[keeps,]
    meta <- meta[keeps,]
  }


  if(c.dup){
    #check for any duplicates, which need to be removed!
    cat("Checking for duplicates...\n")
    dups <- which(duplicated(x) | duplicated(x, fromLast=TRUE))
    if(length(dups) > 0){
      cat("Duplicates detected, indices:", dups, "\nRemoving all of these!\n")
      x <- x[-dups,]
      meta <- meta[-dups,]
    }
  }

  #get perplexity if not provided
  if(!perplex){
    cat("Estimating perplexity...\n")
    perp <- mmtsne::hbeta(x, beta = 1)
    perplex <- perp$H
  }

  #run the tSNE
  cat("Running tSNE...\n")
  tsne.out <- Rtsne::Rtsne(x, dims, initial_dims, perplex,
                           gravity, iter, check_duplicates = FALSE,
                           verbose=TRUE, ...)
  tsne_plot <- cbind(meta, as.data.frame(tsne.out$Y))

  if(!is.character(do.plot)){
    return(list(raw = tsne.out, pca = tsne_plot))
  }

  #############################
  #plot the result
  cat("Preparing plot...\n")

  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  #Categories (pops, fathers, mothers, ect.) are given in do.plot argument. Supports up to two!
  #make the base plot, then add categories as color and fill.
  out <- ggplot2::ggplot(tsne_plot, ggplot2::aes(V1, V2)) + ggplot2::theme_bw() +
    ggplot2::theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          panel.grid = element_blank()) #initialize plot

  #add variables.
  if(do.plot[1] != FALSE){

    v1 <- tsne_plot[,which(colnames(tsne_plot) == do.plot[1])] #get the factors
    v1u <- length(unique(v1)) #number of categories

    # add geoms to plot
    if(length(do.plot) == 1){
      out <- out + ggplot2::geom_point(ggplot2::aes(color = v1))#add the factor
    }
    else{
      # if two plotting variables, prepare the second and add it as well.
      v2 <- tsne_plot[,which(colnames(tsne_plot) == do.plot[2])]
      v2u <- length(unique(v2))
      out <- out + ggplot2::geom_point(ggplot2::aes(color = v1, fill = v2), pch = 21, size = 2.5, stroke = 1.25)
    }


    # change the color scales for the variables
    ## for the first variable
    if(v1u <= 8){#are there too many categories to color with the cbb palette?
      out <- out + ggplot2::scale_color_manual(values = cbbPalette, name = do.plot[1])
    }
    else{
      # if the data is likely continuous:
      if(is.numeric(v1)){
        warning("Plotting ", do.plot[1], " as a continuous variable.\n")
        out <- out + ggplot2::scale_color_viridis_c(name = do.plot[1])
      }
      out <- out + ggplot2::scale_color_viridis_d(name = do.plot[1])
    }

    ## for the second variable if defined
    if(length(do.plot) == 2){
      if(v2u <= 8){#are there too many categories to color with the cbb palette?
        out <- out + ggplot2::scale_fill_manual(values = cbbPalette, name = do.plot[2])
      }
      else{
        if(is.numeric(v2)){
          warning("Plotting ", do.plot[2], " as a continuous variable.\n")
          out <- out + ggplot2::scale_fill_viridis_c(name = do.plot[2])
        }
        else{
          out <- out + ggplot2::scale_fill_viridis_d(name = do.plot[2])
        }
      }
    }
  }

  return(list(tSNE = tsne.out, plot = out))
}
