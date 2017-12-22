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
LD_pairwise_heatmap <- function(x, r = NULL,
                                l.text = "rsq", colors = c("white", "black"),
                                title = "Pairwise LD", t.sizes = c(16, 13, 10, 12, 10),
                                background = "white"){

  #Created in part by Nick Sard

  #remove columns and rows with no data
  x <- x[!apply(x, 1, function(y)all(is.na(y))), !apply(x, 2, function(y)all(is.na(y)))]

  #if first column doesn't contain positions, take the row names and set as first column.
  if(any(is.na(x[,1]))){
    x <- cbind(V1 = row.names(x), x)
  }


  #melting the df so its computer reable and fixing the names of the columns
  heatmap_x <- melt(x, id.vars = "V1")
  names(heatmap_x) <- c("SNPa", "SNPb", "value")

  #getting rid of all the zeros from snps being compared to themselves
  heatmap_x$SNPa <- as.numeric(as.character(heatmap_x$SNPa))
  heatmap_x$SNPb <- as.numeric(as.character(heatmap_x$SNPb))
  heatmap_x <- heatmap_x[!(heatmap_x$SNPa == heatmap_x$SNPb),]

  #removing NA values, since no Dups.  Note, also grabbing only section of interest
  heatmap_x <- heatmap_x[!is.na(heatmap_x$value),]

  #remove any other NAs
  heatmap_x <- heatmap_x[!is.na(heatmap_x$SNPa) & !is.na(heatmap_x$SNPb),]

  #get site names
  ms <- unique(c(colnames(x), x[,1]))
  ms <- c(x[1,1], colnames(x)[-1])
  ms <- as.numeric(as.character(ms))

  #subset down to the desired r if requested
  if(!is.null(r)){
    r <- r*1000000
    heatmap_x <- heatmap_x[heatmap_x$SNPa >= r[1] & heatmap_x$SNPa <= r[2] &
                             heatmap_x$SNPb >= r[1] & heatmap_x$SNPb <= r[2],]
    ms <- ms[ms >= r[1] & ms <= r[2]]
  }

  #finish messing with site names
  ms <- ms/1000000
  ms <- as.factor(ms)
  ms <- sort(ms)

  #set site positions to mb, convert to factor
  heatmap_x$SNPa <- heatmap_x$SNPa/1000000
  heatmap_x$SNPb <- heatmap_x$SNPb/1000000
  heatmap_x$SNPa <- as.factor(heatmap_x$SNPa)
  heatmap_x$SNPb <- as.factor(heatmap_x$SNPb)

  #reordering based on factors
  heatmap_x[["SNPa"]]<-factor(heatmap_x[["SNPa"]],levels= ms,ordered=T)
  heatmap_x[["SNPb"]]<-factor(heatmap_x[["SNPb"]],levels=rev(ms),ordered=T)

  #the plot
  out <- ggplot(heatmap_x, aes(x = SNPa, y=SNPb, fill=value))+
    geom_tile(color = "white")+
    scale_fill_gradient(low = colors[1], high = colors[2]) +
    theme_bw()+
    labs(x = "",y="", fill=l.text)+
    theme(legend.title= element_text(size = t.sizes[2]),
          axis.text = element_text(size = t.sizes[5]),
          panel.grid.major = element_line(color = background),
          plot.title = element_text(size = t.sizes[1], hjust = 0.5),
          axis.title = element_text(size = t.sizes[4]),
          legend.text = element_text(size = t.sizes[3]),
          panel.background = element_rect(fill = background, colour = background)) +
    scale_x_discrete(breaks = levels(heatmap_x$SNPa)[c(T, rep(F, 20))], label = abbreviate) +
    scale_y_discrete(breaks = levels(heatmap_x$SNPb)[c(T, rep(F, 20))], label = abbreviate) +
    ggtitle(title) + ylab("Position (Mb)") + xlab("Position (Mb)")
  return(out)
}
