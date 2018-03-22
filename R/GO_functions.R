#Takes an input file of genes with start/end position noted and maps back the average
#statitstic of choice for that gene from another provided data set with positions. Requires the "zoo"
#package

#' Interpolates stats at genes.
#'
#' \code{gene_ave_stat} takes an input file of genes with start/end position noted and interpolates back the average statitstic of choice for that gene. Requires the "zoo" package.
#'
#'Description of gene_data:
#'    Requires columns titled "group", "start", "end", and "probeID", containing gene linkage group/chromosome, start and end positions, and name, respectively.
#'
#'Description of stat_data:
#'    Requires columns titled "group", "position", and one matching the stat argument, containing the linkage group/chromosome, position, and the value of the observed statistic at that position. Typically produced via windowed gaussian smoothing (from smoothed_ave).
#'
#'Uses spline interpolation from the "zoo" package.
#'
#' @param x Input gene data, with columns named "start", "end" and "probeID".
#' @param y Input stat data, typically from smoothed windows or (not recommended) raw statistics from SNPs.
#' @param stat Name of the statistic to interpolate.
#' @examples
#' gene_ave_stat(stickleGO, randSMOOTHed[randSMOOTHed$pop == "A",], "smoothed_pi")
#'
gene_ave_stat <- function(x, y, stat){

  gdat <- x
  sdat <- y

  if(nrow(y) == 0){
    warning("No stat data provided.\n")
    return()
  }

  #prepare data
  gdat$mid <- rowMeans(gdat[,c(2:3)])
  gdat$probeID <- as.character(gdat$probeID)
  cdf <- data.frame(mid = c(gdat$mid, sdat$position),
                    stat = c(rep(NA, nrow(gdat)), sdat[,stat]),
                    probeID = c(gdat$probeID, rep("snp", nrow(sdat))))

  #use zoo to interplate
  cdf <- dplyr::arrange(cdf, mid)
  cdf <- zoo::zoo(cdf)
  zoo::index(cdf) <- cdf$mid
  cdf$stat <- zoo::na.spline(cdf$stat)

  #clean up
  cdf <- as.data.frame(cdf, stringsAsFactors = F)
  cdf <- cdf[cdf$probeID != "snp",]
  cdf$stat <- as.numeric(cdf$stat)
  cdf$probeID <- as.character(cdf$probeID)
  cdf$mid <- as.numeric(as.character(cdf$mid))

  #remerge
  gdat <- merge(gdat, cdf, by = c("probeID", "mid"))
  gdat <- gdat[,colnames(gdat) != "mid"]
  colnames(gdat)[which(colnames(gdat) == "stat")] <- stat

  return(gdat)
}


