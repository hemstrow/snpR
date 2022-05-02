import_wrapper <- function(geno, snp, samp, mDat){
  if(!isFALSE(snp) & !isFALSE(samp)){
    x <- import.snpR.data(geno$datapath, snp$datapath, samp$datapath, mDat)
  }
  else if(isFALSE(snp)){
    x <- import.snpR.data(geno$datapath, FALSE, samp$datapath, mDat)
  }
  else if(isFALSE(samp)){
    x <- import.snpR.data(geno$datapath, snp$datapath, FALSE, mDat)
  }
  
  x <- calc_fis(x)
  return(x)
}

validation <- function(x){
  maf_spectra <- ggplot2::ggplot(get.snpR.stats(x, ".base", stats = "maf")$single,
                                 ggplot2::aes(x = maf)) +
    ggplot2::geom_density() + 
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean(maf))) +
    ggplot2::theme_bw() + ggplot2::ggtitle("Maf Density")
  
  fis_spectra <- ggplot2::ggplot(get.snpR.stats(x, ".base", stats = "fis")$single,
                                 ggplot2::aes(x = fis)) +
    ggplot2::geom_density() + 
    ggplot2::geom_vline(ggplot2::aes(xintercept = mean(fis))) +
    ggplot2::theme_bw() + ggplot2::ggtitle("Fis Density")
  
  spectra <- gridExtra::arrangeGrob(maf_spectra, fis_spectra, ncol = 1)
  
  facet_options <- list(snp = colnames(x@sample.meta),
                        sample = colnames(x@snp.meta))
  
  return(list(spec = spectra, 
              nsnps = nsnps(x), 
              nsamps = nsamps(x),
              facet_options = facet_options))
}

renderValidation <- function(valid){
  vp <- renderPlot({
    gridExtra::grid.arrange(valid$spec)
  })
  vn <- renderText({
    paste0("Number of SNPs: ", valid$nsnps, "\nNumber of Samples: ", valid$nsamps)
  })
  
  return(list(vp = vp, vn = vn))
}

filter_wrapper <- function(x, maf, maf_facet, hwe, hwe_facet, min_ind, min_loci){
  
  maf_facet <- paste0(maf_facet, collapse = ".")
  hwe_facet <- paste0(hwe_facet, collapse = ".")
  
  x <- filter_snps(x,
                   maf = maf, 
                   maf_facets = maf_facet,
                   hwe = hwe,
                   hwe_facets = hwe_facet,
                   min_ind = min_ind,
                   min_loci = min_loci)
  
  x <- calc_fis(x)
  return(x)
}

facet_check <- function(x){
  return(list(snp = colnames(snp.meta(x)),
              sample = colnames(sample.meta(x))))
}