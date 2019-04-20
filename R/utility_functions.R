#' Import genotype and metadata into a snpRdata object.
#'
#' \code{import.snpR.data} converts genotype and meta data to the snpRdata class, which stores raw genotype data, sample and locus specific metadata, useful data summaries, repeatedly internally used tables, calculated summary statistics, and smoothed statistic data.
#'
#' The snpRdata class is built to contain SNP genotype data for use by functions in the snpR package. It inherits from the S3 class data.frame, in which the genotypes are stored, and can be manipulated identically. It also stores sample and locus specific metadata, genomic summary information, and any results from most snpR functions. The raw data for each of these latter objects is accessable via the at operator.
#' Genotypes are stored in the "character" format, as output by format_snps(). Missing data is noted with "NN".
#'
import.snpR.data <- function(genotypes, snp.meta, sample.meta, mDat){
  # prepare things for addition to data
  if(any(colnames(snp.meta) == "position")){
    snp.meta$position <- as.numeric(as.character(snp.meta$position))
    genotypes <- genotypes[order(snp.meta$position),]
    snp.meta <- dplyr::arrange(snp.meta, position)
  }

  if(any(colnames(snp.meta) == ".snp.id")){
    if(any(duplicated(snp.meta$.snp.id))){stop("Duplicated .snp.id entries found in snp.meta.\n")}
  }
  else{
    snp.meta <- cbind(snp.meta, .snp.id = 1:nrow(snp.meta))
  }
  if(any(colnames(sample.meta) == ".sample.id")){
    if(any(duplicated(sample.meta$.sample.id))){stop("Duplicated .sample.id entries found in sample.meta.\n")}
  }
  else{
    sample.meta <- cbind(sample.meta, .sample.id = 1:nrow(sample.meta))
  }

  rownames(genotypes) <- 1:nrow(genotypes)
  rownames(snp.meta) <- 1:nrow(snp.meta)

  gs <- tabulate_genotypes(genotypes, mDat = mDat, verbose = T)

  fm <- data.frame(facet = rep(".base", nrow(gs$gs)),
                   subfacet = rep(".base", nrow(gs$gs)),
                   facet.type = rep(".base", nrow(gs$gs)))

  fm <- cbind(fm, snp.meta)

  x <- new("snpRdata", .Data = genotypes, sample.meta = sample.meta, snp.meta = snp.meta,
           facet.meta = fm,
           geno.tables = gs,
           mDat = mDat,
           stats = fm,
           snp.form = nchar(genotypes[1,1]), row.names = rownames(genotypes),
           facets = ".base",
           facet.type = ".base")

  # run essential filters (np, bi-al), since otherwise many of the downstream applications, including ac formatting, will be screwy.
  cat("Imput data will be filtered to remove non bi-allelic data.\n")
  invisible(capture.output(x <- filter_snps(x, non_poly = F)))

  # add ac
  invisible(capture.output(x@ac <- format_snps(x, "ac")[,c("n_total", "n_alleles", "ni1", "ni2")]))

  return(x)
}


# function to add a new facet to snpRdata, generating gs, as, and wmat tables, and ac formatted data.
# need to add the ac part.
add.facets.snpR.data <- function(x, facets = NULL){
  if(is.null(facets[1])){return(x)}
  facets <- check.snpR.facet.request(x, facets)
  #===========================turn into list========
  # need to fix any multivariate facets (those with a .)
  comp.facets <- grep("(?<!^)\\.", facets, perl = T)
  if(length(comp.facets) != 0){
    run.facets <- as.list(facets[-c(comp.facets)])
    facet.list <- c(run.facets, strsplit(facets[comp.facets], split = "(?<!^)\\.", perl = T))
  }
  else{
    facet.list <- list(facets)
  }

  #===========================sanity checks=========
  # check that we haven't done these facets before, and remove any that we have.
  all.facets <- character(length = length(facet.list))
  for(i in 1:length(facet.list)){
    all.facets[i] <- paste0(facet.list[[i]], collapse = ".")
  }
  if(all(all.facets %in% x@facets)){return(x)}
  else if (any(all.facets %in% x@facets)){
    already.added <- all.facets[which(all.facets %in% x@facets)]
    warning(paste0("Some facets already present in data:\n\t", paste0(already.added, collapse = "\n\t"), "\n"))
    facet.list <- facet.list[-which(all.facets %in% x@facets)]
  }

  # check that all of the facet columns actually exist
  all.facets <- unlist(unlist(facet.list))
  opts <- c(colnames(x@sample.meta), colnames(x@snp.meta))
  used <- opts[which(opts %in% all.facets)]
  if(!all(all.facets %in% opts)){
    bad.facets <- which(!(all.facets %in% c(colnames(x@sample.meta), colnames(x@snp.meta))))
    stop(paste0("Facets ", paste0(all.facets[bad.facets], collapse = " "), " not found in sample or snp metadata.\n"))
  }
  if(any(duplicated(used))){
    stop(paste0("Facets ", paste0(used[which(duplicated(used))], collapse = " "), " are duplicated in the sample and or snp metadata.\n"))
  }

  #===========================process each facet.===================
  # for sample acets, prep summary data and geno tables.
  # For snp facets, nothing to do here besides add them to the facet list and facet type list.
  # For complex facets, prep summary tables for the sample facet portion if they don't exist and add to facet list.
  added.facets <- character(0)
  oac <- cbind(data.table::as.data.table(x@facet.meta[,c("facet", "subfacet", ".snp.id")]),
               data.table::as.data.table(x@ac)) # grab original ac with meta for later
  for(k in 1:length(facet.list)){
    facets <- facet.list[[k]] # column levels for this facet.
    #=========================figure out unique levels for the facet==========
    # figure out what kind of facets we are working with.

    # save info and get the unique options for each facet
    # save info
    if(length(facets) > 1){
      x@facets <- c(x@facets, paste0(facets, collapse = "."))
      added.facets <- c(added.facets, paste0(facets, collapse = "."))
    }
    else{
      x@facets <- c(x@facets, facets)
      added.facets <- c(added.facets, facets)
    }
    x@facet.type <- c(x@facet.type, "sample")


    # get unique options for this facet
    sample.meta <- x@sample.meta[colnames(x@sample.meta) %in% facets]
    sample.meta <- sample.meta[,sort(colnames(sample.meta))]
    if(!is.data.frame(sample.meta)){
      sample.meta <- as.data.frame(sample.meta)
      colnames(sample.meta) <- colnames(x@sample.meta)[colnames(x@sample.meta) %in% facets]
    }
    sample.opts <- unique(sample.meta)
    if(!is.data.frame(sample.opts)){
      sample.opts <- as.data.frame(sample.opts, stringsAsFactors = F)
      colnames(sample.opts) <- facets[which(facets %in% colnames(x@sample.meta))]
    }
    sample.opts <- dplyr::arrange_all(sample.opts)


    gs <- x@geno.tables
    #=========================get gs matrices==========
    for(i in 1:nrow(sample.opts)){
      matches <- which(apply(sample.meta, 1, function(x) identical(as.character(x), as.character(sample.opts[i,]))))
      t.x <- x[,matches]
      tgs <- tabulate_genotypes(t.x, x@mDat)
      gs$gs <- plyr::rbind.fill.matrix(gs$gs, tgs$gs)
      gs$as <- plyr::rbind.fill.matrix(gs$as, tgs$as)
      gs$wm <- plyr::rbind.fill.matrix(gs$wm, tgs$wm)
      x@facet.meta <- rbind(x@facet.meta,
                            cbind(data.frame(facet = rep(paste0(facets, collapse = "."), nrow(tgs$gs)),
                                             subfacet = rep(paste0(sample.opts[i,], collapse = "."), nrow(tgs$gs)),
                                             facet.type = rep("sample", nrow(tgs$gs)), stringsAsFactors = F),
                                  x@snp.meta))
    }


    #=========================sort, pack, and return==========
    # sort
    x@facet.meta$.reorder <- 1:nrow(x@facet.meta)
    x@facet.meta <- dplyr::arrange(x@facet.meta, .snp.id, facet, subfacet)
    gs$gs <- gs$gs[x@facet.meta$.reorder,]
    gs$as <- gs$as[x@facet.meta$.reorder,]
    gs$wm <- gs$wm[x@facet.meta$.reorder,]
    x@facet.meta <- x@facet.meta[,-ncol(x@facet.meta)]

    # output
    x@geno.tables <- gs
  }
  # add and sort ac formated data.
  invisible(capture.output(nac <- format_snps(x, output = "ac", facets = added.facets)))
  nac <- data.table::as.data.table(nac)
  nac <- rbind(oac, nac[,c("facet", "subfacet", ".snp.id", "n_total","n_alleles", "ni1", "ni2")])
  nac <- dplyr::arrange(nac, .snp.id, facet, subfacet)
  x@ac <- nac[,-c(1:3)]

  # add in dummy rows to stats
  sm <- data.table::as.data.table(x@facet.meta[x@facet.meta$facet %in% added.facets, c("facet", "subfacet", "facet.type", colnames(x@snp.meta))])
  if(ncol(x@stats) - ncol(sm) > 0){
    sm <- cbind(sm, matrix(NA, nrow(sm), ncol(x@stats) - ncol(sm)))
  }
  colnames(sm) <- colnames(x@stats)
  os <- data.table::as.data.table(x@stats)
  if(ncol(os) - ncol(sm) > 0){
    os <- rbind(os, cbind(sm, matrix(NA, nrow(sm), ncol(os) - ncol(sm))))
  }
  else{
    os <- rbind(os, sm)
  }
  x@stats <- dplyr::arrange(os, .snp.id, facet, subfacet)

  return(x)
}

# function to list existing facets.
find.snpR.facets <- function(x){
  facets <- vector("list", length(x@facets))
  names(facets) <- x@facets
  for(i in 1:length(facets)){
    facets[[i]] <- unique(x@facet.meta[x@facet.meta$facet == names(facets)[i],]$subfacet)
  }
  return(facets)
}

# function to pull stats for a given facet
get.snpR.stats <- function(x, facets = NULL){
  if(!is.null(facets)){
    if(facets[1] == "all"){
      facets <- x@facets
    }
  }
  else {
    facets <- ".base"
  }

  for(i in 1:length(facets)){
    if(any(facets %in% colnames(x@snp.meta))){
      if(any(facets %in% colnames(x@sample.meta))){
        samp.facets <- which(facets %in% colnames(x@sample.meta))
        facets <- facets[samp.facets]
      }
      else{
        next() # snp level facet, skip.
      }
    }
  }

  return(x@stats[which(x@stats$facet %in% facets), -which(colnames(x@stats) %in% c(".snp.id", "facet.type"))])
}

# function to apply a function across selected facets
# req: which part of the snpR.data object is required and should be pulled out?
#        gs: genotype tables
#        ac: ac formatted data
#        meta.gs: facet, .snp.id metadata cbound genotype tables.
#        ac.stats: ac data cbound to stats
# facets: if NULL, run without facets.
# cases: ps: per snp.
#        ps.pf: per snp, but split per facet (such as for private alleles, comparisons only exist within a facet!)
# fun: which function should be applied?
apply.snpR.facets <- function(x, facets = NULL, req, fun, case = "ps", ...){
  if(!is.null(facets)){
    if(facets[1] == "all"){
      facets <- x@facets
    }
  }
  else {
    facets <- ".base"
  }


  if(case == "ps"){
    facets <- check.snpR.facet.request(x, facets)
    if(req == "gs"){
      # subset
      gs <- plyr::llply(x@geno.tables, function(y) subset(y, x@facet.meta$facet %in% facets))

      # run the function indicated
      out <- data.table::as.data.table(fun(gs, ...))

      # bind metadata
      out <- cbind(data.table::as.data.table(x@facet.meta[x@facet.meta$facet %in% facets,]), out)
      #out <- cbind(x@facet.meta[x@facet.meta$facet %in% facets,], out)

      # return
      return(out)
    }
    else if(req == "meta.gs"){
      # gs with metadata on facets attached.
      # subset
      gs <- plyr::llply(x@geno.tables, function(y) subset(y, x@facet.meta$facet %in% facets))
      gs <- plyr::llply(gs, function(y) cbind(x@facet.meta[x@facet.meta$facet %in% facets, c("facet", "subfacet", ".snp.id")], y))

      # run the function indicated
      out <- data.table::as.data.table(fun(gs, ...))

      # bind metadata
      out <- cbind(data.table::as.data.table(x@facet.meta[x@facet.meta$facet %in% facets,]), out)

      # return
      return(out)
    }
    else if(req == "ac"){ # really simple, just here for consistancy across functions. you just... run it on the thing.
      # subset
      ac <- x@ac[x@facet.meta$facet %in% facets,]

      # run the function indicated
      out <- data.table::as.data.table(fun(ac, ...))

      # bind metadata
      out <- cbind(data.table::as.data.table(x@facet.meta[x@facet.meta$facet %in% facets,]), out)

      # return
      return(out)
    }
  }
  else if(case == "ps.pf"){
    out <- data.frame() # initialize
    if(req == "meta.gs"){
      # gs with metadata on facets attached.
      # loop through each facet
      for(i in 1:length(facets)){
        # subset
        gs <- plyr::llply(x@geno.tables, function(y) subset(y, x@facet.meta$facet == facets[i]))
        gs <- plyr::llply(gs, function(y) cbind(x@facet.meta[x@facet.meta$facet == facets[i], c("facet", "subfacet", ".snp.id")], y))

        # run the function indicated
        this.out <- data.table::as.data.table(fun(gs, ...))

        # bind metadata
        this.out <- data.table::as.data.table(cbind(x@facet.meta[x@facet.meta$facet == facets[i],], this.out))

        # bind to full data
        out <- rbind(out, this.out)
      }

      # return
      return(out)
    }
  }
  else if(case == "facet.pairwise"){
    if(req == "snpRdata"){
      # initialize
      out1 <- data.frame()
      out2 <- data.frame()


      # in this case, the whole snpRdata object is requested, but should be passed several times with the correct facet noted.
      for(i in 1:length(facets)){
        out <- fun(x, facets[i], ...)
        # if this is a genepop fst output, we should expect a list of size two. Need to return the correct parts!
        if(length(out) == 2){
          suppressWarnings(out1 <- rbind(out1, cbind(data.table::as.data.table(x@snp.meta), facet = facets[i], data.table::as.data.table(out[[1]]))))
          out2 <- rbind(out2, data.frame(comparison = names(out[[2]]), overall_fst = out[[2]]))
        }
      }

      # return
      if(nrow(out2) > 0){
        return(list(out1, out2))
      }
      else{
        return(out1)
      }
    }


    else if (req == "ac.stats"){
      out <- data.frame()
      # need to loop through each facet, since they are all internally compared!
      for(i in 1:length(facets)){
        # subset
        ac <- x@ac[x@facet.meta$facet == facets[i],]
        stats <- x@stats[x@facet.meta$facet == facets[i],]

        # run the function indicated
        this.out <- data.table::as.data.table(fun(cbind(stats, ac), ...))

        # bind metadata and combine with full output
        suppressWarnings(out <- rbind(out, cbind(data.table::as.data.table(x@snp.meta), facet = facets[i], this.out)))

      }

      # return
      return(out)

    }
  }
}

# merge newly calculated stats with a snpRdata object.
merge.snpR.stats <- function(x, stats, type = "stats"){
  if(type == "stats"){
    o.s <- data.table::as.data.table(x@stats)
    n.s<- data.table::as.data.table(stats)
    # if new stats, easy to merge
    if(!all(colnames(n.s) %in% colnames(o.s))){
      o.s <- merge(o.s, n.s,
                       by = colnames(o.s)[which(colnames(o.s) %in% colnames(n.s))], # note here that we are only merging by columns that already exist in both datasets
                       all = T, sort = F)
      o.s <- data.table::setorder(o.s, .snp.id, facet, subfacet)
    }

    # otherwise need to overwrite
    else{
      add.rows <- which(o.s$facet %in% n.s$facet)
      y <- o.s[add.rows,]

      # sort identically
      #y <- dplyr::arrange(y, .snp.id, facet, subfacet)
      n.s <- data.table::setorder(n.s, .snp.id, facet, subfacet)
      y[,which(colnames(y) %in% colnames(n.s))] <- n.s

      #replace and add
      o.s[add.rows,] <- y # add back
    }

    # sort and return
    x@stats <- o.s
  }
  else if(type == "pairwise"){
    o.s <- data.table::as.data.table(x@pairwise.stats)
    n.s <- data.table::as.data.table(stats)
    # since pairwise.stats data.frames aren't automatically initizalized, need to add new rows if they don't already exist.
    missing.rows <- which(!(paste0(n.s$facet, "__", n.s$comparison) %in% paste0(o.s$facet, "__", o.s$comparison)))
    if(length(missing.rows) > 0){
      new.rows <- n.s[,1:which(colnames(n.s) == "comparison")]
      if(ncol(o.s) > ncol(new.rows)){
        new.rows <- data.table::as.data.table(cbind(new.rows, matrix(NA, nrow = nrow(new.rows), ncol = (ncol(o.s) - ncol(new.rows)))))
        colnames(new.rows) <- colnames(o.s)
      }
      o.s <- rbind(o.s, new.rows)
      o.s <- data.table::setorder(o.s, .snp.id, facet, comparison)
    }

    # if new stats, easy to merge
    if(!all(colnames(n.s) %in% colnames(o.s))){
      o.s <- merge(o.s, n.s,
                       by = colnames(o.s)[which(colnames(o.s) %in% colnames(n.s))], # note here that we are only merging by columns that already exist in both datasets
                       all = T, sort = F)
      o.s <- data.table::setorder(o.s, .snp.id, facet, comparison)
    }

    # otherwise need to overwrite
    else{
      add.rows <- which(o.s$facet %in% n.s$facet)
      y <- o.s[add.rows,]

      # sort identically
      # y <- data.table::setorder(y, .snp.id, facet, comparison) # should already be sorted, and this can be slow.
      n.s <- data.table::setorder(n.s, .snp.id, facet, comparison)
      y[,which(colnames(y) %in% colnames(n.s))] <- n.s

      #replace and add
      o.s <- data.table::set(o.s, i = add.rows, j = 1:ncol(o.s), value = y)
    }

    # sort and return -- should be able to skip this.
    # x@pairwise.stats <- data.table::setorder(o.s, .snp.id, facet, comparison)
    x@pairwise.stats <- o.s
  }
  else if(type == "LD"){

    if(length(x@pairwise.LD) == 0){
      x@pairwise.LD <- stats
      return(x)
    }
    else{
      # deal with prox
      prox.check_x <- do.call(paste, as.data.frame(x@pairwise.LD$prox[,which(!colnames(x@pairwise.LD$prox) %in% c("rsq", "proximity", "Dprime", "pval"))]))
      prox.check_stats <- do.call(paste, as.data.frame(stats$prox[,which(!colnames(stats$prox) %in% c("rsq", "proximity", "Dprime", "pval"))]))
      prox.check <- which(!(prox.check_stats %in% prox.check_x))
      if(length(prox.check) != 0){
        x@pairwise.LD$prox <- rbind(x@pairwise.LD$prox, stats$prox[prox.check,])
      }

      # deal with matrices. Overwrite any matrices that already exist and add any new ones.
      facets <- names(x@pairwise.LD$LD_matrices)
      overwrite.facets <- which(facets %in% names(stats$LD_matrices)) # which existing facets need to be overwritten?
      add.facets <- which(!(names(stats$LD_matrices) %in% facets)) # which LD facets need to be added?

      ## overwrite
      if(length(overwrite.facets) != 0){
        x@pairwise.LD$LD_matrices[overwrite.facets] <- stats$LD_matrices[which(names(stats$LD_matrices) %in% facets)]
        names(x@pairwise.LD$LD_matrices)[overwrite.facets] <- names(stats$LD_matrices[which(names(stats$LD_matrices) %in% facets)])
      }

      ## add
      if(length(add.facets) > 0){
        x@pairwise.LD$LD_matrices <- c(x@pairwise.LD$LD_matrices, stats$LD_matrices[add.facets])
      }
    }
  }
  return(x)
}

# subsets snpR data
subset.snpR.data <- function(x, snps = 1:nrow(x), samps = 1:ncol(x), facets = NULL, subfacets = NULL, snp.facets = NULL, snp.subfacets = NULL){
  # if subfacets or snp.subfacets were selected, figure out which samples and loci to keep
  if(!(is.null(snp.facets[1])) & !(is.null(snp.subfacets[1])) | !(is.null(facets[1])) & !(is.null(subfacets[1]))){

    # if snp.subfacets are requested
    if(!(is.null(snp.facets[1])) & !(is.null(snp.subfacets[1]))){
      t.snp.meta <- x@snp.meta

      # check for and get info on complex facets
      complex.snp.facets <- snp.facets[grep("\\.", snp.facets)]
      if(length(complex.snp.facets) > 0){
        for(i in 1:length(complex.snp.facets)){
          tfacets <- unlist(strsplit(complex.snp.facets[i], "\\."))
          tcols <- t.snp.meta[colnames(t.snp.meta) %in% tfacets]
          tcols <- tcols[,match(colnames(tcols), tfacets)]
          t.snp.meta <- cbind(t.snp.meta, do.call(paste, c(tcols, sep=".")))
        }
        colnames(t.snp.meta)[(ncol(t.snp.meta) - length(complex.snp.facets) + 1):ncol(t.snp.meta)] <- complex.snp.facets
      }

      # get the snps to keep
      t.snp.meta <- t.snp.meta[,colnames(t.snp.meta) %in% snp.facets]
      fsnps <- which(as.logical(rowSums(matrix(as.matrix(t.snp.meta) %in% snp.subfacets, nrow(x@snp.meta))))) # here's the snps to keep, those where at least one subfacet level is present in the provided snp.subfacets.
      snps <- snps[snps %in% fsnps]
    }

    # if sample subfacets are requested
    if(!(is.null(facets[1])) & !(is.null(subfacets[1]))){
      t.samp.meta <- x@sample.meta

      # check for and get info on complex facets
      complex.samp.facets <- facets[grep("\\.", facets)]
      if(length(complex.samp.facets) > 0){
        for(i in 1:length(complex.samp.facets)){
          tfacets <- unlist(strsplit(complex.samp.facets[i], "\\."))
          tcols <- t.samp.meta[colnames(t.samp.meta) %in% tfacets]
          tcols <- tcols[,match(colnames(tcols), tfacets)]
          t.samp.meta <- cbind(t.samp.meta, do.call(paste, c(tcols, sep=".")))
        }
        colnames(t.samp.meta)[(ncol(t.samp.meta) - length(complex.samp.facets) + 1):ncol(t.samp.meta)] <- complex.samp.facets
      }

      # get the samples to keep
      t.samp.meta <- t.samp.meta[,colnames(t.samp.meta) %in% facets]
      fsamps <- which(as.logical(rowSums(matrix(as.matrix(t.samp.meta) %in% subfacets, nrow(x@sample.meta))))) # here's the samples to keep, those where at least one subfacet level is present in the provided snp.subfacets.
      samps <- samps[samps %in% fsamps]
    }
  }

  # subset
  if(!identical(samps, 1:ncol(x))){
    dat <- x[snps, samps]
    dat <- import.snpR.data(dat, snp.meta = x@snp.meta[snps,], sample.meta = x@sample.meta[samps,], mDat = x@mDat)
    dat <- add.facets.snpR.data(dat, x@facets[-which(x@facets == ".base")])
    warning("Since samples were subset, any stats will need to be recalculated.\n")
    return(dat)
  }
  else{
    x <- snpRdata(.Data = x[snps,],
                  sample.meta = x@sample.meta,
                  snp.meta = x@snp.meta[snps,],
                  facet.meta = x@facet.meta[x@facet.meta$.snp.id %in% snps,],
                  mDat = x@mDat,
                  snp.form = x@snp.form,
                  geno.tables = list(gs = x@geno.tables$gs[x@facet.meta$.snp.id %in% snps,],
                                     as = x@geno.tables$as[x@facet.meta$.snp.id %in% snps,],
                                     wm = x@geno.tables$wm[x@facet.meta$.snp.id %in% snps,]),
                  ac = x@ac[x@facet.meta$.snp.id %in% snps,],
                  facets = x@facets,
                  facet.type = x@facet.type,
                  stats = x@stats[x@stats$.snp.id %in% snps,],
                  names = x@names,
                  row.names = x@row.names[snps])
    warning("Any window stats will need to be recalculated.\n")
    return(x)
  }
}

# checks requested facets, sorts, remove duplicates, and may remove a type of facet (usually snp based)
# remove.types: snp, sample, complex (both), simple (only one type)
check.snpR.facet.request <- function(x, facets, remove.type = "snp", return.type = F){
  if(any(facets == "all")){
    facets <- x@facets
  }
  if(is.null(facets)){
    facets <- ".base"
    if(return.type){
      return(list(facets, ".base"))
    }
    else{
      return(".base")
    }
  }


  # remove the facet parts as requested.
  facets <- strsplit(facets, "(?<!^)\\.", perl = T)
  to.remove <- logical(length(facets))
  missing.facets <- character(0)
  facet.types <- character(length(facets))
  for(i in 1:length(facets)){
    if(identical(facets[[i]], ".base")){next()}
    facets[[i]] <- sort(facets[[i]])
    samp.facets <- which(facets[[i]] %in% colnames(x@sample.meta))
    snp.facets <- which(facets[[i]] %in% colnames(x@snp.meta))
    missing.facets <- c(missing.facets, facets[[i]][which(!((1:length(facets[[i]])) %in% c(samp.facets, snp.facets)))])

    if(remove.type == "snp"){
      if(length(snp.facets) > 0){
        if(length(snp.facets) == length(facets[[i]])){
          to.remove[i] <- T
        }
        else{
          facets[[i]] <- facets[[i]][-snp.facets]
        }
      }
    }

    else if(remove.type == "sample"){
      if(length(samp.facets) > 0){
        if(length(samp.facets) == length(facets[[i]])){
          to.remove[i] <- T
        }
        else{
          facets[[i]] <- facets[[i]][-samp.facets]
        }
      }
    }

    else if(remove.type == "complex"){
      if(length(snp.facets) > 0 & length(samp.facets) > 0){
        to.remove[i] <- T
      }
    }

    else if(remove.type == "simple"){
      if(!(length(snp.facets) > 0 & length(samp.facets) > 0)){
        to.remove[i] <- T
      }
    }

    if(return.type){
      samp.facets <- which(facets[[i]] %in% colnames(x@sample.meta))
      snp.facets <- which(facets[[i]] %in% colnames(x@snp.meta))
      missing.facets <- c(missing.facets, facets[[i]][which(!((1:length(facets[[i]])) %in% c(samp.facets, snp.facets)))])

      if(length(samp.facets) > 0 & length(snp.facets) > 0){
        facet.types[i] <- "complex"
      }
      else if(length(samp.facets) > 0){
        facet.types[i] <- "sample"
      }
      else if(length(snp.facets) > 0){
        facet.types[i] <- "snp"
      }
    }
    facets[[i]] <- paste(facets[[i]], collapse = ".")
  }

  if(length(missing.facets) > 0){
    dups <- which(duplicated(missing.facets))
    if(length(dups) > 0){
      stop("Facet(s) ", paste(missing.facets[-dups], collapse = ", "), " not found in x metadata.\n")
    }
    stop("Facet(s) ", paste(missing.facets, collapse = ", "), " not found in x metadata.\n")
  }

  if(any(to.remove)){
    facets <- facets[-which(to.remove)]
    facet.types <- facet.types[-which(to.remove)]
  }

  # remove duplicates and return
  dups <- duplicated(facets)
  if(any(dups)){
    facets <- facets[-which(dups)]
    facet.types <- facet.types[-which(to.remove)]
  }
  if(return.type){
    return(list(unlist(facets), facet.types))
  }
  else{
      return(unlist(facets))

  }
}

#' Tabulate allele and genotype counts at each locus.
#'
#' \code{tabulate_genotypes} creates matricies containing counts of observed alleles and genotypes at each locus.
#'
#' This function is pirmarily used interally in several other funcitons, but may occasionally be useful.
#'
#' @param x Input genotype data, where columns are individuals and rows are snps. No metadata.
#' @param mDat Character string. How are missing \emph{genotypes} noted?
#' @param verbose Logical. Should the function report progress?
#'
#' @return A list of matrices. gs is the genotype matrix, as is the allele matrix, and wm is the genotype matrix with missing genotypes.
#'
#' @examples
#' tabulate_genotypes(stickSNPs[,-c(1:3)], "NN")
#'
tabulate_genotypes <- function(x, mDat, verbose = F){

  # get a genotype table
  x <- data.table::setDT(x)

  snp_form <- nchar(x[1,1])   # get information on data format
  x <- data.table::melt(t(x)) # transpose and melt

  gmat <- data.table::dcast(data.table::setDT(x), Var2 ~ value, value.var='value', length)
  gmat <- gmat[,-1]
  mis.cols <- -which(colnames(gmat) == mDat)
  tmat <- gmat[,..mis.cols] # remove missing data

  #get matrix of allele counts
  #initialize
  hs <- substr(colnames(tmat),1,snp_form/2) != substr(colnames(tmat), (snp_form/2 + 1), snp_form*2) # identify heterozygotes.
  if(verbose){cat("Getting allele table...\n")}
  as <- unique(unlist(strsplit(paste0(colnames(tmat)), "")))
  amat <- data.table::as.data.table(matrix(0, nrow(gmat), length(as)))
  colnames(amat) <- as

  #fill in
  for(i in 1:length(as)){
    as[i]
    b <- grep(as[i], colnames(tmat))
    hom <- which(colnames(tmat) == paste0(as[i], as[i]))
    if(length(hom) == 0){
      het <- b
      amat[,..i] <- rowSums(tmat[,..het])
    }
    else{
      het <- b[b != hom]
      if(length(het) > 0){
        if(data.table::is.data.table(tmat[,..het])){
          set(amat, j = i, value = (tmat[,..hom] * 2) + rowSums(tmat[,..het]))
        }
        else{
          amat[,i] <- (tmat[,hom] * 2) + tmat[,het]
        }
      }
      else{
        set(amat, j = i, value = (tmat[,..hom] * 2))
      }
    }
  }
  return(list(gs = as.matrix(tmat), as = amat, wm = as.matrix(gmat)))
}


#'Filter SNP data.
#'
#'\code{filter_snps} filters SNP data to remove loci which violate any of several assumptions and/or individuals which are sequenced at too few SNP loci.
#'
#'Description of x:
#'    Contains metadata in columns 1:ecs. Remainder of columns contain genotype calls for each individual. Each row is a different SNP, as given by format_snps output options 4 or 6.
#'
#'Possible filters:
#'\itemize{
#'    \item{maf, minor allele frequency: }{removes SNPs where the minor allele frequency is too low. Can look for mafs below #'provided either globally or search each population individually.}
#'    \item{hf_hets, high observed heterozygosity: }{removes SNPs where the observed heterozygosity is too high.}
#'    \item{min_ind, minimum individuals: }{removes SNPs that were genotyped in too few individuals.}
#'    \item{min_loci, minimum loci: }{removes individuals sequenced at too few loci.}
#'    \item{non_poly, non-polymorphic SNPs: }{removes SNPs that are not polymorphic (not true SNPs).}
#'    \item{bi_al, non-biallelic SNPs: }{removes SNPs that have more than two observed alleles.}
#'}
#'
#'Note that filtering out poorly sequenced individuals creates a possible conflict with the loci filters, since after individuals are removed, some loci may no longer pass filters. For example, if a portion of individuals in one population all carry the only instances of a rare minor allele that still passes the maf threshold, removing those individuals may cause the loci to no longer be polymorphic in the sample.
#'
#'To counter this, the "re_run" argument can be used to pass the data through a second filtering step after individuals are removed. By default, the "partial" re-run option is used, which re-runs only the non-polymorphic filter (if it was originally set), since these may cause downstream analysis errors. The "full" option re-runs all set filters. Note that re-running any of these filters may cause individuals to fail the individual filter after loci removal, and so subsequent tertiary re-running of the individual filters, followed by the loci filters, and so on, could be justified. This function stops after the second loci re-filtering, since that step is likely to be the most important to prevent downstream analytical errors.
#'
#'Via the "pop" argument, this function can filter by minor allele frequencies in either \emph{all} samples or \emph{in each population and the entire sample}. The latter should be used in instances where populaiton sizes are very different or there are \emph{many} populations, and thus common alleles of interest in one population might be otherwise filtered out. With very small populations, however, this may leave noise in the sample! In most cases, filtering the entire sample is sufficient.
#'
#' @param x data.frame. Input data, in the numeric or character format as given by format_snps options 4 or 6.
#' @param ecs Integer. Number of metadata columns at the start of x.
#' @param maf FALSE or numeric between 0 and 1, default FALSE. Minimum acceptable minor allele frequency
#' @param hf_hets FALSE or numeric between 0 and 1, default FALSE. Maximum acceptable heterozygote frequency.
#' @param min_ind FALSE or integer, default FALSE. Minimum number of individuals in which a loci must be sequenced.
#' @param min_loci FALSE or numeric between 0 and 1, default FALSE. Minimum proportion of SNPs an individual must be genotyped at.
#' @param re_run FALSE, "partial", or "full", default "partial". How should loci be re_filtered after individuals are filtered?
#' @param pop FALSE or table, default FALSE. A table with population information for individuals. Individuals must be sorted in input data in the population order given in this table.
#' @param non_poly boolean, default TRUE. Should non-polymorphic loci be removed?
#' @param bi_al boolean, default TRUE. Should non-biallelic SNPs be removed?
#' @param mDat character variable, default "NN". Format of missing \emph{genotypes}. Overall data format is infered from this. Can be either "NN" or "0000".
#' @param in.tab. FALSE or list. Option to provide tables of snp and genotype counts at each loci, used in many reformatting and filtering steps. Used internally.
#' @param out.tab. FALSE or list. Option to return tables of snp and genotype counts at each loci, used in many reformatting and filtering steps. Used internally.
#'
#' @return A data.frame in the same format as the input, with SNPs and individuals not passing the filters removed.
#'
#' @examples
#' #Strict filtering for missing individuals and unsequenced loci with partial re-run:
#' filter_snps(stickSNPs, 3, 0.05, 0.55, 250, .75)
#'
#' #Strict maf filtering with pops.
#' ##prep pop info
#' l <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' ##filter
#' filter_snps(stickSNPs, 3, 0.05, 0.55, 250, .75, "full", pop = l)
#'
filter_snps <- function(x, maf = FALSE, hf_hets = FALSE, min_ind = FALSE,
                        min_loci = FALSE, re_run = "partial", maf.facets = NULL,
                        non_poly = TRUE, bi_al = TRUE){

  #==============do sanity checks====================

  if(maf){
    if(!is.numeric(maf)){
      stop("maf must be a numeric value.")
    }
    if(length(maf) != 1){
      stop("maf must be a numeric vector of length 1.")
    }
  }

  if(hf_hets){
    if(!is.numeric(hf_hets)){
      stop("hf_hets must be a numeric value.")
    }
    if(length(hf_hets) != 1){
      stop("hf_hets must be a numeric vector of length 1.")
    }
  }

  if(min_ind){
    if(!is.numeric(min_ind)){
      stop("min_ind must be a numeric value.")
    }
    if(length(min_ind) != 1){
      stop("min_ind must be a numeric vector of length 1.")
    }
    if(min_ind > 1 | min_ind < 0){
      stop("min_ind is the minimum proportion of sequenced individuals, and so must be between 0 and 1.\n")
    }
  }

  if(min_loci){
    if(!is.numeric(min_loci) | (min_loci <= 0 | min_loci >= 1) | length(min_loci) != 1){
      stop("min_loci must be a numeric value between but not equal to 0 and 1.")
    }
  }

  if(re_run != FALSE){
    if(re_run != "partial" & re_run != "full"){
      cat("re_run must be set to partial or full if not FALSE.\n")
    }
  }

  if(!is.null(maf.facets[1])){
    # add any needed facets...
    miss.facets <- maf.facets[which(!(maf.facets %in% x@facets))]
    if(length(miss.facets) != 0){
      cat("Adding missing facets...\n")
      # need to fix any multivariate facets (those with a .)
      x <- add.facets.snpR.data(x, miss.facets)
    }

    # check for bad facets to remove (those that don't just consider samples)
    if(any(x@facet.type[x@facets %in% maf.facets] != "sample")){
      vio.facets <- x@facets[match(maf.facets, x@facets)]
      vio.facets <- vio.facets[which(x@facet.type[x@facets %in% maf.facets] != "sample")]
      warning(paste0("Facets over which to maf.filter must be sample specific facets, not snp specific facets! Removing non-sample facets: \n", paste0(vio.facets, collapse = " "), ".\n"))
      maf.facets <- maf.facets[-which(maf.facets %in% vio.facets)]
    }
  }

  #==============set up, get values used later, clean up data a bit,define subfunctions==========
  cat("Initializing...\n")

  #get headers
  headers <- x@snp.meta
  snp_form <- x@snp.form
  mDat <- x@mDat


  #function to filter by loci, to be called before and after min ind filtering (if that is requested.)
  filt_by_loci <- function(){
    #==========================run filters========================
    vio.snps <- logical(nrow(x)) #vector to track status

    amat <- x@geno.tables$as[x@facet.meta$facet == ".base",]
    gmat <- x@geno.tables$gs[x@facet.meta$facet == ".base",]
    wmat <- x@geno.tables$wm[x@facet.meta$facet == ".base",]

    # non-biallelic and non-polymorphic loci
    if(bi_al | non_poly){
      bimat <- ifelse(amat, TRUE, FALSE)

      if(bi_al){
        cat("Filtering non-biallelic loci...\n")
        bi <- ifelse(rowSums(bimat) > 2, T, F) # if true, should keep the allele
        bi.vio <- x@snp.meta$.snp.id[which(bi)] # IDs of violating snps.
        cat(paste0("\t", length(bi.vio), " bad loci\n"))
        vio.snps[which(bi)] <- T
      }

      if(non_poly){
        cat("Filtering non_polymorphic loci...\n")
        np <- ifelse(rowSums(bimat) < 2, T, F) # if true, should keep the allele
        np.vio <- x@snp.meta$.snp.id[which(np)]
        cat(paste0("\t", length(np.vio), " bad loci\n"))
        vio.snps[which(np)] <- T
      }

      #some tests require this, so subset the matrices and redefine things if true and some are multi-allelic
      # if((maf | hf_hets) &  sum(bi) != 0){
      #   x <- x[bi,]
      #   xv <-   xv <- as.vector(t(x))
      #   gs <- unique(xv)
      #   hs <- substr(gs,1,snp_form/2) != substr(gs, (snp_form/2 + 1), snp_form*2)
      #   mpos <- which(gs == mDat)
      #   as <- unique(c(substr(gs,1,snp_form/2), substr(gs, (snp_form/2 + 1), snp_form*2)))
      #   as <- as[as != substr(mDat, 1, snp_form/2)] #that aren't N?
      #   gmat <- gmat[bi,]
      #   tmat <- tmat[bi,]
      #   amat <- amat[bi,]
      #   headers <- headers[bi,]
      # }
      # else{
      #   keep <- keep + bi
      # }
    }

    #========min inds=======
    if(min_ind){
      cat("Filtering loci sequenced in few individuals...\n")
      mi <- wmat[,colnames(wmat) == mDat]
      mi <- (nrow(x@sample.meta) - mi)/nrow(x@sample.meta) < min_ind
      mi.vio <- x@snp.meta$.snp.id[which(mi)]
      cat(paste0("\t", length(mi.vio), " bad loci\n"))
    }

    #========minor allele frequency, both total and by pop. Should only run if bi_al = TRUE.=========
    if(maf){
      #if not filtering with multiple pops
      if(is.null(maf.facets)){
        cat("Filtering low minor allele frequencies, no pops...\n")

        # check to see if we need to calculate mafs:
        if(any(colnames(x@stats) == "maf")){ # check that mafs have been calculated, the all facet must exist
          if(any(is.na(x@stats$maf[x@stats$facet == ".base"]))){ # check that mafs have been calculated for the all facet
            mafs <- 1 - matrixStats::rowMaxs(amat)/rowSums(amat)
          }
          else{
            mafs <- x@stats$maf[x@stats$facet == ".base"]
          }
        }
        else{
          mafs <- 1 - matrixStats::rowMaxs(amat)/rowSums(amat)
        }


        mafs <- mafs < maf #less than required, set to true and reject.
        mafs[is.na(mafs)] <- TRUE
        maf.vio <- x@snp.meta$.snp.id[which(mafs)]
        vio.snps[which(mafs)] <- T
      }
      else{
        cat("Filtering low minor allele frequencies by facet.\n")
        # pmafs <- logical(nrow(x))

        # see if we need to calculate mafs
        if(any(colnames(x@stats) == "maf")){ # mafs have been caluclated
          # get mafs for any uncalculated facets

          # if mafs have been calculated, but not for of the requested any facets...
          if(!any(x@stats$facet %in% maf.facets)){
            x <- calc_maf(x, facets = maf.facets)
          }
          #if mafs have been calculated for some but not all of our facets
          else if(any(is.na(x@stats$maf[x@stats$facet %in% maf.facets]))){
            run.facets <- unique(x@stats$facet[which(is.na(x@stats$maf[x@stats$facet %in% maf.facets]))])
            x <- calc_maf(x, facets = run.facets)
          }
        }
        else{
          x <- calc_maf(x, facets = maf.facets)
        }

        # grab mafs
        mafs <- x@stats[x@stats$facet %in% maf.facets,]
        mafs$maf <- mafs$maf < maf

        # now, figure out in how many subfacets the maf is too low. If all, the loci violates the filter
        cmafs <- reshape2::dcast(mafs[,c(".snp.id", "maf")], fun.aggregate = sum,
                                 formula = ... ~ "maf", value.var = "maf")

        # add in the overall maf, since differential fixation would otherwise be removed.
        if(any(is.na(x@stats$maf[x@stats$facet == ".base"]))){ # check that mafs have been calculated for the all facet
          a.mafs <- 1 - matrixStats::rowMaxs(amat)/rowSums(amat)
        }
        else{
          a.mafs <- x@stats$maf[x@stats$facet == ".base"]
        }

        a.mafs <- a.mafs < maf
        cmafs$maf <- cmafs$maf + a.mafs

        # check vio and report
        maf.vio <- x@snp.meta$.snp.id[which(cmafs$maf == (1 + length(unique(mafs$subfacet))))]
        cat(paste0("\t", length(maf.vio), " bad loci\n"))
        vio.snps[which(cmafs$maf == (1 + length(unique(mafs$subfacet))))] <- T
      }
    }

    #========hf_hets. Should only run if bi_al = TRUE.==========
    # working here
    if(hf_hets){
      cat("Filtering high frequency heterozygote loci...\n")

      # get heterozygote frequency
      hs <- which(substr(colnames(gmat), 1, snp_form/2) != substr(colnames(gmat), (snp_form/2) + 1, snp_form))
      het_f <- rowSums(gmat[,hs])/rowSums(gmat)

      # check violation
      het_f <- het_f > hf_hets #if false, heterozygote frequency is lower than cut-off, keep locus
      het_f.vio <- x@snp.meta$.snp.id[which(het_f)]
      cat(paste0("\t", length(het_f.vio), " bad loci\n"))
      vio.snps[which(het_f)] <- T
    }


    #==========remove violating loci==================
    if(any(vio.snps)){
      vio.ids <- x@snp.meta$.snp.id[which(vio.snps)]
      ngs <- x@geno.tables$gs[-which(x@facet.meta$.snp.id %in% x@snp.meta$.snp.id[which(vio.snps)]),]
      nas <- x@geno.tables$as[-which(x@facet.meta$.snp.id %in% x@snp.meta$.snp.id[which(vio.snps)]),]
      nwm <- x@geno.tables$wm[-which(x@facet.meta$.snp.id %in% x@snp.meta$.snp.id[which(vio.snps)]),]
      ngs <- list(gs = ngs, as = nas, wm = nwm)
      rm(nas, nwm)

      x <- snpRdata(.Data = x[-which(vio.snps),],
                    sample.meta = x@sample.meta,
                    snp.meta = x@snp.meta[-which(vio.snps),],
                    facet.meta = x@facet.meta[-which(x@facet.meta$.snp.id %in% vio.ids),],
                    geno.tables = ngs,
                    ac = x@ac[-which(x@facet.meta$.snp.id %in% vio.ids),],
                    stats = x@stats[-which(x@stats$.snp.id %in% vio.ids),],
                    window.stats = x@window.stats,
                    facets = x@facets,
                    facet.type = x@facet.type,
                    row.names = x@row.names[-which(vio.snps)])
    }
    return(x)
  }

  #funciton to filter by individuals.
  min_loci_filt <- function(){
    cat("Filtering out individuals sequenced in few kept loci...\n")
    mcounts <- colSums(ifelse(x == mDat, 1, 0))
    rejects <- which(mcounts/nrow(x) >= (1 - min_loci))
    if(length(rejects) > 0){
      old.facets <- x@facets
      invisible(capture.output(x <- import.snpR.data(x[,-rejects],
                                                     snp.meta = x@snp.meta,
                                                     sample.meta = x@sample.meta[-rejects,],
                                                     mDat = mDat)))
      cat("Re-calculating and adding facets.\n")
      x <- add.facets.snpR.data(x, old.facets[-which(old.facets == ".base")])
      warning("Any calculated stats will be removed, since individuals were filtered out!\n")
    }
    return(list(x = x, rejects = rejects))
  }

  #==========================call the functions as requested.==================
  if(any(c(non_poly, bi_al, maf, hf_hets, min_ind) != FALSE)){
    cat("Filtering loci. Starting loci:", nrow(x), "\n")

    # run the filter
    x <- filt_by_loci()

    if(nrow(x) == 1){
      stop("No loci remain after filters.")
    }

    cat("\tEnding loci:", nrow(x), "\n")
  }

  # run the minimum sequenced loci filter
  if(min_loci){
    cat("Filtering individuals. Starting individuals:", ncol(x), "\n")
    x <- min_loci_filt()
    if(length(x$rejects) == 0){
      cat("No individuals removed.\n")
    }
    else{
      x <- x$x
      cat("\tEnding individuals:", ncol(x), "\n")
      if(re_run != FALSE){
        cat("Re-filtering loci...\n")

        if(re_run == "partial"){
          maf <- FALSE
          hf_hets <- FALSE
          min_ind <- FALSE
          bi_al <- FALSE
        }
        if(any(c(non_poly, bi_al, maf, hf_hets, min_ind) != FALSE)){
          x <- filt_by_loci() # re-filter loci to make sure that we don't have any surprise non-polys ect.
          cat("\tFinal loci count:", nrow(x), "\n")
        }
        else{
          cat("\tNo variables to re-fitler.\n")
        }
      }
    }
  }

  # return results
  return(x)
}

#'Re-format genomic data.
#'
#'\code{format_snps} reformats SNP data into a range of different possible formats for use in snpR functions and elsewhere. Supports microsatellite data for conversion to a few formats as well.
#'
#'
#'Format options:
#'\itemize{
#'    \item{ac: }{allele count format, allele counts tabulated for all samples or within populations.}
#'    \item{genepop: }{genepop format, genotypes stored as four numeric characters (e.g. "0101", "0204"), transposed, and formatted for genepop. Rownames are individual IDs in genepop format, colnames are SNP ids, matching the first metadata column in input.}
#'    \item{structure: }{STRUCTURE format, two lines per individual: allele calls stored as single character numeric (e.g. "1", "2"). Allele calls per individual stored on two subsequent lines.}
#'    \item{numeric: }{numeric genotype tab format, genotypes stored as four numeric characters (e.g. "0101", "0204").}
#'    \item{hapmap: }{Migrate-n hapmap, allele counts tabulated within populations, in migrate-n hapmap format. Since this migrate-n implementation is iffy, this probably shouldn't be used much.}
#'    \item{character: }{character genotype tab format, genotypes stored as actual base calls (e.g. "AA", "CT").}
#'    \item{pa: }{allele presence/absence format, presence or absence of each possible allele at each possible genotype noted. Interpolation possible, with missing data substituted with allele freqency in all samples or each population.}
#'    \item{rafm: }{RAFM format, two allele calls at each locus stored in subsequent columns, e.g. locus1.1 locus1.2.}
#'    \item{faststructure: }{fastSTRCTURE format, identical to STRUCTURE format save with the addition of filler columns proceeding data such that exactly 6 columns proceed data. These columns can be filled with metadata if desired.}
#'    \item{dadi: }{dadi format SNP data format, requires two columns named "ref" and "anc" with the flanking bases around the SNP, e.g. "ACT" where the middle location is the A/C snp.}
#'    \item{plink: }{PLINK! binary input format, requires columns named "group", "snp", and "position", and may contain a column named "cM", "cm", or "morgans", containing linkage group/chr, snp ID, position in bp, and distance in cM in order to create .bim extended map file.}
#'}
#'
#'Example datasets in each format are available in \code{\link{stickFORMATs}} in elements named for output options.
#'
#'Input formats: For now, all input formats require at least two metadata columns.
#'\itemize{
#'    \item{NN: }{SNP genotypes stored as actual base calls (e.g. "AA", "CT").}
#'    \item{0000: }{SNP genotypes stored as four numeric characters (e.g. "0101", "0204").}
#'    \item{msat_2: }{Microsatellite genotypes stored as four numeric characters (e.g. "0101", "2740").}
#'    \item{msat_3: }{Microsatellite genotypes stored as six numeric characters (e.g. "120123", "233235").}
#'    \item{snp_tab: }{SNP genotypes stored with genotypes in each cell, but only a single nucleotide noted if homozygote and two nucleotides seperated by a space if heterozygote (e.g. "T", "T G").}
#'    \item{sn: }{SNP genotypes stored with genotypes in each cell as 0 (homozyogous allele 1), 1 (heterozygous), or 2 (homozyogus allele 2).}
#'}
#'
#'Currently, msat_2 and msat_3 only support conversion to output option 7. 2, 3, and 4 are forthcoming.
#'
#'
#' @param x data.frame. Input data, in any of the above listed input formats.
#' @param ecs Integer. Number of extra metadata columns at the start of x.
#' @param output Character. Which of the output format should be used?
#' @param input_form Character, default "NN". Which of the above input formats should be used (e.g. "NN", "msat_2")?
#' @param mDat Character, default "NN". The coding for missing \emph{genotypes} in x (typically "NN" or "0000").
#' @param pop FALSE or table, default FALSE. A table with population information for individuals. Individuals must be sorted in input data in the population order given in this table.
#' @param FALSE or table, default FALSE. A table with population information for individuals. Individuals must be sorted in input data in the population order given in this table.
#' @param n_samp Integer or numeric vector, default NA. For output option 3. How many random loci should be selected? Can either be an integer or a numeric vector of loci to use.
#' @param interp_miss boolean, default TRUE. For output option 7. Should missing data be interpolated? Needed for PCA, ect.
#' @param lnames character vector, default NULL. For output option 7, optional vector of locus names by which to name output columns. If not provided, will use 1:nrow(x).
#' @param outfile character vector, default FALSE. If provided, the path to the output file to write to.
#' @param ped data.frame, default NULL. If provided, the six column header for plink .ped files.
#' @param in.tab FALSE or list. Option to provide tables of snp and genotype counts at each loci, used in many reformatting and filtering steps. Used internally.
#' @param out.tab FALSE or list. Option to return tables of snp and genotype counts at each loci, used in many reformatting and filtering steps. Used internally.
#'
#' @return A data.frame in the specified format.
#'
#' @examples
#' #allele count with pops:
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' format_snps(stickSNPs, 3, "ac", pop = pops)
#'
#' #genepop:
#' format_snps(stickSNPs, 3, "genepop")
#'
#' #STRUCTURE, subsetting out 100 random alleles:
#' format_snps(stickSNPs, 3, "structure", n_samp = 100)
#'
#' #STRUCTURE, subseting out the first 100 alleles:
#' format_snps(stickSNPs, 3, "structure", n_samp = 1:100)
#'
#' #fastSTRUCTURE
#' format_snps(stickSNPs, 3, "faststructure")
#'
#' #numeric:
#' format_snps(stickSNPs, 3, "numeric")
#'
#' #hapmap for migrate-n:
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' format_snps(stickSNPs, 3, "hapmap", pop = pops)
#'
#' #character:
#' num <- format_snps(stickSNPs, 3, 4)
#' format_snps(num, 3, "character", input_form = "0000", miss = "00")
#'
#' #presence/absence, SNP data:
#' format_snps(stickSNPs, 3, "pa")
#'
#' #presence/absence, 3 character microsat data (2 character is very similar):
#' format_snps(sthMSATS[seq(1, 13, by = 4),], 3, "pa", input_form = "msat_3", miss = "000")
#'
#' #RAFM, taking only 100 random snps.
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' format_snps(stickSNPs, 3, "rafm", pop = pops, n_samp = 100)
#'
#' #dadi
#' pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
#' format_snps(cbind(ref = "ATA", anc = "ACT", stickSNPs), 5, "dadi", pop = pops)
#'
#' #PLINK! format
#' format_snps(stickSNPs, 3, "plink", outfile = "plink_out")
#' #from command line, then run plink_out.sh to generate plink_out.bed.
#'
#' #PLINK! format with provided ped
#' ped <- data.frame(fam = c(rep(1, 210), rep("FAM2", 210)), ind = 1:420, mat = 1:420, pat = 1:420, sex = sample(1:2, 420, T), pheno = sample(1:2, 420, T))
#' format_snps(stickSNPs, 3, "plink", outfile = "plink_out", ped = ped)
#' #from command line, then run plink_out.sh to generate plink_out.bed.
#'
format_snps <- function(x, output = "ac", facets = NULL, n_samp = NA,
                        interp_miss = T, outfile = FALSE,
                        ped = NULL, input_format = NULL,
                        input_meta_columns = NULL, input_mDat = NULL,
                        sample.meta = NULL, snp.meta = NULL){

  #======================sanity checks================
  # check that a useable output format is given.
  output <- tolower(output) # convert to lower case.
  pos_outs <- c("ac", "genepop", "structure", "0000", "hapmap", "NN", "pa",
                "rafm", "faststructure", "dadi", "plink", "sn", "snprdata")
  if(!(output %in% pos_outs)){
    stop("Unaccepted output format specified. Check documentation for options.\n")
  }

  # check that provided snpRdata objects are in the correct format
  if(is.null(input_format)){
    if(class(x) != "snpRdata"){
      stop("If x is not a snpRdata object, provide input data format.\n")
    }
  }

  # do checks, print info
  facet.types <- x@facet.type[match(facets, x@facets)]
  snp.facets <- which(facet.types == "snp")
  both.facets <- which(facet.types == "both")
  sample.facets <- which(facet.types == "sample")

  if(output == "ac"){
    cat("Converting to allele count format.\n")
  }
  else if(output == "genepop"){
    cat("Converting to genepop format.\n")
    if(length(c(both.facets, snp.facets)) != 0){
      warning("Removing invalid facet types (snp or snp and sample specific).\n")
      facets <- facets[-c(snp.facets, both.facets)]
    }
  }
  else if(output == "structure" | output == "faststructure"){
    if(output == "structure"){cat("Converting to STRUCTURE format.\n")}
    else{cat("Converting to fastSTRUCTURE format.\n")}
    if(length(n_samp) > 1){
      if(is.integer(n_samp)){
        cat("Number of designated sub-samples to take:", length(n_samp), "\n")
      }
      else{
        stop("Number of sub-samples to take must be a positive integer vector.\n")
      }
    }
    else if (is.numeric(n_samp)){
      if(floor(n_samp) == n_samp & n_samp > 0){
        cat("Number of designated sub-samples to take:", n_samp, "\n")
      }
      else{
        stop("Number of sub-samples to take must be a positive integer vector.\n")
      }
    }
    else if (!is.na(n_samp)){
      stop("Number of sub-samples to take must be a positive integer vector.\n")
    }
    if(length(c(both.facets, snp.facets)) != 0){
      warning("Removing invalid facet types (snp or snp and sample specific).\n")
      facets <- facets[-c(snp.facets, both.facets)]
    }
    if(output == "faststructure"){
      if(length(facets) > 5){
        stop("Too many facets selected for fastSTRUCTURE. Limit 5.\n")
      }
    }
  }
  else if(output == "0000"){
    cat("Converting to numeric 2 character format.\n")
  }
  else if(output == "hapmap"){
    cat("Converting to migrate-N hapmap format.\n")
  }
  else if(output == "NN"){
    cat("Converting to NN format.\n")
  }
  else if(output == "pa"){
    cat("Converting to allele presence/absense format.\n")
    if(length(c(both.facets, snp.facets)) != 0){
      warning("Removing invalid facet types (snp or snp and sample specific).\n")
      facets <- facets[-c(snp.facets, both.facets)]
    }
  }
  else if(output == "rafm"){
    cat("Converting to RAFM format.\n")
    if(length(c(both.facets, snp.facets)) != 0){
      warning("Removing invalid facet types (snp or snp and sample specific).\n")
      facets <- facets[-c(snp.facets, both.facets)]
    }
  }
  else if(output == "dadi"){
    if(is.null(input_format)){
      if(!("ref" %in% colnames(x@snp.meta)) | !("anc" %in% colnames(x@snp.meta))){
        stop("Reference and ancestor/outgroup flanking bases required in snp metadata columns named 'ref' and 'anc', respecitvely")
      }
    }
    else {
      if(!("ref" %in% colnames(x)) | !("anc" %in% colnames(x))){
        stop("Reference and ancestor/outgroup flanking bases required in columns named 'ref' and 'anc', respecitvely")
      }
    }
    cat("Converting to dadi format...\n")
    if(length(c(both.facets, snp.facets)) != 0){
      warning("Removing invalid facet types (snp or snp and sample specific).\n")
      facets <- facets[-c(snp.facets, both.facets)]
    }
  }
  else if(output == "plink"){
    if(is.null(input_format)){
      if(!all(c("position") %in% colnames(x@snp.meta)) | !any(c("group", "chr") %in% colnames(x@snp.meta))){
        stop("Columns named position, group/chr, and snp containing position in bp, chr/linkage group/scaffold, and snp ID must be present in snp metadata!")
      }
    }
    else if(!all(c("position") %in% colnames(x)) | !any(c("group", "chr") %in% colnames(x@snp.meta))){
      stop("Columns named position, group, and snp containing position in bp, chr/linkage group/scaffold, and snp ID must be present in x!")
    }
    if(!is.null(ped)){
      if(!is.data.frame(ped)){
        stop("ped must be a six column data frame containg Family ID, Individual ID, Paternal ID, Maternal ID, Sex, and Phenotype and one row per sample. See plink documentation.\n")
      }
      if(ncol(ped) != 6 | nrow(ped) != ncol(data)){
        stop("ped must be a six column data frame containg Family ID, Individual ID, Paternal ID, Maternal ID, Sex, and Phenotype and one row per sample. See plink documentation.\n")
      }
    }
    cat("Converting to PLINK! binary format.")
    if(length(c(both.facets, snp.facets)) != 0){
      warning("Removing invalid facet types (snp or snp and sample specific).\n")
      facets <- facets[-c(snp.facets, both.facets)]
    }
    if(any(colnames(x@snp.meta) == "chr")){
      colnames(x@snp.meta)[which(colnames(x@snp.meta) == "chr")] <- "group"
    }
  }
  else if(output == "sn"){
    cat("Converting to single character numeric format.\n")
    if(length(c(both.facets, snp.facets)) != 0){
      warning("Removing invalid facet types (snp or snp and sample specific).\n")
      facets <- facets[-c(snp.facets, both.facets)]
    }
  }
  else if(output == "snprdata"){
    if(is.null(input_format)){
      stop("Data already in snpRdata object.\n")
    }
    if(all(is.null(snp.meta), is.null(input_meta_columns)) | is.null(input_mDat)){
      stop("sample meta, snp meta, and input metadata must be provided for conversion to snpRdata object.")
    }
    else if(is.null(snp.meta) & ! is.null(input_meta_columns)){
      cat("Using input metadata columns as snp meta.\n")
    }
    cat("Converting to snpRdata object.\n")
  }

  else{
    stop("Please specify output format.")
  }

  #======================put data into snpRdata object if not in that format to start with================
  if(!is.null(input_format)){
    cat("Converting data to snpRdata, NN format.\n")
    if(!is.null(input_meta_columns)){
      headers <- x[,c(1:input_meta_columns)]
      x <- x[,-c(1:input_meta_columns)]
    }
    else{
      headers <- snp.meta
    }
    if(is.null(sample.meta)){
      sample.meta <- data.frame(ID = 1:ncol(x))
    }

    #====================checks===================
    if(is.null(input_mDat)){
      stop("Input missing data encoding required.\n")
    }

    if(input_format == "NN"){
      cat("Input format: NN\n")
      if(nchar(input_mDat) != 2){
        stop("Missing data format must be two character.\n")
      }
      else{
        cat("Missing data format: ", input_mDat, "\n")
      }
    }
    else if(input_format == "0000"){
      cat("Input format: 0000\n")
      if(nchar(input_mDat) != 4){
        stop("Missing data format must be four characters.\n")
      }
      else{
        cat("Missing data format: ", input_mDat, "\n")
      }
    }
    else if(input_format == "snp_tab"){
      cat("Input format: snp_tab.\n")
      if(nchar(input_mDat) != 2){
        stop("Missing data format must be two characters.\n")
      }
      else{
        cat("Missing data format: ", input_mDat, "\n")
      }
    }
    else if(input_format == "sn"){
      cat("Imput format: sn\n")
      if(input_mDat %in% c(0:2)){
        stop("Missing data format must be other than 0, 1, or 2.\n")
      }
      if(nchar(input_mDat) != 1){
        stop("Missing data format must be a single character.\n")
      }
      else{
        cat("Missing data format: ", input_mDat, "\n")
      }
    }
    else{
      stop("Unsupported input format.")
    }

    #====================convert to snpRdata intermediate=============
    # 0000 put into snpRdata
    if(input_format == "0000"){
      #vectorize and replace
      xv <- as.matrix(x)
      xv <- gsub("01", "A", xv)
      xv <- gsub("02", "C", xv)
      xv <- gsub("03", "G", xv)
      xv <- gsub("04", "T", xv)
      xv <- gsub(substr(input_mDat, 1, 2), "N", xv)
      x <- as.data.frame(xv, stringsAsFactors = F)
      rm(xv)

      input_mDat <- "NN"

      #rebind to matrix and remake data.
      if(output == "NN"){
        rdata <- cbind(header, x) #all done if just converting to NN
      }
      else{
        x <- import.snpR.data(x, sample.meta = sample.meta, snp.meta = headers, mDat = "NN")
      }
    }

    # convert snp_tab to snpRdata
    if(input_format == "snp_tab"){
      xv <- as.matrix(x)
      ptf <- nchar(xv) == 1
      xv[ptf] <- paste0(xv[ptf], xv[ptf]) #double up homozygotes
      remove(ptf)
      xv <- gsub(" ", "", xv) #combine heterozygotes.
      xv[xv == input_mDat] <- "NN" #replace with correct missing data
      x <- as.data.frame(xv, stringsAsFactors = F)
      if(output == "NN"){
        rdata <- cbind(header, x)
      }
      else{
        x <- import.snpR.data(x, sample.meta = sample.meta, snp.meta = headers, mDat = "NN")
      }
    }

    #convert "sn" format to NN
    if(input_format == "sn"){
      # save for plink if that's the destination format!
      if(output == "plink"){
        pl_0g_dat <- x
      }
      xv <- as.matrix(x)
      xv[xv == 0] <- "AA"
      xv[xv == 1] <- "AC"
      xv[xv == 2] <- "CC"
      xv[xv == mDat] <- "NN"
      x <- as.data.frame(xv, stringsAsFactors = F)

      if(output == "NN"){
        rdata <- x #all done if just converting to NN
      }
      else{
        x <- import.snpR.data(x, sample.meta = sample.meta, snp.meta = headers, mDat = "NN")
      }
    }

    if(input_format == "NN"){
      x <- import.snpR.data(x, sample.meta = sample.meta, snp.meta = headers, mDat = "NN")
      if(output == "NN"){
        rdata <- x
      }
    }

    if(!is.null(facets)){
      x <- add.facets.snpR.data(x, facets)
    }

    # if this is the desired output, we're done.
    if(output == "snprdata"){
      return(x)
    }
  }

  #======================prepare outfile==================================================================
  if(outfile != FALSE){
    if(is.character(outfile) & length(outfile) == 1){
      cat("Printing results to:", outfile, "\n")
      if(file.exists(outfile)){
        #ask for confirmation before proceeding, since shit will be overwritten.
        cat("Outfile already exits. ")
        resp <- "empty"
        while(resp != "y" & resp != "n"){
          cat("Overwrite? (y or n)\n")
          resp <- readLines(n = 1)
        }
        if(resp == "n"){
          stop("Please provide acceptable path to file for output.\n")
        }
        else{
          cat("\tProceeding with conversion...\n")
        }
      }
    }
    else{
      stop("Outfile must be a character vector of length 1.\n")
    }
  }

  #======================do conversions===================================================================

  # add missing facets
  if(length(facets) == 0){facets <- NULL}

  if(!is.null(facets[1])){
    if(!all(facets %in% x@facets)){
      cat("Adding missing facets.\n")
      new.facets <- facets[which(!(facets %in% x@facets))]
      x <- add.facets.snpR.data(x, new.facets)
    }
  }

  # set facets to ".base" if not requested
  if(is.null(facets[1])){
    facets <- ".base"
  }

  #============convert to allele count, migrate-n, or dadi format. Migrate-n should ALWAYS have multiple pops (why else would you use it?)===============
  if(output == "ac" | output == "hapmap" | output == "dadi"){
    if(output == "hapmap"){cat("WARNING: Data does not have header or pop spacer rows.\n")}

    #=========basic ac constructor function, to be on each facet==========
    # x: stats slot of a snpR object, filtered to the desired facets
    # maj: major allele identities across all facets at the relevent snps
    # min: minor allele identities across all facets at the relevent snps
    # mis.al: missing allele coding
    get.ac <- function(x, maj, min, mis.al){
      # initialize:
      out <- data.frame(n_total = numeric(nrow(x)),
                        n_alleles = numeric(nrow(x)),
                        ni1 = numeric(nrow(x)),
                        ni2 = numeric(nrow(x)))


      # get the column from as matching the target allele.
      maj.col.match <- match(maj, colnames(x))
      out$ni1 <- t(x)[maj.col.match + seq(0, length(x) - ncol(x), by = ncol(x))]

      # ni1 is the rowsums minus this
      out$ni2 <- rowSums(x) - out$ni1
      out$n_total <- rowSums(x)
      out$n_alleles <- rowSums(ifelse(out[,3:4] != 0, 1, 0))
      out[is.na(out)] <- 0

      return(out)
    }

    #=========apply to requested facets=======
    # get missing maf info
    if(!any(colnames(x@stats) == "maf")){
      if(identical(".base", facets)){
        x <- calc_maf(x, facets = c(facets))
      }
      else{
        x <- calc_maf(x, facets = c(".base", facets))
      }
    }
    else if(any(is.na(x@stats$maf[x@stats$facet %in% c(".base", facets)]))){
      miss.facets <- which(c(".base", facets) %in% unique(x@stats[is.na(x@stats$maf[x@stats$facet %in% c(".base", facets)]),]$facet))
      x <- calc_maf(x, facets = c(".base", facets)[miss.facets])
    }

    # grab the major and minors for each snp in the analysis.
    x@stats <- dplyr::arrange(x@stats, .snp.id, facet, subfacet)
    maj <- x@stats$major[x@stats$facet == ".base"]
    min <- x@stats$minor[x@stats$facet == ".base"]
    maj <- rep(maj, each = length(unique(x@facet.meta[x@facet.meta$facet %in% facets,]$subfacet)))
    min <- rep(min, each = length(unique(x@facet.meta[x@facet.meta$facet %in% facets,]$subfacet)))


    # get the ac file
    ac.dat <- get.ac(x@geno.tables$as[x@facet.meta$facet %in% facets,],
                     maj = maj,
                     min = min,
                     substr(x@mDat, 1, nchar(x@mDat)/2))

    if(output == "dadi"){
      rdata <- cbind(x@stats[x@stats$facet %in% facets, c("facet", "subfacet", ".snp.id")], ac.dat)
      ni1 <- reshape2::dcast(rdata, facet + .snp.id ~ subfacet, value.var = c("ni1"))
      ni2 <- reshape2::dcast(rdata, facet + .snp.id ~ subfacet, value.var = c("ni2"))

      rdata <- cbind(ref = x@snp.meta$ref, # since everything is sorted by .snp.id, this will match.
                     anc = x@snp.meta$anc,
                     Allele1 = x@stats[x@stats$facet == ".base", "major"],
                     ni1[order(ni1$.snp.id), 3:ncol(ni1)],
                     Allele2 = x@stats[x@stats$facet == ".base", "minor"],
                     ni2[order(ni2$.snp.id), 3:ncol(ni2)],
                     x@snp.meta[,which(!(colnames(x@snp.meta) %in% c("ref", "anc", ".snp.id")))])
    }
    else{
      rdata <- cbind(x@facet.meta[x@facet.meta$facet %in% facets,],
                     ac.dat)
      if(output == "hapmap"){
        rdata <- dplyr::arrange(rdata, facet, subfacet, .snp.id)
      }
    }
  }



  ##convert to genepop or numeric format (v)
  if (output == "genepop" | output == "0000"){
    #vectorize and replace
    xv <- as.matrix(x)
    xv <- gsub("A", "01", xv)
    xv <- gsub("C", "02", xv)
    xv <- gsub("G", "03", xv)
    xv <- gsub("T", "04", xv)
    xv <- gsub(substr(x@mDat, 1, 2), "0000", xv)

    if(output == "genepop"){ #convert to genepop
      rdata <- as.data.frame(t(xv), stringsAsFactors = F) #remove extra columns and transpose data
      row.names(rdata) <- paste0(row.names(rdata), " ,") #adding space and comma to row names, as required.
    }
    else {#prepare numeric output, otherwise same format
      rdata <- as.data.frame(xv, stringsAsFactors = F)
    }
  }


  ##convert to structure, fastStructure or RAFM format (v)
  if (output == "structure" | output == "rafm" | output == "faststructure"){
    #subset if requested
    if(all(!is.na(n_samp))){
      cat("Subsetting ")
      if(length(n_samp) > 1){
        cat("designated SNPs.\n")
        x <- subset.snpR.data(x, n_samp)
      }
      else{
        cat(n_samp, " random SNPs.\n")
        x <- subset.snpR.data(x, sample(1:nrow(x)))
      }
    }

    # transpose, since these are sample based, then gsub
    xv <- t(as.matrix(x))
    xv <- gsub("A", 1, xv)
    xv <- gsub("C", 2, xv)
    xv <- gsub("G", 3, xv)
    xv <- gsub("T", 4, xv)
    xv <- gsub(substr(x@mDat, 1, nchar(x@mDat)/2), 0, xv)

    #split into two matrices
    xv1 <- substr(xv, 1, 1)
    xv2 <- substr(xv, 2, 2)

    #create output matrix and bind the data to it, structure format
    if(output == "structure" | output == "faststructure"){
      outm <- matrix(NA, nrow = 2*(ncol(x)), ncol =  nrow(x))

      #fill
      outm[seq(1,nrow(outm),2),] <- xv1
      outm[seq(2,nrow(outm),2),] <- xv2

      #add sample names
      snames <- character(nrow(outm))
      snames[seq(1,nrow(outm),2)] <- x@names
      snames[seq(2,nrow(outm),2)] <- x@names
      if(output == "structure"){ # bind sample and metadata for structure
        rdata <- cbind(ind = snames,
                       x@sample.meta[,colnames(x@sample.meta) %in% facets],
                       as.data.frame(outm, stringsAsFactors = F))
      }
      else{ #add a bunch of filler columns for faststructure and change missing data to -9
        outm[outm == 0] <- -9
        if(length(facets) == 1 & facets[1] == ".base"){
          rdata <- cbind(ind = snames,
                         matrix("filler", nrow(outm), 5),
                         as.data.frame(outm, stringsAsFactors = F))
          colnames(rdata)[2:6] <- paste0("filler", 1:5)
        }
        else{
          n.valid.facets <- sum(colnames(x@sample.meta) %in% facets)
          rdata <- cbind(ind = snames,
                         x@sample.meta[,colnames(x@sample.meta) %in% facets],
                         matrix("filler", nrow(outm), 5 - n.valid.facets),
                         as.data.frame(outm, stringsAsFactors = F))
          colnames(rdata)[(1 + n.valid.facets):6] <- paste0("filler", 1:(5 - n.valid.facets))
        }
      }
    }

    #create output matrix and bind to it, RAFM format.
    else{
      outm <- matrix(NA, ncol(x), nrow(x)*2)

      #fill
      outm[,seq(1,ncol(outm),2)] <- xv1
      outm[,seq(2,ncol(outm),2)] <- xv2

      #replance missings with NA
      outm[outm == 0] <- NA

      #add column names
      colnames(outm) <- paste0("locus", sort(rep(1:nrow(x),2)), rep(c(".1", ".2"), ncol(outm)/2))

      #add subpop numbers, if given
      rdata <- cbind(x@sample.meta[,colnames(x@sample.meta) %in% facets],
                     as.data.frame(outm, stringsAsFactors = F))
    }
  }


  #presence/absence format
  if(output == "pa"){
    xv <- as.matrix(x)

    #function to produce vectorized presence absence (as much as possible, not vectorized for NAs)
    pa_alleles <- function(xv, snp_form, mDat, nsamp, nloci){
      nsamp <- nrow(xv)
      nloci <- ncol(xv)
      #get all possible genotypes
      gs <- unique(as.vector(xv))

      #which genotype is the missing data?
      mpos <- which(gs == mDat)

      #what are the possible alleles at all loci?
      as <- unique(c(substr(gs,1,snp_form/2), substr(gs, (snp_form/2 + 1), snp_form*2)))
      as <- sort(as[as != substr(mDat, 1, snp_form/2)]) #that aren't N?

      #make the table
      cat("Creating presence/absence table...\n")

      #convert genotype vector to allele vector
      xva1 <- substr(xv, 1, snp_form/2)
      xva2 <- substr(xv, (snp_form/2+1),snp_form)

      #initialize
      amat <- matrix(0, nsamp, length(as)*nloci)
      colnames(amat) <- paste0(sort(rep(1:nrow(x),length(as))), "_", as) #initialize all of the columns, with locus number followed by allele. Will remove anything empty after assignment.

      #fill in
      for(i in 1:length(as)){
        pr1 <- grep(as[i], xva1) # unique rows which have the allele in either position one or position two.
        pr2 <- grep(as[i], xva2)
        amat[,grep(paste0("_", as[i]),colnames(amat))][pr1] <- amat[,grep(paste0("_", as[i]),colnames(amat))][pr1] + 1 #set the allele as present in the correct rows. This works because we look only at the G amat columns first, then put set only the individual IDs with a G to one. I think.
        amat[,grep(paste0("_", as[i]),colnames(amat))][pr2] <- amat[,grep(paste0("_", as[i]),colnames(amat))][pr2] + 1
      }

      #remove alleles not seen at loci.
      amat <- amat[,which(colSums(amat) != 0)]

      #get allele counts per loci:
      cat("Filling in missing data with NAs.\n")

      ###########
      #fill in missing data with NAs.
      xmc <- which(xv == mDat) #which samples had missing data?
      adj <- floor(xmc / nsamp) #how many loci over do I need to adjust xmc, since in amat each locus occupies two columns?
      adj[xmc%%nsamp == 0] <- adj[xmc%%nsamp == 0] - 1 #shift over anything that got screwed up by being in the last sample
      xmc <- xmc + (nsamp*adj) #adjust xmc for extra columns.
      if(any(amat[xmc] != 0) | any(amat[xmc + nsamp] != 0)){
        stop("Missing data values were not properly identified for replacement with NAs. This usually happens when SNP data is not completely bi-allelic. Try filtering out non-biallelic and non-polymorphic SNPs using filter_snps.\n")
      }
      amat[xmc] <- NA #make the first allele NA
      amat[xmc + nsamp] <- NA #make the second allele (another column over) NA.

      return(amat)
    }

    amat <- pa_alleles(t(xv), 2, x@mDat)

    if(interp_miss){
      #average number observed in columns
      cat("Interpolating missing data...\n")
      afs <- colMeans(amat, TRUE)
      temp <- which(is.na(amat))/nrow(amat)
      fill_element <- floor(temp) + 1 #get the column for each missing data point
      fill_element[which(temp %% 1 == 0)] <- fill_element[which(temp %% 1 == 0)] - 1 #correct for anything in the last row
      amat[which(is.na(amat))] <- afs[fill_element] #fill with the appropriate allele frequency.
    }
    else{cat("Finished. Warning: Missing data counts are also stored!\n")}
    amat <- cbind(samp = as.character(colnames(x)), as.data.frame(amat, stringsAsFactors = F))
    rdata <- amat
  }


  #PLINK formatt
  if(output == "plink"){
    #============convert the genotypes into a binary string==========
    #======================make a vector containing the bits for each genotype============
    # much easier with sn!
    if(exists("pl_0g_dat")){
      # convert to vector
      xv <- as.vector(t(pl_0g_dat))

      # add a few filler genotype calls to fill in the empty bits later on.
      x.dup <- ncol(pl_0g_dat) %% 4
      if(x.dup  != 0){
        x.dup <- 4 - x.dup
      }
      x.dup <- matrix("FF", nrow(pl_0g_dat), x.dup)
      x.dup <- cbind(pl_0g_dat, x.dup)
      xv.dup <- as.vector(t(x.dup))

      # define the vector that we're going to fill with bits for each genotype
      xvt <- character(length(xv.dup))

      xvt[xv.dup == "0"] <- "00"
      xvt[xv.dup == "1"] <- "01"
      xvt[xv.dup == "2"] <- "11"
      xvt[xv.dup == miss] <- "10"
      xvt[xv.dup == "FF"] <- "00"

      # get allele names for map file down the line
      a.names <- matrix(c(0,1), nrow(header), 2, T)
    }

    # harder for the NN format
    else{
      # convert to vector
      xv <- as.vector(t(x))

      # add a few filler genotype calls to fill in the empty bits later on.
      x.dup <- ncol(x) %% 4
      if(x.dup  != 0){
        x.dup <- 4 - x.dup
      }
      x.dup <- matrix("FF", nrow(x), x.dup)
      x.dup <- cbind(x, x.dup)
      xv.dup <- as.vector(t(x.dup))

      # define the vector that we're going to fill with bits for each genotype
      xvt <- character(length(xv.dup))

      # pull allele table to check for violations
      as <- x@geno.tables$as[x@facet.meta$facet == ".base",]

      # if there are any alleles without any data, warn and remove
      if(any(rowSums(as) == 0)){
        missing.dat <- which(rowSums(as) == 0)
        warning(paste0("Removed ", length(missing.dat), "SNPs with no called genotypes.!\n"),
                "Removed SNPs:\n", missing.dat)
        x <- x[-missing.dat,]
        as <- as[-missing.dat,]
        xv <- as.vector(t(x))
      }

      # get the unique alleles for each loci, fixing for loci with one allele.
      unique.as <- cbind(as, N1 = rep(0, nrow(as))) # put in dummies for rows with missing haplotypes
      unique.as <- ifelse(unique.as != 0, TRUE, FALSE) # convert to logical
      if(any(rowSums(unique.as) != 2)){
        unique.as[rowSums(unique.as) != 2, 5] <- TRUE # add one missing allele note
      }
      unique.as <- which(t(unique.as)) %% ncol(unique.as) # get the column for each TRUE
      unique.as[unique.as == 0] <- 5 # fix for sixth column
      unique.as <- c(colnames(as), "M")[unique.as] # add column names.

      # match to observed genotypes
      ## these lines generate a vector equal in length to xv that has the correct genotype options per loci for all three genotypes.
      v00 <- paste0(unique.as[seq(1, length(unique.as), 2)],
                    unique.as[seq(1, length(unique.as), 2)])
      v00 <- rep(v00, each = ncol(x.dup))
      v11 <- paste0(unique.as[seq(2, length(unique.as), 2)],
                    unique.as[seq(2, length(unique.as), 2)])
      v11 <- rep(v11, each = ncol(x.dup))
      v01a <- paste0(unique.as[seq(1, length(unique.as), 2)],
                     unique.as[seq(2, length(unique.as), 2)])
      v01a <- rep(v01a, each = ncol(x.dup))
      v01b <- paste0(unique.as[seq(2, length(unique.as), 2)],
                     unique.as[seq(1, length(unique.as), 2)])
      v01b <- rep(v01b, each = ncol(x.dup))

      ## these lines create xvt, which has the correct binary doublets.
      xvt[xv.dup == v00] <- "00"
      xvt[xv.dup == v11] <- "11"
      xvt[xv.dup == v01a] <- "01"
      xvt[xv.dup == v01b] <- "01"
      xvt[xv.dup == "NN"] <- "10"
      xvt[xv.dup == "FF"] <- "00"

      # get allele names for map file down the line
      a.names <- matrix(unique.as, nrow(x), 2, T)
    }

    #=====================add magic numbers and reorder vector==============

    # order correctly for plink (weird reversed bytes)
    xvt <- paste0(xvt[seq(4, length(xvt), 4)],
                  xvt[seq(3, length(xvt), 4)],
                  xvt[seq(2, length(xvt), 4)],
                  xvt[seq(1, length(xvt), 4)])

    # add magic number and SNP-major identifier
    bed <- c("01101100", "00011011", "00000001", xvt)


    #===============make a ped file=================
    # make an empty set of ped header columns if not provided
    lower.sample.cols <- tolower(colnames(x@sample.meta))
    if(is.null(ped)){
      if(any(lower.sample.cols == "fam")){
        fam <- x@sample.meta[,which(lower.sample.cols == "fam")]
      }
      else{
        fam <- rep("NA", ncol(x))
      }
      if(any(lower.sample.cols == "patid")){
        PatID <- x@sample.meta[,which(lower.sample.cols == "patid")]
      }
      else{
        PatID <- rep("NA", ncol(x))
      }
      if(any(lower.sample.cols == "sex")){
        Sex <- x@sample.meta[,which(lower.sample.cols == "sex")]
      }
      else{
        Sex <- rep("NA", ncol(x))
      }
      if(any(lower.sample.cols == "phenotype")){
        Phenotype <- x@sample.meta[,which(lower.sample.cols == "phenotype")]
      }
      else{
        Phenotype <- rep("NA", ncol(x))
      }
      if(any(lower.sample.cols == "matid")){
        MatID <- x@sample.meta[,which(lower.sample.cols == "matid")]
      }
      else{
        MatID <- rep("NA", ncol(x))
      }
      ped <- data.frame(fam = fam,
                        ind = colnames(x),
                        PatID = PatID,
                        MatID = MatID,
                        Sex = Sex,
                        Phenotype = Phenotype)
    }

    # save .fam
    fam <- ped

    # change missing data value and add a space between alleles.
    x.np <- as.vector(t(x))
    x.np[x.np == x@mDat] <- "00"
    x.np <- gsub("(.)(.)", "\\1 \\2", x.np)

    # rebind
    ped <- cbind(ped, matrix(x.np, nrow(ped), nrow(x)), stringsAsFactors = F)

    #===============make an extended map file=================
    # with morgans
    if(any(colnames(x@snp.meta) %in% c("cM", "cm", "morgans"))){
      bim <- data.frame(chr = x@snp.meta$group,
                            rs = x@snp.meta$.snp.id,
                            cM = x@snp.meta[,which(colnames(x@snp.meta) %in% c("cM", "cm", "morgans"))],
                            bp = x@snp.meta$position,
                            a1 = a.names[,1],
                            a2 = a.names[,2])
    }
    # without morgans
    else{
      bim <- data.frame(chr = x@snp.meta$group,
                            rs = x@snp.meta$snp,
                            bp = x@snp.meta$position,
                            a1 = a.names[,1],
                            a2 = a.names[,2])
    }

    # recode chr
    bim$chr <- as.numeric(as.factor(bim$chr))

    # grab normal map file
    map <- bim[,-c(ncol(bim)-1, ncol(bim))]

    # name output
    rdata <- list(ped = ped, bed = bed, map = map, bim = bim)
  }


  # single-character numeric format
  if(output == "sn"){

    # grab major/minor info via calc_maf, unless already calculated
    if(!any(colnames(x@stats) == "major")){
      x <- calc_maf(x)
    }
    else if(any(is.na(x@stats$major[x@stats$facet == ".base"]))){
      x <- calc_maf(x)
    }

    min <- x@stats$minor[x@stats$facet == ".base"]
    maj <- x@stats$major[x@stats$facet == ".base"]

    # check to see if each allele is the minor, assign a one if so
    a1 <- substr(unlist(t(x)), 1, 1)
    a2 <- substr(unlist(t(x)), 2, 2)

    a1 <- a1 == rep(min, each = ncol(x))
    a2 <- a2 == rep(min, each = ncol(x))

    # collapse to output
    rdata <- t(a1 + a2)
    rdata[as.matrix(x) == x@mDat] <- -1
    rdata <- cbind(x@snp.meta[,-which(colnames(x@snp.meta) == ".snp.id")], as.data.frame(rdata))
  }

  #======================return the final product, printing an outfile if requested.=============
  if(outfile != FALSE){
    cat("Writing output file...\n")

    if(any(facets == ".base")){
      facets <- facets[-which(facets == ".base")]
    }

    # for genepop
    if(output == "genepop"){
      cat("\tPreparing genepop file...\n")
      # get list of snps
      llist <- paste0("SNP", "_", 1:ncol(rdata), ",")
      llist[length(llist)] <- paste0("SNP_", ncol(rdata))

      # write output
      cat(paste0(unlist(strsplit(outfile, split =  "(?<!^)\\.", perl = T))[1], "_genepop\n"), file = outfile)
      cat(llist, "\nPOP\n", file = outfile, append = T)

      # write the tables, splitting by pop if requested:
      if(length(facets) > 0){
        cat("\tWriting genepop file seperated by populations. First provided facet is treated as pop.\t")

        # sort by pop
        write.facets <- sort(unlist(strsplit(facets, split = "(?<!^)\\.", perl = T)))
        facet.cols <- match(write.facets, colnames(x@sample.meta))
        pop <- do.call(paste, as.data.frame(x@sample.meta[,facet.cols]))
        pop <- gsub(" ", ".", pop)
        rdata$pop <- pop

        pop <- sort(unique(rdata$pop)) # for later


        # sort and remove pop column
        rdata$rnames <- rownames(rdata)
        rdata <- dplyr::arrange(rdata, pop)
        pop.rows <- rdata$pop
        rownames(rdata) <- rdata$rnames
        rdata <- rdata[,-which(colnames(rdata) %in% c("pop", "rnames"))]

        #second loop prints results.
        for (i in 1:(length(pop))){
          cat(pop[i], "\t")
          data.table::fwrite(rdata[pop.rows == pop[i],], outfile, quote = F, sep = "\t", col.names = F, row.names = T, append = T)
          cat("POP\n", file = outfile, append = T)
        }
        cat("\t Done.\n")
      }
      else{
        data.table::fwrite(rdata, outfile, quote = F, sep = "\t", col.names = F, row.names = T, append = T)
      }
    }
    else if(output == "ac"){
      #write the raw output
      data.table::fwrite(rdata, outfile, quote = FALSE, col.names = T, sep = "\t", row.names = F)
      #write a bayescan object if pop list was provided.
      if(length(facets) > 0){
        outfile <- paste0(outfile, ".bayes")
        cat("\tWriting bayescan file seperated by populations. First provided facet is treated as pop.\t")
        pop <- x@sample.meta[,colnames(x@sample.meta) == facets[1]]
        u.pops <- unique(pop)

        trdat <- rdata[rdata$facet == facets[1], c("subfacet", ".snp.id", "n_total", "n_alleles", "ni1", "ni2")]

        #write the header
        cat("[loci]=", nrow(x), "\n\n[populations]=", length(u.pops), "\n\n", file = outfile, sep = "")

        #write the data for each population.
        for(i in 1:length(u.pops)){
          cat("[pop]=", i, "\n", file = outfile, append = T, sep = "") #write header

          tdat <- trdat[trdat$subfacet == u.pops[i],]

          wdat <- cbind(snp = tdat$.snp.id, tdat[,3:6])

          data.table::fwrite(wdat,
                             outfile, col.names = F, row.names = F, quote = F, sep = "\t",
                             append = T) # write the data for this population.

          cat("\n", file = outfile, append = T) # add a line break
        }
      }
    }
    else if(output == "structure" | output == "faststructure"){
      data.table::fwrite(rdata, outfile, quote = FALSE, col.names = F, sep = "\t", row.names = F)
    }
    else if(output == "plink"){
      t2 <- as.logical(as.numeric(unlist(strsplit(bed, "")))) # merge the data and convert it to a logical.
      t2 <- BMS::bin2hex(t2) # convet to hex codes
      t2 <- unlist(strsplit(t2,"")) # unlist
      t2 <- paste0("\\\\x",
                   t2[seq(1,length(t2),2)],
                   t2[seq(2,length(t2),2)],
                   collapse = "")
      writeLines(paste0("#!/bin/bash\n\n",
        "echo -n -e ", t2, " > ", outfile, ".bed"), paste0(outfile, ".sh"))
      warning(paste0("To get PLINK .bed file, run ", outfile, ".sh.\n"))
      data.table::fwrite(map, paste0(outfile, ".map"), quote = F, col.names = F, sep = "\t", row.names = F)
      data.table::fwrite(ped, paste0(outfile, ".ped"), quote = F, col.names = F, sep = "\t", row.names = F)
      data.table::fwrite(bim, paste0(outfile, ".bim"), quote = F, col.names = F, sep = "\t", row.names = F)
    }
    else{
      data.table::fwrite(rdata, outfile, quote = FALSE, col.names = T, sep = "\t", row.names = F)
    }
  }

  #return results
  else{
    return(rdata)
  }
}

