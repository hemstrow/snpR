write_neestimator_inputs <- function(x, facets, chr = NULL, methods = "LD",
                                     temporal_methods = c("Pollak", "Nei", "Jorde"),
                                     temporal_details = NULL,
                                     pcrit = c(0.05, 0.02, 0.01), mating = "random",
                                     outfile = "ne_out.txt", max_ind_per_pop = NULL){
  #=====================initialize===============
  # initialize storage directory
  if(!dir.exists("NeEstimator")){
    dir.create("NeEstimator")
    setwd("NeEstimator")
  }
  else{
    setwd("NeEstimator")
  }

  if(file.exists("snps.gen")){
    file.remove("snps.gen")
  }

  #=====================sanity checks============
  ffacets <- .check.snpR.facet.request(x, facets)
  if("temporal" %in% tolower(methods)){
    if(length(facets) > 1){
      stop("Only one facet can be run at a time for the temporal method. Complex facets ('pop.fam') are permitted.\n")
    }
    if(facets == ".base"){
      stop("The temporal method cannot be used for the base level facet.\n")
    }
    
    if(facets[1] != ffacets[1]){
      stop("For the temporal method, complex facets must be provided in alphabetical order in both the 'facets' and 'temporal_details' arguments (for example, 'fam.pop' not 'pop.fam') to prevent downstream issues.")
    }
  }
  facets <- ffacets
  rm(ffacets)

  msg <- character()
  # check methods
  good.methods <- c("ld", "het", "coan", "temporal") # temporal temporarily
  methods <- tolower(methods)
  if(any(!methods %in% good.methods)){
    msg <- c(msg,paste0("Unaccepted methods. Acceptable methods: ", paste(good.methods, collapse = ", ")))
  }
  
  if("temporal" %in% methods & length(methods) != 1){
    msg <- c(msg, "The temporal method cannot currently be run alongside other methods.\n")
  }

  # check temporal methods
  if("temporal" %in% methods){
    
    # check methods
    temporal_methods <- tolower(temporal_methods)
    good.methods <- c("pollak", "nei", "jorde")
    if(any(!temporal_methods %in% good.methods)){
      msg <- c(msg, paste0("Unaccepted temporal methods. Acceptable temporal methods: ", paste(good.methods, collapse = ", ")))
    }

    # check and adjust gens
    x <- .add.facets.snpR.data(x, facets)
    options <- temporal_details[,1:2]
    if(any(options[,1] == options[,2])){
      msg <- c(msg, "Some facet levels are both the time = 1 and time = 2 for the same comparison!\n")
    }
    options <- c(options[,1], options[,2])
    possible_levels <- unlist(lapply(summarize_facets(x, facets), names))
    
    bad_levels <- !options %in% possible_levels
    if(any(bad_levels)){
      msg <- c(msg, paste0("Some facet levels from 'temporal_details' are not recognized. Bad levels:", 
                           paste0(options[bad_levels], collapse = ", "), "\n"))
    }
    
    if(any(temporal_details[,3] <= 0)){
      msg <- c(msg, "All time (generation) gaps specified in 'temporal_details' must be greater than zero.\n")
    }
  }

  # check mating type
  mating <- tolower(mating)
  good.methods <- c("random", "monogamy")
  if(any(!mating %in% good.methods)){
    msg <- c(msg, paste0("Unaccepted mating type. Acceptable mating types: ", paste(good.methods, collapse = ", ")))
  }
  if(length(mating) != 1){
    msg <- c(msg, "Exactly one mating type must be specified.")
  }

  if(length(msg) > 0){
    stop(paste(msg, collapse = "\n"))
  }

  #====================write the data==============
  # format x into genepop format
  if(methods[1] != "temporal"){
    format_snps(x, "genepop", facets, outfile = "snps.gen")
  }
  else{
    uopts <- unique(options)
    fd <- vector("list", length(uopts))
    names(fd) <- uopts
    
    # format
    for(i in 1:length(uopts)){
      .suppress_specific_warning(fd[[i]] <- format_snps(.subset_snpR_data(x, facets = facets, subfacets = uopts[i]), 
                                            output = "genepop"),
                                 "samples were subset")
      fd[[i]] <- cbind(samp = row.names(fd[[i]]), fd[[i]])
    }
    
    # write, two loops because levels might be duplicated and we can save work
    write("snps.gen", "snps.gen")
    write(paste0(paste0("SNP_", 1:nrow(x)), collapse = ", "), "snps.gen", append = TRUE)
    options <- as.character(t(temporal_details[,1:2]))
    for(i in 1:length(options)){
      write("POP", "snps.gen", append = TRUE)
      write.table(fd[[options[i]]], "snps.gen", append = TRUE, 
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
      
    }
    
  }

  # write a chr input file
  if(!is.null(chr[1])){
    chrs <- cbind(paste0("c", as.numeric(as.factor(x@snp.meta[,chr]))), paste0("SNP_", 1:nrow(x)))
    utils::write.table(chrs, "chrinput.map", quote = F, sep = " ", col.names = F, row.names = F)
  }


  #====================write the info file=========
  # methods:
  method.tab <- data.frame(method = c("ld", "het", "coan", "temporal"),
                           val = c(1, 2, 4, 8))
  method.sum <- sum(method.tab[match(methods, method.tab$method), 2])
  ## if temporal...
  if("temporal" %in% methods){
    temporal.tab <- data.frame(method = c("pollak", "nei", "jorde"),
                               val = c(1, 2, 4))
    k <- sum(method.tab[match(temporal_methods, temporal.tab$method), 2])
    if(k == 7){k <- 0} # If doing all, set to 0
    method.sum <- paste(method.sum, k, "   ")
  }
  ## if only doing molecular coancestry, no pcrits
  if(length(methods) == 1 & methods[1] == "coan"){
    pcrit <- 0
  }

  # organize everthing that needs to be written:
  to.write <- c(method.sum,
                paste0(getwd(), "/"),
                "snps.gen",
                2,
                paste0(getwd(), "/"),
                outfile,
                length(pcrit),
                paste(pcrit, collapse = "  "),
                ifelse(mating == "random", 0, 1))

  # write the info file
  writeLines(to.write, "info")

  # if we are doing temporal, need to add stuff
  if("temporal" %in% methods){
    if(ncol(temporal_details) == 4){ # if plan I
      utils::write.table(data.frame(temporal_details[,4], 0, temporal_details[,3]), "info", col.names = F, row.names = F, sep = " ", append = T, quote = F)
    }
    else{ # if plan II
      utils::write.table(data.frame(0, 0, temporal_details[,3]), "info", col.names = F, row.names = F, sep = " ", append = T, quote = F)
    }
    write("0", "info", append = T)
  }

  #====================write the option file=======
  to.write <- c(paste0(method.sum, "  ", 0, "  ", length(pcrit), "  ", 1),
                ifelse(is.null(max_ind_per_pop), 0, max_ind_per_pop),
                0,
                0,
                1,
                1,
                0,
                0,
                0,
                ifelse(is.null(chr), 0, paste0("2  chrinput.map")))
  writeLines(to.write, "option")
  setwd("..")
}

run_neestimator <- function(NeEstimator_path = "/usr/bin/Ne2-1.exe", data_path = "NeEstimator/", verbose = TRUE){
  #==========run NeEstimator==========
  owd <- getwd()
  setwd(data_path)
  call <- paste(NeEstimator_path, "i:info", "o:option", collapse = " ")
  
  if(Sys.info()["sysname"] == "Windows"){
    system(call, show.output.on.console = verbose)
  }
  else{
    system(call, ignore.stdout = !verbose)
  }
  
  setwd(owd)
}

parse_neestimator <- function(path = "NeEstimator/", pattern = "ne_out", facets = NULL, snpRdat = NULL, temporal_details = NULL, temporal_methods = NULL){
  owd <- getwd()
  setwd(path)
  #==========parse output=============
  files <- list.files(pattern = pattern)
  files <- files[-which(files == paste0(pattern, ".txt"))]
  out <- vector("list", length(files))
  if(length(files) == 0){
    stop("No NeEstimator output files found. Double check your NeEstimator path and look for any errors produced by NeEstimator. If the problem persists, please leave an issue at https://github.com/hemstrow/snpR/issues .")
  }

  for(i in 1:length(files)){

    # figure out which type of file we are using
    if(length(grep("Cn", files[i]) != 0)){type <- "cn"; skip <- 13}
    else if(length(grep("LD", files[i]) != 0)){type <- "LD"; skip <- 14}
    else if(length(grep("Ht", files[i]) != 0)){type <- "Ht"; skip <- 14}
    else if(length(grep("Tp", files[i]) != 0)){type <- "Tp"; skip <- 14}
    
    # read in the data
    suppressWarnings(dat <- as.data.frame(data.table::fread(files[i], skip = skip, header = FALSE, sep = "\t")))
    if(ncol(dat) == 1){
      dat <- unlist(strsplit(dat[1,], "\t"))
      dat <- gsub(" +", "", dat)
      dat <- as.data.frame(matrix(dat, nrow = 1))
    }
    # dat <- dat[1:(nrow(dat) - 2),]

    # add column names, adjust data per type
    if(type == "cn"){
      colnames(dat) <- c("pop", "n",	"hm_n", "f^1", "Neb", "Neb_lCI", "Neb_uCI")
      dat$pcrit <- 0
      dat <- dat[,c(1, 8, 5:7)]
    }
    else if(type == "LD" | type == "Ht"){
      fix.rows <- which(rowSums(is.na(dat)) == 2)
      if(length(fix.rows) > 0){
        dat[fix.rows, 3:ncol(dat)] <- dat[fix.rows, 1:(ncol(dat) - 2)]
        dat[fix.rows, 1:2] <- dat[-fix.rows, 1:2][sort(rep(1:(nrow(dat) - length(fix.rows)), length.out = length(fix.rows))),]
      }
      if(type == "LD"){
        if(ncol(dat) == 13){
          colnames(dat) <- c("pop", "n", "pcrit",	"hm_n",	"n_alleles", "r^2", "Exp(r^2)", "LDNe",
                             "LDNe_lCIp", "LDNe_uCIp",
                             "LDNe_lCIj", "LDNe_uCIj", "Eff.df")
          dat <- dat[,c(1,3,8:12)]
        }
        else{
          # pcrit doesn't get printed if pcrit is only one option for some reason.
          colnames(dat) <- c("pop", "n", "hm_n",	"n_alleles", "r^2", "Exp(r^2)", "LDNe",
                             "LDNe_lCIp", "LDNe_uCIp",
                             "LDNe_lCIj", "LDNe_uCIj", "Eff.df")
          dat <- dat[,c(1,7:11)]
          pcrit <- readLines(files[i], n = 10)
          pcrit <- pcrit[10]
          pcrit <- gsub("Lowest allele frequency used: +", "", pcrit)
          pcrit <- as.numeric(pcrit)
          dat$pcrit <- pcrit
          dat <- dat[,c(1, ncol(dat), 2:(ncol(dat) - 1))]
        }
        
      }
      else{
        if(ncol(dat) == 9){
          colnames(dat) <- c("pop", "n", "pcrit",	"hm_n",	"n_alleles", "D", "He_Ne", "He_lCIp", "He_uCIp")
          dat <- dat[,c(1,3,7:9)]
        }
        else{
          # pcrit doesn't get printed if pcrit is only one option for some reason.
          colnames(dat) <- c("pop", "n", "hm_n",	"n_alleles", "D", "He_Ne", "He_lCIp", "He_uCIp")
          dat <- dat[,c(1,6:8)]
          pcrit <- readLines(files[i], n = 10)
          pcrit <- pcrit[10]
          pcrit <- gsub("Lowest allele frequency used: +", "", pcrit)
          pcrit <- as.numeric(pcrit)
          dat$pcrit <- pcrit
          dat <- dat[,c(1, ncol(dat), 2:(ncol(dat) - 1))]
        }
        
      }
      
      dat$pop <- rep(dat$pop[seq(1, nrow(dat), length(unique(dat$pcrit)))], each = length(unique(dat$pcrit)))
    }
    else if(type == "Tp"){
      
      # unfortunately in this case the nice condensed tabular output doesn't print anything but pcrit = 0, so need to parse the larger file
      res <- readLines("ne_out.txt")
      pcrits <- res[grep("Lowest Allele Frequency Used", res)[1]]
      pcrits <- strsplit(pcrits, " +")[[1]]
      pcrits <- pcrits[-c(1:4)]
      pcrits <- as.numeric(gsub("\\+", "", pcrits))
      
      pop_opts <- paste0(temporal_details[,1], "~", temporal_details[,2])
      dat <- expand.grid(pop = pop_opts,
                        pcrit = pcrits)
      dat <- as.data.table(dat)
      dat$pop <- as.character(dat$pop)
      
      
      # function to parse one method part
      parse_tmp_method_part <- function(val, name, pops, pcrits){
        val <- data.table::as.data.table(t(as.data.frame(val)))
        colnames(val) <- as.character(pcrits)
        val$pop <- pops
        val <- data.table::melt(val, id.vars = "pop")
        val$variable <- as.numeric(as.character(val$variable))
        colnames(val) <- c("pop", "pcrit", name)
        return(val)
      }
      
      # function to parse one method
      parse_tmp_method <- function(res, matches, name, pcrits, pops){
        matches <- matches[seq(2, length(matches), 2)]
        
        # ne
        ne <- strsplit(res[(matches  + 4)], " +")
        ne <- lapply(ne, function(y) y[-c(1:4)])
        ne <- parse_tmp_method_part(ne, "_Ne", pops, pcrits)
        
        # lCIp
        lCIp <- strsplit(res[(matches  + 7)], " +")
        lCIp <- lapply(lCIp, function(y) y[-c(1:3)])
        lCIp <- parse_tmp_method_part(lCIp, "_lCIp", pops, pcrits)
        
        # uCIp
        uCIp <- strsplit(res[(matches  + 8)], " +")
        uCIp <- lapply(uCIp, function(y) y[-1])
        uCIp <- parse_tmp_method_part(uCIp, "_uCIp", pops, pcrits)
        
        # lCIj
        lCIj <- strsplit(res[(matches  + 9)], " +")
        lCIj <- lapply(lCIj, function(y) y[-c(1:5)])
        lCIj <- parse_tmp_method_part(lCIj, "_lCIj", pops, pcrits)
        
        # uCIj
        uCIj <- strsplit(res[(matches  + 10)], " +")
        uCIj <- lapply(uCIj, function(y) y[-1])
        uCIj <- parse_tmp_method_part(uCIj, "_uCIj", pops, pcrits)
        
        
        matches_res <- cbind(ne, lCIp[,3], uCIp[,3], lCIj[,3], uCIj[,3])
        colnames(matches_res)[-c(1:2)] <- paste0(name, colnames(matches_res)[-c(1:2)])
        return(matches_res)
      }
      
      # run for each method
      pol <- grep("Pollak", res)
      if(length(pol) > 0){
        pol <- parse_tmp_method(res, pol, "Pollak", pcrits, pop_opts)
        dat <- merge(dat, pol)
      }
      nei <- grep("Nei", res)
      if(length(nei) > 0){
        nei <- parse_tmp_method(res, nei, "Nei", pcrits, pop_opts)
        dat <- merge(dat, nei)
      }
      jor <- grep("Jorde", res)
      if(length(jor) > 0){
        jor <- parse_tmp_method(res, jor, "Jorde", pcrits, pop_opts)
        dat <- merge(dat, jor)
      }
    }
    out[[i]] <- dat
    names(out)[i] <- type
  }
  setwd(owd)
  
  # fix pop names if possible
  if(!is.null(snpRdat)){
    facets <- .check.snpR.facet.request(snpRdat, facets)
    snpRdat <- .add.facets.snpR.data(snpRdat, facets)
    if(type != "Tp"){
      opts <- .get.task.list(snpRdat, facets)
      opts <- opts[,2] # these are the sorted options.
      index <- unique(out[[1]]$pop)
      if(any(index == "")){index <- index[-which(index == "")]}
      tab <- data.frame(index = index, ref = opts) #err
      out <- lapply(out, function(x){x$pop <- tab$ref[match(x$pop, tab$index)];return(x)})
    }
  }
  
  
  # clean and merge
  out <- purrr::reduce(out, dplyr::full_join)
  out <- as.data.frame(out)
  mc <- which(colnames(out) == "pcrit")
  out <- out[,c(1, mc, (2:ncol(out))[-(mc - 1)])]
  out[out == "Infinite"] <- Inf
  out[,-1] <- dplyr::mutate_all(out[,-1], as.numeric)
  out[,grep("Ne", colnames(out))][out[,grep("Ne", colnames(out))] < 0] <- Inf # negative values mean inf
  
  # cast such that pcrit for each pop is applied across columns
  cout <- vector("list", ncol(out) - 2)
  for(i in 1:length(cout)){
    cout[[i]] <- reshape2::dcast(out, pop~pcrit, value.var = colnames(out)[i + 2])
    
    if(i != 1){
      cout[[i]] <- cout[[i]][,-1, drop = FALSE]
      colnames(cout[[i]]) <- paste0(colnames(out)[i + 2], "_", colnames(cout[[i]]))
    }
    else{
      colnames(cout[[i]])[-1] <- paste0(colnames(out)[i + 2], "_", colnames(cout[[i]])[-1])
    }
    
    
  }
  cout <- dplyr::bind_cols(cout)
  
  return(cout)
  
}
