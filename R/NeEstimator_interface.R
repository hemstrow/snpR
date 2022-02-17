write_neestimator_inputs <- function(x, facets, chr = NULL, methods = "LD",
                                     temporal_methods = c("Pollak", "Nei", "Jorde"),
                                     temporal_gens = NULL,
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
  facets <- .check.snpR.facet.request(x, facets)

  msg <- character()
  # check methods
  good.methods <- c("ld", "het", "coan", "temporal")
  methods <- tolower(methods)
  if(any(!methods %in% good.methods)){
    msg <- c(msg,paste0("Unaccepted methods. Acceptable methods: ", paste(good.methods, collapse = ", ")))
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
    tl <- .get.task.list(x, facets)
    if(nrow(temporal_gens) != nrow(tl)){
      msg <- c(msg, paste0("Not enough generation info provided. The number of rows of temporal_gens must equal the number
                           of unique options for the facets provided (", nrow(tl), ")."))
    }
    else{
      cat("Please ensure that the temporal_gens data.frame is in alphabetical order relative to the levels
          of the facet provided (", paste(tl[,2], collapse = ", "), ").\n")
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
  format_snps(x, "genepop", facets, outfile = "snps.gen")

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
    pcrit <- numeric()
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
    utils::write.table(temporal_gens, "info", col.names = F, row.names = F, sep = " ", append = T, quote = F)
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

run_neestimator <- function(NeEstimator_path = "/usr/bin/Ne2-1.exe", data_path = "NeEstimator/"){
  #==========run NeEstimator==========
  owd <- getwd()
  setwd(data_path)
  call <- paste(NeEstimator_path, "i:info", "o:option", collapse = " ")
  if(Sys.info()["sysname"] == "Windows"){
    shell(call)
  }
  else{
    system(call)
  }
  setwd(owd)
}

parse_neestimator <- function(path = "NeEstimator/", pattern = "ne_out", facets = NULL, snpRdat = NULL){
  owd <- getwd()
  setwd(path)
  #==========parse output=============
  files <- list.files(pattern = pattern)
  files <- files[-which(files == paste0(pattern, ".txt"))]
  out <- vector("list", length(files))

  for(i in 1:length(files)){

    # figure out which type of file we are using
    if(length(grep("Cn", files[i]) != 0)){type <- "cn"; skip <- 13}
    else if(length(grep("LD", files[i]) != 0)){type <- "LD"; skip <- 14}
    else if(length(grep("Ht", files[i]) != 0)){type <- "Ht"; skip <- 14}

    # read in the data
    suppressWarnings(dat <- as.data.frame(readr::read_table2(files[i], skip = skip, col_names = F)))
    dat <- dat[1:(nrow(dat) - 2),]

    # add column names, adjust data per type
    if(type == "cn"){
      colnames(dat) <- c("pop", "n",	"hm_n", "f^1", "Neb", "Neb_lCI", "Neb_uCI")
      dat <- dat[,c(1, 5:7)]
    }
    else if(type == "LD" | type == "Ht"){
      fix.rows <- which(rowSums(is.na(dat)) == 2)
      dat[fix.rows, 3:ncol(dat)] <- dat[fix.rows, 1:(ncol(dat) - 2)]
      dat[fix.rows, 1:2] <- dat[-fix.rows, 1:2][sort(rep(1:(nrow(dat) - length(fix.rows)), length.out = length(fix.rows))),]
      if(type == "LD"){
        colnames(dat) <- c("pop", "n", "pcrit",	"hm_n",	"n_alleles", "r^2", "Exp(r^2)", "LDNe",
                           "LDNe_lCIp", "LDNe_uCIp",
                           "LDNe_lCIj", "LDNe_uCIj", "Eff.df")
        dat <- dat[,c(1,3,8:12)]
      }
      else{
        colnames(dat) <- c("pop", "n", "pcrit",	"hm_n",	"n_alleles", "D", "He_Ne", "He_lCIp", "He_uCIp")
        dat <- dat[,c(1,3,7:9)]
      }
    }
    out[[i]] <- dat
    names(out)[i] <- type
  }
  setwd("..")


  # fix pop names if possible
  if(!is.null(snpRdat)){
    facets <- .check.snpR.facet.request(snpRdat, facets)
    snpRdat <- .add.facets.snpR.data(snpRdat, facets)
    opts <- .get.task.list(snpRdat, facets)
    opts <- opts[,2] # these are the sorted options.
    tab <- data.frame(index = unique(out[[1]]$pop), ref = opts)
    out <- lapply(out, function(x){x$pop <- tab$ref[match(x$pop, tab$index)];return(x)})
  }


  # clean and merge
  out <- purrr::reduce(out, dplyr::left_join)
  out <- as.data.frame(out)
  mc <- which(colnames(out) == "pcrit")
  out <- out[,c(1, mc, (2:ncol(out))[-(mc - 1)])]
  out[out == "Infinite"] <- Inf
  out[,-1] <- dplyr::mutate_all(out[,-1], as.numeric)
  out[,grep("Ne", colnames(out))][out[,grep("Ne", colnames(out))] < 0] <- Inf
  
  # cast such that pcrit for each pop is applied across columns
  cout <- as.data.frame(tidyr::pivot_wider(out, id_cols = "pop", names_from = "pcrit", values_from = colnames(out)[-c(1:2)]))
  
  
  return(cout)

}
