# internals, re-distributed from pophelper (OK given GPL-3 liscence)
.clumppExport <- function(qlist=NULL,prefix=NA,parammode=NA,paramrep=NA,exportpath=NULL,path=NULL) {
  # check input
  .is.qlist(qlist)
  
  if(!is.null(path)) {
    warning("clumppExport: Argument 'path' is deprecated. Use 'exportpath'.")
    exportpath <- path
  }
  if(is.null(exportpath)) stop("clumppExport: Argument 'exportpath' is missing. Specify an output directory. To use the current working directory, set exportpath=getwd().")
  
  if(is.na(prefix)) prefix <- "pop"
  prefix <- paste0(prefix,"_K")
  
  # get tabulated runs
  df1 <- .tabulateQ(qlist)
  df2 <- .summariseQ(df1)
  df1l <- as.list(df1)
  df2l <- as.list(df2)
  
  if(is.null(names(qlist))) names(qlist) <- paste0("sample",1:length(qlist))
  
  # k val duplicated
  if(any(duplicated(df2l$k))) stop("clumppExport: Repeating values of K found.")
  # do ind vary?
  if(!all(df2l$ind[1]==df2l$ind)) warning("clumppExport: Number of individuals vary between runs.")
  
  e <- 1
  p <- 1
  len1 <- length(df2l$k)
  while (e <= len1)
  {
    k <- df2l$k[e]
    ind <- df2l$ind[e]
    runs <- df2l$runs[e]
    
    ldata <- vector("list",length=runs)
    for (f in 1:runs)
    {
      sel <- which(names(qlist)==as.character(df1l$file[p]))
      dframe1 <- qlist[[sel]]
      
      # generate df
      dframe3 <- as.matrix(data.frame(V1=paste0(1:ind,":"),dframe1,last=as.character(rep(1,ind)),stringsAsFactors=FALSE))
      
      # add dataframes to list
      ldata[[f]] <- dframe3
      rm(dframe3)
      p=p+1
    }
    
    if(runs > 1 && k > 1)
    {
      if(as.numeric(file.access(exportpath,2))==-1) stop(paste0("clumppExport: Directory ",exportpath," has no write permission."))
      
      currpath <- paste0(exportpath,"/",prefix,k)
      dir.create(currpath)
      #setwd(paste0(path,"/",prefix,k))
      message(paste0("Folder created: ",basename(currpath)))
      out <- paste0(prefix,k,"-combined.txt")
      
      ## file output block
      
      # make 2 line space
      spacer <- matrix(rep("  ",(k+2)*2),nrow=2)
      
      # write file
      write(x=t(format(ldata[[1]],nsmall=15)),file=file.path(currpath,out),ncolumns=k+2)
      for (i in 2:length(ldata))
      {
        write(x=t(spacer),file=file.path(currpath,out),ncolumns=k+2,append=TRUE)
        write(x=t(format(ldata[[i]],nsmall=15)),file=file.path(currpath,out),ncolumns=k+2,append=TRUE)
      }
      message(paste0(out," exported."))
      
      ## paramfile section
      T1 <- factorial(k)*((length(ldata)*(length(ldata)-1))/2)*k*ind
      if(T1 <= 100000000)
      {
        if(is.na(parammode)) parammode <- 2
        if(is.na(paramrep)) paramrep <- 20
      }else{
        if(is.na(parammode)) parammode <- 3
        if(is.na(paramrep)) paramrep <- 500
      }
      out1 <- gsub(".txt","",out)
      params <- c("DATATYPE 1 ",
                  "INDFILE NOTNEEDED.indfile ",
                  paste0("POPFILE ",out," "),
                  paste0("OUTFILE ",out1,"-merged.txt "),
                  paste0("MISCFILE ",out1,"-miscfile.txt "),
                  paste0("K ",k," "),
                  paste0("C ",ind," "),
                  paste0("R ",length(ldata)," "),
                  paste0("M ",parammode," "),
                  "W 0 ",
                  "S 2 ",
                  "GREEDY_OPTION 2 ",
                  paste0("REPEATS ",paramrep," "),
                  "PERMUTATIONFILE NOTNEEDED.permutationfile ",
                  "PRINT_PERMUTED_DATA 1 ",
                  paste0("PERMUTED_DATAFILE ",out1,"-aligned.txt "),
                  "PRINT_EVERY_PERM 0 ",
                  paste0("EVERY_PERMFILE ",out1,".every_permfile "),
                  "PRINT_RANDOM_INPUTORDER 0 ",
                  paste0("RANDOM_INPUTORDERFILE ",out1,".random_inputorderfile "),
                  "OVERRIDE_WARNINGS 0 ",
                  "ORDER_BY_RUN 0 ")
      
      write(x=params,file=file.path(currpath,"paramfile"))
      message("paramfile exported.")
      
      message("-----------------------")
    }else
    {
      if(k==1) message(paste0(prefix,k," not exported. K less than 2."))
      if(runs < 2) message(paste0(prefix,k," not exported. Repeats less than 2."))
      message("-----------------------")
    }
    e <- e + 1
  }
  
  message("Run completed.")
}

.collectClumppOutput <- function(prefix="pop",filetype="aligned",runsdir=NULL,newdir=NULL) {
  
  if(!is.character(prefix)) stop("collectClumppOutput: Argument 'prefix' must be a character.")
  if(!is.character(filetype)) stop("collectClumppOutput: Argument 'filetype' must be a character.")
  
  # check imgoutput
  if(tolower(filetype)!="aligned" && tolower(filetype)!="merged" && tolower(filetype)!="both") stop("collectClumppOutput: Argument 'filetype' set incorrectly. Set as 'aligned', 'merged' or 'both'.")
  
  if(is.null(runsdir)) stop("collectClumppOutput: Argument 'runsdir' missing. Specify a path/directory with CLUMPP results.")
  if(is.null(newdir)) stop("collectClumppOutput: Argument 'newdir' missing. Specify an output path/directory.")
  dirs <- list.dirs(path=runsdir,full.names=TRUE,recursive=FALSE)
  dirs <- dirs[grep(prefix,basename(dirs))]
  if(length(dirs)==0) stop("collectClumppOutput: No directories found with the specified prefix.")
  
  if(!dir.exists(newdir)) {
    dir.create(newdir)
    message(paste0("Directory ",newdir," created."))
  }
  
  if(as.numeric(file.access(newdir,2))==-1) stop(paste0("collectClumppOutput: Directory (",newdir,") has no write permission."))
  
  k <- 0
  l <- 0
  i <- 1
  for (i in seq_along(dirs)) {
    files <- list.files(path=dirs[i])
    sel1 <- grep("aligned",files)
    sel2 <- grep("merged",files)
    if(tolower(filetype)=="aligned") sel3 <- sel1
    if(tolower(filetype)=="merged") sel3 <- sel2
    if(tolower(filetype)=="both") sel3 <- c(sel1,sel2)
    if(length(sel3)==0) {
      warning(paste0("No suitable file found in directory: ",basename(dirs[i])))
    }
    if(length(sel3) != 0) {
      file.copy(from=file.path(dirs[i],files[sel3]),to=newdir)
      k=k+1
      l=l+length(sel3)
    }
  }
  
  message(paste0("Directories processed: ",k,"\nFiles copied: ",l,"\n"))
  
  return(c(k,l))
}

.checkQ <- function(files=NULL,warn=FALSE) {
  
  if(is.null(files)) stop("checkQ: Input is empty.")
  if(class(files) != "list" && class(files) != "character") stop("checkQ: Input is not a character or list datatype.")
  
  len1 <- length(files)
  
  checkvec <- rep("UNIDENTIFIED",length=len1)
  subtype <- rep(NA,length=len1)
  for(i in seq_along(files))
  {
    chk <- FALSE
    
    if(class(files)=="list")
    {
      if(class(files[[i]])=="data.frame") 
      {
        chk <- TRUE
        checkvec[i] <- "data.frame"
      }
    }
    
    if(!chk)
    {
      
      read1 <- readLines(files[i],n=7,warn=FALSE)
      
      # read BAPS file
      if(!chk)
      {
        chk <- any(grepl("RESULTS OF ADMIXTURE ANALYSIS BASED",toupper(read1)))
        if(chk) checkvec[i] <- "BAPS"
      }
      
      # read TESS file
      if(!chk)
      {
        chk <- grepl("ESTIMATED CLUSTERING PROBABILITIES",toupper(read1)[1])
        if(chk) checkvec[i] <- "TESS"
      }
      
      # read STRUCTURE file
      if(!chk)
      {
        chk <- grepl("STRUCTURE BY PRITCHARD",toupper(read1)[4])
        if(chk) checkvec[i] <- "STRUCTURE"
      }
      
      # read BASIC files
      
      rm(read1)
      
      if(!chk)
      {
        seps <- c("","\t",",")
        subtypes <- c("SPACE","TAB","COMMA")
        k=1
        while(!chk)
        {
          if(class(try(suppressWarnings(read.table(files[i],header=FALSE,sep=seps[k],nrows=1,quote="",stringsAsFactors=FALSE))))!="try-error")
          {
            df <- read.table(files[i],header=FALSE,sep=seps[k],nrows=1,quote="",stringsAsFactors=FALSE)
            if(all(sapply(df,is.numeric))) {
              checkvec[i] <- "BASIC"
              subtype[i] <- subtypes[k]
              chk <- TRUE
            }else{
              if((ncol(df) > 2) && (is.character(df[,1])))
              {
                checkvec[i] <- "CLUMPP"
                chk <- TRUE
              }
            }
          }
          k=k+1
          if(k>3) break
        }
      }
    }
    if((!chk) && warn) warning(paste0("checkQ: ",files[i]," is not a STRUCTURE, TESS, BAPS, BASIC or CLUMPP file.\n"))
  }
  
  return(list(type=checkvec,subtype=subtype))
}

.is.qlist <- function(qlist=NULL) {
  
  if(is.null(qlist)) stop("is.qlist: Input is empty.")
  
  # is it a list?
  if(!is.list(qlist)) stop("is.qlist: Input is not a list object.")
  
  # does it have dataframes?
  if(!all(sapply(qlist,is.data.frame))) stop("is.qlist: One or more list elements are not data.frame datatype.")
  
  # any NA?
  if(any(is.na(qlist))) stop("is.qlist: One or more list elements are NA.")
  
  # any NAs in dataframe?
  if(any(unlist(lapply(qlist,sapply,is.na)))) warning("is.qlist: One or more qlist dataframes have NA in data.")
  
  # are dataframes all numeric?
  if(!all(unlist(lapply(qlist,sapply,is.numeric)))) stop("is.qlist: One or more qlist dataframes have non-numeric columns.")
  
  # are there duplicated rownames?
  if(any(unlist(lapply(lapply(qlist,rownames),duplicated)))) stop("is.qlist: One or more qlist dataframes have duplicated row names.")
  
  # are there duplicated column names?
  if(any(unlist(lapply(lapply(qlist,colnames),duplicated)))) stop("is.qlist: One or more qlist dataframes have duplicated column names.")
  
  # are dataframe lists named?
  if(is.null(names(qlist))) stop("is.qlist: List elements are not named. Add names or run as.qlist().")
  
  # are dataframe lists named?
  if(any(is.na(names(qlist)))) stop("is.qlist: One or more list element name is NA. Add names or run as.qlist().")
  
  # are all dataframe lists named?
  if(any(nchar(names(qlist))==0)) stop("is.qlist: One or more list elements are not named. Add names or run as.qlist().")
  
  # are there duplicated list names?
  if(any(duplicated(names(qlist)))) stop("is.qlist: One or more list element names are duplicated. Correct the names or run as.qlist().")
  
  # do column names have format Cluster?
  if(!all(grepl("Cluster[0-9]+$",unlist(lapply(qlist,colnames))))) warning("is.qlist: One or more qlist dataframes have column names that do not conform to format 'ClusterNumber' like 'Cluster1', 'Cluster12' etc. Correct columns or run as.qlist().\n")
  
  # are attributes present?
  a <- lapply(qlist,function(x) as.matrix(unlist(attributes(x))))
  if(any(!sapply(a,function(x) any(grepl("ind",rownames(x)))))) {
    warning("is.qlist: Attribute 'ind' is missing in one or more runs. Run as.qlist().\n")
  }
  if(any(!sapply(a,function(x) any(grepl("k",rownames(x)))))) {
    warning("is.qlist: Attribute 'k' is missing in one or more runs. Run as.qlist().\n")
  }
}

.readQ <- function(files=NULL,filetype="auto",indlabfromfile=FALSE,readci=FALSE) {
  
  if(is.null(files) || (length(files)==0)) stop("readQ: No input files.")
  if(!is.character(files)) stop("readQ: Argument 'files' is not a character datatype.")
  flen <- length(files)
  
  len <- length(files)
  dlist <- vector("list")
  for (i in seq_along(files))
  {
    # check file
    if(filetype=="auto") 
    {
      chk <- tolower(.checkQ(files[i])$type)
      
      if(chk %in% c("structure","tess","baps","basic","clumpp")) 
      {
        if(chk=="structure") dfr <- .readQStructure(files[i],indlabfromfile=indlabfromfile,readci=readci)
        if(chk=="tess") dfr <- .readQTess(files[i])
        if(chk=="basic") dfr <- .readQBasic(files[i])
        if(chk=="clumpp") dfr <- .readQClumpp(files[i])
        if(chk=="baps") dfr <- .readQBaps(files[i])
        dlist <- append(dlist,dfr)
      }else{
        warning(paste0("readQ: Input file ",files[i]," was not identified as a STRUCTURE, TESS, BAPS, BASIC or CLUMPP filetype. Specify 'filetype' manually or check input.\n"))
      }
    }else{
      if(filetype=="structure") dfr <- .readQStructure(files[i],indlabfromfile=indlabfromfile,readci=readci)
      if(filetype=="tess") dfr <- .readQTess(files[i])
      if(filetype=="basic") dfr <- .readQBasic(files[i])
      if(filetype=="clumpp") dfr <- .readQClumpp(files[i])
      if(filetype=="baps") dfr <- .readQBaps(files[i])
      dlist <- append(dlist,dfr)
    }
  }
  return(dlist)
}

.tabulateQ <- function(qlist,sorttable=TRUE,writetable=FALSE,exportpath=NULL) {
  
  # check input
  .is.qlist(qlist)
  if(!is.logical(writetable)) stop("tabulateQ: Argument 'writetable' not set correctly. Set as TRUE or FALSE.")
  if(writetable) {
    # check exportpath
    if(is.null(exportpath)) stop("tabulateQ: Argument 'exportpath' not set. To use current working directory, set 'exportpath=getwd()'.")
  }
  if(!is.logical(sorttable)) stop("tabulateQ: Argument 'sorttable' not set correctly. Set as TRUE or FALSE.")
  
  # get filenames from selection
  filenames <- names(qlist)
  if(is.null(filenames)) filenames <- paste0("sample",1:length(qlist))
  #number of files selected
  flen <- length(filenames)
  
  # make dataframe container
  main <- data.frame(file=filenames,k=1:flen,ind=1:flen,stringsAsFactors=FALSE)
  
  # loop to make dataframe with filenames and other variables
  # initialise variables
  tq_k <- vector(length=flen,mode="numeric")
  tq_ind <- vector(length=flen,mode="numeric")
  tq_loci <- vector(length=flen,mode="numeric")
  tq_burnin <- vector(length=flen,mode="numeric")
  tq_reps <- vector(length=flen,mode="numeric")
  tq_elpd <- vector(length=flen,mode="numeric")
  tq_mvll <- vector(length=flen,mode="numeric")
  tq_vll <- vector(length=flen,mode="numeric")
  tq_gif <- vector(length=flen,mode="numeric")
  tq_rmse <- vector(length=flen,mode="numeric")
  tq_crossentropy <- vector(length=flen,mode="numeric")
  tq_ploidy <- vector(length=flen,mode="numeric")
  
  for (i in seq_along(qlist))
  {
    # read file & error check
    df1 <- qlist[[i]]
    if(!is.data.frame(df1)) stop(paste0("tabulateQ: List item ",i," is not a data.frame object."))
    if(!any(sapply(df1,is.numeric))) stop(paste0("tabulateQ: List item ",i," has non-numeric columns."))
    
    # get k
    tq_k[i] <- ncol(df1)
    # get ind
    tq_ind[i] <- nrow(df1)
    # loci
    tq_loci[i] <- ifelse(is.null(attr(df1,"loci")),NA,attr(df1,"loci"))
    # burnin
    tq_burnin[i] <- ifelse(is.null(attr(df1,"burnin")),NA,attr(df1,"burnin"))
    # reps
    tq_reps[i] <- ifelse(is.null(attr(df1,"reps")),NA,attr(df1,"reps"))
    # elpd
    tq_elpd[i] <- ifelse(is.null(attr(df1,"elpd")),NA,attr(df1,"elpd"))
    # mvll
    tq_mvll[i] <- ifelse(is.null(attr(df1,"mvll")),NA,attr(df1,"mvll"))
    # vll
    tq_vll[i] <- ifelse(is.null(attr(df1,"vll")),NA,attr(df1,"vll"))
    # gif
    tq_gif[i] <- ifelse(is.null(attr(df1,"gif")),NA,attr(df1,"gif"))
    # rmse
    tq_rmse[i] <- ifelse(is.null(attr(df1,"rmse")),NA,attr(df1,"rmse"))
    # crossentropy
    tq_crossentropy[i] <- ifelse(is.null(attr(df1,"crossentropy")),NA,attr(df1,"crossentropy"))
    # ploidy
    tq_ploidy[i] <- ifelse(is.null(attr(df1,"ploidy")),NA,attr(df1,"ploidy"))
  }
  
  # create dataframe
  main <- data.frame(file=filenames,k=tq_k,ind=tq_ind,stringsAsFactors=FALSE)
  if(all(!is.na(tq_loci))) main$loci <- tq_loci
  if(all(!is.na(tq_burnin))) main$burnin <- tq_burnin
  if(all(!is.na(tq_reps))) main$reps <- tq_reps
  if(all(!is.na(tq_elpd))) main$elpd <- tq_elpd
  if(all(!is.na(tq_mvll))) main$mvll <- tq_mvll
  if(all(!is.na(tq_vll))) main$vll <- tq_vll
  if(all(!is.na(tq_gif))) main$gif <- tq_gif
  if(all(!is.na(tq_rmse))) main$rmse <- tq_rmse
  if(all(!is.na(tq_crossentropy))) main$crossentropy <- tq_crossentropy
  if(all(!is.na(tq_ploidy))) main$ploidy <- tq_ploidy
  
  # sort table on K
  if(sorttable) main <- main[with(main,order(ind,k)),]
  
  # write table if opted
  if(writetable) {
    if(as.numeric(file.access(exportpath,2))==-1) stop(paste0("tabulateQ: Directory ",exportpath," has no write permission."))
    write.table(main,file.path(exportpath,"tabulateQ.txt"),quote=FALSE,row.names=FALSE,sep="\t",dec=".")
    message(file.path(exportpath,"tabulateQ.txt exported."))
  }
  
  return(main)
}

.summariseQ <- summarizeQ <- function(data=NULL,writetable=FALSE,exportpath=NULL) {
  
  # does df data contain any data?
  if(is.null(data) || length(data)==0) stop("summariseQ: No input files.")
  if(!is.logical(writetable)) stop("summariseQ: Argument 'writetable' not set correctly. Set as TRUE or FALSE.")
  if(writetable) {
    # check exportpath
    if(is.null(exportpath)) stop("summariseQ: Argument 'exportpath' not set. To use current working directory, set 'exportpath=getwd()'.")
  }
  
  # make sure dataframe
  if(class(data) != "data.frame") stop("summariseQ: Input is not a dataframe.")
  # convert column names to lowercase
  colnames(data) <- tolower(colnames(data))
  # is column k available?
  if(length(grep("k",colnames(data)))==0) stop("summariseQ: Column k not available.")
  # is column ind available?
  if(length(grep("ind",colnames(data)))==0) stop("summariseQ: Column ind not available.")
  
  # check
  #if(nrow(data) < 2) stop("summariseQ: At least 2 runs are required for this function.")
  
  if(all(c("k","ind","loci","elpd") %in% colnames(data)))
  {
    dframe1 <- aggregate(elpd ~ loci + ind + k,data=data,length)
    colnames(dframe1)[4] <- "runs"
    dframe2 <- aggregate(elpd ~ loci + ind + k,data=data,FUN=function(x) c(elpdmean=mean(x,na.rm=TRUE),elpdsd=sd(x,na.rm=TRUE),elpdmin=mean(x,na.rm=TRUE)-sd(x,na.rm=TRUE),elpdmax=mean(x,na.rm=TRUE)+sd(x,na.rm=TRUE)))[,-c(1:3)]
    dframe1 <- cbind(dframe1,dframe2)
  }else{
    dframe1 <- aggregate(file ~ ind + k,data=data[,c("file","k","ind")],length)
    colnames(dframe1)[3] <- "runs"
  }
  
  # write table if opted
  if(writetable) {
    if(as.numeric(file.access(exportpath,2))==-1) stop(paste0("summariseQ: Directory ",exportpath," has no write permission."))
    write.table(dframe1,file.path(exportpath,"summariseQ.txt"),quote=FALSE,row.names=FALSE,sep="\t",dec=".")
    message(file.path(exportpath,"summariseQ.txt exported."))
  }
  
  return(dframe1)
}

.readQStructure <- function(files=NULL,indlabfromfile=FALSE,readci=FALSE) {
  
  if(is.null(files) || (length(files)==0)) stop("readQStructure: No input files.")
  # number of files selected
  flen <- length(files)
  
  #check file
  if(any(.checkQ(files)$type != "STRUCTURE")) stop("readQStructure: Input may be in incorrect format.")
  
  i <- 1
  dlist <- vector("list",length=flen)
  len1 <- length(files)
  for (i in seq_along(files))
  {
    fname <- basename(files[i])
    file1 <- readLines(as.character(files[i]),warn=FALSE)
    
    # find individuals and get number of individuals
    ind <- as.numeric(as.character(gsub("\\D","",grep("\\d individuals",file1,perl=TRUE,ignore.case=TRUE,value=TRUE)[1])))
    if(is.na(ind)) warning(paste0("Number of individuals is NA in file: ",fname,"\n"))
    
    # get value of k & error check
    k <- as.numeric(as.character(gsub("\\D","",grep("\\d populations assumed",file1,perl=TRUE,ignore.case=TRUE,value=TRUE)[1])))
    if(is.na(k)) warning(paste0("Value of K is NA in file: ",fname,"\n"))
    
    # get number of loci & error check
    loci <- as.numeric(gsub("\\D","",grep("\\d loci",file1,perl=TRUE,ignore.case=TRUE,value=TRUE)[1]))
    if(is.na(loci)) warning(paste0("Number of Loci is NA in file: ",files[i],"\n"))
    
    # get burn-in value & error check
    burnin <- as.numeric(gsub("\\D","",grep("\\d Burn-in period",file1,perl=TRUE,ignore.case=TRUE,value=TRUE)[1]))
    if(is.na(burnin)) warning(paste0("Burn-in value is NA in file: ",files[i],"\n"))
    
    # get burn-in value & error check
    reps <- as.numeric(gsub("\\D","",grep("\\d Reps",file1,perl=TRUE,ignore.case=TRUE,value=TRUE)[1]))
    if(is.na(reps)) warning(paste0("Reps value is NA in file: ",files[i],"\n"))
    
    # get est ln prob of data & error check
    elpd <- as.numeric(gsub("=","",gsub("Estimated Ln Prob of Data","",grep("Estimated Ln Prob of Data",file1,perl=TRUE,ignore.case=TRUE,value=TRUE)[1])))
    if(is.na(elpd)) warning(paste0("Estimated Ln Prob of Data is NA in file: ",files[i],"\n"))
    
    # get mn value of ln likelihood & error check
    mvll <- as.numeric(gsub("=","",gsub("Mean value of ln likelihood","",grep("Mean value of ln likelihood",file1,perl=TRUE,ignore.case=TRUE,value=TRUE)[1])))
    if(is.na(mvll)) warning(paste0("Mean value of ln likelihood is NA in file: ",files[i],"\n"))
    
    # get variance of ln likelihood else NA
    vll <- as.numeric(gsub("=","",gsub("Variance of ln likelihood","",grep("Variance of ln likelihood",file1,perl=TRUE,ignore.case=TRUE,value=TRUE)[1])))
    if(is.na(vll)) warning(paste0("Variance of ln likelihood is NA in file: ",files[i],"\n"))
    
    file1 <- file1[grep(".+\\(\\d+\\).+\\:.+",file1)]
    if(length(file1)==0)
    {
      cstart <- charmatch("Inferred ancestry of individuals",file1)
      cend <- charmatch("Estimated Allele Frequencies in each",file1)
      file1 <- file1[(cstart+2):(cend-1)]
    }
    
    file_a <- file1[file1 != ""]
    rm(file1)
    
    # error check
    tc_file_a <- textConnection(file_a)
    file_b <- read.delim(tc_file_a,header=FALSE,sep="",stringsAsFactors=FALSE)
    close(tc_file_a)
    
    suppressWarnings(
      errorcheck <- try(
        file_b[1,as.integer(grep(":",file_b[1,])+1):as.integer(max(grep("^[0-9]|[.]+$",file_b[1,]))),drop=FALSE],
        silent=TRUE)
    )
    rm(file_b)
    
    if(class(errorcheck)=="try-error")
    {
      # using manual substring
      file_a <- gsub("\\([0-9.,]+\\)","",file_a)
      file_b <- gsub(":  ","",substr(file_a,regexpr(":\\W+\\d\\.\\d+",file_a),nchar(file_a)-1))
      file_b <- sub("\\s+$","",sub("^\\s+","",file_b))
      rm(file_a)
      file_c <- as.vector(as.numeric(as.character(unlist(strsplit(file_b," ")))))
      rm(file_b)
      dframe <- as.data.frame(matrix(file_c,nrow=ind,byrow=TRUE),stringsAsFactors=FALSE)
    }else{
      # using textconnection
      tc_file_a <- textConnection(file_a)
      file_b <- read.delim(tc_file_a,header=FALSE,sep="",stringsAsFactors=FALSE)
      close(tc_file_a)
      dframe <- file_b[,as.integer(grep(":",file_b[1,])+1):as.integer(max(grep("^[0-9]|[.]+$",file_b[1,]))),drop=FALSE]
    }
    
    dframe <- as.data.frame(sapply(dframe,as.numeric),stringsAsFactors=FALSE)
    colnames(dframe) <- paste0("Cluster",1:ncol(dframe))
    row.names(dframe) <- 1:nrow(dframe)
    #row.names(dframe) <- sprintf(paste0("%",paste0(rep(0,nchar(nrow(dframe))),collapse=""),nchar(nrow(dframe)),"d"),1:nrow(dframe))
    
    # add labels
    if(indlabfromfile)
    {
      labeldf <- file_b[,(grep("[0-9]",file_b[1,])[1]+1):(grep("[(]",file_b[1,])[1]-1),drop=FALSE]
      
      if(ncol(labeldf) > 1) labeldf <- data.frame(V2=do.call(paste,c(labeldf,sep="_")),stringsAsFactors=FALSE)
      if(nrow(labeldf)==nrow(dframe))
      {
        if(any(duplicated(labeldf[,1]))) 
        {
          warning(paste0("readQStructure: Individual names in file ",fname," not used due to presence of duplicate names.\n"))
        }else{
          row.names(dframe) <- as.character(labeldf[,1])
        }
      }else{
        warning(paste0("readQStructure: Individual names in file ",fname," not used due to incorrect length.\n"))
      }
    }
    
    attr(dframe,"ind") <- nrow(dframe)
    attr(dframe,"k") <- ncol(dframe)
    attr(dframe,"loci") <- loci
    attr(dframe,"burnin") <- burnin
    attr(dframe,"reps") <- reps
    attr(dframe,"elpd") <- elpd
    attr(dframe,"mvll") <- mvll
    attr(dframe,"vll") <- vll
    
    # confidence intervals
    if(readci) {
      cichk <- grep("([0-9.]+,[0-9.]+)",file_b[1,])
      if(length(cichk)!=0) {
        file_b <- apply(file_b[,cichk,drop=FALSE],1,paste0,collapse="")
        file_b <- gsub("[()]","",gsub(")(",",",file_b,fixed=TRUE))
        cframe <- as.data.frame(matrix(as.numeric(unlist(strsplit(file_b,","))),ncol=ncol(dframe)*2,byrow=TRUE))
        colnames(cframe) <- as.vector(t(outer(paste0("Cluster",1:ncol(dframe)),c("L","H"),paste,sep="")))
        row.names(cframe) <- row.names(dframe)
        attr(dframe,"ci") <- cframe
      }else{
        warning("plotQStructure: Confidence intervals could not be read.\n")
      }
    }
    
    dlist[[i]] <- dframe
  }
  
  fnames <- sub(".txt","",basename(files))
  names(dlist) <- fnames
  return(dlist)
}

.readQTess <- function(files=NULL) {
  
  if(is.null(files) || (length(files)==0)) stop("readQTess: No input files.")
  # number of files selected
  flen <- length(files)
  
  # check file
  if(any(.checkQ(files)$type != "TESS")) warning("readQTess: Input may contain incorrect input format.\n")
  
  i <- 1
  dlist <- vector("list",length=flen)
  len1 <- length(files)
  for (i in seq_along(files))
  {
    # read whole file in
    file1 <- readLines(files[i],warn=FALSE)
    
    # extract the cluster table part
    file1 <- file1[3:c(grep("Estimated Allele Frequencies",file1)-1)]
    
    # remove empty lines
    file1 <- file1[file1 != ""]
    
    # create a text connection
    tc_file1 <- textConnection(file1)
    
    # read as a table
    file2 <- read.delim(tc_file1,header=FALSE,sep="\t",stringsAsFactors=FALSE)
    
    # close text connection
    close(tc_file1)
    
    # choose columns 2 to numofcols-2
    dframe <- file2[,2:(ncol(file2)-2)]
    
    # remove temporary files
    rm(file1,file2)
    
    # convert all columns to numeric
    dframe <- as.data.frame(sapply(dframe,as.numeric),stringsAsFactors=FALSE)
    
    # add column names
    colnames(dframe) <- paste0("Cluster",1:ncol(dframe))
    
    # add attributes
    attr(dframe,"ind") <- nrow(dframe)
    attr(dframe,"k") <- ncol(dframe)
    
    # add to list
    dlist[[i]] <- dframe
  }
  
  # add file names as qlist names
  fnames <- sub(".txt","",basename(files))
  names(dlist) <- fnames
  
  return(dlist)
}

.readQBasic <- function(files=NULL) {
  
  if(is.null(files) || (length(files)==0)) stop("readQBasic: No input files.")
  
  # number of files selected
  flen <- length(files)
  
  # check input file type
  chk <- .checkQ(files)
  if(any(chk$type != "BASIC")) warning("readQBasic: Input may be in incorrect format.\n")
  if(any(is.na(chk$subtype))) warning("readQBasic: Input may be in incorrect format.\n")
  
  i <- 1
  dlist <- vector("list",length=flen)
  len1 <- length(files)
  for (i in seq_along(files))
  {
    # read in delimited files
    if(chk$subtype[i]=="SPACE") dframe <- read.delim(files[i],header=FALSE,sep="",dec=".",stringsAsFactors=FALSE)
    if(chk$subtype[i]=="TAB") dframe <- read.delim(files[i],header=FALSE,sep="\t",dec=".",stringsAsFactors=FALSE)
    if(chk$subtype[i]=="COMMA") dframe <- read.delim(files[i],header=FALSE,sep=",",dec=".",stringsAsFactors=FALSE)
    
    # error if columns contain non-numeric
    if(!all(sapply(dframe,is.numeric))) stop("readQBasic: One or more columns are not numeric.")
    colnames(dframe) <- paste0("Cluster",1:ncol(dframe))
    
    # add attributes
    attr(dframe,"ind") <- nrow(dframe)
    attr(dframe,"k") <- ncol(dframe)
    
    dlist[[i]] <- dframe
  }
  
  # remove file name extensions
  fnames <- sub(".txt","",basename(files))
  fnames <- sub(".tsv","",basename(files))
  fnames <- sub(".csv","",basename(files))
  fnames <- sub(".meanQ","",basename(files))
  
  # add file names to qlist
  names(dlist) <- fnames
  
  return(dlist)
}

.readQClumpp <- function(files=NULL) {
  
  if(is.null(files) || (length(files)==0)) stop("readQClumpp: No input files.")
  
  # number of files selected
  flen <- length(files)
  
  # check file
  chk <- .checkQ(files)
  if(any(chk$type != "CLUMPP")) warning("readQClumpp: Input may be in incorrect format.\n")
  
  i <- 1
  k <- 1
  dlist <- vector("list")
  snames <- vector()
  for (i in seq_along(files))
  {
    fname <- gsub(".txt","",basename(files[i]))
    
    df1 <- read.table(files[i],header=FALSE,sep="",dec=".",quote="",stringsAsFactors=FALSE)
    if(class(df1)!="data.frame") stop("readQClumpp: Read error. Check input format.")
    
    df1[,1] <- factor(df1[ ,1])
    indlev <- levels(df1[,1])
    
    # error check
    if((nrow(df1) %% length(indlev)) != 0) stop("readQClumpp: Number of individuals is not a multiple of the total number of rows.")
    
    Ind <- as.numeric(as.character(length(indlev)))
    tempb <- as.numeric(nrow(df1))
    numruns <- as.numeric(tempb/Ind)
    numk <- ncol(df1) - 2
    
    df2 <- data.frame(Num=factor(rep(1:numruns,1,each=Ind)),
                      Ind=factor(rep(1:Ind,numruns)),
                      df1[,2:(numk+1)],stringsAsFactors=FALSE)
    colnames(df2)[3:ncol(df2)] <- paste0("Cluster",1:(ncol(df2)-2))
    
    for(j in 1:numruns)
    {
      dframe <- subset(df2,df2$Num==j)
      dframe$Num <- NULL
      dframe$Ind <- NULL
      snames <- c(snames,paste0(fname,"-",j))
      
      if(!all(sapply(dframe,is.numeric))) stop("readQClumpp: One or more columns are not numeric.")
      
      attr(dframe,"ind") <- nrow(dframe)
      attr(dframe,"k") <- ncol(dframe)
      
      dlist[[k]] <- dframe
      k <- k+1
    }
  }
  
  snames <- sub(".txt","",snames)
  snames <- sub(".tsv","",snames)
  snames <- sub(".csv","",snames)
  snames <- sub(".meanQ","",snames)
  
  names(dlist) <- snames
  return(dlist)
}

.readQTess3 <- function(t3list=NULL,progressbar=FALSE) {
  
  if(is.null(t3list)) stop("readQTess3: Input is empty.")
  if(!any("tess3" %in% class(t3list))) warning("readQTess3: Input cannot be identified as a valid tess3 class object.\n")
  length(t3list)
  
  # initialise loop variables
  len <- length(t3list)
  qlist <- vector("list",length=len)
  if(progressbar) pb <- txtProgressBar(min=0,max=len,style=3)
  
  # loop to read in data
  for(i in seq_along(t3list))
  {
    if(progressbar) setTxtProgressBar(pb,i)
    if(!"tess3.run" %in% names(t3list[[i]])) stop("readQTess3: 'tess3.run' slot not found in list item ",i,".")
    dlist <- t3list[[i]]$tess3.run[[1]]
    dframe <- as.data.frame(dlist$Q,stringsAsFactors=FALSE)
    colnames(dframe) <- paste0("Cluster",1:ncol(dframe))
    
    # attribute values as added if available else set to NA
    if("n" %in% names(dlist)) {attr(dframe,"ind") <- dlist$n} else {attr(dframe,"ind") <-NA}
    if("K" %in% names(dlist)) {attr(dframe,"k") <- dlist$K} else {attr(dframe,"k") <- NA}
    if("L" %in% names(dlist)) {attr(dframe,"loci") <- dlist$L} else {attr(dframe,"loci") <- NA}
    if("gif" %in% names(dlist)) {attr(dframe,"gif") <- dlist$gif} else {attr(dframe,"gif") <- NA}
    if("rmse" %in% names(dlist)) {attr(dframe,"rmse") <- dlist$rmse} else {attr(dframe,"rmse") <- NA}
    if("crossentropy" %in% names(dlist)) {attr(dframe,"crossentropy") <- dlist$crossentropy} else {attr(dframe,"crossentropy") <- NA}
    if("ploidy" %in% names(dlist)) {attr(dframe,"ploidy") <- dlist$ploidy} else {attr(dframe,"ploidy") <- NA}
    
    qlist[[i]] <- dframe
  }
  if(progressbar) close(pb)
  # list items are labelled
  names(qlist) <- paste0("sample",1:len)
  
  return(qlist)
}

.readQBaps <- function(files=NULL) {
  
  if(is.null(files) || (length(files)==0)) stop("readQBaps: No input files.")
  # number of files selected
  flen <- length(files)
  
  # check if file type is BAPS
  if(any(.checkQ(files)$type != "BAPS")) warning("readQBaps: Input may be in incorrect format.\n")
  
  i <- 1
  dlist <- vector("list",length=flen)
  for (i in seq_along(files))
  {
    # read in all lines from file
    file1 <- readLines(files[i],warn=FALSE)
    
    # extract the cluster table part
    file1 <- file1[grep("^1:",file1):length(file1)]
    
    # read table using delimiter : and use column V2
    tc_file1 <- textConnection(file1)
    file2 <- read.delim(tc_file1,sep=":",header=FALSE,stringsAsFactors=FALSE)$V2
    
    # read table using delimiter space
    tc_file2 <- textConnection(file2)
    dframe <- read.delim(tc_file2,sep="",header=FALSE,stringsAsFactors=FALSE)
    
    # close text connections
    close(tc_file1,tc_file2)
    
    # remove temporary objects
    rm(file1,file2)
    
    # convert all columns to numeric
    dframe <- as.data.frame(sapply(dframe,as.numeric),stringsAsFactors=FALSE)
    
    # create valid column names
    colnames(dframe) <- paste0("Cluster",1:ncol(dframe))
    
    # attach attributes to dataframe
    attr(dframe,"ind") <- nrow(dframe)
    attr(dframe,"k") <- ncol(dframe)
    
    # place dataframe in a list
    dlist[[i]] <- dframe
  }
  
  # remove .txt in all file names
  fnames <- sub(".txt","",basename(files))
  # label qlist objects with file names
  names(dlist) <- fnames
  
  return(dlist)
}