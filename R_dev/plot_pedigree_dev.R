#'Make a pedigree plot 
#' 
#'Prepares a pedigree plot from pedigree output generated from run_sequoia 
#'function.
#' 
#'This is uses the packages VisPedigree and Kinship2 to make pedigree plots
#'
#'@param x list produced from run_sequoia command or data input formatted for 
#'visPedigree. Must contain individual id, dam, sire, and sex.
#'@param plot.type "visped" or "kinship2, default "visped". Which pedigree 
#'plotting package should be used?  
#'@param facets character, default NULL. Categorical metadata variables by which
#'  to break up plots. See \code{\link{Facets_in_snpR}} for more details.
#' 
#'@return 
#' 
#'@export
#'@author William Hemstrom
#'@author Melissa Jones
#'
#'@references Sheng, Luan. (2018) visPedigree: A package for tidying and drawing animal pedigree. URL https://github.com/luansheng/visPedigree.
#'@references Sinwnwell, J.P., Therneau, T.M., and Schaid, D.J. (2014) The kinship2 R Package for Pedigree Data. Hum Hered., 78, 91-93. 
#' 
#'@examples
#'# output from run_sequoia 
#'# output from COLONY Best Cluster
#'# example pedigree stickPED

plot_pedigree <- function(x, facets = NULL, plot.type = "visped", run_dupcheck = TRUE, na.omit = TRUE,...){  
  #===================sanity checks===============
  # check that data is in the correct format
  # snpr object then run sequoia
  # list - then extract and proceed
  # option for filtering based on LLR scores? or just have people filter the output first and then go through with the plotting?
  
  # ensure necessary packages are installed
  check.installed("visPedigree", install.type = "github", "luansheng/visPedigree") #luansheng/visPedigree
  check.installed("kinship2") #available via cran
  if(is.snpRdata(x)){
    check.installed("sequoia")
  }
  
  msg <- character(0)
  
  if(!all(c("sire", "dam", "id", "sex") %in% colnames(sample.meta(x)))) {
    warning("Need columns sire, dam, id, and sex in the dataset. \n")
  }
  
  if(x$id %in% c("0", "NA", "*", " ", ",")) {
    stop("Sample id cannot be 0, NA, asterisk, blank space, or comma. \n")
  }
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  #===================prep=========================
  
  # well the data needs to be either 1) output from run_sequoia or 2) independent pedigree construction (individual id order needs to be 1 id, 2 sire, 3 dam)

  err.msg <- "Input data (x) in unrecognized format. See documentation for accepted formats.\n"
  if(is.snpRdata(x)){
    # run run_sequoia, grab the output as x.
    x <- run_sequoia(x, facets = facets, run_parents = TRUE, run_pedigree = TRUE, run_dupcheck = run_dupcheck,
                     min_maf = 0.3, min_ind = 0.5, ...) #this can take some time
  }
  else if(is.list(x)){
    if(!sum(unlist(lapply(lapply(x, names), function(x) "pedigree" %in% x))) == length(x) & length(x) > 0){
      stop(err.msg)
    }
  }
  else if(is.data.frame(x) | is.matrix(x)){
    if(all(c("id", "dam", "sire") %in% colnames(x))){
      # code to to make this look like the output from run_sequoia
      x <- list(.base_.base_.base_.base = list(pedigree = list(Pedigree = x)))
    }
    else{
      stop(err.msg)
    }
  }
  else if(is.character(x[[1]]) & length(x) == 1){
    if(file.exists(x)){
      # coerce something like stick ped up to the same format that run_sequoia returns
      x <- data.table::fread(x)
      
      # make this look like stickPED, since that's easier
      
      # re-run from the top
      return(plot_pedigree(x, plot.type = plot.type, ...))
    }
    else{
      stop(err.msg)
    }
  }
  
  # to grab the parts of the nested pedigree list we want, can use purrr::map
  # e.x. purrr::map(x, c("pedigree", "Pedigree"))
  # this will collapse the list down a bit
  

  x <- purrr::map(x, c("pedigree", "Pedigree")) #but if data is from colony it would probably be a data.frame    
  #if it was run with facets there will be sub-lists
  
  
  # initialize
  out <- vector("list", length(x)) #what do we actually want to have returned? a list would be ok, and I think the forma that run_sequoia will return anyways ...
  names(out) <- names(x)
  
  #==================run===========================       
  for(i in 1:length(out)){
    #grab the pedigree information from the appropriate list part - and currently without filtering inds (via LLR -seq or prob - col)
    x2 <- x[[i]]$pedigree$Pedigree  # basically want to grab the pedigree results from the output for each facet and then run it through and save in something later but this isn't quite the right way to do this
    #how to grab the first elements of the list?
    data <- x2[,c(1,3,2)] #need to find (id, dam, and sire) and rearrange the order to feed into visped - so that the order is id, sire, dam instead of the sequoia output format id, dam, sire
    bads <- which(rowSums(is.na(data[,-1])) == 2)
    if(length(bads) > 0){
      if(length(bads) == nrow(data)){
        # could throw a warning
        next
      }
      
      # clean data of NA inds if na.omit = TRUE.
      else if(na.omit){
        data <- data[-bads,]
      }
    }
    
    
    tp <- visPedigree::tidyped(data) #but want to save each version so they can be accessed somehow later
    visPedigree::visped(tp) #want to somehow add plotting options (colors, pruning)
    tplot <- grDevices::recordPlot()
    
    #but need to add in options for colors ... and more 
    
    ###somehow specify plotting visPedigree vs kinship2 and respond differently - and different if sequoia or colony based input
    # after prepping the plot, add to output
    out[[i]] <- list(plot = tplot, data = tp)
  }
  return(out)
}

### need to remember to have option to remove individuals without family assigned (mum or dad NA) when plotting so that it looks cleaner   

###now there is a fake pedigree


###Will's from 8 Aug 2021

# pedigree <- ped$.base_.base_.base_.base$pedigree$Pedigree
# pedigree <- pedigree[,c(1,3,2)]
# pedigree$Sex <- NA
# pedigree$Sex[1:nrow(ped$.base_.base_.base_.base$pedigree$LifeHistSib)] <- ped$.base_.base_.base_.base$pedigree$LifeHistSib$Sex
# pedigree$Sex[(nrow(ped$.base_.base_.base_.base$pedigree$LifeHistSib)+ 1):nrow(pedigree)] <-
#   substr(pedigree$id[(nrow(ped$.base_.base_.base_.base$pedigree$LifeHistSib)+ 1):nrow(pedigree)], 1, 1)
# pedigree$Sex[pedigree$Sex == "F"] <- 1
# pedigree$Sex[pedigree$Sex == "M"] <- 2
# pedigree$Sex[pedigree$Sex == "U"] <- 3
# pedigree$Sex <- as.numeric(pedigree$Sex)
# 
# newinds <- pedigree[(nrow(ped$.base_.base_.base_.base$pedigree$LifeHistSib)+ 1):nrow(pedigree),]
# oldinds <-  pedigree[1:nrow(ped$.base_.base_.base_.base$pedigree$LifeHistSib),]
# 
# bads <- rowSums(is.na(oldinds[,2:3]))
# bads <- which(bads == 2)
# tp <- rbind(oldinds[-bads,], newinds)
# tp <- visPedigree::tidyped(pedigree[-bads,])
# tp$Sex <- ifelse(tp$Sex == 1, "female", ifelse(tp$Sex == 2, "male", "unknown"))
# visPedigree::visped(tp)
# 
# 
# bads <- rowSums(is.na(tp[,2:3]))
# bads <- which(bads == 1)
# kn <- tp[-bads,]
# 
# kinship2::plot.pedigree(kinship2::pedigree(kn$Ind, kn$Sire, kn$Dam, kn$Sex))


# cod <- data.table::fread("../N_H_etal_LCT/colony/str_SLK_2015/colony/colony_input.BestCluster", header = T)
# cod$Sex <- sample(c("male", "female", "unknown"), nrow(cod), T)
# 
# newmales <- data.frame(id = unique(cod$FatherID), sire = NA, dam = NA, sex = "male")
# newfemales <- data.frame(id = unique(cod$MotherID), sire = NA, dam = NA, sex = "female")
# 
# cod <- rbind(cod[, c(3:6)], newmales, newfemales, use.names = FALSE)
# 
# kinship2::plot.pedigree(kinship2::pedigree(id = cod$OffspringID, dadid = cod$FatherID, momid = cod$MotherID, sex = cod$Sex))