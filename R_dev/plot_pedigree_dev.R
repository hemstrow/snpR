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
#'# output from other pedigree construction programs

plot_pedigree <- function(x, plot.type = "visped", facets = ".base", ...){  
  #===================sanity checks===============
  # check that data is in the correct format
  # snpr object then run sequoia
  # list - then extract and proceed
  
  
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

  if(is.snpRdata(x)){
    # run run_sequoia, grab the output as x.
    
    x <- run_sequoia(x, run_parents = TRUE, run_pedigree = TRUE, run_dupcheck = TRUE, min_maf = 0.3, min_ind = 0.5, ...) #this can take some time
  }
  
  # after the if above, x should be a run_sequoia output object even if it was originally a snpRdata object. 
  ## shouldn't it get saved to something new, instead of overwriting the snpR object? eg. x2, but then it makes keeping track of the facet options difficult. maybe later
  
  # to grab the parts of the nested pedigree list we want, can use purrr::map
  # e.x. purrr::map(x, c("pedigree", "Pedigree"))
  # this will collapse the list down a bit
  
  if(!is.snpRdata(x)){ #actually run_sequoia will result in a list - so then x will be a list
    if(class(x) == "list"){

            x <- purrr::map(x, c("pedigree", "Pedigree")) #but if fed in from colony it would be a data.frame..?    
    #if it was run with facets there will be sublists
            n.tasks <- length(x) # is there some way to grab the names to make specific named tasks? <- yes in sequoia_interface
            #how to access the list elements ? maybe more purrr maps?

#==================run===========================       
for(i in 1:n.tasks){
  #grab the pedigree information from the appropriate list part - and currently without filtering inds (via LLR -seq or prob - col)
  x2 <- purrr::map(.x = x, .f = [[i]])  # basically want to grab the pedigree results from the output for each facet and then run it through and save in something later but this isn't quite the right way to do this
  #how to grab the first elements of the list?
  data <- x2[] #need to find (id, dam, and sire) and rearrange the order to feed into visped - so that the order is id, sire, dam 
  tp <- visPedigree::tidyped(data) #but want to save each version so they can be accessed somehow later
  visPedigree::visped(tp) #want to somehow add plotting options (colors, pruning)
  
###somehow specify plotting visPedigree vs kinship2 and respond differently - and different if sequoia or colony based input
  
}
            
            }
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
  
}