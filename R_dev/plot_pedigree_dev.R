#'Make a pedigree plot 
#' 
#'Prepares a pedigree plot from pedigree output generated from run_sequoia 
#'function.
#' 
#'This is uses the packages VisPedigree and Kinship2 to make pedigree plots
#'
#'@param x list produced from run_sequoia command or data input formatted for 
#'visPedigree. Must contain individual id, dam, sire, and sex.
#'@param tidyped character, default TRUE. Cleans the pedigree for plotting with
#'visPedigree package.
#'@param 
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

plot_pedigree <- function(x, tidyped = TRUE, ...){  
  #===================sanity checks===============
  # check that data is in the correct format
  
  # id names cannot be "0", "NA", "*", " ", ","
  
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
  
  if(length(msg) > 0){
    stop(msg)
  }
  
  #===================prep=========================
  
  # well the data needs to be either 1) output from run_sequoia or 2) independent pedigree construction (individual id order needs to be 1 id, )
  # needs to contain id, dam, sires
  # might contain LLR, metadata, etc
  
  #vispedigree has a good cleaning and sorting function 
  
  if(is.snpRdata(x)){
    # run run_sequoia, grab the output as x.
    
    x <- run_sequoia(x) # with whatever arguments
  }
  
  # after the if above, x should be a run_sequioa output object even if it was originally a snpRdata object.
  
  # to grab the parts of the nested pedigree list we want, can use purrr::map
  # e.x. purrr::map(x, c("pedigree", "Pedigree"))
  # this will collapse the list down a bit
  
  
  
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