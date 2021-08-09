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

plot_pedigree <- function(x, tidyped = TRUE, ...)

{  
#===================sanity checks===============
# check that data is in the correct format

# id names cannot be "0", "NA", "*", " ", ","

# ensure necessary packages are installed
check.installed(visPedigree) #luansheng/visPedigree
check.installed(kinship2) #available via cran

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
}