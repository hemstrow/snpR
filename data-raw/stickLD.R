##############################################################
#Get a subset of stickle data to use as an LD sample

#import
stickLD <- read.table("~/Stickleback/Full Data/2017_reruns/PlotData/LD/full_LDIX_ASP_rsq.txt", header = T, check.names = F)

#pick sites
sites <- unique(c(as.character(stickLD[,1]), as.character(colnames(stickLD))))
sites <- sites[sites != "V1"]
sites <- sites[sites != "NA"]
sites <- sites[seq(1,length(sites),10)]

stickLD <- stickLD[stickLD[,1] %in% sites, c(1, which(colnames(stickLD) %in% sites))]

devtools::use_data(stickLD)
