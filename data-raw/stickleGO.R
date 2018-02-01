############################################################
#prepare sample GO data

#import some real data
stickleGO <- read.table("~/Stickleback/Data/StickleGO_genes_positions.txt", header = T, stringsAsFactors = F)

#subset and rename to be compatible with randPI
stickleGO <- stickleGO[stickleGO$group %in% unique(stickleGO$group)[1:4],] #first four groups
stickleGO <- stickleGO[stickleGO$end <= max(randSMOOTHed$position),] #short enough to be on the mock chromosomes
stickleGO$group <- ifelse(stickleGO$group == "groupI", "chr1",
                          ifelse(stickleGO$group == "groupII", "chr2",
                                 ifelse(stickleGO$group == "groupIII", "chr3", "chr4")))

#save
devtools::use_data(stickleGO, overwrite = T)
