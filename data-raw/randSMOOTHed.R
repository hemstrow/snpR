#generate mock smoothed data from two pops and two chromosomes.
###########################################################################
#set metadata

#positions, groups, pops
groups <- c("chr1", "chr2", "chr3", "chr4")
pos <- sample(10000000, 2000)
pos <- data.frame(position = rep(pos, 4),
                  pop = c(rep("A", 2000), rep("B", 2000), rep("C", 2000), rep("D", 2000)),
                  group = c(rep(groups, 2000)), stringsAsFactors = F)

#sort these:
pos <- dplyr::arrange(pos, pop, group, position)
pos$nk <- floor(rnorm(100, 40, 5))
pos$nk[pos$nk <= 0] <- 1


##########################################################################
#create pi

#background pi, 500 for each chr
backpi <- data.frame(group = c(rep(groups[1], 500),
                  rep(groups[2], 500),
                  rep(groups[3], 500),
                  rep(groups[4], 500)),
                pi = rnorm(2000, .2, .1), stringsAsFactors = F)


#randomize for each pop within this
pos$pi <- numeric(nrow(pos))
for(i in 1:length(groups)){
  pos$pi[pos$group == groups[i]] <- rep(backpi[backpi[,1] == groups[i],2], 4) + rnorm(2000, .05, .025)
}


#randomly adjust up in pop A
pos$pi[pos$pop == "A"] <- pos$pi[pos$pop == "A"] + rnorm(2000, .025, .005)

#randomly reduce pi between 200000 and 400000 on chr 3 in pop B
pos$pi[pos$pop == "B" & pos$position >= 2000000 & pos$position <= 2750000 & pos$group == "chr3"] <- pos$pi[pos$pop == "B" & pos$position >= 2000000 & pos$position <= 2750000 & pos$group == "chr3"] - rnorm(length(pos$pi[pos$pop == "B" & pos$position >= 2000000 & pos$position <= 2750000 & pos$group == "chr3"]), .075, .005)

#fix negative pi
pos$pi[pos$pi < 0] <- 0

##########################################################################
#generate smoothed pi values and window snps

randSMOOTHed <- run_gp(pos, smoothed_ave, parameter = "pi", sigma = 200, nk_weight = TRUE, fixed_window = 50)

##########################################################################
#generate null distribution
randPIBOOTS <- resample_long(randPI[randPI$pop == "A",], "pi", 100, 200, TRUE, randSMOOTHed[randSMOOTHed$pop == "A",], TRUE, 10)

##########################################################################
#save randPI, randSMOOTHed, and randPIBOOTs
randPI <- pos

devtools::use_data(randPI, overwrite = T)
devtools::use_data(randSMOOTHed, overwrite = T)
devtools::use_data(randPIBOOTS, overwrite = T)
