#prepare output with a bunch of statistics for examples.

#get pop info from stickSNPs
pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3))
l <- list(c(names(pops)), as.numeric(pops))

#get alternative formats from stickSNPs
ac <- format_snps(stickSNPs, 3, pop = l)

#pi
pi <- calc_pi(ac)
pi <- cbind(ac[,1:4], pi)

#ho
ho <- calc_Ho(stickSNPs, 3, pop = l)
ho <- ho[,4:ncol(ho)]
ho <- reshape2::melt(ho)
colnames(ho) <- c("pop", "Ho")

#pa
pa <- check_private(ac)
colnames(pa) <- l[[1]]
pa <- reshape2::melt(pa)
colnames(pa) <- c("pop","pa")

#bind together and save
stickSTATs <- cbind(pi, Ho = ho[,2], pa = pa[,2], nk = ac$n_total)
devtools::use_data(stickSTATs, overwrite = T)
