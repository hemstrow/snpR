## ----include=FALSE-------------------------------------------------------
library(snpR); library(ggplot2)

## ----echo=FALSE----------------------------------------------------------
print(head(stickSNPs)[,1:6], row.names = F)

## ------------------------------------------------------------------------
pops <- table(substr(colnames(stickSNPs[,4:ncol(stickSNPs)]), 1, 3)) # this produces a table of the number of individuals in each population. Individuals must be sorted in the order they are given in this table!

pops



## ---- results = "hide"---------------------------------------------------
filt_snps <- filter_snps(x = stickSNPs,
                         ecs = 3, # how many metadata columns are there?
                         maf = 0.05, # what is the minimum minor allele frequency to keep?
                         hf_hets = 0.55, # what is the maximum heterozygote frequency to keep?
                         min_ind =  (ncol(stickSNPs) - 3)/2, # remove all loci sequenced at less than half of the individuals.
                         min_loci = .5, # remove all individuals sequenced at less than half of the loci
                         re_run = "partial", # tell filter_snps to recheck that all loci are polymorphic after removing individuals. Could set to "full" to re-check maf and heterozygote frequency, or FALSE to skip.
                         pop = pops, # tell filter_snps that we want to keep all loci with a maf above 0.1 in at least one population.
                         non_poly = T, # we want to remove non-polymorphic loci
                         bi_al = T, # we want to remove non bi-allelic loci
                         mDat = "NN") # missing genotypes are listed as "NN"

## ---- results="hide"-----------------------------------------------------
# Convert to allele presence/absence
pa <- format_snps(filt_snps, ecs = 3, output = "pa") 

## ---- echo=FALSE---------------------------------------------------------
head(pa)[,1:10]

## ------------------------------------------------------------------------
# add population info
pa <- cbind(pop = substr(pa$samp, 1, 3), pa)

head(pa)[,1:10]


## ------------------------------------------------------------------------
# make a PCA plot
pca <- PCAfromPA(x = pa, ecs = 2, 
                 plot.vars = "pop", # name of the columns to plot in different colors. Can give up to 2.
                 c.dup = T, # should the dataset be checked for any duplicate individuals? Usually not needed, but a good check to run.
                 do.plot = T)

# view the plot
pca$plot

## ------------------------------------------------------------------------
# make a vector of the counts of sequenced loci per individual
counts <- colSums(ifelse(filt_snps[,-c(1:3)] == "NN", 0, 1))


## ------------------------------------------------------------------------
p <- PCAfromPA(x = pa, ecs = 2, plot.vars = "pop", c.dup = T, 
                 mc = 600, # minimum count of sequenced loci to plot
                 counts = counts) # provide the vector of counts.

# view the plot, should be the same!
pca$plot

## ---- results="hide"-----------------------------------------------------
# do a t-SNE plot

tsne <- tSNEfromPA(x = pa, ecs = 2, plot.vars = "pop", c.dup = T,
                   initial_dims = 100, # how many initial PCA dimensions should the t-SNE consider?
                   iter = 500) # how many optimization iterations should be run?


## ---- echo=FALSE---------------------------------------------------------
tsne$plot

## ---- results="hide"-----------------------------------------------------
# convert to ac. Note that we need to give it the pop information from earlier!
ac <- format_snps(filt_snps, 3, output = "ac", pop = pops)

## ---- echo=FALSE---------------------------------------------------------
# view:
head(ac)


## ------------------------------------------------------------------------
# calculate pi
pi <- calc_pi(ac, ecs = 3)

# the result:
head(pi)


## ------------------------------------------------------------------------
# calculate Ho, note that we again provide pop info.
ho <- calc_Ho(x = filt_snps, ecs = 3, pop = pops)

# the result:
head(ho)

## ------------------------------------------------------------------------
# check for private alleles, uses the ac format:
pr.a <- check_private(x = ac, ecs = 3)

# view the result:
head(pr.a)


## ------------------------------------------------------------------------
pairwise.fst <- calc_pairwise_Fst(x = ac, ecs = 3, do.nk = T, # should we also count up the number of sequenced individuals for each locus in each pair of populations?
                                  method = "WC", # specifying that we want to use the Wier and Cockerham 1984 method.
                                  char.dat = filt_snps, pop = pops, c.d.ecs = 3) # we can re-provide this info so that the function can call calc_Ho for us. Alternatively, we could bind a Ho column to the dataset, in which case this can be skipped.

# View the Fst results:
head(pairwise.fst$FST)

# View the number of sequenced individuals per locus (nk) in each pair of populations
head(pairwise.fst$nk)

## ---- results="hide"-----------------------------------------------------
# collect pi, ho, and nk data.

# first, we'll melt the observed heterozygosity data into long form
mho <- reshape2::melt(ho, id.vars = colnames(ho)[1:3])
colnames(mho)[4:5] <- c("pop", "Ho")

# then we'll bind things together. Everything should be in the same order, but the data could be sorted or merged to make sure if wanted. nk  can be taken from our ac data:
s.dat <- cbind(mho, pi = pi$pi, nk = ac$n_total)


# finally, we can smooth
s_piho <- s_ave_multi(s.dat, parms = c("Ho", "pi"), # which variables to we smooth? 
                    sigma = 200, # how big are our windows, in kb?
                    ws = 50, # how much do we slide each window, in kb?
                    nk = T, # are we weighting by nk?
                    levs = c("group", "pop")) # which levels are we splitting the smoothing by?


## ------------------------------------------------------------------------
head(s_piho)

## ---- results="hide"-----------------------------------------------------
# First, melt the data down:

mfst <- reshape2::melt(pairwise.fst$FST, id.vars = colnames(pairwise.fst$FST)[1:3])
colnames(mfst)[4:5] <- c("pop", "FST")

# next, melt the nk data down:
mnk <- reshape2::melt(pairwise.fst$nk, id.vars = colnames(pairwise.fst$nk)[1:3])
colnames(mnk)[4:5] <- c("pop", "nk")

# bind them together:
mfst <- cbind(mfst, nk = mnk$nk)

# return pop to a character variable:
mfst$pop <- as.character(mfst$pop)

# and smooth
s_FST <- s_ave_multi(mfst, parms = "FST", sigma = 200, ws = 50, nk = T, levs = c("group", "pop"))

## ------------------------------------------------------------------------
head(s_piho)

## ------------------------------------------------------------------------
# Plot Ho, pi works the same way!
ggplot(s_piho, aes(x = position, y = smoothed_Ho, color = pop)) + geom_line() + theme_bw() +
  facet_wrap(~group, scales = "free_x") # facet wrap on group to split by linkage group/chromosome.

## ------------------------------------------------------------------------

ggplot(s_FST[s_FST$group == "groupIX",], aes(x = position, y = smoothed_FST, color = pop)) + geom_line() + theme_bw()

## ---- results="hide"-----------------------------------------------------
# subset out the OPL data
OPL <- cbind(filt_snps[,1:3], filt_snps[,grepl("OPL", colnames(filt_snps))])

# run the LD function.
LD <- LD_full_pairwise(OPL, ecs = 3, levels = "group")

## ------------------------------------------------------------------------
# View the proximity table (using tail instead of head, since head has a bunch of rows with missing LD scores).
tail(LD$prox)

ggplot(LD$prox, aes(x = proximity, y = rsq)) + geom_point() + theme_bw()

## ------------------------------------------------------------------------
# plot LD across LG IX:
LD_plot <- LD_pairwise_heatmap(LD$groupIX$rsq)

# View the result
LD_plot$plot


