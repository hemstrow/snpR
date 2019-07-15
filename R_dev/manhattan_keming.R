

# we want to make a manhattan plot, here's an example:
fsts <- dat3@pairwise.stats

fsts <- as.data.frame(fsts)

fsts$group <- as.numeric(as.factor(fsts$group))
fsts <- na.omit(fsts)
fsts <- fsts[fsts$comparison == "OPL~SMR",]

qqman::manhattan(fsts, chr = "group", bp = "position", snp = "snp", p = "fst", logp = FALSE)


# we could use the above funciton if we can make it work well, otherwise ggplot2 might be a good way to go:
fsts <- dat3@pairwise.stats

fsts <- as.data.frame(fsts)
library(ggplot2)
ggplot(fsts, aes(x = position, y = fst)) + geom_point() +
  theme_bw() + facet_grid(group~comparison)



# we will want to make these plots for multiple different statistics. They are stored in snpRdata objects.
# making an object.
dat3 <- stickSNPs
dat3 <- import.snpR.data(genotypes = dat3[,-c(1:3)], snp.meta = dat3[,1:3], sample.meta =
                         data.frame(pop = substr(colnames(dat3)[-c(1:3)], 1, 3), # just getting sample meta data
                                    fam = rep(c("A", "B", "C", "D"), length = ncol(dat3) - 3),
                                    stringsAsFactors = F), mDat = "NN")


# snpRdata objects store data in sockets, accessed by the @ symbol:
dat3 <- calc_ho(dat3, "pop") # add observed heterozygosity, splitting by population
dat3 <- calc_pairwise_fst(dat3, "pop") # add pairwise FST, splitting by population

head(dat3@stats) # see ho
head(dat3@pairwise.stats) # see fst
# all of the stats we'll want to plot with a manhattan plot will be in these two slots!
# one note is the these are data.table objects, from the data.table package. You can convert these back to
# data frames with as.data.frame()

tplot <- make_manhattan(dat3, "pop", "fst")
tplot <- make_manhattan(dat3, "pop.fam", "ho")
