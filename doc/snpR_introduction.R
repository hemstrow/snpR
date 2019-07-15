## ----include=FALSE-------------------------------------------------------
library(snpR); library(ggplot2)

## ----echo=FALSE----------------------------------------------------------
print(head(stickFORMATs$NN)[,1:6], row.names = F)

## ------------------------------------------------------------------------
# some example data is present in stickRAW
# the first three characters of the column names are the population, so let's grab that:
# since the first three columns contain snp.metadata, we need to remove those!
sample_meta <- data.frame(pop = substr(colnames(stickRAW)[-c(1:3)], 1, 3), stringsAsFactors = F)

# grab our sample metadata
snp_meta <- stickRAW[,1:3]

# import, remember to remove metadata from the genotypes!
my.dat <- import.snpR.data(stickRAW[,-c(1:3)], 
                           snp.meta = snp_meta, 
                           sample.meta = sample_meta)


## ------------------------------------------------------------------------
# Import numeric formatted data.
head(stickFORMATs$`0000`)[,1:6]

# get sample metadata, as before.
sample_meta <- data.frame(pop = substr(colnames(stickFORMATs$`0000`)[-c(1:4)], 1, 3), stringsAsFactors = F)

# import
my.dat2 <- format_snps(stickFORMATs$`0000`, input_format = "0000", input_meta_columns = 4, input_mDat = "0000", sample.meta = sample_meta)
#

## ------------------------------------------------------------------------
head(my.dat@geno.tables$as)

## ------------------------------------------------------------------------
# calculate observed heterozygosity, overwritting the my.dat object
my.dat <- calc_ho(my.dat)


## ------------------------------------------------------------------------
# fetch calculated statistics
ho <- get.snpR.stats(my.dat)
head(ho)


## ---- results = "hide"---------------------------------------------------
my.dat <- filter_snps(x = my.dat,
                         maf = 0.05, # what is the minimum minor allele frequency to keep?
                         hf_hets = 0.55, # what is the maximum heterozygote frequency to keep?
                         min_ind =  .5, # remove all loci sequenced at less than half of the individuals.
                         min_loci = .5, # remove all individuals sequenced at less than half of the loci
                         re_run = "partial", # tell filter_snps to recheck that all loci are polymorphic after removing individuals. Could set to "full" to re-check maf and heterozygote frequency, or FALSE to skip.
                         maf.facets = "pop", # we want to keep any SNPs that pass the maf filter in ANY population!
                         non_poly = T, # we want to remove non-polymorphic loci
                         bi_al = T) # we want to remove non bi-allelic loci") # missing genotypes are listed as "NN"

## ------------------------------------------------------------------------
head(stickSNPs@sample.meta)

## ------------------------------------------------------------------------
# generate both tSNE and PCA plots
p <- plot_clusters(stickSNPs, facets = c("pop.fam"))

## ------------------------------------------------------------------------
# View the PCA
p$plots$pca


## ------------------------------------------------------------------------
# view the tSNE
p$plots$tsne

## ------------------------------------------------------------------------
# calculate minor allele frequencies, pi, ho, check for private alleles and HWE.
my.dat <- calc_maf(my.dat, "group.pop")
my.dat <- calc_pi(my.dat, "group.pop")
my.dat <- calc_ho(my.dat, "group.pop")
my.dat <- calc_private(my.dat, "group.pop")
my.dat <- calc_hwe(my.dat, "group.pop")


## ------------------------------------------------------------------------
stats <- get.snpR.stats(my.dat, "all") # using the "all" facet requests statistics for all facets that have been run.

head(stats)

## ------------------------------------------------------------------------
my.dat2 <- calc_ho(stickSNPs, "pop.fam")

stats <- get.snpR.stats(my.dat2, "pop.fam")

head(stats)

## ------------------------------------------------------------------------
my.dat <- calc_pairwise_fst(my.dat, "pop.group")

## ------------------------------------------------------------------------
stats <- get.snpR.stats(my.dat, "pop.group", type = "pairwise")

head(stats)

## ---- results='hide'-----------------------------------------------------

my.dat <- calc_smoothed_averages(my.dat, "pop.group",
                                 sigma = 200, # using a window size of 200
                                 step = 50) # using a step size of 50 between windows



## ------------------------------------------------------------------------
stats <- get.snpR.stats(my.dat, "pop.group", "single.window")

head(stats)


## ------------------------------------------------------------------------
stats <- get.snpR.stats(my.dat, "pop.group", "pairwise.window")

head(stats)


## ------------------------------------------------------------------------
stats <- get.snpR.stats(my.dat, "pop.group", "single.window")

# Plot pi, pi works the same way!
ggplot(stats, aes(x = position, y = pi, color = subfacet)) + geom_line() + theme_bw() +
  facet_wrap(~snp.subfacet, scales = "free_x") # facet wrap on group to split by linkage group/chromosome.

## ------------------------------------------------------------------------

stats <- get.snpR.stats(my.dat, "pop.group", "pairwise.window")


ggplot(stats[stats$snp.subfacet == "groupXIX",], aes(x = position, y = fst, color = subfacet)) + geom_line() + theme_bw()

## ---- results="hide"-----------------------------------------------------
list(pop = "ASP", group = "groupIX")


## ------------------------------------------------------------------------
list(pop = c("SMR", "OPL"), group = c("groupIX", "groupIV"))


## ------------------------------------------------------------------------
list(fam.pop = "A.ASP")

## ------------------------------------------------------------------------
my.dat <- calc_pairwise_ld(my.dat, facets = "group.pop", 
                           subfacets =  list(pop = c("ASP", "OPL"), group = "groupIX"))


## ------------------------------------------------------------------------
stats <- get.snpR.stats(my.dat, "group.pop", "LD")

str(stat)

## ------------------------------------------------------------------------
head(stats$prox)

## ------------------------------------------------------------------------
stats$matrices$ASP$groupIX$rsq[1:10,1:10]

## ------------------------------------------------------------------------
p <- plot_pairwise_LD_heatmap(my.dat, "group.pop")

p$plot

## ------------------------------------------------------------------------
p <- plot_pairwise_LD_heatmap(my.dat, "group.pop", sample.subfacet = "ASP")

p$plot

## ------------------------------------------------------------------------
my.dat <- do_bootstraps(my.dat, "group.pop", boots = 1000, sigma = 200, step = 50, statistics = "all")

stats <- get.snpR.stats(my.dat, "group.pop", "single.window")
head(stats)

## ---- results='hide'-----------------------------------------------------
my.dat2 <- calc_basic_snp_stats(stickSNPs, "group.pop", sigma = 200, step = 50)


## ---- results='hide'-----------------------------------------------------
my.dat2 <- calc_basic_snp_stats(stickSNPs, "group.pop", sigma = 200, step = 50, par = 4)

