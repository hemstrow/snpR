#===========================statistic_index================
#=============possible calculated statistics===============
single.stats <- list(stat = c("ho", "pi", "he", "maf", "private", "association", "hwe", "random_forest", "genomic_prediction", "fis"),
                     col_pattern = list("ho", 
                                "pi",
                                "he",
                                c("maf", "major", "minor", "maj.count", "min.count"),
                                "pa",
                                c("chi_", "p_armitage_", "log_odds_ratio_", "se_", "associated_allele_", "gmmat_"),
                                c("pHWE"),
                                c("RF_importance", "RF_importance_pvals"),
                                c("gp_effect"),
                                "fis"))
window.stats <- list(stat = "tajimas_d",
                  col_pattern = list(c("ws.theta", "ts.theta", "D", "n_snps")))
pairwise.stats <- list(stat = "fst",
                    col_pattern = list(c("fst", "nk")))
sample.stats <- list(stat = c("het_hom_ratio", "hs"),
                  col_pattern = list("Het/Hom",
                                     "hs"))
pop.stats <- list(stat = "ne",
               col_pattern = list(c("LDNe")))
other.stats <- list(stat = c("ld", "genetic_distances", "isolation_by_distance", "geographic_distance"),
                    types = list(c("LD"),
                                 c("genetic_distances"),
                                 c("ibd"),
                                 c("geo_dist")),
                    col_patttern = c(NA, NA, NA, NA))


#=============build lists for each stat====================
statistic_index <- vector("list", length(c(single.stats$stat, window.stats$stat, pairwise.stats$stat, sample.stats$stat, pop.stats$stat, other.stats$stat)))

tracker <- 1
for(i in 1:length(single.stats[[1]])){
  names(statistic_index)[tracker] <- single.stats[[1]][i]
  statistic_index[[tracker]] <- list(category = "single", types = c("single", "single.window", "bootstraps", "weighted.means"),
                                     col_pattern = single.stats[[2]][[i]])
  tracker <- tracker + 1
}

for(i in 1:length(window.stats[[1]])){
  names(statistic_index)[tracker] <- window.stats[[1]][i]
  statistic_index[[tracker]] <- list(category = "window", types = c("single.window", "bootstraps", "weighted.means"),
                                     col_pattern = window.stats[[2]][[i]])
  tracker <- tracker + 1
}

for(i in 1:length(pairwise.stats[[1]])){
  names(statistic_index)[tracker] <- pairwise.stats[[1]][i]
  statistic_index[[tracker]] <- list(category = "pairwise", types = c("pairwise", "pairwise.window", "weighted.means", "fst.matrix", "bootstraps"),
                                     col_pattern = pairwise.stats[[2]][[i]])
  tracker <- tracker + 1
}

for(i in 1:length(sample.stats[[1]])){
  names(statistic_index)[tracker] <- sample.stats[[1]][i]
  statistic_index[[tracker]] <- list(category = "sample", types = c("sample", "weighted.means"),
                                     col_pattern = sample.stats[[2]][[i]])
  tracker <- tracker + 1
}

for(i in 1:length(pop.stats[[1]])){
  names(statistic_index)[tracker] <- pop.stats[[1]][i]
  statistic_index[[tracker]] <- list(category = "pop", types = c("pop"),
                                     col_pattern = pop.stats[[2]][[i]])
  tracker <- tracker + 1
}

for(i in 1:length(other.stats[[1]])){
  names(statistic_index)[tracker] <- other.stats[[1]][i]
  statistic_index[[tracker]] <- list(category = "other", types = other.stats[[2]][[i]],
                                     col_pattern = other.stats[[3]][i])
  tracker <- tracker + 1
}




#============================test_snps=====================
set.seed(1212)
test_snps <- subset_snpR_data(stickSNPs, 1:12, sample(nsamps(stickSNPs), 10, F))
test_snps <- filter_snps(test_snps)

sample.meta(test_snps)$pop <- rep(c("ASP", "PAL"), 5)
sample.meta(test_snps)$fam <- rep(c("A", "B"), each = 5)

#============================save==========================
.internal.data <- list(test_snps = test_snps, statistic_index = statistic_index)
usethis::use_data(.internal.data, internal = T, overwrite = T)
