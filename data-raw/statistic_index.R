#=============possible calculated statistics===============
single.stats <- list(stat = c("ho", "pi", "maf", "private", "association", "hwe"),
                     col_pattern = list("ho", 
                                "pi", 
                                c("maf", "major", "minor", "maj.count", "min.count"),
                                "pa",
                                c("chi_", "p_armitage_", "log_odds_ratio_", "se_", "associated_allele_", "gmmat_"),
                                c("pHWE")))
window.stats <- list(stat = "tajimas_d",
                  col_pattern = list(c("ws.theta", "ts.theta", "D", "n_snps")))
pairwise.stats <- list(stat = "fst",
                    col_pattern = list(c("fst", "nk")))
sample.stats <- list(stat = "het_hom_ratio",
                  col_pattern = list(c("Het/Hom")))
pop.stats <- list(stat = "ne",
               col_pattern = list(c("LDNe")))
other.stats <- list(stat = c("ld", "genetic_distances", "isolation_by_distance", "geographic_distance"),
                    types = list(c("LD"),
                                 c("genetic_distances"),
                                 c("ibd"),
                                 c("geo_dist")),
                    col_patttern = c(NA, NA, NA, NA))


#=============build lists for each stat====================
statistic_index <- vector("list", length(c(single.stats, window.stats, pairwise.stats, sample.stats, pop.stats,
                                           other.stats[[1]])))

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
  statistic_index[[tracker]] <- list(category = "pairwise", types = c("pairwise", "pairwise.window", "weighted.means", "fst.matrix"),
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

usethis::use_data(statistic_index, internal = T, overwrite = T)
