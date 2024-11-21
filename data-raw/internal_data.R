#=============statistic_index================
#=============possible calculated statistics===============
single.stats <- list(stat = c("ho", "pi", "he", "maf", "private", "association", "hwe", "random_forest", "genomic_prediction", 
                              "fis", "allelic_richness", "seg_sites"),
                     col_pattern = list("ho", 
                                "pi",
                                "he",
                                c("maf", "major", "minor", "maj.count", "min.count"),
                                c("pa", "^g$"),
                                c("chi_association", "p_armitage_association", "log_odds_ratio_association", "se_association", "associated_allele", "gmmat_association"),
                                c("pHWE"),
                                c("RF_importance", "RF_importance_pvals"),
                                c("gp_effect"),
                                c("fis", "var_comp_b", "var_comp_c", "nk"),
                                c("richness", "^g$"),
                                c("seg_sites", "g_prob_seg", "prob_seg_var", "prob_seg")),
                     col_types = list("numeric",
                                   "numeric",
                                   "numeric",
                                   c("numeric", "character", "character", "integer", "integer"),
                                   c("numeric", "integer"),
                                   c("numeric", "numeric", "numeric", "numeric", "character", "numeric"),
                                   "numeric",
                                   c("numeric", "numeric"),
                                   "numeric",
                                   c("numeric", "numeric", "numeric", "integer"),
                                   c("numeric", "integer"),
                                   c("numeric", "integer", "numeric", "numeric")))
window.stats <- list(stat = "tajimas_d",
                  col_pattern = list(c("ws.theta", "ts.theta", "num_seg", "D", "n_snps")),
                  col_types = list(c("numeric", "numeric", "integer", "numeric", "integer")))
pairwise.stats <- list(stat = c("fst", "abba_baba"),
                    col_pattern = list(c("fst","var_comp_a", "var_comp_b", "var_comp_c", "zfst", "fst_id", "nk"),
                                       c("D_abba_baba", "abba", "baba", "nk")),
                    col_types = list(c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "integer"),
                                     c("numeric", "numeric", "numeric", "integer")))
sample.stats <- list(stat = c("het_hom_ratio", "hs"),
                  col_pattern = list("Het/Hom",
                                     "hs"),
                  col_types = list("numeric",
                                   "numeric"))
pop.stats <- list(stat = "ne",
               col_pattern = list(c("LDNe_", "Neb_", "He_", "Ne", "CI")),
               col_types = list(c("numeric", "numeric", "numeric", "numeric", "numeric")))
other.stats <- list(stat = c("ld", "genetic_distances", "isolation_by_distance", "geographic_distance", "prop_poly", "allele_frequency_matrix"),
                    types = list(c("LD", "single.window"),
                                 c("genetic_distances"),
                                 c("ibd"),
                                 c("geo_dist"),
                                 c("weighted.means"),
                                 c("allele_frequency_matrix")),
                    col_patttern = list(c("CLD", "rsq", "^pval$", "Dprime"), NA, NA, NA, "prop_poly", NA))


#=============build lists for each stat====================
statistic_index <- vector("list", length(c(single.stats$stat, window.stats$stat, pairwise.stats$stat, sample.stats$stat, pop.stats$stat, other.stats$stat)))

tracker <- 1
for(i in 1:length(single.stats[[1]])){
  names(statistic_index)[tracker] <- single.stats[[1]][i]
  statistic_index[[tracker]] <- list(category = "single", types = c("single", "single.window", "bootstraps", "weighted.means"),
                                     col_pattern = single.stats[[2]][[i]],
                                     col_types = single.stats[[3]][[i]])
  tracker <- tracker + 1
}

for(i in 1:length(window.stats[[1]])){
  names(statistic_index)[tracker] <- window.stats[[1]][i]
  statistic_index[[tracker]] <- list(category = "window", types = c("single.window", "bootstraps", "weighted.means"),
                                     col_pattern = window.stats[[2]][[i]],
                                     col_types = window.stats[[3]][[i]])
  tracker <- tracker + 1
}

for(i in 1:length(pairwise.stats[[1]])){
  names(statistic_index)[tracker] <- pairwise.stats[[1]][i]
  if(pairwise.stats[[1]][i] == "fst"){
    statistic_index[[tracker]] <- list(category = "pairwise", types = c("pairwise", "pairwise.window", "weighted.means", "fst.matrix", "bootstraps"),
                                       col_pattern = pairwise.stats[[2]][[i]],
                                       col_types = pairwise.stats[[3]][[i]])
  }
  else{
    statistic_index[[tracker]] <- list(category = "pairwise", types = c("pairwise", "pairwise.window", "weighted.means"),
                                       col_pattern = pairwise.stats[[2]][[i]],
                                       col_types = pairwise.stats[[3]][[i]])
  }
  
  tracker <- tracker + 1
}

for(i in 1:length(sample.stats[[1]])){
  names(statistic_index)[tracker] <- sample.stats[[1]][i]
  statistic_index[[tracker]] <- list(category = "sample", types = c("sample", "weighted.means"),
                                     col_pattern = sample.stats[[2]][[i]],
                                     col_types = sample.stats[[3]][[i]])
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


#=============stats storage colnames========
pos_col_names <- vector("list", 5)
names(pos_col_names) <- c("stats", "pairwise.stats", "window.stats", "pairwise.window.stats", "sample.stats")
pos_col_types <- pos_col_names
for(i in 1:length(statistic_index)){
  if("single" %in% statistic_index[[i]]$category){
    pos_col_names$stats <- c(pos_col_names$stats, statistic_index[[i]]$col_pattern)
    pos_col_types$stats <- c(pos_col_types$stats, statistic_index[[i]]$col_types)
    if("single.window" %in% statistic_index[[i]]$types){
      pos_col_names$window.stats <- c(pos_col_names$window.stats, statistic_index[[i]]$col_pattern)
      pos_col_types$window.stats <- c(pos_col_types$window.stats, statistic_index[[i]]$col_types)
      
    }
  }
  
  if("pairwise" %in% statistic_index[[i]]$category){
    pos_col_names$pairwise.stats <- c(pos_col_names$pairwise.stats, statistic_index[[i]]$col_pattern)
    pos_col_types$pairwise.stats <- c(pos_col_types$pairwise.stats, statistic_index[[i]]$col_types)
    
    if("pairwise.window" %in% statistic_index[[i]]$types){
      pos_col_names$pairwise.window.stats <- c(pos_col_names$pairwise.window.stats, statistic_index[[i]]$col_pattern)
      pos_col_types$pairwise.window.stats  <- c(pos_col_types$pairwise.window.stats , statistic_index[[i]]$col_types)
    }
  }
  
  if("sample" %in% statistic_index[[i]]$category){
    pos_col_names$sample.stats <- c(pos_col_names$sample.stats, statistic_index[[i]]$col_pattern)
    pos_col_types$sample.stats <- c(pos_col_types$sample.stats, statistic_index[[i]]$col_types)
  }
}

for(i in 1:length(pos_col_names)){
  dups <- which(duplicated(pos_col_names[[i]]))
  if(any(dups)){
    pos_col_names[[i]] <- pos_col_names[[i]][-dups]
    pos_col_types[[i]] <- pos_col_types[[i]][-dups]
  }
}

pos_col_names <- lapply(pos_col_names, function(x) gsub("\\^","", gsub("\\$", "", x))) # fix the exact colname symbols

#=============test_snps=====================
set.seed(1212)
test_snps <- subset_snpR_data(stickSNPs, 1:12, sample(nsamps(stickSNPs), 10, F))
test_snps <- filter_snps(test_snps)

sample.meta(test_snps)$pop <- rep(c("ASP", "PAL"), 5)
sample.meta(test_snps)$fam <- rep(c("A", "B"), each = 5)



#============================save==========================
.internal.data <- list(test_snps = test_snps, statistic_index = statistic_index, pos_col_names = pos_col_names, pos_col_types = pos_col_types)
usethis::use_data(.internal.data, internal = T, overwrite = T)
