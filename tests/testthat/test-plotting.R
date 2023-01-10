#===================plot_structure================
test_that("structure",{
  skip_on_cran(); skip_on_ci()
  
  str_path <- "C://usr/bin/structure.exe"
  clumpp_path <- "C://usr/bin/CLUMPP.exe"
  skip_if(!file.exists(str_path))
  skip_if(!file.exists(clumpp_path))
  
  p <- plot_structure(stickSNPs[1:10, pop = c("ASP", "PAL")], "pop", k = 2:3, method = "structure", structure_path = str_path, clumpp = FALSE)
  
  expect_true(ggplot2::is.ggplot(p$plot))
  expect_equal(as.character(unique(p$plot_data$K)), c("K = 2", "K = 3"))
  expect_true(!any(is.na(p$plot_data)))
  expect_true(max(p$plot_data$Percentage) <= 1)
  expect_equal(as.character(unique(p$plot_data$pop)), c("ASP", "PAL"))
  # not internally calced, just a check for proper prep and parsing. Note that the K plot details were all checked against structure harvester
  
  expect_error(plot_structure(stickSNPs[1:10, pop = c("ASP", "PAL")], "pop", k = 2:3, method = "structure", structure_path = str_path, clumpp = FALSE, iterations = 1),
               regexp = "one or fewer iterations")
  
  
  # stripping axis text
  p_cat <- plot_structure(stickSNPs[1:10, pop = c("ASP", "PAL")], "pop", k = 2:3, 
                          method = "structure", structure_path = str_path, clumpp = FALSE, strip_col_names = "p$")
  expect_equal(p_cat$plot$labels$x, "po")
  
  
  # clumpp
  p2 <-  plot_structure(stickSNPs[1:10, pop = c("ASP", "PAL")], "pop", k = 2:4, reps = 2, method = "structure", 
                        structure_path = str_path, clumpp = TRUE, clumpp_path = clumpp_path)
  expect_true(ggplot2::is.ggplot(p2$plot))
  expect_true(all(c("r_1", "r_2", "clumpp") %in% names(p2$data$K_2)))
  
  # evanno is there and OK?
  expect_true(names(p2$K_plot) == c("raw", "evanno"))
  if(names(p2$K_plot) == c("raw", "evano")){
    expect_true(colnames(p2$K_plot$evanno) == c("K", "mean_est_ln_prob", "lnpK", "lnppK", "deltaK", "sd_est_ln_prob"))
    expect_identical(round(p2$K_plot$evanno$deltaK, 4), c(NA, round(3.889087, 4), NA))
  }
})


test_that("snmf",{
  skip_on_cran();
  skip_if_not_installed("LEA")

  .make_it_quiet(p <- plot_structure(stickSNPs[pop = c("ASP", "PAL")], "pop", k = 2:3, clumpp = FALSE))
  
  expect_true(ggplot2::is.ggplot(p$plot))
  expect_equal(as.character(unique(p$plot_data$K)), c("K = 2", "K = 3"))
  expect_true(!any(is.na(p$plot_data)))
  expect_true(max(p$plot_data$Percentage) <= 1)
  expect_equal(as.character(unique(p$plot_data$pop)), c("ASP", "PAL"))
  # not internally calced, just a check for proper prep and parsing. Note that the K plot details were all checked against structure harvester
})

test_that("snapclust",{
  skip_on_cran();
  skip_if_not_installed("adegenet")

  
  expect_warning(p <- plot_structure(stickSNPs[pop = c("ASP", "PAL")], "pop", k = 2:3, method = "snapclust", clumpp = FALSE), "adegenet maintainers do not")
  
  expect_true(ggplot2::is.ggplot(p$plot))
  expect_equal(as.character(unique(p$plot_data$K)), c("K = 2", "K = 3"))
  expect_true(!any(is.na(p$plot_data)))
  expect_true(max(p$plot_data$Percentage) <= 1)
  expect_equal(as.character(unique(p$plot_data$pop)), c("ASP", "PAL"))
  # not internally calced, just a check for proper prep and parsing. Note that the K plot details were all checked against structure harvester
})

#===================plot_structure_map===================
# test_that("structure map",{
#   skip_on_cran(); skip_on_ci()
#   skip_if_not_installed(c("LEA", "ggrepel", "sf", "ggsn", "scatterpie", "maps", "sf"))
# 
#   lat_long <- data.frame(SMR = c(44.365931, -121.140420), CLF = c(44.267718, -121.255805), OPL = c(44.485958, -121.298360), ASP = c(43.891693, -121.448360), UPD = c(43.891755, -121.451600), PAL = c(43.714114, -121.272797)) # coords for point
#   lat_long <- t(lat_long)
#   colnames(lat_long) <- c("lat", "long")
#   lat_long <- as.data.frame(lat_long)
#   lat_long$pop <- rownames(lat_long)
#   psf <- sf::st_as_sf(as.data.frame(lat_long), coords = c("long", "lat"))
#   psf <- sf::`st_crs<-`(psf, "EPSG:4326")
# 
#   # get the assignments
#   .make_it_quiet(assignments <- plot_structure(stickSNPs, "pop", alpha = 10, k = 3, clumpp = FALSE)) # get structure-like results
# 
#   # get a map of oregon as a background from the maps package. Note that this map is a bit odd as an sf, but works as an example.
#   background <- maps::map("state", "oregon", plot = FALSE)
#   background <- sf::st_as_sf(background)
# 
#   p2 <- plot_structure_map(assignments, k = 3, facet = "pop", pop_coordinates = psf, sf = list(background), radius_scale = .2, scale_bar = list(dist = 40, dist_unit = "km", transform = T), compass = list(symbol = 16, scale = 0.2), ask = FALSE)
#   expect_true(ggplot2::is.ggplot(p2))
# })

#==================plot_clusters=====================

test_that("pca",{
  local_edition(3)
  skip_on_cran();
  
  set.seed(1212)
  .make_it_quiet(p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop"))
  expect_true(ggplot2::is.ggplot(p$plots$pca))
  expect_snapshot_value(p$data$pca[,c("PC1", "PC2")], style = "serialize") # run entirely via R's prcomp function, shouldn't change with a set seed.
})

test_that("tsne",{
  local_edition(3)
  skip_on_cran();
  
  skip_if_not_installed(c("Rtsne", "mmtsne"))
  set.seed(1212)
  .make_it_quiet(p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop", plot_type = "tsne"))
  
  expect_true(ggplot2::is.ggplot(p$plots$tsne))
  expect_snapshot_value(p$data$tsne[,c("PC1", "PC2")], style = "serialize") # run entirely via R's prcomp function, shouldn't change with a set seed.
})

test_that("umap",{
  local_edition(3)
  skip_on_cran();
  skip_if_not_installed(c("umap"))
  set.seed(1212)
  .make_it_quiet(p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop", plot_type = "umap"))
  
  expect_true(ggplot2::is.ggplot(p$plots$umap))
  expect_snapshot_value(p$data$umap[,c("PC1", "PC2")], style = "serialize") # run entirely via R's prcomp function, shouldn't change with a set seed.
})

#==================plot_manhattan==========
test_that("manhattan plots", {
  
  # with snpRdata
  x <- stickSNPs
  sample.meta(x)$phenotype <- sample(c("case", "control"), nsamps(stickSNPs), TRUE)
  x <- calc_association(x, response = "phenotype", method = "armitage")
  p <- plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
                      log.p = TRUE)
  expect_true(ggplot2::is.ggplot(p$plot))
  
  
  # with data.frame
  y <- get.snpR.stats(x, stats = "association")
  p2 <- plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
                       log.p = TRUE)
  expect_true(ggplot2::is.ggplot(p2$plot))
  expect_equivalent(p$data, p2$data)
  
  
  
  
  # with a rug
  rug_data <- data.frame(chr = c("groupX", "groupVIII"), start = c(0, 1000000),
                         end = c(5000000, 6000000), gene = c("A", "B"))

  # point style
  p3 <- plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
                       log.p = TRUE, rug_data = rug_data)
  expect_true(is(p3$plot$layers[[2]]$geom, "GeomRug"))

  # ribbon style
  p4 <- plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
                       log.p = TRUE, rug_data = rug_data, rug_style = "ribbon")
  expect_true(is(p4$plot$layers[[2]]$geom, "GeomSegment"))
  
  # with rug labels
  ## point
  p3 <- plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
                       log.p = TRUE, rug_data = rug_data, rug_label = "gene")
  expect_true(all(c("label", "position") %in% names(p3$plot$layers[[2]]$mapping)))
  ## ribbon
  p4 <- plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
                       log.p = TRUE, rug_data = rug_data, rug_style = "ribbon", rug_label = "gene")
  expect_true(all(c("label", "start_position", "end_position") %in% names(p4$plot$layers[[2]]$mapping)))
  
  
  # with tajimas D
  skip_on_cran(); # slower
  x <- calc_tajimas_d(x, facets = "pop.chr", sigma = 100, step = 50)
  expect_true(ggplot2::is.ggplot(plot_manhattan(x, "D", TRUE, "pop.chr")$p))
})

#=================qq=====================
test_that("qq plots",{
  # with snpRdata 
  asso <- stickSNPs
  sample.meta(asso)$phenotype <- sample(c("case", "control"), nsamps(stickSNPs), TRUE)
  asso <- calc_association(asso, response = "phenotype")

  p <- plot_qq(asso, "gmmat_pval_phenotype")
  ggplot2::is.ggplot(p$.base)

  # multiple facets
  x <- stickSNPs
  sample.meta(x)$phenotype <- sample(c("case", "control"), nsamps(stickSNPs), TRUE)
  x <- calc_association(x, c("pop.fam", "pop", ".base"), "phenotype",
                        method = "armitage")
  p <- plot_qq(x, "p_armitage_phenotype", c("pop.fam", "pop", ".base"))
  expect_true(all(unlist(lapply(p, ggplot2::is.ggplot))))
  out <- lapply(p, ggplot2::ggplot_build)
  expect_equal(length(levels(out$fam.pop$data[[1]]$PANEL)), 24)
  expect_equal(length(levels(out$pop$data[[1]]$PANEL)), 6)
  expect_equal(length(levels(out$.base$data[[1]]$PANEL)), 1)
  
})

#=================LD======================
test_that("LD heatmap", {
  ld <- calc_pairwise_ld(stickSNPs, "pop.chr", subfacets = list(pop = c("ASP", "PAL"), chr = c("groupXIX", "groupIV")))
  p <- plot_pairwise_ld_heatmap(ld, "pop.chr")
  
  expect_true(ggplot2::is.ggplot(p$plot))
  expect_equal(unique(p$plot$data$snp.subfacet), c("groupXIX", "groupIV"))
  expect_equal(unique(p$plot$data$var), c("ASP", "PAL"))
  
  p2 <- plot_pairwise_ld_heatmap(ld, "pop.chr", snp.subfacet = "groupIV", sample.subfacet = "ASP")
  expect_true(ggplot2::is.ggplot(p2$plot))
  expect_equal(unique(p2$plot$data$snp.subfacet), c("groupIV"))
  expect_equal(unique(p2$plot$data$var), c("ASP"))
})

#==============fst======================
test_that("FST heatmap",{
  fst <- calc_pairwise_fst(stickSNPs, "pop")
  p <- plot_pairwise_fst_heatmap(fst, "pop")
  
  expect_true(ggplot2::is.ggplot(p))
  expect_equal(unique(p$data$facet), "pop")
  expect_true("label" %in% names(p$labels))
  
  p2 <- plot_pairwise_fst_heatmap(fst, "pop", print_fst = FALSE)
  expect_true(!"label" %in% names(p2$labels))
})


#============tree====================
test_that("tree generation",{
  skip_if_not_installed("ape")
  if("ggtree" %in% installed.packages()){
    skip_if_not(utils::packageVersion("ggtree") >= numeric_version("3.1.2"))
  }
  
  
  .make_it_quiet(tree <- calc_tree(stickSNPs, update_bib = FALSE))
  
  expect_true(is(tree$.base$.base, "phylo"))
  expect_true(all(colnames(stickSNPs) %in% tree$.base$.base$tip.label))
  
  .make_it_quiet(tree <- calc_tree(stickSNPs, facets = "pop", update_bib = FALSE))
  expect_true(is(tree$pop$.base, "phylo"))
  expect_true(all(unique(sample.meta(stickSNPs)$pop) %in%  tree$pop$.base$tip.label))
  
  .make_it_quiet(tree <- calc_tree(stickSNPs, "pop.chr", update_bib = FALSE))
  expect_true(all(unlist(lapply(unlist(tree, recursive = F), function(x) is(x, "phylo")))))
  expect_equal(names(tree$chr.pop), unique(snp.meta(stickSNPs)$chr))
  expect_true(all(unique(sample.meta(stickSNPs)$pop) %in% tree$chr.pop$groupV$tip.label))
  
  .make_it_quiet(tree <- calc_tree(stickSNPs, "pop", boot = 3, update_bib = FALSE))
  expect_true(is(tree$pop$.base, "phylo"))
  vals <- as.numeric(gsub("%", "", tree$pop$.base$node.label))
  expect_true(all(vals <= 100 & vals >= 0 & !is.na(vals)))
  
})

#============sfs========
test_that("sfs plot",{
  # 1D
  .make_it_quiet(sfs <- plot_sfs(stickSNPs, projection = 10))
  expect_true(ggplot2::is.ggplot(sfs))
  expect_true(sum(sfs$data$N, na.rm = T) <= nsnps(stickSNPs))
  expect_equal(unique(sfs$data$p1), 0)
  expect_true(max(sfs$data$p2[which(!is.na(sfs$data$N))]) <= 100/2) # folded

  # 2D
  .make_it_quiet(sfs2 <- plot_sfs(stickSNPs, facet = "pop", pops =  c("ASP", "UPD"), projection = c(10, 10)))
  expect_true(ggplot2::is.ggplot(sfs2))
  expect_true(sum(sfs2$data$N, na.rm = T) <= nsnps(stickSNPs))
  expect_true(max(sfs2$data$p2[which(!is.na(sfs2$data$N))] + sfs2$data$p1[which(!is.na(sfs2$data$N))]) <= 10) # folded
  
  
  # works with sfs provided
  .make_it_quiet(sfs_provided <- calc_sfs(stickSNPs, projection = 10))
  sfs_provided_plot <- plot_sfs(sfs_provided)
  expect_identical(sfs$data, sfs_provided_plot$data)
  
  .make_it_quiet(sfs_provided2 <- calc_sfs(stickSNPs, facet = "pop", pops =  c("ASP", "UPD"), projection = c(10, 10)))
  sfs_provided2_plot <- plot_sfs(sfs_provided2)
  expect_identical(sfs2$data, sfs_provided2_plot$data)
  
  # sanity checks
  attr(sfs_provided, "pop") <- NULL
  expect_error(plot_sfs(sfs_provided), regexp = "1D SFS must have a pop attribute with length 1")
  attr(sfs_provided, "pop") <- c("check1", "check2")
  expect_error(plot_sfs(sfs_provided), regexp = "1D SFS pop attribute must be length 1")
  
  attr(sfs_provided2, "pop") <- NULL
  expect_error(plot_sfs(sfs_provided2), regexp = "2D SFS must have a pop attribute with length 2")
  attr(sfs_provided2, "pop") <- c("check1")
  expect_error(plot_sfs(sfs_provided2), regexp = "2D SFS pop attribute must be length 2")
  
  expect_error(plot_sfs(sfs2), "x must be either a snpRdata object, a matrix containing a 2D SFS, or a numeric vector containing a 1D sfs")
  
  attr(sfs_provided, "pop") <- "check"
  p1 <- plot_sfs(sfs_provided)
  attr(sfs_provided, "pops") <- "check"
  attr(sfs_provided, "pop") <- NULL
  p2 <- plot_sfs(sfs_provided)
  expect_equal(p1, p2)
  
  expect_warning(.make_it_quiet(sfs <- plot_sfs(stickSNPs, projection = 50, fold = FALSE)), "Without ancestral and derived character states, unfolded spectra will be misleading")
  
})

#=======diagnostic=======
test_that("diagnostic plots",{
  set.seed(1234122)
  # correct plots and axis names
  expect_warning(dp <- plot_diagnostic(.internal.data$test_snps), "Few remaining SNPs after filtering")
  expect_equal(names(dp), c("fis", "sfs", "maf", "pca", "missingness"))
  expect_equal(c(dp$fis$labels$x,
                 dp$sfs$labels$x,
                 dp$maf$labels$x,
                 dp$pca$labels$x,
                 dp$missingness$labels$x),
               c("fis", "Minor Allele Count", "Minor Allele Frequency", "PC1 (33.06%)", "Individual"))
  expect_equal(c(dp$fis$labels$y,
                 dp$sfs$labels$y,
                 dp$maf$labels$y,
                 dp$pca$labels$y,
                 dp$missingness$labels$y),
               c("density", "log10(N)", "density", "PC2 (18.76%)", "Proportion of loci with missing data"))
  expect_false("colour" %in% names(dp$pca$labels))
  expect_false("colour" %in% names(dp$missingness$labels))
  
  # if unfolded sfs
  expect_warning(dp2 <- plot_diagnostic(.internal.data$test_snps, fold_sfs = FALSE), "Few remaining SNPs after filtering")
  expect_equal(dp2$sfs$labels$x, "Derived Allele Count")
  
  # colored by pop
  expect_warning(dp3 <- plot_diagnostic(.internal.data$test_snps, facet = "pop"), "Few remaining SNPs after filtering")
  expect_true("colour" %in% names(dp3$pca$labels))
  expect_true("colour" %in% names(dp3$missingness$labels))
})

