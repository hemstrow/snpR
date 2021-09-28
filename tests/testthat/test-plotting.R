context("plots")

#===================plot_structure================
test_that("structure",{
  skip_on_cran()

  
  str_path <- "C://usr/bin/structure.exe"
  skip_if(!file.exists(str_path))
  
  p <- plot_structure(stickSNPs[pop = c("ASP", "PAL")], "pop", k = 2:3, method = "structure", structure_path = str_path, clumpp = FALSE)
  
  expect_true(ggplot2::is.ggplot(p$plot))
  expect_equal(as.character(unique(p$plot_data$K)), c("K = 2", "K = 3"))
  expect_true(!any(is.na(p$plot_data)))
  expect_true(max(p$plot_data$Percentage) <= 1)
  expect_equal(as.character(unique(p$plot_data$pop)), c("ASP", "PAL"))
  # not internally calced, just a check for proper prep and parsing. Note that the K plot details were all checked against structure harvester
  
  expect_error(plot_structure(stickSNPs[pop = c("ASP", "PAL")], "pop", k = 2:3, method = "structure", structure_path = str_path, clumpp = FALSE, iterations = 1),
               regexp = "one or fewer iterations")
})


test_that("snmf",{
  skip_if_not_installed("LEA")
  
  p <- plot_structure(stickSNPs[pop = c("ASP", "PAL")], "pop", k = 2:3)
  
  expect_true(ggplot2::is.ggplot(p$plot))
  expect_equal(as.character(unique(p$plot_data$K)), c("K = 2", "K = 3"))
  expect_true(!any(is.na(p$plot_data)))
  expect_true(max(p$plot_data$Percentage) <= 1)
  expect_equal(as.character(unique(p$plot_data$pop)), c("ASP", "PAL"))
  # not internally calced, just a check for proper prep and parsing. Note that the K plot details were all checked against structure harvester
})

test_that("snapclust",{
  skip_if_not_installed("adegenet")
  
  p <- plot_structure(stickSNPs[pop = c("ASP", "PAL")], "pop", k = 2:3, method = "snapclust")
  
  expect_true(ggplot2::is.ggplot(p$plot))
  expect_equal(as.character(unique(p$plot_data$K)), c("K = 2", "K = 3"))
  expect_true(!any(is.na(p$plot_data)))
  expect_true(max(p$plot_data$Percentage) <= 1)
  expect_equal(as.character(unique(p$plot_data$pop)), c("ASP", "PAL"))
  # not internally calced, just a check for proper prep and parsing. Note that the K plot details were all checked against structure harvester
})

#===================plot_structure_map===================
test_that("structure map",{
  skip_if_not_installed(c("LEA", "ggrepel", "sf", "ggsn", "scatterpie", "maps"))
  
  lat_long <- data.frame(SMR = c(44.365931, -121.140420), CLF = c(44.267718, -121.255805), OPL = c(44.485958, -121.298360), ASP = c(43.891693, -121.448360), UPD = c(43.891755, -121.451600), PAL = c(43.714114, -121.272797)) # coords for point
  lat_long <- t(lat_long)
  colnames(lat_long) <- c("lat", "long")
  lat_long <- as.data.frame(lat_long)
  lat_long$pop <- rownames(lat_long)
  psf <- sf::st_as_sf(as.data.frame(lat_long), coords = c("long", "lat"))
  psf <- sf::`st_crs<-`(psf, "EPSG:4326")

  # get the assignments
  assignments <- plot_structure(stickSNPs, "pop", alpha = 10, k = 3) # get structure-like results

  # get a map of oregon as a background from the maps package. Note that this map is a bit odd as an sf, but works as an example.
  background <- maps::map("state", "oregon")
  background <- sf::st_as_sf(background)
  
  p2 <- plot_structure_map(assignments, k = 3, facet = "pop", pop_coordinates = psf, sf = list(background), radius_scale = .2, scale_bar = list(dist = 40, dist_unit = "km", transform = T), compass = list(symbol = 16, scale = 0.2))
  expect_true(ggplot2::is.ggplot(p2))
})

#==================plot_clusters=====================
test_that("pca",{
  local_edition(3)
  skip_on_cran()
  
  set.seed(1212)
  p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop")
  expect_true(ggplot2::is.ggplot(p$plots$pca))
  expect_snapshot(p$data$pca) # run entirely via R's prcomp function, shouldn't change with a set seed.
})

test_that("tsne",{
  local_edition(3)
  skip_on_cran()
  
  skip_if_not_installed(c("Rtsne", "mmtsne"))
  set.seed(1212)
  p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop", plot_type = "tsne")
  
  expect_true(ggplot2::is.ggplot(p$plots$tsne))
  expect_snapshot(p$data$tsne) # run entirely via R's prcomp function, shouldn't change with a set seed.
})

test_that("umap",{
  local_edition(3)
  skip_on_cran()
  skip_if_not_installed(c("umap"))
  set.seed(1212)
  p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop", plot_type = "umap")
  
  expect_true(ggplot2::is.ggplot(p$plots$umap))
  expect_snapshot(p$data$umap) # run entirely via R's prcomp function, shouldn't change with a set seed.
})

#==================plot_manhattan==========
test_that("manhattan plots", {
  skip_if_not_installed(c("GMMAT", "AGHmatrix"))
  asso <- stickSNPs[pop = "ASP"]
  set.seed(1212)
  sample.meta(asso)$cat_phenotype <- sample(c("A", "B"), ncol(asso), replace = TRUE)
  asso <- calc_association(asso, response = "cat_phenotype")
  
  p <- plot_manhattan(asso, "gmmat_pval_cat_phenotype", chr = "group", log.p = TRUE)
  expect_true(ggplot2::is.ggplot(p$plot))
})

#=================qq=====================
# test_that("qq plots",{
#   set.seed(1212)
#   sample.meta(asso)$cat_phenotype <- sample(c("A", "B"), ncol(asso), replace = TRUE)
#   asso <- calc_association(asso, response = "cat_phenotype")
#   
#   p <- plot_qq(asso, )
#   
#   
# })

#=================LD======================
test_that("LD heatmap", {
  ld <- calc_pairwise_ld(stickSNPs, "pop.group", subfacets = list(pop = c("ASP", "PAL"), group = c("groupXIX", "groupIV")))
  p <- plot_pairwise_LD_heatmap(ld, "pop.group")
  
  expect_true(ggplot2::is.ggplot(p$plot))
  expect_equal(unique(p$plot$data$snp.subfacet), c("groupXIX", "groupIV"))
  expect_equal(unique(p$plot$data$var), c("ASP", "PAL"))
  
  p2 <- plot_pairwise_LD_heatmap(ld, "pop.group", snp.subfacet = "groupIV", sample.subfacet = "ASP")
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
test_that("tree plot",{
  skip_if_not_installed(c("ggtree", "ape"))
  if("ggtree" %in% installed.packages()){
    skip_if_not(utils::packageVersion("ggtree") >= numeric_version("3.1.2"))
  }
  
  
  tree <- plot_tree(stickSNPs)
  
  expect_true(ggplot2::is.ggplot(tree$.base$.base$plot))
  expect_true(all(colnames(stickSNPs) %in% tree$.base$.base$plot$data$label))
  
  tree <- plot_tree(stickSNPs, facets = "pop")
  expect_true(ggplot2::is.ggplot(tree$pop$.base$plot))
  expect_true(all(unique(sample.meta(stickSNPs)$pop) %in% tree$pop$.base$plot$data$label))
  
  tree <- plot_tree(stickSNPs, "pop.group")
  expect_true(all(unlist(lapply(purrr::map(tree$group.pop, "plot"), ggplot2::is.ggplot))))
  expect_equal(names(tree$group.pop), unique(snp.meta(stickSNPs)$group))
  expect_true(all(unique(sample.meta(stickSNPs)$pop) %in% tree$group.pop$groupV$plot$data$label))
  
  tree <- plot_tree(stickSNPs, "pop", boot = 3)
  expect_true(ggplot2::is.ggplot(tree$pop$.base$plot))
  vals <- as.numeric(gsub("%", "", tree$pop$.base$plot$data$label[-which(tree$pop$.base$plot$data$label %in% sample.meta(stickSNPs)$pop)]))
  expect_true(all(vals <= 100 & vals >= 0 & !is.na(vals)))
  
})

#============sfs========
test_that("sfs plot",{
  # 1D
  sfs <- plot_sfs(stickSNPs, projection = 100)
  expect_true(ggplot2::is.ggplot(sfs))
  expect_true(sum(sfs$data$N, na.rm = T) <= nsnps(stickSNPs))
  expect_equal(unique(sfs$data$p1), 0)
  expect_true(max(sfs$data$p2[which(!is.na(sfs$data$N))]) <= 100/2) # folded

  # 2D
  sfs2 <- plot_sfs(stickSNPs, facet = "pop", pops =  c("ASP", "UPD"), projection = c(50, 50))
  expect_true(ggplot2::is.ggplot(sfs2))
  expect_true(sum(sfs2$data$N, na.rm = T) <= nsnps(stickSNPs))
  expect_true(max(sfs2$data$p2[which(!is.na(sfs2$data$N))] + sfs2$data$p1[which(!is.na(sfs2$data$N))]) <= 50) # folded
  
  
  # works with sfs provided
  sfsp <- calc_sfs(stickSNPs, projection = 100)
  sfsp <- plot_sfs(sfs = sfsp)
  expect_identical(sfs$data, sfsp$data)
  
  sfsp2 <- calc_sfs(stickSNPs, facet = "pop", pops =  c("ASP", "UPD"), projection = c(50, 50))
  sfsp2 <- plot_sfs(sfs = sfsp2)
  expect_identical(sfs2$data, sfsp2$data)
})

