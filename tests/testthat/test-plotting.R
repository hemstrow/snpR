#===================plot_structure================
test_that("structure",{
  skip_on_cran(); skip_on_ci()
  
  str_path <- "C://usr/bin/structure.exe"
  clumpp_path <- "C://usr/bin/CLUMPP.exe"
  skip_if(!file.exists(str_path))
  skip_if(!file.exists(clumpp_path))
  
  p <- plot_structure(stickSNPs[1:10, pop = c("ASP", "PAL")], "pop", k = 2:3, method = "structure", structure_path = str_path, clumpp = FALSE)
  
  expect_true(ggplot2::is_ggplot(p$plot))
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
  expect_true(ggplot2::is_ggplot(p2$plot))
  expect_true(all(c("r_1", "r_2", "clumpp") %in% names(p2$data$K_2)))
  
  # evanno is there and OK?
  expect_true(all(names(p2$K_plot) == c("raw", "evanno")))
  if(all(names(p2$K_plot) == c("raw", "evano"))){
    expect_true(all(colnames(p2$K_plot$evanno) == c("K", "mean_est_ln_prob", "lnpK", "lnppK", "deltaK", "sd_est_ln_prob")))
    expect_identical(round(p2$K_plot$evanno$deltaK, 4), c(NA, round(3.889087, 4), NA))
  }
})


test_that("snmf",{
  skip_on_cran();
  skip_if_not_installed("LEA")

  .make_it_quiet(p <- plot_structure(stickSNPs[pop = c("ASP", "PAL")], "pop", k = 2:3, clumpp = FALSE))
  
  expect_true(ggplot2::is_ggplot(p$plot))
  expect_equal(as.character(unique(p$plot_data$K)), c("K = 2", "K = 3"))
  expect_true(!any(is.na(p$plot_data)))
  expect_true(max(p$plot_data$Percentage) <= 1)
  expect_equal(as.character(unique(p$plot_data$pop)), c("ASP", "PAL"))
  # not internally calced, just a check for proper prep and parsing. Note that the K plot details were all checked against structure harvester
  
  # check that it works OK with ".base"
  
  .make_it_quiet(p <- plot_structure(stickSNPs[pop = c("ASP", "PAL")], ".base", k = 2:3, clumpp = FALSE))
})

test_that("snapclust",{
  skip_on_cran();
  skip_if_not_installed("adegenet")

  
  expect_warning(p <- plot_structure(stickSNPs[pop = c("ASP", "PAL")], "pop", k = 2:3, method = "snapclust", clumpp = FALSE), "adegenet maintainers do not")
  
  expect_true(ggplot2::is_ggplot(p$plot))
  expect_equal(as.character(unique(p$plot_data$K)), c("K = 2", "K = 3"))
  expect_true(!any(is.na(p$plot_data)))
  expect_true(max(p$plot_data$Percentage) <= 1)
  expect_equal(as.character(unique(p$plot_data$pop)), c("ASP", "PAL"))
  # not internally calced, just a check for proper prep and parsing. Note that the K plot details were all checked against structure harvester
})

#===================plot_structure_map===================
# test_that("structure map",{
#   skip_on_cran(); skip_on_ci()
#   skip_if_not_installed(c("LEA", "ggrepel", "ggspatial", "scatterpie", "sf"))
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
#   background <- rnaturalearth::ne_states(iso_a2 = "US", returnclass = "sp")
#   background <- sf::st_as_sf(background)
#   background <- background[background$name %in% "Oregon",]
#   
#   p1 <- plot_structure_map(assignments, k = 3, facet = "pop", 
#                            pop_coordinates = psf,
#                            radius_scale = .2,
#                            ask = FALSE)
#   expect_true(ggplot2::is_ggplot(p1))
#   expect_equal(length(p1$layers), 4)
# 
#   p2 <- plot_structure_map(assignments, k = 3, facet = "pop", 
#                            pop_coordinates = psf, 
#                            layers = list(ggplot2::geom_sf(data = background)),
#                            radius_scale = .2,
#                            ask = FALSE)
#   expect_true(ggplot2::is_ggplot(p2))
#   expect_true(length(p2$layers), 5) # does it have the additional background layer?
# })

#==================plot_clusters=====================

test_that("pca",{
  local_edition(3)
  skip_on_cran();
  
  set.seed(1212)
  .make_it_quiet(p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop", "pca", simplify_output = TRUE))
  expect_true(ggplot2::is_ggplot(p$pca))
  # expect_snapshot_value(p$data$pca[,c("PC1", "PC2")], style = "serialize") # run entirely via R's prcomp function, shouldn't change with a set seed.
  
  # check that loadings are returned if requested
  .make_it_quiet(p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop", "pca", simplify_output = FALSE))
  expect_true("pca_loadings" %in% names(p))
})

test_that("tsne",{
  local_edition(3)
  skip_on_cran();
  
  skip_if_not_installed(c("Rtsne", "mmtsne"))
  set.seed(1212)
  .make_it_quiet(p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop", plot_type = "tsne", simplify_output = TRUE))
  
  expect_true(ggplot2::is_ggplot(p$tsne))
  # expect_snapshot_value(p$data$tsne[,c("PC1", "PC2")], style = "serialize") # run entirely via R's prcomp function, shouldn't change with a set seed.
})

test_that("umap",{
  local_edition(3)
  skip_on_cran();
  skip_if_not_installed(c("umap"))
  set.seed(1212)
  .make_it_quiet(p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop", plot_type = "umap", simplify_output = TRUE))
  
  expect_true(ggplot2::is_ggplot(p$umap))
  # expect_snapshot_value(p$data$umap[,c("PC1", "PC2")], style = "serialize") # run entirely via R's prcomp function, shouldn't change with a set seed.
})

test_that("dapc",{
  skip_on_cran();
  skip_if_not_installed("adegenet")
  
  set.seed(1212)
  expect_error(p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop", "dapc", dapc_clustering_max_n_clust = NULL, simplify_output = TRUE), "please supply all dapc")
  expect_error(p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop", "dapc", 
                                  dapc_clustering_max_n_clust = NULL, dapc_clustering_nclust = 5, simplify_output = TRUE), 
               "supply both dapc_clustering_npca and dapc_clustering_ndisc")
  expect_error(p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop", "dapc",
                                  dapc_clustering_max_n_clust = NULL,
                                  dapc_npca = 5, dapc_ndisc = NULL, simplify_output = TRUE), 
               "pply both dapc_npca and dapc_ndisc")
  
  .make_it_quiet(p <- plot_clusters(stickSNPs, "pop", "dapc", 
                                    dapc_clustering_max_n_clust = NULL, 
                                    dapc_clustering_npca = 20, 
                                    dapc_clustering_nclust = 5, 
                                    dapc_npca = 20, 
                                    dapc_ndisc = 4, simplify_output = TRUE))
  expect_true(ggplot2::is_ggplot(p$dapc))
  
  # case with one discriminant
  .make_it_quiet(p <- plot_clusters(stickSNPs, "pop", "dapc", 
                                    dapc_clustering_max_n_clust = NULL, 
                                    dapc_clustering_npca = 20, 
                                    dapc_clustering_nclust = 2, 
                                    dapc_npca = 20, 
                                    dapc_ndisc = 1, simplify_output = TRUE))
  expect_true(ggplot2::is_ggplot(p$dapc))
  mapping <- capture_output(p$dapc$mapping, print = TRUE)
  expect_true(grepl("`x` -> `\\.cluster`", mapping)) # check that it xluster is mapped to x
  
})

#==================plot_manhattan==========
test_that("manhattan plots", {
  
  # with snpRdata
  x <- stickSNPs
  sample.meta(x)$phenotype <- sample(c("case", "control"), nsamps(stickSNPs), TRUE)
  x <- calc_association(x, response = "phenotype", method = "armitage")
  p <- plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
                      log.p = TRUE, simplify_output = TRUE)
  expect_true(ggplot2::is_ggplot(p))
  
  
  # with data.frame
  y <- get.snpR.stats(x, stats = "association")
  p2 <- plot_manhattan(y$single, "p_armitage_phenotype", chr = "chr",
                       log.p = TRUE, simplify_output = TRUE)
  ggplot2::ggplot_build(p2)$layout$panel_params[[1]]$x$get_labels()
  expect_true(ggplot2::is_ggplot(p2))
  expect_equivalent(ggplot2::ggplot_build(p2)$layout$panel_params[[1]]$x$get_labels(),
                    ggplot2::ggplot_build(p)$layout$panel_params[[1]]$x$get_labels())
  expect_equivalent(ggplot2::ggplot_build(p2)$data[[1]]$y,
                    ggplot2::ggplot_build(p)$data[[1]]$y)
  
  ## simple facets
  x <- calc_ho(x, "pop")
  y <- get.snpR.stats(x, "pop", "ho")
  p1 <- plot_manhattan(y$single, chr = "chr", bp = "position", plot_var = "ho", facets = "subfacet", simplify_output = TRUE)
  p2 <- plot_manhattan(x, chr = "chr", bp = "position", plot_var = "ho", facets = "pop", simplify_output = TRUE)
  expect_equivalent(ggplot2::ggplot_build(p2)$layout$layout$subfacet,
                    ggplot2::ggplot_build(p1)$layout$layout$subfacet)
  
  ## multilevel facets
  x <- calc_ho(x, "pop.fam")
  y <- get.snpR.stats(x, "pop.fam", "ho")
  p1 <- plot_manhattan(y$single, chr = "chr", bp = "position", plot_var = "ho", facets = "subfacet", simplify_output = TRUE)
  p2 <- plot_manhattan(x, chr = "chr", bp = "position", plot_var = "ho", facets = "pop.fam", simplify_output = TRUE)
  expect_equivalent(ggplot2::ggplot_build(p2)$layout$panel_params[[1]]$x$get_labels(),
                    ggplot2::ggplot_build(p1)$layout$panel_params[[1]]$x$get_labels())
  expect_equivalent(ggplot2::ggplot_build(p2)$data[[1]]$y,
                    ggplot2::ggplot_build(p1)$data[[1]]$y)
  expect_equivalent(ggplot2::ggplot_build(p2)$layout$layout$subfacet,
                    ggplot2::ggplot_build(p1)$layout$layout$subfacet)
  y$single$fam <- gsub("\\..+", "", y$single$subfacet)
  y$single$pop <- gsub(".+\\.", "", y$single$subfacet)
  
  ### snpR style facet call
  p3 <- plot_manhattan(y$single, chr = "chr", bp = "position", plot_var = "ho", facets = "fam.pop", simplify_output = TRUE)
  expect_equivalent(ggplot2::ggplot_build(p2)$layout$layout$subfacet,
                    ggplot2::ggplot_build(p3)$layout$layout$subfacet)
  
  ### not snpR style facet call
  p4 <- plot_manhattan(y$single, chr = "chr", bp = "position", plot_var = "ho", facets = c("fam", "pop"), simplify_output = TRUE)
  expect_equivalent(ggplot2::ggplot_build(p3)$layout$layout$subfacet,
                    ggplot2::ggplot_build(p4)$layout$layout$subfacet)
  
  ### errors
  expect_error(plot_manhattan(y$single, chr = "chr", bp = "position", plot_var = "ho", facets = c("fam", "hi")),
               "Some facets not found.+hi.+")
  expect_error(plot_manhattan(y$single, chr = "chr", bp = "position", plot_var = "ho", facets = "fam.hi"),
               "Some facets not found.+hi.+")
  expect_error(plot_manhattan(y$single, chr = "chr", bp = "position", plot_var = "ho", facets = "hi"),
               "Facet not found in")
  
  
  
  # with lambda correction
  pl <- plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
                       log.p = TRUE, lambda_gc_correction = TRUE, simplify_output = TRUE)
  expect_true(".lam" %in% colnames(pl$data))
  if(".lam" %in% colnames(pl$data)){
    expect_true(length(unique(pl$data$.lam)) == 1)
  }
  
  
  
  # with facets
  x <- calc_association(x, response = "phenotype", method = "armitage", facets = "pop")
  p <- plot_manhattan(x, "p_armitage_phenotype", chr = "chr", facets = "pop",
                      log.p = TRUE, simplify_output = TRUE)
  p1 <- ggplot2::ggplot_build(p)
  expect_true(all(p1$layout$layout$subfacet == c("ASP", "CLF", "OPL", "PAL", "SMR", "UPD")))
  ## also with lambda correction
  pl <- plot_manhattan(x, "p_armitage_phenotype", chr = "chr", facets = "pop",
                       log.p = TRUE, lambda_gc_correction = TRUE, simplify_output = TRUE)
  expect_true(".lam" %in% colnames(pl$data))
  if(".lam" %in% colnames(pl$data)){
    expect_true(length(unique(pl$data$.lam)) == 6)
  }
  
  
  
  
  
  
  
  # with a rug
  rug_data <- data.frame(chr = c("groupX", "groupVIII"), start = c(0, 1000000),
                         end = c(5000000, 6000000), gene = c("A", "B"))

  # point style
  p3 <- plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
                       log.p = TRUE, rug_data = rug_data, simplify_output = TRUE)
  expect_true(is(p3$layers[[2]]$geom, "GeomRug"))

  # ribbon style
  p4 <- plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
                       log.p = TRUE, rug_data = rug_data, rug_style = "ribbon", simplify_output = TRUE)
  expect_true(methods::is(p4$layers[[2]]$geom, "GeomSegment"))
  
  # with rug labels
  ## point
  p3 <- plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
                       log.p = TRUE, rug_data = rug_data, rug_label = "gene", simplify_output = TRUE)
  expect_true(all(c("label", "position") %in% names(p3$layers[[2]]$mapping)))
  ## ribbon
  p4 <- plot_manhattan(x, "p_armitage_phenotype", chr = "chr",
                       log.p = TRUE, rug_data = rug_data, rug_style = "ribbon", rug_label = "gene", simplify_output = TRUE)
  expect_true(all(c("label", "start_position", "end_position") %in% names(p4$layers[[2]]$mapping)))
  
  
  # with tajimas D
  skip_on_cran(); # slower
  x <- calc_tajimas_d(x, facets = "pop.chr", sigma = 100, step = 50)
  expect_true(ggplot2::is_ggplot(plot_manhattan(x, "D", TRUE, "pop.chr", simplify_output = TRUE)))
  
  # mean line
  p <- plot_manhattan(x, "ho", facets = "pop", median_line = "red", simplify_output = TRUE)
  p <- ggplot2::ggplot_build(p)
  expect_true(methods::is(p$plot$layers[[2]]$geom, "GeomHline"))
  x <- calc_ho(x)
  p <- plot_manhattan(x, "ho" ,median_line = "red", simplify_output = TRUE)
  p <- ggplot2::ggplot_build(p)
  expect_true(methods::is(p$plot$layers[[2]]$geom, "GeomHline"))
  
  # color by val
  x <- calc_he(x, "pop")
  p <- plot_manhattan(x, "ho", facets = "pop", color_var = "he", 
                      colors = ggplot2::scale_color_viridis_c(option = "A"),
                      simplify_output = TRUE)
  p <- ggplot2::ggplot_build(p)
  expect_true(ggplot2::as_label(p$plot$mapping$colour) == "colvar")
  
  # vlines
  p <- plot_manhattan(x, "ho", facets = "pop", vlines = "red",
                      simplify_output = TRUE)
  p <- ggplot2::ggplot_build(p)
  expect_true(methods::is(p$plot$layers[[2]]$geom, "GeomVline"))
})

#=================qq=====================
test_that("qq plots",{
  set.seed(12211)
  # with snpRdata 
  asso <- stickSNPs
  sample.meta(asso)$phenotype <- sample(c("case", "control"), nsamps(stickSNPs), TRUE)
  asso <- calc_association(asso, response = "phenotype")

  p <- plot_qq(asso, "gmmat_pval_phenotype")
  expect_true(ggplot2::is_ggplot(p$.base))
  
  # lambda correction
  p <- plot_qq(asso, "gmmat_pval_phenotype", lambda_gc_correction = TRUE)
  expect_true(".lam" %in% colnames(p$.base$data))
  if(".lam" %in% colnames(p$.base$data)){
    expect_true(length(unique(p$.base$data$.lam)) == 1)
  }

  # multiple facets
  x <- stickSNPs
  sample.meta(x)$phenotype <- sample(c("case", "control"), nsamps(stickSNPs), TRUE)
  x <- calc_association(x, c("pop.fam", "pop", ".base"), "phenotype",
                        method = "armitage")
  p <- plot_qq(x, "p_armitage_phenotype", c("pop.fam", "pop", ".base"))
  expect_true(all(unlist(lapply(p, ggplot2::is_ggplot))))
  out <- lapply(p, ggplot2::ggplot_build)
  expect_equal(length(levels(out$fam.pop$data[[1]]$PANEL)), 24)
  expect_equal(length(levels(out$pop$data[[1]]$PANEL)), 6)
  expect_equal(length(levels(out$.base$data[[1]]$PANEL)), 1)
  
  # lambda correction
  p <- plot_qq(x, "p_armitage_phenotype", c("pop.fam", "pop", ".base"), lambda_gc_correction = TRUE)
  good <- all(unlist(lapply(p, function(i) ".lam" %in% colnames(i$data))))
  expect_true(good)
  if(good){
    expect_equivalent(unlist(lapply(p, function(i) {
      return(length(unique(i$data$.lam)) > 1)
    })), c(TRUE, TRUE, FALSE))
    
  }
  
  
})

#=================LD======================
test_that("LD heatmap", {
  ld <- calc_pairwise_ld(stickSNPs, "pop.chr", subfacets = list(pop = c("ASP", "PAL"), chr = c("groupXIX", "groupIV")))
  p <- plot_pairwise_ld_heatmap(ld, "pop.chr", simplify_output = TRUE)
  
  expect_true(ggplot2::is_ggplot(p))
  expect_equal(unique(p$data$snp.subfacet), c("groupXIX", "groupIV"))
  expect_equal(unique(p$data$var), c("ASP", "PAL"))
  
  p2 <- plot_pairwise_ld_heatmap(ld, "pop.chr", snp.subfacet = "groupIV", sample.subfacet = "ASP", simplify_output = TRUE)
  expect_true(ggplot2::is_ggplot(p2))
  expect_equal(unique(p2$data$snp.subfacet), c("groupIV"))
  expect_equal(unique(p2$data$var), c("ASP"))
})

#==============fst======================
test_that("FST heatmap",{
  fst <- calc_pairwise_fst(stickSNPs, "pop")
  p <- plot_pairwise_fst_heatmap(fst, "pop")
  
  expect_true(ggplot2::is_ggplot(p))
  expect_equal(unique(p$data$facet), "pop")
  expect_true("label" %in% names(p$labels))
  
  p2 <- plot_pairwise_fst_heatmap(fst, "pop", print_fst = FALSE)
  expect_true(!"label" %in% names(p2$labels))
  
  # lab ordering
  fst <- calc_pairwise_fst(stickSNPs, c("pop", "fam"))
  p3 <- plot_pairwise_fst_heatmap(fst, c("pop", "fam"), 
                                  list(pop = c("PAL", "ASP", "UPD", "CLF", "SMR", "OPL"),
                                       fam = c("A", "B", "C", "D")),print_fst = FALSE)
  p3 <- lapply(p3, ggplot2::ggplot_build)
  expect_equal(p3$pop$layout$panel_params[[1]]$y$get_labels(), c("ASP", "UPD", "CLF", "SMR", "OPL"))
  expect_equal(p3$pop$layout$panel_params[[1]]$x$get_labels(), c("PAL", "ASP", "UPD", "CLF", "SMR"))
  expect_equal(p3$fam$layout$panel_params[[1]]$y$get_labels(), c("B", "C", "D"))
  expect_equal(p3$fam$layout$panel_params[[1]]$x$get_labels(), c("A", "B", "C"))
  
  # print below
  fst <- calc_pairwise_fst(stickSNPs, c("pop", "fam"))
  p3 <- plot_pairwise_fst_heatmap(fst, c("pop", "fam"), 
                                  list(pop = c("PAL", "ASP", "UPD", "CLF", "SMR", "OPL"),
                                       fam = c("A", "B", "C", "D")), lab_lower = TRUE)
  
  p3 <- lapply(p3, ggplot2::ggplot_build)
  expect_equal(p3$pop$layout$panel_params[[1]]$y$get_labels(), c("PAL", "ASP", "UPD", "CLF", "SMR", "OPL"))
  expect_equal(p3$pop$layout$panel_params[[1]]$x$get_labels(), c("PAL", "ASP", "UPD", "CLF", "SMR", "OPL"))
  expect_equal(p3$fam$layout$panel_params[[1]]$y$get_labels(), c("A", "B", "C", "D"))
  expect_equal(p3$fam$layout$panel_params[[1]]$x$get_labels(), c("A", "B", "C", "D"))
  
  # errs
  expect_error(.suppress_specific_warning(plot_pairwise_fst_heatmap(fst, c("pop", "fam"), 
                                                                    list(pop = c("PAL", "ASP", "UPD", "CLF", "SMR", "OPL"),
                                                                         fam = c("A", "B", "D")),print_fst = TRUE, lab_lower = TRUE), "longer object length"),
               "Subfacets in provided facet.order do not exactly match all of those in the provided data for facet: fam.")
  
  expect_error(plot_pairwise_fst_heatmap(fst, c("pop", "fam"), 
                            c("PAL", "ASP", "UPD", "CLF", "SMR", "OPL"), print_fst = TRUE, lab_lower = TRUE),
               "If more than one facet is requested and a facet.order is provided, an order for each facet must be included using a named list, see documentation.")
  
  expect_error(plot_pairwise_fst_heatmap(fst, c("pop"), 
                                         list(c("PAL", "ASP", "UPD", "CLF", "SMR", "OPL"), c("A", "B", "C", "D")), print_fst = TRUE, lab_lower = TRUE),
               "If more than one facet is requested and a facet.order is provided, an order for each facet must be included using a named list with no extra elements, see documentation.\n")
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
  expect_true(ggplot2::is_ggplot(sfs))
  expect_true(sum(sfs$data$N, na.rm = T) <= nsnps(stickSNPs))
  expect_equal(unique(sfs$data$p1), 0)
  expect_true(max(sfs$data$p2[which(!is.na(sfs$data$N))]) <= 100/2) # folded

  # 2D
  .make_it_quiet(sfs2 <- plot_sfs(stickSNPs, facet = "pop", pops =  c("ASP", "UPD"), projection = c(10, 10)))
  expect_true(ggplot2::is_ggplot(sfs2))
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
  expect_warning(dp <- plot_diagnostic(.internal.data$test_snps,
                                       plots = c("fis", "sfs", "maf", "pca", "missingness", "heho")), "Few remaining SNPs after filtering")
  expect_equal(names(dp), c("fis", "sfs", "maf", "pca", "missingness", "heho"))
  expect_equal(c(dp$fis$labels$x,
                 dp$sfs$labels$x,
                 dp$maf$labels$x,
                 dp$pca$labels$x,
                 dp$missingness$labels$x,
                 dp$heho$labels$x),
               c("fis", "Minor Allele Count", "Minor Allele Frequency", "PC1 (36.64%)", "Individual", "Expected Heterozygosity"))
  expect_equal(c(dp$fis$labels$y,
                 dp$sfs$labels$y,
                 dp$maf$labels$y,
                 dp$pca$labels$y,
                 dp$missingness$labels$y,
                 dp$heho$labels$y),
               c("density", "log10(N)", "density", "PC2 (16.26%)", "Proportion of loci with missing data", "Observed Heterozygosity"))
  expect_false("colour" %in% names(dp$pca$labels))
  expect_false("colour" %in% names(dp$missingness$labels))
  
  # if unfolded sfs
  expect_warning(dp2 <- plot_diagnostic(.internal.data$test_snps, fold_sfs = FALSE, plots = "sfs"), "Without ancestral and derived")
  expect_equal(dp2$sfs$labels$x, "Derived Allele Count")
  
  # colored/faceted by pop
  expect_warning(dp3 <- plot_diagnostic(.internal.data$test_snps, facet = "pop", 
                                        plots = c("pca", "missingness", "heho")), "Few remaining SNPs after filtering")
  expect_true("colour" %in% names(dp3$pca$labels))
  expect_true("colour" %in% names(dp3$missingness$labels))
  p1 <- ggplot2::ggplot_build(dp3$heho)
  expect_identical(p1$layout$layout$subfacet, c("ASP", "PAL"))
  
  # maf doesn't error if run alone
  expect_no_error(plot_diagnostic(.internal.data$test_snps, facet = "pop", plots = "maf"))
})


#====roh========
test_that("plot_roh",{
  source_file <-  paste0(tempfile(), ".vcf")
  utils::download.file("https://raw.githubusercontent.com/hemstrow/snpR/refs/heads/dev/R_dev/roh_test.vcf", 
                       destfile = source_file, quiet = TRUE)
  expect_warning(.make_it_quiet(d <- read_vcf(source_file)), "Some levels are duplicated")
  file.remove(source_file)
  expect_warning(
    sample.meta(d) <- cbind(pop = sample(LETTERS[1:4], ncol(d), replace = TRUE),
                            sample.meta(d)), 
    "Some levels are duplicated")
  
  # generate ROH data
  d1 <- calc_roh(d, c("CHROM.pop", "CHROM"), verbose = FALSE)
  
  res <- get.snpR.stats(d1, facets = c("CHROM", "CHROM.pop"), "roh")
  
  # tests--generic errors
  expect_error(plot_roh(d1, facet = c("pop", "sampID")), "Only one")
  expect_error(plot_roh(d1, facet = c("sampID")), "not calculated")
  expect_error(plot_roh(d1, chr = "chr"), "not calculated")
  
  # with snpRobj
  ## with a sampID column
  p1 <- plot_roh(d1, chr = "CHROM")
  expect_identical(ggplot2::ggplot_build(p1)$layout$panel_params[[1]]$y$get_labels(),
                   unique(sample.meta(d1)$sampID))
  d2 <- d1
  ## without
  sample.meta(d2)$sampID <- NULL
  p2 <- plot_roh(d2, chr = "CHROM")
  expect_identical(ggplot2::ggplot_build(p2)$layout$panel_params[[1]]$y$get_labels(),
                   unique(sample.meta(d1)$.sample.id))
  ## with facet
  p1 <- plot_roh(d1, facet = "pop", chr = "CHROM")
  expect_identical(ggplot2::ggplot_build(p1)$layout$panel_params[[1]]$y$get_labels(),
                   sort(unique(sample.meta(d1)$pop)))
  ## with facet, no sampID column
  p2 <- plot_roh(d2, facet = "pop", chr = "CHROM")
  expect_identical(ggplot2::ggplot_build(p2)$layout$panel_params[[1]]$y$get_labels(),
                   sort(unique(sample.meta(d2)$pop)))
  
  
  # with a data.frame
  ## errors
  rroh <- res$roh
  expect_error(p1 <- plot_roh(rroh, chr = "CHROM"), "with column names")
  
  chr_lens <- tapply(snp.meta(d1)$position, snp.meta(d1)$CHROM, max)
  rroh$chr_len <- chr_lens[match(rroh$snp.subfacet, names(chr_lens))] - 100000 # check end error
  expect_error(p1 <- plot_roh(rroh, chr = "snp.subfacet"),"have end positions greater")
  
  ## with no facets
  rroh$chr_len <- chr_lens[match(rroh$snp.subfacet, names(chr_lens))]
  p1 <- plot_roh(rroh, chr = "snp.subfacet")
  expect_identical(ggplot2::ggplot_build(p1)$layout$panel_params[[1]]$y$get_labels(),
                   unique(rroh$sampID))
  ## with facet
  p2 <- plot_roh(rroh, facet = "pop", chr = "snp.subfacet")
  expect_identical(ggplot2::ggplot_build(p2)$layout$panel_params[[1]]$y$get_labels(),
                   sort(unique(rroh$pop)))
  
  # othering chrs
  p1 <- plot_roh(rroh, chr = "snp.subfacet", lab_as_other = c(unique(rroh$snp.subfacet)))
  expect_identical(ggplot2::ggplot_build(p1)$layout$panel_params[[1]]$x$get_labels(),
                   "other")
  p1 <- plot_roh(rroh, chr = "snp.subfacet", lab_as_other = c(unique(rroh$snp.subfacet)[2]))
  expect_identical(ggplot2::ggplot_build(p1)$layout$panel_params[[1]]$x$get_labels(),
                   c(unique(rroh$snp.subfacet)[1], "other"))
})
