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
  set.seed(1212)
  p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop")
  expect_true(ggplot2::is.ggplot(p$plots$pca))
  expect_snapshot(p$data$pca) # run entirely via R's prcomp function, shouldn't change with a set seed.
})

test_that("tsne"){
  skip_if_not_installed(c("Rtsne", "mmtsne"))
  set.seed(1212)
  p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop", plot_type = "tsne")
  
  expect_true(ggplot2::is.ggplot(p$plots$tsne))
  expect_snapshot(p$data$tsne) # run entirely via R's prcomp function, shouldn't change with a set seed.
}

test_that("tsne"){
  skip_if_not_installed(c("umap"))
  set.seed(1212)
  p <- plot_clusters(stickSNPs[pop = c("ASP", "PAL")], "pop", plot_type = "umap")
  
  expect_true(ggplot2::is.ggplot(p$plots$umap))
  expect_snapshot(p$data$umap) # run entirely via R's prcomp function, shouldn't change with a set seed.
}

#==================plot_manhattan==========
test_that("manhattan plots"){
  skip_if_not_installed(c("GMMAT", "AGHmatrix"))
  asso <- stickSNPs[pop = "ASP"]
  set.seed(1212)
  sample.meta(asso)$cat_phenotype <- sample(c("A", "B"), ncol(asso), replace = TRUE)
  asso <- calc_association(asso, response = "cat_phenotype")
  
  p <- plot_manhattan(asso, "gmmat_pval_cat_phenotype", chr = "group", log.p = TRUE)
  expect_true(ggplot2::is.ggplot(p$plot))
}
