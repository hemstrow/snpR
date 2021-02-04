context("facet creation")



test_that("single sample facet",{
  tdat <- add.facets.snpR.data(stickSNPs, "pop")
  tdwm <- cbind(tdat@facet.meta, tdat@geno.tables$wm)
  tdwm$n_tot <- rowSums(tdat@geno.tables$wm)
  
  expect_equal(unique(tdwm[,c("subfacet", "n_tot")])$n_tot, 
               c(nrow(tdat@sample.meta), as.numeric(table(tdat@sample.meta$pop)))) # correct number of samples in each batch? Also checks order in facet meta vs ac
  expect_equal(tdat@ac$n_total, rowSums(tdat@geno.tables$as)) # everything lines up between ac and geno tables
})

test_that("complex sample facet",{
  tdat <- add.facets.snpR.data(stickSNPs, "pop.fam")
  tdwm <- cbind(tdat@facet.meta, tdat@geno.tables$wm)
  tdwm$n_tot <- rowSums(tdat@geno.tables$wm)
  
  ncomp <-  as.data.frame(table(tdat@sample.meta[,c("fam", "pop")]))
  ncomp <- ncomp[order(ncomp$fam),]$Freq
  expect_equal(unique(tdwm[,c("subfacet", "n_tot")])$n_tot, 
               c(nrow(tdat@sample.meta), ncomp)) # correct number of samples in each batch? Also checks order in facet meta vs ac
  
  expect_equal(tdat@ac$n_total, rowSums(tdat@geno.tables$as)) # everything lines up between ac and geno tables
})



# snp facets not actually added, so not checked here.

