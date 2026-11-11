# Generate qq plots from p values

Generate qq (quantile-quantile) plots for p-values, splitting by
provided facets, using either a snpRdata object or a data.frame.
Produces a
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html) object,
which is modifiable as usual.

## Usage

``` r
plot_qq(x, plot_var, facets = NULL, lambda_gc_correction = FALSE)
```

## Arguments

- x:

  snpRdata object or data.frame containing p-values to plot.

- plot_var:

  character, name of the p-values to plot. For a snpRdata object, the
  name should refer to the name of the column produced when fetching
  data via
  [`get.snpR.stats`](https://hemstrow.github.io/snpR/reference/get.snpR.stats.md).
  For a data.frame, it should refer to the column name directly.

- facets:

  character, default NULL. Facets by which to split the plot. See
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).

- lambda_gc_correction:

  Correct for inflated significance due to population and/or family
  structure using the \\\gamma\_{GC}\\ approach described in Price et al
  2010.

## References

Price, A., Zaitlen, N., Reich, D. et al. New approaches to population
stratification in genome-wide association studies. Nat Rev Genet 11,
459–463 (2010). https://doi.org/10.1038/nrg2813

## Author

William Hemstrom

## Examples

``` r
if (FALSE) { # \dontrun{
# from a snpRdata object directly
x <- stickSNPs
sample.meta(x)$phenotype <- sample(c("case", "control"), nsamps(stickSNPs), TRUE)
x <- calc_association(x, c("pop.fam", "pop", ".base"), "phenotype",
                      method = "armitage")
p <- plot_qq(x, "p_armitage_phenotype", c("pop.fam", "pop", ".base"))
                      
                      
# from a data.frame
y <- get.snpR.stats(x, c("pop.fam", "pop", ".base"), "association")
y <- y$single


# with facet/subfacet columns:
p <- plot_qq(y, "p_armitage_phenotype", c("fam.pop", "pop", ".base"))


# from raw data with only one faceting column
z <- y[y$facet == "pop",]
z <- z[,c(1, 5)]
colnames(z)[1] <- "pop"
p <- plot_qq(z, "p_armitage_phenotype", "pop")
} # }
```
