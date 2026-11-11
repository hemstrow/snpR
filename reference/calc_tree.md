# Generate phylogenetic-like clustering trees.

Generate dendritic trees in the style of a phylogenetic tree for
individuals or groups of individuals from snpR data. Note that this
function is not overwrite safe.

## Usage

``` r
calc_tree(
  x,
  facets = NULL,
  distance_method = "Edwards",
  interpolate = "bernoulli",
  tree_method = "nj",
  root = FALSE,
  boot = FALSE,
  boot_par = FALSE,
  update_bib = FALSE
)
```

## Arguments

- x:

  snpRdata object.

- facets:

  character or NULL, default NULL. Facets for which to calculate genetic
  distances, as described in
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).
  If snp or base level facets are requested, distances will be between
  individuals. Otherwise, distances will be between the levels of the
  sample facets.

- distance_method:

  character, default "Edwards". Name of the method to use. Options:

  - Edwards: Angular distance as described in Edwards 1971.

  See
  [`calc_genetic_distances`](https://hemstrow.github.io/snpR/reference/calc_genetic_distances.md).

- interpolate:

  character, default "bernoulli". Missing data interpolation method,
  solely for individual/individual distances. Options detailed in
  documentation for
  [`format_snps`](https://hemstrow.github.io/snpR/reference/format_snps.md).

- tree_method:

  character, default nj. Method by which the tree is constructed from
  genetic distances. Options:

  - nj: Neighbor-joining trees, via
    [`nj`](https://rdrr.io/pkg/ape/man/nj.html).

  - bionj: BIONJ trees, according to Gascuel 1997, via
    [`bionj`](https://rdrr.io/pkg/ape/man/bionj.html).

  - upgma: UPGMA trees, via
    [`hclust`](https://rdrr.io/r/stats/hclust.html).

- root:

  character or FALSE, default FALSE. A vector containing the requested
  roots for each facet. Roots are specified by a string matching either
  the individual sample or sample facet level by which to root. If FALSE
  for a given facet, trees will be unrooted. Note that all UPGMA trees
  are automatically rooted, so this argument is ignored for that tree
  type.

- boot:

  numeric or FALSE, default FALSE. The number of bootstraps to do for
  each facet. See details.

- boot_par:

  numeric or FALSE, default FALSE. If a number, bootstraps will be
  processed in parallel using the supplied number of cores.

- update_bib:

  character or FALSE, default FALSE. If a file path to an existing .bib
  library or to a valid path for a new one, will update or create a .bib
  file including any new citations for methods used. Useful given that
  this function does not return a snpRdata object, so a
  [`citations`](https://hemstrow.github.io/snpR/reference/citations.md)
  cannot be used to fetch references.

## Value

A nested, named list containing plots, trees, and bootstraps for each
facet and facet level.

## Details

Trees are generated via the [`nj`](https://rdrr.io/pkg/ape/man/nj.html)
or [`bionj`](https://rdrr.io/pkg/ape/man/bionj.html) ape package for nj
or bionj trees. Plots of the resulting trees can be produced using other
plotting tools, such as the `geom_tree` function from `ggtree`. These
are not produced automatically because `ggtree` can have unexpected
outputs, but the process is straightforward and examples are provided in
the *examples* section of this documentation. For more information, see
the documentation for those functions and packages.

Bootstraps are conducted by re-sampling SNPs with replacement, according
to Felsenstein (1985). If no snp level facets are provides, loci are
resampled without restraint. If a snp level facet is provided, loci are
only resampled within the levels of that facet (e.g. within
chromosomes).

The genetic distances used to make the trees are calculated using
[`calc_genetic_distances`](https://hemstrow.github.io/snpR/reference/calc_genetic_distances.md).
If a sample facet is used, that function uses code derived from
[`adegenet`](https://rdrr.io/pkg/adegenet/man/adegenet.package.html).
Please cite them and the actual method (e.g. Edwards, A. W. F. (1971))
alongside the tree-building approach.

Bootstrapping is done via the boot.phylo function in the ape package,
and as such does not support parallel runs on Windows machines.

## References

Felsenstein, J. (1985). Confidence Limits on Phylogenies: An Approach
Using the Bootstrap. Evolution, 39(4), 783–791.
https://doi.org/10.2307/2408678

Gascuel, O. (1997).BIONJ: an improved version of the NJ algorithm based
on a simple model of sequence data. Molecular Biology and Evolution,
14(7), 685–695. https://doi.org/10.1093/oxfordjournals.molbev.a025808

Paradis, E., Claude, J. and Strimmer, K. (2004). APE: analyses of
phylogenetics and evolution in R language. Bioinformatics, 20, 289–290.

## Author

William Hemstrom

## Examples

``` r
# Calculate nj trees for the base facet, each chromosome, 
# and for each population.

# note: some versions of ggtree/ggplot2 give an error relating to unable
# to use 'fortify()' when attempting to plot that is potentially due to
# ape masking something when loading snpR. Saving the tree object as
# a .RDS object with 'saveRDS()', restarting R, loading the object in,
# then plotting without loading snpR seems to fix the issue.

if (FALSE) { # \dontrun{

# make the trees
tp <- plot_tree(stickSNPs, c(".base", "pop", "chr"), 
                root = c(FALSE, "PAL", FALSE))

# plot using the ggtree package. Not done internally due to unpredictable
# ggtree behavior
library(ggtree)
ggplot(tp$pop$.base, aes(x, y)) + 
  geom_tree() +
  geom_tiplab() + theme_tree()

# Calculate bionj trees for pop with bootstrapping
tp <- plot_tree(stickSNPs, "pop", root = "PAL", boot = 5)
## plot
ggplot(tp$pop$.base, aes(x, y)) + 
  geom_tree() +
  geom_tiplab() + theme_tree() +
  geom_text2(ggplot2::aes(subset = !isTip, label = label),
                          hjust = -.3)
} # }
```
