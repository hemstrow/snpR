# Create a heatmap from pairwise linkage data.

Prepares a ggplot2 heatmap from pairwise linkage disequilibrium data
stored in a snpRdata object.

## Usage

``` r
plot_pairwise_ld_heatmap(
  x,
  facets = NULL,
  snp.subfacet = NULL,
  sample.subfacet = NULL,
  LD_measure = "CLD",
  r = NULL,
  l.text = "CLD",
  viridis.option = "inferno",
  gradient_colors = NULL,
  title = NULL,
  t.sizes = c(16, 13, 10, 12, 10),
  background = "white",
  simplify_output = FALSE
)
```

## Arguments

- x:

  snpRdata object.

- facets:

  character, default NULL. Categorical metadata variables by which to
  break up plots. Must match facets for which LD data has been
  previously calculated, and only a single facet can be plotted at once.
  See
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md)
  for more details.

- snp.subfacet:

  character, default NULL. Specific snp-specific levels of the provided
  facet to plot. See examples.

- sample.subfacet:

  character, default NULL. Specific sample-specific levels of the
  provided facet to plot. See examples.

- LD_measure:

  character, default rsq. LD metric to plot. Must be present in the
  calculated LD data.

- r:

  Numeric. Region of the chromosome to subset and plot. Given in kb in
  the format numeric vector c(lower, upper).

- l.text:

  character, default "CLD". Legend title.

- viridis.option:

  character, default "inferno". Viridis color scale option to use. Other
  color scales may be substituted by appending the
  scale_color_continuous and scale_fill_continuous ggplot functions to
  the produced plot using the '+' operator. See `scale_gradient` for
  details.

- gradient_colors:

  character vector of length 2 or NULL, default NULL. If provided, the
  low and high LD colors for the color scale. Provided as an alternative
  to the viridis scales. Take care to set a `background` color not on
  this scale. For example, `c("white", "black")` is an excellent choice,
  but will make missing data inaccurately appear appear as if it has a
  very low LD unless a different `background` color is selected.

- title:

  character. Plot title.

- t.sizes:

  numeric, default c(16, 13, 10, 12, 10). Text sizes, given as c(title,
  legend.title, legend.ticks, axis, axis.ticks).

- background:

  character, default "white". Background color for plot.

- simplify_output:

  If TRUE, only the ggplot object will be return. This is optimal, since
  the data is already returned in that object, but is not the default
  due to backwards consistency with old code.

## Value

A list containing:

- plot: A pairwise LD heatmap as a ggplot object.

- dat: Data used to generate the ggplot object. If `simplify_output` is
  `FALSE`, only the ggplot object is returned.

## Details

Since the output is a ggplot object, options can be added or changed by
adding "+ function()" to the end of `LD_pairwise_heatmap`. Some common
options are also built into this function as arguments, but can be
overwritten freely.

Specific facets and facet levels can be requested as long as LD data for
those facets has been previously calculated. Facets should be provided
as described in
[`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).
Only one facet can be plotted at once. Plots for individual SNP or
sample-specific facet levels for the specified facet can also be
requested. See examples. For facets with many categories, this is
strongly recommended, since plots with many facet levels become
computationally challenging to produce and difficult to interpret. A
facets argument of NULL will plot the base level facet, a facet argument
of "all" is not accepted.

Note that NA LD values will be colored white in the resulting plot.

## Author

William Hemstrom

Nicholas Sard

## Examples

``` r
if (FALSE) { # \dontrun{
# get LD data
dat <- calc_pairwise_ld(stickSNPs, c("pop.chr"))

# produce plots for linkage group IX in the ASP and CLF populations.
plot_pairwise_ld_heatmap(dat, c("pop.chr"), "groupIX", c("ASP", "CLF"))

# produce plots for every population for linkage group IV
plot_pairwise_ld_heatmap(dat, c("pop.chr"), "groupIV")
} # }
```
