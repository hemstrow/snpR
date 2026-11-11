# Plot Runs of Homozygosity

Generate a plot displaying the locations of ROHs across the genome.
Takes either a
[`snpRdata`](https://hemstrow.github.io/snpR/reference/snpRdata-class.md)
object for which ROH information has been calculated or a `data.frame`
with ROH info.

## Usage

``` r
plot_roh(x, facet = NULL, chr = "chr", bp = "position", lab_as_other = NULL)
```

## Arguments

- x:

  snpRdata object or an object coercible to a data.frame. If a snpRdata
  object, needs to have had ROH info calculated with
  [`calc_roh`](https://hemstrow.github.io/snpR/reference/calc_roh.md).
  If present in the sample metadata or in the provided data.table, the
  column "sampID" will be taken as sample IDs for plotting. data.frames
  should have columns:

  - start_position: start positions for ROHs

  - end_position: end positions for ROHs

  - sampID: sample IDs

  - chr_len: total chromosome lengths for the chromosomes on which ROHs
    were detected.

  - chr: a column containing chromosome info *named to match the `chr`
    argument!* Note that for data.frame inputs, chromosomes or
    individuals with no ROHs will not be plotted.

- facet:

  character. *sample* facet by which to color individuals. Note that
  only one facet is plotable at a time.

- chr:

  character, default "chr". Column in either snp metadata or x (for
  snpRdata or data.frame objects, respectively) which defines the
  "chromosome" by which SNP positions will be concatenated along the
  x-axis.

- bp:

  character, default "bp". Column in either snp metadata which defines
  the position in bp of each SNP. Ignored for `data.table` inputs.

- lab_as_other:

  character vector, default NULL. Vector of chromosome/scaffold names
  which will be plotted grouped together, labelled as "other". Useful
  for small scaffolds/etc.

## Value

A `ggplot2` plot showing the positions of ROHs (white areas) against
colored chromosomal backdrops. Note that the data used to plot this can
be returned by calling `$data` on the resulting object.

## Author

William Hemstrom
