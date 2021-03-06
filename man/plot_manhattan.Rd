% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_manhattan}
\alias{plot_manhattan}
\title{Generate a manhattan plot from snpRdata or a data.frame.}
\usage{
plot_manhattan(
  x,
  plot_var,
  window = FALSE,
  facets = NULL,
  chr = "chr",
  bp = "position",
  snp = NULL,
  chr.subfacet = NULL,
  sample.subfacet = NULL,
  significant = NULL,
  suggestive = NULL,
  highlight = "significant",
  sig_below = FALSE,
  log.p = FALSE,
  abs = FALSE,
  viridis.option = "plasma",
  viridis.hue = c(0.2, 0.5),
  t.sizes = c(16, 12, 10),
  colors = c("black", "slategray3")
)
}
\arguments{
\item{x}{snpRdata or data.frame object containing the data to be plotted.}

\item{plot_var}{character. A character string naming the statistic to be
plotted. For snpRdata, these names correspond to any previously calculated
statistics.}

\item{window}{logical, default FALSE. If TRUE, sliding window averages will
instead be plotted. These averages must have first been calculated with
calc_smoothed_averags. Ignored if x is a data.frame.}

\item{facets}{character or NULL, default NULL. Facets by which to break
plots, as described in \code{\link{Facets_in_snpR}}. For non-window stats,
the any snp.specific facets will be ignored. Ignored if x is a data.frame.}

\item{chr}{character, default "chr". Column in either snp metadata or x (for
snpRdata or data.frame objects, respectively) which defines the
"chromosome" by which SNP positions will be concatenated along the x-axis.
If window = TRUE and a snpRdata object, this will be ignored in favor of
the SNP specific facet provided to the facets argument.}

\item{bp}{character, default "bp". Column in either snp metadata or x (for
snpRdata or data.frame objects, respectively) which defines the position in
bp of each SNP.}

\item{snp}{character, default NULL. Column in either snp metadata or x (for
snpRdata or data.frame objects, respectively) containing snpIDs to use for
highlighting. Ignored if no highlighting is requested.}

\item{chr.subfacet}{character, default NULL. Specific chromosomes to plot.
See examples.}

\item{sample.subfacet}{character, default NULL. Specific sample-specific
levels of the provided facet to plot. If x is a data.frame, this can refer
to levels of a column titled "subfacet". See examples.}

\item{significant}{numeric, default NULL. Value at which a line will be drawn
designating significant SNPs. If highlight = "significant", SNPs above this
level will also be labeled.}

\item{suggestive}{numeric, default NULL. Value at which a line will be drawn
designating suggestive SNPs. If highlight = "suggestive", SNPs above this
level will also be labeled.}

\item{highlight}{character, numeric, or FALSE, default "significant".
Controls SNP highlighting. If either "significant" or "suggestive", SNPs
above those respetive values will be highlighted. If a numeric vector, SNPs
corresponding to vector entries will be highlighted. See details.}

\item{sig_below}{logical, default FALSE. If TRUE, treats values lower than
the significance threshold as significant.}

\item{log.p}{logical, default FALSE. If TRUE, plot variables and thresholds
will be transformed to -log.}

\item{abs}{logical, default FALSE. If TRUE, converts the plot variable to
it's absolute value.}

\item{viridis.option}{character, default "plasma". Viridis color scale option
to use for significance lines and SNP labels. See
\code{\link[ggplot2]{scale_gradient}} for details.}

\item{viridis.hue}{numeric, default c(0.2, 0.5). Two values between 0 and 1
listing the hues at which to start and stop on the viridis palette defined
by the viridis.option argument. Lower numbers are darker.}

\item{t.sizes}{numeric, default c(16, 12, 10). Text sizes, given as
c(strip.title, axis, axis.ticks).}

\item{colors}{character, default c("black", "slategray3"). Colors to
alternate across chromosomes.}
}
\value{
A list containing \itemize{\item{plot: } A ggplot manhattan plot.
  \item{data: } Raw plot data.}
}
\description{
Creates a ggplot-based manhattan plot, where chromosomes/scaffolds/ect are
concatenated along the x-axis. Can optionally highlight requested SNPs or
those that pass an arbitrary significance threshold and facet plots by
defined sample-specific variables such as population.
}
\details{
Unlike most snpR functions, this function works with either a snpRdata object
or a data.frame. For snpRdata objects snp-specific or sliding window
statistics can be plotted. In both cases, the facet argument can be used to
define facets to plot, as described in \code{\link{Facets_in_snpR}}. For
typical stats, name of the snp meta-data column containing
chromosome/scaffold information must be supplied to the "chr" argument. For
windowed stats, chr is instead inferred from the snp-specific facet used to
create the smoothed windows. In both cases, the requested facets must exactly
match those used to calculate statistics! If x is a data frame, the "chr"
argument must also be given, and the "facets" argument will be ignored.

A column defining the position of the SNP within the chromsome must be
provided, and is "position" by default.

Specific snp and chr levels can also be requested using the chr.subfacet and
sample.subfacet arguments. See examples. For data.frames, sample.subfacets
levels must refer to a column in x titled "subfacet".

Specific snps can be highlighted and annotated. If a significance level is
requested, SNPs above this level will be highlighted by default. SNPs above
the suggestive line can also be highlighted by providing "suggestive" to the
highlight argument. Alternatively, individual SNPs can be highlighted by
providing a numeric vector. For snpR data, this will correspond to the SNP's
row in the snpRdata object. For data.frames, it will correspond to a
".snp.id" column if it exists, and the row number if not. The label for
highlighted SNPs will be either chr_bp by default or given in the column
named by the "snp" argument.
}
\examples{
# make some data
x <- calc_basic_snp_stats(stickSNPs, "pop.group", sigma = 200, step = 50)

# plot pi, breaking apart by population, keeping only the groupIX and
# groupIV chromosomes and the ASP, PAL, and SMR populations, with
# significant and suggestive lines plotted and SNPs
# with pi below the significance level labeled.
plot_manhattan(x, "pi", facets = "pop",
chr = "group", chr.subfacet = c("groupIX", "groupIV"),
sample.subfacet = c("ASP", "OPL", "SMR"),
significant = 0.05, suggestive = 0.15, sig_below = TRUE)

# plot FST for the ASP/PAL comparison across all chromosomes,
# labeling the first 10 SNPs in x (by row) with their ID
plot_manhattan(x, "fst", facets = "pop.group",
sample.subfacet = "ASP~PAL", highlight = 1:20,
chr = "group", snp = ".snp.id")

# plot sliding-window FST between ASP and CLF
# and between OPL and SMR
plot_manhattan(x, "fst", window = TRUE, facets = c("pop.group"),
chr = "group", sample.subfacet = c("ASP~CLF", "OPL~SMR"),
significant = .29, suggestive = .2)

# plot using a data.frame,
# using log-transformed p-values
## grab data
y <- get.snpR.stats(x, "pop")
## plot
plot_manhattan(y, "pHWE", facets = "pop", chr = "group",
significant = 0.0001, suggestive = 0.001,
log.p = TRUE, highlight = FALSE)

}
\author{
William Hemstrom
}
