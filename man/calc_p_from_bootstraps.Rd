% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrapping_functions.R
\name{calc_p_from_bootstraps}
\alias{calc_p_from_bootstraps}
\title{Calculate p-values from bootstrapped distributions.}
\usage{
calc_p_from_bootstraps(
  x,
  facets = "all",
  statistics = "all",
  alt = "two-sided",
  par = FALSE,
  fwe_method = "BY",
  fwe_case = c("by_facet", "overall")
)
}
\arguments{
\item{x}{snpRdata object.}

\item{facets}{character, default "all". Facets to use.}

\item{statistics}{character, default "all". Vector naming the statistics to
calculate p-values for. By default calculates p-values for all possible
stats.}

\item{alt}{character, default "two-sided". Specifies the alternative
hypothesis to be used. Options: \itemize{ \item "less": probability that a
bootstrapped value is as small or smaller than observed. \item "greater":
probability that a bootstrapped value is as large or larger than observed.
\item "two-sided": probability that a bootstrapped value is as or more
extreme than observed. }}

\item{par}{numeric or FALSE, default FALSE. If numeric, the number of cores to
use for parallel processing.}

\item{fwe_method}{character, default c("bonferroni", "holm", "BH", "BY"). Type
of Family-Wise Error correction (multiple testing correction) to use. For
details and options, see \code{\link[stats]{p.adjust}}.}

\item{fwe_case}{character, default c("by_facet", "by_subfacet", "overall").
How should Family-Wise Error correction (multiple testing correction) be
applied? \itemize{\item{"by_facet":} Each facet supplied (such as pop or
pop.fam) is treated as a set of tests. \item{"by_subfacet":} Each level of
each subfacet is treated as a separate set of tests. \item{"overall":} All
tests are treated as a set.}}
}
\value{
snpRdata object, with p-values merged into the stats or pairwise.stats
 sockets.
}
\description{
\code{calc_p_from_bootstraps} finds p-values for observed \emph{smoothed
window} statistics from bootstrapped distributions, such as produced by
\code{\link{do_bootstraps}}.
}
\details{
Calculates p-values for smoothed values of a statistic based upon a
bootstrapped null distribution of that statistic using an empirical continuous
distribution function.

p-values can be generated for specific snp or sample metadata categories
using the facets argument, as described in \code{\link{Facets_in_snpR}}. Only
facets for which bootstrap data and raw statistical data have both been
calculated will be run. "all" and NULL follow the typical facet rules.

Likewise, p-values can be generated for specific statistics using the
statistics argument. Only statistics for which bootstrap data and raw
statistical data have both been calculated will be run. By default, all stats
for which a bootstrap null distribution has been generated will be run.
}
\examples{
\dontrun{
# add statistics and generate bootstraps
x <- calc_basic_snp_stats(stickSNPs, c("chr.pop"), sigma = 200, step = 150)
x <- do_bootstraps(x, facets = c("chr.pop"), boots = 1000, sigma = 200, step = 150)
x <- calc_p_from_bootstraps(x)
get.snpR.stats(x, "chr.pop", "single.window") # pi, ho, etc
get.snpR.stats(x, "chr.pop", "pairwise.window") # fst
}
}
\seealso{
ecdf
}
\author{
William Hemstrom
}
