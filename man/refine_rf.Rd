% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/association_functions.R
\name{refine_rf}
\alias{refine_rf}
\title{Refine a random forest model via sequential removal of uninformative SNPs.}
\usage{
refine_rf(
  rf,
  response,
  facets = NULL,
  subfacet = NULL,
  formula = NULL,
  num.trees = 10000,
  trim = 0.5,
  trim_cuttoffs = NULL,
  importance = "impurity_corrected",
  interpolate = "bernoulli",
  par = FALSE,
  ...
)
}
\arguments{
\item{rf}{The random forest model to be refined. List containg snpRdata
object (named $data) and a \code{\link[ranger]{ranger}} model, named
$models$x$model, where x is the facet/subfacet. Identical to ojects created
by \code{\link{run_random_forest}}.}

\item{response}{character. Name of the column containing the response
variable of interest. Must match a column name in sample metadata. Response
must be categorical, with only two categories.}

\item{facets}{Character, default NULL. Facet to run. Only a single facet and
facet level (subfacet) may be provided at once, and must match a facet and
level in the provided input model. If NULL, runs the base level facet.}

\item{subfacet}{Character, default NULL. Facet level (subfacet) to run. Only
a single facet and facet level (subfacet) may be provided at once, and must
match a facet and level in the provided input model. If NULL, runs the base
level facet.}

\item{formula}{charcter, default NULL. Model for the response variable, as
described in \code{\link[stats]{formula}}. If NULL, the model will be
equivalent to response ~ 1.}

\item{num.trees}{numeric, default 10000. Number of trees to grow. Higher
numbers will increase model accuracy, but increase calculation time. See
\code{\link[ranger]{ranger}} for details.}

\item{trim}{numeric, default 0.5. Percentages of SNPs to be trimmed between
model iterations. Multiple trim levels can be provided corresponding to
different trim cuttoffs. If less trim levels are provided than needed to
describe every trim_cuttoff interval, will trim a single SNP below the
final cuttoff. See details.}

\item{trim_cuttoffs}{numeric, default NULL. Specifies the number of SNPs
below which to change trim percentages. If NULL, the default, trims at the
given level untill 1 SNP remails. See details.}

\item{importance}{character, default "impurity_corrected". The method by
which SNP importance is determined. Options: \itemize{\item{impurity}
\item{impurity_corrected} \item{permutation}}. See
\code{\link[ranger]{ranger}} for details.}

\item{interpolate}{character, default "bernoulli". Interpolation method for
missing data. Options: \itemize{\item{bernoulli: }binomial draws for the
minor allele. \item{af: } insertion of the average allele frequency}.}

\item{par}{numeric, default FALSE. Number of parallel computing cores to use
for computing RFs across multiple facet levels or within a single facet if
only a single category is run (either a one-category facet or no facet).}

\item{...}{Additional arguments passed to \code{\link[ranger]{ranger}}.}
}
\value{
A list containing: \itemize{\item{error_delta: } A data.frame noting
  the number of SNPs and corresponding prediction_error in each model
  iteration. \item{confusion_matrices: } An array containing confusion
  matrices for categorical responses. The third subscript denotes model
  iteration ('[,,1]' would reference model 1.) \item{best_model: } The model
  with the lowest prediction error from the provided dataset, in the format
  provided by \code{\link{run_random_forest}}}
}
\description{
Improves the prediction accuracy of a random forest model via iterative
removal of uninformative SNPs. In each step, the SNPs with the lowest
absolute value importance are removed from the model. Depending on the
provided arguments, multiple trim percentages can be provided for different
SNP number cuttoffs.
}
\details{
Random Forest models can fail to predict well with "noisy" data, where most
explanitory variables are uniformative. Since most whole-genome sequence data
is like this, it can be useful to "trim" a data to remove uniformative snps.
Since even models constructed on "noisy" data tend to pick out the most
important SNPs with a decent degree of accuracy, an initial random forest
model can be used to select SNPs to remove. Sequential removal of unimportant
SNPs in this way can radically improve prediction accuracy, although SNP
p-values stop being informative.

Multiple trim levels can be specified, which will determine the percentage of
SNPs removed according to different cuttoff levels of remaining SNPs. For
example, providing trim levels of 0.9, 0.5, and 0.1 with cuttoffs of 1000 and
100 will trim 90% of SNPs untill 1000 remain, then trim 50% untill 100
remain, then trim 10% thereafter.

If less trim levels are provided than needed for each cuttoff, a single SNP
will be removed each step below the final cuttoff. For example, providing
trim levels of 0.9 and 0.5 and cuttoffs of 1000 and 100 will trim 90% of SNPs
untill 1000 remain, then 50% untill 100 remain, then one at a time.

If the data contains informative SNPs, prediction accuracy should improve on
average untill informative SNPs begin to be removed, at which point accuracy
will decrease.

mtry will be set to the number of SNPs in each run.

As usual, facets can be requested. However, in this case, only a single facet
and facet level (subfacet) may be provided at once, and must match a facet
and level in the provided input model.
}
\references{
Wright, Marvin N and Ziegler, Andreas. (2017). ranger: A Fast
  Implementation of Random Forests for High Dimensional Data in C++ and R.
  \emph{Journal of Statistical Software}.

Goldstein et al. (2011). Random forests for genetic association
  studies. \emph{Statistical Applications in Genetics and Molecular Biology}.
}
\seealso{
\code{\link[ranger]{ranger}} \code{\link[ranger]{predict.ranger}}.

\code{\link{run_random_forest}}
}
