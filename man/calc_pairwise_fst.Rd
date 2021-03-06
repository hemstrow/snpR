% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat_functions.R
\name{calc_pairwise_fst}
\alias{calc_pairwise_fst}
\title{Pairwise FST from SNP data.}
\usage{
calc_pairwise_fst(x, facets, method = "WC")
}
\arguments{
\item{x}{snpRdata. Input SNP data.}

\item{facets}{character. Categorical metadata variables by which to break up
analysis. See \code{\link{Facets_in_snpR}} for more details.}

\item{method}{character, default "WC". Defines the FST estimator to use.
Options: \itemize{ \item{WC: } Wier and Cockerham (1984). \item{Wier: } Wier
(1990) \item{Hohenlohe: } Hohenlohe et al (2010), identical to the STACKS
package. \item{Genepop: } Rousset (2008), uses the genepop package. }}
}
\value{
A snpRdata object with pairwise FST as well as the number of total
 observations at each SNP in each comparison merged in to the pairwise.stats
 socket.
}
\description{
\code{calc_pairwise_fst} calculates pairwise FST for each SNP for each
possible pairwise combination of populations.
}
\details{
Calculates FST according to either Wier and Cockerham 1984, Wier 1990,
Hohenlohe et al 2010, or using the \code{\link[genepop]{Fst}} function from
the genepop package (see references).

If the genpop option is used, several intermediate files will be created in
the current working directory, which are not cleaned automatically.
\code{calc_pairwise_fst} will ask for permission to continue if these files
already exist.

The Wier and Cockerham (1984), Wier (1990), and genepop methods tend to
produce very similar results. Generally, either of the two former options are
prefered for computational efficieny.

The data can be broken up categorically by either SNP and/or sample metadata,
as described in \code{\link{Facets_in_snpR}}. Since this is a pairwise
statistic, at least a single sample level facet must be provided.

Method Options: \itemize{ \item{"WC": }{Wier and Cockerham 1984.}
\item{"Wier": }{Wier 1990.} \item{"Genepop": }{As used in genepop, Rousset
2008.} \item{"Hohenlohe": }{Hohenlohe 2010.} }
}
\examples{
# Using Wier and Cockerham 1984's method
x <- calc_pairwise_fst(stickSNPs, "pop")
get.snpR.stats(x, "pop", "pairwise")

# Using genepop, note that the overall value is part 2 of the returned list
x <- calc_pairwise_fst(stickSNPs, "pop", "genepop")
get.snpR.stats(x[[1]], "pop", "pairwise")
x[[2]] # overall fst
}
\references{
Wier and Cockerham (1984). \emph{Evolution}

Wier (1990). Genetic data analysis. Sinauer,  Sunderland, MA

Hohenlohe et al. (2010). \emph{PLOS Genetics}

Rousset (2008). \emph{Molecular Ecology Resources}
}
\author{
William Hemstrom
}
