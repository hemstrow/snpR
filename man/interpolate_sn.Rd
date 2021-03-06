% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internals.R
\name{interpolate_sn}
\alias{interpolate_sn}
\title{Interpolate sn formatted data.}
\usage{
interpolate_sn(sn, method = "bernoulli", ncp = NULL, ncp.max = 5)
}
\arguments{
\item{sn}{data.frame. Input sn formatted data, as produced by
\code{\link{format_snps}}. Note that \code{\link{format_snps}} has an option
to automatically call this function during formatting.}

\item{method}{character, default "bernoulli". Method to used for
interpolation, either bernoulli or af. See details.}

\item{ncp}{numeric or NULL, default NULL. Number of components to consider for
iPCA sn format interpolations of missing data. If null, the optimum number
will be estimated, with the maximum specified by ncp.max. This can be very
slow.}

\item{ncp.max}{numeric, default 5. Maximum number of components to check for
when determining the optimum number of components to use when interpolating
sn data using the iPCA approach.}
}
\description{
An internal function to interpolate sn formatted data using either the
bernoulli or expected minor allele count approaches. Typically entirely
internal, called via format_snps.
}
\details{
Interpolating missing data in sn formatted data is useful for PCA, genomic
prediction, tSNE, and other methods. Specify method = "af" to insert the
expected number of minor alleles given SNP allele frequency or "bernoulli" to
do binomial draws to determine the number of minor alleles at each missing
data point, where the probability of drawing a minor allele is equal to the
minor allele frequency. The expected number of minor alleles based on the
later method is equal to the interpolated value from the former, but the later
allows for multiple runs to determine the impact of stochastic draws and is
generally prefered and required for some downstream analysis. It is therefore
the default. As a slower but more accurate alternative to "af" interpolation,
"iPCA" may be selected. This an iterative PCA approach to interpolate based on
SNP/SNP covariance via \code{\link[missMDA]{imputePCA}}. If the ncp arugment
is not defined, the number of components used for interpolation will be
estimated using \code{\link[missMDA]{estim_ncpPCA}}. In this case, this method
is much slower than the other methods, especially for large datasets. Setting
an ncp of 2-5 generally results in reasonable inpterpolations without the time
constraint.
}
\author{
William Hemstrom
}
