% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataset_documentation.R
\docType{data}
\name{stickSNPs}
\alias{stickSNPs}
\title{Deschutes stickleback example snpRdata.}
\format{
An object of class \code{snpRdata} with 100 rows and 100 columns.
}
\usage{
stickSNPs
}
\description{
A dataset containing genotypes at 100 SNP loci for 100 Three-spined
Stickleback collected in 2015 from several locations in the Deschutes River,
Oregon, USA. For more information, see Hemstrom et al (in prep).
}
\details{
Stored in the snpRdata object class. See \code{\link{import.snpR.data}} for
details.

Genotypes noted as two characters representing bases sequenced at a specific
site (e.g. "CC" or "GT"). Each SNP occupies one row, each sample occupies on
column. Missing data coded as "NN".

Sample and snp metadata can be found in the sample.meta and snp.meta sockets,
respectively.
}
\references{
Hemstrom et al (in prep)
}
\keyword{datasets}
