% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequoia_interface.R
\name{run_sequoia}
\alias{run_sequoia}
\title{Run Sequoia pedigree/parentage assignment with snpR}
\usage{
run_sequoia(
  x,
  facets = NULL,
  run_dupcheck = FALSE,
  run_parents = FALSE,
  run_pedigree = FALSE,
  run_relatives = FALSE,
  min_maf = 0.3,
  min_ind = 0.5,
  ...
)
}
\arguments{
\item{x}{snpRdata object.}

\item{facets}{character, default NULL. Sample-specific facets over which the
sequoia is called to run. See \code{\link{Facets_in_snpR}}.}

\item{run_dupcheck}{FALSE or TRUE, default FALSE. Uses sequoia to check for 
duplicate samples in the dataset. Duplicate samples should not be included 
for parentage and pedigree construction.}

\item{run_parents}{FALSE or TRUE, default FALSE. Runs parentage assignments 
for the samples. This runs quickly and is required before using the 
run_pedigree command.}

\item{run_pedigree}{FALSE or TRUE, default FALSE. Runs pedigree construction 
for the samples. This process can take a long time.}

\item{run_relatives}{FALSE or TRUE, default FALSE. Runs retrieval of other
relatives which did not pass thresholds for assignment in the main pedigree 
construction model.}

\item{min_maf}{numeric in 0.25:0.5, default 0.3. Minimum allele frequency
cutoff for analysis. Sequoia requires high minor allele frequencies for
parentage and pedigree construction.}

\item{min_ind}{numeric in 0.5:1, default 0.5. Removes loci sequenced in less
than this proportion of individuals. Note that \emph{individuals} with
genotypes for fewer than half of the loci will be automatically removed by
sequoia.}

\item{...}{Additional arguments passed to\code{sequoia}
(during parentage and pedigree reconstruction).}
}
\value{
A nested list with each facet specified containing sequoia output 
  summary information.
}
\description{
Runs the parentage assignment and pedigree construction tool from the
\code{sequoia} package. Note that this function \emph{is not overwrite
safe!}.
}
\details{
This is a limited integration of the program and package written by Jisca 
Huisman. Note that there are many more Sequoia specific arguments that can 
be added to change from the default settings (eg. ErrorM, Tassign, Tfilt, 
etc.) See documentation for \code{sequoia}. These can be passed to the 
pedigree and parentage reconstructions using the ... argument in run_sequoia.
The sequoia package has many features, and snpR facilitates the use of a 
fraction of them. snpR users are encouraged to use the sequoia R package.
}
\examples{
# to follow an example using the stickSNPs example dataset you need 
# to add some variables that don't exist in the actual dataset.
a <- 2013:2015 #create a vector of possible birthyears
b <- c("M", "F", "U") #create a vector of possible sexes
stk <- stickSNPs
set.seed(4865)
sample.meta(stk)$BirthYear <- sample(x = a, size = nsamps(stickSNPs), 
                                     replace = TRUE) #create birth years
sample.meta(stk)$ID <- 1:nsamps(stk) #create unique sampleID
sample.meta(stk)$Sex <- sample(x= b, size = nsamps(stk), 
                               replace = TRUE) # create sexes

# slow, so not run here
\dontrun{
dup <- run_sequoia(x = stk, run_dupcheck = TRUE, run_parents = FALSE, 
                   run_pedigree = FALSE, run_relatives = FALSE)
ped <- run_sequoia(x = stk, run_dupcheck = FALSE, run_parents = TRUE, 
                   run_pedigree = TRUE, run_relatives = FALSE)
rel <- run_sequoia(x = stk, run_dupcheck = FALSE, run_parents = FALSE, 
                   run_pedigree = FALSE, run_relatives = TRUE)
                   
}
}
\references{
Huisman,J. (2017) Pedigree reconstruction from SNP data:
  parentage assignment, sibship clustering and beyond. Mol. Ecol. Resour.,
  17, 1009–1024.
}
\author{
William Hemstrom

Melissa Jones
}
