% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/colony_interface.R
\name{colony_interface}
\alias{colony_interface}
\alias{write_colony_input}
\alias{call_colony}
\alias{parse_colony}
\alias{run_colony}
\title{Interfaces with the COLONY pedigree assignment program.}
\usage{
write_colony_input(
  x,
  outfile = "colony_input",
  method = "FPLS",
  run_length = 2,
  sampleIDs = NULL,
  sibship_prior = 0,
  paternal_sib_size = NULL,
  maternal_sib_size = NULL,
  nruns = 1,
  seed = NULL,
  maternal_genotypes = NULL,
  paternal_genotypes = NULL,
  maternal_inclusion_prob = 0,
  paternal_inclusion_prob = 0,
  update_af = TRUE,
  dioecious = TRUE,
  inbreeding = TRUE,
  male_monogamous = FALSE,
  female_monogamous = FALSE,
  clone_inference = FALSE,
  sibship_scaling = TRUE,
  known_af = FALSE,
  precision = 2,
  dropout = 0.01,
  genotyping_error = 0.01,
  known_maternal_dyads = NULL,
  known_paternal_dyads = NULL,
  known_maternal_max_mismatches = 0,
  known_paternal_max_mismatches = 0,
  known_maternal_sibships = NULL,
  known_paternal_sibships = NULL,
  maternal_exclusions = NULL,
  paternal_exclusions = NULL,
  excluded_maternal_siblings = NULL,
  excluded_paternal_siblings = NULL,
  update_bib = FALSE
)

call_colony(infile, colony_path, update_bib = FALSE, verbose = TRUE)

parse_colony(prefix, x, path = "colony/", sampleIDs = NULL, cleanup = FALSE)

run_colony(
  x,
  colony_path,
  outfile = "colony_input",
  method = "FPLS",
  run_length = 2,
  sampleIDs = NULL,
  sibship_prior = 0,
  paternal_sib_size = NULL,
  maternal_sib_size = NULL,
  nruns = 1,
  seed = NULL,
  maternal_genotypes = NULL,
  paternal_genotypes = NULL,
  maternal_inclusion_prob = 0,
  paternal_inclusion_prob = 0,
  update_af = TRUE,
  dioecious = TRUE,
  inbreeding = TRUE,
  male_monogamous = FALSE,
  female_monogamous = FALSE,
  clone_inference = FALSE,
  sibship_scaling = TRUE,
  known_af = FALSE,
  precision = 2,
  dropout = 0.01,
  genotyping_error = 0.01,
  known_maternal_dyads = NULL,
  known_paternal_dyads = NULL,
  known_maternal_max_mismatches = 0,
  known_paternal_max_mismatches = 0,
  known_maternal_sibships = NULL,
  known_paternal_sibships = NULL,
  maternal_exclusions = NULL,
  paternal_exclusions = NULL,
  excluded_maternal_siblings = NULL,
  excluded_paternal_siblings = NULL,
  update_bib = FALSE,
  cleanup = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{x}{snpRdata object from which metadata for colony results can be found.}

\item{outfile}{character, default "colony_input". Output file name. A file
path may be provided (e.g. "colony/colony_run_1.txt").}

\item{method}{character, default "FPLS". Pedigree reconstruction method. For
more details see the Colony User Guide. Options: \itemize{\item{"FPLS": }
{Pure pairwise likelihood method, combines the full likelihood and
pairwise likelihood methods. A good compromise between speed and
accuracy.}\item{"FL": }{Full Likelihood. More accurate than PLS but more
computationally intensive and slow to run, especially with large complex
datasets.}\item{"PLS": }{Pairwise likelihood score. Less accurate but less
computationally intensive than FL.}}}

\item{run_length}{numeric in c(0,1,2,3), default 2. Length of run:
short/medium/long/very long.}

\item{sampleIDs}{character, default NULL. Name of a column in the sample
metadata that designates sample identifications/"names". Each name must be
unique!}

\item{sibship_prior}{numeric in c(0, 1, 2, 3, 4), default 0. Strength the
sibship size prior (no prior, weak, medium, strong, OR determine from known
prior values). Values other than 0 require additional parameters. Option
for supplying value of 4, currently not implemented in snpR. See Colony
User Guide for more details.}

\item{paternal_sib_size}{numeric, default NULL. Minimum value is 0. The
number of offspring in the candidate pool that are known to share the same
father. If this value is not zero, then you must include a file with the
known paternal sibship/paternity.}

\item{maternal_sib_size}{numeric, default NULL. Minimum value is 0. The
number of offspring in the candidate pool that are known to share the same
mother. If this value is not zero, then you must include a file with the
known paternal sibship/maternity.}

\item{nruns}{integer, default 1. A number of replicate runs for the dataset.}

\item{seed}{integer, default NULL. Supply a four digit integer (eg: 1234,
9876) as a starting point for the algorithm.}

\item{maternal_genotypes}{snpRdata object containing maternal genotypes.}

\item{paternal_genotypes}{snpRdata object containing paternal genotypes.}

\item{maternal_inclusion_prob}{numeric in 0:1, default 0. Probability the
mother is in the dataset ranging from 0 to 1.}

\item{paternal_inclusion_prob}{numeric in 0:1, default 0. Probability the
father is in the dataset ranging from 0 to 1.}

\item{update_af}{character, default TRUE. Should Colony update the allele
frequencies used in the calculations?}

\item{dioecious}{character, default TRUE. Is this species diploid/dioecious?
FALSE = haploid/monoecious. Colony does not work with more ploidy.}

\item{inbreeding}{character, default TRUE. Should Colony assume inbreeding in
the calculations?}

\item{male_monogamous}{character, default FALSE. Should Colony assume males
are monogamous?}

\item{female_monogamous}{character, default FALSE. Should Colony assume
females are monogamous?}

\item{clone_inference}{character, default FALSE. Should Colony infer clones
in the sample set?}

\item{sibship_scaling}{character, default TRUE. Should Colony scale sibling
groups?}

\item{known_af}{character, default FALSE. If TRUE snpR will calculate and
supply mafs, else, supply a numeric vector containing the known maf for
each locus.}

\item{precision}{integer in c(0,1,2,3), default 2. Low/Medium/High/Very High
for calculating the maximum likelihood.}

\item{dropout}{numeric vector where each value is in 0:1, default 0.01.
Supply a flatrate value for all markers, or a vector corresponding to the
allelic droput rate for each marker.}

\item{genotyping_error}{numeric vector where each value is in 0:1, default
0.01. Supply a flatrate value for all markers, or a vector corresponding to
the genotyping error rate for each marker.}

\item{known_maternal_dyads}{character, default NULL. Supply matrix or
dataframe with known maternal-offspring dyads. Offspring ID in column 1,
Maternal ID in column 2.}

\item{known_paternal_dyads}{character, default NULL. Supply matrix or
dataframe with known paternal-offspring dyads. Offspring ID in column 1,
Paternal ID in column 2.}

\item{known_maternal_max_mismatches}{integer in c(0,1,2:nsample), default 0.}

\item{known_paternal_max_mismatches}{integer in c(0,1,2:nsample), default 0.}

\item{known_maternal_sibships}{character, default NULL. Data frame or matrix
with sibship size followed by single column containing all of the sibling
IDs separated by spaces.}

\item{known_paternal_sibships}{character, default NULL. Data frame or matrix
with sibship size followed by single column containing all of the sibling
IDs separated by spaces.}

\item{maternal_exclusions}{character, default NULL. Data.frame or matrix with
column 1 the offspring ID, column 2 the number of excluded females, column
3 the IDs of excluded females separated by spaces.}

\item{paternal_exclusions}{character, default NULL. Data.frame or matrix with
column 1 the offspring ID, column 2 the number of excluded males, column 3
the IDs of excluded males separated by spaces.}

\item{excluded_maternal_siblings}{character, default NULL. Data.frame or
matrix with column 1 the offspring ID, column 2 the number of excluded
siblings, column 3 the IDs of excluded siblings separated by spaces.}

\item{excluded_paternal_siblings}{character, default NULL. Data.frame or
matrix with column 1 the offspring ID, column 2 the number of excluded
siblings, column 3 the IDs of excluded siblings separated by spaces.}

\item{update_bib}{character or FALSE, default FALSE. If a file path to an
existing .bib library or to a valid path for a new one, will update or
create a .bib file including any new citations for methods used. Useful
given that this function does not return a snpRdata object, so a
\code{\link{citations}} cannot be used to fetch references.}

\item{infile}{character. Path to the pre-written colony input file to be run.}

\item{colony_path}{character. Path to the colony executable.}

\item{verbose}{Logical, default TRUE. If TRUE, colony progress will be
reported to console.}

\item{prefix}{character. The prefix for the colony files to be parsed.}

\item{path}{character. Path to the directory containing colony results.}

\item{cleanup}{logical, default FALSE. If TRUE, colony files will be removed
after parsing.}
}
\description{
Interfaces with the command-line version of the COLONY pedigree program.
Requires a snpRdata object containing offspring genotypes and can optionally
take snpRdata objects containing maternal or paternal genotypes, or both. No
facet support.

Create a COLONY infile using snpRdata sets containing offspring and possibly paternal genotypes given
specified parameters. Requires that the COLONY program is installed locally. Input and output files
will be stored in a colony folder created in the current working directory.
}
\details{
Requires that the COLONY program is installed locally if running within snpR.
Text files exported from write_colony_input command can also be imported to
command line versions of COLONY (eg. on a compute cluster). Input and output 
files will be stored in a colony folder created in the current working 
directory. The functions documented here can write input files, call them
using COLONY, and parse some parts of the results given the original
snpRdata object. Note that no facet support is currently available here due
to the complexity and number of possible input parameters that are difficult
to handle across multiple facets and facet levels. Facet support for the
basic, default operation with few additional options may be added in the
future. For now, facets can be handled during parentage and pedigree creation
in snpR using the \code{\link{run_sequoia}} function, which runs a notably
simpler (with respect to implementation) pedigree toolkit.

These functions include many commonly used options but not all possible
parameters for COLONY inputs. See COLONY User Guide for using extra features.

This is still in development. The defaults and the most commonly used options
have been tested and work well, but some of the more esoteric options haven't
been fully tested yet.
}
\section{Functions}{
\itemize{
\item \code{write_colony_input()}: Create a colony input file

\item \code{call_colony()}: Call a colony executable to run a prepared colony input file.

\item \code{parse_colony()}: Parse a previously run colony analysis.

\item \code{run_colony()}: run colony on a snpRdata object, start to finish

}}
\examples{

# A simple example for running all individuals in the snpR object as siblings
# in colony. Not run to avoid clutter.
\dontrun{
  write_colony_input(x = stickSNPs, outfile = "stk.col")
 }

# A more complex example requires 1) creating and adding variables to the
# stickSNPs example dataset and 2) creating subset snpR objects.
a <- 2013:2015 #create a vector of possible birthyears
b <- c("M", "F", "U") #create a vector of possible sexes
stk <- stickSNPs
sample.meta(stk)$BirthYear <-
    sample(x = a, size = nsamps(stk), replace = TRUE) #create birthyears
sample.meta(stk)$ID <- 1:nsamps(stk) #create unique sampleID
sample.meta(stk)$Sex <-
    sample(x= b, size = nsamps(stk), replace = TRUE) # create sexes
#generating snpR objects for male and female potential parents and offspring
# (no U sexes in the potential parents in for this example)
# get list of samples which are now "M" for subsetting
lsir <- which(stk@sample.meta$Sex =="M" &
              stk@sample.meta$BirthYear == "2013") #list sires
ldam <- which(stk@sample.meta$Sex =="F" &
              stk@sample.meta$BirthYear == "2013") #list dams
loff <- which(stk@sample.meta$BirthYear \%in\% c("2014","2015")) #list offspr


#creating new snpR objects for Colony formatting
sir <- subset_snpR_data(x = stk, .samps = lsir)
dam <- subset_snpR_data(x = stk, .samps = ldam)
off <- subset_snpR_data(x = stk, .samps = loff)
# not run to avoid clutter
\dontrun{
  write_colony_input(x = off, outfile = "parents_example.col",
                     maternal_genotypes = dam, paternal_genotypes = sir)
}

# running a simple model
\dontrun{
  ## intentionally shorter run, with a small subset of the samples
  asp <- which(stickSNPs@sample.meta$pop == "ASP")
  test_dat <- subset_snpR_data(stickSNPs, .samps = asp)
  run_colony(x = test_dat, colony_path = "/usr/bin/colony2s.exe",
             method = "PLS", run_length = 1)
}
}
\references{
Jones,O.R. and Wang,J. (2010) COLONY: a program for parentage and
  sibship inference from multilocus genotype data. Mol. Ecol. Resour., 10,
  551–555.
}
\author{
William Hemstrom

Melissa Jones
}
