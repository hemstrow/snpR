#'Deschutes stickleback example snpRdata.
#'
#'A dataset containing genotypes at 100 SNP loci for 100 Three-spined
#'Stickleback collected in 2015 from several locations in the Deschutes River,
#'Oregon, USA. For more information, see Hemstrom et al (in prep).
#'
#'Stored in the snpRdata object class. See \code{\link{import.snpR.data}} for
#'details.
#'
#'Genotypes noted as two characters representing bases sequenced at a specific
#'site (e.g. "CC" or "GT"). Each SNP occupies one row, each sample occupies on
#'column. Missing data coded as "NN".
#'
#'Sample and snp metadata can be found in the sample.meta and snp.meta sockets,
#'respectively.
#'
#'@references Hemstrom et al (in prep)
"stickSNPs"

#'Deschutes stickleback raw data.
#'
#'A dataset containing genotypes at 100 SNP loci for 100 Three-spined
#'Stickleback collected in 2015 from several locations in the Deschutes River,
#'Oregon, USA. For more information, see Hemstrom et al (in prep).
#'
#'Genotypes noted as two characters representing bases sequenced at a specific
#'site (e.g. "CC" or "GT"). Each SNP occupies one row, each sample occupies on
#'column. Missing data coded as "NN".
#'
#'Column names refer to the population of origin (first three characters) and
#'the sample ID.
#'
#'The first three columns contain snp specific meta data:
#'* snp: SNP ID. 
#'* chr:  Linkage group, essentially chromosome. 
#'* position: Position of the SNP along the linkage group, in base pairs.
#'
#'@references Hemstrom et al (in prep)
"stickRAW"



#' An artificial pedigree from the stickleback dataset.
#' 
#' Contains an artificially constructed pedigree based on the 
#' \code{\link{stickSNPs}} dataset. Here samples named 1 through 420 represent 
#' the sequenced individuals in the \code{\link{stickSNPs}} dataset and 
#' individuals inferred via pedigree construction are affixed with "f" or "m". 
"stickPED"

#' Contains example data for 13 microsatellite loci from 1573 steelhead from the
#' Siletz river, Oregon, USA. Samples come from both the summer and winter run
#' wild and hatchery populations. For more details, see Hemstrom et al 2018.
#'
#' The first column contains population info, each remaining column contains
#' information on one microsatellite. Note that this data is transposed with
#' respect to the SNP example data.
#'
#' Individual genotypes are coded as "0000", where the first two numbers are an
#' identifer for the first allele and the third and fourth identify the second
#' allele.
"steelRAW"

#' Siletz Steelhead microsatellite snpRdata
#'
#' Contains example data for 13 microsatellite loci from 1573 steelhead from the
#' Siletz river, Oregon, USA. Samples come from both the summer and winter run
#' wild and hatchery populations. For more details, see Hemstrom et al 2018.
#'
#' Sample metadata can be found using the \code{sample.meta} function.
#' Individuals have been randomly sorted into "families" for example purposes.
"steelMSATs"