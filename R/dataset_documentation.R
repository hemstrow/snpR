#'Deschutes stickleback example snpRdata.
#'
#'A dataset containing genotypes at 833 SNP loci for 420 Three-spined
#'Stickleback collected in 2015 from several locations in the Deschutes River,
#'Oregon, USA. For more information, see Hemstrom et al (in prep).
#'
#'Stored in the snpRdata object class. See \code{\link{import.snpR.data}} for
#'details.
#'
#'Genotypes noted as two characters repesenting bases sequenced at a specific
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
#'A dataset containing genotypes at 833 SNP loci for 420 Three-spined
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
#'The first three columns contain snp specific meta data: \itemize{ \item{snp: }
#'SNP ID. \item{chr: } Linkage group, essentially chromosome. \item{position:
#'} Position of the SNP along the linkage group, in base pairs. }
#'
#'@references Hemstrom et al (in prep)
"stickRAW"



#' Contains an artificially constructed pedigree based on the 
#' \code{\link{stickSNPs}} dataset. Here samples named 1 through 420 represent 
#' the sequenced individuals in the \code{\link{stickSNPs}} dataset and 
#' individuals inferred via pedigree construction are affixed with "f" or "m". 
"stickPED"