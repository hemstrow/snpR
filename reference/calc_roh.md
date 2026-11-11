# Runs of Homozygosity and fROH.

Identifies Runs of Homozygosity and calculates ROH-based genomic
inbreeding (FROH) using (approximately) "PLINK"'s methodology.

## Usage

``` r
calc_roh(
  x,
  facets = NULL,
  window_snps = 50,
  min_snps = 100,
  min_length = 1000,
  max_hets = 1,
  max_missing = 5,
  gap = 1000,
  min_density = 50,
  roh_threshold = 0.05,
  extend = TRUE,
  final_roh_max_het = Inf,
  window_start_method = "conservative",
  genome_length = NULL,
  par = FALSE,
  verbose = 1
)
```

## Arguments

- x:

  snpRdata object

- facets:

  facets over which to calculate ROHs. See
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md)
  for details. Only one snp facet allowed, sample facets are ignored
  except for calculating mean fROH.

- window_snps:

  numeric, default 50. Number of SNPs included in each sliding window
  used to identify possible ROH candidate SNPs. Windows slide by one SNP
  across chromosomes.

- min_snps:

  numeric, default 100. Minimum number of SNPs needed to call a final
  ROH.

- min_length:

  numeric, default 1000. Number of *kb* an ROH must span to be called a
  final ROH.

- max_hets:

  numeric, default 1. Maximum number of heterozygotes a a scanning
  window can contain to be mark containing SNPs as candidates for an
  ROH.

- max_missing:

  numeric, default 1. Maximum number of mising genotypes a a scanning
  window can contain to be mark containing SNPs as candidates for an
  ROH.

- gap:

  numeric or FALSE, default 1000. If not FALSE, final ROHs will be split
  wherever a gap of more than this many base pairs between SNPs is
  detected.

- min_density:

  numeric, default 50. An ROH must contain at least one SNP every
  `min_density` *kb* on average in order to be kept.

- roh_threshold:

  numeric, default 0.05. The percentage of scanning windows for which a
  SNP must pass ROH thresholds in order to be considered part of a final
  ROH. This number is low by default in order to correctly detect the
  boundaries of ROHs.

- extend:

  logical, default TRUE. If TRUE, final ROHs will be extended to include
  any adjacent homozygous loci which were previously not considered part
  of the ROH. Extension stops when either a missing genotype or
  heterozygous loci is detected, or where the distance between loci is
  longer than `gap` (if `gap` is not FALSE).

- final_roh_max_het:

  numeric, default Inf. The number of allowable heterozygotes in the
  final ROHs. ROHs with more than this number of hets after final
  filtering and extension will be removed.

- window_start_method:

  character, default "powerful". Controls how the initial sliding window
  start and end positions are determined. Options:

  - conservative: Windows start and end such that the fist and last SNPs
    of each chromosome have evenly fewer windows which cover them than
    the rest of the genome, but each window has exactly `window_snps`
    SNPs. This generally produces results more similar to PLINK's
    "–homozyg" function.

  - powerfull: Windows start and end such that the start and end of
    chromosomes have the same number of windows covering them as the
    rest of the genome, but those windows will have fewer SNPs within
    them. This will often contain more windows which get close to the
    chromosomal edges that might otherwise be dropped, even if they also
    extend well into the genome.

- genome_length:

  numeric, default NULL. The length of the genome, for use in
  calculating FROH. If NULL, this is calculated by summing the distance
  between the first and last sequenced SNPs of each
  chromosome/scaffold/etc.

- par:

  numeric or FALSE, default FALSE. If numeric, the number of cores to
  use for parallel processing. Parallel processing is done across
  chromosomes/ scaffolds/etc, as provided to `facets`.

- verbose:

  1, 2, or FALSE, default 1. If FALSE, no messages are displayed during
  run. If 1, chromosomes/scaffolds/etc run beginnings are reported; if
  2, progress withing chromosomes/scaffolds/etc are also reported.

## Value

A snpRdata object with ROH and FROH info merged in, fetchable with
[`get.snpR.stats`](https://hemstrow.github.io/snpR/reference/get.snpR.stats.md)

## Details

ROHs are detected in a two-pass approach. First, the genome is scanned
in sliding windows across the genome (windows slid SNP-by-SNP). A window
is considered a candidate ROH if it contains less than `max_hets`
heterozygotes. Windows with more than `max_missing` missing genotype
calls in an individual are discarded. Windows are defined in length by
numbers of SNPs (`window_snps`), with some SNPs ommited at the start and
end of chromosomes if the number of SNPs in the chromosome is not evenly
divisible by window size. SNPs within individuals are given a
probability of being in an ROH given by number of windows they were
included in which were passed the maximum heterozygote count check
divided by the number of SNPs per window, which for most SNPs is the
same as the number of windows they were included in. SNPs at the very
beginning of a chr will not be in as many windows and will thus be
slightly less likely to be detected as in an ROH.

Second, the genome is then slid across SNP-by-SNP to identify ROHs. SNPs
which pass the `roh_threshold` are considered as in ROHs, and chains of
SNPs in ROH are considered as ROHs. Final ROHs are filtered by base-pair
length, number of SNPs, and SNP density, optionally split if gaps
between SNPs are too large, and then optionally extended to cover any
additional homozygous loci which were initially omitted that are not too
far away.

In-all, this approach is generally very similar to that used by PLINK
and is based on the description of PLINK's "–homozyg" function. There
may be some minor differences in window starting and ending locations
due to variation in computation order and extension process, and a few
other details, but the results should usually be very comparable.

FROH is calculated by dividing the length of individual genomes in ROHs
by the total sequenced genome size.

See PLINK's –homozyg documentation for more details.

## References

Purcell, S., Neale, B., Todd-Brown, K., Thomas, L., Ferreira, M. A. R.,
Bender, D., Maller, J., Sklar, P., de Bakker, P. I. W., Daly, M. J., &
Sham, P. C. (2007). PLINK: A Tool Set for Whole-Genome Association and
Population-Based Linkage Analyses. The American Journal of Human
Genetics, 81(3), 559–575. https://doi.org/10.1086/519795

Kardos, M., Luikart, G., & Allendorf, F. W. (2015). Measuring individual
inbreeding in the age of genomics: Marker-based measures are better than
pedigrees. Heredity, 115(1), 63–72. https://doi.org/10.1038/hdy.2015.17

## Author

William Hemstrom

## Examples

``` r
if (FALSE) { # \dontrun{
# not run--the example dataset (stickSNPs) has no ROHs. This is an
# example which could be run with a different dataset, "d".
d <- calc_roh(d, "chr") # run with defaults
res <- get.snpR.stats(d, "chr", "roh")
} # }

```
