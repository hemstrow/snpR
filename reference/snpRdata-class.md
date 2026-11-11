# Storage class for snpR data and calculated statistics.

The snpRdata class stores both raw genotype data, sample and locus
specific metadata, useful data summaries, repeatedly internally used
tables, calculated summary statistics, and smoothed statistic data. Used
by most snpR functions.

## Details

The snpRdata class is built to contain SNP genotype data for use by
functions in the snpR package. It also stores sample and locus specific
metadata, genomic summary information, and any results from most snpR
functions. Results can be gathered using
[`get.snpR.stats`](https://hemstrow.github.io/snpR/reference/get.snpR.stats.md).
Genotypes are stored in the "character" format, as output by
format_snps().

For more information, see
[`import.snpR.data`](https://hemstrow.github.io/snpR/reference/import.snpR.data.md),
the constructor function for this object class.

## Slots

- `sample.meta`:

  data.frame containing sample metadata

- `snp.meta`:

  data.frame containing snp metadata

- `facet.meta`:

  data.frame containing snp metadata for the geno.tables (in the order
  of snps in that data)

- `mDat`:

  character, missing data key

- `snp.form`:

  numeric, number of characters per genotype (not really used)

- `geno.tables`:

  list containing three matrices: gs (genotype counts), as (allele
  counts), and wm (missing counts)

- `facets`:

  character, vector of tabulated facets

- `facet.type`:

  character, types of each tabulated facet.

- `stats`:

  data.frame, all calculated single-snp non-pairwise stats.

- `window.stats`:

  data.frame, all calculated window non-pairwise stats.

- `pairwise.stats`:

  data.frame, all calculated pairwise single snp stats.

- `pairwise.window.stats`:

  data.frame, all calculated window pairwise stats.

- `sample.stats`:

  data.frame, all calculated single sample stats.

- `pop.stats`:

  data.frame, all calculated population level summary stats.

- `pairwise.LD`:

  nested list of matrices, all calculated pairwise LD stats.

- `window.bootstraps`:

  data.frame, all calculated bootstraps

- `sn`:

  list, contains sn formatted data if calculated

- `calced_stats`:

  list, contains information on what statistics have been calculated for
  which facets

- `allele_frequency_matrices`:

  list of matrices, allele frequency matrices if calculated

- `genetic_distances`:

  list of matrices, genetic distance matrices if calculated

- `other`:

  list, contains other miscellaneous calculated statistics that do not
  fit cleanly elsewhere.

## Author

William Hemstrom
