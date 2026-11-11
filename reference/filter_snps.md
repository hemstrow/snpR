# Filter SNPs in snpRdata objects.

`filter_snps` filters snpRdata objects to remove SNPs or individuals
which fail to pass user defined thresholds for several statistics. Since
this function removes all calculated statistics, etc. from the snpRdata
object, this should usually be the first step in an analysis. See
details for filters.

## Usage

``` r
filter_snps(
  x,
  maf = FALSE,
  mac = 0,
  mgc = 0,
  hf_hets = FALSE,
  hwe = FALSE,
  fwe_method = "none",
  hwe_excess_side = "both",
  singletons = FALSE,
  min_ind = FALSE,
  min_loci = FALSE,
  inds_first = FALSE,
  remove_garbage = FALSE,
  re_run = "partial",
  maf_facets = NULL,
  hwe_facets = NULL,
  non_poly = TRUE,
  bi_al = TRUE,
  LD_prune_sigma = FALSE,
  LD_prune_step = LD_prune_sigma,
  LD_prune_r = FALSE,
  LD_prune_method = "CLD",
  LD_prune_facet = NULL,
  LD_prune_par = FALSE,
  LD_prune_use_ME = FALSE,
  LD_prune_ME_sigma = 1e-04,
  verbose = TRUE
)
```

## Arguments

- x:

  snpRdata object.

- maf:

  numeric between 0 and 1 or FALSE, default FALSE. Minimum acceptable
  minor allele frequency. Either maf or mac filtering are allowed, not
  both.

- mac:

  integer, between 0 and 1 minus the total number of individuals. Loci
  with less *or equal to* this many minor alleles will be removed.
  Either maf or mac filtering are allowed, not both. `mac = 1` is
  therefore a singleton filter.

- mgc:

  integer, between 0 and the total number of individuals/2. Loci where
  the minor allele is present in less *or equal to* this many
  individuals will be removed. Either mac or mgc filtering are allowed,
  not both. `mgc = 1` is analagous to a singleton filter, but will also
  remove loci with one homozygous minor individual.

- hf_hets:

  numeric between 0 and 1 or FALSE, default FALSE. Maximum acceptable
  heterozygote frequency.

- hwe:

  numeric between 0 and 1 or FALSE, default FALSE. Minimum acceptable
  HWE p-value.

- fwe_method:

  character, default "none". Option to use for Family-Wise Error rate
  correction for HWE filtering. If requested, only SNPs with p-values
  below the alpha provided to the `hwe` argument *after FWE correction*
  will be removed. See
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) for information on
  method options.

- hwe_excess_side:

  character, default "both". Options:

  - heterozygote: only loci with a heterozygote excess are removed.

  - homozygote: only loci with a homozygote excess are removed.

  - both: loci with either a heterozygote or homozygote excess are
    removed.

- singletons:

  logical, default FALSE. Depricated, use `mac = 1` to remove
  singletons. If TRUE, removes singletons (loci where there is only a
  single minor allele). If population sizes are reasonably high, this is
  more or less redundant if `maf` is also set.

- min_ind:

  numeric between 0 and 1 or FALSE, default FALSE. Minimum proportion of
  individuals in which a loci must be sequenced.

- min_loci:

  numeric between 0 and 1 or FALSE, default FALSE. Minimum proportion of
  SNPs at which an individual must be genotyped.

- inds_first:

  logical, default FALSE. If TRUE, individuals will be filtered out for
  missing data if that option is selected prior to loci being filtered.
  Otherwise loci are filtered first.

- remove_garbage:

  numeric between 0 and 1 or FALSE, default FALSE. Optionally do a
  filter to remove very poorly sequenced individuals and loci jointly
  before applying other filters. This can be used to reduce biases
  caused by very bad loci/individuals prior to full filtering. This
  number should be lower than either the `min_ind` or `min_loci`
  parameters if those arguments are used and should generally be quite
  permissive to remove only truly bad loci or individuals.

- re_run:

  character or FALSE, default "partial". When individuals are removed
  via min_ind, it is possible that some SNPs that initially passed
  filtering steps will now violate some filters. SNP filters can be
  re-run automatically via several methods:

  - partial: Re-filters for non-polymorphic loci (non_poly) only, if
    that filter was requested initially.

  - full: Re-runs the full filtering scheme (save for min_loci).

  Note: if `inds_first = TRUE`, all re-run options other than FALSE will
  re-run the individual filter for missing data again after loci
  filtering.

- maf_facets:

  character or NULL, default NULL. Defines a sample facet over which the
  minor allele frequency can be checked. SNPs will only fail the maf
  filter if they fail in every level of every provided facet.

- hwe_facets:

  character or NULL, default NULL. Defines a sample facet over which the
  hwe filter can be checked. SNPs will fail the hwe filter if they fail
  in any level of any provided facet.

- non_poly:

  logical, default TRUE. If TRUE, non-polymorphic loci will be removed.

- bi_al:

  logical, default TRUE. If TRUE, loci with more than two alleles will
  be removed. Note that this is mostly an internal argument and should
  rarely be used directly, since import.snpR.data and other snpRdata
  object creation functions all pass SNPs through this filter because
  many snpR functions will fail to work if there are more than two
  alleles at a locus.

- LD_prune_sigma:

  numeric or FALSE, default FALSE. If LD pruning, the window size in kb
  across which to consider LD loci pairs. Windows will be equal to two
  times `LD_prune_sigma`.

- LD_prune_step:

  numeric, default equal to `LD_prune_sigma`. Step size between LD
  windows, in kb. The default ensures that each SNP is in exactly two
  windows, allowing for smoother filtering.

- LD_prune_r:

  numeric, default FALSE. The LD filter cutt-off value. Pairs above this
  value will be greedily pruned. Must be between 0 and 1 if LD prunning
  is conducted.

- LD_prune_method:

  character, default "CLD". The method to use for LD calculations and
  prunning. The options are:

  - CLD

  - Dprime

  - rsq

  See
  [`calc_pairwise_ld`](https://hemstrow.github.io/snpR/reference/calc_pairwise_ld.md)
  for details.

- LD_prune_facet:

  character, default NULL. Defines a facet over which the minor allele
  frequency can be checked. Windows will be calculated within any SNP
  parts of the facet, LD values will be calculated only with levels of
  any sample parts. For example 'chr.pop' would calculate LD values
  across windows, then filter SNPs above `LD_prune_r` in any
  populations.

- LD_prune_par:

  numeric, default FALSE. If numeric, LD calculation will be conducted
  in parallel using `LD_prune_par` threads.

- LD_prune_use_ME:

  logical, default FALSE. If TRUE, uses Minimization-Expectation to do
  LD calculations for all `LD_prune_method` options other than "CLD".
  This can be very slow. See
  [`calc_pairwise_ld`](https://hemstrow.github.io/snpR/reference/calc_pairwise_ld.md)
  for details.

- LD_prune_ME_sigma:

  numeric, default 0.0001. The cutt-off value for difference in
  haplotype frequencies to use when conducting ME. See
  [`calc_pairwise_ld`](https://hemstrow.github.io/snpR/reference/calc_pairwise_ld.md)
  for details.

- verbose:

  Logical, default TRUE. If TRUE, some progress updates and filtering
  notes will be printed to the console.

## Value

A data.frame in the same format as the input, with SNPs and individuals
not passing the filters removed.

## Details

Possible filters:

- maf, minor allele frequency: removes SNPs where the minor allele
  frequency is too low. Can look for mafs below provided either globally
  or search each population individually.

- hf_hets, high observed heterozygosity: removes SNPs where the observed
  heterozygosity is too high.

- min_ind, minimum individuals: removes SNPs that were genotyped in too
  few individuals.

- min_loci, minimum loci: removes individuals sequenced at too few loci.

- non_poly, non-polymorphic SNPs: removes SNPs that are not polymorphic
  (not true SNPs).

- bi_al, non-biallelic SNPs: removes SNPs that have more than two
  observed alleles. This is mostly an internal argument, since the
  various snpRdata import options use it automatically to prevent
  downstream errors in other snpR functions.

- remove_garbage: Quickly removes any loci or samples that are jointly
  poorly genotyped prior to other filtering. This can be useful if some
  individuals or loci are of poor enough quality that they will bias
  other filters.

- LD_prune_sigma: Sliding-window based Linkage Disequilibrium filtering
  to remove non-independant loci.

## Filter order

Note that filtering out poorly sequenced individuals creates a possible
conflict with the loci filters, since after individuals are removed,
some loci may no longer pass filters. For example, if a portion of
individuals in one population all carry the only instances of a rare
minor allele that still passes the maf threshold, removing those
individuals may cause the loci to no longer be polymorphic in the
sample.

To counter this, the `re_run` argument can be used to pass the data
through a second filtering step after individuals are removed. By
default, the "partial" re-run option is used, which re-runs only the
non-polymorphic filter (if it was originally set). The "full" option
re-runs all set filters. Note that re-running any of these filters may
cause individuals to fail the individual filter after loci removal, and
so subsequent tertiary re-running of the individual filters, followed by
the loci filters, and so on, could be justified. This is not done
automatically here, however, if the `inds_first` option is selected to
filter poorly sequenced individuals prior to applying loci filters, any
re-run option will re-run individuals filters after the loci filters
have been applied.

Filter order (snps vs samples first) can be controlled with the
`inds_fist` argument. If individuals are filtered first, poorly
sequenced individuals will be re-filtered following locus removal if
re-running is requested.

## LD pruning

LD pruning can be conducted using any of the methods described in
[`calc_pairwise_ld`](https://hemstrow.github.io/snpR/reference/calc_pairwise_ld.md)
using the `LD_prune_` argument family. Loci are removed greedily:
starting from the first locus pair with elevated LD, the more poorly
sequenced loci is removed. If each locus is equally well sequenced, the
second loci (by position) is removed. If any elevated locus pairs
remain, the next pair is addressed the same way until no pairs remain.

## Facet-based filters (within-group filtering)

Via the `maf_facets` argument, this function can filter by minor allele
frequencies in either *all* samples or *each level of a supplied sample
specific facet and the entire dataset*. In the latter case, any SNPs
that pass the maf filter in *any* facet level are considered to pass the
filter. The latter should be used in instances where population sizes
are very different, there are *many* populations, and/or allele
frequencies are very different between populations and thus common
alleles of interest in one population might be otherwise filtered out.

The `hwe_facets` argument is the inverse of this: loci will be removed
if they fail the provided hwe filter in any facet level. In both cases,
Facets should be provided as described in
[`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).

LD filtering can likewise be conducted within facets. SNP-level facets
dictate the levels within which LD values are calculated (chromosome,
etc.), and loci will be removed if they display elevated LD within
sample-level facets.

## Author

William Hemstrom

## Examples

``` r
# Filter with a minor allele frequency of 0.05, maximum heterozygote 
# frequency of 0.55, 50% minimum individuals, and at least 75% of loci 
# sequenced per individual.
filter_snps(stickSNPs, maf = 0.05, hf_hets = 0.55, 
            min_ind = 0.5, min_loci = 0.75)
#> Initializing...
#> Filtering loci. Starting loci: 100 
#> Filtering non-biallelic loci...
#>  0 bad loci
#> Filtering non_polymorphic loci...
#>  0 bad loci
#> Filtering loci sequenced in few individuals...
#>  0 bad loci
#> Filtering low minor allele frequencies, no pops...
#>  10 bad loci
#> Filtering high frequency heterozygote loci...
#>  2 bad loci
#>  Ending loci: 88 
#> Filtering individuals. Starting individuals: 100 
#> Filtering out individuals sequenced in few kept loci...
#> Re-calculating and adding facets.
#>  Ending individuals: 87 
#> Re-filtering loci...
#> Filtering non_polymorphic loci...
#>  0 bad loci
#>  Final loci count: 88 
#>    ASP.r2.141 ASP.r2.144 ASP.r2.146 ASP.r2.196 ASP.r2.216 ASP.r2.221 ASP.r2.227
#> 1          AG         GG         AA         NN         AG         NN         GG
#> 2          CC         CC         CC         CC         CC         CT         CC
#> 3          AA         AA         AA         AG         AA         AA         AA
#> 4          GG         NN         NN         GG         AG         AG         AG
#> 5          AC         AC         AA         AA         AC         AA         AC
#> 6          TT         TT         TT         GT         TT         TT         TT
#> 7          GG         GG         AG         AG         GG         GG         AG
#> 8          CC         AC         CC         CC         CC         CC         CC
#> 9          AA         AA         AG         AA         AA         AA         AA
#> 10         NN         TT         TT         TT         CT         CC         CC
#> 11         AG         AG         AG         AA         AA         AG         AA
#> 12         TT         TT         CT         CT         NN         TT         TT
#> 13         GG         GT         GT         GG         GG         GG         GG
#> 14         AA         AG         AG         AA         AA         AA         AA
#> 15         NN         GG         GT         GG         GT         GT         GG
#> 16         AA         AA         AA         AA         NN         AC         AC
#> 17         AC         NN         CC         AC         AC         CC         AC
#> 18         TT         TT         TT         TT         TT         TT         TT
#> 19         NN         NN         CC         AC         AC         CC         CC
#> 20         CT         NN         TT         CT         NN         CC         TT
#> 21         AA         AA         AA         AA         AA         AA         AA
#> 22         NN         AG         AG         GG         AG         AG         GG
#> 23         CT         CT         CT         NN         TT         NN         TT
#> 24         AA         AA         AA         AA         AA         AA         AA
#> 25         TT         NN         NN         CT         CT         TT         TT
#> 26         AA         AC         AA         AA         AA         AA         NN
#> 27         GG         GG         CG         GG         GG         GG         GG
#> 28         CG         CG         CC         CC         CC         CG         CC
#> 29         AA         AA         AA         AA         AA         AA         AC
#> 30         CC         CT         CC         CT         CT         CT         CC
#> 31         CG         CG         CG         CC         NN         CC         CG
#> 32         AA         AA         AA         AA         NN         AC         AA
#> 33         AA         AA         AA         AA         AG         AA         AG
#> 34         CG         GG         GG         GG         GG         NN         GG
#> 35         GG         TT         GG         GT         TT         GG         GT
#> 36         TT         NN         TT         AT         TT         TT         TT
#> 37         CG         CC         CG         CG         GG         CG         GG
#> 38         AG         AA         AG         AA         AA         AG         AG
#> 39         AG         NN         NN         GG         GG         AG         AG
#> 40         AG         AG         AA         AA         AA         AA         AA
#> 41         AT         TT         TT         NN         TT         TT         NN
#> 42         NN         CC         AA         NN         CC         NN         AA
#> 43         CC         CC         CC         CG         CC         CG         CG
#> 44         CT         TT         TT         TT         TT         CT         CC
#> 45         NN         CC         CT         CC         CT         CC         CC
#> 46         TT         NN         CT         TT         CT         CT         NN
#> 47         AA         AG         AG         AG         AG         AA         AG
#> 48         AA         AA         CC         AA         AC         AA         AA
#> 49         GG         NN         AG         AG         AG         AG         AG
#> 50         AG         AG         AG         AA         NN         AG         AA
#> 51         TT         CT         CC         TT         CC         NN         CT
#> 52         AA         AG         GG         GG         AA         GG         AG
#> 53         AA         AA         AC         AC         AC         AC         NN
#> 54         CG         GG         GG         GG         GG         GG         GG
#> 55         AC         NN         NN         AA         AC         AC         AC
#> 56         TT         AT         AT         AT         AT         AA         AT
#> 57         CT         CC         CT         CC         CT         CT         CC
#> 58         AC         AC         NN         AC         CC         CC         AC
#> 59         GG         GG         GG         GG         GG         GG         GG
#> 60         AA         AT         AT         AA         NN         AT         TT
#> 61         AA         AA         NN         AA         AA         AA         AA
#> 62         AC         AC         AA         AC         NN         AC         AC
#> 63         AA         AA         AA         AA         AA         AA         AA
#> 64         GG         GG         GG         AG         NN         AG         AG
#> 65         TT         GT         TT         GT         TT         GT         NN
#> 66         GT         GG         GT         GG         GG         GT         GG
#> 67         CC         CT         CT         NN         NN         TT         NN
#> 68         AT         AT         AT         AA         AA         AT         NN
#> 69         GG         GG         GG         GG         GG         GG         GG
#> 70         TT         TT         CT         CT         CT         CT         TT
#> 71         CT         CT         TT         TT         CT         TT         TT
#> 72         CT         TT         CT         TT         TT         NN         TT
#> 73         AT         TT         AT         TT         AT         TT         TT
#> 74         CC         CC         CT         CC         NN         CT         CC
#> 75         AG         AG         GG         AG         NN         AA         AA
#> 76         AG         AG         AG         AA         AA         AG         NN
#> 77         CT         CC         CT         TT         CT         CC         CT
#> 78         AG         NN         AG         GG         AG         AA         AG
#> 79         NN         CG         CG         CC         CG         CG         NN
#> 80         CC         CC         NN         CC         CT         CC         CT
#> 81         AA         NN         GG         GG         AA         AG         AG
#> 82         TT         TT         TT         CT         TT         TT         TT
#> 83         CC         CC         CC         CC         NN         CC         CC
#> 84         GG         GG         GG         GT         GG         GT         GT
#> 85         CC         NN         CC         CC         CC         CC         CC
#> 86         GG         GG         GG         GG         GT         GG         GT
#> 87         AA         GG         AA         AG         AG         AG         AG
#> 88         AT         NN         AT         TT         AT         TT         TT
#>    ASP.r2.230 ASP.r2.232 ASP.r2.235 CLF.r2.040 CLF.r2.048 CLF.r2.061 CLF.r2.064
#> 1          AA         GG         AG         NN         AG         AG         AA
#> 2          CC         CT         CC         CC         CC         CC         CC
#> 3          AA         AG         AG         AA         AG         AA         NN
#> 4          NN         GG         NN         GG         NN         NN         NN
#> 5          AA         AC         NN         AA         AC         AC         AA
#> 6          TT         TT         TT         GG         TT         TT         TT
#> 7          GG         AG         GG         GG         AG         AG         AA
#> 8          AC         CC         AA         CC         AC         CC         AC
#> 9          AG         AA         AA         AA         AA         AA         AA
#> 10         CT         CT         CT         CT         CT         CT         CT
#> 11         AG         AA         AA         AA         AA         AG         AA
#> 12         NN         TT         TT         CT         CT         TT         CT
#> 13         GT         GG         GT         GG         GG         GG         GT
#> 14         AG         AG         AG         AA         AA         AA         AA
#> 15         GG         GG         GG         GT         GT         GG         NN
#> 16         AA         AC         AA         AA         AA         AA         AA
#> 17         AC         CC         AC         CC         AC         AC         NN
#> 18         TT         TT         TT         GT         TT         TT         GT
#> 19         AC         CC         CC         NN         AC         AC         NN
#> 20         CT         CC         CT         CC         CC         NN         TT
#> 21         AA         AA         AG         AA         AA         AA         AA
#> 22         AG         AG         GG         GG         GG         AG         NN
#> 23         CT         TT         CT         TT         TT         TT         TT
#> 24         NN         AC         AA         AA         AC         AA         AA
#> 25         TT         TT         TT         CT         TT         CT         TT
#> 26         AA         AA         AA         AA         NN         NN         AA
#> 27         GG         GG         CG         GG         GG         CG         GG
#> 28         CG         CG         CG         CC         CG         CC         CC
#> 29         AC         AC         AA         AC         AC         NN         AA
#> 30         CC         CT         CT         CT         CC         CT         NN
#> 31         CG         CG         NN         CC         CG         GG         CG
#> 32         AA         AA         AA         AA         NN         AA         AA
#> 33         AA         AG         NN         AG         AA         AA         AG
#> 34         CC         GG         CG         GG         CG         CG         NN
#> 35         NN         GG         TT         GG         GG         GG         GG
#> 36         TT         AT         AT         TT         NN         NN         AT
#> 37         CC         GG         CG         CG         GG         NN         CG
#> 38         GG         AG         AG         AG         GG         GG         AG
#> 39         AG         GG         GG         NN         GG         AG         GG
#> 40         AA         AA         AA         AG         NN         AA         AA
#> 41         TT         AA         TT         TT         TT         TT         TT
#> 42         CC         CC         AC         AC         CC         CC         NN
#> 43         CC         CG         CC         CC         CC         CC         CC
#> 44         NN         TT         CT         CT         CT         TT         CT
#> 45         CT         CC         CC         CC         CC         CC         CC
#> 46         CT         CC         TT         TT         NN         TT         TT
#> 47         AG         AG         AG         GG         AG         AA         AA
#> 48         AC         AA         AA         AC         NN         AA         NN
#> 49         AG         AA         NN         GG         GG         GG         AG
#> 50         AG         AA         NN         AA         AA         AG         AA
#> 51         CC         CT         CT         CC         CC         CC         CC
#> 52         GG         AA         AG         AG         AG         AG         AG
#> 53         NN         CC         NN         AC         AC         CC         CC
#> 54         CG         GG         GG         GG         GG         GG         GG
#> 55         CC         AA         AA         CC         AC         AC         AA
#> 56         TT         TT         TT         TT         TT         TT         TT
#> 57         CC         CC         CC         TT         CC         CT         CT
#> 58         AA         AA         AC         AC         AC         CC         AA
#> 59         NN         GG         GG         GG         GG         GG         GG
#> 60         AA         AA         AA         AT         AA         AA         AT
#> 61         AA         AA         AA         NN         AA         AA         AA
#> 62         AA         AA         AC         CC         AA         AC         NN
#> 63         AA         AA         AA         AA         AA         AA         AA
#> 64         AA         GG         AA         AG         AA         AG         GG
#> 65         GT         GT         GT         GT         TT         GT         TT
#> 66         GG         GG         NN         GG         GT         GG         GG
#> 67         NN         CC         NN         CT         NN         NN         NN
#> 68         AT         NN         AA         AA         AA         AT         AA
#> 69         NN         GG         GG         TT         GG         GG         GG
#> 70         TT         TT         TT         TT         TT         TT         TT
#> 71         CC         TT         TT         CT         CT         TT         TT
#> 72         TT         TT         TT         TT         NN         NN         TT
#> 73         TT         AT         TT         AT         AT         TT         AT
#> 74         CT         CC         TT         TT         CT         CC         CC
#> 75         AG         AA         AG         AG         AA         NN         AG
#> 76         AG         NN         NN         AA         AA         AG         AA
#> 77         TT         TT         CC         TT         CT         TT         TT
#> 78         AG         AG         NN         AA         AA         AG         AA
#> 79         CC         CG         CC         CC         NN         GG         CC
#> 80         CC         NN         CC         CC         NN         TT         CC
#> 81         NN         NN         AG         AG         NN         GG         GG
#> 82         TT         CT         CT         CT         TT         CT         TT
#> 83         CC         NN         NN         CC         CC         CC         CC
#> 84         GG         GT         GG         GT         TT         GG         GG
#> 85         AC         CC         NN         CC         CC         CC         AC
#> 86         NN         GG         TT         GG         GG         GG         NN
#> 87         NN         GG         AG         AG         AA         GG         AA
#> 88         TT         AA         AA         AT         AT         NN         AA
#>    CLF.r2.068 CLF.r2.069 OPL.r1.126 OPL.r1.128 OPL.r1.130 OPL.r1.139 OPL.r1.146
#> 1          AG         AA         AG         AG         AA         AA         NN
#> 2          CC         CT         CC         CC         CC         CC         CC
#> 3          AA         AA         AG         AG         AA         AG         AG
#> 4          GG         GG         NN         GG         NN         NN         AG
#> 5          AC         AA         AA         AA         AA         AA         AA
#> 6          TT         TT         TT         TT         TT         NN         TT
#> 7          AG         GG         NN         AG         GG         GG         AG
#> 8          CC         CC         NN         CC         CC         CC         CC
#> 9          AG         AA         AA         AA         AA         AA         AA
#> 10         TT         CC         NN         CT         TT         CT         CT
#> 11         AA         AA         AA         AA         NN         AG         AG
#> 12         TT         CC         TT         TT         CC         TT         CT
#> 13         GT         GG         GG         GG         GG         NN         GG
#> 14         AA         AA         NN         AA         AG         AA         AG
#> 15         GG         GT         TT         GG         GG         GG         GT
#> 16         AA         AA         AA         AC         AA         AA         AA
#> 17         AC         AC         AC         AC         AA         AC         CC
#> 18         TT         GT         TT         NN         NN         GT         TT
#> 19         AC         CC         AC         AC         AC         NN         CC
#> 20         CT         CT         NN         CC         CT         CC         NN
#> 21         AA         AA         AA         AG         AA         AA         AG
#> 22         GG         GG         AG         AG         GG         GG         AG
#> 23         TT         TT         NN         TT         TT         TT         TT
#> 24         AA         AA         AC         AA         AA         AC         AC
#> 25         TT         TT         CT         NN         CC         CC         CC
#> 26         AC         AA         AA         AA         AC         AA         AA
#> 27         CG         NN         NN         GG         GG         GG         GG
#> 28         CC         CC         CG         CG         CC         CG         GG
#> 29         AA         AA         AC         AC         NN         NN         AA
#> 30         CT         NN         CT         CT         CT         CC         CC
#> 31         GG         CC         CG         CG         CC         CC         NN
#> 32         AA         AA         AA         AC         AA         NN         NN
#> 33         AA         AA         AA         AA         AA         AG         AA
#> 34         GG         GG         GG         CG         GG         GG         GG
#> 35         GG         GG         GT         GG         GG         GG         GT
#> 36         NN         AT         AT         TT         TT         TT         NN
#> 37         GG         GG         GG         CG         NN         NN         GG
#> 38         NN         AG         GG         AG         GG         AG         GG
#> 39         AG         AA         AG         GG         AA         NN         AG
#> 40         AG         AA         AA         AA         AA         AA         AA
#> 41         TT         TT         TT         TT         TT         TT         TT
#> 42         CC         AC         AC         CC         AC         NN         NN
#> 43         CC         CC         CG         CG         CC         CG         NN
#> 44         CT         CC         CT         CT         CT         TT         CT
#> 45         CT         NN         NN         CC         CC         CC         CC
#> 46         TT         TT         TT         TT         TT         NN         CT
#> 47         AA         AG         AG         AG         AA         NN         AG
#> 48         AC         AC         AC         AC         CC         AC         AA
#> 49         AG         AG         AG         AG         GG         NN         GG
#> 50         AA         AA         AA         AA         NN         AG         AG
#> 51         CC         CC         CC         CC         CT         CT         CT
#> 52         GG         GG         NN         AG         AG         NN         AG
#> 53         CC         AC         AC         AA         CC         AA         NN
#> 54         GG         GG         GG         CG         CG         GG         CG
#> 55         NN         AA         AC         AC         AA         NN         AA
#> 56         TT         AA         AT         TT         TT         TT         NN
#> 57         CT         CC         CC         TT         CT         TT         NN
#> 58         AC         AA         NN         AA         AA         AC         AA
#> 59         GG         GG         GG         GG         GG         GG         GG
#> 60         AA         AT         AA         AT         AA         AA         AT
#> 61         AG         AA         AA         AA         AG         AG         AG
#> 62         AC         CC         NN         AA         AC         AC         AC
#> 63         AA         AA         AA         AA         AA         AA         AA
#> 64         GG         AG         AA         AA         AG         AA         AG
#> 65         GT         TT         TT         TT         GT         TT         GT
#> 66         NN         GG         GT         GG         GT         GG         NN
#> 67         CC         CC         CC         TT         TT         NN         CT
#> 68         AA         AA         AA         AA         AA         AT         AT
#> 69         GT         GG         GG         GG         GT         GT         GT
#> 70         TT         TT         TT         TT         TT         TT         TT
#> 71         TT         TT         CT         TT         TT         NN         NN
#> 72         NN         TT         TT         TT         CT         TT         TT
#> 73         AT         AT         TT         TT         TT         TT         TT
#> 74         CT         CC         NN         CC         CT         CC         NN
#> 75         AG         AG         GG         AG         NN         AG         AG
#> 76         AA         AA         NN         AA         AA         AG         AG
#> 77         TT         TT         TT         TT         TT         TT         CT
#> 78         AG         AG         AA         AA         NN         AG         AA
#> 79         NN         NN         CG         CG         GG         CC         NN
#> 80         NN         CT         CT         NN         CC         TT         TT
#> 81         GG         AG         NN         GG         NN         NN         AG
#> 82         CT         TT         TT         TT         TT         TT         TT
#> 83         NN         CC         CC         CC         CC         CC         NN
#> 84         GT         GG         GG         GG         GG         GG         GG
#> 85         CC         CC         CC         AC         CC         CC         CC
#> 86         GT         TT         GG         GG         GG         GT         GG
#> 87         GG         GG         AA         AA         AG         GG         NN
#> 88         NN         AT         NN         AT         AA         NN         TT
#>    OPL.r1.147 OPL.r1.149 OPL.r2.007 OPL.r2.012 OPL.r2.022 OPL.r2.023 OPL.r2.028
#> 1          AG         GG         AA         AG         AA         NN         NN
#> 2          CC         CC         CC         CC         CC         CC         CC
#> 3          AG         AG         AG         AG         NN         AA         AG
#> 4          NN         NN         GG         AG         GG         NN         NN
#> 5          AC         AA         AA         AA         AA         AC         AA
#> 6          TT         NN         GT         TT         GT         GT         TT
#> 7          NN         AG         GG         AA         GG         GG         GG
#> 8          CC         CC         AC         AC         CC         CC         CC
#> 9          AA         AA         AA         AA         AA         AA         AG
#> 10         CT         CC         TT         TT         CT         CT         CC
#> 11         AA         AG         AG         AA         AA         NN         NN
#> 12         CT         CC         TT         CT         CT         TT         CT
#> 13         GG         GG         GG         GG         GG         GG         GG
#> 14         AA         NN         NN         AA         AA         AG         GG
#> 15         GT         GG         GT         GT         GG         GT         NN
#> 16         AA         AA         AA         AA         AA         AA         AA
#> 17         CC         CC         CC         CC         AC         CC         CC
#> 18         TT         GT         TT         TT         TT         TT         NN
#> 19         CC         AC         AC         NN         CC         AC         AC
#> 20         TT         CT         CC         CT         CT         TT         CC
#> 21         AA         AA         AA         AA         AA         AA         AA
#> 22         GG         GG         GG         GG         GG         AG         NN
#> 23         TT         CT         TT         NN         TT         CT         CT
#> 24         AC         AA         NN         AA         AA         AC         AA
#> 25         TT         CT         CC         NN         CT         NN         CC
#> 26         AA         AC         AC         AC         AA         AA         CC
#> 27         GG         GG         GG         GG         GG         GG         GG
#> 28         CC         CC         CG         CG         CG         CC         CC
#> 29         AA         AA         AA         AA         AA         AA         AA
#> 30         NN         CT         NN         NN         CC         TT         CC
#> 31         CC         CC         CG         CG         CG         CG         GG
#> 32         AA         AC         AA         AA         AA         AA         AA
#> 33         AA         AG         AA         AG         GG         AG         AA
#> 34         CG         GG         GG         GG         CG         CG         NN
#> 35         GG         GG         GG         GG         GT         GT         GG
#> 36         TT         NN         AT         TT         TT         NN         NN
#> 37         GG         GG         GG         GG         GG         GG         CG
#> 38         GG         GG         GG         NN         AG         AG         GG
#> 39         AG         AA         AA         AG         AA         GG         GG
#> 40         AA         AA         NN         AA         AG         AA         AG
#> 41         TT         TT         TT         TT         TT         TT         AT
#> 42         AC         CC         NN         AC         AC         AC         AA
#> 43         GG         CG         NN         CG         GG         CG         GG
#> 44         CT         CC         CC         CT         CT         TT         CC
#> 45         CC         CT         CC         CC         CC         CT         CC
#> 46         CT         TT         CC         NN         NN         TT         TT
#> 47         AA         AA         AG         AG         AG         AA         AA
#> 48         AC         AC         AC         AA         AA         NN         NN
#> 49         GG         GG         AG         NN         GG         GG         AG
#> 50         AA         AG         AG         AA         AA         AA         AA
#> 51         NN         CC         CC         CT         CT         NN         CC
#> 52         NN         AG         AG         GG         AG         AG         AG
#> 53         AC         AC         AC         AA         AA         AC         NN
#> 54         CG         GG         CG         CG         GG         GG         GG
#> 55         AC         NN         AA         NN         AA         NN         AC
#> 56         AT         TT         AT         NN         AT         AT         TT
#> 57         CT         TT         TT         NN         CT         CT         NN
#> 58         AC         AC         AA         AC         AC         NN         NN
#> 59         GG         GG         GG         GG         GG         GG         GG
#> 60         AT         AA         AA         AA         AT         NN         TT
#> 61         AA         AG         GG         AA         AG         AG         NN
#> 62         AA         AC         NN         AA         AC         AC         AC
#> 63         AA         AA         AA         AA         AA         AA         AA
#> 64         AA         AG         AA         AA         AG         AG         AG
#> 65         TT         NN         TT         GT         TT         TT         TT
#> 66         GT         GG         GT         NN         NN         GG         NN
#> 67         NN         CC         NN         NN         CC         TT         NN
#> 68         AA         AT         AT         NN         AA         AA         AA
#> 69         GG         GT         TT         GG         GT         GT         GG
#> 70         TT         TT         TT         TT         CT         CT         TT
#> 71         TT         CC         TT         TT         TT         CT         TT
#> 72         TT         TT         CT         TT         CT         TT         TT
#> 73         TT         TT         TT         TT         TT         TT         TT
#> 74         TT         CT         CC         CT         TT         CT         CT
#> 75         AG         GG         AG         AG         AG         GG         GG
#> 76         AA         AG         AG         AA         NN         NN         AA
#> 77         CT         TT         CT         TT         TT         CC         TT
#> 78         AA         AA         AG         NN         AG         GG         AA
#> 79         CC         CG         CC         CC         CC         CG         GG
#> 80         CT         NN         TT         CT         CT         TT         CT
#> 81         AG         AG         GG         NN         AG         NN         NN
#> 82         TT         TT         TT         TT         TT         TT         CT
#> 83         CC         CC         NN         CC         CC         NN         CC
#> 84         GG         GG         TT         GG         GG         GG         GG
#> 85         CC         CC         AC         CC         AC         AC         AC
#> 86         GG         GG         GG         GG         GG         GT         GG
#> 87         AA         AG         AG         AG         AG         AG         AA
#> 88         AA         AT         AA         AT         AA         AA         AA
#>    OPL.r2.031 OPL.r2.033 OPL.r2.079 OPL.r2.080 OPL.r2.094 PAL.r2.203 PAL.r2.237
#> 1          AG         AG         AA         GG         AG         AG         GG
#> 2          CC         CC         CC         CC         CC         CC         CC
#> 3          AG         NN         AG         AG         AA         AA         AA
#> 4          NN         NN         NN         AG         NN         AG         NN
#> 5          AA         AA         AA         AA         AA         AA         AC
#> 6          GT         TT         TT         TT         NN         TT         TT
#> 7          GG         GG         GG         GG         GG         NN         GG
#> 8          CC         AC         CC         CC         CC         CC         CC
#> 9          AA         AA         AA         AA         AA         AA         AA
#> 10         CT         CT         NN         CT         CT         TT         TT
#> 11         AA         NN         AA         AA         NN         AG         AG
#> 12         TT         TT         CT         TT         CT         CT         TT
#> 13         GG         GG         GG         GG         TT         GG         GG
#> 14         NN         AA         AG         AA         AA         GG         GG
#> 15         GG         GT         GT         GG         NN         GT         NN
#> 16         AC         AA         AA         AA         AA         CC         AA
#> 17         AC         CC         AC         AC         AC         NN         AA
#> 18         TT         TT         TT         TT         TT         TT         TT
#> 19         NN         AC         NN         AC         AA         CC         CC
#> 20         CT         CC         CT         CC         NN         CC         CC
#> 21         AA         AA         AA         GG         AA         AA         AA
#> 22         AG         GG         NN         AG         AG         GG         GG
#> 23         TT         TT         CT         CT         NN         TT         CT
#> 24         AA         AA         AC         AA         AA         AA         AA
#> 25         CT         NN         CC         CT         CC         TT         TT
#> 26         CC         NN         AA         AC         AC         AC         AA
#> 27         GG         NN         NN         GG         GG         GG         GG
#> 28         CC         CC         CG         GG         CG         CC         CC
#> 29         AA         AA         AA         AA         AA         AA         AC
#> 30         NN         CT         CC         CT         CT         CC         CT
#> 31         GG         CC         NN         NN         CC         CC         CC
#> 32         NN         AA         AA         AA         AA         AA         AA
#> 33         AA         AA         AA         AA         AG         AA         AA
#> 34         GG         CG         CG         CG         NN         GG         GG
#> 35         GG         GG         GT         GG         GG         GG         GT
#> 36         TT         TT         TT         TT         TT         TT         AT
#> 37         GG         NN         GG         GG         CG         GG         CG
#> 38         AG         AG         NN         GG         AG         AA         AA
#> 39         NN         AG         AG         NN         GG         AG         GG
#> 40         AG         AA         AG         AA         AG         AA         AA
#> 41         TT         TT         TT         TT         AT         AT         AT
#> 42         CC         NN         AC         CC         AA         CC         CC
#> 43         CG         GG         CG         CC         NN         CC         CC
#> 44         CC         CT         NN         CT         NN         CT         CC
#> 45         CC         NN         CC         CC         CT         CC         CT
#> 46         NN         CT         CT         TT         CT         NN         TT
#> 47         GG         AG         NN         GG         AA         AG         AG
#> 48         AC         NN         AC         AA         AA         AA         AA
#> 49         AG         AA         NN         GG         GG         AA         AG
#> 50         AA         NN         AA         NN         NN         AG         AG
#> 51         CC         CT         CT         CT         CC         CC         TT
#> 52         GG         GG         GG         GG         AG         GG         GG
#> 53         AC         AA         AA         AC         AC         AA         AA
#> 54         CG         CG         GG         GG         GG         GG         GG
#> 55         AA         NN         AC         AA         NN         NN         NN
#> 56         AT         TT         NN         TT         AT         AT         AA
#> 57         TT         CC         CT         TT         CC         CC         CC
#> 58         AC         NN         CC         AC         AC         CC         AC
#> 59         GG         NN         GG         GG         GG         AG         AA
#> 60         AA         AA         AA         AA         AA         TT         TT
#> 61         AG         AA         AA         GG         AA         AA         AA
#> 62         AC         AA         AA         NN         AA         AA         AA
#> 63         AA         AA         AA         AA         AA         AA         AT
#> 64         AA         AA         AG         AG         AG         GG         AG
#> 65         TT         GT         NN         TT         GT         TT         TT
#> 66         GT         GG         GT         GT         GT         GG         GG
#> 67         CC         CC         CC         CC         CC         NN         NN
#> 68         AA         AT         NN         AA         AA         AT         AT
#> 69         GT         GG         GG         TT         GG         GG         GG
#> 70         CT         TT         CT         TT         TT         TT         TT
#> 71         CT         CT         CT         CT         CT         TT         CC
#> 72         TT         TT         TT         TT         TT         TT         TT
#> 73         TT         TT         TT         TT         NN         TT         TT
#> 74         CC         NN         CT         CC         CT         NN         NN
#> 75         GG         GG         NN         AG         NN         AA         AA
#> 76         NN         AG         AA         AA         AA         AG         NN
#> 77         TT         CC         TT         CT         CT         TT         TT
#> 78         AG         NN         AG         AG         AA         GG         NN
#> 79         CC         CG         NN         CG         NN         CC         NN
#> 80         CC         CT         CT         NN         CC         CC         CC
#> 81         NN         NN         AG         AG         GG         AG         NN
#> 82         TT         TT         TT         TT         TT         NN         CT
#> 83         CC         CC         NN         CG         CC         CG         CG
#> 84         GT         GG         TT         GT         GG         GG         GG
#> 85         CC         NN         AC         CC         CC         CC         CC
#> 86         GG         GG         GG         GG         GG         GG         GG
#> 87         GG         NN         AG         AG         AG         GG         GG
#> 88         TT         AA         AA         AT         AA         NN         NN
#>    PAL.r2.238 PAL.r2.239 PAL.r2.240 PAL.r2.241 PAL.r2.242 PAL.r2.249 PAL.r2.253
#> 1          GG         GG         GG         GG         AG         GG         GG
#> 2          CC         CC         CC         CC         CC         CC         CC
#> 3          AA         AA         AA         AA         AA         NN         AA
#> 4          GG         GG         GG         GG         AG         NN         AG
#> 5          AC         NN         AC         AC         AC         AC         AA
#> 6          TT         TT         TT         TT         TT         TT         TT
#> 7          GG         GG         GG         GG         GG         GG         GG
#> 8          CC         AC         CC         CC         CC         CC         CC
#> 9          AA         AA         AA         AA         AA         AA         AA
#> 10         TT         TT         TT         TT         TT         NN         TT
#> 11         AG         AA         AG         AG         AG         NN         AA
#> 12         TT         TT         TT         NN         TT         TT         TT
#> 13         GG         GG         GG         GG         GG         GG         GG
#> 14         AG         AG         GG         GG         GG         GG         NN
#> 15         TT         NN         GT         GT         GT         GG         TT
#> 16         AA         AA         AA         AC         AA         AA         NN
#> 17         AC         AA         AA         AC         AA         AC         AC
#> 18         TT         TT         TT         TT         TT         TT         TT
#> 19         CC         CC         CC         CC         CC         NN         CC
#> 20         CC         CC         CC         CC         CC         CC         CC
#> 21         AA         AA         AA         AA         AA         AA         NN
#> 22         GG         GG         GG         GG         GG         GG         GG
#> 23         CT         CC         CC         CT         TT         CT         CC
#> 24         AA         AA         AA         AA         AA         AA         AA
#> 25         TT         TT         TT         TT         CT         TT         CT
#> 26         AC         AC         AA         AA         AA         AA         AA
#> 27         GG         GG         GG         GG         GG         GG         GG
#> 28         CC         CC         CC         CC         CC         CC         CC
#> 29         AA         AC         AC         AC         AC         AC         AA
#> 30         CC         NN         CC         CT         CC         CT         CC
#> 31         CC         CC         CC         CC         CC         CC         CC
#> 32         AA         AA         AA         AA         AA         AA         AA
#> 33         AA         AA         AA         AA         AA         AA         AA
#> 34         GG         GG         NN         GG         GG         CG         GG
#> 35         GT         GT         GT         GT         GT         GG         GT
#> 36         TT         NN         NN         TT         AT         NN         NN
#> 37         GG         CG         GG         CG         GG         NN         GG
#> 38         AA         AA         AA         AA         NN         AA         AA
#> 39         GG         GG         GG         GG         GG         NN         GG
#> 40         AA         AA         AA         AA         AA         AA         AA
#> 41         AA         AT         AT         AT         AT         TT         AT
#> 42         CC         CC         CC         CC         CC         CC         CC
#> 43         CC         CC         CC         NN         CG         CG         NN
#> 44         CT         CT         TT         CT         TT         NN         CC
#> 45         CC         CC         CC         CC         CC         CC         CC
#> 46         CT         CT         NN         NN         TT         TT         TT
#> 47         NN         AA         GG         AA         AG         NN         AG
#> 48         AA         AA         AA         AA         AA         AA         AA
#> 49         AG         AG         NN         AG         GG         AG         GG
#> 50         AG         AA         AG         AG         AG         AG         AA
#> 51         TT         CT         NN         CT         CT         CT         CT
#> 52         GG         GG         GG         GG         GG         GG         GG
#> 53         AA         AA         AA         AA         AA         AA         AA
#> 54         GG         GG         GG         GG         GG         NN         GG
#> 55         AC         CC         AC         CC         NN         NN         NN
#> 56         NN         AT         AA         AT         AA         AT         AA
#> 57         CT         CC         CC         CC         CC         CC         CT
#> 58         AC         AA         AA         CC         AC         CC         CC
#> 59         GG         NN         GG         AG         GG         AG         GG
#> 60         NN         NN         TT         AT         AT         AA         TT
#> 61         AA         AA         AA         AA         AA         AA         AA
#> 62         AA         AA         AA         NN         AA         AA         AA
#> 63         AA         AT         NN         AA         AA         AT         AA
#> 64         AG         GG         GG         AG         AA         NN         GG
#> 65         GT         GT         TT         GT         GT         TT         NN
#> 66         GG         GG         GG         GG         GG         GG         GG
#> 67         NN         TT         TT         TT         TT         NN         NN
#> 68         AT         AA         AT         AT         AT         AT         AA
#> 69         GG         GG         GG         GG         GG         GG         GG
#> 70         TT         CT         CT         TT         TT         CT         CT
#> 71         CT         NN         CT         CT         CT         TT         CT
#> 72         TT         TT         TT         TT         TT         TT         TT
#> 73         TT         TT         TT         TT         TT         TT         TT
#> 74         TT         TT         TT         TT         TT         TT         TT
#> 75         AA         AA         AA         AA         NN         AG         AA
#> 76         AG         NN         NN         AG         AG         AG         NN
#> 77         NN         TT         TT         TT         TT         TT         TT
#> 78         GG         GG         GG         GG         GG         GG         GG
#> 79         CC         NN         GG         CC         CG         GG         NN
#> 80         CC         CC         CC         CC         CC         CC         CT
#> 81         GG         AG         NN         AA         AA         NN         AG
#> 82         TT         TT         TT         CT         TT         TT         TT
#> 83         CC         NN         CG         CC         NN         NN         CG
#> 84         GG         GG         GG         GG         GG         GG         GG
#> 85         CC         CC         CC         CC         CC         CC         CC
#> 86         GG         GT         GT         GT         GG         NN         GT
#> 87         GG         GG         GG         GG         NN         GG         GG
#> 88         AT         NN         TT         AA         NN         AA         NN
#>    PAL.r2.275 PAL.r2.276 PAL.r2.278 PAL.r2.281 PAL.r2.285 PAL.r2.287 PAL.r2.290
#> 1          GG         AG         GG         AG         GG         GG         AG
#> 2          CC         CC         CC         NN         CT         CC         CC
#> 3          AA         AA         AA         AA         AA         AA         AA
#> 4          NN         AG         GG         NN         AG         NN         GG
#> 5          AC         AC         AC         AC         CC         AA         AA
#> 6          TT         TT         TT         TT         TT         TT         TT
#> 7          NN         GG         GG         GG         NN         GG         GG
#> 8          CC         CC         CC         CC         CC         CC         NN
#> 9          AA         AA         AA         AA         AA         AA         AA
#> 10         TT         TT         TT         TT         TT         TT         TT
#> 11         AG         AG         AA         AA         AA         AG         AA
#> 12         NN         TT         TT         TT         TT         TT         TT
#> 13         GG         GG         GG         GG         GG         GG         GG
#> 14         GG         AG         GG         AG         GG         GG         GG
#> 15         GG         GT         GT         GG         GT         GT         GT
#> 16         AC         AA         AA         AA         AA         AC         AC
#> 17         CC         NN         AA         NN         AA         AC         AA
#> 18         TT         TT         TT         TT         TT         TT         TT
#> 19         CC         CC         CC         CC         CC         CC         CC
#> 20         CC         CC         CC         CC         CC         CC         CC
#> 21         AA         AA         AA         AA         AA         AA         AA
#> 22         GG         NN         GG         GG         GG         GG         GG
#> 23         TT         CC         NN         CT         CT         CT         NN
#> 24         AA         AA         AA         AA         AA         AA         AA
#> 25         TT         TT         TT         NN         NN         TT         CT
#> 26         AC         AA         AA         NN         NN         AA         AC
#> 27         GG         GG         GG         GG         GG         GG         GG
#> 28         CC         CC         CC         CC         CC         CC         CC
#> 29         AC         AC         AC         AC         NN         AC         NN
#> 30         CC         CT         CT         CT         CC         CC         CT
#> 31         CC         CC         CC         CC         CC         CC         CC
#> 32         AA         AA         NN         AA         AA         AA         AA
#> 33         AA         AA         AA         AA         AA         AA         NN
#> 34         GG         GG         GG         GG         GG         GG         GG
#> 35         GT         GG         GT         GT         GG         NN         GT
#> 36         TT         TT         NN         TT         NN         TT         AT
#> 37         CG         CG         GG         NN         GG         GG         GG
#> 38         NN         NN         AA         NN         AA         AA         AA
#> 39         GG         AG         GG         GG         GG         GG         GG
#> 40         AA         AA         AA         AA         AA         AA         AA
#> 41         AA         NN         TT         TT         NN         AT         TT
#> 42         CC         CC         CC         CC         CC         CC         CC
#> 43         CG         CC         CC         CC         CC         CC         CC
#> 44         TT         CT         NN         NN         TT         TT         TT
#> 45         CC         CC         CC         NN         CC         CC         CC
#> 46         CT         TT         NN         NN         TT         NN         TT
#> 47         NN         NN         AG         AG         AG         NN         AG
#> 48         AA         AA         AA         AA         NN         AA         AA
#> 49         AG         GG         AG         AG         NN         AG         AG
#> 50         AG         AG         AA         NN         AA         AG         AA
#> 51         NN         CC         CC         CC         NN         CC         NN
#> 52         GG         GG         GG         GG         AG         GG         AG
#> 53         AA         AA         AA         AA         AA         AA         AA
#> 54         GG         GG         GG         GG         GG         GG         GG
#> 55         AC         CC         CC         NN         NN         CC         AC
#> 56         TT         AA         TT         AT         TT         NN         NN
#> 57         CT         CC         CT         CT         TT         CT         CC
#> 58         CC         CC         AC         NN         AC         AC         CC
#> 59         AG         GG         NN         AG         GG         GG         AG
#> 60         TT         NN         AT         AT         TT         AT         NN
#> 61         AA         AA         AA         AA         AA         AA         AA
#> 62         AA         AA         AA         AA         AA         AA         AA
#> 63         AA         AA         AA         AT         AA         AA         AT
#> 64         GG         NN         GG         AA         AG         AG         AG
#> 65         GT         TT         GT         GG         TT         NN         GT
#> 66         GG         GG         GG         GG         GG         GG         GG
#> 67         NN         NN         TT         NN         TT         TT         NN
#> 68         AT         AT         NN         AA         NN         AT         AA
#> 69         GG         GG         GG         GG         GG         GG         GG
#> 70         TT         TT         TT         CT         CT         CT         TT
#> 71         CT         CT         CT         CC         CT         CT         NN
#> 72         TT         TT         TT         TT         TT         TT         TT
#> 73         TT         TT         TT         TT         TT         NN         TT
#> 74         TT         TT         TT         TT         TT         TT         TT
#> 75         AA         AG         NN         AA         AA         AA         AG
#> 76         AG         AG         NN         AA         NN         AG         AA
#> 77         TT         TT         TT         TT         TT         TT         TT
#> 78         GG         GG         GG         GG         GG         GG         NN
#> 79         CC         CG         NN         CC         GG         NN         GG
#> 80         CC         CC         CT         CC         CT         CC         CC
#> 81         NN         AA         AG         NN         AG         AG         NN
#> 82         TT         NN         TT         CT         TT         TT         CT
#> 83         CG         CC         GG         CG         NN         CG         CC
#> 84         GG         GG         GG         GG         GG         GG         GG
#> 85         CC         CC         CC         CC         CC         CC         CC
#> 86         GT         GT         GG         GT         GT         TT         GT
#> 87         GG         NN         GG         GG         GG         GG         GG
#> 88         AT         AT         AT         AA         AT         AA         NN
#>    PAL.r2.291 SMR.r1.004 SMR.r1.008 SMR.r1.012 SMR.r1.013 SMR.r1.023 SMR.r1.028
#> 1          AA         GG         AA         AG         NN         AG         AA
#> 2          CC         CC         CT         CC         CC         CC         NN
#> 3          AA         AA         AA         AG         NN         AG         AA
#> 4          GG         NN         NN         AG         GG         NN         NN
#> 5          AC         AA         AC         NN         AA         AA         NN
#> 6          TT         TT         TT         NN         GT         GT         TT
#> 7          GG         AG         GG         NN         GG         AG         AG
#> 8          CC         AC         CC         NN         AC         AC         AC
#> 9          AA         AA         AA         AA         NN         AA         AG
#> 10         TT         CT         CC         CT         TT         CT         NN
#> 11         AG         AA         AG         AA         AA         AA         NN
#> 12         TT         TT         TT         TT         TT         TT         CT
#> 13         GG         GG         GG         NN         GG         GG         GT
#> 14         GG         AA         AG         AG         AA         AG         AA
#> 15         GG         GG         GT         GG         NN         GT         GT
#> 16         AA         AA         AC         AC         NN         AA         AA
#> 17         AC         CC         CC         CC         NN         CC         CC
#> 18         TT         TT         TT         TT         TT         TT         TT
#> 19         CC         AC         CC         AC         AC         CC         CC
#> 20         CC         CC         CT         CC         NN         NN         CT
#> 21         AA         AA         AG         AA         AA         AA         AA
#> 22         GG         AA         AG         AA         AG         NN         GG
#> 23         CT         NN         TT         TT         TT         CT         CT
#> 24         AA         NN         AA         AA         AC         NN         AC
#> 25         TT         TT         CC         CC         CT         CC         NN
#> 26         AA         NN         AC         AA         AC         AA         NN
#> 27         GG         GG         GG         CG         CG         GG         GG
#> 28         CC         CC         GG         CG         NN         CG         CC
#> 29         AC         AC         AC         AA         AA         AA         AC
#> 30         NN         CT         CT         CC         CC         CT         CC
#> 31         CC         CG         CC         CC         CG         NN         CC
#> 32         AA         AA         AA         AA         NN         AA         AC
#> 33         AA         AA         AG         AA         AG         AA         AA
#> 34         GG         GG         CG         GG         CG         CG         CG
#> 35         TT         GT         TT         GT         NN         NN         NN
#> 36         TT         NN         AT         TT         TT         TT         NN
#> 37         GG         GG         GG         GG         GG         CG         NN
#> 38         AA         AG         AA         NN         AG         NN         NN
#> 39         GG         AA         AG         GG         AG         GG         NN
#> 40         AA         AA         AA         AA         AA         AA         AG
#> 41         TT         TT         TT         NN         AA         NN         TT
#> 42         NN         NN         AA         CC         CC         AC         CC
#> 43         CC         CC         CG         CC         CC         CG         CG
#> 44         CC         TT         CC         CT         CT         TT         CT
#> 45         CC         CC         CC         CC         CT         NN         CC
#> 46         CT         NN         CT         CT         NN         NN         CT
#> 47         AG         AG         AA         AG         AA         AA         GG
#> 48         AA         AA         AA         AA         AA         AA         AC
#> 49         AA         GG         AA         GG         AG         GG         GG
#> 50         AG         NN         AG         AA         AA         NN         AA
#> 51         TT         CC         CC         CC         CC         NN         CC
#> 52         NN         AG         NN         NN         AG         AG         GG
#> 53         AA         NN         AC         AC         NN         AA         AC
#> 54         GG         GG         CG         GG         CG         GG         GG
#> 55         CC         AC         AA         AA         AC         AC         NN
#> 56         AT         AT         TT         AA         TT         AT         TT
#> 57         CT         CT         CT         CT         NN         CT         CC
#> 58         AC         AA         AC         AC         AC         NN         NN
#> 59         GG         GG         GG         GG         GG         GG         GG
#> 60         TT         AA         AT         AA         AT         AT         AA
#> 61         AA         AA         AA         AA         AA         AG         AA
#> 62         AA         AA         AA         NN         AC         AC         AA
#> 63         AT         AA         AA         AA         AA         AA         AT
#> 64         AG         AA         AG         AA         AA         NN         AA
#> 65         TT         TT         TT         NN         TT         GT         TT
#> 66         GG         GG         GT         GG         NN         GT         GT
#> 67         TT         NN         CC         CC         NN         CC         CC
#> 68         AT         NN         AT         AA         AA         AA         AA
#> 69         GG         GG         GG         GG         NN         GT         GG
#> 70         CT         TT         CT         TT         NN         CT         TT
#> 71         CT         TT         TT         CT         TT         CT         TT
#> 72         TT         NN         CT         TT         TT         TT         TT
#> 73         TT         TT         TT         TT         TT         TT         TT
#> 74         TT         TT         CC         CT         CT         NN         NN
#> 75         NN         GG         AG         GG         GG         NN         GG
#> 76         AG         AA         AG         AA         AA         AA         AA
#> 77         TT         TT         TT         NN         CT         NN         CT
#> 78         NN         AA         AG         AG         AA         AG         NN
#> 79         GG         NN         CC         NN         CC         CG         NN
#> 80         NN         TT         CC         CT         TT         NN         CC
#> 81         AA         NN         AG         GG         AG         NN         NN
#> 82         NN         TT         TT         NN         TT         NN         NN
#> 83         CC         CC         CC         CC         CC         CC         CC
#> 84         GG         GT         GT         GG         GG         GT         NN
#> 85         CC         CC         CC         CC         CC         AC         CC
#> 86         TT         GG         GG         GT         GG         NN         GT
#> 87         NN         GG         AG         AG         AG         AA         NN
#> 88         AA         AT         AT         AA         AT         NN         NN
#>    SMR.r1.034 SMR.r1.049 SMR.r1.052 SMR.r1.056 SMR.r1.072 SMR.r1.074 SMR.r2.077
#> 1          GG         AG         AG         GG         AA         AG         AG
#> 2          CC         CC         CC         CC         CC         CC         CC
#> 3          AG         GG         AG         AA         AA         AG         AA
#> 4          AG         AG         GG         GG         NN         NN         NN
#> 5          AA         AA         AA         AA         AA         AC         AA
#> 6          TT         GT         NN         TT         TT         TT         TT
#> 7          GG         AG         NN         AA         AG         AG         AG
#> 8          AC         CC         CC         AC         NN         AC         CC
#> 9          NN         AA         AG         AA         AA         AA         AA
#> 10         CT         CC         CT         TT         CT         CT         NN
#> 11         AG         AA         AA         NN         AA         NN         AA
#> 12         CT         CT         CT         TT         CT         TT         TT
#> 13         GG         GG         GT         GG         GT         GG         GT
#> 14         AA         AG         AG         NN         AA         AA         AA
#> 15         GG         TT         GT         NN         GT         GG         GT
#> 16         AA         AA         AA         NN         AC         AA         AA
#> 17         CC         CC         CC         CC         AC         CC         AC
#> 18         TT         TT         TT         TT         GT         TT         NN
#> 19         AC         AC         NN         CC         AC         AC         AC
#> 20         CC         CC         CC         CC         NN         CT         CC
#> 21         AA         AA         AA         AG         AA         AA         NN
#> 22         NN         NN         NN         AA         AG         GG         NN
#> 23         TT         CT         CT         CT         CT         CT         CT
#> 24         NN         AA         AA         AA         AC         AA         AA
#> 25         CT         CT         TT         NN         NN         CT         CT
#> 26         NN         AC         NN         NN         AA         AA         NN
#> 27         NN         GG         GG         CG         GG         GG         CG
#> 28         CC         GG         NN         CG         CC         CG         CG
#> 29         NN         AC         AA         AA         AA         AA         AA
#> 30         CT         CC         NN         CC         CC         TT         CC
#> 31         CG         CC         CC         GG         CC         CC         CC
#> 32         AA         AA         AA         AC         NN         AC         AA
#> 33         AA         AG         AA         AA         AG         AA         AA
#> 34         GG         CG         GG         NN         CG         GG         CG
#> 35         GG         GG         GT         GT         GT         GT         GG
#> 36         TT         TT         NN         NN         TT         NN         TT
#> 37         GG         GG         GG         GG         CG         GG         GG
#> 38         GG         AG         AG         GG         AG         AG         NN
#> 39         AG         AG         AG         AG         AA         AG         AG
#> 40         AA         AA         AA         AG         AA         AA         AG
#> 41         TT         TT         TT         TT         TT         TT         TT
#> 42         NN         AC         NN         CC         AC         NN         CC
#> 43         CC         CC         CC         CC         CG         CC         NN
#> 44         NN         CT         NN         TT         CT         CT         NN
#> 45         CC         CC         CC         NN         NN         NN         CC
#> 46         NN         NN         TT         CT         NN         CT         TT
#> 47         AA         AG         AG         AA         AA         AA         AA
#> 48         NN         AC         AA         CC         AA         AA         AA
#> 49         GG         AG         GG         GG         GG         AG         NN
#> 50         AG         AA         AA         NN         NN         AG         AA
#> 51         NN         CT         NN         CC         CC         CT         CC
#> 52         GG         AG         AG         AG         AA         AG         GG
#> 53         AA         AA         AC         AC         CC         AA         AC
#> 54         GG         GG         GG         GG         GG         GG         GG
#> 55         AC         AC         NN         NN         AA         AA         NN
#> 56         AT         TT         AT         AA         AT         AA         NN
#> 57         TT         CC         CT         CC         CC         NN         CC
#> 58         AA         NN         AC         AA         NN         AA         AC
#> 59         GG         GG         GG         GG         GG         GG         GG
#> 60         AA         NN         AA         NN         AT         AT         AA
#> 61         AA         AA         AA         NN         AA         AA         AA
#> 62         AC         AA         AA         AA         AA         AA         AA
#> 63         AA         AT         AA         AA         AA         AA         AA
#> 64         AG         AA         NN         AG         AA         AG         AA
#> 65         GT         TT         TT         TT         TT         GT         TT
#> 66         NN         GT         NN         GG         GG         NN         NN
#> 67         NN         NN         NN         CC         CT         NN         NN
#> 68         AT         AA         NN         AA         AA         AT         NN
#> 69         GG         GG         GG         GG         GG         GG         GG
#> 70         TT         TT         TT         TT         TT         NN         TT
#> 71         CT         TT         CT         CT         CT         CT         CT
#> 72         TT         TT         TT         TT         TT         NN         CT
#> 73         TT         TT         TT         TT         TT         TT         TT
#> 74         CC         NN         NN         CC         CT         CT         CC
#> 75         NN         GG         AG         AG         AG         GG         GG
#> 76         AG         AA         NN         AA         AA         AG         NN
#> 77         CT         CT         TT         TT         TT         TT         TT
#> 78         NN         AG         AG         AG         AG         NN         NN
#> 79         CC         CC         CG         NN         CG         CG         CG
#> 80         CT         CC         CC         CT         TT         CT         NN
#> 81         NN         GG         GG         NN         NN         NN         AG
#> 82         CT         TT         TT         TT         CT         TT         CT
#> 83         CC         CC         NN         NN         CC         CC         CC
#> 84         GG         GG         GG         NN         GT         GT         GT
#> 85         CC         CC         CC         CC         CC         CC         CC
#> 86         GG         NN         GG         GT         GG         GG         GT
#> 87         AG         AA         GG         AG         AA         NN         AG
#> 88         AT         AA         AA         AT         AT         NN         AT
#>    SMR.r2.086 SMR.r2.091 SMR.r2.097 SMR.r2.105 SMR.r2.139 SMR.r2.175 SMR.r2.256
#> 1          AA         AA         NN         AA         AG         AG         AG
#> 2          CC         CC         CT         CC         CC         CC         CT
#> 3          AG         AA         AA         GG         AA         AG         AG
#> 4          GG         GG         GG         GG         NN         AA         AG
#> 5          AA         NN         AA         AA         AA         AA         AA
#> 6          GT         TT         TT         TT         GT         TT         TT
#> 7          NN         GG         GG         AA         NN         GG         GG
#> 8          AC         NN         AC         NN         CC         AC         CC
#> 9          AA         AA         AA         AA         AA         AA         AA
#> 10         TT         CT         CT         CT         NN         CC         CC
#> 11         AG         AG         AG         AA         NN         NN         AG
#> 12         CC         NN         CT         TT         TT         NN         NN
#> 13         GT         NN         GG         GT         GG         GT         GG
#> 14         AA         AA         AA         AG         AG         NN         AA
#> 15         GT         NN         GG         NN         GT         GT         GG
#> 16         AC         AA         AA         AA         AA         AC         NN
#> 17         AA         AC         AC         CC         CC         CC         CC
#> 18         TT         TT         TT         TT         GT         TT         TT
#> 19         CC         AC         AC         AC         NN         NN         CC
#> 20         NN         CC         CT         CC         CT         CC         CT
#> 21         AA         AA         AA         AA         NN         NN         AA
#> 22         AG         AG         AA         AG         GG         GG         AG
#> 23         TT         CT         TT         TT         CT         TT         TT
#> 24         AA         AA         NN         AA         AA         AA         AC
#> 25         NN         CT         CT         CT         TT         CT         CT
#> 26         AA         AC         CC         AA         NN         AC         AA
#> 27         GG         GG         GG         GG         NN         CG         GG
#> 28         CG         NN         CC         CG         CG         NN         CC
#> 29         NN         AC         AA         AA         AA         AA         NN
#> 30         CC         TT         CT         TT         CC         TT         CC
#> 31         CG         CG         CC         CG         GG         CC         CC
#> 32         AA         AA         AA         AA         NN         AA         AA
#> 33         AG         AA         NN         NN         AG         AA         AA
#> 34         GG         CC         CG         CG         NN         GG         GG
#> 35         GT         GG         NN         GT         GT         GT         GT
#> 36         TT         AT         TT         TT         TT         NN         AT
#> 37         NN         GG         CG         GG         NN         NN         CG
#> 38         GG         GG         AG         GG         NN         GG         NN
#> 39         AA         AG         AG         AG         NN         AG         AA
#> 40         AA         AA         AA         AA         AA         AG         AA
#> 41         TT         TT         AT         AT         TT         AT         TT
#> 42         CC         CC         NN         CC         CC         NN         CC
#> 43         CG         CG         GG         CC         CC         CC         CC
#> 44         CC         CT         CT         TT         CT         CT         NN
#> 45         NN         CC         CC         CC         CC         CT         CC
#> 46         CT         TT         NN         CC         TT         TT         NN
#> 47         AA         AA         AA         AG         AA         AA         AA
#> 48         AA         AA         AC         AC         AA         AA         AC
#> 49         GG         GG         NN         AG         GG         NN         GG
#> 50         AG         AG         AG         AA         AA         AA         AG
#> 51         CC         CC         CT         CT         CC         CC         CT
#> 52         AG         AG         GG         AG         AG         NN         AA
#> 53         AA         AC         AC         AA         AC         AC         NN
#> 54         GG         GG         GG         GG         CG         GG         GG
#> 55         AC         NN         CC         AA         AC         AA         AA
#> 56         AT         NN         TT         AT         AT         AT         AT
#> 57         CT         CT         CC         TT         CT         CC         CC
#> 58         CC         AA         AA         AC         NN         NN         AA
#> 59         GG         GG         GG         GG         GG         GG         GG
#> 60         AT         AT         TT         AA         AT         NN         NN
#> 61         AG         AA         AA         NN         AG         AA         AA
#> 62         AC         AC         AA         AC         AC         NN         AA
#> 63         AT         AA         AT         AT         AA         AA         AA
#> 64         AG         AA         AA         AG         AA         AG         AA
#> 65         TT         TT         NN         GT         TT         TT         TT
#> 66         GG         GT         GT         GG         NN         GT         GT
#> 67         CC         CC         CC         NN         CC         CT         NN
#> 68         AT         AT         AT         AA         AA         NN         AT
#> 69         GT         GG         GG         NN         GT         GG         GG
#> 70         TT         TT         TT         TT         TT         CC         TT
#> 71         CT         CT         NN         TT         TT         TT         TT
#> 72         TT         TT         TT         TT         TT         NN         TT
#> 73         TT         TT         TT         TT         TT         TT         TT
#> 74         CT         CT         CC         CC         CC         NN         CC
#> 75         AG         AG         GG         GG         AG         AG         GG
#> 76         AG         AG         AG         AA         NN         AA         AG
#> 77         CT         TT         CT         CT         CT         CT         TT
#> 78         GG         AA         AA         NN         AG         AG         AG
#> 79         CG         GG         CC         CC         CG         CG         NN
#> 80         CT         CT         TT         TT         TT         CT         NN
#> 81         NN         AG         GG         GG         NN         NN         GG
#> 82         TT         TT         TT         CT         NN         NN         CT
#> 83         CC         CC         CC         CC         CC         CG         CC
#> 84         GG         GG         GG         GG         NN         GG         GG
#> 85         AC         CC         NN         CC         CC         CC         CC
#> 86         GT         GG         GG         GG         GT         GG         GG
#> 87         GG         AG         AA         AG         GG         NN         NN
#> 88         AT         AA         AT         AT         AA         AT         AT
#>    SMR.r2.263 SMR.r2.296 UPD.r2.120 UPD.r2.121 UPD.r2.122 UPD.r2.123 UPD.r2.124
#> 1          AA         GG         AA         AA         AG         AG         AG
#> 2          CC         CC         CT         CT         CC         CC         CC
#> 3          GG         AA         AG         GG         AG         AG         AG
#> 4          GG         AG         GG         NN         AA         AG         AG
#> 5          AA         AA         AC         AC         AC         AC         CC
#> 6          GT         GT         TT         TT         TT         TT         TT
#> 7          AG         GG         NN         GG         AG         AG         NN
#> 8          CC         CC         CC         CC         AC         CC         NN
#> 9          AA         AA         AA         AA         AA         AA         AG
#> 10         CT         TT         CT         CT         TT         TT         CT
#> 11         AA         AA         AG         AA         AG         AG         AA
#> 12         TT         TT         TT         CT         TT         TT         TT
#> 13         GT         GG         GG         GT         GT         GT         GT
#> 14         NN         AG         AG         AG         GG         NN         AG
#> 15         GT         GT         GG         GG         GT         GT         GG
#> 16         AC         AA         AC         AA         AA         AC         AA
#> 17         CC         CC         AA         CC         AA         CC         CC
#> 18         TT         TT         GT         TT         TT         TT         TT
#> 19         CC         AC         CC         AC         AC         NN         CC
#> 20         CC         CC         CC         TT         CT         CT         CT
#> 21         AA         AA         AA         AA         AA         AA         AA
#> 22         AG         AG         GG         AG         AG         AA         GG
#> 23         NN         CT         CT         CT         CT         TT         CT
#> 24         AC         AA         AA         AA         AC         AA         AA
#> 25         CC         CC         TT         CT         CT         TT         TT
#> 26         AA         AA         AA         AA         AA         AC         AA
#> 27         GG         GG         GG         GG         GG         NN         GG
#> 28         CC         CG         GG         CG         CG         CC         CC
#> 29         NN         AC         AC         AA         AA         AC         AC
#> 30         CT         CT         CC         TT         CC         TT         CC
#> 31         CG         GG         CG         CC         CG         CC         CC
#> 32         AA         AA         AC         AA         AA         AA         AA
#> 33         AG         AG         AA         AA         AA         AG         AA
#> 34         CG         CC         GG         GG         GG         GG         CG
#> 35         TT         NN         GG         GG         GT         GG         GT
#> 36         TT         TT         TT         AT         TT         TT         TT
#> 37         CG         CG         CG         CG         GG         CG         GG
#> 38         GG         AG         AA         NN         AA         NN         GG
#> 39         AG         AG         AG         AG         GG         GG         AG
#> 40         AA         AG         NN         AA         AG         AG         NN
#> 41         TT         TT         TT         AT         TT         NN         TT
#> 42         CC         AA         CC         AC         CC         NN         CC
#> 43         CG         GG         CC         CC         CG         CC         CC
#> 44         CC         NN         CT         NN         TT         CT         CT
#> 45         CC         CC         CC         CC         CT         CC         CC
#> 46         NN         TT         CC         TT         CC         CT         TT
#> 47         AA         AA         GG         GG         AA         AG         AG
#> 48         AC         AA         AA         CC         AA         AC         NN
#> 49         GG         AG         AA         AG         AG         AG         GG
#> 50         NN         AA         AG         AA         AG         AG         AA
#> 51         CC         NN         NN         TT         CT         NN         NN
#> 52         AG         AA         GG         AG         GG         AG         GG
#> 53         AA         AA         AC         AA         NN         AC         AC
#> 54         GG         GG         GG         NN         GG         GG         GG
#> 55         AC         NN         AC         AC         AC         AA         CC
#> 56         TT         TT         AA         TT         TT         AT         AT
#> 57         TT         CC         TT         CC         TT         CC         CT
#> 58         CC         AC         AC         AC         AA         AA         AC
#> 59         NN         GG         GG         GG         GG         GG         GG
#> 60         NN         AA         AA         TT         AT         AA         AT
#> 61         AG         NN         AA         AA         AA         AA         AA
#> 62         AC         AC         AA         AC         NN         AC         AA
#> 63         AT         AA         AA         AA         AA         AA         AA
#> 64         AA         AG         AG         GG         AG         AG         NN
#> 65         TT         NN         GT         TT         TT         GT         GT
#> 66         GG         GG         GT         GG         GG         GG         GG
#> 67         CT         NN         CT         CC         CC         CT         CT
#> 68         AA         AA         AT         AA         AT         AT         AA
#> 69         GT         GT         GG         GG         GG         GG         GG
#> 70         TT         TT         CT         TT         TT         TT         NN
#> 71         CT         TT         TT         CT         NN         CT         CT
#> 72         TT         TT         TT         TT         TT         CT         TT
#> 73         TT         AT         TT         TT         TT         TT         TT
#> 74         CT         CT         CT         TT         CC         CT         CT
#> 75         AA         NN         GG         GG         AG         GG         AG
#> 76         AA         NN         AG         AA         AG         AG         NN
#> 77         TT         CT         CT         TT         NN         TT         TT
#> 78         AG         NN         NN         AG         AG         AG         AG
#> 79         GG         CG         GG         CC         CC         CG         NN
#> 80         CT         CT         CT         CT         CC         TT         TT
#> 81         NN         GG         NN         GG         GG         GG         AG
#> 82         TT         CT         TT         TT         TT         CT         CT
#> 83         CC         CC         CC         CC         CC         CC         CG
#> 84         GG         GG         GG         GT         GG         GG         GG
#> 85         AC         CC         CC         CC         CC         NN         CC
#> 86         GT         GT         GT         GG         GG         GT         GT
#> 87         GG         AA         AG         GG         AA         NN         NN
#> 88         AA         AA         AA         AT         AA         AT         AT
#>    UPD.r2.129 UPD.r2.130 UPD.r2.131 UPD.r2.133 UPD.r2.135 UPD.r2.160 UPD.r2.162
#> 1          AA         AG         AG         GG         AG         NN         GG
#> 2          CC         CC         CC         CC         CC         CC         CC
#> 3          AG         AG         NN         NN         AG         AG         AG
#> 4          NN         NN         NN         AA         AG         AG         NN
#> 5          AC         NN         AA         AC         AC         AA         AC
#> 6          TT         TT         TT         TT         TT         TT         NN
#> 7          GG         NN         NN         GG         AG         AG         NN
#> 8          CC         AC         CC         AC         AC         CC         NN
#> 9          NN         AG         AA         AA         AA         AA         AA
#> 10         NN         CT         CT         CT         CT         CT         CT
#> 11         NN         AA         AG         NN         NN         NN         AG
#> 12         TT         TT         TT         CT         CT         CT         TT
#> 13         GT         GG         GG         GT         GG         TT         GT
#> 14         AG         AG         NN         AG         AG         AG         AG
#> 15         GT         NN         GG         NN         GG         GG         GG
#> 16         NN         AA         AA         AA         AC         NN         AA
#> 17         CC         AA         CC         AC         AC         CC         AC
#> 18         TT         TT         TT         GT         TT         TT         TT
#> 19         NN         AC         CC         AC         CC         AC         NN
#> 20         CT         CC         NN         CC         CT         CT         CT
#> 21         AA         AA         AG         AG         AA         NN         AG
#> 22         AG         AG         GG         AA         GG         AG         AG
#> 23         TT         CT         CT         TT         CT         NN         NN
#> 24         AC         AA         AA         AA         AA         AA         AA
#> 25         TT         CT         CT         NN         TT         TT         CT
#> 26         NN         AA         AA         NN         AA         AA         AA
#> 27         GG         CG         CG         NN         GG         NN         GG
#> 28         CC         CG         CG         GG         CC         GG         CG
#> 29         AA         AA         AC         AC         AC         AA         NN
#> 30         CT         CC         NN         CT         CC         CC         CC
#> 31         NN         CC         CC         CC         CC         NN         CG
#> 32         AA         AA         AA         AA         AA         AA         NN
#> 33         AA         AA         AG         AA         NN         AA         NN
#> 34         CG         GG         CG         CG         CG         CG         NN
#> 35         GG         GG         TT         TT         GG         GT         GG
#> 36         NN         NN         AT         NN         TT         AT         TT
#> 37         NN         CG         CG         NN         CC         GG         GG
#> 38         NN         AG         AA         AA         AG         NN         NN
#> 39         GG         AG         AA         GG         NN         AG         AG
#> 40         AG         AA         AG         AA         AG         AG         AA
#> 41         AT         AT         NN         TT         NN         TT         TT
#> 42         CC         AC         AA         AC         CC         AC         CC
#> 43         CC         CC         CG         CC         CC         CC         GG
#> 44         CT         TT         CT         CT         CT         TT         CT
#> 45         CC         CT         CC         CC         CC         NN         CC
#> 46         CT         CT         CT         TT         CT         NN         CT
#> 47         NN         AG         AG         AA         AA         AG         AG
#> 48         NN         AA         AC         AC         AC         AC         NN
#> 49         GG         AG         AA         GG         NN         NN         AA
#> 50         AG         AA         AG         AA         AA         AA         AG
#> 51         CT         CT         CC         CC         CC         NN         CC
#> 52         GG         GG         AG         AG         AG         AA         AG
#> 53         AC         AC         AA         AA         AA         AC         AC
#> 54         GG         CG         CG         CG         GG         GG         GG
#> 55         AA         AC         AC         AA         NN         AA         AA
#> 56         TT         TT         AT         AT         AT         AT         TT
#> 57         CC         CT         NN         CC         CC         CC         CT
#> 58         NN         AA         AA         AA         CC         CC         AC
#> 59         GG         GG         GG         NN         NN         AG         GG
#> 60         AT         AA         AA         AT         AT         NN         TT
#> 61         AA         AA         AA         AA         AA         AA         AA
#> 62         AA         AC         AC         AA         AA         AC         NN
#> 63         AA         AA         AA         AT         AA         AA         AA
#> 64         AG         GG         NN         GG         AG         GG         GG
#> 65         TT         GG         TT         TT         GT         TT         TT
#> 66         GG         NN         GG         GG         GG         GT         GG
#> 67         CT         CT         CC         CC         TT         CT         TT
#> 68         NN         AA         AT         AA         AA         AA         AT
#> 69         GG         GG         GG         GG         GG         GG         GG
#> 70         TT         TT         TT         NN         NN         TT         CT
#> 71         CT         CT         TT         CT         TT         NN         CT
#> 72         TT         NN         CT         TT         CT         CT         TT
#> 73         TT         TT         AT         TT         TT         TT         TT
#> 74         CT         CT         CT         TT         CT         NN         CC
#> 75         AA         AG         AG         AA         AG         AG         AG
#> 76         AG         AA         AG         AA         NN         AA         AG
#> 77         CT         CT         CC         TT         CT         CT         CT
#> 78         AA         AA         AG         AA         NN         AG         NN
#> 79         GG         CG         NN         CG         NN         CG         GG
#> 80         CT         CC         CC         CC         CC         CT         NN
#> 81         NN         NN         NN         AG         GG         AG         GG
#> 82         TT         TT         CT         TT         TT         CT         NN
#> 83         CC         CG         CC         CC         NN         CG         CC
#> 84         GT         GT         GT         GG         GG         NN         GG
#> 85         CC         CC         CC         NN         AC         CC         CC
#> 86         NN         GG         GG         GT         NN         GT         GG
#> 87         GG         AG         AG         GG         NN         NN         AA
#> 88         AT         TT         TT         TT         NN         TT         TT
#>    UPD.r2.174 UPD.r2.189 UPD.r2.201
#> 1          NN         GG         AG
#> 2          CC         CC         CC
#> 3          AA         AA         AA
#> 4          GG         NN         GG
#> 5          AA         AA         AA
#> 6          TT         TT         TT
#> 7          AG         GG         GG
#> 8          CC         AC         AA
#> 9          AG         AA         AG
#> 10         CT         CT         TT
#> 11         AA         AG         NN
#> 12         TT         NN         NN
#> 13         GG         GT         GT
#> 14         AG         AA         AA
#> 15         TT         NN         GG
#> 16         AC         AA         AA
#> 17         AA         AC         AC
#> 18         TT         TT         TT
#> 19         AC         AA         CC
#> 20         CT         CC         CC
#> 21         AA         AA         AA
#> 22         AA         AG         GG
#> 23         NN         TT         TT
#> 24         NN         AA         AA
#> 25         CT         CT         CT
#> 26         AA         AC         AC
#> 27         CG         GG         CG
#> 28         CG         CC         CG
#> 29         AA         NN         AC
#> 30         CT         NN         TT
#> 31         NN         CC         NN
#> 32         AA         AA         AC
#> 33         AA         AA         AA
#> 34         CG         NN         GG
#> 35         GG         GG         NN
#> 36         NN         TT         NN
#> 37         NN         CC         NN
#> 38         GG         AG         AG
#> 39         GG         AG         AG
#> 40         AA         AA         AA
#> 41         AT         AT         TT
#> 42         AA         NN         AC
#> 43         CC         CC         CG
#> 44         NN         CC         TT
#> 45         CC         CC         CC
#> 46         CT         TT         TT
#> 47         AA         AA         AA
#> 48         AA         AC         CC
#> 49         AG         AG         GG
#> 50         AA         AG         AA
#> 51         CC         CT         CC
#> 52         AA         NN         NN
#> 53         AC         AA         AC
#> 54         GG         GG         GG
#> 55         AA         NN         AA
#> 56         AT         AT         AT
#> 57         CT         CT         TT
#> 58         CC         NN         AC
#> 59         GG         GG         GG
#> 60         AA         AA         AT
#> 61         AA         AA         AA
#> 62         NN         AC         AC
#> 63         AT         AA         AA
#> 64         GG         AG         AG
#> 65         TT         NN         GT
#> 66         GT         NN         GG
#> 67         CC         NN         TT
#> 68         AA         AT         AA
#> 69         GG         GG         GG
#> 70         CT         TT         TT
#> 71         TT         TT         TT
#> 72         TT         NN         CT
#> 73         TT         NN         TT
#> 74         CT         NN         CC
#> 75         GG         AG         AG
#> 76         AA         AG         AA
#> 77         CT         CC         CT
#> 78         AG         AA         AG
#> 79         GG         GG         CC
#> 80         CC         CT         CT
#> 81         NN         AG         GG
#> 82         CT         CT         CT
#> 83         CC         CC         CC
#> 84         NN         NN         GG
#> 85         CC         AC         CC
#> 86         GT         GT         TT
#> 87         AA         AG         AG
#> 88         AT         AT         AT

# The same filters, but with minor allele frequency considered per-population
# and a full re-run of loci filters after individual removal.
filter_snps(stickSNPs, maf = 0.05, hf_hets = 0.55, min_ind = 0.5, 
            min_loci = 0.75, re_run = "full", maf_facets = "pop")
#> Adding missing facets...
#> Initializing...
#> Filtering loci. Starting loci: 100 
#> Filtering non-biallelic loci...
#>  0 bad loci
#> Filtering non_polymorphic loci...
#>  0 bad loci
#> Filtering loci sequenced in few individuals...
#>  0 bad loci
#> Filtering low minor allele frequencies by facet.
#>  0 bad loci
#> Filtering high frequency heterozygote loci...
#>  2 bad loci
#>  Ending loci: 98 
#> Filtering individuals. Starting individuals: 100 
#> Filtering out individuals sequenced in few kept loci...
#> Re-calculating and adding facets.
#>  Ending individuals: 88 
#> Re-filtering loci...
#> Filtering non-biallelic loci...
#>  0 bad loci
#> Filtering non_polymorphic loci...
#>  0 bad loci
#> Filtering loci sequenced in few individuals...
#>  0 bad loci
#> Filtering low minor allele frequencies by facet.
#>  0 bad loci
#> Filtering high frequency heterozygote loci...
#>  0 bad loci
#>  Final loci count: 98 
#>    ASP.r2.141 ASP.r2.144 ASP.r2.146 ASP.r2.196 ASP.r2.216 ASP.r2.221 ASP.r2.227
#> 1          AG         GG         AA         NN         AG         NN         GG
#> 2          GG         GG         GG         GT         GG         GG         NN
#> 3          CC         CC         CC         CC         CC         CT         CC
#> 4          AA         AA         AA         AG         AA         AA         AA
#> 5          GG         NN         NN         GG         AG         AG         AG
#> 6          AC         AC         AA         AA         AC         AA         AC
#> 7          TT         TT         TT         GT         TT         TT         TT
#> 8          GG         GG         AG         AG         GG         GG         AG
#> 9          CC         AC         CC         CC         CC         CC         CC
#> 10         AA         AA         AG         AA         AA         AA         AA
#> 11         NN         TT         TT         TT         CT         CC         CC
#> 12         AG         AG         AG         AA         AA         AG         AA
#> 13         TT         TT         CT         CT         NN         TT         TT
#> 14         GG         GT         GT         GG         GG         GG         GG
#> 15         AA         AG         AG         AA         AA         AA         AA
#> 16         NN         GG         GT         GG         GT         GT         GG
#> 17         AA         AA         AA         AA         NN         AC         AC
#> 18         AC         NN         CC         AC         AC         CC         AC
#> 19         TT         TT         TT         TT         TT         TT         TT
#> 20         GG         AG         GG         NN         GG         GG         GG
#> 21         NN         NN         CC         AC         AC         CC         CC
#> 22         CT         NN         TT         CT         NN         CC         TT
#> 23         GG         GG         GG         GG         GG         GG         AG
#> 24         AA         AA         AA         AA         AA         AA         AA
#> 25         NN         AG         AG         GG         AG         AG         GG
#> 26         CT         CT         CT         NN         TT         NN         TT
#> 27         AA         AA         AA         AA         AA         AA         AA
#> 28         TT         NN         NN         CT         CT         TT         TT
#> 29         AA         AC         AA         AA         AA         AA         NN
#> 30         GG         GG         CG         GG         GG         GG         GG
#> 31         CG         CG         CC         CC         CC         CG         CC
#> 32         AA         AA         AA         AA         AA         AA         AC
#> 33         CC         CT         CC         CT         CT         CT         CC
#> 34         CG         CG         CG         CC         NN         CC         CG
#> 35         AA         AA         AA         AA         NN         AC         AA
#> 36         AA         AA         AA         AA         AG         AA         AG
#> 37         CG         GG         GG         GG         GG         NN         GG
#> 38         GG         TT         GG         GT         TT         GG         GT
#> 39         TT         NN         TT         AT         TT         TT         TT
#> 40         CG         CC         CG         CG         GG         CG         GG
#> 41         AG         AA         AG         AA         AA         AG         AG
#> 42         AG         NN         NN         GG         GG         AG         AG
#> 43         AG         AG         AA         AA         AA         AA         AA
#> 44         AT         TT         TT         NN         TT         TT         NN
#> 45         GG         GG         GG         GG         GG         GG         GG
#> 46         NN         CC         AA         NN         CC         NN         AA
#> 47         CC         CC         CC         CG         CC         CG         CG
#> 48         CT         TT         TT         TT         TT         CT         CC
#> 49         GG         GG         GG         NN         GG         GG         GG
#> 50         NN         CC         CT         CC         CT         CC         CC
#> 51         TT         NN         CT         TT         CT         CT         NN
#> 52         AA         AG         AG         AG         AG         AA         AG
#> 53         AA         AA         CC         AA         AC         AA         AA
#> 54         GG         GG         GG         GG         GG         GG         GG
#> 55         GG         NN         AG         AG         AG         AG         AG
#> 56         AG         AG         AG         AA         NN         AG         AA
#> 57         GG         NN         GG         GG         GG         GG         AG
#> 58         TT         CT         CC         TT         CC         NN         CT
#> 59         AA         AG         GG         GG         AA         GG         AG
#> 60         AA         AA         AC         AC         AC         AC         NN
#> 61         TT         NN         TT         TT         TT         TT         NN
#> 62         CG         GG         GG         GG         GG         GG         GG
#> 63         AC         NN         NN         AA         AC         AC         AC
#> 64         TT         AT         AT         AT         AT         AA         AT
#> 65         CT         CC         CT         CC         CT         CT         CC
#> 66         AC         AC         NN         AC         CC         CC         AC
#> 67         GG         GG         GG         GG         GG         GG         GG
#> 68         AA         AT         AT         AA         NN         AT         TT
#> 69         AA         AA         NN         AA         AA         AA         AA
#> 70         AC         AC         AA         AC         NN         AC         AC
#> 71         AA         AA         AA         AA         AA         AA         AA
#> 72         GG         GG         GG         AG         NN         AG         AG
#> 73         TT         GT         TT         GT         TT         GT         NN
#> 74         GT         GG         GT         GG         GG         GT         GG
#> 75         CC         CT         CT         NN         NN         TT         NN
#> 76         AT         AT         AT         AA         AA         AT         NN
#> 77         GG         GG         GG         GG         GG         GG         GG
#> 78         TT         TT         CT         CT         CT         CT         TT
#> 79         CT         CT         TT         TT         CT         TT         TT
#> 80         CT         TT         CT         TT         TT         NN         TT
#> 81         AT         TT         AT         TT         AT         TT         TT
#> 82         CC         CC         CT         CC         NN         CT         CC
#> 83         AG         AG         GG         AG         NN         AA         AA
#> 84         AG         AG         AG         AA         AA         AG         NN
#> 85         CT         CC         CT         TT         CT         CC         CT
#> 86         AG         NN         AG         GG         AG         AA         AG
#> 87         NN         CG         CG         CC         CG         CG         NN
#> 88         CT         CT         CC         CC         CT         CC         CC
#> 89         CC         CC         NN         CC         CT         CC         CT
#> 90         AA         NN         GG         GG         AA         AG         AG
#> 91         TT         TT         TT         CT         TT         TT         TT
#> 92         CC         CC         CC         CC         NN         CC         CC
#> 93         GG         GG         GG         GT         GG         GT         GT
#> 94         CC         NN         CC         CC         CC         CC         CC
#> 95         GG         GG         GG         GG         GG         GG         GG
#> 96         GG         GG         GG         GG         GT         GG         GT
#> 97         AA         GG         AA         AG         AG         AG         AG
#> 98         AT         NN         AT         TT         AT         TT         TT
#>    ASP.r2.230 ASP.r2.232 ASP.r2.235 CLF.r2.040 CLF.r2.048 CLF.r2.061 CLF.r2.064
#> 1          AA         GG         AG         NN         AG         AG         AA
#> 2          GG         GG         GG         GG         GG         GG         GG
#> 3          CC         CT         CC         CC         CC         CC         CC
#> 4          AA         AG         AG         AA         AG         AA         NN
#> 5          NN         GG         NN         GG         NN         NN         NN
#> 6          AA         AC         NN         AA         AC         AC         AA
#> 7          TT         TT         TT         GG         TT         TT         TT
#> 8          GG         AG         GG         GG         AG         AG         AA
#> 9          AC         CC         AA         CC         AC         CC         AC
#> 10         AG         AA         AA         AA         AA         AA         AA
#> 11         CT         CT         CT         CT         CT         CT         CT
#> 12         AG         AA         AA         AA         AA         AG         AA
#> 13         NN         TT         TT         CT         CT         TT         CT
#> 14         GT         GG         GT         GG         GG         GG         GT
#> 15         AG         AG         AG         AA         AA         AA         AA
#> 16         GG         GG         GG         GT         GT         GG         NN
#> 17         AA         AC         AA         AA         AA         AA         AA
#> 18         AC         CC         AC         CC         AC         AC         NN
#> 19         TT         TT         TT         GT         TT         TT         GT
#> 20         GG         NN         GG         GG         GG         GG         GG
#> 21         AC         CC         CC         NN         AC         AC         NN
#> 22         CT         CC         CT         CC         CC         NN         TT
#> 23         GG         GG         NN         GG         GG         GG         GG
#> 24         AA         AA         AG         AA         AA         AA         AA
#> 25         AG         AG         GG         GG         GG         AG         NN
#> 26         CT         TT         CT         TT         TT         TT         TT
#> 27         NN         AC         AA         AA         AC         AA         AA
#> 28         TT         TT         TT         CT         TT         CT         TT
#> 29         AA         AA         AA         AA         NN         NN         AA
#> 30         GG         GG         CG         GG         GG         CG         GG
#> 31         CG         CG         CG         CC         CG         CC         CC
#> 32         AC         AC         AA         AC         AC         NN         AA
#> 33         CC         CT         CT         CT         CC         CT         NN
#> 34         CG         CG         NN         CC         CG         GG         CG
#> 35         AA         AA         AA         AA         NN         AA         AA
#> 36         AA         AG         NN         AG         AA         AA         AG
#> 37         CC         GG         CG         GG         CG         CG         NN
#> 38         NN         GG         TT         GG         GG         GG         GG
#> 39         TT         AT         AT         TT         NN         NN         AT
#> 40         CC         GG         CG         CG         GG         NN         CG
#> 41         GG         AG         AG         AG         GG         GG         AG
#> 42         AG         GG         GG         NN         GG         AG         GG
#> 43         AA         AA         AA         AG         NN         AA         AA
#> 44         TT         AA         TT         TT         TT         TT         TT
#> 45         GG         GG         GG         GG         GG         GG         GG
#> 46         CC         CC         AC         AC         CC         CC         NN
#> 47         CC         CG         CC         CC         CC         CC         CC
#> 48         NN         TT         CT         CT         CT         TT         CT
#> 49         GG         GG         GG         GG         GG         CG         GG
#> 50         CT         CC         CC         CC         CC         CC         CC
#> 51         CT         CC         TT         TT         NN         TT         TT
#> 52         AG         AG         AG         GG         AG         AA         AA
#> 53         AC         AA         AA         AC         NN         AA         NN
#> 54         GG         GG         GG         GG         GG         GG         GT
#> 55         AG         AA         NN         GG         GG         GG         AG
#> 56         AG         AA         NN         AA         AA         AG         AA
#> 57         NN         GG         GG         NN         GG         GG         GG
#> 58         CC         CT         CT         CC         CC         CC         CC
#> 59         GG         AA         AG         AG         AG         AG         AG
#> 60         NN         CC         NN         AC         AC         CC         CC
#> 61         TT         TT         TT         TT         NN         NN         CT
#> 62         CG         GG         GG         GG         GG         GG         GG
#> 63         CC         AA         AA         CC         AC         AC         AA
#> 64         TT         TT         TT         TT         TT         TT         TT
#> 65         CC         CC         CC         TT         CC         CT         CT
#> 66         AA         AA         AC         AC         AC         CC         AA
#> 67         NN         GG         GG         GG         GG         GG         GG
#> 68         AA         AA         AA         AT         AA         AA         AT
#> 69         AA         AA         AA         NN         AA         AA         AA
#> 70         AA         AA         AC         CC         AA         AC         NN
#> 71         AA         AA         AA         AA         AA         AA         AA
#> 72         AA         GG         AA         AG         AA         AG         GG
#> 73         GT         GT         GT         GT         TT         GT         TT
#> 74         GG         GG         NN         GG         GT         GG         GG
#> 75         NN         CC         NN         CT         NN         NN         NN
#> 76         AT         NN         AA         AA         AA         AT         AA
#> 77         NN         GG         GG         TT         GG         GG         GG
#> 78         TT         TT         TT         TT         TT         TT         TT
#> 79         CC         TT         TT         CT         CT         TT         TT
#> 80         TT         TT         TT         TT         NN         NN         TT
#> 81         TT         AT         TT         AT         AT         TT         AT
#> 82         CT         CC         TT         TT         CT         CC         CC
#> 83         AG         AA         AG         AG         AA         NN         AG
#> 84         AG         NN         NN         AA         AA         AG         AA
#> 85         TT         TT         CC         TT         CT         TT         TT
#> 86         AG         AG         NN         AA         AA         AG         AA
#> 87         CC         CG         CC         CC         NN         GG         CC
#> 88         CC         CC         CC         CC         CC         CC         CC
#> 89         CC         NN         CC         CC         NN         TT         CC
#> 90         NN         NN         AG         AG         NN         GG         GG
#> 91         TT         CT         CT         CT         TT         CT         TT
#> 92         CC         NN         NN         CC         CC         CC         CC
#> 93         GG         GT         GG         GT         TT         GG         GG
#> 94         AC         CC         NN         CC         CC         CC         AC
#> 95         GG         GG         GG         GG         GG         GG         GG
#> 96         NN         GG         TT         GG         GG         GG         NN
#> 97         NN         GG         AG         AG         AA         GG         AA
#> 98         TT         AA         AA         AT         AT         NN         AA
#>    CLF.r2.068 CLF.r2.069 OPL.r1.125 OPL.r1.126 OPL.r1.128 OPL.r1.130 OPL.r1.139
#> 1          AG         AA         NN         AG         AG         AA         AA
#> 2          GG         GG         GG         GG         GG         GT         GG
#> 3          CC         CT         NN         CC         CC         CC         CC
#> 4          AA         AA         AG         AG         AG         AA         AG
#> 5          GG         GG         GG         NN         GG         NN         NN
#> 6          AC         AA         AA         AA         AA         AA         AA
#> 7          TT         TT         TT         TT         TT         TT         NN
#> 8          AG         GG         AG         NN         AG         GG         GG
#> 9          CC         CC         CC         NN         CC         CC         CC
#> 10         AG         AA         NN         AA         AA         AA         AA
#> 11         TT         CC         CC         NN         CT         TT         CT
#> 12         AA         AA         AA         AA         AA         NN         AG
#> 13         TT         CC         TT         TT         TT         CC         TT
#> 14         GT         GG         NN         GG         GG         GG         NN
#> 15         AA         AA         NN         NN         AA         AG         AA
#> 16         GG         GT         GT         TT         GG         GG         GG
#> 17         AA         AA         AA         AA         AC         AA         AA
#> 18         AC         AC         CC         AC         AC         AA         AC
#> 19         TT         GT         TT         TT         NN         NN         GT
#> 20         GG         AG         GG         NN         AG         GG         GG
#> 21         AC         CC         AC         AC         AC         AC         NN
#> 22         CT         CT         NN         NN         CC         CT         CC
#> 23         GG         GG         GG         GG         GG         GG         GG
#> 24         AA         AA         NN         AA         AG         AA         AA
#> 25         GG         GG         NN         AG         AG         GG         GG
#> 26         TT         TT         TT         NN         TT         TT         TT
#> 27         AA         AA         NN         AC         AA         AA         AC
#> 28         TT         TT         TT         CT         NN         CC         CC
#> 29         AC         AA         NN         AA         AA         AC         AA
#> 30         CG         NN         CG         NN         GG         GG         GG
#> 31         CC         CC         NN         CG         CG         CC         CG
#> 32         AA         AA         AC         AC         AC         NN         NN
#> 33         CT         NN         CC         CT         CT         CT         CC
#> 34         GG         CC         NN         CG         CG         CC         CC
#> 35         AA         AA         AA         AA         AC         AA         NN
#> 36         AA         AA         AA         AA         AA         AA         AG
#> 37         GG         GG         NN         GG         CG         GG         GG
#> 38         GG         GG         GT         GT         GG         GG         GG
#> 39         NN         AT         NN         AT         TT         TT         TT
#> 40         GG         GG         CG         GG         CG         NN         NN
#> 41         NN         AG         GG         GG         AG         GG         AG
#> 42         AG         AA         GG         AG         GG         AA         NN
#> 43         AG         AA         AG         AA         AA         AA         AA
#> 44         TT         TT         TT         TT         TT         TT         TT
#> 45         GG         GG         GG         GG         GG         NN         GG
#> 46         CC         AC         NN         AC         CC         AC         NN
#> 47         CC         CC         CC         CG         CG         CC         CG
#> 48         CT         CC         TT         CT         CT         CT         TT
#> 49         CG         GG         GG         GG         GG         GG         GG
#> 50         CT         NN         CC         NN         CC         CC         CC
#> 51         TT         TT         CT         TT         TT         TT         NN
#> 52         AA         AG         AA         AG         AG         AA         NN
#> 53         AC         AC         AA         AC         AC         CC         AC
#> 54         GG         GG         GG         GG         GG         GG         GG
#> 55         AG         AG         AG         AG         AG         GG         NN
#> 56         AA         AA         AA         AA         AA         NN         AG
#> 57         GG         GG         GG         GG         GG         AG         AG
#> 58         CC         CC         CC         CC         CC         CT         CT
#> 59         GG         GG         AG         NN         AG         AG         NN
#> 60         CC         AC         AA         AC         AA         CC         AA
#> 61         TT         NN         NN         NN         TT         TT         CT
#> 62         GG         GG         NN         GG         CG         CG         GG
#> 63         NN         AA         AA         AC         AC         AA         NN
#> 64         TT         AA         AT         AT         TT         TT         TT
#> 65         CT         CC         CT         CC         TT         CT         TT
#> 66         AC         AA         NN         NN         AA         AA         AC
#> 67         GG         GG         GG         GG         GG         GG         GG
#> 68         AA         AT         AT         AA         AT         AA         AA
#> 69         AG         AA         AA         AA         AA         AG         AG
#> 70         AC         CC         AA         NN         AA         AC         AC
#> 71         AA         AA         AA         AA         AA         AA         AA
#> 72         GG         AG         NN         AA         AA         AG         AA
#> 73         GT         TT         TT         TT         TT         GT         TT
#> 74         NN         GG         NN         GT         GG         GT         GG
#> 75         CC         CC         NN         CC         TT         TT         NN
#> 76         AA         AA         AA         AA         AA         AA         AT
#> 77         GT         GG         GG         GG         GG         GT         GT
#> 78         TT         TT         TT         TT         TT         TT         TT
#> 79         TT         TT         TT         CT         TT         TT         NN
#> 80         NN         TT         TT         TT         TT         CT         TT
#> 81         AT         AT         TT         TT         TT         TT         TT
#> 82         CT         CC         CT         NN         CC         CT         CC
#> 83         AG         AG         AG         GG         AG         NN         AG
#> 84         AA         AA         AA         NN         AA         AA         AG
#> 85         TT         TT         TT         TT         TT         TT         TT
#> 86         AG         AG         AA         AA         AA         NN         AG
#> 87         NN         NN         NN         CG         CG         GG         CC
#> 88         CC         CC         CC         CC         CC         CC         CC
#> 89         NN         CT         CT         CT         NN         CC         TT
#> 90         GG         AG         NN         NN         GG         NN         NN
#> 91         CT         TT         TT         TT         TT         TT         TT
#> 92         NN         CC         CC         CC         CC         CC         CC
#> 93         GT         GG         GG         GG         GG         GG         GG
#> 94         CC         CC         CC         CC         AC         CC         CC
#> 95         AG         GG         GG         GG         GG         GG         GG
#> 96         GT         TT         NN         GG         GG         GG         GT
#> 97         GG         GG         AG         AA         AA         AG         GG
#> 98         NN         AT         AA         NN         AT         AA         NN
#>    OPL.r1.146 OPL.r1.147 OPL.r1.149 OPL.r2.007 OPL.r2.012 OPL.r2.022 OPL.r2.023
#> 1          NN         AG         GG         AA         AG         AA         NN
#> 2          NN         GG         NN         GG         GG         GG         NN
#> 3          CC         CC         CC         CC         CC         CC         CC
#> 4          AG         AG         AG         AG         AG         NN         AA
#> 5          AG         NN         NN         GG         AG         GG         NN
#> 6          AA         AC         AA         AA         AA         AA         AC
#> 7          TT         TT         NN         GT         TT         GT         GT
#> 8          AG         NN         AG         GG         AA         GG         GG
#> 9          CC         CC         CC         AC         AC         CC         CC
#> 10         AA         AA         AA         AA         AA         AA         AA
#> 11         CT         CT         CC         TT         TT         CT         CT
#> 12         AG         AA         AG         AG         AA         AA         NN
#> 13         CT         CT         CC         TT         CT         CT         TT
#> 14         GG         GG         GG         GG         GG         GG         GG
#> 15         AG         AA         NN         NN         AA         AA         AG
#> 16         GT         GT         GG         GT         GT         GG         GT
#> 17         AA         AA         AA         AA         AA         AA         AA
#> 18         CC         CC         CC         CC         CC         AC         CC
#> 19         TT         TT         GT         TT         TT         TT         TT
#> 20         GG         GG         GG         GG         GG         GG         GG
#> 21         CC         CC         AC         AC         NN         CC         AC
#> 22         NN         TT         CT         CC         CT         CT         TT
#> 23         GG         GG         GG         GG         GG         GG         GG
#> 24         AG         AA         AA         AA         AA         AA         AA
#> 25         AG         GG         GG         GG         GG         GG         AG
#> 26         TT         TT         CT         TT         NN         TT         CT
#> 27         AC         AC         AA         NN         AA         AA         AC
#> 28         CC         TT         CT         CC         NN         CT         NN
#> 29         AA         AA         AC         AC         AC         AA         AA
#> 30         GG         GG         GG         GG         GG         GG         GG
#> 31         GG         CC         CC         CG         CG         CG         CC
#> 32         AA         AA         AA         AA         AA         AA         AA
#> 33         CC         NN         CT         NN         NN         CC         TT
#> 34         NN         CC         CC         CG         CG         CG         CG
#> 35         NN         AA         AC         AA         AA         AA         AA
#> 36         AA         AA         AG         AA         AG         GG         AG
#> 37         GG         CG         GG         GG         GG         CG         CG
#> 38         GT         GG         GG         GG         GG         GT         GT
#> 39         NN         TT         NN         AT         TT         TT         NN
#> 40         GG         GG         GG         GG         GG         GG         GG
#> 41         GG         GG         GG         GG         NN         AG         AG
#> 42         AG         AG         AA         AA         AG         AA         GG
#> 43         AA         AA         AA         NN         AA         AG         AA
#> 44         TT         TT         TT         TT         TT         TT         TT
#> 45         GG         NN         GG         GG         NN         GG         GG
#> 46         NN         AC         CC         NN         AC         AC         AC
#> 47         NN         GG         CG         NN         CG         GG         CG
#> 48         CT         CT         CC         CC         CT         CT         TT
#> 49         GG         GG         GG         GG         GG         GG         GG
#> 50         CC         CC         CT         CC         CC         CC         CT
#> 51         CT         CT         TT         CC         NN         NN         TT
#> 52         AG         AA         AA         AG         AG         AG         AA
#> 53         AA         AC         AC         AC         AA         AA         NN
#> 54         GG         GG         GG         GG         GG         GT         GG
#> 55         GG         GG         GG         AG         NN         GG         GG
#> 56         AG         AA         AG         AG         AA         AA         AA
#> 57         AG         GG         GG         GG         GG         GG         GG
#> 58         CT         NN         CC         CC         CT         CT         NN
#> 59         AG         NN         AG         AG         GG         AG         AG
#> 60         NN         AC         AC         AC         AA         AA         AC
#> 61         TT         TT         TT         NN         TT         NN         TT
#> 62         CG         CG         GG         CG         CG         GG         GG
#> 63         AA         AC         NN         AA         NN         AA         NN
#> 64         NN         AT         TT         AT         NN         AT         AT
#> 65         NN         CT         TT         TT         NN         CT         CT
#> 66         AA         AC         AC         AA         AC         AC         NN
#> 67         GG         GG         GG         GG         GG         GG         GG
#> 68         AT         AT         AA         AA         AA         AT         NN
#> 69         AG         AA         AG         GG         AA         AG         AG
#> 70         AC         AA         AC         NN         AA         AC         AC
#> 71         AA         AA         AA         AA         AA         AA         AA
#> 72         AG         AA         AG         AA         AA         AG         AG
#> 73         GT         TT         NN         TT         GT         TT         TT
#> 74         NN         GT         GG         GT         NN         NN         GG
#> 75         CT         NN         CC         NN         NN         CC         TT
#> 76         AT         AA         AT         AT         NN         AA         AA
#> 77         GT         GG         GT         TT         GG         GT         GT
#> 78         TT         TT         TT         TT         TT         CT         CT
#> 79         NN         TT         CC         TT         TT         TT         CT
#> 80         TT         TT         TT         CT         TT         CT         TT
#> 81         TT         TT         TT         TT         TT         TT         TT
#> 82         NN         TT         CT         CC         CT         TT         CT
#> 83         AG         AG         GG         AG         AG         AG         GG
#> 84         AG         AA         AG         AG         AA         NN         NN
#> 85         CT         CT         TT         CT         TT         TT         CC
#> 86         AA         AA         AA         AG         NN         AG         GG
#> 87         NN         CC         CG         CC         CC         CC         CG
#> 88         CC         CC         CC         CC         CC         CC         CC
#> 89         TT         CT         NN         TT         CT         CT         TT
#> 90         AG         AG         AG         GG         NN         AG         NN
#> 91         TT         TT         TT         TT         TT         TT         TT
#> 92         NN         CC         CC         NN         CC         CC         NN
#> 93         GG         GG         GG         TT         GG         GG         GG
#> 94         CC         CC         CC         AC         CC         AC         AC
#> 95         GG         AG         NN         GG         GG         GG         GG
#> 96         GG         GG         GG         GG         GG         GG         GT
#> 97         NN         AA         AG         AG         AG         AG         AG
#> 98         TT         AA         AT         AA         AT         AA         AA
#>    OPL.r2.028 OPL.r2.031 OPL.r2.033 OPL.r2.079 OPL.r2.080 OPL.r2.094 PAL.r2.203
#> 1          NN         AG         AG         AA         GG         AG         AG
#> 2          GG         GG         GG         GG         GG         GG         GG
#> 3          CC         CC         CC         CC         CC         CC         CC
#> 4          AG         AG         NN         AG         AG         AA         AA
#> 5          NN         NN         NN         NN         AG         NN         AG
#> 6          AA         AA         AA         AA         AA         AA         AA
#> 7          TT         GT         TT         TT         TT         NN         TT
#> 8          GG         GG         GG         GG         GG         GG         NN
#> 9          CC         CC         AC         CC         CC         CC         CC
#> 10         AG         AA         AA         AA         AA         AA         AA
#> 11         CC         CT         CT         NN         CT         CT         TT
#> 12         NN         AA         NN         AA         AA         NN         AG
#> 13         CT         TT         TT         CT         TT         CT         CT
#> 14         GG         GG         GG         GG         GG         TT         GG
#> 15         GG         NN         AA         AG         AA         AA         GG
#> 16         NN         GG         GT         GT         GG         NN         GT
#> 17         AA         AC         AA         AA         AA         AA         CC
#> 18         CC         AC         CC         AC         AC         AC         NN
#> 19         NN         TT         TT         TT         TT         TT         TT
#> 20         GG         GG         GG         GG         GG         GG         NN
#> 21         AC         NN         AC         NN         AC         AA         CC
#> 22         CC         CT         CC         CT         CC         NN         CC
#> 23         GG         GG         GG         NN         GG         NN         GG
#> 24         AA         AA         AA         AA         GG         AA         AA
#> 25         NN         AG         GG         NN         AG         AG         GG
#> 26         CT         TT         TT         CT         CT         NN         TT
#> 27         AA         AA         AA         AC         AA         AA         AA
#> 28         CC         CT         NN         CC         CT         CC         TT
#> 29         CC         CC         NN         AA         AC         AC         AC
#> 30         GG         GG         NN         NN         GG         GG         GG
#> 31         CC         CC         CC         CG         GG         CG         CC
#> 32         AA         AA         AA         AA         AA         AA         AA
#> 33         CC         NN         CT         CC         CT         CT         CC
#> 34         GG         GG         CC         NN         NN         CC         CC
#> 35         AA         NN         AA         AA         AA         AA         AA
#> 36         AA         AA         AA         AA         AA         AG         AA
#> 37         NN         GG         CG         CG         CG         NN         GG
#> 38         GG         GG         GG         GT         GG         GG         GG
#> 39         NN         TT         TT         TT         TT         TT         TT
#> 40         CG         GG         NN         GG         GG         CG         GG
#> 41         GG         AG         AG         NN         GG         AG         AA
#> 42         GG         NN         AG         AG         NN         GG         AG
#> 43         AG         AG         AA         AG         AA         AG         AA
#> 44         AT         TT         TT         TT         TT         AT         AT
#> 45         NN         GG         GG         NN         GG         GG         GG
#> 46         AA         CC         NN         AC         CC         AA         CC
#> 47         GG         CG         GG         CG         CC         NN         CC
#> 48         CC         CC         CT         NN         CT         NN         CT
#> 49         GG         GG         GG         GG         GG         GG         GG
#> 50         CC         CC         NN         CC         CC         CT         CC
#> 51         TT         NN         CT         CT         TT         CT         NN
#> 52         AA         GG         AG         NN         GG         AA         AG
#> 53         NN         AC         NN         AC         AA         AA         AA
#> 54         GG         GG         GG         GG         GG         GG         GG
#> 55         AG         AG         AA         NN         GG         GG         AA
#> 56         AA         AA         NN         AA         NN         NN         AG
#> 57         NN         GG         NN         NN         GG         GG         NN
#> 58         CC         CC         CT         CT         CT         CC         CC
#> 59         AG         GG         GG         GG         GG         AG         GG
#> 60         NN         AC         AA         AA         AC         AC         AA
#> 61         NN         NN         TT         TT         TT         NN         NN
#> 62         GG         CG         CG         GG         GG         GG         GG
#> 63         AC         AA         NN         AC         AA         NN         NN
#> 64         TT         AT         TT         NN         TT         AT         AT
#> 65         NN         TT         CC         CT         TT         CC         CC
#> 66         NN         AC         NN         CC         AC         AC         CC
#> 67         GG         GG         NN         GG         GG         GG         AG
#> 68         TT         AA         AA         AA         AA         AA         TT
#> 69         NN         AG         AA         AA         GG         AA         AA
#> 70         AC         AC         AA         AA         NN         AA         AA
#> 71         AA         AA         AA         AA         AA         AA         AA
#> 72         AG         AA         AA         AG         AG         AG         GG
#> 73         TT         TT         GT         NN         TT         GT         TT
#> 74         NN         GT         GG         GT         GT         GT         GG
#> 75         NN         CC         CC         CC         CC         CC         NN
#> 76         AA         AA         AT         NN         AA         AA         AT
#> 77         GG         GT         GG         GG         TT         GG         GG
#> 78         TT         CT         TT         CT         TT         TT         TT
#> 79         TT         CT         CT         CT         CT         CT         TT
#> 80         TT         TT         TT         TT         TT         TT         TT
#> 81         TT         TT         TT         TT         TT         NN         TT
#> 82         CT         CC         NN         CT         CC         CT         NN
#> 83         GG         GG         GG         NN         AG         NN         AA
#> 84         AA         NN         AG         AA         AA         AA         AG
#> 85         TT         TT         CC         TT         CT         CT         TT
#> 86         AA         AG         NN         AG         AG         AA         GG
#> 87         GG         CC         CG         NN         CG         NN         CC
#> 88         CC         CC         CC         CC         CC         CC         CC
#> 89         CT         CC         CT         CT         NN         CC         CC
#> 90         NN         NN         NN         AG         AG         GG         AG
#> 91         CT         TT         TT         TT         TT         TT         NN
#> 92         CC         CC         CC         NN         CG         CC         CG
#> 93         GG         GT         GG         TT         GT         GG         GG
#> 94         AC         CC         NN         AC         CC         CC         CC
#> 95         NN         GG         GG         GG         GG         GG         GG
#> 96         GG         GG         GG         GG         GG         GG         GG
#> 97         AA         GG         NN         AG         AG         AG         GG
#> 98         AA         TT         AA         AA         AT         AA         NN
#>    PAL.r2.237 PAL.r2.238 PAL.r2.239 PAL.r2.240 PAL.r2.241 PAL.r2.242 PAL.r2.249
#> 1          GG         GG         GG         GG         GG         AG         GG
#> 2          GG         GG         GG         GG         GG         GG         GG
#> 3          CC         CC         CC         CC         CC         CC         CC
#> 4          AA         AA         AA         AA         AA         AA         NN
#> 5          NN         GG         GG         GG         GG         AG         NN
#> 6          AC         AC         NN         AC         AC         AC         AC
#> 7          TT         TT         TT         TT         TT         TT         TT
#> 8          GG         GG         GG         GG         GG         GG         GG
#> 9          CC         CC         AC         CC         CC         CC         CC
#> 10         AA         AA         AA         AA         AA         AA         AA
#> 11         TT         TT         TT         TT         TT         TT         NN
#> 12         AG         AG         AA         AG         AG         AG         NN
#> 13         TT         TT         TT         TT         NN         TT         TT
#> 14         GG         GG         GG         GG         GG         GG         GG
#> 15         GG         AG         AG         GG         GG         GG         GG
#> 16         NN         TT         NN         GT         GT         GT         GG
#> 17         AA         AA         AA         AA         AC         AA         AA
#> 18         AA         AC         AA         AA         AC         AA         AC
#> 19         TT         TT         TT         TT         TT         TT         TT
#> 20         GG         GG         GG         GG         GG         GG         GG
#> 21         CC         CC         CC         CC         CC         CC         NN
#> 22         CC         CC         CC         CC         CC         CC         CC
#> 23         GG         GG         GG         GG         GG         NN         GG
#> 24         AA         AA         AA         AA         AA         AA         AA
#> 25         GG         GG         GG         GG         GG         GG         GG
#> 26         CT         CT         CC         CC         CT         TT         CT
#> 27         AA         AA         AA         AA         AA         AA         AA
#> 28         TT         TT         TT         TT         TT         CT         TT
#> 29         AA         AC         AC         AA         AA         AA         AA
#> 30         GG         GG         GG         GG         GG         GG         GG
#> 31         CC         CC         CC         CC         CC         CC         CC
#> 32         AC         AA         AC         AC         AC         AC         AC
#> 33         CT         CC         NN         CC         CT         CC         CT
#> 34         CC         CC         CC         CC         CC         CC         CC
#> 35         AA         AA         AA         AA         AA         AA         AA
#> 36         AA         AA         AA         AA         AA         AA         AA
#> 37         GG         GG         GG         NN         GG         GG         CG
#> 38         GT         GT         GT         GT         GT         GT         GG
#> 39         AT         TT         NN         NN         TT         AT         NN
#> 40         CG         GG         CG         GG         CG         GG         NN
#> 41         AA         AA         AA         AA         AA         NN         AA
#> 42         GG         GG         GG         GG         GG         GG         NN
#> 43         AA         AA         AA         AA         AA         AA         AA
#> 44         AT         AA         AT         AT         AT         AT         TT
#> 45         GG         GG         GG         GG         GG         GG         GG
#> 46         CC         CC         CC         CC         CC         CC         CC
#> 47         CC         CC         CC         CC         NN         CG         CG
#> 48         CC         CT         CT         TT         CT         TT         NN
#> 49         GG         GG         GG         GG         GG         GG         GG
#> 50         CT         CC         CC         CC         CC         CC         CC
#> 51         TT         CT         CT         NN         NN         TT         TT
#> 52         AG         NN         AA         GG         AA         AG         NN
#> 53         AA         AA         AA         AA         AA         AA         AA
#> 54         GG         GG         GG         GG         GG         GG         GG
#> 55         AG         AG         AG         NN         AG         GG         AG
#> 56         AG         AG         AA         AG         AG         AG         AG
#> 57         NN         GG         GG         GG         GG         NN         GG
#> 58         TT         TT         CT         NN         CT         CT         CT
#> 59         GG         GG         GG         GG         GG         GG         GG
#> 60         AA         AA         AA         AA         AA         AA         AA
#> 61         TT         TT         NN         NN         TT         TT         TT
#> 62         GG         GG         GG         GG         GG         GG         NN
#> 63         NN         AC         CC         AC         CC         NN         NN
#> 64         AA         NN         AT         AA         AT         AA         AT
#> 65         CC         CT         CC         CC         CC         CC         CC
#> 66         AC         AC         AA         AA         CC         AC         CC
#> 67         AA         GG         NN         GG         AG         GG         AG
#> 68         TT         NN         NN         TT         AT         AT         AA
#> 69         AA         AA         AA         AA         AA         AA         AA
#> 70         AA         AA         AA         AA         NN         AA         AA
#> 71         AT         AA         AT         NN         AA         AA         AT
#> 72         AG         AG         GG         GG         AG         AA         NN
#> 73         TT         GT         GT         TT         GT         GT         TT
#> 74         GG         GG         GG         GG         GG         GG         GG
#> 75         NN         NN         TT         TT         TT         TT         NN
#> 76         AT         AT         AA         AT         AT         AT         AT
#> 77         GG         GG         GG         GG         GG         GG         GG
#> 78         TT         TT         CT         CT         TT         TT         CT
#> 79         CC         CT         NN         CT         CT         CT         TT
#> 80         TT         TT         TT         TT         TT         TT         TT
#> 81         TT         TT         TT         TT         TT         TT         TT
#> 82         NN         TT         TT         TT         TT         TT         TT
#> 83         AA         AA         AA         AA         AA         NN         AG
#> 84         NN         AG         NN         NN         AG         AG         AG
#> 85         TT         NN         TT         TT         TT         TT         TT
#> 86         NN         GG         GG         GG         GG         GG         GG
#> 87         NN         CC         NN         GG         CC         CG         GG
#> 88         CC         CC         CC         CC         CC         CC         CC
#> 89         CC         CC         CC         CC         CC         CC         CC
#> 90         NN         GG         AG         NN         AA         AA         NN
#> 91         CT         TT         TT         TT         CT         TT         TT
#> 92         CG         CC         NN         CG         CC         NN         NN
#> 93         GG         GG         GG         GG         GG         GG         GG
#> 94         CC         CC         CC         CC         CC         CC         CC
#> 95         GG         GG         GG         GG         GG         GG         GG
#> 96         GG         GG         GT         GT         GT         GG         NN
#> 97         GG         GG         GG         GG         GG         NN         GG
#> 98         NN         AT         NN         TT         AA         NN         AA
#>    PAL.r2.253 PAL.r2.275 PAL.r2.276 PAL.r2.278 PAL.r2.281 PAL.r2.285 PAL.r2.287
#> 1          GG         GG         AG         GG         AG         GG         GG
#> 2          GG         GG         GG         GG         GG         GG         GG
#> 3          CC         CC         CC         CC         NN         CT         CC
#> 4          AA         AA         AA         AA         AA         AA         AA
#> 5          AG         NN         AG         GG         NN         AG         NN
#> 6          AA         AC         AC         AC         AC         CC         AA
#> 7          TT         TT         TT         TT         TT         TT         TT
#> 8          GG         NN         GG         GG         GG         NN         GG
#> 9          CC         CC         CC         CC         CC         CC         CC
#> 10         AA         AA         AA         AA         AA         AA         AA
#> 11         TT         TT         TT         TT         TT         TT         TT
#> 12         AA         AG         AG         AA         AA         AA         AG
#> 13         TT         NN         TT         TT         TT         TT         TT
#> 14         GG         GG         GG         GG         GG         GG         GG
#> 15         NN         GG         AG         GG         AG         GG         GG
#> 16         TT         GG         GT         GT         GG         GT         GT
#> 17         NN         AC         AA         AA         AA         AA         AC
#> 18         AC         CC         NN         AA         NN         AA         AC
#> 19         TT         TT         TT         TT         TT         TT         TT
#> 20         GG         GG         GG         NN         GG         GG         GG
#> 21         CC         CC         CC         CC         CC         CC         CC
#> 22         CC         CC         CC         CC         CC         CC         CC
#> 23         GG         GG         GG         GG         GG         GG         GG
#> 24         NN         AA         AA         AA         AA         AA         AA
#> 25         GG         GG         NN         GG         GG         GG         GG
#> 26         CC         TT         CC         NN         CT         CT         CT
#> 27         AA         AA         AA         AA         AA         AA         AA
#> 28         CT         TT         TT         TT         NN         NN         TT
#> 29         AA         AC         AA         AA         NN         NN         AA
#> 30         GG         GG         GG         GG         GG         GG         GG
#> 31         CC         CC         CC         CC         CC         CC         CC
#> 32         AA         AC         AC         AC         AC         NN         AC
#> 33         CC         CC         CT         CT         CT         CC         CC
#> 34         CC         CC         CC         CC         CC         CC         CC
#> 35         AA         AA         AA         NN         AA         AA         AA
#> 36         AA         AA         AA         AA         AA         AA         AA
#> 37         GG         GG         GG         GG         GG         GG         GG
#> 38         GT         GT         GG         GT         GT         GG         NN
#> 39         NN         TT         TT         NN         TT         NN         TT
#> 40         GG         CG         CG         GG         NN         GG         GG
#> 41         AA         NN         NN         AA         NN         AA         AA
#> 42         GG         GG         AG         GG         GG         GG         GG
#> 43         AA         AA         AA         AA         AA         AA         AA
#> 44         AT         AA         NN         TT         TT         NN         AT
#> 45         GG         GG         GG         GG         GG         GG         GG
#> 46         CC         CC         CC         CC         CC         CC         CC
#> 47         NN         CG         CC         CC         CC         CC         CC
#> 48         CC         TT         CT         NN         NN         TT         TT
#> 49         GG         GG         GG         GG         GG         GG         GG
#> 50         CC         CC         CC         CC         NN         CC         CC
#> 51         TT         CT         TT         NN         NN         TT         NN
#> 52         AG         NN         NN         AG         AG         AG         NN
#> 53         AA         AA         AA         AA         AA         NN         AA
#> 54         GG         GG         GG         GG         GG         GG         GG
#> 55         GG         AG         GG         AG         AG         NN         AG
#> 56         AA         AG         AG         AA         NN         AA         AG
#> 57         GG         GG         GG         GG         GG         GG         NN
#> 58         CT         NN         CC         CC         CC         NN         CC
#> 59         GG         GG         GG         GG         GG         AG         GG
#> 60         AA         AA         AA         AA         AA         AA         AA
#> 61         NN         TT         TT         TT         TT         TT         TT
#> 62         GG         GG         GG         GG         GG         GG         GG
#> 63         NN         AC         CC         CC         NN         NN         CC
#> 64         AA         TT         AA         TT         AT         TT         NN
#> 65         CT         CT         CC         CT         CT         TT         CT
#> 66         CC         CC         CC         AC         NN         AC         AC
#> 67         GG         AG         GG         NN         AG         GG         GG
#> 68         TT         TT         NN         AT         AT         TT         AT
#> 69         AA         AA         AA         AA         AA         AA         AA
#> 70         AA         AA         AA         AA         AA         AA         AA
#> 71         AA         AA         AA         AA         AT         AA         AA
#> 72         GG         GG         NN         GG         AA         AG         AG
#> 73         NN         GT         TT         GT         GG         TT         NN
#> 74         GG         GG         GG         GG         GG         GG         GG
#> 75         NN         NN         NN         TT         NN         TT         TT
#> 76         AA         AT         AT         NN         AA         NN         AT
#> 77         GG         GG         GG         GG         GG         GG         GG
#> 78         CT         TT         TT         TT         CT         CT         CT
#> 79         CT         CT         CT         CT         CC         CT         CT
#> 80         TT         TT         TT         TT         TT         TT         TT
#> 81         TT         TT         TT         TT         TT         TT         NN
#> 82         TT         TT         TT         TT         TT         TT         TT
#> 83         AA         AA         AG         NN         AA         AA         AA
#> 84         NN         AG         AG         NN         AA         NN         AG
#> 85         TT         TT         TT         TT         TT         TT         TT
#> 86         GG         GG         GG         GG         GG         GG         GG
#> 87         NN         CC         CG         NN         CC         GG         NN
#> 88         CC         CC         CC         CC         CC         CC         CC
#> 89         CT         CC         CC         CT         CC         CT         CC
#> 90         AG         NN         AA         AG         NN         AG         AG
#> 91         TT         TT         NN         TT         CT         TT         TT
#> 92         CG         CG         CC         GG         CG         NN         CG
#> 93         GG         GG         GG         GG         GG         GG         GG
#> 94         CC         CC         CC         CC         CC         CC         CC
#> 95         GG         GG         GG         GG         GG         GG         GG
#> 96         GT         GT         GT         GG         GT         GT         TT
#> 97         GG         GG         NN         GG         GG         GG         GG
#> 98         NN         AT         AT         AT         AA         AT         AA
#>    PAL.r2.290 PAL.r2.291 SMR.r1.004 SMR.r1.008 SMR.r1.012 SMR.r1.013 SMR.r1.023
#> 1          AG         AA         GG         AA         AG         NN         AG
#> 2          GG         GG         GT         GG         GG         GG         GG
#> 3          CC         CC         CC         CT         CC         CC         CC
#> 4          AA         AA         AA         AA         AG         NN         AG
#> 5          GG         GG         NN         NN         AG         GG         NN
#> 6          AA         AC         AA         AC         NN         AA         AA
#> 7          TT         TT         TT         TT         NN         GT         GT
#> 8          GG         GG         AG         GG         NN         GG         AG
#> 9          NN         CC         AC         CC         NN         AC         AC
#> 10         AA         AA         AA         AA         AA         NN         AA
#> 11         TT         TT         CT         CC         CT         TT         CT
#> 12         AA         AG         AA         AG         AA         AA         AA
#> 13         TT         TT         TT         TT         TT         TT         TT
#> 14         GG         GG         GG         GG         NN         GG         GG
#> 15         GG         GG         AA         AG         AG         AA         AG
#> 16         GT         GG         GG         GT         GG         NN         GT
#> 17         AC         AA         AA         AC         AC         NN         AA
#> 18         AA         AC         CC         CC         CC         NN         CC
#> 19         TT         TT         TT         TT         TT         TT         TT
#> 20         GG         GG         GG         AG         AG         GG         GG
#> 21         CC         CC         AC         CC         AC         AC         CC
#> 22         CC         CC         CC         CT         CC         NN         NN
#> 23         GG         GG         AG         GG         GG         GG         AG
#> 24         AA         AA         AA         AG         AA         AA         AA
#> 25         GG         GG         AA         AG         AA         AG         NN
#> 26         NN         CT         NN         TT         TT         TT         CT
#> 27         AA         AA         NN         AA         AA         AC         NN
#> 28         CT         TT         TT         CC         CC         CT         CC
#> 29         AC         AA         NN         AC         AA         AC         AA
#> 30         GG         GG         GG         GG         CG         CG         GG
#> 31         CC         CC         CC         GG         CG         NN         CG
#> 32         NN         AC         AC         AC         AA         AA         AA
#> 33         CT         NN         CT         CT         CC         CC         CT
#> 34         CC         CC         CG         CC         CC         CG         NN
#> 35         AA         AA         AA         AA         AA         NN         AA
#> 36         NN         AA         AA         AG         AA         AG         AA
#> 37         GG         GG         GG         CG         GG         CG         CG
#> 38         GT         TT         GT         TT         GT         NN         NN
#> 39         AT         TT         NN         AT         TT         TT         TT
#> 40         GG         GG         GG         GG         GG         GG         CG
#> 41         AA         AA         AG         AA         NN         AG         NN
#> 42         GG         GG         AA         AG         GG         AG         GG
#> 43         AA         AA         AA         AA         AA         AA         AA
#> 44         TT         TT         TT         TT         NN         AA         NN
#> 45         GG         GG         GG         GG         GG         GG         GG
#> 46         CC         NN         NN         AA         CC         CC         AC
#> 47         CC         CC         CC         CG         CC         CC         CG
#> 48         TT         CC         TT         CC         CT         CT         TT
#> 49         GG         GG         GG         GG         GG         GG         GG
#> 50         CC         CC         CC         CC         CC         CT         NN
#> 51         TT         CT         NN         CT         CT         NN         NN
#> 52         AG         AG         AG         AA         AG         AA         AA
#> 53         AA         AA         AA         AA         AA         AA         AA
#> 54         NN         GG         GG         GG         GG         GG         GG
#> 55         AG         AA         GG         AA         GG         AG         GG
#> 56         AA         AG         NN         AG         AA         AA         NN
#> 57         NN         GG         GG         GG         GG         GG         GG
#> 58         NN         TT         CC         CC         CC         CC         NN
#> 59         AG         NN         AG         NN         NN         AG         AG
#> 60         AA         AA         NN         AC         AC         NN         AA
#> 61         TT         TT         TT         TT         TT         TT         TT
#> 62         GG         GG         GG         CG         GG         CG         GG
#> 63         AC         CC         AC         AA         AA         AC         AC
#> 64         NN         AT         AT         TT         AA         TT         AT
#> 65         CC         CT         CT         CT         CT         NN         CT
#> 66         CC         AC         AA         AC         AC         AC         NN
#> 67         AG         GG         GG         GG         GG         GG         GG
#> 68         NN         TT         AA         AT         AA         AT         AT
#> 69         AA         AA         AA         AA         AA         AA         AG
#> 70         AA         AA         AA         AA         NN         AC         AC
#> 71         AT         AT         AA         AA         AA         AA         AA
#> 72         AG         AG         AA         AG         AA         AA         NN
#> 73         GT         TT         TT         TT         NN         TT         GT
#> 74         GG         GG         GG         GT         GG         NN         GT
#> 75         NN         TT         NN         CC         CC         NN         CC
#> 76         AA         AT         NN         AT         AA         AA         AA
#> 77         GG         GG         GG         GG         GG         NN         GT
#> 78         TT         CT         TT         CT         TT         NN         CT
#> 79         NN         CT         TT         TT         CT         TT         CT
#> 80         TT         TT         NN         CT         TT         TT         TT
#> 81         TT         TT         TT         TT         TT         TT         TT
#> 82         TT         TT         TT         CC         CT         CT         NN
#> 83         AG         NN         GG         AG         GG         GG         NN
#> 84         AA         AG         AA         AG         AA         AA         AA
#> 85         TT         TT         TT         TT         NN         CT         NN
#> 86         NN         NN         AA         AG         AG         AA         AG
#> 87         GG         GG         NN         CC         NN         CC         CG
#> 88         CC         CC         CC         CC         CC         CC         CC
#> 89         CC         NN         TT         CC         CT         TT         NN
#> 90         NN         AA         NN         AG         GG         AG         NN
#> 91         CT         NN         TT         TT         NN         TT         NN
#> 92         CC         CC         CC         CC         CC         CC         CC
#> 93         GG         GG         GT         GT         GG         GG         GT
#> 94         CC         CC         CC         CC         CC         CC         AC
#> 95         GG         GG         GG         GG         GG         GG         GG
#> 96         GT         TT         GG         GG         GT         GG         NN
#> 97         GG         NN         GG         AG         AG         AG         AA
#> 98         NN         AA         AT         AT         AA         AT         NN
#>    SMR.r1.028 SMR.r1.034 SMR.r1.049 SMR.r1.052 SMR.r1.056 SMR.r1.072 SMR.r1.074
#> 1          AA         GG         AG         AG         GG         AA         AG
#> 2          GG         GG         NN         GG         GG         GG         GG
#> 3          NN         CC         CC         CC         CC         CC         CC
#> 4          AA         AG         GG         AG         AA         AA         AG
#> 5          NN         AG         AG         GG         GG         NN         NN
#> 6          NN         AA         AA         AA         AA         AA         AC
#> 7          TT         TT         GT         NN         TT         TT         TT
#> 8          AG         GG         AG         NN         AA         AG         AG
#> 9          AC         AC         CC         CC         AC         NN         AC
#> 10         AG         NN         AA         AG         AA         AA         AA
#> 11         NN         CT         CC         CT         TT         CT         CT
#> 12         NN         AG         AA         AA         NN         AA         NN
#> 13         CT         CT         CT         CT         TT         CT         TT
#> 14         GT         GG         GG         GT         GG         GT         GG
#> 15         AA         AA         AG         AG         NN         AA         AA
#> 16         GT         GG         TT         GT         NN         GT         GG
#> 17         AA         AA         AA         AA         NN         AC         AA
#> 18         CC         CC         CC         CC         CC         AC         CC
#> 19         TT         TT         TT         TT         TT         GT         TT
#> 20         GG         GG         GG         GG         AG         GG         GG
#> 21         CC         AC         AC         NN         CC         AC         AC
#> 22         CT         CC         CC         CC         CC         NN         CT
#> 23         GG         GG         NN         AG         GG         NN         GG
#> 24         AA         AA         AA         AA         AG         AA         AA
#> 25         GG         NN         NN         NN         AA         AG         GG
#> 26         CT         TT         CT         CT         CT         CT         CT
#> 27         AC         NN         AA         AA         AA         AC         AA
#> 28         NN         CT         CT         TT         NN         NN         CT
#> 29         NN         NN         AC         NN         NN         AA         AA
#> 30         GG         NN         GG         GG         CG         GG         GG
#> 31         CC         CC         GG         NN         CG         CC         CG
#> 32         AC         NN         AC         AA         AA         AA         AA
#> 33         CC         CT         CC         NN         CC         CC         TT
#> 34         CC         CG         CC         CC         GG         CC         CC
#> 35         AC         AA         AA         AA         AC         NN         AC
#> 36         AA         AA         AG         AA         AA         AG         AA
#> 37         CG         GG         CG         GG         NN         CG         GG
#> 38         NN         GG         GG         GT         GT         GT         GT
#> 39         NN         TT         TT         NN         NN         TT         NN
#> 40         NN         GG         GG         GG         GG         CG         GG
#> 41         NN         GG         AG         AG         GG         AG         AG
#> 42         NN         AG         AG         AG         AG         AA         AG
#> 43         AG         AA         AA         AA         AG         AA         AA
#> 44         TT         TT         TT         TT         TT         TT         TT
#> 45         GG         GG         GG         GG         GG         GG         GG
#> 46         CC         NN         AC         NN         CC         AC         NN
#> 47         CG         CC         CC         CC         CC         CG         CC
#> 48         CT         NN         CT         NN         TT         CT         CT
#> 49         GG         GG         GG         GG         GG         GG         GG
#> 50         CC         CC         CC         CC         NN         NN         NN
#> 51         CT         NN         NN         TT         CT         NN         CT
#> 52         GG         AA         AG         AG         AA         AA         AA
#> 53         AC         NN         AC         AA         CC         AA         AA
#> 54         GG         GG         GG         GG         GG         GG         GG
#> 55         GG         GG         AG         GG         GG         GG         AG
#> 56         AA         AG         AA         AA         NN         NN         AG
#> 57         NN         GG         NN         GG         NN         GG         GG
#> 58         CC         NN         CT         NN         CC         CC         CT
#> 59         GG         GG         AG         AG         AG         AA         AG
#> 60         AC         AA         AA         AC         AC         CC         AA
#> 61         NN         NN         TT         CT         TT         TT         TT
#> 62         GG         GG         GG         GG         GG         GG         GG
#> 63         NN         AC         AC         NN         NN         AA         AA
#> 64         TT         AT         TT         AT         AA         AT         AA
#> 65         CC         TT         CC         CT         CC         CC         NN
#> 66         NN         AA         NN         AC         AA         NN         AA
#> 67         GG         GG         GG         GG         GG         GG         GG
#> 68         AA         AA         NN         AA         NN         AT         AT
#> 69         AA         AA         AA         AA         NN         AA         AA
#> 70         AA         AC         AA         AA         AA         AA         AA
#> 71         AT         AA         AT         AA         AA         AA         AA
#> 72         AA         AG         AA         NN         AG         AA         AG
#> 73         TT         GT         TT         TT         TT         TT         GT
#> 74         GT         NN         GT         NN         GG         GG         NN
#> 75         CC         NN         NN         NN         CC         CT         NN
#> 76         AA         AT         AA         NN         AA         AA         AT
#> 77         GG         GG         GG         GG         GG         GG         GG
#> 78         TT         TT         TT         TT         TT         TT         NN
#> 79         TT         CT         TT         CT         CT         CT         CT
#> 80         TT         TT         TT         TT         TT         TT         NN
#> 81         TT         TT         TT         TT         TT         TT         TT
#> 82         NN         CC         NN         NN         CC         CT         CT
#> 83         GG         NN         GG         AG         AG         AG         GG
#> 84         AA         AG         AA         NN         AA         AA         AG
#> 85         CT         CT         CT         TT         TT         TT         TT
#> 86         NN         NN         AG         AG         AG         AG         NN
#> 87         NN         CC         CC         CG         NN         CG         CG
#> 88         CC         CC         CC         CC         CC         CC         CC
#> 89         CC         CT         CC         CC         CT         TT         CT
#> 90         NN         NN         GG         GG         NN         NN         NN
#> 91         NN         CT         TT         TT         TT         CT         TT
#> 92         CC         CC         CC         NN         NN         CC         CC
#> 93         NN         GG         GG         GG         NN         GT         GT
#> 94         CC         CC         CC         CC         CC         CC         CC
#> 95         GG         GG         GG         GG         NN         GG         NN
#> 96         GT         GG         NN         GG         GT         GG         GG
#> 97         NN         AG         AA         GG         AG         AA         NN
#> 98         NN         AT         AA         AA         AT         AT         NN
#>    SMR.r2.077 SMR.r2.086 SMR.r2.091 SMR.r2.097 SMR.r2.105 SMR.r2.139 SMR.r2.175
#> 1          AG         AA         AA         NN         AA         AG         AG
#> 2          NN         GG         GT         GG         GG         GG         GG
#> 3          CC         CC         CC         CT         CC         CC         CC
#> 4          AA         AG         AA         AA         GG         AA         AG
#> 5          NN         GG         GG         GG         GG         NN         AA
#> 6          AA         AA         NN         AA         AA         AA         AA
#> 7          TT         GT         TT         TT         TT         GT         TT
#> 8          AG         NN         GG         GG         AA         NN         GG
#> 9          CC         AC         NN         AC         NN         CC         AC
#> 10         AA         AA         AA         AA         AA         AA         AA
#> 11         NN         TT         CT         CT         CT         NN         CC
#> 12         AA         AG         AG         AG         AA         NN         NN
#> 13         TT         CC         NN         CT         TT         TT         NN
#> 14         GT         GT         NN         GG         GT         GG         GT
#> 15         AA         AA         AA         AA         AG         AG         NN
#> 16         GT         GT         NN         GG         NN         GT         GT
#> 17         AA         AC         AA         AA         AA         AA         AC
#> 18         AC         AA         AC         AC         CC         CC         CC
#> 19         NN         TT         TT         TT         TT         GT         TT
#> 20         GG         GG         GG         GG         GG         GG         GG
#> 21         AC         CC         AC         AC         AC         NN         NN
#> 22         CC         NN         CC         CT         CC         CT         CC
#> 23         AG         GG         GG         GG         GG         GG         GG
#> 24         NN         AA         AA         AA         AA         NN         NN
#> 25         NN         AG         AG         AA         AG         GG         GG
#> 26         CT         TT         CT         TT         TT         CT         TT
#> 27         AA         AA         AA         NN         AA         AA         AA
#> 28         CT         NN         CT         CT         CT         TT         CT
#> 29         NN         AA         AC         CC         AA         NN         AC
#> 30         CG         GG         GG         GG         GG         NN         CG
#> 31         CG         CG         NN         CC         CG         CG         NN
#> 32         AA         NN         AC         AA         AA         AA         AA
#> 33         CC         CC         TT         CT         TT         CC         TT
#> 34         CC         CG         CG         CC         CG         GG         CC
#> 35         AA         AA         AA         AA         AA         NN         AA
#> 36         AA         AG         AA         NN         NN         AG         AA
#> 37         CG         GG         CC         CG         CG         NN         GG
#> 38         GG         GT         GG         NN         GT         GT         GT
#> 39         TT         TT         AT         TT         TT         TT         NN
#> 40         GG         NN         GG         CG         GG         NN         NN
#> 41         NN         GG         GG         AG         GG         NN         GG
#> 42         AG         AA         AG         AG         AG         NN         AG
#> 43         AG         AA         AA         AA         AA         AA         AG
#> 44         TT         TT         TT         AT         AT         TT         AT
#> 45         NN         GG         NN         GG         GG         GG         GG
#> 46         CC         CC         CC         NN         CC         CC         NN
#> 47         NN         CG         CG         GG         CC         CC         CC
#> 48         NN         CC         CT         CT         TT         CT         CT
#> 49         GG         GG         GG         GG         GG         GG         GG
#> 50         CC         NN         CC         CC         CC         CC         CT
#> 51         TT         CT         TT         NN         CC         TT         TT
#> 52         AA         AA         AA         AA         AG         AA         AA
#> 53         AA         AA         AA         AC         AC         AA         AA
#> 54         GG         GT         GG         GG         GG         GG         GG
#> 55         NN         GG         GG         NN         AG         GG         NN
#> 56         AA         AG         AG         AG         AA         AA         AA
#> 57         GG         GG         GG         GG         GG         NN         GG
#> 58         CC         CC         CC         CT         CT         CC         CC
#> 59         GG         AG         AG         GG         AG         AG         NN
#> 60         AC         AA         AC         AC         AA         AC         AC
#> 61         TT         TT         TT         TT         NN         TT         TT
#> 62         GG         GG         GG         GG         GG         CG         GG
#> 63         NN         AC         NN         CC         AA         AC         AA
#> 64         NN         AT         NN         TT         AT         AT         AT
#> 65         CC         CT         CT         CC         TT         CT         CC
#> 66         AC         CC         AA         AA         AC         NN         NN
#> 67         GG         GG         GG         GG         GG         GG         GG
#> 68         AA         AT         AT         TT         AA         AT         NN
#> 69         AA         AG         AA         AA         NN         AG         AA
#> 70         AA         AC         AC         AA         AC         AC         NN
#> 71         AA         AT         AA         AT         AT         AA         AA
#> 72         AA         AG         AA         AA         AG         AA         AG
#> 73         TT         TT         TT         NN         GT         TT         TT
#> 74         NN         GG         GT         GT         GG         NN         GT
#> 75         NN         CC         CC         CC         NN         CC         CT
#> 76         NN         AT         AT         AT         AA         AA         NN
#> 77         GG         GT         GG         GG         NN         GT         GG
#> 78         TT         TT         TT         TT         TT         TT         CC
#> 79         CT         CT         CT         NN         TT         TT         TT
#> 80         CT         TT         TT         TT         TT         TT         NN
#> 81         TT         TT         TT         TT         TT         TT         TT
#> 82         CC         CT         CT         CC         CC         CC         NN
#> 83         GG         AG         AG         GG         GG         AG         AG
#> 84         NN         AG         AG         AG         AA         NN         AA
#> 85         TT         CT         TT         CT         CT         CT         CT
#> 86         NN         GG         AA         AA         NN         AG         AG
#> 87         CG         CG         GG         CC         CC         CG         CG
#> 88         CC         CC         CC         CC         CT         CC         CC
#> 89         NN         CT         CT         TT         TT         TT         CT
#> 90         AG         NN         AG         GG         GG         NN         NN
#> 91         CT         TT         TT         TT         CT         NN         NN
#> 92         CC         CC         CC         CC         CC         CC         CG
#> 93         GT         GG         GG         GG         GG         NN         GG
#> 94         CC         AC         CC         NN         CC         CC         CC
#> 95         GG         GG         AG         NN         GG         GG         GG
#> 96         GT         GT         GG         GG         GG         GT         GG
#> 97         AG         GG         AG         AA         AG         GG         NN
#> 98         AT         AT         AA         AT         AT         AA         AT
#>    SMR.r2.256 SMR.r2.263 SMR.r2.296 UPD.r2.120 UPD.r2.121 UPD.r2.122 UPD.r2.123
#> 1          AG         AA         GG         AA         AA         AG         AG
#> 2          GG         GG         GG         GG         GG         GG         GG
#> 3          CT         CC         CC         CT         CT         CC         CC
#> 4          AG         GG         AA         AG         GG         AG         AG
#> 5          AG         GG         AG         GG         NN         AA         AG
#> 6          AA         AA         AA         AC         AC         AC         AC
#> 7          TT         GT         GT         TT         TT         TT         TT
#> 8          GG         AG         GG         NN         GG         AG         AG
#> 9          CC         CC         CC         CC         CC         AC         CC
#> 10         AA         AA         AA         AA         AA         AA         AA
#> 11         CC         CT         TT         CT         CT         TT         TT
#> 12         AG         AA         AA         AG         AA         AG         AG
#> 13         NN         TT         TT         TT         CT         TT         TT
#> 14         GG         GT         GG         GG         GT         GT         GT
#> 15         AA         NN         AG         AG         AG         GG         NN
#> 16         GG         GT         GT         GG         GG         GT         GT
#> 17         NN         AC         AA         AC         AA         AA         AC
#> 18         CC         CC         CC         AA         CC         AA         CC
#> 19         TT         TT         TT         GT         TT         TT         TT
#> 20         GG         GG         GG         GG         GG         NN         AG
#> 21         CC         CC         AC         CC         AC         AC         NN
#> 22         CT         CC         CC         CC         TT         CT         CT
#> 23         NN         GG         GG         GG         GG         GG         GG
#> 24         AA         AA         AA         AA         AA         AA         AA
#> 25         AG         AG         AG         GG         AG         AG         AA
#> 26         TT         NN         CT         CT         CT         CT         TT
#> 27         AC         AC         AA         AA         AA         AC         AA
#> 28         CT         CC         CC         TT         CT         CT         TT
#> 29         AA         AA         AA         AA         AA         AA         AC
#> 30         GG         GG         GG         GG         GG         GG         NN
#> 31         CC         CC         CG         GG         CG         CG         CC
#> 32         NN         NN         AC         AC         AA         AA         AC
#> 33         CC         CT         CT         CC         TT         CC         TT
#> 34         CC         CG         GG         CG         CC         CG         CC
#> 35         AA         AA         AA         AC         AA         AA         AA
#> 36         AA         AG         AG         AA         AA         AA         AG
#> 37         GG         CG         CC         GG         GG         GG         GG
#> 38         GT         TT         NN         GG         GG         GT         GG
#> 39         AT         TT         TT         TT         AT         TT         TT
#> 40         CG         CG         CG         CG         CG         GG         CG
#> 41         NN         GG         AG         AA         NN         AA         NN
#> 42         AA         AG         AG         AG         AG         GG         GG
#> 43         AA         AA         AG         NN         AA         AG         AG
#> 44         TT         TT         TT         TT         AT         TT         NN
#> 45         GG         GG         NN         GG         GG         GG         GG
#> 46         CC         CC         AA         CC         AC         CC         NN
#> 47         CC         CG         GG         CC         CC         CG         CC
#> 48         NN         CC         NN         CT         NN         TT         CT
#> 49         GG         GG         GG         GG         GG         GG         GG
#> 50         CC         CC         CC         CC         CC         CT         CC
#> 51         NN         NN         TT         CC         TT         CC         CT
#> 52         AA         AA         AA         GG         GG         AA         AG
#> 53         AC         AC         AA         AA         CC         AA         AC
#> 54         GG         GG         GG         GG         GG         GG         GG
#> 55         GG         GG         AG         AA         AG         AG         AG
#> 56         AG         NN         AA         AG         AA         AG         AG
#> 57         GG         GG         NN         GG         GG         GG         GG
#> 58         CT         CC         NN         NN         TT         CT         NN
#> 59         AA         AG         AA         GG         AG         GG         AG
#> 60         NN         AA         AA         AC         AA         NN         AC
#> 61         CT         TT         TT         TT         CT         TT         TT
#> 62         GG         GG         GG         GG         NN         GG         GG
#> 63         AA         AC         NN         AC         AC         AC         AA
#> 64         AT         TT         TT         AA         TT         TT         AT
#> 65         CC         TT         CC         TT         CC         TT         CC
#> 66         AA         CC         AC         AC         AC         AA         AA
#> 67         GG         NN         GG         GG         GG         GG         GG
#> 68         NN         NN         AA         AA         TT         AT         AA
#> 69         AA         AG         NN         AA         AA         AA         AA
#> 70         AA         AC         AC         AA         AC         NN         AC
#> 71         AA         AT         AA         AA         AA         AA         AA
#> 72         AA         AA         AG         AG         GG         AG         AG
#> 73         TT         TT         NN         GT         TT         TT         GT
#> 74         GT         GG         GG         GT         GG         GG         GG
#> 75         NN         CT         NN         CT         CC         CC         CT
#> 76         AT         AA         AA         AT         AA         AT         AT
#> 77         GG         GT         GT         GG         GG         GG         GG
#> 78         TT         TT         TT         CT         TT         TT         TT
#> 79         TT         CT         TT         TT         CT         NN         CT
#> 80         TT         TT         TT         TT         TT         TT         CT
#> 81         TT         TT         AT         TT         TT         TT         TT
#> 82         CC         CT         CT         CT         TT         CC         CT
#> 83         GG         AA         NN         GG         GG         AG         GG
#> 84         AG         AA         NN         AG         AA         AG         AG
#> 85         TT         TT         CT         CT         TT         NN         TT
#> 86         AG         AG         NN         NN         AG         AG         AG
#> 87         NN         GG         CG         GG         CC         CC         CG
#> 88         CC         CC         NN         CC         CC         CC         CT
#> 89         NN         CT         CT         CT         CT         CC         TT
#> 90         GG         NN         GG         NN         GG         GG         GG
#> 91         CT         TT         CT         TT         TT         TT         CT
#> 92         CC         CC         CC         CC         CC         CC         CC
#> 93         GG         GG         GG         GG         GT         GG         GG
#> 94         CC         AC         CC         CC         CC         CC         NN
#> 95         GG         GG         AG         GG         GG         GG         GG
#> 96         GG         GT         GT         GT         GG         GG         GT
#> 97         NN         GG         AA         AG         GG         AA         NN
#> 98         AT         AA         AA         AA         AT         AA         AT
#>    UPD.r2.124 UPD.r2.129 UPD.r2.130 UPD.r2.131 UPD.r2.133 UPD.r2.135 UPD.r2.160
#> 1          AG         AA         AG         AG         GG         AG         NN
#> 2          GG         GG         GG         GG         GG         GG         GG
#> 3          CC         CC         CC         CC         CC         CC         CC
#> 4          AG         AG         AG         NN         NN         AG         AG
#> 5          AG         NN         NN         NN         AA         AG         AG
#> 6          CC         AC         NN         AA         AC         AC         AA
#> 7          TT         TT         TT         TT         TT         TT         TT
#> 8          NN         GG         NN         NN         GG         AG         AG
#> 9          NN         CC         AC         CC         AC         AC         CC
#> 10         AG         NN         AG         AA         AA         AA         AA
#> 11         CT         NN         CT         CT         CT         CT         CT
#> 12         AA         NN         AA         AG         NN         NN         NN
#> 13         TT         TT         TT         TT         CT         CT         CT
#> 14         GT         GT         GG         GG         GT         GG         TT
#> 15         AG         AG         AG         NN         AG         AG         AG
#> 16         GG         GT         NN         GG         NN         GG         GG
#> 17         AA         NN         AA         AA         AA         AC         NN
#> 18         CC         CC         AA         CC         AC         AC         CC
#> 19         TT         TT         TT         TT         GT         TT         TT
#> 20         GG         GG         GG         GG         GG         GG         GG
#> 21         CC         NN         AC         CC         AC         CC         AC
#> 22         CT         CT         CC         NN         CC         CT         CT
#> 23         AG         GG         GG         GG         GG         GG         GG
#> 24         AA         AA         AA         AG         AG         AA         NN
#> 25         GG         AG         AG         GG         AA         GG         AG
#> 26         CT         TT         CT         CT         TT         CT         NN
#> 27         AA         AC         AA         AA         AA         AA         AA
#> 28         TT         TT         CT         CT         NN         TT         TT
#> 29         AA         NN         AA         AA         NN         AA         AA
#> 30         GG         GG         CG         CG         NN         GG         NN
#> 31         CC         CC         CG         CG         GG         CC         GG
#> 32         AC         AA         AA         AC         AC         AC         AA
#> 33         CC         CT         CC         NN         CT         CC         CC
#> 34         CC         NN         CC         CC         CC         CC         NN
#> 35         AA         AA         AA         AA         AA         AA         AA
#> 36         AA         AA         AA         AG         AA         NN         AA
#> 37         CG         CG         GG         CG         CG         CG         CG
#> 38         GT         GG         GG         TT         TT         GG         GT
#> 39         TT         NN         NN         AT         NN         TT         AT
#> 40         GG         NN         CG         CG         NN         CC         GG
#> 41         GG         NN         AG         AA         AA         AG         NN
#> 42         AG         GG         AG         AA         GG         NN         AG
#> 43         NN         AG         AA         AG         AA         AG         AG
#> 44         TT         AT         AT         NN         TT         NN         TT
#> 45         GG         GG         GG         GG         GG         GG         GG
#> 46         CC         CC         AC         AA         AC         CC         AC
#> 47         CC         CC         CC         CG         CC         CC         CC
#> 48         CT         CT         TT         CT         CT         CT         TT
#> 49         GG         GG         NN         CG         GG         GG         NN
#> 50         CC         CC         CT         CC         CC         CC         NN
#> 51         TT         CT         CT         CT         TT         CT         NN
#> 52         AG         NN         AG         AG         AA         AA         AG
#> 53         NN         NN         AA         AC         AC         AC         AC
#> 54         GG         GG         GG         GG         GG         GG         GG
#> 55         GG         GG         AG         AA         GG         NN         NN
#> 56         AA         AG         AA         AG         AA         AA         AA
#> 57         GG         NN         GG         NN         GG         NN         GG
#> 58         NN         CT         CT         CC         CC         CC         NN
#> 59         GG         GG         GG         AG         AG         AG         AA
#> 60         AC         AC         AC         AA         AA         AA         AC
#> 61         TT         NN         TT         TT         TT         NN         TT
#> 62         GG         GG         CG         CG         CG         GG         GG
#> 63         CC         AA         AC         AC         AA         NN         AA
#> 64         AT         TT         TT         AT         AT         AT         AT
#> 65         CT         CC         CT         NN         CC         CC         CC
#> 66         AC         NN         AA         AA         AA         CC         CC
#> 67         GG         GG         GG         GG         NN         NN         AG
#> 68         AT         AT         AA         AA         AT         AT         NN
#> 69         AA         AA         AA         AA         AA         AA         AA
#> 70         AA         AA         AC         AC         AA         AA         AC
#> 71         AA         AA         AA         AA         AT         AA         AA
#> 72         NN         AG         GG         NN         GG         AG         GG
#> 73         GT         TT         GG         TT         TT         GT         TT
#> 74         GG         GG         NN         GG         GG         GG         GT
#> 75         CT         CT         CT         CC         CC         TT         CT
#> 76         AA         NN         AA         AT         AA         AA         AA
#> 77         GG         GG         GG         GG         GG         GG         GG
#> 78         NN         TT         TT         TT         NN         NN         TT
#> 79         CT         CT         CT         TT         CT         TT         NN
#> 80         TT         TT         NN         CT         TT         CT         CT
#> 81         TT         TT         TT         AT         TT         TT         TT
#> 82         CT         CT         CT         CT         TT         CT         NN
#> 83         AG         AA         AG         AG         AA         AG         AG
#> 84         NN         AG         AA         AG         AA         NN         AA
#> 85         TT         CT         CT         CC         TT         CT         CT
#> 86         AG         AA         AA         AG         AA         NN         AG
#> 87         NN         GG         CG         NN         CG         NN         CG
#> 88         NN         CC         CC         CT         CC         CC         CC
#> 89         TT         CT         CC         CC         CC         CC         CT
#> 90         AG         NN         NN         NN         AG         GG         AG
#> 91         CT         TT         TT         CT         TT         TT         CT
#> 92         CG         CC         CG         CC         CC         NN         CG
#> 93         GG         GT         GT         GT         GG         GG         NN
#> 94         CC         CC         CC         CC         NN         AC         CC
#> 95         GG         GG         GG         AG         GG         GG         GG
#> 96         GT         NN         GG         GG         GT         NN         GT
#> 97         NN         GG         AG         AG         GG         NN         NN
#> 98         AT         AT         TT         TT         TT         NN         TT
#>    UPD.r2.162 UPD.r2.174 UPD.r2.189 UPD.r2.201
#> 1          GG         NN         GG         AG
#> 2          GG         GG         GG         GG
#> 3          CC         CC         CC         CC
#> 4          AG         AA         AA         AA
#> 5          NN         GG         NN         GG
#> 6          AC         AA         AA         AA
#> 7          NN         TT         TT         TT
#> 8          NN         AG         GG         GG
#> 9          NN         CC         AC         AA
#> 10         AA         AG         AA         AG
#> 11         CT         CT         CT         TT
#> 12         AG         AA         AG         NN
#> 13         TT         TT         NN         NN
#> 14         GT         GG         GT         GT
#> 15         AG         AG         AA         AA
#> 16         GG         TT         NN         GG
#> 17         AA         AC         AA         AA
#> 18         AC         AA         AC         AC
#> 19         TT         TT         TT         TT
#> 20         GG         GG         AG         GG
#> 21         NN         AC         AA         CC
#> 22         CT         CT         CC         CC
#> 23         GG         GG         GG         GG
#> 24         AG         AA         AA         AA
#> 25         AG         AA         AG         GG
#> 26         NN         NN         TT         TT
#> 27         AA         NN         AA         AA
#> 28         CT         CT         CT         CT
#> 29         AA         AA         AC         AC
#> 30         GG         CG         GG         CG
#> 31         CG         CG         CC         CG
#> 32         NN         AA         NN         AC
#> 33         CC         CT         NN         TT
#> 34         CG         NN         CC         NN
#> 35         NN         AA         AA         AC
#> 36         NN         AA         AA         AA
#> 37         NN         CG         NN         GG
#> 38         GG         GG         GG         NN
#> 39         TT         NN         TT         NN
#> 40         GG         NN         CC         NN
#> 41         NN         GG         AG         AG
#> 42         AG         GG         AG         AG
#> 43         AA         AA         AA         AA
#> 44         TT         AT         AT         TT
#> 45         GG         GT         GG         GT
#> 46         CC         AA         NN         AC
#> 47         GG         CC         CC         CG
#> 48         CT         NN         CC         TT
#> 49         GG         CG         GG         GG
#> 50         CC         CC         CC         CC
#> 51         CT         CT         TT         TT
#> 52         AG         AA         AA         AA
#> 53         NN         AA         AC         CC
#> 54         GG         GG         GG         GG
#> 55         AA         AG         AG         GG
#> 56         AG         AA         AG         AA
#> 57         GG         NN         NN         GG
#> 58         CC         CC         CT         CC
#> 59         AG         AA         NN         NN
#> 60         AC         AC         AA         AC
#> 61         TT         CT         TT         NN
#> 62         GG         GG         GG         GG
#> 63         AA         AA         NN         AA
#> 64         TT         AT         AT         AT
#> 65         CT         CT         CT         TT
#> 66         AC         CC         NN         AC
#> 67         GG         GG         GG         GG
#> 68         TT         AA         AA         AT
#> 69         AA         AA         AA         AA
#> 70         NN         NN         AC         AC
#> 71         AA         AT         AA         AA
#> 72         GG         GG         AG         AG
#> 73         TT         TT         NN         GT
#> 74         GG         GT         NN         GG
#> 75         TT         CC         NN         TT
#> 76         AT         AA         AT         AA
#> 77         GG         GG         GG         GG
#> 78         CT         CT         TT         TT
#> 79         CT         TT         TT         TT
#> 80         TT         TT         NN         CT
#> 81         TT         TT         NN         TT
#> 82         CC         CT         NN         CC
#> 83         AG         GG         AG         AG
#> 84         AG         AA         AG         AA
#> 85         CT         CT         CC         CT
#> 86         NN         AG         AA         AG
#> 87         GG         GG         GG         CC
#> 88         CC         CC         NN         CC
#> 89         NN         CC         CT         CT
#> 90         GG         NN         AG         GG
#> 91         NN         CT         CT         CT
#> 92         CC         CC         CC         CC
#> 93         GG         NN         NN         GG
#> 94         CC         CC         AC         CC
#> 95         GG         GG         GG         GG
#> 96         GG         GT         GT         TT
#> 97         AA         AA         AG         AG
#> 98         TT         AT         AT         AT
```
