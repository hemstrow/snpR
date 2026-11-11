# Pairwise LD from SNP data.

`calc_pairwise_ld` calculates LD between each pair of SNPs.

## Usage

``` r
calc_pairwise_ld(
  x,
  facets = NULL,
  subfacets = NULL,
  ss = FALSE,
  par = FALSE,
  CLD = "only",
  use.ME = FALSE,
  sigma = 1e-04,
  window_sigma = NULL,
  window_step = window_sigma * 2,
  window_gaussian = TRUE,
  window_triple_sigma = TRUE,
  verbose = FALSE,
  .prox_only = FALSE
)
```

## Arguments

- x:

  snpRdata. Input SNP data. Note that a SNP column containing snp
  position in base pairs named 'position' is required.

- facets:

  character. Categorical metadata variables by which to break up
  analysis. See
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md)
  for more details.

- subfacets:

  character, default NULL. Subsets the facet levels to run. Given as a
  named list: list(fam = A) will run only fam A, list(fam = c("A", "B"),
  chr = 1) will run only fams A and B on chromosome 1. list(fam = "A",
  pop = "ASP") will run samples in either fam A or pop ASP, list(fam.pop
  = "A.ASP") will run only samples in fam A and pop ASP.

- ss:

  numeric, default NULL. Number of snps to subsample.

- par:

  numeric or FALSE, default FALSE. If numeric, the number of cores to
  use for parallel processing.

- CLD:

  TRUE, FALSE, or "only", default "only". Specifies if the CLD method
  should be used either in addition to or instead of default methods.
  See details.

- use.ME:

  logical, default FALSE. Specifies if the Minimization-Expectation
  haplotype estimation should be used. See details.

- sigma:

  numeric, default 0.0001. If the ME method is used, specifies the
  minimum difference required between steps before haplotype frequencies
  are accepted.

- window_sigma:

  numeric, default NULL. Size of windows in kb within which to calculate
  ld values, if requested. Windows will be two times `window_sigma` in
  size unless `window_triple_sigma` is true, in which case they will be
  six times `window_sigma`.

- window_step:

  numeric or NULL, default two times `window_sigma` (non-overlapping
  windows if `window_triple_sigma` is `FALSE`). Size of the steps
  between windows, in kb.

- window_gaussian:

  logical, default TRUE. If TRUE, windows will be gaussian-smoothed.
  Otherwise, raw averages will be returned. See
  [`calc_smoothed_averages`](https://hemstrow.github.io/snpR/reference/calc_smoothed_averages.md)
  for details.

- window_triple_sigma:

  logical, default TRUE. If TRUE, `window_sigma` values will be tripled
  prior to averaging.

- verbose:

  Logical, default FALSE. If TRUE, some progress updates will be
  reported.

- .prox_only:

  Logical, default FALSE. Primarily for internal use. if TRUE returns
  *ONLY* a proximity table of LD values, not a
  [`snpRdata`](https://hemstrow.github.io/snpR/reference/snpRdata-class.md)
  object.

## Value

a snpRdata object with linkage results stored in the pairwise.LD slot.
Specifically, this slot will contain a list containing any LD matrices
in a nested list broken down facet then by facet levels and a data.frame
containing all pairwise comparisons, their metadata, and calculated
statistics in long format for easy plotting.

## Details

Calculates pairwise linkage disequilibrium between pairs of SNPs using
several different methods. By default uses the Burrow's Composite
Linkage Disequilibrium method.

If cld is not "only", haplotypes are estimated either via direct count
after removing all "0101" double heterozygote haplotypes (if use.ME is
FALSE) or via the Minimization-Expectation method described in
Excoffier, L., and Slatkin, M. (1995). Note that while the latter method
is likely more accurate, it can be *very* slow and often produces
qualitatively equivalent results, and so is not preferred during casual
or preliminary analysis. Either method will calculate D', r-squared, and
the p-value for that r-squared.

Since this process involves many pairwise comparisons, it can be very
slow. As an alternative, average LD values can be calculated within
sliding windows using the `window_` family of arguments. This will be
substantially faster, but individual snp/snp LD values will not be
returned. See
[`calc_smoothed_averages`](https://hemstrow.github.io/snpR/reference/calc_smoothed_averages.md)
for details.

In contrast, Burrow's Composite Linkage Disequilibrium (CLD) can be
calculated very quickly via the
[`cor`](https://rdrr.io/r/stats/cor.html) function from base R.
`calc_pairwise_ld` will perform this method alongside the other methods
if cld = TRUE and by itself if cld = "only". For most analyses, this
will be sufficient and much faster than the other methods. This is the
default behavior.

The data can be broken up categorically by either SNP and/or sample
metadata, as described in
[`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).

Heatmaps of the resulting data can be easily plotted using
[`plot_pairwise_ld_heatmap`](https://hemstrow.github.io/snpR/reference/plot_pairwise_ld_heatmap.md).

## References

Dimitri Zaykin (2004). *Genetic Epidemiology*

Excoffier, L., and Slatkin, M. (1995). *Molecular Biology and Evolution*

Lewontin (1964). *Genetics*

## Author

William Hemstrom

Keming Su

## Examples

``` r
if (FALSE) { # \dontrun{
# not run, slow
## CLD
x <- calc_pairwise_ld(stickSNPs, facets = "chr.pop")
get.snpR.stats(x, "chr.pop", "LD")

## standard haplotype frequency estimation
x <- calc_pairwise_ld(stickSNPs, facets = "chr.pop", CLD = FALSE)
get.snpR.stats(x, "chr.pop", "LD")
} # }

# subset for specific subfacets (ASP and OPL, chromosome IX)
x <- calc_pairwise_ld(stickSNPs, facets = "chr.pop",
                      subfacets = list(pop = c("ASP", "OPL"), 
                                       chr = "groupIX"))
get.snpR.stats(x, "chr.pop", "LD")
#> $prox
#>       s1_chr s1_position s1_.snp.id  s2_chr s2_position s2_.snp.id proximity
#>       <char>       <num>      <int>  <char>       <num>      <int>     <num>
#>   1: groupIX     2075079          7 groupIX     5287711         22   3212632
#>   2: groupIX     2075079          7 groupIX     5960933         24   3885854
#>   3: groupIX     5287711         22 groupIX     5960933         24    673222
#>   4: groupIX     2075079          7 groupIX     6515796         28   4440717
#>   5: groupIX     5287711         22 groupIX     6515796         28   1228085
#>  ---                                                                        
#> 117: groupIX    14164215         65 groupIX    18339767         88   4175552
#> 118: groupIX    14868533         69 groupIX    18339767         88   3471234
#> 119: groupIX    15190693         70 groupIX    18339767         88   3149074
#> 120: groupIX    16306966         78 groupIX    18339767         88   2032801
#> 121: groupIX    17929225         87 groupIX    18339767         88    410542
#>              CLD sample.facet sample.subfacet     S
#>            <num>       <char>          <char> <num>
#>   1: 0.002840909          pop             ASP     9
#>   2: 0.012345679          pop             ASP    10
#>   3: 0.000000000          pop             ASP     8
#>   4: 0.428571429          pop             ASP     8
#>   5: 0.008333333          pop             ASP     7
#>  ---                                               
#> 117: 0.082206796          pop             OPL    16
#> 118: 0.374400056          pop             OPL    19
#> 119: 0.084912281          pop             OPL    16
#> 120: 0.337214543          pop             OPL    17
#> 121: 0.166982204          pop             OPL    20
#> 
#> $matrices
#> $matrices$ASP
#> $matrices$ASP$groupIX
#> $matrices$ASP$groupIX$CLD
#>          2075079     5287711    5960933     6515796     8373442    9284427
#> 2075079       NA 0.002840909 0.01234568 0.428571429 0.047619048 0.20250000
#> 5287711       NA          NA 0.00000000 0.008333333 0.008547009 0.00000000
#> 5960933       NA          NA         NA 0.047619048          NA 0.01234568
#> 6515796       NA          NA         NA          NA 0.008333333 0.71014493
#> 8373442       NA          NA         NA          NA          NA 0.02116402
#> 9284427       NA          NA         NA          NA          NA         NA
#> 11904651      NA          NA         NA          NA          NA         NA
#> 14164215      NA          NA         NA          NA          NA         NA
#> 14868533      NA          NA         NA          NA          NA         NA
#> 15190693      NA          NA         NA          NA          NA         NA
#> 16306966      NA          NA         NA          NA          NA         NA
#> 17929225      NA          NA         NA          NA          NA         NA
#> 18339767      NA          NA         NA          NA          NA         NA
#>          11904651    14164215 14868533    15190693 16306966     17929225
#> 2075079        NA 0.057142857       NA 0.047619048       NA 0.1833333333
#> 5287711        NA 0.011363636       NA 0.011363636       NA 0.0005165289
#> 5960933        NA 0.074074074       NA 0.062500000       NA 0.1851851852
#> 6515796        NA 0.022222222       NA 0.066666667       NA 0.0769230769
#> 8373442        NA 0.007936508       NA 0.035714286       NA 0.0476190476
#> 9284427        NA 0.008035714       NA 0.277504105       NA 0.0458333333
#> 11904651       NA          NA       NA          NA       NA           NA
#> 14164215       NA          NA       NA 0.002267574       NA 0.0654761905
#> 14868533       NA          NA       NA          NA       NA           NA
#> 15190693       NA          NA       NA          NA       NA 0.3174603175
#> 16306966       NA          NA       NA          NA       NA           NA
#> 17929225       NA          NA       NA          NA       NA           NA
#> 18339767       NA          NA       NA          NA       NA           NA
#>           18339767
#> 2075079  0.5714286
#> 5287711  0.1250000
#> 5960933         NA
#> 6515796  0.3500000
#> 8373442  0.0000000
#> 9284427  0.1739130
#> 11904651        NA
#> 14164215 0.2500000
#> 14868533        NA
#> 15190693 0.0000000
#> 16306966        NA
#> 17929225 0.5714286
#> 18339767        NA
#> 
#> $matrices$ASP$groupIX$S
#>     1  2  3  4  5  6  7  8  9 10 11 12 13
#> 1  NA  9 10  8 10 11 NA 11 NA 10 NA 11  8
#> 2  NA NA  8  7  8  9 NA  9 NA  9 NA  9  7
#> 3  NA NA NA  8 NA 10 NA 10 NA  9 NA 10 NA
#> 4  NA NA NA NA  7  8 NA  8 NA  7 NA  8  7
#> 5  NA NA NA NA NA 10 NA 10 NA  9 NA 10  8
#> 6  NA NA NA NA NA NA NA 11 NA 10 NA 11  8
#> 7  NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 8  NA NA NA NA NA NA NA NA NA 10 NA 11  8
#> 9  NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 10 NA NA NA NA NA NA NA NA NA NA NA 10  7
#> 11 NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 12 NA NA NA NA NA NA NA NA NA NA NA NA  8
#> 13 NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 
#> 
#> 
#> $matrices$OPL
#> $matrices$OPL$groupIX
#> $matrices$OPL$groupIX$CLD
#>          2075079    5287711    5960933      6515796      8373442    9284427
#> 2075079       NA 0.04081863 0.09316770 0.0028409091 1.358696e-01 0.22321429
#> 5287711       NA         NA 0.08427527 0.2610797744 1.084071e-01 0.01282051
#> 5960933       NA         NA         NA 0.0009832842 8.656806e-02 0.06576087
#> 6515796       NA         NA         NA           NA 6.459531e-05 0.01481481
#> 8373442       NA         NA         NA           NA           NA 0.22123894
#> 9284427       NA         NA         NA           NA           NA         NA
#> 11904651      NA         NA         NA           NA           NA         NA
#> 14164215      NA         NA         NA           NA           NA         NA
#> 14868533      NA         NA         NA           NA           NA         NA
#> 15190693      NA         NA         NA           NA           NA         NA
#> 16306966      NA         NA         NA           NA           NA         NA
#> 17929225      NA         NA         NA           NA           NA         NA
#> 18339767      NA         NA         NA           NA           NA         NA
#>             11904651     14164215     14868533    15190693     16306966
#> 2075079  0.122807018 0.0785256410 0.0865384615 0.318181818 0.2900361384
#> 5287711  0.015081564 0.0660676533 0.1101937851 0.003906250 0.0625000000
#> 5960933  0.006944444 0.1085420681 0.1008064516 0.001122334 0.1358317215
#> 6515796  0.012987013 0.0004681648 0.0768344218 0.064935065 0.0260145682
#> 8373442  0.361368496 0.0053960717 0.0001194458 0.050000000 0.0006684492
#> 9284427  0.058823529 0.0142269253 0.0149390244 0.000000000 0.0131293403
#> 11904651          NA 0.0075483092 0.0220120065 0.045454545 0.0188153310
#> 14164215          NA           NA 0.2875000000 0.210069444 0.3928571429
#> 14868533          NA           NA           NA 0.818181818 1.0000000000
#> 15190693          NA           NA           NA          NA 0.8100000000
#> 16306966          NA           NA           NA          NA           NA
#> 17929225          NA           NA           NA          NA           NA
#> 18339767          NA           NA           NA          NA           NA
#>             17929225     18339767
#> 2075079  0.013035382 0.2355769231
#> 5287711  0.008264766 0.0001002406
#> 5960933  0.037762238 0.0129870130
#> 6515796  0.010546875 0.0496222527
#> 8373442  0.009708477 0.0305157156
#> 9284427  0.035840108 0.2129629630
#> 11904651 0.022012007 0.0282186949
#> 14164215 0.134101624 0.0822067958
#> 14868533 0.073107004 0.3744000565
#> 15190693 0.009977827 0.0849122807
#> 16306966 0.029239766 0.3372145433
#> 17929225          NA 0.1669822037
#> 18339767          NA           NA
#> 
#> $matrices$OPL$groupIX$S
#>     1  2  3  4  5  6  7  8  9 10 11 12 13
#> 1  NA 17 20 13 18 15 20 15 19 15 18 20 15
#> 2  NA NA 21 14 18 14 20 15 20 17 18 19 15
#> 3  NA NA NA 18 23 18 25 18 24 20 21 24 19
#> 4  NA NA NA NA 18 14 19 14 18 15 16 17 15
#> 5  NA NA NA NA NA 18 24 18 22 19 20 23 18
#> 6  NA NA NA NA NA NA 18 16 18 15 17 19 17
#> 7  NA NA NA NA NA NA NA 19 24 21 22 24 19
#> 8  NA NA NA NA NA NA NA NA 19 15 17 19 16
#> 9  NA NA NA NA NA NA NA NA NA 20 21 23 19
#> 10 NA NA NA NA NA NA NA NA NA NA 19 19 16
#> 11 NA NA NA NA NA NA NA NA NA NA NA 20 17
#> 12 NA NA NA NA NA NA NA NA NA NA NA NA 20
#> 13 NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 
#> 
#> 
#> 

if (FALSE) { # \dontrun{
## not run, really slow
# ME haplotype estimation
x <- calc_pairwise_ld(stickSNPs, facets = "chr.pop", 
                      CLD = FALSE, use.ME = TRUE,
                      subfacets = list(pop = c("ASP", "OPL"), 
                                       chr = "groupIX"))
get.snpR.stats(x, "chr.pop", "LD")
} # }
```
