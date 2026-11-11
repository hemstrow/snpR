# FST from SNP data.

`calc_pairwise_fst` calculates pairwise FST for each SNP for each
possible pairwise combination of populations. `calc_global_fst`
calculates FST for each facet globally across all subfacet levels.

## Usage

``` r
calc_pairwise_fst(
  x,
  facets,
  method = "wc",
  boot = FALSE,
  boot_par = FALSE,
  zfst = FALSE,
  fst_over_one_minus_fst = FALSE,
  keep_components = FALSE,
  cleanup = TRUE,
  verbose = FALSE
)

calc_global_fst(
  x,
  facets,
  boot = FALSE,
  boot_par = FALSE,
  zfst = FALSE,
  fst_over_one_minus_fst = FALSE,
  keep_components = FALSE,
  verbose = FALSE
)
```

## Arguments

- x:

  snpRdata. Input SNP data.

- facets:

  character. Categorical metadata variables by which to break up
  analysis. See
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md)
  for more details.

- method:

  character, default "wc". Defines the FST estimator to use. Options:

  - wc: Weir and Cockerham (1984).

  - genepop: Rousset (2008), uses the genepop package.

- boot:

  numeric or FALSE, default FALSE. The number of bootstraps to do. See
  details.

- boot_par:

  numeric or FALSE, default FALSE. If a number, bootstraps will be
  processed in parallel using the supplied number of cores.

- zfst:

  logical, default FALSE. If TRUE, z-distributed Fst scores (zFST) will
  be calculated, equal to (fst - mean(fst))/sd(fst) within each group.
  The resulting values will be in the column "zfst", accessible using
  the usual
  [`get.snpR.stats`](https://hemstrow.github.io/snpR/reference/get.snpR.stats.md)
  method.

- fst_over_one_minus_fst:

  logical, default FALSE. If TRUE, fst/(1-fst) will be calculated, and
  will be in the column "fst_id" accessible using the usual
  [`get.snpR.stats`](https://hemstrow.github.io/snpR/reference/get.snpR.stats.md)
  method.

- keep_components:

  logical, default FALSE. If TRUE, the variance components "a", "b", and
  "c" will be held and accessible from the `$pairwise` element (named
  "var_comp_a", "var_comp_b", and "var_comp_c", respectively) using the
  usual
  [`get.snpR.stats`](https://hemstrow.github.io/snpR/reference/get.snpR.stats.md)
  method. This may be useful if working with very large datasets that
  need to be run with separate objects for each chromosome, etc. for
  memory purposes. Weighted averages can be generated identically to
  those from snpR by taking the weighted mean (via the
  [`weighted.mean`](https://rdrr.io/r/stats/weighted.mean.html)) of "a"
  divided by the sum of the weighted means of "a", "b", and "c" using
  the number of SNPs called in a comparison (returned in the "nk" column
  from
  [`get.snpR.stats`](https://hemstrow.github.io/snpR/reference/get.snpR.stats.md))
  as weights within each population comparison. Note that this is
  different than taking the weighted mean of a/(a + b + c)!

- cleanup:

  logical, default TRUE. If TRUE, any new files created during FST
  calculation will be automatically removed.

- verbose:

  Logical, default FALSE. If TRUE, some progress updates will be
  reported.

## Value

A snpRdata object with pairwise FST as well as the number of total
observations at each SNP in each comparison merged in to the
pairwise.stats slot.

## Details

Calculates FST according to either Weir and Cockerham 1984 or using the
[`Fst`](https://rdrr.io/pkg/genepop/man/Fst.html) function from the
genepop package (see references). Genepop is not supported for global
FST.

If the genepop option is used, several intermediate files will be
created in the default temporary directory (see
[`tempfile`](https://rdrr.io/r/base/tempfile.html)).

The Weir and Cockerham (1984) and genepop methods tend to produce very
similar results both per SNP and per population. Generally, the former
option may be preferred for computational efficiency.

P-values for group level comparisons can be calculated via bootstrapping
using the boot option. Bootstraps are performed via randomly mixing
individuals amongst different levels of the supplied facet, and thus the
null hypothesis is that all groups are panmictic. P-values are
calculated according to
[`randtest`](https://adeverse.github.io/ade4/reference/randtest.html),
although that function is not directly called.

The data can be broken up categorically by either SNP and/or sample
metadata, as described in
[`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).
Since this is a pairwise statistic, at least a single sample level facet
must be provided.

Method Options:

- "wc": Weir and Cockerham 1984.

- "Genepop": As used in genepop, Rousset 2008.

## Functions

- `calc_pairwise_fst()`: Calculate FST across each pair of pairwise
  subfacet comparisons.

- `calc_global_fst()`: Calculate FST globally across all subfacet
  levels.

## References

Weir and Cockerham (1984). *Evolution*

Weir (1990). Genetic data analysis. Sinauer, Sunderland, MA

Rousset (2008). *Molecular Ecology Resources*

## Author

William Hemstrom

## Examples

``` r
# Using Weir and Cockerham 1984's method
x <- calc_pairwise_fst(stickSNPs, "pop")
get.snpR.stats(x, "pop", "fst")
#> $pairwise
#>      facet comparison        chr position .snp.id           fst  nk
#> 1      pop    ASP~CLF    groupVI   212436       1 -1.313321e-02  28
#> 2      pop    ASP~CLF groupXVIII   350070       2 -2.898551e-02  32
#> 3      pop    ASP~CLF    groupXV   854285       3 -6.589963e-02  34
#> 4      pop    ASP~CLF  groupXIII   907279       4 -6.519082e-02  32
#> 5      pop    ASP~CLF    groupVI  1068620       5  1.736111e-01  20
#> 6      pop    ASP~CLF     groupV  1392753       6 -5.095541e-02  32
#> 7      pop    ASP~CLF    groupIX  2075079       7 -2.103960e-02  34
#> 8      pop    ASP~CLF  groupXIII  2106000       8  2.773344e-02  34
#> 9      pop    ASP~CLF   groupXII  2136587       9 -6.838784e-02  34
#> 10     pop    ASP~CLF  groupXVII  2560839      10 -6.589963e-02  34
#> 11     pop    ASP~CLF   groupIII  2750695      11 -6.104129e-02  32
#> 12     pop    ASP~CLF   groupXIX  3008524      12  3.738318e-02  32
#> 13     pop    ASP~CLF   groupXVI  3624417      13  1.724744e-01  30
#> 14     pop    ASP~CLF   groupVII  3667528      14 -5.479452e-02  32
#> 15     pop    ASP~CLF   groupVII  4189120      15  1.437318e-01  34
#> 16     pop    ASP~CLF    groupXX  4538074      16 -7.814045e-02  30
#> 17     pop    ASP~CLF  groupXIII  4706089      17  6.896552e-02  32
#> 18     pop    ASP~CLF   groupIII  4746808      18 -1.888952e-02  30
#> 19     pop    ASP~CLF  groupXVII  5041995      19  2.885321e-01  34
#> 20     pop    ASP~CLF     groupX  5161637      20 -6.770933e-02  30
#> 21     pop    ASP~CLF   groupXII  5178037      21  5.531915e-02  26
#> 22     pop    ASP~CLF    groupIX  5287711      22 -9.831651e-02  28
#> 23     pop    ASP~CLF  groupVIII  5462135      23 -2.898551e-02  32
#> 24     pop    ASP~CLF    groupIX  5960933      24 -2.898551e-02  32
#> 25     pop    ASP~CLF  groupXIII  6076868      25  5.698779e-02  30
#> 26     pop    ASP~CLF    groupXI  6177332      26  2.661382e-01  30
#> 27     pop    ASP~CLF    groupIV  6291287      27 -6.091371e-02  32
#> 28     pop    ASP~CLF    groupIX  6515796      28 -6.307368e-02  28
#> 29     pop    ASP~CLF   groupXII  6718079      29 -6.042399e-02  26
#> 30     pop    ASP~CLF   groupVII  6967477      30 -1.267996e-02  32
#> 31     pop    ASP~CLF   groupIII  7194069      31  5.963901e-02  34
#> 32     pop    ASP~CLF   groupXVI  7199722      32 -5.969459e-02  30
#> 33     pop    ASP~CLF  groupVIII  7224065      33 -2.927616e-02  30
#> 34     pop    ASP~CLF groupXVIII  7229426      34 -3.491821e-02  30
#> 35     pop    ASP~CLF groupXVIII  8087022      35  5.964215e-03  30
#> 36     pop    ASP~CLF    groupIX  8373442      36 -6.115515e-02  32
#> 37     pop    ASP~CLF    groupXV  8815538      37 -1.021070e-01  30
#> 38     pop    ASP~CLF    groupXI  9019877      38  3.327496e-01  32
#> 39     pop    ASP~CLF   groupXXI  9067876      39 -5.449591e-03  24
#> 40     pop    ASP~CLF  groupXVII  9108184      40  5.513849e-02  30
#> 41     pop    ASP~CLF    groupIX  9284427      41  1.092357e-01  32
#> 42     pop    ASP~CLF     groupV  9767769      42 -2.971250e-02  26
#> 43     pop    ASP~CLF   groupXII  9881312      43 -1.267996e-02  32
#> 44     pop    ASP~CLF    groupXX 10296516      44  3.311258e-02  30
#> 45     pop    ASP~CLF     groupI 10593829      45           NaN  34
#> 46     pop    ASP~CLF   groupXIX 10765968      46 -6.758544e-03  26
#> 47     pop    ASP~CLF    groupII 10829169      47  1.437318e-01  34
#> 48     pop    ASP~CLF    groupXX 10913115      48  3.302895e-02  30
#> 49     pop    ASP~CLF groupXVIII 10974594      49  1.628959e-01  32
#> 50     pop    ASP~CLF  groupXVII 11256661      50 -6.170094e-02  30
#> 51     pop    ASP~CLF  groupVIII 11438870      51  2.817912e-01  28
#> 52     pop    ASP~CLF   groupVII 11477223      52 -4.999555e-02  34
#> 53     pop    ASP~CLF    groupII 11505084      53 -1.793523e-02  28
#> 54     pop    ASP~CLF    groupIX 11904651      54  5.560098e-02  34
#> 55     pop    ASP~CLF  groupVIII 12108945      55  8.108767e-02  28
#> 56     pop    ASP~CLF   groupXIX 12485182      56  1.215500e-01  30
#> 57     pop    ASP~CLF    groupXI 12833733      57 -3.378573e-02  26
#> 58     pop    ASP~CLF groupXVIII 12937786      58  3.478512e-01  30
#> 59     pop    ASP~CLF   groupIII 13177027      59 -4.228420e-02  34
#> 60     pop    ASP~CLF    groupVI 13335577      60  1.345785e-01  26
#> 61     pop    ASP~CLF     groupX 13718292      61  1.895245e-01  22
#> 62     pop    ASP~CLF    groupVI 13895716      62  1.951220e-02  32
#> 63     pop    ASP~CLF   groupXIV 14009990      63 -9.742524e-02  26
#> 64     pop    ASP~CLF   groupXII 14158350      64 -2.724605e-02  34
#> 65     pop    ASP~CLF    groupIX 14164215      65  7.684590e-02  34
#> 66     pop    ASP~CLF    groupVI 14236496      66 -5.912523e-02  30
#> 67     pop    ASP~CLF     groupX 14273796      67           NaN  30
#> 68     pop    ASP~CLF  groupVIII 14445389      68 -7.655502e-02  32
#> 69     pop    ASP~CLF    groupIX 14868533      69  7.749077e-02  30
#> 70     pop    ASP~CLF    groupIX 15190693      70  6.548933e-02  30
#> 71     pop    ASP~CLF   groupXII 15306362      71 -3.173122e-02  34
#> 72     pop    ASP~CLF    groupXX 15460430      72 -8.132147e-02  32
#> 73     pop    ASP~CLF   groupVII 15493025      73 -3.973510e-02  32
#> 74     pop    ASP~CLF    groupXV 15547073      74 -6.170094e-02  30
#> 75     pop    ASP~CLF  groupXIII 15697285      75 -5.760123e-02  28
#> 76     pop    ASP~CLF   groupIII 15983994      76  0.000000e+00  18
#> 77     pop    ASP~CLF   groupXIX 16233699      77  1.215500e-01  30
#> 78     pop    ASP~CLF    groupIX 16306966      78  2.268041e-01  32
#> 79     pop    ASP~CLF     groupI 16498302      79  9.852110e-02  34
#> 80     pop    ASP~CLF  groupXIII 16520929      80 -6.838784e-02  34
#> 81     pop    ASP~CLF   groupXVI 16561891      81  1.005424e-01  28
#> 82     pop    ASP~CLF    groupII 16598364      82 -4.001296e-02  26
#> 83     pop    ASP~CLF    groupIV 17248473      83  9.725061e-02  34
#> 84     pop    ASP~CLF    groupXX 17748944      84 -8.344924e-02  32
#> 85     pop    ASP~CLF  groupXIII 17842228      85 -5.569410e-02  30
#> 86     pop    ASP~CLF   groupXIX 17858769      86  1.695502e-01  28
#> 87     pop    ASP~CLF    groupIX 17929225      87  2.630258e-01  34
#> 88     pop    ASP~CLF    groupIX 18339767      88  8.108767e-02  28
#> 89     pop    ASP~CLF  groupXIII 18836558      89 -1.201717e-01  24
#> 90     pop    ASP~CLF    groupII 19244005      90  5.422518e-02  34
#> 91     pop    ASP~CLF    groupXX 19317273      91  1.181051e-02  26
#> 92     pop    ASP~CLF     groupI 19759331      92  8.842339e-02  24
#> 93     pop    ASP~CLF   groupVII 20165240      93 -9.762100e-03  34
#> 94     pop    ASP~CLF    groupII 20496069      94           NaN  26
#> 95     pop    ASP~CLF     groupI 21284266      95 -2.488336e-02  32
#> 96     pop    ASP~CLF    groupIV 22227311      96 -7.497566e-02  28
#> 97     pop    ASP~CLF     groupI 22746027      97  5.560098e-02  34
#> 98     pop    ASP~CLF   groupVII 25239495      98 -1.123500e-01  28
#> 99     pop    ASP~CLF   groupVII 27599144      99 -5.640423e-02  32
#> 100    pop    ASP~CLF    groupIV 32321895     100  8.330824e-03  28
#> 101    pop    ASP~OPL    groupVI   212436       1 -2.644876e-02  58
#> 102    pop    ASP~OPL groupXVIII   350070       2 -3.376814e-02  66
#> 103    pop    ASP~OPL    groupXV   854285       3 -3.083328e-02  72
#> 104    pop    ASP~OPL  groupXIII   907279       4  1.610105e-01  62
#> 105    pop    ASP~OPL    groupVI  1068620       5 -5.372451e-02  40
#> 106    pop    ASP~OPL     groupV  1392753       6  1.451491e-01  72
#> 107    pop    ASP~OPL    groupIX  2075079       7  1.750041e-02  64
#> 108    pop    ASP~OPL  groupXIII  2106000       8 -2.121558e-02  70
#> 109    pop    ASP~OPL   groupXII  2136587       9  3.153783e-02  72
#> 110    pop    ASP~OPL  groupXVII  2560839      10  3.069926e-02  72
#> 111    pop    ASP~OPL   groupIII  2750695      11  6.144007e-03  62
#> 112    pop    ASP~OPL   groupXIX  3008524      12 -2.088796e-02  56
#> 113    pop    ASP~OPL   groupXVI  3624417      13  1.473334e-02  66
#> 114    pop    ASP~OPL   groupVII  3667528      14  1.936424e-02  70
#> 115    pop    ASP~OPL   groupVII  4189120      15 -3.273141e-02  62
#> 116    pop    ASP~OPL    groupXX  4538074      16 -2.176638e-02  60
#> 117    pop    ASP~OPL  groupXIII  4706089      17 -3.068047e-02  70
#> 118    pop    ASP~OPL   groupIII  4746808      18 -3.959117e-03  74
#> 119    pop    ASP~OPL  groupXVII  5041995      19  2.663816e-02  62
#> 120    pop    ASP~OPL     groupX  5161637      20 -1.813233e-02  66
#> 121    pop    ASP~OPL   groupXII  5178037      21  1.830885e-01  56
#> 122    pop    ASP~OPL    groupIX  5287711      22 -3.055093e-03  60
#> 123    pop    ASP~OPL  groupVIII  5462135      23 -3.440223e-02  70
#> 124    pop    ASP~OPL    groupIX  5960933      24 -3.231976e-02  72
#> 125    pop    ASP~OPL  groupXIII  6076868      25 -1.990941e-02  66
#> 126    pop    ASP~OPL    groupXI  6177332      26  8.059901e-02  60
#> 127    pop    ASP~OPL    groupIV  6291287      27  6.190370e-02  68
#> 128    pop    ASP~OPL    groupIX  6515796      28  4.027542e-01  54
#> 129    pop    ASP~OPL   groupXII  6718079      29  8.244069e-02  60
#> 130    pop    ASP~OPL   groupVII  6967477      30 -3.114541e-02  64
#> 131    pop    ASP~OPL   groupIII  7194069      31 -2.308622e-02  72
#> 132    pop    ASP~OPL   groupXVI  7199722      32 -3.060081e-02  66
#> 133    pop    ASP~OPL  groupVIII  7224065      33 -1.505836e-02  60
#> 134    pop    ASP~OPL groupXVIII  7229426      34 -3.049257e-02  54
#> 135    pop    ASP~OPL groupXVIII  8087022      35 -2.826493e-02  64
#> 136    pop    ASP~OPL    groupIX  8373442      36 -3.230430e-02  68
#> 137    pop    ASP~OPL    groupXV  8815538      37 -1.714656e-02  68
#> 138    pop    ASP~OPL    groupXI  9019877      38  1.166414e-01  66
#> 139    pop    ASP~OPL   groupXXI  9067876      39 -2.015113e-02  48
#> 140    pop    ASP~OPL  groupXVII  9108184      40  1.458560e-01  60
#> 141    pop    ASP~OPL    groupIX  9284427      41  1.825902e-01  60
#> 142    pop    ASP~OPL     groupV  9767769      42  4.488137e-02  62
#> 143    pop    ASP~OPL   groupXII  9881312      43  8.598402e-03  64
#> 144    pop    ASP~OPL    groupXX 10296516      44  2.671541e-02  70
#> 145    pop    ASP~OPL     groupI 10593829      45           NaN  60
#> 146    pop    ASP~OPL   groupXIX 10765968      46 -5.883266e-02  52
#> 147    pop    ASP~OPL    groupII 10829169      47  9.894342e-02  68
#> 148    pop    ASP~OPL    groupXX 10913115      48  1.266249e-01  58
#> 149    pop    ASP~OPL groupXVIII 10974594      49 -2.248356e-02  74
#> 150    pop    ASP~OPL  groupXVII 11256661      50 -2.466147e-02  66
#> 151    pop    ASP~OPL  groupVIII 11438870      51 -3.481320e-02  58
#> 152    pop    ASP~OPL   groupVII 11477223      52 -2.368449e-02  70
#> 153    pop    ASP~OPL    groupII 11505084      53 -6.083345e-03  62
#> 154    pop    ASP~OPL    groupIX 11904651      54 -1.940484e-02  74
#> 155    pop    ASP~OPL  groupVIII 12108945      55  1.557943e-02  58
#> 156    pop    ASP~OPL   groupXIX 12485182      56  2.035997e-02  58
#> 157    pop    ASP~OPL    groupXI 12833733      57 -3.360486e-02  56
#> 158    pop    ASP~OPL groupXVIII 12937786      58  4.694890e-02  68
#> 159    pop    ASP~OPL   groupIII 13177027      59 -2.367062e-02  70
#> 160    pop    ASP~OPL    groupVI 13335577      60 -2.082329e-02  60
#> 161    pop    ASP~OPL     groupX 13718292      61 -2.574469e-02  50
#> 162    pop    ASP~OPL    groupVI 13895716      62  1.441818e-02  72
#> 163    pop    ASP~OPL   groupXIV 14009990      63  2.075538e-02  52
#> 164    pop    ASP~OPL   groupXII 14158350      64 -2.540336e-02  64
#> 165    pop    ASP~OPL    groupIX 14164215      65  3.084491e-01  60
#> 166    pop    ASP~OPL    groupVI 14236496      66  1.488719e-02  58
#> 167    pop    ASP~OPL     groupX 14273796      67           NaN  68
#> 168    pop    ASP~OPL  groupVIII 14445389      68 -3.732624e-02  68
#> 169    pop    ASP~OPL    groupIX 14868533      69  1.993300e-01  70
#> 170    pop    ASP~OPL    groupIX 15190693      70 -3.798855e-03  62
#> 171    pop    ASP~OPL   groupXII 15306362      71 -1.825265e-02  76
#> 172    pop    ASP~OPL    groupXX 15460430      72  1.799863e-01  68
#> 173    pop    ASP~OPL   groupVII 15493025      73  6.661526e-02  66
#> 174    pop    ASP~OPL    groupXV 15547073      74  4.594539e-02  58
#> 175    pop    ASP~OPL  groupXIII 15697285      75 -4.477040e-04  64
#> 176    pop    ASP~OPL   groupIII 15983994      76 -7.772155e-02  44
#> 177    pop    ASP~OPL   groupXIX 16233699      77  4.538012e-02  64
#> 178    pop    ASP~OPL    groupIX 16306966      78  2.043377e-01  64
#> 179    pop    ASP~OPL     groupI 16498302      79 -9.368485e-03  74
#> 180    pop    ASP~OPL  groupXIII 16520929      80 -1.145346e-02  66
#> 181    pop    ASP~OPL   groupXVI 16561891      81 -1.687033e-02  60
#> 182    pop    ASP~OPL    groupII 16598364      82 -3.282814e-02  72
#> 183    pop    ASP~OPL    groupIV 17248473      83  2.535073e-01  74
#> 184    pop    ASP~OPL    groupXX 17748944      84 -3.237916e-02  60
#> 185    pop    ASP~OPL  groupXIII 17842228      85  1.895493e-01  62
#> 186    pop    ASP~OPL   groupXIX 17858769      86  2.886101e-02  56
#> 187    pop    ASP~OPL    groupIX 17929225      87  1.102540e-01  72
#> 188    pop    ASP~OPL    groupIX 18339767      88  2.661748e-02  56
#> 189    pop    ASP~OPL  groupXIII 18836558      89 -3.352789e-02  54
#> 190    pop    ASP~OPL    groupII 19244005      90  1.877812e-01  74
#> 191    pop    ASP~OPL    groupXX 19317273      91  1.787955e-01  60
#> 192    pop    ASP~OPL     groupI 19759331      92  2.187386e-03  42
#> 193    pop    ASP~OPL   groupVII 20165240      93  1.132886e-02  74
#> 194    pop    ASP~OPL    groupII 20496069      94  7.276497e-03  60
#> 195    pop    ASP~OPL     groupI 21284266      95 -4.488946e-02  64
#> 196    pop    ASP~OPL    groupIV 22227311      96  1.404551e-02  62
#> 197    pop    ASP~OPL     groupI 22746027      97  2.900875e-03  64
#> 198    pop    ASP~OPL   groupVII 25239495      98  2.819589e-02  70
#> 199    pop    ASP~OPL   groupVII 27599144      99 -2.884968e-02  62
#> 200    pop    ASP~OPL    groupIV 32321895     100  1.339648e-01  60
#> 201    pop    ASP~PAL    groupVI   212436       1  1.380180e-01  52
#> 202    pop    ASP~PAL groupXVIII   350070       2  2.848134e-02  54
#> 203    pop    ASP~PAL    groupXV   854285       3 -3.650796e-03  54
#> 204    pop    ASP~PAL  groupXIII   907279       4  1.332397e-01  54
#> 205    pop    ASP~PAL    groupVI  1068620       5  2.606590e-03  38
#> 206    pop    ASP~PAL     groupV  1392753       6  2.725543e-02  52
#> 207    pop    ASP~PAL    groupIX  2075079       7  2.097863e-02  56
#> 208    pop    ASP~PAL  groupXIII  2106000       8  2.272388e-01  50
#> 209    pop    ASP~PAL   groupXII  2136587       9  1.369602e-01  54
#> 210    pop    ASP~PAL  groupXVII  2560839      10  8.083761e-02  56
#> 211    pop    ASP~PAL   groupIII  2750695      11  4.948636e-01  52
#> 212    pop    ASP~PAL   groupXIX  3008524      12 -1.619716e-02  52
#> 213    pop    ASP~PAL   groupXVI  3624417      13  8.504872e-03  48
#> 214    pop    ASP~PAL   groupVII  3667528      14  2.286598e-01  54
#> 215    pop    ASP~PAL   groupVII  4189120      15  5.936104e-01  54
#> 216    pop    ASP~PAL    groupXX  4538074      16  3.156297e-02  50
#> 217    pop    ASP~PAL  groupXIII  4706089      17 -4.150376e-02  52
#> 218    pop    ASP~PAL   groupIII  4746808      18  2.676383e-01  48
#> 219    pop    ASP~PAL  groupXVII  5041995      19           NaN  56
#> 220    pop    ASP~PAL     groupX  5161637      20  3.049387e-02  48
#> 221    pop    ASP~PAL   groupXII  5178037      21  1.882905e-01  50
#> 222    pop    ASP~PAL    groupIX  5287711      22  5.137384e-01  52
#> 223    pop    ASP~PAL  groupVIII  5462135      23  2.506749e-02  52
#> 224    pop    ASP~PAL    groupIX  5960933      24  2.506749e-02  52
#> 225    pop    ASP~PAL  groupXIII  6076868      25  3.419417e-01  52
#> 226    pop    ASP~PAL    groupXI  6177332      26  4.249079e-02  48
#> 227    pop    ASP~PAL    groupIV  6291287      27  2.848134e-02  54
#> 228    pop    ASP~PAL    groupIX  6515796      28 -4.280760e-02  46
#> 229    pop    ASP~PAL   groupXII  6718079      29  1.697666e-02  48
#> 230    pop    ASP~PAL   groupVII  6967477      30  8.083761e-02  56
#> 231    pop    ASP~PAL   groupIII  7194069      31  3.062614e-01  56
#> 232    pop    ASP~PAL   groupXVI  7199722      32  1.201923e-01  50
#> 233    pop    ASP~PAL  groupVIII  7224065      33 -2.410547e-02  52
#> 234    pop    ASP~PAL groupXVIII  7229426      34  4.706156e-01  52
#> 235    pop    ASP~PAL groupXVIII  8087022      35  9.219858e-02  52
#> 236    pop    ASP~PAL    groupIX  8373442      36  1.573693e-01  52
#> 237    pop    ASP~PAL    groupXV  8815538      37  1.032312e-01  52
#> 238    pop    ASP~PAL    groupXI  9019877      38 -3.539273e-02  52
#> 239    pop    ASP~PAL   groupXXI  9067876      39 -4.235510e-02  40
#> 240    pop    ASP~PAL  groupXVII  9108184      40  1.470972e-01  50
#> 241    pop    ASP~PAL    groupIX  9284427      41  4.124992e-01  48
#> 242    pop    ASP~PAL     groupV  9767769      42  1.144353e-01  48
#> 243    pop    ASP~PAL   groupXII  9881312      43  8.083761e-02  56
#> 244    pop    ASP~PAL    groupXX 10296516      44  7.311282e-02  48
#> 245    pop    ASP~PAL     groupI 10593829      45           NaN  56
#> 246    pop    ASP~PAL   groupXIX 10765968      46  4.995400e-01  48
#> 247    pop    ASP~PAL    groupII 10829169      47  2.899037e-02  52
#> 248    pop    ASP~PAL    groupXX 10913115      48 -3.068139e-02  46
#> 249    pop    ASP~PAL groupXVIII 10974594      49           NaN  54
#> 250    pop    ASP~PAL  groupXVII 11256661      50  6.006774e-02  52
#> 251    pop    ASP~PAL  groupVIII 11438870      51  6.129905e-02  40
#> 252    pop    ASP~PAL   groupVII 11477223      52 -3.286702e-03  46
#> 253    pop    ASP~PAL    groupII 11505084      53  2.000985e-01  52
#> 254    pop    ASP~PAL    groupIX 11904651      54           NaN  54
#> 255    pop    ASP~PAL  groupVIII 12108945      55 -2.952449e-02  46
#> 256    pop    ASP~PAL   groupXIX 12485182      56 -2.265904e-02  50
#> 257    pop    ASP~PAL    groupXI 12833733      57  2.711961e-02  40
#> 258    pop    ASP~PAL groupXVIII 12937786      58 -5.942235e-02  44
#> 259    pop    ASP~PAL   groupIII 13177027      59  3.306167e-01  54
#> 260    pop    ASP~PAL    groupVI 13335577      60  5.443417e-01  48
#> 261    pop    ASP~PAL     groupX 13718292      61           NaN  42
#> 262    pop    ASP~PAL    groupVI 13895716      62  9.219858e-02  52
#> 263    pop    ASP~PAL   groupXIV 14009990      63  2.822391e-01  36
#> 264    pop    ASP~PAL   groupXII 14158350      64  8.033976e-02  50
#> 265    pop    ASP~PAL    groupIX 14164215      65 -1.605372e-02  56
#> 266    pop    ASP~PAL    groupVI 14236496      66  5.759466e-03  50
#> 267    pop    ASP~PAL     groupX 14273796      67  1.890808e-01  48
#> 268    pop    ASP~PAL  groupVIII 14445389      68  3.418201e-01  46
#> 269    pop    ASP~PAL    groupIX 14868533      69           NaN  54
#> 270    pop    ASP~PAL    groupIX 15190693      70  4.000707e-01  52
#> 271    pop    ASP~PAL   groupXII 15306362      71  5.524675e-02  54
#> 272    pop    ASP~PAL    groupXX 15460430      72 -5.072464e-02  50
#> 273    pop    ASP~PAL   groupVII 15493025      73 -3.455425e-02  50
#> 274    pop    ASP~PAL    groupXV 15547073      74  1.641935e-01  54
#> 275    pop    ASP~PAL  groupXIII 15697285      75  9.111469e-02  50
#> 276    pop    ASP~PAL   groupIII 15983994      76  5.892487e-01  28
#> 277    pop    ASP~PAL   groupXIX 16233699      77 -1.935852e-02  48
#> 278    pop    ASP~PAL    groupIX 16306966      78           NaN  54
#> 279    pop    ASP~PAL     groupI 16498302      79 -2.043572e-02  56
#> 280    pop    ASP~PAL  groupXIII 16520929      80  1.175007e-01  52
#> 281    pop    ASP~PAL   groupXVI 16561891      81 -7.786513e-03  50
#> 282    pop    ASP~PAL    groupII 16598364      82  9.750624e-02  54
#> 283    pop    ASP~PAL    groupIV 17248473      83  1.891374e-01  54
#> 284    pop    ASP~PAL    groupXX 17748944      84  7.329220e-01  50
#> 285    pop    ASP~PAL  groupXIII 17842228      85  1.289552e-01  48
#> 286    pop    ASP~PAL   groupXIX 17858769      86 -1.839488e-02  38
#> 287    pop    ASP~PAL    groupIX 17929225      87  5.311023e-01  54
#> 288    pop    ASP~PAL    groupIX 18339767      88  5.680557e-01  44
#> 289    pop    ASP~PAL  groupXIII 18836558      89 -5.809129e-03  42
#> 290    pop    ASP~PAL    groupII 19244005      90  1.392444e-01  56
#> 291    pop    ASP~PAL    groupXX 19317273      91 -1.432482e-02  50
#> 292    pop    ASP~PAL     groupI 19759331      92 -2.283582e-02  36
#> 293    pop    ASP~PAL   groupVII 20165240      93 -3.676579e-02  50
#> 294    pop    ASP~PAL    groupII 20496069      94  2.683315e-01  42
#> 295    pop    ASP~PAL     groupI 21284266      95  2.286598e-01  54
#> 296    pop    ASP~PAL    groupIV 22227311      96  5.240276e-02  50
#> 297    pop    ASP~PAL     groupI 22746027      97           NaN  56
#> 298    pop    ASP~PAL   groupVII 25239495      98  2.954025e-02  50
#> 299    pop    ASP~PAL   groupVII 27599144      99  5.761323e-01  48
#> 300    pop    ASP~PAL    groupIV 32321895     100  1.007767e-01  42
#> 301    pop    ASP~SMR    groupVI   212436       1 -3.055961e-02  62
#> 302    pop    ASP~SMR groupXVIII   350070       2 -3.682745e-02  64
#> 303    pop    ASP~SMR    groupXV   854285       3 -2.790522e-02  68
#> 304    pop    ASP~SMR  groupXIII   907279       4  5.377446e-02  68
#> 305    pop    ASP~SMR    groupVI  1068620       5 -2.154820e-02  44
#> 306    pop    ASP~SMR     groupV  1392753       6  1.051405e-01  62
#> 307    pop    ASP~SMR    groupIX  2075079       7  5.647609e-02  64
#> 308    pop    ASP~SMR  groupXIII  2106000       8  4.138290e-03  62
#> 309    pop    ASP~SMR   groupXII  2136587       9 -3.086545e-02  62
#> 310    pop    ASP~SMR  groupXVII  2560839      10 -1.512214e-02  66
#> 311    pop    ASP~SMR   groupIII  2750695      11 -3.152378e-02  58
#> 312    pop    ASP~SMR   groupXIX  3008524      12 -1.448195e-02  54
#> 313    pop    ASP~SMR   groupXVI  3624417      13 -1.583938e-02  58
#> 314    pop    ASP~SMR   groupVII  3667528      14 -2.881863e-02  64
#> 315    pop    ASP~SMR   groupVII  4189120      15 -2.614829e-02  64
#> 316    pop    ASP~SMR    groupXX  4538074      16 -2.176638e-02  60
#> 317    pop    ASP~SMR  groupXIII  4706089      17 -3.090728e-02  62
#> 318    pop    ASP~SMR   groupIII  4746808      18  3.399574e-02  66
#> 319    pop    ASP~SMR  groupXVII  5041995      19  1.752962e-02  68
#> 320    pop    ASP~SMR     groupX  5161637      20 -3.791346e-02  64
#> 321    pop    ASP~SMR   groupXII  5178037      21  3.709081e-02  58
#> 322    pop    ASP~SMR    groupIX  5287711      22  1.767216e-01  58
#> 323    pop    ASP~SMR  groupVIII  5462135      23 -1.947863e-02  60
#> 324    pop    ASP~SMR    groupIX  5960933      24 -3.881905e-02  58
#> 325    pop    ASP~SMR  groupXIII  6076868      25  4.951514e-02  58
#> 326    pop    ASP~SMR    groupXI  6177332      26 -1.990307e-02  62
#> 327    pop    ASP~SMR    groupIV  6291287      27  1.502858e-02  60
#> 328    pop    ASP~SMR    groupIX  6515796      28  2.474559e-01  56
#> 329    pop    ASP~SMR   groupXII  6718079      29  6.205223e-02  48
#> 330    pop    ASP~SMR   groupVII  6967477      30 -2.923221e-02  66
#> 331    pop    ASP~SMR   groupIII  7194069      31 -9.335851e-03  60
#> 332    pop    ASP~SMR   groupXVI  7199722      32 -2.988160e-02  60
#> 333    pop    ASP~SMR  groupVIII  7224065      33 -2.852615e-02  66
#> 334    pop    ASP~SMR groupXVIII  7229426      34 -1.968875e-02  62
#> 335    pop    ASP~SMR groupXVIII  8087022      35 -3.234620e-02  60
#> 336    pop    ASP~SMR    groupIX  8373442      36 -2.712438e-02  64
#> 337    pop    ASP~SMR    groupXV  8815538      37  2.007391e-02  64
#> 338    pop    ASP~SMR    groupXI  9019877      38 -3.708065e-02  58
#> 339    pop    ASP~SMR   groupXXI  9067876      39 -8.522660e-03  52
#> 340    pop    ASP~SMR  groupXVII  9108184      40  1.692207e-01  58
#> 341    pop    ASP~SMR    groupIX  9284427      41  1.097465e-01  56
#> 342    pop    ASP~SMR     groupV  9767769      42  1.579337e-01  60
#> 343    pop    ASP~SMR   groupXII  9881312      43 -3.035462e-02  70
#> 344    pop    ASP~SMR    groupXX 10296516      44 -5.054177e-02  60
#> 345    pop    ASP~SMR     groupI 10593829      45           NaN  62
#> 346    pop    ASP~SMR   groupXIX 10765968      46  5.514256e-02  50
#> 347    pop    ASP~SMR    groupII 10829169      47 -3.613596e-02  68
#> 348    pop    ASP~SMR    groupXX 10913115      48  1.592115e-02  54
#> 349    pop    ASP~SMR groupXVIII 10974594      49           NaN  68
#> 350    pop    ASP~SMR  groupXVII 11256661      50  2.455825e-02  58
#> 351    pop    ASP~SMR  groupVIII 11438870      51 -3.187121e-02  46
#> 352    pop    ASP~SMR   groupVII 11477223      52  1.055515e-01  70
#> 353    pop    ASP~SMR    groupII 11505084      53 -4.313915e-02  66
#> 354    pop    ASP~SMR    groupIX 11904651      54 -1.825776e-02  68
#> 355    pop    ASP~SMR  groupVIII 12108945      55  2.018139e-01  56
#> 356    pop    ASP~SMR   groupXIX 12485182      56  3.297867e-02  56
#> 357    pop    ASP~SMR    groupXI 12833733      57 -2.989004e-02  52
#> 358    pop    ASP~SMR groupXVIII 12937786      58  1.904408e-01  54
#> 359    pop    ASP~SMR   groupIII 13177027      59 -3.304979e-02  62
#> 360    pop    ASP~SMR    groupVI 13335577      60 -2.618844e-02  56
#> 361    pop    ASP~SMR     groupX 13718292      61 -7.576503e-03  56
#> 362    pop    ASP~SMR    groupVI 13895716      62 -2.399451e-02  68
#> 363    pop    ASP~SMR   groupXIV 14009990      63 -3.555733e-02  48
#> 364    pop    ASP~SMR   groupXII 14158350      64 -3.160370e-02  66
#> 365    pop    ASP~SMR    groupIX 14164215      65  2.524323e-02  66
#> 366    pop    ASP~SMR    groupVI 14236496      66 -1.879288e-02  52
#> 367    pop    ASP~SMR     groupX 14273796      67           NaN  62
#> 368    pop    ASP~SMR  groupVIII 14445389      68 -3.986270e-02  56
#> 369    pop    ASP~SMR    groupIX 14868533      69  3.950332e-02  62
#> 370    pop    ASP~SMR    groupIX 15190693      70  1.927747e-03  64
#> 371    pop    ASP~SMR   groupXII 15306362      71  4.893582e-03  70
#> 372    pop    ASP~SMR    groupXX 15460430      72  2.938878e-01  64
#> 373    pop    ASP~SMR   groupVII 15493025      73  8.031009e-02  62
#> 374    pop    ASP~SMR    groupXV 15547073      74 -8.506282e-03  54
#> 375    pop    ASP~SMR  groupXIII 15697285      75 -2.229121e-02  58
#> 376    pop    ASP~SMR   groupIII 15983994      76  1.668981e-01  36
#> 377    pop    ASP~SMR   groupXIX 16233699      77  3.297867e-02  56
#> 378    pop    ASP~SMR    groupIX 16306966      78  5.459426e-02  64
#> 379    pop    ASP~SMR     groupI 16498302      79 -1.537988e-03  66
#> 380    pop    ASP~SMR  groupXIII 16520929      80 -2.634475e-02  68
#> 381    pop    ASP~SMR   groupXVI 16561891      81 -2.767784e-02  62
#> 382    pop    ASP~SMR    groupII 16598364      82 -1.971113e-02  58
#> 383    pop    ASP~SMR    groupIV 17248473      83  1.524508e-01  70
#> 384    pop    ASP~SMR    groupXX 17748944      84 -3.606714e-02  58
#> 385    pop    ASP~SMR  groupXIII 17842228      85  2.442944e-01  62
#> 386    pop    ASP~SMR   groupXIX 17858769      86  6.722066e-02  54
#> 387    pop    ASP~SMR    groupIX 17929225      87  8.113347e-02  66
#> 388    pop    ASP~SMR    groupIX 18339767      88 -9.942970e-03  48
#> 389    pop    ASP~SMR  groupXIII 18836558      89 -3.352789e-02  54
#> 390    pop    ASP~SMR    groupII 19244005      90  8.363214e-02  68
#> 391    pop    ASP~SMR    groupXX 19317273      91  2.261224e-01  58
#> 392    pop    ASP~SMR     groupI 19759331      92  1.623520e-01  36
#> 393    pop    ASP~SMR   groupVII 20165240      93 -2.217297e-02  60
#> 394    pop    ASP~SMR    groupII 20496069      94 -2.856392e-02  60
#> 395    pop    ASP~SMR     groupI 21284266      95 -2.599566e-02  58
#> 396    pop    ASP~SMR    groupIV 22227311      96 -3.716582e-02  62
#> 397    pop    ASP~SMR     groupI 22746027      97  4.928458e-03  62
#> 398    pop    ASP~SMR   groupVII 25239495      98 -3.814923e-02  62
#> 399    pop    ASP~SMR   groupVII 27599144      99 -4.155723e-02  60
#> 400    pop    ASP~SMR    groupIV 32321895     100  1.177925e-01  60
#> 401    pop    ASP~UPD    groupVI   212436       1 -5.651565e-02  44
#> 402    pop    ASP~UPD groupXVIII   350070       2  2.150538e-02  50
#> 403    pop    ASP~UPD    groupXV   854285       3 -3.476799e-02  52
#> 404    pop    ASP~UPD  groupXIII   907279       4  1.562174e-01  48
#> 405    pop    ASP~UPD    groupVI  1068620       5 -6.197070e-02  32
#> 406    pop    ASP~UPD     groupV  1392753       6 -6.554640e-03  48
#> 407    pop    ASP~UPD    groupIX  2075079       7  1.148516e-02  50
#> 408    pop    ASP~UPD  groupXIII  2106000       8 -3.474590e-02  42
#> 409    pop    ASP~UPD   groupXII  2136587       9 -4.585797e-02  48
#> 410    pop    ASP~UPD  groupXVII  2560839      10 -2.461878e-02  50
#> 411    pop    ASP~UPD   groupIII  2750695      11 -2.230625e-02  48
#> 412    pop    ASP~UPD   groupXIX  3008524      12 -2.777778e-02  40
#> 413    pop    ASP~UPD   groupXVI  3624417      13 -3.531233e-02  44
#> 414    pop    ASP~UPD   groupVII  3667528      14  3.340013e-02  50
#> 415    pop    ASP~UPD   groupVII  4189120      15  9.093247e-02  48
#> 416    pop    ASP~UPD    groupXX  4538074      16 -5.638445e-02  44
#> 417    pop    ASP~UPD  groupXIII  4706089      17 -3.971392e-02  46
#> 418    pop    ASP~UPD   groupIII  4746808      18 -8.763977e-03  50
#> 419    pop    ASP~UPD  groupXVII  5041995      19  2.099626e-02  52
#> 420    pop    ASP~UPD     groupX  5161637      20 -4.437853e-02  46
#> 421    pop    ASP~UPD   groupXII  5178037      21  2.677330e-02  42
#> 422    pop    ASP~UPD    groupIX  5287711      22 -2.828237e-02  46
#> 423    pop    ASP~UPD  groupVIII  5462135      23 -3.962264e-02  50
#> 424    pop    ASP~UPD    groupIX  5960933      24 -1.996232e-02  48
#> 425    pop    ASP~UPD  groupXIII  6076868      25 -1.208459e-03  50
#> 426    pop    ASP~UPD    groupXI  6177332      26 -2.533035e-02  42
#> 427    pop    ASP~UPD    groupIV  6291287      27 -3.950803e-02  48
#> 428    pop    ASP~UPD    groupIX  6515796      28  3.769117e-02  44
#> 429    pop    ASP~UPD   groupXII  6718079      29 -2.357133e-02  44
#> 430    pop    ASP~UPD   groupVII  6967477      30 -1.387292e-02  46
#> 431    pop    ASP~UPD   groupIII  7194069      31  1.809686e-02  52
#> 432    pop    ASP~UPD   groupXVI  7199722      32  8.603927e-03  46
#> 433    pop    ASP~UPD  groupVIII  7224065      33 -3.825216e-02  48
#> 434    pop    ASP~UPD groupXVIII  7229426      34  1.308837e-01  40
#> 435    pop    ASP~UPD groupXVIII  8087022      35 -3.685891e-02  48
#> 436    pop    ASP~UPD    groupIX  8373442      36 -1.348276e-02  46
#> 437    pop    ASP~UPD    groupXV  8815538      37 -3.081270e-02  46
#> 438    pop    ASP~UPD    groupXI  9019877      38  2.573787e-02  48
#> 439    pop    ASP~UPD   groupXXI  9067876      39 -4.680257e-02  38
#> 440    pop    ASP~UPD  groupXVII  9108184      40 -5.053949e-02  42
#> 441    pop    ASP~UPD    groupIX  9284427      41 -4.797931e-02  42
#> 442    pop    ASP~UPD     groupV  9767769      42 -1.220866e-02  44
#> 443    pop    ASP~UPD   groupXII  9881312      43  3.574309e-02  48
#> 444    pop    ASP~UPD    groupXX 10296516      44 -5.166488e-02  42
#> 445    pop    ASP~UPD     groupI 10593829      45  2.099626e-02  52
#> 446    pop    ASP~UPD   groupXIX 10765968      46 -5.887443e-02  42
#> 447    pop    ASP~UPD    groupII 10829169      47 -3.061228e-02  52
#> 448    pop    ASP~UPD    groupXX 10913115      48 -1.899344e-02  44
#> 449    pop    ASP~UPD groupXVIII 10974594      49  2.731844e-02  46
#> 450    pop    ASP~UPD  groupXVII 11256661      50 -6.629591e-03  48
#> 451    pop    ASP~UPD  groupVIII 11438870      51 -4.528235e-02  46
#> 452    pop    ASP~UPD   groupVII 11477223      52 -3.453697e-02  50
#> 453    pop    ASP~UPD    groupII 11505084      53  5.220178e-02  44
#> 454    pop    ASP~UPD    groupIX 11904651      54           NaN  52
#> 455    pop    ASP~UPD  groupVIII 12108945      55 -4.460788e-02  42
#> 456    pop    ASP~UPD   groupXIX 12485182      56 -3.960646e-03  48
#> 457    pop    ASP~UPD    groupXI 12833733      57  1.483564e-02  36
#> 458    pop    ASP~UPD groupXVIII 12937786      58  5.334968e-03  40
#> 459    pop    ASP~UPD   groupIII 13177027      59 -4.159848e-02  48
#> 460    pop    ASP~UPD    groupVI 13335577      60 -1.216216e-02  42
#> 461    pop    ASP~UPD     groupX 13718292      61  2.156017e-02  40
#> 462    pop    ASP~UPD    groupVI 13895716      62 -4.107238e-02  48
#> 463    pop    ASP~UPD   groupXIV 14009990      63 -2.779762e-02  42
#> 464    pop    ASP~UPD   groupXII 14158350      64 -2.912437e-02  52
#> 465    pop    ASP~UPD    groupIX 14164215      65  5.719799e-02  50
#> 466    pop    ASP~UPD    groupVI 14236496      66 -4.330137e-02  44
#> 467    pop    ASP~UPD     groupX 14273796      67 -1.530207e-02  44
#> 468    pop    ASP~UPD  groupVIII 14445389      68 -2.447217e-02  48
#> 469    pop    ASP~UPD    groupIX 14868533      69           NaN  50
#> 470    pop    ASP~UPD    groupIX 15190693      70 -1.858476e-02  44
#> 471    pop    ASP~UPD   groupXII 15306362      71 -3.571120e-02  52
#> 472    pop    ASP~UPD    groupXX 15460430      72 -3.078876e-02  46
#> 473    pop    ASP~UPD   groupVII 15493025      73 -3.256412e-02  48
#> 474    pop    ASP~UPD    groupXV 15547073      74 -3.557225e-02  46
#> 475    pop    ASP~UPD  groupXIII 15697285      75  1.658078e-02  48
#> 476    pop    ASP~UPD   groupIII 15983994      76 -7.188384e-02  40
#> 477    pop    ASP~UPD   groupXIX 16233699      77  6.052414e-03  46
#> 478    pop    ASP~UPD    groupIX 16306966      78           NaN  50
#> 479    pop    ASP~UPD     groupI 16498302      79 -2.587880e-02  46
#> 480    pop    ASP~UPD  groupXIII 16520929      80 -3.570423e-02  48
#> 481    pop    ASP~UPD   groupXVI 16561891      81 -2.198651e-02  46
#> 482    pop    ASP~UPD    groupII 16598364      82 -5.072295e-03  46
#> 483    pop    ASP~UPD    groupIV 17248473      83  7.977563e-02  50
#> 484    pop    ASP~UPD    groupXX 17748944      84  1.319861e-02  46
#> 485    pop    ASP~UPD  groupXIII 17842228      85  5.570953e-02  50
#> 486    pop    ASP~UPD   groupXIX 17858769      86 -3.194300e-03  42
#> 487    pop    ASP~UPD    groupIX 17929225      87 -3.119333e-02  50
#> 488    pop    ASP~UPD    groupIX 18339767      88  2.890173e-02  40
#> 489    pop    ASP~UPD  groupXIII 18836558      89  7.122694e-02  42
#> 490    pop    ASP~UPD    groupII 19244005      90 -2.046878e-02  48
#> 491    pop    ASP~UPD    groupXX 19317273      91  4.021480e-02  46
#> 492    pop    ASP~UPD     groupI 19759331      92  1.326170e-01  34
#> 493    pop    ASP~UPD   groupVII 20165240      93  8.846148e-03  50
#> 494    pop    ASP~UPD    groupII 20496069      94  4.286042e-02  44
#> 495    pop    ASP~UPD     groupI 21284266      95 -3.523669e-02  44
#> 496    pop    ASP~UPD    groupIV 22227311      96 -4.995615e-02  42
#> 497    pop    ASP~UPD     groupI 22746027      97 -1.100203e-02  52
#> 498    pop    ASP~UPD   groupVII 25239495      98 -1.363384e-02  44
#> 499    pop    ASP~UPD   groupVII 27599144      99 -4.975510e-02  42
#> 500    pop    ASP~UPD    groupIV 32321895     100 -4.865647e-02  48
#> 501    pop    CLF~OPL    groupVI   212436       1 -1.125113e-02  50
#> 502    pop    CLF~OPL groupXVIII   350070       2 -8.292015e-03  58
#> 503    pop    CLF~OPL    groupXV   854285       3 -5.127545e-02  62
#> 504    pop    CLF~OPL  groupXIII   907279       4  1.823409e-01  50
#> 505    pop    CLF~OPL    groupVI  1068620       5  1.348370e-01  32
#> 506    pop    CLF~OPL     groupV  1392753       6  1.549377e-01  64
#> 507    pop    CLF~OPL    groupIX  2075079       7 -6.361556e-02  54
#> 508    pop    CLF~OPL  groupXIII  2106000       8  1.161945e-01  60
#> 509    pop    CLF~OPL   groupXII  2136587       9 -2.531270e-02  62
#> 510    pop    CLF~OPL  groupXVII  2560839      10  1.275269e-02  62
#> 511    pop    CLF~OPL   groupIII  2750695      11 -3.079070e-02  54
#> 512    pop    CLF~OPL   groupXIX  3008524      12 -1.988072e-03  48
#> 513    pop    CLF~OPL   groupXVI  3624417      13  9.380863e-03  60
#> 514    pop    CLF~OPL   groupVII  3667528      14 -4.580117e-02  62
#> 515    pop    CLF~OPL   groupVII  4189120      15  8.751304e-02  52
#> 516    pop    CLF~OPL    groupXX  4538074      16 -4.773997e-02  50
#> 517    pop    CLF~OPL  groupXIII  4706089      17  5.047365e-02  62
#> 518    pop    CLF~OPL   groupIII  4746808      18  5.316289e-02  64
#> 519    pop    CLF~OPL  groupXVII  5041995      19  9.710716e-02  52
#> 520    pop    CLF~OPL     groupX  5161637      20  8.703015e-03  60
#> 521    pop    CLF~OPL   groupXII  5178037      21 -1.739255e-02  46
#> 522    pop    CLF~OPL    groupIX  5287711      22 -5.095599e-02  52
#> 523    pop    CLF~OPL  groupVIII  5462135      23 -2.693503e-02  62
#> 524    pop    CLF~OPL    groupIX  5960933      24 -7.351806e-03  64
#> 525    pop    CLF~OPL  groupXIII  6076868      25  2.868659e-02  56
#> 526    pop    CLF~OPL    groupXI  6177332      26  5.395095e-02  54
#> 527    pop    CLF~OPL    groupIV  6291287      27  1.055024e-02  60
#> 528    pop    CLF~OPL    groupIX  6515796      28  3.397651e-01  50
#> 529    pop    CLF~OPL   groupXII  6718079      29 -3.732072e-02  50
#> 530    pop    CLF~OPL   groupVII  6967477      30  3.747478e-02  52
#> 531    pop    CLF~OPL   groupIII  7194069      31  8.647607e-02  62
#> 532    pop    CLF~OPL   groupXVI  7199722      32 -3.439133e-02  56
#> 533    pop    CLF~OPL  groupVIII  7224065      33 -8.343857e-02  46
#> 534    pop    CLF~OPL groupXVIII  7229426      34 -3.227666e-02  48
#> 535    pop    CLF~OPL groupXVIII  8087022      35 -1.619057e-02  54
#> 536    pop    CLF~OPL    groupIX  8373442      36 -5.525124e-02  60
#> 537    pop    CLF~OPL    groupXV  8815538      37 -2.893487e-02  58
#> 538    pop    CLF~OPL    groupXI  9019877      38  9.717161e-02  58
#> 539    pop    CLF~OPL   groupXXI  9067876      39  1.306584e-01  36
#> 540    pop    CLF~OPL  groupXVII  9108184      40 -5.172414e-02  50
#> 541    pop    CLF~OPL    groupIX  9284427      41 -4.021332e-02  48
#> 542    pop    CLF~OPL     groupV  9767769      42 -6.436382e-02  56
#> 543    pop    CLF~OPL   groupXII  9881312      43 -5.178319e-02  52
#> 544    pop    CLF~OPL    groupXX 10296516      44 -1.343846e-02  64
#> 545    pop    CLF~OPL     groupI 10593829      45           NaN  50
#> 546    pop    CLF~OPL   groupXIX 10765968      46  3.576055e-02  46
#> 547    pop    CLF~OPL    groupII 10829169      47  3.416168e-01  58
#> 548    pop    CLF~OPL    groupXX 10913115      48 -3.710870e-02  52
#> 549    pop    CLF~OPL groupXVIII 10974594      49  1.805513e-01  66
#> 550    pop    CLF~OPL  groupXVII 11256661      50 -5.865680e-02  56
#> 551    pop    CLF~OPL  groupVIII 11438870      51  2.094502e-01  50
#> 552    pop    CLF~OPL   groupVII 11477223      52 -5.404483e-02  60
#> 553    pop    CLF~OPL    groupII 11505084      53 -4.833496e-02  50
#> 554    pop    CLF~OPL    groupIX 11904651      54  1.673456e-02  64
#> 555    pop    CLF~OPL  groupVIII 12108945      55 -4.469861e-02  54
#> 556    pop    CLF~OPL   groupXIX 12485182      56  3.078898e-03  52
#> 557    pop    CLF~OPL    groupXI 12833733      57  8.645533e-03  50
#> 558    pop    CLF~OPL groupXVIII 12937786      58  1.511259e-01  62
#> 559    pop    CLF~OPL   groupIII 13177027      59 -2.851026e-02  60
#> 560    pop    CLF~OPL    groupVI 13335577      60  2.450636e-01  58
#> 561    pop    CLF~OPL     groupX 13718292      61  9.624992e-02  40
#> 562    pop    CLF~OPL    groupVI 13895716      62  1.088139e-01  64
#> 563    pop    CLF~OPL   groupXIV 14009990      63  2.394575e-02  46
#> 564    pop    CLF~OPL   groupXII 14158350      64  6.484988e-03  54
#> 565    pop    CLF~OPL    groupIX 14164215      65  3.215433e-02  50
#> 566    pop    CLF~OPL    groupVI 14236496      66 -3.460570e-02  52
#> 567    pop    CLF~OPL     groupX 14273796      67           NaN  62
#> 568    pop    CLF~OPL  groupVIII 14445389      68 -4.878500e-02  60
#> 569    pop    CLF~OPL    groupIX 14868533      69  3.082672e-02  60
#> 570    pop    CLF~OPL    groupIX 15190693      70  1.831090e-01  52
#> 571    pop    CLF~OPL   groupXII 15306362      71 -4.207938e-02  66
#> 572    pop    CLF~OPL    groupXX 15460430      72  1.062387e-01  60
#> 573    pop    CLF~OPL   groupVII 15493025      73  1.026056e-02  58
#> 574    pop    CLF~OPL    groupXV 15547073      74  7.499092e-02  48
#> 575    pop    CLF~OPL  groupXIII 15697285      75 -4.226495e-02  56
#> 576    pop    CLF~OPL   groupIII 15983994      76 -9.261921e-02  38
#> 577    pop    CLF~OPL   groupXIX 16233699      77 -1.353566e-02  58
#> 578    pop    CLF~OPL    groupIX 16306966      78 -6.004586e-02  56
#> 579    pop    CLF~OPL     groupI 16498302      79  3.080154e-02  64
#> 580    pop    CLF~OPL  groupXIII 16520929      80  1.382196e-02  56
#> 581    pop    CLF~OPL   groupXVI 16561891      81  3.396994e-02  52
#> 582    pop    CLF~OPL    groupII 16598364      82 -2.978852e-02  58
#> 583    pop    CLF~OPL    groupIV 17248473      83  6.517270e-01  64
#> 584    pop    CLF~OPL    groupXX 17748944      84 -5.783207e-02  52
#> 585    pop    CLF~OPL  groupXIII 17842228      85  1.346862e-01  52
#> 586    pop    CLF~OPL   groupXIX 17858769      86  2.194043e-02  52
#> 587    pop    CLF~OPL    groupIX 17929225      87  1.412771e-02  62
#> 588    pop    CLF~OPL    groupIX 18339767      88 -4.031628e-02  52
#> 589    pop    CLF~OPL  groupXIII 18836558      89 -1.196934e-01  42
#> 590    pop    CLF~OPL    groupII 19244005      90           NaN  64
#> 591    pop    CLF~OPL    groupXX 19317273      91 -3.872480e-02  50
#> 592    pop    CLF~OPL     groupI 19759331      92  1.915934e-02  38
#> 593    pop    CLF~OPL   groupVII 20165240      93  1.549377e-01  64
#> 594    pop    CLF~OPL    groupII 20496069      94 -1.619057e-02  54
#> 595    pop    CLF~OPL     groupI 21284266      95 -6.936493e-03  56
#> 596    pop    CLF~OPL    groupIV 22227311      96 -1.353566e-02  58
#> 597    pop    CLF~OPL     groupI 22746027      97 -4.217984e-02  54
#> 598    pop    CLF~OPL   groupVII 25239495      98  1.074453e-01  62
#> 599    pop    CLF~OPL   groupVII 27599144      99 -9.646663e-04  54
#> 600    pop    CLF~OPL    groupIV 32321895     100 -6.423344e-02  48
#> 601    pop    CLF~PAL    groupVI   212436       1  3.767321e-01  44
#> 602    pop    CLF~PAL groupXVIII   350070       2           NaN  46
#> 603    pop    CLF~PAL    groupXV   854285       3 -2.701329e-02  44
#> 604    pop    CLF~PAL  groupXIII   907279       4  1.389056e-01  42
#> 605    pop    CLF~PAL    groupVI  1068620       5  6.702413e-02  30
#> 606    pop    CLF~PAL     groupV  1392753       6  1.562281e-02  44
#> 607    pop    CLF~PAL    groupIX  2075079       7  1.917591e-01  46
#> 608    pop    CLF~PAL  groupXIII  2106000       8  5.165719e-01  40
#> 609    pop    CLF~PAL   groupXII  2136587       9  8.650549e-02  44
#> 610    pop    CLF~PAL  groupXVII  2560839      10  1.016993e-01  46
#> 611    pop    CLF~PAL   groupIII  2750695      11  6.341733e-01  44
#> 612    pop    CLF~PAL   groupXIX  3008524      12  9.750201e-02  44
#> 613    pop    CLF~PAL   groupXVI  3624417      13  4.105650e-01  42
#> 614    pop    CLF~PAL   groupVII  3667528      14  2.433904e-01  46
#> 615    pop    CLF~PAL   groupVII  4189120      15  8.202386e-01  44
#> 616    pop    CLF~PAL    groupXX  4538074      16 -1.630158e-02  40
#> 617    pop    CLF~PAL  groupXIII  4706089      17  7.375307e-02  44
#> 618    pop    CLF~PAL   groupIII  4746808      18  1.420691e-01  38
#> 619    pop    CLF~PAL  groupXVII  5041995      19  3.668227e-01  46
#> 620    pop    CLF~PAL     groupX  5161637      20  8.750000e-02  42
#> 621    pop    CLF~PAL   groupXII  5178037      21  5.890564e-01  40
#> 622    pop    CLF~PAL    groupIX  5287711      22  5.666051e-01  44
#> 623    pop    CLF~PAL  groupVIII  5462135      23           NaN  44
#> 624    pop    CLF~PAL    groupIX  5960933      24           NaN  44
#> 625    pop    CLF~PAL  groupXIII  6076868      25  1.389056e-01  42
#> 626    pop    CLF~PAL    groupXI  6177332      26  4.137116e-01  42
#> 627    pop    CLF~PAL    groupIV  6291287      27  1.016993e-01  46
#> 628    pop    CLF~PAL    groupIX  6515796      28 -3.337306e-02  42
#> 629    pop    CLF~PAL   groupXII  6718079      29 -6.627743e-02  38
#> 630    pop    CLF~PAL   groupVII  6967477      30  3.237155e-01  44
#> 631    pop    CLF~PAL   groupIII  7194069      31  1.016993e-01  46
#> 632    pop    CLF~PAL   groupXVI  7199722      32  5.594406e-02  40
#> 633    pop    CLF~PAL  groupVIII  7224065      33 -2.051749e-03  38
#> 634    pop    CLF~PAL groupXVIII  7229426      34  6.330049e-01  46
#> 635    pop    CLF~PAL groupXVIII  8087022      35           NaN  42
#> 636    pop    CLF~PAL    groupIX  8373442      36  2.331797e-01  44
#> 637    pop    CLF~PAL    groupXV  8815538      37  1.417854e-01  42
#> 638    pop    CLF~PAL    groupXI  9019877      38  2.967814e-01  44
#> 639    pop    CLF~PAL   groupXXI  9067876      39  4.381461e-02  28
#> 640    pop    CLF~PAL  groupXVII  9108184      40 -5.496646e-02  40
#> 641    pop    CLF~PAL    groupIX  9284427      41  8.004869e-01  36
#> 642    pop    CLF~PAL     groupV  9767769      42  3.135426e-01  42
#> 643    pop    CLF~PAL   groupXII  9881312      43  3.237155e-01  44
#> 644    pop    CLF~PAL    groupXX 10296516      44  2.811245e-01  42
#> 645    pop    CLF~PAL     groupI 10593829      45           NaN  46
#> 646    pop    CLF~PAL   groupXIX 10765968      46  3.116406e-01  42
#> 647    pop    CLF~PAL    groupII 10829169      47  1.917808e-02  42
#> 648    pop    CLF~PAL    groupXX 10913115      48 -4.647665e-02  40
#> 649    pop    CLF~PAL groupXVIII 10974594      49  2.433904e-01  46
#> 650    pop    CLF~PAL  groupXVII 11256661      50 -1.384488e-02  42
#> 651    pop    CLF~PAL  groupVIII 11438870      51  8.413672e-02  32
#> 652    pop    CLF~PAL   groupVII 11477223      52 -1.786108e-02  36
#> 653    pop    CLF~PAL    groupII 11505084      53  5.890564e-01  40
#> 654    pop    CLF~PAL    groupIX 11904651      54  9.471095e-02  44
#> 655    pop    CLF~PAL  groupVIII 12108945      55  5.373947e-02  42
#> 656    pop    CLF~PAL   groupXIX 12485182      56  1.277526e-01  44
#> 657    pop    CLF~PAL    groupXI 12833733      57           NaN  34
#> 658    pop    CLF~PAL groupXVIII 12937786      58  2.973689e-01  38
#> 659    pop    CLF~PAL   groupIII 13177027      59  2.321090e-01  44
#> 660    pop    CLF~PAL    groupVI 13335577      60  8.469352e-01  46
#> 661    pop    CLF~PAL     groupX 13718292      61  2.936858e-01  32
#> 662    pop    CLF~PAL    groupVI 13895716      62           NaN  44
#> 663    pop    CLF~PAL   groupXIV 14009990      63  2.452830e-01  30
#> 664    pop    CLF~PAL   groupXII 14158350      64  2.179117e-01  40
#> 665    pop    CLF~PAL    groupIX 14164215      65 -2.897306e-03  46
#> 666    pop    CLF~PAL    groupVI 14236496      66  5.464318e-02  44
#> 667    pop    CLF~PAL     groupX 14273796      67  1.540960e-01  42
#> 668    pop    CLF~PAL  groupVIII 14445389      68  3.363476e-01  38
#> 669    pop    CLF~PAL    groupIX 14868533      69  1.478242e-01  44
#> 670    pop    CLF~PAL    groupIX 15190693      70  7.428012e-01  42
#> 671    pop    CLF~PAL   groupXII 15306362      71  9.591908e-02  44
#> 672    pop    CLF~PAL    groupXX 15460430      72 -5.877624e-02  42
#> 673    pop    CLF~PAL   groupVII 15493025      73 -4.757492e-02  42
#> 674    pop    CLF~PAL    groupXV 15547073      74  1.478242e-01  44
#> 675    pop    CLF~PAL  groupXIII 15697285      75  3.212541e-02  42
#> 676    pop    CLF~PAL   groupIII 15983994      76  8.926518e-01  22
#> 677    pop    CLF~PAL   groupXIX 16233699      77  1.515012e-01  42
#> 678    pop    CLF~PAL    groupIX 16306966      78  3.357287e-01  46
#> 679    pop    CLF~PAL     groupI 16498302      79  1.380217e-01  46
#> 680    pop    CLF~PAL  groupXIII 16520929      80  1.739531e-01  42
#> 681    pop    CLF~PAL   groupXVI 16561891      81  7.424944e-03  42
#> 682    pop    CLF~PAL    groupII 16598364      82           NaN  40
#> 683    pop    CLF~PAL    groupIV 17248473      83  5.601162e-01  44
#> 684    pop    CLF~PAL    groupXX 17748944      84  7.623762e-01  42
#> 685    pop    CLF~PAL  groupXIII 17842228      85  2.089001e-01  38
#> 686    pop    CLF~PAL   groupXIX 17858769      86  2.023825e-01  34
#> 687    pop    CLF~PAL    groupIX 17929225      87  9.471095e-02  44
#> 688    pop    CLF~PAL    groupIX 18339767      88  8.269570e-01  40
#> 689    pop    CLF~PAL  groupXIII 18836558      89 -1.535381e-01  30
#> 690    pop    CLF~PAL    groupII 19244005      90           NaN  46
#> 691    pop    CLF~PAL    groupXX 19317273      91  1.881924e-01  40
#> 692    pop    CLF~PAL     groupI 19759331      92  2.686713e-01  32
#> 693    pop    CLF~PAL   groupVII 20165240      93 -9.992908e-03  40
#> 694    pop    CLF~PAL    groupII 20496069      94  2.240643e-01  36
#> 695    pop    CLF~PAL     groupI 21284266      95  4.540452e-01  46
#> 696    pop    CLF~PAL    groupIV 22227311      96  1.016993e-01  46
#> 697    pop    CLF~PAL     groupI 22746027      97  1.016993e-01  46
#> 698    pop    CLF~PAL   groupVII 25239495      98 -4.683338e-02  42
#> 699    pop    CLF~PAL   groupVII 27599144      99  4.974259e-01  40
#> 700    pop    CLF~PAL    groupIV 32321895     100 -7.175214e-02  30
#> 701    pop    CLF~SMR    groupVI   212436       1 -3.969203e-02  54
#> 702    pop    CLF~SMR groupXVIII   350070       2 -2.378864e-02  56
#> 703    pop    CLF~SMR    groupXV   854285       3 -5.041468e-02  58
#> 704    pop    CLF~SMR  groupXIII   907279       4  5.311598e-02  56
#> 705    pop    CLF~SMR    groupVI  1068620       5  5.528053e-02  36
#> 706    pop    CLF~SMR     groupV  1392753       6  1.073850e-01  54
#> 707    pop    CLF~SMR    groupIX  2075079       7 -5.793631e-02  54
#> 708    pop    CLF~SMR  groupXIII  2106000       8 -4.452982e-02  52
#> 709    pop    CLF~SMR   groupXII  2136587       9 -2.068171e-02  52
#> 710    pop    CLF~SMR  groupXVII  2560839      10 -3.950041e-02  56
#> 711    pop    CLF~SMR   groupIII  2750695      11 -4.829611e-02  50
#> 712    pop    CLF~SMR   groupXIX  3008524      12 -1.529660e-02  46
#> 713    pop    CLF~SMR   groupXVI  3624417      13  6.316841e-02  52
#> 714    pop    CLF~SMR   groupVII  3667528      14 -4.450782e-02  56
#> 715    pop    CLF~SMR   groupVII  4189120      15  1.146369e-01  54
#> 716    pop    CLF~SMR    groupXX  4538074      16 -4.773997e-02  50
#> 717    pop    CLF~SMR  groupXIII  4706089      17  7.384615e-02  54
#> 718    pop    CLF~SMR   groupIII  4746808      18  1.235888e-01  56
#> 719    pop    CLF~SMR  groupXVII  5041995      19  1.271359e-01  58
#> 720    pop    CLF~SMR     groupX  5161637      20 -5.041468e-02  58
#> 721    pop    CLF~SMR   groupXII  5178037      21 -3.552597e-02  48
#> 722    pop    CLF~SMR    groupIX  5287711      22  1.222494e-01  50
#> 723    pop    CLF~SMR  groupVIII  5462135      23  1.897810e-02  52
#> 724    pop    CLF~SMR    groupIX  5960933      24 -1.931944e-02  50
#> 725    pop    CLF~SMR  groupXIII  6076868      25  2.304230e-01  48
#> 726    pop    CLF~SMR    groupXI  6177332      26  1.860632e-01  56
#> 727    pop    CLF~SMR    groupIV  6291287      27 -2.973633e-02  52
#> 728    pop    CLF~SMR    groupIX  6515796      28  1.847354e-01  52
#> 729    pop    CLF~SMR   groupXII  6718079      29 -5.590572e-02  38
#> 730    pop    CLF~SMR   groupVII  6967477      30 -2.199886e-02  54
#> 731    pop    CLF~SMR   groupIII  7194069      31  1.308029e-01  50
#> 732    pop    CLF~SMR   groupXVI  7199722      32 -5.172414e-02  50
#> 733    pop    CLF~SMR  groupVIII  7224065      33 -9.330806e-02  52
#> 734    pop    CLF~SMR groupXVIII  7229426      34  2.274241e-02  56
#> 735    pop    CLF~SMR groupXVIII  8087022      35 -1.142857e-02  50
#> 736    pop    CLF~SMR    groupIX  8373442      36 -4.450782e-02  56
#> 737    pop    CLF~SMR    groupXV  8815538      37 -3.555329e-03  54
#> 738    pop    CLF~SMR    groupXI  9019877      38  2.975113e-01  50
#> 739    pop    CLF~SMR   groupXXI  9067876      39  1.684841e-01  40
#> 740    pop    CLF~SMR  groupXVII  9108184      40 -4.915409e-02  48
#> 741    pop    CLF~SMR    groupIX  9284427      41 -5.231280e-02  44
#> 742    pop    CLF~SMR     groupV  9767769      42  1.030723e-02  54
#> 743    pop    CLF~SMR   groupXII  9881312      43 -1.102864e-02  58
#> 744    pop    CLF~SMR    groupXX 10296516      44  9.934641e-03  54
#> 745    pop    CLF~SMR     groupI 10593829      45           NaN  52
#> 746    pop    CLF~SMR   groupXIX 10765968      46 -9.514655e-02  44
#> 747    pop    CLF~SMR    groupII 10829169      47  1.126899e-01  58
#> 748    pop    CLF~SMR    groupXX 10913115      48 -4.388715e-02  48
#> 749    pop    CLF~SMR groupXVIII 10974594      49  3.068592e-01  60
#> 750    pop    CLF~SMR  groupXVII 11256661      50 -4.468026e-02  48
#> 751    pop    CLF~SMR  groupVIII 11438870      51  1.964007e-01  38
#> 752    pop    CLF~SMR   groupVII 11477223      52  5.059906e-02  60
#> 753    pop    CLF~SMR    groupII 11505084      53  4.970854e-02  54
#> 754    pop    CLF~SMR    groupIX 11904651      54  4.580115e-03  58
#> 755    pop    CLF~SMR  groupVIII 12108945      55 -4.059208e-02  52
#> 756    pop    CLF~SMR   groupXIX 12485182      56 -8.638269e-03  50
#> 757    pop    CLF~SMR    groupXI 12833733      57 -4.989305e-02  46
#> 758    pop    CLF~SMR groupXVIII 12937786      58  5.178571e-02  48
#> 759    pop    CLF~SMR   groupIII 13177027      59  2.744085e-03  52
#> 760    pop    CLF~SMR    groupVI 13335577      60  2.547836e-01  54
#> 761    pop    CLF~SMR     groupX 13718292      61  1.864180e-02  46
#> 762    pop    CLF~SMR    groupVI 13895716      62 -1.018330e-02  60
#> 763    pop    CLF~SMR   groupXIV 14009990      63 -4.957841e-02  42
#> 764    pop    CLF~SMR   groupXII 14158350      64  2.398302e-02  56
#> 765    pop    CLF~SMR    groupIX 14164215      65 -4.849238e-02  56
#> 766    pop    CLF~SMR    groupVI 14236496      66 -6.544590e-02  46
#> 767    pop    CLF~SMR     groupX 14273796      67           NaN  56
#> 768    pop    CLF~SMR  groupVIII 14445389      68 -4.940711e-02  48
#> 769    pop    CLF~SMR    groupIX 14868533      69 -6.117595e-02  52
#> 770    pop    CLF~SMR    groupIX 15190693      70  2.006950e-01  54
#> 771    pop    CLF~SMR   groupXII 15306362      71  3.866746e-02  60
#> 772    pop    CLF~SMR    groupXX 15460430      72  2.191363e-01  56
#> 773    pop    CLF~SMR   groupVII 15493025      73  2.218671e-02  54
#> 774    pop    CLF~SMR    groupXV 15547073      74  6.057727e-03  44
#> 775    pop    CLF~SMR  groupXIII 15697285      75 -5.958204e-02  50
#> 776    pop    CLF~SMR   groupIII 15983994      76 -9.663565e-02  30
#> 777    pop    CLF~SMR   groupXIX 16233699      77 -8.638269e-03  50
#> 778    pop    CLF~SMR    groupIX 16306966      78  1.417195e-02  56
#> 779    pop    CLF~SMR     groupI 16498302      79 -1.849562e-02  56
#> 780    pop    CLF~SMR  groupXIII 16520929      80 -1.330014e-02  58
#> 781    pop    CLF~SMR   groupXVI 16561891      81  1.167964e-01  54
#> 782    pop    CLF~SMR    groupII 16598364      82 -6.953016e-02  44
#> 783    pop    CLF~SMR    groupIV 17248473      83  5.374235e-01  60
#> 784    pop    CLF~SMR    groupXX 17748944      84 -4.771802e-02  50
#> 785    pop    CLF~SMR  groupXIII 17842228      85  1.831090e-01  52
#> 786    pop    CLF~SMR   groupXIX 17858769      86 -8.638269e-03  50
#> 787    pop    CLF~SMR    groupIX 17929225      87  6.127254e-02  56
#> 788    pop    CLF~SMR    groupIX 18339767      88  1.562281e-02  44
#> 789    pop    CLF~SMR  groupXIII 18836558      89 -1.196934e-01  42
#> 790    pop    CLF~SMR    groupII 19244005      90 -4.136898e-02  58
#> 791    pop    CLF~SMR    groupXX 19317273      91 -2.420997e-02  48
#> 792    pop    CLF~SMR     groupI 19759331      92 -6.330289e-02  32
#> 793    pop    CLF~SMR   groupVII 20165240      93 -3.068238e-02  50
#> 794    pop    CLF~SMR    groupII 20496069      94 -5.104232e-02  54
#> 795    pop    CLF~SMR     groupI 21284266      95  3.897612e-02  50
#> 796    pop    CLF~SMR    groupIV 22227311      96 -5.192697e-02  58
#> 797    pop    CLF~SMR     groupI 22746027      97 -4.485815e-02  52
#> 798    pop    CLF~SMR   groupVII 25239495      98 -3.431744e-02  54
#> 799    pop    CLF~SMR   groupVII 27599144      99 -4.803705e-02  52
#> 800    pop    CLF~SMR    groupIV 32321895     100 -3.552597e-02  48
#> 801    pop    CLF~UPD    groupVI   212436       1  1.729034e-02  36
#> 802    pop    CLF~UPD groupXVIII   350070       2           NaN  42
#> 803    pop    CLF~UPD    groupXV   854285       3 -5.786963e-02  42
#> 804    pop    CLF~UPD  groupXIII   907279       4  1.688338e-01  36
#> 805    pop    CLF~UPD    groupVI  1068620       5  2.583519e-01  24
#> 806    pop    CLF~UPD     groupV  1392753       6 -2.158765e-02  40
#> 807    pop    CLF~UPD    groupIX  2075079       7  1.532258e-01  40
#> 808    pop    CLF~UPD  groupXIII  2106000       8  5.376344e-03  32
#> 809    pop    CLF~UPD   groupXII  2136587       9 -3.567208e-02  38
#> 810    pop    CLF~UPD  groupXVII  2560839      10 -3.957404e-02  40
#> 811    pop    CLF~UPD   groupIII  2750695      11 -5.939835e-03  40
#> 812    pop    CLF~UPD   groupXIX  3008524      12  8.624898e-02  32
#> 813    pop    CLF~UPD   groupXVI  3624417      13  1.231294e-01  38
#> 814    pop    CLF~UPD   groupVII  3667528      14  4.638357e-02  42
#> 815    pop    CLF~UPD   groupVII  4189120      15  3.671470e-01  38
#> 816    pop    CLF~UPD    groupXX  4538074      16 -5.691229e-02  34
#> 817    pop    CLF~UPD  groupXIII  4706089      17  6.858579e-02  38
#> 818    pop    CLF~UPD   groupIII  4746808      18 -7.901171e-02  40
#> 819    pop    CLF~UPD  groupXVII  5041995      19  1.039604e-01  42
#> 820    pop    CLF~UPD     groupX  5161637      20 -6.024648e-02  40
#> 821    pop    CLF~UPD   groupXII  5178037      21 -6.867672e-02  32
#> 822    pop    CLF~UPD    groupIX  5287711      22 -6.249991e-02  38
#> 823    pop    CLF~UPD  groupVIII  5462135      23 -3.763441e-02  42
#> 824    pop    CLF~UPD    groupIX  5960933      24  2.556036e-02  40
#> 825    pop    CLF~UPD  groupXIII  6076868      25  1.539197e-01  40
#> 826    pop    CLF~UPD    groupXI  6177332      26  2.050740e-01  36
#> 827    pop    CLF~UPD    groupIV  6291287      27 -6.024648e-02  40
#> 828    pop    CLF~UPD    groupIX  6515796      28 -4.136334e-03  40
#> 829    pop    CLF~UPD   groupXII  6718079      29 -8.161516e-02  34
#> 830    pop    CLF~UPD   groupVII  6967477      30 -5.949852e-02  34
#> 831    pop    CLF~UPD   groupIII  7194069      31  1.878467e-01  42
#> 832    pop    CLF~UPD   groupXVI  7199722      32 -3.843720e-02  36
#> 833    pop    CLF~UPD  groupVIII  7224065      33 -1.152712e-01  34
#> 834    pop    CLF~UPD groupXVIII  7229426      34  2.303626e-01  34
#> 835    pop    CLF~UPD groupXVIII  8087022      35 -1.718170e-02  38
#> 836    pop    CLF~UPD    groupIX  8373442      36 -1.454314e-02  38
#> 837    pop    CLF~UPD    groupXV  8815538      37 -3.843720e-02  36
#> 838    pop    CLF~UPD    groupXI  9019877      38  1.174289e-01  40
#> 839    pop    CLF~UPD   groupXXI  9067876      39  1.961489e-02  26
#> 840    pop    CLF~UPD  groupXVII  9108184      40  1.741718e-02  32
#> 841    pop    CLF~UPD    groupIX  9284427      41  9.304511e-02  30
#> 842    pop    CLF~UPD     groupV  9767769      42 -6.249991e-02  38
#> 843    pop    CLF~UPD   groupXII  9881312      43 -5.321763e-02  36
#> 844    pop    CLF~UPD    groupXX 10296516      44  1.221532e-01  36
#> 845    pop    CLF~UPD     groupI 10593829      45 -9.497965e-03  42
#> 846    pop    CLF~UPD   groupXIX 10765968      46 -3.048902e-02  36
#> 847    pop    CLF~UPD    groupII 10829169      47  5.109489e-02  42
#> 848    pop    CLF~UPD    groupXX 10913115      48 -1.886441e-02  38
#> 849    pop    CLF~UPD groupXVIII 10974594      49 -1.454314e-02  38
#> 850    pop    CLF~UPD  groupXVII 11256661      50 -6.416943e-02  38
#> 851    pop    CLF~UPD  groupVIII 11438870      51  2.578474e-01  38
#> 852    pop    CLF~UPD   groupVII 11477223      52 -7.363647e-02  40
#> 853    pop    CLF~UPD    groupII 11505084      53 -7.692308e-02  32
#> 854    pop    CLF~UPD    groupIX 11904651      54  8.750000e-02  42
#> 855    pop    CLF~UPD  groupVIII 12108945      55  2.920006e-02  38
#> 856    pop    CLF~UPD   groupXIX 12485182      56  2.657068e-02  42
#> 857    pop    CLF~UPD    groupXI 12833733      57           NaN  30
#> 858    pop    CLF~UPD groupXVIII 12937786      58  1.642418e-01  34
#> 859    pop    CLF~UPD   groupIII 13177027      59 -5.342057e-02  38
#> 860    pop    CLF~UPD    groupVI 13335577      60  2.805905e-01  40
#> 861    pop    CLF~UPD     groupX 13718292      61 -6.602641e-02  30
#> 862    pop    CLF~UPD    groupVI 13895716      62  2.556036e-02  40
#> 863    pop    CLF~UPD   groupXIV 14009990      63 -4.125680e-02  36
#> 864    pop    CLF~UPD   groupXII 14158350      64  3.274701e-02  42
#> 865    pop    CLF~UPD    groupIX 14164215      65 -7.680869e-02  40
#> 866    pop    CLF~UPD    groupVI 14236496      66 -7.806955e-02  38
#> 867    pop    CLF~UPD     groupX 14273796      67 -3.538803e-02  38
#> 868    pop    CLF~UPD  groupVIII 14445389      68 -3.663710e-02  40
#> 869    pop    CLF~UPD    groupIX 14868533      69  1.296928e-01  40
#> 870    pop    CLF~UPD    groupIX 15190693      70  1.300773e-01  34
#> 871    pop    CLF~UPD   groupXII 15306362      71 -9.497965e-03  42
#> 872    pop    CLF~UPD    groupXX 15460430      72  3.284691e-05  38
#> 873    pop    CLF~UPD   groupVII 15493025      73 -6.234414e-02  40
#> 874    pop    CLF~UPD    groupXV 15547073      74 -6.773818e-02  36
#> 875    pop    CLF~UPD  groupXIII 15697285      75  4.056162e-02  40
#> 876    pop    CLF~UPD   groupIII 15983994      76  2.424261e-02  34
#> 877    pop    CLF~UPD   groupXIX 16233699      77  1.030534e-02  40
#> 878    pop    CLF~UPD    groupIX 16306966      78  3.086697e-01  42
#> 879    pop    CLF~UPD     groupI 16498302      79  4.259635e-02  36
#> 880    pop    CLF~UPD  groupXIII 16520929      80 -1.596205e-02  38
#> 881    pop    CLF~UPD   groupXVI 16561891      81  4.325237e-02  38
#> 882    pop    CLF~UPD    groupII 16598364      82  5.146713e-02  32
#> 883    pop    CLF~UPD    groupIV 17248473      83  4.118326e-01  40
#> 884    pop    CLF~UPD    groupXX 17748944      84 -2.666018e-02  38
#> 885    pop    CLF~UPD  groupXIII 17842228      85  7.019242e-03  40
#> 886    pop    CLF~UPD   groupXIX 17858769      86  5.713171e-02  38
#> 887    pop    CLF~UPD    groupIX 17929225      87  1.933667e-01  40
#> 888    pop    CLF~UPD    groupIX 18339767      88 -2.188295e-02  36
#> 889    pop    CLF~UPD  groupXIII 18836558      89 -3.976559e-02  30
#> 890    pop    CLF~UPD    groupII 19244005      90 -1.418228e-03  38
#> 891    pop    CLF~UPD    groupXX 19317273      91 -1.072921e-01  36
#> 892    pop    CLF~UPD     groupI 19759331      92 -6.537530e-02  30
#> 893    pop    CLF~UPD   groupVII 20165240      93 -4.440497e-02  40
#> 894    pop    CLF~UPD    groupII 20496069      94  1.377376e-02  38
#> 895    pop    CLF~UPD     groupI 21284266      95  1.210654e-02  36
#> 896    pop    CLF~UPD    groupIV 22227311      96 -6.244608e-02  38
#> 897    pop    CLF~UPD     groupI 22746027      97 -3.206596e-02  42
#> 898    pop    CLF~UPD   groupVII 25239495      98 -7.371027e-02  36
#> 899    pop    CLF~UPD   groupVII 27599144      99 -7.943925e-02  34
#> 900    pop    CLF~UPD    groupIV 32321895     100  3.512244e-02  36
#> 901    pop    OPL~PAL    groupVI   212436       1  2.295653e-01  74
#> 902    pop    OPL~PAL groupXVIII   350070       2  3.360431e-02  80
#> 903    pop    OPL~PAL    groupXV   854285       3 -3.866956e-03  82
#> 904    pop    OPL~PAL  groupXIII   907279       4  3.948354e-01  72
#> 905    pop    OPL~PAL    groupVI  1068620       5 -1.039708e-02  50
#> 906    pop    OPL~PAL     groupV  1392753       6  3.142882e-01  84
#> 907    pop    OPL~PAL    groupIX  2075079       7  1.107191e-01  76
#> 908    pop    OPL~PAL  groupXIII  2106000       8  1.085621e-01  76
#> 909    pop    OPL~PAL   groupXII  2136587       9  1.074538e-02  82
#> 910    pop    OPL~PAL  groupXVII  2560839      10 -8.136920e-03  84
#> 911    pop    OPL~PAL   groupIII  2750695      11  5.529110e-01  74
#> 912    pop    OPL~PAL   groupXIX  3008524      12  1.710363e-02  68
#> 913    pop    OPL~PAL   groupXVI  3624417      13  1.289137e-01  78
#> 914    pop    OPL~PAL   groupVII  3667528      14  1.975531e-02  84
#> 915    pop    OPL~PAL   groupVII  4189120      15  6.125984e-01  72
#> 916    pop    OPL~PAL    groupXX  4538074      16  1.860843e-03  70
#> 917    pop    OPL~PAL  groupXIII  4706089      17 -1.805866e-02  82
#> 918    pop    OPL~PAL   groupIII  4746808      18  4.096351e-01  82
#> 919    pop    OPL~PAL  groupXVII  5041995      19  4.519003e-02  74
#> 920    pop    OPL~PAL     groupX  5161637      20 -1.048699e-02  78
#> 921    pop    OPL~PAL   groupXII  5178037      21  4.694503e-01  70
#> 922    pop    OPL~PAL    groupIX  5287711      22  2.686799e-01  76
#> 923    pop    OPL~PAL  groupVIII  5462135      23  8.083884e-03  82
#> 924    pop    OPL~PAL    groupIX  5960933      24  4.660226e-02  84
#> 925    pop    OPL~PAL  groupXIII  6076868      25  2.165036e-01  78
#> 926    pop    OPL~PAL    groupXI  6177332      26  2.846790e-01  72
#> 927    pop    OPL~PAL    groupIV  6291287      27  1.677159e-01  82
#> 928    pop    OPL~PAL    groupIX  6515796      28  4.695050e-01  68
#> 929    pop    OPL~PAL   groupXII  6718079      29 -4.267805e-03  72
#> 930    pop    OPL~PAL   groupVII  6967477      30  4.083826e-02  76
#> 931    pop    OPL~PAL   groupIII  7194069      31  2.667886e-01  84
#> 932    pop    OPL~PAL   groupXVI  7199722      32  1.689286e-01  76
#> 933    pop    OPL~PAL  groupVIII  7224065      33  1.317330e-02  68
#> 934    pop    OPL~PAL groupXVIII  7229426      34  3.666471e-01  70
#> 935    pop    OPL~PAL groupXVIII  8087022      35  3.462418e-02  76
#> 936    pop    OPL~PAL    groupIX  8373442      36  1.374487e-01  80
#> 937    pop    OPL~PAL    groupXV  8815538      37  1.654300e-01  80
#> 938    pop    OPL~PAL    groupXI  9019877      38  8.710480e-02  78
#> 939    pop    OPL~PAL   groupXXI  9067876      39 -3.028891e-02  52
#> 940    pop    OPL~PAL  groupXVII  9108184      40 -2.429140e-02  70
#> 941    pop    OPL~PAL    groupIX  9284427      41  6.957421e-01  64
#> 942    pop    OPL~PAL     groupV  9767769      42  2.936126e-01  78
#> 943    pop    OPL~PAL   groupXII  9881312      43  1.578409e-01  76
#> 944    pop    OPL~PAL    groupXX 10296516      44  3.074312e-01  82
#> 945    pop    OPL~PAL     groupI 10593829      45           NaN  72
#> 946    pop    OPL~PAL   groupXIX 10765968      46  3.846769e-01  68
#> 947    pop    OPL~PAL    groupII 10829169      47  2.605517e-01  76
#> 948    pop    OPL~PAL    groupXX 10913115      48  3.317210e-02  68
#> 949    pop    OPL~PAL groupXVIII 10974594      49 -9.130023e-03  88
#> 950    pop    OPL~PAL  groupXVII 11256661      50  1.741091e-02  78
#> 951    pop    OPL~PAL  groupVIII 11438870      51  3.573088e-02  62
#> 952    pop    OPL~PAL   groupVII 11477223      52  9.168409e-03  72
#> 953    pop    OPL~PAL    groupII 11505084      53  2.670269e-01  74
#> 954    pop    OPL~PAL    groupIX 11904651      54 -1.000197e-02  84
#> 955    pop    OPL~PAL  groupVIII 12108945      55  9.644132e-03  72
#> 956    pop    OPL~PAL   groupXIX 12485182      56  3.472223e-02  72
#> 957    pop    OPL~PAL    groupXI 12833733      57  5.270332e-02  64
#> 958    pop    OPL~PAL groupXVIII 12937786      58  3.429212e-02  76
#> 959    pop    OPL~PAL   groupIII 13177027      59  2.348957e-01  80
#> 960    pop    OPL~PAL    groupVI 13335577      60  3.048067e-01  80
#> 961    pop    OPL~PAL     groupX 13718292      61 -8.304062e-03  60
#> 962    pop    OPL~PAL    groupVI 13895716      62  1.623117e-01  84
#> 963    pop    OPL~PAL   groupXIV 14009990      63  4.833309e-01  56
#> 964    pop    OPL~PAL   groupXII 14158350      64  1.098335e-01  70
#> 965    pop    OPL~PAL    groupIX 14164215      65  2.155158e-01  72
#> 966    pop    OPL~PAL    groupVI 14236496      66  1.508516e-01  72
#> 967    pop    OPL~PAL     groupX 14273796      67  3.084120e-01  80
#> 968    pop    OPL~PAL  groupVIII 14445389      68  3.890263e-01  74
#> 969    pop    OPL~PAL    groupIX 14868533      69  2.430988e-01  84
#> 970    pop    OPL~PAL    groupIX 15190693      70  2.239049e-01  74
#> 971    pop    OPL~PAL   groupXII 15306362      71  1.561348e-01  86
#> 972    pop    OPL~PAL    groupXX 15460430      72  1.691857e-01  78
#> 973    pop    OPL~PAL   groupVII 15493025      73  6.357983e-02  76
#> 974    pop    OPL~PAL    groupXV 15547073      74  2.931822e-01  72
#> 975    pop    OPL~PAL  groupXIII 15697285      75  1.269403e-02  78
#> 976    pop    OPL~PAL   groupIII 15983994      76  5.575318e-01  48
#> 977    pop    OPL~PAL   groupXIX 16233699      77  7.845122e-02  76
#> 978    pop    OPL~PAL    groupIX 16306966      78  2.506815e-01  78
#> 979    pop    OPL~PAL     groupI 16498302      79  3.261721e-02  86
#> 980    pop    OPL~PAL  groupXIII 16520929      80  4.899504e-02  74
#> 981    pop    OPL~PAL   groupXVI 16561891      81 -1.876845e-02  74
#> 982    pop    OPL~PAL    groupII 16598364      82  5.996148e-02  86
#> 983    pop    OPL~PAL    groupIV 17248473      83           NaN  84
#> 984    pop    OPL~PAL    groupXX 17748944      84  6.085605e-01  70
#> 985    pop    OPL~PAL  groupXIII 17842228      85  4.958126e-01  70
#> 986    pop    OPL~PAL   groupXIX 17858769      86  6.023879e-02  62
#> 987    pop    OPL~PAL    groupIX 17929225      87  1.818702e-01  82
#> 988    pop    OPL~PAL    groupIX 18339767      88  6.272441e-01  68
#> 989    pop    OPL~PAL  groupXIII 18836558      89 -5.146809e-03  60
#> 990    pop    OPL~PAL    groupII 19244005      90           NaN  86
#> 991    pop    OPL~PAL    groupXX 19317273      91  3.003514e-01  74
#> 992    pop    OPL~PAL     groupI 19759331      92  1.217007e-01  50
#> 993    pop    OPL~PAL   groupVII 20165240      93  2.001948e-02  80
#> 994    pop    OPL~PAL    groupII 20496069      94  2.168740e-01  70
#> 995    pop    OPL~PAL     groupI 21284266      95  1.327420e-01  78
#> 996    pop    OPL~PAL    groupIV 22227311      96  1.368291e-01  80
#> 997    pop    OPL~PAL     groupI 22746027      97  1.775555e-02  76
#> 998    pop    OPL~PAL   groupVII 25239495      98  2.312944e-01  84
#> 999    pop    OPL~PAL   groupVII 27599144      99  5.169606e-01  70
#> 1000   pop    OPL~PAL    groupIV 32321895     100 -3.746387e-02  62
#> 1001   pop    OPL~SMR    groupVI   212436       1 -2.134817e-02  84
#> 1002   pop    OPL~SMR groupXVIII   350070       2 -1.805676e-02  90
#> 1003   pop    OPL~SMR    groupXV   854285       3 -1.839844e-02  96
#> 1004   pop    OPL~SMR  groupXIII   907279       4  4.101438e-03  86
#> 1005   pop    OPL~SMR    groupVI  1068620       5 -2.379241e-02  56
#> 1006   pop    OPL~SMR     groupV  1392753       6 -1.931936e-02  94
#> 1007   pop    OPL~SMR    groupIX  2075079       7 -1.158537e-02  84
#> 1008   pop    OPL~SMR  groupXIII  2106000       8  6.414193e-02  88
#> 1009   pop    OPL~SMR   groupXII  2136587       9  6.153640e-02  90
#> 1010   pop    OPL~SMR  groupXVII  2560839      10 -1.080298e-02  94
#> 1011   pop    OPL~SMR   groupIII  2750695      11 -4.602397e-03  80
#> 1012   pop    OPL~SMR   groupXIX  3008524      12 -2.229892e-02  70
#> 1013   pop    OPL~SMR   groupXVI  3624417      13 -1.891660e-02  88
#> 1014   pop    OPL~SMR   groupVII  3667528      14  1.955153e-02  94
#> 1015   pop    OPL~SMR   groupVII  4189120      15 -2.211058e-02  82
#> 1016   pop    OPL~SMR    groupXX  4538074      16 -1.964312e-02  80
#> 1017   pop    OPL~SMR  groupXIII  4706089      17 -1.595170e-02  92
#> 1018   pop    OPL~SMR   groupIII  4746808      18 -1.280859e-02 100
#> 1019   pop    OPL~SMR  groupXVII  5041995      19 -2.191036e-02  86
#> 1020   pop    OPL~SMR     groupX  5161637      20  3.427234e-03  94
#> 1021   pop    OPL~SMR   groupXII  5178037      21  4.669184e-02  78
#> 1022   pop    OPL~SMR    groupIX  5287711      22  4.668986e-02  82
#> 1023   pop    OPL~SMR  groupVIII  5462135      23  7.627025e-03  90
#> 1024   pop    OPL~SMR    groupIX  5960933      24 -1.560056e-02  90
#> 1025   pop    OPL~SMR  groupXIII  6076868      25  9.895353e-02  84
#> 1026   pop    OPL~SMR    groupXI  6177332      26  5.011583e-02  86
#> 1027   pop    OPL~SMR    groupIV  6291287      27 -6.924540e-03  88
#> 1028   pop    OPL~SMR    groupIX  6515796      28  1.184023e-02  78
#> 1029   pop    OPL~SMR   groupXII  6718079      29 -3.136669e-02  72
#> 1030   pop    OPL~SMR   groupVII  6967477      30 -1.102725e-02  86
#> 1031   pop    OPL~SMR   groupIII  7194069      31 -1.614973e-02  88
#> 1032   pop    OPL~SMR   groupXVI  7199722      32 -1.217079e-02  86
#> 1033   pop    OPL~SMR  groupVIII  7224065      33 -2.921057e-02  82
#> 1034   pop    OPL~SMR groupXVIII  7229426      34 -7.790107e-03  80
#> 1035   pop    OPL~SMR groupXVIII  8087022      35 -2.276736e-02  84
#> 1036   pop    OPL~SMR    groupIX  8373442      36 -2.078523e-02  92
#> 1037   pop    OPL~SMR    groupXV  8815538      37  1.222255e-03  92
#> 1038   pop    OPL~SMR    groupXI  9019877      38  9.826265e-02  84
#> 1039   pop    OPL~SMR   groupXXI  9067876      39 -2.915206e-02  64
#> 1040   pop    OPL~SMR  groupXVII  9108184      40 -2.054352e-02  78
#> 1041   pop    OPL~SMR    groupIX  9284427      41 -1.218797e-02  72
#> 1042   pop    OPL~SMR     groupV  9767769      42  5.821202e-03  90
#> 1043   pop    OPL~SMR   groupXII  9881312      43  1.125535e-02  90
#> 1044   pop    OPL~SMR    groupXX 10296516      44 -1.795388e-03  94
#> 1045   pop    OPL~SMR     groupI 10593829      45           NaN  78
#> 1046   pop    OPL~SMR   groupXIX 10765968      46  6.625757e-02  70
#> 1047   pop    OPL~SMR    groupII 10829169      47  9.811852e-02  92
#> 1048   pop    OPL~SMR    groupXX 10913115      48  9.989025e-03  76
#> 1049   pop    OPL~SMR groupXVIII 10974594      49 -2.233281e-03 102
#> 1050   pop    OPL~SMR  groupXVII 11256661      50 -1.828809e-03  84
#> 1051   pop    OPL~SMR  groupVIII 11438870      51 -2.660168e-02  68
#> 1052   pop    OPL~SMR   groupVII 11477223      52  7.272727e-02  96
#> 1053   pop    OPL~SMR    groupII 11505084      53  2.873714e-02  88
#> 1054   pop    OPL~SMR    groupIX 11904651      54 -2.075100e-02  98
#> 1055   pop    OPL~SMR  groupVIII 12108945      55  3.632464e-02  82
#> 1056   pop    OPL~SMR   groupXIX 12485182      56 -1.977367e-02  78
#> 1057   pop    OPL~SMR    groupXI 12833733      57  1.735109e-02  76
#> 1058   pop    OPL~SMR groupXVIII 12937786      58  2.638591e-02  86
#> 1059   pop    OPL~SMR   groupIII 13177027      59 -3.781219e-03  88
#> 1060   pop    OPL~SMR    groupVI 13335577      60 -1.827587e-02  88
#> 1061   pop    OPL~SMR     groupX 13718292      61 -2.188483e-02  74
#> 1062   pop    OPL~SMR    groupVI 13895716      62  7.331152e-02 100
#> 1063   pop    OPL~SMR   groupXIV 14009990      63 -3.183970e-03  68
#> 1064   pop    OPL~SMR   groupXII 14158350      64 -1.285293e-02  86
#> 1065   pop    OPL~SMR    groupIX 14164215      65  1.323447e-01  82
#> 1066   pop    OPL~SMR    groupVI 14236496      66 -2.331006e-02  74
#> 1067   pop    OPL~SMR     groupX 14273796      67           NaN  94
#> 1068   pop    OPL~SMR  groupVIII 14445389      68 -1.698986e-02  84
#> 1069   pop    OPL~SMR    groupIX 14868533      69  9.845411e-02  92
#> 1070   pop    OPL~SMR    groupIX 15190693      70 -1.563068e-02  86
#> 1071   pop    OPL~SMR   groupXII 15306362      71  6.849076e-02 102
#> 1072   pop    OPL~SMR    groupXX 15460430      72  1.975920e-03  92
#> 1073   pop    OPL~SMR   groupVII 15493025      73 -1.983848e-02  88
#> 1074   pop    OPL~SMR    groupXV 15547073      74 -1.519478e-03  72
#> 1075   pop    OPL~SMR  groupXIII 15697285      75 -1.690210e-02  86
#> 1076   pop    OPL~SMR   groupIII 15983994      76  7.182142e-02  56
#> 1077   pop    OPL~SMR   groupXIX 16233699      77 -1.937375e-02  84
#> 1078   pop    OPL~SMR    groupIX 16306966      78  7.581643e-02  88
#> 1079   pop    OPL~SMR     groupI 16498302      79 -2.092312e-02  96
#> 1080   pop    OPL~SMR  groupXIII 16520929      80 -8.118933e-03  90
#> 1081   pop    OPL~SMR   groupXVI 16561891      81 -1.368125e-03  86
#> 1082   pop    OPL~SMR    groupII 16598364      82 -8.260023e-03  90
#> 1083   pop    OPL~SMR    groupIV 17248473      83  1.702692e-03 100
#> 1084   pop    OPL~SMR    groupXX 17748944      84 -6.553020e-03  78
#> 1085   pop    OPL~SMR  groupXIII 17842228      85 -1.218837e-02  84
#> 1086   pop    OPL~SMR   groupXIX 17858769      86 -1.478802e-02  78
#> 1087   pop    OPL~SMR    groupIX 17929225      87 -1.675519e-02  94
#> 1088   pop    OPL~SMR    groupIX 18339767      88 -8.914159e-03  72
#> 1089   pop    OPL~SMR  groupXIII 18836558      89 -2.695259e-02  72
#> 1090   pop    OPL~SMR    groupII 19244005      90  2.725145e-03  98
#> 1091   pop    OPL~SMR    groupXX 19317273      91 -1.727529e-02  82
#> 1092   pop    OPL~SMR     groupI 19759331      92  5.065203e-02  50
#> 1093   pop    OPL~SMR   groupVII 20165240      93  5.961912e-02  90
#> 1094   pop    OPL~SMR    groupII 20496069      94  1.680672e-03  88
#> 1095   pop    OPL~SMR     groupI 21284266      95 -2.725701e-02  82
#> 1096   pop    OPL~SMR    groupIV 22227311      96  1.453055e-02  92
#> 1097   pop    OPL~SMR     groupI 22746027      97 -2.427754e-02  82
#> 1098   pop    OPL~SMR   groupVII 25239495      98  2.933719e-02  96
#> 1099   pop    OPL~SMR   groupVII 27599144      99 -1.774818e-02  82
#> 1100   pop    OPL~SMR    groupIV 32321895     100 -1.897185e-02  80
#> 1101   pop    OPL~UPD    groupVI   212436       1 -1.449303e-02  66
#> 1102   pop    OPL~UPD groupXVIII   350070       2  2.896136e-02  76
#> 1103   pop    OPL~UPD    groupXV   854285       3 -2.450199e-02  80
#> 1104   pop    OPL~UPD  groupXIII   907279       4 -1.307925e-02  66
#> 1105   pop    OPL~UPD    groupVI  1068620       5 -1.001233e-02  44
#> 1106   pop    OPL~UPD     groupV  1392753       6  2.627840e-01  80
#> 1107   pop    OPL~UPD    groupIX  2075079       7  9.894581e-02  70
#> 1108   pop    OPL~UPD  groupXIII  2106000       8 -1.264578e-02  68
#> 1109   pop    OPL~UPD   groupXII  2136587       9  7.526476e-02  76
#> 1110   pop    OPL~UPD  groupXVII  2560839      10  9.553899e-02  78
#> 1111   pop    OPL~UPD   groupIII  2750695      11  5.860181e-02  70
#> 1112   pop    OPL~UPD   groupXIX  3008524      12  3.189927e-03  56
#> 1113   pop    OPL~UPD   groupXVI  3624417      13 -4.579663e-03  74
#> 1114   pop    OPL~UPD   groupVII  3667528      14  2.114714e-01  80
#> 1115   pop    OPL~UPD   groupVII  4189120      15  1.277044e-01  66
#> 1116   pop    OPL~UPD    groupXX  4538074      16  2.119074e-03  64
#> 1117   pop    OPL~UPD  groupXIII  4706089      17 -2.489023e-02  76
#> 1118   pop    OPL~UPD   groupIII  4746808      18  9.386845e-02  84
#> 1119   pop    OPL~UPD  groupXVII  5041995      19 -2.802425e-02  70
#> 1120   pop    OPL~UPD     groupX  5161637      20  5.980812e-03  76
#> 1121   pop    OPL~UPD   groupXII  5178037      21  3.021005e-02  62
#> 1122   pop    OPL~UPD    groupIX  5287711      22 -2.305222e-02  70
#> 1123   pop    OPL~UPD  groupVIII  5462135      23 -2.638047e-02  80
#> 1124   pop    OPL~UPD    groupIX  5960933      24 -3.313811e-02  80
#> 1125   pop    OPL~UPD  groupXIII  6076868      25  4.257685e-02  76
#> 1126   pop    OPL~UPD    groupXI  6177332      26  4.436733e-02  66
#> 1127   pop    OPL~UPD    groupIV  6291287      27  4.650952e-02  76
#> 1128   pop    OPL~UPD    groupIX  6515796      28  2.225418e-01  66
#> 1129   pop    OPL~UPD   groupXII  6718079      29  3.030616e-02  68
#> 1130   pop    OPL~UPD   groupVII  6967477      30  1.823181e-02  66
#> 1131   pop    OPL~UPD   groupIII  7194069      31  1.213942e-03  80
#> 1132   pop    OPL~UPD   groupXVI  7199722      32  4.080943e-02  72
#> 1133   pop    OPL~UPD  groupVIII  7224065      33 -4.223192e-02  64
#> 1134   pop    OPL~UPD groupXVIII  7229426      34  1.143649e-01  58
#> 1135   pop    OPL~UPD groupXVIII  8087022      35 -2.858595e-02  72
#> 1136   pop    OPL~UPD    groupIX  8373442      36  1.673293e-02  74
#> 1137   pop    OPL~UPD    groupXV  8815538      37 -1.976128e-02  74
#> 1138   pop    OPL~UPD    groupXI  9019877      38 -2.265115e-02  74
#> 1139   pop    OPL~UPD   groupXXI  9067876      39 -2.663934e-02  50
#> 1140   pop    OPL~UPD  groupXVII  9108184      40  1.030729e-01  62
#> 1141   pop    OPL~UPD    groupIX  9284427      41  1.859359e-01  58
#> 1142   pop    OPL~UPD     groupV  9767769      42 -7.769257e-03  74
#> 1143   pop    OPL~UPD   groupXII  9881312      43 -1.935996e-02  68
#> 1144   pop    OPL~UPD    groupXX 10296516      44  8.963142e-02  76
#> 1145   pop    OPL~UPD     groupI 10593829      45  4.768635e-02  68
#> 1146   pop    OPL~UPD   groupXIX 10765968      46 -2.615653e-02  62
#> 1147   pop    OPL~UPD    groupII 10829169      47  1.695025e-01  76
#> 1148   pop    OPL~UPD    groupXX 10913115      48  4.453684e-02  66
#> 1149   pop    OPL~UPD groupXVIII 10974594      49  1.923061e-02  80
#> 1150   pop    OPL~UPD  groupXVII 11256661      50 -1.880508e-02  74
#> 1151   pop    OPL~UPD  groupVIII 11438870      51 -2.572183e-02  68
#> 1152   pop    OPL~UPD   groupVII 11477223      52 -2.702379e-02  76
#> 1153   pop    OPL~UPD    groupII 11505084      53 -4.021618e-03  66
#> 1154   pop    OPL~UPD    groupIX 11904651      54 -1.147305e-02  82
#> 1155   pop    OPL~UPD  groupVIII 12108945      55 -1.412387e-03  68
#> 1156   pop    OPL~UPD   groupXIX 12485182      56 -1.907466e-02  70
#> 1157   pop    OPL~UPD    groupXI 12833733      57  4.411765e-02  60
#> 1158   pop    OPL~UPD groupXVIII 12937786      58 -2.685983e-02  72
#> 1159   pop    OPL~UPD   groupIII 13177027      59 -2.514262e-02  74
#> 1160   pop    OPL~UPD    groupVI 13335577      60 -1.694580e-02  74
#> 1161   pop    OPL~UPD     groupX 13718292      61 -5.277616e-03  58
#> 1162   pop    OPL~UPD    groupVI 13895716      62  1.566110e-02  80
#> 1163   pop    OPL~UPD   groupXIV 14009990      63 -2.404408e-02  62
#> 1164   pop    OPL~UPD   groupXII 14158350      64 -1.130772e-02  72
#> 1165   pop    OPL~UPD    groupIX 14164215      65  7.390055e-02  66
#> 1166   pop    OPL~UPD    groupVI 14236496      66 -1.878215e-02  66
#> 1167   pop    OPL~UPD     groupX 14273796      67  2.732432e-02  76
#> 1168   pop    OPL~UPD  groupVIII 14445389      68  1.162612e-02  76
#> 1169   pop    OPL~UPD    groupIX 14868533      69  2.317757e-01  80
#> 1170   pop    OPL~UPD    groupIX 15190693      70 -1.948458e-02  66
#> 1171   pop    OPL~UPD   groupXII 15306362      71  8.153391e-03  84
#> 1172   pop    OPL~UPD    groupXX 15460430      72  2.804125e-01  74
#> 1173   pop    OPL~UPD   groupVII 15493025      73  2.069906e-02  74
#> 1174   pop    OPL~UPD    groupXV 15547073      74  8.670536e-02  64
#> 1175   pop    OPL~UPD  groupXIII 15697285      75  1.144074e-01  76
#> 1176   pop    OPL~UPD   groupIII 15983994      76 -3.297755e-02  60
#> 1177   pop    OPL~UPD   groupXIX 16233699      77 -1.779194e-02  74
#> 1178   pop    OPL~UPD    groupIX 16306966      78  2.386076e-01  74
#> 1179   pop    OPL~UPD     groupI 16498302      79 -2.746462e-02  76
#> 1180   pop    OPL~UPD  groupXIII 16520929      80 -1.705253e-02  70
#> 1181   pop    OPL~UPD   groupXVI 16561891      81 -2.150657e-02  70
#> 1182   pop    OPL~UPD    groupII 16598364      82  1.621622e-02  78
#> 1183   pop    OPL~UPD    groupIV 17248473      83  2.386002e-02  80
#> 1184   pop    OPL~UPD    groupXX 17748944      84 -2.051314e-03  66
#> 1185   pop    OPL~UPD  groupXIII 17842228      85  1.353638e-02  72
#> 1186   pop    OPL~UPD   groupXIX 17858769      86 -1.744033e-02  66
#> 1187   pop    OPL~UPD    groupIX 17929225      87  5.300169e-02  78
#> 1188   pop    OPL~UPD    groupIX 18339767      88 -2.815583e-02  64
#> 1189   pop    OPL~UPD  groupXIII 18836558      89  5.854168e-02  60
#> 1190   pop    OPL~UPD    groupII 19244005      90  8.657465e-02  78
#> 1191   pop    OPL~UPD    groupXX 19317273      91  1.391902e-02  70
#> 1192   pop    OPL~UPD     groupI 19759331      92  3.395840e-02  48
#> 1193   pop    OPL~UPD   groupVII 20165240      93  1.364050e-01  80
#> 1194   pop    OPL~UPD    groupII 20496069      94 -1.783053e-02  72
#> 1195   pop    OPL~UPD     groupI 21284266      95 -4.003460e-02  68
#> 1196   pop    OPL~UPD    groupIV 22227311      96  1.349458e-02  72
#> 1197   pop    OPL~UPD     groupI 22746027      97 -2.629273e-02  72
#> 1198   pop    OPL~UPD   groupVII 25239495      98  1.693543e-01  78
#> 1199   pop    OPL~UPD   groupVII 27599144      99 -1.893491e-02  64
#> 1200   pop    OPL~UPD    groupIV 32321895     100  1.487928e-01  68
#> 1201   pop    PAL~SMR    groupVI   212436       1  2.414190e-01  78
#> 1202   pop    PAL~SMR groupXVIII   350070       2  1.539447e-02  78
#> 1203   pop    PAL~SMR    groupXV   854285       3 -1.416444e-02  78
#> 1204   pop    PAL~SMR  groupXIII   907279       4  2.691158e-01  78
#> 1205   pop    PAL~SMR    groupVI  1068620       5 -3.571065e-02  54
#> 1206   pop    PAL~SMR     groupV  1392753       6  2.692345e-01  74
#> 1207   pop    PAL~SMR    groupIX  2075079       7  1.578409e-01  76
#> 1208   pop    PAL~SMR  groupXIII  2106000       8  2.942282e-01  68
#> 1209   pop    PAL~SMR   groupXII  2136587       9  1.517133e-01  72
#> 1210   pop    PAL~SMR  groupXVII  2560839      10  1.539447e-02  78
#> 1211   pop    PAL~SMR   groupIII  2750695      11  4.661160e-01  70
#> 1212   pop    PAL~SMR   groupXIX  3008524      12  2.930820e-02  66
#> 1213   pop    PAL~SMR   groupXVI  3624417      13  8.795249e-02  70
#> 1214   pop    PAL~SMR   groupVII  3667528      14  1.467447e-01  78
#> 1215   pop    PAL~SMR   groupVII  4189120      15  5.956436e-01  74
#> 1216   pop    PAL~SMR    groupXX  4538074      16  1.860843e-03  70
#> 1217   pop    PAL~SMR  groupXIII  4706089      17 -2.678421e-02  74
#> 1218   pop    PAL~SMR   groupIII  4746808      18  4.853537e-01  74
#> 1219   pop    PAL~SMR  groupXVII  5041995      19  3.360431e-02  80
#> 1220   pop    PAL~SMR     groupX  5161637      20  2.896136e-02  76
#> 1221   pop    PAL~SMR   groupXII  5178037      21  2.915914e-01  72
#> 1222   pop    PAL~SMR    groupIX  5287711      22  1.199833e-01  74
#> 1223   pop    PAL~SMR  groupVIII  5462135      23  6.686478e-02  72
#> 1224   pop    PAL~SMR    groupIX  5960933      24  2.125513e-02  70
#> 1225   pop    PAL~SMR  groupXIII  6076868      25  4.661160e-01  70
#> 1226   pop    PAL~SMR    groupXI  6177332      26  9.251653e-02  74
#> 1227   pop    PAL~SMR    groupIV  6291287      27  1.199833e-01  74
#> 1228   pop    PAL~SMR    groupIX  6515796      28  3.144214e-01  70
#> 1229   pop    PAL~SMR   groupXII  6718079      29 -1.970443e-02  60
#> 1230   pop    PAL~SMR   groupVII  6967477      30  8.055799e-02  78
#> 1231   pop    PAL~SMR   groupIII  7194069      31  3.394360e-01  72
#> 1232   pop    PAL~SMR   groupXVI  7199722      32  1.056051e-01  70
#> 1233   pop    PAL~SMR  groupVIII  7224065      33 -3.889503e-03  74
#> 1234   pop    PAL~SMR groupXVIII  7229426      34  2.463234e-01  78
#> 1235   pop    PAL~SMR groupXVIII  8087022      35  4.251969e-02  72
#> 1236   pop    PAL~SMR    groupIX  8373442      36  1.424934e-01  76
#> 1237   pop    PAL~SMR    groupXV  8815538      37  2.526803e-01  76
#> 1238   pop    PAL~SMR    groupXI  9019877      38 -1.925145e-02  70
#> 1239   pop    PAL~SMR   groupXXI  9067876      39 -2.279507e-02  56
#> 1240   pop    PAL~SMR  groupXVII  9108184      40 -2.528596e-02  68
#> 1241   pop    PAL~SMR    groupIX  9284427      41  6.385074e-01  60
#> 1242   pop    PAL~SMR     groupV  9767769      42  4.235395e-01  76
#> 1243   pop    PAL~SMR   groupXII  9881312      43  6.915983e-02  82
#> 1244   pop    PAL~SMR    groupXX 10296516      44  1.712556e-01  72
#> 1245   pop    PAL~SMR     groupI 10593829      45           NaN  74
#> 1246   pop    PAL~SMR   groupXIX 10765968      46  1.633790e-01  66
#> 1247   pop    PAL~SMR    groupII 10829169      47  3.246450e-02  76
#> 1248   pop    PAL~SMR    groupXX 10913115      48 -2.948129e-02  64
#> 1249   pop    PAL~SMR groupXVIII 10974594      49           NaN  82
#> 1250   pop    PAL~SMR  groupXVII 11256661      50 -2.335848e-02  70
#> 1251   pop    PAL~SMR  groupVIII 11438870      51  1.469926e-02  50
#> 1252   pop    PAL~SMR   groupVII 11477223      52  2.068441e-01  72
#> 1253   pop    PAL~SMR    groupII 11505084      53  1.256814e-01  78
#> 1254   pop    PAL~SMR    groupIX 11904651      54 -8.316165e-03  78
#> 1255   pop    PAL~SMR  groupVIII 12108945      55  1.617460e-01  70
#> 1256   pop    PAL~SMR   groupXIX 12485182      56  4.710754e-02  70
#> 1257   pop    PAL~SMR    groupXI 12833733      57 -1.206022e-02  60
#> 1258   pop    PAL~SMR groupXVIII 12937786      58  1.600208e-01  62
#> 1259   pop    PAL~SMR   groupIII 13177027      59  3.264191e-01  72
#> 1260   pop    PAL~SMR    groupVI 13335577      60  2.935753e-01  76
#> 1261   pop    PAL~SMR     groupX 13718292      61  1.096078e-02  66
#> 1262   pop    PAL~SMR    groupVI 13895716      62  2.842091e-02  80
#> 1263   pop    PAL~SMR   groupXIV 14009990      63  3.585089e-01  52
#> 1264   pop    PAL~SMR   groupXII 14158350      64  5.359708e-02  72
#> 1265   pop    PAL~SMR    groupIX 14164215      65 -1.440475e-02  78
#> 1266   pop    PAL~SMR    groupVI 14236496      66  1.096306e-01  66
#> 1267   pop    PAL~SMR     groupX 14273796      67  2.898385e-01  74
#> 1268   pop    PAL~SMR  groupVIII 14445389      68  3.182718e-01  62
#> 1269   pop    PAL~SMR    groupIX 14868533      69  6.402594e-02  76
#> 1270   pop    PAL~SMR    groupIX 15190693      70  2.088595e-01  76
#> 1271   pop    PAL~SMR   groupXII 15306362      71 -6.804900e-03  80
#> 1272   pop    PAL~SMR    groupXX 15460430      72  2.759021e-01  74
#> 1273   pop    PAL~SMR   groupVII 15493025      73  7.529012e-02  72
#> 1274   pop    PAL~SMR    groupXV 15547073      74  2.187500e-01  68
#> 1275   pop    PAL~SMR  groupXIII 15697285      75  3.216610e-02  72
#> 1276   pop    PAL~SMR   groupIII 15983994      76  8.466454e-01  40
#> 1277   pop    PAL~SMR   groupXIX 16233699      77  6.510422e-02  68
#> 1278   pop    PAL~SMR    groupIX 16306966      78  8.055799e-02  78
#> 1279   pop    PAL~SMR     groupI 16498302      79  5.223701e-02  78
#> 1280   pop    PAL~SMR  groupXIII 16520929      80  1.019909e-01  76
#> 1281   pop    PAL~SMR   groupXVI 16561891      81  9.758486e-03  76
#> 1282   pop    PAL~SMR    groupII 16598364      82  2.353133e-02  72
#> 1283   pop    PAL~SMR    groupIV 17248473      83 -8.949470e-03  80
#> 1284   pop    PAL~SMR    groupXX 17748944      84  7.065301e-01  68
#> 1285   pop    PAL~SMR  groupXIII 17842228      85  5.514233e-01  70
#> 1286   pop    PAL~SMR   groupXIX 17858769      86  1.024989e-01  60
#> 1287   pop    PAL~SMR    groupIX 17929225      87  2.312635e-01  76
#> 1288   pop    PAL~SMR    groupIX 18339767      88  5.694140e-01  60
#> 1289   pop    PAL~SMR  groupXIII 18836558      89 -5.146809e-03  60
#> 1290   pop    PAL~SMR    groupII 19244005      90 -6.870895e-03  80
#> 1291   pop    PAL~SMR    groupXX 19317273      91  3.545566e-01  72
#> 1292   pop    PAL~SMR     groupI 19759331      92  3.260563e-01  44
#> 1293   pop    PAL~SMR   groupVII 20165240      93 -2.000689e-02  66
#> 1294   pop    PAL~SMR    groupII 20496069      94  3.306548e-01  70
#> 1295   pop    PAL~SMR     groupI 21284266      95  1.304909e-01  72
#> 1296   pop    PAL~SMR    groupIV 22227311      96  5.401452e-02  80
#> 1297   pop    PAL~SMR     groupI 22746027      97  2.044570e-02  74
#> 1298   pop    PAL~SMR   groupVII 25239495      98  7.377922e-02  76
#> 1299   pop    PAL~SMR   groupVII 27599144      99  4.664053e-01  68
#> 1300   pop    PAL~SMR    groupIV 32321895     100 -2.522607e-02  62
#> 1301   pop    PAL~UPD    groupVI   212436       1  1.471069e-01  60
#> 1302   pop    PAL~UPD groupXVIII   350070       2           NaN  64
#> 1303   pop    PAL~UPD    groupXV   854285       3 -1.854742e-02  62
#> 1304   pop    PAL~UPD  groupXIII   907279       4  4.380040e-01  58
#> 1305   pop    PAL~UPD    groupVI  1068620       5  7.855150e-02  42
#> 1306   pop    PAL~UPD     groupV  1392753       6 -1.836297e-02  60
#> 1307   pop    PAL~UPD    groupIX  2075079       7           NaN  62
#> 1308   pop    PAL~UPD  groupXIII  2106000       8  2.636139e-01  48
#> 1309   pop    PAL~UPD   groupXII  2136587       9  1.827384e-01  58
#> 1310   pop    PAL~UPD  groupXVII  2560839      10  1.313204e-01  62
#> 1311   pop    PAL~UPD   groupIII  2750695      11  4.019652e-01  60
#> 1312   pop    PAL~UPD   groupXIX  3008524      12 -2.409490e-02  52
#> 1313   pop    PAL~UPD   groupXVI  3624417      13  5.539322e-02  56
#> 1314   pop    PAL~UPD   groupVII  3667528      14  3.671497e-01  64
#> 1315   pop    PAL~UPD   groupVII  4189120      15  3.174836e-01  58
#> 1316   pop    PAL~UPD    groupXX  4538074      16  7.374345e-02  54
#> 1317   pop    PAL~UPD  groupXIII  4706089      17 -3.451487e-02  58
#> 1318   pop    PAL~UPD   groupIII  4746808      18  1.123580e-01  58
#> 1319   pop    PAL~UPD  groupXVII  5041995      19  4.192734e-02  64
#> 1320   pop    PAL~UPD     groupX  5161637      20  4.208999e-02  58
#> 1321   pop    PAL~UPD   groupXII  5178037      21  3.468768e-01  56
#> 1322   pop    PAL~UPD    groupIX  5287711      22  3.646070e-01  62
#> 1323   pop    PAL~UPD  groupVIII  5462135      23  2.224614e-03  62
#> 1324   pop    PAL~UPD    groupIX  5960933      24  8.596231e-02  60
#> 1325   pop    PAL~UPD  groupXIII  6076868      25  4.211690e-01  62
#> 1326   pop    PAL~UPD    groupXI  6177332      26  8.415762e-02  54
#> 1327   pop    PAL~UPD    groupIV  6291287      27  4.888510e-02  62
#> 1328   pop    PAL~UPD    groupIX  6515796      28  8.328603e-02  58
#> 1329   pop    PAL~UPD   groupXII  6718079      29 -2.095439e-02  56
#> 1330   pop    PAL~UPD   groupVII  6967477      30  1.701614e-01  58
#> 1331   pop    PAL~UPD   groupIII  7194069      31  4.299287e-01  64
#> 1332   pop    PAL~UPD   groupXVI  7199722      32  2.041242e-02  56
#> 1333   pop    PAL~UPD  groupVIII  7224065      33 -1.093791e-02  56
#> 1334   pop    PAL~UPD groupXVIII  7229426      34  1.392444e-01  56
#> 1335   pop    PAL~UPD groupXVIII  8087022      35  4.555551e-02  60
#> 1336   pop    PAL~UPD    groupIX  8373442      36  5.366268e-02  58
#> 1337   pop    PAL~UPD    groupXV  8815538      37  1.924122e-01  58
#> 1338   pop    PAL~UPD    groupXI  9019877      38  2.220405e-02  60
#> 1339   pop    PAL~UPD   groupXXI  9067876      39 -4.314795e-02  42
#> 1340   pop    PAL~UPD  groupXVII  9108184      40  1.040741e-01  52
#> 1341   pop    PAL~UPD    groupIX  9284427      41  4.021904e-01  46
#> 1342   pop    PAL~UPD     groupV  9767769      42  2.170551e-01  60
#> 1343   pop    PAL~UPD   groupXII  9881312      43  2.383154e-01  60
#> 1344   pop    PAL~UPD    groupXX 10296516      44  5.059878e-02  54
#> 1345   pop    PAL~UPD     groupI 10593829      45  4.192734e-02  64
#> 1346   pop    PAL~UPD   groupXIX 10765968      46  3.442890e-01  58
#> 1347   pop    PAL~UPD    groupII 10829169      47 -2.040816e-02  60
#> 1348   pop    PAL~UPD    groupXX 10913115      48 -3.938544e-02  54
#> 1349   pop    PAL~UPD groupXVIII 10974594      49  5.734340e-02  60
#> 1350   pop    PAL~UPD  groupXVII 11256661      50 -1.624966e-02  60
#> 1351   pop    PAL~UPD  groupVIII 11438870      51  6.516920e-02  50
#> 1352   pop    PAL~UPD   groupVII 11477223      52 -1.133840e-02  52
#> 1353   pop    PAL~UPD    groupII 11505084      53  4.339255e-01  56
#> 1354   pop    PAL~UPD    groupIX 11904651      54           NaN  62
#> 1355   pop    PAL~UPD  groupVIII 12108945      55 -3.285927e-02  56
#> 1356   pop    PAL~UPD   groupXIX 12485182      56  1.006081e-02  62
#> 1357   pop    PAL~UPD    groupXI 12833733      57           NaN  44
#> 1358   pop    PAL~UPD groupXVIII 12937786      58 -2.307332e-03  48
#> 1359   pop    PAL~UPD   groupIII 13177027      59  2.453490e-01  58
#> 1360   pop    PAL~UPD    groupVI 13335577      60  3.320393e-01  62
#> 1361   pop    PAL~UPD     groupX 13718292      61  5.034373e-02  50
#> 1362   pop    PAL~UPD    groupVI 13895716      62  8.596231e-02  60
#> 1363   pop    PAL~UPD   groupXIV 14009990      63  4.131442e-01  46
#> 1364   pop    PAL~UPD   groupXII 14158350      64  4.966423e-02  58
#> 1365   pop    PAL~UPD    groupIX 14164215      65  8.409845e-05  62
#> 1366   pop    PAL~UPD    groupVI 14236496      66  6.674981e-02  58
#> 1367   pop    PAL~UPD     groupX 14273796      67  1.449408e-01  56
#> 1368   pop    PAL~UPD  groupVIII 14445389      68  2.145274e-01  54
#> 1369   pop    PAL~UPD    groupIX 14868533      69           NaN  64
#> 1370   pop    PAL~UPD    groupIX 15190693      70  3.090006e-01  56
#> 1371   pop    PAL~UPD   groupXII 15306362      71  3.660561e-02  62
#> 1372   pop    PAL~UPD    groupXX 15460430      72 -1.029548e-02  56
#> 1373   pop    PAL~UPD   groupVII 15493025      73 -2.917744e-02  58
#> 1374   pop    PAL~UPD    groupXV 15547073      74  1.033825e-01  60
#> 1375   pop    PAL~UPD  groupXIII 15697285      75  2.468444e-01  62
#> 1376   pop    PAL~UPD   groupIII 15983994      76  4.830157e-01  44
#> 1377   pop    PAL~UPD   groupXIX 16233699      77  3.541123e-02  58
#> 1378   pop    PAL~UPD    groupIX 16306966      78           NaN  64
#> 1379   pop    PAL~UPD     groupI 16498302      79  1.157825e-02  58
#> 1380   pop    PAL~UPD  groupXIII 16520929      80  8.589752e-02  56
#> 1381   pop    PAL~UPD   groupXVI 16561891      81 -2.072245e-02  60
#> 1382   pop    PAL~UPD    groupII 16598364      82  1.938593e-01  60
#> 1383   pop    PAL~UPD    groupIV 17248473      83  4.942115e-03  60
#> 1384   pop    PAL~UPD    groupXX 17748944      84  5.438008e-01  56
#> 1385   pop    PAL~UPD  groupXIII 17842228      85  3.601860e-01  58
#> 1386   pop    PAL~UPD   groupXIX 17858769      86  2.262814e-02  48
#> 1387   pop    PAL~UPD    groupIX 17929225      87  4.290008e-01  60
#> 1388   pop    PAL~UPD    groupIX 18339767      88  6.776190e-01  52
#> 1389   pop    PAL~UPD  groupXIII 18836558      89 -5.303030e-02  48
#> 1390   pop    PAL~UPD    groupII 19244005      90  5.734340e-02  60
#> 1391   pop    PAL~UPD    groupXX 19317273      91  1.556803e-01  60
#> 1392   pop    PAL~UPD     groupI 19759331      92  2.971388e-01  42
#> 1393   pop    PAL~UPD   groupVII 20165240      93  8.547009e-03  56
#> 1394   pop    PAL~UPD    groupII 20496069      94  1.243876e-01  54
#> 1395   pop    PAL~UPD     groupI 21284266      95  1.701614e-01  58
#> 1396   pop    PAL~UPD    groupIV 22227311      96  5.734340e-02  60
#> 1397   pop    PAL~UPD     groupI 22746027      97  4.313381e-03  64
#> 1398   pop    PAL~UPD   groupVII 25239495      98 -2.297344e-02  58
#> 1399   pop    PAL~UPD   groupVII 27599144      99  5.095541e-01  50
#> 1400   pop    PAL~UPD    groupIV 32321895     100  1.187208e-01  50
#> 1401   pop    SMR~UPD    groupVI   212436       1 -1.559760e-02  70
#> 1402   pop    SMR~UPD groupXVIII   350070       2  1.137818e-02  74
#> 1403   pop    SMR~UPD    groupXV   854285       3 -2.703369e-02  76
#> 1404   pop    SMR~UPD  groupXIII   907279       4 -8.147852e-03  72
#> 1405   pop    SMR~UPD    groupVI  1068620       5  5.112321e-02  48
#> 1406   pop    SMR~UPD     groupV  1392753       6  2.174303e-01  70
#> 1407   pop    SMR~UPD    groupIX  2075079       7  1.438692e-01  70
#> 1408   pop    SMR~UPD  groupXIII  2106000       8 -9.915159e-03  60
#> 1409   pop    SMR~UPD   groupXII  2136587       9 -2.643267e-02  66
#> 1410   pop    SMR~UPD  groupXVII  2560839      10  3.471202e-02  72
#> 1411   pop    SMR~UPD   groupIII  2750695      11  2.228159e-03  66
#> 1412   pop    SMR~UPD   groupXIX  3008524      12  1.495582e-02  54
#> 1413   pop    SMR~UPD   groupXVI  3624417      13 -2.570359e-02  66
#> 1414   pop    SMR~UPD   groupVII  3667528      14  6.550385e-02  74
#> 1415   pop    SMR~UPD   groupVII  4189120      15  1.161601e-01  68
#> 1416   pop    SMR~UPD    groupXX  4538074      16  2.119074e-03  64
#> 1417   pop    SMR~UPD  groupXIII  4706089      17 -2.610254e-02  68
#> 1418   pop    SMR~UPD   groupIII  4746808      18  1.548358e-01  76
#> 1419   pop    SMR~UPD  groupXVII  5041995      19 -2.703369e-02  76
#> 1420   pop    SMR~UPD     groupX  5161637      20 -2.789467e-02  74
#> 1421   pop    SMR~UPD   groupXII  5178037      21 -2.310005e-02  64
#> 1422   pop    SMR~UPD    groupIX  5287711      22  8.989927e-02  68
#> 1423   pop    SMR~UPD  groupVIII  5462135      23  5.867770e-03  70
#> 1424   pop    SMR~UPD    groupIX  5960933      24 -8.321345e-03  66
#> 1425   pop    SMR~UPD  groupXIII  6076868      25 -2.042553e-02  68
#> 1426   pop    SMR~UPD    groupXI  6177332      26 -1.995385e-02  68
#> 1427   pop    SMR~UPD    groupIV  6291287      27  2.672193e-03  68
#> 1428   pop    SMR~UPD    groupIX  6515796      28  8.780811e-02  68
#> 1429   pop    SMR~UPD   groupXII  6718079      29  9.356884e-03  56
#> 1430   pop    SMR~UPD   groupVII  6967477      30 -1.629569e-02  68
#> 1431   pop    SMR~UPD   groupIII  7194069      31 -2.112524e-02  68
#> 1432   pop    SMR~UPD   groupXVI  7199722      32  2.703582e-03  66
#> 1433   pop    SMR~UPD  groupVIII  7224065      33 -4.438867e-02  70
#> 1434   pop    SMR~UPD groupXVIII  7229426      34  2.755675e-02  66
#> 1435   pop    SMR~UPD groupXVIII  8087022      35 -2.962400e-02  68
#> 1436   pop    SMR~UPD    groupIX  8373442      36  1.883080e-02  70
#> 1437   pop    SMR~UPD    groupXV  8815538      37 -8.291738e-03  70
#> 1438   pop    SMR~UPD    groupXI  9019877      38  3.248330e-02  66
#> 1439   pop    SMR~UPD   groupXXI  9067876      39 -1.730116e-02  54
#> 1440   pop    SMR~UPD  groupXVII  9108184      40  1.240839e-01  60
#> 1441   pop    SMR~UPD    groupIX  9284427      41  1.101035e-01  54
#> 1442   pop    SMR~UPD     groupV  9767769      42  6.640920e-02  72
#> 1443   pop    SMR~UPD   groupXII  9881312      43  3.803663e-02  74
#> 1444   pop    SMR~UPD    groupXX 10296516      44 -4.769632e-03  66
#> 1445   pop    SMR~UPD     groupI 10593829      45  5.043522e-02  70
#> 1446   pop    SMR~UPD   groupXIX 10765968      46  4.073942e-03  60
#> 1447   pop    SMR~UPD    groupII 10829169      47 -1.875999e-02  76
#> 1448   pop    SMR~UPD    groupXX 10913115      48 -2.092817e-02  62
#> 1449   pop    SMR~UPD groupXVIII 10974594      49  8.051604e-02  74
#> 1450   pop    SMR~UPD  groupXVII 11256661      50 -2.765849e-02  66
#> 1451   pop    SMR~UPD  groupVIII 11438870      51 -2.337442e-02  56
#> 1452   pop    SMR~UPD   groupVII 11477223      52  9.125126e-02  76
#> 1453   pop    SMR~UPD    groupII 11505084      53  1.126632e-01  70
#> 1454   pop    SMR~UPD    groupIX 11904651      54 -9.893466e-03  76
#> 1455   pop    SMR~UPD  groupVIII 12108945      55  1.510401e-01  66
#> 1456   pop    SMR~UPD   groupXIX 12485182      56 -1.580878e-02  68
#> 1457   pop    SMR~UPD    groupXI 12833733      57 -1.811907e-02  56
#> 1458   pop    SMR~UPD groupXVIII 12937786      58  2.088791e-02  58
#> 1459   pop    SMR~UPD   groupIII 13177027      59 -1.088195e-02  66
#> 1460   pop    SMR~UPD    groupVI 13335577      60 -2.398310e-02  70
#> 1461   pop    SMR~UPD     groupX 13718292      61 -2.331810e-02  64
#> 1462   pop    SMR~UPD    groupVI 13895716      62 -1.312992e-02  76
#> 1463   pop    SMR~UPD   groupXIV 14009990      63 -2.963386e-02  58
#> 1464   pop    SMR~UPD   groupXII 14158350      64 -2.643735e-02  74
#> 1465   pop    SMR~UPD    groupIX 14164215      65 -2.974365e-02  72
#> 1466   pop    SMR~UPD    groupVI 14236496      66 -3.880328e-02  60
#> 1467   pop    SMR~UPD     groupX 14273796      67  2.157804e-02  70
#> 1468   pop    SMR~UPD  groupVIII 14445389      68 -1.797445e-02  64
#> 1469   pop    SMR~UPD    groupIX 14868533      69  5.797101e-02  72
#> 1470   pop    SMR~UPD    groupIX 15190693      70 -1.728534e-02  68
#> 1471   pop    SMR~UPD   groupXII 15306362      71 -6.248221e-03  78
#> 1472   pop    SMR~UPD    groupXX 15460430      72  3.959555e-01  70
#> 1473   pop    SMR~UPD   groupVII 15493025      73  2.990638e-02  70
#> 1474   pop    SMR~UPD    groupXV 15547073      74  2.102938e-02  60
#> 1475   pop    SMR~UPD  groupXIII 15697285      75  7.670340e-02  70
#> 1476   pop    SMR~UPD   groupIII 15983994      76  1.665917e-01  52
#> 1477   pop    SMR~UPD   groupXIX 16233699      77 -2.194665e-02  66
#> 1478   pop    SMR~UPD    groupIX 16306966      78  7.405996e-02  74
#> 1479   pop    SMR~UPD     groupI 16498302      79 -3.469437e-02  68
#> 1480   pop    SMR~UPD  groupXIII 16520929      80 -2.029031e-02  72
#> 1481   pop    SMR~UPD   groupXVI 16561891      81 -9.352440e-03  72
#> 1482   pop    SMR~UPD    groupII 16598364      82  6.962853e-02  64
#> 1483   pop    SMR~UPD    groupIV 17248473      83 -2.458925e-02  76
#> 1484   pop    SMR~UPD    groupXX 17748944      84  5.681729e-02  64
#> 1485   pop    SMR~UPD  groupXIII 17842228      85  4.025618e-02  72
#> 1486   pop    SMR~UPD   groupXIX 17858769      86 -2.990996e-03  64
#> 1487   pop    SMR~UPD    groupIX 17929225      87  3.216132e-02  72
#> 1488   pop    SMR~UPD    groupIX 18339767      88 -9.856423e-03  56
#> 1489   pop    SMR~UPD  groupXIII 18836558      89  5.854168e-02  60
#> 1490   pop    SMR~UPD    groupII 19244005      90  8.610997e-03  72
#> 1491   pop    SMR~UPD    groupXX 19317273      91  4.138068e-02  68
#> 1492   pop    SMR~UPD     groupI 19759331      92 -3.919331e-02  42
#> 1493   pop    SMR~UPD   groupVII 20165240      93 -1.069526e-02  66
#> 1494   pop    SMR~UPD    groupII 20496069      94  3.943194e-02  72
#> 1495   pop    SMR~UPD     groupI 21284266      95 -2.909282e-02  62
#> 1496   pop    SMR~UPD    groupIV 22227311      96 -2.834855e-02  72
#> 1497   pop    SMR~UPD     groupI 22746027      97 -2.603145e-02  70
#> 1498   pop    SMR~UPD   groupVII 25239495      98  2.729033e-02  70
#> 1499   pop    SMR~UPD   groupVII 27599144      99 -4.008314e-02  62
#> 1500   pop    SMR~UPD    groupIV 32321895     100  1.297443e-01  68
#> 
#> $weighted.means
#>    facet subfacet snp.facet snp.subfacet weighted_mean_fst     mean_fst
#> 1    pop  ASP~CLF     .base        .base       0.030144775  0.030079042
#> 2    pop  ASP~OPL     .base        .base       0.044866239  0.044622295
#> 3    pop  ASP~PAL     .base        .base       0.166795773  0.166960267
#> 4    pop  ASP~SMR     .base        .base       0.034677072  0.036034847
#> 5    pop  ASP~UPD     .base        .base      -0.006846049 -0.006969351
#> 6    pop  CLF~OPL     .base        .base       0.037613548  0.035375265
#> 7    pop  CLF~PAL     .base        .base       0.275578922  0.275888985
#> 8    pop  CLF~SMR     .base        .base       0.037197014  0.034231385
#> 9    pop  CLF~UPD     .base        .base       0.033804249  0.033489757
#> 10   pop  OPL~PAL     .base        .base       0.234274249  0.237995142
#> 11   pop  OPL~SMR     .base        .base       0.009675088  0.009257510
#> 12   pop  OPL~UPD     .base        .base       0.040626811  0.039056581
#> 13   pop  PAL~SMR     .base        .base       0.224948683  0.233910701
#> 14   pop  PAL~UPD     .base        .base       0.167063505  0.168489328
#> 15   pop  SMR~UPD     .base        .base       0.030598899  0.030039013
#> 
#> $fst.matrix
#> $fst.matrix$pop
#> Key: <p1>
#>        p1        CLF        OPL       PAL         SMR          UPD
#>    <fctr>      <num>      <num>     <num>       <num>        <num>
#> 1:    ASP 0.03014477 0.04486624 0.1667958 0.034677072 -0.006846049
#> 2:    CLF         NA 0.03761355 0.2755789 0.037197014  0.033804249
#> 3:    OPL         NA         NA 0.2342742 0.009675088  0.040626811
#> 4:    PAL         NA         NA        NA 0.224948683  0.167063505
#> 5:    SMR         NA         NA        NA          NA  0.030598899
#> 
#> 

if (FALSE) { # \dontrun{
# Using genepop
x <- calc_pairwise_fst(stickSNPs, "pop", "genepop")
get.snpR.stats(x, "pop", "fst")

# bootstrap p-values for overall pairwise-Fst values
x <- calc_pairwise_fst(stickSNPs, "pop", boot = 5)
get.snpR.stats(x, "pop", "fst")
} # }
```
