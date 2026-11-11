# Calculate basic SNP statistics

Automatically calculate most basic statistics from snpRdata. Calculates
maf, pi, ho, he, pairwise Fst, Fis, HWE divergence, finds private
alleles, and uses Gaussian smoothing to produce per-window averages of
all of these.

## Usage

``` r
calc_basic_snp_stats(
  x,
  facets = NULL,
  fst.method = "WC",
  sigma = NULL,
  step = NULL,
  par = FALSE,
  nk = TRUE
)
```

## Arguments

- x:

  snpRdata object.

- facets:

  character. Categorical metadata variables by which to break up
  analysis. See
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md)
  for more details.

- fst.method:

  character, default "WC". Defines the FST estimator to use. Options:

  - WC: Weir and Cockerham (1984).

  - Weir: Weir (1990)

  - Hohenlohe: Hohenlohe et al (2010), identical to the STACKS package.

  - Genepop: Rousset (2008), uses the genepop package.

- sigma:

  numeric. Designates the width of windows in kilobases. Full window
  size is 6\*sigma.

- step:

  numeric or NULL, default NULL. Designates the number of kilobases
  between each window centroid. If NULL, windows are centered on each
  SNP.

- par:

  numeric or FALSE, default FALSE. If numeric, the number of cores to
  use for parallel processing.

- nk:

  logical, default TRUE. If TRUE, weights SNP contribution to window
  averages by the number of observations at those SNPs.

## Value

A snpRdata object with all of the described statistics merged into the
appropriate sockets.

## Details

The data can be broken up categorically by sample or SNP metadata, as
described in
[`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).
Note that Fst and private allele calculations require a sample specific
contrast (the pairwise part of pairwise Fst), and thus will not be
calculated unless a facet with a sample meta data variable included is
specified. The other stats, in contrast, are snp-specific and thus
ignore any snp meta data variables included in facets. Providing NULL or
"all" to the facets argument works as described in
[`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).

## See also

calc_single_stats calc_pairwise_fst calc_smoothed_averages

## Author

William Hemstrom

## Examples

``` r
x <- calc_basic_snp_stats(stickSNPs, "pop")
get.snpR.stats(x, "pop", stats = "single") # view basic stats
#>     facet subfacet        chr position .snp.id        maf maj.count major
#> 1     pop      ASP    groupVI   212436       1 0.50000000         9     G
#> 2     pop      CLF    groupVI   212436       1 0.70000000         3     G
#> 3     pop      OPL    groupVI   212436       1 0.57500000        17     G
#> 4     pop      PAL    groupVI   212436       1 0.20588235        27     G
#> 5     pop      SMR    groupVI   212436       1 0.59090909        18     G
#> 6     pop      UPD    groupVI   212436       1 0.50000000        13     G
#> 7     pop      ASP groupXVIII   350070       2 0.05000000        19     G
#> 8     pop      CLF groupXVIII   350070       2 0.00000000        12     G
#> 9     pop      OPL groupXVIII   350070       2 0.06521739        43     G
#> 10    pop      PAL groupXVIII   350070       2 0.00000000        34     G
#> 11    pop      SMR groupXVIII   350070       2 0.04545455        42     G
#> 12    pop      UPD groupXVIII   350070       2 0.00000000        30     G
#> 13    pop      ASP    groupXV   854285       3 0.09090909        20     C
#> 14    pop      CLF    groupXV   854285       3 0.08333333        11     C
#> 15    pop      OPL    groupXV   854285       3 0.08000000        46     C
#> 16    pop      PAL    groupXV   854285       3 0.03125000        31     C
#> 17    pop      SMR    groupXV   854285       3 0.06521739        43     C
#> 18    pop      UPD    groupXV   854285       3 0.06666667        28     C
#> 19    pop      ASP  groupXIII   907279       4 0.13636364        19     A
#> 20    pop      CLF  groupXIII   907279       4 0.10000000         9     A
#> 21    pop      OPL  groupXIII   907279       4 0.42500000        23     A
#> 22    pop      PAL  groupXIII   907279       4 0.00000000        32     A
#> 23    pop      SMR  groupXIII   907279       4 0.32608696        31     A
#> 24    pop      UPD  groupXIII   907279       4 0.42307692        15     A
#> 25    pop      ASP    groupVI  1068620       5 0.35714286         9     G
#> 26    pop      CLF    groupVI  1068620       5 0.00000000         6     G
#> 27    pop      OPL    groupVI  1068620       5 0.30769231        18     G
#> 28    pop      PAL    groupVI  1068620       5 0.20833333        19     G
#> 29    pop      SMR    groupVI  1068620       5 0.23333333        23     G
#> 30    pop      UPD    groupVI  1068620       5 0.44444444        10     G
#> 31    pop      ASP     groupV  1392753       6 0.25000000        15     A
#> 32    pop      CLF     groupV  1392753       6 0.25000000         9     A
#> 33    pop      OPL     groupV  1392753       6 0.05769231        49     A
#> 34    pop      PAL     groupV  1392753       6 0.40625000        19     A
#> 35    pop      SMR     groupV  1392753       6 0.07142857        39     A
#> 36    pop      UPD     groupV  1392753       6 0.35714286        18     A
#> 37    pop      ASP    groupIX  2075079       7 0.04545455        21     T
#> 38    pop      CLF    groupIX  2075079       7 0.16666667        10     T
#> 39    pop      OPL    groupIX  2075079       7 0.14285714        36     T
#> 40    pop      PAL    groupIX  2075079       7 0.00000000        34     T
#> 41    pop      SMR    groupIX  2075079       7 0.19047619        34     T
#> 42    pop      UPD    groupIX  2075079       7 0.00000000        28     T
#> 43    pop      ASP  groupXIII  2106000       8 0.22727273        17     G
#> 44    pop      CLF  groupXIII  2106000       8 0.41666667         7     G
#> 45    pop      OPL  groupXIII  2106000       8 0.16666667        40     G
#> 46    pop      PAL  groupXIII  2106000       8 0.00000000        28     G
#> 47    pop      SMR  groupXIII  2106000       8 0.35000000        26     G
#> 48    pop      UPD  groupXIII  2106000       8 0.25000000        15     G
#> 49    pop      ASP   groupXII  2136587       9 0.22727273        17     C
#> 50    pop      CLF   groupXII  2136587       9 0.16666667        10     C
#> 51    pop      OPL   groupXII  2136587       9 0.10000000        45     C
#> 52    pop      PAL   groupXII  2136587       9 0.03125000        31     C
#> 53    pop      SMR   groupXII  2136587       9 0.25000000        30     C
#> 54    pop      UPD   groupXII  2136587       9 0.26923077        19     C
#> 55    pop      ASP  groupXVII  2560839      10 0.09090909        20     A
#> 56    pop      CLF  groupXVII  2560839      10 0.08333333        11     A
#> 57    pop      OPL  groupXVII  2560839      10 0.02000000        49     A
#> 58    pop      PAL  groupXVII  2560839      10 0.00000000        34     A
#> 59    pop      SMR  groupXVII  2560839      10 0.04545455        42     A
#> 60    pop      UPD  groupXVII  2560839      10 0.14285714        24     A
#> 61    pop      ASP   groupIII  2750695      11 0.45000000        11     T
#> 62    pop      CLF   groupIII  2750695      11 0.50000000         6     T
#> 63    pop      OPL   groupIII  2750695      11 0.59523810        17     T
#> 64    pop      PAL   groupIII  2750695      11 0.00000000        32     T
#> 65    pop      SMR   groupIII  2750695      11 0.50000000        19     T
#> 66    pop      UPD   groupIII  2750695      11 0.39285714        17     T
#> 67    pop      ASP   groupXIX  3008524      12 0.25000000        15     A
#> 68    pop      CLF   groupXIX  3008524      12 0.08333333        11     A
#> 69    pop      OPL   groupXIX  3008524      12 0.19444444        29     A
#> 70    pop      PAL   groupXIX  3008524      12 0.31250000        22     A
#> 71    pop      SMR   groupXIX  3008524      12 0.17647059        28     A
#> 72    pop      UPD   groupXIX  3008524      12 0.30000000        14     A
#> 73    pop      ASP   groupXVI  3624417      13 0.11111111        16     T
#> 74    pop      CLF   groupXVI  3624417      13 0.41666667         7     T
#> 75    pop      OPL   groupXVI  3624417      13 0.25000000        36     T
#> 76    pop      PAL   groupXVI  3624417      13 0.03333333        29     T
#> 77    pop      SMR   groupXVI  3624417      13 0.20000000        32     T
#> 78    pop      UPD   groupXVI  3624417      13 0.15384615        22     T
#> 79    pop      ASP   groupVII  3667528      14 0.20000000        16     G
#> 80    pop      CLF   groupVII  3667528      14 0.16666667        10     G
#> 81    pop      OPL   groupVII  3667528      14 0.08000000        46     G
#> 82    pop      PAL   groupVII  3667528      14 0.00000000        34     G
#> 83    pop      SMR   groupVII  3667528      14 0.18181818        36     G
#> 84    pop      UPD   groupVII  3667528      14 0.36666667        19     G
#> 85    pop      ASP   groupVII  4189120      15 0.22727273        17     A
#> 86    pop      CLF   groupVII  4189120      15 0.00000000        12     A
#> 87    pop      OPL   groupVII  4189120      15 0.20000000        32     A
#> 88    pop      PAL   groupVII  4189120      15 0.87500000         4     A
#> 89    pop      SMR   groupVII  4189120      15 0.21428571        33     A
#> 90    pop      UPD   groupVII  4189120      15 0.46153846        14     A
#> 91    pop      ASP    groupXX  4538074      16 0.25000000        15     G
#> 92    pop      CLF    groupXX  4538074      16 0.30000000         7     G
#> 93    pop      OPL    groupXX  4538074      16 0.32500000        27     G
#> 94    pop      PAL    groupXX  4538074      16 0.43333333        17     G
#> 95    pop      SMR    groupXX  4538074      16 0.32500000        27     G
#> 96    pop      UPD    groupXX  4538074      16 0.20833333        19     G
#> 97    pop      ASP  groupXIII  4706089      17 0.15000000        17     A
#> 98    pop      CLF  groupXIII  4706089      17 0.00000000        12     A
#> 99    pop      OPL  groupXIII  4706089      17 0.14000000        43     A
#> 100   pop      PAL  groupXIII  4706089      17 0.18750000        26     A
#> 101   pop      SMR  groupXIII  4706089      17 0.16666667        35     A
#> 102   pop      UPD  groupXIII  4706089      17 0.15384615        22     A
#> 103   pop      ASP   groupIII  4746808      18 0.30000000        14     C
#> 104   pop      CLF   groupIII  4746808      18 0.40000000         6     C
#> 105   pop      OPL   groupIII  4746808      18 0.20370370        43     C
#> 106   pop      PAL   groupIII  4746808      18 0.71428571         8     C
#> 107   pop      SMR   groupIII  4746808      18 0.15217391        39     C
#> 108   pop      UPD   groupIII  4746808      18 0.43333333        17     C
#> 109   pop      ASP  groupXVII  5041995      19 0.00000000        22     T
#> 110   pop      CLF  groupXVII  5041995      19 0.25000000         9     T
#> 111   pop      OPL  groupXVII  5041995      19 0.07500000        37     T
#> 112   pop      PAL  groupXVII  5041995      19 0.00000000        34     T
#> 113   pop      SMR  groupXVII  5041995      19 0.06521739        43     T
#> 114   pop      UPD  groupXVII  5041995      19 0.06666667        28     T
#> 115   pop      ASP     groupX  5161637      20 0.05555556        17     G
#> 116   pop      CLF     groupX  5161637      20 0.08333333        11     G
#> 117   pop      OPL     groupX  5161637      20 0.02083333        47     G
#> 118   pop      PAL     groupX  5161637      20 0.00000000        30     G
#> 119   pop      SMR     groupX  5161637      20 0.06521739        43     G
#> 120   pop      UPD     groupX  5161637      20 0.07142857        26     G
#> 121   pop      ASP   groupXII  5178037      21 0.16666667        15     C
#> 122   pop      CLF   groupXII  5178037      21 0.37500000         5     C
#> 123   pop      OPL   groupXII  5178037      21 0.50000000        19     C
#> 124   pop      PAL   groupXII  5178037      21 0.00000000        32     C
#> 125   pop      SMR   groupXII  5178037      21 0.32500000        27     C
#> 126   pop      UPD   groupXII  5178037      21 0.33333333        16     C
#> 127   pop      ASP    groupIX  5287711      22 0.44444444        10     C
#> 128   pop      CLF    groupIX  5287711      22 0.40000000         6     C
#> 129   pop      OPL    groupIX  5287711      22 0.30952381        29     C
#> 130   pop      PAL    groupIX  5287711      22 0.00000000        34     C
#> 131   pop      SMR    groupIX  5287711      22 0.15000000        34     C
#> 132   pop      UPD    groupIX  5287711      22 0.35714286        18     C
#> 133   pop      ASP  groupVIII  5462135      23 0.05000000        19     G
#> 134   pop      CLF  groupVIII  5462135      23 0.00000000        12     G
#> 135   pop      OPL  groupVIII  5462135      23 0.04000000        48     G
#> 136   pop      PAL  groupVIII  5462135      23 0.00000000        32     G
#> 137   pop      SMR  groupVIII  5462135      23 0.10000000        36     G
#> 138   pop      UPD  groupVIII  5462135      23 0.03333333        29     G
#> 139   pop      ASP    groupIX  5960933      24 0.05000000        19     A
#> 140   pop      CLF    groupIX  5960933      24 0.00000000        12     A
#> 141   pop      OPL    groupIX  5960933      24 0.09615385        47     A
#> 142   pop      PAL    groupIX  5960933      24 0.00000000        32     A
#> 143   pop      SMR    groupIX  5960933      24 0.05263158        36     A
#> 144   pop      UPD    groupIX  5960933      24 0.10714286        25     A
#> 145   pop      ASP  groupXIII  6076868      25 0.30000000        14     G
#> 146   pop      CLF  groupXIII  6076868      25 0.10000000         9     G
#> 147   pop      OPL  groupXIII  6076868      25 0.26086957        34     G
#> 148   pop      PAL  groupXIII  6076868      25 0.00000000        32     G
#> 149   pop      SMR  groupXIII  6076868      25 0.50000000        19     G
#> 150   pop      UPD  groupXIII  6076868      25 0.43333333        17     G
#> 151   pop      ASP    groupXI  6177332      26 0.33333333        12     T
#> 152   pop      CLF    groupXI  6177332      26 0.00000000        12     T
#> 153   pop      OPL    groupXI  6177332      26 0.14285714        36     T
#> 154   pop      PAL    groupXI  6177332      26 0.53333333        14     T
#> 155   pop      SMR    groupXI  6177332      26 0.29545455        31     T
#> 156   pop      UPD    groupXI  6177332      26 0.29166667        17     T
#> 157   pop      ASP    groupIV  6291287      27 0.05000000        19     A
#> 158   pop      CLF    groupIV  6291287      27 0.08333333        11     A
#> 159   pop      OPL    groupIV  6291287      27 0.20833333        38     A
#> 160   pop      PAL    groupIV  6291287      27 0.00000000        34     A
#> 161   pop      SMR    groupIV  6291287      27 0.15000000        34     A
#> 162   pop      UPD    groupIV  6291287      27 0.07142857        26     A
#> 163   pop      ASP    groupIX  6515796      28 0.12500000        14     T
#> 164   pop      CLF    groupIX  6515796      28 0.16666667        10     T
#> 165   pop      OPL    groupIX  6515796      28 0.65789474        13     T
#> 166   pop      PAL    groupIX  6515796      28 0.10000000        27     T
#> 167   pop      SMR    groupIX  6515796      28 0.52500000        19     T
#> 168   pop      UPD    groupIX  6515796      28 0.28571429        20     T
#> 169   pop      ASP   groupXII  6718079      29 0.05555556        17     A
#> 170   pop      CLF   groupXII  6718079      29 0.12500000         7     A
#> 171   pop      OPL   groupXII  6718079      29 0.26190476        31     A
#> 172   pop      PAL   groupXII  6718079      29 0.16666667        25     A
#> 173   pop      SMR   groupXII  6718079      29 0.23333333        23     A
#> 174   pop      UPD   groupXII  6718079      29 0.11538462        23     A
#> 175   pop      ASP   groupVII  6967477      30 0.09090909        20     G
#> 176   pop      CLF   groupVII  6967477      30 0.20000000         8     G
#> 177   pop      OPL   groupVII  6967477      30 0.07142857        39     G
#> 178   pop      PAL   groupVII  6967477      30 0.00000000        34     G
#> 179   pop      SMR   groupVII  6967477      30 0.11363636        39     G
#> 180   pop      UPD   groupVII  6967477      30 0.16666667        20     G
#> 181   pop      ASP   groupIII  7194069      31 0.27272727        16     C
#> 182   pop      CLF   groupIII  7194069      31 0.08333333        11     C
#> 183   pop      OPL   groupIII  7194069      31 0.32000000        34     C
#> 184   pop      PAL   groupIII  7194069      31 0.00000000        34     C
#> 185   pop      SMR   groupIII  7194069      31 0.36842105        24     C
#> 186   pop      UPD   groupIII  7194069      31 0.43333333        17     C
#> 187   pop      ASP   groupXVI  7199722      32 0.15000000        17     A
#> 188   pop      CLF   groupXVI  7199722      32 0.20000000         8     A
#> 189   pop      OPL   groupXVI  7199722      32 0.13043478        40     A
#> 190   pop      PAL   groupXVI  7199722      32 0.40000000        18     A
#> 191   pop      SMR   groupXVI  7199722      32 0.17500000        33     A
#> 192   pop      UPD   groupXVI  7199722      32 0.26923077        19     A
#> 193   pop      ASP  groupVIII  7224065      33 0.27272727        16     C
#> 194   pop      CLF  groupVIII  7224065      33 0.37500000         5     C
#> 195   pop      OPL  groupVIII  7224065      33 0.36842105        24     C
#> 196   pop      PAL  groupVIII  7224065      33 0.23333333        23     C
#> 197   pop      SMR  groupVIII  7224065      33 0.34090909        29     C
#> 198   pop      UPD  groupVIII  7224065      33 0.34615385        17     C
#> 199   pop      ASP groupXVIII  7229426      34 0.38888889        11     C
#> 200   pop      CLF groupXVIII  7229426      34 0.50000000         6     C
#> 201   pop      OPL groupXVIII  7229426      34 0.38888889        22     C
#> 202   pop      PAL groupXVIII  7229426      34 0.00000000        34     C
#> 203   pop      SMR groupXVIII  7229426      34 0.29545455        31     C
#> 204   pop      UPD groupXVIII  7229426      34 0.13636364        19     C
#> 205   pop      ASP groupXVIII  8087022      35 0.10000000        18     A
#> 206   pop      CLF groupXVIII  8087022      35 0.00000000        10     A
#> 207   pop      OPL groupXVIII  8087022      35 0.06818182        41     A
#> 208   pop      PAL groupXVIII  8087022      35 0.00000000        32     A
#> 209   pop      SMR groupXVIII  8087022      35 0.07500000        37     A
#> 210   pop      UPD groupXVIII  8087022      35 0.07142857        26     A
#> 211   pop      ASP    groupIX  8373442      36 0.15000000        17     A
#> 212   pop      CLF    groupIX  8373442      36 0.16666667        10     A
#> 213   pop      OPL    groupIX  8373442      36 0.18750000        39     A
#> 214   pop      PAL    groupIX  8373442      36 0.00000000        32     A
#> 215   pop      SMR    groupIX  8373442      36 0.18181818        36     A
#> 216   pop      UPD    groupIX  8373442      36 0.07692308        24     A
#> 217   pop      ASP    groupXV  8815538      37 0.20000000        16     G
#> 218   pop      CLF    groupXV  8815538      37 0.20000000         8     G
#> 219   pop      OPL    groupXV  8815538      37 0.27083333        35     G
#> 220   pop      PAL    groupXV  8815538      37 0.03125000        31     G
#> 221   pop      SMR    groupXV  8815538      37 0.36363636        28     G
#> 222   pop      UPD    groupXV  8815538      37 0.26923077        19     G
#> 223   pop      ASP    groupXI  9019877      38 0.45000000        11     G
#> 224   pop      CLF    groupXI  9019877      38 0.00000000        12     G
#> 225   pop      OPL    groupXI  9019877      38 0.19565217        37     G
#> 226   pop      PAL    groupXI  9019877      38 0.40625000        19     G
#> 227   pop      SMR    groupXI  9019877      38 0.42105263        22     G
#> 228   pop      UPD    groupXI  9019877      38 0.25000000        21     G
#> 229   pop      ASP   groupXXI  9067876      39 0.16666667        15     T
#> 230   pop      CLF   groupXXI  9067876      39 0.33333333         4     T
#> 231   pop      OPL   groupXXI  9067876      39 0.10000000        27     T
#> 232   pop      PAL   groupXXI  9067876      39 0.13636364        19     T
#> 233   pop      SMR   groupXXI  9067876      39 0.08823529        31     T
#> 234   pop      UPD   groupXXI  9067876      39 0.15000000        17     T
#> 235   pop      ASP  groupXVII  9108184      40 0.45000000        11     G
#> 236   pop      CLF  groupXVII  9108184      40 0.20000000         8     G
#> 237   pop      OPL  groupXVII  9108184      40 0.17500000        33     G
#> 238   pop      PAL  groupXVII  9108184      40 0.16666667        25     G
#> 239   pop      SMR  groupXVII  9108184      40 0.15789474        32     G
#> 240   pop      UPD  groupXVII  9108184      40 0.40909091        13     G
#> 241   pop      ASP    groupIX  9284427      41 0.59090909         9     G
#> 242   pop      CLF    groupIX  9284427      41 0.30000000         7     G
#> 243   pop      OPL    groupIX  9284427      41 0.26315789        28     G
#> 244   pop      PAL    groupIX  9284427      41 1.00000000         0     G
#> 245   pop      SMR    groupIX  9284427      41 0.32352941        23     G
#> 246   pop      UPD    groupIX  9284427      41 0.60000000         8     G
#> 247   pop      ASP     groupV  9767769      42 0.25000000        12     G
#> 248   pop      CLF     groupV  9767769      42 0.40000000         6     G
#> 249   pop      OPL     groupV  9767769      42 0.45652174        25     G
#> 250   pop      PAL     groupV  9767769      42 0.06250000        30     G
#> 251   pop      SMR     groupV  9767769      42 0.56818182        19     G
#> 252   pop      UPD     groupV  9767769      42 0.35714286        18     G
#> 253   pop      ASP   groupXII  9881312      43 0.09090909        20     A
#> 254   pop      CLF   groupXII  9881312      43 0.20000000         8     A
#> 255   pop      OPL   groupXII  9881312      43 0.19047619        34     A
#> 256   pop      PAL   groupXII  9881312      43 0.00000000        34     A
#> 257   pop      SMR   groupXII  9881312      43 0.10416667        43     A
#> 258   pop      UPD   groupXII  9881312      43 0.23076923        20     A
#> 259   pop      ASP    groupXX 10296516      44 0.16666667        15     T
#> 260   pop      CLF    groupXX 10296516      44 0.00000000        12     T
#> 261   pop      OPL    groupXX 10296516      44 0.05769231        49     T
#> 262   pop      PAL    groupXX 10296516      44 0.40000000        18     T
#> 263   pop      SMR    groupXX 10296516      44 0.11904762        37     T
#> 264   pop      UPD    groupXX 10296516      44 0.20833333        19     T
#> 265   pop      ASP     groupI 10593829      45 0.00000000        22     G
#> 266   pop      CLF     groupI 10593829      45 0.00000000        12     G
#> 267   pop      OPL     groupI 10593829      45 0.00000000        38     G
#> 268   pop      PAL     groupI 10593829      45 0.00000000        34     G
#> 269   pop      SMR     groupI 10593829      45 0.00000000        40     G
#> 270   pop      UPD     groupI 10593829      45 0.06666667        28     G
#> 271   pop      ASP   groupXIX 10765968      46 0.43750000         9     C
#> 272   pop      CLF   groupXIX 10765968      46 0.20000000         8     C
#> 273   pop      OPL   groupXIX 10765968      46 0.41666667        21     C
#> 274   pop      PAL   groupXIX 10765968      46 0.00000000        32     C
#> 275   pop      SMR   groupXIX 10765968      46 0.20588235        27     C
#> 276   pop      UPD   groupXIX 10765968      46 0.34615385        17     C
#> 277   pop      ASP    groupII 10829169      47 0.22727273        17     C
#> 278   pop      CLF    groupII 10829169      47 0.00000000        12     C
#> 279   pop      OPL    groupII 10829169      47 0.47826087        24     C
#> 280   pop      PAL    groupII 10829169      47 0.10000000        27     C
#> 281   pop      SMR    groupII 10829169      47 0.23913043        35     C
#> 282   pop      UPD    groupII 10829169      47 0.16666667        25     C
#> 283   pop      ASP    groupXX 10913115      48 0.27777778        13     T
#> 284   pop      CLF    groupXX 10913115      48 0.50000000         6     T
#> 285   pop      OPL    groupXX 10913115      48 0.57500000        17     T
#> 286   pop      PAL    groupXX 10913115      48 0.39285714        17     T
#> 287   pop      SMR    groupXX 10913115      48 0.44444444        20     T
#> 288   pop      UPD    groupXX 10913115      48 0.38461538        16     T
#> 289   pop      ASP groupXVIII 10974594      49 0.00000000        20     G
#> 290   pop      CLF groupXVIII 10974594      49 0.16666667        10     G
#> 291   pop      OPL groupXVIII 10974594      49 0.01851852        53     G
#> 292   pop      PAL groupXVIII 10974594      49 0.00000000        34     G
#> 293   pop      SMR groupXVIII 10974594      49 0.00000000        48     G
#> 294   pop      UPD groupXVIII 10974594      49 0.07692308        24     G
#> 295   pop      ASP  groupXVII 11256661      50 0.15000000        17     C
#> 296   pop      CLF  groupXVII 11256661      50 0.10000000         9     C
#> 297   pop      OPL  groupXVII 11256661      50 0.10869565        41     C
#> 298   pop      PAL  groupXVII 11256661      50 0.03125000        31     C
#> 299   pop      SMR  groupXVII 11256661      50 0.05263158        36     C
#> 300   pop      UPD  groupXVII 11256661      50 0.07142857        26     C
#> 301   pop      ASP  groupVIII 11438870      51 0.38888889        11     T
#> 302   pop      CLF  groupVIII 11438870      51 0.00000000        10     T
#> 303   pop      OPL  groupVIII 11438870      51 0.35000000        26     T
#> 304   pop      PAL  groupVIII 11438870      51 0.18181818        18     T
#> 305   pop      SMR  groupVIII 11438870      51 0.32142857        19     T
#> 306   pop      UPD  groupVIII 11438870      51 0.39285714        17     T
#> 307   pop      ASP   groupVII 11477223      52 0.36363636        14     A
#> 308   pop      CLF   groupVII 11477223      52 0.33333333         8     A
#> 309   pop      OPL   groupVII 11477223      52 0.33333333        32     A
#> 310   pop      PAL   groupVII 11477223      52 0.45833333        13     A
#> 311   pop      SMR   groupVII 11477223      52 0.14583333        41     A
#> 312   pop      UPD   groupVII 11477223      52 0.35714286        18     A
#> 313   pop      ASP    groupII 11505084      53 0.20000000        16     A
#> 314   pop      CLF    groupII 11505084      53 0.37500000         5     A
#> 315   pop      OPL    groupII 11505084      53 0.30952381        29     A
#> 316   pop      PAL    groupII 11505084      53 0.00000000        32     A
#> 317   pop      SMR    groupII 11505084      53 0.17391304        38     A
#> 318   pop      UPD    groupII 11505084      53 0.41666667        14     A
#> 319   pop      ASP    groupIX 11904651      54 0.00000000        22     G
#> 320   pop      CLF    groupIX 11904651      54 0.08333333        11     G
#> 321   pop      OPL    groupIX 11904651      54 0.01923077        51     G
#> 322   pop      PAL    groupIX 11904651      54 0.00000000        32     G
#> 323   pop      SMR    groupIX 11904651      54 0.02173913        45     G
#> 324   pop      UPD    groupIX 11904651      54 0.00000000        30     G
#> 325   pop      ASP  groupVIII 12108945      55 0.50000000         8     G
#> 326   pop      CLF  groupVIII 12108945      55 0.25000000         9     G
#> 327   pop      OPL  groupVIII 12108945      55 0.33333333        28     G
#> 328   pop      PAL  groupVIII 12108945      55 0.46666667        16     G
#> 329   pop      SMR  groupVIII 12108945      55 0.17500000        33     G
#> 330   pop      UPD  groupVIII 12108945      55 0.46153846        14     G
#> 331   pop      ASP   groupXIX 12485182      56 0.33333333        12     A
#> 332   pop      CLF   groupXIX 12485182      56 0.08333333        11     A
#> 333   pop      OPL   groupXIX 12485182      56 0.20000000        32     A
#> 334   pop      PAL   groupXIX 12485182      56 0.34375000        21     A
#> 335   pop      SMR   groupXIX 12485182      56 0.18421053        31     A
#> 336   pop      UPD   groupXIX 12485182      56 0.23333333        23     A
#> 337   pop      ASP    groupXI 12833733      57 0.06250000        15     G
#> 338   pop      CLF    groupXI 12833733      57 0.00000000        10     G
#> 339   pop      OPL    groupXI 12833733      57 0.10000000        36     G
#> 340   pop      PAL    groupXI 12833733      57 0.00000000        24     G
#> 341   pop      SMR    groupXI 12833733      57 0.02777778        35     G
#> 342   pop      UPD    groupXI 12833733      57 0.00000000        20     G
#> 343   pop      ASP groupXVIII 12937786      58 0.44444444        10     C
#> 344   pop      CLF groupXVIII 12937786      58 0.00000000        12     C
#> 345   pop      OPL groupXVIII 12937786      58 0.26000000        37     C
#> 346   pop      PAL groupXVIII 12937786      58 0.42307692        15     C
#> 347   pop      SMR groupXVIII 12937786      58 0.13888889        31     C
#> 348   pop      UPD groupXVIII 12937786      58 0.27272727        16     C
#> 349   pop      ASP   groupIII 13177027      59 0.45454545        12     G
#> 350   pop      CLF   groupIII 13177027      59 0.33333333         8     G
#> 351   pop      OPL   groupIII 13177027      59 0.39583333        29     G
#> 352   pop      PAL   groupIII 13177027      59 0.06250000        30     G
#> 353   pop      SMR   groupIII 13177027      59 0.47500000        21     G
#> 354   pop      UPD   groupIII 13177027      59 0.38461538        16     G
#> 355   pop      ASP    groupVI 13335577      60 0.42857143         8     A
#> 356   pop      CLF    groupVI 13335577      60 0.75000000         3     A
#> 357   pop      OPL    groupVI 13335577      60 0.34782609        30     A
#> 358   pop      PAL    groupVI 13335577      60 0.00000000        34     A
#> 359   pop      SMR    groupVI 13335577      60 0.33333333        28     A
#> 360   pop      UPD    groupVI 13335577      60 0.32142857        19     A
#> 361   pop      ASP     groupX 13718292      61 0.00000000        16     T
#> 362   pop      CLF     groupX 13718292      61 0.16666667         5     T
#> 363   pop      OPL     groupX 13718292      61 0.02941176        33     T
#> 364   pop      PAL     groupX 13718292      61 0.00000000        26     T
#> 365   pop      SMR     groupX 13718292      61 0.05000000        38     T
#> 366   pop      UPD     groupX 13718292      61 0.08333333        22     T
#> 367   pop      ASP    groupVI 13895716      62 0.10000000        18     G
#> 368   pop      CLF    groupVI 13895716      62 0.00000000        12     G
#> 369   pop      OPL    groupVI 13895716      62 0.21153846        41     G
#> 370   pop      PAL    groupVI 13895716      62 0.00000000        32     G
#> 371   pop      SMR    groupVI 13895716      62 0.06250000        45     G
#> 372   pop      UPD    groupVI 13895716      62 0.10714286        25     G
#> 373   pop      ASP   groupXIV 14009990      63 0.37500000        10     A
#> 374   pop      CLF   groupXIV 14009990      63 0.40000000         6     A
#> 375   pop      OPL   groupXIV 14009990      63 0.22222222        28     A
#> 376   pop      PAL   groupXIV 14009990      63 0.80000000         4     A
#> 377   pop      SMR   groupXIV 14009990      63 0.31250000        22     A
#> 378   pop      UPD   groupXIV 14009990      63 0.26923077        19     A
#> 379   pop      ASP   groupXII 14158350      64 0.31818182        15     T
#> 380   pop      CLF   groupXII 14158350      64 0.16666667        10     T
#> 381   pop      OPL   groupXII 14158350      64 0.30952381        29     T
#> 382   pop      PAL   groupXII 14158350      64 0.57142857        12     T
#> 383   pop      SMR   groupXII 14158350      64 0.36363636        28     T
#> 384   pop      UPD   groupXII 14158350      64 0.36666667        19     T
#> 385   pop      ASP    groupIX 14164215      65 0.18181818        18     C
#> 386   pop      CLF    groupIX 14164215      65 0.41666667         7     C
#> 387   pop      OPL    groupIX 14164215      65 0.63157895        14     C
#> 388   pop      PAL    groupIX 14164215      65 0.26470588        25     C
#> 389   pop      SMR    groupIX 14164215      65 0.34090909        29     C
#> 390   pop      UPD    groupIX 14164215      65 0.39285714        17     C
#> 391   pop      ASP    groupVI 14236496      66 0.50000000         9     A
#> 392   pop      CLF    groupVI 14236496      66 0.41666667         7     A
#> 393   pop      OPL    groupVI 14236496      66 0.35000000        26     A
#> 394   pop      PAL    groupVI 14236496      66 0.65625000        11     A
#> 395   pop      SMR    groupVI 14236496      66 0.38235294        21     A
#> 396   pop      UPD    groupVI 14236496      66 0.42307692        15     A
#> 397   pop      ASP     groupX 14273796      67 0.00000000        18     G
#> 398   pop      CLF     groupX 14273796      67 0.00000000        12     G
#> 399   pop      OPL     groupX 14273796      67 0.00000000        50     G
#> 400   pop      PAL     groupX 14273796      67 0.26666667        22     G
#> 401   pop      SMR     groupX 14273796      67 0.00000000        44     G
#> 402   pop      UPD     groupX 14273796      67 0.03846154        25     G
#> 403   pop      ASP  groupVIII 14445389      68 0.25000000        15     A
#> 404   pop      CLF  groupVIII 14445389      68 0.25000000         9     A
#> 405   pop      OPL  groupVIII 14445389      68 0.22916667        37     A
#> 406   pop      PAL  groupVIII 14445389      68 0.73076923         7     A
#> 407   pop      SMR  groupVIII 14445389      68 0.27777778        26     A
#> 408   pop      UPD  groupVIII 14445389      68 0.35714286        18     A
#> 409   pop      ASP    groupIX 14868533      69 0.00000000        20     A
#> 410   pop      CLF    groupIX 14868533      69 0.10000000         9     A
#> 411   pop      OPL    groupIX 14868533      69 0.30000000        35     A
#> 412   pop      PAL    groupIX 14868533      69 0.00000000        34     A
#> 413   pop      SMR    groupIX 14868533      69 0.09523810        38     A
#> 414   pop      UPD    groupIX 14868533      69 0.00000000        30     A
#> 415   pop      ASP    groupIX 15190693      70 0.35000000        13     A
#> 416   pop      CLF    groupIX 15190693      70 0.60000000         4     A
#> 417   pop      OPL    groupIX 15190693      70 0.26190476        31     A
#> 418   pop      PAL    groupIX 15190693      70 0.00000000        32     A
#> 419   pop      SMR    groupIX 15190693      70 0.25000000        33     A
#> 420   pop      UPD    groupIX 15190693      70 0.29166667        17     A
#> 421   pop      ASP   groupXII 15306362      71 0.04545455        21     A
#> 422   pop      CLF   groupXII 15306362      71 0.00000000        12     A
#> 423   pop      OPL   groupXII 15306362      71 0.01851852        53     A
#> 424   pop      PAL   groupXII 15306362      71 0.18750000        26     A
#> 425   pop      SMR   groupXII 15306362      71 0.12500000        42     A
#> 426   pop      UPD   groupXII 15306362      71 0.06666667        28     A
#> 427   pop      ASP    groupXX 15460430      72 0.65000000         7     A
#> 428   pop      CLF    groupXX 15460430      72 0.58333333         5     A
#> 429   pop      OPL    groupXX 15460430      72 0.31250000        33     A
#> 430   pop      PAL    groupXX 15460430      72 0.63333333        11     A
#> 431   pop      SMR    groupXX 15460430      72 0.22727273        34     A
#> 432   pop      UPD    groupXX 15460430      72 0.73076923         7     A
#> 433   pop      ASP   groupVII 15493025      73 0.30000000        14     T
#> 434   pop      CLF   groupVII 15493025      73 0.25000000         9     T
#> 435   pop      OPL   groupVII 15493025      73 0.13043478        40     T
#> 436   pop      PAL   groupVII 15493025      73 0.30000000        21     T
#> 437   pop      SMR   groupVII 15493025      73 0.11904762        37     T
#> 438   pop      UPD   groupVII 15493025      73 0.25000000        21     T
#> 439   pop      ASP    groupXV 15547073      74 0.15000000        17     G
#> 440   pop      CLF    groupXV 15547073      74 0.10000000         9     G
#> 441   pop      OPL    groupXV 15547073      74 0.31578947        26     G
#> 442   pop      PAL    groupXV 15547073      74 0.00000000        34     G
#> 443   pop      SMR    groupXV 15547073      74 0.23529412        26     G
#> 444   pop      UPD    groupXV 15547073      74 0.11538462        23     G
#> 445   pop      ASP  groupXIII 15697285      75 0.55555556         8     A
#> 446   pop      CLF  groupXIII 15697285      75 0.50000000         5     A
#> 447   pop      OPL  groupXIII 15697285      75 0.43478261        26     A
#> 448   pop      PAL  groupXIII 15697285      75 0.31250000        22     A
#> 449   pop      SMR  groupXIII 15697285      75 0.47500000        21     A
#> 450   pop      UPD  groupXIII 15697285      75 0.70000000         9     A
#> 451   pop      ASP   groupIII 15983994      76 0.41666667         7     C
#> 452   pop      CLF   groupIII 15983994      76 0.16666667         5     C
#> 453   pop      OPL   groupIII 15983994      76 0.34375000        21     C
#> 454   pop      PAL   groupIII 15983994      76 1.00000000         0     C
#> 455   pop      SMR   groupIII 15983994      76 0.12500000        21     C
#> 456   pop      UPD   groupIII 15983994      76 0.42857143        16     C
#> 457   pop      ASP   groupXIX 16233699      77 0.33333333        12     A
#> 458   pop      CLF   groupXIX 16233699      77 0.08333333        11     A
#> 459   pop      OPL   groupXIX 16233699      77 0.17391304        38     A
#> 460   pop      PAL   groupXIX 16233699      77 0.36666667        19     A
#> 461   pop      SMR   groupXIX 16233699      77 0.18421053        31     A
#> 462   pop      UPD   groupXIX 16233699      77 0.21428571        22     A
#> 463   pop      ASP    groupIX 16306966      78 0.00000000        20     G
#> 464   pop      CLF    groupIX 16306966      78 0.25000000         9     G
#> 465   pop      OPL    groupIX 16306966      78 0.29545455        31     G
#> 466   pop      PAL    groupIX 16306966      78 0.00000000        34     G
#> 467   pop      SMR    groupIX 16306966      78 0.11363636        39     G
#> 468   pop      UPD    groupIX 16306966      78 0.00000000        30     G
#> 469   pop      ASP     groupI 16498302      79 0.18181818        18     T
#> 470   pop      CLF     groupI 16498302      79 0.00000000        12     T
#> 471   pop      OPL     groupI 16498302      79 0.11538462        46     T
#> 472   pop      PAL     groupI 16498302      79 0.23529412        26     T
#> 473   pop      SMR     groupI 16498302      79 0.09090909        40     T
#> 474   pop      UPD     groupI 16498302      79 0.12500000        21     T
#> 475   pop      ASP  groupXIII 16520929      80 0.22727273        17     T
#> 476   pop      CLF  groupXIII 16520929      80 0.16666667        10     T
#> 477   pop      OPL  groupXIII 16520929      80 0.31818182        30     T
#> 478   pop      PAL  groupXIII 16520929      80 0.50000000        15     T
#> 479   pop      SMR  groupXIII 16520929      80 0.26086957        34     T
#> 480   pop      UPD  groupXIII 16520929      80 0.26923077        19     T
#> 481   pop      ASP   groupXVI 16561891      81 0.33333333        12     A
#> 482   pop      CLF   groupXVI 16561891      81 0.60000000         4     A
#> 483   pop      OPL   groupXVI 16561891      81 0.40476190        25     A
#> 484   pop      PAL   groupXVI 16561891      81 0.43750000        18     A
#> 485   pop      SMR   groupXVI 16561891      81 0.31818182        30     A
#> 486   pop      UPD   groupXVI 16561891      81 0.39285714        17     A
#> 487   pop      ASP    groupII 16598364      82 0.10000000        18     T
#> 488   pop      CLF    groupII 16598364      82 0.00000000         6     T
#> 489   pop      OPL    groupII 16598364      82 0.09615385        47     T
#> 490   pop      PAL    groupII 16598364      82 0.00000000        34     T
#> 491   pop      SMR    groupII 16598364      82 0.05263158        36     T
#> 492   pop      UPD    groupII 16598364      82 0.19230769        21     T
#> 493   pop      ASP    groupIV 17248473      83 0.18181818        18     T
#> 494   pop      CLF    groupIV 17248473      83 0.41666667         7     T
#> 495   pop      OPL    groupIV 17248473      83 0.00000000        52     T
#> 496   pop      PAL    groupIV 17248473      83 0.00000000        32     T
#> 497   pop      SMR    groupIV 17248473      83 0.02083333        47     T
#> 498   pop      UPD    groupIV 17248473      83 0.03571429        27     T
#> 499   pop      ASP    groupXX 17748944      84 0.30000000        14     C
#> 500   pop      CLF    groupXX 17748944      84 0.33333333         8     C
#> 501   pop      OPL    groupXX 17748944      84 0.35000000        26     C
#> 502   pop      PAL    groupXX 17748944      84 1.00000000         0     C
#> 503   pop      SMR    groupXX 17748944      84 0.26315789        28     C
#> 504   pop      UPD    groupXX 17748944      84 0.46153846        14     C
#> 505   pop      ASP  groupXIII 17842228      85 0.65000000         7     G
#> 506   pop      CLF  groupXIII 17842228      85 0.60000000         4     G
#> 507   pop      OPL  groupXIII 17842228      85 0.30952381        29     G
#> 508   pop      PAL  groupXIII 17842228      85 0.89285714         3     G
#> 509   pop      SMR  groupXIII 17842228      85 0.26190476        31     G
#> 510   pop      UPD  groupXIII 17842228      85 0.43333333        17     G
#> 511   pop      ASP   groupXIX 17858769      86 0.37500000        10     A
#> 512   pop      CLF   groupXIX 17858769      86 0.08333333        11     A
#> 513   pop      OPL   groupXIX 17858769      86 0.22500000        31     A
#> 514   pop      PAL   groupXIX 17858769      86 0.40909091        13     A
#> 515   pop      SMR   groupXIX 17858769      86 0.18421053        31     A
#> 516   pop      UPD   groupXIX 17858769      86 0.26923077        19     A
#> 517   pop      ASP    groupIX 17929225      87 0.50000000        11     T
#> 518   pop      CLF    groupIX 17929225      87 0.08333333        11     T
#> 519   pop      OPL    groupIX 17929225      87 0.24000000        38     T
#> 520   pop      PAL    groupIX 17929225      87 0.00000000        32     T
#> 521   pop      SMR    groupIX 17929225      87 0.27272727        32     T
#> 522   pop      UPD    groupIX 17929225      87 0.42857143        16     T
#> 523   pop      ASP    groupIX 18339767      88 0.50000000         8     A
#> 524   pop      CLF    groupIX 18339767      88 0.25000000         9     A
#> 525   pop      OPL    groupIX 18339767      88 0.32500000        27     A
#> 526   pop      PAL    groupIX 18339767      88 1.00000000         0     A
#> 527   pop      SMR    groupIX 18339767      88 0.40625000        19     A
#> 528   pop      UPD    groupIX 18339767      88 0.33333333        16     A
#> 529   pop      ASP  groupXIII 18836558      89 0.33333333        12     C
#> 530   pop      CLF  groupXIII 18836558      89 0.33333333         4     C
#> 531   pop      OPL  groupXIII 18836558      89 0.36111111        23     C
#> 532   pop      PAL  groupXIII 18836558      89 0.50000000        12     C
#> 533   pop      SMR  groupXIII 18836558      89 0.36111111        23     C
#> 534   pop      UPD  groupXIII 18836558      89 0.58333333        10     C
#> 535   pop      ASP    groupII 19244005      90 0.13636364        19     C
#> 536   pop      CLF    groupII 19244005      90 0.00000000        12     C
#> 537   pop      OPL    groupII 19244005      90 0.00000000        52     C
#> 538   pop      PAL    groupII 19244005      90 0.00000000        34     C
#> 539   pop      SMR    groupII 19244005      90 0.02173913        45     C
#> 540   pop      UPD    groupII 19244005      90 0.07692308        24     C
#> 541   pop      ASP    groupXX 19317273      91 0.16666667        15     C
#> 542   pop      CLF    groupXX 19317273      91 0.37500000         5     C
#> 543   pop      OPL    groupXX 19317273      91 0.50000000        21     C
#> 544   pop      PAL    groupXX 19317273      91 0.09375000        29     C
#> 545   pop      SMR    groupXX 19317273      91 0.55000000        18     C
#> 546   pop      UPD    groupXX 19317273      91 0.35714286        18     C
#> 547   pop      ASP     groupI 19759331      92 0.50000000         7     G
#> 548   pop      CLF     groupI 19759331      92 0.20000000         8     G
#> 549   pop      OPL     groupI 19759331      92 0.35714286        18     G
#> 550   pop      PAL     groupI 19759331      92 0.63636364         8     G
#> 551   pop      SMR     groupI 19759331      92 0.18181818        18     G
#> 552   pop      UPD     groupI 19759331      92 0.20000000        16     G
#> 553   pop      ASP   groupVII 20165240      93 0.13636364        19     T
#> 554   pop      CLF   groupVII 20165240      93 0.25000000         9     T
#> 555   pop      OPL   groupVII 20165240      93 0.05769231        49     T
#> 556   pop      PAL   groupVII 20165240      93 0.14285714        24     T
#> 557   pop      SMR   groupVII 20165240      93 0.18421053        31     T
#> 558   pop      UPD   groupVII 20165240      93 0.25000000        21     T
#> 559   pop      ASP    groupII 20496069      94 0.00000000        16     C
#> 560   pop      CLF    groupII 20496069      94 0.00000000        10     C
#> 561   pop      OPL    groupII 20496069      94 0.06818182        41     C
#> 562   pop      PAL    groupII 20496069      94 0.34615385        17     C
#> 563   pop      SMR    groupII 20496069      94 0.02272727        43     C
#> 564   pop      UPD    groupII 20496069      94 0.10714286        25     C
#> 565   pop      ASP     groupI 21284266      95 0.20000000        16     G
#> 566   pop      CLF     groupI 21284266      95 0.33333333         8     G
#> 567   pop      OPL     groupI 21284266      95 0.18181818        36     G
#> 568   pop      PAL     groupI 21284266      95 0.00000000        34     G
#> 569   pop      SMR     groupI 21284266      95 0.15789474        32     G
#> 570   pop      UPD     groupI 21284266      95 0.16666667        20     G
#> 571   pop      ASP    groupIV 22227311      96 0.06250000        15     C
#> 572   pop      CLF    groupIV 22227311      96 0.08333333        11     C
#> 573   pop      OPL    groupIV 22227311      96 0.17391304        38     C
#> 574   pop      PAL    groupIV 22227311      96 0.00000000        34     C
#> 575   pop      SMR    groupIV 22227311      96 0.08695652        42     C
#> 576   pop      UPD    groupIV 22227311      96 0.07692308        24     C
#> 577   pop      ASP     groupI 22746027      97 0.00000000        22     G
#> 578   pop      CLF     groupI 22746027      97 0.08333333        11     G
#> 579   pop      OPL     groupI 22746027      97 0.04761905        40     G
#> 580   pop      PAL     groupI 22746027      97 0.00000000        34     G
#> 581   pop      SMR     groupI 22746027      97 0.05000000        38     G
#> 582   pop      UPD     groupI 22746027      97 0.03333333        29     G
#> 583   pop      ASP   groupVII 25239495      98 0.22222222        14     G
#> 584   pop      CLF   groupVII 25239495      98 0.30000000         7     G
#> 585   pop      OPL   groupVII 25239495      98 0.09615385        47     G
#> 586   pop      PAL   groupVII 25239495      98 0.40625000        19     G
#> 587   pop      SMR   groupVII 25239495      98 0.20454545        35     G
#> 588   pop      UPD   groupVII 25239495      98 0.34615385        17     G
#> 589   pop      ASP   groupVII 27599144      99 0.55000000         9     G
#> 590   pop      CLF   groupVII 27599144      99 0.41666667         7     G
#> 591   pop      OPL   groupVII 27599144      99 0.57142857        18     G
#> 592   pop      PAL   groupVII 27599144      99 0.00000000        28     G
#> 593   pop      SMR   groupVII 27599144      99 0.52500000        19     G
#> 594   pop      UPD   groupVII 27599144      99 0.50000000        11     G
#> 595   pop      ASP    groupIV 32321895     100 0.60000000         8     A
#> 596   pop      CLF    groupIV 32321895     100 0.37500000         5     A
#> 597   pop      OPL    groupIV 32321895     100 0.30000000        28     A
#> 598   pop      PAL    groupIV 32321895     100 0.31818182        15     A
#> 599   pop      SMR    groupIV 32321895     100 0.32500000        27     A
#> 600   pop      UPD    groupIV 32321895     100 0.60714286        11     A
#>     min.count minor         pi        pHWE pHWE_overall_BY pHWE_byfacet_BY
#> 1           9     A 0.52941176 0.502426985               1               1
#> 2           7     A 0.46666667 1.000000000               1               1
#> 3          23     A 0.50128205 0.354542084               1               1
#> 4           7     A 0.33689840 0.537263626               1               1
#> 5          26     A 0.49471459 1.000000000               1               1
#> 6          13     A 0.52000000 1.000000000               1               1
#> 7           1     T 0.10000000 1.000000000               1               1
#> 8           0     T 0.00000000 1.000000000               1               1
#> 9           3     T 0.12463768 1.000000000               1               1
#> 10          0     T 0.00000000 1.000000000               1               1
#> 11          2     T 0.08879493 1.000000000               1               1
#> 12          0     T 0.00000000 1.000000000               1               1
#> 13          2     T 0.17316017 1.000000000               1               1
#> 14          1     T 0.16666667 1.000000000               1               1
#> 15          4     T 0.15020408 1.000000000               1               1
#> 16          1     T 0.06250000 1.000000000               1               1
#> 17          3     T 0.12463768 1.000000000               1               1
#> 18          2     T 0.12873563 1.000000000               1               1
#> 19          3     G 0.24675325 1.000000000               1               1
#> 20          1     G 0.20000000 1.000000000               1               1
#> 21         17     G 0.50128205 0.001718035               1               1
#> 22          0     G 0.00000000 1.000000000               1               1
#> 23         15     G 0.44927536 0.642869188               1               1
#> 24         11     G 0.50769231 0.275057208               1               1
#> 25          5     A 0.49450549 1.000000000               1               1
#> 26          0     A 0.00000000 1.000000000               1               1
#> 27          8     A 0.44307692 1.000000000               1               1
#> 28          5     A 0.34420290 1.000000000               1               1
#> 29          7     A 0.37011494 1.000000000               1               1
#> 30          8     A 0.52287582 1.000000000               1               1
#> 31          5     C 0.39473684 1.000000000               1               1
#> 32          3     C 0.40909091 1.000000000               1               1
#> 33          3     C 0.11085973 1.000000000               1               1
#> 34         13     C 0.49798387 0.159026938               1               1
#> 35          3     C 0.13588850 1.000000000               1               1
#> 36         10     C 0.47619048 0.589931350               1               1
#> 37          1     G 0.09090909 1.000000000               1               1
#> 38          2     G 0.30303030 0.090909091               1               1
#> 39          6     G 0.25087108 1.000000000               1               1
#> 40          0     G 0.00000000 1.000000000               1               1
#> 41          8     G 0.31591173 1.000000000               1               1
#> 42          0     G 0.00000000 1.000000000               1               1
#> 43          5     A 0.36796537 1.000000000               1               1
#> 44          5     A 0.53030303 1.000000000               1               1
#> 45          8     A 0.28368794 0.501043917               1               1
#> 46          0     A 0.00000000 1.000000000               1               1
#> 47         14     A 0.46666667 1.000000000               1               1
#> 48          5     A 0.39473684 1.000000000               1               1
#> 49          5     A 0.36796537 0.438596491               1               1
#> 50          2     A 0.30303030 1.000000000               1               1
#> 51          5     A 0.18367347 1.000000000               1               1
#> 52          1     A 0.06250000 1.000000000               1               1
#> 53         10     A 0.38461538 0.277165729               1               1
#> 54          7     A 0.40923077 1.000000000               1               1
#> 55          2     G 0.17316017 1.000000000               1               1
#> 56          1     G 0.16666667 1.000000000               1               1
#> 57          1     G 0.04000000 1.000000000               1               1
#> 58          0     G 0.00000000 1.000000000               1               1
#> 59          2     G 0.08879493 1.000000000               1               1
#> 60          4     G 0.25396825 1.000000000               1               1
#> 61          9     C 0.52105263 1.000000000               1               1
#> 62          6     C 0.54545455 1.000000000               1               1
#> 63         25     C 0.49361208 1.000000000               1               1
#> 64          0     C 0.00000000 1.000000000               1               1
#> 65         19     C 0.51351351 0.662784266               1               1
#> 66         11     C 0.49470899 0.036392915               1               1
#> 67          5     G 0.39473684 1.000000000               1               1
#> 68          1     G 0.16666667 1.000000000               1               1
#> 69          7     G 0.32222222 1.000000000               1               1
#> 70         10     G 0.44354839 0.233931421               1               1
#> 71          6     G 0.29946524 1.000000000               1               1
#> 72          6     G 0.44210526 0.479876161               1               1
#> 73          2     C 0.20915033 1.000000000               1               1
#> 74          5     C 0.53030303 1.000000000               1               1
#> 75         12     C 0.38297872 0.596423748               1               1
#> 76          1     C 0.06666667 1.000000000               1               1
#> 77          8     C 0.32820513 1.000000000               1               1
#> 78          4     C 0.27076923 1.000000000               1               1
#> 79          4     T 0.33684211 1.000000000               1               1
#> 80          2     T 0.30303030 1.000000000               1               1
#> 81          4     T 0.15020408 0.001302649               1               1
#> 82          0     T 0.00000000 1.000000000               1               1
#> 83          8     T 0.30443975 1.000000000               1               1
#> 84         11     T 0.48045977 0.577811094               1               1
#> 85          5     G 0.36796537 1.000000000               1               1
#> 86          0     G 0.00000000 1.000000000               1               1
#> 87          8     G 0.32820513 1.000000000               1               1
#> 88         28     G 0.22580645 1.000000000               1               1
#> 89          9     G 0.34494774 0.532680899               1               1
#> 90         12     G 0.51692308 0.111589716               1               1
#> 91          5     T 0.39473684 0.479876161               1               1
#> 92          3     T 0.46666667 1.000000000               1               1
#> 93         13     T 0.45000000 0.606942076               1               1
#> 94         13     T 0.50804598 0.614845209               1               1
#> 95         13     T 0.45000000 0.606942076               1               1
#> 96          5     T 0.34420290 0.403726708               1               1
#> 97          3     C 0.26842105 1.000000000               1               1
#> 98          0     C 0.00000000 1.000000000               1               1
#> 99          7     C 0.24571429 1.000000000               1               1
#> 100         6     C 0.31451613 0.434433321               1               1
#> 101         7     C 0.28455285 1.000000000               1               1
#> 102         4     C 0.27076923 1.000000000               1               1
#> 103         6     A 0.44210526 0.479876161               1               1
#> 104         4     A 0.53333333 0.428571429               1               1
#> 105        11     A 0.33053809 1.000000000               1               1
#> 106        20     A 0.42328042 1.000000000               1               1
#> 107         7     A 0.26376812 0.413726602               1               1
#> 108        13     A 0.50804598 0.293882883               1               1
#> 109         0     G 0.00000000 1.000000000               1               1
#> 110         3     G 0.40909091 1.000000000               1               1
#> 111         3     G 0.14230769 1.000000000               1               1
#> 112         0     G 0.00000000 1.000000000               1               1
#> 113         3     G 0.12463768 1.000000000               1               1
#> 114         2     G 0.12873563 1.000000000               1               1
#> 115         1     A 0.11111111 1.000000000               1               1
#> 116         1     A 0.16666667 1.000000000               1               1
#> 117         1     A 0.04166667 1.000000000               1               1
#> 118         0     A 0.00000000 1.000000000               1               1
#> 119         3     A 0.12463768 1.000000000               1               1
#> 120         2     A 0.13756614 1.000000000               1               1
#> 121         3     A 0.29411765 1.000000000               1               1
#> 122         3     A 0.53571429 1.000000000               1               1
#> 123        19     A 0.51351351 0.187616641               1               1
#> 124         0     A 0.00000000 1.000000000               1               1
#> 125        13     A 0.45000000 0.057554297               1               1
#> 126         8     A 0.46376812 1.000000000               1               1
#> 127         8     T 0.52287582 1.000000000               1               1
#> 128         4     T 0.53333333 1.000000000               1               1
#> 129        13     T 0.43786295 1.000000000               1               1
#> 130         0     T 0.00000000 1.000000000               1               1
#> 131         6     T 0.26153846 1.000000000               1               1
#> 132        10     T 0.47619048 0.589931350               1               1
#> 133         1     A 0.10000000 1.000000000               1               1
#> 134         0     A 0.00000000 1.000000000               1               1
#> 135         2     A 0.07836735 1.000000000               1               1
#> 136         0     A 0.00000000 1.000000000               1               1
#> 137         4     A 0.18461538 1.000000000               1               1
#> 138         1     A 0.06666667 1.000000000               1               1
#> 139         1     G 0.10000000 1.000000000               1               1
#> 140         0     G 0.00000000 1.000000000               1               1
#> 141         5     G 0.17722474 0.190076030               1               1
#> 142         0     G 0.00000000 1.000000000               1               1
#> 143         2     G 0.10241821 1.000000000               1               1
#> 144         3     G 0.19841270 1.000000000               1               1
#> 145         6     A 0.44210526 0.479876161               1               1
#> 146         1     A 0.20000000 1.000000000               1               1
#> 147        12     A 0.39420290 0.269879313               1               1
#> 148         0     A 0.00000000 1.000000000               1               1
#> 149        19     A 0.51351351 0.662784266               1               1
#> 150        13     A 0.50804598 1.000000000               1               1
#> 151         6     C 0.47058824 0.457013575               1               1
#> 152         0     C 0.00000000 1.000000000               1               1
#> 153         6     C 0.25087108 1.000000000               1               1
#> 154        16     C 0.51494253 1.000000000               1               1
#> 155        13     C 0.42600423 0.119927642               1               1
#> 156         7     C 0.43115942 0.487414188               1               1
#> 157         1     C 0.10000000 1.000000000               1               1
#> 158         1     C 0.16666667 1.000000000               1               1
#> 159        10     C 0.33687943 0.539425154               1               1
#> 160         0     C 0.00000000 1.000000000               1               1
#> 161         6     C 0.26153846 1.000000000               1               1
#> 162         2     C 0.13756614 1.000000000               1               1
#> 163         2     C 0.23333333 1.000000000               1               1
#> 164         2     C 0.30303030 1.000000000               1               1
#> 165        25     C 0.46230441 1.000000000               1               1
#> 166         3     C 0.18620690 1.000000000               1               1
#> 167        21     C 0.51153846 1.000000000               1               1
#> 168         8     C 0.42328042 0.505314010               1               1
#> 169         1     C 0.11111111 1.000000000               1               1
#> 170         1     C 0.25000000 1.000000000               1               1
#> 171        11     C 0.39605110 0.578114700               1               1
#> 172         5     C 0.28735632 1.000000000               1               1
#> 173         7     C 0.37011494 1.000000000               1               1
#> 174         3     C 0.21230769 1.000000000               1               1
#> 175         2     C 0.17316017 1.000000000               1               1
#> 176         2     C 0.35555556 1.000000000               1               1
#> 177         3     C 0.13588850 1.000000000               1               1
#> 178         0     C 0.00000000 1.000000000               1               1
#> 179         5     C 0.20613108 1.000000000               1               1
#> 180         4     C 0.28985507 1.000000000               1               1
#> 181         6     G 0.41558442 0.504643963               1               1
#> 182         1     G 0.16666667 1.000000000               1               1
#> 183        16     G 0.44408163 1.000000000               1               1
#> 184         0     G 0.00000000 1.000000000               1               1
#> 185        14     G 0.47795164 1.000000000               1               1
#> 186        13     G 0.50804598 1.000000000               1               1
#> 187         3     C 0.26842105 1.000000000               1               1
#> 188         2     C 0.35555556 1.000000000               1               1
#> 189         6     C 0.23188406 1.000000000               1               1
#> 190        12     C 0.49655172 0.027854494               1               1
#> 191         7     C 0.29615385 1.000000000               1               1
#> 192         7     C 0.40923077 0.499130435               1               1
#> 193         6     T 0.41558442 0.504643963               1               1
#> 194         3     T 0.53571429 1.000000000               1               1
#> 195        14     T 0.47795164 0.647820100               1               1
#> 196         7     T 0.37011494 0.527969349               1               1
#> 197        15     T 0.45983087 0.178067424               1               1
#> 198         9     T 0.47076923 0.091990847               1               1
#> 199         7     G 0.50326797 0.176470588               1               1
#> 200         6     G 0.54545455 0.480519481               1               1
#> 201        14     G 0.48888889 0.645906079               1               1
#> 202         0     G 0.00000000 1.000000000               1               1
#> 203        13     G 0.42600423 0.311247720               1               1
#> 204         3     G 0.24675325 1.000000000               1               1
#> 205         2     C 0.18947368 1.000000000               1               1
#> 206         0     C 0.00000000 1.000000000               1               1
#> 207         3     C 0.13002114 1.000000000               1               1
#> 208         0     C 0.00000000 1.000000000               1               1
#> 209         3     C 0.14230769 1.000000000               1               1
#> 210         2     C 0.13756614 1.000000000               1               1
#> 211         3     G 0.26842105 1.000000000               1               1
#> 212         2     G 0.30303030 1.000000000               1               1
#> 213         9     G 0.31117021 1.000000000               1               1
#> 214         0     G 0.00000000 1.000000000               1               1
#> 215         8     G 0.30443975 1.000000000               1               1
#> 216         2     G 0.14769231 1.000000000               1               1
#> 217         4     C 0.33684211 0.306501548               1               1
#> 218         2     C 0.35555556 1.000000000               1               1
#> 219        13     C 0.40336879 0.135048156               1               1
#> 220         1     C 0.06250000 1.000000000               1               1
#> 221        16     C 0.47357294 1.000000000               1               1
#> 222         7     C 0.40923077 0.499130435               1               1
#> 223         9     T 0.52105263 0.245534651               1               1
#> 224         0     T 0.00000000 1.000000000               1               1
#> 225         9     T 0.32173913 0.544273310               1               1
#> 226        13     T 0.49798387 0.159026938               1               1
#> 227        16     T 0.50071124 0.355612865               1               1
#> 228         7     T 0.38888889 0.142028986               1               1
#> 229         3     A 0.29411765 1.000000000               1               1
#> 230         2     A 0.53333333 1.000000000               1               1
#> 231         3     A 0.18620690 1.000000000               1               1
#> 232         3     A 0.24675325 1.000000000               1               1
#> 233         3     A 0.16577540 1.000000000               1               1
#> 234         3     A 0.26842105 1.000000000               1               1
#> 235         9     C 0.52105263 1.000000000               1               1
#> 236         2     C 0.35555556 1.000000000               1               1
#> 237         7     C 0.29615385 1.000000000               1               1
#> 238         5     C 0.28735632 1.000000000               1               1
#> 239         6     C 0.27311522 1.000000000               1               1
#> 240         9     C 0.50649351 1.000000000               1               1
#> 241        13     A 0.50649351 0.554179567               1               1
#> 242         3     A 0.46666667 1.000000000               1               1
#> 243        10     A 0.39829303 0.263260024               1               1
#> 244        26     A 0.00000000 1.000000000               1               1
#> 245        11     A 0.45098039 0.608453838               1               1
#> 246        12     A 0.50526316 0.573231722               1               1
#> 247         4     A 0.40000000 1.000000000               1               1
#> 248         4     A 0.53333333 1.000000000               1               1
#> 249        21     A 0.50724638 1.000000000               1               1
#> 250         2     A 0.12096774 1.000000000               1               1
#> 251        25     A 0.50211416 0.185978072               1               1
#> 252        10     A 0.47619048 0.589931350               1               1
#> 253         2     G 0.17316017 1.000000000               1               1
#> 254         2     G 0.35555556 1.000000000               1               1
#> 255         8     G 0.31591173 1.000000000               1               1
#> 256         0     G 0.00000000 1.000000000               1               1
#> 257         5     G 0.19060284 1.000000000               1               1
#> 258         6     G 0.36923077 1.000000000               1               1
#> 259         3     A 0.29411765 0.176470588               1               1
#> 260         0     A 0.00000000 1.000000000               1               1
#> 261         3     A 0.11085973 1.000000000               1               1
#> 262        12     A 0.49655172 1.000000000               1               1
#> 263         5     A 0.21486643 0.234521576               1               1
#> 264         5     A 0.34420290 1.000000000               1               1
#> 265         0     T 0.00000000 1.000000000               1               1
#> 266         0     T 0.00000000 1.000000000               1               1
#> 267         0     T 0.00000000 1.000000000               1               1
#> 268         0     T 0.00000000 1.000000000               1               1
#> 269         0     T 0.00000000 1.000000000               1               1
#> 270         2     T 0.12873563 1.000000000               1               1
#> 271         7     A 0.52500000 0.048951049               1               1
#> 272         2     A 0.35555556 1.000000000               1               1
#> 273        15     A 0.50000000 1.000000000               1               1
#> 274         0     A 0.00000000 1.000000000               1               1
#> 275         7     A 0.33689840 0.095560724               1               1
#> 276         9     A 0.47076923 0.578215103               1               1
#> 277         5     G 0.36796537 1.000000000               1               1
#> 278         0     G 0.00000000 1.000000000               1               1
#> 279        22     G 0.51014493 0.420947138               1               1
#> 280         3     G 0.18620690 1.000000000               1               1
#> 281        11     G 0.37198068 0.560923985               1               1
#> 282         5     G 0.28735632 0.325670498               1               1
#> 283         5     C 0.42483660 1.000000000               1               1
#> 284         6     C 0.54545455 1.000000000               1               1
#> 285        23     C 0.50128205 1.000000000               1               1
#> 286        11     C 0.49470899 0.570404272               1               1
#> 287        16     C 0.50793651 1.000000000               1               1
#> 288        10     C 0.49230769 0.565809665               1               1
#> 289         0     C 0.00000000 1.000000000               1               1
#> 290         2     C 0.30303030 1.000000000               1               1
#> 291         1     C 0.03703704 1.000000000               1               1
#> 292         0     C 0.00000000 1.000000000               1               1
#> 293         0     C 0.00000000 1.000000000               1               1
#> 294         2     C 0.14769231 1.000000000               1               1
#> 295         3     T 0.26842105 1.000000000               1               1
#> 296         1     T 0.20000000 1.000000000               1               1
#> 297         5     T 0.19806763 1.000000000               1               1
#> 298         1     T 0.06250000 1.000000000               1               1
#> 299         2     T 0.10241821 1.000000000               1               1
#> 300         2     T 0.13756614 1.000000000               1               1
#> 301         7     C 0.50326797 1.000000000               1               1
#> 302         0     C 0.00000000 1.000000000               1               1
#> 303        14     C 0.46666667 1.000000000               1               1
#> 304         4     C 0.31168831 1.000000000               1               1
#> 305         9     C 0.45238095 1.000000000               1               1
#> 306        11     C 0.49470899 1.000000000               1               1
#> 307         8     G 0.48484848 0.190918473               1               1
#> 308         4     G 0.48484848 1.000000000               1               1
#> 309        16     G 0.45390071 1.000000000               1               1
#> 310        11     G 0.51811594 0.238505471               1               1
#> 311         7     G 0.25443262 0.398317665               1               1
#> 312        10     G 0.47619048 1.000000000               1               1
#> 313         4     C 0.33684211 0.306501548               1               1
#> 314         3     C 0.53571429 1.000000000               1               1
#> 315        13     C 0.43786295 0.610776788               1               1
#> 316         0     C 0.00000000 1.000000000               1               1
#> 317         8     C 0.29371981 0.518955161               1               1
#> 318        10     C 0.50724638 1.000000000               1               1
#> 319         0     T 0.00000000 1.000000000               1               1
#> 320         1     T 0.16666667 1.000000000               1               1
#> 321         1     T 0.03846154 1.000000000               1               1
#> 322         0     T 0.00000000 1.000000000               1               1
#> 323         1     T 0.04347826 1.000000000               1               1
#> 324         0     T 0.00000000 1.000000000               1               1
#> 325         8     A 0.53333333 0.477855478               1               1
#> 326         3     A 0.40909091 1.000000000               1               1
#> 327        14     A 0.45528455 0.624198278               1               1
#> 328        14     A 0.51494253 0.325979116               1               1
#> 329         7     A 0.29615385 0.467775468               1               1
#> 330        12     A 0.51692308 1.000000000               1               1
#> 331         6     G 0.47058824 0.457013575               1               1
#> 332         1     G 0.16666667 1.000000000               1               1
#> 333         8     G 0.32820513 0.548415548               1               1
#> 334        11     G 0.46572581 0.094904161               1               1
#> 335         7     G 0.30867710 1.000000000               1               1
#> 336         7     G 0.37011494 0.527969349               1               1
#> 337         1     A 0.12500000 1.000000000               1               1
#> 338         0     A 0.00000000 1.000000000               1               1
#> 339         4     A 0.18461538 1.000000000               1               1
#> 340         0     A 0.00000000 1.000000000               1               1
#> 341         1     A 0.05555556 1.000000000               1               1
#> 342         0     A 0.00000000 1.000000000               1               1
#> 343         8     T 0.52287582 1.000000000               1               1
#> 344         0     T 0.00000000 1.000000000               1               1
#> 345        13     T 0.39265306 0.144646715               1               1
#> 346        11     T 0.50769231 0.573563064               1               1
#> 347         5     T 0.24603175 1.000000000               1               1
#> 348         6     T 0.41558442 1.000000000               1               1
#> 349        10     A 0.51948052 0.542748273               1               1
#> 350         4     A 0.48484848 1.000000000               1               1
#> 351        19     A 0.48847518 0.216808603               1               1
#> 352         2     A 0.12096774 1.000000000               1               1
#> 353        19     A 0.51153846 0.367229217               1               1
#> 354        10     A 0.49230769 1.000000000               1               1
#> 355         6     C 0.52747253 1.000000000               1               1
#> 356         9     C 0.40909091 1.000000000               1               1
#> 357        16     C 0.46376812 0.182055237               1               1
#> 358         0     C 0.00000000 1.000000000               1               1
#> 359        14     C 0.45528455 1.000000000               1               1
#> 360         9     C 0.45238095 0.220869565               1               1
#> 361         0     C 0.00000000 1.000000000               1               1
#> 362         1     C 0.33333333 1.000000000               1               1
#> 363         1     C 0.05882353 1.000000000               1               1
#> 364         0     C 0.00000000 1.000000000               1               1
#> 365         2     C 0.09743590 1.000000000               1               1
#> 366         2     C 0.15942029 1.000000000               1               1
#> 367         2     C 0.18947368 1.000000000               1               1
#> 368         0     C 0.00000000 1.000000000               1               1
#> 369        11     C 0.34012066 0.549761518               1               1
#> 370         0     C 0.00000000 1.000000000               1               1
#> 371         3     C 0.11968085 1.000000000               1               1
#> 372         3     C 0.19841270 1.000000000               1               1
#> 373         6     C 0.50000000 1.000000000               1               1
#> 374         4     C 0.53333333 1.000000000               1               1
#> 375         8     C 0.35555556 0.528850238               1               1
#> 376        16     C 0.33684211 1.000000000               1               1
#> 377        10     C 0.44354839 1.000000000               1               1
#> 378         7     C 0.40923077 1.000000000               1               1
#> 379         7     A 0.45454545 1.000000000               1               1
#> 380         2     A 0.30303030 0.090909091               1               1
#> 381        13     A 0.43786295 0.115401790               1               1
#> 382        16     A 0.50793651 0.620944945               1               1
#> 383        16     A 0.47357294 1.000000000               1               1
#> 384        11     A 0.48045977 0.577811094               1               1
#> 385         4     T 0.31168831 1.000000000               1               1
#> 386         5     T 0.53030303 1.000000000               1               1
#> 387        24     T 0.47795164 0.647820100               1               1
#> 388         9     T 0.40106952 1.000000000               1               1
#> 389        15     T 0.45983087 0.654389537               1               1
#> 390        11     T 0.49470899 0.570404272               1               1
#> 391         9     C 0.52941176 1.000000000               1               1
#> 392         5     C 0.53030303 1.000000000               1               1
#> 393        14     C 0.46666667 0.327434220               1               1
#> 394        21     C 0.46572581 1.000000000               1               1
#> 395        13     C 0.48663102 0.624445842               1               1
#> 396        11     C 0.50769231 0.573563064               1               1
#> 397         0     A 0.00000000 1.000000000               1               1
#> 398         0     A 0.00000000 1.000000000               1               1
#> 399         0     A 0.00000000 1.000000000               1               1
#> 400         8     A 0.40459770 1.000000000               1               1
#> 401         0     A 0.00000000 1.000000000               1               1
#> 402         1     A 0.07692308 1.000000000               1               1
#> 403         5     T 0.39473684 0.479876161               1               1
#> 404         3     T 0.40909091 1.000000000               1               1
#> 405        11     T 0.36081560 1.000000000               1               1
#> 406        19     T 0.40923077 1.000000000               1               1
#> 407        10     T 0.41269841 1.000000000               1               1
#> 408        10     T 0.47619048 1.000000000               1               1
#> 409         0     G 0.00000000 1.000000000               1               1
#> 410         1     G 0.20000000 1.000000000               1               1
#> 411        15     G 0.42857143 0.630928644               1               1
#> 412         0     G 0.00000000 1.000000000               1               1
#> 413         4     G 0.17653891 1.000000000               1               1
#> 414         0     G 0.00000000 1.000000000               1               1
#> 415         7     C 0.47894737 0.219814241               1               1
#> 416         6     C 0.53333333 1.000000000               1               1
#> 417        11     C 0.39605110 0.261700725               1               1
#> 418         0     C 0.00000000 1.000000000               1               1
#> 419        11     C 0.38372093 0.269436547               1               1
#> 420         7     C 0.43115942 0.487414188               1               1
#> 421         1     T 0.09090909 1.000000000               1               1
#> 422         0     T 0.00000000 1.000000000               1               1
#> 423         1     T 0.03703704 1.000000000               1               1
#> 424         6     T 0.31451613 1.000000000               1               1
#> 425         6     T 0.22340426 1.000000000               1               1
#> 426         2     T 0.12873563 1.000000000               1               1
#> 427        13     G 0.47894737 0.479876161               1               1
#> 428         7     G 0.53030303 1.000000000               1               1
#> 429        15     G 0.43882979 0.356657970               1               1
#> 430        19     G 0.48045977 1.000000000               1               1
#> 431        10     G 0.35940803 0.538114842               1               1
#> 432        19     G 0.40923077 0.499130435               1               1
#> 433         6     G 0.44210526 0.479876161               1               1
#> 434         3     G 0.40909091 1.000000000               1               1
#> 435         6     G 0.23188406 1.000000000               1               1
#> 436         9     G 0.43448276 1.000000000               1               1
#> 437         5     G 0.21486643 1.000000000               1               1
#> 438         7     G 0.38888889 1.000000000               1               1
#> 439         3     T 0.26842105 1.000000000               1               1
#> 440         1     T 0.20000000 1.000000000               1               1
#> 441        12     T 0.44381223 0.109068706               1               1
#> 442         0     T 0.00000000 1.000000000               1               1
#> 443         8     T 0.37076649 0.520125242               1               1
#> 444         3     T 0.21230769 1.000000000               1               1
#> 445        10     C 0.52287582 0.539284245               1               1
#> 446         5     C 0.55555556 1.000000000               1               1
#> 447        20     C 0.50241546 0.405280123               1               1
#> 448        10     C 0.44354839 0.233931421               1               1
#> 449        19     C 0.51153846 1.000000000               1               1
#> 450        21     C 0.43448276 0.237181409               1               1
#> 451         5     T 0.53030303 1.000000000               1               1
#> 452         1     T 0.33333333 1.000000000               1               1
#> 453        11     T 0.46572581 0.025571085               1               1
#> 454        16     T 0.00000000 1.000000000               1               1
#> 455         3     T 0.22826087 1.000000000               1               1
#> 456        12     T 0.50793651 0.620944945               1               1
#> 457         6     T 0.47058824 0.457013575               1               1
#> 458         1     T 0.16666667 1.000000000               1               1
#> 459         8     T 0.29371981 1.000000000               1               1
#> 460        11     T 0.48045977 0.085257371               1               1
#> 461         7     T 0.30867710 1.000000000               1               1
#> 462         6     T 0.34920635 1.000000000               1               1
#> 463         0     T 0.00000000 1.000000000               1               1
#> 464         3     T 0.40909091 0.272727273               1               1
#> 465        13     T 0.42600423 1.000000000               1               1
#> 466         0     T 0.00000000 1.000000000               1               1
#> 467         5     T 0.20613108 1.000000000               1               1
#> 468         0     T 0.00000000 1.000000000               1               1
#> 469         4     C 0.31168831 1.000000000               1               1
#> 470         0     C 0.00000000 1.000000000               1               1
#> 471         6     C 0.20814480 1.000000000               1               1
#> 472         8     C 0.37076649 0.520125242               1               1
#> 473         4     C 0.16913319 0.137833239               1               1
#> 474         3     C 0.22826087 1.000000000               1               1
#> 475         5     C 0.36796537 0.438596491               1               1
#> 476         2     C 0.30303030 1.000000000               1               1
#> 477        14     C 0.44397463 0.360620643               1               1
#> 478        15     C 0.51724138 0.141779781               1               1
#> 479        12     C 0.39420290 0.269879313               1               1
#> 480         7     C 0.40923077 0.499130435               1               1
#> 481         6     T 0.47058824 0.457013575               1               1
#> 482         6     T 0.53333333 0.428571429               1               1
#> 483        17     T 0.49361208 0.367093836               1               1
#> 484        14     T 0.50806452 0.608633035               1               1
#> 485        14     T 0.44397463 0.360620643               1               1
#> 486        11     T 0.49470899 0.319806763               1               1
#> 487         2     C 0.18947368 1.000000000               1               1
#> 488         0     C 0.00000000 1.000000000               1               1
#> 489         5     C 0.17722474 1.000000000               1               1
#> 490         0     C 0.00000000 1.000000000               1               1
#> 491         2     C 0.10241821 1.000000000               1               1
#> 492         5     C 0.32307692 1.000000000               1               1
#> 493         4     A 0.31168831 1.000000000               1               1
#> 494         5     A 0.53030303 0.393939394               1               1
#> 495         0     A 0.00000000 1.000000000               1               1
#> 496         0     A 0.00000000 1.000000000               1               1
#> 497         1     A 0.04166667 1.000000000               1               1
#> 498         1     A 0.07142857 1.000000000               1               1
#> 499         6     T 0.44210526 1.000000000               1               1
#> 500         4     T 0.48484848 1.000000000               1               1
#> 501        14     T 0.46666667 1.000000000               1               1
#> 502        30     T 0.00000000 1.000000000               1               1
#> 503        10     T 0.39829303 1.000000000               1               1
#> 504        12     T 0.51692308 0.601992193               1               1
#> 505        13     A 0.47894737 1.000000000               1               1
#> 506         6     A 0.53333333 0.428571429               1               1
#> 507        13     A 0.43786295 0.115401790               1               1
#> 508        25     A 0.19841270 1.000000000               1               1
#> 509        11     A 0.39605110 1.000000000               1               1
#> 510        13     A 0.50804598 0.614845209               1               1
#> 511         6     G 0.50000000 0.440559441               1               1
#> 512         1     G 0.16666667 1.000000000               1               1
#> 513         9     G 0.35769231 0.528255528               1               1
#> 514         9     G 0.50649351 0.065900044               1               1
#> 515         7     G 0.30867710 1.000000000               1               1
#> 516         7     G 0.40923077 0.499130435               1               1
#> 517        11     C 0.52380952 1.000000000               1               1
#> 518         1     C 0.16666667 1.000000000               1               1
#> 519        12     C 0.37224490 0.586422587               1               1
#> 520         0     C 0.00000000 1.000000000               1               1
#> 521        12     C 0.40591966 0.143248314               1               1
#> 522        12     C 0.50793651 1.000000000               1               1
#> 523         8     G 0.53333333 0.477855478               1               1
#> 524         3     G 0.40909091 1.000000000               1               1
#> 525        13     G 0.45000000 1.000000000               1               1
#> 526        28     G 0.00000000 1.000000000               1               1
#> 527        13     G 0.49798387 0.159026938               1               1
#> 528         8     G 0.46376812 0.216045228               1               1
#> 529         6     G 0.47058824 0.457013575               1               1
#> 530         2     G 0.53333333 0.200000000               1               1
#> 531        13     G 0.47460317 1.000000000               1               1
#> 532        12     G 0.52173913 0.026458533               1               1
#> 533        13     G 0.47460317 1.000000000               1               1
#> 534        14     G 0.50724638 0.289271773               1               1
#> 535         3     T 0.24675325 1.000000000               1               1
#> 536         0     T 0.00000000 1.000000000               1               1
#> 537         0     T 0.00000000 1.000000000               1               1
#> 538         0     T 0.00000000 1.000000000               1               1
#> 539         1     T 0.04347826 1.000000000               1               1
#> 540         2     T 0.14769231 1.000000000               1               1
#> 541         3     T 0.29411765 1.000000000               1               1
#> 542         3     T 0.53571429 0.428571429               1               1
#> 543        21     T 0.51219512 0.403464236               1               1
#> 544         3     T 0.17540323 1.000000000               1               1
#> 545        22     T 0.50769231 1.000000000               1               1
#> 546        10     T 0.47619048 1.000000000               1               1
#> 547         7     A 0.53846154 1.000000000               1               1
#> 548         2     A 0.35555556 1.000000000               1               1
#> 549        10     A 0.47619048 0.091990847               1               1
#> 550        14     A 0.48484848 1.000000000               1               1
#> 551         4     A 0.31168831 1.000000000               1               1
#> 552         4     A 0.33684211 1.000000000               1               1
#> 553         3     C 0.24675325 1.000000000               1               1
#> 554         3     C 0.40909091 1.000000000               1               1
#> 555         3     C 0.11085973 1.000000000               1               1
#> 556         4     C 0.25396825 1.000000000               1               1
#> 557         7     C 0.30867710 1.000000000               1               1
#> 558         7     C 0.38888889 0.513043478               1               1
#> 559         0     G 0.00000000 1.000000000               1               1
#> 560         0     G 0.00000000 1.000000000               1               1
#> 561         3     G 0.13002114 1.000000000               1               1
#> 562         9     G 0.47076923 1.000000000               1               1
#> 563         1     G 0.04545455 1.000000000               1               1
#> 564         3     G 0.19841270 1.000000000               1               1
#> 565         4     T 0.33684211 1.000000000               1               1
#> 566         4     T 0.48484848 1.000000000               1               1
#> 567         8     T 0.30443975 0.107022027               1               1
#> 568         0     T 0.00000000 1.000000000               1               1
#> 569         6     T 0.27311522 1.000000000               1               1
#> 570         4     T 0.28985507 1.000000000               1               1
#> 571         1     A 0.12500000 1.000000000               1               1
#> 572         1     A 0.16666667 1.000000000               1               1
#> 573         8     A 0.29371981 1.000000000               1               1
#> 574         0     A 0.00000000 1.000000000               1               1
#> 575         4     A 0.16231884 1.000000000               1               1
#> 576         2     A 0.14769231 1.000000000               1               1
#> 577         0     A 0.00000000 1.000000000               1               1
#> 578         1     A 0.16666667 1.000000000               1               1
#> 579         2     A 0.09291521 1.000000000               1               1
#> 580         0     A 0.00000000 1.000000000               1               1
#> 581         2     A 0.09743590 1.000000000               1               1
#> 582         1     A 0.06666667 1.000000000               1               1
#> 583         4     T 0.36601307 0.341176471               1               1
#> 584         3     T 0.46666667 0.333333333               1               1
#> 585         5     T 0.17722474 1.000000000               1               1
#> 586        13     T 0.49798387 1.000000000               1               1
#> 587         9     T 0.33298097 0.538114842               1               1
#> 588         9     T 0.47076923 1.000000000               1               1
#> 589        11     A 0.52105263 1.000000000               1               1
#> 590         5     A 0.53030303 0.151515152               1               1
#> 591        24     A 0.50174216 0.187539953               1               1
#> 592         0     A 0.00000000 1.000000000               1               1
#> 593        21     A 0.51153846 0.669858722               1               1
#> 594        11     A 0.52380952 1.000000000               1               1
#> 595        12     T 0.50526316 0.573231722               1               1
#> 596         3     T 0.53571429 1.000000000               1               1
#> 597        12     T 0.43076923 1.000000000               1               1
#> 598         7     T 0.45454545 1.000000000               1               1
#> 599        13     T 0.45000000 0.057554297               1               1
#> 600        17     T 0.49470899 1.000000000               1               1
#>             ho         he           fis pa_corrected  g
#> 1   0.33333333 0.50000000  3.846154e-01 0.000000e+00 10
#> 2   0.60000000 0.42000000 -3.333333e-01 0.000000e+00 10
#> 3   0.65000000 0.48875000 -3.068783e-01 0.000000e+00 10
#> 4   0.29411765 0.32698962  1.304348e-01 0.000000e+00 10
#> 5   0.45454545 0.48347107  8.296943e-02 0.000000e+00 10
#> 6   0.53846154 0.50000000 -3.703704e-02 0.000000e+00 10
#> 7   0.10000000 0.09500000  5.551115e-16 1.240114e-01 12
#> 8   0.00000000 0.00000000           NaN 0.000000e+00 12
#> 9   0.13043478 0.12192816 -4.761905e-02 1.270509e-01 12
#> 10  0.00000000 0.00000000           NaN 0.000000e+00 12
#> 11  0.09090909 0.08677687 -2.439024e-02 7.500689e-02 12
#> 12  0.00000000 0.00000000           NaN 0.000000e+00 12
#> 13  0.18181818 0.16528925 -5.263158e-02 0.000000e+00 12
#> 14  0.16666667 0.15277777  3.330669e-16 5.410697e-03 12
#> 15  0.16000000 0.14720000 -6.666667e-02 0.000000e+00 12
#> 16  0.06250000 0.06054688  0.000000e+00 0.000000e+00 12
#> 17  0.13043478 0.12192816 -4.761905e-02 0.000000e+00 12
#> 18  0.13333333 0.12444445 -3.703704e-02 0.000000e+00 12
#> 19  0.27272727 0.23553720 -1.111111e-01 0.000000e+00 10
#> 20  0.20000000 0.18000000 -2.220446e-16 1.186024e-09 10
#> 21  0.85000000 0.48875000 -7.272727e-01 0.000000e+00 10
#> 22  0.00000000 0.00000000           NaN 0.000000e+00 10
#> 23  0.39130435 0.43950851  1.315789e-01 0.000000e+00 10
#> 24  0.69230769 0.48816568 -3.846154e-01 0.000000e+00 10
#> 25  0.42857143 0.45918368  1.428571e-01 3.038486e-05  6
#> 26  0.00000000 0.00000000           NaN 0.000000e+00  6
#> 27  0.46153846 0.42603550 -4.347826e-02 9.969699e-06  6
#> 28  0.41666667 0.32986111 -2.222222e-01 3.463249e-06  6
#> 29  0.33333333 0.35777777  1.025641e-01 4.268787e-06  6
#> 30  0.44444444 0.49382716  1.578947e-01 7.642121e-05  6
#> 31  0.50000000 0.37500000 -2.857143e-01 0.000000e+00 12
#> 32  0.50000000 0.37500000 -2.500000e-01 7.776736e-11 12
#> 33  0.11538462 0.10872781 -4.166667e-02 0.000000e+00 12
#> 34  0.68750000 0.48242188 -3.983051e-01 0.000000e+00 12
#> 35  0.14285714 0.13265306 -5.263158e-02 0.000000e+00 12
#> 36  0.57142857 0.45918368 -2.093023e-01 0.000000e+00 12
#> 37  0.09090909 0.08677687 -4.440892e-16 0.000000e+00 12
#> 38  0.00000000 0.27777778  1.000000e+00 2.551341e-03 12
#> 39  0.28571429 0.24489796 -1.428571e-01 0.000000e+00 12
#> 40  0.00000000 0.00000000           NaN 0.000000e+00 12
#> 41  0.38095238 0.30839002 -2.121212e-01 0.000000e+00 12
#> 42  0.00000000 0.00000000           NaN 0.000000e+00 12
#> 43  0.45454545 0.35123967 -2.500000e-01 0.000000e+00 12
#> 44  0.50000000 0.48611111  6.250000e-02 4.791429e-09 12
#> 45  0.25000000 0.27777778  1.210191e-01 0.000000e+00 12
#> 46  0.00000000 0.00000000           NaN 0.000000e+00 12
#> 47  0.50000000 0.45500000 -7.344633e-02 0.000000e+00 12
#> 48  0.50000000 0.37500000 -2.857143e-01 0.000000e+00 12
#> 49  0.27272727 0.35123967  2.682927e-01 0.000000e+00 12
#> 50  0.33333333 0.27777778 -1.111111e-01 1.144468e-07 12
#> 51  0.20000000 0.18000000 -9.090909e-02 0.000000e+00 12
#> 52  0.06250000 0.06054688  0.000000e+00 0.000000e+00 12
#> 53  0.50000000 0.37500000 -3.103448e-01 0.000000e+00 12
#> 54  0.38461538 0.39349112  6.250000e-02 0.000000e+00 12
#> 55  0.18181818 0.16528925 -5.263158e-02 0.000000e+00 12
#> 56  0.16666667 0.15277777  3.330669e-16 6.900049e-03 12
#> 57  0.04000000 0.03920000  5.551115e-16 0.000000e+00 12
#> 58  0.00000000 0.00000000           NaN 0.000000e+00 12
#> 59  0.09090909 0.08677687 -2.439024e-02 0.000000e+00 12
#> 60  0.28571429 0.24489796 -1.304348e-01 0.000000e+00 12
#> 61  0.50000000 0.49500000  4.255319e-02 0.000000e+00 12
#> 62  0.66666667 0.50000000 -2.500000e-01 0.000000e+00 12
#> 63  0.52380952 0.48185941 -6.280193e-02 0.000000e+00 12
#> 64  0.00000000 0.00000000           NaN 0.000000e+00 12
#> 65  0.57894737 0.50000000 -1.314286e-01 0.000000e+00 12
#> 66  0.78571429 0.47704082 -6.250000e-01 0.000000e+00 12
#> 67  0.50000000 0.37500000 -2.857143e-01 0.000000e+00 12
#> 68  0.16666667 0.15277777  3.330669e-16 1.718854e-11 12
#> 69  0.38888889 0.31327160 -2.142857e-01 0.000000e+00 12
#> 70  0.62500000 0.42968750 -4.285714e-01 0.000000e+00 12
#> 71  0.35294118 0.29065744 -1.851852e-01 0.000000e+00 12
#> 72  0.60000000 0.42000000 -3.846154e-01 0.000000e+00 12
#> 73  0.22222222 0.19753086 -6.666667e-02 0.000000e+00 12
#> 74  0.50000000 0.48611111  6.250000e-02 2.859846e-06 12
#> 75  0.33333333 0.37500000  1.320755e-01 0.000000e+00 12
#> 76  0.06666667 0.06444444 -2.220446e-16 0.000000e+00 12
#> 77  0.30000000 0.32000000  8.800000e-02 0.000000e+00 12
#> 78  0.30769231 0.26035502 -1.428571e-01 0.000000e+00 12
#> 79  0.40000000 0.32000000 -2.000000e-01 0.000000e+00 12
#> 80  0.33333333 0.27777778 -1.111111e-01 1.601039e-07 12
#> 81  0.00000000 0.14720000  1.000000e+00 0.000000e+00 12
#> 82  0.00000000 0.00000000           NaN 0.000000e+00 12
#> 83  0.36363636 0.29752066 -2.000000e-01 0.000000e+00 12
#> 84  0.60000000 0.46444445 -2.600000e-01 0.000000e+00 12
#> 85  0.45454545 0.35123967 -2.500000e-01 0.000000e+00 12
#> 86  0.00000000 0.00000000           NaN 0.000000e+00 12
#> 87  0.30000000 0.32000000  8.800000e-02 0.000000e+00 12
#> 88  0.25000000 0.21875000 -1.111111e-01 1.169278e-10 12
#> 89  0.42857143 0.33673469 -2.500000e-01 0.000000e+00 12
#> 90  0.76923077 0.49704142 -5.189873e-01 0.000000e+00 12
#> 91  0.30000000 0.37500000  2.500000e-01 0.000000e+00 10
#> 92  0.60000000 0.42000000 -3.333333e-01 4.908538e-11 10
#> 93  0.55000000 0.43875000 -2.294118e-01 0.000000e+00 10
#> 94  0.60000000 0.49111111 -1.886792e-01 0.000000e+00 10
#> 95  0.55000000 0.43875000 -2.294118e-01 0.000000e+00 10
#> 96  0.25000000 0.32986111  2.826087e-01 0.000000e+00 10
#> 97  0.30000000 0.25500000 -1.250000e-01 2.596399e-05 12
#> 98  0.00000000 0.00000000           NaN 0.000000e+00 12
#> 99  0.28000000 0.24080000 -1.428571e-01 9.274657e-06 12
#> 100 0.25000000 0.30468750  2.105263e-01 3.001801e-05 12
#> 101 0.33333333 0.27777778 -1.764706e-01 1.643371e-05 12
#> 102 0.30769231 0.26035502 -1.428571e-01 1.869129e-05 12
#> 103 0.60000000 0.42000000 -3.846154e-01 0.000000e+00 10
#> 104 0.80000000 0.48000000 -6.000000e-01 0.000000e+00 10
#> 105 0.33333333 0.32441701 -8.620690e-03 0.000000e+00 10
#> 106 0.42857143 0.40816327 -1.298701e-02 0.000000e+00 10
#> 107 0.21739130 0.25803402  1.791045e-01 0.000000e+00 10
#> 108 0.33333333 0.49111111  3.518519e-01 0.000000e+00 10
#> 109 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 110 0.50000000 0.37500000 -2.500000e-01 4.597364e-02 12
#> 111 0.15000000 0.13875000 -5.555556e-02 0.000000e+00 12
#> 112 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 113 0.13043478 0.12192816 -4.761905e-02 0.000000e+00 12
#> 114 0.13333333 0.12444445 -3.703704e-02 0.000000e+00 12
#> 115 0.11111111 0.10493828  2.220446e-16 0.000000e+00 12
#> 116 0.16666667 0.15277777  3.330669e-16 3.128594e-02 12
#> 117 0.04166667 0.04079860  7.771561e-16 0.000000e+00 12
#> 118 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 119 0.13043478 0.12192816 -4.761905e-02 0.000000e+00 12
#> 120 0.14285714 0.13265306 -4.000000e-02 0.000000e+00 12
#> 121 0.33333333 0.27777778 -1.428571e-01 0.000000e+00  8
#> 122 0.75000000 0.46875000 -5.000000e-01 1.148147e-07  8
#> 123 0.68421053 0.50000000 -3.448276e-01 0.000000e+00  8
#> 124 0.00000000 0.00000000           NaN 0.000000e+00  8
#> 125 0.65000000 0.43875000 -4.615385e-01 0.000000e+00  8
#> 126 0.50000000 0.44444444 -8.196721e-02 0.000000e+00  8
#> 127 0.44444444 0.49382716  1.578947e-01 0.000000e+00 10
#> 128 0.40000000 0.48000000  2.727273e-01 1.604630e-10 10
#> 129 0.42857143 0.42743764  2.173913e-02 0.000000e+00 10
#> 130 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 131 0.30000000 0.25500000 -1.515152e-01 0.000000e+00 10
#> 132 0.57142857 0.45918368 -2.093023e-01 0.000000e+00 10
#> 133 0.10000000 0.09500000  5.551115e-16 4.628571e-02 12
#> 134 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 135 0.08000000 0.07680000 -2.127660e-02 2.291242e-02 12
#> 136 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 137 0.20000000 0.18000000 -8.571429e-02 1.068735e-01 12
#> 138 0.06666667 0.06444444 -2.220446e-16 2.057143e-02 12
#> 139 0.10000000 0.09500000  5.551115e-16 1.200480e-02 12
#> 140 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 141 0.11538462 0.17381657  3.534483e-01 2.360736e-02 12
#> 142 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 143 0.10526316 0.09972299 -2.857143e-02 9.308339e-03 12
#> 144 0.21428571 0.19132654 -8.333333e-02 3.881553e-02 12
#> 145 0.60000000 0.42000000 -3.846154e-01 0.000000e+00 10
#> 146 0.20000000 0.18000000 -2.220446e-16 2.204522e-11 10
#> 147 0.52173913 0.38563327 -3.333333e-01 0.000000e+00 10
#> 148 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 149 0.57894737 0.50000000 -1.314286e-01 0.000000e+00 10
#> 150 0.46666667 0.49111111  8.411215e-02 0.000000e+00 10
#> 151 0.66666667 0.44444444 -4.545455e-01 1.823330e-12 12
#> 152 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 153 0.28571429 0.24489796 -1.428571e-01 7.695485e-16 12
#> 154 0.53333333 0.49777778 -3.703704e-02 9.335927e-11 12
#> 155 0.59090909 0.41632232 -4.000000e-01 1.458149e-14 12
#> 156 0.58333333 0.41319445 -3.750000e-01 4.282562e-14 12
#> 157 0.10000000 0.09500000  5.551115e-16 0.000000e+00 12
#> 158 0.16666667 0.15277777  3.330669e-16 4.843631e-04 12
#> 159 0.41666667 0.32986111 -2.432432e-01 0.000000e+00 12
#> 160 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 161 0.30000000 0.25500000 -1.515152e-01 0.000000e+00 12
#> 162 0.14285714 0.13265306 -4.000000e-02 0.000000e+00 12
#> 163 0.25000000 0.21875000 -7.692308e-02 0.000000e+00 12
#> 164 0.33333333 0.27777778 -1.111111e-01 1.802004e-18 12
#> 165 0.47368421 0.45013850 -2.531646e-02 0.000000e+00 12
#> 166 0.20000000 0.18000000 -7.692308e-02 0.000000e+00 12
#> 167 0.55000000 0.49875000 -7.731959e-02 0.000000e+00 12
#> 168 0.57142857 0.40816327 -3.684211e-01 0.000000e+00 12
#> 169 0.11111111 0.10493828  2.220446e-16 0.000000e+00  8
#> 170 0.25000000 0.21875000  0.000000e+00 1.804027e-04  8
#> 171 0.33333333 0.38662131  1.616766e-01 0.000000e+00  8
#> 172 0.33333333 0.27777778 -1.666667e-01 0.000000e+00  8
#> 173 0.33333333 0.35777777  1.025641e-01 0.000000e+00  8
#> 174 0.23076923 0.20414202 -9.090909e-02 0.000000e+00  8
#> 175 0.18181818 0.16528925 -5.263158e-02 0.000000e+00 10
#> 176 0.40000000 0.32000000 -1.428571e-01 2.979528e-03 10
#> 177 0.14285714 0.13265306 -5.263158e-02 0.000000e+00 10
#> 178 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 179 0.22727273 0.20144628 -1.052632e-01 0.000000e+00 10
#> 180 0.33333333 0.27777778 -1.578947e-01 0.000000e+00 10
#> 181 0.54545455 0.39669421 -3.333333e-01 0.000000e+00 12
#> 182 0.16666667 0.15277777  3.330669e-16 9.084141e-13 12
#> 183 0.48000000 0.43520000 -8.270677e-02 0.000000e+00 12
#> 184 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 185 0.52631579 0.46537396 -1.042945e-01 0.000000e+00 12
#> 186 0.46666667 0.49111111  8.411215e-02 0.000000e+00 12
#> 187 0.30000000 0.25500000 -1.250000e-01 0.000000e+00 10
#> 188 0.40000000 0.32000000 -1.428571e-01 6.054110e-08 10
#> 189 0.26086957 0.22684310 -1.282051e-01 0.000000e+00 10
#> 190 0.80000000 0.48000000 -6.470588e-01 0.000000e+00 10
#> 191 0.35000000 0.28875000 -1.875000e-01 0.000000e+00 10
#> 192 0.53846154 0.39349112 -3.333333e-01 0.000000e+00 10
#> 193 0.54545455 0.39669421 -3.333333e-01 0.000000e+00  8
#> 194 0.75000000 0.46875000 -5.000000e-01 1.910848e-08  8
#> 195 0.42105263 0.46537396  1.219512e-01 0.000000e+00  8
#> 196 0.46666667 0.35777777 -2.727273e-01 0.000000e+00  8
#> 197 0.31818182 0.44938016  3.130841e-01 0.000000e+00  8
#> 198 0.23076923 0.45266272  5.200000e-01 0.000000e+00  8
#> 199 0.77777778 0.47530864 -6.000000e-01 0.000000e+00 12
#> 200 0.33333333 0.50000000  4.117647e-01 0.000000e+00 12
#> 201 0.55555556 0.47530864 -1.409396e-01 0.000000e+00 12
#> 202 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 203 0.31818182 0.41632232  2.575758e-01 0.000000e+00 12
#> 204 0.27272727 0.23553720 -1.111111e-01 0.000000e+00 12
#> 205 0.20000000 0.18000000 -5.882353e-02 5.735292e-02 10
#> 206 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 207 0.13636364 0.12706612 -5.000000e-02 2.159460e-02 10
#> 208 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 209 0.15000000 0.13875000 -5.555556e-02 2.551508e-02 10
#> 210 0.14285714 0.13265306 -4.000000e-02 2.617527e-02 10
#> 211 0.30000000 0.25500000 -1.250000e-01 0.000000e+00 12
#> 212 0.33333333 0.27777778 -1.111111e-01 4.582176e-05 12
#> 213 0.29166667 0.30468750  6.395349e-02 0.000000e+00 12
#> 214 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 215 0.36363636 0.29752066 -2.000000e-01 0.000000e+00 12
#> 216 0.15384615 0.14201184 -4.347826e-02 0.000000e+00 12
#> 217 0.20000000 0.32000000  4.193548e-01 0.000000e+00 10
#> 218 0.40000000 0.32000000 -1.428571e-01 7.692971e-08 10
#> 219 0.54166667 0.39496527 -3.529412e-01 0.000000e+00 10
#> 220 0.06250000 0.06054688  0.000000e+00 0.000000e+00 10
#> 221 0.45454545 0.46280992  4.109589e-02 0.000000e+00 10
#> 222 0.53846154 0.39349112 -3.333333e-01 0.000000e+00 10
#> 223 0.30000000 0.49500000  4.375000e-01 2.451686e-11 12
#> 224 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 225 0.39130435 0.31474480 -2.222222e-01 0.000000e+00 12
#> 226 0.68750000 0.48242188 -3.983051e-01 0.000000e+00 12
#> 227 0.63157895 0.48753463 -2.705882e-01 0.000000e+00 12
#> 228 0.21428571 0.37500000  4.583333e-01 0.000000e+00 12
#> 229 0.33333333 0.27777778 -1.428571e-01 0.000000e+00  6
#> 230 0.66666667 0.44444444 -3.333333e-01 8.543447e-03  6
#> 231 0.20000000 0.18000000 -7.692308e-02 0.000000e+00  6
#> 232 0.27272727 0.23553720 -1.111111e-01 0.000000e+00  6
#> 233 0.17647059 0.16089965 -6.666667e-02 0.000000e+00  6
#> 234 0.30000000 0.25500000 -1.250000e-01 0.000000e+00  6
#> 235 0.50000000 0.49500000  4.255319e-02 0.000000e+00 10
#> 236 0.40000000 0.32000000 -1.428571e-01 4.269096e-11 10
#> 237 0.35000000 0.28875000 -1.875000e-01 0.000000e+00 10
#> 238 0.33333333 0.27777778 -1.666667e-01 0.000000e+00 10
#> 239 0.31578947 0.26592798 -1.612903e-01 0.000000e+00 10
#> 240 0.45454545 0.48347107  1.071429e-01 0.000000e+00 10
#> 241 0.63636364 0.48347107 -2.727273e-01 0.000000e+00 10
#> 242 0.60000000 0.42000000 -3.333333e-01 2.803658e-23 10
#> 243 0.52631579 0.38781163 -3.333333e-01 0.000000e+00 10
#> 244 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 245 0.52941176 0.43771626 -1.803279e-01 0.000000e+00 10
#> 246 0.40000000 0.48000000  2.173913e-01 0.000000e+00 10
#> 247 0.50000000 0.37500000 -2.727273e-01 0.000000e+00 10
#> 248 0.40000000 0.48000000  2.727273e-01 3.821019e-13 10
#> 249 0.47826087 0.49621928  5.836576e-02 0.000000e+00 10
#> 250 0.12500000 0.11718750 -3.448276e-02 0.000000e+00 10
#> 251 0.68181818 0.49070248 -3.695652e-01 0.000000e+00 10
#> 252 0.57142857 0.45918368 -2.093023e-01 0.000000e+00 10
#> 253 0.18181818 0.16528925 -5.263158e-02 0.000000e+00 10
#> 254 0.40000000 0.32000000 -1.428571e-01 2.596087e-04 10
#> 255 0.38095238 0.30839002 -2.121212e-01 0.000000e+00 10
#> 256 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 257 0.20833333 0.18663195 -9.523810e-02 0.000000e+00 10
#> 258 0.46153846 0.35502958 -2.631579e-01 0.000000e+00 10
#> 259 0.11111111 0.27777778  6.363636e-01 2.921766e-07 12
#> 260 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 261 0.11538462 0.10872781 -4.166667e-02 9.079808e-09 12
#> 262 0.53333333 0.48000000 -7.692308e-02 3.419635e-05 12
#> 263 0.14285714 0.20975057  3.406593e-01 3.648059e-08 12
#> 264 0.41666667 0.32986111 -2.222222e-01 3.866324e-07 12
#> 265 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 266 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 267 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 268 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 269 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 270 0.13333333 0.12444445 -3.703704e-02 6.482759e-01 12
#> 271 0.12500000 0.49218750  7.741935e-01 0.000000e+00 10
#> 272 0.40000000 0.32000000 -1.428571e-01 0.000000e+00 10
#> 273 0.50000000 0.48611111  0.000000e+00 0.000000e+00 10
#> 274 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 275 0.17647059 0.32698962  4.838710e-01 0.000000e+00 10
#> 276 0.38461538 0.45266272  1.891892e-01 0.000000e+00 10
#> 277 0.45454545 0.35123967 -2.500000e-01 1.783739e-08 12
#> 278 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 279 0.60869565 0.49905482 -1.984436e-01 2.479688e-06 12
#> 280 0.20000000 0.18000000 -7.692308e-02 6.851436e-10 12
#> 281 0.30434783 0.36389413  1.851852e-01 7.863988e-09 12
#> 282 0.20000000 0.27777778  3.114754e-01 2.694109e-09 12
#> 283 0.33333333 0.40123457  2.258065e-01 0.000000e+00 12
#> 284 0.66666667 0.50000000 -2.500000e-01 2.992217e-21 12
#> 285 0.55000000 0.48875000 -1.000000e-01 0.000000e+00 12
#> 286 0.35714286 0.47704082  2.857143e-01 0.000000e+00 12
#> 287 0.55555556 0.49382716 -9.677419e-02 0.000000e+00 12
#> 288 0.61538462 0.47337278 -2.631579e-01 0.000000e+00 12
#> 289 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 290 0.33333333 0.27777778 -1.111111e-01 2.177778e-01 12
#> 291 0.03703704 0.03635117 -6.661338e-16 0.000000e+00 12
#> 292 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 293 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 294 0.15384615 0.14201184 -4.347826e-02 0.000000e+00 12
#> 295 0.30000000 0.25500000 -1.250000e-01 0.000000e+00 10
#> 296 0.20000000 0.18000000 -2.220446e-16 4.331694e-03 10
#> 297 0.21739130 0.19376181 -1.000000e-01 0.000000e+00 10
#> 298 0.06250000 0.06054688  0.000000e+00 0.000000e+00 10
#> 299 0.10526316 0.09972299 -2.857143e-02 0.000000e+00 10
#> 300 0.14285714 0.13265306 -4.000000e-02 0.000000e+00 10
#> 301 0.55555556 0.47530864 -1.111111e-01 4.422482e-09 10
#> 302 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 303 0.50000000 0.45500000 -7.344633e-02 1.763462e-10 10
#> 304 0.36363636 0.29752066 -1.764706e-01 1.532109e-11 10
#> 305 0.50000000 0.43622449 -1.097561e-01 1.568595e-10 10
#> 306 0.50000000 0.47704082 -1.111111e-02 7.492525e-10 10
#> 307 0.72727273 0.46280992 -5.384615e-01 0.000000e+00 12
#> 308 0.33333333 0.44444444  3.333333e-01 1.516923e-16 12
#> 309 0.50000000 0.44444444 -1.040000e-01 0.000000e+00 12
#> 310 0.75000000 0.49652778 -4.776119e-01 0.000000e+00 12
#> 311 0.20833333 0.24913194  1.843972e-01 0.000000e+00 12
#> 312 0.42857143 0.45918368  1.034483e-01 0.000000e+00 12
#> 313 0.20000000 0.32000000  4.193548e-01 0.000000e+00  8
#> 314 0.75000000 0.46875000 -5.000000e-01 2.843110e-06  8
#> 315 0.52380952 0.42743764 -2.021858e-01 0.000000e+00  8
#> 316 0.00000000 0.00000000           NaN 0.000000e+00  8
#> 317 0.26086957 0.28733459  1.140940e-01 0.000000e+00  8
#> 318 0.50000000 0.48611111  1.492537e-02 0.000000e+00  8
#> 319 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 320 0.16666667 0.15277777  3.330669e-16 5.685619e-01 12
#> 321 0.03846154 0.03772189  1.110223e-15 0.000000e+00 12
#> 322 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 323 0.04347826 0.04253308 -2.220446e-16 0.000000e+00 12
#> 324 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 325 0.75000000 0.50000000 -4.482759e-01 0.000000e+00 12
#> 326 0.50000000 0.37500000 -2.500000e-01 0.000000e+00 12
#> 327 0.38095238 0.44444444  1.666667e-01 0.000000e+00 12
#> 328 0.66666667 0.49777778 -3.084112e-01 0.000000e+00 12
#> 329 0.25000000 0.28875000  1.592920e-01 0.000000e+00 12
#> 330 0.46153846 0.49704142  1.111111e-01 0.000000e+00 12
#> 331 0.66666667 0.44444444 -4.545455e-01 0.000000e+00 12
#> 332 0.16666667 0.15277777  3.330669e-16 2.309141e-12 12
#> 333 0.40000000 0.32000000 -2.258065e-01 0.000000e+00 12
#> 334 0.68750000 0.45117188 -5.000000e-01 0.000000e+00 12
#> 335 0.36842105 0.30055402 -2.000000e-01 0.000000e+00 12
#> 336 0.46666667 0.35777777 -2.727273e-01 0.000000e+00 12
#> 337 0.12500000 0.11718750  0.000000e+00 1.353574e-01 10
#> 338 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 339 0.20000000 0.18000000 -8.571429e-02 1.896189e-01 10
#> 340 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 341 0.05555556 0.05401235  2.220446e-16 3.123632e-02 10
#> 342 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 343 0.44444444 0.49382716  1.578947e-01 2.281271e-10 12
#> 344 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 345 0.52000000 0.38480000 -3.333333e-01 0.000000e+00 12
#> 346 0.38461538 0.48816568  2.500000e-01 0.000000e+00 12
#> 347 0.27777778 0.23919753 -1.333333e-01 0.000000e+00 12
#> 348 0.36363636 0.39669421  1.304348e-01 0.000000e+00 12
#> 349 0.36363636 0.49586777  3.103448e-01 0.000000e+00 12
#> 350 0.66666667 0.44444444 -4.285714e-01 4.374982e-18 12
#> 351 0.62500000 0.47829861 -2.873134e-01 0.000000e+00 12
#> 352 0.12500000 0.11718750 -3.448276e-02 0.000000e+00 12
#> 353 0.65000000 0.49875000 -2.797927e-01 0.000000e+00 12
#> 354 0.46153846 0.47337278  6.493506e-02 0.000000e+00 12
#> 355 0.57142857 0.48979592 -9.090909e-02 0.000000e+00 12
#> 356 0.50000000 0.37500000 -2.500000e-01 0.000000e+00 12
#> 357 0.60869565 0.45368620 -3.218884e-01 0.000000e+00 12
#> 358 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 359 0.47619048 0.44444444 -4.712042e-02 0.000000e+00 12
#> 360 0.64285714 0.43622449 -4.444444e-01 0.000000e+00 12
#> 361 0.00000000 0.00000000           NaN 0.000000e+00  6
#> 362 0.33333333 0.27777778 -2.220446e-16 3.283445e-01  6
#> 363 0.05882353 0.05709342  1.110223e-16 0.000000e+00  6
#> 364 0.00000000 0.00000000           NaN 0.000000e+00  6
#> 365 0.10000000 0.09500000 -2.702703e-02 0.000000e+00  6
#> 366 0.16666667 0.15277777 -4.761905e-02 0.000000e+00  6
#> 367 0.20000000 0.18000000 -5.882353e-02 2.302740e-03 12
#> 368 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 369 0.42307692 0.33357988 -2.500000e-01 1.000122e-02 12
#> 370 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 371 0.12500000 0.11718750 -4.545455e-02 5.661250e-04 12
#> 372 0.21428571 0.19132654 -8.333333e-02 1.930321e-03 12
#> 373 0.50000000 0.46875000  0.000000e+00 0.000000e+00 10
#> 374 0.40000000 0.48000000  2.727273e-01 0.000000e+00 10
#> 375 0.44444444 0.34567901 -2.592593e-01 0.000000e+00 10
#> 376 0.40000000 0.32000000 -2.000000e-01 0.000000e+00 10
#> 377 0.50000000 0.42968750 -1.320755e-01 0.000000e+00 10
#> 378 0.38461538 0.39349112  6.250000e-02 0.000000e+00 10
#> 379 0.45454545 0.43388430  1.110223e-16 0.000000e+00 12
#> 380 0.00000000 0.27777778  1.000000e+00 9.121187e-20 12
#> 381 0.61904762 0.42743764 -4.285714e-01 0.000000e+00 12
#> 382 0.42857143 0.48979592  1.612903e-01 0.000000e+00 12
#> 383 0.45454545 0.46280992  4.109589e-02 0.000000e+00 12
#> 384 0.60000000 0.46444445 -2.600000e-01 0.000000e+00 12
#> 385 0.36363636 0.29752066 -1.764706e-01 0.000000e+00 12
#> 386 0.50000000 0.48611111  6.250000e-02 4.579921e-18 12
#> 387 0.42105263 0.46537396  1.219512e-01 0.000000e+00 12
#> 388 0.41176471 0.38927335 -2.752294e-02 0.000000e+00 12
#> 389 0.40909091 0.44938016  1.126761e-01 0.000000e+00 12
#> 390 0.35714286 0.47704082  2.857143e-01 0.000000e+00 12
#> 391 0.55555556 0.50000000 -5.263158e-02 0.000000e+00 12
#> 392 0.50000000 0.48611111  6.250000e-02 0.000000e+00 12
#> 393 0.60000000 0.45500000 -2.954545e-01 0.000000e+00 12
#> 394 0.43750000 0.45117188  6.250000e-02 0.000000e+00 12
#> 395 0.41176471 0.47231834  1.578947e-01 0.000000e+00 12
#> 396 0.38461538 0.48816568  2.500000e-01 0.000000e+00 12
#> 397 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 398 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 399 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 400 0.40000000 0.39111111  1.176471e-02 5.344359e-01 12
#> 401 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 402 0.07692308 0.07396450 -4.440892e-16 3.450582e-03 12
#> 403 0.30000000 0.37500000  2.500000e-01 0.000000e+00 12
#> 404 0.50000000 0.37500000 -2.500000e-01 0.000000e+00 12
#> 405 0.37500000 0.35329861 -4.020101e-02 0.000000e+00 12
#> 406 0.38461538 0.39349112  6.250000e-02 0.000000e+00 12
#> 407 0.44444444 0.40123457 -7.936508e-02 0.000000e+00 12
#> 408 0.42857143 0.45918368  1.034483e-01 0.000000e+00 12
#> 409 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 410 0.20000000 0.18000000 -2.220446e-16 5.741566e-03 10
#> 411 0.36000000 0.42000000  1.627907e-01 0.000000e+00 10
#> 412 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 413 0.19047619 0.17233561 -8.108108e-02 0.000000e+00 10
#> 414 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 415 0.70000000 0.45500000 -5.000000e-01 0.000000e+00 10
#> 416 0.40000000 0.48000000  2.727273e-01 1.725981e-08 10
#> 417 0.52380952 0.38662131 -3.333333e-01 0.000000e+00 10
#> 418 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 419 0.50000000 0.37500000 -3.125000e-01 0.000000e+00 10
#> 420 0.58333333 0.41319445 -3.750000e-01 0.000000e+00 10
#> 421 0.09090909 0.08677687 -4.440892e-16 1.013035e-03 12
#> 422 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 423 0.03703704 0.03635117 -6.661338e-16 2.411987e-04 12
#> 424 0.37500000 0.30468750 -2.000000e-01 1.889273e-02 12
#> 425 0.25000000 0.21875000 -1.219512e-01 4.474420e-03 12
#> 426 0.13333333 0.12444445 -3.703704e-02 1.555968e-03 12
#> 427 0.30000000 0.45500000  3.863636e-01 0.000000e+00 12
#> 428 0.50000000 0.48611111  6.250000e-02 0.000000e+00 12
#> 429 0.54166667 0.42968750 -2.406639e-01 0.000000e+00 12
#> 430 0.46666667 0.46444445  2.970297e-02 0.000000e+00 12
#> 431 0.45454545 0.35123967 -2.727273e-01 0.000000e+00 12
#> 432 0.53846154 0.39349112 -3.333333e-01 0.000000e+00 12
#> 433 0.60000000 0.42000000 -3.846154e-01 0.000000e+00 12
#> 434 0.50000000 0.37500000 -2.500000e-01 5.705162e-10 12
#> 435 0.26086957 0.22684310 -1.282051e-01 0.000000e+00 12
#> 436 0.46666667 0.42000000 -7.692308e-02 0.000000e+00 12
#> 437 0.23809524 0.20975057 -1.111111e-01 0.000000e+00 12
#> 438 0.35714286 0.37500000  8.450704e-02 0.000000e+00 12
#> 439 0.30000000 0.25500000 -1.250000e-01 0.000000e+00 10
#> 440 0.20000000 0.18000000 -2.220446e-16 1.031933e-05 10
#> 441 0.63157895 0.43213296 -4.400000e-01 0.000000e+00 10
#> 442 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 443 0.47058824 0.35986159 -2.800000e-01 0.000000e+00 10
#> 444 0.23076923 0.20414202 -9.090909e-02 0.000000e+00 10
#> 445 0.66666667 0.49382716 -2.972973e-01 0.000000e+00 10
#> 446 0.60000000 0.50000000 -9.090909e-02 2.054125e-23 10
#> 447 0.60869565 0.49149338 -2.173913e-01 0.000000e+00 10
#> 448 0.62500000 0.42968750 -4.285714e-01 0.000000e+00 10
#> 449 0.55000000 0.49875000 -7.731959e-02 0.000000e+00 10
#> 450 0.60000000 0.42000000 -4.000000e-01 0.000000e+00 10
#> 451 0.50000000 0.48611111  6.250000e-02 0.000000e+00  6
#> 452 0.33333333 0.27777778 -2.220446e-16 0.000000e+00  6
#> 453 0.18750000 0.45117188  6.052632e-01 0.000000e+00  6
#> 454 0.00000000 0.00000000           NaN 0.000000e+00  6
#> 455 0.25000000 0.21875000 -1.000000e-01 0.000000e+00  6
#> 456 0.42857143 0.48979592  1.612903e-01 0.000000e+00  6
#> 457 0.66666667 0.44444444 -4.545455e-01 0.000000e+00 12
#> 458 0.16666667 0.15277777  3.330669e-16 2.419234e-12 12
#> 459 0.34782609 0.28733459 -1.891892e-01 0.000000e+00 12
#> 460 0.73333333 0.46444445 -5.555556e-01 0.000000e+00 12
#> 461 0.36842105 0.30055402 -2.000000e-01 0.000000e+00 12
#> 462 0.42857143 0.33673469 -2.380952e-01 0.000000e+00 12
#> 463 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 464 0.16666667 0.37500000  6.153846e-01 1.240721e-03 12
#> 465 0.40909091 0.41632232  4.060914e-02 0.000000e+00 12
#> 466 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 467 0.22727273 0.20144628 -1.052632e-01 0.000000e+00 12
#> 468 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 469 0.36363636 0.29752066 -1.764706e-01 9.286520e-05 12
#> 470 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 471 0.23076923 0.20414202 -1.111111e-01 1.181337e-05 12
#> 472 0.47058824 0.35986159 -2.800000e-01 1.531012e-04 12
#> 473 0.09090909 0.16528925  4.683544e-01 7.616929e-06 12
#> 474 0.25000000 0.21875000 -1.000000e-01 2.250723e-05 12
#> 475 0.27272727 0.35123967  2.682927e-01 0.000000e+00 12
#> 476 0.33333333 0.27777778 -1.111111e-01 1.517927e-14 12
#> 477 0.54545455 0.43388430 -2.352941e-01 0.000000e+00 12
#> 478 0.73333333 0.50000000 -4.392523e-01 0.000000e+00 12
#> 479 0.52173913 0.38563327 -3.333333e-01 0.000000e+00 12
#> 480 0.53846154 0.39349112 -3.333333e-01 0.000000e+00 12
#> 481 0.66666667 0.44444444 -4.545455e-01 0.000000e+00 10
#> 482 0.80000000 0.48000000 -6.000000e-01 4.078311e-14 10
#> 483 0.61904762 0.48185941 -2.621359e-01 0.000000e+00 10
#> 484 0.62500000 0.49218750 -2.396694e-01 0.000000e+00 10
#> 485 0.54545455 0.43388430 -2.352941e-01 0.000000e+00 10
#> 486 0.64285714 0.47704082 -3.146067e-01 0.000000e+00 10
#> 487 0.20000000 0.18000000 -5.882353e-02 4.570016e-02  6
#> 488 0.00000000 0.00000000           NaN 0.000000e+00  6
#> 489 0.19230769 0.17381657 -8.695652e-02 3.763879e-02  6
#> 490 0.00000000 0.00000000           NaN 0.000000e+00  6
#> 491 0.10526316 0.09972299 -2.857143e-02 1.753124e-02  6
#> 492 0.38461538 0.31065088 -2.000000e-01 1.362200e-01  6
#> 493 0.36363636 0.29752066 -1.764706e-01 0.000000e+00 12
#> 494 0.83333333 0.48611111 -6.666667e-01 1.230349e-02 12
#> 495 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 496 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 497 0.04166667 0.04079860  7.771561e-16 0.000000e+00 12
#> 498 0.07142857 0.06887756 -2.220446e-16 0.000000e+00 12
#> 499 0.40000000 0.42000000  1.000000e-01 0.000000e+00 12
#> 500 0.33333333 0.44444444  3.333333e-01 0.000000e+00 12
#> 501 0.50000000 0.45500000 -7.344633e-02 0.000000e+00 12
#> 502 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 503 0.42105263 0.38781163 -5.882353e-02 0.000000e+00 12
#> 504 0.61538462 0.49704142 -2.000000e-01 0.000000e+00 12
#> 505 0.50000000 0.45500000 -4.651163e-02 0.000000e+00 10
#> 506 0.80000000 0.48000000 -6.000000e-01 5.333094e-24 10
#> 507 0.61904762 0.42743764 -4.285714e-01 0.000000e+00 10
#> 508 0.21428571 0.19132654 -8.333333e-02 0.000000e+00 10
#> 509 0.42857143 0.38662131 -8.433735e-02 0.000000e+00 10
#> 510 0.60000000 0.49111111 -1.886792e-01 0.000000e+00 10
#> 511 0.75000000 0.46875000 -5.555556e-01 0.000000e+00 12
#> 512 0.16666667 0.15277777  3.330669e-16 0.000000e+00 12
#> 513 0.45000000 0.34875000 -2.666667e-01 0.000000e+00 12
#> 514 0.81818182 0.48347107 -6.666667e-01 0.000000e+00 12
#> 515 0.36842105 0.30055402 -2.000000e-01 0.000000e+00 12
#> 516 0.53846154 0.39349112 -3.333333e-01 0.000000e+00 12
#> 517 0.45454545 0.50000000  1.379310e-01 0.000000e+00 12
#> 518 0.16666667 0.15277777  3.330669e-16 0.000000e+00 12
#> 519 0.32000000 0.36480000  1.428571e-01 0.000000e+00 12
#> 520 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 521 0.54545455 0.39669421 -3.548387e-01 0.000000e+00 12
#> 522 0.57142857 0.48979592 -1.304348e-01 0.000000e+00 12
#> 523 0.75000000 0.50000000 -4.482759e-01 0.000000e+00 12
#> 524 0.50000000 0.37500000 -2.500000e-01 0.000000e+00 12
#> 525 0.45000000 0.43875000  1.110223e-16 0.000000e+00 12
#> 526 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 527 0.68750000 0.48242188 -3.983051e-01 0.000000e+00 12
#> 528 0.66666667 0.44444444 -4.666667e-01 0.000000e+00 12
#> 529 0.66666667 0.44444444 -4.545455e-01 0.000000e+00  6
#> 530 0.00000000 0.44444444  1.000000e+00 1.431960e-09  6
#> 531 0.50000000 0.46141975 -5.517241e-02 0.000000e+00  6
#> 532 0.16666667 0.50000000  6.901408e-01 0.000000e+00  6
#> 533 0.50000000 0.46141975 -5.517241e-02 0.000000e+00  6
#> 534 0.33333333 0.48611111  3.529412e-01 0.000000e+00  6
#> 535 0.27272727 0.23553720 -1.111111e-01 1.908300e-01 12
#> 536 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 537 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 538 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 539 0.04347826 0.04253308 -2.220446e-16 5.691700e-03 12
#> 540 0.15384615 0.14201184 -4.347826e-02 4.146810e-02 12
#> 541 0.33333333 0.27777778 -1.428571e-01 0.000000e+00  8
#> 542 0.25000000 0.46875000  5.714286e-01 8.287744e-10  8
#> 543 0.61904762 0.50000000 -2.149533e-01 0.000000e+00  8
#> 544 0.18750000 0.16992188 -7.142857e-02 0.000000e+00  8
#> 545 0.50000000 0.49500000  1.554404e-02 0.000000e+00  8
#> 546 0.42857143 0.45918368  1.034483e-01 0.000000e+00  8
#> 547 0.42857143 0.50000000  2.173913e-01 0.000000e+00 10
#> 548 0.40000000 0.32000000 -1.428571e-01 0.000000e+00 10
#> 549 0.71428571 0.45918368 -5.294118e-01 0.000000e+00 10
#> 550 0.54545455 0.46280992 -1.320755e-01 0.000000e+00 10
#> 551 0.36363636 0.29752066 -1.764706e-01 0.000000e+00 10
#> 552 0.40000000 0.32000000 -2.000000e-01 0.000000e+00 10
#> 553 0.27272727 0.23553720 -1.111111e-01 0.000000e+00 12
#> 554 0.50000000 0.37500000 -2.500000e-01 1.559402e-06 12
#> 555 0.11538462 0.10872781 -4.166667e-02 0.000000e+00 12
#> 556 0.28571429 0.24489796 -1.304348e-01 0.000000e+00 12
#> 557 0.36842105 0.30055402 -2.000000e-01 0.000000e+00 12
#> 558 0.50000000 0.37500000 -3.000000e-01 0.000000e+00 12
#> 559 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 560 0.00000000 0.00000000           NaN 0.000000e+00 10
#> 561 0.13636364 0.12706612 -5.000000e-02 3.863034e-04 10
#> 562 0.53846154 0.45266272 -1.506849e-01 8.664668e-02 10
#> 563 0.04545455 0.04442148 -4.440892e-16 9.364930e-05 10
#> 564 0.21428571 0.19132654 -8.333333e-02 9.599054e-04 10
#> 565 0.40000000 0.32000000 -2.000000e-01 0.000000e+00 12
#> 566 0.33333333 0.44444444  3.333333e-01 3.331096e-06 12
#> 567 0.18181818 0.29752066  4.084507e-01 0.000000e+00 12
#> 568 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 569 0.31578947 0.26592798 -1.612903e-01 0.000000e+00 12
#> 570 0.33333333 0.27777778 -1.578947e-01 0.000000e+00 12
#> 571 0.12500000 0.11718750  0.000000e+00 0.000000e+00 12
#> 572 0.16666667 0.15277777  3.330669e-16 1.384227e-03 12
#> 573 0.34782609 0.28733459 -1.891892e-01 0.000000e+00 12
#> 574 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 575 0.17391304 0.15879017 -7.317073e-02 0.000000e+00 12
#> 576 0.15384615 0.14201184 -4.347826e-02 0.000000e+00 12
#> 577 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 578 0.16666667 0.15277777  3.330669e-16 1.469043e-01 12
#> 579 0.09523810 0.09070295 -2.564103e-02 0.000000e+00 12
#> 580 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 581 0.10000000 0.09500000 -2.702703e-02 0.000000e+00 12
#> 582 0.06666667 0.06444444 -2.220446e-16 0.000000e+00 12
#> 583 0.22222222 0.34567901  4.074074e-01 0.000000e+00 10
#> 584 0.20000000 0.42000000  6.000000e-01 2.904386e-09 10
#> 585 0.19230769 0.17381657 -8.695652e-02 0.000000e+00 10
#> 586 0.56250000 0.48242188 -1.344538e-01 0.000000e+00 10
#> 587 0.40909091 0.32541322 -2.352941e-01 0.000000e+00 10
#> 588 0.53846154 0.45266272 -1.506849e-01 0.000000e+00 10
#> 589 0.50000000 0.49500000  4.255319e-02 0.000000e+00 12
#> 590 0.16666667 0.48611111  7.058824e-01 0.000000e+00 12
#> 591 0.66666667 0.48979592 -3.397129e-01 0.000000e+00 12
#> 592 0.00000000 0.00000000           NaN 0.000000e+00 12
#> 593 0.45000000 0.49875000  1.230769e-01 0.000000e+00 12
#> 594 0.45454545 0.50000000  1.379310e-01 0.000000e+00 12
#> 595 0.40000000 0.48000000  2.173913e-01 0.000000e+00  8
#> 596 0.75000000 0.46875000 -5.000000e-01 9.894405e-15  8
#> 597 0.40000000 0.42000000  7.317073e-02 0.000000e+00  8
#> 598 0.45454545 0.43388430  1.110223e-16 0.000000e+00  8
#> 599 0.65000000 0.43875000 -4.615385e-01 0.000000e+00  8
#> 600 0.50000000 0.47704082 -1.111111e-02 0.000000e+00  8
get.snpR.stats(x, "pop", stats = "fst") # view fst
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

```
