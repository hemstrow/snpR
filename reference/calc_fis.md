# Calculate FIS for individual populations.

Calculates FIS for each individual level in each provided sample facet
according to Weir and Cockerham (1984).

## Usage

``` r
calc_fis(
  x,
  facets = NULL,
  boot = FALSE,
  boot_par = FALSE,
  boot_alt = "two-sided",
  boot_confidence = 0.95,
  keep_components = FALSE
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

- boot:

  numeric or `FALSE`, default `FALSE`. The number of bootstraps to do.
  See details.

- boot_par:

  numeric or `FALSE`, default `FALSE`. If a number, bootstraps will be
  processed in parallel using the supplied number of cores.

- boot_alt:

  character, default "two-sided". The type of *t*-test to conduct on the
  bootstrapped FIS values. Options:

  - two-sided: A two-sided *t*-test, with the alternative hypothesis
    that FIS is not equal to 0.

  - greater: A one-sided *t*-test, with the alternative hypothesis that
    FIS is greater than 0.

  - less: A one-sided *t*-test, with the alternative hypothesis that FIS
    is less than 0.

- boot_confidence:

  numeric, default 0.95. The confidence level to use for generating
  confidence levels for bootstraps. Note that a value of 0.95 will
  result in a call of [`qt`](https://rdrr.io/r/stats/TDist.html) with
  `q = .975` and and so on.

- keep_components:

  logical, default `FALSE`. If TRUE, the variance components "b" and "c"
  will be held and accessible from the `$single` element (named
  "var_comp_b" and "var_comp_c", respectively) using the usual
  [`get.snpR.stats`](https://hemstrow.github.io/snpR/reference/get.snpR.stats.md)
  method. This may be useful if working with very large datasets that
  need to be run with separate objects for each chromosome, etc. for
  memory purposes. Weighted averages can be generated identically to
  those from snpR by taking one minus the weighted mean (via the
  [`weighted.mean`](https://rdrr.io/r/stats/weighted.mean.html)) of "c"
  divided by the sum of the weighted mean of "b" + "c" using the number
  of SNPs called in a comparison (returned in the "nk" column from
  [`get.snpR.stats`](https://hemstrow.github.io/snpR/reference/get.snpR.stats.md))
  as weights within each population comparison. Note that this is
  different than taking the weighted mean of 1 - c/(b + c)!

## Details

Note that FIS is calculated by considering *only data from individual
sample levels*! This means that individual and sub-population variances
are only considered within each sub-population. If snp facets are
provided, weighted means will be provided for each snp facet level,
although raw values are calculated on a per-snp basis and thus ignore
these levels.

If the base facet (facets = NULL or facets = ".base") is requested, FIS
will compare individual to total variance across all samples instead
(equivalent to overall FIT).

Bootstrapping across loci can be done to assess FIS significance. This
is done by re-drawing loci randomly with replacement for each facet
level, calculating the resulting FIS values, and doing one sample
*t*-test with the null hypothesis that FIS = 0 to calculate p-values and
confidence intervals.

## References

Weir and Cockerham (1984). *Evolution*

## Author

William Hemstrom

## Examples

``` r
x <- calc_fis(stickSNPs, c("pop", "pop.chr"))
get.snpR.stats(x, c("pop", "pop.chr"), "fis")
#> $single
#>     facet subfacet        chr position .snp.id           fis
#> 1     pop      ASP    groupVI   212436       1  3.846154e-01
#> 2     pop      CLF    groupVI   212436       1 -3.333333e-01
#> 3     pop      OPL    groupVI   212436       1 -3.068783e-01
#> 4     pop      PAL    groupVI   212436       1  1.304348e-01
#> 5     pop      SMR    groupVI   212436       1  8.296943e-02
#> 6     pop      UPD    groupVI   212436       1 -3.703704e-02
#> 7     pop      ASP groupXVIII   350070       2  5.551115e-16
#> 8     pop      OPL groupXVIII   350070       2 -4.761905e-02
#> 9     pop      SMR groupXVIII   350070       2 -2.439024e-02
#> 10    pop      ASP    groupXV   854285       3 -5.263158e-02
#> 11    pop      CLF    groupXV   854285       3  3.330669e-16
#> 12    pop      OPL    groupXV   854285       3 -6.666667e-02
#> 13    pop      PAL    groupXV   854285       3  0.000000e+00
#> 14    pop      SMR    groupXV   854285       3 -4.761905e-02
#> 15    pop      UPD    groupXV   854285       3 -3.703704e-02
#> 16    pop      ASP  groupXIII   907279       4 -1.111111e-01
#> 17    pop      CLF  groupXIII   907279       4 -2.220446e-16
#> 18    pop      OPL  groupXIII   907279       4 -7.272727e-01
#> 19    pop      SMR  groupXIII   907279       4  1.315789e-01
#> 20    pop      UPD  groupXIII   907279       4 -3.846154e-01
#> 21    pop      ASP    groupVI  1068620       5  1.428571e-01
#> 22    pop      OPL    groupVI  1068620       5 -4.347826e-02
#> 23    pop      PAL    groupVI  1068620       5 -2.222222e-01
#> 24    pop      SMR    groupVI  1068620       5  1.025641e-01
#> 25    pop      UPD    groupVI  1068620       5  1.578947e-01
#> 26    pop      ASP     groupV  1392753       6 -2.857143e-01
#> 27    pop      CLF     groupV  1392753       6 -2.500000e-01
#> 28    pop      OPL     groupV  1392753       6 -4.166667e-02
#> 29    pop      PAL     groupV  1392753       6 -3.983051e-01
#> 30    pop      SMR     groupV  1392753       6 -5.263158e-02
#> 31    pop      UPD     groupV  1392753       6 -2.093023e-01
#> 32    pop      ASP    groupIX  2075079       7 -4.440892e-16
#> 33    pop      CLF    groupIX  2075079       7  1.000000e+00
#> 34    pop      OPL    groupIX  2075079       7 -1.428571e-01
#> 35    pop      SMR    groupIX  2075079       7 -2.121212e-01
#> 36    pop      ASP  groupXIII  2106000       8 -2.500000e-01
#> 37    pop      CLF  groupXIII  2106000       8  6.250000e-02
#> 38    pop      OPL  groupXIII  2106000       8  1.210191e-01
#> 39    pop      SMR  groupXIII  2106000       8 -7.344633e-02
#> 40    pop      UPD  groupXIII  2106000       8 -2.857143e-01
#> 41    pop      ASP   groupXII  2136587       9  2.682927e-01
#> 42    pop      CLF   groupXII  2136587       9 -1.111111e-01
#> 43    pop      OPL   groupXII  2136587       9 -9.090909e-02
#> 44    pop      PAL   groupXII  2136587       9  0.000000e+00
#> 45    pop      SMR   groupXII  2136587       9 -3.103448e-01
#> 46    pop      UPD   groupXII  2136587       9  6.250000e-02
#> 47    pop      ASP  groupXVII  2560839      10 -5.263158e-02
#> 48    pop      CLF  groupXVII  2560839      10  3.330669e-16
#> 49    pop      OPL  groupXVII  2560839      10  5.551115e-16
#> 50    pop      SMR  groupXVII  2560839      10 -2.439024e-02
#> 51    pop      UPD  groupXVII  2560839      10 -1.304348e-01
#> 52    pop      ASP   groupIII  2750695      11  4.255319e-02
#> 53    pop      CLF   groupIII  2750695      11 -2.500000e-01
#> 54    pop      OPL   groupIII  2750695      11 -6.280193e-02
#> 55    pop      SMR   groupIII  2750695      11 -1.314286e-01
#> 56    pop      UPD   groupIII  2750695      11 -6.250000e-01
#> 57    pop      ASP   groupXIX  3008524      12 -2.857143e-01
#> 58    pop      CLF   groupXIX  3008524      12  3.330669e-16
#> 59    pop      OPL   groupXIX  3008524      12 -2.142857e-01
#> 60    pop      PAL   groupXIX  3008524      12 -4.285714e-01
#> 61    pop      SMR   groupXIX  3008524      12 -1.851852e-01
#> 62    pop      UPD   groupXIX  3008524      12 -3.846154e-01
#> 63    pop      ASP   groupXVI  3624417      13 -6.666667e-02
#> 64    pop      CLF   groupXVI  3624417      13  6.250000e-02
#> 65    pop      OPL   groupXVI  3624417      13  1.320755e-01
#> 66    pop      PAL   groupXVI  3624417      13 -2.220446e-16
#> 67    pop      SMR   groupXVI  3624417      13  8.800000e-02
#> 68    pop      UPD   groupXVI  3624417      13 -1.428571e-01
#> 69    pop      ASP   groupVII  3667528      14 -2.000000e-01
#> 70    pop      CLF   groupVII  3667528      14 -1.111111e-01
#> 71    pop      OPL   groupVII  3667528      14  1.000000e+00
#> 72    pop      SMR   groupVII  3667528      14 -2.000000e-01
#> 73    pop      UPD   groupVII  3667528      14 -2.600000e-01
#> 74    pop      ASP   groupVII  4189120      15 -2.500000e-01
#> 75    pop      OPL   groupVII  4189120      15  8.800000e-02
#> 76    pop      PAL   groupVII  4189120      15 -1.111111e-01
#> 77    pop      SMR   groupVII  4189120      15 -2.500000e-01
#> 78    pop      UPD   groupVII  4189120      15 -5.189873e-01
#> 79    pop      ASP    groupXX  4538074      16  2.500000e-01
#> 80    pop      CLF    groupXX  4538074      16 -3.333333e-01
#> 81    pop      OPL    groupXX  4538074      16 -2.294118e-01
#> 82    pop      PAL    groupXX  4538074      16 -1.886792e-01
#> 83    pop      SMR    groupXX  4538074      16 -2.294118e-01
#> 84    pop      UPD    groupXX  4538074      16  2.826087e-01
#> 85    pop      ASP  groupXIII  4706089      17 -1.250000e-01
#> 86    pop      OPL  groupXIII  4706089      17 -1.428571e-01
#> 87    pop      PAL  groupXIII  4706089      17  2.105263e-01
#> 88    pop      SMR  groupXIII  4706089      17 -1.764706e-01
#> 89    pop      UPD  groupXIII  4706089      17 -1.428571e-01
#> 90    pop      ASP   groupIII  4746808      18 -3.846154e-01
#> 91    pop      CLF   groupIII  4746808      18 -6.000000e-01
#> 92    pop      OPL   groupIII  4746808      18 -8.620690e-03
#> 93    pop      PAL   groupIII  4746808      18 -1.298701e-02
#> 94    pop      SMR   groupIII  4746808      18  1.791045e-01
#> 95    pop      UPD   groupIII  4746808      18  3.518519e-01
#> 96    pop      CLF  groupXVII  5041995      19 -2.500000e-01
#> 97    pop      OPL  groupXVII  5041995      19 -5.555556e-02
#> 98    pop      SMR  groupXVII  5041995      19 -4.761905e-02
#> 99    pop      UPD  groupXVII  5041995      19 -3.703704e-02
#> 100   pop      ASP     groupX  5161637      20  2.220446e-16
#> 101   pop      CLF     groupX  5161637      20  3.330669e-16
#> 102   pop      OPL     groupX  5161637      20  7.771561e-16
#> 103   pop      SMR     groupX  5161637      20 -4.761905e-02
#> 104   pop      UPD     groupX  5161637      20 -4.000000e-02
#> 105   pop      ASP   groupXII  5178037      21 -1.428571e-01
#> 106   pop      CLF   groupXII  5178037      21 -5.000000e-01
#> 107   pop      OPL   groupXII  5178037      21 -3.448276e-01
#> 108   pop      SMR   groupXII  5178037      21 -4.615385e-01
#> 109   pop      UPD   groupXII  5178037      21 -8.196721e-02
#> 110   pop      ASP    groupIX  5287711      22  1.578947e-01
#> 111   pop      CLF    groupIX  5287711      22  2.727273e-01
#> 112   pop      OPL    groupIX  5287711      22  2.173913e-02
#> 113   pop      SMR    groupIX  5287711      22 -1.515152e-01
#> 114   pop      UPD    groupIX  5287711      22 -2.093023e-01
#> 115   pop      ASP  groupVIII  5462135      23  5.551115e-16
#> 116   pop      OPL  groupVIII  5462135      23 -2.127660e-02
#> 117   pop      SMR  groupVIII  5462135      23 -8.571429e-02
#> 118   pop      UPD  groupVIII  5462135      23 -2.220446e-16
#> 119   pop      ASP    groupIX  5960933      24  5.551115e-16
#> 120   pop      OPL    groupIX  5960933      24  3.534483e-01
#> 121   pop      SMR    groupIX  5960933      24 -2.857143e-02
#> 122   pop      UPD    groupIX  5960933      24 -8.333333e-02
#> 123   pop      ASP  groupXIII  6076868      25 -3.846154e-01
#> 124   pop      CLF  groupXIII  6076868      25 -2.220446e-16
#> 125   pop      OPL  groupXIII  6076868      25 -3.333333e-01
#> 126   pop      SMR  groupXIII  6076868      25 -1.314286e-01
#> 127   pop      UPD  groupXIII  6076868      25  8.411215e-02
#> 128   pop      ASP    groupXI  6177332      26 -4.545455e-01
#> 129   pop      OPL    groupXI  6177332      26 -1.428571e-01
#> 130   pop      PAL    groupXI  6177332      26 -3.703704e-02
#> 131   pop      SMR    groupXI  6177332      26 -4.000000e-01
#> 132   pop      UPD    groupXI  6177332      26 -3.750000e-01
#> 133   pop      ASP    groupIV  6291287      27  5.551115e-16
#> 134   pop      CLF    groupIV  6291287      27  3.330669e-16
#> 135   pop      OPL    groupIV  6291287      27 -2.432432e-01
#> 136   pop      SMR    groupIV  6291287      27 -1.515152e-01
#> 137   pop      UPD    groupIV  6291287      27 -4.000000e-02
#> 138   pop      ASP    groupIX  6515796      28 -7.692308e-02
#> 139   pop      CLF    groupIX  6515796      28 -1.111111e-01
#> 140   pop      OPL    groupIX  6515796      28 -2.531646e-02
#> 141   pop      PAL    groupIX  6515796      28 -7.692308e-02
#> 142   pop      SMR    groupIX  6515796      28 -7.731959e-02
#> 143   pop      UPD    groupIX  6515796      28 -3.684211e-01
#> 144   pop      ASP   groupXII  6718079      29  2.220446e-16
#> 145   pop      CLF   groupXII  6718079      29  0.000000e+00
#> 146   pop      OPL   groupXII  6718079      29  1.616766e-01
#> 147   pop      PAL   groupXII  6718079      29 -1.666667e-01
#> 148   pop      SMR   groupXII  6718079      29  1.025641e-01
#> 149   pop      UPD   groupXII  6718079      29 -9.090909e-02
#> 150   pop      ASP   groupVII  6967477      30 -5.263158e-02
#> 151   pop      CLF   groupVII  6967477      30 -1.428571e-01
#> 152   pop      OPL   groupVII  6967477      30 -5.263158e-02
#> 153   pop      SMR   groupVII  6967477      30 -1.052632e-01
#> 154   pop      UPD   groupVII  6967477      30 -1.578947e-01
#> 155   pop      ASP   groupIII  7194069      31 -3.333333e-01
#> 156   pop      CLF   groupIII  7194069      31  3.330669e-16
#> 157   pop      OPL   groupIII  7194069      31 -8.270677e-02
#> 158   pop      SMR   groupIII  7194069      31 -1.042945e-01
#> 159   pop      UPD   groupIII  7194069      31  8.411215e-02
#> 160   pop      ASP   groupXVI  7199722      32 -1.250000e-01
#> 161   pop      CLF   groupXVI  7199722      32 -1.428571e-01
#> 162   pop      OPL   groupXVI  7199722      32 -1.282051e-01
#> 163   pop      PAL   groupXVI  7199722      32 -6.470588e-01
#> 164   pop      SMR   groupXVI  7199722      32 -1.875000e-01
#> 165   pop      UPD   groupXVI  7199722      32 -3.333333e-01
#> 166   pop      ASP  groupVIII  7224065      33 -3.333333e-01
#> 167   pop      CLF  groupVIII  7224065      33 -5.000000e-01
#> 168   pop      OPL  groupVIII  7224065      33  1.219512e-01
#> 169   pop      PAL  groupVIII  7224065      33 -2.727273e-01
#> 170   pop      SMR  groupVIII  7224065      33  3.130841e-01
#> 171   pop      UPD  groupVIII  7224065      33  5.200000e-01
#> 172   pop      ASP groupXVIII  7229426      34 -6.000000e-01
#> 173   pop      CLF groupXVIII  7229426      34  4.117647e-01
#> 174   pop      OPL groupXVIII  7229426      34 -1.409396e-01
#> 175   pop      SMR groupXVIII  7229426      34  2.575758e-01
#> 176   pop      UPD groupXVIII  7229426      34 -1.111111e-01
#> 177   pop      ASP groupXVIII  8087022      35 -5.882353e-02
#> 178   pop      OPL groupXVIII  8087022      35 -5.000000e-02
#> 179   pop      SMR groupXVIII  8087022      35 -5.555556e-02
#> 180   pop      UPD groupXVIII  8087022      35 -4.000000e-02
#> 181   pop      ASP    groupIX  8373442      36 -1.250000e-01
#> 182   pop      CLF    groupIX  8373442      36 -1.111111e-01
#> 183   pop      OPL    groupIX  8373442      36  6.395349e-02
#> 184   pop      SMR    groupIX  8373442      36 -2.000000e-01
#> 185   pop      UPD    groupIX  8373442      36 -4.347826e-02
#> 186   pop      ASP    groupXV  8815538      37  4.193548e-01
#> 187   pop      CLF    groupXV  8815538      37 -1.428571e-01
#> 188   pop      OPL    groupXV  8815538      37 -3.529412e-01
#> 189   pop      PAL    groupXV  8815538      37  0.000000e+00
#> 190   pop      SMR    groupXV  8815538      37  4.109589e-02
#> 191   pop      UPD    groupXV  8815538      37 -3.333333e-01
#> 192   pop      ASP    groupXI  9019877      38  4.375000e-01
#> 193   pop      OPL    groupXI  9019877      38 -2.222222e-01
#> 194   pop      PAL    groupXI  9019877      38 -3.983051e-01
#> 195   pop      SMR    groupXI  9019877      38 -2.705882e-01
#> 196   pop      UPD    groupXI  9019877      38  4.583333e-01
#> 197   pop      ASP   groupXXI  9067876      39 -1.428571e-01
#> 198   pop      CLF   groupXXI  9067876      39 -3.333333e-01
#> 199   pop      OPL   groupXXI  9067876      39 -7.692308e-02
#> 200   pop      PAL   groupXXI  9067876      39 -1.111111e-01
#> 201   pop      SMR   groupXXI  9067876      39 -6.666667e-02
#> 202   pop      UPD   groupXXI  9067876      39 -1.250000e-01
#> 203   pop      ASP  groupXVII  9108184      40  4.255319e-02
#> 204   pop      CLF  groupXVII  9108184      40 -1.428571e-01
#> 205   pop      OPL  groupXVII  9108184      40 -1.875000e-01
#> 206   pop      PAL  groupXVII  9108184      40 -1.666667e-01
#> 207   pop      SMR  groupXVII  9108184      40 -1.612903e-01
#> 208   pop      UPD  groupXVII  9108184      40  1.071429e-01
#> 209   pop      ASP    groupIX  9284427      41 -2.727273e-01
#> 210   pop      CLF    groupIX  9284427      41 -3.333333e-01
#> 211   pop      OPL    groupIX  9284427      41 -3.333333e-01
#> 212   pop      SMR    groupIX  9284427      41 -1.803279e-01
#> 213   pop      UPD    groupIX  9284427      41  2.173913e-01
#> 214   pop      ASP     groupV  9767769      42 -2.727273e-01
#> 215   pop      CLF     groupV  9767769      42  2.727273e-01
#> 216   pop      OPL     groupV  9767769      42  5.836576e-02
#> 217   pop      PAL     groupV  9767769      42 -3.448276e-02
#> 218   pop      SMR     groupV  9767769      42 -3.695652e-01
#> 219   pop      UPD     groupV  9767769      42 -2.093023e-01
#> 220   pop      ASP   groupXII  9881312      43 -5.263158e-02
#> 221   pop      CLF   groupXII  9881312      43 -1.428571e-01
#> 222   pop      OPL   groupXII  9881312      43 -2.121212e-01
#> 223   pop      SMR   groupXII  9881312      43 -9.523810e-02
#> 224   pop      UPD   groupXII  9881312      43 -2.631579e-01
#> 225   pop      ASP    groupXX 10296516      44  6.363636e-01
#> 226   pop      OPL    groupXX 10296516      44 -4.166667e-02
#> 227   pop      PAL    groupXX 10296516      44 -7.692308e-02
#> 228   pop      SMR    groupXX 10296516      44  3.406593e-01
#> 229   pop      UPD    groupXX 10296516      44 -2.222222e-01
#> 230   pop      UPD     groupI 10593829      45 -3.703704e-02
#> 231   pop      ASP   groupXIX 10765968      46  7.741935e-01
#> 232   pop      CLF   groupXIX 10765968      46 -1.428571e-01
#> 233   pop      OPL   groupXIX 10765968      46  0.000000e+00
#> 234   pop      SMR   groupXIX 10765968      46  4.838710e-01
#> 235   pop      UPD   groupXIX 10765968      46  1.891892e-01
#> 236   pop      ASP    groupII 10829169      47 -2.500000e-01
#> 237   pop      OPL    groupII 10829169      47 -1.984436e-01
#> 238   pop      PAL    groupII 10829169      47 -7.692308e-02
#> 239   pop      SMR    groupII 10829169      47  1.851852e-01
#> 240   pop      UPD    groupII 10829169      47  3.114754e-01
#> 241   pop      ASP    groupXX 10913115      48  2.258065e-01
#> 242   pop      CLF    groupXX 10913115      48 -2.500000e-01
#> 243   pop      OPL    groupXX 10913115      48 -1.000000e-01
#> 244   pop      PAL    groupXX 10913115      48  2.857143e-01
#> 245   pop      SMR    groupXX 10913115      48 -9.677419e-02
#> 246   pop      UPD    groupXX 10913115      48 -2.631579e-01
#> 247   pop      CLF groupXVIII 10974594      49 -1.111111e-01
#> 248   pop      OPL groupXVIII 10974594      49 -6.661338e-16
#> 249   pop      UPD groupXVIII 10974594      49 -4.347826e-02
#> 250   pop      ASP  groupXVII 11256661      50 -1.250000e-01
#> 251   pop      CLF  groupXVII 11256661      50 -2.220446e-16
#> 252   pop      OPL  groupXVII 11256661      50 -1.000000e-01
#> 253   pop      PAL  groupXVII 11256661      50  0.000000e+00
#> 254   pop      SMR  groupXVII 11256661      50 -2.857143e-02
#> 255   pop      UPD  groupXVII 11256661      50 -4.000000e-02
#> 256   pop      ASP  groupVIII 11438870      51 -1.111111e-01
#> 257   pop      OPL  groupVIII 11438870      51 -7.344633e-02
#> 258   pop      PAL  groupVIII 11438870      51 -1.764706e-01
#> 259   pop      SMR  groupVIII 11438870      51 -1.097561e-01
#> 260   pop      UPD  groupVIII 11438870      51 -1.111111e-02
#> 261   pop      ASP   groupVII 11477223      52 -5.384615e-01
#> 262   pop      CLF   groupVII 11477223      52  3.333333e-01
#> 263   pop      OPL   groupVII 11477223      52 -1.040000e-01
#> 264   pop      PAL   groupVII 11477223      52 -4.776119e-01
#> 265   pop      SMR   groupVII 11477223      52  1.843972e-01
#> 266   pop      UPD   groupVII 11477223      52  1.034483e-01
#> 267   pop      ASP    groupII 11505084      53  4.193548e-01
#> 268   pop      CLF    groupII 11505084      53 -5.000000e-01
#> 269   pop      OPL    groupII 11505084      53 -2.021858e-01
#> 270   pop      SMR    groupII 11505084      53  1.140940e-01
#> 271   pop      UPD    groupII 11505084      53  1.492537e-02
#> 272   pop      CLF    groupIX 11904651      54  3.330669e-16
#> 273   pop      OPL    groupIX 11904651      54  1.110223e-15
#> 274   pop      SMR    groupIX 11904651      54 -2.220446e-16
#> 275   pop      ASP  groupVIII 12108945      55 -4.482759e-01
#> 276   pop      CLF  groupVIII 12108945      55 -2.500000e-01
#> 277   pop      OPL  groupVIII 12108945      55  1.666667e-01
#> 278   pop      PAL  groupVIII 12108945      55 -3.084112e-01
#> 279   pop      SMR  groupVIII 12108945      55  1.592920e-01
#> 280   pop      UPD  groupVIII 12108945      55  1.111111e-01
#> 281   pop      ASP   groupXIX 12485182      56 -4.545455e-01
#> 282   pop      CLF   groupXIX 12485182      56  3.330669e-16
#> 283   pop      OPL   groupXIX 12485182      56 -2.258065e-01
#> 284   pop      PAL   groupXIX 12485182      56 -5.000000e-01
#> 285   pop      SMR   groupXIX 12485182      56 -2.000000e-01
#> 286   pop      UPD   groupXIX 12485182      56 -2.727273e-01
#> 287   pop      ASP    groupXI 12833733      57  0.000000e+00
#> 288   pop      OPL    groupXI 12833733      57 -8.571429e-02
#> 289   pop      SMR    groupXI 12833733      57  2.220446e-16
#> 290   pop      ASP groupXVIII 12937786      58  1.578947e-01
#> 291   pop      OPL groupXVIII 12937786      58 -3.333333e-01
#> 292   pop      PAL groupXVIII 12937786      58  2.500000e-01
#> 293   pop      SMR groupXVIII 12937786      58 -1.333333e-01
#> 294   pop      UPD groupXVIII 12937786      58  1.304348e-01
#> 295   pop      ASP   groupIII 13177027      59  3.103448e-01
#> 296   pop      CLF   groupIII 13177027      59 -4.285714e-01
#> 297   pop      OPL   groupIII 13177027      59 -2.873134e-01
#> 298   pop      PAL   groupIII 13177027      59 -3.448276e-02
#> 299   pop      SMR   groupIII 13177027      59 -2.797927e-01
#> 300   pop      UPD   groupIII 13177027      59  6.493506e-02
#> 301   pop      ASP    groupVI 13335577      60 -9.090909e-02
#> 302   pop      CLF    groupVI 13335577      60 -2.500000e-01
#> 303   pop      OPL    groupVI 13335577      60 -3.218884e-01
#> 304   pop      SMR    groupVI 13335577      60 -4.712042e-02
#> 305   pop      UPD    groupVI 13335577      60 -4.444444e-01
#> 306   pop      CLF     groupX 13718292      61 -2.220446e-16
#> 307   pop      OPL     groupX 13718292      61  1.110223e-16
#> 308   pop      SMR     groupX 13718292      61 -2.702703e-02
#> 309   pop      UPD     groupX 13718292      61 -4.761905e-02
#> 310   pop      ASP    groupVI 13895716      62 -5.882353e-02
#> 311   pop      OPL    groupVI 13895716      62 -2.500000e-01
#> 312   pop      SMR    groupVI 13895716      62 -4.545455e-02
#> 313   pop      UPD    groupVI 13895716      62 -8.333333e-02
#> 314   pop      ASP   groupXIV 14009990      63  0.000000e+00
#> 315   pop      CLF   groupXIV 14009990      63  2.727273e-01
#> 316   pop      OPL   groupXIV 14009990      63 -2.592593e-01
#> 317   pop      PAL   groupXIV 14009990      63 -2.000000e-01
#> 318   pop      SMR   groupXIV 14009990      63 -1.320755e-01
#> 319   pop      UPD   groupXIV 14009990      63  6.250000e-02
#> 320   pop      ASP   groupXII 14158350      64  1.110223e-16
#> 321   pop      CLF   groupXII 14158350      64  1.000000e+00
#> 322   pop      OPL   groupXII 14158350      64 -4.285714e-01
#> 323   pop      PAL   groupXII 14158350      64  1.612903e-01
#> 324   pop      SMR   groupXII 14158350      64  4.109589e-02
#> 325   pop      UPD   groupXII 14158350      64 -2.600000e-01
#> 326   pop      ASP    groupIX 14164215      65 -1.764706e-01
#> 327   pop      CLF    groupIX 14164215      65  6.250000e-02
#> 328   pop      OPL    groupIX 14164215      65  1.219512e-01
#> 329   pop      PAL    groupIX 14164215      65 -2.752294e-02
#> 330   pop      SMR    groupIX 14164215      65  1.126761e-01
#> 331   pop      UPD    groupIX 14164215      65  2.857143e-01
#> 332   pop      ASP    groupVI 14236496      66 -5.263158e-02
#> 333   pop      CLF    groupVI 14236496      66  6.250000e-02
#> 334   pop      OPL    groupVI 14236496      66 -2.954545e-01
#> 335   pop      PAL    groupVI 14236496      66  6.250000e-02
#> 336   pop      SMR    groupVI 14236496      66  1.578947e-01
#> 337   pop      UPD    groupVI 14236496      66  2.500000e-01
#> 338   pop      PAL     groupX 14273796      67  1.176471e-02
#> 339   pop      UPD     groupX 14273796      67 -4.440892e-16
#> 340   pop      ASP  groupVIII 14445389      68  2.500000e-01
#> 341   pop      CLF  groupVIII 14445389      68 -2.500000e-01
#> 342   pop      OPL  groupVIII 14445389      68 -4.020101e-02
#> 343   pop      PAL  groupVIII 14445389      68  6.250000e-02
#> 344   pop      SMR  groupVIII 14445389      68 -7.936508e-02
#> 345   pop      UPD  groupVIII 14445389      68  1.034483e-01
#> 346   pop      CLF    groupIX 14868533      69 -2.220446e-16
#> 347   pop      OPL    groupIX 14868533      69  1.627907e-01
#> 348   pop      SMR    groupIX 14868533      69 -8.108108e-02
#> 349   pop      ASP    groupIX 15190693      70 -5.000000e-01
#> 350   pop      CLF    groupIX 15190693      70  2.727273e-01
#> 351   pop      OPL    groupIX 15190693      70 -3.333333e-01
#> 352   pop      SMR    groupIX 15190693      70 -3.125000e-01
#> 353   pop      UPD    groupIX 15190693      70 -3.750000e-01
#> 354   pop      ASP   groupXII 15306362      71 -4.440892e-16
#> 355   pop      OPL   groupXII 15306362      71 -6.661338e-16
#> 356   pop      PAL   groupXII 15306362      71 -2.000000e-01
#> 357   pop      SMR   groupXII 15306362      71 -1.219512e-01
#> 358   pop      UPD   groupXII 15306362      71 -3.703704e-02
#> 359   pop      ASP    groupXX 15460430      72  3.863636e-01
#> 360   pop      CLF    groupXX 15460430      72  6.250000e-02
#> 361   pop      OPL    groupXX 15460430      72 -2.406639e-01
#> 362   pop      PAL    groupXX 15460430      72  2.970297e-02
#> 363   pop      SMR    groupXX 15460430      72 -2.727273e-01
#> 364   pop      UPD    groupXX 15460430      72 -3.333333e-01
#> 365   pop      ASP   groupVII 15493025      73 -3.846154e-01
#> 366   pop      CLF   groupVII 15493025      73 -2.500000e-01
#> 367   pop      OPL   groupVII 15493025      73 -1.282051e-01
#> 368   pop      PAL   groupVII 15493025      73 -7.692308e-02
#> 369   pop      SMR   groupVII 15493025      73 -1.111111e-01
#> 370   pop      UPD   groupVII 15493025      73  8.450704e-02
#> 371   pop      ASP    groupXV 15547073      74 -1.250000e-01
#> 372   pop      CLF    groupXV 15547073      74 -2.220446e-16
#> 373   pop      OPL    groupXV 15547073      74 -4.400000e-01
#> 374   pop      SMR    groupXV 15547073      74 -2.800000e-01
#> 375   pop      UPD    groupXV 15547073      74 -9.090909e-02
#> 376   pop      ASP  groupXIII 15697285      75 -2.972973e-01
#> 377   pop      CLF  groupXIII 15697285      75 -9.090909e-02
#> 378   pop      OPL  groupXIII 15697285      75 -2.173913e-01
#> 379   pop      PAL  groupXIII 15697285      75 -4.285714e-01
#> 380   pop      SMR  groupXIII 15697285      75 -7.731959e-02
#> 381   pop      UPD  groupXIII 15697285      75 -4.000000e-01
#> 382   pop      ASP   groupIII 15983994      76  6.250000e-02
#> 383   pop      CLF   groupIII 15983994      76 -2.220446e-16
#> 384   pop      OPL   groupIII 15983994      76  6.052632e-01
#> 385   pop      SMR   groupIII 15983994      76 -1.000000e-01
#> 386   pop      UPD   groupIII 15983994      76  1.612903e-01
#> 387   pop      ASP   groupXIX 16233699      77 -4.545455e-01
#> 388   pop      CLF   groupXIX 16233699      77  3.330669e-16
#> 389   pop      OPL   groupXIX 16233699      77 -1.891892e-01
#> 390   pop      PAL   groupXIX 16233699      77 -5.555556e-01
#> 391   pop      SMR   groupXIX 16233699      77 -2.000000e-01
#> 392   pop      UPD   groupXIX 16233699      77 -2.380952e-01
#> 393   pop      CLF    groupIX 16306966      78  6.153846e-01
#> 394   pop      OPL    groupIX 16306966      78  4.060914e-02
#> 395   pop      SMR    groupIX 16306966      78 -1.052632e-01
#> 396   pop      ASP     groupI 16498302      79 -1.764706e-01
#> 397   pop      OPL     groupI 16498302      79 -1.111111e-01
#> 398   pop      PAL     groupI 16498302      79 -2.800000e-01
#> 399   pop      SMR     groupI 16498302      79  4.683544e-01
#> 400   pop      UPD     groupI 16498302      79 -1.000000e-01
#> 401   pop      ASP  groupXIII 16520929      80  2.682927e-01
#> 402   pop      CLF  groupXIII 16520929      80 -1.111111e-01
#> 403   pop      OPL  groupXIII 16520929      80 -2.352941e-01
#> 404   pop      PAL  groupXIII 16520929      80 -4.392523e-01
#> 405   pop      SMR  groupXIII 16520929      80 -3.333333e-01
#> 406   pop      UPD  groupXIII 16520929      80 -3.333333e-01
#> 407   pop      ASP   groupXVI 16561891      81 -4.545455e-01
#> 408   pop      CLF   groupXVI 16561891      81 -6.000000e-01
#> 409   pop      OPL   groupXVI 16561891      81 -2.621359e-01
#> 410   pop      PAL   groupXVI 16561891      81 -2.396694e-01
#> 411   pop      SMR   groupXVI 16561891      81 -2.352941e-01
#> 412   pop      UPD   groupXVI 16561891      81 -3.146067e-01
#> 413   pop      ASP    groupII 16598364      82 -5.882353e-02
#> 414   pop      OPL    groupII 16598364      82 -8.695652e-02
#> 415   pop      SMR    groupII 16598364      82 -2.857143e-02
#> 416   pop      UPD    groupII 16598364      82 -2.000000e-01
#> 417   pop      ASP    groupIV 17248473      83 -1.764706e-01
#> 418   pop      CLF    groupIV 17248473      83 -6.666667e-01
#> 419   pop      SMR    groupIV 17248473      83  7.771561e-16
#> 420   pop      UPD    groupIV 17248473      83 -2.220446e-16
#> 421   pop      ASP    groupXX 17748944      84  1.000000e-01
#> 422   pop      CLF    groupXX 17748944      84  3.333333e-01
#> 423   pop      OPL    groupXX 17748944      84 -7.344633e-02
#> 424   pop      SMR    groupXX 17748944      84 -5.882353e-02
#> 425   pop      UPD    groupXX 17748944      84 -2.000000e-01
#> 426   pop      ASP  groupXIII 17842228      85 -4.651163e-02
#> 427   pop      CLF  groupXIII 17842228      85 -6.000000e-01
#> 428   pop      OPL  groupXIII 17842228      85 -4.285714e-01
#> 429   pop      PAL  groupXIII 17842228      85 -8.333333e-02
#> 430   pop      SMR  groupXIII 17842228      85 -8.433735e-02
#> 431   pop      UPD  groupXIII 17842228      85 -1.886792e-01
#> 432   pop      ASP   groupXIX 17858769      86 -5.555556e-01
#> 433   pop      CLF   groupXIX 17858769      86  3.330669e-16
#> 434   pop      OPL   groupXIX 17858769      86 -2.666667e-01
#> 435   pop      PAL   groupXIX 17858769      86 -6.666667e-01
#> 436   pop      SMR   groupXIX 17858769      86 -2.000000e-01
#> 437   pop      UPD   groupXIX 17858769      86 -3.333333e-01
#> 438   pop      ASP    groupIX 17929225      87  1.379310e-01
#> 439   pop      CLF    groupIX 17929225      87  3.330669e-16
#> 440   pop      OPL    groupIX 17929225      87  1.428571e-01
#> 441   pop      SMR    groupIX 17929225      87 -3.548387e-01
#> 442   pop      UPD    groupIX 17929225      87 -1.304348e-01
#> 443   pop      ASP    groupIX 18339767      88 -4.482759e-01
#> 444   pop      CLF    groupIX 18339767      88 -2.500000e-01
#> 445   pop      OPL    groupIX 18339767      88  1.110223e-16
#> 446   pop      SMR    groupIX 18339767      88 -3.983051e-01
#> 447   pop      UPD    groupIX 18339767      88 -4.666667e-01
#> 448   pop      ASP  groupXIII 18836558      89 -4.545455e-01
#> 449   pop      CLF  groupXIII 18836558      89  1.000000e+00
#> 450   pop      OPL  groupXIII 18836558      89 -5.517241e-02
#> 451   pop      PAL  groupXIII 18836558      89  6.901408e-01
#> 452   pop      SMR  groupXIII 18836558      89 -5.517241e-02
#> 453   pop      UPD  groupXIII 18836558      89  3.529412e-01
#> 454   pop      ASP    groupII 19244005      90 -1.111111e-01
#> 455   pop      SMR    groupII 19244005      90 -2.220446e-16
#> 456   pop      UPD    groupII 19244005      90 -4.347826e-02
#> 457   pop      ASP    groupXX 19317273      91 -1.428571e-01
#> 458   pop      CLF    groupXX 19317273      91  5.714286e-01
#> 459   pop      OPL    groupXX 19317273      91 -2.149533e-01
#> 460   pop      PAL    groupXX 19317273      91 -7.142857e-02
#> 461   pop      SMR    groupXX 19317273      91  1.554404e-02
#> 462   pop      UPD    groupXX 19317273      91  1.034483e-01
#> 463   pop      ASP     groupI 19759331      92  2.173913e-01
#> 464   pop      CLF     groupI 19759331      92 -1.428571e-01
#> 465   pop      OPL     groupI 19759331      92 -5.294118e-01
#> 466   pop      PAL     groupI 19759331      92 -1.320755e-01
#> 467   pop      SMR     groupI 19759331      92 -1.764706e-01
#> 468   pop      UPD     groupI 19759331      92 -2.000000e-01
#> 469   pop      ASP   groupVII 20165240      93 -1.111111e-01
#> 470   pop      CLF   groupVII 20165240      93 -2.500000e-01
#> 471   pop      OPL   groupVII 20165240      93 -4.166667e-02
#> 472   pop      PAL   groupVII 20165240      93 -1.304348e-01
#> 473   pop      SMR   groupVII 20165240      93 -2.000000e-01
#> 474   pop      UPD   groupVII 20165240      93 -3.000000e-01
#> 475   pop      OPL    groupII 20496069      94 -5.000000e-02
#> 476   pop      PAL    groupII 20496069      94 -1.506849e-01
#> 477   pop      SMR    groupII 20496069      94 -4.440892e-16
#> 478   pop      UPD    groupII 20496069      94 -8.333333e-02
#> 479   pop      ASP     groupI 21284266      95 -2.000000e-01
#> 480   pop      CLF     groupI 21284266      95  3.333333e-01
#> 481   pop      OPL     groupI 21284266      95  4.084507e-01
#> 482   pop      SMR     groupI 21284266      95 -1.612903e-01
#> 483   pop      UPD     groupI 21284266      95 -1.578947e-01
#> 484   pop      ASP    groupIV 22227311      96  0.000000e+00
#> 485   pop      CLF    groupIV 22227311      96  3.330669e-16
#> 486   pop      OPL    groupIV 22227311      96 -1.891892e-01
#> 487   pop      SMR    groupIV 22227311      96 -7.317073e-02
#> 488   pop      UPD    groupIV 22227311      96 -4.347826e-02
#> 489   pop      CLF     groupI 22746027      97  3.330669e-16
#> 490   pop      OPL     groupI 22746027      97 -2.564103e-02
#> 491   pop      SMR     groupI 22746027      97 -2.702703e-02
#> 492   pop      UPD     groupI 22746027      97 -2.220446e-16
#> 493   pop      ASP   groupVII 25239495      98  4.074074e-01
#> 494   pop      CLF   groupVII 25239495      98  6.000000e-01
#> 495   pop      OPL   groupVII 25239495      98 -8.695652e-02
#> 496   pop      PAL   groupVII 25239495      98 -1.344538e-01
#> 497   pop      SMR   groupVII 25239495      98 -2.352941e-01
#> 498   pop      UPD   groupVII 25239495      98 -1.506849e-01
#> 499   pop      ASP   groupVII 27599144      99  4.255319e-02
#> 500   pop      CLF   groupVII 27599144      99  7.058824e-01
#> 501   pop      OPL   groupVII 27599144      99 -3.397129e-01
#> 502   pop      SMR   groupVII 27599144      99  1.230769e-01
#> 503   pop      UPD   groupVII 27599144      99  1.379310e-01
#> 504   pop      ASP    groupIV 32321895     100  2.173913e-01
#> 505   pop      CLF    groupIV 32321895     100 -5.000000e-01
#> 506   pop      OPL    groupIV 32321895     100  7.317073e-02
#> 507   pop      PAL    groupIV 32321895     100  1.110223e-16
#> 508   pop      SMR    groupIV 32321895     100 -4.615385e-01
#> 509   pop      UPD    groupIV 32321895     100 -1.111111e-02
#> 
#> $weighted.means
#>     facet subfacet snp.facet snp.subfacet weighted_mean_fis
#> 1     pop      ASP     .base        .base     -6.617460e-02
#> 2     pop      ASP       chr       groupI     -4.100946e-02
#> 3     pop      ASP       chr      groupII      2.770083e-03
#> 4     pop      ASP       chr     groupIII     -4.096386e-02
#> 5     pop      ASP       chr      groupIV      4.862579e-02
#> 6     pop      ASP       chr      groupIX     -1.499760e-01
#> 7     pop      ASP       chr       groupV     -2.799097e-01
#> 8     pop      ASP       chr      groupVI      8.996283e-02
#> 9     pop      ASP       chr     groupVII     -1.641874e-01
#> 10    pop      ASP       chr    groupVIII     -1.574803e-01
#> 11    pop      ASP       chr       groupX      2.220446e-16
#> 12    pop      ASP       chr      groupXI      4.382470e-02
#> 13    pop      ASP       chr     groupXII      4.000000e-02
#> 14    pop      ASP       chr    groupXIII     -1.810605e-01
#> 15    pop      ASP       chr     groupXIV      0.000000e+00
#> 16    pop      ASP       chr     groupXIX     -1.750534e-01
#> 17    pop      ASP       chr      groupXV      1.262136e-01
#> 18    pop      ASP       chr     groupXVI     -2.692308e-01
#> 19    pop      ASP       chr    groupXVII     -2.156640e-02
#> 20    pop      ASP       chr   groupXVIII     -1.653179e-01
#> 21    pop      ASP       chr      groupXX      2.472325e-01
#> 22    pop      ASP       chr     groupXXI     -1.428571e-01
#> 23    pop      CLF     .base        .base     -1.780634e-02
#> 24    pop      CLF       chr       groupI      1.304348e-01
#> 25    pop      CLF       chr      groupII     -5.000000e-01
#> 26    pop      CLF       chr     groupIII     -3.333333e-01
#> 27    pop      CLF       chr      groupIV     -4.285714e-01
#> 28    pop      CLF       chr      groupIX      1.446029e-01
#> 29    pop      CLF       chr       groupV      2.912621e-02
#> 30    pop      CLF       chr      groupVI     -1.464968e-01
#> 31    pop      CLF       chr     groupVII      1.884058e-01
#> 32    pop      CLF       chr    groupVIII     -3.235294e-01
#> 33    pop      CLF       chr       groupX      0.000000e+00
#> 34    pop      CLF       chr     groupXII      6.432749e-02
#> 35    pop      CLF       chr    groupXIII      1.754386e-02
#> 36    pop      CLF       chr     groupXIV      2.727273e-01
#> 37    pop      CLF       chr     groupXIX     -4.347826e-02
#> 38    pop      CLF       chr      groupXV     -6.666667e-02
#> 39    pop      CLF       chr     groupXVI     -2.080537e-01
#> 40    pop      CLF       chr    groupXVII     -1.382114e-01
#> 41    pop      CLF       chr   groupXVIII      2.307692e-01
#> 42    pop      CLF       chr      groupXX      7.032181e-02
#> 43    pop      CLF       chr     groupXXI     -3.333333e-01
#> 44    pop      OPL     .base        .base     -1.204339e-01
#> 45    pop      OPL       chr       groupI     -6.535640e-02
#> 46    pop      OPL       chr      groupII     -1.665328e-01
#> 47    pop      OPL       chr     groupIII     -7.890083e-03
#> 48    pop      OPL       chr      groupIV     -1.110005e-01
#> 49    pop      OPL       chr      groupIX      6.736488e-03
#> 50    pop      OPL       chr       groupV      3.858160e-02
#> 51    pop      OPL       chr      groupVI     -2.626423e-01
#> 52    pop      OPL       chr     groupVII     -3.419104e-02
#> 53    pop      OPL       chr    groupVIII      4.237882e-02
#> 54    pop      OPL       chr       groupX      4.440892e-16
#> 55    pop      OPL       chr      groupXI     -1.657982e-01
#> 56    pop      OPL       chr     groupXII     -1.959288e-01
#> 57    pop      OPL       chr    groupXIII     -2.741218e-01
#> 58    pop      OPL       chr     groupXIV     -2.592593e-01
#> 59    pop      OPL       chr     groupXIX     -1.673077e-01
#> 60    pop      OPL       chr      groupXV     -3.368966e-01
#> 61    pop      OPL       chr     groupXVI     -8.715027e-02
#> 62    pop      OPL       chr    groupXVII     -1.202680e-01
#> 63    pop      OPL       chr   groupXVIII     -1.888022e-01
#> 64    pop      OPL       chr      groupXX     -1.658705e-01
#> 65    pop      OPL       chr     groupXXI     -7.692308e-02
#> 66    pop      PAL     .base        .base     -1.496183e-01
#> 67    pop      PAL       chr       groupI     -2.121212e-01
#> 68    pop      PAL       chr      groupII     -1.275168e-01
#> 69    pop      PAL       chr     groupIII     -1.827676e-02
#> 70    pop      PAL       chr      groupIV      0.000000e+00
#> 71    pop      PAL       chr      groupIX     -4.186047e-02
#> 72    pop      PAL       chr       groupV     -3.265306e-01
#> 73    pop      PAL       chr      groupVI      1.777116e-02
#> 74    pop      PAL       chr     groupVII     -1.929419e-01
#> 75    pop      PAL       chr    groupVIII     -1.885067e-01
#> 76    pop      PAL       chr       groupX      1.176471e-02
#> 77    pop      PAL       chr      groupXI     -2.194377e-01
#> 78    pop      PAL       chr     groupXII     -3.213611e-02
#> 79    pop      PAL       chr    groupXIII     -3.686755e-02
#> 80    pop      PAL       chr     groupXIV     -2.000000e-01
#> 81    pop      PAL       chr     groupXIX     -5.295790e-01
#> 82    pop      PAL       chr      groupXV      0.000000e+00
#> 83    pop      PAL       chr     groupXVI     -4.065230e-01
#> 84    pop      PAL       chr    groupXVII     -1.351351e-01
#> 85    pop      PAL       chr   groupXVIII      2.500000e-01
#> 86    pop      PAL       chr      groupXX      4.462294e-04
#> 87    pop      PAL       chr     groupXXI     -1.111111e-01
#> 88    pop      SMR     .base        .base     -9.068452e-02
#> 89    pop      SMR       chr       groupI      1.932901e-02
#> 90    pop      SMR       chr      groupII      1.195396e-01
#> 91    pop      SMR       chr     groupIII     -1.124692e-01
#> 92    pop      SMR       chr      groupIV     -2.743902e-01
#> 93    pop      SMR       chr      groupIX     -1.725088e-01
#> 94    pop      SMR       chr       groupV     -3.041228e-01
#> 95    pop      SMR       chr      groupVI      5.934854e-02
#> 96    pop      SMR       chr     groupVII     -8.772157e-02
#> 97    pop      SMR       chr    groupVIII      7.577214e-02
#> 98    pop      SMR       chr       groupX     -3.928394e-02
#> 99    pop      SMR       chr      groupXI     -3.177796e-01
#> 100   pop      SMR       chr     groupXII     -1.545886e-01
#> 101   pop      SMR       chr    groupXIII     -9.256572e-02
#> 102   pop      SMR       chr     groupXIV     -1.320755e-01
#> 103   pop      SMR       chr     groupXIX     -5.726872e-02
#> 104   pop      SMR       chr      groupXV     -7.459711e-02
#> 105   pop      SMR       chr     groupXVI     -1.267367e-01
#> 106   pop      SMR       chr    groupXVII     -8.995221e-02
#> 107   pop      SMR       chr   groupXVIII      8.772943e-02
#> 108   pop      SMR       chr      groupXX     -7.916039e-02
#> 109   pop      SMR       chr     groupXXI     -6.666667e-02
#> 110   pop      UPD     .base        .base     -8.369860e-02
#> 111   pop      UPD       chr       groupI     -1.250652e-01
#> 112   pop      UPD       chr      groupII      1.557546e-02
#> 113   pop      UPD       chr     groupIII      1.952233e-02
#> 114   pop      UPD       chr      groupIV     -2.016349e-02
#> 115   pop      UPD       chr      groupIX     -1.305238e-01
#> 116   pop      UPD       chr       groupV     -2.093023e-01
#> 117   pop      UPD       chr      groupVI     -3.067485e-02
#> 118   pop      UPD       chr     groupVII     -1.372773e-01
#> 119   pop      UPD       chr    groupVIII      1.671309e-01
#> 120   pop      UPD       chr       groupX     -3.473227e-02
#> 121   pop      UPD       chr      groupXI      5.921053e-02
#> 122   pop      UPD       chr     groupXII     -1.293171e-01
#> 123   pop      UPD       chr    groupXIII     -1.498276e-01
#> 124   pop      UPD       chr     groupXIV      6.250000e-02
#> 125   pop      UPD       chr     groupXIX     -1.879145e-01
#> 126   pop      UPD       chr      groupXV     -2.086331e-01
#> 127   pop      UPD       chr     groupXVI     -2.823674e-01
#> 128   pop      UPD       chr    groupXVII     -7.613569e-04
#> 129   pop      UPD       chr   groupXVIII      1.254460e-02
#> 130   pop      UPD       chr      groupXX     -1.139534e-01
#> 131   pop      UPD       chr     groupXXI     -1.250000e-01
#> 
```
