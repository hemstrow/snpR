# Run case/control or quantitative association tests on SNP data.

Runs several different association tests on SNP data. The response
variable must have only two different categories (as in case/control)
for most test types, although the "gmmat.score" method supports
quantitative traits. Tests may be broken up by sample-specific facets.

## Usage

``` r
calc_association(
  x,
  facets = NULL,
  response,
  method = "gmmat.score",
  w = c(0, 1, 2),
  formula = NULL,
  family.override = FALSE,
  maxiter = 500,
  sampleID = NULL,
  Gmaf = 0,
  par = FALSE,
  cleanup = TRUE,
  verbose = FALSE
)
```

## Arguments

- x:

  snpRdata object

- facets:

  character, default NULL. Categorical metadata variables by which to
  break up analysis. See
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md)
  for more details.

- response:

  character. Name of the column containing the response variable of
  interest. Must match a column name in sample metadata. Response must
  be categorical, with only two categories.

- method:

  character, default "gmmat.score". Specifies association method.
  Options:

  - gmmat.score: Population/family structure corrected mlm approach,
    based on Chen et al (2016).

  - armitage: Armitage association test, based on Armitage (1955).

  - odds_ratio: Log odds ratio test.

  - chisq: Chi-squared test. See description for more details.

- w:

  numeric, default c(0, 1, 2). Weight variable for each genotype for the
  Armitage association method. See description for details.

- formula:

  character, default set to response ~ 1. Null formula for the response
  variable, as described in
  [`formula`](https://rdrr.io/r/stats/formula.html).

- family.override:

  character, default NULL. Provides an alternative model family object
  to use for GMMAT GWAS regression. By default, uses
  [`gaussian`](https://rdrr.io/r/stats/family.html), link = "identity"
  for a quantitative phenotype and
  [`binomial`](https://rdrr.io/r/stats/family.html), link = "logit" for
  a categorical phenotype.

- maxiter:

  numeric, default 500. Maximum iterations to use when fitting the glmm
  when using the gmmat.score option.

- sampleID:

  character, default NULL. Optional, the name of a column in the sample
  metadata to use as a sampleID when using the gmmat.score option.

- Gmaf:

  numeric, default NULL. If using the "GMMAT" option, can provide and
  optional minor allele frequency filter used when constructing the
  relatedness matrix.

- par:

  numeric or FALSE, default FALSE. Number of parallel cores to use for
  computation.

- cleanup:

  logical, default TRUE. If TRUE, files produced by and for the "gmmat"
  option will be removed after completion.

- verbose:

  Logical, default FALSE. If TRUE, some progress updates will be printed
  to the console.

## Value

A snpRdata object with the resulting association test results merged
into the stats socket.

## Details

Several methods can be used: Armitage, chi-squared, and odds ratio. For
The Armitage approach weights should be provided to the "w" argument,
which specifies the weight for each possible genotype (homozygote 1,
heterozygote, homozygote 2). The default, c(0,1,2), specifies an
additive model. The "gmmat.score" method uses the approach described in
Chen et al. (2016) and implemented in the
[`glmmkin`](https://rdrr.io/pkg/GMMAT/man/glmmkin.html) and
[`glmm.score`](https://rdrr.io/pkg/GMMAT/man/glmm.score.html) functions.
For this method, a 'G' genetic relatedness matrix is first created using
the [`Gmatrix`](https://rdrr.io/pkg/AGHmatrix/man/Gmatrix.html) function
according to Yang et al 2010.

Multi-category data is currently not supported, but is under
development.

Facets are specified as described in
[`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).
NULL and "all" facet specifications function as described.

Note that if all individuals in one facet level have one of the two
possible phenotypes, each association test will return NA or NaN. As a
result, if the response variable overlaps perfectly with the facet
variable, all results will be NA or NaN.

If the chisq method is used, a column will also be returned that
specifies which allele (major or minor) is associated with which
phenotype.

## References

Armitage (1955). Tests for Linear Trends in Proportions and Frequencies.
*Biometrics*.

Chen et al. (2016). Control for Population Structure and Relatedness for
Binary Traits in Genetic Association Studies via Logistic Mixed Models.
*American Journal of Human Genetics*.

Yang et al. (2010). Common SNPs explain a large proportion of the
heritability for human height. *Nature Genetics*.

## Author

William Hemstrom

Keming Su

Avani Chitre

## Examples

``` r
  # add a dummy phenotype and run an association test.
  x <- stickSNPs
  sample.meta(x)$phenotype <- sample(c("A", "B"), nsamps(stickSNPs), TRUE)
#> Warning: Some levels are duplicated across multiple sample meta facets.
#> This will cause issues if those sample facets are run during analysis.
#> Issues:
#> 
#> Level: A in facets: fam, phenotype
#> Level: B in facets: fam, phenotype
  x <- calc_association(x, facets = c("pop"), response = "phenotype", method = "armitage")
  get.snpR.stats(x, "pop", "association")
#> $single
#>     facet subfacet        chr position .snp.id p_armitage_phenotype
#> 1     pop      ASP    groupVI   212436       1           1.00000000
#> 2     pop      CLF    groupVI   212436       1           0.36131043
#> 3     pop      OPL    groupVI   212436       1           0.28901174
#> 4     pop      PAL    groupVI   212436       1           0.35854174
#> 5     pop      SMR    groupVI   212436       1           1.00000000
#> 6     pop      UPD    groupVI   212436       1           0.40138717
#> 7     pop      ASP groupXVIII   350070       2           0.38942370
#> 8     pop      OPL groupXVIII   350070       2           0.29464379
#> 9     pop      SMR groupXVIII   350070       2           0.07464618
#> 10    pop      ASP    groupXV   854285       3           0.65761918
#> 11    pop      CLF    groupXV   854285       3           0.62420611
#> 12    pop      OPL    groupXV   854285       3           0.17294598
#> 13    pop      PAL    groupXV   854285       3           0.36237965
#> 14    pop      SMR    groupXV   854285       3           0.07534771
#> 15    pop      UPD    groupXV   854285       3           0.21484695
#> 16    pop      ASP  groupXIII   907279       4           0.20075068
#> 17    pop      CLF  groupXIII   907279       4           0.57615012
#> 18    pop      OPL  groupXIII   907279       4           0.41324709
#> 19    pop      SMR  groupXIII   907279       4           0.09118157
#> 20    pop      UPD  groupXIII   907279       4           0.80502027
#> 21    pop      ASP    groupVI  1068620       5           0.60840780
#> 22    pop      OPL    groupVI  1068620       5           0.15990034
#> 23    pop      PAL    groupVI  1068620       5           0.92115865
#> 24    pop      SMR    groupVI  1068620       5           0.86461805
#> 25    pop      UPD    groupVI  1068620       5           0.18861420
#> 26    pop      ASP     groupV  1392753       6           0.49015296
#> 27    pop      CLF     groupV  1392753       6           0.27332168
#> 28    pop      OPL     groupV  1392753       6           0.04669009
#> 29    pop      PAL     groupV  1392753       6           0.76490886
#> 30    pop      SMR     groupV  1392753       6           0.36806551
#> 31    pop      UPD     groupV  1392753       6           0.23845478
#> 32    pop      ASP    groupIX  2075079       7           0.42787981
#> 33    pop      CLF    groupIX  2075079       7           0.62420611
#> 34    pop      OPL    groupIX  2075079       7           0.89010618
#> 35    pop      SMR    groupIX  2075079       7           0.60384825
#> 36    pop      ASP  groupXIII  2106000       8           0.81897085
#> 37    pop      CLF  groupXIII  2106000       8           0.79048223
#> 38    pop      OPL  groupXIII  2106000       8           0.62124649
#> 39    pop      SMR  groupXIII  2106000       8           0.16245550
#> 40    pop      UPD  groupXIII  2106000       8           0.05777957
#> 41    pop      ASP   groupXII  2136587       9           0.86199249
#> 42    pop      CLF   groupXII  2136587       9           0.43857803
#> 43    pop      OPL   groupXII  2136587       9           0.42034494
#> 44    pop      PAL   groupXII  2136587       9           0.36237965
#> 45    pop      SMR   groupXII  2136587       9           0.36131043
#> 46    pop      UPD   groupXII  2136587       9           0.83960282
#> 47    pop      ASP  groupXVII  2560839      10           0.65761918
#> 48    pop      CLF  groupXVII  2560839      10           0.01430588
#> 49    pop      OPL  groupXVII  2560839      10           0.36563304
#> 50    pop      SMR  groupXVII  2560839      10           0.13801074
#> 51    pop      UPD  groupXVII  2560839      10           0.04042398
#> 52    pop      ASP   groupIII  2750695      11           0.58006941
#> 53    pop      CLF   groupIII  2750695      11           0.05777957
#> 54    pop      OPL   groupIII  2750695      11           0.55115709
#> 55    pop      SMR   groupIII  2750695      11           1.00000000
#> 56    pop      UPD   groupIII  2750695      11           0.34715253
#> 57    pop      ASP   groupXIX  3008524      12           0.19670560
#> 58    pop      CLF   groupXIX  3008524      12           0.62420611
#> 59    pop      OPL   groupXIX  3008524      12           0.78293508
#> 60    pop      PAL   groupXIX  3008524      12           0.09072987
#> 61    pop      SMR   groupXIX  3008524      12           0.02689534
#> 62    pop      UPD   groupXIX  3008524      12           0.03501498
#> 63    pop      ASP   groupXVI  3624417      13           0.57075039
#> 64    pop      CLF   groupXVI  3624417      13           0.18403863
#> 65    pop      OPL   groupXVI  3624417      13           0.75099205
#> 66    pop      PAL   groupXVI  3624417      13           0.33292161
#> 67    pop      SMR   groupXVI  3624417      13           0.02139757
#> 68    pop      UPD   groupXVI  3624417      13           0.56861821
#> 69    pop      ASP   groupVII  3667528      14           0.59816145
#> 70    pop      CLF   groupVII  3667528      14           0.12133525
#> 71    pop      OPL   groupVII  3667528      14           0.19123538
#> 72    pop      SMR   groupVII  3667528      14           0.80580254
#> 73    pop      UPD   groupVII  3667528      14           0.58135198
#> 74    pop      ASP   groupVII  4189120      15           0.13684387
#> 75    pop      OPL   groupVII  4189120      15           0.64372475
#> 76    pop      PAL   groupVII  4189120      15           0.38273309
#> 77    pop      SMR   groupVII  4189120      15           0.52826533
#> 78    pop      UPD   groupVII  4189120      15           0.69658700
#> 79    pop      ASP    groupXX  4538074      16           1.00000000
#> 80    pop      CLF    groupXX  4538074      16           0.17090352
#> 81    pop      OPL    groupXX  4538074      16           0.90621802
#> 82    pop      PAL    groupXX  4538074      16           0.86461805
#> 83    pop      SMR    groupXX  4538074      16           0.90621802
#> 84    pop      UPD    groupXX  4538074      16           0.93923546
#> 85    pop      ASP  groupXIII  4706089      17           0.88030337
#> 86    pop      OPL  groupXIII  4706089      17           0.33247354
#> 87    pop      PAL  groupXIII  4706089      17           0.17192075
#> 88    pop      SMR  groupXIII  4706089      17           0.53663178
#> 89    pop      UPD  groupXIII  4706089      17           0.56861821
#> 90    pop      ASP   groupIII  4746808      18           0.59816145
#> 91    pop      CLF   groupIII  4746808      18           0.57615012
#> 92    pop      OPL   groupIII  4746808      18           0.53993935
#> 93    pop      PAL   groupIII  4746808      18           0.21535047
#> 94    pop      SMR   groupIII  4746808      18           0.79049531
#> 95    pop      UPD   groupIII  4746808      18           0.43237390
#> 96    pop      CLF  groupXVII  5041995      19           0.27332168
#> 97    pop      OPL  groupXVII  5041995      19           0.53116784
#> 98    pop      SMR  groupXVII  5041995      19           0.10329969
#> 99    pop      UPD  groupXVII  5041995      19           0.21484695
#> 100   pop      ASP     groupX  5161637      20           0.45325470
#> 101   pop      CLF     groupX  5161637      20           0.62420611
#> 102   pop      OPL     groupX  5161637      20           0.26678531
#> 103   pop      SMR     groupX  5161637      20           0.58996173
#> 104   pop      UPD     groupX  5161637      20           0.04042398
#> 105   pop      ASP   groupXII  5178037      21           0.63525630
#> 106   pop      CLF   groupXII  5178037      21           0.50498508
#> 107   pop      OPL   groupXII  5178037      21           1.00000000
#> 108   pop      SMR   groupXII  5178037      21           0.08127681
#> 109   pop      UPD   groupXII  5178037      21           0.53133551
#> 110   pop      ASP    groupIX  5287711      22           0.74911913
#> 111   pop      CLF    groupIX  5287711      22           0.76508719
#> 112   pop      OPL    groupIX  5287711      22           0.42564705
#> 113   pop      SMR    groupIX  5287711      22           0.05096194
#> 114   pop      UPD    groupIX  5287711      22           0.68486037
#> 115   pop      ASP  groupVIII  5462135      23           0.49015296
#> 116   pop      OPL  groupVIII  5462135      23           0.95293298
#> 117   pop      SMR  groupVIII  5462135      23           0.26355248
#> 118   pop      UPD  groupVIII  5462135      23           0.39802472
#> 119   pop      ASP    groupIX  5960933      24           0.19670560
#> 120   pop      OPL    groupIX  5960933      24           0.46621194
#> 121   pop      SMR    groupIX  5960933      24           0.15608367
#> 122   pop      UPD    groupIX  5960933      24           0.92265025
#> 123   pop      ASP  groupXIII  6076868      25           0.42919530
#> 124   pop      CLF  groupXIII  6076868      25           0.57615012
#> 125   pop      OPL  groupXIII  6076868      25           0.30531926
#> 126   pop      SMR  groupXIII  6076868      25           0.47393612
#> 127   pop      UPD  groupXIII  6076868      25           0.88328589
#> 128   pop      ASP    groupXI  6177332      26           0.34278171
#> 129   pop      OPL    groupXI  6177332      26           0.67570585
#> 130   pop      PAL    groupXI  6177332      26           0.64183818
#> 131   pop      SMR    groupXI  6177332      26           0.13799635
#> 132   pop      UPD    groupXI  6177332      26           0.92115865
#> 133   pop      ASP    groupIV  6291287      27           0.19670560
#> 134   pop      CLF    groupIV  6291287      27           0.62420611
#> 135   pop      OPL    groupIV  6291287      27           0.32718688
#> 136   pop      SMR    groupIV  6291287      27           0.09543728
#> 137   pop      UPD    groupIV  6291287      27           0.82549791
#> 138   pop      ASP    groupIX  6515796      28           0.67328998
#> 139   pop      CLF    groupIX  6515796      28           0.43857803
#> 140   pop      OPL    groupIX  6515796      28           0.07229079
#> 141   pop      PAL    groupIX  6515796      28           0.60477285
#> 142   pop      SMR    groupIX  6515796      28           0.32993413
#> 143   pop      UPD    groupIX  6515796      28           0.63999401
#> 144   pop      ASP   groupXII  6718079      29           0.34278171
#> 145   pop      CLF   groupXII  6718079      29           0.04550026
#> 146   pop      OPL   groupXII  6718079      29           0.84934528
#> 147   pop      PAL   groupXII  6718079      29           0.01041499
#> 148   pop      SMR   groupXII  6718079      29           0.14677211
#> 149   pop      UPD   groupXII  6718079      29           0.41644789
#> 150   pop      ASP   groupVII  6967477      30           0.03861380
#> 151   pop      CLF   groupVII  6967477      30           0.17090352
#> 152   pop      OPL   groupVII  6967477      30           0.47553266
#> 153   pop      SMR   groupVII  6967477      30           0.61092956
#> 154   pop      UPD   groupVII  6967477      30           0.66500554
#> 155   pop      ASP   groupIII  7194069      31           0.81897085
#> 156   pop      CLF   groupIII  7194069      31           0.62420611
#> 157   pop      OPL   groupIII  7194069      31           0.18834797
#> 158   pop      SMR   groupIII  7194069      31           0.32282339
#> 159   pop      UPD   groupIII  7194069      31           0.37840966
#> 160   pop      ASP   groupXVI  7199722      32           0.77815969
#> 161   pop      CLF   groupXVI  7199722      32           0.36131043
#> 162   pop      OPL   groupXVI  7199722      32           0.90132060
#> 163   pop      PAL   groupXVI  7199722      32           0.07007554
#> 164   pop      SMR   groupXVI  7199722      32           0.63920743
#> 165   pop      UPD   groupXVI  7199722      32           0.72493893
#> 166   pop      ASP  groupVIII  7224065      33           0.30305233
#> 167   pop      CLF  groupVIII  7224065      33           0.50498508
#> 168   pop      OPL  groupVIII  7224065      33           0.47191958
#> 169   pop      PAL  groupVIII  7224065      33           0.83266211
#> 170   pop      SMR  groupVIII  7224065      33           0.16180287
#> 171   pop      UPD  groupVIII  7224065      33           0.57355880
#> 172   pop      ASP groupXVIII  7229426      34           0.39136594
#> 173   pop      CLF groupXVIII  7229426      34           0.17971249
#> 174   pop      OPL groupXVIII  7229426      34           0.09353251
#> 175   pop      SMR groupXVIII  7229426      34           0.95672093
#> 176   pop      UPD groupXVIII  7229426      34           0.62101437
#> 177   pop      ASP groupXVIII  8087022      35           0.30062299
#> 178   pop      OPL groupXVIII  8087022      35           0.65004138
#> 179   pop      SMR groupXVIII  8087022      35           0.53116784
#> 180   pop      UPD groupXVIII  8087022      35           0.18587673
#> 181   pop      ASP    groupIX  8373442      36           0.09763454
#> 182   pop      CLF    groupIX  8373442      36           0.43857803
#> 183   pop      OPL    groupIX  8373442      36           0.41354806
#> 184   pop      SMR    groupIX  8373442      36           0.57110727
#> 185   pop      UPD    groupIX  8373442      36           0.05181233
#> 186   pop      ASP    groupXV  8815538      37           0.55930495
#> 187   pop      CLF    groupXV  8815538      37           0.36131043
#> 188   pop      OPL    groupXV  8815538      37           0.39174065
#> 189   pop      PAL    groupXV  8815538      37           0.42371080
#> 190   pop      SMR    groupXV  8815538      37           0.06235229
#> 191   pop      UPD    groupXV  8815538      37           0.72493893
#> 192   pop      ASP    groupXI  9019877      38           0.75592611
#> 193   pop      OPL    groupXI  9019877      38           0.26458852
#> 194   pop      PAL    groupXI  9019877      38           0.26996961
#> 195   pop      SMR    groupXI  9019877      38           0.64991795
#> 196   pop      UPD    groupXI  9019877      38           0.46059662
#> 197   pop      ASP   groupXXI  9067876      39           0.34278171
#> 198   pop      OPL   groupXXI  9067876      39           0.79214739
#> 199   pop      PAL   groupXXI  9067876      39           0.38691632
#> 200   pop      SMR   groupXXI  9067876      39           0.59969049
#> 201   pop      UPD   groupXXI  9067876      39           0.77815969
#> 202   pop      ASP  groupXVII  9108184      40           0.14009906
#> 203   pop      CLF  groupXVII  9108184      40           0.36131043
#> 204   pop      OPL  groupXVII  9108184      40           0.44393875
#> 205   pop      PAL  groupXVII  9108184      40           0.71439304
#> 206   pop      SMR  groupXVII  9108184      40           0.63584819
#> 207   pop      UPD  groupXVII  9108184      40           0.35609646
#> 208   pop      ASP    groupIX  9284427      41           0.76623033
#> 209   pop      OPL    groupIX  9284427      41           0.10998447
#> 210   pop      SMR    groupIX  9284427      41           0.88409996
#> 211   pop      UPD    groupIX  9284427      41           0.30062299
#> 212   pop      ASP     groupV  9767769      42           0.02845974
#> 213   pop      CLF     groupV  9767769      42           0.76508719
#> 214   pop      OPL     groupV  9767769      42           0.97980768
#> 215   pop      PAL     groupV  9767769      42           0.84893600
#> 216   pop      SMR     groupV  9767769      42           0.05146414
#> 217   pop      UPD     groupV  9767769      42           0.68486037
#> 218   pop      ASP   groupXII  9881312      43           0.23725583
#> 219   pop      CLF   groupXII  9881312      43           0.17090352
#> 220   pop      OPL   groupXII  9881312      43           0.04874057
#> 221   pop      SMR   groupXII  9881312      43           0.76858895
#> 222   pop      UPD   groupXII  9881312      43           0.01279059
#> 223   pop      ASP    groupXX 10296516      44           0.50233495
#> 224   pop      OPL    groupXX 10296516      44           0.08820496
#> 225   pop      PAL    groupXX 10296516      44           0.51860502
#> 226   pop      SMR    groupXX 10296516      44           0.43953590
#> 227   pop      UPD    groupXX 10296516      44           0.73531669
#> 228   pop      UPD     groupI 10593829      45           0.21484695
#> 229   pop      ASP   groupXIX 10765968      46           0.70291756
#> 230   pop      CLF   groupXIX 10765968      46           0.36131043
#> 231   pop      OPL   groupXIX 10765968      46           0.19708221
#> 232   pop      SMR   groupXIX 10765968      46           0.83615187
#> 233   pop      UPD   groupXIX 10765968      46           0.84760393
#> 234   pop      ASP    groupII 10829169      47           0.30305233
#> 235   pop      OPL    groupII 10829169      47           0.30871849
#> 236   pop      PAL    groupII 10829169      47           0.60477285
#> 237   pop      SMR    groupII 10829169      47           0.88824913
#> 238   pop      UPD    groupII 10829169      47           1.00000000
#> 239   pop      ASP    groupXX 10913115      48           0.23129260
#> 240   pop      CLF    groupXX 10913115      48           1.00000000
#> 241   pop      OPL    groupXX 10913115      48           0.73235632
#> 242   pop      PAL    groupXX 10913115      48           0.84170356
#> 243   pop      SMR    groupXX 10913115      48           0.15149399
#> 244   pop      UPD    groupXX 10913115      48           0.71009681
#> 245   pop      CLF groupXVIII 10974594      49           0.12133525
#> 246   pop      OPL groupXVIII 10974594      49           0.36205019
#> 247   pop      UPD groupXVIII 10974594      49           0.71538812
#> 248   pop      ASP  groupXVII 11256661      50           0.25979646
#> 249   pop      CLF  groupXVII 11256661      50           0.02534732
#> 250   pop      OPL  groupXVII 11256661      50           0.85923355
#> 251   pop      PAL  groupXVII 11256661      50           0.24156659
#> 252   pop      SMR  groupXVII 11256661      50           0.15608367
#> 253   pop      UPD  groupXVII 11256661      50           0.64880708
#> 254   pop      ASP  groupVIII 11438870      51           0.34278171
#> 255   pop      OPL  groupVIII 11438870      51           1.00000000
#> 256   pop      PAL  groupVIII 11438870      51           0.55368023
#> 257   pop      SMR  groupVIII 11438870      51           0.89940035
#> 258   pop      UPD  groupVIII 11438870      51           0.37515543
#> 259   pop      ASP   groupVII 11477223      52           0.12471041
#> 260   pop      CLF   groupVII 11477223      52           0.32718688
#> 261   pop      OPL   groupVII 11477223      52           0.51269076
#> 262   pop      PAL   groupVII 11477223      52           0.09245960
#> 263   pop      SMR   groupVII 11477223      52           0.54695590
#> 264   pop      UPD   groupVII 11477223      52           0.32112107
#> 265   pop      ASP    groupII 11505084      53           0.17308015
#> 266   pop      CLF    groupII 11505084      53           0.50498508
#> 267   pop      OPL    groupII 11505084      53           0.74257454
#> 268   pop      SMR    groupII 11505084      53           0.53817923
#> 269   pop      UPD    groupII 11505084      53           1.00000000
#> 270   pop      CLF    groupIX 11904651      54           0.62420611
#> 271   pop      OPL    groupIX 11904651      54           0.34509047
#> 272   pop      SMR    groupIX 11904651      54           0.32760715
#> 273   pop      ASP  groupVIII 12108945      55           0.14412703
#> 274   pop      CLF  groupVIII 12108945      55           0.27332168
#> 275   pop      OPL  groupVIII 12108945      55           0.30688510
#> 276   pop      PAL  groupVIII 12108945      55           0.58135198
#> 277   pop      SMR  groupVIII 12108945      55           0.69599510
#> 278   pop      UPD  groupVIII 12108945      55           0.28155517
#> 279   pop      ASP   groupXIX 12485182      56           0.57075039
#> 280   pop      CLF   groupXIX 12485182      56           0.62420611
#> 281   pop      OPL   groupXIX 12485182      56           0.85217893
#> 282   pop      PAL   groupXIX 12485182      56           0.19666685
#> 283   pop      SMR   groupXIX 12485182      56           0.36146637
#> 284   pop      UPD   groupXIX 12485182      56           0.20489389
#> 285   pop      ASP    groupXI 12833733      57           0.28504941
#> 286   pop      OPL    groupXI 12833733      57           0.36868827
#> 287   pop      SMR    groupXI 12833733      57           0.35738572
#> 288   pop      ASP groupXVIII 12937786      58           0.15682830
#> 289   pop      OPL groupXVIII 12937786      58           0.87017980
#> 290   pop      PAL groupXVIII 12937786      58           0.95563269
#> 291   pop      SMR groupXVIII 12937786      58           0.41011699
#> 292   pop      UPD groupXVIII 12937786      58           0.86199249
#> 293   pop      ASP   groupIII 13177027      59           0.28083102
#> 294   pop      CLF   groupIII 13177027      59           0.12133525
#> 295   pop      OPL   groupIII 13177027      59           0.35813199
#> 296   pop      PAL   groupIII 13177027      59           0.24156659
#> 297   pop      SMR   groupIII 13177027      59           0.73151622
#> 298   pop      UPD   groupIII 13177027      59           0.34499857
#> 299   pop      ASP    groupVI 13335577      60           0.08766880
#> 300   pop      CLF    groupVI 13335577      60           0.27332168
#> 301   pop      OPL    groupVI 13335577      60           0.20697738
#> 302   pop      SMR    groupVI 13335577      60           0.01263700
#> 303   pop      UPD    groupVI 13335577      60           0.15750435
#> 304   pop      CLF     groupX 13718292      61           0.38647623
#> 305   pop      OPL     groupX 13718292      61           0.38846226
#> 306   pop      SMR     groupX 13718292      61           0.13603713
#> 307   pop      UPD     groupX 13718292      61           0.19043026
#> 308   pop      ASP    groupVI 13895716      62           0.19670560
#> 309   pop      OPL    groupVI 13895716      62           0.27942577
#> 310   pop      SMR    groupVI 13895716      62           0.08851919
#> 311   pop      UPD    groupVI 13895716      62           0.70687919
#> 312   pop      ASP   groupXIV 14009990      63           0.16754628
#> 313   pop      OPL   groupXIV 14009990      63           0.59588309
#> 314   pop      PAL   groupXIV 14009990      63           0.19670560
#> 315   pop      SMR   groupXIV 14009990      63           0.75257628
#> 316   pop      UPD   groupXIV 14009990      63           0.88412209
#> 317   pop      ASP   groupXII 14158350      64           0.59483463
#> 318   pop      CLF   groupXII 14158350      64           0.62420611
#> 319   pop      OPL   groupXII 14158350      64           0.10350435
#> 320   pop      PAL   groupXII 14158350      64           0.11899072
#> 321   pop      SMR   groupXII 14158350      64           0.77397612
#> 322   pop      UPD   groupXII 14158350      64           0.58135198
#> 323   pop      ASP    groupIX 14164215      65           0.47726687
#> 324   pop      CLF    groupIX 14164215      65           0.79048223
#> 325   pop      OPL    groupIX 14164215      65           0.17062194
#> 326   pop      PAL    groupIX 14164215      65           0.56570671
#> 327   pop      SMR    groupIX 14164215      65           0.46959643
#> 328   pop      UPD    groupIX 14164215      65           0.50262824
#> 329   pop      ASP    groupVI 14236496      66           1.00000000
#> 330   pop      CLF    groupVI 14236496      66           0.79048223
#> 331   pop      OPL    groupVI 14236496      66           0.80864089
#> 332   pop      PAL    groupVI 14236496      66           0.10589897
#> 333   pop      SMR    groupVI 14236496      66           0.93757994
#> 334   pop      UPD    groupVI 14236496      66           0.86420308
#> 335   pop      PAL     groupX 14273796      67           0.30631186
#> 336   pop      UPD     groupX 14273796      67           0.18798576
#> 337   pop      ASP  groupVIII 14445389      68           0.60701115
#> 338   pop      CLF  groupVIII 14445389      68           0.27332168
#> 339   pop      OPL  groupVIII 14445389      68           0.45864925
#> 340   pop      PAL  groupVIII 14445389      68           0.83960282
#> 341   pop      SMR  groupVIII 14445389      68           0.25217626
#> 342   pop      UPD  groupVIII 14445389      68           0.21042229
#> 343   pop      CLF    groupIX 14868533      69           0.02534732
#> 344   pop      OPL    groupIX 14868533      69           0.64390193
#> 345   pop      SMR    groupIX 14868533      69           0.74832825
#> 346   pop      ASP    groupIX 15190693      70           0.09763454
#> 347   pop      CLF    groupIX 15190693      70           0.76508719
#> 348   pop      OPL    groupIX 15190693      70           0.52826533
#> 349   pop      SMR    groupIX 15190693      70           0.39180485
#> 350   pop      UPD    groupIX 15190693      70           0.19821062
#> 351   pop      ASP   groupXII 15306362      71           0.42787981
#> 352   pop      OPL   groupXII 15306362      71           0.25456479
#> 353   pop      PAL   groupXII 15306362      71           0.78972582
#> 354   pop      SMR   groupXII 15306362      71           0.81302534
#> 355   pop      UPD   groupXII 15306362      71           0.21484695
#> 356   pop      ASP    groupXX 15460430      72           0.92959545
#> 357   pop      CLF    groupXX 15460430      72           0.18403863
#> 358   pop      OPL    groupXX 15460430      72           0.58168850
#> 359   pop      PAL    groupXX 15460430      72           0.64183818
#> 360   pop      SMR    groupXX 15460430      72           0.63903992
#> 361   pop      UPD    groupXX 15460430      72           0.13480237
#> 362   pop      ASP   groupVII 15493025      73           0.59816145
#> 363   pop      CLF   groupVII 15493025      73           0.27332168
#> 364   pop      OPL   groupVII 15493025      73           0.40845145
#> 365   pop      PAL   groupVII 15493025      73           0.16754628
#> 366   pop      SMR   groupVII 15493025      73           0.37485714
#> 367   pop      UPD   groupVII 15493025      73           0.65635646
#> 368   pop      ASP    groupXV 15547073      74           0.88030337
#> 369   pop      OPL    groupXV 15547073      74           0.21009563
#> 370   pop      SMR    groupXV 15547073      74           0.45660451
#> 371   pop      UPD    groupXV 15547073      74           0.83509889
#> 372   pop      ASP  groupXIII 15697285      75           0.27102765
#> 373   pop      CLF  groupXIII 15697285      75           1.00000000
#> 374   pop      OPL  groupXIII 15697285      75           0.69975960
#> 375   pop      PAL  groupXIII 15697285      75           0.51530470
#> 376   pop      SMR  groupXIII 15697285      75           0.73818101
#> 377   pop      UPD  groupXIII 15697285      75           0.66695453
#> 378   pop      ASP   groupIII 15983994      76           0.40081417
#> 379   pop      CLF   groupIII 15983994      76           0.38647623
#> 380   pop      OPL   groupIII 15983994      76           0.19228004
#> 381   pop      SMR   groupIII 15983994      76           1.00000000
#> 382   pop      UPD   groupIII 15983994      76           0.83001131
#> 383   pop      ASP   groupXIX 16233699      77           0.13361440
#> 384   pop      CLF   groupXIX 16233699      77           0.62420611
#> 385   pop      OPL   groupXIX 16233699      77           0.87885448
#> 386   pop      PAL   groupXIX 16233699      77           0.05653028
#> 387   pop      SMR   groupXIX 16233699      77           0.21009563
#> 388   pop      UPD   groupXIX 16233699      77           0.11899072
#> 389   pop      CLF    groupIX 16306966      78           0.47328947
#> 390   pop      OPL    groupIX 16306966      78           0.32571645
#> 391   pop      SMR    groupIX 16306966      78           0.45743774
#> 392   pop      ASP     groupI 16498302      79           0.47726687
#> 393   pop      OPL     groupI 16498302      79           0.09854293
#> 394   pop      PAL     groupI 16498302      79           0.77151915
#> 395   pop      SMR     groupI 16498302      79           0.47424459
#> 396   pop      UPD     groupI 16498302      79           0.09096895
#> 397   pop      ASP  groupXIII 16520929      80           0.43405434
#> 398   pop      CLF  groupXIII 16520929      80           0.43857803
#> 399   pop      OPL  groupXIII 16520929      80           0.78389165
#> 400   pop      PAL  groupXIII 16520929      80           0.31623269
#> 401   pop      SMR  groupXIII 16520929      80           0.82743508
#> 402   pop      UPD  groupXIII 16520929      80           0.16426437
#> 403   pop      ASP   groupXVI 16561891      81           1.00000000
#> 404   pop      CLF   groupXVI 16561891      81           0.57615012
#> 405   pop      OPL   groupXVI 16561891      81           0.83008015
#> 406   pop      PAL   groupXVI 16561891      81           0.05260189
#> 407   pop      SMR   groupXVI 16561891      81           0.57855968
#> 408   pop      UPD   groupXVI 16561891      81           0.94306703
#> 409   pop      ASP    groupII 16598364      82           0.74688563
#> 410   pop      OPL    groupII 16598364      82           0.90748021
#> 411   pop      SMR    groupII 16598364      82           0.81105464
#> 412   pop      UPD    groupII 16598364      82           0.01494343
#> 413   pop      ASP    groupIV 17248473      83           0.04404471
#> 414   pop      CLF    groupIV 17248473      83           0.62420611
#> 415   pop      SMR    groupIV 17248473      83           0.26678531
#> 416   pop      UPD    groupIV 17248473      83           0.16383568
#> 417   pop      ASP    groupXX 17748944      84           0.21189355
#> 418   pop      CLF    groupXX 17748944      84           0.62420611
#> 419   pop      OPL    groupXX 17748944      84           1.00000000
#> 420   pop      SMR    groupXX 17748944      84           0.08075872
#> 421   pop      UPD    groupXX 17748944      84           0.49901470
#> 422   pop      ASP  groupXIII 17842228      85           0.33208078
#> 423   pop      CLF  groupXIII 17842228      85           0.57615012
#> 424   pop      OPL  groupXIII 17842228      85           0.86392607
#> 425   pop      PAL  groupXIII 17842228      85           0.83681320
#> 426   pop      SMR  groupXIII 17842228      85           0.19788781
#> 427   pop      UPD  groupXIII 17842228      85           0.86461805
#> 428   pop      ASP   groupXIX 17858769      86           0.34577859
#> 429   pop      CLF   groupXIX 17858769      86           0.62420611
#> 430   pop      OPL   groupXIX 17858769      86           0.39073306
#> 431   pop      PAL   groupXIX 17858769      86           0.23725583
#> 432   pop      SMR   groupXIX 17858769      86           0.21009563
#> 433   pop      UPD   groupXIX 17858769      86           0.13480237
#> 434   pop      ASP    groupIX 17929225      87           1.00000000
#> 435   pop      CLF    groupIX 17929225      87           0.62420611
#> 436   pop      OPL    groupIX 17929225      87           0.27888628
#> 437   pop      SMR    groupIX 17929225      87           0.42853867
#> 438   pop      UPD    groupIX 17929225      87           0.13448098
#> 439   pop      ASP    groupIX 18339767      88           1.00000000
#> 440   pop      CLF    groupIX 18339767      88           0.27332168
#> 441   pop      OPL    groupIX 18339767      88           0.88896184
#> 442   pop      SMR    groupIX 18339767      88           0.39089689
#> 443   pop      UPD    groupIX 18339767      88           0.40762595
#> 444   pop      ASP  groupXIII 18836558      89           0.63525630
#> 445   pop      OPL  groupXIII 18836558      89           0.57046936
#> 446   pop      PAL  groupXIII 18836558      89           0.17971249
#> 447   pop      SMR  groupXIII 18836558      89           0.10501127
#> 448   pop      UPD  groupXIII 18836558      89           0.60952338
#> 449   pop      ASP    groupII 19244005      90           0.89819521
#> 450   pop      SMR    groupII 19244005      90           0.24369510
#> 451   pop      UPD    groupII 19244005      90           0.05181233
#> 452   pop      ASP    groupXX 19317273      91           0.57075039
#> 453   pop      OPL    groupXX 19317273      91           0.45325470
#> 454   pop      PAL    groupXX 19317273      91           0.86864517
#> 455   pop      SMR    groupXX 19317273      91           0.94880340
#> 456   pop      UPD    groupXX 19317273      91           0.82549791
#> 457   pop      ASP     groupI 19759331      92           1.00000000
#> 458   pop      CLF     groupI 19759331      92           0.36131043
#> 459   pop      OPL     groupI 19759331      92           1.00000000
#> 460   pop      PAL     groupI 19759331      92           0.10804347
#> 461   pop      SMR     groupI 19759331      92           0.81897085
#> 462   pop      UPD     groupI 19759331      92           1.00000000
#> 463   pop      ASP   groupVII 20165240      93           0.20075068
#> 464   pop      CLF   groupVII 20165240      93           0.27332168
#> 465   pop      OPL   groupVII 20165240      93           0.11479352
#> 466   pop      PAL   groupVII 20165240      93           0.59670122
#> 467   pop      SMR   groupVII 20165240      93           0.51458355
#> 468   pop      UPD   groupVII 20165240      93           0.28008721
#> 469   pop      OPL    groupII 20496069      94           0.65004138
#> 470   pop      PAL    groupII 20496069      94           0.43702909
#> 471   pop      SMR    groupII 20496069      94           0.30605608
#> 472   pop      UPD    groupII 20496069      94           0.92265025
#> 473   pop      ASP     groupI 21284266      95           0.42919530
#> 474   pop      CLF     groupI 21284266      95           0.62420611
#> 475   pop      OPL     groupI 21284266      95           0.85403445
#> 476   pop      SMR     groupI 21284266      95           0.40518736
#> 477   pop      UPD     groupI 21284266      95           0.66500554
#> 478   pop      ASP    groupIV 22227311      96           0.40762595
#> 479   pop      CLF    groupIV 22227311      96           0.62420611
#> 480   pop      OPL    groupIV 22227311      96           0.64496464
#> 481   pop      SMR    groupIV 22227311      96           0.23128935
#> 482   pop      UPD    groupIV 22227311      96           0.05181233
#> 483   pop      CLF     groupI 22746027      97           0.01430588
#> 484   pop      OPL     groupI 22746027      97           0.94349355
#> 485   pop      SMR     groupI 22746027      97           0.09934362
#> 486   pop      UPD     groupI 22746027      97           0.20489389
#> 487   pop      ASP   groupVII 25239495      98           0.23129260
#> 488   pop      CLF   groupVII 25239495      98           0.57615012
#> 489   pop      OPL   groupVII 25239495      98           0.90748021
#> 490   pop      PAL   groupVII 25239495      98           0.80391805
#> 491   pop      SMR   groupVII 25239495      98           0.24499927
#> 492   pop      UPD   groupVII 25239495      98           0.66399272
#> 493   pop      ASP   groupVII 27599144      99           0.71223539
#> 494   pop      CLF   groupVII 27599144      99           0.30910812
#> 495   pop      OPL   groupVII 27599144      99           0.65535968
#> 496   pop      SMR   groupVII 27599144      99           0.76250041
#> 497   pop      UPD   groupVII 27599144      99           0.35931734
#> 498   pop      ASP    groupIV 32321895     100           0.01572530
#> 499   pop      OPL    groupIV 32321895     100           0.27829645
#> 500   pop      PAL    groupIV 32321895     100           0.44087650
#> 501   pop      SMR    groupIV 32321895     100           0.03526631
#> 502   pop      UPD    groupIV 32321895     100           0.95285346
#> 
```
