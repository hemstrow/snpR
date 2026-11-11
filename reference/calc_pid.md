# Probability of Identity

Calculates the probability of sampling two individuals from the same
population with identical genotypes.

## Usage

``` r
calc_pid(x, facets = NULL)
```

## Arguments

- x:

  snpRdata object

- facets:

  facets over which to calculate PID. See
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md)
  for details. SNP facets are ignored.

## Value

A snpRdata object with PID info merged in, fetchable with
[`get.snpR.stats`](https://hemstrow.github.io/snpR/reference/get.snpR.stats.md)

## Details

Both per-locus and across all locus probabilities are calculated. Two
statistics are calculated:

1.  PID: unbiased probability of randomly sampling two identical
    individuals from a population, corrected for sample size.

2.  PIDsib: probability of randomly sampling two genetically identical
    siblings from a population.

## References

Paetkau, D., Waits, L. P., Clarkson, P. L., Craighead, L., Vyse, E.,
Ward, R., & Strobeck, C. (1998). Variation in Genetic Diversity across
the Range of North American Brown Bears. Conservation Biology, 12(2),
418–429. https://doi.org/10.1111/j.1523-1739.1998.96457.x

## Author

William Hemstrom

## Examples

``` r
x <- calc_pid(stickSNPs, "pop")
get.snpR.stats(x, "pop", "pid")
#> $single
#>     facet subfacet        chr position .snp.id        PID    PIDsib
#> 1     pop      ASP    groupVI   212436       1  0.3671875 0.5937500
#> 2     pop      CLF    groupVI   212436       1  0.2781250 0.6461500
#> 3     pop      OPL    groupVI   212436       1  0.3682566 0.6008287
#> 4     pop      PAL    groupVI   212436       1  0.4602482 0.7131062
#> 5     pop      SMR    groupVI   212436       1  0.3708817 0.6041830
#> 6     pop      UPD    groupVI   212436       1  0.3635417 0.5937500
#> 7     pop      ASP groupXVIII   350070       2  0.7648562 0.9083844
#> 8     pop      CLF groupXVIII   350070       2  1.0000000 1.0000000
#> 9     pop      OPL groupXVIII   350070       2  0.7543919 0.8836468
#> 10    pop      PAL groupXVIII   350070       2  1.0000000 1.0000000
#> 11    pop      SMR groupXVIII   350070       2  0.8181818 0.9160470
#> 12    pop      UPD groupXVIII   350070       2  1.0000000 1.0000000
#> 13    pop      ASP    groupXV   854285       3  0.6363636 0.8449559
#> 14    pop      CLF    groupXV   854285       3  0.5302083 0.8559751
#> 15    pop      OPL    groupXV   854285       3  0.7133333 0.8609254
#> 16    pop      PAL    groupXV   854285       3  0.8631153 0.9408278
#> 17    pop      SMR    groupXV   854285       3  0.7543919 0.8836468
#> 18    pop      UPD    groupXV   854285       3  0.7333333 0.8813630
#> 19    pop      ASP  groupXIII   907279       4  0.5276989 0.7852670
#> 20    pop      CLF  groupXIII   907279       4  0.3531250 0.8321500
#> 21    pop      OPL  groupXIII   907279       4  0.3682566 0.6008287
#> 22    pop      PAL  groupXIII   907279       4  1.0000000 1.0000000
#> 23    pop      SMR  groupXIII   907279       4  0.3903603 0.6329294
#> 24    pop      UPD  groupXIII   907279       4  0.3641827 0.6011990
#> 25    pop      ASP    groupVI  1068620       5  0.3504464 0.6198849
#> 26    pop      CLF    groupVI  1068620       5        NaN 1.0000000
#> 27    pop      OPL    groupVI  1068620       5  0.3790210 0.6420293
#> 28    pop      PAL    groupVI  1068620       5  0.4332386 0.7109420
#> 29    pop      SMR    groupVI  1068620       5  0.4268430 0.6902241
#> 30    pop      UPD    groupVI  1068620       5  0.3650794 0.5976223
#> 31    pop      ASP     groupV  1392753       6  0.3843006 0.6777344
#> 32    pop      CLF     groupV  1392753       6  0.2968750 0.6777344
#> 33    pop      OPL     groupV  1392753       6  0.7810671 0.8957053
#> 34    pop      PAL     groupV  1392753       6  0.3673735 0.6048522
#> 35    pop      SMR     groupV  1392753       6  0.7327303 0.8739458
#> 36    pop      UPD     groupV  1392753       6  0.3706294 0.6198849
#> 37    pop      ASP    groupIX  2075079       7  0.7901988 0.9160470
#> 38    pop      CLF    groupIX  2075079       7  0.3333333 0.7511574
#> 39    pop      OPL    groupIX  2075079       7  0.5624060 0.7775927
#> 40    pop      PAL    groupIX  2075079       7  1.0000000 1.0000000
#> 41    pop      SMR    groupIX  2075079       7  0.4887218 0.7272741
#> 42    pop      UPD    groupIX  2075079       7  1.0000000 1.0000000
#> 43    pop      ASP  groupXIII  2106000       8  0.4083807 0.6950238
#> 44    pop      CLF  groupXIII  2106000       8  0.3802083 0.6025029
#> 45    pop      OPL  groupXIII  2106000       8  0.5275739 0.7511574
#> 46    pop      PAL  groupXIII  2106000       8  1.0000000 1.0000000
#> 47    pop      SMR  groupXIII  2106000       8  0.3801858 0.6226344
#> 48    pop      UPD  groupXIII  2106000       8  0.3843006 0.6777344
#> 49    pop      ASP   groupXII  2136587       9  0.4083807 0.6950238
#> 50    pop      CLF   groupXII  2136587       9  0.3333333 0.7511574
#> 51    pop      OPL   groupXII  2136587       9  0.6608510 0.8321500
#> 52    pop      PAL   groupXII  2136587       9  0.8631153 0.9408278
#> 53    pop      SMR   groupXII  2136587       9  0.4272446 0.6777344
#> 54    pop      UPD   groupXII  2136587       9  0.3944857 0.6645721
#> 55    pop      ASP  groupXVII  2560839      10  0.6363636 0.8449559
#> 56    pop      CLF  groupXVII  2560839      10  0.5302083 0.8559751
#> 57    pop      OPL  groupXVII  2560839      10  0.9154887 0.9613762
#> 58    pop      PAL  groupXVII  2560839      10  1.0000000 1.0000000
#> 59    pop      SMR  groupXVII  2560839      10  0.8181818 0.9160470
#> 60    pop      UPD  groupXVII  2560839      10  0.5384615 0.7775927
#> 61    pop      ASP   groupIII  2750695      11  0.3640625 0.5968844
#> 62    pop      CLF   groupIII  2750695      11  0.4000000 0.5937500
#> 63    pop      OPL   groupIII  2750695      11  0.3708882 0.6052113
#> 64    pop      PAL   groupIII  2750695      11  1.0000000 1.0000000
#> 65    pop      SMR   groupIII  2750695      11  0.3650174 0.5937500
#> 66    pop      UPD   groupIII  2750695      11  0.3665865 0.6082972
#> 67    pop      ASP   groupXIX  3008524      12  0.3843006 0.6777344
#> 68    pop      CLF   groupXIX  3008524      12  0.5302083 0.8559751
#> 69    pop      OPL   groupXIX  3008524      12  0.4765370 0.7235306
#> 70    pop      PAL   groupXIX  3008524      12  0.3855311 0.6395493
#> 71    pop      SMR   groupXIX  3008524      12  0.4970588 0.7410232
#> 72    pop      UPD   groupXIX  3008524      12  0.3666667 0.6461500
#> 73    pop      ASP   groupXVI  3624417      13  0.5555556 0.8171011
#> 74    pop      CLF   groupXVI  3624417      13  0.3802083 0.6025029
#> 75    pop      OPL   groupXVI  3624417      13  0.4333710 0.6777344
#> 76    pop      PAL   groupXVI  3624417      13  0.8529419 0.9371130
#> 77    pop      SMR   groupXVI  3624417      13  0.4749226 0.7184000
#> 78    pop      UPD   groupXVI  3624417      13  0.5128205 0.7650643
#> 79    pop      ASP   groupVII  3667528      14  0.4222222 0.7184000
#> 80    pop      CLF   groupVII  3667528      14  0.3333333 0.7511574
#> 81    pop      OPL   groupVII  3667528      14  0.7133333 0.8609254
#> 82    pop      PAL   groupVII  3667528      14  1.0000000 1.0000000
#> 83    pop      SMR   groupVII  3667528      14  0.5021189 0.7356738
#> 84    pop      UPD   groupVII  3667528      14  0.3707990 0.6164463
#> 85    pop      ASP   groupVII  4189120      15  0.4083807 0.6950238
#> 86    pop      CLF   groupVII  4189120      15  1.0000000 1.0000000
#> 87    pop      OPL   groupVII  4189120      15  0.4749226 0.7184000
#> 88    pop      PAL   groupVII  4189120      15  0.5833333 0.7991943
#> 89    pop      SMR   groupVII  4189120      15  0.4611137 0.7057867
#> 90    pop      UPD   groupVII  4189120      15  0.3636364 0.5956024
#> 91    pop      ASP    groupXX  4538074      16  0.3843006 0.6777344
#> 92    pop      CLF    groupXX  4538074      16  0.2781250 0.6461500
#> 93    pop      OPL    groupXX  4538074      16  0.3876064 0.6334381
#> 94    pop      PAL    groupXX  4538074      16  0.3650298 0.5993352
#> 95    pop      SMR    groupXX  4538074      16  0.3876064 0.6334381
#> 96    pop      UPD    groupXX  4538074      16  0.4332386 0.7109420
#> 97    pop      ASP  groupXIII  4706089      17  0.4902530 0.7693844
#> 98    pop      CLF  groupXIII  4706089      17  1.0000000 1.0000000
#> 99    pop      OPL  groupXIII  4706089      17  0.5745500 0.7809442
#> 100   pop      PAL  groupXIII  4706089      17  0.4785714 0.7301254
#> 101   pop      SMR  groupXIII  4706089      17  0.5222039 0.7511574
#> 102   pop      UPD  groupXIII  4706089      17  0.5128205 0.7650643
#> 103   pop      ASP   groupIII  4746808      18  0.3666667 0.6461500
#> 104   pop      CLF   groupIII  4746808      18  0.4000000 0.6064000
#> 105   pop      OPL   groupIII  4746808      18  0.4816035 0.7150504
#> 106   pop      PAL   groupIII  4746808      18  0.3906094 0.6543107
#> 107   pop      SMR   groupIII  4746808      18  0.5495924 0.7669341
#> 108   pop      UPD   groupIII  4746808      18  0.3650298 0.5993352
#> 109   pop      ASP  groupXVII  5041995      19  1.0000000 1.0000000
#> 110   pop      CLF  groupXVII  5041995      19  0.2968750 0.6777344
#> 111   pop      OPL  groupXVII  5041995      19  0.7204238 0.8684693
#> 112   pop      PAL  groupXVII  5041995      19  1.0000000 1.0000000
#> 113   pop      SMR  groupXVII  5041995      19  0.7543919 0.8836468
#> 114   pop      UPD  groupXVII  5041995      19  0.7333333 0.8813630
#> 115   pop      ASP     groupX  5161637      20  0.7322668 0.8991912
#> 116   pop      CLF     groupX  5161637      20  0.5302083 0.8559751
#> 117   pop      OPL     groupX  5161637      20  0.9117451 0.9598256
#> 118   pop      PAL     groupX  5161637      20  1.0000000 1.0000000
#> 119   pop      SMR     groupX  5161637      20  0.7543919 0.8836468
#> 120   pop      UPD     groupX  5161637      20  0.7142857 0.8739458
#> 121   pop      ASP   groupXII  5178037      21  0.4475446 0.7511574
#> 122   pop      CLF   groupXII  5178037      21  0.4531250 0.6136475
#> 123   pop      OPL   groupXII  5178037      21  0.3650174 0.5937500
#> 124   pop      PAL   groupXII  5178037      21  1.0000000 1.0000000
#> 125   pop      SMR   groupXII  5178037      21  0.3876064 0.6334381
#> 126   pop      UPD   groupXII  5178037      21  0.3696970 0.6296296
#> 127   pop      ASP    groupIX  5287711      22  0.3650794 0.5976223
#> 128   pop      CLF    groupIX  5287711      22  0.4000000 0.6064000
#> 129   pop      OPL    groupIX  5287711      22  0.3946977 0.6410760
#> 130   pop      PAL    groupIX  5287711      22  1.0000000 1.0000000
#> 131   pop      SMR    groupIX  5287711      22  0.5473684 0.7693844
#> 132   pop      UPD    groupIX  5287711      22  0.3706294 0.6198849
#> 133   pop      ASP  groupVIII  5462135      23  0.7648562 0.9083844
#> 134   pop      CLF  groupVIII  5462135      23  1.0000000 1.0000000
#> 135   pop      OPL  groupVIII  5462135      23  0.8400000 0.9254118
#> 136   pop      PAL  groupVIII  5462135      23  1.0000000 1.0000000
#> 137   pop      SMR  groupVIII  5462135      23  0.6526316 0.8321500
#> 138   pop      UPD  groupVIII  5462135      23  0.8529419 0.9371130
#> 139   pop      ASP    groupIX  5960933      24  0.7648562 0.9083844
#> 140   pop      CLF    groupIX  5960933      24  1.0000000 1.0000000
#> 141   pop      OPL    groupIX  5960933      24  0.6716189 0.8375130
#> 142   pop      PAL    groupIX  5960933      24  1.0000000 1.0000000
#> 143   pop      SMR    groupIX  5960933      24  0.7894737 0.9040063
#> 144   pop      UPD    groupIX  5960933      24  0.6148383 0.8224007
#> 145   pop      ASP  groupXIII  6076868      25  0.3666667 0.6461500
#> 146   pop      CLF  groupXIII  6076868      25  0.3531250 0.8321500
#> 147   pop      OPL  groupXIII  6076868      25  0.4240542 0.6701341
#> 148   pop      PAL  groupXIII  6076868      25  1.0000000 1.0000000
#> 149   pop      SMR  groupXIII  6076868      25  0.3650174 0.5937500
#> 150   pop      UPD  groupXIII  6076868      25  0.3650298 0.5993352
#> 151   pop      ASP    groupXI  6177332      26  0.3571429 0.6296296
#> 152   pop      CLF    groupXI  6177332      26  1.0000000 1.0000000
#> 153   pop      OPL    groupXI  6177332      26  0.5624060 0.7775927
#> 154   pop      PAL    groupXI  6177332      26  0.3641026 0.5951407
#> 155   pop      SMR    groupXI  6177332      26  0.4023062 0.6486743
#> 156   pop      UPD    groupXI  6177332      26  0.3802083 0.6508292
#> 157   pop      ASP    groupIV  6291287      27  0.7648562 0.9083844
#> 158   pop      CLF    groupIV  6291287      27  0.5302083 0.8559751
#> 159   pop      OPL    groupIV  6291287      27  0.4725202 0.7109420
#> 160   pop      PAL    groupIV  6291287      27  1.0000000 1.0000000
#> 161   pop      SMR    groupIV  6291287      27  0.5473684 0.7693844
#> 162   pop      UPD    groupIV  6291287      27  0.7142857 0.8739458
#> 163   pop      ASP    groupIX  6515796      28  0.5000000 0.7991943
#> 164   pop      CLF    groupIX  6515796      28  0.3333333 0.7511574
#> 165   pop      OPL    groupIX  6515796      28  0.3811745 0.6258457
#> 166   pop      PAL    groupIX  6515796      28  0.6374657 0.8321500
#> 167   pop      SMR    groupIX  6515796      28  0.3656250 0.5945318
#> 168   pop      UPD    groupIX  6515796      28  0.3906094 0.6543107
#> 169   pop      ASP   groupXII  6718079      29  0.7322668 0.8991912
#> 170   pop      CLF   groupXII  6718079      29 -0.1302083 0.7991943
#> 171   pop      OPL   groupXII  6718079      29  0.4204496 0.6694322
#> 172   pop      PAL   groupXII  6718079      29  0.5034913 0.7511574
#> 173   pop      SMR   groupXII  6718079      29  0.4268430 0.6902241
#> 174   pop      UPD   groupXII  6718079      29  0.5893575 0.8114857
#> 175   pop      ASP   groupVII  6967477      30  0.6363636 0.8449559
#> 176   pop      CLF   groupVII  6967477      30  0.2000000 0.7184000
#> 177   pop      OPL   groupVII  6967477      30  0.7327303 0.8739458
#> 178   pop      PAL   groupVII  6967477      30  1.0000000 1.0000000
#> 179   pop      SMR   groupVII  6967477      30  0.6241563 0.8137715
#> 180   pop      UPD   groupVII  6967477      30  0.4848485 0.7511574
#> 181   pop      ASP   groupIII  7194069      31  0.3818182 0.6623181
#> 182   pop      CLF   groupIII  7194069      31  0.5302083 0.8559751
#> 183   pop      OPL   groupIII  7194069      31  0.3943610 0.6358246
#> 184   pop      PAL   groupIII  7194069      31  1.0000000 1.0000000
#> 185   pop      SMR   groupIII  7194069      31  0.3751290 0.6158409
#> 186   pop      UPD   groupIII  7194069      31  0.3650298 0.5993352
#> 187   pop      ASP   groupXVI  7199722      32  0.4902530 0.7693844
#> 188   pop      CLF   groupXVI  7199722      32  0.2000000 0.7184000
#> 189   pop      OPL   groupXVI  7199722      32  0.5900621 0.7924536
#> 190   pop      PAL   groupXVI  7199722      32  0.3670330 0.6064000
#> 191   pop      SMR   groupXVI  7199722      32  0.5075755 0.7425162
#> 192   pop      UPD   groupXVI  7199722      32  0.3944857 0.6645721
#> 193   pop      ASP  groupVIII  7224065      33  0.3818182 0.6623181
#> 194   pop      CLF  groupVIII  7224065      33  0.4531250 0.6136475
#> 195   pop      OPL  groupVIII  7224065      33  0.3751290 0.6158409
#> 196   pop      PAL  groupVIII  7224065      33  0.4268430 0.6902241
#> 197   pop      SMR  groupVIII  7224065      33  0.3845801 0.6263483
#> 198   pop      UPD  groupVIII  7224065      33  0.3703016 0.6241761
#> 199   pop      ASP groupXVIII  7229426      34  0.3602431 0.6094107
#> 200   pop      CLF groupXVIII  7229426      34  0.4000000 0.5937500
#> 201   pop      OPL groupXVIII  7229426      34  0.3709150 0.6094107
#> 202   pop      PAL groupXVIII  7229426      34  1.0000000 1.0000000
#> 203   pop      SMR groupXVIII  7229426      34  0.4023062 0.6486743
#> 204   pop      UPD groupXVIII  7229426      34  0.5276989 0.7852670
#> 205   pop      ASP groupXVIII  8087022      35  0.6000000 0.8321500
#> 206   pop      CLF groupXVIII  8087022      35  1.0000000 1.0000000
#> 207   pop      OPL groupXVIII  8087022      35  0.7440127 0.8789886
#> 208   pop      PAL groupXVIII  8087022      35  1.0000000 1.0000000
#> 209   pop      SMR groupXVIII  8087022      35  0.7204238 0.8684693
#> 210   pop      UPD groupXVIII  8087022      35  0.7142857 0.8739458
#> 211   pop      ASP    groupIX  8373442      36  0.4902530 0.7693844
#> 212   pop      CLF    groupIX  8373442      36  0.3333333 0.7511574
#> 213   pop      OPL    groupIX  8373442      36  0.4978076 0.7301254
#> 214   pop      PAL    groupIX  8373442      36  1.0000000 1.0000000
#> 215   pop      SMR    groupIX  8373442      36  0.5021189 0.7356738
#> 216   pop      UPD    groupIX  8373442      36  0.6923077 0.8655509
#> 217   pop      ASP    groupXV  8815538      37  0.4222222 0.7184000
#> 218   pop      CLF    groupXV  8815538      37  0.2000000 0.7184000
#> 219   pop      OPL    groupXV  8815538      37  0.4185915 0.6635338
#> 220   pop      PAL    groupXV  8815538      37  0.8631153 0.9408278
#> 221   pop      SMR    groupXV  8815538      37  0.3786284 0.6175125
#> 222   pop      UPD    groupXV  8815538      37  0.3944857 0.6645721
#> 223   pop      ASP    groupXI  9019877      38  0.3640625 0.5968844
#> 224   pop      CLF    groupXI  9019877      38  1.0000000 1.0000000
#> 225   pop      OPL    groupXI  9019877      38  0.4858572 0.7224043
#> 226   pop      PAL    groupXI  9019877      38  0.3673735 0.6048522
#> 227   pop      SMR    groupXI  9019877      38  0.3680771 0.6015991
#> 228   pop      UPD    groupXI  9019877      38  0.4102928 0.6777344
#> 229   pop      ASP   groupXXI  9067876      39  0.4475446 0.7511574
#> 230   pop      CLF   groupXXI  9067876      39       -Inf 0.6296296
#> 231   pop      OPL   groupXXI  9067876      39  0.6374657 0.8321500
#> 232   pop      PAL   groupXXI  9067876      39  0.5276989 0.7852670
#> 233   pop      SMR   groupXXI  9067876      39  0.6757944 0.8488086
#> 234   pop      UPD   groupXXI  9067876      39  0.4902530 0.7693844
#> 235   pop      ASP  groupXVII  9108184      40  0.3640625 0.5968844
#> 236   pop      CLF  groupXVII  9108184      40  0.2000000 0.7184000
#> 237   pop      OPL  groupXVII  9108184      40  0.5075755 0.7425162
#> 238   pop      PAL  groupXVII  9108184      40  0.5034913 0.7511574
#> 239   pop      SMR  groupXVII  9108184      40  0.5314757 0.7605912
#> 240   pop      UPD  groupXVII  9108184      40  0.3629261 0.6041830
#> 241   pop      ASP    groupIX  9284427      41  0.3629261 0.6041830
#> 242   pop      CLF    groupIX  9284427      41  0.2781250 0.6461500
#> 243   pop      OPL    groupIX  9284427      41  0.4160647 0.6685876
#> 244   pop      PAL    groupIX  9284427      41  1.0000000 1.0000000
#> 245   pop      SMR    groupIX  9284427      41  0.3837776 0.6341321
#> 246   pop      UPD    groupIX  9284427      41  0.3619048 0.6064000
#> 247   pop      ASP     groupV  9767769      42  0.3571429 0.6777344
#> 248   pop      CLF     groupV  9767769      42  0.4000000 0.6064000
#> 249   pop      OPL     groupV  9767769      42  0.3672801 0.5961183
#> 250   pop      PAL     groupV  9767769      42  0.7500000 0.8879623
#> 251   pop      SMR     groupV  9767769      42  0.3685741 0.5995934
#> 252   pop      UPD     groupV  9767769      42  0.3706294 0.6198849
#> 253   pop      ASP   groupXII  9881312      43  0.6363636 0.8449559
#> 254   pop      CLF   groupXII  9881312      43  0.2000000 0.7184000
#> 255   pop      OPL   groupXII  9881312      43  0.4887218 0.7272741
#> 256   pop      PAL   groupXII  9881312      43  1.0000000 1.0000000
#> 257   pop      SMR   groupXII  9881312      43  0.6493933 0.8264299
#> 258   pop      UPD   groupXII  9881312      43  0.4195804 0.6922377
#> 259   pop      ASP    groupXX 10296516      44  0.4475446 0.7511574
#> 260   pop      CLF    groupXX 10296516      44  1.0000000 1.0000000
#> 261   pop      OPL    groupXX 10296516      44  0.7810671 0.8957053
#> 262   pop      PAL    groupXX 10296516      44  0.3670330 0.6064000
#> 263   pop      SMR    groupXX 10296516      44  0.6102365 0.8067477
#> 264   pop      UPD    groupXX 10296516      44  0.4332386 0.7109420
#> 265   pop      ASP     groupI 10593829      45  1.0000000 1.0000000
#> 266   pop      CLF     groupI 10593829      45  1.0000000 1.0000000
#> 267   pop      OPL     groupI 10593829      45  1.0000000 1.0000000
#> 268   pop      PAL     groupI 10593829      45  1.0000000 1.0000000
#> 269   pop      SMR     groupI 10593829      45  1.0000000 1.0000000
#> 270   pop      UPD     groupI 10593829      45  0.7333333 0.8813630
#> 271   pop      ASP   groupXIX 10765968      46  0.3671875 0.5986557
#> 272   pop      CLF   groupXIX 10765968      46  0.2000000 0.7184000
#> 273   pop      OPL   groupXIX 10765968      46  0.3678768 0.6025029
#> 274   pop      PAL   groupXIX 10765968      46  1.0000000 1.0000000
#> 275   pop      SMR   groupXIX 10765968      46  0.4602482 0.7131062
#> 276   pop      UPD   groupXIX 10765968      46  0.3703016 0.6241761
#> 277   pop      ASP    groupII 10829169      47  0.4083807 0.6950238
#> 278   pop      CLF    groupII 10829169      47  1.0000000 1.0000000
#> 279   pop      OPL    groupII 10829169      47  0.3664596 0.5943411
#> 280   pop      PAL    groupII 10829169      47  0.6374657 0.8321500
#> 281   pop      SMR    groupII 10829169      47  0.4408967 0.6857630
#> 282   pop      UPD    groupII 10829169      47  0.5034913 0.7511574
#> 283   pop      ASP    groupXX 10913115      48  0.3632192 0.6591364
#> 284   pop      CLF    groupXX 10913115      48  0.4000000 0.5937500
#> 285   pop      OPL    groupXX 10913115      48  0.3682566 0.6008287
#> 286   pop      PAL    groupXX 10913115      48  0.3665865 0.6082972
#> 287   pop      SMR    groupXX 10913115      48  0.3660131 0.5976223
#> 288   pop      UPD    groupXX 10913115      48  0.3659674 0.6106579
#> 289   pop      ASP groupXVIII 10974594      49  1.0000000 1.0000000
#> 290   pop      CLF groupXVIII 10974594      49  0.3333333 0.7511574
#> 291   pop      OPL groupXVIII 10974594      49  0.9220949 0.9641444
#> 292   pop      PAL groupXVIII 10974594      49  1.0000000 1.0000000
#> 293   pop      SMR groupXVIII 10974594      49  1.0000000 1.0000000
#> 294   pop      UPD groupXVIII 10974594      49  0.6923077 0.8655509
#> 295   pop      ASP  groupXVII 11256661      50  0.4902530 0.7693844
#> 296   pop      CLF  groupXVII 11256661      50  0.3531250 0.8321500
#> 297   pop      OPL  groupXVII 11256661      50  0.6371842 0.8203171
#> 298   pop      PAL  groupXVII 11256661      50  0.8631153 0.9408278
#> 299   pop      SMR  groupXVII 11256661      50  0.7894737 0.9040063
#> 300   pop      UPD  groupXVII 11256661      50  0.7142857 0.8739458
#> 301   pop      ASP  groupVIII 11438870      51  0.3602431 0.6094107
#> 302   pop      CLF  groupVIII 11438870      51  1.0000000 1.0000000
#> 303   pop      OPL  groupVIII 11438870      51  0.3801858 0.6226344
#> 304   pop      PAL  groupVIII 11438870      51  0.4545455 0.7356738
#> 305   pop      SMR  groupVIII 11438870      51  0.3780751 0.6351349
#> 306   pop      UPD  groupVIII 11438870      51  0.3665865 0.6082972
#> 307   pop      ASP   groupVII 11477223      52  0.3636364 0.6175125
#> 308   pop      CLF   groupVII 11477223      52  0.3333333 0.6296296
#> 309   pop      OPL   groupVII 11477223      52  0.3886693 0.6296296
#> 310   pop      PAL   groupVII 11477223      52  0.3635417 0.5959247
#> 311   pop      SMR   groupVII 11477223      52  0.5623662 0.7741431
#> 312   pop      UPD   groupVII 11477223      52  0.3706294 0.6198849
#> 313   pop      ASP    groupII 11505084      53  0.4222222 0.7184000
#> 314   pop      CLF    groupII 11505084      53  0.4531250 0.6136475
#> 315   pop      OPL    groupII 11505084      53  0.3946977 0.6410760
#> 316   pop      PAL    groupII 11505084      53  1.0000000 1.0000000
#> 317   pop      SMR    groupII 11505084      53  0.5150762 0.7436258
#> 318   pop      UPD    groupII 11505084      53  0.3636364 0.6025029
#> 319   pop      ASP    groupIX 11904651      54  1.0000000 1.0000000
#> 320   pop      CLF    groupIX 11904651      54  0.5302083 0.8559751
#> 321   pop      OPL    groupIX 11904651      54  0.9189266 0.9628117
#> 322   pop      PAL    groupIX 11904651      54  1.0000000 1.0000000
#> 323   pop      SMR    groupIX 11904651      54  0.9076528 0.9581453
#> 324   pop      UPD    groupIX 11904651      54  1.0000000 1.0000000
#> 325   pop      ASP  groupVIII 12108945      55  0.3714286 0.5937500
#> 326   pop      CLF  groupVIII 12108945      55  0.2968750 0.6777344
#> 327   pop      OPL  groupVIII 12108945      55  0.3859649 0.6296296
#> 328   pop      PAL  groupVIII 12108945      55  0.3641026 0.5951407
#> 329   pop      SMR  groupVIII 12108945      55  0.5075755 0.7425162
#> 330   pop      UPD  groupVIII 12108945      55  0.3636364 0.5956024
#> 331   pop      ASP   groupXIX 12485182      56  0.3571429 0.6296296
#> 332   pop      CLF   groupXIX 12485182      56  0.5302083 0.8559751
#> 333   pop      OPL   groupXIX 12485182      56  0.4749226 0.7184000
#> 334   pop      PAL   groupXIX 12485182      56  0.3766226 0.6251616
#> 335   pop      SMR   groupXIX 12485182      56  0.4923393 0.7333207
#> 336   pop      UPD   groupXIX 12485182      56  0.4268430 0.6902241
#> 337   pop      ASP    groupXI 12833733      57  0.6886161 0.8879623
#> 338   pop      CLF    groupXI 12833733      57  1.0000000 1.0000000
#> 339   pop      OPL    groupXI 12833733      57  0.6526316 0.8321500
#> 340   pop      PAL    groupXI 12833733      57  1.0000000 1.0000000
#> 341   pop      SMR    groupXI 12833733      57  0.8797232 0.9470817
#> 342   pop      UPD    groupXI 12833733      57  1.0000000 1.0000000
#> 343   pop      ASP groupXVIII 12937786      58  0.3650794 0.5976223
#> 344   pop      CLF groupXVIII 12937786      58  1.0000000 1.0000000
#> 345   pop      OPL groupXVIII 12937786      58  0.4270697 0.6707266
#> 346   pop      PAL groupXVIII 12937786      58  0.3641827 0.6011990
#> 347   pop      SMR groupXVIII 12937786      58  0.5623213 0.7822583
#> 348   pop      UPD groupXVIII 12937786      58  0.3818182 0.6623181
#> 349   pop      ASP   groupIII 13177027      59  0.3636364 0.5963390
#> 350   pop      CLF   groupIII 13177027      59  0.3333333 0.6296296
#> 351   pop      OPL   groupIII 13177027      59  0.3737016 0.6074900
#> 352   pop      PAL   groupIII 13177027      59  0.7500000 0.8879623
#> 353   pop      SMR   groupIII 13177027      59  0.3656250 0.5945318
#> 354   pop      UPD   groupIII 13177027      59  0.3659674 0.6106579
#> 355   pop      ASP    groupVI 13335577      60  0.3714286 0.6001666
#> 356   pop      CLF    groupVI 13335577      60  0.2968750 0.6777344
#> 357   pop      OPL    groupVI 13335577      60  0.3833992 0.6235005
#> 358   pop      PAL    groupVI 13335577      60  1.0000000 1.0000000
#> 359   pop      SMR    groupVI 13335577      60  0.3859649 0.6296296
#> 360   pop      UPD    groupVI 13335577      60  0.3780751 0.6351349
#> 361   pop      ASP     groupX 13718292      61  1.0000000 1.0000000
#> 362   pop      CLF     groupX 13718292      61       -Inf 0.7511574
#> 363   pop      OPL     groupX 13718292      61  0.8719604 0.9441290
#> 364   pop      PAL     groupX 13718292      61  1.0000000 1.0000000
#> 365   pop      SMR     groupX 13718292      61  0.8000000 0.9083844
#> 366   pop      UPD     groupX 13718292      61  0.6666667 0.8559751
#> 367   pop      ASP    groupVI 13895716      62  0.6000000 0.8321500
#> 368   pop      CLF    groupVI 13895716      62  1.0000000 1.0000000
#> 369   pop      OPL    groupVI 13895716      62  0.4716022 0.7081484
#> 370   pop      PAL    groupVI 13895716      62  1.0000000 1.0000000
#> 371   pop      SMR    groupVI 13895716      62  0.7639707 0.8879623
#> 372   pop      UPD    groupVI 13895716      62  0.6148383 0.8224007
#> 373   pop      ASP   groupXIV 14009990      63  0.3571429 0.6136475
#> 374   pop      CLF   groupXIV 14009990      63  0.4000000 0.6064000
#> 375   pop      OPL   groupXIV 14009990      63  0.4464052 0.6991312
#> 376   pop      PAL   groupXIV 14009990      63  0.4222222 0.7184000
#> 377   pop      SMR   groupXIV 14009990      63  0.3855311 0.6395493
#> 378   pop      UPD   groupXIV 14009990      63  0.3944857 0.6645721
#> 379   pop      ASP   groupXII 14158350      64  0.3686080 0.6367115
#> 380   pop      CLF   groupXII 14158350      64  0.3333333 0.7511574
#> 381   pop      OPL   groupXII 14158350      64  0.3946977 0.6410760
#> 382   pop      PAL   groupXII 14158350      64  0.3646354 0.6001666
#> 383   pop      SMR   groupXII 14158350      64  0.3786284 0.6175125
#> 384   pop      UPD   groupXII 14158350      64  0.3707990 0.6164463
#> 385   pop      ASP    groupIX 14164215      65  0.4545455 0.7356738
#> 386   pop      CLF    groupIX 14164215      65  0.3802083 0.6025029
#> 387   pop      OPL    groupIX 14164215      65  0.3751290 0.6158409
#> 388   pop      PAL    groupIX 14164215      65  0.4106158 0.6675518
#> 389   pop      SMR    groupIX 14164215      65  0.3845801 0.6263483
#> 390   pop      UPD    groupIX 14164215      65  0.3665865 0.6082972
#> 391   pop      ASP    groupVI 14236496      66  0.3671875 0.5937500
#> 392   pop      CLF    groupVI 14236496      66  0.3802083 0.6025029
#> 393   pop      OPL    groupVI 14236496      66  0.3801858 0.6226344
#> 394   pop      PAL    groupVI 14236496      66  0.3766226 0.6251616
#> 395   pop      SMR    groupVI 14236496      66  0.3709099 0.6113384
#> 396   pop      UPD    groupVI 14236496      66  0.3641827 0.6011990
#> 397   pop      ASP     groupX 14273796      67  1.0000000 1.0000000
#> 398   pop      CLF     groupX 14273796      67  1.0000000 1.0000000
#> 399   pop      OPL     groupX 14273796      67  1.0000000 1.0000000
#> 400   pop      PAL     groupX 14273796      67  0.4036630 0.6662518
#> 401   pop      SMR     groupX 14273796      67  1.0000000 1.0000000
#> 402   pop      UPD     groupX 14273796      67  0.8271780 0.9280870
#> 403   pop      ASP  groupVIII 14445389      68  0.3843006 0.6777344
#> 404   pop      CLF  groupVIII 14445389      68  0.2968750 0.6777344
#> 405   pop      OPL  groupVIII 14445389      68  0.4512002 0.6935089
#> 406   pop      PAL  groupVIII 14445389      68  0.3944857 0.6645721
#> 407   pop      SMR  groupVIII 14445389      68  0.4052288 0.6591364
#> 408   pop      UPD  groupVIII 14445389      68  0.3706294 0.6198849
#> 409   pop      ASP    groupIX 14868533      69  1.0000000 1.0000000
#> 410   pop      CLF    groupIX 14868533      69  0.3531250 0.8321500
#> 411   pop      OPL    groupIX 14868533      69  0.4031929 0.6461500
#> 412   pop      PAL    groupIX 14868533      69  1.0000000 1.0000000
#> 413   pop      SMR    groupIX 14868533      69  0.6666667 0.8388017
#> 414   pop      UPD    groupIX 14868533      69  1.0000000 1.0000000
#> 415   pop      ASP    groupIX 15190693      70  0.3612847 0.6226344
#> 416   pop      CLF    groupIX 15190693      70  0.4000000 0.6064000
#> 417   pop      OPL    groupIX 15190693      70  0.4204496 0.6694322
#> 418   pop      PAL    groupIX 15190693      70  1.0000000 1.0000000
#> 419   pop      SMR    groupIX 15190693      70  0.4306156 0.6777344
#> 420   pop      UPD    groupIX 15190693      70  0.3802083 0.6508292
#> 421   pop      ASP   groupXII 15306362      71  0.7901988 0.9160470
#> 422   pop      CLF   groupXII 15306362      71  1.0000000 1.0000000
#> 423   pop      OPL   groupXII 15306362      71  0.9220949 0.9641444
#> 424   pop      PAL   groupXII 15306362      71  0.4785714 0.7301254
#> 425   pop      SMR   groupXII 15306362      71  0.6027668 0.7991943
#> 426   pop      UPD   groupXII 15306362      71  0.7333333 0.8813630
#> 427   pop      ASP    groupXX 15460430      72  0.3612847 0.6226344
#> 428   pop      CLF    groupXX 15460430      72  0.3802083 0.6025029
#> 429   pop      OPL    groupXX 15460430      72  0.3965936 0.6395493
#> 430   pop      PAL    groupXX 15460430      72  0.3707990 0.6164463
#> 431   pop      SMR    groupXX 15460430      72  0.4499886 0.6950238
#> 432   pop      UPD    groupXX 15460430      72  0.3944857 0.6645721
#> 433   pop      ASP   groupVII 15493025      73  0.3666667 0.6461500
#> 434   pop      CLF   groupVII 15493025      73  0.2968750 0.6777344
#> 435   pop      OPL   groupVII 15493025      73  0.5900621 0.7924536
#> 436   pop      PAL   groupVII 15493025      73  0.3877404 0.6461500
#> 437   pop      SMR   groupVII 15493025      73  0.6102365 0.8067477
#> 438   pop      UPD   groupVII 15493025      73  0.4102928 0.6777344
#> 439   pop      ASP    groupXV 15547073      74  0.4902530 0.7693844
#> 440   pop      CLF    groupXV 15547073      74  0.3531250 0.8321500
#> 441   pop      OPL    groupXV 15547073      74  0.3895769 0.6378941
#> 442   pop      PAL    groupXV 15547073      74  1.0000000 1.0000000
#> 443   pop      SMR    groupXV 15547073      74  0.4319328 0.6887010
#> 444   pop      UPD    groupXV 15547073      74  0.5893575 0.8114857
#> 445   pop      ASP  groupXIII 15697285      75  0.3650794 0.5976223
#> 446   pop      CLF  groupXIII 15697285      75  0.4531250 0.5937500
#> 447   pop      OPL  groupXIII 15697285      75  0.3687182 0.5990938
#> 448   pop      PAL  groupXIII 15697285      75  0.3855311 0.6395493
#> 449   pop      SMR  groupXIII 15697285      75  0.3656250 0.5945318
#> 450   pop      UPD  groupXIII 15697285      75  0.3877404 0.6461500
#> 451   pop      ASP   groupIII 15983994      76  0.3802083 0.6025029
#> 452   pop      CLF   groupIII 15983994      76       -Inf 0.7511574
#> 453   pop      OPL   groupIII 15983994      76  0.3766226 0.6251616
#> 454   pop      PAL   groupIII 15983994      76  1.0000000 1.0000000
#> 455   pop      SMR   groupIII 15983994      76  0.5605114 0.7991943
#> 456   pop      UPD   groupIII 15983994      76  0.3646354 0.6001666
#> 457   pop      ASP   groupXIX 16233699      77  0.3571429 0.6296296
#> 458   pop      CLF   groupXIX 16233699      77  0.5302083 0.8559751
#> 459   pop      OPL   groupXIX 16233699      77  0.5150762 0.7436258
#> 460   pop      PAL   groupXIX 16233699      77  0.3707990 0.6164463
#> 461   pop      SMR   groupXIX 16233699      77  0.4923393 0.7333207
#> 462   pop      UPD   groupXIX 16233699      77  0.4395604 0.7057867
#> 463   pop      ASP    groupIX 16306966      78  1.0000000 1.0000000
#> 464   pop      CLF    groupIX 16306966      78  0.2968750 0.6777344
#> 465   pop      OPL    groupIX 16306966      78  0.4023062 0.6486743
#> 466   pop      PAL    groupIX 16306966      78  1.0000000 1.0000000
#> 467   pop      SMR    groupIX 16306966      78  0.6241563 0.8137715
#> 468   pop      UPD    groupIX 16306966      78  1.0000000 1.0000000
#> 469   pop      ASP     groupI 16498302      79  0.4545455 0.7356738
#> 470   pop      CLF     groupI 16498302      79  1.0000000 1.0000000
#> 471   pop      OPL     groupI 16498302      79  0.6261538 0.8114857
#> 472   pop      PAL     groupI 16498302      79  0.4319328 0.6887010
#> 473   pop      SMR     groupI 16498302      79  0.6796537 0.8449559
#> 474   pop      UPD     groupI 16498302      79  0.5605114 0.7991943
#> 475   pop      ASP  groupXIII 16520929      80  0.4083807 0.6950238
#> 476   pop      CLF  groupXIII 16520929      80  0.3333333 0.7511574
#> 477   pop      OPL  groupXIII 16520929      80  0.3923445 0.6367115
#> 478   pop      PAL  groupXIII 16520929      80  0.3638393 0.5937500
#> 479   pop      SMR  groupXIII 16520929      80  0.4240542 0.6701341
#> 480   pop      UPD  groupXIII 16520929      80  0.3944857 0.6645721
#> 481   pop      ASP   groupXVI 16561891      81  0.3571429 0.6296296
#> 482   pop      CLF   groupXVI 16561891      81  0.4000000 0.6064000
#> 483   pop      OPL   groupXVI 16561891      81  0.3708882 0.6052113
#> 484   pop      PAL   groupXVI 16561891      81  0.3653846 0.5986557
#> 485   pop      SMR   groupXVI 16561891      81  0.3923445 0.6367115
#> 486   pop      UPD   groupXVI 16561891      81  0.3665865 0.6082972
#> 487   pop      ASP    groupII 16598364      82  0.6000000 0.8321500
#> 488   pop      CLF    groupII 16598364      82        NaN 1.0000000
#> 489   pop      OPL    groupII 16598364      82  0.6716189 0.8375130
#> 490   pop      PAL    groupII 16598364      82  1.0000000 1.0000000
#> 491   pop      SMR    groupII 16598364      82  0.7894737 0.9040063
#> 492   pop      UPD    groupII 16598364      82  0.4577142 0.7255381
#> 493   pop      ASP    groupIV 17248473      83  0.4545455 0.7356738
#> 494   pop      CLF    groupIV 17248473      83  0.3802083 0.6025029
#> 495   pop      OPL    groupIV 17248473      83  1.0000000 1.0000000
#> 496   pop      PAL    groupIV 17248473      83  1.0000000 1.0000000
#> 497   pop      SMR    groupIV 17248473      83  0.9117451 0.9598256
#> 498   pop      UPD    groupIV 17248473      83  0.8411120 0.9329015
#> 499   pop      ASP    groupXX 17748944      84  0.3666667 0.6461500
#> 500   pop      CLF    groupXX 17748944      84  0.3333333 0.6296296
#> 501   pop      OPL    groupXX 17748944      84  0.3801858 0.6226344
#> 502   pop      PAL    groupXX 17748944      84  1.0000000 1.0000000
#> 503   pop      SMR    groupXX 17748944      84  0.4160647 0.6685876
#> 504   pop      UPD    groupXX 17748944      84  0.3636364 0.5956024
#> 505   pop      ASP  groupXIII 17842228      85  0.3612847 0.6226344
#> 506   pop      CLF  groupXIII 17842228      85  0.4000000 0.6064000
#> 507   pop      OPL  groupXIII 17842228      85  0.3946977 0.6410760
#> 508   pop      PAL  groupXIII 17842228      85  0.6148383 0.8224007
#> 509   pop      SMR  groupXIII 17842228      85  0.4204496 0.6694322
#> 510   pop      UPD  groupXIII 17842228      85  0.3650298 0.5993352
#> 511   pop      ASP   groupXIX 17858769      86  0.3571429 0.6136475
#> 512   pop      CLF   groupXIX 17858769      86  0.5302083 0.8559751
#> 513   pop      OPL   groupXIX 17858769      86  0.4484423 0.6968600
#> 514   pop      PAL   groupXIX 17858769      86  0.3629261 0.6041830
#> 515   pop      SMR   groupXIX 17858769      86  0.4923393 0.7333207
#> 516   pop      UPD   groupXIX 17858769      86  0.3944857 0.6645721
#> 517   pop      ASP    groupIX 17929225      87  0.3640625 0.5937500
#> 518   pop      CLF    groupIX 17929225      87  0.5302083 0.8559751
#> 519   pop      OPL    groupIX 17929225      87  0.4427668 0.6851046
#> 520   pop      PAL    groupIX 17929225      87  1.0000000 1.0000000
#> 521   pop      SMR    groupIX 17929225      87  0.4149009 0.6623181
#> 522   pop      UPD    groupIX 17929225      87  0.3646354 0.6001666
#> 523   pop      ASP    groupIX 18339767      88  0.3714286 0.5937500
#> 524   pop      CLF    groupIX 18339767      88  0.2968750 0.6777344
#> 525   pop      OPL    groupIX 18339767      88  0.3876064 0.6334381
#> 526   pop      PAL    groupIX 18339767      88  1.0000000 1.0000000
#> 527   pop      SMR    groupIX 18339767      88  0.3673735 0.6048522
#> 528   pop      UPD    groupIX 18339767      88  0.3696970 0.6296296
#> 529   pop      ASP  groupXIII 18836558      89  0.3571429 0.6296296
#> 530   pop      CLF  groupXIII 18836558      89       -Inf 0.6296296
#> 531   pop      OPL  groupXIII 18836558      89  0.3755566 0.6184208
#> 532   pop      PAL  groupXIII 18836558      89  0.3636364 0.5937500
#> 533   pop      SMR  groupXIII 18836558      89  0.3755566 0.6184208
#> 534   pop      UPD  groupXIII 18836558      89  0.3636364 0.6025029
#> 535   pop      ASP    groupII 19244005      90  0.5276989 0.7852670
#> 536   pop      CLF    groupII 19244005      90  1.0000000 1.0000000
#> 537   pop      OPL    groupII 19244005      90  1.0000000 1.0000000
#> 538   pop      PAL    groupII 19244005      90  1.0000000 1.0000000
#> 539   pop      SMR    groupII 19244005      90  0.9076528 0.9581453
#> 540   pop      UPD    groupII 19244005      90  0.6923077 0.8655509
#> 541   pop      ASP    groupXX 19317273      91  0.4475446 0.7511574
#> 542   pop      CLF    groupXX 19317273      91  0.4531250 0.6136475
#> 543   pop      OPL    groupXX 19317273      91  0.3656250 0.5937500
#> 544   pop      PAL    groupXX 19317273      91  0.6576666 0.8409057
#> 545   pop      SMR    groupXX 19317273      91  0.3665635 0.5968844
#> 546   pop      UPD    groupXX 19317273      91  0.3706294 0.6198849
#> 547   pop      ASP     groupI 19759331      92  0.3802083 0.5937500
#> 548   pop      CLF     groupI 19759331      92  0.2000000 0.7184000
#> 549   pop      OPL     groupI 19759331      92  0.3706294 0.6198849
#> 550   pop      PAL     groupI 19759331      92  0.3636364 0.6175125
#> 551   pop      SMR     groupI 19759331      92  0.4545455 0.7356738
#> 552   pop      UPD     groupI 19759331      92  0.4222222 0.7184000
#> 553   pop      ASP   groupVII 20165240      93  0.5276989 0.7852670
#> 554   pop      CLF   groupVII 20165240      93  0.2968750 0.6777344
#> 555   pop      OPL   groupVII 20165240      93  0.7810671 0.8957053
#> 556   pop      PAL   groupVII 20165240      93  0.5384615 0.7775927
#> 557   pop      SMR   groupVII 20165240      93  0.4923393 0.7333207
#> 558   pop      UPD   groupVII 20165240      93  0.4102928 0.6777344
#> 559   pop      ASP    groupII 20496069      94  1.0000000 1.0000000
#> 560   pop      CLF    groupII 20496069      94  1.0000000 1.0000000
#> 561   pop      OPL    groupII 20496069      94  0.7440127 0.8789886
#> 562   pop      PAL    groupII 20496069      94  0.3703016 0.6241761
#> 563   pop      SMR    groupII 20496069      94  0.9031606 0.9563185
#> 564   pop      UPD    groupII 20496069      94  0.6148383 0.8224007
#> 565   pop      ASP     groupI 21284266      95  0.4222222 0.7184000
#> 566   pop      CLF     groupI 21284266      95  0.3333333 0.6296296
#> 567   pop      OPL     groupI 21284266      95  0.5021189 0.7356738
#> 568   pop      PAL     groupI 21284266      95  1.0000000 1.0000000
#> 569   pop      SMR     groupI 21284266      95  0.5314757 0.7605912
#> 570   pop      UPD     groupI 21284266      95  0.4848485 0.7511574
#> 571   pop      ASP    groupIV 22227311      96  0.6886161 0.8879623
#> 572   pop      CLF    groupIV 22227311      96  0.5302083 0.8559751
#> 573   pop      OPL    groupIV 22227311      96  0.5150762 0.7436258
#> 574   pop      PAL    groupIV 22227311      96  1.0000000 1.0000000
#> 575   pop      SMR    groupIV 22227311      96  0.6916996 0.8506652
#> 576   pop      UPD    groupIV 22227311      96  0.6923077 0.8655509
#> 577   pop      ASP     groupI 22746027      97  1.0000000 1.0000000
#> 578   pop      CLF     groupI 22746027      97  0.5302083 0.8559751
#> 579   pop      OPL     groupI 22746027      97  0.8095238 0.9123822
#> 580   pop      PAL     groupI 22746027      97  1.0000000 1.0000000
#> 581   pop      SMR     groupI 22746027      97  0.8000000 0.9083844
#> 582   pop      UPD     groupI 22746027      97  0.8529419 0.9371130
#> 583   pop      ASP   groupVII 25239495      98  0.3888889 0.6991312
#> 584   pop      CLF   groupVII 25239495      98  0.2781250 0.6461500
#> 585   pop      OPL   groupVII 25239495      98  0.6716189 0.8375130
#> 586   pop      PAL   groupVII 25239495      98  0.3673735 0.6048522
#> 587   pop      SMR   groupVII 25239495      98  0.4736095 0.7142969
#> 588   pop      UPD   groupVII 25239495      98  0.3703016 0.6241761
#> 589   pop      ASP   groupVII 27599144      99  0.3640625 0.5968844
#> 590   pop      CLF   groupVII 27599144      99  0.3802083 0.6025029
#> 591   pop      OPL   groupVII 27599144      99  0.3684211 0.6001666
#> 592   pop      PAL   groupVII 27599144      99  1.0000000 1.0000000
#> 593   pop      SMR   groupVII 27599144      99  0.3656250 0.5945318
#> 594   pop      UPD   groupVII 27599144      99  0.3640625 0.5937500
#> 595   pop      ASP    groupIV 32321895     100  0.3619048 0.6064000
#> 596   pop      CLF    groupIV 32321895     100  0.4531250 0.6136475
#> 597   pop      OPL    groupIV 32321895     100  0.3975232 0.6461500
#> 598   pop      PAL    groupIV 32321895     100  0.3686080 0.6367115
#> 599   pop      SMR    groupIV 32321895     100  0.3876064 0.6334381
#> 600   pop      UPD    groupIV 32321895     100  0.3665865 0.6082972
#> 
#> $weighted.means
#>   facet subfacet snp.facet snp.subfacet          PID       PIDsib
#> 1   pop      ASP     .base        .base 1.018934e-32 1.150169e-14
#> 2   pop      CLF     .base        .base 3.867094e-33 7.344235e-13
#> 3   pop      OPL     .base        .base 3.305705e-30 1.046881e-14
#> 4   pop      PAL     .base        .base 4.247995e-19 2.781777e-09
#> 5   pop      SMR     .base        .base 2.118404e-29 4.263310e-14
#> 6   pop      UPD     .base        .base 4.142719e-33 1.471158e-15
#> 

```
