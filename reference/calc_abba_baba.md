# Conduct an ABBA/BABA test for gene flow.

Estimate D according to Green et al (2010) via an ABBA/BABA test for a
specific set of three different populations.

## Usage

``` r
calc_abba_baba(
  x,
  facet,
  p1,
  p2,
  p3,
  jackknife = FALSE,
  jackknife_par = FALSE,
  sigma = NULL
)
```

## Arguments

- x:

  snpRdata. Input SNP data.

- facet:

  character. Categorical metadata variables by which to break up
  analysis. Must contain a sample facet. See
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md)
  for more details.

- p1:

  character. Name of population 1, must match a category present in the
  provided facet.

- p2:

  character. Name of population 2, must match a category present in the
  provided facet.

- p3:

  character. Name of population 3, must match a category present in the
  provided facet.

- jackknife:

  logical, default FALSE. If TRUE, block-jackknifed significance for D
  will be calculated with window size sigma according to any SNP facet
  levels. See details.

- jackknife_par:

  numeric or FALSE, default FALSE. If numeric, jackknifes per SNP levels
  will be run with the requested number of processing threads.

- sigma:

  numeric or NULL, default NULL. If jackknifes are requested, the size
  of the windows to block by, in kb. Should be large enough to account
  for linkage!

## Details

An ABBA/BABA test according to Green et al (2010) tests the hypothesis
that there are an equal number of loci where population 1 and population
2 are more closely related to population 3. The ratio of these two
scenarios is given as ABBA/BABA, where: ABBA = (1 - p1)p2p3 and BABA =
p1(1 - p2)p3. where p1, p2, and p3 are the derived allele frequencies in
populations 1 through 3, respectively. *D* values are provided for both
the overall comparison and within any levels of provided snp facets.

*p*-values for *D* can be calculated by a block jackknifing approach
according to Maier et. al (2022) by removing SNPs from each of *n*
genomic windows (blocks) and then calculating *D* for all other sites.
Non-overlapping windows are "blocked" by reference to any provided SNP
metadata facets (usually chromosome), each with a length equal to
*sigma*, so providing *sigma = 100* will use 100kb windows.

## References

Maier, R. et al (2022). On the limits of fitting complex models of
population history to genetic data. BioRxiv. doi:
10.1101/2022.05.08.491072

Green, ER, et al (2010). A Draft Sequence of the Neandertal Genome.
Science, 328(5979), 710–722. doi: 10.1126/science.1188021

## Author

William Hemstrom

## Examples

``` r

# add the ref and anc columns
# note: here these results are meaningless since they are arbitrary.
x <- stickSNPs
maf <- get.snpR.stats(x, ".base", stats = "maf")$single
snp.meta(x)$ref <- maf$major
#> Warning: Some levels are duplicated across multiple sample meta facets.
#> This will cause issues if those sample facets are run during analysis.
#> Issues:
#> 
#> Level: A in facets: fam, ref
#> Level: C in facets: fam, ref
snp.meta(x)$anc <- maf$minor
#> Warning: Some levels are duplicated across multiple sample meta facets.
#> This will cause issues if those sample facets are run during analysis.
#> Issues:
#> 
#> Level: A in facets: fam, ref, anc
#> Level: C in facets: fam, ref, anc
#> Level: G in facets: ref, anc
#> Level: T in facets: ref, anc

# run with jackknifing, 1000kb windows!
# Test if ASP or UPD have more geneflow with PAL...
x <- calc_abba_baba(x, "pop.chr", "ASP", "UPD", "PAL", TRUE, sigma = 1000)
get.snpR.stats(x, "pop.chr", "abba_baba") # gets the per chr results
#> $pairwise
#>     facet  comparison        chr position ref anc .snp.id  D_abba_baba
#> 1     pop ASP~UPD~PAL    groupVI   212436   G   A       1  0.000000000
#> 2     pop ASP~UPD~PAL groupXVIII   350070   G   T       2  1.000000000
#> 3     pop ASP~UPD~PAL    groupXV   854285   C   T       3  0.166666667
#> 4     pop ASP~UPD~PAL  groupXIII   907279   A   G       4 -0.645669291
#> 5     pop ASP~UPD~PAL    groupVI  1068620   G   A       5 -0.180327869
#> 6     pop ASP~UPD~PAL     groupV  1392753   A   C       6 -0.250000000
#> 7     pop ASP~UPD~PAL    groupIX  2075079   T   G       7  1.000000000
#> 8     pop ASP~UPD~PAL  groupXIII  2106000   G   A       8 -0.062500000
#> 9     pop ASP~UPD~PAL   groupXII  2136587   C   A       9 -0.112149533
#> 10    pop ASP~UPD~PAL  groupXVII  2560839   A   G      10 -0.250000000
#> 11    pop ASP~UPD~PAL   groupIII  2750695   T   C      11  0.116788321
#> 12    pop ASP~UPD~PAL   groupXIX  3008524   A   G      12 -0.125000000
#> 13    pop ASP~UPD~PAL   groupXVI  3624417   T   C      13 -0.185185185
#> 14    pop ASP~UPD~PAL   groupVII  3667528   G   T      14 -0.396825397
#> 15    pop ASP~UPD~PAL   groupVII  4189120   A   G      15 -0.489051095
#> 16    pop ASP~UPD~PAL    groupXX  4538074   G   T      16  0.117647059
#> 17    pop ASP~UPD~PAL  groupXIII  4706089   A   C      17 -0.014925373
#> 18    pop ASP~UPD~PAL   groupIII  4746808   C   A      18 -0.281690141
#> 19    pop ASP~UPD~PAL  groupXVII  5041995   T   G      19 -1.000000000
#> 20    pop ASP~UPD~PAL     groupX  5161637   G   A      20 -0.133333333
#> 21    pop ASP~UPD~PAL   groupXII  5178037   C   A      21 -0.428571429
#> 22    pop ASP~UPD~PAL    groupIX  5287711   C   T      22  0.180327869
#> 23    pop ASP~UPD~PAL  groupVIII  5462135   G   A      23  0.208333333
#> 24    pop ASP~UPD~PAL    groupIX  5960933   A   G      24 -0.390243902
#> 25    pop ASP~UPD~PAL  groupXIII  6076868   G   A      25 -0.281690141
#> 26    pop ASP~UPD~PAL    groupXI  6177332   T   C      26  0.096774194
#> 27    pop ASP~UPD~PAL    groupIV  6291287   A   C      27 -0.187500000
#> 28    pop ASP~UPD~PAL    groupIX  6515796   T   C      28 -0.473684211
#> 29    pop ASP~UPD~PAL   groupXII  6718079   A   C      29 -0.378378378
#> 30    pop ASP~UPD~PAL   groupVII  6967477   G   C      30 -0.333333333
#> 31    pop ASP~UPD~PAL   groupIII  7194069   C   G      31 -0.341935484
#> 32    pop ASP~UPD~PAL   groupXVI  7199722   A   C      32 -0.352272727
#> 33    pop ASP~UPD~PAL  groupVIII  7224065   C   T      33 -0.170731707
#> 34    pop ASP~UPD~PAL groupXVIII  7229426   C   G      34  0.602409639
#> 35    pop ASP~UPD~PAL groupXVIII  8087022   A   C      35  0.181818182
#> 36    pop ASP~UPD~PAL    groupIX  8373442   A   G      36  0.358490566
#> 37    pop ASP~UPD~PAL    groupXV  8815538   G   C      37 -0.191489362
#> 38    pop ASP~UPD~PAL    groupXI  9019877   G   T      38  0.421052632
#> 39    pop ASP~UPD~PAL   groupXXI  9067876   T   A      39  0.062500000
#> 40    pop ASP~UPD~PAL  groupXVII  9108184   G   C      40  0.083333333
#> 41    pop ASP~UPD~PAL    groupIX  9284427   G   A      41          NaN
#> 42    pop ASP~UPD~PAL     groupV  9767769   G   A      42 -0.250000000
#> 43    pop ASP~UPD~PAL   groupXII  9881312   A   G      43 -0.500000000
#> 44    pop ASP~UPD~PAL    groupXX 10296516   T   A      44 -0.136363636
#> 45    pop ASP~UPD~PAL     groupI 10593829   G   T      45 -1.000000000
#> 46    pop ASP~UPD~PAL   groupXIX 10765968   C   A      46  0.190000000
#> 47    pop ASP~UPD~PAL    groupII 10829169   C   G      47  0.190476190
#> 48    pop ASP~UPD~PAL    groupXX 10913115   T   C      48 -0.238095238
#> 49    pop ASP~UPD~PAL groupXVIII 10974594   G   C      49 -1.000000000
#> 50    pop ASP~UPD~PAL  groupXVII 11256661   C   T      50  0.392857143
#> 51    pop ASP~UPD~PAL  groupVIII 11438870   T   C      51 -0.008333333
#> 52    pop ASP~UPD~PAL   groupVII 11477223   A   G      52  0.014084507
#> 53    pop ASP~UPD~PAL    groupII 11505084   A   C      53 -0.481481481
#> 54    pop ASP~UPD~PAL    groupIX 11904651   G   T      54          NaN
#> 55    pop ASP~UPD~PAL  groupVIII 12108945   G   A      55  0.076923077
#> 56    pop ASP~UPD~PAL   groupXIX 12485182   A   G      56  0.243243243
#> 57    pop ASP~UPD~PAL    groupXI 12833733   G   A      57  1.000000000
#> 58    pop ASP~UPD~PAL groupXVIII 12937786   C   T      58  0.361702128
#> 59    pop ASP~UPD~PAL   groupIII 13177027   G   A      59  0.142857143
#> 60    pop ASP~UPD~PAL    groupVI 13335577   A   C      60  0.225806452
#> 61    pop ASP~UPD~PAL     groupX 13718292   T   C      61 -1.000000000
#> 62    pop ASP~UPD~PAL    groupVI 13895716   G   C      62 -0.038461538
#> 63    pop ASP~UPD~PAL   groupXIV 14009990   A   C      63  0.239130435
#> 64    pop ASP~UPD~PAL   groupXII 14158350   T   A      64 -0.107382550
#> 65    pop ASP~UPD~PAL    groupIX 14164215   C   T      65 -0.488721805
#> 66    pop ASP~UPD~PAL    groupVI 14236496   A   C      66  0.153846154
#> 67    pop ASP~UPD~PAL     groupX 14273796   G   A      67 -1.000000000
#> 68    pop ASP~UPD~PAL  groupVIII 14445389   A   T      68 -0.250000000
#> 69    pop ASP~UPD~PAL    groupIX 14868533   A   G      69          NaN
#> 70    pop ASP~UPD~PAL    groupIX 15190693   A   C      70  0.133333333
#> 71    pop ASP~UPD~PAL   groupXII 15306362   A   T      71 -0.200000000
#> 72    pop ASP~UPD~PAL    groupXX 15460430   A   G      72 -0.187500000
#> 73    pop ASP~UPD~PAL   groupVII 15493025   T   G      73  0.125000000
#> 74    pop ASP~UPD~PAL    groupXV 15547073   G   T      74  0.150000000
#> 75    pop ASP~UPD~PAL  groupXIII 15697285   A   C      75 -0.302325581
#> 76    pop ASP~UPD~PAL   groupIII 15983994   C   T      76          NaN
#> 77    pop ASP~UPD~PAL   groupXIX 16233699   A   T      77  0.294117647
#> 78    pop ASP~UPD~PAL    groupIX 16306966   G   T      78          NaN
#> 79    pop ASP~UPD~PAL     groupI 16498302   T   C      79  0.217391304
#> 80    pop ASP~UPD~PAL  groupXIII 16520929   T   C      80 -0.112149533
#> 81    pop ASP~UPD~PAL   groupXVI 16561891   A   T      81 -0.128205128
#> 82    pop ASP~UPD~PAL    groupII 16598364   T   C      82 -0.363636364
#> 83    pop ASP~UPD~PAL    groupIV 17248473   T   A      83  0.714285714
#> 84    pop ASP~UPD~PAL    groupXX 17748944   C   T      84          NaN
#> 85    pop ASP~UPD~PAL  groupXIII 17842228   G   A      85  0.416666667
#> 86    pop ASP~UPD~PAL   groupXIX 17858769   A   G      86  0.239130435
#> 87    pop ASP~UPD~PAL    groupIX 17929225   T   C      87  0.142857143
#> 88    pop ASP~UPD~PAL    groupIX 18339767   A   G      88          NaN
#> 89    pop ASP~UPD~PAL  groupXIII 18836558   C   G      89 -0.473684211
#> 90    pop ASP~UPD~PAL    groupII 19244005   C   T      90  0.309090909
#> 91    pop ASP~UPD~PAL    groupXX 19317273   C   T      91 -0.470588235
#> 92    pop ASP~UPD~PAL     groupI 19759331   G   A      92  0.600000000
#> 93    pop ASP~UPD~PAL   groupVII 20165240   T   C      93 -0.357142857
#> 94    pop ASP~UPD~PAL    groupII 20496069   C   G      94 -1.000000000
#> 95    pop ASP~UPD~PAL     groupI 21284266   G   T      95  0.111111111
#> 96    pop ASP~UPD~PAL    groupIV 22227311   C   A      96 -0.111111111
#> 97    pop ASP~UPD~PAL     groupI 22746027   G   A      97 -1.000000000
#> 98    pop ASP~UPD~PAL   groupVII 25239495   G   T      98 -0.298969072
#> 99    pop ASP~UPD~PAL   groupVII 27599144   G   A      99  0.100000000
#> 100   pop ASP~UPD~PAL    groupIV 32321895   A   T     100 -0.014925373
#>           abba       baba nk
#> 1   0.19852941 0.19852941 78
#> 2   0.05000000 0.00000000 84
#> 3   0.08219697 0.05871212 84
#> 4   0.07867133 0.36538462 80
#> 5   0.15707672 0.22619048 56
#> 6   0.09542411 0.15904018 80
#> 7   0.04545455 0.00000000 84
#> 8   0.17045455 0.19318182 70
#> 9   0.16089379 0.20154065 80
#> 10  0.07792208 0.12987013 84
#> 11  0.27321429 0.21607143 80
#> 12  0.12031250 0.15468750 72
#> 13  0.09088319 0.13219373 74
#> 14  0.12666667 0.29333333 84
#> 15  0.01529720 0.04458042 80
#> 16  0.11215278 0.08854167 74
#> 17  0.10312500 0.10625000 78
#> 18  0.04857143 0.08666667 78
#> 19  0.00000000 0.06666667 86
#> 20  0.05158730 0.06746032 76
#> 21  0.11111111 0.27777778 74
#> 22  0.28571429 0.19841270 80
#> 23  0.04833333 0.03166667 82
#> 24  0.04464286 0.10178571 80
#> 25  0.17000000 0.30333333 82
#> 26  0.11018519 0.09074074 72
#> 27  0.04642857 0.06785714 82
#> 28  0.08035714 0.22500000 74
#> 29  0.04095442 0.09081197 74
#> 30  0.07575758 0.15151515 80
#> 31  0.15454545 0.31515152 86
#> 32  0.06576923 0.13730769 76
#> 33  0.13671329 0.19300699 78
#> 34  0.33585859 0.08333333 74
#> 35  0.09285714 0.06428571 80
#> 36  0.13846154 0.06538462 78
#> 37  0.14158654 0.20865385 78
#> 38  0.20039062 0.08164063 80
#> 39  0.12234848 0.10795455 60
#> 40  0.22159091 0.18750000 72
#> 41  0.00000000 0.00000000 68
#> 42  0.15066964 0.25111607 76
#> 43  0.06993007 0.20979021 82
#> 44  0.07916667 0.10416667 72
#> 45  0.00000000 0.06666667 86
#> 46  0.28605769 0.19471154 74
#> 47  0.17045455 0.11590909 82
#> 48  0.10378510 0.16865079 72
#> 49  0.00000000 0.07692308 80
#> 50  0.13493304 0.05881696 80
#> 51  0.19318182 0.19642857 68
#> 52  0.12662338 0.12310606 74
#> 53  0.11666667 0.33333333 76
#> 54  0.00000000 0.00000000 84
#> 55  0.14358974 0.12307692 72
#> 56  0.16770833 0.10208333 80
#> 57  0.06250000 0.00000000 60
#> 58  0.18648019 0.08741259 66
#> 59  0.26223776 0.19667832 80
#> 60  0.29081633 0.18367347 76
#> 61  0.00000000 0.08333333 66
#> 62  0.08928571 0.09642857 80
#> 63  0.05480769 0.03365385 62
#> 64  0.08636364 0.10714286 80
#> 65  0.08116883 0.23634454 84
#> 66  0.09915865 0.07271635 76
#> 67  0.00000000 0.02820513 74
#> 68  0.04326923 0.07211538 74
#> 69  0.00000000 0.00000000 84
#> 70  0.24791667 0.18958333 76
#> 71  0.03446970 0.05170455 84
#> 72  0.06416667 0.09378205 76
#> 73  0.15750000 0.12250000 78
#> 74  0.13269231 0.09807692 80
#> 75  0.11458333 0.21388889 80
#> 76  0.00000000 0.00000000 56
#> 77  0.16587302 0.09047619 76
#> 78  0.00000000 0.00000000 84
#> 79  0.12165775 0.07820856 80
#> 80  0.08304196 0.10402098 78
#> 81  0.11383929 0.14732143 78
#> 82  0.08076923 0.17307692 80
#> 83  0.17532468 0.02922078 82
#> 84  0.00000000 0.00000000 76
#> 85  0.03946429 0.01625000 78
#> 86  0.16193182 0.09943182 64
#> 87  0.28571429 0.21428571 82
#> 88  0.00000000 0.00000000 68
#> 89  0.06944444 0.19444444 66
#> 90  0.12587413 0.06643357 82
#> 91  0.09709821 0.26971726 78
#> 92  0.14545455 0.03636364 56
#> 93  0.08766234 0.18506494 78
#> 94  0.00000000 0.07005495 70
#> 95  0.16666667 0.13333333 78
#> 96  0.05769231 0.07211538 76
#> 97  0.00000000 0.03333333 86
#> 98  0.08627137 0.15985577 76
#> 99  0.27500000 0.22500000 70
#> 100 0.16071429 0.16558442 70
#> 
#> $weighted.means
#>    facet    subfacet snp.facet snp.subfacet D_abba_baba_jackknife
#> 1    pop ASP~UPD~PAL       chr       groupI                    NA
#> 2    pop ASP~UPD~PAL       chr      groupII                    NA
#> 3    pop ASP~UPD~PAL       chr     groupIII                    NA
#> 4    pop ASP~UPD~PAL       chr      groupIV                    NA
#> 5    pop ASP~UPD~PAL       chr      groupIX                    NA
#> 6    pop ASP~UPD~PAL       chr       groupV                    NA
#> 7    pop ASP~UPD~PAL       chr      groupVI                    NA
#> 8    pop ASP~UPD~PAL       chr     groupVII                    NA
#> 9    pop ASP~UPD~PAL       chr    groupVIII                    NA
#> 10   pop ASP~UPD~PAL       chr       groupX                    NA
#> 11   pop ASP~UPD~PAL       chr      groupXI                    NA
#> 12   pop ASP~UPD~PAL       chr     groupXII                    NA
#> 13   pop ASP~UPD~PAL       chr    groupXIII                    NA
#> 14   pop ASP~UPD~PAL       chr     groupXIV                    NA
#> 15   pop ASP~UPD~PAL       chr     groupXIX                    NA
#> 16   pop ASP~UPD~PAL       chr      groupXV                    NA
#> 17   pop ASP~UPD~PAL       chr     groupXVI                    NA
#> 18   pop ASP~UPD~PAL       chr    groupXVII                    NA
#> 19   pop ASP~UPD~PAL       chr   groupXVIII                    NA
#> 20   pop ASP~UPD~PAL       chr      groupXX                    NA
#> 21   pop ASP~UPD~PAL       chr     groupXXI                    NA
#>    D_abba_baba_p_value  D_abba_baba
#> 1                   NA  0.109856907
#> 2                   NA -0.211599173
#> 3                   NA -0.048932585
#> 4                   NA  0.135987882
#> 5                   NA -0.008755932
#> 6                   NA -0.250000000
#> 7                   NA  0.035554682
#> 8                   NA -0.157011913
#> 9                   NA -0.043345106
#> 10                  NA -0.552554938
#> 11                  NA  0.367938040
#> 12                  NA -0.301593120
#> 13                  NA -0.287231988
#> 14                  NA  0.239130435
#> 15                  NA  0.168792465
#> 16                  NA -0.012421169
#> 17                  NA -0.212902729
#> 18                  NA -0.009583655
#> 19                  NA  0.361501281
#> 20                  NA -0.227296543
#> 21                  NA  0.062500000
#> 
get.snpR.stats(x, "pop", "abba_baba") # gets the overall results
#> $pairwise
#>     facet  comparison        chr position ref anc .snp.id  D_abba_baba
#> 1     pop ASP~UPD~PAL    groupVI   212436   G   A       1  0.000000000
#> 2     pop ASP~UPD~PAL groupXVIII   350070   G   T       2  1.000000000
#> 3     pop ASP~UPD~PAL    groupXV   854285   C   T       3  0.166666667
#> 4     pop ASP~UPD~PAL  groupXIII   907279   A   G       4 -0.645669291
#> 5     pop ASP~UPD~PAL    groupVI  1068620   G   A       5 -0.180327869
#> 6     pop ASP~UPD~PAL     groupV  1392753   A   C       6 -0.250000000
#> 7     pop ASP~UPD~PAL    groupIX  2075079   T   G       7  1.000000000
#> 8     pop ASP~UPD~PAL  groupXIII  2106000   G   A       8 -0.062500000
#> 9     pop ASP~UPD~PAL   groupXII  2136587   C   A       9 -0.112149533
#> 10    pop ASP~UPD~PAL  groupXVII  2560839   A   G      10 -0.250000000
#> 11    pop ASP~UPD~PAL   groupIII  2750695   T   C      11  0.116788321
#> 12    pop ASP~UPD~PAL   groupXIX  3008524   A   G      12 -0.125000000
#> 13    pop ASP~UPD~PAL   groupXVI  3624417   T   C      13 -0.185185185
#> 14    pop ASP~UPD~PAL   groupVII  3667528   G   T      14 -0.396825397
#> 15    pop ASP~UPD~PAL   groupVII  4189120   A   G      15 -0.489051095
#> 16    pop ASP~UPD~PAL    groupXX  4538074   G   T      16  0.117647059
#> 17    pop ASP~UPD~PAL  groupXIII  4706089   A   C      17 -0.014925373
#> 18    pop ASP~UPD~PAL   groupIII  4746808   C   A      18 -0.281690141
#> 19    pop ASP~UPD~PAL  groupXVII  5041995   T   G      19 -1.000000000
#> 20    pop ASP~UPD~PAL     groupX  5161637   G   A      20 -0.133333333
#> 21    pop ASP~UPD~PAL   groupXII  5178037   C   A      21 -0.428571429
#> 22    pop ASP~UPD~PAL    groupIX  5287711   C   T      22  0.180327869
#> 23    pop ASP~UPD~PAL  groupVIII  5462135   G   A      23  0.208333333
#> 24    pop ASP~UPD~PAL    groupIX  5960933   A   G      24 -0.390243902
#> 25    pop ASP~UPD~PAL  groupXIII  6076868   G   A      25 -0.281690141
#> 26    pop ASP~UPD~PAL    groupXI  6177332   T   C      26  0.096774194
#> 27    pop ASP~UPD~PAL    groupIV  6291287   A   C      27 -0.187500000
#> 28    pop ASP~UPD~PAL    groupIX  6515796   T   C      28 -0.473684211
#> 29    pop ASP~UPD~PAL   groupXII  6718079   A   C      29 -0.378378378
#> 30    pop ASP~UPD~PAL   groupVII  6967477   G   C      30 -0.333333333
#> 31    pop ASP~UPD~PAL   groupIII  7194069   C   G      31 -0.341935484
#> 32    pop ASP~UPD~PAL   groupXVI  7199722   A   C      32 -0.352272727
#> 33    pop ASP~UPD~PAL  groupVIII  7224065   C   T      33 -0.170731707
#> 34    pop ASP~UPD~PAL groupXVIII  7229426   C   G      34  0.602409639
#> 35    pop ASP~UPD~PAL groupXVIII  8087022   A   C      35  0.181818182
#> 36    pop ASP~UPD~PAL    groupIX  8373442   A   G      36  0.358490566
#> 37    pop ASP~UPD~PAL    groupXV  8815538   G   C      37 -0.191489362
#> 38    pop ASP~UPD~PAL    groupXI  9019877   G   T      38  0.421052632
#> 39    pop ASP~UPD~PAL   groupXXI  9067876   T   A      39  0.062500000
#> 40    pop ASP~UPD~PAL  groupXVII  9108184   G   C      40  0.083333333
#> 41    pop ASP~UPD~PAL    groupIX  9284427   G   A      41          NaN
#> 42    pop ASP~UPD~PAL     groupV  9767769   G   A      42 -0.250000000
#> 43    pop ASP~UPD~PAL   groupXII  9881312   A   G      43 -0.500000000
#> 44    pop ASP~UPD~PAL    groupXX 10296516   T   A      44 -0.136363636
#> 45    pop ASP~UPD~PAL     groupI 10593829   G   T      45 -1.000000000
#> 46    pop ASP~UPD~PAL   groupXIX 10765968   C   A      46  0.190000000
#> 47    pop ASP~UPD~PAL    groupII 10829169   C   G      47  0.190476190
#> 48    pop ASP~UPD~PAL    groupXX 10913115   T   C      48 -0.238095238
#> 49    pop ASP~UPD~PAL groupXVIII 10974594   G   C      49 -1.000000000
#> 50    pop ASP~UPD~PAL  groupXVII 11256661   C   T      50  0.392857143
#> 51    pop ASP~UPD~PAL  groupVIII 11438870   T   C      51 -0.008333333
#> 52    pop ASP~UPD~PAL   groupVII 11477223   A   G      52  0.014084507
#> 53    pop ASP~UPD~PAL    groupII 11505084   A   C      53 -0.481481481
#> 54    pop ASP~UPD~PAL    groupIX 11904651   G   T      54          NaN
#> 55    pop ASP~UPD~PAL  groupVIII 12108945   G   A      55  0.076923077
#> 56    pop ASP~UPD~PAL   groupXIX 12485182   A   G      56  0.243243243
#> 57    pop ASP~UPD~PAL    groupXI 12833733   G   A      57  1.000000000
#> 58    pop ASP~UPD~PAL groupXVIII 12937786   C   T      58  0.361702128
#> 59    pop ASP~UPD~PAL   groupIII 13177027   G   A      59  0.142857143
#> 60    pop ASP~UPD~PAL    groupVI 13335577   A   C      60  0.225806452
#> 61    pop ASP~UPD~PAL     groupX 13718292   T   C      61 -1.000000000
#> 62    pop ASP~UPD~PAL    groupVI 13895716   G   C      62 -0.038461538
#> 63    pop ASP~UPD~PAL   groupXIV 14009990   A   C      63  0.239130435
#> 64    pop ASP~UPD~PAL   groupXII 14158350   T   A      64 -0.107382550
#> 65    pop ASP~UPD~PAL    groupIX 14164215   C   T      65 -0.488721805
#> 66    pop ASP~UPD~PAL    groupVI 14236496   A   C      66  0.153846154
#> 67    pop ASP~UPD~PAL     groupX 14273796   G   A      67 -1.000000000
#> 68    pop ASP~UPD~PAL  groupVIII 14445389   A   T      68 -0.250000000
#> 69    pop ASP~UPD~PAL    groupIX 14868533   A   G      69          NaN
#> 70    pop ASP~UPD~PAL    groupIX 15190693   A   C      70  0.133333333
#> 71    pop ASP~UPD~PAL   groupXII 15306362   A   T      71 -0.200000000
#> 72    pop ASP~UPD~PAL    groupXX 15460430   A   G      72 -0.187500000
#> 73    pop ASP~UPD~PAL   groupVII 15493025   T   G      73  0.125000000
#> 74    pop ASP~UPD~PAL    groupXV 15547073   G   T      74  0.150000000
#> 75    pop ASP~UPD~PAL  groupXIII 15697285   A   C      75 -0.302325581
#> 76    pop ASP~UPD~PAL   groupIII 15983994   C   T      76          NaN
#> 77    pop ASP~UPD~PAL   groupXIX 16233699   A   T      77  0.294117647
#> 78    pop ASP~UPD~PAL    groupIX 16306966   G   T      78          NaN
#> 79    pop ASP~UPD~PAL     groupI 16498302   T   C      79  0.217391304
#> 80    pop ASP~UPD~PAL  groupXIII 16520929   T   C      80 -0.112149533
#> 81    pop ASP~UPD~PAL   groupXVI 16561891   A   T      81 -0.128205128
#> 82    pop ASP~UPD~PAL    groupII 16598364   T   C      82 -0.363636364
#> 83    pop ASP~UPD~PAL    groupIV 17248473   T   A      83  0.714285714
#> 84    pop ASP~UPD~PAL    groupXX 17748944   C   T      84          NaN
#> 85    pop ASP~UPD~PAL  groupXIII 17842228   G   A      85  0.416666667
#> 86    pop ASP~UPD~PAL   groupXIX 17858769   A   G      86  0.239130435
#> 87    pop ASP~UPD~PAL    groupIX 17929225   T   C      87  0.142857143
#> 88    pop ASP~UPD~PAL    groupIX 18339767   A   G      88          NaN
#> 89    pop ASP~UPD~PAL  groupXIII 18836558   C   G      89 -0.473684211
#> 90    pop ASP~UPD~PAL    groupII 19244005   C   T      90  0.309090909
#> 91    pop ASP~UPD~PAL    groupXX 19317273   C   T      91 -0.470588235
#> 92    pop ASP~UPD~PAL     groupI 19759331   G   A      92  0.600000000
#> 93    pop ASP~UPD~PAL   groupVII 20165240   T   C      93 -0.357142857
#> 94    pop ASP~UPD~PAL    groupII 20496069   C   G      94 -1.000000000
#> 95    pop ASP~UPD~PAL     groupI 21284266   G   T      95  0.111111111
#> 96    pop ASP~UPD~PAL    groupIV 22227311   C   A      96 -0.111111111
#> 97    pop ASP~UPD~PAL     groupI 22746027   G   A      97 -1.000000000
#> 98    pop ASP~UPD~PAL   groupVII 25239495   G   T      98 -0.298969072
#> 99    pop ASP~UPD~PAL   groupVII 27599144   G   A      99  0.100000000
#> 100   pop ASP~UPD~PAL    groupIV 32321895   A   T     100 -0.014925373
#>           abba       baba nk
#> 1   0.19852941 0.19852941 78
#> 2   0.05000000 0.00000000 84
#> 3   0.08219697 0.05871212 84
#> 4   0.07867133 0.36538462 80
#> 5   0.15707672 0.22619048 56
#> 6   0.09542411 0.15904018 80
#> 7   0.04545455 0.00000000 84
#> 8   0.17045455 0.19318182 70
#> 9   0.16089379 0.20154065 80
#> 10  0.07792208 0.12987013 84
#> 11  0.27321429 0.21607143 80
#> 12  0.12031250 0.15468750 72
#> 13  0.09088319 0.13219373 74
#> 14  0.12666667 0.29333333 84
#> 15  0.01529720 0.04458042 80
#> 16  0.11215278 0.08854167 74
#> 17  0.10312500 0.10625000 78
#> 18  0.04857143 0.08666667 78
#> 19  0.00000000 0.06666667 86
#> 20  0.05158730 0.06746032 76
#> 21  0.11111111 0.27777778 74
#> 22  0.28571429 0.19841270 80
#> 23  0.04833333 0.03166667 82
#> 24  0.04464286 0.10178571 80
#> 25  0.17000000 0.30333333 82
#> 26  0.11018519 0.09074074 72
#> 27  0.04642857 0.06785714 82
#> 28  0.08035714 0.22500000 74
#> 29  0.04095442 0.09081197 74
#> 30  0.07575758 0.15151515 80
#> 31  0.15454545 0.31515152 86
#> 32  0.06576923 0.13730769 76
#> 33  0.13671329 0.19300699 78
#> 34  0.33585859 0.08333333 74
#> 35  0.09285714 0.06428571 80
#> 36  0.13846154 0.06538462 78
#> 37  0.14158654 0.20865385 78
#> 38  0.20039062 0.08164063 80
#> 39  0.12234848 0.10795455 60
#> 40  0.22159091 0.18750000 72
#> 41  0.00000000 0.00000000 68
#> 42  0.15066964 0.25111607 76
#> 43  0.06993007 0.20979021 82
#> 44  0.07916667 0.10416667 72
#> 45  0.00000000 0.06666667 86
#> 46  0.28605769 0.19471154 74
#> 47  0.17045455 0.11590909 82
#> 48  0.10378510 0.16865079 72
#> 49  0.00000000 0.07692308 80
#> 50  0.13493304 0.05881696 80
#> 51  0.19318182 0.19642857 68
#> 52  0.12662338 0.12310606 74
#> 53  0.11666667 0.33333333 76
#> 54  0.00000000 0.00000000 84
#> 55  0.14358974 0.12307692 72
#> 56  0.16770833 0.10208333 80
#> 57  0.06250000 0.00000000 60
#> 58  0.18648019 0.08741259 66
#> 59  0.26223776 0.19667832 80
#> 60  0.29081633 0.18367347 76
#> 61  0.00000000 0.08333333 66
#> 62  0.08928571 0.09642857 80
#> 63  0.05480769 0.03365385 62
#> 64  0.08636364 0.10714286 80
#> 65  0.08116883 0.23634454 84
#> 66  0.09915865 0.07271635 76
#> 67  0.00000000 0.02820513 74
#> 68  0.04326923 0.07211538 74
#> 69  0.00000000 0.00000000 84
#> 70  0.24791667 0.18958333 76
#> 71  0.03446970 0.05170455 84
#> 72  0.06416667 0.09378205 76
#> 73  0.15750000 0.12250000 78
#> 74  0.13269231 0.09807692 80
#> 75  0.11458333 0.21388889 80
#> 76  0.00000000 0.00000000 56
#> 77  0.16587302 0.09047619 76
#> 78  0.00000000 0.00000000 84
#> 79  0.12165775 0.07820856 80
#> 80  0.08304196 0.10402098 78
#> 81  0.11383929 0.14732143 78
#> 82  0.08076923 0.17307692 80
#> 83  0.17532468 0.02922078 82
#> 84  0.00000000 0.00000000 76
#> 85  0.03946429 0.01625000 78
#> 86  0.16193182 0.09943182 64
#> 87  0.28571429 0.21428571 82
#> 88  0.00000000 0.00000000 68
#> 89  0.06944444 0.19444444 66
#> 90  0.12587413 0.06643357 82
#> 91  0.09709821 0.26971726 78
#> 92  0.14545455 0.03636364 56
#> 93  0.08766234 0.18506494 78
#> 94  0.00000000 0.07005495 70
#> 95  0.16666667 0.13333333 78
#> 96  0.05769231 0.07211538 76
#> 97  0.00000000 0.03333333 86
#> 98  0.08627137 0.15985577 76
#> 99  0.27500000 0.22500000 70
#> 100 0.16071429 0.16558442 70
#> 
#> $weighted.means
#>   facet    subfacet snp.facet snp.subfacet D_abba_baba_jackknife
#> 1   pop ASP~UPD~PAL     .base        .base           -0.06393297
#>   D_abba_baba_p_value D_abba_baba
#> 1          0.05981081 -0.06403927
#> 

# smoothed windowed averages
x <- calc_smoothed_averages(x, "pop.chr", sigma = 200, step = 200, 
   nk = TRUE, stats.type = "pairwise")
get.snpR.stats(x, "pop.chr", "abba_baba")
#> $pairwise
#>     facet  comparison        chr position ref anc .snp.id  D_abba_baba
#> 1     pop ASP~UPD~PAL    groupVI   212436   G   A       1  0.000000000
#> 2     pop ASP~UPD~PAL groupXVIII   350070   G   T       2  1.000000000
#> 3     pop ASP~UPD~PAL    groupXV   854285   C   T       3  0.166666667
#> 4     pop ASP~UPD~PAL  groupXIII   907279   A   G       4 -0.645669291
#> 5     pop ASP~UPD~PAL    groupVI  1068620   G   A       5 -0.180327869
#> 6     pop ASP~UPD~PAL     groupV  1392753   A   C       6 -0.250000000
#> 7     pop ASP~UPD~PAL    groupIX  2075079   T   G       7  1.000000000
#> 8     pop ASP~UPD~PAL  groupXIII  2106000   G   A       8 -0.062500000
#> 9     pop ASP~UPD~PAL   groupXII  2136587   C   A       9 -0.112149533
#> 10    pop ASP~UPD~PAL  groupXVII  2560839   A   G      10 -0.250000000
#> 11    pop ASP~UPD~PAL   groupIII  2750695   T   C      11  0.116788321
#> 12    pop ASP~UPD~PAL   groupXIX  3008524   A   G      12 -0.125000000
#> 13    pop ASP~UPD~PAL   groupXVI  3624417   T   C      13 -0.185185185
#> 14    pop ASP~UPD~PAL   groupVII  3667528   G   T      14 -0.396825397
#> 15    pop ASP~UPD~PAL   groupVII  4189120   A   G      15 -0.489051095
#> 16    pop ASP~UPD~PAL    groupXX  4538074   G   T      16  0.117647059
#> 17    pop ASP~UPD~PAL  groupXIII  4706089   A   C      17 -0.014925373
#> 18    pop ASP~UPD~PAL   groupIII  4746808   C   A      18 -0.281690141
#> 19    pop ASP~UPD~PAL  groupXVII  5041995   T   G      19 -1.000000000
#> 20    pop ASP~UPD~PAL     groupX  5161637   G   A      20 -0.133333333
#> 21    pop ASP~UPD~PAL   groupXII  5178037   C   A      21 -0.428571429
#> 22    pop ASP~UPD~PAL    groupIX  5287711   C   T      22  0.180327869
#> 23    pop ASP~UPD~PAL  groupVIII  5462135   G   A      23  0.208333333
#> 24    pop ASP~UPD~PAL    groupIX  5960933   A   G      24 -0.390243902
#> 25    pop ASP~UPD~PAL  groupXIII  6076868   G   A      25 -0.281690141
#> 26    pop ASP~UPD~PAL    groupXI  6177332   T   C      26  0.096774194
#> 27    pop ASP~UPD~PAL    groupIV  6291287   A   C      27 -0.187500000
#> 28    pop ASP~UPD~PAL    groupIX  6515796   T   C      28 -0.473684211
#> 29    pop ASP~UPD~PAL   groupXII  6718079   A   C      29 -0.378378378
#> 30    pop ASP~UPD~PAL   groupVII  6967477   G   C      30 -0.333333333
#> 31    pop ASP~UPD~PAL   groupIII  7194069   C   G      31 -0.341935484
#> 32    pop ASP~UPD~PAL   groupXVI  7199722   A   C      32 -0.352272727
#> 33    pop ASP~UPD~PAL  groupVIII  7224065   C   T      33 -0.170731707
#> 34    pop ASP~UPD~PAL groupXVIII  7229426   C   G      34  0.602409639
#> 35    pop ASP~UPD~PAL groupXVIII  8087022   A   C      35  0.181818182
#> 36    pop ASP~UPD~PAL    groupIX  8373442   A   G      36  0.358490566
#> 37    pop ASP~UPD~PAL    groupXV  8815538   G   C      37 -0.191489362
#> 38    pop ASP~UPD~PAL    groupXI  9019877   G   T      38  0.421052632
#> 39    pop ASP~UPD~PAL   groupXXI  9067876   T   A      39  0.062500000
#> 40    pop ASP~UPD~PAL  groupXVII  9108184   G   C      40  0.083333333
#> 41    pop ASP~UPD~PAL    groupIX  9284427   G   A      41          NaN
#> 42    pop ASP~UPD~PAL     groupV  9767769   G   A      42 -0.250000000
#> 43    pop ASP~UPD~PAL   groupXII  9881312   A   G      43 -0.500000000
#> 44    pop ASP~UPD~PAL    groupXX 10296516   T   A      44 -0.136363636
#> 45    pop ASP~UPD~PAL     groupI 10593829   G   T      45 -1.000000000
#> 46    pop ASP~UPD~PAL   groupXIX 10765968   C   A      46  0.190000000
#> 47    pop ASP~UPD~PAL    groupII 10829169   C   G      47  0.190476190
#> 48    pop ASP~UPD~PAL    groupXX 10913115   T   C      48 -0.238095238
#> 49    pop ASP~UPD~PAL groupXVIII 10974594   G   C      49 -1.000000000
#> 50    pop ASP~UPD~PAL  groupXVII 11256661   C   T      50  0.392857143
#> 51    pop ASP~UPD~PAL  groupVIII 11438870   T   C      51 -0.008333333
#> 52    pop ASP~UPD~PAL   groupVII 11477223   A   G      52  0.014084507
#> 53    pop ASP~UPD~PAL    groupII 11505084   A   C      53 -0.481481481
#> 54    pop ASP~UPD~PAL    groupIX 11904651   G   T      54          NaN
#> 55    pop ASP~UPD~PAL  groupVIII 12108945   G   A      55  0.076923077
#> 56    pop ASP~UPD~PAL   groupXIX 12485182   A   G      56  0.243243243
#> 57    pop ASP~UPD~PAL    groupXI 12833733   G   A      57  1.000000000
#> 58    pop ASP~UPD~PAL groupXVIII 12937786   C   T      58  0.361702128
#> 59    pop ASP~UPD~PAL   groupIII 13177027   G   A      59  0.142857143
#> 60    pop ASP~UPD~PAL    groupVI 13335577   A   C      60  0.225806452
#> 61    pop ASP~UPD~PAL     groupX 13718292   T   C      61 -1.000000000
#> 62    pop ASP~UPD~PAL    groupVI 13895716   G   C      62 -0.038461538
#> 63    pop ASP~UPD~PAL   groupXIV 14009990   A   C      63  0.239130435
#> 64    pop ASP~UPD~PAL   groupXII 14158350   T   A      64 -0.107382550
#> 65    pop ASP~UPD~PAL    groupIX 14164215   C   T      65 -0.488721805
#> 66    pop ASP~UPD~PAL    groupVI 14236496   A   C      66  0.153846154
#> 67    pop ASP~UPD~PAL     groupX 14273796   G   A      67 -1.000000000
#> 68    pop ASP~UPD~PAL  groupVIII 14445389   A   T      68 -0.250000000
#> 69    pop ASP~UPD~PAL    groupIX 14868533   A   G      69          NaN
#> 70    pop ASP~UPD~PAL    groupIX 15190693   A   C      70  0.133333333
#> 71    pop ASP~UPD~PAL   groupXII 15306362   A   T      71 -0.200000000
#> 72    pop ASP~UPD~PAL    groupXX 15460430   A   G      72 -0.187500000
#> 73    pop ASP~UPD~PAL   groupVII 15493025   T   G      73  0.125000000
#> 74    pop ASP~UPD~PAL    groupXV 15547073   G   T      74  0.150000000
#> 75    pop ASP~UPD~PAL  groupXIII 15697285   A   C      75 -0.302325581
#> 76    pop ASP~UPD~PAL   groupIII 15983994   C   T      76          NaN
#> 77    pop ASP~UPD~PAL   groupXIX 16233699   A   T      77  0.294117647
#> 78    pop ASP~UPD~PAL    groupIX 16306966   G   T      78          NaN
#> 79    pop ASP~UPD~PAL     groupI 16498302   T   C      79  0.217391304
#> 80    pop ASP~UPD~PAL  groupXIII 16520929   T   C      80 -0.112149533
#> 81    pop ASP~UPD~PAL   groupXVI 16561891   A   T      81 -0.128205128
#> 82    pop ASP~UPD~PAL    groupII 16598364   T   C      82 -0.363636364
#> 83    pop ASP~UPD~PAL    groupIV 17248473   T   A      83  0.714285714
#> 84    pop ASP~UPD~PAL    groupXX 17748944   C   T      84          NaN
#> 85    pop ASP~UPD~PAL  groupXIII 17842228   G   A      85  0.416666667
#> 86    pop ASP~UPD~PAL   groupXIX 17858769   A   G      86  0.239130435
#> 87    pop ASP~UPD~PAL    groupIX 17929225   T   C      87  0.142857143
#> 88    pop ASP~UPD~PAL    groupIX 18339767   A   G      88          NaN
#> 89    pop ASP~UPD~PAL  groupXIII 18836558   C   G      89 -0.473684211
#> 90    pop ASP~UPD~PAL    groupII 19244005   C   T      90  0.309090909
#> 91    pop ASP~UPD~PAL    groupXX 19317273   C   T      91 -0.470588235
#> 92    pop ASP~UPD~PAL     groupI 19759331   G   A      92  0.600000000
#> 93    pop ASP~UPD~PAL   groupVII 20165240   T   C      93 -0.357142857
#> 94    pop ASP~UPD~PAL    groupII 20496069   C   G      94 -1.000000000
#> 95    pop ASP~UPD~PAL     groupI 21284266   G   T      95  0.111111111
#> 96    pop ASP~UPD~PAL    groupIV 22227311   C   A      96 -0.111111111
#> 97    pop ASP~UPD~PAL     groupI 22746027   G   A      97 -1.000000000
#> 98    pop ASP~UPD~PAL   groupVII 25239495   G   T      98 -0.298969072
#> 99    pop ASP~UPD~PAL   groupVII 27599144   G   A      99  0.100000000
#> 100   pop ASP~UPD~PAL    groupIV 32321895   A   T     100 -0.014925373
#>           abba       baba nk
#> 1   0.19852941 0.19852941 78
#> 2   0.05000000 0.00000000 84
#> 3   0.08219697 0.05871212 84
#> 4   0.07867133 0.36538462 80
#> 5   0.15707672 0.22619048 56
#> 6   0.09542411 0.15904018 80
#> 7   0.04545455 0.00000000 84
#> 8   0.17045455 0.19318182 70
#> 9   0.16089379 0.20154065 80
#> 10  0.07792208 0.12987013 84
#> 11  0.27321429 0.21607143 80
#> 12  0.12031250 0.15468750 72
#> 13  0.09088319 0.13219373 74
#> 14  0.12666667 0.29333333 84
#> 15  0.01529720 0.04458042 80
#> 16  0.11215278 0.08854167 74
#> 17  0.10312500 0.10625000 78
#> 18  0.04857143 0.08666667 78
#> 19  0.00000000 0.06666667 86
#> 20  0.05158730 0.06746032 76
#> 21  0.11111111 0.27777778 74
#> 22  0.28571429 0.19841270 80
#> 23  0.04833333 0.03166667 82
#> 24  0.04464286 0.10178571 80
#> 25  0.17000000 0.30333333 82
#> 26  0.11018519 0.09074074 72
#> 27  0.04642857 0.06785714 82
#> 28  0.08035714 0.22500000 74
#> 29  0.04095442 0.09081197 74
#> 30  0.07575758 0.15151515 80
#> 31  0.15454545 0.31515152 86
#> 32  0.06576923 0.13730769 76
#> 33  0.13671329 0.19300699 78
#> 34  0.33585859 0.08333333 74
#> 35  0.09285714 0.06428571 80
#> 36  0.13846154 0.06538462 78
#> 37  0.14158654 0.20865385 78
#> 38  0.20039062 0.08164063 80
#> 39  0.12234848 0.10795455 60
#> 40  0.22159091 0.18750000 72
#> 41  0.00000000 0.00000000 68
#> 42  0.15066964 0.25111607 76
#> 43  0.06993007 0.20979021 82
#> 44  0.07916667 0.10416667 72
#> 45  0.00000000 0.06666667 86
#> 46  0.28605769 0.19471154 74
#> 47  0.17045455 0.11590909 82
#> 48  0.10378510 0.16865079 72
#> 49  0.00000000 0.07692308 80
#> 50  0.13493304 0.05881696 80
#> 51  0.19318182 0.19642857 68
#> 52  0.12662338 0.12310606 74
#> 53  0.11666667 0.33333333 76
#> 54  0.00000000 0.00000000 84
#> 55  0.14358974 0.12307692 72
#> 56  0.16770833 0.10208333 80
#> 57  0.06250000 0.00000000 60
#> 58  0.18648019 0.08741259 66
#> 59  0.26223776 0.19667832 80
#> 60  0.29081633 0.18367347 76
#> 61  0.00000000 0.08333333 66
#> 62  0.08928571 0.09642857 80
#> 63  0.05480769 0.03365385 62
#> 64  0.08636364 0.10714286 80
#> 65  0.08116883 0.23634454 84
#> 66  0.09915865 0.07271635 76
#> 67  0.00000000 0.02820513 74
#> 68  0.04326923 0.07211538 74
#> 69  0.00000000 0.00000000 84
#> 70  0.24791667 0.18958333 76
#> 71  0.03446970 0.05170455 84
#> 72  0.06416667 0.09378205 76
#> 73  0.15750000 0.12250000 78
#> 74  0.13269231 0.09807692 80
#> 75  0.11458333 0.21388889 80
#> 76  0.00000000 0.00000000 56
#> 77  0.16587302 0.09047619 76
#> 78  0.00000000 0.00000000 84
#> 79  0.12165775 0.07820856 80
#> 80  0.08304196 0.10402098 78
#> 81  0.11383929 0.14732143 78
#> 82  0.08076923 0.17307692 80
#> 83  0.17532468 0.02922078 82
#> 84  0.00000000 0.00000000 76
#> 85  0.03946429 0.01625000 78
#> 86  0.16193182 0.09943182 64
#> 87  0.28571429 0.21428571 82
#> 88  0.00000000 0.00000000 68
#> 89  0.06944444 0.19444444 66
#> 90  0.12587413 0.06643357 82
#> 91  0.09709821 0.26971726 78
#> 92  0.14545455 0.03636364 56
#> 93  0.08766234 0.18506494 78
#> 94  0.00000000 0.07005495 70
#> 95  0.16666667 0.13333333 78
#> 96  0.05769231 0.07211538 76
#> 97  0.00000000 0.03333333 86
#> 98  0.08627137 0.15985577 76
#> 99  0.27500000 0.22500000 70
#> 100 0.16071429 0.16558442 70
#> 
#> $pairwise.window
#>     facet    subfacet snp.facet snp.subfacet position sigma step nk.status
#> 1     pop ASP~UPD~PAL       chr      groupVI   200000   200  200      TRUE
#> 2     pop ASP~UPD~PAL       chr   groupXVIII   200000   200  200      TRUE
#> 3     pop ASP~UPD~PAL       chr      groupVI   400000   200  200      TRUE
#> 4     pop ASP~UPD~PAL       chr   groupXVIII   400000   200  200      TRUE
#> 5     pop ASP~UPD~PAL       chr      groupXV   800000   200  200      TRUE
#> 6     pop ASP~UPD~PAL       chr    groupXIII   800000   200  200      TRUE
#> 7     pop ASP~UPD~PAL       chr      groupVI  1000000   200  200      TRUE
#> 8     pop ASP~UPD~PAL       chr      groupXV  1000000   200  200      TRUE
#> 9     pop ASP~UPD~PAL       chr    groupXIII  1000000   200  200      TRUE
#> 10    pop ASP~UPD~PAL       chr      groupVI  1200000   200  200      TRUE
#> 11    pop ASP~UPD~PAL       chr       groupV  1200000   200  200      TRUE
#> 12    pop ASP~UPD~PAL       chr       groupV  1400000   200  200      TRUE
#> 13    pop ASP~UPD~PAL       chr    groupXIII  2000000   200  200      TRUE
#> 14    pop ASP~UPD~PAL       chr      groupIX  2000000   200  200      TRUE
#> 15    pop ASP~UPD~PAL       chr     groupXII  2000000   200  200      TRUE
#> 16    pop ASP~UPD~PAL       chr    groupXIII  2200000   200  200      TRUE
#> 17    pop ASP~UPD~PAL       chr      groupIX  2200000   200  200      TRUE
#> 18    pop ASP~UPD~PAL       chr     groupXII  2200000   200  200      TRUE
#> 19    pop ASP~UPD~PAL       chr    groupXVII  2400000   200  200      TRUE
#> 20    pop ASP~UPD~PAL       chr    groupXVII  2600000   200  200      TRUE
#> 21    pop ASP~UPD~PAL       chr     groupIII  2600000   200  200      TRUE
#> 22    pop ASP~UPD~PAL       chr     groupIII  2800000   200  200      TRUE
#> 23    pop ASP~UPD~PAL       chr     groupXIX  3000000   200  200      TRUE
#> 24    pop ASP~UPD~PAL       chr     groupXIX  3200000   200  200      TRUE
#> 25    pop ASP~UPD~PAL       chr     groupXVI  3600000   200  200      TRUE
#> 26    pop ASP~UPD~PAL       chr     groupVII  3600000   200  200      TRUE
#> 27    pop ASP~UPD~PAL       chr     groupXVI  3800000   200  200      TRUE
#> 28    pop ASP~UPD~PAL       chr     groupVII  3800000   200  200      TRUE
#> 29    pop ASP~UPD~PAL       chr     groupVII  4000000   200  200      TRUE
#> 30    pop ASP~UPD~PAL       chr     groupVII  4200000   200  200      TRUE
#> 31    pop ASP~UPD~PAL       chr      groupXX  4400000   200  200      TRUE
#> 32    pop ASP~UPD~PAL       chr    groupXIII  4600000   200  200      TRUE
#> 33    pop ASP~UPD~PAL       chr     groupIII  4600000   200  200      TRUE
#> 34    pop ASP~UPD~PAL       chr      groupXX  4600000   200  200      TRUE
#> 35    pop ASP~UPD~PAL       chr    groupXIII  4800000   200  200      TRUE
#> 36    pop ASP~UPD~PAL       chr     groupIII  4800000   200  200      TRUE
#> 37    pop ASP~UPD~PAL       chr     groupXII  5000000   200  200      TRUE
#> 38    pop ASP~UPD~PAL       chr    groupXVII  5000000   200  200      TRUE
#> 39    pop ASP~UPD~PAL       chr       groupX  5000000   200  200      TRUE
#> 40    pop ASP~UPD~PAL       chr      groupIX  5200000   200  200      TRUE
#> 41    pop ASP~UPD~PAL       chr     groupXII  5200000   200  200      TRUE
#> 42    pop ASP~UPD~PAL       chr    groupXVII  5200000   200  200      TRUE
#> 43    pop ASP~UPD~PAL       chr       groupX  5200000   200  200      TRUE
#> 44    pop ASP~UPD~PAL       chr      groupIX  5400000   200  200      TRUE
#> 45    pop ASP~UPD~PAL       chr    groupVIII  5400000   200  200      TRUE
#> 46    pop ASP~UPD~PAL       chr    groupVIII  5600000   200  200      TRUE
#> 47    pop ASP~UPD~PAL       chr      groupIX  5800000   200  200      TRUE
#> 48    pop ASP~UPD~PAL       chr    groupXIII  6000000   200  200      TRUE
#> 49    pop ASP~UPD~PAL       chr      groupIX  6000000   200  200      TRUE
#> 50    pop ASP~UPD~PAL       chr      groupXI  6000000   200  200      TRUE
#> 51    pop ASP~UPD~PAL       chr    groupXIII  6200000   200  200      TRUE
#> 52    pop ASP~UPD~PAL       chr      groupXI  6200000   200  200      TRUE
#> 53    pop ASP~UPD~PAL       chr      groupIV  6200000   200  200      TRUE
#> 54    pop ASP~UPD~PAL       chr      groupIX  6400000   200  200      TRUE
#> 55    pop ASP~UPD~PAL       chr      groupIV  6400000   200  200      TRUE
#> 56    pop ASP~UPD~PAL       chr      groupIX  6600000   200  200      TRUE
#> 57    pop ASP~UPD~PAL       chr     groupXII  6600000   200  200      TRUE
#> 58    pop ASP~UPD~PAL       chr     groupXII  6800000   200  200      TRUE
#> 59    pop ASP~UPD~PAL       chr     groupVII  6800000   200  200      TRUE
#> 60    pop ASP~UPD~PAL       chr     groupIII  7000000   200  200      TRUE
#> 61    pop ASP~UPD~PAL       chr     groupXVI  7000000   200  200      TRUE
#> 62    pop ASP~UPD~PAL       chr     groupVII  7000000   200  200      TRUE
#> 63    pop ASP~UPD~PAL       chr   groupXVIII  7200000   200  200      TRUE
#> 64    pop ASP~UPD~PAL       chr     groupIII  7200000   200  200      TRUE
#> 65    pop ASP~UPD~PAL       chr     groupXVI  7200000   200  200      TRUE
#> 66    pop ASP~UPD~PAL       chr    groupVIII  7200000   200  200      TRUE
#> 67    pop ASP~UPD~PAL       chr   groupXVIII  7400000   200  200      TRUE
#> 68    pop ASP~UPD~PAL       chr    groupVIII  7400000   200  200      TRUE
#> 69    pop ASP~UPD~PAL       chr   groupXVIII  8000000   200  200      TRUE
#> 70    pop ASP~UPD~PAL       chr   groupXVIII  8200000   200  200      TRUE
#> 71    pop ASP~UPD~PAL       chr      groupIX  8200000   200  200      TRUE
#> 72    pop ASP~UPD~PAL       chr      groupIX  8400000   200  200      TRUE
#> 73    pop ASP~UPD~PAL       chr      groupXV  8800000   200  200      TRUE
#> 74    pop ASP~UPD~PAL       chr      groupXV  9000000   200  200      TRUE
#> 75    pop ASP~UPD~PAL       chr    groupXVII  9000000   200  200      TRUE
#> 76    pop ASP~UPD~PAL       chr      groupXI  9000000   200  200      TRUE
#> 77    pop ASP~UPD~PAL       chr     groupXXI  9000000   200  200      TRUE
#> 78    pop ASP~UPD~PAL       chr      groupIX  9200000   200  200      TRUE
#> 79    pop ASP~UPD~PAL       chr    groupXVII  9200000   200  200      TRUE
#> 80    pop ASP~UPD~PAL       chr      groupXI  9200000   200  200      TRUE
#> 81    pop ASP~UPD~PAL       chr     groupXXI  9200000   200  200      TRUE
#> 82    pop ASP~UPD~PAL       chr      groupIX  9400000   200  200      TRUE
#> 83    pop ASP~UPD~PAL       chr       groupV  9600000   200  200      TRUE
#> 84    pop ASP~UPD~PAL       chr       groupV  9800000   200  200      TRUE
#> 85    pop ASP~UPD~PAL       chr     groupXII  9800000   200  200      TRUE
#> 86    pop ASP~UPD~PAL       chr     groupXII 10000000   200  200      TRUE
#> 87    pop ASP~UPD~PAL       chr      groupXX 10200000   200  200      TRUE
#> 88    pop ASP~UPD~PAL       chr      groupXX 10400000   200  200      TRUE
#> 89    pop ASP~UPD~PAL       chr       groupI 10400000   200  200      TRUE
#> 90    pop ASP~UPD~PAL       chr     groupXIX 10600000   200  200      TRUE
#> 91    pop ASP~UPD~PAL       chr       groupI 10600000   200  200      TRUE
#> 92    pop ASP~UPD~PAL       chr   groupXVIII 10800000   200  200      TRUE
#> 93    pop ASP~UPD~PAL       chr     groupXIX 10800000   200  200      TRUE
#> 94    pop ASP~UPD~PAL       chr      groupXX 10800000   200  200      TRUE
#> 95    pop ASP~UPD~PAL       chr      groupII 10800000   200  200      TRUE
#> 96    pop ASP~UPD~PAL       chr   groupXVIII 11000000   200  200      TRUE
#> 97    pop ASP~UPD~PAL       chr      groupXX 11000000   200  200      TRUE
#> 98    pop ASP~UPD~PAL       chr      groupII 11000000   200  200      TRUE
#> 99    pop ASP~UPD~PAL       chr    groupXVII 11200000   200  200      TRUE
#> 100   pop ASP~UPD~PAL       chr    groupXVII 11400000   200  200      TRUE
#> 101   pop ASP~UPD~PAL       chr     groupVII 11400000   200  200      TRUE
#> 102   pop ASP~UPD~PAL       chr    groupVIII 11400000   200  200      TRUE
#> 103   pop ASP~UPD~PAL       chr      groupII 11400000   200  200      TRUE
#> 104   pop ASP~UPD~PAL       chr     groupVII 11600000   200  200      TRUE
#> 105   pop ASP~UPD~PAL       chr    groupVIII 11600000   200  200      TRUE
#> 106   pop ASP~UPD~PAL       chr      groupII 11600000   200  200      TRUE
#> 107   pop ASP~UPD~PAL       chr      groupIX 11800000   200  200      TRUE
#> 108   pop ASP~UPD~PAL       chr      groupIX 12000000   200  200      TRUE
#> 109   pop ASP~UPD~PAL       chr    groupVIII 12000000   200  200      TRUE
#> 110   pop ASP~UPD~PAL       chr    groupVIII 12200000   200  200      TRUE
#> 111   pop ASP~UPD~PAL       chr     groupXIX 12400000   200  200      TRUE
#> 112   pop ASP~UPD~PAL       chr     groupXIX 12600000   200  200      TRUE
#> 113   pop ASP~UPD~PAL       chr   groupXVIII 12800000   200  200      TRUE
#> 114   pop ASP~UPD~PAL       chr      groupXI 12800000   200  200      TRUE
#> 115   pop ASP~UPD~PAL       chr   groupXVIII 13000000   200  200      TRUE
#> 116   pop ASP~UPD~PAL       chr     groupIII 13000000   200  200      TRUE
#> 117   pop ASP~UPD~PAL       chr      groupXI 13000000   200  200      TRUE
#> 118   pop ASP~UPD~PAL       chr      groupVI 13200000   200  200      TRUE
#> 119   pop ASP~UPD~PAL       chr     groupIII 13200000   200  200      TRUE
#> 120   pop ASP~UPD~PAL       chr      groupVI 13400000   200  200      TRUE
#> 121   pop ASP~UPD~PAL       chr       groupX 13600000   200  200      TRUE
#> 122   pop ASP~UPD~PAL       chr      groupVI 13800000   200  200      TRUE
#> 123   pop ASP~UPD~PAL       chr       groupX 13800000   200  200      TRUE
#> 124   pop ASP~UPD~PAL       chr      groupVI 14000000   200  200      TRUE
#> 125   pop ASP~UPD~PAL       chr      groupIX 14000000   200  200      TRUE
#> 126   pop ASP~UPD~PAL       chr     groupXII 14000000   200  200      TRUE
#> 127   pop ASP~UPD~PAL       chr     groupXIV 14000000   200  200      TRUE
#> 128   pop ASP~UPD~PAL       chr      groupVI 14200000   200  200      TRUE
#> 129   pop ASP~UPD~PAL       chr      groupIX 14200000   200  200      TRUE
#> 130   pop ASP~UPD~PAL       chr     groupXII 14200000   200  200      TRUE
#> 131   pop ASP~UPD~PAL       chr       groupX 14200000   200  200      TRUE
#> 132   pop ASP~UPD~PAL       chr     groupXIV 14200000   200  200      TRUE
#> 133   pop ASP~UPD~PAL       chr      groupVI 14400000   200  200      TRUE
#> 134   pop ASP~UPD~PAL       chr       groupX 14400000   200  200      TRUE
#> 135   pop ASP~UPD~PAL       chr    groupVIII 14400000   200  200      TRUE
#> 136   pop ASP~UPD~PAL       chr    groupVIII 14600000   200  200      TRUE
#> 137   pop ASP~UPD~PAL       chr      groupIX 14800000   200  200      TRUE
#> 138   pop ASP~UPD~PAL       chr      groupIX 15000000   200  200      TRUE
#> 139   pop ASP~UPD~PAL       chr      groupIX 15200000   200  200      TRUE
#> 140   pop ASP~UPD~PAL       chr     groupXII 15200000   200  200      TRUE
#> 141   pop ASP~UPD~PAL       chr      groupXV 15400000   200  200      TRUE
#> 142   pop ASP~UPD~PAL       chr     groupXII 15400000   200  200      TRUE
#> 143   pop ASP~UPD~PAL       chr     groupVII 15400000   200  200      TRUE
#> 144   pop ASP~UPD~PAL       chr      groupXX 15400000   200  200      TRUE
#> 145   pop ASP~UPD~PAL       chr      groupXV 15600000   200  200      TRUE
#> 146   pop ASP~UPD~PAL       chr    groupXIII 15600000   200  200      TRUE
#> 147   pop ASP~UPD~PAL       chr     groupVII 15600000   200  200      TRUE
#> 148   pop ASP~UPD~PAL       chr      groupXX 15600000   200  200      TRUE
#> 149   pop ASP~UPD~PAL       chr    groupXIII 15800000   200  200      TRUE
#> 150   pop ASP~UPD~PAL       chr     groupIII 15800000   200  200      TRUE
#> 151   pop ASP~UPD~PAL       chr     groupIII 16000000   200  200      TRUE
#> 152   pop ASP~UPD~PAL       chr      groupIX 16200000   200  200      TRUE
#> 153   pop ASP~UPD~PAL       chr     groupXIX 16200000   200  200      TRUE
#> 154   pop ASP~UPD~PAL       chr    groupXIII 16400000   200  200      TRUE
#> 155   pop ASP~UPD~PAL       chr      groupIX 16400000   200  200      TRUE
#> 156   pop ASP~UPD~PAL       chr     groupXIX 16400000   200  200      TRUE
#> 157   pop ASP~UPD~PAL       chr     groupXVI 16400000   200  200      TRUE
#> 158   pop ASP~UPD~PAL       chr       groupI 16400000   200  200      TRUE
#> 159   pop ASP~UPD~PAL       chr      groupII 16400000   200  200      TRUE
#> 160   pop ASP~UPD~PAL       chr    groupXIII 16600000   200  200      TRUE
#> 161   pop ASP~UPD~PAL       chr     groupXVI 16600000   200  200      TRUE
#> 162   pop ASP~UPD~PAL       chr       groupI 16600000   200  200      TRUE
#> 163   pop ASP~UPD~PAL       chr      groupII 16600000   200  200      TRUE
#> 164   pop ASP~UPD~PAL       chr      groupIV 17200000   200  200      TRUE
#> 165   pop ASP~UPD~PAL       chr      groupIV 17400000   200  200      TRUE
#> 166   pop ASP~UPD~PAL       chr      groupXX 17600000   200  200      TRUE
#> 167   pop ASP~UPD~PAL       chr    groupXIII 17800000   200  200      TRUE
#> 168   pop ASP~UPD~PAL       chr      groupIX 17800000   200  200      TRUE
#> 169   pop ASP~UPD~PAL       chr     groupXIX 17800000   200  200      TRUE
#> 170   pop ASP~UPD~PAL       chr      groupXX 17800000   200  200      TRUE
#> 171   pop ASP~UPD~PAL       chr    groupXIII 18000000   200  200      TRUE
#> 172   pop ASP~UPD~PAL       chr      groupIX 18000000   200  200      TRUE
#> 173   pop ASP~UPD~PAL       chr     groupXIX 18000000   200  200      TRUE
#> 174   pop ASP~UPD~PAL       chr      groupIX 18200000   200  200      TRUE
#> 175   pop ASP~UPD~PAL       chr      groupIX 18400000   200  200      TRUE
#> 176   pop ASP~UPD~PAL       chr    groupXIII 18800000   200  200      TRUE
#> 177   pop ASP~UPD~PAL       chr    groupXIII 19000000   200  200      TRUE
#> 178   pop ASP~UPD~PAL       chr      groupXX 19200000   200  200      TRUE
#> 179   pop ASP~UPD~PAL       chr      groupII 19200000   200  200      TRUE
#> 180   pop ASP~UPD~PAL       chr      groupXX 19400000   200  200      TRUE
#> 181   pop ASP~UPD~PAL       chr      groupII 19400000   200  200      TRUE
#> 182   pop ASP~UPD~PAL       chr       groupI 19600000   200  200      TRUE
#> 183   pop ASP~UPD~PAL       chr       groupI 19800000   200  200      TRUE
#> 184   pop ASP~UPD~PAL       chr     groupVII 20000000   200  200      TRUE
#> 185   pop ASP~UPD~PAL       chr     groupVII 20200000   200  200      TRUE
#> 186   pop ASP~UPD~PAL       chr      groupII 20400000   200  200      TRUE
#> 187   pop ASP~UPD~PAL       chr      groupII 20600000   200  200      TRUE
#> 188   pop ASP~UPD~PAL       chr       groupI 21200000   200  200      TRUE
#> 189   pop ASP~UPD~PAL       chr       groupI 21400000   200  200      TRUE
#> 190   pop ASP~UPD~PAL       chr      groupIV 22200000   200  200      TRUE
#> 191   pop ASP~UPD~PAL       chr      groupIV 22400000   200  200      TRUE
#> 192   pop ASP~UPD~PAL       chr       groupI 22600000   200  200      TRUE
#> 193   pop ASP~UPD~PAL       chr       groupI 22800000   200  200      TRUE
#> 194   pop ASP~UPD~PAL       chr     groupVII 25200000   200  200      TRUE
#> 195   pop ASP~UPD~PAL       chr     groupVII 25400000   200  200      TRUE
#> 196   pop ASP~UPD~PAL       chr     groupVII 27400000   200  200      TRUE
#> 197   pop ASP~UPD~PAL       chr     groupVII 27600000   200  200      TRUE
#> 198   pop ASP~UPD~PAL       chr      groupIV 32200000   200  200      TRUE
#> 199   pop ASP~UPD~PAL       chr      groupIV 32400000   200  200      TRUE
#>     gaussian n_snps triple_sigma  D_abba_baba       abba       baba
#> 1       TRUE      1        FALSE  0.000000000 0.19852941 0.19852941
#> 2       TRUE      1        FALSE  1.000000000 0.05000000 0.00000000
#> 3       TRUE      1        FALSE  0.000000000 0.19852941 0.19852941
#> 4       TRUE      1        FALSE  1.000000000 0.05000000 0.00000000
#> 5       TRUE      1        FALSE  0.166666667 0.08219697 0.05871212
#> 6       TRUE      1        FALSE -0.645669291 0.07867133 0.36538462
#> 7       TRUE      1        FALSE -0.180327869 0.15707672 0.22619048
#> 8       TRUE      1        FALSE  0.166666667 0.08219697 0.05871212
#> 9       TRUE      1        FALSE -0.645669291 0.07867133 0.36538462
#> 10      TRUE      1        FALSE -0.180327869 0.15707672 0.22619048
#> 11      TRUE      1        FALSE -0.250000000 0.09542411 0.15904018
#> 12      TRUE      1        FALSE -0.250000000 0.09542411 0.15904018
#> 13      TRUE      1        FALSE -0.062500000 0.17045455 0.19318182
#> 14      TRUE      1        FALSE  1.000000000 0.04545455 0.00000000
#> 15      TRUE      1        FALSE -0.112149533 0.16089379 0.20154065
#> 16      TRUE      1        FALSE -0.062500000 0.17045455 0.19318182
#> 17      TRUE      1        FALSE  1.000000000 0.04545455 0.00000000
#> 18      TRUE      1        FALSE -0.112149533 0.16089379 0.20154065
#> 19      TRUE      1        FALSE -0.250000000 0.07792208 0.12987013
#> 20      TRUE      1        FALSE -0.250000000 0.07792208 0.12987013
#> 21      TRUE      1        FALSE  0.116788321 0.27321429 0.21607143
#> 22      TRUE      1        FALSE  0.116788321 0.27321429 0.21607143
#> 23      TRUE      1        FALSE -0.125000000 0.12031250 0.15468750
#> 24      TRUE      1        FALSE -0.125000000 0.12031250 0.15468750
#> 25      TRUE      1        FALSE -0.185185185 0.09088319 0.13219373
#> 26      TRUE      1        FALSE -0.396825397 0.12666667 0.29333333
#> 27      TRUE      1        FALSE -0.185185185 0.09088319 0.13219373
#> 28      TRUE      1        FALSE -0.396825397 0.12666667 0.29333333
#> 29      TRUE      1        FALSE -0.489051095 0.01529720 0.04458042
#> 30      TRUE      1        FALSE -0.489051095 0.01529720 0.04458042
#> 31      TRUE      1        FALSE  0.117647059 0.11215278 0.08854167
#> 32      TRUE      1        FALSE -0.014925373 0.10312500 0.10625000
#> 33      TRUE      1        FALSE -0.281690141 0.04857143 0.08666667
#> 34      TRUE      1        FALSE  0.117647059 0.11215278 0.08854167
#> 35      TRUE      1        FALSE -0.014925373 0.10312500 0.10625000
#> 36      TRUE      1        FALSE -0.281690141 0.04857143 0.08666667
#> 37      TRUE      1        FALSE -0.428571429 0.11111111 0.27777778
#> 38      TRUE      1        FALSE -1.000000000 0.00000000 0.06666667
#> 39      TRUE      1        FALSE -0.133333333 0.05158730 0.06746032
#> 40      TRUE      1        FALSE  0.180327869 0.28571429 0.19841270
#> 41      TRUE      1        FALSE -0.428571429 0.11111111 0.27777778
#> 42      TRUE      1        FALSE -1.000000000 0.00000000 0.06666667
#> 43      TRUE      1        FALSE -0.133333333 0.05158730 0.06746032
#> 44      TRUE      1        FALSE  0.180327869 0.28571429 0.19841270
#> 45      TRUE      1        FALSE  0.208333333 0.04833333 0.03166667
#> 46      TRUE      1        FALSE  0.208333333 0.04833333 0.03166667
#> 47      TRUE      1        FALSE -0.390243902 0.04464286 0.10178571
#> 48      TRUE      1        FALSE -0.281690141 0.17000000 0.30333333
#> 49      TRUE      1        FALSE -0.390243902 0.04464286 0.10178571
#> 50      TRUE      1        FALSE  0.096774194 0.11018519 0.09074074
#> 51      TRUE      1        FALSE -0.281690141 0.17000000 0.30333333
#> 52      TRUE      1        FALSE  0.096774194 0.11018519 0.09074074
#> 53      TRUE      1        FALSE -0.187500000 0.04642857 0.06785714
#> 54      TRUE      1        FALSE -0.473684211 0.08035714 0.22500000
#> 55      TRUE      1        FALSE -0.187500000 0.04642857 0.06785714
#> 56      TRUE      1        FALSE -0.473684211 0.08035714 0.22500000
#> 57      TRUE      1        FALSE -0.378378378 0.04095442 0.09081197
#> 58      TRUE      1        FALSE -0.378378378 0.04095442 0.09081197
#> 59      TRUE      1        FALSE -0.333333333 0.07575758 0.15151515
#> 60      TRUE      1        FALSE -0.341935484 0.15454545 0.31515152
#> 61      TRUE      1        FALSE -0.352272727 0.06576923 0.13730769
#> 62      TRUE      1        FALSE -0.333333333 0.07575758 0.15151515
#> 63      TRUE      1        FALSE  0.602409639 0.33585859 0.08333333
#> 64      TRUE      1        FALSE -0.341935484 0.15454545 0.31515152
#> 65      TRUE      1        FALSE -0.352272727 0.06576923 0.13730769
#> 66      TRUE      1        FALSE -0.170731707 0.13671329 0.19300699
#> 67      TRUE      1        FALSE  0.602409639 0.33585859 0.08333333
#> 68      TRUE      1        FALSE -0.170731707 0.13671329 0.19300699
#> 69      TRUE      1        FALSE  0.181818182 0.09285714 0.06428571
#> 70      TRUE      1        FALSE  0.181818182 0.09285714 0.06428571
#> 71      TRUE      1        FALSE  0.358490566 0.13846154 0.06538462
#> 72      TRUE      1        FALSE  0.358490566 0.13846154 0.06538462
#> 73      TRUE      1        FALSE -0.191489362 0.14158654 0.20865385
#> 74      TRUE      1        FALSE -0.191489362 0.14158654 0.20865385
#> 75      TRUE      1        FALSE  0.083333333 0.22159091 0.18750000
#> 76      TRUE      1        FALSE  0.421052632 0.20039062 0.08164063
#> 77      TRUE      1        FALSE  0.062500000 0.12234848 0.10795455
#> 78      TRUE      1        FALSE          NaN 0.00000000 0.00000000
#> 79      TRUE      1        FALSE  0.083333333 0.22159091 0.18750000
#> 80      TRUE      1        FALSE  0.421052632 0.20039062 0.08164062
#> 81      TRUE      1        FALSE  0.062500000 0.12234848 0.10795455
#> 82      TRUE      1        FALSE          NaN 0.00000000 0.00000000
#> 83      TRUE      1        FALSE -0.250000000 0.15066964 0.25111607
#> 84      TRUE      1        FALSE -0.250000000 0.15066964 0.25111607
#> 85      TRUE      1        FALSE -0.500000000 0.06993007 0.20979021
#> 86      TRUE      1        FALSE -0.500000000 0.06993007 0.20979021
#> 87      TRUE      1        FALSE -0.136363636 0.07916667 0.10416667
#> 88      TRUE      1        FALSE -0.136363636 0.07916667 0.10416667
#> 89      TRUE      1        FALSE -1.000000000 0.00000000 0.06666667
#> 90      TRUE      1        FALSE  0.190000000 0.28605769 0.19471154
#> 91      TRUE      1        FALSE -1.000000000 0.00000000 0.06666667
#> 92      TRUE      1        FALSE -1.000000000 0.00000000 0.07692308
#> 93      TRUE      1        FALSE  0.190000000 0.28605769 0.19471154
#> 94      TRUE      1        FALSE -0.238095238 0.10378510 0.16865079
#> 95      TRUE      1        FALSE  0.190476190 0.17045455 0.11590909
#> 96      TRUE      1        FALSE -1.000000000 0.00000000 0.07692308
#> 97      TRUE      1        FALSE -0.238095238 0.10378510 0.16865079
#> 98      TRUE      1        FALSE  0.190476190 0.17045455 0.11590909
#> 99      TRUE      1        FALSE  0.392857143 0.13493304 0.05881696
#> 100     TRUE      1        FALSE  0.392857143 0.13493304 0.05881696
#> 101     TRUE      1        FALSE  0.014084507 0.12662338 0.12310606
#> 102     TRUE      1        FALSE -0.008333333 0.19318182 0.19642857
#> 103     TRUE      1        FALSE -0.481481481 0.11666667 0.33333333
#> 104     TRUE      1        FALSE  0.014084507 0.12662338 0.12310606
#> 105     TRUE      1        FALSE -0.008333333 0.19318182 0.19642857
#> 106     TRUE      1        FALSE -0.481481481 0.11666667 0.33333333
#> 107     TRUE      1        FALSE          NaN 0.00000000 0.00000000
#> 108     TRUE      1        FALSE          NaN 0.00000000 0.00000000
#> 109     TRUE      1        FALSE  0.076923077 0.14358974 0.12307692
#> 110     TRUE      1        FALSE  0.076923077 0.14358974 0.12307692
#> 111     TRUE      1        FALSE  0.243243243 0.16770833 0.10208333
#> 112     TRUE      1        FALSE  0.243243243 0.16770833 0.10208333
#> 113     TRUE      1        FALSE  0.361702128 0.18648019 0.08741259
#> 114     TRUE      1        FALSE  1.000000000 0.06250000 0.00000000
#> 115     TRUE      1        FALSE  0.361702128 0.18648019 0.08741259
#> 116     TRUE      1        FALSE  0.142857143 0.26223776 0.19667832
#> 117     TRUE      1        FALSE  1.000000000 0.06250000 0.00000000
#> 118     TRUE      1        FALSE  0.225806452 0.29081633 0.18367347
#> 119     TRUE      1        FALSE  0.142857143 0.26223776 0.19667832
#> 120     TRUE      1        FALSE  0.225806452 0.29081633 0.18367347
#> 121     TRUE      1        FALSE -1.000000000 0.00000000 0.08333333
#> 122     TRUE      1        FALSE -0.038461538 0.08928571 0.09642857
#> 123     TRUE      1        FALSE -1.000000000 0.00000000 0.08333333
#> 124     TRUE      1        FALSE -0.038461538 0.08928571 0.09642857
#> 125     TRUE      1        FALSE -0.488721805 0.08116883 0.23634454
#> 126     TRUE      1        FALSE -0.107382550 0.08636364 0.10714286
#> 127     TRUE      1        FALSE  0.239130435 0.05480769 0.03365385
#> 128     TRUE      1        FALSE  0.153846154 0.09915865 0.07271635
#> 129     TRUE      1        FALSE -0.488721805 0.08116883 0.23634454
#> 130     TRUE      1        FALSE -0.107382550 0.08636364 0.10714286
#> 131     TRUE      1        FALSE -1.000000000 0.00000000 0.02820513
#> 132     TRUE      1        FALSE  0.239130435 0.05480769 0.03365385
#> 133     TRUE      1        FALSE  0.153846154 0.09915865 0.07271635
#> 134     TRUE      1        FALSE -1.000000000 0.00000000 0.02820513
#> 135     TRUE      1        FALSE -0.250000000 0.04326923 0.07211538
#> 136     TRUE      1        FALSE -0.250000000 0.04326923 0.07211538
#> 137     TRUE      1        FALSE          NaN 0.00000000 0.00000000
#> 138     TRUE      2        FALSE  0.133333333 0.10309480 0.07883720
#> 139     TRUE      1        FALSE  0.133333333 0.24791667 0.18958333
#> 140     TRUE      1        FALSE -0.200000000 0.03446970 0.05170455
#> 141     TRUE      1        FALSE  0.150000000 0.13269231 0.09807692
#> 142     TRUE      1        FALSE -0.200000000 0.03446970 0.05170455
#> 143     TRUE      1        FALSE  0.125000000 0.15750000 0.12250000
#> 144     TRUE      1        FALSE -0.187500000 0.06416667 0.09378205
#> 145     TRUE      1        FALSE  0.150000000 0.13269231 0.09807692
#> 146     TRUE      1        FALSE -0.302325581 0.11458333 0.21388889
#> 147     TRUE      1        FALSE  0.125000000 0.15750000 0.12250000
#> 148     TRUE      1        FALSE -0.187500000 0.06416667 0.09378205
#> 149     TRUE      1        FALSE -0.302325581 0.11458333 0.21388889
#> 150     TRUE      1        FALSE          NaN 0.00000000 0.00000000
#> 151     TRUE      1        FALSE          NaN 0.00000000 0.00000000
#> 152     TRUE      1        FALSE          NaN 0.00000000 0.00000000
#> 153     TRUE      1        FALSE  0.294117647 0.16587302 0.09047619
#> 154     TRUE      1        FALSE -0.112149533 0.08304196 0.10402098
#> 155     TRUE      1        FALSE          NaN 0.00000000 0.00000000
#> 156     TRUE      1        FALSE  0.294117647 0.16587302 0.09047619
#> 157     TRUE      1        FALSE -0.128205128 0.11383929 0.14732143
#> 158     TRUE      1        FALSE  0.217391304 0.12165775 0.07820856
#> 159     TRUE      1        FALSE -0.363636364 0.08076923 0.17307692
#> 160     TRUE      1        FALSE -0.112149533 0.08304196 0.10402098
#> 161     TRUE      1        FALSE -0.128205128 0.11383929 0.14732143
#> 162     TRUE      1        FALSE  0.217391304 0.12165775 0.07820856
#> 163     TRUE      1        FALSE -0.363636364 0.08076923 0.17307692
#> 164     TRUE      1        FALSE  0.714285714 0.17532468 0.02922078
#> 165     TRUE      1        FALSE  0.714285714 0.17532468 0.02922078
#> 166     TRUE      1        FALSE          NaN 0.00000000 0.00000000
#> 167     TRUE      1        FALSE  0.416666667 0.03946429 0.01625000
#> 168     TRUE      1        FALSE  0.142857143 0.28571429 0.21428571
#> 169     TRUE      1        FALSE  0.239130435 0.16193182 0.09943182
#> 170     TRUE      1        FALSE          NaN 0.00000000 0.00000000
#> 171     TRUE      1        FALSE  0.416666667 0.03946429 0.01625000
#> 172     TRUE      1        FALSE  0.142857143 0.28571429 0.21428571
#> 173     TRUE      1        FALSE  0.239130435 0.16193182 0.09943182
#> 174     TRUE      1        FALSE          NaN 0.00000000 0.00000000
#> 175     TRUE      1        FALSE          NaN 0.00000000 0.00000000
#> 176     TRUE      1        FALSE -0.473684211 0.06944444 0.19444444
#> 177     TRUE      1        FALSE -0.473684211 0.06944444 0.19444444
#> 178     TRUE      1        FALSE -0.470588235 0.09709821 0.26971726
#> 179     TRUE      1        FALSE  0.309090909 0.12587413 0.06643357
#> 180     TRUE      1        FALSE -0.470588235 0.09709821 0.26971726
#> 181     TRUE      1        FALSE  0.309090909 0.12587413 0.06643357
#> 182     TRUE      1        FALSE  0.600000000 0.14545455 0.03636364
#> 183     TRUE      1        FALSE  0.600000000 0.14545455 0.03636364
#> 184     TRUE      1        FALSE -0.357142857 0.08766234 0.18506494
#> 185     TRUE      1        FALSE -0.357142857 0.08766234 0.18506494
#> 186     TRUE      1        FALSE -1.000000000 0.00000000 0.07005495
#> 187     TRUE      1        FALSE -1.000000000 0.00000000 0.07005495
#> 188     TRUE      1        FALSE  0.111111111 0.16666667 0.13333333
#> 189     TRUE      1        FALSE  0.111111111 0.16666667 0.13333333
#> 190     TRUE      1        FALSE -0.111111111 0.05769231 0.07211538
#> 191     TRUE      1        FALSE -0.111111111 0.05769231 0.07211538
#> 192     TRUE      1        FALSE -1.000000000 0.00000000 0.03333333
#> 193     TRUE      1        FALSE -1.000000000 0.00000000 0.03333333
#> 194     TRUE      1        FALSE -0.298969072 0.08627137 0.15985577
#> 195     TRUE      1        FALSE -0.298969072 0.08627137 0.15985577
#> 196     TRUE      1        FALSE  0.100000000 0.27500000 0.22500000
#> 197     TRUE      1        FALSE  0.100000000 0.27500000 0.22500000
#> 198     TRUE      1        FALSE -0.014925373 0.16071429 0.16558442
#> 199     TRUE      1        FALSE -0.014925373 0.16071429 0.16558442
#> 
#> $weighted.means
#>    facet    subfacet snp.facet snp.subfacet D_abba_baba_jackknife
#> 1    pop ASP~UPD~PAL       chr       groupI                    NA
#> 2    pop ASP~UPD~PAL       chr      groupII                    NA
#> 3    pop ASP~UPD~PAL       chr     groupIII                    NA
#> 4    pop ASP~UPD~PAL       chr      groupIV                    NA
#> 5    pop ASP~UPD~PAL       chr      groupIX                    NA
#> 6    pop ASP~UPD~PAL       chr       groupV                    NA
#> 7    pop ASP~UPD~PAL       chr      groupVI                    NA
#> 8    pop ASP~UPD~PAL       chr     groupVII                    NA
#> 9    pop ASP~UPD~PAL       chr    groupVIII                    NA
#> 10   pop ASP~UPD~PAL       chr       groupX                    NA
#> 11   pop ASP~UPD~PAL       chr      groupXI                    NA
#> 12   pop ASP~UPD~PAL       chr     groupXII                    NA
#> 13   pop ASP~UPD~PAL       chr    groupXIII                    NA
#> 14   pop ASP~UPD~PAL       chr     groupXIV                    NA
#> 15   pop ASP~UPD~PAL       chr     groupXIX                    NA
#> 16   pop ASP~UPD~PAL       chr      groupXV                    NA
#> 17   pop ASP~UPD~PAL       chr     groupXVI                    NA
#> 18   pop ASP~UPD~PAL       chr    groupXVII                    NA
#> 19   pop ASP~UPD~PAL       chr   groupXVIII                    NA
#> 20   pop ASP~UPD~PAL       chr      groupXX                    NA
#> 21   pop ASP~UPD~PAL       chr     groupXXI                    NA
#>    D_abba_baba_p_value  D_abba_baba
#> 1                   NA  0.109856907
#> 2                   NA -0.211599173
#> 3                   NA -0.048932585
#> 4                   NA  0.135987882
#> 5                   NA -0.008755932
#> 6                   NA -0.250000000
#> 7                   NA  0.035554682
#> 8                   NA -0.157011913
#> 9                   NA -0.043345106
#> 10                  NA -0.552554938
#> 11                  NA  0.367938040
#> 12                  NA -0.301593120
#> 13                  NA -0.287231988
#> 14                  NA  0.239130435
#> 15                  NA  0.168792465
#> 16                  NA -0.012421169
#> 17                  NA -0.212902729
#> 18                  NA -0.009583655
#> 19                  NA  0.361501281
#> 20                  NA -0.227296543
#> 21                  NA  0.062500000
#> 
```
