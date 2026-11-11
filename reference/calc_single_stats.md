# Calculate standard single-SNP genetic statistics

These functions calculate a basic genetic statistics from SNP data
contained in snpRdata objects, splitting samples by any number of
provided facets. Since these statistics all use only a single SNP at a
time, they ignore any SNP specific facet levels.

## Usage

``` r
calc_pi(x, facets = NULL)

calc_maf(x, facets = NULL)

calc_ho(x, facets = NULL)

calc_private(x, facets = NULL, rarefaction = TRUE, g = 0)

calc_hwe(
  x,
  facets = NULL,
  method = "exact",
  fwe_method = "BY",
  fwe_case = c("by_facet", "overall")
)

calc_he(x, facets = NULL)

calc_allelic_richness(x, facets = NULL, g = 0)
```

## Arguments

- x:

  snpRdata. Input SNP data.

- facets:

  character. Categorical metadata variables by which to break up
  analysis. See
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md)
  for more details.

- rarefaction:

  logical, default TRUE. Should the number of segregating sites be
  estimated via rarefaction? See details.

- g:

  numeric, default 0. If doing rarefaction, controls the number of
  *alleles/gene copies* to rarefact to. If 0, this will rarefact to the
  smallest sample size per locus. If g \< 0, this will rarefact to to
  the smallest sample size per locus minus the absolute value of g. If
  positive, this will rarefact to g, and any loci where the smallest
  sample size is less than g will be dropped from the calculation.

- method:

  character, default "exact". Defines the method to use for calculating
  p-values for HWE divergence. Options:

  - exact: Uses the exact test as described in Wigginton et al (2005).

  - chisq: Uses a chi-squared test.

  See details

- fwe_method:

  character, default "BY". Type of Family-Wise Error correction
  (multiple testing correction) to use. For details and options, see
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). If no correction
  is desired, set this argument to "none".

- fwe_case:

  character, default c("by_facet", "by_subfacet", "overall"). How should
  Family-Wise Error correction (multiple testing correction) be applied?

  - "by_facet": Each facet supplied (such as pop or pop.fam) is treated
    as a set of tests.

  - "by_subfacet": Each level of each subfacet is treated as a separate
    set of tests.

  - "overall": All tests are treated as a set.

## Value

snpRdata object with requested stats merged into the stats socket

## Details

The data can be broken up categorically by sample metadata, as described
in
[`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).

## Functions

- `calc_pi()`: \\\pi\\ (nucleotide diversity/average number of pairwise
  differences)

- `calc_maf()`: minor allele frequency

- `calc_ho()`: observed heterozygosity

- `calc_private()`: find private alleles

- `calc_hwe()`: p-values for Hardy-Weinberg Equilibrium divergence

- `calc_he()`: expected heterozygosity

- `calc_allelic_richness()`: allelic richness (standardized number of
  alleles per locus via rarefaction)

## \\\pi\\

Calculates \\\pi\\ (nucleotide diversity/average number of pairwise
differences) according to Hohenlohe et al. (2010).

## H_(E)

Calculates traditional expected heterozygosity \\2pq\\. Note that this
will produce results almost identical to \\\pi\\.

## H_(O)

Calculates observed heterozygosity.

## maf

Calculates minor allele frequencies and note identities and counts of
major and minor alleles.

## private alleles

Determines if each SNP is a private allele across all levels in each
sample facet. Will return an error if no sample facets are provided. If
rarefaction is requested, the estimated number of private alleles will
be calculated according to Smith and Grassle (1977). Note that the
standardized sample size (*g*) will vary across loci due to differences
in sequencing coverage at those loci, equal to the smallest number of
alleles sequenced in any population at that locus minus one. Instead of
weighted averages, the value stored in the `$weighted.means` slot in the
returned value is the total number of private alleles per population.

## hwe

Calculates a p-value for the null hypothesis that a population is in HWE
at a given locus. Several methods available:

- "exact": Exact test according to Wigginton, JE, Cutler, DJ, and
  Abecasis, GR (2005). Slightly slower.

- "chisq": Chi-squared test. May produce poor results when sample sizes
  for any observed or expected genotypes are small.

For the exact test, code derived from
<http://csg.sph.umich.edu/abecasis/Exact/snp_hwe.r>

## allelic richness

Calculates the allelic richness, the estimated number of alleles per
locus standardized via rarefaction for sample size according to Hurlburt
(1971). Note that the standardized sample size (*g*) will vary across
loci due to differences in sequencing coverage at those loci, equal to
the smallest number of alleles sequenced in any population at that locus
minus one. Weighted averages are weighted by *g*.

## References

Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005). *American Journal of
Human Genetics*

Hohenlohe et al. (2010). *PLOS Genetics*.

Hurlburt (1971). *Ecology*.

Smith and Grassle (1977). *Biometrics*

## Author

William Hemstrom

## Examples

``` r
# base facet
x <- calc_pi(stickSNPs)
get.snpR.stats(x)
#>     facet subfacet        chr position .snp.id major minor maj.count min.count
#> 1   .base    .base    groupVI   212436       1     G     A        87        85
#> 2   .base    .base groupXVIII   350070       2     G     T       180         6
#> 3   .base    .base    groupXV   854285       3     C     T       179        13
#> 4   .base    .base  groupXIII   907279       4     A     G       129        47
#> 5   .base    .base    groupVI  1068620       5     G     A        85        33
#> 6   .base    .base     groupV  1392753       6     A     C       149        37
#> 7   .base    .base    groupIX  2075079       7     T     G       163        17
#> 8   .base    .base  groupXIII  2106000       8     G     A       133        37
#> 9   .base    .base   groupXII  2136587       9     C     A       152        30
#> 10  .base    .base  groupXVII  2560839      10     A     G       180        10
#> 11  .base    .base   groupIII  2750695      11     T     C       102        70
#> 12  .base    .base   groupXIX  3008524      12     A     G       119        35
#> 13  .base    .base   groupXVI  3624417      13     T     C       142        32
#> 14  .base    .base   groupVII  3667528      14     G     T       161        29
#> 15  .base    .base   groupVII  4189120      15     A     G       112        62
#> 16  .base    .base    groupXX  4538074      16     G     T       112        52
#> 17  .base    .base  groupXIII  4706089      17     A     C       155        27
#> 18  .base    .base   groupIII  4746808      18     C     A       127        61
#> 19  .base    .base  groupXVII  5041995      19     T     G       173        11
#> 20  .base    .base     groupX  5161637      20     G     A       174         8
#> 21  .base    .base   groupXII  5178037      21     C     A       114        46
#> 22  .base    .base    groupIX  5287711      22     C     T       131        41
#> 23  .base    .base  groupVIII  5462135      23     G     A       176         8
#> 24  .base    .base    groupIX  5960933      24     A     G       171        11
#> 25  .base    .base  groupXIII  6076868      25     G     A       125        51
#> 26  .base    .base    groupXI  6177332      26     T     C       122        48
#> 27  .base    .base    groupIV  6291287      27     A     C       162        20
#> 28  .base    .base    groupIX  6515796      28     T     C       103        61
#> 29  .base    .base   groupXII  6718079      29     A     C       126        28
#> 30  .base    .base   groupVII  6967477      30     G     C       160        16
#> 31  .base    .base   groupIII  7194069      31     C     G       136        50
#> 32  .base    .base   groupXVI  7199722      32     A     C       135        37
#> 33  .base    .base  groupVIII  7224065      33     C     T       114        54
#> 34  .base    .base groupXVIII  7229426      34     C     G       123        43
#> 35  .base    .base groupXVIII  8087022      35     A     C       164        10
#> 36  .base    .base    groupIX  8373442      36     A     G       158        24
#> 37  .base    .base    groupXV  8815538      37     G     C       137        43
#> 38  .base    .base    groupXI  9019877      38     G     T       122        54
#> 39  .base    .base   groupXXI  9067876      39     T     A       113        17
#> 40  .base    .base  groupXVII  9108184      40     G     C       122        38
#> 41  .base    .base    groupIX  9284427      41     G     A        75        75
#> 42  .base    .base     groupV  9767769      42     G     A       110        66
#> 43  .base    .base   groupXII  9881312      43     A     G       159        23
#> 44  .base    .base    groupXX 10296516      44     T     A       150        28
#> 45  .base    .base     groupI 10593829      45     G     T       174         2
#> 46  .base    .base   groupXIX 10765968      46     C     A       114        40
#> 47  .base    .base    groupII 10829169      47     C     G       140        46
#> 48  .base    .base    groupXX 10913115      48     T     C        89        71
#> 49  .base    .base groupXVIII 10974594      49     G     C       189         5
#> 50  .base    .base  groupXVII 11256661      50     C     T       160        14
#> 51  .base    .base  groupVIII 11438870      51     T     C       101        45
#> 52  .base    .base   groupVII 11477223      52     A     G       126        56
#> 53  .base    .base    groupII 11505084      53     A     C       134        38
#> 54  .base    .base    groupIX 11904651      54     G     T       191         3
#> 55  .base    .base  groupVIII 12108945      55     G     A       108        58
#> 56  .base    .base   groupXIX 12485182      56     A     G       130        40
#> 57  .base    .base    groupXI 12833733      57     G     A       140         6
#> 58  .base    .base groupXVIII 12937786      58     C     T       121        43
#> 59  .base    .base   groupIII 13177027      59     G     A       116        64
#> 60  .base    .base    groupVI 13335577      60     A     C       122        54
#> 61  .base    .base     groupX 13718292      61     T     C       140         6
#> 62  .base    .base    groupVI 13895716      62     G     C       173        19
#> 63  .base    .base   groupXIV 14009990      63     A     C        89        51
#> 64  .base    .base   groupXII 14158350      64     T     A       113        65
#> 65  .base    .base    groupIX 14164215      65     C     T       110        68
#> 66  .base    .base    groupVI 14236496      66     A     C        89        73
#> 67  .base    .base     groupX 14273796      67     G     A       171         9
#> 68  .base    .base  groupVIII 14445389      68     A     T       112        58
#> 69  .base    .base    groupIX 14868533      69     A     G       166        20
#> 70  .base    .base    groupIX 15190693      70     A     C       130        42
#> 71  .base    .base   groupXII 15306362      71     A     T       182        16
#> 72  .base    .base    groupXX 15460430      72     A     G        97        83
#> 73  .base    .base   groupVII 15493025      73     T     G       142        36
#> 74  .base    .base    groupXV 15547073      74     G     T       135        27
#> 75  .base    .base  groupXIII 15697285      75     A     C        91        85
#> 76  .base    .base   groupIII 15983994      76     C     T        70        48
#> 77  .base    .base   groupXIX 16233699      77     A     T       133        39
#> 78  .base    .base    groupIX 16306966      78     G     T       163        21
#> 79  .base    .base     groupI 16498302      79     T     C       163        25
#> 80  .base    .base  groupXIII 16520929      80     T     C       125        55
#> 81  .base    .base   groupXVI 16561891      81     A     T       106        68
#> 82  .base    .base    groupII 16598364      82     T     C       162        14
#> 83  .base    .base    groupIV 17248473      83     T     A       183        11
#> 84  .base    .base    groupXX 17748944      84     C     T        90        76
#> 85  .base    .base  groupXIII 17842228      85     G     A        91        81
#> 86  .base    .base   groupXIX 17858769      86     A     G       115        39
#> 87  .base    .base    groupIX 17929225      87     T     C       140        48
#> 88  .base    .base    groupIX 18339767      88     A     G        79        73
#> 89  .base    .base  groupXIII 18836558      89     C     G        84        60
#> 90  .base    .base    groupII 19244005      90     C     T       186         6
#> 91  .base    .base    groupXX 19317273      91     C     T       106        62
#> 92  .base    .base     groupI 19759331      92     G     A        75        41
#> 93  .base    .base   groupVII 20165240      93     T     C       153        27
#> 94  .base    .base    groupII 20496069      94     C     G       152        16
#> 95  .base    .base     groupI 21284266      95     G     T       146        26
#> 96  .base    .base    groupIV 22227311      96     C     A       164        16
#> 97  .base    .base     groupI 22746027      97     G     A       174         6
#> 98  .base    .base   groupVII 25239495      98     G     T       139        43
#> 99  .base    .base   groupVII 27599144      99     G     A        92        72
#> 100 .base    .base    groupIV 32321895     100     A     T        94        64
#>            maf         pi
#> 1   0.49418605 0.50285598
#> 2   0.03225806 0.06277245
#> 3   0.06770833 0.12690881
#> 4   0.26704545 0.39370130
#> 5   0.27966102 0.40634507
#> 6   0.19892473 0.32043011
#> 7   0.09444444 0.17200497
#> 8   0.21764706 0.34256874
#> 9   0.16483516 0.27685022
#> 10  0.05263158 0.10025063
#> 11  0.40697674 0.48551612
#> 12  0.22727273 0.35353535
#> 13  0.18390805 0.30190685
#> 14  0.15263158 0.26003899
#> 15  0.35632184 0.46136469
#> 16  0.31707317 0.43573246
#> 17  0.14835165 0.25408293
#> 18  0.32446809 0.44072136
#> 19  0.05978261 0.11303160
#> 20  0.04395604 0.08451217
#> 21  0.28750000 0.41226415
#> 22  0.23837209 0.36522508
#> 23  0.04347826 0.08363032
#> 24  0.06043956 0.11420072
#> 25  0.28977273 0.41396104
#> 26  0.28235294 0.40765750
#> 27  0.10989011 0.19670937
#> 28  0.37195122 0.47007332
#> 29  0.18181818 0.29946524
#> 30  0.09090909 0.16623377
#> 31  0.26881720 0.39523394
#> 32  0.21511628 0.33965728
#> 33  0.32142857 0.43883661
#> 34  0.25903614 0.38619934
#> 35  0.05747126 0.10896286
#> 36  0.13186813 0.23022282
#> 37  0.23888889 0.36567349
#> 38  0.30681818 0.42779221
#> 39  0.13076923 0.22909958
#> 40  0.23750000 0.36446541
#> 41  0.50000000 0.50335570
#> 42  0.37500000 0.47142857
#> 43  0.12637363 0.22202659
#> 44  0.15730337 0.26661588
#> 45  0.01136364 0.02259740
#> 46  0.25974026 0.38706392
#> 47  0.24731183 0.37430979
#> 48  0.44375000 0.49677673
#> 49  0.02577320 0.05047807
#> 50  0.08045977 0.14882732
#> 51  0.30821918 0.42938120
#> 52  0.30769231 0.42838929
#> 53  0.22093023 0.34625323
#> 54  0.01546392 0.03060734
#> 55  0.34939759 0.45739321
#> 56  0.23529412 0.36199095
#> 57  0.04109589 0.07935758
#> 58  0.26219512 0.38927129
#> 59  0.35555556 0.46083178
#> 60  0.30681818 0.42779221
#> 61  0.04109589 0.07935758
#> 62  0.09895833 0.17926483
#> 63  0.36428571 0.46649538
#> 64  0.36516854 0.46626039
#> 65  0.38202247 0.47483019
#> 66  0.45061728 0.49819799
#> 67  0.05000000 0.09553073
#> 68  0.34117647 0.45221023
#> 69  0.10752688 0.19296716
#> 70  0.24418605 0.37127703
#> 71  0.08080808 0.14931036
#> 72  0.46111111 0.49975171
#> 73  0.20224719 0.32450962
#> 74  0.16666667 0.27950311
#> 75  0.48295455 0.50227273
#> 76  0.40677966 0.48674489
#> 77  0.22674419 0.35271318
#> 78  0.11413043 0.20331433
#> 79  0.13297872 0.23182387
#> 80  0.30555556 0.42675357
#> 81  0.39080460 0.47890506
#> 82  0.07954545 0.14727273
#> 83  0.05670103 0.10752631
#> 84  0.45783133 0.49945235
#> 85  0.47093023 0.50122399
#> 86  0.25324675 0.38069773
#> 87  0.25531915 0.38229605
#> 88  0.48026316 0.50252701
#> 89  0.41666667 0.48951049
#> 90  0.03125000 0.06086387
#> 91  0.36904762 0.46849159
#> 92  0.35344828 0.46101949
#> 93  0.15000000 0.25642458
#> 94  0.09523810 0.17336755
#> 95  0.15116279 0.25812593
#> 96  0.08888889 0.16288020
#> 97  0.03333333 0.06480447
#> 98  0.23626374 0.36288021
#> 99  0.43902439 0.49558581
#> 100 0.40506329 0.48504394

# multiple facets
x <- calc_pi(stickSNPs, facets = c("pop", "pop.fam"))
get.snpR.stats(x, c("pop", "pop.fam"))
#>        facet subfacet        chr position .snp.id         pi        maf
#> 1    fam.pop    A.ASP    groupVI   212436       1 0.50000000 0.75000000
#> 2    fam.pop    A.CLF    groupVI   212436       1 0.50000000 0.75000000
#> 3    fam.pop    A.OPL    groupVI   212436       1 0.53030303 0.58333333
#> 4    fam.pop    A.PAL    groupVI   212436       1 0.53333333 0.40000000
#> 5    fam.pop    A.SMR    groupVI   212436       1 0.54545455 0.50000000
#> 6    fam.pop    A.UPD    groupVI   212436       1 0.53333333 0.33333333
#> 7    fam.pop    B.ASP    groupVI   212436       1 0.33333333 0.16666667
#> 8    fam.pop    B.CLF    groupVI   212436       1 1.00000000 0.50000000
#> 9    fam.pop    B.OPL    groupVI   212436       1 0.53333333 0.60000000
#> 10   fam.pop    B.PAL    groupVI   212436       1 0.00000000 0.00000000
#> 11   fam.pop    B.SMR    groupVI   212436       1 0.46666667 0.70000000
#> 12   fam.pop    B.UPD    groupVI   212436       1 0.60000000 0.50000000
#> 13   fam.pop    C.ASP    groupVI   212436       1 0.50000000 0.75000000
#> 14   fam.pop    C.CLF    groupVI   212436       1 0.00000000 1.00000000
#> 15   fam.pop    C.OPL    groupVI   212436       1 0.53571429 0.62500000
#> 16   fam.pop    C.PAL    groupVI   212436       1 0.42857143 0.25000000
#> 17   fam.pop    C.SMR    groupVI   212436       1 0.53333333 0.40000000
#> 18   fam.pop    C.UPD    groupVI   212436       1 0.53571429 0.62500000
#> 19   fam.pop    D.ASP    groupVI   212436       1 0.66666667 0.50000000
#> 20   fam.pop    D.CLF    groupVI   212436       1 1.00000000 0.50000000
#> 21   fam.pop    D.OPL    groupVI   212436       1 0.55555556 0.50000000
#> 22   fam.pop    D.PAL    groupVI   212436       1 0.25000000 0.12500000
#> 23   fam.pop    D.SMR    groupVI   212436       1 0.40909091 0.75000000
#> 24   fam.pop    D.UPD    groupVI   212436       1 0.60000000 0.50000000
#> 25       pop      ASP    groupVI   212436       1 0.52941176 0.50000000
#> 26       pop      CLF    groupVI   212436       1 0.46666667 0.70000000
#> 27       pop      OPL    groupVI   212436       1 0.50128205 0.57500000
#> 28       pop      PAL    groupVI   212436       1 0.33689840 0.20588235
#> 29       pop      SMR    groupVI   212436       1 0.49471459 0.59090909
#> 30       pop      UPD    groupVI   212436       1 0.52000000 0.50000000
#> 31   fam.pop    A.ASP groupXVIII   350070       2 0.33333333 0.16666667
#> 32   fam.pop    A.CLF groupXVIII   350070       2 0.00000000 0.00000000
#> 33   fam.pop    A.OPL groupXVIII   350070       2 0.00000000 0.00000000
#> 34   fam.pop    A.PAL groupXVIII   350070       2 0.00000000 0.00000000
#> 35   fam.pop    A.SMR groupXVIII   350070       2 0.00000000 0.00000000
#> 36   fam.pop    A.UPD groupXVIII   350070       2 0.00000000 0.00000000
#> 37   fam.pop    B.ASP groupXVIII   350070       2 0.00000000 0.00000000
#> 38   fam.pop    B.CLF groupXVIII   350070       2 0.00000000 0.00000000
#> 39   fam.pop    B.OPL groupXVIII   350070       2 0.20000000 0.10000000
#> 40   fam.pop    B.PAL groupXVIII   350070       2 0.00000000 0.00000000
#> 41   fam.pop    B.SMR groupXVIII   350070       2 0.20000000 0.10000000
#> 42   fam.pop    B.UPD groupXVIII   350070       2 0.00000000 0.00000000
#> 43   fam.pop    C.ASP groupXVIII   350070       2 0.00000000 0.00000000
#> 44   fam.pop    C.CLF groupXVIII   350070       2 0.00000000 0.00000000
#> 45   fam.pop    C.OPL groupXVIII   350070       2 0.00000000 0.00000000
#> 46   fam.pop    C.PAL groupXVIII   350070       2 0.00000000 0.00000000
#> 47   fam.pop    C.SMR groupXVIII   350070       2 0.16666667 0.08333333
#> 48   fam.pop    C.UPD groupXVIII   350070       2 0.00000000 0.00000000
#> 49   fam.pop    D.ASP groupXVIII   350070       2 0.00000000 0.00000000
#> 50   fam.pop    D.CLF groupXVIII   350070       2 0.00000000 0.00000000
#> 51   fam.pop    D.OPL groupXVIII   350070       2 0.30303030 0.16666667
#> 52   fam.pop    D.PAL groupXVIII   350070       2 0.00000000 0.00000000
#> 53   fam.pop    D.SMR groupXVIII   350070       2 0.00000000 0.00000000
#> 54   fam.pop    D.UPD groupXVIII   350070       2 0.00000000 0.00000000
#> 55       pop      ASP groupXVIII   350070       2 0.10000000 0.05000000
#> 56       pop      CLF groupXVIII   350070       2 0.00000000 0.00000000
#> 57       pop      OPL groupXVIII   350070       2 0.12463768 0.06521739
#> 58       pop      PAL groupXVIII   350070       2 0.00000000 0.00000000
#> 59       pop      SMR groupXVIII   350070       2 0.08879493 0.04545455
#> 60       pop      UPD groupXVIII   350070       2 0.00000000 0.00000000
#> 61   fam.pop    A.ASP    groupXV   854285       3 0.00000000 0.00000000
#> 62   fam.pop    A.CLF    groupXV   854285       3 0.50000000 0.25000000
#> 63   fam.pop    A.OPL    groupXV   854285       3 0.00000000 0.00000000
#> 64   fam.pop    A.PAL    groupXV   854285       3 0.00000000 0.00000000
#> 65   fam.pop    A.SMR    groupXV   854285       3 0.00000000 0.00000000
#> 66   fam.pop    A.UPD    groupXV   854285       3 0.00000000 0.00000000
#> 67   fam.pop    B.ASP    groupXV   854285       3 0.33333333 0.16666667
#> 68   fam.pop    B.CLF    groupXV   854285       3 0.00000000 0.00000000
#> 69   fam.pop    B.OPL    groupXV   854285       3 0.14285714 0.07142857
#> 70   fam.pop    B.PAL    groupXV   854285       3 0.25000000 0.12500000
#> 71   fam.pop    B.SMR    groupXV   854285       3 0.00000000 0.00000000
#> 72   fam.pop    B.UPD    groupXV   854285       3 0.25000000 0.12500000
#> 73   fam.pop    C.ASP    groupXV   854285       3 0.33333333 0.16666667
#> 74   fam.pop    C.CLF    groupXV   854285       3 0.00000000 0.00000000
#> 75   fam.pop    C.OPL    groupXV   854285       3 0.30303030 0.16666667
#> 76   fam.pop    C.PAL    groupXV   854285       3 0.00000000 0.00000000
#> 77   fam.pop    C.SMR    groupXV   854285       3 0.30303030 0.16666667
#> 78   fam.pop    C.UPD    groupXV   854285       3 0.25000000 0.12500000
#> 79   fam.pop    D.ASP    groupXV   854285       3 0.00000000 0.00000000
#> 80   fam.pop    D.CLF    groupXV   854285       3 0.00000000 0.00000000
#> 81   fam.pop    D.OPL    groupXV   854285       3 0.14285714 0.07142857
#> 82   fam.pop    D.PAL    groupXV   854285       3 0.00000000 0.00000000
#> 83   fam.pop    D.SMR    groupXV   854285       3 0.20000000 0.10000000
#> 84   fam.pop    D.UPD    groupXV   854285       3 0.00000000 0.00000000
#> 85       pop      ASP    groupXV   854285       3 0.17316017 0.09090909
#> 86       pop      CLF    groupXV   854285       3 0.16666667 0.08333333
#> 87       pop      OPL    groupXV   854285       3 0.15020408 0.08000000
#> 88       pop      PAL    groupXV   854285       3 0.06250000 0.03125000
#> 89       pop      SMR    groupXV   854285       3 0.12463768 0.06521739
#> 90       pop      UPD    groupXV   854285       3 0.12873563 0.06666667
#> 91   fam.pop    A.ASP  groupXIII   907279       4 0.33333333 0.16666667
#> 92   fam.pop    A.CLF  groupXIII   907279       4 0.50000000 0.25000000
#> 93   fam.pop    A.OPL  groupXIII   907279       4 0.57142857 0.50000000
#> 94   fam.pop    A.PAL  groupXIII   907279       4 0.00000000 0.00000000
#> 95   fam.pop    A.SMR  groupXIII   907279       4 0.40909091 0.25000000
#> 96   fam.pop    A.UPD  groupXIII   907279       4 0.66666667 0.50000000
#> 97   fam.pop    B.ASP  groupXIII   907279       4 0.33333333 0.16666667
#> 98   fam.pop    B.CLF  groupXIII   907279       4 0.00000000 0.00000000
#> 99   fam.pop    B.OPL  groupXIII   907279       4 0.46666667 0.30000000
#> 100  fam.pop    B.PAL  groupXIII   907279       4 0.00000000 0.00000000
#> 101  fam.pop    B.SMR  groupXIII   907279       4 0.55555556 0.50000000
#> 102  fam.pop    B.UPD  groupXIII   907279       4 0.53333333 0.33333333
#> 103  fam.pop    C.ASP  groupXIII   907279       4 0.33333333 0.16666667
#> 104  fam.pop    C.CLF  groupXIII   907279       4        NaN         NA
#> 105  fam.pop    C.OPL  groupXIII   907279       4 0.54545455 0.50000000
#> 106  fam.pop    C.PAL  groupXIII   907279       4 0.00000000 0.00000000
#> 107  fam.pop    C.SMR  groupXIII   907279       4 0.40909091 0.25000000
#> 108  fam.pop    C.UPD  groupXIII   907279       4 0.57142857 0.50000000
#> 109  fam.pop    D.ASP  groupXIII   907279       4 0.00000000 0.00000000
#> 110  fam.pop    D.CLF  groupXIII   907279       4 0.00000000 0.00000000
#> 111  fam.pop    D.OPL  groupXIII   907279       4 0.53333333 0.40000000
#> 112  fam.pop    D.PAL  groupXIII   907279       4 0.00000000 0.00000000
#> 113  fam.pop    D.SMR  groupXIII   907279       4 0.48484848 0.33333333
#> 114  fam.pop    D.UPD  groupXIII   907279       4 0.53571429 0.37500000
#> 115      pop      ASP  groupXIII   907279       4 0.24675325 0.13636364
#> 116      pop      CLF  groupXIII   907279       4 0.20000000 0.10000000
#> 117      pop      OPL  groupXIII   907279       4 0.50128205 0.42500000
#> 118      pop      PAL  groupXIII   907279       4 0.00000000 0.00000000
#> 119      pop      SMR  groupXIII   907279       4 0.44927536 0.32608696
#> 120      pop      UPD  groupXIII   907279       4 0.50769231 0.42307692
#> 121  fam.pop    A.ASP    groupVI  1068620       5 0.00000000 0.00000000
#> 122  fam.pop    A.CLF    groupVI  1068620       5 0.00000000 0.00000000
#> 123  fam.pop    A.OPL    groupVI  1068620       5 0.25000000 0.12500000
#> 124  fam.pop    A.PAL    groupVI  1068620       5 0.42857143 0.25000000
#> 125  fam.pop    A.SMR    groupVI  1068620       5 0.53571429 0.37500000
#> 126  fam.pop    A.UPD    groupVI  1068620       5 1.00000000 0.50000000
#> 127  fam.pop    B.ASP    groupVI  1068620       5 0.50000000 0.25000000
#> 128  fam.pop    B.CLF    groupVI  1068620       5        NaN         NA
#> 129  fam.pop    B.OPL    groupVI  1068620       5 0.53333333 0.66666667
#> 130  fam.pop    B.PAL    groupVI  1068620       5 0.50000000 0.25000000
#> 131  fam.pop    B.SMR    groupVI  1068620       5 0.53571429 0.37500000
#> 132  fam.pop    B.UPD    groupVI  1068620       5 0.53571429 0.37500000
#> 133  fam.pop    C.ASP    groupVI  1068620       5 1.00000000 0.50000000
#> 134  fam.pop    C.CLF    groupVI  1068620       5        NaN         NA
#> 135  fam.pop    C.OPL    groupVI  1068620       5 0.33333333 0.16666667
#> 136  fam.pop    C.PAL    groupVI  1068620       5 0.53333333 0.33333333
#> 137  fam.pop    C.SMR    groupVI  1068620       5 0.25000000 0.12500000
#> 138  fam.pop    C.UPD    groupVI  1068620       5 1.00000000 0.50000000
#> 139  fam.pop    D.ASP    groupVI  1068620       5 0.50000000 0.75000000
#> 140  fam.pop    D.CLF    groupVI  1068620       5 0.00000000 0.00000000
#> 141  fam.pop    D.OPL    groupVI  1068620       5 0.53333333 0.33333333
#> 142  fam.pop    D.PAL    groupVI  1068620       5 0.00000000 0.00000000
#> 143  fam.pop    D.SMR    groupVI  1068620       5 0.00000000 0.00000000
#> 144  fam.pop    D.UPD    groupVI  1068620       5 0.60000000 0.50000000
#> 145      pop      ASP    groupVI  1068620       5 0.49450549 0.35714286
#> 146      pop      CLF    groupVI  1068620       5 0.00000000 0.00000000
#> 147      pop      OPL    groupVI  1068620       5 0.44307692 0.30769231
#> 148      pop      PAL    groupVI  1068620       5 0.34420290 0.20833333
#> 149      pop      SMR    groupVI  1068620       5 0.37011494 0.23333333
#> 150      pop      UPD    groupVI  1068620       5 0.52287582 0.44444444
#> 151  fam.pop    A.ASP     groupV  1392753       6 0.33333333 0.16666667
#> 152  fam.pop    A.CLF     groupV  1392753       6 0.50000000 0.25000000
#> 153  fam.pop    A.OPL     groupV  1392753       6 0.16666667 0.08333333
#> 154  fam.pop    A.PAL     groupV  1392753       6 0.46666667 0.30000000
#> 155  fam.pop    A.SMR     groupV  1392753       6 0.00000000 0.00000000
#> 156  fam.pop    A.UPD     groupV  1392753       6 0.53333333 0.33333333
#> 157  fam.pop    B.ASP     groupV  1392753       6 0.60000000 0.50000000
#> 158  fam.pop    B.CLF     groupV  1392753       6 1.00000000 0.50000000
#> 159  fam.pop    B.OPL     groupV  1392753       6 0.14285714 0.07142857
#> 160  fam.pop    B.PAL     groupV  1392753       6 0.53571429 0.62500000
#> 161  fam.pop    B.SMR     groupV  1392753       6 0.35555556 0.20000000
#> 162  fam.pop    B.UPD     groupV  1392753       6 0.57142857 0.50000000
#> 163  fam.pop    C.ASP     groupV  1392753       6 0.00000000 0.00000000
#> 164  fam.pop    C.CLF     groupV  1392753       6 0.00000000 0.00000000
#> 165  fam.pop    C.OPL     groupV  1392753       6 0.00000000 0.00000000
#> 166  fam.pop    C.PAL     groupV  1392753       6 0.53571429 0.37500000
#> 167  fam.pop    C.SMR     groupV  1392753       6 0.00000000 0.00000000
#> 168  fam.pop    C.UPD     groupV  1392753       6 0.53571429 0.37500000
#> 169  fam.pop    D.ASP     groupV  1392753       6 0.50000000 0.25000000
#> 170  fam.pop    D.CLF     groupV  1392753       6 0.50000000 0.25000000
#> 171  fam.pop    D.OPL     groupV  1392753       6 0.16666667 0.08333333
#> 172  fam.pop    D.PAL     groupV  1392753       6 0.53333333 0.33333333
#> 173  fam.pop    D.SMR     groupV  1392753       6 0.20000000 0.10000000
#> 174  fam.pop    D.UPD     groupV  1392753       6 0.33333333 0.16666667
#> 175      pop      ASP     groupV  1392753       6 0.39473684 0.25000000
#> 176      pop      CLF     groupV  1392753       6 0.40909091 0.25000000
#> 177      pop      OPL     groupV  1392753       6 0.11085973 0.05769231
#> 178      pop      PAL     groupV  1392753       6 0.49798387 0.40625000
#> 179      pop      SMR     groupV  1392753       6 0.13588850 0.07142857
#> 180      pop      UPD     groupV  1392753       6 0.47619048 0.35714286
#> 181  fam.pop    A.ASP    groupIX  2075079       7 0.33333333 0.16666667
#> 182  fam.pop    A.CLF    groupIX  2075079       7 0.00000000 0.00000000
#> 183  fam.pop    A.OPL    groupIX  2075079       7 0.20000000 0.10000000
#> 184  fam.pop    A.PAL    groupIX  2075079       7 0.00000000 0.00000000
#> 185  fam.pop    A.SMR    groupIX  2075079       7 0.46666667 0.30000000
#> 186  fam.pop    A.UPD    groupIX  2075079       7 0.00000000 0.00000000
#> 187  fam.pop    B.ASP    groupIX  2075079       7 0.00000000 0.00000000
#> 188  fam.pop    B.CLF    groupIX  2075079       7 0.00000000 0.00000000
#> 189  fam.pop    B.OPL    groupIX  2075079       7 0.30303030 0.16666667
#> 190  fam.pop    B.PAL    groupIX  2075079       7 0.00000000 0.00000000
#> 191  fam.pop    B.SMR    groupIX  2075079       7 0.35555556 0.20000000
#> 192  fam.pop    B.UPD    groupIX  2075079       7 0.00000000 0.00000000
#> 193  fam.pop    C.ASP    groupIX  2075079       7 0.00000000 0.00000000
#> 194  fam.pop    C.CLF    groupIX  2075079       7 0.00000000 0.00000000
#> 195  fam.pop    C.OPL    groupIX  2075079       7 0.20000000 0.10000000
#> 196  fam.pop    C.PAL    groupIX  2075079       7 0.00000000 0.00000000
#> 197  fam.pop    C.SMR    groupIX  2075079       7 0.35555556 0.20000000
#> 198  fam.pop    C.UPD    groupIX  2075079       7 0.00000000 0.00000000
#> 199  fam.pop    D.ASP    groupIX  2075079       7 0.00000000 0.00000000
#> 200  fam.pop    D.CLF    groupIX  2075079       7 0.66666667 0.50000000
#> 201  fam.pop    D.OPL    groupIX  2075079       7 0.35555556 0.20000000
#> 202  fam.pop    D.PAL    groupIX  2075079       7 0.00000000 0.00000000
#> 203  fam.pop    D.SMR    groupIX  2075079       7 0.16666667 0.08333333
#> 204  fam.pop    D.UPD    groupIX  2075079       7 0.00000000 0.00000000
#> 205      pop      ASP    groupIX  2075079       7 0.09090909 0.04545455
#> 206      pop      CLF    groupIX  2075079       7 0.30303030 0.16666667
#> 207      pop      OPL    groupIX  2075079       7 0.25087108 0.14285714
#> 208      pop      PAL    groupIX  2075079       7 0.00000000 0.00000000
#> 209      pop      SMR    groupIX  2075079       7 0.31591173 0.19047619
#> 210      pop      UPD    groupIX  2075079       7 0.00000000 0.00000000
#> 211  fam.pop    A.ASP  groupXIII  2106000       8 0.33333333 0.16666667
#> 212  fam.pop    A.CLF  groupXIII  2106000       8 0.50000000 0.25000000
#> 213  fam.pop    A.OPL  groupXIII  2106000       8 0.25000000 0.12500000
#> 214  fam.pop    A.PAL  groupXIII  2106000       8 0.00000000 0.00000000
#> 215  fam.pop    A.SMR  groupXIII  2106000       8 0.33333333 0.16666667
#> 216  fam.pop    A.UPD  groupXIII  2106000       8 1.00000000 0.50000000
#> 217  fam.pop    B.ASP  groupXIII  2106000       8 0.33333333 0.16666667
#> 218  fam.pop    B.CLF  groupXIII  2106000       8 1.00000000 0.50000000
#> 219  fam.pop    B.OPL  groupXIII  2106000       8 0.26373626 0.14285714
#> 220  fam.pop    B.PAL  groupXIII  2106000       8 0.00000000 0.00000000
#> 221  fam.pop    B.SMR  groupXIII  2106000       8 0.40909091 0.25000000
#> 222  fam.pop    B.UPD  groupXIII  2106000       8 0.50000000 0.25000000
#> 223  fam.pop    C.ASP  groupXIII  2106000       8 0.33333333 0.16666667
#> 224  fam.pop    C.CLF  groupXIII  2106000       8 0.00000000 1.00000000
#> 225  fam.pop    C.OPL  groupXIII  2106000       8 0.14285714 0.07142857
#> 226  fam.pop    C.PAL  groupXIII  2106000       8 0.00000000 0.00000000
#> 227  fam.pop    C.SMR  groupXIII  2106000       8 0.46666667 0.30000000
#> 228  fam.pop    C.UPD  groupXIII  2106000       8 0.25000000 0.12500000
#> 229  fam.pop    D.ASP  groupXIII  2106000       8 0.66666667 0.50000000
#> 230  fam.pop    D.CLF  groupXIII  2106000       8 0.50000000 0.25000000
#> 231  fam.pop    D.OPL  groupXIII  2106000       8 0.48484848 0.33333333
#> 232  fam.pop    D.PAL  groupXIII  2106000       8 0.00000000 0.00000000
#> 233  fam.pop    D.SMR  groupXIII  2106000       8 0.53030303 0.58333333
#> 234  fam.pop    D.UPD  groupXIII  2106000       8 0.53333333 0.33333333
#> 235      pop      ASP  groupXIII  2106000       8 0.36796537 0.22727273
#> 236      pop      CLF  groupXIII  2106000       8 0.53030303 0.41666667
#> 237      pop      OPL  groupXIII  2106000       8 0.28368794 0.16666667
#> 238      pop      PAL  groupXIII  2106000       8 0.00000000 0.00000000
#> 239      pop      SMR  groupXIII  2106000       8 0.46666667 0.35000000
#> 240      pop      UPD  groupXIII  2106000       8 0.39473684 0.25000000
#> 241  fam.pop    A.ASP   groupXII  2136587       9 0.33333333 0.16666667
#> 242  fam.pop    A.CLF   groupXII  2136587       9 0.50000000 0.25000000
#> 243  fam.pop    A.OPL   groupXII  2136587       9 0.16666667 0.08333333
#> 244  fam.pop    A.PAL   groupXII  2136587       9 0.00000000 0.00000000
#> 245  fam.pop    A.SMR   groupXII  2136587       9 0.42857143 0.25000000
#> 246  fam.pop    A.UPD   groupXII  2136587       9 0.00000000 0.00000000
#> 247  fam.pop    B.ASP   groupXII  2136587       9 0.33333333 0.16666667
#> 248  fam.pop    B.CLF   groupXII  2136587       9 0.00000000 0.00000000
#> 249  fam.pop    B.OPL   groupXII  2136587       9 0.14285714 0.07142857
#> 250  fam.pop    B.PAL   groupXII  2136587       9 0.00000000 0.00000000
#> 251  fam.pop    B.SMR   groupXII  2136587       9 0.46666667 0.30000000
#> 252  fam.pop    B.UPD   groupXII  2136587       9 0.33333333 0.16666667
#> 253  fam.pop    C.ASP   groupXII  2136587       9 0.53333333 0.33333333
#> 254  fam.pop    C.CLF   groupXII  2136587       9 1.00000000 0.50000000
#> 255  fam.pop    C.OPL   groupXII  2136587       9 0.30303030 0.16666667
#> 256  fam.pop    C.PAL   groupXII  2136587       9 0.00000000 0.00000000
#> 257  fam.pop    C.SMR   groupXII  2136587       9 0.40909091 0.25000000
#> 258  fam.pop    C.UPD   groupXII  2136587       9 0.42857143 0.25000000
#> 259  fam.pop    D.ASP   groupXII  2136587       9 0.50000000 0.25000000
#> 260  fam.pop    D.CLF   groupXII  2136587       9 0.00000000 0.00000000
#> 261  fam.pop    D.OPL   groupXII  2136587       9 0.16666667 0.08333333
#> 262  fam.pop    D.PAL   groupXII  2136587       9 0.33333333 0.16666667
#> 263  fam.pop    D.SMR   groupXII  2136587       9 0.35555556 0.20000000
#> 264  fam.pop    D.UPD   groupXII  2136587       9 0.57142857 0.50000000
#> 265      pop      ASP   groupXII  2136587       9 0.36796537 0.22727273
#> 266      pop      CLF   groupXII  2136587       9 0.30303030 0.16666667
#> 267      pop      OPL   groupXII  2136587       9 0.18367347 0.10000000
#> 268      pop      PAL   groupXII  2136587       9 0.06250000 0.03125000
#> 269      pop      SMR   groupXII  2136587       9 0.38461538 0.25000000
#> 270      pop      UPD   groupXII  2136587       9 0.40923077 0.26923077
#> 271  fam.pop    A.ASP  groupXVII  2560839      10 0.33333333 0.16666667
#> 272  fam.pop    A.CLF  groupXVII  2560839      10 0.00000000 0.00000000
#> 273  fam.pop    A.OPL  groupXVII  2560839      10 0.00000000 0.00000000
#> 274  fam.pop    A.PAL  groupXVII  2560839      10 0.00000000 0.00000000
#> 275  fam.pop    A.SMR  groupXVII  2560839      10 0.00000000 0.00000000
#> 276  fam.pop    A.UPD  groupXVII  2560839      10 0.00000000 0.00000000
#> 277  fam.pop    B.ASP  groupXVII  2560839      10 0.00000000 0.00000000
#> 278  fam.pop    B.CLF  groupXVII  2560839      10 0.00000000 0.00000000
#> 279  fam.pop    B.OPL  groupXVII  2560839      10 0.00000000 0.00000000
#> 280  fam.pop    B.PAL  groupXVII  2560839      10 0.00000000 0.00000000
#> 281  fam.pop    B.SMR  groupXVII  2560839      10 0.00000000 0.00000000
#> 282  fam.pop    B.UPD  groupXVII  2560839      10 0.42857143 0.25000000
#> 283  fam.pop    C.ASP  groupXVII  2560839      10 0.33333333 0.16666667
#> 284  fam.pop    C.CLF  groupXVII  2560839      10 0.00000000 0.00000000
#> 285  fam.pop    C.OPL  groupXVII  2560839      10 0.16666667 0.08333333
#> 286  fam.pop    C.PAL  groupXVII  2560839      10 0.00000000 0.00000000
#> 287  fam.pop    C.SMR  groupXVII  2560839      10 0.16666667 0.08333333
#> 288  fam.pop    C.UPD  groupXVII  2560839      10 0.00000000 0.00000000
#> 289  fam.pop    D.ASP  groupXVII  2560839      10 0.00000000 0.00000000
#> 290  fam.pop    D.CLF  groupXVII  2560839      10 0.50000000 0.25000000
#> 291  fam.pop    D.OPL  groupXVII  2560839      10 0.00000000 0.00000000
#> 292  fam.pop    D.PAL  groupXVII  2560839      10 0.00000000 0.00000000
#> 293  fam.pop    D.SMR  groupXVII  2560839      10 0.16666667 0.08333333
#> 294  fam.pop    D.UPD  groupXVII  2560839      10 0.42857143 0.25000000
#> 295      pop      ASP  groupXVII  2560839      10 0.17316017 0.09090909
#> 296      pop      CLF  groupXVII  2560839      10 0.16666667 0.08333333
#> 297      pop      OPL  groupXVII  2560839      10 0.04000000 0.02000000
#> 298      pop      PAL  groupXVII  2560839      10 0.00000000 0.00000000
#> 299      pop      SMR  groupXVII  2560839      10 0.08879493 0.04545455
#> 300      pop      UPD  groupXVII  2560839      10 0.25396825 0.14285714
#> 301  fam.pop    A.ASP   groupIII  2750695      11 0.50000000 0.25000000
#> 302  fam.pop    A.CLF   groupIII  2750695      11 0.50000000 0.75000000
#> 303  fam.pop    A.OPL   groupIII  2750695      11 0.53333333 0.60000000
#> 304  fam.pop    A.PAL   groupIII  2750695      11 0.00000000 0.00000000
#> 305  fam.pop    A.SMR   groupIII  2750695      11 0.46666667 0.30000000
#> 306  fam.pop    A.UPD   groupIII  2750695      11 0.53333333 0.33333333
#> 307  fam.pop    B.ASP   groupIII  2750695      11 0.53333333 0.33333333
#> 308  fam.pop    B.CLF   groupIII  2750695      11 1.00000000 0.50000000
#> 309  fam.pop    B.OPL   groupIII  2750695      11 0.53333333 0.60000000
#> 310  fam.pop    B.PAL   groupIII  2750695      11 0.00000000 0.00000000
#> 311  fam.pop    B.SMR   groupIII  2750695      11 0.53333333 0.60000000
#> 312  fam.pop    B.UPD   groupIII  2750695      11 0.57142857 0.50000000
#> 313  fam.pop    C.ASP   groupIII  2750695      11 0.60000000 0.50000000
#> 314  fam.pop    C.CLF   groupIII  2750695      11 1.00000000 0.50000000
#> 315  fam.pop    C.OPL   groupIII  2750695      11 0.48484848 0.66666667
#> 316  fam.pop    C.PAL   groupIII  2750695      11 0.00000000 0.00000000
#> 317  fam.pop    C.SMR   groupIII  2750695      11 0.53333333 0.60000000
#> 318  fam.pop    C.UPD   groupIII  2750695      11 0.60000000 0.50000000
#> 319  fam.pop    D.ASP   groupIII  2750695      11 0.50000000 0.75000000
#> 320  fam.pop    D.CLF   groupIII  2750695      11 0.50000000 0.25000000
#> 321  fam.pop    D.OPL   groupIII  2750695      11 0.55555556 0.50000000
#> 322  fam.pop    D.PAL   groupIII  2750695      11 0.00000000 0.00000000
#> 323  fam.pop    D.SMR   groupIII  2750695      11 0.57142857 0.50000000
#> 324  fam.pop    D.UPD   groupIII  2750695      11 0.42857143 0.25000000
#> 325      pop      ASP   groupIII  2750695      11 0.52105263 0.45000000
#> 326      pop      CLF   groupIII  2750695      11 0.54545455 0.50000000
#> 327      pop      OPL   groupIII  2750695      11 0.49361208 0.59523810
#> 328      pop      PAL   groupIII  2750695      11 0.00000000 0.00000000
#> 329      pop      SMR   groupIII  2750695      11 0.51351351 0.50000000
#> 330      pop      UPD   groupIII  2750695      11 0.49470899 0.39285714
#> 331  fam.pop    A.ASP   groupXIX  3008524      12 0.53333333 0.33333333
#> 332  fam.pop    A.CLF   groupXIX  3008524      12 0.00000000 0.00000000
#> 333  fam.pop    A.OPL   groupXIX  3008524      12 0.00000000 0.00000000
#> 334  fam.pop    A.PAL   groupXIX  3008524      12 0.46666667 0.30000000
#> 335  fam.pop    A.SMR   groupXIX  3008524      12 0.35555556 0.20000000
#> 336  fam.pop    A.UPD   groupXIX  3008524      12 0.60000000 0.50000000
#> 337  fam.pop    B.ASP   groupXIX  3008524      12 0.33333333 0.16666667
#> 338  fam.pop    B.CLF   groupXIX  3008524      12 1.00000000 0.50000000
#> 339  fam.pop    B.OPL   groupXIX  3008524      12 0.42857143 0.25000000
#> 340  fam.pop    B.PAL   groupXIX  3008524      12 0.53571429 0.37500000
#> 341  fam.pop    B.SMR   groupXIX  3008524      12 0.33333333 0.16666667
#> 342  fam.pop    B.UPD   groupXIX  3008524      12 0.33333333 0.16666667
#> 343  fam.pop    C.ASP   groupXIX  3008524      12 0.53333333 0.33333333
#> 344  fam.pop    C.CLF   groupXIX  3008524      12 0.00000000 0.00000000
#> 345  fam.pop    C.OPL   groupXIX  3008524      12 0.46666667 0.30000000
#> 346  fam.pop    C.PAL   groupXIX  3008524      12 0.57142857 0.50000000
#> 347  fam.pop    C.SMR   groupXIX  3008524      12 0.35555556 0.20000000
#> 348  fam.pop    C.UPD   groupXIX  3008524      12 0.50000000 0.25000000
#> 349  fam.pop    D.ASP   groupXIX  3008524      12 0.00000000 0.00000000
#> 350  fam.pop    D.CLF   groupXIX  3008524      12 0.00000000 0.00000000
#> 351  fam.pop    D.OPL   groupXIX  3008524      12 0.35555556 0.20000000
#> 352  fam.pop    D.PAL   groupXIX  3008524      12 0.00000000 0.00000000
#> 353  fam.pop    D.SMR   groupXIX  3008524      12 0.25000000 0.12500000
#> 354  fam.pop    D.UPD   groupXIX  3008524      12 0.50000000 0.25000000
#> 355      pop      ASP   groupXIX  3008524      12 0.39473684 0.25000000
#> 356      pop      CLF   groupXIX  3008524      12 0.16666667 0.08333333
#> 357      pop      OPL   groupXIX  3008524      12 0.32222222 0.19444444
#> 358      pop      PAL   groupXIX  3008524      12 0.44354839 0.31250000
#> 359      pop      SMR   groupXIX  3008524      12 0.29946524 0.17647059
#> 360      pop      UPD   groupXIX  3008524      12 0.44210526 0.30000000
#> 361  fam.pop    A.ASP   groupXVI  3624417      13 0.50000000 0.25000000
#> 362  fam.pop    A.CLF   groupXVI  3624417      13 0.50000000 0.75000000
#> 363  fam.pop    A.OPL   groupXVI  3624417      13 0.35555556 0.20000000
#> 364  fam.pop    A.PAL   groupXVI  3624417      13 0.20000000 0.10000000
#> 365  fam.pop    A.SMR   groupXVI  3624417      13 0.48484848 0.33333333
#> 366  fam.pop    A.UPD   groupXVI  3624417      13 0.00000000 0.00000000
#> 367  fam.pop    B.ASP   groupXVI  3624417      13 0.00000000 0.00000000
#> 368  fam.pop    B.CLF   groupXVI  3624417      13 0.00000000 0.00000000
#> 369  fam.pop    B.OPL   groupXVI  3624417      13 0.53030303 0.41666667
#> 370  fam.pop    B.PAL   groupXVI  3624417      13 0.00000000 0.00000000
#> 371  fam.pop    B.SMR   groupXVI  3624417      13 0.25000000 0.12500000
#> 372  fam.pop    B.UPD   groupXVI  3624417      13 0.25000000 0.12500000
#> 373  fam.pop    C.ASP   groupXVI  3624417      13 0.33333333 0.16666667
#> 374  fam.pop    C.CLF   groupXVI  3624417      13 1.00000000 0.50000000
#> 375  fam.pop    C.OPL   groupXVI  3624417      13 0.14285714 0.07142857
#> 376  fam.pop    C.PAL   groupXVI  3624417      13 0.00000000 0.00000000
#> 377  fam.pop    C.SMR   groupXVI  3624417      13 0.42857143 0.25000000
#> 378  fam.pop    C.UPD   groupXVI  3624417      13 0.53333333 0.33333333
#> 379  fam.pop    D.ASP   groupXVI  3624417      13 0.00000000 0.00000000
#> 380  fam.pop    D.CLF   groupXVI  3624417      13 0.50000000 0.25000000
#> 381  fam.pop    D.OPL   groupXVI  3624417      13 0.48484848 0.33333333
#> 382  fam.pop    D.PAL   groupXVI  3624417      13 0.00000000 0.00000000
#> 383  fam.pop    D.SMR   groupXVI  3624417      13 0.16666667 0.08333333
#> 384  fam.pop    D.UPD   groupXVI  3624417      13 0.33333333 0.16666667
#> 385      pop      ASP   groupXVI  3624417      13 0.20915033 0.11111111
#> 386      pop      CLF   groupXVI  3624417      13 0.53030303 0.41666667
#> 387      pop      OPL   groupXVI  3624417      13 0.38297872 0.25000000
#> 388      pop      PAL   groupXVI  3624417      13 0.06666667 0.03333333
#> 389      pop      SMR   groupXVI  3624417      13 0.32820513 0.20000000
#> 390      pop      UPD   groupXVI  3624417      13 0.27076923 0.15384615
#> 391  fam.pop    A.ASP   groupVII  3667528      14 0.33333333 0.16666667
#> 392  fam.pop    A.CLF   groupVII  3667528      14 0.00000000 0.00000000
#> 393  fam.pop    A.OPL   groupVII  3667528      14 0.30303030 0.16666667
#> 394  fam.pop    A.PAL   groupVII  3667528      14 0.00000000 0.00000000
#> 395  fam.pop    A.SMR   groupVII  3667528      14 0.35555556 0.20000000
#> 396  fam.pop    A.UPD   groupVII  3667528      14 0.53333333 0.33333333
#> 397  fam.pop    B.ASP   groupVII  3667528      14 0.33333333 0.16666667
#> 398  fam.pop    B.CLF   groupVII  3667528      14 0.00000000 0.00000000
#> 399  fam.pop    B.OPL   groupVII  3667528      14 0.00000000 0.00000000
#> 400  fam.pop    B.PAL   groupVII  3667528      14 0.00000000 0.00000000
#> 401  fam.pop    B.SMR   groupVII  3667528      14 0.20000000 0.10000000
#> 402  fam.pop    B.UPD   groupVII  3667528      14 0.42857143 0.25000000
#> 403  fam.pop    C.ASP   groupVII  3667528      14 0.53333333 0.33333333
#> 404  fam.pop    C.CLF   groupVII  3667528      14 1.00000000 0.50000000
#> 405  fam.pop    C.OPL   groupVII  3667528      14 0.00000000 0.00000000
#> 406  fam.pop    C.PAL   groupVII  3667528      14 0.00000000 0.00000000
#> 407  fam.pop    C.SMR   groupVII  3667528      14 0.16666667 0.08333333
#> 408  fam.pop    C.UPD   groupVII  3667528      14 0.53571429 0.37500000
#> 409  fam.pop    D.ASP   groupVII  3667528      14 0.00000000 0.00000000
#> 410  fam.pop    D.CLF   groupVII  3667528      14 0.50000000 0.25000000
#> 411  fam.pop    D.OPL   groupVII  3667528      14 0.26373626 0.14285714
#> 412  fam.pop    D.PAL   groupVII  3667528      14 0.00000000 0.00000000
#> 413  fam.pop    D.SMR   groupVII  3667528      14 0.48484848 0.33333333
#> 414  fam.pop    D.UPD   groupVII  3667528      14 0.57142857 0.50000000
#> 415      pop      ASP   groupVII  3667528      14 0.33684211 0.20000000
#> 416      pop      CLF   groupVII  3667528      14 0.30303030 0.16666667
#> 417      pop      OPL   groupVII  3667528      14 0.15020408 0.08000000
#> 418      pop      PAL   groupVII  3667528      14 0.00000000 0.00000000
#> 419      pop      SMR   groupVII  3667528      14 0.30443975 0.18181818
#> 420      pop      UPD   groupVII  3667528      14 0.48045977 0.36666667
#> 421  fam.pop    A.ASP   groupVII  4189120      15 0.33333333 0.16666667
#> 422  fam.pop    A.CLF   groupVII  4189120      15 0.00000000 0.00000000
#> 423  fam.pop    A.OPL   groupVII  4189120      15 0.00000000 0.00000000
#> 424  fam.pop    A.PAL   groupVII  4189120      15 0.25000000 0.87500000
#> 425  fam.pop    A.SMR   groupVII  4189120      15 0.40909091 0.25000000
#> 426  fam.pop    A.UPD   groupVII  4189120      15 1.00000000 0.50000000
#> 427  fam.pop    B.ASP   groupVII  4189120      15 0.53333333 0.33333333
#> 428  fam.pop    B.CLF   groupVII  4189120      15 0.00000000 0.00000000
#> 429  fam.pop    B.OPL   groupVII  4189120      15 0.53333333 0.40000000
#> 430  fam.pop    B.PAL   groupVII  4189120      15 0.00000000 1.00000000
#> 431  fam.pop    B.SMR   groupVII  4189120      15 0.20000000 0.10000000
#> 432  fam.pop    B.UPD   groupVII  4189120      15 0.57142857 0.50000000
#> 433  fam.pop    C.ASP   groupVII  4189120      15 0.53333333 0.33333333
#> 434  fam.pop    C.CLF   groupVII  4189120      15 0.00000000 0.00000000
#> 435  fam.pop    C.OPL   groupVII  4189120      15 0.46666667 0.30000000
#> 436  fam.pop    C.PAL   groupVII  4189120      15 0.42857143 0.75000000
#> 437  fam.pop    C.SMR   groupVII  4189120      15 0.40909091 0.25000000
#> 438  fam.pop    C.UPD   groupVII  4189120      15 0.53571429 0.37500000
#> 439  fam.pop    D.ASP   groupVII  4189120      15 0.00000000 0.00000000
#> 440  fam.pop    D.CLF   groupVII  4189120      15 0.00000000 0.00000000
#> 441  fam.pop    D.OPL   groupVII  4189120      15 0.25000000 0.12500000
#> 442  fam.pop    D.PAL   groupVII  4189120      15 0.25000000 0.87500000
#> 443  fam.pop    D.SMR   groupVII  4189120      15 0.42857143 0.25000000
#> 444  fam.pop    D.UPD   groupVII  4189120      15 0.57142857 0.50000000
#> 445      pop      ASP   groupVII  4189120      15 0.36796537 0.22727273
#> 446      pop      CLF   groupVII  4189120      15 0.00000000 0.00000000
#> 447      pop      OPL   groupVII  4189120      15 0.32820513 0.20000000
#> 448      pop      PAL   groupVII  4189120      15 0.22580645 0.87500000
#> 449      pop      SMR   groupVII  4189120      15 0.34494774 0.21428571
#> 450      pop      UPD   groupVII  4189120      15 0.51692308 0.46153846
#> 451  fam.pop    A.ASP    groupXX  4538074      16 0.00000000 0.00000000
#> 452  fam.pop    A.CLF    groupXX  4538074      16 0.66666667 0.50000000
#> 453  fam.pop    A.OPL    groupXX  4538074      16 0.46666667 0.30000000
#> 454  fam.pop    A.PAL    groupXX  4538074      16 0.53333333 0.40000000
#> 455  fam.pop    A.SMR    groupXX  4538074      16 0.48484848 0.33333333
#> 456  fam.pop    A.UPD    groupXX  4538074      16 0.33333333 0.16666667
#> 457  fam.pop    B.ASP    groupXX  4538074      16 0.33333333 0.16666667
#> 458  fam.pop    B.CLF    groupXX  4538074      16 0.00000000 0.00000000
#> 459  fam.pop    B.OPL    groupXX  4538074      16 0.35555556 0.20000000
#> 460  fam.pop    B.PAL    groupXX  4538074      16 0.53333333 0.33333333
#> 461  fam.pop    B.SMR    groupXX  4538074      16 0.53571429 0.37500000
#> 462  fam.pop    B.UPD    groupXX  4538074      16 0.53333333 0.33333333
#> 463  fam.pop    C.ASP    groupXX  4538074      16 0.53333333 0.33333333
#> 464  fam.pop    C.CLF    groupXX  4538074      16        NaN         NA
#> 465  fam.pop    C.OPL    groupXX  4538074      16 0.48484848 0.33333333
#> 466  fam.pop    C.PAL    groupXX  4538074      16 0.53571429 0.62500000
#> 467  fam.pop    C.SMR    groupXX  4538074      16 0.30303030 0.16666667
#> 468  fam.pop    C.UPD    groupXX  4538074      16 0.33333333 0.16666667
#> 469  fam.pop    D.ASP    groupXX  4538074      16 0.66666667 0.50000000
#> 470  fam.pop    D.CLF    groupXX  4538074      16 0.50000000 0.25000000
#> 471  fam.pop    D.OPL    groupXX  4538074      16 0.57142857 0.50000000
#> 472  fam.pop    D.PAL    groupXX  4538074      16 0.53333333 0.33333333
#> 473  fam.pop    D.SMR    groupXX  4538074      16 0.57142857 0.50000000
#> 474  fam.pop    D.UPD    groupXX  4538074      16 0.33333333 0.16666667
#> 475      pop      ASP    groupXX  4538074      16 0.39473684 0.25000000
#> 476      pop      CLF    groupXX  4538074      16 0.46666667 0.30000000
#> 477      pop      OPL    groupXX  4538074      16 0.45000000 0.32500000
#> 478      pop      PAL    groupXX  4538074      16 0.50804598 0.43333333
#> 479      pop      SMR    groupXX  4538074      16 0.45000000 0.32500000
#> 480      pop      UPD    groupXX  4538074      16 0.34420290 0.20833333
#> 481  fam.pop    A.ASP  groupXIII  4706089      17 0.00000000 0.00000000
#> 482  fam.pop    A.CLF  groupXIII  4706089      17 0.00000000 0.00000000
#> 483  fam.pop    A.OPL  groupXIII  4706089      17 0.35555556 0.20000000
#> 484  fam.pop    A.PAL  groupXIII  4706089      17 0.42857143 0.25000000
#> 485  fam.pop    A.SMR  groupXIII  4706089      17 0.40909091 0.25000000
#> 486  fam.pop    A.UPD  groupXIII  4706089      17 0.33333333 0.16666667
#> 487  fam.pop    B.ASP  groupXIII  4706089      17 0.50000000 0.25000000
#> 488  fam.pop    B.CLF  groupXIII  4706089      17 0.00000000 0.00000000
#> 489  fam.pop    B.OPL  groupXIII  4706089      17 0.36263736 0.21428571
#> 490  fam.pop    B.PAL  groupXIII  4706089      17 0.42857143 0.25000000
#> 491  fam.pop    B.SMR  groupXIII  4706089      17 0.35555556 0.20000000
#> 492  fam.pop    B.UPD  groupXIII  4706089      17 0.42857143 0.25000000
#> 493  fam.pop    C.ASP  groupXIII  4706089      17 0.33333333 0.16666667
#> 494  fam.pop    C.CLF  groupXIII  4706089      17 0.00000000 0.00000000
#> 495  fam.pop    C.OPL  groupXIII  4706089      17 0.00000000 0.00000000
#> 496  fam.pop    C.PAL  groupXIII  4706089      17 0.25000000 0.12500000
#> 497  fam.pop    C.SMR  groupXIII  4706089      17 0.00000000 0.00000000
#> 498  fam.pop    C.UPD  groupXIII  4706089      17 0.33333333 0.16666667
#> 499  fam.pop    D.ASP  groupXIII  4706089      17 0.50000000 0.25000000
#> 500  fam.pop    D.CLF  groupXIII  4706089      17 0.00000000 0.00000000
#> 501  fam.pop    D.OPL  groupXIII  4706089      17 0.26373626 0.14285714
#> 502  fam.pop    D.PAL  groupXIII  4706089      17 0.25000000 0.12500000
#> 503  fam.pop    D.SMR  groupXIII  4706089      17 0.35555556 0.20000000
#> 504  fam.pop    D.UPD  groupXIII  4706089      17 0.00000000 0.00000000
#> 505      pop      ASP  groupXIII  4706089      17 0.26842105 0.15000000
#> 506      pop      CLF  groupXIII  4706089      17 0.00000000 0.00000000
#> 507      pop      OPL  groupXIII  4706089      17 0.24571429 0.14000000
#> 508      pop      PAL  groupXIII  4706089      17 0.31451613 0.18750000
#> 509      pop      SMR  groupXIII  4706089      17 0.28455285 0.16666667
#> 510      pop      UPD  groupXIII  4706089      17 0.27076923 0.15384615
#> 511  fam.pop    A.ASP   groupIII  4746808      18 0.60000000 0.50000000
#> 512  fam.pop    A.CLF   groupIII  4746808      18 0.66666667 0.50000000
#> 513  fam.pop    A.OPL   groupIII  4746808      18 0.40909091 0.25000000
#> 514  fam.pop    A.PAL   groupIII  4746808      18 0.53333333 0.66666667
#> 515  fam.pop    A.SMR   groupIII  4746808      18 0.40909091 0.25000000
#> 516  fam.pop    A.UPD   groupIII  4746808      18 0.33333333 0.16666667
#> 517  fam.pop    B.ASP   groupIII  4746808      18 0.50000000 0.25000000
#> 518  fam.pop    B.CLF   groupIII  4746808      18 1.00000000 0.50000000
#> 519  fam.pop    B.OPL   groupIII  4746808      18 0.36263736 0.21428571
#> 520  fam.pop    B.PAL   groupIII  4746808      18 0.53571429 0.62500000
#> 521  fam.pop    B.SMR   groupIII  4746808      18 0.35555556 0.20000000
#> 522  fam.pop    B.UPD   groupIII  4746808      18 0.53571429 0.62500000
#> 523  fam.pop    C.ASP   groupIII  4746808      18 0.33333333 0.16666667
#> 524  fam.pop    C.CLF   groupIII  4746808      18        NaN         NA
#> 525  fam.pop    C.OPL   groupIII  4746808      18 0.26373626 0.14285714
#> 526  fam.pop    C.PAL   groupIII  4746808      18 0.53333333 0.66666667
#> 527  fam.pop    C.SMR   groupIII  4746808      18 0.16666667 0.08333333
#> 528  fam.pop    C.UPD   groupIII  4746808      18 0.42857143 0.25000000
#> 529  fam.pop    D.ASP   groupIII  4746808      18 0.50000000 0.25000000
#> 530  fam.pop    D.CLF   groupIII  4746808      18 0.50000000 0.25000000
#> 531  fam.pop    D.OPL   groupIII  4746808      18 0.36263736 0.21428571
#> 532  fam.pop    D.PAL   groupIII  4746808      18 0.25000000 0.87500000
#> 533  fam.pop    D.SMR   groupIII  4746808      18 0.16666667 0.08333333
#> 534  fam.pop    D.UPD   groupIII  4746808      18 0.53571429 0.62500000
#> 535      pop      ASP   groupIII  4746808      18 0.44210526 0.30000000
#> 536      pop      CLF   groupIII  4746808      18 0.53333333 0.40000000
#> 537      pop      OPL   groupIII  4746808      18 0.33053809 0.20370370
#> 538      pop      PAL   groupIII  4746808      18 0.42328042 0.71428571
#> 539      pop      SMR   groupIII  4746808      18 0.26376812 0.15217391
#> 540      pop      UPD   groupIII  4746808      18 0.50804598 0.43333333
#> 541  fam.pop    A.ASP  groupXVII  5041995      19 0.00000000 0.00000000
#> 542  fam.pop    A.CLF  groupXVII  5041995      19 0.50000000 0.25000000
#> 543  fam.pop    A.OPL  groupXVII  5041995      19 0.20000000 0.10000000
#> 544  fam.pop    A.PAL  groupXVII  5041995      19 0.00000000 0.00000000
#> 545  fam.pop    A.SMR  groupXVII  5041995      19 0.30303030 0.16666667
#> 546  fam.pop    A.UPD  groupXVII  5041995      19 0.00000000 0.00000000
#> 547  fam.pop    B.ASP  groupXVII  5041995      19 0.00000000 0.00000000
#> 548  fam.pop    B.CLF  groupXVII  5041995      19 0.00000000 0.00000000
#> 549  fam.pop    B.OPL  groupXVII  5041995      19 0.33333333 0.16666667
#> 550  fam.pop    B.PAL  groupXVII  5041995      19 0.00000000 0.00000000
#> 551  fam.pop    B.SMR  groupXVII  5041995      19 0.00000000 0.00000000
#> 552  fam.pop    B.UPD  groupXVII  5041995      19 0.42857143 0.25000000
#> 553  fam.pop    C.ASP  groupXVII  5041995      19 0.00000000 0.00000000
#> 554  fam.pop    C.CLF  groupXVII  5041995      19 1.00000000 0.50000000
#> 555  fam.pop    C.OPL  groupXVII  5041995      19 0.16666667 0.08333333
#> 556  fam.pop    C.PAL  groupXVII  5041995      19 0.00000000 0.00000000
#> 557  fam.pop    C.SMR  groupXVII  5041995      19 0.16666667 0.08333333
#> 558  fam.pop    C.UPD  groupXVII  5041995      19 0.00000000 0.00000000
#> 559  fam.pop    D.ASP  groupXVII  5041995      19 0.00000000 0.00000000
#> 560  fam.pop    D.CLF  groupXVII  5041995      19 0.50000000 0.25000000
#> 561  fam.pop    D.OPL  groupXVII  5041995      19 0.00000000 0.00000000
#> 562  fam.pop    D.PAL  groupXVII  5041995      19 0.00000000 0.00000000
#> 563  fam.pop    D.SMR  groupXVII  5041995      19 0.00000000 0.00000000
#> 564  fam.pop    D.UPD  groupXVII  5041995      19 0.00000000 0.00000000
#> 565      pop      ASP  groupXVII  5041995      19 0.00000000 0.00000000
#> 566      pop      CLF  groupXVII  5041995      19 0.40909091 0.25000000
#> 567      pop      OPL  groupXVII  5041995      19 0.14230769 0.07500000
#> 568      pop      PAL  groupXVII  5041995      19 0.00000000 0.00000000
#> 569      pop      SMR  groupXVII  5041995      19 0.12463768 0.06521739
#> 570      pop      UPD  groupXVII  5041995      19 0.12873563 0.06666667
#> 571  fam.pop    A.ASP     groupX  5161637      20 0.00000000 0.00000000
#> 572  fam.pop    A.CLF     groupX  5161637      20 0.50000000 0.25000000
#> 573  fam.pop    A.OPL     groupX  5161637      20 0.16666667 0.08333333
#> 574  fam.pop    A.PAL     groupX  5161637      20 0.00000000 0.00000000
#> 575  fam.pop    A.SMR     groupX  5161637      20 0.16666667 0.08333333
#> 576  fam.pop    A.UPD     groupX  5161637      20 0.33333333 0.16666667
#> 577  fam.pop    B.ASP     groupX  5161637      20 0.50000000 0.25000000
#> 578  fam.pop    B.CLF     groupX  5161637      20 0.00000000 0.00000000
#> 579  fam.pop    B.OPL     groupX  5161637      20 0.00000000 0.00000000
#> 580  fam.pop    B.PAL     groupX  5161637      20 0.00000000 0.00000000
#> 581  fam.pop    B.SMR     groupX  5161637      20 0.00000000 0.00000000
#> 582  fam.pop    B.UPD     groupX  5161637      20 0.00000000 0.00000000
#> 583  fam.pop    C.ASP     groupX  5161637      20 0.00000000 0.00000000
#> 584  fam.pop    C.CLF     groupX  5161637      20 0.00000000 0.00000000
#> 585  fam.pop    C.OPL     groupX  5161637      20 0.00000000 0.00000000
#> 586  fam.pop    C.PAL     groupX  5161637      20 0.00000000 0.00000000
#> 587  fam.pop    C.SMR     groupX  5161637      20 0.00000000 0.00000000
#> 588  fam.pop    C.UPD     groupX  5161637      20 0.25000000 0.12500000
#> 589  fam.pop    D.ASP     groupX  5161637      20 0.00000000 0.00000000
#> 590  fam.pop    D.CLF     groupX  5161637      20 0.00000000 0.00000000
#> 591  fam.pop    D.OPL     groupX  5161637      20 0.00000000 0.00000000
#> 592  fam.pop    D.PAL     groupX  5161637      20 0.00000000 0.00000000
#> 593  fam.pop    D.SMR     groupX  5161637      20 0.30303030 0.16666667
#> 594  fam.pop    D.UPD     groupX  5161637      20 0.00000000 0.00000000
#> 595      pop      ASP     groupX  5161637      20 0.11111111 0.05555556
#> 596      pop      CLF     groupX  5161637      20 0.16666667 0.08333333
#> 597      pop      OPL     groupX  5161637      20 0.04166667 0.02083333
#> 598      pop      PAL     groupX  5161637      20 0.00000000 0.00000000
#> 599      pop      SMR     groupX  5161637      20 0.12463768 0.06521739
#> 600      pop      UPD     groupX  5161637      20 0.13756614 0.07142857
#> 601  fam.pop    A.ASP   groupXII  5178037      21 0.66666667 0.50000000
#> 602  fam.pop    A.CLF   groupXII  5178037      21 0.50000000 0.25000000
#> 603  fam.pop    A.OPL   groupXII  5178037      21 0.53333333 0.40000000
#> 604  fam.pop    A.PAL   groupXII  5178037      21 0.00000000 0.00000000
#> 605  fam.pop    A.SMR   groupXII  5178037      21 0.53333333 0.40000000
#> 606  fam.pop    A.UPD   groupXII  5178037      21 0.00000000 0.00000000
#> 607  fam.pop    B.ASP   groupXII  5178037      21 0.50000000 0.25000000
#> 608  fam.pop    B.CLF   groupXII  5178037      21 1.00000000 0.50000000
#> 609  fam.pop    B.OPL   groupXII  5178037      21 0.57142857 0.50000000
#> 610  fam.pop    B.PAL   groupXII  5178037      21 0.00000000 0.00000000
#> 611  fam.pop    B.SMR   groupXII  5178037      21 0.55555556 0.50000000
#> 612  fam.pop    B.UPD   groupXII  5178037      21 0.42857143 0.25000000
#> 613  fam.pop    C.ASP   groupXII  5178037      21 0.00000000 0.00000000
#> 614  fam.pop    C.CLF   groupXII  5178037      21        NaN         NA
#> 615  fam.pop    C.OPL   groupXII  5178037      21 0.54545455 0.50000000
#> 616  fam.pop    C.PAL   groupXII  5178037      21 0.00000000 0.00000000
#> 617  fam.pop    C.SMR   groupXII  5178037      21 0.42857143 0.25000000
#> 618  fam.pop    C.UPD   groupXII  5178037      21 0.60000000 0.50000000
#> 619  fam.pop    D.ASP   groupXII  5178037      21 0.00000000 0.00000000
#> 620  fam.pop    D.CLF   groupXII  5178037      21 1.00000000 0.50000000
#> 621  fam.pop    D.OPL   groupXII  5178037      21 0.53571429 0.62500000
#> 622  fam.pop    D.PAL   groupXII  5178037      21 0.00000000 0.00000000
#> 623  fam.pop    D.SMR   groupXII  5178037      21 0.30303030 0.16666667
#> 624  fam.pop    D.UPD   groupXII  5178037      21 0.53571429 0.37500000
#> 625      pop      ASP   groupXII  5178037      21 0.29411765 0.16666667
#> 626      pop      CLF   groupXII  5178037      21 0.53571429 0.37500000
#> 627      pop      OPL   groupXII  5178037      21 0.51351351 0.50000000
#> 628      pop      PAL   groupXII  5178037      21 0.00000000 0.00000000
#> 629      pop      SMR   groupXII  5178037      21 0.45000000 0.32500000
#> 630      pop      UPD   groupXII  5178037      21 0.46376812 0.33333333
#> 631  fam.pop    A.ASP    groupIX  5287711      22 0.60000000 0.50000000
#> 632  fam.pop    A.CLF    groupIX  5287711      22 0.50000000 0.25000000
#> 633  fam.pop    A.OPL    groupIX  5287711      22 0.46666667 0.30000000
#> 634  fam.pop    A.PAL    groupIX  5287711      22 0.00000000 0.00000000
#> 635  fam.pop    A.SMR    groupIX  5287711      22 0.25000000 0.12500000
#> 636  fam.pop    A.UPD    groupIX  5287711      22 0.66666667 0.50000000
#> 637  fam.pop    B.ASP    groupIX  5287711      22 0.00000000 0.00000000
#> 638  fam.pop    B.CLF    groupIX  5287711      22        NaN         NA
#> 639  fam.pop    B.OPL    groupIX  5287711      22 0.52747253 0.42857143
#> 640  fam.pop    B.PAL    groupIX  5287711      22 0.00000000 0.00000000
#> 641  fam.pop    B.SMR    groupIX  5287711      22 0.20000000 0.10000000
#> 642  fam.pop    B.UPD    groupIX  5287711      22 0.42857143 0.25000000
#> 643  fam.pop    C.ASP    groupIX  5287711      22 0.60000000 0.50000000
#> 644  fam.pop    C.CLF    groupIX  5287711      22 0.00000000 1.00000000
#> 645  fam.pop    C.OPL    groupIX  5287711      22 0.16666667 0.08333333
#> 646  fam.pop    C.PAL    groupIX  5287711      22 0.00000000 0.00000000
#> 647  fam.pop    C.SMR    groupIX  5287711      22 0.35555556 0.20000000
#> 648  fam.pop    C.UPD    groupIX  5287711      22 0.57142857 0.50000000
#> 649  fam.pop    D.ASP    groupIX  5287711      22 0.66666667 0.50000000
#> 650  fam.pop    D.CLF    groupIX  5287711      22 0.50000000 0.25000000
#> 651  fam.pop    D.OPL    groupIX  5287711      22 0.60000000 0.50000000
#> 652  fam.pop    D.PAL    groupIX  5287711      22 0.00000000 0.00000000
#> 653  fam.pop    D.SMR    groupIX  5287711      22 0.30303030 0.16666667
#> 654  fam.pop    D.UPD    groupIX  5287711      22 0.42857143 0.25000000
#> 655      pop      ASP    groupIX  5287711      22 0.52287582 0.44444444
#> 656      pop      CLF    groupIX  5287711      22 0.53333333 0.40000000
#> 657      pop      OPL    groupIX  5287711      22 0.43786295 0.30952381
#> 658      pop      PAL    groupIX  5287711      22 0.00000000 0.00000000
#> 659      pop      SMR    groupIX  5287711      22 0.26153846 0.15000000
#> 660      pop      UPD    groupIX  5287711      22 0.47619048 0.35714286
#> 661  fam.pop    A.ASP  groupVIII  5462135      23 0.00000000 0.00000000
#> 662  fam.pop    A.CLF  groupVIII  5462135      23 0.00000000 0.00000000
#> 663  fam.pop    A.OPL  groupVIII  5462135      23 0.00000000 0.00000000
#> 664  fam.pop    A.PAL  groupVIII  5462135      23 0.00000000 0.00000000
#> 665  fam.pop    A.SMR  groupVIII  5462135      23 0.00000000 0.00000000
#> 666  fam.pop    A.UPD  groupVIII  5462135      23 0.00000000 0.00000000
#> 667  fam.pop    B.ASP  groupVIII  5462135      23 0.00000000 0.00000000
#> 668  fam.pop    B.CLF  groupVIII  5462135      23 0.00000000 0.00000000
#> 669  fam.pop    B.OPL  groupVIII  5462135      23 0.16666667 0.08333333
#> 670  fam.pop    B.PAL  groupVIII  5462135      23 0.00000000 0.00000000
#> 671  fam.pop    B.SMR  groupVIII  5462135      23 0.00000000 0.00000000
#> 672  fam.pop    B.UPD  groupVIII  5462135      23 0.25000000 0.12500000
#> 673  fam.pop    C.ASP  groupVIII  5462135      23 0.00000000 0.00000000
#> 674  fam.pop    C.CLF  groupVIII  5462135      23 0.00000000 0.00000000
#> 675  fam.pop    C.OPL  groupVIII  5462135      23 0.00000000 0.00000000
#> 676  fam.pop    C.PAL  groupVIII  5462135      23 0.00000000 0.00000000
#> 677  fam.pop    C.SMR  groupVIII  5462135      23 0.53571429 0.37500000
#> 678  fam.pop    C.UPD  groupVIII  5462135      23 0.00000000 0.00000000
#> 679  fam.pop    D.ASP  groupVIII  5462135      23 0.50000000 0.25000000
#> 680  fam.pop    D.CLF  groupVIII  5462135      23 0.00000000 0.00000000
#> 681  fam.pop    D.OPL  groupVIII  5462135      23 0.16666667 0.08333333
#> 682  fam.pop    D.PAL  groupVIII  5462135      23 0.00000000 0.00000000
#> 683  fam.pop    D.SMR  groupVIII  5462135      23 0.16666667 0.08333333
#> 684  fam.pop    D.UPD  groupVIII  5462135      23 0.00000000 0.00000000
#> 685      pop      ASP  groupVIII  5462135      23 0.10000000 0.05000000
#> 686      pop      CLF  groupVIII  5462135      23 0.00000000 0.00000000
#> 687      pop      OPL  groupVIII  5462135      23 0.07836735 0.04000000
#> 688      pop      PAL  groupVIII  5462135      23 0.00000000 0.00000000
#> 689      pop      SMR  groupVIII  5462135      23 0.18461538 0.10000000
#> 690      pop      UPD  groupVIII  5462135      23 0.06666667 0.03333333
#> 691  fam.pop    A.ASP    groupIX  5960933      24 0.00000000 0.00000000
#> 692  fam.pop    A.CLF    groupIX  5960933      24 0.00000000 0.00000000
#> 693  fam.pop    A.OPL    groupIX  5960933      24 0.30303030 0.16666667
#> 694  fam.pop    A.PAL    groupIX  5960933      24 0.00000000 0.00000000
#> 695  fam.pop    A.SMR    groupIX  5960933      24 0.00000000 0.00000000
#> 696  fam.pop    A.UPD    groupIX  5960933      24 0.53333333 0.33333333
#> 697  fam.pop    B.ASP    groupIX  5960933      24 0.00000000 0.00000000
#> 698  fam.pop    B.CLF    groupIX  5960933      24 0.00000000 0.00000000
#> 699  fam.pop    B.OPL    groupIX  5960933      24 0.00000000 0.00000000
#> 700  fam.pop    B.PAL    groupIX  5960933      24 0.00000000 0.00000000
#> 701  fam.pop    B.SMR    groupIX  5960933      24 0.00000000 0.00000000
#> 702  fam.pop    B.UPD    groupIX  5960933      24 0.25000000 0.12500000
#> 703  fam.pop    C.ASP    groupIX  5960933      24 0.33333333 0.16666667
#> 704  fam.pop    C.CLF    groupIX  5960933      24 0.00000000 0.00000000
#> 705  fam.pop    C.OPL    groupIX  5960933      24 0.30303030 0.16666667
#> 706  fam.pop    C.PAL    groupIX  5960933      24 0.00000000 0.00000000
#> 707  fam.pop    C.SMR    groupIX  5960933      24 0.00000000 0.00000000
#> 708  fam.pop    C.UPD    groupIX  5960933      24 0.00000000 0.00000000
#> 709  fam.pop    D.ASP    groupIX  5960933      24 0.00000000 0.00000000
#> 710  fam.pop    D.CLF    groupIX  5960933      24 0.00000000 0.00000000
#> 711  fam.pop    D.OPL    groupIX  5960933      24 0.14285714 0.07142857
#> 712  fam.pop    D.PAL    groupIX  5960933      24 0.00000000 0.00000000
#> 713  fam.pop    D.SMR    groupIX  5960933      24 0.35555556 0.20000000
#> 714  fam.pop    D.UPD    groupIX  5960933      24 0.00000000 0.00000000
#> 715      pop      ASP    groupIX  5960933      24 0.10000000 0.05000000
#> 716      pop      CLF    groupIX  5960933      24 0.00000000 0.00000000
#> 717      pop      OPL    groupIX  5960933      24 0.17722474 0.09615385
#> 718      pop      PAL    groupIX  5960933      24 0.00000000 0.00000000
#> 719      pop      SMR    groupIX  5960933      24 0.10241821 0.05263158
#> 720      pop      UPD    groupIX  5960933      24 0.19841270 0.10714286
#> 721  fam.pop    A.ASP  groupXIII  6076868      25 0.50000000 0.25000000
#> 722  fam.pop    A.CLF  groupXIII  6076868      25 0.00000000 0.00000000
#> 723  fam.pop    A.OPL  groupXIII  6076868      25 0.30303030 0.16666667
#> 724  fam.pop    A.PAL  groupXIII  6076868      25 0.00000000 0.00000000
#> 725  fam.pop    A.SMR  groupXIII  6076868      25 0.55555556 0.50000000
#> 726  fam.pop    A.UPD  groupXIII  6076868      25 0.60000000 0.50000000
#> 727  fam.pop    B.ASP  groupXIII  6076868      25 0.60000000 0.50000000
#> 728  fam.pop    B.CLF  groupXIII  6076868      25 1.00000000 0.50000000
#> 729  fam.pop    B.OPL  groupXIII  6076868      25 0.40909091 0.25000000
#> 730  fam.pop    B.PAL  groupXIII  6076868      25 0.00000000 0.00000000
#> 731  fam.pop    B.SMR  groupXIII  6076868      25 0.46666667 0.30000000
#> 732  fam.pop    B.UPD  groupXIII  6076868      25 0.57142857 0.50000000
#> 733  fam.pop    C.ASP  groupXIII  6076868      25 0.53333333 0.33333333
#> 734  fam.pop    C.CLF  groupXIII  6076868      25        NaN         NA
#> 735  fam.pop    C.OPL  groupXIII  6076868      25 0.25000000 0.12500000
#> 736  fam.pop    C.PAL  groupXIII  6076868      25 0.00000000 0.00000000
#> 737  fam.pop    C.SMR  groupXIII  6076868      25 0.42857143 0.75000000
#> 738  fam.pop    C.UPD  groupXIII  6076868      25 0.53571429 0.37500000
#> 739  fam.pop    D.ASP  groupXIII  6076868      25 0.00000000 0.00000000
#> 740  fam.pop    D.CLF  groupXIII  6076868      25 0.00000000 0.00000000
#> 741  fam.pop    D.OPL  groupXIII  6076868      25 0.52747253 0.42857143
#> 742  fam.pop    D.PAL  groupXIII  6076868      25 0.00000000 0.00000000
#> 743  fam.pop    D.SMR  groupXIII  6076868      25 0.55555556 0.50000000
#> 744  fam.pop    D.UPD  groupXIII  6076868      25 0.53571429 0.37500000
#> 745      pop      ASP  groupXIII  6076868      25 0.44210526 0.30000000
#> 746      pop      CLF  groupXIII  6076868      25 0.20000000 0.10000000
#> 747      pop      OPL  groupXIII  6076868      25 0.39420290 0.26086957
#> 748      pop      PAL  groupXIII  6076868      25 0.00000000 0.00000000
#> 749      pop      SMR  groupXIII  6076868      25 0.51351351 0.50000000
#> 750      pop      UPD  groupXIII  6076868      25 0.50804598 0.43333333
#> 751  fam.pop    A.ASP    groupXI  6177332      26 0.66666667 0.50000000
#> 752  fam.pop    A.CLF    groupXI  6177332      26 0.00000000 0.00000000
#> 753  fam.pop    A.OPL    groupXI  6177332      26 0.00000000 0.00000000
#> 754  fam.pop    A.PAL    groupXI  6177332      26 0.53333333 0.60000000
#> 755  fam.pop    A.SMR    groupXI  6177332      26 0.40909091 0.25000000
#> 756  fam.pop    A.UPD    groupXI  6177332      26 0.50000000 0.25000000
#> 757  fam.pop    B.ASP    groupXI  6177332      26 0.33333333 0.16666667
#> 758  fam.pop    B.CLF    groupXI  6177332      26 0.00000000 0.00000000
#> 759  fam.pop    B.OPL    groupXI  6177332      26 0.53333333 0.40000000
#> 760  fam.pop    B.PAL    groupXI  6177332      26 0.53571429 0.37500000
#> 761  fam.pop    B.SMR    groupXI  6177332      26 0.48484848 0.33333333
#> 762  fam.pop    B.UPD    groupXI  6177332      26 0.53333333 0.33333333
#> 763  fam.pop    C.ASP    groupXI  6177332      26 0.66666667 0.50000000
#> 764  fam.pop    C.CLF    groupXI  6177332      26 0.00000000 0.00000000
#> 765  fam.pop    C.OPL    groupXI  6177332      26 0.26373626 0.14285714
#> 766  fam.pop    C.PAL    groupXI  6177332      26 0.57142857 0.50000000
#> 767  fam.pop    C.SMR    groupXI  6177332      26 0.46666667 0.30000000
#> 768  fam.pop    C.UPD    groupXI  6177332      26 0.42857143 0.25000000
#> 769  fam.pop    D.ASP    groupXI  6177332      26 0.50000000 0.25000000
#> 770  fam.pop    D.CLF    groupXI  6177332      26 0.00000000 0.00000000
#> 771  fam.pop    D.OPL    groupXI  6177332      26 0.00000000 0.00000000
#> 772  fam.pop    D.PAL    groupXI  6177332      26 0.50000000 0.75000000
#> 773  fam.pop    D.SMR    groupXI  6177332      26 0.46666667 0.30000000
#> 774  fam.pop    D.UPD    groupXI  6177332      26 0.53333333 0.33333333
#> 775      pop      ASP    groupXI  6177332      26 0.47058824 0.33333333
#> 776      pop      CLF    groupXI  6177332      26 0.00000000 0.00000000
#> 777      pop      OPL    groupXI  6177332      26 0.25087108 0.14285714
#> 778      pop      PAL    groupXI  6177332      26 0.51494253 0.53333333
#> 779      pop      SMR    groupXI  6177332      26 0.42600423 0.29545455
#> 780      pop      UPD    groupXI  6177332      26 0.43115942 0.29166667
#> 781  fam.pop    A.ASP    groupIV  6291287      27 0.00000000 0.00000000
#> 782  fam.pop    A.CLF    groupIV  6291287      27 0.50000000 0.25000000
#> 783  fam.pop    A.OPL    groupIV  6291287      27 0.40909091 0.25000000
#> 784  fam.pop    A.PAL    groupIV  6291287      27 0.00000000 0.00000000
#> 785  fam.pop    A.SMR    groupIV  6291287      27 0.20000000 0.10000000
#> 786  fam.pop    A.UPD    groupIV  6291287      27 0.00000000 0.00000000
#> 787  fam.pop    B.ASP    groupIV  6291287      27 0.33333333 0.16666667
#> 788  fam.pop    B.CLF    groupIV  6291287      27 0.00000000 0.00000000
#> 789  fam.pop    B.OPL    groupIV  6291287      27 0.36263736 0.21428571
#> 790  fam.pop    B.PAL    groupIV  6291287      27 0.00000000 0.00000000
#> 791  fam.pop    B.SMR    groupIV  6291287      27 0.16666667 0.08333333
#> 792  fam.pop    B.UPD    groupIV  6291287      27 0.00000000 0.00000000
#> 793  fam.pop    C.ASP    groupIV  6291287      27 0.00000000 0.00000000
#> 794  fam.pop    C.CLF    groupIV  6291287      27 0.00000000 0.00000000
#> 795  fam.pop    C.OPL    groupIV  6291287      27 0.35555556 0.20000000
#> 796  fam.pop    C.PAL    groupIV  6291287      27 0.00000000 0.00000000
#> 797  fam.pop    C.SMR    groupIV  6291287      27 0.53333333 0.33333333
#> 798  fam.pop    C.UPD    groupIV  6291287      27 0.25000000 0.12500000
#> 799  fam.pop    D.ASP    groupIV  6291287      27 0.00000000 0.00000000
#> 800  fam.pop    D.CLF    groupIV  6291287      27 0.00000000 0.00000000
#> 801  fam.pop    D.OPL    groupIV  6291287      27 0.30303030 0.16666667
#> 802  fam.pop    D.PAL    groupIV  6291287      27 0.00000000 0.00000000
#> 803  fam.pop    D.SMR    groupIV  6291287      27 0.30303030 0.16666667
#> 804  fam.pop    D.UPD    groupIV  6291287      27 0.25000000 0.12500000
#> 805      pop      ASP    groupIV  6291287      27 0.10000000 0.05000000
#> 806      pop      CLF    groupIV  6291287      27 0.16666667 0.08333333
#> 807      pop      OPL    groupIV  6291287      27 0.33687943 0.20833333
#> 808      pop      PAL    groupIV  6291287      27 0.00000000 0.00000000
#> 809      pop      SMR    groupIV  6291287      27 0.26153846 0.15000000
#> 810      pop      UPD    groupIV  6291287      27 0.13756614 0.07142857
#> 811  fam.pop    A.ASP    groupIX  6515796      28 0.33333333 0.16666667
#> 812  fam.pop    A.CLF    groupIX  6515796      28 0.00000000 0.00000000
#> 813  fam.pop    A.OPL    groupIX  6515796      28 0.53333333 0.33333333
#> 814  fam.pop    A.PAL    groupIX  6515796      28 0.25000000 0.12500000
#> 815  fam.pop    A.SMR    groupIX  6515796      28 0.53571429 0.62500000
#> 816  fam.pop    A.UPD    groupIX  6515796      28 0.53333333 0.33333333
#> 817  fam.pop    B.ASP    groupIX  6515796      28 0.50000000 0.25000000
#> 818  fam.pop    B.CLF    groupIX  6515796      28 1.00000000 0.50000000
#> 819  fam.pop    B.OPL    groupIX  6515796      28 0.25000000 0.87500000
#> 820  fam.pop    B.PAL    groupIX  6515796      28 0.00000000 0.00000000
#> 821  fam.pop    B.SMR    groupIX  6515796      28 0.53030303 0.41666667
#> 822  fam.pop    B.UPD    groupIX  6515796      28 0.33333333 0.16666667
#> 823  fam.pop    C.ASP    groupIX  6515796      28 0.00000000 0.00000000
#> 824  fam.pop    C.CLF    groupIX  6515796      28 0.00000000 0.00000000
#> 825  fam.pop    C.OPL    groupIX  6515796      28 0.48484848 0.66666667
#> 826  fam.pop    C.PAL    groupIX  6515796      28 0.25000000 0.12500000
#> 827  fam.pop    C.SMR    groupIX  6515796      28 0.53030303 0.41666667
#> 828  fam.pop    C.UPD    groupIX  6515796      28 0.42857143 0.25000000
#> 829  fam.pop    D.ASP    groupIX  6515796      28 0.00000000 0.00000000
#> 830  fam.pop    D.CLF    groupIX  6515796      28 0.50000000 0.25000000
#> 831  fam.pop    D.OPL    groupIX  6515796      28 0.48484848 0.66666667
#> 832  fam.pop    D.PAL    groupIX  6515796      28 0.25000000 0.12500000
#> 833  fam.pop    D.SMR    groupIX  6515796      28 0.42857143 0.75000000
#> 834  fam.pop    D.UPD    groupIX  6515796      28 0.53571429 0.37500000
#> 835      pop      ASP    groupIX  6515796      28 0.23333333 0.12500000
#> 836      pop      CLF    groupIX  6515796      28 0.30303030 0.16666667
#> 837      pop      OPL    groupIX  6515796      28 0.46230441 0.65789474
#> 838      pop      PAL    groupIX  6515796      28 0.18620690 0.10000000
#> 839      pop      SMR    groupIX  6515796      28 0.51153846 0.52500000
#> 840      pop      UPD    groupIX  6515796      28 0.42328042 0.28571429
#> 841  fam.pop    A.ASP   groupXII  6718079      29 0.00000000 0.00000000
#> 842  fam.pop    A.CLF   groupXII  6718079      29 0.00000000 0.00000000
#> 843  fam.pop    A.OPL   groupXII  6718079      29 0.00000000 0.00000000
#> 844  fam.pop    A.PAL   groupXII  6718079      29 0.25000000 0.12500000
#> 845  fam.pop    A.SMR   groupXII  6718079      29 0.00000000 0.00000000
#> 846  fam.pop    A.UPD   groupXII  6718079      29 0.33333333 0.16666667
#> 847  fam.pop    B.ASP   groupXII  6718079      29 0.33333333 0.16666667
#> 848  fam.pop    B.CLF   groupXII  6718079      29        NaN         NA
#> 849  fam.pop    B.OPL   groupXII  6718079      29 0.30303030 0.16666667
#> 850  fam.pop    B.PAL   groupXII  6718079      29 0.33333333 0.16666667
#> 851  fam.pop    B.SMR   groupXII  6718079      29 0.53333333 0.40000000
#> 852  fam.pop    B.UPD   groupXII  6718079      29 0.00000000 0.00000000
#> 853  fam.pop    C.ASP   groupXII  6718079      29 0.00000000 0.00000000
#> 854  fam.pop    C.CLF   groupXII  6718079      29 0.00000000 0.00000000
#> 855  fam.pop    C.OPL   groupXII  6718079      29 0.53030303 0.41666667
#> 856  fam.pop    C.PAL   groupXII  6718079      29 0.25000000 0.12500000
#> 857  fam.pop    C.SMR   groupXII  6718079      29 0.53333333 0.33333333
#> 858  fam.pop    C.UPD   groupXII  6718079      29 0.33333333 0.16666667
#> 859  fam.pop    D.ASP   groupXII  6718079      29        NaN         NA
#> 860  fam.pop    D.CLF   groupXII  6718079      29 0.50000000 0.25000000
#> 861  fam.pop    D.OPL   groupXII  6718079      29 0.53333333 0.40000000
#> 862  fam.pop    D.PAL   groupXII  6718079      29 0.42857143 0.25000000
#> 863  fam.pop    D.SMR   groupXII  6718079      29 0.33333333 0.16666667
#> 864  fam.pop    D.UPD   groupXII  6718079      29 0.25000000 0.12500000
#> 865      pop      ASP   groupXII  6718079      29 0.11111111 0.05555556
#> 866      pop      CLF   groupXII  6718079      29 0.25000000 0.12500000
#> 867      pop      OPL   groupXII  6718079      29 0.39605110 0.26190476
#> 868      pop      PAL   groupXII  6718079      29 0.28735632 0.16666667
#> 869      pop      SMR   groupXII  6718079      29 0.37011494 0.23333333
#> 870      pop      UPD   groupXII  6718079      29 0.21230769 0.11538462
#> 871  fam.pop    A.ASP   groupVII  6967477      30 0.00000000 0.00000000
#> 872  fam.pop    A.CLF   groupVII  6967477      30 0.00000000 0.00000000
#> 873  fam.pop    A.OPL   groupVII  6967477      30 0.25000000 0.12500000
#> 874  fam.pop    A.PAL   groupVII  6967477      30 0.00000000 0.00000000
#> 875  fam.pop    A.SMR   groupVII  6967477      30 0.25000000 0.12500000
#> 876  fam.pop    A.UPD   groupVII  6967477      30 0.50000000 0.25000000
#> 877  fam.pop    B.ASP   groupVII  6967477      30 0.00000000 0.00000000
#> 878  fam.pop    B.CLF   groupVII  6967477      30 1.00000000 0.50000000
#> 879  fam.pop    B.OPL   groupVII  6967477      30 0.00000000 0.00000000
#> 880  fam.pop    B.PAL   groupVII  6967477      30 0.00000000 0.00000000
#> 881  fam.pop    B.SMR   groupVII  6967477      30 0.30303030 0.16666667
#> 882  fam.pop    B.UPD   groupVII  6967477      30 0.33333333 0.16666667
#> 883  fam.pop    C.ASP   groupVII  6967477      30 0.53333333 0.33333333
#> 884  fam.pop    C.CLF   groupVII  6967477      30 0.00000000 0.00000000
#> 885  fam.pop    C.OPL   groupVII  6967477      30 0.30303030 0.16666667
#> 886  fam.pop    C.PAL   groupVII  6967477      30 0.00000000 0.00000000
#> 887  fam.pop    C.SMR   groupVII  6967477      30 0.00000000 0.00000000
#> 888  fam.pop    C.UPD   groupVII  6967477      30 0.00000000 0.00000000
#> 889  fam.pop    D.ASP   groupVII  6967477      30 0.00000000 0.00000000
#> 890  fam.pop    D.CLF   groupVII  6967477      30 0.50000000 0.25000000
#> 891  fam.pop    D.OPL   groupVII  6967477      30 0.00000000 0.00000000
#> 892  fam.pop    D.PAL   groupVII  6967477      30 0.00000000 0.00000000
#> 893  fam.pop    D.SMR   groupVII  6967477      30 0.30303030 0.16666667
#> 894  fam.pop    D.UPD   groupVII  6967477      30 0.53333333 0.33333333
#> 895      pop      ASP   groupVII  6967477      30 0.17316017 0.09090909
#> 896      pop      CLF   groupVII  6967477      30 0.35555556 0.20000000
#> 897      pop      OPL   groupVII  6967477      30 0.13588850 0.07142857
#> 898      pop      PAL   groupVII  6967477      30 0.00000000 0.00000000
#> 899      pop      SMR   groupVII  6967477      30 0.20613108 0.11363636
#> 900      pop      UPD   groupVII  6967477      30 0.28985507 0.16666667
#> 901  fam.pop    A.ASP   groupIII  7194069      31 0.53333333 0.33333333
#> 902  fam.pop    A.CLF   groupIII  7194069      31 0.50000000 0.25000000
#> 903  fam.pop    A.OPL   groupIII  7194069      31 0.35555556 0.20000000
#> 904  fam.pop    A.PAL   groupIII  7194069      31 0.00000000 0.00000000
#> 905  fam.pop    A.SMR   groupIII  7194069      31 0.48484848 0.33333333
#> 906  fam.pop    A.UPD   groupIII  7194069      31 0.53333333 0.33333333
#> 907  fam.pop    B.ASP   groupIII  7194069      31 0.53333333 0.33333333
#> 908  fam.pop    B.CLF   groupIII  7194069      31 0.00000000 0.00000000
#> 909  fam.pop    B.OPL   groupIII  7194069      31 0.43956044 0.28571429
#> 910  fam.pop    B.PAL   groupIII  7194069      31 0.00000000 0.00000000
#> 911  fam.pop    B.SMR   groupIII  7194069      31 0.50000000 0.75000000
#> 912  fam.pop    B.UPD   groupIII  7194069      31 0.53571429 0.62500000
#> 913  fam.pop    C.ASP   groupIII  7194069      31 0.53333333 0.33333333
#> 914  fam.pop    C.CLF   groupIII  7194069      31 0.00000000 0.00000000
#> 915  fam.pop    C.OPL   groupIII  7194069      31 0.48484848 0.33333333
#> 916  fam.pop    C.PAL   groupIII  7194069      31 0.00000000 0.00000000
#> 917  fam.pop    C.SMR   groupIII  7194069      31 0.35555556 0.20000000
#> 918  fam.pop    C.UPD   groupIII  7194069      31 0.25000000 0.12500000
#> 919  fam.pop    D.ASP   groupIII  7194069      31 0.00000000 0.00000000
#> 920  fam.pop    D.CLF   groupIII  7194069      31 0.00000000 0.00000000
#> 921  fam.pop    D.OPL   groupIII  7194069      31 0.52747253 0.42857143
#> 922  fam.pop    D.PAL   groupIII  7194069      31 0.00000000 0.00000000
#> 923  fam.pop    D.SMR   groupIII  7194069      31 0.53030303 0.41666667
#> 924  fam.pop    D.UPD   groupIII  7194069      31 0.53571429 0.62500000
#> 925      pop      ASP   groupIII  7194069      31 0.41558442 0.27272727
#> 926      pop      CLF   groupIII  7194069      31 0.16666667 0.08333333
#> 927      pop      OPL   groupIII  7194069      31 0.44408163 0.32000000
#> 928      pop      PAL   groupIII  7194069      31 0.00000000 0.00000000
#> 929      pop      SMR   groupIII  7194069      31 0.47795164 0.36842105
#> 930      pop      UPD   groupIII  7194069      31 0.50804598 0.43333333
#> 931  fam.pop    A.ASP   groupXVI  7199722      32 0.33333333 0.16666667
#> 932  fam.pop    A.CLF   groupXVI  7199722      32 0.50000000 0.25000000
#> 933  fam.pop    A.OPL   groupXVI  7199722      32 0.20000000 0.10000000
#> 934  fam.pop    A.PAL   groupXVI  7199722      32 0.46666667 0.30000000
#> 935  fam.pop    A.SMR   groupXVI  7199722      32 0.25000000 0.12500000
#> 936  fam.pop    A.UPD   groupXVI  7199722      32 0.66666667 0.50000000
#> 937  fam.pop    B.ASP   groupXVI  7199722      32 0.33333333 0.16666667
#> 938  fam.pop    B.CLF   groupXVI  7199722      32        NaN         NA
#> 939  fam.pop    B.OPL   groupXVI  7199722      32 0.20000000 0.10000000
#> 940  fam.pop    B.PAL   groupXVI  7199722      32 0.60000000 0.50000000
#> 941  fam.pop    B.SMR   groupXVI  7199722      32 0.30303030 0.16666667
#> 942  fam.pop    B.UPD   groupXVI  7199722      32 0.53571429 0.37500000
#> 943  fam.pop    C.ASP   groupXVI  7199722      32 0.00000000 0.00000000
#> 944  fam.pop    C.CLF   groupXVI  7199722      32 0.00000000 0.00000000
#> 945  fam.pop    C.OPL   groupXVI  7199722      32 0.16666667 0.08333333
#> 946  fam.pop    C.PAL   groupXVI  7199722      32 0.53571429 0.37500000
#> 947  fam.pop    C.SMR   groupXVI  7199722      32 0.35555556 0.20000000
#> 948  fam.pop    C.UPD   groupXVI  7199722      32 0.33333333 0.16666667
#> 949  fam.pop    D.ASP   groupXVI  7199722      32 1.00000000 0.50000000
#> 950  fam.pop    D.CLF   groupXVI  7199722      32 0.50000000 0.25000000
#> 951  fam.pop    D.OPL   groupXVI  7199722      32 0.36263736 0.21428571
#> 952  fam.pop    D.PAL   groupXVI  7199722      32 0.60000000 0.50000000
#> 953  fam.pop    D.SMR   groupXVI  7199722      32 0.35555556 0.20000000
#> 954  fam.pop    D.UPD   groupXVI  7199722      32 0.25000000 0.12500000
#> 955      pop      ASP   groupXVI  7199722      32 0.26842105 0.15000000
#> 956      pop      CLF   groupXVI  7199722      32 0.35555556 0.20000000
#> 957      pop      OPL   groupXVI  7199722      32 0.23188406 0.13043478
#> 958      pop      PAL   groupXVI  7199722      32 0.49655172 0.40000000
#> 959      pop      SMR   groupXVI  7199722      32 0.29615385 0.17500000
#> 960      pop      UPD   groupXVI  7199722      32 0.40923077 0.26923077
#> 961  fam.pop    A.ASP  groupVIII  7224065      33 0.33333333 0.16666667
#> 962  fam.pop    A.CLF  groupVIII  7224065      33 0.00000000 0.00000000
#> 963  fam.pop    A.OPL  groupVIII  7224065      33 0.55555556 0.50000000
#> 964  fam.pop    A.PAL  groupVIII  7224065      33 0.25000000 0.12500000
#> 965  fam.pop    A.SMR  groupVIII  7224065      33 0.30303030 0.16666667
#> 966  fam.pop    A.UPD  groupVIII  7224065      33 0.66666667 0.50000000
#> 967  fam.pop    B.ASP  groupVIII  7224065      33 0.60000000 0.50000000
#> 968  fam.pop    B.CLF  groupVIII  7224065      33 1.00000000 0.50000000
#> 969  fam.pop    B.OPL  groupVIII  7224065      33 0.53333333 0.40000000
#> 970  fam.pop    B.PAL  groupVIII  7224065      33 0.42857143 0.25000000
#> 971  fam.pop    B.SMR  groupVIII  7224065      33 0.54545455 0.50000000
#> 972  fam.pop    B.UPD  groupVIII  7224065      33 0.42857143 0.25000000
#> 973  fam.pop    C.ASP  groupVIII  7224065      33 0.53333333 0.33333333
#> 974  fam.pop    C.CLF  groupVIII  7224065      33        NaN         NA
#> 975  fam.pop    C.OPL  groupVIII  7224065      33 0.25000000 0.12500000
#> 976  fam.pop    C.PAL  groupVIII  7224065      33 0.25000000 0.12500000
#> 977  fam.pop    C.SMR  groupVIII  7224065      33 0.53571429 0.37500000
#> 978  fam.pop    C.UPD  groupVIII  7224065      33 0.60000000 0.50000000
#> 979  fam.pop    D.ASP  groupVIII  7224065      33 0.00000000 0.00000000
#> 980  fam.pop    D.CLF  groupVIII  7224065      33 0.66666667 0.50000000
#> 981  fam.pop    D.OPL  groupVIII  7224065      33 0.53333333 0.40000000
#> 982  fam.pop    D.PAL  groupVIII  7224065      33 0.60000000 0.50000000
#> 983  fam.pop    D.SMR  groupVIII  7224065      33 0.48484848 0.33333333
#> 984  fam.pop    D.UPD  groupVIII  7224065      33 0.42857143 0.25000000
#> 985      pop      ASP  groupVIII  7224065      33 0.41558442 0.27272727
#> 986      pop      CLF  groupVIII  7224065      33 0.53571429 0.37500000
#> 987      pop      OPL  groupVIII  7224065      33 0.47795164 0.36842105
#> 988      pop      PAL  groupVIII  7224065      33 0.37011494 0.23333333
#> 989      pop      SMR  groupVIII  7224065      33 0.45983087 0.34090909
#> 990      pop      UPD  groupVIII  7224065      33 0.47076923 0.34615385
#> 991  fam.pop    A.ASP groupXVIII  7229426      34 0.53333333 0.33333333
#> 992  fam.pop    A.CLF groupXVIII  7229426      34 0.50000000 0.25000000
#> 993  fam.pop    A.OPL groupXVIII  7229426      34 0.46666667 0.30000000
#> 994  fam.pop    A.PAL groupXVIII  7229426      34 0.00000000 0.00000000
#> 995  fam.pop    A.SMR groupXVIII  7229426      34 0.54545455 0.50000000
#> 996  fam.pop    A.UPD groupXVIII  7229426      34 0.33333333 0.16666667
#> 997  fam.pop    B.ASP groupXVIII  7229426      34 0.66666667 0.50000000
#> 998  fam.pop    B.CLF groupXVIII  7229426      34 0.00000000 1.00000000
#> 999  fam.pop    B.OPL groupXVIII  7229426      34 0.46666667 0.30000000
#> 1000 fam.pop    B.PAL groupXVIII  7229426      34 0.00000000 0.00000000
#> 1001 fam.pop    B.SMR groupXVIII  7229426      34 0.30303030 0.16666667
#> 1002 fam.pop    B.UPD groupXVIII  7229426      34 0.33333333 0.16666667
#> 1003 fam.pop    C.ASP groupXVIII  7229426      34 0.50000000 0.25000000
#> 1004 fam.pop    C.CLF groupXVIII  7229426      34 1.00000000 0.50000000
#> 1005 fam.pop    C.OPL groupXVIII  7229426      34 0.60000000 0.50000000
#> 1006 fam.pop    C.PAL groupXVIII  7229426      34 0.00000000 0.00000000
#> 1007 fam.pop    C.SMR groupXVIII  7229426      34 0.25000000 0.12500000
#> 1008 fam.pop    C.UPD groupXVIII  7229426      34 0.00000000 0.00000000
#> 1009 fam.pop    D.ASP groupXVIII  7229426      34 0.66666667 0.50000000
#> 1010 fam.pop    D.CLF groupXVIII  7229426      34 0.66666667 0.50000000
#> 1011 fam.pop    D.OPL groupXVIII  7229426      34 0.55555556 0.50000000
#> 1012 fam.pop    D.PAL groupXVIII  7229426      34 0.00000000 0.00000000
#> 1013 fam.pop    D.SMR groupXVIII  7229426      34 0.48484848 0.33333333
#> 1014 fam.pop    D.UPD groupXVIII  7229426      34 0.50000000 0.25000000
#> 1015     pop      ASP groupXVIII  7229426      34 0.50326797 0.38888889
#> 1016     pop      CLF groupXVIII  7229426      34 0.54545455 0.50000000
#> 1017     pop      OPL groupXVIII  7229426      34 0.48888889 0.38888889
#> 1018     pop      PAL groupXVIII  7229426      34 0.00000000 0.00000000
#> 1019     pop      SMR groupXVIII  7229426      34 0.42600423 0.29545455
#> 1020     pop      UPD groupXVIII  7229426      34 0.24675325 0.13636364
#> 1021 fam.pop    A.ASP groupXVIII  8087022      35 0.00000000 0.00000000
#> 1022 fam.pop    A.CLF groupXVIII  8087022      35 0.00000000 0.00000000
#> 1023 fam.pop    A.OPL groupXVIII  8087022      35 0.35555556 0.20000000
#> 1024 fam.pop    A.PAL groupXVIII  8087022      35 0.00000000 0.00000000
#> 1025 fam.pop    A.SMR groupXVIII  8087022      35 0.00000000 0.00000000
#> 1026 fam.pop    A.UPD groupXVIII  8087022      35 0.00000000 0.00000000
#> 1027 fam.pop    B.ASP groupXVIII  8087022      35 0.00000000 0.00000000
#> 1028 fam.pop    B.CLF groupXVIII  8087022      35 0.00000000 0.00000000
#> 1029 fam.pop    B.OPL groupXVIII  8087022      35 0.16666667 0.08333333
#> 1030 fam.pop    B.PAL groupXVIII  8087022      35 0.00000000 0.00000000
#> 1031 fam.pop    B.SMR groupXVIII  8087022      35 0.25000000 0.12500000
#> 1032 fam.pop    B.UPD groupXVIII  8087022      35 0.25000000 0.12500000
#> 1033 fam.pop    C.ASP groupXVIII  8087022      35 0.33333333 0.16666667
#> 1034 fam.pop    C.CLF groupXVIII  8087022      35 0.00000000 0.00000000
#> 1035 fam.pop    C.OPL groupXVIII  8087022      35 0.00000000 0.00000000
#> 1036 fam.pop    C.PAL groupXVIII  8087022      35 0.00000000 0.00000000
#> 1037 fam.pop    C.SMR groupXVIII  8087022      35 0.00000000 0.00000000
#> 1038 fam.pop    C.UPD groupXVIII  8087022      35 0.00000000 0.00000000
#> 1039 fam.pop    D.ASP groupXVIII  8087022      35 0.50000000 0.25000000
#> 1040 fam.pop    D.CLF groupXVIII  8087022      35 0.00000000 0.00000000
#> 1041 fam.pop    D.OPL groupXVIII  8087022      35 0.00000000 0.00000000
#> 1042 fam.pop    D.PAL groupXVIII  8087022      35 0.00000000 0.00000000
#> 1043 fam.pop    D.SMR groupXVIII  8087022      35 0.30303030 0.16666667
#> 1044 fam.pop    D.UPD groupXVIII  8087022      35 0.25000000 0.12500000
#> 1045     pop      ASP groupXVIII  8087022      35 0.18947368 0.10000000
#> 1046     pop      CLF groupXVIII  8087022      35 0.00000000 0.00000000
#> 1047     pop      OPL groupXVIII  8087022      35 0.13002114 0.06818182
#> 1048     pop      PAL groupXVIII  8087022      35 0.00000000 0.00000000
#> 1049     pop      SMR groupXVIII  8087022      35 0.14230769 0.07500000
#> 1050     pop      UPD groupXVIII  8087022      35 0.13756614 0.07142857
#> 1051 fam.pop    A.ASP    groupIX  8373442      36 0.00000000 0.00000000
#> 1052 fam.pop    A.CLF    groupIX  8373442      36 0.00000000 0.00000000
#> 1053 fam.pop    A.OPL    groupIX  8373442      36 0.30303030 0.16666667
#> 1054 fam.pop    A.PAL    groupIX  8373442      36 0.00000000 0.00000000
#> 1055 fam.pop    A.SMR    groupIX  8373442      36 0.48484848 0.33333333
#> 1056 fam.pop    A.UPD    groupIX  8373442      36 0.66666667 0.50000000
#> 1057 fam.pop    B.ASP    groupIX  8373442      36 0.53333333 0.33333333
#> 1058 fam.pop    B.CLF    groupIX  8373442      36 0.00000000 0.00000000
#> 1059 fam.pop    B.OPL    groupIX  8373442      36 0.35555556 0.20000000
#> 1060 fam.pop    B.PAL    groupIX  8373442      36 0.00000000 0.00000000
#> 1061 fam.pop    B.SMR    groupIX  8373442      36 0.30303030 0.16666667
#> 1062 fam.pop    B.UPD    groupIX  8373442      36 0.00000000 0.00000000
#> 1063 fam.pop    C.ASP    groupIX  8373442      36 0.00000000 0.00000000
#> 1064 fam.pop    C.CLF    groupIX  8373442      36 1.00000000 0.50000000
#> 1065 fam.pop    C.OPL    groupIX  8373442      36 0.30303030 0.16666667
#> 1066 fam.pop    C.PAL    groupIX  8373442      36 0.00000000 0.00000000
#> 1067 fam.pop    C.SMR    groupIX  8373442      36 0.00000000 0.00000000
#> 1068 fam.pop    C.UPD    groupIX  8373442      36 0.00000000 0.00000000
#> 1069 fam.pop    D.ASP    groupIX  8373442      36 0.50000000 0.25000000
#> 1070 fam.pop    D.CLF    groupIX  8373442      36 0.50000000 0.25000000
#> 1071 fam.pop    D.OPL    groupIX  8373442      36 0.36263736 0.21428571
#> 1072 fam.pop    D.PAL    groupIX  8373442      36 0.00000000 0.00000000
#> 1073 fam.pop    D.SMR    groupIX  8373442      36 0.35555556 0.20000000
#> 1074 fam.pop    D.UPD    groupIX  8373442      36 0.00000000 0.00000000
#> 1075     pop      ASP    groupIX  8373442      36 0.26842105 0.15000000
#> 1076     pop      CLF    groupIX  8373442      36 0.30303030 0.16666667
#> 1077     pop      OPL    groupIX  8373442      36 0.31117021 0.18750000
#> 1078     pop      PAL    groupIX  8373442      36 0.00000000 0.00000000
#> 1079     pop      SMR    groupIX  8373442      36 0.30443975 0.18181818
#> 1080     pop      UPD    groupIX  8373442      36 0.14769231 0.07692308
#> 1081 fam.pop    A.ASP    groupXV  8815538      37 0.60000000 0.50000000
#> 1082 fam.pop    A.CLF    groupXV  8815538      37 0.50000000 0.25000000
#> 1083 fam.pop    A.OPL    groupXV  8815538      37 0.53030303 0.41666667
#> 1084 fam.pop    A.PAL    groupXV  8815538      37 0.00000000 0.00000000
#> 1085 fam.pop    A.SMR    groupXV  8815538      37 0.46666667 0.30000000
#> 1086 fam.pop    A.UPD    groupXV  8815538      37 0.50000000 0.25000000
#> 1087 fam.pop    B.ASP    groupXV  8815538      37 0.00000000 0.00000000
#> 1088 fam.pop    B.CLF    groupXV  8815538      37 1.00000000 0.50000000
#> 1089 fam.pop    B.OPL    groupXV  8815538      37 0.43956044 0.28571429
#> 1090 fam.pop    B.PAL    groupXV  8815538      37 0.00000000 0.00000000
#> 1091 fam.pop    B.SMR    groupXV  8815538      37 0.54545455 0.50000000
#> 1092 fam.pop    B.UPD    groupXV  8815538      37 0.53571429 0.37500000
#> 1093 fam.pop    C.ASP    groupXV  8815538      37 0.50000000 0.25000000
#> 1094 fam.pop    C.CLF    groupXV  8815538      37        NaN         NA
#> 1095 fam.pop    C.OPL    groupXV  8815538      37 0.46666667 0.30000000
#> 1096 fam.pop    C.PAL    groupXV  8815538      37 0.00000000 0.00000000
#> 1097 fam.pop    C.SMR    groupXV  8815538      37 0.30303030 0.16666667
#> 1098 fam.pop    C.UPD    groupXV  8815538      37 0.53333333 0.33333333
#> 1099 fam.pop    D.ASP    groupXV  8815538      37 0.00000000 0.00000000
#> 1100 fam.pop    D.CLF    groupXV  8815538      37 0.00000000 0.00000000
#> 1101 fam.pop    D.OPL    groupXV  8815538      37 0.16666667 0.08333333
#> 1102 fam.pop    D.PAL    groupXV  8815538      37 0.25000000 0.12500000
#> 1103 fam.pop    D.SMR    groupXV  8815538      37 0.55555556 0.50000000
#> 1104 fam.pop    D.UPD    groupXV  8815538      37 0.25000000 0.12500000
#> 1105     pop      ASP    groupXV  8815538      37 0.33684211 0.20000000
#> 1106     pop      CLF    groupXV  8815538      37 0.35555556 0.20000000
#> 1107     pop      OPL    groupXV  8815538      37 0.40336879 0.27083333
#> 1108     pop      PAL    groupXV  8815538      37 0.06250000 0.03125000
#> 1109     pop      SMR    groupXV  8815538      37 0.47357294 0.36363636
#> 1110     pop      UPD    groupXV  8815538      37 0.40923077 0.26923077
#> 1111 fam.pop    A.ASP    groupXI  9019877      38 0.50000000 0.25000000
#> 1112 fam.pop    A.CLF    groupXI  9019877      38 0.00000000 0.00000000
#> 1113 fam.pop    A.OPL    groupXI  9019877      38 0.25000000 0.12500000
#> 1114 fam.pop    A.PAL    groupXI  9019877      38 0.55555556 0.50000000
#> 1115 fam.pop    A.SMR    groupXI  9019877      38 0.53333333 0.40000000
#> 1116 fam.pop    A.UPD    groupXI  9019877      38 0.53333333 0.33333333
#> 1117 fam.pop    B.ASP    groupXI  9019877      38 0.53333333 0.66666667
#> 1118 fam.pop    B.CLF    groupXI  9019877      38 0.00000000 0.00000000
#> 1119 fam.pop    B.OPL    groupXI  9019877      38 0.40909091 0.25000000
#> 1120 fam.pop    B.PAL    groupXI  9019877      38 0.53571429 0.37500000
#> 1121 fam.pop    B.SMR    groupXI  9019877      38 0.46666667 0.30000000
#> 1122 fam.pop    B.UPD    groupXI  9019877      38 0.53571429 0.37500000
#> 1123 fam.pop    C.ASP    groupXI  9019877      38 0.53333333 0.33333333
#> 1124 fam.pop    C.CLF    groupXI  9019877      38 0.00000000 0.00000000
#> 1125 fam.pop    C.OPL    groupXI  9019877      38 0.26373626 0.14285714
#> 1126 fam.pop    C.PAL    groupXI  9019877      38 0.53333333 0.33333333
#> 1127 fam.pop    C.SMR    groupXI  9019877      38 0.53571429 0.37500000
#> 1128 fam.pop    C.UPD    groupXI  9019877      38 0.00000000 0.00000000
#> 1129 fam.pop    D.ASP    groupXI  9019877      38 0.66666667 0.50000000
#> 1130 fam.pop    D.CLF    groupXI  9019877      38 0.00000000 0.00000000
#> 1131 fam.pop    D.OPL    groupXI  9019877      38 0.40909091 0.25000000
#> 1132 fam.pop    D.PAL    groupXI  9019877      38 0.53571429 0.37500000
#> 1133 fam.pop    D.SMR    groupXI  9019877      38 0.53333333 0.60000000
#> 1134 fam.pop    D.UPD    groupXI  9019877      38 0.53333333 0.33333333
#> 1135     pop      ASP    groupXI  9019877      38 0.52105263 0.45000000
#> 1136     pop      CLF    groupXI  9019877      38 0.00000000 0.00000000
#> 1137     pop      OPL    groupXI  9019877      38 0.32173913 0.19565217
#> 1138     pop      PAL    groupXI  9019877      38 0.49798387 0.40625000
#> 1139     pop      SMR    groupXI  9019877      38 0.50071124 0.42105263
#> 1140     pop      UPD    groupXI  9019877      38 0.38888889 0.25000000
#> 1141 fam.pop    A.ASP   groupXXI  9067876      39 0.33333333 0.16666667
#> 1142 fam.pop    A.CLF   groupXXI  9067876      39 1.00000000 0.50000000
#> 1143 fam.pop    A.OPL   groupXXI  9067876      39 0.00000000 0.00000000
#> 1144 fam.pop    A.PAL   groupXXI  9067876      39 0.00000000 0.00000000
#> 1145 fam.pop    A.SMR   groupXXI  9067876      39 0.00000000 0.00000000
#> 1146 fam.pop    A.UPD   groupXXI  9067876      39 0.33333333 0.16666667
#> 1147 fam.pop    B.ASP   groupXXI  9067876      39 0.50000000 0.25000000
#> 1148 fam.pop    B.CLF   groupXXI  9067876      39        NaN         NA
#> 1149 fam.pop    B.OPL   groupXXI  9067876      39 0.00000000 0.00000000
#> 1150 fam.pop    B.PAL   groupXXI  9067876      39 0.33333333 0.16666667
#> 1151 fam.pop    B.SMR   groupXXI  9067876      39 0.33333333 0.16666667
#> 1152 fam.pop    B.UPD   groupXXI  9067876      39 0.00000000 0.00000000
#> 1153 fam.pop    C.ASP   groupXXI  9067876      39 0.33333333 0.16666667
#> 1154 fam.pop    C.CLF   groupXXI  9067876      39 1.00000000 0.50000000
#> 1155 fam.pop    C.OPL   groupXXI  9067876      39 0.33333333 0.16666667
#> 1156 fam.pop    C.PAL   groupXXI  9067876      39 0.25000000 0.12500000
#> 1157 fam.pop    C.SMR   groupXXI  9067876      39 0.25000000 0.12500000
#> 1158 fam.pop    C.UPD   groupXXI  9067876      39 0.33333333 0.16666667
#> 1159 fam.pop    D.ASP   groupXXI  9067876      39 0.00000000 0.00000000
#> 1160 fam.pop    D.CLF   groupXXI  9067876      39 0.00000000 0.00000000
#> 1161 fam.pop    D.OPL   groupXXI  9067876      39 0.35555556 0.20000000
#> 1162 fam.pop    D.PAL   groupXXI  9067876      39 1.00000000 0.50000000
#> 1163 fam.pop    D.SMR   groupXXI  9067876      39 0.25000000 0.12500000
#> 1164 fam.pop    D.UPD   groupXXI  9067876      39 0.50000000 0.25000000
#> 1165     pop      ASP   groupXXI  9067876      39 0.29411765 0.16666667
#> 1166     pop      CLF   groupXXI  9067876      39 0.53333333 0.33333333
#> 1167     pop      OPL   groupXXI  9067876      39 0.18620690 0.10000000
#> 1168     pop      PAL   groupXXI  9067876      39 0.24675325 0.13636364
#> 1169     pop      SMR   groupXXI  9067876      39 0.16577540 0.08823529
#> 1170     pop      UPD   groupXXI  9067876      39 0.26842105 0.15000000
#> 1171 fam.pop    A.ASP  groupXVII  9108184      40 0.53333333 0.66666667
#> 1172 fam.pop    A.CLF  groupXVII  9108184      40 0.00000000 0.00000000
#> 1173 fam.pop    A.OPL  groupXVII  9108184      40 0.25000000 0.12500000
#> 1174 fam.pop    A.PAL  groupXVII  9108184      40 0.00000000 0.00000000
#> 1175 fam.pop    A.SMR  groupXVII  9108184      40 0.42857143 0.25000000
#> 1176 fam.pop    A.UPD  groupXVII  9108184      40 0.53333333 0.33333333
#> 1177 fam.pop    B.ASP  groupXVII  9108184      40 0.53333333 0.33333333
#> 1178 fam.pop    B.CLF  groupXVII  9108184      40        NaN         NA
#> 1179 fam.pop    B.OPL  groupXVII  9108184      40 0.35555556 0.20000000
#> 1180 fam.pop    B.PAL  groupXVII  9108184      40 0.53571429 0.37500000
#> 1181 fam.pop    B.SMR  groupXVII  9108184      40 0.00000000 0.00000000
#> 1182 fam.pop    B.UPD  groupXVII  9108184      40 0.50000000 0.25000000
#> 1183 fam.pop    C.ASP  groupXVII  9108184      40 0.60000000 0.50000000
#> 1184 fam.pop    C.CLF  groupXVII  9108184      40 1.00000000 0.50000000
#> 1185 fam.pop    C.OPL  groupXVII  9108184      40 0.42857143 0.25000000
#> 1186 fam.pop    C.PAL  groupXVII  9108184      40 0.25000000 0.12500000
#> 1187 fam.pop    C.SMR  groupXVII  9108184      40 0.46666667 0.30000000
#> 1188 fam.pop    C.UPD  groupXVII  9108184      40 0.33333333 0.83333333
#> 1189 fam.pop    D.ASP  groupXVII  9108184      40 0.00000000 0.00000000
#> 1190 fam.pop    D.CLF  groupXVII  9108184      40 0.50000000 0.25000000
#> 1191 fam.pop    D.OPL  groupXVII  9108184      40 0.26373626 0.14285714
#> 1192 fam.pop    D.PAL  groupXVII  9108184      40 0.33333333 0.16666667
#> 1193 fam.pop    D.SMR  groupXVII  9108184      40 0.20000000 0.10000000
#> 1194 fam.pop    D.UPD  groupXVII  9108184      40 0.33333333 0.16666667
#> 1195     pop      ASP  groupXVII  9108184      40 0.52105263 0.45000000
#> 1196     pop      CLF  groupXVII  9108184      40 0.35555556 0.20000000
#> 1197     pop      OPL  groupXVII  9108184      40 0.29615385 0.17500000
#> 1198     pop      PAL  groupXVII  9108184      40 0.28735632 0.16666667
#> 1199     pop      SMR  groupXVII  9108184      40 0.27311522 0.15789474
#> 1200     pop      UPD  groupXVII  9108184      40 0.50649351 0.40909091
#> 1201 fam.pop    A.ASP    groupIX  9284427      41 0.60000000 0.50000000
#> 1202 fam.pop    A.CLF    groupIX  9284427      41 0.50000000 0.25000000
#> 1203 fam.pop    A.OPL    groupIX  9284427      41 0.53571429 0.37500000
#> 1204 fam.pop    A.PAL    groupIX  9284427      41 0.00000000 1.00000000
#> 1205 fam.pop    A.SMR    groupIX  9284427      41 0.42857143 0.25000000
#> 1206 fam.pop    A.UPD    groupIX  9284427      41 0.00000000 1.00000000
#> 1207 fam.pop    B.ASP    groupIX  9284427      41 0.33333333 0.83333333
#> 1208 fam.pop    B.CLF    groupIX  9284427      41 0.00000000 0.00000000
#> 1209 fam.pop    B.OPL    groupIX  9284427      41 0.46666667 0.30000000
#> 1210 fam.pop    B.PAL    groupIX  9284427      41 0.00000000 1.00000000
#> 1211 fam.pop    B.SMR    groupIX  9284427      41 0.48484848 0.33333333
#> 1212 fam.pop    B.UPD    groupIX  9284427      41 0.57142857 0.50000000
#> 1213 fam.pop    C.ASP    groupIX  9284427      41 0.60000000 0.50000000
#> 1214 fam.pop    C.CLF    groupIX  9284427      41 1.00000000 0.50000000
#> 1215 fam.pop    C.OPL    groupIX  9284427      41 0.20000000 0.10000000
#> 1216 fam.pop    C.PAL    groupIX  9284427      41 0.00000000 1.00000000
#> 1217 fam.pop    C.SMR    groupIX  9284427      41 0.60000000 0.50000000
#> 1218 fam.pop    C.UPD    groupIX  9284427      41 0.66666667 0.50000000
#> 1219 fam.pop    D.ASP    groupIX  9284427      41 0.66666667 0.50000000
#> 1220 fam.pop    D.CLF    groupIX  9284427      41 1.00000000 0.50000000
#> 1221 fam.pop    D.OPL    groupIX  9284427      41 0.46666667 0.30000000
#> 1222 fam.pop    D.PAL    groupIX  9284427      41 0.00000000 1.00000000
#> 1223 fam.pop    D.SMR    groupIX  9284427      41 0.42857143 0.25000000
#> 1224 fam.pop    D.UPD    groupIX  9284427      41 0.53333333 0.66666667
#> 1225     pop      ASP    groupIX  9284427      41 0.50649351 0.59090909
#> 1226     pop      CLF    groupIX  9284427      41 0.46666667 0.30000000
#> 1227     pop      OPL    groupIX  9284427      41 0.39829303 0.26315789
#> 1228     pop      PAL    groupIX  9284427      41 0.00000000 1.00000000
#> 1229     pop      SMR    groupIX  9284427      41 0.45098039 0.32352941
#> 1230     pop      UPD    groupIX  9284427      41 0.50526316 0.60000000
#> 1231 fam.pop    A.ASP     groupV  9767769      42 0.53333333 0.33333333
#> 1232 fam.pop    A.CLF     groupV  9767769      42 0.66666667 0.50000000
#> 1233 fam.pop    A.OPL     groupV  9767769      42 0.53030303 0.58333333
#> 1234 fam.pop    A.PAL     groupV  9767769      42 0.20000000 0.10000000
#> 1235 fam.pop    A.SMR     groupV  9767769      42 0.53333333 0.60000000
#> 1236 fam.pop    A.UPD     groupV  9767769      42 0.60000000 0.50000000
#> 1237 fam.pop    B.ASP     groupV  9767769      42 0.00000000 0.00000000
#> 1238 fam.pop    B.CLF     groupV  9767769      42 1.00000000 0.50000000
#> 1239 fam.pop    B.OPL     groupV  9767769      42 0.54545455 0.50000000
#> 1240 fam.pop    B.PAL     groupV  9767769      42 0.00000000 0.00000000
#> 1241 fam.pop    B.SMR     groupV  9767769      42 0.54545455 0.50000000
#> 1242 fam.pop    B.UPD     groupV  9767769      42 0.42857143 0.25000000
#> 1243 fam.pop    C.ASP     groupV  9767769      42 0.50000000 0.25000000
#> 1244 fam.pop    C.CLF     groupV  9767769      42 0.00000000 0.00000000
#> 1245 fam.pop    C.OPL     groupV  9767769      42 0.46666667 0.30000000
#> 1246 fam.pop    C.PAL     groupV  9767769      42 0.25000000 0.12500000
#> 1247 fam.pop    C.SMR     groupV  9767769      42 0.48484848 0.66666667
#> 1248 fam.pop    C.UPD     groupV  9767769      42 0.53333333 0.33333333
#> 1249 fam.pop    D.ASP     groupV  9767769      42 1.00000000 0.50000000
#> 1250 fam.pop    D.CLF     groupV  9767769      42 1.00000000 0.50000000
#> 1251 fam.pop    D.OPL     groupV  9767769      42 0.53030303 0.41666667
#> 1252 fam.pop    D.PAL     groupV  9767769      42 0.00000000 0.00000000
#> 1253 fam.pop    D.SMR     groupV  9767769      42 0.55555556 0.50000000
#> 1254 fam.pop    D.UPD     groupV  9767769      42 0.53571429 0.37500000
#> 1255     pop      ASP     groupV  9767769      42 0.40000000 0.25000000
#> 1256     pop      CLF     groupV  9767769      42 0.53333333 0.40000000
#> 1257     pop      OPL     groupV  9767769      42 0.50724638 0.45652174
#> 1258     pop      PAL     groupV  9767769      42 0.12096774 0.06250000
#> 1259     pop      SMR     groupV  9767769      42 0.50211416 0.56818182
#> 1260     pop      UPD     groupV  9767769      42 0.47619048 0.35714286
#> 1261 fam.pop    A.ASP   groupXII  9881312      43 0.33333333 0.16666667
#> 1262 fam.pop    A.CLF   groupXII  9881312      43 0.00000000 0.00000000
#> 1263 fam.pop    A.OPL   groupXII  9881312      43 0.25000000 0.12500000
#> 1264 fam.pop    A.PAL   groupXII  9881312      43 0.00000000 0.00000000
#> 1265 fam.pop    A.SMR   groupXII  9881312      43 0.16666667 0.08333333
#> 1266 fam.pop    A.UPD   groupXII  9881312      43 0.53333333 0.33333333
#> 1267 fam.pop    B.ASP   groupXII  9881312      43 0.33333333 0.16666667
#> 1268 fam.pop    B.CLF   groupXII  9881312      43 0.00000000 0.00000000
#> 1269 fam.pop    B.OPL   groupXII  9881312      43 0.25000000 0.12500000
#> 1270 fam.pop    B.PAL   groupXII  9881312      43 0.00000000 0.00000000
#> 1271 fam.pop    B.SMR   groupXII  9881312      43 0.16666667 0.08333333
#> 1272 fam.pop    B.UPD   groupXII  9881312      43 0.00000000 0.00000000
#> 1273 fam.pop    C.ASP   groupXII  9881312      43 0.00000000 0.00000000
#> 1274 fam.pop    C.CLF   groupXII  9881312      43 0.00000000 0.00000000
#> 1275 fam.pop    C.OPL   groupXII  9881312      43 0.40909091 0.25000000
#> 1276 fam.pop    C.PAL   groupXII  9881312      43 0.00000000 0.00000000
#> 1277 fam.pop    C.SMR   groupXII  9881312      43 0.00000000 0.00000000
#> 1278 fam.pop    C.UPD   groupXII  9881312      43 0.42857143 0.25000000
#> 1279 fam.pop    D.ASP   groupXII  9881312      43 0.00000000 0.00000000
#> 1280 fam.pop    D.CLF   groupXII  9881312      43 0.66666667 0.50000000
#> 1281 fam.pop    D.OPL   groupXII  9881312      43 0.36263736 0.21428571
#> 1282 fam.pop    D.PAL   groupXII  9881312      43 0.00000000 0.00000000
#> 1283 fam.pop    D.SMR   groupXII  9881312      43 0.40909091 0.25000000
#> 1284 fam.pop    D.UPD   groupXII  9881312      43 0.42857143 0.25000000
#> 1285     pop      ASP   groupXII  9881312      43 0.17316017 0.09090909
#> 1286     pop      CLF   groupXII  9881312      43 0.35555556 0.20000000
#> 1287     pop      OPL   groupXII  9881312      43 0.31591173 0.19047619
#> 1288     pop      PAL   groupXII  9881312      43 0.00000000 0.00000000
#> 1289     pop      SMR   groupXII  9881312      43 0.19060284 0.10416667
#> 1290     pop      UPD   groupXII  9881312      43 0.36923077 0.23076923
#> 1291 fam.pop    A.ASP    groupXX 10296516      44 0.50000000 0.25000000
#> 1292 fam.pop    A.CLF    groupXX 10296516      44 0.00000000 0.00000000
#> 1293 fam.pop    A.OPL    groupXX 10296516      44 0.20000000 0.10000000
#> 1294 fam.pop    A.PAL    groupXX 10296516      44 0.46666667 0.30000000
#> 1295 fam.pop    A.SMR    groupXX 10296516      44 0.00000000 0.00000000
#> 1296 fam.pop    A.UPD    groupXX 10296516      44 0.00000000 0.00000000
#> 1297 fam.pop    B.ASP    groupXX 10296516      44 0.53333333 0.33333333
#> 1298 fam.pop    B.CLF    groupXX 10296516      44 0.00000000 0.00000000
#> 1299 fam.pop    B.OPL    groupXX 10296516      44 0.00000000 0.00000000
#> 1300 fam.pop    B.PAL    groupXX 10296516      44 0.53333333 0.66666667
#> 1301 fam.pop    B.SMR    groupXX 10296516      44 0.46666667 0.30000000
#> 1302 fam.pop    B.UPD    groupXX 10296516      44 0.25000000 0.12500000
#> 1303 fam.pop    C.ASP    groupXX 10296516      44 0.00000000 0.00000000
#> 1304 fam.pop    C.CLF    groupXX 10296516      44 0.00000000 0.00000000
#> 1305 fam.pop    C.OPL    groupXX 10296516      44 0.14285714 0.07142857
#> 1306 fam.pop    C.PAL    groupXX 10296516      44 0.53333333 0.66666667
#> 1307 fam.pop    C.SMR    groupXX 10296516      44 0.20000000 0.10000000
#> 1308 fam.pop    C.UPD    groupXX 10296516      44 0.60000000 0.50000000
#> 1309 fam.pop    D.ASP    groupXX 10296516      44 0.00000000 0.00000000
#> 1310 fam.pop    D.CLF    groupXX 10296516      44 0.00000000 0.00000000
#> 1311 fam.pop    D.OPL    groupXX 10296516      44 0.14285714 0.07142857
#> 1312 fam.pop    D.PAL    groupXX 10296516      44 0.25000000 0.12500000
#> 1313 fam.pop    D.SMR    groupXX 10296516      44 0.16666667 0.08333333
#> 1314 fam.pop    D.UPD    groupXX 10296516      44 0.25000000 0.12500000
#> 1315     pop      ASP    groupXX 10296516      44 0.29411765 0.16666667
#> 1316     pop      CLF    groupXX 10296516      44 0.00000000 0.00000000
#> 1317     pop      OPL    groupXX 10296516      44 0.11085973 0.05769231
#> 1318     pop      PAL    groupXX 10296516      44 0.49655172 0.40000000
#> 1319     pop      SMR    groupXX 10296516      44 0.21486643 0.11904762
#> 1320     pop      UPD    groupXX 10296516      44 0.34420290 0.20833333
#> 1321 fam.pop    A.ASP     groupI 10593829      45 0.00000000 0.00000000
#> 1322 fam.pop    A.CLF     groupI 10593829      45 0.00000000 0.00000000
#> 1323 fam.pop    A.OPL     groupI 10593829      45 0.00000000 0.00000000
#> 1324 fam.pop    A.PAL     groupI 10593829      45 0.00000000 0.00000000
#> 1325 fam.pop    A.SMR     groupI 10593829      45 0.00000000 0.00000000
#> 1326 fam.pop    A.UPD     groupI 10593829      45 0.00000000 0.00000000
#> 1327 fam.pop    B.ASP     groupI 10593829      45 0.00000000 0.00000000
#> 1328 fam.pop    B.CLF     groupI 10593829      45 0.00000000 0.00000000
#> 1329 fam.pop    B.OPL     groupI 10593829      45 0.00000000 0.00000000
#> 1330 fam.pop    B.PAL     groupI 10593829      45 0.00000000 0.00000000
#> 1331 fam.pop    B.SMR     groupI 10593829      45 0.00000000 0.00000000
#> 1332 fam.pop    B.UPD     groupI 10593829      45 0.25000000 0.12500000
#> 1333 fam.pop    C.ASP     groupI 10593829      45 0.00000000 0.00000000
#> 1334 fam.pop    C.CLF     groupI 10593829      45 0.00000000 0.00000000
#> 1335 fam.pop    C.OPL     groupI 10593829      45 0.00000000 0.00000000
#> 1336 fam.pop    C.PAL     groupI 10593829      45 0.00000000 0.00000000
#> 1337 fam.pop    C.SMR     groupI 10593829      45 0.00000000 0.00000000
#> 1338 fam.pop    C.UPD     groupI 10593829      45 0.00000000 0.00000000
#> 1339 fam.pop    D.ASP     groupI 10593829      45 0.00000000 0.00000000
#> 1340 fam.pop    D.CLF     groupI 10593829      45 0.00000000 0.00000000
#> 1341 fam.pop    D.OPL     groupI 10593829      45 0.00000000 0.00000000
#> 1342 fam.pop    D.PAL     groupI 10593829      45 0.00000000 0.00000000
#> 1343 fam.pop    D.SMR     groupI 10593829      45 0.00000000 0.00000000
#> 1344 fam.pop    D.UPD     groupI 10593829      45 0.25000000 0.12500000
#> 1345     pop      ASP     groupI 10593829      45 0.00000000 0.00000000
#> 1346     pop      CLF     groupI 10593829      45 0.00000000 0.00000000
#> 1347     pop      OPL     groupI 10593829      45 0.00000000 0.00000000
#> 1348     pop      PAL     groupI 10593829      45 0.00000000 0.00000000
#> 1349     pop      SMR     groupI 10593829      45 0.00000000 0.00000000
#> 1350     pop      UPD     groupI 10593829      45 0.12873563 0.06666667
#> 1351 fam.pop    A.ASP   groupXIX 10765968      46 0.00000000 0.00000000
#> 1352 fam.pop    A.CLF   groupXIX 10765968      46 0.50000000 0.25000000
#> 1353 fam.pop    A.OPL   groupXIX 10765968      46 0.42857143 0.25000000
#> 1354 fam.pop    A.PAL   groupXIX 10765968      46 0.00000000 0.00000000
#> 1355 fam.pop    A.SMR   groupXIX 10765968      46 0.46666667 0.30000000
#> 1356 fam.pop    A.UPD   groupXIX 10765968      46 0.66666667 0.50000000
#> 1357 fam.pop    B.ASP   groupXIX 10765968      46 0.00000000 0.00000000
#> 1358 fam.pop    B.CLF   groupXIX 10765968      46 0.00000000 0.00000000
#> 1359 fam.pop    B.OPL   groupXIX 10765968      46 0.48484848 0.33333333
#> 1360 fam.pop    B.PAL   groupXIX 10765968      46 0.00000000 0.00000000
#> 1361 fam.pop    B.SMR   groupXIX 10765968      46 0.25000000 0.12500000
#> 1362 fam.pop    B.UPD   groupXIX 10765968      46 0.53571429 0.37500000
#> 1363 fam.pop    C.ASP   groupXIX 10765968      46 0.50000000 0.75000000
#> 1364 fam.pop    C.CLF   groupXIX 10765968      46        NaN         NA
#> 1365 fam.pop    C.OPL   groupXIX 10765968      46 0.60000000 0.50000000
#> 1366 fam.pop    C.PAL   groupXIX 10765968      46 0.00000000 0.00000000
#> 1367 fam.pop    C.SMR   groupXIX 10765968      46 0.50000000 0.25000000
#> 1368 fam.pop    C.UPD   groupXIX 10765968      46 0.33333333 0.16666667
#> 1369 fam.pop    D.ASP   groupXIX 10765968      46 0.00000000 1.00000000
#> 1370 fam.pop    D.CLF   groupXIX 10765968      46 0.50000000 0.25000000
#> 1371 fam.pop    D.OPL   groupXIX 10765968      46 0.53333333 0.60000000
#> 1372 fam.pop    D.PAL   groupXIX 10765968      46 0.00000000 0.00000000
#> 1373 fam.pop    D.SMR   groupXIX 10765968      46 0.30303030 0.16666667
#> 1374 fam.pop    D.UPD   groupXIX 10765968      46 0.53571429 0.37500000
#> 1375     pop      ASP   groupXIX 10765968      46 0.52500000 0.43750000
#> 1376     pop      CLF   groupXIX 10765968      46 0.35555556 0.20000000
#> 1377     pop      OPL   groupXIX 10765968      46 0.50000000 0.41666667
#> 1378     pop      PAL   groupXIX 10765968      46 0.00000000 0.00000000
#> 1379     pop      SMR   groupXIX 10765968      46 0.33689840 0.20588235
#> 1380     pop      UPD   groupXIX 10765968      46 0.47076923 0.34615385
#> 1381 fam.pop    A.ASP    groupII 10829169      47 0.33333333 0.16666667
#> 1382 fam.pop    A.CLF    groupII 10829169      47 0.00000000 0.00000000
#> 1383 fam.pop    A.OPL    groupII 10829169      47 0.40909091 0.75000000
#> 1384 fam.pop    A.PAL    groupII 10829169      47 0.00000000 0.00000000
#> 1385 fam.pop    A.SMR    groupII 10829169      47 0.48484848 0.33333333
#> 1386 fam.pop    A.UPD    groupII 10829169      47 0.60000000 0.50000000
#> 1387 fam.pop    B.ASP    groupII 10829169      47 0.33333333 0.16666667
#> 1388 fam.pop    B.CLF    groupII 10829169      47 0.00000000 0.00000000
#> 1389 fam.pop    B.OPL    groupII 10829169      47 0.48484848 0.33333333
#> 1390 fam.pop    B.PAL    groupII 10829169      47 0.33333333 0.16666667
#> 1391 fam.pop    B.SMR    groupII 10829169      47 0.16666667 0.08333333
#> 1392 fam.pop    B.UPD    groupII 10829169      47 0.00000000 0.00000000
#> 1393 fam.pop    C.ASP    groupII 10829169      47 0.33333333 0.16666667
#> 1394 fam.pop    C.CLF    groupII 10829169      47 0.00000000 0.00000000
#> 1395 fam.pop    C.OPL    groupII 10829169      47 0.53030303 0.41666667
#> 1396 fam.pop    C.PAL    groupII 10829169      47 0.25000000 0.12500000
#> 1397 fam.pop    C.SMR    groupII 10829169      47 0.40909091 0.25000000
#> 1398 fam.pop    C.UPD    groupII 10829169      47 0.00000000 0.00000000
#> 1399 fam.pop    D.ASP    groupII 10829169      47 0.66666667 0.50000000
#> 1400 fam.pop    D.CLF    groupII 10829169      47 0.00000000 0.00000000
#> 1401 fam.pop    D.OPL    groupII 10829169      47 0.53333333 0.40000000
#> 1402 fam.pop    D.PAL    groupII 10829169      47 0.25000000 0.12500000
#> 1403 fam.pop    D.SMR    groupII 10829169      47 0.46666667 0.30000000
#> 1404 fam.pop    D.UPD    groupII 10829169      47 0.42857143 0.25000000
#> 1405     pop      ASP    groupII 10829169      47 0.36796537 0.22727273
#> 1406     pop      CLF    groupII 10829169      47 0.00000000 0.00000000
#> 1407     pop      OPL    groupII 10829169      47 0.51014493 0.47826087
#> 1408     pop      PAL    groupII 10829169      47 0.18620690 0.10000000
#> 1409     pop      SMR    groupII 10829169      47 0.37198068 0.23913043
#> 1410     pop      UPD    groupII 10829169      47 0.28735632 0.16666667
#> 1411 fam.pop    A.ASP    groupXX 10913115      48 0.50000000 0.25000000
#> 1412 fam.pop    A.CLF    groupXX 10913115      48 0.50000000 0.75000000
#> 1413 fam.pop    A.OPL    groupXX 10913115      48 0.55555556 0.50000000
#> 1414 fam.pop    A.PAL    groupXX 10913115      48 0.53571429 0.62500000
#> 1415 fam.pop    A.SMR    groupXX 10913115      48 0.53571429 0.62500000
#> 1416 fam.pop    A.UPD    groupXX 10913115      48 0.60000000 0.50000000
#> 1417 fam.pop    B.ASP    groupXX 10913115      48 0.00000000 0.00000000
#> 1418 fam.pop    B.CLF    groupXX 10913115      48 0.00000000 0.00000000
#> 1419 fam.pop    B.OPL    groupXX 10913115      48 0.53333333 0.60000000
#> 1420 fam.pop    B.PAL    groupXX 10913115      48 0.53571429 0.37500000
#> 1421 fam.pop    B.SMR    groupXX 10913115      48 0.53030303 0.41666667
#> 1422 fam.pop    B.UPD    groupXX 10913115      48 0.60000000 0.50000000
#> 1423 fam.pop    C.ASP    groupXX 10913115      48 0.53333333 0.33333333
#> 1424 fam.pop    C.CLF    groupXX 10913115      48 1.00000000 0.50000000
#> 1425 fam.pop    C.OPL    groupXX 10913115      48 0.53030303 0.58333333
#> 1426 fam.pop    C.PAL    groupXX 10913115      48 0.42857143 0.25000000
#> 1427 fam.pop    C.SMR    groupXX 10913115      48 0.33333333 0.16666667
#> 1428 fam.pop    C.UPD    groupXX 10913115      48 0.53333333 0.66666667
#> 1429 fam.pop    D.ASP    groupXX 10913115      48 0.00000000 1.00000000
#> 1430 fam.pop    D.CLF    groupXX 10913115      48 0.66666667 0.50000000
#> 1431 fam.pop    D.OPL    groupXX 10913115      48 0.53571429 0.62500000
#> 1432 fam.pop    D.PAL    groupXX 10913115      48 0.50000000 0.25000000
#> 1433 fam.pop    D.SMR    groupXX 10913115      48 0.55555556 0.50000000
#> 1434 fam.pop    D.UPD    groupXX 10913115      48 0.00000000 0.00000000
#> 1435     pop      ASP    groupXX 10913115      48 0.42483660 0.27777778
#> 1436     pop      CLF    groupXX 10913115      48 0.54545455 0.50000000
#> 1437     pop      OPL    groupXX 10913115      48 0.50128205 0.57500000
#> 1438     pop      PAL    groupXX 10913115      48 0.49470899 0.39285714
#> 1439     pop      SMR    groupXX 10913115      48 0.50793651 0.44444444
#> 1440     pop      UPD    groupXX 10913115      48 0.49230769 0.38461538
#> 1441 fam.pop    A.ASP groupXVIII 10974594      49 0.00000000 0.00000000
#> 1442 fam.pop    A.CLF groupXVIII 10974594      49 0.00000000 0.00000000
#> 1443 fam.pop    A.OPL groupXVIII 10974594      49 0.00000000 0.00000000
#> 1444 fam.pop    A.PAL groupXVIII 10974594      49 0.00000000 0.00000000
#> 1445 fam.pop    A.SMR groupXVIII 10974594      49 0.00000000 0.00000000
#> 1446 fam.pop    A.UPD groupXVIII 10974594      49 0.33333333 0.16666667
#> 1447 fam.pop    B.ASP groupXVIII 10974594      49 0.00000000 0.00000000
#> 1448 fam.pop    B.CLF groupXVIII 10974594      49 1.00000000 0.50000000
#> 1449 fam.pop    B.OPL groupXVIII 10974594      49 0.14285714 0.07142857
#> 1450 fam.pop    B.PAL groupXVIII 10974594      49 0.00000000 0.00000000
#> 1451 fam.pop    B.SMR groupXVIII 10974594      49 0.00000000 0.00000000
#> 1452 fam.pop    B.UPD groupXVIII 10974594      49 0.25000000 0.12500000
#> 1453 fam.pop    C.ASP groupXVIII 10974594      49 0.00000000 0.00000000
#> 1454 fam.pop    C.CLF groupXVIII 10974594      49 0.00000000 0.00000000
#> 1455 fam.pop    C.OPL groupXVIII 10974594      49 0.00000000 0.00000000
#> 1456 fam.pop    C.PAL groupXVIII 10974594      49 0.00000000 0.00000000
#> 1457 fam.pop    C.SMR groupXVIII 10974594      49 0.00000000 0.00000000
#> 1458 fam.pop    C.UPD groupXVIII 10974594      49 0.00000000 0.00000000
#> 1459 fam.pop    D.ASP groupXVIII 10974594      49 0.00000000 0.00000000
#> 1460 fam.pop    D.CLF groupXVIII 10974594      49 0.50000000 0.25000000
#> 1461 fam.pop    D.OPL groupXVIII 10974594      49 0.00000000 0.00000000
#> 1462 fam.pop    D.PAL groupXVIII 10974594      49 0.00000000 0.00000000
#> 1463 fam.pop    D.SMR groupXVIII 10974594      49 0.00000000 0.00000000
#> 1464 fam.pop    D.UPD groupXVIII 10974594      49 0.00000000 0.00000000
#> 1465     pop      ASP groupXVIII 10974594      49 0.00000000 0.00000000
#> 1466     pop      CLF groupXVIII 10974594      49 0.30303030 0.16666667
#> 1467     pop      OPL groupXVIII 10974594      49 0.03703704 0.01851852
#> 1468     pop      PAL groupXVIII 10974594      49 0.00000000 0.00000000
#> 1469     pop      SMR groupXVIII 10974594      49 0.00000000 0.00000000
#> 1470     pop      UPD groupXVIII 10974594      49 0.14769231 0.07692308
#> 1471 fam.pop    A.ASP  groupXVII 11256661      50 0.50000000 0.25000000
#> 1472 fam.pop    A.CLF  groupXVII 11256661      50 0.00000000 0.00000000
#> 1473 fam.pop    A.OPL  groupXVII 11256661      50 0.20000000 0.10000000
#> 1474 fam.pop    A.PAL  groupXVII 11256661      50 0.00000000 0.00000000
#> 1475 fam.pop    A.SMR  groupXVII 11256661      50 0.00000000 0.00000000
#> 1476 fam.pop    A.UPD  groupXVII 11256661      50 0.00000000 0.00000000
#> 1477 fam.pop    B.ASP  groupXVII 11256661      50 0.33333333 0.16666667
#> 1478 fam.pop    B.CLF  groupXVII 11256661      50 0.00000000 0.00000000
#> 1479 fam.pop    B.OPL  groupXVII 11256661      50 0.40909091 0.25000000
#> 1480 fam.pop    B.PAL  groupXVII 11256661      50 0.25000000 0.12500000
#> 1481 fam.pop    B.SMR  groupXVII 11256661      50 0.35555556 0.20000000
#> 1482 fam.pop    B.UPD  groupXVII 11256661      50 0.00000000 0.00000000
#> 1483 fam.pop    C.ASP  groupXVII 11256661      50 0.33333333 0.16666667
#> 1484 fam.pop    C.CLF  groupXVII 11256661      50 0.00000000 0.00000000
#> 1485 fam.pop    C.OPL  groupXVII 11256661      50 0.00000000 0.00000000
#> 1486 fam.pop    C.PAL  groupXVII 11256661      50 0.00000000 0.00000000
#> 1487 fam.pop    C.SMR  groupXVII 11256661      50 0.00000000 0.00000000
#> 1488 fam.pop    C.UPD  groupXVII 11256661      50 0.00000000 0.00000000
#> 1489 fam.pop    D.ASP  groupXVII 11256661      50 0.00000000 0.00000000
#> 1490 fam.pop    D.CLF  groupXVII 11256661      50 0.50000000 0.25000000
#> 1491 fam.pop    D.OPL  groupXVII 11256661      50 0.16666667 0.08333333
#> 1492 fam.pop    D.PAL  groupXVII 11256661      50 0.00000000 0.00000000
#> 1493 fam.pop    D.SMR  groupXVII 11256661      50 0.00000000 0.00000000
#> 1494 fam.pop    D.UPD  groupXVII 11256661      50 0.53333333 0.33333333
#> 1495     pop      ASP  groupXVII 11256661      50 0.26842105 0.15000000
#> 1496     pop      CLF  groupXVII 11256661      50 0.20000000 0.10000000
#> 1497     pop      OPL  groupXVII 11256661      50 0.19806763 0.10869565
#> 1498     pop      PAL  groupXVII 11256661      50 0.06250000 0.03125000
#> 1499     pop      SMR  groupXVII 11256661      50 0.10241821 0.05263158
#> 1500     pop      UPD  groupXVII 11256661      50 0.13756614 0.07142857
#> 1501 fam.pop    A.ASP  groupVIII 11438870      51 0.33333333 0.16666667
#> 1502 fam.pop    A.CLF  groupVIII 11438870      51 0.00000000 0.00000000
#> 1503 fam.pop    A.OPL  groupVIII 11438870      51 0.46666667 0.30000000
#> 1504 fam.pop    A.PAL  groupVIII 11438870      51 0.50000000 0.25000000
#> 1505 fam.pop    A.SMR  groupVIII 11438870      51 0.42857143 0.25000000
#> 1506 fam.pop    A.UPD  groupVIII 11438870      51 0.60000000 0.50000000
#> 1507 fam.pop    B.ASP  groupVIII 11438870      51 0.50000000 0.75000000
#> 1508 fam.pop    B.CLF  groupVIII 11438870      51 0.00000000 0.00000000
#> 1509 fam.pop    B.OPL  groupVIII 11438870      51 0.40909091 0.25000000
#> 1510 fam.pop    B.PAL  groupVIII 11438870      51 0.33333333 0.16666667
#> 1511 fam.pop    B.SMR  groupVIII 11438870      51 0.33333333 0.16666667
#> 1512 fam.pop    B.UPD  groupVIII 11438870      51 0.53571429 0.37500000
#> 1513 fam.pop    C.ASP  groupVIII 11438870      51 0.53333333 0.33333333
#> 1514 fam.pop    C.CLF  groupVIII 11438870      51 0.00000000 0.00000000
#> 1515 fam.pop    C.OPL  groupVIII 11438870      51 0.55555556 0.50000000
#> 1516 fam.pop    C.PAL  groupVIII 11438870      51 0.33333333 0.16666667
#> 1517 fam.pop    C.SMR  groupVIII 11438870      51 0.50000000 0.25000000
#> 1518 fam.pop    C.UPD  groupVIII 11438870      51 0.42857143 0.25000000
#> 1519 fam.pop    D.ASP  groupVIII 11438870      51 1.00000000 0.50000000
#> 1520 fam.pop    D.CLF  groupVIII 11438870      51 0.00000000 0.00000000
#> 1521 fam.pop    D.OPL  groupVIII 11438870      51 0.53571429 0.37500000
#> 1522 fam.pop    D.PAL  groupVIII 11438870      51 0.33333333 0.16666667
#> 1523 fam.pop    D.SMR  groupVIII 11438870      51 0.55555556 0.50000000
#> 1524 fam.pop    D.UPD  groupVIII 11438870      51 0.60000000 0.50000000
#> 1525     pop      ASP  groupVIII 11438870      51 0.50326797 0.38888889
#> 1526     pop      CLF  groupVIII 11438870      51 0.00000000 0.00000000
#> 1527     pop      OPL  groupVIII 11438870      51 0.46666667 0.35000000
#> 1528     pop      PAL  groupVIII 11438870      51 0.31168831 0.18181818
#> 1529     pop      SMR  groupVIII 11438870      51 0.45238095 0.32142857
#> 1530     pop      UPD  groupVIII 11438870      51 0.49470899 0.39285714
#> 1531 fam.pop    A.ASP   groupVII 11477223      52 0.53333333 0.33333333
#> 1532 fam.pop    A.CLF   groupVII 11477223      52 0.66666667 0.50000000
#> 1533 fam.pop    A.OPL   groupVII 11477223      52 0.46666667 0.30000000
#> 1534 fam.pop    A.PAL   groupVII 11477223      52 0.53333333 0.60000000
#> 1535 fam.pop    A.SMR   groupVII 11477223      52 0.16666667 0.08333333
#> 1536 fam.pop    A.UPD   groupVII 11477223      52 0.60000000 0.50000000
#> 1537 fam.pop    B.ASP   groupVII 11477223      52 0.60000000 0.50000000
#> 1538 fam.pop    B.CLF   groupVII 11477223      52 0.00000000 0.00000000
#> 1539 fam.pop    B.OPL   groupVII 11477223      52 0.40909091 0.25000000
#> 1540 fam.pop    B.PAL   groupVII 11477223      52 0.53333333 0.33333333
#> 1541 fam.pop    B.SMR   groupVII 11477223      52 0.16666667 0.08333333
#> 1542 fam.pop    B.UPD   groupVII 11477223      52 0.53571429 0.37500000
#> 1543 fam.pop    C.ASP   groupVII 11477223      52 0.53333333 0.33333333
#> 1544 fam.pop    C.CLF   groupVII 11477223      52 0.00000000 0.00000000
#> 1545 fam.pop    C.OPL   groupVII 11477223      52 0.53030303 0.41666667
#> 1546 fam.pop    C.PAL   groupVII 11477223      52 1.00000000 0.50000000
#> 1547 fam.pop    C.SMR   groupVII 11477223      52 0.30303030 0.16666667
#> 1548 fam.pop    C.UPD   groupVII 11477223      52 0.53333333 0.33333333
#> 1549 fam.pop    D.ASP   groupVII 11477223      52 0.50000000 0.25000000
#> 1550 fam.pop    D.CLF   groupVII 11477223      52 0.66666667 0.50000000
#> 1551 fam.pop    D.OPL   groupVII 11477223      52 0.49450549 0.35714286
#> 1552 fam.pop    D.PAL   groupVII 11477223      52 0.53333333 0.33333333
#> 1553 fam.pop    D.SMR   groupVII 11477223      52 0.40909091 0.25000000
#> 1554 fam.pop    D.UPD   groupVII 11477223      52 0.42857143 0.25000000
#> 1555     pop      ASP   groupVII 11477223      52 0.48484848 0.36363636
#> 1556     pop      CLF   groupVII 11477223      52 0.48484848 0.33333333
#> 1557     pop      OPL   groupVII 11477223      52 0.45390071 0.33333333
#> 1558     pop      PAL   groupVII 11477223      52 0.51811594 0.45833333
#> 1559     pop      SMR   groupVII 11477223      52 0.25443262 0.14583333
#> 1560     pop      UPD   groupVII 11477223      52 0.47619048 0.35714286
#> 1561 fam.pop    A.ASP    groupII 11505084      53 0.33333333 0.16666667
#> 1562 fam.pop    A.CLF    groupII 11505084      53 1.00000000 0.50000000
#> 1563 fam.pop    A.OPL    groupII 11505084      53 0.46666667 0.30000000
#> 1564 fam.pop    A.PAL    groupII 11505084      53 0.00000000 0.00000000
#> 1565 fam.pop    A.SMR    groupII 11505084      53 0.00000000 0.00000000
#> 1566 fam.pop    A.UPD    groupII 11505084      53 0.66666667 0.50000000
#> 1567 fam.pop    B.ASP    groupII 11505084      53 0.33333333 0.16666667
#> 1568 fam.pop    B.CLF    groupII 11505084      53 0.00000000 0.00000000
#> 1569 fam.pop    B.OPL    groupII 11505084      53 0.53571429 0.62500000
#> 1570 fam.pop    B.PAL    groupII 11505084      53 0.00000000 0.00000000
#> 1571 fam.pop    B.SMR    groupII 11505084      53 0.16666667 0.08333333
#> 1572 fam.pop    B.UPD    groupII 11505084      53 0.33333333 0.16666667
#> 1573 fam.pop    C.ASP    groupII 11505084      53 0.53333333 0.33333333
#> 1574 fam.pop    C.CLF    groupII 11505084      53        NaN         NA
#> 1575 fam.pop    C.OPL    groupII 11505084      53 0.35555556 0.20000000
#> 1576 fam.pop    C.PAL    groupII 11505084      53 0.00000000 0.00000000
#> 1577 fam.pop    C.SMR    groupII 11505084      53 0.30303030 0.16666667
#> 1578 fam.pop    C.UPD    groupII 11505084      53 0.53333333 0.66666667
#> 1579 fam.pop    D.ASP    groupII 11505084      53 0.00000000 0.00000000
#> 1580 fam.pop    D.CLF    groupII 11505084      53 0.66666667 0.50000000
#> 1581 fam.pop    D.OPL    groupII 11505084      53 0.36263736 0.21428571
#> 1582 fam.pop    D.PAL    groupII 11505084      53 0.00000000 0.00000000
#> 1583 fam.pop    D.SMR    groupII 11505084      53 0.53030303 0.41666667
#> 1584 fam.pop    D.UPD    groupII 11505084      53 0.53571429 0.37500000
#> 1585     pop      ASP    groupII 11505084      53 0.33684211 0.20000000
#> 1586     pop      CLF    groupII 11505084      53 0.53571429 0.37500000
#> 1587     pop      OPL    groupII 11505084      53 0.43786295 0.30952381
#> 1588     pop      PAL    groupII 11505084      53 0.00000000 0.00000000
#> 1589     pop      SMR    groupII 11505084      53 0.29371981 0.17391304
#> 1590     pop      UPD    groupII 11505084      53 0.50724638 0.41666667
#> 1591 fam.pop    A.ASP    groupIX 11904651      54 0.00000000 0.00000000
#> 1592 fam.pop    A.CLF    groupIX 11904651      54 0.00000000 0.00000000
#> 1593 fam.pop    A.OPL    groupIX 11904651      54 0.16666667 0.08333333
#> 1594 fam.pop    A.PAL    groupIX 11904651      54 0.00000000 0.00000000
#> 1595 fam.pop    A.SMR    groupIX 11904651      54 0.16666667 0.08333333
#> 1596 fam.pop    A.UPD    groupIX 11904651      54 0.00000000 0.00000000
#> 1597 fam.pop    B.ASP    groupIX 11904651      54 0.00000000 0.00000000
#> 1598 fam.pop    B.CLF    groupIX 11904651      54 0.00000000 0.00000000
#> 1599 fam.pop    B.OPL    groupIX 11904651      54 0.00000000 0.00000000
#> 1600 fam.pop    B.PAL    groupIX 11904651      54 0.00000000 0.00000000
#> 1601 fam.pop    B.SMR    groupIX 11904651      54 0.00000000 0.00000000
#> 1602 fam.pop    B.UPD    groupIX 11904651      54 0.00000000 0.00000000
#> 1603 fam.pop    C.ASP    groupIX 11904651      54 0.00000000 0.00000000
#> 1604 fam.pop    C.CLF    groupIX 11904651      54 1.00000000 0.50000000
#> 1605 fam.pop    C.OPL    groupIX 11904651      54 0.00000000 0.00000000
#> 1606 fam.pop    C.PAL    groupIX 11904651      54 0.00000000 0.00000000
#> 1607 fam.pop    C.SMR    groupIX 11904651      54 0.00000000 0.00000000
#> 1608 fam.pop    C.UPD    groupIX 11904651      54 0.00000000 0.00000000
#> 1609 fam.pop    D.ASP    groupIX 11904651      54 0.00000000 0.00000000
#> 1610 fam.pop    D.CLF    groupIX 11904651      54 0.00000000 0.00000000
#> 1611 fam.pop    D.OPL    groupIX 11904651      54 0.00000000 0.00000000
#> 1612 fam.pop    D.PAL    groupIX 11904651      54 0.00000000 0.00000000
#> 1613 fam.pop    D.SMR    groupIX 11904651      54 0.00000000 0.00000000
#> 1614 fam.pop    D.UPD    groupIX 11904651      54 0.00000000 0.00000000
#> 1615     pop      ASP    groupIX 11904651      54 0.00000000 0.00000000
#> 1616     pop      CLF    groupIX 11904651      54 0.16666667 0.08333333
#> 1617     pop      OPL    groupIX 11904651      54 0.03846154 0.01923077
#> 1618     pop      PAL    groupIX 11904651      54 0.00000000 0.00000000
#> 1619     pop      SMR    groupIX 11904651      54 0.04347826 0.02173913
#> 1620     pop      UPD    groupIX 11904651      54 0.00000000 0.00000000
#> 1621 fam.pop    A.ASP  groupVIII 12108945      55 0.53333333 0.33333333
#> 1622 fam.pop    A.CLF  groupVIII 12108945      55 0.50000000 0.25000000
#> 1623 fam.pop    A.OPL  groupVIII 12108945      55 0.55555556 0.50000000
#> 1624 fam.pop    A.PAL  groupVIII 12108945      55 0.53571429 0.62500000
#> 1625 fam.pop    A.SMR  groupVIII 12108945      55 0.16666667 0.08333333
#> 1626 fam.pop    A.UPD  groupVIII 12108945      55 0.33333333 0.83333333
#> 1627 fam.pop    B.ASP  groupVIII 12108945      55 0.50000000 0.75000000
#> 1628 fam.pop    B.CLF  groupVIII 12108945      55 0.00000000 0.00000000
#> 1629 fam.pop    B.OPL  groupVIII 12108945      55 0.40909091 0.25000000
#> 1630 fam.pop    B.PAL  groupVIII 12108945      55 0.60000000 0.50000000
#> 1631 fam.pop    B.SMR  groupVIII 12108945      55 0.46666667 0.30000000
#> 1632 fam.pop    B.UPD  groupVIII 12108945      55 0.53571429 0.37500000
#> 1633 fam.pop    C.ASP  groupVIII 12108945      55 0.66666667 0.50000000
#> 1634 fam.pop    C.CLF  groupVIII 12108945      55 1.00000000 0.50000000
#> 1635 fam.pop    C.OPL  groupVIII 12108945      55 0.46666667 0.30000000
#> 1636 fam.pop    C.PAL  groupVIII 12108945      55 0.42857143 0.25000000
#> 1637 fam.pop    C.SMR  groupVIII 12108945      55 0.00000000 0.00000000
#> 1638 fam.pop    C.UPD  groupVIII 12108945      55 0.53333333 0.33333333
#> 1639 fam.pop    D.ASP  groupVIII 12108945      55 1.00000000 0.50000000
#> 1640 fam.pop    D.CLF  groupVIII 12108945      55 0.50000000 0.25000000
#> 1641 fam.pop    D.OPL  groupVIII 12108945      55 0.46666667 0.30000000
#> 1642 fam.pop    D.PAL  groupVIII 12108945      55 0.57142857 0.50000000
#> 1643 fam.pop    D.SMR  groupVIII 12108945      55 0.46666667 0.30000000
#> 1644 fam.pop    D.UPD  groupVIII 12108945      55 0.53333333 0.33333333
#> 1645     pop      ASP  groupVIII 12108945      55 0.53333333 0.50000000
#> 1646     pop      CLF  groupVIII 12108945      55 0.40909091 0.25000000
#> 1647     pop      OPL  groupVIII 12108945      55 0.45528455 0.33333333
#> 1648     pop      PAL  groupVIII 12108945      55 0.51494253 0.46666667
#> 1649     pop      SMR  groupVIII 12108945      55 0.29615385 0.17500000
#> 1650     pop      UPD  groupVIII 12108945      55 0.51692308 0.46153846
#> 1651 fam.pop    A.ASP   groupXIX 12485182      56 0.53333333 0.33333333
#> 1652 fam.pop    A.CLF   groupXIX 12485182      56 0.00000000 0.00000000
#> 1653 fam.pop    A.OPL   groupXIX 12485182      56 0.25000000 0.12500000
#> 1654 fam.pop    A.PAL   groupXIX 12485182      56 0.53571429 0.37500000
#> 1655 fam.pop    A.SMR   groupXIX 12485182      56 0.35555556 0.20000000
#> 1656 fam.pop    A.UPD   groupXIX 12485182      56 0.60000000 0.50000000
#> 1657 fam.pop    B.ASP   groupXIX 12485182      56 0.50000000 0.25000000
#> 1658 fam.pop    B.CLF   groupXIX 12485182      56 1.00000000 0.50000000
#> 1659 fam.pop    B.OPL   groupXIX 12485182      56 0.35555556 0.20000000
#> 1660 fam.pop    B.PAL   groupXIX 12485182      56 0.53571429 0.37500000
#> 1661 fam.pop    B.SMR   groupXIX 12485182      56 0.30303030 0.16666667
#> 1662 fam.pop    B.UPD   groupXIX 12485182      56 0.25000000 0.12500000
#> 1663 fam.pop    C.ASP   groupXIX 12485182      56 0.66666667 0.50000000
#> 1664 fam.pop    C.CLF   groupXIX 12485182      56 0.00000000 0.00000000
#> 1665 fam.pop    C.OPL   groupXIX 12485182      56 0.40909091 0.25000000
#> 1666 fam.pop    C.PAL   groupXIX 12485182      56 0.57142857 0.50000000
#> 1667 fam.pop    C.SMR   groupXIX 12485182      56 0.42857143 0.25000000
#> 1668 fam.pop    C.UPD   groupXIX 12485182      56 0.42857143 0.25000000
#> 1669 fam.pop    D.ASP   groupXIX 12485182      56 0.50000000 0.25000000
#> 1670 fam.pop    D.CLF   groupXIX 12485182      56 0.00000000 0.00000000
#> 1671 fam.pop    D.OPL   groupXIX 12485182      56 0.35555556 0.20000000
#> 1672 fam.pop    D.PAL   groupXIX 12485182      56 0.25000000 0.12500000
#> 1673 fam.pop    D.SMR   groupXIX 12485182      56 0.25000000 0.12500000
#> 1674 fam.pop    D.UPD   groupXIX 12485182      56 0.25000000 0.12500000
#> 1675     pop      ASP   groupXIX 12485182      56 0.47058824 0.33333333
#> 1676     pop      CLF   groupXIX 12485182      56 0.16666667 0.08333333
#> 1677     pop      OPL   groupXIX 12485182      56 0.32820513 0.20000000
#> 1678     pop      PAL   groupXIX 12485182      56 0.46572581 0.34375000
#> 1679     pop      SMR   groupXIX 12485182      56 0.30867710 0.18421053
#> 1680     pop      UPD   groupXIX 12485182      56 0.37011494 0.23333333
#> 1681 fam.pop    A.ASP    groupXI 12833733      57 0.00000000 0.00000000
#> 1682 fam.pop    A.CLF    groupXI 12833733      57 0.00000000 0.00000000
#> 1683 fam.pop    A.OPL    groupXI 12833733      57 0.00000000 0.00000000
#> 1684 fam.pop    A.PAL    groupXI 12833733      57 0.00000000 0.00000000
#> 1685 fam.pop    A.SMR    groupXI 12833733      57 0.00000000 0.00000000
#> 1686 fam.pop    A.UPD    groupXI 12833733      57 0.00000000 0.00000000
#> 1687 fam.pop    B.ASP    groupXI 12833733      57 0.00000000 0.00000000
#> 1688 fam.pop    B.CLF    groupXI 12833733      57 0.00000000 0.00000000
#> 1689 fam.pop    B.OPL    groupXI 12833733      57 0.30303030 0.16666667
#> 1690 fam.pop    B.PAL    groupXI 12833733      57 0.00000000 0.00000000
#> 1691 fam.pop    B.SMR    groupXI 12833733      57 0.00000000 0.00000000
#> 1692 fam.pop    B.UPD    groupXI 12833733      57 0.00000000 0.00000000
#> 1693 fam.pop    C.ASP    groupXI 12833733      57 0.00000000 0.00000000
#> 1694 fam.pop    C.CLF    groupXI 12833733      57 0.00000000 0.00000000
#> 1695 fam.pop    C.OPL    groupXI 12833733      57 0.20000000 0.10000000
#> 1696 fam.pop    C.PAL    groupXI 12833733      57 0.00000000 0.00000000
#> 1697 fam.pop    C.SMR    groupXI 12833733      57 0.16666667 0.08333333
#> 1698 fam.pop    C.UPD    groupXI 12833733      57 0.00000000 0.00000000
#> 1699 fam.pop    D.ASP    groupXI 12833733      57 1.00000000 0.50000000
#> 1700 fam.pop    D.CLF    groupXI 12833733      57 0.00000000 0.00000000
#> 1701 fam.pop    D.OPL    groupXI 12833733      57 0.16666667 0.08333333
#> 1702 fam.pop    D.PAL    groupXI 12833733      57 0.00000000 0.00000000
#> 1703 fam.pop    D.SMR    groupXI 12833733      57 0.00000000 0.00000000
#> 1704 fam.pop    D.UPD    groupXI 12833733      57 0.00000000 0.00000000
#> 1705     pop      ASP    groupXI 12833733      57 0.12500000 0.06250000
#> 1706     pop      CLF    groupXI 12833733      57 0.00000000 0.00000000
#> 1707     pop      OPL    groupXI 12833733      57 0.18461538 0.10000000
#> 1708     pop      PAL    groupXI 12833733      57 0.00000000 0.00000000
#> 1709     pop      SMR    groupXI 12833733      57 0.05555556 0.02777778
#> 1710     pop      UPD    groupXI 12833733      57 0.00000000 0.00000000
#> 1711 fam.pop    A.ASP groupXVIII 12937786      58 0.53333333 0.66666667
#> 1712 fam.pop    A.CLF groupXVIII 12937786      58 0.00000000 0.00000000
#> 1713 fam.pop    A.OPL groupXVIII 12937786      58 0.46666667 0.30000000
#> 1714 fam.pop    A.PAL groupXVIII 12937786      58 0.53571429 0.37500000
#> 1715 fam.pop    A.SMR groupXVIII 12937786      58 0.00000000 0.00000000
#> 1716 fam.pop    A.UPD groupXVIII 12937786      58 0.00000000 0.00000000
#> 1717 fam.pop    B.ASP groupXVIII 12937786      58 0.53333333 0.33333333
#> 1718 fam.pop    B.CLF groupXVIII 12937786      58 0.00000000 0.00000000
#> 1719 fam.pop    B.OPL groupXVIII 12937786      58 0.40909091 0.25000000
#> 1720 fam.pop    B.PAL groupXVIII 12937786      58 0.50000000 0.75000000
#> 1721 fam.pop    B.SMR groupXVIII 12937786      58 0.35555556 0.20000000
#> 1722 fam.pop    B.UPD groupXVIII 12937786      58 0.00000000 0.00000000
#> 1723 fam.pop    C.ASP groupXVIII 12937786      58 0.50000000 0.25000000
#> 1724 fam.pop    C.CLF groupXVIII 12937786      58 0.00000000 0.00000000
#> 1725 fam.pop    C.OPL groupXVIII 12937786      58 0.43956044 0.28571429
#> 1726 fam.pop    C.PAL groupXVIII 12937786      58 0.53571429 0.37500000
#> 1727 fam.pop    C.SMR groupXVIII 12937786      58 0.53333333 0.33333333
#> 1728 fam.pop    C.UPD groupXVIII 12937786      58 0.57142857 0.50000000
#> 1729 fam.pop    D.ASP groupXVIII 12937786      58 1.00000000 0.50000000
#> 1730 fam.pop    D.CLF groupXVIII 12937786      58 0.00000000 0.00000000
#> 1731 fam.pop    D.OPL groupXVIII 12937786      58 0.36263736 0.21428571
#> 1732 fam.pop    D.PAL groupXVIII 12937786      58 0.53333333 0.33333333
#> 1733 fam.pop    D.SMR groupXVIII 12937786      58 0.16666667 0.08333333
#> 1734 fam.pop    D.UPD groupXVIII 12937786      58 0.53333333 0.33333333
#> 1735     pop      ASP groupXVIII 12937786      58 0.52287582 0.44444444
#> 1736     pop      CLF groupXVIII 12937786      58 0.00000000 0.00000000
#> 1737     pop      OPL groupXVIII 12937786      58 0.39265306 0.26000000
#> 1738     pop      PAL groupXVIII 12937786      58 0.50769231 0.42307692
#> 1739     pop      SMR groupXVIII 12937786      58 0.24603175 0.13888889
#> 1740     pop      UPD groupXVIII 12937786      58 0.41558442 0.27272727
#> 1741 fam.pop    A.ASP   groupIII 13177027      59 0.53333333 0.33333333
#> 1742 fam.pop    A.CLF   groupIII 13177027      59 0.50000000 0.25000000
#> 1743 fam.pop    A.OPL   groupIII 13177027      59 0.46666667 0.30000000
#> 1744 fam.pop    A.PAL   groupIII 13177027      59 0.00000000 0.00000000
#> 1745 fam.pop    A.SMR   groupIII 13177027      59 0.53333333 0.60000000
#> 1746 fam.pop    A.UPD   groupIII 13177027      59 0.60000000 0.50000000
#> 1747 fam.pop    B.ASP   groupIII 13177027      59 0.33333333 0.83333333
#> 1748 fam.pop    B.CLF   groupIII 13177027      59 1.00000000 0.50000000
#> 1749 fam.pop    B.OPL   groupIII 13177027      59 0.53846154 0.50000000
#> 1750 fam.pop    B.PAL   groupIII 13177027      59 0.25000000 0.12500000
#> 1751 fam.pop    B.SMR   groupIII 13177027      59 0.57142857 0.50000000
#> 1752 fam.pop    B.UPD   groupIII 13177027      59 0.53571429 0.37500000
#> 1753 fam.pop    C.ASP   groupIII 13177027      59 0.33333333 0.16666667
#> 1754 fam.pop    C.CLF   groupIII 13177027      59 1.00000000 0.50000000
#> 1755 fam.pop    C.OPL   groupIII 13177027      59 0.54545455 0.50000000
#> 1756 fam.pop    C.PAL   groupIII 13177027      59 0.00000000 0.00000000
#> 1757 fam.pop    C.SMR   groupIII 13177027      59 0.54545455 0.50000000
#> 1758 fam.pop    C.UPD   groupIII 13177027      59 0.53333333 0.33333333
#> 1759 fam.pop    D.ASP   groupIII 13177027      59 0.66666667 0.50000000
#> 1760 fam.pop    D.CLF   groupIII 13177027      59 0.50000000 0.25000000
#> 1761 fam.pop    D.OPL   groupIII 13177027      59 0.40909091 0.25000000
#> 1762 fam.pop    D.PAL   groupIII 13177027      59 0.25000000 0.12500000
#> 1763 fam.pop    D.SMR   groupIII 13177027      59 0.46666667 0.30000000
#> 1764 fam.pop    D.UPD   groupIII 13177027      59 0.53333333 0.33333333
#> 1765     pop      ASP   groupIII 13177027      59 0.51948052 0.45454545
#> 1766     pop      CLF   groupIII 13177027      59 0.48484848 0.33333333
#> 1767     pop      OPL   groupIII 13177027      59 0.48847518 0.39583333
#> 1768     pop      PAL   groupIII 13177027      59 0.12096774 0.06250000
#> 1769     pop      SMR   groupIII 13177027      59 0.51153846 0.47500000
#> 1770     pop      UPD   groupIII 13177027      59 0.49230769 0.38461538
#> 1771 fam.pop    A.ASP    groupVI 13335577      60 0.50000000 0.25000000
#> 1772 fam.pop    A.CLF    groupVI 13335577      60 0.66666667 0.50000000
#> 1773 fam.pop    A.OPL    groupVI 13335577      60 0.35555556 0.20000000
#> 1774 fam.pop    A.PAL    groupVI 13335577      60 0.00000000 0.00000000
#> 1775 fam.pop    A.SMR    groupVI 13335577      60 0.48484848 0.33333333
#> 1776 fam.pop    A.UPD    groupVI 13335577      60 0.53333333 0.33333333
#> 1777 fam.pop    B.ASP    groupVI 13335577      60 0.60000000 0.50000000
#> 1778 fam.pop    B.CLF    groupVI 13335577      60 0.00000000 1.00000000
#> 1779 fam.pop    B.OPL    groupVI 13335577      60 0.53030303 0.41666667
#> 1780 fam.pop    B.PAL    groupVI 13335577      60 0.00000000 0.00000000
#> 1781 fam.pop    B.SMR    groupVI 13335577      60 0.35555556 0.20000000
#> 1782 fam.pop    B.UPD    groupVI 13335577      60 0.53571429 0.37500000
#> 1783 fam.pop    C.ASP    groupVI 13335577      60 0.66666667 0.50000000
#> 1784 fam.pop    C.CLF    groupVI 13335577      60 0.00000000 1.00000000
#> 1785 fam.pop    C.OPL    groupVI 13335577      60 0.48484848 0.33333333
#> 1786 fam.pop    C.PAL    groupVI 13335577      60 0.00000000 0.00000000
#> 1787 fam.pop    C.SMR    groupVI 13335577      60 0.57142857 0.50000000
#> 1788 fam.pop    C.UPD    groupVI 13335577      60 0.25000000 0.12500000
#> 1789 fam.pop    D.ASP    groupVI 13335577      60        NaN         NA
#> 1790 fam.pop    D.CLF    groupVI 13335577      60 0.50000000 0.75000000
#> 1791 fam.pop    D.OPL    groupVI 13335577      60 0.53030303 0.41666667
#> 1792 fam.pop    D.PAL    groupVI 13335577      60 0.00000000 0.00000000
#> 1793 fam.pop    D.SMR    groupVI 13335577      60 0.48484848 0.33333333
#> 1794 fam.pop    D.UPD    groupVI 13335577      60 0.60000000 0.50000000
#> 1795     pop      ASP    groupVI 13335577      60 0.52747253 0.42857143
#> 1796     pop      CLF    groupVI 13335577      60 0.40909091 0.75000000
#> 1797     pop      OPL    groupVI 13335577      60 0.46376812 0.34782609
#> 1798     pop      PAL    groupVI 13335577      60 0.00000000 0.00000000
#> 1799     pop      SMR    groupVI 13335577      60 0.45528455 0.33333333
#> 1800     pop      UPD    groupVI 13335577      60 0.45238095 0.32142857
#> 1801 fam.pop    A.ASP     groupX 13718292      61 0.00000000 0.00000000
#> 1802 fam.pop    A.CLF     groupX 13718292      61        NaN         NA
#> 1803 fam.pop    A.OPL     groupX 13718292      61 0.00000000 0.00000000
#> 1804 fam.pop    A.PAL     groupX 13718292      61 0.00000000 0.00000000
#> 1805 fam.pop    A.SMR     groupX 13718292      61 0.00000000 0.00000000
#> 1806 fam.pop    A.UPD     groupX 13718292      61 0.00000000 0.00000000
#> 1807 fam.pop    B.ASP     groupX 13718292      61 0.00000000 0.00000000
#> 1808 fam.pop    B.CLF     groupX 13718292      61        NaN         NA
#> 1809 fam.pop    B.OPL     groupX 13718292      61 0.00000000 0.00000000
#> 1810 fam.pop    B.PAL     groupX 13718292      61 0.00000000 0.00000000
#> 1811 fam.pop    B.SMR     groupX 13718292      61 0.00000000 0.00000000
#> 1812 fam.pop    B.UPD     groupX 13718292      61 0.25000000 0.12500000
#> 1813 fam.pop    C.ASP     groupX 13718292      61 0.00000000 0.00000000
#> 1814 fam.pop    C.CLF     groupX 13718292      61 1.00000000 0.50000000
#> 1815 fam.pop    C.OPL     groupX 13718292      61 0.50000000 0.25000000
#> 1816 fam.pop    C.PAL     groupX 13718292      61 0.00000000 0.00000000
#> 1817 fam.pop    C.SMR     groupX 13718292      61 0.35555556 0.20000000
#> 1818 fam.pop    C.UPD     groupX 13718292      61 0.50000000 0.25000000
#> 1819 fam.pop    D.ASP     groupX 13718292      61        NaN         NA
#> 1820 fam.pop    D.CLF     groupX 13718292      61 0.00000000 0.00000000
#> 1821 fam.pop    D.OPL     groupX 13718292      61 0.00000000 0.00000000
#> 1822 fam.pop    D.PAL     groupX 13718292      61 0.00000000 0.00000000
#> 1823 fam.pop    D.SMR     groupX 13718292      61 0.00000000 0.00000000
#> 1824 fam.pop    D.UPD     groupX 13718292      61 0.00000000 0.00000000
#> 1825     pop      ASP     groupX 13718292      61 0.00000000 0.00000000
#> 1826     pop      CLF     groupX 13718292      61 0.33333333 0.16666667
#> 1827     pop      OPL     groupX 13718292      61 0.05882353 0.02941176
#> 1828     pop      PAL     groupX 13718292      61 0.00000000 0.00000000
#> 1829     pop      SMR     groupX 13718292      61 0.09743590 0.05000000
#> 1830     pop      UPD     groupX 13718292      61 0.15942029 0.08333333
#> 1831 fam.pop    A.ASP    groupVI 13895716      62 0.53333333 0.33333333
#> 1832 fam.pop    A.CLF    groupVI 13895716      62 0.00000000 0.00000000
#> 1833 fam.pop    A.OPL    groupVI 13895716      62 0.40909091 0.25000000
#> 1834 fam.pop    A.PAL    groupVI 13895716      62 0.00000000 0.00000000
#> 1835 fam.pop    A.SMR    groupVI 13895716      62 0.16666667 0.08333333
#> 1836 fam.pop    A.UPD    groupVI 13895716      62 0.33333333 0.16666667
#> 1837 fam.pop    B.ASP    groupVI 13895716      62 0.00000000 0.00000000
#> 1838 fam.pop    B.CLF    groupVI 13895716      62 0.00000000 0.00000000
#> 1839 fam.pop    B.OPL    groupVI 13895716      62 0.36263736 0.21428571
#> 1840 fam.pop    B.PAL    groupVI 13895716      62 0.00000000 0.00000000
#> 1841 fam.pop    B.SMR    groupVI 13895716      62 0.16666667 0.08333333
#> 1842 fam.pop    B.UPD    groupVI 13895716      62 0.25000000 0.12500000
#> 1843 fam.pop    C.ASP    groupVI 13895716      62 0.00000000 0.00000000
#> 1844 fam.pop    C.CLF    groupVI 13895716      62 0.00000000 0.00000000
#> 1845 fam.pop    C.OPL    groupVI 13895716      62 0.16666667 0.08333333
#> 1846 fam.pop    C.PAL    groupVI 13895716      62 0.00000000 0.00000000
#> 1847 fam.pop    C.SMR    groupVI 13895716      62 0.00000000 0.00000000
#> 1848 fam.pop    C.UPD    groupVI 13895716      62 0.00000000 0.00000000
#> 1849 fam.pop    D.ASP    groupVI 13895716      62 0.00000000 0.00000000
#> 1850 fam.pop    D.CLF    groupVI 13895716      62 0.00000000 0.00000000
#> 1851 fam.pop    D.OPL    groupVI 13895716      62 0.43956044 0.28571429
#> 1852 fam.pop    D.PAL    groupVI 13895716      62 0.00000000 0.00000000
#> 1853 fam.pop    D.SMR    groupVI 13895716      62 0.16666667 0.08333333
#> 1854 fam.pop    D.UPD    groupVI 13895716      62 0.25000000 0.12500000
#> 1855     pop      ASP    groupVI 13895716      62 0.18947368 0.10000000
#> 1856     pop      CLF    groupVI 13895716      62 0.00000000 0.00000000
#> 1857     pop      OPL    groupVI 13895716      62 0.34012066 0.21153846
#> 1858     pop      PAL    groupVI 13895716      62 0.00000000 0.00000000
#> 1859     pop      SMR    groupVI 13895716      62 0.11968085 0.06250000
#> 1860     pop      UPD    groupVI 13895716      62 0.19841270 0.10714286
#> 1861 fam.pop    A.ASP   groupXIV 14009990      63 0.60000000 0.50000000
#> 1862 fam.pop    A.CLF   groupXIV 14009990      63 0.50000000 0.25000000
#> 1863 fam.pop    A.OPL   groupXIV 14009990      63 0.42857143 0.25000000
#> 1864 fam.pop    A.PAL   groupXIV 14009990      63 0.50000000 0.75000000
#> 1865 fam.pop    A.SMR   groupXIV 14009990      63 0.46666667 0.30000000
#> 1866 fam.pop    A.UPD   groupXIV 14009990      63 0.33333333 0.16666667
#> 1867 fam.pop    B.ASP   groupXIV 14009990      63 0.50000000 0.25000000
#> 1868 fam.pop    B.CLF   groupXIV 14009990      63 1.00000000 0.50000000
#> 1869 fam.pop    B.OPL   groupXIV 14009990      63 0.53571429 0.37500000
#> 1870 fam.pop    B.PAL   groupXIV 14009990      63 0.50000000 0.75000000
#> 1871 fam.pop    B.SMR   groupXIV 14009990      63 0.42857143 0.25000000
#> 1872 fam.pop    B.UPD   groupXIV 14009990      63 0.53571429 0.37500000
#> 1873 fam.pop    C.ASP   groupXIV 14009990      63 0.50000000 0.25000000
#> 1874 fam.pop    C.CLF   groupXIV 14009990      63 0.00000000 0.00000000
#> 1875 fam.pop    C.OPL   groupXIV 14009990      63 0.30303030 0.16666667
#> 1876 fam.pop    C.PAL   groupXIV 14009990      63 0.33333333 0.83333333
#> 1877 fam.pop    C.SMR   groupXIV 14009990      63 0.57142857 0.50000000
#> 1878 fam.pop    C.UPD   groupXIV 14009990      63 0.50000000 0.25000000
#> 1879 fam.pop    D.ASP   groupXIV 14009990      63 1.00000000 0.50000000
#> 1880 fam.pop    D.CLF   groupXIV 14009990      63 0.00000000 1.00000000
#> 1881 fam.pop    D.OPL   groupXIV 14009990      63 0.25000000 0.12500000
#> 1882 fam.pop    D.PAL   groupXIV 14009990      63 0.33333333 0.83333333
#> 1883 fam.pop    D.SMR   groupXIV 14009990      63 0.33333333 0.16666667
#> 1884 fam.pop    D.UPD   groupXIV 14009990      63 0.42857143 0.25000000
#> 1885     pop      ASP   groupXIV 14009990      63 0.50000000 0.37500000
#> 1886     pop      CLF   groupXIV 14009990      63 0.53333333 0.40000000
#> 1887     pop      OPL   groupXIV 14009990      63 0.35555556 0.22222222
#> 1888     pop      PAL   groupXIV 14009990      63 0.33684211 0.80000000
#> 1889     pop      SMR   groupXIV 14009990      63 0.44354839 0.31250000
#> 1890     pop      UPD   groupXIV 14009990      63 0.40923077 0.26923077
#> 1891 fam.pop    A.ASP   groupXII 14158350      64 0.33333333 0.16666667
#> 1892 fam.pop    A.CLF   groupXII 14158350      64 0.66666667 0.50000000
#> 1893 fam.pop    A.OPL   groupXII 14158350      64 0.42857143 0.25000000
#> 1894 fam.pop    A.PAL   groupXII 14158350      64 0.46666667 0.70000000
#> 1895 fam.pop    A.SMR   groupXII 14158350      64 0.54545455 0.50000000
#> 1896 fam.pop    A.UPD   groupXII 14158350      64 0.53333333 0.33333333
#> 1897 fam.pop    B.ASP   groupXII 14158350      64 0.53333333 0.33333333
#> 1898 fam.pop    B.CLF   groupXII 14158350      64 0.00000000 0.00000000
#> 1899 fam.pop    B.OPL   groupXII 14158350      64 0.48484848 0.33333333
#> 1900 fam.pop    B.PAL   groupXII 14158350      64 0.53571429 0.37500000
#> 1901 fam.pop    B.SMR   groupXII 14158350      64 0.46666667 0.30000000
#> 1902 fam.pop    B.UPD   groupXII 14158350      64 0.53571429 0.62500000
#> 1903 fam.pop    C.ASP   groupXII 14158350      64 0.60000000 0.50000000
#> 1904 fam.pop    C.CLF   groupXII 14158350      64 0.00000000 0.00000000
#> 1905 fam.pop    C.OPL   groupXII 14158350      64 0.30303030 0.16666667
#> 1906 fam.pop    C.PAL   groupXII 14158350      64 0.00000000 1.00000000
#> 1907 fam.pop    C.SMR   groupXII 14158350      64 0.48484848 0.33333333
#> 1908 fam.pop    C.UPD   groupXII 14158350      64 0.42857143 0.25000000
#> 1909 fam.pop    D.ASP   groupXII 14158350      64 0.50000000 0.25000000
#> 1910 fam.pop    D.CLF   groupXII 14158350      64 0.00000000 0.00000000
#> 1911 fam.pop    D.OPL   groupXII 14158350      64 0.55555556 0.50000000
#> 1912 fam.pop    D.PAL   groupXII 14158350      64 0.53333333 0.33333333
#> 1913 fam.pop    D.SMR   groupXII 14158350      64 0.46666667 0.30000000
#> 1914 fam.pop    D.UPD   groupXII 14158350      64 0.42857143 0.25000000
#> 1915     pop      ASP   groupXII 14158350      64 0.45454545 0.31818182
#> 1916     pop      CLF   groupXII 14158350      64 0.30303030 0.16666667
#> 1917     pop      OPL   groupXII 14158350      64 0.43786295 0.30952381
#> 1918     pop      PAL   groupXII 14158350      64 0.50793651 0.57142857
#> 1919     pop      SMR   groupXII 14158350      64 0.47357294 0.36363636
#> 1920     pop      UPD   groupXII 14158350      64 0.48045977 0.36666667
#> 1921 fam.pop    A.ASP    groupIX 14164215      65 0.33333333 0.16666667
#> 1922 fam.pop    A.CLF    groupIX 14164215      65 0.00000000 0.00000000
#> 1923 fam.pop    A.OPL    groupIX 14164215      65 0.57142857 0.50000000
#> 1924 fam.pop    A.PAL    groupIX 14164215      65 0.46666667 0.30000000
#> 1925 fam.pop    A.SMR    groupIX 14164215      65 0.53030303 0.41666667
#> 1926 fam.pop    A.UPD    groupIX 14164215      65 0.50000000 0.25000000
#> 1927 fam.pop    B.ASP    groupIX 14164215      65 0.33333333 0.16666667
#> 1928 fam.pop    B.CLF    groupIX 14164215      65 1.00000000 0.50000000
#> 1929 fam.pop    B.OPL    groupIX 14164215      65 0.53333333 0.60000000
#> 1930 fam.pop    B.PAL    groupIX 14164215      65 0.53571429 0.37500000
#> 1931 fam.pop    B.SMR    groupIX 14164215      65 0.25000000 0.12500000
#> 1932 fam.pop    B.UPD    groupIX 14164215      65 0.57142857 0.50000000
#> 1933 fam.pop    C.ASP    groupIX 14164215      65 0.53333333 0.33333333
#> 1934 fam.pop    C.CLF    groupIX 14164215      65 1.00000000 0.50000000
#> 1935 fam.pop    C.OPL    groupIX 14164215      65 0.20000000 0.90000000
#> 1936 fam.pop    C.PAL    groupIX 14164215      65 0.42857143 0.25000000
#> 1937 fam.pop    C.SMR    groupIX 14164215      65 0.48484848 0.33333333
#> 1938 fam.pop    C.UPD    groupIX 14164215      65 0.25000000 0.12500000
#> 1939 fam.pop    D.ASP    groupIX 14164215      65 0.00000000 0.00000000
#> 1940 fam.pop    D.CLF    groupIX 14164215      65 0.50000000 0.75000000
#> 1941 fam.pop    D.OPL    groupIX 14164215      65 0.55555556 0.50000000
#> 1942 fam.pop    D.PAL    groupIX 14164215      65 0.25000000 0.12500000
#> 1943 fam.pop    D.SMR    groupIX 14164215      65 0.53030303 0.41666667
#> 1944 fam.pop    D.UPD    groupIX 14164215      65 0.53571429 0.62500000
#> 1945     pop      ASP    groupIX 14164215      65 0.31168831 0.18181818
#> 1946     pop      CLF    groupIX 14164215      65 0.53030303 0.41666667
#> 1947     pop      OPL    groupIX 14164215      65 0.47795164 0.63157895
#> 1948     pop      PAL    groupIX 14164215      65 0.40106952 0.26470588
#> 1949     pop      SMR    groupIX 14164215      65 0.45983087 0.34090909
#> 1950     pop      UPD    groupIX 14164215      65 0.49470899 0.39285714
#> 1951 fam.pop    A.ASP    groupVI 14236496      66 0.53333333 0.33333333
#> 1952 fam.pop    A.CLF    groupVI 14236496      66 0.50000000 0.25000000
#> 1953 fam.pop    A.OPL    groupVI 14236496      66 0.46666667 0.30000000
#> 1954 fam.pop    A.PAL    groupVI 14236496      66 0.53571429 0.62500000
#> 1955 fam.pop    A.SMR    groupVI 14236496      66 0.57142857 0.50000000
#> 1956 fam.pop    A.UPD    groupVI 14236496      66 0.33333333 0.16666667
#> 1957 fam.pop    B.ASP    groupVI 14236496      66 0.60000000 0.50000000
#> 1958 fam.pop    B.CLF    groupVI 14236496      66 0.00000000 1.00000000
#> 1959 fam.pop    B.OPL    groupVI 14236496      66 0.46666667 0.30000000
#> 1960 fam.pop    B.PAL    groupVI 14236496      66 0.42857143 0.75000000
#> 1961 fam.pop    B.SMR    groupVI 14236496      66 0.33333333 0.16666667
#> 1962 fam.pop    B.UPD    groupVI 14236496      66 0.57142857 0.50000000
#> 1963 fam.pop    C.ASP    groupVI 14236496      66 0.50000000 0.75000000
#> 1964 fam.pop    C.CLF    groupVI 14236496      66 0.00000000 0.00000000
#> 1965 fam.pop    C.OPL    groupVI 14236496      66 0.53333333 0.40000000
#> 1966 fam.pop    C.PAL    groupVI 14236496      66 0.53571429 0.62500000
#> 1967 fam.pop    C.SMR    groupVI 14236496      66 0.46666667 0.30000000
#> 1968 fam.pop    C.UPD    groupVI 14236496      66 0.50000000 0.75000000
#> 1969 fam.pop    D.ASP    groupVI 14236496      66 1.00000000 0.50000000
#> 1970 fam.pop    D.CLF    groupVI 14236496      66 0.66666667 0.50000000
#> 1971 fam.pop    D.OPL    groupVI 14236496      66 0.53333333 0.40000000
#> 1972 fam.pop    D.PAL    groupVI 14236496      66 0.53571429 0.62500000
#> 1973 fam.pop    D.SMR    groupVI 14236496      66 0.55555556 0.50000000
#> 1974 fam.pop    D.UPD    groupVI 14236496      66 0.53571429 0.37500000
#> 1975     pop      ASP    groupVI 14236496      66 0.52941176 0.50000000
#> 1976     pop      CLF    groupVI 14236496      66 0.53030303 0.41666667
#> 1977     pop      OPL    groupVI 14236496      66 0.46666667 0.35000000
#> 1978     pop      PAL    groupVI 14236496      66 0.46572581 0.65625000
#> 1979     pop      SMR    groupVI 14236496      66 0.48663102 0.38235294
#> 1980     pop      UPD    groupVI 14236496      66 0.50769231 0.42307692
#> 1981 fam.pop    A.ASP     groupX 14273796      67 0.00000000 0.00000000
#> 1982 fam.pop    A.CLF     groupX 14273796      67 0.00000000 0.00000000
#> 1983 fam.pop    A.OPL     groupX 14273796      67 0.00000000 0.00000000
#> 1984 fam.pop    A.PAL     groupX 14273796      67 0.35555556 0.20000000
#> 1985 fam.pop    A.SMR     groupX 14273796      67 0.00000000 0.00000000
#> 1986 fam.pop    A.UPD     groupX 14273796      67 0.00000000 0.00000000
#> 1987 fam.pop    B.ASP     groupX 14273796      67 0.00000000 0.00000000
#> 1988 fam.pop    B.CLF     groupX 14273796      67 0.00000000 0.00000000
#> 1989 fam.pop    B.OPL     groupX 14273796      67 0.00000000 0.00000000
#> 1990 fam.pop    B.PAL     groupX 14273796      67 0.57142857 0.50000000
#> 1991 fam.pop    B.SMR     groupX 14273796      67 0.00000000 0.00000000
#> 1992 fam.pop    B.UPD     groupX 14273796      67 0.00000000 0.00000000
#> 1993 fam.pop    C.ASP     groupX 14273796      67 0.00000000 0.00000000
#> 1994 fam.pop    C.CLF     groupX 14273796      67 0.00000000 0.00000000
#> 1995 fam.pop    C.OPL     groupX 14273796      67 0.00000000 0.00000000
#> 1996 fam.pop    C.PAL     groupX 14273796      67 0.00000000 0.00000000
#> 1997 fam.pop    C.SMR     groupX 14273796      67 0.00000000 0.00000000
#> 1998 fam.pop    C.UPD     groupX 14273796      67 0.00000000 0.00000000
#> 1999 fam.pop    D.ASP     groupX 14273796      67 0.00000000 0.00000000
#> 2000 fam.pop    D.CLF     groupX 14273796      67 0.00000000 0.00000000
#> 2001 fam.pop    D.OPL     groupX 14273796      67 0.00000000 0.00000000
#> 2002 fam.pop    D.PAL     groupX 14273796      67 0.66666667 0.50000000
#> 2003 fam.pop    D.SMR     groupX 14273796      67 0.00000000 0.00000000
#> 2004 fam.pop    D.UPD     groupX 14273796      67 0.25000000 0.12500000
#> 2005     pop      ASP     groupX 14273796      67 0.00000000 0.00000000
#> 2006     pop      CLF     groupX 14273796      67 0.00000000 0.00000000
#> 2007     pop      OPL     groupX 14273796      67 0.00000000 0.00000000
#> 2008     pop      PAL     groupX 14273796      67 0.40459770 0.26666667
#> 2009     pop      SMR     groupX 14273796      67 0.00000000 0.00000000
#> 2010     pop      UPD     groupX 14273796      67 0.07692308 0.03846154
#> 2011 fam.pop    A.ASP  groupVIII 14445389      68 0.00000000 0.00000000
#> 2012 fam.pop    A.CLF  groupVIII 14445389      68 0.50000000 0.25000000
#> 2013 fam.pop    A.OPL  groupVIII 14445389      68 0.53030303 0.41666667
#> 2014 fam.pop    A.PAL  groupVIII 14445389      68 0.20000000 0.90000000
#> 2015 fam.pop    A.SMR  groupVIII 14445389      68 0.40909091 0.25000000
#> 2016 fam.pop    A.UPD  groupVIII 14445389      68 0.53333333 0.33333333
#> 2017 fam.pop    B.ASP  groupVIII 14445389      68 0.50000000 0.25000000
#> 2018 fam.pop    B.CLF  groupVIII 14445389      68 0.00000000 0.00000000
#> 2019 fam.pop    B.OPL  groupVIII 14445389      68 0.20000000 0.10000000
#> 2020 fam.pop    B.PAL  groupVIII 14445389      68 0.25000000 0.87500000
#> 2021 fam.pop    B.SMR  groupVIII 14445389      68 0.60000000 0.50000000
#> 2022 fam.pop    B.UPD  groupVIII 14445389      68 0.42857143 0.25000000
#> 2023 fam.pop    C.ASP  groupVIII 14445389      68 0.53333333 0.33333333
#> 2024 fam.pop    C.CLF  groupVIII 14445389      68 1.00000000 0.50000000
#> 2025 fam.pop    C.OPL  groupVIII 14445389      68 0.43956044 0.28571429
#> 2026 fam.pop    C.PAL  groupVIII 14445389      68 0.66666667 0.50000000
#> 2027 fam.pop    C.SMR  groupVIII 14445389      68 0.46666667 0.30000000
#> 2028 fam.pop    C.UPD  groupVIII 14445389      68 0.57142857 0.50000000
#> 2029 fam.pop    D.ASP  groupVIII 14445389      68 0.66666667 0.50000000
#> 2030 fam.pop    D.CLF  groupVIII 14445389      68 0.50000000 0.25000000
#> 2031 fam.pop    D.OPL  groupVIII 14445389      68 0.16666667 0.08333333
#> 2032 fam.pop    D.PAL  groupVIII 14445389      68 0.50000000 0.25000000
#> 2033 fam.pop    D.SMR  groupVIII 14445389      68 0.25000000 0.12500000
#> 2034 fam.pop    D.UPD  groupVIII 14445389      68 0.53333333 0.33333333
#> 2035     pop      ASP  groupVIII 14445389      68 0.39473684 0.25000000
#> 2036     pop      CLF  groupVIII 14445389      68 0.40909091 0.25000000
#> 2037     pop      OPL  groupVIII 14445389      68 0.36081560 0.22916667
#> 2038     pop      PAL  groupVIII 14445389      68 0.40923077 0.73076923
#> 2039     pop      SMR  groupVIII 14445389      68 0.41269841 0.27777778
#> 2040     pop      UPD  groupVIII 14445389      68 0.47619048 0.35714286
#> 2041 fam.pop    A.ASP    groupIX 14868533      69 0.00000000 0.00000000
#> 2042 fam.pop    A.CLF    groupIX 14868533      69 0.00000000 0.00000000
#> 2043 fam.pop    A.OPL    groupIX 14868533      69 0.20000000 0.10000000
#> 2044 fam.pop    A.PAL    groupIX 14868533      69 0.00000000 0.00000000
#> 2045 fam.pop    A.SMR    groupIX 14868533      69 0.35555556 0.20000000
#> 2046 fam.pop    A.UPD    groupIX 14868533      69 0.00000000 0.00000000
#> 2047 fam.pop    B.ASP    groupIX 14868533      69 0.00000000 0.00000000
#> 2048 fam.pop    B.CLF    groupIX 14868533      69 0.00000000 0.00000000
#> 2049 fam.pop    B.OPL    groupIX 14868533      69 0.52747253 0.42857143
#> 2050 fam.pop    B.PAL    groupIX 14868533      69 0.00000000 0.00000000
#> 2051 fam.pop    B.SMR    groupIX 14868533      69 0.00000000 0.00000000
#> 2052 fam.pop    B.UPD    groupIX 14868533      69 0.00000000 0.00000000
#> 2053 fam.pop    C.ASP    groupIX 14868533      69 0.00000000 0.00000000
#> 2054 fam.pop    C.CLF    groupIX 14868533      69 0.00000000 0.00000000
#> 2055 fam.pop    C.OPL    groupIX 14868533      69 0.54545455 0.50000000
#> 2056 fam.pop    C.PAL    groupIX 14868533      69 0.00000000 0.00000000
#> 2057 fam.pop    C.SMR    groupIX 14868533      69 0.16666667 0.08333333
#> 2058 fam.pop    C.UPD    groupIX 14868533      69 0.00000000 0.00000000
#> 2059 fam.pop    D.ASP    groupIX 14868533      69 0.00000000 0.00000000
#> 2060 fam.pop    D.CLF    groupIX 14868533      69 1.00000000 0.50000000
#> 2061 fam.pop    D.OPL    groupIX 14868533      69 0.26373626 0.14285714
#> 2062 fam.pop    D.PAL    groupIX 14868533      69 0.00000000 0.00000000
#> 2063 fam.pop    D.SMR    groupIX 14868533      69 0.25000000 0.12500000
#> 2064 fam.pop    D.UPD    groupIX 14868533      69 0.00000000 0.00000000
#> 2065     pop      ASP    groupIX 14868533      69 0.00000000 0.00000000
#> 2066     pop      CLF    groupIX 14868533      69 0.20000000 0.10000000
#> 2067     pop      OPL    groupIX 14868533      69 0.42857143 0.30000000
#> 2068     pop      PAL    groupIX 14868533      69 0.00000000 0.00000000
#> 2069     pop      SMR    groupIX 14868533      69 0.17653891 0.09523810
#> 2070     pop      UPD    groupIX 14868533      69 0.00000000 0.00000000
#> 2071 fam.pop    A.ASP    groupIX 15190693      70 0.53333333 0.33333333
#> 2072 fam.pop    A.CLF    groupIX 15190693      70 0.66666667 0.50000000
#> 2073 fam.pop    A.OPL    groupIX 15190693      70 0.20000000 0.10000000
#> 2074 fam.pop    A.PAL    groupIX 15190693      70 0.00000000 0.00000000
#> 2075 fam.pop    A.SMR    groupIX 15190693      70 0.53333333 0.40000000
#> 2076 fam.pop    A.UPD    groupIX 15190693      70 0.66666667 0.50000000
#> 2077 fam.pop    B.ASP    groupIX 15190693      70 0.50000000 0.25000000
#> 2078 fam.pop    B.CLF    groupIX 15190693      70 1.00000000 0.50000000
#> 2079 fam.pop    B.OPL    groupIX 15190693      70 0.48484848 0.33333333
#> 2080 fam.pop    B.PAL    groupIX 15190693      70 0.00000000 0.00000000
#> 2081 fam.pop    B.SMR    groupIX 15190693      70 0.46666667 0.30000000
#> 2082 fam.pop    B.UPD    groupIX 15190693      70 0.00000000 0.00000000
#> 2083 fam.pop    C.ASP    groupIX 15190693      70 0.53333333 0.33333333
#> 2084 fam.pop    C.CLF    groupIX 15190693      70        NaN         NA
#> 2085 fam.pop    C.OPL    groupIX 15190693      70 0.46666667 0.30000000
#> 2086 fam.pop    C.PAL    groupIX 15190693      70 0.00000000 0.00000000
#> 2087 fam.pop    C.SMR    groupIX 15190693      70 0.30303030 0.16666667
#> 2088 fam.pop    C.UPD    groupIX 15190693      70 0.42857143 0.25000000
#> 2089 fam.pop    D.ASP    groupIX 15190693      70 0.66666667 0.50000000
#> 2090 fam.pop    D.CLF    groupIX 15190693      70 0.50000000 0.75000000
#> 2091 fam.pop    D.OPL    groupIX 15190693      70 0.46666667 0.30000000
#> 2092 fam.pop    D.PAL    groupIX 15190693      70 0.00000000 0.00000000
#> 2093 fam.pop    D.SMR    groupIX 15190693      70 0.30303030 0.16666667
#> 2094 fam.pop    D.UPD    groupIX 15190693      70 0.60000000 0.50000000
#> 2095     pop      ASP    groupIX 15190693      70 0.47894737 0.35000000
#> 2096     pop      CLF    groupIX 15190693      70 0.53333333 0.60000000
#> 2097     pop      OPL    groupIX 15190693      70 0.39605110 0.26190476
#> 2098     pop      PAL    groupIX 15190693      70 0.00000000 0.00000000
#> 2099     pop      SMR    groupIX 15190693      70 0.38372093 0.25000000
#> 2100     pop      UPD    groupIX 15190693      70 0.43115942 0.29166667
#> 2101 fam.pop    A.ASP   groupXII 15306362      71 0.00000000 0.00000000
#> 2102 fam.pop    A.CLF   groupXII 15306362      71 0.00000000 0.00000000
#> 2103 fam.pop    A.OPL   groupXII 15306362      71 0.00000000 0.00000000
#> 2104 fam.pop    A.PAL   groupXII 15306362      71 0.42857143 0.25000000
#> 2105 fam.pop    A.SMR   groupXII 15306362      71 0.16666667 0.08333333
#> 2106 fam.pop    A.UPD   groupXII 15306362      71 0.00000000 0.00000000
#> 2107 fam.pop    B.ASP   groupXII 15306362      71 0.00000000 0.00000000
#> 2108 fam.pop    B.CLF   groupXII 15306362      71 0.00000000 0.00000000
#> 2109 fam.pop    B.OPL   groupXII 15306362      71 0.00000000 0.00000000
#> 2110 fam.pop    B.PAL   groupXII 15306362      71 0.25000000 0.12500000
#> 2111 fam.pop    B.SMR   groupXII 15306362      71 0.16666667 0.08333333
#> 2112 fam.pop    B.UPD   groupXII 15306362      71 0.42857143 0.25000000
#> 2113 fam.pop    C.ASP   groupXII 15306362      71 0.00000000 0.00000000
#> 2114 fam.pop    C.CLF   groupXII 15306362      71 0.00000000 0.00000000
#> 2115 fam.pop    C.OPL   groupXII 15306362      71 0.14285714 0.07142857
#> 2116 fam.pop    C.PAL   groupXII 15306362      71 0.00000000 0.00000000
#> 2117 fam.pop    C.SMR   groupXII 15306362      71 0.16666667 0.08333333
#> 2118 fam.pop    C.UPD   groupXII 15306362      71 0.00000000 0.00000000
#> 2119 fam.pop    D.ASP   groupXII 15306362      71 0.50000000 0.25000000
#> 2120 fam.pop    D.CLF   groupXII 15306362      71 0.00000000 0.00000000
#> 2121 fam.pop    D.OPL   groupXII 15306362      71 0.00000000 0.00000000
#> 2122 fam.pop    D.PAL   groupXII 15306362      71 0.53571429 0.37500000
#> 2123 fam.pop    D.SMR   groupXII 15306362      71 0.40909091 0.25000000
#> 2124 fam.pop    D.UPD   groupXII 15306362      71 0.00000000 0.00000000
#> 2125     pop      ASP   groupXII 15306362      71 0.09090909 0.04545455
#> 2126     pop      CLF   groupXII 15306362      71 0.00000000 0.00000000
#> 2127     pop      OPL   groupXII 15306362      71 0.03703704 0.01851852
#> 2128     pop      PAL   groupXII 15306362      71 0.31451613 0.18750000
#> 2129     pop      SMR   groupXII 15306362      71 0.22340426 0.12500000
#> 2130     pop      UPD   groupXII 15306362      71 0.12873563 0.06666667
#> 2131 fam.pop    A.ASP    groupXX 15460430      72 0.60000000 0.50000000
#> 2132 fam.pop    A.CLF    groupXX 15460430      72 0.50000000 0.25000000
#> 2133 fam.pop    A.OPL    groupXX 15460430      72 0.40909091 0.25000000
#> 2134 fam.pop    A.PAL    groupXX 15460430      72 0.46666667 0.70000000
#> 2135 fam.pop    A.SMR    groupXX 15460430      72 0.40909091 0.25000000
#> 2136 fam.pop    A.UPD    groupXX 15460430      72 0.50000000 0.75000000
#> 2137 fam.pop    B.ASP    groupXX 15460430      72 0.00000000 1.00000000
#> 2138 fam.pop    B.CLF    groupXX 15460430      72 1.00000000 0.50000000
#> 2139 fam.pop    B.OPL    groupXX 15460430      72 0.53030303 0.41666667
#> 2140 fam.pop    B.PAL    groupXX 15460430      72 0.53571429 0.62500000
#> 2141 fam.pop    B.SMR    groupXX 15460430      72 0.40909091 0.25000000
#> 2142 fam.pop    B.UPD    groupXX 15460430      72 0.33333333 0.83333333
#> 2143 fam.pop    C.ASP    groupXX 15460430      72 0.60000000 0.50000000
#> 2144 fam.pop    C.CLF    groupXX 15460430      72 0.00000000 1.00000000
#> 2145 fam.pop    C.OPL    groupXX 15460430      72 0.53333333 0.40000000
#> 2146 fam.pop    C.PAL    groupXX 15460430      72 0.53333333 0.33333333
#> 2147 fam.pop    C.SMR    groupXX 15460430      72 0.25000000 0.12500000
#> 2148 fam.pop    C.UPD    groupXX 15460430      72 0.53571429 0.62500000
#> 2149 fam.pop    D.ASP    groupXX 15460430      72 0.50000000 0.75000000
#> 2150 fam.pop    D.CLF    groupXX 15460430      72 0.50000000 0.75000000
#> 2151 fam.pop    D.OPL    groupXX 15460430      72 0.36263736 0.21428571
#> 2152 fam.pop    D.PAL    groupXX 15460430      72 0.33333333 0.83333333
#> 2153 fam.pop    D.SMR    groupXX 15460430      72 0.40909091 0.25000000
#> 2154 fam.pop    D.UPD    groupXX 15460430      72 0.42857143 0.75000000
#> 2155     pop      ASP    groupXX 15460430      72 0.47894737 0.65000000
#> 2156     pop      CLF    groupXX 15460430      72 0.53030303 0.58333333
#> 2157     pop      OPL    groupXX 15460430      72 0.43882979 0.31250000
#> 2158     pop      PAL    groupXX 15460430      72 0.48045977 0.63333333
#> 2159     pop      SMR    groupXX 15460430      72 0.35940803 0.22727273
#> 2160     pop      UPD    groupXX 15460430      72 0.40923077 0.73076923
#> 2161 fam.pop    A.ASP   groupVII 15493025      73 0.53333333 0.33333333
#> 2162 fam.pop    A.CLF   groupVII 15493025      73 0.00000000 0.00000000
#> 2163 fam.pop    A.OPL   groupVII 15493025      73 0.20000000 0.10000000
#> 2164 fam.pop    A.PAL   groupVII 15493025      73 0.42857143 0.25000000
#> 2165 fam.pop    A.SMR   groupVII 15493025      73 0.25000000 0.12500000
#> 2166 fam.pop    A.UPD   groupVII 15493025      73 0.33333333 0.16666667
#> 2167 fam.pop    B.ASP   groupVII 15493025      73 0.53333333 0.33333333
#> 2168 fam.pop    B.CLF   groupVII 15493025      73 1.00000000 0.50000000
#> 2169 fam.pop    B.OPL   groupVII 15493025      73 0.42857143 0.25000000
#> 2170 fam.pop    B.PAL   groupVII 15493025      73 0.42857143 0.25000000
#> 2171 fam.pop    B.SMR   groupVII 15493025      73 0.16666667 0.08333333
#> 2172 fam.pop    B.UPD   groupVII 15493025      73 0.42857143 0.25000000
#> 2173 fam.pop    C.ASP   groupVII 15493025      73 0.53333333 0.33333333
#> 2174 fam.pop    C.CLF   groupVII 15493025      73 0.00000000 0.00000000
#> 2175 fam.pop    C.OPL   groupVII 15493025      73 0.00000000 0.00000000
#> 2176 fam.pop    C.PAL   groupVII 15493025      73 0.53333333 0.33333333
#> 2177 fam.pop    C.SMR   groupVII 15493025      73 0.35555556 0.20000000
#> 2178 fam.pop    C.UPD   groupVII 15493025      73 0.33333333 0.16666667
#> 2179 fam.pop    D.ASP   groupVII 15493025      73 0.00000000 0.00000000
#> 2180 fam.pop    D.CLF   groupVII 15493025      73 0.66666667 0.50000000
#> 2181 fam.pop    D.OPL   groupVII 15493025      73 0.36263736 0.21428571
#> 2182 fam.pop    D.PAL   groupVII 15493025      73 0.53571429 0.37500000
#> 2183 fam.pop    D.SMR   groupVII 15493025      73 0.16666667 0.08333333
#> 2184 fam.pop    D.UPD   groupVII 15493025      73 0.53571429 0.37500000
#> 2185     pop      ASP   groupVII 15493025      73 0.44210526 0.30000000
#> 2186     pop      CLF   groupVII 15493025      73 0.40909091 0.25000000
#> 2187     pop      OPL   groupVII 15493025      73 0.23188406 0.13043478
#> 2188     pop      PAL   groupVII 15493025      73 0.43448276 0.30000000
#> 2189     pop      SMR   groupVII 15493025      73 0.21486643 0.11904762
#> 2190     pop      UPD   groupVII 15493025      73 0.38888889 0.25000000
#> 2191 fam.pop    A.ASP    groupXV 15547073      74 0.33333333 0.16666667
#> 2192 fam.pop    A.CLF    groupXV 15547073      74 0.50000000 0.25000000
#> 2193 fam.pop    A.OPL    groupXV 15547073      74 0.33333333 0.16666667
#> 2194 fam.pop    A.PAL    groupXV 15547073      74 0.00000000 0.00000000
#> 2195 fam.pop    A.SMR    groupXV 15547073      74 0.00000000 0.00000000
#> 2196 fam.pop    A.UPD    groupXV 15547073      74 0.00000000 0.00000000
#> 2197 fam.pop    B.ASP    groupXV 15547073      74 0.00000000 0.00000000
#> 2198 fam.pop    B.CLF    groupXV 15547073      74 0.00000000 0.00000000
#> 2199 fam.pop    B.OPL    groupXV 15547073      74 0.49450549 0.35714286
#> 2200 fam.pop    B.PAL    groupXV 15547073      74 0.00000000 0.00000000
#> 2201 fam.pop    B.SMR    groupXV 15547073      74 0.60000000 0.50000000
#> 2202 fam.pop    B.UPD    groupXV 15547073      74 0.42857143 0.25000000
#> 2203 fam.pop    C.ASP    groupXV 15547073      74 0.66666667 0.50000000
#> 2204 fam.pop    C.CLF    groupXV 15547073      74 0.00000000 0.00000000
#> 2205 fam.pop    C.OPL    groupXV 15547073      74 0.53571429 0.37500000
#> 2206 fam.pop    C.PAL    groupXV 15547073      74 0.00000000 0.00000000
#> 2207 fam.pop    C.SMR    groupXV 15547073      74 0.46666667 0.30000000
#> 2208 fam.pop    C.UPD    groupXV 15547073      74 0.00000000 0.00000000
#> 2209 fam.pop    D.ASP    groupXV 15547073      74 0.00000000 0.00000000
#> 2210 fam.pop    D.CLF    groupXV 15547073      74 0.00000000 0.00000000
#> 2211 fam.pop    D.OPL    groupXV 15547073      74 0.46666667 0.30000000
#> 2212 fam.pop    D.PAL    groupXV 15547073      74 0.00000000 0.00000000
#> 2213 fam.pop    D.SMR    groupXV 15547073      74 0.35555556 0.20000000
#> 2214 fam.pop    D.UPD    groupXV 15547073      74 0.33333333 0.16666667
#> 2215     pop      ASP    groupXV 15547073      74 0.26842105 0.15000000
#> 2216     pop      CLF    groupXV 15547073      74 0.20000000 0.10000000
#> 2217     pop      OPL    groupXV 15547073      74 0.44381223 0.31578947
#> 2218     pop      PAL    groupXV 15547073      74 0.00000000 0.00000000
#> 2219     pop      SMR    groupXV 15547073      74 0.37076649 0.23529412
#> 2220     pop      UPD    groupXV 15547073      74 0.21230769 0.11538462
#> 2221 fam.pop    A.ASP  groupXIII 15697285      75 0.60000000 0.50000000
#> 2222 fam.pop    A.CLF  groupXIII 15697285      75 1.00000000 0.50000000
#> 2223 fam.pop    A.OPL  groupXIII 15697285      75 0.53333333 0.40000000
#> 2224 fam.pop    A.PAL  groupXIII 15697285      75 0.42857143 0.25000000
#> 2225 fam.pop    A.SMR  groupXIII 15697285      75 0.55555556 0.50000000
#> 2226 fam.pop    A.UPD  groupXIII 15697285      75 0.33333333 0.83333333
#> 2227 fam.pop    B.ASP  groupXIII 15697285      75 0.53333333 0.66666667
#> 2228 fam.pop    B.CLF  groupXIII 15697285      75 1.00000000 0.50000000
#> 2229 fam.pop    B.OPL  groupXIII 15697285      75 0.53571429 0.37500000
#> 2230 fam.pop    B.PAL  groupXIII 15697285      75 0.57142857 0.50000000
#> 2231 fam.pop    B.SMR  groupXIII 15697285      75 0.53333333 0.40000000
#> 2232 fam.pop    B.UPD  groupXIII 15697285      75 0.42857143 0.75000000
#> 2233 fam.pop    C.ASP  groupXIII 15697285      75 1.00000000 0.50000000
#> 2234 fam.pop    C.CLF  groupXIII 15697285      75 0.00000000 1.00000000
#> 2235 fam.pop    C.OPL  groupXIII 15697285      75 0.52747253 0.42857143
#> 2236 fam.pop    C.PAL  groupXIII 15697285      75 0.42857143 0.25000000
#> 2237 fam.pop    C.SMR  groupXIII 15697285      75 0.53333333 0.60000000
#> 2238 fam.pop    C.UPD  groupXIII 15697285      75 0.57142857 0.50000000
#> 2239 fam.pop    D.ASP  groupXIII 15697285      75 0.66666667 0.50000000
#> 2240 fam.pop    D.CLF  groupXIII 15697285      75 0.50000000 0.25000000
#> 2241 fam.pop    D.OPL  groupXIII 15697285      75 0.53846154 0.50000000
#> 2242 fam.pop    D.PAL  groupXIII 15697285      75 0.42857143 0.25000000
#> 2243 fam.pop    D.SMR  groupXIII 15697285      75 0.53333333 0.40000000
#> 2244 fam.pop    D.UPD  groupXIII 15697285      75 0.42857143 0.75000000
#> 2245     pop      ASP  groupXIII 15697285      75 0.52287582 0.55555556
#> 2246     pop      CLF  groupXIII 15697285      75 0.55555556 0.50000000
#> 2247     pop      OPL  groupXIII 15697285      75 0.50241546 0.43478261
#> 2248     pop      PAL  groupXIII 15697285      75 0.44354839 0.31250000
#> 2249     pop      SMR  groupXIII 15697285      75 0.51153846 0.47500000
#> 2250     pop      UPD  groupXIII 15697285      75 0.43448276 0.70000000
#> 2251 fam.pop    A.ASP   groupIII 15983994      76 0.00000000 0.00000000
#> 2252 fam.pop    A.CLF   groupIII 15983994      76 0.00000000 0.00000000
#> 2253 fam.pop    A.OPL   groupIII 15983994      76 0.53333333 0.40000000
#> 2254 fam.pop    A.PAL   groupIII 15983994      76 0.00000000 1.00000000
#> 2255 fam.pop    A.SMR   groupIII 15983994      76 0.25000000 0.12500000
#> 2256 fam.pop    A.UPD   groupIII 15983994      76 0.60000000 0.50000000
#> 2257 fam.pop    B.ASP   groupIII 15983994      76 0.50000000 0.25000000
#> 2258 fam.pop    B.CLF   groupIII 15983994      76        NaN         NA
#> 2259 fam.pop    B.OPL   groupIII 15983994      76 0.55555556 0.50000000
#> 2260 fam.pop    B.PAL   groupIII 15983994      76 0.00000000 1.00000000
#> 2261 fam.pop    B.SMR   groupIII 15983994      76 0.50000000 0.25000000
#> 2262 fam.pop    B.UPD   groupIII 15983994      76 0.42857143 0.25000000
#> 2263 fam.pop    C.ASP   groupIII 15983994      76 0.50000000 0.75000000
#> 2264 fam.pop    C.CLF   groupIII 15983994      76        NaN         NA
#> 2265 fam.pop    C.OPL   groupIII 15983994      76 0.50000000 0.25000000
#> 2266 fam.pop    C.PAL   groupIII 15983994      76 0.00000000 1.00000000
#> 2267 fam.pop    C.SMR   groupIII 15983994      76 0.00000000 0.00000000
#> 2268 fam.pop    C.UPD   groupIII 15983994      76 0.60000000 0.50000000
#> 2269 fam.pop    D.ASP   groupIII 15983994      76 1.00000000 0.50000000
#> 2270 fam.pop    D.CLF   groupIII 15983994      76 0.50000000 0.25000000
#> 2271 fam.pop    D.OPL   groupIII 15983994      76 0.25000000 0.12500000
#> 2272 fam.pop    D.PAL   groupIII 15983994      76 0.00000000 1.00000000
#> 2273 fam.pop    D.SMR   groupIII 15983994      76 0.25000000 0.12500000
#> 2274 fam.pop    D.UPD   groupIII 15983994      76 0.57142857 0.50000000
#> 2275     pop      ASP   groupIII 15983994      76 0.53030303 0.41666667
#> 2276     pop      CLF   groupIII 15983994      76 0.33333333 0.16666667
#> 2277     pop      OPL   groupIII 15983994      76 0.46572581 0.34375000
#> 2278     pop      PAL   groupIII 15983994      76 0.00000000 1.00000000
#> 2279     pop      SMR   groupIII 15983994      76 0.22826087 0.12500000
#> 2280     pop      UPD   groupIII 15983994      76 0.50793651 0.42857143
#> 2281 fam.pop    A.ASP   groupXIX 16233699      77 0.53333333 0.33333333
#> 2282 fam.pop    A.CLF   groupXIX 16233699      77 0.00000000 0.00000000
#> 2283 fam.pop    A.OPL   groupXIX 16233699      77 0.25000000 0.12500000
#> 2284 fam.pop    A.PAL   groupXIX 16233699      77 0.46666667 0.30000000
#> 2285 fam.pop    A.SMR   groupXIX 16233699      77 0.30303030 0.16666667
#> 2286 fam.pop    A.UPD   groupXIX 16233699      77 0.60000000 0.50000000
#> 2287 fam.pop    B.ASP   groupXIX 16233699      77 0.50000000 0.25000000
#> 2288 fam.pop    B.CLF   groupXIX 16233699      77 1.00000000 0.50000000
#> 2289 fam.pop    B.OPL   groupXIX 16233699      77 0.30303030 0.16666667
#> 2290 fam.pop    B.PAL   groupXIX 16233699      77 0.60000000 0.50000000
#> 2291 fam.pop    B.SMR   groupXIX 16233699      77 0.35555556 0.20000000
#> 2292 fam.pop    B.UPD   groupXIX 16233699      77 0.25000000 0.12500000
#> 2293 fam.pop    C.ASP   groupXIX 16233699      77 0.53333333 0.33333333
#> 2294 fam.pop    C.CLF   groupXIX 16233699      77 0.00000000 0.00000000
#> 2295 fam.pop    C.OPL   groupXIX 16233699      77 0.36263736 0.21428571
#> 2296 fam.pop    C.PAL   groupXIX 16233699      77 0.57142857 0.50000000
#> 2297 fam.pop    C.SMR   groupXIX 16233699      77 0.53333333 0.33333333
#> 2298 fam.pop    C.UPD   groupXIX 16233699      77 0.33333333 0.16666667
#> 2299 fam.pop    D.ASP   groupXIX 16233699      77 1.00000000 0.50000000
#> 2300 fam.pop    D.CLF   groupXIX 16233699      77 0.00000000 0.00000000
#> 2301 fam.pop    D.OPL   groupXIX 16233699      77 0.30303030 0.16666667
#> 2302 fam.pop    D.PAL   groupXIX 16233699      77 0.33333333 0.16666667
#> 2303 fam.pop    D.SMR   groupXIX 16233699      77 0.20000000 0.10000000
#> 2304 fam.pop    D.UPD   groupXIX 16233699      77 0.25000000 0.12500000
#> 2305     pop      ASP   groupXIX 16233699      77 0.47058824 0.33333333
#> 2306     pop      CLF   groupXIX 16233699      77 0.16666667 0.08333333
#> 2307     pop      OPL   groupXIX 16233699      77 0.29371981 0.17391304
#> 2308     pop      PAL   groupXIX 16233699      77 0.48045977 0.36666667
#> 2309     pop      SMR   groupXIX 16233699      77 0.30867710 0.18421053
#> 2310     pop      UPD   groupXIX 16233699      77 0.34920635 0.21428571
#> 2311 fam.pop    A.ASP    groupIX 16306966      78 0.00000000 0.00000000
#> 2312 fam.pop    A.CLF    groupIX 16306966      78 0.00000000 0.00000000
#> 2313 fam.pop    A.OPL    groupIX 16306966      78 0.25000000 0.12500000
#> 2314 fam.pop    A.PAL    groupIX 16306966      78 0.00000000 0.00000000
#> 2315 fam.pop    A.SMR    groupIX 16306966      78 0.40909091 0.25000000
#> 2316 fam.pop    A.UPD    groupIX 16306966      78 0.00000000 0.00000000
#> 2317 fam.pop    B.ASP    groupIX 16306966      78 0.00000000 0.00000000
#> 2318 fam.pop    B.CLF    groupIX 16306966      78 0.00000000 0.00000000
#> 2319 fam.pop    B.OPL    groupIX 16306966      78 0.48484848 0.33333333
#> 2320 fam.pop    B.PAL    groupIX 16306966      78 0.00000000 0.00000000
#> 2321 fam.pop    B.SMR    groupIX 16306966      78 0.00000000 0.00000000
#> 2322 fam.pop    B.UPD    groupIX 16306966      78 0.00000000 0.00000000
#> 2323 fam.pop    C.ASP    groupIX 16306966      78 0.00000000 0.00000000
#> 2324 fam.pop    C.CLF    groupIX 16306966      78 0.00000000 0.00000000
#> 2325 fam.pop    C.OPL    groupIX 16306966      78 0.52747253 0.42857143
#> 2326 fam.pop    C.PAL    groupIX 16306966      78 0.00000000 0.00000000
#> 2327 fam.pop    C.SMR    groupIX 16306966      78 0.16666667 0.08333333
#> 2328 fam.pop    C.UPD    groupIX 16306966      78 0.00000000 0.00000000
#> 2329 fam.pop    D.ASP    groupIX 16306966      78 0.00000000 0.00000000
#> 2330 fam.pop    D.CLF    groupIX 16306966      78 0.50000000 0.75000000
#> 2331 fam.pop    D.OPL    groupIX 16306966      78 0.35555556 0.20000000
#> 2332 fam.pop    D.PAL    groupIX 16306966      78 0.00000000 0.00000000
#> 2333 fam.pop    D.SMR    groupIX 16306966      78 0.20000000 0.10000000
#> 2334 fam.pop    D.UPD    groupIX 16306966      78 0.00000000 0.00000000
#> 2335     pop      ASP    groupIX 16306966      78 0.00000000 0.00000000
#> 2336     pop      CLF    groupIX 16306966      78 0.40909091 0.25000000
#> 2337     pop      OPL    groupIX 16306966      78 0.42600423 0.29545455
#> 2338     pop      PAL    groupIX 16306966      78 0.00000000 0.00000000
#> 2339     pop      SMR    groupIX 16306966      78 0.20613108 0.11363636
#> 2340     pop      UPD    groupIX 16306966      78 0.00000000 0.00000000
#> 2341 fam.pop    A.ASP     groupI 16498302      79 0.33333333 0.16666667
#> 2342 fam.pop    A.CLF     groupI 16498302      79 0.00000000 0.00000000
#> 2343 fam.pop    A.OPL     groupI 16498302      79 0.20000000 0.10000000
#> 2344 fam.pop    A.PAL     groupI 16498302      79 0.53333333 0.40000000
#> 2345 fam.pop    A.SMR     groupI 16498302      79 0.00000000 0.00000000
#> 2346 fam.pop    A.UPD     groupI 16498302      79 0.33333333 0.16666667
#> 2347 fam.pop    B.ASP     groupI 16498302      79 0.33333333 0.16666667
#> 2348 fam.pop    B.CLF     groupI 16498302      79 0.00000000 0.00000000
#> 2349 fam.pop    B.OPL     groupI 16498302      79 0.36263736 0.21428571
#> 2350 fam.pop    B.PAL     groupI 16498302      79 0.25000000 0.12500000
#> 2351 fam.pop    B.SMR     groupI 16498302      79 0.42857143 0.25000000
#> 2352 fam.pop    B.UPD     groupI 16498302      79 0.66666667 0.50000000
#> 2353 fam.pop    C.ASP     groupI 16498302      79 0.53333333 0.33333333
#> 2354 fam.pop    C.CLF     groupI 16498302      79 0.00000000 0.00000000
#> 2355 fam.pop    C.OPL     groupI 16498302      79 0.00000000 0.00000000
#> 2356 fam.pop    C.PAL     groupI 16498302      79 0.25000000 0.12500000
#> 2357 fam.pop    C.SMR     groupI 16498302      79 0.16666667 0.08333333
#> 2358 fam.pop    C.UPD     groupI 16498302      79 0.00000000 0.00000000
#> 2359 fam.pop    D.ASP     groupI 16498302      79 0.00000000 0.00000000
#> 2360 fam.pop    D.CLF     groupI 16498302      79 0.00000000 0.00000000
#> 2361 fam.pop    D.OPL     groupI 16498302      79 0.26373626 0.14285714
#> 2362 fam.pop    D.PAL     groupI 16498302      79 0.42857143 0.25000000
#> 2363 fam.pop    D.SMR     groupI 16498302      79 0.16666667 0.08333333
#> 2364 fam.pop    D.UPD     groupI 16498302      79 0.00000000 0.00000000
#> 2365     pop      ASP     groupI 16498302      79 0.31168831 0.18181818
#> 2366     pop      CLF     groupI 16498302      79 0.00000000 0.00000000
#> 2367     pop      OPL     groupI 16498302      79 0.20814480 0.11538462
#> 2368     pop      PAL     groupI 16498302      79 0.37076649 0.23529412
#> 2369     pop      SMR     groupI 16498302      79 0.16913319 0.09090909
#> 2370     pop      UPD     groupI 16498302      79 0.22826087 0.12500000
#> 2371 fam.pop    A.ASP  groupXIII 16520929      80 0.60000000 0.50000000
#> 2372 fam.pop    A.CLF  groupXIII 16520929      80 0.50000000 0.25000000
#> 2373 fam.pop    A.OPL  groupXIII 16520929      80 0.35555556 0.20000000
#> 2374 fam.pop    A.PAL  groupXIII 16520929      80 0.55555556 0.50000000
#> 2375 fam.pop    A.SMR  groupXIII 16520929      80 0.48484848 0.33333333
#> 2376 fam.pop    A.UPD  groupXIII 16520929      80 0.53333333 0.33333333
#> 2377 fam.pop    B.ASP  groupXIII 16520929      80 0.53333333 0.33333333
#> 2378 fam.pop    B.CLF  groupXIII 16520929      80 0.00000000 0.00000000
#> 2379 fam.pop    B.OPL  groupXIII 16520929      80 0.55555556 0.50000000
#> 2380 fam.pop    B.PAL  groupXIII 16520929      80 0.53571429 0.62500000
#> 2381 fam.pop    B.SMR  groupXIII 16520929      80 0.30303030 0.16666667
#> 2382 fam.pop    B.UPD  groupXIII 16520929      80 0.42857143 0.25000000
#> 2383 fam.pop    C.ASP  groupXIII 16520929      80 0.00000000 0.00000000
#> 2384 fam.pop    C.CLF  groupXIII 16520929      80 0.00000000 0.00000000
#> 2385 fam.pop    C.OPL  groupXIII 16520929      80 0.40909091 0.25000000
#> 2386 fam.pop    C.PAL  groupXIII 16520929      80 0.57142857 0.50000000
#> 2387 fam.pop    C.SMR  groupXIII 16520929      80 0.46666667 0.30000000
#> 2388 fam.pop    C.UPD  groupXIII 16520929      80 0.42857143 0.25000000
#> 2389 fam.pop    D.ASP  groupXIII 16520929      80 0.00000000 0.00000000
#> 2390 fam.pop    D.CLF  groupXIII 16520929      80 0.50000000 0.25000000
#> 2391 fam.pop    D.OPL  groupXIII 16520929      80 0.48484848 0.33333333
#> 2392 fam.pop    D.PAL  groupXIII 16520929      80 0.50000000 0.25000000
#> 2393 fam.pop    D.SMR  groupXIII 16520929      80 0.40909091 0.25000000
#> 2394 fam.pop    D.UPD  groupXIII 16520929      80 0.50000000 0.25000000
#> 2395     pop      ASP  groupXIII 16520929      80 0.36796537 0.22727273
#> 2396     pop      CLF  groupXIII 16520929      80 0.30303030 0.16666667
#> 2397     pop      OPL  groupXIII 16520929      80 0.44397463 0.31818182
#> 2398     pop      PAL  groupXIII 16520929      80 0.51724138 0.50000000
#> 2399     pop      SMR  groupXIII 16520929      80 0.39420290 0.26086957
#> 2400     pop      UPD  groupXIII 16520929      80 0.40923077 0.26923077
#> 2401 fam.pop    A.ASP   groupXVI 16561891      81 0.53333333 0.33333333
#> 2402 fam.pop    A.CLF   groupXVI 16561891      81 1.00000000 0.50000000
#> 2403 fam.pop    A.OPL   groupXVI 16561891      81 0.53333333 0.40000000
#> 2404 fam.pop    A.PAL   groupXVI 16561891      81 0.53333333 0.60000000
#> 2405 fam.pop    A.SMR   groupXVI 16561891      81 0.48484848 0.33333333
#> 2406 fam.pop    A.UPD   groupXVI 16561891      81 0.66666667 0.50000000
#> 2407 fam.pop    B.ASP   groupXVI 16561891      81 1.00000000 0.50000000
#> 2408 fam.pop    B.CLF   groupXVI 16561891      81 1.00000000 0.50000000
#> 2409 fam.pop    B.OPL   groupXVI 16561891      81 0.46666667 0.30000000
#> 2410 fam.pop    B.PAL   groupXVI 16561891      81 0.53333333 0.33333333
#> 2411 fam.pop    B.SMR   groupXVI 16561891      81 0.35555556 0.20000000
#> 2412 fam.pop    B.UPD   groupXVI 16561891      81 0.42857143 0.25000000
#> 2413 fam.pop    C.ASP   groupXVI 16561891      81 0.53333333 0.33333333
#> 2414 fam.pop    C.CLF   groupXVI 16561891      81 1.00000000 0.50000000
#> 2415 fam.pop    C.OPL   groupXVI 16561891      81 0.53333333 0.60000000
#> 2416 fam.pop    C.PAL   groupXVI 16561891      81 0.42857143 0.25000000
#> 2417 fam.pop    C.SMR   groupXVI 16561891      81 0.46666667 0.30000000
#> 2418 fam.pop    C.UPD   groupXVI 16561891      81 0.53571429 0.37500000
#> 2419 fam.pop    D.ASP   groupXVI 16561891      81 0.50000000 0.25000000
#> 2420 fam.pop    D.CLF   groupXVI 16561891      81 0.50000000 0.75000000
#> 2421 fam.pop    D.OPL   groupXVI 16561891      81 0.48484848 0.33333333
#> 2422 fam.pop    D.PAL   groupXVI 16561891      81 0.57142857 0.50000000
#> 2423 fam.pop    D.SMR   groupXVI 16561891      81 0.53030303 0.41666667
#> 2424 fam.pop    D.UPD   groupXVI 16561891      81 0.57142857 0.50000000
#> 2425     pop      ASP   groupXVI 16561891      81 0.47058824 0.33333333
#> 2426     pop      CLF   groupXVI 16561891      81 0.53333333 0.60000000
#> 2427     pop      OPL   groupXVI 16561891      81 0.49361208 0.40476190
#> 2428     pop      PAL   groupXVI 16561891      81 0.50806452 0.43750000
#> 2429     pop      SMR   groupXVI 16561891      81 0.44397463 0.31818182
#> 2430     pop      UPD   groupXVI 16561891      81 0.49470899 0.39285714
#> 2431 fam.pop    A.ASP    groupII 16598364      82 0.33333333 0.16666667
#> 2432 fam.pop    A.CLF    groupII 16598364      82 0.00000000 0.00000000
#> 2433 fam.pop    A.OPL    groupII 16598364      82 0.16666667 0.08333333
#> 2434 fam.pop    A.PAL    groupII 16598364      82 0.00000000 0.00000000
#> 2435 fam.pop    A.SMR    groupII 16598364      82 0.00000000 0.00000000
#> 2436 fam.pop    A.UPD    groupII 16598364      82 0.53333333 0.33333333
#> 2437 fam.pop    B.ASP    groupII 16598364      82 0.00000000 0.00000000
#> 2438 fam.pop    B.CLF    groupII 16598364      82        NaN         NA
#> 2439 fam.pop    B.OPL    groupII 16598364      82 0.36263736 0.21428571
#> 2440 fam.pop    B.PAL    groupII 16598364      82 0.00000000 0.00000000
#> 2441 fam.pop    B.SMR    groupII 16598364      82 0.00000000 0.00000000
#> 2442 fam.pop    B.UPD    groupII 16598364      82 0.00000000 0.00000000
#> 2443 fam.pop    C.ASP    groupII 16598364      82 0.50000000 0.25000000
#> 2444 fam.pop    C.CLF    groupII 16598364      82 0.00000000 0.00000000
#> 2445 fam.pop    C.OPL    groupII 16598364      82 0.16666667 0.08333333
#> 2446 fam.pop    C.PAL    groupII 16598364      82 0.00000000 0.00000000
#> 2447 fam.pop    C.SMR    groupII 16598364      82 0.00000000 0.00000000
#> 2448 fam.pop    C.UPD    groupII 16598364      82 0.33333333 0.16666667
#> 2449 fam.pop    D.ASP    groupII 16598364      82 0.00000000 0.00000000
#> 2450 fam.pop    D.CLF    groupII 16598364      82 0.00000000 0.00000000
#> 2451 fam.pop    D.OPL    groupII 16598364      82 0.00000000 0.00000000
#> 2452 fam.pop    D.PAL    groupII 16598364      82 0.00000000 0.00000000
#> 2453 fam.pop    D.SMR    groupII 16598364      82 0.30303030 0.16666667
#> 2454 fam.pop    D.UPD    groupII 16598364      82 0.53333333 0.33333333
#> 2455     pop      ASP    groupII 16598364      82 0.18947368 0.10000000
#> 2456     pop      CLF    groupII 16598364      82 0.00000000 0.00000000
#> 2457     pop      OPL    groupII 16598364      82 0.17722474 0.09615385
#> 2458     pop      PAL    groupII 16598364      82 0.00000000 0.00000000
#> 2459     pop      SMR    groupII 16598364      82 0.10241821 0.05263158
#> 2460     pop      UPD    groupII 16598364      82 0.32307692 0.19230769
#> 2461 fam.pop    A.ASP    groupIV 17248473      83 0.33333333 0.16666667
#> 2462 fam.pop    A.CLF    groupIV 17248473      83 0.66666667 0.50000000
#> 2463 fam.pop    A.OPL    groupIV 17248473      83 0.00000000 0.00000000
#> 2464 fam.pop    A.PAL    groupIV 17248473      83 0.00000000 0.00000000
#> 2465 fam.pop    A.SMR    groupIV 17248473      83 0.16666667 0.08333333
#> 2466 fam.pop    A.UPD    groupIV 17248473      83 0.33333333 0.16666667
#> 2467 fam.pop    B.ASP    groupIV 17248473      83 0.53333333 0.33333333
#> 2468 fam.pop    B.CLF    groupIV 17248473      83 0.00000000 0.00000000
#> 2469 fam.pop    B.OPL    groupIV 17248473      83 0.00000000 0.00000000
#> 2470 fam.pop    B.PAL    groupIV 17248473      83 0.00000000 0.00000000
#> 2471 fam.pop    B.SMR    groupIV 17248473      83 0.00000000 0.00000000
#> 2472 fam.pop    B.UPD    groupIV 17248473      83 0.00000000 0.00000000
#> 2473 fam.pop    C.ASP    groupIV 17248473      83 0.33333333 0.16666667
#> 2474 fam.pop    C.CLF    groupIV 17248473      83 1.00000000 0.50000000
#> 2475 fam.pop    C.OPL    groupIV 17248473      83 0.00000000 0.00000000
#> 2476 fam.pop    C.PAL    groupIV 17248473      83 0.00000000 0.00000000
#> 2477 fam.pop    C.SMR    groupIV 17248473      83 0.00000000 0.00000000
#> 2478 fam.pop    C.UPD    groupIV 17248473      83 0.00000000 0.00000000
#> 2479 fam.pop    D.ASP    groupIV 17248473      83 0.00000000 0.00000000
#> 2480 fam.pop    D.CLF    groupIV 17248473      83 0.66666667 0.50000000
#> 2481 fam.pop    D.OPL    groupIV 17248473      83 0.00000000 0.00000000
#> 2482 fam.pop    D.PAL    groupIV 17248473      83 0.00000000 0.00000000
#> 2483 fam.pop    D.SMR    groupIV 17248473      83 0.00000000 0.00000000
#> 2484 fam.pop    D.UPD    groupIV 17248473      83 0.00000000 0.00000000
#> 2485     pop      ASP    groupIV 17248473      83 0.31168831 0.18181818
#> 2486     pop      CLF    groupIV 17248473      83 0.53030303 0.41666667
#> 2487     pop      OPL    groupIV 17248473      83 0.00000000 0.00000000
#> 2488     pop      PAL    groupIV 17248473      83 0.00000000 0.00000000
#> 2489     pop      SMR    groupIV 17248473      83 0.04166667 0.02083333
#> 2490     pop      UPD    groupIV 17248473      83 0.07142857 0.03571429
#> 2491 fam.pop    A.ASP    groupXX 17748944      84 0.33333333 0.16666667
#> 2492 fam.pop    A.CLF    groupXX 17748944      84 0.50000000 0.25000000
#> 2493 fam.pop    A.OPL    groupXX 17748944      84 0.53571429 0.62500000
#> 2494 fam.pop    A.PAL    groupXX 17748944      84 0.00000000 1.00000000
#> 2495 fam.pop    A.SMR    groupXX 17748944      84 0.48484848 0.33333333
#> 2496 fam.pop    A.UPD    groupXX 17748944      84 0.53333333 0.33333333
#> 2497 fam.pop    B.ASP    groupXX 17748944      84 0.00000000 0.00000000
#> 2498 fam.pop    B.CLF    groupXX 17748944      84 0.00000000 0.00000000
#> 2499 fam.pop    B.OPL    groupXX 17748944      84 0.48484848 0.33333333
#> 2500 fam.pop    B.PAL    groupXX 17748944      84 0.00000000 1.00000000
#> 2501 fam.pop    B.SMR    groupXX 17748944      84 0.53571429 0.37500000
#> 2502 fam.pop    B.UPD    groupXX 17748944      84 0.53571429 0.62500000
#> 2503 fam.pop    C.ASP    groupXX 17748944      84 0.53333333 0.66666667
#> 2504 fam.pop    C.CLF    groupXX 17748944      84 0.00000000 0.00000000
#> 2505 fam.pop    C.OPL    groupXX 17748944      84 0.36263736 0.21428571
#> 2506 fam.pop    C.PAL    groupXX 17748944      84 0.00000000 1.00000000
#> 2507 fam.pop    C.SMR    groupXX 17748944      84 0.42857143 0.25000000
#> 2508 fam.pop    C.UPD    groupXX 17748944      84 0.53333333 0.66666667
#> 2509 fam.pop    D.ASP    groupXX 17748944      84 0.50000000 0.25000000
#> 2510 fam.pop    D.CLF    groupXX 17748944      84 0.50000000 0.75000000
#> 2511 fam.pop    D.OPL    groupXX 17748944      84 0.53333333 0.33333333
#> 2512 fam.pop    D.PAL    groupXX 17748944      84 0.00000000 1.00000000
#> 2513 fam.pop    D.SMR    groupXX 17748944      84 0.20000000 0.10000000
#> 2514 fam.pop    D.UPD    groupXX 17748944      84 0.33333333 0.16666667
#> 2515     pop      ASP    groupXX 17748944      84 0.44210526 0.30000000
#> 2516     pop      CLF    groupXX 17748944      84 0.48484848 0.33333333
#> 2517     pop      OPL    groupXX 17748944      84 0.46666667 0.35000000
#> 2518     pop      PAL    groupXX 17748944      84 0.00000000 1.00000000
#> 2519     pop      SMR    groupXX 17748944      84 0.39829303 0.26315789
#> 2520     pop      UPD    groupXX 17748944      84 0.51692308 0.46153846
#> 2521 fam.pop    A.ASP  groupXIII 17842228      85 0.60000000 0.50000000
#> 2522 fam.pop    A.CLF  groupXIII 17842228      85 0.50000000 0.75000000
#> 2523 fam.pop    A.OPL  groupXIII 17842228      85 0.53030303 0.41666667
#> 2524 fam.pop    A.PAL  groupXIII 17842228      85 0.00000000 1.00000000
#> 2525 fam.pop    A.SMR  groupXIII 17842228      85 0.53571429 0.37500000
#> 2526 fam.pop    A.UPD  groupXIII 17842228      85 0.53333333 0.33333333
#> 2527 fam.pop    B.ASP  groupXIII 17842228      85 0.50000000 0.75000000
#> 2528 fam.pop    B.CLF  groupXIII 17842228      85        NaN         NA
#> 2529 fam.pop    B.OPL  groupXIII 17842228      85 0.33333333 0.16666667
#> 2530 fam.pop    B.PAL  groupXIII 17842228      85 0.00000000 1.00000000
#> 2531 fam.pop    B.SMR  groupXIII 17842228      85 0.30303030 0.16666667
#> 2532 fam.pop    B.UPD  groupXIII 17842228      85 0.53571429 0.37500000
#> 2533 fam.pop    C.ASP  groupXIII 17842228      85 0.60000000 0.50000000
#> 2534 fam.pop    C.CLF  groupXIII 17842228      85 1.00000000 0.50000000
#> 2535 fam.pop    C.OPL  groupXIII 17842228      85 0.49450549 0.35714286
#> 2536 fam.pop    C.PAL  groupXIII 17842228      85 0.33333333 0.83333333
#> 2537 fam.pop    C.SMR  groupXIII 17842228      85 0.35555556 0.20000000
#> 2538 fam.pop    C.UPD  groupXIII 17842228      85 0.57142857 0.50000000
#> 2539 fam.pop    D.ASP  groupXIII 17842228      85 0.00000000 1.00000000
#> 2540 fam.pop    D.CLF  groupXIII 17842228      85 0.66666667 0.50000000
#> 2541 fam.pop    D.OPL  groupXIII 17842228      85 0.35555556 0.20000000
#> 2542 fam.pop    D.PAL  groupXIII 17842228      85 0.53333333 0.66666667
#> 2543 fam.pop    D.SMR  groupXIII 17842228      85 0.48484848 0.33333333
#> 2544 fam.pop    D.UPD  groupXIII 17842228      85 0.57142857 0.50000000
#> 2545     pop      ASP  groupXIII 17842228      85 0.47894737 0.65000000
#> 2546     pop      CLF  groupXIII 17842228      85 0.53333333 0.60000000
#> 2547     pop      OPL  groupXIII 17842228      85 0.43786295 0.30952381
#> 2548     pop      PAL  groupXIII 17842228      85 0.19841270 0.89285714
#> 2549     pop      SMR  groupXIII 17842228      85 0.39605110 0.26190476
#> 2550     pop      UPD  groupXIII 17842228      85 0.50804598 0.43333333
#> 2551 fam.pop    A.ASP   groupXIX 17858769      86 0.53333333 0.33333333
#> 2552 fam.pop    A.CLF   groupXIX 17858769      86 0.00000000 0.00000000
#> 2553 fam.pop    A.OPL   groupXIX 17858769      86 0.42857143 0.25000000
#> 2554 fam.pop    A.PAL   groupXIX 17858769      86 0.53333333 0.33333333
#> 2555 fam.pop    A.SMR   groupXIX 17858769      86 0.42857143 0.25000000
#> 2556 fam.pop    A.UPD   groupXIX 17858769      86 0.60000000 0.50000000
#> 2557 fam.pop    B.ASP   groupXIX 17858769      86 0.50000000 0.25000000
#> 2558 fam.pop    B.CLF   groupXIX 17858769      86 1.00000000 0.50000000
#> 2559 fam.pop    B.OPL   groupXIX 17858769      86 0.42857143 0.25000000
#> 2560 fam.pop    B.PAL   groupXIX 17858769      86 0.66666667 0.50000000
#> 2561 fam.pop    B.SMR   groupXIX 17858769      86 0.30303030 0.16666667
#> 2562 fam.pop    B.UPD   groupXIX 17858769      86 0.33333333 0.16666667
#> 2563 fam.pop    C.ASP   groupXIX 17858769      86 0.66666667 0.50000000
#> 2564 fam.pop    C.CLF   groupXIX 17858769      86 0.00000000 0.00000000
#> 2565 fam.pop    C.OPL   groupXIX 17858769      86 0.36263736 0.21428571
#> 2566 fam.pop    C.PAL   groupXIX 17858769      86 0.57142857 0.50000000
#> 2567 fam.pop    C.SMR   groupXIX 17858769      86 0.42857143 0.25000000
#> 2568 fam.pop    C.UPD   groupXIX 17858769      86 0.53333333 0.33333333
#> 2569 fam.pop    D.ASP   groupXIX 17858769      86 1.00000000 0.50000000
#> 2570 fam.pop    D.CLF   groupXIX 17858769      86 0.00000000 0.00000000
#> 2571 fam.pop    D.OPL   groupXIX 17858769      86 0.35555556 0.20000000
#> 2572 fam.pop    D.PAL   groupXIX 17858769      86 0.50000000 0.25000000
#> 2573 fam.pop    D.SMR   groupXIX 17858769      86 0.20000000 0.10000000
#> 2574 fam.pop    D.UPD   groupXIX 17858769      86 0.25000000 0.12500000
#> 2575     pop      ASP   groupXIX 17858769      86 0.50000000 0.37500000
#> 2576     pop      CLF   groupXIX 17858769      86 0.16666667 0.08333333
#> 2577     pop      OPL   groupXIX 17858769      86 0.35769231 0.22500000
#> 2578     pop      PAL   groupXIX 17858769      86 0.50649351 0.40909091
#> 2579     pop      SMR   groupXIX 17858769      86 0.30867710 0.18421053
#> 2580     pop      UPD   groupXIX 17858769      86 0.40923077 0.26923077
#> 2581 fam.pop    A.ASP    groupIX 17929225      87 0.33333333 0.16666667
#> 2582 fam.pop    A.CLF    groupIX 17929225      87 0.50000000 0.25000000
#> 2583 fam.pop    A.OPL    groupIX 17929225      87 0.48484848 0.33333333
#> 2584 fam.pop    A.PAL    groupIX 17929225      87 0.00000000 0.00000000
#> 2585 fam.pop    A.SMR    groupIX 17929225      87 0.53333333 0.40000000
#> 2586 fam.pop    A.UPD    groupIX 17929225      87 0.60000000 0.50000000
#> 2587 fam.pop    B.ASP    groupIX 17929225      87 0.60000000 0.50000000
#> 2588 fam.pop    B.CLF    groupIX 17929225      87 0.00000000 0.00000000
#> 2589 fam.pop    B.OPL    groupIX 17929225      87 0.48484848 0.33333333
#> 2590 fam.pop    B.PAL    groupIX 17929225      87 0.00000000 0.00000000
#> 2591 fam.pop    B.SMR    groupIX 17929225      87 0.48484848 0.33333333
#> 2592 fam.pop    B.UPD    groupIX 17929225      87 0.42857143 0.25000000
#> 2593 fam.pop    C.ASP    groupIX 17929225      87 0.33333333 0.83333333
#> 2594 fam.pop    C.CLF    groupIX 17929225      87 0.00000000 0.00000000
#> 2595 fam.pop    C.OPL    groupIX 17929225      87 0.30303030 0.16666667
#> 2596 fam.pop    C.PAL    groupIX 17929225      87 0.00000000 0.00000000
#> 2597 fam.pop    C.SMR    groupIX 17929225      87 0.35555556 0.20000000
#> 2598 fam.pop    C.UPD    groupIX 17929225      87 0.57142857 0.50000000
#> 2599 fam.pop    D.ASP    groupIX 17929225      87 0.66666667 0.50000000
#> 2600 fam.pop    D.CLF    groupIX 17929225      87 0.00000000 0.00000000
#> 2601 fam.pop    D.OPL    groupIX 17929225      87 0.26373626 0.14285714
#> 2602 fam.pop    D.PAL    groupIX 17929225      87 0.00000000 0.00000000
#> 2603 fam.pop    D.SMR    groupIX 17929225      87 0.30303030 0.16666667
#> 2604 fam.pop    D.UPD    groupIX 17929225      87 0.60000000 0.50000000
#> 2605     pop      ASP    groupIX 17929225      87 0.52380952 0.50000000
#> 2606     pop      CLF    groupIX 17929225      87 0.16666667 0.08333333
#> 2607     pop      OPL    groupIX 17929225      87 0.37224490 0.24000000
#> 2608     pop      PAL    groupIX 17929225      87 0.00000000 0.00000000
#> 2609     pop      SMR    groupIX 17929225      87 0.40591966 0.27272727
#> 2610     pop      UPD    groupIX 17929225      87 0.50793651 0.42857143
#> 2611 fam.pop    A.ASP    groupIX 18339767      88 0.53333333 0.66666667
#> 2612 fam.pop    A.CLF    groupIX 18339767      88 0.50000000 0.25000000
#> 2613 fam.pop    A.OPL    groupIX 18339767      88 0.42857143 0.25000000
#> 2614 fam.pop    A.PAL    groupIX 18339767      88 0.00000000 1.00000000
#> 2615 fam.pop    A.SMR    groupIX 18339767      88 0.53571429 0.62500000
#> 2616 fam.pop    A.UPD    groupIX 18339767      88 0.66666667 0.50000000
#> 2617 fam.pop    B.ASP    groupIX 18339767      88 0.66666667 0.50000000
#> 2618 fam.pop    B.CLF    groupIX 18339767      88 1.00000000 0.50000000
#> 2619 fam.pop    B.OPL    groupIX 18339767      88 0.53333333 0.60000000
#> 2620 fam.pop    B.PAL    groupIX 18339767      88 0.00000000 1.00000000
#> 2621 fam.pop    B.SMR    groupIX 18339767      88 0.42857143 0.25000000
#> 2622 fam.pop    B.UPD    groupIX 18339767      88 0.53333333 0.33333333
#> 2623 fam.pop    C.ASP    groupIX 18339767      88 0.50000000 0.25000000
#> 2624 fam.pop    C.CLF    groupIX 18339767      88 0.00000000 0.00000000
#> 2625 fam.pop    C.OPL    groupIX 18339767      88 0.48484848 0.33333333
#> 2626 fam.pop    C.PAL    groupIX 18339767      88 0.00000000 1.00000000
#> 2627 fam.pop    C.SMR    groupIX 18339767      88 0.46666667 0.30000000
#> 2628 fam.pop    C.UPD    groupIX 18339767      88 0.33333333 0.16666667
#> 2629 fam.pop    D.ASP    groupIX 18339767      88 1.00000000 0.50000000
#> 2630 fam.pop    D.CLF    groupIX 18339767      88 0.50000000 0.25000000
#> 2631 fam.pop    D.OPL    groupIX 18339767      88 0.20000000 0.10000000
#> 2632 fam.pop    D.PAL    groupIX 18339767      88 0.00000000 1.00000000
#> 2633 fam.pop    D.SMR    groupIX 18339767      88 0.60000000 0.50000000
#> 2634 fam.pop    D.UPD    groupIX 18339767      88 0.53571429 0.37500000
#> 2635     pop      ASP    groupIX 18339767      88 0.53333333 0.50000000
#> 2636     pop      CLF    groupIX 18339767      88 0.40909091 0.25000000
#> 2637     pop      OPL    groupIX 18339767      88 0.45000000 0.32500000
#> 2638     pop      PAL    groupIX 18339767      88 0.00000000 1.00000000
#> 2639     pop      SMR    groupIX 18339767      88 0.49798387 0.40625000
#> 2640     pop      UPD    groupIX 18339767      88 0.46376812 0.33333333
#> 2641 fam.pop    A.ASP  groupXIII 18836558      89 0.00000000 0.00000000
#> 2642 fam.pop    A.CLF  groupXIII 18836558      89        NaN         NA
#> 2643 fam.pop    A.OPL  groupXIII 18836558      89 0.46666667 0.30000000
#> 2644 fam.pop    A.PAL  groupXIII 18836558      89 0.57142857 0.50000000
#> 2645 fam.pop    A.SMR  groupXIII 18836558      89 0.53333333 0.40000000
#> 2646 fam.pop    A.UPD  groupXIII 18836558      89 0.50000000 0.75000000
#> 2647 fam.pop    B.ASP  groupXIII 18836558      89 0.60000000 0.50000000
#> 2648 fam.pop    B.CLF  groupXIII 18836558      89 0.00000000 1.00000000
#> 2649 fam.pop    B.OPL  groupXIII 18836558      89 0.53571429 0.62500000
#> 2650 fam.pop    B.PAL  groupXIII 18836558      89 0.53333333 0.33333333
#> 2651 fam.pop    B.SMR  groupXIII 18836558      89 0.53333333 0.40000000
#> 2652 fam.pop    B.UPD  groupXIII 18836558      89 0.33333333 0.83333333
#> 2653 fam.pop    C.ASP  groupXIII 18836558      89 0.53333333 0.33333333
#> 2654 fam.pop    C.CLF  groupXIII 18836558      89 0.00000000 0.00000000
#> 2655 fam.pop    C.OPL  groupXIII 18836558      89 0.53333333 0.40000000
#> 2656 fam.pop    C.PAL  groupXIII 18836558      89 0.53333333 0.33333333
#> 2657 fam.pop    C.SMR  groupXIII 18836558      89 0.42857143 0.25000000
#> 2658 fam.pop    C.UPD  groupXIII 18836558      89 0.53333333 0.66666667
#> 2659 fam.pop    D.ASP  groupXIII 18836558      89 1.00000000 0.50000000
#> 2660 fam.pop    D.CLF  groupXIII 18836558      89 0.00000000 0.00000000
#> 2661 fam.pop    D.OPL  groupXIII 18836558      89 0.25000000 0.12500000
#> 2662 fam.pop    D.PAL  groupXIII 18836558      89 0.00000000 1.00000000
#> 2663 fam.pop    D.SMR  groupXIII 18836558      89 0.53571429 0.37500000
#> 2664 fam.pop    D.UPD  groupXIII 18836558      89 0.42857143 0.25000000
#> 2665     pop      ASP  groupXIII 18836558      89 0.47058824 0.33333333
#> 2666     pop      CLF  groupXIII 18836558      89 0.53333333 0.33333333
#> 2667     pop      OPL  groupXIII 18836558      89 0.47460317 0.36111111
#> 2668     pop      PAL  groupXIII 18836558      89 0.52173913 0.50000000
#> 2669     pop      SMR  groupXIII 18836558      89 0.47460317 0.36111111
#> 2670     pop      UPD  groupXIII 18836558      89 0.50724638 0.58333333
#> 2671 fam.pop    A.ASP    groupII 19244005      90 0.33333333 0.16666667
#> 2672 fam.pop    A.CLF    groupII 19244005      90 0.00000000 0.00000000
#> 2673 fam.pop    A.OPL    groupII 19244005      90 0.00000000 0.00000000
#> 2674 fam.pop    A.PAL    groupII 19244005      90 0.00000000 0.00000000
#> 2675 fam.pop    A.SMR    groupII 19244005      90 0.00000000 0.00000000
#> 2676 fam.pop    A.UPD    groupII 19244005      90 0.53333333 0.33333333
#> 2677 fam.pop    B.ASP    groupII 19244005      90 0.53333333 0.33333333
#> 2678 fam.pop    B.CLF    groupII 19244005      90 0.00000000 0.00000000
#> 2679 fam.pop    B.OPL    groupII 19244005      90 0.00000000 0.00000000
#> 2680 fam.pop    B.PAL    groupII 19244005      90 0.00000000 0.00000000
#> 2681 fam.pop    B.SMR    groupII 19244005      90 0.00000000 0.00000000
#> 2682 fam.pop    B.UPD    groupII 19244005      90 0.00000000 0.00000000
#> 2683 fam.pop    C.ASP    groupII 19244005      90 0.00000000 0.00000000
#> 2684 fam.pop    C.CLF    groupII 19244005      90 0.00000000 0.00000000
#> 2685 fam.pop    C.OPL    groupII 19244005      90 0.00000000 0.00000000
#> 2686 fam.pop    C.PAL    groupII 19244005      90 0.00000000 0.00000000
#> 2687 fam.pop    C.SMR    groupII 19244005      90 0.00000000 0.00000000
#> 2688 fam.pop    C.UPD    groupII 19244005      90 0.00000000 0.00000000
#> 2689 fam.pop    D.ASP    groupII 19244005      90 0.00000000 0.00000000
#> 2690 fam.pop    D.CLF    groupII 19244005      90 0.00000000 0.00000000
#> 2691 fam.pop    D.OPL    groupII 19244005      90 0.00000000 0.00000000
#> 2692 fam.pop    D.PAL    groupII 19244005      90 0.00000000 0.00000000
#> 2693 fam.pop    D.SMR    groupII 19244005      90 0.16666667 0.08333333
#> 2694 fam.pop    D.UPD    groupII 19244005      90 0.00000000 0.00000000
#> 2695     pop      ASP    groupII 19244005      90 0.24675325 0.13636364
#> 2696     pop      CLF    groupII 19244005      90 0.00000000 0.00000000
#> 2697     pop      OPL    groupII 19244005      90 0.00000000 0.00000000
#> 2698     pop      PAL    groupII 19244005      90 0.00000000 0.00000000
#> 2699     pop      SMR    groupII 19244005      90 0.04347826 0.02173913
#> 2700     pop      UPD    groupII 19244005      90 0.14769231 0.07692308
#> 2701 fam.pop    A.ASP    groupXX 19317273      91 0.00000000 0.00000000
#> 2702 fam.pop    A.CLF    groupXX 19317273      91 1.00000000 0.50000000
#> 2703 fam.pop    A.OPL    groupXX 19317273      91 0.55555556 0.50000000
#> 2704 fam.pop    A.PAL    groupXX 19317273      91 0.25000000 0.12500000
#> 2705 fam.pop    A.SMR    groupXX 19317273      91 0.48484848 0.66666667
#> 2706 fam.pop    A.UPD    groupXX 19317273      91 0.66666667 0.50000000
#> 2707 fam.pop    B.ASP    groupXX 19317273      91 0.50000000 0.25000000
#> 2708 fam.pop    B.CLF    groupXX 19317273      91 0.00000000 1.00000000
#> 2709 fam.pop    B.OPL    groupXX 19317273      91 0.53030303 0.41666667
#> 2710 fam.pop    B.PAL    groupXX 19317273      91 0.25000000 0.12500000
#> 2711 fam.pop    B.SMR    groupXX 19317273      91 0.54545455 0.50000000
#> 2712 fam.pop    B.UPD    groupXX 19317273      91 0.53571429 0.37500000
#> 2713 fam.pop    C.ASP    groupXX 19317273      91 0.00000000 0.00000000
#> 2714 fam.pop    C.CLF    groupXX 19317273      91 0.00000000 0.00000000
#> 2715 fam.pop    C.OPL    groupXX 19317273      91 0.42857143 0.75000000
#> 2716 fam.pop    C.PAL    groupXX 19317273      91 0.00000000 0.00000000
#> 2717 fam.pop    C.SMR    groupXX 19317273      91 0.53333333 0.66666667
#> 2718 fam.pop    C.UPD    groupXX 19317273      91 0.53571429 0.37500000
#> 2719 fam.pop    D.ASP    groupXX 19317273      91 0.66666667 0.50000000
#> 2720 fam.pop    D.CLF    groupXX 19317273      91 0.00000000 0.00000000
#> 2721 fam.pop    D.OPL    groupXX 19317273      91 0.53030303 0.41666667
#> 2722 fam.pop    D.PAL    groupXX 19317273      91 0.25000000 0.12500000
#> 2723 fam.pop    D.SMR    groupXX 19317273      91 0.53333333 0.40000000
#> 2724 fam.pop    D.UPD    groupXX 19317273      91 0.42857143 0.25000000
#> 2725     pop      ASP    groupXX 19317273      91 0.29411765 0.16666667
#> 2726     pop      CLF    groupXX 19317273      91 0.53571429 0.37500000
#> 2727     pop      OPL    groupXX 19317273      91 0.51219512 0.50000000
#> 2728     pop      PAL    groupXX 19317273      91 0.17540323 0.09375000
#> 2729     pop      SMR    groupXX 19317273      91 0.50769231 0.55000000
#> 2730     pop      UPD    groupXX 19317273      91 0.47619048 0.35714286
#> 2731 fam.pop    A.ASP     groupI 19759331      92 0.66666667 0.50000000
#> 2732 fam.pop    A.CLF     groupI 19759331      92 1.00000000 0.50000000
#> 2733 fam.pop    A.OPL     groupI 19759331      92 0.46666667 0.30000000
#> 2734 fam.pop    A.PAL     groupI 19759331      92 0.53333333 0.66666667
#> 2735 fam.pop    A.SMR     groupI 19759331      92 0.00000000 0.00000000
#> 2736 fam.pop    A.UPD     groupI 19759331      92 0.00000000 0.00000000
#> 2737 fam.pop    B.ASP     groupI 19759331      92 0.00000000 1.00000000
#> 2738 fam.pop    B.CLF     groupI 19759331      92 0.00000000 0.00000000
#> 2739 fam.pop    B.OPL     groupI 19759331      92 0.60000000 0.50000000
#> 2740 fam.pop    B.PAL     groupI 19759331      92 0.50000000 0.75000000
#> 2741 fam.pop    B.SMR     groupI 19759331      92 0.53333333 0.33333333
#> 2742 fam.pop    B.UPD     groupI 19759331      92 0.66666667 0.50000000
#> 2743 fam.pop    C.ASP     groupI 19759331      92 0.53333333 0.33333333
#> 2744 fam.pop    C.CLF     groupI 19759331      92 0.00000000 0.00000000
#> 2745 fam.pop    C.OPL     groupI 19759331      92 0.53333333 0.33333333
#> 2746 fam.pop    C.PAL     groupI 19759331      92 0.53571429 0.62500000
#> 2747 fam.pop    C.SMR     groupI 19759331      92 0.00000000 0.00000000
#> 2748 fam.pop    C.UPD     groupI 19759331      92 0.33333333 0.16666667
#> 2749 fam.pop    D.ASP     groupI 19759331      92 1.00000000 0.50000000
#> 2750 fam.pop    D.CLF     groupI 19759331      92 0.50000000 0.25000000
#> 2751 fam.pop    D.OPL     groupI 19759331      92 0.53333333 0.33333333
#> 2752 fam.pop    D.PAL     groupI 19759331      92 0.66666667 0.50000000
#> 2753 fam.pop    D.SMR     groupI 19759331      92 0.53333333 0.33333333
#> 2754 fam.pop    D.UPD     groupI 19759331      92 0.33333333 0.16666667
#> 2755     pop      ASP     groupI 19759331      92 0.53846154 0.50000000
#> 2756     pop      CLF     groupI 19759331      92 0.35555556 0.20000000
#> 2757     pop      OPL     groupI 19759331      92 0.47619048 0.35714286
#> 2758     pop      PAL     groupI 19759331      92 0.48484848 0.63636364
#> 2759     pop      SMR     groupI 19759331      92 0.31168831 0.18181818
#> 2760     pop      UPD     groupI 19759331      92 0.33684211 0.20000000
#> 2761 fam.pop    A.ASP   groupVII 20165240      93 0.33333333 0.16666667
#> 2762 fam.pop    A.CLF   groupVII 20165240      93 0.00000000 0.00000000
#> 2763 fam.pop    A.OPL   groupVII 20165240      93 0.00000000 0.00000000
#> 2764 fam.pop    A.PAL   groupVII 20165240      93 0.33333333 0.16666667
#> 2765 fam.pop    A.SMR   groupVII 20165240      93 0.53571429 0.37500000
#> 2766 fam.pop    A.UPD   groupVII 20165240      93 0.66666667 0.50000000
#> 2767 fam.pop    B.ASP   groupVII 20165240      93 0.33333333 0.16666667
#> 2768 fam.pop    B.CLF   groupVII 20165240      93 1.00000000 0.50000000
#> 2769 fam.pop    B.OPL   groupVII 20165240      93 0.26373626 0.14285714
#> 2770 fam.pop    B.PAL   groupVII 20165240      93 0.42857143 0.25000000
#> 2771 fam.pop    B.SMR   groupVII 20165240      93 0.20000000 0.10000000
#> 2772 fam.pop    B.UPD   groupVII 20165240      93 0.42857143 0.25000000
#> 2773 fam.pop    C.ASP   groupVII 20165240      93 0.33333333 0.16666667
#> 2774 fam.pop    C.CLF   groupVII 20165240      93 0.00000000 0.00000000
#> 2775 fam.pop    C.OPL   groupVII 20165240      93 0.16666667 0.08333333
#> 2776 fam.pop    C.PAL   groupVII 20165240      93 0.00000000 0.00000000
#> 2777 fam.pop    C.SMR   groupVII 20165240      93 0.20000000 0.10000000
#> 2778 fam.pop    C.UPD   groupVII 20165240      93 0.25000000 0.12500000
#> 2779 fam.pop    D.ASP   groupVII 20165240      93 0.00000000 0.00000000
#> 2780 fam.pop    D.CLF   groupVII 20165240      93 0.66666667 0.50000000
#> 2781 fam.pop    D.OPL   groupVII 20165240      93 0.00000000 0.00000000
#> 2782 fam.pop    D.PAL   groupVII 20165240      93 0.25000000 0.12500000
#> 2783 fam.pop    D.SMR   groupVII 20165240      93 0.35555556 0.20000000
#> 2784 fam.pop    D.UPD   groupVII 20165240      93 0.42857143 0.25000000
#> 2785     pop      ASP   groupVII 20165240      93 0.24675325 0.13636364
#> 2786     pop      CLF   groupVII 20165240      93 0.40909091 0.25000000
#> 2787     pop      OPL   groupVII 20165240      93 0.11085973 0.05769231
#> 2788     pop      PAL   groupVII 20165240      93 0.25396825 0.14285714
#> 2789     pop      SMR   groupVII 20165240      93 0.30867710 0.18421053
#> 2790     pop      UPD   groupVII 20165240      93 0.38888889 0.25000000
#> 2791 fam.pop    A.ASP    groupII 20496069      94 0.00000000 0.00000000
#> 2792 fam.pop    A.CLF    groupII 20496069      94 0.00000000 0.00000000
#> 2793 fam.pop    A.OPL    groupII 20496069      94 0.00000000 0.00000000
#> 2794 fam.pop    A.PAL    groupII 20496069      94 0.53333333 0.40000000
#> 2795 fam.pop    A.SMR    groupII 20496069      94 0.00000000 0.00000000
#> 2796 fam.pop    A.UPD    groupII 20496069      94 0.00000000 0.00000000
#> 2797 fam.pop    B.ASP    groupII 20496069      94 0.00000000 0.00000000
#> 2798 fam.pop    B.CLF    groupII 20496069      94 0.00000000 0.00000000
#> 2799 fam.pop    B.OPL    groupII 20496069      94 0.20000000 0.10000000
#> 2800 fam.pop    B.PAL    groupII 20496069      94 0.53333333 0.33333333
#> 2801 fam.pop    B.SMR    groupII 20496069      94 0.16666667 0.08333333
#> 2802 fam.pop    B.UPD    groupII 20496069      94 0.25000000 0.12500000
#> 2803 fam.pop    C.ASP    groupII 20496069      94 0.00000000 0.00000000
#> 2804 fam.pop    C.CLF    groupII 20496069      94 0.00000000 0.00000000
#> 2805 fam.pop    C.OPL    groupII 20496069      94 0.16666667 0.08333333
#> 2806 fam.pop    C.PAL    groupII 20496069      94 0.33333333 0.16666667
#> 2807 fam.pop    C.SMR    groupII 20496069      94 0.00000000 0.00000000
#> 2808 fam.pop    C.UPD    groupII 20496069      94 0.00000000 0.00000000
#> 2809 fam.pop    D.ASP    groupII 20496069      94 0.00000000 0.00000000
#> 2810 fam.pop    D.CLF    groupII 20496069      94 0.00000000 0.00000000
#> 2811 fam.pop    D.OPL    groupII 20496069      94 0.16666667 0.08333333
#> 2812 fam.pop    D.PAL    groupII 20496069      94 0.66666667 0.50000000
#> 2813 fam.pop    D.SMR    groupII 20496069      94 0.00000000 0.00000000
#> 2814 fam.pop    D.UPD    groupII 20496069      94 0.42857143 0.25000000
#> 2815     pop      ASP    groupII 20496069      94 0.00000000 0.00000000
#> 2816     pop      CLF    groupII 20496069      94 0.00000000 0.00000000
#> 2817     pop      OPL    groupII 20496069      94 0.13002114 0.06818182
#> 2818     pop      PAL    groupII 20496069      94 0.47076923 0.34615385
#> 2819     pop      SMR    groupII 20496069      94 0.04545455 0.02272727
#> 2820     pop      UPD    groupII 20496069      94 0.19841270 0.10714286
#> 2821 fam.pop    A.ASP     groupI 21284266      95 0.33333333 0.16666667
#> 2822 fam.pop    A.CLF     groupI 21284266      95 0.66666667 0.50000000
#> 2823 fam.pop    A.OPL     groupI 21284266      95 0.16666667 0.08333333
#> 2824 fam.pop    A.PAL     groupI 21284266      95 0.00000000 0.00000000
#> 2825 fam.pop    A.SMR     groupI 21284266      95 0.20000000 0.10000000
#> 2826 fam.pop    A.UPD     groupI 21284266      95 0.33333333 0.16666667
#> 2827 fam.pop    B.ASP     groupI 21284266      95 0.33333333 0.16666667
#> 2828 fam.pop    B.CLF     groupI 21284266      95 0.00000000 0.00000000
#> 2829 fam.pop    B.OPL     groupI 21284266      95 0.40909091 0.25000000
#> 2830 fam.pop    B.PAL     groupI 21284266      95 0.00000000 0.00000000
#> 2831 fam.pop    B.SMR     groupI 21284266      95 0.20000000 0.10000000
#> 2832 fam.pop    B.UPD     groupI 21284266      95 0.00000000 0.00000000
#> 2833 fam.pop    C.ASP     groupI 21284266      95 0.33333333 0.16666667
#> 2834 fam.pop    C.CLF     groupI 21284266      95 0.00000000 0.00000000
#> 2835 fam.pop    C.OPL     groupI 21284266      95 0.46666667 0.30000000
#> 2836 fam.pop    C.PAL     groupI 21284266      95 0.00000000 0.00000000
#> 2837 fam.pop    C.SMR     groupI 21284266      95 0.35555556 0.20000000
#> 2838 fam.pop    C.UPD     groupI 21284266      95 0.53333333 0.33333333
#> 2839 fam.pop    D.ASP     groupI 21284266      95 1.00000000 0.50000000
#> 2840 fam.pop    D.CLF     groupI 21284266      95 0.66666667 0.50000000
#> 2841 fam.pop    D.OPL     groupI 21284266      95 0.20000000 0.10000000
#> 2842 fam.pop    D.PAL     groupI 21284266      95 0.00000000 0.00000000
#> 2843 fam.pop    D.SMR     groupI 21284266      95 0.42857143 0.25000000
#> 2844 fam.pop    D.UPD     groupI 21284266      95 0.33333333 0.16666667
#> 2845     pop      ASP     groupI 21284266      95 0.33684211 0.20000000
#> 2846     pop      CLF     groupI 21284266      95 0.48484848 0.33333333
#> 2847     pop      OPL     groupI 21284266      95 0.30443975 0.18181818
#> 2848     pop      PAL     groupI 21284266      95 0.00000000 0.00000000
#> 2849     pop      SMR     groupI 21284266      95 0.27311522 0.15789474
#> 2850     pop      UPD     groupI 21284266      95 0.28985507 0.16666667
#> 2851 fam.pop    A.ASP    groupIV 22227311      96 0.33333333 0.16666667
#> 2852 fam.pop    A.CLF    groupIV 22227311      96 0.00000000 0.00000000
#> 2853 fam.pop    A.OPL    groupIV 22227311      96 0.42857143 0.25000000
#> 2854 fam.pop    A.PAL    groupIV 22227311      96 0.00000000 0.00000000
#> 2855 fam.pop    A.SMR    groupIV 22227311      96 0.16666667 0.08333333
#> 2856 fam.pop    A.UPD    groupIV 22227311      96 0.00000000 0.00000000
#> 2857 fam.pop    B.ASP    groupIV 22227311      96 0.00000000 0.00000000
#> 2858 fam.pop    B.CLF    groupIV 22227311      96 0.00000000 0.00000000
#> 2859 fam.pop    B.OPL    groupIV 22227311      96 0.40909091 0.25000000
#> 2860 fam.pop    B.PAL    groupIV 22227311      96 0.00000000 0.00000000
#> 2861 fam.pop    B.SMR    groupIV 22227311      96 0.16666667 0.08333333
#> 2862 fam.pop    B.UPD    groupIV 22227311      96 0.00000000 0.00000000
#> 2863 fam.pop    C.ASP    groupIV 22227311      96 0.00000000 0.00000000
#> 2864 fam.pop    C.CLF    groupIV 22227311      96 1.00000000 0.50000000
#> 2865 fam.pop    C.OPL    groupIV 22227311      96 0.30303030 0.16666667
#> 2866 fam.pop    C.PAL    groupIV 22227311      96 0.00000000 0.00000000
#> 2867 fam.pop    C.SMR    groupIV 22227311      96 0.20000000 0.10000000
#> 2868 fam.pop    C.UPD    groupIV 22227311      96 0.42857143 0.25000000
#> 2869 fam.pop    D.ASP    groupIV 22227311      96 0.00000000 0.00000000
#> 2870 fam.pop    D.CLF    groupIV 22227311      96 0.00000000 0.00000000
#> 2871 fam.pop    D.OPL    groupIV 22227311      96 0.14285714 0.07142857
#> 2872 fam.pop    D.PAL    groupIV 22227311      96 0.00000000 0.00000000
#> 2873 fam.pop    D.SMR    groupIV 22227311      96 0.16666667 0.08333333
#> 2874 fam.pop    D.UPD    groupIV 22227311      96 0.00000000 0.00000000
#> 2875     pop      ASP    groupIV 22227311      96 0.12500000 0.06250000
#> 2876     pop      CLF    groupIV 22227311      96 0.16666667 0.08333333
#> 2877     pop      OPL    groupIV 22227311      96 0.29371981 0.17391304
#> 2878     pop      PAL    groupIV 22227311      96 0.00000000 0.00000000
#> 2879     pop      SMR    groupIV 22227311      96 0.16231884 0.08695652
#> 2880     pop      UPD    groupIV 22227311      96 0.14769231 0.07692308
#> 2881 fam.pop    A.ASP     groupI 22746027      97 0.00000000 0.00000000
#> 2882 fam.pop    A.CLF     groupI 22746027      97 0.00000000 0.00000000
#> 2883 fam.pop    A.OPL     groupI 22746027      97 0.20000000 0.10000000
#> 2884 fam.pop    A.PAL     groupI 22746027      97 0.00000000 0.00000000
#> 2885 fam.pop    A.SMR     groupI 22746027      97 0.16666667 0.08333333
#> 2886 fam.pop    A.UPD     groupI 22746027      97 0.33333333 0.16666667
#> 2887 fam.pop    B.ASP     groupI 22746027      97 0.00000000 0.00000000
#> 2888 fam.pop    B.CLF     groupI 22746027      97 0.00000000 0.00000000
#> 2889 fam.pop    B.OPL     groupI 22746027      97 0.25000000 0.12500000
#> 2890 fam.pop    B.PAL     groupI 22746027      97 0.00000000 0.00000000
#> 2891 fam.pop    B.SMR     groupI 22746027      97 0.25000000 0.12500000
#> 2892 fam.pop    B.UPD     groupI 22746027      97 0.00000000 0.00000000
#> 2893 fam.pop    C.ASP     groupI 22746027      97 0.00000000 0.00000000
#> 2894 fam.pop    C.CLF     groupI 22746027      97 0.00000000 0.00000000
#> 2895 fam.pop    C.OPL     groupI 22746027      97 0.00000000 0.00000000
#> 2896 fam.pop    C.PAL     groupI 22746027      97 0.00000000 0.00000000
#> 2897 fam.pop    C.SMR     groupI 22746027      97 0.00000000 0.00000000
#> 2898 fam.pop    C.UPD     groupI 22746027      97 0.00000000 0.00000000
#> 2899 fam.pop    D.ASP     groupI 22746027      97 0.00000000 0.00000000
#> 2900 fam.pop    D.CLF     groupI 22746027      97 0.50000000 0.25000000
#> 2901 fam.pop    D.OPL     groupI 22746027      97 0.00000000 0.00000000
#> 2902 fam.pop    D.PAL     groupI 22746027      97 0.00000000 0.00000000
#> 2903 fam.pop    D.SMR     groupI 22746027      97 0.00000000 0.00000000
#> 2904 fam.pop    D.UPD     groupI 22746027      97 0.00000000 0.00000000
#> 2905     pop      ASP     groupI 22746027      97 0.00000000 0.00000000
#> 2906     pop      CLF     groupI 22746027      97 0.16666667 0.08333333
#> 2907     pop      OPL     groupI 22746027      97 0.09291521 0.04761905
#> 2908     pop      PAL     groupI 22746027      97 0.00000000 0.00000000
#> 2909     pop      SMR     groupI 22746027      97 0.09743590 0.05000000
#> 2910     pop      UPD     groupI 22746027      97 0.06666667 0.03333333
#> 2911 fam.pop    A.ASP   groupVII 25239495      98 0.00000000 0.00000000
#> 2912 fam.pop    A.CLF   groupVII 25239495      98 0.66666667 0.50000000
#> 2913 fam.pop    A.OPL   groupVII 25239495      98 0.16666667 0.08333333
#> 2914 fam.pop    A.PAL   groupVII 25239495      98 0.55555556 0.50000000
#> 2915 fam.pop    A.SMR   groupVII 25239495      98 0.48484848 0.33333333
#> 2916 fam.pop    A.UPD   groupVII 25239495      98 0.33333333 0.16666667
#> 2917 fam.pop    B.ASP   groupVII 25239495      98 0.33333333 0.16666667
#> 2918 fam.pop    B.CLF   groupVII 25239495      98 0.00000000 0.00000000
#> 2919 fam.pop    B.OPL   groupVII 25239495      98 0.14285714 0.07142857
#> 2920 fam.pop    B.PAL   groupVII 25239495      98 0.53571429 0.37500000
#> 2921 fam.pop    B.SMR   groupVII 25239495      98 0.20000000 0.10000000
#> 2922 fam.pop    B.UPD   groupVII 25239495      98 0.57142857 0.50000000
#> 2923 fam.pop    C.ASP   groupVII 25239495      98 0.53333333 0.33333333
#> 2924 fam.pop    C.CLF   groupVII 25239495      98        NaN         NA
#> 2925 fam.pop    C.OPL   groupVII 25239495      98 0.16666667 0.08333333
#> 2926 fam.pop    C.PAL   groupVII 25239495      98 0.53571429 0.37500000
#> 2927 fam.pop    C.SMR   groupVII 25239495      98 0.00000000 0.00000000
#> 2928 fam.pop    C.UPD   groupVII 25239495      98 0.50000000 0.25000000
#> 2929 fam.pop    D.ASP   groupVII 25239495      98 1.00000000 0.50000000
#> 2930 fam.pop    D.CLF   groupVII 25239495      98 0.50000000 0.25000000
#> 2931 fam.pop    D.OPL   groupVII 25239495      98 0.26373626 0.14285714
#> 2932 fam.pop    D.PAL   groupVII 25239495      98 0.53333333 0.33333333
#> 2933 fam.pop    D.SMR   groupVII 25239495      98 0.48484848 0.33333333
#> 2934 fam.pop    D.UPD   groupVII 25239495      98 0.53571429 0.37500000
#> 2935     pop      ASP   groupVII 25239495      98 0.36601307 0.22222222
#> 2936     pop      CLF   groupVII 25239495      98 0.46666667 0.30000000
#> 2937     pop      OPL   groupVII 25239495      98 0.17722474 0.09615385
#> 2938     pop      PAL   groupVII 25239495      98 0.49798387 0.40625000
#> 2939     pop      SMR   groupVII 25239495      98 0.33298097 0.20454545
#> 2940     pop      UPD   groupVII 25239495      98 0.47076923 0.34615385
#> 2941 fam.pop    A.ASP   groupVII 27599144      99 0.50000000 0.75000000
#> 2942 fam.pop    A.CLF   groupVII 27599144      99 0.66666667 0.50000000
#> 2943 fam.pop    A.OPL   groupVII 27599144      99 0.25000000 0.87500000
#> 2944 fam.pop    A.PAL   groupVII 27599144      99 0.00000000 0.00000000
#> 2945 fam.pop    A.SMR   groupVII 27599144      99 0.54545455 0.50000000
#> 2946 fam.pop    A.UPD   groupVII 27599144      99 0.50000000 0.75000000
#> 2947 fam.pop    B.ASP   groupVII 27599144      99 0.33333333 0.16666667
#> 2948 fam.pop    B.CLF   groupVII 27599144      99 0.00000000 0.00000000
#> 2949 fam.pop    B.OPL   groupVII 27599144      99 0.55555556 0.50000000
#> 2950 fam.pop    B.PAL   groupVII 27599144      99 0.00000000 0.00000000
#> 2951 fam.pop    B.SMR   groupVII 27599144      99 0.53571429 0.62500000
#> 2952 fam.pop    B.UPD   groupVII 27599144      99 0.60000000 0.50000000
#> 2953 fam.pop    C.ASP   groupVII 27599144      99 0.53333333 0.66666667
#> 2954 fam.pop    C.CLF   groupVII 27599144      99 0.00000000 1.00000000
#> 2955 fam.pop    C.OPL   groupVII 27599144      99 0.54545455 0.50000000
#> 2956 fam.pop    C.PAL   groupVII 27599144      99 0.00000000 0.00000000
#> 2957 fam.pop    C.SMR   groupVII 27599144      99 0.53333333 0.60000000
#> 2958 fam.pop    C.UPD   groupVII 27599144      99 0.33333333 0.16666667
#> 2959 fam.pop    D.ASP   groupVII 27599144      99 0.50000000 0.75000000
#> 2960 fam.pop    D.CLF   groupVII 27599144      99 0.50000000 0.25000000
#> 2961 fam.pop    D.OPL   groupVII 27599144      99 0.54545455 0.50000000
#> 2962 fam.pop    D.PAL   groupVII 27599144      99 0.00000000 0.00000000
#> 2963 fam.pop    D.SMR   groupVII 27599144      99 0.53333333 0.40000000
#> 2964 fam.pop    D.UPD   groupVII 27599144      99 0.53333333 0.66666667
#> 2965     pop      ASP   groupVII 27599144      99 0.52105263 0.55000000
#> 2966     pop      CLF   groupVII 27599144      99 0.53030303 0.41666667
#> 2967     pop      OPL   groupVII 27599144      99 0.50174216 0.57142857
#> 2968     pop      PAL   groupVII 27599144      99 0.00000000 0.00000000
#> 2969     pop      SMR   groupVII 27599144      99 0.51153846 0.52500000
#> 2970     pop      UPD   groupVII 27599144      99 0.52380952 0.50000000
#> 2971 fam.pop    A.ASP    groupIV 32321895     100 0.33333333 0.83333333
#> 2972 fam.pop    A.CLF    groupIV 32321895     100 0.66666667 0.50000000
#> 2973 fam.pop    A.OPL    groupIV 32321895     100 0.35555556 0.20000000
#> 2974 fam.pop    A.PAL    groupIV 32321895     100 0.53333333 0.33333333
#> 2975 fam.pop    A.SMR    groupIV 32321895     100 0.40909091 0.25000000
#> 2976 fam.pop    A.UPD    groupIV 32321895     100 0.33333333 0.83333333
#> 2977 fam.pop    B.ASP    groupIV 32321895     100 0.50000000 0.25000000
#> 2978 fam.pop    B.CLF    groupIV 32321895     100        NaN         NA
#> 2979 fam.pop    B.OPL    groupIV 32321895     100 0.35555556 0.20000000
#> 2980 fam.pop    B.PAL    groupIV 32321895     100 0.53333333 0.33333333
#> 2981 fam.pop    B.SMR    groupIV 32321895     100 0.42857143 0.25000000
#> 2982 fam.pop    B.UPD    groupIV 32321895     100 0.57142857 0.50000000
#> 2983 fam.pop    C.ASP    groupIV 32321895     100 0.60000000 0.50000000
#> 2984 fam.pop    C.CLF    groupIV 32321895     100 0.00000000 0.00000000
#> 2985 fam.pop    C.OPL    groupIV 32321895     100 0.35555556 0.20000000
#> 2986 fam.pop    C.PAL    groupIV 32321895     100 0.53333333 0.33333333
#> 2987 fam.pop    C.SMR    groupIV 32321895     100 0.53333333 0.40000000
#> 2988 fam.pop    C.UPD    groupIV 32321895     100 0.60000000 0.50000000
#> 2989 fam.pop    D.ASP    groupIV 32321895     100 0.50000000 0.75000000
#> 2990 fam.pop    D.CLF    groupIV 32321895     100 1.00000000 0.50000000
#> 2991 fam.pop    D.OPL    groupIV 32321895     100 0.53333333 0.60000000
#> 2992 fam.pop    D.PAL    groupIV 32321895     100 0.50000000 0.25000000
#> 2993 fam.pop    D.SMR    groupIV 32321895     100 0.53333333 0.40000000
#> 2994 fam.pop    D.UPD    groupIV 32321895     100 0.53571429 0.62500000
#> 2995     pop      ASP    groupIV 32321895     100 0.50526316 0.60000000
#> 2996     pop      CLF    groupIV 32321895     100 0.53571429 0.37500000
#> 2997     pop      OPL    groupIV 32321895     100 0.43076923 0.30000000
#> 2998     pop      PAL    groupIV 32321895     100 0.45454545 0.31818182
#> 2999     pop      SMR    groupIV 32321895     100 0.45000000 0.32500000
#> 3000     pop      UPD    groupIV 32321895     100 0.49470899 0.60714286
#>      maj.count major min.count minor
#> 1            1     G         3     A
#> 2            1     G         3     A
#> 3            5     G         7     A
#> 4            6     G         4     A
#> 5            6     G         6     A
#> 6            4     G         2     A
#> 7            5     G         1     A
#> 8            1     G         1     A
#> 9            4     G         6     A
#> 10           8     G         0     A
#> 11           3     G         7     A
#> 12           3     G         3     A
#> 13           1     G         3     A
#> 14           0     G         2     A
#> 15           3     G         5     A
#> 16           6     G         2     A
#> 17           6     G         4     A
#> 18           3     G         5     A
#> 19           2     G         2     A
#> 20           1     G         1     A
#> 21           5     G         5     A
#> 22           7     G         1     A
#> 23           3     G         9     A
#> 24           3     G         3     A
#> 25           9     G         9     A
#> 26           3     G         7     A
#> 27          17     G        23     A
#> 28          27     G         7     A
#> 29          18     G        26     A
#> 30          13     G        13     A
#> 31           5     G         1     T
#> 32           4     G         0     T
#> 33          12     G         0     T
#> 34          10     G         0     T
#> 35          12     G         0     T
#> 36           6     G         0     T
#> 37           6     G         0     T
#> 38           2     G         0     T
#> 39           9     G         1     T
#> 40           8     G         0     T
#> 41           9     G         1     T
#> 42           8     G         0     T
#> 43           6     G         0     T
#> 44           2     G         0     T
#> 45          12     G         0     T
#> 46           8     G         0     T
#> 47          11     G         1     T
#> 48           8     G         0     T
#> 49           2     G         0     T
#> 50           4     G         0     T
#> 51          10     G         2     T
#> 52           8     G         0     T
#> 53          10     G         0     T
#> 54           8     G         0     T
#> 55          19     G         1     T
#> 56          12     G         0     T
#> 57          43     G         3     T
#> 58          34     G         0     T
#> 59          42     G         2     T
#> 60          30     G         0     T
#> 61           6     C         0     T
#> 62           3     C         1     T
#> 63          10     C         0     T
#> 64           8     C         0     T
#> 65          12     C         0     T
#> 66           6     C         0     T
#> 67           5     C         1     T
#> 68           2     C         0     T
#> 69          13     C         1     T
#> 70           7     C         1     T
#> 71          12     C         0     T
#> 72           7     C         1     T
#> 73           5     C         1     T
#> 74           2     C         0     T
#> 75          10     C         2     T
#> 76           8     C         0     T
#> 77          10     C         2     T
#> 78           7     C         1     T
#> 79           4     C         0     T
#> 80           4     C         0     T
#> 81          13     C         1     T
#> 82           8     C         0     T
#> 83           9     C         1     T
#> 84           8     C         0     T
#> 85          20     C         2     T
#> 86          11     C         1     T
#> 87          46     C         4     T
#> 88          31     C         1     T
#> 89          43     C         3     T
#> 90          28     C         2     T
#> 91           5     A         1     G
#> 92           3     A         1     G
#> 93           4     A         4     G
#> 94          10     A         0     G
#> 95           9     A         3     G
#> 96           2     A         2     G
#> 97           5     A         1     G
#> 98           2     A         0     G
#> 99           7     A         3     G
#> 100          8     A         0     G
#> 101          5     A         5     G
#> 102          4     A         2     G
#> 103          5     A         1     G
#> 104          0     A         0     G
#> 105          6     A         6     G
#> 106          8     A         0     G
#> 107          9     A         3     G
#> 108          4     A         4     G
#> 109          4     A         0     G
#> 110          4     A         0     G
#> 111          6     A         4     G
#> 112          6     A         0     G
#> 113          8     A         4     G
#> 114          5     A         3     G
#> 115         19     A         3     G
#> 116          9     A         1     G
#> 117         23     A        17     G
#> 118         32     A         0     G
#> 119         31     A        15     G
#> 120         15     A        11     G
#> 121          4     G         0     A
#> 122          2     G         0     A
#> 123          7     G         1     A
#> 124          6     G         2     A
#> 125          5     G         3     A
#> 126          1     G         1     A
#> 127          3     G         1     A
#> 128          0     G         0     A
#> 129          2     G         4     A
#> 130          3     G         1     A
#> 131          5     G         3     A
#> 132          5     G         3     A
#> 133          1     G         1     A
#> 134          0     G         0     A
#> 135          5     G         1     A
#> 136          4     G         2     A
#> 137          7     G         1     A
#> 138          1     G         1     A
#> 139          1     G         3     A
#> 140          4     G         0     A
#> 141          4     G         2     A
#> 142          6     G         0     A
#> 143          6     G         0     A
#> 144          3     G         3     A
#> 145          9     G         5     A
#> 146          6     G         0     A
#> 147         18     G         8     A
#> 148         19     G         5     A
#> 149         23     G         7     A
#> 150         10     G         8     A
#> 151          5     A         1     C
#> 152          3     A         1     C
#> 153         11     A         1     C
#> 154          7     A         3     C
#> 155         10     A         0     C
#> 156          4     A         2     C
#> 157          3     A         3     C
#> 158          1     A         1     C
#> 159         13     A         1     C
#> 160          3     A         5     C
#> 161          8     A         2     C
#> 162          4     A         4     C
#> 163          4     A         0     C
#> 164          2     A         0     C
#> 165         14     A         0     C
#> 166          5     A         3     C
#> 167         12     A         0     C
#> 168          5     A         3     C
#> 169          3     A         1     C
#> 170          3     A         1     C
#> 171         11     A         1     C
#> 172          4     A         2     C
#> 173          9     A         1     C
#> 174          5     A         1     C
#> 175         15     A         5     C
#> 176          9     A         3     C
#> 177         49     A         3     C
#> 178         19     A        13     C
#> 179         39     A         3     C
#> 180         18     A        10     C
#> 181          5     T         1     G
#> 182          4     T         0     G
#> 183          9     T         1     G
#> 184         10     T         0     G
#> 185          7     T         3     G
#> 186          4     T         0     G
#> 187          6     T         0     G
#> 188          2     T         0     G
#> 189         10     T         2     G
#> 190          8     T         0     G
#> 191          8     T         2     G
#> 192          8     T         0     G
#> 193          6     T         0     G
#> 194          2     T         0     G
#> 195          9     T         1     G
#> 196          8     T         0     G
#> 197          8     T         2     G
#> 198          8     T         0     G
#> 199          4     T         0     G
#> 200          2     T         2     G
#> 201          8     T         2     G
#> 202          8     T         0     G
#> 203         11     T         1     G
#> 204          8     T         0     G
#> 205         21     T         1     G
#> 206         10     T         2     G
#> 207         36     T         6     G
#> 208         34     T         0     G
#> 209         34     T         8     G
#> 210         28     T         0     G
#> 211          5     G         1     A
#> 212          3     G         1     A
#> 213          7     G         1     A
#> 214          8     G         0     A
#> 215          5     G         1     A
#> 216          1     G         1     A
#> 217          5     G         1     A
#> 218          1     G         1     A
#> 219         12     G         2     A
#> 220          4     G         0     A
#> 221          9     G         3     A
#> 222          3     G         1     A
#> 223          5     G         1     A
#> 224          0     G         2     A
#> 225         13     G         1     A
#> 226          8     G         0     A
#> 227          7     G         3     A
#> 228          7     G         1     A
#> 229          2     G         2     A
#> 230          3     G         1     A
#> 231          8     G         4     A
#> 232          8     G         0     A
#> 233          5     G         7     A
#> 234          4     G         2     A
#> 235         17     G         5     A
#> 236          7     G         5     A
#> 237         40     G         8     A
#> 238         28     G         0     A
#> 239         26     G        14     A
#> 240         15     G         5     A
#> 241          5     C         1     A
#> 242          3     C         1     A
#> 243         11     C         1     A
#> 244         10     C         0     A
#> 245          6     C         2     A
#> 246          4     C         0     A
#> 247          5     C         1     A
#> 248          2     C         0     A
#> 249         13     C         1     A
#> 250          8     C         0     A
#> 251          7     C         3     A
#> 252          5     C         1     A
#> 253          4     C         2     A
#> 254          1     C         1     A
#> 255         10     C         2     A
#> 256          8     C         0     A
#> 257          9     C         3     A
#> 258          6     C         2     A
#> 259          3     C         1     A
#> 260          4     C         0     A
#> 261         11     C         1     A
#> 262          5     C         1     A
#> 263          8     C         2     A
#> 264          4     C         4     A
#> 265         17     C         5     A
#> 266         10     C         2     A
#> 267         45     C         5     A
#> 268         31     C         1     A
#> 269         30     C        10     A
#> 270         19     C         7     A
#> 271          5     A         1     G
#> 272          4     A         0     G
#> 273         12     A         0     G
#> 274         10     A         0     G
#> 275         10     A         0     G
#> 276          6     A         0     G
#> 277          6     A         0     G
#> 278          2     A         0     G
#> 279         12     A         0     G
#> 280          8     A         0     G
#> 281         10     A         0     G
#> 282          6     A         2     G
#> 283          5     A         1     G
#> 284          2     A         0     G
#> 285         11     A         1     G
#> 286          8     A         0     G
#> 287         11     A         1     G
#> 288          6     A         0     G
#> 289          4     A         0     G
#> 290          3     A         1     G
#> 291         14     A         0     G
#> 292          8     A         0     G
#> 293         11     A         1     G
#> 294          6     A         2     G
#> 295         20     A         2     G
#> 296         11     A         1     G
#> 297         49     A         1     G
#> 298         34     A         0     G
#> 299         42     A         2     G
#> 300         24     A         4     G
#> 301          3     T         1     C
#> 302          1     T         3     C
#> 303          4     T         6     C
#> 304         10     T         0     C
#> 305          7     T         3     C
#> 306          4     T         2     C
#> 307          4     T         2     C
#> 308          1     T         1     C
#> 309          4     T         6     C
#> 310          8     T         0     C
#> 311          4     T         6     C
#> 312          4     T         4     C
#> 313          3     T         3     C
#> 314          1     T         1     C
#> 315          4     T         8     C
#> 316          8     T         0     C
#> 317          4     T         6     C
#> 318          3     T         3     C
#> 319          1     T         3     C
#> 320          3     T         1     C
#> 321          5     T         5     C
#> 322          6     T         0     C
#> 323          4     T         4     C
#> 324          6     T         2     C
#> 325         11     T         9     C
#> 326          6     T         6     C
#> 327         17     T        25     C
#> 328         32     T         0     C
#> 329         19     T        19     C
#> 330         17     T        11     C
#> 331          4     A         2     G
#> 332          4     A         0     G
#> 333          8     A         0     G
#> 334          7     A         3     G
#> 335          8     A         2     G
#> 336          3     A         3     G
#> 337          5     A         1     G
#> 338          1     A         1     G
#> 339          6     A         2     G
#> 340          5     A         3     G
#> 341          5     A         1     G
#> 342          5     A         1     G
#> 343          4     A         2     G
#> 344          2     A         0     G
#> 345          7     A         3     G
#> 346          4     A         4     G
#> 347          8     A         2     G
#> 348          3     A         1     G
#> 349          2     A         0     G
#> 350          4     A         0     G
#> 351          8     A         2     G
#> 352          6     A         0     G
#> 353          7     A         1     G
#> 354          3     A         1     G
#> 355         15     A         5     G
#> 356         11     A         1     G
#> 357         29     A         7     G
#> 358         22     A        10     G
#> 359         28     A         6     G
#> 360         14     A         6     G
#> 361          3     T         1     C
#> 362          1     T         3     C
#> 363          8     T         2     C
#> 364          9     T         1     C
#> 365          8     T         4     C
#> 366          6     T         0     C
#> 367          4     T         0     C
#> 368          2     T         0     C
#> 369          7     T         5     C
#> 370          4     T         0     C
#> 371          7     T         1     C
#> 372          7     T         1     C
#> 373          5     T         1     C
#> 374          1     T         1     C
#> 375         13     T         1     C
#> 376          8     T         0     C
#> 377          6     T         2     C
#> 378          4     T         2     C
#> 379          4     T         0     C
#> 380          3     T         1     C
#> 381          8     T         4     C
#> 382          8     T         0     C
#> 383         11     T         1     C
#> 384          5     T         1     C
#> 385         16     T         2     C
#> 386          7     T         5     C
#> 387         36     T        12     C
#> 388         29     T         1     C
#> 389         32     T         8     C
#> 390         22     T         4     C
#> 391          5     G         1     T
#> 392          4     G         0     T
#> 393         10     G         2     T
#> 394         10     G         0     T
#> 395          8     G         2     T
#> 396          4     G         2     T
#> 397          5     G         1     T
#> 398          2     G         0     T
#> 399         14     G         0     T
#> 400          8     G         0     T
#> 401          9     G         1     T
#> 402          6     G         2     T
#> 403          4     G         2     T
#> 404          1     G         1     T
#> 405         10     G         0     T
#> 406          8     G         0     T
#> 407         11     G         1     T
#> 408          5     G         3     T
#> 409          2     G         0     T
#> 410          3     G         1     T
#> 411         12     G         2     T
#> 412          8     G         0     T
#> 413          8     G         4     T
#> 414          4     G         4     T
#> 415         16     G         4     T
#> 416         10     G         2     T
#> 417         46     G         4     T
#> 418         34     G         0     T
#> 419         36     G         8     T
#> 420         19     G        11     T
#> 421          5     A         1     G
#> 422          4     A         0     G
#> 423         12     A         0     G
#> 424          1     A         7     G
#> 425          9     A         3     G
#> 426          1     A         1     G
#> 427          4     A         2     G
#> 428          2     A         0     G
#> 429          6     A         4     G
#> 430          0     A         8     G
#> 431          9     A         1     G
#> 432          4     A         4     G
#> 433          4     A         2     G
#> 434          2     A         0     G
#> 435          7     A         3     G
#> 436          2     A         6     G
#> 437          9     A         3     G
#> 438          5     A         3     G
#> 439          4     A         0     G
#> 440          4     A         0     G
#> 441          7     A         1     G
#> 442          1     A         7     G
#> 443          6     A         2     G
#> 444          4     A         4     G
#> 445         17     A         5     G
#> 446         12     A         0     G
#> 447         32     A         8     G
#> 448          4     A        28     G
#> 449         33     A         9     G
#> 450         14     A        12     G
#> 451          4     G         0     T
#> 452          2     G         2     T
#> 453          7     G         3     T
#> 454          6     G         4     T
#> 455          8     G         4     T
#> 456          5     G         1     T
#> 457          5     G         1     T
#> 458          2     G         0     T
#> 459          8     G         2     T
#> 460          4     G         2     T
#> 461          5     G         3     T
#> 462          4     G         2     T
#> 463          4     G         2     T
#> 464          0     G         0     T
#> 465          8     G         4     T
#> 466          3     G         5     T
#> 467         10     G         2     T
#> 468          5     G         1     T
#> 469          2     G         2     T
#> 470          3     G         1     T
#> 471          4     G         4     T
#> 472          4     G         2     T
#> 473          4     G         4     T
#> 474          5     G         1     T
#> 475         15     G         5     T
#> 476          7     G         3     T
#> 477         27     G        13     T
#> 478         17     G        13     T
#> 479         27     G        13     T
#> 480         19     G         5     T
#> 481          6     A         0     C
#> 482          4     A         0     C
#> 483          8     A         2     C
#> 484          6     A         2     C
#> 485          9     A         3     C
#> 486          5     A         1     C
#> 487          3     A         1     C
#> 488          2     A         0     C
#> 489         11     A         3     C
#> 490          6     A         2     C
#> 491          8     A         2     C
#> 492          6     A         2     C
#> 493          5     A         1     C
#> 494          2     A         0     C
#> 495         12     A         0     C
#> 496          7     A         1     C
#> 497         10     A         0     C
#> 498          5     A         1     C
#> 499          3     A         1     C
#> 500          4     A         0     C
#> 501         12     A         2     C
#> 502          7     A         1     C
#> 503          8     A         2     C
#> 504          6     A         0     C
#> 505         17     A         3     C
#> 506         12     A         0     C
#> 507         43     A         7     C
#> 508         26     A         6     C
#> 509         35     A         7     C
#> 510         22     A         4     C
#> 511          3     C         3     A
#> 512          2     C         2     A
#> 513          9     C         3     A
#> 514          2     C         4     A
#> 515          9     C         3     A
#> 516          5     C         1     A
#> 517          3     C         1     A
#> 518          1     C         1     A
#> 519         11     C         3     A
#> 520          3     C         5     A
#> 521          8     C         2     A
#> 522          3     C         5     A
#> 523          5     C         1     A
#> 524          0     C         0     A
#> 525         12     C         2     A
#> 526          2     C         4     A
#> 527         11     C         1     A
#> 528          6     C         2     A
#> 529          3     C         1     A
#> 530          3     C         1     A
#> 531         11     C         3     A
#> 532          1     C         7     A
#> 533         11     C         1     A
#> 534          3     C         5     A
#> 535         14     C         6     A
#> 536          6     C         4     A
#> 537         43     C        11     A
#> 538          8     C        20     A
#> 539         39     C         7     A
#> 540         17     C        13     A
#> 541          6     T         0     G
#> 542          3     T         1     G
#> 543          9     T         1     G
#> 544         10     T         0     G
#> 545         10     T         2     G
#> 546          6     T         0     G
#> 547          6     T         0     G
#> 548          2     T         0     G
#> 549          5     T         1     G
#> 550          8     T         0     G
#> 551         12     T         0     G
#> 552          6     T         2     G
#> 553          6     T         0     G
#> 554          1     T         1     G
#> 555         11     T         1     G
#> 556          8     T         0     G
#> 557         11     T         1     G
#> 558          8     T         0     G
#> 559          4     T         0     G
#> 560          3     T         1     G
#> 561         12     T         0     G
#> 562          8     T         0     G
#> 563         10     T         0     G
#> 564          8     T         0     G
#> 565         22     T         0     G
#> 566          9     T         3     G
#> 567         37     T         3     G
#> 568         34     T         0     G
#> 569         43     T         3     G
#> 570         28     T         2     G
#> 571          4     G         0     A
#> 572          3     G         1     A
#> 573         11     G         1     A
#> 574          8     G         0     A
#> 575         11     G         1     A
#> 576          5     G         1     A
#> 577          3     G         1     A
#> 578          2     G         0     A
#> 579         12     G         0     A
#> 580          8     G         0     A
#> 581         12     G         0     A
#> 582          8     G         0     A
#> 583          6     G         0     A
#> 584          2     G         0     A
#> 585         14     G         0     A
#> 586          8     G         0     A
#> 587         10     G         0     A
#> 588          7     G         1     A
#> 589          4     G         0     A
#> 590          4     G         0     A
#> 591         10     G         0     A
#> 592          6     G         0     A
#> 593         10     G         2     A
#> 594          6     G         0     A
#> 595         17     G         1     A
#> 596         11     G         1     A
#> 597         47     G         1     A
#> 598         30     G         0     A
#> 599         43     G         3     A
#> 600         26     G         2     A
#> 601          2     C         2     A
#> 602          3     C         1     A
#> 603          6     C         4     A
#> 604         10     C         0     A
#> 605          6     C         4     A
#> 606          2     C         0     A
#> 607          3     C         1     A
#> 608          1     C         1     A
#> 609          4     C         4     A
#> 610          8     C         0     A
#> 611          5     C         5     A
#> 612          6     C         2     A
#> 613          6     C         0     A
#> 614          0     C         0     A
#> 615          6     C         6     A
#> 616          8     C         0     A
#> 617          6     C         2     A
#> 618          3     C         3     A
#> 619          4     C         0     A
#> 620          1     C         1     A
#> 621          3     C         5     A
#> 622          6     C         0     A
#> 623         10     C         2     A
#> 624          5     C         3     A
#> 625         15     C         3     A
#> 626          5     C         3     A
#> 627         19     C        19     A
#> 628         32     C         0     A
#> 629         27     C        13     A
#> 630         16     C         8     A
#> 631          3     C         3     T
#> 632          3     C         1     T
#> 633          7     C         3     T
#> 634         10     C         0     T
#> 635          7     C         1     T
#> 636          2     C         2     T
#> 637          2     C         0     T
#> 638          0     C         0     T
#> 639          8     C         6     T
#> 640          8     C         0     T
#> 641          9     C         1     T
#> 642          6     C         2     T
#> 643          3     C         3     T
#> 644          0     C         2     T
#> 645         11     C         1     T
#> 646          8     C         0     T
#> 647          8     C         2     T
#> 648          4     C         4     T
#> 649          2     C         2     T
#> 650          3     C         1     T
#> 651          3     C         3     T
#> 652          8     C         0     T
#> 653         10     C         2     T
#> 654          6     C         2     T
#> 655         10     C         8     T
#> 656          6     C         4     T
#> 657         29     C        13     T
#> 658         34     C         0     T
#> 659         34     C         6     T
#> 660         18     C        10     T
#> 661          6     G         0     A
#> 662          4     G         0     A
#> 663         12     G         0     A
#> 664         10     G         0     A
#> 665         10     G         0     A
#> 666          6     G         0     A
#> 667          6     G         0     A
#> 668          2     G         0     A
#> 669         11     G         1     A
#> 670          8     G         0     A
#> 671         10     G         0     A
#> 672          7     G         1     A
#> 673          4     G         0     A
#> 674          2     G         0     A
#> 675         14     G         0     A
#> 676          6     G         0     A
#> 677          5     G         3     A
#> 678          8     G         0     A
#> 679          3     G         1     A
#> 680          4     G         0     A
#> 681         11     G         1     A
#> 682          8     G         0     A
#> 683         11     G         1     A
#> 684          8     G         0     A
#> 685         19     G         1     A
#> 686         12     G         0     A
#> 687         48     G         2     A
#> 688         32     G         0     A
#> 689         36     G         4     A
#> 690         29     G         1     A
#> 691          6     A         0     G
#> 692          4     A         0     G
#> 693         10     A         2     G
#> 694          8     A         0     G
#> 695         10     A         0     G
#> 696          4     A         2     G
#> 697          6     A         0     G
#> 698          2     A         0     G
#> 699         14     A         0     G
#> 700          8     A         0     G
#> 701          8     A         0     G
#> 702          7     A         1     G
#> 703          5     A         1     G
#> 704          2     A         0     G
#> 705         10     A         2     G
#> 706          8     A         0     G
#> 707         10     A         0     G
#> 708          8     A         0     G
#> 709          2     A         0     G
#> 710          4     A         0     G
#> 711         13     A         1     G
#> 712          8     A         0     G
#> 713          8     A         2     G
#> 714          6     A         0     G
#> 715         19     A         1     G
#> 716         12     A         0     G
#> 717         47     A         5     G
#> 718         32     A         0     G
#> 719         36     A         2     G
#> 720         25     A         3     G
#> 721          3     G         1     A
#> 722          4     G         0     A
#> 723         10     G         2     A
#> 724         10     G         0     A
#> 725          5     G         5     A
#> 726          3     G         3     A
#> 727          3     G         3     A
#> 728          1     G         1     A
#> 729          9     G         3     A
#> 730          8     G         0     A
#> 731          7     G         3     A
#> 732          4     G         4     A
#> 733          4     G         2     A
#> 734          0     G         0     A
#> 735          7     G         1     A
#> 736          6     G         0     A
#> 737          2     G         6     A
#> 738          5     G         3     A
#> 739          4     G         0     A
#> 740          4     G         0     A
#> 741          8     G         6     A
#> 742          8     G         0     A
#> 743          5     G         5     A
#> 744          5     G         3     A
#> 745         14     G         6     A
#> 746          9     G         1     A
#> 747         34     G        12     A
#> 748         32     G         0     A
#> 749         19     G        19     A
#> 750         17     G        13     A
#> 751          2     T         2     C
#> 752          4     T         0     C
#> 753         10     T         0     C
#> 754          4     T         6     C
#> 755          9     T         3     C
#> 756          3     T         1     C
#> 757          5     T         1     C
#> 758          2     T         0     C
#> 759          6     T         4     C
#> 760          5     T         3     C
#> 761          8     T         4     C
#> 762          4     T         2     C
#> 763          2     T         2     C
#> 764          2     T         0     C
#> 765         12     T         2     C
#> 766          4     T         4     C
#> 767          7     T         3     C
#> 768          6     T         2     C
#> 769          3     T         1     C
#> 770          4     T         0     C
#> 771          8     T         0     C
#> 772          1     T         3     C
#> 773          7     T         3     C
#> 774          4     T         2     C
#> 775         12     T         6     C
#> 776         12     T         0     C
#> 777         36     T         6     C
#> 778         14     T        16     C
#> 779         31     T        13     C
#> 780         17     T         7     C
#> 781          4     A         0     C
#> 782          3     A         1     C
#> 783          9     A         3     C
#> 784         10     A         0     C
#> 785          9     A         1     C
#> 786          6     A         0     C
#> 787          5     A         1     C
#> 788          2     A         0     C
#> 789         11     A         3     C
#> 790          8     A         0     C
#> 791         11     A         1     C
#> 792          6     A         0     C
#> 793          6     A         0     C
#> 794          2     A         0     C
#> 795          8     A         2     C
#> 796          8     A         0     C
#> 797          4     A         2     C
#> 798          7     A         1     C
#> 799          4     A         0     C
#> 800          4     A         0     C
#> 801         10     A         2     C
#> 802          8     A         0     C
#> 803         10     A         2     C
#> 804          7     A         1     C
#> 805         19     A         1     C
#> 806         11     A         1     C
#> 807         38     A        10     C
#> 808         34     A         0     C
#> 809         34     A         6     C
#> 810         26     A         2     C
#> 811          5     T         1     C
#> 812          4     T         0     C
#> 813          4     T         2     C
#> 814          7     T         1     C
#> 815          3     T         5     C
#> 816          4     T         2     C
#> 817          3     T         1     C
#> 818          1     T         1     C
#> 819          1     T         7     C
#> 820          6     T         0     C
#> 821          7     T         5     C
#> 822          5     T         1     C
#> 823          4     T         0     C
#> 824          2     T         0     C
#> 825          4     T         8     C
#> 826          7     T         1     C
#> 827          7     T         5     C
#> 828          6     T         2     C
#> 829          2     T         0     C
#> 830          3     T         1     C
#> 831          4     T         8     C
#> 832          7     T         1     C
#> 833          2     T         6     C
#> 834          5     T         3     C
#> 835         14     T         2     C
#> 836         10     T         2     C
#> 837         13     T        25     C
#> 838         27     T         3     C
#> 839         19     T        21     C
#> 840         20     T         8     C
#> 841          6     A         0     C
#> 842          2     A         0     C
#> 843          8     A         0     C
#> 844          7     A         1     C
#> 845          8     A         0     C
#> 846          5     A         1     C
#> 847          5     A         1     C
#> 848          0     A         0     C
#> 849         10     A         2     C
#> 850          5     A         1     C
#> 851          6     A         4     C
#> 852          6     A         0     C
#> 853          6     A         0     C
#> 854          2     A         0     C
#> 855          7     A         5     C
#> 856          7     A         1     C
#> 857          4     A         2     C
#> 858          5     A         1     C
#> 859          0     A         0     C
#> 860          3     A         1     C
#> 861          6     A         4     C
#> 862          6     A         2     C
#> 863          5     A         1     C
#> 864          7     A         1     C
#> 865         17     A         1     C
#> 866          7     A         1     C
#> 867         31     A        11     C
#> 868         25     A         5     C
#> 869         23     A         7     C
#> 870         23     A         3     C
#> 871          6     G         0     C
#> 872          2     G         0     C
#> 873          7     G         1     C
#> 874         10     G         0     C
#> 875          7     G         1     C
#> 876          3     G         1     C
#> 877          6     G         0     C
#> 878          1     G         1     C
#> 879         12     G         0     C
#> 880          8     G         0     C
#> 881         10     G         2     C
#> 882          5     G         1     C
#> 883          4     G         2     C
#> 884          2     G         0     C
#> 885         10     G         2     C
#> 886          8     G         0     C
#> 887         12     G         0     C
#> 888          8     G         0     C
#> 889          4     G         0     C
#> 890          3     G         1     C
#> 891         10     G         0     C
#> 892          8     G         0     C
#> 893         10     G         2     C
#> 894          4     G         2     C
#> 895         20     G         2     C
#> 896          8     G         2     C
#> 897         39     G         3     C
#> 898         34     G         0     C
#> 899         39     G         5     C
#> 900         20     G         4     C
#> 901          4     C         2     G
#> 902          3     C         1     G
#> 903          8     C         2     G
#> 904         10     C         0     G
#> 905          8     C         4     G
#> 906          4     C         2     G
#> 907          4     C         2     G
#> 908          2     C         0     G
#> 909         10     C         4     G
#> 910          8     C         0     G
#> 911          1     C         3     G
#> 912          3     C         5     G
#> 913          4     C         2     G
#> 914          2     C         0     G
#> 915          8     C         4     G
#> 916          8     C         0     G
#> 917          8     C         2     G
#> 918          7     C         1     G
#> 919          4     C         0     G
#> 920          4     C         0     G
#> 921          8     C         6     G
#> 922          8     C         0     G
#> 923          7     C         5     G
#> 924          3     C         5     G
#> 925         16     C         6     G
#> 926         11     C         1     G
#> 927         34     C        16     G
#> 928         34     C         0     G
#> 929         24     C        14     G
#> 930         17     C        13     G
#> 931          5     A         1     C
#> 932          3     A         1     C
#> 933          9     A         1     C
#> 934          7     A         3     C
#> 935          7     A         1     C
#> 936          2     A         2     C
#> 937          5     A         1     C
#> 938          0     A         0     C
#> 939          9     A         1     C
#> 940          3     A         3     C
#> 941         10     A         2     C
#> 942          5     A         3     C
#> 943          6     A         0     C
#> 944          2     A         0     C
#> 945         11     A         1     C
#> 946          5     A         3     C
#> 947          8     A         2     C
#> 948          5     A         1     C
#> 949          1     A         1     C
#> 950          3     A         1     C
#> 951         11     A         3     C
#> 952          3     A         3     C
#> 953          8     A         2     C
#> 954          7     A         1     C
#> 955         17     A         3     C
#> 956          8     A         2     C
#> 957         40     A         6     C
#> 958         18     A        12     C
#> 959         33     A         7     C
#> 960         19     A         7     C
#> 961          5     C         1     T
#> 962          2     C         0     T
#> 963          5     C         5     T
#> 964          7     C         1     T
#> 965         10     C         2     T
#> 966          2     C         2     T
#> 967          3     C         3     T
#> 968          1     C         1     T
#> 969          6     C         4     T
#> 970          6     C         2     T
#> 971          6     C         6     T
#> 972          6     C         2     T
#> 973          4     C         2     T
#> 974          0     C         0     T
#> 975          7     C         1     T
#> 976          7     C         1     T
#> 977          5     C         3     T
#> 978          3     C         3     T
#> 979          4     C         0     T
#> 980          2     C         2     T
#> 981          6     C         4     T
#> 982          3     C         3     T
#> 983          8     C         4     T
#> 984          6     C         2     T
#> 985         16     C         6     T
#> 986          5     C         3     T
#> 987         24     C        14     T
#> 988         23     C         7     T
#> 989         29     C        15     T
#> 990         17     C         9     T
#> 991          4     C         2     G
#> 992          3     C         1     G
#> 993          7     C         3     G
#> 994         10     C         0     G
#> 995          6     C         6     G
#> 996          5     C         1     G
#> 997          2     C         2     G
#> 998          0     C         2     G
#> 999          7     C         3     G
#> 1000         8     C         0     G
#> 1001        10     C         2     G
#> 1002         5     C         1     G
#> 1003         3     C         1     G
#> 1004         1     C         1     G
#> 1005         3     C         3     G
#> 1006         8     C         0     G
#> 1007         7     C         1     G
#> 1008         6     C         0     G
#> 1009         2     C         2     G
#> 1010         2     C         2     G
#> 1011         5     C         5     G
#> 1012         8     C         0     G
#> 1013         8     C         4     G
#> 1014         3     C         1     G
#> 1015        11     C         7     G
#> 1016         6     C         6     G
#> 1017        22     C        14     G
#> 1018        34     C         0     G
#> 1019        31     C        13     G
#> 1020        19     C         3     G
#> 1021         6     A         0     C
#> 1022         2     A         0     C
#> 1023         8     A         2     C
#> 1024        10     A         0     C
#> 1025         8     A         0     C
#> 1026         4     A         0     C
#> 1027         4     A         0     C
#> 1028         2     A         0     C
#> 1029        11     A         1     C
#> 1030         8     A         0     C
#> 1031         7     A         1     C
#> 1032         7     A         1     C
#> 1033         5     A         1     C
#> 1034         2     A         0     C
#> 1035        12     A         0     C
#> 1036         8     A         0     C
#> 1037        12     A         0     C
#> 1038         8     A         0     C
#> 1039         3     A         1     C
#> 1040         4     A         0     C
#> 1041        10     A         0     C
#> 1042         6     A         0     C
#> 1043        10     A         2     C
#> 1044         7     A         1     C
#> 1045        18     A         2     C
#> 1046        10     A         0     C
#> 1047        41     A         3     C
#> 1048        32     A         0     C
#> 1049        37     A         3     C
#> 1050        26     A         2     C
#> 1051         6     A         0     G
#> 1052         4     A         0     G
#> 1053        10     A         2     G
#> 1054        10     A         0     G
#> 1055         8     A         4     G
#> 1056         2     A         2     G
#> 1057         4     A         2     G
#> 1058         2     A         0     G
#> 1059         8     A         2     G
#> 1060         8     A         0     G
#> 1061        10     A         2     G
#> 1062         8     A         0     G
#> 1063         4     A         0     G
#> 1064         1     A         1     G
#> 1065        10     A         2     G
#> 1066         8     A         0     G
#> 1067        10     A         0     G
#> 1068         6     A         0     G
#> 1069         3     A         1     G
#> 1070         3     A         1     G
#> 1071        11     A         3     G
#> 1072         6     A         0     G
#> 1073         8     A         2     G
#> 1074         8     A         0     G
#> 1075        17     A         3     G
#> 1076        10     A         2     G
#> 1077        39     A         9     G
#> 1078        32     A         0     G
#> 1079        36     A         8     G
#> 1080        24     A         2     G
#> 1081         3     G         3     C
#> 1082         3     G         1     C
#> 1083         7     G         5     C
#> 1084         8     G         0     C
#> 1085         7     G         3     C
#> 1086         3     G         1     C
#> 1087         6     G         0     C
#> 1088         1     G         1     C
#> 1089        10     G         4     C
#> 1090         8     G         0     C
#> 1091         6     G         6     C
#> 1092         5     G         3     C
#> 1093         3     G         1     C
#> 1094         0     G         0     C
#> 1095         7     G         3     C
#> 1096         8     G         0     C
#> 1097        10     G         2     C
#> 1098         4     G         2     C
#> 1099         4     G         0     C
#> 1100         4     G         0     C
#> 1101        11     G         1     C
#> 1102         7     G         1     C
#> 1103         5     G         5     C
#> 1104         7     G         1     C
#> 1105        16     G         4     C
#> 1106         8     G         2     C
#> 1107        35     G        13     C
#> 1108        31     G         1     C
#> 1109        28     G        16     C
#> 1110        19     G         7     C
#> 1111         3     G         1     T
#> 1112         4     G         0     T
#> 1113         7     G         1     T
#> 1114         5     G         5     T
#> 1115         6     G         4     T
#> 1116         4     G         2     T
#> 1117         2     G         4     T
#> 1118         2     G         0     T
#> 1119         9     G         3     T
#> 1120         5     G         3     T
#> 1121         7     G         3     T
#> 1122         5     G         3     T
#> 1123         4     G         2     T
#> 1124         2     G         0     T
#> 1125        12     G         2     T
#> 1126         4     G         2     T
#> 1127         5     G         3     T
#> 1128         8     G         0     T
#> 1129         2     G         2     T
#> 1130         4     G         0     T
#> 1131         9     G         3     T
#> 1132         5     G         3     T
#> 1133         4     G         6     T
#> 1134         4     G         2     T
#> 1135        11     G         9     T
#> 1136        12     G         0     T
#> 1137        37     G         9     T
#> 1138        19     G        13     T
#> 1139        22     G        16     T
#> 1140        21     G         7     T
#> 1141         5     T         1     A
#> 1142         1     T         1     A
#> 1143        10     T         0     A
#> 1144         6     T         0     A
#> 1145        12     T         0     A
#> 1146         5     T         1     A
#> 1147         3     T         1     A
#> 1148         0     T         0     A
#> 1149         4     T         0     A
#> 1150         5     T         1     A
#> 1151         5     T         1     A
#> 1152         4     T         0     A
#> 1153         5     T         1     A
#> 1154         1     T         1     A
#> 1155         5     T         1     A
#> 1156         7     T         1     A
#> 1157         7     T         1     A
#> 1158         5     T         1     A
#> 1159         2     T         0     A
#> 1160         2     T         0     A
#> 1161         8     T         2     A
#> 1162         1     T         1     A
#> 1163         7     T         1     A
#> 1164         3     T         1     A
#> 1165        15     T         3     A
#> 1166         4     T         2     A
#> 1167        27     T         3     A
#> 1168        19     T         3     A
#> 1169        31     T         3     A
#> 1170        17     T         3     A
#> 1171         2     G         4     C
#> 1172         4     G         0     C
#> 1173         7     G         1     C
#> 1174         8     G         0     C
#> 1175         6     G         2     C
#> 1176         4     G         2     C
#> 1177         4     G         2     C
#> 1178         0     G         0     C
#> 1179         8     G         2     C
#> 1180         5     G         3     C
#> 1181        10     G         0     C
#> 1182         3     G         1     C
#> 1183         3     G         3     C
#> 1184         1     G         1     C
#> 1185         6     G         2     C
#> 1186         7     G         1     C
#> 1187         7     G         3     C
#> 1188         1     G         5     C
#> 1189         2     G         0     C
#> 1190         3     G         1     C
#> 1191        12     G         2     C
#> 1192         5     G         1     C
#> 1193         9     G         1     C
#> 1194         5     G         1     C
#> 1195        11     G         9     C
#> 1196         8     G         2     C
#> 1197        33     G         7     C
#> 1198        25     G         5     C
#> 1199        32     G         6     C
#> 1200        13     G         9     C
#> 1201         3     G         3     A
#> 1202         3     G         1     A
#> 1203         5     G         3     A
#> 1204         0     G         8     A
#> 1205         6     G         2     A
#> 1206         0     G         2     A
#> 1207         1     G         5     A
#> 1208         2     G         0     A
#> 1209         7     G         3     A
#> 1210         0     G         6     A
#> 1211         8     G         4     A
#> 1212         4     G         4     A
#> 1213         3     G         3     A
#> 1214         1     G         1     A
#> 1215         9     G         1     A
#> 1216         0     G         4     A
#> 1217         3     G         3     A
#> 1218         2     G         2     A
#> 1219         2     G         2     A
#> 1220         1     G         1     A
#> 1221         7     G         3     A
#> 1222         0     G         8     A
#> 1223         6     G         2     A
#> 1224         2     G         4     A
#> 1225         9     G        13     A
#> 1226         7     G         3     A
#> 1227        28     G        10     A
#> 1228         0     G        26     A
#> 1229        23     G        11     A
#> 1230         8     G        12     A
#> 1231         4     G         2     A
#> 1232         2     G         2     A
#> 1233         5     G         7     A
#> 1234         9     G         1     A
#> 1235         4     G         6     A
#> 1236         3     G         3     A
#> 1237         4     G         0     A
#> 1238         1     G         1     A
#> 1239         6     G         6     A
#> 1240         8     G         0     A
#> 1241         6     G         6     A
#> 1242         6     G         2     A
#> 1243         3     G         1     A
#> 1244         2     G         0     A
#> 1245         7     G         3     A
#> 1246         7     G         1     A
#> 1247         4     G         8     A
#> 1248         4     G         2     A
#> 1249         1     G         1     A
#> 1250         1     G         1     A
#> 1251         7     G         5     A
#> 1252         6     G         0     A
#> 1253         5     G         5     A
#> 1254         5     G         3     A
#> 1255        12     G         4     A
#> 1256         6     G         4     A
#> 1257        25     G        21     A
#> 1258        30     G         2     A
#> 1259        19     G        25     A
#> 1260        18     G        10     A
#> 1261         5     A         1     G
#> 1262         2     A         0     G
#> 1263         7     A         1     G
#> 1264        10     A         0     G
#> 1265        11     A         1     G
#> 1266         4     A         2     G
#> 1267         5     A         1     G
#> 1268         2     A         0     G
#> 1269         7     A         1     G
#> 1270         8     A         0     G
#> 1271        11     A         1     G
#> 1272         4     A         0     G
#> 1273         6     A         0     G
#> 1274         2     A         0     G
#> 1275         9     A         3     G
#> 1276         8     A         0     G
#> 1277        12     A         0     G
#> 1278         6     A         2     G
#> 1279         4     A         0     G
#> 1280         2     A         2     G
#> 1281        11     A         3     G
#> 1282         8     A         0     G
#> 1283         9     A         3     G
#> 1284         6     A         2     G
#> 1285        20     A         2     G
#> 1286         8     A         2     G
#> 1287        34     A         8     G
#> 1288        34     A         0     G
#> 1289        43     A         5     G
#> 1290        20     A         6     G
#> 1291         3     T         1     A
#> 1292         4     T         0     A
#> 1293         9     T         1     A
#> 1294         7     T         3     A
#> 1295        10     T         0     A
#> 1296         2     T         0     A
#> 1297         4     T         2     A
#> 1298         2     T         0     A
#> 1299        14     T         0     A
#> 1300         2     T         4     A
#> 1301         7     T         3     A
#> 1302         7     T         1     A
#> 1303         6     T         0     A
#> 1304         2     T         0     A
#> 1305        13     T         1     A
#> 1306         2     T         4     A
#> 1307         9     T         1     A
#> 1308         3     T         3     A
#> 1309         2     T         0     A
#> 1310         4     T         0     A
#> 1311        13     T         1     A
#> 1312         7     T         1     A
#> 1313        11     T         1     A
#> 1314         7     T         1     A
#> 1315        15     T         3     A
#> 1316        12     T         0     A
#> 1317        49     T         3     A
#> 1318        18     T        12     A
#> 1319        37     T         5     A
#> 1320        19     T         5     A
#> 1321         6     G         0     T
#> 1322         4     G         0     T
#> 1323        10     G         0     T
#> 1324        10     G         0     T
#> 1325        10     G         0     T
#> 1326         6     G         0     T
#> 1327         6     G         0     T
#> 1328         2     G         0     T
#> 1329         6     G         0     T
#> 1330         8     G         0     T
#> 1331        10     G         0     T
#> 1332         7     G         1     T
#> 1333         6     G         0     T
#> 1334         2     G         0     T
#> 1335        12     G         0     T
#> 1336         8     G         0     T
#> 1337        10     G         0     T
#> 1338         8     G         0     T
#> 1339         4     G         0     T
#> 1340         4     G         0     T
#> 1341        10     G         0     T
#> 1342         8     G         0     T
#> 1343        10     G         0     T
#> 1344         7     G         1     T
#> 1345        22     G         0     T
#> 1346        12     G         0     T
#> 1347        38     G         0     T
#> 1348        34     G         0     T
#> 1349        40     G         0     T
#> 1350        28     G         2     T
#> 1351         2     C         0     A
#> 1352         3     C         1     A
#> 1353         6     C         2     A
#> 1354         8     C         0     A
#> 1355         7     C         3     A
#> 1356         2     C         2     A
#> 1357         6     C         0     A
#> 1358         2     C         0     A
#> 1359         8     C         4     A
#> 1360         8     C         0     A
#> 1361         7     C         1     A
#> 1362         5     C         3     A
#> 1363         1     C         3     A
#> 1364         0     C         0     A
#> 1365         3     C         3     A
#> 1366         8     C         0     A
#> 1367         3     C         1     A
#> 1368         5     C         1     A
#> 1369         0     C         4     A
#> 1370         3     C         1     A
#> 1371         4     C         6     A
#> 1372         8     C         0     A
#> 1373        10     C         2     A
#> 1374         5     C         3     A
#> 1375         9     C         7     A
#> 1376         8     C         2     A
#> 1377        21     C        15     A
#> 1378        32     C         0     A
#> 1379        27     C         7     A
#> 1380        17     C         9     A
#> 1381         5     C         1     G
#> 1382         4     C         0     G
#> 1383         3     C         9     G
#> 1384         8     C         0     G
#> 1385         8     C         4     G
#> 1386         3     C         3     G
#> 1387         5     C         1     G
#> 1388         2     C         0     G
#> 1389         8     C         4     G
#> 1390         5     C         1     G
#> 1391        11     C         1     G
#> 1392         8     C         0     G
#> 1393         5     C         1     G
#> 1394         2     C         0     G
#> 1395         7     C         5     G
#> 1396         7     C         1     G
#> 1397         9     C         3     G
#> 1398         8     C         0     G
#> 1399         2     C         2     G
#> 1400         4     C         0     G
#> 1401         6     C         4     G
#> 1402         7     C         1     G
#> 1403         7     C         3     G
#> 1404         6     C         2     G
#> 1405        17     C         5     G
#> 1406        12     C         0     G
#> 1407        24     C        22     G
#> 1408        27     C         3     G
#> 1409        35     C        11     G
#> 1410        25     C         5     G
#> 1411         3     T         1     C
#> 1412         1     T         3     C
#> 1413         5     T         5     C
#> 1414         3     T         5     C
#> 1415         3     T         5     C
#> 1416         3     T         3     C
#> 1417         6     T         0     C
#> 1418         2     T         0     C
#> 1419         4     T         6     C
#> 1420         5     T         3     C
#> 1421         7     T         5     C
#> 1422         3     T         3     C
#> 1423         4     T         2     C
#> 1424         1     T         1     C
#> 1425         5     T         7     C
#> 1426         6     T         2     C
#> 1427         5     T         1     C
#> 1428         2     T         4     C
#> 1429         0     T         2     C
#> 1430         2     T         2     C
#> 1431         3     T         5     C
#> 1432         3     T         1     C
#> 1433         5     T         5     C
#> 1434         8     T         0     C
#> 1435        13     T         5     C
#> 1436         6     T         6     C
#> 1437        17     T        23     C
#> 1438        17     T        11     C
#> 1439        20     T        16     C
#> 1440        16     T        10     C
#> 1441         4     G         0     C
#> 1442         4     G         0     C
#> 1443        12     G         0     C
#> 1444        10     G         0     C
#> 1445        12     G         0     C
#> 1446         5     G         1     C
#> 1447         6     G         0     C
#> 1448         1     G         1     C
#> 1449        13     G         1     C
#> 1450         8     G         0     C
#> 1451        12     G         0     C
#> 1452         7     G         1     C
#> 1453         6     G         0     C
#> 1454         2     G         0     C
#> 1455        14     G         0     C
#> 1456         8     G         0     C
#> 1457        12     G         0     C
#> 1458         8     G         0     C
#> 1459         4     G         0     C
#> 1460         3     G         1     C
#> 1461        14     G         0     C
#> 1462         8     G         0     C
#> 1463        12     G         0     C
#> 1464         4     G         0     C
#> 1465        20     G         0     C
#> 1466        10     G         2     C
#> 1467        53     G         1     C
#> 1468        34     G         0     C
#> 1469        48     G         0     C
#> 1470        24     G         2     C
#> 1471         3     C         1     T
#> 1472         2     C         0     T
#> 1473         9     C         1     T
#> 1474         8     C         0     T
#> 1475         8     C         0     T
#> 1476         6     C         0     T
#> 1477         5     C         1     T
#> 1478         2     C         0     T
#> 1479         9     C         3     T
#> 1480         7     C         1     T
#> 1481         8     C         2     T
#> 1482         8     C         0     T
#> 1483         5     C         1     T
#> 1484         2     C         0     T
#> 1485        12     C         0     T
#> 1486         8     C         0     T
#> 1487        10     C         0     T
#> 1488         8     C         0     T
#> 1489         4     C         0     T
#> 1490         3     C         1     T
#> 1491        11     C         1     T
#> 1492         8     C         0     T
#> 1493        10     C         0     T
#> 1494         4     C         2     T
#> 1495        17     C         3     T
#> 1496         9     C         1     T
#> 1497        41     C         5     T
#> 1498        31     C         1     T
#> 1499        36     C         2     T
#> 1500        26     C         2     T
#> 1501         5     T         1     C
#> 1502         2     T         0     C
#> 1503         7     T         3     C
#> 1504         3     T         1     C
#> 1505         6     T         2     C
#> 1506         3     T         3     C
#> 1507         1     T         3     C
#> 1508         2     T         0     C
#> 1509         9     T         3     C
#> 1510         5     T         1     C
#> 1511         5     T         1     C
#> 1512         5     T         3     C
#> 1513         4     T         2     C
#> 1514         2     T         0     C
#> 1515         5     T         5     C
#> 1516         5     T         1     C
#> 1517         3     T         1     C
#> 1518         6     T         2     C
#> 1519         1     T         1     C
#> 1520         4     T         0     C
#> 1521         5     T         3     C
#> 1522         5     T         1     C
#> 1523         5     T         5     C
#> 1524         3     T         3     C
#> 1525        11     T         7     C
#> 1526        10     T         0     C
#> 1527        26     T        14     C
#> 1528        18     T         4     C
#> 1529        19     T         9     C
#> 1530        17     T        11     C
#> 1531         4     A         2     G
#> 1532         2     A         2     G
#> 1533         7     A         3     G
#> 1534         4     A         6     G
#> 1535        11     A         1     G
#> 1536         3     A         3     G
#> 1537         3     A         3     G
#> 1538         2     A         0     G
#> 1539         9     A         3     G
#> 1540         4     A         2     G
#> 1541        11     A         1     G
#> 1542         5     A         3     G
#> 1543         4     A         2     G
#> 1544         2     A         0     G
#> 1545         7     A         5     G
#> 1546         1     A         1     G
#> 1547        10     A         2     G
#> 1548         4     A         2     G
#> 1549         3     A         1     G
#> 1550         2     A         2     G
#> 1551         9     A         5     G
#> 1552         4     A         2     G
#> 1553         9     A         3     G
#> 1554         6     A         2     G
#> 1555        14     A         8     G
#> 1556         8     A         4     G
#> 1557        32     A        16     G
#> 1558        13     A        11     G
#> 1559        41     A         7     G
#> 1560        18     A        10     G
#> 1561         5     A         1     C
#> 1562         1     A         1     C
#> 1563         7     A         3     C
#> 1564        10     A         0     C
#> 1565        10     A         0     C
#> 1566         2     A         2     C
#> 1567         5     A         1     C
#> 1568         2     A         0     C
#> 1569         3     A         5     C
#> 1570         6     A         0     C
#> 1571        11     A         1     C
#> 1572         5     A         1     C
#> 1573         4     A         2     C
#> 1574         0     A         0     C
#> 1575         8     A         2     C
#> 1576         8     A         0     C
#> 1577        10     A         2     C
#> 1578         2     A         4     C
#> 1579         2     A         0     C
#> 1580         2     A         2     C
#> 1581        11     A         3     C
#> 1582         8     A         0     C
#> 1583         7     A         5     C
#> 1584         5     A         3     C
#> 1585        16     A         4     C
#> 1586         5     A         3     C
#> 1587        29     A        13     C
#> 1588        32     A         0     C
#> 1589        38     A         8     C
#> 1590        14     A        10     C
#> 1591         6     G         0     T
#> 1592         4     G         0     T
#> 1593        11     G         1     T
#> 1594        10     G         0     T
#> 1595        11     G         1     T
#> 1596         6     G         0     T
#> 1597         6     G         0     T
#> 1598         2     G         0     T
#> 1599        12     G         0     T
#> 1600         8     G         0     T
#> 1601        12     G         0     T
#> 1602         8     G         0     T
#> 1603         6     G         0     T
#> 1604         1     G         1     T
#> 1605        14     G         0     T
#> 1606         8     G         0     T
#> 1607        10     G         0     T
#> 1608         8     G         0     T
#> 1609         4     G         0     T
#> 1610         4     G         0     T
#> 1611        14     G         0     T
#> 1612         6     G         0     T
#> 1613        12     G         0     T
#> 1614         8     G         0     T
#> 1615        22     G         0     T
#> 1616        11     G         1     T
#> 1617        51     G         1     T
#> 1618        32     G         0     T
#> 1619        45     G         1     T
#> 1620        30     G         0     T
#> 1621         4     G         2     A
#> 1622         3     G         1     A
#> 1623         5     G         5     A
#> 1624         3     G         5     A
#> 1625        11     G         1     A
#> 1626         1     G         5     A
#> 1627         1     G         3     A
#> 1628         2     G         0     A
#> 1629         9     G         3     A
#> 1630         3     G         3     A
#> 1631         7     G         3     A
#> 1632         5     G         3     A
#> 1633         2     G         2     A
#> 1634         1     G         1     A
#> 1635         7     G         3     A
#> 1636         6     G         2     A
#> 1637         8     G         0     A
#> 1638         4     G         2     A
#> 1639         1     G         1     A
#> 1640         3     G         1     A
#> 1641         7     G         3     A
#> 1642         4     G         4     A
#> 1643         7     G         3     A
#> 1644         4     G         2     A
#> 1645         8     G         8     A
#> 1646         9     G         3     A
#> 1647        28     G        14     A
#> 1648        16     G        14     A
#> 1649        33     G         7     A
#> 1650        14     G        12     A
#> 1651         4     A         2     G
#> 1652         4     A         0     G
#> 1653         7     A         1     G
#> 1654         5     A         3     G
#> 1655         8     A         2     G
#> 1656         3     A         3     G
#> 1657         3     A         1     G
#> 1658         1     A         1     G
#> 1659         8     A         2     G
#> 1660         5     A         3     G
#> 1661        10     A         2     G
#> 1662         7     A         1     G
#> 1663         2     A         2     G
#> 1664         2     A         0     G
#> 1665         9     A         3     G
#> 1666         4     A         4     G
#> 1667         6     A         2     G
#> 1668         6     A         2     G
#> 1669         3     A         1     G
#> 1670         4     A         0     G
#> 1671         8     A         2     G
#> 1672         7     A         1     G
#> 1673         7     A         1     G
#> 1674         7     A         1     G
#> 1675        12     A         6     G
#> 1676        11     A         1     G
#> 1677        32     A         8     G
#> 1678        21     A        11     G
#> 1679        31     A         7     G
#> 1680        23     A         7     G
#> 1681         4     G         0     A
#> 1682         4     G         0     A
#> 1683         6     G         0     A
#> 1684         8     G         0     A
#> 1685         8     G         0     A
#> 1686         4     G         0     A
#> 1687         4     G         0     A
#> 1688         2     G         0     A
#> 1689        10     G         2     A
#> 1690         6     G         0     A
#> 1691         8     G         0     A
#> 1692         6     G         0     A
#> 1693         6     G         0     A
#> 1694         2     G         0     A
#> 1695         9     G         1     A
#> 1696         4     G         0     A
#> 1697        11     G         1     A
#> 1698         2     G         0     A
#> 1699         1     G         1     A
#> 1700         2     G         0     A
#> 1701        11     G         1     A
#> 1702         6     G         0     A
#> 1703         8     G         0     A
#> 1704         8     G         0     A
#> 1705        15     G         1     A
#> 1706        10     G         0     A
#> 1707        36     G         4     A
#> 1708        24     G         0     A
#> 1709        35     G         1     A
#> 1710        20     G         0     A
#> 1711         2     C         4     T
#> 1712         4     C         0     T
#> 1713         7     C         3     T
#> 1714         5     C         3     T
#> 1715         8     C         0     T
#> 1716         4     C         0     T
#> 1717         4     C         2     T
#> 1718         2     C         0     T
#> 1719         9     C         3     T
#> 1720         1     C         3     T
#> 1721         8     C         2     T
#> 1722         4     C         0     T
#> 1723         3     C         1     T
#> 1724         2     C         0     T
#> 1725        10     C         4     T
#> 1726         5     C         3     T
#> 1727         4     C         2     T
#> 1728         4     C         4     T
#> 1729         1     C         1     T
#> 1730         4     C         0     T
#> 1731        11     C         3     T
#> 1732         4     C         2     T
#> 1733        11     C         1     T
#> 1734         4     C         2     T
#> 1735        10     C         8     T
#> 1736        12     C         0     T
#> 1737        37     C        13     T
#> 1738        15     C        11     T
#> 1739        31     C         5     T
#> 1740        16     C         6     T
#> 1741         4     G         2     A
#> 1742         3     G         1     A
#> 1743         7     G         3     A
#> 1744         8     G         0     A
#> 1745         4     G         6     A
#> 1746         3     G         3     A
#> 1747         1     G         5     A
#> 1748         1     G         1     A
#> 1749         7     G         7     A
#> 1750         7     G         1     A
#> 1751         4     G         4     A
#> 1752         5     G         3     A
#> 1753         5     G         1     A
#> 1754         1     G         1     A
#> 1755         6     G         6     A
#> 1756         8     G         0     A
#> 1757         6     G         6     A
#> 1758         4     G         2     A
#> 1759         2     G         2     A
#> 1760         3     G         1     A
#> 1761         9     G         3     A
#> 1762         7     G         1     A
#> 1763         7     G         3     A
#> 1764         4     G         2     A
#> 1765        12     G        10     A
#> 1766         8     G         4     A
#> 1767        29     G        19     A
#> 1768        30     G         2     A
#> 1769        21     G        19     A
#> 1770        16     G        10     A
#> 1771         3     A         1     C
#> 1772         2     A         2     C
#> 1773         8     A         2     C
#> 1774        10     A         0     C
#> 1775         8     A         4     C
#> 1776         4     A         2     C
#> 1777         3     A         3     C
#> 1778         0     A         2     C
#> 1779         7     A         5     C
#> 1780         8     A         0     C
#> 1781         8     A         2     C
#> 1782         5     A         3     C
#> 1783         2     A         2     C
#> 1784         0     A         2     C
#> 1785         8     A         4     C
#> 1786         8     A         0     C
#> 1787         4     A         4     C
#> 1788         7     A         1     C
#> 1789         0     A         0     C
#> 1790         1     A         3     C
#> 1791         7     A         5     C
#> 1792         8     A         0     C
#> 1793         8     A         4     C
#> 1794         3     A         3     C
#> 1795         8     A         6     C
#> 1796         3     A         9     C
#> 1797        30     A        16     C
#> 1798        34     A         0     C
#> 1799        28     A        14     C
#> 1800        19     A         9     C
#> 1801         6     T         0     C
#> 1802         0     T         0     C
#> 1803        10     T         0     C
#> 1804         4     T         0     C
#> 1805        10     T         0     C
#> 1806         6     T         0     C
#> 1807         4     T         0     C
#> 1808         0     T         0     C
#> 1809        14     T         0     C
#> 1810         8     T         0     C
#> 1811        12     T         0     C
#> 1812         7     T         1     C
#> 1813         6     T         0     C
#> 1814         1     T         1     C
#> 1815         3     T         1     C
#> 1816         8     T         0     C
#> 1817         8     T         2     C
#> 1818         3     T         1     C
#> 1819         0     T         0     C
#> 1820         4     T         0     C
#> 1821         6     T         0     C
#> 1822         6     T         0     C
#> 1823         8     T         0     C
#> 1824         6     T         0     C
#> 1825        16     T         0     C
#> 1826         5     T         1     C
#> 1827        33     T         1     C
#> 1828        26     T         0     C
#> 1829        38     T         2     C
#> 1830        22     T         2     C
#> 1831         4     G         2     C
#> 1832         4     G         0     C
#> 1833         9     G         3     C
#> 1834        10     G         0     C
#> 1835        11     G         1     C
#> 1836         5     G         1     C
#> 1837         6     G         0     C
#> 1838         2     G         0     C
#> 1839        11     G         3     C
#> 1840         8     G         0     C
#> 1841        11     G         1     C
#> 1842         7     G         1     C
#> 1843         6     G         0     C
#> 1844         2     G         0     C
#> 1845        11     G         1     C
#> 1846         8     G         0     C
#> 1847        12     G         0     C
#> 1848         6     G         0     C
#> 1849         2     G         0     C
#> 1850         4     G         0     C
#> 1851        10     G         4     C
#> 1852         6     G         0     C
#> 1853        11     G         1     C
#> 1854         7     G         1     C
#> 1855        18     G         2     C
#> 1856        12     G         0     C
#> 1857        41     G        11     C
#> 1858        32     G         0     C
#> 1859        45     G         3     C
#> 1860        25     G         3     C
#> 1861         3     A         3     C
#> 1862         3     A         1     C
#> 1863         6     A         2     C
#> 1864         1     A         3     C
#> 1865         7     A         3     C
#> 1866         5     A         1     C
#> 1867         3     A         1     C
#> 1868         1     A         1     C
#> 1869         5     A         3     C
#> 1870         1     A         3     C
#> 1871         6     A         2     C
#> 1872         5     A         3     C
#> 1873         3     A         1     C
#> 1874         2     A         0     C
#> 1875        10     A         2     C
#> 1876         1     A         5     C
#> 1877         4     A         4     C
#> 1878         3     A         1     C
#> 1879         1     A         1     C
#> 1880         0     A         2     C
#> 1881         7     A         1     C
#> 1882         1     A         5     C
#> 1883         5     A         1     C
#> 1884         6     A         2     C
#> 1885        10     A         6     C
#> 1886         6     A         4     C
#> 1887        28     A         8     C
#> 1888         4     A        16     C
#> 1889        22     A        10     C
#> 1890        19     A         7     C
#> 1891         5     T         1     A
#> 1892         2     T         2     A
#> 1893         6     T         2     A
#> 1894         3     T         7     A
#> 1895         6     T         6     A
#> 1896         4     T         2     A
#> 1897         4     T         2     A
#> 1898         2     T         0     A
#> 1899         8     T         4     A
#> 1900         5     T         3     A
#> 1901         7     T         3     A
#> 1902         3     T         5     A
#> 1903         3     T         3     A
#> 1904         2     T         0     A
#> 1905        10     T         2     A
#> 1906         0     T         4     A
#> 1907         8     T         4     A
#> 1908         6     T         2     A
#> 1909         3     T         1     A
#> 1910         4     T         0     A
#> 1911         5     T         5     A
#> 1912         4     T         2     A
#> 1913         7     T         3     A
#> 1914         6     T         2     A
#> 1915        15     T         7     A
#> 1916        10     T         2     A
#> 1917        29     T        13     A
#> 1918        12     T        16     A
#> 1919        28     T        16     A
#> 1920        19     T        11     A
#> 1921         5     C         1     T
#> 1922         4     C         0     T
#> 1923         4     C         4     T
#> 1924         7     C         3     T
#> 1925         7     C         5     T
#> 1926         3     C         1     T
#> 1927         5     C         1     T
#> 1928         1     C         1     T
#> 1929         4     C         6     T
#> 1930         5     C         3     T
#> 1931         7     C         1     T
#> 1932         4     C         4     T
#> 1933         4     C         2     T
#> 1934         1     C         1     T
#> 1935         1     C         9     T
#> 1936         6     C         2     T
#> 1937         8     C         4     T
#> 1938         7     C         1     T
#> 1939         4     C         0     T
#> 1940         1     C         3     T
#> 1941         5     C         5     T
#> 1942         7     C         1     T
#> 1943         7     C         5     T
#> 1944         3     C         5     T
#> 1945        18     C         4     T
#> 1946         7     C         5     T
#> 1947        14     C        24     T
#> 1948        25     C         9     T
#> 1949        29     C        15     T
#> 1950        17     C        11     T
#> 1951         4     A         2     C
#> 1952         3     A         1     C
#> 1953         7     A         3     C
#> 1954         3     A         5     C
#> 1955         4     A         4     C
#> 1956         5     A         1     C
#> 1957         3     A         3     C
#> 1958         0     A         2     C
#> 1959         7     A         3     C
#> 1960         2     A         6     C
#> 1961         5     A         1     C
#> 1962         4     A         4     C
#> 1963         1     A         3     C
#> 1964         2     A         0     C
#> 1965         6     A         4     C
#> 1966         3     A         5     C
#> 1967         7     A         3     C
#> 1968         1     A         3     C
#> 1969         1     A         1     C
#> 1970         2     A         2     C
#> 1971         6     A         4     C
#> 1972         3     A         5     C
#> 1973         5     A         5     C
#> 1974         5     A         3     C
#> 1975         9     A         9     C
#> 1976         7     A         5     C
#> 1977        26     A        14     C
#> 1978        11     A        21     C
#> 1979        21     A        13     C
#> 1980        15     A        11     C
#> 1981         4     G         0     A
#> 1982         4     G         0     A
#> 1983        10     G         0     A
#> 1984         8     G         2     A
#> 1985        12     G         0     A
#> 1986         6     G         0     A
#> 1987         6     G         0     A
#> 1988         2     G         0     A
#> 1989        14     G         0     A
#> 1990         4     G         4     A
#> 1991        10     G         0     A
#> 1992         6     G         0     A
#> 1993         6     G         0     A
#> 1994         2     G         0     A
#> 1995        12     G         0     A
#> 1996         8     G         0     A
#> 1997        12     G         0     A
#> 1998         6     G         0     A
#> 1999         2     G         0     A
#> 2000         4     G         0     A
#> 2001        14     G         0     A
#> 2002         2     G         2     A
#> 2003        10     G         0     A
#> 2004         7     G         1     A
#> 2005        18     G         0     A
#> 2006        12     G         0     A
#> 2007        50     G         0     A
#> 2008        22     G         8     A
#> 2009        44     G         0     A
#> 2010        25     G         1     A
#> 2011         6     A         0     T
#> 2012         3     A         1     T
#> 2013         7     A         5     T
#> 2014         1     A         9     T
#> 2015         9     A         3     T
#> 2016         4     A         2     T
#> 2017         3     A         1     T
#> 2018         2     A         0     T
#> 2019         9     A         1     T
#> 2020         1     A         7     T
#> 2021         3     A         3     T
#> 2022         6     A         2     T
#> 2023         4     A         2     T
#> 2024         1     A         1     T
#> 2025        10     A         4     T
#> 2026         2     A         2     T
#> 2027         7     A         3     T
#> 2028         4     A         4     T
#> 2029         2     A         2     T
#> 2030         3     A         1     T
#> 2031        11     A         1     T
#> 2032         3     A         1     T
#> 2033         7     A         1     T
#> 2034         4     A         2     T
#> 2035        15     A         5     T
#> 2036         9     A         3     T
#> 2037        37     A        11     T
#> 2038         7     A        19     T
#> 2039        26     A        10     T
#> 2040        18     A        10     T
#> 2041         6     A         0     G
#> 2042         4     A         0     G
#> 2043         9     A         1     G
#> 2044        10     A         0     G
#> 2045         8     A         2     G
#> 2046         6     A         0     G
#> 2047         6     A         0     G
#> 2048         2     A         0     G
#> 2049         8     A         6     G
#> 2050         8     A         0     G
#> 2051        12     A         0     G
#> 2052         8     A         0     G
#> 2053         4     A         0     G
#> 2054         2     A         0     G
#> 2055         6     A         6     G
#> 2056         8     A         0     G
#> 2057        11     A         1     G
#> 2058         8     A         0     G
#> 2059         4     A         0     G
#> 2060         1     A         1     G
#> 2061        12     A         2     G
#> 2062         8     A         0     G
#> 2063         7     A         1     G
#> 2064         8     A         0     G
#> 2065        20     A         0     G
#> 2066         9     A         1     G
#> 2067        35     A        15     G
#> 2068        34     A         0     G
#> 2069        38     A         4     G
#> 2070        30     A         0     G
#> 2071         4     A         2     C
#> 2072         2     A         2     C
#> 2073         9     A         1     C
#> 2074        10     A         0     C
#> 2075         6     A         4     C
#> 2076         2     A         2     C
#> 2077         3     A         1     C
#> 2078         1     A         1     C
#> 2079         8     A         4     C
#> 2080         6     A         0     C
#> 2081         7     A         3     C
#> 2082         6     A         0     C
#> 2083         4     A         2     C
#> 2084         0     A         0     C
#> 2085         7     A         3     C
#> 2086         8     A         0     C
#> 2087        10     A         2     C
#> 2088         6     A         2     C
#> 2089         2     A         2     C
#> 2090         1     A         3     C
#> 2091         7     A         3     C
#> 2092         8     A         0     C
#> 2093        10     A         2     C
#> 2094         3     A         3     C
#> 2095        13     A         7     C
#> 2096         4     A         6     C
#> 2097        31     A        11     C
#> 2098        32     A         0     C
#> 2099        33     A        11     C
#> 2100        17     A         7     C
#> 2101         6     A         0     T
#> 2102         4     A         0     T
#> 2103        12     A         0     T
#> 2104         6     A         2     T
#> 2105        11     A         1     T
#> 2106         6     A         0     T
#> 2107         6     A         0     T
#> 2108         2     A         0     T
#> 2109        14     A         0     T
#> 2110         7     A         1     T
#> 2111        11     A         1     T
#> 2112         6     A         2     T
#> 2113         6     A         0     T
#> 2114         2     A         0     T
#> 2115        13     A         1     T
#> 2116         8     A         0     T
#> 2117        11     A         1     T
#> 2118         8     A         0     T
#> 2119         3     A         1     T
#> 2120         4     A         0     T
#> 2121        14     A         0     T
#> 2122         5     A         3     T
#> 2123         9     A         3     T
#> 2124         8     A         0     T
#> 2125        21     A         1     T
#> 2126        12     A         0     T
#> 2127        53     A         1     T
#> 2128        26     A         6     T
#> 2129        42     A         6     T
#> 2130        28     A         2     T
#> 2131         3     A         3     G
#> 2132         3     A         1     G
#> 2133         9     A         3     G
#> 2134         3     A         7     G
#> 2135         9     A         3     G
#> 2136         1     A         3     G
#> 2137         0     A         4     G
#> 2138         1     A         1     G
#> 2139         7     A         5     G
#> 2140         3     A         5     G
#> 2141         9     A         3     G
#> 2142         1     A         5     G
#> 2143         3     A         3     G
#> 2144         0     A         2     G
#> 2145         6     A         4     G
#> 2146         4     A         2     G
#> 2147         7     A         1     G
#> 2148         3     A         5     G
#> 2149         1     A         3     G
#> 2150         1     A         3     G
#> 2151        11     A         3     G
#> 2152         1     A         5     G
#> 2153         9     A         3     G
#> 2154         2     A         6     G
#> 2155         7     A        13     G
#> 2156         5     A         7     G
#> 2157        33     A        15     G
#> 2158        11     A        19     G
#> 2159        34     A        10     G
#> 2160         7     A        19     G
#> 2161         4     T         2     G
#> 2162         4     T         0     G
#> 2163         9     T         1     G
#> 2164         6     T         2     G
#> 2165         7     T         1     G
#> 2166         5     T         1     G
#> 2167         4     T         2     G
#> 2168         1     T         1     G
#> 2169         6     T         2     G
#> 2170         6     T         2     G
#> 2171        11     T         1     G
#> 2172         6     T         2     G
#> 2173         4     T         2     G
#> 2174         2     T         0     G
#> 2175        14     T         0     G
#> 2176         4     T         2     G
#> 2177         8     T         2     G
#> 2178         5     T         1     G
#> 2179         2     T         0     G
#> 2180         2     T         2     G
#> 2181        11     T         3     G
#> 2182         5     T         3     G
#> 2183        11     T         1     G
#> 2184         5     T         3     G
#> 2185        14     T         6     G
#> 2186         9     T         3     G
#> 2187        40     T         6     G
#> 2188        21     T         9     G
#> 2189        37     T         5     G
#> 2190        21     T         7     G
#> 2191         5     G         1     T
#> 2192         3     G         1     T
#> 2193         5     G         1     T
#> 2194        10     G         0     T
#> 2195         8     G         0     T
#> 2196         6     G         0     T
#> 2197         6     G         0     T
#> 2198         2     G         0     T
#> 2199         9     G         5     T
#> 2200         8     G         0     T
#> 2201         3     G         3     T
#> 2202         6     G         2     T
#> 2203         2     G         2     T
#> 2204         2     G         0     T
#> 2205         5     G         3     T
#> 2206         8     G         0     T
#> 2207         7     G         3     T
#> 2208         6     G         0     T
#> 2209         4     G         0     T
#> 2210         2     G         0     T
#> 2211         7     G         3     T
#> 2212         8     G         0     T
#> 2213         8     G         2     T
#> 2214         5     G         1     T
#> 2215        17     G         3     T
#> 2216         9     G         1     T
#> 2217        26     G        12     T
#> 2218        34     G         0     T
#> 2219        26     G         8     T
#> 2220        23     G         3     T
#> 2221         3     A         3     C
#> 2222         1     A         1     C
#> 2223         6     A         4     C
#> 2224         6     A         2     C
#> 2225         5     A         5     C
#> 2226         1     A         5     C
#> 2227         2     A         4     C
#> 2228         1     A         1     C
#> 2229         5     A         3     C
#> 2230         4     A         4     C
#> 2231         6     A         4     C
#> 2232         2     A         6     C
#> 2233         1     A         1     C
#> 2234         0     A         2     C
#> 2235         8     A         6     C
#> 2236         6     A         2     C
#> 2237         4     A         6     C
#> 2238         4     A         4     C
#> 2239         2     A         2     C
#> 2240         3     A         1     C
#> 2241         7     A         7     C
#> 2242         6     A         2     C
#> 2243         6     A         4     C
#> 2244         2     A         6     C
#> 2245         8     A        10     C
#> 2246         5     A         5     C
#> 2247        26     A        20     C
#> 2248        22     A        10     C
#> 2249        21     A        19     C
#> 2250         9     A        21     C
#> 2251         2     C         0     T
#> 2252         2     C         0     T
#> 2253         6     C         4     T
#> 2254         0     C         4     T
#> 2255         7     C         1     T
#> 2256         3     C         3     T
#> 2257         3     C         1     T
#> 2258         0     C         0     T
#> 2259         5     C         5     T
#> 2260         0     C         4     T
#> 2261         3     C         1     T
#> 2262         6     C         2     T
#> 2263         1     C         3     T
#> 2264         0     C         0     T
#> 2265         3     C         1     T
#> 2266         0     C         4     T
#> 2267         4     C         0     T
#> 2268         3     C         3     T
#> 2269         1     C         1     T
#> 2270         3     C         1     T
#> 2271         7     C         1     T
#> 2272         0     C         4     T
#> 2273         7     C         1     T
#> 2274         4     C         4     T
#> 2275         7     C         5     T
#> 2276         5     C         1     T
#> 2277        21     C        11     T
#> 2278         0     C        16     T
#> 2279        21     C         3     T
#> 2280        16     C        12     T
#> 2281         4     A         2     T
#> 2282         4     A         0     T
#> 2283         7     A         1     T
#> 2284         7     A         3     T
#> 2285        10     A         2     T
#> 2286         3     A         3     T
#> 2287         3     A         1     T
#> 2288         1     A         1     T
#> 2289        10     A         2     T
#> 2290         3     A         3     T
#> 2291         8     A         2     T
#> 2292         7     A         1     T
#> 2293         4     A         2     T
#> 2294         2     A         0     T
#> 2295        11     A         3     T
#> 2296         4     A         4     T
#> 2297         4     A         2     T
#> 2298         5     A         1     T
#> 2299         1     A         1     T
#> 2300         4     A         0     T
#> 2301        10     A         2     T
#> 2302         5     A         1     T
#> 2303         9     A         1     T
#> 2304         7     A         1     T
#> 2305        12     A         6     T
#> 2306        11     A         1     T
#> 2307        38     A         8     T
#> 2308        19     A        11     T
#> 2309        31     A         7     T
#> 2310        22     A         6     T
#> 2311         4     G         0     T
#> 2312         4     G         0     T
#> 2313         7     G         1     T
#> 2314        10     G         0     T
#> 2315         9     G         3     T
#> 2316         6     G         0     T
#> 2317         6     G         0     T
#> 2318         2     G         0     T
#> 2319         8     G         4     T
#> 2320         8     G         0     T
#> 2321        10     G         0     T
#> 2322         8     G         0     T
#> 2323         6     G         0     T
#> 2324         2     G         0     T
#> 2325         8     G         6     T
#> 2326         8     G         0     T
#> 2327        11     G         1     T
#> 2328         8     G         0     T
#> 2329         4     G         0     T
#> 2330         1     G         3     T
#> 2331         8     G         2     T
#> 2332         8     G         0     T
#> 2333         9     G         1     T
#> 2334         8     G         0     T
#> 2335        20     G         0     T
#> 2336         9     G         3     T
#> 2337        31     G        13     T
#> 2338        34     G         0     T
#> 2339        39     G         5     T
#> 2340        30     G         0     T
#> 2341         5     T         1     C
#> 2342         4     T         0     C
#> 2343         9     T         1     C
#> 2344         6     T         4     C
#> 2345        12     T         0     C
#> 2346         5     T         1     C
#> 2347         5     T         1     C
#> 2348         2     T         0     C
#> 2349        11     T         3     C
#> 2350         7     T         1     C
#> 2351         6     T         2     C
#> 2352         2     T         2     C
#> 2353         4     T         2     C
#> 2354         2     T         0     C
#> 2355        14     T         0     C
#> 2356         7     T         1     C
#> 2357        11     T         1     C
#> 2358         6     T         0     C
#> 2359         4     T         0     C
#> 2360         4     T         0     C
#> 2361        12     T         2     C
#> 2362         6     T         2     C
#> 2363        11     T         1     C
#> 2364         8     T         0     C
#> 2365        18     T         4     C
#> 2366        12     T         0     C
#> 2367        46     T         6     C
#> 2368        26     T         8     C
#> 2369        40     T         4     C
#> 2370        21     T         3     C
#> 2371         3     T         3     C
#> 2372         3     T         1     C
#> 2373         8     T         2     C
#> 2374         5     T         5     C
#> 2375         8     T         4     C
#> 2376         4     T         2     C
#> 2377         4     T         2     C
#> 2378         2     T         0     C
#> 2379         5     T         5     C
#> 2380         3     T         5     C
#> 2381        10     T         2     C
#> 2382         6     T         2     C
#> 2383         6     T         0     C
#> 2384         2     T         0     C
#> 2385         9     T         3     C
#> 2386         4     T         4     C
#> 2387         7     T         3     C
#> 2388         6     T         2     C
#> 2389         4     T         0     C
#> 2390         3     T         1     C
#> 2391         8     T         4     C
#> 2392         3     T         1     C
#> 2393         9     T         3     C
#> 2394         3     T         1     C
#> 2395        17     T         5     C
#> 2396        10     T         2     C
#> 2397        30     T        14     C
#> 2398        15     T        15     C
#> 2399        34     T        12     C
#> 2400        19     T         7     C
#> 2401         4     A         2     T
#> 2402         1     A         1     T
#> 2403         6     A         4     T
#> 2404         4     A         6     T
#> 2405         8     A         4     T
#> 2406         2     A         2     T
#> 2407         1     A         1     T
#> 2408         1     A         1     T
#> 2409         7     A         3     T
#> 2410         4     A         2     T
#> 2411         8     A         2     T
#> 2412         6     A         2     T
#> 2413         4     A         2     T
#> 2414         1     A         1     T
#> 2415         4     A         6     T
#> 2416         6     A         2     T
#> 2417         7     A         3     T
#> 2418         5     A         3     T
#> 2419         3     A         1     T
#> 2420         1     A         3     T
#> 2421         8     A         4     T
#> 2422         4     A         4     T
#> 2423         7     A         5     T
#> 2424         4     A         4     T
#> 2425        12     A         6     T
#> 2426         4     A         6     T
#> 2427        25     A        17     T
#> 2428        18     A        14     T
#> 2429        30     A        14     T
#> 2430        17     A        11     T
#> 2431         5     T         1     C
#> 2432         2     T         0     C
#> 2433        11     T         1     C
#> 2434        10     T         0     C
#> 2435        12     T         0     C
#> 2436         4     T         2     C
#> 2437         6     T         0     C
#> 2438         0     T         0     C
#> 2439        11     T         3     C
#> 2440         8     T         0     C
#> 2441         6     T         0     C
#> 2442         8     T         0     C
#> 2443         3     T         1     C
#> 2444         2     T         0     C
#> 2445        11     T         1     C
#> 2446         8     T         0     C
#> 2447         8     T         0     C
#> 2448         5     T         1     C
#> 2449         4     T         0     C
#> 2450         2     T         0     C
#> 2451        14     T         0     C
#> 2452         8     T         0     C
#> 2453        10     T         2     C
#> 2454         4     T         2     C
#> 2455        18     T         2     C
#> 2456         6     T         0     C
#> 2457        47     T         5     C
#> 2458        34     T         0     C
#> 2459        36     T         2     C
#> 2460        21     T         5     C
#> 2461         5     T         1     A
#> 2462         2     T         2     A
#> 2463        12     T         0     A
#> 2464        10     T         0     A
#> 2465        11     T         1     A
#> 2466         5     T         1     A
#> 2467         4     T         2     A
#> 2468         2     T         0     A
#> 2469        14     T         0     A
#> 2470         8     T         0     A
#> 2471        12     T         0     A
#> 2472         8     T         0     A
#> 2473         5     T         1     A
#> 2474         1     T         1     A
#> 2475        14     T         0     A
#> 2476         6     T         0     A
#> 2477        12     T         0     A
#> 2478         6     T         0     A
#> 2479         4     T         0     A
#> 2480         2     T         2     A
#> 2481        12     T         0     A
#> 2482         8     T         0     A
#> 2483        12     T         0     A
#> 2484         8     T         0     A
#> 2485        18     T         4     A
#> 2486         7     T         5     A
#> 2487        52     T         0     A
#> 2488        32     T         0     A
#> 2489        47     T         1     A
#> 2490        27     T         1     A
#> 2491         5     C         1     T
#> 2492         3     C         1     T
#> 2493         3     C         5     T
#> 2494         0     C         8     T
#> 2495         8     C         4     T
#> 2496         4     C         2     T
#> 2497         4     C         0     T
#> 2498         2     C         0     T
#> 2499         8     C         4     T
#> 2500         0     C         6     T
#> 2501         5     C         3     T
#> 2502         3     C         5     T
#> 2503         2     C         4     T
#> 2504         2     C         0     T
#> 2505        11     C         3     T
#> 2506         0     C         8     T
#> 2507         6     C         2     T
#> 2508         2     C         4     T
#> 2509         3     C         1     T
#> 2510         1     C         3     T
#> 2511         4     C         2     T
#> 2512         0     C         8     T
#> 2513         9     C         1     T
#> 2514         5     C         1     T
#> 2515        14     C         6     T
#> 2516         8     C         4     T
#> 2517        26     C        14     T
#> 2518         0     C        30     T
#> 2519        28     C        10     T
#> 2520        14     C        12     T
#> 2521         3     G         3     A
#> 2522         1     G         3     A
#> 2523         7     G         5     A
#> 2524         0     G         8     A
#> 2525         5     G         3     A
#> 2526         4     G         2     A
#> 2527         1     G         3     A
#> 2528         0     G         0     A
#> 2529         5     G         1     A
#> 2530         0     G         8     A
#> 2531        10     G         2     A
#> 2532         5     G         3     A
#> 2533         3     G         3     A
#> 2534         1     G         1     A
#> 2535         9     G         5     A
#> 2536         1     G         5     A
#> 2537         8     G         2     A
#> 2538         4     G         4     A
#> 2539         0     G         4     A
#> 2540         2     G         2     A
#> 2541         8     G         2     A
#> 2542         2     G         4     A
#> 2543         8     G         4     A
#> 2544         4     G         4     A
#> 2545         7     G        13     A
#> 2546         4     G         6     A
#> 2547        29     G        13     A
#> 2548         3     G        25     A
#> 2549        31     G        11     A
#> 2550        17     G        13     A
#> 2551         4     A         2     G
#> 2552         4     A         0     G
#> 2553         6     A         2     G
#> 2554         4     A         2     G
#> 2555         6     A         2     G
#> 2556         3     A         3     G
#> 2557         3     A         1     G
#> 2558         1     A         1     G
#> 2559         6     A         2     G
#> 2560         2     A         2     G
#> 2561        10     A         2     G
#> 2562         5     A         1     G
#> 2563         2     A         2     G
#> 2564         2     A         0     G
#> 2565        11     A         3     G
#> 2566         4     A         4     G
#> 2567         6     A         2     G
#> 2568         4     A         2     G
#> 2569         1     A         1     G
#> 2570         4     A         0     G
#> 2571         8     A         2     G
#> 2572         3     A         1     G
#> 2573         9     A         1     G
#> 2574         7     A         1     G
#> 2575        10     A         6     G
#> 2576        11     A         1     G
#> 2577        31     A         9     G
#> 2578        13     A         9     G
#> 2579        31     A         7     G
#> 2580        19     A         7     G
#> 2581         5     T         1     C
#> 2582         3     T         1     C
#> 2583         8     T         4     C
#> 2584        10     T         0     C
#> 2585         6     T         4     C
#> 2586         3     T         3     C
#> 2587         3     T         3     C
#> 2588         2     T         0     C
#> 2589         8     T         4     C
#> 2590         8     T         0     C
#> 2591         8     T         4     C
#> 2592         6     T         2     C
#> 2593         1     T         5     C
#> 2594         2     T         0     C
#> 2595        10     T         2     C
#> 2596         6     T         0     C
#> 2597         8     T         2     C
#> 2598         4     T         4     C
#> 2599         2     T         2     C
#> 2600         4     T         0     C
#> 2601        12     T         2     C
#> 2602         8     T         0     C
#> 2603        10     T         2     C
#> 2604         3     T         3     C
#> 2605        11     T        11     C
#> 2606        11     T         1     C
#> 2607        38     T        12     C
#> 2608        32     T         0     C
#> 2609        32     T        12     C
#> 2610        16     T        12     C
#> 2611         2     A         4     G
#> 2612         3     A         1     G
#> 2613         6     A         2     G
#> 2614         0     A         8     G
#> 2615         3     A         5     G
#> 2616         2     A         2     G
#> 2617         2     A         2     G
#> 2618         1     A         1     G
#> 2619         4     A         6     G
#> 2620         0     A         6     G
#> 2621         6     A         2     G
#> 2622         4     A         2     G
#> 2623         3     A         1     G
#> 2624         2     A         0     G
#> 2625         8     A         4     G
#> 2626         0     A         8     G
#> 2627         7     A         3     G
#> 2628         5     A         1     G
#> 2629         1     A         1     G
#> 2630         3     A         1     G
#> 2631         9     A         1     G
#> 2632         0     A         6     G
#> 2633         3     A         3     G
#> 2634         5     A         3     G
#> 2635         8     A         8     G
#> 2636         9     A         3     G
#> 2637        27     A        13     G
#> 2638         0     A        28     G
#> 2639        19     A        13     G
#> 2640        16     A         8     G
#> 2641         4     C         0     G
#> 2642         0     C         0     G
#> 2643         7     C         3     G
#> 2644         4     C         4     G
#> 2645         6     C         4     G
#> 2646         1     C         3     G
#> 2647         3     C         3     G
#> 2648         0     C         2     G
#> 2649         3     C         5     G
#> 2650         4     C         2     G
#> 2651         6     C         4     G
#> 2652         1     C         5     G
#> 2653         4     C         2     G
#> 2654         2     C         0     G
#> 2655         6     C         4     G
#> 2656         4     C         2     G
#> 2657         6     C         2     G
#> 2658         2     C         4     G
#> 2659         1     C         1     G
#> 2660         2     C         0     G
#> 2661         7     C         1     G
#> 2662         0     C         4     G
#> 2663         5     C         3     G
#> 2664         6     C         2     G
#> 2665        12     C         6     G
#> 2666         4     C         2     G
#> 2667        23     C        13     G
#> 2668        12     C        12     G
#> 2669        23     C        13     G
#> 2670        10     C        14     G
#> 2671         5     C         1     T
#> 2672         4     C         0     T
#> 2673        10     C         0     T
#> 2674        10     C         0     T
#> 2675        10     C         0     T
#> 2676         4     C         2     T
#> 2677         4     C         2     T
#> 2678         2     C         0     T
#> 2679        14     C         0     T
#> 2680         8     C         0     T
#> 2681        12     C         0     T
#> 2682         6     C         0     T
#> 2683         6     C         0     T
#> 2684         2     C         0     T
#> 2685        14     C         0     T
#> 2686         8     C         0     T
#> 2687        12     C         0     T
#> 2688         6     C         0     T
#> 2689         4     C         0     T
#> 2690         4     C         0     T
#> 2691        14     C         0     T
#> 2692         8     C         0     T
#> 2693        11     C         1     T
#> 2694         8     C         0     T
#> 2695        19     C         3     T
#> 2696        12     C         0     T
#> 2697        52     C         0     T
#> 2698        34     C         0     T
#> 2699        45     C         1     T
#> 2700        24     C         2     T
#> 2701         6     C         0     T
#> 2702         1     C         1     T
#> 2703         5     C         5     T
#> 2704         7     C         1     T
#> 2705         4     C         8     T
#> 2706         2     C         2     T
#> 2707         3     C         1     T
#> 2708         0     C         2     T
#> 2709         7     C         5     T
#> 2710         7     C         1     T
#> 2711         6     C         6     T
#> 2712         5     C         3     T
#> 2713         4     C         0     T
#> 2714         2     C         0     T
#> 2715         2     C         6     T
#> 2716         8     C         0     T
#> 2717         2     C         4     T
#> 2718         5     C         3     T
#> 2719         2     C         2     T
#> 2720         2     C         0     T
#> 2721         7     C         5     T
#> 2722         7     C         1     T
#> 2723         6     C         4     T
#> 2724         6     C         2     T
#> 2725        15     C         3     T
#> 2726         5     C         3     T
#> 2727        21     C        21     T
#> 2728        29     C         3     T
#> 2729        18     C        22     T
#> 2730        18     C        10     T
#> 2731         2     G         2     A
#> 2732         1     G         1     A
#> 2733         7     G         3     A
#> 2734         2     G         4     A
#> 2735         4     G         0     A
#> 2736         4     G         0     A
#> 2737         0     G         2     A
#> 2738         2     G         0     A
#> 2739         3     G         3     A
#> 2740         1     G         3     A
#> 2741         4     G         2     A
#> 2742         2     G         2     A
#> 2743         4     G         2     A
#> 2744         2     G         0     A
#> 2745         4     G         2     A
#> 2746         3     G         5     A
#> 2747         6     G         0     A
#> 2748         5     G         1     A
#> 2749         1     G         1     A
#> 2750         3     G         1     A
#> 2751         4     G         2     A
#> 2752         2     G         2     A
#> 2753         4     G         2     A
#> 2754         5     G         1     A
#> 2755         7     G         7     A
#> 2756         8     G         2     A
#> 2757        18     G        10     A
#> 2758         8     G        14     A
#> 2759        18     G         4     A
#> 2760        16     G         4     A
#> 2761         5     T         1     C
#> 2762         4     T         0     C
#> 2763        12     T         0     C
#> 2764         5     T         1     C
#> 2765         5     T         3     C
#> 2766         2     T         2     C
#> 2767         5     T         1     C
#> 2768         1     T         1     C
#> 2769        12     T         2     C
#> 2770         6     T         2     C
#> 2771         9     T         1     C
#> 2772         6     T         2     C
#> 2773         5     T         1     C
#> 2774         2     T         0     C
#> 2775        11     T         1     C
#> 2776         6     T         0     C
#> 2777         9     T         1     C
#> 2778         7     T         1     C
#> 2779         4     T         0     C
#> 2780         2     T         2     C
#> 2781        14     T         0     C
#> 2782         7     T         1     C
#> 2783         8     T         2     C
#> 2784         6     T         2     C
#> 2785        19     T         3     C
#> 2786         9     T         3     C
#> 2787        49     T         3     C
#> 2788        24     T         4     C
#> 2789        31     T         7     C
#> 2790        21     T         7     C
#> 2791         6     C         0     G
#> 2792         4     C         0     G
#> 2793        10     C         0     G
#> 2794         6     C         4     G
#> 2795        12     C         0     G
#> 2796         6     C         0     G
#> 2797         2     C         0     G
#> 2798         2     C         0     G
#> 2799         9     C         1     G
#> 2800         4     C         2     G
#> 2801        11     C         1     G
#> 2802         7     C         1     G
#> 2803         4     C         0     G
#> 2804         2     C         0     G
#> 2805        11     C         1     G
#> 2806         5     C         1     G
#> 2807        10     C         0     G
#> 2808         6     C         0     G
#> 2809         4     C         0     G
#> 2810         2     C         0     G
#> 2811        11     C         1     G
#> 2812         2     C         2     G
#> 2813        10     C         0     G
#> 2814         6     C         2     G
#> 2815        16     C         0     G
#> 2816        10     C         0     G
#> 2817        41     C         3     G
#> 2818        17     C         9     G
#> 2819        43     C         1     G
#> 2820        25     C         3     G
#> 2821         5     G         1     T
#> 2822         2     G         2     T
#> 2823        11     G         1     T
#> 2824        10     G         0     T
#> 2825         9     G         1     T
#> 2826         5     G         1     T
#> 2827         5     G         1     T
#> 2828         2     G         0     T
#> 2829         9     G         3     T
#> 2830         8     G         0     T
#> 2831         9     G         1     T
#> 2832         6     G         0     T
#> 2833         5     G         1     T
#> 2834         2     G         0     T
#> 2835         7     G         3     T
#> 2836         8     G         0     T
#> 2837         8     G         2     T
#> 2838         4     G         2     T
#> 2839         1     G         1     T
#> 2840         2     G         2     T
#> 2841         9     G         1     T
#> 2842         8     G         0     T
#> 2843         6     G         2     T
#> 2844         5     G         1     T
#> 2845        16     G         4     T
#> 2846         8     G         4     T
#> 2847        36     G         8     T
#> 2848        34     G         0     T
#> 2849        32     G         6     T
#> 2850        20     G         4     T
#> 2851         5     C         1     A
#> 2852         4     C         0     A
#> 2853         6     C         2     A
#> 2854        10     C         0     A
#> 2855        11     C         1     A
#> 2856         4     C         0     A
#> 2857         4     C         0     A
#> 2858         2     C         0     A
#> 2859         9     C         3     A
#> 2860         8     C         0     A
#> 2861        11     C         1     A
#> 2862         6     C         0     A
#> 2863         4     C         0     A
#> 2864         1     C         1     A
#> 2865        10     C         2     A
#> 2866         8     C         0     A
#> 2867         9     C         1     A
#> 2868         6     C         2     A
#> 2869         2     C         0     A
#> 2870         4     C         0     A
#> 2871        13     C         1     A
#> 2872         8     C         0     A
#> 2873        11     C         1     A
#> 2874         8     C         0     A
#> 2875        15     C         1     A
#> 2876        11     C         1     A
#> 2877        38     C         8     A
#> 2878        34     C         0     A
#> 2879        42     C         4     A
#> 2880        24     C         2     A
#> 2881         6     G         0     A
#> 2882         4     G         0     A
#> 2883         9     G         1     A
#> 2884        10     G         0     A
#> 2885        11     G         1     A
#> 2886         5     G         1     A
#> 2887         6     G         0     A
#> 2888         2     G         0     A
#> 2889         7     G         1     A
#> 2890         8     G         0     A
#> 2891         7     G         1     A
#> 2892         8     G         0     A
#> 2893         6     G         0     A
#> 2894         2     G         0     A
#> 2895        12     G         0     A
#> 2896         8     G         0     A
#> 2897        10     G         0     A
#> 2898         8     G         0     A
#> 2899         4     G         0     A
#> 2900         3     G         1     A
#> 2901        12     G         0     A
#> 2902         8     G         0     A
#> 2903        10     G         0     A
#> 2904         8     G         0     A
#> 2905        22     G         0     A
#> 2906        11     G         1     A
#> 2907        40     G         2     A
#> 2908        34     G         0     A
#> 2909        38     G         2     A
#> 2910        29     G         1     A
#> 2911         4     G         0     T
#> 2912         2     G         2     T
#> 2913        11     G         1     T
#> 2914         5     G         5     T
#> 2915         8     G         4     T
#> 2916         5     G         1     T
#> 2917         5     G         1     T
#> 2918         2     G         0     T
#> 2919        13     G         1     T
#> 2920         5     G         3     T
#> 2921         9     G         1     T
#> 2922         4     G         4     T
#> 2923         4     G         2     T
#> 2924         0     G         0     T
#> 2925        11     G         1     T
#> 2926         5     G         3     T
#> 2927        10     G         0     T
#> 2928         3     G         1     T
#> 2929         1     G         1     T
#> 2930         3     G         1     T
#> 2931        12     G         2     T
#> 2932         4     G         2     T
#> 2933         8     G         4     T
#> 2934         5     G         3     T
#> 2935        14     G         4     T
#> 2936         7     G         3     T
#> 2937        47     G         5     T
#> 2938        19     G        13     T
#> 2939        35     G         9     T
#> 2940        17     G         9     T
#> 2941         1     G         3     A
#> 2942         2     G         2     A
#> 2943         1     G         7     A
#> 2944         8     G         0     A
#> 2945         6     G         6     A
#> 2946         1     G         3     A
#> 2947         5     G         1     A
#> 2948         2     G         0     A
#> 2949         5     G         5     A
#> 2950         8     G         0     A
#> 2951         3     G         5     A
#> 2952         3     G         3     A
#> 2953         2     G         4     A
#> 2954         0     G         2     A
#> 2955         6     G         6     A
#> 2956         4     G         0     A
#> 2957         4     G         6     A
#> 2958         5     G         1     A
#> 2959         1     G         3     A
#> 2960         3     G         1     A
#> 2961         6     G         6     A
#> 2962         8     G         0     A
#> 2963         6     G         4     A
#> 2964         2     G         4     A
#> 2965         9     G        11     A
#> 2966         7     G         5     A
#> 2967        18     G        24     A
#> 2968        28     G         0     A
#> 2969        19     G        21     A
#> 2970        11     G        11     A
#> 2971         1     A         5     T
#> 2972         2     A         2     T
#> 2973         8     A         2     T
#> 2974         4     A         2     T
#> 2975         9     A         3     T
#> 2976         1     A         5     T
#> 2977         3     A         1     T
#> 2978         0     A         0     T
#> 2979         8     A         2     T
#> 2980         4     A         2     T
#> 2981         6     A         2     T
#> 2982         4     A         4     T
#> 2983         3     A         3     T
#> 2984         2     A         0     T
#> 2985         8     A         2     T
#> 2986         4     A         2     T
#> 2987         6     A         4     T
#> 2988         3     A         3     T
#> 2989         1     A         3     T
#> 2990         1     A         1     T
#> 2991         4     A         6     T
#> 2992         3     A         1     T
#> 2993         6     A         4     T
#> 2994         3     A         5     T
#> 2995         8     A        12     T
#> 2996         5     A         3     T
#> 2997        28     A        12     T
#> 2998        15     A         7     T
#> 2999        27     A        13     T
#> 3000        11     A        17     T

# HWE with family-wise error correction
if (FALSE) { # \dontrun{
x <- calc_hwe(stickSNPs, facets = c("pop", "pop.fam"))
get.snpR.stats(x, c("pop", "pop.fam"))
} # }
```
