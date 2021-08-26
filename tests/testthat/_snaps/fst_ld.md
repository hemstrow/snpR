# correct genepop

    Code
      tdfst
    Output
      $pairwise
             group position facet comparison     fst
      1   groupXIX    67921   pop    ASP~PAL  0.3672
      2    groupIX   100382   pop    ASP~PAL  0.2189
      3     groupX   101821   pop    ASP~PAL -0.1111
      4     groupI   182629   pop    ASP~PAL  0.2704
      5    groupIV   194493   pop    ASP~PAL -0.1798
      6   groupVII   208605   pop    ASP~PAL  0.1667
      7    groupVI   212436   pop    ASP~PAL -0.1296
      8   groupXXI   239621   pop    ASP~PAL -0.0353
      9  groupXVII   250065   pop    ASP~PAL  0.1667
      10    groupX   267302   pop    ASP~PAL -0.0817
      
      $weighted.means
        facet subfacet snp.facet snp.subfacet weighted_mean_fst
      1   pop  ASP~PAL     .base        .base            0.0494
      
      $fst.matrix
      $fst.matrix$pop
          p1    PAL
      1: ASP 0.0494
      
      

# correct traditional ld

    Code
      tdld
    Output
      $LD
      $LD$prox
           s1_group s1_position s1_.snp.id  s2_group s2_position s2_.snp.id proximity         rsq     Dprime      pval sample.facet
      1    groupXIX       67921          1   groupIX      100382          2     32461         NaN 0.00000000       NaN        .base
      1.1  groupXIX       67921          1    groupX      101821          3     33900 0.064935065 1.00000000 0.6102987        .base
      1.2  groupXIX       67921          1    groupI      182629          4    114708 0.042735043 1.00000000 0.6792776        .base
      1.3  groupXIX       67921          1   groupIV      194493          5    126572 0.151515152 1.00000000 0.4362749        .base
      1.4  groupXIX       67921          1  groupVII      208605          6    140684 0.100000000 1.00000000 0.5270893        .base
      1.5  groupXIX       67921          1   groupVI      212436          7    144515 0.107142857 1.00000000 0.5126908        .base
      1.6  groupXIX       67921          1  groupXXI      239621          8    171700 0.002267574 0.04761905 0.9241257        .base
      1.7  groupXIX       67921          1 groupXVII      250065          9    182144 0.030769231 1.00000000 0.7257210        .base
      1.8  groupXIX       67921          1    groupX      267302         10    199381 0.066666667 0.33333333 0.6055766        .base
      2     groupIX      100382          2    groupX      101821          3      1439 0.025000000 1.00000000 0.7518296        .base
      2.1   groupIX      100382          2    groupI      182629          4     82247 0.020979021 1.00000000 0.7720590        .base
      2.2   groupIX      100382          2   groupIV      194493          5     94111 0.076923077 1.00000000 0.5790997        .base
      2.3   groupIX      100382          2  groupVII      208605          6    108223 0.027777778 1.00000000 0.7388827        .base
      2.4   groupIX      100382          2   groupVI      212436          7    112054         NaN 0.00000000       NaN        .base
      2.5   groupIX      100382          2  groupXXI      239621          8    139239 0.127272727 1.00000000 0.4755327        .base
      2.6   groupIX      100382          2 groupXVII      250065          9    149683 0.020408163 1.00000000 0.7750970        .base
      2.7   groupIX      100382          2    groupX      267302         10    166920 0.020979021 1.00000000 0.7720590        .base
      3      groupX      101821          3    groupI      182629          4     80808 0.005917160 1.00000000 0.8777310        .base
      3.1    groupX      101821          3   groupIV      194493          5     92672 0.047619048 0.33333333 0.6625206        .base
      3.2    groupX      101821          3  groupVII      208605          6    106784 0.009523810 1.00000000 0.8452520        .base
      3.3    groupX      101821          3   groupVI      212436          7    110615 0.062500000 1.00000000 0.6170751        .base
      3.4    groupX      101821          3  groupXXI      239621          8    137800         NaN 0.00000000       NaN        .base
      3.5    groupX      101821          3 groupXVII      250065          9    148244 0.020408163 1.00000000 0.7750970        .base
      3.6    groupX      101821          3    groupX      267302         10    165481 0.022222222 1.00000000 0.7655945        .base
      4      groupI      182629          4   groupIV      194493          5     11864          NA 0.00000000        NA        .base
      4.1    groupI      182629          4  groupVII      208605          6     25976          NA 0.00000000        NA        .base
      4.2    groupI      182629          4   groupVI      212436          7     29807 0.102564103 1.00000000 0.5218394        .base
      4.3    groupI      182629          4  groupXXI      239621          8     56992 0.042735043 1.00000000 0.6792776        .base
      4.4    groupI      182629          4 groupXVII      250065          9     67436 0.008264463 1.00000000 0.8557254        .base
      4.5    groupI      182629          4    groupX      267302         10     84673 0.018181818 1.00000000 0.7874065        .base
      5     groupIV      194493          5  groupVII      208605          6     14112          NA 0.00000000        NA        .base
      5.1   groupIV      194493          5   groupVI      212436          7     17943 0.074380165 1.00000000 0.5854409        .base
      5.2   groupIV      194493          5  groupXXI      239621          8     45128 0.000000000 0.00000000 1.0000000        .base
      5.3   groupIV      194493          5 groupXVII      250065          9     55572          NA 0.00000000        NA        .base
      5.4   groupIV      194493          5    groupX      267302         10     72809 0.111111111 1.00000000 0.5049851        .base
      6    groupVII      208605          6   groupVI      212436          7      3831 0.181818182 1.00000000 0.3937686        .base
      6.1  groupVII      208605          6  groupXXI      239621          8     31016 0.085714286 1.00000000 0.5581846        .base
      6.2  groupVII      208605          6 groupXVII      250065          9     41460 0.008264463 1.00000000 0.8557254        .base
      6.3  groupVII      208605          6    groupX      267302         10     58697 0.047619048 1.00000000 0.6625206        .base
      7     groupVI      212436          7  groupXXI      239621          8     27185 0.057142857 0.40000000 0.6325851        .base
      7.1   groupVI      212436          7 groupXVII      250065          9     37629 0.030769231 1.00000000 0.7257210        .base
      7.2   groupVI      212436          7    groupX      267302         10     54866 0.109090909 1.00000000 0.5088828        .base
      8    groupXXI      239621          8 groupXVII      250065          9     10444 0.003472222 0.12500000 0.9061856        .base
      8.1  groupXXI      239621          8    groupX      267302         10     27681 0.022222222 0.20000000 0.7655945        .base
      9   groupXVII      250065          9    groupX      267302         10     17237 0.020979021 1.00000000 0.7720590        .base
          sample.subfacet
      1             .base
      1.1           .base
      1.2           .base
      1.3           .base
      1.4           .base
      1.5           .base
      1.6           .base
      1.7           .base
      1.8           .base
      2             .base
      2.1           .base
      2.2           .base
      2.3           .base
      2.4           .base
      2.5           .base
      2.6           .base
      2.7           .base
      3             .base
      3.1           .base
      3.2           .base
      3.3           .base
      3.4           .base
      3.5           .base
      3.6           .base
      4             .base
      4.1           .base
      4.2           .base
      4.3           .base
      4.4           .base
      4.5           .base
      5             .base
      5.1           .base
      5.2           .base
      5.3           .base
      5.4           .base
      6             .base
      6.1           .base
      6.2           .base
      6.3           .base
      7             .base
      7.1           .base
      7.2           .base
      8             .base
      8.1           .base
      9             .base
      
      $LD$matrices
      $LD$matrices$.base
      $LD$matrices$.base$.base
      $LD$matrices$.base$.base$Dprime
             67921 100382 101821 182629    194493 208605 212436     239621 250065    267302
      67921     NA      0      1      1 1.0000000      1      1 0.04761905  1.000 0.3333333
      100382    NA     NA      1      1 1.0000000      1      0 1.00000000  1.000 1.0000000
      101821    NA     NA     NA      1 0.3333333      1      1 0.00000000  1.000 1.0000000
      182629    NA     NA     NA     NA 0.0000000      0      1 1.00000000  1.000 1.0000000
      194493    NA     NA     NA     NA        NA      0      1 0.00000000  0.000 1.0000000
      208605    NA     NA     NA     NA        NA     NA      1 1.00000000  1.000 1.0000000
      212436    NA     NA     NA     NA        NA     NA     NA 0.40000000  1.000 1.0000000
      239621    NA     NA     NA     NA        NA     NA     NA         NA  0.125 0.2000000
      250065    NA     NA     NA     NA        NA     NA     NA         NA     NA 1.0000000
      267302    NA     NA     NA     NA        NA     NA     NA         NA     NA        NA
      
      $LD$matrices$.base$.base$rsq
             67921 100382     101821     182629     194493     208605     212436      239621      250065     267302
      67921     NA    NaN 0.06493506 0.04273504 0.15151515 0.10000000 0.10714286 0.002267574 0.030769231 0.06666667
      100382    NA     NA 0.02500000 0.02097902 0.07692308 0.02777778        NaN 0.127272727 0.020408163 0.02097902
      101821    NA     NA         NA 0.00591716 0.04761905 0.00952381 0.06250000         NaN 0.020408163 0.02222222
      182629    NA     NA         NA         NA         NA         NA 0.10256410 0.042735043 0.008264463 0.01818182
      194493    NA     NA         NA         NA         NA         NA 0.07438017 0.000000000          NA 0.11111111
      208605    NA     NA         NA         NA         NA         NA 0.18181818 0.085714286 0.008264463 0.04761905
      212436    NA     NA         NA         NA         NA         NA         NA 0.057142857 0.030769231 0.10909091
      239621    NA     NA         NA         NA         NA         NA         NA          NA 0.003472222 0.02222222
      250065    NA     NA         NA         NA         NA         NA         NA          NA          NA 0.02097902
      267302    NA     NA         NA         NA         NA         NA         NA          NA          NA         NA
      
      $LD$matrices$.base$.base$pval
             67921 100382    101821    182629    194493    208605    212436    239621    250065    267302
      67921     NA    NaN 0.6102987 0.6792776 0.4362749 0.5270893 0.5126908 0.9241257 0.7257210 0.6055766
      100382    NA     NA 0.7518296 0.7720590 0.5790997 0.7388827       NaN 0.4755327 0.7750970 0.7720590
      101821    NA     NA        NA 0.8777310 0.6625206 0.8452520 0.6170751       NaN 0.7750970 0.7655945
      182629    NA     NA        NA        NA        NA        NA 0.5218394 0.6792776 0.8557254 0.7874065
      194493    NA     NA        NA        NA        NA        NA 0.5854409 1.0000000        NA 0.5049851
      208605    NA     NA        NA        NA        NA        NA 0.3937686 0.5581846 0.8557254 0.6625206
      212436    NA     NA        NA        NA        NA        NA        NA 0.6325851 0.7257210 0.5088828
      239621    NA     NA        NA        NA        NA        NA        NA        NA 0.9061856 0.7655945
      250065    NA     NA        NA        NA        NA        NA        NA        NA        NA 0.7720590
      267302    NA     NA        NA        NA        NA        NA        NA        NA        NA        NA
      
      
      
      
      

# correct ME ld

    Code
      tdldme
    Output
      $LD
      $LD$prox
           s1_group s1_position s1_.snp.id  s2_group s2_position s2_.snp.id proximity         rsq    Dprime      pval sample.facet
      1    groupXIX       67921          1   groupIX      100382          2     32461 0.353994002 0.9999389 0.2340669        .base
      1.1  groupXIX       67921          1    groupX      101821          3     33900 0.064935065 1.0000000 0.6102987        .base
      1.2  groupXIX       67921          1    groupI      182629          4    114708 0.042735043 1.0000000 0.6792776        .base
      1.3  groupXIX       67921          1   groupIV      194493          5    126572 0.151515152 1.0000000 0.4362749        .base
      1.4  groupXIX       67921          1  groupVII      208605          6    140684 0.100000000 1.0000000 0.5270893        .base
      1.5  groupXIX       67921          1   groupVI      212436          7    144515 0.194966106 0.9997940 0.3771826        .base
      1.6  groupXIX       67921          1  groupXXI      239621          8    171700 0.011083805 0.1052797 0.8332312        .base
      1.7  groupXIX       67921          1 groupXVII      250065          9    182144 0.030769231 1.0000000 0.7257210        .base
      1.8  groupXIX       67921          1    groupX      267302         10    199381 0.108329998 0.4094221 0.5103644        .base
      2     groupIX      100382          2    groupX      101821          3      1439 0.025000000 1.0000000 0.7518296        .base
      2.1   groupIX      100382          2    groupI      182629          4     82247 0.020979021 1.0000000 0.7720590        .base
      2.2   groupIX      100382          2   groupIV      194493          5     94111 0.076923077 1.0000000 0.5790997        .base
      2.3   groupIX      100382          2  groupVII      208605          6    108223 0.027777778 1.0000000 0.7388827        .base
      2.4   groupIX      100382          2   groupVI      212436          7    112054 0.259192408 0.9998711 0.3085740        .base
      2.5   groupIX      100382          2  groupXXI      239621          8    139239 0.222194584 0.9999378 0.3458086        .base
      2.6   groupIX      100382          2 groupXVII      250065          9    149683 0.020408163 1.0000000 0.7750970        .base
      2.7   groupIX      100382          2    groupX      267302         10    166920 0.032582005 0.9815785 0.7180922        .base
      3      groupX      101821          3    groupI      182629          4     80808 0.005917160 1.0000000 0.8777310        .base
      3.1    groupX      101821          3   groupIV      194493          5     92672 0.047619048 0.3333333 0.6625206        .base
      3.2    groupX      101821          3  groupVII      208605          6    106784 0.009523810 1.0000000 0.8452520        .base
      3.3    groupX      101821          3   groupVI      212436          7    110615 0.062500000 1.0000000 0.6170751        .base
      3.4    groupX      101821          3  groupXXI      239621          8    137800 0.048813740 0.9985965 0.6585785        .base
      3.5    groupX      101821          3 groupXVII      250065          9    148244 0.020408163 1.0000000 0.7750970        .base
      3.6    groupX      101821          3    groupX      267302         10    165481 0.034539801 0.9957009 0.7101179        .base
      4      groupI      182629          4   groupIV      194493          5     11864          NA 0.0000000        NA        .base
      4.1    groupI      182629          4  groupVII      208605          6     25976          NA 0.0000000        NA        .base
      4.2    groupI      182629          4   groupVI      212436          7     29807 0.102564103 1.0000000 0.5218394        .base
      4.3    groupI      182629          4  groupXXI      239621          8     56992 0.042735043 1.0000000 0.6792776        .base
      4.4    groupI      182629          4 groupXVII      250065          9     67436 0.008264463 1.0000000 0.8557254        .base
      4.5    groupI      182629          4    groupX      267302         10     84673 0.018181818 1.0000000 0.7874065        .base
      5     groupIV      194493          5  groupVII      208605          6     14112          NA 0.0000000        NA        .base
      5.1   groupIV      194493          5   groupVI      212436          7     17943 0.092592006 0.9998079 0.5428037        .base
      5.2   groupIV      194493          5  groupXXI      239621          8     45128 0.000000000 0.0000000 1.0000000        .base
      5.3   groupIV      194493          5 groupXVII      250065          9     55572          NA 0.0000000        NA        .base
      5.4   groupIV      194493          5    groupX      267302         10     72809 0.135706586 0.9998987 0.4612638        .base
      6    groupVII      208605          6   groupVI      212436          7      3831 0.246370641 0.9999875 0.3208490        .base
      6.1  groupVII      208605          6  groupXXI      239621          8     31016 0.085714286 1.0000000 0.5581846        .base
      6.2  groupVII      208605          6 groupXVII      250065          9     41460 0.024023016 0.1549936 0.7565707        .base
      6.3  groupVII      208605          6    groupX      267302         10     58697 0.047619048 1.0000000 0.6625206        .base
      7     groupVI      212436          7  groupXXI      239621          8     27185 0.076875619 0.4295364 0.5792168        .base
      7.1   groupVI      212436          7 groupXVII      250065          9     37629 0.047496855 0.9987162 0.6629271        .base
      7.2   groupVI      212436          7    groupX      267302         10     54866 0.130421519 0.9999492 0.4701233        .base
      8    groupXXI      239621          8 groupXVII      250065          9     10444 0.003472222 0.1250000 0.9061856        .base
      8.1  groupXXI      239621          8    groupX      267302         10     27681 0.022222222 0.2000000 0.7655945        .base
      9   groupXVII      250065          9    groupX      267302         10     17237 0.032582005 0.9815785 0.7180922        .base
          sample.subfacet
      1             .base
      1.1           .base
      1.2           .base
      1.3           .base
      1.4           .base
      1.5           .base
      1.6           .base
      1.7           .base
      1.8           .base
      2             .base
      2.1           .base
      2.2           .base
      2.3           .base
      2.4           .base
      2.5           .base
      2.6           .base
      2.7           .base
      3             .base
      3.1           .base
      3.2           .base
      3.3           .base
      3.4           .base
      3.5           .base
      3.6           .base
      4             .base
      4.1           .base
      4.2           .base
      4.3           .base
      4.4           .base
      4.5           .base
      5             .base
      5.1           .base
      5.2           .base
      5.3           .base
      5.4           .base
      6             .base
      6.1           .base
      6.2           .base
      6.3           .base
      7             .base
      7.1           .base
      7.2           .base
      8             .base
      8.1           .base
      9             .base
      
      $LD$matrices
      $LD$matrices$.base
      $LD$matrices$.base$.base
      $LD$matrices$.base$.base$Dprime
             67921    100382 101821 182629    194493 208605    212436    239621    250065    267302
      67921     NA 0.9999389      1      1 1.0000000      1 0.9997940 0.1052797 1.0000000 0.4094221
      100382    NA        NA      1      1 1.0000000      1 0.9998711 0.9999378 1.0000000 0.9815785
      101821    NA        NA     NA      1 0.3333333      1 1.0000000 0.9985965 1.0000000 0.9957009
      182629    NA        NA     NA     NA 0.0000000      0 1.0000000 1.0000000 1.0000000 1.0000000
      194493    NA        NA     NA     NA        NA      0 0.9998079 0.0000000 0.0000000 0.9998987
      208605    NA        NA     NA     NA        NA     NA 0.9999875 1.0000000 0.1549936 1.0000000
      212436    NA        NA     NA     NA        NA     NA        NA 0.4295364 0.9987162 0.9999492
      239621    NA        NA     NA     NA        NA     NA        NA        NA 0.1250000 0.2000000
      250065    NA        NA     NA     NA        NA     NA        NA        NA        NA 0.9815785
      267302    NA        NA     NA     NA        NA     NA        NA        NA        NA        NA
      
      $LD$matrices$.base$.base$rsq
             67921   100382     101821     182629     194493     208605     212436     239621      250065     267302
      67921     NA 0.353994 0.06493506 0.04273504 0.15151515 0.10000000 0.19496611 0.01108381 0.030769231 0.10833000
      100382    NA       NA 0.02500000 0.02097902 0.07692308 0.02777778 0.25919241 0.22219458 0.020408163 0.03258200
      101821    NA       NA         NA 0.00591716 0.04761905 0.00952381 0.06250000 0.04881374 0.020408163 0.03453980
      182629    NA       NA         NA         NA         NA         NA 0.10256410 0.04273504 0.008264463 0.01818182
      194493    NA       NA         NA         NA         NA         NA 0.09259201 0.00000000          NA 0.13570659
      208605    NA       NA         NA         NA         NA         NA 0.24637064 0.08571429 0.024023016 0.04761905
      212436    NA       NA         NA         NA         NA         NA         NA 0.07687562 0.047496855 0.13042152
      239621    NA       NA         NA         NA         NA         NA         NA         NA 0.003472222 0.02222222
      250065    NA       NA         NA         NA         NA         NA         NA         NA          NA 0.03258200
      267302    NA       NA         NA         NA         NA         NA         NA         NA          NA         NA
      
      $LD$matrices$.base$.base$pval
             67921    100382    101821    182629    194493    208605    212436    239621    250065    267302
      67921     NA 0.2340669 0.6102987 0.6792776 0.4362749 0.5270893 0.3771826 0.8332312 0.7257210 0.5103644
      100382    NA        NA 0.7518296 0.7720590 0.5790997 0.7388827 0.3085740 0.3458086 0.7750970 0.7180922
      101821    NA        NA        NA 0.8777310 0.6625206 0.8452520 0.6170751 0.6585785 0.7750970 0.7101179
      182629    NA        NA        NA        NA        NA        NA 0.5218394 0.6792776 0.8557254 0.7874065
      194493    NA        NA        NA        NA        NA        NA 0.5428037 1.0000000        NA 0.4612638
      208605    NA        NA        NA        NA        NA        NA 0.3208490 0.5581846 0.7565707 0.6625206
      212436    NA        NA        NA        NA        NA        NA        NA 0.5792168 0.6629271 0.4701233
      239621    NA        NA        NA        NA        NA        NA        NA        NA 0.9061856 0.7655945
      250065    NA        NA        NA        NA        NA        NA        NA        NA        NA 0.7180922
      267302    NA        NA        NA        NA        NA        NA        NA        NA        NA        NA
      
      
      
      
      

