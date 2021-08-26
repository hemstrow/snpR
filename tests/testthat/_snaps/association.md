# correct gmmat

    Code
      asgmmat
    Output
      $single
         facet subfacet     group position gmmat_n_phenotype gmmat_af_phenotype gmmat_score_phenotype gmmat_var_phenotype
      1  .base    .base  groupXIX    67921                 8          0.3125000            -0.7080890            1.214100
      2  .base    .base   groupIX   100382                 9          0.1666670             0.3038430            0.707075
      3  .base    .base    groupX   101821                10          0.1000000            -1.2219700            0.571736
      4  .base    .base    groupI   182629                 7          0.0714286             1.1361300            0.480942
      5  .base    .base   groupIV   194493                 8          0.2500000            -0.5835270            1.448930
      6  .base    .base  groupVII   208605                 8          0.1250000             0.2930060            0.520488
      7  .base    .base   groupVI   212436                 9          0.3333330             1.9072500            1.313800
      8  .base    .base  groupXXI   239621                 9          0.3888890            -0.2357000            2.005750
      9  .base    .base groupXVII   250065                 8          0.1250000             0.5114610            0.766673
      10 .base    .base    groupX   267302                 9          0.2777780             0.0388775            0.949583
         gmmat_pval_phenotype
      1             0.5204650
      2             0.7178450
      3             0.1060780
      4             0.1013690
      5             0.6278380
      6             0.6846430
      7             0.0961199
      8             0.8678220
      9             0.5591350
      10            0.9681760
      

# correct odds

    Code
      asodds
    Output
      $single
         subfacet facet     group position log_odds_ratio_cat_phenotype se_cat_phenotype
      1     .base .base  groupXIX    67921                   -0.8450980      0.107799900
      2     .base .base   groupIX   100382                   -0.4771213      0.124938737
      3     .base .base    groupX   101821                   -0.1962946      0.174518861
      4     .base .base    groupI   182629                          Inf              Inf
      5     .base .base   groupIV   194493                    0.0000000      0.062469368
      6     .base .base  groupVII   208605                          Inf              Inf
      7     .base .base   groupVI   212436                    0.3010300      0.017381053
      8     .base .base  groupXXI   239621                   -0.3679768     -0.005232717
      9     .base .base groupXVII   250065                          Inf              Inf
      10    .base .base    groupX   267302                   -0.1760913      0.038582977
      

# correct chisq

    Code
      aschi
    Output
      $single
         subfacet facet     group position chi_stat_cat_phenotype chi_p_cat_phenotype associated_allele_cat_phenotype
      1     .base .base  groupXIX    67921             2.61818182           0.1056454                         minor_A
      2     .base .base   groupIX   100382             0.72000000           0.3961439                         minor_A
      3     .base .base    groupX   101821             0.09259259           0.7609067                         minor_A
      4     .base .base    groupI   182629             1.43589744           0.2308044                         major_A
      5     .base .base   groupIV   194493             0.00000000           1.0000000                            NA_A
      6     .base .base  groupVII   208605             1.37142857           0.2415666                         major_A
      7     .base .base   groupVI   212436             0.45000000           0.5023350                         major_A
      8     .base .base  groupXXI   239621             0.74805195           0.3870937                         minor_A
      9     .base .base groupXVII   250065             1.37142857           0.2415666                         major_A
      10    .base .base    groupX   267302             0.13846154           0.7098153                         minor_A
      

