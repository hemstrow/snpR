# NeEstimator

    Code
      ne
    Output
      $pop
        facet pop LDNe_0.05 LDNe_0.02 LDNe_0.01 LDNe_lCIp_0.05 LDNe_lCIp_0.02 LDNe_lCIp_0.01 LDNe_uCIp_0.05 LDNe_uCIp_0.02 LDNe_uCIp_0.01
      1   pop ASP     662.2     890.6     935.5          485.5          639.7          705.8          885.8         1391.4         1688.7
        LDNe_lCIj_0.05 LDNe_lCIj_0.02 LDNe_lCIj_0.01 LDNe_uCIj_0.05 LDNe_uCIj_0.02 LDNe_uCIj_0.01
      1          257.4          319.5          334.2            Inf            Inf            Inf
      

# colony

    Code
      col$dyads[col$dyads$Probability > 0.5, -3]
    Output
            Sample1    Sample2    type
      2  ASP.r2.143 ASP.r2.188 HalfSib
      3  ASP.r2.143 ASP.r2.195 HalfSib
      4  ASP.r2.143 ASP.r2.227 HalfSib
      8  ASP.r2.182 ASP.r2.232 HalfSib
      10 ASP.r2.186 ASP.r2.215 HalfSib
      12 ASP.r2.187 ASP.r2.216 HalfSib
      14 ASP.r2.188 ASP.r2.195 HalfSib
      15 ASP.r2.188 ASP.r2.232 HalfSib

---

    Code
      col$clusters[col$clusters$ClusterProbability > 0.5, 3]
    Output
      [38;5;246m# A tibble: 40 x 1[39m
         OffspringID
         [3m[38;5;246m<chr>[39m[23m      
      [38;5;250m 1[39m ASP.r2.140 
      [38;5;250m 2[39m ASP.r2.141 
      [38;5;250m 3[39m ASP.r2.142 
      [38;5;250m 4[39m ASP.r2.144 
      [38;5;250m 5[39m ASP.r2.145 
      [38;5;250m 6[39m ASP.r2.146 
      [38;5;250m 7[39m ASP.r2.147 
      [38;5;250m 8[39m ASP.r2.148 
      [38;5;250m 9[39m ASP.r2.149 
      [38;5;250m10[39m ASP.r2.150 
      [38;5;246m# ... with 30 more rows[39m

