# Summarize possible snpRdata object facet options

List either all of the possible SNP and sample facets (if called with no
facets) or all of the categories for each requested facet.

## Usage

``` r
summarize_facets(x, facets = NULL)
```

## Arguments

- x:

  snpRdata object

- facets:

  character. Categorical metadata variables by which to break up
  analysis. See
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md)
  for more details. If NULL, the possible SNP and sample facets will be
  listed. If facets are instead provided, the *categories* for each
  facet will instead be listed.

## Value

A named list containing either the possible SNP and sample facets or the
categories for all of the requested facets.

## Details

If called with sample facets, returns a table of counts of each category
otherwise returns a vector of categories.

## Author

William Hemstrom

## Examples

``` r
# list available facets
summarize_facets(stickSNPs)
#> Returning list of facets. For more details, try asking for information about a specific facet (such as 'pop' or 'pop.chr')!
#> $SNP
#> [1] "chr"      "position" ".snp.id" 
#> 
#> $sample
#> [1] "pop"        "fam"        ".sample.id"
#> 

# return details for a few facets
summarize_facets(stickSNPs, c("pop", "chr.pop", "fam.pop"))
#> $pop
#> 
#> ASP CLF OPL PAL SMR UPD 
#>  11   6  27  17  24  15 
#> 
#> $chr.pop
#>   [1] "groupVI.ASP"    "groupVI.CLF"    "groupVI.OPL"    "groupVI.PAL"   
#>   [5] "groupVI.SMR"    "groupVI.UPD"    "groupXVIII.ASP" "groupXVIII.CLF"
#>   [9] "groupXVIII.OPL" "groupXVIII.PAL" "groupXVIII.SMR" "groupXVIII.UPD"
#>  [13] "groupXV.ASP"    "groupXV.CLF"    "groupXV.OPL"    "groupXV.PAL"   
#>  [17] "groupXV.SMR"    "groupXV.UPD"    "groupXIII.ASP"  "groupXIII.CLF" 
#>  [21] "groupXIII.OPL"  "groupXIII.PAL"  "groupXIII.SMR"  "groupXIII.UPD" 
#>  [25] "groupV.ASP"     "groupV.CLF"     "groupV.OPL"     "groupV.PAL"    
#>  [29] "groupV.SMR"     "groupV.UPD"     "groupIX.ASP"    "groupIX.CLF"   
#>  [33] "groupIX.OPL"    "groupIX.PAL"    "groupIX.SMR"    "groupIX.UPD"   
#>  [37] "groupXII.ASP"   "groupXII.CLF"   "groupXII.OPL"   "groupXII.PAL"  
#>  [41] "groupXII.SMR"   "groupXII.UPD"   "groupXVII.ASP"  "groupXVII.CLF" 
#>  [45] "groupXVII.OPL"  "groupXVII.PAL"  "groupXVII.SMR"  "groupXVII.UPD" 
#>  [49] "groupIII.ASP"   "groupIII.CLF"   "groupIII.OPL"   "groupIII.PAL"  
#>  [53] "groupIII.SMR"   "groupIII.UPD"   "groupXIX.ASP"   "groupXIX.CLF"  
#>  [57] "groupXIX.OPL"   "groupXIX.PAL"   "groupXIX.SMR"   "groupXIX.UPD"  
#>  [61] "groupXVI.ASP"   "groupXVI.CLF"   "groupXVI.OPL"   "groupXVI.PAL"  
#>  [65] "groupXVI.SMR"   "groupXVI.UPD"   "groupVII.ASP"   "groupVII.CLF"  
#>  [69] "groupVII.OPL"   "groupVII.PAL"   "groupVII.SMR"   "groupVII.UPD"  
#>  [73] "groupXX.ASP"    "groupXX.CLF"    "groupXX.OPL"    "groupXX.PAL"   
#>  [77] "groupXX.SMR"    "groupXX.UPD"    "groupX.ASP"     "groupX.CLF"    
#>  [81] "groupX.OPL"     "groupX.PAL"     "groupX.SMR"     "groupX.UPD"    
#>  [85] "groupVIII.ASP"  "groupVIII.CLF"  "groupVIII.OPL"  "groupVIII.PAL" 
#>  [89] "groupVIII.SMR"  "groupVIII.UPD"  "groupXI.ASP"    "groupXI.CLF"   
#>  [93] "groupXI.OPL"    "groupXI.PAL"    "groupXI.SMR"    "groupXI.UPD"   
#>  [97] "groupIV.ASP"    "groupIV.CLF"    "groupIV.OPL"    "groupIV.PAL"   
#> [101] "groupIV.SMR"    "groupIV.UPD"    "groupXXI.ASP"   "groupXXI.CLF"  
#> [105] "groupXXI.OPL"   "groupXXI.PAL"   "groupXXI.SMR"   "groupXXI.UPD"  
#> [109] "groupI.ASP"     "groupI.CLF"     "groupI.OPL"     "groupI.PAL"    
#> [113] "groupI.SMR"     "groupI.UPD"     "groupII.ASP"    "groupII.CLF"   
#> [117] "groupII.OPL"    "groupII.PAL"    "groupII.SMR"    "groupII.UPD"   
#> [121] "groupXIV.ASP"   "groupXIV.CLF"   "groupXIV.OPL"   "groupXIV.PAL"  
#> [125] "groupXIV.SMR"   "groupXIV.UPD"  
#> 
#> $fam.pop
#> 
#> A.ASP A.CLF A.OPL A.PAL A.SMR A.UPD B.ASP B.CLF B.OPL B.PAL B.SMR B.UPD C.ASP 
#>     3     2     6     5     6     3     3     1     7     4     6     4     3 
#> C.CLF C.OPL C.PAL C.SMR C.UPD D.ASP D.CLF D.OPL D.PAL D.SMR D.UPD 
#>     1     7     4     6     4     2     2     7     4     6     4 
#> 
```
