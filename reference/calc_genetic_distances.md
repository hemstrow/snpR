# Calculate genetic distances between individuals or groups of samples.

Calculates the genetic distances between either individuals or the
levels of any sample specific facets, broken apart by any requested snp
level facets. For details on methods, see details.

## Usage

``` r
calc_genetic_distances(
  x,
  facets = NULL,
  method = "Edwards",
  interpolate = "bernoulli"
)
```

## Arguments

- x:

  snpRdata object.

- facets:

  character or NULL, default NULL. Facets for which to calculate genetic
  distances, as described in
  [`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).
  If snp or base level facets are requested, distances will be between
  individuals. Otherwise, distances will be between the levels of the
  sample facets.

- method:

  character, default "Edwards". Name of the method to use. Options:

  - Edwards: Angular distance as described in Edwards

  1.  

  See details.

- interpolate:

  character, default "bernoulli". Missing data interpolation method,
  solely for individual/individual distances. Options detailed in
  documentation for
  [`format_snps`](https://hemstrow.github.io/snpR/reference/format_snps.md).

## Value

An overwrite-safe snpRdata object with genetic distance information (a
named, nested list containing distance measures named according to
facets and facet levels) added.

## Details

If a sample facet is requested, distances are calculated via code
adapted from code derived from
[`adegenet`](https://rdrr.io/pkg/adegenet/man/adegenet.package.html).
Please cite them alongside the tree-building and distance methods.
Available methods:

- Edwards: Angular distance as described in Edwards 1971.

- Nei: Nei's (1978) genetic distance measure.

Otherwise, genetic distances are calculated via the
[`dist`](https://rdrr.io/r/stats/dist.html) function using genotypes
formatted numerically (the "sn" option in
[`format_snps`](https://hemstrow.github.io/snpR/reference/format_snps.md)).
In all cases, missing genotypes are first interpolated according to the
method provided to the 'interpolate' argument.

## References

Edwards, A. W. F. (1971). Distances between populations on the basis of
gene frequencies. Biometrics, 873-881.

Nei, M. (1978). Estimation of average heterozygosity and genetic
distance from a small number of individuals. Genetics, 23, 341-369.

Jombart, T. (2008). adegenet: a R package for the multivariate analysis
of genetic markers Bioinformatics 24: 1403-1405.

## Author

William Hemstrom

## Examples

``` r
# by pop:
y <- calc_genetic_distances(stickSNPs, facets = "pop", method = "Edwards")
get.snpR.stats(y, "pop", "genetic_distance")
#> $pop
#> $pop$.base
#> $pop$.base$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.19567350                                            
#> OPL 0.14557964 0.17974134                                 
#> PAL 0.26273601 0.32767445 0.30325697                      
#> SMR 0.13011878 0.17413296 0.08655796 0.29362542           
#> UPD 0.10316254 0.19275330 0.13913664 0.26700420 0.12429385
#> 
#> 
#> 

# by chr and pop jointly
y <- calc_genetic_distances(stickSNPs, facets = "pop.chr", 
                            method = "Edwards")
get.snpR.stats(y, "pop.chr", "genetic_distance")
#> $chr.pop
#> $chr.pop$groupVI
#> $chr.pop$groupVI$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.25660019                                            
#> OPL 0.07924888 0.27433887                                 
#> PAL 0.27234262 0.39289021 0.29461931                      
#> SMR 0.07444861 0.22585586 0.07692285 0.25894041           
#> UPD 0.05132643 0.29256798 0.07277946 0.26101352 0.08212160
#> 
#> 
#> $chr.pop$groupXVIII
#> $chr.pop$groupXVIII$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.29143816                                            
#> OPL 0.07794949 0.22536849                                 
#> PAL 0.24295419 0.35226915 0.24914322                      
#> SMR 0.11471006 0.22023260 0.07335531 0.23435337           
#> UPD 0.15830277 0.23453873 0.13186523 0.17824321 0.13819295
#> 
#> 
#> $chr.pop$groupXV
#> $chr.pop$groupXV$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.03148000                                            
#> OPL 0.08838262 0.11703470                                 
#> PAL 0.20567053 0.18123586 0.28653542                      
#> SMR 0.08920877 0.10711398 0.05618791 0.28107228           
#> UPD 0.04349273 0.03724235 0.10239992 0.20802840 0.07730373
#> 
#> 
#> $chr.pop$groupXIII
#> $chr.pop$groupXIII$Edwards
#>           ASP       CLF       OPL       PAL       SMR
#> CLF 0.1320022                                        
#> OPL 0.1281510 0.1841410                              
#> PAL 0.2467409 0.2732171 0.2976303                    
#> SMR 0.1324830 0.1947162 0.0891242 0.3457030          
#> UPD 0.1285920 0.2007071 0.1080649 0.3287302 0.1008935
#> 
#> 
#> $chr.pop$groupV
#> $chr.pop$groupV$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.08047310                                            
#> OPL 0.17742141 0.14293958                                 
#> PAL 0.15885276 0.23008695 0.32878298                      
#> SMR 0.20720407 0.15183114 0.05765847 0.36223246           
#> UPD 0.08263034 0.06246718 0.20405474 0.19435913 0.21250749
#> 
#> 
#> $chr.pop$groupIX
#> $chr.pop$groupIX$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.20377523                                            
#> OPL 0.24116034 0.16352755                                 
#> PAL 0.32788690 0.39299334 0.40856673                      
#> SMR 0.17373485 0.14807071 0.10549804 0.35237666           
#> UPD 0.09069563 0.21100244 0.21569149 0.32600023 0.16503609
#> 
#> 
#> $chr.pop$groupXII
#> $chr.pop$groupXII$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.12255862                                            
#> OPL 0.15234040 0.09306658                                 
#> PAL 0.20851237 0.29957791 0.28689252                      
#> SMR 0.10404523 0.13912768 0.10890449 0.23105171           
#> UPD 0.08880017 0.10784281 0.10666556 0.26290491 0.07362240
#> 
#> 
#> $chr.pop$groupXVII
#> $chr.pop$groupXVII$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.20834023                                            
#> OPL 0.15762391 0.10263765                                 
#> PAL 0.17304121 0.21669826 0.12332919                      
#> SMR 0.16124832 0.10443489 0.04626631 0.12006859           
#> UPD 0.10740061 0.12888597 0.12896446 0.19335759 0.11841234
#> 
#> 
#> $chr.pop$groupIII
#> $chr.pop$groupIII$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.13153361                                            
#> OPL 0.06745803 0.14025208                                 
#> PAL 0.43961563 0.45689146 0.48673370                      
#> SMR 0.12653100 0.15243293 0.09622599 0.53831107           
#> UPD 0.07503234 0.16717941 0.11188195 0.43360338 0.15729501
#> 
#> 
#> $chr.pop$groupXIX
#> $chr.pop$groupXIX$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.21419025                                            
#> OPL 0.09448917 0.13096405                                 
#> PAL 0.22525640 0.26386578 0.24609677                      
#> SMR 0.13279287 0.09374108 0.07528015 0.19532331           
#> UPD 0.07426999 0.15886045 0.05221378 0.21176943 0.07854861
#> 
#> 
#> $chr.pop$groupXVI
#> $chr.pop$groupXVI$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.18557512                                            
#> OPL 0.08160442 0.11472564                                 
#> PAL 0.14004047 0.23725605 0.18891101                      
#> SMR 0.05278950 0.15225383 0.05097274 0.16168462           
#> UPD 0.07029203 0.15198356 0.08714485 0.10751838 0.06155318
#> 
#> 
#> $chr.pop$groupVII
#> $chr.pop$groupVII$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.14130459                                            
#> OPL 0.08951102 0.16991230                                 
#> PAL 0.30602675 0.37091117 0.31592877                      
#> SMR 0.08801032 0.14819459 0.09536588 0.32238416           
#> UPD 0.09868662 0.19315639 0.16554847 0.29734660 0.12329106
#> 
#> 
#> $chr.pop$groupXX
#> $chr.pop$groupXX$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.15620696                                            
#> OPL 0.17869189 0.11432505                                 
#> PAL 0.29425158 0.34754383 0.33935078                      
#> SMR 0.18360918 0.15825102 0.06443307 0.35870792           
#> UPD 0.09220821 0.15467737 0.16462493 0.26826140 0.18024172
#> 
#> 
#> $chr.pop$groupX
#> $chr.pop$groupX$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.17188736                                            
#> OPL 0.07988761 0.11771682                                 
#> PAL 0.23932255 0.30185710 0.23730310                      
#> SMR 0.09224339 0.08073638 0.05110458 0.25957026           
#> UPD 0.14436508 0.09627963 0.10733183 0.21444169 0.08516823
#> 
#> 
#> $chr.pop$groupVIII
#> $chr.pop$groupVIII$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.23818390                                            
#> OPL 0.06493169 0.20905381                                 
#> PAL 0.18820517 0.22655641 0.19705705                      
#> SMR 0.12023738 0.21583316 0.07284321 0.21492700           
#> UPD 0.04820242 0.23225900 0.06325434 0.15852107 0.11424897
#> 
#> 
#> $chr.pop$groupXI
#> $chr.pop$groupXI$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.39734021                                            
#> OPL 0.14876304 0.27606895                                 
#> PAL 0.13331858 0.42673887 0.23793114                      
#> SMR 0.04036849 0.37136179 0.14130976 0.12074936           
#> UPD 0.13551542 0.31216939 0.15290215 0.12199585 0.10100492
#> 
#> 
#> $chr.pop$groupIV
#> $chr.pop$groupIV$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.12514022                                            
#> OPL 0.21706730 0.25747130                                 
#> PAL 0.21987285 0.28423904 0.22446812                      
#> SMR 0.15657608 0.19871718 0.07465067 0.18235309           
#> UPD 0.09031265 0.19726779 0.15696764 0.18476523 0.11191674
#> 
#> 
#> $chr.pop$groupXXI
#> $chr.pop$groupXXI$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.13762901                                            
#> OPL 0.06982229 0.20695203                                 
#> PAL 0.02991273 0.16736897 0.03993042                      
#> SMR 0.08405276 0.22103950 0.01425140 0.05417410           
#> UPD 0.01614638 0.15368978 0.05369105 0.01376801 0.06792944
#> 
#> 
#> $chr.pop$groupI
#> $chr.pop$groupI$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.20042950                                            
#> OPL 0.08860246 0.13636896                                 
#> PAL 0.15314579 0.30182404 0.18552196                      
#> SMR 0.13763393 0.11864775 0.06522413 0.22000868           
#> UPD 0.14566557 0.15709630 0.10047934 0.22460525 0.08560579
#> 
#> 
#> $chr.pop$groupII
#> $chr.pop$groupII$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.22890137                                            
#> OPL 0.17246596 0.26983474                                 
#> PAL 0.29470893 0.30071240 0.27576095                      
#> SMR 0.09224878 0.20141061 0.11428097 0.22737207           
#> UPD 0.14117168 0.23817635 0.15215231 0.29154144 0.13471703
#> 
#> 
#> $chr.pop$groupXIV
#> $chr.pop$groupXIV$Edwards
#>            ASP        CLF        OPL        PAL        SMR
#> CLF 0.01814469                                            
#> OPL 0.11877785 0.13684865                                 
#> PAL 0.31420384 0.29648678 0.42890286                      
#> SMR 0.04656041 0.06469143 0.07231756 0.35943021           
#> UPD 0.08023063 0.09833949 0.03863940 0.39192320 0.03370171
#> 
#> 
#> 

if (FALSE) { # \dontrun{
# by pop and fam separately
y <- calc_genetic_distances(stickSNPs, facets = c("pop", "fam"), 
                            method = "Edwards")
get.snpR.stats(y, c("pop", "chr"), "genetic_distance")

# individuals across all snps + plot
y <- calc_genetic_distances(stickSNPs)
dat <- get.snpR.stats(y, stats = "genetic_distance")$.base$.base$Edwards
heatmap(as.matrix(dat))
} # }
```
