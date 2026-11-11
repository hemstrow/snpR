# Subset snpRdata objects

Subsets snpRdata objects by specific snps, samples, facets, subfacets,
etc. Statistics will need to be recalculated after subsetting to avoid
confusion. The bracket operators can be alternatively, following the
same syntax.

## Usage

``` r
# S4 method for class 'snpRdata,ANY,ANY,ANY'
x[i, j, ..., drop = FALSE]

subset_snpR_data(
  x,
  .snps = 1:nsnps(x),
  .samps = 1:nsamps(x),
  ...,
  .facets = NULL,
  .subfacets = NULL,
  .snp.facets = NULL,
  .snp.subfacets = NULL
)
```

## Arguments

- x:

  snpRdata object.

- i:

  numeric, default `1:nrow(x)`. Row numbers corresponding to SNPs to
  keep. Negative subscripts allowed.

- j:

  numeric, default `1:nrow(x)`. Row numbers corresponding to SNPs to
  keep. Negative subscripts allowed.

- ...:

  Facet subsetting may be specified by providing the facet as an
  argument and then providing the levels to keep or remove. Setting pop
  = "ASP", for example, would keep only samples in the "ASP" level of
  the pop facet, and setting pop.fam = "ASP.A" would keep only samples
  the ASP pop and the A fam. Negative subscripts are allowed: pop =
  -c("ASP", "PAL") would remove samples in the ASP or PAL pops.
  Subsetting by multiple facets is supported, although negative and
  positive subscripts cannot be mixed across sample or SNP facets. They
  may be mixed between the two.

- drop:

  logical, default FALSE. Deprecated.

- .snps:

  numeric, default `1:nrow(x)`. Row numbers corresponding to SNPs to
  keep. Negative subscripts allowed.

- .samps:

  numeric, default `1:ncol(x)`. Column numbers corresponding to samples
  to keep. Negative subscripts allowed.

- .facets:

  Character, default NULL. Sample facet to select by. Alternative to
  direct specification as in `pop = "ASP"`.

- .subfacets:

  Character, default NULL. Sample facet level matching a level of the
  provided `.facets`. Alternative to direct specification as in
  `pop = "ASP"`.

- .snp.facets:

  Character, default NULL. SNP facet to select by. Alternative to direct
  specification as in `chr = "groupIX"`.

- .snp.subfacets:

  Character, default NULL. SNP facet level matching a level of the
  provided `.snp.facets`. Alternative to direct specification as in
  `chr = "groupIX"`.

## Details

Sample and snp facets to subset over can be provided. Filtering by facet
categories is performed by naming the facet as an argument name, then
providing the levels to keep or remove to that argument. See examples.
Facets designated as described in
[`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).

## Functions

- `x[i`: extraction operator

- `subset_snpR_data()`: subset_snpR_data

## Examples

``` r
# Keep only individuals in the ASP and PAL populations 
# and on the LGIX or LGIV chromosome.
subset_snpR_data(stickSNPs, pop = c("ASP", "PAL"), 
                 chr = c("groupIX", "groupIV"))
#>    ASP.r2.141 ASP.r2.144 ASP.r2.146 ASP.r2.147 ASP.r2.196 ASP.r2.216 ASP.r2.221
#> 1          TT         TT         TT         TT         GT         TT         TT
#> 2          CT         NN         TT         CC         CT         NN         CC
#> 3          AA         AA         AA         NN         AA         AA         AA
#> 4          AA         AA         AA         AA         AA         AA         AA
#> 5          TT         NN         NN         NN         CT         CT         TT
#> 6          AA         AA         AA         AA         AA         AG         AA
#> 7          AG         AA         AG         AG         AA         AA         AG
#> 8          GG         GG         GG         GG         GG         GG         GG
#> 9          CT         CC         CT         CC         CC         CT         CT
#> 10         AA         AA         NN         AA         AA         AA         AA
#> 11         AC         AC         AA         AC         AC         NN         AC
#> 12         GG         GG         GG         GG         GG         GG         GG
#> 13         AT         TT         AT         TT         TT         AT         TT
#> 14         CT         CC         CT         CT         TT         CT         CC
#> 15         AG         NN         AG         NN         GG         AG         AA
#> 16         CC         NN         CC         NN         CC         CC         CC
#> 17         AT         NN         AT         AT         TT         AT         TT
#>    ASP.r2.227 ASP.r2.230 ASP.r2.232 ASP.r2.235 PAL.r2.203 PAL.r2.237 PAL.r2.238
#> 1          TT         TT         TT         TT         TT         TT         TT
#> 2          TT         CT         CC         CT         CC         CC         CC
#> 3          AA         AA         AA         AG         AA         AA         AA
#> 4          AA         NN         AC         AA         AA         AA         AA
#> 5          TT         TT         TT         TT         TT         TT         TT
#> 6          AG         AA         AG         NN         AA         AA         AA
#> 7          AG         GG         AG         AG         AA         AA         AA
#> 8          GG         GG         GG         GG         GG         GG         GG
#> 9          CC         CC         CC         CC         CC         CC         CT
#> 10         AA         AA         AA         AA         AA         AA         AA
#> 11         AC         AA         AA         AC         AA         AA         AA
#> 12         GG         NN         GG         GG         GG         GG         GG
#> 13         TT         TT         AT         TT         TT         TT         TT
#> 14         CT         TT         TT         CC         TT         TT         NN
#> 15         AG         AG         AG         NN         GG         NN         GG
#> 16         CC         AC         CC         NN         CC         CC         CC
#> 17         TT         TT         AA         AA         NN         NN         AT
#>    PAL.r2.239 PAL.r2.240 PAL.r2.241 PAL.r2.242 PAL.r2.249 PAL.r2.253 PAL.r2.275
#> 1          TT         TT         TT         TT         TT         TT         TT
#> 2          CC         CC         CC         CC         CC         CC         CC
#> 3          AA         AA         AA         AA         AA         NN         AA
#> 4          AA         AA         AA         AA         AA         AA         AA
#> 5          TT         TT         TT         CT         TT         CT         TT
#> 6          AA         AA         AA         AA         AA         AA         AA
#> 7          AA         AA         AA         NN         AA         AA         NN
#> 8          GG         GG         GG         GG         GG         GG         GG
#> 9          CC         CC         CC         CC         CC         CT         CT
#> 10         AA         AA         AA         AA         AA         AA         AA
#> 11         AA         AA         NN         AA         AA         AA         AA
#> 12         GG         GG         GG         GG         GG         GG         GG
#> 13         TT         TT         TT         TT         TT         TT         TT
#> 14         TT         TT         TT         TT         TT         TT         TT
#> 15         GG         GG         GG         GG         GG         GG         GG
#> 16         CC         CC         CC         CC         CC         CC         CC
#> 17         NN         TT         AA         NN         AA         NN         AT
#>    PAL.r2.276 PAL.r2.278 PAL.r2.281 PAL.r2.285 PAL.r2.287 PAL.r2.290 PAL.r2.291
#> 1          TT         TT         TT         TT         TT         TT         TT
#> 2          CC         CC         CC         CC         CC         CC         CC
#> 3          AA         AA         AA         AA         AA         AA         AA
#> 4          AA         AA         AA         AA         AA         AA         AA
#> 5          TT         TT         NN         NN         TT         CT         TT
#> 6          AA         AA         AA         AA         AA         NN         AA
#> 7          NN         AA         NN         AA         AA         AA         AA
#> 8          GG         GG         GG         GG         GG         NN         GG
#> 9          CC         CT         CT         TT         CT         CC         CT
#> 10         AA         AA         AA         AA         AA         AA         AA
#> 11         AA         AA         AA         AA         AA         AA         AA
#> 12         GG         GG         GG         GG         GG         GG         GG
#> 13         TT         TT         TT         TT         NN         TT         TT
#> 14         TT         TT         TT         TT         TT         TT         TT
#> 15         GG         GG         GG         GG         GG         NN         NN
#> 16         CC         CC         CC         CC         CC         CC         CC
#> 17         AT         AT         AA         AT         AA         NN         AA
                 
# keep individuals/SNPs in the first 10 columns/rows.
subset_snpR_data(stickSNPs, 1:10, 1:10)
#>    ASP.r2.141 ASP.r2.144 ASP.r2.146 ASP.r2.147 ASP.r2.196 ASP.r2.216 ASP.r2.221
#> 1          AG         GG         AA         AA         NN         AG         NN
#> 2          GG         GG         GG         GG         GT         GG         GG
#> 3          CC         CC         CC         CC         CC         CC         CT
#> 4          AA         AA         AA         AA         AG         AA         AA
#> 5          GG         NN         NN         AA         GG         AG         AG
#> 6          AC         AC         AA         AA         AA         AC         AA
#> 7          TT         TT         TT         TT         GT         TT         TT
#> 8          GG         GG         AG         AG         AG         GG         GG
#> 9          CC         AC         CC         AC         CC         CC         CC
#> 10         AA         AA         AG         AA         AA         AA         AA
#>    ASP.r2.227 ASP.r2.230 ASP.r2.232
#> 1          GG         AA         GG
#> 2          NN         GG         GG
#> 3          CC         CC         CT
#> 4          AA         AA         AG
#> 5          AG         NN         GG
#> 6          AC         AA         AC
#> 7          TT         TT         TT
#> 8          AG         GG         AG
#> 9          CC         AC         CC
#> 10         AA         AG         AA

# negative subscripts: remove individuals in family B
subset_snpR_data(stickSNPs, fam = -"B")
#>     ASP.r2.141 ASP.r2.146 ASP.r2.147 ASP.r2.196 ASP.r2.221 ASP.r2.227
#> 1           AG         AA         AA         NN         NN         GG
#> 2           GG         GG         GG         GT         GG         NN
#> 3           CC         CC         CC         CC         CT         CC
#> 4           AA         AA         AA         AG         AA         AA
#> 5           GG         NN         AA         GG         AG         AG
#> 6           AC         AA         AA         AA         AA         AC
#> 7           TT         TT         TT         GT         TT         TT
#> 8           GG         AG         AG         AG         GG         AG
#> 9           CC         CC         AC         CC         CC         CC
#> 10          AA         AG         AA         AA         AA         AA
#> 11          NN         TT         CT         TT         CC         CC
#> 12          AG         AG         NN         AA         AG         AA
#> 13          TT         CT         TT         CT         TT         TT
#> 14          GG         GT         NN         GG         GG         GG
#> 15          AA         AG         AA         AA         AA         AA
#> 16          NN         GT         TT         GG         GT         GG
#> 17          AA         AA         AA         AA         AC         AC
#> 18          AC         CC         CC         AC         CC         AC
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         GG         GG         NN         GG         GG
#> 21          NN         CC         CC         AC         CC         CC
#> 22          CT         TT         CC         CT         CC         TT
#> 23          GG         GG         GG         GG         GG         AG
#> 24          AA         AA         NN         AA         AA         AA
#> 25          NN         AG         GG         GG         AG         GG
#> 26          CT         CT         CT         NN         NN         TT
#> 27          AA         AA         AA         AA         AA         AA
#> 28          TT         NN         NN         CT         TT         TT
#> 29          AA         AA         NN         AA         AA         NN
#> 30          GG         CG         GG         GG         GG         GG
#> 31          CG         CC         CC         CC         CG         CC
#> 32          AA         AA         NN         AA         AA         AC
#> 33          CC         CC         CC         CT         CT         CC
#> 34          CG         CG         CG         CC         CC         CG
#> 35          AA         AA         AC         AA         AC         AA
#> 36          AA         AA         AA         AA         AA         AG
#> 37          CG         GG         GG         GG         NN         GG
#> 38          GG         GG         GT         GT         GG         GT
#> 39          TT         TT         NN         AT         TT         TT
#> 40          CG         CG         NN         CG         CG         GG
#> 41          AG         AG         AG         AA         AG         AG
#> 42          AG         NN         NN         GG         AG         AG
#> 43          AG         AA         AA         AA         AA         AA
#> 44          AT         TT         TT         NN         TT         NN
#> 45          GG         GG         GG         GG         GG         GG
#> 46          NN         AA         AA         NN         NN         AA
#> 47          CC         CC         CG         CG         CG         CG
#> 48          CT         TT         NN         TT         CT         CC
#> 49          GG         GG         GG         NN         GG         GG
#> 50          NN         CT         CC         CC         CC         CC
#> 51          TT         CT         CT         TT         CT         NN
#> 52          AA         AG         AA         AG         AA         AG
#> 53          AA         CC         NN         AA         AA         AA
#> 54          GG         GG         GG         GG         GG         GG
#> 55          GG         AG         NN         AG         AG         AG
#> 56          AG         AG         AG         AA         AG         AA
#> 57          GG         GG         NN         GG         GG         AG
#> 58          TT         CC         NN         TT         NN         CT
#> 59          AA         GG         AG         GG         GG         AG
#> 60          AA         AC         NN         AC         AC         NN
#> 61          TT         TT         NN         TT         TT         NN
#> 62          CG         GG         NN         GG         GG         GG
#> 63          AC         NN         NN         AA         AC         AC
#> 64          TT         AT         TT         AT         AA         AT
#> 65          CT         CT         CC         CC         CT         CC
#> 66          AC         NN         NN         AC         CC         AC
#> 67          GG         GG         NN         GG         GG         GG
#> 68          AA         AT         AA         AA         AT         TT
#> 69          AA         NN         AA         AA         AA         AA
#> 70          AC         AA         AC         AC         AC         AC
#> 71          AA         AA         AT         AA         AA         AA
#> 72          GG         GG         GG         AG         AG         AG
#> 73          TT         TT         TT         GT         GT         NN
#> 74          GT         GT         GG         GG         GT         GG
#> 75          AC         NN         AC         CC         AC         AC
#> 76          CC         CT         CT         NN         TT         NN
#> 77          AT         AT         AT         AA         AT         NN
#> 78          GG         GG         GG         GG         GG         GG
#> 79          TT         CT         TT         CT         CT         TT
#> 80          CT         TT         TT         TT         TT         TT
#> 81          AA         AT         AT         AT         AT         AA
#> 82          CT         CT         TT         TT         NN         TT
#> 83          AT         AT         TT         TT         TT         TT
#> 84          CC         CT         CT         CC         CT         CC
#> 85          AG         GG         AA         AG         AA         AA
#> 86          AG         AG         AG         AA         AG         NN
#> 87          CT         CT         CT         TT         CC         CT
#> 88          AG         AG         NN         GG         AA         AG
#> 89          NN         CG         CG         CC         CG         NN
#> 90          CT         CC         CC         CC         CC         CC
#> 91          CC         NN         CT         CC         CC         CT
#> 92          AA         GG         NN         GG         AG         AG
#> 93          TT         TT         TT         CT         TT         TT
#> 94          CC         CC         CC         CC         CC         CC
#> 95          GG         GG         NN         GT         GT         GT
#> 96          CC         CC         NN         CC         CC         CC
#> 97          GG         GG         GG         GG         GG         GG
#> 98          GG         GG         NN         GG         GG         GT
#> 99          AA         AA         AA         AG         AG         AG
#> 100         AT         AT         AT         TT         TT         TT
#>     ASP.r2.230 ASP.r2.235 CLF.r2.040 CLF.r2.048 CLF.r2.064 CLF.r2.068
#> 1           AA         AG         NN         AG         AA         AG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AA         AG         AA         AG         NN         AA
#> 5           NN         NN         GG         NN         NN         GG
#> 6           AA         NN         AA         AC         AA         AC
#> 7           TT         TT         GG         TT         TT         TT
#> 8           GG         GG         GG         AG         AA         AG
#> 9           AC         AA         CC         AC         AC         CC
#> 10          AG         AA         AA         AA         AA         AG
#> 11          CT         CT         CT         CT         CT         TT
#> 12          AG         AA         AA         AA         AA         AA
#> 13          NN         TT         CT         CT         CT         TT
#> 14          GT         GT         GG         GG         GT         GT
#> 15          AG         AG         AA         AA         AA         AA
#> 16          GG         GG         GT         GT         NN         GG
#> 17          AA         AA         AA         AA         AA         AA
#> 18          AC         AC         CC         AC         NN         AC
#> 19          TT         TT         GT         TT         GT         TT
#> 20          GG         GG         GG         GG         GG         GG
#> 21          AC         CC         NN         AC         NN         AC
#> 22          CT         CT         CC         CC         TT         CT
#> 23          GG         NN         GG         GG         GG         GG
#> 24          AA         AG         AA         AA         AA         AA
#> 25          AG         GG         GG         GG         NN         GG
#> 26          CT         CT         TT         TT         TT         TT
#> 27          NN         AA         AA         AC         AA         AA
#> 28          TT         TT         CT         TT         TT         TT
#> 29          AA         AA         AA         NN         AA         AC
#> 30          GG         CG         GG         GG         GG         CG
#> 31          CG         CG         CC         CG         CC         CC
#> 32          AC         AA         AC         AC         AA         AA
#> 33          CC         CT         CT         CC         NN         CT
#> 34          CG         NN         CC         CG         CG         GG
#> 35          AA         AA         AA         NN         AA         AA
#> 36          AA         NN         AG         AA         AG         AA
#> 37          CC         CG         GG         CG         NN         GG
#> 38          NN         TT         GG         GG         GG         GG
#> 39          TT         AT         TT         NN         AT         NN
#> 40          CC         CG         CG         GG         CG         GG
#> 41          GG         AG         AG         GG         AG         NN
#> 42          AG         GG         NN         GG         GG         AG
#> 43          AA         AA         AG         NN         AA         AG
#> 44          TT         TT         TT         TT         TT         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          CC         AC         AC         CC         NN         CC
#> 47          CC         CC         CC         CC         CC         CC
#> 48          NN         CT         CT         CT         CT         CT
#> 49          GG         GG         GG         GG         GG         CG
#> 50          CT         CC         CC         CC         CC         CT
#> 51          CT         TT         TT         NN         TT         TT
#> 52          AG         AG         GG         AG         AA         AA
#> 53          AC         AA         AC         NN         NN         AC
#> 54          GG         GG         GG         GG         GT         GG
#> 55          AG         NN         GG         GG         AG         AG
#> 56          AG         NN         AA         AA         AA         AA
#> 57          NN         GG         NN         GG         GG         GG
#> 58          CC         CT         CC         CC         CC         CC
#> 59          GG         AG         AG         AG         AG         GG
#> 60          NN         NN         AC         AC         CC         CC
#> 61          TT         TT         TT         NN         CT         TT
#> 62          CG         GG         GG         GG         GG         GG
#> 63          CC         AA         CC         AC         AA         NN
#> 64          TT         TT         TT         TT         TT         TT
#> 65          CC         CC         TT         CC         CT         CT
#> 66          AA         AC         AC         AC         AA         AC
#> 67          NN         GG         GG         GG         GG         GG
#> 68          AA         AA         AT         AA         AT         AA
#> 69          AA         AA         NN         AA         AA         AG
#> 70          AA         AC         CC         AA         NN         AC
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AA         AA         AG         AA         GG         GG
#> 73          GT         GT         GT         TT         TT         GT
#> 74          GG         NN         GG         GT         GG         NN
#> 75          AA         NN         AA         AC         CC         AC
#> 76          NN         NN         CT         NN         NN         CC
#> 77          AT         AA         AA         AA         AA         AA
#> 78          NN         GG         TT         GG         GG         GT
#> 79          TT         TT         TT         TT         TT         TT
#> 80          CC         TT         CT         CT         TT         TT
#> 81          AT         AA         TT         NN         AT         AT
#> 82          TT         TT         TT         NN         TT         NN
#> 83          TT         TT         AT         AT         AT         AT
#> 84          CT         TT         TT         CT         CC         CT
#> 85          AG         AG         AG         AA         AG         AG
#> 86          AG         NN         AA         AA         AA         AA
#> 87          TT         CC         TT         CT         TT         TT
#> 88          AG         NN         AA         AA         AA         AG
#> 89          CC         CC         CC         NN         CC         NN
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CC         CC         CC         NN         CC         NN
#> 92          NN         AG         AG         NN         GG         GG
#> 93          TT         CT         CT         TT         TT         CT
#> 94          CC         NN         CC         CC         CC         NN
#> 95          GG         GG         GT         TT         GG         GT
#> 96          AC         NN         CC         CC         AC         CC
#> 97          GG         GG         GG         GG         GG         AG
#> 98          NN         TT         GG         GG         NN         GT
#> 99          NN         AG         AG         AA         AA         GG
#> 100         TT         AA         AT         AT         AA         NN
#>     CLF.r2.069 OPL.r1.080 OPL.r1.083 OPL.r1.086 OPL.r1.088 OPL.r1.109
#> 1           AA         AG         AG         AG         NN         NN
#> 2           GG         NN         GT         GG         GG         GT
#> 3           CT         CT         CT         CC         CT         CC
#> 4           AA         AG         NN         AG         NN         NN
#> 5           GG         NN         NN         AG         NN         GG
#> 6           AA         AA         AC         AA         AA         NN
#> 7           TT         NN         GT         TT         TT         NN
#> 8           GG         GG         GG         GG         GG         AG
#> 9           CC         AC         CC         CC         NN         CC
#> 10          AA         AA         AA         AA         AA         AA
#> 11          CC         NN         NN         CC         CC         CC
#> 12          AA         NN         NN         AA         AG         AG
#> 13          CC         TT         NN         TT         TT         CT
#> 14          GG         GG         GG         GG         GG         GG
#> 15          AA         AG         AA         AA         AA         NN
#> 16          GT         GT         NN         NN         GT         NN
#> 17          AA         AA         AC         NN         NN         AA
#> 18          AC         CC         CC         CC         CC         CC
#> 19          GT         TT         TT         GT         TT         NN
#> 20          AG         GG         NN         GG         GG         GG
#> 21          CC         AC         AA         NN         AC         NN
#> 22          CT         CT         CT         NN         CC         NN
#> 23          GG         GG         AG         GG         GG         GG
#> 24          AA         AA         AA         AG         AA         AA
#> 25          GG         NN         AG         GG         GG         AG
#> 26          TT         TT         TT         NN         TT         TT
#> 27          AA         AC         AA         AC         AA         NN
#> 28          TT         CT         CT         NN         NN         CT
#> 29          AA         AA         NN         AA         AC         NN
#> 30          NN         NN         GG         NN         CG         NN
#> 31          CC         CC         CG         CC         CC         CC
#> 32          AA         AA         AC         NN         AA         AC
#> 33          NN         NN         CC         CT         NN         TT
#> 34          CC         NN         CG         NN         NN         NN
#> 35          AA         AA         AA         NN         AA         AA
#> 36          AA         AG         AA         AA         NN         AG
#> 37          GG         CG         GG         CG         CG         CG
#> 38          GG         GG         NN         NN         GT         GT
#> 39          AT         NN         NN         TT         NN         AT
#> 40          GG         NN         GG         GG         NN         CG
#> 41          AG         NN         NN         NN         NN         AG
#> 42          AA         AG         AG         AG         GG         AG
#> 43          AA         AA         AA         NN         AG         AG
#> 44          TT         TT         TT         NN         TT         TT
#> 45          GG         GG         NN         GG         GG         GG
#> 46          AC         AC         NN         NN         NN         AA
#> 47          CC         CG         CC         CG         CG         CG
#> 48          CC         CC         NN         CT         NN         NN
#> 49          GG         GG         GG         GG         GG         GG
#> 50          NN         CC         CC         CT         NN         CC
#> 51          TT         CC         NN         TT         NN         CT
#> 52          AG         AG         AA         NN         AG         AA
#> 53          AC         AA         AA         AC         NN         AC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          AG         NN         NN         NN         GG         AG
#> 56          AA         AA         NN         NN         AG         AG
#> 57          GG         NN         NN         NN         GG         GG
#> 58          CC         CT         CC         CC         CT         CT
#> 59          GG         AG         AG         AG         AA         GG
#> 60          AC         AC         AC         NN         AC         AC
#> 61          NN         NN         NN         TT         NN         TT
#> 62          GG         GG         CG         GG         GG         GG
#> 63          AA         AC         NN         AA         AA         AA
#> 64          AA         NN         AT         NN         TT         AT
#> 65          CC         NN         TT         NN         TT         CT
#> 66          AA         AC         AC         AC         AC         NN
#> 67          GG         NN         GG         GG         GG         GG
#> 68          AT         AT         AA         AT         AA         NN
#> 69          AA         AA         AA         NN         AG         AA
#> 70          CC         AA         NN         NN         AC         AC
#> 71          AA         AT         AA         AA         AA         AA
#> 72          AG         GG         AG         AG         NN         AA
#> 73          TT         TT         TT         TT         TT         TT
#> 74          GG         GT         GG         NN         NN         GG
#> 75          NN         AC         AC         AC         AC         CC
#> 76          CC         CT         NN         CC         NN         NN
#> 77          AA         AA         AA         NN         AT         AT
#> 78          GG         GG         NN         NN         GT         NN
#> 79          TT         TT         TT         TT         TT         CT
#> 80          TT         CT         CT         NN         CT         TT
#> 81          AT         NN         AA         NN         NN         AT
#> 82          TT         TT         TT         TT         NN         TT
#> 83          AT         TT         TT         TT         TT         TT
#> 84          CC         CT         NN         NN         CC         NN
#> 85          AG         GG         NN         AG         AG         GG
#> 86          AA         AA         AA         NN         AG         AG
#> 87          TT         NN         TT         CT         TT         TT
#> 88          AG         NN         NN         NN         AG         AA
#> 89          NN         NN         CC         CG         CG         NN
#> 90          CC         CC         CC         NN         CC         CC
#> 91          CT         NN         CT         CT         NN         NN
#> 92          AG         NN         AG         AG         AG         NN
#> 93          TT         TT         TT         TT         NN         TT
#> 94          CC         CC         CC         NN         CC         CG
#> 95          GG         NN         NN         GG         NN         NN
#> 96          CC         CC         AC         NN         NN         CC
#> 97          GG         GG         GG         NN         GG         NN
#> 98          TT         GG         GT         GG         GG         GT
#> 99          GG         NN         AG         AA         AG         AG
#> 100         AT         NN         AT         AT         AT         NN
#>     OPL.r1.110 OPL.r1.125 OPL.r1.126 OPL.r1.128 OPL.r1.139 OPL.r1.146
#> 1           AG         NN         AG         AG         AA         NN
#> 2           GG         GG         GG         GG         GG         NN
#> 3           NN         NN         CC         CC         CC         CC
#> 4           AG         AG         AG         AG         AG         AG
#> 5           GG         GG         NN         GG         NN         AG
#> 6           AA         AA         AA         AA         AA         AA
#> 7           NN         TT         TT         TT         NN         TT
#> 8           NN         AG         NN         AG         GG         AG
#> 9           CC         CC         NN         CC         CC         CC
#> 10          AA         NN         AA         AA         AA         AA
#> 11          NN         CC         NN         CT         CT         CT
#> 12          NN         AA         AA         AA         AG         AG
#> 13          NN         TT         TT         TT         TT         CT
#> 14          TT         NN         GG         GG         NN         GG
#> 15          AA         NN         NN         AA         AA         AG
#> 16          GT         GT         TT         GG         GG         GT
#> 17          AC         AA         AA         AC         AA         AA
#> 18          AC         CC         AC         AC         AC         CC
#> 19          TT         TT         TT         NN         GT         TT
#> 20          GG         GG         NN         AG         GG         GG
#> 21          AA         AC         AC         AC         NN         CC
#> 22          CC         NN         NN         CC         CC         NN
#> 23          GG         GG         GG         GG         GG         GG
#> 24          AA         NN         AA         AG         AA         AG
#> 25          AG         NN         AG         AG         GG         AG
#> 26          TT         TT         NN         TT         TT         TT
#> 27          AC         NN         AC         AA         AC         AC
#> 28          CT         TT         CT         NN         CC         CC
#> 29          NN         NN         AA         AA         AA         AA
#> 30          CG         CG         NN         GG         GG         GG
#> 31          NN         NN         CG         CG         CG         GG
#> 32          AA         AC         AC         AC         NN         AA
#> 33          TT         CC         CT         CT         CC         CC
#> 34          CG         NN         CG         CG         CC         NN
#> 35          AC         AA         AA         AC         NN         NN
#> 36          AA         AA         AA         AA         AG         AA
#> 37          GG         NN         GG         CG         GG         GG
#> 38          NN         GT         GT         GG         GG         GT
#> 39          NN         NN         AT         TT         TT         NN
#> 40          NN         CG         GG         CG         NN         GG
#> 41          NN         GG         GG         AG         AG         GG
#> 42          AA         GG         AG         GG         NN         AG
#> 43          NN         AG         AA         AA         AA         AA
#> 44          AT         TT         TT         TT         TT         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          CC         NN         AC         CC         NN         NN
#> 47          CG         CC         CG         CG         CG         NN
#> 48          NN         TT         CT         CT         TT         CT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         NN         CC         CC         CC
#> 51          CT         CT         TT         TT         NN         CT
#> 52          AA         AA         AG         AG         NN         AG
#> 53          AA         AA         AC         AC         AC         AA
#> 54          GG         GG         GG         GG         GG         GG
#> 55          AA         AG         AG         AG         NN         GG
#> 56          AG         AA         AA         AA         AG         AG
#> 57          NN         GG         GG         GG         AG         AG
#> 58          CT         CC         CC         CC         CT         CT
#> 59          GG         AG         NN         AG         NN         AG
#> 60          AC         AA         AC         AA         AA         NN
#> 61          TT         NN         NN         TT         CT         TT
#> 62          GG         NN         GG         CG         GG         CG
#> 63          NN         AA         AC         AC         NN         AA
#> 64          NN         AT         AT         TT         TT         NN
#> 65          NN         CT         CC         TT         TT         NN
#> 66          AA         NN         NN         AA         AC         AA
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AT         AT         AA         AT         AA         AT
#> 69          AA         AA         AA         AA         AG         AG
#> 70          AA         AA         NN         AA         AC         AC
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AG         NN         AA         AA         AA         AG
#> 73          NN         TT         TT         TT         TT         GT
#> 74          NN         NN         GT         GG         GG         NN
#> 75          NN         AC         AA         AC         AA         AC
#> 76          TT         NN         CC         TT         NN         CT
#> 77          NN         AA         AA         AA         AT         AT
#> 78          NN         GG         GG         GG         GT         GT
#> 79          NN         TT         TT         TT         TT         TT
#> 80          CT         TT         CT         TT         NN         NN
#> 81          AT         AT         AT         AT         AT         AA
#> 82          TT         TT         TT         TT         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          CT         CT         NN         CC         CC         NN
#> 85          AG         AG         GG         AG         AG         AG
#> 86          AG         AA         NN         AA         AG         AG
#> 87          TT         TT         TT         TT         TT         CT
#> 88          AG         AA         AA         AA         AG         AA
#> 89          NN         NN         CG         CG         CC         NN
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CT         CT         CT         NN         TT         TT
#> 92          GG         NN         NN         GG         NN         AG
#> 93          TT         TT         TT         TT         TT         TT
#> 94          CC         CC         CC         CC         CC         NN
#> 95          GT         GG         GG         GG         GG         GG
#> 96          CC         CC         CC         AC         CC         CC
#> 97          GG         GG         GG         GG         GG         GG
#> 98          GT         NN         GG         GG         GT         GG
#> 99          NN         AG         AA         AA         GG         NN
#> 100         NN         AA         NN         AT         NN         TT
#>     OPL.r1.147 OPL.r2.007 OPL.r2.012 OPL.r2.022 OPL.r2.028 OPL.r2.031
#> 1           AG         AA         AG         AA         NN         AG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AG         AG         AG         NN         AG         AG
#> 5           NN         GG         AG         GG         NN         NN
#> 6           AC         AA         AA         AA         AA         AA
#> 7           TT         GT         TT         GT         TT         GT
#> 8           NN         GG         AA         GG         GG         GG
#> 9           CC         AC         AC         CC         CC         CC
#> 10          AA         AA         AA         AA         AG         AA
#> 11          CT         TT         TT         CT         CC         CT
#> 12          AA         AG         AA         AA         NN         AA
#> 13          CT         TT         CT         CT         CT         TT
#> 14          GG         GG         GG         GG         GG         GG
#> 15          AA         NN         AA         AA         GG         NN
#> 16          GT         GT         GT         GG         NN         GG
#> 17          AA         AA         AA         AA         AA         AC
#> 18          CC         CC         CC         AC         CC         AC
#> 19          TT         TT         TT         TT         NN         TT
#> 20          GG         GG         GG         GG         GG         GG
#> 21          CC         AC         NN         CC         AC         NN
#> 22          TT         CC         CT         CT         CC         CT
#> 23          GG         GG         GG         GG         GG         GG
#> 24          AA         AA         AA         AA         AA         AA
#> 25          GG         GG         GG         GG         NN         AG
#> 26          TT         TT         NN         TT         CT         TT
#> 27          AC         NN         AA         AA         AA         AA
#> 28          TT         CC         NN         CT         CC         CT
#> 29          AA         AC         AC         AA         CC         CC
#> 30          GG         GG         GG         GG         GG         GG
#> 31          CC         CG         CG         CG         CC         CC
#> 32          AA         AA         AA         AA         AA         AA
#> 33          NN         NN         NN         CC         CC         NN
#> 34          CC         CG         CG         CG         GG         GG
#> 35          AA         AA         AA         AA         AA         NN
#> 36          AA         AA         AG         GG         AA         AA
#> 37          CG         GG         GG         CG         NN         GG
#> 38          GG         GG         GG         GT         GG         GG
#> 39          TT         AT         TT         TT         NN         TT
#> 40          GG         GG         GG         GG         CG         GG
#> 41          GG         GG         NN         AG         GG         AG
#> 42          AG         AA         AG         AA         GG         NN
#> 43          AA         NN         AA         AG         AG         AG
#> 44          TT         TT         TT         TT         AT         TT
#> 45          NN         GG         NN         GG         NN         GG
#> 46          AC         NN         AC         AC         AA         CC
#> 47          GG         NN         CG         GG         GG         CG
#> 48          CT         CC         CT         CT         CC         CC
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CC         CC         CC         CC
#> 51          CT         CC         NN         NN         TT         NN
#> 52          AA         AG         AG         AG         AA         GG
#> 53          AC         AC         AA         AA         NN         AC
#> 54          GG         GG         GG         GT         GG         GG
#> 55          GG         AG         NN         GG         AG         AG
#> 56          AA         AG         AA         AA         AA         AA
#> 57          GG         GG         GG         GG         NN         GG
#> 58          NN         CC         CT         CT         CC         CC
#> 59          NN         AG         GG         AG         AG         GG
#> 60          AC         AC         AA         AA         NN         AC
#> 61          TT         NN         TT         NN         NN         NN
#> 62          CG         CG         CG         GG         GG         CG
#> 63          AC         AA         NN         AA         AC         AA
#> 64          AT         AT         NN         AT         TT         AT
#> 65          CT         TT         NN         CT         NN         TT
#> 66          AC         AA         AC         AC         NN         AC
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AT         AA         AA         AT         TT         AA
#> 69          AA         GG         AA         AG         NN         AG
#> 70          AA         NN         AA         AC         AC         AC
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AA         AA         AA         AG         AG         AA
#> 73          TT         TT         GT         TT         TT         TT
#> 74          GT         GT         NN         NN         NN         GT
#> 75          AA         AA         AA         AC         AC         CC
#> 76          NN         NN         NN         CC         NN         CC
#> 77          AA         AT         NN         AA         AA         AA
#> 78          GG         TT         GG         GT         GG         GT
#> 79          TT         TT         TT         CT         TT         CT
#> 80          TT         TT         TT         TT         TT         CT
#> 81          AT         AT         AT         AT         TT         AT
#> 82          TT         CT         TT         CT         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          TT         CC         CT         TT         CT         CC
#> 85          AG         AG         AG         AG         GG         GG
#> 86          AA         AG         AA         NN         AA         NN
#> 87          CT         CT         TT         TT         TT         TT
#> 88          AA         AG         NN         AG         AA         AG
#> 89          CC         CC         CC         CC         GG         CC
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CT         TT         CT         CT         CT         CC
#> 92          AG         GG         NN         AG         NN         NN
#> 93          TT         TT         TT         TT         CT         TT
#> 94          CC         NN         CC         CC         CC         CC
#> 95          GG         TT         GG         GG         GG         GT
#> 96          CC         AC         CC         AC         AC         CC
#> 97          AG         GG         GG         GG         NN         GG
#> 98          GG         GG         GG         GG         GG         GG
#> 99          AA         AG         AG         AG         AA         GG
#> 100         AA         AA         AT         AA         AA         TT
#>     OPL.r2.033 OPL.r2.080 OPL.r2.094 PAL.r2.203 PAL.r2.238 PAL.r2.239
#> 1           AG         GG         AG         AG         GG         GG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           NN         AG         AA         AA         AA         AA
#> 5           NN         AG         NN         AG         GG         GG
#> 6           AA         AA         AA         AA         AC         NN
#> 7           TT         TT         NN         TT         TT         TT
#> 8           GG         GG         GG         NN         GG         GG
#> 9           AC         CC         CC         CC         CC         AC
#> 10          AA         AA         AA         AA         AA         AA
#> 11          CT         CT         CT         TT         TT         TT
#> 12          NN         AA         NN         AG         AG         AA
#> 13          TT         TT         CT         CT         TT         TT
#> 14          GG         GG         TT         GG         GG         GG
#> 15          AA         AA         AA         GG         AG         AG
#> 16          GT         GG         NN         GT         TT         NN
#> 17          AA         AA         AA         CC         AA         AA
#> 18          CC         AC         AC         NN         AC         AA
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         GG         GG         NN         GG         GG
#> 21          AC         AC         AA         CC         CC         CC
#> 22          CC         CC         NN         CC         CC         CC
#> 23          GG         GG         NN         GG         GG         GG
#> 24          AA         GG         AA         AA         AA         AA
#> 25          GG         AG         AG         GG         GG         GG
#> 26          TT         CT         NN         TT         CT         CC
#> 27          AA         AA         AA         AA         AA         AA
#> 28          NN         CT         CC         TT         TT         TT
#> 29          NN         AC         AC         AC         AC         AC
#> 30          NN         GG         GG         GG         GG         GG
#> 31          CC         GG         CG         CC         CC         CC
#> 32          AA         AA         AA         AA         AA         AC
#> 33          CT         CT         CT         CC         CC         NN
#> 34          CC         NN         CC         CC         CC         CC
#> 35          AA         AA         AA         AA         AA         AA
#> 36          AA         AA         AG         AA         AA         AA
#> 37          CG         CG         NN         GG         GG         GG
#> 38          GG         GG         GG         GG         GT         GT
#> 39          TT         TT         TT         TT         TT         NN
#> 40          NN         GG         CG         GG         GG         CG
#> 41          AG         GG         AG         AA         AA         AA
#> 42          AG         NN         GG         AG         GG         GG
#> 43          AA         AA         AG         AA         AA         AA
#> 44          TT         TT         AT         AT         AA         AT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          NN         CC         AA         CC         CC         CC
#> 47          GG         CC         NN         CC         CC         CC
#> 48          CT         CT         NN         CT         CT         CT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          NN         CC         CT         CC         CC         CC
#> 51          CT         TT         CT         NN         CT         CT
#> 52          AG         GG         AA         AG         NN         AA
#> 53          NN         AA         AA         AA         AA         AA
#> 54          GG         GG         GG         GG         GG         GG
#> 55          AA         GG         GG         AA         AG         AG
#> 56          NN         NN         NN         AG         AG         AA
#> 57          NN         GG         GG         NN         GG         GG
#> 58          CT         CT         CC         CC         TT         CT
#> 59          GG         GG         AG         GG         GG         GG
#> 60          AA         AC         AC         AA         AA         AA
#> 61          TT         TT         NN         NN         TT         NN
#> 62          CG         GG         GG         GG         GG         GG
#> 63          NN         AA         NN         NN         AC         CC
#> 64          TT         TT         AT         AT         NN         AT
#> 65          CC         TT         CC         CC         CT         CC
#> 66          NN         AC         AC         CC         AC         AA
#> 67          NN         GG         GG         AG         GG         NN
#> 68          AA         AA         AA         TT         NN         NN
#> 69          AA         GG         AA         AA         AA         AA
#> 70          AA         NN         AA         AA         AA         AA
#> 71          AA         AA         AA         AA         AA         AT
#> 72          AA         AG         AG         GG         AG         GG
#> 73          GT         TT         GT         TT         GT         GT
#> 74          GG         GT         GT         GG         GG         GG
#> 75          AC         CC         AC         AA         AC         AA
#> 76          CC         CC         CC         NN         NN         TT
#> 77          AT         AA         AA         AT         AT         AA
#> 78          GG         TT         GG         GG         GG         GG
#> 79          TT         TT         TT         TT         TT         CT
#> 80          CT         CT         CT         TT         CT         NN
#> 81          AA         AT         NN         TT         AT         AT
#> 82          TT         TT         TT         TT         TT         TT
#> 83          TT         TT         NN         TT         TT         TT
#> 84          NN         CC         CT         NN         TT         TT
#> 85          GG         AG         NN         AA         AA         AA
#> 86          AG         AA         AA         AG         AG         NN
#> 87          CC         CT         CT         TT         NN         TT
#> 88          NN         AG         AA         GG         GG         GG
#> 89          CG         CG         NN         CC         CC         NN
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CT         NN         CC         CC         CC         CC
#> 92          NN         AG         GG         AG         GG         AG
#> 93          TT         TT         TT         NN         TT         TT
#> 94          CC         CG         CC         CG         CC         NN
#> 95          GG         GT         GG         GG         GG         GG
#> 96          NN         CC         CC         CC         CC         CC
#> 97          GG         GG         GG         GG         GG         GG
#> 98          GG         GG         GG         GG         GG         GT
#> 99          NN         AG         AG         GG         GG         GG
#> 100         AA         AT         AA         NN         AT         NN
#>     PAL.r2.240 PAL.r2.242 PAL.r2.249 PAL.r2.253 PAL.r2.276 PAL.r2.278
#> 1           GG         AG         GG         GG         AG         GG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AA         AA         NN         AA         AA         AA
#> 5           GG         AG         NN         AG         AG         GG
#> 6           AC         AC         AC         AA         AC         AC
#> 7           TT         TT         TT         TT         TT         TT
#> 8           GG         GG         GG         GG         GG         GG
#> 9           CC         CC         CC         CC         CC         CC
#> 10          AA         AA         AA         AA         AA         AA
#> 11          TT         TT         NN         TT         TT         TT
#> 12          AG         AG         NN         AA         AG         AA
#> 13          TT         TT         TT         TT         TT         TT
#> 14          GG         GG         GG         GG         GG         GG
#> 15          GG         GG         GG         NN         AG         GG
#> 16          GT         GT         GG         TT         GT         GT
#> 17          AA         AA         AA         NN         AA         AA
#> 18          AA         AA         AC         AC         NN         AA
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         GG         GG         GG         GG         NN
#> 21          CC         CC         NN         CC         CC         CC
#> 22          CC         CC         CC         CC         CC         CC
#> 23          GG         NN         GG         GG         GG         GG
#> 24          AA         AA         AA         NN         AA         AA
#> 25          GG         GG         GG         GG         NN         GG
#> 26          CC         TT         CT         CC         CC         NN
#> 27          AA         AA         AA         AA         AA         AA
#> 28          TT         CT         TT         CT         TT         TT
#> 29          AA         AA         AA         AA         AA         AA
#> 30          GG         GG         GG         GG         GG         GG
#> 31          CC         CC         CC         CC         CC         CC
#> 32          AC         AC         AC         AA         AC         AC
#> 33          CC         CC         CT         CC         CT         CT
#> 34          CC         CC         CC         CC         CC         CC
#> 35          AA         AA         AA         AA         AA         NN
#> 36          AA         AA         AA         AA         AA         AA
#> 37          NN         GG         CG         GG         GG         GG
#> 38          GT         GT         GG         GT         GG         GT
#> 39          NN         AT         NN         NN         TT         NN
#> 40          GG         GG         NN         GG         CG         GG
#> 41          AA         NN         AA         AA         NN         AA
#> 42          GG         GG         NN         GG         AG         GG
#> 43          AA         AA         AA         AA         AA         AA
#> 44          AT         AT         TT         AT         NN         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          CC         CC         CC         CC         CC         CC
#> 47          CC         CG         CG         NN         CC         CC
#> 48          TT         TT         NN         CC         CT         NN
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CC         CC         CC         CC
#> 51          NN         TT         TT         TT         TT         NN
#> 52          GG         AG         NN         AG         NN         AG
#> 53          AA         AA         AA         AA         AA         AA
#> 54          GG         GG         GG         GG         GG         GG
#> 55          NN         GG         AG         GG         GG         AG
#> 56          AG         AG         AG         AA         AG         AA
#> 57          GG         NN         GG         GG         GG         GG
#> 58          NN         CT         CT         CT         CC         CC
#> 59          GG         GG         GG         GG         GG         GG
#> 60          AA         AA         AA         AA         AA         AA
#> 61          NN         TT         TT         NN         TT         TT
#> 62          GG         GG         NN         GG         GG         GG
#> 63          AC         NN         NN         NN         CC         CC
#> 64          AA         AA         AT         AA         AA         TT
#> 65          CC         CC         CC         CT         CC         CT
#> 66          AA         AC         CC         CC         CC         AC
#> 67          GG         GG         AG         GG         GG         NN
#> 68          TT         AT         AA         TT         NN         AT
#> 69          AA         AA         AA         AA         AA         AA
#> 70          AA         AA         AA         AA         AA         AA
#> 71          NN         AA         AT         AA         AA         AA
#> 72          GG         AA         NN         GG         NN         GG
#> 73          TT         GT         TT         NN         TT         GT
#> 74          GG         GG         GG         GG         GG         GG
#> 75          AC         AA         AA         AA         AA         AC
#> 76          TT         TT         NN         NN         NN         TT
#> 77          AT         AT         AT         AA         AT         NN
#> 78          GG         GG         GG         GG         GG         GG
#> 79          CT         TT         CT         CT         TT         TT
#> 80          CT         CT         TT         CT         CT         CT
#> 81          AA         AA         AT         AT         AT         AT
#> 82          TT         TT         TT         TT         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          TT         TT         TT         TT         TT         TT
#> 85          AA         NN         AG         AA         AG         NN
#> 86          NN         AG         AG         NN         AG         NN
#> 87          TT         TT         TT         TT         TT         TT
#> 88          GG         GG         GG         GG         GG         GG
#> 89          GG         CG         GG         NN         CG         NN
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CC         CC         CC         CT         CC         CT
#> 92          NN         AA         NN         AG         AA         AG
#> 93          TT         TT         TT         TT         NN         TT
#> 94          CG         NN         NN         CG         CC         GG
#> 95          GG         GG         GG         GG         GG         GG
#> 96          CC         CC         CC         CC         CC         CC
#> 97          GG         GG         GG         GG         GG         GG
#> 98          GT         GG         NN         GT         GT         GG
#> 99          GG         NN         GG         GG         NN         GG
#> 100         TT         NN         AA         NN         AT         AT
#>     PAL.r2.281 PAL.r2.287 PAL.r2.290 PAL.r2.291 SMR.r1.004 SMR.r1.008
#> 1           AG         GG         AG         AA         GG         AA
#> 2           GG         GG         GG         GG         GT         GG
#> 3           NN         CC         CC         CC         CC         CT
#> 4           AA         AA         AA         AA         AA         AA
#> 5           NN         NN         GG         GG         NN         NN
#> 6           AC         AA         AA         AC         AA         AC
#> 7           TT         TT         TT         TT         TT         TT
#> 8           GG         GG         GG         GG         AG         GG
#> 9           CC         CC         NN         CC         AC         CC
#> 10          AA         AA         AA         AA         AA         AA
#> 11          TT         TT         TT         TT         CT         CC
#> 12          AA         AG         AA         AG         AA         AG
#> 13          TT         TT         TT         TT         TT         TT
#> 14          GG         GG         GG         GG         GG         GG
#> 15          AG         GG         GG         GG         AA         AG
#> 16          GG         GT         GT         GG         GG         GT
#> 17          AA         AC         AC         AA         AA         AC
#> 18          NN         AC         AA         AC         CC         CC
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         GG         GG         GG         GG         AG
#> 21          CC         CC         CC         CC         AC         CC
#> 22          CC         CC         CC         CC         CC         CT
#> 23          GG         GG         GG         GG         AG         GG
#> 24          AA         AA         AA         AA         AA         AG
#> 25          GG         GG         GG         GG         AA         AG
#> 26          CT         CT         NN         CT         NN         TT
#> 27          AA         AA         AA         AA         NN         AA
#> 28          NN         TT         CT         TT         TT         CC
#> 29          NN         AA         AC         AA         NN         AC
#> 30          GG         GG         GG         GG         GG         GG
#> 31          CC         CC         CC         CC         CC         GG
#> 32          AC         AC         NN         AC         AC         AC
#> 33          CT         CC         CT         NN         CT         CT
#> 34          CC         CC         CC         CC         CG         CC
#> 35          AA         AA         AA         AA         AA         AA
#> 36          AA         AA         NN         AA         AA         AG
#> 37          GG         GG         GG         GG         GG         CG
#> 38          GT         NN         GT         TT         GT         TT
#> 39          TT         TT         AT         TT         NN         AT
#> 40          NN         GG         GG         GG         GG         GG
#> 41          NN         AA         AA         AA         AG         AA
#> 42          GG         GG         GG         GG         AA         AG
#> 43          AA         AA         AA         AA         AA         AA
#> 44          TT         AT         TT         TT         TT         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          CC         CC         CC         NN         NN         AA
#> 47          CC         CC         CC         CC         CC         CG
#> 48          NN         TT         TT         CC         TT         CC
#> 49          GG         GG         GG         GG         GG         GG
#> 50          NN         CC         CC         CC         CC         CC
#> 51          NN         NN         TT         CT         NN         CT
#> 52          AG         NN         AG         AG         AG         AA
#> 53          AA         AA         AA         AA         AA         AA
#> 54          GG         GG         NN         GG         GG         GG
#> 55          AG         AG         AG         AA         GG         AA
#> 56          NN         AG         AA         AG         NN         AG
#> 57          GG         NN         NN         GG         GG         GG
#> 58          CC         CC         NN         TT         CC         CC
#> 59          GG         GG         AG         NN         AG         NN
#> 60          AA         AA         AA         AA         NN         AC
#> 61          TT         TT         TT         TT         TT         TT
#> 62          GG         GG         GG         GG         GG         CG
#> 63          NN         CC         AC         CC         AC         AA
#> 64          AT         NN         NN         AT         AT         TT
#> 65          CT         CT         CC         CT         CT         CT
#> 66          NN         AC         CC         AC         AA         AC
#> 67          AG         GG         AG         GG         GG         GG
#> 68          AT         AT         NN         TT         AA         AT
#> 69          AA         AA         AA         AA         AA         AA
#> 70          AA         AA         AA         AA         AA         AA
#> 71          AT         AA         AT         AT         AA         AA
#> 72          AA         AG         AG         AG         AA         AG
#> 73          GG         NN         GT         TT         TT         TT
#> 74          GG         GG         GG         GG         GG         GT
#> 75          NN         AC         AC         AC         CC         AA
#> 76          NN         TT         NN         TT         NN         CC
#> 77          AA         AT         AA         AT         NN         AT
#> 78          GG         GG         GG         GG         GG         GG
#> 79          CT         CT         TT         CT         TT         CT
#> 80          CC         CT         NN         CT         TT         TT
#> 81          TT         AA         AT         AT         AA         AT
#> 82          TT         TT         TT         TT         NN         CT
#> 83          TT         NN         TT         TT         TT         TT
#> 84          TT         TT         TT         TT         TT         CC
#> 85          AA         AA         AG         NN         GG         AG
#> 86          AA         AG         AA         AG         AA         AG
#> 87          TT         TT         TT         TT         TT         TT
#> 88          GG         GG         NN         NN         AA         AG
#> 89          CC         NN         GG         GG         NN         CC
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CC         CC         CC         NN         TT         CC
#> 92          NN         AG         NN         AA         NN         AG
#> 93          CT         TT         CT         NN         TT         TT
#> 94          CG         CG         CC         CC         CC         CC
#> 95          GG         GG         GG         GG         GT         GT
#> 96          CC         CC         CC         CC         CC         CC
#> 97          GG         GG         GG         GG         GG         GG
#> 98          GT         TT         GT         TT         GG         GG
#> 99          GG         GG         GG         NN         GG         AG
#> 100         AA         AA         NN         AA         AT         AT
#>     SMR.r1.012 SMR.r1.023 SMR.r1.028 SMR.r1.034 SMR.r1.052 SMR.r1.056
#> 1           AG         AG         AA         GG         AG         GG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         NN         CC         CC         CC
#> 4           AG         AG         AA         AG         AG         AA
#> 5           AG         NN         NN         AG         GG         GG
#> 6           NN         AA         NN         AA         AA         AA
#> 7           NN         GT         TT         TT         NN         TT
#> 8           NN         AG         AG         GG         NN         AA
#> 9           NN         AC         AC         AC         CC         AC
#> 10          AA         AA         AG         NN         AG         AA
#> 11          CT         CT         NN         CT         CT         TT
#> 12          AA         AA         NN         AG         AA         NN
#> 13          TT         TT         CT         CT         CT         TT
#> 14          NN         GG         GT         GG         GT         GG
#> 15          AG         AG         AA         AA         AG         NN
#> 16          GG         GT         GT         GG         GT         NN
#> 17          AC         AA         AA         AA         AA         NN
#> 18          CC         CC         CC         CC         CC         CC
#> 19          TT         TT         TT         TT         TT         TT
#> 20          AG         GG         GG         GG         GG         AG
#> 21          AC         CC         CC         AC         NN         CC
#> 22          CC         NN         CT         CC         CC         CC
#> 23          GG         AG         GG         GG         AG         GG
#> 24          AA         AA         AA         AA         AA         AG
#> 25          AA         NN         GG         NN         NN         AA
#> 26          TT         CT         CT         TT         CT         CT
#> 27          AA         NN         AC         NN         AA         AA
#> 28          CC         CC         NN         CT         TT         NN
#> 29          AA         AA         NN         NN         NN         NN
#> 30          CG         GG         GG         NN         GG         CG
#> 31          CG         CG         CC         CC         NN         CG
#> 32          AA         AA         AC         NN         AA         AA
#> 33          CC         CT         CC         CT         NN         CC
#> 34          CC         NN         CC         CG         CC         GG
#> 35          AA         AA         AC         AA         AA         AC
#> 36          AA         AA         AA         AA         AA         AA
#> 37          GG         CG         CG         GG         GG         NN
#> 38          GT         NN         NN         GG         GT         GT
#> 39          TT         TT         NN         TT         NN         NN
#> 40          GG         CG         NN         GG         GG         GG
#> 41          NN         NN         NN         GG         AG         GG
#> 42          GG         GG         NN         AG         AG         AG
#> 43          AA         AA         AG         AA         AA         AG
#> 44          NN         NN         TT         TT         TT         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          CC         AC         CC         NN         NN         CC
#> 47          CC         CG         CG         CC         CC         CC
#> 48          CT         TT         CT         NN         NN         TT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         NN         CC         CC         CC         NN
#> 51          CT         NN         CT         NN         TT         CT
#> 52          AG         AA         GG         AA         AG         AA
#> 53          AA         AA         AC         NN         AA         CC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          GG         GG         GG         GG         GG         GG
#> 56          AA         NN         AA         AG         AA         NN
#> 57          GG         GG         NN         GG         GG         NN
#> 58          CC         NN         CC         NN         NN         CC
#> 59          NN         AG         GG         GG         AG         AG
#> 60          AC         AA         AC         AA         AC         AC
#> 61          TT         TT         NN         NN         CT         TT
#> 62          GG         GG         GG         GG         GG         GG
#> 63          AA         AC         NN         AC         NN         NN
#> 64          AA         AT         TT         AT         AT         AA
#> 65          CT         CT         CC         TT         CT         CC
#> 66          AC         NN         NN         AA         AC         AA
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AA         AT         AA         AA         AA         NN
#> 69          AA         AG         AA         AA         AA         NN
#> 70          NN         AC         AA         AC         AA         AA
#> 71          AA         AA         AT         AA         AA         AA
#> 72          AA         NN         AA         AG         NN         AG
#> 73          NN         GT         TT         GT         TT         TT
#> 74          GG         GT         GT         NN         NN         GG
#> 75          NN         AC         AC         AC         AC         NN
#> 76          CC         CC         CC         NN         NN         CC
#> 77          AA         AA         AA         AT         NN         AA
#> 78          GG         GT         GG         GG         GG         GG
#> 79          TT         CT         TT         TT         TT         TT
#> 80          CT         CT         TT         CT         CT         CT
#> 81          AA         NN         AT         AT         AT         AT
#> 82          TT         TT         TT         TT         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          CT         NN         NN         CC         NN         CC
#> 85          GG         NN         GG         NN         AG         AG
#> 86          AA         AA         AA         AG         NN         AA
#> 87          NN         NN         CT         CT         TT         TT
#> 88          AG         AG         NN         NN         AG         AG
#> 89          NN         CG         NN         CC         CG         NN
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CT         NN         CC         CT         CC         CT
#> 92          GG         NN         NN         NN         GG         NN
#> 93          NN         NN         NN         CT         TT         TT
#> 94          CC         CC         CC         CC         NN         NN
#> 95          GG         GT         NN         GG         GG         NN
#> 96          CC         AC         CC         CC         CC         CC
#> 97          GG         GG         GG         GG         GG         NN
#> 98          GT         NN         GT         GG         GG         GT
#> 99          AG         AA         NN         AG         GG         AG
#> 100         AA         NN         NN         AT         AA         AT
#>     SMR.r1.072 SMR.r2.001 SMR.r2.077 SMR.r2.086 SMR.r2.097 SMR.r2.105
#> 1           AA         AG         AG         AA         NN         AA
#> 2           GG         GG         NN         GG         GG         GG
#> 3           CC         CC         CC         CC         CT         CC
#> 4           AA         AA         AA         AG         AA         GG
#> 5           NN         GG         NN         GG         GG         GG
#> 6           AA         AA         AA         AA         AA         AA
#> 7           TT         GT         TT         GT         TT         TT
#> 8           AG         AG         AG         NN         GG         AA
#> 9           NN         CC         CC         AC         AC         NN
#> 10          AA         AA         AA         AA         AA         AA
#> 11          CT         NN         NN         TT         CT         CT
#> 12          AA         NN         AA         AG         AG         AA
#> 13          CT         NN         TT         CC         CT         TT
#> 14          GT         GG         GT         GT         GG         GT
#> 15          AA         AG         AA         AA         AA         AG
#> 16          GT         GG         GT         GT         GG         NN
#> 17          AC         AA         AA         AC         AA         AA
#> 18          AC         CC         AC         AA         AC         CC
#> 19          GT         GT         NN         TT         TT         TT
#> 20          GG         NN         GG         GG         GG         GG
#> 21          AC         NN         AC         CC         AC         AC
#> 22          NN         CC         CC         NN         CT         CC
#> 23          NN         NN         AG         GG         GG         GG
#> 24          AA         NN         NN         AA         AA         AA
#> 25          AG         AG         NN         AG         AA         AG
#> 26          CT         CT         CT         TT         TT         TT
#> 27          AC         AC         AA         AA         NN         AA
#> 28          NN         CT         CT         NN         CT         CT
#> 29          AA         NN         NN         AA         CC         AA
#> 30          GG         GG         CG         GG         GG         GG
#> 31          CC         CG         CG         CG         CC         CG
#> 32          AA         AC         AA         NN         AA         AA
#> 33          CC         NN         CC         CC         CT         TT
#> 34          CC         NN         CC         CG         CC         CG
#> 35          NN         AA         AA         AA         AA         AA
#> 36          AG         AA         AA         AG         NN         NN
#> 37          CG         GG         CG         GG         CG         CG
#> 38          GT         GG         GG         GT         NN         GT
#> 39          TT         TT         TT         TT         TT         TT
#> 40          CG         NN         GG         NN         CG         GG
#> 41          AG         NN         NN         GG         AG         GG
#> 42          AA         AA         AG         AA         AG         AG
#> 43          AA         AA         AG         AA         AA         AA
#> 44          TT         TT         TT         TT         AT         AT
#> 45          GG         NN         NN         GG         GG         GG
#> 46          AC         NN         CC         CC         NN         CC
#> 47          CG         CC         NN         CG         GG         CC
#> 48          CT         NN         NN         CC         CT         TT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          NN         CC         CC         NN         CC         CC
#> 51          NN         CT         TT         CT         NN         CC
#> 52          AA         AA         AA         AA         AA         AG
#> 53          AA         AA         AA         AA         AC         AC
#> 54          GG         NN         GG         GT         GG         GG
#> 55          GG         NN         NN         GG         NN         AG
#> 56          NN         AA         AA         AG         AG         AA
#> 57          GG         AG         GG         GG         GG         GG
#> 58          CC         NN         CC         CC         CT         CT
#> 59          AA         AG         GG         AG         GG         AG
#> 60          CC         CC         AC         AA         AC         AA
#> 61          TT         NN         TT         TT         TT         NN
#> 62          GG         GG         GG         GG         GG         GG
#> 63          AA         NN         NN         AC         CC         AA
#> 64          AT         TT         NN         AT         TT         AT
#> 65          CC         CT         CC         CT         CC         TT
#> 66          NN         CC         AC         CC         AA         AC
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AT         AA         AA         AT         TT         AA
#> 69          AA         AA         AA         AG         AA         NN
#> 70          AA         AC         AA         AC         AA         AC
#> 71          AA         AA         AA         AT         AT         AT
#> 72          AA         AG         AA         AG         AA         AG
#> 73          TT         GT         TT         TT         NN         GT
#> 74          GG         GG         NN         GG         GT         GG
#> 75          AC         CC         AA         AC         AA         CC
#> 76          CT         NN         NN         CC         CC         NN
#> 77          AA         NN         NN         AT         AT         AA
#> 78          GG         GG         GG         GT         GG         NN
#> 79          TT         TT         TT         TT         TT         TT
#> 80          CT         CT         CT         CT         NN         TT
#> 81          AA         AA         AT         AT         AT         AA
#> 82          TT         NN         CT         TT         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          CT         CC         CC         CT         CC         CC
#> 85          AG         AG         GG         AG         GG         GG
#> 86          AA         NN         NN         AG         AG         AA
#> 87          TT         CT         TT         CT         CT         CT
#> 88          AG         NN         NN         GG         AA         NN
#> 89          CG         CC         CG         CG         CC         CC
#> 90          CC         CC         CC         CC         CC         CT
#> 91          TT         NN         NN         CT         TT         TT
#> 92          NN         NN         AG         NN         GG         GG
#> 93          CT         TT         CT         TT         TT         CT
#> 94          CC         CC         CC         CC         CC         CC
#> 95          GT         NN         GT         GG         GG         GG
#> 96          CC         CC         CC         AC         NN         CC
#> 97          GG         GG         GG         GG         NN         GG
#> 98          GG         GG         GT         GT         GG         GG
#> 99          AA         AA         AG         GG         AA         AG
#> 100         AT         AT         AT         AT         AT         AT
#>     SMR.r2.139 SMR.r2.256 SMR.r2.263 SMR.r2.296 UPD.r2.121 UPD.r2.122
#> 1           AG         AG         AA         GG         AA         AG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CT         CC         CC         CT         CC
#> 4           AA         AG         GG         AA         GG         AG
#> 5           NN         AG         GG         AG         NN         AA
#> 6           AA         AA         AA         AA         AC         AC
#> 7           GT         TT         GT         GT         TT         TT
#> 8           NN         GG         AG         GG         GG         AG
#> 9           CC         CC         CC         CC         CC         AC
#> 10          AA         AA         AA         AA         AA         AA
#> 11          NN         CC         CT         TT         CT         TT
#> 12          NN         AG         AA         AA         AA         AG
#> 13          TT         NN         TT         TT         CT         TT
#> 14          GG         GG         GT         GG         GT         GT
#> 15          AG         AA         NN         AG         AG         GG
#> 16          GT         GG         GT         GT         GG         GT
#> 17          AA         NN         AC         AA         AA         AA
#> 18          CC         CC         CC         CC         CC         AA
#> 19          GT         TT         TT         TT         TT         TT
#> 20          GG         GG         GG         GG         GG         NN
#> 21          NN         CC         CC         AC         AC         AC
#> 22          CT         CT         CC         CC         TT         CT
#> 23          GG         NN         GG         GG         GG         GG
#> 24          NN         AA         AA         AA         AA         AA
#> 25          GG         AG         AG         AG         AG         AG
#> 26          CT         TT         NN         CT         CT         CT
#> 27          AA         AC         AC         AA         AA         AC
#> 28          TT         CT         CC         CC         CT         CT
#> 29          NN         AA         AA         AA         AA         AA
#> 30          NN         GG         GG         GG         GG         GG
#> 31          CG         CC         CC         CG         CG         CG
#> 32          AA         NN         NN         AC         AA         AA
#> 33          CC         CC         CT         CT         TT         CC
#> 34          GG         CC         CG         GG         CC         CG
#> 35          NN         AA         AA         AA         AA         AA
#> 36          AG         AA         AG         AG         AA         AA
#> 37          NN         GG         CG         CC         GG         GG
#> 38          GT         GT         TT         NN         GG         GT
#> 39          TT         AT         TT         TT         AT         TT
#> 40          NN         CG         CG         CG         CG         GG
#> 41          NN         NN         GG         AG         NN         AA
#> 42          NN         AA         AG         AG         AG         GG
#> 43          AA         AA         AA         AG         AA         AG
#> 44          TT         TT         TT         TT         AT         TT
#> 45          GG         GG         GG         NN         GG         GG
#> 46          CC         CC         CC         AA         AC         CC
#> 47          CC         CC         CG         GG         CC         CG
#> 48          CT         NN         CC         NN         NN         TT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CC         CC         CC         CT
#> 51          TT         NN         NN         TT         TT         CC
#> 52          AA         AA         AA         AA         GG         AA
#> 53          AA         AC         AC         AA         CC         AA
#> 54          GG         GG         GG         GG         GG         GG
#> 55          GG         GG         GG         AG         AG         AG
#> 56          AA         AG         NN         AA         AA         AG
#> 57          NN         GG         GG         NN         GG         GG
#> 58          CC         CT         CC         NN         TT         CT
#> 59          AG         AA         AG         AA         AG         GG
#> 60          AC         NN         AA         AA         AA         NN
#> 61          TT         CT         TT         TT         CT         TT
#> 62          CG         GG         GG         GG         NN         GG
#> 63          AC         AA         AC         NN         AC         AC
#> 64          AT         AT         TT         TT         TT         TT
#> 65          CT         CC         TT         CC         CC         TT
#> 66          NN         AA         CC         AC         AC         AA
#> 67          GG         GG         NN         GG         GG         GG
#> 68          AT         NN         NN         AA         TT         AT
#> 69          AG         AA         AG         NN         AA         AA
#> 70          AC         AA         AC         AC         AC         NN
#> 71          AA         AA         AT         AA         AA         AA
#> 72          AA         AA         AA         AG         GG         AG
#> 73          TT         TT         TT         NN         TT         TT
#> 74          NN         GT         GG         GG         GG         GG
#> 75          AC         NN         AC         AC         AC         AC
#> 76          CC         NN         CT         NN         CC         CC
#> 77          AA         AT         AA         AA         AA         AT
#> 78          GT         GG         GT         GT         GG         GG
#> 79          TT         TT         TT         TT         TT         TT
#> 80          TT         TT         CT         TT         CT         NN
#> 81          AT         AT         AT         AT         AA         AT
#> 82          TT         TT         TT         TT         TT         TT
#> 83          TT         TT         TT         AT         TT         TT
#> 84          CC         CC         CT         CT         TT         CC
#> 85          AG         GG         AA         NN         GG         AG
#> 86          NN         AG         AA         NN         AA         AG
#> 87          CT         TT         TT         CT         TT         NN
#> 88          AG         AG         AG         NN         AG         AG
#> 89          CG         NN         GG         CG         CC         CC
#> 90          CC         CC         CC         NN         CC         CC
#> 91          TT         NN         CT         CT         CT         CC
#> 92          NN         GG         NN         GG         GG         GG
#> 93          NN         CT         TT         CT         TT         TT
#> 94          CC         CC         CC         CC         CC         CC
#> 95          NN         GG         GG         GG         GT         GG
#> 96          CC         CC         AC         CC         CC         CC
#> 97          GG         GG         GG         AG         GG         GG
#> 98          GT         GG         GT         GT         GG         GG
#> 99          GG         NN         GG         AA         GG         AA
#> 100         AA         AT         AA         AA         AT         AA
#>     UPD.r2.123 UPD.r2.129 UPD.r2.130 UPD.r2.131 UPD.r2.135 UPD.r2.160
#> 1           AG         AA         AG         AG         AG         NN
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AG         AG         AG         NN         AG         AG
#> 5           AG         NN         NN         NN         AG         AG
#> 6           AC         AC         NN         AA         AC         AA
#> 7           TT         TT         TT         TT         TT         TT
#> 8           AG         GG         NN         NN         AG         AG
#> 9           CC         CC         AC         CC         AC         CC
#> 10          AA         NN         AG         AA         AA         AA
#> 11          TT         NN         CT         CT         CT         CT
#> 12          AG         NN         AA         AG         NN         NN
#> 13          TT         TT         TT         TT         CT         CT
#> 14          GT         GT         GG         GG         GG         TT
#> 15          NN         AG         AG         NN         AG         AG
#> 16          GT         GT         NN         GG         GG         GG
#> 17          AC         NN         AA         AA         AC         NN
#> 18          CC         CC         AA         CC         AC         CC
#> 19          TT         TT         TT         TT         TT         TT
#> 20          AG         GG         GG         GG         GG         GG
#> 21          NN         NN         AC         CC         CC         AC
#> 22          CT         CT         CC         NN         CT         CT
#> 23          GG         GG         GG         GG         GG         GG
#> 24          AA         AA         AA         AG         AA         NN
#> 25          AA         AG         AG         GG         GG         AG
#> 26          TT         TT         CT         CT         CT         NN
#> 27          AA         AC         AA         AA         AA         AA
#> 28          TT         TT         CT         CT         TT         TT
#> 29          AC         NN         AA         AA         AA         AA
#> 30          NN         GG         CG         CG         GG         NN
#> 31          CC         CC         CG         CG         CC         GG
#> 32          AC         AA         AA         AC         AC         AA
#> 33          TT         CT         CC         NN         CC         CC
#> 34          CC         NN         CC         CC         CC         NN
#> 35          AA         AA         AA         AA         AA         AA
#> 36          AG         AA         AA         AG         NN         AA
#> 37          GG         CG         GG         CG         CG         CG
#> 38          GG         GG         GG         TT         GG         GT
#> 39          TT         NN         NN         AT         TT         AT
#> 40          CG         NN         CG         CG         CC         GG
#> 41          NN         NN         AG         AA         AG         NN
#> 42          GG         GG         AG         AA         NN         AG
#> 43          AG         AG         AA         AG         AG         AG
#> 44          NN         AT         AT         NN         NN         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          NN         CC         AC         AA         CC         AC
#> 47          CC         CC         CC         CG         CC         CC
#> 48          CT         CT         TT         CT         CT         TT
#> 49          GG         GG         NN         CG         GG         NN
#> 50          CC         CC         CT         CC         CC         NN
#> 51          CT         CT         CT         CT         CT         NN
#> 52          AG         NN         AG         AG         AA         AG
#> 53          AC         NN         AA         AC         AC         AC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          AG         GG         AG         AA         NN         NN
#> 56          AG         AG         AA         AG         AA         AA
#> 57          GG         NN         GG         NN         NN         GG
#> 58          NN         CT         CT         CC         CC         NN
#> 59          AG         GG         GG         AG         AG         AA
#> 60          AC         AC         AC         AA         AA         AC
#> 61          TT         NN         TT         TT         NN         TT
#> 62          GG         GG         CG         CG         GG         GG
#> 63          AA         AA         AC         AC         NN         AA
#> 64          AT         TT         TT         AT         AT         AT
#> 65          CC         CC         CT         NN         CC         CC
#> 66          AA         NN         AA         AA         CC         CC
#> 67          GG         GG         GG         GG         NN         AG
#> 68          AA         AT         AA         AA         AT         NN
#> 69          AA         AA         AA         AA         AA         AA
#> 70          AC         AA         AC         AC         AA         AC
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AG         AG         GG         NN         AG         GG
#> 73          GT         TT         GG         TT         GT         TT
#> 74          GG         GG         NN         GG         GG         GT
#> 75          AC         AC         AC         CC         AC         CC
#> 76          CT         CT         CT         CC         TT         CT
#> 77          AT         NN         AA         AT         AA         AA
#> 78          GG         GG         GG         GG         GG         GG
#> 79          TT         TT         TT         TT         NN         TT
#> 80          CT         CT         CT         TT         TT         NN
#> 81          NN         AT         AT         AA         AT         AT
#> 82          CT         TT         NN         CT         CT         CT
#> 83          TT         TT         TT         AT         TT         TT
#> 84          CT         CT         CT         CT         CT         NN
#> 85          GG         AA         AG         AG         AG         AG
#> 86          AG         AG         AA         AG         NN         AA
#> 87          TT         CT         CT         CC         CT         CT
#> 88          AG         AA         AA         AG         NN         AG
#> 89          CG         GG         CG         NN         NN         CG
#> 90          CT         CC         CC         CT         CC         CC
#> 91          TT         CT         CC         CC         CC         CT
#> 92          GG         NN         NN         NN         GG         AG
#> 93          CT         TT         TT         CT         TT         CT
#> 94          CC         CC         CG         CC         NN         CG
#> 95          GG         GT         GT         GT         GG         NN
#> 96          NN         CC         CC         CC         AC         CC
#> 97          GG         GG         GG         AG         GG         GG
#> 98          GT         NN         GG         GG         NN         GT
#> 99          NN         GG         AG         AG         NN         NN
#> 100         AT         AT         TT         TT         NN         TT
#>     UPD.r2.162 UPD.r2.189 UPD.r2.201
#> 1           GG         GG         AG
#> 2           GG         GG         GG
#> 3           CC         CC         CC
#> 4           AG         AA         AA
#> 5           NN         NN         GG
#> 6           AC         AA         AA
#> 7           NN         TT         TT
#> 8           NN         GG         GG
#> 9           NN         AC         AA
#> 10          AA         AA         AG
#> 11          CT         CT         TT
#> 12          AG         AG         NN
#> 13          TT         NN         NN
#> 14          GT         GT         GT
#> 15          AG         AA         AA
#> 16          GG         NN         GG
#> 17          AA         AA         AA
#> 18          AC         AC         AC
#> 19          TT         TT         TT
#> 20          GG         AG         GG
#> 21          NN         AA         CC
#> 22          CT         CC         CC
#> 23          GG         GG         GG
#> 24          AG         AA         AA
#> 25          AG         AG         GG
#> 26          NN         TT         TT
#> 27          AA         AA         AA
#> 28          CT         CT         CT
#> 29          AA         AC         AC
#> 30          GG         GG         CG
#> 31          CG         CC         CG
#> 32          NN         NN         AC
#> 33          CC         NN         TT
#> 34          CG         CC         NN
#> 35          NN         AA         AC
#> 36          NN         AA         AA
#> 37          NN         NN         GG
#> 38          GG         GG         NN
#> 39          TT         TT         NN
#> 40          GG         CC         NN
#> 41          NN         AG         AG
#> 42          AG         AG         AG
#> 43          AA         AA         AA
#> 44          TT         AT         TT
#> 45          GG         GG         GT
#> 46          CC         NN         AC
#> 47          GG         CC         CG
#> 48          CT         CC         TT
#> 49          GG         GG         GG
#> 50          CC         CC         CC
#> 51          CT         TT         TT
#> 52          AG         AA         AA
#> 53          NN         AC         CC
#> 54          GG         GG         GG
#> 55          AA         AG         GG
#> 56          AG         AG         AA
#> 57          GG         NN         GG
#> 58          CC         CT         CC
#> 59          AG         NN         NN
#> 60          AC         AA         AC
#> 61          TT         TT         NN
#> 62          GG         GG         GG
#> 63          AA         NN         AA
#> 64          TT         AT         AT
#> 65          CT         CT         TT
#> 66          AC         NN         AC
#> 67          GG         GG         GG
#> 68          TT         AA         AT
#> 69          AA         AA         AA
#> 70          NN         AC         AC
#> 71          AA         AA         AA
#> 72          GG         AG         AG
#> 73          TT         NN         GT
#> 74          GG         NN         GG
#> 75          CC         AC         CC
#> 76          TT         NN         TT
#> 77          AT         AT         AA
#> 78          GG         GG         GG
#> 79          CT         TT         TT
#> 80          CT         TT         TT
#> 81          TT         AT         AT
#> 82          TT         NN         CT
#> 83          TT         NN         TT
#> 84          CC         NN         CC
#> 85          AG         AG         AG
#> 86          AG         AG         AA
#> 87          CT         CC         CT
#> 88          NN         AA         AG
#> 89          GG         GG         CC
#> 90          CC         NN         CC
#> 91          NN         CT         CT
#> 92          GG         AG         GG
#> 93          NN         CT         CT
#> 94          CC         CC         CC
#> 95          GG         NN         GG
#> 96          CC         AC         CC
#> 97          GG         GG         GG
#> 98          GG         GT         TT
#> 99          AA         AG         AG
#> 100         TT         AT         AT

# negative subscripts: remove individuals in pop PAL AND fam B or C, 
# and keep only SNPs on LGIV
subset_snpR_data(stickSNPs, pop.fam = -c("PAL.B", "PAL.C"),
                 chr = "groupIV")
#>   ASP.r2.141 ASP.r2.144 ASP.r2.146 ASP.r2.147 ASP.r2.196 ASP.r2.216 ASP.r2.221
#> 1         AA         AA         AA         AA         AA         AA         AA
#> 2         AT         TT         AT         TT         TT         AT         TT
#> 3         CC         NN         CC         NN         CC         CC         CC
#> 4         AT         NN         AT         AT         TT         AT         TT
#>   ASP.r2.227 ASP.r2.230 ASP.r2.232 ASP.r2.235 CLF.r2.040 CLF.r2.048 CLF.r2.061
#> 1         AA         NN         AC         AA         AA         AC         AA
#> 2         TT         TT         AT         TT         AT         AT         TT
#> 3         CC         AC         CC         NN         CC         CC         CC
#> 4         TT         TT         AA         AA         AT         AT         NN
#>   CLF.r2.064 CLF.r2.068 CLF.r2.069 OPL.r1.078 OPL.r1.080 OPL.r1.083 OPL.r1.086
#> 1         AA         AA         AA         AC         AC         AA         AC
#> 2         AT         AT         AT         TT         TT         TT         TT
#> 3         AC         CC         CC         AC         CC         AC         NN
#> 4         AA         NN         AT         NN         NN         AT         AT
#>   OPL.r1.087 OPL.r1.088 OPL.r1.109 OPL.r1.110 OPL.r1.116 OPL.r1.125 OPL.r1.126
#> 1         AA         AA         NN         AC         AA         NN         AC
#> 2         TT         TT         TT         TT         TT         TT         TT
#> 3         NN         NN         CC         CC         CC         CC         CC
#> 4         AT         AT         NN         NN         NN         AA         NN
#>   OPL.r1.128 OPL.r1.130 OPL.r1.139 OPL.r1.146 OPL.r1.147 OPL.r1.149 OPL.r2.007
#> 1         AA         AA         AC         AC         AC         AA         NN
#> 2         TT         TT         TT         TT         TT         TT         TT
#> 3         AC         CC         CC         CC         CC         CC         AC
#> 4         AT         AA         NN         TT         AA         AT         AA
#>   OPL.r2.012 OPL.r2.022 OPL.r2.023 OPL.r2.028 OPL.r2.031 OPL.r2.033 OPL.r2.079
#> 1         AA         AA         AC         AA         AA         AA         AC
#> 2         TT         TT         TT         TT         TT         TT         TT
#> 3         CC         AC         AC         AC         CC         NN         AC
#> 4         AT         AA         AA         AA         TT         AA         AA
#>   OPL.r2.080 OPL.r2.094 PAL.r2.203 PAL.r2.239 PAL.r2.240 PAL.r2.249 PAL.r2.253
#> 1         AA         AA         AA         AA         AA         AA         AA
#> 2         TT         NN         TT         TT         TT         TT         TT
#> 3         CC         CC         CC         CC         CC         CC         CC
#> 4         AT         AA         NN         NN         TT         AA         NN
#>   PAL.r2.278 PAL.r2.281 PAL.r2.290 PAL.r2.291 SMR.r1.003 SMR.r1.004 SMR.r1.008
#> 1         AA         AA         AA         AA         AA         NN         AA
#> 2         TT         TT         TT         TT         TT         TT         TT
#> 3         CC         CC         CC         CC         AC         CC         CC
#> 4         AT         AA         NN         AA         NN         AT         AT
#>   SMR.r1.012 SMR.r1.013 SMR.r1.023 SMR.r1.028 SMR.r1.034 SMR.r1.049 SMR.r1.052
#> 1         AA         AC         NN         AC         NN         AA         AA
#> 2         TT         TT         TT         TT         TT         TT         TT
#> 3         CC         CC         AC         CC         CC         CC         CC
#> 4         AA         AT         NN         NN         AT         AA         AA
#>   SMR.r1.056 SMR.r1.072 SMR.r1.074 SMR.r2.001 SMR.r2.077 SMR.r2.086 SMR.r2.091
#> 1         AA         AC         AA         AC         AA         AA         AA
#> 2         TT         TT         TT         TT         TT         TT         TT
#> 3         CC         CC         CC         CC         CC         AC         CC
#> 4         AT         AT         NN         AT         AT         AT         AA
#>   SMR.r2.097 SMR.r2.105 SMR.r2.139 SMR.r2.175 SMR.r2.256 SMR.r2.263 SMR.r2.296
#> 1         NN         AA         AA         AA         AC         AC         AA
#> 2         TT         TT         TT         TT         TT         TT         AT
#> 3         NN         CC         CC         CC         CC         AC         CC
#> 4         AT         AT         AA         AT         AT         AA         AA
#>   UPD.r2.120 UPD.r2.121 UPD.r2.122 UPD.r2.123 UPD.r2.124 UPD.r2.129 UPD.r2.130
#> 1         AA         AA         AC         AA         AA         AC         AA
#> 2         TT         TT         TT         TT         TT         TT         TT
#> 3         CC         CC         CC         NN         CC         CC         CC
#> 4         AA         AT         AA         AT         AT         AT         TT
#>   UPD.r2.131 UPD.r2.133 UPD.r2.135 UPD.r2.160 UPD.r2.162 UPD.r2.174 UPD.r2.189
#> 1         AA         AA         AA         AA         AA         NN         AA
#> 2         AT         TT         TT         TT         TT         TT         NN
#> 3         CC         NN         AC         CC         CC         CC         AC
#> 4         TT         TT         NN         TT         TT         AT         AT
#>   UPD.r2.201
#> 1         AA
#> 2         TT
#> 3         CC
#> 4         AT
                 
# using the bracket operator, same as example 1
stickSNPs[pop = c("ASP", "PAL"), chr = c("groupIX", "groupIV")]
#>    ASP.r2.141 ASP.r2.144 ASP.r2.146 ASP.r2.147 ASP.r2.196 ASP.r2.216 ASP.r2.221
#> 1          TT         TT         TT         TT         GT         TT         TT
#> 2          CT         NN         TT         CC         CT         NN         CC
#> 3          AA         AA         AA         NN         AA         AA         AA
#> 4          AA         AA         AA         AA         AA         AA         AA
#> 5          TT         NN         NN         NN         CT         CT         TT
#> 6          AA         AA         AA         AA         AA         AG         AA
#> 7          AG         AA         AG         AG         AA         AA         AG
#> 8          GG         GG         GG         GG         GG         GG         GG
#> 9          CT         CC         CT         CC         CC         CT         CT
#> 10         AA         AA         NN         AA         AA         AA         AA
#> 11         AC         AC         AA         AC         AC         NN         AC
#> 12         GG         GG         GG         GG         GG         GG         GG
#> 13         AT         TT         AT         TT         TT         AT         TT
#> 14         CT         CC         CT         CT         TT         CT         CC
#> 15         AG         NN         AG         NN         GG         AG         AA
#> 16         CC         NN         CC         NN         CC         CC         CC
#> 17         AT         NN         AT         AT         TT         AT         TT
#>    ASP.r2.227 ASP.r2.230 ASP.r2.232 ASP.r2.235 PAL.r2.203 PAL.r2.237 PAL.r2.238
#> 1          TT         TT         TT         TT         TT         TT         TT
#> 2          TT         CT         CC         CT         CC         CC         CC
#> 3          AA         AA         AA         AG         AA         AA         AA
#> 4          AA         NN         AC         AA         AA         AA         AA
#> 5          TT         TT         TT         TT         TT         TT         TT
#> 6          AG         AA         AG         NN         AA         AA         AA
#> 7          AG         GG         AG         AG         AA         AA         AA
#> 8          GG         GG         GG         GG         GG         GG         GG
#> 9          CC         CC         CC         CC         CC         CC         CT
#> 10         AA         AA         AA         AA         AA         AA         AA
#> 11         AC         AA         AA         AC         AA         AA         AA
#> 12         GG         NN         GG         GG         GG         GG         GG
#> 13         TT         TT         AT         TT         TT         TT         TT
#> 14         CT         TT         TT         CC         TT         TT         NN
#> 15         AG         AG         AG         NN         GG         NN         GG
#> 16         CC         AC         CC         NN         CC         CC         CC
#> 17         TT         TT         AA         AA         NN         NN         AT
#>    PAL.r2.239 PAL.r2.240 PAL.r2.241 PAL.r2.242 PAL.r2.249 PAL.r2.253 PAL.r2.275
#> 1          TT         TT         TT         TT         TT         TT         TT
#> 2          CC         CC         CC         CC         CC         CC         CC
#> 3          AA         AA         AA         AA         AA         NN         AA
#> 4          AA         AA         AA         AA         AA         AA         AA
#> 5          TT         TT         TT         CT         TT         CT         TT
#> 6          AA         AA         AA         AA         AA         AA         AA
#> 7          AA         AA         AA         NN         AA         AA         NN
#> 8          GG         GG         GG         GG         GG         GG         GG
#> 9          CC         CC         CC         CC         CC         CT         CT
#> 10         AA         AA         AA         AA         AA         AA         AA
#> 11         AA         AA         NN         AA         AA         AA         AA
#> 12         GG         GG         GG         GG         GG         GG         GG
#> 13         TT         TT         TT         TT         TT         TT         TT
#> 14         TT         TT         TT         TT         TT         TT         TT
#> 15         GG         GG         GG         GG         GG         GG         GG
#> 16         CC         CC         CC         CC         CC         CC         CC
#> 17         NN         TT         AA         NN         AA         NN         AT
#>    PAL.r2.276 PAL.r2.278 PAL.r2.281 PAL.r2.285 PAL.r2.287 PAL.r2.290 PAL.r2.291
#> 1          TT         TT         TT         TT         TT         TT         TT
#> 2          CC         CC         CC         CC         CC         CC         CC
#> 3          AA         AA         AA         AA         AA         AA         AA
#> 4          AA         AA         AA         AA         AA         AA         AA
#> 5          TT         TT         NN         NN         TT         CT         TT
#> 6          AA         AA         AA         AA         AA         NN         AA
#> 7          NN         AA         NN         AA         AA         AA         AA
#> 8          GG         GG         GG         GG         GG         NN         GG
#> 9          CC         CT         CT         TT         CT         CC         CT
#> 10         AA         AA         AA         AA         AA         AA         AA
#> 11         AA         AA         AA         AA         AA         AA         AA
#> 12         GG         GG         GG         GG         GG         GG         GG
#> 13         TT         TT         TT         TT         NN         TT         TT
#> 14         TT         TT         TT         TT         TT         TT         TT
#> 15         GG         GG         GG         GG         GG         NN         NN
#> 16         CC         CC         CC         CC         CC         CC         CC
#> 17         AT         AT         AA         AT         AA         NN         AA

# bracket operator, excluding first ten SNPs and only keeping pop = PAL,
# fam = B
stickSNPs[-c(1:10), pop.fam = c("PAL.B")]
#>    PAL.r2.237 PAL.r2.241 PAL.r2.275 PAL.r2.285
#> 1          TT         TT         TT         TT
#> 2          AG         AG         AG         AA
#> 3          TT         NN         NN         TT
#> 4          GG         GG         GG         GG
#> 5          GG         GG         GG         GG
#> 6          NN         GT         GG         GT
#> 7          AA         AC         AC         AA
#> 8          AA         AC         CC         AA
#> 9          TT         TT         TT         TT
#> 10         GG         GG         GG         GG
#> 11         CC         CC         CC         CC
#> 12         CC         CC         CC         CC
#> 13         GG         GG         GG         GG
#> 14         AA         AA         AA         AA
#> 15         GG         GG         GG         GG
#> 16         CT         CT         TT         CT
#> 17         AA         AA         AA         AA
#> 18         TT         TT         TT         NN
#> 19         AA         AA         AC         NN
#> 20         GG         GG         GG         GG
#> 21         CC         CC         CC         CC
#> 22         AC         AC         AC         NN
#> 23         CT         CT         CC         CC
#> 24         CC         CC         CC         CC
#> 25         AA         AA         AA         AA
#> 26         AA         AA         AA         AA
#> 27         GG         GG         GG         GG
#> 28         GT         GT         GT         GG
#> 29         AT         TT         TT         NN
#> 30         CG         CG         CG         GG
#> 31         AA         AA         NN         AA
#> 32         GG         GG         GG         GG
#> 33         AA         AA         AA         AA
#> 34         AT         AT         AA         NN
#> 35         GG         GG         GG         GG
#> 36         CC         CC         CC         CC
#> 37         CC         NN         CG         CC
#> 38         CC         CT         TT         TT
#> 39         GG         GG         GG         GG
#> 40         CT         CC         CC         CC
#> 41         TT         NN         CT         TT
#> 42         AG         AA         NN         AG
#> 43         AA         AA         AA         NN
#> 44         GG         GG         GG         GG
#> 45         AG         AG         AG         NN
#> 46         AG         AG         AG         AA
#> 47         NN         GG         GG         GG
#> 48         TT         CT         NN         NN
#> 49         GG         GG         GG         AG
#> 50         AA         AA         AA         AA
#> 51         TT         TT         TT         TT
#> 52         GG         GG         GG         GG
#> 53         NN         CC         AC         NN
#> 54         AA         AT         TT         TT
#> 55         CC         CC         CT         TT
#> 56         AC         CC         CC         AC
#> 57         AA         AG         AG         GG
#> 58         TT         AT         TT         TT
#> 59         AA         AA         AA         AA
#> 60         AA         NN         AA         AA
#> 61         AT         AA         AA         AA
#> 62         AG         AG         GG         AG
#> 63         TT         GT         GT         TT
#> 64         GG         GG         GG         GG
#> 65         AC         AC         AC         AC
#> 66         NN         TT         NN         TT
#> 67         AT         AT         AT         NN
#> 68         GG         GG         GG         GG
#> 69         TT         TT         TT         CT
#> 70         CC         CT         CT         CT
#> 71         NN         AA         AT         AT
#> 72         TT         TT         TT         TT
#> 73         TT         TT         TT         TT
#> 74         NN         TT         TT         TT
#> 75         AA         AA         AA         AA
#> 76         NN         AG         AG         NN
#> 77         TT         TT         TT         TT
#> 78         NN         GG         GG         GG
#> 79         NN         CC         CC         GG
#> 80         CC         CC         CC         CC
#> 81         CC         CC         CC         CT
#> 82         NN         AA         NN         AG
#> 83         CT         CT         TT         TT
#> 84         CG         CC         CG         NN
#> 85         GG         GG         GG         GG
#> 86         CC         CC         CC         CC
#> 87         GG         GG         GG         GG
#> 88         GG         GT         GT         GT
#> 89         GG         GG         GG         GG
#> 90         NN         AA         AT         AT

# using the .facet etc arguments, useful when the facet is stored as an
# object
target_facet <- "pop"
target_subfacet <- "ASP"
subset_snpR_data(stickSNPs, 
                 .facets = target_facet, 
                 .subfacets = target_subfacet)
#> Warning: Since samples were subset, any stats will need to be recalculated.
#>     ASP.r2.141 ASP.r2.144 ASP.r2.146 ASP.r2.147 ASP.r2.196 ASP.r2.216
#> 1           AG         GG         AA         AA         NN         AG
#> 2           GG         GG         GG         GG         GT         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AA         AA         AA         AA         AG         AA
#> 5           GG         NN         NN         AA         GG         AG
#> 6           AC         AC         AA         AA         AA         AC
#> 7           TT         TT         TT         TT         GT         TT
#> 8           GG         GG         AG         AG         AG         GG
#> 9           CC         AC         CC         AC         CC         CC
#> 10          AA         AA         AG         AA         AA         AA
#> 11          NN         TT         TT         CT         TT         CT
#> 12          AG         AG         AG         NN         AA         AA
#> 13          TT         TT         CT         TT         CT         NN
#> 14          GG         GT         GT         NN         GG         GG
#> 15          AA         AG         AG         AA         AA         AA
#> 16          NN         GG         GT         TT         GG         GT
#> 17          AA         AA         AA         AA         AA         NN
#> 18          AC         NN         CC         CC         AC         AC
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         AG         GG         GG         NN         GG
#> 21          NN         NN         CC         CC         AC         AC
#> 22          CT         NN         TT         CC         CT         NN
#> 23          GG         GG         GG         GG         GG         GG
#> 24          AA         AA         AA         NN         AA         AA
#> 25          NN         AG         AG         GG         GG         AG
#> 26          CT         CT         CT         CT         NN         TT
#> 27          AA         AA         AA         AA         AA         AA
#> 28          TT         NN         NN         NN         CT         CT
#> 29          AA         AC         AA         NN         AA         AA
#> 30          GG         GG         CG         GG         GG         GG
#> 31          CG         CG         CC         CC         CC         CC
#> 32          AA         AA         AA         NN         AA         AA
#> 33          CC         CT         CC         CC         CT         CT
#> 34          CG         CG         CG         CG         CC         NN
#> 35          AA         AA         AA         AC         AA         NN
#> 36          AA         AA         AA         AA         AA         AG
#> 37          CG         GG         GG         GG         GG         GG
#> 38          GG         TT         GG         GT         GT         TT
#> 39          TT         NN         TT         NN         AT         TT
#> 40          CG         CC         CG         NN         CG         GG
#> 41          AG         AA         AG         AG         AA         AA
#> 42          AG         NN         NN         NN         GG         GG
#> 43          AG         AG         AA         AA         AA         AA
#> 44          AT         TT         TT         TT         NN         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          NN         CC         AA         AA         NN         CC
#> 47          CC         CC         CC         CG         CG         CC
#> 48          CT         TT         TT         NN         TT         TT
#> 49          GG         GG         GG         GG         NN         GG
#> 50          NN         CC         CT         CC         CC         CT
#> 51          TT         NN         CT         CT         TT         CT
#> 52          AA         AG         AG         AA         AG         AG
#> 53          AA         AA         CC         NN         AA         AC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          GG         NN         AG         NN         AG         AG
#> 56          AG         AG         AG         AG         AA         NN
#> 57          GG         NN         GG         NN         GG         GG
#> 58          TT         CT         CC         NN         TT         CC
#> 59          AA         AG         GG         AG         GG         AA
#> 60          AA         AA         AC         NN         AC         AC
#> 61          TT         NN         TT         NN         TT         TT
#> 62          CG         GG         GG         NN         GG         GG
#> 63          AC         NN         NN         NN         AA         AC
#> 64          TT         AT         AT         TT         AT         AT
#> 65          CT         CC         CT         CC         CC         CT
#> 66          AC         AC         NN         NN         AC         CC
#> 67          GG         GG         GG         NN         GG         GG
#> 68          AA         AT         AT         AA         AA         NN
#> 69          AA         AA         NN         AA         AA         AA
#> 70          AC         AC         AA         AC         AC         NN
#> 71          AA         AA         AA         AT         AA         AA
#> 72          GG         GG         GG         GG         AG         NN
#> 73          TT         GT         TT         TT         GT         TT
#> 74          GT         GG         GT         GG         GG         GG
#> 75          AC         AC         NN         AC         CC         AC
#> 76          CC         CT         CT         CT         NN         NN
#> 77          AT         AT         AT         AT         AA         AA
#> 78          GG         GG         GG         GG         GG         GG
#> 79          TT         TT         CT         TT         CT         CT
#> 80          CT         CT         TT         TT         TT         CT
#> 81          AA         NN         AT         AT         AT         AT
#> 82          CT         TT         CT         TT         TT         TT
#> 83          AT         TT         AT         TT         TT         AT
#> 84          CC         CC         CT         CT         CC         NN
#> 85          AG         AG         GG         AA         AG         NN
#> 86          AG         AG         AG         AG         AA         AA
#> 87          CT         CC         CT         CT         TT         CT
#> 88          AG         NN         AG         NN         GG         AG
#> 89          NN         CG         CG         CG         CC         CG
#> 90          CT         CT         CC         CC         CC         CT
#> 91          CC         CC         NN         CT         CC         CT
#> 92          AA         NN         GG         NN         GG         AA
#> 93          TT         TT         TT         TT         CT         TT
#> 94          CC         CC         CC         CC         CC         NN
#> 95          GG         GG         GG         NN         GT         GG
#> 96          CC         NN         CC         NN         CC         CC
#> 97          GG         GG         GG         GG         GG         GG
#> 98          GG         GG         GG         NN         GG         GT
#> 99          AA         GG         AA         AA         AG         AG
#> 100         AT         NN         AT         AT         TT         AT
#>     ASP.r2.221 ASP.r2.227 ASP.r2.230 ASP.r2.232 ASP.r2.235
#> 1           NN         GG         AA         GG         AG
#> 2           GG         NN         GG         GG         GG
#> 3           CT         CC         CC         CT         CC
#> 4           AA         AA         AA         AG         AG
#> 5           AG         AG         NN         GG         NN
#> 6           AA         AC         AA         AC         NN
#> 7           TT         TT         TT         TT         TT
#> 8           GG         AG         GG         AG         GG
#> 9           CC         CC         AC         CC         AA
#> 10          AA         AA         AG         AA         AA
#> 11          CC         CC         CT         CT         CT
#> 12          AG         AA         AG         AA         AA
#> 13          TT         TT         NN         TT         TT
#> 14          GG         GG         GT         GG         GT
#> 15          AA         AA         AG         AG         AG
#> 16          GT         GG         GG         GG         GG
#> 17          AC         AC         AA         AC         AA
#> 18          CC         AC         AC         CC         AC
#> 19          TT         TT         TT         TT         TT
#> 20          GG         GG         GG         NN         GG
#> 21          CC         CC         AC         CC         CC
#> 22          CC         TT         CT         CC         CT
#> 23          GG         AG         GG         GG         NN
#> 24          AA         AA         AA         AA         AG
#> 25          AG         GG         AG         AG         GG
#> 26          NN         TT         CT         TT         CT
#> 27          AA         AA         NN         AC         AA
#> 28          TT         TT         TT         TT         TT
#> 29          AA         NN         AA         AA         AA
#> 30          GG         GG         GG         GG         CG
#> 31          CG         CC         CG         CG         CG
#> 32          AA         AC         AC         AC         AA
#> 33          CT         CC         CC         CT         CT
#> 34          CC         CG         CG         CG         NN
#> 35          AC         AA         AA         AA         AA
#> 36          AA         AG         AA         AG         NN
#> 37          NN         GG         CC         GG         CG
#> 38          GG         GT         NN         GG         TT
#> 39          TT         TT         TT         AT         AT
#> 40          CG         GG         CC         GG         CG
#> 41          AG         AG         GG         AG         AG
#> 42          AG         AG         AG         GG         GG
#> 43          AA         AA         AA         AA         AA
#> 44          TT         NN         TT         AA         TT
#> 45          GG         GG         GG         GG         GG
#> 46          NN         AA         CC         CC         AC
#> 47          CG         CG         CC         CG         CC
#> 48          CT         CC         NN         TT         CT
#> 49          GG         GG         GG         GG         GG
#> 50          CC         CC         CT         CC         CC
#> 51          CT         NN         CT         CC         TT
#> 52          AA         AG         AG         AG         AG
#> 53          AA         AA         AC         AA         AA
#> 54          GG         GG         GG         GG         GG
#> 55          AG         AG         AG         AA         NN
#> 56          AG         AA         AG         AA         NN
#> 57          GG         AG         NN         GG         GG
#> 58          NN         CT         CC         CT         CT
#> 59          GG         AG         GG         AA         AG
#> 60          AC         NN         NN         CC         NN
#> 61          TT         NN         TT         TT         TT
#> 62          GG         GG         CG         GG         GG
#> 63          AC         AC         CC         AA         AA
#> 64          AA         AT         TT         TT         TT
#> 65          CT         CC         CC         CC         CC
#> 66          CC         AC         AA         AA         AC
#> 67          GG         GG         NN         GG         GG
#> 68          AT         TT         AA         AA         AA
#> 69          AA         AA         AA         AA         AA
#> 70          AC         AC         AA         AA         AC
#> 71          AA         AA         AA         AA         AA
#> 72          AG         AG         AA         GG         AA
#> 73          GT         NN         GT         GT         GT
#> 74          GT         GG         GG         GG         NN
#> 75          AC         AC         AA         CC         NN
#> 76          TT         NN         NN         CC         NN
#> 77          AT         NN         AT         NN         AA
#> 78          GG         GG         NN         GG         GG
#> 79          CT         TT         TT         TT         TT
#> 80          TT         TT         CC         TT         TT
#> 81          AT         AA         AT         NN         AA
#> 82          NN         TT         TT         TT         TT
#> 83          TT         TT         TT         AT         TT
#> 84          CT         CC         CT         CC         TT
#> 85          AA         AA         AG         AA         AG
#> 86          AG         NN         AG         NN         NN
#> 87          CC         CT         TT         TT         CC
#> 88          AA         AG         AG         AG         NN
#> 89          CG         NN         CC         CG         CC
#> 90          CC         CC         CC         CC         CC
#> 91          CC         CT         CC         NN         CC
#> 92          AG         AG         NN         NN         AG
#> 93          TT         TT         TT         CT         CT
#> 94          CC         CC         CC         NN         NN
#> 95          GT         GT         GG         GT         GG
#> 96          CC         CC         AC         CC         NN
#> 97          GG         GG         GG         GG         GG
#> 98          GG         GT         NN         GG         TT
#> 99          AG         AG         NN         GG         AG
#> 100         TT         TT         TT         AA         AA
```
