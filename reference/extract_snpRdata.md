# Get from or overwrite components of a snpRdata object

Fetch or overwrite the major parts of a snpRdata object (genotypes, snp
meta, or sample meta). If overwritten, any calculated stats will be
removed, since their values may be dependant upon changes in metadata.

## Usage

``` r
# S4 method for class 'snpRdata'
genotypes(x)

# S4 method for class 'snpRdata'
genotypes(x) <- value

# S4 method for class 'snpRdata'
snp.meta(x, value)

# S4 method for class 'snpRdata'
snp.meta(x) <- value

# S4 method for class 'snpRdata'
sample.meta(x, value)

# S4 method for class 'snpRdata'
sample.meta(x) <- value
```

## Arguments

- x:

  snpRdata object to get genotype data from.

- value:

  Genotypes, snp metadata, or sample metadata

## Functions

- `genotypes(snpRdata)`: view genotypes

- `genotypes(snpRdata) <- value`: set genotypes

- `snp.meta(snpRdata)`: view snp meta

- `snp.meta(snpRdata) <- value`: set snp meta

- `sample.meta(snpRdata)`: view sample meta

- `sample.meta(snpRdata) <- value`: set sample meta

## Examples

``` r
# copy test data
test <- stickSNPs

# show genotypes
genotypes(test)
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
#>     ASP.r2.221 ASP.r2.227 ASP.r2.230 ASP.r2.232 ASP.r2.235 CLF.r2.040
#> 1           NN         GG         AA         GG         AG         NN
#> 2           GG         NN         GG         GG         GG         GG
#> 3           CT         CC         CC         CT         CC         CC
#> 4           AA         AA         AA         AG         AG         AA
#> 5           AG         AG         NN         GG         NN         GG
#> 6           AA         AC         AA         AC         NN         AA
#> 7           TT         TT         TT         TT         TT         GG
#> 8           GG         AG         GG         AG         GG         GG
#> 9           CC         CC         AC         CC         AA         CC
#> 10          AA         AA         AG         AA         AA         AA
#> 11          CC         CC         CT         CT         CT         CT
#> 12          AG         AA         AG         AA         AA         AA
#> 13          TT         TT         NN         TT         TT         CT
#> 14          GG         GG         GT         GG         GT         GG
#> 15          AA         AA         AG         AG         AG         AA
#> 16          GT         GG         GG         GG         GG         GT
#> 17          AC         AC         AA         AC         AA         AA
#> 18          CC         AC         AC         CC         AC         CC
#> 19          TT         TT         TT         TT         TT         GT
#> 20          GG         GG         GG         NN         GG         GG
#> 21          CC         CC         AC         CC         CC         NN
#> 22          CC         TT         CT         CC         CT         CC
#> 23          GG         AG         GG         GG         NN         GG
#> 24          AA         AA         AA         AA         AG         AA
#> 25          AG         GG         AG         AG         GG         GG
#> 26          NN         TT         CT         TT         CT         TT
#> 27          AA         AA         NN         AC         AA         AA
#> 28          TT         TT         TT         TT         TT         CT
#> 29          AA         NN         AA         AA         AA         AA
#> 30          GG         GG         GG         GG         CG         GG
#> 31          CG         CC         CG         CG         CG         CC
#> 32          AA         AC         AC         AC         AA         AC
#> 33          CT         CC         CC         CT         CT         CT
#> 34          CC         CG         CG         CG         NN         CC
#> 35          AC         AA         AA         AA         AA         AA
#> 36          AA         AG         AA         AG         NN         AG
#> 37          NN         GG         CC         GG         CG         GG
#> 38          GG         GT         NN         GG         TT         GG
#> 39          TT         TT         TT         AT         AT         TT
#> 40          CG         GG         CC         GG         CG         CG
#> 41          AG         AG         GG         AG         AG         AG
#> 42          AG         AG         AG         GG         GG         NN
#> 43          AA         AA         AA         AA         AA         AG
#> 44          TT         NN         TT         AA         TT         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          NN         AA         CC         CC         AC         AC
#> 47          CG         CG         CC         CG         CC         CC
#> 48          CT         CC         NN         TT         CT         CT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CT         CC         CC         CC
#> 51          CT         NN         CT         CC         TT         TT
#> 52          AA         AG         AG         AG         AG         GG
#> 53          AA         AA         AC         AA         AA         AC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          AG         AG         AG         AA         NN         GG
#> 56          AG         AA         AG         AA         NN         AA
#> 57          GG         AG         NN         GG         GG         NN
#> 58          NN         CT         CC         CT         CT         CC
#> 59          GG         AG         GG         AA         AG         AG
#> 60          AC         NN         NN         CC         NN         AC
#> 61          TT         NN         TT         TT         TT         TT
#> 62          GG         GG         CG         GG         GG         GG
#> 63          AC         AC         CC         AA         AA         CC
#> 64          AA         AT         TT         TT         TT         TT
#> 65          CT         CC         CC         CC         CC         TT
#> 66          CC         AC         AA         AA         AC         AC
#> 67          GG         GG         NN         GG         GG         GG
#> 68          AT         TT         AA         AA         AA         AT
#> 69          AA         AA         AA         AA         AA         NN
#> 70          AC         AC         AA         AA         AC         CC
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AG         AG         AA         GG         AA         AG
#> 73          GT         NN         GT         GT         GT         GT
#> 74          GT         GG         GG         GG         NN         GG
#> 75          AC         AC         AA         CC         NN         AA
#> 76          TT         NN         NN         CC         NN         CT
#> 77          AT         NN         AT         NN         AA         AA
#> 78          GG         GG         NN         GG         GG         TT
#> 79          CT         TT         TT         TT         TT         TT
#> 80          TT         TT         CC         TT         TT         CT
#> 81          AT         AA         AT         NN         AA         TT
#> 82          NN         TT         TT         TT         TT         TT
#> 83          TT         TT         TT         AT         TT         AT
#> 84          CT         CC         CT         CC         TT         TT
#> 85          AA         AA         AG         AA         AG         AG
#> 86          AG         NN         AG         NN         NN         AA
#> 87          CC         CT         TT         TT         CC         TT
#> 88          AA         AG         AG         AG         NN         AA
#> 89          CG         NN         CC         CG         CC         CC
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CC         CT         CC         NN         CC         CC
#> 92          AG         AG         NN         NN         AG         AG
#> 93          TT         TT         TT         CT         CT         CT
#> 94          CC         CC         CC         NN         NN         CC
#> 95          GT         GT         GG         GT         GG         GT
#> 96          CC         CC         AC         CC         NN         CC
#> 97          GG         GG         GG         GG         GG         GG
#> 98          GG         GT         NN         GG         TT         GG
#> 99          AG         AG         NN         GG         AG         AG
#> 100         TT         TT         TT         AA         AA         AT
#>     CLF.r2.048 CLF.r2.061 CLF.r2.064 CLF.r2.068 CLF.r2.069 OPL.r1.078
#> 1           AG         AG         AA         AG         AA         AG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         CC         CC         CT         CC
#> 4           AG         AA         NN         AA         AA         AG
#> 5           NN         NN         NN         GG         GG         AA
#> 6           AC         AC         AA         AC         AA         AA
#> 7           TT         TT         TT         TT         TT         TT
#> 8           AG         AG         AA         AG         GG         GG
#> 9           AC         CC         AC         CC         CC         CC
#> 10          AA         AA         AA         AG         AA         AA
#> 11          CT         CT         CT         TT         CC         CC
#> 12          AA         AG         AA         AA         AA         AA
#> 13          CT         TT         CT         TT         CC         NN
#> 14          GG         GG         GT         GT         GG         GG
#> 15          AA         AA         AA         AA         AA         NN
#> 16          GT         GG         NN         GG         GT         NN
#> 17          AA         AA         AA         AA         AA         AC
#> 18          AC         AC         NN         AC         AC         CC
#> 19          TT         TT         GT         TT         GT         NN
#> 20          GG         GG         GG         GG         AG         GG
#> 21          AC         AC         NN         AC         CC         NN
#> 22          CC         NN         TT         CT         CT         CC
#> 23          GG         GG         GG         GG         GG         GG
#> 24          AA         AA         AA         AA         AA         AA
#> 25          GG         AG         NN         GG         GG         AG
#> 26          TT         TT         TT         TT         TT         CT
#> 27          AC         AA         AA         AA         AA         AC
#> 28          TT         CT         TT         TT         TT         NN
#> 29          NN         NN         AA         AC         AA         AA
#> 30          GG         CG         GG         CG         NN         GG
#> 31          CG         CC         CC         CC         CC         CG
#> 32          AC         NN         AA         AA         AA         AA
#> 33          CC         CT         NN         CT         NN         NN
#> 34          CG         GG         CG         GG         CC         NN
#> 35          NN         AA         AA         AA         AA         AA
#> 36          AA         AA         AG         AA         AA         NN
#> 37          CG         CG         NN         GG         GG         GG
#> 38          GG         GG         GG         GG         GG         GT
#> 39          NN         NN         AT         NN         AT         NN
#> 40          GG         NN         CG         GG         GG         CG
#> 41          GG         GG         AG         NN         AG         AG
#> 42          GG         AG         GG         AG         AA         GG
#> 43          NN         AA         AA         AG         AA         NN
#> 44          TT         TT         TT         TT         TT         TT
#> 45          GG         GG         GG         GG         GG         NN
#> 46          CC         CC         NN         CC         AC         AC
#> 47          CC         CC         CC         CC         CC         NN
#> 48          CT         TT         CT         CT         CC         CT
#> 49          GG         CG         GG         CG         GG         GG
#> 50          CC         CC         CC         CT         NN         CT
#> 51          NN         TT         TT         TT         TT         CT
#> 52          AG         AA         AA         AA         AG         AG
#> 53          NN         AA         NN         AC         AC         NN
#> 54          GG         GG         GT         GG         GG         NN
#> 55          GG         GG         AG         AG         AG         AA
#> 56          AA         AG         AA         AA         AA         NN
#> 57          GG         GG         GG         GG         GG         AG
#> 58          CC         CC         CC         CC         CC         CT
#> 59          AG         AG         AG         GG         GG         AG
#> 60          AC         CC         CC         CC         AC         AC
#> 61          NN         NN         CT         TT         NN         TT
#> 62          GG         GG         GG         GG         GG         GG
#> 63          AC         AC         AA         NN         AA         NN
#> 64          TT         TT         TT         TT         AA         AT
#> 65          CC         CT         CT         CT         CC         NN
#> 66          AC         CC         AA         AC         AA         NN
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AA         AA         AT         AA         AT         AT
#> 69          AA         AA         AA         AG         AA         GG
#> 70          AA         AC         NN         AC         CC         NN
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AA         AG         GG         GG         AG         NN
#> 73          TT         GT         TT         GT         TT         NN
#> 74          GT         GG         GG         NN         GG         GT
#> 75          AC         AC         CC         AC         NN         AC
#> 76          NN         NN         NN         CC         CC         NN
#> 77          AA         AT         AA         AA         AA         AA
#> 78          GG         GG         GG         GT         GG         NN
#> 79          TT         TT         TT         TT         TT         TT
#> 80          CT         TT         TT         TT         TT         CT
#> 81          NN         AT         AT         AT         AT         AA
#> 82          NN         NN         TT         NN         TT         TT
#> 83          AT         TT         AT         AT         AT         TT
#> 84          CT         CC         CC         CT         CC         CC
#> 85          AA         NN         AG         AG         AG         NN
#> 86          AA         AG         AA         AA         AA         NN
#> 87          CT         TT         TT         TT         TT         CT
#> 88          AA         AG         AA         AG         AG         GG
#> 89          NN         GG         CC         NN         NN         CG
#> 90          CC         CC         CC         CC         CC         CC
#> 91          NN         TT         CC         NN         CT         CT
#> 92          NN         GG         GG         GG         AG         AG
#> 93          TT         CT         TT         CT         TT         TT
#> 94          CC         CC         CC         NN         CC         CG
#> 95          TT         GG         GG         GT         GG         GG
#> 96          CC         CC         AC         CC         CC         AC
#> 97          GG         GG         GG         AG         GG         AG
#> 98          GG         GG         NN         GT         TT         GG
#> 99          AA         GG         AA         GG         GG         NN
#> 100         AT         NN         AA         NN         AT         NN
#>     OPL.r1.080 OPL.r1.083 OPL.r1.086 OPL.r1.087 OPL.r1.088 OPL.r1.109
#> 1           AG         AG         AG         NN         NN         NN
#> 2           NN         GT         GG         GG         GG         GT
#> 3           CT         CT         CC         CC         CT         CC
#> 4           AG         NN         AG         NN         NN         NN
#> 5           NN         NN         AG         AG         NN         GG
#> 6           AA         AC         AA         AA         AA         NN
#> 7           NN         GT         TT         GT         TT         NN
#> 8           GG         GG         GG         AG         GG         AG
#> 9           AC         CC         CC         AC         NN         CC
#> 10          AA         AA         AA         AA         AA         AA
#> 11          NN         NN         CC         CT         CC         CC
#> 12          NN         NN         AA         AG         AG         AG
#> 13          TT         NN         TT         TT         TT         CT
#> 14          GG         GG         GG         GG         GG         GG
#> 15          AG         AA         AA         AA         AA         NN
#> 16          GT         NN         NN         NN         GT         NN
#> 17          AA         AC         NN         AC         NN         AA
#> 18          CC         CC         CC         CC         CC         CC
#> 19          TT         TT         GT         NN         TT         NN
#> 20          GG         NN         GG         GG         GG         GG
#> 21          AC         AA         NN         NN         AC         NN
#> 22          CT         CT         NN         CC         CC         NN
#> 23          GG         AG         GG         GG         GG         GG
#> 24          AA         AA         AG         AA         AA         AA
#> 25          NN         AG         GG         GG         GG         AG
#> 26          TT         TT         NN         NN         TT         TT
#> 27          AC         AA         AC         AA         AA         NN
#> 28          CT         CT         NN         NN         NN         CT
#> 29          AA         NN         AA         NN         AC         NN
#> 30          NN         GG         NN         GG         CG         NN
#> 31          CC         CG         CC         CG         CC         CC
#> 32          AA         AC         NN         NN         AA         AC
#> 33          NN         CC         CT         NN         NN         TT
#> 34          NN         CG         NN         CG         NN         NN
#> 35          AA         AA         NN         NN         AA         AA
#> 36          AG         AA         AA         AA         NN         AG
#> 37          CG         GG         CG         CG         CG         CG
#> 38          GG         NN         NN         NN         GT         GT
#> 39          NN         NN         TT         NN         NN         AT
#> 40          NN         GG         GG         NN         NN         CG
#> 41          NN         NN         NN         AG         NN         AG
#> 42          AG         AG         AG         AG         GG         AG
#> 43          AA         AA         NN         NN         AG         AG
#> 44          TT         TT         NN         TT         TT         TT
#> 45          GG         NN         GG         NN         GG         GG
#> 46          AC         NN         NN         NN         NN         AA
#> 47          CG         CC         CG         CC         CG         CG
#> 48          CC         NN         CT         NN         NN         NN
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CT         NN         NN         CC
#> 51          CC         NN         TT         NN         NN         CT
#> 52          AG         AA         NN         AG         AG         AA
#> 53          AA         AA         AC         AC         NN         AC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          NN         NN         NN         AG         GG         AG
#> 56          AA         NN         NN         AG         AG         AG
#> 57          NN         NN         NN         GG         GG         GG
#> 58          CT         CC         CC         CC         CT         CT
#> 59          AG         AG         AG         AG         AA         GG
#> 60          AC         AC         NN         NN         AC         AC
#> 61          NN         NN         TT         TT         NN         TT
#> 62          GG         CG         GG         CG         GG         GG
#> 63          AC         NN         AA         AC         AA         AA
#> 64          NN         AT         NN         AT         TT         AT
#> 65          NN         TT         NN         CT         TT         CT
#> 66          AC         AC         AC         AA         AC         NN
#> 67          NN         GG         GG         GG         GG         GG
#> 68          AT         AA         AT         NN         AA         NN
#> 69          AA         AA         NN         AG         AG         AA
#> 70          AA         NN         NN         AC         AC         AC
#> 71          AT         AA         AA         AA         AA         AA
#> 72          GG         AG         AG         AG         NN         AA
#> 73          TT         TT         TT         TT         TT         TT
#> 74          GT         GG         NN         GT         NN         GG
#> 75          AC         AC         AC         NN         AC         CC
#> 76          CT         NN         CC         NN         NN         NN
#> 77          AA         AA         NN         AT         AT         AT
#> 78          GG         NN         NN         GT         GT         NN
#> 79          TT         TT         TT         TT         TT         CT
#> 80          CT         CT         NN         NN         CT         TT
#> 81          NN         AA         NN         NN         NN         AT
#> 82          TT         TT         TT         CT         NN         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          CT         NN         NN         CC         CC         NN
#> 85          GG         NN         AG         NN         AG         GG
#> 86          AA         AA         NN         AG         AG         AG
#> 87          NN         TT         CT         CT         TT         TT
#> 88          NN         NN         NN         AG         AG         AA
#> 89          NN         CC         CG         NN         CG         NN
#> 90          CC         CC         NN         CC         CC         CC
#> 91          NN         CT         CT         CT         NN         NN
#> 92          NN         AG         AG         NN         AG         NN
#> 93          TT         TT         TT         CT         NN         TT
#> 94          CC         CC         NN         CC         CC         CG
#> 95          NN         NN         GG         GT         NN         NN
#> 96          CC         AC         NN         NN         NN         CC
#> 97          GG         GG         NN         NN         GG         NN
#> 98          GG         GT         GG         GG         GG         GT
#> 99          NN         AG         AA         AG         AG         AG
#> 100         NN         AT         AT         AT         AT         NN
#>     OPL.r1.110 OPL.r1.116 OPL.r1.125 OPL.r1.126 OPL.r1.128 OPL.r1.130
#> 1           AG         AG         NN         AG         AG         AA
#> 2           GG         GG         GG         GG         GG         GT
#> 3           NN         CT         NN         CC         CC         CC
#> 4           AG         NN         AG         AG         AG         AA
#> 5           GG         AG         GG         NN         GG         NN
#> 6           AA         AA         AA         AA         AA         AA
#> 7           NN         TT         TT         TT         TT         TT
#> 8           NN         GG         AG         NN         AG         GG
#> 9           CC         CC         CC         NN         CC         CC
#> 10          AA         NN         NN         AA         AA         AA
#> 11          NN         NN         CC         NN         CT         TT
#> 12          NN         NN         AA         AA         AA         NN
#> 13          NN         TT         TT         TT         TT         CC
#> 14          TT         GG         NN         GG         GG         GG
#> 15          AA         AG         NN         NN         AA         AG
#> 16          GT         GG         GT         TT         GG         GG
#> 17          AC         AC         AA         AA         AC         AA
#> 18          AC         CC         CC         AC         AC         AA
#> 19          TT         NN         TT         TT         NN         NN
#> 20          GG         NN         GG         NN         AG         GG
#> 21          AA         AC         AC         AC         AC         AC
#> 22          CC         CT         NN         NN         CC         CT
#> 23          GG         AG         GG         GG         GG         GG
#> 24          AA         AA         NN         AA         AG         AA
#> 25          AG         AG         NN         AG         AG         GG
#> 26          TT         NN         TT         NN         TT         TT
#> 27          AC         AA         NN         AC         AA         AA
#> 28          CT         CC         TT         CT         NN         CC
#> 29          NN         AA         NN         AA         AA         AC
#> 30          CG         GG         CG         NN         GG         GG
#> 31          NN         CG         NN         CG         CG         CC
#> 32          AA         AC         AC         AC         AC         NN
#> 33          TT         CC         CC         CT         CT         CT
#> 34          CG         CG         NN         CG         CG         CC
#> 35          AC         AA         AA         AA         AC         AA
#> 36          AA         NN         AA         AA         AA         AA
#> 37          GG         CG         NN         GG         CG         GG
#> 38          NN         GG         GT         GT         GG         GG
#> 39          NN         NN         NN         AT         TT         TT
#> 40          NN         CG         CG         GG         CG         NN
#> 41          NN         NN         GG         GG         AG         GG
#> 42          AA         NN         GG         AG         GG         AA
#> 43          NN         NN         AG         AA         AA         AA
#> 44          AT         TT         TT         TT         TT         TT
#> 45          GG         GG         GG         GG         GG         NN
#> 46          CC         CC         NN         AC         CC         AC
#> 47          CG         CG         CC         CG         CG         CC
#> 48          NN         CC         TT         CT         CT         CT
#> 49          GG         CG         GG         GG         GG         GG
#> 50          CC         CC         CC         NN         CC         CC
#> 51          CT         CT         CT         TT         TT         TT
#> 52          AA         AG         AA         AG         AG         AA
#> 53          AA         NN         AA         AC         AC         CC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          AA         GG         AG         AG         AG         GG
#> 56          AG         AA         AA         AA         AA         NN
#> 57          NN         GG         GG         GG         GG         AG
#> 58          CT         CC         CC         CC         CC         CT
#> 59          GG         AA         AG         NN         AG         AG
#> 60          AC         AA         AA         AC         AA         CC
#> 61          TT         TT         NN         NN         TT         TT
#> 62          GG         CG         NN         GG         CG         CG
#> 63          NN         AC         AA         AC         AC         AA
#> 64          NN         AT         AT         AT         TT         TT
#> 65          NN         NN         CT         CC         TT         CT
#> 66          AA         AA         NN         NN         AA         AA
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AT         AA         AT         AA         AT         AA
#> 69          AA         AA         AA         AA         AA         AG
#> 70          AA         AA         AA         NN         AA         AC
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AG         AA         NN         AA         AA         AG
#> 73          NN         GT         TT         TT         TT         GT
#> 74          NN         GT         NN         GT         GG         GT
#> 75          NN         NN         AC         AA         AC         AA
#> 76          TT         CT         NN         CC         TT         TT
#> 77          NN         AA         AA         AA         AA         AA
#> 78          NN         GG         GG         GG         GG         GT
#> 79          NN         CT         TT         TT         TT         TT
#> 80          CT         NN         TT         CT         TT         TT
#> 81          AT         AA         AT         AT         AT         TT
#> 82          TT         CT         TT         TT         TT         CT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          CT         NN         CT         NN         CC         CT
#> 85          AG         AG         AG         GG         AG         NN
#> 86          AG         NN         AA         NN         AA         AA
#> 87          TT         NN         TT         TT         TT         TT
#> 88          AG         NN         AA         AA         AA         NN
#> 89          NN         NN         NN         CG         CG         GG
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CT         CC         CT         CT         NN         CC
#> 92          GG         NN         NN         NN         GG         NN
#> 93          TT         CT         TT         TT         TT         TT
#> 94          CC         CC         CC         CC         CC         CC
#> 95          GT         NN         GG         GG         GG         GG
#> 96          CC         CC         CC         CC         AC         CC
#> 97          GG         NN         GG         GG         GG         GG
#> 98          GT         GG         NN         GG         GG         GG
#> 99          NN         NN         AG         AA         AA         AG
#> 100         NN         NN         AA         NN         AT         AA
#>     OPL.r1.139 OPL.r1.146 OPL.r1.147 OPL.r1.149 OPL.r2.007 OPL.r2.012
#> 1           AA         NN         AG         GG         AA         AG
#> 2           GG         NN         GG         NN         GG         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AG         AG         AG         AG         AG         AG
#> 5           NN         AG         NN         NN         GG         AG
#> 6           AA         AA         AC         AA         AA         AA
#> 7           NN         TT         TT         NN         GT         TT
#> 8           GG         AG         NN         AG         GG         AA
#> 9           CC         CC         CC         CC         AC         AC
#> 10          AA         AA         AA         AA         AA         AA
#> 11          CT         CT         CT         CC         TT         TT
#> 12          AG         AG         AA         AG         AG         AA
#> 13          TT         CT         CT         CC         TT         CT
#> 14          NN         GG         GG         GG         GG         GG
#> 15          AA         AG         AA         NN         NN         AA
#> 16          GG         GT         GT         GG         GT         GT
#> 17          AA         AA         AA         AA         AA         AA
#> 18          AC         CC         CC         CC         CC         CC
#> 19          GT         TT         TT         GT         TT         TT
#> 20          GG         GG         GG         GG         GG         GG
#> 21          NN         CC         CC         AC         AC         NN
#> 22          CC         NN         TT         CT         CC         CT
#> 23          GG         GG         GG         GG         GG         GG
#> 24          AA         AG         AA         AA         AA         AA
#> 25          GG         AG         GG         GG         GG         GG
#> 26          TT         TT         TT         CT         TT         NN
#> 27          AC         AC         AC         AA         NN         AA
#> 28          CC         CC         TT         CT         CC         NN
#> 29          AA         AA         AA         AC         AC         AC
#> 30          GG         GG         GG         GG         GG         GG
#> 31          CG         GG         CC         CC         CG         CG
#> 32          NN         AA         AA         AA         AA         AA
#> 33          CC         CC         NN         CT         NN         NN
#> 34          CC         NN         CC         CC         CG         CG
#> 35          NN         NN         AA         AC         AA         AA
#> 36          AG         AA         AA         AG         AA         AG
#> 37          GG         GG         CG         GG         GG         GG
#> 38          GG         GT         GG         GG         GG         GG
#> 39          TT         NN         TT         NN         AT         TT
#> 40          NN         GG         GG         GG         GG         GG
#> 41          AG         GG         GG         GG         GG         NN
#> 42          NN         AG         AG         AA         AA         AG
#> 43          AA         AA         AA         AA         NN         AA
#> 44          TT         TT         TT         TT         TT         TT
#> 45          GG         GG         NN         GG         GG         NN
#> 46          NN         NN         AC         CC         NN         AC
#> 47          CG         NN         GG         CG         NN         CG
#> 48          TT         CT         CT         CC         CC         CT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CC         CT         CC         CC
#> 51          NN         CT         CT         TT         CC         NN
#> 52          NN         AG         AA         AA         AG         AG
#> 53          AC         AA         AC         AC         AC         AA
#> 54          GG         GG         GG         GG         GG         GG
#> 55          NN         GG         GG         GG         AG         NN
#> 56          AG         AG         AA         AG         AG         AA
#> 57          AG         AG         GG         GG         GG         GG
#> 58          CT         CT         NN         CC         CC         CT
#> 59          NN         AG         NN         AG         AG         GG
#> 60          AA         NN         AC         AC         AC         AA
#> 61          CT         TT         TT         TT         NN         TT
#> 62          GG         CG         CG         GG         CG         CG
#> 63          NN         AA         AC         NN         AA         NN
#> 64          TT         NN         AT         TT         AT         NN
#> 65          TT         NN         CT         TT         TT         NN
#> 66          AC         AA         AC         AC         AA         AC
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AA         AT         AT         AA         AA         AA
#> 69          AG         AG         AA         AG         GG         AA
#> 70          AC         AC         AA         AC         NN         AA
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AA         AG         AA         AG         AA         AA
#> 73          TT         GT         TT         NN         TT         GT
#> 74          GG         NN         GT         GG         GT         NN
#> 75          AA         AC         AA         NN         AA         AA
#> 76          NN         CT         NN         CC         NN         NN
#> 77          AT         AT         AA         AT         AT         NN
#> 78          GT         GT         GG         GT         TT         GG
#> 79          TT         TT         TT         TT         TT         TT
#> 80          NN         NN         TT         CC         TT         TT
#> 81          AT         AA         AT         NN         AT         AT
#> 82          TT         TT         TT         TT         CT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          CC         NN         TT         CT         CC         CT
#> 85          AG         AG         AG         GG         AG         AG
#> 86          AG         AG         AA         AG         AG         AA
#> 87          TT         CT         CT         TT         CT         TT
#> 88          AG         AA         AA         AA         AG         NN
#> 89          CC         NN         CC         CG         CC         CC
#> 90          CC         CC         CC         CC         CC         CC
#> 91          TT         TT         CT         NN         TT         CT
#> 92          NN         AG         AG         AG         GG         NN
#> 93          TT         TT         TT         TT         TT         TT
#> 94          CC         NN         CC         CC         NN         CC
#> 95          GG         GG         GG         GG         TT         GG
#> 96          CC         CC         CC         CC         AC         CC
#> 97          GG         GG         AG         NN         GG         GG
#> 98          GT         GG         GG         GG         GG         GG
#> 99          GG         NN         AA         AG         AG         AG
#> 100         NN         TT         AA         AT         AA         AT
#>     OPL.r2.022 OPL.r2.023 OPL.r2.028 OPL.r2.031 OPL.r2.033 OPL.r2.079
#> 1           AA         NN         NN         AG         AG         AA
#> 2           GG         NN         GG         GG         GG         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           NN         AA         AG         AG         NN         AG
#> 5           GG         NN         NN         NN         NN         NN
#> 6           AA         AC         AA         AA         AA         AA
#> 7           GT         GT         TT         GT         TT         TT
#> 8           GG         GG         GG         GG         GG         GG
#> 9           CC         CC         CC         CC         AC         CC
#> 10          AA         AA         AG         AA         AA         AA
#> 11          CT         CT         CC         CT         CT         NN
#> 12          AA         NN         NN         AA         NN         AA
#> 13          CT         TT         CT         TT         TT         CT
#> 14          GG         GG         GG         GG         GG         GG
#> 15          AA         AG         GG         NN         AA         AG
#> 16          GG         GT         NN         GG         GT         GT
#> 17          AA         AA         AA         AC         AA         AA
#> 18          AC         CC         CC         AC         CC         AC
#> 19          TT         TT         NN         TT         TT         TT
#> 20          GG         GG         GG         GG         GG         GG
#> 21          CC         AC         AC         NN         AC         NN
#> 22          CT         TT         CC         CT         CC         CT
#> 23          GG         GG         GG         GG         GG         NN
#> 24          AA         AA         AA         AA         AA         AA
#> 25          GG         AG         NN         AG         GG         NN
#> 26          TT         CT         CT         TT         TT         CT
#> 27          AA         AC         AA         AA         AA         AC
#> 28          CT         NN         CC         CT         NN         CC
#> 29          AA         AA         CC         CC         NN         AA
#> 30          GG         GG         GG         GG         NN         NN
#> 31          CG         CC         CC         CC         CC         CG
#> 32          AA         AA         AA         AA         AA         AA
#> 33          CC         TT         CC         NN         CT         CC
#> 34          CG         CG         GG         GG         CC         NN
#> 35          AA         AA         AA         NN         AA         AA
#> 36          GG         AG         AA         AA         AA         AA
#> 37          CG         CG         NN         GG         CG         CG
#> 38          GT         GT         GG         GG         GG         GT
#> 39          TT         NN         NN         TT         TT         TT
#> 40          GG         GG         CG         GG         NN         GG
#> 41          AG         AG         GG         AG         AG         NN
#> 42          AA         GG         GG         NN         AG         AG
#> 43          AG         AA         AG         AG         AA         AG
#> 44          TT         TT         AT         TT         TT         TT
#> 45          GG         GG         NN         GG         GG         NN
#> 46          AC         AC         AA         CC         NN         AC
#> 47          GG         CG         GG         CG         GG         CG
#> 48          CT         TT         CC         CC         CT         NN
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CT         CC         CC         NN         CC
#> 51          NN         TT         TT         NN         CT         CT
#> 52          AG         AA         AA         GG         AG         NN
#> 53          AA         NN         NN         AC         NN         AC
#> 54          GT         GG         GG         GG         GG         GG
#> 55          GG         GG         AG         AG         AA         NN
#> 56          AA         AA         AA         AA         NN         AA
#> 57          GG         GG         NN         GG         NN         NN
#> 58          CT         NN         CC         CC         CT         CT
#> 59          AG         AG         AG         GG         GG         GG
#> 60          AA         AC         NN         AC         AA         AA
#> 61          NN         TT         NN         NN         TT         TT
#> 62          GG         GG         GG         CG         CG         GG
#> 63          AA         NN         AC         AA         NN         AC
#> 64          AT         AT         TT         AT         TT         NN
#> 65          CT         CT         NN         TT         CC         CT
#> 66          AC         NN         NN         AC         NN         CC
#> 67          GG         GG         GG         GG         NN         GG
#> 68          AT         NN         TT         AA         AA         AA
#> 69          AG         AG         NN         AG         AA         AA
#> 70          AC         AC         AC         AC         AA         AA
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AG         AG         AG         AA         AA         AG
#> 73          TT         TT         TT         TT         GT         NN
#> 74          NN         GG         NN         GT         GG         GT
#> 75          AC         AC         AC         CC         AC         AC
#> 76          CC         TT         NN         CC         CC         CC
#> 77          AA         AA         AA         AA         AT         NN
#> 78          GT         GT         GG         GT         GG         GG
#> 79          CT         CT         TT         CT         TT         CT
#> 80          TT         CT         TT         CT         CT         CT
#> 81          AT         AA         TT         AT         AA         AT
#> 82          CT         TT         TT         TT         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          TT         CT         CT         CC         NN         CT
#> 85          AG         GG         GG         GG         GG         NN
#> 86          NN         NN         AA         NN         AG         AA
#> 87          TT         CC         TT         TT         CC         TT
#> 88          AG         GG         AA         AG         NN         AG
#> 89          CC         CG         GG         CC         CG         NN
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CT         TT         CT         CC         CT         CT
#> 92          AG         NN         NN         NN         NN         AG
#> 93          TT         TT         CT         TT         TT         TT
#> 94          CC         NN         CC         CC         CC         NN
#> 95          GG         GG         GG         GT         GG         TT
#> 96          AC         AC         AC         CC         NN         AC
#> 97          GG         GG         NN         GG         GG         GG
#> 98          GG         GT         GG         GG         GG         GG
#> 99          AG         AG         AA         GG         NN         AG
#> 100         AA         AA         AA         TT         AA         AA
#>     OPL.r2.080 OPL.r2.094 PAL.r2.203 PAL.r2.237 PAL.r2.238 PAL.r2.239
#> 1           GG         AG         AG         GG         GG         GG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AG         AA         AA         AA         AA         AA
#> 5           AG         NN         AG         NN         GG         GG
#> 6           AA         AA         AA         AC         AC         NN
#> 7           TT         NN         TT         TT         TT         TT
#> 8           GG         GG         NN         GG         GG         GG
#> 9           CC         CC         CC         CC         CC         AC
#> 10          AA         AA         AA         AA         AA         AA
#> 11          CT         CT         TT         TT         TT         TT
#> 12          AA         NN         AG         AG         AG         AA
#> 13          TT         CT         CT         TT         TT         TT
#> 14          GG         TT         GG         GG         GG         GG
#> 15          AA         AA         GG         GG         AG         AG
#> 16          GG         NN         GT         NN         TT         NN
#> 17          AA         AA         CC         AA         AA         AA
#> 18          AC         AC         NN         AA         AC         AA
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         GG         NN         GG         GG         GG
#> 21          AC         AA         CC         CC         CC         CC
#> 22          CC         NN         CC         CC         CC         CC
#> 23          GG         NN         GG         GG         GG         GG
#> 24          GG         AA         AA         AA         AA         AA
#> 25          AG         AG         GG         GG         GG         GG
#> 26          CT         NN         TT         CT         CT         CC
#> 27          AA         AA         AA         AA         AA         AA
#> 28          CT         CC         TT         TT         TT         TT
#> 29          AC         AC         AC         AA         AC         AC
#> 30          GG         GG         GG         GG         GG         GG
#> 31          GG         CG         CC         CC         CC         CC
#> 32          AA         AA         AA         AC         AA         AC
#> 33          CT         CT         CC         CT         CC         NN
#> 34          NN         CC         CC         CC         CC         CC
#> 35          AA         AA         AA         AA         AA         AA
#> 36          AA         AG         AA         AA         AA         AA
#> 37          CG         NN         GG         GG         GG         GG
#> 38          GG         GG         GG         GT         GT         GT
#> 39          TT         TT         TT         AT         TT         NN
#> 40          GG         CG         GG         CG         GG         CG
#> 41          GG         AG         AA         AA         AA         AA
#> 42          NN         GG         AG         GG         GG         GG
#> 43          AA         AG         AA         AA         AA         AA
#> 44          TT         AT         AT         AT         AA         AT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          CC         AA         CC         CC         CC         CC
#> 47          CC         NN         CC         CC         CC         CC
#> 48          CT         NN         CT         CC         CT         CT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CT         CC         CT         CC         CC
#> 51          TT         CT         NN         TT         CT         CT
#> 52          GG         AA         AG         AG         NN         AA
#> 53          AA         AA         AA         AA         AA         AA
#> 54          GG         GG         GG         GG         GG         GG
#> 55          GG         GG         AA         AG         AG         AG
#> 56          NN         NN         AG         AG         AG         AA
#> 57          GG         GG         NN         NN         GG         GG
#> 58          CT         CC         CC         TT         TT         CT
#> 59          GG         AG         GG         GG         GG         GG
#> 60          AC         AC         AA         AA         AA         AA
#> 61          TT         NN         NN         TT         TT         NN
#> 62          GG         GG         GG         GG         GG         GG
#> 63          AA         NN         NN         NN         AC         CC
#> 64          TT         AT         AT         AA         NN         AT
#> 65          TT         CC         CC         CC         CT         CC
#> 66          AC         AC         CC         AC         AC         AA
#> 67          GG         GG         AG         AA         GG         NN
#> 68          AA         AA         TT         TT         NN         NN
#> 69          GG         AA         AA         AA         AA         AA
#> 70          NN         AA         AA         AA         AA         AA
#> 71          AA         AA         AA         AT         AA         AT
#> 72          AG         AG         GG         AG         AG         GG
#> 73          TT         GT         TT         TT         GT         GT
#> 74          GT         GT         GG         GG         GG         GG
#> 75          CC         AC         AA         AC         AC         AA
#> 76          CC         CC         NN         NN         NN         TT
#> 77          AA         AA         AT         AT         AT         AA
#> 78          TT         GG         GG         GG         GG         GG
#> 79          TT         TT         TT         TT         TT         CT
#> 80          CT         CT         TT         CC         CT         NN
#> 81          AT         NN         TT         NN         AT         AT
#> 82          TT         TT         TT         TT         TT         TT
#> 83          TT         NN         TT         TT         TT         TT
#> 84          CC         CT         NN         NN         TT         TT
#> 85          AG         NN         AA         AA         AA         AA
#> 86          AA         AA         AG         NN         AG         NN
#> 87          CT         CT         TT         TT         NN         TT
#> 88          AG         AA         GG         NN         GG         GG
#> 89          CG         NN         CC         NN         CC         NN
#> 90          CC         CC         CC         CC         CC         CC
#> 91          NN         CC         CC         CC         CC         CC
#> 92          AG         GG         AG         NN         GG         AG
#> 93          TT         TT         NN         CT         TT         TT
#> 94          CG         CC         CG         CG         CC         NN
#> 95          GT         GG         GG         GG         GG         GG
#> 96          CC         CC         CC         CC         CC         CC
#> 97          GG         GG         GG         GG         GG         GG
#> 98          GG         GG         GG         GG         GG         GT
#> 99          AG         AG         GG         GG         GG         GG
#> 100         AT         AA         NN         NN         AT         NN
#>     PAL.r2.240 PAL.r2.241 PAL.r2.242 PAL.r2.249 PAL.r2.253 PAL.r2.275
#> 1           GG         GG         AG         GG         GG         GG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AA         AA         AA         NN         AA         AA
#> 5           GG         GG         AG         NN         AG         NN
#> 6           AC         AC         AC         AC         AA         AC
#> 7           TT         TT         TT         TT         TT         TT
#> 8           GG         GG         GG         GG         GG         NN
#> 9           CC         CC         CC         CC         CC         CC
#> 10          AA         AA         AA         AA         AA         AA
#> 11          TT         TT         TT         NN         TT         TT
#> 12          AG         AG         AG         NN         AA         AG
#> 13          TT         NN         TT         TT         TT         NN
#> 14          GG         GG         GG         GG         GG         GG
#> 15          GG         GG         GG         GG         NN         GG
#> 16          GT         GT         GT         GG         TT         GG
#> 17          AA         AC         AA         AA         NN         AC
#> 18          AA         AC         AA         AC         AC         CC
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         GG         GG         GG         GG         GG
#> 21          CC         CC         CC         NN         CC         CC
#> 22          CC         CC         CC         CC         CC         CC
#> 23          GG         GG         NN         GG         GG         GG
#> 24          AA         AA         AA         AA         NN         AA
#> 25          GG         GG         GG         GG         GG         GG
#> 26          CC         CT         TT         CT         CC         TT
#> 27          AA         AA         AA         AA         AA         AA
#> 28          TT         TT         CT         TT         CT         TT
#> 29          AA         AA         AA         AA         AA         AC
#> 30          GG         GG         GG         GG         GG         GG
#> 31          CC         CC         CC         CC         CC         CC
#> 32          AC         AC         AC         AC         AA         AC
#> 33          CC         CT         CC         CT         CC         CC
#> 34          CC         CC         CC         CC         CC         CC
#> 35          AA         AA         AA         AA         AA         AA
#> 36          AA         AA         AA         AA         AA         AA
#> 37          NN         GG         GG         CG         GG         GG
#> 38          GT         GT         GT         GG         GT         GT
#> 39          NN         TT         AT         NN         NN         TT
#> 40          GG         CG         GG         NN         GG         CG
#> 41          AA         AA         NN         AA         AA         NN
#> 42          GG         GG         GG         NN         GG         GG
#> 43          AA         AA         AA         AA         AA         AA
#> 44          AT         AT         AT         TT         AT         AA
#> 45          GG         GG         GG         GG         GG         GG
#> 46          CC         CC         CC         CC         CC         CC
#> 47          CC         NN         CG         CG         NN         CG
#> 48          TT         CT         TT         NN         CC         TT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CC         CC         CC         CC
#> 51          NN         NN         TT         TT         TT         CT
#> 52          GG         AA         AG         NN         AG         NN
#> 53          AA         AA         AA         AA         AA         AA
#> 54          GG         GG         GG         GG         GG         GG
#> 55          NN         AG         GG         AG         GG         AG
#> 56          AG         AG         AG         AG         AA         AG
#> 57          GG         GG         NN         GG         GG         GG
#> 58          NN         CT         CT         CT         CT         NN
#> 59          GG         GG         GG         GG         GG         GG
#> 60          AA         AA         AA         AA         AA         AA
#> 61          NN         TT         TT         TT         NN         TT
#> 62          GG         GG         GG         NN         GG         GG
#> 63          AC         CC         NN         NN         NN         AC
#> 64          AA         AT         AA         AT         AA         TT
#> 65          CC         CC         CC         CC         CT         CT
#> 66          AA         CC         AC         CC         CC         CC
#> 67          GG         AG         GG         AG         GG         AG
#> 68          TT         AT         AT         AA         TT         TT
#> 69          AA         AA         AA         AA         AA         AA
#> 70          AA         NN         AA         AA         AA         AA
#> 71          NN         AA         AA         AT         AA         AA
#> 72          GG         AG         AA         NN         GG         GG
#> 73          TT         GT         GT         TT         NN         GT
#> 74          GG         GG         GG         GG         GG         GG
#> 75          AC         AC         AA         AA         AA         AC
#> 76          TT         TT         TT         NN         NN         NN
#> 77          AT         AT         AT         AT         AA         AT
#> 78          GG         GG         GG         GG         GG         GG
#> 79          CT         TT         TT         CT         CT         TT
#> 80          CT         CT         CT         TT         CT         CT
#> 81          AA         AA         AA         AT         AT         AT
#> 82          TT         TT         TT         TT         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          TT         TT         TT         TT         TT         TT
#> 85          AA         AA         NN         AG         AA         AA
#> 86          NN         AG         AG         AG         NN         AG
#> 87          TT         TT         TT         TT         TT         TT
#> 88          GG         GG         GG         GG         GG         GG
#> 89          GG         CC         CG         GG         NN         CC
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CC         CC         CC         CC         CT         CC
#> 92          NN         AA         AA         NN         AG         NN
#> 93          TT         CT         TT         TT         TT         TT
#> 94          CG         CC         NN         NN         CG         CG
#> 95          GG         GG         GG         GG         GG         GG
#> 96          CC         CC         CC         CC         CC         CC
#> 97          GG         GG         GG         GG         GG         GG
#> 98          GT         GT         GG         NN         GT         GT
#> 99          GG         GG         NN         GG         GG         GG
#> 100         TT         AA         NN         AA         NN         AT
#>     PAL.r2.276 PAL.r2.278 PAL.r2.281 PAL.r2.285 PAL.r2.287 PAL.r2.290
#> 1           AG         GG         AG         GG         GG         AG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         NN         CT         CC         CC
#> 4           AA         AA         AA         AA         AA         AA
#> 5           AG         GG         NN         AG         NN         GG
#> 6           AC         AC         AC         CC         AA         AA
#> 7           TT         TT         TT         TT         TT         TT
#> 8           GG         GG         GG         NN         GG         GG
#> 9           CC         CC         CC         CC         CC         NN
#> 10          AA         AA         AA         AA         AA         AA
#> 11          TT         TT         TT         TT         TT         TT
#> 12          AG         AA         AA         AA         AG         AA
#> 13          TT         TT         TT         TT         TT         TT
#> 14          GG         GG         GG         GG         GG         GG
#> 15          AG         GG         AG         GG         GG         GG
#> 16          GT         GT         GG         GT         GT         GT
#> 17          AA         AA         AA         AA         AC         AC
#> 18          NN         AA         NN         AA         AC         AA
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         NN         GG         GG         GG         GG
#> 21          CC         CC         CC         CC         CC         CC
#> 22          CC         CC         CC         CC         CC         CC
#> 23          GG         GG         GG         GG         GG         GG
#> 24          AA         AA         AA         AA         AA         AA
#> 25          NN         GG         GG         GG         GG         GG
#> 26          CC         NN         CT         CT         CT         NN
#> 27          AA         AA         AA         AA         AA         AA
#> 28          TT         TT         NN         NN         TT         CT
#> 29          AA         AA         NN         NN         AA         AC
#> 30          GG         GG         GG         GG         GG         GG
#> 31          CC         CC         CC         CC         CC         CC
#> 32          AC         AC         AC         NN         AC         NN
#> 33          CT         CT         CT         CC         CC         CT
#> 34          CC         CC         CC         CC         CC         CC
#> 35          AA         NN         AA         AA         AA         AA
#> 36          AA         AA         AA         AA         AA         NN
#> 37          GG         GG         GG         GG         GG         GG
#> 38          GG         GT         GT         GG         NN         GT
#> 39          TT         NN         TT         NN         TT         AT
#> 40          CG         GG         NN         GG         GG         GG
#> 41          NN         AA         NN         AA         AA         AA
#> 42          AG         GG         GG         GG         GG         GG
#> 43          AA         AA         AA         AA         AA         AA
#> 44          NN         TT         TT         NN         AT         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          CC         CC         CC         CC         CC         CC
#> 47          CC         CC         CC         CC         CC         CC
#> 48          CT         NN         NN         TT         TT         TT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         NN         CC         CC         CC
#> 51          TT         NN         NN         TT         NN         TT
#> 52          NN         AG         AG         AG         NN         AG
#> 53          AA         AA         AA         NN         AA         AA
#> 54          GG         GG         GG         GG         GG         NN
#> 55          GG         AG         AG         NN         AG         AG
#> 56          AG         AA         NN         AA         AG         AA
#> 57          GG         GG         GG         GG         NN         NN
#> 58          CC         CC         CC         NN         CC         NN
#> 59          GG         GG         GG         AG         GG         AG
#> 60          AA         AA         AA         AA         AA         AA
#> 61          TT         TT         TT         TT         TT         TT
#> 62          GG         GG         GG         GG         GG         GG
#> 63          CC         CC         NN         NN         CC         AC
#> 64          AA         TT         AT         TT         NN         NN
#> 65          CC         CT         CT         TT         CT         CC
#> 66          CC         AC         NN         AC         AC         CC
#> 67          GG         NN         AG         GG         GG         AG
#> 68          NN         AT         AT         TT         AT         NN
#> 69          AA         AA         AA         AA         AA         AA
#> 70          AA         AA         AA         AA         AA         AA
#> 71          AA         AA         AT         AA         AA         AT
#> 72          NN         GG         AA         AG         AG         AG
#> 73          TT         GT         GG         TT         NN         GT
#> 74          GG         GG         GG         GG         GG         GG
#> 75          AA         AC         NN         AC         AC         AC
#> 76          NN         TT         NN         TT         TT         NN
#> 77          AT         NN         AA         NN         AT         AA
#> 78          GG         GG         GG         GG         GG         GG
#> 79          TT         TT         CT         CT         CT         TT
#> 80          CT         CT         CC         CT         CT         NN
#> 81          AT         AT         TT         AT         AA         AT
#> 82          TT         TT         TT         TT         TT         TT
#> 83          TT         TT         TT         TT         NN         TT
#> 84          TT         TT         TT         TT         TT         TT
#> 85          AG         NN         AA         AA         AA         AG
#> 86          AG         NN         AA         NN         AG         AA
#> 87          TT         TT         TT         TT         TT         TT
#> 88          GG         GG         GG         GG         GG         NN
#> 89          CG         NN         CC         GG         NN         GG
#> 90          CC         CC         CC         CC         CC         CC
#> 91          CC         CT         CC         CT         CC         CC
#> 92          AA         AG         NN         AG         AG         NN
#> 93          NN         TT         CT         TT         TT         CT
#> 94          CC         GG         CG         NN         CG         CC
#> 95          GG         GG         GG         GG         GG         GG
#> 96          CC         CC         CC         CC         CC         CC
#> 97          GG         GG         GG         GG         GG         GG
#> 98          GT         GG         GT         GT         TT         GT
#> 99          NN         GG         GG         GG         GG         GG
#> 100         AT         AT         AA         AT         AA         NN
#>     PAL.r2.291 SMR.r1.003 SMR.r1.004 SMR.r1.008 SMR.r1.012 SMR.r1.013
#> 1           AA         AA         GG         AA         AG         NN
#> 2           GG         GG         GT         GG         GG         GG
#> 3           CC         CC         CC         CT         CC         CC
#> 4           AA         AG         AA         AA         AG         NN
#> 5           GG         NN         NN         NN         AG         GG
#> 6           AC         AC         AA         AC         NN         AA
#> 7           TT         NN         TT         TT         NN         GT
#> 8           GG         AG         AG         GG         NN         GG
#> 9           CC         CC         AC         CC         NN         AC
#> 10          AA         AA         AA         AA         AA         NN
#> 11          TT         NN         CT         CC         CT         TT
#> 12          AG         NN         AA         AG         AA         AA
#> 13          TT         TT         TT         TT         TT         TT
#> 14          GG         GG         GG         GG         NN         GG
#> 15          GG         AA         AA         AG         AG         AA
#> 16          GG         GG         GG         GT         GG         NN
#> 17          AA         AC         AA         AC         AC         NN
#> 18          AC         AC         CC         CC         CC         NN
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         GG         GG         AG         AG         GG
#> 21          CC         AC         AC         CC         AC         AC
#> 22          CC         CC         CC         CT         CC         NN
#> 23          GG         GG         AG         GG         GG         GG
#> 24          AA         NN         AA         AG         AA         AA
#> 25          GG         AG         AA         AG         AA         AG
#> 26          CT         CT         NN         TT         TT         TT
#> 27          AA         AA         NN         AA         AA         AC
#> 28          TT         TT         TT         CC         CC         CT
#> 29          AA         NN         NN         AC         AA         AC
#> 30          GG         GG         GG         GG         CG         CG
#> 31          CC         NN         CC         GG         CG         NN
#> 32          AC         AA         AC         AC         AA         AA
#> 33          NN         CC         CT         CT         CC         CC
#> 34          CC         CC         CG         CC         CC         CG
#> 35          AA         NN         AA         AA         AA         NN
#> 36          AA         AA         AA         AG         AA         AG
#> 37          GG         CC         GG         CG         GG         CG
#> 38          TT         GT         GT         TT         GT         NN
#> 39          TT         NN         NN         AT         TT         TT
#> 40          GG         GG         GG         GG         GG         GG
#> 41          AA         AG         AG         AA         NN         AG
#> 42          GG         AG         AA         AG         GG         AG
#> 43          AA         AA         AA         AA         AA         AA
#> 44          TT         NN         TT         TT         NN         AA
#> 45          GG         GG         GG         GG         GG         GG
#> 46          NN         CC         NN         AA         CC         CC
#> 47          CC         CC         CC         CG         CC         CC
#> 48          CC         TT         TT         CC         CT         CT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CC         CC         CC         CT
#> 51          CT         NN         NN         CT         CT         NN
#> 52          AG         AA         AG         AA         AG         AA
#> 53          AA         AA         AA         AA         AA         AA
#> 54          GG         GG         GG         GG         GG         GG
#> 55          AA         GG         GG         AA         GG         AG
#> 56          AG         AA         NN         AG         AA         AA
#> 57          GG         NN         GG         GG         GG         GG
#> 58          TT         NN         CC         CC         CC         CC
#> 59          NN         NN         AG         NN         NN         AG
#> 60          AA         AA         NN         AC         AC         NN
#> 61          TT         TT         TT         TT         TT         TT
#> 62          GG         GG         GG         CG         GG         CG
#> 63          CC         NN         AC         AA         AA         AC
#> 64          AT         TT         AT         TT         AA         TT
#> 65          CT         CC         CT         CT         CT         NN
#> 66          AC         NN         AA         AC         AC         AC
#> 67          GG         NN         GG         GG         GG         GG
#> 68          TT         NN         AA         AT         AA         AT
#> 69          AA         AA         AA         AA         AA         AA
#> 70          AA         AC         AA         AA         NN         AC
#> 71          AT         AA         AA         AA         AA         AA
#> 72          AG         AG         AA         AG         AA         AA
#> 73          TT         TT         TT         TT         NN         TT
#> 74          GG         NN         GG         GT         GG         NN
#> 75          AC         AA         CC         AA         NN         CC
#> 76          TT         NN         NN         CC         CC         NN
#> 77          AT         AA         NN         AT         AA         AA
#> 78          GG         GG         GG         GG         GG         NN
#> 79          CT         TT         TT         CT         TT         NN
#> 80          CT         TT         TT         TT         CT         TT
#> 81          AT         TT         AA         AT         AA         AA
#> 82          TT         NN         NN         CT         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          TT         CC         TT         CC         CT         CT
#> 85          NN         GG         GG         AG         GG         GG
#> 86          AG         AA         AA         AG         AA         AA
#> 87          TT         CT         TT         TT         NN         CT
#> 88          NN         NN         AA         AG         AG         AA
#> 89          GG         NN         NN         CC         NN         CC
#> 90          CC         CC         CC         CC         CC         CC
#> 91          NN         CT         TT         CC         CT         TT
#> 92          AA         NN         NN         AG         GG         AG
#> 93          NN         CT         TT         TT         NN         TT
#> 94          CC         CC         CC         CC         CC         CC
#> 95          GG         NN         GT         GT         GG         GG
#> 96          CC         AC         CC         CC         CC         CC
#> 97          GG         NN         GG         GG         GG         GG
#> 98          TT         GT         GG         GG         GT         GG
#> 99          NN         AG         GG         AG         AG         AG
#> 100         AA         NN         AT         AT         AA         AT
#>     SMR.r1.023 SMR.r1.028 SMR.r1.034 SMR.r1.049 SMR.r1.052 SMR.r1.056
#> 1           AG         AA         GG         AG         AG         GG
#> 2           GG         GG         GG         NN         GG         GG
#> 3           CC         NN         CC         CC         CC         CC
#> 4           AG         AA         AG         GG         AG         AA
#> 5           NN         NN         AG         AG         GG         GG
#> 6           AA         NN         AA         AA         AA         AA
#> 7           GT         TT         TT         GT         NN         TT
#> 8           AG         AG         GG         AG         NN         AA
#> 9           AC         AC         AC         CC         CC         AC
#> 10          AA         AG         NN         AA         AG         AA
#> 11          CT         NN         CT         CC         CT         TT
#> 12          AA         NN         AG         AA         AA         NN
#> 13          TT         CT         CT         CT         CT         TT
#> 14          GG         GT         GG         GG         GT         GG
#> 15          AG         AA         AA         AG         AG         NN
#> 16          GT         GT         GG         TT         GT         NN
#> 17          AA         AA         AA         AA         AA         NN
#> 18          CC         CC         CC         CC         CC         CC
#> 19          TT         TT         TT         TT         TT         TT
#> 20          GG         GG         GG         GG         GG         AG
#> 21          CC         CC         AC         AC         NN         CC
#> 22          NN         CT         CC         CC         CC         CC
#> 23          AG         GG         GG         NN         AG         GG
#> 24          AA         AA         AA         AA         AA         AG
#> 25          NN         GG         NN         NN         NN         AA
#> 26          CT         CT         TT         CT         CT         CT
#> 27          NN         AC         NN         AA         AA         AA
#> 28          CC         NN         CT         CT         TT         NN
#> 29          AA         NN         NN         AC         NN         NN
#> 30          GG         GG         NN         GG         GG         CG
#> 31          CG         CC         CC         GG         NN         CG
#> 32          AA         AC         NN         AC         AA         AA
#> 33          CT         CC         CT         CC         NN         CC
#> 34          NN         CC         CG         CC         CC         GG
#> 35          AA         AC         AA         AA         AA         AC
#> 36          AA         AA         AA         AG         AA         AA
#> 37          CG         CG         GG         CG         GG         NN
#> 38          NN         NN         GG         GG         GT         GT
#> 39          TT         NN         TT         TT         NN         NN
#> 40          CG         NN         GG         GG         GG         GG
#> 41          NN         NN         GG         AG         AG         GG
#> 42          GG         NN         AG         AG         AG         AG
#> 43          AA         AG         AA         AA         AA         AG
#> 44          NN         TT         TT         TT         TT         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          AC         CC         NN         AC         NN         CC
#> 47          CG         CG         CC         CC         CC         CC
#> 48          TT         CT         NN         CT         NN         TT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          NN         CC         CC         CC         CC         NN
#> 51          NN         CT         NN         NN         TT         CT
#> 52          AA         GG         AA         AG         AG         AA
#> 53          AA         AC         NN         AC         AA         CC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          GG         GG         GG         AG         GG         GG
#> 56          NN         AA         AG         AA         AA         NN
#> 57          GG         NN         GG         NN         GG         NN
#> 58          NN         CC         NN         CT         NN         CC
#> 59          AG         GG         GG         AG         AG         AG
#> 60          AA         AC         AA         AA         AC         AC
#> 61          TT         NN         NN         TT         CT         TT
#> 62          GG         GG         GG         GG         GG         GG
#> 63          AC         NN         AC         AC         NN         NN
#> 64          AT         TT         AT         TT         AT         AA
#> 65          CT         CC         TT         CC         CT         CC
#> 66          NN         NN         AA         NN         AC         AA
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AT         AA         AA         NN         AA         NN
#> 69          AG         AA         AA         AA         AA         NN
#> 70          AC         AA         AC         AA         AA         AA
#> 71          AA         AT         AA         AT         AA         AA
#> 72          NN         AA         AG         AA         NN         AG
#> 73          GT         TT         GT         TT         TT         TT
#> 74          GT         GT         NN         GT         NN         GG
#> 75          AC         AC         AC         AA         AC         NN
#> 76          CC         CC         NN         NN         NN         CC
#> 77          AA         AA         AT         AA         NN         AA
#> 78          GT         GG         GG         GG         GG         GG
#> 79          CT         TT         TT         TT         TT         TT
#> 80          CT         TT         CT         TT         CT         CT
#> 81          NN         AT         AT         AA         AT         AT
#> 82          TT         TT         TT         TT         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          NN         NN         CC         NN         NN         CC
#> 85          NN         GG         NN         GG         AG         AG
#> 86          AA         AA         AG         AA         NN         AA
#> 87          NN         CT         CT         CT         TT         TT
#> 88          AG         NN         NN         AG         AG         AG
#> 89          CG         NN         CC         CC         CG         NN
#> 90          CC         CC         CC         CC         CC         CC
#> 91          NN         CC         CT         CC         CC         CT
#> 92          NN         NN         NN         GG         GG         NN
#> 93          NN         NN         CT         TT         TT         TT
#> 94          CC         CC         CC         CC         NN         NN
#> 95          GT         NN         GG         GG         GG         NN
#> 96          AC         CC         CC         CC         CC         CC
#> 97          GG         GG         GG         GG         GG         NN
#> 98          NN         GT         GG         NN         GG         GT
#> 99          AA         NN         AG         AA         GG         AG
#> 100         NN         NN         AT         AA         AA         AT
#>     SMR.r1.072 SMR.r1.074 SMR.r2.001 SMR.r2.077 SMR.r2.086 SMR.r2.091
#> 1           AA         AG         AG         AG         AA         AA
#> 2           GG         GG         GG         NN         GG         GT
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AA         AG         AA         AA         AG         AA
#> 5           NN         NN         GG         NN         GG         GG
#> 6           AA         AC         AA         AA         AA         NN
#> 7           TT         TT         GT         TT         GT         TT
#> 8           AG         AG         AG         AG         NN         GG
#> 9           NN         AC         CC         CC         AC         NN
#> 10          AA         AA         AA         AA         AA         AA
#> 11          CT         CT         NN         NN         TT         CT
#> 12          AA         NN         NN         AA         AG         AG
#> 13          CT         TT         NN         TT         CC         NN
#> 14          GT         GG         GG         GT         GT         NN
#> 15          AA         AA         AG         AA         AA         AA
#> 16          GT         GG         GG         GT         GT         NN
#> 17          AC         AA         AA         AA         AC         AA
#> 18          AC         CC         CC         AC         AA         AC
#> 19          GT         TT         GT         NN         TT         TT
#> 20          GG         GG         NN         GG         GG         GG
#> 21          AC         AC         NN         AC         CC         AC
#> 22          NN         CT         CC         CC         NN         CC
#> 23          NN         GG         NN         AG         GG         GG
#> 24          AA         AA         NN         NN         AA         AA
#> 25          AG         GG         AG         NN         AG         AG
#> 26          CT         CT         CT         CT         TT         CT
#> 27          AC         AA         AC         AA         AA         AA
#> 28          NN         CT         CT         CT         NN         CT
#> 29          AA         AA         NN         NN         AA         AC
#> 30          GG         GG         GG         CG         GG         GG
#> 31          CC         CG         CG         CG         CG         NN
#> 32          AA         AA         AC         AA         NN         AC
#> 33          CC         TT         NN         CC         CC         TT
#> 34          CC         CC         NN         CC         CG         CG
#> 35          NN         AC         AA         AA         AA         AA
#> 36          AG         AA         AA         AA         AG         AA
#> 37          CG         GG         GG         CG         GG         CC
#> 38          GT         GT         GG         GG         GT         GG
#> 39          TT         NN         TT         TT         TT         AT
#> 40          CG         GG         NN         GG         NN         GG
#> 41          AG         AG         NN         NN         GG         GG
#> 42          AA         AG         AA         AG         AA         AG
#> 43          AA         AA         AA         AG         AA         AA
#> 44          TT         TT         TT         TT         TT         TT
#> 45          GG         GG         NN         NN         GG         NN
#> 46          AC         NN         NN         CC         CC         CC
#> 47          CG         CC         CC         NN         CG         CG
#> 48          CT         CT         NN         NN         CC         CT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          NN         NN         CC         CC         NN         CC
#> 51          NN         CT         CT         TT         CT         TT
#> 52          AA         AA         AA         AA         AA         AA
#> 53          AA         AA         AA         AA         AA         AA
#> 54          GG         GG         NN         GG         GT         GG
#> 55          GG         AG         NN         NN         GG         GG
#> 56          NN         AG         AA         AA         AG         AG
#> 57          GG         GG         AG         GG         GG         GG
#> 58          CC         CT         NN         CC         CC         CC
#> 59          AA         AG         AG         GG         AG         AG
#> 60          CC         AA         CC         AC         AA         AC
#> 61          TT         TT         NN         TT         TT         TT
#> 62          GG         GG         GG         GG         GG         GG
#> 63          AA         AA         NN         NN         AC         NN
#> 64          AT         AA         TT         NN         AT         NN
#> 65          CC         NN         CT         CC         CT         CT
#> 66          NN         AA         CC         AC         CC         AA
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AT         AT         AA         AA         AT         AT
#> 69          AA         AA         AA         AA         AG         AA
#> 70          AA         AA         AC         AA         AC         AC
#> 71          AA         AA         AA         AA         AT         AA
#> 72          AA         AG         AG         AA         AG         AA
#> 73          TT         GT         GT         TT         TT         TT
#> 74          GG         NN         GG         NN         GG         GT
#> 75          AC         AC         CC         AA         AC         AC
#> 76          CT         NN         NN         NN         CC         CC
#> 77          AA         AT         NN         NN         AT         AT
#> 78          GG         GG         GG         GG         GT         GG
#> 79          TT         NN         TT         TT         TT         TT
#> 80          CT         CT         CT         CT         CT         CT
#> 81          AA         NN         AA         AT         AT         AA
#> 82          TT         NN         NN         CT         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          CT         CT         CC         CC         CT         CT
#> 85          AG         GG         AG         GG         AG         AG
#> 86          AA         AG         NN         NN         AG         AG
#> 87          TT         TT         CT         TT         CT         TT
#> 88          AG         NN         NN         NN         GG         AA
#> 89          CG         CG         CC         CG         CG         GG
#> 90          CC         CC         CC         CC         CC         CC
#> 91          TT         CT         NN         NN         CT         CT
#> 92          NN         NN         NN         AG         NN         AG
#> 93          CT         TT         TT         CT         TT         TT
#> 94          CC         CC         CC         CC         CC         CC
#> 95          GT         GT         NN         GT         GG         GG
#> 96          CC         CC         CC         CC         AC         CC
#> 97          GG         NN         GG         GG         GG         AG
#> 98          GG         GG         GG         GT         GT         GG
#> 99          AA         NN         AA         AG         GG         AG
#> 100         AT         NN         AT         AT         AT         AA
#>     SMR.r2.097 SMR.r2.105 SMR.r2.139 SMR.r2.175 SMR.r2.256 SMR.r2.263
#> 1           NN         AA         AG         AG         AG         AA
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CT         CC         CC         CC         CT         CC
#> 4           AA         GG         AA         AG         AG         GG
#> 5           GG         GG         NN         AA         AG         GG
#> 6           AA         AA         AA         AA         AA         AA
#> 7           TT         TT         GT         TT         TT         GT
#> 8           GG         AA         NN         GG         GG         AG
#> 9           AC         NN         CC         AC         CC         CC
#> 10          AA         AA         AA         AA         AA         AA
#> 11          CT         CT         NN         CC         CC         CT
#> 12          AG         AA         NN         NN         AG         AA
#> 13          CT         TT         TT         NN         NN         TT
#> 14          GG         GT         GG         GT         GG         GT
#> 15          AA         AG         AG         NN         AA         NN
#> 16          GG         NN         GT         GT         GG         GT
#> 17          AA         AA         AA         AC         NN         AC
#> 18          AC         CC         CC         CC         CC         CC
#> 19          TT         TT         GT         TT         TT         TT
#> 20          GG         GG         GG         GG         GG         GG
#> 21          AC         AC         NN         NN         CC         CC
#> 22          CT         CC         CT         CC         CT         CC
#> 23          GG         GG         GG         GG         NN         GG
#> 24          AA         AA         NN         NN         AA         AA
#> 25          AA         AG         GG         GG         AG         AG
#> 26          TT         TT         CT         TT         TT         NN
#> 27          NN         AA         AA         AA         AC         AC
#> 28          CT         CT         TT         CT         CT         CC
#> 29          CC         AA         NN         AC         AA         AA
#> 30          GG         GG         NN         CG         GG         GG
#> 31          CC         CG         CG         NN         CC         CC
#> 32          AA         AA         AA         AA         NN         NN
#> 33          CT         TT         CC         TT         CC         CT
#> 34          CC         CG         GG         CC         CC         CG
#> 35          AA         AA         NN         AA         AA         AA
#> 36          NN         NN         AG         AA         AA         AG
#> 37          CG         CG         NN         GG         GG         CG
#> 38          NN         GT         GT         GT         GT         TT
#> 39          TT         TT         TT         NN         AT         TT
#> 40          CG         GG         NN         NN         CG         CG
#> 41          AG         GG         NN         GG         NN         GG
#> 42          AG         AG         NN         AG         AA         AG
#> 43          AA         AA         AA         AG         AA         AA
#> 44          AT         AT         TT         AT         TT         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          NN         CC         CC         NN         CC         CC
#> 47          GG         CC         CC         CC         CC         CG
#> 48          CT         TT         CT         CT         NN         CC
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CC         CT         CC         CC
#> 51          NN         CC         TT         TT         NN         NN
#> 52          AA         AG         AA         AA         AA         AA
#> 53          AC         AC         AA         AA         AC         AC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          NN         AG         GG         NN         GG         GG
#> 56          AG         AA         AA         AA         AG         NN
#> 57          GG         GG         NN         GG         GG         GG
#> 58          CT         CT         CC         CC         CT         CC
#> 59          GG         AG         AG         NN         AA         AG
#> 60          AC         AA         AC         AC         NN         AA
#> 61          TT         NN         TT         TT         CT         TT
#> 62          GG         GG         CG         GG         GG         GG
#> 63          CC         AA         AC         AA         AA         AC
#> 64          TT         AT         AT         AT         AT         TT
#> 65          CC         TT         CT         CC         CC         TT
#> 66          AA         AC         NN         NN         AA         CC
#> 67          GG         GG         GG         GG         GG         NN
#> 68          TT         AA         AT         NN         NN         NN
#> 69          AA         NN         AG         AA         AA         AG
#> 70          AA         AC         AC         NN         AA         AC
#> 71          AT         AT         AA         AA         AA         AT
#> 72          AA         AG         AA         AG         AA         AA
#> 73          NN         GT         TT         TT         TT         TT
#> 74          GT         GG         NN         GT         GT         GG
#> 75          AA         CC         AC         NN         NN         AC
#> 76          CC         NN         CC         CT         NN         CT
#> 77          AT         AA         AA         NN         AT         AA
#> 78          GG         NN         GT         GG         GG         GT
#> 79          TT         TT         TT         CC         TT         TT
#> 80          NN         TT         TT         TT         TT         CT
#> 81          AT         AA         AT         AA         AT         AT
#> 82          TT         TT         TT         NN         TT         TT
#> 83          TT         TT         TT         TT         TT         TT
#> 84          CC         CC         CC         NN         CC         CT
#> 85          GG         GG         AG         AG         GG         AA
#> 86          AG         AA         NN         AA         AG         AA
#> 87          CT         CT         CT         CT         TT         TT
#> 88          AA         NN         AG         AG         AG         AG
#> 89          CC         CC         CG         CG         NN         GG
#> 90          CC         CT         CC         CC         CC         CC
#> 91          TT         TT         TT         CT         NN         CT
#> 92          GG         GG         NN         NN         GG         NN
#> 93          TT         CT         NN         NN         CT         TT
#> 94          CC         CC         CC         CG         CC         CC
#> 95          GG         GG         NN         GG         GG         GG
#> 96          NN         CC         CC         CC         CC         AC
#> 97          NN         GG         GG         GG         GG         GG
#> 98          GG         GG         GT         GG         GG         GT
#> 99          AA         AG         GG         NN         NN         GG
#> 100         AT         AT         AA         AT         AT         AA
#>     SMR.r2.296 UPD.r2.120 UPD.r2.121 UPD.r2.122 UPD.r2.123 UPD.r2.124
#> 1           GG         AA         AA         AG         AG         AG
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CT         CT         CC         CC         CC
#> 4           AA         AG         GG         AG         AG         AG
#> 5           AG         GG         NN         AA         AG         AG
#> 6           AA         AC         AC         AC         AC         CC
#> 7           GT         TT         TT         TT         TT         TT
#> 8           GG         NN         GG         AG         AG         NN
#> 9           CC         CC         CC         AC         CC         NN
#> 10          AA         AA         AA         AA         AA         AG
#> 11          TT         CT         CT         TT         TT         CT
#> 12          AA         AG         AA         AG         AG         AA
#> 13          TT         TT         CT         TT         TT         TT
#> 14          GG         GG         GT         GT         GT         GT
#> 15          AG         AG         AG         GG         NN         AG
#> 16          GT         GG         GG         GT         GT         GG
#> 17          AA         AC         AA         AA         AC         AA
#> 18          CC         AA         CC         AA         CC         CC
#> 19          TT         GT         TT         TT         TT         TT
#> 20          GG         GG         GG         NN         AG         GG
#> 21          AC         CC         AC         AC         NN         CC
#> 22          CC         CC         TT         CT         CT         CT
#> 23          GG         GG         GG         GG         GG         AG
#> 24          AA         AA         AA         AA         AA         AA
#> 25          AG         GG         AG         AG         AA         GG
#> 26          CT         CT         CT         CT         TT         CT
#> 27          AA         AA         AA         AC         AA         AA
#> 28          CC         TT         CT         CT         TT         TT
#> 29          AA         AA         AA         AA         AC         AA
#> 30          GG         GG         GG         GG         NN         GG
#> 31          CG         GG         CG         CG         CC         CC
#> 32          AC         AC         AA         AA         AC         AC
#> 33          CT         CC         TT         CC         TT         CC
#> 34          GG         CG         CC         CG         CC         CC
#> 35          AA         AC         AA         AA         AA         AA
#> 36          AG         AA         AA         AA         AG         AA
#> 37          CC         GG         GG         GG         GG         CG
#> 38          NN         GG         GG         GT         GG         GT
#> 39          TT         TT         AT         TT         TT         TT
#> 40          CG         CG         CG         GG         CG         GG
#> 41          AG         AA         NN         AA         NN         GG
#> 42          AG         AG         AG         GG         GG         AG
#> 43          AG         NN         AA         AG         AG         NN
#> 44          TT         TT         AT         TT         NN         TT
#> 45          NN         GG         GG         GG         GG         GG
#> 46          AA         CC         AC         CC         NN         CC
#> 47          GG         CC         CC         CG         CC         CC
#> 48          NN         CT         NN         TT         CT         CT
#> 49          GG         GG         GG         GG         GG         GG
#> 50          CC         CC         CC         CT         CC         CC
#> 51          TT         CC         TT         CC         CT         TT
#> 52          AA         GG         GG         AA         AG         AG
#> 53          AA         AA         CC         AA         AC         NN
#> 54          GG         GG         GG         GG         GG         GG
#> 55          AG         AA         AG         AG         AG         GG
#> 56          AA         AG         AA         AG         AG         AA
#> 57          NN         GG         GG         GG         GG         GG
#> 58          NN         NN         TT         CT         NN         NN
#> 59          AA         GG         AG         GG         AG         GG
#> 60          AA         AC         AA         NN         AC         AC
#> 61          TT         TT         CT         TT         TT         TT
#> 62          GG         GG         NN         GG         GG         GG
#> 63          NN         AC         AC         AC         AA         CC
#> 64          TT         AA         TT         TT         AT         AT
#> 65          CC         TT         CC         TT         CC         CT
#> 66          AC         AC         AC         AA         AA         AC
#> 67          GG         GG         GG         GG         GG         GG
#> 68          AA         AA         TT         AT         AA         AT
#> 69          NN         AA         AA         AA         AA         AA
#> 70          AC         AA         AC         NN         AC         AA
#> 71          AA         AA         AA         AA         AA         AA
#> 72          AG         AG         GG         AG         AG         NN
#> 73          NN         GT         TT         TT         GT         GT
#> 74          GG         GT         GG         GG         GG         GG
#> 75          AC         AC         AC         AC         AC         AC
#> 76          NN         CT         CC         CC         CT         CT
#> 77          AA         AT         AA         AT         AT         AA
#> 78          GT         GG         GG         GG         GG         GG
#> 79          TT         CT         TT         TT         TT         NN
#> 80          TT         TT         CT         NN         CT         CT
#> 81          AT         AT         AA         AT         NN         AT
#> 82          TT         TT         TT         TT         CT         TT
#> 83          AT         TT         TT         TT         TT         TT
#> 84          CT         CT         TT         CC         CT         CT
#> 85          NN         GG         GG         AG         GG         AG
#> 86          NN         AG         AA         AG         AG         NN
#> 87          CT         CT         TT         NN         TT         TT
#> 88          NN         NN         AG         AG         AG         AG
#> 89          CG         GG         CC         CC         CG         NN
#> 90          NN         CC         CC         CC         CT         NN
#> 91          CT         CT         CT         CC         TT         TT
#> 92          GG         NN         GG         GG         GG         AG
#> 93          CT         TT         TT         TT         CT         CT
#> 94          CC         CC         CC         CC         CC         CG
#> 95          GG         GG         GT         GG         GG         GG
#> 96          CC         CC         CC         CC         NN         CC
#> 97          AG         GG         GG         GG         GG         GG
#> 98          GT         GT         GG         GG         GT         GT
#> 99          AA         AG         GG         AA         NN         NN
#> 100         AA         AA         AT         AA         AT         AT
#>     UPD.r2.129 UPD.r2.130 UPD.r2.131 UPD.r2.133 UPD.r2.135 UPD.r2.160
#> 1           AA         AG         AG         GG         AG         NN
#> 2           GG         GG         GG         GG         GG         GG
#> 3           CC         CC         CC         CC         CC         CC
#> 4           AG         AG         NN         NN         AG         AG
#> 5           NN         NN         NN         AA         AG         AG
#> 6           AC         NN         AA         AC         AC         AA
#> 7           TT         TT         TT         TT         TT         TT
#> 8           GG         NN         NN         GG         AG         AG
#> 9           CC         AC         CC         AC         AC         CC
#> 10          NN         AG         AA         AA         AA         AA
#> 11          NN         CT         CT         CT         CT         CT
#> 12          NN         AA         AG         NN         NN         NN
#> 13          TT         TT         TT         CT         CT         CT
#> 14          GT         GG         GG         GT         GG         TT
#> 15          AG         AG         NN         AG         AG         AG
#> 16          GT         NN         GG         NN         GG         GG
#> 17          NN         AA         AA         AA         AC         NN
#> 18          CC         AA         CC         AC         AC         CC
#> 19          TT         TT         TT         GT         TT         TT
#> 20          GG         GG         GG         GG         GG         GG
#> 21          NN         AC         CC         AC         CC         AC
#> 22          CT         CC         NN         CC         CT         CT
#> 23          GG         GG         GG         GG         GG         GG
#> 24          AA         AA         AG         AG         AA         NN
#> 25          AG         AG         GG         AA         GG         AG
#> 26          TT         CT         CT         TT         CT         NN
#> 27          AC         AA         AA         AA         AA         AA
#> 28          TT         CT         CT         NN         TT         TT
#> 29          NN         AA         AA         NN         AA         AA
#> 30          GG         CG         CG         NN         GG         NN
#> 31          CC         CG         CG         GG         CC         GG
#> 32          AA         AA         AC         AC         AC         AA
#> 33          CT         CC         NN         CT         CC         CC
#> 34          NN         CC         CC         CC         CC         NN
#> 35          AA         AA         AA         AA         AA         AA
#> 36          AA         AA         AG         AA         NN         AA
#> 37          CG         GG         CG         CG         CG         CG
#> 38          GG         GG         TT         TT         GG         GT
#> 39          NN         NN         AT         NN         TT         AT
#> 40          NN         CG         CG         NN         CC         GG
#> 41          NN         AG         AA         AA         AG         NN
#> 42          GG         AG         AA         GG         NN         AG
#> 43          AG         AA         AG         AA         AG         AG
#> 44          AT         AT         NN         TT         NN         TT
#> 45          GG         GG         GG         GG         GG         GG
#> 46          CC         AC         AA         AC         CC         AC
#> 47          CC         CC         CG         CC         CC         CC
#> 48          CT         TT         CT         CT         CT         TT
#> 49          GG         NN         CG         GG         GG         NN
#> 50          CC         CT         CC         CC         CC         NN
#> 51          CT         CT         CT         TT         CT         NN
#> 52          NN         AG         AG         AA         AA         AG
#> 53          NN         AA         AC         AC         AC         AC
#> 54          GG         GG         GG         GG         GG         GG
#> 55          GG         AG         AA         GG         NN         NN
#> 56          AG         AA         AG         AA         AA         AA
#> 57          NN         GG         NN         GG         NN         GG
#> 58          CT         CT         CC         CC         CC         NN
#> 59          GG         GG         AG         AG         AG         AA
#> 60          AC         AC         AA         AA         AA         AC
#> 61          NN         TT         TT         TT         NN         TT
#> 62          GG         CG         CG         CG         GG         GG
#> 63          AA         AC         AC         AA         NN         AA
#> 64          TT         TT         AT         AT         AT         AT
#> 65          CC         CT         NN         CC         CC         CC
#> 66          NN         AA         AA         AA         CC         CC
#> 67          GG         GG         GG         NN         NN         AG
#> 68          AT         AA         AA         AT         AT         NN
#> 69          AA         AA         AA         AA         AA         AA
#> 70          AA         AC         AC         AA         AA         AC
#> 71          AA         AA         AA         AT         AA         AA
#> 72          AG         GG         NN         GG         AG         GG
#> 73          TT         GG         TT         TT         GT         TT
#> 74          GG         NN         GG         GG         GG         GT
#> 75          AC         AC         CC         CC         AC         CC
#> 76          CT         CT         CC         CC         TT         CT
#> 77          NN         AA         AT         AA         AA         AA
#> 78          GG         GG         GG         GG         GG         GG
#> 79          TT         TT         TT         NN         NN         TT
#> 80          CT         CT         TT         CT         TT         NN
#> 81          AT         AT         AA         AA         AT         AT
#> 82          TT         NN         CT         TT         CT         CT
#> 83          TT         TT         AT         TT         TT         TT
#> 84          CT         CT         CT         TT         CT         NN
#> 85          AA         AG         AG         AA         AG         AG
#> 86          AG         AA         AG         AA         NN         AA
#> 87          CT         CT         CC         TT         CT         CT
#> 88          AA         AA         AG         AA         NN         AG
#> 89          GG         CG         NN         CG         NN         CG
#> 90          CC         CC         CT         CC         CC         CC
#> 91          CT         CC         CC         CC         CC         CT
#> 92          NN         NN         NN         AG         GG         AG
#> 93          TT         TT         CT         TT         TT         CT
#> 94          CC         CG         CC         CC         NN         CG
#> 95          GT         GT         GT         GG         GG         NN
#> 96          CC         CC         CC         NN         AC         CC
#> 97          GG         GG         AG         GG         GG         GG
#> 98          NN         GG         GG         GT         NN         GT
#> 99          GG         AG         AG         GG         NN         NN
#> 100         AT         TT         TT         TT         NN         TT
#>     UPD.r2.162 UPD.r2.174 UPD.r2.189 UPD.r2.201
#> 1           GG         NN         GG         AG
#> 2           GG         GG         GG         GG
#> 3           CC         CC         CC         CC
#> 4           AG         AA         AA         AA
#> 5           NN         GG         NN         GG
#> 6           AC         AA         AA         AA
#> 7           NN         TT         TT         TT
#> 8           NN         AG         GG         GG
#> 9           NN         CC         AC         AA
#> 10          AA         AG         AA         AG
#> 11          CT         CT         CT         TT
#> 12          AG         AA         AG         NN
#> 13          TT         TT         NN         NN
#> 14          GT         GG         GT         GT
#> 15          AG         AG         AA         AA
#> 16          GG         TT         NN         GG
#> 17          AA         AC         AA         AA
#> 18          AC         AA         AC         AC
#> 19          TT         TT         TT         TT
#> 20          GG         GG         AG         GG
#> 21          NN         AC         AA         CC
#> 22          CT         CT         CC         CC
#> 23          GG         GG         GG         GG
#> 24          AG         AA         AA         AA
#> 25          AG         AA         AG         GG
#> 26          NN         NN         TT         TT
#> 27          AA         NN         AA         AA
#> 28          CT         CT         CT         CT
#> 29          AA         AA         AC         AC
#> 30          GG         CG         GG         CG
#> 31          CG         CG         CC         CG
#> 32          NN         AA         NN         AC
#> 33          CC         CT         NN         TT
#> 34          CG         NN         CC         NN
#> 35          NN         AA         AA         AC
#> 36          NN         AA         AA         AA
#> 37          NN         CG         NN         GG
#> 38          GG         GG         GG         NN
#> 39          TT         NN         TT         NN
#> 40          GG         NN         CC         NN
#> 41          NN         GG         AG         AG
#> 42          AG         GG         AG         AG
#> 43          AA         AA         AA         AA
#> 44          TT         AT         AT         TT
#> 45          GG         GT         GG         GT
#> 46          CC         AA         NN         AC
#> 47          GG         CC         CC         CG
#> 48          CT         NN         CC         TT
#> 49          GG         CG         GG         GG
#> 50          CC         CC         CC         CC
#> 51          CT         CT         TT         TT
#> 52          AG         AA         AA         AA
#> 53          NN         AA         AC         CC
#> 54          GG         GG         GG         GG
#> 55          AA         AG         AG         GG
#> 56          AG         AA         AG         AA
#> 57          GG         NN         NN         GG
#> 58          CC         CC         CT         CC
#> 59          AG         AA         NN         NN
#> 60          AC         AC         AA         AC
#> 61          TT         CT         TT         NN
#> 62          GG         GG         GG         GG
#> 63          AA         AA         NN         AA
#> 64          TT         AT         AT         AT
#> 65          CT         CT         CT         TT
#> 66          AC         CC         NN         AC
#> 67          GG         GG         GG         GG
#> 68          TT         AA         AA         AT
#> 69          AA         AA         AA         AA
#> 70          NN         NN         AC         AC
#> 71          AA         AT         AA         AA
#> 72          GG         GG         AG         AG
#> 73          TT         TT         NN         GT
#> 74          GG         GT         NN         GG
#> 75          CC         CC         AC         CC
#> 76          TT         CC         NN         TT
#> 77          AT         AA         AT         AA
#> 78          GG         GG         GG         GG
#> 79          CT         CT         TT         TT
#> 80          CT         TT         TT         TT
#> 81          TT         AA         AT         AT
#> 82          TT         TT         NN         CT
#> 83          TT         TT         NN         TT
#> 84          CC         CT         NN         CC
#> 85          AG         GG         AG         AG
#> 86          AG         AA         AG         AA
#> 87          CT         CT         CC         CT
#> 88          NN         AG         AA         AG
#> 89          GG         GG         GG         CC
#> 90          CC         CC         NN         CC
#> 91          NN         CC         CT         CT
#> 92          GG         NN         AG         GG
#> 93          NN         CT         CT         CT
#> 94          CC         CC         CC         CC
#> 95          GG         NN         NN         GG
#> 96          CC         CC         AC         CC
#> 97          GG         GG         GG         GG
#> 98          GG         GT         GT         TT
#> 99          AA         AA         AG         AG
#> 100         TT         AT         AT         AT

# show or overwrite snp meta
snp.meta(test)
#>            chr position .snp.id
#> 1      groupVI   212436       1
#> 2   groupXVIII   350070       2
#> 3      groupXV   854285       3
#> 4    groupXIII   907279       4
#> 5      groupVI  1068620       5
#> 6       groupV  1392753       6
#> 7      groupIX  2075079       7
#> 8    groupXIII  2106000       8
#> 9     groupXII  2136587       9
#> 10   groupXVII  2560839      10
#> 11    groupIII  2750695      11
#> 12    groupXIX  3008524      12
#> 13    groupXVI  3624417      13
#> 14    groupVII  3667528      14
#> 15    groupVII  4189120      15
#> 16     groupXX  4538074      16
#> 17   groupXIII  4706089      17
#> 18    groupIII  4746808      18
#> 19   groupXVII  5041995      19
#> 20      groupX  5161637      20
#> 21    groupXII  5178037      21
#> 22     groupIX  5287711      22
#> 23   groupVIII  5462135      23
#> 24     groupIX  5960933      24
#> 25   groupXIII  6076868      25
#> 26     groupXI  6177332      26
#> 27     groupIV  6291287      27
#> 28     groupIX  6515796      28
#> 29    groupXII  6718079      29
#> 30    groupVII  6967477      30
#> 31    groupIII  7194069      31
#> 32    groupXVI  7199722      32
#> 33   groupVIII  7224065      33
#> 34  groupXVIII  7229426      34
#> 35  groupXVIII  8087022      35
#> 36     groupIX  8373442      36
#> 37     groupXV  8815538      37
#> 38     groupXI  9019877      38
#> 39    groupXXI  9067876      39
#> 40   groupXVII  9108184      40
#> 41     groupIX  9284427      41
#> 42      groupV  9767769      42
#> 43    groupXII  9881312      43
#> 44     groupXX 10296516      44
#> 45      groupI 10593829      45
#> 46    groupXIX 10765968      46
#> 47     groupII 10829169      47
#> 48     groupXX 10913115      48
#> 49  groupXVIII 10974594      49
#> 50   groupXVII 11256661      50
#> 51   groupVIII 11438870      51
#> 52    groupVII 11477223      52
#> 53     groupII 11505084      53
#> 54     groupIX 11904651      54
#> 55   groupVIII 12108945      55
#> 56    groupXIX 12485182      56
#> 57     groupXI 12833733      57
#> 58  groupXVIII 12937786      58
#> 59    groupIII 13177027      59
#> 60     groupVI 13335577      60
#> 61      groupX 13718292      61
#> 62     groupVI 13895716      62
#> 63    groupXIV 14009990      63
#> 64    groupXII 14158350      64
#> 65     groupIX 14164215      65
#> 66     groupVI 14236496      66
#> 67      groupX 14273796      67
#> 68   groupVIII 14445389      68
#> 69     groupIX 14868533      69
#> 70     groupIX 15190693      70
#> 71    groupXII 15306362      71
#> 72     groupXX 15460430      72
#> 73    groupVII 15493025      73
#> 74     groupXV 15547073      74
#> 75   groupXIII 15697285      75
#> 76    groupIII 15983994      76
#> 77    groupXIX 16233699      77
#> 78     groupIX 16306966      78
#> 79      groupI 16498302      79
#> 80   groupXIII 16520929      80
#> 81    groupXVI 16561891      81
#> 82     groupII 16598364      82
#> 83     groupIV 17248473      83
#> 84     groupXX 17748944      84
#> 85   groupXIII 17842228      85
#> 86    groupXIX 17858769      86
#> 87     groupIX 17929225      87
#> 88     groupIX 18339767      88
#> 89   groupXIII 18836558      89
#> 90     groupII 19244005      90
#> 91     groupXX 19317273      91
#> 92      groupI 19759331      92
#> 93    groupVII 20165240      93
#> 94     groupII 20496069      94
#> 95      groupI 21284266      95
#> 96     groupIV 22227311      96
#> 97      groupI 22746027      97
#> 98    groupVII 25239495      98
#> 99    groupVII 27599144      99
#> 100    groupIV 32321895     100
snp.meta(test) <- data.frame(pos = sample(10000, nrow(test), replace = TRUE),
                             chr = sample(LETTERS[1:4], nrow(test), 
                                          replace = TRUE))
#> Warning: Some levels are duplicated across multiple sample meta facets.
#> This will cause issues if those sample facets are run during analysis.
#> Issues:
#> 
#> Level: A in facets: fam, chr
#> Level: B in facets: fam, chr
#> Level: C in facets: fam, chr
#> Level: D in facets: fam, chr

#show or overwrite sample meta
sample.meta(test)
#>     pop fam .sample.id
#> 1   ASP   A          1
#> 2   ASP   B          2
#> 3   ASP   C          3
#> 4   ASP   D          4
#> 5   ASP   A          5
#> 6   ASP   B          6
#> 7   ASP   C          7
#> 8   ASP   D          8
#> 9   ASP   A          9
#> 10  ASP   B         10
#> 11  ASP   C         11
#> 12  CLF   D         12
#> 13  CLF   A         13
#> 14  CLF   B         14
#> 15  CLF   C         15
#> 16  CLF   D         16
#> 17  CLF   A         17
#> 18  OPL   B         18
#> 19  OPL   C         19
#> 20  OPL   D         20
#> 21  OPL   A         21
#> 22  OPL   B         22
#> 23  OPL   C         23
#> 24  OPL   D         24
#> 25  OPL   A         25
#> 26  OPL   B         26
#> 27  OPL   C         27
#> 28  OPL   D         28
#> 29  OPL   A         29
#> 30  OPL   B         30
#> 31  OPL   C         31
#> 32  OPL   D         32
#> 33  OPL   A         33
#> 34  OPL   B         34
#> 35  OPL   C         35
#> 36  OPL   D         36
#> 37  OPL   A         37
#> 38  OPL   B         38
#> 39  OPL   C         39
#> 40  OPL   D         40
#> 41  OPL   A         41
#> 42  OPL   B         42
#> 43  OPL   C         43
#> 44  OPL   D         44
#> 45  PAL   A         45
#> 46  PAL   B         46
#> 47  PAL   C         47
#> 48  PAL   D         48
#> 49  PAL   A         49
#> 50  PAL   B         50
#> 51  PAL   C         51
#> 52  PAL   D         52
#> 53  PAL   A         53
#> 54  PAL   B         54
#> 55  PAL   C         55
#> 56  PAL   D         56
#> 57  PAL   A         57
#> 58  PAL   B         58
#> 59  PAL   C         59
#> 60  PAL   D         60
#> 61  PAL   A         61
#> 62  SMR   B         62
#> 63  SMR   C         63
#> 64  SMR   D         64
#> 65  SMR   A         65
#> 66  SMR   B         66
#> 67  SMR   C         67
#> 68  SMR   D         68
#> 69  SMR   A         69
#> 70  SMR   B         70
#> 71  SMR   C         71
#> 72  SMR   D         72
#> 73  SMR   A         73
#> 74  SMR   B         74
#> 75  SMR   C         75
#> 76  SMR   D         76
#> 77  SMR   A         77
#> 78  SMR   B         78
#> 79  SMR   C         79
#> 80  SMR   D         80
#> 81  SMR   A         81
#> 82  SMR   B         82
#> 83  SMR   C         83
#> 84  SMR   D         84
#> 85  SMR   A         85
#> 86  UPD   B         86
#> 87  UPD   C         87
#> 88  UPD   D         88
#> 89  UPD   A         89
#> 90  UPD   B         90
#> 91  UPD   C         91
#> 92  UPD   D         92
#> 93  UPD   A         93
#> 94  UPD   B         94
#> 95  UPD   C         95
#> 96  UPD   D         96
#> 97  UPD   A         97
#> 98  UPD   B         98
#> 99  UPD   C         99
#> 100 UPD   D        100
sample.meta(test) <- data.frame(fam = sample(LETTERS[1:4], ncol(test), 
                                             replace = TRUE), 
                                pop = sample(LETTERS[5:8], ncol(test), 
                                             replace = TRUE))
#> Warning: Some levels are duplicated across multiple sample meta facets.
#> This will cause issues if those sample facets are run during analysis.
#> Issues:
#> 
#> Level: A in facets: fam, chr
#> Level: B in facets: fam, chr
#> Level: C in facets: fam, chr
#> Level: D in facets: fam, chr
```
