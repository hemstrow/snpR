# Merge two snpRdata objects

Merge two snpRdata objects using sample and SNP metadata. Functions much
like base R's [`merge`](https://rdrr.io/r/base/merge.html) function, but
the 'by' and 'all' options can be specified at the SNP and sample level.

## Usage

``` r
merge_snpRdata(
  x,
  y,
  by.sample = intersect(names(sample.meta(x)), names(sample.meta(y))),
  by.sample.x = by.sample,
  by.sample.y = by.sample,
  by.snp = intersect(names(snp.meta(x)), names(snp.meta(y))),
  by.snp.x = by.snp,
  by.snp.y = by.snp,
  all = TRUE,
  all.x.snps = all,
  all.y.snps = all,
  all.x.samples = all,
  all.y.samples = all,
  resolve_conflicts = "error"
)
```

## Arguments

- x, y:

  `snpRdata` objects to merge

- by.sample, by.sample.x, by.sample.y:

  Columns of sample metadata by which to merge across samples–function
  identically to the `by`, `by.x`, and `by.y` arguments to
  [`merge`](https://rdrr.io/r/base/merge.html), see documentation there
  for details.

- by.snp, by.snp.x, by.snp.y:

  Columns of SNP metadata by which to merge across SNPs–function
  idetically to the `by`, `by.x`, and `by.y` arguments to
  [`merge`](https://rdrr.io/r/base/merge.html), see documentation there
  for details.

- all:

  logical, default TRUE. If TRUE, all samples and SNPs will be
  maintained in the output `snpRdata` object, with missing data matching
  the missing data format of `x` added where genotypes are not in either
  `x` or `y`.

- all.x.snps, all.y.snps:

  logical, default `all`. Keep SNPs in the data even if they are only
  present in `x` or `y`, respectively.

- all.x.samples, all.y.samples:

  logical, default `all`. Keep samples in the data even if they are only
  present in `x` or `y`, respectively.

- resolve_conflicts:

  character, default 'error'. Controls how conflicting genotypic
  information in `x` and `y` is handled. See 'Details' for options and
  explanation.

## Value

A merged `snpRdata` object.

## Details

While this function can be used essentially identically to how one might
use base R's [`merge`](https://rdrr.io/r/base/merge.html) function,
there are a few differences to note.

First, samples that are genotyped at identical loci in both data sets
can be handled several ways, controlled by the `resolve_conflicts`
argument:

- warning:Return a harsh warning and a data frame with more information
  on genotypes at identical samples/SNPs are different between `x` and
  `y`.

- error: The default, return an error when conflicts are detected.

- x: Use genotypes from `x` to resolve conflicts.

- y: Use genotypes from `y` to resolve conflicts.

- random: Randomly sample (non-missing) genotypes from `x` and `y` to
  resolve conflicts.

Note that called genotypes are always taken over un-called genotypes
when there are merge conflicts, and missing data in one but not the
other data set will not trigger an error or a warning if those options
are selected.

Secondly, the `by` and `all` arugment families from
[`merge`](https://rdrr.io/r/base/merge.html) are extended to refer to
either samples or SNPs, such that all samples can be maintained but not
all SNPs, for example.

Lastly, all of the `all` family of arguments default to `TRUE` instead
of `FALSE`, since purely overlapping genotypes/SNPs is unlikely to be
desired. `FALSE` values provided to any specific `all` argument will
sill override `all = TRUE`, as in
[`merge`](https://rdrr.io/r/base/merge.html).

At present, `merge_snpRdata` is not maximally efficient in that it will
remove all tabulated statistics and re-tabulate all internal summaries.
Improvements are in development.

## Author

William Hemstrom

## Examples

``` r
# create data to merge in
y <- data.frame(s1 = c("GG", "NN"),
                s2 = c("GG", "TG"),
                s3 = c("GG", "TT"),
                s4 = c("GA", "TT"),
                s5 = c("GG", "GT"),
                s6 = c("NN", "GG"))
                
snp.y <- data.frame(chr = c("groupVI", "test_chr"),
                    position = c(212436, 10))
                   
samp.y <- data.frame(pop = c("ASP", "ASP", "ASP", "test1", "test2", "test3"),
                     ID = c(1, 2, 3, "A1", "A2", "A3"),
                     fam = c("A", "B", "C", "T", "T", "T"))
y <- import.snpR.data(y, snp.y, samp.y)

x <- stickSNPs
sample.meta(x)$ID <- 1:ncol(x)

if (FALSE) { # \dontrun{
# Not run, will error due to conflicts
z <- merge_snpRdata(x, y)

# Not run, will return a warning and report mismatches
z <- merge_snpRdata(x, y, resolve_conflicts = "warning")
} # }

# take a random genotype in the case of conflicts
z <- merge_snpRdata(x, y, resolve_conflicts = "random")
#> Warning: Some levels are duplicated across multiple sample meta facets.
#> This will cause issues if those sample facets are run during analysis.
#> Issues:
#> 
#> Level: 10    in facets: ID, position
z
#>     V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 V20 V21
#> 1   NN TG TT NN NN NN NN NN NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN
#> 2   AG GG AA AG AA AG AA NN AG  NN  GG  AA  GG  AG  NN  AG  NN  AG  AG  AA  AG
#> 3   GG GG GG GG GG GG GG GT GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 4   CC CC CC CC CC CC CC CC CC  CT  CC  CC  CT  CC  CC  CC  NN  CC  CC  CT  CC
#> 5   AA AA AA AA NN AA AA AG AA  AA  AA  AA  AG  AG  AA  AG  AG  AG  AG  AA  AG
#> 6   GG NN NN NN NN GG AA GG AG  AG  AG  NN  GG  NN  GG  NN  GG  NN  GG  GG  AA
#> 7   AC AC AA AC AA AC AA AA AC  AA  AC  AA  AC  NN  AA  AC  AA  AA  AA  AA  AA
#> 8   TT TT TT TT TT TT TT GT TT  TT  TT  TT  TT  TT  GG  TT  TT  TT  TT  TT  TT
#> 9   GG GG AG AG AA AG AG AG GG  GG  AG  GG  AG  GG  GG  AG  AG  NN  AG  GG  GG
#> 10  CC AC CC CC AC CC AC CC CC  CC  CC  AC  CC  AA  CC  AC  CC  NN  CC  CC  CC
#> 11  AA AA AG AA AA AG AA AA AA  AA  AA  AG  AA  AA  AA  AA  NN  AA  AA  AA  AA
#> 12  NN TT TT CT CT TT CT TT CT  CC  CC  CT  CT  CT  CT  CT  CC  NN  CT  CC  CC
#> 13  AG AG AG AG AA AA NN AA AA  AG  AA  AG  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 14  TT TT CT TT CT TT TT CT NN  TT  TT  NN  TT  TT  CT  CT  TT  TT  TT  CC  NN
#> 15  GG GT GT GG GT GT NN GG GG  GG  GG  GT  GG  GT  GG  GG  NN  GG  GG  GG  GG
#> 16  AA AG AG AA AA AA AA AA AA  AA  AA  AG  AG  AG  AA  AA  NN  NN  AA  AA  NN
#> 17  NN GG GT GG NN GG TT GG GT  GT  GG  GG  GG  GG  GT  GT  GT  TT  GG  GT  NN
#> 18  AA AA AA AA AA AA AA AA NN  AC  AC  AA  AC  AA  AA  AA  AA  AA  AC  AA  AC
#> 19  AC NN CC AC NN AC CC AC AC  CC  AC  AC  CC  AC  CC  AC  CC  AC  AC  AC  CC
#> 20  TT TT TT TT GT TT TT TT TT  TT  TT  TT  TT  TT  GT  TT  TT  TT  NN  GT  NN
#> 21  GG AG GG GG GG GG GG NN GG  GG  GG  GG  NN  GG  GG  GG  GG  NN  AG  AG  GG
#> 22  NN NN CC AC NN AC CC AC AC  CC  CC  AC  CC  CC  NN  AC  AC  AC  AC  CC  NN
#> 23  CT NN TT NN TT CT CC CT NN  CC  TT  CT  CC  CT  CC  CC  NN  NN  CC  CT  CC
#> 24  GG GG GG GG GG GG GG GG GG  GG  AG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG
#> 25  AA AA AA AA AA AA NN AA AA  AA  AA  AA  AA  AG  AA  AA  NN  AA  AG  AA  AA
#> 26  NN AG AG AG NN GG GG GG AG  AG  GG  AG  AG  GG  GG  GG  NN  AG  AG  GG  AG
#> 27  CT CT CT TT TT TT CT NN TT  NN  TT  CT  TT  CT  TT  TT  TT  NN  TT  TT  CT
#> 28  AA AA AA AA AA AA AA AA AA  AA  AA  NN  AC  AA  AA  AC  NN  AC  AA  AA  AC
#> 29  TT NN NN CT TT TT NN CT CT  TT  TT  TT  TT  TT  CT  TT  TT  CT  NN  TT  NN
#> 30  AA AC AA NN AA AC NN AA AA  AA  NN  AA  AA  AA  AA  NN  NN  AA  AA  AA  AA
#> 31  GG GG CG CG GG CG GG GG GG  GG  GG  GG  GG  CG  GG  GG  CG  NN  GG  NN  GG
#> 32  CG CG CC CC CC CC CC CC CC  CG  CC  CG  CG  CG  CC  CG  NN  CG  CG  CC  CG
#> 33  AA AA AA NN AA AA NN AA AA  AA  AC  AC  AC  AA  AC  AC  AC  AC  AC  AA  AA
#> 34  CC CT CC CT NN CT CC CT CT  CT  CC  CC  CT  CT  CT  CC  CC  CT  CT  NN  NN
#> 35  CG CG CG GG CG GG CG CC NN  CC  CG  CG  CG  NN  CC  CG  NN  CG  CG  CC  NN
#> 36  AA AA AA AA AA AA AC AA NN  AC  AA  AA  AA  AA  AA  NN  AA  AA  AC  AA  AA
#> 37  AA AA AA AA AG AA AA AA AG  AA  AG  AA  AG  NN  AG  AA  AA  AA  AA  AA  NN
#> 38  CG GG GG CG NN GG GG GG GG  NN  GG  CC  GG  CG  GG  CG  NN  GG  CG  GG  GG
#> 39  GG TT GG GG GG GG GT GT TT  GG  GT  NN  GG  TT  GG  GG  GT  GT  GG  GG  GT
#> 40  TT NN TT NN AT NN NN AT TT  TT  TT  TT  AT  AT  TT  NN  NN  AT  TT  AT  NN
#> 41  CG CC CG NN CG GG NN CG GG  CG  GG  CC  GG  CG  CG  GG  CG  GG  CG  GG  CG
#> 42  AG AA AG GG AG NN AG AA AA  AG  AG  GG  AG  AG  AG  GG  GG  GG  AG  AG  AG
#> 43  AG NN NN AG GG AG NN GG GG  AG  AG  AG  GG  GG  NN  GG  GG  AG  GG  AA  GG
#> 44  AG AG AA AA AA AG AA AA AA  AA  AA  AA  AA  AA  AG  NN  AG  AA  AA  AA  NN
#> 45  AT TT TT TT TT TT TT NN TT  TT  NN  TT  AA  TT  TT  TT  TT  TT  TT  TT  TT
#> 46  GG GG GG GG GG GG GG GG GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN
#> 47  NN CC AA CC NN CC AA NN CC  NN  AA  CC  CC  AC  AC  CC  NN  AC  CC  AC  AC
#> 48  CC CC CC CC CC CC CG CG CC  CG  CG  CC  CG  CC  CC  CC  CC  CG  CG  CC  NN
#> 49  CT TT TT TT CT CT NN TT TT  CT  CC  NN  TT  CT  CT  CT  TT  CT  CT  CC  CT
#> 50  GG GG GG CG GG CG GG NN GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 51  NN CC CT CC CC CT CC CC CT  CC  CC  CT  CC  CC  CC  CC  CC  NN  CC  NN  CT
#> 52  TT NN CT TT TT TT CT TT CT  CT  NN  CT  CC  TT  TT  NN  CT  TT  TT  TT  CT
#> 53  AA AG AG AA AA AA AA AG AG  AA  AG  AG  AG  AG  GG  AG  AA  AG  AG  AG  AG
#> 54  AA AA CC AA NN AC NN AA AC  AA  AA  AC  AA  AA  AC  NN  AA  AC  AC  AC  NN
#> 55  GG GG GG GG GT GG GG GG GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN
#> 56  GG NN AG GG AG AG NN AG AG  AG  AG  AG  AA  NN  GG  GG  AG  AG  AG  AG  AA
#> 57  AG AG AG AG AA AA AG AA NN  AG  AA  AG  AA  NN  AA  AA  AA  AA  AA  AA  NN
#> 58  GG NN GG GG GG GG NN GG GG  GG  AG  NN  GG  GG  NN  GG  GG  GG  GG  GG  AG
#> 59  TT CT CC CC CC CC NN TT CC  NN  CT  CC  CT  CT  CC  CC  CC  CC  CC  CC  CT
#> 60  AA AG GG AG AG GG AG GG AA  GG  AG  GG  AA  AG  AG  AG  AG  NN  AG  GG  AG
#> 61  AA AA AC CC CC CC NN AC AC  AC  NN  NN  CC  NN  AC  AC  AA  AC  AA  AC  AC
#> 62  TT NN TT NN CT TT NN TT TT  TT  NN  TT  TT  TT  TT  NN  NN  NN  TT  NN  TT
#> 63  CG GG GG GG GG GG NN GG GG  GG  GG  CG  GG  GG  GG  GG  NN  GG  CG  GG  GG
#> 64  AC NN NN AC AA NN NN AA AC  AC  AC  CC  AA  AA  CC  AC  AA  AC  AC  AA  NN
#> 65  TT AT AT TT TT TT TT AT AT  AA  AT  TT  TT  TT  TT  TT  AT  AT  TT  AA  AT
#> 66  CT CC CT CT CT CT CC CC CT  CT  CC  CC  CC  CC  TT  CC  CT  CC  TT  CC  NN
#> 67  AC AC NN CC AA AC NN AC CC  CC  AC  AA  AA  AC  AC  AC  NN  NN  AA  AA  NN
#> 68  GG GG GG GG GG GG NN GG GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 69  AA AT AT AA AT AA AA AA NN  AT  TT  AA  AA  AA  AT  AA  AT  AA  AT  AT  AT
#> 70  AA AA NN AA AA AG AA AA AA  AA  AA  AA  AA  AA  NN  AA  AA  AA  AA  AA  GG
#> 71  AC AC AA AC NN AC AC AC NN  AC  AC  AA  AA  AC  CC  AA  AA  NN  AA  CC  NN
#> 72  AA AA AA AA AA AA AT AA AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 73  GG GG GG AG GG GG GG AG NN  AG  AG  AA  GG  AA  AG  AA  NN  AA  AA  AG  NN
#> 74  TT GT TT GT TT GT TT GT TT  GT  NN  GT  GT  GT  GT  TT  TT  TT  TT  TT  NN
#> 75  GT GG GT GG GG NN GG GG GG  GT  GG  GG  GG  NN  GG  GT  NN  GT  GG  GG  GT
#> 76  AC AC NN AC CC AC AC CC AC  AC  AC  AA  CC  NN  AA  AC  AC  AA  AC  NN  AC
#> 77  CC CT CT NN NN CC CT NN NN  TT  NN  NN  CC  NN  CT  NN  NN  CC  TT  CC  NN
#> 78  AT AT AT AT AA AA AT AA AA  AT  NN  AT  NN  AA  AA  AA  AA  AA  AA  AA  AA
#> 79  GG GG GG GG GG GT GG GG GG  GG  GG  NN  GG  GG  TT  GG  GG  GG  GG  GG  NN
#> 80  TT TT CT TT TT TT TT CT CT  CT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 81  CT CT TT TT TT TT TT TT CT  TT  TT  CC  TT  TT  CT  CT  TT  CT  TT  TT  CT
#> 82  AA NN AT AT AT AT AT AT AT  AT  AA  AT  NN  AA  TT  NN  AT  AT  AT  AT  AA
#> 83  CT TT CT NN TT NN TT TT TT  NN  TT  TT  TT  TT  TT  NN  TT  TT  TT  TT  TT
#> 84  AT TT AT TT AT AT TT TT AT  TT  TT  TT  AT  TT  AT  AT  TT  TT  TT  AT  TT
#> 85  CC CC CT CC CC CT CT CC NN  CT  CC  CT  CC  TT  TT  CT  CT  NN  CC  CC  CC
#> 86  AG AG GG NN AG AG AA AG NN  AA  AA  AG  AA  AG  AG  AA  AG  GG  AG  AG  NN
#> 87  AG AG AG AG AA AA AG AA AA  AG  NN  AG  NN  NN  AA  AA  AA  NN  AA  AA  NN
#> 88  CT CC CT TT TT TT CT TT CT  CC  CT  TT  TT  CC  TT  CT  TT  TT  TT  TT  CT
#> 89  AG NN AG AG AA AG NN GG AG  AA  AG  AG  AG  NN  AA  AA  AA  AA  AA  AG  GG
#> 90  NN CG CG GG CC NN CG CC CG  CG  NN  CC  CG  CC  CC  NN  NN  CG  CG  NN  CG
#> 91  CT CT CC CC CC CC CC CC CT  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 92  CC CC NN TT CC NN CT CC CT  CC  CT  CC  NN  CC  CC  NN  CT  CT  NN  CT  CT
#> 93  AA NN GG GG GG GG NN GG AA  AG  AG  NN  NN  AG  AG  NN  NN  NN  GG  AG  AG
#> 94  TT TT TT CT TT CT TT CT TT  TT  TT  TT  CT  CT  CT  TT  TT  TT  TT  TT  TT
#> 95  CC CC CC CC CC NN CC CC NN  CC  CC  CC  NN  NN  CC  CC  CC  CC  CC  CC  CG
#> 96  GG GG GG GG GG GT NN GT GG  GT  GT  GG  GT  GG  GT  TT  GG  GG  GG  GG  GG
#> 97  CC NN CC CC AC CC NN CC CC  CC  CC  AC  CC  NN  CC  CC  CC  CC  AC  CC  AC
#> 98  GG GG GG GG GG AG GG GG GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  AG
#> 99  GG GG GG GG NN GT NN GG GT  GG  GT  NN  GG  TT  GG  GG  NN  GG  GG  TT  GG
#> 100 AA GG AA GG AA GG AA AG AG  AG  AG  NN  GG  AG  AG  AA  AG  AA  AA  GG  NN
#> 101 AT NN AT NN AA NN AT TT AT  TT  TT  TT  AA  AA  AT  AT  AA  NN  AT  AT  NN
#>     V22 V23 V24 V25 V26 V27 V28 V29 V30 V31 V32 V33 V34 V35 V36 V37 V38 V39 V40
#> 1    NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN
#> 2    AG  AG  AG  NN  NN  NN  AG  AG  AG  AG  AA  AA  AA  NN  AG  GG  AA  AG  AA
#> 3    NN  GT  GG  GG  GG  GT  GG  GG  GG  GG  GG  GT  GG  NN  GG  NN  GG  GG  GG
#> 4    CT  CT  CC  CC  CT  CC  NN  CT  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 5    AG  NN  AG  NN  NN  NN  AG  NN  AG  NN  AG  AA  AG  AG  AG  AG  AG  AG  NN
#> 6    NN  NN  AG  AG  NN  GG  GG  AG  NN  NN  NN  NN  NN  AG  NN  NN  GG  AG  GG
#> 7    AA  AC  AA  AA  AA  NN  AA  AA  AA  AA  AA  AA  AA  AA  AC  AA  AA  AA  AA
#> 8    NN  GT  TT  GT  TT  NN  NN  TT  GT  TT  TT  TT  NN  TT  TT  NN  GT  TT  GT
#> 9    GG  GG  GG  AG  GG  AG  NN  GG  GG  GG  GG  GG  GG  AG  NN  AG  GG  AA  GG
#> 10   AC  CC  CC  AC  NN  CC  CC  CC  CC  AC  CC  CC  CC  CC  CC  CC  AC  AC  CC
#> 11   AA  AA  AA  AA  AA  AA  AA  NN  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 12   NN  NN  CC  CT  CC  CC  NN  NN  CT  CT  NN  TT  CT  CT  CT  CC  TT  TT  CT
#> 13   NN  NN  AA  AG  AG  AG  NN  NN  AA  NN  AA  NN  AG  AG  AA  AG  AG  AA  AA
#> 14   TT  NN  TT  TT  TT  CT  NN  TT  TT  TT  CT  CC  TT  CT  CT  CC  TT  CT  CT
#> 15   GG  GG  GG  GG  GG  GG  TT  GG  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG
#> 16   AG  AA  AA  AA  AA  NN  AA  AG  NN  AA  AG  AG  AA  AG  AA  NN  NN  AA  AA
#> 17   GT  NN  NN  NN  GT  NN  GT  GG  GG  GT  GT  GG  GG  GT  GT  GG  GT  GT  GG
#> 18   AA  AC  NN  AC  NN  AA  AC  AC  AC  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 19   CC  CC  CC  CC  CC  CC  AC  CC  AC  CC  AC  AA  AC  CC  CC  CC  CC  CC  AC
#> 20   TT  TT  GT  NN  TT  NN  TT  NN  TT  TT  TT  NN  GT  TT  TT  GT  TT  TT  TT
#> 21   GG  NN  GG  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 22   AC  AA  NN  NN  AC  NN  AA  AC  NN  AC  NN  AC  NN  CC  CC  AC  AC  NN  CC
#> 23   CT  CT  NN  CC  CC  NN  CC  CT  CT  CC  CT  CT  CC  NN  TT  CT  CC  CT  CT
#> 24   GG  AG  GG  GG  GG  GG  GG  AG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG
#> 25   AA  AA  AG  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AG  AA  AA  AA  AA  AA
#> 26   NN  AG  GG  GG  GG  AG  AG  AG  AG  GG  NN  GG  GG  AG  GG  GG  GG  GG  GG
#> 27   TT  TT  NN  NN  TT  TT  TT  NN  TT  TT  CT  TT  TT  TT  TT  CT  TT  NN  TT
#> 28   AC  AA  AC  AA  AA  NN  AC  AA  AA  AA  AC  AA  AC  AC  AC  AA  NN  AA  AA
#> 29   CT  CT  NN  NN  NN  CT  CT  CC  CT  NN  CC  CC  CC  CC  TT  CT  CC  NN  CT
#> 30   AA  NN  AA  NN  AC  NN  NN  AA  CC  NN  AA  AC  AA  AA  AA  AC  AC  AC  AA
#> 31   NN  GG  NN  GG  CG  NN  CG  GG  GG  NN  NN  GG  GG  GG  GG  GG  GG  GG  GG
#> 32   CC  CG  CC  CG  CC  CC  NN  CG  CC  CC  CG  CC  CG  GG  CC  CC  CG  CG  CG
#> 33   AA  AC  NN  NN  AA  AC  AA  AC  AA  AA  AA  NN  NN  AA  AA  AA  AA  AA  AA
#> 34   NN  CC  CT  NN  NN  TT  TT  CC  NN  CT  CC  CT  CC  CC  NN  CT  NN  NN  CC
#> 35   NN  CG  NN  CG  NN  NN  CG  CG  GG  CC  NN  CC  CC  NN  CC  CC  CG  CG  CG
#> 36   AA  AA  NN  NN  AA  AA  AC  AA  NN  AA  AA  AA  NN  NN  AA  AC  AA  AA  AA
#> 37   AG  AA  AA  AA  NN  AG  AA  NN  AA  AA  AA  AA  AG  AA  AA  AG  AA  AG  GG
#> 38   CG  GG  CG  CG  CG  CG  GG  CG  GG  CG  CG  GG  GG  GG  CG  GG  GG  GG  CG
#> 39   GG  NN  NN  NN  GT  GT  NN  GG  GG  GG  GT  GG  GG  GT  GG  GG  GG  GG  GT
#> 40   NN  NN  TT  NN  NN  AT  NN  NN  TT  TT  TT  TT  TT  NN  TT  NN  AT  TT  TT
#> 41   NN  GG  GG  NN  NN  CG  NN  CG  GG  NN  GG  NN  NN  GG  GG  GG  GG  GG  GG
#> 42   NN  NN  NN  AG  NN  AG  NN  NN  AG  AG  NN  GG  AG  GG  GG  GG  GG  NN  AG
#> 43   AG  AG  AG  AG  GG  AG  AA  NN  NN  AG  AG  AA  NN  AG  AG  AA  AA  AG  AA
#> 44   AA  AA  NN  NN  AG  AG  NN  NN  AG  AA  AG  AA  AA  AA  AA  AA  NN  AA  AG
#> 45   TT  TT  NN  TT  TT  TT  AT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 46   GG  NN  GG  NN  GG  GG  GG  GG  GG  GG  NN  NN  GG  GG  NN  GG  GG  NN  GG
#> 47   AC  NN  NN  NN  NN  AA  CC  CC  CC  NN  AC  AC  NN  NN  AC  CC  NN  AC  AC
#> 48   CG  CC  CG  CC  CG  CG  CG  CG  CG  GG  CG  CC  CG  NN  GG  CG  NN  CG  GG
#> 49   CC  NN  CT  NN  NN  NN  NN  CC  CC  CT  NN  CT  TT  CT  CT  CC  CC  CT  CT
#> 50   GG  GG  GG  GG  GG  GG  GG  CG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 51   CC  CC  CT  NN  NN  CC  CC  CC  CC  NN  CC  CC  CC  CC  CC  CT  CC  CC  CC
#> 52   CC  NN  TT  NN  NN  CT  CT  CT  NN  CT  CT  TT  NN  CT  CT  TT  CC  NN  NN
#> 53   AG  AA  NN  AG  AG  AA  AA  AG  GG  AG  NN  AA  NN  AG  AA  AA  AG  AG  AG
#> 54   AA  AA  AC  AC  NN  AC  AA  NN  AC  NN  AC  CC  AC  AA  AC  AC  AC  AA  AA
#> 55   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GT
#> 56   NN  NN  NN  AG  GG  AG  AA  GG  AG  AA  NN  GG  NN  GG  GG  GG  AG  NN  GG
#> 57   AA  NN  NN  AG  AG  AG  AG  AA  AA  NN  AA  NN  AG  AG  AA  AG  AG  AA  AA
#> 58   NN  NN  NN  GG  GG  GG  NN  GG  GG  NN  NN  AG  AG  AG  GG  GG  GG  GG  GG
#> 59   CT  CC  CC  CC  CT  CT  CT  CC  CC  CT  CT  CT  CT  CT  NN  CC  CC  CT  CT
#> 60   AG  AG  AG  AG  AA  GG  GG  AA  GG  GG  GG  AG  NN  AG  NN  AG  AG  GG  AG
#> 61   AC  AC  NN  NN  AC  AC  AC  AA  AC  AA  AA  CC  AA  NN  AC  AC  AC  AA  AA
#> 62   NN  NN  TT  TT  NN  TT  TT  TT  NN  TT  TT  TT  CT  TT  TT  TT  NN  TT  NN
#> 63   GG  CG  GG  CG  GG  GG  GG  CG  CG  CG  GG  CG  GG  CG  CG  GG  CG  CG  GG
#> 64   AC  NN  AA  AC  AA  AA  NN  AC  AA  NN  AC  AA  NN  AA  AC  NN  AA  NN  AA
#> 65   NN  AT  NN  AT  TT  AT  NN  AT  AT  TT  NN  TT  TT  NN  AT  TT  AT  NN  AT
#> 66   NN  TT  NN  CT  TT  CT  NN  NN  TT  CC  CT  CT  TT  NN  CT  TT  TT  NN  CT
#> 67   AC  AC  AC  AA  AC  NN  AA  AA  AC  NN  CC  AA  AC  AA  AC  AC  AA  AC  AC
#> 68   NN  GG  GG  GG  GG  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 69   AT  AA  AT  NN  AA  NN  AT  AA  AA  AA  AA  AA  AA  AT  AT  AA  AA  AA  AT
#> 70   AA  AA  NN  AG  AG  AA  AA  AA  AG  AA  AA  AG  AG  AG  AA  AG  GG  AA  AG
#> 71   AA  NN  NN  AC  AC  AC  AA  AA  AC  AA  AA  AC  AC  AC  AA  AC  NN  AA  AC
#> 72   AT  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 73   GG  AG  AG  AG  NN  AA  AG  AA  AA  AA  AG  AG  AA  AG  AA  AG  AA  AA  AG
#> 74   TT  TT  TT  TT  TT  TT  NN  GT  TT  GT  NN  GT  TT  GT  TT  NN  TT  GT  TT
#> 75   GT  GG  NN  GT  NN  GG  NN  GT  GT  GG  GT  GT  GG  NN  GT  GG  GT  NN  NN
#> 76   AC  AC  AC  NN  AC  CC  NN  NN  CC  AC  AC  AA  AA  AC  AA  NN  AA  AA  AC
#> 77   CT  NN  CC  NN  NN  NN  TT  CT  CC  CC  CC  TT  NN  CT  NN  CC  NN  NN  CC
#> 78   AA  AA  NN  AT  AT  AT  NN  AA  AA  AT  NN  AA  AT  AT  AA  AT  AT  NN  AA
#> 79   GG  NN  NN  GT  GT  NN  NN  GG  GT  GG  GG  GT  GT  GT  GG  GT  TT  GG  GT
#> 80   TT  TT  TT  TT  TT  CT  NN  CT  CT  TT  CT  TT  TT  TT  TT  TT  TT  TT  CT
#> 81   CT  CT  NN  NN  CT  TT  CT  NN  CT  CT  CT  TT  NN  NN  TT  CC  TT  TT  TT
#> 82   NN  AA  NN  NN  NN  AT  AT  AA  AT  AA  AT  TT  AT  AA  AT  NN  AT  AT  AT
#> 83   TT  TT  TT  CT  NN  TT  TT  CT  TT  TT  TT  CT  TT  TT  TT  TT  CT  TT  CT
#> 84   TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 85   CT  NN  NN  CC  CC  NN  CT  NN  CC  NN  CT  CT  CC  NN  TT  CT  CC  CT  TT
#> 86   GG  NN  AG  NN  AG  GG  AG  AG  GG  GG  NN  NN  AG  AG  AG  GG  AG  AG  AG
#> 87   AA  AA  NN  AG  AG  AG  AG  NN  NN  AG  AA  AA  AG  AG  AA  AG  AG  AA  NN
#> 88   NN  TT  CT  CT  TT  TT  TT  NN  TT  CC  TT  TT  TT  CT  CT  TT  CT  TT  TT
#> 89   NN  NN  NN  AG  AG  AA  AG  NN  AG  NN  AG  NN  AG  AA  AA  AA  AG  NN  AG
#> 90   NN  CC  CG  NN  CG  NN  NN  NN  CC  CG  NN  GG  CC  NN  CC  CG  CC  CC  CC
#> 91   CC  CC  NN  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 92   NN  CT  CT  CT  NN  NN  CT  CC  CC  CT  CT  CC  TT  TT  CT  NN  TT  CT  CT
#> 93   NN  AG  AG  NN  AG  NN  GG  NN  NN  NN  AG  NN  NN  AG  AG  AG  GG  NN  AG
#> 94   TT  TT  TT  CT  NN  TT  TT  CT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 95   CC  CC  NN  CC  CC  CG  CC  CC  CC  CC  NN  CC  CC  NN  CC  CC  NN  CC  CC
#> 96   NN  NN  GG  GT  NN  NN  GT  NN  GT  GG  TT  GG  GG  GG  GG  GG  TT  GG  GG
#> 97   CC  AC  NN  NN  NN  CC  CC  CC  CC  NN  AC  CC  CC  CC  CC  CC  AC  CC  AC
#> 98   GG  GG  NN  NN  GG  NN  GG  NN  GG  GG  GG  GG  GG  GG  AG  NN  GG  GG  GG
#> 99   GG  GT  GG  GG  GG  GT  GT  GG  GG  GG  GG  GG  GT  GG  GG  GG  GG  GG  GG
#> 100  NN  AG  AA  AG  AG  AG  NN  NN  GG  NN  AG  AG  GG  NN  AA  AG  AG  AG  AG
#> 101  NN  AT  AT  AT  AT  NN  NN  NN  TT  AA  AA  AA  NN  TT  AA  AT  AA  AT  AA
#>     V41 V42 V43 V44 V45 V46 V47 V48 V49 V50 V51 V52 V53 V54 V55 V56 V57 V58 V59
#> 1    NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN
#> 2    NN  NN  GG  GG  AG  GG  AG  AG  GG  GG  GG  GG  GG  AG  GG  NN  AG  AA  GG
#> 3    NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 4    CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  NN  CC
#> 5    AA  AG  AA  AA  AA  AG  AA  AA  AA  AA  AA  AA  AA  AA  NN  NN  AG  AA  AA
#> 6    NN  NN  AG  NN  AG  AG  NN  AG  NN  GG  GG  GG  GG  AG  NN  GG  NN  NN  GG
#> 7    AC  AA  AA  AC  AC  AA  AA  AA  AC  AC  NN  AC  AC  AC  AC  AA  AA  NN  AC
#> 8    GT  TT  TT  TT  TT  TT  NN  TT  TT  TT  TT  TT  TT  TT  TT  GT  GT  TT  TT
#> 9    GG  GG  GG  NN  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  AG  AG  GG
#> 10   CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  AC  CC  CC  CC  CC  AC  AC  AC  CC
#> 11   AA  AG  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  NN  AA  AG  AA
#> 12   CT  CC  TT  TT  TT  CT  CT  TT  TT  TT  TT  TT  TT  TT  NN  TT  CT  NN  TT
#> 13   NN  NN  AA  AG  AG  AA  NN  AG  AG  AG  AA  AG  AG  AG  NN  AA  AA  NN  AA
#> 14   TT  CT  TT  NN  TT  TT  CT  CT  TT  TT  TT  TT  NN  TT  TT  TT  TT  CT  TT
#> 15   GG  GG  GG  GG  GG  GG  TT  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GT  GG
#> 16   AG  GG  NN  GG  AG  AA  AA  GG  GG  AG  AG  GG  GG  GG  GG  AA  AG  AA  GG
#> 17   GT  NN  TT  GG  GT  GG  NN  GT  NN  TT  NN  GT  GT  GT  GG  NN  GT  GT  GT
#> 18   AA  AA  NN  AC  AA  AA  AA  CC  AA  AA  AA  AA  AC  AA  AA  NN  AA  AA  AA
#> 19   CC  CC  AC  CC  NN  AC  AC  NN  AA  AC  AA  AA  AC  AA  AC  NN  CC  CC  AA
#> 20   TT  NN  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 21   GG  GG  GG  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN
#> 22   AC  AC  CC  CC  CC  AC  AA  CC  CC  CC  CC  CC  CC  CC  NN  AC  CC  CC  CC
#> 23   TT  CC  CC  CC  CC  CC  NN  CC  CC  CC  CC  CC  CC  CC  CC  NN  NN  CT  CC
#> 24   GG  GG  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  NN  GG  GG  AG  GG  GG
#> 25   AA  AA  NN  AA  AA  GG  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 26   AG  NN  GG  GG  NN  AG  AG  GG  GG  GG  GG  GG  GG  GG  GG  AG  NN  GG  GG
#> 27   CT  CT  CC  TT  CC  CT  NN  TT  CT  CT  CC  CC  CT  TT  CT  TT  CT  CT  NN
#> 28   AC  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AC  NN  AC  AA
#> 29   NN  CC  CT  TT  TT  CT  CC  TT  TT  TT  TT  TT  TT  CT  TT  CT  CC  NN  TT
#> 30   AA  CC  AA  AC  AA  AC  AC  AC  AA  AC  AC  AA  AA  AA  AA  AC  AA  NN  AA
#> 31   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  CG  GG  GG  GG
#> 32   CC  CC  CC  CC  CC  GG  CG  CC  CC  CC  CC  CC  CC  CC  CC  NN  CG  CC  CC
#> 33   AA  AA  AA  AC  AC  AA  AA  AA  AC  AA  AC  AC  AC  AC  AC  AA  AA  AC  AC
#> 34   TT  CC  CC  CC  CT  CT  CT  CC  CT  CC  NN  CC  CT  CC  CT  CC  CT  CC  CT
#> 35   CG  GG  CC  CC  CC  NN  CC  CC  CC  CC  CC  CC  CC  CC  CC  CG  NN  CC  CC
#> 36   AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  NN  AA  AC  NN
#> 37   AG  AA  AA  AA  AA  AA  AG  AA  AA  AA  AA  AA  AA  AA  AA  AG  AA  AA  AA
#> 38   CG  NN  GG  GG  GG  CG  NN  GG  GG  GG  GG  NN  GG  GG  CG  CG  CG  CG  GG
#> 39   GT  GG  GT  GT  GG  GG  GG  GG  GT  GT  GT  GT  GT  GT  GG  NN  NN  NN  GT
#> 40   NN  NN  NN  TT  TT  TT  TT  TT  AT  TT  NN  NN  TT  AT  NN  TT  TT  NN  NN
#> 41   GG  CG  GG  CG  CG  GG  CG  GG  CG  GG  CG  GG  CG  GG  NN  GG  CG  NN  GG
#> 42   AG  GG  AA  NN  NN  GG  AG  AA  AA  AA  AA  AA  AA  NN  AA  AG  NN  NN  AA
#> 43   GG  GG  GG  GG  AG  NN  GG  AG  GG  GG  GG  GG  GG  GG  NN  AG  GG  NN  GG
#> 44   AA  AG  AA  AA  AA  AA  AG  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AG  AA
#> 45   TT  AT  AT  AA  NN  TT  AT  AT  AT  AA  AT  AT  AT  AT  TT  AA  NN  TT  TT
#> 46   GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 47   AC  AA  CC  CC  CC  CC  AA  CC  CC  CC  CC  CC  CC  CC  CC  CC  AC  CC  CC
#> 48   CG  GG  NN  CG  CC  CC  NN  CC  CC  CC  CC  CC  NN  CG  CG  CC  CG  CG  CC
#> 49   TT  CC  CC  TT  CT  CT  NN  CT  CC  CT  CT  TT  CT  TT  NN  CT  TT  CT  NN
#> 50   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 51   CT  CC  CC  CC  CC  CC  CT  CC  CT  CC  CC  CC  CC  CC  CC  CT  NN  CC  CC
#> 52   TT  TT  TT  CT  TT  TT  CT  NN  TT  CT  CT  NN  NN  TT  TT  NN  NN  CT  NN
#> 53   AA  AA  AG  NN  NN  GG  AA  AG  AG  NN  AA  GG  AA  AG  NN  AA  AA  GG  AG
#> 54   NN  NN  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AC  AA
#> 55   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 56   GG  AG  GG  AG  GG  GG  GG  AA  AG  AG  AG  NN  AG  GG  AG  AG  GG  GG  AG
#> 57   AA  AA  AA  AG  AG  NN  NN  AG  AG  AG  AA  AG  AG  AG  AG  AA  NN  AA  AA
#> 58   GG  NN  GG  GG  GG  GG  GG  NN  NN  GG  GG  GG  GG  NN  GG  GG  GG  NN  GG
#> 59   NN  CC  CT  NN  CC  CT  CC  CC  TT  TT  CT  NN  CT  CT  CT  CC  NN  CC  CC
#> 60   AG  AG  GG  GG  GG  GG  AG  GG  GG  GG  GG  GG  GG  GG  GG  AG  AG  GG  GG
#> 61   AC  NN  AA  AA  AA  AC  AC  AA  AA  AA  AA  AA  AA  AA  AA  NN  AA  AC  AA
#> 62   TT  NN  NN  TT  TT  TT  NN  NN  TT  TT  NN  NN  TT  TT  TT  TT  TT  NN  TT
#> 63   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  CG  GG  GG  GG
#> 64   NN  AC  NN  AC  CC  AA  NN  NN  NN  AC  CC  AC  CC  NN  NN  AC  AC  NN  CC
#> 65   AT  TT  AA  TT  AA  TT  AT  AT  AA  NN  AT  AA  AT  AA  AT  TT  AT  TT  TT
#> 66   CT  NN  CT  CT  CC  TT  CC  CC  CC  CT  CC  CC  CC  CC  CC  NN  CT  CC  CT
#> 67   NN  NN  CC  CC  CC  AC  AC  CC  AC  AC  AA  AA  CC  AC  CC  AC  NN  NN  AC
#> 68   GG  GG  GG  AG  GG  GG  GG  AG  AA  GG  NN  GG  AG  GG  AG  GG  GG  GG  NN
#> 69   NN  TT  TT  TT  NN  AA  AA  TT  TT  NN  NN  TT  AT  AT  AA  AT  AT  AA  AT
#> 70   AG  NN  AA  AA  AA  GG  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AG  AA  AA
#> 71   AC  AC  AA  AA  AA  NN  AA  AA  AA  AA  AA  AA  NN  AA  AA  AC  AC  AA  AA
#> 72   AA  AA  AA  AA  AA  AA  AA  AA  AT  AA  AT  NN  AA  AA  AT  AA  AA  AT  AA
#> 73   AG  AG  GG  GG  NN  AG  AG  GG  AG  AG  GG  GG  AG  AA  NN  AA  NN  AA  GG
#> 74   TT  TT  NN  GT  TT  TT  GT  TT  TT  GT  GT  TT  GT  GT  TT  TT  GT  TT  GT
#> 75   GG  NN  GG  GG  GG  GT  GT  GG  GG  GG  GG  GG  GG  GG  GG  NN  GT  GT  GG
#> 76   AC  AC  AA  AC  AA  CC  AC  AA  AC  AC  AA  AC  AC  AA  AA  CC  AC  AC  AC
#> 77   TT  NN  NN  NN  NN  CC  CC  NN  NN  NN  TT  TT  TT  TT  NN  NN  CC  CC  TT
#> 78   AA  AA  AA  AT  AT  AA  AA  AT  AT  AT  AA  AT  AT  AT  AT  AA  AA  AA  NN
#> 79   GT  GG  GG  GG  GG  TT  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  GT  GG  GG
#> 80   CT  TT  CT  TT  TT  TT  TT  TT  TT  TT  CT  CT  TT  TT  CT  NN  CT  TT  TT
#> 81   CT  TT  CT  CT  CT  CT  CT  TT  CC  CT  NN  CT  CT  CT  TT  TT  CT  TT  CT
#> 82   AA  TT  AT  AT  AT  AT  NN  TT  NN  AT  AT  AA  AA  AA  AT  AA  NN  AT  AT
#> 83   TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 84   TT  TT  TT  TT  TT  TT  NN  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 85   CT  CT  TT  TT  TT  CC  CT  NN  NN  TT  TT  TT  TT  TT  TT  CT  NN  NN  TT
#> 86   GG  GG  AA  AA  AG  AG  NN  AA  AA  AA  AA  AA  AA  NN  AG  GG  NN  GG  NN
#> 87   NN  AA  NN  AG  AG  AA  AA  AG  NN  AG  NN  NN  AG  AG  AG  AA  AA  AA  NN
#> 88   CC  TT  TT  TT  TT  CT  CT  TT  TT  NN  TT  TT  TT  TT  TT  CT  NN  CT  TT
#> 89   GG  AA  GG  GG  GG  AG  AA  GG  NN  GG  GG  GG  GG  GG  GG  AA  AG  NN  GG
#> 90   CG  GG  NN  CC  CG  CG  NN  CC  NN  CC  NN  GG  CC  CG  GG  CC  CG  NN  NN
#> 91   CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 92   TT  CT  CT  CC  CC  NN  CC  CC  CC  CC  CC  CC  CC  CC  CC  TT  NN  CC  CT
#> 93   NN  NN  AG  NN  AA  AG  GG  AG  NN  GG  AG  NN  AA  AA  NN  AG  NN  NN  AG
#> 94   TT  CT  TT  TT  NN  TT  TT  NN  CT  TT  TT  TT  CT  TT  TT  TT  NN  NN  TT
#> 95   NN  CC  CG  CG  CC  CG  CC  CG  CG  CC  NN  CG  CC  NN  NN  CC  CC  CC  GG
#> 96   GG  GG  GG  GG  GG  GT  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GT  NN  GG
#> 97   AC  AC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  AC  CC  CC
#> 98   GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 99   GT  GG  GT  GT  GT  GG  GG  GG  GG  GG  GT  GT  GT  GG  NN  GG  NN  GT  GG
#> 100  AG  AA  GG  GG  NN  AG  AG  GG  GG  GG  GG  GG  GG  NN  GG  AG  AA  NN  GG
#> 101  AA  AA  NN  AT  AT  AT  AA  NN  NN  AT  NN  TT  AA  NN  AA  AT  NN  NN  AT
#>     V60 V61 V62 V63 V64 V65 V66 V67 V68 V69 V70 V71 V72 V73 V74 V75 V76 V77 V78
#> 1    NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN
#> 2    AG  GG  GG  AG  AA  AA  GG  AA  AG  NN  AA  AG  GG  AG  AG  GG  AA  AG  AG
#> 3    GG  GG  GG  GG  GG  GG  GT  GG  GG  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG
#> 4    NN  CT  CC  CC  CC  CC  CC  CT  CC  CT  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 5    AA  AA  AA  AA  AA  AG  AA  AA  AG  AA  GG  AA  AG  GG  AG  AA  AA  AG  AA
#> 6    NN  AG  NN  GG  GG  NN  NN  NN  AG  GG  GG  NN  AG  AG  GG  GG  NN  NN  GG
#> 7    AC  CC  AA  AA  AC  AC  AA  AC  NN  AA  AA  AA  AA  AA  AA  AA  AA  AC  AA
#> 8    TT  TT  TT  TT  TT  NN  TT  TT  NN  TT  TT  GT  TT  GT  NN  TT  TT  TT  GT
#> 9    GG  NN  GG  GG  GG  AG  AG  GG  NN  GG  AA  NN  GG  AG  NN  AA  AG  AG  AG
#> 10   CC  CC  CC  NN  CC  CC  AC  CC  NN  AC  NN  CC  AC  CC  CC  AC  NN  AC  CC
#> 11   AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  NN  AA  AG  AA  AA  AA  AA
#> 12   TT  TT  TT  TT  TT  NN  CT  CC  CT  CT  CT  NN  CT  CC  CT  TT  CT  CT  NN
#> 13   AA  AA  AG  AA  AG  NN  AA  AG  AA  AG  AA  NN  AG  AA  AA  NN  AA  NN  NN
#> 14   TT  TT  TT  TT  TT  TT  TT  TT  TT  CT  TT  TT  CT  CT  CT  TT  CT  TT  NN
#> 15   GG  GG  GG  GG  GG  GG  GG  GG  NN  GG  GT  GG  GG  GG  GT  GG  GT  GG  GG
#> 16   AG  GG  GG  GG  GG  AA  AA  AG  AG  AA  AG  AG  AA  AG  AG  NN  AA  AA  AG
#> 17   GG  GT  GT  GT  GG  GG  GG  GT  GG  GG  NN  GT  GG  TT  GT  NN  GT  GG  GG
#> 18   AA  AA  AC  AC  AA  AC  AA  AC  AC  AA  AA  AA  AA  AA  AA  NN  AC  AA  AA
#> 19   NN  AA  AC  AA  AC  AC  CC  CC  CC  AC  CC  CC  CC  CC  CC  CC  AC  CC  CC
#> 20   TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  GT  TT  TT  TT  TT  GT  TT  GT
#> 21   GG  GG  GG  GG  GG  GG  GG  AG  AG  GG  GG  GG  GG  GG  GG  AG  GG  GG  NN
#> 22   CC  CC  CC  CC  CC  AC  AC  CC  AC  AC  AC  NN  AC  AC  NN  CC  AC  AC  NN
#> 23   CC  CC  CC  CC  CC  CC  CC  CT  CC  CT  CC  CT  CC  CC  CC  CC  NN  CT  CC
#> 24   GG  GG  GG  GG  GG  GG  AG  GG  GG  GG  GG  GG  GG  NN  AG  GG  NN  GG  NN
#> 25   AA  AA  AA  AA  AA  NN  AA  AG  AA  AA  AA  NN  AA  AA  AA  AG  AA  AA  NN
#> 26   GG  GG  GG  GG  GG  AG  AA  AG  AA  AA  AG  GG  NN  NN  NN  AA  AG  GG  AG
#> 27   CT  CT  CT  NN  CT  CT  NN  TT  TT  TT  TT  CT  TT  CT  CT  CT  CT  CT  CT
#> 28   AA  AA  AA  AA  AA  AA  NN  AA  AA  NN  AA  AA  NN  AA  AA  AA  AC  AA  AC
#> 29   NN  NN  TT  CT  TT  TT  TT  CC  CC  CT  CT  TT  CT  CT  TT  NN  NN  CT  CT
#> 30   NN  NN  AA  AC  AA  NN  NN  AC  AA  CC  AA  NN  NN  AC  NN  NN  AA  AA  NN
#> 31   GG  GG  GG  GG  GG  GG  GG  GG  CG  GG  GG  NN  NN  GG  GG  CG  GG  GG  GG
#> 32   CC  CC  CC  CC  CC  NN  CC  GG  CG  CC  CG  CG  CC  GG  NN  CG  CC  CG  CG
#> 33   AC  NN  AC  NN  AC  AA  AC  AC  AA  AA  AA  AA  NN  AC  AA  AA  AA  AA  AC
#> 34   CT  CC  CC  CT  NN  CC  CT  CT  CC  CT  TT  CC  CT  CC  NN  CC  CC  TT  NN
#> 35   CC  CC  CC  CC  CC  CC  CG  CC  CC  CC  CG  GG  CG  CC  CC  GG  CC  CC  NN
#> 36   AA  AA  AA  AA  AA  NN  AA  AA  AA  AA  AA  NN  AA  AA  AA  AC  NN  AC  AA
#> 37   AA  AA  AA  NN  AA  AA  AA  AG  AA  NN  NN  AG  AA  AG  AA  AA  AG  AA  AA
#> 38   GG  GG  GG  GG  GG  CC  GG  CG  GG  CG  CG  NN  GG  CG  GG  NN  CG  GG  GG
#> 39   GT  GG  NN  GT  TT  GT  GT  TT  GT  NN  GT  GT  GG  GG  GT  GT  GT  GT  GG
#> 40   TT  NN  TT  AT  TT  NN  NN  AT  TT  TT  TT  TT  TT  TT  NN  NN  TT  NN  TT
#> 41   NN  GG  GG  GG  GG  GG  GG  GG  GG  CG  GG  NN  GG  GG  GG  GG  CG  GG  NN
#> 42   NN  AA  AA  AA  AA  AG  AG  AA  NN  AG  GG  NN  GG  AG  AG  GG  AG  AG  NN
#> 43   GG  GG  GG  GG  GG  AG  AA  AG  GG  AG  AG  NN  AG  AG  AG  AG  AA  AG  AA
#> 44   AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AG  AA  AA  AA
#> 45   TT  NN  AT  TT  TT  NN  TT  TT  NN  AT  AT  TT  TT  TT  TT  TT  TT  TT  TT
#> 46   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN
#> 47   CC  CC  CC  CC  NN  CC  NN  AA  CC  NN  CC  CC  NN  AC  NN  CC  AC  NN  NN
#> 48   CC  CC  CC  CC  CC  CC  CC  CG  CC  GG  CC  CC  CC  CC  CC  CC  CG  CC  CC
#> 49   NN  TT  TT  TT  CC  TT  TT  CC  CT  CT  TT  CT  NN  CT  NN  TT  CT  CT  NN
#> 50   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 51   NN  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  NN  NN  NN  CC
#> 52   NN  TT  NN  TT  CT  NN  NN  CT  CT  NN  CC  TT  NN  NN  TT  CT  NN  CT  CT
#> 53   AG  AG  NN  AG  AG  AA  AG  AA  AG  AA  AG  AA  AA  AG  AG  AA  AA  AA  AA
#> 54   AA  NN  AA  AA  AA  AA  AA  AA  AA  AC  AC  AA  NN  AC  AA  CC  AA  AA  AA
#> 55   GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN
#> 56   AG  NN  AG  AG  AA  GG  GG  AA  GG  NN  AG  GG  GG  AG  GG  GG  GG  AG  NN
#> 57   NN  AA  AG  AA  AG  AA  NN  AG  AA  AG  AA  AA  AG  AA  AA  NN  NN  AG  AA
#> 58   GG  GG  NN  NN  GG  NN  GG  GG  GG  GG  GG  NN  GG  NN  GG  NN  GG  GG  AG
#> 59   CC  NN  CC  NN  TT  NN  CC  CC  CC  CT  CT  CC  NN  CT  NN  CC  CC  CT  NN
#> 60   GG  AG  GG  AG  NN  NN  AG  NN  NN  GG  AG  AG  GG  AG  AG  AG  AA  AG  AG
#> 61   AA  AA  AA  AA  AA  AA  NN  AC  AC  AC  AA  AC  AA  AA  AC  AC  CC  AA  CC
#> 62   TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  NN  TT  NN  TT  CT  TT  TT  TT  NN
#> 63   GG  GG  GG  GG  GG  GG  GG  CG  GG  GG  GG  CG  GG  GG  GG  GG  GG  GG  GG
#> 64   NN  NN  CC  AC  CC  NN  AC  AA  AA  CC  AA  AC  AC  AC  NN  NN  AA  AA  NN
#> 65   AT  TT  NN  NN  AT  TT  AT  TT  AA  TT  AT  AT  AT  TT  AT  AA  AT  AA  TT
#> 66   CT  TT  CT  CC  CT  CC  CT  CT  CT  CC  TT  CT  TT  CC  CT  CC  CC  NN  CT
#> 67   NN  AC  AC  CC  AC  NN  AA  AC  AC  AA  AC  NN  AA  NN  AC  AA  NN  AA  CC
#> 68   AG  GG  GG  AG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 69   AT  TT  AT  NN  TT  NN  AA  AT  AA  TT  AA  AT  AA  NN  AA  NN  AT  AT  AA
#> 70   AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  NN  AG  AA  AA  AA  NN  AA  AA  AA
#> 71   AA  AA  AA  AA  AA  AC  AA  AA  NN  AA  AC  AC  AC  AA  AA  AA  AA  AA  AC
#> 72   AT  AA  AA  AT  AT  AA  AA  AA  AA  AT  AT  AA  AA  AT  AA  AA  AA  AA  AA
#> 73   AA  AG  AG  AG  AG  AG  AA  AG  AA  AA  AG  AA  AG  AA  NN  AG  AA  AG  AG
#> 74   GG  TT  NN  GT  TT  TT  TT  TT  NN  NN  GT  TT  GT  TT  TT  TT  TT  GT  GT
#> 75   GG  GG  GG  GG  GG  NN  GG  GT  GG  GT  GG  NN  NN  GT  NN  GG  GG  NN  GG
#> 76   NN  AC  AC  AC  AC  AA  CC  AA  NN  AA  CC  AC  AC  AA  AC  NN  AC  AC  CC
#> 77   NN  TT  TT  NN  TT  NN  NN  CC  CC  CC  NN  CC  NN  NN  NN  CC  CT  NN  NN
#> 78   AA  NN  AT  AA  AT  AA  NN  AT  AA  AT  AA  AA  AT  AA  NN  AA  AA  AT  NN
#> 79   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  GT  GG  GG  GG  GG  GG  GG  GG
#> 80   CT  CT  CT  TT  CT  TT  TT  CT  TT  TT  TT  TT  TT  TT  TT  TT  TT  NN  TT
#> 81   CC  CT  CT  NN  CT  TT  TT  TT  CT  NN  TT  TT  CT  TT  CT  CT  CT  CT  CT
#> 82   TT  AT  AA  AT  AT  TT  AA  AT  AA  AT  AA  AT  AT  AA  AT  AT  AA  NN  AA
#> 83   TT  TT  TT  TT  TT  NN  NN  CT  TT  TT  TT  TT  TT  TT  TT  TT  TT  NN  NN
#> 84   TT  TT  NN  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 85   TT  TT  TT  TT  TT  CC  TT  CC  CT  CC  CC  CC  CC  NN  NN  CC  CT  CT  CC
#> 86   AA  AA  AA  AG  NN  GG  GG  AG  GG  GG  GG  AG  NN  GG  AG  AG  AG  GG  AG
#> 87   AA  NN  AG  AA  AG  AA  AA  AG  AA  AG  AA  NN  AG  AA  NN  AA  AA  AG  NN
#> 88   TT  TT  TT  TT  TT  CT  TT  TT  NN  CT  CT  CT  CT  CT  TT  TT  TT  TT  CT
#> 89   GG  GG  GG  NN  NN  NN  AA  AG  AG  AA  NN  AG  NN  AG  AG  AG  AG  NN  NN
#> 90   CC  GG  NN  GG  GG  NN  NN  CC  NN  CC  CC  CG  CC  CC  CG  NN  CG  CG  CC
#> 91   CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CT  CC  CC  CC  CC  CC  CC  CC  CC
#> 92   CC  CT  CC  CC  NN  CT  TT  CC  CT  TT  TT  TT  CT  CC  CC  CT  TT  CT  NN
#> 93   NN  AG  AG  NN  AA  NN  NN  AG  GG  GG  GG  NN  NN  GG  GG  NN  NN  NN  NN
#> 94   CT  TT  TT  CT  NN  CT  TT  TT  NN  TT  CT  NN  CT  TT  TT  TT  CT  TT  TT
#> 95   CG  NN  CG  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  CC  NN  NN  CC  CC  CC
#> 96   GG  GG  GG  GG  GG  NN  GT  GT  GG  GG  GG  NN  GG  GG  GG  NN  GT  GT  NN
#> 97   CC  CC  CC  CC  CC  AC  CC  CC  CC  NN  CC  CC  CC  CC  CC  CC  CC  CC  CC
#> 98   GG  GG  GG  GG  GG  NN  GG  GG  GG  NN  GG  GG  GG  GG  GG  NN  GG  NN  GG
#> 99   GT  GT  TT  GT  TT  GT  GG  GG  GT  GG  GG  GT  GG  NN  GG  GT  GG  GG  GG
#> 100  GG  GG  GG  GG  NN  AG  GG  AG  AG  AA  AG  GG  AG  AA  GG  AG  AA  NN  AA
#> 101  AA  AT  AA  NN  AA  NN  AT  AT  AA  AT  AT  AA  AT  AA  AA  AT  AT  NN  AT
#>     V79 V80 V81 V82 V83 V84 V85 V86 V87 V88 V89 V90 V91 V92 V93 V94 V95 V96 V97
#> 1    NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN  NN
#> 2    AG  AA  AA  AG  AG  GG  AG  AG  AA  GG  AA  AA  AG  AG  AG  AA  NN  GG  NN
#> 3    NN  GG  GT  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 4    CC  CC  CC  CC  CC  CC  CC  CT  CC  CC  CT  CT  CC  CC  CC  CC  CC  CC  CC
#> 5    AA  AG  AA  AG  NN  NN  AG  AG  GG  AA  AG  GG  AG  AG  AG  AG  AG  AG  AA
#> 6    NN  GG  GG  NN  NN  AA  AA  AG  GG  AG  GG  NN  AA  AG  AG  NN  AG  NN  GG
#> 7    AA  AA  NN  NN  AA  AC  AA  AA  AA  AA  AC  AC  AC  AC  CC  AC  AA  AC  AA
#> 8    TT  GT  TT  TT  TT  TT  TT  TT  GT  GT  TT  TT  TT  TT  TT  TT  TT  NN  TT
#> 9    AG  NN  GG  NN  NN  GG  GG  GG  AG  GG  NN  GG  AG  AG  NN  GG  AG  NN  AG
#> 10   CC  AC  NN  AC  CC  AC  AC  CC  CC  CC  CC  CC  AC  CC  NN  CC  CC  NN  CC
#> 11   AA  AA  AA  AG  AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AG  NN  AA  AA  AG
#> 12   NN  TT  CT  CT  CT  CT  CC  CC  CT  TT  CT  CT  TT  TT  CT  NN  CT  CT  CT
#> 13   AA  AG  AG  AA  AG  NN  NN  AG  AA  AA  AG  AA  AG  AG  AA  NN  NN  AG  AA
#> 14   TT  CC  NN  TT  TT  CT  NN  NN  TT  TT  TT  CT  TT  TT  TT  TT  CT  TT  TT
#> 15   GT  GT  NN  GG  GG  GT  GT  GG  GT  GG  GG  GT  GT  GT  GT  GT  TT  GT  GG
#> 16   AA  AA  AA  AG  NN  AG  NN  AA  NN  AG  AG  AG  GG  NN  AG  AG  AG  AG  AG
#> 17   GT  GT  NN  NN  GG  NN  GT  GG  GT  GT  GG  GG  GT  GT  GG  GT  GG  GG  TT
#> 18   AA  AC  AA  AA  AA  AA  AC  NN  AC  AA  AC  AA  AA  AC  AA  NN  NN  AA  AC
#> 19   AC  AA  AC  AA  CC  AC  CC  CC  CC  CC  AA  CC  AA  CC  CC  CC  CC  AC  AA
#> 20   NN  TT  TT  TT  TT  GT  TT  TT  TT  TT  GT  TT  TT  TT  TT  TT  TT  TT  TT
#> 21   GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  AG  GG  GG  GG  GG  GG
#> 22   AC  CC  AC  AC  CC  AC  NN  CC  CC  AC  CC  AC  AC  NN  CC  NN  AC  NN  AC
#> 23   CC  NN  CC  CC  NN  CC  CC  CT  CC  CC  CC  TT  CT  CT  CT  CT  CT  CT  CT
#> 24   AG  GG  GG  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  AG  GG  GG  GG  GG
#> 25   NN  AA  AA  AA  AG  AG  NN  AA  AA  AA  AA  AA  AA  AA  AA  AA  NN  AG  AA
#> 26   NN  AG  AG  AG  GG  AA  GG  AG  AG  AG  GG  AG  AG  AA  GG  AG  AG  AG  AA
#> 27   CT  TT  CT  CT  CT  TT  TT  TT  NN  CT  CT  CT  CT  TT  CT  TT  NN  NN  NN
#> 28   AA  AA  AA  AA  AA  AA  AA  AC  AC  AA  AA  AA  AC  AA  AA  AC  AA  AA  NN
#> 29   CT  NN  CT  CT  CT  NN  CT  CT  CC  CC  TT  CT  CT  TT  TT  TT  TT  CT  CT
#> 30   NN  AA  AC  AA  AA  NN  AC  AA  AA  AA  AA  AA  AA  AC  AA  NN  AA  AA  AA
#> 31   CG  GG  GG  CG  CG  NN  CG  GG  GG  GG  GG  GG  GG  NN  GG  GG  NN  GG  CG
#> 32   CG  CG  NN  CG  CG  GG  NN  CC  CC  CG  GG  CG  CG  CC  CC  CC  GG  CG  CG
#> 33   AA  NN  AC  AA  AC  AC  AA  NN  NN  AC  AC  AA  AA  AC  AC  AA  AA  NN  AA
#> 34   CC  CC  TT  CC  NN  CT  TT  CC  CT  CT  CC  TT  CC  TT  CC  CT  CC  CC  CT
#> 35   CC  CG  CG  CC  CC  CC  CC  CC  CG  GG  CG  CC  CG  CC  CC  NN  NN  CG  NN
#> 36   AA  AA  AA  AA  AA  AA  AA  AA  AA  AA  AC  AA  AA  AA  AA  AA  AA  NN  AA
#> 37   AA  AG  AA  AA  AG  AA  AA  AA  AG  AG  AA  AA  AA  AG  AA  AA  AA  NN  AA
#> 38   CG  GG  CC  GG  CG  CG  GG  GG  CG  CC  GG  GG  GG  GG  CG  CG  CG  NN  CG
#> 39   GG  GT  GG  GG  TT  TT  GT  GT  TT  NN  GG  GG  GT  GG  GT  GG  GT  GG  GG
#> 40   TT  TT  AT  NN  AT  NN  NN  AT  TT  TT  TT  AT  TT  TT  TT  NN  AT  TT  NN
#> 41   GG  NN  GG  CG  CG  NN  NN  CG  CG  CG  CG  CG  GG  CG  GG  NN  GG  GG  NN
#> 42   NN  GG  GG  AG  AA  AA  GG  NN  GG  AG  AA  NN  AA  NN  GG  NN  NN  NN  GG
#> 43   AG  AA  AG  AG  AA  GG  AG  AA  AG  AG  AG  AG  GG  GG  AG  GG  AG  AG  GG
#> 44   AG  AA  AA  AA  AG  AA  AG  AA  AA  AG  NN  AA  AG  AG  NN  AG  AG  AA  AA
#> 45   TT  TT  TT  AT  NN  TT  AT  TT  TT  TT  TT  AT  TT  NN  TT  AT  TT  TT  AT
#> 46   NN  GG  NN  GG  GG  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  GG  GT
#> 47   CC  CC  CC  AC  AA  AC  NN  CC  CC  AA  CC  AC  CC  NN  CC  CC  AC  CC  AA
#> 48   NN  CG  CG  CC  CG  CC  CC  CC  CG  GG  CC  CC  CG  CC  CC  CC  CC  GG  CC
#> 49   NN  CC  CT  TT  CT  CT  CT  NN  CC  NN  CT  NN  TT  CT  CT  CT  TT  CT  NN
#> 50   GG  GG  GG  NN  CG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  NN  GG  CG
#> 51   CC  NN  CC  CT  CC  CC  CT  CC  CC  CC  CC  CC  CT  CC  CC  CC  NN  CC  CC
#> 52   TT  CT  TT  CT  CT  TT  TT  NN  NN  TT  CC  TT  CC  CT  TT  CT  NN  CT  CT
#> 53   AA  AA  AA  AG  AG  AA  AA  AA  AA  AA  GG  GG  AA  AG  AG  NN  AG  AG  AA
#> 54   AA  AA  AA  AA  AC  AC  AA  AC  AC  AA  AA  CC  AA  AC  NN  NN  AC  NN  AA
#> 55   GG  GT  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 56   NN  GG  GG  AG  AA  GG  NN  GG  GG  AG  AA  AG  AG  AG  GG  GG  NN  AA  AG
#> 57   AA  AG  AG  AA  AG  AA  AA  AG  NN  AA  AG  AA  AG  AG  AA  AG  AA  AG  AA
#> 58   GG  GG  GG  GG  NN  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  NN  GG  GG  NN
#> 59   CC  CC  CC  CT  CC  CC  CC  CT  CC  NN  NN  TT  CT  NN  NN  CT  NN  CC  CC
#> 60   GG  AG  AG  GG  AG  AG  NN  AA  AG  AA  GG  AG  GG  AG  GG  GG  AA  AG  AA
#> 61   AC  AA  AC  AC  AA  AA  AC  NN  AA  AA  AC  AA  NN  AC  AC  AC  AC  AC  AC
#> 62   TT  TT  TT  TT  TT  TT  TT  CT  TT  TT  TT  CT  TT  TT  TT  NN  TT  TT  CT
#> 63   GG  GG  GG  CG  CG  CG  GG  GG  GG  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG
#> 64   NN  AC  NN  AC  AC  AA  AA  AA  AC  NN  AC  AC  AC  AA  CC  AA  AA  AA  AA
#> 65   NN  AT  NN  TT  AT  AT  AT  AT  TT  TT  AA  TT  TT  AT  AT  TT  AT  TT  AT
#> 66   CC  CT  CT  CT  NN  CC  CC  CC  TT  CC  TT  CC  TT  CC  CT  CC  CC  CT  CT
#> 67   AC  CC  AA  AA  AA  AA  NN  AA  CC  AC  AC  AC  AA  AA  AC  NN  CC  AC  CC
#> 68   GG  GG  GG  GG  GG  NN  GG  GG  NN  GG  GG  GG  GG  GG  GG  GG  AG  GG  GG
#> 69   AA  AT  AT  AA  AA  AT  NN  NN  NN  AA  AA  TT  AT  AA  AT  AT  NN  TT  AA
#> 70   AA  AG  AA  AA  AA  AA  AA  AA  AG  NN  AA  AA  AA  AA  AA  AA  AA  AA  AA
#> 71   AA  AC  AC  AC  AC  AA  NN  AA  AC  AC  AA  AC  NN  AC  AA  AA  AC  NN  NN
#> 72   AA  AT  AA  AA  AA  AT  AA  AA  AT  AA  AA  AA  AA  AA  AA  AA  AA  AA  AT
#> 73   AA  AG  AA  GG  NN  GG  AG  AA  AA  AG  AG  GG  AG  AG  NN  AG  GG  GG  GG
#> 74   TT  TT  TT  GG  TT  TT  TT  TT  TT  NN  GT  TT  TT  GT  GT  TT  TT  TT  TT
#> 75   NN  GG  GT  NN  GG  GG  GT  GT  GG  GG  GT  GG  GG  GG  GG  GG  GT  GG  GT
#> 76   AA  AC  AC  AC  CC  CC  NN  NN  AC  AC  AC  AC  AC  AC  AC  AC  CC  CC  CC
#> 77   NN  CC  CC  CT  CC  CC  CT  NN  CT  NN  CT  CC  CC  CT  CT  CT  CT  TT  CC
#> 78   NN  AT  AT  AA  AT  AA  NN  AT  AA  AA  AT  AA  AT  AT  AA  NN  AA  AT  AA
#> 79   GG  GT  GG  GG  GG  GG  GG  GG  GT  GT  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 80   TT  TT  TT  TT  TT  NN  CC  TT  TT  TT  CT  TT  TT  TT  NN  TT  TT  CT  CT
#> 81   CT  CT  CT  CT  TT  CT  TT  TT  CT  TT  TT  CT  NN  CT  CT  CT  NN  CT  TT
#> 82   AT  AT  AA  AT  AA  AA  AA  AT  AT  AT  AT  AA  AT  NN  AT  AT  AT  TT  AA
#> 83   CT  TT  TT  NN  CT  TT  NN  TT  TT  TT  TT  TT  TT  CT  TT  TT  CT  TT  TT
#> 84   TT  TT  TT  TT  AT  TT  TT  TT  TT  AT  TT  TT  TT  TT  TT  TT  TT  TT  TT
#> 85   CC  CT  CT  CT  CT  TT  NN  CC  CT  CT  CT  TT  CC  CT  CT  CT  NN  CC  CT
#> 86   GG  AG  AG  AG  AG  AA  AG  GG  AA  NN  GG  GG  AG  GG  AG  AA  AG  AG  GG
#> 87   NN  AG  AG  AA  AG  AA  AA  AG  AA  NN  AG  AA  AG  AG  NN  AG  AA  AG  AA
#> 88   TT  CT  TT  CT  CC  TT  CT  TT  TT  CT  CT  TT  NN  TT  TT  CT  CT  CT  CT
#> 89   NN  GG  AA  AA  AG  AA  AG  AG  AG  NN  NN  AG  AG  AG  AG  AA  AG  NN  AG
#> 90   CG  CG  GG  CG  NN  CG  CG  NN  GG  CG  GG  CC  CC  CG  NN  GG  CG  GG  GG
#> 91   CC  CC  CC  CC  CT  CC  CC  CC  CC  NN  CC  CC  CC  CT  NN  CC  CC  CC  CC
#> 92   NN  CT  CT  CC  CC  CC  CT  NN  CT  CT  CT  CT  CC  TT  TT  CT  CT  NN  CC
#> 93   AG  NN  AG  NN  NN  AG  NN  GG  NN  GG  NN  GG  GG  GG  AG  NN  AG  GG  NN
#> 94   CT  TT  TT  TT  CT  TT  NN  CT  TT  CT  TT  TT  TT  CT  CT  TT  CT  NN  CT
#> 95   CC  CC  CC  CG  CC  CC  CG  CC  CC  CC  CC  CC  CC  CC  CG  CC  CG  CC  CC
#> 96   GT  GG  GG  GT  GT  GG  GG  GG  GG  GG  GG  GT  GG  GG  GG  GT  NN  GG  NN
#> 97   CC  AC  CC  CC  CC  NN  CC  CC  AC  CC  CC  CC  CC  NN  CC  CC  CC  CC  CC
#> 98   GG  GG  AG  GG  AG  GG  GG  GG  GG  AG  GG  GG  GG  GG  GG  GG  GG  GG  GG
#> 99   GT  GT  GG  GG  GG  GT  GG  GG  GT  GT  GT  GG  GG  GT  GT  NN  GT  GG  GT
#> 100  AG  GG  AG  AG  AG  GG  NN  NN  GG  AA  AG  GG  AA  NN  NN  GG  NN  AA  AA
#> 101  AT  AT  AA  TT  TT  TT  AT  AT  AA  AA  AA  AT  AA  AT  AT  AT  TT  TT  AT
#>     V98 V99 V100 V101 V102 V103
#> 1    NN  NN   NN   TT   GT   GG
#> 2    AG  AG   GG   GA   GG   NN
#> 3    GG  GG   GG   NN   NN   NN
#> 4    CC  CC   CC   NN   NN   NN
#> 5    AG  AA   AA   NN   NN   NN
#> 6    AG  GG   NN   NN   NN   NN
#> 7    AC  AA   AA   NN   NN   NN
#> 8    TT  TT   TT   NN   NN   NN
#> 9    AG  GG   GG   NN   NN   NN
#> 10   AC  AA   AC   NN   NN   NN
#> 11   AA  AG   AA   NN   NN   NN
#> 12   CT  TT   CT   NN   NN   NN
#> 13   NN  NN   AG   NN   NN   NN
#> 14   CT  NN   NN   NN   NN   NN
#> 15   GG  GT   GT   NN   NN   NN
#> 16   AG  AA   AA   NN   NN   NN
#> 17   GG  GG   NN   NN   NN   NN
#> 18   AC  AA   AA   NN   NN   NN
#> 19   AC  AC   AC   NN   NN   NN
#> 20   TT  TT   TT   NN   NN   NN
#> 21   GG  GG   AG   NN   NN   NN
#> 22   CC  CC   AA   NN   NN   NN
#> 23   CT  CC   CC   NN   NN   NN
#> 24   GG  GG   GG   NN   NN   NN
#> 25   AA  AA   AA   NN   NN   NN
#> 26   GG  GG   AG   NN   NN   NN
#> 27   CT  TT   TT   NN   NN   NN
#> 28   AA  AA   AA   NN   NN   NN
#> 29   TT  CT   CT   NN   NN   NN
#> 30   AA  AC   AC   NN   NN   NN
#> 31   GG  CG   GG   NN   NN   NN
#> 32   CC  CG   CC   NN   NN   NN
#> 33   AC  AC   NN   NN   NN   NN
#> 34   CC  TT   NN   NN   NN   NN
#> 35   CC  NN   CC   NN   NN   NN
#> 36   AA  AC   AA   NN   NN   NN
#> 37   NN  AA   AA   NN   NN   NN
#> 38   CG  GG   NN   NN   NN   NN
#> 39   GG  NN   GG   NN   NN   NN
#> 40   TT  NN   TT   NN   NN   NN
#> 41   CC  NN   CC   NN   NN   NN
#> 42   AG  AG   AG   NN   NN   NN
#> 43   NN  AG   AG   NN   NN   NN
#> 44   AG  AA   AA   NN   NN   NN
#> 45   NN  TT   AT   NN   NN   NN
#> 46   GG  GT   GG   NN   NN   NN
#> 47   CC  AC   NN   NN   NN   NN
#> 48   CC  CG   CC   NN   NN   NN
#> 49   CT  TT   CC   NN   NN   NN
#> 50   GG  GG   GG   NN   NN   NN
#> 51   CC  CC   CC   NN   NN   NN
#> 52   CT  TT   TT   NN   NN   NN
#> 53   AA  AA   AA   NN   NN   NN
#> 54   AC  CC   AC   NN   NN   NN
#> 55   GG  GG   GG   NN   NN   NN
#> 56   NN  GG   AG   NN   NN   NN
#> 57   AA  AA   AG   NN   NN   NN
#> 58   NN  GG   NN   NN   NN   NN
#> 59   CC  CC   CT   NN   NN   NN
#> 60   AG  NN   NN   NN   NN   NN
#> 61   AA  AC   AA   NN   NN   NN
#> 62   NN  NN   TT   NN   NN   NN
#> 63   GG  GG   GG   NN   NN   NN
#> 64   NN  AA   NN   NN   NN   NN
#> 65   AT  AT   AT   NN   NN   NN
#> 66   CC  TT   CT   NN   NN   NN
#> 67   CC  AC   NN   NN   NN   NN
#> 68   NN  GG   GG   NN   NN   NN
#> 69   AT  AT   AA   NN   NN   NN
#> 70   AA  AA   AA   NN   NN   NN
#> 71   AA  AC   AC   NN   NN   NN
#> 72   AA  AA   AA   NN   NN   NN
#> 73   AG  AG   AG   NN   NN   NN
#> 74   GT  GT   NN   NN   NN   NN
#> 75   GG  GG   NN   NN   NN   NN
#> 76   AC  CC   AC   NN   NN   NN
#> 77   TT  TT   NN   NN   NN   NN
#> 78   AA  AA   AT   NN   NN   NN
#> 79   GG  GG   GG   NN   NN   NN
#> 80   NN  TT   TT   NN   NN   NN
#> 81   TT  TT   TT   NN   NN   NN
#> 82   AT  AT   AT   NN   NN   NN
#> 83   CT  CT   NN   NN   NN   NN
#> 84   TT  TT   NN   NN   NN   NN
#> 85   CT  CC   NN   NN   NN   NN
#> 86   AG  AG   AG   NN   NN   NN
#> 87   NN  AA   AG   NN   NN   NN
#> 88   CT  CT   CC   NN   NN   NN
#> 89   NN  AG   AA   NN   NN   NN
#> 90   NN  CC   GG   NN   NN   NN
#> 91   CC  CC   NN   NN   NN   NN
#> 92   CC  CT   CT   NN   NN   NN
#> 93   GG  GG   AG   NN   NN   NN
#> 94   TT  CT   CT   NN   NN   NN
#> 95   NN  CC   CC   NN   NN   NN
#> 96   GG  GG   NN   NN   NN   NN
#> 97   AC  CC   AC   NN   NN   NN
#> 98   GG  GG   GG   NN   NN   NN
#> 99   NN  TT   GT   NN   NN   NN
#> 100  NN  AG   AG   NN   NN   NN
#> 101  NN  AT   AT   NN   NN   NN
```
