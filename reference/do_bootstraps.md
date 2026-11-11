# Generate bootstrapped windowed statistics

`do_bootstraps` creates a distribution of bootstrapped smoothed values
for requested statistics.

## Usage

``` r
do_bootstraps(
  x,
  facets = NULL,
  boots,
  sigma,
  step = 2 * sigma,
  statistics = "all",
  nk = TRUE,
  par = FALSE,
  triple_sigma = FALSE,
  gaussian = TRUE,
  do.p = TRUE,
  p.alt = "two-sided"
)
```

## Arguments

- x:

  snpRdata object.

- facets:

  character or NULL, default NULL. Categories by which to break up
  bootstraps.

- boots:

  numeric. Number of bootstraps to generate.

- sigma:

  numeric. Designates the width of windows in kilobases. Full window
  size is 6\*sigma.

- step:

  numeric or NULL, default default `2*sigma` (non-overlapping windows).
  Designates the number of kilobases between each window centroid. If
  NULL, windows are centered on each SNP.

- statistics:

  character. Designates the statistic(s) to smooth, typically "all",
  "single", or "pairwise". See details for options.

- nk:

  logical, default TRUE. If TRUE, weights SNP contribution to window
  averages by the number of observations at those SNPs.

- par:

  numeric or FALSE, default FALSE. If numeric, the number of cores to
  use for parallel processing.

- triple_sigma:

  Logical, default FALSE. If TRUE, sigma will be tripled to create
  windows of 6\*sigma total.

- gaussian:

  Logical, default TRUE. If TRUE, windows will be gaussian smoothed. If
  not, windows will be raw averages.

- do.p:

  logical, default TRUE. Determines if p-values should be calculated for
  sliding windows.

- p.alt:

  character, default "two-sided". Specifies the alternative hypothesis
  to be used. Options:

  - "less": probability that a bootstrapped value is as small or smaller
    than observed.

  - "greater": probability that a bootstrapped value is as large or
    larger than observed.

  - "two-sided": probability that a bootstrapped value is as or more
    extreme than observed.

## Value

A snpRdata object with bootstrapped windows merged in to the
window.bootstraps slot. If do.p is TRUE, it will also merge p values in
for bootstrapped statistics into the stats or window.stats sockets.

## Details

Bootstraps are conducted as described by Hohenlohe et al. (2010). For
each bootstrap, this function draws random window position, then draws
random statistics from all provided SNPs to fill each observed position
on that window and calculates a smoothed statistic for that window using
Gaussian-smoothing. Note that a "position" column *must* be present in
the snp metadata of the snpRdata object to do any window calculations.

Bootstraps for multiple statistics can be calculated at once. If the
statistics argument is set to "all", all calculated stats will be run.
If it is set to "single", then all non-pairwise statistics will be
bootstrapped, if it is set to "pairwise", then all pairwise statistics
will be bootstrapped. Individual statistics can also be requested by
name ("pi", "ho", etc.). All statistics bootstrapped at the same time
will be calculated using *the same randomly filled windows*, and thus do
not represent independent observations between statistics. This is still
computationally more efficient, however.

The data can be broken up categorically by snp or sample metadata, as
described in
[`Facets_in_snpR`](https://hemstrow.github.io/snpR/reference/Facets_in_snpR.md).

As described in Hohenlohe et al. (2010), the contribution of individual
SNPs to window averages can be weighted by the number of observations
per SNP by setting the nk argument to TRUE, as is the default. For
bootstraps, nk values are randomly drawn for each SNP in each window.

Possible centers for windows can either SNPs (if no step size is
provided), or every step kilobases from the 0 position of each snp level
facet category (chromosome, etc.).

If do.p is TRUE, calculates p-values for smoothed values of a statistic
based upon the bootstrapped null distribution of that statistic using an
empirical continuous distribution function. Note that in this case, the
minimum possible p-value for a window depends upon the number of
bootstraps calculated (if only a 1000 bootstraps were performed, the
minimum possible p-value is about .001, or one in a thousand.)

Note that this function will return an error if equivalent windowed
statistics have not first been calculated for the designated facets (if
the "chromosome" facet is requested with a sigma of 200 and a step of
50, do_bootstraps will error unless
[`calc_smoothed_averages`](https://hemstrow.github.io/snpR/reference/calc_smoothed_averages.md)
has not yet been run with the same facet, sigma, and step values).

## References

Hohenlohe et al. (2010). *PLOS Genetics*

## Author

William Hemstrom

## Examples

``` r
# add statistics
dat <- calc_basic_snp_stats(stickSNPs, "chr", sigma = 200, step = 150)

# do bootstraps
dat <- do_bootstraps(dat, facets = "chr", boots = 100, 
                     sigma = 200, step = 150)
#> Sigma: 200 
#> Total number of input snps:  100 
#> Number of boots: 100 
#> Unique snp level facets: chr 
#> Bootstrap facet:
#>  chr 
#> Bootstrap Progress:
#> 0%
#> 5% 
#> 10% 
#> 16% 
#> 22% 
#> 27% 
#> 33% 
#> 39% 
#> 45% 
#> 51% 
#> 56% 
#> 62% 
#> 67% 
#> 73% 
#> 78% 
#> 84% 
#> 89% 
#> 95% 
#> 100% 

# fetch results, bootstraps and then p-values on original stats
get.snpR.stats(dat, "chr", "bootstraps")
#>     facet subfacet snp.facet triple_sigma gaussian stat sigma   nk step
#> 1   .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 2   .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 3   .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 4   .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 5   .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 6   .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 7   .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 8   .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 9   .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 10  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 11  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 12  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 13  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 14  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 15  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 16  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 17  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 18  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 19  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 20  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 21  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 22  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 23  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 24  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 25  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 26  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 27  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 28  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 29  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 30  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 31  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 32  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 33  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 34  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 35  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 36  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 37  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 38  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 39  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 40  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 41  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 42  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 43  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 44  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 45  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 46  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 47  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 48  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 49  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 50  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 51  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 52  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 53  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 54  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 55  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 56  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 57  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 58  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 59  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 60  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 61  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 62  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 63  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 64  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 65  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 66  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 67  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 68  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 69  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 70  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 71  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 72  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 73  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 74  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 75  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 76  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 77  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 78  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 79  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 80  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 81  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 82  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 83  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 84  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 85  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 86  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 87  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 88  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 89  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 90  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 91  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 92  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 93  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 94  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 95  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 96  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 97  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 98  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 99  .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 100 .base    .base       chr        FALSE     TRUE   pi   200 TRUE  150
#> 101 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 102 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 103 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 104 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 105 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 106 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 107 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 108 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 109 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 110 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 111 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 112 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 113 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 114 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 115 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 116 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 117 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 118 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 119 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 120 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 121 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 122 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 123 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 124 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 125 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 126 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 127 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 128 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 129 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 130 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 131 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 132 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 133 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 134 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 135 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 136 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 137 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 138 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 139 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 140 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 141 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 142 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 143 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 144 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 145 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 146 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 147 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 148 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 149 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 150 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 151 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 152 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 153 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 154 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 155 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 156 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 157 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 158 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 159 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 160 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 161 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 162 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 163 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 164 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 165 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 166 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 167 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 168 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 169 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 170 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 171 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 172 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 173 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 174 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 175 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 176 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 177 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 178 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 179 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 180 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 181 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 182 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 183 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 184 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 185 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 186 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 187 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 188 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 189 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 190 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 191 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 192 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 193 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 194 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 195 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 196 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 197 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 198 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 199 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 200 .base    .base       chr        FALSE     TRUE   ho   200 TRUE  150
#> 201 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 202 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 203 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 204 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 205 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 206 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 207 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 208 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 209 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 210 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 211 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 212 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 213 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 214 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 215 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 216 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 217 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 218 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 219 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 220 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 221 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 222 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 223 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 224 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 225 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 226 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 227 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 228 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 229 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 230 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 231 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 232 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 233 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 234 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 235 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 236 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 237 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 238 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 239 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 240 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 241 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 242 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 243 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 244 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 245 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 246 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 247 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 248 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 249 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 250 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 251 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 252 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 253 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 254 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 255 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 256 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 257 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 258 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 259 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 260 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 261 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 262 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 263 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 264 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 265 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 266 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 267 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 268 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 269 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 270 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 271 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 272 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 273 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 274 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 275 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 276 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 277 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 278 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 279 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 280 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 281 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 282 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 283 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 284 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 285 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 286 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 287 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 288 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 289 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 290 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 291 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 292 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 293 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 294 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 295 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 296 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 297 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 298 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 299 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 300 .base    .base       chr        FALSE     TRUE pHWE   200 TRUE  150
#> 301 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 302 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 303 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 304 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 305 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 306 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 307 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 308 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 309 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 310 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 311 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 312 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 313 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 314 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 315 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 316 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 317 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 318 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 319 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 320 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 321 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 322 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 323 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 324 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 325 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 326 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 327 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 328 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 329 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 330 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 331 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 332 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 333 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 334 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 335 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 336 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 337 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 338 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 339 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 340 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 341 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 342 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 343 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 344 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 345 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 346 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 347 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 348 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 349 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 350 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 351 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 352 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 353 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 354 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 355 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 356 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 357 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 358 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 359 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 360 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 361 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 362 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 363 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 364 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 365 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 366 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 367 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 368 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 369 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 370 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 371 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 372 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 373 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 374 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 375 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 376 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 377 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 378 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 379 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 380 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 381 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 382 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 383 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 384 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 385 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 386 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 387 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 388 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 389 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 390 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 391 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 392 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 393 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 394 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 395 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 396 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 397 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 398 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 399 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 400 .base    .base       chr        FALSE     TRUE  fis   200 TRUE  150
#> 401 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 402 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 403 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 404 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 405 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 406 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 407 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 408 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 409 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 410 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 411 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 412 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 413 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 414 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 415 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 416 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 417 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 418 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 419 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 420 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 421 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 422 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 423 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 424 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 425 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 426 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 427 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 428 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 429 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 430 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 431 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 432 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 433 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 434 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 435 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 436 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 437 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 438 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 439 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 440 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 441 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 442 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 443 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 444 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 445 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 446 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 447 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 448 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 449 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 450 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 451 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 452 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 453 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 454 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 455 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 456 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 457 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 458 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 459 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 460 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 461 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 462 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 463 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 464 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 465 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 466 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 467 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 468 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 469 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 470 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 471 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 472 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 473 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 474 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 475 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 476 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 477 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 478 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 479 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 480 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 481 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 482 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 483 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 484 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 485 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 486 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 487 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 488 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 489 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 490 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 491 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 492 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 493 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 494 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 495 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 496 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 497 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 498 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 499 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#> 500 .base    .base       chr        FALSE     TRUE   he   200 TRUE  150
#>             value
#> 1    0.4388366125
#> 2    0.3822960519
#> 3    0.3396572827
#> 4    0.4994523549
#> 5    0.3656734947
#> 6    0.0955307263
#> 7    0.2581259350
#> 8    0.1472727273
#> 9    0.0225974026
#> 10   0.3628802137
#> 11   0.3937012987
#> 12   0.4684915883
#> 13   0.4277922078
#> 14   0.1269088133
#> 15   0.3396572827
#> 16   0.1792648342
#> 17   0.1733675506
#> 18   0.3644654088
#> 19   0.1142007164
#> 20   0.5028559772
#> 21   0.3425687435
#> 22   0.4139610390
#> 23   0.3712770298
#> 24   0.3619909502
#> 25   0.1089628596
#> 26   0.0504780728
#> 27   0.4388366125
#> 28   0.5033557047
#> 29   0.4967767296
#> 30   0.3806977336
#> 31   0.0504780728
#> 32   0.0955307263
#> 33   0.3644654088
#> 34   0.4684915883
#> 35   0.0306073394
#> 36   0.0306073394
#> 37   0.4994523549
#> 38   0.4714285714
#> 39   0.3652250782
#> 40   0.3527131783
#> 41   0.4277922078
#> 42   0.4388366125
#> 43   0.2220265922
#> 44   0.2795031056
#> 45   0.0793575815
#> 46   0.3204301075
#> 47   0.1792648342
#> 48   0.4613646934
#> 49   0.3019068500
#> 50   0.5012239902
#> 51   0.4357324555
#> 52   0.4994523549
#> 53   0.4700733204
#> 54   0.3019068500
#> 55   0.1130315990
#> 56   0.3656734947
#> 57   0.4662603948
#> 58   0.1720049659
#> 59   0.3656734947
#> 60   0.2600389864
#> 61   0.4277922078
#> 62   0.4277922078
#> 63   0.4981979910
#> 64   0.0648044693
#> 65   0.0648044693
#> 66   0.3245096172
#> 67   0.5012239902
#> 68   0.4388366125
#> 69   0.2768502216
#> 70   0.1967093680
#> 71   0.4789050561
#> 72   0.3656734947
#> 73   0.0627724499
#> 74   0.4407213562
#> 75   0.4522102332
#> 76   0.4967767296
#> 77   0.1628801986
#> 78   0.4981979910
#> 79   0.3870639165
#> 80   0.5033557047
#> 81   0.5033557047
#> 82   0.4850439410
#> 83   0.0627724499
#> 84   0.5033557047
#> 85   0.1142007164
#> 86   0.3204301075
#> 87   0.1142007164
#> 88   0.4063450674
#> 89   0.0627724499
#> 90   0.3937012987
#> 91   0.4267535692
#> 92   0.3527131783
#> 93   0.0306073394
#> 94   0.1967093680
#> 95   0.0836303160
#> 96   0.3952339436
#> 97   0.4981979910
#> 98   0.2564245810
#> 99   0.0648044693
#> 100  0.2220265922
#> 101  0.4047619048
#> 102  0.3617021277
#> 103  0.4302325581
#> 104  0.3855421687
#> 105  0.3888888889
#> 106  0.0777777778
#> 107  0.2325581395
#> 108  0.1590909091
#> 109  0.0227272727
#> 110  0.3626373626
#> 111  0.4431818182
#> 112  0.4285714286
#> 113  0.4318181818
#> 114  0.1354166667
#> 115  0.4302325581
#> 116  0.1979166667
#> 117  0.1666666667
#> 118  0.3750000000
#> 119  0.0989010989
#> 120  0.4767441860
#> 121  0.3411764706
#> 122  0.4204545455
#> 123  0.4418604651
#> 124  0.4705882353
#> 125  0.1149425287
#> 126  0.0515463918
#> 127  0.4047619048
#> 128  0.4400000000
#> 129  0.5125000000
#> 130  0.5064935065
#> 131  0.0515463918
#> 132  0.0777777778
#> 133  0.3750000000
#> 134  0.4285714286
#> 135  0.0309278351
#> 136  0.0309278351
#> 137  0.3855421687
#> 138  0.4772727273
#> 139  0.3372093023
#> 140  0.4534883721
#> 141  0.4318181818
#> 142  0.4047619048
#> 143  0.2527472527
#> 144  0.3333333333
#> 145  0.0821917808
#> 146  0.3548387097
#> 147  0.1979166667
#> 148  0.3908045977
#> 149  0.2758620690
#> 150  0.5000000000
#> 151  0.4878048780
#> 152  0.3855421687
#> 153  0.4268292683
#> 154  0.2758620690
#> 155  0.1195652174
#> 156  0.3888888889
#> 157  0.4831460674
#> 158  0.1666666667
#> 159  0.3888888889
#> 160  0.2421052632
#> 161  0.4545454545
#> 162  0.4545454545
#> 163  0.4814814815
#> 164  0.0666666667
#> 165  0.0666666667
#> 166  0.3595505618
#> 167  0.5000000000
#> 168  0.4047619048
#> 169  0.2857142857
#> 170  0.2197802198
#> 171  0.6206896552
#> 172  0.3888888889
#> 173  0.0645161290
#> 174  0.3723404255
#> 175  0.4000000000
#> 176  0.5125000000
#> 177  0.1777777778
#> 178  0.4814814815
#> 179  0.2597402597
#> 180  0.4400000000
#> 181  0.4400000000
#> 182  0.5063291139
#> 183  0.0645161290
#> 184  0.4400000000
#> 185  0.0989010989
#> 186  0.3548387097
#> 187  0.0989010989
#> 188  0.3898305085
#> 189  0.0645161290
#> 190  0.4431818182
#> 191  0.5222222222
#> 192  0.4534883721
#> 193  0.0309278351
#> 194  0.2197802198
#> 195  0.0869565217
#> 196  0.3870967742
#> 197  0.4814814815
#> 198  0.3000000000
#> 199  0.0666666667
#> 200  0.2527472527
#> 201  0.6157488726
#> 202  0.5952381578
#> 203  0.0093652672
#> 204  0.0468486225
#> 205  0.7710853254
#> 206  0.1894128961
#> 207  0.3946036410
#> 208  1.0000000000
#> 209  1.0000000000
#> 210  1.0000000000
#> 211  0.2832921808
#> 212  0.4849435852
#> 213  1.0000000000
#> 214  1.0000000000
#> 215  0.0093652672
#> 216  0.5916578006
#> 217  0.5453537899
#> 218  1.0000000000
#> 219  0.2745021501
#> 220  0.6696626548
#> 221  1.0000000000
#> 222  1.0000000000
#> 223  0.0856641117
#> 224  0.0045545792
#> 225  1.0000000000
#> 226  1.0000000000
#> 227  0.6157488726
#> 228  0.3549029963
#> 229  0.8230104487
#> 230  0.0019224541
#> 231  1.0000000000
#> 232  0.1894128961
#> 233  1.0000000000
#> 234  0.4849435852
#> 235  1.0000000000
#> 236  1.0000000000
#> 237  0.0468486225
#> 238  1.0000000000
#> 239  0.5518748234
#> 240  0.0047465338
#> 241  1.0000000000
#> 242  0.6157488726
#> 243  0.3467137289
#> 244  0.1108728613
#> 245  1.0000000000
#> 246  0.5113073257
#> 247  0.5916578006
#> 248  0.1663864888
#> 249  0.4727813865
#> 250  1.0000000000
#> 251  0.3157466665
#> 252  0.0468486225
#> 253  0.4788401118
#> 254  0.4727813865
#> 255  1.0000000000
#> 256  0.7710853254
#> 257  0.8202605039
#> 258  0.5653493973
#> 259  0.7710853254
#> 260  0.4455942033
#> 261  0.6217596253
#> 262  0.6217596253
#> 263  0.8238416373
#> 264  1.0000000000
#> 265  1.0000000000
#> 266  0.5078101487
#> 267  1.0000000000
#> 268  0.6157488726
#> 269  1.0000000000
#> 270  0.5913226546
#> 271  0.0067076617
#> 272  0.7710853254
#> 273  1.0000000000
#> 274  0.1586180018
#> 275  0.3358125606
#> 276  0.8230104487
#> 277  1.0000000000
#> 278  0.8238416373
#> 279  0.0063540720
#> 280  0.3549029963
#> 281  0.3549029963
#> 282  0.8156303048
#> 283  1.0000000000
#> 284  0.3549029963
#> 285  0.2745021501
#> 286  0.5113073257
#> 287  0.2745021501
#> 288  0.7539182309
#> 289  1.0000000000
#> 290  0.2832921808
#> 291  0.0450618612
#> 292  0.0047465338
#> 293  1.0000000000
#> 294  0.5913226546
#> 295  1.0000000000
#> 296  1.0000000000
#> 297  0.8238416373
#> 298  0.2036723562
#> 299  1.0000000000
#> 300  0.3467137289
#> 301  0.0780790591
#> 302  0.0541429853
#> 303 -0.2686567164
#> 304  0.2291421857
#> 305 -0.0638661202
#> 306  0.1866840731
#> 307  0.0995762712
#> 308 -0.0807453416
#> 309 -0.0057803468
#> 310  0.0006729475
#> 311 -0.1264940239
#> 312  0.0856793146
#> 313 -0.0094656489
#> 314 -0.0674157303
#> 315 -0.2686567164
#> 316 -0.1046511628
#> 317  0.0388751034
#> 318 -0.0290924881
#> 319  0.1346153846
#> 320  0.0522164808
#> 321  0.0040883074
#> 322 -0.0157778479
#> 323 -0.1914422722
#> 324 -0.3023255814
#> 325 -0.0552147239
#> 326 -0.0212765957
#> 327  0.0780790591
#> 328  0.1266094421
#> 329 -0.0318572794
#> 330 -0.3333333333
#> 331 -0.0212765957
#> 332  0.1866840731
#> 333 -0.0290924881
#> 334  0.0856793146
#> 335 -0.0105263158
#> 336 -0.0105263158
#> 337  0.2291421857
#> 338 -0.0124688279
#> 339  0.0771246724
#> 340 -0.2878787879
#> 341 -0.0094656489
#> 342  0.0780790591
#> 343 -0.1392405063
#> 344 -0.1940298507
#> 345 -0.0359712230
#> 346 -0.1080291971
#> 347 -0.1046511628
#> 348  0.1536903039
#> 349  0.0867256637
#> 350  0.0024563319
#> 351 -0.1203319502
#> 352  0.2291421857
#> 353  0.0925096031
#> 354  0.0867256637
#> 355 -0.0581395349
#> 356 -0.0638661202
#> 357 -0.0364283758
#> 358  0.0312046444
#> 359 -0.0638661202
#> 360  0.0693069307
#> 361 -0.0629199756
#> 362 -0.0629199756
#> 363  0.0337565810
#> 364 -0.0289017341
#> 365 -0.0289017341
#> 366 -0.1086614173
#> 367  0.0024563319
#> 368  0.0780790591
#> 369 -0.0322011469
#> 370 -0.1180124224
#> 371 -0.2982946603
#> 372 -0.0638661202
#> 373 -0.0279329609
#> 374  0.1558609959
#> 375  0.1160631383
#> 376 -0.0318572794
#> 377 -0.0920245399
#> 378  0.0337565810
#> 379  0.3303964758
#> 380  0.1266094421
#> 381  0.1266094421
#> 382 -0.0441767068
#> 383 -0.0279329609
#> 384  0.1266094421
#> 385  0.1346153846
#> 386 -0.1080291971
#> 387  0.1346153846
#> 388  0.0409777139
#> 389 -0.0279329609
#> 390 -0.1264940239
#> 391 -0.2252489748
#> 392 -0.2878787879
#> 393 -0.0105263158
#> 394 -0.1180124224
#> 395 -0.0400000000
#> 396  0.0206978119
#> 397  0.0337565810
#> 398 -0.1710526316
#> 399 -0.0289017341
#> 400 -0.1392405063
#> 401  0.4362244888
#> 402  0.3802625633
#> 403  0.3376825322
#> 404  0.4964436065
#> 405  0.3636419765
#> 406  0.0950000000
#> 407  0.2566252018
#> 408  0.1464359428
#> 409  0.0224690154
#> 410  0.3608863703
#> 411  0.3914643553
#> 412  0.4657029483
#> 413  0.4253615688
#> 414  0.1262478241
#> 415  0.3376825322
#> 416  0.1783311578
#> 417  0.1723356086
#> 418  0.3621875000
#> 419  0.1135732392
#> 420  0.4999323960
#> 421  0.3405536345
#> 422  0.4116089899
#> 423  0.3691184460
#> 424  0.3598615942
#> 425  0.1083366285
#> 426  0.0502178843
#> 427  0.4362244888
#> 428  0.5000000000
#> 429  0.4936718750
#> 430  0.3782256672
#> 431  0.0502178843
#> 432  0.0950000000
#> 433  0.3621875000
#> 434  0.4657029483
#> 435  0.0304495744
#> 436  0.0304495744
#> 437  0.4964436065
#> 438  0.4687500000
#> 439  0.3631016734
#> 440  0.3506625246
#> 441  0.4253615688
#> 442  0.4362244888
#> 443  0.2208066713
#> 444  0.2777777822
#> 445  0.0788140357
#> 446  0.3187073636
#> 447  0.1783311578
#> 448  0.4587131727
#> 449  0.3001717583
#> 450  0.4983098969
#> 451  0.4330755497
#> 452  0.4964436065
#> 453  0.4672070199
#> 454  0.3001717583
#> 455  0.1124172991
#> 456  0.3636419765
#> 457  0.4636409548
#> 458  0.1710493755
#> 459  0.3636419765
#> 460  0.2586703616
#> 461  0.4253615688
#> 462  0.4253615688
#> 463  0.4951226939
#> 464  0.0644444382
#> 465  0.0644444382
#> 466  0.3226865283
#> 467  0.4983098969
#> 468  0.4362244888
#> 469  0.2753290601
#> 470  0.1956285474
#> 471  0.4761527292
#> 472  0.3636419765
#> 473  0.0624349551
#> 474  0.4383770971
#> 475  0.4495501726
#> 476  0.4936718750
#> 477  0.1619753105
#> 478  0.4951226939
#> 479  0.3845505147
#> 480  0.5000000000
#> 481  0.5000000000
#> 482  0.4819740422
#> 483  0.0624349551
#> 484  0.5000000000
#> 485  0.1135732392
#> 486  0.3187073636
#> 487  0.1135732392
#> 488  0.4029014678
#> 489  0.0624349551
#> 490  0.3914643553
#> 491  0.4243827195
#> 492  0.3506625246
#> 493  0.0304495744
#> 494  0.1956285474
#> 495  0.0831758018
#> 496  0.3931090260
#> 497  0.4951226939
#> 498  0.2550000000
#> 499  0.0644444382
#> 500  0.2208066713
get.snpR.stats(dat, "chr", "single.window")
#>     facet subfacet snp.facet snp.subfacet position sigma n_snps step nk.status
#> 1   .base    .base       chr       groupI 10500000   200      1  150      TRUE
#> 2   .base    .base       chr       groupI 10650000   200      1  150      TRUE
#> 3   .base    .base       chr       groupI 16350000   200      1  150      TRUE
#> 4   .base    .base       chr       groupI 16500000   200      1  150      TRUE
#> 5   .base    .base       chr       groupI 16650000   200      1  150      TRUE
#> 6   .base    .base       chr       groupI 19650000   200      1  150      TRUE
#> 7   .base    .base       chr       groupI 19800000   200      1  150      TRUE
#> 8   .base    .base       chr       groupI 19950000   200      1  150      TRUE
#> 9   .base    .base       chr       groupI 21150000   200      1  150      TRUE
#> 10  .base    .base       chr       groupI 21300000   200      1  150      TRUE
#> 11  .base    .base       chr       groupI 21450000   200      1  150      TRUE
#> 12  .base    .base       chr       groupI 22650000   200      1  150      TRUE
#> 13  .base    .base       chr       groupI 22800000   200      1  150      TRUE
#> 14  .base    .base       chr      groupII 10650000   200      1  150      TRUE
#> 15  .base    .base       chr      groupII 10800000   200      1  150      TRUE
#> 16  .base    .base       chr      groupII 10950000   200      1  150      TRUE
#> 17  .base    .base       chr      groupII 11400000   200      1  150      TRUE
#> 18  .base    .base       chr      groupII 11550000   200      1  150      TRUE
#> 19  .base    .base       chr      groupII 11700000   200      1  150      TRUE
#> 20  .base    .base       chr      groupII 16500000   200      1  150      TRUE
#> 21  .base    .base       chr      groupII 16650000   200      1  150      TRUE
#> 22  .base    .base       chr      groupII 19050000   200      1  150      TRUE
#> 23  .base    .base       chr      groupII 19200000   200      1  150      TRUE
#> 24  .base    .base       chr      groupII 19350000   200      1  150      TRUE
#> 25  .base    .base       chr      groupII 20400000   200      1  150      TRUE
#> 26  .base    .base       chr      groupII 20550000   200      1  150      TRUE
#> 27  .base    .base       chr     groupIII  2700000   200      1  150      TRUE
#> 28  .base    .base       chr     groupIII  2850000   200      1  150      TRUE
#> 29  .base    .base       chr     groupIII  4650000   200      1  150      TRUE
#> 30  .base    .base       chr     groupIII  4800000   200      1  150      TRUE
#> 31  .base    .base       chr     groupIII  7050000   200      1  150      TRUE
#> 32  .base    .base       chr     groupIII  7200000   200      1  150      TRUE
#> 33  .base    .base       chr     groupIII  7350000   200      1  150      TRUE
#> 34  .base    .base       chr     groupIII 13050000   200      1  150      TRUE
#> 35  .base    .base       chr     groupIII 13200000   200      1  150      TRUE
#> 36  .base    .base       chr     groupIII 13350000   200      1  150      TRUE
#> 37  .base    .base       chr     groupIII 15900000   200      1  150      TRUE
#> 38  .base    .base       chr     groupIII 16050000   200      1  150      TRUE
#> 39  .base    .base       chr      groupIV  6150000   200      1  150      TRUE
#> 40  .base    .base       chr      groupIV  6300000   200      1  150      TRUE
#> 41  .base    .base       chr      groupIV  6450000   200      1  150      TRUE
#> 42  .base    .base       chr      groupIV 17100000   200      1  150      TRUE
#> 43  .base    .base       chr      groupIV 17250000   200      1  150      TRUE
#> 44  .base    .base       chr      groupIV 17400000   200      1  150      TRUE
#> 45  .base    .base       chr      groupIV 22050000   200      1  150      TRUE
#> 46  .base    .base       chr      groupIV 22200000   200      1  150      TRUE
#> 47  .base    .base       chr      groupIV 22350000   200      1  150      TRUE
#> 48  .base    .base       chr      groupIV 32250000   200      1  150      TRUE
#> 49  .base    .base       chr      groupIV 32400000   200      1  150      TRUE
#> 50  .base    .base       chr      groupIX  1950000   200      1  150      TRUE
#> 51  .base    .base       chr      groupIX  2100000   200      1  150      TRUE
#> 52  .base    .base       chr      groupIX  2250000   200      1  150      TRUE
#> 53  .base    .base       chr      groupIX  5100000   200      1  150      TRUE
#> 54  .base    .base       chr      groupIX  5250000   200      1  150      TRUE
#> 55  .base    .base       chr      groupIX  5400000   200      1  150      TRUE
#> 56  .base    .base       chr      groupIX  5850000   200      1  150      TRUE
#> 57  .base    .base       chr      groupIX  6000000   200      1  150      TRUE
#> 58  .base    .base       chr      groupIX  6150000   200      1  150      TRUE
#> 59  .base    .base       chr      groupIX  6450000   200      1  150      TRUE
#> 60  .base    .base       chr      groupIX  6600000   200      1  150      TRUE
#> 61  .base    .base       chr      groupIX  8250000   200      1  150      TRUE
#> 62  .base    .base       chr      groupIX  8400000   200      1  150      TRUE
#> 63  .base    .base       chr      groupIX  8550000   200      1  150      TRUE
#> 64  .base    .base       chr      groupIX  9150000   200      1  150      TRUE
#> 65  .base    .base       chr      groupIX  9300000   200      1  150      TRUE
#> 66  .base    .base       chr      groupIX  9450000   200      1  150      TRUE
#> 67  .base    .base       chr      groupIX 11850000   200      1  150      TRUE
#> 68  .base    .base       chr      groupIX 12000000   200      1  150      TRUE
#> 69  .base    .base       chr      groupIX 14100000   200      1  150      TRUE
#> 70  .base    .base       chr      groupIX 14250000   200      1  150      TRUE
#> 71  .base    .base       chr      groupIX 14700000   200      1  150      TRUE
#> 72  .base    .base       chr      groupIX 14850000   200      1  150      TRUE
#> 73  .base    .base       chr      groupIX 15000000   200      2  150      TRUE
#> 74  .base    .base       chr      groupIX 15150000   200      1  150      TRUE
#> 75  .base    .base       chr      groupIX 15300000   200      1  150      TRUE
#> 76  .base    .base       chr      groupIX 16200000   200      1  150      TRUE
#> 77  .base    .base       chr      groupIX 16350000   200      1  150      TRUE
#> 78  .base    .base       chr      groupIX 16500000   200      1  150      TRUE
#> 79  .base    .base       chr      groupIX 17850000   200      1  150      TRUE
#> 80  .base    .base       chr      groupIX 18000000   200      1  150      TRUE
#> 81  .base    .base       chr      groupIX 18150000   200      1  150      TRUE
#> 82  .base    .base       chr      groupIX 18300000   200      1  150      TRUE
#> 83  .base    .base       chr      groupIX 18450000   200      1  150      TRUE
#> 84  .base    .base       chr       groupV  1200000   200      1  150      TRUE
#> 85  .base    .base       chr       groupV  1350000   200      1  150      TRUE
#> 86  .base    .base       chr       groupV  1500000   200      1  150      TRUE
#> 87  .base    .base       chr       groupV  9600000   200      1  150      TRUE
#> 88  .base    .base       chr       groupV  9750000   200      1  150      TRUE
#> 89  .base    .base       chr       groupV  9900000   200      1  150      TRUE
#> 90  .base    .base       chr      groupVI   150000   200      1  150      TRUE
#> 91  .base    .base       chr      groupVI   300000   200      1  150      TRUE
#> 92  .base    .base       chr      groupVI   900000   200      1  150      TRUE
#> 93  .base    .base       chr      groupVI  1050000   200      1  150      TRUE
#> 94  .base    .base       chr      groupVI  1200000   200      1  150      TRUE
#> 95  .base    .base       chr      groupVI 13200000   200      1  150      TRUE
#> 96  .base    .base       chr      groupVI 13350000   200      1  150      TRUE
#> 97  .base    .base       chr      groupVI 13500000   200      1  150      TRUE
#> 98  .base    .base       chr      groupVI 13800000   200      1  150      TRUE
#> 99  .base    .base       chr      groupVI 13950000   200      1  150      TRUE
#> 100 .base    .base       chr      groupVI 14100000   200      1  150      TRUE
#> 101 .base    .base       chr      groupVI 14250000   200      1  150      TRUE
#> 102 .base    .base       chr     groupVII  3600000   200      1  150      TRUE
#> 103 .base    .base       chr     groupVII  3750000   200      1  150      TRUE
#> 104 .base    .base       chr     groupVII  4050000   200      1  150      TRUE
#> 105 .base    .base       chr     groupVII  4200000   200      1  150      TRUE
#> 106 .base    .base       chr     groupVII  4350000   200      1  150      TRUE
#> 107 .base    .base       chr     groupVII  6900000   200      1  150      TRUE
#> 108 .base    .base       chr     groupVII  7050000   200      1  150      TRUE
#> 109 .base    .base       chr     groupVII 11400000   200      1  150      TRUE
#> 110 .base    .base       chr     groupVII 11550000   200      1  150      TRUE
#> 111 .base    .base       chr     groupVII 15300000   200      1  150      TRUE
#> 112 .base    .base       chr     groupVII 15450000   200      1  150      TRUE
#> 113 .base    .base       chr     groupVII 15600000   200      1  150      TRUE
#> 114 .base    .base       chr     groupVII 20100000   200      1  150      TRUE
#> 115 .base    .base       chr     groupVII 20250000   200      1  150      TRUE
#> 116 .base    .base       chr     groupVII 25050000   200      1  150      TRUE
#> 117 .base    .base       chr     groupVII 25200000   200      1  150      TRUE
#> 118 .base    .base       chr     groupVII 25350000   200      1  150      TRUE
#> 119 .base    .base       chr     groupVII 27450000   200      1  150      TRUE
#> 120 .base    .base       chr     groupVII 27600000   200      1  150      TRUE
#> 121 .base    .base       chr    groupVIII  5400000   200      1  150      TRUE
#> 122 .base    .base       chr    groupVIII  5550000   200      1  150      TRUE
#> 123 .base    .base       chr    groupVIII  7050000   200      1  150      TRUE
#> 124 .base    .base       chr    groupVIII  7200000   200      1  150      TRUE
#> 125 .base    .base       chr    groupVIII  7350000   200      1  150      TRUE
#> 126 .base    .base       chr    groupVIII 11250000   200      1  150      TRUE
#> 127 .base    .base       chr    groupVIII 11400000   200      1  150      TRUE
#> 128 .base    .base       chr    groupVIII 11550000   200      1  150      TRUE
#> 129 .base    .base       chr    groupVIII 12000000   200      1  150      TRUE
#> 130 .base    .base       chr    groupVIII 12150000   200      1  150      TRUE
#> 131 .base    .base       chr    groupVIII 12300000   200      1  150      TRUE
#> 132 .base    .base       chr    groupVIII 14250000   200      1  150      TRUE
#> 133 .base    .base       chr    groupVIII 14400000   200      1  150      TRUE
#> 134 .base    .base       chr    groupVIII 14550000   200      1  150      TRUE
#> 135 .base    .base       chr       groupX  5100000   200      1  150      TRUE
#> 136 .base    .base       chr       groupX  5250000   200      1  150      TRUE
#> 137 .base    .base       chr       groupX 13650000   200      1  150      TRUE
#> 138 .base    .base       chr       groupX 13800000   200      1  150      TRUE
#> 139 .base    .base       chr       groupX 14100000   200      1  150      TRUE
#> 140 .base    .base       chr       groupX 14250000   200      1  150      TRUE
#> 141 .base    .base       chr       groupX 14400000   200      1  150      TRUE
#> 142 .base    .base       chr      groupXI  6000000   200      1  150      TRUE
#> 143 .base    .base       chr      groupXI  6150000   200      1  150      TRUE
#> 144 .base    .base       chr      groupXI  6300000   200      1  150      TRUE
#> 145 .base    .base       chr      groupXI  8850000   200      1  150      TRUE
#> 146 .base    .base       chr      groupXI  9000000   200      1  150      TRUE
#> 147 .base    .base       chr      groupXI  9150000   200      1  150      TRUE
#> 148 .base    .base       chr      groupXI 12750000   200      1  150      TRUE
#> 149 .base    .base       chr      groupXI 12900000   200      1  150      TRUE
#> 150 .base    .base       chr     groupXII  1950000   200      1  150      TRUE
#> 151 .base    .base       chr     groupXII  2100000   200      1  150      TRUE
#> 152 .base    .base       chr     groupXII  2250000   200      1  150      TRUE
#> 153 .base    .base       chr     groupXII  5100000   200      1  150      TRUE
#> 154 .base    .base       chr     groupXII  5250000   200      1  150      TRUE
#> 155 .base    .base       chr     groupXII  6600000   200      1  150      TRUE
#> 156 .base    .base       chr     groupXII  6750000   200      1  150      TRUE
#> 157 .base    .base       chr     groupXII  6900000   200      1  150      TRUE
#> 158 .base    .base       chr     groupXII  9750000   200      1  150      TRUE
#> 159 .base    .base       chr     groupXII  9900000   200      1  150      TRUE
#> 160 .base    .base       chr     groupXII 10050000   200      1  150      TRUE
#> 161 .base    .base       chr     groupXII 14100000   200      1  150      TRUE
#> 162 .base    .base       chr     groupXII 14250000   200      1  150      TRUE
#> 163 .base    .base       chr     groupXII 15150000   200      1  150      TRUE
#> 164 .base    .base       chr     groupXII 15300000   200      1  150      TRUE
#> 165 .base    .base       chr     groupXII 15450000   200      1  150      TRUE
#> 166 .base    .base       chr    groupXIII   750000   200      1  150      TRUE
#> 167 .base    .base       chr    groupXIII   900000   200      1  150      TRUE
#> 168 .base    .base       chr    groupXIII  1050000   200      1  150      TRUE
#> 169 .base    .base       chr    groupXIII  1950000   200      1  150      TRUE
#> 170 .base    .base       chr    groupXIII  2100000   200      1  150      TRUE
#> 171 .base    .base       chr    groupXIII  2250000   200      1  150      TRUE
#> 172 .base    .base       chr    groupXIII  4650000   200      1  150      TRUE
#> 173 .base    .base       chr    groupXIII  4800000   200      1  150      TRUE
#> 174 .base    .base       chr    groupXIII  6000000   200      1  150      TRUE
#> 175 .base    .base       chr    groupXIII  6150000   200      1  150      TRUE
#> 176 .base    .base       chr    groupXIII 15600000   200      1  150      TRUE
#> 177 .base    .base       chr    groupXIII 15750000   200      1  150      TRUE
#> 178 .base    .base       chr    groupXIII 16350000   200      1  150      TRUE
#> 179 .base    .base       chr    groupXIII 16500000   200      1  150      TRUE
#> 180 .base    .base       chr    groupXIII 16650000   200      1  150      TRUE
#> 181 .base    .base       chr    groupXIII 17700000   200      1  150      TRUE
#> 182 .base    .base       chr    groupXIII 17850000   200      1  150      TRUE
#> 183 .base    .base       chr    groupXIII 18000000   200      1  150      TRUE
#> 184 .base    .base       chr    groupXIII 18750000   200      1  150      TRUE
#> 185 .base    .base       chr    groupXIII 18900000   200      1  150      TRUE
#> 186 .base    .base       chr     groupXIV 13950000   200      1  150      TRUE
#> 187 .base    .base       chr     groupXIV 14100000   200      1  150      TRUE
#> 188 .base    .base       chr     groupXIX  2850000   200      1  150      TRUE
#> 189 .base    .base       chr     groupXIX  3000000   200      1  150      TRUE
#> 190 .base    .base       chr     groupXIX  3150000   200      1  150      TRUE
#> 191 .base    .base       chr     groupXIX 10650000   200      1  150      TRUE
#> 192 .base    .base       chr     groupXIX 10800000   200      1  150      TRUE
#> 193 .base    .base       chr     groupXIX 10950000   200      1  150      TRUE
#> 194 .base    .base       chr     groupXIX 12300000   200      1  150      TRUE
#> 195 .base    .base       chr     groupXIX 12450000   200      1  150      TRUE
#> 196 .base    .base       chr     groupXIX 12600000   200      1  150      TRUE
#> 197 .base    .base       chr     groupXIX 16050000   200      1  150      TRUE
#> 198 .base    .base       chr     groupXIX 16200000   200      1  150      TRUE
#> 199 .base    .base       chr     groupXIX 16350000   200      1  150      TRUE
#> 200 .base    .base       chr     groupXIX 17700000   200      1  150      TRUE
#> 201 .base    .base       chr     groupXIX 17850000   200      1  150      TRUE
#> 202 .base    .base       chr     groupXIX 18000000   200      1  150      TRUE
#> 203 .base    .base       chr      groupXV   750000   200      1  150      TRUE
#> 204 .base    .base       chr      groupXV   900000   200      1  150      TRUE
#> 205 .base    .base       chr      groupXV  1050000   200      1  150      TRUE
#> 206 .base    .base       chr      groupXV  8700000   200      1  150      TRUE
#> 207 .base    .base       chr      groupXV  8850000   200      1  150      TRUE
#> 208 .base    .base       chr      groupXV  9000000   200      1  150      TRUE
#> 209 .base    .base       chr      groupXV 15450000   200      1  150      TRUE
#> 210 .base    .base       chr      groupXV 15600000   200      1  150      TRUE
#> 211 .base    .base       chr     groupXVI  3450000   200      1  150      TRUE
#> 212 .base    .base       chr     groupXVI  3600000   200      1  150      TRUE
#> 213 .base    .base       chr     groupXVI  3750000   200      1  150      TRUE
#> 214 .base    .base       chr     groupXVI  7050000   200      1  150      TRUE
#> 215 .base    .base       chr     groupXVI  7200000   200      1  150      TRUE
#> 216 .base    .base       chr     groupXVI  7350000   200      1  150      TRUE
#> 217 .base    .base       chr     groupXVI 16500000   200      1  150      TRUE
#> 218 .base    .base       chr     groupXVI 16650000   200      1  150      TRUE
#> 219 .base    .base       chr    groupXVII  2400000   200      1  150      TRUE
#> 220 .base    .base       chr    groupXVII  2550000   200      1  150      TRUE
#> 221 .base    .base       chr    groupXVII  2700000   200      1  150      TRUE
#> 222 .base    .base       chr    groupXVII  4950000   200      1  150      TRUE
#> 223 .base    .base       chr    groupXVII  5100000   200      1  150      TRUE
#> 224 .base    .base       chr    groupXVII  9000000   200      1  150      TRUE
#> 225 .base    .base       chr    groupXVII  9150000   200      1  150      TRUE
#> 226 .base    .base       chr    groupXVII  9300000   200      1  150      TRUE
#> 227 .base    .base       chr    groupXVII 11100000   200      1  150      TRUE
#> 228 .base    .base       chr    groupXVII 11250000   200      1  150      TRUE
#> 229 .base    .base       chr    groupXVII 11400000   200      1  150      TRUE
#> 230 .base    .base       chr   groupXVIII   300000   200      1  150      TRUE
#> 231 .base    .base       chr   groupXVIII   450000   200      1  150      TRUE
#> 232 .base    .base       chr   groupXVIII  7050000   200      1  150      TRUE
#> 233 .base    .base       chr   groupXVIII  7200000   200      1  150      TRUE
#> 234 .base    .base       chr   groupXVIII  7350000   200      1  150      TRUE
#> 235 .base    .base       chr   groupXVIII  7950000   200      1  150      TRUE
#> 236 .base    .base       chr   groupXVIII  8100000   200      1  150      TRUE
#> 237 .base    .base       chr   groupXVIII  8250000   200      1  150      TRUE
#> 238 .base    .base       chr   groupXVIII 10800000   200      1  150      TRUE
#> 239 .base    .base       chr   groupXVIII 10950000   200      1  150      TRUE
#> 240 .base    .base       chr   groupXVIII 11100000   200      1  150      TRUE
#> 241 .base    .base       chr   groupXVIII 12750000   200      1  150      TRUE
#> 242 .base    .base       chr   groupXVIII 12900000   200      1  150      TRUE
#> 243 .base    .base       chr   groupXVIII 13050000   200      1  150      TRUE
#> 244 .base    .base       chr      groupXX  4350000   200      1  150      TRUE
#> 245 .base    .base       chr      groupXX  4500000   200      1  150      TRUE
#> 246 .base    .base       chr      groupXX  4650000   200      1  150      TRUE
#> 247 .base    .base       chr      groupXX 10200000   200      1  150      TRUE
#> 248 .base    .base       chr      groupXX 10350000   200      1  150      TRUE
#> 249 .base    .base       chr      groupXX 10800000   200      1  150      TRUE
#> 250 .base    .base       chr      groupXX 10950000   200      1  150      TRUE
#> 251 .base    .base       chr      groupXX 11100000   200      1  150      TRUE
#> 252 .base    .base       chr      groupXX 15300000   200      1  150      TRUE
#> 253 .base    .base       chr      groupXX 15450000   200      1  150      TRUE
#> 254 .base    .base       chr      groupXX 15600000   200      1  150      TRUE
#> 255 .base    .base       chr      groupXX 17550000   200      1  150      TRUE
#> 256 .base    .base       chr      groupXX 17700000   200      1  150      TRUE
#> 257 .base    .base       chr      groupXX 17850000   200      1  150      TRUE
#> 258 .base    .base       chr      groupXX 19200000   200      1  150      TRUE
#> 259 .base    .base       chr      groupXX 19350000   200      1  150      TRUE
#> 260 .base    .base       chr     groupXXI  9000000   200      1  150      TRUE
#> 261 .base    .base       chr     groupXXI  9150000   200      1  150      TRUE
#>     gaussian triple_sigma        maf         pi         pHWE pHWE_overall_BY
#> 1       TRUE        FALSE 0.01136364 0.02259740 1.0000000000       1.0000000
#> 2       TRUE        FALSE 0.01136364 0.02259740 1.0000000000       1.0000000
#> 3       TRUE        FALSE 0.13297872 0.23182387 1.0000000000       1.0000000
#> 4       TRUE        FALSE 0.13297872 0.23182387 1.0000000000       1.0000000
#> 5       TRUE        FALSE 0.13297872 0.23182387 1.0000000000       1.0000000
#> 6       TRUE        FALSE 0.35344828 0.46101949 0.5745139787       1.0000000
#> 7       TRUE        FALSE 0.35344828 0.46101949 0.5745139787       1.0000000
#> 8       TRUE        FALSE 0.35344828 0.46101949 0.5745139787       1.0000000
#> 9       TRUE        FALSE 0.15116279 0.25812593 0.3946036410       1.0000000
#> 10      TRUE        FALSE 0.15116279 0.25812593 0.3946036410       1.0000000
#> 11      TRUE        FALSE 0.15116279 0.25812593 0.3946036410       1.0000000
#> 12      TRUE        FALSE 0.03333333 0.06480447 1.0000000000       1.0000000
#> 13      TRUE        FALSE 0.03333333 0.06480447 1.0000000000       1.0000000
#> 14      TRUE        FALSE 0.24731183 0.37430979 0.4176680333       1.0000000
#> 15      TRUE        FALSE 0.24731183 0.37430979 0.4176680333       1.0000000
#> 16      TRUE        FALSE 0.24731183 0.37430979 0.4176680333       1.0000000
#> 17      TRUE        FALSE 0.22093023 0.34625323 0.5451131188       1.0000000
#> 18      TRUE        FALSE 0.22093023 0.34625323 0.5451131188       1.0000000
#> 19      TRUE        FALSE 0.22093023 0.34625323 0.5451131188       1.0000000
#> 20      TRUE        FALSE 0.07954545 0.14727273 1.0000000000       1.0000000
#> 21      TRUE        FALSE 0.07954545 0.14727273 1.0000000000       1.0000000
#> 22      TRUE        FALSE 0.03125000 0.06086387 1.0000000000       1.0000000
#> 23      TRUE        FALSE 0.03125000 0.06086387 1.0000000000       1.0000000
#> 24      TRUE        FALSE 0.03125000 0.06086387 1.0000000000       1.0000000
#> 25      TRUE        FALSE 0.09523810 0.17336755 0.5453537899       1.0000000
#> 26      TRUE        FALSE 0.09523810 0.17336755 0.5453537899       1.0000000
#> 27      TRUE        FALSE 0.40697674 0.48551612 1.0000000000       1.0000000
#> 28      TRUE        FALSE 0.40697674 0.48551612 1.0000000000       1.0000000
#> 29      TRUE        FALSE 0.32446809 0.44072136 0.1586180018       1.0000000
#> 30      TRUE        FALSE 0.32446809 0.44072136 0.1586180018       1.0000000
#> 31      TRUE        FALSE 0.26881720 0.39523394 1.0000000000       1.0000000
#> 32      TRUE        FALSE 0.26881720 0.39523394 1.0000000000       1.0000000
#> 33      TRUE        FALSE 0.26881720 0.39523394 1.0000000000       1.0000000
#> 34      TRUE        FALSE 0.35555556 0.46083178 0.6473190987       1.0000000
#> 35      TRUE        FALSE 0.35555556 0.46083178 0.6473190987       1.0000000
#> 36      TRUE        FALSE 0.35555556 0.46083178 0.6473190987       1.0000000
#> 37      TRUE        FALSE 0.40677966 0.48674489 0.0009766236       0.4986248
#> 38      TRUE        FALSE 0.40677966 0.48674489 0.0009766236       0.4986248
#> 39      TRUE        FALSE 0.10989011 0.19670937 0.5913226546       1.0000000
#> 40      TRUE        FALSE 0.10989011 0.19670937 0.5913226546       1.0000000
#> 41      TRUE        FALSE 0.10989011 0.19670937 0.5913226546       1.0000000
#> 42      TRUE        FALSE 0.05670103 0.10752631 1.0000000000       1.0000000
#> 43      TRUE        FALSE 0.05670103 0.10752631 1.0000000000       1.0000000
#> 44      TRUE        FALSE 0.05670103 0.10752631 1.0000000000       1.0000000
#> 45      TRUE        FALSE 0.08888889 0.16288020 1.0000000000       1.0000000
#> 46      TRUE        FALSE 0.08888889 0.16288020 1.0000000000       1.0000000
#> 47      TRUE        FALSE 0.08888889 0.16288020 1.0000000000       1.0000000
#> 48      TRUE        FALSE 0.40506329 0.48504394 0.8156303048       1.0000000
#> 49      TRUE        FALSE 0.40506329 0.48504394 0.8156303048       1.0000000
#> 50      TRUE        FALSE 0.09444444 0.17200497 0.5653493973       1.0000000
#> 51      TRUE        FALSE 0.09444444 0.17200497 0.5653493973       1.0000000
#> 52      TRUE        FALSE 0.09444444 0.17200497 0.5653493973       1.0000000
#> 53      TRUE        FALSE 0.23837209 0.36522508 0.5518748234       1.0000000
#> 54      TRUE        FALSE 0.23837209 0.36522508 0.5518748234       1.0000000
#> 55      TRUE        FALSE 0.23837209 0.36522508 0.5518748234       1.0000000
#> 56      TRUE        FALSE 0.06043956 0.11420072 0.2745021501       1.0000000
#> 57      TRUE        FALSE 0.06043956 0.11420072 0.2745021501       1.0000000
#> 58      TRUE        FALSE 0.06043956 0.11420072 0.2745021501       1.0000000
#> 59      TRUE        FALSE 0.37195122 0.47007332 0.4788401118       1.0000000
#> 60      TRUE        FALSE 0.37195122 0.47007332 0.4788401118       1.0000000
#> 61      TRUE        FALSE 0.13186813 0.23022282 1.0000000000       1.0000000
#> 62      TRUE        FALSE 0.13186813 0.23022282 1.0000000000       1.0000000
#> 63      TRUE        FALSE 0.13186813 0.23022282 1.0000000000       1.0000000
#> 64      TRUE        FALSE 0.50000000 0.50335570 0.3549029963       1.0000000
#> 65      TRUE        FALSE 0.50000000 0.50335570 0.3549029963       1.0000000
#> 66      TRUE        FALSE 0.50000000 0.50335570 0.3549029963       1.0000000
#> 67      TRUE        FALSE 0.01546392 0.03060734 1.0000000000       1.0000000
#> 68      TRUE        FALSE 0.01546392 0.03060734 1.0000000000       1.0000000
#> 69      TRUE        FALSE 0.38202247 0.47483019 0.1816157055       1.0000000
#> 70      TRUE        FALSE 0.38202247 0.47483019 0.1816157055       1.0000000
#> 71      TRUE        FALSE 0.10752688 0.19296716 0.0637270727       1.0000000
#> 72      TRUE        FALSE 0.10752688 0.19296716 0.0637270727       1.0000000
#> 73      TRUE        FALSE 0.16510951 0.26809970 0.0729704505       1.0000000
#> 74      TRUE        FALSE 0.24418605 0.37127703 0.0856641117       1.0000000
#> 75      TRUE        FALSE 0.24418605 0.37127703 0.0856641117       1.0000000
#> 76      TRUE        FALSE 0.11413043 0.20331433 0.0861747630       1.0000000
#> 77      TRUE        FALSE 0.11413043 0.20331433 0.0861747630       1.0000000
#> 78      TRUE        FALSE 0.11413043 0.20331433 0.0861747630       1.0000000
#> 79      TRUE        FALSE 0.25531915 0.38229605 0.5952381578       1.0000000
#> 80      TRUE        FALSE 0.25531915 0.38229605 0.5952381578       1.0000000
#> 81      TRUE        FALSE 0.48026316 0.50252701 0.8209148938       1.0000000
#> 82      TRUE        FALSE 0.48026316 0.50252701 0.8209148938       1.0000000
#> 83      TRUE        FALSE 0.48026316 0.50252701 0.8209148938       1.0000000
#> 84      TRUE        FALSE 0.19892473 0.32043011 0.5113073257       1.0000000
#> 85      TRUE        FALSE 0.19892473 0.32043011 0.5113073257       1.0000000
#> 86      TRUE        FALSE 0.19892473 0.32043011 0.5113073257       1.0000000
#> 87      TRUE        FALSE 0.37500000 0.47142857 1.0000000000       1.0000000
#> 88      TRUE        FALSE 0.37500000 0.47142857 1.0000000000       1.0000000
#> 89      TRUE        FALSE 0.37500000 0.47142857 1.0000000000       1.0000000
#> 90      TRUE        FALSE 0.49418605 0.50285598 0.6696626548       1.0000000
#> 91      TRUE        FALSE 0.49418605 0.50285598 0.6696626548       1.0000000
#> 92      TRUE        FALSE 0.27966102 0.40634507 0.7539182309       1.0000000
#> 93      TRUE        FALSE 0.27966102 0.40634507 0.7539182309       1.0000000
#> 94      TRUE        FALSE 0.27966102 0.40634507 0.7539182309       1.0000000
#> 95      TRUE        FALSE 0.30681818 0.42779221 0.6217596253       1.0000000
#> 96      TRUE        FALSE 0.30681818 0.42779221 0.6217596253       1.0000000
#> 97      TRUE        FALSE 0.30681818 0.42779221 0.6217596253       1.0000000
#> 98      TRUE        FALSE 0.09895833 0.17926483 0.5916578006       1.0000000
#> 99      TRUE        FALSE 0.09895833 0.17926483 0.5916578006       1.0000000
#> 100     TRUE        FALSE 0.45061728 0.49819799 0.8238416373       1.0000000
#> 101     TRUE        FALSE 0.45061728 0.49819799 0.8238416373       1.0000000
#> 102     TRUE        FALSE 0.15263158 0.26003899 0.4455942033       1.0000000
#> 103     TRUE        FALSE 0.15263158 0.26003899 0.4455942033       1.0000000
#> 104     TRUE        FALSE 0.35632184 0.46136469 0.1663864888       1.0000000
#> 105     TRUE        FALSE 0.35632184 0.46136469 0.1663864888       1.0000000
#> 106     TRUE        FALSE 0.35632184 0.46136469 0.1663864888       1.0000000
#> 107     TRUE        FALSE 0.09090909 0.16623377 1.0000000000       1.0000000
#> 108     TRUE        FALSE 0.09090909 0.16623377 1.0000000000       1.0000000
#> 109     TRUE        FALSE 0.30769231 0.42838929 0.6216568634       1.0000000
#> 110     TRUE        FALSE 0.30769231 0.42838929 0.6216568634       1.0000000
#> 111     TRUE        FALSE 0.20224719 0.32450962 0.5078101487       1.0000000
#> 112     TRUE        FALSE 0.20224719 0.32450962 0.5078101487       1.0000000
#> 113     TRUE        FALSE 0.20224719 0.32450962 0.5078101487       1.0000000
#> 114     TRUE        FALSE 0.15000000 0.25642458 0.2036723562       1.0000000
#> 115     TRUE        FALSE 0.15000000 0.25642458 0.2036723562       1.0000000
#> 116     TRUE        FALSE 0.23626374 0.36288021 1.0000000000       1.0000000
#> 117     TRUE        FALSE 0.23626374 0.36288021 1.0000000000       1.0000000
#> 118     TRUE        FALSE 0.23626374 0.36288021 1.0000000000       1.0000000
#> 119     TRUE        FALSE 0.43902439 0.49558581 0.1781610827       1.0000000
#> 120     TRUE        FALSE 0.43902439 0.49558581 0.1781610827       1.0000000
#> 121     TRUE        FALSE 0.04347826 0.08363032 1.0000000000       1.0000000
#> 122     TRUE        FALSE 0.04347826 0.08363032 1.0000000000       1.0000000
#> 123     TRUE        FALSE 0.32142857 0.43883661 0.6157488726       1.0000000
#> 124     TRUE        FALSE 0.32142857 0.43883661 0.6157488726       1.0000000
#> 125     TRUE        FALSE 0.32142857 0.43883661 0.6157488726       1.0000000
#> 126     TRUE        FALSE 0.30821918 0.42938120 0.7848884601       1.0000000
#> 127     TRUE        FALSE 0.30821918 0.42938120 0.7848884601       1.0000000
#> 128     TRUE        FALSE 0.30821918 0.42938120 0.7848884601       1.0000000
#> 129     TRUE        FALSE 0.34939759 0.45739321 1.0000000000       1.0000000
#> 130     TRUE        FALSE 0.34939759 0.45739321 1.0000000000       1.0000000
#> 131     TRUE        FALSE 0.34939759 0.45739321 1.0000000000       1.0000000
#> 132     TRUE        FALSE 0.34117647 0.45221023 0.3358125606       1.0000000
#> 133     TRUE        FALSE 0.34117647 0.45221023 0.3358125606       1.0000000
#> 134     TRUE        FALSE 0.34117647 0.45221023 0.3358125606       1.0000000
#> 135     TRUE        FALSE 0.04395604 0.08451217 1.0000000000       1.0000000
#> 136     TRUE        FALSE 0.04395604 0.08451217 1.0000000000       1.0000000
#> 137     TRUE        FALSE 0.04109589 0.07935758 1.0000000000       1.0000000
#> 138     TRUE        FALSE 0.04109589 0.07935758 1.0000000000       1.0000000
#> 139     TRUE        FALSE 0.05000000 0.09553073 0.1894128961       1.0000000
#> 140     TRUE        FALSE 0.05000000 0.09553073 0.1894128961       1.0000000
#> 141     TRUE        FALSE 0.05000000 0.09553073 0.1894128961       1.0000000
#> 142     TRUE        FALSE 0.28235294 0.40765750 0.1854042068       1.0000000
#> 143     TRUE        FALSE 0.28235294 0.40765750 0.1854042068       1.0000000
#> 144     TRUE        FALSE 0.28235294 0.40765750 0.1854042068       1.0000000
#> 145     TRUE        FALSE 0.30681818 0.42779221 1.0000000000       1.0000000
#> 146     TRUE        FALSE 0.30681818 0.42779221 1.0000000000       1.0000000
#> 147     TRUE        FALSE 0.30681818 0.42779221 1.0000000000       1.0000000
#> 148     TRUE        FALSE 0.04109589 0.07935758 1.0000000000       1.0000000
#> 149     TRUE        FALSE 0.04109589 0.07935758 1.0000000000       1.0000000
#> 150     TRUE        FALSE 0.16483516 0.27685022 1.0000000000       1.0000000
#> 151     TRUE        FALSE 0.16483516 0.27685022 1.0000000000       1.0000000
#> 152     TRUE        FALSE 0.16483516 0.27685022 1.0000000000       1.0000000
#> 153     TRUE        FALSE 0.28750000 0.41226415 0.2716893594       1.0000000
#> 154     TRUE        FALSE 0.28750000 0.41226415 0.2716893594       1.0000000
#> 155     TRUE        FALSE 0.18181818 0.29946524 0.7037786543       1.0000000
#> 156     TRUE        FALSE 0.18181818 0.29946524 0.7037786543       1.0000000
#> 157     TRUE        FALSE 0.18181818 0.29946524 0.7037786543       1.0000000
#> 158     TRUE        FALSE 0.12637363 0.22202659 0.3467137289       1.0000000
#> 159     TRUE        FALSE 0.12637363 0.22202659 0.3467137289       1.0000000
#> 160     TRUE        FALSE 0.12637363 0.22202659 0.3467137289       1.0000000
#> 161     TRUE        FALSE 0.36516854 0.46626039 0.8202605039       1.0000000
#> 162     TRUE        FALSE 0.36516854 0.46626039 0.8202605039       1.0000000
#> 163     TRUE        FALSE 0.08080808 0.14931036 1.0000000000       1.0000000
#> 164     TRUE        FALSE 0.08080808 0.14931036 1.0000000000       1.0000000
#> 165     TRUE        FALSE 0.08080808 0.14931036 1.0000000000       1.0000000
#> 166     TRUE        FALSE 0.26704545 0.39370130 0.2832921808       1.0000000
#> 167     TRUE        FALSE 0.26704545 0.39370130 0.2832921808       1.0000000
#> 168     TRUE        FALSE 0.26704545 0.39370130 0.2832921808       1.0000000
#> 169     TRUE        FALSE 0.21764706 0.34256874 1.0000000000       1.0000000
#> 170     TRUE        FALSE 0.21764706 0.34256874 1.0000000000       1.0000000
#> 171     TRUE        FALSE 0.21764706 0.34256874 1.0000000000       1.0000000
#> 172     TRUE        FALSE 0.14835165 0.25408293 0.6824678761       1.0000000
#> 173     TRUE        FALSE 0.14835165 0.25408293 0.6824678761       1.0000000
#> 174     TRUE        FALSE 0.28977273 0.41396104 1.0000000000       1.0000000
#> 175     TRUE        FALSE 0.28977273 0.41396104 1.0000000000       1.0000000
#> 176     TRUE        FALSE 0.48295455 0.50227273 0.0865209741       1.0000000
#> 177     TRUE        FALSE 0.48295455 0.50227273 0.0865209741       1.0000000
#> 178     TRUE        FALSE 0.30555556 0.42675357 0.0450618612       1.0000000
#> 179     TRUE        FALSE 0.30555556 0.42675357 0.0450618612       1.0000000
#> 180     TRUE        FALSE 0.30555556 0.42675357 0.0450618612       1.0000000
#> 181     TRUE        FALSE 0.47093023 0.50122399 1.0000000000       1.0000000
#> 182     TRUE        FALSE 0.47093023 0.50122399 1.0000000000       1.0000000
#> 183     TRUE        FALSE 0.47093023 0.50122399 1.0000000000       1.0000000
#> 184     TRUE        FALSE 0.41666667 0.48951049 0.2308441699       1.0000000
#> 185     TRUE        FALSE 0.41666667 0.48951049 0.2308441699       1.0000000
#> 186     TRUE        FALSE 0.36428571 0.46649538 0.7962799890       1.0000000
#> 187     TRUE        FALSE 0.36428571 0.46649538 0.7962799890       1.0000000
#> 188     TRUE        FALSE 0.22727273 0.35353535 0.0084317296       0.6072647
#> 189     TRUE        FALSE 0.22727273 0.35353535 0.0084317296       0.6072647
#> 190     TRUE        FALSE 0.22727273 0.35353535 0.0084317296       0.6072647
#> 191     TRUE        FALSE 0.25974026 0.38706392 0.0063540720       0.5799196
#> 192     TRUE        FALSE 0.25974026 0.38706392 0.0063540720       0.5799196
#> 193     TRUE        FALSE 0.25974026 0.38706392 0.0063540720       0.5799196
#> 194     TRUE        FALSE 0.23529412 0.36199095 0.0045545792       0.5799196
#> 195     TRUE        FALSE 0.23529412 0.36199095 0.0045545792       0.5799196
#> 196     TRUE        FALSE 0.23529412 0.36199095 0.0045545792       0.5799196
#> 197     TRUE        FALSE 0.22674419 0.35271318 0.0047465338       0.5799196
#> 198     TRUE        FALSE 0.22674419 0.35271318 0.0047465338       0.5799196
#> 199     TRUE        FALSE 0.22674419 0.35271318 0.0047465338       0.5799196
#> 200     TRUE        FALSE 0.25324675 0.38069773 0.0019224541       0.4986248
#> 201     TRUE        FALSE 0.25324675 0.38069773 0.0019224541       0.4986248
#> 202     TRUE        FALSE 0.25324675 0.38069773 0.0019224541       0.4986248
#> 203     TRUE        FALSE 0.06770833 0.12690881 1.0000000000       1.0000000
#> 204     TRUE        FALSE 0.06770833 0.12690881 1.0000000000       1.0000000
#> 205     TRUE        FALSE 0.06770833 0.12690881 1.0000000000       1.0000000
#> 206     TRUE        FALSE 0.23888889 0.36567349 0.7710853254       1.0000000
#> 207     TRUE        FALSE 0.23888889 0.36567349 0.7710853254       1.0000000
#> 208     TRUE        FALSE 0.23888889 0.36567349 0.7710853254       1.0000000
#> 209     TRUE        FALSE 0.16666667 0.27950311 0.1108728613       1.0000000
#> 210     TRUE        FALSE 0.16666667 0.27950311 0.1108728613       1.0000000
#> 211     TRUE        FALSE 0.18390805 0.30190685 0.4727813865       1.0000000
#> 212     TRUE        FALSE 0.18390805 0.30190685 0.4727813865       1.0000000
#> 213     TRUE        FALSE 0.18390805 0.30190685 0.4727813865       1.0000000
#> 214     TRUE        FALSE 0.21511628 0.33965728 0.0093652672       0.6072647
#> 215     TRUE        FALSE 0.21511628 0.33965728 0.0093652672       0.6072647
#> 216     TRUE        FALSE 0.21511628 0.33965728 0.0093652672       0.6072647
#> 217     TRUE        FALSE 0.39080460 0.47890506 0.0067076617       0.5799196
#> 218     TRUE        FALSE 0.39080460 0.47890506 0.0067076617       0.5799196
#> 219     TRUE        FALSE 0.05263158 0.10025063 1.0000000000       1.0000000
#> 220     TRUE        FALSE 0.05263158 0.10025063 1.0000000000       1.0000000
#> 221     TRUE        FALSE 0.05263158 0.10025063 1.0000000000       1.0000000
#> 222     TRUE        FALSE 0.05978261 0.11303160 1.0000000000       1.0000000
#> 223     TRUE        FALSE 0.05978261 0.11303160 1.0000000000       1.0000000
#> 224     TRUE        FALSE 0.23750000 0.36446541 1.0000000000       1.0000000
#> 225     TRUE        FALSE 0.23750000 0.36446541 1.0000000000       1.0000000
#> 226     TRUE        FALSE 0.23750000 0.36446541 1.0000000000       1.0000000
#> 227     TRUE        FALSE 0.08045977 0.14882732 1.0000000000       1.0000000
#> 228     TRUE        FALSE 0.08045977 0.14882732 1.0000000000       1.0000000
#> 229     TRUE        FALSE 0.08045977 0.14882732 1.0000000000       1.0000000
#> 230     TRUE        FALSE 0.03225806 0.06277245 1.0000000000       1.0000000
#> 231     TRUE        FALSE 0.03225806 0.06277245 1.0000000000       1.0000000
#> 232     TRUE        FALSE 0.25903614 0.38619934 0.3988157088       1.0000000
#> 233     TRUE        FALSE 0.25903614 0.38619934 0.3988157088       1.0000000
#> 234     TRUE        FALSE 0.25903614 0.38619934 0.3988157088       1.0000000
#> 235     TRUE        FALSE 0.05747126 0.10896286 1.0000000000       1.0000000
#> 236     TRUE        FALSE 0.05747126 0.10896286 1.0000000000       1.0000000
#> 237     TRUE        FALSE 0.05747126 0.10896286 1.0000000000       1.0000000
#> 238     TRUE        FALSE 0.02577320 0.05047807 1.0000000000       1.0000000
#> 239     TRUE        FALSE 0.02577320 0.05047807 1.0000000000       1.0000000
#> 240     TRUE        FALSE 0.02577320 0.05047807 1.0000000000       1.0000000
#> 241     TRUE        FALSE 0.26219512 0.38927129 0.7806955345       1.0000000
#> 242     TRUE        FALSE 0.26219512 0.38927129 0.7806955345       1.0000000
#> 243     TRUE        FALSE 0.26219512 0.38927129 0.7806955345       1.0000000
#> 244     TRUE        FALSE 0.31707317 0.43573246 0.3157466665       1.0000000
#> 245     TRUE        FALSE 0.31707317 0.43573246 0.3157466665       1.0000000
#> 246     TRUE        FALSE 0.31707317 0.43573246 0.3157466665       1.0000000
#> 247     TRUE        FALSE 0.15730337 0.26661588 0.2155797046       1.0000000
#> 248     TRUE        FALSE 0.15730337 0.26661588 0.2155797046       1.0000000
#> 249     TRUE        FALSE 0.44375000 0.49677673 0.8230104487       1.0000000
#> 250     TRUE        FALSE 0.44375000 0.49677673 0.8230104487       1.0000000
#> 251     TRUE        FALSE 0.44375000 0.49677673 0.8230104487       1.0000000
#> 252     TRUE        FALSE 0.46111111 0.49975171 0.8317187491       1.0000000
#> 253     TRUE        FALSE 0.46111111 0.49975171 0.8317187491       1.0000000
#> 254     TRUE        FALSE 0.46111111 0.49975171 0.8317187491       1.0000000
#> 255     TRUE        FALSE 0.45783133 0.49945235 0.0468486225       1.0000000
#> 256     TRUE        FALSE 0.45783133 0.49945235 0.0468486225       1.0000000
#> 257     TRUE        FALSE 0.45783133 0.49945235 0.0468486225       1.0000000
#> 258     TRUE        FALSE 0.36904762 0.46849159 0.4849435852       1.0000000
#> 259     TRUE        FALSE 0.36904762 0.46849159 0.4849435852       1.0000000
#> 260     TRUE        FALSE 0.13076923 0.22909958 0.5830459866       1.0000000
#> 261     TRUE        FALSE 0.13076923 0.22909958 0.5830459866       1.0000000
#>             ho         he           fis p_fis p_he p_ho p_pHWE p_pi
#> 1   0.02272727 0.02246902 -0.0057803468  0.84 0.02 0.02   0.00 0.02
#> 2   0.02272727 0.02246902 -0.0057803468  0.84 0.02 0.02   0.00 0.02
#> 3   0.24468085 0.23059076 -0.0557749260  0.66 0.64 0.64   0.00 0.64
#> 4   0.24468085 0.23059076 -0.0557749260  0.66 0.64 0.64   0.00 0.64
#> 5   0.24468085 0.23059076 -0.0557749260  0.66 0.64 0.64   0.00 0.64
#> 6   0.50000000 0.45704519 -0.0853578464  0.48 0.46 0.12   0.90 0.46
#> 7   0.50000000 0.45704519 -0.0853578464  0.48 0.46 0.12   0.90 0.46
#> 8   0.50000000 0.45704519 -0.0853578464  0.48 0.46 0.12   0.90 0.46
#> 9   0.23255814 0.25662520  0.0995762712  0.32 0.68 0.60   0.66 0.66
#> 10  0.23255814 0.25662520  0.0995762712  0.32 0.68 0.60   0.66 0.68
#> 11  0.23255814 0.25662520  0.0995762712  0.32 0.66 0.60   0.64 0.68
#> 12  0.06666667 0.06444444 -0.0289017341  0.88 0.24 0.24   0.00 0.24
#> 13  0.06666667 0.06444444 -0.0289017341  0.92 0.18 0.24   0.00 0.24
#> 14  0.34408602 0.37229738  0.0811485643  0.44 0.86 0.84   0.66 0.86
#> 15  0.34408602 0.37229738  0.0811485643  0.44 0.86 0.84   0.66 0.86
#> 16  0.34408602 0.37229738  0.0811485643  0.44 0.86 0.84   0.66 0.86
#> 17  0.32558140 0.34424013  0.0600315956  0.56 0.90 0.78   0.84 0.90
#> 18  0.32558140 0.34424013  0.0600315956  0.56 0.90 0.78   0.84 0.90
#> 19  0.32558140 0.34424013  0.0600315956  0.56 0.90 0.78   0.84 0.90
#> 20  0.15909091 0.14643594 -0.0807453416  0.50 0.46 0.46   0.00 0.44
#> 21  0.15909091 0.14643594 -0.0807453416  0.50 0.46 0.46   0.00 0.46
#> 22  0.06250000 0.06054688 -0.0270270270  0.98 0.12 0.12   0.00 0.12
#> 23  0.06250000 0.06054688 -0.0270270270  0.98 0.12 0.12   0.00 0.12
#> 24  0.06250000 0.06054687 -0.0270270270  0.98 0.12 0.12   0.00 0.12
#> 25  0.16666667 0.17233561  0.0388751034  0.62 0.52 0.50   0.84 0.52
#> 26  0.16666667 0.17233561  0.0388751034  0.62 0.52 0.50   0.86 0.52
#> 27  0.48837209 0.48269335 -0.0059171598  0.84 0.30 0.16   0.58 0.30
#> 28  0.48837209 0.48269335 -0.0059171598  0.84 0.30 0.16   0.58 0.30
#> 29  0.37234043 0.43837710  0.1558609959  0.12 0.48 0.96   0.30 0.48
#> 30  0.37234043 0.43837710  0.1558609959  0.12 0.48 0.96   0.30 0.48
#> 31  0.38709677 0.39310903  0.0206978119  0.72 0.74 0.92   0.00 0.74
#> 32  0.38709677 0.39310903  0.0206978119  0.72 0.74 0.92   0.00 0.76
#> 33  0.38709677 0.39310903  0.0206978119  0.72 0.74 0.92   0.00 0.74
#> 34  0.48888889 0.45827161 -0.0612466125  0.64 0.46 0.16   0.88 0.46
#> 35  0.48888889 0.45827161 -0.0612466125  0.64 0.46 0.16   0.88 0.46
#> 36  0.48888889 0.45827161 -0.0612466125  0.64 0.46 0.16   0.88 0.46
#> 37  0.27118644 0.48261994  0.4449760766  0.00 0.30 0.70   0.00 0.30
#> 38  0.27118644 0.48261994  0.4449760766  0.00 0.30 0.70   0.00 0.30
#> 39  0.21978022 0.19562855 -0.1180124224  0.36 0.60 0.60   0.94 0.60
#> 40  0.21978022 0.19562855 -0.1180124224  0.36 0.60 0.60   0.94 0.60
#> 41  0.21978022 0.19562855 -0.1180124224  0.32 0.60 0.60   0.94 0.60
#> 42  0.11340206 0.10697205 -0.0549450549  0.68 0.32 0.38   0.00 0.32
#> 43  0.11340206 0.10697205 -0.0549450549  0.68 0.32 0.38   0.00 0.32
#> 44  0.11340206 0.10697205 -0.0549450549  0.68 0.32 0.38   0.00 0.32
#> 45  0.17777778 0.16197531 -0.0920245399  0.48 0.48 0.52   0.00 0.48
#> 46  0.17777778 0.16197531 -0.0920245399  0.46 0.48 0.52   0.00 0.48
#> 47  0.17777778 0.16197531 -0.0920245399  0.48 0.48 0.52   0.00 0.48
#> 48  0.50632911 0.48197404 -0.0441767068  0.68 0.30 0.10   0.74 0.30
#> 49  0.50632911 0.48197404 -0.0441767068  0.68 0.30 0.10   0.74 0.30
#> 50  0.16666667 0.17104938  0.0312046444  0.72 0.48 0.46   0.90 0.50
#> 51  0.16666667 0.17104938  0.0312046444  0.72 0.48 0.50   0.90 0.50
#> 52  0.16666667 0.17104938  0.0312046444  0.70 0.50 0.50   0.90 0.50
#> 53  0.33720930 0.36310167  0.0771246724  0.54 0.96 0.80   0.86 0.96
#> 54  0.33720930 0.36310167  0.0771246724  0.52 0.98 0.82   0.88 0.96
#> 55  0.33720930 0.36310167  0.0771246724  0.52 0.96 0.82   0.88 0.96
#> 56  0.09890110 0.11357324  0.1346153846  0.16 0.42 0.38   0.42 0.42
#> 57  0.09890110 0.11357324  0.1346153846  0.16 0.42 0.38   0.44 0.40
#> 58  0.09890110 0.11357324  0.1346153846  0.16 0.42 0.38   0.42 0.42
#> 59  0.42682927 0.46720702  0.0925096031  0.34 0.36 0.66   0.74 0.36
#> 60  0.42682927 0.46720702  0.0925096031  0.34 0.36 0.66   0.72 0.38
#> 61  0.24175824 0.22895785 -0.0503978780  0.68 0.64 0.62   0.00 0.64
#> 62  0.24175824 0.22895785 -0.0503978780  0.68 0.64 0.62   0.00 0.64
#> 63  0.24175824 0.22895785 -0.0503978780  0.68 0.64 0.62   0.00 0.64
#> 64  0.44000000 0.50000000  0.1266094421  0.22 0.00 0.46   0.62 0.00
#> 65  0.44000000 0.50000000  0.1266094421  0.22 0.00 0.46   0.64 0.00
#> 66  0.44000000 0.50000000  0.1266094421  0.22 0.00 0.46   0.64 0.00
#> 67  0.03092784 0.03044957 -0.0105263158  0.88 0.06 0.04   0.00 0.08
#> 68  0.03092784 0.03044957 -0.0105263158  0.88 0.08 0.08   0.00 0.08
#> 69  0.40449438 0.47216260  0.1488447071  0.16 0.34 0.78   0.32 0.34
#> 70  0.40449438 0.47216260  0.1488447071  0.16 0.34 0.78   0.32 0.34
#> 71  0.15053763 0.19192970  0.2208106473  0.08 0.56 0.44   0.24 0.56
#> 72  0.15053763 0.19192970  0.2208106473  0.08 0.56 0.44   0.24 0.56
#> 73  0.27328925 0.26658985  0.0471040007  0.60 0.70 0.70   0.24 0.70
#> 74  0.44186047 0.36911845 -0.1914422722  0.20 0.86 0.44   0.26 0.86
#> 75  0.44186047 0.36911845 -0.1914422722  0.20 0.86 0.44   0.26 0.86
#> 76  0.16304348 0.20220935  0.1989436620  0.08 0.60 0.46   0.26 0.60
#> 77  0.16304348 0.20220935  0.1989436620  0.08 0.60 0.46   0.26 0.60
#> 78  0.16304348 0.20220935  0.1989436620  0.08 0.60 0.46   0.26 0.60
#> 79  0.36170213 0.38026256  0.0541429853  0.56 0.82 0.92   1.00 0.82
#> 80  0.36170213 0.38026256  0.0541429853  0.56 0.84 0.92   0.98 0.82
#> 81  0.48684211 0.49922091  0.0314136126  0.70 0.10 0.18   0.72 0.10
#> 82  0.48684211 0.49922091  0.0314136126  0.70 0.10 0.18   0.72 0.10
#> 83  0.48684211 0.49922091  0.0314136126  0.70 0.10 0.18   0.72 0.10
#> 84  0.35483871 0.31870736 -0.1080291971  0.38 0.78 0.88   0.84 0.82
#> 85  0.35483871 0.31870736 -0.1080291971  0.38 0.78 0.88   0.84 0.82
#> 86  0.35483871 0.31870736 -0.1080291971  0.42 0.78 0.88   0.84 0.78
#> 87  0.47727273 0.46875000 -0.0124688279  0.94 0.34 0.26   0.00 0.36
#> 88  0.47727273 0.46875000 -0.0124688279  0.94 0.34 0.26   0.00 0.36
#> 89  0.47727273 0.46875000 -0.0124688279  0.94 0.34 0.26   0.00 0.34
#> 90  0.47674419 0.49993240  0.0522164808  0.58 0.08 0.28   0.86 0.08
#> 91  0.47674419 0.49993240  0.0522164808  0.58 0.08 0.28   0.86 0.08
#> 92  0.38983051 0.40290147  0.0409777139  0.60 0.72 0.82   0.84 0.72
#> 93  0.38983051 0.40290147  0.0409777139  0.60 0.72 0.82   0.86 0.72
#> 94  0.38983051 0.40290147  0.0409777139  0.60 0.72 0.82   0.84 0.72
#> 95  0.45454545 0.42536157 -0.0629199756  0.64 0.68 0.32   0.92 0.60
#> 96  0.45454545 0.42536157 -0.0629199756  0.64 0.60 0.32   0.88 0.60
#> 97  0.45454545 0.42536157 -0.0629199756  0.64 0.60 0.32   0.88 0.60
#> 98  0.19791667 0.17833116 -0.1046511628  0.44 0.52 0.56   0.98 0.56
#> 99  0.19791667 0.17833116 -0.1046511628  0.46 0.56 0.54   0.98 0.56
#> 100 0.48148148 0.49512269  0.0337565810  0.70 0.26 0.22   0.62 0.20
#> 101 0.48148148 0.49512269  0.0337565810  0.64 0.26 0.20   0.62 0.20
#> 102 0.24210526 0.25867036  0.0693069307  0.54 0.70 0.62   0.66 0.70
#> 103 0.24210526 0.25867036  0.0693069307  0.54 0.70 0.62   0.66 0.68
#> 104 0.39080460 0.45871317  0.1536903039  0.14 0.44 0.80   0.32 0.44
#> 105 0.39080460 0.45871317  0.1536903039  0.14 0.46 0.80   0.32 0.44
#> 106 0.39080460 0.45871317  0.1536903039  0.14 0.46 0.80   0.30 0.46
#> 107 0.18181818 0.16528925 -0.0943396226  0.46 0.48 0.52   0.00 0.48
#> 108 0.18181818 0.16528925 -0.0943396226  0.46 0.48 0.52   0.00 0.48
#> 109 0.46153846 0.42603550 -0.0778443114  0.50 0.60 0.32   0.92 0.60
#> 110 0.46153846 0.42603550 -0.0778443114  0.50 0.60 0.32   0.92 0.60
#> 111 0.35955056 0.32268653 -0.1086614173  0.38 0.84 0.90   0.80 0.84
#> 112 0.35955056 0.32268653 -0.1086614173  0.38 0.84 0.90   0.80 0.84
#> 113 0.35955056 0.32268653 -0.1086614173  0.38 0.84 0.90   0.78 0.82
#> 114 0.30000000 0.25500000 -0.1710526316  0.22 0.66 0.78   0.38 0.66
#> 115 0.30000000 0.25500000 -0.1710526316  0.22 0.66 0.78   0.38 0.66
#> 116 0.36263736 0.36088637  0.0006729475  0.80 0.96 0.94   0.00 0.98
#> 117 0.36263736 0.36088637  0.0006729475  0.82 0.96 0.94   0.00 0.96
#> 118 0.36263736 0.36088637  0.0006729475  0.80 0.96 0.94   0.00 0.96
#> 119 0.41463415 0.49256395  0.1641881639  0.12 0.30 0.70   0.32 0.30
#> 120 0.41463415 0.49256395  0.1641881639  0.12 0.30 0.70   0.32 0.30
#> 121 0.08695652 0.08317580 -0.0400000000  0.72 0.28 0.32   0.50 0.28
#> 122 0.08695652 0.08317580 -0.0400000000  0.70 0.28 0.32   0.00 0.28
#> 123 0.40476190 0.43622449  0.0780790591  0.44 0.50 0.72   0.94 0.50
#> 124 0.40476190 0.43622449  0.0780790591  0.44 0.50 0.70   0.94 0.50
#> 125 0.40476190 0.43622449  0.0780790591  0.44 0.56 0.72   1.00 0.58
#> 126 0.45205479 0.42644023 -0.0531914894  0.68 0.60 0.40   0.76 0.60
#> 127 0.45205479 0.42644023 -0.0531914894  0.68 0.60 0.40   0.76 0.60
#> 128 0.45205479 0.42644023 -0.0531914894  0.68 0.60 0.40   0.76 0.60
#> 129 0.45783133 0.45463783 -0.0009637006  0.82 0.46 0.32   0.52 0.46
#> 130 0.45783133 0.45463783 -0.0009637006  0.82 0.46 0.32   0.52 0.46
#> 131 0.45783133 0.45463783 -0.0009637006  0.82 0.46 0.32   0.52 0.46
#> 132 0.40000000 0.44955017  0.1160631383  0.30 0.46 0.78   0.50 0.46
#> 133 0.40000000 0.44955017  0.1160631383  0.30 0.46 0.78   0.52 0.46
#> 134 0.40000000 0.44955017  0.1160631383  0.30 0.46 0.78   0.52 0.46
#> 135 0.08791209 0.08404781 -0.0404624277  0.70 0.28 0.32   0.50 0.28
#> 136 0.08791209 0.08404781 -0.0404624277  0.70 0.28 0.32   0.52 0.28
#> 137 0.08219178 0.07881404 -0.0359712230  0.76 0.26 0.30   0.00 0.26
#> 138 0.08219178 0.07881404 -0.0359712230  0.76 0.24 0.28   0.00 0.26
#> 139 0.07777778 0.09500000  0.1866840731  0.08 0.32 0.28   0.32 0.32
#> 140 0.07777778 0.09500000  0.1866840731  0.08 0.30 0.28   0.32 0.32
#> 141 0.07777778 0.09500000  0.1866840731  0.08 0.30 0.28   0.32 0.28
#> 142 0.47058824 0.40525951 -0.1554332875  0.22 0.72 0.30   0.32 0.72
#> 143 0.47058824 0.40525951 -0.1554332875  0.22 0.72 0.30   0.32 0.72
#> 144 0.47058824 0.40525951 -0.1554332875  0.22 0.72 0.30   0.32 0.72
#> 145 0.43181818 0.42536157 -0.0094656489  0.84 0.60 0.54   0.00 0.60
#> 146 0.43181818 0.42536157 -0.0094656489  0.84 0.60 0.58   0.00 0.60
#> 147 0.43181818 0.42536157 -0.0094656489  0.84 0.60 0.54   0.00 0.60
#> 148 0.08219178 0.07881404 -0.0359712230  0.76 0.26 0.30   0.00 0.26
#> 149 0.08219178 0.07881404 -0.0359712230  0.76 0.24 0.28   0.00 0.26
#> 150 0.28571429 0.27532906 -0.0322011469  0.78 0.72 0.76   0.52 0.70
#> 151 0.28571429 0.27532906 -0.0322011469  0.78 0.72 0.76   0.50 0.72
#> 152 0.28571429 0.27532906 -0.0322011469  0.78 0.72 0.76   0.52 0.72
#> 153 0.47500000 0.40968750 -0.1532846715  0.22 0.72 0.30   0.38 0.72
#> 154 0.47500000 0.40968750 -0.1532846715  0.22 0.72 0.30   0.38 0.72
#> 155 0.28571429 0.29752066  0.0462065031  0.60 0.74 0.76   0.86 0.74
#> 156 0.28571429 0.29752066  0.0462065031  0.60 0.74 0.76   0.86 0.74
#> 157 0.28571429 0.29752066  0.0462065031  0.60 0.74 0.76   0.86 0.74
#> 158 0.25274725 0.22080667 -0.1392405063  0.22 0.64 0.68   0.56 0.62
#> 159 0.25274725 0.22080667 -0.1392405063  0.26 0.64 0.68   0.56 0.64
#> 160 0.25274725 0.22080667 -0.1392405063  0.26 0.64 0.68   0.56 0.62
#> 161 0.48314607 0.46364095 -0.0364283758  0.74 0.44 0.18   0.72 0.42
#> 162 0.48314607 0.46364095 -0.0364283758  0.74 0.44 0.20   0.72 0.44
#> 163 0.16161616 0.14855627 -0.0828729282  0.48 0.46 0.46   0.00 0.46
#> 164 0.16161616 0.14855627 -0.0828729282  0.48 0.46 0.46   0.00 0.46
#> 165 0.16161616 0.14855627 -0.0828729282  0.48 0.46 0.46   0.00 0.46
#> 166 0.44318182 0.39146436 -0.1264940239  0.30 0.76 0.40   0.48 0.76
#> 167 0.44318182 0.39146436 -0.1264940239  0.30 0.76 0.40   0.48 0.76
#> 168 0.44318182 0.39146436 -0.1264940239  0.30 0.76 0.40   0.48 0.76
#> 169 0.34117647 0.34055363  0.0040883074  0.74 0.90 0.84   0.00 0.90
#> 170 0.34117647 0.34055363  0.0040883074  0.74 0.90 0.84   0.00 0.90
#> 171 0.34117647 0.34055363  0.0040883074  0.74 0.90 0.84   0.00 0.90
#> 172 0.27472527 0.25268688 -0.0817307692  0.48 0.64 0.70   0.86 0.64
#> 173 0.27472527 0.25268688 -0.0817307692  0.48 0.64 0.70   0.86 0.64
#> 174 0.42045455 0.41160899 -0.0157778479  0.96 0.72 0.70   0.52 0.70
#> 175 0.42045455 0.41160899 -0.0157778479  0.96 0.70 0.68   0.52 0.70
#> 176 0.60227273 0.49941891 -0.2004686280  0.16 0.10 0.02   0.26 0.10
#> 177 0.60227273 0.49941891 -0.2004686280  0.16 0.10 0.02   0.26 0.10
#> 178 0.52222222 0.42438272 -0.2252489748  0.14 0.68 0.02   0.16 0.68
#> 179 0.52222222 0.42438272 -0.2252489748  0.16 0.68 0.02   0.16 0.68
#> 180 0.52222222 0.42438272 -0.2252489748  0.16 0.68 0.02   0.16 0.70
#> 181 0.50000000 0.49830990  0.0024563319  0.76 0.10 0.12   0.58 0.10
#> 182 0.50000000 0.49830990  0.0024563319  0.76 0.12 0.12   0.58 0.10
#> 183 0.50000000 0.49830990  0.0024563319  0.76 0.12 0.12   0.58 0.10
#> 184 0.41666667 0.48611111  0.1497005988  0.16 0.30 0.70   0.38 0.30
#> 185 0.41666667 0.48611111  0.1497005988  0.16 0.30 0.70   0.38 0.30
#> 186 0.44285714 0.46316326  0.0510204082  0.60 0.44 0.44   0.76 0.42
#> 187 0.44285714 0.46316326  0.0510204082  0.60 0.44 0.44   0.76 0.42
#> 188 0.45454545 0.35123967 -0.2881355932  0.06 0.94 0.32   0.12 0.94
#> 189 0.45454545 0.35123967 -0.2881355932  0.06 0.94 0.32   0.12 0.94
#> 190 0.45454545 0.35123967 -0.2881355932  0.06 0.94 0.32   0.12 0.94
#> 191 0.25974026 0.38455051  0.3303964758  0.00 0.80 0.70   0.10 0.80
#> 192 0.25974026 0.38455051  0.3303964758  0.00 0.80 0.70   0.10 0.80
#> 193 0.25974026 0.38455051  0.3303964758  0.00 0.80 0.70   0.10 0.80
#> 194 0.47058824 0.35986159 -0.3023255814  0.04 0.96 0.30   0.02 0.96
#> 195 0.47058824 0.35986159 -0.3023255814  0.04 0.96 0.30   0.04 0.96
#> 196 0.47058824 0.35986159 -0.3023255814  0.04 0.96 0.30   0.04 0.96
#> 197 0.45348837 0.35066252 -0.2878787879  0.10 0.94 0.38   0.08 0.94
#> 198 0.45348837 0.35066252 -0.2878787879  0.10 0.94 0.36   0.08 0.90
#> 199 0.45348837 0.35066252 -0.2878787879  0.10 0.94 0.38   0.08 0.94
#> 200 0.50649351 0.37822567 -0.3333333333  0.02 0.84 0.08   0.02 0.84
#> 201 0.50649351 0.37822567 -0.3333333333  0.02 0.84 0.08   0.02 0.84
#> 202 0.50649351 0.37822567 -0.3333333333  0.02 0.84 0.08   0.02 0.84
#> 203 0.13541667 0.12624782 -0.0674157303  0.52 0.44 0.44   0.52 0.44
#> 204 0.13541667 0.12624782 -0.0674157303  0.52 0.44 0.44   0.50 0.44
#> 205 0.13541667 0.12624782 -0.0674157303  0.52 0.44 0.44   0.50 0.44
#> 206 0.38888889 0.36364198 -0.0638661202  0.60 0.88 0.84   0.76 0.88
#> 207 0.38888889 0.36364198 -0.0638661202  0.60 0.88 0.84   0.76 0.88
#> 208 0.38888889 0.36364198 -0.0638661202  0.60 0.88 0.88   0.78 0.88
#> 209 0.33333333 0.27777778 -0.1940298507  0.18 0.72 0.80   0.26 0.72
#> 210 0.33333333 0.27777778 -0.1940298507  0.18 0.74 0.80   0.28 0.74
#> 211 0.27586207 0.30017176  0.0867256637  0.36 0.78 0.74   0.70 0.78
#> 212 0.27586207 0.30017176  0.0867256637  0.36 0.78 0.74   0.70 0.74
#> 213 0.27586207 0.30017176  0.0867256637  0.40 0.78 0.74   0.70 0.78
#> 214 0.43023256 0.33768253 -0.2686567164  0.14 0.88 0.58   0.16 0.88
#> 215 0.43023256 0.33768253 -0.2686567164  0.14 0.88 0.58   0.16 0.88
#> 216 0.43023256 0.33768253 -0.2686567164  0.14 0.88 0.58   0.16 0.88
#> 217 0.62068966 0.47615273 -0.2982946603  0.06 0.34 0.00   0.12 0.32
#> 218 0.62068966 0.47615273 -0.2982946603  0.04 0.32 0.00   0.12 0.32
#> 219 0.10526316 0.09972299 -0.0502793296  0.68 0.32 0.38   0.00 0.32
#> 220 0.10526316 0.09972299 -0.0502793296  0.68 0.32 0.38   0.00 0.32
#> 221 0.10526316 0.09972299 -0.0502793296  0.68 0.32 0.38   0.00 0.32
#> 222 0.11956522 0.11241730 -0.0581395349  0.64 0.36 0.42   0.00 0.36
#> 223 0.11956522 0.11241730 -0.0581395349  0.66 0.34 0.42   0.00 0.36
#> 224 0.37500000 0.36218750 -0.0290924881  0.86 0.98 1.00   0.00 0.98
#> 225 0.37500000 0.36218750 -0.0290924881  0.82 0.98 1.00   0.00 0.98
#> 226 0.37500000 0.36218750 -0.0290924881  0.82 0.98 1.00   0.00 0.98
#> 227 0.16091954 0.14797199 -0.0817610063  0.48 0.46 0.46   0.50 0.46
#> 228 0.16091954 0.14797199 -0.0817610063  0.48 0.46 0.46   0.50 0.46
#> 229 0.16091954 0.14797199 -0.0817610063  0.48 0.46 0.46   0.52 0.46
#> 230 0.06451613 0.06243496 -0.0279329609  0.94 0.18 0.18   0.00 0.18
#> 231 0.06451613 0.06243496 -0.0279329609  0.98 0.18 0.18   0.00 0.18
#> 232 0.34939759 0.38387284  0.0958174905  0.34 0.82 0.84   0.66 0.82
#> 233 0.34939759 0.38387284  0.0958174905  0.34 0.82 0.84   0.66 0.82
#> 234 0.34939759 0.38387284  0.0958174905  0.34 0.82 0.84   0.66 0.82
#> 235 0.11494253 0.10833663 -0.0552147239  0.66 0.34 0.40   0.00 0.34
#> 236 0.11494253 0.10833663 -0.0552147239  0.68 0.34 0.38   0.00 0.34
#> 237 0.11494253 0.10833663 -0.0552147239  0.66 0.34 0.40   0.00 0.34
#> 238 0.05154639 0.05021788 -0.0212765957  0.98 0.12 0.12   0.00 0.12
#> 239 0.05154639 0.05021788 -0.0212765957  0.98 0.12 0.12   0.00 0.12
#> 240 0.05154639 0.05021788 -0.0212765957  0.98 0.12 0.12   0.00 0.12
#> 241 0.37804878 0.38689768  0.0290023202  0.72 0.80 1.00   0.76 0.80
#> 242 0.37804878 0.38689768  0.0290023202  0.72 0.80 1.00   0.76 0.80
#> 243 0.37804878 0.38689768  0.0290023202  0.72 0.80 1.00   0.76 0.80
#> 244 0.48780488 0.43307555 -0.1203319502  0.32 0.58 0.16   0.50 0.58
#> 245 0.48780488 0.43307555 -0.1203319502  0.32 0.60 0.16   0.50 0.58
#> 246 0.48780488 0.43307555 -0.1203319502  0.32 0.60 0.16   0.50 0.58
#> 247 0.22471910 0.26511804  0.1578947368  0.12 0.70 0.60   0.38 0.70
#> 248 0.22471910 0.26511804  0.1578947368  0.12 0.70 0.60   0.38 0.70
#> 249 0.51250000 0.49367187 -0.0318572794  0.82 0.30 0.04   0.72 0.30
#> 250 0.51250000 0.49367187 -0.0318572794  0.82 0.26 0.04   0.68 0.28
#> 251 0.51250000 0.49367187 -0.0318572794  0.82 0.30 0.04   0.68 0.28
#> 252 0.47777778 0.49697531  0.0442057942  0.60 0.14 0.26   0.62 0.14
#> 253 0.47777778 0.49697531  0.0442057942  0.60 0.14 0.26   0.62 0.14
#> 254 0.47777778 0.49697531  0.0442057942  0.60 0.14 0.26   0.62 0.14
#> 255 0.38554217 0.49644361  0.2291421857  0.02 0.14 0.94   0.24 0.14
#> 256 0.38554217 0.49644361  0.2291421857  0.02 0.14 0.94   0.24 0.14
#> 257 0.38554217 0.49644361  0.2291421857  0.02 0.18 1.00   0.24 0.14
#> 258 0.42857143 0.46570295  0.0856793146  0.40 0.38 0.62   0.74 0.38
#> 259 0.42857143 0.46570295  0.0856793146  0.40 0.38 0.62   0.74 0.38
#> 260 0.26153846 0.22733728 -0.1428571429  0.22 0.64 0.70   0.90 0.64
#> 261 0.26153846 0.22733728 -0.1428571429  0.22 0.64 0.70   0.90 0.64
#>     p_fis_overall_BY p_fis_byfacet_BY p_he_overall_BY p_he_byfacet_BY
#> 1                  1                1               1               1
#> 2                  1                1               1               1
#> 3                  1                1               1               1
#> 4                  1                1               1               1
#> 5                  1                1               1               1
#> 6                  1                1               1               1
#> 7                  1                1               1               1
#> 8                  1                1               1               1
#> 9                  1                1               1               1
#> 10                 1                1               1               1
#> 11                 1                1               1               1
#> 12                 1                1               1               1
#> 13                 1                1               1               1
#> 14                 1                1               1               1
#> 15                 1                1               1               1
#> 16                 1                1               1               1
#> 17                 1                1               1               1
#> 18                 1                1               1               1
#> 19                 1                1               1               1
#> 20                 1                1               1               1
#> 21                 1                1               1               1
#> 22                 1                1               1               1
#> 23                 1                1               1               1
#> 24                 1                1               1               1
#> 25                 1                1               1               1
#> 26                 1                1               1               1
#> 27                 1                1               1               1
#> 28                 1                1               1               1
#> 29                 1                1               1               1
#> 30                 1                1               1               1
#> 31                 1                1               1               1
#> 32                 1                1               1               1
#> 33                 1                1               1               1
#> 34                 1                1               1               1
#> 35                 1                1               1               1
#> 36                 1                1               1               1
#> 37                 0                0               1               1
#> 38                 0                0               1               1
#> 39                 1                1               1               1
#> 40                 1                1               1               1
#> 41                 1                1               1               1
#> 42                 1                1               1               1
#> 43                 1                1               1               1
#> 44                 1                1               1               1
#> 45                 1                1               1               1
#> 46                 1                1               1               1
#> 47                 1                1               1               1
#> 48                 1                1               1               1
#> 49                 1                1               1               1
#> 50                 1                1               1               1
#> 51                 1                1               1               1
#> 52                 1                1               1               1
#> 53                 1                1               1               1
#> 54                 1                1               1               1
#> 55                 1                1               1               1
#> 56                 1                1               1               1
#> 57                 1                1               1               1
#> 58                 1                1               1               1
#> 59                 1                1               1               1
#> 60                 1                1               1               1
#> 61                 1                1               1               1
#> 62                 1                1               1               1
#> 63                 1                1               1               1
#> 64                 1                1               0               0
#> 65                 1                1               0               0
#> 66                 1                1               0               0
#> 67                 1                1               1               1
#> 68                 1                1               1               1
#> 69                 1                1               1               1
#> 70                 1                1               1               1
#> 71                 1                1               1               1
#> 72                 1                1               1               1
#> 73                 1                1               1               1
#> 74                 1                1               1               1
#> 75                 1                1               1               1
#> 76                 1                1               1               1
#> 77                 1                1               1               1
#> 78                 1                1               1               1
#> 79                 1                1               1               1
#> 80                 1                1               1               1
#> 81                 1                1               1               1
#> 82                 1                1               1               1
#> 83                 1                1               1               1
#> 84                 1                1               1               1
#> 85                 1                1               1               1
#> 86                 1                1               1               1
#> 87                 1                1               1               1
#> 88                 1                1               1               1
#> 89                 1                1               1               1
#> 90                 1                1               1               1
#> 91                 1                1               1               1
#> 92                 1                1               1               1
#> 93                 1                1               1               1
#> 94                 1                1               1               1
#> 95                 1                1               1               1
#> 96                 1                1               1               1
#> 97                 1                1               1               1
#> 98                 1                1               1               1
#> 99                 1                1               1               1
#> 100                1                1               1               1
#> 101                1                1               1               1
#> 102                1                1               1               1
#> 103                1                1               1               1
#> 104                1                1               1               1
#> 105                1                1               1               1
#> 106                1                1               1               1
#> 107                1                1               1               1
#> 108                1                1               1               1
#> 109                1                1               1               1
#> 110                1                1               1               1
#> 111                1                1               1               1
#> 112                1                1               1               1
#> 113                1                1               1               1
#> 114                1                1               1               1
#> 115                1                1               1               1
#> 116                1                1               1               1
#> 117                1                1               1               1
#> 118                1                1               1               1
#> 119                1                1               1               1
#> 120                1                1               1               1
#> 121                1                1               1               1
#> 122                1                1               1               1
#> 123                1                1               1               1
#> 124                1                1               1               1
#> 125                1                1               1               1
#> 126                1                1               1               1
#> 127                1                1               1               1
#> 128                1                1               1               1
#> 129                1                1               1               1
#> 130                1                1               1               1
#> 131                1                1               1               1
#> 132                1                1               1               1
#> 133                1                1               1               1
#> 134                1                1               1               1
#> 135                1                1               1               1
#> 136                1                1               1               1
#> 137                1                1               1               1
#> 138                1                1               1               1
#> 139                1                1               1               1
#> 140                1                1               1               1
#> 141                1                1               1               1
#> 142                1                1               1               1
#> 143                1                1               1               1
#> 144                1                1               1               1
#> 145                1                1               1               1
#> 146                1                1               1               1
#> 147                1                1               1               1
#> 148                1                1               1               1
#> 149                1                1               1               1
#> 150                1                1               1               1
#> 151                1                1               1               1
#> 152                1                1               1               1
#> 153                1                1               1               1
#> 154                1                1               1               1
#> 155                1                1               1               1
#> 156                1                1               1               1
#> 157                1                1               1               1
#> 158                1                1               1               1
#> 159                1                1               1               1
#> 160                1                1               1               1
#> 161                1                1               1               1
#> 162                1                1               1               1
#> 163                1                1               1               1
#> 164                1                1               1               1
#> 165                1                1               1               1
#> 166                1                1               1               1
#> 167                1                1               1               1
#> 168                1                1               1               1
#> 169                1                1               1               1
#> 170                1                1               1               1
#> 171                1                1               1               1
#> 172                1                1               1               1
#> 173                1                1               1               1
#> 174                1                1               1               1
#> 175                1                1               1               1
#> 176                1                1               1               1
#> 177                1                1               1               1
#> 178                1                1               1               1
#> 179                1                1               1               1
#> 180                1                1               1               1
#> 181                1                1               1               1
#> 182                1                1               1               1
#> 183                1                1               1               1
#> 184                1                1               1               1
#> 185                1                1               1               1
#> 186                1                1               1               1
#> 187                1                1               1               1
#> 188                1                1               1               1
#> 189                1                1               1               1
#> 190                1                1               1               1
#> 191                0                0               1               1
#> 192                0                0               1               1
#> 193                0                0               1               1
#> 194                1                1               1               1
#> 195                1                1               1               1
#> 196                1                1               1               1
#> 197                1                1               1               1
#> 198                1                1               1               1
#> 199                1                1               1               1
#> 200                1                1               1               1
#> 201                1                1               1               1
#> 202                1                1               1               1
#> 203                1                1               1               1
#> 204                1                1               1               1
#> 205                1                1               1               1
#> 206                1                1               1               1
#> 207                1                1               1               1
#> 208                1                1               1               1
#> 209                1                1               1               1
#> 210                1                1               1               1
#> 211                1                1               1               1
#> 212                1                1               1               1
#> 213                1                1               1               1
#> 214                1                1               1               1
#> 215                1                1               1               1
#> 216                1                1               1               1
#> 217                1                1               1               1
#> 218                1                1               1               1
#> 219                1                1               1               1
#> 220                1                1               1               1
#> 221                1                1               1               1
#> 222                1                1               1               1
#> 223                1                1               1               1
#> 224                1                1               1               1
#> 225                1                1               1               1
#> 226                1                1               1               1
#> 227                1                1               1               1
#> 228                1                1               1               1
#> 229                1                1               1               1
#> 230                1                1               1               1
#> 231                1                1               1               1
#> 232                1                1               1               1
#> 233                1                1               1               1
#> 234                1                1               1               1
#> 235                1                1               1               1
#> 236                1                1               1               1
#> 237                1                1               1               1
#> 238                1                1               1               1
#> 239                1                1               1               1
#> 240                1                1               1               1
#> 241                1                1               1               1
#> 242                1                1               1               1
#> 243                1                1               1               1
#> 244                1                1               1               1
#> 245                1                1               1               1
#> 246                1                1               1               1
#> 247                1                1               1               1
#> 248                1                1               1               1
#> 249                1                1               1               1
#> 250                1                1               1               1
#> 251                1                1               1               1
#> 252                1                1               1               1
#> 253                1                1               1               1
#> 254                1                1               1               1
#> 255                1                1               1               1
#> 256                1                1               1               1
#> 257                1                1               1               1
#> 258                1                1               1               1
#> 259                1                1               1               1
#> 260                1                1               1               1
#> 261                1                1               1               1
#>     p_ho_overall_BY p_ho_byfacet_BY p_pHWE_overall_BY p_pHWE_byfacet_BY
#> 1                 1               1         0.0000000         0.0000000
#> 2                 1               1         0.0000000         0.0000000
#> 3                 1               1         0.0000000         0.0000000
#> 4                 1               1         0.0000000         0.0000000
#> 5                 1               1         0.0000000         0.0000000
#> 6                 1               1         1.0000000         1.0000000
#> 7                 1               1         1.0000000         1.0000000
#> 8                 1               1         1.0000000         1.0000000
#> 9                 1               1         1.0000000         1.0000000
#> 10                1               1         1.0000000         1.0000000
#> 11                1               1         1.0000000         1.0000000
#> 12                1               1         0.0000000         0.0000000
#> 13                1               1         0.0000000         0.0000000
#> 14                1               1         1.0000000         1.0000000
#> 15                1               1         1.0000000         1.0000000
#> 16                1               1         1.0000000         1.0000000
#> 17                1               1         1.0000000         1.0000000
#> 18                1               1         1.0000000         1.0000000
#> 19                1               1         1.0000000         1.0000000
#> 20                1               1         0.0000000         0.0000000
#> 21                1               1         0.0000000         0.0000000
#> 22                1               1         0.0000000         0.0000000
#> 23                1               1         0.0000000         0.0000000
#> 24                1               1         0.0000000         0.0000000
#> 25                1               1         1.0000000         1.0000000
#> 26                1               1         1.0000000         1.0000000
#> 27                1               1         1.0000000         1.0000000
#> 28                1               1         1.0000000         1.0000000
#> 29                1               1         1.0000000         1.0000000
#> 30                1               1         1.0000000         1.0000000
#> 31                1               1         0.0000000         0.0000000
#> 32                1               1         0.0000000         0.0000000
#> 33                1               1         0.0000000         0.0000000
#> 34                1               1         1.0000000         1.0000000
#> 35                1               1         1.0000000         1.0000000
#> 36                1               1         1.0000000         1.0000000
#> 37                1               1         0.0000000         0.0000000
#> 38                1               1         0.0000000         0.0000000
#> 39                1               1         1.0000000         1.0000000
#> 40                1               1         1.0000000         1.0000000
#> 41                1               1         1.0000000         1.0000000
#> 42                1               1         0.0000000         0.0000000
#> 43                1               1         0.0000000         0.0000000
#> 44                1               1         0.0000000         0.0000000
#> 45                1               1         0.0000000         0.0000000
#> 46                1               1         0.0000000         0.0000000
#> 47                1               1         0.0000000         0.0000000
#> 48                1               1         1.0000000         1.0000000
#> 49                1               1         1.0000000         1.0000000
#> 50                1               1         1.0000000         1.0000000
#> 51                1               1         1.0000000         1.0000000
#> 52                1               1         1.0000000         1.0000000
#> 53                1               1         1.0000000         1.0000000
#> 54                1               1         1.0000000         1.0000000
#> 55                1               1         1.0000000         1.0000000
#> 56                1               1         1.0000000         1.0000000
#> 57                1               1         1.0000000         1.0000000
#> 58                1               1         1.0000000         1.0000000
#> 59                1               1         1.0000000         1.0000000
#> 60                1               1         1.0000000         1.0000000
#> 61                1               1         0.0000000         0.0000000
#> 62                1               1         0.0000000         0.0000000
#> 63                1               1         0.0000000         0.0000000
#> 64                1               1         1.0000000         1.0000000
#> 65                1               1         1.0000000         1.0000000
#> 66                1               1         1.0000000         1.0000000
#> 67                1               1         0.0000000         0.0000000
#> 68                1               1         0.0000000         0.0000000
#> 69                1               1         1.0000000         1.0000000
#> 70                1               1         1.0000000         1.0000000
#> 71                1               1         1.0000000         1.0000000
#> 72                1               1         1.0000000         1.0000000
#> 73                1               1         1.0000000         1.0000000
#> 74                1               1         1.0000000         1.0000000
#> 75                1               1         1.0000000         1.0000000
#> 76                1               1         1.0000000         1.0000000
#> 77                1               1         1.0000000         1.0000000
#> 78                1               1         1.0000000         1.0000000
#> 79                1               1         1.0000000         1.0000000
#> 80                1               1         1.0000000         1.0000000
#> 81                1               1         1.0000000         1.0000000
#> 82                1               1         1.0000000         1.0000000
#> 83                1               1         1.0000000         1.0000000
#> 84                1               1         1.0000000         1.0000000
#> 85                1               1         1.0000000         1.0000000
#> 86                1               1         1.0000000         1.0000000
#> 87                1               1         0.0000000         0.0000000
#> 88                1               1         0.0000000         0.0000000
#> 89                1               1         0.0000000         0.0000000
#> 90                1               1         1.0000000         1.0000000
#> 91                1               1         1.0000000         1.0000000
#> 92                1               1         1.0000000         1.0000000
#> 93                1               1         1.0000000         1.0000000
#> 94                1               1         1.0000000         1.0000000
#> 95                1               1         1.0000000         1.0000000
#> 96                1               1         1.0000000         1.0000000
#> 97                1               1         1.0000000         1.0000000
#> 98                1               1         1.0000000         1.0000000
#> 99                1               1         1.0000000         1.0000000
#> 100               1               1         1.0000000         1.0000000
#> 101               1               1         1.0000000         1.0000000
#> 102               1               1         1.0000000         1.0000000
#> 103               1               1         1.0000000         1.0000000
#> 104               1               1         1.0000000         1.0000000
#> 105               1               1         1.0000000         1.0000000
#> 106               1               1         1.0000000         1.0000000
#> 107               1               1         0.0000000         0.0000000
#> 108               1               1         0.0000000         0.0000000
#> 109               1               1         1.0000000         1.0000000
#> 110               1               1         1.0000000         1.0000000
#> 111               1               1         1.0000000         1.0000000
#> 112               1               1         1.0000000         1.0000000
#> 113               1               1         1.0000000         1.0000000
#> 114               1               1         1.0000000         1.0000000
#> 115               1               1         1.0000000         1.0000000
#> 116               1               1         0.0000000         0.0000000
#> 117               1               1         0.0000000         0.0000000
#> 118               1               1         0.0000000         0.0000000
#> 119               1               1         1.0000000         1.0000000
#> 120               1               1         1.0000000         1.0000000
#> 121               1               1         1.0000000         1.0000000
#> 122               1               1         0.0000000         0.0000000
#> 123               1               1         1.0000000         1.0000000
#> 124               1               1         1.0000000         1.0000000
#> 125               1               1         1.0000000         1.0000000
#> 126               1               1         1.0000000         1.0000000
#> 127               1               1         1.0000000         1.0000000
#> 128               1               1         1.0000000         1.0000000
#> 129               1               1         1.0000000         1.0000000
#> 130               1               1         1.0000000         1.0000000
#> 131               1               1         1.0000000         1.0000000
#> 132               1               1         1.0000000         1.0000000
#> 133               1               1         1.0000000         1.0000000
#> 134               1               1         1.0000000         1.0000000
#> 135               1               1         1.0000000         1.0000000
#> 136               1               1         1.0000000         1.0000000
#> 137               1               1         0.0000000         0.0000000
#> 138               1               1         0.0000000         0.0000000
#> 139               1               1         1.0000000         1.0000000
#> 140               1               1         1.0000000         1.0000000
#> 141               1               1         1.0000000         1.0000000
#> 142               1               1         1.0000000         1.0000000
#> 143               1               1         1.0000000         1.0000000
#> 144               1               1         1.0000000         1.0000000
#> 145               1               1         0.0000000         0.0000000
#> 146               1               1         0.0000000         0.0000000
#> 147               1               1         0.0000000         0.0000000
#> 148               1               1         0.0000000         0.0000000
#> 149               1               1         0.0000000         0.0000000
#> 150               1               1         1.0000000         1.0000000
#> 151               1               1         1.0000000         1.0000000
#> 152               1               1         1.0000000         1.0000000
#> 153               1               1         1.0000000         1.0000000
#> 154               1               1         1.0000000         1.0000000
#> 155               1               1         1.0000000         1.0000000
#> 156               1               1         1.0000000         1.0000000
#> 157               1               1         1.0000000         1.0000000
#> 158               1               1         1.0000000         1.0000000
#> 159               1               1         1.0000000         1.0000000
#> 160               1               1         1.0000000         1.0000000
#> 161               1               1         1.0000000         1.0000000
#> 162               1               1         1.0000000         1.0000000
#> 163               1               1         0.0000000         0.0000000
#> 164               1               1         0.0000000         0.0000000
#> 165               1               1         0.0000000         0.0000000
#> 166               1               1         1.0000000         1.0000000
#> 167               1               1         1.0000000         1.0000000
#> 168               1               1         1.0000000         1.0000000
#> 169               1               1         0.0000000         0.0000000
#> 170               1               1         0.0000000         0.0000000
#> 171               1               1         0.0000000         0.0000000
#> 172               1               1         1.0000000         1.0000000
#> 173               1               1         1.0000000         1.0000000
#> 174               1               1         1.0000000         1.0000000
#> 175               1               1         1.0000000         1.0000000
#> 176               1               1         1.0000000         1.0000000
#> 177               1               1         1.0000000         1.0000000
#> 178               1               1         1.0000000         1.0000000
#> 179               1               1         1.0000000         1.0000000
#> 180               1               1         1.0000000         1.0000000
#> 181               1               1         1.0000000         1.0000000
#> 182               1               1         1.0000000         1.0000000
#> 183               1               1         1.0000000         1.0000000
#> 184               1               1         1.0000000         1.0000000
#> 185               1               1         1.0000000         1.0000000
#> 186               1               1         1.0000000         1.0000000
#> 187               1               1         1.0000000         1.0000000
#> 188               1               1         1.0000000         1.0000000
#> 189               1               1         1.0000000         1.0000000
#> 190               1               1         1.0000000         1.0000000
#> 191               1               1         1.0000000         1.0000000
#> 192               1               1         1.0000000         1.0000000
#> 193               1               1         1.0000000         1.0000000
#> 194               1               1         0.4581408         0.4581408
#> 195               1               1         0.8908293         0.8908293
#> 196               1               1         0.8908293         0.8908293
#> 197               1               1         1.0000000         1.0000000
#> 198               1               1         1.0000000         1.0000000
#> 199               1               1         1.0000000         1.0000000
#> 200               1               1         0.4581408         0.4581408
#> 201               1               1         0.4581408         0.4581408
#> 202               1               1         0.4581408         0.4581408
#> 203               1               1         1.0000000         1.0000000
#> 204               1               1         1.0000000         1.0000000
#> 205               1               1         1.0000000         1.0000000
#> 206               1               1         1.0000000         1.0000000
#> 207               1               1         1.0000000         1.0000000
#> 208               1               1         1.0000000         1.0000000
#> 209               1               1         1.0000000         1.0000000
#> 210               1               1         1.0000000         1.0000000
#> 211               1               1         1.0000000         1.0000000
#> 212               1               1         1.0000000         1.0000000
#> 213               1               1         1.0000000         1.0000000
#> 214               1               1         1.0000000         1.0000000
#> 215               1               1         1.0000000         1.0000000
#> 216               1               1         1.0000000         1.0000000
#> 217               0               0         1.0000000         1.0000000
#> 218               0               0         1.0000000         1.0000000
#> 219               1               1         0.0000000         0.0000000
#> 220               1               1         0.0000000         0.0000000
#> 221               1               1         0.0000000         0.0000000
#> 222               1               1         0.0000000         0.0000000
#> 223               1               1         0.0000000         0.0000000
#> 224               1               1         0.0000000         0.0000000
#> 225               1               1         0.0000000         0.0000000
#> 226               1               1         0.0000000         0.0000000
#> 227               1               1         1.0000000         1.0000000
#> 228               1               1         1.0000000         1.0000000
#> 229               1               1         1.0000000         1.0000000
#> 230               1               1         0.0000000         0.0000000
#> 231               1               1         0.0000000         0.0000000
#> 232               1               1         1.0000000         1.0000000
#> 233               1               1         1.0000000         1.0000000
#> 234               1               1         1.0000000         1.0000000
#> 235               1               1         0.0000000         0.0000000
#> 236               1               1         0.0000000         0.0000000
#> 237               1               1         0.0000000         0.0000000
#> 238               1               1         0.0000000         0.0000000
#> 239               1               1         0.0000000         0.0000000
#> 240               1               1         0.0000000         0.0000000
#> 241               1               1         1.0000000         1.0000000
#> 242               1               1         1.0000000         1.0000000
#> 243               1               1         1.0000000         1.0000000
#> 244               1               1         1.0000000         1.0000000
#> 245               1               1         1.0000000         1.0000000
#> 246               1               1         1.0000000         1.0000000
#> 247               1               1         1.0000000         1.0000000
#> 248               1               1         1.0000000         1.0000000
#> 249               1               1         1.0000000         1.0000000
#> 250               1               1         1.0000000         1.0000000
#> 251               1               1         1.0000000         1.0000000
#> 252               1               1         1.0000000         1.0000000
#> 253               1               1         1.0000000         1.0000000
#> 254               1               1         1.0000000         1.0000000
#> 255               1               1         1.0000000         1.0000000
#> 256               1               1         1.0000000         1.0000000
#> 257               1               1         1.0000000         1.0000000
#> 258               1               1         1.0000000         1.0000000
#> 259               1               1         1.0000000         1.0000000
#> 260               1               1         1.0000000         1.0000000
#> 261               1               1         1.0000000         1.0000000
#>     p_pi_overall_BY p_pi_byfacet_BY
#> 1                 1               1
#> 2                 1               1
#> 3                 1               1
#> 4                 1               1
#> 5                 1               1
#> 6                 1               1
#> 7                 1               1
#> 8                 1               1
#> 9                 1               1
#> 10                1               1
#> 11                1               1
#> 12                1               1
#> 13                1               1
#> 14                1               1
#> 15                1               1
#> 16                1               1
#> 17                1               1
#> 18                1               1
#> 19                1               1
#> 20                1               1
#> 21                1               1
#> 22                1               1
#> 23                1               1
#> 24                1               1
#> 25                1               1
#> 26                1               1
#> 27                1               1
#> 28                1               1
#> 29                1               1
#> 30                1               1
#> 31                1               1
#> 32                1               1
#> 33                1               1
#> 34                1               1
#> 35                1               1
#> 36                1               1
#> 37                1               1
#> 38                1               1
#> 39                1               1
#> 40                1               1
#> 41                1               1
#> 42                1               1
#> 43                1               1
#> 44                1               1
#> 45                1               1
#> 46                1               1
#> 47                1               1
#> 48                1               1
#> 49                1               1
#> 50                1               1
#> 51                1               1
#> 52                1               1
#> 53                1               1
#> 54                1               1
#> 55                1               1
#> 56                1               1
#> 57                1               1
#> 58                1               1
#> 59                1               1
#> 60                1               1
#> 61                1               1
#> 62                1               1
#> 63                1               1
#> 64                0               0
#> 65                0               0
#> 66                0               0
#> 67                1               1
#> 68                1               1
#> 69                1               1
#> 70                1               1
#> 71                1               1
#> 72                1               1
#> 73                1               1
#> 74                1               1
#> 75                1               1
#> 76                1               1
#> 77                1               1
#> 78                1               1
#> 79                1               1
#> 80                1               1
#> 81                1               1
#> 82                1               1
#> 83                1               1
#> 84                1               1
#> 85                1               1
#> 86                1               1
#> 87                1               1
#> 88                1               1
#> 89                1               1
#> 90                1               1
#> 91                1               1
#> 92                1               1
#> 93                1               1
#> 94                1               1
#> 95                1               1
#> 96                1               1
#> 97                1               1
#> 98                1               1
#> 99                1               1
#> 100               1               1
#> 101               1               1
#> 102               1               1
#> 103               1               1
#> 104               1               1
#> 105               1               1
#> 106               1               1
#> 107               1               1
#> 108               1               1
#> 109               1               1
#> 110               1               1
#> 111               1               1
#> 112               1               1
#> 113               1               1
#> 114               1               1
#> 115               1               1
#> 116               1               1
#> 117               1               1
#> 118               1               1
#> 119               1               1
#> 120               1               1
#> 121               1               1
#> 122               1               1
#> 123               1               1
#> 124               1               1
#> 125               1               1
#> 126               1               1
#> 127               1               1
#> 128               1               1
#> 129               1               1
#> 130               1               1
#> 131               1               1
#> 132               1               1
#> 133               1               1
#> 134               1               1
#> 135               1               1
#> 136               1               1
#> 137               1               1
#> 138               1               1
#> 139               1               1
#> 140               1               1
#> 141               1               1
#> 142               1               1
#> 143               1               1
#> 144               1               1
#> 145               1               1
#> 146               1               1
#> 147               1               1
#> 148               1               1
#> 149               1               1
#> 150               1               1
#> 151               1               1
#> 152               1               1
#> 153               1               1
#> 154               1               1
#> 155               1               1
#> 156               1               1
#> 157               1               1
#> 158               1               1
#> 159               1               1
#> 160               1               1
#> 161               1               1
#> 162               1               1
#> 163               1               1
#> 164               1               1
#> 165               1               1
#> 166               1               1
#> 167               1               1
#> 168               1               1
#> 169               1               1
#> 170               1               1
#> 171               1               1
#> 172               1               1
#> 173               1               1
#> 174               1               1
#> 175               1               1
#> 176               1               1
#> 177               1               1
#> 178               1               1
#> 179               1               1
#> 180               1               1
#> 181               1               1
#> 182               1               1
#> 183               1               1
#> 184               1               1
#> 185               1               1
#> 186               1               1
#> 187               1               1
#> 188               1               1
#> 189               1               1
#> 190               1               1
#> 191               1               1
#> 192               1               1
#> 193               1               1
#> 194               1               1
#> 195               1               1
#> 196               1               1
#> 197               1               1
#> 198               1               1
#> 199               1               1
#> 200               1               1
#> 201               1               1
#> 202               1               1
#> 203               1               1
#> 204               1               1
#> 205               1               1
#> 206               1               1
#> 207               1               1
#> 208               1               1
#> 209               1               1
#> 210               1               1
#> 211               1               1
#> 212               1               1
#> 213               1               1
#> 214               1               1
#> 215               1               1
#> 216               1               1
#> 217               1               1
#> 218               1               1
#> 219               1               1
#> 220               1               1
#> 221               1               1
#> 222               1               1
#> 223               1               1
#> 224               1               1
#> 225               1               1
#> 226               1               1
#> 227               1               1
#> 228               1               1
#> 229               1               1
#> 230               1               1
#> 231               1               1
#> 232               1               1
#> 233               1               1
#> 234               1               1
#> 235               1               1
#> 236               1               1
#> 237               1               1
#> 238               1               1
#> 239               1               1
#> 240               1               1
#> 241               1               1
#> 242               1               1
#> 243               1               1
#> 244               1               1
#> 245               1               1
#> 246               1               1
#> 247               1               1
#> 248               1               1
#> 249               1               1
#> 250               1               1
#> 251               1               1
#> 252               1               1
#> 253               1               1
#> 254               1               1
#> 255               1               1
#> 256               1               1
#> 257               1               1
#> 258               1               1
#> 259               1               1
#> 260               1               1
#> 261               1               1
```
