# Calculate gaussian weights.

Get the gaussian weight of a statistic given its position relative to
the center and scaling factor sigma.

## Usage

``` r
gaussian_weight(p, c, s)
```

## Arguments

- p:

  numeric. Position of value to smooth.

- c:

  numeric. Position of center of window.

- s:

  numeric. Sigma, scaling factor.

## Value

A numeric value, the relative weight of the point.
