# (Weighted) geometric mean of a vector

(Weighted) geometric mean of a vector

## Usage

``` r
geo_mean(x, p = NULL, max_denom = 1024L, approx = TRUE)
```

## Arguments

- x:

  An Expression (vector)

- p:

  Numeric weight vector (default: uniform)

- max_denom:

  Maximum denominator for rational approximation

- approx:

  If TRUE (default), use SOC approximation. If FALSE, use exact power
  cone.

## Value

A GeoMean or GeoMeanApprox atom
