# (Weighted) geometric mean

The (weighted) geometric mean of the entries of `x`, of all of them or
of each row or column.

## Usage

``` r
geo_mean(
  x,
  p = NULL,
  max_denom = 1024L,
  approx = TRUE,
  axis = NULL,
  keepdims = FALSE
)
```

## Arguments

- x:

  An Expression of any shape. With `axis = NULL` its entries are taken
  in column-major order, so a row vector, a column vector and a matrix
  holding the same entries give the same result.

- p:

  Numeric weight vector (default: uniform), one weight per entry being
  reduced: `length(x)` entries with `axis = NULL`, `ncol(x)` with
  `axis = 1`, `nrow(x)` with `axis = 2`. Entries with zero weight take
  no part.

- max_denom:

  Maximum denominator for rational approximation

- approx:

  If TRUE (default), use SOC approximation. If FALSE, use exact power
  cone.

- axis:

  `NULL` (default) for one geometric mean of all entries, `1` for one
  per row, `2` for one per column.

- keepdims:

  Keep the reduced dimension with length one.

## Value

A GeoMean or GeoMeanApprox atom: a scalar when `axis = NULL`, otherwise
one geometric mean per row or column.
