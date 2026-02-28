# Peak-to-peak (range): max(x) - min(x)

Computes the range of values along an axis: `max(x) - min(x)`. The
result is always nonnegative.

## Usage

``` r
ptp(x, axis = NULL, keepdims = FALSE)
```

## Arguments

- x:

  An Expression or numeric value.

- axis:

  NULL (all), 0 (columns), or 1 (rows).

- keepdims:

  Logical; keep reduced dimension?

## Value

An Expression representing max(x) - min(x).
