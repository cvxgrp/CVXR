# Compute kth Order Differences of an Expression

Takes in an expression and returns an expression with the kth order
differences along the given axis. The output shape is the same as the
input except the size along the specified axis is reduced by k.

## Usage

``` r
cvxr_diff(x, k = 1L, axis = 2L)
```

## Arguments

- x:

  An Expression or numeric value.

- k:

  Integer. The number of times values are differenced. Default is 1.
  (Mapped from R's `lag` argument in diff.default; use `differences` for
  repeated differencing which maps to `k` here.)

- axis:

  Integer. The axis along which the difference is taken. 2 = along
  rows/down columns (default), 1 = along columns/across rows.

## Value

An Expression representing the kth order differences.
