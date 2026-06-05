# Sum of k smallest entries

Sum of k smallest entries

## Usage

``` r
sum_smallest(x, k, axis = NULL, keepdims = FALSE)
```

## Arguments

- x:

  An Expression

- k:

  Number of smallest entries to sum

- axis:

  NULL (all entries), 1 (row-wise), or 2 (column-wise)

- keepdims:

  Logical; keep the reduced dimension as size 1

## Value

An Expression equal to -SumLargest(-x, k)
