# Sum the entries of an expression

Sum the entries of an expression

## Usage

``` r
sum_entries(x, axis = NULL, keepdims = FALSE)
```

## Arguments

- x:

  An Expression or numeric value.

- axis:

  NULL (sum all), 1 (row-wise, like apply(X,1,sum)), or 2 (column-wise,
  like apply(X,2,sum)).

- keepdims:

  Logical: if TRUE, keep the reduced dimension as size 1.

## Value

A SumEntries expression.
