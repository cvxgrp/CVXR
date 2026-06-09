# Attach a label to an expression

CVXPY-parity setter that returns its first argument so calls can be
chained (e.g. `sum_squares(x) |> set_label("cost")`). See
[`format_labeled()`](https://www.cvxgrp.org/CVXR/reference/format_labeled.md)
for the pretty-printer that consumes labels.

## Usage

``` r
set_label(x, value)
```

## Arguments

- x:

  An Expression object.

- value:

  A label (character; coerced via `as.character`). Pass `NULL` to clear
  an existing label.

## Value

`x` with the label updated.

## See also

[`label()`](https://www.cvxgrp.org/CVXR/reference/label.md),
[`format_labeled()`](https://www.cvxgrp.org/CVXR/reference/format_labeled.md)
