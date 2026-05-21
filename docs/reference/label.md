# Get the label of an expression

Returns the human-readable label set via
[`set_label()`](https://www.cvxgrp.org/CVXR/reference/set_label.md) (or
`label(x) <- ...`), or `NULL` if no label has been set.

## Usage

``` r
label(x)
```

## Arguments

- x:

  An Expression object.

## Value

A length-1 character string, or `NULL`.

## See also

[`set_label()`](https://www.cvxgrp.org/CVXR/reference/set_label.md),
[`format_labeled()`](https://www.cvxgrp.org/CVXR/reference/format_labeled.md)
