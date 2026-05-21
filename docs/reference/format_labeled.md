# Pretty-print an expression with labels substituted

Recursive analogue of
[`expr_name()`](https://www.cvxgrp.org/CVXR/reference/expr_name.md) that
substitutes user-supplied labels (see
[`set_label()`](https://www.cvxgrp.org/CVXR/reference/set_label.md)) for
sub-expressions wherever they are set, falling back to the structural
name on unlabelled nodes. Mirrors CVXPY's `Expression.format_labeled`.

## Usage

``` r
format_labeled(x)
```

## Arguments

- x:

  An Expression object.

## Value

A character string.

## See also

[`set_label()`](https://www.cvxgrp.org/CVXR/reference/set_label.md),
[`label()`](https://www.cvxgrp.org/CVXR/reference/label.md)
