# Set the label of an expression

R replacement form of
[`set_label()`](https://www.cvxgrp.org/CVXR/reference/set_label.md).
`label(x) <- value` stores `value` (coerced to character) in `x`'s
internal label slot; setting to `NULL` clears the label. Equivalent to
`x <- set_label(x, value)`.

## Usage

``` r
label(x) <- value
```

## Arguments

- x:

  An Expression object.

- value:

  A character string, or `NULL` to clear.

## Value

`x`, invisibly, with the label updated.

## See also

[`set_label()`](https://www.cvxgrp.org/CVXR/reference/set_label.md),
[`format_labeled()`](https://www.cvxgrp.org/CVXR/reference/format_labeled.md)
