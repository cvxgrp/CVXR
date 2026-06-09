# Sign of a product of two expressions

Determines whether the product of two expressions is nonnegative,
nonpositive, or unknown, using the sign multiplication table.

## Usage

``` r
mul_sign(lh_expr, rh_expr)
```

## Arguments

- lh_expr:

  An Expression object (left-hand operand)

- rh_expr:

  An Expression object (right-hand operand)

## Value

Named logical vector c(is_nonneg, is_nonpos)
