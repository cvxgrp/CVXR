# Indicator function for constraints

Creates an expression that equals 0 if all constraints are satisfied and
+Inf otherwise. Use this to embed constraints into the objective.

## Usage

``` r
indicator(constraints, err_tol = 0.001)
```

## Arguments

- constraints:

  A list of constraint objects

- err_tol:

  Numeric tolerance for checking constraint satisfaction (default 1e-3)

## Value

An Indicator expression
