# Create a Positive Semidefinite Constraint

Constrains a square matrix expression to be positive semidefinite (PSD):
\\X \succeq 0\\. The expression must be square.

## Usage

``` r
PSD(expr, constr_id = NULL)
```

## Arguments

- expr:

  A CVXR expression representing a square matrix.

- constr_id:

  Optional integer constraint ID.

## Value

A `PSD` constraint object.
