# Check DPP Compliance

Determines whether an expression or problem satisfies the rules of
Disciplined Parameterized Programming (DPP). A DPP-compliant problem
enables caching the compilation across parameter value changes.

## Usage

``` r
is_dpp(x, context = "dcp")
```

## Arguments

- x:

  An expression, constraint, or problem object.

- context:

  Either `"dcp"` (default) or `"dgp"`: which discipline to check
  parameterization against. Mirrors CVXPY's `is_dpp(context=...)`.

## Value

Logical scalar.
