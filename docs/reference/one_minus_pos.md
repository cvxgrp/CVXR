# The difference 1 - x with domain (0, 1)

Log-log concave atom for DGP. Solve with `psolve(problem, gp = TRUE)`.

## Usage

``` r
one_minus_pos(x)
```

## Arguments

- x:

  An Expression (elementwise in (0, 1))

## Value

A OneMInusPos atom

## Examples

``` r
x <- Variable(pos = TRUE)
prob <- Problem(Maximize(one_minus_pos(x)), list(x >= 0.1, x <= 0.5))
if (FALSE) psolve(prob, gp = TRUE) # \dontrun{}
```
