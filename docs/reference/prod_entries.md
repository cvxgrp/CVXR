# Product of entries along an axis

Used in DGP (geometric programming) context. Solve with
`psolve(problem, gp = TRUE)`.

## Usage

``` r
prod_entries(x, axis = NULL, keepdims = FALSE)
```

## Arguments

- x:

  An Expression

- axis:

  NULL (all), 1 (row-wise), or 2 (column-wise)

- keepdims:

  Whether to keep reduced dimensions

## Value

A Prod atom

## Examples

``` r
x <- Variable(3, pos = TRUE)
prob <- Problem(Minimize(prod_entries(x)), list(x >= 2))
if (FALSE) psolve(prob, gp = TRUE) # \dontrun{}
```
