# Convert CVXR Object to LaTeX

Renders a CVXR `Problem`, `Expression`, or `Constraint` as a LaTeX
string. Problem-level output uses the `optidef` package (`mini*`/`maxi*`
environments) and atom macros from `dcp.sty` (shipped as
`system.file("sty", "dcp.sty", package = "CVXR")`).

## Usage

``` r
to_latex(x, ...)
```

## Arguments

- x:

  A `Problem`, `Expression`, `Constraint`, or `Objective`.

- ...:

  Reserved for future options.

## Value

A character string containing LaTeX code.

## Examples

``` r
x <- Variable(3, name = "x")
cat(to_latex(p_norm(x, 2)))
#> \cvxnorm{x}_{2}
# \cvxnorm{x}_2
```
