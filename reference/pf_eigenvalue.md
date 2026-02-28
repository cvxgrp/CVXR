# Perron-Frobenius eigenvalue of a positive matrix

Log-log convex atom for DGP. Solve with `psolve(problem, gp = TRUE)`.

## Usage

``` r
pf_eigenvalue(X)
```

## Arguments

- X:

  An Expression (positive square matrix)

## Value

A PfEigenvalue atom (scalar)

## Examples

``` r
X <- Variable(c(2, 2), pos = TRUE)
prob <- Problem(Minimize(pf_eigenvalue(X)),
                list(X[1,1] >= 0.1, X[2,2] >= 0.1))
if (FALSE) psolve(prob, gp = TRUE) # \dontrun{}
```
