# Unity resolvent (I - X) inverse for positive square matrix X

Log-log convex atom for DGP. Solve with `psolve(problem, gp = TRUE)`.

## Usage

``` r
eye_minus_inv(X)
```

## Arguments

- X:

  An Expression (positive square matrix with spectral radius \< 1)

## Value

An EyeMinusInv atom

## Examples

``` r
X <- Variable(c(2, 2), pos = TRUE)
prob <- Problem(Minimize(sum(eye_minus_inv(X))), list(X <= 0.4))
if (FALSE) psolve(prob, gp = TRUE, solver = "SCS") # \dontrun{}
```
