# Geometric matrix multiplication A diamond X

Computes the geometric matrix product where (A diamond X)\_ij = prod_k
X_kj^A_ik. Log-log affine atom for DGP. Solve with
`psolve(problem, gp = TRUE)`.

## Usage

``` r
gmatmul(A, X)
```

## Arguments

- A:

  A constant matrix

- X:

  An Expression (positive matrix)

## Value

A Gmatmul atom

## Examples

``` r
x <- Variable(2, pos = TRUE)
A <- matrix(c(1, 0, 0, 1), 2, 2)
prob <- Problem(Minimize(sum(gmatmul(A, x))), list(x >= 0.5))
if (FALSE) psolve(prob, gp = TRUE) # \dontrun{}
```
