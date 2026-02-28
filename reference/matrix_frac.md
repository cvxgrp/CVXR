# Matrix fractional function

Computes \\\mathrm{trace}(X^T P^{-1} X)\\. If P is a constant matrix,
uses a QuadForm shortcut for efficiency.

## Usage

``` r
matrix_frac(X, P)
```

## Arguments

- X:

  A matrix expression (n by m)

- P:

  A square matrix expression (n by n), must be PSD

## Value

An expression representing \\\mathrm{trace}(X^T P^{-1} X)\\
