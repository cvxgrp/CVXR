# Quadratic form x^T P x

When `x` is constant, returns `t(Conj(x)) %*% P %*% x` (affine in `P`).
When `P` is constant, returns a `QuadForm` atom (quadratic in `x`). At
least one of `x` or `P` must be constant.

## Usage

``` r
quad_form(x, P, assume_PSD = FALSE)
```

## Arguments

- x:

  An Expression (vector)

- P:

  An Expression (square matrix, symmetric/Hermitian)

- assume_PSD:

  If TRUE, assume P is PSD without checking (only when P is constant).

## Value

A QuadForm atom or an affine Expression
