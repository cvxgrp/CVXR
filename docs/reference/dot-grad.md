# Per-atom Subgradient Hook (private)

R counterpart of CVXPY's `Atom._grad(self, values)`. Returns a list of
one Jacobian per argument, shape `(prod(arg.shape), prod(self.shape))`.
Called by the chain-rule walker `grad(x)`. The leading dot follows R's
private-name convention (cf. `.Machine`, `.libPaths`, `.cvxr_*`); the
name maps one-to-one onto CVXPY's `_grad`.

## Usage

``` r
.grad(x, values, ...)
```

## Arguments

- x:

  An atom expression.

- values:

  A list of numeric values, one per argument.

- ...:

  Reserved for future use; method dispatch ignores it.

## Value

A list of sparse Jacobians (one per argument).
