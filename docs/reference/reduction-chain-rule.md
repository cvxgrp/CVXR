# Reduction chain-rule hooks (dict-in / dict-out)

Walk the solving chain during `Problem$backward()` /
`Problem$derivative()`. Each hook takes and returns a named list keyed
by `as.character(leaf@id)` whose values are shaped arrays (the leaf's
`dim`). The base `Reduction` methods are identity pass-throughs.

## Usage

``` r
var_backward(x, del_vars)

var_forward(x, dvars)

param_backward(x, dparams)

param_forward(x, param_deltas)
```

## Arguments

- x:

  A `Reduction`.

- del_vars:

  Named list `var-id -> gradient array` (outer representation).

- dvars:

  Named list `var-id -> delta array` (inner representation).

- dparams:

  Named list `param-id -> gradient array` (inner representation).

- param_deltas:

  Named list `param-id -> delta array` (outer representation).

## Value

For `var_backward`: the same map in the inner (reduced) representation.

For `var_forward`: the same map in the outer (original) representation.

For `param_backward`: the same map in the outer (original)
representation.

For `param_forward`: the same map in the inner (transformed)
representation.
