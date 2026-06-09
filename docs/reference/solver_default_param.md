# Standard Solver Parameter Mappings

Returns a named list mapping standard CVXR parameter names (`reltol`,
`abstol`, `feastol`, `num_iter`) to solver-specific parameter names and
their default values. Used internally by
[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md) to translate
standard parameters into solver-native names.

## Usage

``` r
solver_default_param()
```

## Value

A named list keyed by solver name (e.g. `"CLARABEL"`, `"OSQP"`). Each
element is a list of standard parameter mappings, where each mapping has
`name` (solver-native parameter name) and `value` (default value).

## See also

[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)
