# Standard Solver Parameter Mappings

Returns a named list mapping standard CVXR parameter names (`reltol`,
`abstol`, `feastol`, `num_iter`) to solver-specific parameter names.
Used internally by
[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md) to translate
standard parameters into solver-native names. A standard parameter with
no entry for the chosen solver does not reach that solver;
[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md) warns when
that happens.

## Usage

``` r
solver_default_param()
```

## Value

A named list keyed by solver name (e.g. `"CLARABEL"`, `"OSQP"`). Each
element is a list of standard parameter mappings, where each mapping has
`name` (the solver-native parameter name the standard name translates
to) and `value` (the solver's own default for that parameter, recorded
for reference when the mapping was written; CVXR does not apply these
values – when you do not set a parameter, the solver uses its own
current default).

## See also

[`psolve`](https://www.cvxgrp.org/CVXR/reference/psolve.md)
