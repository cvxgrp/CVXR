# Reduction chain-rule hooks (backward / forward, parameter / variable)

Reduction chain-rule hooks (backward / forward, parameter / variable)

## Usage

``` r
param_backward(x, param, dparams)

param_forward(x, param, delta)

var_backward(x, var, value)

var_forward(x, var, value)
```

## Arguments

- x:

  A `Reduction`.

- param:

  A `Parameter` from the original problem.

- dparams:

  Named list of param-id -\> gradient arrays for the reduced
  (transformed) parameters.

- delta:

  Numeric array of perturbations to the original `param`.

- var:

  A `Variable` from the original problem.

- value:

  Numeric array (gradient or delta) attached to `var`.

## Value

A gradient array for `param`, or `NULL` if this reduction does not touch
`param`.

Named list of transformed-param-id -\> delta arrays, or `NULL` if this
reduction does not touch `param`.
