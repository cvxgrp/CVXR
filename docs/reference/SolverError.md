# Solver Error condition

Creates a custom R condition of class "SolverError" for solver failures.

## Usage

``` r
SolverError(message, call = sys.call(-1L))
```

## Arguments

- message:

  Character error message

- call:

  The call to include in the condition (default: caller's call)

## Value

A condition object of class c("SolverError", "error", "condition")
