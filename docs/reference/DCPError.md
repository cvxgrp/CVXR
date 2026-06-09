# DCP Error condition

Creates a custom R condition of class "DCPError" for disciplined convex
programming violations.

## Usage

``` r
DCPError(message, call = sys.call(-1L))
```

## Arguments

- message:

  Character error message

- call:

  The call to include in the condition (default: caller's call)

## Value

A condition object of class c("DCPError", "error", "condition")
