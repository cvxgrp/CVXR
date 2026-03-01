# Partial trace of a tensor product expression

Assumes `expr` is a 2D square matrix representing a Kronecker product of
`length(dims)` subsystems. Returns the partial trace over the subsystem
at index `axis` (1-indexed).

## Usage

``` r
partial_trace(expr, dims, axis = 1L)
```

## Arguments

- expr:

  An Expression (2D square matrix)

- dims:

  Integer vector of subsystem dimensions

- axis:

  Integer (1-indexed) subsystem to trace out

## Value

An Expression representing the partial trace
