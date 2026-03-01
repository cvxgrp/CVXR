# Partial transpose of a tensor product expression

Assumes `expr` is a 2D square matrix representing a Kronecker product of
`length(dims)` subsystems. Returns the partial transpose with the
transpose applied to the subsystem at index `axis` (1-indexed).

## Usage

``` r
partial_transpose(expr, dims, axis = 1L)
```

## Arguments

- expr:

  An Expression (2D square matrix)

- dims:

  Integer vector of subsystem dimensions

- axis:

  Integer (1-indexed) subsystem to transpose

## Value

An Expression representing the partial transpose
