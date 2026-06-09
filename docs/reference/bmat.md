# Construct a Block Matrix

Takes a list of lists. Each internal list is stacked horizontally. The
internal lists are stacked vertically.

## Usage

``` r
bmat(block_lists)
```

## Arguments

- block_lists:

  A list of lists of Expression objects (or numerics). Each inner list
  forms one block row.

## Value

An Expression representing the block matrix.
