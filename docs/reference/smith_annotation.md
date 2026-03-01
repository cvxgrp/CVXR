# Smith Form Annotation for an Expression Node

Returns LaTeX-math annotation data for visualizing the canonicalization
pipeline. Each atom class can override this to provide a custom LaTeX
name, definition, and conic form. The default stub auto-generates from
class metadata.

## Usage

``` r
smith_annotation(expr, aux_var = "t", child_vars = character(0), ...)
```

## Arguments

- expr:

  An Expression, Atom, or Leaf.

- aux_var:

  Character: the auxiliary variable name assigned to this node (e.g.,
  "t_3").

- child_vars:

  Character vector: auxiliary variable names of the children.

- ...:

  Reserved for future use.

## Value

A list with components: latex_name, latex_definition, conic, doc_topic,
developer.
