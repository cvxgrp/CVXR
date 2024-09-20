# Notes on code.

## LinOp.R

CVXcanon.LinOp: R object backing CPP object

## LinOpVector.R

CVXcanon.LinOp: R object backing CPP LinOpVector object

## lin_op.R

Defines: Singleton named `LINOP_TYPES` character array matching C++ LinOp types
implementations

`LinOp`: LinOp (pure R) factory method
`LinConstr`: LinConstraint (pure R) factory method
