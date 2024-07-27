## CVXPY SOURCE: cvxpy/lin_ops/lin_op.py

# Create lists to represent linear operators and constraints
LinOp <- function(type, dim, args = list(), data = NULL) {
  if (is.null(dim)) dim <- c(1L ,1L )   # TODO: Get rid of this with proper multi-dimensional handling.
  if (!is.character(type)) stop("type must be a character string")
  if (!is.numeric(dim)) stop("dim must be a numeric vector")
  if (!is.list(args)) stop("args must be a list of arguments")
  list(type = type, dim = dim, args = args, data = data, class = "LinOp")
}

## Since lin_op often gets imported in python as lo, we will use that naming convention
# A variable.
# Data: var id.
lo.VARIABLE <- "variable"
# Promoting a scalar expression.
# Data: None
lo.PROMOTE <- "promote"
# Multiplying an expression by a constant.
# Data: LinOp evaluating to the left hand multiple.
lo.MUL <- "mul"
# Multiplying an expression by a constant on the right.
# Data: LinOp evaluating to the right hand multiple.
lo.RMUL <- "rmul"
# Multiplying an expression elementwise by a constant.
# Data: LinOp evaluating to the left hand multiple.
lo.MUL_ELEM <- "mul_elem"
# Dividing an expression by a scalar constant.
# Data: LinOp evaluating to the divisor.
lo.DIV <- "div"
# Summing expressions.
lo.SUM <- "sum"
# Negating an expression.
lo.NEG <- "neg"
# An index/slice into an expression.
# Data: (row slice, col slice).
lo.INDEX <- "index"
# The transpose of an expression.
# Data: None.
lo.TRANSPOSE <- "transpose"
# The sum of the entries of an expression.
# Data: None
lo.SUM_ENTRIES <- "sum_entries"
# The sum of the diagonal entries of an expression.
# Data: None
lo.TRACE <- "trace"
# An expression cast into a different shape.
# Data: None
lo.RESHAPE <- "reshape"
# Converts a vector to a diagonal matrix.
# Data: None
lo.DIAG_VEC <- "diag_vec"
# Converts the diagonal of a matrix to a vector.
# Data: None
lo.DIAG_MAT <- "diag_mat"
# Vectorized upper triangular portion of a matrix.
# Data: None
lo.UPPER_TRI <- "upper_tri"
# The 1D discrete convolution of two vectors.
# Data: LinOp evaluating to the left hand term.
lo.CONV <- "conv"
# The Kronecker product of two matrices.
# Data: LinOp evaluating to the left hand term (variable in the right-hand term).
lo.KRON_R <- "kron_r"
# Data: LinOp evaluating to the right hand term (variable in the left-hand term).
lo.KRON_L <- "kron_l"
# Horizontally concatenating operators.
# Data: None
lo.HSTACK <- "hstack"
# Vertically concatenating operators.
# Data: None
lo.VSTACK <- "vstack"
# A scalar constant.
# Data: Python float.
lo.SCALAR_CONST <- "scalar_const"
# A dense matrix/vector constant.
# Data: NumPy matrix.
lo.DENSE_CONST <- "dense_const"
# A sparse matrix constant.
# Data: SciPy sparse matrix.
lo.SPARSE_CONST <- "sparse_const"
# Some function of parameters.
# Data: CVXPY expression.
lo.PARAM <- "param"
# An expression with no variables.
# Data: None
lo.NO_OP <- "no_op"
# ID in coefficients for constants.
lo.CONSTANT_ID <- "constant_id"

