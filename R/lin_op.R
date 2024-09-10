## Analog of cvxpy/lin_ops/lin_op.py + cvxpy/lin_ops/lin_constraints.py

# The types of linear operators.
# A variable.
# Data: var id.
VARIABLE = "variable"
# Some function of parameters.
# Data: CVXPY expression.
PARAM = "param"
# Promoting a scalar expression.
# Data: None
PROMOTE = "promote"
# Multiplying an expression by a constant.
# Data: LinOp evaluating to the left hand multiple.
MUL = "mul"
# Multiplying an expression by a constant on the right.
# Data: LinOp evaluating to the right hand multiple.
RMUL = "rmul"
# Multiplying an expression elementwise by a constant.
# Data: LinOp evaluating to the left hand multiple.
MUL_ELEM = "mul_elem"
# Dividing an expression by a scalar constant.
# Data: LinOp evaluating to the divisor.
DIV = "div"
# Summing expressions.
SUM = "sum"
# Negating an expression.
NEG = "neg"
# An index/slice into an expression.
# Data: (row slice, col slice).
INDEX = "index"
# The transpose of an expression.
# Data: None.
TRANSPOSE = "transpose"
# The sum of the entries of an expression.
# Data: None
SUM_ENTRIES = "sum_entries"
# The sum of the diagonal entries of an expression.
# Data: None
TRACE = "trace"
# An expression cast into a different shape.
# Data: None
RESHAPE = "reshape"
# Converts a vector to a diagonal matrix.
# Data: None
DIAG_VEC = "diag_vec"
# Converts the diagonal of a matrix to a vector.
# Data: None
DIAG_MAT = "diag_mat"
# Vectorized upper triangular portion of a matrix.
# Data: None
UPPER_TRI = "upper_tri"
# The 1D discrete convolution of two vectors.
# Data: LinOp evaluating to the left hand term.
CONV = "conv"
# Horizontally concatenating operators.
# Data: None
HSTACK = "hstack"
# Vertically concatenating operators.
# Data: None
VSTACK = "vstack"
# A scalar constant.
# Data: Python float.
SCALAR_CONST = "scalar_const"
# A dense matrix/vector constant.
# Data: NumPy matrix.
DENSE_CONST = "dense_const"
# A sparse matrix constant.
# Data: SciPy sparse matrix.
SPARSE_CONST = "sparse_const"
# An expression with no variables.
# Data: None
NO_OP = "no_op"
# The Kronecker product of two matrices.
# Data: LinOp evaluating to the left hand term (variable in the right-hand term).
KRON = "kron" ## equivalent to KRON_R, here for backwards compatibility.
# The Kronecker product of two matrices.
# Data: LinOp evaluating to the left hand term (variable in the right-hand term).
KRON_R = "kron_r"
# Data: LinOp evaluating to the right hand term (variable in the left-hand term).
KRON_L = "kron_l"

# ID in coefficients for constants.
CONSTANT_ID = "constant_id"

# Types of linear operators
LINOP_TYPES <-
  c( ## from LinOp.hpp
    VARIABLE,
    PARAM,
    PROMOTE,
    MUL,
    RMUL,
    MUL_ELEM,
    DIV,
    SUM,
    NEG,
    INDEX,
    TRANSPOSE,
    SUM_ENTRIES,
    TRACE,
    RESHAPE,
    DIAG_VEC,
    DIAG_MAT,
    UPPER_TRI,
    CONV,
    HSTACK,
    VSTACK,
    SCALAR_CONST,
    DENSE_CONST,
    SPARSE_CONST,
    NO_OP,
    KRON,
    KRON_R,
    KRON_L
  )
names(LINOP_TYPES) <- toupper(LINOP_TYPES)

#' An R LinOp class as distinguished from a C++ class of the same name. Always has a unique id
#' @param type the type of LinOp, one of the types above
#' @param dim the shape of the LinOp, a tuple, so for us a vector of integers
#' @param args the arguments of the LinOp
#' @param data the data for the LinOp, which is later set to C++ LinOp objects' linOp_data_ field
#' @importFrom uuid UUIDgenerate
#' @return an object of class "LinOp"
LinOp <- function(type, dim, args, data = NULL, id = uuid::UUIDgenerate()) {
  get_type <- function() type
  set_type <- function(what) type <<- what
  get_dim <- function() dim
  set_dim <- function(what) dim <<- what
  get_args <- function() args
  set_args <- function(what) args <<- what
  result <- list(get_type = get_type, set_type = set_type,
                 get_dim = get_dim, set_dim = set_dim,
                 get_args = get_args, set_args = set_args,
                 get_id = function() id)
  result$self <- result
  class(result) <- "LinOp"
  result
}

#' @method print LinOp
print.LinOp <- function(x, ...) {
  sprintf("LinOp(%s, dim = [%s])",
          x$self$get_type(),
          paste0(x$self$get_dim(), collapse = ", ")
          )
}

#' Make a (R) Linear Constraint
#' @param expr the expression
#' @param constr_id the constaint id
#' @param dim the shape

#' @param data the data for the LinOp, which is later set to C++ LinOp objects' linOp_data_ field
#' @return an object of class "LinOp"
make_lin_constraint <- function(expr, constr_id, dim, class = "LinConstr") {
  result <- list(get_expr = function() expr
                 get_constr_id = function() constr_id,
                 get_dim = function() dim)
  result$self <- result
  class(result) <- class
  result
}

LinConstr <- make_lin_constraint
LinEqConstr <- function(expr, constr_id, dim) { LinConstr(expr, constr_id, dim, class = "LinEqConstr") }
LinLeqConstr <- function(expr, constr_id, dim) { LinConstr(expr, constr_id, dim, class = "LinLeqConstr") }

print_lin_constraint <- function(x, ...) {
  sprintf("%s(%s, dim = [%s])",
          class(x),
          x$self$get_expr(),
          paste0(x$self$get_dim(), collapse = ", ")
          )
}

#' @method print LinConstr
print.LinConstr <- print_lin_constraint

#' @method print LinEqConstr
print.LinEqConstr <- print_lin_constraint

#' @method print LinLeqConstr
print.LinLeqConstr <- print_lin_constraint






