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

#' An R LinOp class as distinguished from a C++ class of the same name
#' @param type the type of LinOp, one of the types above
#' @param dim the shape of the LinOp, a tuple, so for us a vector of integers
#' @param args the arguments of the LinOp
#' @param data the data for the LinOp, which is later set to C++ LinOp objects' linOp_data_ field
#' @importFrom uuid UUIDgenerate
#' @return an object of class "LinOp"
LinOp <- function(type, dim, args, data = NULL, uuid = uuid::UUIDgenerate()) {
  ## if (is.na(match(type, LINOP_TYPES))) stop(sprintf("LinOp: Unknown type: %s", type))
  self <- environment()
  class(self) <- c("LinOp", class(self))
  self$self <- self  
  ## These checks were in the earlier version. TODO: Need some validation
  ## if(is.null(dim)) dim <- c(1,1)   # TODO: Get rid of this with proper multi-dimensional handling.
  ## if(!is.character(type)) stop("type must be a character string")
  ## if(!is.numeric(dim)) stop("dim must be a numeric vector")
  ## if(!is.list(args)) stop("args must be a list of arguments")
  ## list(type = type, dim = dim, args = args, data = data, class = "LinOp")
  
  
  ## Read THIS NOW!!!
  ## TODO: Not sure if C++ object really needs to be created yet.
  ## Which means the methods below that iterface to the C++ object may not be needed.
  ## That would mean a much simpler interface...
  self
}

#' @method print LinOp
print.LinOp <- function(x, ...) {
  comps <- ls(x)
  methods <- comps[sapply(comps, function(name) is.function(x[[name]]))]
  others <- setdiff(comps, methods)
  
  out <- sprintf("LinOp(%s)", x$type)
  if (length(others) > 0) {
    out <- c(
      out,
      " " = paste0("Slots: ", paste(others, collapse = ", "))
    )
  }
  if (length(methods) > 0) {
    out <- c(
      out,
      " " = paste0("Methods: ", paste(comps[methods], collapse = ", "))
    )
  }
  cli::cli_bullets(out)
}

#' Make a (R) Linear Constraint
#' @param expr the expression
#' @param constr_id the constaint id
#' @param dim the shape
#' @param class the class to set this object to
#' @param data the data for the LinOp, which is later set to C++ LinOp objects' linOp_data_ field
#' @return an object of class "LinOp"
make_lin_constraint <- function(expr, constr_id, dim, class = "LinConstr") {
  self <- environment()
  if(is.null(dim)) dim <- rep(1L, 2L)   # TODO: Get rid of this with proper multi-dimensional handling.
  if(!is.integer(dim)) stop("dim must be a integer vector")
  class(self) <- c(class, class(self))
  rm("class", envir = self)
  self$self <- self  
  self
}

LinConstr <- make_lin_constraint
LinEqConstr <- function(expr, constr_id, dim) { LinConstr(expr, constr_id, dim, class = "LinEqConstr") }
LinLeqConstr <- function(expr, constr_id, dim) { LinConstr(expr, constr_id, dim, class = "LinLeqConstr") }

print_lin_constraint <- function(x, ...) {
  comps <- ls(x)
  methods <- comps[sapply(comps, function(name) is.function(x[[name]]))]
  others <- setdiff(comps, methods)
  
  out <- sprintf("%s", class(x)[[1L]])
  if (length(others) > 0) {
    out <- c(
      out,
      " " = paste0("Slots: ", paste(others, collapse = ", "))
    )
  }
  if (length(methods) > 0) {
    out <- c(
      out,
      " " = paste0("Methods: ", paste(comps[methods], collapse = ", "))
    )
  }
  cli::cli_bullets(out)
}

#' @method print LinConstr
print.LinConstr <- print_lin_constraint

#' @method print LinEqConstr
print.LinEqConstr <- print_lin_constraint

#' @method print LinLeqConstr
print.LinLeqConstr <- print_lin_constraint






