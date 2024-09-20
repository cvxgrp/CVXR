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


# Create lists to represent linear operators and constraints
LinOp <- function(type, dim, args = list(), data = NULL, class = "LinOp", id = uuid::UUIDgenerate()) {
  if(is.null(dim)) dim <- c(1,1)   # TODO: Get rid of this with proper multi-dimensional handling.
  if(!is.character(type)) stop("type must be a character string")
  if(!is.numeric(dim)) stop("dim must be a numeric vector")
  if(!is.list(args)) stop("args must be a list of arguments")
  list(type = type, dim = dim, args = args, data = data, id = id, class = "LinOp")
}

LinConstr <- function(expr, constr_id, dim, class = "LinConstr", id = uuid::UUIDgenerate()) {
  if(is.null(dim)) dim <- c(1,1)   # TODO: Get rid of this with proper multi-dimensional handling.
    ##if(!is.character(constr_id)) stop("constr_id must be a character string")
  if(!is.integer(constr_id)) stop("constr_id must be an integer")
  if(!is.numeric(dim)) stop("dim must be a numeric vector")
  list(expr = expr, constr_id = constr_id, dim = dim, id = id, class = class)
}

LinEqConstr <- function(expr, constr_id, dim) { LinConstr(expr, constr_id, dim, class = "LinEqConstr") }
LinLeqConstr <- function(expr, constr_id, dim) { LinConstr(expr, constr_id, dim, class = "LinLeqConstr") }

## get_id <- function() {
##   # sample.int(.Machine$integer.max, 1)
##     uuid::UUIDgenerate()
## }

create_var <- function(dim, var_id = get_id()) {
  LinOp(VARIABLE, dim, list(), var_id)
}

create_param <- function(value, dim) {
  LinOp(PARAM, dim, list(), value)
}

create_const <- function(value, dim, sparse = FALSE) {
  if(all(dim == c(1,1)))
    op_type <- SCALAR_CONST
  else if(sparse)
    op_type <- SPARSE_CONST
  else
    op_type <- DENSE_CONST
  LinOp(op_type, dim, list(), value)
}

lo.is_scalar <- function(operator) {
  length(operator$dim) == 0 || as.integer(prod(operator$dim)) == 1
}

lo.is_const <- function(operator) {
  operator$type %in% c(SCALAR_CONST, SPARSE_CONST, DENSE_CONST)
}

lo.sum_expr <- function(operators) {
  LinOp(SUM, operators[[1]]$dim, operators)
}

lo.neg_expr <- function(operator) {
  LinOp(NEG, operator$dim, list(operator))
}

lo.sub_expr <- function(lh_op, rh_op) {
  lo.sum_expr(list(lh_op, lo.neg_expr(rh_op)))
}

lo.promote_lin_ops_for_mul <- function(lh_op, rh_op) {
  tmp <- mul_dims_promote(lh_op$dim, rh_op$dim)
  lh_dim <- tmp[[1]]
  rh_dim <- tmp[[2]]
  new_dim <- tmp[[3]]
  lh_op <- LinOp(lh_op$type, lh_dim, lh_op$args, lh_op$data)
  rh_op <- LinOp(rh_op$type, rh_dim, rh_op$args, rh_op$data)
  list(lh_op, rh_op, new_dim)
}

lo.mul_expr <- function(lh_op, rh_op, dim) {
  LinOp(MUL_EXPR, dim, list(rh_op), lh_op)
}

lo.rmul_expr <- function(lh_op, rh_op, dim) {
  LinOp(RMUL_EXPR, dim, list(lh_op), rh_op)
}

lo.multiply <- function(lh_op, rh_op) {
  LinOp(MUL_ELEM, lh_op$dim, list(rh_op), lh_op)
}

lo.kron <- function(lh_op, rh_op, dim) {
  LinOp(KRON, dim, list(rh_op), lh_op)
}

lo.div_expr <- function(lh_op, rh_op) {
  LinOp(DIV, lh_op$dim, list(lh_op), rh_op)
}

lo.promote <- function(operator, dim) {
  LinOp(PROMOTE, dim, list(operator))
}

lo.sum_entries <- function(operator, dim) {
  LinOp(SUM_ENTRIES, dim, list(operator))
}

lo.trace <- function(operator) {
  LinOp(TRACE, c(1,1), list(operator))
}

lo.index <- function(operator, dim, keys) {
  LinOp(INDEX, dim, list(operator), keys)
}

lo.conv <- function(lh_op, rh_op, dim) {
  LinOp(CONV, dim, list(rh_op), lh_op)
}

lo.transpose <- function(operator) {
  if(length(operator$dim) < 2)
    operator
  else if(length(operator$dim) > 2)
    stop("Unimplemented")
  else {
    new_dim = c(operator$dim[2], operator$dim[1])
    LinOp(TRANSPOSE, new_dim, list(operator))
  }
}

lo.reshape <- function(operator, dim) {
  LinOp(RESHAPE_EXPR, dim, list(operator))
}

lo.diag_vec <- function(operator) {
  new_dim <- c(operator$dim[1], operator$dim[1])
  LinOp(DIAG_VEC, new_dim, list(operator))
}

lo.diag_mat <- function(operator) {
  new_dim = c(operator$dim[1], 1)
  LinOp(DIAG_MAT, new_dim, list(operator))
}

lo.upper_tri <- function(operator) {
  entries <- operator$dim[1] * operator$dim[2]
  new_dim <- c(floor((entries - operator$dim[1])/2), 1)
  LinOp(UPPER_TRI, new_dim, list(operator))
}

lo.hstack <- function(operators, dim) {
  LinOp(HSTACK, dim, operators)
}

lo.vstack <- function(operators, dim) {
  LinOp(VSTACK, dim, operators)
}

get_constr_expr <- function(lh_op, rh_op) {
  if(missing(rh_op))
    lh_op
  else
    lo.sum_expr(list(lh_op, lo.neg_expr(rh_op)))
}

create_eq <- function(lh_op, rh_op, constr_id = get_id()) {
  expr <- get_constr_expr(lh_op, rh_op)
  LinEqConstr(expr, constr_id, lh_op$dim)
}

create_leq <- function(lh_op, rh_op, constr_id = get_id()) {
  expr <- get_constr_expr(lh_op, rh_op)
  LinLeqConstr(expr, constr_id, lh_op$dim)
}

create_geq <- function(lh_op, rh_op, constr_id = get_id()) {
  if(!missing(rh_op))
    rh_op <- lo.neg_expr(rh_op)
  create_leq(lo.neg_expr(lh_op), rh_op, constr_id)
}

get_expr_vars <- function(operator) {
  if(operator$type == VARIABLE)
    list(list(operator$data, operator$dim))
  else {
    vars_ <- list()
    for(arg in operator$args)
      vars_ <- c(vars_, get_expr_vars(arg))
    vars_
  }
}

get_expr_params <- function(operator) {
  if(operator$type == PARAM)
    parameters(operator$data)
  else {
    params <- list()
    for(arg in operator$args)
      params <- c(params, get_expr_params(arg))
    if(is(operator$data, "LinOp"))
      params <- c(params, get_expr_params(operator$data))
    params
  }
}

copy_constr <- function(constr, func) {
  expr <- func(constr$expr)
  new(class(constr), expr, constr$constr_id, constr$dim)
}

replace_new_vars <- function(expr, id_to_new_var) {
  if(expr$type == VARIABLE && expr$data %in% id_to_new_var)
    id_to_new_var[expr$data]
  else {
    new_args <- list()
    for(arg in expr$args)
      new_args <- c(new_args, replace_new_vars(arg, id_to_new_var))
    LinOp(expr$type, expr$dim, new_args, expr$data)
  }
}

replace_params_with_consts <- function(expr) {
  if(expr$type == PARAM)
    create_const(expr$data$value, expr$dim)
  else {
    new_args <- list()
    for(arg in expr$args)
      new_args <- c(new_args, replace_params_with_consts(arg))
    # Data could also be a parameter
    if(is(expr$data, "LinOp") && expr$data$type == PARAM) {
      data_lin_op <- expr$data
      data <- create_const(data_lin_op$data$value, data_lin_op$dim)
    } else
      data <- expr$data
    LinOp(expr$type, expr$dim, new_args, data)
  }
}
