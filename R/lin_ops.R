# Types of linear operators
VARIABLE <- "variable"
PARAM <- "param"
PROMOTE <- "promote"
MUL_EXPR <- "mul_expr"
RMUL_EXPR <- "rmul_expr"
MUL_ELEM <- "mul_elem"
DIV <- "div"
SUM <- "sum"
NEG <- "neg"
INDEX <- "index"
TRANSPOSE <- "transpose"
SUM_ENTRIES <- "sum_entries"
TRACE <- "trace"
RESHAPE_EXPR <- "reshape_expr"
DIAG_VEC <- "diag_vec"
DIAG_MAT <- "diag_mat"
UPPER_TRI <- "upper_tri"
CONV <- "conv"
HSTACK <- "hstack"
VSTACK <- "vstack"
SCALAR_CONST <- "scalar_const"
DENSE_CONST <- "dense_const"
SPARSE_CONST <- "sparse_const"
NO_OP <- "no_op"
KRON <- "kron_r"
KRON_R <- "kron_r"
KRON_L <- "kron_l"

CONSTANT_ID <- "constant_id"

LINOP_TYPES <- c(VARIABLE = "VARIABLE",
                 PARAM = "PARAM",
                 PROMOTE = "PROMOTE",
                 MUL_EXPR = "MUL_EXPR",
                 RMUL_EXPR = "RMUL_EXPR",
                 MUL_ELEM = "MUL_ELEM",
                 DIV = "DIV",
                 SUM = "SUM",
                 NEG = "NEG",
                 INDEX = "INDEX",
                 TRANSPOSE = "TRANSPOSE",
                 SUM_ENTRIES = "SUM_ENTRIES",
                 TRACE = "TRACE",
                 RESHAPE_EXPR = "RESHAPE_EXPR",
                 DIAG_VEC = "DIAG_VEC",
                 DIAG_MAT = "DIAG_MAT",
                 UPPER_TRI = "UPPER_TRI",
                 CONV = "CONV",
                 KRON = "KRON",
                 HSTACK = "HSTACK",
                 VSTACK = "VSTACK",
                 SCALAR_CONST = "SCALAR_CONST",
                 DENSE_CONST = "DENSE_CONST",
                 SPARSE_CONST = "SPARSE_CONST",

                 NO_OP = "NO_OP",
                 CONSTANT_ID = "CONSTANT_ID")

# Create lists to represent linear operators and constraints
LinOp <- function(type, dim, args = list(), data = NULL, class = "LinOp") {
  if(is.null(dim)) dim <- c(1,1)   # TODO: Get rid of this with proper multi-dimensional handling.
  if(!is.character(type)) stop("type must be a character string")
  if(!is.numeric(dim)) stop("dim must be a numeric vector")
  if(!is.list(args)) stop("args must be a list of arguments")
  list(type = type, dim = dim, args = args, data = data, class = "LinOp")
}

LinConstr <- function(expr, constr_id, dim, class = "LinConstr") {
  if(is.null(dim)) dim <- c(1,1)   # TODO: Get rid of this with proper multi-dimensional handling.
    ##if(!is.character(constr_id)) stop("constr_id must be a character string")
  if(!is.integer(constr_id)) stop("constr_id must be an integer")
  if(!is.numeric(dim)) stop("dim must be a numeric vector")
  list(expr = expr, constr_id = constr_id, dim = dim, class = class)
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
