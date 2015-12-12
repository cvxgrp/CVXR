# Types of linear operators
VARIABLE = "variable"
PROMOTE = "promote"
MUL = "mul"
RMUL = "rmul"
MUL_ELEM = "mul_elem"
DIV = "div"
SUM = "sum"
NEG = "neg"
INDEX = "index"
TRANSPOSE = "transpose"
SUM_ENTRIES = "sum_entries"
TRACE = "trace"
RESHAPE = "reshape"
DIAG_VEC = "diag_vec"
DIAG_MAT = "diag_mat"
UPPER_TRI = "upper_tri"
CONV = "conv"
KRON = "kron"
HSTACK = "hstack"
VSTACK = "vstack"
SCALAR_CONST = "scalar_const"
DENSE_CONST = "dense_const"
SPARSE_CONST = "sparse_const"
PARAM = "param"
NO_OP = "no_op"
CONSTANT_ID = "constant_id"

# Create lists to represent linear operators and constraints
LinOp <- function(type, size, args = list(), data = NA_real_) {
  if(!is.character(type)) stop("type must be a character string")
  if(!is.numeric(size)) stop("size must be a numeric vector")
  if(!is.list(args)) stop("args must be a list of arguments")
  list(type = type, size = size, args = args, data = data)
}

LinConstr <- function(expr, constr_id, size) {
  if(!is.character(constr_id)) stop("constr_id must be a character string")
  if(!is.numeric(size)) stop("size must be a numeric vector")
  list(expr = expr, constr_id = constr_id, size = size)
}

LinEqConstr <- function(expr, constr_id, size) { LinConstr(expr, constr_id, size) }
LinLeqConstr <- function(expr, constr_id, size) { LinConstr(expr, constr_id, size) }

get_id <- function() {
  sample.int(.Machine$integer.max, 1)
}

create_var <- function(size, var_id = get_id()) {
  LinOp(VARIABLE, size, list(), var_id)
}

create_param <- function(value, size) {
  LinOp(PARAM, size, list(), value)
}

create_const <- function(value, size, sparse = FALSE) {
  if(all(size == c(1,1)))
    op_type <- SCALAR_CONST
  else if(sparse)
    op_type <- SPARSE_CONST
  else
    op_type <- DENSE_CONST
  LinOp(op_type, size, list(), value)
}

sum_expr <- function(operators) {
  LinOp(SUM, operators[[1]]$size, operators)
}

neg_expr <- function(operator) {
  LinOp(NEG, operator$size, list(operator))  
}

sub_expr <- function(lh_op, rh_op) {
  sum_expr(list(lh_op, neg_expr(rh_op)))
}

mul_expr <- function(lh_op, rh_op, size) {
  LinOp(MUL, size, list(rh_op), lh_op)
}

rmul_expr <- function(lh_op, rh_op, size) {
  LinOp(RMUL, size, list(lh_op), rh_op)
}

mul_elemwise <- function(lh_op, rh_op) {
  LinOp(MUL_ELEM, lh_op$size, list(rh_op), lh_op)
}

kron <- function(lh_op, rh_op, size) {
  LinOp(KRON, size, list(rh_op), lh_op)
}

div_expr <- function(lh_op, rh_op) {
  LinOp(DIV, lh_op$size, list(lh_op), rh_op)
}

promote <- function(operator, size) {
  LinOp(PROMOTE, size, list(operator))
}

sum_entries <- function(operator) {
  LinOp(SUM_ENTRIES, c(1,1), list(operator))
}

trace <- function(operator) {
  LinOp(TRACE, c(1,1), list(operator))
}

index <- function(operator, size, keys) {
  LinOp(INDEX, size, list(operator), keys)
}

conv <- function(lh_op, rh_op, size) {
  LinOp(CONV, size, list(rh_op), lh_op)
}

transpose <- function(operator) {
  size = c(operator$size[2], operator$size[1])
  LinOp(TRANSPOSE, size, list(operator))
}

reshape <- function(operator, size) {
  LinOp(RESHAPE, size, list(operator))
}

diag_vec <- function(operator) {
  size <- c(operator$size[1], operator$size[1])
  LinOp(DIAG_VEC, size, list(operator))
}

diag_mat <- function(operator) {
  size = c(operator$size[1], 1)
  LinOp(DIAG_MAT, size, list(operator))
}

upper_tri <- function(operator) {
  entries <- operator$size[1] * operator$size[2]
  size <- c(floor((entries - operator$size[1])/2), 1)
  LinOp(UPPER_TRI, size, list(operator))
}

hstack <- function(operators, size) {
  LinOp(HSTACK, size, operators)
}

vstack <- function(operators, size) {
  LinOp(VSTACK, size, operators)
}

get_constr_expr <- function(lh_op, rh_op) {
  if(missing(rh_op))
    lh_op
  else
    sum_expr(list(lh_op, neg_expr(rh_op)))
}

create_eq <- function(lh_op, rh_op, constr_id = get_id()) {
  expr <- get_constr_expr(lh_op, rh_op)
  LinEqConstr(expr, constr_id, lh_op$size)
}

create_leq <- function(lh_op, rh_op, constr_id = get_id()) {
  expr <- get_constr_expr(lh_op, rh_op)
  LinLeqConstr(expr, constr_id, lh_op$size)
}

create_geq <- function(lh_op, rh_op, constr_id) {
  if(!missing(rh_op))
    rh_op <- neg_expr(rh_op)
  create_leq(neg_expr(lh_op), rh_op, constr_id)
}

get_expr_vars <- function(operator) {
  if(operator$type == VARIABLE)
    list(list(operator$data, operator$size))
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
  new(class(constr), expr, constr$constr_id, constr$size)
}

replace_new_vars <- function(expr, id_to_new_var) {
  if(expr$type == VARIABLE && expr$data %in% id_to_new_var)
    id_to_new_var[expr$data]
  else {
    new_args <- list()
    for(arg in expr$args)
      new_args <- c(new_args, replace_new_vars(arg, id_to_new_var))
    LinOp(expr$type, expr$size, new_args, expr$data)
  }
}

replace_params_with_consts <- function(expr) {
  if(expr$type == PARAM)
    create_const(expr$data$value, expr$size)
  else {
    new_args <- list()
    for(arg in expr$args)
      new_args <- c(new_args, replace_params_with_consts(arg))
    # Data could also be a parameter
    if(is(expr$data, "LinOp") && expr$data$type == PARAM) {
      data_lin_op <- expr$data
      data <- create_const(data_lin_op$data$value, data_lin_op$size)
    } else
      data <- expr$data
    LinOp(expr$type, expr$size, new_args, data)
  }
}
