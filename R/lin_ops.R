LinOp <- setClass("LinOp", representation(type = "character", size = "numeric", args = "list", data = "ConstValORExpr"))

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

create_const <- function(value, size, sparse = FALSE) {
  if(all(size == c(1,1)))
    op_type <- SCALAR_CONST
  else if(sparse)
    op_type <- SPARSE_CONST
  else
    op_type <- DENSE_CONST
  LinOp(op_type, size, list(), value)
}

create_param <- function(value, size) {
  LinOp(PARAM, size, list(), value)
}

create_var <- function(size, var_id) {
  if(missing(var_id))
    var_id <- 12345    # TODO: Get unique ID
  LinOp(VARIABLE, size, list(), var_id)
}

sum_expr <- function(operators) {
  LinOp(SUM, size(operators[1]), operators, NA_real_)
}

neg_expr <- function(operator) {
  LinOp(NEG, size(operator), list(operator), NA_real_)  
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
  LinOp(MUL_ELEM, size(lh_op), list(rh_op), lh_op)
}

kron <- function(lh_op, rh_op, size) {
  LinOp(KRON, size, list(rh_op), lh_op)
}

div_expr <- function(lh_op, rh_op) {
  LinOp(DIV, size(lh_op), list(lh_op), rh_op)
}

promote <- function(operator, size) {
  LinOp(PROMOTE, size, list(operator), NA_real_)
}

diag_vec <- function(operator) {
  size <- c(size(operator)[1], size(operator)[1])
  LinOp(DIAG_VEC, size, list(operator), NA_real_)
}