## Analog of cvxpy/lin_ops/lin_utils.py
## NEEDS: lin_op.R

#' Counter class for ids
#' @param type the type of LinOp, one of the types above
#' @param dim the dim of the LinOp, a tuple, so for us a vector of integers
#' @param args the arguments of the LinOp
#' @param data the data for the LinOp, which is later set to C++ LinOp objects' linOp_data_ field
#' @return an object of class "LinOp"
Counter <- function() {
  self <- environment()
  self$count <- 1L
  class(self) <- c("Counter", class(self))
  self$self <- self
  self
}

ID_COUNTER <- Counter()

#' Returns a new id and updates the id counter
#' @return a new id
get_id <- function() {
  new_id <- ID_COUNTER$count
  ID_COUNTER$count <- ID_COUNTER$count + 1L
  new_id
}

#' Create a new internal variable
#' @param dim a tuple, the (rows, cols) dimensions of the variable
#' @param var_id the id of the variable
#' @return a LinOp representing the new variable
lo.create_var <- function(dim, var_id = get_id()) {
  LinOp(VARIABLE, dim, list(), var_id)
}

#' Wrap a parameter.  (This is the v1.3.0 version and the 1.5.2 version is different!)
#' @param dim a tuple, the (rows, cols) dimensions of the variable
#' @param param_id the id of the variable
#' @return a LinOp wrapping the parameter
lo.create_param <- function(dim, param_id = get_id()) {
  LinOp(PARAM, dim, list(), param_id)
}

#' Wrap a constant
#' @param value a scalar, matrix or sparse matrix to wrap
#' @param dim a tuple, the (rows, cols) dimensions of the variable
#' @param sparse a boolean flag indicating if sparse matrix or not
#' @return a LinOp wrapping the constant
lo.create_const <- function(value, dim, sparse = FALSE) {
  if(all(dim == 1L)) {
    op_type <- SCALAR_CONST
    if (length(dim) > 1L) value <- value[1L, 1L]
  } else if(sparse) {
    op_type <- SPARSE_CONST
  } else {
    op_type <- DENSE_CONST
  }
  LinOp(op_type, dim, list(), value)
}

#' Is the LinOp a scalar?
#' @param operator a LinOp
#' @return TRUE or FALSE
lo.is_scalar <- function(operator) {
  length(operator$dim) == 0L || prod(operator$dim) == 1
}

#' Is the LinOp a constant?
#' @param operator a LinOp
#' @return TRUE or FALSE
lo.is_const <- function(operator) {
  operator$type %in% c(SCALAR_CONST, SPARSE_CONST, DENSE_CONST)
}

#' Sum a list of LinOps
#' @param operators a list of LinOps
#' @return a LinOp representing the sum
lo.sum_expr <- function(operators) {
  LinOp(SUM, operators[[1]]$dim, operators, NULL)
}

#' Negate a LinOp
#' @param operator a LinOp
#' @return the negated LinOp
lo.neg_expr <- function(operator) {
  LinOp(NEG, operator$dim, list(operator), NULL)
}

#' Subtract one LinOp from another
#' @param lh_op a LinOp
#' @param rh_op a LinOp
#' @return a LinOp representing the subtracted result 
lo.sub_expr <- function(lh_op, rh_op) {
  lo.sum_expr(list(lh_op, lo.neg_expr(rh_op)))
}

#' Promote LinOp arguments for multiplication
#' @param lh_op a LinOp
#' @param rh_op a LinOp
#' @return a list of promoted LinOps respectively and the new dimension
lo.promote_lin_ops_for_mul <- function(lh_op, rh_op) {
  promotion <- mul_dims_promote(lh_op$dim, rh_op$dim)
  lh_op <- LinOp(lh_op$type, promotion[[1L]], lh_op$args, lh_op$data)
  rh_op <- LinOp(rh_op$type, promotion[[2L]], rh_op$args, rh_op$data)
  list(lh_op, rh_op, promotion[[3L]])
}


#' Multiply two LinOps, with constant on left
#' @param lh_op a LinOp
#' @param rh_op a LinOp
#' @param dim the dim of the result
#' @return a LinOp representing the product
lo.mul_expr <- function(lh_op, rh_op, dim) {
  LinOp(MUL_EXPR, dim, list(rh_op), lh_op)
}

#' Multiply two LinOps, with constant on right
#' @param lh_op a LinOp
#' @param rh_op a LinOp
#' @param dim the dim of the result
#' @return a LinOp representing the product
lo.rmul_expr <- function(lh_op, rh_op, dim) {
  LinOp(RMUL_EXPR, dim, list(lh_op), rh_op)
}

#' Multiply two LinOps elementwise
#' @param lh_op a LinOp
#' @param rh_op a LinOp
#' @return a LinOp representing the product
lo.multiply <- function(lh_op, rh_op) {
  dim <- pmax(lh_op$dim, rh_op$dim)
  LinOp(MUL_ELEM, dim, list(rh_op), lh_op)
}

#' Construct Kronecker product of two matrices, where the right operand is a Variable
#' @param lh_op a LinOp
#' @param rh_op a LinOp
#' @param dim the dimension of the result
#' @return a LinOp representing the kronecker product
lo.kron_r<- function(lh_op, rh_op, dim) {
  LinOp(KRON_R, dim, list(rh_op), lh_op)
}

#' Construct Kronecker product of two matrices, where the left operand is a Variable
#' @param lh_op a LinOp
#' @param rh_op a LinOp
#' @param dim the dimension of the result
#' @return a LinOp representing the kronecker product
lo.kron_l<- function(lh_op, rh_op, dim) {
  LinOp(KRON_L, dim, list(lh_op), rh_op)
}

#' Divide one LinOp by another
#' @param lh_op a LinOp
#' @param rh_op a LinOp
#' @return a LinOp representing the product
lo.div_expr <- function(lh_op, rh_op) {
  LinOp(DIV, lh_op$dim, list(lh_op), rh_op)
}

#' Promote a scalar operator to a specified shape
#' @param operator a LinOp
#' @param dim the desired shape
#' @return a LinOp representing the promotion
lo.promote <- function(operator, dim) {
  LinOp(PROMOTE, dim, list(operator), NULL)
}

#' Sum the entries of an operator
#' @param operator a LinOp
#' @param dim the dimensions of the sum
#' @return a LinOp representing the sum
lo.sum_entries <- function(operator, dim) {
  LinOp(SUM_ENTRIES, dim, list(operator), NULL)
}

#' Sum the diagonal entries of an operator
#' @param operator a LinOp
#' @return a LinOp representing the sum
lo.trace <- function(operator) {
  LinOp(TRACE, c(1L, 1L), list(operator), NULL)
}

#' Indexes/slices an operator
#' @param operator a LinOp
#' @param dim the shape of the LinOp after indexing/slicing
#' @param keys the (row slice, column slice) 
#' @return a LinOp representing the index/slice
lo.index <- function(operator, dim, keys) {
  LinOp(INDEX, dim, list(operator), keys)
}


#' Discrete convolution of two vectors
#' @param lh_op the left-hand operator in the convolution
#' @param rh_op the right-hand operator in the convolution
#' @param dim the shape of the LinOp after convolution
#' @return a LinOp representing the convolution
lo.conv <- function(lh_op, rh_op, dim) {
  LinOp(CONV, dim, list(rh_op), lh_op)
}

#' Trnaspose an operator
#' @param operator a LinOp
#' @return a LinOp representing the transpose
lo.transpose <- function(operator) {
  if(length(operator$dim) < 2)
    operator
  else if(length(operator$dim) > 2)
    stop("Unimplemented")
  else {
    new_dim = c(operator$dim[2L], operator$dim[1L])
    LinOp(TRANSPOSE, new_dim, list(operator))
  }
}

#' Reshape an operator
#' @param operator the operator to reshape
#' @param dim the (rows, cols) of the reshaped operator
#' @return a LinOp representing the reshaped operator
lo.redim <- function(operator, dim) {
  LinOp(REDIM_EXPR, dim, list(operator))
}
lo.reshape <- lo.redim

#' Convert a LinOp vector into a diagonal matrix. (Note version different for cvxpy v1.5.2 and this is v1.3.0 version)
#' @param operator the operator to convert
#' @return a LinOp representing the matrix diagonal
lo.diag_vec <- function(operator) {
  new_dim <- rep(operator$dim[1L], 2L)
  LinOp(DIAG_VEC, new_dim, list(operator), NULL)
}

#' Convert the diagonal of a LinOp matrix to a vector
#' @param operator the operator to convert
#' @return a LinOp representing the matrix diagonal
lo.diag_mat <- function(operator) {
  new_dim = c(operator$dim[1L], 1L)
  LinOp(DIAG_MAT, new_dim, list(operator), NULL)
}


#' Create vectorized upper triangular portion of a square matrix (excluding diagonal)
#' @param operator the operator to convert
#' @return a LinOp representing the vectorized upper triangle
lo.upper_tri <- function(operator) {
  entries <- operator$dim[1L] * operator$dim[2L]
  new_dim <- c(as.integer(floor(entries - operator$dim[1L]) / 2), 1L)
  LinOp(UPPER_TRI, new_dim, list(operator), NULL)
}

#' Concatenate operators horizontally, (like `cbind`)
#' @param operators the operators to horizontally stack
#' @param dim the (rows, cols) of the stacked operators
#' @return a LinOp representing the stacked expression
lo.hstack <- function(operators, dim) {
  LinOp(HSTACK, dim, operators, NULL)
}

#' Concatenate operators vertically, (like `rbind`)
#' @param operators the operators to vertically stack
#' @param dim the (rows, cols) of the stacked operators
#' @return a LinOp representing the stacked expression
lo.vstack <- function(operators, dim) {
  LinOp(VSTACK, dim, operators, NULL)
}

#' Return the operator in the constraint left op right, i.e. rewrite as: left - right op 0
#' @param lh_op the operator
#' @param rh_op the right operator
#' @return the operator in the constraint
get_constr_expr <- function(lh_op, rh_op = NULL) {
  if(is.null(rh_op))
    lh_op
  else
    lo.sum_expr(list(lh_op, lo.neg_expr(rh_op)))
}

#' Create an internal equality constraint
#' @param lh_op the lhs operator
#' @param rh_op the rhs operator
#' @param constr_id the id of the constraint
#' @return the operator in the constraint
create_eq <- function(lh_op, rh_op = NULL, constr_id = get_id()) {
  expr <- get_constr_expr(lh_op, rh_op)
  LinEqConstr(expr, constr_id, lh_op$dim)
}

#' Create an internal less than or equal constraint
#' @param lh_op the lhs operator
#' @param rh_op the rhs operator
#' @param constr_id the id of the constraint
#' @return the operator in the constraint
create_leq <- function(lh_op, rh_op = NULL, constr_id = get_id()) {
  expr <- get_constr_expr(lh_op, rh_op)
  LinLeqConstr(expr, constr_id, lh_op$dim)
}

#' Create an internal greater than or equal constraint
#' @param lh_op the lhs operator
#' @param rh_op the rhs operator
#' @param constr_id the id of the constraint
#' @return the operator in the constraint
create_geq <- function(lh_op, rh_op = NULL, constr_id = get_id()) {
  if(!is.null(rh_op))
    rh_op <- lo.neg_expr(rh_op)
  create_leq(lo.neg_expr(lh_op), rh_op, constr_id)
}


## NARAS_NOTE:
## These routines below are recursive and this should be done in C++ for efficiency,
## also the ids and the dims can be returned as separate parallel lists as
## almost surely they need to be processed in parallel.
## For get_expr_vars, it also seems that a single named vector of dims with the names being var_ids
## would be more succint? Perhaps similar considerations apply to others.
## Questions:
## The recursive call will fail if there are no variables. Is it guaranteed that there will be some?
##  ACTUALLY THIS CODE IS ONLY USED IN TESTS, not in CVXR.
## Try a find-grep-dired for get_expr_vars...

#' Get a list of the variables in the operator and their shapes
#' @param operator the lhs operator
#' @return a list of var id, var dim pairs
get_expr_vars <- function(operator) {
  if(operator$type == VARIABLE)
    list(list(id = operator$data, dim = operator$dim))
  else {
    vars_ <- list()
    for(arg in operator$args)
      vars_ <- c(vars_, get_expr_vars(arg))
    vars_
  }
}

## NEVER_USED_IN_CVXR
#' Get a list of the parameters in the operator and their shapes
#' @param operator the lhs operator
#' @return a list of var id, var dim pairs
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

## NEVER_USED_IN_CVXR_ONLY_TEST
copy_constr <- function(constr, func) {
  expr <- func(constr$expr)
  new(class(constr), expr, constr$constr_id, constr$dim)
}

## NEVER_USED_IN_CVXR
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

## THIS IS USED IN CVXR in reductions in cvxpy/reductions/eval_params.py
replace_params_with_consts <- function(expr) {
  if(expr$type == PARAM)
    create_const(expr$data$value, expr$dim)
  else {
    new_args <- list()
    for(arg in expr$args)
      new_args <- c(new_args, replace_params_with_consts(arg))
    # Data could also be a parameter
    ## NARAS replaced "is" by "inherits" which is faster
    if(inherits(expr$data, "LinOp") && expr$data$type == PARAM) {
      data_lin_op <- expr$data
      data <- create_const(data_lin_op$data$value, data_lin_op$dim)
    } else
      data <- expr$data
    LinOp(expr$type, expr$dim, new_args, data)
  }
}


