## CVXPY SOURCE: cvxpy/atoms/affine/cumsum.py
#
# Difference Matrix
#
# Returns a sparse matrix representation of the first order difference operator.
#
# @param dim The length of the matrix dimensions.
# @param axis The axis to take the difference along.
# @return A square matrix representing the first order difference.
## TODO: Recode in C++
get_diff_mat <- function(dim, axis) {
  # Construct a sparse matrix representation
  val_arr <- c()
  row_arr <- c()
  col_arr <- c()
  
  for(i in 1:dim) {
    val_arr <- c(val_arr, 1)
    row_arr <- c(row_arr, i)
    col_arr <- c(col_arr, i)
    if(i > 1) {
      val_arr <- c(val_arr, -1)
      row_arr <- c(row_arr, i)
      col_arr <- c(col_arr, i-1)
    }
  }
  
  mat <- sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = c(dim, dim))
  
  if(axis == 2)
    mat
  else
    t(mat)
}

#'
#' The CumSum class.
#'
#' This class represents the cumulative sum.
#'
#' @slot expr An \linkS4class{Expression} to be summed.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, and \code{2} indicates columns. The default is \code{2}.
#' @name CumSum-class
#' @aliases CumSum
#' @rdname CumSum-class
.CumSum <- setClass("CumSum", prototype = prototype(axis = 2), contains = c("AffAtom", "AxisAtom"))

#' @param expr An \linkS4class{Expression} to be summed.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, and \code{2} indicates columns. The default is \code{2}.
#' @rdname CumSum-class
CumSum <- function(expr, axis = 2) { .CumSum(expr = expr, axis = axis) }

#' @param object A \linkS4class{CumSum} object.
#' @param values A list of arguments to the atom.
#' @describeIn CumSum The cumulative sum of the values along the specified axis.
setMethod("to_numeric", "CumSum", function(object, values) {
  # apply(values[[1]], object@axis, base::cumsum)
  if(object@axis == 1)
    do.call(rbind, lapply(seq_len(nrow(values[[1]])), function(i) { base::cumsum(values[[1]][i,]) }))
  else if(object@axis == 2)
    do.call(cbind, lapply(seq_len(ncol(values[[1]])), function(j) { base::cumsum(values[[1]][,j]) }))
  else
    base::cumsum(values[[1]])
})

#' @describeIn CumSum The dimensions of the atom.
setMethod("dim_from_args", "CumSum", function(object) { dim(object@args[[1]]) })

#' @describeIn CumSum Returns the axis along which the cumulative sum is taken.
setMethod("get_data", "CumSum", function(object) { list(object@axis) })

#' @param values A list of numeric values for the arguments
#' @describeIn CumSum Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "CumSum", function(object, values) {
  # TODO: This is inefficient
  val_dim <- dim(values[[1]])
  collapse <- setdiff(1:length(val_dim), object@axis)
  dim <- val_dim[collapse]
  mat <- matrix(0, nrow = dim, ncol = dim)
  mat[lower.tri(mat, diag = TRUE)] <- 1

  # var <- Variable(dim(object@args[[1]]))
  var <- new("Variable", dim = dim(object@args[[1]]))
  if(object@axis == 2)
    grad <- .grad(new("MulExpression", lh_exp = mat, rh_exp = var), values)[[2]]
  else
    grad <- .grad(new("MulExpression", lh_exp = var, rh_exp = t(mat)), values)[[1]]
  list(grad)
})

CumSum.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  # Implicit O(n) definition:
  # X = Y[:1,:] - Y[1:,:]
  Y <- create_var(dim)
  axis <- data[[1]]
  collapse <- setdiff(1:length(dim), axis)
  new_dim <- dim[collapse]
  diff_mat <- get_diff_mat(new_dim, axis)
  diff_mat <- create_const(diff_mat, c(new_dim, new_dim), sparse = TRUE)

  if(axis == 2)
    diff <- lo.mul_expr(diff_mat, Y)
  else
    diff <- lo.rmul_expr(Y, diff_mat)
  list(Y, list(create_eq(arg_objs[[1]], diff)))
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn CumSum The graph implementation of the atom.
setMethod("graph_implementation", "CumSum", function(object, arg_objs, dim, data = NA_real_) {
  CumSum.graph_implementation(arg_objs, dim, data)
})

