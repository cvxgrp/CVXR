## CVXPY SOURCE: cvxpy/atoms/affine/sum.py

#'
#' The SumEntries class.
#'
#' This class represents the sum of all entries in a vector or matrix.
#'
#' @slot expr An \linkS4class{Expression} representing a vector or matrix.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @slot keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @name SumEntries-class
#' @aliases SumEntries
#' @rdname SumEntries-class
.SumEntries <- setClass("SumEntries", contains = c("AxisAtom", "AffAtom"))

#' @param expr An \linkS4class{Expression} representing a vector or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @param keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @rdname SumEntries-class
SumEntries <- function(expr, axis = NA_real_, keepdims = FALSE) { .SumEntries(expr = expr, axis = axis, keepdims = keepdims) }

#' @param object A \linkS4class{SumEntries} object.
#' @param values A list of arguments to the atom.
#' @describeIn SumEntries Sum the entries along the specified axis.
setMethod("to_numeric", "SumEntries", function(object, values) {
  apply_with_keepdims(values[[1]], sum, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn SumEntries Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "SumEntries", function(object) { TRUE })

#' @describeIn SumEntries Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "SumEntries", function(object) { FALSE })

SumEntries.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  # TODO: Handle keepdims properly by setting a 1-D vector's dimension to c(len, NA_integer_).
  axis <- data[[1]]
  keepdims <- data[[2]]
  if(is.na(axis))
    obj <- lu.sum_entries(arg_objs[[1]], dim)
  else if(axis == 1) {
    # if(keepdims)
    #  const_dim <- c(arg_objs[[1]]$dim[2], 1)
    # else
    #  const_dim <- c(arg_objs[[1]]$dim[2], NA_integer_)
    
    # Always treat result as a column vector.
    const_dim <- c(arg_objs[[1]]$dim[2], 1)
    ones <- lu.create_const(array(1, dim = const_dim), const_dim)
    obj <- lu.rmul_expr(arg_objs[[1]], ones, dim)
  } else {   # axis == 2
    # if(keepdims)
    #  const_dim <- c(1, arg_objs[[1]]$dim[1])
    # else
    #  const_dim <- c(arg_objs[[1]]$dim[1], NA_integer_)
    # ones <- lu.create_const(array(1, dim = const_dim), const_dim)
    # obj <- lu.mul_expr(ones, arg_objs[[1]], dim)
    
    if(keepdims) {
      # Keep result as a row vector.
      const_dim <- c(1, arg_objs[[1]]$dim[1])
      ones <- lu.create_const(array(1, dim = const_dim), const_dim)
      obj <- lu.mul_expr(ones, arg_objs[[1]], dim)
    } else {
      # Treat collapsed 1-D vector as a column vector.
      const_dim <- c(arg_objs[[1]]$dim[1], 1)
      ones <- lu.create_const(array(1, dim = const_dim), const_dim)
      obj <- lu.rmul_expr(lo.transpose(arg_objs[[1]]), ones, dim)
    }
  }
  list(obj, list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn SumEntries The graph implementation of the atom.
setMethod("graph_implementation", "SumEntries", function(object, arg_objs, dim, data = NA_real_) {
  SumEntries.graph_implementation(arg_objs, dim, data)
})
