## CVXPY SOURCE: cvxpy/atoms/min.py
#'
#' The MinEntries class.
#'
#' The minimum of an expression.
#'
#' @slot x An \linkS4class{Expression} representing a vector or matrix.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @slot keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @name MinEntries-class
#' @aliases MinEntries
#' @rdname MinEntries-class
.MinEntries <- setClass("MinEntries", contains = "AxisAtom")

#' @param x An \linkS4class{Expression} representing a vector or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @param keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @rdname MinEntries-class
MinEntries <- function(x, axis = NA_real_, keepdims = FALSE) { .MinEntries(expr = x, axis = axis, keepdims = keepdims) }

#' @param object A \linkS4class{MinEntries} object.
#' @param values A list of arguments to the atom.
#' @describeIn MinEntries The largest entry in \code{x}.
setMethod("to_numeric", "MinEntries", function(object, values) {
  apply_with_keepdims(values[[1]], min, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn MinEntries The sign of the atom.
setMethod("sign_from_args",  "MinEntries", function(object) { c(is_nonneg(object@args[[1]]), is_nonpos(object@args[[1]])) })

#' @describeIn MinEntries The atom is not convex.
setMethod("is_atom_convex", "MinEntries", function(object) { FALSE })

#' @describeIn MinEntries The atom is concave.
setMethod("is_atom_concave", "MinEntries", function(object) { TRUE })

#' @describeIn MinEntries Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "MinEntries", function(object) { FALSE })

#' @describeIn MinEntries Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "MinEntries", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn MinEntries The atom is weakly increasing in every argument.
setMethod("is_incr", "MinEntries", function(object, idx) { TRUE })

#' @param idx An index into the atom.
#' @describeIn MinEntries The atom is not weakly decreasing in any argument.
setMethod("is_decr", "MinEntries", function(object, idx) { FALSE })

#' @describeIn MinEntries Is \code{x} piecewise linear?
setMethod("is_pwl", "MinEntries", function(object) { is_pwl(object@args[[1]]) })

#' @param values A list of numeric values for the arguments
#' @describeIn MinEntries Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "MinEntries", function(object, values) { .axis_grad(object, values) })

#' @param value A numeric value
#' @describeIn MinEntries Gives the (sub/super)gradient of the atom w.r.t. each column variable
setMethod(".column_grad", "MinEntries", function(object, value) {
  # Grad: 1 for a largest index
  value <- as.vector(value)
  idx <- (value == min(value))
  D <- matrix(0, nrow = length(value), ncol = 1)
  D[idx,1] <- 1
  D
})

