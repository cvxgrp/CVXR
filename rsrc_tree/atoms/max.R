## CVXPY SOURCE: cvxpy/atoms/max.py
#'
#' The MaxEntries class.
#'
#' The maximum of an expression.
#'
#' @slot x An \linkS4class{Expression} representing a vector or matrix.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @slot keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @name MaxEntries-class
#' @aliases MaxEntries
#' @rdname MaxEntries-class
.MaxEntries <- setClass("MaxEntries", contains = "AxisAtom")

#' @param x An \linkS4class{Expression} representing a vector or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @param keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @rdname MaxEntries-class
MaxEntries <- function(x, axis = NA_real_, keepdims = FALSE) { .MaxEntries(expr = x, axis = axis, keepdims = keepdims) }

#' @param object A \linkS4class{MaxEntries} object.
#' @param values A list of arguments to the atom.
#' @describeIn MaxEntries The largest entry in \code{x}.
setMethod("to_numeric", "MaxEntries", function(object, values) {
  apply_with_keepdims(values[[1]], max, axis = object@axis, keepdims = object@keepdims)
})

#' @describeIn MaxEntries The sign of the atom.
setMethod("sign_from_args",  "MaxEntries", function(object) { c(is_nonneg(object@args[[1]]), is_nonpos(object@args[[1]])) })

#' @describeIn MaxEntries The atom is convex.
setMethod("is_atom_convex", "MaxEntries", function(object) { TRUE })

#' @describeIn MaxEntries The atom is not concave.
setMethod("is_atom_concave", "MaxEntries", function(object) { FALSE })

#' @describeIn MaxEntries Is the atom log-log convex.
setMethod("is_atom_log_log_convex", "MaxEntries", function(object) { TRUE })

#' @describeIn MaxEntries Is the atom log-log concave.
setMethod("is_atom_log_log_concave", "MaxEntries", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn MaxEntries The atom is weakly increasing in every argument.
setMethod("is_incr", "MaxEntries", function(object, idx) { TRUE })

#' @param idx An index into the atom.
#' @describeIn MaxEntries The atom is not weakly decreasing in any argument.
setMethod("is_decr", "MaxEntries", function(object, idx) { FALSE })

#' @describeIn MaxEntries Is \code{x} piecewise linear?
setMethod("is_pwl", "MaxEntries", function(object) { is_pwl(object@args[[1]]) })

#' @param values A list of numeric values for the arguments
#' @describeIn MaxEntries Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "MaxEntries", function(object, values) { .axis_grad(object, values) })

#' @param value A numeric value
#' @describeIn MaxEntries Gives the (sub/super)gradient of the atom w.r.t. each column variable
setMethod(".column_grad", "MaxEntries", function(object, value) {
  # Grad: 1 for a largest index
  value <- as.vector(value)
  idx <- (value == max(value))
  D <- matrix(0, nrow = length(value), ncol = 1)
  D[idx,1] <- 1
  D
})

