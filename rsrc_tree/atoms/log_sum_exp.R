## CVXPY SOURCE: cvxpy/atoms/log_sum_exp.py
#'
#' The LogSumExp class.
#'
#' The natural logarithm of the sum of the elementwise exponential, \eqn{\log\sum_{i=1}^n e^{x_i}}.
#'
#' @slot x An \linkS4class{Expression} representing a vector or matrix.
#' @slot axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @slot keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @name LogSumExp-class
#' @aliases LogSumExp
#' @rdname LogSumExp-class
.LogSumExp <- setClass("LogSumExp", contains = "AxisAtom")

#' @param x An \linkS4class{Expression} representing a vector or matrix.
#' @param axis (Optional) The dimension across which to apply the function: \code{1} indicates rows, \code{2} indicates columns, and \code{NA} indicates rows and columns. The default is \code{NA}.
#' @param keepdims (Optional) Should dimensions be maintained when applying the atom along an axis? If \code{FALSE}, result will be collapsed into an \eqn{n x 1} column vector. The default is \code{FALSE}.
#' @rdname LogSumExp-class
LogSumExp <- function(x, axis = NA_real_, keepdims = FALSE) { .LogSumExp(expr = x, axis = axis, keepdims = keepdims) }

#' @param object A \linkS4class{LogSumExp} object.
#' @param values A list of arguments to the atom.
#' @describeIn LogSumExp Evaluates \eqn{e^x} elementwise, sums, and takes the natural log.
setMethod("to_numeric", "LogSumExp", function(object, values) {
  if(is.na(object@axis))
    log(sum(exp(values[[1]])))
  else
    # log(apply(exp(values[[1]]), object@axis, sum))
    log(apply_with_keepdims(exp(values[[1]]), sum, axis = object@axis, keepdims = object@keepdims))
})

#' @param values A list of numeric values.
#' @describeIn LogSumExp Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "LogSumExp", function(object, values) { .axis_grad(object, values) })

#' @param value A numeric value.
#' @describeIn LogSumExp Gives the (sub/super)gradient of the atom w.r.t. each column variable.
setMethod(".column_grad", "LogSumExp", function(object, value) {
  denom <- sum(exp(value))
  nom <- exp(value)
  D <- nom/denom
  D
})

#' @describeIn LogSumExp Returns sign (is positive, is negative) of the atom.
setMethod("sign_from_args",  "LogSumExp", function(object) { c(FALSE, FALSE) })

#' @describeIn LogSumExp The atom is convex.
setMethod("is_atom_convex", "LogSumExp", function(object) { TRUE })

#' @describeIn LogSumExp The atom is not concave.
setMethod("is_atom_concave", "LogSumExp", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn LogSumExp The atom is weakly increasing in the index.
setMethod("is_incr", "LogSumExp", function(object, idx) { TRUE })

#' @param idx An index into the atom.
#' @describeIn LogSumExp The atom is not weakly decreasing in the index.
setMethod("is_decr", "LogSumExp", function(object, idx) { FALSE })
