## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/log.py

#'
#' The Log class.
#'
#' This class represents the elementwise natural logarithm \eqn{\log(x)}.
#'
#' @slot x An \linkS4class{Expression} or numeric constant.
#' @name Log-class
#' @aliases Log
#' @rdname Log-class
.Log <- setClass("Log", representation(x = "ConstValORExpr"), contains = "Elementwise")

#' @param x An \linkS4class{Expression} or numeric constant.
#' @rdname Log-class
Log <- function(x) { .Log(x = x) }

setMethod("initialize", "Log", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., atom_args = list(.Object@x))
})

#' @param object A \linkS4class{Log} object.
#' @param values A list of arguments to the atom.
#' @describeIn Log The elementwise natural logarithm of the input value.
setMethod("to_numeric", "Log", function(object, values) { log(values[[1]]) })

#' @describeIn Log The sign of the atom is unknown.
setMethod("sign_from_args", "Log", function(object) { c(FALSE, FALSE) })

#' @describeIn Log The atom is not convex.
setMethod("is_atom_convex", "Log", function(object) { FALSE })

#' @describeIn Log The atom is concave.
setMethod("is_atom_concave", "Log", function(object) { TRUE })

#' @describeIn Log Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Log", function(object) { FALSE })

#' @describeIn Log Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Log", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn Log The atom is weakly increasing.
setMethod("is_incr", "Log", function(object, idx) { TRUE })

#' @describeIn Log The atom is not weakly decreasing.
setMethod("is_decr", "Log", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn Log Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "Log", function(object, values) {
  rows <- size(object@args[[1]])
  cols <- size(object)

  # Outside domain or on boundary
  if(min(values[[1]]) <= 0)
    return(list(NA_real_))   # Non-differentiable
  else {
    grad_vals <- 1.0/values[[1]]
    return(list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols)))
  }
})

#' @describeIn Log Returns constraints describng the domain of the node
setMethod(".domain", "Log", function(object) { list(object@args[[1]] >= 0) })

