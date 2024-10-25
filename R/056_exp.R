## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/exp.py
#'
#' The Exp class.
#'
#' This class represents the elementwise natural exponential \eqn{e^x}.
#'
#' @slot x An \linkS4class{Expression} object.
#' @name Exp-class
#' @aliases Exp
#' @rdname Exp-class
.Exp <- setClass("Exp", representation(x = "Expression"), contains = "Elementwise")

#' @param x An \linkS4class{Expression} object.
#' @rdname Exp-class
Exp <- function(x) { .Exp(x = x) }

setMethod("initialize", "Exp", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., atom_args = list(.Object@x))
})

#' @param object An \linkS4class{Exp} object.
#' @param values A list of arguments to the atom.
#' @describeIn Exp The matrix with each element exponentiated.
setMethod("to_numeric", "Exp", function(object, values) { exp(values[[1]]) })

#' @describeIn Exp The atom is positive.
setMethod("sign_from_args", "Exp", function(object) { c(TRUE, FALSE) })

#' @describeIn Exp The atom is convex.
setMethod("is_atom_convex", "Exp", function(object) { TRUE })

#' @describeIn Exp The atom is not concave.
setMethod("is_atom_concave", "Exp", function(object) { FALSE })

#' @describeIn Exp Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Exp", function(object) { TRUE })

#' @describeIn Exp Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Exp", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn Exp The atom is weakly increasing.
setMethod("is_incr", "Exp", function(object, idx) { TRUE })

#' @describeIn Exp The atom is not weakly decreasing.
setMethod("is_decr", "Exp", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn Exp Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "Exp", function(object, values) {
  rows <- size(object@args[[1]])
  cols <- size(object)
  grad_vals <- exp(values[[1]])
  list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols))
})
