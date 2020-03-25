#'
#' The DGumbel class.
#'
#' This class represents the probability density function of the Gumbel distribution, i.e., 
#' \eqn{\exp(-(x + e^{-x}))} applied elementwise.
#'
#' @slot x An \linkS4class{Expression} or numeric constant.
#' @name DGumbel-class
#' @aliases DGumbel
#' @rdname DGumbel-class
.DGumbel <- setClass("DGumbel", representation(x = "ConstValORExpr"), contains = "Elementwise")

#' @param x An \linkS4class{Expression} or numeric constant.
#' @rdname DGumbel-class
DGumbel <- function(x) { .DGumbel(x = x) }

setMethod("initialize", "DGumbel", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., atom_args = list(.Object@x))
})

#' @param object A \linkS4class{DGumbel} object.
#' @param values A list of arguments to the atom.
#' @describeIn DGumbel The Gumbel probability density function evaluated elementwise.
setMethod("to_numeric", "DGumbel", function(object, values) { exp(-(values[[1]] + exp(-values[[1]]))) })

#' @describeIn DGumbel The atom is always positive.
setMethod("sign_from_args", "DGumbel", function(object) { c(TRUE, FALSE) })

#' @describeIn DGumbel The atom is not convex.
setMethod("is_atom_convex", "DGumbel", function(object) { FALSE })

#' @describeIn DGumbel The atom is not concave.
setMethod("is_atom_concave", "DGumbel", function(object) { FALSE })

#' @describeIn DGumbel Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "DGumbel", function(object) { FALSE })

#' @describeIn DGumbel Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "DGumbel", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn DGumbel The atom is monotonically increasing for \eqn{x \leq 0}.
setMethod("is_incr", "DGumbel", function(object, idx) { is_nonpos(object@args[[1]]) })

#' @describeIn DGumbel The atom is monotonically decreasing for \eqn{x \geq 0}.
setMethod("is_decr", "DGumbel", function(object, idx) { is_nonneg(object@args[[1]]) })

#' @param values A list of numeric values for the arguments
#' @describeIn DGumbel Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "DGumbel", function(object, values) {
  rows <- size(object@args[[1]])
  cols <- size(object)
  exp_val <- exp(-values[[1]])
  grad_vals <- exp(-(values[[1]] + exp_val))*(exp_val - 1)
  return(list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols)))
})
