## CVXPY SOURCE: cvxpy/atoms/affine/elementwise/xexp.py
#'
#' The XExp class.
#'
#' This class represents elementwise \eqn{x\exp^x}.
#'
#' @slot x An \linkS4class{Expression}.
#' @name XExp-class
#' @aliases XExp
#' @rdname XExp-class
.XExp <- setClass("XExp", representation(x = "ConstValORExpr"), contains = "Elementwise")

#' @param x An \linkS4class{Expression}.
#' @rdname XExp-class
XExp <- function(x) { .XExp(x = x) }

setMethod("initialize", "XExp", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., atom_args = list(.Object@x))
})

#' @param object A \linkS4class{XExp} object.
#' @param values A list of arguments to the atom.
#' @describeIn XExp The numeric value of the expression.
setMethod("to_numeric", "XExp", function(object, values) {
  values[[1]]*base::exp(values[[1]])
})

#' @describeIn XExp The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "XExp", function(object) {
  # Depends on the sign of x.
  c(is_nonneg(object@args[[1]]), is_nonpos(object@args[[1]]))
})

#' @describeIn XExp Is the atom convex?
setMethod("is_atom_convex", "XExp", function(object) { is_nonneg(object@args[[1]]) })

#' @describeIn XExp Is the atom concave?
setMethod("is_atom_concave", "XExp", function(object) { FALSE })

#' @describeIn XExp Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "XExp", function(object) { TRUE })

#' @describeIn XExp Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "XExp", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn XExp Is the atom weakly increasing in the index?
setMethod("is_incr", "XExp", function(object, idx) { TRUE })

#' @describeIn XExp Is the atom weakly decreasing in the index?
setMethod("is_decr", "XExp", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn XExp Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "XExp", function(object, values) {
  rows <- size(object@args[[1]])
  cols <- size(object)
  grad_vals <- base::exp(values[[1]])*(1 + values[[1]])
  return(list(Elementwise.elemwise_grad_to_diag(grad_vals, rows, cols)))
})

#' @describeIn XExp Returns constraints describing the domain of the node
setMethod(".domain", "XExp", function(object) { list(object@args[[1]] >= 0) })
