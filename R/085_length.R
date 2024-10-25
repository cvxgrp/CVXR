## CVXPY SOURCE: cvxpy/atoms/length.py
#'
#' The VecLength class.
#'
#' This class represents the length of a vector (index of last nonzero, ones-based).
#'
#' @slot x An \linkS4class{Expression} representing a vector.
#' @name VecLength-class
#' @aliases VecLength
#' @rdname VecLength-class
.VecLength <- setClass("VecLength", representation(x = "ConstValORExpr"),
                         validity = function(object) {
                           if(!is_vector(object@args[[1]]))
                             stop("[VecLength: x] The argument x must be a vector.")
                           return(TRUE)
                         }, contains = "Atom")

#' @param x An \linkS4class{Expression} representing a vector.
#' @rdname VecLength-class
VecLength <- function(x) { .Length(x = x) }

setMethod("initialize", "VecLength", function(.Object, ..., x) {
  .Object@x <- x
  callNextMethod(.Object, ..., atom_args = list(.Object@x))
})

#' @param object A \linkS4class{VecLength} object.
#' @param values A list of arguments to the atom.
#' @describeIn VecLength The length of the vector.
setMethod("to_numeric", "VecLength", function(object, values) {
  outside_tol <- abs(values[[1]]) > ATOM_EVAL_TOL
  return(base::max(which(outside_tol)))
})

#' @describeIn VecLength The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "VecLength", function(object) { c(1,1) })

#' @describeIn VecLength The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "VecLength", function(object) { c(TRUE, FALSE) })

#' @describeIn VecLength Is the atom convex?
setMethod("is_atom_convex", "VecLength", function(object) { FALSE })

#' @describeIn VecLength Is the atom concave?
setMethod("is_atom_concave", "VecLength", function(object) { FALSE })

#' @describeIn VecLength Is the atom quasiconvex?
setMethod("is_atom_quasiconvex", "VecLength", function(object) { TRUE })

#' @describeIn VecLength Is the atom quasiconcave?
setMethod("is_atom_quasiconcave", "VecLength", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn VecLength Is the atom weakly increasing in the index?
setMethod("is_incr", "VecLength", function(object, idx) { FALSE })

#' @describeIn VecLength Is the atom weakly decreasing in the index?
setMethod("is_decr", "VecLength", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn VecLength Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "VecLength", function(object, values) { NA_real_ })
