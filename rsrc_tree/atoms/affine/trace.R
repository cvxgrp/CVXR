## CVXPY SOURCE: cvxpy/atoms/affine/trace.py

#'
#' The Trace class.
#'
#' This class represents the sum of the diagonal entries in a matrix.
#'
#' @slot expr An \linkS4class{Expression} representing a matrix.
#' @name Trace-class
#' @aliases Trace
#' @rdname Trace-class
.Trace <- setClass("Trace", representation(expr = "Expression"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} representing a matrix.
#' @rdname Trace-class
Trace <- function(expr) { .Trace(expr = expr) }

setMethod("initialize", "Trace", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object A \linkS4class{Trace} object.
#' @param values A list of arguments to the atom.
#' @describeIn Trace Sum the diagonal entries.
setMethod("to_numeric", "Trace", function(object, values) { sum(diag(values[[1]])) })

#' @describeIn Trace Check the argument is a square matrix.
setMethod("validate_args", "Trace", function(object) {
  arg_dim <- dim(object@args[[1]])
  if(arg_dim[1] != arg_dim[2])
    stop("Argument to Trace must be a square matrix")
})

#' @describeIn Trace The atom is a scalar.
setMethod("dim_from_args", "Trace", function(object){ c(1,1) })

#' @describeIn Trace The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "Trace", function(object) {
  is_nonneg <- is_nonneg(object@args[[1]]) || is_psd(object@args[[1]])
  is_nonpos <- is_nonpos(object@args[[1]]) || is_nsd(object@args[[1]])
  c(is_nonneg, is_nonpos)
})

#' @describeIn Trace A logical value indicating whether the atom is real.
setMethod("is_real", "Trace", function(object) {
  is_real(object@args[[1]]) || is_hermitian(object@args[[1]])
})

#' @describeIn Trace A logical value indicating whether the atom is complex.
setMethod("is_complex", "Trace", function(object) { !is_real(object) })

#' @describeIn Trace Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Trace", function(object) { TRUE })

#' @describeIn Trace Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Trace", function(object) { FALSE })

Trace.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.trace(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Trace The graph implementation of the atom.
setMethod("graph_implementation", "Trace", function(object, arg_objs, dim, data = NA_real_) {
  Trace.graph_implementation(arg_objs, dim, data)
})

