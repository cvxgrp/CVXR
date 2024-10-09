## CVXPY SOURCE: cvxpy/atoms/affine/transpose.py

#'
#' The Transpose class.
#'
#' This class represents the matrix transpose.
#'
#' @name Transpose-class
#' @aliases Transpose
#' @rdname Transpose-class
.Transpose <- setClass("Transpose", representation(expr = "Expression", axes = "ConstValORNULL"), prototype(axes = NULL), contains = "AffAtom")

Transpose <- function(expr, axes = NULL) { .Transpose(expr = expr, axes = axes) }

setMethod("initialize", "Transpose", function(.Object, ..., expr, axes = NULL) {
  .Object@expr <- expr
  .Object@axes <- axes
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object A \linkS4class{Transpose} object.
#' @param values A list of arguments to the atom.
#' @describeIn Transpose The transpose of the given value.
setMethod("to_numeric", "Transpose", function(object, values) {
  if(is.vector(values[[1]]))
    return(t(values[[1]]))
  else if(is(values[[1]], "Matrix")) {
    if(!is.null(object@axes))
      stop("Cannot permute Matrix object axes to (", paste(object@axes, collapse = ","), ")")
    return(t(values[[1]]))
  } else
    return(aperm(values[[1]], perm = object@axes))
})

#' @describeIn Transpose Is the expression symmetric?
setMethod("is_symmetric", "Transpose", function(object) { is_symmetric(object@args[[1]]) })

#' @describeIn Transpose Is the expression skew symmetric?
setMethod("is_skew_symmetric", "Transpose", function(object) { is_skew_symmetric(object@args[[1]]) })

#' @describeIn Transpose Is the expression hermitian?
setMethod("is_hermitian", "Transpose", function(object) { is_hermitian(object@args[[1]]) })

#' @describeIn Transpose The dimensions of the atom.
setMethod("dim_from_args", "Transpose", function(object) {
  if(is.null(object@axes))
    rev(dim(object@args[[1]]))
  else
    dim(object@args[[1]])[object@axes]
})

#' @describeIn Transpose Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Transpose", function(object) { TRUE })

#' @describeIn Transpose Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Transpose", function(object) { TRUE })

#' @describeIn Transpose Returns the axes for transposition.
setMethod("get_data", "Transpose", function(object) { list(object@axes) })

Transpose.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lu.transpose(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Transpose The graph implementation of the atom.
setMethod("graph_implementation", "Transpose", function(object, arg_objs, dim, data = NA_real_) {
  Transpose.graph_implementation(arg_objs, dim, data)
})
