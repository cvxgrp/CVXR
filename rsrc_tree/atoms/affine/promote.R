## CVXPY SOURCE: cvxpy/atoms/affine/promote.py

#'
#' The Promote class.
#' 
#' This class represents the promotion of a scalar expression into a vector/matrix.
#'
#' @slot expr An \linkS4class{Expression} or numeric constant.
#' @slot promoted_dim The desired dimensions.
#' @name Promote-class
#' @aliases Promote
#' @rdname Promote-class
.Promote <- setClass("Promote", representation(expr = "Expression", promoted_dim = "numeric"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} or numeric constant.
#' @param promoted_dim The desired dimensions.
#' @rdname Promote-class
Promote <- function(expr, promoted_dim) { .Promote(expr = expr, promoted_dim = promoted_dim) } 

promote <- function(expr, promoted_dim) { 
  expr <- as.Constant(expr)
  if(!all(dim(expr) == promoted_dim)) {
    if(!is_scalar(expr))
      stop("Only scalars may be promoted.")
    return(Promote(expr = expr, promoted_dim = promoted_dim))
  } else
    return(expr)
}

setMethod("initialize", "Promote", function(.Object, ..., expr, promoted_dim) {
  .Object@expr <- expr
  .Object@promoted_dim <- promoted_dim
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object A \linkS4class{Promote} object.
#' @param values A list containing the value to promote.
#' @describeIn Promote Promotes the value to the new dimensions.
setMethod("to_numeric", "Promote", function(object, values) {
  array(1, dim = object@promoted_dim) * as.vector(values[[1]])[1]
})

#' @describeIn Promote Is the expression symmetric?
setMethod("is_symmetric", "Promote", function(object) {
  ndim(object) == 2 && dim(object)[1] == dim(object)[2]
})

#' @describeIn Promote Returns the (row, col) dimensions of the expression.
setMethod("dim_from_args", "Promote", function(object) { object@promoted_dim })

#' @describeIn Promote Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "Promote", function(object) { TRUE })

#' @describeIn Promote Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "Promote", function(object) { TRUE })

#' @describeIn Promote Returns information needed to reconstruct the expression besides the args.
setMethod("get_data", "Promote", function(object) { list(object@promoted_dim) })

Promote.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lu.promote(arg_objs[[1]], dim), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Promote The graph implementation of the atom.
setMethod("graph_implementation", "Promote", function(object, arg_objs, dim, data = NA_real_) {
  Promote.graph_implementation(arg_objs, dim, data)
})
