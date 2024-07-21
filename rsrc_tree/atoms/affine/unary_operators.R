## CVXPY SOURCE: cvxpy/atoms/affine/unary_operators.py

#'
#' The UnaryOperator class.
#'
#' This base class represents expressions involving unary operators.
#'
#' @slot expr The \linkS4class{Expression} that is being operated upon.
#' @slot op_name A \code{character} string indicating the unary operation.
#' @name UnaryOperator-class
#' @aliases UnaryOperator
#' @rdname UnaryOperator-class
UnaryOperator <- setClass("UnaryOperator", representation(expr = "Expression"), contains = "AffAtom")

setMethod("initialize", "UnaryOperator", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

setMethod("op_name", "UnaryOperator", function(object) { stop("Unimplemented") })
setMethod("op_func", "UnaryOperator", function(object) { stop("Unimplemented") })

#' @param x,object A \linkS4class{UnaryOperator} object.
#' @describeIn UnaryOperator Returns the expression in string form.
setMethod("name", "UnaryOperator", function(x) {
  paste(op_name(x), name(x@args[[1]]), sep = "")
})

#' @param values A list of arguments to the atom.
#' @describeIn UnaryOperator Applies the unary operator to the value.
setMethod("to_numeric", "UnaryOperator", function(object, values) {
  op_func(object)(values[[1]])
})

#'
#' The NegExpression class.
#'
#' This class represents the negation of an affine expression.
#'
#' @name NegExpression-class
#' @aliases NegExpression
#' @rdname NegExpression-class
.NegExpression <- setClass("NegExpression", contains = "UnaryOperator")
NegExpression <- function(expr) { .NegExpression(expr = expr) }

setMethod("op_name", "NegExpression", function(object) { "-" })
setMethod("op_func", "NegExpression", function(object) { function(x) { -x } })

#' @param object A \linkS4class{NegExpression} object.
#' @describeIn NegExpression The (row, col) dimensions of the expression.
setMethod("dim_from_args", "NegExpression", function(object) { dim(object@args[[1]]) })

#' @describeIn NegExpression The (is positive, is negative) sign of the expression.
setMethod("sign_from_args", "NegExpression", function(object) { c(is_nonpos(object@args[[1]]), is_nonneg(object@args[[1]])) })

#' @param idx An index into the atom.
#' @describeIn NegExpression The expression is not weakly increasing in any argument.
setMethod("is_incr", "NegExpression", function(object, idx) { FALSE })

#' @describeIn NegExpression The expression is weakly decreasing in every argument.
setMethod("is_decr", "NegExpression", function(object, idx) { TRUE })

#' @describeIn NegExpression Is the expression symmetric?
setMethod("is_symmetric", "NegExpression", function(object) { is_symmetric(object@args[[1]]) })

#' @describeIn NegExpression Is the expression Hermitian?
setMethod("is_hermitian", "NegExpression", function(object) { is_hermitian(object@args[[1]]) })

NegExpression.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.neg_expr(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn NegExpression The graph implementation of the expression.
setMethod("graph_implementation", "NegExpression", function(object, arg_objs, dim, data = NA_real_) {
  NegExpression.graph_implementation(arg_objs, dim, data)
})
