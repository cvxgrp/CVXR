#'
#' The AffAtom class.
#'
#' This virtual class represents an affine atomic expression.
#'
#' @aliases AffAtom
#' @export
AffAtom <- setClass("AffAtom", contains = c("VIRTUAL", "Atom"))

setMethod("func_curvature", "AffAtom", function(object) { Curvature(curvature = CURV_AFFINE_KEY) })

setMethod("sign_from_args", "AffAtom", function(object) { 
  arg_signs <- lapply(object@.args, function(arg) { arg@dcp_attr@sign })
  Reduce("+", arg_signs)
})

setMethod("monotonicity", "AffAtom", function(object) { rep(INCREASING, length(object@.args)) })

#'
#' The AddExpression class.
#'
#' This class represents the sum of any number of expressions.
#'
#' @slot arg_groups A \code{list} of \S4class{Expression}s and numeric data.frame, matrix, or vector objects.
#' @aliases AddExpression
#' @export
AddExpression <- setClass("AddExpression", representation(arg_groups = "list"), prototype(arg_groups = list()), contains = "AffAtom")

setMethod("init_dcp_attr", "AddExpression", function(object) {
  arg_dcp <- lapply(object@.args, function(arg) { arg@dcp_attr })
  Reduce("+", arg_dcp)
})

setMethod("initialize", "AddExpression", function(.Object, ..., arg_groups = list()) {
  .Object@arg_groups <- arg_groups
  .Object <- callNextMethod(.Object, ..., .args = arg_groups)   # Casts R values to Constant objects
  .Object@.args <- lapply(arg_groups, function(group) { if(is(group,"AddExpression")) group@.args else group })
  .Object@.args <- flatten_list(.Object@.args)   # Need to flatten list of expressions
  return(.Object)
})

AddExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  arg_objs <- lapply(arg_objs, function(arg) { if(size(arg) != size) promote(arg, size) else arg })
  list(sum_expr(arg_objs), list())
}

#'
#' The UnaryOperator class.
#'
#' This base class represents expressions involving unary operators.
#'
#' @slot expr The \S4class{Expression} that is being operated upon.
#' @slot op_name A \code{character} string indicating the unary operation.
#' @aliases UnaryOperator
#' @export
UnaryOperator <- setClass("UnaryOperator", representation(expr = "Expression", op_name = "character"), contains = "AffAtom")

setMethod("init_dcp_attr", "UnaryOperator", function(object) {
  .Primitive(object@op_name)(object@.args[[1]]@dcp_attr)
})

setMethod("initialize", "UnaryOperator", function(.Object, ..., expr, op_name) {
  .Object@expr = expr
  .Object@op_name = op_name
  callNextMethod(.Object, ..., .args = list(.Object@expr))
})

#'
#' The NegExpression class.
#'
#' This class represents the negation of an affine expression.
#'
#' @aliases NegExpression
#' @export
NegExpression <- setClass("NegExpression", contains = "UnaryOperator")

setMethod("initialize", "NegExpression", function(.Object, ...) {
  callNextMethod(.Object, ..., op_name = "-")
})

NegExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(neg_expr(arg_objs[[1]]), list())
}

#'
#' The BinaryOperator class.
#'
#' This base class represents expressions involving binary operators.
#'
#' @slot lh_exp The \S4class{Expression} on the left-hand side of the operator.
#' @slot rh_exp The \S4class{Expression} on the right-hand side of the operator.
#' @slot op_name A \code{character} string indicating the binary operation.
#' @aliases BinaryOperator
#' @export
BinaryOperator <- setClass("BinaryOperator", representation(lh_exp = "Expression", rh_exp = "Expression", op_name = "character"), contains = "AffAtom")

setMethod("init_dcp_attr", "BinaryOperator", function(object) {
  .Primitive(object@op_name)(object@.args[[1]]@dcp_attr, object@.args[[2]]@dcp_attr)
})

setMethod("validate_args", "BinaryOperator", function(object) {
  .Primitive(object@op_name)(object@.args[[1]]@dcp_attr@shape, object@.args[[2]]@dcp_attr@shape)
})

setMethod("initialize", "BinaryOperator", function(.Object, ..., lh_exp, rh_exp, op_name) {
  .Object@lh_exp = lh_exp
  .Object@rh_exp = rh_exp
  .Object@op_name = op_name
  callNextMethod(.Object, ..., .args = list(.Object@lh_exp, .Object@rh_exp))
})

#'
#' The MulExpression class.
#'
#' This class represents the product of two linear expressions.
#'
#' @aliases MulExpression
#' @export
MulExpression <- setClass("MulExpression", contains = "BinaryOperator")

setMethod("initialize", "MulExpression", function(.Object, ...) {
  callNextMethod(.Object, ..., op_name = "*")
})

MulExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  if(size[2] != 1 && all(size(arg_objs[[2]]) == c(1,1))) {
    arg <- promote(arg_objs[[2]], list(size[2], 1))
    arg_objs[2] <- diag_vec(arg)
  }
  list(mul_expr(arg_objs[[1]], arg_objs[[2]], size), list())
}

#'
#' The RMulExpression class.
#'
#' This class represents product of an expression with a constant on the right.
#'
#' @aliases RMulExpression
#' @export
RMulExpression <- setClass("RMulExpression", contains = "MulExpression")

RMulExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  if(size[1] != 1 && all(size(arg_objs[[1]]) == c(1,1))) {
    arg <- promote(arg_objs[[1]], list(size[1], 1))
    arg_objs[[1]] <- diag_vec(arg)
  }
  list(rmul_expr(arg_objs[[1]], arg_objs[[2]], size), list())
}

#'
#' The DivExpression class.
#'
#' This class represents one expression divided by another expression.
#'
#' @aliases DivExpression
#' @export
DivExpression <- setClass("DivExpression", contains = "BinaryOperator")

setMethod("initialize", "DivExpression", function(.Object, ...) {
  .Object@op_name = "/"
  callNextMethod(.Object, ..., op_name = "/")
})

DivExpression.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(div_expr(arg_objs[[1]], arg_objs[[2]]), list())
}

#'
#' The Reshape class.
#'
#' This class represents the reshaping of an expression. The operator vectorizes the expression,
#' then unvectorizes it into the new shape. Entries are stored in column-major order.
#'
#' @slot expr The \S4class{Expression} to reshape.
#' @slot rows The new number of rows.
#' @slot cols The new number of columns.
#' @aliases Reshape
#' @export
Reshape <- setClass("Reshape", representation(expr = "Expression", rows = "numeric", cols = "numeric"), contains = "AffAtom")

setMethod("validate_args", "Reshape", function(object) {
  old_len <- size(object@.args[1])[1] * size(object@.args[1])[2]
  new_len <- object@rows * object@cols
  if(old_len != new_len)
    stop(sprintf("Invalid reshape dimensions (%i, %i)", object@rows, object@cols))
})

setMethod("initialize", "Reshape", function(.Object, ..., expr, rows, cols) {
  .Object@rows <- rows
  .Object@cols <- cols
  .Object@expr <- expr
  callNextMethod(.Object, ..., .args = list(.Object@expr))
})

setMethod("shape_from_args", "Reshape", function(object) {
  Shape(rows = object@rows, cols = object@cols)
})

Reshape.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(reshape(arg_objs[1], size), list())
}

#'
#' The SumEntries class.
#'
#' This class represents sum of all the entries in an expression.
#'
#' @slot expr The \S4class{Expression} to sum the entries of.
#' @aliases SumEntries
#' @export
SumEntries <- setClass("SumEntries", representation(expr = "Expression"), contains = "AffAtom")

setMethod("initialize", "SumEntries", function(.Object, ..., expr) {
  .Object@expr = expr
  callNextMethod(.Object, ..., .args = list(.Object@expr))
})

setMethod("shape_from_args", "SumEntries", function(object){
  Shape(rows = 1, cols = 1)
})

SumEntries.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(sum_entries(arg_objs[[1]]), list())
}

#'
#' The Trace class.
#'
#' This class represents the sum of the diagonal entries in a matrix.
#'
#' @slot expr The \S4class{Expression} to sum the diagonal of.
#' @aliases Trace
#' @export
Trace <- setClass("Trace", representation(expr = "Expression"), contains = "AffAtom")

setMethod("validate_args", "Trace", function(object) {
  size <- size(object@.args[1])
  if(size[1] != size[2])
    stop("Argument to trace must be a square matrix")
})

setMethod("initialize", "Trace", function(.Object, ..., expr) {
  .Object@expr = expr
  callNextMethod(.Object, ..., .args = list(.Object@expr))
})

setMethod("shape_from_args", "Trace", function(object){
  Shape(rows = 1, cols = 1)
})

Trace.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(trace(arg_objs[[1]]), list())
}

#'
#' The Transpose class.
#'
#' This class represents the matrix transpose.
#'
#' @aliases Transpose
#' @export
Transpose <- setClass("Transpose", contains = "AffAtom")

setMethod("shape_from_args", "Transpose", function(object) {
  obj_size = size(object@.args[[1]])
  Shape(rows = obj_size[2], cols = obj_size[1])
})

Transpose.graph_implementation <- function(arg_objs, size, data = NA_real_) {
  list(transpose(arg_objs[[1]]), list())  
}
