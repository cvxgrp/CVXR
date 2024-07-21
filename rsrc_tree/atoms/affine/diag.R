## CVXPY SOURCE: cvxpy/atoms/affine/diag.py
#'
#' The DiagVec class.
#'
#' This class represents the conversion of a vector into a diagonal matrix.
#'
#' @slot expr An \linkS4class{Expression} representing the vector to convert.
#' @name DiagVec-class
#' @aliases DiagVec
#' @rdname DiagVec-class
.DiagVec <- setClass("DiagVec", representation(expr = "Expression"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} representing the vector to convert.
#' @rdname DiagVec-class
DiagVec <- function(expr) { .DiagVec(expr = expr) }

setMethod("initialize", "DiagVec", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object A \linkS4class{DiagVec} object.
#' @param values A list of arguments to the atom.
#' @describeIn DiagVec Convert the vector constant into a diagonal matrix.
setMethod("to_numeric", "DiagVec", function(object, values) { diag(as.vector(values[[1]])) })

#' @describeIn DiagVec The dimensions of the atom.
setMethod("dim_from_args", "DiagVec", function(object) {
  rows <- dim(object@args[[1]])[1]
  c(rows, rows)
})

#' @describeIn DiagVec Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "DiagVec", function(object) { TRUE })

#' @describeIn DiagVec Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "DiagVec", function(object) { TRUE })

#' @describeIn DiagVec Is the expression symmetric?
setMethod("is_symmetric", "DiagVec", function(object) { TRUE })

#' @describeIn DiagVec Is the expression hermitian?
setMethod("is_hermitian", "DiagVec", function(object) { TRUE })

#' @describeIn DiagVec Is the atom positive semidefinite?
setMethod("is_psd", "DiagVec", function(object) { is_nonneg(object) })

#' @describeIn DiagVec Is the atom negative semidefinite?
setMethod("is_nsd", "DiagVec", function(object) { is_nonpos(object) })

DiagVec.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.diag_vec(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn DiagVec The graph implementation of the atom.
setMethod("graph_implementation", "DiagVec", function(object, arg_objs, dim, data = NA_real_) {
  DiagVec.graph_implementation(arg_objs, dim, data)
})

#'
#' The DiagMat class.
#'
#' This class represents the extraction of the diagonal from a square matrix.
#'
#' @slot expr An \linkS4class{Expression} representing the matrix whose diagonal we are interested in.
#' @name DiagMat-class
#' @aliases DiagMat
#' @rdname DiagMat-class
.DiagMat <- setClass("DiagMat", representation(expr = "Expression"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} representing the matrix whose diagonal we are interested in.
#' @rdname DiagMat-class
DiagMat <- function(expr) { .DiagMat(expr = expr) }

setMethod("initialize", "DiagMat", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object A \linkS4class{DiagMat} object.
#' @param values A list of arguments to the atom.
#' @describeIn DiagMat Extract the diagonal from a square matrix constant.
setMethod("to_numeric", "DiagMat", function(object, values) { diag(values[[1]]) })

#' @describeIn DiagMat The size of the atom.
setMethod("dim_from_args", "DiagMat", function(object) {
  rows <- dim(object@args[[1]])[1]
  c(rows, 1)
})

#' @describeIn DiagMat Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "DiagMat", function(object) { TRUE })

#' @describeIn DiagMat Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "DiagMat", function(object) { TRUE })

#' @describeIn DiagMat A logical value indicating whether the atom is nonnegative.
setMethod("is_nonneg", "DiagMat", function(object) {
  is_nonneg(object@args[[1]]) || is_psd(object@args[[1]])
})

DiagMat.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.diag_mat(arg_objs[[1]]), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn DiagMat The graph implementation of the atom.
setMethod("graph_implementation", "DiagMat", function(object, arg_objs, dim, data = NA_real_) {
  DiagMat.graph_implementation(arg_objs, dim, data)
})

#' 
#' The Diag atom.
#' 
#' Turns an expression into a DiagVec object
#' 
#' @param expr An \linkS4class{Expression} that represents a vector or square matrix.
#' @return An \linkS4class{Expression} representing the diagonal vector/matrix.
#' @rdname Diag-int
Diag <- function(expr) {
  expr <- as.Constant(expr)
  if(is_vector(expr))
    return(DiagVec(Vec(expr)))
  else if(ndim(expr) == 2 && nrow(expr) == ncol(expr))
    return(DiagMat(expr = expr))
  else
    stop("Argument to Diag must be a vector or square matrix.")
}

