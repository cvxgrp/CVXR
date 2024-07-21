## CVXPY SOURCE: cvxpy/atoms/affine/real.py

#'
#' The Real class.
#'
#' This class represents the real part of an expression.
#'
#' @slot expr An \linkS4class{Expression} representing a vector or matrix.
#' @name Real-class
#' @aliases Real
#' @rdname Real-class
.Real <- setClass("Real", representation(expr = "Expression"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} representing a vector or matrix.
#' @rdname Real-class
Real <- function(expr) { .Real(expr = expr) }

setMethod("initialize", "Real", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object An \linkS4class{Real} object.
#' @param values A list of arguments to the atom.
#' @describeIn Real The imaginary part of the given value.
setMethod("to_numeric", "Real", function(object, values) { Re(values[[1]]) })

#' @describeIn Real The dimensions of the atom.
setMethod("dim_from_args", "Real", function(object) { dim(object@args[[1]]) })

#' @describeIn Real Is the atom imaginary?
setMethod("is_imag", "Real", function(object) { FALSE })

#' @describeIn Real Is the atom complex valued?
setMethod("is_complex", "Real", function(object) { FALSE })

#' @describeIn Real Is the atom symmetric?
setMethod("is_symmetric", "Real", function(object) { is_hermitian(object@args[[1]]) })
