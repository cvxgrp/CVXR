## CVXPY SOURCE: cvxpy/atoms/affine/conj.py

#'
#' The Conjugate class.
#' 
#' This class represents the complex conjugate of an expression.
#' 
#' @slot expr An \linkS4class{Expression} or R numeric data.
#' @name Conjugate-class
#' @aliases Conjugate
#' @rdname Conjugate-class
.Conjugate <- setClass("Conjugate", representation(expr = "ConstValORExpr"), contains = "AffAtom")

#' @param expr An \linkS4class{Expression} or R numeric data.
#' @rdname Conjugate-class
Conjugate <- function(expr) { .Conjugate(expr = expr) }

setMethod("initialize", "Conjugate", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., atom_args = list(.Object@expr))
})

#' @param object A \linkS4class{Conjugate} object.
#' @param values A list of arguments to the atom.
#' @describeIn Conjugate Elementwise complex conjugate of the constant.
setMethod("to_numeric", "Conjugate", function(object, values) { Conj(values[[1]]) })

#' @describeIn Conjugate The (row, col) dimensions of the expression.
setMethod("dim_from_args", "Conjugate", function(object) { dim(object@args[[1]]) })

#' @param idx An index into the atom.
#' @describeIn Conjugate Is the composition weakly increasing in argument idx?
setMethod("is_incr", "Conjugate", function(object, idx) { FALSE })

#' @describeIn Conjugate Is the composition weakly decreasing in argument idx?
setMethod("is_decr", "Conjugate", function(object, idx) { FALSE })

#' @describeIn Conjugate Is the expression symmetric?
setMethod("is_symmetric", "Conjugate", function(object) { is_symmetric(object@args[[1]]) })

#' @describeIn Conjugate Is the expression hermitian?
setMethod("is_hermitian", "Conjugate", function(object) { is_hermitian(object@args[[1]]) })

Conjugate.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(arg_objs[[1]], list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Conv The graph implementation of the atom.
setMethod("graph_implementation", "Conjugate", function(object, arg_objs, dim, data = NA_real_) {
  Conjugate.graph_implementation(arg_objs, dim, data)
})

#'
#' The Conv class.
#'
#' This class represents the 1-D discrete convolution of two vectors.
#'
#' @slot lh_exp An \linkS4class{Expression} or R numeric data representing the left-hand vector.
#' @slot rh_exp An \linkS4class{Expression} or R numeric data representing the right-hand vector.
#' @name Conv-class
#' @aliases Conv
#' @rdname Conv-class
.Conv <- setClass("Conv", representation(lh_exp = "ConstValORExpr", rh_exp = "ConstValORExpr"), contains = "AffAtom")

#' @param lh_exp An \linkS4class{Expression} or R numeric data representing the left-hand vector.
#' @param rh_exp An \linkS4class{Expression} or R numeric data representing the right-hand vector.
#' @rdname Conv-class
Conv <- function(lh_exp, rh_exp) { .Conv(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("initialize", "Conv", function(.Object, ..., lh_exp, rh_exp) {
  .Object@lh_exp <- lh_exp
  .Object@rh_exp <- rh_exp
  callNextMethod(.Object, ..., atom_args = list(.Object@lh_exp, .Object@rh_exp))
})

#' @param object A \linkS4class{Conv} object.
#' @param values A list of arguments to the atom.
#' @describeIn Conv The convolution of the two values.
setMethod("to_numeric", "Conv", function(object, values) {
  .Call('_CVXR_cpp_convolve', PACKAGE = 'CVXR', as.vector(values[[1]]), as.vector(values[[2]]))
})

#' @describeIn Conv Check both arguments are vectors and the first is a constant.
setMethod("validate_args", "Conv", function(object) {
  if(!is_vector(object@args[[1]]) || !is_vector(object@args[[2]]))
    stop("The arguments to Conv must resolve to vectors.")
  if(!is_constant(object@args[[1]]))
    stop("The first argument to Conv must be constant.")
})

#' @describeIn Conv The dimensions of the atom.
setMethod("dim_from_args", "Conv", function(object) {
  lh_length <- dim(object@args[[1]])[1]
  rh_length <- dim(object@args[[2]])[1]
  c(lh_length + rh_length - 1, 1)
})

#' @describeIn Conv The sign of the atom.
setMethod("sign_from_args", "Conv", function(object) { mul_sign(object@args[[1]], object@args[[2]]) })

#' @param idx An index into the atom.
#' @describeIn Conv Is the left-hand expression positive?
setMethod("is_incr", "Conv", function(object, idx) { is_nonneg(object@args[[1]]) })

#' @param idx An index into the atom.
#' @describeIn Conv Is the left-hand expression negative?
setMethod("is_decr", "Conv", function(object, idx) { is_nonpos(object@args[[1]]) })

Conv.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  list(lo.conv(arg_objs[[1]], arg_objs[[2]], dim), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector representing the dimensions of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Conv The graph implementation of the atom.
setMethod("graph_implementation", "Conv", function(object, arg_objs, dim, data = NA_real_) {
  Conv.graph_implementation(arg_objs, dim, data)
})
