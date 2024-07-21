## CVXPY SOURCE: cvxpy/atoms/quad_over_lin.py
#'
#' The QuadOverLin class.
#'
#' This class represents the sum of squared entries in X divided by a scalar y, \eqn{\sum_{i,j} X_{i,j}^2/y}.
#'
#' @slot x An \linkS4class{Expression} or numeric matrix.
#' @slot y A scalar \linkS4class{Expression} or numeric constant.
#' @name QuadOverLin-class
#' @aliases QuadOverLin
#' @rdname QuadOverLin-class
.QuadOverLin <- setClass("QuadOverLin", representation(x = "ConstValORExpr", y = "ConstValORExpr"), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric matrix.
#' @param y A scalar \linkS4class{Expression} or numeric constant.
#' @rdname QuadOverLin-class
QuadOverLin <- function(x, y) { .QuadOverLin(x = x, y = y) }

setMethod("initialize", "QuadOverLin", function(.Object, ..., x, y) {
  .Object@x <- x
  .Object@y <- y
  callNextMethod(.Object, ..., atom_args = list(.Object@x, .Object@y))
})

#' @describeIn QuadOverLin Does the atom handle complex numbers?
setMethod("allow_complex", "QuadOverLin", function(object) { TRUE })

#' @param object A \linkS4class{QuadOverLin} object.
#' @param values A list of arguments to the atom.
#' @describeIn QuadOverLin The sum of the entries of \code{x} squared over \code{y}.
setMethod("to_numeric", "QuadOverLin", function(object, values) { sum(Mod(values[[1]])^2) / values[[2]] })

#' @describeIn QuadOverLin Check the dimensions of the arguments.
setMethod("validate_args",   "QuadOverLin", function(object) {
  if(!is_scalar(object@args[[2]]))
    stop("The second argument to QuadOverLin must be a scalar.")
  if(is_complex(object@args[[2]]))
    stop("The second argument to QuadOverLin cannot be complex.")
  callNextMethod()
})

#' @describeIn QuadOverLin The atom is a scalar.
setMethod("dim_from_args", "QuadOverLin", function(object) { c(1,1) })

#' @describeIn QuadOverLin The atom is positive.
setMethod("sign_from_args",  "QuadOverLin", function(object) { c(TRUE, FALSE) })

#' @describeIn QuadOverLin The atom is convex.
setMethod("is_atom_convex", "QuadOverLin", function(object) { TRUE })

#' @describeIn QuadOverLin The atom is not concave.
setMethod("is_atom_concave", "QuadOverLin", function(object) { FALSE })

#' @describeIn QuadOverLin Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "QuadOverLin", function(object) { TRUE })

#' @describeIn QuadOverLin Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "QuadOverLin", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn QuadOverLin A logical value indicating whether the atom is weakly increasing in argument \code{idx}.
setMethod("is_incr", "QuadOverLin", function(object, idx) { (idx == 1) && is_nonneg(object@args[[idx]]) })

#' @describeIn QuadOverLin A logical value indicating whether the atom is weakly decreasing in argument \code{idx}.
setMethod("is_decr", "QuadOverLin", function(object, idx) { ((idx == 1) && is_nonpos(object@args[[idx]])) || (idx == 2) })

#' @describeIn QuadOverLin Quadratic if \code{x} is affine and \code{y} is constant.
setMethod("is_quadratic", "QuadOverLin", function(object) { is_affine(object@args[[1]]) && is_constant(object@args[[2]]) })

#' @describeIn QuadOverLin Quadratic of piecewise affine if \code{x} is piecewise linear and \code{y} is constant.
setMethod("is_qpwa", "QuadOverLin", function(object) { is_pwl(object@args[[1]]) && is_constant(object@args[[2]]) })

#' @describeIn QuadOverLin Returns constraints describing the domain of the node
setMethod(".domain", "QuadOverLin", function(object) { list(object@args[[2]] >= 0) })

#' @param values A list of numeric values for the arguments
#' @describeIn QuadOverLin Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "QuadOverLin", function(object, values) {
  X <- values[[1]]
  y <- as.vector(values[[2]])
  if(y <= 0)
    return(list(NA_real_, NA_real_))
  else {
    # DX = 2X/y, Dy = -||X||^2_2/y^2
    Dy <- -sum(Mod(X)^2)/y^2
    Dy <- Matrix(Dy, sparse = TRUE)
    DX <- 2.0*X/y
    DX <- Matrix(as.vector(t(DX)), sparse = TRUE)
    return(list(DX, Dy))
  }
})

