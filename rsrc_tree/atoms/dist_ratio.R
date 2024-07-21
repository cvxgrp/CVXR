## CVXPY SOURCE: cvxpy/atoms/dist_ratio.py

#'
#' The DistRatio class.
#'
#' This class represents the ratio of the \eqn{l_2}-norm distance from \eqn{x} to two points \eqn{a} and \eqn{b}, i.e.,
#'
#' \deqn{||x - a||_2 / ||x - b||_2},
#'
#' where \eqn{a} and \eqn{b} are constants with \eqn{norm(x - a)_2 \leq norm(x - b)}.
#'
#' @slot x An \linkS4class{Expression}.
#' @name DistRatio-class
#' @aliases DistRatio
#' @rdname DistRatio-class
.DistRatio <- setClass("DistRatio", representation(x = "ConstValORExpr", a = "numeric", b = "numeric"),
                       validity = function(object) {
                         if(!is_constant(object@args[[2]]))
                           stop("[DistRatio: a] The argument a must be a constant.")
                         if(!is_constant(object@args[[3]]))
                           stop("[DistRatio: b] The argument b must be a constant.")
                         return(TRUE)
                       }, contains = "Atom")

#' @param x An \linkS4class{Expression}.
#' @param a A numeric value.
#' @param b A numeric value.
#' @rdname DistRatio-class
DistRatio <- function(x, a, b) { .DistRatio(x = x, a = a, b = b) }

setMethod("initialize", "DistRatio", function(.Object, ..., x, a, b) {
  .Object@x <- x
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@x, a, b))
  .Object@a <- value(.Object@args[[2]])
  .Object@b <- value(.Object@args[[3]])
  .Object
})

#' @param object A \linkS4class{DistRatio} object.
#' @param values A list of arguments to the atom.
#' @describeIn DistRatio The distance ratio of the expression.
setMethod("to_numeric", "DistRatio", function(object, values) {
  sqrt(sum((values[[1]] - object@a)^2)) / sqrt(sum((values[[1]] - object@b)^2))
})

#' @describeIn DistRatio The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "DistRatio", function(object) { c(1,1) })

#' @describeIn DistRatio The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "DistRatio", function(object) { c(TRUE, FALSE) })

#' @describeIn DistRatio Is the atom convex?
setMethod("is_atom_convex", "DistRatio", function(object) { FALSE })

#' @describeIn DistRatio Is the atom concave?
setMethod("is_atom_concave", "DistRatio", function(object) { FALSE })

#' @describeIn DistRatio Is the atom quasiconvex?
setMethod("is_atom_quasiconvex", "DistRatio", function(object) { TRUE })

#' @describeIn DistRatio Is the atom quasiconcave?
setMethod("is_atom_quasiconcave", "DistRatio", function(object) { FALSE })

#' @param idx An index into the atom.
#' @describeIn DistRatio Is the atom weakly increasing in the index?
setMethod("is_incr", "DistRatio", function(object, idx) { FALSE })

#' @describeIn DistRatio Is the atom weakly decreasing in the index?
setMethod("is_decr", "DistRatio", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn DistRatio Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "DistRatio", function(object, values) { return(NA_real_) })

