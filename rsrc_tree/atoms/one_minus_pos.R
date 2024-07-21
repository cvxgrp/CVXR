## CVXPY SOURCE: cvxpy/atoms/one_minus_pos.py

#'
#' The OneMinusPos class.
#'
#' This class represents the difference \eqn{1 - x} with domain \eqn{\{x : 0 < x < 1}\}
#'
#' @slot x An \linkS4class{Expression} or numeric matrix.
#' @name OneMinusPos-class
#' @aliases OneMinusPos
#' @rdname OneMinusPos-class
.OneMinusPos <- setClass("OneMinusPos", representation(x = "ConstValORExpr", .ones = "ConstVal"), prototype(.ones = NA_real_), contains = "Atom")

#' @param x An \linkS4class{Expression} or numeric matrix.
#' @rdname OneMinusPos-class
OneMinusPos <- function(x) { .OneMinusPos(x = x) }

setMethod("initialize", "OneMinusPos", function(.Object, ..., x) {
  .Object@x <- x
  .Object@.ones <- matrix(1, nrow = nrow(x), ncol = ncol(x))
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@x))
  .Object@args[[1]] <- x
  .Object
})

#' @describeIn OneMinusPos The name and arguments of the atom.
setMethod("name", "OneMinusPos", function(x) { paste(class(x), x@args[[1]]) })

#' @param object A \linkS4class{OneMinusPos} object.
#' @param values A list of arguments to the atom.
#' @describeIn OneMinusPos Returns one minus the value.
setMethod("to_numeric", "OneMinusPos", function(object, values) { object@.ones - values[[1]] })

#' @describeIn OneMinusPos The dimensions of the atom.
setMethod("dim_from_args", "OneMinusPos", function(object) { dim(object@args[[1]]) })

#' @describeIn OneMinusPos Returns the sign (is positive, is negative) of the atom.
setMethod("sign_from_args", "OneMinusPos", function(object) { c(TRUE, FALSE) })

#' @describeIn OneMinusPos Is the atom convex?
setMethod("is_atom_convex", "OneMinusPos", function(object) { FALSE })

#' @describeIn OneMinusPos Is the atom concave?
setMethod("is_atom_concave", "OneMinusPos", function(object) { FALSE })

#' @describeIn OneMinusPos Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "OneMinusPos", function(object) { FALSE })

#' @describeIn OneMinusPos Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "OneMinusPos", function(object) { TRUE })

#' @param idx An index into the atom.
#' @describeIn OneMinusPos Is the atom weakly increasing in the argument \code{idx}?
setMethod("is_incr", "OneMinusPos", function(object, idx) { FALSE })

#' @param idx An index into the atom.
#' @describeIn OneMinusPos Is the atom weakly decreasing in the argument \code{idx}?
setMethod("is_decr", "OneMinusPos", function(object, idx) { TRUE })

#' @param values A list of numeric values for the arguments
#' @describeIn OneMinusPos Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "OneMinusPos", function(object, values) { Matrix(-object@.ones, sparse = TRUE) })

#'
#' The DiffPos atom.
#'
#' The difference between expressions, \eqn{x - y}, where \eqn{x > y > 0}.
#'
#' @param x An \linkS4class{Expression}
#' @param y An \linkS4class{Expression}
#' @return The difference x - y with domain {x,y: x > y > 0}.
DiffPos <- function(x, y) {
  x * OneMinusPos(y/x)
}

