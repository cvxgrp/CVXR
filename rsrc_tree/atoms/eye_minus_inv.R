## CVXPY SOURCE: cvxpy/atoms/eye_minus_inv.py

#'
#' The EyeMinusInv class.
#'
#' This class represents the unity resolvent of an elementwise positive matrix \eqn{X}, i.e., \eqn{(I - X)^{-1}},
#' and it enforces the constraint that the spectral radius of \eqn{X} is at most 1.
#' This atom is log-log convex.
#'
#' @slot X An \linkS4class{Expression} or numeric matrix.
#' @name EyeMinusInv-class
#' @aliases EyeMinusInv
#' @rdname EyeMinusInv-class
.EyeMinusInv <- setClass("EyeMinusInv", representation(X = "ConstValORExpr"),
                         validity = function(object) {
                           if(length(dim(object@X)) != 2 || nrow(object@X) != ncol(object@X))
                             stop("[EyeMinusInv: X] The argument X must be a square matrix.")
                           return(TRUE)
                          }, contains = "Atom")

#' @param X An \linkS4class{Expression} or numeric matrix.
#' @rdname EyeMinusInv-class
EyeMinusInv <- function(X) { .EyeMinusInv(X = X) }

setMethod("initialize", "EyeMinusInv", function(.Object, ..., X) {
  .Object@X <- X
  .Object <- callNextMethod(.Object, ..., atom_args = list(.Object@X))
  .Object@args[[1]] <- X
  .Object
})

#' @param object,x An \linkS4class{EyeMinusInv} object.
#' @param values A list of arguments to the atom.
#' @describeIn EyeMinusInv The unity resolvent of the matrix.
setMethod("to_numeric", "EyeMinusInv", function(object, values) {
  base::solve(diag(nrow(object@args[[1]])) - values[[1]])
})

#' @describeIn EyeMinusInv The name and arguments of the atom.
setMethod("name", "EyeMinusInv", function(x) { paste(class(x), x@args[[1]]) })

#' @describeIn EyeMinusInv The dimensions of the atom determined from its arguments.
setMethod("dim_from_args", "EyeMinusInv", function(object) { dim(object@args[[1]]) })

#' @describeIn EyeMinusInv The (is positive, is negative) sign of the atom.
setMethod("sign_from_args", "EyeMinusInv", function(object) { c(TRUE, FALSE) })

#' @describeIn EyeMinusInv Is the atom convex?
setMethod("is_atom_convex", "EyeMinusInv", function(object) { FALSE })

#' @describeIn EyeMinusInv Is the atom concave?
setMethod("is_atom_concave", "EyeMinusInv", function(object) { FALSE })

#' @describeIn EyeMinusInv Is the atom log-log convex?
setMethod("is_atom_log_log_convex", "EyeMinusInv", function(object) { TRUE })

#' @describeIn EyeMinusInv Is the atom log-log concave?
setMethod("is_atom_log_log_concave", "EyeMinusInv", function(object) { FALSE })

# TODO: Figure out monotonicity.
#' @param idx An index into the atom.
#' @describeIn EyeMinusInv Is the atom weakly increasing in the index?
setMethod("is_incr", "EyeMinusInv", function(object, idx) { FALSE })

#' @describeIn EyeMinusInv Is the atom weakly decreasing in the index?
setMethod("is_decr", "EyeMinusInv", function(object, idx) { FALSE })

#' @param values A list of numeric values for the arguments
#' @describeIn EyeMinusInv Gives the (sub/super)gradient of the atom w.r.t. each variable
setMethod(".grad", "EyeMinusInv", function(object, values) { NA_real_ })

# The resolvent of a positive matrix, (sI - X)^(-1).
# For an elementwise positive matrix X and a positive scalar s, this atom computes
# (sI - X)^(-1), and it enforces the constraint that the spectral radius of X/s is
# at most 1.
# This atom is log-log convex.
Resolvent <- function(X, s) {
  1.0 / (s * EyeMinusInv(X / s))
}

