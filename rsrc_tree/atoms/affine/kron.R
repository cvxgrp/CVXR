## CVXPY SOURCE: cvxpy/atoms/affine/kron.py
#'
#' The Kron class.
#'
#' This class represents the kronecker product.
#'
#' @slot lh_exp An \linkS4class{Expression} or numeric constant representing the left-hand matrix.
#' @slot rh_exp An \linkS4class{Expression} or numeric constant representing the right-hand matrix.
#' @name Kron-class
#' @aliases Kron
#' @rdname Kron-class
.Kron <- setClass("Kron", representation(lh_exp = "ConstValORExpr", rh_exp = "ConstValORExpr"), contains = "AffAtom")

#' @param lh_exp An \linkS4class{Expression} or numeric constant representing the left-hand matrix.
#' @param rh_exp An \linkS4class{Expression} or numeric constant representing the right-hand matrix.
#' @rdname Kron-class
Kron <- function(lh_exp, rh_exp) { .Kron(lh_exp = lh_exp, rh_exp = rh_exp) }

setMethod("initialize", "Kron", function(.Object, ..., lh_exp, rh_exp) {
  .Object@lh_exp <- lh_exp
  .Object@rh_exp <- rh_exp
  callNextMethod(.Object, ..., atom_args = list(.Object@lh_exp, .Object@rh_exp))
})

#' @param object A \linkS4class{Kron} object.
#' @param values A list of arguments to the atom.
#' @describeIn Kron The kronecker product of the two values.
setMethod("to_numeric", "Kron", function(object, values) {
  base::kronecker(values[[1]], values[[2]])
})

#' @describeIn Kron Check both arguments are vectors and the first is a constant.
setMethod("validate_args", "Kron", function(object) {
  if(!(is_constant(object@args[[1]]) || is_constant(object@args[[2]])))
    stop("The first argument to Kron must be constant.")
  else if(ndim(object@args[[1]]) != 2 || ndim(object@args[[2]]) != 2)
    stop("Kron requires matrix arguments.")
})

#' @describeIn Kron The dimensions of the atom.
setMethod("dim_from_args", "Kron", function(object) {
  rows <- dim(object@args[[1]])[1] * dim(object@args[[2]])[1]
  cols <- dim(object@args[[1]])[2] * dim(object@args[[2]])[2]
  return(c(rows, cols))
})

#' @describeIn Kron Is the atom convex?
setMethod("is_atom_convex", "Kron", function(object) {
  if(dpp_scope_active()) {
    # Kron is not DPP if any parameters are present
    x <- object@args[[1]]
    y <- object@args[[2]]
    return((is_constant(x) || is_constant(y)) && (is_param_free(x) && is_param_free(y)))
  } else
    return(is_constant(object@args[[1]]) || is_constant(object@args[[2]]))
})

#' @describeIn Kron Is the atom concave?
setMethod("is_atom_concave", "Kron", function(object) {
  return(is_atom_convex(object))
})

#' @describeIn Kron The sign of the atom.
setMethod("sign_from_args", "Kron", function(object) { mul_sign(object@args[[1]], object@args[[2]]) })

#' @param idx An index into the atom.
#' @describeIn Kron Is the composition non-decreasing in argument \code{idx}?
setMethod("is_incr", "Kron", function(object, idx) {
  cst_loc <- ifelse(is_constant(object@args[[1]]), 1, 2)
  is_nonneg(object@args[[cst_loc]])
})

#' @describeIn Kron Is the composition non-increasing in argument \code{idx}?
setMethod("is_decr", "Kron", function(object, idx) {
  cst_loc <- ifelse(is_constant(object@args[[1]]), 1, 2)
  is_nonpos(object@args[[1]])
})

#' @describeIn Kron Is the atom a positive semidefinite matrix?
setMethod("is_psd", "Kron", function(object) {
  # Check a *sufficient condition* that the expression is PSD, 
  # by checking if both arguments are PSD or both are NSD.
  case1 <- is_psd(object@args[[1]]) && is_psd(object@args[[2]])
  case2 <- is_nsd(object@args[[1]]) && is_nsd(object@args[[2]])
  return(case1 || case2)
})

#' @describeIn Kron Is the atom a negative semidefinite matrix?
setMethod("is_nsd", "Kron", function(object) {
  # Check a *sufficient condition* that the expression is NSD, 
  # by checking if one argument is PSD and the other is NSD.
  case1 <- is_psd(object@args[[1]]) && is_nsd(object@args[[2]])
  case2 <- is_nsd(object@args[[1]]) && is_psd(object@args[[2]])
  return(case1 || case2)
})

Kron.graph_implementation <- function(arg_objs, dim, data = NA_real_) {
  if(is_constant(object@args[[1]]))
    list(lo.kron_r(arg_objs[[1]], arg_objs[[2]], dim), list())
  else
    list(lo.kron_l(arg_objs[[1]], arg_objs[[2]], dim), list())
}

#' @param arg_objs A list of linear expressions for each argument.
#' @param dim A vector with two elements representing the size of the resulting expression.
#' @param data A list of additional data required by the atom.
#' @describeIn Kron The graph implementation of the atom.
setMethod("graph_implementation", "Kron", function(object, arg_objs, dim, data = NA_real_) {
  Kron.graph_implementation(arg_objs, dim, data)
})
