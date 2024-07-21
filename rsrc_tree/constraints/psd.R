## CVXPY SOURCE: cvxpy/expression/consraints/psd.py
#'
#' The PSDConstraint class.
#'
#' This class represents the positive semidefinite constraint, \eqn{\frac{1}{2}(X + X^T) \succeq 0}, i.e. \eqn{z^T(X + X^T)z \geq 0} for all \eqn{z}.
#'
#' @slot expr An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{X}.
#' @name PSDConstraint-class
#' @aliases PSDConstraint
#' @rdname PSDConstraint-class
.PSDConstraint <- setClass("PSDConstraint", representation(expr = "ConstValORExpr"),
                           validity = function(object) {
                             expr_dim <- dim(object@expr)
                             if(length(expr_dim) != 2 || expr_dim[1] != expr_dim[2])
                               stop("Non-square matrix in positive definite constraint.")
                             return(TRUE)
                           }, contains = "Constraint")

#' @param expr An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{X}.
#' @param constr_id (Optional) A numeric value representing the constraint ID.
#' @rdname PSDConstraint-class
PSDConstraint <- function(expr, constr_id = NA_integer_) { .PSDConstraint(expr = expr, constr_id = constr_id) }

setMethod("initialize", "PSDConstraint", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., args = list(expr))
})

#' @param x,object A \linkS4class{PSDConstraint} object.
#' @describeIn PSDConstraint The string representation of the constraint.
setMethod("name", "PSDConstraint", function(x) {
  # paste(as.character(x@args[[1]]), ">> 0")
  paste(name(x@args[[1]]), ">> 0")
})

#' @describeIn PSDConstraint The constraint is DCP if the constrained expression is affine.
setMethod("is_dcp", "PSDConstraint", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping
    return(is_affine(object@args[[1]]))
  }
  return(is_affine(object@args[[1]]))
})

#' @describeIn PSDConstraint Is the constraint DGP?
setMethod("is_dgp", "PSDConstraint", function(object, dpp = FALSE) { FALSE })

#' @describeIn PSDConstraint Is the constraint DQCP?
setMethod("is_dqcp", "PSDConstraint", function(object) { is_dcp(object) })

#' @describeIn PSDConstraint A \linkS4class{Expression} representing the residual of the constraint.
setMethod("residual", "PSDConstraint", function(object) {
  val <- value(expr(object))
  if(any(is.na(val)))
    return(NA_real_)
  min_eig <- LambdaMin(object@args[[1]] + t(object@args[[1]]))/2
  return(value(Neg(min_eig)))
})
