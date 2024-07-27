## CVXPY SOURCE: cvxpy/expression/constraints/nonpos.py

#'
#' The NonPosConstraint class
#'
#' A constraint of the form \eqn{x \leq 0}.
#'
#' The preferred way of creating a NonPosConstraint constraint is through
#' operator overloading. To constrain an expression x to be non-positive,
#' simply write \code{x <= 0}; to constrain x to be non-negative, write
#' \code{x >= 0}. The former creates a NonPosConstraint constraint with x
#' as its argument, while the latter creates one with -x as its argument.
#' Strict inequalities are not supported, as they do not make sense in a
#' numerical setting.
#'
#' @rdname NonPosConstraint-class
.NonPosConstraint <- setClass("NonPosConstraint", representation(expr = "Expression"), contains = "Constraint")
NonPosConstraint <- function(expr, constr_id = NA_integer_) { .NonPosConstraint(expr = expr, constr_id = constr_id) }

setMethod("initialize", "NonPosConstraint", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., args = list(expr))
  if(!is_real(.Object@args[[1]]))
    stop("Input to NonPosConstraint must be real")
})

#' @param x,object A \linkS4class{NonPosConstraint} object.
#' @describeIn NonPosConstraint The string representation of the constraint.
setMethod("name", "NonPosConstraint", function(x) {
  # paste(as.character(x@args[[1]]), "<= 0")
  paste(name(x@args[[1]]), "<= 0")
})

#' @describeIn NonPosConstraint A non-positive constraint is DCP if its argument is convex.
setMethod("is_dcp", "NonPosConstraint", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping.
    return(is_convex(object@args[[1]]))
  }
  return(is_convex(object@args[[1]]))
})

#' @describeIn NonPosConstraint Is the constraint DGP?
setMethod("is_dgp", "NonPosConstraint", function(object, dpp = FALSE) { FALSE })

#' @describeIn NonPosConstraint Is the constraint DQCP?
setMethod("is_dqcp", "NonPosConstraint", function(object) { is_quasiconvex(object@args[[1]]) })

#' @describeIn NonPosConstraint The residual of the constraint.
setMethod("residual", "NonPosConstraint", function(object) {
  val <- value(expr(object))
  if(any(is.na(val)))
    return(NA_real_)
  return(pmax(val, 0))
})

#' @describeIn NonPosConstraint The violation of the constraint.
setMethod("violation", "NonPosConstraint", function(object) {
  resid <- residual(object)
  if(any(is.na(resid)))
    stop("Cannot compute the violation of a constraint whose expression is NA-valued.")
  viol <- base::norm(resid, type = "2")
  return(viol)
})

#'
#' The NonNegConstraint class
#'
#' A constraint of the form \eqn{x \geq 0}.
#'
#' This class was created to account for the fact that the ConicSolver interface
#' returns matrix data stated with respect to the nonnegative orthant, rather than
#' the nonpositive orthant. This class can be removed if the behavior of ConicSolver is
#' changed. However the current behavior of ConicSolver means CVXR's dual variable
#' and Lagrangian convention follows the most common convention in the literature.
#'
#' @rdname NonNegConstraint-class
.NonNegConstraint <- setClass("NonNegConstraint", representation(expr = "Expression"), contains = "Constraint")
NonNegConstraint <- function(expr, constr_id = NA_integer_) { .NonNegConstraint(expr = expr, constr_id = constr_id) }

setMethod("initialize", "NonNegConstraint", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., args = list(expr))
  if(!is_real(.Object@args[[1]]))
    stop("Input to NonNegConstraint must be real")
})

#' @param x,object A \linkS4class{NonNegConstraint} object.
#' @describeIn NonNegConstraint The string representation of the constraint.
setMethod("name", "NonNegConstraint", function(x) {
  # paste("0 <=", as.character(x@args[[1]]))
  paste("0 <=", name(x@args[[1]]))
})

#' @describeIn NonNegConstraint A non-positive constraint is DCP if its argument is convex.
setMethod("is_dcp", "NonNegConstraint", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping.
    return(is_concave(object@args[[1]]))
  }
  return(is_concave(object@args[[1]]))
})

#' @describeIn NonNegConstraint Is the constraint DGP?
setMethod("is_dgp", "NonNegConstraint", function(object, dpp = FALSE) { FALSE })

#' @describeIn NonNegConstraint Is the constraint DQCP?
setMethod("is_dqcp", "NonNegConstraint", function(object) { is_quasiconcave(object@args[[1]]) })

#' @describeIn NonNegConstraint The residual of the constraint.
setMethod("residual", "NonNegConstraint", function(object) {
  val <- value(expr(object))
  if(any(is.na(val)))
    return(NA_real_)
  return(pmin(val, 0))
})

#' @describeIn NonNegConstraint The violation of the constraint.
setMethod("violation", "NonNegConstraint", function(object) {
  resid <- residual(object)
  if(any(is.na(resid)))
    stop("Cannot compute the violation of a constraint whose expression is NA-valued.")
  viol <- base::norm(resid, type = "2")
  return(viol)
})

#'
#' The IneqConstraint class
#'
#' A constraint of the form \eqn{x \leq y}.
#'
#' @slot lhs The expression to be upper-bounded by rhs.
#' @slot rhs The expression to be lower-bounded by lhs.
#' @rdname IneqConstraint-class
.IneqConstraint <- setClass("IneqConstraint", representation(lhs = "ConstValORExpr", rhs = "ConstValORExpr", expr = "ConstValORExpr"),
                            prototype(expr = NA_real_), contains = "Constraint")

IneqConstraint <- function(lhs, rhs, constr_id = NA_integer_) { .IneqConstraint(lhs = lhs, rhs = rhs, constr_id = constr_id) }

setMethod("initialize", "IneqConstraint", function(.Object, ..., lhs, rhs, expr = NA_real_) {
  .Object@lhs <- lhs
  .Object@rhs <- rhs
  .Object@expr <- lhs - rhs
  if(is_complex(.Object@expr))   # TODO: Remove this restriction when implemented.
    stop("Inequality constraints cannot be complex.")
  callNextMethod(.Object, ..., args = list(lhs, rhs))
})

#' @param x,object A \linkS4class{IneqConstraint} object.
#' @describeIn IneqConstraint The string representation of the constraint.
setMethod("name", "IneqConstraint", function(x) {
  # paste(as.character(x@args[[1]]), "<=", as.character(x@args[[2]]))
  paste(name(x@args[[1]]), "<=", name(x@args[[2]]))
})

#' @describeIn IneqConstraint The dimensions of the constrained expression.
setMethod("dim", "IneqConstraint", function(x) { dim(x@expr) })

#' @describeIn IneqConstraint The size of the constrained expression.
setMethod("size", "IneqConstraint", function(object) { size(object@expr) })

#' @describeIn IneqConstraint The expression to constrain.
setMethod("expr", "IneqConstraint", function(object) { object@expr })

#' @describeIn IneqConstraint A non-positive constraint is DCP if its argument is convex.
setMethod("is_dcp", "IneqConstraint", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping.
    return(is_convex(object@expr))
  }
  return(is_convex(object@expr))
})

#' @describeIn IneqConstraint Is the constraint DGP?
setMethod("is_dgp", "IneqConstraint", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping.
    return(is_log_log_convex(object@args[[1]]) && is_log_log_concave(object@args[[2]]))
  }
  return(is_log_log_convex(object@args[[1]]) && is_log_log_concave(object@args[[2]]))
})

#' @param context Must be either 'dcp' (disciplined convex program) or 'dgp' (disciplined geometric program).
#' @describeIn Constraint is the constraint DPP in the given context?
setMethod("is_dpp", "IneqConstraint", function(object, context = "dcp") {
  if(tolower(context) == "dcp")
    return(is_dcp(object, dpp = TRUE))
  else if(tolower(context) == "dgp")
    return(is_dgp(object, dpp = TRUE))
  else
    stop("Unsupported context ", context)
})

#' @describeIn IneqConstraint Is the constraing DQCP?
setMethod("is_dqcp", "IneqConstraint", function(object) {
  is_dcp(object) ||
  (is_quasiconvex(object@args[[1]]) && is_constant(object@args[[2]])) ||
  (is_constant(object@args[[1]]) && is_quasiconcave(object@args[[2]]))
})

#' @describeIn IneqConstraint The residual of the constraint.
setMethod("residual", "IneqConstraint", function(object) {
  val <- value(object@expr)
  if(any(is.na(val)))
    return(NA_real_)
  return(pmax(val, 0))
})
