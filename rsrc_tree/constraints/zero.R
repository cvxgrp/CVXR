## CVXPY SOURCE: cvxpy/expression/consraints/zero.py
#'
#' The ZeroConstraint class
#'
#' @rdname ZeroConstraint-class
.ZeroConstraint <- setClass("ZeroConstraint", representation(expr = "Expression"), contains = "Constraint")
ZeroConstraint <- function(expr, constr_id = NA_integer_) { .ZeroConstraint(expr = expr, constr_id = constr_id) }

setMethod("initialize", "ZeroConstraint", function(.Object, ..., expr) {
  .Object@expr <- expr
  callNextMethod(.Object, ..., args = list(expr))
})

#' @param x,object A \linkS4class{ZeroConstraint} object.
#' @describeIn ZeroConstraint The string representation of the constraint.
setMethod("name", "ZeroConstraint", function(x) {
  # paste(as.character(x@args[[1]]), "== 0")
  paste(name(x@args[[1]]), "== 0")
})

setMethod("show", "ZeroConstraint", function(object) {
  print(paste("ZeroConstraint(", as.character(object@args[[1]]), ")", sep = ""))
})

#' @rdname ZeroConstraint-class
setMethod("as.character", "ZeroConstraint", function(x) {
  paste(as.character(x@args[[1]]), "== 0")
})

#' @describeIn ZeroConstraint The dimensions of the constrained expression.
setMethod("dim", "ZeroConstraint", function(x) { dim(x@args[[1]]) })

#' @describeIn Constraint The size of the constrained expression.
setMethod("size", "ZeroConstraint", function(object) { size(object@args[[1]]) })

#' @describeIn ZeroConstraint Is the constraint DCP?
setMethod("is_dcp", "ZeroConstraint", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping.
    return(is_affine(object@args[[1]]))
  }
  return(is_affine(object@args[[1]]))
})

#' @describeIn ZeroConstraint Is the constraint DGP?
setMethod("is_dgp", "ZeroConstraint", function(object, dpp = FALSE) { FALSE })

#' @describeIn ZeroConstraint Is the constraint DQCP?
setMethod("is_dqcp", "ZeroConstraint", function(object) { is_dcp(object) })

#' @describeIn ZeroConstraint The residual of a constraint
setMethod("residual", "ZeroConstraint", function(object) {
  val <- value(expr(object))
  if(any(is.na(val)))
    return(NA_real_)
  return(abs(val))
})

#' @describeIn ZeroConstraint The dual values of a constraint.
setMethod("dual_value", "ZeroConstraint", function(object) {
  return(value(object@dual_variables[[1]]))
})

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn ZeroConstraint Replaces the dual values of a constraint.
setReplaceMethod("dual_value", "ZeroConstraint", function(object, value) {
  value(object@dual_variables[[1]]) <- value
  return(object)
})

#'
#' The EqConstraint class
#'
#' @rdname EqConstraint-class
.EqConstraint <- setClass("EqConstraint", representation(lhs = "ConstValORExpr", rhs = "ConstValORExpr", expr = "ConstValORExpr"), prototype(expr = NA_real_), contains = "Constraint")
EqConstraint <- function(lhs, rhs, constr_id = NA_integer_) { .EqConstraint(lhs = lhs, rhs = rhs, constr_id = constr_id) }

setMethod("initialize", "EqConstraint", function(.Object, ..., lhs, rhs, expr = NA_real_) {
  .Object@lhs <- lhs
  .Object@rhs <- rhs
  .Object@expr <- lhs - rhs
  callNextMethod(.Object, ..., args = list(lhs, rhs))
})

setMethod(".construct_dual_variables", "EqConstraint", function(object, args) {
  callNextMethod(object, list(object@expr))
})

#' @param x,object A \linkS4class{EqConstraint} object.
#' @describeIn EqConstraint The string representation of the constraint.
setMethod("name", "EqConstraint", function(x) {
  # paste(as.character(x@args[[1]]), "==", as.character(x@args[[2]]))
  paste(name(x@args[[1]]), "==", name(x@args[[2]]))
})

setMethod("show", "EqConstraint", function(object) {
  print(paste("EqConstraint(", as.character(object@args[[1]]), ", ", as.character(object@args[[2]]), ")", sep = ""))
})

#' @rdname EqConstraint-class
setMethod("as.character", "EqConstraint", function(x) {
  paste(as.character(x@args[[1]]), "==", as.character(x@args[[2]]))
})

#' @describeIn EqConstraint The dimensions of the constrained expression.
setMethod("dim", "EqConstraint", function(x) { dim(x@expr) })

#' @describeIn EqConstraint The size of the constrained expression.
setMethod("size", "EqConstraint", function(object) { size(object@expr) })

#' @describeIn EqConstraint The expression to constrain.
setMethod("expr", "EqConstraint", function(object) { object@expr })

#' @describeIn EqConstraint Is the constraint DCP?
setMethod("is_dcp", "EqConstraint", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()
    return(is_affine(object@expr))
  }
  return(is_affine(object@expr))
})

#' @describeIn EqConstraint Is the constraint DGP?
setMethod("is_dgp", "EqConstraint", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()
    return(is_log_log_affine(object@args[[1]]) && is_log_log_affine(object@args[[2]]))
  }
  return(is_log_log_affine(object@args[[1]]) && is_log_log_affine(object@args[[2]]))
})

setMethod("is_dqcp", "EqConstraint", function(object) {
  return(is_dcp(object))
})

#' @describeIn EqConstraint The residual of the constraint..
setMethod("residual", "EqConstraint", function(object) {
  val <- value(object@expr)
  if(any(is.na(val)))
    return(NA_real_)
  return(abs(val))
})

#' @describeIn EqConstraint The dual values of a constraint.
setMethod("dual_value", "EqConstraint", function(object) {
  return(value(object@dual_variables[[1]]))
})

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn EqConstraint Replaces the dual values of a constraint.
setReplaceMethod("dual_value", "EqConstraint", function(object, value) {
  value(object@dual_variables[[1]]) <- value
  return(object)
})
