#'
#' The Constraint class.
#'
#' This virtual class represents a mathematical constraint.
#'
#' @name Constraint-class
#' @aliases Constraint
#' @rdname Constraint-class
setClass("Constraint", representation(dual_variables = "list", constr_id = "integer"), prototype(dual_variables = list(), constr_id = NA_integer_), contains = "Canonical")

setMethod("initialize", "Constraint", function(.Object, ..., dual_variables = list(), constr_id = NA_integer_) {
  .Object <- callNextMethod(.Object, ...)
  # TODO: Cast constants
  # .Object@args <- lapply(args, as.Constant)
  .Object@constr_id <- ifelse(is.na(constr_id), get_id(), constr_id)
  # .Object@dual_variables <- lapply(.Object@args, function(arg) { Variable(dim(arg)) })
  .Object@dual_variables <- lapply(.Object@args, function(arg) { new("Variable", dim = dim(arg)) })
  return(.Object)
})

#' @param x,object A \linkS4class{Constraint} object.
#' @rdname Constraint-class
setMethod("as.character", "Constraint", function(x) { name(x) })

setMethod("show", "Constraint", function(object) {
  print(paste(class(object), "(", as.character(object@args[[1]]), ")", sep = ""))
})

#' @describeIn Constraint The dimensions of the constrained expression.
setMethod("dim", "Constraint", function(x) { dim(x@args[[1]]) })

#' @describeIn Constraint The size of the constrained expression.
setMethod("size", "Constraint", function(object) { size(object@args[[1]]) })

#' @describeIn Constraint Is the constraint real?
setMethod("is_real", "Constraint", function(object) { !is_complex(object) })

#' @describeIn Constraint Is the constraint imaginary?
setMethod("is_imag", "Constraint", function(object) { all(sapply(object@args, is_imag)) })

#' @describeIn Constraint Is the constraint complex?
setMethod("is_complex", "Constraint", function(object) { any(sapply(object@args, is_complex)) })

#' @param dpp A logical value indicating whether we are solving a disciplined parameterized program (DPP).
#' @describeIn Constraint Is the constraint DCP?
setMethod("is_dcp", "Constraint", function(object, dpp = FALSE) { stop("Unimplemented") })

#' @describeIn Constraint Is the constraint DGP?
setMethod("is_dgp", "Constraint", function(object, dpp = FALSE) { stop("Unimplemented") })

#' @param context Must be either 'dcp' (disciplined convex program) or 'dgp' (disciplined geometric program).
#' @describeIn Constraint is the constraint DPP in the given context?
setMethod("is_dpp", "Constraint", function(object, context = "dcp") {
  if(tolower(context) == "dcp")
    return(is_dcp(object, dpp = TRUE))
  else if(tolower(context) == "dgp")
    return(is_dpp(object, dpp = TRUE))
  else
    stop("Unsupported context ", context)
})

#' @describeIn Constraint The residual of a constraint
setMethod("residual", "Constraint", function(object) { stop("Unimplemented") })

#' @describeIn Constraint The numeric residual of a constraint. The violation is defined as the
#' distance between the constrained expression's value and its projection onto the domain of the
#' constraint: ||\\Pi(v) - v||_2^2, where `v` is the value of the constrained expression and
#' `\\Pi` is the projection operator onto the constraint's domain.
setMethod("violation", "Constraint", function(object) {
  resid <- residual(object)
  if(any(is.na(resid)))
    stop("Cannot compute the violation of a constraint whose expression is NA-valued.")
  return(resid)
})

#' @param tolerance The tolerance for checking if the constraint is violated.
#' @describeIn Constraint Checks whether the constraint violation is less than a tolerance.
setMethod("constr_value", "Constraint", function(object, tolerance = 1e-8) {
  resid <- residual(object)
  if(any(is.na(resid)))
    stop("Cannot compute the value of a constraint whose expression is NA-valued.")
  return(all(resid <= tolerance))
})


#'
#' A Class Union of List and Constraint
#'
#' @name ListORConstr-class
#' @rdname ListORConstr-class
setClassUnion("ListORConstr", c("list", "Constraint"))

# Helper function since syntax is different for LinOp (list) vs. Constraint object
#' @param object A list or \linkS4class{Constraint} object.
#' @describeIn ListORConstr Returns the ID associated with the list or constraint.
setMethod("id", "ListORConstr", function(object) {
  if(is.list(object))
    object$constr_id
  else
    object@constr_id
})

#' @describeIn Constraint Information needed to reconstruct the object aside from the args.
setMethod("get_data", "Constraint", function(object) { list(id(object)) })

#' @describeIn Constraint The dual values of a constraint.
setMethod("dual_value", "Constraint", function(object) {
  dual_vals <- lapply(object@dual_variables, value)
  if(length(dual_vals) == 1)
    return(dual_vals[[1]])
  else
    return(dual_vals)
})

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn Constraint Replaces the dual values of a constraint.
setReplaceMethod("dual_value", "Constraint", function(object, value) {
  value(object@dual_variables[[1]]) <- value
  return(object)
})

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
    return(is_dpp(object, dpp = TRUE))
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

#'
#' The FiniteSet class.
#'
#' This class represents a constraint that each entry of an Expression to take a value in a given set of real numbers.
#'
#' @slot expre The given expression to be constrained. This Expression must be affine.
#' If expre has multiple elements, then the constraint is applied separately to
#' each element, i.e., after solving a problem with this constraint, we should have:
#' \code{for(e in flatten(expre)) { print(value(e) %in% vec) # => TRUE }
#' @slot vec The finite collection of values to which each entry of expre is to be constrained.
#' @slot ineq_form A logical value controlling how this constraint is canonicalized into mixed-integer linear constraints.
#' If TRUE, then we use a formulation with (size(vec) - 1) inequality constraints,
#' one equality constraint, and (size(vec) - 1) binary variables for each element
#' of expre. If FALSE, then we use a formulation with size(vec) binary variables and two
#' equality constraints for each element of expre. Defaults to FALSE. The case ineq_form = TRUE may speed up some mixed-integer
#' solvers that use simple branch and bound methods.
#' @name FiniteSet-class
#' @aliases FiniteSet
#' @rdname FiniteSet-class
.FiniteSet <- setClass("FiniteSet", representation(expre = "ConstValORExpr", vec = "list", ineq_form = "logical"),
                       prototype(ineq_form = FALSE), contains = "Constraint")

#' @param expre An affine Expression object.
#' @param vec The finite collection of values to which each entry of expre is to be constrained.
#' @param ineq_form A logical value controlling how this constraint is canonicalized.
#' @param constr_id (Optional) An integer representing the unique ID of the constraint.
#' @rdname FiniteSet-class
FiniteSet <- function(expre, vec, ineq_form = FALSE, constr_id = NA_integer_) { .FiniteSet(expre = expre, vec = vec, ineq_form = ineq_form, constr_id = constr_id) }

setMethod("initialize", "FiniteSet", function(.Object, ..., expre, vec, ineq_form = FALSE) {
  vec <- flatten(as.Constant(vec))
  if(!is_affine(expre))
    stop("Provided Expression must be affine, but had curvature ", curvature(expre))

  # Note: we use the term "expre" rather than "expr" since
  # "expr" is already a property used by all Constraint classes.
  .Object@expre <- expre
  .Object@vec <- vec
  .Object@ineq_form <- ineq_form
  callNextMethod(.Object, ..., args = list(expre, vec))
})

#' @param x,object A \linkS4class{FiniteSet} object.
#' @describeIn FiniteSet The string representation of the constraint.
setMethod("name", "FiniteSet", function(x) {
  paste("FiniteSet(", as.character(x@args[[1]]), ", ", as.character(x@args[[2]]), ")", sep = "")
})

#' @describeIn FiniteSet Information needed to reconstruct the object aside from the args.
setMethod("get_data", "FiniteSet", function(object) { list(object@ineq_form, id(object)) })

#' @describeIn FiniteSet The constraint is DCP if the constrained expression is affine.
setMethod("is_dcp", "FiniteSet", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping
    return(is_affine(object@args[[1]]))
  }
  return(is_affine(object@args[[1]]))
})

#' @describeIn FiniteSet Is the constraint DGP?
setMethod("is_dgp", "FiniteSet", function(object, dpp = FALSE) { FALSE })

#' @describeIn FiniteSet Is the constraint DQCP?
setMethod("is_dqcp", "FiniteSet", function(object) { is_dcp(object) })

#' @describeIn FiniteSet The size of the constrained expression.
setMethod("size", "FiniteSet", function(object) { size(object@expre) })

setMethod("dim", "FiniteSet", function(x) { dim(object@expre) })

#' @describeIn FiniteSet The residual of the constraint.
setMethod("residual", "FiniteSet", function(object) {
  expr_val <- as.vector(value(object@expre))
  vec_val <- value(object@vec)
  resids <- sapply(expr_val, min(abs(val - vec_val)))
  res <- max(resids)
  return(res)
})

#'
#' The ExpCone class.
#'
#' This class represents a reformulated exponential cone constraint operating elementwise on \eqn{a, b, c}.
#'
#' Original cone:
#' \deqn{
#' K = \{(x,y,z) | y > 0, ye^{x/y} \leq z\} \cup \{(x,y,z) | x \leq 0, y = 0, z \geq 0\}
#' }
#' Reformulated cone:
#' \deqn{
#' K = \{(x,y,z) | y, z > 0, y\log(y) + x \leq y\log(z)\} \cup \{(x,y,z) | x \leq 0, y = 0, z \geq 0\}
#' }
#'
#' @slot x The variable \eqn{x} in the exponential cone.
#' @slot y The variable \eqn{y} in the exponential cone.
#' @slot z The variable \eqn{z} in the exponential cone.
#' @name ExpCone-class
#' @aliases ExpCone
#' @rdname ExpCone-class
.ExpCone <- setClass("ExpCone", representation(x = "ConstValORExpr", y = "ConstValORExpr", z = "ConstValORExpr"), contains = "Constraint")

#' @param x The variable \eqn{x} in the exponential cone.
#' @param y The variable \eqn{y} in the exponential cone.
#' @param z The variable \eqn{z} in the exponential cone.
#' @param constr_id (Optional) A numeric value representing the constraint ID.
#' @rdname ExpCone-class
## #' @export
ExpCone <- function(x, y, z, constr_id = NA_integer_) { .ExpCone(x = x, y = y, z = z, constr_id = constr_id) }

setMethod("initialize", "ExpCone", function(.Object, ..., x, y, z) {
  .Object@x <- as.Constant(x)
  .Object@y <- as.Constant(y)
  .Object@z <- as.Constant(z)
  args <- list(.Object@x, .Object@y, .Object@z)
  for(val in args) {
    if(!(is_affine(val) && is_real(val)))
      stop("All arguments must be affine and real")
  }
  xs <- dim(.Object@x)
  ys <- dim(.Object@y)
  zs <- dim(.Object@z)
  if(!all(xs == ys) || !all(xs == zs))
    stop("All arguments must have the same shapes. Provided arguments have shapes ", xs, ", ", ys, ", and ", zs)
  callNextMethod(.Object, ..., args = list(.Object@x, .Object@y, .Object@z))
})

setMethod("show", "ExpCone", function(object) {
  print(paste("ExpCone(", as.character(object@x), ", ", as.character(object@y), ", ", as.character(object@z), ")", sep = ""))
})

#' @rdname ExpCone-class
setMethod("as.character", "ExpCone", function(x) {
  paste("ExpCone(", as.character(x@x), ", ", as.character(x@y), ", ", as.character(x@z), ")", sep = "")
})

#' @param object A \linkS4class{ExpCone} object.
#' @describeIn ExpCone The size of the \code{x} argument.
setMethod("residual", "ExpCone", function(object) {
  # TODO: The projection should be implemented directly.
  if(any(is.na(value(object@x))) || any(is.na(value(object@y))) || any(is.na(value(object@z))))
    return(NA_real_)
  # x <- Variable(dim(object@x))
  # y <- Variable(dim(object@y))
  # z <- Variable(dim(object@z))
  x <- new("Variable", dim = dim(object@x))
  y <- new("Variable", dim = dim(object@y))
  z <- new("Variable", dim = dim(object@z))
  constr <- list(ExpCone(x, y, z))
  obj <- Minimize(Norm2(HStack(x, y, z) - HStack(value(object@x), value(object@y), value(object@z))))
  prob <- Problem(obj, constr)
  result <- solve(prob)
  return(result$value)
})

#' @describeIn ExpCone The number of entries in the combined cones.
setMethod("size", "ExpCone", function(object) { 3 * num_cones(object) })

#' @describeIn ExpCone The number of elementwise cones.
setMethod("num_cones", "ExpCone", function(object) { size(object@x) })

setMethod("as_quad_approx", "ExpCone", function(object, m, k) {
  if(as.integer(m) != m)
    stop("m must be an integer")
  if(as.integer(k) != k)
    stop("k must be an integer")
  return(RelEntrConeQuad(object@y, object@z, -object@x, m, k))
})

#' @describeIn ExpCone The dimensions of the exponential cones.
setMethod("cone_sizes", "ExpCone", function(object) { rep(3, num_cones(object)) })

#' @describeIn ExpCone An exponential constraint is DCP if each argument is affine.
setMethod("is_dcp", "ExpCone", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping
    return(all(sapply(object@args, is_affine)))
  }
  return(all(sapply(object@args, is_affine)))
})

#' @describeIn ExpCone Is the constraint DGP?
setMethod("is_dgp", "ExpCone", function(object, dpp = FALSE) { FALSE })

#' @describeIn ExpCone Is the constraint DQCP?
setMethod("is_dqcp", "ExpCone", function(object) { is_dcp(object) })

setMethod("dim", "ExpCone", function(x) { c(3, dim(x@x)) })

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn ExpCone Replaces the dual values of an exponential cone constraint.
setReplaceMethod("dual_value", "ExpCone", function(object, value) {
  # TODO: verify that reshaping below works correctly
  value <- t(matrix(t(value), nrow = 3, byrow = FALSE))
  dv0 <- matrix(value[, 1], nrow = nrow(object@x), ncol = ncol(object@x))
  dv1 <- matrix(value[, 2], nrow = nrow(object@y), ncol = ncol(object@y))
  dv2 <- matrix(value[, 3], nrow = nrow(object@z), ncol = ncol(object@z))

  value(object@dual_variables[[1]]) <- dv0
  value(object@dual_variables[[2]]) <- dv1
  value(object@dual_variables[[3]]) <- dv2
  return(object)
})

#'
#' The RelEntrConeQuad class.
#'
#' This class represents an approximate construction of the scalar relative entropy cone.
#'
#' \deqn{
#'  K_{re}=\\text{cl}\\{(x,y,z)\\in\\mathbb{R}_{++} \\times \\mathbb{R}_{++}\\times\\mathbb{R}_{++}\\:x\\log(x/y)\\leq z\\}
#' }
#'
#' Since the above definition is very similar to the ExpCone, we provide a conversion method.
#'
#' More details on the approximation can be found in Theorem-3 on page-10 in the paper:
#' Semidefinite Approximations of the Matrix Logarithm.
#'
#' @slot x The variable \eqn{x} in the (approximate) scalar relative entropy cone.
#' @slot y The variable \eqn{y} in the (approximate) scalar relative entropy cone.
#' @slot z The variable \eqn{z} in the (approximate) scalar relative entropy cone.
#' @slot m An integer directly related to the number of generated nodes for the quadrature approximation used in the algorithm.
#' @slot k An integer controlling the approximation.
#' @name RelEntrConeQuad-class
#' @aliases RelEntrConeQuad
#' @rdname RelEntrConeQuad-class
.RelEntrConeQuad <- setClass("RelEntrConeQuad", representation(x = "ConstValORExpr", y = "ConstValORExpr", z = "ConstValORExpr", m = "numeric", k = "numeric"),
                             validity = function(object) {
                               if(as.integer(object@m) != object@m)
                                 stop("[RelEntrConeQuad: m] The argument m must be an integer")
                               if(as.integer(object@k) != object@k)
                                 stop("[RelEntrConeQuad: k] The argument k must be an integer")
                               return(TRUE)
                              }, contains = "Constraint")

#' @param x The variable \eqn{x} in the (approximate) scalar relative entropy cone.
#' @param y The variable \eqn{y} in the (approximate) scalar relative entropy cone.
#' @param z The variable \eqn{z} in the (approximate) scalar relative entropy cone.
#' @param m An integer directly related to the number of generated nodes for the quadrature approximation used in the algorithm.
#' @param k An integer controlling the approximation.
#' @param constr_id (Optional) A numeric value representing the constraint ID.
#' @rdname RelEntrConeQuad-class
## #' @export
RelEntrConeQuad <- function(x, y, z, m, k, constr_id = NA_integer_) { .RelEntrConeQuad(x = x, y = y, z = z, m = m, k = k, constr_id = constr_id) }

setMethod("initialize", "RelEntrConeQuad", function(.Object, ..., x, y, z, m, k) {
  .Object@x <- as.Constant(x)
  .Object@y <- as.Constant(y)
  .Object@z <- as.Constant(z)
  args <- list(.Object@x, .Object@y, .Object@z)
  for(val in args) {
    if(!(is_affine(val) && is_real(val)))
      stop("All arguments must be affine and real")
  }
  .Object@m <- m
  .Object@k <- k
  xs <- dim(.Object@x)
  ys <- dim(.Object@y)
  zs <- dim(.Object@z)
  if(!all(xs == ys) || !all(xs == zs))
    stop("All arguments must have the same shapes. Provided arguments have shapes ", xs, ", ", ys, ", and ", zs)
  callNextMethod(.Object, ..., args = list(.Object@x, .Object@y, .Object@z))
})

#' @describeIn RelEntrConeQuad Information needed to reconstruct the object aside from the args.
setMethod("get_data", "RelEntrConeQuad", function(object) { list(object@m, object@k, id(object)) })

setMethod("show", "RelEntrConeQuad", function(object) {
  print(paste("RelEntrConeQuad(", as.character(object@x), ", ", as.character(object@y), ", ", as.character(object@z), ", ", object@m, ", ", object@k, ")", sep = ""))
})

#' @rdname RelEntrConeQuad-class
setMethod("as.character", "RelEntrConeQuad", function(x) {
  paste("RelEntrConeQuad(", as.character(x@x), ", ", as.character(x@y), ", ", as.character(x@z), ", ", x@m, ", ", x@k, ")", sep = "")
})

setMethod("residual", "RelEntrConeQuad", function(object) {
  # TODO: The projection should be implemented directly.
  if(is.na(value(object@x)) || is.na(value(object@y)) || is.na(value(object@z)))
    return(NA_real_)

  x <- new("Variable", dim = dim(object@x))
  y <- new("Variable", dim = dim(object@y))
  z <- new("Variable", dim = dim(object@z))
  constr <- list(RelEntrConeQuad(x, y, z, object@m, object@k))
  obj <- Minimize(Norm2(HStack(x, y, z) - HStack(value(object@x), value(object@y), value(object@z))))
  problem <- Problem(obj, constr)
  result <- solve(problem)
  return(result$value)
})

#' @describeIn RelEntrConeQuad The number of entries in the combined cones.
setMethod("size", "RelEntrConeQuad", function(object) { 3*num_cones(object) })

#' @describeIn RelEntrConeQuad The number of elementwise cones.
setMethod("num_cones", "RelEntrConeQuad", function(object) { size(object@x) })

#' @describeIn RelEntrConeQuad The dimensions of the exponential cones.
setMethod("cone_sizes", "RelEntrConeQuad", function(object) { rep(3, num_cones(object)) })

#' @describeIn RelEntrConeQuad An exponential constraint is DCP if each argument is affine.
setMethod("is_dcp", "RelEntrConeQuad", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping
    return(all(sapply(object@args, is_affine)))
  }
  return(all(sapply(object@args, is_affine)))
})

#' @describeIn RelEntrConeQuad Is the constraint DGP?
setMethod("is_dgp", "RelEntrConeQuad", function(object, dpp = FALSE) { FALSE })

#' @describeIn RelEntrConeQuad Is the constraint DQCP?
setMethod("is_dqcp", "RelEntrConeQuad", function(object) { is_dcp(object) })

setMethod("dim", "RelEntrConeQuad", function(x) { c(3, dim(x@x)) })

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn RelEntrConeQuad Replaces the dual values of an exponential cone constraint.
setReplaceMethod("dual_value", "RelEntrConeQuad", function(object, value) {
  stop("Unimplemented")
})

#'
#' The OpRelEntrConeQuad class.
#'
#' This class represents an approximate construction of the scalar relative entropy cone.
#'
#' \deqn{
#'  K_{re}^n=\\text{cl}\\{(X,Y,T)\\in\\mathbb{H}^n_{++} \\times \\mathbb{H}^n_{++}\\times\\mathbb{H}^n_{++}\\:D_{\\text{op}}\\succeq T\\}
#' }
#'
#' This approximation uses \eqn{m + k} semidefinite constraints.
#'
#' More details on the approximation can be found in Theorem-3 on page-10 in the paper:
#' Semidefinite Approximations of the Matrix Logarithm.
#'
#' @slot X The variable \eqn{X} in the (approximate) operator relative entropy cone.
#' @slot Y The variable \eqn{Y} in the (approximate) operator relative entropy cone.
#' @slot Z The variable \eqn{Z} in the (approximate) operator relative entropy cone.
#' @slot m A positive integer that controls the number of quadrature nodes used in a local approximation of the matrix logarithm. Increasing this value results in better local approximations, but does not significantly expand the region of inputs for which the approximation is effective.
#' @slot k A positive integer that sets the number of scaling points about which the quadrature approximation is performed. Increasing this value will expand the region of inputs over which the approximation is effective.
#' @name OpRelEntrConeQuad-class
#' @aliases OpRelEntrConeQuad
#' @rdname OpRelEntrConeQuad-class
.OpRelEntrConeQuad <- setClass("OpRelEntrConeQuad", representation(X = "ConstValORExpr", Y = "ConstValORExpr", Z = "ConstValORExpr", m = "numeric", k = "numeric"),
                               validity = function(object) {
                                 if(as.integer(object@m) != object@m)
                                   stop("[OpRelEntrConeQuad: m] The argument m must be an integer")
                                 if(object@m <= 0)
                                   stop("[OpRelEntrConeQuad: m] The argument m must be positive")
                                 if(as.integer(object@k) != object@k)
                                   stop("[OpRelEntrConeQuad: k] The argument k must be an integer")
                                 if(object@k <= 0)
                                   stop("[OpRelEntrConeQuad: k] The argument k must be positive")
                                 return(TRUE)
                               }, contains = "Constraint")

#' @param X The variable \eqn{X} in the (approximate) operator relative entropy cone.
#' @param Y The variable \eqn{y} in the (approximate) operator relative entropy cone.
#' @param Z The variable \eqn{z} in the (approximate) operator relative entropy cone.
#' @param m A positive integer that controls the number of quadrature nodes used in a local approximation of the matrix logarithm. Increasing this value results in better local approximations, but does not significantly expand the region of inputs for which the approximation is effective.
#' @param k A positive integer that sets the number of scaling points about which the quadrature approximation is performed. Increasing this value will expand the region of inputs over which the approximation is effective.
#' @param constr_id (Optional) A numeric value representing the constraint ID.
#' @rdname OpRelEntrConeQuad-class
## #' @export
OpRelEntrConeQuad <- function(X, Y, Z, m, k, constr_id = NA_integer_) { .OpRelEntrConeQuad(X = X, Y = Y, Z = Z, m = m, k = k, constr_id = constr_id) }

setMethod("initialize", "OpRelEntrConeQuad", function(.Object, ..., X, Y, Z, m, k) {
  .Object@X <- as.Constant(X)
  .Object@Y <- as.Constant(Y)
  .Object@Z <- as.Constant(Z)
  if(!is_hermitian(X) || !is_hermitian(Y) || !is_hermitian(Z))
    stop("One of the input matrices has not explicitly been declared as symmetric or"
         "Hermitian. If the inputs are Variable objects, try declaring them with the"
         "symmetric=True or Hermitian=True properties. If the inputs are general "
         "Expression objects that are known to be symmetric or Hermitian, then you"
         "can wrap them with the symmetric_wrap and hermitian_wrap atoms. Failure to"
         "do one of these things will cause this function to impose a symmetry or"
         "conjugate-symmetry constraint internally, in a way that is very"
         "inefficient.")
  .Object@m <- m
  .Object@k <- k
  Xs <- dim(.Object@X)
  Ys <- dim(.Object@Y)
  Zs <- dim(.Object@Z)
  if(!all(Xs == Ys) || !all(Xs == Zs))
    stop("All arguments must have the same shapes. Provided arguments have shapes ", Xs, ", ", Ys, ", and ", Zs)
  callNextMethod(.Object, ..., args = list(.Object@X, .Object@Y, .Object@Z))
})

#' @describeIn OpRelEntrConeQuad Information needed to reconstruct the object aside from the args.
setMethod("get_data", "OpRelEntrConeQuad", function(object) { list(object@m, object@k, id(object)) })

setMethod("show", "OpRelEntrConeQuad", function(object) {
  print(paste("OpRelEntrConeQuad(", as.character(object@X), ", ", as.character(object@Y), ", ", as.character(object@Z), ", ", object@m, ", ", object@k, ")", sep = ""))
})

#' @rdname OpRelEntrConeQuad-class
setMethod("as.character", "OpRelEntrConeQuad", function(x) {
  paste("OpRelEntrConeQuad(", as.character(x@X), ", ", as.character(x@Y), ", ", as.character(x@Z), ", ", x@m, ", ", x@k, ")", sep = "")
})

setMethod("residual", "OpRelEntrConeQuad", function(object) {
  stop("Unimplemented")
})

#' @describeIn OpRelEntrConeQuad The number of entries in the combined cones.
setMethod("size", "OpRelEntrConeQuad", function(object) { 3*num_cones(object) })

#' @describeIn OpRelEntrConeQuad The number of elementwise cones.
setMethod("num_cones", "OpRelEntrConeQuad", function(object) { size(object@X) })

#' @describeIn OpRelEntrConeQuad The dimensions of the exponential cones.
setMethod("cone_sizes", "OpRelEntrConeQuad", function(object) { rep(3, num_cones(object)) })

#' @describeIn OpRelEntrConeQuad An exponential constraint is DCP if each argument is affine.
setMethod("is_dcp", "OpRelEntrConeQuad", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping
    return(all(sapply(object@args, is_affine)))
  }
  return(all(sapply(object@args, is_affine)))
})

#' @describeIn OpRelEntrConeQuad Is the constraint DGP?
setMethod("is_dgp", "OpRelEntrConeQuad", function(object, dpp = FALSE) { FALSE })

#' @describeIn OpRelEntrConeQuad Is the constraint DQCP?
setMethod("is_dqcp", "OpRelEntrConeQuad", function(object) { is_dcp(object) })

setMethod("dim", "OpRelEntrConeQuad", function(x) { c(3, dim(x@X)) })

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn OpRelEntrConeQuad Replaces the dual values of an exponential cone constraint.
setReplaceMethod("dual_value", "OpRelEntrConeQuad", function(object, value) {
  stop("Unimplemented")
})

#'
#' The PowCone3D class.
#'
#' This class represents a collection of 3D power cone constraints
#'
#' \deqn{x[i]^alpha[i] * y[i]^(1-alpha[i]) \geq |z[i]| \text{for all} i, x \geq 0, y \geq 0}
#'
#' If the parameter alpha is a scalar, it will be promoted to a vector matching
#' the (common) sizes of x, y, z. The numeric value of alpha (or its components,
#' in the vector case) must be a number in the open interval (0, 1).
#'
#' We store flattened representations of the arguments (x, y, z, and alpha) as
#' Expression objects. We construct dual variables with respect to these
#' flattened representations.
#'
#' @slot x An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{x}.
#' @slot y An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{y}.
#' @slot z An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{z}.
#' @slot alpha An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{\alpha}. Must be in the open interval (0, 1).
#' @name PowCone3D-class
#' @aliases PowCone3D
#' @rdname PowCone3D-class
.PowCone3D <- setClass("PowCone3D", representation(x = "ConstValORExpr", y = "ConstValORExpr", z = "ConstValORExpr", alpha = "ConstValORExpr"), contains = "Constraint")

#' @param x An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{x}.
#' @param y An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{y}.
#' @param z An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{z}.
#' @param alpha An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{\alpha}. Must be in the open interval (0, 1).
#' @rdname PowCone3D-class
PowCone3D <- function(x, y, z, alpha, constr_id = NA_integer_) { .PowCone3D(x = x, y = y, z = z, alpha = alpha, constr_id = constr_id) }

setMethod("initialize", "PowCone3D", function(.Object, ..., x, y, z, alpha) {
  .Object@x <- as.Constant(x)
  .Object@y <- as.Constant(y)
  .Object@z <- as.Constant(z)
  for(val in list(.Object@x, .Object@y, .Object@z)) {
    if(!(is_affine(val) && is_real(val)))
      stop("All arguments must be affine and real")
  }

  alpha <- as.Constant(alpha)
  if(is_scalar(alpha))
    alpha <- Promote(alpha, dim(.Object@x))
  .Object@alpha <- alpha

  alpha_val <- value(.Object@alpha)
  if(any(alpha_val <= 0) || any(alpha_val >= 1))
    stop("alpha must have entries in the open interval (0, 1)")
  arg_dims <- list(dim(.Object@x), dim(.Object@y), dim(.Object@z), dim(.Object@alpha))
  for(i in 2:length(arg_dims)) {
    s <- arg_dims[[i]]
    if(any(arg_dims[[1]] != s))
      stop("All arguments must have the same dimensions")
  }
  callNextMethod(.Object, ..., args = list(.Object@x, .Object@y, .Object@z))
})

#' @rdname PowCone3D-class
setMethod("as.character", "PowCone3D", function(x) {
  paste("PowCone3D(", as.character(x@x), ", ", as.character(x@y), ", ", as.character(x@z), ", ", as.character(x@alpha), ")", sep = "")
})

#' @describeIn PowCone3D A \linkS4class{Expression} representing the residual of the constraint.
setMethod("residual", "PowCone3D", function(object) {
  # TODO: The projection should be implemented directly.
  if(is.na(value(object@x)) || is.na(value(object@y)) || is.na(value(object@z)))
    return(NA_real_)

  x <- new("Variable", dim = dim(object@x))
  y <- new("Variable", dim = dim(object@y))
  z <- new("Variable", dim = dim(object@z))
  constr <- list(PowCone3D(x, y, z, object@alpha))
  obj <- Minimize(Norm2(HStack(x, y, z) - HStack(value(object@x), value(object@y), value(object@z))))
  problem <- Problem(obj, constr)
  result <- solve(problem, solver = "SCS", eps = 1e-8)
  return(result$value)
})

#' @describeIn PowCone3D Information needed to reconstruct the object aside from the args.
setMethod("get_data", "PowCone3D", function(object) { list(object@alpha, id(object)) })

#' @describeIn PowCone3D A logical value indicating whether the constraint is imaginary.
setMethod("is_imag", "PowCone3D", function(object) { FALSE })

#' @describeIn PowCone3D A logical value indicating whether the constraint is complex.
setMethod("is_complex", "PowCone3D", function(object) { FALSE })

#' @describeIn PowCone3D The number of entries in the combined cones.
setMethod("size", "SOC", function(object) { 3*num_cones(object) })

#' @describeIn PowCone3D The number of elementwise cones.
setMethod("num_cones", "PowCone3D", function(object) { size(object@x) })

#' @describeIn PowCone3D The dimensions of the second-order cones.
setMethod("cone_sizes", "PowCone3D", function(object) { rep(3, num_cones(object)) })

#' @describeIn PowCone3D The constraint is DCP if the constrained expression is affine.
setMethod("is_dcp", "PowCone3D", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping
    args_ok <- all(sapply(object@args, is_affine))
    exps_ok <- !is(object@alpha, "Parameter")
    return(args_ok && exps_ok)
  }
  return(all(sapply(object@args, is_affine)))
})

#' @describeIn PowCone3D Is the constraint DGP?
setMethod("is_dgp", "PowCone3D", function(object, dpp = FALSE) { FALSE })

#' @describeIn PowCone3D Is the constraint DQCP?
setMethod("is_dqcp", "PowCone3D", function(object) { is_dcp(object) })

#' @describeIn PowCone3D The dimensions of the constrained expression.
setMethod("dim", "PowCone3D", function(x) {
  s <- c(3, dim(x@x))
  # Note: this can be a 3-tuple of ndim(x) == 2.
  return(s)
})

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn PowCone3D Replaces the dual values of a second-order cone constraint.
setReplaceMethod("dual_value", "PowCone3D", function(object, value) {
  # TODO: figure out why the reshaping has to be done differently,
  # compared to ExpCone constraints.
  value <- t(matrix(t(value), ncol = 3, byrow = FALSE))
  dv0 <- matrix(value[1,], nrow = nrow(object@x), ncol = ncol(object@x))
  dv1 <- matrix(value[2,], nrow = nrow(object@y), ncol = ncol(object@y))
  dv2 <- matrix(value[3,], nrow = nrow(object@z), ncol = ncol(object@z))

  value(object@dual_variables[[1]]) <- dv0
  value(object@dual_variables[[2]]) <- dv1
  value(object@dual_variables[[3]]) <- dv2
  return(object)
})

#'
#' The PowConeND class.
#'
#' This class represents a collection of N-dimensional power cone constraints
#' that is mathematically equivalent to the following code snippet:
#'
#' \code{apply(W^alpha, axis, prod) >= abs(z)},
#' W >= 0
#'
#' All arguments must be Expression-like, and z must satisfy ndim(z) <= 1. The
#' rows (resp. columns) of alpha must sum to 1 when axis = 1 (resp. axis = 2).
#'
#' Note: unlike PowCone3D, we make no attempt to promote alpha to the
#' appropriate shape. The dimensions of W and alpha must match exactly.
#'
#' Note: Dual variables are not currently implemented for this type of constraint.
#'
#' @slot W An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{W}.
#' @slot z An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{z}.
#' @slot alpha An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{\alpha}. Must be in the open interval (0, 1).
#' @slot axis The dimension along which to constrain: \code{1} indicates rows, and \code{2} indicates columns. The default is \code{2}.
#' @name PowConeND-class
#' @aliases PowConeND
#' @rdname PowConeND-class
.PowConeND <- setClass("PowConeND", representation(W = "ConstValORExpr", z = "ConstValORExpr", alpha = "ConstValORExpr", axis = "numeric"), prototype(axis = 2),
                       validity = function(object) {
                                      if(length(axis) > 1 || (axis != 1 && axis != 2))
                                        stop("[PowConeND: axis] axis must be either 1 (rows) or 2 (columns)")
                                      return(TRUE)
                                    }, contains = "Constraint")

#' @param W An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{W}.
#' @param z An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{z}.
#' @param alpha An \linkS4class{Expression}, numeric element, vector, or matrix representing \eqn{\alpha}. Must be in the open interval (0, 1).
#' @param axis The dimension along which to constrain: \code{1} indicates rows, and \code{2} indicates columns. The default is \code{2}.
#' @rdname PowConeND-class
PowConeND <- function(W, z, alpha, axis = 2, constr_id = NA_integer_) { .PowConeND(W = W, z = z, alpha = alpha, axis = axis, constr_id = constr_id) }

setMethod("initialize", "PowConeND", function(.Object, ..., W, z, alpha, axis = 2) {
  W <- as.Constant(W)
  if(!(is_real(W) && is_affine(W)))
    stop("Invalid first argument; W must be affine and real.")

  z <- as.Constant(z)
  # if(!(ndim(z) <= 1 || (ndim(z) == 2 && ncol(z) == 1)) || !(is_real(z) && is_affine(z)))
  if(ndim(z) > 1 || !(is_real(z) && is_affine(z)))
    stop("Invalid second argument. z must be affine, real, and have at most one ndim(z) <= 1.")

  # Check z has one entry per cone.
  if((ndim(W) <= 1 && size(z) > 1) ||
     (ndim(W) == 2 && size(z) != dim(W)[axis]) ||
     (ndim(W) == 1 && axis == 1))
    stop("Argument dimensions and axis are incompatible")

  axis_opp <- ifelse(axis == 1, 2, 1)
  if(ndim(W) == 2 && dim(W)[axis_opp] <= 1)
    stop("PowConeND requires left-hand-side to have at least two terms.")

  alpha <- as.Constant(alpha)
  if(any(dim(alpha) != dim(W)))
    stop("W and alpha dimensions must be equal")
  if(any(value(alpha) <= 0))
    stop("Argument alpha must be entry-wise positive.")
  if(any(abs(1 - apply(value(alpha), axis_opp, sum)) > 1e-6))
    stop("Argument alpha must sum to 1 along specified axis.")

  .Object@W <- W
  .Object@z <- z
  .Object@alpha <- alpha
  .Object@axis <- axis
  if(ndim(z) == 0)
    z <- flatten(z)

  callNextMethod(.Object, ..., args = list(W, z))
})

#' @rdname PowConeND-class
setMethod("as.character", "PowConeND", function(x) {
  paste("PowConeND(", as.character(x@x), ", ", as.character(x@W), ", ", as.character(x@z), ", ", as.character(x@alpha), ")", sep = "")
})

#' @describeIn PowConeND A logical value indicating whether the constraint is imaginary.
setMethod("is_imag", "PowConeND", function(object) { FALSE })

#' @describeIn PowConeND A logical value indicating whether the constraint is complex.
setMethod("is_complex", "PowConeND", function(object) { FALSE })

#' @describeIn PowConeND Information needed to reconstruct the object aside from the args.
setMethod("get_data", "PowConeND", function(object) { list(object@alpha, object@axis, id(object)) })

#' @describeIn PowConeND A \linkS4class{Expression} representing the residual of the constraint.
setMethod("residual", "PowConeND", function(object) {
  # TODO: The projection should be implemented directly.
  if(is.na(value(object@W)) || is.na(value(object@z)))
    return(NA_real_)

  W <- new("Variable", dim = dim(object@W))
  z <- new("Variable", dim = dim(object@z))
  constr <- list(PowConeND(W, z, object@alpha, axis = object@axis))
  obj <- Minimize(Norm2(HStack(flatten(W), flatten(z)) -
                        HStack(value(flatten(object@W)), value(flatten(object@z)))))
  problem <- Problem(obj, constr)
  result <- solve(problem, solver = "SCS", eps = 1e-8)
  return(result$value)
})

#' @describeIn PowConeND The number of elementwise cones.
setMethod("num_cones", "PowConeND", function(object) { size(object@z) })

#' @describeIn PowConeND The number of entries in the combined cones.
setMethod("size", "PowConeND", function(object) {
  axis_opp <- ifelse(object@axis == 1, 2, 1)
  cone_size <- 1 + dim(object@args[[1]])[axis_opp]
  return(cone_size * num_cones(object))
})

#' @describeIn PowConeND The dimensions of the second-order cones.
setMethod("cone_sizes", "PowConeND", function(object) {
  axis_opp <- ifelse(object@axis == 1, 2, 1)
  cone_size <- 1 + dim(object@args[[1]])[axis_opp]
  rep(cone_size, num_cones(object))
})

#' @describeIn PowConeND The constraint is DCP if the constrained expression is affine.
setMethod("is_dcp", "PowConeND", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping
    args_ok <- (is_affine(object@args[[1]]) && is_affine(object@args[[2]]))
    exps_ok <- !is(object@alpha, "Parameter")
    return(args_ok && exps_ok)
  }
  return(TRUE)
})

#' @describeIn PowConeND Is the constraint DGP?
setMethod("is_dgp", "PowConeND", function(object, dpp = FALSE) { FALSE })

#' @describeIn PowConeND Is the constraint DQCP?
setMethod("is_dqcp", "PowConeND", function(object) { is_dcp(object) })

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn PowConeND Replaces the dual values of a second-order cone constraint.
setReplaceMethod("dual_value", "PowConeND", function(object, value) {
  stop("Unimplemented")
})

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

#'
#' The SOC class.
#'
#' This class represents a second-order cone constraint for each row/column, i.e. \eqn{\|x\|_2 \leq t}.
#' It assumes t is a vector the same length as X's rows (resp. columns) for axis == 1 (resp. 2).
#'
#' @slot t The scalar part of the second-order constraint.
#' @slot X A matrix whose rows/columns are each a cone.
#' @slot axis The dimension along which to constrain: \code{1} indicates rows, and \code{2} indicates columns. The default is \code{2}.
#' @name SOC-class
#' @aliases SOC
#' @rdname SOC-class
.SOC <- setClass("SOC", representation(t = "ConstValORExpr", X = "ConstValORExpr", axis = "numeric"),
                        prototype(t = NA_real_, X = NA_real_, axis = 2),
                        validity = function(object) {
                          if(length(object@axis) > 1 || (object@axis != 1 && object@axis != 2))
                            stop("[SOC: axis] axis must be either 1 (rows) or 2 (columns)")
                          return(TRUE)
                        }, contains = "Constraint")

#' @param t The scalar part of the second-order constraint.
#' @param X A matrix whose rows/columns are each a cone.
#' @param axis The dimension along which to slice: \code{1} indicates rows, and \code{2} indicates columns. The default is \code{2}.
#' @param constr_id (Optional) A numeric value representing the constraint ID.
#' @rdname SOC-class
## #' @export
SOC <- function(t, X, axis = 2, constr_id = NA_integer_) { .SOC(t = t, X = X, axis = axis, constr_id = constr_id) }

setMethod("initialize", "SOC", function(.Object, ..., t, X, axis = 2) {
  t_dim <- dim(t)
  t_size <- size(t)
  X_dim <- dim(X)
  if(length(t_dim) >= 2 || !is_real(t))
    stop("Invalid first argument")

  # Check t has one entry per cone.
  if((length(X_dim) <= 1 and t_size > 1) ||
     (length(X_dim) == 2 && t_size != X_dim[axis]) ||
     (length(X_dim) == 1 && axis == 1))
    stop("Argument dimensions and axis are incompatible")

  if(is.null(t_dim) || length(t_dim) == 0)
    t <- flatten(t)

  .Object@t <- t
  .Object@X <- X
  .Object@axis <- axis
  callNextMethod(.Object, ..., args = list(t, X))
})

#' @param x,object A \linkS4class{SOC} object.
#' @rdname SOC-class
setMethod("as.character", "SOC", function(x) {
  paste("SOC(", as.character(x@t), ", ", as.character(x@X), ")", sep = "")
})

#' @describeIn SOC The residual of the second-order constraint.
setMethod("residual", "SOC", function(object) {
  # For each cone, returns:
  #
  #      ||(t,X) - proj(t,X)||
  #      with
  #      proj(t,X) = (t,X)                       if t >= ||x||
  #                  0.5*(t/||x|| + 1)(||x||,x)  if -||x|| < t < ||x||
  #                  0                           if t <= -||x||
  #      References:
  #           https://docs.mosek.com/modeling-cookbook/practical.html#distance-to-a-cone
  #           https://math.stackexchange.com/questions/2509986/projection-onto-the-second-order-cone

  t <- value(object@args[[1]])
  X <- value(object@args[[2]])
  if(is.na(t) || is.na(X))
    return(NA)

  # Reduce axis = 2 to axis = 1.
  if(object@axis == 2)
    X <- t(X)

  promoted <- (ndim(X) == 1)
  X <- matrix(X)

  # Initializing with zeros makes "0 if t <= -||x||" the default case for the projection
  if(is.null(dim(t)))
    t_proj <- matrix(0)
  else
    t_proj <- matrix(0, nrow = nrow(t), ncol = ncol(t))
  if(is.null(dim(X)))
    X_proj <- matrix(0)
  else
    X_proj <- matrix(0, nrow = nrow(X), ncol = ncol(X))

  norms <- apply(X, 1, function(row) { norm(row, "2") })

  # 1. proj(t,X) = (t,X) if t >= ||x||
  t_geq_x_norm <- (t >= norms)
  t_proj[t_geq_x_norm] <- t[t_geq_x_norm]
  X_proj[t_geq_x_norm] <- X[t_geq_x_norm]

  # 2. proj(t,X) = 0.5*(t/||x|| + 1)(||x||,x)  if -||x|| < t < ||x||
  abs_t_less_x_norm <- (abs(t) < norms)
  avg_coeff <- 0.5 * (1 + t/norms)
  X_proj[abs_t_less_x_norm] <- diag(avg_coeff[abs_t_less_x_norm]) %*% X[abs_t_less_x_norm]
  t_proj[abs_t_less_x_norm] <- avg_coeff[abs_t_less_x_norm] * norms[abs_t_less_x_norm]

  Xt <- cbind(X, t)
  Xt_proj <- cbind(X_proj, t_proj)
  resid <- apply(Xt - Xt_proj, 1, function(row) { norm(row, "2") })

  # Demote back to 1D.
  if(promoted)
    return(resid[1])
  else
    return(resid)
})

#' @describeIn SOC Information needed to reconstruct the object aside from the args.
setMethod("get_data", "SOC", function(object) { list(object@axis, id(object)) })

#' @describeIn SOC The number of elementwise cones.
setMethod("num_cones", "SOC", function(object) { size(object@args[[1]]) })

#' @describeIn SOC The number of entries in the combined cones.
setMethod("size", "SOC", function(object) {
  if(object@axis == 2)   # Collapse columns.
    idx <- 1
  else if(object@axis == 1)   # Collapse rows.
    idx <- 2
  else
    stop("Unimplemented")
  cone_size <- 1 + dim(object@args[[2]])[idx]
  return(cone_size * num_cones(object))
})

#' @describeIn SOC The dimensions of the second-order cones.
setMethod("cone_sizes", "SOC", function(object) {
  if(object@axis == 2)   # Collapse columns.
    idx <- 1
  else if(object@axis == 1)   # Collapse rows.
    idx <- 2
  else
    stop("Unimplemented")
  cone_size <- 1 + dim(object@args[[2]])[idx]
  return(rep(cone_size, num_cones(object)))
})

#' @describeIn SOC An SOC constraint is DCP if each of its arguments is affine.
setMethod("is_dcp", "SOC", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping.
    return(all(sapply(object@args, is_affine)))
  }
  return(all(sapply(object@args, is_affine)))
})

#' @describeIn SOC Is the constraint DGP?
setMethod("is_dgp", "SOC", function(object, dpp = FALSE) { FALSE })

#' @describeIn SOC Is the constraint DQCP?
setMethod("is_dqcp", "SOC", function(object) { is_dcp(object) })

#' @param value A numeric scalar, vector, or matrix.
#' @describeIn SOC Replaces the dual values of a second-order cone constraint.
setReplaceMethod("dual_value", "SOC", function(object, value) {
  if(object@axis == 2)   # Collapse columns.
    idx <- 1
  else if(object@axis == 1)   # Collapse rows.
    idx <- 2
  else
    stop("Unimplemented")

  cone_size <- 1 + dim(object@args[[2]])[idx]
  value <- t(matrix(t(value), nrow = cone_size, byrow = FALSE))

  t <- value[,1]
  X <- value[,2:ncol(value)]
  if(object@axis == 2)
    X <- t(X)

  value(object@dual_variables[[1]]) <- t
  value(object@dual_variables[[2]]) <- X
  return(object)
})

## Deleted duplicated and not-in-effect functions:
## format_axis, format_elemwise, get_spacing_matrix
