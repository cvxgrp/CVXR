## CVXPY SOURCE: cvxpy/problems/objective.py

#'
#' The Objective class.
#'
#' This class represents an optimization objective.
#'
#' @slot expr A scalar \linkS4class{Expression} to optimize.
#' @slot NAME (Internal) Name used for debug messages.
#' @name Objective-class
#' @aliases Objective
#' @rdname Objective-class
.Objective <- setClass("Objective",
                       representation(expr = "ConstValORExpr", NAME = "character"),
                       contains = "Canonical")

#' @param expr A scalar \linkS4class{Expression} to optimize.
#' @rdname Objective-class
Objective <- function(expr) { .Objective(expr = expr) }

setMethod("initialize", "Objective", function(.Object, ..., expr) {
  .Object@expr <- expr
  .Object@NAME  <- "objective"
  .Object@args <- list(as.Constant(expr))
  # .Object <- callNextMethod(.Object, ..., args = list(as.Constant(expr)))

  # Validate that the objective resolves to a scalar.
  if(!is_scalar(.Object@args[[1]]))
    stop("The objective must resolve to a scalar")
  if(!is_real(.Object@args[[1]]))
    stop("The objective must be real valued")
  return(.Object)
})

setMethod("show", "Objective", function(object) {
  cat("Objective(", name(object@args[[1]]), ")", sep = "")
})

setMethod("as.character", "Objective", function(x) {
  paste("Objective", name(x@args[[1]]))
})

#' @param object An \linkS4class{Objective} object.
#' @describeIn Objective The value of the objective expression.
setMethod("value", "Objective", function(object) {
  v <- value(object@args[[1]])
  if(is.na(v))
    return(NA_real_)
  else
    return(intf_scalar_value(v))
})

#' @describeIn Objective Is the objective a quadratic function?
setMethod("is_quadratic", "Objective", function(object) {
  is_quadratic(object@args[[1]])
})

#' @describeIn Objective Is the objective a quadratic of piecewise affine function?
setMethod("is_qpwa", "Objective", function(object) {
  is_qpwa(object@args[[1]])
})

setClassUnion("ListORNULL", c("list", "NULL"))

#'
#' Arithmetic Operations on Objectives
#'
#' Add, subtract, multiply, or divide optimization objectives.
#'
#' @param e1 The left-hand \linkS4class{Minimize}, \linkS4class{Maximize}, or numeric value.
#' @param e2 The right-hand \linkS4class{Minimize}, \linkS4class{Maximize}, or numeric value.
#' @return A \linkS4class{Minimize} or \linkS4class{Maximize} object.
#' @name Objective-arith
NULL

#============================#
# Objective Class Arithmetic
#============================#
#' @rdname Objective-arith
setMethod("+", signature(e1 = "Objective", e2 = "numeric"), function(e1, e2) { if(length(e2) == 1 && e2 == 0) e1 else stop("Unimplemented") })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "numeric", e2 = "Objective"), function(e1, e2) { e2 + e1 })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Objective", e2 = "numeric"), function(e1, e2) { if(length(e2) == 1 && e2 == 0) e1 else stop("Unimplemented") })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "numeric", e2 = "Objective"), function(e1, e2) { if(length(e1) == 0 && e1 == 0) -e2 else stop("Unimplemented") })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Objective", e2 = "Minimize"), function(e1, e2) { e1 + (-e2) })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Objective", e2 = "Maximize"), function(e1, e2) { e1 + (-e2) })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Minimize", e2 = "Objective"), function(e1, e2) { (-e2) + e1 })

#' @rdname Objective-arith
setMethod("-", signature(e1 = "Maximize", e2 = "Objective"), function(e1, e2) { (-e2) + e1 })

#' @rdname Objective-arith
setMethod("*", signature(e1 = "Minimize", e2 = "numeric"), function(e1, e2) {
  if(e2 >= 0) Minimize(expr = e1@args[[1]] * e2) else Maximize(expr = e1@args[[1]] * e2)
})

#' @rdname Objective-arith
setMethod("*", signature(e1 = "Maximize", e2 = "numeric"), function(e1, e2) {
  if(e2 < 0) Minimize(expr = e1@args[[1]] * e2) else Maximize(expr = e1@args[[1]] * e2)
})

#' @rdname Objective-arith
setMethod("*", signature(e1 = "numeric", e2 = "Minimize"), function(e1, e2) { e2 * e1 })

#' @rdname Objective-arith
setMethod("*", signature(e1 = "numeric", e2 = "Maximize"), function(e1, e2) { e2 * e1 })

#' @rdname Objective-arith
setMethod("/", signature(e1 = "Objective", e2 = "numeric"), function(e1, e2) { e1 * (1.0/e2) })



#'
#' The Minimize class.
#'
#' This class represents an optimization objective for minimization.
#'
#' @slot expr A scalar \linkS4class{Expression} to minimize.
#' @name Minimize-class
#' @aliases Minimize
#' @rdname Minimize-class
.Minimize <- setClass("Minimize", representation(expr = "ConstValORExpr"), contains = "Objective")

#' @param expr A scalar \linkS4class{Expression} to minimize.
#' @rdname Minimize-class
#' @export
Minimize <- function(expr) { .Minimize(expr = expr) }

#' @param object A \linkS4class{Minimize} object.
#' @describeIn Minimize Pass on the target expression's objective and constraints.
setMethod("canonicalize", "Minimize", function(object) { canonical_form(object@args[[1]]) })

#' @describeIn Minimize A logical value indicating whether the objective is convex.
setMethod("is_dcp", "Minimize", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping.
    return(is_convex(object@args[[1]]))
  }
  return(is_convex(object@args[[1]]))
})

#' @describeIn Minimize A logical value indicating whether the objective is log-log convex.
setMethod("is_dgp", "Minimize", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping.
    return(is_log_log_convex(object@args[[1]]))
  }
  return(is_log_log_convex(object@args[[1]]))
})

#' @describeIn Minimize Is the objective DPP?
setMethod("is_dpp", "Minimize", function(object, context = "dcp") {
  dpp_scope()   # TODO: Implement DPP scoping.
  if(tolower(context) == "dcp")
    return(is_dcp(object, dpp = TRUE))
  else if(tolower(context) == "dgp")
    return(is_dgp(object, dpp = TRUE))
  else
    stop("Unsupported context ", context)
})

#' @describeIn Minimize A logical value indicating whether the objective is quasiconvex.
setMethod("is_dqcp", "Minimize", function(object) {
  is_quasiconvex(object@args[[1]])
})

# The value of the objective given the solver primal value.
Minimize.primal_to_result <- function(result) { result }

#===========================#
# Minimize Class Arithmetic
#===========================#
#' @rdname Objective-arith
setMethod("-", signature(e1 = "Minimize", e2 = "missing"), function(e1, e2) { Maximize(expr = -e1@args[[1]]) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Minimize", e2 = "Minimize"), function(e1, e2) { Minimize(e1@args[[1]] + e2@args[[1]]) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Minimize", e2 = "Maximize"), function(e1, e2) { stop("Problem does not follow DCP rules") })


#'
#' The Maximize class.
#'
#' This class represents an optimization objective for maximization.
#'
#' @slot expr A scalar \linkS4class{Expression} to maximize.
#' @examples
#' x <- Variable(3)
#' alpha <- c(0.8,1.0,1.2)
#' obj <- sum(log(alpha + x))
#' constr <- list(x >= 0, sum(x) == 1)
#' prob <- Problem(Maximize(obj), constr)
#' result <- solve(prob)
#' result$value
#' result$getValue(x)
#' @name Maximize-class
#' @aliases Maximize
#' @rdname Maximize-class
.Maximize <- setClass("Maximize", contains = "Objective")

#' @param expr A scalar \linkS4class{Expression} to maximize.
#' @rdname Maximize-class
#' @export
Maximize <- function(expr) { .Maximize(expr = expr) }

#' @param object A \linkS4class{Maximize} object.
#' @describeIn Maximize Negates the target expression's objective.
setMethod("canonicalize", "Maximize", function(object) {
  canon <- canonical_form(object@args[[1]])
  list(lu.neg_expr(canon[[1]]), canon[[2]])
})

#' @describeIn Maximize A logical value indicating whether the objective is concave.
setMethod("is_dcp", "Maximize", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping.
    return(is_concave(object@args[[1]]))
  }
  return(is_concave(object@args[[1]]))
})

#' @describeIn Maximize A logical value indicating whether the objective is log-log concave.
setMethod("is_dgp", "Maximize", function(object, dpp = FALSE) {
  if(dpp) {
    dpp_scope()   # TODO: Implement DPP scoping.
    return(is_log_log_concave(object@args[[1]]))
  }
  return(is_log_log_concave(object@args[[1]]))
})

#' @describeIn Maximize Is the objective DPP?
setMethod("is_dpp", "Maximize", function(object, context = "dcp") {
  dpp_scope()   # TODO: Implement DPP scoping.
  if(tolower(context) == "dcp")
    return(is_dcp(object, dpp = TRUE))
  else if(tolower(context) == "dgp")
    return(is_dgp(object, dpp = TRUE))
  else
    stop("Unsupported context ", context)
})

#' @describeIn Maximize A logical value indicating whether the objective is quasiconcave.
setMethod("is_dqcp", "Maximize", function(object) { is_quasiconcave(object@args[[1]]) })

# The value of the objective given the solver primal value.
Maximize.primal_to_result <- function(result) { -result }

#===========================#
# Maximize Class Arithmetic
#===========================#
#' @rdname Objective-arith
setMethod("-", signature(e1 = "Maximize", e2 = "missing"), function(e1, e2) { Minimize(expr = -e1@args[[1]]) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Maximize", e2 = "Maximize"), function(e1, e2) { Maximize(expr = e1@args[[1]] + e2@args[[1]]) })

#' @rdname Objective-arith
setMethod("+", signature(e1 = "Maximize", e2 = "Minimize"), function(e1, e2) { stop("Problem does not follow DCP rules") })

