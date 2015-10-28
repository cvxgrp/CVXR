#'
#' The Minimize class.
#'
#' This class represents a minimization problem.
#'
#' @slot expr The expression to minimize.
#' @aliases Minimize
#' @export
Minimize <- setClass("Minimize", representation(expr = "ConstValORExpr"))

setMethod("initialize", "Minimize", function(.Object, expr) {
    .Object@expr <- as.Constant(expr)
    if(!all(size(.Object@expr) == c(1,1)))
      stop("The objective must resolve to a scalar")
    return(.Object)
})

setMethod("canonicalize", "Minimize", function(object) { canonical_form(object@expr) })
setMethod("variables", "Minimize", function(object) { variables(object@expr) })
setMethod("parameters", "Minimize", function(object) { parameters(object@expr) })
setMethod("is_dcp", "Minimize", function(object) { is_convex(object@expr) })

#'
#' The Maximize class.
#'
#' This class represents a maximization problem.
#'
#' @slot expr The expression to maximize.
#' @aliases Maximize
#' @export
Maximize <- setClass("Maximize", contains = "Minimize")

setMethod("-", signature(e1 = "Minimize", e2 = "missing"), function(e1, e2) { Maximize(expr = - e1@expr) })
setMethod("+", signature(e1 = "Minimize", e2 = "Minimize"), function(e1, e2) { Minimize(e1@expr + e2@expr) })
setMethod("+", signature(e1 = "Minimize", e2 = "Maximize"), function(e1, e2) { stop("Problem does not follow DCP rules") })
setMethod("-", signature(e1 = "Minimize", e2 = "Minimize"), function(e1, e2) { e1 + -e2 })
setMethod("-", signature(e1 = "Minimize", e2 = "Maximize"), function(e1, e2) { e1 + -e2 })
setMethod("*", signature(e1 = "Minimize", e2 = "numeric"), function(e1, e2) {
  ifelse(e2 >= 0, Minimize(expr = e1@expr * e2), Maximize(expr = e1@expr * e2))
})
setMethod("*", signature(e1 = "Maximize", e2 = "numeric"), function(e1, e2) {
  ifelse(e2 < 0, Minimize(expr = e1@expr * e2), Maximize(expr = e1@expr * e2))
})
setMethod("/", signature(e1 = "Minimize", e2 = "numeric"), function(e1, e2) { e1 * (1.0/e2) })

setMethod("-", signature(e1 = "Maximize", e2 = "missing"), function(e1, e2) { Minimize(-e1@expr) })
setMethod("+", signature(e1 = "Maximize", e2 = "Maximize"), function(e1, e2) { Maximize(expr = e1@expr + e2@expr) })
setMethod("+", signature(e1 = "Maximize", e2 = "Minimize"), function(e1, e2) { stop("Problem does not follow DCP rules") })

setMethod("canonicalize", "Maximize", function(object) {
  canon <- callNextMethod(object)
  obj <- canon[[1]]
  constraints <- canon[[2]]
  list(neg_expr(obj), constraints)
})

setMethod("is_dcp", "Maximize", function(object) { is_concave(object@expr) })

#'
#' The Problem class.
#'
#' This class represents the convex optimization problem.
#'
#' @slot objective The expression to minimize or maximize.
#' @slot constraints (Optional) A list of constraints on the problem variables.
#' @aliases Problem
#' @export
.Problem <- setClass("Problem", representation(objective = "Minimize", constraints = "list"),
                    prototype(constraints = list()),
                    validity = function(object) {
                      if(!(class(object@objective) %in% c("Minimize", "Maximize")))
                        stop("Problem objective must be Minimize or Maximize")
                      return(TRUE)
                    })

Problem <- function(objective, constraints = list(), ...) {
  .Problem(objective = objective, constraints = constraints, ...)
}

setMethod("is_dcp", "Problem", function(object) {
  is_dcp_list <- lapply(c(object@objective, object@constraints), is_dcp)
  all(unlist(is_dcp_list))
})

setMethod("+", signature(e1 = "Problem", e2 = "missing"), function(e1, e2) { Problem(objective = e1@objective, constraints = e1@constraints) })
setMethod("-", signature(e1 = "Problem", e2 = "missing"), function(e1, e2) { Problem(objective = -e1@objective, constraints = e1@constraints) })
setMethod("+", signature(e1 = "Problem", e2 = "Problem"), function(e1, e2) {
  Problem(objective = e1@objective + e2@objective, constraints = unique(c(e1@constraints, e2@constraints)))
})
setMethod("-", signature(e1 = "Problem", e2 = "Problem"), function(e1, e2) {
  Problem(objective = e1@objective - e2@objective, constraints = unique(c(e1@constraints, e2@constraints)))
})
setMethod("*", signature(e1 = "Problem", e2 = "numeric"), function(e1, e2) {
  Problem(objective = e1@objective * e2, constraints = e1@constraints)
})
setMethod("*", signature(e1 = "numeric", e2 = "Problem"), function(e1, e2) { 
  Problem(objective = e2 * e1@objective, constraints = e1@constraints)
})
setMethod("/", signature(e1 = "Problem", e2 = "numeric"), function(e1, e2) {
  Problem(objective = e1@objective * (1.0/e2), constraints = e1@constraints)
})

setMethod("canonicalize", "Problem", function(object) {
  obj_canon <- canonical_form(object@objective)
  canon_constr <- lapply(object@constraints, function(x) { canonical_form(x)[2] })
  list(obj_canon[1], c(obj_canon[2], canon_constr))
})

setMethod("variables", "Problem", function(object) {
  vars_ <- variables(object@objective)
  constrs_ <- lapply(object@constraints, function(constr) { variables(constr) })
  unique(flatten_list(c(vars_, constrs_)))
})

setMethod("parameters", "Problem", function(object) {
  params <- parameters(object@objective)
  constrs_ <- lapply(object@constraints, function(constr) { parameters(constr) })
  unique(flatten_list(c(params, constrs_)))
})
