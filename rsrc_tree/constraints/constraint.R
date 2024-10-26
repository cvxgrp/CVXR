## CVXPY SOURCE: cvxpy/expression/constraints/.py
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
setMethod("is_dpp", "Constraint", function(object, context = c("dcp", "dgp")) {
  context <- match.arg(context)
  if(tolower(context) == "dcp")
    return(is_dcp(object, dpp = TRUE))
  else
    return(is_dgp(object, dpp = TRUE))

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


## Begin R-specific Code Section
## Update ListORConstr to include the Constraint class just defined above
setIs("Constraint", "ListORConstr")
## End R-specific Code Section

# Helper function since syntax is different for LinOp (list) vs. Constraint object
# @param object A list or \linkS4class{Constraint} object.
# Returns the ID associated with the list or constraint.
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
