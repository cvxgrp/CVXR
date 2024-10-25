## CVXPY SOURCE: cvxpy/expression/constraints/finite_set.py

#'
#' The FiniteSet class.
#'
#' This class represents a constraint that each entry of an Expression to take a value in a given set of real numbers.
#'
#' @slot expre The given expression, which must be affine, to be constrained.
#' If expre has multiple elements, then the constraint is applied separately to
#' each element, i.e., after solving a problem with this constraint, we should have:
#' `for(e in flatten(expre)) print(value(e) %in% vec) # => TRUE `
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
    saved_scope <- .CVXR_options$dpp_scope_active
    .CVXR_options$dpp_scope_active <- TRUE
    on.exit({
      .CVXR_options$dpp_scope_active <- saved_scope
    })
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
