## CVXPY SOURCE: cvxpy/reductions/qp2quad_form/qp2symbolic_qp.py

## Uses definitions in dcpcanon.R which needs to precede this when sourcing!

#'
#' The Qp2SymbolicQp class.
#'
#' This class reduces a quadratic problem to a problem that consists of affine
#' expressions and symbolic quadratic forms.
#'
#' @rdname Qp2SymbolicQp-class
.Qp2SymbolicQp <- setClass("Qp2SymbolicQp", contains = "Canonicalization")
Qp2SymbolicQp <- function(problem = NULL) { .Qp2SymbolicQp(problem = problem) }

setMethod("initialize", "Qp2SymbolicQp", function(.Object, ...) {
  callNextMethod(.Object, ..., canon_methods = Qp2QuadForm.CANON_METHODS)
})

Qp2SymbolicQp.accepts <- function(problem) {
  is_qpwa(expr(problem@objective)) &&
  length(intersect(c("PSD", "NSD"), convex_attributes(variables(problem)))) == 0 &&
  all(sapply(problem@constraints, function(c) {
        (inherits(c, c("NonPosConstraint", "NonNegConstraint", "IneqConstraint")) && is_pwl(expr(c))) ||
        (inherits(c, c("ZeroConstraint", "EqConstraint")) && are_args_affine(list(c)))
  }))
}

# Problems with quadratic, piecewise affine objectives, piecewise-linear constraints, inequality constraints,
# and affine equality constraints are accepted.
setMethod("accepts", signature(object = "Qp2SymbolicQp", problem = "Problem"), function(object, problem) {
    Qp2SymbolicQp.accepts(problem)
})

# Converts a QP to an even more symbolic form.
setMethod("perform", signature(object = "Qp2SymbolicQp", problem = "Problem"), function(object, problem) {
  if(!accepts(object, problem))
    stop("Cannot reduce problem to symbolic QP")
  callNextMethod(object, problem)
})
