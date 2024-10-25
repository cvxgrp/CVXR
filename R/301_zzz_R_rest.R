
# QPSolver requires objectives to be stuffed in the following way.
#'
#' Is the QP objective stuffed?
#'
#' @param objective A \linkS4class{Minimize} or \linkS4class{Maximize} object representing the optimization objective.
#' @return Is the objective a stuffed QP?
is_stuffed_qp_objective <- function(objective) {
  expr <- expr(objective)
  return(inherits(expr, "AddExpression") && length(expr@args) == 2 && inherits(expr@args[[1]], "QuadForm") && inherits(expr@args[[2]], "MulExpression") && is_affine(expr@args[[2]]))
}
