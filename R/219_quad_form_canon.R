## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/quad_form_canon.py

#'
#' Dcp2Cone canonicalizer for the quadratic form atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a quadratic form atom,
#' where the objective function consists of the scaled objective function
#' from the quadratic over linear canonicalization and same with the
#' constraints.
Dcp2Cone.quad_form_canon <- function(expr, args) {
  # TODO: This doesn't work with parameters!
  decomp <- .decomp_quad(value(args[[2]]))
  scale <- decomp[[1]]
  M1 <- decomp[[2]]
  M2 <- decomp[[3]]
  M1_dim <- dim(M1)
  M2_dim <- dim(M2)

  # Special case where P == 0.
  if((is.null(M1_dim) || size(M1) == 0) && (is.null(M2_dim) && size(M2) == 0))
    return(list(Constant(0), list()))

  if(!is.null(M1_dim) && prod(M1_dim) > 0)
    expr <- sum_squares(Constant(t(M1)) %*% args[[1]])
  else if(!is.null(M2_dim) && prod(M2_dim) > 0) {
    scale <- -scale
    expr <- sum_squares(Constant(t(M2)) %*% args[[1]])
  }
  canon <- Dcp2Cone.quad_over_lin_canon(expr, expr@args)
  obj <- canon[[1]]
  constr <- canon[[2]]
  return(list(scale * obj, constr))
}
