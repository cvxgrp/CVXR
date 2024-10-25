## CVXPY SOURCE: cvxpy/reductions/dgp2dcp/atom_canonicalizers/mul_canon.py
#'
#' Dgp2Dcp canonicalizer for the multiplication atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the multiplication atom of a DGP expression,
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.mul_canon <- function(expr, args) {
  # expr <- NULL
  return(list(AddExpression(args), list()))
}

