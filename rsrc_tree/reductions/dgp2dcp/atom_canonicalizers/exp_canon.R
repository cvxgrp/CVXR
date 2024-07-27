## CVXPY SOURCE: cvxpy/reductions/dgp2dcp/atom_canonicalizers/exp_canon.py

#'
#' Dgp2Dcp canonicalizer for the exp atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the exp atom of a DGP expression,
#' where the returned expression is the transformed DCP equivalent.
Dgp2Dcp.exp_canon <- function(expr, args) {
  # expr <- NULL
  return(list(Exp(args[[1]]), list()))
}

