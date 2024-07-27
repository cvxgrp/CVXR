## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/lambda_max_canon.py
#'
#' Dcp2Cone canonicalizer for the lambda maximization atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a lambda maximization atom
#' where t is the objective function and a PSD constraint and a
#' constraint requiring I*t to be symmetric.
Dcp2Cone.lambda_max_canon <- function(expr, args) {
  A <- args[[1]]
  n <- nrow(A)
  t <- Variable()
  prom_t <- promote(t, c(n,1))
  # Constrain I*t - A to be PSD; note this expression must be symmetric
  tmp_expr <- DiagVec(prom_t) - A
  constr <- list(PSDConstraint(tmp_expr))
  if(!is_symmetric(A)) {
    ut <- UpperTri(A)
    lt <- UpperTri(t(A))
    constr <- c(constr, ut == lt)
  }
  return(list(t, constr))
}

