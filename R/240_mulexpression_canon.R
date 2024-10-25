## CVXPY SOURCE: cvxpy/reductions/dgp2dcp/atom_canonicalizers/mulexpression_canon.py
#'
#' Dgp2Dcp canonicalizer for the multiplication expression atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of values for the expr variable
#' @return A canonicalization of the multiplication expression atom
#' of a DGP expression, where the returned expression is the transformed
#' DCP equivalent.
Dgp2Dcp.mulexpression_canon <- function(expr, args) {
  lhs <- args[[1]]
  rhs <- args[[2]]
  dims <- mul_dims_promote(dim(lhs), dim(rhs))
  lhs_dim <- dims[[1]]
  rhs_dim <- dims[[2]]
  lhs <- reshape_expr(lhs, lhs_dim)
  rhs <- reshape_expr(rhs, rhs_dim)
  rows <- list()

  # TODO: Parallelize this for large matrices.
  for(i in seq_len(nrow(lhs))) {
    row <- list()
    for(j in seq_len(ncol(rhs))) {
      hstack_args <- lapply(seq_len(ncol(lhs)), function(k) { lhs[i,k] + rhs[k,j] })
      row <- c(row, list(log_sum_exp(do.call("HStack", hstack_args))))
    }
    rows <- c(rows, list(row))
  }
  mat <- bmat(rows)
  if(!all(dim(mat) == dim(expr)))
    mat <- reshape_expr(mat, dim(expr))
  return(list(mat, list()))
}

