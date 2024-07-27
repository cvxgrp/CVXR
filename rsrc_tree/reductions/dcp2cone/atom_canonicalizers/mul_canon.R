## CVXPY SOURCE: cvxpy/reductions/dcp2cone/atom_canonicalizers/mul_canon.py

#'
#' Dcp2Cone canonicalizer for the multiplication of expressions atom
#'
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from the multiplication of expressions atom.
Dcp2Cone.mul_canon <- function(expr, args) {
  # TODO(akshayka): expose as a reduction for user's convenience

  # Only allow param * var (not var * param). Associate right to left.
  # TODO: Only descend if both sides have parameters
  lhs <- args[[1]]
  rhs <- args[[2]]
  lhs_parms <- parameters(lhs)
  rhs_parms <- parameters(rhs)
  if((is.null(lhs_parms) || length(lhs_parms) == 0) && (is.null(rhs_parms) && length(rhs_parms) == 0))
    return(list(copy(expr, args), list()))

  op_type <- class(expr)
  if(length(variables(lhs)) > 0) {
    if(dpp_scope()) {   # TODO: Implement DPP scoping with global variables.
      if(!is_affine(rhs))
        stop("rhs must be affine if DPP")
    }
    t <- new("Variable", dim = dim(lhs))
    return(list(do.call(op_type, list(t, rhs)), list(t == lhs)))
  } else if(length(variables(rhs))) {
    if(dpp_scope()) {
      if(!is_affine(lhs))
        stop("lhs must be affine if DPP")
    }
    t <- new("Variable", dim = dim(rhs))
    return(list(do.call(op_type, list(lhs, t)), list(t == rhs)))
  }

  # Neither side has variables. One side must be affine in parameters.
  lhs_affine <- FALSE
  rhs_affine <- FALSE
  if(dpp_scope()) {
    lhs_affine <- is_affine(lhs)
    rhs_affine <- is_affine(rhs)
  }
  if(!(lhs_affine || rhs_affine))
    stop("Either lhs or rhs must be affine in parameters")

  if(lhs_affine) {
    t <- new("Variable", dim = dim(rhs))
    return(list(lhs %*% t, list(t == rhs)))
  } else {
    t <- new("Variable", dim = dim(lhs))
    return(list(t %*% rhs, list(t == lhs)))
  }
}
