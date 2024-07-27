## CVXPY SOURCE: cvxpy/reductions/qp2quad_form/atom_canonicalizers/power_canon.py

Qp2QuadForm.power_canon <- function(expr, args) {
  affine_expr <- args[[1]]
  if(is.numeric(expr@p))
    p <- expr@p
  else
    p <- value(expr@p)

  if(is_constant(expr))
    return(list(Constant(value(expr)), list()))
  else if(p == 0)
    return(list(matrix(1, nrow = nrow(affine_expr), ncol = ncol(affine_expr)), list()))
  else if(p == 1)
    return(list(affine_expr, list()))
  else if(p == 2) {
    if(is(affine_expr, "Variable")) {
      affine_expr_size <- size(affine_expr)
      speye <- sparseMatrix(1:affine_expr_size, 1:affine_expr_size, x = rep(1, affine_expr_size))
      return(list(SymbolicQuadForm(affine_expr, speye, expr), list()))
    } else {
      # t <- Variable(dim(affine_expr))
      t <- new("Variable", dim = dim(affine_expr))
      t_size <- size(t)
      speye <- sparseMatrix(1:t_size, 1:t_size, x = rep(1, t_size))
      return(list(SymbolicQuadForm(t, speye, expr), list(affine_expr == t)))
    }
  }
  stop("Non-constant quadratic forms cannot be raised to a power greater than 2.")
}
