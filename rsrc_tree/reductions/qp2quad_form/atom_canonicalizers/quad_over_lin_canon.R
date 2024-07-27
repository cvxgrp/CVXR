## CVXPY SOURCE: cvxpy/reductions/qp2quad_form/atom_canonicalizers/quad_over_lin_canon.py
Qp2QuadForm.quad_over_lin_canon <- function(expr, args) {
  affine_expr <- args[[1]]
  y <- args[[2]]
  # Simplify if y has no parameters
  if(length(parameters(y)) == 0) {
    affine_expr_size <- size(affine_expr)
    speye <- sparseMatrix(1:affine_expr_size, 1:affine_expr_size, x = rep(1, affine_expr_size))
    quad_mat <- speye/value(y)
  } else {
    # TODO: This code path produces an intermediate dense matrix, but it should be sparse the whole time.
    affine_expr_size <- size(affine_expr)
    speye <- sparseMatrix(1:affine_expr_size, 1:affine_expr_size, x = rep(1, affine_expr_size))
    quad_mat <- speye/y
  }

  if(is(affine_expr, "Variable"))
    return(list(SymbolicQuadForm(affine_expr, quad_mat, expr), list()))
  else {
    # t <- Variable(dim(affine_expr))
    t <- new("Variable", dim = dim(affine_expr))
    return(list(SymbolicQuadForm(t, quad_mat, expr), list(affine_expr == t)))
  }
}
