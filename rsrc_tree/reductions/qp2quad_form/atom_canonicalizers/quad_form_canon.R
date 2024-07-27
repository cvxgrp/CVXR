## CVXPY SOURCE: cvxpy/reductions/qp2quad_form/atom_canonicalizers/quad_form_canon.py

Qp2QuadForm.quad_form_canon <- function(expr, args) {
  affine_expr <- expr@args[[1]]
  P <- expr@args[[2]]
  if(is(affine_expr, "Variable"))
    return(list(SymbolicQuadForm(affine_expr, P, expr), list()))
  else {
    # t <- Variable(dim(affine_expr))
    t <- new("Variable", dim = dim(affine_expr))
    return(list(SymbolicQuadForm(t, P, expr), list(affine_expr == t)))
  }
}

