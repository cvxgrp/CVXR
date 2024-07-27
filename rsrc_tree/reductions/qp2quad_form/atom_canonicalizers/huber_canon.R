## CVXPY SOURCE: cvxpy/reductions/qp2quad_form/atom_canonicalizers/huber_canon.py
# Atom canonicalizers
Qp2QuadForm.huber_canon <- function(expr, args) {
  M <- expr@M
  x <- args[[1]]
  expr_dim <- dim(expr)
  # n <- Variable(expr_dim)
  # s <- Variable(expr_dim)
  n <- new("Variable", dim = expr_dim)
  s <- new("Variable", dim = expr_dim)

  # n^2 + 2*M*|s|
  # TODO: Make use of recursion inherent to canonicalization process and just return a power / abs expression for readability's sake.
  power_expr <- Power(n, 2)
  canon <- Qp2QuadForm.power_canon(power_expr, power_expr@args)
  n2 <- canon[[1]]
  constr_sq <- canon[[2]]

  abs_expr <- abs(s)
  canon <- EliminatePwl.abs_canon(abs_expr, abs_expr@args)
  abs_s <- canon[[1]]
  constr_abs <- canon[[2]]
  obj <- n2 + 2*M*abs_s

  constraints <- c(constr_sq, constr_abs)
  constraints <- c(constraints, list(x == s + n))
  return(list(obj, constraints))
}
