#'
#' The Qp2SymbolicQp class.
#'
#' This class reduces a quadratic problem to a problem that consists of affine
#' expressions and symbolic quadratic forms.
#'
#' @rdname Qp2SymbolicQp-class
.Qp2SymbolicQp <- setClass("Qp2SymbolicQp", contains = "Canonicalization")
Qp2SymbolicQp <- function(problem = NULL) { .Qp2SymbolicQp(problem = problem) }

setMethod("initialize", "Qp2SymbolicQp", function(.Object, ...) {
  callNextMethod(.Object, ..., canon_methods = Qp2QuadForm.CANON_METHODS)
})

Qp2SymbolicQp.accepts <- function(problem) {
  is_qpwa(expr(problem@objective)) &&
  length(intersect(c("PSD", "NSD"), convex_attributes(variables(problem)))) == 0 &&
  all(sapply(problem@constraints, function(c) {
        (inherits(c, c("NonPosConstraint", "IneqConstraint")) && is_pwl(expr(c))) ||
        (inherits(c, c("ZeroConstraint", "EqConstraint")) && are_args_affine(list(c)))
  }))
}

# Problems with quadratic, piecewise affine objectives, piecewise-linear constraints, inequality constraints,
# and affine equality constraints are accepted.
setMethod("accepts", signature(object = "Qp2SymbolicQp", problem = "Problem"), function(object, problem) {
    Qp2SymbolicQp.accepts(problem)
})

# Converts a QP to an even more symbolic form.
setMethod("perform", signature(object = "Qp2SymbolicQp", problem = "Problem"), function(object, problem) {
  if(!accepts(object, problem))
    stop("Cannot reduce problem to symbolic QP")
  callNextMethod(object, problem)
})

#'
#' The QpMatrixStuffing class.
#'
#' This class fills in numeric values for the problem instance and
#' outputs a DCP-compliant minimization problem with an objective
#' of the form
#'
#' QuadForm(x, p) + t(q) %*% x
#'
#' and Zero/NonPos constraints, both of which exclusively carry
#' affine arguments
#'
#' @rdname QpMatrixStuffing-class
QpMatrixStuffing <- setClass("QpMatrixStuffing", contains = "MatrixStuffing")

setMethod("accepts", signature(object = "QpMatrixStuffing", problem = "Problem"), function(object, problem) {
    inherits(problem@objective, "Minimize") &&
        is_quadratic(problem@objective) &&
        is_dcp(problem) &&
        length(convex_attributes(variables(problem))) == 0 &&
        are_args_affine(problem@constraints) &&
        all(sapply(problem@constraints, inherits, what = c("ZeroConstraint", "NonPosConstraint", "EqConstraint", "IneqConstraint") ))
})

setMethod("stuffed_objective", signature(object = "QpMatrixStuffing", problem = "Problem", extractor = "CoeffExtractor"), function(object, problem, extractor) {
  # Extract to t(x) %*% P %*% x + t(q) %*% x, and store r
  expr <- copy(expr(problem@objective))     # TODO: Need to copy objective?
  Pqr <- coeff_quad_form(extractor, expr)
  P <- Pqr[[1]]
  q <- Pqr[[2]]
  r <- Pqr[[3]]

  # Concatenate all variables in one vector
  boolint <- extract_mip_idx(variables(problem))
  boolean <- boolint[[1]]
  integer <- boolint[[2]]
  # x <- Variable(extractor@N, boolean = boolean, integer = integer)
  #   new_obj <- quad_form(x, P) + t(q) %*% x
  x <- Variable(extractor@N, 1, boolean = boolean, integer = integer)
  new_obj <- new("QuadForm", x = x, P = P) + t(q) %*% x
  return(list(new_obj, x, r))
})

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

Qp2QuadForm.power_canon <- function(expr, args) {
  affine_expr <- args[[1]]
  p <- expr@p
  if(is_constant(expr))
    return(list(Constant(value(expr)), list()))
  else if(p == 0)
    return(list(matrix(1, nrow = nrow(affine_expr), ncol = ncol(affine_expr)), list()))
  else if(p == 1)
    return(list(affine_expr, list()))
  else if(p == 2) {
    if(is(affine_expr, "Variable"))
      return(list(SymbolicQuadForm(affine_expr, make_sparse_diagonal_matrix(size(affine_expr)), expr), list()))
    else {
      # t <- Variable(dim(affine_expr))
      t <- new("Variable", dim = dim(affine_expr))
      return(list(SymbolicQuadForm(t, make_sparse_diagonal_matrix(size(t)), expr), list(affine_expr == t)))
    }
  }
  stop("Non-constant quadratic forms cannot be raised to a power greater than 2.")
}

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

Qp2QuadForm.quad_over_lin_canon <- function(expr, args) {
  affine_expr <- args[[1]]
  y <- args[[2]]
  if(is(affine_expr, "Variable"))
    return(list(SymbolicQuadForm(affine_expr, diag(size(affine_expr))/y, expr), list()))
  else {
    # t <- Variable(dim(affine_expr))
    t <- new("Variable", dim = dim(affine_expr))
    return(list(SymbolicQuadForm(t, diag(size(affine_expr))/y, expr), list(affine_expr == t)))
  }
}

Qp2QuadForm.CANON_METHODS <- list(
  # Reuse cone canonicalization methods.
  Abs = Dcp2Cone.CANON_METHODS$Abs,
  CumSum = Dcp2Cone.CANON_METHODS$CumSum,
  MaxElemwise = Dcp2Cone.CANON_METHODS$MaxElemwise,
  MinElemwise = Dcp2Cone.CANON_METHODS$MinElemwise,
  SumLargest = Dcp2Cone.CANON_METHODS$SumLargest,
  MaxEntries = Dcp2Cone.CANON_METHODS$MaxEntries,
  MinEntries = Dcp2Cone.CANON_METHODS$MinEntries,
  Norm1 = Dcp2Cone.CANON_METHODS$Norm1,
  NormInf = Dcp2Cone.CANON_METHODS$NormInf,
  Indicator = Dcp2Cone.CANON_METHODS$Indicator,
  SpecialIndex = Dcp2Cone.CANON_METHODS$SpecialIndex,
  
  # Canonicalizations that are different for QPs.
  QuadOverLin = Qp2QuadForm.quad_over_lin_canon,
  Power = Qp2QuadForm.power_canon,
  Huber = Qp2QuadForm.huber_canon,
  QuadForm = Qp2QuadForm.quad_form_canon)
