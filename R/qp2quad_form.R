#'
#' The Qp2SymbolicQp class.
#'
#' This class reduces a quadratic problem to a problem that consists of affine
#' expressions and symbolic quadratic forms.
#'
#' @rdname Qp2SymbolicQp-class
setClass("Qp2SymbolicQp", contains = "Canonicalization")

Qp2SymbolicQp.accepts <- function(problem) {
  is_qpwa(problem@objective@expr) &&
    length(intersect(c("PSD", "NSD"), convex_attributes(variables(problem)))) == 0 &&
    all(sapply(problem@constraints, function(c) {
      ((class(c) == "NonPosConstraint" || class(c) == "IneqConstraint")) && is_pwl(c@args[[1]])) ||
        ((class(c) == "ZeroConstraint" || class(c) == "EqConstraint") && are_args_affine(list(c)))
    }))
}

# Problems with quadratic, piecewise affine objectives, piecewise-linear constraints, inequality constraints,
# and affine equality constraints are accepted.
setMethod("accepts", signature(object = "Qp2SymbolicQp", problem = "Problem"), function(object, problem) {
    Qp2SymbolicQp.accepts(problem)
})

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
setClass("QpMatrixStuffing", contains = "MatrixStuffing")

setMethod("accepts", signature(object = "QpMatrixStuffing", problem = "Problem"), function(object, problem) {
    class(problem@objective) == "Minimize" &&
        is_quadratic(problem@objective) &&
        is_dcp(problem) &&
        length(convex_attributes(variables(problem))) == 0 &&
        are_args_affine(problem@constraints) &&
        all(sapply(problem@constraints, function(c) { class(c) %in% c("ZeroConstraint", "NonPosConstraint", "EqConstraint", "IneqConstraint") }))
})

setMethod("stuffed_objective", signature(object = "QpMatrixStuffing", problem = "Problem", extractor = "CoeffExtractor"), function(object, problem, extractor) {
  # Extract to t(x) %*% P %*% x + t(q) %*% x, and store r
  # TODO: Need to copy objective?
  Pqr <- quad_form(extractor, problem_copy@objective@expr)
  P <- Pqr[[1]]
  q <- Pqr[[2]]
  r <- Pqr[[3]]

  # Concatenate all variables in one vector
  boolint <- extract_mip_idx(variables(problem))
  boolean <- boolint[[1]]
  integer <- boolint[[2]]
  x <- Variable(extractor@N, boolean = boolean, integer = integer)
  new_obj <- QuadForm(x, P) + t(q) %*% x

  inverse_data@r <- r
  return(list(new_obj, x, r))
})

# Converts a QP to an even more symbolic form.
setMethod("perform", signature(object = "Qp2SymbolicQp", problem = "Problem"), function(object, problem) {
  if(!(accepts(object, problem)))
    stop("Cannot reduce problem to symbolic QP")
  return(perform(Canonicalization(qp_canon_methods), problem))
})

# Atom canonicalizers
Qp2QuadForm.CANON_METHODS <- list(# Reuse cone canonicalization methods.
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
                                  QuadOverLin = quad_over_lin_canon,
                                  Power = power_canon,
                                  Huber = huber_canon,
                                  QuadForm = quad_form_canon)

huber_canon <- function(expr, args) {
  M <- expr@M
  x <- args[[1]]
  shape <- shape(expr)
  n <- Variable(shape)
  s <- Variable(shape)
  
  # n^2 + 2*M*|s|
  # TODO: Make use of recursion inherent to canonicalization process and just return a power / abs expression for readability's sake.
  power_expr <- Power(n, 2)
  canon <- power_canon(power_expr, power_expr@args)
  n2 <- canon[[1]]
  constr_sq <- canon[[2]]
  
  abs_expr <- abs(s)
  canon <- abs_canon(abs_expr, abs_expr@args)
  abs_s <- canon[[1]]
  constr_abs <- canon[[2]]
  obj <- n2 + 2*M*abs_s
  
  constraints <- c(constr_sq, constr_abs)
  constraints <- c(constraints, list(x == s + n))
  return(list(obj, constraints))
}

power_canon <- function(expr, args) {
  affine_expr <- args[[1]]
  p <- expr@p
  if(is_constant(expr))
    return(list(Constant(value(expr)), list()))
  else if(p == 0)
    return(list(matrix(1, nrow = shape(affine_expr)[1], ncol = shape(affine_expr)[2]), list()))
  else if(p == 1)
    return(list(affine_expr, list()))
  else if(p == 2) {
    if(is(affine_expr, "Variable"))
      return(list(SymbolicQuadForm(affine_expr, diag(size(affine_expr)), expr), list()))
    else {
      t <- Variable(shape(affine_expr))
      return(list(SymbolicQuadForm(t, diag(size(t)), expr), list(affine_expr == t)))
    }
  }
  stop("Non-constant quadratic forms cannot be raised to a power greater than 2.")
}

quad_form_canon <- function(expr, args) {
  affine_expr <- expr@args[[1]]
  P <- expr@args[[2]]
  if(is(affine_expr, "Variable"))
    return(list(SymbolicQuadForm(affine_expr, P, expr), list()))
  else {
    t <- Variable(shape(affine_expr))
    return(list(SymbolicQuadForm(t, P, expr), list(affine_expr == t)))
  }
}

quad_over_lin_canon <- function(expr, args) {
  affine_expr <- args[[1]]
  y <- args[[2]]
  if(is(affine_expr, "Variable"))
    return(list(SymbolicQuadForm(affine_expr, diag(size(affine_expr))/y, expr), list()))
  else {
    t <- Variable(shape(affine_expr))
    return(list(SymbolicQuadForm(t, diag(size(affine_expr)/y), expr), list(affine_expr == t)))
  }
}
