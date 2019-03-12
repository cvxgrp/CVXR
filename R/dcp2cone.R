#'
#' Reduce DCP Problem to Conic Form
#'
#' This reduction takes as input (minimization) DCP problems and converts them into problems
#' with affine objectives and conic constraints whose arguments are affine.
#'
#' @rdname Dcp2Cone-class
.Dcp2Cone <- setClass("Dcp2Cone", contains = "Canonicalization")
Dcp2Cone <- function(problem = NULL) { .Dcp2Cone(problem = problem) }

setMethod("initialize", "Dcp2Cone", function(.Object, ...) {
  callNextMethod(.Object, ..., canon_methods = Dcp2Cone.CANON_METHODS)
})

setMethod("accepts", signature(object = "Dcp2Cone", problem = "Problem"), function(object, problem) {
  class(problem@objective) == "Minimize" && is_dcp(problem)
})

setMethod("perform", signature(object = "Dcp2Cone", problem = "Problem"), function(object, problem) {
  if(!accepts(object, problem))
    stop("Cannot reduce problem to cone program")
  callNextMethod(object, problem)
})

#'
#' Construct Matrices for Linear Cone Problems
#'
#' Linear cone problems are assumed to have a linear objective and cone constraints,
#' which may have zero or more arguments, all of which must be affine.
#'
#' minimize c^Tx
#' subject to cone_constr1(A_1*x + b_1, ...)
#'            ...
#'            cone_constrK(A_K*x + b_K, ...)
#'
#' @rdname ConeMatrixStuffing-class
ConeMatrixStuffing <- setClass("ConeMatrixStuffing", contains = "MatrixStuffing")

setMethod("accepts", signature(object = "ConeMatrixStuffing", problem = "Problem"), function(object, problem) {
    class(problem@objective) == "Minimize" &&
        is_affine(problem@objective@expr) &&
        length(convex_attributes(variables(problem))) == 0 &&
        are_args_affine(problem@constraints)
})

setMethod("stuffed_objective", signature(object = "ConeMatrixStuffing", problem = "Problem", extractor = "CoeffExtractor"), function(object, problem, extractor) {
  # Extract to t(c) %*% x, store in r
  CR <- affine(extractor, problem@objective@expr)
  C <- CR[[1]]
  R <- CR[[2]]

  c <- matrix(C, ncol = 1)   # TODO: Check if converted to dense matrix and flattened like in CVXPY
  boolint <- extract_mip_idx(variables(problem))
  boolean <- boolint[[1]]
  integer <- boolint[[2]]
  x <- Variable(extractor@N, boolean = boolean, integer = integer)

  new_obj <- t(c) %*% x + 0

  return(list(new_obj, x, R[[1]]))
})

# Atom canonicalizers.
Dcp2Cone.cumsum_canon <- function(expr, args) {
  X <- args[[1]]
  axis <- expr@axis

  # Implicit O(n) definition:
  # X = Y[1,:] - Y[2:nrow(Y),:]
  Y <- Variable(dim(expr))
  if(axis == 1)   # Cumulative sum on each column
    constr <- list(X[2:nrow(X),] == Y[2:nrow(Y),] - Y[1:(nrow(Y)-1),], Y[1,] == X[1,])   # TODO: Check corner cases
  else   # Cumulative sum on each row
    constr <- list(X[,2:ncol(X)] == Y[,2:ncol(Y)] - Y[,1:(ncol(Y)-1)], Y[,1] == X[,1])
  return(list(Y, constr))
}

Dcp2Cone.entr_canon <- function(expr, args) {
  x <- args[[1]]
  expr_dim <- dim(expr)
  t <- Variable(expr_dim)

  # -x*log(x) >= t is equivalent to x/exp(t/x) <= 1
  # TODO: ExpCone requires each of its inputs to be a Variable; is this something we want to change?
  ones <- Constant(matrix(1, nrow = expr_dim[1], ncol = expr_dim[2]))
  constraints <- list(ExpCone(t, x, ones))
  return(list(t, constraints))
}

Dcp2Cone.exp_canon <- function(expr, args) {
  expr_dim <- dim(expr)
  x <- promote(args[[1]], expr_dim)
  t <- Variable(expr_dim)
  ones <- Constant(matrix(1, nrow = expr_dim[1], ncol = expr_dim[2]))
  constraints <- list(ExpCone(x, ones, t))
  return(list(t, constraints))
}

Dcp2Cone.geo_mean_canon <- function(expr, args) {
  x <- args[[1]]
  w <- expr@w
  expr_dim <- dim(expr)
  t <- Variable(expr_dim)

  x_list <- lapply(1:length(w), function(i) { x[i] })

  # TODO: Catch cases where we have (0,0,1)?
  # TODO: What about curvature case (should be affine) in trivial case of (0,0,1)?
  # Should this behavior match with what we do in power?
  return(list(t, gm_constrs(t, x_list, w)))
}

Dcp2Cone.huber_canon <- function(expr, args) {
  M <- expr@M
  x <- args[[1]]
  expr_dim <- dim(expr)
  n <- Variable(expr_dim)
  s <- Variable(expr_dim)

  # n^2 + 2*M*|s|
  # TODO: Make use of recursion inherent to canonicalization process and just return a
  # power/abs expression for readiability's sake
  power_expr <- power(n,2)
  canon <- Dcp2Cone.power_canon(power_expr, power_expr@args)
  n2 <- canon[[1]]
  constr_sq <- canon[[2]]

  abs_expr <- abs(s)
  canon <- EliminatePwl.abs_canon(abs_expr, abs_expr@args)
  abs_s <- canon[[1]]
  constr_abs <- canon[[2]]

  obj <- n2 + 2*M*abs_s

  # x == s + n
  constraints <- c(constr_sq, constr_abs)
  constraints <- c(constraints, x == s + n)
  return(list(obj, constraints))
}

Dcp2Cone.indicator_canon <- function(expr, args) {
  return(list(0, args))
}

Dcp2Cone.kl_div_canon <- function(expr, args) {
  expr_dim <- dim(expr)
  x <- promote(args[[1]], expr_dim)
  y <- promote(args[[2]], expr_dim)
  t <- Variable(expr_dim)
  constraints <- list(ExpCone(t, x, y), y >= 0)
  obj <- y - x - t
  return(list(obj, constraints))
}

Dcp2Cone.lambda_max_canon <- function(expr, args) {
  A <- args[[1]]
  n <- nrow(A)
  t <- Variable()
  prom_t <- promote(t, c(n,1))
  # Constraint I*t - A to be PSD; note this expression must be symmetric
  tmp_expr <- DiagVec(prom_t) - A
  constr <- list(tmp_expr == t(tmp_expr), PSDConstraint(tmp_expr))
  return(list(t, constr))
}

Dcp2Cone.lambda_sum_largest_canon <- function(expr, args) {
  # S_k(X) denotes lambda_sum_largest(X, k)
  # t >= k S_k(X - Z) + trace(Z), Z is PSD
  # implies
  # t >= ks + trace(Z)
  # Z is PSD
  # sI >= X - Z (PSD sense)
  # which implies
  # t >= ks + trace(Z) >= S_k(sI + Z) >= S_k(X)
  # We use the fact that
  # S_k(X) = sup_{sets of k orthonormal vectors u_i}\sum_{i}u_i^T X u_i
  # and if Z >= X in PSD sense then
  # \sum_{i}u_i^T Z u_i >= \sum_{i}u_i^T X u_i
  #
  # We have equality when s = lambda_k and Z diagonal
  #  with Z_{ii} = (lambda_i - lambda_k)_+

  X <- expr@args[[1]]
  k <- expr@k
  Z <- Variable(c(nrow(X), nrow(X)), PSD = TRUE)
  canon <- Dcp2Cone.lambda_max_canon(expr, list(X - Z))
  obj <- canon[[1]]
  constr <- canon[[2]]
  obj <- k*obj + trace(Z)
  return(list(obj, constr))
}

Dcp2Cone.log1p_canon <- function(expr, args) {
  return(Dcp2Cone.log_canon(expr, list(args[[1]] + 1)))
}

Dcp2Cone.log_canon <- function(expr, args) {
  x <- args[[1]]
  expr_dim <- dim(expr)
  t <- Variable(expr_dim)
  ones <- Constant(matrix(1, nrow = expr_dim[1], ncol = expr_dim[2]))
  # TODO: ExpCone requires each of its inputs to be a Variable; is this something that we want to change?
  constraints <- list(ExpCone(t, ones, x))
  return(list(t, constraints))
}

Dcp2Cone.log_det_canon <- function(expr, args) {
  # Reduces the atom to an affine expression and list of constraints.
  #
  # Creates the equivalent problem::
  #
  # maximize    sum(log(D[i, i]))
  # subject to: D diagonal
  # diag(D) = diag(Z)
  # Z is upper triangular.
  # [D Z; t(Z) A] is positive semidefinite
  #
  # The problem computes the LDL factorization:
  #
  # A = (Z^TD^{-1})D(D^{-1}Z)
  #
  # This follows from the inequality:
  #
  # \det(A) >= \det(D) + \det([D, Z; Z^T, A])/\det(D) >= \det(D)
  #
  # because (Z^TD^{-1})D(D^{-1}Z) is a feasible D, Z that achieves
  # det(A) = det(D) and the objective maximizes det(D).
  #
  # Parameters
  # ----------
  # expr : log_det
  # args : list of arguments for the expression
  #
  # Returns
  # -------
  # (Variable for objective, list of constraints)

  A <- args[[1]]   # n by n matrix
  n <- nrow(A)
  # Require that X and A are PSD.
  X <- Variable(c(2*n, 2*n), PSD = TRUE)
  constraints <- list(PSDConstraint(A))

  # Fix Z as upper triangular
  # TODO: Represent Z as upper triangular vector
  Z <- Variable(c(n,n))
  Z_lower_tri <- upper_tri(t(Z))
  constraints <- list(Z_lower_tri == 0)

  # Fix diag(D) = Diag(Z): D[i,i] = Z[i,i]
  D <- Variable(n)
  constraints <- c(constraints, D == DiagMat(Z))
  # Fix X using the fact that A must be affine by the DCP rules
  # X[1:n, 1:n] == D
  constraints <- c(constraints, X[1:n, 1:n] == DiagVec(D))
  # X[1:n, (n+1):(2*n)] == Z
  constraints <- c(constraints, X[1:n, (n+1):(2*n)] == Z)
  # X[(n+1):(2*n),  (n+1):(2*n)] == A
  constraints <- c(constraints, X[(n+1):(2*n), (n+1):(2*n)] == A)
  # Add the objective sum(log(D[i,i]))
  log_expr <- log(D)
  canon <- Dcp2Cone.log_canon(log_expr, log_expr@args)
  obj <- canon[[1]]
  constr <- canon[[2]]
  constraints <- c(constraints, constr)
  return(list(sum(obj), constraints))
}

Dcp2Cone.log_sum_exp_canon <- function(expr, args) {
  x <- args[[1]]
  x_dim <- dim(x)
  expr_dim <- dim(expr)
  axis <- expr@axis
  keepdims <- expr@keepdims
  t <- Variable(expr_dim)
  
  # log(sum(exp(x))) <= t <=> sum(exp(x-t)) <= 1.
  if(is.na(axis))   # shape = c(1,1)
    promoted_t <- promote(t, x_dim)
  else if(axis == 2)   # shape = c(1,n)
    promoted_t <- Constant(matrix(1, nrow = x_dim[1], ncol = 1) %*% reshape_expr(t, c(1 + x_dim[2], x_dim[3:length(x_dim)])))
  else   # shape = c(m,1)
    promoted_t <- reshape_expr(t, c(1 + x_dim[1], x_dim[2:(length(x_dim)-1)])) %*% Constant(matrix(1, nrow = 1, ncol = x_dim[2]))
  
  exp_expr <- exp(x - promoted_t)
  canon <- Dcp2Cone.exp_canon(exp_expr, exp_expr@args)
  obj <- sum_entries(canon[[1]], axis = axis, keepdims = keepdims)
  ones <- Constant(matrix(1, nrow = expr_dim[1], ncol = expr_dim[2]))
  constraints <- c(canon[[2]], obj <= ones)
  return(list(t, constraints))
}

Dcp2Cone.logistic_canon <- function(expr, args) {
  x <- args[[1]]
  expr_dim <- dim(expr)
  # log(1 + exp(x)) <= t is equivalent to exp(-t) + exp(x - t) <= 1
  t0 <- Variable(expr_dim)
  canon1 <- Dcp2Cone.exp_canon(expr, list(-t0))
  canon2 <- Dcp2Cone.exp_canon(expr, list(x - t0))

  t1 <- canon1[[1]]
  constr1 <- canon1[[2]]
  t2 <- canon2[[1]]
  constr2 <- canon2[[2]]

  ones <- Constant(matrix(1, nrow = expr_dim[1], ncol = expr_dim[2]))
  constraints <- c(constr1, constr2, list(t1 + t2 <= ones))
  return(list(t0, constraints))
}

Dcp2Cone.matrix_frac_canon <- function(expr, args) {
  X <- args[[1]]   # n by m matrix
  P <- args[[2]]   # n by n matrix

  if(length(dim(X)) == 1)
    X <- reshape_expr(X, c(nrow(X), 1))
  X_dim <- dim(X)
  n <- X_dim[1]
  m <- X_dim[2]

  # Create a matrix with Schur complement Tvar - t(X) %*% inv(P) %*% X
  M <- Variable(c(n+m, n+m), PSD = TRUE)
  Tvar <- Variable(c(m,m), symmetric = TRUE)
  constraints <- list()
  # Fix M using the fact that P must be affine by the DCP rules.
  # M[1:n, 1:n] == P
  constraints <- c(constraints, M[1:n, 1:n] == P)
  # M[1:n, (n+1):(n+m)] == X
  constraints <- c(constraints, M[1:n, (n+1):(n+m)] == X)
  # M[(n+1):(n+m), (n+1):(n+m)] == Tvar
  constraints <- c(constraints, M[(n+1):(n+m), (n+1):(n+m)] == Tvar)
  return(list(trace(Tvar), constraints))
}

Dcp2Cone.normNuc_canon <- function(expr, args) {
  A <- args[[1]]
  A_dim <- dim(A)
  m <- A_dim[1]
  n <- A_dim[2]

  # Create the equivalent problem:
  #   minimize (trace(U) + trace(V))/2
  #   subject to:
  #            [U A; t(A) V] is positive semidefinite
  X <- Variable(c(m+n, m+n), PSD = TRUE)
  constraints <- list()

  # Fix X using the fact that A must be affine by the DCP rules.
  # X[1:rows, (rows+1):(rows+cols)] == A
  constraints <- c(constraints, X[1:m, (m+1):(m+n)] == A)
  trace_value <- 0.5*trace(X)
  return(list(trace_value, constraints))
}

Dcp2Cone.pnorm_canon <- function(expr, args) {
  x <- args[[1]]
  p <- expr@p
  axis <- expr@axis
  expr_dim <- dim(expr)
  t <- Variable(expr_dim)

  if(p == 2) {
    if(is.na(axis)) {
      if(!is.null(expr_dim))
        stop("Dimensions should be NULL")
      return(list(t, list(SOC(t, vec(x)))))
    } else
      return(list(t, list(SOC(vec(t), x, axis))))
  }

  # We need an absolute value constraint for the symmetric convex branches (p > 1)
  constraints <- list()
  if(p > 1) {
    # TODO: Express this more naturally (recursively) in terms of the other atoms
    abs_expr <- abs(x)
    canon <- EliminatePwl.abs_canon(abs_expr, abs_expr@args)
    x <- canon[[1]]
    abs_constraints <- canon[[2]]
    constraints <- c(constraints, abs_constraints)
  }

  # Now, we take care of the remaining convex and concave branches to create the
  # rational powers. We need a new variable, r, and the constraint sum(r) == t
  r <- Variable(dim(x))
  constraints <- c(constraints, list(sum(r) == t))

  # TODO: No need to run gm_constr to form the tree each time.
  # We only need to form the tree once.
  promoted_t <- Constant(matrix(1, nrow = nrow(x), ncol = ncol(x))) %*% t
  p <- gmp::as.bigq(p)
  if(p < 0)
    constraints <- c(constraints, gm_constrs(promoted_t, list(x, r), c(-p/(1-p), 1/(1-p))))
  else if(p > 0 && p < 1)
    constraints <- c(constraints, gm_constrs(r, list(x, promoted_t), c(p, 1-p)))
  else if(p > 1)
    constraints <- c(constraints, gm_constrs(x, list(r, promoted_t), c(1/p, 1-1/p)))
  return(list(t, constraints))
}

Dcp2Cone.power_canon <- function(expr, args) {
  x <- args[[1]]
  p <- expr@p
  w <- expr@w

  if(p == 1)
    return(list(x, list()))

  expr_dim <- dim(expr)
  ones <- Constant(matrix(1, nrow = expr_dim[1], ncol = expr_dim[2]))
  if(p == 0)
    return(list(ones, list()))
  else {
    t <- Variable(expr_dim)
    # TODO: gm_constrs requires each of its inputs to be a Variable; is this something that we want to change?
    if(p > 0 && p < 1)
      return(list(t, gm_constrs(t, list(x, ones), w)))
    else if(p > 1)
      return(list(t, gm_constrs(x, list(t, ones), w)))
    else if(p < 0)
      return(list(t, gm_constrs(ones, list(x, t), w)))
    else
      stop("This power is not yet supported")
  }
}

Dcp2Cone.quad_form_canon <- function(expr, args) {
  decomp <- .decomp_quad(value(args[[2]]))
  scale <- decomp[[1]]
  M1 <- decomp[[2]]
  M2 <- decomp[[3]]

  if(size(M1) > 0)
    expr <- sum_squares(Constant(t(M1)) %*% args[[1]])
  else if(size(M2)> 0) {
    scale <- -scale
    expr <- sum_squares(Constant(t(M2)) %*% args[[1]])
  }
  canon <- Dcp2Cone.quad_over_lin_canon(expr, expr@args)
  obj <- canon[[1]]
  constr <- canon[[2]]
  return(list(scale * obj, constr))
}

Dcp2Cone.quad_over_lin_canon <- function(expr, args) {
  # quad_over_lin := sum_{ij} X^2_{ij} / y
  x <- args[[1]]
  y <- flatten(args[[2]])
  
  # Pre-condition: dim = c()
  t <- Variable(1)
  
  # (y+t, y-t, 2*x) must lie in the second-order cone, where y+t is the scalar part
  # of the second-order cone constraint
  constraints <- list(SOC(t = y+t, X = hstack(list(y-t, 2*flatten(x))), axis = 2))
  return(list(t, constraints))
}

Dcp2Cone.sigma_max_canon <- function(expr, args) {
  A <- args[[1]]
  A_dim <- dim(A)
  n <- A_dim[1]
  m <- A_dim[2]
  X <- Variable(c(n+m, n+m), PSD = TRUE)

  expr_dim <- dim(expr)
  t <- Variable(expr_dim)
  constraints <- list()

  # Fix X using the fact that A must be affine by the DCP rules.
  # X[1:n, 1:n] == I_n*t
  constraints <- c(constraints, X[1:n, 1:n] == Constant(sparseMatrix(i = 1:n, j = 1:n, x = 1) %*% t))

  # X[1:n, (n+1):(n+m)] == A
  constraints <- c(constraints, X[1:n, (n+1):(n+m)] == A)

  # X[(n+1):(n+m), (n+1):(n+m)] == I_m*t
  constraints <- c(constraints, X[(n+1):(n+m), (n+1):(n+m)] == Constant(sparseMatrix(i = 1:m, j = 1:m, x = 1) %*% t))
  return(list(t, constraints))
}

# TODO: Remove pwl canonicalize methods and use EliminatePwl reduction instead.
Dcp2Cone.CANON_METHODS <- list(CumSum = Dcp2Cone.cumsum_canon,
                               GeoMean = Dcp2Cone.geo_mean_canon,
                               LambdaMax = Dcp2Cone.lambda_max_canon,
                               LambdaSumLargest = Dcp2Cone.lambda_sum_largest_canon,
                               LogDet = Dcp2Cone.log_det_canon,
                               LogSumExp = Dcp2Cone.log_sum_exp_canon,
                               MatrixFrac = Dcp2Cone.matrix_frac_canon,
                               MaxEntries = EliminatePwl.CANON_METHODS$MaxEntries,
                               MinEntries = EliminatePwl.CANON_METHODS$MinEntries,
                               Norm1 = EliminatePwl.CANON_METHODS$Norm1,
                               NormNuc = Dcp2Cone.normNuc_canon,
                               NormInf = EliminatePwl.CANON_METHODS$NormInf,
                               Pnorm = Dcp2Cone.pnorm_canon,
                               QuadForm = Dcp2Cone.quad_form_canon,
                               QuadOverLin = Dcp2Cone.quad_over_lin_canon,
                               SigmaMax = Dcp2Cone.sigma_max_canon,
                               SumLargest = EliminatePwl.CANON_METHODS$SumLargest,
                               Abs = EliminatePwl.CANON_METHODS$Abs,
                               Entr = Dcp2Cone.entr_canon,
                               Exp = Dcp2Cone.exp_canon,
                               Huber = Dcp2Cone.huber_canon,
                               KLDiv = Dcp2Cone.kl_div_canon,
                               Log = Dcp2Cone.log_canon,
                               Log1p = Dcp2Cone.log1p_canon,
                               Logistic = Dcp2Cone.logistic_canon,
                               MaxElemwise = EliminatePwl.CANON_METHODS$MaxElemwise,
                               MinElemwise = EliminatePwl.CANON_METHODS$MinElemwise,
                               Power = Dcp2Cone.power_canon,
                               Indicator = Dcp2Cone.indicator_canon,
                               SpecialIndex = special_index_canon)
