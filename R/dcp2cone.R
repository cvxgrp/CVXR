setClass("Dcp2Cone", contains = "Canonicalization")

setMethod("accepts", signature(object = "Dcp2Cone", problem = "Problem"), function(object, problem) {
  class(problem@objective) == "Minimize" && is_dcp(problem)
})

setMethod("apply", signature(object = "Dcp2Cone", problem = "Problem"), function(object, problem) {
  if(!accepts(object, problem))
    stop("Cannot reduce problem to cone program")
  callNextMethod(object, problem)
})

setClass("ConeMatrixStuffing", contains = "MatrixStuffing")

setMethod("accepts", signature(object = "ConeMatrixStuffing", problem = "Problem"), function(object, problem) {
    class(problem@objective) == "Minimize" &&
        is_affine(problem@objective@expr) &&
        length(convex_attributes(variables(problem))) == 0 &&
        are_args_affine(problem@constraints)
})

setMethod("stuffed_objective", signature(object = "ConeMatrixStuffing", problem = "Problem", inverse_data = "InverseData"), function(object, problem, inverse_data) {
  extractor <- CoeffExtractor(inverse_data)

  # Extract to t(c) %*% x, store in r
  CR <- get_coeffs(extractor, problem@objective@expr)
  C <- CR[[1]]
  R <- CR[[2]]

  c <- matrix(C, ncol = 1)   # TODO: Check if converted to dense matrix and flattened like in CVXPY
  boolint <- extract_mip_idx(variables(problem))
  boolean <- boolint[[1]]
  integer <- boolint[[2]]
  x <- Variable(inverse_data@x_length, boolean = boolean, integer = integer)

  new_obj <- t(c) %*% x + 0

  inverse_data@r <- R[[1]]
  return(list(new_obj, x))
})

cumsum_canon <- function(expr, args) {
  X <- args[[1]]
  axis <- expr@axis

  # Implicit O(n) definition:
  # X = Y[1,:] - Y[2:nrow(Y),:]
  Y <- Variable(shape(expr))
  if(axis == 1)   # Cumulative sum on each column
    constr <- list(X[2:nrow(X),] == Y[2:nrow(Y),] - Y[1:(nrow(Y)-1),], Y[1,] == X[1,])   # TODO: Check corner cases
  else   # Cumulative sum on each row
    constr <- list(X[,2:ncol(X)] == Y[,2:ncol(Y)] - Y[,1:(ncol(Y)-1)], Y[,1] == X[,1])
  return(list(Y, constr))
}

entr_canon <- function(expr, args) {
  x <- args[[1]]
  shape <- shape(expr)
  t <- Variable(shape)

  # -x*log(x) >= t is equivalent to x/exp(t/x) <= 1
  # TODO: ExpCone requires each of its inputs to be a Variable; is this something we want to change?
  ones <- Constant(matrix(1, nrow = shape[1], ncol = shape[2]))
  constraints <- list(ExpCone(t, x, ones))
  return(list(t, constraints))
}

exp_canon <- function(expr, args) {
  shape <- shape(expr)
  x <- promote(args[[1]], shape)
  t <- Variable(shape)
  ones <- Constant(matrix(1, nrow = shape[1], ncol = shape[2]))
  constraints <- list(ExpCone(x, ones, t))
  return(list(t, constraints))
}

geo_mean_canon <- function(expr, args) {
  x <- args[[1]]
  w <- expr@w
  shape <- shape(expr)
  t <- Variable(shape)

  x_list <- lapply(1:length(w), function(i) { x[i] })

  # TODO: Catch cases where we have (0,0,1)?
  # TODO: What about curvature case (should be affine) in trivial case of (0,0,1)?
  # Should this behavior match with what we do in power?
  return(list(t, gm_constrs(t, x_list, w)))
}

huber_canon <- function(expr, args) {
  M <- expr@M
  x <- args[[1]]
  shape <- shape(expr)
  n <- Variable(shape)
  s <- Variable(shape)

  # n^2 + 2*M*|s|
  # TODO: Make use of recursion inherent to canonicalization process and just return a
  # power/abs expression for readiability's sake
  power_expr <- power(n,2)
  canon <- power_canon(power_expr, power_expr@args)
  n2 <- canon[[1]]
  constr_sq <- canon[[2]]

  abs_expr <- abs(s)
  canon <- abs_canon(abs_expr, abs_expr@args)
  abs_s <- canon[[1]]
  constr_abs <- canon[[2]]

  obj <- n2 + 2*M*abs_s

  # x == s + n
  constraints <- c(constr_sq, constr_abs)
  constraints <- c(constraints, x == s + n)
  return(list(obj, constraints))
}

indicator_canon <- function(expr, args) {
  return(list(0, args))
}

kl_div_canon <- function(expr, args) {
  shape <- shape(expr)
  x <- promote(args[[1]], shape)
  y <- promote(args[[2]], shape)
  t <- Variable(shape)
  constraints <- list(ExpCone(t, x, y), y >= 0)
  obj <- y - x - t
  return(list(obj, constraints))
}

lambda_max_canon <- function(expr, args) {
  A <- args[[1]]
  n <- shape(A)[1]
  t <- Variable()
  prom_t <- promote(t, c(n,1))
  # Constraint I*t - A to be PSD; note this expression must be symmetric
  tmp_expr <- diag_vec(prom_t) - A
  constr <- list(tmp_expr == t(tmp_expr), PSD(tmp_expr))
  return(list(t, constr))
}

lambda_sum_largest_canon <- function(expr, args) {
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
  Z <- Variable(c(shape(X)[1], shape(X)[1]), PSD = TRUE)
  canon <- lambda_max_canon(expr, list(X - Z))
  obj <- canon[[1]]
  constr <- canon[[2]]
  obj <- k*obj + trace(Z)
  return(list(obj, constr))
}

log1p_canon <- function(expr, args) {
  return(log_canon(expr, list(args[[1]] + 1)))
}

log_canon <- function(expr, args) {
  x <- args[[1]]
  shape <- shape(expr)
  t <- Variable(shape)
  ones <- Constant(matrix(1, nrow = shape[1], ncol = shape[2]))
  # TODO: ExpCone requires each of its inputs to be a Variable; is this something that we want to change?
  constraints <- list(ExpCone(t, ones, x))
  return(list(t, constraints))
}

log_det_canon <- function(expr, args) {
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
  n <- shape(A)[1]
  # Require that X and A are PSD.
  X <- Variable(c(2*n, 2*n), PSD = TRUE)
  constraints <- list(PSD(A))

  # Fix Z as upper triangular
  # TODO: Represent Z as upper triangular vector
  Z <- Variable(c(n,n))
  Z_lower_tri <- upper_tri(transpose(Z))
  constraints <- list(Z_lower_tri == 0)

  # Fix diag(D) = Diag(Z): D[i,i] = Z[i,i]
  D <- Variable(n)
  constraints <- c(constraints, D == diag_mat(Z))
  # Fix X using the fact that A must be affine by the DCP rules
  # X[1:n, 1:n] == D
  constraints <- c(constraints, X[1:n, 1:n] == diag_vec(D))
  # X[1:n, (n+1):(2*n)] == Z
  constraints <- c(constraints, X[1:n, (n+1):(2*n)] == Z)
  # X[(n+1):(2*n),  (n+1):(2*n)] == A
  constraints <- c(constraints, X[(n+1):(2*n), (n+1):(2*n)] == A)
  # Add the objective sum(log(D[i,i]))
  log_expr <- log(D)
  canon <- log_canon(log_expr, log_expr@args)
  obj <- canon[[1]]
  constr <- canon[[2]]
  constraints <- c(constraints, constr)
  return(list(sum(obj), constraints))
}

logistic_canon <- function(expr, args) {
  x <- args[[1]]
  shape <- shape(expr)
  # log(1 + exp(x)) <= t is equivalent to exp(-t) + exp(x - t) <= 1
  t0 <- Variable(shape)
  canon1 <- exp_canon(expr, list(-t0))
  canon2 <- exp_canon(expr, list(x - t0))

  t1 <- canon1[[1]]
  constr1 <- canon1[[2]]
  t2 <- canon2[[1]]
  constr2 <- canon2[[2]]

  ones <- Constant(matrix(1, nrow = shape[1], ncol = shape[2]))
  constraints <- c(constr1, constr2, list(t1 + t2 <= ones))
  return(list(t0, constraints))
}

matrix_frac_canon <- function(expr, args) {
  X <- args[[1]]   # n by m matrix
  P <- args[[2]]   # n by n matrix

  if(length(shape(X)) == 1)
    X <- reshape(X, c(shape(X)[1], 1))
  shape <- shape(X)
  n <- shape[1]
  m <- shape[2]

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

normNuc_canon <- function(expr, args) {
  A <- args[[1]]
  shape <- shape(A)
  m <- shape[1]
  n <- shape[2]

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

pnorm_canon <- function(expr, args) {
  x <- args[[1]]
  p <- expr@p
  axis <- expr@axis
  shape <- shape(expr)
  t <- Variable(shape)

  if(p == 2) {
    if(is.na(axis)) {
      if(!is.null(shape))
        stop("shape should be NULL")
      return(list(t, list(SOC(t, vec(x)))))
    } else
      return(list(t, list(SOC(vec(t), x, axis))))
  }

  # We need an absolute value constraint for the symmetric convex branches (p > 1)
  constraints <- list()
  if(p > 1) {
    # TODO: Express this more naturally (recursively) in terms of the other atoms
    abs_expr <- abs(x)
    canon <- abs_canon(abs_expr, abs_expr@args)
    x <- canon[[1]]
    abs_constraints <- canon[[2]]
    constraints <- c(constraints, abs_constraints)
  }

  # Now, we take care of the remaining convex and concave branches to create the
  # rational powers. We need a new variable, r, and the constraint sum(r) == t
  r <- Variable(shape(x))
  constraints <- c(constraints, list(sum(r) == t))

  # TODO: No need to run gm_constr to form the tree each time.
  # We only need to form the tree once.
  promoted_t <- Constant(matrix(1, nrow = shape(x)[1], ncol = shape(x)[2])) %*% t
  p <- Fraction(p)
  if(p < 0)
    constraints <- c(constraints, gm_constrs(promoted_t, list(x, r), c(-p/(1-p), 1/(1-p))))
  else if(p > 0 && p < 1)
    constraints <- c(constraints, gm_constrs(r, list(x, promoted_t), c(p, 1-p)))
  else if(p > 1)
    constraints <- c(constraints, gm_constrs(x, list(r, promoted_t), c(1/p, 1-1/p)))
  return(list(t, constraints))
}

power_canon <- function(expr, args) {
  x <- args[[1]]
  p <- expr@p
  w <- expr@w

  if(p == 1)
    return(list(x, list()))

  shape <- shape(expr)
  ones <- Constant(matrix(1, nrow = shape[1], ncol = shape[2]))
  if(p == 0)
    return(list(ones, list()))
  else {
    t <- Variable(shape)
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

quad_form_canon <- function(expr, args) {
  decomp <- decomp_quad(value(args[[1]]))
  scale <- decomp[[1]]
  M1 <- decomp[[2]]
  M2 <- decomp[[3]]

  if(size(M1) > 0)
    expr <- sum_squares(Constant(t(M1)) %*% args[[1]])
  else if(size(M2)> 0) {
    scale <- -scale
    expr <- sum_squares(Constant(t(M2)) %*% args[[1]])
  }
  canon <- quad_over_lin_canon(expr, expr@args)
  obj <- canon[[1]]
  constr <- canon[[2]]
  return(list(scale * obj, constr))
}

quad_over_lin_canon <- function(expr, args) {
  # quad_over_lin := sum_{ij} X^2_{ij} / y
  x <- args[[1]]
  y <- matrix(args[[1]], ncol = 1)
  # Pre-condition: shape = ()
  t <- Variable(1)
  # (y+t, y-t, 2*x) must lie in the second-order cone, where y+t is the scalar part
  # of the second-order cone constraint
  constraints <- list(SOC(t = y+t, X = hstack(list(y-t, 2*matrix(x, ncol = 1))), axis = 1))
  return(list(t, constraints))
}

sigma_max_canon <- function(expr, args) {
  A <- args[[1]]
  shape <- shape(A)
  n <- shape[1]
  m <- shape[2]
  X <- Variable(c(n+m, n+m), PSD = TRUE)

  shape <- shape(expr)
  t <- Variable(shape)
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
