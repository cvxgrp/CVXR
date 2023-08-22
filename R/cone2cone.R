APPROX_CONES <- list("RelEntrConeQuad" = list("SOC"),
                     "OpRelEntrConeQuad" = list("PSD"))

# Helper function for returning the weights and nodes for an n-point Gauss-Legendre quadrature on [0, 1].
gauss_legendre <- function(n) {
  beta <- 0.5/sqrt(rep(1, n-1) - (2*seq(1, n-1))^(-2))
  
  Tmat <- matrix(0, nrow = n, ncol = n)
  for(i in seq(1, n - 1)) {
    Tmat[i,i+1] <- beta[i]
    Tmat[i+1,i] <- beta[i]
  }
  
  Tmat_eig <- eigen(Tmat, only.values = FALSE)
  D <- Tmat_eig$values
  V <- Tmat_eig$vectors
  
  x <- D
  x <- sort(x)
  i <- order(x)
  w <- 2*(V[1,i])^2
  x <- (x + 1)/2
  w <- w/2
  return(list(w, x))
}

# For each i, enforce a constraint that (X[i,], y[i], z[i]) belongs to the rotated quadratic cone
# { (x,y,z) : || x ||^2 <= y z, 0 <= (y, z) }. This implementation doesn't enforce (x,y) >= 0!
# That should be imposed by the calling function.
rotated_quad_cone <- function(X, y, z) {
  m <- size(y)
  if(size(z) != m)
    stop("z must have size ", m)
  if(nrow(X) != m)
    stop("X must have ", m, " rows")
  if(length(dim(X)) < 2)
    X <- Reshape(X, c(m,1))
  
  #####################################
  # Comments from Dcp2Cone.quad_over_lin_canon:
  #   quad_over_lin := sum_{i} x^2_{i} / y
  #   t = Variable(1,) is the epigraph variable.
  #   Becomes a constraint
  #   SOC(t=y + t, X=[y - t, 2*x])
  ####################################
  soc_X_col0 <- Reshape(y - z, c(m,1))
  soc_X <- HStack(soc_X_col0, 2*X)
  soc_t <- y + z
  con - SOC(t = soc_t, X = soc_X, axis = 1)
  return(con)
}

#'
#' Use linear and SOC constraints to approximately enforce
#' con@x * log(con@x / con@y) <= con@z.
#'
#' We rely on an SOC characterization of 2-by-2 PSD matrices.
#' Namely, a matrix
#'      [ a, b ]
#'      [ b, c ]
#' is PSD if and only if (a, c) >= 0 and a*c >= b^2.
#' That system of constraints can be expressed as
#' a >= quad_over_lin(b, c).
#'
#' Note: constraint canonicalization in CVXR uses a return format
#' list(lead_con, con_list) where lead_con is a Constraint that might be
#' used in dual variable recovery and con_list consists of extra
#' Constraint objects as needed.
#'
RelEntrConeQuad_canon <- function(con, args) {
  k <- con@k
  m <- con@m
  x <- con@x
  y <- con@y
  n <- size(x)
  
  # Z has been declared so as to allow for proper vectorization.
  Z <- Variable(k+1, n)
  gauss <- gauss_legendre(m)
  w <- gauss[[1]]
  t <- gauss[[2]]
  Tvar <- Variable(m, n)
  lead_con <- Zero(w %*% Tvar + con@z/2^k)
  constrs <- list(Zero(Z[1] - y))
  
  for(i in seq(k)) {
    # The following matrix needs to be PSD.
    #     [Z[i]  , Z[i+1]]
    #     [Z[i+1], x     ]
    # The below recipe for imposing a 2x2 matrix as PSD follows from pg. 35, Ex. 2.6
    # of Boyd's convex optimization. Where the constraint simply becomes a rotated 
    # quadratic cone, see Dcp2Cone.quad_over_lin_canon for the very similar scalar case.
    epi <- Z[i,]
    stackedZ <- Z[i+1,]
    cons <- rotated_quad_cone(stackedZ, epi, x)
    constrs <- c(constrs, cons)
    constrs <- c(constrs, list(epi >= 0, x >= 0))
  }
  
  for(i in seq(m)) {
    off_diag <- -(t[i]^0.5)*Tvar[i,]
    # The following matrix needs to be PSD.
    #     [ Z[k] - x - T[i] , off_diag      ]
    #     [ off_diag        , x - t[i]*T[i] ]
    epi <- (Z[k,] - x - Tvar[i,])
    cons <- rotated_quad_cone(off_diag, epi, x - t[i]*Tvar[i,])
    constrs <- c(constrs, cons)
    constrs <- c(constrs, list(epi >= 0, x - t[i]*Tvar[i,] >= 0))
  }
  
  return(list(lead_con, constrs))
}

OpRelEntrConeQuad_canon <- function(con, args) {
  k <- con@k
  m <- conm
  X <- con@X
  Y <- con@Y
  
  if(!is_real(X))
    stop("X must be real")
  if(!is_real(Y))
    stop("Y must be real")
  if(!is_real(con@Z))
    stop("Z must be real")
  
  X_dim <- dim(X)
  Zs <- lapply(seq(k+1), function(i) { Variable(X_dim[1], X_dim[2], symmetric = TRUE) })
  Ts <- lapply(seq(m+1), function(i) { Variable(X_dim[1], X_dim[2], symmetric = TRUE) })
  constrs <- list(Zeros(Zs[[1]] - Y))
  
  if(!is_symmetric(X)) {
    ut <- upper_tri(X)
    lt <- upper_tri(t(X))
    constrs <- c(constrs, ut == lt)
  }
  if(!is_symmetric(Y)) {
    ut <- upper_tri(Y)
    lt <- upper_tri(t(Y))
    constrs <- c(constrs, ut == lt)
  }
  if(!is_symmetric(con@Z)) {
    ut <- upper_tri(con@Z)
    lt <- upper_tri(t(con@Z))
    constrs <- c(constrs, ut == lt)
  }
  
  gauss <- gauss_legendre(m)
  w <- gauss[[1]]
  t <- gauss[[2]]
  sum_list <- lapply(seq(m), function(i) { w[i] * Ts[i] })
  lead_con <- Zero(AddExpression(sum_list) + con@Z/2^k)
  
  for(i in seq(k)) {
    #    [Z[i],   Z[i+1]]
    #    [Z[i+1], x     ]
    constrs <- c(constrs, Bmat(list(list(Zs[i], Zs[i+1]), list(t(Zs[i+1]), X))) %>>% 0)
  }
  
  for(i in seq(m)) {
    off_diag <- -(t[i]^0.5) * Ts[i]
    # The following matrix needs to be PSD.
    #    [Z[k] - x - T[i], off_diag]
    #    [off_diag,        x - t[i]*T[i]]
    constrs <- c(constrs, Bmat(list(list(Zs[k] - X - Ts[i], off_diag), list(t(off_diag), X - t[i] * Ts[i]))) %>>% 0)
  }
  
  return(list(lead_con, constrs))
}

von_neumann_entr_QuadApprox <- function(expr, args) {
  m <- expr@quad_approx[[1]]
  k <- expr@quad_approx[[2]]
  canon <- von_neumann_entr_canon(expr, args)
  epi <- canon[[1]]
  initial_cons <- canon[[2]]
  
  cons <- list()
  for(con in initial_cons) {
    if(is(con, "ExpCone")) {   # Should only hit this once.
      qa_con <- as_quad_approx(con, m, k)
      canon <- RelEntrConeQuad_canon(qa_con, NULL)
      qa_con_canon_lead <- canon[[1]]
      qa_con_canon <- canon[[2]]
      cons <- c(cons, qa_con_canon_lead)
      cons <- c(cons, qa_con_canon)
    } else
      cons <- c(cons, con)
  }
  return(list(epi, cons))
}

von_neumann_entr_canon_dispatch <- function(expr, args) {
  if(!is.null(expr@quad_approx))
    von_neumann_entr_QuadApprox(expr, args)
  else
    von_neumann_entr_canon(expr, args)
}

#'
#' The QuadApprox class.
#'
#' This class represents a quadratic approximation.
#'
#' @rdname QuadApprox-class
.QuadApprox <- setClass("QuadApprox", contains = "Canonicalization")

QuadApprox <- function(problem = NULL) { .QuadApprox(problem = problem) }

QuadApprox.CANON_METHODS <- list("RelEntrConeQuad" = RelEntrConeQuad_canon, 
                                 "OpRelEntrConeQuad" = OpRelEntrConeQuad_canon)

setMethod("initialize", "QuadApprox", function(.Object, ...) {
  callNextMethod(.Object, ..., canon_methods = QuadApprox.CANON_METHODS)
})

# An "exotic" cone is defined as any cone that isn't supported by ParamConeProg. 
# If ParamConeProg is updated to support more cones, then it may be necessary to change this file.
EXOTIC_CONES <- list("PowConeND" = list("PowCone3D"))

# con : PowConeND. We can extract metadata from this, e.g., con@alpha and con@axis.
# args : tuple of length two with W = args[[1]], z = args[[2]].
pow_nd_canon <- function(con, args) {
  data <- get_data(con)
  alpha <- data[[1]]
  axis <- data[[2]]
  alpha <- value(alpha)
  W <- args[[1]]
  z <- args[[2]]
  
  if(axis == 1) {
    W <- t(W)
    alpha <- t(alpha)
  }
  
  if(ndim(W) == 1) {
    W <- Reshape(W, c(size(W), 1))
    alpha <- Reshape(alpha, c(size(W), 1))
  }
  
  W_dim <- dim(W)
  n <- W_dim[1]
  k <- W_dim[2]
  if(n == 2)
    can_canon <- PowCone3D(W[1,], W[2,], z, alpha[1,])
  else {
    Tvar <- Variable(n-2, k)
    
    x_3d <- list()
    y_3d <- list()
    z_3d <- list()
    alpha_3d <- list()
    for(j in seq_len(k)) {
      x_3d <- c(x_3d, W[seq_len(n-1),j])
      y_3d <- c(y_3d, Tvar[,j])
      y_3d <- c(y_3d, W[n,j])
      z_3d <- c(z_3d, z[j])
      z_3d <- c(z_3d, Tvar[,j])
      
      r_nums <- alpha[,j]
      r_dens <- rev(base::cumsum(rev(r_nums)))   # Equivalent to sapply(seq(n), function(i) { sum(alpha[seq(i,n), j]) }).
      r <- r_nums / r_dens
      alpha_3d <- c(alpha_3d, r[seq_len(n-1)])
    }

    # TODO: Ideally, we should construct x, y, z, alpha_p3d by applying suitable sparse matrices to W, z, Tvar,
    # rather than using the HStack atom. (HStack will probably result in longer compile times).    
    x_3d <- do.call(HStack, x_3d)
    y_3d <- do.call(HStack, y_3d)
    z_3d <- do.call(HStack, z_3d)
    alpha_p3d <- do.call(HStack, alpha_3d)

    can_con <- PowCone3D(x_3d, y_3d, z_3d, alpha_p3d)
  }
  
  # Return a single PowCone3D constraint defined over all auxiliary variables needed for the reduction to go through. There are no "auxiliary constraints" beyond this one.
  return(list(can_con, list()))
}

#'
#' The Exotic2Common class.
#'
#' This class represents a reduction of an exotic cone to a common cone constraint.
#'
#' @rdname Exotic2Common-class
.Exotic2Common <- setClass("Exotic2Common", contains = "Canonicalization")

Exotic2Common <- function(problem = NULL) { .Exotic2Common(problem = problem) }

Exotic2Common.CANON_METHODS <- list("PowConeND" = pow_nd_canon)

setMethod("initialize", "Exotic2Common", function(.Object, ...) {
  callNextMethod(.Object, ..., canon_methods = Exotic2Common.CANON_METHODS)
})
