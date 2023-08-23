Cone2Cone.FREE <- "fr"
Cone2Cone.ZERO <- "0"
Cone2Cone.NONNEG <- "+"
Cone2Cone.EXP <- "e"
Cone2Cone.DUAL_EXP <- "de"
Cone2Cone.SOC <- "q"
Cone2Cone.PSD <- "s"
Cone2Cone.POW3D <- "pp3"
Cone2Cone.DUAL_POW3D <- "dp3"

#'
#'    CVXR represents cone programs as
#'
#'        (P-Opt) min{ t(c) %*% x : A %*% x + b in K } + d,
#'
#'    where the corresponding dual is
#'
#'        (D-Opt) max{ -b %*% y : c = t(A) %*% y, y in K^* } + d.
#'
#'    For some solvers, it is much easier to specify a problem of the form (D-Opt) than it
#'    is to specify a problem of the form (P-Opt). The purpose of this reduction is to handle
#'    mapping between (P-Opt) and (D-Opt) so that a solver interface can pretend the original
#'    problem was stated in terms (D-Opt).
#'
#'    Usage
#'    -----
#'    Dualize applies to ParamConeProg problems. It accesses (P-Opt) data by calling
#'    ``c, d, A, b = apply_parameters(problem)``. It assumes the solver interface
#'    has already executed its ``format_constraints`` function on the ParamConeProg problem.
#'
#'    A solver interface is responsible for calling both Dualize.perform and Dualize.invert.
#'    The call to Dualize.perform should be one of the first things that happens, and the
#'    call to Dualize.invert should be one of the last things that happens.
#'
#'    The "data" dict returned by Dualize.perform is keyed by A_KEY, B_KEY, C_KEY, and 'K_dir',
#'    which respectively provide the dual constraint matrix (t(A)), the dual constraint
#'    right-hand-side (c), the dual objective vector (-b), and the dual cones (K^*).
#'    The solver interface should interpret this data is a new primal problem, just with a
#'    maximization objective. Given a numerical solution, the solver interface should first
#'    construct a CVXR Solution object where `y` is a primal variable, divided into
#'    several blocks according to the structure of elementary cones appearing in K^*. The only
#'    dual variable we use is that corresponding to the equality constraint `c = A^T y`.
#'    No attempt should be made to map unbounded / infeasible status codes for (D-Opt) back
#'    to unbounded / infeasible status codes for (P-Opt); all such mappings are handled in
#'    Dualize.invert. Refer to Dualize.invert for detailed documentation.
#'
#'    Assumptions
#'    -----------
#'    The problem has no integer or boolean constraints. This is necessary because strong
#'    duality does not hold for problems with discrete constraints.
#'
#'    Dualize.perform assumes "SOLVER.format_constraints()" has already been called. This
#'    assumption allows flexibility in how a solver interface chooses to vectorize a
#'    feasible set (e.g. how to order conic constraints, or how to vectorize the PSD cone).
#'
#'    Additional notes
#'    ----------------
#'
#'    Dualize.invert is written in a way which is agnostic to how a solver formats constraints,
#'    but it also imposes specific requirements on the input. Providing correct input to
#'    Dualize.invert requires consideration to the effect of ``SOLVER.format_constraints`` and
#'    the output of ``apply_parameters(problem)``.
#'    
Dualize.perform <- function(problem) {
  tmp <- apply_parameters(problem)
  c <- tmp[[1]]
  d <- tmp[[2]]
  A <- tmp[[3]]
  b <- tmp[[4]]
  
  Kp <- cone_dims(problem)    # zero, nonneg, exp, soc, psd
  
  Kd <- list()
  Kd[[FREE]] <- Kp@zero       # length of block of unconstrained variables.
  Kd[[NONNEG]] <- Kp@nonneg   # length of block of nonneg variables.
  Kd[[SOC]] <- Kp@soc         # lengths of blocks of soc-constrained variables.
  Kd[[PSD]] <- Kp@psd         # "orders" of PSD variables.
  Kd[DUAL_EXP]] <- Kp@exp     # number of length-3 blocks of dual exp cone variables.
  Kd[[DUAL_POW3D]] <- Kp@p3d  # scale parameters for dual 3d power cones.
  
  data <- list()
  data[[A_KEY]] <- t(A)
  data[[B_KEY]] <- c
  data[[C_KEY]] <- -b
  data$K_dir <- Kd
  data$dualized <- TRUE
  
  inv_data <- list()
  inv_data[[OBJ_OFFSET]] <- d
  inv_data$constr_map <- problem@constr_map
  inv_data$x_id <- id(problem@x)
  inv_data$K_dir <- Kd
  inv_data$dualized <- TRUE
  
  return(list(data, inv_data))
}

#'
#'  ``solution`` is a CVXR Solution object, formatted where
#'
#'    (D-Opt) max{ -b %*% y : c = t(A) %*% y, y in K^* } + d
#'
#'    is the primal problem from the solver's perspective. The purpose of this function
#'    is to map such a solution back to the format
#'
#'    (P-Opt) min{ t(c) %*% x : A %*% x + b in K } + d.
#'
#'    This function handles mapping of primal and dual variables, and solver status codes.
#'    The variable "x" in (P-Opt) is trivially populated from the dual variables to the
#'    constraint "c = t(A) %*% y" in (D-Opt). Status codes also map back in a simple way.
#'
#'    Details on required formatting of solution.primal_vars
#'    ------------------------------------------------------
#'
#'    We assume the dict solution@primal_vars is keyed by string-enums Cone2Cone.FREE ('fr'), 
#'    Cone2Cone.NONNEG ('+'), Cone2Cone.SOC ('s'), Cone2Cone.PSD ('p'), and Cone2Cone.DUAL_EXP ('de'). 
#'    The corresponding values are described below.
#'
#'    solution@primal_vars[[Cone2Cone.FREE]] should be a single vector. It corresponds to the 
#'    (possibly concatenated) components of "y" which are subject to no conic constraints. 
#'    We map these variables back to dual variables for equality constraints in (P-Opt).
#'
#'    solution@primal_vars[[Cone2Cone.NONNEG]] should also be a single vector, this time giving 
#'    the possibly concatenated components of "y" which must be >= 0. We map these variables
#'    back to dual variables for inequality constraints in (P-Opt).
#'
#'    solution@primal_vars[[Cone2Cone.SOC]] is a list of vectors specifying blocks of "y" which 
#'    belong to the second-order-cone under the CVXR standard ({ z : z[1] >= || z[1 + seq_len(length(z) - 1)] || }).
#'    We map these variables back to dual variables for SOC constraints in (P-Opt).
#'
#'    solution@primal_vars[[Cone2Cone.PSD]] is a list of symmetric positive semidefinite 
#'    matrices which result by lifting the vectorized PSD blocks of "y" back into matrix form.
#'    We assign these as dual variables to PSD constraints appearing in (P-Opt).
#'
#'    solution@primal_vars[[Cone2Cone.DUAL_EXP]] is a vector of concatenated length-3 slices 
#'    of y, where each constituent length-3 slice belongs to dual exponential cone as implied 
#'    by the CVXR standard of the primal exponential cone (see constraints.R:ExpCone).
#'    We map these back to dual variables for exponential cone constraints in (P-Opt).
#'
Dualize.invert <- function(solution, inv_data) {
  status <- solution@status
  prob_attr <- solution@attr
  primal_vars <- NULL
  dual_vars <- NULL
  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution@opt_val + inv_data[[OBJ_OFFSET]]
    primal_vars <- list()
    primal_vars[[inv_data$x_id]] <- solution@dual_vars[[EQ_DUAL]]
    dual_vars <- list()
    direct_prims <- solution@primal_vars
    constr_map <- inv_data$constr_map
    
    i <- 1
    for(con in constr_map$Zero) {
      dv <- direct_prims[[Cone2Cone.FREE]][seq(i, i + size(con) - 1)]
      if(size(dv) > 1)
        dual_vars[[id(con)]] <- dv
      else
        dual_vars[[id(con)]] <- as.numeric(dv)   # TODO: Is this same as dv.item()?
      i <- i + size(con)
    }
    
    i <- 1
    for(con in constr_map$Nonneg) {
      dv <- direct_prims[[Cone2Cone.NONNEG]][seq(i, i + size(con) - 1)]
      if(size(dv) > 1)
        dual_vars[[id(con)]] <- dv
      else
        dual_vars[[id(con)]] <- as.numeric(dv)
      i <- i + size(con)
    }
    
    i <- 1
    for(con in constr_map$SOC) {
      block_len <- nrow(con)
      dv <- Reduce(base::c, direct_prims[[Cone2Cone.SOC]][seq(i, i + size(con) - 1)])
      dual_vars[[id(con)]] <- dv
      i <- i + block_len
    }
    
    i <- 1
    for(i in seq(length(constr_map$PSD))) {
      con <- constr_map$PSD[[i]]
      dv <- direct_prims[[Cone2Cone.PSD]][i]
      dual_vars[[id(con)]] <- dv
    }
    
    i <- 1
    for(con in constr_map$ExpCone) {
      dv <- direct_prims[[Cone2Cone.DUAL_EXP]][seq(i, i + size(con) - 1)]
      dual_vars[[id(con)]] <- dv
      i <- i + size(con)
    }
    
    i <- 1   # TODO: Should we reset i = 1 here?
    for(con in constr_map$PowCone) {
      dv <- direct_prims[[Cone2Cone.DUAL_POW3D]][seq_i, i + size(con) - 1]
      dual_vars[[id(con)]] <- dv
      i <- i + size(con)
    }
  } else if(status == INFEASIBLE) {
    status <- UNBOUNDED
    opt_val <- -Inf
  } else if(status == INFEASIBLE_INACCURATE) {
    status <- UNBOUNDED_INACCURATE
    opt_val <- -Inf
  } else if(status == UNBOUNDED) {
    status <- INFEASIBLE
    opt_val <- Inf
  } else if(status == UNBOUNDED_INACCURATE) {
    status <- INFEASIBLE_INACCURATE
    opt_val <- Inf
  } else {
    status <- SOLVER_ERROR
    opt_val <- NA_real_
  }
  sol <- Solution(status, opt_val, primal_vars, dual_vars, prob_attr)
  return(sol)
}

# TODO: Convert affine2direct.py

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
