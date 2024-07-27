## CVXPY SOURCE: cvxpy/reductions/cone2cone/affine2direct.py
## Remember "apply" is "perform" for us

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
  Kd[[Cone2Cone.FREE]] <- Kp@zero       # length of block of unconstrained variables.
  Kd[[Cone2Cone.NONNEG]] <- Kp@nonneg   # length of block of nonneg variables.
  Kd[[Cone2Cone.SOC]] <- Kp@soc         # lengths of blocks of soc-constrained variables.
  Kd[[Cone2Cone.PSD]] <- Kp@psd         # "orders" of PSD variables.
  Kd[[Cone2Cone.DUAL_EXP]] <- Kp@exp     # number of length-3 blocks of dual exp cone variables.
  Kd[[Cone2Cone.DUAL_POW3D]] <- Kp@p3d  # scale parameters for dual 3d power cones.

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
    for(con in constr_map$PowCone3D) {
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

#'
#' CVXR represents mixed-integer cone programs as
#'
#'      (Aff)   min{ t(c) %*% x : A %*% x + b in K,
#'                                x[bools] in {0, 1}, x[ints] in Z } + d.
#'
#' Some solvers do not accept input in the form (Aff). A general pattern we find
#' across solver types is that the feasible set is represented by
#'
#'      (Dir)   min{ f %*% y : G %*% y <=_{K_aff} h, y in K_dir
#'                             y[bools] in {0, 1}, y[ints] in Z } + d,
#'
#'    where K_aff is built from a list convex cones which includes the zero cone (ZERO),
#'    and K_dir is built from a list of convex cones which includes the free cone (FREE).
#'
#'    This reduction handles mapping back and forth between problems stated in terms
#'    of (Aff) and (Dir), by way of adding slack variables.
#'
#'    Notes
#'    -----
#'    Support for semidefinite constraints has not yet been implemented in this
#'    reduction.
#'
#'    If the problem has no integer constraints, then the Dualize reduction should be
#'    used instead.
#'
#'    Because this reduction is only intended for mixed-integer problems, this reduction
#'    makes no attempt to recover dual variables when mapping between (Aff) and (Dir).
#'
#' ======================================================================================
#'
#' "prob" is a ParamConeProg which represents
#'
#'    (Aff)   min{ t(c) %*% x : A %*% x + b in K,
#'                           x[bools] in {0, 1}, x[ints] in Z } + d.
#'
#' We return data for an equivalent problem
#'
#'    (Dir)   min{ f %*% y : G %*% y <=_{K_aff} h, y in K_dir
#'                           y[bools] in {0, 1}, y[ints] in Z } + d,
#'
#' where
#'
#'    (1) K_aff is built from cone types specified in "affine" (a list of strings),
#'    (2) a primal solution for (Dir) can be mapped back to a primal solution
#'        for (Aff) by selecting the leading ``size(c)`` block of y's components.
#'
#' In the returned dict "data", data[[A_KEY]] = G, data[[B_KEY]] = h, data[[C_KEY]] = f,
#' data[['K_aff']] = K_aff, data[['K_dir']] = K_dir, data[[BOOL_IDX]] = bools,
#' and data[[INT_IDX]] = ints. The rows of G are ordered according to Cone2Cone.ZERO,
#' then (as applicable) Cone2Cone.NONNEG, Cone2Cone.SOC, and Cone2Cone.EXP. If  "c"
#' is the objective vector in (Aff), then ``y[seq_len(size(c) - 1)]`` should contain
#' the optimal solution to (Aff). The columns of G correspond first to variables in
#' cones Cone2Cone.FREE, then Cone2Cone.NONNEG, then Cone2Cone.SOC, then Cone2Cone.EXP.
#' The length of the free cone is equal to ``size(c)``.
#'
#' Assumptions
#' -----------
#' The function call ``cdAb = apply_parameters(prob)`` returns (A,b) with
#' rows formatted first for the zero cone, then for the nonnegative orthant, then
#' second order cones, then the exponential cone. Removing this assumption will
#' require adding additional data to ParamConeProg objects.
#'
Slacks.perform <- function(prob, affine) {
  tmp <- apply_parameters(prob)   # A %*% x + b in K
  c <- tmp[[1]]
  d <- tmp[[2]]
  A <- tmp[[3]]
  b <- tmp[[4]]
  A <- -A   # A %*% x <=_K b.
  cone_dims <- cone_dims(prob)
  if(!is.null(cone_dims@psd)) {
    # This will need to account for different conventions: does order-n PSD constraint
    # give rise to n^2 rows in A, or floor(n*(n-1)/2) rows?
    stop("Unimplemented")
  }

  for(val in affine) {
    if(!(val %in% c(Cone2Cone.ZERO, Cone2Cone.NONNEG, Cone2Cone.EXP, Cone2Cone.SOC, Cone2Cone.POW3D)))
      stop("Unimplemented")
  }
  if(!(Cone2Cone.ZERO %in% affine))
    affine <- c(affine, Cone2Cone.ZERO)

  cone_lens <- list()
  cone_lens[[Cone2Cone.ZERO]] <- cone_dims@zero
  cone_lens[[Cone2Cone.NONNEG]] <- cone_dims@nonneg
  cone_lens[[Cone2Cone.SOC]] <- sum(cone_dims@soc)
  cone_lens[[Cone2Cone.EXP]] <- 3*cone_dims@exp
  cone_lens[[Cone2Cone.POW3D]] <- 3*length(cone_dims@p3d)

  # If the rows of A are formatted in an order different from
  # zero -> nonneg -> soc -> exp -> pow, then the below block of code should
  # change. Right now there isn't enough data in (c, d, A, b, cone_dims,
  # constr_map) which allows us to figure out the ordering of these rows.
  row_offsets <- list()
  row_offsets[[Cone2Cone.ZERO]] <- 0
  row_offsets[[Cone2Cone.NONNEG]] <- cone_lens[[Cone2Cone.ZERO]]
  row_offsets[[Cone2Cone.SOC]] <- cone_lens[[Cone2Cone.ZERO]] + cone_lens[[Cone2Cone.NONNEG]]
  row_offsets[[Cone2Cone.EXP]] <- cone_lens[[Cone2Cone.ZERO]] + cone_lens[[Cone2Cone.NONNEG]] + cone_lens[[Cone2Cone.SOC]]
  row_offsets[[Cone2Cone.POW3D]] <- cone_lens[[Cone2Cone.ZERO]] + cone_lens[[Cone2Cone.NONNEG]] + cone_lens[[Cone2Cone.SOC]] + cone_lens[[Cone2Cone.EXP]]

  A_aff <- list()
  b_aff <- list()
  A_slk <- list()
  b_slk <- list()
  total_slack <- 0
  for(co_type in c(Cone2Cone.ZERO, Cone2Cone.NONNEG, Cone2Cone.SOC, Cone2Cone.EXP, Cone2Cone.POW3D)) {
    # The order of the list in the for loop means that the matrix "G" in "G %*% z <=_{K_aff} h"
    # will always have rows ordered by the zero cone, then the nonnegative orthant,
    # then second order cones, and finally exponential cones. Changing the order
    # of items in this list would change the order of row blocks in "G".
    #
    # If the order is changed, then this affects which columns of the final matrix
    # "G" correspond to which types of cones. For example, c(Cone2Cone.ZERO,
    # Cone2Cone.SOC, Cone2Cone.EXP, Cone2Cone.POW3D, Cone2Cone.NONNEG) and Cone2Cone.NONNEG
    # is not in "affine", then the columns of G with nonnegative variables occur
    # after all free variables, soc variables, exp variables, and pow3d variables.

    co_dim <- cone_lens[[co_type]]
    if(co_dim > 0) {
      r <- row_offsets[[co_type]]
      A_temp <- A[seq(r, r + co_dim - 1),]
      b_temp <- b[seq(r, r + co_dim - 1)]
      if(co_type %in% affine) {
        A_aff <- c(A_aff, list(A_temp))
        b_aff <- c(b_aff, list(b_temp))
      } else {
        total_slack <- total_slack + length(b_temp)
        A_slk <- c(A_slk, list(A_temp))
        b_slk <- c(b_slk, list(b_temp))
      }
    }
  }

  K_dir <- list()
  K_dir[[Cone2Cone.FREE]] <- size(prob@x)
  if(Cone2Cone.NONNEG %in% affine)
    K_dir[[Cone2Cone.NONNEG]] <- 0
  else
    K_dir[[Cone2Cone.NONNEG]] <- cone_dims@nonneg
  if(Cone2Cone.SOC %in% affine)
    K_dir[[Cone2Cone.SOC]] <- list()
  else
    K_dir[[Cone2Cone.SOC]] <- cone_dims@soc
  if(Cone2Cone.EXP %in% affine)
    K_dir[[Cone2Cone.EXP]] <- 0
  else
    K_dir[[Cone2Cone.EXP]] <- cone_dims@exp
  K_dir[[Cone2Cone.PSD]] <- list()   # Not currently supported in this reduction
  K_dir[[Cone2Cone.DUAL_EXP]] <- 0  # Not currently supported in CVXR
  if(Cone2Cone.POW3D %in% affine)
    K_dir[[Cone2Cone.POW3D]] <- list()
  else
    K_dir[[Cone2Cone.POW3D]] <- cone_dims@p3d
  K_dir[[Cone2Cone.DUAL_POW3D]] <- list()   # Not currently supported in CVXR

  K_aff <- list()
  if(Cone2Cone.NONNEG %in% affine)
    K_aff[[Cone2Cone.NONNEG]] <- cone_dims@nonneg
  else
    K_aff[[Cone2Cone.NONNEG]] <- 0
  if(ConeCone.SOC %in% affine)
    K_aff[[Cone2Cone.SOC]] <- cone_dims@soc
  else
    K_aff[[Cone2Cone.SOC]] <- list()
  if(Cone2Cone.EXP %in% affine)
    K_aff[[Cone2Cone.EXP]] <- cone_dims@exp
  else
    K_aff[[Cone2Cone.EXP]] <- 0
  K_aff[[Cone2Cone.PSD]] <- list()   # Currently not supported in this reduction
  K_aff[[Cone2Cone.ZERO]] <- cone_dims@zero + total_slack
  if(Cone2Cone.POW3D %in% affine)
    K_aff[[Cone2Cone.POW3D]] <- cone_dims@p3d
  else
    K_aff[[Cone2Cone.POW3D]] <- list()

  data <- list()
  if(length(A_slk) > 0) {
    # We need to introduce slack variables.
    A_slk <- Matrix(do.call("rbind", A_slk), sparse = TRUE)
    eye <- sparseMatrix(seq(total_slack), seq(total_slack), x = rep(1, total_slack))
    if(length(A_aff) > 0) {
      A_aff <- Matrix(do.call("rbind", A_aff), sparse = TRUE)
      G <- Matrix(rbind(cbind(A_slk, eye), cbind(A_aff, matrix(0, nrow = total_slack, ncol = total_slack))), sparse = TRUE)
      h <- Reduce("c", c(b_slk, b_aff))   # Concatenate lists, then turn to vector
    } else {
      G <- cbind(A_slk, eye)
      h <- Reduce("c", b_slk)
    }
    f <- c(c, rep(0, total_slack))
  } else if(length(A_aff) > 0) {
    # No slack variables were introduced.
    G <- Matrix(do.call("rbind", A_aff), sparse = TRUE)
    h <- Reduce("c", b_aff)
    f <- c
  } else
    stop("Must have at least one slack or affine variable")

  data[[A_KEY]] <- G
  data[[B_KEY]] <- h
  data[[C_KEY]] <- f
  data[[BOOL_IDX]] <- sapply(prob@x@boolean_idx, function(t) { as.integer(t[1]) })
  data[[INT_IDX]] <- sapply(prob@x@integer_idx, function(t) { as.integer(t[1]) })
  data$K_dir <- K_dir
  data$K_aff <- K_aff

  inv_data <- list()
  inv_data$x_id <- id(prob@x)
  inv_data$K_dir <- K_dir
  inv_data$K_aff <- K_aff
  inv_data[[OBJ_OFFSET]] <- d

  return(list(data, inv_data))
}

Slack.invert <- function(solution, inv_data) {
  if(solution@status %in% SOLUTION_PRESENT) {
    prim_vars <- solution@primal_vars
    x <- prim_vars[[Cone2Cone.FREE]]
    x_id <- id(x)
    for(i in seq(length(solution@primal_vars))) {
      if(id(solution@primal_vars[[i]]) == x_id) {
        solution@primal_vars[[i]] <- NULL
        break
      }
    }
    prim_vars[[inv_data$x_id]] <- x
  }
  solution@opt_val <- solution@opt_val + inv_data[[OBJ_OFFSET]]
  return(solution)
}
