#' Expands n*floor((n+1)/2) lower triangular to full matrix and scales
#' off-diagonal by 1/sqrt(2), as per the SCS specification.
#' @param lower_tri A matrix representing the lower triangular part of
#'   the matrix, stacked in column-major order
#' @param n The number of rows (columns) in the full square matrix.
#' @return A matrix that is the scaled expansion of the lower
#'   triangular matrix.
tri_to_full <- function(lower_tri, n) {
  full <- matrix(0, nrow = n, ncol = n)
  full[upper.tri(full, diag = TRUE)] <- lower_tri
  full[lower.tri(full, diag = TRUE)] <- lower_tri

  unscaled_diag <- diag(full)
  full <- full/sqrt(2)
  diag(full) <- unscaled_diag

  matrix(full, nrow = n*n, byrow = FALSE)
}

#' Return `V` so that _vec[indices] belongs to the SCS-standard PSD
#' cone_ can be written in natural cvxpy syntax as _V >> 0_.
#' @param vec list of affine expressions
#' @param indices contains nonnegative integers, which can index into
#'   `vec`
#' @details This function is similar to `tri_to_full`, which is also
#'   found in this file. The difference is that this function works
#'   without indexed assignment `mat[i,j] = expr`. Such indexed
#'   assignment cannot be used, because this function builds a CVXR
#'   expression, rather than a numeric matrix.
SCS.psdvec_to_psdmat <- function(vec, indices) {
  n <- floor(sqrt(2 * length(indices)))
  zero_mat <- matrix(0, nrow = n, ncol = n)
  indseq <- seq_len(n)
  i_seq <- unlist(lapply(indseq, seq.int, to = dim))
  j_seq <- unlist(lapply(indseq, function(i) rep(i, times = dim - i + 1)))

  mats <- lapply(seq_along(indices),
                 function(i) {
                   mat <- zero_mat
                   idx <- indices[i]
                   r <- i_seq[i]; c <- j_seq[i];
                   if (r == c) {
                     mat[r, r] <- 1
                   } else {
                     mat[r, c] <- mat[c, r] <- 1 / sqrt(2)
                   }
                   vec[idx] * mat
                 })
  do.call(sum, mats)
}


#' An interface for the SCS solver
#'
#' @name SCS-class
#' @aliases SCS
#' @rdname SCS-class
#' @export
SCS <- setClass("SCS", contains = "ConicSolver",
                prototype = list(
                    EXP_CONE_ORDER = c(0L, 1L, 2L),
                    SUPPORTED_CONSTRAINTS = c(ConicSolver()@SUPPORTED_CONSTRAINTS, "SOC", "ExpCone", "PSDConstraint"),
                    REQUIRES_CONSTR = TRUE)
                )
## Not needed usually
## #' @rdname SCS-class
## #' @export
## SCS <- function() { new("SCS") }

# Solver capabilities.
#' @describeIn SCS Can the solver handle mixed-integer programs?
setMethod("mip_capable", "SCS", function(solver) { FALSE })
##setMethod("requires_constr", "SCS", function(solver) { TRUE })
##setMethod("supported_constraints", "SCS", function(solver) { c(supported_constraints(ConicSolver()), "SOC", "ExpCone", "PSDConstraint") })

# Map of SCS status to CVXR status.
#' @param solver,object,x A \linkS4class{SCS} object.
#' @param status A status code returned by the solver.
#' @describeIn SCS Converts status returned by SCS solver to its respective CVXPY status.
setMethod("status_map", "SCS", function(solver, status) {
  if(status == 1)
    return(OPTIMAL)
  else if(status == 2)
    return(OPTIMAL_INACCURATE)
  else if(status == -1)
    return(UNBOUNDED)
  else if(status == -6)
    return(UNBOUNDED_INACCURATE)
  else if(status == -2)
    return(INFEASIBLE)
  else if(status == -7)
    return(INFEASIBLE_INACCURATE)
  else if(status %in% c(-3, -4, -5))
    return(SOLVER_ERROR)
  else
    stop("SCS status unrecognized: ", status)
})

#' @describeIn SCS Returns the name of the solver
setMethod("name", "SCS", function(x) { SCS_NAME })

#' @describeIn SCS Imports the solver
##setMethod("import_solver", "SCS", function(solver) { requireNamespace("scs", quietly = TRUE) })
setMethod("import_solver", "SCS", function(solver) { TRUE }) ## we require scs

#' @describeIn SCS supports quadratic objective?
setMethod("supports_quad_obj", "SCS", function(solver) { TRUE }) ## we require scs > 3.0

#' Construct a linear operator to multiply by PSD constraint coefficients.
#' @details
#' Special cases PSD constraints, as SCS expects constraints to be
#' imposed on solely the lower triangular part of the variable matrix.
#' Moreover, it requires the off-diagonal coefficients to be scaled by
#' sqrt(2), and applies to the symmetric part of the constrained expression.
#' @param constr a list of constraints
#' @return a linear operator to multiply by PSD constraint coefficients.
#' @importFrom Matrix sparseMatrix
setMethod("psd_format_mat", "SCS",  function(solver, constr) {
  rows <- cols <- nrow(constr)
  ind_seq <- seq_len(rows)
  entries <- rows * (cols + 1) %/% 2
  val_array <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                    dims = c(rows, cols))
  val_tri <- lower.tri(val_array, diag = TRUE)
  col_arr <- which(val_tri)
  val_array[val_tri] <- sqrt(2.0)
  diag(val_array) <- 1.0
  n <- rows * cols
  scaled_lower_tri <- Matrix::sparseMatrix(i = seq_len(entries), j = col_arr,
                                           x = val_array@x, dims = c(entries, n))
  idx <- seq_len(n)
  val_symm <- rep(0.5 , 2 * n)
  tmp <- ind_seq - 1L
  res <- rows * tmp
  row_symm <- c(idx, idx)
  col_symm <- c(idx, unlist(lapply(ind_seq, `+`, res)))
  symm_matrix <- Matrix::sparseMatrix(i = row_symm, j = col_symm, x = val_symm,
                                      dims = c(n, n))
  scaled_lower_tri %*% symm_matrix
})

## THIS METHOD now dispatches to the super class
## #' @describeIn SCS Returns a new problem and data for inverting the new solution
## setMethod("perform", signature(object = "SCS", problem = "Problem"), function(object, problem) {
##   # Returns a new problem and data for inverting the new solution.
##   data <- list()
##   inv_data <- list()
##   inv_data[[object@var_id]] <- id(variables(problem)[[1]])

##   # Parse the coefficient vector from the objective.
##   offsets <- ConicSolver.get_coeff_offset(problem@objective@args[[1]])
##   data[[C_KEY]] <- offsets[[1]]
##   data[[OFFSET]] <- offsets[[2]]
##   data[[C_KEY]] <- as.vector(data[[C_KEY]])
##   inv_data[[OFFSET]] <- data[[OFFSET]][1]

##   # Order and group nonlinear constraints.
##   constr_map <- group_constraints(problem@constraints)
##   data[[ConicSolver()@dims]] <- ConeDims(constr_map)
##   inv_data[[ConicSolver()@dims]] <- data[[ConicSolver()@dims]]

##   # SCS requires constraints to be specified in the following order:
##   # 1) Zero cone.
##   # 2) Non-negative orthant.
##   # 3) SOC.
##   # 4) PSD.
##   # 5) Exponential.
##   zero_constr <- constr_map$Zero
##   neq_constr <- c(constr_map$NonPos, constr_map$SOC, constr_map$PSD, constr_map$ExpCone)
##   inv_data[[object@eq_constr]] <- zero_constr
##   inv_data[[object@neq_constr]] <- neq_constr

##   # Obtain A, b such that Ax + s = b, s \in cones.
##   # Note that SCS mandates that the cones MUST be ordered with
##   # zero cones first, then non-negative orthant, then SOC, then
##   # PSD, then exponential.
##   offsets <- group_coeff_offset(object, problem, c(zero_constr, neq_constr), object@exp_cone_order)
##   data[[A_KEY]] <- offsets[[1]]
##   data[[B_KEY]] <- offsets[[2]]
##   return(list(object, data, inv_data))
## })

#'
#' Extracts the dual value for constraint starting at offset.
#'
#' Special cases PSD constraints, as per the SCS specification.
#'
#' @param result_vec The vector to extract dual values from.
#' @param offset The starting point of the vector to extract from.
#' @param constraint A \linkS4class{Constraint} object.
#' @return The dual values for the corresponding PSD constraints
SCS.extract_dual_value <- function(result_vec, offset, constraint) {
  if(is(constraint, "PSDConstraint")) {
    dim <- nrow(constraint)
    lower_tri_dim <- (dim * (dim + 1)) %/% 2
    new_offset <- offset + lower_tri_dim
    lower_tri <- result_vec[(offset + 1):new_offset]
    full <- tri_to_full(lower_tri, dim)
    return(list(full, new_offset))
  } else
    return(extract_dual_value(result_vec, offset, constraint))
}

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A list containing data necessary for the inversion.
#' @describeIn SCS Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "SCS", solution = "list", inverse_data = "list"), function(object, solution, inverse_data) {
  # Returns the solution to the original problem given the inverse_data.
  status <- status_map(object, solution$info$status)

  attr <- list()
  attr[[SOLVE_TIME]] <- solution$info$solveTime
  attr[[SETUP_TIME]] <- solution$info$setupTime
  attr[[NUM_ITERS]] <- solution$info$iter

  if(status %in% SOLUTION_PRESENT) {
    primal_val <- solution$info$pobj
    opt_val <- primal_val + inverse_data[[OFFSET]]
    primal_vars <- list()
    var_id <- inverse_data[[object@var_id]]
    primal_vars[[as.character(var_id)]] <- as.matrix(solution$x)

    num_zero <- inverse_data[[object@dims]]@zero
    eq_idx <- seq_len(num_zero)
    ineq_idx <- seq(num_zero + 1, length.out = length(solution$y) - num_zero)
    eq_dual_vars <- get_dual_values(solution$y[eq_idx], SCS.extract_dual_value, inverse_data[[object@eq_constr]])
    ineq_dual_vars <- get_dual_values(solution$y[ineq_idx], SCS.extract_dual_value, inverse_data[[object@neq_constr]])

    dual_vars <- list()
    dual_vars <- utils::modifyList(dual_vars, eq_dual_vars)
    dual_vars <- utils::modifyList(dual_vars, ineq_dual_vars)
    return(Solution(status, opt_val, primal_vars, dual_vars, attr))
  } else
    return(failure_solution(status))
})

#' @param data Data generated via an apply call.
#' @param warm_start A boolean of whether to warm start the solver.
#' @param verbose A boolean of whether to enable solver verbosity.
#' @param feastol The feasible tolerance on the primal and dual residual.
#' @param reltol The relative tolerance on the duality gap.
#' @param abstol The absolute tolerance on the duality gap.
#' @param num_iter The maximum number of iterations.
#' @param solver_opts A list of Solver specific options
#' @param solver_cache Cache for the solver.
#' @describeIn SCS Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "SCS", function(object, data, warm_start, verbose, feastol, reltol, abstol,
                                            num_iter, solver_opts, solver_cache = new.env(parent = emptyenv()) ) {
  A  <- data[[A_KEY]]
  if (!inherits(A, "dgCMatrix")) A  <- as(as(A, "CsparseMatrix"), "generalMatrix")

  args <- list(A = A, b = data[[B_KEY]], c = data[[C_KEY]])
  if(warm_start && !is.null(solver_cache) && length(solver_cache) > 0 && name(object) %in% names(solver_cache)) {
    obj_name <- name(object)
    args$x <- solver_cache[[obj_name]]$x
    args$y <- solver_cache[[obj_name]]$y
    args$s <- solver_cache[[obj_name]]$s
  }
  cones <- SCS.dims_to_solver_dict(data[[ConicSolver()@dims]])

  ## if(!all(c(is.null(feastol), is.null(reltol), is.null(abstol)))) {
  ##   warning("Ignoring inapplicable parameter feastol/reltol/abstol for SCS.")
  ## }
  solver_defaults  <- SOLVER_DEFAULT_PARAM$SCS

  if(is.null(num_iter)) {
    num_iter <- solver_defaults$max_iters
  }
  if (is.null(reltol)) {
    reltol <- solver_defaults$eps_rel
  }
  if (is.null(abstol)) {
    abstol  <- solver_defaults$eps_abs
  }
  if (is.null(feastol)) {
    feastol  <- solver_defaults$eps_infeas
  }

  control =  scs::scs_control(max_iters = num_iter, verbose = verbose, eps_rel = reltol, eps_abs = abstol, eps_infeas = feastol)

  #Fill in parameter values
  control[names(solver_opts)] <- solver_opts

  # Returns the result of the call to the solver.
  results <- scs::scs(A = args$A, b = args$b, obj = args$c, cone = cones, control = control)
  if(!is.null(solver_cache) && length(solver_cache) > 0)
    solver_cache[[name(object)]] <- results
  return(results)
})
