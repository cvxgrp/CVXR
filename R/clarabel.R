#' An interface for the CLARABEL solver
#'
#' @name CLARABEL-class
#' @aliases CLARABEL
#' @rdname CLARABEL-class
#' @export
setClass("CLARABEL", representation(exp_cone_order = "numeric"),   # Order of exponential cone arguments for solver. Internal only!
                prototype(exp_cone_order = c(0, 1, 2)), contains = "ConicSolver")

#' @rdname CLARABEL-class
#' @export
CLARABEL <- function() { new("CLARABEL") }

# Solver capabilities.
#' @describeIn CLARABEL Can the solver handle mixed-integer programs?
setMethod("mip_capable", "CLARABEL", function(solver) { FALSE })
setMethod("requires_constr", "CLARABEL", function(solver) { TRUE })
setMethod("supported_constraints", "CLARABEL", function(solver) { c(supported_constraints(ConicSolver()), "SOC", "ExpCone", "PSDConstraint") })
## Clarabel also supports power cone, but we have not implemented power cones yet in CVXR!

# Map of CLARABEL status to CVXR status.
#' @param solver,object,x A \linkS4class{CLARABEL} object.
#' @param status A status code returned by the solver.
#' @describeIn CLARABEL Converts status returned by CLARABEL solver to its respective CVXPY status.
setMethod("status_map", "CLARABEL", function(solver, status) {
  clarabel_status_map <- clarabel::solver_status_descriptions()
  if(status == 2L)
    return(OPTIMAL)
  else if(status == 3L)
    return(INFEASIBLE)
  else if(status == 4L)
    return(UNBOUNDED)
  else if(status == 5L)
    return (OPTIMAL_INACCURATE)
  else if(status == 6L)
    return(INFEASIBLE_INACCURATE)
  else if(status == 7L)
    return(UNBOUNDED_INACCURATE)
  else if(status == 8L || status == 9L)
    return(USER_LIMIT)
  else if(status == 10L || status == 11L)
    return(SOLVER_ERROR)
  else
    stop("CLARABEL status unrecognized: ", status)
})

#' @describeIn CLARABEL Returns the name of the solver
setMethod("name", "CLARABEL", function(x) { CLARABEL_NAME })

#' @describeIn CLARABEL Imports the solver
## Since CLARABEL is now optional, we check if it is available
setMethod("import_solver", "CLARABEL", function(solver) {
    requireNamespace("clarabel", quietly = TRUE)
})

#' @param problem A \linkS4class{Problem} object.
#' @param constr A \linkS4class{Constraint} to format.
#' @param exp_cone_order A list indicating how the exponential cone arguments are ordered.
#' @describeIn CLARABEL Return a linear operator to multiply by PSD constraint coefficients.
setMethod("reduction_format_constr", "CLARABEL", function(object, problem, constr, exp_cone_order) {
  # Extract coefficient and offset vector from constraint.
  # Special cases PSD constraints, as CLARABEL expects constraints to be
  # imposed on solely the lower triangular part of the variable matrix.
  # Moreover, it requires the off-diagonal coefficients to be scaled by
  # sqrt(2).
  if(is(constr, "PSDConstraint")) {
    expr <- expr(constr)
    triangularized_expr <- scaled_upper_tri(expr + t(expr))/2
    extractor <- CoeffExtractor(InverseData(problem))
    Ab <- affine(extractor, triangularized_expr)
    A_prime <- Ab[[1]]
    b_prime <- Ab[[2]]

    # CLARABEL requests constraints to be formatted as Ax + s = b,
    # where s is constrained to reside in some cone. Here, however,
    # we are formatting the constraint as A"x + b" = -Ax + b; h ence,
    # A = -A", b = b".
    return(list(-1*A_prime, b_prime))
  } else
    callNextMethod(object, problem, constr, exp_cone_order)
})

#' @describeIn CLARABEL Returns a new problem and data for inverting the new solution
setMethod("perform", signature(object = "CLARABEL", problem = "Problem"), function(object, problem) {
  # Returns a new problem and data for inverting the new solution.
  data <- list()
  inv_data <- list()
  inv_data[[object@var_id]] <- id(variables(problem)[[1]])

  # Parse the coefficient vector from the objective.
  offsets <- ConicSolver.get_coeff_offset(problem@objective@args[[1]])
  data[[C_KEY]] <- offsets[[1]]
  data[[OFFSET]] <- offsets[[2]]
  data[[C_KEY]] <- as.vector(data[[C_KEY]])
  inv_data[[OFFSET]] <- data[[OFFSET]][1]

  # Order and group nonlinear constraints.
  constr_map <- group_constraints(problem@constraints)
  data[[ConicSolver()@dims]] <- ConeDims(constr_map)
  inv_data[[ConicSolver()@dims]] <- data[[ConicSolver()@dims]]

  # CLARABEL requires constraints to be specified in the following order:
  # 1) Zero cone.
  # 2) Non-negative orthant.
  # 3) SOC.
  # 4) PSD.
  # 5) Exponential
  # NOT YET IMPLEMENTED 6) Power cone
  zero_constr <- constr_map$ZeroConstraint
  neq_constr <- c(constr_map$NonPosConstraint, constr_map$SOC, constr_map$PSDConstraint, constr_map$ExpCone)
  inv_data[[object@eq_constr]] <- zero_constr
  inv_data[[object@neq_constr]] <- neq_constr

  # Obtain A, b such that Ax + s = b, s \in cones.
  # Note that CLARABEL mandates that the cones MUST be ordered with
  # zero cones first, then non-negative orthant, then SOC, then
  # PSD, then exponential.
  offsets <- group_coeff_offset(object, problem, c(zero_constr, neq_constr), object@exp_cone_order)
  data[[A_KEY]] <- offsets[[1]]
  data[[B_KEY]] <- offsets[[2]]
  return(list(object, data, inv_data))
})

#'
#' Extracts the dual value for constraint starting at offset.
#'
#' Special cases PSD constraints, as per the CLARABEL specification.
#'
#' @param result_vec The vector to extract dual values from.
#' @param offset The starting point of the vector to extract from.
#' @param constraint A \linkS4class{Constraint} object.
#' @return The dual values for the corresponding PSD constraints
CLARABEL.extract_dual_value <- function(result_vec, offset, constraint) {
  if(is(constraint, "PSDConstraint")) {
    dim <- nrow(constraint)
    upper_tri_dim <- floor(dim*(dim+1)/2)
    new_offset <- offset + upper_tri_dim
    upper_tri <- result_vec[(offset + 1):new_offset]
    full <- triu_to_full(upper_tri, dim)
    return(list(full, new_offset))
  } else
    return(extract_dual_value(result_vec, offset, constraint))
}

#' @param solution The raw solution returned by the solver.
#' @param inverse_data A list containing data necessary for the inversion.
#' @describeIn CLARABEL Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "CLARABEL", solution = "list", inverse_data = "list"), function(object, solution, inverse_data) {
  # Returns the solution to the original problem given the inverse_data.
  status <- status_map(object, solution$status)

  attr <- list()
  attr[[SOLVE_TIME]] <- solution$solve_time
  attr[[SETUP_TIME]] <- solution$i
  attr[[NUM_ITERS]] <- solution$iterations

  if(status %in% SOLUTION_PRESENT) {
    primal_val <- solution$obj_val
    opt_val <- primal_val + inverse_data[[OFFSET]]
    primal_vars <- list()
    var_id <- inverse_data[[object@var_id]]
    primal_vars[[as.character(var_id)]] <- as.matrix(solution$x)

    num_zero <- inverse_data[[ConicSolver()@dims]]@zero
    eq_idx <- seq_len(num_zero)
    ineq_idx <- seq(num_zero + 1, length.out = length(solution$z) - num_zero)
    eq_dual_vars <- get_dual_values(solution$z[eq_idx], CLARABEL.extract_dual_value, inverse_data[[object@eq_constr]])
    ineq_dual_vars <- get_dual_values(solution$z[ineq_idx], CLARABEL.extract_dual_value, inverse_data[[object@neq_constr]])

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
#' @describeIn CLARABEL Solve a problem represented by data returned from apply.
setMethod("solve_via_data", "CLARABEL", function(object, data, warm_start, verbose, feastol, reltol, abstol,
                                            num_iter, solver_opts, solver_cache) {
  if (missing(solver_cache)) solver_cache <- new.env(parent=emptyenv())
  # TODO: Cast A to dense because scs R package rejects sparse matrices?
  ## Fix until scs::scs can handle sparse symmetric matrices
  A  <- data[[A_KEY]]
  ## Fix for Matrix version 1.3
  ## if (inherits(A, "dsCMatrix")) A <- as(A, "dgCMatrix")
  ## if (!inherits(A, "dgCMatrix")) A  <- as(as(A, "CsparseMatrix"), "dgCMatrix")
  ## Matrix 1.5 change!
  if (!inherits(A, "dgCMatrix")) A  <- as(as(A, "CsparseMatrix"), "generalMatrix")

  args <- list(A = A, b = data[[B_KEY]], c = data[[C_KEY]])
  if(warm_start && !is.null(solver_cache) && length(solver_cache) > 0 && name(object) %in% names(solver_cache)) {
    args$x <- solver_cache[[name(object)]]$x
    args$y <- solver_cache[[name(object)]]$y
    args$s <- solver_cache[[name(object)]]$s
  }
  cones <- CLARABEL.dims_to_solver_dict(data[[ConicSolver()@dims]])

  ## if(!all(c(is.null(feastol), is.null(reltol), is.null(abstol)))) {
  ##   warning("Ignoring inapplicable parameter feastol/reltol/abstol for CLARABEL.")
  ## }
  solver_defaults  <- SOLVER_DEFAULT_PARAM$CLARABEL

  if(is.null(num_iter)) {
    num_iter <- solver_defaults$max_iter
  }
  if (is.null(reltol)) {
    reltol <- solver_defaults$tol_gap_rel
  }
  if (is.null(abstol)) {
    abstol  <- solver_defaults$tol_gap_abs
  }
  if (is.null(feastol)) {
    feastol  <- solver_defaults$tol_feas
  }

  if (is.null(verbose)) {
    verbose  <- solver_defaults$verbose
  }
  
  control <- clarabel::clarabel_control(max_iter = num_iter, verbose = verbose, tol_gap_rel = reltol, tol_gap_abs = abstol, tol_feas = feastol)

  #Fill in parameter values
  control[names(solver_opts)] <- solver_opts

  # Returns the result of the call to the solver.
  results <- clarabel::clarabel(A = args$A, b = args$b, q = args$c, cones = cones, control = control)
  if(!is.null(solver_cache) && length(solver_cache) > 0)
    solver_cache[[name(object)]] <- results
  return(results)
})


## #' Return an vectorization of symmetric matrix using the upper triangular part,
## #' still in column order.
## #' @param S a symmetric matrix
## #' @return vector of values
## vec_u <- function(S) {
##   n <- nrow(S)
##   sqrt2 <- sqrt(2.0)
##   upper_tri <- upper.tri(S, diag = FALSE)
##   S[upper_tri] <- S[upper_tri] * sqrt2
##   S[upper.tri(S, diag = TRUE)]
## }

## #' Return the symmetric matrix from the [vec] vectorization
## #' @param v a vector
## #' @return a symmetric matrix
## mat_u <- function(v) {
##   n <- (sqrt(8 * length(v) + 1) - 1) / 2
##   sqrt2 <- sqrt(2.0)
##   S <- matrix(0, n, n)
##   upper_tri <- upper.tri(S, diag = TRUE)
##   S[upper_tri] <- v / sqrt2
##   S <- S + t(S)
##   diag(S) <- diag(S) / sqrt(2)
##   S
## }


#'
#' Utility methods for special handling of semidefinite constraints.
#'
#' @param matrix The matrix to get the lower triangular matrix for
#' @return The lower triangular part of the matrix, stacked in column-major order
updated_scaled_lower_tri <- function(matrix) {
  rows <- cols <- nrow(matrix)
  entries <- floor(rows * (cols + 1)/2)
  
  row_arr <- seq_len(entries)
  
  col_arr <- matrix(seq_len(rows * cols), nrow = rows, ncol = cols)
  col_arr <- col_arr[lower.tri(col_arr, diag = TRUE)]
  val_arr <- matrix(sqrt(2), nrow = rows, ncol = cols)
  diag(val_arr) <- 1
  val_arr <- val_arr[lower.tri(val_arr, diag = T)]
  coeff <- Constant(sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = c(entries, rows*cols)))
  vectorized_matrix <- reshape_expr(matrix, c(rows * cols, 1))
  coeff %*% vectorized_matrix
}

#'
#' Utility methods for special handling of semidefinite constraints.
#'
#' @param matrix The matrix to get the upper triangular matrix for
#' @return The upper triangular part of the matrix, stacked in column-major order
scaled_upper_tri <- function(matrix) {
  rows <- cols <- nrow(matrix)
  entries <- floor(rows * (cols + 1)/2)
  
  row_arr <- seq_len(entries)
  
  col_arr <- matrix(seq_len(rows * cols), nrow = rows, ncol = cols)
  col_arr <- col_arr[upper.tri(col_arr, diag = TRUE)]
  val_arr <- matrix(sqrt(2), nrow = rows, ncol = cols)
  diag(val_arr) <- 1
  val_arr <- val_arr[upper.tri(val_arr, diag = T)]
  coeff <- Constant(sparseMatrix(i = row_arr, j = col_arr, x = val_arr, dims = c(entries, rows*cols)))
  vectorized_matrix <- reshape_expr(matrix, c(rows*cols, 1))
  coeff %*% vectorized_matrix
}


#'
#' Expands upper triangular to full matrix.
#'
#' @param upper_tri A matrix representing the uppertriangular part of the matrix,
#' stacked in column-major order
#' @param n The number of rows (columns) in the full square matrix.
#' @return A matrix that is the scaled expansion of the upper triangular matrix.
triu_to_full <- function(upper_tri, n) {
  # Expands n*floor((n+1)/2) upper triangular to full matrix.
  # Scales off-diagonal by 1/sqrt(2), as per the SCS specification.
  sqrt2 <- sqrt(2.0)
  result <- matrix(0, nrow = n, ncol = n)
  result[upper.tri(result, diag = TRUE)] <- upper_tri
  result <- result + t(result)
  diag(result) <- diag(result) / 2

  result_upper_tri <- upper.tri(result, diag = FALSE)
  result_lower_tri <- lower.tri(result, diag = FALSE)

  result[result_upper_tri] <- result[result_upper_tri] / sqrt2
  result[result_lower_tri] <- result[result_lower_tri] / sqrt2
  matrix(result, nrow = n * n, byrow = FALSE)
}

#' Utility method for formatting a ConeDims instance into a dictionary
#' that can be supplied to Clarabel
#' @param cone_dims A \linkS4class{ConeDims} instance.
#' @return The dimensions of the cones.
CLARABEL.dims_to_solver_dict <- function(cone_dims) {
  cones <- list(z = as.integer(cone_dims@zero),
                l = as.integer(cone_dims@nonpos),
                q = sapply(cone_dims@soc, as.integer),
                ep = as.integer(cone_dims@exp),
                s = sapply(cone_dims@psd, as.integer))
  return(cones)
}
