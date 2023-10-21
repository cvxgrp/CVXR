## Return a diagonal sparse matrix of given size
sp.eye <- function(n, repr = c("C", "R", "T"), negate = FALSE) {
  repr <- match.arg(repr)
  ind <- seq.int(0, length.out = n)
  vals <- if (negate) rep(-1.0, n) else rep(1.0, n)
  Matrix::sparseMatrix(i = ind, j = ind, x = vals, index1 = FALSE, dims = c(n, n), repr = repr)
}

#' A wrapper for linear operators
#' @param linear_op a sparse matrix or a function that takes one argument X
#' @param dim the dimensions
#' @return a function that will perform the linear operation
as_block_diag_linear_operator <- function(matrices) {
  ncols <- lapply(matrices, ncol)
  col_indices <- c(0, cumsum(ncols))
  index_seq <- seq_along(matrices)
  do.call(rbind,
          lapply(index_seq,
                 function(i) matrices[[i]] %*% X[, col_indices[i]:col_indices[i + 1]]))
}

# Utility method for formatting a ConeDims instance into a dictionary
# that can be supplied to solvers.
dims_to_solver_dict <- function(cone_dims)  {
  list(f =  cone_dims@zero,
       l = cone_dims@nonneg,
       q = cone_dims@soc,
       ep = cone_dims@exp,
       s = cone_dims@psd,
       p = cone_dims@p3d
  }
}

#' The ConicSolver class.
#'
#' Conic solver class with reduction semantics.
#'
#' @rdname ConicSolver-class
ConicSolver <- setClass("ConicSolver",
                        slots = list(
                          # The key that maps to ConeDims in the data returned by apply()
                          ## DIMS = "character",  Already in super class ReductionSolver
                          # Every conic solver must support Zero and NonNeg constraints.
                          SUPPORTED_CONSTRAINTS = "character",
                          # Some solvers cannot solve problems that do not have constraints.
                          # For such solvers, REQUIRES_CONSTR should be set to True.
                          REQUIRES_CONSTR = "logical",
                          # Does it support quadratic objective?
                          SUPPORTS_QUAD_OBJECTIVE = "logical",
                          # If a solver supports exponential cones, it must specify the corresponding order
                          # The cvxpy standard for the exponential cone is:
                          #     K_e = closure{(x,y,z) |  z >= y * exp(x/y), y>0}.
                          # Whenever a solver uses this convention, EXP_CONE_ORDER should be [0, 1, 2].
                          EXP_CONE_ORDER = "integer"),
                        prototype =
                          list(SUPPORTED_CONSTRAINTS = c("ZeroConstraint", "NonNegConstraint"),
                               REQUIRES_CONSTR = FALSE,
                               SUPPORTS_QUAD_OBJECTIVE = FALSE,
                               EXP_CONE_ORDER = integer(0)),
                        contains = "ReductionSolver")

# Every conic solver must support Zero and NonNeg constraints.
setMethod("supported_constraints", "ConicSolver", function(solver) { solver@SUPPORTED_CONSTRAINTS })

# Some solvers cannot solve problems that do not have constraints.
# For such solvers, requires_constr should return TRUE.
setMethod("requires_constr", "ConicSolver", function(solver) { solver@REQUIRES_CONSTR })

# Does this support quadratic objective? By default no.
setMethod("supports_quad_obj", "ConicSolver", function(solver) { solver@SUPPORTS_QUAD_OBJECTIVE })

#' @param object A \linkS4class{ConicSolver} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn ConicSolver Can the problem be solved with a conic solver?
setMethod("accepts", signature(object = "ConicSolver", problem = "Problem"), function(object, problem) {
  is(problem, "ParamConeProg") &&
    (object@MIP_CAPABLE || !is_mixed_integer(problem)) &&
    length(convex_attributes(variables(problem))) == 0 &&
    (length(problem@constraints) > 0 || !object@REQUIRES_CONSTR) &&
    all(as.logical(lapply(problem@constraints, inherits, what = object@SUPPORTED_CONSTRAINTS)))
})

#' Create a spacing matrix in compressed column format
#'
#' Construct a sparse matrix in compressed column format with a
#' specific block structure.  The matrix has a series of diagonal
#' blocks of ones separated by blocks of zeros. Note that wrong inputs
#' can cause this routine to error out.
#'
#' @param dim the dimensions of the matrix, a vector of num_row and num_column
#' @param n the number of columns in the matrix
#' @param offset number of leading rows of zeros
#' @param streak the number of ones along the diagonal for each block
#' @param spacing the number of rows filled with zeros between each block
#' @param num_blocks the number of blocks
#' @return a sparse matrix in compressed column format with the specified block structure
#' @examples
#' mat <- create_spacing_matrix(m = 20, n = 15, spacing = 1, streak = 3, num_blocks = 4, offset = 2)
#' print(mat)
#' @importFrom Matrix sparseMatrix
#' @noRd
ConicSolver.get_spacing_matrix <- function(dim, spacing, streak, num_blocks, offset) {
  # Calculate number of non-zero entries
  num_values <- num_blocks * streak
  values <- rep(1.0, num_values)
  streak_plus_spacing <- streak + spacing
  row_arr <- offset +
    as.integer(
      matrix(seq.int(from = 0, length.out = num_blocks * streak_plus_spacing),
             nrow = streak_plus_spacing, ncol = num_blocks)[seq_len(streak), ]
    )
  col_arr <- seq.int(0, length.out = num_values)
  Matrix::sparseMatrix(i = row_arr, j = col_arr, x = values, dims = dim, index1 = FALSE)
}

#' Return a matrix to multiply by PSD constraint coefficients, default is identity
#'
#' @param constr the constraints
#' @importFrom Matrix sparseMatrix
#' @noRd
ConicSolver.psd_format_mat <- function(constr) {
  n <- length(constr)
  ind <- seq.int(0, length.out = n)
  Matrix::sparseMatrix(i = ind, j = ind, x = rep(1, n), index1 = FALSE, dims = c(n, n))
}

#' Format constraints and return a `ParamConeProg`
#' @details
#' The `ParamConeProg` will have problem data tensors that will yield the
#' coefficient `A` and offset `b` for the constraint in the following formats:
#'   - Linear equations: `(A, b)` such that `A * x + b == 0`,
#'   - Linear inequalities: `(A, b)` such that `A * x + b >= 0`,
#'   - Second order cone: `(A, b)` such that `A * x + b` in SOC,
#'   - Exponential cone: `(A, b)` such that `A * x + b` in EXP,
#'   - Semidefinite cone: `(A, b)` such that `A * x + b` in PSD,
#' The `CVXR` standard for the exponential cone is:
#'     K_e = closure{(x,y,z) |  z >= y * exp(x/y), y>0}.
#' Whenever a solver uses this convention, exp_cone_order should be
#' [0, 1, 2].
#' The `CVXR` standard for the second order cone is:
#'     SOC(n) = { x : x[0] >= norm(x[1:n], 2)  }.
#' All currently supported solvers use this convention.
#' Args:
#' @param object the solver object
#' @param problem a `ParamConeProb` object that is problem that is the provenance of the constraint
#' @param exp_cone_order a list indicating how the exponential cone arguments are ordered
#' @return a `ParamConeProg` with structured A.
#' @describeIn ConicSolver Return a list representing a cone program whose problem data tensors
setMethod("format_constr", "ConicSolver", function(object, problem, exp_cone_order) {
  constr <- problem@constraints
  restruct_mat <- list()
  for (ct in constr) {
    n <- length(ct)
    total_height <- sum(unlist(lapply(ct@args, function(c) dim(c)[1L])))
    if (inherits(ct, "ZeroConstraint")) {
      restruct_mat.append(list(sp.eye(n = n, negate = TRUE))) ## Why csr format? For perf in python? Ignoring
    } else if (inherits(ct, "NonnegConstraint")) {
      restruct_mat.append(list(sp.eye(n = n)))  ## Why csr format? For perf in python? Ignoring
    } else if (inherits(ct, "SOC")) {
      stopifnot(constr@axis == 0) ## SOC must be lowered to axis == 0
      t_spacer <- ConicSolver.get_spacing_matrix(dim = c(total_height, size(ct@args[1L])),
                                                 spacing = constr.args[2L][1L],
                                                 streak = 1,
                                                 num_blocks = size(ct@args[1L]),
                                                 offset = 0
                                                 )
      X_spacer <- ConicSolver.get_spacing_matrix(dim = c(total_height, size(ct@args[1L])),
                                                 spacing = 1,
                                                 streak = constr.args[2L][1L],
                                                 num_blocks = size(ct@args[1L]),
                                                 offset = 1
                                                 )
      restruct_mat.append(list(cbind2(t_spacer, X_spacer)))
    } else if (inherits(ct, "ExpCone")) {
      l <- length(exp_cone_order)
      args <- constr@args
      arg_mats <- list()
      for (i in seq_along(args)) {
        arg <- args[[i]]
        space_mat <- ConicSolver.get_spacing_matrix( dim = c(total_height, size(arg)),
                                                    spacing = l - 1,
                                                    streak = 1,
                                                    num_blocks = size(arg),
                                                    offset = exp_cone_order[[i]]
                                                    )
        arg_mats <- arg_mats.append(space_mat)
      }
      restruct_mat.append(list(do.call(cbind, arg_mats)))
      } else if (inherits(ct, "PSD")) {
        restruct_mat.append(list(ConicSolver.psd_format_mat(constr)))
      } else {
        stop("Unsupported constraint type.")
      }
  }
  # Form new ParamConeProg
  ## Much simplified version in R but keeping comments
  if (length(restruct_mat) > 0) {
    # TODO(akshayka): profile to see whether using linear operators
    # or bmat is faster

    restruct_mat <- as_block_diag_linear_operator(restruct_mat)
    a_size<- prod(dim(problem@A))
    n <- ncol(restruct_mat)
    # this is equivalent to but _much_ faster than:
    #    restruct_mat_rep = sp.block_diag([restruct_mat]*(problem.x.size + 1))
    #    restruct_A = restruct_mat_rep * problem.A
    unspecified <- a_size %/% n ; remainder <- a_size %% n;
    reshaped_A <- problem@A; dim(reshaped_A) <- c(ncol(restruct_mat), unspecified)
    # Because of a bug in scipy versions <  1.20, `reshape`
    # can overflow if indices are int32s.
  } else {
    restructured_A <- problem@A
  }
  ## Return new formatted parametrized cone program
  new("ParamConeProg",
      c = problem@c,
      x = problem@x,
      A = restructured_A,
      variables = problem@variables,
      var_id_to_col = problem@var_id_to_col,
      constraints = problem@constraints,
      parameters = problem@parameters,
      param_id_to_col = problem@param_id_to_col,
      P = problem@P,
      formatted = TRUE)
})

#' @param solution A \linkS4class{Solution} object to invert.
#' @param inverse_data A \linkS4class{InverseData} object containing data necessary for the inversion.
#' @describeIn ConicSolver Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "ConicSolver", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  # Returns the solution to the original problem given the inverse_data.
  status <- solution$status

  if(status %in% SOLUTION_PRESENT) {
    opt_val <- solution$value
    primal_vars <- list()
    primal_vars[[inverse_data[[object@var_id]]]] <- solution$primal
    eq_dual <- get_dual_values(solution$eq_dual, extract_dual_value, inverse_data[object@eq_constr])
    leq_dual <- get_dual_values(solution$ineq_dual, extract_dual_value, inverse_data[object@neq_constr])
    eq_dual <- utils::modifyList(eq_dual, leq_dual)
    dual_vars <- eq_dual
  } else {
    primal_vars <- list()
    primal_vars[[inverse_data[[object@var_id]]]] <- NA_real_
    dual_vars <- NA
    if(status == INFEASIBLE)
      opt_val <- Inf
    else if(status == UNBOUNDED)
      opt_val <- -Inf
    else
      opt_val <- NA_real_
  }
  return(Solution(status, opt_val, primal_vars, dual_vars, list()))
})

ConicSolver.prepare_data_and_inv_data <- function(object, problem) {
  inv_data <- list()
  inv_data[[object@VAR_ID]] <- id(problem@x)

  # Format constraints
  #
  # By default cvxpy follows the SCS convention, which requires
  # constraints to be specified in the following order:
  # 1. zero cone
  # 2. non-negative orthant
  # 3. soc
  # 4. psd
  # 5. exponential
  # 6. three-dimensional power cones
  if (!problem@formatted) {
    problem <- reduction_format_constr(object, problem)  ## Exp cone order is a characteristic of the specific solver
  }
  data <- list()
  data[[PARAM_PROB]] <- problem
  data[[DIMS]] <- inv_data[[DIMS]] <- problem@cone_dims
  constr_map <- problem@constr_map
  inv_data[[EQ_CONSTR]] <- constr_map$Zero
  ## For why I do this below, see function group_constraints in reductions.R
  k <- match(c("NonNeg", "SOC", "PSD", "ExpCone", "PowCone3D"),
             names(constr_map))
  inv_data[[NEQ_CONSTR]] <- constr_map[ k[!is.na(k)] ]

  list(problem = problem, data = data, inv_data = inv_data)
}

#' @describeIn ConicSolver Return problem data and data for inverting the solution
setMethod("perform", signature(object = "ConicSolver", problem = "Problem"), function(object, problem) {
  # This is a reference implementation following SCS conventions
  # Implementations for other solvers may amend or override the implementation entirely
  tmp_dat <- prepare_data_and_inv_data(object, problem)
  problem <- tmp_dat$problem
  data <- tmp_dat$data
  inv_data <- tmp_dat$inv_data

  ## Since problem is ParamConeProg, we can use named args
  ## See method "apply_params" for ParamConeProg in dcp2cone.R
  if (!is.na(problem@P) {
    tmp_parm <- apply_parameters(problem, quad_obj = TRUE)
  } else {
    tmp_parm <- apply_parameters(problem)
  }
  data[[P_KEY]] <- tmp_parm$P
  data[[C_KEY]] <- tmp_param$c
  inv_data[[OFFSET]] <- tmp_parm$d
  data[[A_KEY]] <- -tmp_parm$A
  data[[B_KEY]] <- tmp_parm$b
  return(list(data = data, inv_data = inv_data))
}

