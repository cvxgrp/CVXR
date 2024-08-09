# QPSolver requires objectives to be stuffed in the following way.
#'
#' Is the QP objective stuffed?
#'
#' @param objective A \linkS4class{Minimize} or \linkS4class{Maximize} object representing the optimization objective.
#' @return Is the objective a stuffed QP?
is_stuffed_qp_objective <- function(objective) {
  expr <- expr(objective)
  return(inherits(expr, "AddExpression") && length(expr@args) == 2 && inherits(expr@args[[1]], "QuadForm") && inherits(expr@args[[2]], "MulExpression") && is_affine(expr@args[[2]]))
}

## File cvxpy/reductions/solvers/qp_solvers/

#'
#' A QP solver interface.
#'
# Slots IS_MIP, REQUIRES_CONSTR, and SUPPORTED_CONSTRAINTS are for internal use only!
setClass("QpSolver",
         slots = list(
           IS_MIP = "character",
           REQUIRES_CONSTR = "logical",
           SUPPORTED_CONSTRAINTS = "character"),
         prototype = list(
           IS_MIP = "IS_MIP",
           REQUIRES_CONSTR = FALSE,
           SUPPORTED_CONSTRAINTS = c("ZeroConstraint", "NonPosConstraint")),
         contains = "ReductionSolver")

#' @param object A \linkS4class{QpSolver} object.
#' @describeIn QpSolver What classes of constraints does the solver support?
setMethod("supported_constraints", "QpSolver", function(solver) { object@SUPPORTED_CONSTRAINTS })

#' @describeIn QPSolver Can the solver solve problems that do not have constraints?
setMethod("requires_constr", "QpSolver", function(solver) { object@REQUIRES_CONSTR })

#' @param object A \linkS4class{QpSolver} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn QpSolver Is this a QP problem?
setMethod("accepts", signature(object = "QpSolver", problem = "Problem"), function(object, problem) {
  is(problem, "ParamQuadProg") &&
    (mip_capable(object) || !is_mixed_integer(problem)) &&
    length(convex_attributes(list(problem@x))) == 0 &&
    (length(problem@constraints) > 0 || !requires_constr(object)) &&
    ## TO FIX: the next statement should not use class(x) %in% ... REPLACE
    ##all(sapply(problem@constraints, function(c) { class(c) %in% supported_constraints(object) }))
    do.call(all, lapply(problem@constraints, inherits, what = object@SUPPORTED_CONSTRAINTS))
})

QpSolver.prepare_data_and_inv_data <- function(object, problem) {
  inv_data <- list()
  inv_data[[object@VAR_ID]] <- id(problem@x)

  constr_map <- group_constraints(problem@constraints)
  inv_data[[object@DIMS]] <- data[[object@DIMS]] <- ConeDims(constr_map)

  # Add information about integer variables.
  inv_data[[object@IS_MIP]] <- is_mixed_integer(problem)

  data <- list()
  data[[PARAM_PROB]] <- problem
  list(problem = problem, data = data, inv_data = inv_data)
}

#' @describeIn QpSolver Constructs a QP problem data stored in a list
setMethod("perform", signature(object = "QpSolver", problem = "Problem"), function(object, problem) {
  # Construct QP problem data stored in a dictionary.
  # The QP has the following form
  #
  #    minimize 1/2 x' P x + q' x
  #    subject to A x = b
  #               F x <= g

  tmp_dat <- QPSolver.prepare_data_and_inv_data(object, problem)
  problem <- tmp_dat$problem
  data <- tmp_dat$data
  inv_data <- tmp_dat$inv_data

  tmp_parm <- apply_parameters(problem)
  ## Since problem is ParamQuadProg, we can use named args
  ## See method "apply_params" for ParamQuadProg in qp2quad_form.R
  P <- tmp_parm$P ## tmp_parm[[1L]]
  q <- tmp_parm$q ## tmp_parm[[2L]]
  d <- tmp_parm$d ## tmp_parm[[3L]]
  AF <- tmp_parm$AF ## tmp_parm[[4L]]
  bg <- tmp_parm$b  ## tmp_parm[[5L]]
  inv_data[[OFFSET]] <- d

  # Get number of variables.
  n <- size(problem@x)
  len_eq <- data[[object@DIMS]]@zero
  len_leq <- data[[object@DIMS]]@nonpos

  if(len_eq > 0) {
    A <- AF[1:len_eq,]
    b <- -bg[1:len_eq]
  } else {
    A <- Matrix(0, nrow = 0, ncol = n, sparse = TRUE)
    b <- matrix(0, nrow = 0, ncol = 1)
  }

  if(len_leq > 0 && len_leq < nrow(AF)) {
    Fmat <- AF[(len_eq + 1):nrow(AF),]
    g <- -bg[(len_eq + 1):nrow(AF)]
  } else {
    Fmat <- Matrix(0, nrow = 0, ncol = n, sparse = TRUE)
    g <- -matrix(0, nrow = 0, ncol = 1)
  }

  # Create dictionary with problem data.
  data[[P_KEY]] <- Matrix(P, sparse = TRUE)
  data[[Q_KEY]] <- q
  data[[A_KEY]] <- Matrix(A, sparse = TRUE)
  data[[B_KEY]] <- b
  data[[F_KEY]] <- Matrix(Fmat, sparse = TRUE)
  data[[G_KEY]] <- g
  data[[BOOL_IDX]] <- sapply(problem@x@boolean_idx, `[[`, 1L)
  data[[INT_IDX]] <- sapply(problem@x@integer_idx, `[[`, 1L)
  data$n_var <- n
  data$n_eq <- nrow(A)
  data$n_ineq <- nrow(Fmat)

  return(list(data, inv_data))
})
