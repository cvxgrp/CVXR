## CVXPY SOURCE: cvxpy/reductions/solvers/utilities.py
######################################################################
#                      Solver utility functions.
######################################################################

#'
#' Stacks the values of the given variables.
#'
#' @param variables a list of variables.
#' @param default value to use when variable value is NA
#' @param byrow logical value indicating whether to fill the matrix by rows. Defaults to FALSE.
stack_vals <- function(variables, default, byrow = FALSE) {
  value <- list()
  for(variable in variables) {
    val <- value(variable)
    if(is.na(val))   # Unknown values.
      mat <- matrix(rep(default, size(variable)), ncol = 1)
    else
      mat <- matrix(val, ncol = 1, byrow = byrow)
    value <- c(value, list(mat))
  }
  # return(do.call("c", value))
  return(do.call("rbind", value))
}

expcone_permutor <- function(n_cones, exp_cone_order) {
  if(n_cones == 0)
    return(c())
  order <- rep(exp_cone_order, n_cones)   # e.g., c(1,0,2, 1,0,2, 1, 0,2, ...)
  offsets <- 3*do.call("c", lapply(seq(n_cones), function(i) { rep(i-1, 3) }))   # e.g., c(0,0,0, 3,3,3, 6,6,6, ...)
  perm <- order + offsets
  return(perm)
}

#'
#' Gets a specified value of a dual variable.
#'
#' @param result_vec A vector containing the dual variable values.
#' @param offset An offset to get correct index of dual values.
#' @param constraint A list of the constraints in the problem.
#' @return A list of a dual variable value and its offset.
extract_dual_value <- function(result_vec, offset, constraint) {
  value <- result_vec[seq(offset + 1, length.out = size(constraint))]
  if(size(constraint) == 1)
    value <- as.numeric(value)
  offset <- offset + size(constraint)
  return(list(value, offset))
}

#'
#' Gets the values of the dual variables.
#'
#' @param result_vec A vector containing the dual variable values.
#' @param parse_func Function handle for the parser.
#' @param constraints A list of the constraints in the problem.
#' @return A map of constraint ID to dual variable value.
get_dual_values <- function(result_vec, parse_func, constraints) {
  dual_vars <- list()
  offset <- 0
  for(constr in constraints) {
    # TODO: Reshape based on dual variable size.
    parsed <- parse_func(result_vec, offset, constr)
    dual_vars[[as.character(id(constr))]] <- parsed[[1]]
    offset <- parsed[[2]]
  }
  return(dual_vars)
}

