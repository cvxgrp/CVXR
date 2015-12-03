.SymData <- setClass("SymData", representation(objective = "list", constraints = "list", solver = "Solver", .constr_map = "list", .dims = "numeric", .var_offsets = "list", .var_sizes = "list", .x_length = "numeric", .presolve_status = "character"),
                     prototype(.constr_map = NA, .dims = NA_real_, .var_offsets = NA, .var_sizes = NA, .x_length = NA_real_, .presolve_status = NA_character_),
                     validity = function(object) {
                       if(!is.na(object@.constr_map))
                         stop("[Validation: SymData] .constr_map is an internal variable that should not be set by user")
                       if(!is.na(object@.dims))
                         stop("[Validation: SymData] .dims is an internal variable that should not be set by user")
                       if(!is.na(object@.var_offsets))
                         stop("[Validation: SymData] .var_offsets is an internal variable that should not be set by user")
                       if(!is.na(object@.var_sizes))
                         stop("[Validation: SymData] .var_sizes is an internal variable that should not be set by user")
                       if(!is.na(object@.x_length))
                         stop("[Validation: SymData] .x_length is an internal variable that should not be set by user")
                       if(!is.na(object@.presolve_status))
                         stop("[Validation: SymData] .presolve_status is an internal variable that should not be set by user")
                      })

SymData <- function(objective, constraints, solver) {
  .SymData(objective = objective, constraints = constraints, solver = solver)
}

setMethod("initialize", "SymData", function(.Object, objective, constraints, solver, .constr_map = NA, .dims = NA_real_, .var_offsets = NA, .var_sizes = NA, .x_length = NA_real_, .presolve_status = NA_character_) {
  .Object@objective <- objective
  .Object@constraints <- constraints
  .Object@.constr_map <- SymData.filter_constraints(constraints)
  .Object@.presolve_status <- SymData.presolve(objective, .Objective@.constr_map)
  .Object@.dims <- SymData.format_for_solver(.Object@.constr_map, solver)
  
  all_ineq <- c(.Object@.constr_map[[EQ]], .Object@.constr_map[[LEQ]])
  # CVXOPT can have variables that only live in NonLinearConstraints.
  if(name(solver) == CVXOPT)
    nonlinear <- .Object@.constr_map[[EXP]]
  else
    nonlinear <- list()
  var_data <- SymData.get_var_offsets(objective, all_ineq, nonlinear)
  .Object@.var_offsets <- var_data[[1]]
  .Object@.var_sizes <- var_data[[2]]
  .Object@.x_length <- var_data[[3]]
  .Object
})

SymData.filter_constraints <- function(constraints) {
  constr_map <- list()
  constr_map[[EQ]] <- constraints[sapply(constraints, function(c) { is(c, "LinEqConstr") })]
  constr_map[[LEQ]] <- constraints[sapply(constraints, function(c) { is(c, "LinLeqConstr") })]
  constr_map[[SOC]] <- constraints[sapply(constraints, function(c) { is(c, "SOC") })]
  constr_map[[SDP]] <- constraints[sapply(constraints, function(c) { is(c, "SDP") })]
  constr_map[[EXP]] <- constraints[sapply(constraints, function(c) { is(c, "ExpCone") })]
  constr_map[[BOOL]] <- constraints[sapply(constraints, function(c) { is(c, "BoolConstr") })]
  constr_map[[INT]] <- constraints[sapply(constraints, function(c) { is(c, "IntConstr") })]
  constr_map
}

SymData.presolve <- function(objective, constr_map) {
  # Remove redundant constraints
  constr_map <- lapply(constr_map, function(constraints) { 
      constr_ids <- sapply(constraints, function(c) { c@constr_id })
       constraints[!duplicated(constr_ids)]
    })
  
  # If there are no constraints, the problem is unbounded if any of the coefficients are non-zero.
  # If all the coefficients are zero, then return the constant term and set all variables to zero.
  # TODO: Deal with the case when constr_maps has no values
  
  # Remove constraints with no variables or parameters.
  for(key in c(EQ, LEQ)) {
    new_constraints <- list()
    for(constr in constr_map[[key]]) {
      vars_ <- get_expr_vars(constr@expr)
      if(length(vars_) == 0 && is.na(get_expr_params(constr@expr))) {
        prob <- get_problem_matrix(list(constr))   # TODO: Call to canonInterface
        V <- prob[[1]]
        I <- prob[[2]]
        J <- prob[[3]]
        coeff <- prob[[4]]
        sign <- sign(coeff)
        
        # For equality constraint, coeff must be zero.
        # For inequality (i.e. <= 0) constraint, coeff must be negative.
        if((key == EQ && !is_zero(sign)) || (key == LEQ && !is_negative(sign)))
          return(INFEASIBLE)
      } else
        new_constraints <- c(new_constraints, constr)
    }
    constr_map[[key]] <- new_constraints
  }
  return(NA)
}

SymData.format_for_solver <- function(constr_map, solver) {
  dims <- list()
  dims[[EQ_DIM]] <- sum(sapply(constr_map[[EQ]], function(c) { size(c)[1] * size(c)[2] }))
  dims[[LEQ_DIM]] <- sum(sapply(constr_map[[LEQ]], function(c) { size(c)[1] * size(c)[2] }))
  dims[[SOC_DIM]] <- list()
  dims[[SDP_DIM]] <- list()
  dims[[EXP_DIM]] <- 0
  dims[[BOOL_IDS]] <- list()
  dims[[INT_IDS]] <- list()
  
  # Formats nonlinear constraints for the solver
  for(constr_type in names(dims)) {
    if(!(constr_type %in% c(EQ, LEQ))) {
      for(constr in constr_map[[constr_type]])
        constr_format(constr, constr_map[[EQ]], constr_map[[LEQ]], dims, solver)
    }
  }
  dims
}

SymData.get_var_offsets <- function(objective, constraints, nonlinear) {
  vars_ <- get_expr_vars(objective)
  vars_ <- c(vars_, lapply(constraints, function(constr) { get_expr_vars(constr@expr) }))
  
  # If CVXOPT is the solver, some of the variables are in NonLinearConstraints.
  for(constr in nonlinear)
    vars_ <- c(vars_, lapply(variables(constr), function(nonlin_var) { get_expr_vars(nonlin_var) }))

  # TODO: Create var_names, var_sizes, and vert_offset
}

#'
#' The ProblemData class.
#'
#' This class represents the symbolic and numeric data for a problem.
#'
#' @slot sym_data A \S4class{SymData} object representing the symbolic data for the problem.
#' @slot matrix_data A \S4class{MatrixData} object representing the numeric data for the problem.
#' @slot prev_result A \code{list} representing the result of the last solve
#' @aliases ProblemData
#' @export
ProblemData <- setClass("ProblemData", representation(sym_data = "SymData", matrix_data = "MatrixData", prev_result = "list"))
