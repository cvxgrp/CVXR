#'
#' The SymData class.
#'
#' This class represents the symbolic information for the conic form convex optimization problem.
#'
#' @slot objective A \code{list} representing the objective.
#' @slot constraints A \code{list} of canonicalized constraints.
#' @slot constr_map A \code{list} mapping constraint type to a list of constraints.
#' @slot dims A \code{list} of the dimensions of the cones.
#' @slot var_offsets A \code{list} mapping variable ID to horizontal offset.
#' @slot var_sizes A \code{list} mapping variable ID to variable dimensions.
#' @slot x_length The length of the x vector.
#' @slot presolve_status A \code{character} string indicating the status of the pre-solver. May be NA if pre-solve has failed.
#' @aliases SymData
#' @export
.SymData <- setClass("SymData", representation(objective = "list", constraints = "list", solver = "Solver", .constr_map = "list", .dims = "numeric", .var_offsets = "list", .var_sizes = "list", .x_length = "numeric", .presolve_status = "character"),
                     prototype(.constr_map = NULL, .dims = NA_real_, .var_offsets = NULL, .var_sizes = NULL, .x_length = NA_real_, .presolve_status = NA_character_),
                     validity = function(object) {
                       if(!is.null(object@.constr_map))
                         stop("[Validation: SymData] .constr_map is an internal variable that should not be set by user")
                       if(!is.na(object@.dims))
                         stop("[Validation: SymData] .dims is an internal variable that should not be set by user")
                       if(!is.null(object@.var_offsets))
                         stop("[Validation: SymData] .var_offsets is an internal variable that should not be set by user")
                       if(!is.null(object@.var_sizes))
                         stop("[Validation: SymData] .var_sizes is an internal variable that should not be set by user")
                       if(!is.na(object@.x_length))
                         stop("[Validation: SymData] .x_length is an internal variable that should not be set by user")
                       if(!is.na(object@.presolve_status))
                         stop("[Validation: SymData] .presolve_status is an internal variable that should not be set by user")
                      })

SymData <- function(objective, constraints, solver) {
  .SymData(objective = objective, constraints = constraints, solver = solver)
}

setMethod("initialize", "SymData", function(.Object, objective, constraints, solver, .constr_map = NULL, .dims = NA_real_, .var_offsets = NULL, .var_sizes = NULL, .x_length = NA_real_, .presolve_status = NA_character_) {
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
  constr_map[[EQ]]   <- constraints[sapply(constraints, function(c) { is(c, "LinEqConstr") })]
  constr_map[[LEQ]]  <- constraints[sapply(constraints, function(c) { is(c, "LinLeqConstr") })]
  constr_map[[SOC]]  <- constraints[sapply(constraints, function(c) { is(c, "SOC") })]
  constr_map[[SDP]]  <- constraints[sapply(constraints, function(c) { is(c, "SDP") })]
  constr_map[[EXP]]  <- constraints[sapply(constraints, function(c) { is(c, "ExpCone") })]
  constr_map[[BOOL]] <- constraints[sapply(constraints, function(c) { is(c, "BoolConstr") })]
  constr_map[[INT]]  <- constraints[sapply(constraints, function(c) { is(c, "IntConstr") })]
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
  return(NA)   # TODO: Return objective and constr_map as well?
}

SymData.format_for_solver <- function(constr_map, solver) {
  dims <- list()
  dims[[EQ_DIM]]   <- sum(sapply(constr_map[[EQ]],  function(c) { size(c)[1] * size(c)[2] }))
  dims[[LEQ_DIM]]  <- sum(sapply(constr_map[[LEQ]], function(c) { size(c)[1] * size(c)[2] }))
  dims[[SOC_DIM]]  <- list()
  dims[[SDP_DIM]]  <- list()
  dims[[EXP_DIM]]  <- 0
  dims[[BOOL_IDS]] <- list()
  dims[[INT_IDS]]  <- list()
  
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

  # TODO: Create var_sizes sorted by unique var_id
  # var_sizes <- sapply(vars_, function(var_names) { var_names[[2]] })
  # size_prods <- sapply(var_sizes, function(var_size) { var_size[1]*var_size[2] })
  # var_offsets <- cumsum(c(0, head(size_prods, n = -1)))
  # vert_offset <- sum(size_prods)
  # list(var_offsets = var_offsets, var_sizes = var_sizes, vert_offset = vert_offset)
}

.MatrixCache <- setClass("MatrixCache", representation(coo_tup = "list", const_vec = "numeric", constraints = "list", x_length = "numeric", .size = "numeric", .param_coo_tup = "list"),
                         prototype(.size = NA_real_, .param_coo_tup = NULL), validity = function(object) {
                           if(!is.na(object@.size))
                             stop("[Validation: MatrixCache] .size is an internal variable that should not be set by user")
                           if(!is.null(object@.param_coo_tup))
                             stop("[Validation: MatrixCache] .param_coo_tup is an internal variable that should not be set by user")
                         })

MatrixCache <- function(coo_tup, const_vec, constraints, x_length) {
  .MatrixCache(coo_tup = coo_tup, const_vec = const_vec, constraints = constraints, x_length = x_length)
}

setMethod("initialize", "MatrixCache", function(.Object, coo_tup, const_vec, constraints, x_length, .size = NA_real_, .param_coo_tup = NULL) {
  .Object@coo_tup <- coo_tup
  .Object@const_vec <- const_vec
  .Object@constraints <- constraints
  
  rows <- sum(sapply(constraints, function(c) { size(c)[1] * size(c)[2] }))
  cols <- x_length
  .Object@.size <- c(rows, cols)
  .Object@.param_coo_tup <- list(list(), list(), list())
  .Object
})

setMethod("reset_param_data", "MatrixCache", function(object) {
  object@.param_coo_tup <- list(list(), list(), list())
  object
})

.MatrixData <- setClass("MatrixData", representation(sym_data = "SymData", solver = "Solver", .obj_cache = "MatrixCache", .eq_cache = "MatrixCache", .ineq_cache = "MatrixCache"),
                        prototype(.obj_cache = NULL, .eq_cache = NULL, .ineq_cache = NULL), validity = function(object) {
                          if(!is.null(object@.obj_cache))
                            stop("[Validation: MatrixData] .obj_cache is an internal variable that should not be set by user")
                          if(!is.null(object@.eq_cache))
                            stop("[Validation: MatrixData] .eq_cache is an internal variable that should not be set by user")
                          if(!is.null(object@.ineq_cache))
                            stop("[Validation: MatrixData] .ineq_cache is an internal variable that should not be set by user")
                        })

MatrixData <- function(sym_data, solver) { .MatrixData(sym_data = sym_data, solver = solver) }

setMethod("initialize", "MatrixData", function(.Object, sym_data, solver, .obj_cache = NULL, .eq_cache = NULL, .ineq_cache = NULL) {
  .Object@sym_data <- sym_data
  
  # Cache everything possible.
  .Object@obj_cache <- .init_matrix_cache(.dummy_constr(.Object), .Object@sym_data@x_length)
  .Object@obj_cache <- .lin_matrix(.Object, .Object@obj_cache, caching = TRUE)
  
  # Separate constraints based on the solver being used.
  constr_types <- split_constr(solver, .Object@sym_data@constr_map)
  eq_constr <- constr_types[[1]]
  ineq_constr <- constr_types[[2]]
  nonlin_constr <- constr_types[[3]]
  
  # Equality constraints.
  .Object@.eq_cache <- .init_matrix_cache(eq_constr, .Object@sym_data@x_length)
  .Object@.eq_cache <- .lin_matrix(.Object, .Object@eq_cache, cachine = TRUE)
  
  # Inequality constraints.
  .Object@.ineq_cache <- .init_matrix_cache(ineq_constr, .Object@sym_data@x_length)
  .Object@.ineq_cache <- .lin_matrix(.Object, .Object@.ineq_cache, caching = TRUE)
  
  # TODO: Nonlinear constraints require returning an oracle (function), which R does not support. Need an alternative way.
  .Object
})

setMethod(".dummy_constr", "MatrixData", function(object) {
  list(create_eq(object@sym_data@objective))
})

setMethod("get_objective", "MatrixData", function(object) {
  mat <- .cache_to_matrix(object, object@.mat_cache)
  c <- mat[[1]]
  offset <- mat[[2]]
  c <- as.vector(t(c))
  offset <- as.numeric(offset[[1]])
  # Negate offset because was negated before.
  list(c, -offset)
})

setMethod("get_eq_constr", "MatrixData", function(object) { .cache_to_matrix(object, object@.eq_cache) })
setMethod("get_ineq_constr", "MatrixData", function(object) { .cache_to_matrix(object, object@.ineq_cache) })

setMethod(".init_matrix_cache", "MatrixData", function(object, constraints, x_length) {
  rows <- sum(sapply(constraints, function(c) { size(c)[1] * size(c)[2] }))
  COO <- list(list(), list(), list())
  const_vec <- rep(0, rows)
  MatrixCache(COO, const_vec, constraints, x_length)
})

setMethod(".lin_matrix", signature(object = "MatrixData", mat_cache = "MatrixCache", caching = "logical"), function(object, mat_cache, caching = FALSE) {
  active_constr <- list()
  constr_offsets <- list()
  vert_offset <- 0
  
  for(constr in mat_cache@constraints) {
    # Process the constraint if it has a parameter and not caching or it doesn't have a parameter and caching.
    has_param <- length(get_expr_params(constr@expr)) > 0
    if((has_param && !caching) || (!has_param && caching)) {
      # If parameterized, convert the parameters into constant nodes.
      if(has_param)
        constr <- copy_constr(constr, replace_params_with_consts)   # TODO: Don't think I need to copy constraint since R automatically creates copies
      active_constr <- c(active_constr, constr)
      constr_offsets <- c(constr_offsets, vert_offset)
    }
    vert_offset <- vert_offset + size(constr)[1]*size(constr)[2]
  }
  
  # Convert the constraints into a matrix and vector offset and add them to the matrix cache.
  if(length(active_constr) > 0) {
    mat <- get_problem_matrix(active_constr, object@sym_data@var_offsets, constr_offsets)  # TODO: Need to implement for canonInterface
    V <- mat[[1]]
    I <- mat[[2]]
    J <- mat[[3]]
    const_vec <- mat[[4]]
    
    mat_cache@const_vec <- mat_cache@const_vec + const_vec
    mat_cache@coo_tup[[1]] <- c(mat_cache@coo_tup[[1]], V)
    mat_cache@coo_tup[[2]] <- c(mat_cache@coo_tup[[2]], I)
    mat_cache@coo_tup[[3]] <- c(mat_cache@coo_tup[[3]], J)
  }
  mat_cache
})

setMethod(".cache_to_matrix", signature(object = "MatrixData", mat_cache = "MatrixCache"), function(object, mat_cache) {
  # Get parameter values.
  param_cache <- .init_matrix_cache(mat_cache@constraints, mat_cache@.size[1])
  param_cache <- .lin_matrix(object, param_cache)   # TODO: Need to re-assign this
  mat_size <- mat_cache@.size
  rows <- mat_size[1]
  cols <- mat_size[2]
  
  # Create the constraints matrix. Combine the cached data with the parameter data.
  mat_coo_tup <- mat_cache@coo_tup
  V <- mat_coo_tup[[1]]
  I <- mat_coo_tup[[2]]
  J <- mat_coo_tup[[3]]
  
  param_coo_tup <- param_cache@coo_tup
  Vp <- param_coo_tup[[1]]
  Ip <- param_coo_tup[[2]]
  Jp <- param_coo_tup[[3]]
  
  if(length(V) + length(Vp) > 0) {
    mat <- sparseMatrix(i = I + Ip, j = J + Jp, x = V + Vp, dims = c(rows, cols))
    mat <- as.matrix(mat)
  } else   # Empty matrix
    mat <- matrix(0, rows = rows, cols = cols)
  
  const_vec <- mat_cache@const_vec + param_cache@const_vec
  list(mat, -const_vec)
})

# TODO: Fill this out when nonlinear constraints are finished
# setMethod(".nonlin_matrix", "MatrixData", function(object, nonlin_constr) {
#  rows <- sum(sapply(nonlin_constr, function(c) { size(c)[1] * size(c)[2] }))
#  cols <- object@sym_data@x_length
#  var_offsets <- object@sym_data@var_offsets
#  
#  big_x <- rep(0, cols)
# })

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
.ProblemData <- setClass("ProblemData", representation(sym_data = "SymData", matrix_data = "MatrixData", prev_result = "list"),
                                        prototype(sym_data = NULL, matrix_data = NULL, prev_result = NULL))
ProblemData <- function(sym_data, matrix_data, prev_result) { 
  .ProblemData(sym_data = sym_data, matrix_data = matrix_data, prev_result = prev_result)
}
