#
# The SymData class.
#
# This class represents the symbolic information for the conic form convex optimization problem.
#
# @slot objective A \code{list} representing the objective.
# @slot constraints A \code{list} of canonicalized constraints.
# @slot .constr_map (Internal) A \code{list} mapping constraint type to a list of constraints.
# @slot .dims (Internal) A \code{list} of the dimensions of the cones.
# @slot .var_offsets (Internal) A \code{numeric} vector mapping variable ID to horizontal offset.
# @slot .var_dims (Internal) A \code{list} mapping variable ID to variable dimensions.
# @slot .x_length (Internal) The length of the x vector.
# @slot .presolve_status (Internal) A \code{character} string indicating the status of the pre-solver. May be \code{NA} if pre-solve has failed.
# @rdname SymData-class
.SymData <- setClass("SymData", representation(objective = "list", constraints = "list", .constr_map = "list", .dims = "list", .var_offsets = "numeric", .var_dims = "list", .x_length = "numeric", .presolve_status = "character"),
                     prototype(.constr_map = list(), .dims = list(), .var_offsets = NA_integer_, .var_dims = list(), .x_length = NA_real_, .presolve_status = NA_character_))

#
# Symbolic Data Constructor
#
# Construct a \linkS4class{SymData} object.
#
# @param objective A \code{list} representing the objective.
# @param constraints A \code{list} of canonicalized constraints.
# @param solver A \linkS4class{Solver} for which to format the data.
# @return A \linkS4class{SymData} object.
# @docType methods
# @rdname SymData
SymData <- function(objective, constraints, solver) {
  constr_map <- SymData.filter_constraints(constraints)
  tmp <- SymData.presolve(objective, constr_map)
  constr_map <- tmp$constr_map
  presolve_status <- tmp$status

  tmp <- SymData.format_for_solver(constr_map, solver)
  constr_map <- tmp$constr_map
  dims <- tmp$dims

  all_ineq <- c(constr_map[[EQ_MAP]], constr_map[[LEQ_MAP]])
  # CVXOPT can have variables that only live in NonLinearConstraints.
  if(name(solver) == CVXOPT_NAME)
    nonlinear <- constr_map[[EXP_MAP]]
  else
    nonlinear <- list()
  var_data <- SymData.get_var_offsets(objective, all_ineq, nonlinear)
  .SymData(objective = objective, constraints = constraints, .constr_map = constr_map, .dims = dims, .var_offsets = var_data[[1]], .var_dims = var_data[[2]], .x_length = var_data[[3]], .presolve_status = presolve_status)
}

#
# Filter Constraints
#
# Separate the constraints by type.
#
# @param constraints A list of \linkS4class{Constraint} objects.
# @return A list of type to an ordered set of constraints. The types are linear equality (1), linear \eqn{leq} (2), SOC (3), SDP (4), exponential cone (5), boolean (6), and integer (7).
# @rdname SymData-filter_constraints
SymData.filter_constraints <- function(constraints) {
  ## constr_map <- list()
  ## constr_map[[EQ_MAP]]   <- if(length(constraints) == 0) list() else constraints[sapply(constraints, function(c) { is.list(c) && c$class == "LinEqConstr" })]
  ## constr_map[[LEQ_MAP]]  <- if(length(constraints) == 0) list() else constraints[sapply(constraints, function(c) { is.list(c) && c$class == "LinLeqConstr" })]
  ## constr_map[[SOC_MAP]]  <- if(length(constraints) == 0) list() else constraints[sapply(constraints, function(c) { is(c, "SOC") })]
  ## constr_map[[PSD_MAP]]  <- if(length(constraints) == 0) list() else constraints[sapply(constraints, function(c) { is(c, "SDP") })]
  ## constr_map[[EXP_MAP]]  <- if(length(constraints) == 0) list() else constraints[sapply(constraints, function(c) { is(c, "ExpCone") })]
  ## constr_map[[BOOL_MAP]] <- if(length(constraints) == 0) list() else constraints[sapply(constraints, function(c) { is(c, "BoolConstr") })]
  ## constr_map[[INT_MAP]]  <- if(length(constraints) == 0) list() else constraints[sapply(constraints, function(c) { is(c, "IntConstr") })]
  ## constr_map
    EQ_MAP_TAG <- 1L
    LEQ_MAP_TAG <- 2L
    SOC_MAP_TAG <- 3L
    SOC_EW_MAP_TAG <- 4L
    PSD_MAP_TAG <- 5L
    EXP_MAP_TAG <- 6L
    BOOL_MAP_TAG <- 7L
    INT_MAP_TAG <- 8L

    constr_tags <- sapply(constraints,
                          function(c) {
                              if (is.list(c) && c$class == "LinEqConstr" ) {
                                  EQ_MAP_TAG
                              } else if (is.list(c) && c$class == "LinLeqConstr") {
                                  LEQ_MAP_TAG
                              } else if (is(c, "SOC")) {
                                  SOC_MAP_TAG
                              } else if (is(c, "SDP")) {
                                  PSD_MAP_TAG
                              } else if (is(c, "ExpCone")) {
                                  EXP_MAP_TAG
                              } else if (is(c, "BoolConstr")) {
                                  BOOL_MAP_TAG
                              } else if (is(c, "IntConstr")) {
                                  INT_MAP_TAG
                              } else {
                                  stop("SymData.filter_constraints: Unimplemented Constraint")
                              }
                          })

    constr_map <- list(constraints[constr_tags == EQ_MAP_TAG],
                       constraints[constr_tags == LEQ_MAP_TAG],
                       constraints[constr_tags == SOC_MAP_TAG],
                       constraints[constr_tags == PSD_MAP_TAG],
                       constraints[constr_tags == EXP_MAP_TAG],
                       constraints[constr_tags == BOOL_MAP_TAG],
                       constraints[constr_tags == INT_MAP_TAG])
    names(constr_map) <- c(EQ_MAP, LEQ_MAP, SOC_MAP, PSD_MAP, EXP_MAP, BOOL_MAP, INT_MAP)
    constr_map
}

#
# Pre-Solve
#
# Eliminates unnecessary constraints and short circuits the solver if possible.
#
# @param objective A list representing the canonicalized objective.
# @param constr_map a list mapping constraint type to a list of constraints.
# @return A list containing the constraint map and feasibility status of the problem.
# @rdname SymData-presolve
SymData.presolve <- function(objective, constr_map) {
  # Remove redundant constraints
  constr_map <- lapply(constr_map, function(constraints) {
      constr_ids <- sapply(constraints, function(c) {
        if(is(c, "Constraint"))
          c@constr_id
        else if(is.list(c) && !is.null(c$constr_id))
          c$constr_id
        else
          stop("Invalid constraint class ", class(c))
        })
      constraints[!duplicated(constr_ids)]
  })

  # If there are no constraints, the problem is unbounded if any of the coefficients are non-zero.
  # If all the coefficients are zero, then return the constant term and set all variables to zero.

  # TODO: Deal with the case when constr_map has no values
  if(length(constr_map) == 0) {
    warning("Handling of empty constraint map is unfinished")
    print(objective)
  }

  # Remove constraints with no variables or parameters.
  for(key in c(EQ_MAP, LEQ_MAP)) {
    new_constraints <- list()
    for(constr in constr_map[[key]]) {
      vars_ <- get_expr_vars(constr$expr)
      if(length(vars_) == 0 && length(get_expr_params(constr$expr)) == 0) {
        prob <- get_problem_matrix(list(constr))
        V <- prob[[1]]
        I <- prob[[2]]
        J <- prob[[3]]

        coeff <- prob[[4]]
        sign <- intf_sign(coeff)
        is_pos <- sign[1]
        is_neg <- sign[2]

        # For equality constraint, coeff must be zero.
        # For inequality (i.e. <= 0) constraint, coeff must be negative.
        if((key == EQ_MAP && !(is_pos && is_neg)) || (key == LEQ_MAP && !is_neg))
          return(list(constr_map = constr_map, status = INFEASIBLE))
      } else
        new_constraints <- c(new_constraints, list(constr))
    }
    constr_map[[key]] <- new_constraints
  }
  list(constr_map = constr_map, status = NA_character_)
}

#
# Format for Solver
#
# Formats the problem for the solver.
#
# @param constr_map a list mapping constraint type to a list of constraints.
# @param solver A string indicating the solver being targeted.
# @return A list containing the constraint map and dimensions of the cones.
# @rdname SymData-format_for_solver
SymData.format_for_solver <- function(constr_map, solver) {
  dims <- list()
  if(length(constr_map[[EQ_MAP]]) > 0)
    dims[[EQ_DIM]]   <- sum(sapply(constr_map[[EQ_MAP]], function(c) { prod(dim(c)) }))
  else
    dims[[EQ_DIM]] <- 0
  if(length(constr_map[[LEQ_MAP]]) > 0)
    dims[[LEQ_DIM]]  <- sum(sapply(constr_map[[LEQ_MAP]], function(c) { prod(dim(c)) }))
  else
    dims[[LEQ_DIM]] <- 0
  dims[[SOC_DIM]]  <- c()
  dims[[PSD_DIM]]  <- c()
  dims[[EXP_DIM]]  <- 0
  dims[[BOOL_IDS]] <- c()
  dims[[INT_IDS]]  <- c()

  # Formats nonlinear constraints for the solver
  for(constr_type in names(constr_map)) {
    if(!(constr_type %in% c(EQ_MAP, LEQ_MAP))) {
      for(constr in constr_map[[constr_type]]) {
        tmp <- format_constr(constr, constr_map[[EQ_MAP]], constr_map[[LEQ_MAP]], dims, solver)
        constr_map[[EQ_MAP]] <- tmp$eq_constr
        constr_map[[LEQ_MAP]] <- tmp$leq_constr
        dims <- tmp$dims
      }
    }
  }
  list(constr_map = constr_map, dims = dims)
}

#
# Get Variable Offsets
#
# Maps each variable to a horizontal offset.
#
# @param objective A list representing the canonicalized objective.
# @param constraints A list of canonicalized constraints.
# @param nonlinear A list of nonlinear constraints.
# @return A list containing variable offsets, variable dimensions, and vertical offset.
# @rdname SymData-get_var_offsets
SymData.get_var_offsets <- function(objective, constraints, nonlinear) {
  vars_ <- get_expr_vars(objective)
  for(constr in constraints)
    vars_ <- c(vars_, get_expr_vars(constr$expr))

  # If CVXOPT is the solver, some of the variables are in NonLinearConstraints.
  for(constr in nonlinear) {
    for(nonlin_var in variables(constr))
      vars_ <- c(vars_, get_expr_vars(nonlin_var))
  }

  # Ensure variables are always in same order for same problem.
  var_names <- unique(vars_)
  var_ids <- sapply(var_names, function(id_and_dim) { id_and_dim[[1]] })
  names(var_names) <- var_ids
  var_names <- var_names[order(var_ids)]

  # Map variable IDs to offsets and dimensions.
  var_dims <- lapply(var_names, function(var) { var[[2]] })
  dim_prods <- sapply(var_dims, function(var_dim) { prod(var_dim) })
  if(length(dim_prods) == 0) {
    var_offsets <- integer(0)
    vert_offset <- 0
  } else {
    var_offsets <- base::cumsum(c(0, head(dim_prods, n = -1)))
    names(var_offsets) <- names(var_names)
    vert_offset <- sum(dim_prods)
  }
  list(var_offsets = var_offsets, var_dims = var_dims, vert_offset = vert_offset)
}

#
# The MatrixCache class.
#
# This class represents a cached version of the matrix and vector pair in an affine constraint.
#
# @slot coo_tup A \code{(V, I, J)} triplet for the COO (coordinate format) matrix. \code{I} are the row indices, \code{J} are the column indices, and \code{V} are the corresponding values.
# @slot .param_coo_tup (Internal) A \code{(V, I, J)} triplet for the parameterized COO matrix.
# @slot const_vec The vector offset.
# @slot constraints A list of constraints in the matrix.
# @slot .dim (Internal) The \code{c(rows, cols)} dimensions of the matrix.
# @rdname MatrixCache-class
.MatrixCache <- setClass("MatrixCache", representation(coo_tup = "list", .param_coo_tup = "list", const_vec = "numeric", constraints = "list", .dim = "numeric"),
                         prototype(.dim = NA_real_, .param_coo_tup = list(c(), c(), c())))

#
# Matrix Cache Constructor
#
# Construct a \linkS4class{MatrixCache} object.
#
# @param coo_tup A \code{(V, I, J)} triplet for the COO (coordinate format) matrix. \code{I} are the row indices, \code{J} are the column indices, and \code{V} are the corresponding values.
# @param const_vec The vector offset.
# @param constraints A list of constraints in the matrix.
# @param x_length The number of columns in the matrix.
# @return A \linkS4class{MatrixCache} object.
# @rdname MatrixCache
MatrixCache <- function(coo_tup, const_vec, constraints, x_length) {
  if(length(constraints) == 0)
    rows <- 0
  else
    rows <- sum(sapply(constraints, function(c) { prod(dim(c)) }))
  cols <- x_length
  .MatrixCache(coo_tup = coo_tup, const_vec = const_vec, constraints = constraints, .dim = c(rows, cols), .param_coo_tup = list(c(), c(), c()))
}

reset_param_data <- function(object) {
  object@.param_coo_tup <- list(c(), c(), c())
  object
}

#
# The MatrixData class.
#
# This class represents the matrices for the conic form convex optimization problem.
#
# @slot sym_data The \linkS4class{SymData} for the conic form problem.
# @slot solver A \linkS4class{Solver} for which to format the data.
# @slot nonlin A logical value indicating whether nonlinear constraints are needed.
# @slot .obj_cache (Internal) The \linkS4class{MatrixCache} for the objective to be used internally.
# @slot .eq_cache (Internal) The \linkS4class{MatrixCache} for the equality constraints to be used internally.
# @slot .ineq_cache (Internal) The \linkS4class{MatrixCache} for the inequality constraints to be used internally.
# @rdname MatrixData-class
.MatrixData <- setClass("MatrixData", representation(sym_data = "SymData", solver = "Solver", nonlin = "logical", .obj_cache = "MatrixCache", .eq_cache = "MatrixCache", .ineq_cache = "MatrixCache"),
                        prototype(.obj_cache = NULL, .eq_cache = NULL, .ineq_cache = NULL), validity = function(object) {
                          if(!is.null(object@.obj_cache))
                            stop("[Validation: MatrixData] .obj_cache is an internal variable that should not be set by user")
                          if(!is.null(object@.eq_cache))
                            stop("[Validation: MatrixData] .eq_cache is an internal variable that should not be set by user")
                          if(!is.null(object@.ineq_cache))
                            stop("[Validation: MatrixData] .ineq_cache is an internal variable that should not be set by user")
                          return(TRUE)
                        })

#
# Matrix Data Constructor
#
# Construct a \linkS4class{MatrixData} object.
#
# @param sym_data The \linkS4class{SymData} for the conic form problem.
# @param solver A \linkS4class{Solver} for which to format the data.
# @param nonlin A logical value indicating whether nonlinear constraints are needed.
# @return A \linkS4class{MatrixData} object.
# @docType methods
# @rdname MatrixData
MatrixData <- function(sym_data, solver, nonlin = FALSE) { .MatrixData(sym_data = sym_data, solver = solver, nonlin = nonlin) }

setMethod("initialize", "MatrixData", function(.Object, sym_data, solver, nonlin, .obj_cache = NULL, .eq_cache = NULL, .ineq_cache = NULL) {
    ##browser()
  .Object@sym_data <- sym_data

  # Cache everything possible.
  .Object@.obj_cache <- .init_matrix_cache(.Object, .dummy_constr(.Object), .Object@sym_data@.x_length)
  .Object@.obj_cache <- .lin_matrix(.Object, .Object@.obj_cache, caching = TRUE)

  # Separate constraints based on the solver being used.
  constr_types <- split_constr(solver, .Object@sym_data@.constr_map)
  eq_constr <- constr_types[[1]]
  ineq_constr <- constr_types[[2]]
  nonlin_constr <- constr_types[[3]]

  # Equality constraints.
  .Object@.eq_cache <- .init_matrix_cache(.Object, eq_constr, .Object@sym_data@.x_length)
  .Object@.eq_cache <- .lin_matrix(.Object, .Object@.eq_cache, caching = TRUE)

  # Inequality constraints.
  .Object@.ineq_cache <- .init_matrix_cache(.Object, ineq_constr, .Object@sym_data@.x_length)
  .Object@.ineq_cache <- .lin_matrix(.Object, .Object@.ineq_cache, caching = TRUE)

  # TODO: Nonlinear constraints require returning an oracle (function), which R does not support. Need an alternative way.
  .Object@nonlin <- nonlin
  if(nonlin) stop("Nonlinear constraints are unimplemented")
  # if(nonlin)
  #  .Object@F <- .nonlin_matrix(nonlin_constr)
  # else
  #  .Object@F <- NULL
  .Object
})

# Returns a dummy constraint for the objective.
.dummy_constr <- function(object) { list(create_eq(object@sym_data@objective)) }

# Returns the linear objective and scalar offset.
setMethod("get_objective", "MatrixData", function(object) {
  mat <- .cache_to_matrix(object, object@.obj_cache)
  c <- mat[[1]]
  offset <- mat[[2]]
  c <- as.numeric(t(c))
  offset <- as.numeric(offset[[1]])
  # Negate offset because was negated before.
  list(c, -offset)
})

# Returns the matrix and vector for the equality constraint.
setMethod("get_eq_constr", "MatrixData", function(object) { .cache_to_matrix(object, object@.eq_cache) })

# Returns the matrix and vector for the inequality constraint.
setMethod("get_ineq_constr", "MatrixData", function(object) { .cache_to_matrix(object, object@.ineq_cache) })

# Returns the oracle function for the nonlinear constraints.
setMethod("get_nonlin_constr", "MatrixData", function(object) { object@F} )

#
# Initialize Matrix Cache
#
# Initializes the data structures for the cached matrix.
#
# @param constraints A list of constraints in the matrix.
# @param x_length The number of columns in the matrix.
# @return A \linkS4class{MatrixCache} object.
# @rdname init_matrix_cache-int
.init_matrix_cache <- function(object, constraints, x_length) {
  if(length(constraints) == 0)
    rows <- 0
  else
    rows <- sum(sapply(constraints, function(c) { prod(dim(c)) }))
  COO <- list(c(), c(), c())
  const_vec <- rep(0, rows)
  MatrixCache(COO, const_vec, constraints, x_length)
}

#
# Constraints to Matrix/Vector
#
# Computes a matrix and vector representing a list of constraints. In the matrix, each constraint is given a block of rows. Each variable coefficient is inserted as a block with upper left corner at matrix[variable offset, constraint offset]. The constant term in the constraint is added to the vector.
#
# @param mat_cache A \linkS4class{MatrixCache} object representing the cached version of the matrix-vector pair.
# @param caching A logical value indicating whether the data should be cached.
# @return A \linkS4class{MatrixCache} object
# @rdname lin_matrix-int
.lin_matrix <- function(object, mat_cache, caching = FALSE) {
  active_constr <- list()
  constr_offsets <- c()
  vert_offset <- 0

  for(constr in mat_cache@constraints) {
    # Process the constraint if it has a parameter and not caching or it doesn't have a parameter and caching.
    has_param <- length(get_expr_params(constr$expr)) > 0
    if((has_param && !caching) || (!has_param && caching)) {
      # If parameterized, convert the parameters into constant nodes.
      if(has_param)
        constr <- copy_constr(constr, replace_params_with_consts)   # TODO: Don't think I need to copy constraint since R automatically creates copies
      active_constr <- c(active_constr, list(constr))
      constr_offsets <- c(constr_offsets, vert_offset)
    }
    vert_offset <- vert_offset + prod(dim(constr))
  }
  storage.mode(constr_offsets) <- "integer"

  # Convert the constraints into a matrix and vector offset and add them to the matrix cache.
  expr_list <- lapply(active_constr, function(con) { con$expr })
  if(length(active_constr) > 0) {
    mat <- get_problem_matrix(expr_list, object@sym_data@.var_offsets, constr_offsets)
    V <- mat[[1]]
    I <- mat[[2]]
    J <- mat[[3]]
    const_vec <- mat[[4]]

    # Convert the constant offset to the correct data type.
    conv_vec <- as.matrix(const_vec)

    mat_cache@const_vec[1:length(conv_vec)] <- mat_cache@const_vec[1:length(conv_vec)] + conv_vec
    mat_cache@coo_tup[[1]] <- c(mat_cache@coo_tup[[1]], V)
    mat_cache@coo_tup[[2]] <- c(mat_cache@coo_tup[[2]], I)
    mat_cache@coo_tup[[3]] <- c(mat_cache@coo_tup[[3]], J)
  }
  mat_cache
}

#
# Cache to Matrix
#
# Converts the cached representation of the constraints matrix.
#
# @param mat_cache A \linkS4class{MatrixCache} object representing the cached version of the matrix-vector pair.
# @return A list of \code{c(matrix, vector)}.
# @rdname cache_to_matrix-int
.cache_to_matrix <- function(object, mat_cache) {
  # Get parameter values.
  param_cache <- .init_matrix_cache(object, mat_cache@constraints, mat_cache@.dim[1])
  param_cache <- .lin_matrix(object, param_cache)
  mat_dim <- mat_cache@.dim
  rows <- mat_dim[1]
  cols <- mat_dim[2]

  # Create the constraints matrix. Combine the cached data with the parameter data.
  mat_coo_tup <- mat_cache@coo_tup
  V <- mat_coo_tup[[1]]
  I <- mat_coo_tup[[2]]
  J <- mat_coo_tup[[3]]

  param_coo_tup <- param_cache@coo_tup
  Vp <- param_coo_tup[[1]]
  Ip <- param_coo_tup[[2]]
  Jp <- param_coo_tup[[3]]

  if((length(V) + length(Vp)) > 0) {
    # TODO: R uses 1-indexing, but get_problem_matrix in canonInterface returns with 0-indexing
    mat <- sparseMatrix(i = c(I, Ip) + 1, j = c(J, Jp) + 1, x = c(V, Vp), dims = c(rows, cols))
    # mat <- as.matrix(mat)
  } else   # Empty matrix
    # mat <- matrix(0, nrow = rows, ncol = cols)
    mat <- sparseMatrix(i = c(), j = c(), dims = c(rows, cols))

  # Convert 2D ND arrays to 1D
  combo_vec <- mat_cache@const_vec + param_cache@const_vec
  const_vec <- as.vector(combo_vec)
  list(mat, -const_vec)
}

# TODO: Double-check scoping of the returned oracle function
#
# Nonlinear Constraint Matrix
#
# Returns an oracle for the nonlinear constraints. The oracle computes the combined function value, gradient, and Hessian.
#
# @param nonlin_constr A list of nonlinear constraints represented as oracle functions.
# @return An oracle function.
# @rdname nonlin_matrix-int
.nonlin_matrix <- function(object, nonlin_constr) {
  rows <- sum(sapply(nonlin_constr, function(c) { prod(dim(c)) }))
  cols <- object@sym_data@.x_length
  var_offsets <- object@sym_data@.var_offsets

  big_x <- matrix(0, nrow = cols, ncol = 1)
  for(constr in nonlin_constr)
    constr <- place_x0(constr, big_x, var_offsets)

  # Oracle for function value, gradient, and Hessian
  oracle <- function(x = NA, z = NA) {
    if(is.na(x))
      return(list(rows, big_x))
    big_f <- matrix(0, nrow = rows, ncol = 1)
    big_Df <- sparseMatrix(i = c(), j = c(), dims = c(rows, cols))
    if(!is.na(z))
      big_H <- sparseMatrix(i = c(), j = c(), dims = c(cols, cols))
    offset <- 0

    for(constr in nonlin_constr) {
      constr_entries <- prod(dim(constr))
      local_x <- extract_variables(constr, x, var_offsets)
      if(!is.na(z)) {
        tmp <- constr@f(local_x, z[offset:(offset + constr_entries)])
        f <- tmp[[1]]
        Df <- tmp[[2]]
        H <- tmp[[3]]
      } else {
        result <- constr@f(local_x)
        if(!is.na(result)) {
          f <- result[[1]]
          Df <- result[[2]]
        } else
          return(NA)
      }

      big_f[offset:(offset + constr_entries)] <- f
      big_Df <- place_Df(constr, big_Df, Df, var_offsets, offset)
      if(!is.na(z))
        big_H <- place_H(constr, big_H, H, var_offsets)
      offset <- offset + constr_entries
    }

    if(is.na(z))
      return(list(big_f, big_Df))
    return(list(big_f, big_Df, big_H))
  }

  return(oracle)
}

setClassUnion("SymDataORNull", c("SymData", "NULL"))
setClassUnion("MatrixDataORNull", c("MatrixData", "NULL"))

#
# The ProblemData class.
#
# This class is a wrapper for the symbolic and numeric data of a problem.
#
# @slot sym_data A \linkS4class{SymData} object representing the symbolic data for the problem.
# @slot matrix_data A \linkS4class{MatrixData} object representing the numeric data for the problem.
# @slot prev_result A \code{list} containing the result of the last solve.
# @rdname ProblemData-class
.ProblemData <- setClass("ProblemData", representation(sym_data = "SymDataORNull", matrix_data = "MatrixDataORNull", prev_result = "list"),
                                        prototype(sym_data = NULL, matrix_data = NULL, prev_result = list()))

#
# Problem Data Constructor
#
# Construct a \linkS4class{ProblemData} object.
#
# @param sym_data A \linkS4class{SymData} object representing the symbolic data for the problem.
# @param matrix_data A \linkS4class{MatrixData} object representing the numeric data for the problem.
# @param prev_result A \code{list} containing the result of the last solve.
# @return A \linkS4class{ProblemData} object.
# @docType methods
# @rdname ProblemData
ProblemData <- function(sym_data = NULL, matrix_data = NULL, prev_result = list()) {
  .ProblemData(sym_data = sym_data, matrix_data = matrix_data, prev_result = prev_result)
}
