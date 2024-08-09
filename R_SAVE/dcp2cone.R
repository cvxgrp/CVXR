## Uses definitions in dcpcanon.R which needs to precede this when sourcing!

#'
#' Reduce DCP Problem to Conic Form
#'
#' This reduction takes as input (minimization) DCP problems and converts them into problems
#' with affine objectives and conic constraints whose arguments are affine.
#'
#' @rdname Dcp2Cone-class
.Dcp2Cone <- setClass("Dcp2Cone", representation(quad_obj = "logical", cone_canon_methods = "list", quad_canon_methods = "list"),
                                  prototype(quad_obj = FALSE, cone_canon_methods = Dcp2Cone.CANON_METHODS, quad_canon_methods = Qp2QuadForm.QUAD_CANON_METHODS), contains = "Canonicalization")
Dcp2Cone <- function(problem = NULL, quad_obj = FALSE) { .Dcp2Cone(problem = problem, quad_obj = quad_obj) }

setMethod("initialize", "Dcp2Cone", function(.Object, ..., quad_obj = FALSE, cone_canon_methods = Dcp2Cone.CANON_METHODS, quad_canon_methods = Qp2QuadForm.QUAD_CANON_METHODS) {
  .Object <- callNextMethod(.Object, ...)
  .Object@cone_canon_methods <- cone_canon_methods
  .Object@quad_canon_methods <- quad_canon_methods
  .Object@quad_obj <- quad_obj
  .Object
})

#' @param object A \linkS4class{Dcp2Cone} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn Dcp2Cone A problem is accepted if it is a minimization and is DCP.
setMethod("accepts", signature(object = "Dcp2Cone", problem = "Problem"), function(object, problem) {
  inherits(problem@objective, "Minimize") && is_dcp(problem)
})

#' @describeIn Dcp2Cone Converts a DCP problem to a conic form.
setMethod("perform", signature(object = "Dcp2Cone", problem = "Problem"), function(object, problem) {
  if(!accepts(object, problem))
    stop("Cannot reduce problem to cone program")

  inverse_data <- InverseData(problem)

  canon <- canonicalize_tree(object, problem@objective, TRUE)
  canon_objective <- canon[[1]]
  canon_constraints <- canon[[2]]

  for(constraint in problem@constraints) {
    # canon_constr is the constraint re-expressed in terms of
    # its canonicalized arguments, and aux_constr are the constraints
    # generated while canonicalizing the arguments of the original
    # constraint
    canon <- canonicalize_tree(object, constraints, FALSE)
    canon_constr <- canon[[1]]
    aux_constr <- canon[[2]]
    canon_constraints <- c(canon_constraints, aux_constr, list(canon_constr))
    inverse_data@cons_id_map[[as.character(id(constraint))]] <- id(canon_constr)
  }

  new_problem <- Problem(canon_objective, canon_constraints)
  return(list(new_problem, inverse_data))
})

#' @describeIn Dcp2Cone Recursively canonicalize an Expression.
setMethod("canonicalize_tree", signature(object = "Dcp2Cone", expr = "Expression", affine_above = "logical"), function(object, expr, affine_above) {
  # TODO: don't copy affine expressions?
  if(inherits(expr, "PartialProblem")) {
    canon <- canonicalize_tree(object, expr@args[[1]]@objective@expr, FALSE)
    canon_expr <- canon[[1]]
    constrs <- canon[[2]]
    for(constr in expr@args[[1]]@constraints) {
      canon <- canonicalize_tree(object, constr, FALSE)
      canon_constr <- canon[[1]]
      aux_constr <- canon[[2]]
      constrs <- c(constrs, list(canon_constr), aux_constr)
    }
  } else {
    affine_atom <- !(class(expr) %in% names(object@cone_canon_methods))
    canon_args <- list()
    constrs <- list()
    for(arg in expr@args) {
      tmp <- canonicalize_tree(object, arg, affine_atom && affine_above)
      canon_arg <- tmp[[1]]
      c <- tmp[[2]]
      canon_args <- c(canon_args, list(canon_arg))
      constrs <- c(constrs, c)
    }
    tmp <- canonicalize_expr(object, expr, canon_args, affine_above)
    canon_expr <- tmp[[1]]
    c <- tmp[[2]]
    constrs <- c(constrs, c)
  }
  return(list(canon_expr, constrs))
})

#' @describeIn Dcp2Cone Canonicalize an expression wrt canonicalized arguments.
setMethod("canonicalize_expr", "Dcp2Cone", function(object, expr, args, affine_above) {
  # Constant trees are collapsed, but parameter trees are preserved.
  if(is(expr, "Expression") && (is_constant(expr) && (is.null(parameters(expr)) || length(parameters(expr)) == 0)))
    return(list(expr, list()))
  else if(object@quad_obj && affine_above && class(expr) %in% names(object@quad_canon_methods))
    return(object@quad_canon_methods[[class(expr)]](expr, args))
  else if(class(expr) %in% names(object@cone_canon_methods))
    return(object@cone_canon_methods[[class(expr)]](expr, args))
  else
    return(list(copy(expr, args), list()))
})

#'
#' Summary of Cone Dimensions Present in Constraints
#'
#' Constraints must be formatted as dictionary that maps from
#' constraint type to a list of constraints of that type.
#'
#' @slot zero The dimension of the zero cone.
#' @slot nonneg The dimension of the non-negative cone.
#' @slot exp The number of 3-dimensional exponential cones.
#' @slot soc A vector of the second-order cone dimensions.
#' @slot psd A vector of the positive semidefinite cone dimensions, where the dimension of the PSD cone of k by k matrices is k.
#' @rdname ConeDims-class
ConeDims <- setClass("ConeDims",
                     representation(zero = "integer", nonneg = "integer", exp = "integer",
                                    soc = "integer", psd = "integer", p3d = "numeric"),
                     prototype(zero = NA_integer_, nonneg = NA_integer_, exp = NA_integer_,
                               soc = NA_integer_, psd = NA_integer_, p3d = NA_real_))

setMethod("initialize", "ConeDims",
          function(.Object, constr_map) {
            .Object$zero <- Reduce(sum, lapply(constr_map$Zero, size))
            .Object$nonneg <- Reduce(sum, lapply(constr_map$NonNeg, size))
            .Object$exp <- Reduce(sum, lapply(constr_map$ExpCone, num_cones))
            .Object$soc <- unlist(lapply(constr_map$SOC, cone_sizes))
            .Object$psd <- unlist(lapply(constr_map$PSD, nrow))
            .Object$p3d <- unlist(lapply(constr_map$PowCone3D, function(x) value(x@alpha)))
            .Object
          })

setMethod("show", "ConeDims", function(object) {
  cat(
    sprintf(
      "zero: %i, nonneg: %i, exp: %i, soc: [%s], psd: [%s], p3d: [%s]\n",
      object@zero, object@nonneg, object@exp, paste(object@soc, collapse = ","),
      paste(object@psd, collapse = ","), paste(object@p3d, collapse = ",")
    )
  )
})

setMethod("as.character", "ConeDims", function(x) {
  sprintf(
    "%i equalities, %i inequalities, %i exponential cones, SOC constraints: [%s], PSD constraints: [%s], Power Cones: [%s]",
    object@zero, object@nonneg, object@exp, paste(object@soc, collapse = ","),
    paste(object@psd, collapse = ","), paste(object@p3d, collapse = ",")
  )
})

setMethod("$", signature(x = "ConeDims"), function(x, name) {
  if (name == EQ_DIM)
    return(object@zero)
  else if (name == LEQ_DIM)
    return(object@nonneg)
  else if (name == EXP_DIM)
    return(object@exp)
  else if (name == SOC_DIM)
    return(object@soc)
  else if (name == PSD_DIM)
    return(object@psd)
  else if (name == "p3")   # P3D_DIM = "p3"
    return(object@p3d)
  else
    stop("Unknown key: ", name)
})

# TODO(akshayka): unit tests
#'
#' Parametrized Cone Program
#'
#' minimize   c'x  + d + [(1/2)x'Px]
#' subject to cone_constr1(A_1*x + b_1, ...)
#' ...
#' cone_constrK(A_i*x + b_i, ...)
#'
#' The constant offsets d and b are the last column of c and A.
#'
#' @rdname ParamConeProg-class
ParamConeProg <- setClass("ParamConeProg", representation(c = "numeric", x = "Variable", A = "numeric", variables = "list", var_id_to_col = "list",
                                                          constraints = "list", parameters = "list", param_id_to_col = "list", P = "numeric", formatted = "logical",
                                                          reduced_A = "numeric", reduced_P = "numeric", constr_size = "numeric", constr_map = "list", cone_dims = "S4ORNULL",
                                                          id_to_param = "list", param_id_to_size = "list", total_param_size = "numeric", id_to_var = "list"),
                          prototype(P = NA_real_, formatted = FALSE, reduced_A = NA_real_, reduced_P = NA_real_, constr_size = NA_real_, constr_map = list(),
                                    cone_dims = NULL, id_to_param = list(), param_id_to_size = list(), total_param_size = NA_real_, id_to_var = list()), contains = "ParamProb")

setMethod("initialize", "ParamConeProg", function(.Object, ..., c, x, A, variables, var_id_to_col, constraints, parameters, param_id_to_col, P = NA_real_, formatted = FALSE,
                                                  reduced_A = NA_real_, reduced_P = NA_real_, constr_size = NA_real_, constr_map = list(), cone_dims = NULL, id_to_param = list(),
                                                  param_id_to_size = list(), total_param_size = NA_real_, id_to_var = list()) {
  # The problem data tensors; c is for the constraint, and A for
  # the problem data matrix
  .Object@c <- c
  .Object@A <- A
  .Object@P <- P
  # The variable
  .Object@x <- x

  # Form a reduced representation of A and P, for faster application
  # of parameters.
  .Object@reduced_A = ReducedMat(.Object@A, size(.Object@x))
  .Object@reduced_P = ReducedMat(.Object@P, size(.Object@x), quad_form = TRUE)

  .Object@constraints <- constraints
  .Object@constr_size <- sum(sapply(constraints, size))
  .Object@constr_map <- group_constraints(constraints)
  .Object@cone_dims <- ConeDims(.Object@constr_map)
  .Object@parameters <- parameters
  .Object@param_id_to_col <- param_id_to_col

  .Object@id_to_param <- list()
  .Object@param_id_to_size <- list()
  .Object@total_param_size <- 0
  for(p in .Object@parameters) {
    p_id_char <- as.character(id(p))
    p_size <- size(p)
    .Object@id_to_param[[p_id_char]] <- p
    .Object@param_id_to_size[[p_id_char]] <- p_size
    .Object@total_param_size <- .Object@total_param_size + p_size
  }

  # TODO technically part of inverse data.
  .Object@variables <- variables
  .Object@var_id_to_col <- var_id_to_col
  .Object@id_to_var <- list()
  for(v in .Object@variables)
    .Object@id_to_var[[as.character(id(v))]] <- v

  # whether this param cone prog has been formatted for a solver
  .Object@formatted <- formatted
  .Object
})

## Add ParamConeProg to class union ParamProgORNULL
setIs("ParamConeProg", "ParamProgORNULL")

#' @param object A \linkS4class{ParamConeProg} object.
#' @describeIn ParamConeProg Is the problem mixed-integer?
setMethod("is_mixed_integer", "ParamConeProg", function(object) {
  return(object@x@attributes$boolean || object@x@attributes$integer)
})

#' @param id_to_param_value: (optional) dict mapping parameter ids to values.
#' @param zero_offset: (optional) if True, zero out the constant offset in the parameter vector.
#' @param keep_zeros: (optional) if True, store explicit zeros in A where parameters are affected.
#' @param quad_obj: (optional) if True, include quadratic objective term.
#' @describeIn ParamConeProg Returns A, b after applying parameters (and reshaping)
setMethod("apply_parameters", "ParamConeProg", function(object, id_to_param_value = NULL, zero_offset = FALSE, keep_zeros = FALSE, quad_obj = FALSE) {
  cache(object@reduced_A, keep_zeros)

  param_value <- function(idx) {
    if(is.null(id_to_param_value))
      return(value(object@id_to_param[[idx]]))
    else
      return(id_to_param_value[[idx]])
  }

  param_vec <- canonInterface.get_parameter_vector(object@total_param_size, object@param_id_to_col, object@param_id_to_size, param_value, zero_offset = zero_offset)
  cd = canonInterface.get_matrix_from_tensor(object@c, param_vec, size(object@x), with_offset = TRUE)
  c <- cd[[1]]
  d <- cd[[2]]
  c <- as.vector(c)

  Ab <- get_matrix_from_tensor(object@reduced_A, param_vec, with_offset = TRUE)
  A <- Ab[[1]]
  b <- Ab[[2]]

  if(quad_obj) {
    cache(object@reduced_P, keep_zeros)
    P <- get_matrix_from_tensor(object@reduced_P, param_vec, with_offset = FALSE)[[1]]
    return(list(P = p, c = c, d = d, A = A, b = matrix(b)))
  } else
    return(list(c = c, d = d, A = A, b = matrix(b)))
})

#' @describeIn ParamConeProg Multiplies by Jacobian of parameter mapping. Assumes delA is sparse. Returns a mapping of param id to dparam.
setMethod("apply_param_jac", "ParamConeProg", function(object, delc, delA, delb, active_params = NULL) {
  if(!(is.na(object@P) || is.null(object@P)))
    stop("Can't apply Jacobian with a quadratic objective")

  if(is.null(active_params))
    active_params <- sapply(object@parameters, id)

  if(length(object@c) == 1)
    del_param_vec <- c()
  else
    del_param_vec <- delc %*% object@c[1:(length(object@c) - 1)]

  flatdelA <- delA
  dim(flatdelA) <- c(prod(dim(delA)), 1)
  sparsedelb <- Matrix(delb, sparse = TRUE)
  delAb <- rbind(flatdelA, sparsedelb)

  one_gig_of_doubles <- 125000000
  if(nrow(delAb) < one_gig_of_doubles) {
    # fast path: if delAb is small enough, just materialize it
    # in memory because sparse-matrix @ dense vector is much faster
    # than sparse @ sparse
    del_param_vec <- del_param_vec + t(object@A) %*% matrix(delAb)
  } else {
    # slow path.
    # TODO: make this faster by intelligently operating on the
    # sparse matrix data / making use of reduced_A
    del_param_vec <- del_param_vec + matrix(t(delAb) %*% object@A)
  }
  del_param_vec <- as.vector(del_param_vec)

  param_id_to_delta_param <- list()
  for(param_id in names(object@param_id_to_col)) {
    col <- object@param_id_to_col[[param_id]]
    if(param_id %in% names(active_params)) {
      param <- object@id_to_param[[param_id]]
      delta <- del_param_vec[col:(col + size(param))]
      param_id_to_delta_param[[param_id]] <- matrix(delta, nrow = nrow(param), ncol = ncol(param), byrow = FALSE)
    }
  }

  return(param_id_to_delta_param)
})

#' @describeIn ParamConeProg Splits the solution into individual variables.
setMethod("split_solution", "ParamConeProg", function(object, sltn, active_vars = NULL) {
  if(is.null(active_vars))
    active_vars <- sapply(object@variables, function(v) { as.character(id(v)) })
  # var id to solution.
  sltn_dict <- list()
  for(var_id in names(object@var_id_to_col)) {
    col <- object@var_id_to_col[[var_id]]
    if(var_id %in% active_vars) {
      var <- object@id_to_var[[var_id]]
      value <- sltn[col:(size(var) + col)]
      if(attributes_were_lowered(var)) {
        orig_var <- variable_of_provenance(var)
        value <- recover_value_for_variable(orig_var, value, project = FALSE)
        sltn_dict[[as.character(id(orig_var))]] = matrix(value, nrow = nrow(orig_var), ncol = ncol(orig_var), byrow = FALSE)
      } else
        sltn_dict[[var_id]] = matrix(value, nrow = nrow(var), ncol = ncol(var), byrow = FALSE)
    }
  }
  return(sltn_dict)
})

#' @describeIn ParamConeProg Adjoint of split solution.
setMethod("split_adjoint", "ParamConeProg", function(object, del_vars = NULL) {
  var_vec <- rep(0, size(object@x))
  if(is.null(del_vars))
    return(var_vec)

  for(var_id in names(del_vars)) {
    delta <- del_vars[[var_id]]
    var <- object@id_to_var[[var_id]]
    col <- object@var_id_to_col[[var_id]]
    if(attributes_were_lowered(var)) {
      orig_var <- variable_of_provenance(var)
      if(attributes_present(list(orig_var), SYMMETRIC_ATTRIBUTES))
        delta <- delta + t(delta) - diag(diag(delta))
      delta <- lower_value(orig_var, delta)
    }
    var_vec[col:(col + size(var))] <- as.vector(delta)
  }
  return(var_vec)
})

#'
#' Construct Matrices for Linear Cone Problems
#'
#' Linear cone problems are assumed to have a linear objective and cone constraints,
#' which may have zero or more arguments, all of which must be affine.
#'
#' minimize c^Tx
#' subject to cone_constr1(A_1*x + b_1, ...)
#'            ...
#'            cone_constrK(A_K*x + b_K, ...)
#'
#' @rdname ConeMatrixStuffing-class
ConeMatrixStuffing <- setClass("ConeMatrixStuffing", representation(quad_obj = "logical", canon_backend = "character"),
                               prototype(quad_obj = FALSE, canon_backend = NA_character_), contains = "MatrixStuffing")

#' @param object A \linkS4class{ConeMatrixStuffing} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn ConeMatrixStuffing Is the solver accepted?
setMethod("accepts", signature(object = "ConeMatrixStuffing", problem = "Problem"), function(object, problem) {
  valid_obj_curv <- (object@quad_obj && is_quadratic(problem@objective@expr)) || is_affine(problem@objective@expr)
  return(inherits(problem@objective, "Minimize") &&
         valid_obj_curv &&
         length(convex_attributes(variables(problem))) == 0 &&
         are_args_affine(problem@constraints) &&
         is_dpp(problem))
})

#' @param extractor Used to extract the affine coefficients of the objective.
#' @describeIn ConeMatrixStuffing Returns a list of the stuffed matrices
setMethod("stuffed_objective", signature(object = "ConeMatrixStuffing", problem = "Problem", extractor = "CoeffExtractor"), function(object, problem, extractor) {
  # concatenate all variables in one vector
  boolint <- extract_mip_idx(variables(problem))
  boolean <- boolint[[1]]
  integer <- boolint[[2]]
  # x <- Variable(extractor@x_length, boolean = boolean, integer = integer)
  x <- Variable(extractor@x_length, 1, boolean = boolean, integer = integer)

  if(object@quad_obj) {
    # extract to 0.5 %*% t(x) %*% P %*% x + t(q) %*% x + r
    expr <- problem@objective@expr
    tmp <- quad_form(extractor, expr)
    params_to_P <- tmp[[1]]
    params_to_c <- tmp[[2]]
    # Handle 0.5 factor.
    params_to_P <- 2*params_to_P
  } else {
    # Extract to t(c) %*% x + r; c is represented by a ma
    params_to_c <- affine(extractor, problem@objective@expr)
    params_to_P <- NULL
  }
  return(list(params_to_P, params_to_c, x))
})

#' @describeIn ConeMatrixStuffing Constructs matrices and returns a ParamConeProg (parametrized cone program)
setMethod("perform", signature(object = "ConeMatrixStuffing", problem = "Problem"), function(object, problem) {
  inverse_data <- InverseData(problem)
  # Form the constraints
  extractor <- CoeffExtractor(inverse_data, object@canon_backend)
  stuffed <- stuffed_objective(object, problem, extractor)
  params_to_P <- stuffed[[1]]
  params_to_c <- stuffed[[2]]
  flattened_variable <- stuffed[[3]]

  # Lower equality and inequality to Zero and NonNeg.
  cons <- list()
  for(con in problem@constraints) {
    if(is(con, "EqConstraint"))
      con <- lower_equality(con)
    else if(is(con, "IneqConstraint"))
      con <- lower_ineq_to_nonneg(con)
    else if(is(con, "NonPosConstraint"))
      con <- nonpos2nonneg(con)
    else if(is(con, "SOC") && con@axis == 1)
      con <- SOC(con@args[[1]], t(con@args[[2]]), axis = 2, constr_id = con@constr_id)
    else if(is(con, "PowCone3D") && ndim(con@args[[1]]) > 1) {
      xyz <- con@args
      x <- xyz[[1]]
      y <- xyz[[2]]
      z <- xyz[[3]]
      alpha <- con@alpha
      con <- PowCone3D(flatten(x), flatten(y), flatten(z), flatten(alpha), constr_id = con@constr_id)
    } else if(is(con, "ExpCone") && ndim(con@args[[1]]) > 1) {
      xyz <- con@args
      x <- xyz[[1]]
      y <- xyz[[2]]
      z <- xyz[[3]]
      con <- ExpCone(flatten(x), flatten(y), flatten(z), constr_id = con@constr_id)
    }
    cons <- c(cons, con)
  }

  # Reorder constraints to Zero, NonNeg, SOC, PSD, EXP, PowCone3D
  constr_map <- group_constraints(cons)
  ordered_cons <- c(constr_map$Zero, constr_map$NonNeg, constr_map$SOC, constr_map$PSD,
                    constr_map$ExpCone, constr_map$PowCone3D)
  inverse_data@cons_id_map <- list()
  for(con in ordered_cons) {
    con_id <- id(con)
    inverse_data@cons_id_map[[as.character(con_id)]] <- con_id
  }

  inverse_data@constraints <- ordered_cons
  # Batch expressions together, then split apart.
  for(c in ordered_cons) {
    for(arg in c@args)
      expr_list <- c(expr_list, arg)
  }
  params_to_problem_data <- affine(extractor, expr_list)

  inverse_data@minimize <- (inherits(problem@objective, "Minimize"))
  new_prob <- ParamConeProg(params_to_c,
                            flattened_variable,
                            params_to_problem_data,
                            variables(problem),
                            inverse_data@var_offsets,
                            ordered_cons,
                            parameters(problem),
                            inverse_data@param_id_map,
                            P = params_to_P)
  return(list(new_prob, inverse_data))
})

#' @param solution A \linkS4class{Solution} object to invert.
#' @param inverse_data A \linkS4class{InverseData} object containing data necessary for the inversion.
#' @describeIn ConeMatrixStuffing Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "ConeMatrixStuffing", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  var_map <- inverse_data@var_offsets
  con_map <- inverse_data@cons_id_map
  # Flip sign of opt val if maximize.
  opt_val <- solution@opt_val
  if(!(solution@status %in% ERROR) && !inverse_data@minimize)
    opt_val <- -solution@opt_val

  primal_vars <- list()
  dual_vars <- list()
  if(!(solution@status %in% SOLUTION_PRESENT))
    return(Solution(solution@status, opt_val, primal_vars, dual_vars, solution@attr))

  # Split vectorized variable into components.
  x_opt <- values(solution@primal_vars)[[1]]
  for(var_id in names(var_map)) {
    offset <- var_map[[var_id]]
    var_dim <- inverse_data@var_shapes[[var_id]]
    var_size <- as.integer(prod(var_dim))
    primal_vars[[var_id]] = matrix(x_opt[offset:(offset + var_size)], nrow = var_dim[1], ncol = var_dim[2], byrow = FALSE)
  }

  if(!is.null(solution@dual_vars) && length(solution@dual_vars) > 0) {
    for(old_con in names(con_map)) {
      new_con <- con_map[[old_con]]
      con_obj <- inverse_data@id2cons[[old_con]]
      con_obj_dim <- dim(con_obj)
      # TODO rationalize Exponential.
      if(all(con_obj_dim == c(1,1)) || is(con_obj, "ExpCone") || is(con_obj, "SOC"))
        dual_vars[[old_con]] <- solution@dual_vars[[new_con]]
      else
        dual_vars[[old_con]] <- matrix(solution@dual_vars[[new_con]], nrow = con_obj_dim[1], ncol = con_obj_dim[2], byrow = FALSE)
    }
  }

  return(Solution(solution@status, opt_val, primal_vars, dual_vars, solution@attr))
})

# TODO: ConeMatrixStuffing invert function.

