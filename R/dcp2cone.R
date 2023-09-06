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
ConeDims <- setClass("ConeDims", representation(zero = "numeric", nonneg = "numeric", exp = "numeric", soc = "numeric", psd = "numeric", p3d = "numeric"),
                                 prototype(zero = NA_integer_, nonneg = NA_integer_, exp = NA_integer_, soc = NA_integer_, psd = NA_integer_, p3d = NA_real_))

setMethod("initialize", "ConeDims", function(.Object, zero = NA_integer_, nonneg = NA_integer_, exp = NA_integer_, soc = NA_integer_, psd = NA_integer_, p3d = list()) {
  .Object@zero <- as.integer(sum(sapply(constr_map$Zero, size)))
  .Object@nonneg <- as.integer(sum(sapply(constr_map$NonNeg, size)))
  .Object@exp <- as.integer(sum(sapply(constr_map$ExpCone, num_cones)))
  
  .Object@soc <- c()
  for(con in constr_map$SOC)
    .Object@soc <- c(.Object@soc, cone_sizes(con))
  .Object@soc <- as.integer(.Object@soc)
  .Object@psd <- as.integer(sapply(constr_map$PSD, nrow))
  
  p3d <- c()
  if(!is.null(constr_map$PowCone3D)) {
    for(con in constr_map$PowCone3D)
      p3d <- c(p3d, value(con@alpha))
  }
  .Object@p3d <- p3d
  .Object
})

setMethod("show", "ConeDims", function(object) {
  cat(paste("zero:", object@zero), paste("nonneg:", object@nonneg), paste("exp:", object@exp),
      paste("soc:", object@soc), paste("psd:", object@psd), paste("p3d:", object@p3d), sep = ", ")
})

setMethod("as.character", "ConeDims", function(x) {
  sprintf("%i equalities, %i inequalities, %i exponential cones,\nSOC constraints: %s, PSD constraints: %s,\n3d power cones %s.", 
          x@zero, x@nonneg, x@exp, 
          paste("(", paste(x@soc, sep = ", "), ")", sep = ""), 
          paste("(", paste(x@psd, sep = ", "), ")", sep = ""), 
          paste("(", paste(x@p3d, sep = ", "), ")", sep = ""))
})

setMethod("$", signature(x = "ConeDims"), function(x, name) {
  if(name == EQ_DIM)
    return(object@zero)
  else if(name == LEQ_DIM)
    return(object@nonneg)
  else if(name == EXP_DIM)
    return(object@exp)
  else if(name == SOC_DIM)
    return(object@soc)
  else if(name == PSD_DIM)
    return(object@psd)
  else if(name == "p3")   # P3D_DIM = "p3"
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
  c, d = canonInterface.get_matrix_from_tensor(object@c, param_vec, size(object@x), with_offset = TRUE)
  c <- cd[[1]]
  d <- cd[[2]]
  c <- as.vector(c)
  
  Ab <- get_matrix_from_tensor(object@reduced_A, param_vec, with_offset = TRUE)
  A <- Ab[[1]]
  b <- Ab[[2]]
  
  if(quad_obj) {
    cache(object@reduced_P, keep_zeros)
    P <- get_matrix_from_tensor(object@reduced_P, param_vec, with_offset = FALSE)[[1]]
    return(list(P, c, d, A, matrix(b)))
  } else
    return(list(c, d, A, matrix(b)))
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

# Atom canonicalizers.
#' 
#' Dcp2Cone canonicalizer for the entropy atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from an entropy atom where 
#' the objective function is just the variable t with an ExpCone constraint.
Dcp2Cone.entr_canon <- function(expr, args) {
  x <- args[[1]]
  expr_dim <- dim(expr)
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = expr_dim)

  # -x*log(x) >= t is equivalent to x*exp(t/x) <= 1
  # TODO: ExpCone requires each of its inputs to be a Variable; is this something we want to change?
  if(is.null(expr_dim))
    ones <- Constant(1)
  else
    ones <- Constant(matrix(1, nrow = expr_dim[1], ncol = expr_dim[2]))
  constraints <- list(ExpCone(t, x, ones))
  return(list(t, constraints))
}

#' 
#' Dcp2Cone canonicalizer for the exponential atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from an exponential atom 
#' where the objective function is the variable t with an ExpCone constraint.
Dcp2Cone.exp_canon <- function(expr, args) {
  expr_dim <- dim(expr)
  x <- promote(args[[1]], expr_dim)
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = expr_dim)
  if(is.null(expr_dim))
    ones <- Constant(1)
  else
    ones <- Constant(matrix(1, nrow = expr_dim[1], ncol = expr_dim[2]))
  constraints <- list(ExpCone(x, ones, t))
  return(list(t, constraints))
}

#' 
#' Dcp2Cone canonicalizer for the geometric mean atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a geometric mean atom 
#' where the objective function is the variable t with geometric mean constraints
Dcp2Cone.geo_mean_canon <- function(expr, args) {
  x <- args[[1]]
  w <- expr@w
  expr_dim <- dim(expr)
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = expr_dim)

  x_list <- lapply(seq_len(w), function(i) { x[i] })

  # TODO: Catch cases where we have (0,0,1)?
  # TODO: What about curvature case (should be affine) in trivial case of (0,0,1)?
  # Should this behavior match with what we do in power?
  return(list(t, gm_constrs(t, x_list, w)))
}

#' 
#' Dcp2Cone canonicalizer for the huber atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a huber atom where the objective 
#' function is the variable t with square and absolute constraints
Dcp2Cone.huber_canon <- function(expr, args) {
  M <- expr@M
  x <- args[[1]]
  expr_dim <- dim(expr)
  # n <- Variable(expr_dim)
  # s <- Variable(expr_dim)
  n <- new("Variable", dim = expr_dim)
  s <- new("Variable", dim = expr_dim)

  # n^2 + 2*M*|s|
  # TODO: Make use of recursion inherent to canonicalization process and just return a
  # power/abs expression for readiability's sake
  power_expr <- power(n,2)
  canon <- Dcp2Cone.power_canon(power_expr, power_expr@args)
  n2 <- canon[[1]]
  constr_sq <- canon[[2]]

  abs_expr <- abs(s)
  canon <- EliminatePwl.abs_canon(abs_expr, abs_expr@args)
  abs_s <- canon[[1]]
  constr_abs <- canon[[2]]

  obj <- n2 + 2*M*abs_s

  # x == s + n
  constraints <- c(constr_sq, constr_abs)
  constraints <- c(constraints, x == s + n)
  return(list(obj, constraints))
}

#' 
#' Dcp2Cone canonicalizer for the indicator atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from an indicator atom and
#' where 0 is the objective function with the given constraints
#' in the function.
Dcp2Cone.indicator_canon <- function(expr, args) {
  return(list(Constant(0), args))
}

#' 
#' Dcp2Cone canonicalizer for the KL Divergence atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a KL divergence atom
#' where t is the objective function with the ExpCone constraints.
Dcp2Cone.kl_div_canon <- function(expr, args) {
  expr_dim <- dim(expr)
  x <- promote(args[[1]], expr_dim)
  y <- promote(args[[2]], expr_dim)
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = expr_dim)
  constraints <- list(ExpCone(t, x, y))
  obj <- y - x - t
  return(list(obj, constraints))
}

#' 
#' Dcp2Cone canonicalizer for the lambda maximization atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a lambda maximization atom
#' where t is the objective function and a PSD constraint and a
#' constraint requiring I*t to be symmetric.
Dcp2Cone.lambda_max_canon <- function(expr, args) {
  A <- args[[1]]
  n <- nrow(A)
  t <- Variable()
  prom_t <- promote(t, c(n,1))
  # Constrain I*t - A to be PSD; note this expression must be symmetric
  tmp_expr <- DiagVec(prom_t) - A
  constr <- list(PSDConstraint(tmp_expr))
  if(!is_symmetric(A)) {
    ut <- UpperTri(A)
    lt <- UpperTri(t(A))
    constr <- c(constr, ut == lt)
  }
  return(list(t, constr))
}

#' 
#' Dcp2Cone canonicalizer for the largest lambda sum atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a lambda sum of the k
#' largest elements atom where k*t + trace(Z) is the objective function.
#' t denotes the variable subject to constraints and Z is a PSD matrix variable
#' whose dimensions consist of the length of the vector at hand. The constraints
#' require the the diagonal matrix of the vector to be symmetric and PSD.
Dcp2Cone.lambda_sum_largest_canon <- function(expr, args) {
  # S_k(X) denotes lambda_sum_largest(X, k)
  # t >= k S_k(X - Z) + trace(Z), Z is PSD
  # implies
  # t >= ks + trace(Z)
  # Z is PSD
  # sI >= X - Z (PSD sense)
  # which implies
  # t >= ks + trace(Z) >= S_k(sI + Z) >= S_k(X)
  # We use the fact that
  # S_k(X) = sup_{sets of k orthonormal vectors u_i}\sum_{i}u_i^T X u_i
  # and if Z >= X in PSD sense then
  # \sum_{i}u_i^T Z u_i >= \sum_{i}u_i^T X u_i
  #
  # We have equality when s = lambda_k and Z diagonal
  #  with Z_{ii} = (lambda_i - lambda_k)_+

  X <- expr@args[[1]]
  k <- expr@k
  # Z <- Variable(c(nrow(X), nrow(X)), PSD = TRUE)
  Z <- Variable(nrow(X), ncol(X), PSD = TRUE)
  canon <- Dcp2Cone.lambda_max_canon(expr, list(X - Z))
  obj <- canon[[1]]
  constr <- canon[[2]]
  obj <- k*obj + matrix_trace(Z)
  return(list(obj, constr))
}

#' 
#' Dcp2Cone canonicalizer for the log 1p atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a log 1p atom where
#' t is the objective function and the constraints consist of
#' ExpCone constraints + 1.
Dcp2Cone.log1p_canon <- function(expr, args) {
  return(Dcp2Cone.log_canon(expr, list(args[[1]] + 1)))
}

#' 
#' Dcp2Cone canonicalizer for the logistic function atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from the logistic atom
#' where the objective function is given by t0 and the 
#' constraints consist of the ExpCone constraints.
Dcp2Cone.logistic_canon <- function(expr, args) {
  x <- args[[1]]
  expr_dim <- dim(expr)
  # log(1 + exp(x)) <= t is equivalent to exp(-t) + exp(x - t) <= 1
  # t0 <- Variable(expr_dim)
  t0 <- new("Variable", dim = expr_dim)
  canon1 <- Dcp2Cone.exp_canon(expr, list(-t0))
  canon2 <- Dcp2Cone.exp_canon(expr, list(x - t0))
  
  t1 <- canon1[[1]]
  constr1 <- canon1[[2]]
  t2 <- canon2[[1]]
  constr2 <- canon2[[2]]
  
  if(is.null(expr_dim))
    ones <- Constant(1)
  else
    ones <- Constant(matrix(1, nrow = expr_dim[1], ncol = expr_dim[2]))
  constraints <- c(constr1, constr2, list(t1 + t2 <= ones))
  return(list(t0, constraints))
}

#' 
#' Dcp2Cone canonicalizer for the log atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a log atom where
#' t is the objective function and the constraints consist of
#' ExpCone constraints
Dcp2Cone.log_canon <- function(expr, args) {
  x <- args[[1]]
  expr_dim <- dim(expr)
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = expr_dim)
  if(is.null(expr_dim))
    ones <- Constant(1)
  else
    ones <- Constant(matrix(1, nrow = expr_dim[1], ncol = expr_dim[2]))
  # TODO: ExpCone requires each of its inputs to be a Variable; is this something that we want to change?
  constraints <- list(ExpCone(t, ones, x))
  return(list(t, constraints))
}

#' 
#' Dcp2Cone canonicalizer for the log determinant atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a log determinant atom where
#' the objective function is the sum of the log of the vector D
#' and the constraints consist of requiring the matrix Z to be
#' diagonal and the diagonal Z to equal D, Z to be upper triangular
#' and DZ; t(Z)A to be positive semidefinite, where A is a n by n
#' matrix.
Dcp2Cone.log_det_canon <- function(expr, args) {
  # Reduces the atom to an affine expression and list of constraints.
  #
  # Creates the equivalent problem::
  #
  # maximize    sum(log(D[i, i]))
  # subject to: D diagonal
  # diag(D) = diag(Z)
  # Z is upper triangular.
  # [D Z; t(Z) A] is positive semidefinite
  #
  # The problem computes the LDL factorization:
  #
  # A = (Z^TD^{-1})D(D^{-1}Z)
  #
  # This follows from the inequality:
  #
  # \det(A) >= \det(D) + \det([D, Z; Z^T, A])/\det(D) >= \det(D)
  #
  # because (Z^TD^{-1})D(D^{-1}Z) is a feasible D, Z that achieves
  # det(A) = det(D) and the objective maximizes det(D).
  #
  # Parameters
  # ----------
  # expr : log_det
  # args : list of arguments for the expression
  #
  # Returns
  # -------
  # (Variable for objective, list of constraints)
  
  A <- args[[1]]  # n by n matrix.
  n <- nrow(A)
  
  z_small <- Variable(floor(n*(n+1)/2))
  Z <- UpperTri.vec_to_upper_tri(z_small, strict = FALSE)
  d_small <- DiagMat(Z)  # a vector
  D <- DiagVec(d_small)  # a matrix
  X <- bmat(list(list(D, Z),
                 list(t(Z), A)))
  
  constraints <- list(PSDConstraint(X))
  log_expr <- log(d)
  canon <- Dcp2Cone.log_canon(log_expr, log_expr@args)
  obj <- canon[[1]]
  constr <- canon[[2]]
  constraints <- c(constraints, constr)
  return(list(sum(obj), constraints))
}

#' 
#' Dcp2Cone canonicalizer for the log sum of the exp atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from the log sum
#' of the exp atom where the objective is the t variable
#' and the constraints consist of the ExpCone constraints and
#' requiring t to be less than a matrix of ones of the same size.
Dcp2Cone.log_sum_exp_canon <- function(expr, args) {
  x <- args[[1]]
  x_dim <- dim(x)
  expr_dim <- dim(expr)
  axis <- expr@axis
  keepdims <- expr@keepdims
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = expr_dim)
  
  # log(sum(exp(x))) <= t <=> sum(exp(x-t)) <= 1.
  if(is.na(axis))   # shape = c(1,1)
    promoted_t <- promote(t, x_dim)
  else if(axis == 2)   # shape = c(1,n)
    promoted_t <- Constant(matrix(1, nrow = x_dim[1], ncol = 1)) %*% reshape_expr(t, c(1, x_dim[2:length(x_dim)])))
  else   # shape = c(m,1)
    promoted_t <- reshape_expr(t, c(x_dim[1:(length(x_dim)-1)], 1)) %*% Constant(matrix(1, nrow = 1, ncol = x_dim[2]))
  
  exp_expr <- Exp(x - promoted_t)
  canon <- Dcp2Cone.exp_canon(exp_expr, exp_expr@args)
  obj <- sum_entries(canon[[1]], axis = axis, keepdims = keepdims)
  if(is.null(expr_dim))
    ones <- Constant(1)
  else
    ones <- Constant(matrix(1, nrow = expr_dim[1], ncol = expr_dim[2]))
  constraints <- c(canon[[2]], obj <= ones)
  return(list(t, constraints))
}

#' 
#' Dcp2Cone canonicalizer for the matrix fraction atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from the matrix fraction
#' atom, where the objective function is the trace of Tvar, a 
#' m by m matrix where the constraints consist of the matrix of
#' the Schur complement of Tvar to consist of P, an n by n, given
#' matrix, X, an n by m given matrix, and Tvar.
Dcp2Cone.matrix_frac_canon <- function(expr, args) {
  X <- args[[1]]   # n by m matrix
  P <- args[[2]]   # n by n matrix

  if(length(dim(X)) == 1)
    X <- reshape_expr(X, c(nrow(X), 1))
  X_dim <- dim(X)
  n <- X_dim[1]
  m <- X_dim[2]
  
  Tvar <- Variable(m, m, symmetric = TRUE)
  # A matrix with Schur complement Tvar - t(X) %*% inv(P) %*% X.
  M <- bmat(list(list(P, X),
                 list(t(X), Tvar)))
  constraints <- list(PSDConstraint(M))
  
  if(!is_symmetric(P)) {
    ut <- upper_tri(P)
    lt <- upper_tri(t(P))
    constraints <- c(constraints, ut == lt)
  }
  return(list(matrix_trace(Tvar), constraints))
}

#' 
#' Dcp2Cone canonicalizer for the multiplication of expressions atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from the multiplication of expressions atom.
Dcp2Cone.mul_canon <- function(expr, args) {
  # TODO(akshayka): expose as a reduction for user's convenience
  
  # Only allow param * var (not var * param). Associate right to left.
  # TODO: Only descend if both sides have parameters
  lhs <- args[[1]]
  rhs <- args[[2]]
  lhs_parms <- parameters(lhs)
  rhs_parms <- parameters(rhs)
  if((is.null(lhs_parms) || length(lhs_parms) == 0) && (is.null(rhs_parms) && length(rhs_parms) == 0))
    return(list(copy(expr, args), list()))
  
  op_type <- class(expr)
  if(length(variables(lhs)) > 0) {
    if(dpp_scope()) {   # TODO: Implement DPP scoping with global variables.
      if(!is_affine(rhs))
        stop("rhs must be affine if DPP")
    }
    t <- new("Variable", dim = dim(lhs))
    return(list(do.call(op_type, list(t, rhs)), list(t == lhs)))
  } else if(length(variables(rhs))) {
    if(dpp_scope()) {
      if(!is_affine(lhs))
        stop("lhs must be affine if DPP")
    }
    t <- new("Variable", dim = dim(rhs))
    return(list(do.call(op_type, list(lhs, t)), list(t == rhs)))
  }
  
  # Neither side has variables. One side must be affine in parameters.
  lhs_affine <- FALSE
  rhs_affine <- FALSE
  if(dpp_scope()) {
    lhs_affine <- is_affine(lhs)
    rhs_affine <- is_affine(rhs)
  }
  if(!(lhs_affine || rhs_affine))
    stop("Either lhs or rhs must be affine in parameters")
  
  if(lhs_affine) {
    t <- new("Variable", dim = dim(rhs))
    return(list(lhs %*% t, list(t == rhs)))
  } else {
    t <- new("Variable", dim = dim(lhs))
    return(list(t %*% rhs, list(t == lhs)))
  }
}

#' 
#' Dcp2Cone canonicalizer for the nuclear norm atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a nuclear norm atom,
#' where the objective function consists of .5 times the trace of
#' a matrix X of size m+n by m+n where the constraint consist of
#' the top right corner of the matrix being the original matrix.
Dcp2Cone.normNuc_canon <- function(expr, args) {
  A <- args[[1]]
  A_dim <- dim(A)
  m <- A_dim[1]
  n <- A_dim[2]

  # Create the equivalent problem:
  #   minimize (trace(U) + trace(V))/2
  #   subject to:
  #            [U A; t(A) V] is positive semidefinite
  U <- Variable(m, m, symmetric = TRUE)
  V <- Variable(n, n, symmetric = TRUE)
  X <- bmat(list(list(U, A),
                 list(t(A), V)))
  constraints <- list(X %>>% 0)
  trace_value <- 0.5 * (matrix_trace(U) + matrix_trace(V))
  return(list(trace_value, constraints))
}

Dcp2Cone.perspective_canon <- function(expr, args) {
  # Only working for minimization right now.
  
  if(is_convex(expr@f))
    aux_prob <- Problem(Minimize(expr@f))
  else
    aux_prob <- Problem(Maximize(expr@f))
  # Does numerical solution value of epigraph t coincide with expr.f numerical
  # value at opt?
  solver_opts <- list("use_quad_obj" = FALSE)
  chain <- .construct_chain(aux_prob, solver_opts=solver_opts, ignore_dpp=TRUE)
  num_red <- length(chain@reductions)
  if(num_red > 1)
    chain@reductions = chain@reductions[[seq(num_red - 1)]]  # skip solver reduction
  prob_canon <- apply(chain, aux_prob)[[1]]  # grab problem instance
  # get cone representation of c, A, and b for some problem.
  
  c <- as.vector(prob_canon@c)
  c_len <- length(c)
  c <- c[seq_len(c_len - 1)]
  d <- as.vector(prob_canon@c)
  d <- d[c_len]
  
  Ab <- matrix(prob_canon@A, ncol = length(c) + 1, byrow = FALSE)
  Ab_ncol <- ncol(Ab)
  A <- Ab[, seq_len(Ab_ncol - 1)]
  b <- Ab[, Ab_ncol]
  
  # given f in epigraph form, aka epi f = \{(x,t) | f(x) \leq t\}
  # = \{(x,t) | Fx +tg + e \in K} for K a cone, the epigraph of the
  # perspective, \{(x,s,t) | sf(x/s) \leq t} = \{(x,s,t) | Fx + tg + se \in K\}
  # If I have the problem "minimize f(x)" written in the CVXPY compatible
  # "c^Tx, Ax+b \in K" form, I can re-write this in the graph form above via
  # x,t \in \epi f iff Ax + b \in K and t-c^Tx \in R_+ which I can further write
  # with block matrices as Fx + tg + e \in K \times R_+
  # with F = [A ], g = [0], e = [b]
  #          [-c]      [1]      [-d]
  
  # Actually, all we need is Ax + 0*t + sb \in K, -c^Tx + t - ds >= 0
  
  t <- Variable()
  s <- args[[1]]
  x_canon <- prob_canon@x
  constraints <- list()
  
  if(!is.null(dim(A))) {
    # Rules out the case where f is affine and requires no additional
    # constraints.
    x_pers <- A@x_canon + s*b
  
    i <- 1
    for(con in prob_canon@constraints) {
      sz <- size(con)
      var_slice <- x_pers[seq(i, i + sz - 1)]
      pers_constraint <- form_cone_constraint(var_slice, con)
      constraints <- c(constraints, list(pers_constraint))
      i <- i + sz
    }
  }
  
  constraints <- c(constraints, list(-c@x_canon + t - s*d >= 0))
    
  # recover initial variables
  
  end_inds <- c(sort(unlist(prob_canon@var_id_to_col), decreasing = FALSE), nrow(x_canon) + 1)
  
  for(var in variables(expr@f)) {
    start_ind <- prob_canon@var_id_to_col[[as.character(id(var))]]
    end_ind <- end_inds[which(end_inds == start_ind) + 1]
    
    if(var@attributes$diag)  # checking for diagonal first because diagonal is also symmetric
      constraints <- c(constraints, list(diag(var) == x_canon[start_ind:(end_ind - 1)]))
    else if(is_symmetric(var) && size(var) > 1) {
      n <- nrow(var)
      inds <- which(upper.tri(matrix(1, nrow = n, ncol = n), diag = TRUE), arr.ind = TRUE)  # includes diagonal
      constraints <- c(constraints, list(var[inds] == x_canon[start_ind:(end_ind - 1)]))
    } else
      constraints <- c(constraints, list(vec(var) == x_canon[start_ind:(end_ind - 1)]))
  }
  
  if(is_convex(expr@f))
    return(list(t, constraints))
  else
    return(list(-t, constraints))
}

#' 
#' Dcp2Cone canonicalizer for the p norm atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a pnorm atom, where
#' the objective is a variable t of dimension of the original
#' vector in the problem and the constraints consist of geometric
#' mean constraints.
Dcp2Cone.pnorm_canon <- function(expr, args) {
  x <- args[[1]]
  p <- expr@p
  axis <- expr@axis
  expr_dim <- dim(expr)
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = expr_dim)

  if(p == 2) {
    if(is.na(axis)) {
      # if(!is.null(expr_dim))
      #  stop("Dimensions should be NULL")
      if(!all(expr_dim == c(1,1)))
        stop("Dimensions should be c(1,1)")
      return(list(t, list(SOC(t, vec(x)))))
    } else
      return(list(t, list(SOC(vec(t), x, axis))))
  }

  # We need an absolute value constraint for the symmetric convex branches (p > 1)
  constraints <- list()
  if(p > 1) {
    # TODO: Express this more naturally (recursively) in terms of the other atoms
    abs_expr <- abs(x)
    canon <- EliminatePwl.abs_canon(abs_expr, abs_expr@args)
    x <- canon[[1]]
    abs_constraints <- canon[[2]]
    constraints <- c(constraints, abs_constraints)
  }

  # Now, we take care of the remaining convex and concave branches to create the
  # rational powers. We need a new variable, r, and the constraint sum(r) == t
  # r <- Variable(dim(x))
  r <- new("Variable", dim = dim(x))
  constraints <- c(constraints, list(sum(r) == t))

  # TODO: No need to run gm_constr to form the tree each time.
  # We only need to form the tree once.
  promoted_t <- Constant(matrix(1, nrow = nrow(x), ncol = ncol(x))) %*% t
  p <- gmp::as.bigq(p)
  if(p < 0)
    constraints <- c(constraints, gm_constrs(promoted_t, list(x, r), c(-p/(1-p), 1/(1-p))))
  else if(p > 0 && p < 1)
    constraints <- c(constraints, gm_constrs(r, list(x, promoted_t), c(p, 1-p)))
  else if(p > 1)
    constraints <- c(constraints, gm_constrs(x, list(r, promoted_t), c(1/p, 1-1/p)))
  return(list(t, constraints))
}

#' 
#' Dcp2Cone canonicalizer for the power atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a power atom, where
#' the objective function consists of the variable t which is
#' of the dimension of the original vector from the power atom
#' and the constraints consists of geometric mean constraints.
Dcp2Cone.power_canon <- function(expr, args) {
  x <- args[[1]]
  p <- expr@p_rational
  w <- expr@w

  if(p == 1)
    return(list(x, list()))

  expr_dim <- dim(expr)
  if(is.null(expr_dim))
    ones <- Constant(1)
  else
    ones <- Constant(matrix(1, nrow = expr_dim[1], ncol = expr_dim[2]))
  if(p == 0)
    return(list(ones, list()))
  else {
    # t <- Variable(expr_dim)
    t <- new("Variable", dim = expr_dim)
    # TODO: gm_constrs requires each of its inputs to be a Variable; is this something that we want to change?
    if(p > 0 && p < 1)
      return(list(t, gm_constrs(t, list(x, ones), w)))
    else if(p > 1)
      return(list(t, gm_constrs(x, list(t, ones), w)))
    else if(p < 0)
      return(list(t, gm_constrs(ones, list(x, t), w)))
    else
      stop("This power is not yet supported")
  }
}

#' 
#' Dcp2Cone canonicalizer for the quadratic form atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a quadratic form atom,
#' where the objective function consists of the scaled objective function
#' from the quadratic over linear canonicalization and same with the
#' constraints.
Dcp2Cone.quad_form_canon <- function(expr, args) {
  # TODO: This doesn't work with parameters!
  decomp <- .decomp_quad(value(args[[2]]))
  scale <- decomp[[1]]
  M1 <- decomp[[2]]
  M2 <- decomp[[3]]
  M1_dim <- dim(M1)
  M2_dim <- dim(M2)
  
  # Special case where P == 0.
  if((is.null(M1_dim) || size(M1) == 0) && (is.null(M2_dim) && size(M2) == 0))
    return(list(Constant(0), list()))

  if(!is.null(M1_dim) && prod(M1_dim) > 0)
    expr <- sum_squares(Constant(t(M1)) %*% args[[1]])
  else if(!is.null(M2_dim) && prod(M2_dim) > 0) {
    scale <- -scale
    expr <- sum_squares(Constant(t(M2)) %*% args[[1]])
  }
  canon <- Dcp2Cone.quad_over_lin_canon(expr, expr@args)
  obj <- canon[[1]]
  constr <- canon[[2]]
  return(list(scale * obj, constr))
}

#' 
#' Dcp2Cone canonicalizer for the quadratic over linear term atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a quadratic over linear
#' term atom where the objective function consists of a one
#' dimensional variable t with SOC constraints.
Dcp2Cone.quad_over_lin_canon <- function(expr, args) {
  # quad_over_lin := sum_{ij} X^2_{ij} / y
  x <- args[[1]]
  y <- flatten(args[[2]])
  
  # Pre-condition: dim = c()
  t <- Variable(1)
  
  # (y+t, y-t, 2*x) must lie in the second-order cone, where y+t is the scalar part
  # of the second-order cone constraint
  # BUG: In Python, flatten produces single dimension (n,), but in R, we always treat 
  # these as column vectors with dimension (n,1), necessitating the use of VStack.
  # constraints <- list(SOC(t = y+t, X = HStack(y-t, 2*flatten(x)), axis = 2))
  constraints <- list(SOC(t = y+t, X = VStack(y-t, 2*flatten(x)), axis = 2))
  return(list(t, constraints))
}

#' 
#' Dcp2Cone canonicalizer for the relative entropy atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a relative entropy atom
#' where the objective function consists of the variable t
#' that is of the same dimension as the original expression
#' with specified constraints in the function.
Dcp2Cone.rel_entr_canon <- function(expr, args) {
  expr_dim <- dim(expr)
  x <- promote(args[[1]], expr_dim)
  y <- promote(args[[2]], expr_dim)
  t <- new("Variable", dim = expr_dim)
  constraints <- list(ExpCone(t, x, y))
  obj <- -t
  return(list(obj, constraints))
}

#' 
#' Dcp2Cone canonicalizer for the sigma max atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a sigma max atom
#' where the objective function consists of the variable t
#' that is of the same dimension as the original expression
#' with specified constraints in the function.
Dcp2Cone.sigma_max_canon <- function(expr, args) {
  A <- args[[1]]
  A_dim <- dim(A)
  n <- A_dim[1]
  m <- A_dim[2]
  expr_dim <- dim(expr)
  if(prod(expr_dim) != 1)
    stop("Invalid shape of expr in sigma_max canonicalization")
  
  t <- new("Variable", dim = expr_dim)
  tI_n <- sparseMatrix(seq(n), seq(n), x = rep(1, n)) * t
  tI_m <- sparseMatrix(seq(m), seq(m), x = rep(1, m)) * t
  X <- bmat(list(list(tI_n, A), 
                 list(t(A), tI_m)))
  constraints <- list(PSDConstraint(X))
  return(list(t, constraints))
}

#' 
#' Dcp2Cone canonicalizer for the support function atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a support function atom
#' where the objective function consists of the variable t
#' that is of the same dimension as the original expression
#' with specified constraints in the function.
Dcp2Cone.suppfunc_canon <- function(expr, args) {
  # The user-supplied argument to the support function:
  y <- flatten(args[[1]])
  parent <- expr@.parent
  
  # This Defines the set "X" associated with this support function:
  conic <- conic_repr_of_set(parent)
  A <- conic[[1]]
  b <- conic[[2]]
  K_self <- conic[[3]]
  
  # The main part of the duality trick for representing the epigraph
  # of this support function:
  eta <- Variable(size(b))
  expr@.eta <- eta

  n <- ncol(A)
  n0 <- size(y)
  if(n > n0) {
    # The description of the set "X" used in this support
    # function included n - n0 > 0 auxiliary variables.
    # We can pretend these variables were user-defined
    # by appending a suitable number of zeros to y.
    y_lift <- HStack(list(y, rep(0, n - n0)))
  } else
    y_lift <- y
  local_cons <- list(t(A) %*% eta + y_lift == 0)
  
  # Now, the conic constraints on eta.
  #   nonneg, exp, soc, psd
  nonnegsel <- K_sels$nonneg
  if(size(nonnegsel) > 0) {
    temp_expr <- eta[nonnegsel]
    local_cons <- c(local_cons, list(temp_expr >= 0))
  }
  
  socsels <- K_sels$soc
  socsels_len <- length(socsels)
  for(socsel in socsels) {
    tempsca <- eta[socsel[1]]
    tempvec <- eta[socsel[seq(2, socsels_len)]]
    soccon <- SOC(tempsca, tempvec)
    local_cons <- c(local_cons, list(soccon))
  }
  
  psdsels <- K_sels$psd
  for(psdsel in psdsels) {
    curmat <- scs_psdvec_to_psdmat(eta, psdsel)
    local_cons <- c(local_cons, list(curmat %>>% 0))
  }
  
  expsel <- K_sels$exp
  if(size(expsel) > 0) {
    matexpsel <- matrix(expsel, ncol = 3, byrow = TRUE)
    curr_u <- eta[matexpsel[,1]]
    curr_v <- eta[matexpsel[,2]]
    curr_w <- eta[matexpsel[,3]]
    # (curr_u, curr_v, curr_w) needs to belong to the dual
    # exponential cone, as used by the SCS solver. We map
    # this to a primal exponential cone as follows.
    ec <- ExpCone(-curr_v, -curr_u, base::exp(1) * curr_w)
    local_cons.append(ec)
  }
  
  epigraph <- b %*% eta
  return(list(epigraph, local_cons))
}

#' 
#' Dcp2Cone canonicalizer for the trinv atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a trinv atom
#' where the objective function consists of the variable t
#' that is of the same dimension as the original expression
#' with specified constraints in the function.
Dcp2Cone.tr_inv_canon <- function(expr, args) {
  X <- args[[1]]
  n <- nrow(X)
  su <- NULL
  
  constraints <- list()
  for(i in seq_len(n)) {
    ei <- matrix(0, nrow = n, ncol = 1)
    ei[i] <- 1.0
    ui <- Variable(1, 1)
    R <- bmat(list(list(X, ei),
                   list(t(ei), ui)))
    constraints <- c(constraints, list(R %>>% 0))
    if(is.null(su))
      su <- ui
    else
      su <- su + ui
  }
  return(list(su, constraints))
}

#' 
#' Dcp2Cone canonicalizer for the von Neumann entropy atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a von Neumann entropy atom
#' where the objective function consists of the variable t
#' that is of the same dimension as the original expression
#' with specified constraints in the function.
Dcp2Cone.von_neumann_entr_canon <- function(expr, args) {
  N <- args[[1]]
  if(!is_real(N))
    stop("N must be real")
  
  n <- nrow(N)
  x <- Variable(n)
  t <- Variable()
  
  # START code that applies to all spectral functions
  constrs <- list()
  if(n > 1) {
    for(r in seq(2, n)) {
      # lambda_sum_largest(N, r) <= sum(x[1:(r-1)])
      expr_r <- lambda_sum_largest(N, r)
      canon <- Dcp2Cone.lambda_sum_largest_canon(expr_r, expr_r@args)
      epi <- canon[[1]]
      cons <- canon[[2]]
      constrs <- c(constrs, cons)
      con <- NonPosConstraint(epi - sum(x[1:(r-1)]))
      constrs <- c(constrs, list(con))
    }
  }
  
  # trace(N) \leq sum(x)
  con <- (matrix_trace(N) == sum(x))
  constrs <- c(constrs, list(con))
  
  # trace(N) == sum(x)
  con <- ZeroConstraint(matrix_trace(N) - sum(x))
  constrs <- c(constrs, list(con))
  
  # x[1:(n-1)] >= x[2:n]
  #   x[1] >= x[2],  x[2] >= x[3], ...
  if(n > 1) {
    con <- NonPosConstraint(x[2:n] - x[1:(n - 1)])
    constrs <- c(constrs, list(con))
  }
  
  # END code that applies to all spectral functions
  
  # sum(entr(x)) >= t
  canon <- Dcp2Cone.entr_canon(x, list(x))
  hypos <- canon[[1]]
  entr_cons <- canon[[2]]
  constrs <- c(constrs, entr_cons)
  con <- NonPosConstraint(t - sum(hypos))
  constrs <- c(constrs, list(con))
  
  return(list(t, constrs))
}

#' 
#' Dcp2Cone canonicalizer for the xexp atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A cone program constructed from a xexp atom
#' where the objective function consists of the variable t
#' that is of the same dimension as the original expression
#' with specified constraints in the function.
Dcp2Cone.xexp_canon <- function(expr, args) {
  x <- args[[1]]
  u <- new("Variable", dim = dim(expr), nonneg = TRUE)
  t <- new("Variable", dim = dim(expr), nonneg = TRUE)
  power_expr <- Power(x, 2)
  canon <- Dcp2Cone.power_canon(power_expr, power_expr@args)
  power_obj <- canon[[1]]
  constraints <- canon[[2]]
  
  constraints <- c(constraints, list(ExpCone(u, x, t), u >= power_obj, x >= 0))
  return(list(t, constraints))
}

# TODO: Remove pwl canonicalize methods and use EliminatePwl reduction instead.
Dcp2Cone.CANON_METHODS <- list(CumMax = EliminatePwl.CANON_METHODS$CumMax,
                               CumSum = EliminatePwl.CANON_METHODS$CumSum,
                               GeoMean = Dcp2Cone.geo_mean_canon,
                               LambdaMax = Dcp2Cone.lambda_max_canon,
                               LambdaSumLargest = Dcp2Cone.lambda_sum_largest_canon,
                               LogDet = Dcp2Cone.log_det_canon,
                               LogSumExp = Dcp2Cone.log_sum_exp_canon,
                               MatrixFrac = Dcp2Cone.matrix_frac_canon,
                               MaxEntries = EliminatePwl.CANON_METHODS$MaxEntries,
                               MinEntries = EliminatePwl.CANON_METHODS$MinEntries,
                               Norm1 = EliminatePwl.CANON_METHODS$Norm1,
                               NormNuc = Dcp2Cone.normNuc_canon,
                               NormInf = EliminatePwl.CANON_METHODS$NormInf,
                               Pnorm = Dcp2Cone.pnorm_canon,
                               QuadForm = Dcp2Cone.quad_form_canon,
                               QuadOverLin = Dcp2Cone.quad_over_lin_canon,
                               SigmaMax = Dcp2Cone.sigma_max_canon,
                               SumLargest = EliminatePwl.CANON_METHODS$SumLargest,
                               Abs = EliminatePwl.CANON_METHODS$Abs,
                               Entr = Dcp2Cone.entr_canon,
                               Exp = Dcp2Cone.exp_canon,
                               Huber = Dcp2Cone.huber_canon,
                               KLDiv = Dcp2Cone.kl_div_canon,
                               Log = Dcp2Cone.log_canon,
                               Log1p = Dcp2Cone.log1p_canon,
                               Logistic = Dcp2Cone.logistic_canon,
                               MaxElemwise = EliminatePwl.CANON_METHODS$MaxElemwise,
                               MinElemwise = EliminatePwl.CANON_METHODS$MinElemwise,
                               Power = Dcp2Cone.power_canon,
                               Indicator = Dcp2Cone.indicator_canon,
                               SpecialIndex = special_index_canon)
