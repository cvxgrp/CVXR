# Reduction utility functions.
lower_inequality <- function(inequality) {
  lhs <- inequality@args[[1]]
  rhs <- inequality@args[[2]]
  NonPosConstraint(lhs - rhs, constr_id = inequality@constr_id)
}

lower_equality <- function(equality) {
  lhs <- equality@args[[1]]
  rhs <- equality@args[[2]]
  ZeroConstraint(lhs - rhs, constr_id = equality@constr_id)
}

special_index_canon <- function(expr, args) {
  select_mat <- expr@.select_mat
  final_dim <- dim(expr@.select_mat)
  select_vec <- matrix(select_mat, nrow = size(select_mat))
  
  # Select the chosen entries from expr.
  arg <- args[[1]]
  identity <- diag(size(arg))
  lowered <- Reshape(identity[select_vec]*vec(arg), final_dim)
  list(lowered, list())
}

are_args_affine <- function(constraints) {
  all(sapply(constraints, function(constr) {
    all(sapply(constr@args, function(arg) { is_affine(arg) }))
  }))
}

# Factory function for infeasible or unbounded solutions.
failure_solution <- function(status) {
  if(status == INFEASIBLE)
    opt_val = Inf
  else if(status == UNBOUNDED)
    opt_val = -Inf
  else
    opt_val = NA_real_
  return(Solution(status, opt_val, list(), list(), list()))
}

#'
#' The InverseData class.
#'
#' This class represents the data encoding an optimization problem.
#'
#' @rdname InverseData-class
.InverseData <- setClass("InverseData", representation(problem = "Problem", id_map = "list", var_offsets = "list", x_length = "numeric", var_dims = "list",
                                                       id2var = "list", real2imag = "list", id2cons = "list", cons_id_map = "list"),
                                        prototype(id_map = list(), var_offsets = list(), x_length = NA_real_, var_dims = list(), id2var = list(),
                                                  real2imag = list(), id2cons = list(), cons_id_map = list()))

InverseData <- function(problem) { .InverseData(problem = problem) }

setMethod("initialize", "InverseData", function(.Object, ..., problem, id_map = list(), var_offsets = list(), x_length = NA_real_, var_dims = list(), id2var = list(), real2imag = list(), id2cons = list(), cons_id_map = list()) {
  # Basic variable offset information
  varis <- variables(problem)
  varoffs <- get_var_offsets(.Object, varis)
  .Object@id_map <- varoffs$id_map
  .Object@var_offsets <- varoffs$var_offsets
  .Object@x_length <- varoffs$x_length
  .Object@var_dims <- varoffs$var_dims

  # Map of variable id to variable
  .Object@id2var <- setNames(varis, sapply(varis, function(var) { as.character(id(var)) }))

  # Map of real to imaginary parts of complex variables
  var_comp <- varis[sapply(varis, function(var) { is_complex(var) })]
  .Object@real2imag <- setNames(var_comp, sapply(var_comp, function(var) { as.character(id(var)) }))
  constrs <- constraints(problem)
  constr_comp <- constrs[sapply(constrs, function(cons) { is_complex(cons) })]
  constr_dict <- setNames(constr_comp, sapply(constr_comp, function(cons) { as.character(id(cons)) }))
  .Object@real2imag <- update(.Object@real2imag, constr_dict)

  # Map of constraint id to constraint
  .Object@id2cons <- setNames(constrs, sapply(constrs, function(cons) { id(cons) }))
  .Object@cons_id_map <- list()
  return(.Object)
})

setMethod("get_var_offsets", signature(object = "InverseData", variables = "list"), function(object, variables) {
  var_dims <- list()
  var_offsets <- list()
  id_map <- list()
  vert_offset <- 0
  for(x in variables) {
    var_dims[[as.character(id(x))]] <- dim(x)   # TODO: Redefine Variable class to include dim parameter
    var_offsets[[as.character(id(x))]] <- vert_offset
    id_map[[as.character(id(x))]] <- list(vert_offset, size(x))
    vert_offset <- vert_offset + size(x)
  }
  return(list(id_map = id_map, var_offsets = var_offsets, x_length = vert_offset, var_dims = var_dims))
})

setClassUnion("InverseDataORList", c("InverseData", "list"))

#'
#' The Reduction class.
#'
#' This virtual class represents a reduction, an actor that transforms a problem
#' into an equivalent problem. By equivalent, we mean that there exists a mapping
#' between solutions of either problem: if we reduce a problem \eqn{A} to another
#' problem \eqn{B} and then proceed to find a solution to \eqn{B}, we can convert
#' it to a solution of \eqn{A} with at most a moderate amount of effort.
#'
#' Every reduction supports three methods: accepts, perform, and invert. The accepts
#' method of a particular reduction codifies the types of problems that it is applicable
#' to, the perform method takes a problem and reduces it to a (new) equivalent form,
#' and the invert method maps solutions from reduced-to problems to their problems
#' of provenance.
#'
#' @rdname Reduction-class
setClass("Reduction", representation(problem = "Problem", .emitted_problem = "Problem", .retrieval_data = "InverseDataORList"), 
                      prototype(problem = NULL, .emitted_problem = NULL, .retrieval_data = NULL), contains = "VIRTUAL")

#' @param object A \linkS4class{Reduction} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn Reduction States whether the reduction accepts a problem.
setMethod("accepts", signature(object = "Reduction", problem = "Problem"), function(object, problem) { stop("Unimplemented") })

#' @describeIn Reduction Reduces the owned problem to an equivalent problem.
setMethod("reduce", "Reduction", function(object) {
  if(!is.null(object@.emitted_problem))
    return(object@.emitted_problem)
  
  if(is.null(object@problem))
    stop("The reduction was constructed without a Problem")
  
  tmp <- perform(object, object@problem)
  object@.emitted_problem <- problem
  object@.retrieval_data <- retrieval_data
  return(list(object, object@.emitted_problem))
})

setMethod("retrieve", signature(object = "Reduction", solution = "Solution"), function(object, solution) {
  if(is.null(object@.retrieval_data))
    stop("reduce must be called before retrieve")
  return(invert(object, solution, object@.retrieval_data))
})

#' @describeIn Reduction Performs the reduction on a problem and returns an equivalent problem.
setMethod("perform", signature(object = "Reduction", problem = "Problem"), function(object, problem) { stop("Unimplemented") })

#' @param solution A \linkS4class{Solution} to a problem that generated the inverse data.
#' @param inverse_data The data encoding the original problem.
#' @describeIn Reduction Returns a solution to the original problem given the inverse data.
setMethod("invert", signature(object = "Reduction", solution = "Solution", inverse_data = "list"), function(object, solution, inverse_data) { stop("Unimplemented") })

replace_params_with_consts <- function(expr) {
  if(is.list(expr))
    return(lapply(expr, function(elem) { replace_params_with_consts(elem) }))
  else if(len(parameters(expr)) == 0)
    return(expr)
  else if(is(expr, "Parameter")) {
    if(is.na(value(expr)))
      stop("Problem contains unspecified parameters")
    return(Constant(value(expr)))
  } else {
    new_args <- list()
    for(arg in expr@args)
        51
    new_args <- c(new_args, replace_params_with_consts(arg))
    return(new_args)
  }
}

#'
#' The Canonicalization class.
#'
#' This class represents a canonicalization reduction.
#'
#' @rdname Canonicalization-class
.Canonicalization <- setClass("Canonicalization", representation(problem = "Problem", canon_methods = "list"), prototype(canon_methods = list()), contains = "Reduction")

Canonicalization <- function(problem, canon_methods) { .Canonicalization(problem = problem, canon_methods = canon_methods) }

setMethod("initialize", function(.Object, ..., problem, canon_methods) {
  .Object@canon_methods <- canon_methods
  callNextMethod(.Object, ..., problem = problem)
})

setMethod("perform", signature(object = "Canonicalization", problem = "Problem"), function(object, problem) {
  inverse_data <- InverseData(problem)

  canon <- canonicalize_tree(object, problem@objective)
  canon_objective <- canon[[1]]
  canon_constraints <- canon[[2]]

  for(constraint in problem@constraints) {
    # canon_constr is the constraint re-expressed in terms of its canonicalized arguments,
    # and aux_constr are the constraints generated while canonicalizing the arguments of the original constraint.
    canon <- canonicalize_tree(object, constraint)
    canon_constr <- canon[[1]]
    aux_constr <- canon[[2]]
    canon_constraints <- c(canon_constraints, aux_constr, list(canon_constr))
    inverse_data@cons_id_map[[as.character(constraint@id)]] <- canon_constr@id   # TODO: Check this updates like dict().update in Python
  }
  new_problem <- Problem(canon_objective, canon_constraints)
  return(list(new_problem, inverse_data))
})

setMethod("invert", signature(object = "Canonicalization", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  pvars <- list()
  for(vid in names(inverse_data@id_map)) {
    if(vid %in% names(solution@primal_vars))
      pvars[[as.character(vid)]] <- solution@primal_vars[[vid]]
  }

  dvars <- list()
  for(orig_id in names(inverse_data@cons_id_map)) {
    vid <- inverse_data@cons_id_map[[orig_id]]
    if(as.character(vid) %in% names(solution@dual_vars))
      dvars[[as.character(orig_id)]] <- solution@dual_vars[[as.character(vid)]]
  }
  return(Solution(solution@status, solution@opt_val, pvars, dvars, solution@attr))
})

setMethod("canonicalize_tree", signature(object = "Canonicalization", expr = "Expression"), function(object, expr) {
  # TODO: Don't copy affine expressions?
  if(class(expr) == "PartialProblem") {
    canon <- canonicalize_tree(object, expr@args[[1]]@objective@expr)
    canon_expr <- canon[[1]]
    constrs <- canon[[2]]
    for(constr in expr@args[[1]]@constraints) {
      canon <- canonicalize_tree(constr)
      canon_constr <- canon[[1]]
      aux_constr <- canon[[2]]
      constrs <- c(constrs, list(canon_constr), aux_constr)
    }
  } else {
    canon_args <- list()
    constrs <- list()
    for(arg in expr@args) {
      canon <- canonicalize_tree(arg)
      canon_arg <- canon[[1]]
      c <- canon[[2]]
      canon_args <- c(canon_args, canon_arg)
      constrs <- c(constrs, c)
    }
  }
  return(list(canon_expr, constrs))
})

setMethod("canonicalize_expr", signature(object = "Canonicalization", expr = "Expression", args = "list"), function(object, expr, args) {
  if(is(expr, "Expression") && is_constant(expr)) {
    return(list(expr, list()))
  } else if(class(expr) %in% canon_methods(object))
    return(canon_methods(object)[[class(expr)]](expr, args))
  else
    return(list(copy(expr, args), list()))   # TODO: copy must be defined in Canonical virtual class
})

#'
#' The Chain class.
#'
#' This class represents a reduction that replaces symbolic parameters with
#' their constraint values.
#'
#' @rdname Chain-class
.Chain <- setClass("Chain", representation(problem = "ProblemORNull", reductions = "list"), prototype(problem = NULL, reductions = list()), contains = "Reduction")

Chain <- function(problem, reductions) { .Chain(problem = problem, reductions = reductions) }

setMethod("initialize", function(.Object, ..., problem = NULL, reductions = list()) {
  .Object@reductions <- reductions
  callNextMethod(.Object, ..., problem = problem)
})

setMethod("as.character", "Chain", function(x) { paste(sapply(x@reductions, as.character), collapse = ", ") })
setMethod("show", "Chain", function(object) { paste("Chain(reductions = (", as.character(object@reductions),"))") })
setMethod("accepts", signature(object = "Chain", problem = "Problem"), function(object, problem) {
  for(r in object@reductions) {
    if(!accepts(r, problem))
      return(FALSE)
    problem <- perform(r, problem)[[1]]
  }
  return(TRUE)
})

setMethod("perform", signature(object = "Chain", problem = "Problem"), function(object, problem) {
  inverse_data <- list()
  for(r in object@reductions) {
    res <- perform(r, problem)
    problem <- res[[1]]
    inv <- res[[2]]
    inverse_data <- c(inverse_data, inv)
  }
  return(list(problem, inverse_data))
})

setMethod("invert", signature(object = "Chain", solution = "Solution", inverse_data = "list"), function(object, solution, inverse_data) {
  m <- min(length(object@reductions), inverse_data)
  for(i in rev(1:m)) {
    r <- object@reductions[[i]]
    inv <- inverse_data[[i]]
    solution <- invert(r, solution, inv)
  }
  return(solution)
})

# Returns a list of the relevant attributes present among the variables.
attributes_present <- function(variables, attr_map) {
  attr_map[sapply(attr_map, function(attr) {
    any(sapply(variables, function(v) { !is.null(v@attributes[[attr]]) }))
  })]
}

# Convex attributes that generate constraints.
CONVEX_ATTRIBUTES <- c("nonneg", "nonpos", "symmetric", "diag", "PSD", "NSD")
convex_attributes <- function(variables) {
  # Returns a list of the (constraint-generating) convex attributes present among the variables.
  attributes_present(variables, CONVEX_ATTRIBUTES)
}

# Attributes related to symmetry.
SYMMETRIC_ATTRIBUTES <- c("symmetric", "PSD", "NSD")
symmetric_attributes <- function(variables) {
  # Returns a list of the (constraint-generating) symmetric attributes present among the variables.
  attributes_present(variables, SYMMETRIC_ATTRIBUTES)
}

#'
#' The CvxAttr2Constr class.
#'
#' This class represents a reduction that expands convex variable attributes into constraints.
#'
#' @rdname CvxAttr2Constr-class
setClass("CvxAttr2Constr", contains = "Reduction")

setMethod("accepts", signature(object = "CvxAttr2Constr", problem = "Problem"), function(object, problem) { TRUE })

setMethod("perform", signature(object = "CvxAttr2Constr", problem = "Problem"), function(object, problem) {
  if(length(convex_attributes(variables(problem))) == 0)
    return(list(problem), list())

  # For each unique variable, add constraints.
  id2new_var <- list()
  id2new_obj <- list()
  id2old_var <- list()
  constr <- list()
  for(var in variables(problem)) {
    vid <- as.character(id(var))
    if(!(vid %in% names(id2new_var))) {
      id2old_var[[vid]] <- var
      new_var <- FALSE
      new_attr <- var@attributes
      for(key in CONVEX_ATTRIBUTES) {
        if(!is.null(new_attr[[key]])) {
          new_var <- TRUE
          new_attr[[key]] <- FALSE
        }
      }

      if(symmetric_attributes(list(var))) {
        n <- nrow(var)
        new_dim <- c(floor(n*(n+1)/2), 1)
        upper_tri <- do.call(Variable, c(list(new_dim), new_attr))
        id2new_var[[vid]] <- upper_tri
        fill_coeff <- Constant(upper_tri_to_full(n))
        full_mat <- fill_coeff %*% upper_tri
        obj <- reshape(full_mat, c(n, n))
      } else if(!is.null(var@attributes$diag)) {
        diag_var <- do.call(Variable, c(list(nrow(var)), new_attr))
        id2new_var[[vid]] <- diag_var
        obj <- diag(diag_var)
      } else if(new_var) {
        obj <- do.call(Variable, c(list(dim(var)), new_attr))
        id2new_var[[vid]] <- obj
      } else {
        obj <- var
        id2new_var[[vid]] <- obj
      }

      id2new_obj[[vid]] <- obj
      if(is_nonneg(var) || is_pos(var))
        constr <- c(constr, obj >= 0)
      else if(is_nonpos(var) || is_neg(var))
        constr <- c(constr, obj <= 0)
      else if(is_psd(var))
        constr <- c(constr, obj %>>% 0)
      else if(!is.null(var@attributes$NSD))
        constr <- c(constr, obj %<<% 0)
    }
  }

  # Create new problem.
  obj <- tree_copy(problem@objective, id_objects = id2new_obj)   # TODO: Implement tree_copy in Canonical virtual class
  cons_id_map <- list()
  for(cons in problem@constraints) {
    constr <- c(constr, tree_copy(cons, id_objects = id2new_obj))
    cons_id_map[[as.character(id(cons))]] <- id(constr[length(constr)])
  }
  inverse_data <- list(id2new_var, id2old_var, cons_id_map)
  return(list(Problem(obj, constr), inverse_data))
})

setMethod("invert", signature(object = "CvxAttr2Constr", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  if(is.null(inverse_data) || length(inverse_data) == 0)
    return(solution)

  id2new_var <- inverse_data[[1]]
  id2old_var <- inverse_data[[2]]
  cons_id_map <- inverse_data[[3]]
  pvars <- list()
  for(id in names(id2old_var)) {
    var <- id2old_var[[id]]
    new_var <- id2new_var[[id]]

    # Need to map from constrained to symmetric variable.
    nvid <- as.character(id(new_var))
    if(nvid %in% names(solution@primal_vars)) {
      if(!is.null(var@attributes$diag))
        pvars[[id]] <- Diagonal(x = as.vector(solution@primal_vars[[nvid]]))
      else if(length(symmetric_attributes(list(var))) > 0) {
        n <- nrow(var)
        value <- matrix(0, nrow = n, ncol = n)   # Variable is symmetric
        idxs <- upper.tri(value, diag = TRUE)
        value[idxs] <- as.vector(solution@primal_vars[[nvid]])
        value <- value + t(value) - diag(diag(value))
        pvars[[id]] <- value
      } else
        pvars[[id]] <- project(var, solution@primal_vars[[nvid]])
    }
  }

  dvars <- list()
  for(orig_id in names(cons_id_map)) {
    vid <- cons_id_map[[orig_id]]
    if(vid %in% names(solution@dual_vars))
      dvars[[orig_id]] <- solution@dual_vars[[vid]]
  }
  return(Solution(solution@status, solution@opt_val, pvars, dvars, solution@attr))
})

#'
#' The EvalParams class.
#'
#' This class represents a reduction that replaces symbolic parameters with
#' their constaint values.
#'
#' @rdname EvalParams-class
setClass("EvalParams", contains = "Reduction")

setMethod("accepts", signature(object = "EvalParams", problem = "Problem"), function(object, problem) { TRUE })
setMethod("perform", signature(object = "EvalParams", problem = "Problem"), function(object, problem) {
  # Do not instantiate a new objective if it does not contain parameters.
  if(length(parameters(problem@objective)) > 0) {
    obj_expr <- replace_params_with_consts(problem@objective@expr)
    if(class(problem@objective) == "Maximize")
      objective <- Maximize(obj_expr)
    else
      objective <- Minimize(obj_expr)
  } else
    objective <- problem@objective

  constraints <- list()
  for(c in problem@constraints) {
    args <- list()
    for(arg in c@args)
      args <- c(args, replace_params_with_consts(arg))

    # Do not instantiate a new constraint object if it did not contain parameters.
    id_match <- mapply(function(new, old) { id(new) == id(old) }, args, c@args)
    if(all(unlist(id_match)))
      constraints <- c(constraints, c)
    else {   # Otherwise, create a copy of the constraint.
      data <- get_data(c)
      if(!is.na(data) && length(data) > 0)
        constraints <- c(constraints, do.call(class(c), c(args, data)))
      else
        constraints <- c(constraints, do.call(class(c), args))
    }
  }
  return(list(Problem(objective, constraints)), list())
})

setMethod("invert", signature(object = "EvalParams", solution = "Solution", inverse_data = "list"), function(object, solution, inverse_data) { solution })

#'
#' The FlipObjective class.
#'
#' This class represents a reduction that flips a minimization objective to a
#' maximization and vice versa.
#'
#' @rdname FlipObjective-class
setClass("FlipObjective", contains = "Reduction")

setMethod("accepts", signature(object = "FlipObjective", problem = "Problem"), function(object, problem) { TRUE })
setMethod("perform", signature(object = "FlipObjective", problem = "Problem"), function(object, problem) {
  is_maximize <- class(problem@objective) == "Maximize"
  if(class(problem@objective) == "Maximize")
    objective <- Maximize
  else
    objective <- Minimize
  problem <- Problem(objective(-problem@objective@expr), problem@constraints)
  return(list(problem, list()))
})

setMethod("invert", signature(object = "FlipObjective", solution = "Solution", inverse_data = "list"), function(object, solution, inverse_data) {
  if(!is.null(solution@opt_val))
    solution@opt_val <- -solution@opt_val
  return(solution)
})

extract_mip_idx <- function(variables) {
  # Coalesces bool, int indices for variables.
  # The indexing scheme assumes that the variables will be coalesced into a single
  # one-dimensional variable with each variable being reshaped in Fortran order.

  ravel_multi_index <- function(multi_index, x, vert_offset) {
    # Ravel a multi-index and add a vertical offset to it
    # TODO: I have no idea what the following Python code does
    # ravel_idx <- np.ravel_multi_index(multi_index, max(x.shape, (1,)), order = "F")
    return(sapply(ravel_idx, function(idx) { vert_offset + idx }))
  }

  boolean_idx <- c()
  integer_idx <- c()
  vert_offset <- 0
  for(x in variables) {
    if(!is.null(x@boolean_idx)) {
      multi_index <- x@boolean_idx
      boolean_idx <- c(boolean_idx, ravel_multi_index(multi_index, x, vert_offset))
    }
    if(!is.null(x@integer_idx)) {
      multi_index <- x@integer_idx
      integer_idx <- c(integer_idx, ravel_multi_index(multi_index, x, vert_offset))
    }
    vert_offset <- vert_offset + size(x)
  }
  return(list(boolean_idx, integer_idx))
}

setClass("MatrixStuffing", contains = "Reduction")

setMethod("perform", signature(object = "MatrixStuffing", problem = "Problem"), function(object, problem) {
  inverse_data <- InverseData(problem)

  # Form the constraints
  extractor <- CoeffExtractor(inverse_data)
  stuffed <- stuffed_objective(object, problem, inverse_data)
  new_obj <- stuffed[[1]]
  new_var <- stuffed[[2]]
  inverse_data@r <- stuffed[[3]]
  
  # Lower equality and inequality to ZeroConstraint and NonPosConstraint.
  cons <- list()
  for(con in problem@constraints) {
    if(is(con, "EqConstraint"))
      con <- lower_equality(con)
    else if(is(con, "IneqConstraint"))
      con <- lower_inequality(con)
    else if(is(con, "SOC") && con@axis == 1)
      con <- SOC(con@args[[1]], t(con@args[[2]]), axis = 2, constr_id = con@constr_id)
    cons <- c(cons, con)
  }
    
  # Batch expressions together, then split apart.
  expr_list <- flatten_list(lapply(cons, function(c) { c@args }))
  Abfull <- affine(extractor, expr_list)
  Afull <- Abfull[[1]]
  bfull <- Abfull[[2]]
  new_cons <- list()
  offset <- 0
  
  for(con in cons) {
    arg_list <- list()
    for(arg in con@args) {
      A <- Afull[(offset + 1):(offset + size(arg) + 1),]
      b <- bfull[(offset + 1):(offset + size(arg) + 1)]
      arg_list <- c(arg_list, reshape(A %*% new_var + b, dim(arg)))
      offset <- offset + size(arg)
    }
    new_cons <- c(new_cons, copy(con, arg_list))
    inverse_data@cons_id_map[[as.character(id(con))]] <- id(new_cons[[length(new_cons)]])
  }

  # Map of old constraint id to new constraint id.
  inverse_data@minimize <- class(problem@objective) == "Minimize"
  new_prob <- Problem(Minimize(new_obj), new_cons)
  return(list(new_prob, inverse_data))
})

setMethod("invert", signature(object = "MatrixStuffing", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  var_map <- inverse_data@var_offsets
  con_map <- inverse_data@cons_id_map

  # Flip sign of optimal value if maximization.
  opt_val <- solution@opt_val
  if(!(solution@status %in% ERROR) && !inverse_data@minimize)
    opt_val <- -solution@opt_val

  primal_vars <- list()
  dual_vars <- list()
  if(!(solution@status %in% SOLUTION_PRESENT))
    return(Solution(solution@status, opt_val, primal_vars, dual_vars, solution@attr))

  # Split vectorized variable into components.
  x_opt <- solution@primal_vars[[1]]
  for(var_id in names(var_map)) {
    offset <- var_map[[var_id]]
    var_dim <- inverse_data@var_dims[[var_id]]
    size <- prod(var_dim)
    primal_vars[[var_id]] <- matrix(x_opt[offset:(offset+size)], nrow = var_dim[1], ncol = var_dim[2])
  }

  # Remap dual variables if dual exists (problem is convex).
  if(length(solution@dual_vars) > 0) {
    for(old_con in names(con_map)) {
      new_con <- con_map[[old_con]]
      con_obj <- inverse_data@id2cons[[old_con]]
      obj_dim <- dim(con_obj)
      # TODO: Rationalize Exponential.
      if(length(obj_dim) == 0 || is(con_obj, "ExpCone") || is(con_obj, "SOC"))
        dual_vars[[old_con]] <- solution@dual_vars[[new_con]]
      else
        dual_vars[[old_con]] <- matrix(solution@dual_vars[[new_con]], nrow = obj_dim[1], ncol = obj_dim[2])
    }
  }

  # Add constant part
  if(inverse_data@minimize)
    opt_val <- opt_val + inverse_data@r
  else
    opt_val <- opt_val - inverse_data@r

  return(Solution(solution@status, opt_val, primal_vars, dual_vars, solution@attr))
})

setMethod("stuffed_objective", signature(object = "MatrixStuffing", problem = "Problem", inverse_data = "InverseData"), function(object, problem, inverse_data) {
  stop("Unimplemented")
})
