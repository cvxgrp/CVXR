# Reduction utility functions.
lower_inequality <- function(inequality) {
  lhs <- inequality@args[[1]]
  rhs <- inequality@args[[2]]
  NonPosConstraint(lhs - rhs, id = inequality@id)
}

lower_equality <- function(equality) {
  lhs <- equality@args[[1]]
  rhs <- equality@args[[2]]
  ZeroConstraint(lhs - rhs, id = equality@id)
}

special_index_canon <- function(expr, args) {
  select_mat <- expr@.select_mat
  final_dim <- dim(expr@.select_mat)
  if(is.null(final_dim))
    final_dim <- c(length(select_mat), 1)
  # select_vec <- matrix(select_mat, nrow = prod(final_dim))
  select_vec <- as.vector(select_mat)

  # Select the chosen entries from expr.
  arg <- args[[1]]
  arg_size <- size(arg)
  # identity <- diag(size(arg))
  identity <- sparseMatrix(i = 1:arg_size, j = 1:arg_size, x = rep(1, arg_size))

  v <- vec(arg)
  idmat <- matrix(identity[select_vec,], ncol = arg_size)
  if(is_scalar(v) || is_scalar(as.Constant(idmat)))
    lowered <- Reshape(idmat * v, final_dim)
  else
    lowered <- Reshape(idmat %*% v, final_dim)
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

setClassUnion("ProblemORNULL", c("Problem", "NULL"))

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
setClass("Reduction", representation(problem = "ProblemORNULL", .emitted_problem = "ANY", .retrieval_data = "ANY"),
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
  object@.emitted_problem <- tmp[[1]]
  object@.retrieval_data <- tmp[[2]]
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

#'
#' The Canonicalization class.
#'
#' This class represents a canonicalization reduction.
#'
#' @rdname Canonicalization-class
.Canonicalization <- setClass("Canonicalization", representation(canon_methods = "list"), prototype(canon_methods = list()), contains = "Reduction")

Canonicalization <- function(problem, canon_methods) { .Canonicalization(problem = problem, canon_methods = canon_methods) }

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
    vid <- as.character(inverse_data@cons_id_map[[orig_id]])
    if(vid %in% names(solution@dual_vars))
      dvars[[orig_id]] <- solution@dual_vars[[vid]]
  }
  return(Solution(solution@status, solution@opt_val, pvars, dvars, solution@attr))
})

setMethod("canonicalize_tree", "Canonicalization", function(object, expr) {
  # TODO: Don't copy affine expressions?
  if(class(expr) == "PartialProblem") {
    canon <- canonicalize_tree(object, expr@args[[1]]@objective@expr)
    canon_expr <- canon[[1]]
    constrs <- canon[[2]]
    for(constr in expr@args[[1]]@constraints) {
      canon <- canonicalize_tree(object, constr)
      canon_constr <- canon[[1]]
      aux_constr <- canon[[2]]
      constrs <- c(constrs, list(canon_constr), aux_constr)
    }
  } else {
    canon_args <- list()
    constrs <- list()
    for(arg in expr@args) {
      canon <- canonicalize_tree(object, arg)
      canon_arg <- canon[[1]]
      c <- canon[[2]]
      canon_args <- c(canon_args, list(canon_arg))
      constrs <- c(constrs, c)
    }
    canon <- canonicalize_expr(object, expr, canon_args)
    canon_expr <- canon[[1]]
    c <- canon[[2]]
    constrs <- c(constrs, c)
  }
  return(list(canon_expr, constrs))
})

setMethod("canonicalize_expr", "Canonicalization", function(object, expr, args) {
  if(is(expr, "Expression") && is_constant(expr)) {
    return(list(expr, list()))
  } else if(class(expr) %in% names(object@canon_methods))
    return(object@canon_methods[[class(expr)]](expr, args))
  else
    return(list(copy(expr, args), list()))
})

#'
#' The Chain class.
#'
#' This class represents a reduction that replaces symbolic parameters with
#' their constraint values.
#'
#' @rdname Chain-class
.Chain <- setClass("Chain", representation(reductions = "list"), prototype(reductions = list()), contains = "Reduction")
Chain <- function(problem = NULL, reductions = list()) { .Chain(problem = problem, reductions = reductions) }

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
    inverse_data <- c(inverse_data, list(inv))
  }
  return(list(problem, inverse_data))
})

setMethod("invert", signature(object = "Chain", solution = "SolutionORList", inverse_data = "list"), function(object, solution, inverse_data) {
  m <- min(length(object@reductions), length(inverse_data))
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
    any(sapply(variables, function(v) { !is.null(v@attributes[[attr]]) && v@attributes[[attr]] }))
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
CvxAttr2Constr <- setClass("CvxAttr2Constr", contains = "Reduction")

setMethod("accepts", signature(object = "CvxAttr2Constr", problem = "Problem"), function(object, problem) { TRUE })

setMethod("perform", signature(object = "CvxAttr2Constr", problem = "Problem"), function(object, problem) {
  if(length(convex_attributes(variables(problem))) == 0)
    return(list(problem, list()))

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

      if(length(symmetric_attributes(list(var))) > 0) {
        n <- nrow(var)
        new_dim <- c(floor(n*(n+1)/2), 1)
        # upper_tri <- do.call(Variable, c(list(new_dim), new_attr))
        upper_tri <- do.call(.Variable, c(list(dim = new_dim), new_attr))
        id2new_var[[vid]] <- upper_tri
        fill_coeff <- Constant(upper_tri_to_full(n))
        full_mat <- fill_coeff %*% upper_tri
        obj <- reshape_expr(full_mat, c(n, n))
      } else if(!is.null(var@attributes$diag)) {
        # diag_var <- do.call(Variable, c(list(nrow(var)), new_attr))
        diag_var <- do.call(.Variable, c(list(dim = c(nrow(var), 1)), new_attr))
        id2new_var[[vid]] <- diag_var
        obj <- diag(diag_var)
      } else if(new_var) {
        # obj <- do.call(Variable, c(list(dim(var)), new_attr))
        obj <- do.call(.Variable, c(list(dim = dim(var)), new_attr))
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
  obj <- tree_copy(problem@objective, id_objects = id2new_obj)
  cons_id_map <- list()
  for(cons in problem@constraints) {
    constr <- c(constr, tree_copy(cons, id_objects = id2new_obj))
    cons_id_map[[as.character(cons@id)]] <- constr[[length(constr)]]@id
  }
  inverse_data <- list(id2new_var, id2old_var, cons_id_map)
  return(list(Problem(obj, constr), inverse_data))
})

setMethod("invert", signature(object = "CvxAttr2Constr", solution = "Solution", inverse_data = "list"), function(object, solution, inverse_data) {
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
    vid <- as.character(cons_id_map[[orig_id]])
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
EvalParams <- setClass("EvalParams", contains = "Reduction")

EvalParams.replace_params_with_consts <- function(expr) {
  if(is.list(expr))
    return(lapply(expr, function(elem) { EvalParams.replace_params_with_consts(elem) }))
  else if(length(parameters(expr)) == 0)
    return(expr)
  else if(is(expr, "Parameter")) {
    if(any(is.na(value(expr))))
      stop("Problem contains unspecified parameters")
    return(Constant(value(expr)))
  } else {
    new_args <- list()
    for(arg in expr@args)
      new_args <- c(new_args, EvalParams.replace_params_with_consts(arg))
    return(copy(expr, new_args))
  }
}

setMethod("accepts", signature(object = "EvalParams", problem = "Problem"), function(object, problem) { TRUE })
setMethod("perform", signature(object = "EvalParams", problem = "Problem"), function(object, problem) {
  # Do not instantiate a new objective if it does not contain parameters.
  if(length(parameters(problem@objective)) > 0) {
    obj_expr <- EvalParams.replace_params_with_consts(problem@objective@expr)
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
      args <- c(args, EvalParams.replace_params_with_consts(arg))

    # Do not instantiate a new constraint object if it did not contain parameters.
    id_match <- mapply(function(new, old) { new@id == old@id }, args, c@args)
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
  return(list(Problem(objective, constraints), list()))
})

setMethod("invert", signature(object = "EvalParams", solution = "Solution", inverse_data = "list"), function(object, solution, inverse_data) { solution })

#'
#' The FlipObjective class.
#'
#' This class represents a reduction that flips a minimization objective to a
#' maximization and vice versa.
#'
#' @rdname FlipObjective-class
FlipObjective <- setClass("FlipObjective", contains = "Reduction")

setMethod("accepts", signature(object = "FlipObjective", problem = "Problem"), function(object, problem) { TRUE })
setMethod("perform", signature(object = "FlipObjective", problem = "Problem"), function(object, problem) {
  if(class(problem@objective) == "Maximize")
    objective <- Minimize
  else
    objective <- Maximize
  problem <- Problem(objective(-problem@objective@expr), problem@constraints)
  return(list(problem, list()))
})

setMethod("invert", signature(object = "FlipObjective", solution = "Solution", inverse_data = "list"), function(object, solution, inverse_data) {
  if(!is.null(solution@opt_val))
    solution@opt_val <- -solution@opt_val
  return(solution)
})

#'
#' The MatrixStuffing class.
#'
#' @rdname MatrixStuffing-class
MatrixStuffing <- setClass("MatrixStuffing", contains = "Reduction")

setMethod("perform", signature(object = "MatrixStuffing", problem = "Problem"), function(object, problem) {
  inverse_data <- InverseData(problem)

  # Form the constraints
  extractor <- CoeffExtractor(inverse_data)
  stuffed <- stuffed_objective(object, problem, extractor)
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
      con <- SOC(con@args[[1]], t(con@args[[2]]), axis = 2, id = con@id)
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
      arg_size <- size(arg)
      if(arg_size == 1) {
        A <- matrix(Afull[offset + 1,], nrow = 1)
        b <- bfull[offset + 1]
      } else {
        A <- Afull[(offset + 1):(offset + arg_size),]
        b <- bfull[(offset + 1):(offset + arg_size)]
      }
      arg_list <- c(arg_list, reshape_expr(A %*% new_var + b, dim(arg)))
      offset <- offset + arg_size
    }
    new_cons <- c(new_cons, copy(con, arg_list))
    inverse_data@cons_id_map[[as.character(con@id)]] <- new_cons[[length(new_cons)]]@id
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
    primal_vars[[var_id]] <- matrix(x_opt[(offset+1):(offset+size)], nrow = var_dim[1], ncol = var_dim[2])
  }

  # Remap dual variables if dual exists (problem is convex).
  if(length(solution@dual_vars) > 0) {
    for(old_con in names(con_map)) {
      new_con <- as.character(con_map[[old_con]])
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

setMethod("stuffed_objective", signature(object = "MatrixStuffing", problem = "Problem", extractor = "CoeffExtractor"), function(object, problem, extractor) {
  stop("Unimplemented")
})

extract_mip_idx <- function(variables) {
  # Coalesces bool, int indices for variables.
  # The indexing scheme assumes that the variables will be coalesced into a single
  # one-dimensional variable with each variable being reshaped in Fortran order.

  ravel_multi_index <- function(multi_index, x, vert_offset) {
    # Ravel a multi-index and add a vertical offset to it
    ravel_idx <- array(FALSE, dim(x))
    ravel_idx[multi_index] <- TRUE
    ravel_idx <- which(ravel_idx, arr.ind = FALSE)
    return(sapply(ravel_idx, function(idx) { vert_offset + idx }))
  }

  boolean_idx <- c()
  integer_idx <- c()
  vert_offset <- 0
  for(x in variables) {
    if(nrow(x@boolean_idx) > 0) {
      multi_index <- x@boolean_idx
      boolean_idx <- c(boolean_idx, ravel_multi_index(multi_index, x, vert_offset))
    }
    if(nrow(x@integer_idx) > 0) {
      multi_index <- x@integer_idx
      integer_idx <- c(integer_idx, ravel_multi_index(multi_index, x, vert_offset))
    }
    vert_offset <- vert_offset + size(x)
  }

  if(is.null(boolean_idx))
    boolean_idx <- matrix(0, nrow = 0, ncol = 1)
  else
    boolean_idx <- as.matrix(boolean_idx)
  if(is.null(integer_idx))
    integer_idx <- matrix(0, nrow = 0, ncol = 1)
  else
    integer_idx <- as.matrix(integer_idx)
  return(list(boolean_idx, integer_idx))
}
