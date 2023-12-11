#==============================#
# Reduction utility functions
#==============================#
lower_ineq_to_nonpos <- function(inequality) {
  lhs <- inequality@args[[1]]
  rhs <- inequality@args[[2]]
  NonPosConstraint(lhs - rhs, constr_id = inequality@constr_id)
}

lower_ineq_to_nonneg <- function(inequality) {
  lhs <- inequality@args[[1]]
  rhs <- inequality@args[[2]]
  NonNegConstraint(rhs - lhs, constr_id = inequality@constr_id)
}

lower_equality <- function(equality) {
  lhs <- equality@args[[1]]
  rhs <- equality@args[[2]]
  ZeroConstraint(lhs - rhs, constr_id = equality@constr_id)
}

nonpos2nonneg <- function(nonpos) {
  NonNegConstraint(-nonpos@expr, constr_id = nonpos@constr_id)
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
  identity <- sparseMatrix(i = seq(arg_size), j = seq(arg_size), x = rep(1, arg_size))

  v <- vec(arg)
  idmat <- matrix(identity[select_vec,], ncol = arg_size)
  if(is_scalar(v) || is_scalar(as.Constant(idmat)))
    lowered <- Reshape(idmat * v, final_dim)
  else
    lowered <- Reshape(idmat %*% v, final_dim)
  list(lowered, list())
}

#'
#' Are the arguments affine?
#'
#' @param constraints A list of \linkS4class{Constraint} object.
#' @return All the affine arguments in given constraints.
are_args_affine <- function(constraints) {
  all(sapply(constraints, function(constr) { all(sapply(constr@args, is_affine)) }))
}


## #' REPLACED BELOW by modified version
## #' Organize the constraints into a list keyed by constraint names.
## #'
## #'@param constraints A list of \linkS4class{Constraint} objects
## #'@return A list keyed by constraint types where list[[cone_type]]
## #'   maps to a list of exactly those constraints that are of type
## #'   cone_type.
## group_constraints <- function(constraints) {
##   constr_map <- list()
##   for(c in constraints) {
##     if(class(c) %in% names(constr_map))
##       constr_map[[class(c)]] <- c(constr_map[[class(c)]], c)
##     else
##       constr_map[[class(c)]] <- list(c)
##   }
##   return(constr_map)
## }


#'
#' Organize the constraints into a list keyed by constraint names.
#'
#'@param constraints A list of \linkS4class{Constraint} objects
#'@return A list keyed by constraint types where list[[cone_type]] maps to a list of exactly those constraints that are of type cone_type.
group_constraints <- function(constraints) {
  ## The constr_types list below should match the map named used in ConeDims-class (file dcp2cone.R)
  constr_map <- list()
  constr_names <- character(0)
  for (constr in constraints) {
    cl <- class(constr)
    index <- match(cl, constr_names)
    if (is.na(index)) { ## class not yet appeared in list
      constr_map[[cl]] <- list(constr) ## add new named item to list
      constr_names <- c(constr_names, cl)  ## add name to list of names
    } else {
      constr_map[[cl]] <- c(constr_map[[cl]], list(constr))
    }
  }
  constr_map
}

#'
#' The ReducedMat class.
#'
#' This is a utility class for condensing the mapping from parameters to problem data.
#'
#' For maximum efficiency of representation and application, the mapping from
#' parameters to problem data must be condensed. It begins as a CSC sparse matrix
#' matrix_data, such that multiplying by a parameter vector gives the problem data.
#' The row index array and column pointer array are saved as problem_data_index,
#' and a CSR matrix reduced_mat that when multiplied by a parameter vector gives
#' the values array. The ReducedMat class caches the condensed representation
#' and provides a method for multiplying by a parameter vector.
#'
#' This class consolidates code from ParamConeProg and ParamQuadProg.
#'
#' @rdname ReducedMat-class
ReducedMat <- setClass("ReducedMat", representation(matrix_data = "numeric", var_len = "integer", quad_form = "logical", reduced_mat = "numeric", problem_data_index = "ListORNULL", mapping_nonzero = "numeric"),
                        prototype(quad_form = FALSE, reduced_mat = NA_real_, problem_data_index = NULL, mapping_nonzero = NA_integer_))

#' @param keep_zeros (Optional) If TRUE, store explicit zeros in A where parameters are affected.
#' @describeIn ReducedMat Cache computed attributes if not present.
setMethod("cache", "ReducedMat", function(object, keep_zeros = FALSE) {
  # Short circuit null case.
  if(is.na(object@matrix_data))
    return(object)

  if(is.na(object@reduced_mat)) {
    # Form a reduced representation of the mapping, for faster application of parameters.
    if(!is.null(dim(object@matrix_data)) && prod(dim(object@matrix_data)) != 0) {
      tmp <- canonInterface.reduce_problem_data_tensor(object@matrix_data, object@var_len, object@quad_form)
      reduced_mat <- tmp[[1]]
      indices <- tmp[[2]]
      indptr <- tmp[[3]]
      shape <- tmp[[4]]
      object@reduced_mat <- reduced_mat
      object@problem_data_index <- list(indices, indptr, shape)
    } else {
      object@reduced_mat <- object@matrix_data
      object@problem_data_index <- NULL
    }
  }

  if(keep_zeros && is.na(object@mapping_nonzero)) {
    object@mapping_nonzero <- canonInterface.A_mapping_nonzero_rows(object@matrix_data, object@var_len)
  }
  return(object)
})

#' @param param_vec Flattened parameter vector
#' @param with_offset (Optional) A logical value indicating whether to return offset. Defaults to TRUE.
#' @describeIn ReducedMat Wraps get_matrix_from_tensor in canonInterface
setMethod("get_matrix_from_tensor", "ReducedMat", function(object, param_vec, with_offset = TRUE) {
  canonInterface.get_matrix_from_tensor(object@reduced_mat, param_vec, object@var_len, nonzero_rows = object@mapping_nonzero, with_offset = with_offset, problem_data_index = object@problem_data_index)
})

setClassUnion("ReducedMatORNULL", c("ReducedMat", "NULL"))

# Factory function for infeasible or unbounded solutions.
failure_solution <- function(status, attr = NULL) {
  if(status %in% c(INFEASIBLE, INFEASIBLE_INACCURATE))
    opt_val <- Inf
  else if(status %in% c(UNBOUNDED, UNBOUNDED_INACCURATE))
    opt_val <- -Inf
  else
    opt_val <- NA_real_

  if(is.null(attr))
    attr <- list()
  if(status == INFEASIBLE_OR_UNBOUNDED)
    attr$message <- INF_OR_UNB_MESSAGE

  return(Solution(status, opt_val, list(), list(), attr))
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
#' A reduction that is instantiated with a non-NULL problem offers two key methods:
#' reduce and retrieve. The reduce method converts the problem the reduction
#' was instantiated with to an equivalent problem. The retrieve method takes
#' as an argument a Solution for the equivalent problem and returns a Solution
#' for the problem owned by the reduction.
#'
#' Every reduction supports three low-level methods: accepts, perform, and invert.
#' The accepts method of a particular reduction specifies the types of problems
#' that it is applicable to, the perform method takes a problem and reduces it to
#' an equivalent form, and the invert method maps solutions from reduced-to problems
#' to their problems of provenance.
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
    return(list(object, object@.emitted_problem))

  if(is.null(object@problem))
    stop("The reduction was constructed without a Problem")

  tmp <- perform(object, object@problem)
  object <- tmp[[1]]
  object@.emitted_problem <- tmp[[2]]
  object@.retrieval_data <- tmp[[3]]
  return(list(object, object@.emitted_problem))
})

#' @param object A \linkS4class{Reduction} object.
#' @param solution A \linkS4class{Solution} object.
#' @describeIn Reduction Retrieves a solution to the owned problem.
setMethod("retrieve", signature(object = "Reduction", solution = "Solution"), function(object, solution) {
  if(is.null(object@.retrieval_data))
    stop("reduce must be called before retrieve")
  return(invert(object, solution, object@.retrieval_data))
})

#' @param object A \linkS4class{Reduction} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn Reduction Performs the reduction on a problem and returns an equivalent problem.
setMethod("perform", signature(object = "Reduction", problem = "Problem"), function(object, problem) { stop("Unimplemented") })

#' @param object A \linkS4class{Reduction} object.
#' @param solution A \linkS4class{Solution} to a problem that generated the inverse data.
#' @param inverse_data The data encoding the original problem.
#' @describeIn Reduction Returns a solution to the original problem given the inverse data.
setMethod("invert", signature(object = "Reduction", solution = "Solution", inverse_data = "list"), function(object, solution, inverse_data) { stop("Unimplemented") })

#'
#' The Canonicalization class.
#'
#' This class represents a canonicalization reduction.
#'
#' This reduction recursively canonicalizes every expression tree in a problem,
#' visiting each node. At every node, this reduction first canonicalizes its
#' arguments; it then canonicalizes the node, using the canonicalized arguments.
#'
#' The attribute canon_methods is a list mapping node types to functions that
#' canonicalize them; the signature of these canonicalizing functions must be
#' \code{canon_func(expr, canon_args) --> (new_expr, constraints) }
#' where expr is the Expression (node) to canonicalize, canon_args is a list of
#' the canonicalized arguments of this expression, new_expr is a canonicalized
#' expression, and constraints is a list of constraints introduced while
#' canonicalizing expr.
#'
#' @rdname Canonicalization-class
.Canonicalization <- setClass("Canonicalization", slots = list(canon_methods = "list"),
                              prototype = list(canon_methods = list()), contains = "Reduction")

Canonicalization <- function(problem, canon_methods) { .Canonicalization(problem = problem, canon_methods = canon_methods) }

#' @param object A \linkS4class{Canonicalization} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn Canonicalization Recursively canonicalize the objective and every constraint.
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
    inverse_data@cons_id_map[[as.character(id(constraint))]] <- id(canon_constr)   # TODO: Check this updates like dict().update in Python
  }

  new_problem <- Problem(canon_objective, canon_constraints)
  return(list(object, new_problem, inverse_data))
})


#'
#' The InverseData class.
#'
#' This class stores the data useful for solution retrieval.
#'
#' @rdname InverseData-class
.InverseData <- setClass("InverseData", representation(problem = "Problem", id_map = "list", var_offsets = "list", x_length = "numeric", var_dims = "list",
                                                       param_dims = "list", param_to_size = "list", param_id_map = "list",
                                                       id2var = "list", id2cons = "list", cons_id_map = "list", constraints = "ListORNULL"),
                         prototype(id_map = list(), var_offsets = list(), x_length = NA_real_, var_dims = list(),
                                   param_dims = list(), param_to_size = list(), param_id_map = list(), id2var = list(), id2cons = list(),
                                   cons_id_map = list(), constraints = NULL))

InverseData <- function(problem) { .InverseData(problem = problem) }

## Add InverseData to class union InverseDataORNUL
setIs("InverseData", "InverseDataORNULL")

setMethod("initialize", "InverseData", function(.Object, ..., problem, id_map = list(), var_offsets = list(), x_length = NA_real_, var_dims = list(), id2var = list(), real2imag = list(), id2cons = list(), cons_id_map = list(), r = NA_real_, minimize = NA, sorted_constraints = list(), is_mip = NA) {
  # Basic variable offset information
  varis <- variables(problem)
  varoffs <- InverseData.get_var_offsets(varis)
  .Object@id_map <- varoffs$id_map
  .Object@var_offsets <- varoffs$var_offsets
  .Object@x_length <- varoffs$vert_offset
  .Object@var_dims <- varoffs$var_dims

  .Object@param_dims <- list()
  # Always start with CONSTANT_ID.
  .Object@param_to_size[[CONSTANT_ID]] <- 1
  .Object@param_id_map <- list()
  offset <- 0
  for(param in parameters(problem)) {
    pid <- as.character(id(param))
    .Object@param_dims[[pid]] <- dim(param)
    .Object@param_to_size[[pid]] <- size(param)
    .Object@param_id_map[[pid]] <- offset
    offset <- offset + size(param)
  }

  # Map of variable id to variable
  .Object@id2var <- stats::setNames(varis, sapply(varis, function(var) { as.character(id(var)) }))

  # Map of constraint id to constraint
  constrs <- problem@constraints
  .Object@id2cons <- stats::setNames(constrs, sapply(constrs, function(cons) { as.character(id(cons)) }))
  .Object@cons_id_map <- list()
  .Object@constraints <- NULL
  return(.Object)
})

InverseData.get_var_offsets <- function(variables) {
  var_dims <- list()
  var_offsets <- list()
  id_map <- list()
  vert_offset <- 0
  for(x in variables) {
    xid <- as.character(id(x))
    var_dims[[xid]] <- dim(x)
    var_offsets[[xid]] <- vert_offset
    id_map[[xid]] <- c(vert_offset, size(x))
    vert_offset <- vert_offset + size(x)
  }
  return(list(id_map = id_map, var_offsets = var_offsets, vert_offset = vert_offset, var_dims = var_dims))
}

#' @param solution A \linkS4class{Solution} to a problem that generated the inverse data.
#' @param inverse_data An \linkS4class{InverseData} object that contains the data encoding the original problem.
#' @describeIn Canonicalization Performs the reduction on a problem and returns an equivalent problem.
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

#' @param expr An \linkS4class{Expression} object.
#' @describeIn Canonicalization Recursively canonicalize an Expression.
setMethod("canonicalize_tree", "Canonicalization", function(object, expr) {
  # TODO: Don't copy affine expressions?
  if(inherits(expr, "PartialProblem")) {
    canon <- canonicalize_tree(object, expr(expr@args[[1]]@objective))
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

#' @param args List of arguments to canonicalize the expression.
#' @describeIn Canonicalization Canonicalize an expression, w.r.t. canonicalized arguments.
setMethod("canonicalize_expr", "Canonicalization", function(object, expr, args) {
  expr_parms <- parameters(expr)
  if(is(expr, "Expression") && is_constant(expr) && (is.null(expr_parms) || length(expr_parms) == 0)) {
    return(list(expr, list()))
  } else if(inherits(expr, names(object@canon_methods)))
    return(object@canon_methods[[class(expr)]](expr, args))   # TODO: Not working for DgpCanonMethods.
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

#' @param x,object A \linkS4class{Chain} object.
#' @rdname Chain-class
setMethod("as.character", "Chain", function(x) { paste(sapply(x@reductions, as.character), collapse = ", ") })

setMethod("show", "Chain", function(object) { paste("Chain(reductions = (", as.character(object@reductions),"))") })

# TODO: How to implement this with S4 setMethod? Signature is x = "DgpCanonMethods", i = "character", j = "missing".
'[[.Chain' <- function(x, i, j, ..., exact = TRUE) { do.call("$", list(x, i)) }

#' @param name The type of reduction.
#' @describeIn Chain Returns the reduction of the specified type.
setMethod("$", signature(x = "Chain"), function(x, name) {
  for(reduction in x@reductions) {
    if(is(reduction, name))
      return(reduction)
  }
  stop(name, " reduction not found")
})

#' @param problem A \linkS4class{Problem} object to check.
#' @describeIn Chain A problem is accepted if the sequence of reductions is valid. In particular, the i-th reduction must accept the output of the i-1th
#' reduction, with the first reduction (self.reductions[0]) in the sequence taking as input the supplied problem.
setMethod("accepts", signature(object = "Chain", problem = "Problem"), function(object, problem) {
  for(i in seq_along(object@reductions)) {
    r <- object@reductions[[i]]
    if(!accepts(r, problem))
      return(FALSE)

    tmp <- perform(r, problem)
    object@reductions[[i]] <- tmp[[1]]
    problem <- tmp[[2]]
  }
  return(TRUE)
})

#' @describeIn Chain Applies the chain to a problem and returns an equivalent problem.
setMethod("perform", signature(object = "Chain", problem = "Problem"), function(object, problem, verbose = FALSE) {
  inverse_data <- list()
  for(i in seq_along(object@reductions)) {
    r <- object@reductions[[i]]
    if(verbose)
      print(paste("Applying reduction", class(r)))
    res <- perform(r, problem)

    object@reductions[[i]] <- res[[1]]
    problem <- res[[2]]
    inv <- res[[3]]
    inverse_data <- c(inverse_data, list(inv))
  }
  return(list(object, problem, inverse_data))
})

#' @param solution A \linkS4class{Solution} or list.
#' @param inverse_data A list that contains the data encoding the original problem.
#' @describeIn Chain Performs the reduction on a problem and returns an equivalent problem.
setMethod("invert", signature(object = "Chain", solution = "SolutionORList", inverse_data = "list"), function(object, solution, inverse_data) {
  m <- min(length(object@reductions), length(inverse_data))
  for(i in rev(seq_len(m))) {
    r <- object@reductions[[i]]
    inv <- inverse_data[[i]]
    solution <- invert(r, solution, inv)
  }
  return(solution)
})

#==================================#
# Convex Attributes to Constraints
#==================================#

# Convex attributes that generate constraints.
CONVEX_ATTRIBUTES <- c("nonneg", "nonpos", "pos", "neg", "symmetric", "diag", "PSD", "NSD")

# Attributes related to symmetry.
SYMMETRIC_ATTRIBUTES <- c("symmetric", "PSD", "NSD")

convex_attributes <- function(variables) {
  # Returns a list of the (constraint-generating) convex attributes present among the variables.
  attributes_present(variables, CONVEX_ATTRIBUTES)
}

attributes_present <- function(variables, attr_map) {
  # Returns a list of the relevant attributes present among the variables.
  # attr_map[sapply(attr_map, function(attr) {
  #  any(sapply(variables, function(v) { !is.null(v@attributes[[attr]]) && v@attributes[[attr]] }))
  # })]

  attr_list <- c()
  for(attr in attr_map) {
    if(any(sapply(variables, function(v) { !is.null(v@attributes[[attr]]) && v@attributes[[attr]] })))
      attr_list <- c(attr_list, attr)
  }
  return(attr_list)
}

recover_value_for_variable <- function(variable, lowered_value, project = TRUE) {
  if(variable@attributes$diag) {
    v_flat <- as.vector(lowered_value)
    return(sparseMatrix(i = seq_along(v_flat), j = seq_along(v_flat), x = v_flat))
  } else if(length(attributes_present(list(variable), SYMMETRIC_ATTRIBUTES)) > 0) {
    n <- nrow(variable)
    value <- matrix(0, nrow = nrow(variable), ncol = ncol(variable))
    triu_mask <- upper.tri(matrix(0, nrow =  n, ncol = n), diag = TRUE)
    value[triu_mask] <- as.vector(lowered_value)
    return(value + t(value) - diag(diag(value)))
  } else if(project)
    return(project(variable, lowered_value))
  else
    return(lowered_value)
}

lower_value <- function(variable, value) {
  if(length(attributes_present(list(variable), SYMMETRIC_ATTRIBUTES)) > 0) {
    # return(upper.tri(value, diag = TRUE))
    n <- nrow(variable)
    triu_mask <- upper.tri(matrix(0, nrow = n, ncol = n), diag = TRUE)
    return(value[triu_mask])
  } else if(variable@attributes$diag)
    return(diag(value))
  else
    return(value)
}

#'
#' The CvxAttr2Constr class.
#'
#' This class represents a reduction that expands convex variable attributes into constraints.
#'
#' @rdname CvxAttr2Constr-class
CvxAttr2Constr <- setClass("CvxAttr2Constr", contains = "Reduction")

setMethod("accepts", signature(object = "CvxAttr2Constr", problem = "Problem"), function(object, problem) { TRUE })

#' @param object A \linkS4class{CvxAttr2Constr} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn CvxAttr2Constr Expand convex variable attributes to constraints.
setMethod("perform", signature(object = "CvxAttr2Constr", problem = "Problem"), function(object, problem) {
  if(length(attributes_present(variables(problem), CONVEX_ATTRIBUTES)) == 0)
    return(list(object, problem, list()))

  # For each unique variable, add constraints.
  id2new_var <- list()
  id2new_obj <- list()
  id2old_var <- list()
  constr <- list()
  for(var in variables(problem)) {
    vid <- id(var)
    vid_char <- as.character(id(var))
    if(!(vid_char %in% names(id2new_var))) {
      id2old_var[[vid_char]] <- var
      new_var <- FALSE
      new_attr <- var@attributes
      for(key in CONVEX_ATTRIBUTES) {
        if(new_attr[[key]]) {
          new_var <- TRUE
          new_attr[[key]] <- FALSE
        }
      }

      if(length(attributes_present(list(var), SYMMETRIC_ATTRIBUTES)) > 0) {
        n <- nrow(var)
        new_dim <- c(floor(n*(n+1)/2), 1)
        # upper_tri <- do.call(Variable, c(list(new_dim), new_attr))
        upper_tri <- do.call(.Variable, c(list(dim = new_dim, var_id = vid), new_attr))
        upper_tri <- set_variable_of_provenance(upper_tri, var)
        id2new_var[[vid_char]] <- upper_tri
        fill_coeff <- Constant(upper_tri_to_full(n))
        full_mat <- fill_coeff %*% upper_tri
        obj <- reshape_expr(full_mat, c(n, n))
      } else if(var@attributes$diag) {
        # diag_var <- do.call(Variable, c(list(nrow(var)), new_attr))
        diag_var <- do.call(.Variable, c(list(dim = c(nrow(var), 1), var_id = vid), new_attr))
        diag_var <- set_variable_of_provenance(diag_var, var)
        id2new_var[[vid_char]] <- diag_var
        obj <- diag(diag_var)
      } else if(new_var) {
        # obj <- do.call(Variable, c(list(dim(var)), new_attr))
        obj <- do.call(.Variable, c(list(dim = dim(var), var_id = vid), new_attr))
        obj <- set_variable_of_provenance(obj, var)
        id2new_var[[vid_char]] <- obj
      } else {
        obj <- var
        id2new_var[[vid_char]] <- obj
      }

      vid_new <- id(var)
      vid_new_char <- as.character(vid_new)
      id2new_obj[[vid_new_char]] <- obj
      if(is_pos(var) || is_nonneg(var))
        constr <- c(constr, obj >= 0)
      else if(is_neg(var) || is_nonpos(var))
        constr <- c(constr, obj <= 0)
      else if(is_psd(var))
        constr <- c(constr, obj %>>% 0)
      else if(var@attributes$NSD)
        constr <- c(constr, obj %<<% 0)
    }
  }

  # Create new problem.
  obj <- tree_copy(problem@objective, id_objects = id2new_obj)
  cons_id_map <- list()
  for(cons in problem@constraints) {
    constr <- c(constr, tree_copy(cons, id_objects = id2new_obj))
    cons_id_map[[as.character(id(cons))]] <- id(constr[[length(constr)]])
  }
  inverse_data <- list(id2new_var, id2old_var, cons_id_map)
  return(list(object, Problem(obj, constr), inverse_data))
})

#' @param object A \linkS4class{CvxAttr2Constr} object.
#' @param solution A \linkS4class{Solution} to a problem that generated the inverse data.
#' @param inverse_data The inverse data returned by an invocation to apply.
#' @describeIn CvxAttr2Constr Performs the reduction on a problem and returns an equivalent problem.
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
    nvid <- as.character(id(new_var))
    if(nvid %in% names(solution@primal_vars))
      pvars[[id]] <- recover_value_for_variable(var, solution@primal_vars[[nvid]])
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
    return(lapply(expr, EvalParams.replace_params_with_consts ))
  else if(is.null(parameters(expr)) || length(parameters(expr)) == 0)
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

#' @param object A \linkS4class{EvalParams} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn EvalParams Replace parameters with constant values.
setMethod("perform", signature(object = "EvalParams", problem = "Problem"), function(object, problem) {
  # Do not instantiate a new objective if it does not contain parameters.
  if(length(parameters(problem@objective)) > 0) {
    obj_expr <- EvalParams.replace_params_with_consts(expr(problem@objective))
    objective <- do.call(class(problem@objective), list(obj_expr))
    # if(inherits(problem@objective, "Maximize"))
    #   objective <- Maximize(obj_expr)
    # else
    #   objective <- Minimize(obj_expr)
  } else
    objective <- problem@objective

  constraints <- list()
  for(c in problem@constraints) {
    args <- list()
    for(arg in c@args)
      args <- c(args, EvalParams.replace_params_with_consts(arg))

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
  return(list(object, Problem(objective, constraints), list()))
})

#' @param object A \linkS4class{EvalParams} object.
#' @param solution A \linkS4class{Solution} to a problem that generated the inverse data.
#' @param inverse_data The inverse data returned by an invocation to apply.
#' @describeIn EvalParams Returns a solution to the original problem given the inverse_data.
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

#' @param object A \linkS4class{FlipObjective} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn FlipObjective Flip a minimization objective to a maximization and vice versa: \eqn{\max f(x) = -\min -f(x)}
setMethod("perform", signature(object = "FlipObjective", problem = "Problem"), function(object, problem) {
  if(inherits(problem@objective, "Maximize"))
    objective <- Minimize
  else
    objective <- Maximize
  problem <- Problem(objective(-expr(problem@objective)), problem@constraints)
  return(list(object, problem, list()))
})

#' @param object A \linkS4class{FlipObjective} object.
#' @param solution A \linkS4class{Solution} to a problem that generated the inverse data.
#' @param inverse_data The inverse data returned by an invocation to apply.
#' @describeIn FlipObjective Map the solution of the flipped problem to that of the original.
setMethod("invert", signature(object = "FlipObjective", solution = "Solution", inverse_data = "list"), function(object, solution, inverse_data) {
  if(!is.na(solution@opt_val) && !is.null(solution@opt_val))
    solution@opt_val <- -solution@opt_val
  return(solution)
})

#'
#' The MatrixStuffing class.
#'
#' This class stuffs a problem into a standard form for a family of solvers.
#'
#' @rdname MatrixStuffing-class
MatrixStuffing <- setClass("MatrixStuffing", contains = "Reduction")

#' @param object A \linkS4class{MatrixStuffing} object.
#' @param problem A \linkS4class{Problem} object to stuff; the arguments of every constraint must be affine.
#' @describeIn MatrixStuffing Returns a stuffed problem. The returned problem is a minimization problem in which every
#' constraint in the problem has affine arguments that are expressed in the form A %*% x + b.
setMethod("perform", signature(object = "MatrixStuffing", problem = "Problem"), function(object, problem) {
  stop("Unimplemented")
})

#' @param object A \linkS4class{MatrixStuffing} object.
#' @param solution A \linkS4class{Solution} to a problem that generated the inverse data.
#' @param inverse_data The data encoding the original problem.
#' @describeIn MatrixStuffing Returns the solution to the original problem given the inverse_data.
setMethod("invert", signature(object = "MatrixStuffing", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  stop("Unimplemented")
})

## setMethod("stuffed_objective", signature(object = "MatrixStuffing", problem = "Problem", inverse_data = "InverseData"), function(object, problem, inverse_data) {
##   stop("Unimplemented")
## })

#'
#' Coalesces bool, int indices for variables.
#'
#' @param variables A list of \linkS4class{Variable} objects.
#' @return Coalesces bool, int indices for variables. The indexing scheme assumes that the variables will be coalesced into
#' a single one-dimensional variable, with each variable being reshaped in Fortran order.
extract_mip_idx <- function(variables) {
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
