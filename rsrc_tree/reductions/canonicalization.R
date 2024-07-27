## CVXPY SOURCE: cvxpy/reductions/canonicalization.py
## NOTE how `apply` is renamed to `perform`!
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

