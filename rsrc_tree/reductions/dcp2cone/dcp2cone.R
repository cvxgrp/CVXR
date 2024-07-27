## CVXPY SOURCE: cvxpy/reductions/dcp2cone/dcp2cone.py

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

