## CVXPY SOURCE: cvxpy/reductions/eval_params.py

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

