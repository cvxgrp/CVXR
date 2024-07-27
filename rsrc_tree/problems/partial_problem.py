#'
#' A Partial Optimization Problem
#'
#' @slot opt_vars The variables to optimize over.
#' @slot dont_opt_vars The variables to not optimize over.
#' @name PartialProblem-class
#' @aliases PartialProblem
#' @rdname PartialProblem-class
.PartialProblem <- setClass("PartialProblem", representation(.prob = "Problem", opt_vars = "list", dont_opt_vars = "list", solver = "ANY", .solve_kwargs = "list"),
                                              prototype(opt_vars = list(), dont_opt_vars = list(), solver = NA_character_, .solve_kwargs = list()),
                             contains = "Expression")

PartialProblem <- function(prob, opt_vars, dont_opt_vars, solver, ...) {
  .PartialProblem(.prob = prob, opt_vars = opt_vars, dont_opt_vars = dont_opt_vars, solver = solver, .solve_kwargs = list(...))
}

setMethod("initialize", "PartialProblem", function(.Object, ..., .prob, opt_vars = list(), dont_opt_vars = list(), solver = NA_character_, .solve_kwargs = list()) {
  .Object@.prob <- .prob
  .Object@opt_vars <- opt_vars
  .Object@dont_opt_vars <- dont_opt_vars
  .Object@solver <- solver
  .Object@.solve_kwargs <- .solve_kwargs
  callNextMethod(.Object, ..., args = list(.prob))
})

#' @describeIn PartialProblem Returns info needed to reconstruct the expression besides the args.
setMethod("get_data", "PartialProblem", function(object) { list(object@opt_vars, object@dont_opt_vars, object@solver) })

#' @describeIn PartialProblem Is the expression constant?
setMethod("is_constant", "PartialProblem", function(object) {
  length(variables(object@args[[1]])) == 0
})

#' @describeIn PartialProblem Is the expression convex?
setMethod("is_convex", "PartialProblem", function(object) {
  is_dcp(object@args[[1]]) && inherits(object@args[[1]]@objective, "Minimize")
})

#' @describeIn PartialProblem Is the expression concave?
setMethod("is_concave", "PartialProblem", function(object) {
  is_dcp(object@args[[1]]) && inherits(object@args[[1]]@objective, "Maximize")
})

#' @describeIn PartialProblem Is the expression a disciplined parameterized expression?
setMethod("is_dpp", "PartialProblem", function(object, context = "dcp") {
  if(tolower(context) %in% c("dcp", "dgp"))
    is_dpp(object@args[[1]], context)
  else
    stop("Unsupported context ", context)
})

#' @describeIn PartialProblem Is the expression log-log convex?
setMethod("is_log_log_convex", "PartialProblem", function(object) {
  is_dgp(object@args[[1]]) && inherits(object@args[[1]]@objective, "Minimize")
})

#' @describeIn PartialProblem Is the expression log-log concave?
setMethod("is_log_log_concave", "PartialProblem", function(object) {
  is_dgp(object@args[[1]]) && inherits(object@args[[1]]@objective, "Maximize")
})

#' @describeIn PartialProblem Is the expression nonnegative?
setMethod("is_nonneg", "PartialProblem", function(object) {
  is_nonneg(object@args[[1]]@objective@args[[1]])
})

#' @describeIn PartialProblem Is the expression nonpositive?
setMethod("is_nonpos", "PartialProblem", function(object) {
  is_nonpos(object@args[[1]]@objective@args[[1]])
})

#' @describeIn PartialProblem Is the expression imaginary?
setMethod("is_imag", "PartialProblem", function(object) { FALSE })

#' @describeIn PartialProblem Is the expression complex valued?
setMethod("is_complex", "PartialProblem", function(object) { FALSE })

#' @describeIn PartialProblem Returns the (row, col) dimensions of the expression.
setMethod("dim", "PartialProblem", function(x) { c(1,1) })

setMethod("name", "PartialProblem", function(x) { cat("PartialProblem(", x@args[[1]], ")") })

#' @describeIn PartialProblem Returns the variables in the problem.
setMethod("variables", "PartialProblem", function(object) { variables(object@args[[1]]) })

#' @describeIn PartialProblem Returns the parameters in the problem.
setMethod("parameters", "PartialProblem", function(object) { parameters(object@args[[1]]) })

#' @describeIn PartialProblem Returns the constants in the problem.
setMethod("constants", "PartialProblem", function(object) { constants(object@args[[1]]) })

#' @describeIn PartialProblem Gives the (sub/super)gradient of the expression wrt each variable. Matrix expressions are vectorized, so the gradient is a matrix. NA indicates variable values unknown or outside domain.
setMethod("grad", "PartialProblem", function(object) {
  # Subgrad of g(y) = min f_0(x,y)
  #                   s.t. f_i(x,y) <= 0, i = 1,..,p
  #                        h_i(x,y) == 0, i = 1,...,q
  # Given by Df_0(x^*,y) + \sum_i Df_i(x^*,y) \lambda^*_i
  #          + \sum_i Dh_i(x^*,y) \nu^*_i
  # where x^*, \lambda^*_i, \nu^*_i are optimal primal/dual variables.
  # Add PSD constraints in same way.

  # Short circuit for constant
  if(is_constant(object))
    return(constant_grad(object))

  old_vals <- list()
  for(var in variables(object))
    old_vals[[as.character(id(var))]] <- value(var)

  fix_vars <- list()
  for(var in object@dont_opt_vars) {
    if(is.na(value(var)))
      return(error_grad(object))
    else
      fix_vars <- c(fix_vars, list(var == value(var)))
  }

  prob <- Problem(object@args[[1]]@objective, c(fix_vars, object@args[[1]]@constraints))
  result <- do.call("solve", c(list(a = prob, solver = object@solver), object@.solve_kwargs))

  # Compute gradient.
  if(result$status %in% SOLUTION_PRESENT) {
    sign <- as.numeric(is_convex(object) - is_concave(object))

    # Form Lagrangian.
    lagr <- object@args[[1]]@objective@args[[1]]
    for(constr in object@args[[1]]@constraints) {
      # TODO: better way to get constraint expressions.
      lagr_multiplier <- as.Constant(sign * result$getDualValue(constr))
      lprod <- t(lagr_multiplier) %*% constr@expr
      if(is_scalar(lprod))
        lagr <- lagr + sum(lprod)
      else
        lagr <- lagr + matrix_trace(lprod)
    }

    grad_map <- grad(lagr)   # TODO: After finishing grad implementation, we probably need to update this call by passing in result of solving problem.
    res_grad <- list()
    for(var in object@dont_opt_vars) {
      var_id_char <- as.character(id(var))
      res_grad[[var_id_char]] <- grad_map[[var_id_char]]
    }
  } else   # Unbounded, infeasible, or solver error.
    res_grad <- error_grad(object)

  # Restore the original values to the variables.
  for(var in variables(object))
    value(var) <- old_vals[[as.character(id(var))]]
  return(res_grad)
})

#' @describeIn PartialProblem A list of constraints describing the closure of the region where the expression is finite.
setMethod("domain", "PartialProblem", function(object) {
  # Variables optimized over are replaced in object@args[[1]]
  obj_expr <- object@args[[1]]@objective@args[[1]]
  c(object@args[[1]]@constraints, domain(obj_expr))
})

#' @describeIn PartialProblem Returns the numeric value of the expression.
setMethod("value", "PartialProblem", function(object) {
  old_vals <- list()
  for(var in variables(object))
    old_vals[[as.character(id(var))]] <- value(var)

  fix_vars <- list()
  for(var in object@dont_opt_vars) {
    if(is.na(value(var)))
      return(NA_real_)
    else
      fix_vars <- c(fix_vars, list(var == value(var)))
  }

  prob <- Problem(object@args[[1]]@objective, c(fix_vars, object@args[[1]]@constraints))
  result <- do.call("solve", c(list(a = prob, solver = object@solver), object@.solve_kwargs))

  # Restore the original values to the variables.
  for(var in variables(object))
    value(var) <- old_vals[[as.character(id(var))]]

  # Need to get value returned by solver in case of stacked partial_optimizes.
  return(result$value)
})

#' @describeIn PartialProblem Returns the graph implementation of the object. Chain ids of all the opt_vars.
setMethod("canonicalize", "PartialProblem", function(object) {
  # Canonical form for objective and problem switches from minimize to maximize.
  canon <- canonical_form(object@args[[1]]@objective@args[[1]])
  obj <- canon[[1]]
  constrs <- canon[[2]]
  for(cons in object@args[[1]]@constraints)
    constrs <- c(constrs, list(canonical_form(cons)[[2]]))
  list(obj, constrs)
})

