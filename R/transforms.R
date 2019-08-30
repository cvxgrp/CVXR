#########################################################################################
# TODO: Uncomment once partial optimization and related objects/functions are complete! #
#########################################################################################

# # The indicator I(constraints) = 0 if constraints hold, +Inf otherwise.
# .Indicator <- setClass("Indicator", representation(args = "list", err_tol = "numeric"),
#                                     prototype(args = list(), err_tol = 1e-3), contains = "Expression")
# Indicator <- function(constraints, err_tol = 1e-3) { .Indicator(args = constraints, err_tol = err_tol) }
# indicator <- Indicator
#
# setMethod("initialize", "Indicator", function(.Object, ..., args = list(), err_tol = 1e-3) {
#   .Object@args <- args
#   .Object@err_tol <- err_tol
#   callNextMethod(.Object, ...)
# })
#
# setMethod("is_convex", "Indicator", function(object) { TRUE })
# setMethod("is_concave", "Indicator", function(object) { FALSE })
# setMethod("is_positive", "Indicator", function(object) { TRUE })
# setMethod("is_negative", "Indicator", function(object) { FALSE })
# setMethod("get_data", "Indicator", function(object) { list(object@err_tol) })
# setMethod("size", "Indicator", function(object) { c(1,1) })
# setMethod("name", "Indicator", function(object) { cat("Indicator(", object@args, ")") })
# setMethod("domain", "Indicator", function(object) { object@args })
# setMethod("value", "Indicator", function(object) {
#   vals <- sapply(object@args, function(cons) { value(cons) })
#   if(all(!is.na(vals)))
#     return(0)
#   else
#     return(Inf)
# })
# setMethod("grad", "Indicator", function(object) { stop("Unimplemented") })
# setMethod("canonicalize", "Indicator", function(object) {
#   constraints <- list()
#   for(cons in object@args)
#     constraints <- c(constraints, list(canonical_form(cons)[[2]]))
#   list(create_const(0, c(1,1)), constraints)
# })

#'
#' Affine Approximation to an Expression
#'
#' Gives an elementwise lower (upper) bound for convex (concave) expressions that is tight
#' at the current variable/parameter values. No guarantees for non-DCP expressions.
#' 
#' If f and g are convex, the objective f-g can be (heuristically) minimized using the
#' implementation below of the convex-concave method:
#' 
#' \code{for(iters in 1:N)
#'    solve(Problem(Minimize(f - linearize(g))))}
#'
#' @param expr An \linkS4class{Expression} to linearize.
#' @return An affine expression or \code{NA} if cannot be linearized.
#' @docType methods
#' @rdname linearize
#' @export
linearize <- function(expr) {
  expr <- as.Constant(expr)
  if(is_affine(expr))
    return(expr)
  else {
    tangent <- value(expr)
    if(any(is.na(tangent)))
      stop("Cannot linearize non-affine expression with missing variable values.")
    grad_map <- grad(expr)
    for(var in variables(expr)) {
      grad_var <- grad_map[[as.character(var@id)]]
      if(any(is.na(grad_var)))
        return(NA_real_)
      else if(is_matrix(var)) {
        flattened <- t(Constant(grad_var)) %*% Vec(var - value(var))
        tangent <- tangent + Reshape(flattened, dim(expr))
      } else
        tangent <- tangent + t(Constant(grad_var)) %*% (var - value(var))
    }
  }
  return(tangent)
}

# partial_optimize <- function(prob, opt_vars = list(), dont_opt_vars = list()) {
#   # One of the two arguments must be specified.
#   if(length(opt_vars) == 0 && length(dont_opt_vars) == 0)
#     stop("partial_optimize called with neither opt_vars nor dont_opt_vars")
#   # If opt_vars is not specified, it's the complement of dont_opt_vars.
#   else if(length(opt_vars) == 0) {
#     ids <- sapply(dont_opt_vars, function(var) { id(var) })
#     opt_vars <- lapply(variables(prob), function(var) { if(!(id(var) %in% ids)) var })
#   # If dont_opt_vars is not specified, it's the complement of opt_vars.
#   } else if(length(dont_opt_vars) == 0) {
#     ids <- sapply(opt_vars, function(var) { id(var) })
#     dont_opt_vars <- lapply(variables(prob), function(var) { if(!(id(var) %in% ids)) var })
#   } else if(length(opt_vars) != 0 && length(dont_opt_vars) != 0) {
#     ids <- sapply(c(opt_vars, dont_opt_vars), function(var) { id(var) })
#     for(var in variables(prob)) {
#       if(!(id(var) %in% ids))
#         stop("If opt_vars and dont_opt_vars are both specified, they must contain all variables in the problem.")
#     }
#   }
#   PartialProblem(prob, opt_vars, dont_opt_vars)
# }
#
# .PartialProblem <- setClass("PartialProblem", representation(.prob = "Problem", opt_vars = "list", dont_opt_vars = "list", args = "list"),
#                                               prototype(opt_vars = list(), dont_opt_vars = list(), args = list()),
#                             contains = "Expression")
# PartialProblem <- function(prob, opt_vars, dont_opt_vars) {
#   .PartialProblem(.prob = prob, opt_vars = opt_vars, dont_opt_vars = dont_opt_vars)
# }
#
# setMethod("initialize", "PartialProblem", function(.Object, ..., .prob, opt_vars = list(), dont_opt_vars = list(), args = list()) {
#   .Object@opt_vars <- opt_vars
#   .Object@dont_opt_vars <- dont_opt_vars
#   .Object@args <- list(prob)
#
#   # Replace the opt_vars in prob with new variables.
#   id_to_new_var <- .Object@opt_vars   # TODO: Do I need to make copies like in CVXPY?
#   names(id_to_new_var) <- sapply(.Object@opt_vars, function(var) { id(var) })
#   new_obj <- PartialProblem.replace_new_vars(prob@objective, id_to_new_var)
#   new_constrs <- lapply(prob@constraints, function(con) { PartialProblem.replace_new_vars(con, id_to_new_var) })
#   .Object@.prob <- Problem(new_obj, new_constrs)
#   callNextMethod(.Object, ...)
# })
#
# setMethod("get_data", "PartialProblem", function(object) { c(object@opt_vars, object@dont_opt_vars) })
# setMethod("is_convex", "PartialProblem", function(object) {
#   is_dcp(object@args[[1]]) && class(object@args[[1]]@objective) == "Minimize"
# })
# setMethod("is_concave", "PartialProblem", function(object) {
#   is_dcp(object@args[[1]]) && class(object@args[[1]]@objective) == "Maximize"
# })
# setMethod("is_log_log_convex", "PartialProblem", function(object) {
#   is_dgp(object@args[[1]]) && class(object@args[[1]]@objective) == "Minimize"
# })
# setMethod("is_log_log_concave", "PartialProblem", function(object) {
#   is_dgp(object@args[[1]]) && class(object@args[[1]]@objective) == "Maximize"
# })
# setMethod("is_nonneg", "PartialProblem", function(object) { is_nonneg(object@args[[1]]@objective@args[[1]]) })
# setMethod("is_nonpos", "PartialProblem", function(object) { is_nonpos(object@args[[1]]@objective@args[[1]]) })
# setMethod("is_imag", "PartialProblem", function(object) { FALSE })
# setMethod("is_complex", "PartialProblem", function(object) { FALSE })
# setMethod("dim", "PartialProblem", function(x) { c() })
# setMethod("variables", "PartialProblem", function(object) { variables(object@args[[1]]) })
# setMethod("parameters", "PartialProblem", function(object) { parameters(object@args[[1]]) })
# setMethod("constants", "PartialProblem", function(object) { constants(object@args[[1]]) })
# setMethod("grad", "PartialProblem", function(object) {
#   # Short circuit for constant
#   if(is_constant(object))
#     return(constant_grad(object))
#
#   old_vals <- sapply(variables(object), function(var) { value(var) })
#   names(old_vals) <- sapply(variables(object), function(var) { id(var) })
#   fix_vars <- list()
#   for(var in variables(object)) {
#     if(is.na(value(var)))
#       return(error_grad(object))
#     else
#       fix_vars <- c(fix_vars, list(var == value(var)))
#   }
#   prob <- Problem(object@.prob@objective, c(fix_vars, object@.prob@constraints))
#   sol <- solve(prob)
#
#   # Compute gradient.
#   if(sol$status %in% SOLUTION_PRESENT) {
#     sign <- is_convex(object) - is_concave(object)
#     # Form Lagrangian.
#     lagr <- object@.prob@objective@args[[1]]
#     for(constr in object@.prob@constraints) {
#       # TODO: better way to get constraint expressions.
#       lagr_multiplier <- as.Constant(sign*sol$constr@dual_value)   # TODO: Fix this once result retrieval finished
#       lagr <- lagr + Trace(t(lagr_multiplier) * expr(constr))
#     }
#     grad_map <- grad(lagr)
#     result <- lapply(variables(object), function(var) { grad_map[[id(var)]] })
#     names(result) <- sapply(variables(object), function(var) { id(var) })
#   } else   # Unbounded, infeasible, or solver error.
#     result <- error_grad(object)
#   # Restore the original values to the variables.
#   for(var in variables(object))
#     value(var) <- old_vals[[as.character(id(var))]]
#   result
# })
#
# setMethod("domain", "PartialProblem", function(object) {
#   # Variables optimizd over are replaced in object@.prob
#   obj_expr <- object@.prob@objective@args[[1]]
#   c(object@.prob@constraints, domain(obj_expr))
# })
#
# setMethod("value", "PartialProblem", function(object) {
#   old_vals <- sapply(variables(object), function(var) { value(var) })
#   names(old_vals) <- sapply(variables(object), function(var) { id(var) })
#   for(var in variables(object)) {
#     if(is.na(value(var)))
#       return(NA_real_)
#     else
#       fix_vars <- c(fix_vars, list(var == value(var)))
#   }
#   prob <- Problem(object@args[[1]]@objective, fix_vars)
#   result <- solve(prob)
#
#   # Restore the original values to the variables.
#   for(var in variables(object))
#     value(var) <- old_vals[[as.character(id(var))]]
#   result
# })
#
# PartialProblem.replace_new_vars <- function(obj, id_to_new_var) {
#   if(is(obj, "Variable") && id(obj) %in% names(id_to_new_var))
#     return(id_to_new_var[[id(obj)]])
#   # Leaves outside of optimized variables are preserved.
#   else if(length(obj@args) == 0)
#     return(obj)
#   else if(is(obj, "PartialProblem")) {
#     prob <- obj@args[[1]]
#     new_obj <- PartialProblem.replace_new_vars(prob@objective, id_to_new_var)
#     new_constr <- list()
#     for(constr in prob@constraints)
#       new_constr <- c(new_constr, list(PartialProblem.replace_new_vars(constr, id_to_new_var)))
#     new_args <- list(Problem(new_obj, new_constr))
#     return(copy(obj, new_args))    # TODO: What the heck does the copy function do? Does it even matter in R?
#   # Parent nodes are copied.
#   } else {
#     new_args <- list()
#     for(arg in obj@args)
#       new_args <- c(new_args, list(PartialProblem.replace_new_vars(arg, id_to_new_var)))
#     return(copy(obj, new_args))
#   }
# }
#
# setMethod("canonicalize", "PartialProblem", function(object) {
#   # Canonical form for objective and problem switches from minimize to maximize
#   canon <- canonical_form(object@.prob@objective@args[[1]])
#   obj <- canon[[1]]
#   constrs <- canon[[2]]
#   for(cons in object@.prob@constraints)
#     constrs <- c(constrs, list(canonical_form(cons)[[2]]))
#   list(obj, constrs)
# })
#
# weighted_sum <- function(objective, weights) {
#   num_objs <- length(objective)
#   weighted_objs <- lapply(1:num_objs, function(i) { objectives[[i]] * weights[i] })
#   Reduce("+", weighted_objs)
# }
#
# targets_and_priorities <- function(objectives, priorities, targets, limits = list(), off_target = 1e-5) {
#   num_objs <- length(objectives)
#   new_objs <- list()
#
#   for(i in 1:num_objs) {
#     obj <- objectives[[i]]
#     sign <- ifelse(is_positive(Constant(priorities[[i]])), 1, -1)
#     off_target <- sign * off_target
#     if(class(obj) == "Minimize") {
#       expr <- (priorities[[i]] - off_target)*Pos(obj@args[[1]] - targets[[i]])
#       expr <- expr + off_target*obj@args[[1]]
#       if(length(limits) > 0)
#         expr <- expr + sign*indicator(list(obj@args[[1]] <= limits[[i]]))
#       new_objs <- c(new_objs, expr)
#     } else {   # Maximize
#       expr <- (priorities[[i]] - off_target)*MinElemwise(obj@args[[1]], targets[[i]])
#       expr <- expr + off_target*obj@args[[1]]
#       if(length(limits) > 0)
#         expr <- expr + sign*indicator(list(obj@args[[1]] >= limits[[i]]))
#       new_objs <- c(new_objs, expr)
#     }
#   }
#
#   obj_expr <- Reduce("+", new_objs)
#   if(is_convex(obj_expr))
#     return(Minimize(obj_expr))
#   else
#     return(Maximize(obj_expr))
# }
#
# Scalarize.max <- function(objectives, weights) {
#   num_objs <- length(objectives)
#   expr <- MaxElemwise(lapply(1:num_objs, function(i) { objectives[[i]]*weights[[i]] }))
#   Minimize(expr)
# }
#
# Scalarize.log_sum_exp <- function(objectives, weights, gamma) {
#   num_objs <- length(objectives)
#   terms <- lapply(1:num_objs, function(i) { (objectives[[i]]*weights[[i]])@args[[1]] })
#   expr <- LogSumExp(VStack(terms))
#   Minimize(expr)
# }
#
# get_separable_problems <- function(problem) {
#   # obj_terms contains the terms in the objective functions. We have to
#   # deal with the special case where the objective function is not a sum.
#   if(is(problem@objective@args[[1]], "AddExpression"))
#     obj_terms <- problem@objective@args[[1]]@args
#   else
#     obj_terms <- list(problem@objective@args[[1]])
#
#   # Remove constant terms, which will be appended to the first separable problem.
#   constant_terms <- lapply(obj_terms, function(term) { if(is_constant(term)) term })
#   obj_terms <- lapply(obj_terms, function(term) { if(!is_constant(term)) term })
#
#   constraints <- problem@constraints
#   num_obj_terms <- length(obj_terms)
#   num_terms <- length(obj_terms) + length(constraints)
#
#   # Objective terms and constraints are indexed from 0 to num_terms - 1.
#   # TODO: Finish this section
#   stop("Unimplemented: need set and connected graph operations.")
#
#   # After splitting, construct subproblems from appropriate objective terms and constraints.
#   term_ids_per_subproblem <- rep(list(), num_components)
#   for(i in 1:length(labels)) {
#     label <- labels[[i]]
#     term_ids_per_subproblem[[label]] <- c(term_ids_per_subproblem[[label]], i)
#   }
#   problem_list <- list()
#   for(index in 1:num_components) {
#     terms <- lapply(term_ids_per_subproblem[[index]], function(i) { if(i < num_obj_terms) obj_terms[[i]] })
#     # If we just call sum, we'll have an extra 0 in the objective.
#     if(!is.null(terms))
#       obj <- sum(terms[-1], terms[1])
#     else
#       obj <- Constant(0)
#     constrs <- lapply(term_ids_per_subproblem[[index]], function(i) {
#         if(i >= num_obj_terms) constraints[[i - num_obj_terms]]
#       })
#     # TODO: What should we do in place of copy?
#     problem_list <- c(problem_list, list(Problem(problem@objective.copy(list(obj)), constrs)))
#   }
#
#   # Append constant terms to the first separable problem.
#   if(!is.null(constant_terms)) {
#     # Avoid adding an extra 0 in the objective.
#     sum_constant_terms <- sum(constant_terms[-1], constant_terms[1])
#     if(!is.null(problem_list))
#       problem_list[[1]]@objective@args[[1]] <- problem_list[[1]]@objective@args[[1]] + sum_constant_terms
#     else
#       problem_list <- c(problem_list, list(Problem(problem@objective.copy(list(sum_constant_terms)))))
#   }
#   problem_list
# }
