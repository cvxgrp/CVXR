#'
#' The Indicator class.
#'
#' An expression representing the convex function I(constraints) = 0 if constraints hold, and +Inf otherwise.
#'
#' @slot constraints A list of \linkS4class{Constraint}s.
#' @slot err_tol A numeric tolerance for determing whether the constraints hold.
#' @name Indicator-class
#' @aliases Indicator
#' @rdname Indicator-class
.Indicator <- setClass("Indicator", representation(constraints = "list", err_tol = "numeric"),
                                    prototype(constraints = list(), err_tol = 1e-3), contains = "Expression")
Indicator <- function(constraints, err_tol = 1e-3) { .Indicator(args = constraints, err_tol = err_tol) }
indicator <- Indicator

setMethod("initialize", "Indicator", function(.Object, ..., constraints = list(), err_tol = 1e-3) {
  .Object@constraints <- constraints
  .Object@err_tol <- err_tol
  callNextMethod(.Object, ..., args = constraints)
})

#' @describeIn Indicator Expression is constant if all constraints have constant args.
setMethod("is_constant", "Indicator", function(object) {
  all_args <- do.call("+", lapply(object@args, function(c) { c@args }))
  all(sapply(all_args, is_constant))
})

#' @describeIn Indicator Is the expression convex?
setMethod("is_convex", "Indicator", function(object) { TRUE })

#' @describeIn Indicator Is the expression concave?
setMethod("is_concave", "Indicator", function(object) { FALSE })

#' @describeIn Indicator Is the expression log-log convex?
setMethod("is_log_log_convex", "Indicator", function(object) { FALSE })

#' @describeIn Indicator Is the expression log-log concave?
setMethod("is_log_log_concave", "Indicator", function(object) { FALSE })

#' @describeIn Indicator Is the expression positive?
setMethod("is_nonneg", "Indicator", function(object) { TRUE })

#' @describeIn Indicator Is the expression negative?
setMethod("is_nonpos", "Indicator", function(object) { FALSE })

#' @describeIn Indicator Is the expression imaginary?
setMethod("is_imag", "Indicator", function(object) { FALSE })

#' @describeIn Indicator Is the expression complex valued?
setMethod("is_complex", "Indicator", function(object) { FALSE })

#' @describeIn Indicator Returns information needed to reconstruct the expression besides the args.
setMethod("get_data", "Indicator", function(object) { list(object@err_tol) })

#' @describeIn Indicator Returns the (row, col) dimensions of the expression.
setMethod("dim", "Indicator", function(x) { c(1,1) })

#' @describeIn Indicator Returns the string representation of the expression.
setMethod("name", "Indicator", function(x) { cat("Indicator(", x@args, ")") })

#' @describeIn Indicator A list of constraints describe the closure of the region where the expression is finite.
setMethod("domain", "Indicator", function(object) { object@args })

#' @describeIn Indicator Returns the numeric value of the expression.
setMethod("value", "Indicator", function(object) {
  vals <- sapply(object@args, function(cons) { constr_value(cons, tolerance = object@err_tol) })
  ifelse(all(vals), 0.0, Inf)
})

#' @describeIn Indicator Gives the (sub/super)gradient of the expression wrt each variable. Matrix expressions are vectorized, so the gradient is a matrix. NA indicates variable values unknown or outside domain.
setMethod("grad", "Indicator", function(object) {
  # TODO: Implement gradient.
  stop("Unimplemented")
})

#############
#           #
# Linearize #
#           #
#############

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
      grad_var <- grad_map[[as.character(id(var))]]
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

####################
#                  #
# Partial Optimize #
#                  #
####################

#'
#' Partially optimizes the given problem over the specified variables.
#'
#'    Either opt_vars or dont_opt_vars must be given.
#'    If both are given, they must contain all the variables in the problem.
#'
#'    Partial optimize is useful for two-stage optimization and graph implementations.
#'    For example, we can write
#'
#'        x <- Variable(n)
#'        t <- Variable(n)
#'        abs_x <- partial_optimize(Problem(Minimize(sum(t)),
#'                                  list(-t <= x, x <= t)), opt_vars = list(t))
#'
#'    to define the entrywise absolute value of x.
#'
#' @param prob The problem to partially optimize.
#' @param opt_vars The variables to optimize over.
#' @param dont_opt_vars The variables to not optimize over.
#' @param solver The default solver to use for value and grad.
#' @param ... Additional solver specific keyword arguments.
#' @return An expression representing the partial optimization.
#'         Convex for minimization objectives and concave for maximization objectives.
#'
partial_optimize <- function(prob, opt_vars = list(), dont_opt_vars = list(), solver = NA, ...) {
  # One of the two arguments must be specified.
  if((is.null(opt_vars) || length(opt_vars) == 0) && (is.null(dont_opt_vars) || length(dont_opt_vars) == 0))
    stop("partial_optimize called with neither opt_vars nor dont_opt_vars")
  # If opt_vars is not specified, it's the complement of dont_opt_vars.
  else if(is.null(opt_vars) || length(opt_vars) == 0) {
    ids <- sapply(dont_opt_vars, id)
    opt_vars <- lapply(variables(prob), function(var) { if(!(id(var) %in% ids)) var })
  # If dont_opt_vars is not specified, it's the complement of opt_vars.
  } else if(is.null(dont_opt_vars) || length(dont_opt_vars) == 0) {
    ids <- sapply(opt_vars, id)
    dont_opt_vars <- lapply(variables(prob), function(var) { if(!(id(var) %in% ids)) var })
  } else if(!is.null(opt_vars) && length(opt_vars) != 0 && !is.null(dont_opt_vars) && length(dont_opt_vars) != 0) {
    ids <- sapply(c(opt_vars, dont_opt_vars), id)
    for(var in variables(prob)) {
      if(!(id(var) %in% ids))
        stop("If opt_vars and dont_opt_vars are both specified, they must contain all variables in the problem.")
    }
  }

  # Replace the opt_vars in prob with new variables.
  id_to_new_var <- list()
  for(var in opt_vars)
    id_to_new_var[[as.character(id(var))]] <- do.call(".Variable", c(list(dim = dim(var)), var@attributes))
  new_obj <- tree_copy(prob@objective, id_to_new_var)
  new_constrs <- lapply(prob@constraints, function(con) { tree_copy(con, id_to_new_var) })
  new_var_prob <- Problem(new_obj, new_constrs)
  PartialProblem(new_var_prob, opt_vars, dont_opt_vars, solver, ...)
}


#############
#           #
# Scalarize #
#           #
#############

weighted_sum <- function(objectives, weights) {
  num_objs <- length(objectives)
  weighted_objs <- lapply(seq(num_objs), function(i) { objectives[[i]] * weights[i] })
  Reduce("+", weighted_objs)
}

targets_and_priorities <- function(objectives, priorities, targets, limits = list(), off_target = 1e-5) {
  # Combines objectives with penalties within a range between target and limit.
  #
  # Each Minimize objective i has value
  #
  #   priorities[i]*objectives[i] when objectives[i] >= targets[i]
  #
  #   +infinity when objectives[i] > limits[i]
  #
  # Each Maximize objective i has value
  #
  #   priorities[i]*objectives[i] when objectives[i] <= targets[i]
  #
  #   +infinity when objectives[i] < limits[i]
  #
  # Args:
  #   objectives: A list of Minimize/Maximize objectives.
  #   priorities: The weight within the range.
  #   targets: The start (end) of penalty for Minimize (Maximize).
  #   limits: The hard end (start) of penalty for Minimize (Maximize).
  #   off_target: Penalty outside of target.
  #
  # Returns:
  #   A Minimize/Maximize objective.

  num_objs <- length(objectives)
  new_objs <- list()

  for(i in seq(num_objs)) {
    obj <- objectives[[i]]
    sign <- ifelse(is_nonneg(as.Constant(priorities[[i]])), 1, -1)
    off_target <- sign*off_target
    if(inherits(obj, "Minimize")) {
      expr <- (priorities[[i]] - off_target)*Pos(obj@args[[1]] - targets[[i]])
      expr <- expr + off_target*obj@args[[1]]
      if(!is.null(limits) && length(limits) > 0)
        expr <- expr + sign*indicator(list(obj@args[[1]] <= limits[[i]]))
      new_objs <- c(new_objs, expr)
    } else {   # Maximize
      expr <- (priorities[[i]] - off_target)*MinElemwise(obj@args[[1]], targets[[i]])
      expr <- expr + off_target*obj@args[[1]]
      if(length(limits) > 0)
        expr <- expr + sign*indicator(list(obj@args[[1]] >= limits[[i]]))
      new_objs <- c(new_objs, expr)
    }
  }

  obj_expr <- Reduce("+", new_objs)
  if(is_convex(obj_expr))
    return(Minimize(obj_expr))
  else
    return(Maximize(obj_expr))
}

Scalarize.max <- function(objectives, weights) {
  # Combines objectives as max of weighted terms.
  #
  # Args:
  #   objectives: A list of Minimize/Maximize objectives.
  #   weights: A vector of weights.
  #
  # Returns:
  #   A Minimize objective.
  #

  num_objs <- length(objectives)
  expr <- .MaxElemwise(atom_args = lapply(seq(num_objs), function(i) { (objectives[[i]]*weights[i])@args[[1]] }))
  return(Minimize(expr))
}

Scalarize.log_sum_exp <- function(objectives, weights, gamma) {
  # Combines objectives as log_sum_exp of weighted terms.
  #
  # The objective takes the form
  #   log(sum_{i=1}^n exp(gamma*weights[i]*objectives[i]))/gamma
  # As gamma goes to 0, log_sum_exp approaches weighted_sum. As gamma goes to infinity,
  # log_sum_exp approaches max.
  #
  # Args:
  #   objectives: A list of Minimize/Maximize objectives.
  #   weights: A vector of weights.
  #   gamma: Parameter interpolating between weighted_sum and max.
  #
  # Returns:
  #   A Minimize objective.

  num_objs <- length(objectives)
  terms <- lapply(seq(num_objs), function(i) { (objectives[[i]]*weights[i])@args[[1]] })
  expr <- LogSumExp(gamma*VStack(terms))/gamma
  return(Minimize(expr))
}

#####################
#                   #
# Support Functions #
#                   #
#####################

scs_coniclift <- function(x, constraints) {
  # Return (A, b, K) so that
  #      {x : x satisfies constraints}
  # can be written as
  #      {x : exists y where A @ [x; y] + b in K}.
  #
  #  Parameters
  #  ----------
  #  x: Variable
  #  constraints: list of Constraint objects. Each Constraint must be DCP-compatible.
  #
  #  Notes
  #  -----
  #  This function DOES NOT work when x has attributes, like PSD=TRUE, diag=TRUE,
  #  symmetric=TRUE, etc...

  # The objective value is only used to make sure that "x" participates in the
  # problem. So, if constraints is an empty list, then the support function is
  # the standard support function for R^n.
  prob <- Problem(Minimize(sum(x)), constraints)

  tmp <- get_problem_data(prob, solver = "SCS")
  data <- tmp[[1]]
  chain <- tmp[[2]]
  invdata <- tmp[[3]]
  inv <- invdata[length(invdata)-2]
  x_offset <- inv@var_offsets[as.character(id(x))]
  x_indices <- seq(x_offset, x_offset + size(x)) + 1

  A <- data$A
  x_selector <- as.logical(matrix(0, nrow = ncol(A), ncol = 1))
  x_selector[x_indices] <- TRUE
  A_x <- A[,x_selector]
  A_other <- A[,!x_selector]
  A <- -rbind(A_x, A_other)
  b <- data$b
  K <- data$dims
  return(list(A, b, K))
}

scs_cone_selectors <- function(K) {
  # Parse a ConeDims object, as returned from SCS's perform function.
  #
  # Return a dictionary which gives row-wise information for the affine operator
  # returned from SCS's perform function.
  #
  # Parameters
  #   K: ConeDims
  #
  # Returns
  #   selectors: List keyed by strings, which specify cone types. Values are
  #   R vectors or lists of vectors. The vectors give row indices of the affine
  #   operator (A, b) returned by SCS's perform function.

  if(K@p3d) {
    # TODO: Implement this.
    stop("Unimplemented: SuppFunc doesn't yet support feasible sets represented with power cone constraints")
  }

  idx <- K@zero + 1
  nonneg_idxs <- seq(idx, idx + K@nonneg - 1)
  idx <- idx + K@nonneg

  soc_idxs <- list()
  for(soc in K@soc) {
    idxs <- seq(idx, idx + soc - 1)
    soc_idxs <- c(soc_idxs, list(idxs))
    idx <- idx + soc
  }

  psd_idxs <- list()
  for(psd in K@psd) {
    veclen <- psd*floor((psd + 1)/2)
    psd_idxs <- c(psd_idxs, list(seq(idx, idx + veclen - 1)))
    idx <- idx + veclen
  }

  expsize <- 3*K@exp
  exp_idxs <- seq(idx, idx + expsize - 1)
  selectors <- list(nonneg = nonneg_idxs, exp = exp_idxs, soc = soc_idxs, psd = psd_idxs)
  return(selectors)
}

#'
#' The SuppFunc class
#'
#' Given a list of CVXR Constraint objects constraints involving a real CVXR
#' Variable x, consider the convex set
#'
#'
#'    S = \\{ v : \\text{it's possible to satisfy all } \\texttt{constraints}
#'                \\text{ when } \\texttt{x.value} = v \\}.
#'
#' This object represents the *support function* of :math:`S`.
#' This is the convex function
#'
#'    y \\mapsto \\max\\{ \\langle y, v \\rangle : v \\in S \\}.
#'
#' The support function is a fundamental object in convex analysis.
#' It's extremely useful for expressing dual problems using Fenchel duality.
#'
#'  Notes
#'  -----
#'  You are allowed to use CVXR Variables other than x to define constraints,
#'  but the set S only consists of objects (vectors or matrices) with the same
#'  dimensions as x.
#'
#'  It's possible for the support function to take the value +Inf for a fixed
#'  vector y. This is an important point, and it's one reason why support
#'  functions are actually formally defined with the supremum sup rather than
#'  the maximum max.
#'
#' @slot x This \linkS4class{Variable} object cannot have any attributes, such as \code{PSD = TRUE}, \code{nonneg = TRUE}, \code{symmetric = TRUE}, etc...
#' @slot constraints A list of \linkS4class{Constraint}s. Usually, these are constraints over x, and some number of auxiliary CVXR Variables. It is valid to supply \code{constraints = list()}.
#' @examples
#' # If h = SuppFunc(x, constraints), then you can use h just like any other
#' # scalar-valued atom in CVXR. For example, if x was a CVXR Variable with
#' # ndim(x) == 1, you could do the following:
#'
#' z <- Variable(10)
#' A <- matrix(rnorm(size(x)*10), nrow = size(x), ncol = 10)
#' c <- matrix(runif(10), nrow = 10, ncol = 1)
#' objective <- h(A %*% z) - t(c) %*% z
#' prob <- Problem(Minimize(objective), list())
#' result <- solve(prob)
#' @name SuppFunc-class
#' @aliases SuppFunc
#' @rdname SuppFunc-class
SuppFunc <- setClass("SuppFunc", representation(x = "Variable", constraints = "list", A = "numeric", b = "numeric", K_sels = "numeric"),
                                 prototype(A = NA_real_, b = NA_real_, K_sels = NA_real_))

setMethod("initialize", "SuppFunc", function(.Object, x, constraints, A = NA_real_, b = NA_real_, K_sels = NA_real_) {
  if(!is(x, "Variable"))
    stop("The first argument must be an unmodified CVXR Variable object")
  if(any(sapply(CONVEX_ATTRIBUTES, function(attr) { x@attributes[[attr]] })))
    stop("The first argument cannot have any declared attributes")
  for(con in constraints) {
    con_params <- parameters(con)
    if(length(con_params) > 0)
      stop("Convex sets described with Parameter objects are not allowed")
  }

  .Object@x <- x
  .Object@constraints <- constraints
  .Object@A <- NA_real_
  .Object@b <- NA_real_
  .Object@K_sels <- NA_real_
  .Object <- .compute_conic_repr_of_set(.Object)
  return(.Object)
})

# TODO: Is a callable S4 object possible to create? E.g., f = SuppFunc(x, cons); f() --> call.
call.SuppFunc <- function(object, y) {
  sigma_at_y <- SuppFuncAtom(y, object)
  return(sigma_at_y)
}

.compute_conic_repr_of_set <- function(object) {
  if(length(object@constraints) == 0) {
    dummy <- Variable()
    constrs <- list(dummy == 1)
  } else
    constrs <- object@constraints

  Abk <- scs_coniclift(object@x, constrs)
  A <- Abk[[1]]
  b <- Abk[[2]]
  k <- Abk[[3]]
  K_sels <- scs_cone_selectors(K)

  object@A <- A
  object@b <- b
  object@K_sels <- K_sels
  return(object)
}

conic_repr_of_set.SuppFunc <- function(object) {
  return(list(object@A, object@b, object@K_sels))
}
