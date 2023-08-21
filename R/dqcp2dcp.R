BisectionData <- list("BisectionData" = list("feas_problem", "param", "tighten_lower", "tighten_upper"))

Dqcp2Dcp.get_lazy_and_real_constraints <- function(constraints) {
  lazy_constraints <- list()
  real_constraints <- list()
  for(c in constraints) {
    if(is.function(c))
      lazy_constraints <- c(lazy_constraints, c)
    else
      real_constraints <- c(real_constraints, c)
  }
  return(list(lazy_constraints, real_constraints))
}

#'
#' Reduce DQCP Problem to Parametrized DCP Problem
#'
#' This reduction takes as input a DQCP problem and returns a parameterized
#' DCP problem that can be solved by bisection. Some of the constraints might
#' be lazy, i.e., callables that return a constraint when called. The problem
#' will only be DCP once the lazy constraints are replaced with actual
#' constraints.
#' 
#' Problems emitted by this reduction can be solved with the bisect function.
#'
#' @rdname Dqcp2Cone-class
.Dqcp2Cone <- setClass("Dqcp2Cone", prototype(.bisection_data = NULL), contains = "Canonicalization")
Dqcp2Cone <- function(problem = NULL) { .Dqcp2Cone(problem = problem) }

setMethod("initialize", "Dqcp2Cone", function(.Object, ...) {
  callNextMethod(.Object, ..., canon_methods = Dcp2Cone.CANON_METHODS)
})

#' @param object A \linkS4class{Dqcp2Cone} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn Dqcp2Cone A problem is accepted if it is a minimization and is DQCP.
setMethod("accepts", signature(object = "Dqcp2Cone", problem = "Problem"), function(object, problem) {
  inherits(problem@objective, "Minimize") && is_dqcp(problem)
})

#' @param solution A \linkS4class{Solution} object to invert.
#' @param inverse_data A \linkS4class{InverseData} object containing data necessary for the inversion.
#' @describeIn D2cp2Cone Returns a solution to the original problem given the inverse data.
setMethod("invert", signature(object = "Dqcp2Cone", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  pvars <- list()
  for(vid in inverse_data@id_map) {
    if(vid in solution@primal_vars)
      pvars[[vid]] <- solution@primal_vars[[vid]]
    else
      # Variable was optimized out because it was unconstrained.
      pvars[[vid]] <- 0.0
  }
  return(Solution(solution@status, solution@opt_val, pvars, list(), solution@attr))
}

#' @describeIn Dqcp2Cone Recursively canonicalize the objective and every constraint.
setMethod("perform", signature(object = "Dqcp2Cone", problem = "Problem"), function(object, problem) {
  constraints <- list()
  for(constr in problem@constraints)
    constraints <- c(constraints, Dqcp2Dcp.canonicalize_constraint(object, constr))
  tmp <- Dqcp2Dcp.get_lazy_and_real_constraints(constraints)
  lazy <- tmp[[1]]
  real <- tmp[[2]]
  feas_problem <- Problem(Minimize(0), real)
  feas_problem@.lazy_constraints <- lazy
  
  objective <- problem@objective@expr
  if(is_nonneg(objective))
    t <- Parameter(nonneg = TRUE)
  else if(is_nonpos(objective))
    t <- Parameter(nonpos = TRUE)
  else
    t <- Parameter()
  constraints <- c(constraints, Dqcp2Dcp.canonicalize_constraint(objective <= t))
  
  tmp <- Dqcp2Dcp.get_lazy_and_real_constraints(constraints)
  lazy <- tmp[[1]]
  real <- tmp[[2]]
  param_problem <- Problem(Minimize(0), real)
  param_problem@.lazy_constraints <- lazy
  param_problem@.bisection_data <- BisectionData(feas_problem, t, tighten_fns(objective))
  return(list(param_problem, InverseData(problem)))
})

#' @param expr An \linkS4class{Expression} object.
#' @describeIn "Dqcp2Cone Recursively canonicalize an Expression.
setMethod("canonicalize_tree", "Dqcp2Cone", function(object, expr) {
  tmp <- D2cp2Cone.canon_args(expr)
  canon_args <- tmp[[1]]
  constrs <- tmp[[2]]
  tmp2 <- canonicalize_expr(object, expr, canon_args)
  canon_expr <- tmp2[[1]]
  c <- tmp2[[2]]
  constrs <- c(constrs, c)
  return(list(canon_expr, constrs))
})

#' @describeIn Dqcp2Cone Canonicalize arguments of an expression. Like canonicalize_tree, but preserves signs.
setMethod("canon_args", "Dqcp2Cone", function(object, expr) {
  canon_args <- list()
  constrs <- list()
  for(arg in expr@args) {
    tmp <- canonicalize_tree(object, arg)
    canon_args <- tmp[[1]]
    c <- tmp[[2]]
    if(is(canon_args, "Variable")) {
      if(is_nonneg(arg))
        canon_arg@attributes$nonneg <- TRUE
      else if(is_nonpos(arg))
        canon_arg@attributes$nonpos <- TRUE
    }
    canon_args <- c(canon_args, canon_arg)
    constrs <- c(constrs, c)
  }
  return(list(canon_args, constrs))
})

#' @describeIn Dqcp2Cone Recursively canonicalize a constraint. The DQCP grammar has expressions of the form INCR QCVX DCP and DECR QCCV DCP, i.e., zero or more real/scalar increasing (or decreasing) atoms, composed with a quasiconvex (or quasiconcave) atom, composed with DCP expressions. The monotone functions are inverted by applying their inverses to both sides of a constraint. The QCVX (QCCV) atom is lowered by replacing it with its sublevel (superlevel) set. The DCP expressions are canonicalized via graph implementations.
setMethod("canonicalize_constraint", "Dqcp2Cone", function(object, expr) {
  lhs <- constr@args[[1]]
  rhs <- constr@args[[2]]
  
  if(is(constr, "Inequality")) {
    # Taking inverses can yield +/-; this is handled here.
    lhs_val <- as.vector(value(lhs))
    rhs_val <- as.vector(value(rhs))
    if(all(lhs_val == -Inf) || all(rhs_val == Inf)) {
      # Constraint is redundant.
      return(list(TRUE))
    } else if(any(lhs_val == Inf) || any(rhs_val == -Inf)) {
      # Constraint is infeasible.
      return(list(FALSE))
    }
  }
  
  if(is_dcp(constr)) {
    tmp <- canonicalize_tree(object, constr)
    canon_constr <- tmp[[1]]
    aux_constr <- tmp[[2]]
    constr_list <- c(list(canon_constr), aux_constr)
    return(constr_list)
  }
  
  # Canonicalize lhs <= rhs.
  # Either lhs or rhs is quasiconvex (and not convex).
  if(!is(constr, "Inequality"))
    stop("constr must be of class Inequality")
  
  # Short-circuit zero-valued expressions to simplify inverse logic.
  if(is_zero(lhs))
    return(canonicalize_constraint(object, 0 <= rhs))
  if(is_zero(rhs))
    return(canonicalize_constraint(object, lhs <= 0))
  
  if(is_quasiconvex(lhs) && !is_convex(lhs)) {
    # Quasiconvex <= constant.
    if(!is_constant(rhs))
      stop("rhs must be a constant")
    if(invertible(lhs)) {
      # Apply inverse to both sides of constraint.
      rhs <- (inverse(lhs))(rhs)
      idx <- .non_const_idx(lhs)[1]
      expr <- lhs@args[[idx]]
      if(is_incr(lhs, idx))
        return(canonicalize_constraint(object, expr <= rhs))
      if(!is_decr(lhs, idx))
        stop("lhs must be decreasing")
      return(canonicalize_constraint(object, expr >= rhs))
    } else if(is(lhs, "MaxEntries") || is(lhs, "MaxElemwise")) {
      # Lower maximum.
      res <- list()
      for(arg in lhs@args)
        res <- c(res, canonicalize_constraint(object, arg <= rhs))
      return(res)
    } else {
      # Replace quasiconvex atom with a sublevel set.
      tmp <- canon_args(object, lhs)
      canon_args <- tmp[[1]]
      aux_args_constr <- tmp[[2]]
      sublevel_set <- sublevel(copy(lhs, canon_args), t = rhs)
      return(c(sublevel_set, aux_args_constr))
    }
  }
  
  # Constant <= quasiconcave.
  if(!is_quasiconcave(rhs))
    stop("rhs must be quasiconcave")
  if(!is_constant(lhs))
    stop("lhs must be constant")
  if(invertible(rhs)) {
    # Apply inverse to both sides of constraint.
    lhs <- (inverse(rhs))(lhs)
    idx <- .non_const_idx(rhs)[1]
    expr <- rhs@args[[idx]]
    if(is_incr(rhs, idx))
      return(canonicalize_constraint(object, lhs <= expr))
    if(!is_decr(rhs, idx))
      stop("rhs must be decreasing in index ", idx)
    return(canonicalize_constraint(object, lhs >= expr))
  } else if(is(rhs, "MinEntries") || is(lhs, "MinElemwise")) {
    # Lower minimum.
    res <- list()
    for(arg in rhs@args)
      res <- c(res, canonicalize_constraint(lhs <= expr))
    return(res)
  } else {
    # Replace quasiconcave atom with a superlevel set.
    tmp <- canon_args(object, rhs)
    canon_args <- tmp[[1]]
    aux_args_constr <- tmp[[2]]
    superlevel_set <- superlevel(copy(rhs, canon_args), t = lhs)
    return(c(superlevel_set, aux_args_constr))
  }
})

# These atoms are always invertible. Others (like AddExpression, DivExpression, SumEntries, and CumSum) are only invertible in special cases, checked in the invertible function.
INVERTIBLE_FNS <- c("Ceil", "Floor", "NegExpression", "Exp", "Log", "Log1p", "Logistic", "Power", "Abs")

# Inverses are extended-value functions.
inverse <- function(expr) {
  if(class(expr) == "Ceil")
    return(function(t) { Floor(t) })
  else if(class(expr) == "Floor")
    return(function(t) { Ceil(t) })
  else if(class(expr) == "NegExpression")
    return(function(t) { -t })
  else if(class(expr) == "Exp") {
    exp_inv <- function(t) {
      if(is_nonneg(t))
        return(Log(t))
      else
        return(-Inf)
    }
    return(exp_inv)
  } else if(class(expr) == "Log")
    return(function(t) { Exp(t) })
  else if(class(expr) == "Log1p")
    return(function(t) { Exp(t) - 1 })
  else if(class(expr) == "Logistic") {
    logistic_inv <- function(t) {
      if(is_nonneg(t))
        return(Log(Exp(t) - 1))
      else
        return(-Inf)
    }
    return(logistic_inv)
  } else if(class(expr) == "Power") {
    power_inv <- function(t) {
      if(value(expr@p) == 1)
        return(t)
      else if(is_nonneg(t))
        return(Power(t, 1/value(expr@p)))
      else
        return(Inf)
    }
    return(power_inv)
  } else if(class(expr) == "Multiply") {
    if(is_constant(expr@args[[1]]))
      const <- expr@args[[1]]
    else
      const <- expr@args[[2]]
    return(function(t) { t / const })
  } else if(class(expr) == "DivExpression") {
    # Either const / x <= t or x / const <= t.
    if(is_constant(expr@args[[1]])) {
      # Numerator is constant.
      const <- expr@args[[1]]
      return(function(t) { const / t })
    } else {
      # Denominator is constant.
      const <- expr@args[[2]]
      return(function(t) { const * t })
    } else if(class(expr) == "AddExpression") {
      if(is_constant(expr@args[[1]]))
        const <- expr@args[[1]]
      else
        const <- expr@args[[2]]
      return(function(t) { t - const })
    } else if(class(expr) == "Abs") {
      arg <- expr@args[[1]]
      if(is_nonneg(arg))
        return(function(t) { t })
      else if(is_nonpos(arg))
        return(function(t) { -t })
      else
        stop("Sign of argument must be known")
    } else if(class(expr) %in% c("SumEntries", "CumSum"))
      return(function(t) { t })
    else
      stop("Expression cannot be inverted")
  }
}

invertible <- function(expr) {
  if(is(expr, "Multiply") || is(expr, "DivExpression") || is(expr, "AddExpression"))
    return(length(.non_const_idx(expr)) == 1)
  else if(is(expr, "SumEntries") || is(expr, "CumSum"))
    return(.is_real(expr))
  else
    return(class(expr) %in% INVERTIBLE_FNS)
}

# Sublevel sets for quasiconvex atoms.
# In the following functions, FALSE is a placeholder for an infeasible constraint (one that cannot be represented in a DCP way), and TRUE is a placeholder for the absence of a constraint.

dist_ratio_sub <- function(expr, t) {
  x <- expr@args[[1]]
  a <- expr@a
  b <- expr@b
  
  sublevel_set <- function() {
    if(value(t) > 1)
      return(FALSE)
    tsq <- value(t)^2
    return(((1 - tsq^2)*SumSquares(x) - (2*(a - tsq*b)) %*% x + SumSquares(a) - tsq*SumSquares(b)) <= 0)
  }
  return(list(sublevel_set))
}

mul_sup <- function(expr, t) {
  x <- expr@args[[1]]
  y <- expr@args[[2]]
  if(is_nonneg(x) && is_nonneg(y))
    return(list(x >= t * InvPos(y)))
  else if(is_nonpos(x) && is_nonpos(y))
    return(list(-x >= t * InvPos(-y)))
  else
    stop("Incorrect signs")
}

mul_sub <- function(expr, t) {
  x <- expr@args[[1]]
  y <- expr@args[[2]]
  if(is_nonneg(x) && is_nonpos(y))
    return(list(y <= t * InvPos(x)))
  else if(is_nonpos(x) && is_nonneg(y))
    return(list(x <= t * InvPos(y)))
  else
    stop("Incorrect signs")
}

ratio_sup <- function(expr, t) {
  x <- expr@args[[1]]
  y <- expr@args[[2]]
  if(is_nonneg(y))
    return(list(x >= t * y))
  else if(is_nonpos(y))
    return(list(x <= t * y))
  else
    stop("The denominator's sign must be known")
}

ratio_sub <- function(expr, t) {
  x <- expr@args[[1]]
  y <- expr@args[[2]]
  if(is_nonneg(y))
    return(list(x <= t * y))
  else if(is_nonpos(y))
    return(list(x >= t * y))
  else
    stop("The denominator's sign must be known")
}

length_sub <- function(expr, t) {
  arg <- expr@args[[1]]
  if(is(t, "Parameter")) {
    sublevel_set <- function() {
      if(value(t) < 0)
        return(FALSE)
      if(value(t) >= size(arg))
        return(TRUE)
      idx_start <- as.integer(value(Floor(t))) + 1
      return(arg[seq(idx_start, size(arg))] == 0)   # TODO: Check if we should use nrow(arg) or size(arg).
    }
    return(list(sublevel_set))
  } else {
    idx_start <- as.integer(value(Floor(t))) + 1
    return(list(arg[seq(idx_start, size(arg))] == 0))   # TODO: Check if we should use nrow(arg) or size(arg).
  }
}

sign_sup <- function(expr, t) {
  x <- expr@args[[1]]
  
  superlevel_set <- function() {
    if(value(t) <= -1)
      return(TRUE)
    else if(value(t) <= 1)
      return(x >= 0)
    else
      return(FALSE)
  }
  return(list(superlevel_set))
}

sign_sub <- function(expr, t) {
  x <- expr@args[[1]]
  
  sublevel_set <- function() {
    if(value(t) >= 1)
      return(TRUE)
    else if(value(t) >= -1)
      return(x <= 0)
    else
      return(FALSE)
  }
  return(list(sublevel_set))
}

gen_lambda_max_sub <- function(expr, t) {
  return(list(expr@args[[1]] == t(expr@args[[1]]), expr@args[[2]] %>>% 0, (t * expr@args[[2]] - expr@args[[1]]) %>>% 0))
}

condition_number_sub <- function(expr, t) {
  A <- expr@args[[1]]
  n <- nrow(A)
  u <- Variable(pos = TRUE)
  
  prom_ut <- promote(u * t, c(n, 1))
  prom_u <- promote(u, c(n, 1))
  tmp_expr1 <- A - diag_vec(prom_u)
  tmp_expr2 <- diag_vec(prom_ut) - A
  
  return(list(upper_tri(A) == upper_tri(t(A)), PSD(A), PSD(tmp_expr1), PSD(tmp_expr2)))
}

SUBLEVEL_SETS <- list(
  "Multiply" = mul_sub,
  "DivExpression" = ratio_sub,
  "VecLength" = length_sub, 
  "Sign" = sign_sub,
  "DistRatio" = dist_ratio_sub,
  "GenLambdaMax" = gen_lambda_max_sub,
  "ConditionNumber" = condition_number_sub
)

SUPERLEVEL_SETS <- list(
  "Multiply" = mul_sup,
  "DivExpression" = ratio_sup,
  "Sign" = sign_sup
)

# Return the t-level sublevel set for expr.
# Returned as a constraint phi_t(x) <= 0, where phi_t(x) is convex.
sublevel <- function(expr, t) {
  if(class(expr) %in% names(SUBLEVEL_SETS))
    return(SUBLEVEL_SETS[[class(expr)]](expr, t))
  else
    stop("Unimplemented: The ", class(expr), " atom is not yet supported in DQCP")
}

# Return the t-level superlevel set for expr.
# Returned as a constraint phi_t(x) >= 0, where phi_t(x) is concave.
superlevel <- function(expr, t) {
  if(class(expr) %in% names(SUPERLEVEL_SETS))
    return(SUPERLEVEL_SETS[[class(expr)]](expr, t))
  else
    stop("Unimplemented: The ", class(expr), " atom is not yet supported in DQCP")
}

integer_valued_fns <- c("Ceil", "Floor", "VecLength")

# Tuples fns such that t infeasible implies fns[1](t) infeasible
# (or sup of infeasible set), t feasible implies fns[2](t)
# (or inf of infeasible set).
tighten_fns <- function(expr) {
  if(class(expr) %in% integer_valued_fns)
    list(base::ceiling, base::floor)
  else if(is_nonneg(expr))
    list(function(t) { base::pmax(t, 0) }, function(t) { t })
  else if(is_nonpos(expr))
    list(function(t) { t }, function(t) { base::pmin(t, 0) })
  else
    list(function(t) { t }, function(t) { t })
}
