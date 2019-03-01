setClass("Dgp2Dcp", contains = "Canonicalization")

setMethod("accepts", signature(object = "Dgp2Dcp", problem = "Problem"), function(object, problem) {
  return(is_dgp(problem))
})

setMethod("perform", signature(object = "Dgp2Dcp", problem = "Problem"), function(object, problem) {
  if(!accepts(object, problem))
    stop("The supplied problem is not DGP")
  object@canon_methods <- DgpCanonMethods()
  tmp <- perform(Canonicalization(), problem)
  equiv_problem <- tmp[[1]]
  inverse_data <- tmp[[2]]
  inverse_data@problem <- problem
  return(list(object, equiv_problem, inverse_data))
})

setMethod("canonicalize_expr", signature(object = "Dgp2Dcp", expr = "Expression"), function(object, expr, args) {
  if(class(expr) %in% object@canon_methods)
    return(object@canon_methods[class(expr)](expr, args))
  else
    return(list(copy(expr, args), list()))
})

setMethod("invert", signature(object = "Dgp2Dcp", solution = "Solution", inverse_data = "InverseData"), function(object, solution, inverse_data) {
  solution <- callNextMethod(object, solution, inverse_data)
  if(solution@status == SOLVER_ERROR)
    return(solution)
  for(vid in names(solution@primal_vars))
    solution@primal_vars[vid] <- exp(solution@primal_vars[vid])
  # f(x) = e^{F(u)}
  solution@opt_val <- exp(solution@opt_val)
  return(solution)
})

# Atom canonicalizers
# TODO: Implement sum_largest/sum_smallest.
Dgp2Dcp.CANON_METHODS <- list(AddExpression = add_canon,
                              Constant = constant_canon,
                              DivExpression = div_canon,
                              Exp = exp_canon,
                              EyeMinusInv = eye_minus_inv_canon,
                              GeoMean = geo_mean_canon,
                              Log = log_canon,
                              MulExpression = mulexpression_canon,
                              Multiply = mul_canon,
                              Norm1 = norm1_canon,
                              NormInf = norm_inf_canon,
                              OneMinusPos = one_minus_pos_canon,
                              Parameter = parameter_canon,
                              PfEigenvalue = pf_eigenvalue_canon,
                              Pnorm = pnorm_canon,
                              Power = power_canon,
                              Prod = prod_canon,
                              QuadForm = quad_form_canon,
                              QuadOverLin = quad_over_lin_canon,
                              Trace = trace_canon,
                              SumEntries = sum_canon,
                              Variable = NULL,
                              
                              MaxEntries = EliminatePwl.CANON_METHODS$MaxEntries,
                              MinEntries = EliminatePwl.CANON_METHODS$MinEntries,
                              MaxElemwise = EliminatePwl.CANON_METHODS$MaxElemwise,
                              MinElemwise = EliminatePwl.CANON_METHODS$MinElemwise)

# Canonicalization of DGPs is a stateful procedure, hence the need for a class.
.DgpCanonMethods <- setClass("DgpCanonMethods", representation(.variables = "list"), prototype(.variables = list()), contains = "list")
DgpCanonMethods <- function(...) { .DgpCanonMethods(...) }

setMethod("[", signature(x = "DgpCanonMethods", i = "character", j = "missing", drop = "ANY"), function(x, i, j, ..., drop = TRUE) {
  if(i == "Variable") {
    # TODO: Check scoping of x here is correct.
    variable_canon <- function(variable, args) {
      # Swap out positive variables for unconstrained variables.
      if(id(variable) %in% names(x@.variables))
        return(list(variable, list()))
      else {
        log_variable <- Variable(dim(variable), var_id = id(variable))
        x@.variables[as.character(id(variable))] <- log_variable
        return(list(log_variable), list())
      }
    }
    return(variable_canon)
  } else
    return(Dgp2Dcp.CANON_METHODS[[key]])
})

add_canon <- function(expr, args) {
  if(is_scalar(expr))
    return(list(log_sum_exp(hstack(args)), list()))
  
  rows <- c()
  summands <- sapply(args, function(s) { if(is_scalar(s)) promote(s, dim(expr)) else s })
  if(length(dim(expr)) == 1) {
    for(i in 1:nrow(expr)) {
      summand_args <- lapply(summands, function(summand) { summand[i] })
      rows <- c(rows, log_sum_exp(hstack(summand_args)))
    }
    return(list(reshape(bmat(rows), dim(expr)), list()))
  } else {
    for(i in 1:nrow(expr)) {
      row <- c()
      for(j in 1:ncol(expr)) {
        summand_args <- lapply(summands, function(summand) { summand[i,j] })
        rows <- c(rows, log_sum_exp(hstack(summand_args)))
      }
    }
    return(list(reshape(bmat(rows), dim(expr)), list()))
  }
}

constant_canon <- function(expr, args) {
  # args <- list()
  return(list(Constant(log(value(expr))), list()))
}

div_canon <- function(expr, args) {
  # expr <- NULL
  # x / y == x * y^(-1)
  return(list(args[[1]] - args[[2]], list()))
}

exp_canon <- function(expr, args) {
  # expr <- NULL
  return(list(Exp(args[[1]]), list()))
}

eye_minus_inv_canon <- function(expr, args) {
  X <- args[[1]]
  # (I - X)^(-1) <= T iff there exists 0 <= Y <= T s.t. YX + Y <= Y.
  # Y represents log(Y) here, hence no positivity constraint.
  Y <- Variable(dim(X))
  prod <- matmul(Y, X)
  lhs <- mulexpression_canon(prod, prod@args)[[1]]
  lhs <- lhs + diag(1, nrow(prod))
  return(list(Y, list(lhs <= Y)))
}

geo_mean_canon <- function(expr, args) {
  out <- 0.0
  for(i in 1:length(args[[1]])) {
    x_i <- args[[1]][[i]]
    p_i <- expr@p[i]
    out <- out + p_i * x_i
  }
  return(list((1 / sum(expr@p))*out, list()))
}

log_canon <- function(expr, args) {
  return(list(Log(args[[1]]), list()))
}

mul_canon <- function(expr, args) {
  # expr <- NULL
  return(list(AddExpression(args), list()))
}

mulexpression_canon <- function(expr, args) {
  lhs <- args[[1]]
  rhs <- args[[2]]
  dims <- mul_dims_promote(dim(lhs), dim(rhs))
  lhs_dim <- dims[[1]]
  rhs_dim <- dims[[2]]
  lhs <- reshape(lhs, lhs_dim)
  rhs <- reshape(rhs, rhs_dim)
  rows <- c()
  
  # TODO: Parallelize this for large matrices.
  for(i in 1:nrow(lhs)) {
    row <- c()
    for(j in 1:ncol(rhs)) {
      hstack_args <- lapply(1:ncol(lhs), function(k) { lhs[i,k] + rhs[k,j] })
      row <- c(row, log_sum_exp(hstack(hstack_args)))
    }
    rows <- c(rows, row)
  }
  mat <- bmat(rows)
  if(!all(dim(mat) == dim(expr)))
    mat <- reshape(mat, dim(expr))
  return(list(mat, list()))
}

nonpos_constr_canon <- function(expr, args) {
  if(length(args) != 2)
    stop("Must have exactly 2 arguments")
  return(list(NonPos(args[[1]] - args[[2]], constr_id = id(expr)), list()))
}

norm1_canon <- function(expr, args) {
  if(length(args) != 1)
    stop("Must have exactly 1 argument")
  tmp <- SumEntries(args[[1]], axis = expr@axis, keepdims = expr@keepdims)
  return(sum_canon(tmp, tmp@args))
}

norm_inf_canon <- function(expr, args) {
  if(length(args) != 1)
    stop("Must have exactly 1 argument")
  tmp <- MaxEntries(args[[1]], axis = expr@axis, keepdims = expr@keepdims)
  return(max_canon(tmp, tmp@args))
}

one_minus_pos_canon <- function(expr, args) {
  return(list(Log(expr@ones - Exp(args[[1]])), list()))
}

parameter_canon <- function(expr, args) {
  # args <- list()
  return(list(Parameter(log(value(expr)), name = name(expr)), list()))
}

pf_eigenvalue_canon <- function(expr, args) {
  X <- args[[1]]
  # rho(X) <= lambda iff there exists v s.t. Xv <= lambda v.
  # v and lambda represent log variables, hence no positivity constraints.
  lambd <- Variable()
  v <- Variable(nrow(X))
  lhs <- matmul(X, v)
  rhs <- lambd*v
  lhs <- mulexpression_canon(lhs, lhs@args)[[1]]
  rhs <- mul_canon(rhs, rhs@args)[[1]]
  return(list(lambd, list(lhs <= rhs)))
}

pnorm_canon <- function(expr, args) {
  x <- args[[1]]
  p <- expr@original_p
  if(is.null(dim(x)))
    x <- promote(p, c(1))
  if(is.na(expr@axis) || length(dim(x)) == 1) {
    x <- Vec(x)
    hstack_args <- lapply(x, function(xi) { xi^p })
    return(list((1.0/p) * log_sum_exp(hstack(hstack_args)), list()))
  }
  
  if(expr@axis == 2)
    x <- t(x)
  
  rows <- c()
  for(i in 1:nrow(x)) {
    row <- x[i]
    hstack_args <- lapply(row, function(xi) { xi^p })
    rows <- c(rows, (1.0/p)*log_sum_exp(hstack_args))
  }
  return(list(vstack(rows), list()))
}

power_canon <- function(expr, args) {
  # y = log(x); x^p --> exp(y^p) --> p*log(exp(y)) = p*y.
  return(list(expr@p*args[[1]], list()))
}

prod_canon <- function(expr, args) {
  return(list(SumEntries(args[[1]], axis = expr@axis, keepdims = expr@keepdims), list()))
}

quad_form_canon <- function(expr, args) {
  x <- args[[1]]
  P <- args[[2]]
  elems <- list()
  for(i in 1:nrow(P)) {
    for(j in 1:nrow(P))
      elems <- c(elems, P[i,j] + x[i] + x[j])
  }
  return(list(log_sum_exp(hstack(elems)), list()))
}

quad_over_lin_canon <- function(expr, args) {
  x <- Vec(args[[1]])
  y <- args[[2]]
  numerator <- sum(sapply(x, function(xi) { 2*xi }))
  return(list(numerator - y, list()))
}

sum_canon <- function(expr, args) {
  X <- args[[1]]
  if(is.na(expr@axis)) {
    x <- Vec(X)
    summation <- do.call("sum", args = lapply(x, function(xi) { xi }))
    canon <- add_canon(summation, summation@args)[[1]]
    return(list(reshape(canon, dim(expr)), list()))
  }
  
  if(expr@axis == 2)
    X <- t(X)
  
  rows <- list()
  for(i in 1:nrow(X)) {
    x <- Vec(X[i])
    summation <- do.call("sum", args = lapply(x, function(xi) { xi }))
    canon <- add_canon(summation, summation@args)[[1]]
    rows <- c(rows, canon)
  }
  canon <- hstack(rows)
  return(list(reshape(canon, dim(expr)), list()))
}

trace_canon <- function(expr, args) {
  diag_sum <- sum(Diag(args[[1]]))
  return(add_canon(diag_sum, diag_sum@args))
}

zero_constr_canon <- function(expr, args) {
  if(length(args) != 2)
    stop("Must have exactly 2 arguments")
  return(list(Zero(args[[1]] - args[[2]], constr_id = id(expr)), list()))
}
