#'
#' The EliminatePwl class.
#' 
#' This class eliminates piecewise linear atoms.
#' 
#' @rdname EliminatePwl-class
.EliminatePwl <- setClass("EliminatePwl", contains = "Canonicalization")
EliminatePwl <- function(problem = NULL) { .EliminatePwl(problem = problem) }

setMethod("initialize", "EliminatePwl", function(.Object, ...) {
  callNextMethod(.Object, ..., canon_methods = EliminatePwl.CANON_METHODS)
})

#' @param object An \linkS4class{EliminatePwl} object.
#' @param problem A \linkS4class{Problem} object.
#' @describeIn EliminatePwl Does this problem contain piecewise linear atoms?
setMethod("accepts", signature(object = "EliminatePwl", problem = "Problem"), function(object, problem) {
  atom_types <- sapply(atoms(problem), function(atom) { class(atom) })
  pwl_types <- c("Abs", "MaxElemwise", "SumLargest", "MaxEntries", "Norm1", "NormInf")
  return(any(sapply(atom_types, function(atom) { atom %in% pwl_types })))
})

setMethod("perform", signature(object = "EliminatePwl", problem = "Problem"), function(object, problem) {
  if(!accepts(object, problem))
    stop("Cannot canonicalize away piecewise linear atoms.")
  callNextMethod(object, problem)
})

# Atom canonicalizers.
#' 
#' EliminatePwl canonicalizer for the absolute atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A canonicalization of the picewise-lienar atom
#' constructed from an absolute atom where the objective function 
#' consists of the variable that is of the same dimension as the 
#' original expression and the constraints consist of splitting 
#' the absolute value into two inequalities.
#' 
EliminatePwl.abs_canon <- function(expr, args) {
  x <- args[[1]]
  # t <- Variable(dim(expr))
  t <- new("Variable", dim = dim(expr))
  constraints <- list(t >= x, t >= -x)
  return(list(t, constraints))
}

#' 
#' EliminatePwl canonicalizer for the cumulative max atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A canonicalization of the piecewise-lienar atom
#' constructed from a cumulative max atom where the objective
#' function consists of the variable Y which is of the same
#' dimension as the original expression and the constraints
#' consist of row/column constraints depending on the axis
EliminatePwl.cummax_canon <- function(expr, args) {
  X <- args[[1]]
  axis <- expr@axis
  
  # Implicit O(n) definition:
  # Y_{k} = maximum(Y_{k-1}, X_k)
  # Y <- Variable(dim(expr))
  Y <- new("Variable", dim = dim(expr))
  constr <- list(X <= Y)
  if(axis == 2) {
    if(nrow(Y) > 1)
      constr <- c(constr, list(Y[1:(nrow(Y)-1),] <= Y[2:nrow(Y),]))
  } else {
    if(ncol(Y) > 1)
      constr <- c(constr, list(Y[,1:(ncol(Y)-1)] <= Y[,2:ncol(Y)]))
  }
  return(list(Y, constr))
}

#' 
#' EliminatePwl canonicalizer for the cumulative sum atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A canonicalization of the piecewise-lienar atom
#' constructed from a cumulative sum atom where the objective
#' is Y that is of the same dimension as the matrix of the expression
#' and the constraints consist of various row constraints
EliminatePwl.cumsum_canon <- function(expr, args) {
  X <- args[[1]]
  axis <- expr@axis
  
  # Implicit O(n) definition:
  # X = Y[1,:] - Y[2:nrow(Y),:]
  # Y <- Variable(dim(expr))
  Y <- new("Variable", dim = dim(expr))
  if(axis == 2) {  # Cumulative sum on each column
    constr <- list(Y[1,] == X[1,])
    if(nrow(Y) > 1)
      constr <- c(constr, list(X[2:nrow(X),] == Y[2:nrow(Y),] - Y[1:(nrow(Y)-1),]))
  } else {   # Cumulative sum on each row
    constr <- list(Y[,1] == X[,1])
    if(ncol(Y) > 1)
      constr <- c(constr, list(X[,2:ncol(X)] == Y[,2:ncol(Y)] - Y[,1:(ncol(Y)-1)]))
  }
  return(list(Y, constr))
}

#' 
#' EliminatePwl canonicalizer for the max entries atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A canonicalization of the piecewise-lienar atom
#' constructed from the max entries atom where the objective
#' function consists of the variable t of the same size as
#' the original expression and the constraints consist of
#' a vector multiplied by a vector of 1's.
EliminatePwl.max_entries_canon <- function(expr, args) {
  x <- args[[1]]
  axis <- expr@axis
  # expr_dim <- dim(expr)
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = dim(expr))
  
  if(is.na(axis))   # dim(expr) = c(1,1)
    promoted_t <- promote(t, dim(x))
  else if(axis == 2)   # dim(expr) = c(1,n)
    promoted_t <- Constant(matrix(1, nrow = nrow(x), ncol = 1) %*% reshape_expr(t, c(1, ncol(x))))
  else   # shape = c(m,1)
    promoted_t <- reshape_expr(t, c(nrow(x), 1)) %*% Constant(matrix(1, nrow = 1, ncol = ncol(x)))
  
  constraints <- list(x <= promoted_t)
  return(list(t, constraints))
}

#' 
#' EliminatePwl canonicalizer for the elementwise maximum atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A canonicalization of the piecewise-lienar atom
#' constructed by a elementwise maximum atom where the
#' objective function is the variable t of the same dimension
#' as the expression and the constraints consist of a simple
#' inequality.
EliminatePwl.max_elemwise_canon <- function(expr, args) {
  # expr_dim <- dim(expr)
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = dim(expr))
  constraints <- lapply(args, function(elem) { t >= elem })
  return(list(t, constraints))
}

#' 
#' EliminatePwl canonicalizer for the minimum entries atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A canonicalization of the piecewise-lienar atom
#' constructed by a minimum entries atom where the
#' objective function is the negative of variable 
#' t produced by max_elemwise_canon of the same dimension
#' as the expression and the constraints consist of a simple
#' inequality.
EliminatePwl.min_entries_canon <- function(expr, args) {
  if(length(args) != 1)
    stop("Length of args must be one")
  tmp <- MaxEntries(-args[[1]])
  canon <- EliminatePwl.max_entries_canon(tmp, tmp@args)
  return(list(-canon[[1]], canon[[2]]))
}

#' 
#' EliminatePwl canonicalizer for the elementwise minimum atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A canonicalization of the piecewise-lienar atom
#' constructed by a minimum elementwise atom where the
#' objective function is the negative of variable t
#' t produced by max_elemwise_canon of the same dimension
#' as the expression and the constraints consist of a simple
#' inequality.
EliminatePwl.min_elemwise_canon <- function(expr, args) {
  tmp <- do.call(MaxElemwise, lapply(args, function(arg) { -arg }))
  canon <- EliminatePwl.max_elemwise_canon(tmp, tmp@args)
  return(list(-canon[[1]], canon[[2]]))
}

#' 
#' EliminatePwl canonicalizer for the 1 norm atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A canonicalization of the piecewise-lienar atom
#' constructed by the norm1 atom where the objective functino
#' consists of the sum of the variables created by the
#' abs_canon function and the constraints consist of
#' constraints generated by abs_canon.
EliminatePwl.norm1_canon <- function(expr, args) {
  x <- args[[1]]
  axis <- expr@axis
  
  # We need an absolute value constraint for the symmetric convex branches (p >= 1)
  constraints <- list()
  
  # TODO: Express this more naturally (recursively) in terms of the other atoms
  abs_expr <- abs(x)
  xconstr <- EliminatePwl.abs_canon(abs_expr, abs_expr@args)
  abs_x <- xconstr[[1]]
  abs_constraints <- xconstr[[2]]
  constraints <- c(constraints, abs_constraints)
  return(list(SumEntries(abs_x, axis = axis), constraints))
}

#' 
#' EliminatePwl canonicalizer for the infinite norm atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A canonicalization of the piecewise-lienar atom
#' constructed by the infinite norm atom where the objective
#' function consists variable t of the same dimension as the
#' expression and the constraints consist of a vector
#' constructed by multiplying t to a vector of 1's
EliminatePwl.norm_inf_canon <- function(expr, args) {
  x <- args[[1]]
  axis <- expr@axis
  # expr_dim <- dim(expr)
  # t <- Variable(expr_dim)
  t <- new("Variable", dim = dim(expr))
  
  if(is.na(axis))   # dim(expr) = c(1,1)
    promoted_t <- promote(t, dim(x))
  else if(axis == 2)   # dim(expr) = c(1,n)
    promoted_t <- Constant(matrix(1, nrow = nrow(x), ncol = 1) %*% reshape_expr(t, c(1, ncol(x))))
  else   # shape = c(m,1)
    promoted_t <- reshape_expr(t, c(nrow(x), 1)) %*% Constant(matrix(1, nrow = 1, ncol = ncol(x)))
  
  return(list(t, list(x <= promoted_t, x + promoted_t >= 0)))
}

#' 
#' EliminatePwl canonicalizer for the largest sum atom
#' 
#' @param expr An \linkS4class{Expression} object
#' @param args A list of \linkS4class{Constraint} objects
#' @return A canonicalization of the piecewise-lienar atom
#' constructed by the k largest sums atom where the objective
#' function consists of the sum of variables t that is of
#' the same dimension as the expression plus k
EliminatePwl.sum_largest_canon <- function(expr, args) {
  x <- args[[1]]
  k <- expr@k
  
  # min sum(t) + kq
  # s.t. x <= t + q, 0 <= t
  # t <- Variable(dim(x))
  t <- new("Variable", dim = dim(x))
  q <- Variable()
  obj <- sum(t) + k*q
  constraints <- list(x <= t + q, t >= 0)
  return(list(obj, constraints))
}

EliminatePwl.CANON_METHODS <- list(Abs = EliminatePwl.abs_canon,
                                   CumMax = EliminatePwl.cummax_canon,
                                   CumSum = EliminatePwl.cumsum_canon,
                                   MaxElemwise = EliminatePwl.max_elemwise_canon,
                                   MaxEntries = EliminatePwl.max_entries_canon,
                                   MinElemwise = EliminatePwl.min_elemwise_canon,
                                   MinEntries = EliminatePwl.min_entries_canon,
                                   Norm1 = EliminatePwl.norm1_canon,
                                   NormInf = EliminatePwl.norm_inf_canon,
                                   SumLargest = EliminatePwl.sum_largest_canon)
