# Scalar functions
geo_mean <- GeoMean
harmonic_mean <- HarmonicMean
lambda_max <- LambdaMax
lambda_min <- LambdaMin
lambda_sum_largest <- LambdaSumLargest
lambda_sum_smallest <- LambdaSumSmallest
log_det <- LogDet
log_sum_exp <- LogSumExp
matrix_frac <- MatrixFrac
max_entries <- MaxEntries
min_entries <- MinEntries
mixed_norm <- MixedNorm
norm1 <- Norm1
norm2 <- Norm2
norm_inf <- NormInf
norm_nuc <- NormNuc
pnorm <- Pnorm
quad_form <- QuadForm
quad_over_lin <- QuadOverLin
sum_entries <- SumEntries
sum_largest <- SumLargest
sum_smallest <- SumSmallest
sum_squares <- SumSquares
matrix_trace <- Trace
tv <- TotalVariation

max.Expression <- function(..., na.rm = FALSE) {
  if(na.rm)
    warning("na.rm is unimplemented for Expression objects")
  
  vals <- list(...)
  is_expr <- sapply(vals, function(v) { is(v, "Expression") })
  max_args <- lapply(vals[is_expr], function(expr) { MaxEntries(expr) })
  if(!all(is_expr)) {
    max_num <- max(sapply(vals[!is_expr], function(v) { max(v, na.rm = na.rm) }))
    max_args <- c(max_args, max_num)
  }
  .MaxElemwise(args = max_args)
}

min.Expression <- function(..., na.rm = FALSE) {
  if(na.rm)
    warning("na.rm is unimplemented for Expression objects")
  
  vals <- list(...)
  is_expr <- sapply(vals, function(v) { is(v, "Expression") })
  min_args <- lapply(vals[is_expr], function(expr) { MinEntries(expr) })
  if(!all(is_expr)) {
    min_num <- min(sapply(vals[!is_expr], function(v) { min(v, na.rm = na.rm) }))
    min_args <- c(min_args, min_num)
  }
  min_args <- lapply(min_args, function(arg) { -as.Constant(arg) })
  -.MaxElemwise(args = min_args)
}

setMethod("norm", signature(x = "Expression", type = "character"), function(x, type) {
  x <- as.Constant(x)
  type <- substr(type, 1, 1)
  
  # Norms for scalars same as absolute value
  if(type %in% c("O", "o", "1"))                  # Maximum absolute column sum
    MaxEntries(Pnorm(x = x, p = 1, axis = 2))
  else if(type %in% c("I", "i"))                  # Maximum absolute row sum
    MaxEntries(Pnorm(x = x, p = 1, axis = 1))
  else if(type %in% c("E", "e", "F", "f"))        # Frobenius norm (Euclidean norm if x is treated as a vector)
    Pnorm(x = x, p = 2, axis = NA_real_)
  else if(type %in% c("M", "m"))                  # Maximum modulus (absolute value) of all elements in x
    MaxEntries(Abs(x = x))
  else if(type == "2")                            # Spectral norm (largest singular value of x)
    SigmaMax(A = x)
  else
    stop("argument type[1]='", type, "' must be one of 'M','1','O','I','F' or 'E'")
})

sum.Expression <- function(..., na.rm = FALSE) {
  if(na.rm)
    warning("na.rm is unimplemented for Expression objects")
  
  vals <- list(...)
  is_expr <- sapply(vals, function(v) { is(v, "Expression") })
  sum_expr <- lapply(vals[is_expr], function(expr) { SumEntries(expr = expr) })
  if(all(is_expr))
    Reduce("+", sum_expr)
  else {
    sum_num <- sum(sapply(vals[!is_expr], function(v) { sum(v, na.rm = na.rm) }))
    Reduce("+", sum_expr) + sum_num
  }
}

mean.Expression <- function(x, trim = 0, na.rm = FALSE, ...) {
  if(na.rm)
    stop("na.rm is unimplemented for Expression objects")
  if(trim != 0)
    stop("trim is unimplemented for Expression objects")
  SumEntries(expr = x) / prod(size(x))
}

# Elementwise functions
entr <- Entr
huber <- Huber
inv_pos <- InvPos
kl_div <- KLDiv
logistic <- Logistic
max_elemwise <- MaxElemwise
min_elemwise <- MinElemwise
mul_elemwise <- MulElemwise
neg <- Neg
pos <- Pos
power <- Power
scalene <- Scalene
square <- Square

setMethod("abs", "Expression", function(x) { Abs(x = x) })
setMethod("exp", "Expression", function(x) { Exp(x = x) })
setMethod("log", "Expression", function(x, base = exp(1)) { Log(x = x)/log(base) })
setMethod("log10", "Expression", function(x) { log(x, base = 10) })
setMethod("log2", "Expression", function(x) { log(x, base = 2) })
setMethod("log1p", "Expression", function(x) { Log1p(x = x) })
setMethod("sqrt", "Expression", function(x) { Sqrt(x = x) })

# Matrix/vector operations
bmat <- Bmat
conv <- Conv
hstack <- HStack
kron <- Kron
reshape_expr <- Reshape
sigma_max <- SigmaMax
vec <- Vec
vstack <- VStack

setMethod("cumsum", "Expression", function(x) { CumSum(expr = Vec(x)) })   # Flatten matrix in column-major order to match R's behavior
setMethod("diag", "Expression", function(x, nrow, ncol) {
  if(nargs() == 1L)
    Diag(x)
  else if(is_matrix(x))
    stop("'nrow' or 'ncol' cannot be specified when 'x' is a matrix")
  else {
    expr <- as.Constant(x)
    n <- length(expr)
    if(!missing(nrow))
      n <- nrow
    if(missing(ncol))
      ncol <- n
    expr*(base::diag(n)[1:n, 1:ncol])
  }
})
setMethod("diff", "Expression", function(x, lag = 1, differences = 1, ...) { Diff(x = x, lag = lag, k = differences, ...) })

setMethod("kronecker", signature(X = "Expression", Y = "ANY"), function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
  if(FUN != "*" || make.dimnames)
    stop("Unimplemented")
  Kron(X, Y)
})
setMethod("kronecker", signature(X = "ANY", Y = "Expression"), function(X, Y, FUN = "*", make.dimnames = FALSE, ...) {
  if(FUN != "*" || make.dimnames)
    stop("Unimplemented")
  Kron(X, Y)
})
setMethod("%x%", signature(X = "Expression", Y = "ANY"), function(X, Y) { Kron(lh_exp = X, rh_exp = Y) })
setMethod("%x%", signature(X = "ANY", Y = "Expression"), function(X, Y) { Kron(lh_exp = X, rh_exp = Y) })

setMethod("cbind2", signature(x = "Expression", y = "ANY"), function(x, y, ...) { HStack(x, y) })
setMethod("cbind2", signature(x = "ANY", y = "Expression"), function(x, y, ...) { HStack(x, y) })
setMethod("rbind2", signature(x = "Expression", y = "ANY"), function(x, y, ...) { VStack(x, y) })
setMethod("rbind2", signature(x = "ANY", y = "Expression"), function(x, y, ...) { VStack(x, y) })
