context("test-g01-atoms")
TOL <- 1e-5

CONSTANT <- "CONSTANT"
AFFINE <- "AFFINE"
CONVEX <- "CONVEX"
CONCAVE <- "CONCAVE"
ZERO <- "ZERO"
UNKNOWN <- "UNKNOWN"
NONNEG <- "NONNEGATIVE"
NONPOS <- "NONPOSITIVE"

a <- Variable(name = "a")

x <- Variable(2, name = "x")
y <- Variable(2, name = "y")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("test the norm_inf function", {
  exp <- x + y
  atom <- norm_inf(exp)

  # expect_equal(dim(atom), NULL)
  expect_equal(dim(atom), c(1,1))
  expect_equal(curvature(atom), CONVEX)
  expect_true(is_convex(atom))
  expect_true(is_concave(-atom))
  expect_equal(curvature(norm_inf(atom)), CONVEX)
  expect_equal(curvature(norm_inf(-atom)), CONVEX)
})

test_that("test the norm1 function", {
  exp <- x + y
  atom <- norm1(exp)

  # expect_equal(dim(atom), NULL)
  expect_equal(dim(atom), c(1,1))
  expect_equal(curvature(atom), CONVEX)
  expect_equal(curvature(norm1(atom)), CONVEX)
  expect_equal(curvature(norm1(-atom)), CONVEX)
})

test_that("test quad_form function", {
  P <- Parameter(2,2)
  expr <- quad_form(x, P)
  expect_false(is_dcp(expr))
})

test_that("test the power function", {
  for(a_dim in list(c(1, 1), c(3, 1), c(2, 3))) {
    # x_pow <- Variable(a_dim)
    # y_pow <- Variable(a_dim)
    x_pow <- new("Variable", dim = a_dim)
    y_pow <- new("Variable", dim = a_dim)
    exp <- x_pow + y_pow

    for(p in c(0, 1, 2, 3, 2.7, 0.67, -1, -2.3, 4/5)) {
      atom <- power(exp, p)
      expect_equal(dim(atom), a_dim)

      if(p > 1 || p < 0)
        expect_equal(curvature(atom), CONVEX)
      else if(p == 1)
        expect_equal(curvature(atom), AFFINE)
      else if(p == 0)
        expect_equal(curvature(atom), CONSTANT)
      else
        expect_equal(curvature(atom), CONCAVE)

      if(p != 1)
        expect_equal(sign(atom), NONNEG)
    }
  }
  expect_error(value(power(-1, 3)),
               "Power cannot be applied to negative values", fixed = TRUE)

  expect_error(value(power(0, -1)),
               "Power cannot be applied to negative or zero values", fixed = TRUE)

})

test_that("test the geo_mean function", {
  atom <- geo_mean(x)
  # expect_equal(dim(atom), NULL)
  expect_equal(dim(atom), c(1,1))
  expect_equal(curvature(atom), CONCAVE)
  expect_equal(sign(atom), NONNEG)
})

test_that("test the harmonic_mean function", {
  atom <- harmonic_mean(x)
  # expect_equal(dim(atom), NULL)
  expect_equal(dim(atom), c(1,1))
  expect_equal(curvature(atom), CONCAVE)
  expect_equal(sign(atom), NONNEG)
})

test_that("test the p_norm function", {
  atom <- p_norm(x, p = 1.5)
  # expect_equal(dim(atom), NULL)
  expect_equal(dim(atom), c(1,1))
  expect_equal(curvature(atom), CONVEX)
  expect_equal(sign(atom), NONNEG)

  atom <- p_norm(x, p = 1)
  # expect_equal(dim(atom), NULL)
  expect_equal(dim(atom), c(1,1))
  expect_equal(curvature(atom), CONVEX)
  expect_equal(sign(atom), NONNEG)

  atom <- p_norm(x, p = 2)
  # expect_equal(dim(atom), NULL)
  expect_equal(dim(atom), c(1,1))
  expect_equal(curvature(atom), CONVEX)
  expect_equal(sign(atom), NONNEG)

  atom <- p_norm(x, p = Inf)
  # expect_equal(dim(atom), NULL)
  expect_equal(dim(atom), c(1,1))
  expect_equal(curvature(atom), CONVEX)
  expect_equal(sign(atom), NONNEG)

  atom <- p_norm(x, p = 0.5)
  # expect_equal(dim(atom), NULL)
  expect_equal(dim(atom), c(1,1))
  expect_equal(curvature(atom), CONCAVE)
  expect_equal(sign(atom), NONNEG)

  atom <- p_norm(x, p = 0.7)
  # expect_equal(dim(atom), NULL)
  expect_equal(dim(atom), c(1,1))
  expect_equal(curvature(atom), CONCAVE)
  expect_equal(sign(atom), NONNEG)

  atom <- p_norm(x, p = -0.1)
  # expect_equal(dim(atom), NULL)
  expect_equal(dim(atom), c(1,1))
  expect_equal(curvature(atom), CONCAVE)
  expect_equal(sign(atom), NONNEG)

  atom <- p_norm(x, p = -1)
  # expect_equal(dim(atom), NULL)
  expect_equal(dim(atom), c(1,1))
  expect_equal(curvature(atom), CONCAVE)
  expect_equal(sign(atom), NONNEG)

  atom <- p_norm(x, p = -1.3)
  # expect_equal(dim(atom), NULL)
  expect_equal(dim(atom), c(1,1))
  expect_equal(curvature(atom), CONCAVE)
  expect_equal(sign(atom), NONNEG)
})

test_that("test matrix norms", {
  for(p in c("1", "2", "I", "F", "M")) {
    for(var in c(A, C)) {
      atom <- norm(var, p)
      # expect_equal(dim(atom), NULL)
      expect_equal(dim(atom), c(1,1))
      expect_equal(curvature(atom), CONVEX)
      expect_equal(sign(atom), NONNEG)

      value(var) <- matrix(rnorm(size(var)), nrow = nrow(var), ncol = ncol(var))
      expect_equal(value(norm(var, p)), base:::norm(value(var), type = p))
    }
  }
})

test_that("test the quad_over_lin function", {
  atom <- quad_over_lin(x^2, a)
  expect_equal(curvature(atom), CONVEX)

  atom <- quad_over_lin(-x^2, a)
  expect_equal(curvature(atom), CONVEX)

  atom <- quad_over_lin(sqrt(x), a)
  expect_equal(curvature(atom), UNKNOWN)
  expect_false(is_dcp(atom))

  expect_error(quad_over_lin(x, x),
               "The second argument to QuadOverLin must be a scalar.", fixed = TRUE)
})

test_that("test the arg count for max_elemwise and min_elemwise", {
  expect_error(max_elemwise(1),
               "argument \"arg2\" is missing, with no default", fixed = TRUE)
  expect_error(min_elemwise(1),
               "argument \"arg2\" is missing, with no default", fixed = TRUE)
})

test_that("test the matrix_frac function", {
  atom <- matrix_frac(x, A)
  # expect_equal(dim(atom), NULL)
  expect_equal(dim(atom), c(1,1))
  expect_equal(curvature(atom), CONVEX)

  # Test matrix_frac dim validation
  expect_error(matrix_frac(x, C),
               "The second argument to MatrixFrac must be a square matrix.", fixed = TRUE)
  expect_error(matrix_frac(Variable(3), A),
               "The arguments to MatrixFrac have incompatible dimensions.", fixed = TRUE)
})

test_that("test the sign for max_entries", {
  expect_equal(sign(max_entries(1)), NONNEG)
  expect_equal(sign(max_entries(-2)), NONPOS)
  expect_equal(sign(max_entries(Variable())), UNKNOWN)
  expect_equal(sign(max_entries(0)), ZERO)

  # Test with axis argument
  expect_equal(dim(max_entries(Variable(2), axis = 1)), c(2, 1))
  expect_equal(dim(max_entries(Variable(2), axis = 2, keepdims = TRUE)), c(1, 1))
  # expect_equal(dim(max_entries(Variable(c(2, 3)), axis = 1)), 2)
  # expect_equal(dim(max_entries(Variable(c(2, 3)), axis = 2, keepdims = TRUE)), c(1, 3))
  expect_equal(dim(max_entries(Variable(2, 3), axis = 1)), c(2, 1))
  expect_equal(dim(max_entries(Variable(2, 3), axis = 2, keepdims = TRUE)), c(1, 3))

  # Invalid axis
  expect_error(max_entries(x, axis = 4),
               "Invalid argument for axis. Must be an integer between 1 and 2", fixed = TRUE)
})

test_that("test the sign for min_entries", {
  expect_equal(sign(min_entries(1)), NONNEG)
  expect_equal(sign(min_entries(-2)), NONPOS)
  expect_equal(sign(min_entries(Variable())), UNKNOWN)
  expect_equal(sign(min_entries(0)), ZERO)

  # Test with axis argument
  expect_equal(dim(min_entries(Variable(2), axis = 1)), c(2, 1))
  expect_equal(dim(min_entries(Variable(2), axis = 2, keepdims = TRUE)), c(1, 1))
  # expect_equal(dim(min_entries(Variable(c(2, 3)), axis = 1)), 2)
  # expect_equal(dim(min_entries(Variable(c(2, 3)), axis = 2, keepdims = TRUE)), c(1, 3))
  expect_equal(dim(min_entries(Variable(2, 3), axis = 1)), c(2, 1))
  expect_equal(dim(min_entries(Variable(2, 3), axis = 2, keepdims = TRUE)), c(1, 3))

  # Invalid axis
  expect_error(min_entries(x, axis = 4),
               "Invalid argument for axis. Must be an integer between 1 and 2", fixed = TRUE)
})

test_that("test sign logic for max_elemwise", {
  expect_equal(sign(max_elemwise(1, 2)), NONNEG)
  expect_equal(sign(max_elemwise(1, Variable())), NONNEG)
  expect_equal(sign(max_elemwise(1, -2)), NONNEG)
  expect_equal(sign(max_elemwise(1, 0)), NONNEG)

  expect_equal(sign(max_elemwise(Variable(), 0)), NONNEG)
  expect_equal(sign(max_elemwise(Variable(), Variable())), UNKNOWN)
  expect_equal(sign(max_elemwise(Variable(), -2)), UNKNOWN)

  expect_equal(sign(max_elemwise(0, 0)), ZERO)
  expect_equal(sign(max_elemwise(0, -2)), ZERO)

  expect_equal(sign(max_elemwise(-3, -2)), NONPOS)

  # Many args
  expect_equal(sign(max_elemwise(-2, Variable(), 0, -1, Variable(), -1)), NONNEG)

  # Promotion
  expect_equal(sign(max_elemwise(1, Variable(2))), NONNEG)
  expect_equal(dim(max_elemwise(1, Variable(2))), c(2, 1))
})

test_that("test sign logic for min_elemwise", {
  expect_equal(sign(min_elemwise(1, 2)), NONNEG)
  expect_equal(sign(min_elemwise(1, Variable())), UNKNOWN)
  expect_equal(sign(min_elemwise(1, -2)), NONPOS)
  expect_equal(sign(min_elemwise(1, 0)), ZERO)

  expect_equal(sign(min_elemwise(Variable(), 0)), NONPOS)
  expect_equal(sign(min_elemwise(Variable(), Variable())), UNKNOWN)
  expect_equal(sign(min_elemwise(Variable(), -2)), NONPOS)

  expect_equal(sign(min_elemwise(0, 0)), ZERO)
  expect_equal(sign(min_elemwise(0, -2)), NONPOS)

  expect_equal(sign(min_elemwise(-3, -2)), NONPOS)

  # Many args
  expect_equal(sign(min_elemwise(-2, Variable(), 0, -1, Variable(), 1)), NONPOS)

  # Promotion
  expect_equal(sign(min_elemwise(-1, Variable(2))), NONPOS)
  expect_equal(dim(min_elemwise(-1, Variable(2))), c(2, 1))
})

test_that("test the sum_entries function", {
  expect_equal(sign(sum_entries(1)), NONNEG)
  expect_equal(sign(sum_entries(c(1, -1))), UNKNOWN)
  expect_equal(curvature(sum_entries(c(1, -1))), CONSTANT)
  expect_equal(sign(sum_entries(Variable(2))), UNKNOWN)
  # expect_equal(dim(sum_entries(Variable(2))), NULL)
  expect_equal(dim(sum_entries(Variable(2))), c(1, 1))
  expect_equal(curvature(sum_entries(Variable(2))), AFFINE)
  # expect_equal(dim(sum_entries(Variable(c(2, 1)), keepdims = TRUE)), c(1, 1))
  expect_equal(dim(sum_entries(Variable(2, 1), keepdims = TRUE)), c(1, 1))

  # Mixed curvature
  mat <- matrix(c(1,-1), nrow = 1, ncol = 2)
  expect_equal(curvature(sum_entries( mat %*% Variable(2)^2 )), UNKNOWN)

  # Test with axis argument
  expect_equal(dim(sum_entries(Variable(2), axis = 1)), c(2, 1))
  # expect_equal(dim(sum_entries(Variable(2), axis = 2)), NULL)
  expect_equal(dim(sum_entries(Variable(2), axis = 2)), c(1, 1))
  # expect_equal(dim(sum_entries(Variable(c(2, 3)), axis = 1)), 2)
  # expect_equal(dim(sum_entries(Variable(c(2, 3)), axis = 2, keepdims = TRUE)), c(1, 3))
  # expect_equal(dim(sum_entries(Variable(c(2, 3)), axis = 2, keepdims = FALSE)), 3)
  expect_equal(dim(sum_entries(Variable(2, 3), axis = 1)), c(2, 1))
  expect_equal(dim(sum_entries(Variable(2, 3), axis = 2, keepdims = TRUE)), c(1, 3))
  expect_equal(dim(sum_entries(Variable(2, 3), axis = 2, keepdims = FALSE)), c(3, 1))

  # Invalid axis
  expect_error(sum_entries(x, axis = 4),
               "Invalid argument for axis. Must be an integer between 1 and 2", fixed = TRUE)

  A <- diag(3)
  expect_equal(value(CVXR::sum_entries(A)), 3)

  A <- diag(3)
  expect_equal(value(CVXR::sum_entries(A, axis = 1)), c(1, 1, 1))

})

test_that("test the multiply function", {
  expect_equal(sign(multiply(c(1, -1), x)), UNKNOWN)
  expect_equal(curvature(multiply(c(1, -1), x)), AFFINE)
  expect_equal(dim(multiply(c(1, -1), x)), c(2, 1))
  pos_param <- Parameter(2, nonneg = TRUE)
  neg_param <- Parameter(2, nonpos = TRUE)
  expect_equal(sign(multiply(pos_param, pos_param)), NONNEG)
  expect_equal(sign(multiply(pos_param, neg_param)), NONPOS)
  expect_equal(sign(multiply(neg_param, neg_param)), NONNEG)

  expect_equal(curvature(multiply(neg_param, x^2)), CONCAVE)

  # Test promotion
  expect_equal(dim(multiply(c(1, -1), 1)), c(2, 1))
  expect_equal(dim(multiply(1, C)), dim(C))
  # expect_error(multiply(x, c(1, -1)))

  expect_equal(sign(multiply(x, c(1, -1))), UNKNOWN)
  expect_equal(curvature(multiply(x, c(1, -1))), AFFINE)
  expect_equal(dim(multiply(x, c(1, -1))), c(2,1))
})

test_that("test the vstack function", {
  atom <- vstack(x, y, x)
  expect_equal(dim(atom), c(6, 1))

  atom <- vstack(A, C, B)
  expect_equal(dim(atom), c(7, 2))

  entries <- list()
  for(i in 1:dim(x)[1]) {
   for(j in 1:dim(x)[2]) {
     entries <- c(entries, x[i, j])
   }
  }
  atom <- do.call(vstack, entries)
  # atom <- vstack(x[1,1], x[2,1])

  expect_error(vstack(C, 1),
               "All the input dimensions except for axis 1 must match exactly.", fixed = TRUE)

})

test_that("test the reshape_expr function", {
  expr <- reshape_expr(A, c(4, 1))
  expect_equal(sign(expr), UNKNOWN)
  expect_equal(curvature(expr), AFFINE)
  expect_equal(dim(expr), c(4, 1))

  expr <- reshape_expr(expr, c(2, 2))
  expect_equal(dim(expr), c(2, 2))

  expr <- reshape_expr(x^2, c(1, 2))
  expect_equal(sign(expr), NONNEG)
  expect_equal(curvature(expr), CONVEX)
  expect_equal(dim(expr), c(1, 2))

  expect_error(reshape_expr(C, c(5, 4)),
               "Invalid reshape dimensions (54)", fixed = TRUE)
})

test_that("test the vec function", {
  expr <- vec(C)
  expect_equal(sign(expr), UNKNOWN)
  expect_equal(curvature(expr), AFFINE)
  expect_equal(dim(expr), c(6, 1))

  expr <- vec(x)
  expect_equal(dim(expr), c(2, 1))

  expr <- vec(a^2)
  expect_equal(sign(expr), NONNEG)
  expect_equal(curvature(expr), CONVEX)
  expect_equal(dim(expr), c(1, 1))
})

test_that("test the diag function", {
  expr <- diag(x)
  expect_equal(sign(expr), UNKNOWN)
  expect_equal(curvature(expr), AFFINE)
  expect_equal(dim(expr), c(2,2))

  expr <- diag(A)
  expect_equal(sign(expr), UNKNOWN)
  expect_equal(curvature(expr), AFFINE)
  expect_equal(dim(expr), c(2, 1))

  expr <- diag(t(x))
  expect_equal(sign(expr), UNKNOWN)
  expect_equal(curvature(expr), AFFINE)
  expect_equal(dim(expr), c(2, 2))

  expect_error(diag(C),
               "Argument to Diag must be a vector or square matrix.", fixed = TRUE)
})

test_that("test the matrix_trace function", {
  expr <- matrix_trace(A)
  expect_equal(sign(expr), UNKNOWN)
  expect_equal(curvature(expr), AFFINE)
  # expect_equal(dim(expr), NULL)
  expect_equal(dim(expr), c(1, 1))

  expect_error(matrix_trace(C),
               "Argument to Trace must be a square matrix", fixed = TRUE)
})

test_that("test the log1p function", {
  expr <- log1p(Constant(1))
  expect_equal(sign(expr), NONNEG)
  expect_equal(curvature(expr), CONSTANT)
  expect_equal(dim(expr), c(1, 1))
  expr <- CVXR::log1p(Constant(-0.5))
  expect_equal(sign(expr), NONPOS)
})

test_that("test the upper_tri function", {
  expect_error(upper_tri(C),
               "Argument to UpperTri must be a square matrix.", fixed = TRUE)
})

test_that("test the huber function", {
  huber(x, 1)
  expect_error(huber(x, -1),
               "M must be a non-negative scalar constant", fixed = TRUE)
  expect_error(huber(x, c(1, 1)),
               "M must be a non-negative scalar constant", fixed = TRUE)

  # M parameter
  M <- Parameter(nonneg = TRUE)
  # Valid
  huber(x, M)
  value(M) <- 1
  expect_equal(value(huber(2, M)), 3, tolerance = TOL)
  # Invalid
  M <- Parameter(nonpos = TRUE)
  expect_error(huber(x, M),
               "M must be a non-negative scalar constant", fixed = TRUE)
})

test_that("test the sum_largest function", {
  expect_error(sum_largest(x, -1),
               "[SumLargest: validation] k must be a positive integer", fixed = TRUE)
  expect_error(lambda_sum_largest(x, 2.4),
               "First argument must be a square matrix.", fixed = TRUE)
  expect_error(lambda_sum_largest(Variable(2, 2), 2.4),
               "Second argument must be a positive integer.", fixed = TRUE)
})

test_that("test the sum_smallest function", {
  expect_error(sum_smallest(x, -1),
               "[SumLargest: validation] k must be a positive integer", fixed = TRUE)
  expect_error(lambda_sum_smallest(Variable(2, 2), 2.4),
               "Second argument must be a positive integer.", fixed = TRUE)
})

test_that("test the bmat function", {
  v_np <- matrix(1, nrow = 3, ncol = 1)
  v_00 <- matrix(c(0,0), nrow = 2, ncol = 1)
  v_12 <- matrix(c(1,2), nrow = 2, ncol = 1)
  expr <- bmat(list(list(v_np, v_np), list(v_00, v_12)))
  expect_equal(dim(expr), c(5, 2))
  const <- rbind(cbind(v_np, v_np), cbind(c(0, 0), c(1, 2)))
  expect_equal(value(expr), const)
})

test_that("test the conv function", {
  a <- matrix(1, nrow = 3, ncol = 1)
  b <- Parameter(2, nonneg = TRUE)
  expr <- conv(a, b)
  expect_true(is_nonneg(expr))
  expect_equal(dim(expr), c(4, 1))
  b <- Parameter(2, nonpos = TRUE)
  expr <- conv(a, b)
  expect_true(is_nonpos(expr))
  expect_error(conv(x, -1),
               "The first argument to Conv must be constant.", fixed = TRUE)
  expect_error(conv(cbind(c(0, 1), c(0, 1)), x),
               "The arguments to Conv must resolve to vectors.", fixed = TRUE)
})

test_that("test the kronecker function", {
  a <- matrix(1, nrow = 3, ncol = 2)
  b <- Parameter(2, nonneg = TRUE)
  expr <- kronecker(a, b)
  expect_true(is_nonneg(expr))
  expect_equal(dim(expr), c(6, 2))
  b <- Parameter(2, nonpos = TRUE)
  expr <- kronecker(a, b)
  expect_true(is_nonpos(expr))
  expect_error(kronecker(x, -1),
               "The first argument to Kron must be constant.", fixed = TRUE)
})

# test_that("test DCP properties of partial optimize", {
#   # Evaluate the 1-norm in the usual way (i.e., in epigraph form)
#   dims <- 3
#   x <- Variable(dims)
#   t <- Variable(dims)
#   xval <- matrix(rep(-5, dims), nrow = dims, ncol = 1)
#   p2 <- Problem(Minimize(sum_entries(t)), list(-t <= x, x <= t))
#   g <- partial_optimize(p2, list(t), list(x))
#   expect_equal(curvature(g), CONVEX)
#
#   p2 <- Problem(Maximize(sum_entries(t)), list(-t <= x, x <= t))
#   g <- partial_optimize(p2, list(t), list(x))
#   expect_equal(curvature(g), CONCAVE)
#
#   p2 <- Problem(Maximize(t[1]^2), list(-t <= x, x <= t))
#   g <- partial_optimize(p2, list(t), list(x))
#   expect_false(is_convex(g))
#   expect_false(is_concave(g))
# })
#
# test_that("test the partial_optimize eval 1-norm", {
#   # Evaluate the 1-norm in the usual way (i.e., in epigraph form)
#   dims <- 3
#   x <- Variable(dims)
#   t <- Variable(dims)
#   xval <- matrix(rep(-5, dims), nrow = dims, ncol = 1)
#   p1 <- Problem(Minimize(sum_entries(t)), list(-t <= xval, xval <= t))
#   result1 <- solve(p1)
#
#   # Minimize the 1-norm via partial_optimize
#   p2 <- Problem(Minimize(sum_entries(t)), list(-t <= x, x <= t))
#   # g <- partial_optimize(p2, list(t), list(x))
#   # p3 <- Problem(Minimize(g), list(x == xval))
#   # result3 <- solve(p3)
#   # expect_equal(result1$value, -result3$value)
#
#   # Try leaving out args
#
#   # Minimize the 1-norm via partial_optimize
#   p2 <- Problem(Minimize(sum_entries(t)), list(-t <= x, x <= t))
#   # g <- partial_optimize(p2, opt_vars = list(t))
#   # p3 <- Problem(Minimize(g), list(x == xval))
#   # result3 <- solve(p3)
#   # expect_equal(result1$value, result3$value)
#
#   # Minimize the 1-norm via partial_optimize
#   # g <- partial_optimize(p2, dont_opt_vars = list(x))
#   # p3 <- Problem(Minimize(g), list(x == xval))
#   # result3 <- solve(p3)
#   # expect_equal(result1$value, result3$value)
#
#   # expect_error(partial_optimize(p2))
#   # expect_error(partial_optimize(p2, list(), list(x)))
# })
#
# test_that("test partial_optimize min 1-norm", {
#   # Minimize the 1-norm in the usual way
#   dims <- 3
#   x <- Variable(dims)
#   t <- Variable(dims)
#   p1 <- Problem(Minimize(sum_entries(t)), list(-t <= x, x <= t))
#
#   # Minimize the 1-norm via partial_optimize
#   # g <- partial_optimize(p1, list(t), list(x))
#   # p2 <- Problem(Minimize(g))
#   # result2 <- solve(p2)
#
#   result1 <- solve(p1)
#   # expect_equal(result1$value, result2$value)
# })
#
# test_that("test partial_optimize simple problem", {
#   x <- Variable(1)
#   y <- Variable(1)
#
#   # Solve the (simple) two-stage problem by "combining" the two stages (i.e., by solving a single linear program)
#   p1 <- Problem(Minimize(x+y), list(x+y >= 3, y >= 4, x >= 5))
#   result1 <- solve(p1)
#
#   # Solve the two-stage problem via partial_optimize
#   p2 <- Problem(Minimize(y), list(x+y >= 3, y >= 4))
#   # g <- partial_optimize(p2, list(y), list(x))
#   # p3 <- Problem(Minimize(x+g), list(x >= 5))
#   # result3 <- solve(p3)
#   # expect_equal(result1$value, result3$value)
# })
#
# test_that("test partial_optimize special var", {
#   x <- Bool(1)
#   y <- Int(1)
#
#   # Solve the (simple) two-stage problem by "combining" the two stages (i.e., by solving a single linear program)
#   p1 <- Problem(Minimize(x+y), list(x+y >= 3, y >= 4, x >= 5))
#   # result1 <- solve(p1)
#
#   # Solve the two-stage problem via partial_optimize
#   p2 <- Problem(Minimize(y), list(x+y >= 3, y >= 4))
#   # g <- partial_optimize(p2, list(y), list(x))
#   # p3 <- Problem(Minimize(x+g), list(x >= 5))
#   # result3 <- solve(p3)
#   # expect_equal(result1$value, result3$value)
# })
#
# test_that("test partial_optimize special constr", {
#   x <- Variable(1)
#   y <- Variable(1)
#
#   # Solve the (simple) two-stage problem by "combining" the two stages (i.e., by solving a single linear program)
#   p1 <- Problem(Minimize(x+exp(y)), list(x+y >= 3, y >= 4, x >= 5))
#   result1 <- solve(p1)
#
#   # Solve the two-stage problem via partial_optimize
#   p2 <- Problem(Minimize(exp(y)), list(x+y >= 3, y >= 4))
#   # g <- partial_optimize(p2, list(y), list(x))
#   # p3 <- Problem(Minimize(x+g), list(x >= 5))
#   # result3 <- solve(p3)
#   # expect_equal(result1$value, result3$value)
# })
#
# test_that("test partial_optimize with parameters", {
#   x <- Variable(1)
#   y <- Variable(1)
#   gamma <- Parameter()
#
#   # Solve the (simple) two-stage problem by "combining" the two stages (i.e., by solving a single linear program)
#   p1 <- Problem(Minimize(x+y), list(x+y >= gamma, y >= 4, x >= 5))
#   gamma@value <- 3
#   # result1 <- solve(p1)
#
#   # Solve the two-stage problem via partial_optimize
#   p2 <- Problem(Minimize(y), list(x+y >= gamma, y >= 4))
#   # g <- partial_optimize(p2, list(y), list(x))
#   # p3 <- Problem(Minimize(x+g), list(x >= 5))
#   # result3 <- solve(p3)
#   # expect_equal(result1$value, result3$value)
# })
#
# test_that("test partial_optimize numeric function", {
#   x <- Variable(1)
#   y <- Variable(1)
#   xval <- 4
#
#   # Solve the (simple) two-stage problem by "combining" the two stages (i.e., by solving a single linear program)
#   p1 <- Problem(Minimize(y), list(xval+y >= 3))
#   result1 <- solve(p1)
#
#   # Solve the two-stage problem via partial_optimize
#   constr <- list(y >= -100)
#   p2 <- Problem(Minimize(y), c(x+y >= 3, constr))
#   # g <- partial_optimize(p2, list(y), list(x))
#   # x@value <- xval
#   # y@value <- 42
#   # const[1]@dual_variable@value <- 42
#   # result <- g@value
#   # expect_equal(result, result1$value)
#   # expect_equal(y@value, 42)
#   # expect_equal(constr[1]@dual_value, 42)
#
#   # No variables optimized over
#   p2 <- Problem(Minimize(y), list(x+y >= 3))
#   # g <- partial_optimize(p2, list(), list(x,y))
#   # x@value <- xval
#   # y@value <- 42
#   # p2@constraints[1]@dual_variable@value <- 42
#   # result <- g@value
#   # expect_equal(result, y@value)
#   # expect_equal(y@value, 42)
#   # expect_equal(p2@constraints[1]@dual_value, 42)
# })
#
# test_that("test partial_optimize stacked", {
#   # Minimize the 1-norm in the usual way
#   dims <- 3
#   x <- Variable(dims)
#   t <- Variable(dims)
#   p1 <- Problem(Minimize(sum_entries(t)), list(-t <= x, x <= t))
#
#   # Minimize the 1-norm via partial_optimize
#   # g <- partial_optimize(p1, list(t), list(x))
#   # g2 <- partial_optimize(Problem(Minimize(g)), list(x))
#   # p2 <- Problem(Minimize(g2))
#   # result2 <- solve(p2)
#
#   result1 <- solve(p1)
#   # expect_equal(result1$value, result2$value)
# })

test_that("test the NonNegative Variable class", {
  # x <- NonNegative()
  x <- Variable(nonneg = TRUE)
  p <- Problem(Minimize(5+x), list(x >= 3))
  result <- solve(p)
  expect_equal(result$value, 8, tolerance = TOL)
  expect_equal(result$getValue(x), 3, tolerance = TOL)
})

test_that("test mixed_norm", {
  y <- Variable(5, 5)
  obj <- Minimize(mixed_norm(y, Inf, 1))
  prob <- Problem(obj, list(y == matrix(1, nrow = 5, ncol = 5)))
  result <- solve(prob)
  expect_equal(result$value, 5, tolerance = TOL)
})

# test_that("test whether changing an array constant breaks DCP", {
#   c <- matrix(c(1, 2), nrow = 2, ncol = 1)
#   x@value <- c(1, 1)
#   expr <- t(c) %*% x^2
#   expect_equal(value(expr), 3, tolerance = TOL)
#   expect_true(is_dcp(expr))
#
#   c[1] <- -1
#   expect_equal(value(expr), 3, tolerance = TOL)
#   expect_true(is_dcp(expr))
# })

test_that("test that norm1 and normInf match definition for matrices", {
  A <- rbind(c(1,2), c(3,4))
  print(A)
  X <- Variable(2, 2)
  obj <- Minimize(norm1(X))
  prob <- Problem(obj, list(X == A))
  result <- solve(prob)
  print(result$value)
  expect_equal(result$value, value(norm1(A)), TOL)

  obj <- Minimize(norm_inf(X))
  prob <- Problem(obj, list(X == A))
  result <- solve(prob)
  print(result$value)
  expect_equal(result$value, value(norm_inf(A)), TOL)
})

#DK, uncomment once Indicator in transforms is uncommented
# test_that("test indicator", {
#   x <- Variable()
#   constraints <- list(0 <= x, x <= 1)
#   expr <- Indicator(constraints)
#   x@value <- .5
#   expect_equal(value(expr), 0.0)
#   x@value <- 2
#   expect_equal(value(expr), Inf)
#
# })
