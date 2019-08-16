context("test-g01-expressions")
TOL <- 1e-6

CONSTANT <- "CONSTANT"
AFFINE <- "AFFINE"
CONVEX <- "CONVEX"
CONCAVE <- "CONCAVE"
UNKNOWN <- "UNKNOWN"

ZERO <- "ZERO"
NONNEG <- "NONNEGATIVE"
NONPOS <- "NONPOSITIVE"

a <- Variable(name = "a")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

canonical_form <- CVXR:::canonical_form
save_value <- CVXR:::save_value

test_that("Test the Variable class", {
  x <- Variable(2, name = "x")
  y <- Variable()
  expect_equal(dim(x), c(2,1))
  expect_equal(dim(y), c(1,1))
  expect_equal(curvature(x), AFFINE)
  
  expect_error(Variable(2,2, diag = TRUE, symmetric = TRUE))
  expect_error(Variable(2,0))
  expect_error(Variable(2,0.5))
})

test_that("Test assigning a value to a variable", {
  # Scalar variable
  a <- Variable()
  value(a) <- 1
  expect_equal(value(a), 1)
  expect_error(value(a) <- c(2,1))

  # Test assigning None
  value(a) <- 1
  value(a) <- NA_real_
  expect_true(is.na(value(a)))

  # Vector variable
  x <- Variable(2)
  value(x) <- c(2,1)
  expect_equal(value(x), matrix(c(2,1)), tolerance = TOL)

  # Matrix variable
  A <- Variable(3, 2)
  value(A) <- matrix(1, nrow = 3, ncol = 2)
  expect_equal(value(A), matrix(1, nrow = 3, ncol = 2))

  # Test assigning negative val to non-negative variable
  x <- Variable(nonneg = TRUE)
  expect_error(value(x) <- -2)

  # Small negative values are rounded to zero
  # *** BEGIN BN EDIT
  # This last test does not seem right given the change
  # made in interface.R for intf_sign!!
  # value(x) <- -1e-8
  # expect_equal(value(x), 0)
  # *** END BN EDIT
})

test_that("Test transposing variables", {
  var <- t(a)
  expect_equal(dim(var), c(1,1))

  a <- save_value(a, 2)
  var <- t(a)
  expect_equal(value(var), 2)

  var <- t(x)
  expect_equal(dim(var), c(1,2))

  x <- save_value(x, matrix(c(1, 2), nrow = 2, ncol = 1))
  var <- t(x)
  expect_equal(value(var)[1,1], 1)
  expect_equal(value(var)[1,2], 2)

  var <- t(C)
  expect_equal(dim(var), c(2,3))

  index <- var[2,1]
  expect_equal(dim(index), c(1,1))

  var <- t(t(x))
  expect_equal(dim(var), c(2,1))
})

test_that("Test the Constant class", {
  c <- Constant(2)
  expect_equal(value(c), matrix(2))
  expect_equal(dim(c), c(1,1))
  expect_equal(curvature(c), CONSTANT)
  expect_equal(sign(c), NONNEG)
  expect_equal(sign(Constant(-2)), NONPOS)
  expect_equal(sign(Constant(0)), ZERO)
  expect_equal(canonical_form(c)[[1]]$size, c(1,1))
  expect_equal(canonical_form(c)[[2]], list())

  # Test the sign
  c <- Constant(matrix(2, nrow = 1, ncol = 2))
  expect_equal(dim(c), c(1,2))
  expect_equal(sign(c), NONNEG)
  expect_equal(sign(-c), NONPOS)
  expect_equal(sign(0*c), ZERO)
  c <- Constant(matrix(c(2, -2), nrow = 1, ncol = 2))
  expect_equal(sign(c), UNKNOWN)

  # Test sign of a complex expression
  c <- Constant(matrix(c(1,2), nrow = 2, ncol = 1))
  Acon <- Constant(matrix(1, nrow = 2, ncol = 2))
  exp <- t(c) %*% Acon %*% c
  expect_equal(sign(exp), NONNEG)
  expect_equal(sign(t(c) %*% c), NONNEG)
  exp <- t(t(c))
  expect_equal(sign(exp), NONNEG)
  exp <- t(c) %*% A
  expect_equal(sign(exp), UNKNOWN)
})

test_that("test R vectors as constants", {
  c <- matrix(c(1,2), nrow = 1, ncol = 2)
  p  <- Parameter(2)
  value(p) <- c(1,1)
  expect_equal(value(c %*% p), 3)
  expect_equal(dim(c %*% x), c(1,1))
})

test_that("test the Parameters class", {
  p <- Parameter(name = "p")
  expect_equal(name(p), "p")
  expect_equal(dim(p), c(1,1))

  p <- Parameter(4, 3, nonneg = TRUE)
  expect_error(value(p) <- 1)

  val <- matrix(-1, nrow = 4, ncol = 3)
  val[1,1] <- 2

  p <- Parameter(4, 3, nonneg = TRUE)
  expect_error(value(p) <- val)

  p <- Parameter(4, 3, nonpos = TRUE)
  expect_error(value(p) <- val)

  # No error for unknown sign
  p <- Parameter(4, 3)
  value(p) <- val

  # Initialize a parameter with a value
  p <- Parameter(value = 10)
  expect_equal(value(p), 10)

  # Test assigning NA
  value(p) <- 10
  value(p) <- NA_real_
  expect_true(is.na(value(p)))

  expect_error(p <- Parameter(2, 1, nonpos = TRUE, value = c(2,1)))
  expect_error(p <- Parameter(4, 3, nonneg = TRUE, value = c(1,2)))
})

test_that("test the AddExpression class", {
  # Vectors
  c <- Constant(c(2,2))
  exp <- x + c
  expect_equal(curvature(exp), AFFINE)
  expect_equal(sign(exp), UNKNOWN)
  expect_equal(canonical_form(exp)[[1]]$size, c(2,1))
  expect_equal(canonical_form(exp)[[2]], list())
  expect_equal(dim(exp), c(2,1))

  z <- Variable(2, name = "z")
  exp <- exp + z + x
  expect_error(x + y)

  # Matrices
  exp <- A + B
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(2,2))
  expect_error(A + C)
  expect_error(AddExpression(A, C))

  # Test that sum is flattened
  exp <- x + c + x
  expect_equal(length(exp@args), 3)
})

test_that("test the SubExpression class", {
  # Vectors
  c <- Constant(c(2,2))
  exp <- x - c
  expect_equal(curvature(exp), AFFINE)
  expect_equal(sign(exp), UNKNOWN)
  expect_equal(canonical_form(exp)[[1]]$size, c(2,1))
  expect_equal(canonical_form(exp)[[2]], list())
  expect_equal(dim(exp), c(2,1))

  z <- Variable(2, name = "z")
  exp <- exp - z - x
  expect_error(x - y)

  # Matrices
  exp <- A - B
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(2,2))
  expect_error(A - C)
})

test_that("test the MulExpression class", {
  # Vectors
  c <- Constant(matrix(2, nrow = 1, ncol = 2))
  expr <- c %*% x
  expect_equal(curvature(expr), AFFINE)
  expect_equal(sign(c[1]*x), UNKNOWN)
  expect_equal(canonical_form(expr)[[1]]$size, c(1,1))
  expect_equal(canonical_form(expr)[[2]], list())
  expect_equal(dim(expr), c(1,1))

  expect_error(matrix(c(2,2,3), nrow = 3, ncol = 1) %*% x)

  # Matrices
  expect_error(Constant(cbind(c(2,1), c(2,2))) %*% C)

  # Affine times affine is okay
  expect_warning(q <- A %*% B)
  expect_true(is_quadratic(q))

  # Non-affine times non-constant raises error
  expect_error(expect_warning(A %*% B) %*% A)

  # Constant expressions
  Tmat <- Constant(cbind(c(1,2,3), c(3,5,5)))
  expr <- (Tmat + Tmat) %*% B
  expect_equal(curvature(expr), AFFINE)
  expect_equal(dim(expr), c(3,2))

  # Expression that would break sign multiplication without promotion
  c <- Constant(matrix(c(2, 2, -2), nrow = 1, ncol = 3))
  expr <- matrix(c(1,2), nrow = 1, ncol = 2) + c %*% C
  expect_equal(sign(expr), UNKNOWN)

  # Scalar constants on the right should be moved left
  expr <- C*2
  expect_equivalent(value(expr@args[[1]]), matrix(2))

  # Scalar variables on the left should be moved right
  expr <- a*c(2,1)
  expect_equivalent(value(expr@args[[1]]), matrix(c(2,1)))
})

test_that("test matrix multiplication operator %*%", {
  # Vectors
  c <- Constant(matrix(2, nrow = 1, ncol = 2))
  exp <- c %*% x
  expect_equal(curvature(exp), AFFINE)
  expect_equal(sign(exp), UNKNOWN)
  expect_equal(canonical_form(exp)[[1]]$size, c(1,1))
  expect_equal(canonical_form(exp)[[2]], list())
  expect_equal(dim(exp), c(1,1))

  # expect_error(x %*% 2)    Note: Allow scalars to be multiplied with %*% to distinguish MulExpression from MulElemwise.
  expect_error(x %*% matrix(c(2,2,3), nrow = 3, ncol = 1))

  # Matrices
  expect_error(Constant(cbind(c(2,1), c(2,2))) %*% C)

  # Affine times affine is okay
  expect_warning(q <- A %*% B)
  expect_true(is_quadratic(q))

  # Non-affine times non-constant raises error
  expect_error(expect_warning(A %*% B %*% A))

  # Constant expressions
  Tmat <- Constant(cbind(c(1,2,3), c(3,5,5)))
  exp <- (Tmat + Tmat) %*% B
  expect_equal(curvature(exp), AFFINE)
  expect_equal(sign(exp), UNKNOWN)

  # Expression that would break sign multiplication without promotion
  c <- Constant(matrix(c(2,2,-2), nrow = 1, ncol = 3))
  exp <- matrix(c(1,2), nrow = 1, ncol = 2) + c %*% C
  expect_equal(sign(exp), UNKNOWN)
})

test_that("test the DivExpression class", {
  # Vectors
  exp <- x/2
  expect_equal(curvature(exp), AFFINE)
  expect_equal(sign(exp), UNKNOWN)
  expect_equal(canonical_form(exp)[[1]]$size, c(2,1))
  expect_equal(canonical_form(exp)[[2]], list())
  expect_equal(dim(exp), c(2,1))

  expect_error(x/c(2,2,3))

  # Constant expressions
  c <- Constant(2)
  exp <- c/(3-5)
  expect_equal(curvature(exp), CONSTANT)
  expect_equal(dim(exp), c(1,1))
  expect_equal(sign(exp), NEGATIVE)

  # Parameters
  p <- Parameter(nonneg = TRUE)
  value(p) <- 2
  exp <- 2/p
  expect_equal(value(exp), 1)

  rho <- Parameter(nonneg = TRUE)
  value(rho) <- 1

  expect_equal(sign(rho), NONNEG)
  expect_equal(sign(Constant(2)), NONNEG)
  expect_equal(sign(Constant(2)/Constant(2)), NONNEG)
  expect_equal(sign(Constant(2)*rho), NONNEG)
  expect_equal(sign(rho/2), NONNEG)
})

test_that("test the NegExpression class", {
  # Vectors
  exp <- -x
  expect_equal(curvature(exp), AFFINE)
  expect_true(is_affine(exp))
  expect_equal(sign(exp), UNKNOWN)
  expect_false(is_nonneg(exp))
  expect_equal(canonical_form(exp)[[1]]$size, c(2,1))
  expect_equal(canonical_form(exp)[[2]], list())
  expect_equal(dim(exp), dim(x))

  # Matrices
  exp <- -C
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(3,2))
})

test_that("test promotion of scalar constants", {
  # Vectors
  exp <- x + 2
  expect_equal(curvature(exp), AFFINE)
  expect_true(is_affine(exp))
  expect_equal(sign(exp), UNKNOWN)
  expect_false(is_nonpos(exp))
  expect_equal(canonical_form(exp)[[1]]$size, c(2,1))
  expect_equal(canonical_form(exp)[[2]], list())
  expect_equal(dim(exp), c(2,1))

  expect_equal(dim(4 - x), c(2,1))
  expect_equal(dim(4 * x), c(2,1))
  expect_equal(dim(4 <= x), c(2,1))
  expect_equal(dim(4 == x), c(2,1))
  expect_equal(dim(x >= 4), c(2,1))

  # Matrices
  exp <- (A + 2) + 4
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(3 * A), c(2,2))
  expect_equal(dim(exp), c(2,2))
})

test_that("test indexing expression", {
  # Tuple of integers as key
  exp <- x[2,1]
  expect_equal(curvature(exp), AFFINE)
  expect_true(is_affine(exp))
  expect_equal(dim(exp), c(1,1))
  expect_equal(value(exp), NA)

  exp <- t(x[2,1])
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(1,1))
  expect_error(x[3,1])

  # Slicing
  exp <- C[1:2,2]
  expect_equal(dim(exp), c(2,1))
  exp <- C[1:nrow(C),1:2]
  expect_equal(dim(exp), c(3,2))
  exp <- C[seq(1, nrow(C), 2), seq(1, ncol(C), 2)]
  expect_equal(dim(exp), c(2,1))
  exp <- C[1:3, seq(1,2,2)]
  expect_equal(dim(exp), c(3,1))
  exp <- C[1:nrow(C),1]
  expect_equal(dim(exp), c(3,1))

  c <- Constant(cbind(c(1,-2), c(0,4)))
  exp <- c[2,2]
  expect_equal(curvature(exp), CONSTANT)
  expect_equal(sign(exp), UNKNOWN)
  expect_equal(sign(c[1,2]), UNKNOWN)
  expect_equal(sign(c[2,1]), UNKNOWN)
  expect_equal(dim(exp), c(1,1))
  expect_equal(value(exp), 4)

  c <- Constant(cbind(c(1,-2,3), c(0,4,5), c(7,8,9)))
  exp <- c[1:3, seq(1,4,2)]
  expect_equal(curvature(exp), CONSTANT)
  expect_true(is_constant(exp))
  expect_equal(dim(exp), c(3,2))
  expect_equal(value(exp[1,2]), 7)

  # Slice of transpose
  exp <- t(C)[1:2,2]
  expect_equal(dim(exp), c(2,1))

  # Arithmetic expression indexing
  exp <- (x + z)[2,1]
  expect_equal(curvature(exp), AFFINE)
  expect_equal(sign(exp), UNKNOWN)
  expect_equal(dim(exp), c(1,1))

  exp <- (x + a)[2,1]
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(1,1))

  exp <- (x - z)[2,1]
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(1,1))

  exp <- (x - a)[2,1]
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(1,1))

  exp <- (-x)[2,1]
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(1,1))

  c <- Constant(rbind(c(1,2), c(3,4)))
  exp <- (c %*% x)[2,1]
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(1,1))

  c <- Constant(rbind(c(1,2), c(3,4)))
  exp <- (c*a)[2,1]
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(1,1))
})

test_that("test negative indices", {
  c <- Constant(rbind(c(1,2), c(3,4)))
  exp <- c[-1,-1]
  expect_equal(value(exp), 4)
  expect_equal(dim(exp), c(1,1))
  expect_equal(curvature(exp), CONSTANT)

  # TODO: More testing of R's negative indices (and sequences of negative indices)
})

test_that("test indexing with logical matrices", {
    ##  require(Matrix)
  A <- rbind(1:4, 5:8, 9:12)
  C <- Constant(A)

  # Logical matrix
  expr <- C[A <= 2]
  expect_equal(dim(expr), c(2,1))
  expect_equal(sign(expr), NONNEG)
  expect_equal(matrix(A[A <= 2]), value(expr))

  expr <- C[A %% 2 == 0]
  expect_equal(dim(expr), c(6,1))
  expect_equal(sign(expr), NONNEG)
  expect_equal(matrix(A[A %% 2 == 0]), value(expr))

  # Logical vector for rows, index for columns
  expr <- C[c(TRUE, FALSE, TRUE), 4]
  expect_equal(dim(expr), c(2,1))
  expect_equal(sign(expr), NONNEG)
  expect_equal(matrix(A[c(TRUE, FALSE, TRUE), 4]), value(expr))

  # Index for rows, logical vector for columns
  expr <- C[2, c(TRUE, FALSE, FALSE, TRUE)]
  expect_equal(dim(expr), c(1, 2))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[2, c(TRUE, FALSE, FALSE, TRUE), drop = FALSE], value(expr))

  # Logical vector for rows, slice for columns
  expr <- C[c(TRUE, TRUE, TRUE), 2:3]
  expect_equal(dim(expr), c(3,2))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[c(TRUE, TRUE, TRUE), 2:3], value(expr))

  # Slice for rows, logical vector for columns
  expr <- C[2:(nrow(C)-1), c(TRUE, FALSE, TRUE, TRUE)]
  expect_equal(dim(expr), c(1,3))    # Always cast 1-D arrays as column vectors. Edit: NOT!!
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[2:(nrow(A)-1), c(TRUE, FALSE, TRUE, TRUE), drop = FALSE], value(expr))

  # Logical vectors for rows and columns
  expr <- C[c(TRUE, TRUE, TRUE), c(TRUE, FALSE, TRUE, TRUE)]
  expect_equal(dim(expr), c(3,3))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[c(TRUE, TRUE, TRUE), c(TRUE, FALSE, TRUE, TRUE)], value(expr))
})

test_that("test indexing with vectors/matrices of indices", {
  A <- rbind(1:4, 5:8, 9:12)
  C <- Constant(A)

  # Vector for rows
  expr <- C[c(1,2)]
  expect_equal(dim(expr), c(2,4))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[c(1,2),], value(expr))

  # Vector for rows, index for columns
  expr <- C[c(1,3),4]
  expect_equal(dim(expr), c(2,1))
  expect_equal(sign(expr), NONNEG)
  expect_equal(matrix(A[c(1,3),4]), value(expr))

  # Index for rows, vector for columns
  expr <- C[2,c(1,3)]
  expect_equal(dim(expr), c(1,2))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[2,c(1,3), drop = FALSE], value(expr))

  # Vector for rows, slice for columns
  expr <- C[c(1,3),2:3]
  expect_equal(dim(expr), c(2,2))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[c(1,3), 2:3], value(expr))

  # Vector for rows and columns
  expr <- C[c(1,2), c(2,4)]
  expect_equal(dim(expr), c(2,2))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[c(1,2), c(2,4)], value(expr))

  # Matrix for rows, vector for columns
  expr <- C[matrix(c(1,2)), c(2,4)]
  expect_equal(dim(expr), c(2,2))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[matrix(c(1,2)), c(2,4)], value(expr))

  # Matrix for rows and columns
  expr <- C[matrix(c(1,2)), matrix(c(2,4))]
  expect_equal(dim(expr), c(2,2))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[matrix(c(1,2)), matrix(c(2,4))], value(expr))
})

test_that("test powers", {
  exp <- x^2
  expect_equal(curvature(exp), CONVEX)
  exp <- x^0.5
  expect_equal(curvature(exp), CONCAVE)
  exp <- x^(-1)
  expect_equal(curvature(exp), CONVEX)
})

test_that("test built-in sum (not good usage)", {
  a_copy <- a
  value(a_copy) <- 1
  expr <- sum(a_copy)
  expect_equal(value(expr), 1)

  x_copy <- x
  value(x_copy) <- c(1,2)
  expr <- sum(x_copy)
  expect_equal(value(expr), 3)
})

test_that("test piecewise linear", {
  A <- matrix(stats::rnorm(6), nrow = 2, ncol = 3)
  b <- matrix(stats::rnorm(2))

  expr <- A %*% y - b
  expect_true(is_pwl(expr))

  expr <- max_elemwise(1, 3*y)
  expect_true(is_pwl(expr))

  expr <- abs(y)
  expect_true(is_pwl(expr))

  expr <- p_norm(3*y, 1)
  expect_true(is_pwl(expr))

  expr <- p_norm(3*y^2, 1)
  expect_false(is_pwl(expr))
})
