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
  skip_on_cran()
  xv <- Variable(2)
  yv <- Variable(2)
  expect_false(name(xv) == name(yv))
  
  xv <- Variable(2, name = "x")
  yv <- Variable()
  expect_equal(name(xv), "x")
  expect_equal(dim(xv), c(2,1))
  expect_equal(dim(yv), c(1,1))
  expect_equal(curvature(xv), AFFINE)
  # expect_equal(dim(canonical_form(x)[[1]]), c(2, 1))
  # expect_equal(canonical_form(x)[[2]], list())
  
  expect_equal(as.character(x), "Variable(2, 1)")
  expect_equal(as.character(A), "Variable(2, 2)")
  
  # Test dim provided as vector instead of tuple
  expect_equal(dim(Variable(c(2), integer = TRUE)), c(2,1))
  
  # # Scalar variable
  # coeff <- coefficients(a)
  # expect_equal(coeff[[as.character(id(a))]], c(1))
  
  # # Vector variable.
  # coeffs <- coefficients(xv)
  # expect_equal(names(coeffs), c(as.character(id(xv))))
  # vec <- coeffs[[as.character(id(xv))]][[1]]
  # expect_equal(dim(vec), c(2,2))
  # expect_equal(vec[1,1], 1)
  
  # # Matrix variable.
  # coeffs <- coefficients(A)
  # expect_equal(names(coeffs), c(as.character(id(A))))
  # expect_equal(length(coeffs[[as.character(id(A))]]), 2) || 0 %in% dim(self)
  # mat <- coeffs[[as.character(id(A))]][[2]]
  # expect_equal(dim(mat), c(2,4))
  # expect_equal(mat[1,3], 1)

  expect_error(Variable(2, 2, diag = TRUE, symmetric = TRUE), 
               "Cannot set more than one special attribute.", fixed = TRUE)
  expect_error(Variable(2, 0), "Invalid dimensions (2,0)", fixed = TRUE)
  expect_error(Variable(2, 0.5), "Invalid dimensions (2,0.5)", fixed = TRUE)
  # expect_error(Variable(2, name = 1), "Variable name 1 must be a string", fixed = TRUE)
})

test_that("Test assigning a value to a variable", {
  skip_on_cran()
  # Scalar variable
  a <- Variable()
  value(a) <- 1
  expect_equal(value(a), matrix(1))
  expect_error(value(a) <- c(2,1), "Invalid dimensions (2,1) for Variable value", fixed = TRUE)

  # Test assigning None
  value(a) <- 1
  value(a) <- NA_real_
  expect_true(is.na(value(a)))

  # Vector variable
  x <- Variable(2)
  value(x) <- c(2,1)
  expect_equal(value(x), matrix(c(2,1)))

  # Matrix variable
  A <- Variable(3, 2)
  value(A) <- matrix(1, nrow = 3, ncol = 2)
  expect_equal(value(A), matrix(1, nrow = 3, ncol = 2))

  # Test assigning negative val to non-negative variable
  x <- Variable(nonneg = TRUE)
  expect_error(value(x) <- -2, "Value must be nonnegative", fixed = TRUE)
})

test_that("Test transposing variables", {
  skip_on_cran()
  var <- t(a)
  expect_equal(name(var), "a")
  expect_equal(dim(var), c(1,1))

  a <- save_value(a, 2)
  var <- t(a)
  expect_equal(value(var), matrix(2))
  
  var <- x
  expect_equal(name(var), "x")
  expect_equal(dim(var), c(2,1))

  xv <- Variable(2, 1, name = "x")
  var <- t(xv)
  expect_equal(name(var), "t(x)")
  expect_equal(dim(var), c(1,2))

  x <- save_value(x, matrix(c(1, 2), nrow = 2, ncol = 1))
  var <- t(x)
  expect_equal(value(var)[1,1], 1)
  expect_equal(value(var)[1,2], 2)

  var <- t(C)
  expect_equal(name(var), "t(C)")
  expect_equal(dim(var), c(2,3))
  
  # coeffs <- coefficients(canonical_form(var)[[1]])
  # mat <- values(coeffs)[[1]][[1]]
  # expect_equal(dim(mat), c(2,6))
  # expect_equal(mat[2,4], 1)

  index <- var[2,1]
  expect_equal(index, "t(C)[2,1]")
  expect_equal(dim(index), c(1,1))

  var <- t(t(x))
  expect_equal(name(var), "t(t(x))")
  expect_equal(dim(var), c(2,1))
})

test_that("Test the Constant class", {
  skip_on_cran()
  c <- Constant(2.0)
  expect_equal(name(c), as.character(2.0))
  
  c <- Constant(2)
  expect_equal(value(c), matrix(2))
  expect_equal(dim(c), c(1,1))
  expect_equal(curvature(c), CONSTANT)
  expect_equal(sign(c), NONNEG)
  expect_equal(sign(Constant(-2)), NONPOS)
  expect_equal(sign(Constant(0)), ZERO)
  # expect_equal(canonical_form(c)[[1]]$dim, c(1,1))
  # expect_equal(canonical_form(c)[[2]], list())
  
  # coeffs <- coefficients(c)
  # expect_equal(names(coeffs), c(CONSTANT))
  # expect_equal(coeffs[[CONSTANT]], c(2))

  # Test the sign
  c <- Constant(matrix(2, nrow = 1, ncol = 2))
  expect_equal(dim(c), c(1,2))
  expect_equal(sign(c), NONNEG)
  expect_equal(sign(-c), NONPOS)
  expect_equal(sign(0*c), ZERO)
  c <- Constant(matrix(c(2, -2), nrow = 1, ncol = 2))
  expect_equal(sign(c), UNKNOWN)
  
  c <- Constant(matrix(0, nrow = 2, ncol = 1))
  expect_equal(dim(c), c(2,1))

  # Test sign of a complex expression
  c <- Constant(matrix(c(1,2), nrow = 2, ncol = 1))
  expect_equal(dim(c), c(2,1))
  Acon <- Constant(matrix(1, nrow = 2, ncol = 2))
  exp <- t(c) %*% Acon %*% c
  expect_equal(sign(exp), NONNEG)
  expect_equal(sign(t(c) %*% c), NONNEG)
  exp <- t(t(c))
  expect_equal(sign(exp), NONNEG)
  exp <- t(c) %*% A
  expect_equal(sign(exp), UNKNOWN)
  
  # Test as.character.
  expect_equal(as.character(c), "Constant(CONSTANT, NONNEGATIVE, (2,1))")
})

test_that("test constant PSD and NSD", {
  skip_on_cran()
  
  n <- 5
  set.seed(0)
  U <- matrix(rnorm(n*n), nrow = n, ncol = n)
  U <- U %*% t(U)
  U <- eigen(U)$vectors   # U is now an orthogonal matrix
  
  # Try four indefinite matrices with different eigenvalue
  # spread around the origin.
  v1 <- matrix(c(3, 2, 1, 1e-8, -1))
  P <- Constant(U %*% diag(v1) %*% t(U))
  expect_false(is_psd(P))
  expect_false(is_nsd(P))
  
  v2 <- matrix(c(3, 2, 2, 1e-6, -1))
  P <- Constant(U %*% diag(v2) %*% t(U))
  expect_false(is_psd(P))
  expect_false(is_nsd(P))
  
  v3 <- matrix(c(3, 2, 2, 1e-4, -1e-6))
  P <- Constant(U %*% diag(v3) %*% t(U))
  expect_false(is_psd(P))
  expect_false(is_nsd(P))
  
  v4 <- matrix(c(-1, 3, 0, 0, 0))
  P <- Constant(U %*% diag(v4) %*% t(U))
  expect_false(is_psd(P))
  expect_false(is_nsd(P))
  
  # Try a test case given in GitHub issue 1451.
  # (Should be equivalent to v4 above).
  P <- Constant(rbind(c(1,2), c(2,1)))
  x <- Variable(2)
  expr <- quad_form(x, P)
  expect_false(is_dcp(expr))
  expect_false(is_dcp(-expr))
  expect_false(gershgorin_psd_check(value(P), tolerance = 0.99))
  
  # Useful Gershgorin disc check
  P <- Constant(rbind(c(2,1), c(1,2)))
  expect_true(gershgorin_psd_check(value(P), tolerance = 0.0))
  
  # Verify good behavior for large eigenvalues
  P <- Constant(diag(c(rep(1e-4, 9), -1e4)))
  expect_false(is_psd(P))
  expect_false(is_nsd(P))
  
  # Check a case when the matrix is in fact PSD.
  P <- Constant(matrix(1, nrow = 5, ncol = 5))
  expect_true(is_psd(P))
  expect_false(is_nsd(P))
  
  # Check with sparse inputs
  P <- Constant(Matrix(diag(10), sparse = TRUE))
  expect_true(gershgorin_psd_check(value(P), EIGVAL_TOL))
  expect_true(is_psd(P))
  expect_true(is_nsd(-P))
  Q <- -EIGVAL_TOL/2 * P
  expect_true(gershgorin_psd_check(value(Q), EIGVAL_TOL))
  Q <- -1.1*EIGVAL_TOL*P
  expect_false(gershgorin_psd_check(value(Q), EIGVAL_TOL))
  expect_false(is_psd(Q))
})

test_that("test constant skew symmetric", {
  # Define inputs
  M1_false <- diag(3)
  M2_true <- matrix(0, nrow = 3, ncol = 3)
  M3_true <- rbind(c(0,1), c(-1,0))
  M4_true <- rbind(c(0,-1), c(1,0))
  M5_false <- rbind(c(0,1), c(1,0))
  M6_false <- rbind(c(1,1), c(-1,0))
  M7_false <- rbind(c(0,1), c(-1.1,0))
  
  # Test dense constants
  C <- Constant(M1_false)
  expect_false(is_skew_symmetric(C))
  C <- Constant(M2_true)
  expect_true(is_skew_symmetric(C))
  C <- Constant(M3_true)
  expect_true(is_skew_symmetric(C))
  C <- Constant(M4_true)
  expect_true(is_skew_symmetric(C))
  C <- Constant(M5_false)
  expect_false(is_skew_symmetric(C))
  C <- Constant(M6_false)
  expect_false(is_skew_symmetric(C))
  C <- Constant(M7_false)
  expect_false(is_skew_symmetric(C))
  
  # Test sparse constants
  C <- Constant(Matrix(M1_false, sparse = TRUE))
  expect_false(is_skew_symmetric(C))
  C <- Constant(Matrix(M2_true, sparse = TRUE))
  expect_true(is_skew_symmetric(C))
  C <- Constant(Matrix(M3_true, sparse = TRUE))
  expect_true(is_skew_symmetric(C))
  C <- Constant(Matrix(M4_true, sparse = TRUE))
  expect_true(is_skew_symmetric(C))
  C <- Constant(Matrix(M5_false, sparse = TRUE))
  expect_false(is_skew_symmetric(C))
  C <- Constant(Matrix(M6_false, sparse = TRUE))
  expect_false(is_skew_symmetric(C))
  C <- Constant(Matrix(M7_false, sparse = TRUE))
  expect_false(is_skew_symmetric(C))
  
  # Test complex inputs: never recognized as skew-symmetric.
  C <- Constant(1i * M2_true)
  # From a mathematical standpoint one can argue that this should
  # be true, but I don't think there's precedent for CVXR
  # automatically converting complex expressions with zero imaginary-part
  # into equivalent real expressions.
  expect_false(is_skew_symmetric(C))
  C <- Constant(1i * M3_true)
  expect_false(is_skew_symmetric(C))
  C <-  Constant(1i * M4_true)
  expect_false(is_skew_symmetric(C))
})

test_that("test R vectors as constants", {
  skip_on_cran()
  c <- matrix(c(1,2), nrow = 1, ncol = 2)
  p  <- Parameter(2)
  value(p) <- c(1,1)
  expect_equal(value(c %*% p), matrix(3))
  expect_equal(dim(c %*% x), c(1,1))
})

test_that("test Parameter class on good inputs", {
  skip_on_cran()
  p <- Parameter(name = "p")
  expect_equal(name(p), "p")
  expect_equal(dim(p), c(1,1))
  
  # Entry-wise constraints on parameter values.
  val <- -matrix(1, nrow = 4, ncol = 3)
  val[1, 1] <- 2
  p <- Parameter(4, 3)
  value(p) <- val

  # Initialize a parameter with a value; later set it to NA.
  p <- Parameter(value = 10)
  expect_equal(value(p), matrix(10))
  value(p) <- 10
  value(p) <- NA_real_
  expect_true(is.na(value(p)))
  
  # Test parameter representation.
  p <- Parameter(4, 3, nonpos = TRUE)
  expect_equal(as.character(p), "Parameter(4, 3, nonpos = TRUE)")

  # Test valid diagonal parameter.
  p <- Parameter(2, 2, diag = TRUE)
  value(p) <- Matrix(diag(2), sparse = TRUE)
  expect_equal(as.matrix(value(p)), diag(2), check.attributes = FALSE, tolerance = 1e-10)
})

test_that("test PSD and NSD parameters", {
  skip_on_cran()
  
  # Test valid rank-deficient PSD parameter.
  set.seed(42)
  a <- matrix(rnorm(100*95), nrow = 100, ncol = 95)
  a2 <- a %*% t(a)   # This must be a PSD matrix.
  p <- Parameter(100, 100, PSD = TRUE)
  value(p) <- a2
  expect_equal(value(p), a2, tolerance = 1e-10)

  # Test positive definite matrix with non-distinct eigenvalues
  m <- 10
  n <- 5
  A <- matrix(rnorm(m*n), nrow = m, ncol = n) + 1i * matrix(rnorm(m*n), nrow = m, ncol = n)   # a random complex matrix
  A <- Conj(t(A)) %*% A   # a random Hermitian positive definite matrix
  A <- rbind(cbind(Re(A), -Im(A)), cbind(Im(A), Re(A)))

  p <- Parameter(2*n, 2*n, PSD = TRUE)
  value(p) <- A
  expect_equal(value(p), A, tolerance = TOL)

  # Test invalid PSD and NSD parameters
  n <- 5
  P <- Parameter(n, n, PSD = TRUE)
  N <- Parameter(n, n, NSD = TRUE)
  
  set.seed(0)
  U <- matrix(rnorm(n*n), nrow = n, ncol = n)
  U <- U %*% t(U)
  U <- eigen(U)$vectors   # U is now an orthogonal matrix
  v1 <- matrix(c(3, 2, 1, 1e-8, -1))
  v2 <- matrix(c(3, 2, 2, 1e-6, -1))
  v3 <- matrix(c(3, 2, 2, 1e-4, -1e-6))
  v4 <- matrix(c(-1, 3, 0, 0, 0))
  vs <- list(v1, v2, v3, v4)
  
  for(vi in vs) {
    expect_error(value(P) <- U %*% diag(vi) %*% t(U),
                 "Parameter value must be positive semidefinite.", fixed = TRUE)
    expect_error(value(N) <- -U %*% diag(vi) %*% t(U),
                 "Parameter value must be negative semidefinite.", fixed = TRUE)
  }
})

test_that("test the Parameter class on bad inputs",{
  skip_on_cran()
  p <- Parameter(name = 'p')
  expect_equal(name(p), "p")
  expect_equal(dim(p), c(1,1))

  p <- Parameter(4, 3, nonneg = TRUE)
  expect_error(value(p) <- 1, "Invalid dimensions c(1,1) for Parameter value.", fixed = TRUE)
  # expect_error(value(p) <- c(1,1), "Invalid dimensions (2,1) for Parameter value.", fixed = TRUE)

  val <- -matrix(1, nrow = 4, ncol = 3)
  val[1,1] <- 2

  p <- Parameter(4, 3, nonneg = TRUE)
  expect_error(value(p) <- val, "Parameter value must be nonnegative.", fixed = TRUE)

  p <- Parameter(4, 3, nonpos = TRUE)
  expect_error(value(p) <- val, "Parameter value must be nonpositive.", fixed = TRUE)

  expect_error(p <- Parameter(2, 1, nonpos = TRUE, value = matrix(c(2,1))),
               "Parameter value must be nonpositive", fixed = TRUE)

  expect_error(p <- Parameter(4, 3, nonneg = TRUE, value = matrix(c(1,2))),
               "Invalid dimensions (2,1) for Parameter value.", fixed = TRUE)

  expect_error(p <- Parameter(2, 2, diag = TRUE, symmetric = TRUE),
               "Cannot set more than one special attribute in Parameter.", fixed = TRUE)

  # Boolean
  expect_error(p <- Parameter(2, 2, boolean = TRUE, value = rbind(c(1,1), c(1,-1))),
               "Parameter value must be boolean.", fixed = TRUE)

  # Integer
  expect_error(p <- Parameter(2, 2, integer = TRUE, value = rbind(c(1,1.5), c(1,-1))),
               "Parameter value must be integer.", fixed = TRUE)

  # Diag
  expect_error(p <- Parameter(2, 2, diag = TRUE, value = rbind(c(1,1), c(1,-1))),
               "Parameter value must be diagonal.", fixed = TRUE)

  # Symmetric
  expect_error(p <- Parameter(2, 2, symmetric = TRUE, value = rbind(c(1,1), c(-1,-1))),
               "Parameter value must be symmetric.", fixed = TRUE)
})

test_that("test symmetric variables",{
  skip_on_cran()
  expect_error(v <- Variable(4, 3, symmetric = TRUE),
               "Invalid dimensions (4,3). Must be a square matrix.", fixed = TRUE)

  v <- Variable(2, 2, symmetric = TRUE)
  expect_true(is_symmetric(v))
  v <- Variable(2, 2, PSD = TRUE)
  expect_true(is_symmetric(v))
  v <- Variable(2, 2, NSD = TRUE)
  expect_true(is_symmetric(v))
  v <- Variable(2, 2, diag = TRUE)
  expect_true(is_symmetric(v))
  expect_true(is_symmetric(a))
  expect_false(is_symmetric(A))

  v <- Variable(2, 2, symmetric = TRUE)
  expr <- v + v
  expect_true(is_symmetric(expr))
  expr <- -v
  expect_true(is_symmetric(expr))
  expr <- t(v)
  expect_true(is_symmetric(expr))
  expr <- Re(v)
  expect_true(is_symmetric(expr))
  expr <- Im(v)
  expect_true(is_symmetric(expr))
  expr <- Conj(v)
  expect_true(is_symmetric(expr))
  expr <- Promote(Variable(), c(2,2))
  expect_true(is_symmetric(expr))

})

test_that("test Hermitian variables", {
  skip_on_cran()
  expect_error(v <- Variable(4, 3, hermitian = TRUE),
               "Invalid dimensions (4,3). Must be a square matrix.", fixed = TRUE)

  v <- Variable(2, 2, hermitian = TRUE)
  expect_true(is_hermitian(v))
  # v <- Variable(2, 2, PSD = TRUE)
  # expect_true(is_symmetric(v))
  # v <- Variable(2, 2, NSD = TRUE)
  # expect_true(is_symmetric(v))
  v <- Variable(2, 2, diag = TRUE)
  expect_true(is_hermitian(v))


  v <- Variable(2, 2, hermitian = TRUE)
  expr <- v + v
  expect_true(is_hermitian(expr))
  expr <- -v
  expect_true(is_hermitian(expr))
  expr <- t(v)
  expect_true(is_hermitian(expr))
  expr <- Re(v)
  expect_true(is_hermitian(expr))
  expr <- Im(v)
  expect_true(is_hermitian(expr))
  expr <- Conj(v)
  expect_true(is_hermitian(expr))
  expr <- Promote(Variable(), c(2,2))
  expect_true(is_hermitian(expr))
})

test_that("test rounding for attributes", {
  skip_on_cran()

  # Nonpos
  v <- Variable(1, nonpos = TRUE)
  expect_equal(project(v, 1), 0)
  v <- Variable(2, nonpos = TRUE)
  expect_equal(project(v, c(1, -1)), c(0,-1))

  # Nonneg
  v <- Variable(1, nonneg = TRUE)
  expect_equal(project(v, -1), 0)
  v <- Variable(2, nonneg = TRUE)
  expect_equal(project(v, c(1, -1)), c(1,0))

  # Boolean
  v <- Variable(2, 2, boolean = TRUE)
  expect_equal(project(v, cbind(c(1, -1), c(1, 0))), c(1, 0, 1, 0), check.attributes = FALSE)

  # Integer
  v <- Variable(2, 2, integer = TRUE)
  expect_equal(project(v, cbind(c(1, -1.6), c(1, 0))), c(1, -2, 1, 0), check.attributes = FALSE)

  # Symmetric
  v <- Variable(2, 2, symmetric = TRUE)
  expect_equal(project(v, rbind(c(1, -1), c(1, 0))), c(1, 0, 0, 0), check.attributes = FALSE)

  # PSD
  v <- Variable(2, 2, PSD = TRUE)
  expect_equal(project(v, rbind(c(1, -1), c(1, -1))), c(1, 0, 0, 0), check.attributes = FALSE)

  # NSD
  v <- Variable(2, 2, NSD = TRUE)
  expect_equal(project(v, rbind(c(1, -1), c(1, -1))), c(0, 0, 0, -1), check.attributes = FALSE)

  # diag
  v <- Variable(2, 2, diag = TRUE)
  expect_equal(matrix(project(v, rbind(c(1, -1), c(1, 0)))), c(1, 0, 0, 0), check.attributes = FALSE)

  # Hermitian
  v <- Variable(2, 2, hermitian = TRUE)
  expect_equal(project(v, rbind(c(1, -1i), c(1, 0))), c(1, 0.5 + 0.5i, 0.5 - 0.5i, 0), check.attributes = FALSE)

  A <- Constant(1.0)
  expect_true(is_psd(A))
  expect_false(is_nsd(A))
  A <- Constant(-1.0)
  expect_false(is_psd(A))
  expect_true(is_nsd(A))
  A <- Constant(0.0)
  expect_true(is_psd(A))
  expect_true(is_nsd(A))
})

test_that("test the AddExpression class", {
  skip_on_cran()
  
  # Vectors
  c <- Constant(c(2,2))
  exp <- x + c
  expect_equal(curvature(exp), AFFINE)
  expect_equal(sign(exp), UNKNOWN)
  # expect_equal(canonical_form(exp)[[1]]$dim, c(2,1))
  # expect_equal(canonical_form(exp)[[2]], list())
  # expect_equal(name(exp), paste(name(x), "+", name(c)))
  expect_equal(dim(exp), c(2,1))

  z <- Variable(2, name = "z")
  exp <- exp + z + x
  
  # Incompatible dimensions
  expect_error(x + y)

  # Matrices
  exp <- A + B
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(2,2))
  
  # Incompatible dimensions
  expect_error(A + C)
  
  # Incompatible dimensions
  expect_error(AddExpression(A, C))

  # Test that sum is flattened
  exp <- x + c + x
  expect_equal(length(exp@args), 3)
  
  # Test repr
  expect_equal(as.character(exp), "Expression(AFFINE, UNKNOWN, (2,1))")
})

test_that("test the SubExpression class", {
  skip_on_cran()
  # Vectors
  c <- Constant(c(2,2))
  exp <- x - c
  expect_equal(curvature(exp), AFFINE)
  expect_equal(sign(exp), UNKNOWN)
  # expect_equal(canonical_form(exp)[[1]]$dim, c(2,1))
  # expect_equal(canonical_form(exp)[[2]], list())
  # expect_equal(name(exp), paste(name(x), "-", name(Constant(c(2,2)))))
  expect_equal(dim(exp), c(2,1))

  z <- Variable(2, name = "z")
  exp <- exp - z - x
  
  # Incompatible dimensions
  expect_error(x - y)

  # Matrices
  exp <- A - B
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(2,2))
  
  # Incompatible dimensions
  expect_error(A - C)
  
  # Test repr
  expect_equal(as.character(exp), "Expression(AFFINE, UNKNOWN, (2,1))")
})

test_that("test the MulExpression class", {
  skip_on_cran()
  
  # Vectors
  c <- Constant(matrix(2, nrow = 1, ncol = 2))
  exp <- c %*% x
  expect_equal(curvature(exp), AFFINE)
  expect_equal(sign(c[1]*x), UNKNOWN)
  # expect_equal(canonical_form(exp)[[1]]$dim, c(1,1))
  # expect_equal(canonical_form(exp)[[2]], list())
  # expect_equal(name(exp), paste(name(c), "*", name(x)))
  expect_equal(dim(exp), c(1,1))

  # Incompatible dimensions
  expect_error(matrix(c(2,2,3), nrow = 3, ncol = 1) %*% x)

  # Matrices: incompatible dimensions
  expect_error(Constant(rbind(c(2,1), c(2,2))) %*% C)

  # Affine times affine is okay
  expect_warning(q <- A %*% B)
  expect_true(is_quadratic(q))

  # Constant expressions
  Tmat <- Constant(rbind(c(1,2,3), c(3,5,5)))
  expr <- (Tmat + Tmat) %*% B
  expect_equal(curvature(expr), AFFINE)
  expect_equal(dim(expr), c(3,2))

  # Expression that would break sign multiplication without promotion
  c <- Constant(matrix(c(2, 2, -2), nrow = 1, ncol = 3))
  exp <- matrix(c(1,2), nrow = 1, ncol = 2) + c %*% C
  expect_equal(sign(exp), UNKNOWN)
})

test_that("test matrix multiplication operator %*%", {
  skip_on_cran()
  
  # Vectors
  c <- Constant(matrix(2, nrow = 1, ncol = 2))
  exp <- c %*% x
  expect_equal(curvature(exp), AFFINE)
  expect_equal(sign(exp), UNKNOWN)
  # expect_equal(canonical_form(exp)[[1]]$dim, c(1,1))
  # expect_equal(canonical_form(exp)[[2]], list())
  expect_equal(name(exp), paste(name(c), "%*%", name(x)))
  expect_equal(dim(exp), c(1,1))

  # Note: Allow scalars to be multiplied with %*% to distinguish MulExpression from MulElemwise.
  # expect_error(x %*% 2, "Scalar operands are not allowed, use '*' instead", fixed = TRUE)
  
  # Incompatible dimensions
  expect_error(x %*% matrix(c(2,2,3), nrow = 3, ncol = 1))

  # Matrices
  expect_error(Constant(rbind(c(2,1), c(2,2))) %*% C)

  # Affine times affine is okay
  expect_warning(q <- A %*% B)
  expect_true(is_quadratic(q))

  # Non-affine times non-constant raises error
  # expect_error(expect_warning(A %*% B %*% A))

  # Constant expressions
  Tmat <- Constant(rbind(c(1,2,3), c(3,5,5)))
  exp <- (Tmat + Tmat) %*% B
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(3,2))
  
  # Expression that would break sign multiplication without promotion
  c <- Constant(matrix(c(2,2,-2), nrow = 1, ncol = 3))
  exp <- matrix(c(1,2), nrow = 1, ncol = 2) + c %*% C
  expect_equal(sign(exp), UNKNOWN)

  # Testing dim.
  a <- Parameter(1)
  x <- Variable(1)
  expr <- a %*% x
  expect_equal(dim(expr), c(1,1))

  A <- Parameter(4, 4)
  z <- Variable(4, 1)
  expr <- A %*% z
  expect_equal(dim(expr), c(4,1))

  v <- Variable(1,1)
  col_scalar <- Parameter(1,1)
  expect_true(identical(dim(v), dim(col_scalar), dim(t(col_scalar))))
})

test_that("test the DivExpression class", {
  skip_on_cran()
  
  # Vectors
  exp <- x/2
  expect_equal(curvature(exp), AFFINE)
  expect_equal(sign(exp), UNKNOWN)
  # expect_equal(canonical_form(exp)[[1]]$dim, c(2,1))
  # expect_equal(canonical_form(exp)[[2]], list())
  # expect_equal(name(exp), paste(name(c), "*", name(x)))
  expect_equal(dim(exp), c(2,1))

  expect_error(x/c(2,2,3), "Incompatible dimensions for division", fixed = TRUE)

  # Constant expressions
  c <- Constant(2)
  exp <- c/(3-5)
  expect_equal(curvature(exp), CONSTANT)
  expect_equal(dim(exp), c(1,1))
  expect_equal(sign(exp), NONPOS)

  # Parameters
  p <- Parameter(nonneg = TRUE)
  exp <- 2/p
  value(p) <- 2
  expect_equal(value(exp), 1)

  rho <- Parameter(nonneg = TRUE)
  value(rho) <- 1

  expect_equal(sign(rho), NONNEG)
  expect_equal(sign(Constant(2)), NONNEG)
  expect_equal(sign(Constant(2)/Constant(2)), NONNEG)
  expect_equal(sign(Constant(2)*rho), NONNEG)
  expect_equal(sign(rho/2), NONNEG)
  
  # Broadcasting.
  x <- Variable(3, 3)
  c <- matrix(1:3, nrow = 3, ncol = 1)
  expr <- x / c
  expect_equal(dim(expr), c(3, 3))
  value(x) <- matrix(1, nrow = 3, ncol = 3)
  A <- matrix(1, nrow = 3, ncol = 3) / c
  expect_equal(A, value(expr))
  # expect_error(x / 1:3, "Incompatible dimensions for division.")
})

test_that("test the NegExpression class", {
  skip_on_cran()
  
  # Vectors
  exp <- -x
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(2,1))
  expect_true(is_affine(exp))
  expect_equal(sign(exp), UNKNOWN)
  expect_false(is_nonneg(exp))
  # expect_equal(canonical_form(exp)[[1]]$dim, c(2,1))
  # expect_equal(canonical_form(exp)[[2]], list())
  # expect_equal(name(exp), paste("-", name(x), sep = ""))
  expect_equal(dim(exp), dim(x))

  # Matrices
  exp <- -C
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(3,2))
})

test_that("test promotion of scalar constants", {
  skip_on_cran()
  
  # Vectors
  exp <- x + 2
  expect_equal(curvature(exp), AFFINE)
  expect_true(is_affine(exp))
  expect_equal(sign(exp), UNKNOWN)
  expect_false(is_nonpos(exp))
  # expect_equal(canonical_form(exp)[[1]]$dim, c(2,1))
  # expect_equal(canonical_form(exp)[[2]], list())
  # expect_equal(name(exp), paste(name(x), "+", name(Constant(2))))
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
  skip_on_cran()
  
  # Tuple of integers as key
  exp <- x[2,1]
  # expect_equal(name(exp), "x[2,1]")
  expect_equal(curvature(exp), AFFINE)
  expect_true(is_affine(exp))
  expect_equal(dim(exp), c(1,1))
  # coeff <- coefficients(canonical_form(exp)[[1]])[[as.character(id(x))]][[1]]
  # expect_equal(coeff[1,2], 1)
  expect_equal(value(exp), NA_real_)

  exp <- t(x[2,1])
  # expect_equal(name(exp), "x[2,1]")
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(1,1))
  
  expect_error(x[3,1], "Too many indices for expression", fixed = TRUE)
  expect_error(x[3], "Index 2 is out of bounds for axis 2 with size 2", fixed = TRUE)

  # Slicing
  exp <- C[1:2,2]
  # expect_equal(name(exp), "C[1:2,2]")
  expect_equal(dim(exp), c(2,1))
  exp <- C[1:nrow(C),1:2]
  # expect_equal(name(exp), "C[1:nrow(C),1:2]")
  expect_equal(dim(exp), c(3,2))
  exp <- C[seq(1, nrow(C), 2), seq(1, ncol(C), 2)]
  # expect_equal(name(exp), "C[seq(1::3,1::3)]")
  expect_equal(dim(exp), c(2,1))
  exp <- C[1:3, seq(1,2,2)]
  # expect_equal(name(exp), "C[1:3,1]")
  expect_equal(dim(exp), c(3,1))
  exp <- C[1:nrow(C),1]
  # expect_equal(name(exp), "C[1:,1]")
  expect_equal(dim(exp), c(3,1))

  c <- Constant(rbind(c(1,-2), c(0,4)))
  exp <- c[2,2]
  expect_equal(curvature(exp), CONSTANT)
  expect_equal(sign(exp), UNKNOWN)
  expect_equal(sign(c[1,2]), UNKNOWN)
  expect_equal(sign(c[2,1]), UNKNOWN)
  expect_equal(dim(exp), c(1,1))
  expect_equal(value(exp), matrix(4))

  c <- Constant(rbind(c(1,-2,3), c(0,4,5), c(7,8,9)))
  exp <- c[1:3, seq(1,4,2)]
  expect_equal(curvature(exp), CONSTANT)
  expect_true(is_constant(exp))
  expect_equal(dim(exp), c(3,2))
  expect_equal(value(exp[1,2]), matrix(7))

  # Slice of transpose
  exp <- t(C)[1:2,2]
  expect_equal(dim(exp), c(2,1))

  # Arithmetic expression indexing
  exp <- (x + z)[2,1]
  # expect_equal(name(exp), "x[2,1] + z[2,1]")
  expect_equal(curvature(exp), AFFINE)
  expect_equal(sign(exp), UNKNOWN)
  expect_equal(dim(exp), c(1,1))

  exp <- (x + a)[2,1]
  # expect_equal(name(exp), "x[2,1] + a")
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(1,1))

  exp <- (x - z)[2,1]
  # expect_equal(name(exp), "x[2,1] - z[2,1]")
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(1,1))

  exp <- (x - a)[2,1]
  # expect_equal(name(exp), "x[2,1] - a")
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(1,1))

  exp <- (-x)[2,1]
  # expect_equal(name(exp), "-x[2,1]")
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(1,1))

  c <- Constant(rbind(c(1,2), c(3,4)))
  exp <- (c %*% x)[2,1]
  # expect_equal(name(exp), "(2,4)^T * x[1:,1]")
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(1,1))

  c <- Constant(rbind(c(1,2), c(3,4)))
  exp <- (c*a)[2,1]
  # expect_equal(name(exp), "2*a")
  expect_equal(curvature(exp), AFFINE)
  expect_equal(dim(exp), c(1,1))
})

# Note: Not really applicable to CVXR, but leaving here for records.
# test_that("test NA as idx", {
#   skip_on_cran()
#   expr <- a[NA, NA]
#   expect_equal(dim(expr), c(1, 1))
#
#   expr <- x[, NA]
#   expect_equal(dim(expr), c(2, 1))
#
#   expr <- x[NA, ]
#   expect_equal(dim(expr), c(1, 2))
#
#   expr <- Constant(c(1,2))[NA,]
#   expect_equal(dim(expr), c(1, 2))
#   expect_equal(value(expr), c(1,2))
#
# })

test_that("test out of bounds indices", {
  skip_on_cran()
  # TODO: These tests, especially for negative indices, need to be adapted to CVXR.
  
  expect_error(x[100], "Index 100 is out of bounds for axis 1 with size 2", fixed = TRUE)
  expect_error(x[-100], "Index -100 is out of bounds for axis 1 with size 2", fixed = TRUE)

  exp <- x[-100]
  expect_equal(size(exp), 1)
  expect_equal(value(exp), c())
  
  exp <- C[100:2]
  expect_equal(dim(exp), c(1,2))
  
  exp <- C[,-199:2]
  expect_equal(dim(exp), c(3,2))
  
  exp <- C[,-199:-3]
  expect_equal(dim(exp), c(3,1))
})

test_that("test negative indices", {
  skip_on_cran()
  # TODO: Double check indexing matches Python unit tests.
  
  c <- Constant(rbind(c(1,2), c(3,4)))
  exp <- c[-1,-1]
  expect_equal(value(exp), matrix(4))
  expect_equal(dim(exp), c(1,1))
  expect_equal(curvature(exp), CONSTANT)

  c <- Constant(1:4)
  exp <- c[2:(nrow(c)-1)]
  expect_equal(value(exp), matrix(c(2,3)))
  expect_equal(dim(exp), c(2,1))
  expect_equal(curvature(exp), CONSTANT)

  c <- Constant(1:4)
  exp <- c[seq(4,1,-1)]
  expect_equal(value(exp), matrix(c(4,3,2,1)))
  expect_equal(dim(exp), c(4,1))
  expect_equal(curvature(exp), CONSTANT)

  x <- Variable(4)
  expect_equal(dim(x[seq(4,1,-1)]), c(4,1))
  prob <- Problem(Minimize(0), list(x[seq(4,1,-1)] == c))
  result <- solve(prob, solver = "SCS")
  expect_equal(result$getValue(x), matrix(c(4,3,2,1)))

  x <- Variable(2)
  expect_equal(dim(x[seq(4,1,-1)]), c(2,1))
  
  c <- Constant(rbind(c(1,2), c(3,4)))
  expr <- c[1, seq(2,1,-1)]
  expect_equal(dim(expr), c(1,1))
  expect_equal(value(expr), matrix(3))
  
  expr <- c[1, seq(2,1,-1)]
  expect_equal(dim(expr), c(2,1))
  expect_equal(value(expr), matrix(c(3,1)))
})

test_that("test indexing with logical matrices", {
  skip_on_cran()
  
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
  expr <- C[2, c(TRUE, FALSE, FALSE, TRUE)]   # TODO: Add option for drop = TRUE in overloaded [.
  expect_equal(dim(expr), c(1, 2))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[2, c(TRUE, FALSE, FALSE, TRUE), drop = FALSE], value(expr))

  # Logical vector for rows, slice for columns
  expr <- C[c(TRUE, TRUE, TRUE), 2:3]
  expect_equal(dim(expr), c(3,2))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[c(TRUE, TRUE, TRUE), 2:3], value(expr))

  # Slice for rows, logical vector for columns
  expr <- C[2:(nrow(C)-1), c(TRUE, FALSE, TRUE, TRUE)]   # TODO: Add option for drop = TRUE in overloaded [.
  expect_equal(dim(expr), c(1,3))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[2:(nrow(A)-1), c(TRUE, FALSE, TRUE, TRUE), drop = FALSE], value(expr))

  # Logical vectors for rows and columns
  expr <- C[c(TRUE, TRUE, TRUE), c(TRUE, FALSE, TRUE, TRUE)]
  expect_equal(dim(expr), c(3,3))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[c(TRUE, TRUE, TRUE), c(TRUE, FALSE, TRUE, TRUE)], value(expr))
})

test_that("test indexing with vectors/matrices of indices", {
  skip_on_cran()
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
  expr <- C[2,c(1,3)]   # TODO: Add option for drop = TRUE in overloaded [.
  expect_equal(dim(expr), c(1,2))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[2,c(1,3), drop = FALSE], value(expr))

  # Vector for rows, slice for columns
  expr <- C[c(1,3),2:3]
  expect_equal(dim(expr), c(2,2))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[c(1,3), 2:3], value(expr))

  # Vector for rows and columns
  expr <- C[2:(nrow(C)-1), c(2,4)]
  expect_equal(dim(expr), c(2,2))
  expect_equal(sign(expr), NONNEG)
  expect_equal(A[2:(nrow(C)-1), c(2,4)], value(expr))

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
  skip_on_cran()
  exp <- x^2
  expect_equal(curvature(exp), CONVEX)
  exp <- x^0.5
  expect_equal(curvature(exp), CONCAVE)
  exp <- x^(-1)
  expect_equal(curvature(exp), CONVEX)
})

test_that("test sum_entries function", {
  skip_on_cran()
  a_copy <- a
  value(a_copy) <- 1
  expr <- sum_entries(a_copy)
  expect_equal(value(expr), 1)

  x_copy <- x
  value(x_copy) <- c(1,2)
  expr <- sum_entries(x_copy)
  expect_equal(value(expr), 3)
})

test_that("test copy function for variable types", {
  skip_on_cran()
  x <- Variable(3, 4, name = "x")
  y <- copy(x)
  expect_equal(dim(y), c(3,4))
  expect_equal(name(y), "x")
  
  x <- Variable(5, 5, PSD = TRUE, name = "x")
  y <- copy(x)
  expect_equal(dim(y), c(5,5))
})

test_that("test copy function for Parameters", {
  skip_on_cran()
  x <- Parameter(3, 4, name = "x", nonneg = TRUE)
  y <- copy(x)
  expect_equal(dim(y), c(3,4))
  expect_equal(name(y), "x")
  expect_equal(sign(y), NONNEG)
})

test_that("test copy function for Constants", {
  skip_on_cran()
  x <- Constant(2)
  y <- copy(x)
  expect_equal(dim(y), c(1,1))
  expect_equal(value(y), 2)
})

test_that("test piecewise linear", {
  skip_on_cran()
  A <- matrix(1, nrow = 2, ncol = 3)
  b <- matrix(c(1, 1))

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

test_that("test multiply broadcasting", {
  skip_on_cran()
  
  y <- Parameter(3, 1)
  z <- Variable(1, 3)
  value(y) <- 1:2
  value(z) <- (1:2 - 1)
  expr <- multiply(y, z)
  expect_equal(value(expr), value(y)*value(z))
  
  prob <- Problem(Minimize(sum(expr)), list(z == value(z)))
  result <- solve(prob, solver = "SCS")
  expect_equal(result$getValue(expr), result$getValue(y)*result$getValue(z))
  
  set.seed(0)
  m <- 3
  n <- 4
  A <- matrix(runif(m*n), nrow = m, ncol = n)
  
  col_scale <- Variable(n)
  expect_error(multiply(A, col_scale), "Cannot broadcast dimensions (3,4) (4,1)", fixed = TRUE)
  
  col_scale <- Variable(m, 1)
  R <- multiply(A, row_scale)
  expect_equal(dim(R), c(m,n))
})

test_that("test addition broadcasting", {
  skip_on_cran()
  
  y <- Parameter(3, 1)
  z <- Variable(1, 3)
  value(y) <- 1:2
  value(z) <- 1:2 - 1
  expr <- y + z
  expect_equal(value(expr), value(y) + value(z))

  prob <- Problem(Minimize(sum(expr)), list(z == value(z)))
  result <- solve(prob, solver = "SCS")
  expect_equal(result$getValue(expr), result$getValue(y) + result$getValue(z))
  
  set.seed(0)
  m <- 3
  n <- 4
  A <- matrix(rnorm(m*n), nrow = m, ncol = n)
  
  col_scale <- Variable(n)
  
  expect_error(A + col_scale, "Cannot broadcast dimensions (3, 4) (4, 1)", fixed = TRUE)
  
  col_scale <- Variable(c(1, n))
  C <- A + col_scale
  expect_equal(dim(C), c(m, n))
  
  row_scale <- Variable(c(m, 1))
  R <- A + row_scale
  expect_equal(dim(R), c(m, n))
})

test_that("test curvature string is populated for log-log expressions", {
  skip_on_cran()
  xv <- Variable(pos = TRUE)
  monomial <- x*x*x
  expect_equal(curvature(monomial), LOG_LOG_AFFINE)
  
  posynomial <- x*x*x + x
  expect_equal(curvature(posynomial), LOG_LOG_CONVEX)
  
  llcv <- 1/(x*x*x + x)
  expect_equal(curvature(llcv), LOG_LOG_CONCAVE)
})

test_that("test conversion of native t(x) %*% A %*% x into QuadForms", {
  # Trivial quad form
  xv <- Variable(2)
  Ac <- Constant(rbind(c(1, 0), c(0, -1)))
  expr <- t(xv) %*% Ac %*% xv
  expect_true(is(expr, QuadForm))
  
  # QuadForm inside nested expr: 0.5 * (t(x) %*% A %*% x) + t(x) %*% x
  xv <- Variable(2)
  Ac <- Constant(rbind(c(1, 0), c(0, -1)))
  expr <- (1/2) * (t(xv) %*% Ac %*% xv) + t(xv) %*% xv
  expect_true(is(expr@args[[1]]@args[[2]], QuadForm))
  expect_true(identical(expr@args[[1]]@args[[2]]@args[[1]], x))
  
  # QuadForm inside nested expr: (0.5 * t(c) %*% c) * (t(x) %*% A %*% x) + t(x) %*% x
  xv <- Variable(2)
  Ac <- Constant(rbind(c(1, 0), c(0, -1)))
  cc <- Constant(c(2, -2))
  expr <- (1/2 * t(cc) %*% cc) * (t(xv) %*% Ac %*% xv) + t(xv) %*% xv
  expect_true(is(expr@args[[1]]@args[[2]], QuadForm))
  expect_true(identical(expr@args[[1]]@args[[2]]@args[[1]], x))
  
  # QuadForm with sparse matrices
  xv <- Variable(2)
  Ac <- Constant(Matrix(diag(2), sparse = TRUE))
  expr <- t(xv) %*% Ac %*% xv
  expect_true(is(expr, QuadForm))
  
  # QuadForm with mismatched dimensions raises error
  xv <- Variable(2)
  Ac <- Constant(diag(3))
  expect_error(t(xv) %*% Ac %*% xv)
  
  # QuadForm with PSD-wrapped matrix
  xv <- Variable(2)
  Ac <- Constant(rbind(c(1, 0), c(0, 1)))
  expr <- t(xv) %*% psd_wrap(Ac) %*% xv
  expect_true(is(expr, QuadForm))
  
  # QuadForm with nested subexpr
  xv <- Variable(2)
  Ac <- Constant(rbind(c(2, 0, 0), c(0, 0, 1)))
  M <- Constant(rbind(c(2, 0, 0), c(0, 2, 0), c(0, 0, 2)))
  b <- Constant(c(1, 2, 3))
  
  yv <- Ac %*% xv - b
  expr <- t(yv) %*% M %*% yv
  expect_true(is(expr, QuadForm))
  expect_true(identical(expr@args[[1]], yv))
  expect_true(identical(expr@args[[2]], M))
  
  # QuadForm with parameters
  xv <- Variable(2)
  Ap <- Parameter(2, 2, symmetric = TRUE)
  expr <- t(xv) %*% Ap %*% xv
  expect_true(is(expr, QuadForm))
  
  # Expect error for asymmetric/nonhermitian matrices
  xv <- Variable(2)
  Ac <- Constant(rbind(c(1, 0), c(1, 1)))
  expect_error(t(xv) %*% Ac %*% xv)
  
  xv <- Variable(2)
  Ac <- Constant(rbind(c(1, 1i), c(1i, 1)))
  expect_error(t(xv) %*% Ac %*% xv)
  
  # Not a quad_form because t(x) %*% A %*% y where x, y not necessarily equal
  xv <- Variable(2)
  yv <- Variable(2)
  Ac <- Constant(rbind(c(1, 0), c(0, -1)))
  expr <- t(xv) %*% Ac %*% yv
  expect_false(is(expr, QuadForm))
  
  # Not a quad_form because M is variable
  xv <- Variable(2)
  M <- Variable(2, 2)
  expr <- t(xv) %*% M %*% xv
  expect_false(is(expr, QuadForm))
  
  xc <- Constant(c(1, 0))
  M <- Variable(2, 2)
  expr <- t(xc) %*% M %*% xc
  expect_false(is(expr, QuadForm))
})
