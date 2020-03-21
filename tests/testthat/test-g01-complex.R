context("test-g01-complex")
TOL <- 1e-6

test_that("test the Variable class", {
  skip_on_cran()
  x <- Variable(2, complex = FALSE)
  y <- Variable(2, complex = TRUE)
  z <- Variable(2, imag = TRUE)

  expect_false(is_complex(x))
  expect_false(is_imag(x))
  expect_true(is_complex(y))
  expect_false(is_imag(y))
  expect_true(is_complex(z))
  expect_true(is_imag(z))

  expect_error(value(x) <- c(1i, 0))

  value(y) <- c(1, 0)
  value(y) <- c(1i, 0)

  expect_error(value(z) <- c(1, 0))
})

test_that("test the Parameter class", {
  skip_on_cran()
  x <- Parameter(2, complex = FALSE)
  y <- Parameter(2, complex = TRUE)
  z <- Parameter(2, imag = TRUE)

  expect_false(is_complex(x))
  expect_false(is_imag(x))
  expect_true(is_complex(y))
  expect_false(is_imag(y))
  expect_true(is_complex(z))
  expect_true(is_imag(z))

  expect_error(value(x) <- c(1i, 0))

  value(y) <- c(1, 0)
  value(y) <- c(1i, 0)

  expect_error(value(z) <- c(1, 0))
})

test_that("test the Constant class", {
  skip_on_cran()
  x <- Constant(2)
  y <- Constant(2i + 1)
  z <- Constant(2i)

  expect_false(is_complex(x))
  expect_false(is_imag(x))
  expect_true(is_complex(y))
  expect_false(is_imag(y))
  expect_true(is_complex(z))
  expect_true(is_imag(z))
})

test_that("test objectives", {
  skip_on_cran()
  x <- Variable(complex = TRUE)
  expect_error(Minimize(x))
  expect_error(Maximize(x))
})

test_that("test basic arithmetic expressions", {
  skip_on_cran()
  x <- Variable(complex = TRUE)
  y <- Variable(imag = TRUE)
  z <- Variable()

  expr <- x + z
  expect_true(is_complex(expr))
  expect_false(is_imag(expr))

  expr <- y + z
  expect_true(is_complex(expr))
  expect_false(is_imag(expr))

  expr <- y*z
  expect_true(is_complex(expr))
  expect_true(is_imag(expr))

  expr <- y*y
  expect_false(is_complex(expr))
  expect_false(is_imag(expr))

  expr <- y/2
  expect_true(is_complex(expr))
  expect_true(is_imag(expr))

  expr <- y/1i
  expect_false(is_complex(expr))
  expect_false(is_imag(expr))

  A <- matrix(1, nrow = 2, ncol = 2)
  expr <- A %*% y %*% A
  expect_true(is_complex(expr))
  expect_true(is_imag(expr))
})

test_that("test real function", {
  skip_on_cran()
  A <- matrix(1, nrow = 2, ncol = 2)
  expr <- Constant(A) + 1i*Constant(A)
  expr <- Re(expr)
  expect_true(is_real(expr))
  expect_false(is_complex(expr))
  expect_false(is_imag(expr))
  expect_equal(value(expr), A)

  x <- Variable(complex = TRUE)
  expr <- Im(x) + Re(x)
  expect_true(is_real(expr))
})

test_that("test imag function", {
  skip_on_cran()
  A <- matrix(1, nrow = 2, ncol = 2)
  expr <- Constant(A) + 2i*Constant(A)
  expr <- Im(expr)
  expect_true(is_real(expr))
  expect_false(is_complex(expr))
  expect_false(is_imag(expr))
  expect_equal(value(expr), 2*A)
})

test_that("test conj function", {
  skip_on_cran()
  A <- matrix(1, nrow = 2, ncol = 2)
  expr <- Constant(A) + 1i*Constant(A)
  expr <- Conj(expr)
  expect_false(is_real(expr))
  expect_true(is_complex(expr))
  expect_false(is_imag(expr))
  expect_equal(value(expr), A - 1i*A)
})

test_that("test canonicalization for affine atoms", {
  skip_on_cran()
  # Scalar.
  x <- Variable()
  expr <- Im(x + 1i*x)
  prob <- Problem(Minimize(expr), list(x >= 0))
  result <- solve(prob)
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$getValue(x), 0, tolerance = TOL)

  x <- Variable(imag = TRUE)
  expr <- 1i*x
  prob <- Problem(Minimize(expr), list(Im(x) <= 1))
  result <- solve(prob)
  expect_equal(result$value, -1, tolerance = TOL)
  expect_equal(result$getValue(x), 1i, tolerance = TOL)

  x <- Variable(2)
  expr <- x/1i
  prob <- Problem(Minimize(expr[1]*1i + expr[2]*1i), list(Re(x + 1i) >= 1))
  result <- solve(prob)
  expect_equal(result$value, -Inf)
  prob <- Problem(Minimize(expr[1]*1i + expr[2]*1i), list(Re(x + 1i) <= 1))
  result <- solve(prob)
  expect_equal(result$value, -2, tolerance = TOL)
  expect_equal(result$getValue(x), as.matrix(c(1, 1)), tolerance = TOL)
  prob <- Problem(Minimize(expr[1]*1i + expr[2]*1i), list(Re(x + 1i) >= 1, Conj(x) <= 0))
  result <- solve(prob)   # TODO_NARAS_1: OSQP returns a solver_error, but ECOS is correct.
  expect_equal(result$value, Inf)

  x <- Variable(2,2)
  y <- Variable(3,2, complex = TRUE)
  expr <- vstack(x, y)
  prob <- Problem(Minimize(sum(Im(Conj(expr)))), list(x == 0, Re(y) == 0, Im(y) <= 1))
  result <- solve(prob)
  expect_equal(result$value, -6, tolerance = TOL)
  expect_equal(result$getValue(y), 1i*matrix(1, nrow = 3, ncol = 2), tolerance = TOL)
  expect_equal(result$getValue(x), matrix(0, nrow = 2, ncol = 2), tolerance = TOL)
})

# test_that("test with parameters", {
#  p <- Parameter(imag = TRUE, value = 1i)
#  x <- Variable(2, complex = TRUE)
#  prob <- Problem(Maximize(sum(Im(x) + Re(x))), list(abs(p*x) <= 2))
#  result <- solve(prob)
#  expect_equal(result$value, 4*sqrt(2), tolerance = TOL)
#  val <- matrix(sqrt(2), nrow = 2, ncol = 1)
#  expect_equal(result$getValue(x), val + 1i*val)
# })

test_that("test with absolute value", {
  skip_on_cran()
  x <- Variable(2, complex = TRUE)
  prob <- Problem(Maximize(sum(Im(x) + Re(x))), list(abs(x) <= 2))
  result <- solve(prob)
  expect_equal(result$value, 4*sqrt(2))
  val <- matrix(sqrt(2), nrow = 2, ncol = 1)
  expect_equal(result$getValue(x), val + 1i*val)
})

test_that("test with p-norm", {
  skip_on_cran()
  x <- Variable(1, 2, complex = TRUE)
  prob <- Problem(Maximize(sum(Im(x) + Re(x))), list(norm1(x) <= 2))
  result <- solve(prob)
  expect_equal(result$value, 2*sqrt(2), tolerance = TOL)
  val <- matrix(sqrt(2)/2, nrow = 2, ncol = 1)
  # expect_equal(result$getValue(x), val + 1i*val)

  x <- Variable(2, 2, complex = TRUE)
  prob <- Problem(Maximize(sum(Im(x) + Re(x))), list(p_norm(x, p = 2) <= sqrt(8)))
  result <- solve(prob)
  expect_equal(result$value, 8, tolerance = TOL)
  val <- matrix(1, nrow = 2, ncol = 2)
  expect_equal(result$getValue(x), val + 1i*val, tolerance = TOL)
})

test_that("test matrix norms", {
  skip_on_cran()
  P <- 0:7 - 2i*(0:7)
  P <- matrix(P, nrow = 2, ncol = 4, byrow = TRUE)
  sigma_max <- base:::norm(P, type = "2")
  X <- Variable(2, 4, complex = TRUE)
  prob <- Problem(Minimize(norm(X, "2")), list(X == P))
  result <- solve(prob)
  expect_equal(result$value, sigma_max, tolerance = 1e-3)

  norm_nuc_val <- sum(svd(P)$d)
  X <- Variable(2, 4, complex = TRUE)
  prob <- Problem(Minimize(norm_nuc(X)), list(X == P))
  result <- solve(prob, solver = "SCS", eps = 1e-4)
  expect_equal(result$value, norm_nuc_val, tolerance = 0.1)
})

test_that("test log-determinant", {
  skip_on_cran()
  P <- (0:8) - 2i*(0:8)
  P <- matrix(P, nrow = 3, ncol = 3)
  P <- Conj(t(P)) %*% P/100 + diag(0.1, 3)
  value <- value(log_det(P))
  X <- Variable(3, 3, complex = TRUE)
  prob <- Problem(Maximize(log_det(X)), list(X == P))
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, value, tolerance = 1e-2)
})

test_that("test eigenvalue atoms", {
  skip_on_cran()
  P <- (0:8) - 2i*(0:8)
  P <- matrix(P, nrow = 3, ncol = 3)
  P1 <- Conj(t(P)) %*% P/10 + diag(0.1, 3)
  P2 <- rbind(c(10, 1i, 0), c(-1i, 10, 0), c(0, 0, 1))
  for(P in list(P1, P2)) {
    value <- value(lambda_max(P))
    # X <- Variable(dim(P), complex = TRUE)
    X <- Variable(nrow(P), ncol(P), complex = TRUE)
    prob <- Problem(Minimize(lambda_max(X)), list(X == P))
    result <- solve(prob, solver = "SCS", eps = 1e-5)
    expect_equal(result$value, value, tolerance = 1e-2)

    eigs <- Re(eigen(P, only.values = TRUE)$values)
    value <- value(sum_largest(eigs, 2))
    # X <- Variable(dim(P), complex = TRUE)
    X <- Variable(nrow(P), ncol(P), complex = TRUE)
    prob <- Problem(Minimize(lambda_sum_largest(X, 2)), list(X == P))
    result <- solve(prob, solver = "SCS", eps = 1e-8, verbose = TRUE)
    expect_equal(result$value, value, tolerance = 1e-3)
    expect_equal(result$getValue(X), P, tolerance = 1e-3)

    value <- value(sum_smallest(eigs, 2))
    # X <- Variable(dim(P), complex = TRUE)
    X <- Variable(nrow(P), ncol(P), complex = TRUE)
    prob <- Problem(Maximize(lambda_sum_smallest(X, 2)), list(X == P))
    result <- solve(prob, solver = "SCS", eps = 1e-6)
    expect_equal(result$value, value, tolerance = 1e-3)
  }
})

test_that("test quad_form atom", {
  skip_on_cran()
  # Create a random positive definite Hermitian matrix for all tests.
  set.seed(42)
  P <- matrix(rnorm(9), nrow = 3, ncol = 3) - 1i*matrix(rnorm(9), nrow = 3, ncol = 3)
  P <- Conj(t(P)) %*% P

  # Solve a problem with real variable.
  b <- 0:2
  x <- Variable(3, complex = FALSE)
  value <- value(quad_form(b, P))
  prob <- Problem(Minimize(quad_form(x, P)), list(x == b))
  result <- solve(prob)
  expect_equal(result$value, Re(value[1]))

  # Solve a problem with complex variable.
  b <- (0:2) + 3i*(0:2 + 10)
  x <- Variable(3, complex = TRUE)
  value <- value(quad_form(b, P))
  prob <- Problem(Minimize(quad_form(x, P)), list(x == b))
  result <- solve(prob)
  normalization <- max(abs(result$value), abs(value))
  expect_equal(result$value/normalization, Re(value[1])/normalization, tolerance = 1e-5)

  # Solve a problem with imaginary variable.
  b <- 3i*(0:2 + 10)
  x <- Variable(3, imag = TRUE)
  value <- value(quad_form(b, P))
  expr <- quad_form(x, P)
  prob <- Problem(Minimize(expr), list(x == b))
  result <- solve(prob)
  normalization <- max(abs(result$value), abs(value))
  expect_equal(result$value/normalization, Re(value[1])/normalization, tolerance = TOL)
})

test_that("test matrix_frac atom", {
  skip_on_cran()
  P <- rbind(c(10, 1i), c(-1i, 10))
  Y <- Variable(2, 2, complex = TRUE)
  b <- 0:1
  x <- Variable(2, complex = FALSE)
  value <- value(matrix_frac(b, P))
  expr <- matrix_frac(x, Y)
  prob <- Problem(Minimize(expr), list(x == b, Y == P))
  result <- solve(prob, solver = "SCS", eps = 1e-6, max_iters = 7500, verbose = TRUE)
  expect_equal(result$value, Re(value[1]), tolerance = 1e-3)

  b <- (0:1 + 3i*(0:1 + 10))
  x <- Variable(2, complex = TRUE)
  value <- value(matrix_frac(b, P))
  expr <- matrix_frac(x, Y)
  prob <- Problem(Minimize(expr), list(x == b, Y == P))
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, Re(value[1]), tolerance = 1e-3)

  b <- (0:1 + 10)/10i
  x <- Variable(2, imag = TRUE)
  value <- value(matrix_frac(b, P))
  expr <- matrix_frac(x, Y)
  prob <- Problem(Minimize(expr), list(x == b, Y == P))
  result <- solve(prob, solver = "SCS", eps = 1e-5, max_iters = 7500)
  expect_equal(result$value, Re(value[1]), tolerance = 1e-3)
})

test_that("test Hermitian variables", {
  skip_on_cran()
  X <- Variable(2, 2, hermitian = TRUE)
  prob <- Problem(Minimize(Im(X[2,1])), list(X[1,1] == 2, X[2,2] == 3, X[1,2] == 1+1i))
  result <- solve(prob)
  expect_equal(result$getValue(X), cbind(c(2, 1-1i), c(1+1i, 3)))
})

test_that("test positive semidefinite variables", {
  skip_on_cran()
  X <- Variable(2, 2, hermitian = TRUE)
  prob <- Problem(Minimize(Im(X[2,1])), list(X %>>% 0, X[1,1] == -1))
  result <- solve(prob)
  expect_equal(result$status, "infeasible")
})

test_that("test promotion of complex variables", {
  skip_on_cran()
  v <- Variable(complex = TRUE)
  obj <- Maximize(Re(sum(v * matrix(1, nrow = 2, ncol = 2))))
  con <- list(cvxr_norm(v) <= 1)
  prob <- Problem(obj, con)
  result <- solve(prob)
  expect_equal(result$value, 4.0, tolerance = TOL)
})

# TODO_NARAS_2: Figure out how to handle complex sparse matrices in R.
# test_that("test problem with complex sparse matrix", {
#   # Define sparse matrix [[0, 1i], [-1i, 0]]
#   require(Matrix)
#   A <- rbind(c(0, 1i), c(-1i, 0))
#   A_sparse <- Matrix(A, sparse = TRUE)
#
#   # Feasibility with sparse matrix.
#   rho <- Variable(2, 2, complex = TRUE)
#   Id <- diag(2)
#   obj <- Maximize(0)
#   cons <- list(A_sparse %*% rho == Id)
#   prob <- Problem(obj, cons)
#   result <- solve(prob)
#   rho_sparse <- value(rho)
#   # Infeasible here, which is wrong!
#
#   # Feasibility with R matrices.
#   rho <- Variable(2, 2, complex = TRUE)
#   Id <- diag(2)
#   obj <- Maximize(0)
#   cons <- list(A %*% rho == Id)
#   prob <- Problem(obj, cons)
#   result <- solve(prob)
#   expect_equal(result$getValue(rho), rho_sparse)
# })

test_that("test with special index", {
  skip_on_cran()
  c <- c(0, 1)
  n <- length(c)

  # Create optimization variables.
  f <- Variable(n, n, hermitian = TRUE)

  # Create constraints.
  constraints <- list(f %>>% 0)
  for(k in seq_len(n-1)) {   # TODO: Check indices match range in Python.
    i <- seq(from = n-k, length.out = k)
    indices <- (i*n) + i - (n-k) + 1
    constraints <- c(constraints, sum(vec(f)[indices]) == c[n-k+1])
  }

  # Form objective.
  obj <- Maximize(c[1] - Re(matrix_trace(f)))

  # Form and solve problem.
  prob <- Problem(obj, constraints)
  sol <- solve(prob)
  expect_equal(sol$status, "optimal")
})

test_that("test that complex arguments are rejected", {
  skip_on_cran()
  x <- Variable(complex = TRUE)
  expect_error(x >= 0, "Inequality constraints cannot be complex.")
  expect_error(quad_over_lin(x, x), "The second argument to QuadOverLin cannot be complex.")
  expect_error(sum_largest(x, 2), "Arguments to SumLargest cannot be complex.")

  x <- Variable(2, complex = TRUE)
  for(atom in c("GeoMean", "LogSumExp", "MaxEntries", "Entr", "Exp", "Huber", "Log", "Log1p", "Logistic")) {
    print(atom)
    error_msg <- paste("Arguments to ", atom, " cannot be complex.", sep = "")

    if(atom %in% c("LogSumExp", "MaxEntries"))
      expect_error(new(atom, expr = x), error_msg)
    else
      expect_error(new(atom, x = x), error_msg)
  }

  x <- Variable(2, complex = TRUE)
  for(atom in c("MaxElemwise", "KLDiv")) {
    print(atom)
    error_msg <- paste("Arguments to ", atom, " cannot be complex.", sep = "")

    if(atom == "MaxElemwise")
      expect_error(max_elemwise(x, x), error_msg)
    else
      expect_error(kl_div(x, x), error_msg)
  }

  x <- Variable(2, complex = TRUE)
  for(atom in c(inv_pos, sqrt, function(y) { power(y, 0.2) })) {
    expect_error(atom(x), "Arguments to Power cannot be complex.")
  }

  x <- Variable(2, complex = TRUE)
  for(atom in c(harmonic_mean, function(y) { p_norm(y, 0.2) })) {
    expect_error(atom(x), "Pnorm(x, p) cannot have x complex for p < 1.", fixed = TRUE)
  }
})
