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

  expect_error(value(x) <- c(1i, 0), "Variable value must be real", fixed = TRUE)

  value(y) <- c(1, 0)
  value(y) <- c(1i, 0)

  expect_error(value(z) <- c(1, 0), "Variable value must be imaginary", fixed = TRUE)
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

  expect_error(value(x) <- c(1i, 0), "Parameter value must be real", fixed = TRUE)

  value(y) <- c(1, 0)
  value(y) <- c(1i, 0)

  expect_error(value(z) <- c(1, 0), "Parameter value must be imaginary", fixed = TRUE)
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
  expect_error(Minimize(x), "The Minimize objective must be real valued", fixed = TRUE)
  expect_error(Maximize(x), "The Maximize objective must be real valued", fixed = TRUE)
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
  expr <- (A %*% A) %*% y
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
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 0, tolerance = TOL)
  expect_equal(result$getValue(x), 0, tolerance = TOL)

  x <- Variable(imag = TRUE)
  expr <- 1i*x
  prob <- Problem(Minimize(expr), list(Im(x) <= 1))
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, -1, tolerance = TOL)
  expect_equal(result$getValue(x), 1i, tolerance = TOL)

  x <- Variable(2)
  expr <- x/1i
  prob <- Problem(Minimize(expr[1]*1i + expr[2]*1i), list(Re(x + 1i) >= 1))
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, -Inf)
  prob <- Problem(Minimize(expr[1]*1i + expr[2]*1i), list(Re(x + 1i) <= 1))
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, -2, tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(1, 1)), tolerance = TOL)
  prob <- Problem(Minimize(expr[1]*1i + expr[2]*1i), list(Re(x + 1i) >= 1, Conj(x) <= 0))
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, Inf)

  x <- Variable(2, 2)
  y <- Variable(3, 2, complex = TRUE)
  expr <- vstack(x, y)
  prob <- Problem(Minimize(sum(Im(Conj(expr)))), list(x == 0, Re(y) == 0, Im(y) <= 1))
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, -6, tolerance = TOL)
  expect_equal(result$getValue(y), 1i*matrix(1, nrow = 3, ncol = 2), tolerance = TOL)
  expect_equal(result$getValue(x), matrix(0, nrow = 2, ncol = 2), tolerance = TOL)
  
  x <- Variable(2, 2)
  y <- Variable(3, 2, complex = TRUE)
  expr <- vstack(x, y)
  prob <- Problem(Minimize(sum(Im(H(expr)))), list(x == 0, Re(y) == 0, Im(y) <= 1))
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, -6, tolerance = TOL)
  expect_equal(result$getValue(y), 1i*matrix(1, nrow = 3, ncol = 2), tolerance = TOL)
  expect_equal(result$getValue(x), matrix(0, nrow = 2, ncol = 2), tolerance = TOL)
})

test_that("test with parameters", {
  p <- Parameter(imag = TRUE, value = 1i)
  x <- Variable(2, complex = TRUE)
  prob <- Problem(Maximize(sum(Im(x) + Re(x))), list(abs(p*x) <= 2))
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 4*sqrt(2), tolerance = TOL)
  val <- matrix(sqrt(2), nrow = 2, ncol = 1)
  expect_equal(result$getValue(x), val + 1i*val)
})

test_that("test problems where imaginary is missing", {
  Z <- Variable(2, 2, hermitian = TRUE)
  constraints <- list(matrix_trace(Re(Z)) == 1)
  obj <- Minimize(0)
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  
  Z <- Variable(2, 2, imag = TRUE)
  obj <- Minimize(matrix_trace(Re(Z)))
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, 0, tolerance = TOL)
})

test_that("test with absolute value", {
  skip_on_cran()
  x <- Variable(2, complex = TRUE)
  prob <- Problem(Maximize(sum(Im(x) + Re(x))), list(abs(x) <= 2))
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 4*sqrt(2), tolerance = TOL)
  val <- matrix(sqrt(2), nrow = 2, ncol = 1)
  expect_equal(result$getValue(x), val + 1i*val, tolerance = TOL)
})

test_that("test with SOC", {
  x <- Variable(2, complex = TRUE)
  t <- Variable()
  prob <- Problem(Minimize(t), list(SOC(t, x), x == 2i))
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, 2*sqrt(2), tolerance = TOL)
  expect_equal(result$getValue(x), matrix(c(2i, 2i)), tolerance = TOL)
})

test_that("test complex with p-norm", {
  skip_on_cran()
  x <- Variable(1, 2, complex = TRUE)
  prob <- Problem(Maximize(sum(Im(x) + Re(x))), list(norm1(x) <= 2))
  result <- solve(prob, solver = "ECOS")
  expect_equal(result$value, 2*sqrt(2), tolerance = TOL)
  val <- matrix(sqrt(2)/2, nrow = 2, ncol = 1)
  # expect_equal(result$getValue(x), val + 1i*val, tolerance = TOL)

  x <- Variable(2, 2, complex = TRUE)
  prob <- Problem(Maximize(sum(Im(x) + Re(x))), list(p_norm(x, p = 2) <= sqrt(8)))
  result <- solve(prob, solver = "ECOS")
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
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, sigma_max, tolerance = 1e-1)

  norm_nuc_val <- sum(svd(P)$d)
  X <- Variable(2, 4, complex = TRUE)
  prob <- Problem(Minimize(norm_nuc(X)), list(X == P))
  result <- solve(prob, solver = "SCS", eps = 1e-4)
  expect_equal(result$value, norm_nuc_val, tolerance = 1e-1)
})

test_that("test log-determinant", {
  skip_on_cran()
  P <- (0:8) - 2i*(0:8)
  P <- matrix(P, nrow = 3, ncol = 3, byrow = TRUE)
  P <- Conj(t(P)) %*% P/100 + diag(0.1, 3)
  logdet_value <- value(log_det(P))
  X <- Variable(3, 3, complex = TRUE)
  objective <- log_det(X)
  prob <- Problem(Maximize(objective), list(X == P))
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, logdet_value, tolerance = 1e-2)
  objective_value <- result$getValue(objective)
  ld <- as.numeric(determinant(P, logarithm = TRUE)$modulus)
  expect_equal(objective_value, ld, tolerance = 1e-3)
  
  # Test case for Issue 1816 from CVXPY.
  # The optimal solution is the identity matrix scaled by 3.
  cons <- list(X %>>% 0, Re(matrix_trace(X)) <= 9)
  obj <- log_det(X)
  prob <- Problem(Maximize(obj), cons)
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$getValue(obj), 3*log(3), tolerance = TOL)
})

test_that("test eigenvalue atoms", {
  skip_on_cran()
  P <- (0:8) - 2i*(0:8)
  P <- matrix(P, nrow = 3, ncol = 3, byrow = TRUE)
  P1 <- Conj(t(P)) %*% P/10 + diag(0.1, 3)
  P2 <- rbind(c(10, 1i, 0), c(-1i, 10, 0), c(0, 0, 1))
  for(P in list(P1, P2)) {
    val <- value(lambda_max(P))
    X <- new("Variable", dim = dim(P), complex = TRUE)
    # X <- Variable(nrow(P), ncol(P), complex = TRUE)
    prob <- Problem(Minimize(lambda_max(X)), list(X == P))
    result <- solve(prob, solver = "SCS", eps = 1e-5)
    expect_equal(result$value, val, tolerance = 1e-2)

    eigs <- Re(eigen(P, only.values = TRUE)$values)
    val <- value(sum_largest(eigs, 2))
    X <- new("Variable", dim = dim(P), complex = TRUE)
    # X <- Variable(nrow(P), ncol(P), complex = TRUE)
    prob <- Problem(Minimize(lambda_sum_largest(X, 2)), list(X == P))
    result <- solve(prob, solver = "SCS", eps = 1e-8)
    expect_equal(result$value, val, tolerance = 1e-3)
    expect_equal(result$getValue(X), P, tolerance = 1e-3)

    val <- value(sum_smallest(eigs, 2))
    X <- new("Variable", dim = dim(P), complex = TRUE)
    # X <- Variable(nrow(P), ncol(P), complex = TRUE)
    prob <- Problem(Maximize(lambda_sum_smallest(X, 2)), list(X == P))
    result <- solve(prob, solver = "SCS", eps = 1e-6)
    expect_equal(result$value, val, tolerance = 1e-3)
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
  val <- value(quad_form(b, P))
  prob <- Problem(Minimize(quad_form(x, P)), list(x == b))
  result <- solve(prob, solver = "ECOS")
  expect_equal(result$value, val, tolerance = TOL)

  # Solve a problem with complex variable.
  b <- (0:2) + 3i*(0:2 + 10)
  x <- Variable(3, complex = TRUE)
  val <- value(quad_form(b, P))
  prob <- Problem(Minimize(quad_form(x, P)), list(x == b))
  result <- solve(prob, solver = "ECOS")
  normalization <- max(abs(result$value), abs(val))
  expect_equal(result$value/normalization, val/normalization, tolerance = 1e-5)

  # Solve a problem with imaginary variable.
  b <- 3i*(0:2 + 10)
  x <- Variable(3, imag = TRUE)
  val <- value(quad_form(b, P))
  expr <- quad_form(x, P)
  prob <- Problem(Minimize(expr), list(x == b))
  result <- solve(prob, solver = "ECOS")
  normalization <- max(abs(result$value), abs(value))
  expect_equal(result$value/normalization, val/normalization, tolerance = TOL)
})

test_that("test matrix_frac atom", {
  skip_on_cran()
  P <- rbind(c(10, 1i), c(-1i, 10))
  Y <- Variable(2, 2, complex = TRUE)
  b <- 0:1
  x <- Variable(2, complex = FALSE)
  val <- value(matrix_frac(b, P))
  expr <- matrix_frac(x, Y)
  prob <- Problem(Minimize(expr), list(x == b, Y == P))
  result <- solve(prob, solver = "SCS", eps = 1e-6, max_iters = 7500, verbose = TRUE)
  expect_equal(result$value, val, tolerance = 1e-3)

  b <- (0:1 + 3i*(0:1 + 10))
  x <- Variable(2, complex = TRUE)
  val <- value(matrix_frac(b, P))
  expr <- matrix_frac(x, Y)
  prob <- Problem(Minimize(expr), list(x == b, Y == P))
  result <- solve(prob, solver = "SCS", eps = 1e-6)
  expect_equal(result$value, val, tolerance = 1e-3)

  b <- (0:1 + 10)/10i
  x <- Variable(2, imag = TRUE)
  val <- value(matrix_frac(b, P))
  expr <- matrix_frac(x, Y)
  prob <- Problem(Minimize(expr), list(x == b, Y == P))
  result <- solve(prob, solver = "SCS", eps = 1e-5, max_iters = 7500)
  expect_equal(result$value, val, tolerance = 1e-3)
})

test_that("test quad_over_lin atom", {
  P <- rbind(c(10, 1i), c(-1i, 10))
  X <- Variable(2, 2, complex = TRUE)
  b <- 1
  y <- Variable(complex = FALSE)
  
  val <- value(quad_over_lin(P, b))
  expr <- quad_over_lin(X, y)
  prob <- Problem(Minimize(expr), list(X == P, y == b))
  result <- solve(prob, solver = "SCS", eps = 1e-6, max_iters = 7500, verbose = TRUE)
  expect_equal(result$value, val, tolerance = 1e-3)
  
  expr <- quad_over_lin(X - P, y)
  prob <- Problem(Minimize(expr), list(y == b))
  result <- solve(prob, solver = "SCS", eps = 1e-6, max_iters = 7500, verbose = TRUE)
  expect_equal(result$value, 0, tolerance = 1e-3)
  expect_equal(result$getValue(X), P, tolerance = 1e-3)
})

test_that("test Hermitian variables", {
  skip_on_cran()
  X <- Variable(2, 2, hermitian = TRUE)
  prob <- Problem(Minimize(Im(X[2,1])), list(X[1,1] == 2, X[2,2] == 3, X[1,2] == 1+1i))
  result <- solve(prob, solver =)
  expect_equal(result$getValue(X), cbind(c(2, 1-1i), c(1+1i, 3)))
})

test_that("test positive semidefinite variables", {
  skip_on_cran()
  X <- Variable(2, 2, hermitian = TRUE)
  prob <- Problem(Minimize(Im(X[2,1])), list(X %>>% 0, X[1,1] == -1))
  result <- solve(prob, solver = "SCS")
  expect_equal(result$status, INFEASIBLE)
})

test_that("test promotion of complex variables", {
  skip_on_cran()
  v <- Variable(complex = TRUE)
  obj <- Maximize(Re(sum(v*matrix(1, nrow = 2, ncol = 2))))
  con <- list(cvxr_norm(v) <= 1)
  prob <- Problem(obj, con)
  result <- solve(prob, solver = "ECOS")
  expect_equal(result$value, 4.0, tolerance = TOL)
})

test_that("test problem with complex sparse matrix", {
  # TODO: Currently, Matrix library does NOT support complex sparse matrices.
  skip_on_cran()
  
  # Define sparse matrix [[0, 1i], [-1i, 0]]
  require(Matrix)
  row <- c(1, 2)
  col <- c(2, 1)
  data <- c(1i, -1i)
  A <- sparseMatrix(i = row, j = col, x = data, dims = c(2, 2), repr = "C")
  # A <- rbind(c(0, 1i), c(-1i, 0))

  # Feasibility with sparse matrix.
  rho <- Variable(2, 2, complex = TRUE)
  Id <- diag(2)
  obj <- Maximize(0)
  cons <- list(A %*% rho == Id)
  prob <- Problem(obj, cons)
  result <- solve(prob, solver = "SCS")
  rho_sparse <- result$getValue(rho)
  # Infeasible here, which is wrong!

  # Feasibility with R matrices.
  rho <- Variable(2, 2, complex = TRUE)
  Id <- diag(2)
  obj <- Maximize(0)
  cons <- list(as.matrix(A) %*% rho == Id)
  prob <- Problem(obj, cons)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$getValue(rho), rho_sparse, tolerance = TOL)
})

test_that("test with special index", {
  skip_on_cran()
  c <- c(0, 1)
  n <- length(c)

  # Create optimization variables.
  f <- Variable(n, n, hermitian = TRUE)

  # Create constraints.
  constraints <- list(f %>>% 0)
  for(k in seq(1, n-1)) {
    i <- seq(from = n-k, length.out = k)
    indices <- (i*n) + i - (n-k) + 1
    constraints <- c(constraints, sum(vec(f)[indices]) == c[n-k+1])
  }

  # Form objective.
  obj <- Maximize(c[1] - Re(matrix_trace(f)))

  # Form and solve problem.
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "SCS")
})

test_that("test that complex arguments are rejected", {
  skip_on_cran()
  x <- Variable(complex = TRUE)
  expect_error(x >= 0, "Inequality constraints cannot be complex.", fixed = TRUE)
  expect_error(quad_over_lin(x, x), "The second argument to QuadOverLin cannot be complex.", fixed = TRUE)
  expect_error(sum_largest(x, 2), "Arguments to SumLargest cannot be complex.", fixed = TRUE)
  expect_error(dotsort(x, 2), "Arguments to DotSort cannot be complex.", fixed = TRUE)
  expect_error(dotsort(Variable(2), c(1+2i)), "Arguments to DotSort cannot be complex.", fixed = TRUE)

  x <- Variable(2, complex = TRUE)
  for(atom in c("GeoMean", "LogSumExp", "MaxEntries", "Entr", "Exp", "Huber", "Log", "Log1p", "Logistic")) {
    print(atom)
    error_msg <- paste("Arguments to", atom, "cannot be complex.")

    if(atom %in% c("LogSumExp", "MaxEntries"))
      expect_error(new(atom, expr = x), error_msg, fixed = TRUE)
    else
      expect_error(new(atom, x = x), error_msg, fixed = TRUE)
  }

  x <- Variable(2, complex = TRUE)
  for(atom in c("MaxElemwise", "KLDiv")) {
    print(atom)
    error_msg <- paste("Arguments to", atom, "cannot be complex.")

    if(atom == "MaxElemwise")
      expect_error(max_elemwise(x, x), error_msg, fixed = TRUE)
    else
      expect_error(kl_div(x, x), error_msg, fixed = TRUE)
  }

  x <- Variable(2, complex = TRUE)
  for(atom in c(inv_pos, sqrt, function(y) { power(y, 0.2) }))
    expect_error(atom(x), "Arguments to Power cannot be complex.", fixed = TRUE)

  x <- Variable(2, complex = TRUE)
  for(atom in c(harmonic_mean, function(y) { p_norm(y, 0.2) }))
    expect_error(atom(x), "Pnorm(x, p) cannot have x complex for p < 1.", fixed = TRUE)
})

test_that("test diag of mat and of vector", {
  X <- Variable(2, 2, complex = TRUE)
  obj <- Maximize(matrix_trace(Re(X)))
  cons <- list(diag(X) == 1)
  prob <- Problem(obj, cons)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, 2, tolerance = TOL)
  
  x <- Variable(2, complex = TRUE)
  X <- diag(x)
  obj <- Maximize(matrix_trace(Re(X)))
  cons <- list(diag(X) == 1)
  prob <- Problem(obj, cons)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$value, 2, tolerance = TOL)
})

test_that("test a QP with a complex variable", {
  A0 <- matrix(c(0+1i, 2-1i))
  A1 <- rbind(c(2, -1+1i), c(4-3i, -3+2i))
  Z <- Variable(complex = TRUE)
  X <- Variable(2)
  B <- matrix(c(2+1i, -2i))
  
  objective <- Minimize(sum_squares(A0*Z + A1 %*% X - B))
  prob <- Problem(objective)
  result <- solve(prob, solver = "SCS")
  expect_equal(result$status, OPTIMAL)
  expect_false(is.na(result$getDualValue(constraints[[1]])))
})

test_that("test PSD checking from CVXPY Issue #1491", {
  x <- Variable(2, complex = TRUE)
  P1 <- diag(2)
  P2 <- rbind(c(1+0i, 0+0i), c(0-0i, 1+0i))
  print("P1 is real:", curvature(quad_form(x, P1)))
  print("P2 is complex:", curvature(quad_form(x, P2)))
  expect_true(is_dcp(quad_form(x, P2)))
})

test_that("test bool", {
  # The purpose of this test is to make sure that we don't try to recover dual
  # variables unless they're actually present. See CVXPY issue #1133.
  
  bool_var <- Variable(boolean = TRUE)
  complex_var <- Variable(complex = TRUE)
  
  constraints <- list(Re(complex_var) <= bool_var)
  obj <- Maximize(Re(complex_var))
  prob <- Problem(obj, constraints)
  result <- solve(prob, solver = "ECOS_BB")
  expect_equal(result$value, 1, tolerance = 1e-4)
})

test_that("test partial trace", {
  # Test a problem with partial_trace.
  # rho_ABC = rho_A \otimes rho_B \otimes rho_C.
  # Here \otimes signifies Kronecker product.
  # Each rho_i is normalized, i.e., Tr(rho_i) = 1.
  
  # Set random state.
  set.seed(1)
  
  # Generate test case.
  rho_A <- matrix(runif(4*4), nrow = 4, ncol = 4) + 1i*matrix(runif(4*4), nrow = 4, ncol = 4)
  rho_A <- rho_A/sum(diag(rho_A))
  rho_B <- matrix(runif(3*3), nrow = 3, ncol = 3) + 1i*matrix(runif(3*3), nrow = 3, ncol = 3)
  rho_B <- rho_B/sum(diag(rho_B))
  rho_C <- matrix(runif(2*2), nrow = 2, ncol = 2) + 1i*matrix(runif(2*2), nrow = 2, ncol = 2)
  rho_C <- rho_C/sum(diag(rho_C))
  rho_AB <- kronecker(rho_A, rho_B)
  rho_AC <- kronecker(rho_A, rho_C)
  
  # Construct a CVXR Variable with value equal to rho_A \otimes rho_B \otimes rho_C.
  rho_ABC_val <- kronecker(rho_AB, rho_C)
  rho_ABC <- new("Variable", dim = dim(rho_ABC_val), complex = TRUE)
  cons <- list(rho_ABC_val == rho_ABC, 
               rho_AB == partial_trace(rho_ABC, c(4, 3, 2), axis = 2),   # TODO: Check axis parameter handling matches CVXPY.
               rho_AC == partial_trace(rho_ABC, c(4, 3, 2), axis = 1))
  prob <- Problem(Minimize(0), cons)
  result <- solve(prob)
  
  print(rho_ABC_val)
  expect_true(is.allclose(result$getValue(rho_ABC), rho_ABC_val))
})

test_that("test partial transpose", {
  # Test a problem with partial_transpose.
  # rho_ABC = rho_A \otimes rho_B \otimes rho_C.
  # Here \otimes signifies Kronecker product.
  # Each rho_i is normalized, i.e. Tr(rho_i) = 1.
  
  # Set random state.
  set.seed(1)
  
  # Generate three test cases.
  rho_A <- matrix(runif(8*8), nrow = 8, ncol = 8) + 1i*matrix(runif(8*8), nrow = 8, ncol = 8)
  rho_A <- rho_A/sum(diag(rho_A))
  rho_B <- matrix(runif(6*6), nrow = 6, ncol = 6) + 1i*matrix(runif(6*6), nrow = 6, ncol = 6)
  rho_B <- rho_B/sum(diag(rho_B))
  rho_C <- matrix(runif(4*4), nrow = 4, ncol = 4) + 1i*matrix(runif(4*4), nrow = 4, ncol = 4)
  rho_C <- rho_C/sum(diag(rho_C))
  
  rho_TC <- kronecker(kronecker(rho_A, rho_B), t(rho_C))
  rho_TB <- kronecker(kronecker(rho_A, t(rho_B)), rho_C)
  
  # Construct a CVXR variable with value equal to rho_A \otimes rho_B \otimes rho_C.
  rho_ABC_val <- kronecker(kronecker(rho_A, rho_B), rho_C)
  rho_ABC <- new("Variable", dim = dim(rho_ABC_val), complex = TRUE)
  cons <- list(rho_ABC_val == rho_ABC,
               rho_TC <- partial_transpose(rho_ABC, c(8, 6, 4), axis = 2),   # TODO: Check axis parameter handling matches CVXPY.
               rho_TC <- partial_transpose(rho_ABC, c(8, 6, 4), axis = 1))
  prob <- Problem(Minimize(0), cons)
  result <- solve(prob)
  
  print(rho_ABC_val)
  expect_true(is.allclose(result$getValue(rho_ABC)), rho_ABC_val)
})

test_that("test duals", {
  set.seed(0)
  u_real <- runif(3)
  u_imag <- runif(3)
  u <- u_real + 1i*u_imag
  
  y <- Variable(3)
  helper_objective <- Minimize(norm(y - u_real))
  
  ####################################
  #
  #   Compute reference values
  #
  ####################################
  con_a <- PowCone3D(y[1], y[2], y[3], c(0.25))
  helper_prob_a <- Problem(helper_objective, list(con_a))
  result <- solve(helper_prob_a)
  expect_dual_a <- result$getDualValue(con_a)
  
  con_b <- ExpCone(y[1], y[2], y[3])
  helper_prob_b <- Problem(helper_objective, list(con_b))
  result <- solve(helper_prob_b)
  expect_dual_b <- result$getDualValue(con_b)
  
  con_c <- SOC(y[3], y[1:2])
  helper_prob_c <- Problem(helper_objective, list(con_c))
  result <- solve(helper_prob_c)
  expect_dual_c <- result$getDualValue(con_c)
  
  ####################################
  #
  #   Run tests
  #
  ####################################
  x <- Variable(3, complex = TRUE)
  actual_objective <- Minimize(cvxr_norm(x - u))
  coupling_con <- Re(x) == y
  
  con_a_test <- con_a
  # con_a_test <- copy(con_a)
  prob_a <- Problem(actual_objective, list(coupling_con, con_a_test))
  result <- solve(prob_a)
  actual_dual_a <- result$getDualValue(con_a_test)
  expect_equal(actual_dual_a, expect_dual_a, tolerance = 1e-2)
  
  con_b_test <- con_b
  # con_b_test <- copy(con_b)
  prob_b <- Problem(actual_objective, list(coupling_con, con_b_test))
  result <- solve(prob_b)
  actual_dual_b <- result$getDualValue(con_b_test)
  expect_equal(actual_dual_b, expect_dual_b, tolerance = 1e-2)
  
  con_c_test <- con_c
  # con_c_test <- copy(con_c)
  prob_c <- Problem(actual_objective, list(coupling_con, con_c_test))
  result <- solve(prob_c)
  actual_dual_c <- result$getDualValue(con_c_test)
  expect_equal(actual_dual_c[1], expect_dual_c[1], tolerance = 1e-2)
  expect_equal(actual_dual_c[2], expect_dual_c[2], tolerance = 1e-2)
})

test_that("test illegal complex args", {
  x <- Variable(3, complex = TRUE)
  expect_error(ExpCone(x[1], x[2], x[3]))
  expect_error(PowCone3D(x[1], x[2], c(0.5)))
  expect_error(PowConeND(x[1], x[1:2], c(0.5, 0.5)))
  expect_error(RelEntrConeQuad(x[1], x[2], x[3], 1, 1))
  expect_error(NonNeg(x))
  expect_error(NonPos(x))
})











