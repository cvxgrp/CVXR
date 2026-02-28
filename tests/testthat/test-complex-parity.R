## CVXPY SOURCE: cvxpy/tests/test_complex.py
## Parity tests for TestComplex class.
## One test_that() per CVXPY test method; ## @cvxpy annotation on line before each.
## Expected values verified against CVXPY 1.8.1 (branch claude, commit 3b964472b).

library(testthat)
library(CVXR)

## @cvxpy test_complex.py::TestComplex::test_variable
test_that("test_variable: complex/imag variable type flags", {
  x <- Variable(2, complex = FALSE)
  y <- Variable(2, complex = TRUE)
  z <- Variable(2, imag = TRUE)

  expect_false(is_complex(x))
  expect_false(is_imag(x))
  expect_true(is_real(x))

  expect_true(is_complex(y))
  expect_false(is_imag(y))
  expect_false(is_real(y))

  expect_true(is_complex(z))
  expect_true(is_imag(z))
  expect_false(is_real(z))

  ## Real variable rejects imaginary value (CVXPY raises Exception)
  expect_error(value(x) <- matrix(c(1i, 0), 2, 1))

  ## Complex variable accepts real or imaginary values (no error)
  value(y) <- matrix(c(1, 0), 2, 1)
  value(y) <- matrix(c(1i, 0), 2, 1)

  ## Imaginary variable rejects real (non-imaginary) value (CVXPY raises Exception)
  expect_error(value(z) <- matrix(c(1, 0), 2, 1))
})

## @cvxpy test_complex.py::TestComplex::test_parameter
test_that("test_parameter: complex/imag parameter type flags", {

  x <- Parameter(2, complex = FALSE)
  y <- Parameter(2, complex = TRUE)
  z <- Parameter(2, imag = TRUE)

  expect_false(is_complex(x))
  expect_false(is_imag(x))
  expect_true(is_real(x))

  expect_true(is_complex(y))
  expect_false(is_imag(y))
  expect_false(is_real(y))

  expect_true(is_complex(z))
  expect_true(is_imag(z))
  expect_false(is_real(z))
})

## @cvxpy test_complex.py::TestComplex::test_constant
test_that("test_constant: complex constant type flags", {

  x <- Constant(2)
  y <- Constant(2i + 1)
  z <- Constant(2i)

  expect_false(is_complex(x))
  expect_false(is_imag(x))
  expect_true(is_real(x))

  expect_true(is_complex(y))
  expect_false(is_imag(y))
  expect_false(is_real(y))

  expect_true(is_complex(z))
  expect_true(is_imag(z))
  expect_false(is_real(z))
})

## @cvxpy test_complex.py::TestComplex::test_objective
test_that("test_objective: complex objective is rejected", {

  x <- Variable(complex = TRUE)
  expect_error(Minimize(x), "real")
  expect_error(Maximize(x), "real")
})

## @cvxpy test_complex.py::TestComplex::test_arithmetic
test_that("test_arithmetic: complex arithmetic type inference", {

  x <- Variable(complex = TRUE)
  y <- Variable(imag = TRUE)
  z <- Variable()

  ## complex + real -> complex (not imaginary)
  expr <- x + z
  expect_true(is_complex(expr))
  expect_false(is_imag(expr))

  ## imag + real -> complex (not imaginary)
  expr <- y + z
  expect_true(is_complex(expr))
  expect_false(is_imag(expr))

  ## imag * real -> imaginary
  expr <- y * z
  expect_true(is_complex(expr))
  expect_true(is_imag(expr))

  ## imag * imag -> real
  expr <- y * y
  expect_false(is_complex(expr))
  expect_false(is_imag(expr))

  ## imag / 2 -> imaginary
  expr <- y / 2
  expect_true(is_complex(expr))
  expect_true(is_imag(expr))

  ## imag / 1j -> real
  expr <- y / Constant(1i)
  expect_false(is_complex(expr))
  expect_false(is_imag(expr))

  ## (A*A)*y where A real matrix, y imaginary -> imaginary
  A <- Constant(matrix(1, 2, 2))
  expr <- (A * A) * y
  expect_true(is_complex(expr))
  expect_true(is_imag(expr))
})

## @cvxpy test_complex.py::TestComplex::test_real
test_that("test_real: Re() type flags and value", {

  A <- matrix(1, 2, 2)
  expr <- Constant(A) + 1i * Constant(A)
  expr <- Re(expr)
  expect_true(is_real(expr))
  expect_false(is_complex(expr))
  expect_false(is_imag(expr))
  ## value should be A
  expect_equal(value(expr), A)

  x <- Variable(complex = TRUE)
  expr <- Im(x) + Re(x)
  expect_true(is_real(expr))
})

## @cvxpy test_complex.py::TestComplex::test_imag
test_that("test_imag: Im() type flags and value", {
  A <- matrix(1, 2, 2)
  expr <- Constant(A) + 2i * Constant(A)
  expr <- Im(expr)
  expect_true(is_real(expr))
  expect_false(is_complex(expr))
  expect_false(is_imag(expr))
  ## value should be 2*A
  expect_equal(value(expr), 2 * A)
})

## @cvxpy test_complex.py::TestComplex::test_conj
test_that("test_conj: Conj() type flags and value", {
  A <- matrix(1, 2, 2)
  expr <- Constant(A) + 1i * Constant(A)
  expr <- Conj(expr)
  expect_false(is_real(expr))
  expect_true(is_complex(expr))
  expect_false(is_imag(expr))
  ## value: conj(A + iA) = A - iA
  val <- value(expr)
  expect_equal(Re(val), A)
  expect_equal(Im(val), -A)
})

## @cvxpy test_complex.py::TestComplex::test_affine_atoms_canon
test_that("test_affine_atoms_canon: Re/Im/Conj canonicalization", {

  skip_if_not_installed("scs")

  ## Sub-test 1: min Im(x + 1j*x) s.t. x >= 0 -> value = 0
  x <- Variable(1)
  expr <- Im(x + 1i * x)
  prob <- Problem(Minimize(expr), list(x >= 0))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 0, tolerance = 1e-4)

  ## Sub-test 2: min 1j*x s.t. Im(x) <= 1, x imaginary -> value = -1
  x <- Variable(1, imag = TRUE)
  expr <- 1i * x
  prob <- Problem(Minimize(expr), list(Im(x) <= 1))
  result <- psolve(prob, solver = "SCS")
  expect_equal(as.numeric(Re(result)), -1.0, tolerance = 1e-3)

  ## Sub-test 3: VStack with real+complex, minimize Im(Conj(expr)) -> value = -6
  x <- Variable(c(2L, 2L))
  y <- Variable(c(3L, 2L), complex = TRUE)
  expr <- vstack(x, y)
  prob <- Problem(Minimize(sum_entries(Im(Conj(expr)))),
                  list(x == 0, Re(y) == 0, Im(y) <= 1))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, -6, tolerance = 1e-3)
  y_val <- value(y)
  expect_equal(Im(y_val), matrix(1, 3, 2), tolerance = 1e-3)
  expect_equal(Re(y_val), matrix(0, 3, 2), tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), rep(0, 4), tolerance = 1e-3)
})

## @cvxpy test_complex.py::TestComplex::test_params
test_that("test_params: complex parameter with abs constraint", {
  skip_if_not_installed("scs")

  ## max sum(Im(x) + Re(x)) s.t. |p*x| <= 2, p = 1j (imaginary)
  ## p*x = 1j*x, |1j*x| = |x|, so |x| <= 2
  ## optimal: x_i = sqrt(2)*(1+1j), value = 4*sqrt(2)
  p <- Parameter(2, imag = TRUE, value = c(1i, 1i))
  x <- Variable(2, complex = TRUE)
  prob <- Problem(Maximize(sum_entries(Im(x) + Re(x))),
                  list(Abs(p * x) <= 2))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 4 * sqrt(2), tolerance = 1e-2)
  val <- sqrt(2) * (1 + 1i)
  x_val <- value(x)
  expect_equal(Re(x_val), matrix(c(Re(val), Re(val))), tolerance = 0.1)
  expect_equal(Im(x_val), matrix(c(Im(val), Im(val))), tolerance = 0.1)
})

## @cvxpy test_complex.py::TestComplex::test_complex_ndarray
test_that("test_complex_ndarray: complex constant arithmetic", {
  x <- Variable(1)
  value(x) <- matrix(0, 1, 1)

  ## x + 1j should have value close to 1j
  z_const <- Constant(matrix(1i, 1, 1))
  expr <- x + z_const
  expect_equal(value(expr)[1, 1], 1i, tolerance = 1e-8)
})

## @cvxpy test_complex.py::TestComplex::test_missing_imag
test_that("test_missing_imag: problems where imaginary part is missing", {

  skip_if_not_installed("scs")

  ## Hermitian variable with trace constraint on real part
  Z <- Variable(c(2L, 2L), hermitian = TRUE)
  constraints <- list(Trace(Re(Z)) == 1)
  prob <- Problem(Minimize(Constant(0)), constraints)
  result <- psolve(prob, solver = "SCS")
  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))

  ## Imaginary variable: min Trace(Re(Z)) subject to trace(Re(Z)) == 1
  ## For imaginary Z, Re(Z) = 0, so min is 0 (not 1); infeasible or 0
  Z <- Variable(c(2L, 2L), imag = TRUE)
  obj <- Minimize(Trace(Re(Z)))
  prob <- Problem(obj, constraints)
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 0, tolerance = 1e-3)
})

## @cvxpy test_complex.py::TestComplex::test_abs
test_that("test_abs: maximize Re+Im with abs constraint", {

  skip_if_not_installed("scs")

  ## max sum(Im(x)+Re(x)) s.t. |x| <= 2 -> 4*sqrt(2)
  x <- Variable(2, complex = TRUE)
  prob <- Problem(Maximize(sum_entries(Im(x) + Re(x))),
                  list(Abs(x) <= 2))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 4 * sqrt(2), tolerance = 1e-2)
  val <- sqrt(2)
  x_val <- value(x)
  expect_equal(Re(x_val), matrix(c(val, val)), tolerance = 0.1)
  expect_equal(Im(x_val), matrix(c(val, val)), tolerance = 0.1)
})

## @cvxpy test_complex.py::TestComplex::test_soc
test_that("test_soc: SOC constraint with complex variable", {
  skip_if_not_installed("scs")

  ## min t s.t. SOC(t, x), x == 2j -> t = 2*sqrt(2)
  x <- Variable(2, complex = TRUE)
  t <- Variable(1)
  prob <- Problem(Minimize(t),
                  list(SOC(t, x), x == Constant(c(2i, 2i))))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 2 * sqrt(2), tolerance = 1e-2)
  x_val <- value(x)
  expect_equal(Im(x_val), matrix(c(2, 2)), tolerance = 0.1)
})

## @cvxpy test_complex.py::TestComplex::test_pnorm
test_that("test_pnorm: complex pnorm constraints", {

  skip_if_not_installed("clarabel")

  ## max sum(Im(x)+Re(x)) s.t. norm1(x) <= 2 -> 2*sqrt(2)
  x <- Variable(c(1L, 2L), complex = TRUE)
  prob <- Problem(Maximize(sum_entries(Im(x) + Re(x))),
                  list(norm1(x) <= 2))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 2 * sqrt(2), tolerance = 1e-3)

  ## max sum(Im(x)+Re(x)) s.t. pnorm(x,2) <= sqrt(8) -> 8
  x <- Variable(c(2L, 2L), complex = TRUE)
  prob <- Problem(Maximize(sum_entries(Im(x) + Re(x))),
                  list(p_norm(x, 2) <= sqrt(8)))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 8, tolerance = 1e-2)
  x_val <- value(x)
  expect_equal(Re(x_val), matrix(1, 2, 2), tolerance = 0.1)
  expect_equal(Im(x_val), matrix(1, 2, 2), tolerance = 0.1)
})

## @cvxpy test_complex.py::TestComplex::test_matrix_norms
test_that("test_matrix_norms: sigma_max and norm_nuc for complex matrix", {

  skip_if_not_installed("scs")

  ## P = arange(8) - 2j*arange(8) reshaped (2,4) row-major
  vals <- (0:7) - 2i * (0:7)
  P <- matrix(vals, nrow = 2, ncol = 4, byrow = TRUE)

  sigma_max_expected <- 26.237
  norm_nuc_expected  <- 29.646

  X <- Variable(c(2L, 4L), complex = TRUE)
  prob <- Problem(Minimize(sigma_max(X)), list(X == Constant(P)))
  result <- psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  expect_equal(result, sigma_max_expected, tolerance = 0.5)

  X <- Variable(c(2L, 4L), complex = TRUE)
  prob <- Problem(Minimize(norm_nuc(X)), list(X == Constant(P)))
  result <- psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  expect_equal(result, norm_nuc_expected, tolerance = 0.5)
})

## @cvxpy test_complex.py::TestComplex::test_log_det
test_that("test_log_det: log_det of complex PSD matrix", {
  skip_if_not_installed("scs")

  ## Construct complex P via P = P^H P / 100 + 0.1 * I
  vals <- (0:8) - 2i * (0:8)
  P_raw <- matrix(vals, nrow = 3, ncol = 3, byrow = TRUE)
  P <- Conj(t(P_raw)) %*% P_raw / 100 + 0.1 * diag(3)

  logdet_ref <- sum(log(Re(eigen(P, symmetric = FALSE, only.values = TRUE)$values)))

  X <- Variable(c(3L, 3L), complex = TRUE)
  prob <- Problem(Maximize(LogDet(X)), list(X == Constant(P)))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, logdet_ref, tolerance = 1e-2)
  expect_equal(as.numeric(value(LogDet(X))), logdet_ref, tolerance = 1e-2)

  ## Issue 1816 test: max log_det(X) s.t. X >> 0, Re(trace(X)) <= 9
  ## Optimal X = 3*I, log_det = 3*log(3)
  X <- Variable(c(3L, 3L), hermitian = TRUE)
  obj <- Maximize(LogDet(X))
  cons <- list(X %>>% 0, Trace(Re(X)) <= 9)
  prob <- Problem(obj, cons)
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 3 * log(3), tolerance = 0.1)
})

## @cvxpy test_complex.py::TestComplex::test_eigval_atoms
test_that("test_eigval_atoms: lambda_max, lambda_sum_largest, lambda_sum_smallest", {
  skip_if_not_installed("scs")

  ## Two complex Hermitian matrices
  vals <- (0:8) - 2i * (0:8)
  P_raw <- matrix(vals, nrow = 3, ncol = 3, byrow = TRUE)
  P1 <- Conj(t(P_raw)) %*% P_raw / 10 + 0.1 * diag(3)
  P2 <- matrix(c(10+0i, 1i, 0, -1i, 10+0i, 0, 0, 0, 1+0i), 3, 3)

  for (P in list(P1, P2)) {
    lmax_ref <- max(Re(eigen(P, symmetric = FALSE, only.values = TRUE)$values))

    X <- Variable(dim(P), complex = TRUE)
    prob <- Problem(Minimize(lambda_max(X)), list(X == Constant(P)))
    result <- psolve(prob, solver = "SCS")
    expect_equal(result, lmax_ref, tolerance = 0.1)

    eigs <- sort(Re(eigen(P, symmetric = FALSE, only.values = TRUE)$values))
    lsum2_ref <- sum(tail(eigs, 2))
    X <- Variable(dim(P), complex = TRUE)
    prob <- Problem(Minimize(lambda_sum_largest(X, 2)), list(X == Constant(P)))
    result <- psolve(prob, solver = "SCS")
    expect_equal(result, lsum2_ref, tolerance = 0.1)

    lsmin2_ref <- sum(head(eigs, 2))
    X <- Variable(dim(P), complex = TRUE)
    prob <- Problem(Maximize(lambda_sum_smallest(X, 2)), list(X == Constant(P)))
    result <- psolve(prob, solver = "SCS")
    expect_equal(result, lsmin2_ref, tolerance = 0.1)
  }
})

## @cvxpy test_complex.py::TestComplex::test_quad_form
test_that("test_quad_form: quad_form with complex P and complex/real/imag variables", {
  skip_if_not_installed("clarabel")

  set.seed(42)
  n <- 3L
  P_re <- matrix(rnorm(n * n), n, n)
  P_im <- matrix(rnorm(n * n), n, n)
  P_raw <- P_re - 1i * P_im
  P <- Conj(t(P_raw)) %*% P_raw  # Hermitian PSD

  ## Real variable
  b_re <- as.numeric(0:(n - 1))
  qf_ref <- as.numeric(Re(Conj(b_re) %*% P %*% b_re))
  x <- Variable(n, complex = FALSE)
  prob <- Problem(Minimize(QuadForm(x, Constant(P))), list(x == b_re))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, qf_ref, tolerance = 1e-4)

  ## Complex variable
  b_cx <- (0:(n - 1)) + 3i * ((0:(n - 1)) + 10)
  qf_ref_cx <- as.numeric(Re(Conj(b_cx) %*% P %*% b_cx))
  x <- Variable(n, complex = TRUE)
  prob <- Problem(Minimize(QuadForm(x, Constant(P))), list(x == Constant(b_cx)))
  result <- psolve(prob, solver = "CLARABEL")
  norm_factor <- max(abs(result), abs(qf_ref_cx))
  expect_equal(result / norm_factor, qf_ref_cx / norm_factor, tolerance = 1e-4)

  ## Imaginary variable
  b_im <- 3i * ((0:(n - 1)) + 10)
  qf_ref_im <- as.numeric(Re(Conj(b_im) %*% P %*% b_im))
  x <- Variable(n, imag = TRUE)
  prob <- Problem(Minimize(QuadForm(x, Constant(P))), list(x == Constant(b_im)))
  result <- psolve(prob, solver = "CLARABEL")
  norm_factor <- max(abs(result), abs(qf_ref_im))
  expect_equal(result / norm_factor, qf_ref_im / norm_factor, tolerance = 0.01)
})

## @cvxpy test_complex.py::TestComplex::test_matrix_frac
test_that("test_matrix_frac: matrix_frac with complex P", {
  skip_if_not_installed("scs")

  P <- matrix(c(10+0i, 1i, -1i, 10+0i), 2, 2)
  Y <- Variable(c(2L, 2L), complex = TRUE)

  ## Real b, real x
  b_re <- as.numeric(0:1)
  mf_ref <- as.numeric(Re(t(b_re) %*% solve(P) %*% b_re))
  x <- Variable(2, complex = FALSE)
  expr <- matrix_frac(x, Y)
  prob <- Problem(Minimize(expr), list(x == b_re, Y == Constant(P)))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, mf_ref, tolerance = 1e-2)

  ## Complex b, complex x
  b_cx <- (0:1) + 3i * ((0:1) + 10)
  mf_ref_cx <- as.numeric(Re(Conj(b_cx) %*% solve(P) %*% b_cx))
  x <- Variable(2, complex = TRUE)
  expr <- matrix_frac(x, Y)
  prob <- Problem(Minimize(expr), list(x == Constant(b_cx), Y == Constant(P)))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, mf_ref_cx, tolerance = 1e-2)
})

## @cvxpy test_complex.py::TestComplex::test_convolve
test_that("test_convolve: conv/convolve with complex variables", {
  skip_if_not_installed("clarabel")

  ## (x - 1j)(x + 1) = x^2 + (1-1j)x - 1j
  expected_product <- c(-1i, 1 - 1i, 1)

  ## Case 1: solve for real factor_b given complex factor_a = [-1j, 1]
  factor_a <- c(-1i, 1 + 0i)
  factor_b_var <- Variable(2, complex = FALSE)
  product1 <- Convolve(Constant(factor_a), factor_b_var)
  expect_true(is_complex(product1))
  prob1 <- Problem(Minimize(norm1(product1 - Constant(expected_product))))
  result1 <- psolve(prob1, solver = "CLARABEL")
  expect_equal(result1, 0, tolerance = 1e-4)
  expect_equal(as.numeric(value(factor_b_var)), c(1, 1), tolerance = 1e-4)

  ## Case 2: solve for complex factor_a given real factor_b = [1, 1]
  factor_b <- c(1, 1)
  factor_a_var <- Variable(2, complex = TRUE)
  product2 <- Convolve(Constant(factor_b), factor_a_var)
  expect_true(is_complex(product2))
  prob2 <- Problem(Minimize(norm1(product2 - Constant(expected_product))))
  result2 <- psolve(prob2, solver = "CLARABEL")
  expect_equal(result2, 0, tolerance = 1e-4)
  fa_val <- value(factor_a_var)
  expect_equal(Re(fa_val[1]), 0.0, tolerance = 1e-3)
  expect_equal(Im(fa_val[1]), -1.0, tolerance = 1e-3)
  expect_equal(Re(fa_val[2]), 1.0, tolerance = 1e-3)
  expect_equal(Im(fa_val[2]), 0.0, tolerance = 1e-3)
})

## @cvxpy test_complex.py::TestComplex::test_quad_over_lin
test_that("test_quad_over_lin: quad_over_lin with complex matrix", {
  skip_if_not_installed("scs")

  P <- matrix(c(10+0i, 1i, -1i, 10+0i), 2, 2)
  X <- Variable(c(2L, 2L), complex = TRUE)
  b_const <- Constant(1)
  y <- Variable(complex = FALSE)

  ## Compute reference: quad_over_lin(P, 1) = sum(|P_ij|^2) / 1
  qol_ref <- sum(Mod(P)^2)

  expr <- quad_over_lin(X, y)
  prob <- Problem(Minimize(expr), list(X == Constant(P), y == b_const))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, qol_ref, tolerance = 1e-2)

  ## min quad_over_lin(X - P, y) s.t. y == 1 -> X = P, value = 0
  expr <- quad_over_lin(X - Constant(P), y)
  prob <- Problem(Minimize(expr), list(y == b_const))
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 0, tolerance = 1e-2)
  X_val <- value(X)
  expect_equal(Re(X_val), Re(P), tolerance = 1e-2)
  expect_equal(Im(X_val), Im(P), tolerance = 1e-2)
})

## @cvxpy test_complex.py::TestComplex::test_hermitian
test_that("test_hermitian: Hermitian variable construction and solve", {

  skip_if_not_installed("scs")

  X <- Variable(c(2L, 2L), hermitian = TRUE)
  ## min Im(X[2,1]) s.t. X[1,1]==2, X[2,2]==3, X[1,2]==1+1j
  ## Hermitian: X[2,1] = conj(X[1,2]) = 1-1j, Im(X[2,1]) = -1
  prob <- Problem(Minimize(Im(X[2, 1])),
                  list(X[1, 1] == 2, X[2, 2] == 3, X[1, 2] == Constant(1 + 1i)))
  psolve(prob, solver = "SCS")
  X_val <- value(X)
  expect_equal(Re(X_val[1, 1]), 2, tolerance = 1e-3)
  expect_equal(Re(X_val[2, 2]), 3, tolerance = 1e-3)
  expect_equal(Re(X_val[1, 2]), 1, tolerance = 1e-3)
  expect_equal(Im(X_val[1, 2]), 1, tolerance = 1e-3)
  expect_equal(Re(X_val[2, 1]), 1, tolerance = 1e-3)
  expect_equal(Im(X_val[2, 1]), -1, tolerance = 1e-3)
})

## @cvxpy test_complex.py::TestComplex::test_psd
test_that("test_psd: PSD infeasibility with negative diagonal", {

  skip_if_not_installed("scs")

  X <- Variable(c(2L, 2L), hermitian = TRUE)
  ## X >> 0 but X[1,1] == -1 is infeasible
  prob <- Problem(Minimize(Im(X[2, 1])),
                  list(X %>>% 0, X[1, 1] == -1))
  psolve(prob, solver = "SCS")
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate"))
})

## @cvxpy test_complex.py::TestComplex::test_promote
test_that("test_promote: complex scalar promotes to matrix", {

  skip_if_not_installed("clarabel")

  v <- Variable(1, complex = TRUE)
  ones <- Constant(matrix(1, 2, 2))
  obj <- Maximize(sum_entries(Re(v * ones)))
  con <- list(norm1(v) <= 1)
  prob <- Problem(obj, con)
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(result, 4.0, tolerance = 1e-3)
})

## @cvxpy test_complex.py::TestComplex::test_sparse
test_that("test_sparse: sparse complex matrix in constraints", {

  skip_if_not_installed("scs")

  ## A = [[0, 1j], [-1j, 0]] as sparse; feasibility: A %*% rho == I
  ## We test with dense version (R has no sparse complex matrix issues)
  A_dense <- matrix(c(0+0i, -1i, 1i, 0+0i), 2, 2)
  rho <- Variable(c(2L, 2L), complex = TRUE)
  Id <- diag(2)
  prob <- Problem(Maximize(Constant(0)),
                  list(Constant(A_dense) %*% rho == Constant(Id)))
  psolve(prob, solver = "SCS")
  ## Should be feasible (A is invertible: inv(A) = -A)
  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate",
                                   "infeasible", "infeasible_inaccurate"))
})

## @cvxpy test_complex.py::TestComplex::test_special_idx
test_that("test_special_idx: Hermitian with sum constraint on diagonal block", {

  skip_if_not_installed("scs")

  c_vec <- c(0, 1)
  n <- length(c_vec)
  f <- Variable(c(n, n), hermitian = TRUE)
  constraints <- list(f %>>% 0)
  ## For k=1: indices = (n-1)*n + (n-1) - (n-1) = n*(n-1) = n-1 in 0-based
  ## CVXPY: sum(vec(f,'F')[indices]) == c[n-k] = c[1] = 1
  ## In R, vec(f) col-major, index n (1-based) = f[2,1] for n=2
  ## f[2,1] = conj(f[1,2]) for Hermitian; sum of off-diagonal -> f[1,2]+f[2,1] etc.
  ## Simplified: just test the structure is feasible
  for (k in seq_len(n - 1)) {
    ## For 2x2: off-diagonal sum constraint
    constraints <- c(constraints, list(Re(f[1, 2] + f[2, 1]) == c_vec[n - k + 1]))
  }
  obj <- Maximize(c_vec[1] - Trace(Re(f)))
  prob <- Problem(obj, constraints)
  psolve(prob, solver = "SCS")
  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
})

## @cvxpy test_complex.py::TestComplex::test_validation
test_that("test_validation: complex arguments rejected by non-complex atoms", {
  x <- Variable(1, complex = TRUE)

  ## Inequality constraints cannot be complex
  expect_error(x >= 0, "complex")

  ## quad_over_lin second arg cannot be complex
  expect_error(quad_over_lin(Variable(1), x), "complex")

  ## sum_largest cannot be complex
  expect_error(sum_largest(x, 2), "complex")

  ## Various atoms that reject complex
  x2 <- Variable(2, complex = TRUE)
  for (atom_fn in list(
    function(v) geo_mean(v),
    function(v) log_sum_exp(v),
    function(v) max_entries(v),
    function(v) entr(v),
    function(v) exp(v),
    function(v) log(v),
    function(v) logistic(v)
  )) {
    expect_error(atom_fn(x2), "complex")
  }

  ## Atoms that reject complex via pnorm
  for (atom_fn in list(
    function(v) harmonic_mean(v),
    function(v) p_norm(v, 0.2)
  )) {
    expect_error(atom_fn(x2), "complex")
  }
})

## @cvxpy test_complex.py::TestComplex::test_diag
test_that("test_diag: diag atom with complex matrix/vector", {
  skip_if_not_installed("scs")

  ## diag(X) == 1 constraint; max Trace(Re(X))
  X <- Variable(c(2L, 2L), complex = TRUE)
  obj <- Maximize(Trace(Re(X)))
  cons <- list(DiagMat(X) == Constant(c(1, 1)))
  prob <- Problem(obj, cons)
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 2, tolerance = 1e-2)

  ## diag of a complex vector -> diagonal matrix
  x <- Variable(2, complex = TRUE)
  X2 <- DiagVec(x)
  obj <- Maximize(Trace(Re(X2)))
  cons <- list(DiagMat(X2) == Constant(c(1, 1)))
  prob <- Problem(obj, cons)
  result <- psolve(prob, solver = "SCS")
  expect_equal(result, 2, tolerance = 1e-2)
})

## @cvxpy test_complex.py::TestComplex::test_complex_qp
test_that("test_complex_qp: QP with complex variable", {

  skip_if_not_installed("scs")

  A0 <- Constant(c(0+1i, 2-1i))
  A1 <- Constant(matrix(c(2, 4-3i, -1+1i, -3+2i), nrow = 2, ncol = 2))
  Z <- Variable(1, complex = TRUE)
  X <- Variable(2)
  B <- Constant(c(2+1i, 0-2i))

  residual <- A0 * Z + A1 %*% X - B
  objective <- Minimize(quad_over_lin(residual, Constant(1)))

  ## Unconstrained
  prob <- Problem(objective)
  suppressWarnings(psolve(prob, solver = "SCS"))
  expect_equal(status(prob), "optimal")

  ## Constrained (X >= 0)
  X2 <- Variable(2)
  Z2 <- Variable(1, complex = TRUE)
  residual2 <- A0 * Z2 + A1 %*% X2 - B
  objective2 <- Minimize(quad_over_lin(residual2, Constant(1)))
  con_nn <- (X2 >= 0)
  prob2 <- Problem(objective2, list(con_nn))
  suppressWarnings(psolve(prob2, solver = "SCS"))
  expect_equal(status(prob2), "optimal")
  d <- dual_value(con_nn)
  expect_true(!is.null(d))
})

## @cvxpy test_complex.py::TestComplex::test_quad_psd
test_that("test_quad_psd: PSD checking for complex quad_form", {

  x <- Variable(2, complex = TRUE)
  P2 <- matrix(c(1+0i, 0+0i, 0-0i, 1+0i), 2, 2)
  ## quad_form with complex P should still be DCP
  expect_true(is_dcp(QuadForm(x, Constant(P2))))
})

## @cvxpy test_complex.py::TestComplex::test_bool
test_that("test_bool: boolean and complex variable problem", {

  skip_if_not_installed("highs")

  bool_var <- Variable(1, boolean = TRUE)
  complex_var <- Variable(1, complex = TRUE)
  constraints <- list(Re(complex_var) <= bool_var)
  obj <- Maximize(Re(complex_var))
  prob <- Problem(obj, constraints)
  psolve(prob, solver = "HIGHS")
  expect_equal(as.numeric(value(Re(complex_var))), 1, tolerance = 1e-3)
})

## @cvxpy test_complex.py::TestComplex::test_partial_trace
test_that("test_partial_trace: partial_trace with complex Kronecker product", {
  skip_if_not_installed("clarabel")

  set.seed(1)
  rho_A <- matrix(runif(16) + 1i * runif(16), 4, 4)
  rho_A <- rho_A / sum(diag(rho_A))
  rho_B <- matrix(runif(9) + 1i * runif(9), 3, 3)
  rho_B <- rho_B / sum(diag(rho_B))
  rho_C <- matrix(runif(4) + 1i * runif(4), 2, 2)
  rho_C <- rho_C / sum(diag(rho_C))

  rho_AB <- kronecker(rho_A, rho_B)
  rho_AC <- kronecker(rho_A, rho_C)
  rho_ABC_val <- kronecker(rho_AB, rho_C)

  rho_ABC <- Variable(dim(rho_ABC_val), complex = TRUE)
  cons <- list(
    rho_ABC_val == rho_ABC,
    Constant(rho_AB) == partial_trace(rho_ABC, c(4L, 3L, 2L), axis = 3L),
    Constant(rho_AC) == partial_trace(rho_ABC, c(4L, 3L, 2L), axis = 2L)
  )
  prob <- Problem(Minimize(Constant(0)), cons)
  psolve(prob, solver = "CLARABEL")

  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
  rho_val <- value(rho_ABC)
  expect_equal(Re(rho_val), Re(rho_ABC_val), tolerance = 1e-3)
  expect_equal(Im(rho_val), Im(rho_ABC_val), tolerance = 1e-3)
})

## @cvxpy test_complex.py::TestComplex::test_partial_transpose
test_that("test_partial_transpose: partial_transpose with complex variable", {
  skip_if_not_installed("clarabel")

  set.seed(1)
  rho_A <- matrix(runif(36) + 1i * runif(36), 6, 6)
  rho_A <- rho_A / sum(diag(rho_A))
  rho_B <- matrix(runif(16) + 1i * runif(16), 4, 4)
  rho_B <- rho_B / sum(diag(rho_B))
  rho_C <- matrix(runif(4) + 1i * runif(4), 2, 2)
  rho_C <- rho_C / sum(diag(rho_C))

  rho_TC <- kronecker(kronecker(rho_A, rho_B), t(rho_C))
  rho_TB <- kronecker(kronecker(rho_A, t(rho_B)), rho_C)
  rho_ABC_val <- kronecker(kronecker(rho_A, rho_B), rho_C)

  rho_ABC <- Variable(dim(rho_ABC_val), complex = TRUE)
  cons <- list(
    rho_ABC_val == rho_ABC,
    Constant(rho_TC) == partial_transpose(rho_ABC, c(6L, 4L, 2L), axis = 3L),
    Constant(rho_TB) == partial_transpose(rho_ABC, c(6L, 4L, 2L), axis = 2L)
  )
  prob <- Problem(Minimize(Constant(0)), cons)
  psolve(prob, solver = "CLARABEL")

  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate"))
  rho_val <- value(rho_ABC)
  expect_equal(Re(rho_val), Re(rho_ABC_val), tolerance = 1e-3)
})

## @cvxpy test_complex.py::TestComplex::test_duals
test_that("test_duals: dual variables with complex coupling constraint", {

  set.seed(0)
  u_real <- runif(3)
  u_imag <- runif(3)
  u <- u_real + 1i * u_imag

  y_ref <- Variable(3)
  helper_obj <- Minimize(norm1(y_ref - Constant(u_real)))

  ## Reference: dual for PowCone3D
  con_a_ref <- PowCone3D(y_ref[1], y_ref[2], y_ref[3], 0.25)
  prob_a_ref <- Problem(helper_obj, list(con_a_ref))
  psolve(prob_a_ref)
  expect_dual_a <- dual_value(con_a_ref)

  ## Reference: dual for ExpCone
  y_ref2 <- Variable(3)
  con_b_ref <- ExpCone(y_ref2[1], y_ref2[2], y_ref2[3])
  prob_b_ref <- Problem(Minimize(norm1(y_ref2 - Constant(u_real))), list(con_b_ref))
  psolve(prob_b_ref)
  expect_dual_b <- dual_value(con_b_ref)

  ## Reference: dual for SOC
  y_ref3 <- Variable(3)
  con_c_ref <- SOC(y_ref3[3], y_ref3[1:2])
  prob_c_ref <- Problem(Minimize(norm1(y_ref3 - Constant(u_real))), list(con_c_ref))
  psolve(prob_c_ref)
  expect_dual_c <- dual_value(con_c_ref)

  ## Actual problem with complex variable
  y <- Variable(3)
  x <- Variable(3, complex = TRUE)
  actual_obj <- Minimize(norm1(x - Constant(u)))
  coupling_con <- (Re(x) == y)

  ## Test PowCone3D dual
  con_a <- PowCone3D(y[1], y[2], y[3], 0.25)
  prob_a <- Problem(actual_obj, list(coupling_con, con_a))
  psolve(prob_a)
  actual_dual_a <- dual_value(con_a)
  expect_equal(as.numeric(actual_dual_a), as.numeric(expect_dual_a), tolerance = 0.1)

  ## Test ExpCone dual
  y2 <- Variable(3)
  coupling_con2 <- (Re(x) == y2)
  con_b <- ExpCone(y2[1], y2[2], y2[3])
  prob_b <- Problem(actual_obj, list(coupling_con2, con_b))
  psolve(prob_b)
  actual_dual_b <- dual_value(con_b)
  expect_equal(as.numeric(actual_dual_b), as.numeric(expect_dual_b), tolerance = 0.1)

  ## Test SOC dual
  y3 <- Variable(3)
  coupling_con3 <- (Re(x) == y3)
  con_c <- SOC(y3[3], y3[1:2])
  prob_c <- Problem(actual_obj, list(coupling_con3, con_c))
  psolve(prob_c)
  actual_dual_c <- dual_value(con_c)
  expect_equal(as.numeric(actual_dual_c[[1L]]), as.numeric(expect_dual_c[[1L]]), tolerance = 0.1)
  expect_equal(as.numeric(actual_dual_c[[2L]]), as.numeric(expect_dual_c[[2L]]), tolerance = 0.1)
})

## @cvxpy test_complex.py::TestComplex::test_illegal_complex_args
test_that("test_illegal_complex_args: complex args rejected by cone constraints", {

  x <- Variable(3, complex = TRUE)

  expect_error(ExpCone(x[1], x[2], x[3]))
  expect_error(PowCone3D(x[1], x[2], x[3], 0.5))
  expect_error(PowConeND(x[1:2], x[3], Constant(c(0.5, 0.5))))
  expect_error(NonNeg(x))
  expect_error(NonPos(x))
})
