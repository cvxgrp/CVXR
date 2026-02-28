## Tests for Phase C3: End-to-end complex problem solving
## These tests solve complex problems via the full pipeline.
## CVXPY SOURCE: tests/test_complex.py

# ── Variable/Constant/Objective validation ─────────────────────────

## @cvxpy NONE
test_that("Complex variable type flags are correct", {
  x <- Variable(2)
  expect_true(is_real(x))
  expect_false(is_complex(x))
  expect_false(is_imag(x))

  z <- Variable(2, complex = TRUE)
  expect_false(is_real(z))
  expect_true(is_complex(z))
  expect_false(is_imag(z))

  w <- Variable(2, imag = TRUE)
  expect_false(is_real(w))
  expect_true(is_complex(w))
  expect_true(is_imag(w))
})

## @cvxpy NONE
test_that("Complex constant type flags are correct", {
  expect_false(is_complex(Constant(1)))
  expect_true(is_complex(Constant(1i)))
  expect_true(is_complex(Constant(1 + 1i)))
  expect_true(is_imag(Constant(2i)))
  expect_false(is_imag(Constant(1 + 2i)))
})

## @cvxpy NONE
test_that("Objective rejects complex expression", {
  z <- Variable(2, complex = TRUE)
  expect_error(Minimize(sum_entries(z)), "real")
})

# ── Real/Imag/Conj atoms ──────────────────────────────────────────

## @cvxpy NONE
test_that("Re/Im/Conj dispatch on CVXR expressions", {
  z <- Variable(2, complex = TRUE)
  expect_true(S7_inherits(Re(z), Real_))
  expect_true(S7_inherits(Im(z), Imag_))
  expect_true(S7_inherits(Conj(z), Conj_))
})

# ── Simple complex variable solve ─────────────────────────────────

## @cvxpy NONE
test_that("Complex least squares: min sum|x - b|", {
  skip_if_not_installed("clarabel")
  x <- Variable(2, complex = TRUE)
  b <- Constant(c(1 + 1i, 2 - 1i))
  prob <- Problem(Minimize(sum_entries(Abs(x - b))))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  x_val <- value(x)
  expect_equal(Re(x_val), Re(matrix(c(1 + 1i, 2 - 1i))), tolerance = 1e-4)
  expect_equal(Im(x_val), Im(matrix(c(1 + 1i, 2 - 1i))), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Minimize real part of complex variable", {
  skip_if_not_installed("clarabel")
  z <- Variable(2, complex = TRUE)
  prob <- Problem(Minimize(sum_entries(Re(z))),
                  list(Im(z) == 1, Re(z) >= -5))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(result, -10, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Imaginary variable solve", {
  skip_if_not_installed("clarabel")
  z <- Variable(2, imag = TRUE)
  prob <- Problem(Minimize(sum_entries(-Im(z))),
                  list(Im(z) <= 3))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(result, -6, tolerance = 1e-4)
})

# ── Complex affine atoms via solver ────────────────────────────────

## @cvxpy NONE
test_that("Complex variable with affine constraints", {
  skip_if_not_installed("clarabel")
  z <- Variable(2, complex = TRUE)
  ## min sum(Im(z)) s.t. sum(Re(z)) == 3, sum(Im(z)) == 4, Im(z) >= 0
  prob <- Problem(Minimize(sum_entries(Im(z))),
                  list(sum_entries(Re(z)) == 3,
                       sum_entries(Im(z)) == 4,
                       Im(z) >= 0))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 4, tolerance = 1e-4)
})

# ── Complex abs ───────────────────────────────────────────────────

## @cvxpy NONE
test_that("Maximize with complex abs constraint", {
  skip_if_not_installed("clarabel")
  ## CVXPY test_abs: maximize sum(Re(x) + Im(x)) s.t. |x| <= 2
  ## Optimal: x_i = sqrt(2)*(1+1i), value = 4*sqrt(2)
  x <- Variable(2, complex = TRUE)
  prob <- Problem(Maximize(sum_entries(Re(x) + Im(x))),
                  list(Abs(x) <= 2))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 4 * sqrt(2), tolerance = 1e-3)
})

# ── Complex norms ─────────────────────────────────────────────────

## @cvxpy NONE
test_that("Norm1 with complex variable", {
  skip_if_not_installed("clarabel")
  ## maximize sum(Re(x) + Im(x)) s.t. ||x||_1 <= 2
  x <- Variable(2, complex = TRUE)
  prob <- Problem(Maximize(sum_entries(Re(x) + Im(x))),
                  list(Norm1(x) <= 2))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 2 * sqrt(2), tolerance = 1e-3)
})

# ── Hermitian variable ────────────────────────────────────────────

## @cvxpy NONE
test_that("Hermitian variable construction and solve", {
  skip_if_not_installed("clarabel")
  n <- 2L
  X <- Variable(c(n, n), hermitian = TRUE)
  ## Set entries: X[1,1]=2, X[2,2]=3, X[1,2]=1+1i
  prob <- Problem(Minimize(0),
                  list(X[1, 1] == 2,
                       X[2, 2] == 3,
                       Re(X[1, 2]) == 1,
                       Im(X[1, 2]) == 1))
  psolve(prob)
  expect_equal(status(prob), "optimal")
  X_val <- value(X)
  ## X[2,1] should be conjugate of X[1,2] = 1-1i
  expect_equal(Re(X_val[1, 1]), 2, tolerance = 1e-4)
  expect_equal(Re(X_val[2, 2]), 3, tolerance = 1e-4)
  expect_equal(Re(X_val[1, 2]), 1, tolerance = 1e-4)
  expect_equal(Im(X_val[1, 2]), 1, tolerance = 1e-4)
  expect_equal(Re(X_val[2, 1]), 1, tolerance = 1e-4)
  expect_equal(Im(X_val[2, 1]), -1, tolerance = 1e-4)
})

# ── PSD constraint on Hermitian ───────────────────────────────────

## @cvxpy NONE
test_that("PSD infeasibility for non-PSD Hermitian", {
  skip_if_not_installed("clarabel")
  n <- 2L
  X <- Variable(c(n, n), hermitian = TRUE)
  ## X[1,1] = -1 with X >> 0 should be infeasible
  prob <- Problem(Minimize(0),
                  list(X %>>% 0, X[1, 1] == -1))
  psolve(prob)
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate"))
})

# ── Complex quad_form ─────────────────────────────────────────────

## @cvxpy NONE
test_that("Complex quad_form minimization", {
  skip_if_not_installed("clarabel")
  n <- 2L
  P <- matrix(c(4, 1+1i, 1-1i, 5), n, n)  # Hermitian PSD
  x <- Variable(n, complex = TRUE)
  prob <- Problem(Minimize(QuadForm(x, P)),
                  list(Re(x) >= 1))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  ## Value should be positive since P is PSD and x is constrained
  expect_true(result > 0)
})

# ── Complex equality dual recovery ────────────────────────────────

## @cvxpy NONE
test_that("Dual variables recovered for complex equality", {
  skip_if_not_installed("clarabel")
  z <- Variable(2, complex = TRUE)
  con <- (z == Constant(c(1 + 1i, 2 + 2i)))
  prob <- Problem(Minimize(sum_entries(Re(z) + Im(z))), list(con))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 6, tolerance = 1e-4)
  ## Dual should be non-null
  d <- dual_value(con)
  expect_true(!is.null(d))
})

# ── Log det with Hermitian variable ───────────────────────────────

## @cvxpy NONE
test_that("Maximize log_det of Hermitian PSD variable", {
  skip_if_not_installed("scs")
  n <- 2L
  X <- Variable(c(n, n), hermitian = TRUE)
  prob <- Problem(Maximize(LogDet(X)),
                  list(X %>>% 0, Trace(Re(X)) <= 4))
  result <- psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  ## Optimal: X = 2*I, log_det = log(4) = 2*log(2)
  expect_equal(result, 2 * log(2), tolerance = 0.1)
})

# ── SOC constraint with complex variable ──────────────────────────

## @cvxpy NONE
test_that("SOC constraint with imaginary variable", {
  skip_if_not_installed("clarabel")
  ## minimize t s.t. SOC(t, Im(x)), Im(x) == 2
  x <- Variable(2, imag = TRUE)
  t_var <- Variable(1)
  prob <- Problem(Minimize(t_var),
                  list(SOC(t_var, Im(x)), Im(x) == 2))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  ## ||[2, 2]||_2 = 2*sqrt(2)
  expect_equal(result, 2 * sqrt(2), tolerance = 1e-3)
})

# ── Mixed real/complex constraints ────────────────────────────────

## @cvxpy NONE
test_that("Problem with both real and complex variables", {
  skip_if_not_installed("clarabel")
  x <- Variable(1)
  z <- Variable(1, complex = TRUE)
  prob <- Problem(Minimize(x + Re(z)),
                  list(x >= 1, Im(z) == 2, Re(z) >= 3))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 4, tolerance = 1e-4)
})

# ── Chain insertion order verified ────────────────────────────────

## @cvxpy NONE
test_that("Complex2Real appears in solving chain", {
  z <- Variable(2, complex = TRUE)
  prob <- Problem(Minimize(sum_entries(Re(z))), list(Im(z) == 1))
  chain <- construct_solving_chain(prob)
  red_names <- vapply(chain@reductions, function(r) sub("^.*::", "", class(r)[[1L]]),
                       character(1L))
  ## Complex2Real should appear before FlipObjective and Dcp2Cone
  expect_true("Complex2Real" %in% red_names)
  c2r_pos <- which(red_names == "Complex2Real")
  dcp_pos <- which(red_names == "Dcp2Cone")
  expect_true(c2r_pos < dcp_pos)
})

## @cvxpy NONE
test_that("Real problems do not get Complex2Real", {
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 0))
  chain <- construct_solving_chain(prob)
  red_names <- vapply(chain@reductions, function(r) sub("^.*::", "", class(r)[[1L]]),
                       character(1L))
  expect_false("Complex2Real" %in% red_names)
})

# ── Promotion: complex scalar * matrix ────────────────────────────

## @cvxpy NONE
test_that("Complex scalar times matrix", {
  skip_if_not_installed("clarabel")
  v <- Variable(1, complex = TRUE)
  ones <- Constant(matrix(1, 2, 2))
  prob <- Problem(Maximize(sum_entries(Re(v * ones))),
                  list(Abs(v) <= 1))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 4, tolerance = 1e-3)
})

# ── Complex conjugate arithmetic ──────────────────────────────────

## @cvxpy NONE
test_that("Conj and transpose work together (expr_H)", {
  z <- Variable(c(2, 2), complex = TRUE)
  ## expr_H should be Conj(t(z))
  zh <- expr_H(z)
  expect_true(S7_inherits(zh, Conj_))
})

# ── Validate complex variable with multiple solvers ───────────────

## @cvxpy NONE
test_that("Complex solve works with SCS", {
  skip_if_not_installed("scs")
  z <- Variable(2, complex = TRUE)
  prob <- Problem(Minimize(sum_entries(Re(z))),
                  list(Im(z) == 1, Re(z) >= -3))
  result <- psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  expect_equal(result, -6, tolerance = 0.1)
})

## @cvxpy NONE
test_that("Complex solve works with Clarabel", {
  skip_if_not_installed("clarabel")
  z <- Variable(2, complex = TRUE)
  prob <- Problem(Minimize(sum_entries(Re(z))),
                  list(Im(z) == 1, Re(z) >= -3))
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  expect_equal(result, -6, tolerance = 1e-4)
})
