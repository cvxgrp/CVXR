## Wave 6: Remaining Misc Parity Tests
## Tests from CVXPY test_nonlinear_atoms.py, test_quad_form.py,
## test_kron_canon.py, test_convolution.py gaps
## All expected values verified via `uv run python` against CVXPY 1.8.1

# ═══════════════════════════════════════════════════════════════════════
# Entr edge cases (CVXPY test_nonlinear_atoms.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_entr
test_that("entr(0) evaluates to 0", {
  ## CVXPY: entr(0).value = -0.0 (which is 0)
  ## entr(x) = -x * log(x), limit as x→0+ is 0
  e <- Entr(Constant(0))
  expect_equal(as.numeric(value(e)), 0.0, tolerance = 1e-10)
})

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_entr
test_that("entr(-1) evaluates to -Inf (outside domain)", {
  ## CVXPY: entr(-1).value = -inf
  e <- Entr(Constant(-1))
  val <- as.numeric(value(e))
  expect_true(is.infinite(val) && val < 0)
})

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_entr
test_that("entr(1) evaluates to 0", {
  ## CVXPY: entr(1).value = -0.0 (which is -1*log(1) = 0)
  e <- Entr(Constant(1))
  expect_equal(as.numeric(value(e)), 0.0, tolerance = 1e-10)
})

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_entr
test_that("entr(0.5) evaluates to 0.3466", {
  ## CVXPY: entr(0.5).value = 0.34657359027997264
  ## -0.5 * log(0.5) = 0.5 * log(2) = 0.3466
  e <- Entr(Constant(0.5))
  expect_equal(as.numeric(value(e)), 0.34657359027997264, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════════
# Quad form: singular P (CVXPY test_quad_form.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_quad_form.py::TestNonOptimal::test_singular_quad_form
test_that("quad_form: singular (rank-deficient) PSD matrix", {
  ## CVXPY: P=diag(1,0,1), minimize x^T P x s.t. x >= 1 → value=2.0
  ## x[2] is free since P[2,2]=0 → only x[1]^2 + x[3]^2 minimized
  P <- diag(c(1, 0, 1))
  x <- Variable(3)
  prob <- Problem(Minimize(quad_form(x, P)), list(x >= 1))
  result <- psolve(prob)
  expect_equal(result, 2.0, tolerance = 1e-3)
})

## @cvxpy test_quad_form.py::TestNonOptimal::test_psd_exactly_tolerance
test_that("quad_form: NSD eigenvalue within tolerance is accepted", {
  ## CVXPY: P=diag(-1e-14, 1, 1) → is_convex=TRUE, solves to 2.0
  ## Near-zero negative eigenvalue treated as PSD within tolerance
  P <- diag(c(-1e-14, 1, 1))
  x <- Variable(3)
  qf <- quad_form(x, P)
  expect_true(is_convex(qf))
  prob <- Problem(Minimize(qf), list(x >= 1))
  result <- psolve(prob)
  expect_equal(result, 2.0, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════════
# Quad form: Parameter P (CVXPY test_quad_form.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_quad_form.py::TestNonOptimal::test_param_quad_form
test_that("quad_form: parameter P is DCP", {
  ## CVXPY: quad_form(x, P_param) with PSD=True → is_dcp=TRUE
  P <- Parameter(c(2, 2), PSD = TRUE)
  value(P) <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  x <- Variable(2)
  qf <- quad_form(x, P)
  expect_true(is_dcp(qf))
})

## @cvxpy test_quad_form.py::TestNonOptimal::test_param_quad_form
test_that("quad_form: parameter P solve", {
  ## CVXPY: x=[1,1], P=[[2,0.5],[0.5,1]] → x^T P x = 4.0
  P <- Parameter(c(2, 2), PSD = TRUE)
  value(P) <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  x <- Variable(2)
  prob <- Problem(Minimize(quad_form(x, P)), list(x >= 1))
  result <- psolve(prob)
  expect_equal(result, 4.0, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════════
# Convolution: scalar (CVXPY test_convolution.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_convolution.py::TestConvolution::test_0D_conv
test_that("conv: scalar convolution", {
  ## CVXPY: conv(constant=3, x), minimize s.t. x>=2 → value=6.0
  x <- Variable(1)
  c_val <- c(3.0)
  prob <- Problem(Minimize(conv(Constant(c_val), x)), list(x >= 2))
  result <- psolve(prob)
  expect_equal(result, 6.0, tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════════
# Kron with Parameter (CVXPY test_kron_canon.py)
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_kron_canon.py::TestKronLeftVar::test_gen_kronl_param
test_that("kron: parameter times variable", {
  ## CVXPY: kron(I_2, x) with x all-ones → 4x4 block-diagonal → sum=8
  P <- Parameter(c(2, 2))
  value(P) <- diag(2)
  x <- Variable(c(2, 2))
  expr <- kron(P, x)
  expect_equal(expr@shape, c(4L, 4L))
  expect_true(is_affine(expr))
  prob <- Problem(Minimize(sum_entries(expr)), list(x >= 1))
  result <- psolve(prob)
  ## kron(I_2, ones(2,2)) = block diagonal with two [[1,1],[1,1]]
  ## sum = 4 + 4 = 8
  expect_equal(result, 8.0, tolerance = 1e-3)
})

## @cvxpy test_kron_canon.py::TestKronRightVar::test_gen_kronr_const
test_that("kron: constant times variable", {
  ## Basic kron(constant_matrix, variable) is affine
  A <- matrix(c(1, 0, 0, 1), 2, 2)
  x <- Variable(c(2, 2))
  expr <- kron(Constant(A), x)
  expect_true(is_affine(expr))
  expect_equal(expr@shape, c(4L, 4L))
})

# ═══════════════════════════════════════════════════════════════════════
# Entr solve tests
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_entr_prob
test_that("entr: maximize sum(entr(x)) s.t. sum(x) == 1", {
  ## Maximum entropy distribution → uniform: x_i = 1/n, value = log(n)
  n <- 4L
  x <- Variable(n, nonneg = TRUE)
  prob <- Problem(Maximize(sum_entries(Entr(x))),
                  list(sum_entries(x) == 1))
  result <- psolve(prob)
  expect_equal(result, log(n), tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════════
# KL divergence and relative entropy edge cases
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_kl_div
test_that("kl_div: constant evaluation", {
  ## kl_div(x, y) = x * log(x/y) - x + y
  ## kl_div(1, 1) = 0 - 0 = 0
  e <- KlDiv(Constant(1), Constant(1))
  expect_equal(as.numeric(value(e)), 0.0, tolerance = 1e-10)
})

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_rel_entr
test_that("rel_entr: constant evaluation", {
  ## rel_entr(x, y) = x * log(x/y)
  ## rel_entr(1, 1) = 0
  e <- RelEntr(Constant(1), Constant(1))
  expect_equal(as.numeric(value(e)), 0.0, tolerance = 1e-10)
})

## @cvxpy test_nonlinear_atoms.py::TestNonlinearAtoms::test_rel_entr
test_that("rel_entr: basic solve", {
  ## minimize sum(rel_entr(x, y)) s.t. x == y, sum(x) == 1 → value = 0
  x <- Variable(3, nonneg = TRUE)
  y <- Variable(3, nonneg = TRUE)
  prob <- Problem(Minimize(sum_entries(RelEntr(x, y))),
                  list(x == y, sum_entries(x) == 1))
  result <- psolve(prob)
  expect_equal(result, 0.0, tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════════
# Logistic and log_sum_exp edge cases
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("logistic: constant evaluation", {
  ## logistic(x) = log(1 + exp(x))
  ## logistic(0) = log(2)
  e <- Logistic(Constant(0))
  expect_equal(as.numeric(value(e)), log(2), tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_log_sum_exp
test_that("log_sum_exp: constant evaluation", {
  ## log_sum_exp([0, 0, 0]) = log(3)
  e <- LogSumExp(Constant(c(0, 0, 0)))
  expect_equal(as.numeric(value(e)), log(3), tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_log_sum_exp
test_that("log_sum_exp: basic solve", {
  ## minimize log_sum_exp(x) s.t. x == [1, 2, 3]
  x <- Variable(3)
  prob <- Problem(Minimize(LogSumExp(x)), list(x == c(1, 2, 3)))
  result <- psolve(prob)
  expected <- log(exp(1) + exp(2) + exp(3))
  expect_equal(result, expected, tolerance = 1e-4)
})
