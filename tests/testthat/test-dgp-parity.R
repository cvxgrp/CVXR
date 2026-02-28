## CVXPY parity tests for DGP (Phase DGP-5)
## Ported from cvxpy/tests/test_dgp2dcp.py
## CVXPY reference values verified via `uv run python`

## ── Basic GP ──────────────────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_basic_gp
test_that("CVXPY parity: basic GP (box optimization)", {
  ## test_basic_gp
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  obj <- Minimize(1 / (x * y * z))
  constr <- list(2*x*y + 2*x*z + 2*y*z <= 1, x >= 2*y)
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 15.59, tolerance = 0.02)
})

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_basic_equality_constraint
test_that("CVXPY parity: equality constraint", {
  ## test_basic_equality_constraint
  x <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x), list(x == 1))
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 1.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-3)
})

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_maximum
test_that("CVXPY parity: maximum atom in GP", {
  ## test_maximum
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(max_elemwise(x * sqrt(y), 3 * x * sqrt(y))),
                  list(x == 1, y == 4))
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 6.0, tolerance = 1e-3)
})

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_minimum
test_that("CVXPY parity: minimum atom in GP", {
  ## test_minimum
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Maximize(min_elemwise(x * sqrt(y), 3 * x * sqrt(y),
                                         1 / (x * sqrt(y) + 3 * x * sqrt(y)))),
                  list(x == 1, y == 4))
  val <- psolve(prob, gp = TRUE)
  ## 1 / (2 + 6) = 0.125
  expect_equal(val, 0.125, tolerance = 1e-3)
})

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_div
test_that("CVXPY parity: division constraint", {
  ## test_div
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x * y), list(y / 3 <= x, y >= 1))
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 1/3, tolerance = 1e-3)
  expect_equal(as.numeric(value(y)), 1.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), 1/3, tolerance = 1e-2)
})

## ── Sum / trace ──────────────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_sum_scalar
test_that("CVXPY parity: sum scalar", {
  ## test_sum_scalar
  w <- Variable(pos = TRUE)
  h <- Variable(pos = TRUE)
  prob <- Problem(Minimize(h), list(w * h >= 10, sum(w) <= 5))
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 2.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(w)), 5.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(h)), 2.0, tolerance = 1e-2)
})

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_sum_vector
test_that("CVXPY parity: sum vector", {
  ## test_sum_vector
  w <- Variable(2, pos = TRUE)
  h <- Variable(2, pos = TRUE)
  prob <- Problem(Minimize(sum(h)),
                  list(multiply(w, h) >= 10, sum(w) <= 10))
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 4.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(h)), c(2, 2), tolerance = 1e-2)
  expect_equal(as.numeric(value(w)), c(5, 5), tolerance = 1e-2)
})

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_sum_matrix
test_that("CVXPY parity: sum matrix", {
  ## test_sum_matrix
  w <- Variable(c(2, 2), pos = TRUE)
  h <- Variable(c(2, 2), pos = TRUE)
  prob <- Problem(Minimize(sum(h)),
                  list(multiply(w, h) >= 10, sum(w) <= 20))
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 8.0, tolerance = 1e-2)
})

## ── one_minus_pos ─────────────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_one_minus_pos
test_that("CVXPY parity: one_minus_pos", {
  ## test_one_minus_pos
  x <- Variable(pos = TRUE)
  prob <- Problem(Maximize(x), list(one_minus_pos(x) >= 0.4))
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 0.6, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), 0.6, tolerance = 1e-2)
})

## ── eye_minus_inv ─────────────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_paper_example_eye_minus_inv
test_that("CVXPY parity: eye_minus_inv paper example", {
  ## test_paper_example_eye_minus_inv
  X <- Variable(c(2, 2), pos = TRUE)
  prob <- Problem(Minimize(matrix_trace(eye_minus_inv(X))),
                  list(geo_mean(DiagMat(X)) == 0.1,
                       geo_mean(hstack(X[1,2], X[2,1])) == 0.1))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  expect_equal(val, 2.25, tolerance = 0.05)
  expect_equal(as.numeric(value(X)), rep(0.1, 4), tolerance = 0.01)
})

## ── pf_eigenvalue ─────────────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_pf_matrix_completion
test_that("CVXPY parity: PF matrix completion", {
  ## test_pf_matrix_completion
  X <- Variable(c(3, 3), pos = TRUE)
  prob <- Problem(Minimize(pf_eigenvalue(X)),
                  list(X[1,2] == 1, X[1,3] == 1.2,
                       X[2,1] == 0.7, X[3,1] == 4, X[3,2] == 3.2,
                       X[1,1] * X[2,2] * X[3,3] == 1))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  expect_equal(val, 3.478, tolerance = 0.05)
})

## ── gmatmul ──────────────────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_gmatmul
test_that("CVXPY parity: gmatmul solve", {
  ## test_gmatmul
  x <- Variable(2, pos = TRUE)
  A <- matrix(c(-5, 1, 2, -3), 2, 2)  ## column-major: [[-5,2],[1,-3]]
  b <- c(3, 2)
  prob <- Problem(Minimize(1), list(gmatmul(A, x) == b))
  val <- psolve(prob, gp = TRUE)
  expected_x <- exp(solve(A, log(b)))
  expect_equal(as.numeric(value(x)), expected_x, tolerance = 1e-2)
})

## ── xexp ─────────────────────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_xexp
test_that("CVXPY parity: xexp product", {
  ## test_xexp
  x <- Variable(2, pos = TRUE)
  b <- c(1, 0.5)
  prob <- Problem(Minimize(prod(xexp(x))), list(x >= b))
  val <- psolve(prob, gp = TRUE)
  ## xexp([1, 0.5]) = [e, 0.5*sqrt(e)], product = e * 0.5*sqrt(e) = 0.5*e^1.5
  expected <- 0.5 * exp(1)^1.5
  expect_equal(val, expected, tolerance = 1e-3)
})

## ── pnorm on variable ───────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_pnorm
test_that("CVXPY parity: pnorm on variable", {
  ## test_pnorm (part D)
  x <- Variable(2, pos = TRUE)
  prob <- Problem(Minimize(p_norm(x, 2)), list(x == c(3, 4)))
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 5.0, tolerance = 1e-3)
})

## ── sum_squares ──────────────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_sum_squares_vector
test_that("CVXPY parity: sum_squares vector", {
  ## test_sum_squares_vector
  w <- Variable(2, pos = TRUE)
  h <- Variable(2, pos = TRUE)
  prob <- Problem(Minimize(sum_squares(h)),
                  list(multiply(w, h) >= 10, sum(w) <= 10))
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 8.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(h)), c(2, 2), tolerance = 1e-2)
  expect_equal(as.numeric(value(w)), c(5, 5), tolerance = 1e-2)
})

## ── parameter DGP ────────────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_parameter
test_that("CVXPY parity: parameter in GP", {
  ## test_parameter
  param <- Parameter(pos = TRUE)
  x <- Variable(pos = TRUE)

  value(param) <- 1.0
  prob <- Problem(Minimize(x), list(x == param))
  val1 <- psolve(prob, gp = TRUE)
  expect_equal(val1, 1.0, tolerance = 1e-3)

  value(param) <- 2.0
  val2 <- psolve(prob, gp = TRUE)
  expect_equal(val2, 2.0, tolerance = 1e-3)
})

## ── documentation example ────────────────────────────────────────

## @cvxpy test_dgp2dcp.py::TestDgp2Dcp::test_documentation_prob
test_that("CVXPY parity: documentation problem", {
  ## test_documentation_prob
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  prob <- Problem(Maximize(x * y * z),
                  list(4*x*y*z + 2*x*z <= 10,
                       x <= 2*y, y <= 2*x, z >= 1))
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 2.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(z)), 1.0, tolerance = 1e-2)
})
