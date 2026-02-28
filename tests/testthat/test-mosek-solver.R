## Tests for MOSEK solver integration
## All tests are wrapped in skip_if_not_installed("Rmosek") so they
## pass cleanly when Rmosek is not available.

## @cvxpy NONE
test_that("MOSEK LP", {
  skip_if_not_installed("Rmosek")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(x[1] + x[2] >= 3, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "MOSEK")
  expect_equal(val, 3.0, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(3.0, 0.0), ncol = 1), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("MOSEK QP", {
  skip_if_not_installed("Rmosek")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x) + x[1]),
                  list(x[1] + x[2] == 1, x[1] >= 0, x[2] >= 0))
  val <- psolve(prob, solver = "MOSEK")
  expect_equal(val, 0.875, tolerance = 1e-4)
  expect_equal(value(x), matrix(c(0.25, 0.75), ncol = 1), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("MOSEK SOCP", {
  skip_if_not_installed("Rmosek")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(cvxr_norm(x, 2) <= 1))
  val <- psolve(prob, solver = "MOSEK")
  expect_equal(val, -sqrt(2), tolerance = 1e-4)
  expect_equal(value(x), matrix(c(-1/sqrt(2), -1/sqrt(2)), ncol = 1), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("MOSEK ExpCone via entr", {
  skip_if_not_installed("Rmosek")
  x <- Variable(2)
  prob <- Problem(Maximize(entr(x[1]) + entr(x[2])),
                  list(x[1] + x[2] == 1))
  val <- psolve(prob, solver = "MOSEK")
  expect_equal(val, log(2), tolerance = 1e-4)
  expect_equal(value(x), matrix(c(0.5, 0.5), ncol = 1), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("MOSEK PowCone via pnorm", {
  skip_if_not_installed("Rmosek")
  x <- Variable(3)
  prob <- Problem(Minimize(cvxr_norm(x, 1.5)),
                  list(sum(x) == 1, x >= 0))
  val <- psolve(prob, solver = "MOSEK")
  ## Reference: CVXPY Clarabel gives ~0.6934
  expect_equal(val, 0.6934, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), rep(1/3, 3), tolerance = 1e-2)
})

## @cvxpy NONE
test_that("MOSEK SDP via lambda_max", {
  skip_if_not_installed("Rmosek")
  X <- Variable(shape = c(2, 2), symmetric = TRUE)
  prob <- Problem(Minimize(lambda_max(X)),
                  list(X[1, 1] == 1, X[2, 2] == 2, X[1, 2] >= 0.5))
  val <- psolve(prob, solver = "MOSEK")
  expected <- 1.5 + sqrt(2) / 2  # (3 + sqrt(2)) / 2
  expect_equal(val, expected, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("MOSEK infeasible problem", {
  skip_if_not_installed("Rmosek")
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 1, x <= -1))
  val <- psolve(prob, solver = "MOSEK")
  expect_true(is.infinite(val) && val > 0)
  expect_true(status(prob) %in% c("infeasible", "infeasible_inaccurate"))
})

## @cvxpy NONE
test_that("MOSEK unbounded problem", {
  skip_if_not_installed("Rmosek")
  x <- Variable()
  prob <- Problem(Minimize(x))
  val <- psolve(prob, solver = "MOSEK")
  expect_true(is.infinite(val) && val < 0)
  expect_true(status(prob) %in% c("unbounded", "unbounded_inaccurate"))
})

## @cvxpy NONE
test_that("MOSEK LP matches Clarabel", {
  skip_if_not_installed("Rmosek")
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + 2 * x[2]),
                  list(x[1] + x[2] >= 3, x[1] >= 0, x[2] >= 0))
  val_mosek <- psolve(prob, solver = "MOSEK")
  val_clarabel <- psolve(prob, solver = "CLARABEL")
  expect_equal(val_mosek, val_clarabel, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("MOSEK QP matches Clarabel", {
  skip_if_not_installed("Rmosek")
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x) + x[1]),
                  list(x[1] + x[2] == 1, x[1] >= 0, x[2] >= 0))
  val_mosek <- psolve(prob, solver = "MOSEK")
  val_clarabel <- psolve(prob, solver = "CLARABEL")
  expect_equal(val_mosek, val_clarabel, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("MOSEK SOCP matches Clarabel", {
  skip_if_not_installed("Rmosek")
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  prob <- Problem(Minimize(x[1] + x[2]),
                  list(cvxr_norm(x, 2) <= 1))
  val_mosek <- psolve(prob, solver = "MOSEK")
  val_clarabel <- psolve(prob, solver = "CLARABEL")
  expect_equal(val_mosek, val_clarabel, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("MOSEK ExpCone matches Clarabel", {
  skip_if_not_installed("Rmosek")
  skip_if_not_installed("clarabel")
  x <- Variable(2)
  prob <- Problem(Maximize(entr(x[1]) + entr(x[2])),
                  list(x[1] + x[2] == 1))
  val_mosek <- psolve(prob, solver = "MOSEK")
  val_clarabel <- psolve(prob, solver = "CLARABEL")
  expect_equal(val_mosek, val_clarabel, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("MOSEK SDP matches Clarabel", {
  skip_if_not_installed("Rmosek")
  skip_if_not_installed("clarabel")
  X <- Variable(shape = c(2, 2), symmetric = TRUE)
  prob <- Problem(Minimize(lambda_max(X)),
                  list(X[1, 1] == 1, X[2, 2] == 2, X[1, 2] >= 0.5))
  val_mosek <- psolve(prob, solver = "MOSEK")
  val_clarabel <- psolve(prob, solver = "CLARABEL")
  expect_equal(val_mosek, val_clarabel, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("MOSEK log_det problem", {
  skip_if_not_installed("Rmosek")
  X <- Variable(shape = c(2, 2), symmetric = TRUE)
  prob <- Problem(Maximize(log_det(X)),
                  list(X[1, 1] <= 2, X[2, 2] <= 3, X == t(X)))
  val <- psolve(prob, solver = "MOSEK")
  ## log(det(diag(2,3))) = log(6)
  expect_equal(val, log(6), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("MOSEK norm_nuc problem", {
  skip_if_not_installed("Rmosek")
  X <- Variable(shape = c(2, 2))
  prob <- Problem(Minimize(norm_nuc(X)),
                  list(X[1, 1] == 1, X[2, 2] == 1, X[1, 2] == 0, X[2, 1] == 0))
  val <- psolve(prob, solver = "MOSEK")
  expect_equal(val, 2.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("MOSEK solver name constant is exported", {
  expect_equal(CVXR::MOSEK_SOLVER, "MOSEK")
})
