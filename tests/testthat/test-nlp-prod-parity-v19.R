## CVXPY 1.9 DNLP `prod` parity (nlp_tests/test_prod.py).
##
## `prod` is a smooth atom (atoms/prod.py:73-75), so prod-containing problems
## are DNLP and solvable through the NLP path. CVXPY runs the IPOPT solves
## through a Python DerivativeChecker; CVXR has no such gradient checker, so
## these exercise the same models through psolve(..., nlp = TRUE) on the
## available R NLP backend (IPOPT or UNO via sparsediff) and verify the
## optimizer and objective values. `cp.prod(x, axis=...)` maps to CVXR's
## `prod_entries(x, axis=...)` (CVXPY axis 0 -> CVXR axis 2, axis 1 -> axis 1).

.prod_nlp_solve <- function(prob, ...) {
  skip_if_not_installed("sparsediff")
  if (!requireNamespace("ipopt", quietly = TRUE) &&
      !requireNamespace("Uno", quietly = TRUE)) {
    skip("No R NLP backend (ipopt or Uno) installed")
  }
  output <- capture.output(val <- psolve(prob, nlp = TRUE, ...), type = "output")
  expect_equal(status(prob), "optimal")
  invisible(val)
}

## ---- DNLP detection (no solver) ------------------------------------------

## @cvxpy nlp_tests/test_prod.py::TestProdDNLP::test_prod_is_smooth nlp_tests/test_prod.py::TestProdDNLP::test_prod_is_smooth_convex nlp_tests/test_prod.py::TestProdDNLP::test_prod_is_smooth_concave
test_that("DNLP: prod over positive entries is smooth, linearizable convex and concave", {
  x <- Variable(3, pos = TRUE)
  p <- prod(x)
  expect_true(is_atom_smooth(p))
  expect_true(is_smooth(p))
  expect_true(is_linearizable_convex(p))
  expect_true(is_linearizable_concave(p))
})

## @cvxpy nlp_tests/test_prod.py::TestProdDNLP::test_prod_composition_smooth
test_that("DNLP: log(prod(x)) is smooth", {
  x <- Variable(3, pos = TRUE)
  expr <- log(prod(x))
  expect_true(is_smooth(expr))
})

## @cvxpy nlp_tests/test_prod.py::TestProdDNLP::test_prod_problem_is_dnlp
test_that("DNLP: a problem with prod is DNLP", {
  x <- Variable(3, pos = TRUE)
  prob <- Problem(Maximize(prod(x)), list(sum_entries(x) <= 3))
  expect_true(is_dnlp(prob))
})

## ---- IPOPT/UNO solves -----------------------------------------------------

## @cvxpy nlp_tests/test_prod.py::TestProdIPOPT::test_prod_maximize_positive
test_that("DNLP: maximize prod over positive simplex (AM-GM)", {
  x <- Variable(3, pos = TRUE)
  value(x) <- rep(0.9, 3)
  prob <- Problem(Maximize(prod(x)), list(sum_entries(x) <= 3))
  val <- .prod_nlp_solve(prob)
  expect_equal(as.numeric(value(x)), c(1, 1, 1), tolerance = 1e-4)
  expect_equal(val, 1, tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_prod.py::TestProdIPOPT::test_prod_minimize_squared
test_that("DNLP: minimize prod(x)^2 with mixed-sign variables", {
  x <- Variable(3)
  value(x) <- c(1, -1, 2)
  prob <- Problem(Minimize(prod(x)^2),
                  list(x[1] >= 0.5, x[2] <= -0.5, x[3] >= 1,
                       sum_entries(x) == 2))
  val <- .prod_nlp_solve(prob)
  xv <- as.numeric(value(x))
  expect_true(xv[1] >= 0.5 - 1e-4)
  expect_true(xv[2] <= -0.5 + 1e-4)
  expect_true(xv[3] >= 1 - 1e-4)
  expect_equal(sum(xv), 2, tolerance = 1e-4)
  expect_equal(val, prod(xv)^2, tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_prod.py::TestProdIPOPT::test_prod_with_axis nlp_tests/test_prod.py::TestProdIPOPT::test_prod_with_axis_large_matrix
test_that("DNLP: row-wise prod (axis = 1) maximized per row (AM-GM)", {
  X <- Variable(c(2, 3), pos = TRUE)
  value(X) <- matrix(0.9, 2, 3)
  prob <- Problem(Maximize(sum_entries(prod_entries(X, axis = 1))),
                  list(sum_entries(X, axis = 1) <= 3))
  val <- .prod_nlp_solve(prob)
  expect_equal(value(X), matrix(1, 2, 3), tolerance = 1e-4)
  expect_equal(val, 2, tolerance = 1e-4)

  X2 <- Variable(c(4, 5), pos = TRUE)
  value(X2) <- matrix(0.9, 4, 5)
  prob2 <- Problem(Maximize(sum_entries(prod_entries(X2, axis = 1))),
                   list(sum_entries(X2, axis = 1) <= 5))
  val2 <- .prod_nlp_solve(prob2)
  expect_equal(value(X2), matrix(1, 4, 5), tolerance = 1e-4)
  expect_equal(val2, 4, tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_prod.py::TestProdIPOPT::test_prod_with_axis_zero nlp_tests/test_prod.py::TestProdIPOPT::test_prod_with_axis_zero_large_matrix
test_that("DNLP: column-wise prod (axis = 2) maximized per column (AM-GM)", {
  X <- Variable(c(3, 2), pos = TRUE)
  value(X) <- matrix(0.9, 3, 2)
  prob <- Problem(Maximize(sum_entries(prod_entries(X, axis = 2))),
                  list(sum_entries(X, axis = 2) <= 3))
  val <- .prod_nlp_solve(prob)
  expect_equal(value(X), matrix(1, 3, 2), tolerance = 1e-4)
  expect_equal(val, 2, tolerance = 1e-4)

  X2 <- Variable(c(5, 4), pos = TRUE)
  value(X2) <- matrix(0.9, 5, 4)
  prob2 <- Problem(Maximize(sum_entries(prod_entries(X2, axis = 2))),
                   list(sum_entries(X2, axis = 2) <= 5))
  val2 <- .prod_nlp_solve(prob2)
  expect_equal(value(X2), matrix(1, 5, 4), tolerance = 1e-4)
  expect_equal(val2, 4, tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_prod.py::TestProdIPOPT::test_prod_with_zero_start
test_that("DNLP: minimize (prod(x) - 1)^2 from a near-zero start", {
  x <- Variable(3)
  value(x) <- c(1, 0.1, 2)
  prob <- Problem(Minimize((prod(x) - 1)^2),
                  list(x[1] >= 0.5, x[3] >= 1, sum_entries(x) == 3))
  .prod_nlp_solve(prob)
  expect_equal(prod(as.numeric(value(x))), 1, tolerance = 1e-3)
})

## @cvxpy nlp_tests/test_prod.py::TestProdIPOPT::test_prod_negative_values
test_that("DNLP: minimize (prod(x) + 4)^2 with fixed mixed-sign variables", {
  x <- Variable(2)
  value(x) <- c(2, -2)
  prob <- Problem(Minimize((prod(x) + 4)^2), list(x[1] == 2, x[2] == -2))
  val <- .prod_nlp_solve(prob)
  xv <- as.numeric(value(x))
  expect_equal(xv[1], 2, tolerance = 1e-3)
  expect_equal(xv[2], -2, tolerance = 1e-3)
  expect_equal(prod(xv), -4, tolerance = 1e-3)
  expect_equal(val, 0, tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_prod.py::TestProdIPOPT::test_prod_in_constraint
test_that("DNLP: minimize sum(x) subject to prod(x) >= 8 (AM-GM)", {
  x <- Variable(3, pos = TRUE)
  value(x) <- rep(2.1, 3)
  prob <- Problem(Minimize(sum_entries(x)), list(prod(x) >= 8))
  val <- .prod_nlp_solve(prob)
  expect_equal(as.numeric(value(x)), c(2, 2, 2), tolerance = 1e-3)
  expect_equal(val, 6, tolerance = 1e-3)
})

## @cvxpy nlp_tests/test_prod.py::TestProdIPOPT::test_prod_log_composition
test_that("DNLP: log(prod(x)) equals sum(log(x))", {
  n <- 4
  x <- Variable(n, pos = TRUE)

  value(x) <- rep(0.9, n)
  prob1 <- Problem(Maximize(log(prod(x))), list(sum_entries(x) <= n))
  val1 <- .prod_nlp_solve(prob1)
  x1 <- as.numeric(value(x))

  value(x) <- rep(0.9, n)
  prob2 <- Problem(Maximize(sum_entries(log(x))), list(sum_entries(x) <= n))
  val2 <- .prod_nlp_solve(prob2)
  x2 <- as.numeric(value(x))

  expect_equal(val1, val2, tolerance = 1e-4)
  expect_equal(x1, x2, tolerance = 1e-3)
})

## @cvxpy nlp_tests/test_prod.py::TestProdIPOPT::test_prod_exp_log_relationship
test_that("DNLP: prod(x) equals exp(sum(log(x)))", {
  n <- 3
  x <- Variable(n, pos = TRUE)

  value(x) <- rep(0.9, n)
  prob1 <- Problem(Maximize(prod(x)), list(sum_entries(x) <= n))
  prod_val <- .prod_nlp_solve(prob1)
  x1 <- as.numeric(value(x))

  value(x) <- rep(0.9, n)
  prob2 <- Problem(Maximize(exp(sum_entries(log(x)))), list(sum_entries(x) <= n))
  exp_sum_log_val <- .prod_nlp_solve(prob2)
  x2 <- as.numeric(value(x))

  expect_equal(prod_val, exp_sum_log_val, tolerance = 1e-4)
  expect_equal(x1, x2, tolerance = 1e-3)
})

## @cvxpy nlp_tests/test_prod.py::TestProdIPOPT::test_prod_single_element
test_that("DNLP: maximize prod of a single element", {
  x <- Variable(1, pos = TRUE)
  value(x) <- 1
  prob <- Problem(Maximize(prod(x)), list(x <= 5))
  val <- .prod_nlp_solve(prob)
  expect_equal(as.numeric(value(x)), 5, tolerance = 1e-4)
  expect_equal(val, 5, tolerance = 1e-4)
})
