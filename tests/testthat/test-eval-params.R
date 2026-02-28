## Tests for EvalParams reduction / Parameter support
## CVXPY SOURCE: tests/test_problem.py (parameter-related tests)

## @cvxpy NONE
test_that("scalar parameter in objective", {
  p <- Parameter(value = 2)
  x <- Variable()
  prob <- Problem(Minimize(p * x + 1), list(x >= -5))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(drop(value(x)), -5, tolerance = 1e-4)
  expect_equal(result, -9, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("change parameter and re-solve", {
  p <- Parameter(value = 2)
  x <- Variable()

  prob1 <- Problem(Minimize(p * x), list(x >= -10))
  r1 <- psolve(prob1)
  expect_equal(r1, -20, tolerance = 1e-4)

  ## Change parameter value, re-solve
  value(p) <- 5
  prob2 <- Problem(Minimize(p * x), list(x >= -10))
  r2 <- psolve(prob2)
  expect_equal(r2, -50, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("parameter in constraint", {
  b <- Parameter(value = 3)
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= b))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(result, 3, tolerance = 1e-4)

  ## Change bound
  value(b) <- 7
  prob2 <- Problem(Minimize(x), list(x >= b))
  r2 <- psolve(prob2)
  expect_equal(r2, 7, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("vector parameter", {
  c_param <- Parameter(3, value = c(1, 2, 3))
  x <- Variable(3)
  prob <- Problem(Minimize(sum(c_param * x)), list(x >= 0, x <= 1))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  ## All x should be 0 since c_param > 0
  expect_equal(drop(value(x)), c(0, 0, 0), tolerance = 1e-4)
  expect_equal(result, 0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("parameter in both objective and constraints", {
  a <- Parameter(value = 1)
  b <- Parameter(value = 5)
  x <- Variable()
  prob <- Problem(Minimize(a * x), list(x >= b))
  result <- psolve(prob)
  expect_equal(result, 5, tolerance = 1e-4)

  ## Change both
  value(a) <- -1
  value(b) <- 3
  prob2 <- Problem(Minimize(a * x), list(x <= b))
  r2 <- psolve(prob2)
  expect_equal(r2, -3, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("unspecified parameter errors", {
  p <- Parameter()
  x <- Variable()
  prob <- Problem(Minimize(p * x), list(x >= 0))
  expect_error(psolve(prob), "unspecified parameter")
})

## @cvxpy NONE
test_that("parameter with nonneg attribute", {
  p <- Parameter(nonneg = TRUE, value = 2)
  x <- Variable()
  prob <- Problem(Minimize(p * x), list(x >= -1, x <= 1))
  result <- psolve(prob)
  expect_equal(result, -2, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("matrix parameter", {
  A_val <- matrix(c(1, 0, 0, 2), 2, 2)
  A <- Parameter(c(2, 2), value = A_val)
  x <- Variable(2)
  prob <- Problem(Minimize(sum(A %*% x)), list(x >= 0, x <= 1))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  ## sum(A %*% x) = sum([1*x1, 2*x2]) = x1 + 2*x2; minimize → x = (0,0)
  expect_equal(result, 0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("parameter-free problem is unaffected", {
  ## Ensure EvalParams is not prepended when not needed
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 1))
  result <- psolve(prob)
  expect_equal(result, 1, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("parameter in quadratic objective", {
  p <- Parameter(nonneg = TRUE, value = 2)
  x <- Variable()
  prob <- Problem(Minimize(p * x^2 + x), list(x >= -10))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  ## min 2x^2 + x → x = -1/4, obj = -1/8
  expect_equal(drop(value(x)), -0.25, tolerance = 1e-3)
  expect_equal(result, -0.125, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("parameter in SOCP constraint", {
  ## ||x|| <= p, minimize -sum(x)
  p <- Parameter(value = 1)
  x <- Variable(2)
  prob <- Problem(Minimize(-sum(x)), list(cvxr_norm(x, 2) <= p))
  result <- psolve(prob)
  expect_equal(status(prob), "optimal")
  ## x* = (1/sqrt(2), 1/sqrt(2)), obj = -sqrt(2)
  expect_equal(result, -sqrt(2), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("EvalParams reduction directly", {
  ## Test the reduction interface directly
  p <- Parameter(value = 5)
  x <- Variable()
  prob <- Problem(Minimize(p * x), list(x >= 1))

  ep <- CVXR:::EvalParams()
  expect_true(CVXR:::reduction_accepts(ep, prob))

  result <- CVXR:::reduction_apply(ep, prob)
  new_prob <- result[[1]]

  ## New problem should have no parameters
  expect_equal(length(parameters(new_prob@objective)), 0L)
  expect_equal(length(parameters(new_prob@constraints[[1]])), 0L)

  ## Original problem still has parameter

  expect_equal(length(parameters(prob@objective)), 1L)
})
