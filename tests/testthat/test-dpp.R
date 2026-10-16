## Tests for Disciplined Parameterized Programming (DPP)
## Verifies: is_dpp detection, tensor fast path, parameter re-solve

# ── is_dpp detection ──────────────────────────────────────────────

## @cvxpy NONE
test_that("DPP: simple linear parameter problem is DPP", {
  p <- Parameter(value = 1)
  x <- Variable()
  prob <- Problem(Minimize(p * x), list(x >= 1))
  expect_true(is_dpp(prob))
})

## @cvxpy NONE
test_that("DPP: quadratic-in-params is NOT DPP", {
  p1 <- Parameter(value = 1)
  p2 <- Parameter(value = 2)
  x <- Variable()
  prob <- Problem(Minimize(p1 * p2 * x), list(x >= 1))
  expect_false(is_dpp(prob))
})

## @cvxpy NONE
test_that("DPP: parameter in constraint is DPP", {
  p <- Parameter(value = 3)
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= p))
  expect_true(is_dpp(prob))
})

## @cvxpy NONE
test_that("DPP: no parameters is DPP", {
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= 1))
  expect_true(is_dpp(prob))
})

## @cvxpy NONE
test_that("DPP: kron with parameter arg is NOT DPP", {
  P <- Parameter(c(2, 2), value = diag(2))
  x <- Variable(c(2, 2))
  prob <- Problem(Minimize(sum_entries(kron(P, x))), list(x >= 1))
  expect_false(is_dpp(prob))
})

## @cvxpy NONE
test_that("DPP: conv with parameter kernel is NOT DPP", {
  k <- Parameter(3, value = c(1, 2, 1))
  x <- Variable(5)
  prob <- Problem(Minimize(sum_entries(conv(k, x))), list(x >= 0, x <= 1))
  expect_false(is_dpp(prob))
})

## dotsort.py:104-113 -- under a DPP scope dotsort is convex iff X is constant
## or W is parameter-affine. Measured against CVXPY 1.9.3:
## notes/audit/repro/13_dotsort_param_w.{py,R}.
## @cvxpy NONE
test_that("DPP: dotsort with a parameter-affine W is DPP, with exp(param) W is not", {
  x <- Variable(3)
  p <- Parameter(2, nonneg = TRUE, value = c(1, 1))
  expect_true(is_dpp(Problem(Minimize(dotsort(x, p)), list(x >= c(1, 2, 3)))))
  expect_true(is_dpp(Problem(Minimize(dotsort(x, 2 * p)), list(x >= c(1, 2, 3)))))
  q <- Parameter(2, value = c(0, 0))
  expect_false(is_dpp(Problem(Minimize(dotsort(x, exp(q))), list(x >= c(1, 2, 3)))))
  ## X constant: DPP regardless of W
  r <- Parameter(value = 2)
  expect_true(is_dpp(Problem(Minimize(dotsort(Constant(c(1, 2, 3)), r^2)))))
})

## dotsort_canon.py:30-33 -- a non-Constant W stays symbolic so the compiled
## tensor depends on the parameter and the DPP fast path is used.
## @cvxpy NONE
test_that("DPP: dotsort with a Parameter W caches the tensor and re-solves", {
  x <- Variable(3)
  p <- Parameter(2, nonneg = TRUE, value = c(1, 1))
  prob <- Problem(Minimize(dotsort(x, p)), list(x >= c(1, 2, 3)))
  expect_equal(psolve(prob, solver = "CLARABEL"), 5, tolerance = 1e-5)
  expect_false(is.null(prob@.cache$param_prog))
  value(p) <- c(1, 0.5)
  expect_equal(psolve(prob, solver = "CLARABEL"), 4, tolerance = 1e-5)
  value(p) <- c(0, 0)
  expect_equal(psolve(prob, solver = "CLARABEL"), 0, tolerance = 1e-5)
})

## The warning is emitted when the solving chain is built; CVXPY 1.9.3 caches
## the chain and warns ONCE per problem (measured: 1, 0, 0 warnings over three
## solves). CVXR's chain cache behaves the same, so only the first solve warns.
## @cvxpy NONE
test_that("DPP: dotsort with exp(param) W solves through EvalParams", {
  x <- Variable(3)
  q <- Parameter(2, value = c(0, 0))
  prob <- Problem(Minimize(dotsort(x, exp(q))), list(x >= c(1, 2, 3)))
  expect_warning(v1 <- psolve(prob, solver = "CLARABEL"), class = "DPPWarning")
  expect_equal(v1, 5, tolerance = 1e-5)
  value(q) <- c(0, log(0.5))
  expect_no_warning(v2 <- psolve(prob, solver = "CLARABEL"))
  expect_equal(v2, 4, tolerance = 1e-5)
})

# ── DPP re-solve (fast path) ─────────────────────────────────────

## @cvxpy NONE
test_that("DPP: parameter change re-solve gives correct values", {
  p <- Parameter(value = 1)
  x <- Variable()
  prob <- Problem(Minimize(p * x), list(x >= 1))

  r1 <- psolve(prob)
  expect_equal(r1, 1.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-4)

  value(p) <- 5
  r2 <- psolve(prob)
  expect_equal(r2, 5.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("DPP: parameter in constraint re-solve", {
  b <- Parameter(value = 3)
  x <- Variable()
  prob <- Problem(Minimize(x), list(x >= b))

  r1 <- psolve(prob)
  expect_equal(r1, 3.0, tolerance = 1e-4)

  value(b) <- 7
  r2 <- psolve(prob)
  expect_equal(r2, 7.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("DPP: vector parameter re-solve", {
  n <- 3L
  c_param <- Parameter(n, value = c(1, 2, 3))
  x <- Variable(n)
  prob <- Problem(Minimize(t(c_param) %*% x), list(x >= 0, x <= 1))

  r1 <- psolve(prob)
  ## Minimize c'x with 0<=x<=1: set x[i]=0 where c[i]>0, x[i]=1 where c[i]<0
  expect_equal(r1, 0.0, tolerance = 1e-4)

  value(c_param) <- c(-1, -2, -3)
  r2 <- psolve(prob)
  ## All coefficients negative: x = (1,1,1), value = -1-2-3 = -6
  expect_equal(r2, -6.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("DPP: matrix parameter in constraint", {
  A <- Parameter(c(2, 2), value = diag(2))
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(A %*% x >= 1))

  r1 <- psolve(prob)
  expect_equal(r1, 2.0, tolerance = 1e-4)

  ## Change to [[2,0],[0,2]] -> x >= 0.5 each -> sum = 1.0
  value(A) <- 2 * diag(2)
  r2 <- psolve(prob)
  expect_equal(r2, 1.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("DPP: cached param_prog exists after first solve", {
  p <- Parameter(value = 1)
  x <- Variable()
  prob <- Problem(Minimize(p * x), list(x >= 1))

  psolve(prob)
  ## After solving a DPP problem, the param_prog should be cached
  expect_true(!is.null(prob@.cache$param_prog))
})

## @cvxpy NONE
test_that("DPP: non-DPP problem does NOT cache param_prog", {
  p1 <- Parameter(value = 1)
  p2 <- Parameter(value = 2)
  x <- Variable()
  prob <- Problem(Minimize(p1 * p2 * x), list(x >= 1))

  psolve(prob)
  ## Non-DPP: EvalParams is used, no tensor caching
  expect_null(prob@.cache$param_prog)
})

# ── DPP with different solvers ────────────────────────────────────

## @cvxpy NONE
test_that("DPP: fast path works with SCS", {
  p <- Parameter(value = 2)
  x <- Variable()
  prob <- Problem(Minimize(p * x), list(x >= 1))

  r1 <- psolve(prob, solver = "SCS")
  expect_equal(r1, 2.0, tolerance = 1e-3)

  value(p) <- 10
  r2 <- psolve(prob, solver = "SCS")
  expect_equal(r2, 10.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("DPP: fast path works with Clarabel", {
  p <- Parameter(value = 2)
  x <- Variable()
  prob <- Problem(Minimize(p * x), list(x >= 1))

  r1 <- psolve(prob, solver = "CLARABEL")
  expect_equal(r1, 2.0, tolerance = 1e-4)

  value(p) <- 10
  r2 <- psolve(prob, solver = "CLARABEL")
  expect_equal(r2, 10.0, tolerance = 1e-4)
})

# ── Edge cases ────────────────────────────────────────────────────

## @cvxpy NONE
test_that("DPP: unset parameter errors with informative message", {
  p <- Parameter()
  x <- Variable()
  prob <- Problem(Minimize(p * x), list(x >= 0))
  expect_error(psolve(prob), class = "ParameterError")
})

## @cvxpy NONE
test_that("DPP: multiple parameters", {
  a <- Parameter(value = 1)
  b <- Parameter(value = 2)
  x <- Variable()
  prob <- Problem(Minimize(a * x + b), list(x >= 1))

  r1 <- psolve(prob)
  expect_equal(r1, 3.0, tolerance = 1e-4)  # 1*1 + 2

  value(a) <- 3
  value(b) <- -1
  r2 <- psolve(prob)
  expect_equal(r2, 2.0, tolerance = 1e-4)  # 3*1 + (-1)
})

## @cvxpy NONE
test_that("DPP: parameter affects optimal variable values", {
  p <- Parameter(c(2, 1), value = matrix(c(1, -1), 2, 1))
  x <- Variable(2)
  prob <- Problem(Minimize(t(p) %*% x), list(x >= 0, sum_entries(x) == 1))

  r1 <- psolve(prob)
  ## min x1 - x2 s.t. x>=0, x1+x2=1 -> x=(0,1), val=-1
  expect_equal(r1, -1.0, tolerance = 1e-4)

  value(p) <- matrix(c(-1, 1), 2, 1)
  r2 <- psolve(prob)
  ## min -x1 + x2 s.t. x>=0, x1+x2=1 -> x=(1,0), val=-1
  expect_equal(r2, -1.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(1, 0), tolerance = 1e-3)
})

# ── DPP + QP re-solve ────────────────────────────────────────────────

## @cvxpy NONE
test_that("DPP: QP with parameter does not crash on re-solve", {
  set.seed(1)
  n <- 10; m <- 5
  A <- matrix(rnorm(m * n), m, n)
  b <- rnorm(m)
  x <- Variable(n)
  gamma <- Parameter(nonneg = TRUE)
  obj <- Minimize(sum_squares(A %*% x - b) + gamma * p_norm(x, 1))
  prob <- Problem(obj)

  expect_true(is_dpp(prob))
  value(gamma) <- 1.0
  r1 <- psolve(prob)
  expect_true(is.finite(r1))

  value(gamma) <- 2.0
  r2 <- psolve(prob)        # was crashing
  expect_true(is.finite(r2))
  expect_true(r2 >= r1)     # larger gamma -> more regularization -> larger obj
})
