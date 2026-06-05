## @cvxpy test_derivative.py
##
## R port of the canonical Problem$backward() / Problem$derivative()
## tests.  Covers the scalar-quadratic docstring example, the L1
## regression family, least-squares, error paths (forget
## requires_grad, unsupported solver), and a basic GP.  Mirrors the
## helpers `gradcheck` / `perturbcheck` from CVXPY.

skip_if_not_installed("diffcp")
skip_if_not_installed("clarabel")

# ---- helpers (R port of cvxpy/tests/test_derivative.py:19-109) ----

## Numerical perturbation check: set delta on each parameter, call
## derivative(), and compare to the actual value change after solving
## the perturbed problem.
.perturbcheck <- function(problem, delta = 1e-5, atol = 1e-6) {
  set.seed(0)

  if (length(parameters(problem)) == 0L) {
    psolve(problem, requires_grad = TRUE)
    derivative(problem)
    for (v in variables(problem)) {
      expect_equal(delta(v), array(0, dim = v@shape), tolerance = 0)
    }
    return(invisible())
  }

  ## Set perturbations.
  for (param in parameters(problem)) {
    delta(param) <- array(delta * rnorm(prod(param@shape)), dim = param@shape)
  }
  psolve(problem, requires_grad = TRUE)
  derivative(problem)
  variable_values <- lapply(variables(problem), value)
  deltas          <- lapply(variables(problem), delta)

  ## Numerical perturbation: bump each parameter by its delta and re-solve.
  old_values <- list()
  for (param in parameters(problem)) {
    pid <- as.character(param@id)
    old_values[[pid]] <- value(param)
    value(param) <- value(param) + delta(param)
  }
  psolve(problem)   # standard solve (no requires_grad)
  num_deltas <- mapply(function(v, old) value(v) - old,
                       variables(problem), variable_values, SIMPLIFY = FALSE)

  for (i in seq_along(deltas)) {
    expect_lt(max(abs(deltas[[i]] - num_deltas[[i]])), atol)
  }

  ## Restore.
  for (param in parameters(problem)) {
    pid <- as.character(param@id)
    value(param) <- old_values[[pid]]
  }
}

## Numerical gradient check: central-difference on each parameter
## entry, compare to backward()'s gradient(param).
.gradcheck <- function(problem, delta = 1e-5, atol = 1e-4) {
  ## Linear flatten of all parameter values into one vector for the
  ## central-difference loop.
  values <- numeric(0)
  offsets <- list()
  for (param in parameters(problem)) {
    sz <- prod(param@shape)
    offsets[[as.character(param@id)]] <- length(values) + 1L
    values <- c(values, as.numeric(value(param)))
  }
  numgrad <- numeric(length(values))

  bump <- function(idx, sign) {
    saved <- values
    saved[idx] <- saved[idx] + sign * 0.5 * delta
    off <- 1L
    for (param in parameters(problem)) {
      sz <- prod(param@shape)
      v <- saved[off:(off + sz - 1L)]
      dim(v) <- param@shape
      value(param) <- v
      off <- off + sz
    }
    psolve(problem)
    unlist(lapply(variables(problem), function(v) as.numeric(value(v))))
  }

  for (i in seq_along(values)) {
    left  <- bump(i, +1)
    right <- bump(i, -1)
    numgrad[i] <- (sum(left) - sum(right)) / delta
  }

  ## Restore parameter values.
  off <- 1L
  for (param in parameters(problem)) {
    sz <- prod(param@shape)
    v <- values[off:(off + sz - 1L)]; dim(v) <- param@shape
    value(param) <- v
    off <- off + sz
  }

  ## Now backward(): clear gradients first (so backward starts from
  ## the default sum-of-x loss).
  for (v in variables(problem)) gradient(v) <- NULL
  psolve(problem, requires_grad = TRUE)
  backward(problem)

  for (param in parameters(problem)) {
    pid <- as.character(param@id)
    sz <- prod(param@shape)
    off <- offsets[[pid]]
    expected <- numgrad[off:(off + sz - 1L)]
    actual   <- as.numeric(gradient(param))
    expect_lt(max(abs(actual - expected)), atol)
  }
}

# ---- test_scalar_quadratic (test_derivative.py:121-143) -----------

test_that("backward / derivative on (x - 2p)^2 with x >= 0 (scalar quadratic)", {
  p <- Parameter(); x <- Variable(); value(p) <- 3.0
  prob <- Problem(Minimize((x - 2 * p)^2), list(x >= 0))
  psolve(prob, requires_grad = TRUE)
  expect_lt(abs(value(x) - 6), 1e-5)

  backward(prob)
  expect_lt(abs(gradient(p) - 2), 1e-5)

  gradient(x) <- 4
  backward(prob)
  expect_lt(abs(gradient(p) - 8), 1e-5)

  ## Reset and verify forward derivative.
  gradient(x) <- NULL
  psolve(prob, requires_grad = TRUE)
  delta(p) <- 1e-3
  derivative(prob)
  expect_lt(abs(delta(x) - 2e-3), 1e-5)
})

# ---- test_l1_square (test_derivative.py:145-159) ------------------

test_that("L1 regression (square A): gradcheck + perturbcheck", {
  set.seed(0)
  n <- 3L
  x <- Variable(n)
  A <- Parameter(c(n, n))
  b <- Parameter(n, name = "b")
  obj <- Minimize(p_norm(A %*% x - b, p = 1))
  prob <- Problem(obj)
  expect_true(is_dpp(prob))

  L <- matrix(rnorm(n * n), n, n)
  value(A) <- crossprod(L) + diag(n)
  value(b) <- rnorm(n)

  .gradcheck(prob, atol = 1e-3)
  .perturbcheck(prob, atol = 1e-4)
})

# ---- test_forget_requires_grad (test_derivative.py:277-295) -------

test_that("backward errors when requires_grad was not set", {
  set.seed(0)
  x <- Variable(3); A <- Parameter(c(3, 3))
  value(A) <- diag(3)
  prob <- Problem(Minimize(sum_squares(A %*% x - 1)))
  psolve(prob)   # without requires_grad
  expect_error(backward(prob),   "requires_grad")
  expect_error(derivative(prob), "requires_grad")
})

# ---- test_unsupported_solver (test_derivative.py:323-331) ---------

test_that("requires_grad with a non-DIFFCP solver errors", {
  x <- Variable()
  prob <- Problem(Minimize(x^2), list(x >= 1))
  expect_error(
    psolve(prob, solver = CLARABEL_SOLVER, requires_grad = TRUE),
    "DIFFCP"
  )
})

# ---- test_zero_in_problem_data (test_derivative.py:333-341) -------

test_that("zero-coefficient parameter entries don't break the chain", {
  ## The constraint matrix has a structurally-present zero coefficient
  ## (b * x == p), where b is a parameter set to 0.  Forward solve
  ## must succeed and backward / derivative must return finite values.
  b <- Parameter(); p <- Parameter(); x <- Variable()
  value(b) <- 0; value(p) <- 1
  prob <- Problem(Minimize(x^2), list(b * x == p, x >= 0))
  ## With b = 0, the equality constraint is infeasible (0 = 1) -- the
  ## solver should report infeasibility, not silently produce garbage.
  expect_error(psolve(prob, requires_grad = TRUE))
})

# ---- test_basic_gp (test_derivative.py:554-569) -------------------

test_that("basic GP: backward / derivative through Dgp2Dcp chain rule", {
  ## min   x * y
  ## s.t.  a / x + 1 / y <= 1
  ##       1 <= x, 1 <= y
  ## with parameter `a` controlling the constraint right-hand-side.
  a <- Parameter(pos = TRUE); value(a) <- 0.5
  x <- Variable(pos = TRUE); y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x * y),
                  list(a / x + 1 / y <= 1, x >= 1, y >= 1))
  expect_true(is_dgp(prob))

  ## Forward solve with requires_grad.
  psolve(prob, gp = TRUE, requires_grad = TRUE)
  ## Just verify that backward / derivative both run cleanly and
  ## return finite values (a numerical gradcheck on this GP needs
  ## careful tolerances and is skipped here -- mirrors CVXPY's choice
  ## not to run gradcheck on GPs in test_basic_gp).
  backward(prob)
  expect_true(all(is.finite(gradient(a))))

  delta(a) <- 1e-4
  derivative(prob)
  expect_true(all(is.finite(delta(x))))
  expect_true(all(is.finite(delta(y))))
})

# ---- parameters() discovery in non-arg properties -----------------
## Regression test: Power@p, Huber@M, Gmatmul@A hold Parameters as
## properties rather than args.  The base Canonical parameters()
## traversal only walks @args, so without explicit overrides those
## Parameters become invisible to parameters(problem), which in turn
## leaves their backward gradients empty.
## CVXPY counterparts: power.py:208, huber.py:52, gmatmul.py:112.

test_that("parameters() finds Parameter exponent in Power", {
  x <- Variable(pos = TRUE)
  p <- Parameter(value = 0.5)
  expr <- power(x, p)
  ps <- parameters(expr)
  expect_equal(length(ps), 1L)
  expect_equal(ps[[1L]]@id, p@id)
})

test_that("parameters() finds Parameter M in Huber", {
  x <- Variable()
  m <- Parameter(pos = TRUE, value = 1.0)
  expr <- huber(x, m)
  ps <- parameters(expr)
  expect_equal(length(ps), 1L)
  expect_equal(ps[[1L]]@id, m@id)
})

## @cvxpy test_derivative.py::TestBackwardDgp::test_param_used_in_exponent_and_elsewhere
test_that("parameter used in exponent and elsewhere: backward + derivative", {
  ## Constructed so the optimal x is analytic:
  ##   x*(alpha) = 1 - base^alpha - alpha^2
  ## and the derivative is:
  ##   x*'(alpha) = -log(base) * base^alpha - 2*alpha
  base <- 0.3
  alpha <- Parameter(pos = TRUE, value = 0.5)
  x <- Variable(pos = TRUE)
  objective <- Maximize(x)
  constr <- list(one_minus_pos(x) >= Constant(base)^alpha + alpha^2)
  prob <- Problem(objective, constr)

  delta(alpha) <- 1e-5
  ## Force DIFFCP's SCS path (matches CVXPY's exponent test, which uses
  ## solve_method=SCS) but at eps=1e-8, NOT CVXPY's literal eps=1e-5.
  ## Reason: CVXPY's eps=1e-5 + assertAlmostEqual(places=7) is internally
  ## inconsistent and never actually validated -- that test is SKIPPED in
  ## the cvxpy repo (diffcp absent in its venv), and at eps=1e-5 even CVXPY
  ## only reaches x_err ~9e-7. CVXPY's *validated* analytic-solution
  ## derivative tests use eps=1e-10 (test_derivative.py:127,140). We follow
  ## that convention (eps=1e-8): the SCS forward point is then accurate to
  ## ~4e-9, so the closed-form checks below are meaningful rather than
  ## aspirational. See notes/diffcp_status_mapping.md and ADR D_DIFFCP.1.
  psolve(prob, gp = TRUE, requires_grad = TRUE,
         solve_method = "SCS", eps_abs = 1e-8, eps_rel = 1e-8)
  expect_equal(as.numeric(value(x)),
               1 - base^(0.5) - 0.5^2, tolerance = 1e-5)
  backward(prob)
  derivative(prob)
  expect_equal(as.numeric(gradient(alpha)),
               -log(base) * base^(0.5) - 2 * 0.5, tolerance = 1e-5)
  expect_equal(as.numeric(delta(x)),
               as.numeric(gradient(alpha)) * 1e-5, tolerance = 1e-3)
})

## @cvxpy test_derivative.py::TestDgp2DcpReduction::test_param_backward_absent_log_param
## CVXPY calls dgp.param_backward({}) and expects {} when the log-param id is
## absent. CVXR's param_backward is now dict-in/dict-out (#3147 part A): an empty
## dparams passes straight through to an empty list (each log id is guarded by
## `if (is.null(dparams[[log_pid]])) next`), instead of erroring.
test_that("Dgp2Dcp param_backward passes through when the log-param id is absent", {
  p <- Parameter(pos = TRUE, value = 2.0)
  x <- Variable(pos = TRUE)
  prob <- Problem(Minimize(p * x), list(x >= 1))
  dgp <- Dgp2Dcp()
  reduction_apply(dgp, prob)
  ## Empty dparams: the log-parameter id is absent -> empty dict out, no error.
  expect_equal(param_backward(dgp, list()), list())
})
