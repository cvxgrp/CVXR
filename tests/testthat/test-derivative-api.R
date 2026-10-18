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

## @cvxpy test_derivative.py::TestBackward::test_unsupported_solver
test_that("requires_grad with a non-DIFFCP solver errors", {
  x <- Variable()
  param <- Parameter(value = 1)
  prob <- Problem(Minimize(x), list(x <= param))
  expect_error(
    psolve(prob, solver = CLARABEL_SOLVER, requires_grad = TRUE),
    "When requires_grad is True, the only supported solver is SCS"
  )
})

## @cvxpy NONE
## The rest of the requires_grad guard (problem.py:1061-1078): QCP is refused,
## SCS is accepted (and routed to DIFFCP), and a non-DCP problem is a DPPError
## (not DPP), with no separate DCP check.
test_that("requires_grad guard mirrors CVXPY: QCP, SCS, non-DCP", {
  x <- Variable()
  p <- Parameter(value = 1.5)
  prob <- Problem(Minimize(ceiling(x)), list(x >= p))
  expect_error(psolve(prob, qcp = TRUE, requires_grad = TRUE),
               "Cannot compute gradients of DQCP problems.", fixed = TRUE)
  y <- Variable()
  q <- Parameter(value = 2)
  prob2 <- Problem(Minimize(y), list(y >= q))
  expect_equal(psolve(prob2, solver = "SCS", requires_grad = TRUE), 2, tolerance = 1e-4)
  z <- Variable()
  prob3 <- Problem(Minimize(sqrt(z)), list(z >= 1))
  expect_error(psolve(prob3, requires_grad = TRUE), class = "DPPError")
})

# ---- test_zero_in_problem_data (test_derivative.py:333-340) -------

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

# ---- test_basic_gp (test_derivative.py:668-683) -------------------

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
## CVXPY counterparts: power.py:233-249, huber.py:51-53, gmatmul.py:111-113.

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

## Regression guard, 2026-08-11: DIFFCP extends SCS_Solver and takes the same
## lower-tri sqrt(2)-scaled svec PSD format.  That used to reach DIFFCP by S7
## method inheritance (`solver_psd_format_mat` on SCS_Solver); once the format
## became a class PROPERTY (D_19.6 step A), inheritance no longer supplies it --
## `.fast_new` gives a subclass only what its own constructor names.  A blanket
## NA there stuffs PSD blocks as full n x n matrices while the cone dims still
## declare n(n+1)/2 svec rows, and the solve dies on the size mismatch.
##
## ABSOLUTE ORACLE (see test-psd-dual-recovery.R): for  min tr(CX)  s.t. X >= I,
## the optimum is tr(C) and the dual is exactly C.
## @cvxpy NONE
test_that("DIFFCP solves a PSD problem with the SCS svec format", {
  n <- 3L
  C <- diag(c(2, 3, 4))
  X <- Variable(c(n, n), symmetric = TRUE)
  prob <- Problem(Minimize(matrix_trace(C %*% X)), list((X - diag(n)) %>>% 0))

  ## The stuffed A must have n(n+1)/2 = 6 PSD rows, not n^2 = 9.
  d <- problem_data(prob, solver = "DIFFCP")
  expect_equal(nrow(d[[1L]]$A), 6L)
  expect_equal(d[[1L]]$dims@psd, n)

  prob2 <- Problem(Minimize(matrix_trace(C %*% X)), list((X - diag(n)) %>>% 0))
  expect_equal(psolve(prob2, solver = "DIFFCP"), sum(diag(C)), tolerance = 1e-5)
  expect_equal(status(prob2), "optimal")
})

## @cvxpy NONE
## CVXPY problem.py:1064-1066 refuses requires_grad on a non-DPP problem with
## DPPError before solving. CVXR used to solve it and then crash in backward()
## ("no applicable method for `@` applied to NULL").
test_that("requires_grad on a non-DPP problem is a DPPError", {
  p <- Parameter(value = 2)
  q <- Parameter(value = 3)
  x <- Variable()
  prob <- Problem(Minimize(square(x) - p * q * x))
  expect_false(is_dpp(prob))
  expect_error(psolve(prob, requires_grad = TRUE), class = "DPPError")
  ## A parameter in quad_over_lin's denominator is not DPP either.
  z <- Variable()
  r <- Parameter(value = 1)
  expect_error(psolve(Problem(Minimize(quad_over_lin(z, r) + z)), requires_grad = TRUE),
               class = "DPPError")
  ## The DGP context is checked with gp = TRUE: quad_over_lin over a positive
  ## parameter is DGP-DPP and still differentiates.
  xg <- Variable(2, pos = TRUE)
  yg <- Parameter(pos = TRUE, value = 1)
  pg <- Problem(Minimize(quad_over_lin(xg, yg)), list(xg >= 1))
  expect_equal(psolve(pg, gp = TRUE, requires_grad = TRUE), 2, tolerance = 1e-4)
})

## Derivatives through a Hermitian variable (claims-audit row C0587, repro 18).
## The compact imaginary part is the strict upper triangle in ROW-major order
## (CVXPY `np.triu_indices(n, k=1)`, complex2real.py:128-179; CVXR's canon
## `vec_to_upper_tri`). The derivative helpers used column-major `upper.tri`,
## which coincides for n <= 3 and permutes the (1,4)/(2,3) entries for n = 4.
## minimize ||X - p*C||_F^2 s.t. trace(X) == 1 has X* = p*C + (1 - p*tr(C))/n * I,
## so dX/dp = D = C - tr(C)/n * I and, for a gradient seed G on X,
## dL/dp = sum(Re(G)*Re(D)) + sum_{i<j} Im(G)_ij * Im(D)_ij.
## RECORDED DIVERGENCE (ADR D_PORT.7): CVXPY 1.9.3 seeds `np.asarray(variable.gradient,
## dtype=np.float64)` (problems/problem.py:1259-1260) and so discards the imaginary
## part of a complex gradient (5.7778 for n = 3 where 6.3678 is right); CVXR keeps it.
## @cvxpy NONE
test_that("Hermitian variable: delta and gradient use the row-major compact imaginary order", {
  skip_if_not_installed("scs")
  triu_rm <- function(M) {
    idx <- which(upper.tri(M), arr.ind = TRUE)
    idx <- idx[order(idx[, 1L], idx[, 2L]), , drop = FALSE]
    M[idx]
  }
  hermitian <- function(n, re_seed, im_vals) {
    A <- matrix(seq_len(n * n) / (n * n), n, n, byrow = TRUE) + re_seed
    U <- matrix(0, n, n)
    idx <- which(upper.tri(U), arr.ind = TRUE)
    idx <- idx[order(idx[, 1L], idx[, 2L]), , drop = FALSE]
    U[idx] <- im_vals[seq_len(n * (n - 1L) / 2L)]
    (A + t(A)) / 2 + 1i * (U - t(U))
  }
  for (n in c(3L, 4L)) {
    C <- hermitian(n, 0, (1:6) / 10)
    G <- hermitian(n, 1, c(1.2, 0.7, 1.1, 0.8, 1.0, 0.9))
    D <- C - Re(sum(diag(C))) / n * diag(n)

    X <- Variable(c(n, n), hermitian = TRUE)
    p <- Parameter()
    value(p) <- 0.5
    prob <- Problem(Minimize(sum_squares(X - p * C)), list(matrix_trace(X) == 1))
    psolve(prob, requires_grad = TRUE)

    delta(p) <- 1
    derivative(prob)
    dX <- delta(X)
    expect_equal(Re(dX), Re(D), tolerance = 1e-3)
    expect_equal(Im(dX), Im(D), tolerance = 1e-3)

    gradient(X) <- G
    backward(prob)
    expected <- sum(Re(G) * Re(D)) + sum(triu_rm(Im(G)) * triu_rm(Im(D)))
    expect_equal(as.numeric(gradient(p)), expected, tolerance = 1e-3)
  }
})

## Claims-audit decision 12 (C0622/C0636, repro 21). CvxAttr2Constr.var_backward
## (cvx_attr2constr.py:313-330) lowers a sparsity variable's gradient to its
## pattern entries; CVXR passed the full 3x3 gradient to the 3-long reduced
## variable, so P.gradient came back as recycled garbage.
## @cvxpy NONE
test_that("sparsity variable: backward and derivative stay on the pattern", {
  skip_if_not_installed("scs")
  pattern <- rbind(c(1L, 1L), c(2L, 3L), c(3L, 2L))
  x <- Variable(c(3L, 3L), sparsity = pattern)
  P <- Parameter(c(3L, 3L))
  value(P) <- matrix(1:9 / 10, 3L, 3L, byrow = TRUE)
  prob <- Problem(Minimize(sum_squares(x - P)))
  psolve(prob, requires_grad = TRUE)
  mask <- matrix(0, 3L, 3L)
  mask[pattern] <- 1
  expect_equal(as.matrix(value(x)), value(P) * mask, tolerance = 1e-5)

  G <- matrix(1:9, 3L, 3L, byrow = TRUE)
  gradient(x) <- G
  backward(prob)
  expect_equal(as.matrix(gradient(P)), G * mask, tolerance = 1e-4)

  delta(P) <- matrix(1, 3L, 3L)
  derivative(prob)
  expect_equal(as.matrix(delta(x)), mask, tolerance = 1e-4)
})
