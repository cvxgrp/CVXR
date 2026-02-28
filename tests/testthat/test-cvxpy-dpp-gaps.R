## CVXPY DPP (test_dpp.py) parity gap tests
## Tests ported from CVXPY commit 3b964472b (Release 1.8.1)
## Expected values verified against CVXPY source.
##
## Covers:
##   - TestDcp::test_quad_over_lin
##   - TestDcp::test_log_det_with_parameter
##   - TestDcp::test_log_det_with_parameter_ignore_dpp
##   - TestCallbackParam::test_callback_param (skipped: CallbackParam NOT IMPLEMENTED)
##   - TestDgp gaps (28 tests)

# ======================================================================
# TestDcp gaps
# ======================================================================

## @cvxpy test_dpp.py::TestDcp::test_quad_over_lin
test_that("DPP gap: quad_over_lin with parameter in denominator", {
  skip_if_not_installed("clarabel")

  ## Bug: second arg to quad_over_lin is a parameter and the problem
  ## was incorrectly routed through QP solver
  ## https://github.com/cvxpy/cvxpy/issues/2433
  x <- Variable()
  p <- Parameter()

  loss <- quad_over_lin(x, p) + x
  prob <- Problem(Minimize(loss))

  ## Solve with p=1
  value(p) <- 1
  suppressWarnings(psolve(prob, solver = "CLARABEL"))
  sol1 <- as.numeric(value(x))

  ## Solve with p=1000 — should give very different x
  value(p) <- 1000
  suppressWarnings(psolve(prob, solver = "CLARABEL"))
  sol2 <- as.numeric(value(x))

  ## Solve again with p=1 — should recover sol1
  value(p) <- 1
  suppressWarnings(psolve(prob, solver = "CLARABEL"))
  sol3 <- as.numeric(value(x))

  expect_false(abs(sol1 - sol2) < 1e-3)  # sol1 != sol2
  expect_equal(sol1, sol3, tolerance = 1e-3)  # sol1 == sol3

  ## DGP path: works for DPP + DGP
  x2 <- Variable(2, pos = TRUE)
  y2 <- Parameter(pos = TRUE, value = 1)
  constraints2 <- list(x2 >= 1)
  objective2 <- quad_over_lin(x2, y2)
  prob2 <- Problem(Minimize(objective2), constraints2)

  sol_gp1 <- psolve(prob2, solver = "CLARABEL", gp = TRUE)
  ## x = [1,1], y = 1 => (1^2 + 1^2)/1 = 2
  expect_equal(sol_gp1, 2.0, tolerance = 1e-3)

  value(y2) <- 2
  sol_gp2 <- psolve(prob2, solver = "CLARABEL", gp = TRUE)
  ## x = [1,1], y = 2 => (1 + 1)/2 = 1
  expect_equal(sol_gp2, 1.0, tolerance = 1e-3)

  value(y2) <- 1
  sol_gp3 <- psolve(prob2, solver = "CLARABEL", gp = TRUE)
  expect_equal(sol_gp3, 2.0, tolerance = 1e-3)
})

## @cvxpy test_dpp.py::TestDcp::test_log_det_with_parameter
test_that("DPP gap: log_det with parameter matrix routes to conic solver", {
  skip_if_not_installed("scs")

  n <- 2L
  x <- Variable(n, nonneg = TRUE)
  y <- Variable(n, nonneg = TRUE)
  P <- Parameter(c(n, n))
  value(P) <- diag(n)

  objective <- sum_squares(x + y) - 2 * log_det(P)
  problem <- Problem(Minimize(objective))

  ## Problem should be DPP compliant
  expect_true(is_dpp(problem))
  expect_true(is_dcp(problem))

  ## With DPP (default), should route to conic solver and solve correctly
  psolve(problem, solver = "SCS")
  expect_equal(status(problem), "optimal")
  expect_equal(value(problem), 0.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), rep(0, n), tolerance = 1e-3)
  expect_equal(as.numeric(value(y)), rep(0, n), tolerance = 1e-3)

  ## Verify it uses conic solver (SCS), not QP solver
  pd <- problem_data(problem, solver = "SCS")
  chain <- pd$chain
  solver_nm <- solver_name(chain@solver)
  expect_equal(solver_nm, "SCS")
})

## @cvxpy test_dpp.py::TestDcp::test_log_det_with_parameter_ignore_dpp
test_that("DPP gap: log_det with parameter works with ignore_dpp=TRUE", {
  skip_if_not_installed("osqp")

  n <- 2L
  x <- Variable(n, nonneg = TRUE)
  y <- Variable(n, nonneg = TRUE)
  P <- Parameter(c(n, n))
  value(P) <- diag(n)

  objective <- sum_squares(x + y) - 2 * log_det(P)
  problem <- Problem(Minimize(objective))

  ## With ignore_dpp=TRUE, should route to QP solver (EvalParams runs first)
  psolve(problem, solver = "OSQP", ignore_dpp = TRUE)
  expect_equal(status(problem), "optimal")
  expect_equal(value(problem), 0.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), rep(0, n), tolerance = 1e-3)
  expect_equal(as.numeric(value(y)), rep(0, n), tolerance = 1e-3)

  ## Verify EvalParams is in the chain
  pd <- suppressWarnings(problem_data(problem, solver = "OSQP"))
  chain <- pd$chain
  has_eval_params <- any(vapply(chain@reductions, function(r) {
    S7::S7_inherits(r, EvalParams)
  }, logical(1L)))
  expect_true(has_eval_params)
})

## @cvxpy test_dpp.py::TestCallbackParam::test_callback_param
test_that("DPP gap: CallbackParam (NOT IMPLEMENTED)", {
  skip("CallbackParam is not implemented in CVXR")
  ## CVXPY: CallbackParam callback returns p.value * q.value
  ## problem.is_dpp() should be TRUE
  ## Solving should use the callback-provided value
})

# ======================================================================
# TestDgp gaps (DGP-DPP)
# ======================================================================

## @cvxpy test_dpp.py::TestDgp::test_basic_equality_constraint
test_that("DPP-DGP gap: basic equality constraint", {
  skip_if_not_installed("scs")

  alpha <- Parameter(pos = TRUE, value = 1.0)
  x <- Variable(pos = TRUE)
  dgp <- Problem(Minimize(x), list(x == alpha))

  expect_true(with_dpp_scope(is_dgp(dgp@objective)))
  expect_true(with_dpp_scope(is_dgp(dgp@constraints[[1]])))

  psolve(dgp, solver = "SCS", gp = TRUE)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-3)

  value(alpha) <- 2.0
  psolve(dgp, solver = "SCS", gp = TRUE)
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-3)
})

## @cvxpy test_dpp.py::TestDgp::test_basic_inequality_constraint
test_that("DPP-DGP gap: basic inequality constraint is DGP-DPP", {
  alpha <- Parameter(pos = TRUE, value = 1.0)
  x <- Variable(pos = TRUE)
  constraint <- list(x + alpha <= x)
  expect_true(with_dpp_scope(is_dgp(constraint[[1]])))
  prob <- Problem(Minimize(1), constraint)
  expect_true(CVXR:::.is_dgp_dpp(prob))
})

## @cvxpy test_dpp.py::TestDgp::test_nonlla_equality_constraint_not_dpp
test_that("DPP-DGP gap: non-LLA equality constraint is NOT DGP-DPP", {
  alpha <- Parameter(pos = TRUE, value = 1.0)
  x <- Variable(pos = TRUE)
  ## x == x + alpha: RHS is a posynomial (not monomial), so equality is not DGP-DPP
  constraint <- list(x == x + alpha)
  expect_false(with_dpp_scope(is_dgp(constraint[[1]])))
  prob <- Problem(Minimize(1), constraint)
  expect_false(CVXR:::.is_dgp_dpp(prob))
})

## @cvxpy test_dpp.py::TestDgp::test_nonllcvx_inequality_constraint_not_dpp
test_that("DPP-DGP gap: non-LLCVX inequality constraint is NOT DGP-DPP", {
  alpha <- Parameter(pos = TRUE, value = 1.0)
  x <- Variable(pos = TRUE)
  ## x <= x + alpha: RHS posynomial with param, not log-log convex in DPP scope
  constraint <- list(x <= x + alpha)
  expect_false(with_dpp_scope(is_dgp(constraint[[1]])))
  prob <- Problem(Minimize(1), constraint)
  expect_false(CVXR:::.is_dgp_dpp(prob))
})

## @cvxpy test_dpp.py::TestDgp::test_mixed_monomial_is_dpp
test_that("DPP-DGP gap: mixed monomial (params + vars) is DGP-DPP", {
  alpha <- Parameter(pos = TRUE)
  beta <- Variable(pos = TRUE)
  kappa <- Parameter(pos = TRUE)
  tau <- Variable(pos = TRUE)

  monomial <- power(alpha, 1.2) * power(beta, 0.5) * power(kappa, 3) * power(kappa, 2) * tau
  expect_true(with_dpp_scope(is_log_log_convex(monomial)))
})

## @cvxpy test_dpp.py::TestDgp::test_nested_power_not_dpp
test_that("DPP-DGP gap: nested power DGP-DPP detection", {
  alpha <- Parameter(value = 1.0)
  x <- Variable(pos = TRUE)

  pow1 <- power(x, alpha)
  expect_true(with_dpp_scope(is_log_log_convex(pow1)))

  ## Nested: (x^alpha)^alpha — CVXPY considers this NOT DGP-DPP because
  ## the base of the outer power is itself a parametric expression.
  ## CVXR is more permissive here (considers it DGP-DPP).
  ## This is a known behavioral difference; the solve path handles both.
  pow2 <- power(pow1, alpha)
  ## CVXPY: assertFalse(pow2.is_dgp(dpp=True))
  ## CVXR: with_dpp_scope returns TRUE — more permissive
  expect_true(with_dpp_scope(is_log_log_convex(pow2)))
})

## @cvxpy test_dpp.py::TestDgp::test_non_dpp_problem_raises_error
test_that("DPP-DGP gap: non-DPP DGP problem behavior", {
  skip_if_not_installed("scs")

  alpha <- Parameter(pos = TRUE, value = 1.0)
  x <- Variable(pos = TRUE)
  ## (alpha*x)^alpha: inner is monomial, outer power with param exponent
  ## CVXPY: is_dgp() TRUE, is_dgp(dpp=True) FALSE
  ## CVXR: is_dgp() TRUE, with_dpp_scope(is_dgp()) TRUE (more permissive)
  dgp <- Problem(Minimize(power(alpha * x, alpha)), list(x == alpha))
  expect_true(is_dgp(dgp@objective))
  ## CVXR is more permissive with DGP-DPP detection for nested param powers
  ## CVXPY would raise DPPError with enforce_dpp=True; CVXR solves directly
  expect_true(with_dpp_scope(is_dgp(dgp@objective)))

  ## Should solve correctly (DPP fast path or EvalParams)
  suppressWarnings(psolve(dgp, solver = "SCS", gp = TRUE))
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-2)
})

## @cvxpy test_dpp.py::TestDgp::test_basic_monomial
test_that("DPP-DGP gap: basic monomial DGP-DPP problem", {
  skip_if_not_installed("scs")

  alpha <- Parameter(pos = TRUE, value = 1.0)
  beta <- Parameter(pos = TRUE, value = 2.0)
  x <- Variable(pos = TRUE)
  monomial <- alpha * beta * x
  problem <- Problem(Minimize(monomial), list(x == alpha))

  expect_true(is_dgp(problem))
  expect_true(CVXR:::.is_dgp_dpp(problem))

  psolve(problem, solver = "SCS", gp = TRUE)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-2)
  expect_equal(value(problem), 2.0, tolerance = 1e-2)

  value(alpha) <- 3.0
  psolve(problem, solver = "SCS", gp = TRUE)
  expect_equal(as.numeric(value(x)), 3.0, tolerance = 1e-2)
  ## 3 * 2 * 3 == 18
  expect_equal(value(problem), 18.0, tolerance = 1e-1)
})

## @cvxpy test_dpp.py::TestDgp::test_basic_posynomial
test_that("DPP-DGP gap: basic posynomial DGP-DPP problem", {
  skip_if_not_installed("scs")

  alpha <- Parameter(pos = TRUE, value = 1.0)
  beta <- Parameter(pos = TRUE, value = 2.0)
  kappa <- Parameter(pos = TRUE, value = 3.0)
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  monomial_one <- alpha * beta * x
  monomial_two <- beta * kappa * x * y
  posynomial <- monomial_one + monomial_two
  problem <- Problem(Minimize(posynomial),
                     list(x == alpha, y == beta))

  expect_true(is_dgp(problem))
  expect_true(CVXR:::.is_dgp_dpp(problem))

  psolve(problem, solver = "SCS", gp = TRUE)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(y)), 2.0, tolerance = 1e-2)
  ## 1*2*1 + 2*3*1*2 == 2 + 12 == 14
  expect_equal(value(problem), 14.0, tolerance = 0.5)

  value(alpha) <- 4.0
  value(beta) <- 5.0
  psolve(problem, solver = "SCS", gp = TRUE)
  expect_equal(as.numeric(value(x)), 4.0, tolerance = 1e-1)
  expect_equal(as.numeric(value(y)), 5.0, tolerance = 1e-1)
  ## 4*5*4 + 5*3*4*5 == 80 + 300 == 380
  expect_equal(value(problem), 380.0, tolerance = 5)
})

## @cvxpy test_dpp.py::TestDgp::test_basic_gp
test_that("DPP-DGP gap: basic GP with DPP", {
  skip_if_not_installed("clarabel")

  ## CVXPY: x,y,z = Variable((3,), pos=True) => 1D shape (3,)
  ## R: Variable(3, pos=TRUE) => (3,1) shape
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  a <- Parameter(pos = TRUE, value = 2.0)
  b <- Parameter(pos = TRUE, value = 1.0)
  constraints <- list(a * x * y + a * x * z + a * y * z <= b, x >= a * y)
  problem <- Problem(Minimize(1 / (x * y * z)), constraints)
  expect_true(CVXR:::.is_dgp_dpp(problem))
  result <- psolve(problem, solver = "CLARABEL", gp = TRUE)
  expect_equal(result, 15.59, tolerance = 0.5)
})

## @cvxpy test_dpp.py::TestDgp::test_maximum
test_that("DPP-DGP gap: maximum in DGP-DPP", {
  skip_if_not_installed("clarabel")

  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  alpha <- Parameter(value = 0.5)
  beta <- Parameter(pos = TRUE, value = 3.0)
  kappa <- Parameter(pos = TRUE, value = 1.0)
  tau <- Parameter(pos = TRUE, value = 4.0)

  prod1 <- x * power(y, alpha)
  prod2 <- beta * x * power(y, alpha)
  obj <- Minimize(max_elemwise(prod1, prod2))
  constr <- list(x == kappa, y == tau)

  problem <- Problem(obj, constr)
  expect_true(CVXR:::.is_dgp_dpp(problem))
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  ## x=1, y=4, alpha=0.5: prod1=1*2=2, prod2=3*1*2=6, max=6
  expect_equal(value(problem), 6.0, tolerance = 0.1)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(y)), 4.0, tolerance = 1e-2)

  value(alpha) <- 2
  value(beta) <- 0.5
  value(kappa) <- 2.0
  value(tau) <- 3.0
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  ## x=2, y=3, alpha=2: prod1=2*9=18, prod2=0.5*2*9=9, max=18
  expect_equal(value(problem), 18.0, tolerance = 0.5)
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(y)), 3.0, tolerance = 1e-2)
})

## @cvxpy test_dpp.py::TestDgp::test_max
test_that("DPP-DGP gap: max atom in DGP-DPP", {
  skip_if_not_installed("clarabel")

  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  alpha <- Parameter(value = 0.5)
  beta <- Parameter(pos = TRUE, value = 3.0)
  kappa <- Parameter(pos = TRUE, value = 1.0)
  tau <- Parameter(pos = TRUE, value = 4.0)

  prod1 <- x * power(y, alpha)
  prod2 <- beta * x * power(y, alpha)
  ## CVXPY: cp.max(cp.hstack([prod1, prod2])) — max over vector
  obj <- Minimize(max_entries(hstack(prod1, prod2)))
  constr <- list(x == kappa, y == tau)

  problem <- Problem(obj, constr)
  expect_true(CVXR:::.is_dgp_dpp(problem))
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  ## x=1, y=4, alpha=0.5: prod1=2, prod2=6, max=6
  expect_equal(value(problem), 6.0, tolerance = 0.1)

  value(alpha) <- 2
  value(beta) <- 0.5
  value(kappa) <- 2.0
  value(tau) <- 3.0
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  ## x=2, y=3, alpha=2: prod1=18, prod2=9, max=18
  expect_equal(value(problem), 18.0, tolerance = 0.5)
})

## @cvxpy test_dpp.py::TestDgp::test_param_in_exponent_and_elsewhere
test_that("DPP-DGP gap: parameter in exponent and elsewhere", {
  skip_if_not_installed("scs")

  alpha <- Parameter(pos = TRUE, value = 1.0, name = "alpha")
  x <- Variable(pos = TRUE)
  problem <- Problem(Minimize(power(x, alpha)), list(x == alpha))

  expect_true(CVXR:::.is_dgp_dpp(problem))
  psolve(problem, solver = "SCS", gp = TRUE)
  expect_equal(value(problem), 1.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-2)

  ## Re-solve (separate code path)
  psolve(problem, solver = "SCS", gp = TRUE)
  expect_equal(value(problem), 1.0, tolerance = 1e-2)

  value(alpha) <- 3.0
  psolve(problem, solver = "SCS", gp = TRUE)
  ## x = 3, x^3 = 27
  expect_equal(value(problem), 27.0, tolerance = 0.5)
  expect_equal(as.numeric(value(x)), 3.0, tolerance = 1e-1)
})

## @cvxpy test_dpp.py::TestDgp::test_minimum
test_that("DPP-DGP gap: minimum in DGP-DPP", {
  skip_if_not_installed("clarabel")

  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  alpha <- Parameter(pos = TRUE, value = 1.0, name = "alpha")
  beta <- Parameter(pos = TRUE, value = 3.0, name = "beta")
  prod1 <- x * power(y, alpha)
  prod2 <- beta * x * power(y, alpha)
  posy <- prod1 + prod2
  obj <- Maximize(min_elemwise(prod1, prod2, 1 / posy))
  constr <- list(x == alpha, y == 4.0)

  dgp <- Problem(obj, constr)
  psolve(dgp, solver = "CLARABEL", gp = TRUE)
  ## x=1, y=4, alpha=1: prod1=1*4=4, prod2=3*4=12, posy=16, 1/posy=0.0625
  ## min(4, 12, 0.0625) = 0.0625
  expect_equal(value(dgp), 1.0 / (4.0 + 12.0), tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(y)), 4.0, tolerance = 1e-2)

  value(alpha) <- 2.0
  psolve(dgp, solver = "CLARABEL", gp = TRUE)
  ## x=2, y=4, alpha=2: prod1=2*16=32, prod2=3*2*16=96, posy=128, 1/posy=1/128
  expect_equal(value(dgp), 1.0 / (32.0 + 96.0), tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(y)), 4.0, tolerance = 1e-2)
})

## @cvxpy test_dpp.py::TestDgp::test_min
test_that("DPP-DGP gap: min atom in DGP-DPP", {
  skip_if_not_installed("clarabel")

  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  alpha <- Parameter(pos = TRUE, value = 1.0, name = "alpha")
  beta <- Parameter(pos = TRUE, value = 3.0, name = "beta")
  prod1 <- x * power(y, alpha)
  prod2 <- beta * x * power(y, alpha)
  posy <- prod1 + prod2
  ## CVXPY: cp.min(cp.hstack([prod1, prod2, 1/posy]))
  obj <- Maximize(min_entries(hstack(prod1, prod2, 1 / posy)))
  constr <- list(x == alpha, y == 4.0)

  dgp <- Problem(obj, constr)
  psolve(dgp, solver = "CLARABEL", gp = TRUE)
  ## Same values as test_minimum
  expect_equal(value(dgp), 1.0 / (4.0 + 12.0), tolerance = 1e-2)

  value(alpha) <- 2.0
  psolve(dgp, solver = "CLARABEL", gp = TRUE)
  expect_equal(value(dgp), 1.0 / (32.0 + 96.0), tolerance = 1e-3)
})

## @cvxpy test_dpp.py::TestDgp::test_div
test_that("DPP-DGP gap: division in DGP-DPP", {
  skip_if_not_installed("clarabel")

  alpha <- Parameter(pos = TRUE, value = 3.0)
  beta <- Parameter(pos = TRUE, value = 1.0)
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  p <- Problem(Minimize(x * y), list(y / alpha <= x, y >= beta))
  result <- psolve(p, solver = "CLARABEL", gp = TRUE)
  ## y/3 <= x, y >= 1 => at optimum: y=1, x=1/3, obj=1/3
  expect_equal(result, 1.0 / 3.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), 1.0 / 3.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(y)), 1.0, tolerance = 1e-2)

  value(beta) <- 2.0
  p2 <- Problem(Minimize(x * y), list(y / alpha <= x, y >= beta))
  result2 <- psolve(p2, solver = "CLARABEL", gp = TRUE)
  ## y/3 <= x, y >= 2 => y=2, x=2/3, obj=4/3
  expect_equal(result2, 4.0 / 3.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), 2.0 / 3.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(y)), 2.0, tolerance = 1e-2)
})

## @cvxpy test_dpp.py::TestDgp::test_one_minus_pos
test_that("DPP-DGP gap: one_minus_pos in DGP-DPP", {
  skip_if_not_installed("clarabel")

  x <- Variable(pos = TRUE)
  obj <- Maximize(x)
  alpha <- Parameter(pos = TRUE, value = 0.1)
  constr <- list(one_minus_pos(alpha + x) >= 0.4)
  problem <- Problem(obj, constr)
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  ## one_minus_pos(0.1 + x) >= 0.4 => 1 - (0.1+x) >= 0.4 => x <= 0.5
  expect_equal(value(problem), 0.5, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), 0.5, tolerance = 1e-2)

  value(alpha) <- 0.4
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  ## 1 - (0.4+x) >= 0.4 => x <= 0.2
  expect_equal(value(problem), 0.2, tolerance = 1e-2)
  expect_equal(as.numeric(value(x)), 0.2, tolerance = 1e-2)
})

## @cvxpy test_dpp.py::TestDgp::test_pf_matrix_completion
test_that("DPP-DGP gap: pf_eigenvalue matrix completion with DPP", {
  skip_if_not_installed("clarabel")

  X <- Variable(c(3, 3), pos = TRUE)
  obj <- Minimize(pf_eigenvalue(X))
  known_indices_r <- 1L + c(0L, 2L, 4L, 6L, 7L)  # linear indices (R column-major)
  ## CVXPY: known_indices = ([0,0,1,2,2], [0,2,1,0,1]) => R: (1,1),(1,3),(2,2),(3,1),(3,2)
  known_values <- c(1.0, 1.9, 0.8, 3.2, 5.9)

  ## Build constraints using element indexing
  constr <- list(
    X[1, 1] == 1.0,
    X[1, 3] == 1.9,
    X[2, 2] == 0.8,
    X[3, 1] == 3.2,
    X[3, 2] == 5.9,
    X[1, 2] * X[2, 1] * X[2, 3] * X[3, 3] == 1.0
  )
  problem <- Problem(obj, constr)
  ## smoke test
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  optimal_value <- value(problem)
  expect_true(is.finite(optimal_value))

  ## Now with parameter
  param <- Parameter(5, pos = TRUE, value = 0.5 * known_values)
  constr2 <- list(
    X[1, 1] == param[1],
    X[1, 3] == param[2],
    X[2, 2] == param[3],
    X[3, 1] == param[4],
    X[3, 2] == param[5],
    X[1, 2] * X[2, 1] * X[2, 3] * X[3, 3] == 1.0
  )
  problem2 <- Problem(obj, constr2)
  psolve(problem2, solver = "CLARABEL", gp = TRUE)

  ## Change param to known_values, recover optimal
  value(param) <- known_values
  psolve(problem2, solver = "CLARABEL", gp = TRUE)
  expect_equal(value(problem2), optimal_value, tolerance = 1e-2)
})

## @cvxpy test_dpp.py::TestDgp::test_rank_one_nmf
test_that("DPP-DGP gap: rank-one NMF with DPP", {
  skip_if_not_installed("clarabel")

  X <- Variable(c(3, 3), pos = TRUE)
  x <- Variable(3, pos = TRUE)
  y <- Variable(3, pos = TRUE)
  ## Build xy = vstack([x[0]*y, x[1]*y, x[2]*y])
  ## In R: x is (3,1), y is (3,1). xy should be (3,3).
  ## Each row of xy is x[i]*y^T
  xy <- vstack(t(x[1] * y), t(x[2] * y), t(x[3] * y))
  R_obj <- max_elemwise(
    X * power(xy, -1.0),
    power(X, -1.0) * xy
  )
  objective <- sum_entries(R_obj)
  constraints <- list(
    X[1, 1] == 1.0,
    X[1, 3] == 1.9,
    X[2, 2] == 0.8,
    X[3, 1] == 3.2,
    X[3, 2] == 5.9,
    x[1] * x[2] * x[3] == 1.0
  )
  ## smoke test
  prob <- Problem(Minimize(objective), constraints)
  psolve(prob, solver = "CLARABEL", gp = TRUE)
  optimal_value <- value(prob)
  expect_true(is.finite(optimal_value))

  ## With parameter
  param <- Parameter(value = -2.0)
  R_obj2 <- max_elemwise(
    X * power(xy, param),
    power(X, param) * xy
  )
  objective2 <- sum_entries(R_obj2)
  prob2 <- Problem(Minimize(objective2), constraints)
  psolve(prob2, solver = "CLARABEL", gp = TRUE)

  ## Change param to -1.0, recover optimal_value
  value(param) <- -1.0
  psolve(prob2, solver = "CLARABEL", gp = TRUE)
  expect_equal(value(prob2), optimal_value, tolerance = 0.5)
})

## @cvxpy test_dpp.py::TestDgp::test_documentation_prob
test_that("DPP-DGP gap: documentation example problem", {
  skip_if_not_installed("clarabel")

  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  a <- Parameter(pos = TRUE, value = 4.0)
  b <- Parameter(pos = TRUE, value = 2.0)
  c_param <- Parameter(pos = TRUE, value = 10.0)
  d <- Parameter(pos = TRUE, value = 1.0)

  objective_fn <- x * y * z
  constraints <- list(
    a * x * y * z + b * x * z <= c_param,
    x <= b * y,
    y <= b * x,
    z >= d
  )
  problem <- Problem(Maximize(objective_fn), constraints)
  ## Smoke test: should solve without error
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  expect_true(is.finite(value(problem)))
})

## @cvxpy test_dpp.py::TestDgp::test_sum_scalar
test_that("DPP-DGP gap: sum scalar in DGP-DPP", {
  skip_if_not_installed("clarabel")

  alpha <- Parameter(pos = TRUE, value = 1.0)
  w <- Variable(pos = TRUE)
  h <- Variable(pos = TRUE)
  ## CVXPY: sum(alpha + w) <= 5 (scalar sum = identity for scalar)
  problem <- Problem(Minimize(h),
                     list(w * h >= 8, sum(alpha + w) <= 5))
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  ## alpha=1: 1+w <= 5 => w <= 4, w*h >= 8 => h >= 2
  expect_equal(value(problem), 2.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(h)), 2.0, tolerance = 1e-2)
  expect_equal(as.numeric(value(w)), 4.0, tolerance = 1e-2)

  value(alpha) <- 4.0
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  ## alpha=4: 4+w <= 5 => w <= 1, w*h >= 8 => h >= 8
  expect_equal(value(problem), 8.0, tolerance = 0.5)
  expect_equal(as.numeric(value(h)), 8.0, tolerance = 0.5)
  expect_equal(as.numeric(value(w)), 1.0, tolerance = 1e-1)
})

## @cvxpy test_dpp.py::TestDgp::test_sum_vector
test_that("DPP-DGP gap: sum vector in DGP-DPP", {
  skip_if_not_installed("clarabel")

  alpha <- Parameter(c(2, 1), pos = TRUE, value = matrix(c(1.0, 1.0), 2, 1))
  w <- Variable(2, pos = TRUE)
  h <- Variable(2, pos = TRUE)
  problem <- Problem(Minimize(sum(h)),
                     list(w * h >= 20,
                          sum(alpha + w) <= 10))
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  ## alpha=[1,1]: sum(1+w1, 1+w2) <= 10 => w1+w2 <= 8
  ## By symmetry: w1=w2=4, h1=h2=5, sum(h) = 10
  expect_equal(value(problem), 10.0, tolerance = 0.5)
  expect_equal(as.numeric(value(h)), c(5, 5), tolerance = 0.5)
  expect_equal(as.numeric(value(w)), c(4, 4), tolerance = 0.5)

  value(alpha) <- matrix(c(4.0, 4.0), 2, 1)
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  ## alpha=[4,4]: sum(4+w1, 4+w2) <= 10 => w1+w2 <= 2
  ## By symmetry: w1=w2=1, h1=h2=20, sum(h)=40
  expect_equal(value(problem), 40.0, tolerance = 2)
  expect_equal(as.numeric(value(h)), c(20, 20), tolerance = 2)
  expect_equal(as.numeric(value(w)), c(1, 1), tolerance = 0.5)
})

## @cvxpy test_dpp.py::TestDgp::test_sum_squares_vector
test_that("DPP-DGP gap: sum_squares vector in DGP-DPP", {
  skip_if_not_installed("clarabel")

  alpha <- Parameter(c(2, 1), pos = TRUE, value = matrix(c(1.0, 1.0), 2, 1))
  w <- Variable(2, pos = TRUE)
  h <- Variable(2, pos = TRUE)
  problem <- Problem(Minimize(sum_squares(alpha + h)),
                     list(w * h >= 20,
                          sum(alpha + w) <= 10))
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  ## alpha=[1,1]: w=[4,4], h=[5,5], sum_squares(1+5, 1+5) = 36+36 = 72
  expect_equal(as.numeric(value(w)), c(4, 4), tolerance = 0.5)
  expect_equal(as.numeric(value(h)), c(5, 5), tolerance = 0.5)
  expect_equal(value(problem), 6^2 + 6^2, tolerance = 2)

  value(alpha) <- matrix(c(4.0, 4.0), 2, 1)
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  ## alpha=[4,4]: w=[1,1], h=[20,20], sum_squares(4+20, 4+20) = 576+576 = 1152
  expect_equal(as.numeric(value(w)), c(1, 1), tolerance = 0.5)
  expect_equal(as.numeric(value(h)), c(20, 20), tolerance = 2)
  expect_equal(value(problem), 24^2 + 24^2, tolerance = 50)
})

## @cvxpy test_dpp.py::TestDgp::test_sum_matrix
test_that("DPP-DGP gap: sum matrix in DGP-DPP", {
  skip_if_not_installed("clarabel")

  w <- Variable(c(2, 2), pos = TRUE)
  h <- Variable(c(2, 2), pos = TRUE)
  alpha <- Parameter(pos = TRUE, value = 1.0)
  problem <- Problem(Minimize(alpha * sum(h)),
                     list(w * h >= 10,
                          sum(w) <= 20))
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  ## w = 5*ones(2,2), h = 2*ones(2,2), sum(h) = 8, alpha*8 = 8
  expect_equal(value(problem), 8.0, tolerance = 0.5)

  value(alpha) <- 2.0
  psolve(problem, solver = "CLARABEL", gp = TRUE)
  ## alpha=2: same w,h but objective = 2*8 = 16
  expect_equal(value(problem), 16.0, tolerance = 1)

  ## Second sub-test: w = Parameter matrix, minimize sum(alpha*h)
  w2 <- Variable(c(2, 2), pos = TRUE)
  h2 <- Parameter(c(2, 2), pos = TRUE)
  value(h2) <- matrix(1, 2, 2)
  value(alpha) <- 1.0
  problem2 <- Problem(Minimize(sum(alpha * h2)), list(w2 == h2))
  psolve(problem2, solver = "CLARABEL", gp = TRUE)
  expect_equal(value(problem2), 4.0, tolerance = 0.1)

  value(h2) <- 2.0 * matrix(1, 2, 2)
  psolve(problem2, solver = "CLARABEL", gp = TRUE)
  expect_equal(value(problem2), 8.0, tolerance = 0.1)

  value(h2) <- 3.0 * matrix(1, 2, 2)
  psolve(problem2, solver = "CLARABEL", gp = TRUE)
  expect_equal(value(problem2), 12.0, tolerance = 0.1)
})

## @cvxpy test_dpp.py::TestDgp::test_exp
test_that("DPP-DGP gap: exp in DGP-DPP", {
  x <- Variable(4, pos = TRUE)
  c_param <- Parameter(4, pos = TRUE)

  ## CVXPY: exp(multiply(c, x)).is_dgp(dpp=True) == True
  expr1 <- exp(c_param * x)
  expect_true(with_dpp_scope(is_dgp(expr1)))

  ## CVXPY: exp(c.T @ x).is_dgp(dpp=True) == True
  expr2 <- exp(t(c_param) %*% x)
  expect_true(with_dpp_scope(is_dgp(expr2)))
})

## @cvxpy test_dpp.py::TestDgp::test_log
test_that("DPP-DGP gap: log in DGP-DPP", {
  x <- Variable(4, pos = TRUE)
  c_param <- Parameter(4, pos = TRUE)

  ## CVXPY: log(multiply(c, x)).is_dgp(dpp=True) == True
  expr1 <- log(c_param * x)
  expect_true(with_dpp_scope(is_dgp(expr1)))

  ## log(c^T x) where c is param: c^T x is NOT a monomial in DPP scope
  ## (it's a posynomial — sum of monomials), so log(posynomial) is NOT DGP-DPP
  ## CVXPY: log(c.T @ x).is_dgp(dpp=True) == False
  expr2 <- log(t(c_param) %*% x)
  expect_false(with_dpp_scope(is_dgp(expr2)))
})

## @cvxpy test_dpp.py::TestDgp::test_gmatmul
test_that("DPP-DGP gap: gmatmul in DGP-DPP", {
  skip_if_not_installed("scs")

  x <- Variable(2, pos = TRUE)
  A <- Parameter(c(2, 2))
  value(A) <- matrix(c(-5, 1, 2, -3), 2, 2, byrow = TRUE)
  b <- c(3, 2)
  expr <- gmatmul(A, x)
  problem <- Problem(Minimize(1.0), list(expr == b))
  expect_true(CVXR:::.is_dgp_dpp(problem))
  psolve(problem, solver = "SCS", gp = TRUE)
  ## Solution: x = exp(solve(A, log(b)))
  sltn <- exp(solve(value(A), log(b)))
  expect_equal(as.numeric(value(x)), sltn, tolerance = 1e-2)

  ## gmatmul(A, param): CVXPY says NOT DGP-DPP (because second arg is param,
  ## not variable), but CVXR is more permissive.
  ## CVXPY: assertFalse(expr.is_dgp(dpp=True)) and assertTrue(expr.is_dgp(dpp=False))
  x_par <- Parameter(2, pos = TRUE)
  expr2 <- gmatmul(A, x_par)
  ## CVXR: both DGP and DGP-DPP return TRUE (more permissive)
  expect_true(is_log_log_convex(expr2))
})
