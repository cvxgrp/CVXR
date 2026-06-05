## DQCP CVXPY 1.9.0 parity tests
## Ported from: cvxpy/tests/test_dqcp.py (CVXPY 1.9.0, worktree ~/GitHub/cvxpy-1.9)
##
## Closes the 6 test_dqcp.py gaps new in CVXPY 1.9.0 not covered by
## test-dqcp-cvxpy-parity.R / test-dqcp-parity.R / test-dqcp-properties.R:
##   test_bisection_low_zero, test_dist_ratio_dqcp, test_dqcp_power_infeasible,
##   test_length_verbose, test_sign_verbose, test_dist_ratio_verbose.
##
## The three *_verbose tests are regressions of a CVXPY 1.9.0 fix: bisect()'s
## verbose path must not evaluate lazy constraints before the bisection
## parameter t has a value (CVXPY `_lower_problem()` ran before `t.value` was
## set -> TypeError). CVXR's bisect() (bisection.R) already avoids the
## premature .lower_problem() call in its verbose block; these guard it.

# =====================================================================
# test_bisection_low_zero
# =====================================================================
## @cvxpy test_dqcp.py::TestDqcp::test_bisection_low_zero
test_that("CVXPY 1.9 parity: bisection_low_zero — low=0 must not infinite-loop", {
  skip_if_not_installed("clarabel")

  x <- Variable()
  problem <- Problem(Minimize(ceiling(x)), list(x >= 15, x <= 17))
  result <- psolve(problem, solver = "CLARABEL", qcp = TRUE, low = 0, high = 100)
  expect_equal(status(problem), "optimal")
  expect_equal(as.numeric(result), 15.0, tolerance = 1e-3)
})

# =====================================================================
# test_dist_ratio_dqcp
# =====================================================================
## @cvxpy test_dqcp.py::TestDqcp::test_dist_ratio_dqcp
test_that("CVXPY 1.9 parity: dist_ratio_dqcp — correct optimum at 0", {
  skip_if_not_installed("clarabel")

  x <- Variable(2)
  a <- c(0.0, 0.0)
  b <- c(3.0, 0.0)
  atom <- dist_ratio(x, a, b)
  prob <- Problem(Minimize(atom), list(norm(x) <= 2))
  result <- psolve(prob, solver = "CLARABEL", qcp = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(result), 0.0, tolerance = 1e-3)
})

# =====================================================================
# test_dqcp_power_infeasible
# =====================================================================
## @cvxpy test_dqcp.py::TestDqcp::test_dqcp_power_infeasible
test_that("CVXPY 1.9 parity: dqcp_power_infeasible — power(ceil(x),2) <= -5", {
  skip_if_not_installed("clarabel")

  x <- Variable(nonneg = TRUE)
  problem <- Problem(Minimize(0), list(power(ceiling(x), 2) <= -5))
  psolve(problem, solver = "CLARABEL", qcp = TRUE)
  expect_equal(status(problem), "infeasible")
})

# =====================================================================
# test_length_verbose
# =====================================================================
## @cvxpy test_dqcp.py::TestDqcp::test_length_verbose
test_that("CVXPY 1.9 parity: length_verbose — verbose=TRUE must not crash", {
  skip_if_not_installed("clarabel")

  ## CVXR's length atom is `length_expr` (base R `length()` returns 1 for
  ## the S7 Expression object, so `length()` is NOT the DQCP atom).
  x <- Variable(5)
  problem <- Problem(Minimize(length_expr(x)), list(x[1] == 2.0, x[2] == 1.0))
  result <- psolve(problem, solver = "CLARABEL", qcp = TRUE, verbose = TRUE)
  expect_equal(status(problem), "optimal")
  expect_equal(as.numeric(value(problem@objective)), 2)
})

# =====================================================================
# test_sign_verbose
# =====================================================================
## @cvxpy test_dqcp.py::TestDqcp::test_sign_verbose
test_that("CVXPY 1.9 parity: sign_verbose — verbose=TRUE must not crash", {
  skip_if_not_installed("clarabel")

  x <- Variable()
  problem <- Problem(Minimize(sign(x)), list(-2 <= x, x <= -0.5))
  result <- psolve(problem, solver = "CLARABEL", qcp = TRUE, verbose = TRUE)
  expect_equal(status(problem), "optimal")
  expect_equal(as.numeric(value(problem@objective)), -1)
})

# =====================================================================
# test_dist_ratio_verbose
# =====================================================================
## @cvxpy test_dqcp.py::TestDqcp::test_dist_ratio_verbose
test_that("CVXPY 1.9 parity: dist_ratio_verbose — verbose=TRUE must not crash", {
  skip_if_not_installed("clarabel")

  x <- Variable(2)
  a <- c(0.0, 0.0)
  b <- c(3.0, 0.0)
  problem <- Problem(Minimize(dist_ratio(x, a, b)), list(norm(x) <= 2))
  result <- psolve(problem, solver = "CLARABEL", qcp = TRUE, verbose = TRUE)
  expect_equal(status(problem), "optimal")
  expect_equal(as.numeric(result), 0.0, tolerance = 1e-3)
})
