## Port of cvxpy/tests/test_atoms.py partial_optimize-related tests.

## @cvxpy test_atoms.py::TestAtoms::test_partial_optimize_dcp
test_that("partial_optimize: DCP curvature", {
  dims <- 3L
  x <- Variable(dims)
  t <- Variable(dims)

  p_min <- Problem(Minimize(sum_entries(t)), list(-t <= x, x <= t))
  g <- partial_optimize(p_min, list(t), list(x))
  expect_true(is_convex(g))
  expect_false(is_concave(g))

  p_max <- Problem(Maximize(sum_entries(t)), list(-t <= x, x <= t))
  g <- partial_optimize(p_max, list(t), list(x))
  expect_true(is_concave(g))
  expect_false(is_convex(g))

  ## Non-DCP inner -> neither convex nor concave.
  p_bad <- Problem(Maximize(t[1]^2), list(-t <= x, x <= t))
  g <- partial_optimize(p_bad, list(t), list(x))
  expect_false(is_convex(g))
  expect_false(is_concave(g))
})

## @cvxpy test_atoms.py::TestAtoms::test_partial_optimize_eval_1norm
test_that("partial_optimize: value() evaluates inner solve", {
  skip_if_not_installed("clarabel")
  dims <- 3L
  x <- Variable(dims)
  t <- Variable(dims)
  xval <- rep(-5, dims)

  ## Baseline: direct 1-norm in epigraph form.
  p1 <- Problem(Minimize(sum_entries(t)), list(-t <= xval, xval <= t))
  psolve(p1, solver = "CLARABEL")

  ## partial_optimize with both opt + dont_opt explicit.
  p2 <- Problem(Minimize(sum_entries(t)), list(-t <= x, x <= t))
  g  <- partial_optimize(p2, list(t), list(x))
  p3 <- Problem(Minimize(g), list(x == xval))
  psolve(p3, solver = "CLARABEL")
  expect_equal(value(p1), value(p3), tolerance = 1e-3)

  ## opt_vars only (dont_opt_vars filled by complement).
  g  <- partial_optimize(p2, opt_vars = list(t))
  p3 <- Problem(Minimize(g), list(x == xval))
  psolve(p3, solver = "CLARABEL")
  expect_equal(value(p1), value(p3), tolerance = 1e-3)

  ## dont_opt_vars only.
  g  <- partial_optimize(p2, dont_opt_vars = list(x))
  p3 <- Problem(Minimize(g), list(x == xval))
  psolve(p3, solver = "CLARABEL")
  expect_equal(value(p1), value(p3), tolerance = 1e-3)
})

## @cvxpy test_atoms.py::TestAtoms::test_partial_optimize_eval_1norm (errors)
test_that("partial_optimize: error cases (CVXPY parity)", {
  dims <- 3L
  x <- Variable(dims)
  t <- Variable(dims)
  p2 <- Problem(Minimize(sum_entries(t)), list(-t <= x, x <= t))

  ## Neither opt_vars nor dont_opt_vars.
  expect_error(partial_optimize(p2), regexp = "neither")

  ## Both specified but not covering all vars (here opt_vars empty,
  ## dont_opt_vars = list(x) misses t).
  expect_error(partial_optimize(p2, list(), list(x)),
               regexp = "every variable")
})

## @cvxpy test_atoms.py::TestAtoms::test_partial_optimize_min_1norm
test_that("partial_optimize: embedding inside outer Problem (1-norm)", {
  skip_if_not_installed("clarabel")
  dims <- 3L
  x <- Variable(dims)
  t <- Variable(dims)
  p1 <- Problem(Minimize(sum_entries(t)), list(-t <= x, x <= t))
  v1 <- psolve(p1, solver = "CLARABEL")

  g  <- partial_optimize(p1, list(t), list(x))
  p2 <- Problem(Minimize(g))
  v2 <- psolve(p2, solver = "CLARABEL")
  expect_equal(v1, v2, tolerance = 1e-3)
})

## @cvxpy test_atoms.py::TestAtoms::test_partial_optimize_simple_problem
test_that("partial_optimize: two-stage LP decomposition", {
  skip_if_not_installed("clarabel")
  x <- Variable(1)
  y <- Variable(1)

  ## Combined one-shot.
  p1 <- Problem(Minimize(x + y), list(x + y >= 3, y >= 4, x >= 5))
  v1 <- psolve(p1, solver = "CLARABEL")

  ## Decomposed via partial_optimize.
  p2 <- Problem(Minimize(y), list(x + y >= 3, y >= 4))
  g  <- partial_optimize(p2, list(y), list(x))
  p3 <- Problem(Minimize(x + g), list(x >= 5))
  v3 <- psolve(p3, solver = "CLARABEL")
  expect_equal(v1, v3, tolerance = 1e-3)
})
