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

## ------------------------------------------------------------------
## PartialProblem.grad (partial_optimize.py:225-277), claims-audit
## decision 13, ledger rows C0969/C0968 (CVXR 1.9.2.9148). Before 9148
## the Lagrangian code was registered on the per-atom `.grad` hook, so
## `grad(partial_problem)` aborted "grad() not implemented"; the body
## also keyed the result off the hook's per-argument list, used the
## constraint's first argument instead of lhs - rhs, and gated on
## optimal/optimal_inaccurate only. Expected values measured on CVXPY
## 1.9.3 + Clarabel (scratch po_d13_cvxpy.py, 2026-10-10).
## ------------------------------------------------------------------

## @cvxpy test_grad.py::TestSpecialCases::test_partial_problem
test_that("partial_optimize: grad of a partial minimization/maximization", {
  skip_if_not_installed("clarabel")
  a <- Variable(name = "a")
  x <- Variable(2, name = "x")
  ka <- as.character(a@id)
  kx <- as.character(x@id)
  as_vec <- function(m) as.numeric(as.matrix(m))

  for (obj in list(Minimize(a^-1), Maximize(entr(a)))) {
    prob <- Problem(obj, list(x + a >= c(5, 8)))

    ## Optimize over nothing: no values -> NULL.
    expr <- partial_optimize(prob, dont_opt_vars = list(x, a), solver = "CLARABEL")
    value(a) <- NULL
    value(x) <- NULL
    g <- grad(expr)
    expect_null(g[[ka]])
    expect_null(g[[kx]])

    ## Outside the domain (inner problem infeasible) -> NULL.
    value(a) <- 1
    value(x) <- c(5, 5)
    g <- grad(expr)
    expect_null(g[[ka]])
    expect_null(g[[kx]])

    ## Inside: the gradient is the objective's; x does not enter.
    value(a) <- 1
    value(x) <- c(10, 10)
    g <- grad(expr)
    expect_equal(as_vec(g[[ka]]), as_vec(grad(obj@args[[1L]])[[ka]]),
                 tolerance = 1e-6)
    expect_equal(as_vec(g[[kx]]), c(0, 0), tolerance = 1e-6)

    ## Optimize over x: constraint slack, dual zero.
    expr <- partial_optimize(prob, opt_vars = list(x), solver = "CLARABEL")
    value(a) <- 1
    g <- grad(expr)
    expect_equal(as_vec(g[[ka]]), as_vec(grad(obj@args[[1L]])[[ka]]),
                 tolerance = 1e-6)

    ## Optimize over a: d/dx equals the dual of x + a >= c(5, 8) at x = 0
    ## (C0969: the Lagrangian must use lhs - rhs, which carries x).
    fix_prob <- Problem(obj, list(x + a >= c(5, 8), x == 0))
    psolve(fix_prob, solver = "CLARABEL")
    dual_val <- as.numeric(dual_value(fix_prob@constraints[[1L]]))
    expr <- partial_optimize(prob, opt_vars = list(a), solver = "CLARABEL")
    value(x) <- c(0, 0)
    g <- grad(expr)
    expect_equal(as_vec(g[[kx]]), dual_val, tolerance = 1e-5)

    ## Optimize over everything: nothing left to differentiate.
    expr <- partial_optimize(prob, opt_vars = list(x, a), solver = "CLARABEL")
    expect_length(grad(expr), 0L)
  }
})

## @cvxpy NONE -- C0968: the gate is SOLUTION_PRESENT, which includes user_limit
test_that("partial_optimize: grad is returned when the inner solve hits a user limit", {
  skip_if_not_installed("clarabel")
  a <- Variable(name = "a")
  x <- Variable(2, name = "x")
  prob <- Problem(Maximize(entr(a)), list(x + a >= c(5, 8)))

  ## Clarabel stops after one iteration -> MaxIterations -> user_limit,
  ## with a (rough) primal/dual point present.
  inner <- Problem(Maximize(entr(a)), list(x + a >= c(5, 8), x == c(0, 0)))
  suppressWarnings(psolve(inner, solver = "CLARABEL", max_iter = 1L))
  expect_equal(status(inner), "user_limit")

  expr <- partial_optimize(prob, opt_vars = list(a), solver = "CLARABEL",
                           max_iter = 1L)
  value(x) <- c(0, 0)
  g <- suppressWarnings(grad(expr))
  gx <- g[[as.character(x@id)]]
  expect_false(is.null(gx))
  expect_length(as.numeric(as.matrix(gx)), 2L)
  expect_true(all(is.finite(as.numeric(as.matrix(gx)))))
})

## @cvxpy NONE -- partial_optimize.py:197-203: bounds (-Inf, Inf); a bounded
## solver reads them through get_expr_bounds_if_supported. CVXPY 1.9.3 solves
## this on HIGHS: optimal 1.5, x = (1.5, 0).
test_that("partial_optimize: get_bounds is unbounded and a bounded solver reads it", {
  skip_if_not_installed("highs")
  x <- Variable(2)
  y <- Variable(2)
  pp <- partial_optimize(Problem(Minimize(sum(y)), list(y >= x)), opt_vars = list(y))
  b <- get_bounds(pp)
  expect_equal(as.numeric(b[[1L]]), -Inf)
  expect_equal(as.numeric(b[[2L]]), Inf)
  prob <- Problem(Minimize(max_elemwise(pp, 3 - x[1])), list(x >= 0, x <= 5))
  expect_equal(psolve(prob, solver = "HIGHS"), 1.5, tolerance = 1e-6)
  expect_equal(status(prob), "optimal")
})

## @cvxpy NONE -- problem_form.py:207-228 predicts cones from atoms(), which
## reaches the inner problem through PartialProblem.args = [prob]
## (partial_optimize.py:121). CVXPY 1.9.3 predicts {SOC, Zero} / {ExpCone,
## Zero} here and raises SolverError on HIGHS; CVXR 1.9.2.9201 predicted no
## cones, so HIGHS returned Inf "infeasible" (SOC) or crashed (ExpCone).
test_that("partial_optimize: cone prediction sees the inner problem's atoms", {
  x <- Variable(name = "x")
  y <- Variable(name = "y")
  outer <- function(inner_obj) {
    pp <- partial_optimize(Problem(Minimize(inner_obj), list(y >= 1)),
                           opt_vars = list(y), dont_opt_vars = list(x))
    Problem(Minimize(pp), list(x == 0))
  }
  cone_names <- function(prob)
    sort(vapply(CVXR:::.required_cone_types(prob), function(c) c@name, ""))
  expect_equal(cone_names(outer(y - x)), character(0))
  expect_equal(cone_names(outer(square(y - x))), "SOC")
  expect_equal(cone_names(outer(exp(y - x))), "ExpCone")

  skip_if_not_installed("highs")
  expect_error(psolve(outer(square(y - x)), solver = "HIGHS"),
               "does not support the required cone types")
  expect_error(psolve(outer(exp(y - x)), solver = "HIGHS"),
               "does not support the required cone types")
  expect_equal(psolve(outer(y - x), solver = "HIGHS"), 1, tolerance = 1e-6)
  prob <- outer(square(y - x))
  expect_equal(psolve(prob, solver = "CLARABEL"), 1, tolerance = 1e-6)
  expect_equal(status(prob), "optimal")
})
