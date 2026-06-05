# Tests for the DNLP predicate cluster (CVXPY 1.9.0).
#
# is_dnlp / is_smooth / is_linearizable_convex / is_linearizable_concave /
# is_atom_smooth are new in CVXPY 1.9.0. The classically-differentiable atoms
# (affine, exp, log, entr, logistic, kl_div, rel_entr, xexp, power with a
# constant exponent, geo_mean, log_sum_exp, quad_form, quad_over_lin, prod,
# trig, hyperbolic, normcdf) report is_atom_smooth = TRUE, so a smooth convex
# function such as sum_squares is linearizable in BOTH directions. Piecewise
# atoms (abs, min, max) stay convex/concave-only. These tests pin the machinery
# and the curvature classes. Expected values cross-checked against the CVXPY
# 1.9.0 composition rules (expression.py:349-405, atom.py:191-298,
# leaf.py:267-273).
#
library(CVXR)

## @cvxpy expressions/leaf.py::Leaf.is_linearizable_convex
test_that("a leaf is smooth and DNLP", {
  x <- Variable(2)
  expect_true(is_dnlp(x))
  expect_true(is_smooth(x))
  expect_true(is_linearizable_convex(x))
  expect_true(is_linearizable_concave(x))
})

## @cvxpy atoms/atom.py::Atom.is_linearizable_convex
test_that("an affine expression is smooth and DNLP", {
  x <- Variable(2)
  e <- 2 * x + 1
  expect_true(is_smooth(e))
  expect_true(is_dnlp(e))
  expect_true(is_linearizable_convex(e))
  expect_true(is_linearizable_concave(e))
})

## @cvxpy atoms/atom.py::Atom.is_linearizable_convex
test_that("a convex-only atom is linearizable-convex and DNLP but not smooth", {
  x <- Variable(2)
  ## max is a piecewise-linear convex atom: not is_atom_smooth, so it is
  ## linearizable-convex only (unlike sum_squares, which is smooth).
  cvx <- max_entries(x)
  expect_true(is_linearizable_convex(cvx))
  expect_false(is_linearizable_concave(cvx))
  expect_true(is_dnlp(cvx))
  ## convex-only -> not smooth (smooth needs both directions)
  expect_false(is_smooth(cvx))
})

## @cvxpy atoms/atom.py::Atom.is_linearizable_concave
test_that("a concave-only atom is linearizable-concave and DNLP but not smooth", {
  x <- Variable(2)
  ## min is a piecewise-linear concave atom: not is_atom_smooth, so it is
  ## linearizable-concave only (unlike sqrt = power(., 0.5), which is smooth).
  s <- min_entries(x)
  expect_true(is_linearizable_concave(s))
  expect_false(is_linearizable_convex(s))
  expect_true(is_dnlp(s))
  expect_false(is_smooth(s))
})

## @cvxpy expressions/expression.py::Expression.is_dnlp
test_that("a non-DCP, non-smooth expression is not DNLP", {
  x <- Variable(2); y <- Variable(2)
  ## Difference of two piecewise-linear convex atoms: neither DCP nor smooth.
  nondcp <- sum_entries(abs(x)) - sum_entries(abs(y))
  expect_false(is_dcp(nondcp))
  expect_false(is_dnlp(nondcp))
})

## @cvxpy atoms/atom.py::Atom.is_atom_smooth
test_that("is_atom_smooth is TRUE for smooth atoms, FALSE for piecewise atoms", {
  x <- Variable(2)
  expect_true(is_atom_smooth(sum_squares(x)))
  expect_false(is_atom_smooth(abs(x)))
})

## @cvxpy problems/objective.py::Minimize.is_dnlp
test_that("Minimize objective DNLP delegates to linearizable-convex", {
  x <- Variable(2)
  expect_true(is_dnlp(Minimize(sum_squares(x))))
  ## concave-only objective under Minimize is not linearizable-convex
  expect_false(is_dnlp(Minimize(min_entries(x))))
})

## @cvxpy problems/objective.py::Maximize.is_dnlp
test_that("Maximize objective DNLP delegates to linearizable-concave", {
  x <- Variable(2)
  expect_true(is_dnlp(Maximize(-sum_squares(x))))
  ## convex-only objective under Maximize is not linearizable-concave
  expect_false(is_dnlp(Maximize(max_entries(x))))
})

## @cvxpy constraints/nonpos.py::Inequality.is_dnlp
test_that("constraint DNLP predicates mirror the curvature rules", {
  x <- Variable(2)
  ## Inequality: (lhs - rhs) linearizable-convex
  expect_true(is_dnlp(sum_squares(x) <= 1))
  ## Equality / Zero: argument smooth representable (affine here)
  expect_true(is_dnlp(sum_entries(x) == 0))
  ## A convex-but-not-smooth (piecewise-linear) equality argument is not DNLP
  expect_false(is_dnlp(max_entries(x) == 1))
})

## @cvxpy problems/problem.py::Problem.is_dnlp
test_that("a DCP problem is DNLP; a non-DCP one is not", {
  x <- Variable(2)
  p_dcp <- Problem(Minimize(sum_squares(x)),
                   list(x <= 1, sum_entries(x) == 0))
  expect_true(is_dnlp(p_dcp))

  y <- Variable(2)
  p_nondcp <- Problem(Minimize(sum_entries(abs(x)) - sum_entries(abs(y))))
  expect_false(is_dnlp(p_nondcp))
})

## @cvxpy transforms/indicator.py::indicator.is_linearizable_convex
test_that("Indicator is linearizable-convex only", {
  x <- Variable(2)
  ind <- Indicator(list(x >= 0))
  expect_true(is_linearizable_convex(ind))
  expect_false(is_linearizable_concave(ind))
})
