## Wave 5: Power & Approx Atom Tests
## Tests from CVXPY test_power_atom.py gaps
## All expected values verified via `uv run python` against CVXPY 1.8.1

# ═══════════════════════════════════════════════════════════════════════
# Power: approx vs exact agreement
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("power approx: minimize x^3 s.t. x >= 2", {
  ## CVXPY: prob.value = 8.0 for both exact and approx
  x <- Variable()
  prob <- Problem(Minimize(power(x, 3)), list(x >= 2))
  result <- psolve(prob)
  expect_equal(result, 8.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("power approx: maximize x^0.5 s.t. x <= 4", {
  x <- Variable()
  prob <- Problem(Maximize(power(x, 0.5)), list(x <= 4, x >= 0))
  result <- psolve(prob)
  expect_equal(result, 2.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("power approx: minimize x^2 equivalent to sum_squares", {
  x <- Variable(3)
  vals <- c(1, 2, 3)
  prob1 <- Problem(Minimize(sum_entries(power(x, 2))), list(x >= vals))
  prob2 <- Problem(Minimize(sum_squares(x)), list(x >= vals))
  r1 <- psolve(prob1)
  r2 <- psolve(prob2)
  expect_equal(r1, r2, tolerance = 1e-4)
  expect_equal(r1, 14.0, tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════════
# Pnorm: approx vs exact agreement
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("pnorm approx: minimize ||x||_3 s.t. x >= [1,2,3]", {
  ## CVXPY: value = 3.30192725803897
  x <- Variable(3)
  prob <- Problem(Minimize(p_norm(x, 3)), list(x >= c(1, 2, 3)))
  result <- psolve(prob)
  expect_equal(result, 3.30192725803897, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("pnorm approx: minimize ||x||_1.5 s.t. x >= [1,1,1]", {
  x <- Variable(3)
  prob <- Problem(Minimize(p_norm(x, 1.5)), list(x >= c(1, 1, 1)))
  result <- psolve(prob)
  ## ||[1,1,1]||_1.5 = (1^1.5 + 1^1.5 + 1^1.5)^(1/1.5) = 3^(2/3) ≈ 2.0801
  expect_equal(result, 3^(2 / 3), tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════════
# GeoMean: basic solve
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_problem.py::TestProblem::test_geo_mean
test_that("geo_mean: maximize s.t. sum <= 3", {
  ## CVXPY: value = 1.0, x = [1, 1, 1] (AM-GM equality)
  x <- Variable(3, nonneg = TRUE)
  prob <- Problem(Maximize(geo_mean(x)), list(sum_entries(x) <= 3))
  result <- psolve(prob)
  expect_equal(result, 1.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("geo_mean: weighted with equal weights", {
  ## Equal weights = standard geo_mean
  x <- Variable(2, nonneg = TRUE)
  prob <- Problem(Maximize(geo_mean(x, c(1, 1))), list(sum_entries(x) <= 4))
  result <- psolve(prob)
  ## Maximum of sqrt(x1*x2) s.t. x1+x2<=4 → x1=x2=2, value=2
  expect_equal(result, 2.0, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════════
# GeoMean: single non-zero weight = affine
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("single-weight geo_mean: is affine", {
  ## CVXPY: geo_mean(x, [0,0,1]).is_affine() = True
  x <- Variable(3, nonneg = TRUE)
  gm <- geo_mean(x, c(0, 0, 1))
  expect_true(is_affine(gm))
})

## @cvxpy NONE
test_that("single-weight geo_mean: value equals x[k]", {
  ## CVXPY: geo_mean([10,20,30], [0,0,1]) = 30.0
  x <- Variable(3, nonneg = TRUE)
  gm <- geo_mean(x, c(0, 0, 1))
  save_leaf_value(x, matrix(c(10, 20, 30), 3, 1))
  expect_equal(as.numeric(value(gm)), 30.0, tolerance = 1e-10)
})

## @cvxpy NONE
test_that("single-weight geo_mean: solve concentrates in active entry", {
  ## CVXPY: maximize geo_mean(x, [0,0,1]) s.t. sum(x)<=6 → value=6, x=[0,0,6]
  x <- Variable(3, nonneg = TRUE)
  prob <- Problem(Maximize(geo_mean(x, c(0, 0, 1))), list(sum_entries(x) <= 6))
  result <- psolve(prob)
  expect_equal(result, 6.0, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════════
# inv_prod
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("inv_prod: basic solve", {
  ## inv_prod(x) = 1/prod(x)
  ## minimize inv_prod(x) s.t. x <= [1, 2, 3], x >= 0.5 → value = 1/6
  ## (need upper bounds; with only lower bounds, solver pushes x→∞, value→0)
  x <- Variable(3, nonneg = TRUE)
  prob <- Problem(Minimize(inv_prod(x)), list(x <= c(1, 2, 3), x >= 0.5))
  result <- psolve(prob)
  expect_equal(result, 1 / 6, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════════
# Power: edge cases
# ═══════════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power: p=0 gives constant 1", {
  x <- Variable()
  expr <- power(x, 0)
  save_leaf_value(x, matrix(5, 1, 1))
  expect_equal(as.numeric(value(expr)), 1.0)
})

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power: p=1 is affine", {
  x <- Variable()
  expr <- power(x, 1)
  expect_true(is_affine(expr))
})

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power: p=2 is convex", {
  x <- Variable()
  expr <- power(x, 2)
  expect_true(is_convex(expr))
  expect_false(is_concave(expr))
})

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power: p=0.5 is concave", {
  x <- Variable()
  expr <- power(x, 0.5)
  expect_true(is_concave(expr))
  expect_false(is_convex(expr))
})

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power: p=-1 is convex", {
  x <- Variable()
  expr <- power(x, -1)
  expect_true(is_convex(expr))
})

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("power: p=3 is convex (odd integer > 1)", {
  x <- Variable()
  expr <- power(x, 3)
  expect_true(is_convex(expr))
})
