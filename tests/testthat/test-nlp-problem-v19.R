## CVXPY 1.9 NLP initial-point logic (nlp_tests/test_problem.py::TestProblem).
##
## Tests the internal `.set_nlp_initial_point` (= CVXPY `_set_nlp_initial_point`):
## for each variable without a user value, build an initial point from
## get_bounds() -- midpoint when both bounds are finite, one unit inside a single
## finite bound, else 0; sign attributes (nonneg) fold into the bounds. CVXPY's
## `None` bound is R's `Inf`/`-Inf`.

.set_init <- get(".set_nlp_initial_point", envir = asNamespace("CVXR"))

## @cvxpy nlp_tests/test_problem.py::TestProblem::test_set_initial_point_scalar_bounds
test_that("NLP init point: scalar bounds variants", {
  cases <- list(
    list(bounds = NULL,              exp = 0.0),
    list(bounds = list(-Inf, Inf),   exp = 0.0),
    list(bounds = list(-Inf, 3.5),   exp = 2.5),
    list(bounds = list(3.5, Inf),    exp = 4.5),
    list(bounds = list(3.5, 4.5),    exp = 4.0))
  for (cs in cases) {
    x <- if (is.null(cs$bounds)) Variable(3) else Variable(3, bounds = cs$bounds)
    prob <- Problem(Minimize(sum_entries(x)))
    .set_init(prob)
    expect_equal(as.numeric(value(x)), rep(cs$exp, 3))
  }
})

## @cvxpy nlp_tests/test_problem.py::TestProblem::test_set_initial_point_mixed_inf_and_finite
test_that("NLP init point: per-entry mixed inf/finite bounds", {
  lb <- c(-Inf, 3.5, -Inf, -1.5, 2, 2.5)
  ub <- c(-4, 4.5, Inf, 4.5, Inf, 4.5)
  x <- Variable(6, bounds = list(lb, ub))
  prob <- Problem(Minimize(sum_entries(x)))
  .set_init(prob)
  expect_equal(as.numeric(value(x)), c(-5, 4.0, 0.0, 1.5, 3, 3.5))
})

## @cvxpy nlp_tests/test_problem.py::TestProblem::test_set_initial_point_two_variables
test_that("NLP init point: two variables", {
  x <- Variable(2, bounds = list(-Inf, Inf))
  y <- Variable(2, bounds = list(-3, Inf))
  prob <- Problem(Minimize(sum_entries(x) + sum_entries(y)))
  .set_init(prob)
  expect_equal(as.numeric(value(x)), rep(0, 2))
  expect_equal(as.numeric(value(y)), rep(-2, 2))
})

## @cvxpy nlp_tests/test_problem.py::TestProblem::test_set_initial_point_nonnegative_attributes
test_that("NLP init point: nonneg attribute -> 1", {
  x <- Variable(2, nonneg = TRUE)
  prob <- Problem(Minimize(sum_entries(x)))
  .set_init(prob)
  expect_equal(as.numeric(value(x)), rep(1, 2))
})

## @cvxpy nlp_tests/test_problem.py::TestProblem::test_set_initial_point_nonnegative_attributes_and_bounds
test_that("NLP init point: nonneg attribute + bounds[1, None] -> 2", {
  x <- Variable(2, nonneg = TRUE, bounds = list(1, Inf))
  prob <- Problem(Minimize(sum_entries(x)))
  .set_init(prob)
  expect_equal(as.numeric(value(x)), rep(2, 2))
})
