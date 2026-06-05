# Tests for best_of NLP random restarts + sample_bounds (CVXPY 1.9.0).
#
# psolve(prob, nlp = TRUE, best_of = n) runs the solve from n random initial
# points and keeps the best (lowest-objective) result; the per-run objectives
# are exposed via solver_stats(prob)@extra_stats$all_objs_from_best_of. Random
# initialization draws from each variable's sample_bounds (if set) or its finite
# variable bounds.
#
# CVXPY SOURCE: reductions/solvers/nlp_solving_chain.py:69-238,
#   tests/nlp_tests/test_best_of.py

library(CVXR)

## @cvxpy nlp_tests/test_best_of.py::TestBestOf::test_path_planning_best_of_three
test_that("best_of with finite variable bounds runs N times and keeps the best", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  set.seed(1)
  x <- Variable(bounds = list(-5, 5))
  y <- Variable(bounds = list(-3, 3))
  prob <- Problem(Minimize((x - 1)^2 + (y - 2)^2))
  val <- psolve(prob, solver = "UNO", nlp = TRUE, best_of = 3)
  expect_equal(status(prob), "optimal")
  expect_equal(val, 0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 1, tolerance = 1e-3)
  expect_equal(as.numeric(value(y)), 2, tolerance = 1e-3)

  all_objs <- solver_stats(prob)@extra_stats$all_objs_from_best_of
  expect_equal(length(all_objs), 3L)
  expect_equal(min(all_objs), val, tolerance = 1e-6)
})

## @cvxpy nlp_tests/test_best_of.py::TestBestOf::test_circle_packing_best_of_one
test_that("sample_bounds enables best_of on an otherwise-unbounded variable", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  set.seed(2)
  ## x has NO finite variable bounds, but sample_bounds supplies the sampling
  ## region, so best_of can randomize it without error.
  x <- Variable(2)
  sample_bounds(x) <- c(-5, 5)
  expect_equal(sample_bounds(x), list(-5, 5))   # stored as (low, high); broadcast at sampling
  prob <- Problem(Minimize(sum_squares(x - c(1, 2))))
  val <- psolve(prob, solver = "UNO", nlp = TRUE, best_of = 4)
  expect_equal(status(prob), "optimal")
  expect_equal(val, 0, tolerance = 1e-4)
  expect_equal(length(solver_stats(prob)@extra_stats$all_objs_from_best_of), 4L)
})

## @cvxpy nlp_tests/test_best_of.py::TestBestOf::test_path_planning_best_of_four
test_that("best_of errors on an unbounded variable with no sample_bounds or value", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  x <- Variable(bounds = list(-5, 5))
  y <- Variable(bounds = list(-3, NULL))   # one infinite bound, no sample_bounds
  prob <- Problem(Minimize((x - 1)^2 + (y - 2)^2))
  expect_error(psolve(prob, solver = "UNO", nlp = TRUE, best_of = 3), "non-finite sampling bounds")
})

## @cvxpy nlp_tests/test_best_of.py::TestBestOf::test_path_planning_best_of_five
test_that("best_of does not randomize a variable that has an assigned value", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  set.seed(3)
  x <- Variable(bounds = list(-5, 5))
  y <- Variable()           # unbounded, but...
  value(y) <- 5             # ...assigned a value -> restored each run, not randomized
  prob <- Problem(Minimize((x - 1)^2 + (y - 2)^2))
  val <- psolve(prob, solver = "UNO", nlp = TRUE, best_of = 3)
  expect_equal(status(prob), "optimal")
  expect_equal(val, 0, tolerance = 1e-4)
  expect_equal(length(solver_stats(prob)@extra_stats$all_objs_from_best_of), 3L)
})

## @cvxpy reductions/solvers/nlp_solving_chain.py
test_that("best_of must be a positive integer", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("Uno")
  x <- Variable(bounds = list(-5, 5))
  prob <- Problem(Minimize((x - 1)^2))
  expect_error(psolve(prob, solver = "UNO", nlp = TRUE, best_of = 0), "positive integer")
  expect_error(psolve(prob, solver = "UNO", nlp = TRUE, best_of = -2), "positive integer")
})

## The remaining three CVXPY best_of tests use IPOPT (the default NLP backend):
## they exercise random-restart distinctness, infeasible -> +Inf, and unbounded
## -> -Inf, all of which depend on robust IPOPT status detection.

## @cvxpy nlp_tests/test_best_of.py::TestBestOf::test_path_planning_best_of_two
test_that("best_of randomizes within sample_bounds even when value is set", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("ipopt")
  set.seed(0); n <- 5
  radius <- runif(n, 1.0, 3.0)
  centers <- Variable(c(n, 2L), name = "c")
  cons <- list()
  for (i in seq_len(n - 1L)) {
    cons <- c(cons, list(
      sum_entries(square(centers[i, ] - centers[(i + 1L):n, ]), axis = 1) >=
        (radius[i] + radius[(i + 1L):n])^2))
  }
  obj <- Minimize(max_entries(cvxr_norm(centers, "inf", axis = 1) + matrix(radius, n, 1)))
  prob <- Problem(obj, cons)
  value(centers) <- matrix(runif(2 * n), n, 2)   # value set, but...
  sample_bounds(centers) <- c(-5, 5)             # ...best_of still randomizes within these
  n_runs <- 10
  val <- psolve(prob, solver = "IPOPT", nlp = TRUE, best_of = n_runs, print_level = 0L)
  all_objs <- solver_stats(prob)@extra_stats$all_objs_from_best_of
  expect_equal(max(table(all_objs)), 1L)         # every restart distinct
  expect_equal(length(all_objs), n_runs)
  cc <- matrix(as.numeric(value(centers)), n, 2)
  manual_obj <- max(apply(abs(cc), 1, max) + radius)
  expect_equal(manual_obj, val, tolerance = 1e-6)
  expect_equal(manual_obj, min(all_objs), tolerance = 1e-6)
})

## @cvxpy nlp_tests/test_best_of.py::TestBestOf::test_best_of_infeasible_problem
test_that("best_of returns +Inf objective for an infeasible problem", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("ipopt")
  x <- Variable(bounds = list(-5, 5))
  y <- Variable(bounds = list(-3, 3))
  prob <- Problem(Minimize((x - 1)^2 + (y - 2)^2), list(x + y == 10))  # max x+y = 8
  val <- psolve(prob, solver = "IPOPT", nlp = TRUE, best_of = 20, print_level = 0L)
  expect_equal(val, Inf)
})

## @cvxpy nlp_tests/test_best_of.py::TestBestOf::test_best_of_with_unbounded
test_that("best_of returns -Inf objective for an unbounded problem", {
  skip_if_not_installed("sparsediff")
  skip_if_not_installed("ipopt")
  x <- Variable()
  sample_bounds(x) <- c(-5, 5)
  prob <- Problem(Minimize(x))
  val <- psolve(prob, solver = "IPOPT", nlp = TRUE, best_of = 20, print_level = 0L)
  expect_equal(val, -Inf)
})
