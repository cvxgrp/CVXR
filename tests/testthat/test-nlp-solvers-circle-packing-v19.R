## CVXPY 1.9 NLP circle-packing examples
## (nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_circle_packing_*).
##
## Three formulations of the same min-enclosing-square circle packing. CVXPY
## solves with IPOPT and asserts the numpy-RNG-specific optimal centers/objective;
## R's RNG generates different radii, so those exact values are not reproducible.
## We instead assert OPTIMAL status, packing feasibility (every pair of centers is
## at least the sum of radii apart), and the faithful DerivativeChecker mirror.
## CVXPY axis=0 (per-column) -> CVXR axis=2; radius is a (1, n) row so it
## broadcasts against the (1, n) per-column reductions.

.cp_pack_constraints <- function(centers, radius, n) {
  con <- list()
  for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
    con <- c(con, list(sum_entries(square(centers[, i] - centers[, j])) >=
                         (radius[i] + radius[j])^2))
  }
  con
}

.cp_assert_feasible <- function(centers, radius, n, tol = 1e-6) {
  C <- matrix(as.numeric(value(centers)), 2, n)
  for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
    dist_sq <- sum((C[, i] - C[, j])^2)
    expect_gte(dist_sq - (radius[i] + radius[j])^2, -tol)
  }
}

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_circle_packing_formulation_one
test_that("NLP example: circle packing (epigraph formulation)", {
  set.seed(5); n <- 3
  radius <- runif(n, 1.0, 3.0)
  centers <- Variable(c(2L, n), name = "c")
  value(centers) <- matrix(runif(2 * n, -5.0, 5.0), 2, n)
  t <- Variable(); value(t) <- 10
  con <- c(.cp_pack_constraints(centers, radius, n),
           list(max_entries(cvxr_norm(centers, "inf", axis = 2) + matrix(radius, 1, n)) <= t))
  prob <- Problem(Minimize(t), con)
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  .cp_assert_feasible(centers, radius, n)
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_circle_packing_formulation_two
test_that("NLP example: circle packing (norm_inf objective)", {
  set.seed(5); n <- 3
  radius <- runif(n, 1.0, 3.0)
  centers <- Variable(c(2L, n), name = "c")
  value(centers) <- matrix(runif(2 * n, -5.0, 5.0), 2, n)
  con <- .cp_pack_constraints(centers, radius, n)
  obj <- Minimize(max_entries(cvxr_norm(centers, "inf", axis = 2) + matrix(radius, 1, n)))
  prob <- Problem(obj, con)
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  .cp_assert_feasible(centers, radius, n)
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_nlp_solvers.py::TestNLPExamples::test_circle_packing_formulation_three
test_that("NLP example: circle packing (max max abs objective)", {
  set.seed(5); n <- 3
  radius <- runif(n, 1.0, 3.0)
  centers <- Variable(c(2L, n), name = "c")
  value(centers) <- matrix(runif(2 * n, -5.0, 5.0), 2, n)
  con <- .cp_pack_constraints(centers, radius, n)
  obj <- Minimize(max_entries(max_entries(abs(centers), axis = 2) + matrix(radius, 1, n)))
  prob <- Problem(obj, con)
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  .cp_assert_feasible(centers, radius, n)
  nlp_derivative_check(prob)
})
