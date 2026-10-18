## problem_data() shares the compile cache with psolve().
##
## CVXPY SOURCE: problems/problem.py get_problem_data (1.9.3 lines 849-897).
## CVXPY has one compile path: on a cache hit it calls every reduction's
## update_parameters() and the solver's apply() on the cached param_prog; on a
## miss it applies the chain and caches param_prog. solve() calls it, so a
## user's get_problem_data() between two solves is harmless.
##
## Through 1.9.2.9513 CVXR's problem_data() always re-applied the full chain.
## Dgp2Dcp's apply builds fresh log-space Parameters each time, so after a
## problem_data() call the param_prog cached by the first psolve() pointed at
## orphaned Parameters and the next solve after a parameter edit used the
## first solve's values (status "optimal", value 4.999994 instead of
## 2 sqrt(6) = 4.898979). Found by the lifecycle fuzz
## (notes/lifecycle_fuzz_2026-10-17.md, finding 6.1).

gp_problem <- function() {
  a <- Parameter(1, pos = TRUE, value = 1)
  b <- Parameter(1, pos = TRUE, value = 1)
  x <- Variable(1, pos = TRUE)
  prob <- Problem(Minimize(a * x + b / x), list(x >= 0.1, x <= 10))
  list(a = a, b = b, x = x, prob = prob)
}

## @cvxpy problems/problem.py::Problem.get_problem_data test_dgp_dpp.py::TestDgpDpp::test_fast_path_with_changing_params
test_that("DGP: problem_data() between a parameter edit and a solve does not stale the solve", {
  require_solver("CLARABEL")
  g <- gp_problem()
  v1 <- psolve(g$prob, solver = "CLARABEL", gp = TRUE)
  expect_equal(v1, 2, tolerance = 1e-5)
  value(g$a) <- 2
  value(g$b) <- 3
  pd <- problem_data(g$prob, solver = "CLARABEL", gp = TRUE)
  expect_type(pd, "list")
  v2 <- psolve(g$prob, solver = "CLARABEL", gp = TRUE)
  expect_equal(status(g$prob), "optimal")
  expect_equal(v2, 2 * sqrt(6), tolerance = 1e-5)
  expect_equal(as.numeric(value(g$x)), sqrt(3 / 2), tolerance = 1e-4)
})

## @cvxpy problems/problem.py::Problem.get_problem_data test_dgp_dpp.py::TestDgpDpp::test_fast_path_with_changing_params
test_that("DGP: problem_data() before a parameter edit does not stale the next solve", {
  require_solver("CLARABEL")
  g <- gp_problem()
  psolve(g$prob, solver = "CLARABEL", gp = TRUE)
  invisible(problem_data(g$prob, solver = "CLARABEL", gp = TRUE))
  value(g$a) <- 2
  value(g$b) <- 3
  v2 <- psolve(g$prob, solver = "CLARABEL", gp = TRUE)
  expect_equal(v2, 2 * sqrt(6), tolerance = 1e-5)
  expect_equal(as.numeric(value(g$x)), sqrt(3 / 2), tolerance = 1e-4)
})

## @cvxpy problems/problem.py::Problem.get_problem_data
test_that("DGP: the deprecated get_problem_data() takes the same cached path", {
  require_solver("CLARABEL")
  g <- gp_problem()
  psolve(g$prob, solver = "CLARABEL", gp = TRUE)
  value(g$a) <- 2
  value(g$b) <- 3
  suppressWarnings(invisible(get_problem_data(g$prob, solver = "CLARABEL", gp = TRUE)))
  v2 <- psolve(g$prob, solver = "CLARABEL", gp = TRUE)
  expect_equal(v2, 2 * sqrt(6), tolerance = 1e-5)
})

## @cvxpy problems/problem.py::Problem.get_problem_data
test_that("DGP: problem_data() called first caches, and a later edit is reflected in both data and solve", {
  require_solver("CLARABEL")
  g <- gp_problem()
  pd1 <- problem_data(g$prob, solver = "CLARABEL", gp = TRUE)
  ## The first call populates the cache; the second reuses it (CVXPY's fast
  ## path), and the data must follow the new parameter values.
  value(g$a) <- 2
  value(g$b) <- 3
  pd2 <- problem_data(g$prob, solver = "CLARABEL", gp = TRUE)
  expect_false(isTRUE(all.equal(pd1$data[["b"]], pd2$data[["b"]])) &&
               isTRUE(all.equal(pd1$data[["c"]], pd2$data[["c"]])))
  expect_equal(length(pd1$inverse_data), length(pd2$inverse_data))
  v <- psolve(g$prob, solver = "CLARABEL", gp = TRUE)
  expect_equal(v, 2 * sqrt(6), tolerance = 1e-5)
})

## @cvxpy problems/problem.py::Problem.get_problem_data
test_that("DCP: problem_data() between an edit and a solve stays correct (control)", {
  require_solver("CLARABEL")
  p <- Parameter(2, value = c(1, 1))
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - p)), list(x >= -10))
  psolve(prob, solver = "CLARABEL")
  value(p) <- c(3, 4)
  invisible(problem_data(prob, solver = "CLARABEL"))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), c(3, 4), tolerance = 1e-4)
})

## @cvxpy problems/problem.py::Problem.get_problem_data
test_that("problem_data() after psolve() returns the same data as a fresh problem", {
  require_solver("CLARABEL")
  A <- Parameter(c(2, 3), value = matrix(1:6, 2, 3))
  b <- Parameter(2, value = c(1, 2))
  x <- Variable(3, nonneg = TRUE)
  prob <- Problem(Minimize(sum(x)), list(A %*% x >= b))
  psolve(prob, solver = "CLARABEL")
  value(A) <- matrix(c(2, 0, 1, 3, 0, 5), 2, 3)
  value(b) <- c(4, 3)
  pd <- problem_data(prob, solver = "CLARABEL")
  A2 <- Parameter(c(2, 3), value = matrix(c(2, 0, 1, 3, 0, 5), 2, 3))
  b2 <- Parameter(2, value = c(4, 3))
  x2 <- Variable(3, nonneg = TRUE)
  fresh <- problem_data(Problem(Minimize(sum(x2)), list(A2 %*% x2 >= b2)),
                        solver = "CLARABEL")
  expect_equal(as.numeric(pd$data[["c"]]), as.numeric(fresh$data[["c"]]))
  expect_equal(as.numeric(pd$data[["b"]]), as.numeric(fresh$data[["b"]]))
  expect_equal(as.matrix(pd$data[["A"]]), as.matrix(fresh$data[["A"]]))
})
