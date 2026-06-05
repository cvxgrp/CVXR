## CVXPY 1.9 NLP variable-attribute initialization
## (nlp_tests/test_initialization.py::TestVariableAttributeInit).
##
## Diagonal variables in the NLP path: CvxAttr2Constr lowers a diag (n,n)
## variable to an n-vector, and the original variable's value must propagate to
## that lowered vector so the NLP initial point survives (fixed in
## cvx_attr2constr.R, mirroring CVXPY cvx_attr2constr.py:251-252).

## @cvxpy nlp_tests/test_initialization.py::TestVariableAttributeInit::test_simple_diag_variable
test_that("init: diagonal variable solves without crashing", {
  if (!.de_have_backend()) skip("No R NLP backend")
  skip_if_not_installed("ipopt")
  n <- 3
  D <- Variable(c(n, n), diag = TRUE)
  prob <- Problem(Maximize(sum_entries(log(exp(D)))))
  ## CVXPY only checks that this does not crash (max_iter caps the unbounded max).
  .de_ipopt_solve(prob, max_iter = 10)
  expect_true(status(prob) %in% c("optimal", "optimal_inaccurate", "user_limit",
                                  "unbounded", "solver_error"))
})

## @cvxpy nlp_tests/test_initialization.py::TestVariableAttributeInit::test_diag_variable_value_sparse_init
test_that("init: diagonal variable with sparse value initialization", {
  if (!.de_have_backend()) skip("No R NLP backend")
  skip_if_not_installed("ipopt")
  skip_if_not_installed("Matrix")
  n <- 3
  D <- Variable(c(n, n), diag = TRUE)
  value(D) <- Matrix::Diagonal(x = c(1, 2, 3))   # sparse diagonal initialization
  prob <- Problem(Minimize(sum_squares(D)), list(sum_entries(D) >= 1))
  .de_ipopt_solve(prob)
  expect_equal(status(prob), "optimal")
  expect_equal(sum(as.matrix(value(D))), 1.0, tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_initialization.py::TestVariableAttributeInit::test_advanced_pricing_problem
test_that("init: nonconvex advanced pricing problem with diag variable", {
  if (!.de_have_backend()) skip("No R NLP backend")
  skip_if_not_installed("ipopt")
  set.seed(42); n <- 5; N <- 10; rank <- 2
  D <- matrix(rnorm(n * N), n, N)
  Pitilde <- matrix(rnorm((n + 1) * N), n + 1, N)
  Etilde <- Variable(c(n, n + 1L))
  Ediag <- Variable(c(n, n), diag = TRUE)
  B <- Variable(c(n, rank)); Cmat <- Variable(c(rank, n))
  prob <- Problem(
    Maximize(sum_entries(D * (Etilde %*% Pitilde) - exp(Etilde %*% Pitilde))),
    list(Etilde[, 1:n] == Ediag + B %*% Cmat,
         Ediag <= 0,
         Etilde >= -10, Etilde <= 10,
         B >= -5, B <= 5, Cmat >= -5, Cmat <= 5))
  expect_false(is_dcp(prob))
  value(Etilde) <- matrix(rnorm(n * (n + 1)), n, n + 1)
  value(Ediag) <- diag(-abs(rnorm(n)))
  value(B) <- matrix(rnorm(n * rank), n, rank)
  value(Cmat) <- matrix(rnorm(rank * n), rank, n)
  .de_ipopt_solve(prob, max_iter = 200)
  expect_equal(status(prob), "optimal")
})
