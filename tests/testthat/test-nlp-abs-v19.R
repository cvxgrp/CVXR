## CVXPY 1.9 NLP lasso (abs) stress tests (nlp_tests/test_abs.py::TestAbs).
##
## Lasso sum_squares(A x - b) + lambda * sum(|x|) solved via the DNLP path (abs
## is PWL -> epigraph canonicalization) must reach the same objective as the DCP
## (CLARABEL) solve. DNLP_opt == DCP_opt is a data-independent invariant, so
## R-generated A/b verify the same property CVXPY checks. The factor sweep is
## reduced by default and runs CVXPY's full np.linspace(0.1, 1, 20) under
## CVXR_NLP_FULL=1. Only the small case runs the DerivativeChecker mirror.

.lasso_factors <- function() if (.nlp_full()) seq(0.1, 1, length.out = 20) else c(0.1, 0.55, 1.0)

.lasso_compare <- function(m, n, check_derivative = FALSE) {
  if (!.de_have_backend()) skip("No R NLP backend")
  skip_if_not_installed("ipopt")
  set.seed(0)
  for (factor in .lasso_factors()) {
    b <- rnorm(m)
    A <- matrix(rnorm(m * n), m, n)
    lmbda <- factor * 2 * max(abs(t(A) %*% b))

    x <- Variable(n, name = "x")
    obj <- sum_squares(A %*% x - b) + lmbda * sum_entries(abs(x))
    prob_dcp <- Problem(Minimize(obj))
    psolve(prob_dcp, solver = "CLARABEL")
    obj_dcp <- value(prob_dcp)

    x <- Variable(n, name = "x")
    obj <- sum_squares(A %*% x - b) + lmbda * sum_entries(abs(x))
    prob_nlp <- Problem(Minimize(obj))
    .de_ipopt_solve(prob_nlp, hessian_approximation = "exact")
    obj_nlp <- value(prob_nlp)

    expect_lt(abs(obj_nlp - obj_dcp) / abs(obj_nlp), 1e-4)
    if (check_derivative) nlp_derivative_check(prob_nlp)
  }
}

## @cvxpy nlp_tests/test_abs.py::TestAbs::test_lasso_square_small
test_that("lasso: square small (m=n=10)", {
  .lasso_compare(10, 10, check_derivative = TRUE)
})

## @cvxpy nlp_tests/test_abs.py::TestAbs::test_lasso_square
test_that("lasso: square (m=n=50)", {
  .lasso_compare(50, 50)
})

## @cvxpy nlp_tests/test_abs.py::TestAbs::test_lasso_underdetermined
test_that("lasso: underdetermined (m=100, n=200)", {
  .lasso_compare(100, 200)
})

## @cvxpy nlp_tests/test_abs.py::TestAbs::test_lasso_overdetermined
test_that("lasso: overdetermined (m=200, n=100)", {
  .lasso_compare(200, 100)
})
