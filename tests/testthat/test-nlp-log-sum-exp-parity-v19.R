## CVXPY 1.9 DNLP log_sum_exp parity (nlp_tests/test_log_sum_exp.py).
##
## log_sum_exp is a smooth convex atom, so it solves through the NLP path
## (canonicalized to exp/sum/log by Dnlp2Smooth). CVXPY runs a Python
## DerivativeChecker, mirrored here by nlp_derivative_check. The random
## DCP-equivalence cases (test_two / test_three) assert the data-independent
## invariant DNLP_opt == DCP_opt (CLARABEL); test_three's larger (300x100) size
## runs only under CVXR_NLP_FULL=1 to keep R CMD check fast.

.lse_nlp_solve <- function(prob) {
  skip_if_not_installed("sparsediff")
  if (!requireNamespace("ipopt", quietly = TRUE) &&
      !requireNamespace("Uno", quietly = TRUE)) {
    skip("No R NLP backend (ipopt or Uno) installed")
  }
  output <- capture.output(val <- psolve(prob, nlp = TRUE), type = "output")
  expect_equal(status(prob), "optimal")
  invisible(val)
}

## @cvxpy nlp_tests/test_log_sum_exp.py::TestLogSumExp::test_one
test_that("DNLP log_sum_exp: minimize log_sum_exp(x) with x >= 1", {
  x <- Variable(3)
  value(x) <- rep(1, 3)
  prob <- Problem(Minimize(log_sum_exp(x)), list(x >= 1))
  val <- .lse_nlp_solve(prob)
  expect_equal(val, log(3 * exp(1)), tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_log_sum_exp.py::TestLogSumExp::test_two
test_that("DNLP log_sum_exp: log_sum_exp(A x) matches DCP optimum", {
  skip_if_not_installed("sparsediff")
  if (!requireNamespace("ipopt", quietly = TRUE) &&
      !requireNamespace("Uno", quietly = TRUE)) skip("No R NLP backend")
  set.seed(3); m <- 50; n <- 10
  A <- matrix(rnorm(m * n), m, n)
  x <- Variable(n)
  prob <- Problem(Minimize(log_sum_exp(A %*% x)), list(x >= 0, sum_entries(x) == 1))
  invisible(capture.output(psolve(prob, nlp = TRUE, print_level = 0L), type = "output"))
  expect_equal(status(prob), "optimal")
  dnlp_val <- value(prob)
  psolve(prob, solver = "CLARABEL")
  expect_equal(dnlp_val, value(prob), tolerance = 1e-5)
  value(x) <- as.numeric(value(x)) + 0.1   # perturb off the boundary, as CVXPY does
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_log_sum_exp.py::TestLogSumExp::test_three
test_that("DNLP log_sum_exp: log_sum_exp(square(A x)) matches DCP optimum", {
  skip_if_not_installed("sparsediff")
  if (!requireNamespace("ipopt", quietly = TRUE) &&
      !requireNamespace("Uno", quietly = TRUE)) skip("No R NLP backend")
  sizes <- if (.nlp_full()) list(c(50L, 25L), c(300L, 100L)) else list(c(50L, 25L))
  for (mn in sizes) {
    m <- mn[1L]; n <- mn[2L]
    set.seed(11)
    A <- matrix(rnorm(m * n), m, n)
    x <- Variable(n); y <- Variable(n)
    prob <- Problem(Minimize(log_sum_exp(square(A %*% x))),
                    list(x >= 0, x + y == 1, y >= 0))
    invisible(capture.output(psolve(prob, nlp = TRUE, print_level = 0L), type = "output"))
    expect_equal(status(prob), "optimal")
    dnlp_val <- value(prob)
    psolve(prob, solver = "CLARABEL")
    expect_equal(dnlp_val, value(prob), tolerance = 1e-5)
    nlp_derivative_check(prob)
  }
})
