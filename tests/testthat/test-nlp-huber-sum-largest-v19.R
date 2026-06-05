## CVXPY 1.9 DNLP parity for huber / sum_largest / sum_smallest
## (nlp_tests/test_huber_sum_largest.py).
##
## These PWL atoms are convex/concave (not smooth), so a DNLP that uses them
## is canonicalized to a smooth epigraph by Dnlp2Smooth. The epigraph
## variables must be initialized so the NLP backend has a starting point;
## sum_largest_canon (shared with the conic path) sets q to the (k+1)-th
## largest entry and t to the positive part above q on the top-k entries.
## CVXPY runs these through IPOPT plus a Python DerivativeChecker (dropped
## here); the huber case uses a small deterministic instance instead of the
## upstream 100-variable random sweep to keep the suite fast.

.hsl_nlp_solve <- function(prob) {
  skip_if_not_installed("sparsediff")
  if (!requireNamespace("ipopt", quietly = TRUE) &&
      !requireNamespace("Uno", quietly = TRUE)) {
    skip("No R NLP backend (ipopt or Uno) installed")
  }
  output <- capture.output(val <- psolve(prob, nlp = TRUE), type = "output")
  expect_equal(status(prob), "optimal")
  invisible(val)
}

## @cvxpy nlp_tests/test_huber_sum_largest.py::TestSumLargestCanonWarmStart::test_warmstart_feasibility
test_that("sum_largest_canon initializes a feasible warm start", {
  x <- Variable(6)
  value(x) <- c(1, 10, 6, 8, 3, 9)
  res <- sum_largest_canon(SumLargest(x, 2L), list(x))
  obj <- res[[1L]]
  ## obj = SumEntries(t) + 2 * q
  t <- obj@args[[1L]]@args[[1L]]
  q <- obj@args[[2L]]@args[[2L]]

  ## q is the max of the non-top-k entries (x_{[k+1]} = 8), not the max of the
  ## k smallest -- the regression this test guards.
  expect_equal(as.numeric(value(q)), 8)
  expected_t <- numeric(6)
  expected_t[2] <- 10 - 8  # top-k element
  expected_t[6] <- 9 - 8   # top-k element
  expect_equal(as.numeric(value(t)), expected_t)

  xv <- c(1, 10, 6, 8, 3, 9)
  expect_true(all(xv <= as.numeric(value(t)) + as.numeric(value(q)) + 1e-10))
  expect_true(all(as.numeric(value(t)) >= -1e-10))
})

## @cvxpy nlp_tests/test_huber_sum_largest.py::TestNonsmoothNontrivial::test_sum_largest
test_that("DNLP: sum_largest matches the conic solve", {
  x <- Variable(5)
  value(x) <- rep(0.2, 5)
  w <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  prob <- Problem(Minimize(sum_largest(x * w, 2)), list(sum_entries(x) == 1, x >= 0))
  nlp_value <- .hsl_nlp_solve(prob)
  conic_value <- psolve(prob, solver = "CLARABEL")
  expect_equal(nlp_value, conic_value, tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_huber_sum_largest.py::TestNonsmoothNontrivial::test_sum_smallest
test_that("DNLP: sum_smallest matches the conic solve", {
  x <- Variable(5)
  value(x) <- rep(0.2, 5)
  w <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  prob <- Problem(Maximize(sum_smallest(x * w, 2)), list(sum_entries(x) == 1, x >= 0))
  nlp_value <- .hsl_nlp_solve(prob)
  conic_value <- psolve(prob, solver = "CLARABEL")
  expect_equal(nlp_value, conic_value, tolerance = 1e-4)
})

## @cvxpy nlp_tests/test_huber_sum_largest.py::TestNonsmoothNontrivial::test_huber
test_that("DNLP: huber robust regression matches the conic solve", {
  set.seed(1)
  n <- 6; S <- 9
  beta_true <- 5 * rnorm(n)
  X <- matrix(rnorm(n * S), n, S)
  v <- rnorm(S)
  Y <- as.numeric(crossprod(X, beta_true)) + v
  beta <- Variable(n)
  value(beta) <- rep(0, n)
  prob <- Problem(Minimize(sum_entries(huber(t(X) %*% beta - Y, 1))))
  nlp_value <- .hsl_nlp_solve(prob)
  conic_value <- psolve(prob, solver = "CLARABEL")
  expect_equal(nlp_value, conic_value, tolerance = 1e-4)
})
