## CVXPY 1.9 DNLP hyperbolic-atom parity (nlp_tests/test_hyperbolic.py).
##
## sinh/tanh/asinh/atanh are smooth atoms; composed with logistic they give
## nonconvex-but-smooth DNLPs. CVXPY solves through IPOPT plus a Python
## DerivativeChecker (dropped here, no R equivalent) and only checks the solve
## reaches an optimum. test_sinh / test_atanh equivalents already live in
## test-nlp-solving-chain.R; this covers the tanh and asinh cases.
## CVXPY exposes these via the cp.nlp submodule; CVXR exposes them directly.

.hyp_nlp_solve <- function(prob) {
  skip_if_not_installed("sparsediff")
  if (!requireNamespace("ipopt", quietly = TRUE) &&
      !requireNamespace("Uno", quietly = TRUE)) {
    skip("No R NLP backend (ipopt or Uno) installed")
  }
  output <- capture.output(val <- psolve(prob, nlp = TRUE), type = "output")
  expect_equal(status(prob), "optimal")
  invisible(val)
}

## @cvxpy nlp_tests/test_hyperbolic.py::TestHyperbolic::test_tanh
test_that("DNLP hyperbolic: sum(tanh(logistic(2x))) solves", {
  n <- 10
  x <- Variable(n)
  value(x) <- rep(1, n)
  prob <- Problem(Minimize(sum_entries(tanh(logistic(x * 2)))),
                  list(x >= 0.1, sum_entries(x) == 10))
  .hyp_nlp_solve(prob)
})

## @cvxpy nlp_tests/test_hyperbolic.py::TestHyperbolic::test_asinh
test_that("DNLP hyperbolic: sum(asinh(logistic(3x))) solves", {
  n <- 10
  x <- Variable(n)
  value(x) <- rep(1, n)
  prob <- Problem(Minimize(sum_entries(asinh(logistic(x * 3)))),
                  list(x >= 0.1, sum_entries(x) == 10))
  .hyp_nlp_solve(prob)
})
