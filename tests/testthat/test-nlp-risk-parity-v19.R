## CVXPY 1.9 NLP risk-parity examples (nlp_tests/test_risk_parity.py::TestRiskParity).
##
## Sigma is hardcoded (deterministic), so the assertion -- normalized risk
## contributions equal the target -- is the defining property of a risk-parity
## portfolio and is reproducible verbatim. Solve with IPOPT, assert the risk
## contributions, and run the DerivativeChecker mirror. CVXPY's
## `derivative_test='none'` is an IPOPT-side flag with no R analogue.

.rp_sigma <- function() {
  1e-5 * matrix(c(
    41.16, 22.03, 18.64, -4.74,  6.27, 10.1 , 14.52,  3.18,
    22.03, 58.57, 32.92, -5.04,  4.02,  3.7 , 26.76,  2.17,
    18.64, 32.92, 81.02,  0.53,  6.05,  2.02, 25.52,  1.56,
    -4.74, -5.04,  0.53, 20.6 ,  2.52,  0.57,  0.2 ,  3.6 ,
    6.27,  4.02,  6.05,  2.52, 10.13,  2.59,  4.32,  3.13,
    10.1 ,  3.7 ,  2.02,  0.57,  2.59, 22.89,  3.97,  3.26,
    14.52, 26.76, 25.52,  0.2 ,  4.32,  3.97, 29.91,  3.25,
    3.18,  2.17,  1.56,  3.6 ,  3.13,  3.26,  3.25, 13.63), 8, 8, byrow = TRUE)
}

## @cvxpy nlp_tests/test_risk_parity.py::TestRiskParity::test_vanilla_risk_parity_formulation_one
test_that("risk parity: vanilla formulation one", {
  Sigma <- .rp_sigma(); n <- 8; risk_target <- rep(1 / n, n)
  w <- Variable(n, nonneg = TRUE, name = "w"); t <- Variable(n, name = "t")
  value(w) <- rep(1 / n, n)
  term1 <- sum_entries((w^2) * (t^2)) / quad_form(w, Sigma)
  term2 <- (sqrt(sum(risk_target^2)))^2 * quad_form(w, Sigma)
  term3 <- -2 * sum_entries(risk_target * (w * t))
  prob <- Problem(Minimize(term1 + term2 + term3),
                  list(sum_entries(w) == 1, t == Sigma %*% w))
  .de_ipopt_solve(prob)
  wv <- as.numeric(value(w))
  rc <- wv * as.numeric(Sigma %*% wv); rc <- rc / sum(rc)
  expect_lt(sqrt(sum((rc - risk_target)^2)), 1e-5)
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_risk_parity.py::TestRiskParity::test_vanilla_risk_parity_formulation_two
## CVXPY's formulation two is an unimplemented placeholder (`pass`); mirror it.
test_that("risk parity: vanilla formulation two (placeholder)", {
  succeed("CVXPY test_vanilla_risk_parity_formulation_two is a `pass` stub")
})

## @cvxpy nlp_tests/test_risk_parity.py::TestRiskParity::test_group_risk_parity_formulation_one
test_that("risk parity: group formulation one", {
  Sigma <- .rp_sigma(); n <- 8; b <- c(0.4, 0.6)
  w <- Variable(n, nonneg = TRUE, name = "w"); t <- Variable(n, name = "t")
  value(w) <- rep(1 / n, n)
  groups <- list(c(1L, 2L, 6L), c(4L, 5L, 3L, 7L, 8L))  # CVXPY 0-based + 1
  term1 <- 0; term2 <- 0; term3 <- 0
  for (k in seq_along(groups)) {
    g <- groups[[k]]
    term1 <- term1 + square(sum_entries(w[g] * t[g])) / quad_form(w, Sigma)
    term2 <- term2 + (abs(b[k]))^2 * quad_form(w, Sigma)
    term3 <- term3 - 2 * b[k] * sum_entries(w[g] * t[g])
  }
  prob <- Problem(Minimize(term1 + term2 + term3),
                  list(sum_entries(w) == 1, t == Sigma %*% w))
  .de_ipopt_solve(prob)
  wv <- as.numeric(value(w))
  rc <- wv * as.numeric(Sigma %*% wv); rc <- rc / sum(rc)
  rc_g <- vapply(groups, function(g) sum(rc[g]), numeric(1))
  expect_lt(sqrt(sum((rc_g - b)^2)), 1e-5)
  nlp_derivative_check(prob)
})

## @cvxpy nlp_tests/test_risk_parity.py::TestRiskParity::test_group_risk_parity_formulation_two
test_that("risk parity: group formulation two", {
  Sigma <- .rp_sigma(); n <- 8; b <- c(0.4, 0.6)
  w <- Variable(n, nonneg = TRUE, name = "w"); t <- Variable(n, name = "t")
  value(w) <- rep(1 / n, n)
  groups <- list(c(1L, 2L, 6L), c(4L, 5L, 3L, 7L, 8L))
  obj <- 0
  for (k in seq_along(groups)) {
    g <- groups[[k]]
    obj <- obj + square(sum_entries(w[g] * t[g]) / quad_form(w, Sigma) - b[k])
  }
  prob <- Problem(Minimize(obj), list(sum_entries(w) == 1, t == Sigma %*% w))
  .de_ipopt_solve(prob)
  wv <- as.numeric(value(w))
  rc <- wv * as.numeric(Sigma %*% wv); rc <- rc / sum(rc)
  rc_g <- vapply(groups, function(g) sum(rc[g]), numeric(1))
  expect_lt(sqrt(sum((rc_g - b)^2)), 1e-5)
  nlp_derivative_check(prob)
})
