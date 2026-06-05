## CVXPY 1.9 NLP Gaussian-MLE stress tests (nlp_tests/test_ML_Gaussian_stress.py).
##
## Three equivalent formulations of the Gaussian maximum-likelihood estimate of
## sigma (and mu). The optimum is the closed-form MLE -- sigma_opt =
## ||data||/sqrt(n) (zero mean) or ||data - mean(data)||/sqrt(n) (nonzero mean),
## mu_opt = mean(data) -- a data-independent invariant, so R-generated data
## verifies the same property CVXPY checks. Each (n, method) is solved with IPOPT
## and run through the DerivativeChecker mirror. The size sweep is reduced by
## default and runs CVXPY's full np.arange(2, 100, 5) under CVXR_NLP_FULL=1.

.mle_sizes <- function() if (.nlp_full()) seq(2L, 97L, by = 5L) else c(3L, 22L, 62L)
.mle_tol <- 1e-3
.mle_rel <- function(v, opt) abs(v - opt) / max(1, abs(opt))

## @cvxpy nlp_tests/test_ML_Gaussian_stress.py::TestStressMLE::test_zero_mean
test_that("MLE stress: zero mean, three formulations", {
  if (!.de_have_backend()) skip("No R NLP backend")
  skip_if_not_installed("ipopt")
  for (n in .mle_sizes()) {
    set.seed(n)
    data <- rnorm(n)
    sigma_opt <- (1 / sqrt(n)) * sqrt(sum(data^2))
    res <- sum(data^2)
    for (method in 1:3) {
      if (method == 1L) {
        sigma <- Variable(1, nonneg = TRUE)
        obj <- (n / 2) * log(2 * pi * square(sigma)) + (1 / (2 * square(sigma))) * res
      } else if (method == 2L) {
        sigma2 <- Variable(1, nonneg = TRUE)
        obj <- (n / 2) * log(2 * pi * sigma2) + (1 / (2 * sigma2)) * res
        sigma <- sigma2  # solved value square-rooted below
      } else {
        sigma <- Variable(1, nonneg = TRUE)
        obj <- n * log(sqrt(2 * pi) * sigma) + (1 / (2 * square(sigma))) * res
      }
      prob <- Problem(Minimize(obj))
      .de_ipopt_solve(prob)
      nlp_derivative_check(prob)
      sval <- if (method == 2L) sqrt(as.numeric(value(sigma))) else as.numeric(value(sigma))
      expect_lt(.mle_rel(sval, sigma_opt), .mle_tol)
    }
  }
})

## @cvxpy nlp_tests/test_ML_Gaussian_stress.py::TestStressMLE::test_nonzero_mean
test_that("MLE stress: nonzero mean, three formulations", {
  if (!.de_have_backend()) skip("No R NLP backend")
  skip_if_not_installed("ipopt")
  for (n in .mle_sizes()) {
    set.seed(n)
    data <- rnorm(n)
    sigma_opt <- (1 / sqrt(n)) * sqrt(sum((data - mean(data))^2))
    mu_opt <- mean(data)
    for (method in 1:3) {
      mu <- Variable(1, name = "mu")
      res <- sum_entries(square(data - mu))
      if (method == 1L) {
        sigma <- Variable(1, nonneg = TRUE)
        obj <- (n / 2) * log(2 * pi * square(sigma)) + (1 / (2 * square(sigma))) * res
      } else if (method == 2L) {
        sigma2 <- Variable(1, nonneg = TRUE)
        obj <- (n / 2) * log(2 * pi * sigma2) + (1 / (2 * sigma2)) * res
        sigma <- sigma2
      } else {
        sigma <- Variable(1, nonneg = TRUE)
        obj <- n * log(sqrt(2 * pi) * sigma) + (1 / (2 * square(sigma))) * res
      }
      prob <- Problem(Minimize(obj))
      .de_ipopt_solve(prob)
      nlp_derivative_check(prob)
      sval <- if (method == 2L) sqrt(as.numeric(value(sigma))) else as.numeric(value(sigma))
      expect_lt(.mle_rel(sval, sigma_opt), .mle_tol)
      expect_lt(.mle_rel(as.numeric(value(mu)), mu_opt), .mle_tol)
    }
  }
})
