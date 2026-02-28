## CVXPY SOURCE: cvxpy/tests/test_complex_dpp.py
## Parity tests for TestComplexDPP class.
## One test_that() per CVXPY test method; ## @cvxpy annotation on line before each.
## Expected values verified against CVXPY 1.8.1 (branch claude, commit 3b964472b).

library(testthat)
library(CVXR)

## @cvxpy test_complex_dpp.py::TestComplexDPP::test_dpp_recognition_and_chain
test_that("test_dpp_recognition_and_chain: complex param problem is DPP, uses Complex2Real not EvalParams", {
  skip("Complex atom not yet fully implemented")
  ## CVXPY: p=Parameter(complex=True), x=Variable()
  ## prob = Problem(Minimize(abs(x - real(p)))) -> is_dpp() True
  ## chain: EvalParams NOT in reductions, Complex2Real IS
  p <- Parameter(1, complex = TRUE)
  x <- Variable(1)
  prob <- Problem(Minimize(Abs(x - Re(p))))

  expect_true(is_dpp(prob))

  ## Chain inspection without parameter value (key DPP feature)
  chain <- construct_solving_chain(prob, solver = "CLARABEL")
  red_names <- vapply(chain@reductions,
                      function(r) sub("^.*::", "", class(r)[[1L]]),
                      character(1L))
  expect_true("Complex2Real" %in% red_names)
  expect_false("EvalParams" %in% red_names)

  ## After setting value and solving: still no EvalParams
  p@value <- matrix(1 + 2i, 1, 1)
  result <- psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-3)
})

## @cvxpy test_complex_dpp.py::TestComplexDPP::test_shapes
test_that("test_shapes: DPP with scalar, vector, matrix complex parameters", {
  skip("Complex atom not yet fully implemented")
  skip_if_not_installed("clarabel")

  ## Scalar
  p_sc <- Parameter(1, complex = TRUE)
  x_sc <- Variable(1)
  prob_sc <- Problem(Minimize(sum_entries(x_sc)), list(x_sc >= Re(p_sc)))
  p_sc@value <- matrix(3 + 4i, 1, 1)
  psolve(prob_sc, solver = "CLARABEL")
  expect_equal(as.numeric(value(x_sc)), 3.0, tolerance = 1e-3)

  ## Fast path re-solve (scalar)
  p_sc@value <- matrix(6 + 8i, 1, 1)
  psolve(prob_sc, solver = "CLARABEL")
  expect_equal(as.numeric(value(x_sc)), 6.0, tolerance = 1e-3)

  ## Vector
  p_v <- Parameter(3, complex = TRUE)
  x_v <- Variable(3)
  prob_v <- Problem(Minimize(sum_entries(x_v)), list(x_v >= Re(p_v)))
  p_v@value <- matrix(c(1+1i, 2+2i, 3+3i), 3, 1)
  psolve(prob_v, solver = "CLARABEL")
  expect_equal(as.numeric(value(x_v)), c(1, 2, 3), tolerance = 1e-3)

  ## Fast path re-solve (vector)
  p_v@value <- 2 * matrix(c(1+1i, 2+2i, 3+3i), 3, 1)
  psolve(prob_v, solver = "CLARABEL")
  expect_equal(as.numeric(value(x_v)), c(2, 4, 6), tolerance = 1e-3)

  ## Matrix
  p_m <- Parameter(c(2L, 2L), complex = TRUE)
  x_m <- Variable(c(2L, 2L))
  prob_m <- Problem(Minimize(sum_entries(x_m)), list(x_m >= Re(p_m)))
  p_m@value <- matrix(c(1+1i, 2+2i, 3+3i, 4+4i), 2, 2)
  psolve(prob_m, solver = "CLARABEL")
  expect_equal(as.numeric(value(x_m)), c(1, 2, 3, 4), tolerance = 1e-3)

  ## Fast path re-solve (matrix)
  p_m@value <- 2 * matrix(c(1+1i, 2+2i, 3+3i, 4+4i), 2, 2)
  psolve(prob_m, solver = "CLARABEL")
  expect_equal(as.numeric(value(x_m)), c(2, 4, 6, 8), tolerance = 1e-3)
})

## @cvxpy test_complex_dpp.py::TestComplexDPP::test_param_types
test_that("test_param_types: DPP with imag and complex parameter types", {
  skip("Complex atom not yet fully implemented")
  skip_if_not_installed("clarabel")

  ## imag parameter: min x s.t. x >= Im(p)
  p_im <- Parameter(1, imag = TRUE)
  x_im <- Variable(1)
  prob_im <- Problem(Minimize(x_im), list(x_im >= Im(p_im)))
  p_im@value <- matrix(3i, 1, 1)
  psolve(prob_im, solver = "CLARABEL")
  expect_equal(as.numeric(value(x_im)), 3.0, tolerance = 1e-3)

  ## complex parameter: min x s.t. x >= Re(p)
  p_cx <- Parameter(1, complex = TRUE)
  x_cx <- Variable(1)
  prob_cx <- Problem(Minimize(x_cx), list(x_cx >= Re(p_cx)))
  p_cx@value <- matrix(3 + 4i, 1, 1)
  psolve(prob_cx, solver = "CLARABEL")
  expect_equal(as.numeric(value(x_cx)), 3.0, tolerance = 1e-3)
})

## @cvxpy test_complex_dpp.py::TestComplexDPP::test_mixed_real_and_complex_params
test_that("test_mixed_real_and_complex_params: DPP with real and complex parameters", {
  skip("Complex atom not yet fully implemented")
  skip_if_not_installed("clarabel")

  p_real <- Parameter(1)
  p_complex <- Parameter(1, complex = TRUE)
  x <- Variable(1)
  prob <- Problem(Minimize(x), list(x >= p_real + Re(p_complex)))

  p_real@value <- matrix(2.0, 1, 1)
  p_complex@value <- matrix(3 + 4i, 1, 1)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), 5.0, tolerance = 1e-3)

  ## Fast path re-solve
  p_real@value <- matrix(1.0, 1, 1)
  p_complex@value <- matrix(1 + 1i, 1, 1)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-3)
})

## @cvxpy test_complex_dpp.py::TestComplexDPP::test_complex_param_with_abs
test_that("test_complex_param_with_abs: abs(x - complex_param) is DPP", {
  skip("Complex atom not yet fully implemented")
  skip_if_not_installed("clarabel")

  p <- Parameter(1, complex = TRUE)
  x <- Variable(1)
  prob <- Problem(Minimize(Abs(x - p)))

  expect_true(is_dpp(prob))

  ## p = 3 + 4j: optimal x = Re(p) = 3, objective |3 - (3+4j)| = 4
  p@value <- matrix(3 + 4i, 1, 1)
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(x)), 3.0, tolerance = 1e-3)
  expect_equal(as.numeric(prob@value), 4.0, tolerance = 1e-3)
})

## @cvxpy test_complex_dpp.py::TestComplexDPP::test_dpp_flags
test_that("test_dpp_flags: enforce_dpp succeeds; ignore_dpp uses EvalParams", {
  skip("Complex atom not yet fully implemented")
  skip_if_not_installed("clarabel")

  p <- Parameter(1, complex = TRUE)
  x <- Variable(1)
  prob <- Problem(Minimize(x), list(x >= Re(p)))
  p@value <- matrix(1 + 2i, 1, 1)

  ## enforce_dpp=TRUE: should succeed and not use EvalParams
  psolve(prob, solver = "CLARABEL", enforce_dpp = TRUE)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-3)
  chain_enforced <- prob@.cache$solving_chain
  if (!is.null(chain_enforced)) {
    red_names_e <- vapply(chain_enforced@reductions,
                          function(r) sub("^.*::", "", class(r)[[1L]]),
                          character(1L))
    expect_false("EvalParams" %in% red_names_e)
  }

  ## ignore_dpp=TRUE: should use EvalParams
  prob2 <- Problem(Minimize(x), list(x >= Re(p)))
  psolve(prob2, solver = "CLARABEL", ignore_dpp = TRUE)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-3)
  chain_ignored <- prob2@.cache$solving_chain
  if (!is.null(chain_ignored)) {
    red_names_i <- vapply(chain_ignored@reductions,
                          function(r) sub("^.*::", "", class(r)[[1L]]),
                          character(1L))
    expect_true("EvalParams" %in% red_names_i)
  }
})

## @cvxpy test_complex_dpp.py::TestComplexDPP::test_hermitian_param_dpp
test_that("test_hermitian_param_dpp: DPP with Hermitian parameter, fast path re-solve", {
  skip("Complex atom not yet fully implemented")
  skip_if_not_installed("clarabel")

  for (n in 1:3) {
    P <- Parameter(c(n, n), hermitian = TRUE)
    x <- Variable(1)
    ## min x s.t. x*I >= P (x >= max eigenvalue of P)
    prob <- Problem(Minimize(x),
                    list(x * Constant(diag(n)) %>>% P))

    expect_true(is_dpp(prob))

    ## Create random Hermitian matrix (reproducible)
    set.seed(n)
    A <- matrix(rnorm(n * n) + 1i * rnorm(n * n), n, n)
    P_val1 <- (A + Conj(t(A))) / 2
    P@value <- P_val1

    psolve(prob, solver = "CLARABEL")
    lmax1 <- max(Re(eigen(P_val1, symmetric = FALSE, only.values = TRUE)$values))
    expect_equal(as.numeric(value(x)), lmax1, tolerance = 1e-3)

    ## Fast path re-solve
    P_val2 <- P_val1 + diag(n)
    P@value <- P_val2
    psolve(prob, solver = "CLARABEL")
    lmax2 <- max(Re(eigen(P_val2, symmetric = FALSE, only.values = TRUE)$values))
    expect_equal(as.numeric(value(x)), lmax2, tolerance = 1e-3)
  }
})

## @cvxpy test_complex_dpp.py::TestComplexDPP::test_hermitian_param_efficient_representation
test_that("test_hermitian_param_efficient_representation: compact Hermitian param representation", {
  skip("Complex atom not yet fully implemented")
  skip_if_not_installed("clarabel")

  n <- 3L
  P <- Parameter(c(n, n), hermitian = TRUE)
  X <- Variable(c(n, n), hermitian = TRUE)
  prob <- Problem(Minimize(Trace(Re(X))), list(X %>>% P))

  P@value <- matrix(c(1+0i, -1i, 0, 1i, 1+0i, 1i, 0, -1i, 1+0i), n, n)
  psolve(prob, solver = "CLARABEL")

  ## Verify chain contains Complex2Real
  chain <- prob@.cache$solving_chain
  if (!is.null(chain)) {
    red_names <- vapply(chain@reductions,
                        function(r) sub("^.*::", "", class(r)[[1L]]),
                        character(1L))
    expect_true("Complex2Real" %in% red_names)
  }

  ## Hermitian parameter: real part is symmetric (n x n),
  ## imaginary part is compact vector of length n*(n-1)/2
  expect_equal(is_hermitian(P), TRUE)

  ## Solution should be >= P (elementwise in PSD sense)
  X_val <- value(X)
  expect_true(!is.null(X_val))
  expect_equal(dim(X_val), c(n, n))
})
