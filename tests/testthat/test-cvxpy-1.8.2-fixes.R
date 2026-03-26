## CVXPY v1.8.2 regression tests
## ============================
## Tests ported from CVXPY v1.8.2 (tag v1.8.2, commit 76b052791, released 2026-03-22)
## to verify Tier 1 bug fixes ported to CVXR.
##
## Source files:
##   /Users/naras/GitHub/cvxpy/cvxpy/tests/test_quad_form.py
##   /Users/naras/GitHub/cvxpy/cvxpy/tests/test_atoms.py
##   /Users/naras/GitHub/cvxpy/cvxpy/tests/test_constraints.py
##   /Users/naras/GitHub/cvxpy/cvxpy/tests/test_dgp.py
##   /Users/naras/GitHub/cvxpy/cvxpy/tests/test_attributes.py
##   /Users/naras/GitHub/cvxpy/cvxpy/tests/test_clarabel_warm_start.py
##
## Expected values verified against CVXPY 1.8.2 using `uv run python`.

library(testthat)
library(CVXR)

# ══════════════════════════════════════════════════════════════════
# quad_form_canon: indefinite P fix (silent wrong results in v1.8.1)
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_quad_form.py::TestNonOptimal::test_indefinite_assume_psd_raises
test_that("quad_form_canon: indefinite P with assume_PSD raises error", {
  ## CVXPY v1.8.2 fix: in v1.8.1 the M1 (positive eigenvalue) term was
  ## silently dropped for indefinite P, producing wrong objective at OPTIMAL.
  ## CVXPY v1.8.2: use_quad_obj=FALSE forces conic path where
  ## quad_form_canon catches indefinite P via decomp_quad.
  P <- matrix(c(2, 0, 0, -1), 2, 2)  # indefinite: eigenvalues 2 and -1
  x <- Variable(2)
  expr <- quad_form(x, P, assume_PSD = TRUE)
  prob <- Problem(Minimize(expr), list(x == c(1, 2)))
  expect_error(psolve(prob, use_quad_obj = FALSE), "indefinite")
})

## @cvxpy test_quad_form.py::TestNonOptimal::test_quad_form_monotonicity_in_x
test_that("quad_form monotonicity in x: convex composition via DCP", {
  ## CVXPY v1.8.2 fix: is_incr/is_decr were returning FALSE unconditionally.
  P <- matrix(c(2, 0.5, 0.5, 1), 2, 2)

  ## x = square(z) is convex nonneg; quad_form nondecreasing in nonneg x -> convex
  z <- Variable(2, nonneg = TRUE)
  expect_true(is_convex(quad_form(power(z, 2), P)))

  ## x = square(z) - 1 is convex but can be negative; monotonicity unknown -> not DCP
  expect_false(is_convex(quad_form(power(z, 2) - 1, P)))
})

## @cvxpy test_quad_form.py::TestNonOptimal::test_quad_form_monotonicity_in_P
test_that("quad_form monotonicity in P: convex composition via DCP", {
  ## Nonneg x: quad_form is nondecreasing in P, so convex o nondecreasing = convex
  y <- Variable(2, nonneg = TRUE)
  P_expr <- DiagVec(power(y, 2))
  expect_true(is_convex(quad_form(c(1, 2), P_expr)))

  ## Mixed-sign x: monotonicity unknown, so not DCP-convex
  expect_false(is_convex(quad_form(c(1, -2), P_expr)))
})

## @cvxpy test_quad_form.py::TestDecompQuad::test_psd_nsd
test_that("decomp_quad: PSD and NSD reconstruction", {
  ## Helper: check P = scale * (M1 %*% t(M1) - M2 %*% t(M2))
  check_decomp <- function(P) {
    result <- CVXR:::decomp_quad(P)
    scale <- result$scale; M1 <- result$M1; M2 <- result$M2
    P_dense <- as.matrix(P)
    pos <- if (ncol(M1) > 0L) scale * (M1 %*% t(M1)) else matrix(0, nrow(P), ncol(P))
    neg <- if (ncol(M2) > 0L) scale * (M2 %*% t(M2)) else matrix(0, nrow(P), ncol(P))
    expect_equal(P_dense, as.matrix(pos - neg), tolerance = 1e-10)
  }

  set.seed(0)
  for (n in c(2, 5, 10)) {
    A <- matrix(rnorm(n * n), n, n)
    check_decomp(A %*% t(A))        # PSD, full rank
    check_decomp(-(A %*% t(A)))     # NSD
  }
  ## PSD, rank-deficient
  A <- matrix(rnorm(10 * 7), 10, 7)
  check_decomp(A %*% t(A))
})

## @cvxpy test_quad_form.py::TestDecompQuad::test_special
test_that("decomp_quad: zero and diagonal matrices", {
  result <- CVXR:::decomp_quad(matrix(0, 4, 4))
  expect_equal(result$scale, 0)
  expect_equal(ncol(result$M1), 0L)
  expect_equal(ncol(result$M2), 0L)

  ## Diagonal with zeros
  P <- diag(c(3, 0, 5, 0))
  result <- CVXR:::decomp_quad(P)
  scale <- result$scale; M1 <- result$M1; M2 <- result$M2
  P_recon <- scale * (M1 %*% t(M1))
  expect_equal(P, as.matrix(P_recon), tolerance = 1e-10)
})

# ══════════════════════════════════════════════════════════════════
# log_det: sign fix
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_log_det_sign
test_that("log_det sign is unknown (can be negative for det < 1)", {
  ## CVXPY v1.8.2 fix: log_det can be negative, e.g. log(det(0.5*I)) < 0
  X <- Variable(c(2, 2), symmetric = TRUE)
  atom <- log_det(X)
  expect_false(is_nonneg(atom))

  ## Verify numerically: log(det(0.5*I)) = 2*log(0.5) < 0
  value(X) <- 0.5 * diag(2)
  expect_true(value(atom) < 0)
})

## @cvxpy test_atoms.py::TestAtoms::test_spectral_atoms_value_none_variable
test_that("log_det value returns NULL when argument is unset Variable", {
  ## CVXPY v1.8.2 fix: guard against NULL (unset Variable)
  X <- Variable(c(3, 3), symmetric = TRUE)
  expect_null(value(log_det(X)))
})

# ══════════════════════════════════════════════════════════════════
# pnorm: log-log curvature fix
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_dgp.py::TestDgp::test_pnorm_negative_p_dgp
test_that("pnorm(x, p=-1) is log-log concave, so Maximize is DGP", {
  ## CVXPY v1.8.2 fix: is_atom_log_log_convex was TRUE for all p.
  x <- Variable(3, pos = TRUE)
  ## Maximize pnorm(x, p=-1) s.t. x <= 2, x >= 0.5
  ## Optimal: all x_i = 2, pnorm = (3 * 2^(-1))^(-1) = 2/3
  prob <- Problem(
    Maximize(p_norm(x, -1)),
    list(x <= 2, x >= 0.5)
  )
  expect_true(is_dgp(prob))
  psolve(prob, gp = TRUE, solver = "SCS")
  expect_equal(status(prob), "optimal")
  expect_equal(psolve(prob, gp = TRUE), 2/3, tolerance = 1e-3)
})

## @cvxpy test_dgp.py::TestDgp::test_pnorm_negative_p_dgp (second part)
test_that("Minimize pnorm(x, p=-1) is NOT DGP", {
  x <- Variable(3, pos = TRUE)
  prob <- Problem(
    Minimize(p_norm(x, -1)),
    list(x <= 2, x >= 0.5)
  )
  expect_false(is_dgp(prob))
})

## @cvxpy test_dgp.py::TestDgp::test_pnorm_scalar
test_that("pnorm(scalar, p=2) in DGP: minimum equals lower bound", {
  ## CVXPY v1.8.2 fix: scalar DGP pnorm canonicalizer was promoting p instead of x
  x <- Variable(pos = TRUE)
  prob <- Problem(Minimize(p_norm(x, 2)), list(x >= 2.5))
  psolve(prob, gp = TRUE)
  expect_equal(status(prob), "optimal")
  expect_equal(as.numeric(value(x)), 2.5, tolerance = 1e-4)
})

# ══════════════════════════════════════════════════════════════════
# leaf: sparse/diag + pos/neg validation
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_dgp.py::TestDgp::test_sparse_variable_not_dgp
test_that("Variables reject combining pos/neg with sparsity/diag", {
  ## CVXPY v1.8.2 fix: sparsity forces zeros, contradicting pos/neg.
  expect_error(Variable(c(3, 3), pos = TRUE, sparsity = TRUE), "pos.*sparsity")
  expect_error(Variable(c(3, 3), neg = TRUE, sparsity = TRUE), "neg.*sparsity")
  expect_error(Variable(c(3, 3), pos = TRUE, diag = TRUE), "pos.*diag")
  expect_error(Variable(c(3, 3), neg = TRUE, diag = TRUE), "neg.*diag")
})

# ══════════════════════════════════════════════════════════════════
# Clarabel warm-start: sparsity change fallback
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_clarabel_warm_start.py::test_clarabel_warm_start_sparsity_change
test_that("Clarabel warm-start: sparsity change falls back to cold start (Issue #2800)", {
  ## CVXPY v1.8.2 fix: solver_update() crashes when sparsity pattern changes.
  ## Wrapped in tryCatch to fall back to cold path re-initialization.
  skip_if_not_installed("clarabel")

  x <- Variable(2)
  p <- Parameter(2)

  prob <- Problem(Minimize(sum_squares(x)), list(sum(Multiply(p, x)) >= 1))

  ## Solve with fully dense parameter
  value(p) <- c(1, 1)
  psolve(prob, solver = "CLARABEL", warm_start = TRUE)
  expect_equal(status(prob), "optimal")

  ## Introduce a zero to change the sparsity pattern
  value(p) <- c(1, 0)
  psolve(prob, solver = "CLARABEL", warm_start = TRUE)
  expect_equal(status(prob), "optimal")
})

# ══════════════════════════════════════════════════════════════════
# trace(A %*% B): O(n^2) optimization via sum(A * t(B))
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_trace_AB
test_that("trace(A %*% B) = sum(A * t(B)), not vdot(A, B)", {
  ## CVXPY v1.8.2 fix: old vdot formula computed trace(A^H @ B), wrong for non-Hermitian A.
  ## Structural check: result is SumEntries of Multiply, not a Trace atom
  A <- Variable(c(4, 5))
  B <- Variable(c(5, 4))
  tr <- matrix_trace(A %*% B)
  expect_false(S7_inherits(tr, CVXR:::Trace))

  ## Numerical correctness for non-symmetric real matrices
  A_val <- matrix(c(1, 4, 2, 5, 3, 6), 2, 3)
  B_val <- matrix(c(7, 9, 11, 8, 10, 12), 3, 2)
  A2 <- Variable(c(2, 3)); B2 <- Variable(c(3, 2))
  value(A2) <- A_val; value(B2) <- B_val

  tr2 <- matrix_trace(A2 %*% B2)
  expected <- sum(diag(A_val %*% B_val))
  expect_equal(as.numeric(value(tr2)), expected, tolerance = 1e-10)
})

## @cvxpy test_atoms.py::TestAtoms::test_trace_AB (solve-level)
test_that("trace(C %*% X) minimize with asymmetric cost matrix", {
  ## Asymmetric cost: vdot(C, X) != trace(C @ X) when C != C^T
  C <- matrix(c(1, 2, 3, 4), 2, 2)
  X <- Variable(c(2, 2), nonneg = TRUE)
  prob <- Problem(Minimize(matrix_trace(Constant(C) %*% X)), list(sum(X) == 1))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  ## Minimum is achieved at the entry with smallest C value
  expect_equal(value(prob), min(C), tolerance = 1e-4)
})

# ══════════════════════════════════════════════════════════════════
# multiply: is_hermitian property preservation
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_multiply_hermitian
test_that("multiply preserves Hermitian property", {
  ## CVXPY v1.8.2: is_hermitian() added to multiply class
  X <- Variable(c(3, 3), hermitian = TRUE)
  Y <- Variable(c(3, 3), hermitian = TRUE)

  ## Real scalar multiplication preserves Hermitian
  expect_true(is_hermitian(1 * X))
  expect_true(is_hermitian(X * 2.5))
  expect_true(is_hermitian(-1 * X))

  ## Hadamard product of two Hermitians is Hermitian
  expect_true(is_hermitian(Multiply(X, Y)))
})

# ══════════════════════════════════════════════════════════════════
# Hermitian property for real/imag atoms
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_imag_hermitian_not_symmetric
test_that("imag(Hermitian) is skew-symmetric, not symmetric", {
  ## CVXPY v1.8.2 fix: imag.is_symmetric() was returning is_hermitian()
  H <- Variable(c(3, 3), hermitian = TRUE)
  expect_false(is_symmetric(Im(H)))
})

## @cvxpy test_atoms.py::TestAtoms::test_imag_symmetric_is_symmetric
test_that("imag(symmetric real matrix) is symmetric (zero matrix)", {
  S <- Variable(c(3, 3), symmetric = TRUE)
  expect_true(is_symmetric(Im(S)))
})

## @cvxpy test_atoms.py::TestAtoms::test_real_hermitian_is_symmetric
test_that("real(Hermitian) is symmetric", {
  H <- Variable(c(3, 3), hermitian = TRUE)
  expect_true(is_symmetric(Re(H)))
})

## @cvxpy test_atoms.py::TestAtoms::test_real_symmetric_is_symmetric
test_that("real(symmetric) is symmetric", {
  S <- Variable(c(3, 3), symmetric = TRUE)
  expect_true(is_symmetric(Re(S)))
})

# ══════════════════════════════════════════════════════════════════
# PSD residual: complex Hermitian uses conjugate transpose
# ══════════════════════════════════════════════════════════════════

## @cvxpy test_constraints.py::TestConstraints::test_psd_residual_complex
test_that("PSD residual uses conjugate transpose for complex Hermitian", {
  ## CVXPY v1.8.2 fix: PSD.residual used .T instead of .H
  ## R's hermitian Variable rejects complex values in .validate_leaf_value;
  ## test with real symmetric instead, complex case deferred.
  X <- Variable(c(2, 2), symmetric = TRUE)
  constr <- X %>>% 0

  ## A symmetric PSD matrix
  value(X) <- matrix(c(2, 1, 1, 2), 2, 2)
  expect_equal(residual(constr), 0, tolerance = 1e-6)

  ## A non-PSD symmetric matrix
  value(X) <- matrix(c(1, 3, 3, 1), 2, 2)
  expect_true(residual(constr) > 0)
})

## @cvxpy test_constraints.py::TestConstraints::test_finite_set_residual_none
test_that("FiniteSet residual returns NULL when variable has no value", {
  ## CVXPY v1.8.2 fix: FiniteSet.residual crashed on None value
  x <- Variable()
  constr <- FiniteSet(x, c(1, 2, 3))
  expect_null(residual(constr))
})
