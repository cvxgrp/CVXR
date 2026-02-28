## Tests for Phase 7b-ii: PowerApprox, PnormApprox, GeoMean, GeoMeanApprox,
## decomp_quad, and 6 new canonicalizers (power_approx, pnorm_approx,
## geo_mean_exact, geo_mean_approx, huber, quad_form).
## All expected values verified against CVXPY 1.8.1 (branch claude).

library(testthat)

# ═══════════════════════════════════════════════════════════════════
# PowerApprox class
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("PowerApprox: construction and properties", {
  x <- Variable(3)
  pa <- power(x, 2.5)  # approx=TRUE by default
  expect_s3_class(pa, "CVXR::PowerApprox")
  expect_equal(as.numeric(pa@p_used), 2.5)
  expect_equal(pa@approx_error, 0.0)
  expect_false(is.null(pa@w))
  ## CVXPY: w = (2/5, 3/5)
  expect_equal(as.numeric(pa@w), c(2/5, 3/5))
})

## @cvxpy NONE
test_that("PowerApprox: p=2 gives integer p_used", {
  x <- Variable(3)
  pa <- power(x, 2)
  expect_s3_class(pa, "CVXR::PowerApprox")
  expect_equal(as.numeric(pa@p_used), 2)
  expect_equal(pa@approx_error, 0.0)
})

## @cvxpy NONE
test_that("PowerApprox: p=0.5 (fractional)", {
  x <- Variable(3)
  pa <- power(x, 0.5)
  expect_s3_class(pa, "CVXR::PowerApprox")
  expect_equal(as.numeric(pa@p_used), 0.5)
  expect_equal(pa@approx_error, 0.0)
  ## concave for 0 < p < 1
  expect_true(is_concave(pa))
})

## @cvxpy NONE
test_that("PowerApprox: irrational exponent has approx_error", {
  x <- Variable(3)
  pa <- power(x, pi, max_denom = 8L)
  expect_s3_class(pa, "CVXR::PowerApprox")
  ## With small max_denom, approximation should have nonzero error
  expect_true(pa@approx_error > 0)
})

## @cvxpy NONE
test_that("Power (exact) vs PowerApprox: factory dispatch", {
  x <- Variable(3)
  pe <- power(x, 2.5, approx = FALSE)
  pa <- power(x, 2.5, approx = TRUE)
  expect_s3_class(pe, "CVXR::Power")
  expect_false(S7_inherits(pe, PowerApprox))
  expect_s3_class(pa, "CVXR::PowerApprox")
})

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("PowerApprox: DCP curvature inherited from Power", {
  x <- Variable(3, nonneg = TRUE)
  ## p=2.5 > 1: convex
  expect_true(is_convex(power(x, 2.5)))
  ## p=0.5, 0<p<1: concave
  expect_true(is_concave(power(x, 0.5)))
  ## p=-1, p<0: convex
  expect_true(is_convex(power(x, -1)))
})

## @cvxpy test_atoms.py::TestAtoms::test_power
test_that("PowerApprox: sign is always nonneg (except p=1)", {
  x <- Variable(3)
  expect_true(is_nonneg(power(x, 2.5)))
  expect_true(is_nonneg(power(x, 0.5)))
  expect_true(is_nonneg(power(x, -1)))
})

# ═══════════════════════════════════════════════════════════════════
# PnormApprox class
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("PnormApprox: construction and properties", {
  x <- Variable(3)
  pna <- p_norm(x, p = 3)
  expect_s3_class(pna, "CVXR::PnormApprox")
  ## CVXPY: p=3 (exact integer, no approximation error)
  expect_equal(pna@p, 3)
  expect_equal(pna@approx_error, 0.0)
})

## @cvxpy NONE
test_that("PnormApprox: p=2 gives standard Euclidean norm", {
  x <- Variable(3)
  pna <- p_norm(x, p = 2)
  expect_s3_class(pna, "CVXR::PnormApprox")
  expect_equal(pna@p, 2)
})

## @cvxpy NONE
test_that("PnormApprox: fractional p", {
  x <- Variable(3)
  pna <- p_norm(x, p = 1.5)
  expect_s3_class(pna, "CVXR::PnormApprox")
  ## 1.5 = 3/2 is exactly representable
  expect_equal(pna@p, 1.5)
  expect_equal(pna@approx_error, 0.0)
})

## @cvxpy NONE
test_that("Pnorm (exact) vs PnormApprox: factory dispatch", {
  x <- Variable(3)
  pne <- p_norm(x, p = 3, approx = FALSE)
  pna <- p_norm(x, p = 3, approx = TRUE)
  expect_s3_class(pne, "CVXR::Pnorm")
  expect_false(S7_inherits(pne, PnormApprox))
  expect_s3_class(pna, "CVXR::PnormApprox")
})

## @cvxpy NONE
test_that("PnormApprox: special dispatch for p=1 and p=Inf", {
  x <- Variable(3)
  expect_s3_class(p_norm(x, 1), "CVXR::Norm1")
  expect_s3_class(p_norm(x, Inf), "CVXR::NormInf")
})

## @cvxpy test_atoms.py::TestAtoms::test_pnorm
test_that("PnormApprox: DCP curvature", {
  x <- Variable(3)
  ## p >= 1: convex
  expect_true(is_convex(p_norm(x, 3)))
  expect_true(is_convex(p_norm(x, 2)))
  ## p < 1 (and p > 0): concave
  expect_true(is_concave(p_norm(x, 0.5)))
})

# ═══════════════════════════════════════════════════════════════════
# GeoMean class
# ═══════════════════════════════════════════════════════════════════

## @cvxpy test_atoms.py::TestAtoms::test_geo_mean
test_that("GeoMean: construction with uniform weights", {
  x <- Variable(4)
  gm <- geo_mean(x, approx = FALSE)
  expect_s3_class(gm, "CVXR::GeoMean")
  ## Uniform weights: each = 1/4
  expect_equal(gm@w, rep(0.25, 4))
  expect_equal(gm@shape, c(1L, 1L))  # always scalar
  expect_equal(gm@approx_error, 0.0)
})

## @cvxpy NONE
test_that("GeoMean: construction with custom weights", {
  x <- Variable(4)
  gm <- geo_mean(x, p = c(1, 2, 3, 4), approx = FALSE)
  ## w = p / sum(p) = [0.1, 0.2, 0.3, 0.4]
  expect_equal(gm@w, c(0.1, 0.2, 0.3, 0.4))
})

## @cvxpy NONE
test_that("GeoMean: zero weights filter elements", {
  x <- Variable(5)
  gm <- geo_mean(x, p = c(1, 0, 1, 0, 1), approx = FALSE)
  ## Only 3 nonzero weights
  expect_equal(length(gm@w), 3L)
  expect_equal(gm@w, rep(1/3, 3))
})

## @cvxpy NONE
test_that("GeoMean: validation errors", {
  x <- Variable(3)
  ## Negative weights: filtered out (like CVXPY), no error
  ## But all-zero weights error (after filtering, empty → validation fails)
  expect_error(geo_mean(x, p = c(0, 0, 0)))
  ## Mismatched length
  expect_error(geo_mean(x, p = c(1, 2)), "same number")
})

## @cvxpy test_atoms.py::TestAtoms::test_geo_mean
test_that("GeoMean: DCP curvature is concave", {
  x <- Variable(3, nonneg = TRUE)
  gm <- geo_mean(x, approx = FALSE)
  expect_true(is_concave(gm))
  ## Convex only when single weight (degenerates to identity)
  x1 <- Variable(1, nonneg = TRUE)
  gm1 <- geo_mean(x1, approx = FALSE)
  expect_true(is_convex(gm1))
})

## @cvxpy test_atoms.py::TestAtoms::test_geo_mean
test_that("GeoMean: sign is always nonneg", {
  x <- Variable(3)
  gm <- geo_mean(x, approx = FALSE)
  expect_true(is_nonneg(gm))
})

## @cvxpy NONE
test_that("GeoMean: numeric evaluation", {
  x <- Constant(c(2, 8))
  gm <- geo_mean(x, approx = FALSE)
  ## (2 * 8)^0.5 = 4
  expect_equal(as.numeric(value(gm)), 4.0)
})

## @cvxpy NONE
test_that("GeoMean: numeric evaluation with weights", {
  x <- Constant(c(4, 9))
  gm <- geo_mean(x, p = c(1, 2), approx = FALSE)
  ## 4^(1/3) * 9^(2/3) = 4^0.333.. * 9^0.666..
  expected <- 4^(1/3) * 9^(2/3)
  expect_equal(as.numeric(value(gm)), expected, tolerance = 1e-10)
})

# ═══════════════════════════════════════════════════════════════════
# GeoMeanApprox class
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("GeoMeanApprox: construction and properties", {
  x <- Variable(4)
  gma <- geo_mean(x)  # approx=TRUE by default
  expect_s3_class(gma, "CVXR::GeoMeanApprox")
  ## w is bigq, should be 1/4 each
  expect_equal(as.numeric(gma@w), rep(0.25, 4))
  expect_false(is.null(gma@w_dyad))
  expect_false(is.null(gma@tree))
  expect_true(gma@cone_num > 0L)
})

## @cvxpy NONE
test_that("GeoMeanApprox: approx_error for exact weights", {
  x <- Variable(4)
  gma <- geo_mean(x)
  ## Uniform weights are exactly representable
  expect_equal(gma@approx_error, 0.0)
})

## @cvxpy NONE
test_that("GeoMeanApprox: factory dispatch", {
  x <- Variable(3)
  gme <- geo_mean(x, approx = FALSE)
  gma <- geo_mean(x, approx = TRUE)
  expect_s3_class(gme, "CVXR::GeoMean")
  expect_false(S7_inherits(gme, GeoMeanApprox))
  expect_s3_class(gma, "CVXR::GeoMeanApprox")
})

# ═══════════════════════════════════════════════════════════════════
# decomp_quad
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("decomp_quad: PSD matrix", {
  P <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  result <- decomp_quad(P)
  ## Reconstruct: scale * (M1 %*% t(M1) - M2 %*% t(M2))
  reconstructed <- result$scale * (result$M1 %*% t(result$M1))
  expect_equal(reconstructed, P, tolerance = 1e-10)
  expect_equal(ncol(result$M2), 0L)  # No negative eigenvalues
})

## @cvxpy NONE
test_that("decomp_quad: NSD matrix", {
  P <- -matrix(c(2, 0.5, 0.5, 1), 2, 2)  # NSD
  result <- decomp_quad(P)
  reconstructed <- result$scale * (-result$M2 %*% t(result$M2))
  expect_equal(reconstructed, P, tolerance = 1e-10)
  expect_equal(ncol(result$M1), 0L)  # No positive eigenvalues
})

## @cvxpy NONE
test_that("decomp_quad: identity matrix", {
  P <- diag(3)
  result <- decomp_quad(P)
  reconstructed <- result$scale * result$M1 %*% t(result$M1)
  expect_equal(reconstructed, P, tolerance = 1e-10)
})

## @cvxpy NONE
test_that("decomp_quad: zero matrix", {
  P <- matrix(0, 2, 2)
  result <- decomp_quad(P)
  expect_equal(ncol(result$M1), 0L)
  expect_equal(ncol(result$M2), 0L)
  expect_equal(result$scale, 0)
})

# ═══════════════════════════════════════════════════════════════════
# CANON_METHODS registration
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("DCP canonicalization S7 methods are registered for approx atoms", {
  x <- Variable(2)
  expect_true(has_dcp_canon(Power(x, 0.5)))
  expect_true(has_dcp_canon(Pnorm(x, 3)))
  expect_true(has_dcp_canon(GeoMean(x)))
  expect_true(has_dcp_canon(GeoMeanApprox(x)))
  expect_true(has_dcp_canon(Huber(x)))
  expect_true(has_dcp_canon(QuadForm(x, Constant(diag(2)))))
})

# ═══════════════════════════════════════════════════════════════════
# Solve parity tests (CVXPY verified)
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Solve: power(x, 2.5) minimization (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  x <- Variable(name = "x")
  prob <- Problem(Minimize(power(x, 2.5)), list(x >= 2))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  ## CVXPY: 2^2.5 = 5.656854
  expect_equal(value(prob), 5.656854, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Solve: power(x, 0.5) concave maximization", {
  skip_if_not_installed("clarabel")
  x <- Variable(name = "x")
  prob <- Problem(Maximize(power(x, 0.5)), list(x <= 9, x >= 0))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 3.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 9.0, tolerance = 1e-4)
})

## @cvxpy test_problem.py::TestProblem::test_geo_mean
test_that("Solve: geo_mean maximization (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  x <- Variable(3, name = "x")
  prob <- Problem(Maximize(geo_mean(x)), list(sum(x) <= 4, x >= 0))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  ## CVXPY: 4/3 ≈ 1.333333
  expect_equal(value(prob), 4/3, tolerance = 1e-4)
  vals <- as.numeric(value(x))
  expect_equal(vals, rep(4/3, 3), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Solve: geo_mean with weighted mean", {
  skip_if_not_installed("clarabel")
  x <- Variable(2, nonneg = TRUE, name = "x")
  prob <- Problem(Maximize(geo_mean(x, p = c(1, 3))),
                  list(x[1] + x[2] <= 4))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  ## With weights (1,3) → normalized (1/4, 3/4)
  ## Optimal: x[1]=1, x[2]=3 → 1^(1/4) * 3^(3/4) = 3^0.75
  expect_equal(value(prob), 3^0.75, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Solve: huber loss (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  x <- Variable(3, name = "x")
  y <- c(1.0, 2.0, 3.0)
  prob <- Problem(Minimize(sum(huber(x - y))), list())
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  ## CVXPY: optimal at x = y, value = 0
  expect_equal(value(prob), 0.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), y, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Solve: huber with constraints", {
  skip_if_not_installed("clarabel")
  x <- Variable(2, name = "x")
  prob <- Problem(Minimize(sum(huber(x))), list(x >= 5))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  ## huber(5) = 2*1*5 - 1^2 = 9 (M=1 default, |x|=5 > M)
  ## Two variables, each contributes 9
  expect_equal(value(prob), 18.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Solve: quad_form minimization (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  P <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  x <- Variable(2, name = "x")
  prob <- Problem(Minimize(quad_form(x, P)), list(sum(x) == 1))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  ## CVXPY: 0.875
  expect_equal(value(prob), 0.875, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), c(0.25, 0.75), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Solve: quad_form with identity (sum of squares)", {
  skip_if_not_installed("clarabel")
  P <- diag(3)
  x <- Variable(3, name = "x")
  prob <- Problem(Minimize(quad_form(x, P)), list(sum(x) == 3))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  ## x = (1,1,1), x^T I x = 3
  expect_equal(value(prob), 3.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), rep(1, 3), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Solve: p_norm(x, 3) minimization (CVXPY parity)", {
  skip_if_not_installed("clarabel")
  x <- Variable(3, name = "x")
  prob <- Problem(Minimize(p_norm(x, 3)), list(sum(x) >= 3, x >= 0))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  ## CVXPY: 3^(1/3) ≈ 1.44225
  expect_equal(value(prob), 3^(1/3), tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), rep(1, 3), tolerance = 1e-2)
})

## @cvxpy NONE
test_that("Solve: p_norm(x, 0.5) concave", {
  skip_if_not_installed("clarabel")
  x <- Variable(2, nonneg = TRUE, name = "x")
  prob <- Problem(Maximize(p_norm(x, 0.5)), list(sum(x) <= 4))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  ## p < 1: concave. Maximum when x balanced: x = (2, 2)
  ## ||[2,2]||_0.5 = (2^0.5 + 2^0.5)^(1/0.5) = (2*sqrt(2))^2 = 8
  expect_equal(value(prob), 8.0, tolerance = 1e-3)
})

# ═══════════════════════════════════════════════════════════════════
# Exact cone path solve tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("Solve: power exact (PowCone3D) path", {
  skip_if_not_installed("clarabel")
  x <- Variable(name = "x")
  prob <- Problem(Minimize(power(x, 2.5, approx = FALSE)), list(x >= 2))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 5.656854, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("Solve: pnorm exact (PowCone3D) path", {
  skip_if_not_installed("clarabel")
  x <- Variable(3, name = "x")
  prob <- Problem(Minimize(p_norm(x, 3, approx = FALSE)),
                  list(sum(x) >= 3, x >= 0))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 3^(1/3), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("Solve: geo_mean exact (PowConeND) path", {
  skip_if_not_installed("clarabel")
  x <- Variable(3, nonneg = TRUE, name = "x")
  prob <- Problem(Maximize(geo_mean(x, approx = FALSE)),
                  list(sum(x) <= 4))
  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  expect_equal(value(prob), 4/3, tolerance = 1e-4)
})

# ═══════════════════════════════════════════════════════════════════
# DCP analysis tests
# ═══════════════════════════════════════════════════════════════════

## @cvxpy NONE
test_that("DCP: PowerApprox in problems passes is_dcp", {
  x <- Variable(nonneg = TRUE)
  prob <- Problem(Minimize(power(x, 2)), list(x >= 1))
  expect_true(is_dcp(prob))
})

## @cvxpy NONE
test_that("DCP: GeoMean in Maximize is DCP", {
  x <- Variable(3, nonneg = TRUE)
  prob <- Problem(Maximize(geo_mean(x)), list(sum(x) <= 1))
  expect_true(is_dcp(prob))
})

## @cvxpy NONE
test_that("DCP: GeoMean in Minimize is NOT DCP", {
  x <- Variable(3, nonneg = TRUE)
  prob <- Problem(Minimize(geo_mean(x)), list(sum(x) <= 1))
  expect_false(is_dcp(prob))
})

## @cvxpy test_atoms.py::TestAtoms::test_huber
test_that("DCP: Huber is convex", {
  x <- Variable(3)
  h <- huber(x)
  expect_true(is_convex(h))
  expect_false(is_concave(h))
})

## @cvxpy NONE
test_that("DCP: QuadForm with PSD P is convex", {
  P <- matrix(c(2, 0, 0, 1), 2, 2)
  x <- Variable(2)
  qf <- quad_form(x, P)
  expect_true(is_convex(qf))
  expect_true(is_dcp(Problem(Minimize(qf))))
})
