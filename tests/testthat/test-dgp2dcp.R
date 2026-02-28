## Tests for Dgp2Dcp reduction (Phase DGP-1)
## CVXPY reference values verified via `uv run python`

## @cvxpy NONE
test_that("basic monomial minimization (box optimization)", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(pos = TRUE)
  obj <- Minimize(1 / (x * y * z))
  constr <- list(2*x*y + 2*x*z + 2*y*z <= 1, x >= 2*y)
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 15.5884571697, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 0.5773580361, tolerance = 1e-4)
  expect_equal(as.numeric(value(y)), 0.2886790022, tolerance = 1e-4)
  expect_equal(as.numeric(value(z)), 0.3848898473, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("monomial maximization", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  obj <- Maximize(x * y)
  constr <- list(2*x + y <= 1)
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 0.125, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 0.25, tolerance = 1e-3)
  expect_equal(as.numeric(value(y)), 0.50, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("division", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  obj <- Minimize(x / y)
  constr <- list(x >= 1, y <= 2)
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 0.5, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(y)), 2.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("power atom", {
  x <- Variable(pos = TRUE)
  obj <- Minimize(x^3)
  constr <- list(x >= 2)
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 8.0, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("equality constraint (monomial)", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  obj <- Minimize(x)
  constr <- list(x * y == 4, y <= 2)
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 2.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 2.0, tolerance = 1e-3)
  expect_equal(as.numeric(value(y)), 2.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("scalar posynomial: x + 1/x", {
  x <- Variable(pos = TRUE)
  obj <- Minimize(x + 1/x)
  prob <- Problem(obj)
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 2.0, tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("exp/log atoms in DGP context", {
  x <- Variable(pos = TRUE)
  obj <- Minimize(exp(x))
  constr <- list(x >= 1)
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, exp(1), tolerance = 1e-4)
  expect_equal(as.numeric(value(x)), 1.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("sum of positive variables (vector GP)", {
  x <- Variable(3, pos = TRUE)
  obj <- Minimize(sum(x))
  constr <- list(x[1] * x[2] * x[3] >= 1)
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 3.0, tolerance = 1e-4)
  ## By AM-GM, optimal at x=1,1,1
  expect_equal(as.numeric(value(x)), rep(1, 3), tolerance = 1e-3)
})

## @cvxpy NONE
test_that("trace in DGP (matrix diagonal sum)", {
  A <- Variable(c(2, 2), pos = TRUE)
  obj <- Minimize(Trace(A))
  constr <- list(A[1,1] * A[2,2] >= 4)
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE)
  ## By AM-GM: trace = a11+a22 >= 2*sqrt(a11*a22) >= 2*2 = 4
  expect_equal(val, 4.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("error: non-DGP with gp=TRUE", {
  x <- Variable()
  prob <- Problem(Minimize(x^2), list(x >= -1))
  expect_error(psolve(prob, gp = TRUE), "not DGP compliant")
})

## @cvxpy NONE
test_that("error: DGP without gp=TRUE suggests gp=TRUE", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x + y), list(x * y <= 1))
  expect_error(psolve(prob), "gp = TRUE")
})

## @cvxpy NONE
test_that("Dgp2Dcp reduction_accepts", {
  dgp <- Dgp2Dcp()
  ## DGP problem → TRUE
  x <- Variable(pos = TRUE)
  prob_gp <- Problem(Minimize(x), list(x >= 1))
  expect_true(reduction_accepts(dgp, prob_gp))
  ## Non-DGP → FALSE
  y <- Variable()
  prob_dcp <- Problem(Minimize(y^2), list(y >= 0))
  expect_false(reduction_accepts(dgp, prob_dcp))
})

## @cvxpy NONE
test_that("multiple solvers produce consistent DGP results", {
  ## Test with SCS
  if (requireNamespace("scs", quietly = TRUE)) {
    x <- Variable(pos = TRUE)
    prob <- Problem(Minimize(x + 1/x))
    val_scs <- psolve(prob, solver = "SCS", gp = TRUE)
    expect_equal(val_scs, 2.0, tolerance = 1e-3)
  }
  ## Test with Clarabel (fresh objects to avoid cache leakage)
  if (requireNamespace("clarabel", quietly = TRUE)) {
    x2 <- Variable(pos = TRUE)
    prob2 <- Problem(Minimize(x2 + 1/x2))
    val_clar <- psolve(prob2, solver = "CLARABEL", gp = TRUE)
    expect_equal(val_clar, 2.0, tolerance = 1e-4)
  }
})

## @cvxpy NONE
test_that("constant canonicalizer rejects non-positive values", {
  ## Direct test of .dgp_constant_canon
  c_pos <- Constant(5)
  result <- CVXR:::.dgp_constant_canon(c_pos, list())
  expect_true(S7_inherits(result[[1L]], Constant))
  expect_equal(as.numeric(value(result[[1L]])), log(5), tolerance = 1e-10)

  ## Negative constant should error
  c_neg <- Constant(-1)
  expect_error(CVXR:::.dgp_constant_canon(c_neg, list()), "positive")
})

## @cvxpy NONE
test_that("GP with product of two terms (multiply test)", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  obj <- Minimize(x * y)
  constr <- list(x >= 2, y >= 3)
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 6.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("simple monomial: x", {
  x <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x), list(x >= 5))
  val <- psolve(prob, gp = TRUE)
  expect_equal(val, 5.0, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("GP infeasibility detected", {
  x <- Variable(pos = TRUE)
  prob <- Problem(Minimize(x), list(x >= 2, x <= 1))
  val <- psolve(prob, gp = TRUE)
  expect_true(is.na(val) || is.infinite(val))
})

## ── Phase DGP-2: Extended Canonicalizers ─────────────────────────
## CVXPY reference values verified via `uv run python`

## @cvxpy NONE
test_that("geo_mean in DGP", {
  x <- Variable(3, pos = TRUE)
  obj <- Minimize(1 / geo_mean(x))
  constr <- list(x[1] * x[2] * x[3] <= 1, x[1] >= 0.5)
  prob <- Problem(obj, constr)
  val <- psolve(prob, gp = TRUE)
  ## CVXPY: value=1.0, x=[0.887, 1.062, 1.062]
  expect_equal(val, 1.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("xexp in DGP", {
  y <- Variable(pos = TRUE)
  prob <- Problem(Minimize(xexp(y)), list(y >= 1))
  val <- psolve(prob, gp = TRUE)
  ## xexp(1) = 1 * exp(1) = e
  expect_equal(val, exp(1), tolerance = 1e-4)
  expect_equal(as.numeric(value(y)), 1.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("p_norm in DGP (p=2)", {
  z <- Variable(3, pos = TRUE)
  prob <- Problem(Minimize(p_norm(z, 2)),
                  list(z[1] >= 0.5, z[2] >= 0.3, z[3] >= 0.4))
  val <- psolve(prob, gp = TRUE)
  ## CVXPY: 0.70710678 (= sqrt(0.5))
  expect_equal(val, sqrt(0.5), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("p_norm in DGP (p=3)", {
  z <- Variable(2, pos = TRUE)
  prob <- Problem(Minimize(p_norm(z, 3)), list(z[1] >= 1, z[2] >= 2))
  val <- psolve(prob, gp = TRUE)
  ## CVXPY: 2.08008382 = (1^3 + 2^3)^(1/3)
  expect_equal(val, (1 + 8)^(1/3), tolerance = 1e-4)
})

## @cvxpy NONE
test_that("quad_over_lin in DGP (vector)", {
  a <- Variable(2, pos = TRUE)
  b <- Variable(pos = TRUE)
  prob <- Problem(Minimize(quad_over_lin(a, b)),
                  list(a[1] >= 0.5, a[2] >= 0.3, b <= 2, b >= 1))
  val <- psolve(prob, gp = TRUE)
  ## CVXPY: 0.17 = (0.25 + 0.09) / 2
  expect_equal(val, 0.17, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("quad_over_lin in DGP (scalar)", {
  a <- Variable(pos = TRUE)
  b <- Variable(pos = TRUE)
  prob <- Problem(Minimize(quad_over_lin(reshape_expr(a, c(1L, 1L)), b)),
                  list(a >= 2, b >= 1, b <= 4))
  val <- psolve(prob, gp = TRUE)
  ## CVXPY: 1.0 = 4/4
  expect_equal(val, 1.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("norm_inf in DGP", {
  w <- Variable(3, pos = TRUE)
  prob <- Problem(Minimize(norm_inf(w)),
                  list(w[1] >= 0.5, w[2] >= 0.3, w[3] >= 0.4))
  val <- psolve(prob, gp = TRUE)
  ## CVXPY: 0.5 (max of lower bounds)
  expect_equal(val, 0.5, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("norm_inf in DGP (2 vars)", {
  w <- Variable(2, pos = TRUE)
  prob <- Problem(Minimize(norm_inf(w)), list(w[1] >= 3, w[2] >= 1))
  val <- psolve(prob, gp = TRUE)
  ## CVXPY: 3.0
  expect_equal(val, 3.0, tolerance = 1e-4)
})

## @cvxpy NONE
test_that("quad_form in DGP (solver=SCS)", {
  ## Note: quad_form DGP requires explicit conic solver because the
  ## DGP-transformed problem has log_sum_exp (not QP-solvable).
  ## Note: CVXPY has a bug where quad_form canonicalizer is dead code
  ## (function vs class key mismatch). R implementation is correct.
  x <- Variable(2, pos = TRUE)
  P <- matrix(c(4, 1, 1, 9), 2, 2)
  prob <- Problem(Minimize(quad_form(x, P)),
                  list(x[1] >= 0.5, x[2] >= 0.3))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## Correct value: x^T P x at lower bounds = 0.5^2*4 + 2*0.5*0.3*1 + 0.3^2*9 = 2.11
  expect_equal(val, 2.11, tolerance = 1e-3)
  expect_equal(as.numeric(value(x)), c(0.5, 0.3), tolerance = 1e-2)
})

## @cvxpy NONE
test_that("DGP atoms: is_dgp property checks", {
  x <- Variable(3, pos = TRUE)
  y <- Variable(pos = TRUE)

  ## geo_mean is log-log concave → valid in Maximize objective
  expect_true(is_dgp(Problem(Maximize(geo_mean(x)))))

  ## p_norm is log-log convex → valid in Minimize objective
  expect_true(is_dgp(Problem(Minimize(p_norm(x, 2)))))

  ## xexp is log-log convex → valid in Minimize
  expect_true(is_dgp(Problem(Minimize(xexp(y)))))

  ## quad_form is log-log convex → valid in Minimize
  P <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  x2 <- Variable(2, pos = TRUE)
  expect_true(is_dgp(Problem(Minimize(quad_form(x2, P)))))

  ## quad_over_lin is log-log convex → valid in Minimize
  expect_true(is_dgp(Problem(Minimize(quad_over_lin(x, y)))))

  ## norm_inf is log-log convex → valid in Minimize
  expect_true(is_dgp(Problem(Minimize(norm_inf(x)))))
})

## ── Phase DGP-3: DGP-Specific Atoms ──────────────────────────────
## CVXPY reference values verified via `uv run python`

## @cvxpy NONE
test_that("prod_entries in DGP (minimize)", {
  x <- Variable(3, pos = TRUE)
  prob <- Problem(Minimize(prod_entries(x)),
                  list(x[1] >= 2, x[2] >= 3, x[3] >= 4))
  val <- psolve(prob, gp = TRUE)
  ## CVXPY: 24.0
  expect_equal(val, 24.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("prod(x) dispatch creates Prod in DGP", {
  x <- Variable(3, pos = TRUE)
  prob <- Problem(Maximize(prod(x)),
                  list(x[1] <= 2, x[2] <= 3, x[3] <= 4))
  val <- psolve(prob, gp = TRUE)
  ## CVXPY: 24.0
  expect_equal(val, 24.0, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("cumprod in DGP", {
  x <- Variable(3, pos = TRUE)
  prob <- Problem(Minimize(sum(cumprod(x))),
                  list(x[1] >= 2, x[2] >= 3, x[3] >= 4))
  val <- psolve(prob, gp = TRUE)
  ## CVXPY: cumprod=[2,6,24], sum=32
  expect_equal(val, 32.0, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("one_minus_pos in DGP", {
  y <- Variable(pos = TRUE)
  prob <- Problem(Maximize(one_minus_pos(y)),
                  list(y >= 0.1, y <= 0.5))
  val <- psolve(prob, gp = TRUE)
  ## CVXPY: 0.9 (at y=0.1)
  expect_equal(val, 0.9, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("diff_pos in DGP", {
  a <- Variable(pos = TRUE)
  b <- Variable(pos = TRUE)
  prob <- Problem(Maximize(diff_pos(a, b)),
                  list(a <= 5, b >= 0.5, b <= 1))
  val <- psolve(prob, gp = TRUE)
  ## CVXPY: 4.5 (a=5, b=0.5)
  expect_equal(val, 4.5, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("pf_eigenvalue in DGP", {
  X <- Variable(c(2, 2), pos = TRUE)
  prob <- Problem(Minimize(pf_eigenvalue(X)),
                  list(X[1,1] >= 0.1, X[2,2] >= 0.1,
                       X[1,2] <= 0.01, X[2,1] <= 0.01))
  val <- psolve(prob, gp = TRUE)
  ## CVXPY: ~0.1 (diagonal dominates)
  expect_equal(val, 0.1, tolerance = 1e-2)
})

## @cvxpy NONE
test_that("gmatmul in DGP", {
  A <- matrix(c(1, 3, 2, 4), 2, 2)  ## column-major: [[1,2],[3,4]]
  x <- Variable(2, pos = TRUE)
  prob <- Problem(Minimize(sum(gmatmul(A, x))),
                  list(x[1] >= 0.5, x[2] >= 0.5))
  val <- psolve(prob, gp = TRUE)
  ## CVXPY: 0.132813 (gmatmul=[[0.125, 0.0078125]])
  expect_equal(val, 0.1328125, tolerance = 1e-3)
})

## @cvxpy NONE
test_that("eye_minus_inv in DGP", {
  Y <- Variable(c(2, 2), pos = TRUE)
  prob <- Problem(Minimize(sum(eye_minus_inv(Y))),
                  list(Y <= 0.4))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## CVXPY: ~2.0 (Y→0 ⇒ (I-Y)^{-1} → I, sum(I_{2x2}) = 2)
  expect_equal(val, 2.0, tolerance = 0.05)
})

## @cvxpy NONE
test_that("resolvent in DGP", {
  Z <- Variable(c(2, 2), pos = TRUE)
  prob <- Problem(Minimize(sum(resolvent(Z, 2))),
                  list(Z <= 0.5))
  val <- psolve(prob, gp = TRUE, solver = "SCS")
  ## CVXPY: ~1.0 (Z→0 ⇒ (I - Z/s)^{-1}/s → I/s, sum(I/2) = 1)
  expect_equal(val, 1.0, tolerance = 0.05)
})

## @cvxpy NONE
test_that("DGP Phase 3 atoms: is_dgp checks", {
  x <- Variable(3, pos = TRUE)
  y <- Variable(pos = TRUE)
  X <- Variable(c(2, 2), pos = TRUE)

  ## prod is log-log affine → valid in both min and max
  expect_true(is_dgp(Problem(Minimize(prod_entries(x)))))
  expect_true(is_dgp(Problem(Maximize(prod_entries(x)))))

  ## cumprod is log-log affine
  expect_true(is_dgp(Problem(Minimize(sum(cumprod(x))))))

  ## one_minus_pos is log-log concave → valid in Maximize
  expect_true(is_dgp(Problem(Maximize(one_minus_pos(y)))))

  ## pf_eigenvalue is log-log convex → valid in Minimize
  expect_true(is_dgp(Problem(Minimize(pf_eigenvalue(X)))))

  ## gmatmul is log-log affine → valid in both
  A <- matrix(c(1,0,0,1), 2, 2)
  expect_true(is_dgp(Problem(Minimize(sum(gmatmul(A, X))))))
})
