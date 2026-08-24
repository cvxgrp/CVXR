## @cvxpy cvxpy/utilities/coeff_extractor.py:227-262
##
## A scalar multiplier on a VECTOR-valued square() was silently discarded on the
## QP path, giving the wrong minimizer with no error and no warning.
##
## `power_canon` builds SymbolicQuadForm(x, Diagonal(n), square(x)).  square()
## is elementwise, so the placeholder is a VECTOR and each component carries its
## own affine coefficient.  The extractor took a single `scalar_coeff` and, for a
## non-scalar placeholder, fell back to 1.0 -- discarding every coefficient.
## CVXPY's DIAGONAL path (coeff_extractor.py:245-262) builds the contribution
## from c_part per component instead.
##
## WHY THIS WENT UNNOTICED, and how these tests are built to catch it:
##
##   A PURE quadratic has the SAME MINIMIZER under any positive scale, so
##   `2*sum(square(x))` and `sum(square(x))` agree on x, and the reported
##   objective is recomputed from the original expression and looks right too.
##   Only a QUADRATIC PLUS A LINEAR TERM makes the scale change the answer:
##
##       minimize  a*sum(x^2) + b*sum(x)   =>   x_i = -b/(2a)
##
##   Every block below therefore carries a linear term.  Asserting on the
##   objective value alone is not enough either -- assert on x.
##
## Unaffected by the defect (they take the SCALAR path) and pinned here as
## controls: sum_squares(), quad_form(), p_norm(x,2)^2.

quad_solvers <- c("OSQP", "CLARABEL")

## minimize a*sum(x^2) - 4*sum(x) over a box -> x_i = 2/a
fit <- function(build, solver, n = 2L) {
  x <- Variable(n)
  prob <- Problem(Minimize(build(x)), list(x >= -10, x <= 10))
  val <- psolve(prob, solver = solver)
  list(x = as.numeric(value(x)), obj = drop(val))
}
expected <- function(a, n = 2L) {
  xs <- 2 / a
  list(x = rep(xs, n), obj = n * (a * xs^2 - 4 * xs))
}

test_that("a scalar multiplier on sum(square(x)) is not discarded", {
  for (s in quad_solvers) {
    for (a in c(1, 2, 3, 0.5)) {
      got  <- fit(function(x) a * sum_entries(square(x)) - 4 * sum_entries(x), s)
      want <- expected(a)
      expect_equal(got$x,   want$x,   tolerance = 1e-4,
                   info = paste(s, "a =", a))
      expect_equal(got$obj, want$obj, tolerance = 1e-4,
                   info = paste(s, "a =", a))
    }
  }
})

test_that("the multiplier is honored wherever it sits in the expression", {
  for (s in quad_solvers) {
    want <- expected(2)
    ## Outside the sum, inside the sum, on the right, and via `^`.
    expect_equal(fit(function(x) 2 * sum_entries(square(x)) - 4 * sum_entries(x), s)$x,
                 want$x, tolerance = 1e-4, info = paste(s, "outside"))
    expect_equal(fit(function(x) sum_entries(2 * square(x)) - 4 * sum_entries(x), s)$x,
                 want$x, tolerance = 1e-4, info = paste(s, "inside"))
    expect_equal(fit(function(x) sum_entries(square(x)) * 2 - 4 * sum_entries(x), s)$x,
                 want$x, tolerance = 1e-4, info = paste(s, "right"))
    expect_equal(fit(function(x) 2 * sum_entries(x^2) - 4 * sum_entries(x), s)$x,
                 want$x, tolerance = 1e-4, info = paste(s, "caret"))
  }
})

test_that("per-component weights are kept distinct, not collapsed", {
  ## The sharpest form: the coefficient VECTOR differs per component, so a
  ## single scalar cannot represent it however it is chosen.
  ## minimize sum_i w_i x_i^2 - 4 x_i  ->  x_i = 2/w_i
  w <- c(1, 4)
  for (s in quad_solvers) {
    x <- Variable(2)
    prob <- Problem(Minimize(sum_entries(w * square(x)) - 4 * sum_entries(x)),
                    list(x >= -10, x <= 10))
    psolve(prob, solver = s)
    expect_equal(as.numeric(value(x)), 2 / w, tolerance = 1e-4, info = s)
  }
})

test_that("ridge regression honors lambda", {
  ## The idiom that made this matter: lambda * sum(square(w)).  Before the fix
  ## every lambda produced the SAME fit.  Checked against the closed form
  ## w = (A'A + lambda I)^-1 A'b.
  set.seed(1)
  A <- matrix(rnorm(30), 15, 2)
  b <- rnorm(15)
  for (lam in c(0.1, 1, 10)) {
    w  <- Variable(2)
    prob <- Problem(Minimize(sum_entries(square(A %*% w - b)) +
                             lam * sum_entries(square(w))))
    psolve(prob, solver = "CLARABEL")
    want <- solve(t(A) %*% A + lam * diag(2), t(A) %*% b)
    expect_equal(as.numeric(value(w)), as.numeric(want), tolerance = 1e-4,
                 info = paste("lambda =", lam))
  }
})

test_that("different lambdas give different ridge fits", {
  ## A guard on the guard: if the penalty were dropped again, every lambda
  ## would coincide and the block above could only fail by tolerance.
  set.seed(2)
  A <- matrix(rnorm(30), 15, 2)
  b <- rnorm(15)
  fits <- lapply(c(0.1, 10), function(lam) {
    w <- Variable(2)
    psolve(Problem(Minimize(sum_entries(square(A %*% w - b)) +
                            lam * sum_entries(square(w)))), solver = "CLARABEL")
    as.numeric(value(w))
  })
  expect_gt(max(abs(fits[[1]] - fits[[2]])), 1e-3)
})

test_that("the scalar path is unchanged — controls", {
  ## These take the SCALAR branch and were always correct; pinned so the split
  ## introduced with the fix cannot regress them.
  for (s in quad_solvers) {
    want <- expected(2)
    expect_equal(fit(function(x) 2 * sum_squares(x) - 4 * sum_entries(x), s)$x,
                 want$x, tolerance = 1e-4, info = paste(s, "sum_squares"))
    expect_equal(fit(function(x) 2 * quad_form(x, diag(2)) - 4 * sum_entries(x), s)$x,
                 want$x, tolerance = 1e-4, info = paste(s, "quad_form"))
  }
  ## CLARABEL only: p_norm(x, 2)^2 canonicalizes through an SOC, and OSQP is a
  ## QP-only solver that rejects SOC cones. Running it on both solvers is a bug
  ## in the TEST, not in CVXR -- it cost one spurious ERROR in the pre-fix run.
  expect_equal(fit(function(x) 2 * p_norm(x, 2)^2 - 4 * sum_entries(x), "CLARABEL")$x,
               expected(2)$x, tolerance = 1e-4, info = "CLARABEL p_norm^2")
})

test_that("a non-diagonal P still works through quad_form", {
  ## quad_form with a genuinely non-diagonal P must keep working: it is scalar,
  ## so it takes the scalar path and never meets the diagonal-P guard.
  P <- matrix(c(2, 1, 1, 2), 2, 2)
  for (s in quad_solvers) {
    x <- Variable(2)
    prob <- Problem(Minimize(quad_form(x, P) - 4 * sum_entries(x)),
                    list(x >= -10, x <= 10))
    psolve(prob, solver = s)
    ## minimize x'Px - 4'x  ->  2Px = 4*1  ->  x = 2 * P^-1 1
    expect_equal(as.numeric(value(x)),
                 as.numeric(2 * solve(P, rep(1, 2))), tolerance = 1e-4, info = s)
  }
})

test_that("the unscaled case still solves — no regression at coefficient 1", {
  for (s in quad_solvers) {
    want <- expected(1)
    got <- fit(function(x) sum_entries(square(x)) - 4 * sum_entries(x), s)
    expect_equal(got$x,   want$x,   tolerance = 1e-4, info = s)
    expect_equal(got$obj, want$obj, tolerance = 1e-4, info = s)
  }
})
