## Differentiation regression NET for the dict-diff-chain refactor (ADR D_19.5).
##
## Purpose: a before/after oracle + protective net for CVXPY PR #3147 part (A)
## (dict-in/dict-out diff hooks). Part (A) has LANDED: the complex-Variable and
## symmetric-Variable derivative cases below are now active and pass. It (1)
## protects the working real path (silent-zeroing / mis-length / recycling would
## otherwise pass unnoticed) and (2) keeps the one remaining dormant case
## (complex-PARAMETER derivative, blocked by the complex-sparse boundary, NOT
## #3147) pinned as a skip()-ed `@v19-pending` oracle.
##
## HARD INVARIANT (D_19.5 addendum): R recycling in the diff chain is a
## silent-correctness disaster (exact-multiple lengths warn NOTHING). Every
## backward()/derivative() here runs under `.no_recycle()`, which escalates any
## recycling-class warning to a hard failure. The net can NEVER go green over one.

skip_if_not_installed("diffcp")
skip_if_not_installed("clarabel")

## Escalate recycling-class warnings (and the replacement-length variant) to
## errors so a length mismatch in the diff chain fails loudly instead of
## silently returning recycled garbage.
.no_recycle <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("longer object length|number of items to replace|recycl",
                msg, ignore.case = TRUE)) {
        stop("RECYCLING in diff chain (forbidden, D_19.5): ", msg, call. = FALSE)
      }
    }
  )
}

# ===================================================================
# ACTIVE — the protective net (must stay green THROUGH the refactor)
# ===================================================================

## Real vector, analytic: min ||x - p||^2  =>  x* = p, dx*/dp = I.
##   forward:  delta(x) = delta(p)
##   backward: gradient(p) = (dx/dp)^T gradient(x) = gradient(x)
test_that("net[real-vec]: forward delta and backward grad match dx/dp = I", {
  p <- Parameter(2); value(p) <- c(1, 2)
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - p)))
  .no_recycle({
    psolve(prob, requires_grad = TRUE)
    delta(p) <- c(0.3, -0.4)
    derivative(prob)
    expect_equal(as.numeric(delta(x)), c(0.3, -0.4), tolerance = 1e-5)

    gradient(x) <- c(1, 1)
    backward(prob)
    expect_equal(as.numeric(gradient(p)), c(1, 1), tolerance = 1e-5)
  })
})

## Two parameters, analytic: min ||x - a||^2 + ||x - b||^2  =>  x* = (a+b)/2,
##   dx*/da = dx*/db = 0.5 I.  backward grad(a) = 0.5 * gradient(x).
test_that("net[multi-param]: backward grad splits 0.5/0.5 across two params", {
  a <- Parameter(2); b <- Parameter(2)
  value(a) <- c(1, 1); value(b) <- c(3, 3)
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x - a) + sum_squares(x - b)))
  .no_recycle({
    psolve(prob, requires_grad = TRUE)
    gradient(x) <- c(1, 1)
    backward(prob)
    expect_equal(as.numeric(gradient(a)), c(0.5, 0.5), tolerance = 1e-5)
    expect_equal(as.numeric(gradient(b)), c(0.5, 0.5), tolerance = 1e-5)
  })
})

## SOC, structural guard: parametrized second-order cone. We don't pin an
## analytic value (SOC couples the coordinates) -- the point is that the
## derivative path produces a non-NULL, correctly-shaped, finite delta with NO
## recycling. This is exactly the silent-breakage a chain refactor risks.
test_that("net[SOC]: derivative produces a finite, correctly-shaped delta", {
  p <- Parameter(2); value(p) <- c(1, 2)
  x <- Variable(2); t <- Variable()
  prob <- Problem(Minimize(t), list(SOC(t, x - p), x >= 0.5))
  .no_recycle({
    psolve(prob, requires_grad = TRUE)
    delta(p) <- c(1e-3, 0)
    derivative(prob)
    dv <- delta(x)
    expect_false(is.null(dv))
    expect_length(as.numeric(dv), 2L)
    expect_true(all(is.finite(as.numeric(dv))))
  })
})

# ===================================================================
# POST-(A) — complex-Variable + symmetric-Variable derivatives (now ACTIVE).
# Part (A) wired the Complex2Real var hooks (real/imag split) and the
# symmetric-Variable upper-tri<->full adjoint reshape, so the complex-VARIABLE
# and symmetric-VARIABLE cases below pass. The single remaining dormant case is
# the complex-PARAMETER derivative (first test), blocked by the complex-sparse
# boundary (EvalParams), which is a separate sweep item, NOT #3147.
# ===================================================================

## @cvxpy test_derivative.py::TestBackwardComplex::test_derivative_complex_param_complex_variable
## @v19-pending: complex-PARAMETER derivative [sweep] -- distinct from the dict-API
## fix. A complex Parameter is forced through EvalParams (CVXR has no complex
## sparse; it substitutes the complex param to a constant), so no param_prog is
## cached and derivative() cannot differentiate w.r.t. it. Unblocking needs the
## complex param to flow as a real/imag-split FREE param through param_prog
## (the "carry complex sparse as a real/imag pair" idea) -- tied to the
## complex-sparse boundary, NOT #3147. The complex-VARIABLE cases below DO work.
test_that("net[complex param+var]: dv = dq  (v* = q)", {
  skip("@v19-pending: complex-parameter derivative blocked by EvalParams (no complex sparse); separate from #3147-A")
  q <- Parameter(complex = TRUE); v <- Variable(complex = TRUE)
  prob <- Problem(Minimize(sum_squares(Re(v) - Re(q)) + sum_squares(Im(v) - Im(q))))
  value(q) <- 2 + 3i
  .no_recycle({
    psolve(prob, requires_grad = TRUE)
    delta(q) <- 1 + 0.5i
    derivative(prob)
    expect_equal(as.numeric(Re(delta(v))), 1.0, tolerance = 1e-2)
    expect_equal(as.numeric(Im(delta(v))), 0.5, tolerance = 1e-2)
  })
})

## @cvxpy test_derivative.py::TestBackwardComplex::test_derivative_complex_variable_imag_param
test_that("net[complex var, imag param]: perturb b => dz = 1i", {
  a <- Parameter(); b <- Parameter(); z <- Variable(complex = TRUE)
  prob <- Problem(Minimize(sum_squares(Re(z) - a) + sum_squares(Im(z) - b)))
  value(a) <- 3; value(b) <- 4
  .no_recycle({
    psolve(prob, requires_grad = TRUE)
    delta(a) <- 0; delta(b) <- 1
    derivative(prob)
    expect_equal(as.numeric(Im(delta(z))), 1.0, tolerance = 1e-2)
    expect_equal(as.numeric(Re(delta(z))), 0.0, tolerance = 1e-2)
  })
})

## @cvxpy test_derivative.py::TestBackwardComplex::test_derivative_complex_variable_purely_imaginary
test_that("net[purely imaginary]: perturb p => dw = 1i", {
  p <- Parameter(); w <- Variable(complex = TRUE)
  prob <- Problem(Minimize(sum_squares(Re(w)) + sum_squares(Im(w) - p)))
  value(p) <- 5
  .no_recycle({
    psolve(prob, requires_grad = TRUE)
    delta(p) <- 1
    derivative(prob)
    expect_equal(as.numeric(Im(delta(w))), 1.0, tolerance = 1e-2)
  })
})

## CVXR regression guard for the D_19.1 "concrete correctness follow-up": a
## symmetric/PSD Variable's gradient is reshaped upper-tri<->full in the
## diff chain (CvxAttr2Constr var_forward/var_backward). Analytic:
## min ||X - P||^2 over symmetric X with symmetric P  =>  X* = P, dX/dP = I.
## (Was @v19-pending pre-#3147-A; now active.)
test_that("net[symmetric var]: dX = dP  (X* = P), no recycling", {
  P <- Parameter(c(2, 2), symmetric = TRUE)
  value(P) <- matrix(c(1, 0.5, 0.5, 2), 2, 2)
  X <- Variable(c(2, 2), symmetric = TRUE)
  prob <- Problem(Minimize(sum_squares(X - P)))
  dP <- matrix(c(0.1, 0.05, 0.05, -0.2), 2, 2)
  .no_recycle({
    psolve(prob, requires_grad = TRUE)
    delta(P) <- dP
    derivative(prob)
    expect_equal(as.matrix(delta(X)), dP, tolerance = 1e-4)
  })
})
