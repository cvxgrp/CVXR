## @cvxpy cvxpy/atoms/affine/wraps.py:107-130
## @cvxpy cvxpy/reductions/complex2real/canonicalizers/matrix_canon.py:197-200
## @cvxpy cvxpy/reductions/cvx_attr2constr.py:151
##
## `nsd_wrap` is the seventh Wrap subclass; CVXR had six. Its two upstream
## consumers had been handled differently, one soundly and one not:
##
##   cvx_attr2constr.py:151   CVXR used `-psd_wrap(-expr)` -- curvature- and
##                            sign-equivalent, but two extra nodes and not the
##                            upstream tree.
##   matrix_canon.py:200      CVXR had NO substitute. It applied psd_wrap
##                            UNCONDITIONALLY, with neither the is_psd() guard
##                            nor the NSD arm, so for a complex quad_form with
##                            an NSD P it asserted the real embedding was PSD.
##                            is_dcp() then said TRUE and canonicalization
##                            aborted with "must be a minimization DCP problem"
##                            -- CVXPY solves the same problem.
##
## Ledger + paired CVXPY reproducer:
##   notes/audit/completeness_ledger.md findings 7 and 8
##   notes/audit/repro/04_complex_quad_form_nsd.{py,R}

library(testthat)
library(CVXR)

test_that("nsd_wrap exists and asserts NSD", {
  A <- Constant(matrix(c(1, 0, 0, 1), 2, 2))
  w <- nsd_wrap(A)
  expect_s3_class(w, "CVXR::nsd_wrap")
  expect_s3_class(w, "CVXR::Wrap")
  expect_true(is_nsd(w))
  expect_false(is_psd(w))
  expect_true(is_hermitian(w))
  expect_true(is_symmetric(w))          # real argument
  expect_equal(w@shape, c(2L, 2L))

  ## Mirror image of psd_wrap, so the two cannot drift apart.
  p <- psd_wrap(A)
  expect_true(is_psd(p))
  expect_false(is_nsd(p))
})

test_that("nsd_wrap is a no-op numerically and requires a square matrix", {
  M <- matrix(c(-1, 0, 0, -2), 2, 2)
  expect_equal(value(nsd_wrap(Constant(M))), M)
  expect_error(nsd_wrap(Constant(matrix(1:6, 2, 3))), "square matrix")
  expect_error(psd_wrap(Constant(matrix(1:6, 2, 3))), "square matrix")
})

test_that("nsd_wrap has an NLP diff-engine converter, like its siblings", {
  ## registry.R registers the pass-through converters per CLASS, not on the
  ## Wrap base (mirroring CVXPY's ATOM_CONVERTERS, which lists "nsd_wrap"
  ## explicitly at registry.py:249). A new Wrap subclass without its own line
  ## would abort the diff engine on any DNLP problem containing it.
  conv <- diff_engine_convert
  for (cls in list(nsd_wrap, psd_wrap, nonneg_wrap, nonpos_wrap,
                   hermitian_wrap, symmetric_wrap, skew_symmetric_wrap)) {
    expect_true(tryCatch({ S7::method(conv, cls); TRUE },
                         error = function(e) FALSE))
  }
})

test_that("nsd_wrap of a complex argument is Hermitian but not symmetric", {
  ## CVXPY wraps.py:124-125 -- "symmetric" means REAL symmetric.
  Z <- Variable(c(2, 2), hermitian = TRUE)
  w <- nsd_wrap(Z)
  expect_true(is_hermitian(w))
  expect_false(is_symmetric(w))
})

test_that("an NSD variable attribute lowers through nsd_wrap", {
  ## cvx_attr2constr.py:151. Previously `-psd_wrap(-expr)`; the observable
  ## properties must be unchanged by the switch.
  X <- Variable(c(2, 2), NSD = TRUE)
  prob <- Problem(Maximize(matrix_trace(X)), list(X[1, 1] >= -3))
  val <- psolve(prob, solver = "SCS")
  expect_equal(status(prob), "optimal")
  expect_true(val <= 1e-6)                    # trace of an NSD matrix is <= 0
  ev <- eigen(value(X), symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev <= 1e-6))                # solution really is NSD
})

test_that("complex quad_form with an NSD P is maximizable, as in CVXPY", {
  ## THE finding. Three cells, so the divergence is attributed rather than
  ## merely observed: only the complex+NSD cell exercises the missing branch.
  ## Measured against cvxpy 1.9.2: all three are optimal, at -1, 1 and -1.
  n <- 2L
  Ppsd <- diag(n)
  Pnsd <- -diag(n)

  ## A: real + NSD + Maximize -- never enters complex2real.
  xr <- Variable(n)
  pa <- Problem(Maximize(quad_form(xr, Pnsd, assume_PSD = FALSE)),
                list(xr[1] >= 1))
  expect_equal(psolve(pa, solver = "SCS"), -1, tolerance = 1e-4)

  ## B: complex + PSD + Minimize -- enters it, and psd_wrap is the RIGHT wrap.
  xb <- Variable(n, complex = TRUE)
  pb <- Problem(Minimize(Re(quad_form(xb, Ppsd, assume_PSD = FALSE))),
                list(Re(xb)[1] >= 1))
  expect_equal(psolve(pb, solver = "SCS"), 1, tolerance = 1e-4)

  ## C: complex + NSD + Maximize -- the cell that used to abort.
  xc <- Variable(n, complex = TRUE)
  pc <- Problem(Maximize(Re(quad_form(xc, Pnsd, assume_PSD = FALSE))),
                list(Re(xc)[1] >= 1))
  expect_true(is_dcp(pc))
  expect_equal(psolve(pc, solver = "SCS"), -1, tolerance = 1e-4)
  expect_equal(status(pc), "optimal")
})

test_that("complex quad_form with an INDEFINITE P is not wrapped", {
  ## Upstream's guard has no `else`: an indefinite P gets NO wrap, so nothing
  ## claims it is PSD and the failure surfaces downstream on its own terms.
  ## Before the fix CVXR asserted PSD here too.
  Pind <- matrix(c(1, 0, 0, -1), 2, 2)
  expect_false(is_psd(Constant(Pind)))
  expect_false(is_nsd(Constant(Pind)))

  x <- Variable(2, complex = TRUE)
  q <- quad_form(x, Pind, assume_PSD = FALSE)
  ## An indefinite quadratic form is neither convex nor concave, so neither
  ## direction may be accepted as DCP.
  expect_false(is_dcp(Problem(Minimize(Re(q)), list(Re(x)[1] >= 1))))
  expect_false(is_dcp(Problem(Maximize(Re(q)), list(Re(x)[1] >= 1))))
})
