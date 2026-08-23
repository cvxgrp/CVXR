## Regression guard: PSD constraint dual recovery must not be transposed.
##
## CVXR-specific bug (not a CVXPY parity item), found 2026-08-10 and fixed in
## 893ca4b.  The since-deleted per-solver expanders (`clarabel_triu_to_full` /
## `scs_tri_to_full`, now one shared `tri_to_full` in utilities/psd_utils.R --
## ADR D_19.6) transliterated numpy's
## `triu_indices`/`tril_indices` -- which enumerate ROW-major -- as R's
## `upper.tri()`/`lower.tri()`, which are COLUMN-major.  The two orders coincide
## only at n = 2, so every PSD dual with n >= 3 came back transposed while the
## objective stayed correct.  See cvxpy/utilities/psd_utils.py:52-61.
##
## ABSOLUTE ORACLE (deliberately not a self-referential comparison): for
##     minimize  tr(C X)   subject to   X >= I  (PSD)
## the Lagrangian is tr(CX) - tr(S(X - I)), so stationarity gives S = C and the
## optimum is tr(C).  With C symmetric positive definite the dual S = C is
## feasible and unique.  A transposition of the dual is therefore detectable
## whenever C is NOT symmetric under the scrambling -- which is why C must have
## distinct off-diagonal entries, and why n = 2 alone cannot catch this.

## Symmetric positive definite (diagonally dominant), distinct off-diagonals.
.psd_dual_C <- function(n) {
  stopifnot(n >= 2L)
  C <- matrix(0, n, n)
  diag(C) <- seq(2, length.out = n)
  for (i in seq_len(n - 1L)) {
    C[i, i + 1L] <- i / 2
    C[i + 1L, i] <- i / 2
  }
  C
}

.psd_dual_check <- function(solver, n, tol = 1e-5) {
  C <- .psd_dual_C(n)
  X <- Variable(c(n, n), symmetric = TRUE)
  con <- PSD(X - diag(n))
  val <- psolve(Problem(Minimize(sum(C * X)), list(con)), solver = solver)
  dual <- matrix(as.numeric(dual_value(con)), n, n)
  list(val = val, dual = dual, C = C)
}

## @cvxpy NONE -- CVXR-specific: R column-major vs numpy row-major triangle order
test_that("PSD dual is recovered untransposed (CLARABEL, n = 2..4)", {
  for (n in 2:4) {
    r <- .psd_dual_check("CLARABEL", n)
    expect_equal(r$val, sum(diag(r$C)), tolerance = 1e-6,
                 info = paste("objective, n =", n))
    expect_lt(max(abs(r$dual - r$C)), 1e-5)
    ## the dual must equal C, not its transpose-scramble: assert asymmetry of
    ## the failure mode is actually exercised
    expect_false(isTRUE(all.equal(r$C, t(r$C)[, rev(seq_len(n)), drop = FALSE])))
  }
})

## @cvxpy NONE -- CVXR-specific: R column-major vs numpy row-major triangle order
test_that("PSD dual is recovered untransposed (SCS, n = 2..4)", {
  skip_if_not_installed("scs")
  for (n in 2:4) {
    r <- .psd_dual_check("SCS", n, tol = 1e-4)
    expect_equal(r$val, sum(diag(r$C)), tolerance = 1e-4,
                 info = paste("objective, n =", n))
    expect_lt(max(abs(r$dual - r$C)), 1e-4)
  }
})

## @cvxpy NONE -- CVXR-specific: R column-major vs numpy row-major triangle order
test_that("PSD dual is recovered untransposed (MOSEK, n = 3)", {
  require_solver("MOSEK")
  r <- .psd_dual_check("MOSEK", 3L)
  expect_equal(r$val, sum(diag(r$C)), tolerance = 1e-5)
  expect_lt(max(abs(r$dual - r$C)), 1e-5)
})

## @cvxpy NONE -- CVXR-specific: pins the exact scrambling that shipped in 1.9.1
test_that("the n = 3 dual is not the transposed-triangle scramble", {
  ## Shipped 1.9.1 returned, for C = [[2,1,0],[1,3,1],[0,1,4]]:
  ##   [2, 1, 3/sqrt(2); 1, 0, 1; 3/sqrt(2), 1, 4]
  ## i.e. the (2,2) entry lost to (1,3).  Pin that it does not come back.
  n <- 3L
  C <- matrix(c(2, 1, 0, 1, 3, 1, 0, 1, 4), n, n, byrow = TRUE)
  X <- Variable(c(n, n), symmetric = TRUE)
  con <- PSD(X - diag(n))
  psolve(Problem(Minimize(sum(C * X)), list(con)), solver = "CLARABEL")
  dual <- matrix(as.numeric(dual_value(con)), n, n)
  expect_lt(max(abs(dual - C)), 1e-5)
  expect_gt(abs(dual[2L, 2L]), 2.5)          # was 0 under the bug
  expect_lt(abs(dual[1L, 3L]), 1e-5)         # was 3/sqrt(2) = 2.121 under the bug
})
