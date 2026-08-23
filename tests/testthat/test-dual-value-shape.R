## Dual-value SHAPES, per constraint kind (ADR D_19.6 step C).
##
## `ConeMatrixStuffing.invert` gives each dual the shape of its constraint, in
## F (column-major) order, EXCEPT for an exemption list transliterated from
## CVXPY (cone_matrix_stuffing.py:473-478): scalars, ExpCone, SOC, and
## multi-cone PSD. Both halves matter, so both are asserted here -- a future
## "simplification" that drops an exemption has to fail something.
##
## These are CVXR-side shape guarantees with no single upstream test to map to;
## the behavior they pin is CVXPY's, cited at each block.
##
## The VALUES (not just shapes) are checked against absolute oracles, never
## against a reshaped twin of the same problem -- see the methodological warning
## in test-psd-dual-recovery.R.

skip_if_not_installed("clarabel")

# ---- reshaped: the constraint's own shape --------------------------

## @cvxpy NONE
test_that("NonNeg and Zero duals take the constraint's matrix shape", {
  X <- Variable(c(2L, 3L))
  p <- Problem(Minimize(sum_entries(X)), list(X >= 1))
  psolve(p, solver = "CLARABEL")
  d <- dual_value(p@constraints[[1L]])
  expect_equal(dim(d), c(2L, 3L))
  ## min sum(X) s.t. X >= 1: every entry is active with multiplier 1.
  expect_equal(as.numeric(d), rep(1, 6L), tolerance = 1e-6)

  Y <- Variable(c(2L, 3L))
  q <- Problem(Minimize(sum_entries(square(Y))), list(Y == 2))
  psolve(q, solver = "CLARABEL")
  dq <- dual_value(q@constraints[[1L]])
  expect_equal(dim(dq), c(2L, 3L))
  ## min sum(Y^2) s.t. Y = 2: stationarity gives nu = 2*Y = 4 per entry
  ## (sign convention aside, the magnitude is what pins the shape mapping).
  expect_equal(abs(as.numeric(dq)), rep(4, 6L), tolerance = 1e-5)
})

## @cvxpy NONE
test_that("a vector constraint's dual is an (n, 1) column, not a bare vector", {
  z <- Variable(3L)
  p <- Problem(Minimize(sum_entries(z)), list(z >= 1))
  psolve(p, solver = "CLARABEL")
  d <- dual_value(p@constraints[[1L]])
  expect_equal(dim(d), c(3L, 1L))
})

# ---- exempt: shape () -- CVXPY's `shape == ()` clause ---------------

## @cvxpy NONE
test_that("a scalar constraint's dual stays a scalar", {
  t <- Variable()
  p <- Problem(Minimize(t), list(t >= 1))
  psolve(p, solver = "CLARABEL")
  d <- dual_value(p@constraints[[1L]])
  expect_null(dim(d))
  expect_length(d, 1L)
  expect_equal(as.numeric(d), 1, tolerance = 1e-6)
})

# ---- PSD: (n, n), composing with the triangle fix (893ca4b) ---------

## @cvxpy NONE
test_that("PSD duals are (n, n) and still equal C for n = 3 and n = 4", {
  for (n in c(3L, 4L)) {
    C <- matrix(0, n, n)
    diag(C) <- seq(2, length.out = n)
    ## Distinct off-diagonals: a transposition must be detectable.
    for (i in seq_len(n - 1L)) {
      C[i, i + 1L] <- C[i + 1L, i] <- 0.1 * i
    }
    X <- Variable(c(n, n), symmetric = TRUE)
    p <- Problem(Minimize(matrix_trace(C %*% X)), list((X - diag(n)) %>>% 0))
    val <- psolve(p, solver = "CLARABEL")
    d <- dual_value(p@constraints[[1L]])
    expect_equal(dim(d), c(n, n))
    ## Analytic: optimum is tr(C) and the unique dual is C itself.
    expect_equal(val, sum(diag(C)), tolerance = 1e-5)
    expect_equal(as.matrix(d), C, tolerance = 1e-5)
  }
})

# ---- exempt: SOC and ExpCone keep their per-argument dual lists -----

## @cvxpy NONE
test_that("SOC duals stay a per-argument list (exempt from the reshape)", {
  u <- Variable(3L)
  s <- Variable()
  p <- Problem(Minimize(s), list(SOC(s, u), u >= c(1, 2, 2)))
  psolve(p, solver = "CLARABEL")
  d <- dual_value(p@constraints[[1L]])
  expect_type(d, "list")
  expect_length(d, 2L)
  ## Analytic: at the optimum ||u|| = 3 and the SOC dual is (1, -u/||u||).
  expect_equal(as.numeric(d[[1L]]), 1, tolerance = 1e-5)
  expect_equal(as.numeric(d[[2L]]), -c(1, 2, 2) / 3, tolerance = 1e-5)
})

## @cvxpy NONE
test_that("ExpCone duals stay a per-argument list (exempt from the reshape)", {
  xv <- Variable()
  yv <- Variable()
  zv <- Variable()
  p <- Problem(Minimize(zv), list(ExpCone(xv, yv, zv), xv >= 1, yv == 1))
  psolve(p, solver = "CLARABEL")
  d <- dual_value(p@constraints[[1L]])
  expect_type(d, "list")
  expect_length(d, 3L)
  for (part in d) expect_length(as.numeric(part), 1L)
})

# ---- PowCone3D: reshaped to (3, num_cones), composing with 2c16fca --

## @cvxpy NONE
test_that("PowCone3D duals survive the reshape for num_cones() >= 2", {
  nc <- 2L
  cvec <- c(0.7, 1.3)
  xv <- Variable(nc)
  yv <- Variable(nc)
  zv <- Variable(nc)
  p <- Problem(
    Minimize(sum_entries(xv) + sum_entries(yv)),
    list(PowCone3D(xv, yv, zv, rep(0.5, nc)), zv == cvec, xv >= 0, yv >= 0)
  )
  psolve(p, solver = "CLARABEL")
  d <- dual_value(p@constraints[[1L]])
  expect_type(d, "list")
  expect_length(d, 3L)
  ## Each cone contributes its own (x, y, z) dual; the z-part carries the sign
  ## that a transposed (coord x cone) table would scramble.
  for (part in d) expect_length(as.numeric(part), nc)
})
