## Regression guard: PowCone3D duals must not be transposed across cones.
##
## CVXR-specific bug (not a CVXPY parity item), found 2026-08-10 and fixed in
## 1.9.1-1.  save_dual_value(PowCone3D) filled its (coord x cone) table
## byrow, but format_constraints stuffs PowCone3D INTERLEAVED per cone --
## (x1,y1,z1, x2,y2,z2, ...) -- and CVXR hands save_dual_value the raw flat
## vector (CVXPY pre-shapes it in ConeMatrixStuffing.invert; CVXR does not).
## The result was a transposed dual for every constraint with num_cones() >= 2.
##
## ABSOLUTE ORACLE (deliberately not a comparison against a reshaped twin of
## the same problem -- that is blind to a permutation applied to both halves).
## For
##     minimize  sum(x) + sum(y) - c'z   subject to   PowCone3D(x, y, z, alpha)
## each cone is independent and the KKT conditions force the per-cone dual to
## (1, 1, -c_i).  Distinct c_i make a transposition detectable; nc = 1 cannot
## catch it, because the two fill orders coincide there.

## @cvxpy NONE -- CVXR-specific: flat interleaved dual layout, not a CVXPY port
test_that("PowCone3D duals are per-coordinate, not transposed (nc = 2, the smallest failing case)", {
  alpha <- 0.4
  cc <- c(0.4, 0.9)
  n <- length(cc)
  x <- Variable(n); y <- Variable(n); z <- Variable(n)
  con <- PowCone3D(x, y, z, rep(alpha, n))
  psolve(Problem(Minimize(sum(x) + sum(y) - sum(cc * z)), list(con)),
         solver = "CLARABEL")
  dv <- lapply(con@dual_variables, function(d) as.numeric(value(d)))
  expect_equal(dv[[1L]], rep(1, n), tolerance = 1e-5)
  expect_equal(dv[[2L]], rep(1, n), tolerance = 1e-5)
  expect_equal(dv[[3L]], -cc,       tolerance = 1e-5)
})

## @cvxpy NONE -- CVXR-specific: flat interleaved dual layout, not a CVXPY port
test_that("PowCone3D duals are per-coordinate for nc = 3 and nc = 6", {
  alpha <- 0.4
  for (cc in list(c(0.6, 0.75, 0.9), c(0.5, 0.55, 0.6, 0.7, 0.8, 0.95))) {
    n <- length(cc)
    x <- Variable(n); y <- Variable(n); z <- Variable(n)
    con <- PowCone3D(x, y, z, rep(alpha, n))
    psolve(Problem(Minimize(sum(x) + sum(y) - sum(cc * z)), list(con)),
           solver = "CLARABEL")
    dv <- lapply(con@dual_variables, function(d) as.numeric(value(d)))
    expect_equal(dv[[1L]], rep(1, n), tolerance = 1e-5, info = paste("nc =", n))
    expect_equal(dv[[2L]], rep(1, n), tolerance = 1e-5, info = paste("nc =", n))
    expect_equal(dv[[3L]], -cc,       tolerance = 1e-5, info = paste("nc =", n))
  }
})

## @cvxpy NONE -- CVXR-specific: matrix-shaped args, same interleaved layout
test_that("PowCone3D duals are correct for matrix-shaped arguments", {
  alpha <- 0.4; m <- 2L; k <- 3L
  cm <- matrix(c(0.55, 0.6, 0.7, 0.8, 0.9, 0.95), m, k)
  x <- Variable(c(m, k)); y <- Variable(c(m, k)); z <- Variable(c(m, k))
  con <- PowCone3D(x, y, z, alpha)
  psolve(Problem(Minimize(sum(x) + sum(y) - sum(cm * z)), list(con)),
         solver = "CLARABEL")
  dv <- lapply(con@dual_variables, value)
  expect_equal(as.numeric(dv[[1L]]), rep(1, m * k), tolerance = 1e-5)
  expect_equal(as.numeric(dv[[2L]]), rep(1, m * k), tolerance = 1e-5)
  ## column-major (F-order) flattening, matching how the cone is stuffed
  expect_equal(as.numeric(dv[[3L]]), -as.numeric(cm), tolerance = 1e-5)
  expect_equal(dim(dv[[3L]]), c(m, k))
})

## @cvxpy NONE -- CVXR-specific: nc = 1 must keep working (both fills coincide)
test_that("PowCone3D scalar cone (nc = 1) is unaffected", {
  alpha <- 0.4; cval <- 0.4
  x <- Variable(1); y <- Variable(1); z <- Variable(1)
  con <- PowCone3D(x, y, z, alpha)
  psolve(Problem(Minimize(x + y - cval * z), list(con)), solver = "CLARABEL")
  dv <- lapply(con@dual_variables, function(d) as.numeric(value(d)))
  expect_equal(dv[[1L]], 1,      tolerance = 1e-5)
  expect_equal(dv[[2L]], 1,      tolerance = 1e-5)
  expect_equal(dv[[3L]], -cval,  tolerance = 1e-5)
})
