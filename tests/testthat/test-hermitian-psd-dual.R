## Dual variable of a PSD constraint on a Hermitian matrix.
##
## CVXPY 1.9.2 (complex2real.py:340-349) added a factor of 2 that CVXR was
## missing: the recovered dual was exactly HALF its correct value, silently,
## for every complex PSD constraint.
##
## A Hermitian X >> 0 is solved through the real dilation
##     Y = [[re(X), -im(X)], [im(X), re(X)]] >> 0
## and the dilation doubles the inner product,
##     <D, Y> = 2*Re(<D[1:n,1:n] + 1i*D[(n+1):2n,1:n], X>),
## so the recovery must scale by 2.
##
## These tests do NOT assert the buggy value with a fudge factor -- they check
## the two things that must hold independent of implementation: the KKT
## stationarity condition and strong duality.

## @cvxpy test_complex.py::TestComplex::test_hermitian_psd_dual
test_that("Hermitian PSD dual satisfies KKT and strong duality", {
  ## R fills matrices column-major, so this is [[2, 1-1i], [1+1i, 3]].
  C <- matrix(c(2 + 0i, 1 + 1i, 1 - 1i, 3 + 0i), 2, 2)
  A <- matrix(c(1 + 0i, -0.5i, 0.5i, 1 + 0i), 2, 2)
  X <- Variable(c(2, 2), hermitian = TRUE)
  con <- X %>>% A
  prob <- Problem(Minimize(Re(matrix_trace(C %*% X))), list(con))
  val <- psolve(prob, solver = "CLARABEL")
  expect_equal(val, 4, tolerance = 1e-5)

  ## KKT gives the dual variable Y == C exactly.
  ## Before the fix this returned C/2 = [[1, 0.5+0.5i], [0.5-0.5i, 1.5]].
  Y <- dual_value(con)
  expect_equal(as.vector(Y), as.vector(C), tolerance = 1e-5)

  ## Strong duality: the optimal value equals Re(<Y, A>).
  expect_equal(Re(sum(Conj(Y) * A)), val, tolerance = 1e-5)
})

## @cvxpy NONE -- guards the scaling in the opposite direction.
## A REAL symmetric PSD constraint never goes through the dilation, so its
## dual must be untouched by the factor of 2.
test_that("real symmetric PSD duals are unaffected by the dilation scaling", {
  Cr <- matrix(c(2, 1, 1, 3), 2, 2)
  Ar <- matrix(c(1, 0, 0, 1), 2, 2)
  Xr <- Variable(c(2, 2), symmetric = TRUE)
  conr <- Xr %>>% Ar
  probr <- Problem(Minimize(matrix_trace(Cr %*% Xr)), list(conr))
  valr <- psolve(probr, solver = "CLARABEL")

  Yr <- dual_value(conr)
  expect_equal(as.vector(as.matrix(Yr)), as.vector(Cr), tolerance = 1e-5)
  ## Strong duality for the real case.
  expect_equal(sum(as.matrix(Yr) * Ar), valr, tolerance = 1e-5)
})
