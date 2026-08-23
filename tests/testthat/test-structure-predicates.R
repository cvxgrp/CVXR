## @cvxpy cvxpy/atoms/affine/kron.py:96-102
## @cvxpy cvxpy/atoms/affine/transpose.py:71-75
## @cvxpy cvxpy/atoms/affine/trace.py:100-104
##
## Three structure predicates that CVXR under-reported.  All three were
## CONSERVATIVE -- they lost a valid deduction rather than asserting a false one,
## so they rejected valid problems rather than accepting invalid ones.  Audit
## finding #9.
##
##   Kron.is_nsd               PSD (x) NSD is NSD          returned FALSE
##   Transpose.is_skew_symmetric  t() preserves skew-symmetry  returned FALSE
##   Trace.is_real             trace of Hermitian is real  returned FALSE
##
## These predicates are not exported, so the tests reach them through the
## namespace -- which is deliberate: it pins the INTERNAL contract the
## canonicalizers rely on, not a user-facing API.

ns <- asNamespace("CVXR")
.is_skew_symmetric <- get("is_skew_symmetric", ns)
.is_hermitian      <- get("is_hermitian", ns)
.is_real           <- get("is_real", ns)
.is_complex        <- get("is_complex", ns)

# -- Kron ------------------------------------------------------------

test_that("kron(PSD, NSD) and kron(NSD, PSD) are NSD", {
  I2  <- Constant(diag(2))
  nI2 <- Constant(-diag(2))
  P   <- Variable(c(2, 2), PSD = TRUE)
  N   <- Variable(c(2, 2), NSD = TRUE)

  expect_true(is_nsd(kron(I2, nI2)))     # both constant
  expect_true(is_nsd(kron(nI2, I2)))
  expect_true(is_nsd(kron(I2, N)))       # constant x variable
  expect_true(is_nsd(kron(nI2, P)))
})

test_that("kron is_nsd stays FALSE where the deduction is not valid", {
  ## The condition is SUFFICIENT, not necessary -- it must not over-claim.
  I2    <- Constant(diag(2))
  nI2   <- Constant(-diag(2))
  indef <- Constant(matrix(c(1, 0, 0, -1), 2, 2))
  P     <- Variable(c(2, 2), PSD = TRUE)
  N     <- Variable(c(2, 2), NSD = TRUE)

  expect_false(is_nsd(kron(I2, P)))      # PSD (x) PSD is PSD, not NSD
  expect_false(is_nsd(kron(nI2, N)))     # NSD (x) NSD is PSD, not NSD
  expect_false(is_nsd(kron(I2, indef)))  # nothing deducible
  expect_false(is_nsd(kron(indef, P)))
})

test_that("kron is_psd still works — it was never the gap", {
  ## AffAtom's monotonicity propagation already covered both upstream is_psd
  ## cases, which is why no is_psd method was added.  Pinned so a future
  ## "completeness" pass does not add one and narrow this.
  I2  <- Constant(diag(2))
  nI2 <- Constant(-diag(2))
  N   <- Variable(c(2, 2), NSD = TRUE)
  expect_true(is_psd(kron(I2, I2)))
  expect_true(is_psd(kron(nI2, N)))
})

# -- Transpose -------------------------------------------------------

test_that("t() preserves skew-symmetry", {
  ## t(K) == -K for skew-symmetric K, and -K is skew-symmetric too.
  K <- Constant(matrix(c(0, 1, -2,
                         -1, 0, 3,
                         2, -3, 0), 3, 3, byrow = TRUE))
  expect_true(.is_skew_symmetric(K))          # guard on the guard
  expect_true(.is_skew_symmetric(t(K)))
  expect_true(.is_skew_symmetric(t(t(K))))
})

test_that("t() does not invent skew-symmetry", {
  expect_false(.is_skew_symmetric(t(Variable(c(3, 3)))))
  expect_false(.is_skew_symmetric(t(Constant(diag(3)))))   # symmetric, not skew
})

test_that("t() still preserves symmetry and Hermitian-ness", {
  S <- Variable(c(3, 3), symmetric = TRUE)
  H <- Variable(c(3, 3), hermitian = TRUE)
  expect_true(is_symmetric(t(S)))
  expect_true(.is_hermitian(t(H)))
})

# -- Trace -----------------------------------------------------------

test_that("the trace of a Hermitian matrix is real", {
  H <- Variable(c(2, 2), hermitian = TRUE)
  expect_true(.is_hermitian(H))               # guard on the guard
  expect_false(.is_real(H))                   # the matrix itself is complex
  expect_true(.is_real(matrix_trace(H)))      # its trace is not
  expect_false(.is_complex(matrix_trace(H)))
})

test_that("trace is_real still follows a real argument", {
  S <- Variable(c(2, 2), symmetric = TRUE)
  expect_true(.is_real(matrix_trace(S)))
  expect_false(.is_complex(matrix_trace(S)))
})

test_that("trace does not claim a genuinely complex trace is real", {
  ## A plain complex variable is neither real nor Hermitian, so nothing is
  ## deducible and the trace must stay complex.
  Z <- Variable(c(2, 2), complex = TRUE)
  expect_false(.is_real(Z))
  expect_false(.is_hermitian(Z))
  expect_false(.is_real(matrix_trace(Z)))
  expect_true(.is_complex(matrix_trace(Z)))
})

test_that("a Hermitian trace solves, and returns a real value", {
  ## The point of the predicate: Complex2Real no longer splits this into real
  ## and imaginary parts with a redundant Im(.) == 0 row.
  H <- Variable(c(2, 2), hermitian = TRUE)
  prob <- Problem(Minimize(matrix_trace(H)),
                  list(H %>>% diag(2), matrix_trace(H) <= 10))
  val <- psolve(prob, solver = "SCS")
  expect_true(is.finite(drop(Re(val))))
  expect_equal(drop(Im(as.complex(val))), 0, tolerance = 1e-6)
})
