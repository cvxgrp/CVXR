## Tolerance of the symmetric / Hermitian / skew-symmetric tests on constant
## matrices (claims audit C0487, decision 22, CVXR 1.9.2.9155).
##
## CVXPY SOURCE: interface/matrix_utilities.py:276-312 is_hermitian /
## is_skew_symmetric, :315-366 is_sparse_symmetric, :369-413
## is_sparse_skew_symmetric (CVXPY 1.9.3).
##
## Upstream compares with `np.allclose` (elementwise |a - b| <= 1e-8 + 1e-5 |b|,
## False on NaN).  Dense matrices are compared entry by entry; sparse matrices
## must match STRUCTURALLY (same stored positions after summing duplicates and
## dropping explicit zeros) and then allclose on the paired values.  Through
## CVXR 1.9.2.9154 all three tests used `all.equal` at tolerance 1e-10, which is
## a MEAN relative difference, treats NaN == NaN as equal, and has no absolute
## part.  Every expectation below was measured on CVXPY 1.9.3
## (scratch d22_cvxpy.py / d22_skew.py, recorded in
## notes/audit/claims_audit_item1.md, "Decision 22").

## @cvxpy NONE (R twin of numpy.allclose; numpy's own semantics)
test_that(".np_allclose has numpy's rtol/atol, NaN and infinity semantics", {
  np_allclose <- CVXR:::.np_allclose
  ## rtol 1e-5 relative to |b|
  expect_true(np_allclose(1, 1 + 1e-6))
  expect_false(np_allclose(1, 1 + 1e-4))
  ## atol 1e-8 keeps small numbers comparable
  expect_true(np_allclose(0, 5e-9))
  expect_false(np_allclose(0, 2e-8))
  ## asymmetric in (a, b), exactly as numpy: the tolerance scales with |b|
  expect_true(np_allclose(1e6 + 5, 1e6))
  expect_true(np_allclose(1e6, 1e6 + 5))
  expect_false(np_allclose(1 + 1e-4, 1))
  ## complex: modulus of the difference
  expect_true(np_allclose(1 + 1i, 1 + 1i + 1e-6))
  expect_false(np_allclose(1 + 1i, 1 + 1.001i))
  ## NaN anywhere -> FALSE, never NA
  expect_false(np_allclose(c(1, NaN), c(1, NaN)))
  expect_identical(np_allclose(NaN, NaN), FALSE)
  ## infinities must be equal
  expect_true(np_allclose(c(1, Inf), c(1, Inf)))
  expect_false(np_allclose(Inf, -Inf))
  expect_false(np_allclose(Inf, 1e300))
  ## matrices compared entrywise, scalar b recycled
  expect_true(np_allclose(matrix(1e-9, 2, 2), 0))
  expect_true(np_allclose(numeric(0), numeric(0)))
})

## @cvxpy NONE (tolerance of matrix_utilities.py:276-299 is_hermitian, dense)
test_that("dense symmetry test is np.allclose: rtol 1e-5, atol 1e-8, FALSE on NaN", {
  herm <- CVXR:::intf_is_hermitian
  ## off-diagonal pair differs by 1e-7: symmetric upstream, rejected through 9154
  P1 <- matrix(c(2, 1 + 1e-7, 1, 2), 2)
  expect_true(herm(P1)$is_symmetric)
  expect_true(herm(P1)$is_hermitian)
  ## 1e-4 is outside rtol
  P2 <- matrix(c(2, 1 + 1e-4, 1, 2), 2)
  expect_false(herm(P2)$is_symmetric)
  ## tiny magnitudes: atol decides
  P3 <- matrix(c(1e-9, 2e-9, 1e-9, 1e-9), 2)
  expect_true(herm(P3)$is_symmetric)
  ## a 1e-9 entry against an exact 0 is symmetric when DENSE ...
  D2 <- matrix(c(2, 0, 1e-9, 2), 2)
  expect_true(herm(D2)$is_symmetric)
  ## ... and not when SPARSE: the stored positions differ
  S2 <- Matrix::Matrix(D2, sparse = TRUE)
  expect_false(herm(S2)$is_symmetric)
  ## NaN: np.allclose is False; all.equal called it symmetric and is_psd then
  ## crashed in eigen() with "infinite or missing values in 'x'"
  Pn <- matrix(c(1, NaN, NaN, 1), 2)
  expect_false(herm(Pn)$is_symmetric)
  expect_false(herm(Pn)$is_hermitian)
  Cn <- Constant(Pn)
  expect_false(is_symmetric(Cn))
  expect_false(CVXR:::is_psd(Cn))
  expect_false(CVXR:::is_nsd(Cn))
  ## complex: Hermitian within tolerance; "symmetric" is always FALSE
  P5 <- matrix(c(2 + 0i, 1 - 1i + 1e-7, 1 + 1i, 2 + 0i), 2)
  expect_false(herm(P5)$is_symmetric)
  expect_true(herm(P5)$is_hermitian)
  P6 <- matrix(c(2 + 0i, 1 + 1i, 1 + 1i, 2 + 0i), 2)
  expect_false(herm(P6)$is_symmetric)
  expect_false(herm(P6)$is_hermitian)
})

## @cvxpy NONE (matrix_utilities.py:315-366 is_sparse_symmetric)
test_that("sparse symmetry test is structure-exact, then np.allclose on paired values", {
  herm <- CVXR:::intf_is_hermitian
  S1 <- Matrix::Matrix(matrix(c(2, 1 + 1e-7, 1, 2), 2), sparse = TRUE)
  expect_s4_class(S1, "dgCMatrix")
  expect_true(herm(S1)$is_symmetric)
  S4 <- Matrix::Matrix(matrix(c(4, 1 + 1e-7, 0, 0,  1, 4, 1, 0,
                                0, 1, 4, 1,  0, 0, 1, 4), 4), sparse = TRUE)
  expect_true(herm(S4)$is_symmetric)
  S5 <- Matrix::Matrix(matrix(c(4, 1 + 1e-4, 0, 0,  1, 4, 1, 0,
                                0, 1, 4, 1,  0, 0, 1, 4), 4), sparse = TRUE)
  expect_false(herm(S5)$is_symmetric)
  ## Matrix's symmetric and diagonal storage classes store one triangle only;
  ## the test must expand them before comparing positions
  Ssym <- Matrix::Matrix(matrix(c(2, 1, 1, 2), 2), sparse = TRUE)
  expect_s4_class(Ssym, "dsCMatrix")
  expect_true(herm(Ssym)$is_symmetric)
  expect_true(herm(Matrix::Diagonal(3, 1:3))$is_symmetric)
  ## triplet input with duplicates (summed) and an explicit zero (dropped)
  Tdup <- Matrix::sparseMatrix(i = c(2L, 2L, 1L, 3L), j = c(1L, 1L, 2L, 1L),
                               x = c(2, 3, 5, 0), dims = c(3L, 3L), repr = "T")
  expect_true(herm(Tdup)$is_symmetric)
  ## same multiset of values in the two triangles, different positions
  A <- Matrix::sparseMatrix(i = c(2L, 1L), j = c(1L, 3L), x = c(5, 5), dims = c(3L, 3L))
  expect_false(herm(A)$is_symmetric)
})

## @cvxpy NONE (end to end: quad_form accepts what CVXPY accepts)
test_that("quad_form accepts a matrix symmetric to np.allclose tolerance", {
  require_solver("CLARABEL")
  x <- Variable(2)
  P1 <- matrix(c(2, 1 + 1e-7, 1, 2), 2)
  prob <- Problem(Minimize(quad_form(x, P1)), list(x >= 1))
  val <- psolve(prob, solver = "CLARABEL")
  ## CVXPY 1.9.3 + CLARABEL: optimal, 6.00000010989, x = (1, 1)
  expect_equal(status(prob), "optimal")
  expect_equal(val, 6, tolerance = 1e-6)
  S1 <- Matrix::Matrix(P1, sparse = TRUE)
  prob2 <- Problem(Minimize(quad_form(x, S1)), list(x >= 1))
  val2 <- psolve(prob2, solver = "CLARABEL")
  expect_equal(val2, 6, tolerance = 1e-6)
  ## 1e-4 is still rejected, with the upstream message
  P2 <- matrix(c(2, 1 + 1e-4, 1, 2), 2)
  expect_error(quad_form(x, P2), "must be symmetric/Hermitian")
})

## @cvxpy NONE (tolerance of matrix_utilities.py:302-312 is_skew_symmetric and
## :369-413 is_sparse_skew_symmetric)
test_that("skew-symmetry test is np.allclose(A + t(A), 0): absolute 1e-8 only", {
  skew <- CVXR:::intf_is_skew_symmetric
  ## the reference is 0, so rtol contributes nothing: a 1e12 pair off by 1e-3
  ## is NOT skew-symmetric (all.equal's relative test said it was)
  K3 <- matrix(c(0, -1e12 - 1e-3, 1e12, 0), 2)
  expect_false(skew(K3))
  ## a 1e-9 diagonal entry doubles to 2e-9 <= 1e-8: skew-symmetric
  K4 <- matrix(c(1e-9, -1, 1, 0), 2)
  expect_true(skew(K4))
  K5 <- matrix(c(1e-7, -1, 1, 0), 2)
  expect_false(skew(K5))
  ## off-diagonal pair off by 1e-7: 1e-7 > 1e-8, not skew-symmetric
  K1 <- matrix(c(0, -1 - 1e-7, 1, 0), 2)
  expect_false(skew(K1))
  ## NaN -> FALSE (all.equal said TRUE)
  Kn <- matrix(c(0, NaN, NaN, 0), 2)
  expect_false(skew(Kn))
  expect_false(CVXR:::is_skew_symmetric(Constant(Kn)))
  ## sparse: the diagonal entry mirrors itself; a stored entry with no mirror
  ## fails the structural test even when its value is within atol
  expect_true(skew(Matrix::Matrix(K4, sparse = TRUE)))
  expect_false(skew(Matrix::Matrix(K1, sparse = TRUE)))
  S6 <- Matrix::Matrix(matrix(c(0, 0, 1e-9, 0), 2), sparse = TRUE)
  expect_false(skew(S6))
  ## explicit zero on the diagonal is dropped (test_sparse_symmetry_checks_positions K)
  K <- Matrix::sparseMatrix(i = c(2L, 1L, 3L), j = c(1L, 2L, 3L), x = c(5, -5, 0),
                            dims = c(3L, 3L), repr = "T")
  expect_true(skew(K))
  expect_true(CVXR:::is_skew_symmetric(Constant(K)))
})
