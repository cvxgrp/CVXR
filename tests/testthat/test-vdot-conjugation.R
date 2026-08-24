## @cvxpy cvxpy/atoms/affine/binary_operators.py:590-660
##
## vdot() is CONJUGATE-LINEAR in its first argument.
##
##     x = deep_flatten(x); y = deep_flatten(y)
##     prod = multiply(conj(x), y)          # <- conj
##     return cvxpy_sum(prod)
##
## CVXR computed `sum(vec(x) * vec(y))` with no conjugation, so for complex
## input it returned the bilinear form instead of the Hermitian inner product --
## silently, and only for complex arguments. numpy.vdot conjugates too, so the
## R answer disagreed with both CVXPY and numpy.
##
## The wrap is UNCONDITIONAL, matching upstream. conj() is non-monotone
## (conj.py:41-49), so it drops a real CONVEX argument to unknown curvature in
## CVXPY as well -- measured: cp.vdot(cp.square(z), c) is not DCP upstream.
## Skipping the wrap for real input would leave CVXR more permissive than
## CVXPY, which is the shape of the conv() defect fixed just before this.
##
## Ledger + paired CVXPY reproducer:
##   notes/audit/completeness_ledger.md finding 2
##   notes/audit/repro/02_vdot_conjugation.{py,R}

library(testthat)
library(CVXR)

test_that("vdot conjugates its first argument (complex)", {
  z <- c(1 + 2i, 3 - 1i)
  w <- c(2 + 0i, 1 + 1i)

  ## numpy.vdot(z, w) == sum(Conj(z) * w) == 4+0i.  The un-conjugated product
  ## sum(z * w) is 6+6i -- a DIFFERENT number, which is what makes this test
  ## discriminating rather than a tautology.
  expect_equal(drop(value(vdot(z, w))), sum(Conj(z) * w))
  expect_equal(drop(value(scalar_product(z, w))), sum(Conj(z) * w))
  expect_false(isTRUE(all.equal(drop(value(vdot(z, w))), sum(z * w))))
})

test_that("vdot is conjugate-linear, not bilinear, in argument order", {
  ## <x, y> = conj(<y, x>).  A bilinear (un-conjugated) form would instead be
  ## symmetric, so this pins the conjugation from a second direction.
  ##
  ## NOTE the operands: <z, w> must have a NONZERO IMAGINARY PART or the
  ## asymmetry assertion is vacuous. The pair used in the test above gives
  ## 4+0i, for which <z,w> and <w,z> agree and a bilinear form would pass this
  ## too. These operands give 6+1i.
  z <- c(1 + 2i, 3 - 1i)
  w <- c(2 + 1i, 1 + 1i)
  expect_true(Im(drop(value(vdot(z, w)))) != 0)          # guard the guard
  expect_equal(drop(value(vdot(z, w))), Conj(drop(value(vdot(w, z)))))
  expect_false(isTRUE(all.equal(drop(value(vdot(z, w))),
                                drop(value(vdot(w, z))))))
})

test_that("vdot with a complex Variable conjugates the variable", {
  x <- Variable(2, complex = TRUE)
  cst <- c(2 + 0i, 1 + 1i)
  value(x) <- c(1 + 2i, 3 - 1i)
  expect_equal(drop(value(vdot(x, cst))), sum(Conj(c(1 + 2i, 3 - 1i)) * cst))
})

test_that("vdot is unchanged for real arguments", {
  ## Conjugation is the identity on reals, so every existing real use must
  ## return exactly what it did before.
  x <- Variable(3)
  value(x) <- c(1, 2, 3)
  expect_equal(drop(value(vdot(x, c(1, 2, 3)))), 14)
  expect_equal(drop(value(vdot(c(1, 2, 3), c(1, 2, 3)))), 14)
  expect_equal(drop(value(scalar_product(c(1, 2, 3), c(4, 5, 6)))), 32)
  expect_true(is_affine(vdot(x, c(1, 2, 3))))
})

test_that("vdot of a real convex argument matches CVXPY's unknown curvature", {
  ## Measured upstream: cp.vdot(cp.square(z), c).is_convex() is False and the
  ## minimization is not DCP, because conj() is non-monotone. CVXR must agree
  ## rather than being more permissive.
  z <- Variable(3)
  e <- vdot(square(z), c(1, 2, 3))
  expect_false(is_convex(e))
  expect_false(is_concave(e))
  expect_false(is_dcp(Problem(Minimize(e))))
})

test_that("vdot solves a real least-squares problem correctly", {
  ## End-to-end guard: the extra Conj_ node must not disturb canonicalization
  ## on the purely real path (its graph_implementation is the identity).
  x <- Variable(3)
  a <- c(1, 2, 3)
  prob <- Problem(Minimize(sum_squares(x - a)), list(vdot(x, c(1, 1, 1)) == 3))
  val <- psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
  expect_equal(sum(value(x)), 3, tolerance = 1e-6)
  expect_true(is.finite(val))
})

test_that("vdot accepts nested lists, flattened independently", {
  ## CVXPY SOURCE: binary_operators.py:616-617 (deep_flatten on both args) and
  ## reshape.py:164-184. Measured upstream: cp.vdot([[a],[b]], [1., 2.]) with
  ## a = 3, b = 4 gives 11.0.
  a <- Variable(); b <- Variable()
  value(a) <- 3; value(b) <- 4

  expect_equal(drop(value(vdot(list(a, b), c(1, 2)))), 11)
  expect_equal(drop(value(vdot(list(list(a), list(b)), c(1, 2)))), 11)
  ## Nesting need not match between the two arguments.
  expect_equal(drop(value(vdot(list(list(a), list(b)), list(1, 2)))), 11)
  expect_equal(drop(value(vdot(list(a, b), list(list(1), 2)))), 11)
})

test_that("deep_flatten stacks into a COLUMN, not a row", {
  ## The R-specific trap: CVXPY flattens to 1-D and concatenates with hstack;
  ## CVXR shapes are always 2-D, so the pieces are (n, 1) columns and the
  ## concatenation must be vstack. Getting it backwards yields an (n, k)
  ## matrix that still "works" for some downstream ops.
  a <- Variable(); b <- Variable()
  value(a) <- 3; value(b) <- 4
  f <- deep_flatten(list(a, b))
  expect_equal(f@shape, c(2L, 1L))
  expect_equal(as.vector(value(f)), c(3, 4))

  ## A matrix flattens column-major, matching order = "F".
  M <- Constant(matrix(c(1, 2, 3, 4), 2, 2))
  expect_equal(deep_flatten(M)@shape, c(4L, 1L))
  expect_equal(as.vector(value(deep_flatten(M))), c(1, 2, 3, 4))

  ## Mixed nesting and depth.
  g <- deep_flatten(list(M, list(a, c(5, 6))))
  expect_equal(g@shape, c(7L, 1L))
  expect_equal(as.vector(value(g)), c(1, 2, 3, 4, 3, 5, 6))

  expect_error(deep_flatten(list()), "empty list")
  expect_error(deep_flatten("nope"), "cannot flatten")
})

test_that("cvxr_outer rejects nested lists and matrices, as CVXPY does", {
  ## Upstream's DOCSTRING claims outer() takes nested lists and flattens
  ## non-vectors; its BODY does neither. Measured against cvxpy 1.9.2, both
  ## raise "x must be a 1-d array." CVXR must reject them too -- this test
  ## exists so nobody "fixes" cvxr_outer toward the docstring.
  expect_error(cvxr_outer(Constant(matrix(1:4, 2, 2)), c(1, 2)), "must be a vector")
  expect_error(cvxr_outer(c(1, 2), Constant(matrix(1:4, 2, 2))), "must be a vector")
})

test_that("cvxr_outer is unchanged and matches CVXPY's outer", {
  ## outer() was ported correctly; it moved files with vdot, so pin it.
  z <- c(1 + 2i, 3 - 1i)
  w <- c(2 + 0i, 1 + 1i)
  o <- cvxr_outer(z, w)
  expect_equal(o@shape, c(2L, 2L))
  ## CVXPY's outer does NOT conjugate: outer(z, w)[0,0] == z[0] * w[0].
  expect_equal(value(o)[1, 1], z[1] * w[1])
  expect_equal(value(o), outer(z, w))

  x <- Variable(3)
  expect_error(cvxr_outer(Variable(c(2, 2)), x), "must be a vector")
})
