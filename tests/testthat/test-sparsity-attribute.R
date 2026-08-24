## @cvxpy tests/test_expressions.py, tests/test_problem.py
##
## The `sparsity` leaf attribute must actually constrain the solve.
##
## Why this file exists: it did not. `sparsity` was accepted by the
## constructor, stored in the attribute list, listed in CONVEX_ATTRIBUTES and in
## the reduction's own `reduction_attributes` -- and then consumed by nothing.
## `Variable(c(2,2), sparsity = <diagonal>)` solved as a fully dense variable.
## The header of cvx_attr2constr.R said "Deferred: sparsity" while the
## constructor advertised it, which is the worst of the three possible states:
## a user gets no error, no warning, and a different answer.
##
## Two structural details made it possible, and both are fixed here:
##   * `isTRUE(attrs$sparsity)` was the test at every site. `sparsity` is an
##     INDEX SET, so `isTRUE` is FALSE for every real pattern and each site was
##     told "not set". `.attr_set()` is Python's truthiness rule and is now the
##     test.
##   * the attribute-clearing loop in `reduction_apply` hardcoded four keys and
##     omitted "sparsity", while `reduction_attributes` twenty lines above
##     listed five. Two lists for one concept, disagreeing.
##
## CVXPY SOURCE: leaf.py:148-152 (sparse_idx), :714-728 (_has_dim_reducing_attr,
## _reduced_size), cvx_attr2constr.py:32-43 (CONVEX_ATTRIBUTES), :155-161
## (build_dim_reduced_expression's sparse arm), :93-97 (recover_value_for_leaf).
##
## R DEVIATION, deliberate and consistent with `boolean`/`integer`: CVXPY takes
## a numpy multi-index (a tuple of per-dimension arrays) and RAISES TypeError on
## `sparsity=True`. CVXR accepts the same four spellings it accepts for
## boolean/integer -- TRUE, flat column-major indices, a two-column (row, col)
## matrix, and a logical mask -- all canonicalized by `.mip_idx()` to a 1-based
## flat index. See the block comment above `.mip_idx` in expressions/leaf.R.

.sp_diag2 <- function() cbind(c(1, 2), c(1, 2))   # the 2x2 diagonal

test_that("every spelling canonicalizes to the same .sparse_idx", {
  expect_equal(as.integer(Variable(c(2, 2), sparsity = c(1, 4))@.sparse_idx), c(1L, 4L))
  expect_equal(as.integer(Variable(c(2, 2), sparsity = .sp_diag2())@.sparse_idx), c(1L, 4L))
  expect_equal(as.integer(Variable(c(2, 2),
    sparsity = matrix(c(TRUE, FALSE, FALSE, TRUE), 2, 2))@.sparse_idx), c(1L, 4L))
  ## TRUE means "every entry", as it does for boolean/integer.
  expect_equal(as.integer(Variable(c(2, 2), sparsity = TRUE)@.sparse_idx), 1:4)
  ## Parameters go through the same canonicalizer.
  expect_equal(as.integer(Parameter(c(2, 2), sparsity = c(1, 4))@.sparse_idx), c(1L, 4L))
})

test_that("a bad sparsity index is rejected with a real message", {
  ## Pre-fix EVERY index spelling died on `if (sparsity)` with R's bare
  ## "the condition has length > 1", before reaching any validation at all.
  expect_error(Variable(c(2, 2), sparsity = c(1, 99)), "outside")
  expect_error(Variable(c(2, 2), sparsity = c(1, 1)), "unique")
  expect_error(Variable(c(2, 2), sparsity = c(1.5, 2)), "whole numbers")
  expect_error(Variable(c(2, 2), sparsity = matrix(c(TRUE, FALSE), 1, 2)), "mask")
})

test_that("sparsity constrains the solve on every solver", {
  ## min sum(X) s.t. X >= -1, X sparse on the diagonal. Only the two diagonal
  ## entries are free, so the optimum is -2. Pre-fix CVXR returned -4, the
  ## fully dense answer. CVXPY 1.9.2 returns -1.9999999999430869.
  for (s in intersect(c("CLARABEL", "OSQP", "ECOS", "SCS"), installed_solvers())) {
    X <- Variable(c(2, 2), sparsity = .sp_diag2())
    prob <- Problem(Minimize(sum_entries(X)), list(X >= -1))
    expect_equal(as.numeric(psolve(prob, solver = s)), -2, tolerance = 1e-4,
                 label = paste("objective under", s))
    expect_equal(status(prob), "optimal", label = paste("status under", s))
    ## Off-pattern entries are structural zeros, not merely small.
    v <- as.numeric(value(X))
    expect_equal(v[c(2L, 3L)], c(0, 0), tolerance = 1e-9,
                 label = paste("off-pattern entries under", s))
  }
})

test_that("an off-diagonal pattern is honored, not just the diagonal", {
  ## Pattern (1,2), (2,3), (3,1) -- three entries of a 3x3, none on the
  ## diagonal, so a diag-shaped implementation would get this wrong.
  X <- Variable(c(3, 3), sparsity = cbind(c(1, 2, 3), c(2, 3, 1)))
  prob <- Problem(Minimize(sum_entries(X)), list(X >= -1))
  expect_equal(as.numeric(psolve(prob, solver = "CLARABEL")), -3, tolerance = 1e-6)
  m <- matrix(as.numeric(value(X)), 3, 3)
  expect_equal(c(m[1, 2], m[2, 3], m[3, 1]), c(-1, -1, -1), tolerance = 1e-6)
  ## every other entry is zero
  m[1, 2] <- 0; m[2, 3] <- 0; m[3, 1] <- 0
  expect_equal(as.numeric(m), rep(0, 9), tolerance = 1e-9)
})

test_that("a dense variable is unaffected by the sparsity machinery", {
  Y <- Variable(c(2, 2))
  expect_equal(as.numeric(psolve(Problem(Minimize(sum_entries(Y)), list(Y >= -1)),
                                 solver = "CLARABEL")), -4, tolerance = 1e-6)
})

test_that("the leaf is dimension-reducing and reports the right reduced size", {
  w <- Variable(c(3, 3), sparsity = cbind(1:3, 1:3))
  expect_true(CVXR:::.cvxattr_has_dim_reducing_attr(w))
  expect_equal(CVXR:::.cvxattr_reduced_size(w), 3L)
  expect_equal(CVXR:::convex_attributes(list(w)), "sparsity")
  ## A dense leaf is neither.
  d <- Variable(c(3, 3))
  expect_false(CVXR:::.cvxattr_has_dim_reducing_attr(d))
  expect_equal(CVXR:::.cvxattr_reduced_size(d), 9L)
})

test_that("project() and value<- enforce the pattern", {
  v <- Variable(c(2, 2), sparsity = c(1, 4))
  expect_equal(as.numeric(CVXR:::project(v, matrix(0.7, 2, 2))), c(0.7, 0, 0, 0.7))
  ## An off-pattern nonzero must be rejected, and the message must name the
  ## reason rather than falling through to "value must be real."
  expect_error({ value(v) <- matrix(c(1, 2, 3, 4), 2, 2) }, "sparsity pattern")
  expect_no_error({ value(v) <- matrix(c(1, 0, 0, 4), 2, 2) })
  expect_equal(as.numeric(value(v)), c(1, 0, 0, 4))
})

test_that("a value set on the leaf survives the reduction", {
  ## The lowered nnz-vector must receive the stored entries, and the recovered
  ## solution must scatter back into full shape.
  v <- Variable(c(2, 2), sparsity = c(1, 4))
  value(v) <- matrix(c(3, 0, 0, 5), 2, 2)
  expect_equal(as.numeric(CVXR:::.cvxattr_lower_value(v)), c(3, 5))
  expect_equal(as.numeric(CVXR:::.cvxattr_build_full_value(v, c(3, 5))),
               c(3, 0, 0, 5))
})

test_that("sparsity still cannot be combined with the attributes it contradicts", {
  ## pos/neg force strict signs, sparsity forces zeros. CVXPY leaf.py:161-167.
  expect_error(Variable(c(2, 2), sparsity = c(1), pos = TRUE), "Cannot combine")
  expect_error(Variable(c(2, 2), sparsity = c(1), neg = TRUE), "Cannot combine")
  ## and it is one of the mutually exclusive dimension-reducing attributes.
  expect_error(Variable(c(2, 2), sparsity = c(1), diag = TRUE), "more than one")
  expect_error(Variable(c(2, 2), sparsity = c(1), symmetric = TRUE), "more than one")
  expect_error(Variable(c(2, 2), sparsity = c(1), PSD = TRUE), "more than one")
})

test_that("a sparse variable works inside a larger problem", {
  ## Least squares with a sparsity-constrained coefficient matrix: the zeros
  ## must bind, so the fit is strictly worse than the dense one.
  set.seed(11)
  Amat <- matrix(rnorm(12), 4, 3)
  bvec <- rnorm(4)
  Xs <- Variable(c(3, 1), sparsity = c(1, 3))     # middle coefficient forced to 0
  fit_sparse <- psolve(Problem(Minimize(sum_squares(Amat %*% Xs - bvec))),
                       solver = "CLARABEL")
  expect_equal(as.numeric(value(Xs))[2], 0, tolerance = 1e-9)

  Xd <- Variable(c(3, 1))
  fit_dense <- psolve(Problem(Minimize(sum_squares(Amat %*% Xd - bvec))),
                      solver = "CLARABEL")
  expect_gt(as.numeric(fit_sparse), as.numeric(fit_dense))
})
