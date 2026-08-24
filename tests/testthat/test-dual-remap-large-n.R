## Dual re-keying across MANY constraints, asserted BY INDEX.
##
## Every `reduction_invert` that renames constraints moves the solution's duals
## from new constraint ids back onto the original ones. Those loops are now
## O(n) -- one vectorized `match()` instead of a per-id lookup into a named list
## (see `.remap_by_id_map`, zzz_R_specific/utility.R, and ADR D_PERF.7). A
## mistake there does not crash: it returns plausible-looking duals attached to
## the WRONG constraints.
##
## Before this file, no test in the suite asserted distinct duals across a large
## constraint set -- every multi-constraint dual assertion sat at index 2 or 3 of
## a handful of constraints, which is exactly the regime where an off-by-one or
## an ordering assumption survives. Two specific hazards it exists to catch:
##
##   1. A permuted id map. The oracle below gives constraint i the dual value i,
##      so any permutation fails loudly and names the index.
##   2. Lexicographic key order. The remap loops iterate `ls(env)`, which sorts
##      keys AS STRINGS, so "10" precedes "9". With n >= 200 the ids span three
##      digit widths, so a fix that accidentally depends on numeric order fails
##      here and cannot reach a release.
##
## The values are absolute oracles derived from LP duality, never a comparison
## against a reshaped twin of the same problem -- see the methodological warning
## in test-psd-dual-recovery.R.

skip_if_not_installed("clarabel")

N_CONS <- 250L

# ---- active constraints: dual = objective coefficient ---------------

## @cvxpy NONE
test_that("duals stay attached to their own constraint across 250 constraints", {
  ## minimize sum(i * x_i) s.t. x_i >= b_i, with i > 0, so every lower bound is
  ## ACTIVE at the optimum (x_i = b_i) and its dual is exactly the objective
  ## coefficient i. Distinct per constraint by construction.
  cvec <- as.numeric(seq_len(N_CONS))
  b <- as.numeric((seq_len(N_CONS) %% 7L) - 3L)
  x <- Variable(N_CONS)
  cons <- lapply(seq_len(N_CONS), function(i) x[i] >= b[i])
  p <- Problem(Minimize(sum(cvec * x)), cons)
  psolve(p, solver = "CLARABEL")

  expect_equal(status(p), "optimal")
  expect_equal(as.numeric(value(x)), b, tolerance = 1e-6)

  duals <- vapply(cons, function(cc) as.numeric(dual_value(cc)), numeric(1))
  ## The whole vector at once: a permutation reports every wrong index.
  expect_equal(duals, cvec, tolerance = 1e-6)
})

# ---- active/inactive pattern by index -------------------------------

## @cvxpy NONE
test_that("only the active constraints carry nonzero duals, by index", {
  ## Lower bounds are active (dual = i); upper bounds are slack by 10 and so
  ## are INACTIVE (dual = 0). Both families are present in one problem, so a
  ## remap that shifts duals between families fails even though the multiset of
  ## dual values is unchanged.
  cvec <- as.numeric(seq_len(N_CONS))
  b <- as.numeric((seq_len(N_CONS) %% 7L) - 3L)
  x <- Variable(N_CONS)
  lo <- lapply(seq_len(N_CONS), function(i) x[i] >= b[i])
  hi <- lapply(seq_len(N_CONS), function(i) x[i] <= b[i] + 10)
  p <- Problem(Minimize(sum(cvec * x)), c(lo, hi))
  psolve(p, solver = "CLARABEL")

  expect_equal(status(p), "optimal")
  lo_duals <- vapply(lo, function(cc) as.numeric(dual_value(cc)), numeric(1))
  hi_duals <- vapply(hi, function(cc) as.numeric(dual_value(cc)), numeric(1))
  expect_equal(lo_duals, cvec, tolerance = 1e-6)
  expect_equal(hi_duals, rep(0, N_CONS), tolerance = 1e-6)
})

# ---- the same, with a matrix-shaped constraint in the mix -----------

## @cvxpy NONE
test_that("a matrix constraint keeps its shape while scalar duals stay indexed", {
  ## Exercises the ConeMatrixStuffing reshape branch alongside the bulk remap:
  ## the matrix constraint must come back shaped, and the 250 scalar duals must
  ## still land on their own constraints.
  cvec <- as.numeric(seq_len(N_CONS))
  b <- as.numeric((seq_len(N_CONS) %% 7L) - 3L)
  x <- Variable(N_CONS)
  Y <- Variable(c(2L, 3L))
  cons <- c(lapply(seq_len(N_CONS), function(i) x[i] >= b[i]),
            list(Y >= 1))
  p <- Problem(Minimize(sum(cvec * x) + sum_entries(Y)), cons)
  psolve(p, solver = "CLARABEL")

  expect_equal(status(p), "optimal")
  scal <- vapply(cons[seq_len(N_CONS)],
                 function(cc) as.numeric(dual_value(cc)), numeric(1))
  expect_equal(scal, cvec, tolerance = 1e-6)
  expect_equal(dim(dual_value(cons[[N_CONS + 1L]])), c(2L, 3L))
})

# ---- expr_copy memo: unchanged when a map IS supplied ---------------

## @cvxpy NONE
test_that("expr_copy still memoizes when given an id_objects map", {
  ## `expr_copy` skips the memo entirely when `id_objects` is NULL, because a
  ## freshly allocated environment cannot produce a hit. The path that DOES
  ## supply a map (tree_copy, cvx_attr2constr.R) must be unaffected, so both
  ## halves are pinned: the map is populated, and a second call returns the
  ## identical object rather than rebuilding.
  a <- Variable(3L)
  bexp <- a + a           # AddExpression: has its own expr_copy method
  env <- new.env(hash = TRUE, parent = emptyenv())

  c1 <- expr_copy(bexp, id_objects = env)
  expect_true(exists(as.character(bexp@id), envir = env, inherits = FALSE))
  c2 <- expr_copy(bexp, id_objects = env)
  expect_identical(c1, c2)

  ## With no map, each call builds a fresh copy -- and must still be a copy.
  d1 <- expr_copy(bexp)
  d2 <- expr_copy(bexp)
  expect_false(identical(d1, d2))
  expect_equal(dim(d1), dim(bexp))
})
