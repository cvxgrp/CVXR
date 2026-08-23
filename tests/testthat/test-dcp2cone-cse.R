## Common-subexpression elimination in Dcp2Cone (CVXPY 1.9.2, PR #3355).
##
## When the same subtree appears twice in a problem, the canonicalizer must emit
## ONE set of auxiliary variables and epigraph constraints, not two. These tests
## are the R port of CVXPY's test_dcp2cone_cse.py.

## Auxiliary variables introduced by canonicalization: everything except the
## user's own variables (compared by id, since canonicalization copies nodes).
.aux_vars <- function(new_prob, ...) {
  user_ids <- vapply(list(...), function(v) v@id, integer(1L))
  Filter(function(v) !(v@id %in% user_ids), variables(new_prob))
}

.apply_dcp2cone <- function(prob, quad_obj = FALSE) {
  reduction_apply(Dcp2Cone(quad_obj = quad_obj), prob)[[1L]]
}

## All atoms of a class reachable from an expression (upstream `_find_atoms`).
.find_atoms <- function(expr, cls) {
  found <- list()
  stack <- list(expr)
  while (length(stack) > 0L) {
    e <- stack[[length(stack)]]; stack[[length(stack)]] <- NULL
    if (.s7_is(e, cls)) found[[length(found) + 1L]] <- e
    if (.s7_is(e, Expression) && length(e@args) > 0L) stack <- c(stack, e@args)
  }
  found
}

## @cvxpy test_dcp2cone_cse.py::TestDcp2ConeCSE::test_scalar_norm1_dedup
test_that("a scalar norm1 in objective and constraint yields one epigraph variable", {
  x <- Variable()
  prob <- Problem(Minimize(norm1(x)), list(norm1(x) <= 1))
  new_prob <- .apply_dcp2cone(prob)

  expect_length(variables(new_prob), 2L)
  aux <- .aux_vars(new_prob, x)
  expect_length(aux, 1L)
  expect_equal(prod(aux[[1L]]@shape), 1L)
})

## @cvxpy test_dcp2cone_cse.py::TestDcp2ConeCSE::test_vector_norm1_dedup
test_that("a vector norm1 used twice yields one epigraph variable of its shape", {
  x <- Variable(3)
  prob <- Problem(Minimize(norm1(x)), list(norm1(x) <= 5))
  aux <- .aux_vars(.apply_dcp2cone(prob), x)
  expect_length(aux, 1L)
  expect_equal(prod(aux[[1L]]@shape), 3L)
})

## @cvxpy test_dcp2cone_cse.py::TestDcp2ConeCSE::test_distinct_subtrees_not_merged
test_that("structurally different subtrees are not merged", {
  x <- Variable(3)
  prob <- Problem(Minimize(norm1(x) + norm1(-x)), list(x >= -1))
  ## norm1(x) and norm1(-x) differ below the root, so each needs its own
  ## epigraph variable. A key that ignored child structure would merge them.
  expect_length(.aux_vars(.apply_dcp2cone(prob), x), 2L)
})

## @cvxpy test_dcp2cone_cse.py::TestDcp2ConeCSE::test_solve_matches_unduplicated
test_that("a duplicated subtree solves identically to an explicitly shared one", {
  skip_if_not_installed("clarabel")
  x <- Variable(3)
  prob_dup <- Problem(Minimize(norm1(x) + sum_entries(x)),
                      list(norm1(x) <= 4, sum_entries(x) >= 1))
  val_dup <- psolve(prob_dup, solver = "CLARABEL")

  y <- Variable(3)
  shared <- norm1(y)
  prob_shared <- Problem(Minimize(shared + sum_entries(y)),
                         list(shared <= 4, sum_entries(y) >= 1))
  val_shared <- psolve(prob_shared, solver = "CLARABEL")

  expect_equal(status(prob_dup), "optimal")
  expect_equal(val_dup, val_shared, tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), as.numeric(value(y)), tolerance = 1e-5)
})

## The coefficient extractor must keep the two occurrences as DISTINCT rows even
## though CSE makes them share one canonical SymbolicQuadForm -- if the shared
## node were counted once, the objective would be halved.
## @cvxpy test_dcp2cone_cse.py::TestDcp2ConeCSE::test_shared_quad_form_solves_correctly
test_that("a shared QuadForm reused twice solves correctly", {
  skip_if_not_installed("clarabel")
  set.seed(0)
  A <- matrix(rnorm(25), 5, 5)
  z <- rnorm(5)
  P <- t(A) %*% A
  q <- -2 * (P %*% z)
  w <- Variable(5)
  qf <- quad_form(w, P)
  prob <- Problem(Minimize(0.5 * qf + 0.5 * qf + t(q) %*% w))
  psolve(prob, solver = "CLARABEL")
  expect_equal(as.numeric(value(w)), as.numeric(z), tolerance = 1e-4)
})

## @cvxpy test_dcp2cone_cse.py::TestDcp2ConeCSE::test_quad_obj_shared_subtree_dedup
test_that("with quad_obj, a repeated quad-eligible subtree shares one SymbolicQuadForm", {
  skip_if_not_installed("clarabel")
  x <- Variable(3)
  qol <- quad_over_lin(x, 1)
  prob <- Problem(Minimize(qol + qol + sum_entries(x)))
  new_prob <- .apply_dcp2cone(prob, quad_obj = TRUE)

  sqfs <- .find_atoms(new_prob@objective@args[[1L]], SymbolicQuadForm)
  expect_length(sqfs, 2L)
  ## Both occurrences must be the SAME node, not two equal ones.
  expect_true(all(vapply(sqfs, function(s) s@id == sqfs[[1L]]@id, logical(1L))))

  y <- Variable(3)
  ref <- Problem(Minimize(2 * quad_over_lin(y, 1) + sum_entries(y)))
  expect_equal(psolve(prob, solver = "CLARABEL"),
               psolve(ref, solver = "CLARABEL"), tolerance = 1e-5)
  expect_equal(as.numeric(value(x)), as.numeric(value(y)), tolerance = 1e-5)
})

## The same subtree in the objective (quad branch) and in a constraint (cone
## branch) must NOT share a result. This is what the affine_above component of
## the cache key is for: the structural keys match, the contexts do not.
## @cvxpy test_dcp2cone_cse.py::TestDcp2ConeCSE::test_quad_obj_cross_context_not_merged
test_that("with quad_obj, objective and constraint uses of one subtree stay separate", {
  skip_if_not_installed("clarabel")
  x <- Variable(3)
  qol <- quad_over_lin(x, 1)
  prob <- Problem(Minimize(qol + sum_entries(x)), list(qol <= 5))
  new_prob <- .apply_dcp2cone(prob, quad_obj = TRUE)

  expect_length(.find_atoms(new_prob@objective@args[[1L]], SymbolicQuadForm), 1L)
  expect_equal(sum(vapply(new_prob@constraints,
                          function(c) .s7_is(c, SOC), logical(1L))), 1L)

  psolve(prob, solver = "CLARABEL")
  expect_equal(status(prob), "optimal")
})

## R DIFFERENCE, deliberate. Upstream merges two SEPARATE `Constant(arr)` calls
## on one large array, because numpy's ndarray interface stores the array by
## reference so both wrappers share `id(value)`. R has no such exposed object
## identity, so CVXR keys a large Constant by the wrapper's own id: binding a
## Constant once and reusing it merges (upstream's branch 4, and the common R
## idiom), while two separate constructions do not. Small constants
## (<= 64 entries) are keyed BY VALUE in both, which is the case that matters
## for defaults such as huber()'s M = 0.5.
## @cvxpy test_dcp2cone_cse.py::TestDcp2ConeCSE::test_shared_ndarray_constant_dedup
test_that("a large Constant deduplicates when the same object is reused", {
  set.seed(0)
  arr <- rnorm(100)
  x <- Variable(100)

  shared <- Constant(arr)
  prob_shared <- Problem(Minimize(norm1(shared * x)),
                         list(norm1(shared * x) <= 5, x >= -1))
  aux <- .aux_vars(.apply_dcp2cone(prob_shared), x)
  expect_length(aux, 1L)
  expect_equal(prod(aux[[1L]]@shape), 100L)

  ## Two separate constructions: NOT merged in R (see the note above).
  c1 <- Constant(arr); c2 <- Constant(arr)
  prob_two <- Problem(Minimize(norm1(c1 * x)),
                      list(norm1(c2 * x) <= 5, x >= -1))
  expect_length(.aux_vars(.apply_dcp2cone(prob_two), x), 2L)
})

## @cvxpy test_dcp2cone_cse.py::TestDcp2ConeCSE::test_sparse_constant_key_uses_sparse_contents
test_that("sparse Constants key on contents, not storage order", {
  m_coo <- Matrix::sparseMatrix(i = c(1, 2), j = c(2, 1), x = c(1, 2), dims = c(2, 2))
  ## Same matrix, entries supplied in the other order.
  m_same <- Matrix::sparseMatrix(i = c(2, 1), j = c(1, 2), x = c(2, 1), dims = c(2, 2))
  ## Same matrix, one entry given as two duplicates that sum to it.
  m_dup <- Matrix::sparseMatrix(i = c(1, 1, 2), j = c(2, 2, 1), x = c(0.25, 0.75, 2),
                                dims = c(2, 2))
  m_diff <- Matrix::sparseMatrix(i = c(1, 2), j = c(2, 1), x = c(1, 3), dims = c(2, 2))

  kc <- StructuralKeyCache()
  k      <- expr_key(Constant(m_coo), kc)
  k_same <- expr_key(Constant(m_same), kc)
  k_dup  <- expr_key(Constant(m_dup), kc)
  k_diff <- expr_key(Constant(m_diff), kc)

  expect_equal(k, k_same)
  expect_equal(k, k_dup)
  expect_false(k == k_diff)
})

## @cvxpy test_dcp2cone_cse.py::TestDcp2ConeCSE::test_parameter_subtree_dedup
test_that("subtrees sharing a Parameter deduplicate", {
  p <- Parameter(2, nonneg = TRUE)
  value(p) <- c(1, 2)
  x <- Variable(2)
  prob <- Problem(Minimize(norm1(p * x)),
                  list(norm1(p * x) <= 10, x >= 0))
  aux <- .aux_vars(.apply_dcp2cone(prob), x)
  expect_length(aux, 1L)
  expect_equal(prod(aux[[1L]]@shape), 2L)
})

## Upstream spies on `expr_key` to prove one key cache is threaded through the
## whole apply. R has no equivalent monkey-patch that is safe inside a package
## namespace, so this asserts the OBSERVABLE consequence instead: one cache
## memoizes every node of a deep chain, so interning 21 nested subtrees costs 21
## signatures, and re-keying the same tree adds none.
## @cvxpy test_dcp2cone_cse.py::TestDcp2ConeCSE::test_structural_key_cache_threaded_through_apply
test_that("one structural key cache memoizes every node of a deep chain", {
  x <- Variable()
  expr <- x
  for (i in 1:20) expr <- expr + 1

  kc <- StructuralKeyCache()
  k1 <- expr_key(expr, kc)
  expect_type(k1, "integer")
  ## The memo used to be an environment countable with ls(); it is a C++ map
  ## now, so the same observable comes from the store itself.
  n_after_first <- .CseKeys__memo_size(kc)

  k2 <- expr_key(expr, kc)          # same tree again: pure cache hits
  expect_equal(k1, k2)
  expect_equal(.CseKeys__memo_size(kc), n_after_first)
  ## Re-keying interns no new signatures either.
  n_sig <- .CseKeys__n_signatures(kc)
  expr_key(expr, kc)
  expect_equal(.CseKeys__n_signatures(kc), n_sig)

  ## 21 nodes: the Variable plus 20 additions (each `+ 1` also interns its
  ## Constant, so node count is what is asserted, not signature count).
  expect_gte(n_after_first, 21L)

  ## The same cache is threaded through a real apply (upstream asserts this by
  ## spying; here the observable is that the reduction completes and the chain
  ## canonicalizes to a single affine expression). Deliberately not solved:
  ## `min x` with only an upper bound is unbounded, and upstream does not solve
  ## it either -- it only applies the reduction.
  prob <- Problem(Minimize(expr), list(expr <= 25))
  new_prob <- .apply_dcp2cone(prob)
  expect_length(.aux_vars(new_prob, x), 0L)
})
