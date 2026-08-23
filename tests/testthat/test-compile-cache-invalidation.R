## @cvxpy cvxpy/problems/problem.py:831-841
##
## A compile that ABORTS must leave the cache empty, never stale.
##
## CVXPY:
##     if key != self._cache.key:
##         self._cache.invalidate()                      # clear FIRST
##         solving_chain = self._construct_chain(...)     # may raise
##         self._cache.key = key                          # record on success only
##         self._cache.solving_chain = solving_chain
##
## CVXR's .compile() recorded the key BEFORE building and cleared nothing, so an
## abort from construct_solving_chain() left the NEW key paired with the
## PREVIOUS chain. The next call with the same arguments matched the cache and
## got that stale chain back -- so `psolve(prob, gp = TRUE)` on a non-DGP
## problem correctly aborted the first time and then silently solved it as an
## ordinary DCP problem the second time. With no prior successful compile the
## same path yielded a NULL chain and "no applicable method for `@` applied to
## an object of class NULL".
##
## notes/session_handoff_2026-08-17_completeness_audit_and_fixes.md section 5.1

library(testthat)
library(CVXR)

test_that("a repeated failing compile keeps reporting the real error", {
  ## No prior successful compile: the pre-fix failure mode here was a NULL
  ## chain and an unrelated `@`-on-NULL message on the SECOND call.
  x <- Variable(c(2, 2))
  y <- Variable(c(2, 2))
  bad <- Problem(Minimize(sum_entries(x %*% y)), list(x >= 1, y >= 1))
  expect_false(is_dcp(bad))
  for (i in 1:3) {
    expect_error(psolve(bad), "not DCP compliant")
  }
})

test_that("an aborted compile does not leave a stale chain behind", {
  ## THE serious case: a SUCCESSFUL compile, then a failing one with a
  ## different cache key. Pre-fix, the second gp = TRUE call returned the
  ## non-gp chain's answer.
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  expect_true(is_dcp(prob))
  expect_false(is_dgp(prob))

  expect_equal(psolve(prob, solver = "CLARABEL"), 2, tolerance = 1e-6)

  ## Every subsequent gp = TRUE must abort. Repeat: once was never the problem.
  for (i in 1:3) {
    expect_error(psolve(prob, gp = TRUE), "not DGP compliant")
  }

  ## ... and the original request must still work afterwards, i.e. the
  ## invalidation did not throw away a usable cache permanently.
  expect_equal(psolve(prob, solver = "CLARABEL"), 2, tolerance = 1e-6)
})

test_that("the cache is empty, not stale, after an abort", {
  ## Inspect directly rather than inferring from behavior.
  x <- Variable(2)
  prob <- Problem(Minimize(sum_entries(x)), list(x >= 1))
  invisible(psolve(prob, solver = "CLARABEL"))
  expect_false(is.null(prob@.cache$compile_chain))   # populated by the good solve

  try(psolve(prob, gp = TRUE), silent = TRUE)
  expect_null(prob@.cache$compile_key)
  expect_null(prob@.cache$compile_chain)
  expect_null(prob@.cache$param_prog)
  expect_null(prob@.cache$compile_inverse_data)
})

test_that("a successful re-compile under a new key still swaps the chain", {
  ## Guard the other direction: invalidating must not break ordinary
  ## re-compilation when the key changes and the build succeeds.
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x)), list(x >= 1))
  v1 <- psolve(prob, solver = "CLARABEL")
  v2 <- psolve(prob, solver = "SCS")
  expect_equal(v1, 2, tolerance = 1e-5)
  expect_equal(v2, 2, tolerance = 1e-5)
  expect_false(is.null(prob@.cache$compile_chain))
  ## Back again, to exercise the cache-hit path.
  expect_equal(psolve(prob, solver = "CLARABEL"), 2, tolerance = 1e-5)
})

test_that("repeated identical solves still hit the cache", {
  ## The invalidation must not defeat caching for the normal path: same key
  ## twice must NOT rebuild.
  x <- Variable(2)
  prob <- Problem(Minimize(sum_squares(x)), list(x >= 1))
  invisible(psolve(prob, solver = "CLARABEL"))
  key1 <- prob@.cache$compile_key
  chain1 <- prob@.cache$compile_chain
  invisible(psolve(prob, solver = "CLARABEL"))
  expect_identical(prob@.cache$compile_key, key1)
  expect_identical(prob@.cache$compile_chain, chain1)   # same object, not rebuilt
})
