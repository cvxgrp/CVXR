## CVXPY SOURCE: utilities/performance_utils.py (compute_once / lazyprop)
##
## CVXPY has no test module for performance_utils, so these are R-side unit
## tests of the one helper every member cache goes through. The property
## under test is _cache_key's: a cached member answers per (object, member,
## dpp scope, quad_form_dpp scope), never across scopes.

## @cvxpy NONE
test_that("compute_once evaluates once per object and returns the cached value", {
  x <- Constant(matrix(c(2, 1, 1, 3), 2, 2))
  n <- 0L
  fn <- function(obj) { n <<- n + 1L; "computed" }
  expect_identical(compute_once(x, "probe", fn), "computed")
  expect_identical(compute_once(x, "probe", fn), "computed")
  expect_identical(n, 1L)
  ## a different object has its own cache
  y <- Constant(matrix(c(2, 1, 1, 3), 2, 2))
  expect_identical(compute_once(y, "probe", fn), "computed")
  expect_identical(n, 2L)
  ## a different member name on the same object is a miss
  expect_identical(compute_once(x, "probe2", fn), "computed")
  expect_identical(n, 3L)
})

## @cvxpy NONE
test_that("compute_once keys the cache on the active scopes, as _cache_key does", {
  x <- Constant(1)
  seen <- character(0)
  fn <- function(obj) {
    tag <- paste0(if (dpp_scope_active()) "D" else "-",
                  if (quad_form_dpp_scope_active()) "Q" else "-")
    seen <<- c(seen, tag)
    tag
  }
  expect_identical(compute_once(x, "k", fn), "--")
  expect_identical(with_dpp_scope(compute_once(x, "k", fn)), "D-")
  expect_identical(with_quad_form_dpp_scope(compute_once(x, "k", fn)), "-Q")
  expect_identical(
    with_dpp_scope(with_quad_form_dpp_scope(compute_once(x, "k", fn))), "DQ")
  ## second round: every scope state is now a hit and keeps its own answer
  expect_identical(compute_once(x, "k", fn), "--")
  expect_identical(with_dpp_scope(compute_once(x, "k", fn)), "D-")
  expect_identical(with_quad_form_dpp_scope(compute_once(x, "k", fn)), "-Q")
  expect_identical(
    with_dpp_scope(with_quad_form_dpp_scope(compute_once(x, "k", fn))), "DQ")
  expect_identical(seen, c("--", "D-", "-Q", "DQ"))
})

## @cvxpy NONE
test_that("compute_once caches FALSE and NULL results, not just truthy ones", {
  x <- Constant(1)
  n <- 0L
  fn_false <- function(obj) { n <<- n + 1L; FALSE }
  expect_false(compute_once(x, "f", fn_false))
  expect_false(compute_once(x, "f", fn_false))
  expect_identical(n, 1L)
  m <- 0L
  fn_null <- function(obj) { m <<- m + 1L; NULL }
  expect_null(compute_once(x, "g", fn_null))
  expect_null(compute_once(x, "g", fn_null))
  expect_identical(m, 1L)
})

## `scoped = FALSE` is CVXPY's plain instance attribute set once
## (`if self._x is None: self._x = ...`): one answer per object, whatever scope
## is active, and the body runs once.
## @cvxpy NONE
test_that("compute_once(scoped = FALSE) keeps one answer across scopes", {
  x <- Constant(1)
  n <- 0L
  fn <- function(obj) { n <<- n + 1L; "once" }
  expect_identical(compute_once(x, "u", fn, scoped = FALSE), "once")
  expect_identical(with_dpp_scope(compute_once(x, "u", fn, scoped = FALSE)), "once")
  expect_identical(
    with_quad_form_dpp_scope(compute_once(x, "u", fn, scoped = FALSE)), "once")
  expect_identical(
    with_dpp_scope(with_quad_form_dpp_scope(compute_once(x, "u", fn, scoped = FALSE))),
    "once")
  expect_identical(n, 1L)
  keys <- ls(x@.cache, all.names = TRUE)
  expect_true("u" %in% keys)
  expect_false(any(grepl("^u_", keys)))
})

## Constant stores its value-level tests once per instance, as CVXPY does with
## `_psd_test`, `_symm`/`_herm` and `_nonneg`/`_nonpos` (constant.py 222-290);
## only the decorated wrappers is_psd / is_hermitian are scope-keyed. Before
## 1.9.2.9144 the eigenvalue test, the Hermitian check and the sign scan ran
## once per scope (`qp` mem_alloc +16%, notes/benchmark_results.md
## "Regression B").
## @cvxpy constant.py::Constant::is_psd constant.py::Constant::is_hermitian constant.py::Constant::is_nonneg
test_that("Constant value tests are cached once, wrappers per scope", {
  x <- Constant(matrix(c(2, 1, 1, 3), 2, 2))
  expect_true(is_psd(x))
  expect_true(with_dpp_scope(is_psd(x)))
  expect_true(is_nonneg(x))
  expect_true(with_dpp_scope(is_nonneg(x)))
  expect_true(is_symmetric(x))
  expect_true(with_dpp_scope(is_symmetric(x)))
  keys <- ls(x@.cache, all.names = TRUE)
  ## the decorated wrappers: one entry per scope, as @perf.compute_once
  expect_true(all(c("is_psd", "is_psd_dpp", "is_hermitian", "is_hermitian_dpp",
                    "is_symmetric", "is_symmetric_dpp") %in% keys))
  ## the stored tests: exactly one entry each, no scope variant
  expect_true(all(c("psd_test", "symm_herm", "sign") %in% keys))
  expect_false(any(c("psd_test_dpp", "symm_herm_dpp", "sign_dpp") %in% keys))
  ## negative-definite: the NSD test is stored once too
  y <- Constant(-matrix(c(2, 1, 1, 3), 2, 2))
  expect_true(is_nsd(y))
  expect_true(with_dpp_scope(is_nsd(y)))
  keys_y <- ls(y@.cache, all.names = TRUE)
  expect_true("nsd_test" %in% keys_y)
  expect_false("nsd_test_dpp" %in% keys_y)
})

## Problem.is_dgp is @perf.compute_once in CVXPY (problem.py:315), so a DPP-scope
## query and a plain query are separate cache entries. The nested parametric
## power is DGP but not DPP: the in-scope answer must not leak out of scope.
## @cvxpy test_dpp.py::TestDgp::test_nested_power_not_dpp
test_that("Problem-level DGP cache is keyed on the DPP scope", {
  x <- Variable(pos = TRUE)
  alpha <- Parameter(pos = TRUE, value = 2)
  prob <- Problem(Minimize(power(power(x, alpha), alpha)), list(x >= 1))
  expect_false(with_dpp_scope(is_dgp(prob)))
  expect_true(is_dgp(prob))
  prob2 <- Problem(Minimize(power(power(x, alpha), alpha)), list(x >= 1))
  expect_true(is_dgp(prob2))
  expect_false(with_dpp_scope(is_dgp(prob2)))
})
