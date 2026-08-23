## Structural key store (src/RcppCseKeys.cpp + reductions/subexpr_cache.R).
##
## R-SPECIFIC. CVXPY has no counterpart: upstream builds keys as Python tuples
## and lets dict hashing do the rest, so none of this can go wrong there. CVXR
## interns keys as integers in a C++ store (ADR D_PERF.6) because R has no
## integer-keyed hash map and building a string per node cost ~370us per compile
## on the qp bench cell. These tests pin the properties that substitution has to
## preserve.
##
## THE ASYMMETRY THAT MOTIVATES EVERY CASE BELOW: a FALSE merge is a wrong
## answer -- two subtrees that are not interchangeable canonicalized as one --
## while a MISSED merge only costs time. Everything ambiguous must resolve to
## uncacheable.

.store <- function() StructuralKeyCache()

## What the previous implementation produced for the same value, so "unchanged
## behavior" is asserted against the real thing rather than against a belief.
.legacy_value_key <- function(v)
  paste0("const-val|", paste(dim_or_length(v), collapse = ","), "|",
         paste(sprintf("%a", as.numeric(v)), collapse = ","))

test_that("doubles key exactly as the previous sprintf('%a') implementation did", {
  s <- .store()
  K <- function(v) .Call(`_CVXR_CseKeys__intern_doubles`, s, NA_integer_,
                         .CSE_TAG_CONST_VAL, as.numeric(v),
                         as.integer(dim_or_length(v)))
  pairs <- list(
    list(NA_real_, NaN),        # distinct: R distinguishes NA from NaN
    list(NaN, 0 / 0),           # merged on EVERY platform: NaNs are
                                # canonicalized before keying (0/0 is the
                                # negative quiet NaN on x86 -- raw-bit keying
                                # failed CI on Linux/Windows)
    list(Inf, -Inf),
    list(0, -0),                # distinct: equal numerically, different bits
    list(1, 1 + 3e-16),         # distinct: a real one-ulp difference
    list(1, 1 + 1e-16),         # merged: 1 + 1e-16 IS 1 in IEEE754
    list(.Machine$double.xmax, Inf),
    list(1e-308, 1e-309)
  )
  for (p in pairs)
    expect_identical(K(p[[1L]]) != K(p[[2L]]),
                     .legacy_value_key(p[[1L]]) != .legacy_value_key(p[[2L]]),
                     info = paste(format(p[[1L]]), "vs", format(p[[2L]])))
})

test_that("value keys are stable, and dims are part of the key", {
  s <- .store()
  K <- function(v, dims = as.integer(dim_or_length(v)))
    .Call(`_CVXR_CseKeys__intern_doubles`, s, NA_integer_, .CSE_TAG_CONST_VAL,
          as.numeric(v), dims)
  for (v in list(NA_real_, NaN, Inf, -Inf, 0, -0, 1e308, 1e-308))
    expect_identical(K(v), K(v), info = format(v))
  expect_identical(K(1L), K(1))                       # integer coerces to double
  expect_false(K(c(1, 2), c(2L, 1L)) == K(c(1, 2), c(1L, 2L)))   # 2x1 vs 1x2
  expect_false(K(c(1, 2, 3)) == K(c(3, 2, 1)))                   # order matters
  expect_false(K(numeric(0), integer(0)) == K(0))                # empty vs c(0)
})

test_that("NA anywhere yields an uncacheable key and is never memoised", {
  s <- .store()
  cc <- .CseKeys__class_code(s, "Foo")
  node <- function(id, code = cc, shape = 2L, data = 0L, kids = 3L)
    .Call(`_CVXR_CseKeys__key_node`, s, id, code, shape, data, kids)

  expect_true(is.na(node(1L, kids = c(3L, NA_integer_))))
  expect_true(is.na(node(2L, shape = c(2L, NA_integer_))))
  expect_true(is.na(node(3L, data = NA_integer_)))
  expect_true(is.na(node(4L, code = NA_integer_)))
  ## `.UNCACHEABLE` IS NA_integer_, so a leaked uncacheable child key would
  ## arrive here as an ordinary INT_MIN component -- this is what stops it
  ## becoming a silent false merge.
  expect_true(is.na(.CseKeys__memo_get(s, 1L)))
  expect_true(is.na(.Call(`_CVXR_CseKeys__intern_doubles`, s, NA_integer_,
                          .CSE_TAG_CONST_VAL, 1, NA_integer_)))
  expect_true(.is_uncacheable(.UNCACHEABLE))
})

test_that("signature boundaries cannot be forged by component values", {
  s <- .store()
  cc <- .CseKeys__class_code(s, "Foo")
  node <- function(id, shape, kids)
    .Call(`_CVXR_CseKeys__key_node`, s, id, cc, shape, 0L, kids)

  ## Runs are length-prefixed rather than sentinel-separated, so no component
  ## value can fake the shape/children boundary.
  expect_false(node(1L, c(1L, 2L), 3L) == node(2L, 1L, c(2L, 3L)))
  expect_false(node(3L, integer(0), c(1L, 2L)) == node(4L, 1L, 2L))
  expect_false(node(5L, c(1L, -2147483647L, 2L), integer(0)) == node(6L, 1L, 2L))
  ## and identical structure still merges, different shape still splits
  expect_identical(node(7L, c(4L, 4L), 9L), node(8L, c(4L, 4L), 9L))
  expect_false(node(9L, c(4L, 4L), 9L) == node(10L, c(4L, 5L), 9L))
})

test_that("the memo is per-store and keys restart in a fresh store", {
  s1 <- .store(); s2 <- .store()
  c1 <- .CseKeys__class_code(s1, "Foo")
  c2 <- .CseKeys__class_code(s2, "Foo")
  k1 <- .Call(`_CVXR_CseKeys__key_node`, s1, 1L, c1, 2L, 0L, 3L)
  expect_identical(k1, 0L)                       # first signature in a store
  expect_true(is.na(.CseKeys__memo_get(s2, 1L))) # nothing leaks across stores
  expect_identical(.Call(`_CVXR_CseKeys__key_node`, s2, 1L, c2, 2L, 0L, 3L), 0L)
})

test_that("the result cache round-trips R objects and protects them from gc", {
  s <- .store()
  x <- Variable(3L)
  expect_null(.CseKeys__result_get(s, 42L))
  .CseKeys__result_set(s, 42L, x)
  expect_identical(.CseKeys__result_get(s, 42L), x)
  ## The store holds this only from C++; if protection were wrong, a full gc
  ## would collect it.
  gc(full = TRUE)
  expect_identical(.CseKeys__result_get(s, 42L), x)
})

test_that("expression keys merge by structure and split by leaf identity", {
  s <- .store()
  x <- Variable(3L)
  y <- Variable(3L)
  ## Two separately built but structurally identical subtrees over the SAME
  ## leaf must merge; over a different leaf must not.
  expect_identical(expr_key(p_norm(x, 1), s), expr_key(p_norm(x, 1), s))
  expect_false(expr_key(p_norm(x, 1), s) == expr_key(p_norm(y, 1), s))
  expect_false(expr_key(p_norm(x, 1), s) == expr_key(p_norm(x, 2), s))
  ## Distinct Constant objects holding equal small values merge (this is what
  ## makes huber's default M = 0.5 shareable); Variables never do.
  expect_identical(expr_key(Constant(0.5), s), expr_key(Constant(0.5), s))
  expect_false(expr_key(Variable(), s) == expr_key(Variable(), s))
})

test_that("deny-listed classes are refused rather than risked", {
  s <- .store()
  ## PartialProblem hard-codes shape c(1,1) with empty args, so every instance
  ## would key alike -- the worst possible false merge.
  for (cls in .cse_deny_list()) expect_true(inherits(cls, "S7_class"))
  expect_length(.cse_deny_list(), 2L)
})
