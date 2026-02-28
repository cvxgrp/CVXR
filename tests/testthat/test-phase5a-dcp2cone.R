## Tests for Phase 5a: Dcp2Cone end-to-end integration
## ============================================================

# ── Dcp2Cone accepts ──────────────────────────────────────────────

## @cvxpy NONE
test_that("Dcp2Cone accepts DCP Minimize problems", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(x)), list(x >= 0))
  d2c <- Dcp2Cone()
  expect_true(reduction_accepts(d2c, p))
})

## @cvxpy NONE
test_that("Dcp2Cone rejects Maximize (needs FlipObjective first)", {
  x <- Variable(2)
  p <- Problem(Maximize(sum(x)), list(x <= 1))
  d2c <- Dcp2Cone()
  expect_false(reduction_accepts(d2c, p))
})

## @cvxpy NONE
test_that("Dcp2Cone rejects non-DCP problems", {
  x <- Variable(1)
  p <- Problem(Minimize(-exp(x)))  ## concave, not convex
  d2c <- Dcp2Cone()
  ## -exp(x) is concave → Minimize(-exp(x)) requires concave objective → not DCP
  expect_false(reduction_accepts(d2c, p))
})

# ── Dcp2Cone apply: affine-only problem (passthrough) ─────────────

## @cvxpy NONE
test_that("Dcp2Cone apply on affine problem passes through", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(x)), list(x >= 0))
  d2c <- Dcp2Cone()
  result <- reduction_apply(d2c, p)
  new_p <- result[[1L]]
  inv_data <- result[[2L]]
  ## Affine objective should still be affine
  expect_true(is_affine(new_p@objective@args[[1L]]))
  expect_true(S7_inherits(inv_data, InverseData))
})

# ── Dcp2Cone apply: exp(x) ───────────────────────────────────────

## @cvxpy NONE
test_that("Dcp2Cone canonicalizes exp(x) → affine + ExpCone", {
  x <- Variable(1)
  p <- Problem(Minimize(exp(x)), list(x >= -1))
  d2c <- Dcp2Cone()
  expect_true(reduction_accepts(d2c, p))
  result <- reduction_apply(d2c, p)
  new_p <- result[[1L]]

  ## Objective expression should now be affine (it's the epigraph variable t)
  expect_true(is_affine(new_p@objective@args[[1L]]))

  ## Should have ExpCone constraint(s) in the new problem
  has_expcone <- any(vapply(new_p@constraints, function(c) S7_inherits(c, ExpCone), logical(1)))
  expect_true(has_expcone)
})

# ── Dcp2Cone apply: log(x) ───────────────────────────────────────

## @cvxpy NONE
test_that("Dcp2Cone canonicalizes -log(x) → affine + ExpCone", {
  x <- Variable(1, nonneg = TRUE)
  p <- Problem(Minimize(-log(x)), list(x >= 1))
  d2c <- Dcp2Cone()
  expect_true(reduction_accepts(d2c, p))
  result <- reduction_apply(d2c, p)
  new_p <- result[[1L]]
  expect_true(is_affine(new_p@objective@args[[1L]]))
  has_expcone <- any(vapply(new_p@constraints, function(c) S7_inherits(c, ExpCone), logical(1)))
  expect_true(has_expcone)
})

# ── Dcp2Cone apply: abs(x) ───────────────────────────────────────

## @cvxpy NONE
test_that("Dcp2Cone canonicalizes abs(x) → affine + NonNeg/Inequality", {
  x <- Variable(1)
  p <- Problem(Minimize(abs(x)))
  d2c <- Dcp2Cone()
  expect_true(reduction_accepts(d2c, p))
  result <- reduction_apply(d2c, p)
  new_p <- result[[1L]]
  ## abs replaced by epigraph variable t
  expect_true(is_affine(new_p@objective@args[[1L]]))
  ## Should have >= constraints (t >= x, t >= -x)
  expect_true(length(new_p@constraints) >= 2L)
})

# ── Dcp2Cone apply: max(x) ───────────────────────────────────────

## @cvxpy NONE
test_that("Dcp2Cone canonicalizes max_entries(x)", {
  x <- Variable(3)
  p <- Problem(Minimize(max_entries(x)), list(x >= -1, x <= 1))
  d2c <- Dcp2Cone()
  expect_true(reduction_accepts(d2c, p))
  result <- reduction_apply(d2c, p)
  new_p <- result[[1L]]
  expect_true(is_affine(new_p@objective@args[[1L]]))
})

# ── Dcp2Cone apply: compound expression ──────────────────────────

## @cvxpy NONE
test_that("Dcp2Cone canonicalizes exp(x) + abs(y)", {
  x <- Variable(1)
  y <- Variable(1)
  p <- Problem(Minimize(exp(x) + abs(y)), list(x >= 0, y >= -1))
  d2c <- Dcp2Cone()
  expect_true(reduction_accepts(d2c, p))
  result <- reduction_apply(d2c, p)
  new_p <- result[[1L]]
  ## After canonicalization, objective args should all be affine
  expect_true(is_affine(new_p@objective@args[[1L]]))
})

# ── FlipObjective + Dcp2Cone pipeline ────────────────────────────

## @cvxpy NONE
test_that("FlipObjective + Dcp2Cone handles Maximize problems", {
  x <- Variable(1)
  p <- Problem(Maximize(log(x)), list(x >= 1, x <= 10))

  ## Step 1: Flip
  fo <- FlipObjective()
  r1 <- reduction_apply(fo, p)
  flipped <- r1[[1L]]
  expect_true(S7_inherits(flipped@objective, Minimize))

  ## flipped minimizes -log(x) which is convex → DCP
  expect_true(is_dcp(flipped))

  ## Step 2: Dcp2Cone
  d2c <- Dcp2Cone()
  r2 <- reduction_apply(d2c, flipped)
  conic <- r2[[1L]]
  expect_true(is_affine(conic@objective@args[[1L]]))
})

# ── Constraint canonicalization ──────────────────────────────────

## @cvxpy NONE
test_that("Dcp2Cone canonicalizes constraint args", {
  x <- Variable(1)
  y <- Variable(1)
  ## exp(x) <= y is a valid DCP constraint (convex <= affine)
  p <- Problem(Minimize(y), list(exp(x) <= y, x >= -2))
  d2c <- Dcp2Cone()
  expect_true(reduction_accepts(d2c, p))
  result <- reduction_apply(d2c, p)
  new_p <- result[[1L]]
  ## The exp constraint should be replaced by ExpCone
  has_expcone <- any(vapply(new_p@constraints, function(c) S7_inherits(c, ExpCone), logical(1)))
  expect_true(has_expcone)
})

# ── InverseData constraint mapping ──────────────────────────────

## @cvxpy NONE
test_that("InverseData records constraint ID mapping", {
  x <- Variable(2)
  cons <- list(x >= 0, x <= 5)
  p <- Problem(Minimize(sum(x)), cons)
  d2c <- Dcp2Cone()
  result <- reduction_apply(d2c, p)
  inv <- result[[2L]]
  ## cons_id_map should have entries for original constraints
  ids_stored <- ls(inv@cons_id_map)
  expect_equal(length(ids_stored), 2L)
})

# ── Dcp2Cone: norm1 ────────────────────────────────────────────

## @cvxpy NONE
test_that("Dcp2Cone canonicalizes norm1(x)", {
  x <- Variable(3)
  p <- Problem(Minimize(norm1(x)), list(x >= -1, x <= 1))
  d2c <- Dcp2Cone()
  expect_true(reduction_accepts(d2c, p))
  result <- reduction_apply(d2c, p)
  new_p <- result[[1L]]
  expect_true(is_affine(new_p@objective@args[[1L]]))
})

# ── Dcp2Cone: norm_inf ─────────────────────────────────────────

## @cvxpy NONE
test_that("Dcp2Cone canonicalizes norm_inf(x)", {
  x <- Variable(3)
  p <- Problem(Minimize(norm_inf(x)), list(x >= -1, x <= 1))
  d2c <- Dcp2Cone()
  expect_true(reduction_accepts(d2c, p))
  result <- reduction_apply(d2c, p)
  new_p <- result[[1L]]
  expect_true(is_affine(new_p@objective@args[[1L]]))
})

# ── Dcp2Cone: sum_largest ──────────────────────────────────────

## @cvxpy NONE
test_that("Dcp2Cone canonicalizes sum_largest(x, k)", {
  x <- Variable(5)
  p <- Problem(Minimize(sum_largest(x, 2)), list(x >= -1, x <= 1))
  d2c <- Dcp2Cone()
  expect_true(reduction_accepts(d2c, p))
  result <- reduction_apply(d2c, p)
  new_p <- result[[1L]]
  expect_true(is_affine(new_p@objective@args[[1L]]))
})

# ── reduction_reduce / reduction_retrieve ──────────────────────

## @cvxpy NONE
test_that("reduction_reduce caches result", {
  x <- Variable(2)
  p <- Problem(Minimize(sum(x)), list(x >= 0))
  fo <- FlipObjective()
  new_p <- reduction_reduce(fo, p)
  expect_true(S7_inherits(new_p, Problem))
  ## Second call returns cached
  new_p2 <- reduction_reduce(fo, p)
  expect_identical(new_p, new_p2)
})
