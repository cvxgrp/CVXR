## Tests for Phase 7b-i: power_tools.R — Exact Rational Arithmetic Utilities
## CVXPY SOURCE: utilities/power_tools.py
## Cross-validated against CVXPY output where noted.
## Note: gmp functions (as.bigq, as.bigz, etc.) are available via CVXR's importFrom

# ══════════════════════════════════════════════════════════════════════
# next_pow2
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("next_pow2 returns correct values", {
  ## CVXPY-verified
  expect_equal(next_pow2(0L), 1L)
  expect_equal(next_pow2(1L), 1L)
  expect_equal(next_pow2(3L), 4L)
  expect_equal(next_pow2(8L), 8L)
  expect_equal(next_pow2(1023L), 1024L)
  expect_equal(next_pow2(5L), 8L)
  expect_equal(next_pow2(16L), 16L)
  expect_equal(next_pow2(17L), 32L)
})

# ══════════════════════════════════════════════════════════════════════
# is_power2
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("is_power2 identifies powers of 2 correctly", {
  expect_true(is_power2(1L))
  expect_true(is_power2(2L))
  expect_true(is_power2(4L))
  expect_true(is_power2(1024L))
  expect_true(is_power2(2^10))

  expect_false(is_power2(0L))
  expect_false(is_power2(-4L))
  expect_false(is_power2(3L))
  expect_false(is_power2(6L))
  ## Note: R doesn't distinguish int vs float like Python, so is_power2(1.0) is TRUE in R
  ## (Python's is_power2(1.0) is False because isinstance(1.0, Integral) is False)
})

## @cvxpy NONE
test_that("is_power2 handles bigz input", {
  expect_true(is_power2(as.bigz(4L)))
  expect_true(is_power2(as.bigz(1024L)))
  expect_false(is_power2(as.bigz(0L)))
  expect_false(is_power2(as.bigz(3L)))
})

# ══════════════════════════════════════════════════════════════════════
# .limit_denominator (CPython algorithm)
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that(".limit_denominator matches CPython Fraction.limit_denominator", {
  ## CVXPY-verified: Fraction(7,16).limit_denominator(5) = 2/5
  result <- .limit_denominator(as.bigq(7L, 16L), as.bigz(5L))
  expect_equal(result, as.bigq(2L, 5L))

  ## Fraction(1,3).limit_denominator(1) = 0/1
  result <- .limit_denominator(as.bigq(1L, 3L), as.bigz(1L))
  expect_equal(result, as.bigq(0L, 1L))

  ## Fraction(1,3).limit_denominator(10) = 1/3 (already within limit)
  result <- .limit_denominator(as.bigq(1L, 3L), as.bigz(10L))
  expect_equal(result, as.bigq(1L, 3L))

  ## Fraction(355,113).limit_denominator(100) = 311/99 (pi approx)
  result <- .limit_denominator(as.bigq(355L, 113L), as.bigz(100L))
  expect_equal(result, as.bigq(311L, 99L))
})

## @cvxpy NONE
test_that(".limit_denominator returns unchanged fraction when denom within limit", {
  result <- .limit_denominator(as.bigq(1L, 4L), as.bigz(4L))
  expect_equal(result, as.bigq(1L, 4L))

  result <- .limit_denominator(as.bigq(3L, 7L), as.bigz(100L))
  expect_equal(result, as.bigq(3L, 7L))
})

## @cvxpy NONE
test_that(".limit_denominator handles edge cases", {
  ## Integer fraction (denominator = 1)
  result <- .limit_denominator(as.bigq(5L), as.bigz(10L))
  expect_equal(result, as.bigq(5L))

  ## Zero
  result <- .limit_denominator(as.bigq(0L), as.bigz(5L))
  expect_equal(result, as.bigq(0L))

  ## max_denom < 1 should error
  expect_error(.limit_denominator(as.bigq(1L, 3L), as.bigz(0L)))
})

## @cvxpy NONE
test_that(".limit_denominator handles negative fractions", {
  result <- .limit_denominator(as.bigq(-7L, 16L), as.bigz(5L))
  expect_equal(result, as.bigq(-2L, 5L))
})

## @cvxpy NONE
test_that(".limit_denominator semiconvergent is better than last convergent", {
  ## Key test: 7/16 with max_denom=5
  ## Last convergent would be 1/2 but semiconvergent 2/5 is closer
  result <- .limit_denominator(as.bigq(7L, 16L), as.bigz(5L))
  ## 7/16 = 0.4375
  ## 1/2 = 0.5, |diff| = 0.0625
  ## 2/5 = 0.4, |diff| = 0.0375  <- closer!
  expect_equal(result, as.bigq(2L, 5L))
})

# ══════════════════════════════════════════════════════════════════════
# is_dyad
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("is_dyad identifies dyadic fractions", {
  expect_true(is_dyad(as.bigq(1L, 4L)))
  expect_true(is_dyad(as.bigq(0L)))
  expect_true(is_dyad(as.bigq(1L)))
  expect_true(is_dyad(0L))
  expect_true(is_dyad(1L))

  expect_false(is_dyad(as.bigq(1L, 3L)))
  expect_false(is_dyad(as.bigq(-1L, 4L)))
  expect_false(is_dyad(as.bigq(1L, 6L)))
})

# ══════════════════════════════════════════════════════════════════════
# is_weight / is_dyad_weight
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("is_weight validates weight vectors", {
  expect_true(is_weight(as.bigq(c(1L, 2L), c(3L, 3L))))
  expect_true(is_weight(c(0L, 1L, 0L)))
  expect_true(is_weight(as.bigq(c(1L, 1L), c(2L, 2L))))

  ## Doesn't sum to 1
  expect_false(is_weight(as.bigq(c(2L, 2L), c(3L, 3L))))

  ## Float inputs: 0.1 + 0.9 in IEEE 754 bigq != 1 exactly
  expect_false(is_weight(c(0.1, 0.9)))

  ## Negative element
  expect_false(is_weight(as.bigq(c(-1L, 2L), c(1L, 1L))))
})

## @cvxpy NONE
test_that("is_dyad_weight validates dyadic weight vectors", {
  expect_true(is_dyad_weight(as.bigq(c(1L, 1L), c(2L, 2L))))
  expect_true(is_dyad_weight(c(0L, 1L, 0L)))

  ## Sum to 1 but not dyadic (1/3, 2/3)
  expect_false(is_dyad_weight(as.bigq(c(1L, 2L), c(3L, 3L))))
})

# ══════════════════════════════════════════════════════════════════════
# make_frac
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("make_frac matches CVXPY docstring examples", {
  ## CVXPY-verified: make_frac([.123,.345,.532], 10) = (1/10, 2/5, 1/2)
  result <- make_frac(c(.123, .345, .532), 10L)
  expect_equal(result, as.bigq(c(1L, 4L, 5L), c(10L, 10L, 10L)))

  ## CVXPY-verified: make_frac([.123,.345,.532], 100) = (3/25, 7/20, 53/100)
  result <- make_frac(c(.123, .345, .532), 100L)
  expect_equal(result, as.bigq(c(12L, 35L, 53L), c(100L, 100L, 100L)))
})

## @cvxpy NONE
test_that("make_frac results sum to 1", {
  result <- make_frac(c(1, 2, 3), 8L)
  expect_equal(sum(result), as.bigq(1L))

  result <- make_frac(c(.123, .345, .532), 1024L)
  expect_equal(sum(result), as.bigq(1L))
})

# ══════════════════════════════════════════════════════════════════════
# dyad_completion
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("dyad_completion matches CVXPY docstring examples", {
  ## CVXPY-verified: (1/3, 1/3, 1/3) -> (1/4, 1/4, 1/4, 1/4)
  w <- as.bigq(c(1L, 1L, 1L), c(3L, 3L, 3L))
  result <- dyad_completion(w)
  expected <- as.bigq(c(1L, 1L, 1L, 1L), c(4L, 4L, 4L, 4L))
  expect_equal(result, expected)

  ## CVXPY-verified: (1/3, 1/5, 7/15) -> (5/16, 3/16, 7/16, 1/16)
  w <- as.bigq(c(1L, 1L, 7L), c(3L, 5L, 15L))
  result <- dyad_completion(w)
  expected <- as.bigq(c(5L, 3L, 7L, 1L), c(16L, 16L, 16L, 16L))
  expect_equal(result, expected)
})

## @cvxpy NONE
test_that("dyad_completion returns dyadic input unchanged", {
  w <- as.bigq(c(1L, 0L, 0L), c(1L, 1L, 1L))
  result <- dyad_completion(w)
  expect_equal(result, w)

  w <- as.bigq(c(1L, 1L), c(2L, 2L))
  result <- dyad_completion(w)
  expect_equal(result, w)

  w <- as.bigq(c(1L, 3L), c(4L, 4L))
  result <- dyad_completion(w)
  expect_equal(result, w)
})

## @cvxpy NONE
test_that("dyad_completion result is valid dyadic weight", {
  w <- as.bigq(c(2L, 3L), c(5L, 5L))
  result <- dyad_completion(w)
  expect_true(is_dyad_weight(result))
  expect_equal(sum(result), as.bigq(1L))
})

# ══════════════════════════════════════════════════════════════════════
# fracify
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("fracify matches CVXPY docstring: integer input [1,2,3]", {
  ## CVXPY-verified
  result <- fracify(c(1, 2, 3))
  expect_equal(result$w, as.bigq(c(1L, 2L, 3L), c(6L, 6L, 6L)))
  expect_equal(result$w_dyad,
               as.bigq(c(1L, 2L, 3L, 2L), c(8L, 8L, 8L, 8L)))
  expect_true(is_weight(result$w))
  expect_true(is_dyad_weight(result$w_dyad))
})

## @cvxpy NONE
test_that("fracify matches CVXPY docstring: five ones", {
  ## CVXPY-verified
  result <- fracify(c(1, 1, 1, 1, 1))
  expect_equal(result$w, as.bigq(c(1L, 1L, 1L, 1L, 1L), c(5L, 5L, 5L, 5L, 5L)))
  expect_equal(result$w_dyad,
               as.bigq(c(1L, 1L, 1L, 1L, 1L, 3L), c(8L, 8L, 8L, 8L, 8L, 8L)))
})

## @cvxpy NONE
test_that("fracify matches CVXPY docstring: float input", {
  ## CVXPY-verified: fracify([.23,.56,.87])
  result <- fracify(c(.23, .56, .87))
  expect_equal(result$w, as.bigq(c(23L, 56L, 87L), c(166L, 166L, 166L)))
  expect_true(is_dyad_weight(result$w_dyad))
})

## @cvxpy NONE
test_that("fracify handles standard basis unit vector", {
  ## CVXPY-verified
  result <- fracify(c(0, 0, 1))
  expect_equal(result$w, as.bigq(c(0L, 0L, 1L), c(1L, 1L, 1L)))
  expect_equal(result$w_dyad, as.bigq(c(0L, 0L, 1L), c(1L, 1L, 1L)))
})

## @cvxpy NONE
test_that("fracify already-dyadic input returns itself", {
  ## CVXPY-verified
  a <- as.bigq(c(1L, 1L, 3L), c(2L, 8L, 8L))
  result <- fracify(a)
  expect_equal(result$w, a)
  expect_equal(result$w_dyad, a)
})

## @cvxpy NONE
test_that("fracify with force_dyad", {
  ## CVXPY-verified: fracify([3.4, 8, 3/2], force_dyad=True)
  result <- fracify(c(3.4, 8, 1.5), force_dyad = TRUE)
  expect_equal(result$w, result$w_dyad)  # force_dyad means w == w_dyad
  expect_true(is_dyad_weight(result$w))
  expect_equal(sum(result$w), as.bigq(1L))
})

## @cvxpy NONE
test_that("fracify rejects negative inputs", {
  expect_error(fracify(c(1, -2, 3)), "nonnegative")
})

## @cvxpy NONE
test_that("fracify result components are valid", {
  for (a in list(c(1, 2, 3), c(1, 1), c(3, 5, 7, 11))) {
    result <- fracify(a)
    expect_true(is_weight(result$w))
    expect_true(is_dyad_weight(result$w_dyad))
    expect_true(check_dyad(result$w, result$w_dyad))
  }
})

# ══════════════════════════════════════════════════════════════════════
# .split_dyad
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that(".split_dyad splits (1/2, 1/2) correctly", {
  w <- as.bigq(c(1L, 1L), c(2L, 2L))
  ch <- .split_dyad(w)
  expect_length(ch, 2L)
  ## Children should each be valid dyadic weight vectors
  expect_true(is_dyad_weight(ch[[1L]]))
  expect_true(is_dyad_weight(ch[[2L]]))
  ## d_tup = 1/2 * (child1 + child2)
  expect_equal(w, (ch[[1L]] + ch[[2L]]) / as.bigq(2L))
})

## @cvxpy NONE
test_that(".split_dyad returns empty for basis vectors", {
  w1 <- as.bigq(c(1L, 0L), c(1L, 1L))
  expect_equal(.split_dyad(w1), list())

  w2 <- as.bigq(c(0L, 1L, 0L), c(1L, 1L, 1L))
  expect_equal(.split_dyad(w2), list())

  w3 <- as.bigq(c(0L, 0L, 0L, 1L), c(1L, 1L, 1L, 1L))
  expect_equal(.split_dyad(w3), list())
})

## @cvxpy NONE
test_that(".split_dyad children satisfy parent = 1/2*(child1 + child2)", {
  w <- as.bigq(c(1L, 1L, 1L, 1L), c(4L, 4L, 4L, 4L))
  ch <- .split_dyad(w)
  expect_length(ch, 2L)
  expect_equal(w, (ch[[1L]] + ch[[2L]]) / as.bigq(2L))
})

# ══════════════════════════════════════════════════════════════════════
# decompose
# ══════════════════════════════════════════════════════════════════════
## @cvxpy test_power_tools.py::TestGeoMean::test_multi_step_dyad_completion
test_that("decompose (1/2, 1/2) gives 3 nodes", {
  ## CVXPY-verified: 3 nodes (root + 2 basis vectors)
  w <- as.bigq(c(1L, 1L), c(2L, 2L))
  tree <- decompose(w)
  expect_equal(length(tree$keys), 3L)
  ## Root has 2 children, basis vectors have 0
  expect_length(tree$children[[1L]], 2L)
  expect_length(tree$children[[2L]], 0L)
  expect_length(tree$children[[3L]], 0L)
})

## @cvxpy test_power_tools.py::TestGeoMean::test_multi_step_dyad_completion
test_that("decompose (1/4, 1/4, 1/4, 1/4) gives 7 nodes", {
  ## CVXPY-verified
  w <- as.bigq(c(1L, 1L, 1L, 1L), c(4L, 4L, 4L, 4L))
  tree <- decompose(w)
  expect_equal(length(tree$keys), 7L)
})

## @cvxpy NONE
test_that("decompose (1/8, 1/4, 3/8, 1/4) gives 8 nodes", {
  ## CVXPY-verified
  w <- as.bigq(c(1L, 2L, 3L, 2L), c(8L, 8L, 8L, 8L))
  tree <- decompose(w)
  expect_equal(length(tree$keys), 8L)
})

## @cvxpy NONE
test_that("decompose root is first key", {
  w <- as.bigq(c(1L, 1L, 1L, 1L), c(4L, 4L, 4L, 4L))
  tree <- decompose(w)
  expect_equal(tree$keys[[1L]], w)
})

## @cvxpy NONE
test_that("decompose rejects non-dyadic input", {
  w <- as.bigq(c(1L, 2L), c(3L, 3L))
  expect_error(decompose(w), "dyadic weight")
})

## @cvxpy NONE
test_that("decompose handles basis vector input", {
  w <- as.bigq(c(0L, 1L, 0L), c(1L, 1L, 1L))
  tree <- decompose(w)
  expect_equal(length(tree$keys), 1L)
  expect_length(tree$children[[1L]], 0L)
})

# ══════════════════════════════════════════════════════════════════════
# check_dyad
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("check_dyad validates correct completions", {
  ## CVXPY docstring: (1/3,1/3,1/3) with (1/4,1/4,1/4,1/4)
  w <- as.bigq(c(1L, 1L, 1L), c(3L, 3L, 3L))
  w_dyad <- as.bigq(c(1L, 1L, 1L, 1L), c(4L, 4L, 4L, 4L))
  expect_true(check_dyad(w, w_dyad))
})

## @cvxpy NONE
test_that("check_dyad accepts w == w_dyad when already dyadic", {
  w <- as.bigq(c(1L, 0L, 3L), c(4L, 1L, 4L))
  expect_true(check_dyad(w, w))

  ## Integer input
  expect_true(check_dyad(c(1L, 0L, 0L), c(1L, 0L, 0L)))
})

## @cvxpy NONE
test_that("check_dyad rejects invalid completions", {
  ## w doesn't sum to 1
  w <- as.bigq(c(2L, 1L), c(3L, 1L))
  expect_false(check_dyad(w, w))

  ## Wrong dyadic completion
  w <- as.bigq(c(2L, 3L), c(5L, 5L))
  w_dyad <- as.bigq(c(3L, 4L, 1L), c(8L, 8L, 8L))
  expect_false(check_dyad(w, w_dyad))

  ## Correct completion for (2/5, 3/5)
  w_dyad_correct <- as.bigq(c(2L, 3L, 3L), c(8L, 8L, 8L))
  expect_true(check_dyad(w, w_dyad_correct))
})

# ══════════════════════════════════════════════════════════════════════
# approx_error
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("approx_error is zero for exact weights", {
  e <- approx_error(c(1, 1, 1), as.bigq(c(1L, 1L, 1L), c(3L, 3L, 3L)))
  expect_true(e <= 1e-10)
})

## @cvxpy NONE
test_that("approx_error is positive for inexact weights", {
  e <- approx_error(c(1, 2, 3), as.bigq(c(1L, 1L, 1L), c(3L, 3L, 3L)))
  expect_true(e > 0)
})

# ══════════════════════════════════════════════════════════════════════
# lower_bound / over_bound
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("lower_bound matches CVXPY docstring", {
  expect_equal(lower_bound(c(0L, 1L, 0L)), 0L)
  expect_equal(lower_bound(as.bigq(c(1L, 1L), c(2L, 2L))), 1L)
  expect_equal(lower_bound(as.bigq(c(1L, 1L, 1L, 1L), c(4L, 4L, 4L, 4L))), 3L)
  expect_equal(lower_bound(as.bigq(c(1L, 7L), c(8L, 8L))), 3L)
})

## @cvxpy NONE
test_that("over_bound is nonneg for decompose output", {
  w <- as.bigq(c(1L, 1L), c(2L, 2L))
  tree <- decompose(w)
  expect_true(over_bound(w, tree) >= 0L)

  w <- as.bigq(c(1L, 1L, 1L, 1L), c(4L, 4L, 4L, 4L))
  tree <- decompose(w)
  expect_true(over_bound(w, tree) >= 0L)
})

# ══════════════════════════════════════════════════════════════════════
# prettytuple / prettydict
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("prettytuple formats correctly", {
  tup <- as.bigq(c(1L, 1L), c(2L, 2L))
  expect_equal(prettytuple(tup), "(1/2, 1/2)")

  tup <- as.bigq(c(0L, 1L, 0L), c(1L, 1L, 1L))
  expect_equal(prettytuple(tup), "(0, 1, 0)")
})

## @cvxpy NONE
test_that("prettydict produces non-empty output", {
  w <- as.bigq(c(1L, 1L), c(2L, 2L))
  tree <- decompose(w)
  result <- prettydict(tree)
  expect_true(nchar(result) > 0)
  expect_true(grepl("1/2", result))
})

# ══════════════════════════════════════════════════════════════════════
# pow_high
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("pow_high(2) matches CVXPY: p=2, w=(1/2, 1/2)", {
  result <- pow_high(2)
  expect_equal(result$p, 2L)  # integer since 1/(1/2) = 2
  expect_equal(result$w, as.bigq(c(1L, 1L), c(2L, 2L)))
})

## @cvxpy NONE
test_that("pow_high(3) matches CVXPY: p=3, w=(1/3, 2/3)", {
  result <- pow_high(3)
  expect_equal(result$p, 3L)  # integer
  expect_equal(result$w, as.bigq(c(1L, 2L), c(3L, 3L)))
})

## @cvxpy NONE
test_that("pow_high(2.5) matches CVXPY: p=5/2, w=(2/5, 3/5)", {
  result <- pow_high(2.5)
  expect_equal(result$p, as.bigq(5L, 2L))
  expect_equal(result$w, as.bigq(c(2L, 3L), c(5L, 5L)))
})

## @cvxpy NONE
test_that("pow_high approx=FALSE returns plain numerics", {
  result <- pow_high(2, approx = FALSE)
  expect_equal(result$p, 2)
  expect_equal(result$w, c(0.5, 0.5))
  expect_true(is.numeric(result$p))
  expect_true(is.numeric(result$w))
})

## @cvxpy NONE
test_that("pow_high rejects p <= 1", {
  expect_error(pow_high(1), "p > 1")
  expect_error(pow_high(0.5), "p > 1")
})

# ══════════════════════════════════════════════════════════════════════
# pow_mid
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("pow_mid(0.5) matches CVXPY", {
  result <- pow_mid(0.5)
  expect_equal(result$p, as.bigq(1L, 2L))
  expect_equal(result$w, as.bigq(c(1L, 1L), c(2L, 2L)))
})

## @cvxpy NONE
test_that("pow_mid(0.3) matches CVXPY", {
  result <- pow_mid(0.3)
  expect_equal(result$p, as.bigq(3L, 10L))
  expect_equal(result$w, as.bigq(c(3L, 7L), c(10L, 10L)))
})

## @cvxpy NONE
test_that("pow_mid approx=FALSE returns plain numerics", {
  result <- pow_mid(0.5, approx = FALSE)
  expect_equal(result$p, 0.5)
  expect_equal(result$w, c(0.5, 0.5))
  expect_true(is.numeric(result$p))
})

## @cvxpy NONE
test_that("pow_mid rejects invalid p", {
  expect_error(pow_mid(0), "0 < p < 1")
  expect_error(pow_mid(1), "0 < p < 1")
  expect_error(pow_mid(1.5), "0 < p < 1")
})

# ══════════════════════════════════════════════════════════════════════
# pow_neg
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("pow_neg(-1) matches CVXPY: p=-1, w=(1/2, 1/2)", {
  result <- pow_neg(-1)
  expect_equal(result$p, as.bigq(-1L))
  expect_equal(result$w, as.bigq(c(1L, 1L), c(2L, 2L)))
})

## @cvxpy NONE
test_that("pow_neg(-2) matches CVXPY: p=-2, w=(2/3, 1/3)", {
  result <- pow_neg(-2)
  expect_equal(result$p, as.bigq(-2L))
  expect_equal(result$w, as.bigq(c(2L, 1L), c(3L, 3L)))
})

## @cvxpy NONE
test_that("pow_neg approx=FALSE returns plain numerics", {
  result <- pow_neg(-1, approx = FALSE)
  expect_equal(result$p, -1)
  expect_equal(result$w, c(0.5, 0.5))
  expect_true(is.numeric(result$p))
})

## @cvxpy NONE
test_that("pow_neg rejects p >= 0", {
  expect_error(pow_neg(0), "p < 0")
  expect_error(pow_neg(1), "p < 0")
})

# ══════════════════════════════════════════════════════════════════════
# pow_high/mid/neg weights sum correctly
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("pow_* weights sum to 1 in approx mode", {
  r <- pow_high(2); expect_equal(sum(r$w), as.bigq(1L))
  r <- pow_high(3); expect_equal(sum(r$w), as.bigq(1L))
  r <- pow_high(7); expect_equal(sum(r$w), as.bigq(1L))
  r <- pow_mid(0.3); expect_equal(sum(r$w), as.bigq(1L))
  r <- pow_mid(0.7); expect_equal(sum(r$w), as.bigq(1L))
  r <- pow_neg(-1); expect_equal(sum(r$w), as.bigq(1L))
  r <- pow_neg(-0.5); expect_equal(sum(r$w), as.bigq(1L))
})

## @cvxpy NONE
test_that("pow_* weights sum to 1 in non-approx mode", {
  r <- pow_high(2, approx = FALSE); expect_equal(sum(r$w), 1)
  r <- pow_mid(0.5, approx = FALSE); expect_equal(sum(r$w), 1)
  r <- pow_neg(-1, approx = FALSE); expect_equal(sum(r$w), 1)
})

# ══════════════════════════════════════════════════════════════════════
# .bigq_vec helper
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that(".bigq_vec safely concatenates mixed types", {
  result <- .bigq_vec(as.bigq(1L, 2L), as.bigq(1L, 3L))
  expect_true(is.bigq(result))
  expect_equal(length(result), 2L)

  ## Mixing integer and bigq
  result <- .bigq_vec(3L, as.bigq(1L, 2L))
  expect_true(is.bigq(result))
  expect_equal(result[1], as.bigq(3L))
  expect_equal(result[2], as.bigq(1L, 2L))
})

# ══════════════════════════════════════════════════════════════════════
# .bigq_to_key
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that(".bigq_to_key produces consistent keys", {
  w <- as.bigq(c(1L, 1L), c(2L, 2L))
  k1 <- .bigq_to_key(w)
  k2 <- .bigq_to_key(w)
  expect_equal(k1, k2)

  ## Different vectors produce different keys
  w2 <- as.bigq(c(1L, 1L), c(3L, 3L))
  expect_false(.bigq_to_key(w) == .bigq_to_key(w2))
})

# ══════════════════════════════════════════════════════════════════════
# gm — SOC constraint factory
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("gm creates valid SOC constraint", {
  t_var <- Variable(c(1L, 1L))
  x <- Variable(c(1L, 1L))
  y <- Variable(c(1L, 1L))
  constr <- gm(t_var, x, y)
  expect_true(inherits(constr, "CVXR::SOC"))
})

## @cvxpy NONE
test_that("gm creates SOC with correct shapes for vectors", {
  t_var <- Variable(c(3L, 1L))
  x <- Variable(c(3L, 1L))
  y <- Variable(c(3L, 1L))
  constr <- gm(t_var, x, y)
  expect_true(inherits(constr, "CVXR::SOC"))
  ## SOC shape should match t_var shape
  expect_equal(constr@shape, c(3L, 1L))
})

# ══════════════════════════════════════════════════════════════════════
# gm_constrs — full constraint builder
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("gm_constrs with (1/2, 1/2) produces 1 SOC constraint", {
  t_var <- Variable(c(1L, 1L))
  x1 <- Variable(c(1L, 1L))
  x2 <- Variable(c(1L, 1L))
  p <- as.bigq(c(1L, 1L), c(2L, 2L))

  constrs <- gm_constrs(t_var, list(x1, x2), p)
  ## decompose(1/2,1/2) has 1 interior node -> 1 SOC constraint
  soc_count <- sum(vapply(constrs, function(c) inherits(c, "CVXR::SOC"), logical(1)))
  expect_equal(soc_count, 1L)
})

## @cvxpy NONE
test_that("gm_constrs with (1/4,1/4,1/4,1/4) produces 3 SOC constraints", {
  t_var <- Variable(c(1L, 1L))
  x_list <- lapply(1:4, function(i) Variable(c(1L, 1L)))
  p <- as.bigq(c(1L, 1L, 1L, 1L), c(4L, 4L, 4L, 4L))

  constrs <- gm_constrs(t_var, x_list, p)
  soc_count <- sum(vapply(constrs, function(c) inherits(c, "CVXR::SOC"), logical(1)))
  expect_equal(soc_count, 3L)
})

## @cvxpy NONE
test_that("gm_constrs single variable returns equality", {
  t_var <- Variable(c(1L, 1L))
  x1 <- Variable(c(1L, 1L))
  p <- as.bigq(1L)

  constrs <- gm_constrs(t_var, list(x1), p)
  ## Single-variable: should produce t == x constraint
  expect_true(length(constrs) >= 1L)
})

## @cvxpy NONE
test_that("gm_constrs with non-dyadic weights works (uses dyad_completion)", {
  t_var <- Variable(c(1L, 1L))
  x1 <- Variable(c(1L, 1L))
  x2 <- Variable(c(1L, 1L))
  x3 <- Variable(c(1L, 1L))
  p <- as.bigq(c(1L, 1L, 1L), c(3L, 3L, 3L))

  constrs <- gm_constrs(t_var, list(x1, x2, x3), p)
  ## Should produce SOC constraints (dyad_completion adds dummy variable)
  soc_count <- sum(vapply(constrs, function(c) inherits(c, "CVXR::SOC"), logical(1)))
  expect_true(soc_count >= 1L)
})

## @cvxpy NONE
test_that("gm_constrs rejects invalid weight vector", {
  t_var <- Variable(c(1L, 1L))
  x1 <- Variable(c(1L, 1L))
  p <- as.bigq(c(1L, 1L), c(3L, 3L))  # sum = 2/3 != 1

  expect_error(gm_constrs(t_var, list(x1), p), "sum to 1")
})

# ══════════════════════════════════════════════════════════════════════
# powcone_constrs
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("powcone_constrs creates PowCone3D constraint", {
  t_var <- Variable(c(1L, 1L))
  x1 <- Variable(c(1L, 1L))
  x2 <- Variable(c(1L, 1L))
  constrs <- powcone_constrs(t_var, list(x1, x2), 0.5)
  expect_length(constrs, 1L)
  expect_true(inherits(constrs[[1L]], "CVXR::PowCone3D"))
})

# ══════════════════════════════════════════════════════════════════════
# get_max_denom
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("get_max_denom finds maximum denominator", {
  tup <- as.bigq(c(1L, 3L), c(4L, 8L))
  expect_equal(get_max_denom(tup), 8L)

  tup <- as.bigq(c(1L, 0L), c(1L, 1L))
  expect_equal(get_max_denom(tup), 1L)
})

# ══════════════════════════════════════════════════════════════════════
# Integration: fracify -> decompose -> constraint count
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("fracify + decompose pipeline works end-to-end", {
  result <- fracify(c(1, 2, 3))
  tree <- decompose(result$w_dyad)
  expect_true(length(tree$keys) >= 1L)

  ## Lower bound should be achievable
  lb <- lower_bound(result$w_dyad)
  ob <- over_bound(result$w_dyad, tree)
  expect_true(lb >= 0L)
  expect_true(ob >= 0L)
})

## @cvxpy NONE
test_that("fracify + check_dyad round-trip", {
  for (a in list(c(1, 2, 3), c(1, 1, 1, 1, 1), c(3, 7))) {
    result <- fracify(a)
    expect_true(check_dyad(result$w, result$w_dyad))
  }
})

# ══════════════════════════════════════════════════════════════════════
# CVXPY cross-validation via Python
# ══════════════════════════════════════════════════════════════════════
## @cvxpy NONE
test_that("limit_denominator matches Python for float-originated fractions", {
  ## as.bigq(0.3) should match Python Fraction(0.3).limit_denominator(10)
  ## Python: Fraction(0.3).limit_denominator(10) = Fraction(3, 10)
  result <- .limit_denominator(as.bigq(0.3), as.bigz(10L))
  expect_equal(result, as.bigq(3L, 10L))
})

## @cvxpy NONE
test_that("decompose node count matches CVXPY for various inputs", {
  ## (1/2, 1/2) -> 3 nodes
  w <- as.bigq(c(1L, 1L), c(2L, 2L))
  expect_equal(length(decompose(w)$keys), 3L)

  ## (1/4, 1/4, 1/4, 1/4) -> 7 nodes
  w <- as.bigq(c(1L, 1L, 1L, 1L), c(4L, 4L, 4L, 4L))
  expect_equal(length(decompose(w)$keys), 7L)

  ## (1/8, 1/4, 3/8, 1/4) -> 8 nodes (fracify([1,2,3]) w_dyad)
  w <- as.bigq(c(1L, 2L, 3L, 2L), c(8L, 8L, 8L, 8L))
  expect_equal(length(decompose(w)$keys), 8L)
})
