## Partial `boolean` / `integer` index lists -- the attribute must actually
## reach the solver.
##
## Before this was fixed, `Variable(2, boolean = c(1))` looked continuous to
## `is_mixed_integer()` and contributed nothing to `bool_idx`, so the problem
## was solved as its continuous relaxation and reported "optimal" -- a silent
## wrong answer, the same shape as the two transposed-dual bugs fixed earlier
## on this branch.
##
## Expected values are CVXPY 1.9.2's, measured with
## `uv run --with "cvxpy==1.9.2"`, using its index spelling (a numpy
## multi-index) which differs from CVXR's; see `.mip_idx` in expressions/leaf.R.

# ---- the reproducer, in all four accepted spellings -----------------

## min sum(x) s.t. x >= 0.3, x <= 5 with ONLY x[1] boolean.
## CVXPY: obj 1.3, x = (1, 0.3). Relaxed (the old CVXR answer): obj 0.6.
.partial_bool_problem <- function(v) {
  Problem(Minimize(sum_entries(v)), list(v >= 0.3, v <= 5))
}

## @cvxpy NONE
test_that("a partial boolean attribute reaches the solver (all four spellings)", {
  require_solver("GLPK_MI")
  spellings <- list(
    flat   = Variable(2, boolean = c(1)),
    coords = Variable(2, boolean = cbind(1, 1)),
    mask   = Variable(2, boolean = c(TRUE, FALSE))
  )
  for (nm in names(spellings)) {
    v <- spellings[[nm]]
    p <- .partial_bool_problem(v)
    expect_true(is_mixed_integer(p), info = nm)
    expect_equal(psolve(p, solver = "GLPK_MI"), 1.3, tolerance = 1e-6, info = nm)
    expect_equal(as.numeric(value(v)), c(1, 0.3), tolerance = 1e-6, info = nm)
  }

  ## The whole-variable form is unchanged, and the no-attribute form still
  ## solves as a continuous problem -- the two ends the partial case sits between.
  vt <- Variable(2, boolean = TRUE)
  expect_equal(psolve(.partial_bool_problem(vt), solver = "GLPK_MI"), 2, tolerance = 1e-6)
  vn <- Variable(2)
  pn <- .partial_bool_problem(vn)
  expect_false(is_mixed_integer(pn))
  expect_equal(psolve(pn, solver = "GLPK_MI"), 0.6, tolerance = 1e-6)
})

## The upstream test this ports asserts a 2x2 carrying BOTH attributes, each on
## a different pair of entries. CVXPY writes them as numpy multi-indices --
## `boolean=[(0,0),(0,1)]` is row indices (0,0) and column indices (0,1), i.e.
## entries (0,0) and (0,1) -- so the 1-based CVXR coordinates are as below.
## Two independent attributes on one variable is exactly what the old
## `else if` in `.extract_mip_idx` dropped.
## @cvxpy test_attributes.py::TestMultipleAttributes::test_bool_int_variable
test_that("boolean and integer index lists on the same variable both apply", {
  require_solver("GLPK_MI")
  X <- Variable(c(2, 2),
                boolean = cbind(c(1, 1), c(1, 2)),   # entries (1,1), (1,2)
                integer = cbind(c(2, 2), c(1, 2)))   # entries (2,1), (2,2)
  lb <- matrix(c(0.5, 2.3, -0.5, 3.2), 2, 2)
  p <- Problem(Minimize(sum_entries(X)), list(X >= lb))
  psolve(p, solver = "GLPK_MI")
  expect_equal(matrix(as.numeric(value(X)), 2, 2),
               matrix(c(1, 3, 0, 4), 2, 2), tolerance = 1e-6)
})

# ---- the disjointness invariant -------------------------------------

## An entry named by BOTH attributes is not a conflict: {0,1} is the integers
## intersected with [0,1], so applying both yields the boolean constraint.
## `.extract_mip_idx` enforces that once, by removing boolean indices from the
## integer set, rather than leaving it to the statement order inside eight
## solver interfaces. Measured on CVXPY 1.9.2: boolean wins there too.
## @cvxpy NONE
test_that("an entry in both boolean and integer is boolean", {
  require_solver("GLPK_MI")
  y <- Variable(2, boolean = c(1), integer = c(1))
  p <- Problem(Maximize(y[1]), list(y <= 5, y >= 0))
  expect_equal(psolve(p, solver = "GLPK_MI"), 1, tolerance = 1e-6)  # 5 if integer won
})

# ---- offsets across variables ---------------------------------------

## Each variable's indices are shifted by its offset in the stuffed vector
## (matrix_stuffing.py:113). Two variables are the minimum needed to catch a
## missing `+ offset`, which would otherwise half-work.
## @cvxpy NONE
test_that("index lists are offset by the variable's position in the stuffed vector", {
  require_solver("GLPK_MI")
  a <- Variable(3, boolean = c(2))
  b <- Variable(2, integer = c(1))
  p <- Problem(Minimize(sum_entries(a) + sum_entries(b)),
               list(a >= 0.4, b >= 1.2, a <= 9, b <= 9))
  d <- problem_data(p, solver = "GLPK_MI")
  expect_equal(d[[1L]]$bool_idx, 2L)   # a's 2nd entry
  expect_equal(d[[1L]]$int_idx, 4L)    # b's 1st entry, shifted past a's 3
})

# ---- validation ------------------------------------------------------

## R recycles or silently mis-indexes where numpy raises, and silent
## mis-indexing is the failure mode this whole bug class keeps producing, so
## the index forms are validated at construction.
## @cvxpy NONE
test_that("malformed index lists are rejected at construction", {
  expect_error(Variable(2, boolean = c(0)),            "outside 1:2")
  expect_error(Variable(2, boolean = c(99)),           "outside 1:2")
  expect_error(Variable(2, boolean = c(1.5)),          "whole numbers")
  expect_error(Variable(2, boolean = c(1, 1)),         "unique")
  expect_error(Variable(2, boolean = NA_integer_),     "must not contain")
  expect_error(Variable(2, boolean = "x"),             "must be")
  expect_error(Variable(c(2, 2), boolean = cbind(3, 1)), "outside the variable")
  expect_error(Variable(c(2, 2), boolean = c(TRUE, FALSE)), "logical mask")
  ## integer takes the same path
  expect_error(Variable(2, integer = c(7)),            "outside 1:2")
})

# ---- shape-reducing attributes are refused ---------------------------

## symmetric / PSD / NSD / diag replace the variable with a SMALLER one, so an
## index list into the original denotes different entries afterwards: on a 3x3
## symmetric, `boolean = c(5)` means entry (2,2) but position 5 of the
## 6-element upper-triangle vector is entries (3,2) & (2,3). CVXPY refuses this
## combination too (numpy raises inside its index ravel), so refusing is the
## parity behavior -- with a message that says why.
## @cvxpy NONE
test_that("a partial index list cannot be combined with a shape-reducing attribute", {
  expect_error(Variable(c(3, 3), symmetric = TRUE, boolean = c(5)),
               "cannot be combined with")
  expect_error(Variable(c(3, 3), symmetric = TRUE, boolean = c(5, 9)),
               "cannot be combined with")
  expect_error(Variable(c(3, 3), diag = TRUE, integer = c(1)),
               "cannot be combined with")
  expect_error(Variable(c(2, 2), PSD = TRUE, boolean = c(1)),
               "cannot be combined with")
})

## The WHOLE-variable form still works with those attributes: it survives the
## rebuild because it does not name positions.
## @cvxpy NONE
test_that("a whole-variable boolean attribute survives a shape-reducing rebuild", {
  require_solver("GLPK_MI")
  X <- Variable(c(2, 2), symmetric = TRUE, boolean = TRUE)
  p <- Problem(Minimize(sum_entries(X)), list(X >= 0.3, X <= 5))
  expect_equal(psolve(p, solver = "GLPK_MI"), 4, tolerance = 1e-6)
  expect_equal(as.numeric(value(X)), rep(1, 4), tolerance = 1e-6)
})

## A same-shape rebuild (bounds lowering) keeps the index list, since the
## positions still mean what they meant. This is the branch where a length-1
## list used to be silently promoted to TRUE and a longer one raised
## "'length = 2' in coercion to 'logical(1)'".
## @cvxpy NONE
test_that("a same-shape rebuild preserves the index list", {
  require_solver("GLPK_MI")
  v <- Variable(2, nonneg = TRUE, boolean = c(1))
  p <- .partial_bool_problem(v)
  expect_equal(psolve(p, solver = "GLPK_MI"), 1.3, tolerance = 1e-6)
  expect_equal(as.numeric(value(v)), c(1, 0.3), tolerance = 1e-6)

  v2 <- Variable(3, nonneg = TRUE, boolean = c(1, 3))
  p2 <- Problem(Minimize(sum_entries(v2)), list(v2 >= 0.3, v2 <= 5))
  expect_equal(psolve(p2, solver = "GLPK_MI"), 2.3, tolerance = 1e-6)
  expect_equal(as.numeric(value(v2)), c(1, 0.3, 1), tolerance = 1e-6)
})

# ---- Parameters take the same attribute ------------------------------

## A Parameter's boolean/integer attribute constrains its VALUE rather than a
## solver variable, and the same index forms apply. CVXPY rejects a
## non-conforming value with "Parameter value must be boolean." -- measured on
## 1.9.2 -- which is what the message here matches; before, a partial index list
## fell through to "value must be real", true of 0.7 and unrelated to the reason.
## @cvxpy NONE
test_that("a partial boolean Parameter validates its value on the named entries", {
  p <- Parameter(2, boolean = c(1))
  expect_error(value(p) <- c(0.7, 0.7), "must be boolean")
  value(p) <- c(1, 0.7)                       # entry 1 conforms; entry 2 is free
  expect_equal(as.numeric(value(p)), c(1, 0.7))

  q <- Parameter(2, integer = c(2))
  expect_error(value(q) <- c(0.5, 0.5), "must be integer")
  value(q) <- c(0.5, 3)
  expect_equal(as.numeric(value(q)), c(0.5, 3))

  ## Same construction-time refusal as Variables for a shape-reducing attribute.
  expect_error(Parameter(c(2, 2), symmetric = TRUE, boolean = c(1)),
               "cannot be combined with")
})

# ---- integer-only, and the 2-D spellings ------------------------------

## @cvxpy NONE
test_that("a partial integer attribute reaches the solver", {
  require_solver("GLPK_MI")
  v <- Variable(2, integer = c(1))
  p <- Problem(Minimize(sum_entries(v)), list(v >= 0.3, v <= 5))
  expect_true(is_mixed_integer(p))
  expect_equal(psolve(p, solver = "GLPK_MI"), 1.3, tolerance = 1e-6)
  expect_equal(as.numeric(value(v)), c(1, 0.3), tolerance = 1e-6)
})

## Coordinates and masks on a genuinely 2-D variable, where flat position and
## (row, column) differ: entry (2,1) of a 2x3 is flat position 2, entry (1,2) is
## flat position 3. All three spellings must select the same single entry.
## @cvxpy NONE
test_that("coordinates and masks agree with flat indices on a 2-D variable", {
  require_solver("GLPK_MI")
  target_flat <- 4L                      # entry (2,2) of a 2x3, column-major
  mask <- matrix(FALSE, 2, 3); mask[2, 2] <- TRUE
  spellings <- list(
    flat   = Variable(c(2, 3), boolean = target_flat),
    coords = Variable(c(2, 3), boolean = cbind(2, 2)),
    mask   = Variable(c(2, 3), boolean = mask)
  )
  for (nm in names(spellings)) {
    v <- spellings[[nm]]
    p <- Problem(Minimize(sum_entries(v)), list(v >= 0.3, v <= 5))
    expect_equal(psolve(p, solver = "GLPK_MI"), 0.3 * 5 + 1, tolerance = 1e-6, info = nm)
    got <- matrix(as.numeric(value(v)), 2, 3)
    expect_equal(got[2, 2], 1, tolerance = 1e-6, info = nm)
    expect_equal(sum(abs(got - 0.3) < 1e-6), 5L, info = nm)   # the other five untouched
  }
})
