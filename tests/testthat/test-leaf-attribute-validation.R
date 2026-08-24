## @cvxpy tests/test_expressions.py, tests/test_problem.py
##
## Leaf attribute validation and lowering: the checks CVXPY's `Leaf.__init__`
## and `_ensure_valid_bounds` perform, plus the two callers of `_bound_domain`.
##
## Why this file exists: each of these was a SILENT wrong answer, not an error.
## An attribute CVXR accepted but did not enforce does not announce itself --
## the problem solves, the status is "optimal", and the number is wrong.
##
## CVXPY SOURCE: leaf.py:155-160 (dimension-reducing cardinality),
## leaf.py:839-930 (_ensure_valid_bounds), leaf.py:350-414 (_bound_domain),
## leaf.py:416-430 (domain), leaf.py:465-475 (project), leaf.py:670-671
## (the 'in bounds' arm of the attribute-name cascade).

## ---------------------------------------------------------------------------
## At most ONE dimension-reducing attribute (CVXPY leaf.py:155-160)
##
## Pre-fix CVXR accepted every one of these and SILENTLY DROPPED the second
## attribute, because the reduction's `if (has_sym || has_psd || has_nsd)` /
## `else if (has_diag)` chain assumes an exclusivity nothing enforced.
## Measured pre-fix: symmetric+diag gave -4 for a problem whose diagonal answer
## is -2; PSD+NSD gave 1 for a problem that is infeasible if both hold.
## ---------------------------------------------------------------------------

test_that("more than one dimension-reducing attribute is rejected", {
  expect_error(Variable(c(2, 2), symmetric = TRUE, diag = TRUE), "more than one")
  expect_error(Variable(c(2, 2), PSD = TRUE, NSD = TRUE), "more than one")
  expect_error(Variable(c(2, 2), PSD = TRUE, diag = TRUE), "more than one")
  expect_error(Variable(c(2, 2), symmetric = TRUE, hermitian = TRUE), "more than one")
  expect_error(Variable(c(2, 2), NSD = TRUE, diag = TRUE), "more than one")
  ## Parameters go through the same validator.
  expect_error(Parameter(c(2, 2), symmetric = TRUE, diag = TRUE), "more than one")
})

test_that("exactly one dimension-reducing attribute is still accepted", {
  for (a in c("symmetric", "diag", "PSD", "NSD", "hermitian")) {
    args <- list(c(2, 2)); args[[a]] <- TRUE
    expect_no_error(do.call(Variable, args))
  }
})

## The exclusivity rule is what makes these two answers well defined at all.
test_that("a single dimension-reducing attribute still lowers correctly", {
  X <- Variable(c(2, 2), diag = TRUE)
  expect_equal(as.numeric(psolve(Problem(Minimize(sum_entries(X)), list(X >= -1)),
                                 solver = "CLARABEL")), -2, tolerance = 1e-6)
  Y <- Variable(c(2, 2), symmetric = TRUE)
  expect_equal(as.numeric(psolve(Problem(Minimize(sum_entries(Y)), list(Y >= -1)),
                                 solver = "CLARABEL")), -4, tolerance = 1e-6)
})

## ---------------------------------------------------------------------------
## Numeric array bounds must match the leaf's shape
## (CVXPY leaf.py:909-921 -- `valid_array = ... and val.shape == self.shape`)
##
## Pre-fix, an n x n bounds array survived onto the n(n+1)/2 upper-triangular
## variable that CvxAttr2Constr builds for a symmetric leaf, and
## `.extract_lower_bounds` then reshaped it with `array(b, dim = .shape(v))`,
## TRUNCATING it. Measured pre-fix on the 3x3 below, whose true optimum is
## sum(LB3) = -31:  HIGHS -26, while CLARABEL/SCS/OSQP gave -31. Feasible but
## suboptimal, so nothing downstream noticed. CVXPY raises instead, and only on
## the HiGHS path, because BOUNDED_VARIABLES solvers are the ones that keep the
## bounds attribute (reduce_bounds = FALSE).
## ---------------------------------------------------------------------------

test_that("mis-shaped numeric bounds are rejected at construction", {
  expect_error(Variable(c(3, 3), bounds = list(matrix(-1, 2, 2), matrix(1, 2, 2))),
               "same dimensions")
  expect_error(Variable(c(3, 3), bounds = list(rep(-1, 4), rep(1, 4))),
               "same dimensions")
  ## A dim-carrying bound of the right LENGTH but the wrong shape.
  expect_error(Variable(c(2, 3), bounds = list(matrix(-1, 3, 2), matrix(1, 3, 2))),
               "same dimensions")
})

test_that("the R spellings of a correctly sized bound are still accepted", {
  ## Dim-less vector of the right total length: ordinary R, and CVXR shapes a
  ## vector variable as c(n, 1), so this must not be read as a shape mismatch.
  v <- Variable(3, bounds = list(rep(-1, 3), rep(1, 3)))
  expect_equal(as.numeric(psolve(Problem(Minimize(sum_entries(v))), solver = "CLARABEL")),
               -3, tolerance = 1e-6)
  expect_no_error(Variable(c(2, 3), bounds = list(matrix(-1, 2, 3), matrix(1, 2, 3))))
  expect_no_error(Variable(c(2, 3), bounds = list(-1, 1)))          # scalars
  expect_no_error(Variable(c(2, 3), bounds = list(NULL, NULL)))     # NULLs
})

test_that("per-entry bounds on a symmetric variable agree across solvers", {
  LB3 <- matrix(c(-1, -2, -3, -2, -4, -5, -3, -5, -6), 3, 3)
  UB3 <- matrix(1, 3, 3)
  expect_equal(sum(LB3), -31)
  for (s in intersect(c("CLARABEL", "SCS", "OSQP"), installed_solvers())) {
    X <- Variable(c(3, 3), symmetric = TRUE, bounds = list(LB3, UB3))
    expect_equal(as.numeric(psolve(Problem(Minimize(sum_entries(X))), solver = s)),
                 -31, tolerance = 1e-4, label = paste("objective under", s))
  }
  ## The reduction cannot carry 9 bounds onto a 6-entry variable. Erroring is
  ## what CVXPY does here; the pre-fix behavior was to truncate and return -26.
  skip_if_not("HIGHS" %in% installed_solvers(), "HiGHS not installed")
  X2 <- Variable(c(3, 3), symmetric = TRUE, bounds = list(LB3, UB3))
  expect_error(psolve(Problem(Minimize(sum_entries(X2))), solver = "HIGHS"),
               "same dimensions")
})

## ---------------------------------------------------------------------------
## project() reads the canonical index, not the raw attribute
## (CVXPY leaf.py:465-475, indexing with self.boolean_idx / integer_idx)
##
## `_validate_value` routes every `value<-` through project (leaf.py:608), so a
## projection bug is a validation bug. Pre-fix, of CVXR's four documented
## spellings two broke in OPPOSITE directions: the coordinate-matrix form
## rounded the wrong entries (as.integer(cbind(c(1,2), c(1,2))) is 1,2,1,2) and
## REJECTED a valid value; the logical-mask form matched no branch at all, so
## project() was a no-op and an INVALID value was silently accepted.
## ---------------------------------------------------------------------------

.bool_spellings <- function() list(
  "TRUE"        = TRUE,
  "flat"        = c(1, 4),
  "coords"      = cbind(c(1, 2), c(1, 2)),
  "logical mask" = matrix(c(TRUE, FALSE, FALSE, TRUE), 2, 2)
)

test_that("every boolean spelling yields the same canonical index", {
  sp <- .bool_spellings()
  for (nm in c("flat", "coords", "logical mask")) {
    v <- Variable(c(2, 2), boolean = sp[[nm]])
    expect_equal(as.integer(v@.boolean_idx), c(1L, 4L), label = paste("index for", nm))
  }
  expect_equal(as.integer(Variable(c(2, 2), boolean = TRUE)@.boolean_idx), 1:4)
})

test_that("project() touches exactly the canonical index, for every spelling", {
  sp <- .bool_spellings()
  val <- matrix(0.7, 2, 2)
  for (nm in c("flat", "coords", "logical mask")) {
    v <- Variable(c(2, 2), boolean = sp[[nm]])
    got <- as.numeric(CVXR:::project(v, val))
    ## entries 1 and 4 rounded to 1; entries 2 and 3 untouched
    expect_equal(got, c(1, 0.7, 0.7, 1), label = paste("projection for", nm))
  }
  expect_equal(as.numeric(CVXR:::project(Variable(c(2, 2), boolean = TRUE), val)),
               rep(1, 4))
})

test_that("value<- accepts and rejects identically across boolean spellings", {
  sp <- .bool_spellings()
  ok  <- matrix(c(1, 0.7, 0.7, 0), 2, 2)   # diagonal is 0/1  -> valid
  bad <- matrix(c(0.7, 0, 0, 0.7), 2, 2)   # diagonal is 0.7  -> invalid
  for (nm in c("flat", "coords", "logical mask")) {
    v <- Variable(c(2, 2), boolean = sp[[nm]])
    expect_no_error({ value(v) <- ok })
    w <- Variable(c(2, 2), boolean = sp[[nm]])
    expect_error({ value(w) <- bad }, "boolean", label = paste("rejection for", nm))
  }
})

test_that("integer spellings project through the canonical index too", {
  v <- Variable(c(2, 2), integer = cbind(c(1, 2), c(1, 2)))
  expect_equal(as.integer(v@.integer_idx), c(1L, 4L))
  expect_equal(as.numeric(CVXR:::project(v, matrix(1.4, 2, 2))), c(1, 1.4, 1.4, 1))
})

## ---------------------------------------------------------------------------
## domain(Leaf) -- the second caller of _bound_domain (CVXPY leaf.py:416-430)
## Pre-fix this was `function(x) list()`.
## ---------------------------------------------------------------------------

test_that("domain() of a leaf reports its attribute constraints", {
  expect_length(CVXR:::domain(Variable(2)), 0)
  expect_length(CVXR:::domain(Variable(2, nonneg = TRUE)), 1)
  expect_length(CVXR:::domain(Variable(2, pos = TRUE)), 1)
  expect_length(CVXR:::domain(Variable(2, nonpos = TRUE)), 1)
  expect_length(CVXR:::domain(Variable(2, neg = TRUE)), 1)
  ## lower and upper are separate constraints
  expect_length(CVXR:::domain(Variable(2, bounds = list(-1, 1))), 2)
  ## semidefiniteness, CVXPY leaf.py:426-429
  expect_length(CVXR:::domain(Variable(c(2, 2), PSD = TRUE)), 1)
  expect_length(CVXR:::domain(Variable(c(2, 2), NSD = TRUE)), 1)
  ## Parameters share the implementation; a Constant has no attributes.
  expect_length(CVXR:::domain(Parameter(2, nonneg = TRUE)), 1)
  expect_length(CVXR:::domain(Constant(3)), 0)
})

test_that("the domain constraint is about the right leaf and the right side", {
  x <- Variable(2, nonneg = TRUE)
  d <- CVXR:::domain(x)[[1L]]
  ## It constrains x itself, not a copy or an auxiliary.
  ids <- vapply(variables(d), function(v) as.integer(CVXR:::.id(v)), integer(1))
  expect_true(as.integer(CVXR:::.id(x)) %in% ids)
  ## Satisfied wherever the attribute is satisfied. (An attribute-violating
  ## value cannot be ASSIGNED to test the other side -- `value<-` projects and
  ## rejects it -- which is exactly the invariant domain() is reporting.)
  value(x) <- c(0, 0)
  expect_true(all(as.numeric(violation(d)) <= 1e-9))
  value(x) <- c(3, 4)
  expect_true(all(as.numeric(violation(d)) <= 1e-9))

  ## A bounds leaf gives one constraint per side, and the upper one binds.
  b <- Variable(2, bounds = list(-1, 1))
  db <- CVXR:::domain(b)
  expect_length(db, 2)
  value(b) <- c(1, 1)
  expect_true(all(vapply(db, function(cc) max(as.numeric(violation(cc))), numeric(1)) <= 1e-9))
})

test_that("domain() of a partial_optimize collects its leaves' constraints", {
  x <- Variable(2, nonneg = TRUE)
  y <- Variable(2, nonneg = TRUE)
  po <- partial_optimize(Problem(Minimize(sum_entries(x) + sum_entries(y)),
                                 list(x + y >= 1)), list(y), list(x))
  ## CVXPY 1.9.2 returns 3 here; pre-fix CVXR returned 1, because the two leaf
  ## contributions were lost to the `list()` stub.
  expect_length(CVXR:::domain(po), 3)
})

## ---------------------------------------------------------------------------
## The attribute-name cascade (CVXPY leaf.py:653-676)
## Pre-fix an out-of-bounds assignment said "value must be real."
## ---------------------------------------------------------------------------

test_that("rejecting a value names the attribute that rejected it", {
  v <- Variable(2, bounds = list(-1, 1))
  expect_error({ value(v) <- c(5, 5) }, "in bounds")
  n <- Variable(2, nonneg = TRUE)
  expect_error({ value(n) <- c(-5, -5) }, "nonnegative")
  b <- Variable(2, boolean = TRUE)
  expect_error({ value(b) <- c(0.7, 0.7) }, "boolean")
})
