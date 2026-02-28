## Phase 0: Operator dispatch spike
## Validates technical risks R1 (%*% dispatch) and R2 ([ indexing)
## Uses S7's method() dispatch — the mechanism we'll use in the real package.

## ── Minimal S7 Expression spike for operator testing ────────────────

# A minimal expression class for testing operators
SpikeExpr <- S7::new_class("SpikeExpr",
  package = "CVXR",
  properties = list(
    data = S7::class_any,
    shape = S7::class_integer
  ),
  constructor = function(data, shape = NULL) {
    if (is.null(shape)) {
      if (is.matrix(data)) {
        shape <- dim(data)
      } else if (is.numeric(data)) {
        shape <- c(length(data), 1L)
      } else {
        shape <- c(1L, 1L)
      }
    }
    S7::new_object(S7::S7_object(),
                   data = data,
                   shape = as.integer(shape))
  }
)

# ── Risk R1: Test %*% (matmul) dispatch ──────────────────────────────

# Register %*% methods via S7::method()
S7::method(`%*%`, list(SpikeExpr, S7::class_any)) <- function(x, y) {
  list(op = "matmul", direction = "lhs")
}

S7::method(`%*%`, list(S7::class_any, SpikeExpr)) <- function(x, y) {
  list(op = "matmul", direction = "rhs")
}

## @cvxpy NONE
test_that("Risk R1: %*% works — SpikeExpr %*% matrix (LHS)", {
  e <- SpikeExpr(matrix(1:6, 2, 3))
  m <- matrix(1:12, 3, 4)

  result <- e %*% m
  expect_equal(result$op, "matmul")
  expect_equal(result$direction, "lhs")
})

## @cvxpy NONE
test_that("Risk R1: %*% works — matrix %*% SpikeExpr (RHS, critical)", {
  e <- SpikeExpr(matrix(1:6, 2, 3))
  m <- matrix(1:6, 2, 3)

  # This is the critical test: plain matrix on LHS, our class on RHS
  result <- m %*% e
  expect_equal(result$op, "matmul")
  expect_equal(result$direction, "rhs")
})

## @cvxpy NONE
test_that("Risk R1: %*% works — SpikeExpr %*% SpikeExpr", {
  e1 <- SpikeExpr(matrix(1:6, 2, 3))
  e2 <- SpikeExpr(matrix(1:12, 3, 4))

  result <- e1 %*% e2
  expect_equal(result$op, "matmul")
  expect_equal(result$direction, "lhs")
})

# ── Risk R2: Test [ indexing dispatch ────────────────────────────────

S7::method(`[`, SpikeExpr) <- function(x, i, j, ..., drop = TRUE) {
  has_i <- !missing(i)
  has_j <- !missing(j)
  list(op = "index",
       has_i = has_i, has_j = has_j,
       i = if (has_i) i else NULL,
       j = if (has_j) j else NULL,
       drop = drop)
}

## @cvxpy NONE
test_that("Risk R2: [ indexing works — e[i, j]", {
  e <- SpikeExpr(matrix(1:12, 3, 4))

  result <- e[2, 3]
  expect_equal(result$op, "index")
  expect_equal(result$i, 2)
  expect_equal(result$j, 3)
  expect_true(result$has_i)
  expect_true(result$has_j)
})

## @cvxpy NONE
test_that("Risk R2: [ indexing works — e[i, ] (row selection, missing j)", {
  e <- SpikeExpr(matrix(1:12, 3, 4))

  result <- e[1:2, ]
  expect_equal(result$op, "index")
  expect_equal(result$i, 1:2)
  expect_true(result$has_i)
  expect_false(result$has_j)
  expect_null(result$j)
})

## @cvxpy NONE
test_that("Risk R2: [ indexing works — e[, j] (column selection, missing i)", {
  e <- SpikeExpr(matrix(1:12, 3, 4))

  result <- e[, 2:3]
  expect_equal(result$op, "index")
  expect_false(result$has_i)
  expect_null(result$i)
  expect_true(result$has_j)
  expect_equal(result$j, 2:3)
})

## @cvxpy NONE
test_that("Risk R2: [ indexing works — e[i] (single index)", {
  e <- SpikeExpr(1:10)

  result <- e[3]
  expect_equal(result$op, "index")
  expect_equal(result$i, 3)
})

# ── Test Ops dispatch (+, -, *, /, ==, <=, >=) ──────────────────────

S7::method(`+`, list(SpikeExpr, SpikeExpr)) <- function(e1, e2) {
  list(op = "+", lhs_class = "SpikeExpr", rhs_class = "SpikeExpr")
}

S7::method(`+`, list(SpikeExpr, S7::class_numeric)) <- function(e1, e2) {
  list(op = "+", lhs_class = "SpikeExpr", rhs_class = "numeric")
}

S7::method(`+`, list(S7::class_numeric, SpikeExpr)) <- function(e1, e2) {
  list(op = "+", lhs_class = "numeric", rhs_class = "SpikeExpr")
}

S7::method(`-`, list(SpikeExpr, SpikeExpr)) <- function(e1, e2) {
  list(op = "-", lhs_class = "SpikeExpr", rhs_class = "SpikeExpr")
}

S7::method(`*`, list(SpikeExpr, S7::class_numeric)) <- function(e1, e2) {
  list(op = "*", lhs_class = "SpikeExpr", rhs_class = "numeric")
}

S7::method(`*`, list(S7::class_numeric, SpikeExpr)) <- function(e1, e2) {
  list(op = "*", lhs_class = "numeric", rhs_class = "SpikeExpr")
}

S7::method(`==`, list(SpikeExpr, SpikeExpr)) <- function(e1, e2) {
  list(op = "==", type = "equality_constraint")
}

S7::method(`<=`, list(SpikeExpr, SpikeExpr)) <- function(e1, e2) {
  list(op = "<=", type = "inequality_constraint")
}

S7::method(`>=`, list(SpikeExpr, SpikeExpr)) <- function(e1, e2) {
  list(op = ">=", type = "inequality_constraint")
}

# ── Unary operator workaround ─────────────────────────────────────────
# S7's Ops.S7_object always evaluates both e1 and e2, which fails for
# unary operators (- and +). Workaround: register an S3 Ops method for
# the specific S7 class that intercepts unary ops before S7's handler,
# then delegates binary ops via NextMethod().
.spike_ops_handler <- function(e1, e2) {
  if (missing(e2)) {
    if (.Generic == "-") return(list(op = "unary_neg"))
    if (.Generic == "+") return(e1)
    stop("Unsupported unary operator: ", .Generic)
  }
  NextMethod()  # Delegates to S7::Ops.S7_object for binary ops
}
registerS3method("Ops", "CVXR::SpikeExpr", .spike_ops_handler)

## @cvxpy NONE
test_that("Ops: SpikeExpr + SpikeExpr", {
  e1 <- SpikeExpr(1)
  e2 <- SpikeExpr(2)
  result <- e1 + e2
  expect_equal(result$op, "+")
  expect_equal(result$lhs_class, "SpikeExpr")
  expect_equal(result$rhs_class, "SpikeExpr")
})

## @cvxpy NONE
test_that("Ops: SpikeExpr + numeric", {
  e <- SpikeExpr(1)
  result <- e + 5
  expect_equal(result$op, "+")
  expect_equal(result$rhs_class, "numeric")
})

## @cvxpy NONE
test_that("Ops: numeric + SpikeExpr (reverse dispatch)", {
  e <- SpikeExpr(1)
  result <- 5 + e
  expect_equal(result$op, "+")
  expect_equal(result$lhs_class, "numeric")
  expect_equal(result$rhs_class, "SpikeExpr")
})

## @cvxpy NONE
test_that("Ops: SpikeExpr * numeric and numeric * SpikeExpr", {
  e <- SpikeExpr(1)
  r1 <- e * 3
  expect_equal(r1$op, "*")
  r2 <- 3 * e
  expect_equal(r2$op, "*")
})

## @cvxpy NONE
test_that("Ops: comparison operators return constraint-like objects", {
  e1 <- SpikeExpr(1)
  e2 <- SpikeExpr(2)

  r_eq <- e1 == e2
  expect_equal(r_eq$op, "==")
  expect_equal(r_eq$type, "equality_constraint")

  r_leq <- e1 <= e2
  expect_equal(r_leq$op, "<=")

  r_geq <- e1 >= e2
  expect_equal(r_geq$op, ">=")
})

## @cvxpy NONE
test_that("Ops: unary negation works", {
  e <- SpikeExpr(1)
  result <- -e
  expect_equal(result$op, "unary_neg")
})
