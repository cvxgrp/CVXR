## Tests for DGP log-log curvature annotations and is_dgp methods

## @cvxpy NONE
test_that("log-log curvature: leaf variables and constants", {
 x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  z <- Variable(3, pos = TRUE)

  ## Positive variables are log-log affine (both convex and concave)
  expect_true(is_log_log_convex(x))
  expect_true(is_log_log_concave(x))
  expect_true(is_log_log_affine(x))

  ## Non-positive variable is NOT log-log convex/concave
  w <- Variable()
  expect_false(is_log_log_convex(w))
  expect_false(is_log_log_concave(w))

  ## Positive constants are log-log affine
  c1 <- Constant(5)
  expect_true(is_log_log_convex(c1))
  expect_true(is_log_log_concave(c1))
  expect_true(is_log_log_affine(c1))

  ## Non-positive constant is NOT log-log
  c2 <- Constant(-1)
  expect_false(is_log_log_convex(c2))
})

## @cvxpy NONE
test_that("log-log curvature: AddExpression", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  ## Addition of positive vars: log-log convex (posynomial)
  s <- x + y
  expect_true(is_atom_log_log_convex(s))
  expect_false(is_atom_log_log_concave(s))
  expect_true(is_log_log_convex(s))
  expect_false(is_log_log_concave(s))
  expect_false(is_log_log_affine(s))
})

## @cvxpy NONE
test_that("log-log curvature: MulExpression", {
  x <- Variable(pos = TRUE)
  A <- Constant(matrix(1:4, 2, 2))

  ## Matrix multiply is log-log convex only
  m <- A %*% x
  expect_true(is_atom_log_log_convex(m))
  expect_false(is_atom_log_log_concave(m))
})

## @cvxpy NONE
test_that("log-log curvature: Multiply (elementwise)", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  ## Elementwise multiply is log-log affine
  m <- x * y
  expect_true(is_atom_log_log_convex(m))
  expect_true(is_atom_log_log_concave(m))
  expect_true(is_log_log_affine(m))
})

## @cvxpy NONE
test_that("log-log curvature: DivExpression", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  ## Division is log-log affine
  d <- x / y
  expect_true(is_atom_log_log_convex(d))
  expect_true(is_atom_log_log_concave(d))
  expect_true(is_log_log_affine(d))
})

## @cvxpy NONE
test_that("log-log curvature: Power", {
  x <- Variable(pos = TRUE)

  ## x^2 is log-log affine
  p2 <- power(x, 2)
  expect_true(is_atom_log_log_convex(p2))
  expect_true(is_atom_log_log_concave(p2))
  expect_true(is_log_log_affine(p2))

  ## x^(-1) is log-log affine
  pm1 <- power(x, -1)
  expect_true(is_atom_log_log_convex(pm1))
  expect_true(is_atom_log_log_concave(pm1))

  ## x^0.5 is log-log affine
  phalf <- power(x, 0.5)
  expect_true(is_atom_log_log_convex(phalf))
  expect_true(is_atom_log_log_concave(phalf))
})

## @cvxpy NONE
test_that("log-log curvature: Pnorm", {
  x <- Variable(3, pos = TRUE)

  ## p-norm is log-log convex
  pn <- p_norm(x, 2)
  expect_true(is_atom_log_log_convex(pn))
  expect_false(is_atom_log_log_concave(pn))
  expect_true(is_log_log_convex(pn))
  expect_false(is_log_log_affine(pn))
})

## @cvxpy NONE
test_that("log-log curvature: compositions", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  ## Monomial: 2 * x^2 * y^(-1) is log-log affine
  mono <- 2 * power(x, 2) * power(y, -1)
  expect_true(is_log_log_affine(mono))

  ## Posynomial: x + y is log-log convex but not concave
  posy <- x + y
  expect_true(is_log_log_convex(posy))
  expect_false(is_log_log_concave(posy))

  ## Ratio of monomials: x / y is log-log affine
  ratio <- x / y
  expect_true(is_log_log_affine(ratio))

  ## Sum of monomials (posynomial): x^2 + x*y + y^2
  posy2 <- power(x, 2) + x * y + power(y, 2)
  expect_true(is_log_log_convex(posy2))
  expect_false(is_log_log_concave(posy2))
})

## @cvxpy NONE
test_that("log-log curvature: non-DGP expressions", {
  x <- Variable()  ## not positive
  y <- Variable(pos = TRUE)

  ## Non-positive variable makes expressions non-DGP
  expect_false(is_log_log_convex(x + y))
  expect_false(is_log_log_convex(x * y))
})

## @cvxpy NONE
test_that("is_dgp: Minimize objective", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  ## Log-log convex objective: OK for Minimize
  expect_true(is_dgp(Minimize(x + y)))
  expect_true(is_dgp(Minimize(x * y)))
  expect_true(is_dgp(Minimize(power(x, 2))))

  ## Non-DGP objective
  w <- Variable()
  expect_false(is_dgp(Minimize(w)))
})

## @cvxpy NONE
test_that("is_dgp: Maximize objective", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  ## Log-log concave (affine) objective: OK for Maximize
  expect_true(is_dgp(Maximize(x * y)))
  expect_true(is_dgp(Maximize(x / y)))
  expect_true(is_dgp(Maximize(power(x, -1))))

  ## Posynomial (log-log convex only) is NOT valid for Maximize
  expect_false(is_dgp(Maximize(x + y)))
})

## @cvxpy NONE
test_that("is_dgp: Inequality constraint", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  ## lhs log-log convex, rhs log-log concave (monomial): valid DGP
  ## x + y <= x*y  (posynomial <= monomial)
  con1 <- x + y <= x * y
  expect_true(is_dgp(con1))

  ## Scalar constant as RHS: valid (constant is log-log affine)
  con2 <- x + y <= 1
  expect_true(is_dgp(con2))

  ## Monomial <= monomial: valid
  con3 <- x <= 2 * y
  expect_true(is_dgp(con3))

  ## Non-DGP: non-positive variable
  w <- Variable()
  con_bad <- w <= 1
  expect_false(is_dgp(con_bad))
})

## @cvxpy NONE
test_that("is_dgp: Equality constraint", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  ## Both sides log-log affine: valid
  con1 <- x * y == 1
  expect_true(is_dgp(con1))

  ## Monomial == monomial: valid
  con2 <- x == 2 * y
  expect_true(is_dgp(con2))

  ## Posynomial == constant: NOT valid (posynomial is not log-log affine)
  con3 <- (x + y) == 1
  expect_false(is_dgp(con3))
})

## @cvxpy NONE
test_that("is_dgp: Zero constraint", {
  x <- Variable(pos = TRUE)
  ## Zero constraints are never DGP
  z <- Zero(x - 1)
  expect_false(is_dgp(z))
})

## @cvxpy NONE
test_that("is_dgp: Problem", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)

  ## Valid DGP problem: posynomial objective + monomial constraints
  prob1 <- Problem(Minimize(x + y), list(x * y >= 1))
  expect_true(is_dgp(prob1))

  ## Valid DGP: monomial objective with posynomial constraint
  prob2 <- Problem(Minimize(1 / (x * y)), list(x + y <= 1))
  expect_true(is_dgp(prob2))

  ## Not DGP: non-positive variable
  w <- Variable()
  prob3 <- Problem(Minimize(w), list(w >= 0))
  expect_false(is_dgp(prob3))

  ## Not DGP: constraint violates DGP
  prob4 <- Problem(Minimize(x), list(x + y >= 1))
  ## x + y is log-log convex (not concave), so x + y >= 1 fails
  ## Because >= is rewritten as 1 <= x+y, lhs=1 (ll-affine OK),
  ## rhs=x+y (ll-convex but NOT ll-concave, FAIL for Inequality rhs)
  expect_false(is_dgp(prob4))

  ## Not DGP: equality with posynomial
  prob5 <- Problem(Minimize(x), list(x + y == 1))
  expect_false(is_dgp(prob5))
})

## @cvxpy NONE
test_that("is_dgp: cone constraints are not DGP", {
  x <- Variable(3, pos = TRUE)
  ## Cone constraints are never DGP
  expect_false(is_dgp(SOC(x[1], x[2:3])))
})

# ── is_dgp on expressions ──────────────────────────────────────────

## @cvxpy NONE
test_that("is_dgp works on expressions", {
  x <- Variable(pos = TRUE)
  y <- Variable(pos = TRUE)
  expect_true(is_dgp(x))
  expect_true(is_dgp(x * y))
  expect_true(is_dgp(x + y))
  expect_false(is_dgp(x - y))   # subtraction is not log-log convex/concave
})

## @cvxpy NONE
test_that("is_log_log_* exported and works on expressions", {
  x <- Variable(pos = TRUE)
  expect_true(is_log_log_affine(x))
  expect_true(is_log_log_convex(x))
  expect_true(is_log_log_concave(x))
})
