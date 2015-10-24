a <- Variable(name = "a")
b <- Variable(name = "b")
c <- Variable(name = "c")

x <- Variable(2, name = "x")
y <- Variable(3, name = "y")
z <- Variable(2, name = "z")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("test the variables method", {
  p <- Problem(Minimize(a), list(a <= x, b <= A + 2))
  vars_ <- variables(p)
  ref <- list(a, x, b, A)
  mapply(function(v, r) { expect_equal(v, r) }, vars_, ref)
})

test_that("test the parameters method", {
  p1 <- Parameter()
  p2 <- Parameter(3, sign = "negative")
  p3 <- Parameter(4, 4, sign = "positive")
  p <- Problem(Minimize(p1), list(a + p1 <= p2, b <= p3 + p3 + 2))
  params <- parameters(p)
  ref <- c(p1, p2, p3)
  mapply(function(p, r) { expect_equal(p, r) }, params, ref)
})

test_that("test the is_dcp method", {
  p <- Problem(Minimize(NormInf(a)))
  expect_true(is_dcp(p))
  
  p <- Problem(Maximize(NormInf(a)))
  expect_false(is_dcp(p))
})
