a <- Variable(name = "a")

x <- Variable(2, name = "x")
y <- Variable(2, name = "y")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("Test gradient for AffineProd", {
  expr <- AffineProd(C, A)
  C@value <- rbind(c(1,-2), c(3,4), c(-1,-3))
  A@value <- rbind(c(3,2), c(-5,1))
  
  expect_equal(as.matrix(grad(expr)[[C@id]]), rbind(c(3,0,0,2,0,0), c(0,3,0,0,2,0), c(0,0,3,0,0,2),
                                                    c(-5,0,0,1,0,0), c(0,-5,0,0,1,0), c(0,0,-5,0,0,1)), tolereance = TOL)
  expect_equal(as.matrix(grad(expr)[[A@id]]), rbind(c(1,3,-1,0,0,0), c(-2,4,-3,0,0,0), c(0,0,0,1,3,-1), 
                                                   c(0,0,0,-2,4,-3)), tolerance = TOL)
})

test_that("Test gradient for Pnorm", {
  expr <- Pnorm(x, 1)
  x@value <- c(-1,0)
  expect_equal(as.matrix(grad(expr)[[x@id]]), as.matrix(c(-1,0)), tolerance = TOL)
  
  x@value <- c(0,10)
  expect_equal(as.matrix(grad(expr)[[x@id]]), as.matrix(c(0,1)), tolerance = TOL)
  
  expr <- Pnorm(x, 2)
  x@value <- c(-3,4)
  expect_equal(as.matrix(grad(expr)[[x@id]]), as.matrix(c(-3.0/5),4.0/5), tolerance = TOL)
  
  expr <- Pnorm(x, 0.5)
  x@value <- c(-1,2)
  expect_equal(grad(expr)[[x@id]], NA)
  
  expr <- Pnorm(x, 0.5)
  x@value <- c(0,0)
  expect_equal(grad(expr)[[x@id]], NA)
  
  expr <- Pnorm(x, 2)
  x@value <- c(0,0)
  expect_equal(as.matrix(grad(expr)[[x@id]]), as.matrix(c(0,0)))
  
  expr <- Pnorm(A, 2)
  A@value <- rbind(c(2,-2), c(2,2))
  expect_equal(as.matrix(grad(expr)[[A@id]]), matrix(c(0.5,0.5,-0.5,0.5)))
  
  expr <- Pnorm(A, 2, axis = 0)
  A@value <- rbind(c(3,-3), c(4,4))
  expect_equal(as.matrix(grad(expr)[[A@id]]), rbind(c(0.6,0), c(0.8,0), c(0,-0.6), c(0,0.8)))
  
  expr <- Pnorm(A, 2, axis = 1)
  A@value <- rbind(c(3,-4), c(4,3))
  expect_equal(as.matrix(grad(expr)[[A@id]]), rbind(c(0.6,0), c(0,0.8), c(-0.8,0), c(0,0.6)))
  
  expr <- Pnorm(A, 0.5)
  A@value <- rbind(c(3,-4), c(4,3))
  expect_equal(grad(expr)[[A@id]], NA)
})

test_that("Test gradient for LogSumExp", {
  expr <- LogSumExp(x)
  x@value <- c(0,1)
  e <- exp(1)
  expect_equal(as.matrix(grad(expr)[[x@id]]), c(1.0/(1+e), e/(1+e)))
  
  expr <- LogSumExp(A)
  A@value <- rbind(c(0,1), c(-1,0))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(1.0/(2+e+1.0/e), 1.0/e/(2+e+1.0/e), e/(2+e+1.0/e), 1.0/(2+e+1.0/e)))
  
  expr <- LogSumExp(A, axis = 0)
  A@value <- rbind(c(0,1), c(-1,0))
  expect_equal(as.matrix(grad(expr)[[A@id]]), cbind(c(1.0/(1+1.0/e), 1.0/e/(1+1.0/e), 0, 0), c(0, 0, e/(1+e), 1.0/(1+e))))
})

test_that("Test gradient for GeoMean", {
  expr <- GeoMean(x)
  x@value <- c(1,2)
  expect_equal(as.matrix(grad(expr)[[x@id]]), c(sqrt(2)/2, 1.0/2/sqrt(2)))
  
  x@value <- c(0,2)
  expect_equal(grad(expr)[[x@id]], NA)
  
  expr <- GeoMean(x, c(1,0))
  x@value <- c(1,2)
  expect_equal(as.matrix(grad(expr)[[x@id]]), c(1,0))
  
  # No exception for single weight
  x@value <- c(-1,2)
  expect_equal(grad(expr)[[x@id]], NA)
})

test_that("Test gradient for LambdaMax", {
  expr <- LambdaMax(A)
  A@value <- rbind(c(2,0), c(0,1))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(1,0,0,0))
  
  A@value <- rbind(c(1,0), c(0,2))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(0,0,0,1))
  
  A@value <- rbind(c(1,0), c(0,1))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(0,0,0,1))
})

test_that("Test gradient for MatrixFrac", {
  expr <- MatrixFrac(A, B)
  A@value <- diag(rep(1, 2))
  B@value <- diag(rep(1, 2))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(2,0,0,2))
  expect_equal(as.matrix(grad(expr)[[B@id]]), c(-1,0,0,-1))
  
  B@value <- matrix(0, nrow = 2, ncol = 2)
  expect_equal(grad(expr)[[A@id]], NA)
  expect_equal(grad(expr)[[B@id]], NA)
  
  expr <- MatrixFrac(x, A)
  x@value <- c(2,3)
  A@value <- diag(rep(1,2))
  expect_equal(as.matrix(grad(expr)[[x@id]]), c(4,6))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(-4,-6,-6,-9))
})

test_that("Test gradient for NormNuc", {
  expr <- NormNuc(A)
  A@value <- rbind(c(10,4), c(4,30))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(1,0,0,1))
})

test_that("Test gradient for LogDet", {
  expr <- LogDet(A)
  A@value <- 2*diag(rep(1,2))
  expect_equal(as.matrix(grad(expr)[[A@id]]), 1.0/2*diag(rep(1,2)))
  
  mat <- rbind(c(1,2), c(3,5))
  A@value <- t(mat) %*% mat
  val <- t(solve(A@value))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
  
  A@value <- matrix(0, nrow = 2, ncol = 2)
  expect_equal(grad(expr)[[A@id]], NA)
  
  A@value <- -rbind(c(1,2), c(3,4))
  expect_equal(grad(expr)[[A@id]], NA)
})

test_that("Test gradient for QuadOverLin", {
  expr <- QuadOverLin(x, a)
  x@value <- c(1,2)
  a@value <- 2
  expect_equal(as.matrix(grad(expr)[[x@id]]), c(1,2))
  expect_equal(grad(expr)[[a@id]], -1.25)
  
  a@value <- 0
  expect_equal(grad(expr)[[x@id]], NA)
  expect_equal(grad(expr)[[a@id]], NA)
  
  expr <- QuadOverLin(A, a)
  A@value <- diag(rep(1,2))
  a@value <- 2
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(1,0,0,1))
  expect_equal(grad(expr)[[a@id]], -0.5)
  
  expr <- QuadOverLin(x, a) + QuadOverLin(y, a)
  x@value <- c(1,2)
  a@value <- 2
  y@value <- c(1,2)
  a@value <- 2
  expect_equal(as.matrix(grad(expr)[[x@id]]), c(1,2))
  expect_equal(as.matrix(grad(expr)[[y@id]]), c(1,2))
  expect_equal(grad(expr)[[a@id]], -2.5)
})

test_that("Test gradient for MaxEntries", {
  expr <- MaxEntries(x)
  x@value <- c(2,1)
  expect_equal(as.matrix(grad(expr)[[x@id]]), c(1,0))
  
  expr <- MaxEntries(A)
  A@value <- rbind(c(1,2), c(4,3))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(0,1,0,0))
  
  expr <- MaxEntries(A, axis = 0)
  A@value <- rbind(c(1,2), c(4,3))
  expect_equal(as.matrix(grad(expr)[[A@id]]), rbind(c(0,0), c(1,0), c(0,0), c(0,1)))
  
  expr <- MaxEntries(A, axis = 1)
  A@value <- rbind(c(1,2), c(4,3))
  expect_equal(as.matrix(grad(expr)[[A@id]]), rbind(c(0,0), c(0,1), c(1,0), c(0,0)))
})

test_that("Test SigmaMax", {
  expr <- SigmaMax(A)
  A@value <- rbind(c(1,0), c(0,2))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(0,0,0,1))
  
  A@value <- rbind(c(1,0), c(0,1))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(1,0,0,0))
})

test_that("Test SumLargest", {
  expr <- SumLargest(A, 2)
  
  A@value <- rbind(c(4,3), c(2,1))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(1,0,1,0))
  
  A@value <- rbind(c(1,2), c(3,0.5))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(0,1,1,0))
})

test_that("Test Abs", {
  expr <- Abs(A)
  A@value <- rbind(c(1,2), c(-1,0))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(1,1,-1,0))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
})

test_that("Test linearize method", {
  # Affine
  expr <- (2*x - 5)[1]
  x@value <- c(1,2)
  lin_expr <- linearize(expr)
  x@value <- c(55,22)
  expect_equal(value(lin_expr), value(expr))
  x@value <- c(-1,-5)
  expect_equal(value(lin_expr), value(expr))
  
  # Convex
  expr <- A^2 + 5
  expect_error(linearize(expr))
  
  A@value <- rbind(c(1,2), c(3,4))
  lin_expr <- linearize(expr)
  manual <- value(expr) + 2*Reshape(value(Diag(Vec(A))) * Vec(A - A@value), 2, 2)
  expect_equal(value(lin_expr), value(expr))
  A@value <- rbind(c(-5,-5), c(8.2,4.4))
  expect_true(all(value(lin_expr) <= value(expr)))
  expect_equal(value(lin_expr), value(manual))
  
  # Concave
  expr <- Log(x)/2
  x@value <- c(1,2)
  lin_expr <- linearize(expr)
  manual <- value(expr) + value(Diag(0.5*x^-1))*(x - x@value)
  expect_equal(value(lin_expr), value(expr))
  x@value <- c(3,4.4)
  expect_true(all(value(lin_expr) >= value(expr)))
  expect_equal(value(lin_expr), value(manual))
})

test_that("Test gradient for Log", {
  expr <- Log(a)
  a@value <- 2
  expect_equal(grad(expr)[[a@id]], 1.0/2)
  
  a@value <- 3
  expect_equal(grad(expr)[[a@id]], 1.0/3)
  
  a@value <- -1
  expect_equal(grad(expr)[[a@id]], NA)
  
  expr <- Log(x)
  x@value <- c(3,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1/3,1/4))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Log(x)
  x@value <- c(1e-9,4)
  expect_equal(grad(expr)[[x@id]], NA)
  
  expr <- Log(A)
  A@value <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(1, 1/2, 1/3, 1/4))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
})

test_that("Test domain for Log1p", {
  expr <- Log1p(a)
  a@value <- 2
  expect_equal(grad(expr)[[a@id]], 1.0/3)
  
  a@value <- 3
  expect_equal(grad(expr)[[a@id]], 1.0/4)
  
  a@value <- -1
  expect_equal(grad(expr)[[a@id]], NA)
  
  expr <- Log1p(x)
  x@value <- c(3,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1/4,1/5))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Log1p(x)
  x@value <- c(-1e-9-1,4)
  expect_equal(grad(expr)[[x@id]], NA)
  
  expr <- Log1p(A)
  A@value <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(1/2, 1/3, 1/4, 1/5))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
})

test_that("Test domain for Entr", {
  expr <- Entr(a)
  a@value <- 2
  expect_equal(grad(expr)[[a@id]], -log(2)-1)
  
  a@value <- 3
  expect_equal(grad(expr)[[a@id]], -(log(3)+1))
  
  a@value <- -1
  expect_equal(grad(expr)[[a@id]], NA)
  
  expr <- Entr(x)
  x@value <- c(3,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(-(log(c(3,4)) + 1))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Entr(x)
  x@value <- c(-1e-9,4)
  expect_equal(grad(expr)[[x@id]], NA)
  
  expr <- Entr(A)
  A@value <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(-(log(1:4)+1))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
})

test_that("Test domain for Exp", {
  expr <- Exp(a)
  a@value <- 2
  expect_equal(grad(expr)[[a@id]], exp(2))
  
  a@value <- 3
  expect_equal(grad(expr)[[a@id]], exp(3))
  
  a@value <- -1
  expect_equal(grad(expr)[[a@id]], exp(-1))
  
  expr <- Exp(x)
  x@value <- c(3,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(exp(c(3,4)))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Exp(x)
  x@value <- c(-1e-9,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(exp(c(-1e-9,4)))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Exp(A)
  A@value <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(exp(1:4))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
})

test_that("Test domain for logistic", {
  expr <- Logistic(a)
  a@value <- 2
  expect_equal(grad(expr)[[a@id]], exp(2)/(1+exp(2)))
  
  a@value <- 3
  expect_equal(grad(expr)[[a@id]], exp(3)/(1+exp(3)))
  
  a@value <- -1
  expect_equal(grad(expr)[[a@id]], exp(-1)/(1+exp(-1)))
  
  expr <- Logistic(x)
  x@value <- c(3,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(exp(c(3,4))/(1+exp(3,4)))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Logistic(x)
  x@value <- c(-1e-9,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(exp(c(-1e-9,4))/(1+exp(c(-1e-9,4))))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Logistic(A)
  A@value <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(exp(1:4)/(1+exp(1:4)))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
})

test_that("Test domain for Huber", {
  expr <- Huber(a)
  a@value <- 2
  expect_equal(grad(expr)[[a@id]], 2)
  
  expr <- Huber(a, M = 2)
  a@value <- 3
  expect_equal(grad(expr)[[a@id]], 4)
  
  a@value <- -1
  expect_equal(grad(expr)[[a@id]], -2)
  
  expr <- Huber(x)
  x@value <- c(3,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(2,2))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Huber(x)
  x@value <- c(-1e-9,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(0,2))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Huber(A, M = 3)
  A@value <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(2,4,6,6))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
})

test_that("Test domain for KLDiv", {
  b <- Variable()
  expr <- KLDiv(a, b)
  a@value <- 2
  b@value <- 4
  expect_equal(grad(expr)[[a@id]], log(2/4))
  expect_equal(grad(expr)[[b@id]], 1-(2/4))
  
  a@value <- 3
  b@value <- 0
  expect_equal(grad(expr)[[a@id]], NA)
  expect_equal(grad(expr)[[b@id]], NA)
  
  a@value <- -1
  b@value <- 2
  expect_equal(grad(expr)[[a@id]], NA)
  expect_equal(grad(expr)[[b@id]], NA)
  
  y <- Variable(2)
  expr <- KLDiv(x, y)
  x@value <- c(3,4)
  y@value <- c(5,8)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(log(c(3,4)) - log(c(5,8)))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  val <-matrix(0, nrow = 2, ncol = 2) + diag(c(1-3/5,1-4/8))
  expect_equal(as.matrix(grad(expr)[[y@id]]), val)
  
  expr <- KLDiv(x, y)
  x@value <- c(-1e-9,4)
  y@value <- c(1,2)
  expect_equal(grad(expr)[[x@id]], NA)
  expect_equal(grad(expr)[[y@id]], NA)
  
  expr <- KLDiv(A, B)
  A@value <- rbind(c(1,2), c(3,4))
  B@value <- rbind(c(5,1), c(3.5,2.3))
  div <- as.vector(A@value / B@value)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(log(div))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(1-div)
  expect_equal(as.matrix(grad(expr)[[B@id]]), val)
})

test_that("Test domain for MaxElemwise", {
  b <- Variable()
  expr <- MaxElemwise(a, b)
  a@value <- 2
  b@value <- 4
  expect_equal(grad(expr)[[a@id]], 0)
  expect_equal(grad(expr)[[b@id]], 1)
  
  a@value <- 3
  b@value <- 0
  expect_equal(grad(expr)[[a@id]], 1)
  expect_equal(grad(expr)[[b@id]], 0)
  
  a@value <- -1
  b@value <- 2
  expect_equal(grad(expr)[[a@id]], 0)
  expect_equal(grad(expr)[[b@id]], 1)
  
  y <- Variable(2)
  expr <- MaxElemwise(x, y)
  x@value <- c(3,4)
  y@value <- c(5,-5)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(0,1))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1,0))
  expect_equal(as.matrix(grad(expr)[[y@id]]), val)
  
  expr <- MaxElemwise(x, y)
  x@value <- c(-1e-9,4)
  y@value <- c(1,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(0,1))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1,0))
  expect_equal(as.matrix(grad(expr)[[y@id]]), val)
  
  expr <- MaxElemwise(A, B)
  A@value <- rbind(c(1,2), c(3,4))
  B@value <- rbind(c(5,1), c(3,2.3))
  div <- as.vector(A@value / B@value)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(0,1,1,1))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(1,0,0,0))
  expect_equal(as.matrix(grad(expr)[[B@id]]), val)
})

test_that("Test domain for MinElemwise", {
  b <- Variable()
  expr <- MinElemwise(a, b)
  a@value <- 2
  b@value <- 4
  expect_equal(grad(expr)[[a@id]], 1)
  expect_equal(grad(expr)[[b@id]], 0)
  
  a@value <- 3
  b@value <- 0
  expect_equal(grad(expr)[[a@id]], 0)
  expect_equal(grad(expr)[[b@id]], 1)
  
  a@value <- -1
  b@value <- 2
  expect_equal(grad(expr)[[a@id]], 1)
  expect_equal(grad(expr)[[b@id]], 0)
  
  y <- Variable(2)
  expr <- MinElemwise(x, y)
  x@value <- c(3,4)
  y@value <- c(5,-5)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1,0))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(0,1))
  expect_equal(as.matrix(grad(expr)[[y@id]]), val)
  
  expr <- MinElemwise(x, y)
  x@value <- c(-1e-9,4)
  y@value <- c(1,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1,1))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(0,0))
  expect_equal(as.matrix(grad(expr)[[y@id]]), val)
  
  expr <- MinElemwise(A, B)
  A@value <- rbind(c(1,2), c(3,4))
  B@value <- rbind(c(5,1), c(3,2.3))
  div <- as.vector(A@value / B@value)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(1,0,1,0))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(0,1,0,1))
  expect_equal(as.matrix(grad(expr)[[B@id]]), val)
})

test_that("Test domain for Power", {
  expr <- Sqrt(a)
  a@value <- 2
  expect_equal(grad(expr)[[a@id]], 0.5/sqrt(2))
  
  a@value <- 3
  expect_equal(grad(expr)[[a@id]], 0.5/sqrt(3))
  
  a@value <- -1
  expect_equal(grad(expr)[[a@id]], NA)
  
  expr <- x^3
  x@value <- c(3,4)
  expect_equal(as.matrix(grad(expr)[[x@id]]), rbind(c(27,0), c(0,48)))
  
  expr <- x^3
  x@value <- c(-1e-9,4)
  expect_equal(as.matrix(grad(expr)[[x@id]]), rbind(c(0,0), c(0,48)))
  
  expr <- A^2
  A@value <- rbind(c(1,-2), c(3,4))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(2,-4,6,8))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
  
  # Constant
  expr <- a^0
  expect_equal(grad(expr)[[a@id]], 0)
  
  expr <- x^0
  expect_equal(as.matrix(grad(expr)[[x@id]]), matrix(0, nrow = 2, ncol = 2))
})

test_that("Test grad for partial minimization/maximization problems", {
  for(obj in list(Minimize(a^-1), Maximize(Entr(a)))) {
    prob <- Problem(obj, list(x + a >= c(5,8)))
    
    # Optimize over nothing
    expr <- partial_optimize(prob, dont_opt_vars = list(x, a))
    a@value <- NA
    x@value <- NA
    grad <- grad(expr)
    expect_equal(grad[[a@id]], NA)
    expect_equal(grad[[x@id]], NA)
    
    # Outside domain
    a@value <- 1.0
    x@value <- c(5,5)
    grad <- grad(expr)
    expect_equal(grad[[a@id]], NA)
    expect_equal(grad[[x@id]], NA)
    
    a@value <- 1
    x@value <- c(10,10)
    grad <- grad(expr)
    expect_equal(grad[[a@id]], grad(obj@.args[[1]])[[a@id]])
    expect_equal(as.matrix(grad[[x@id]]), rep(0,4))
    
    # Optimize over x
    expr <- partial_optimize(prob, opt_vars = list(x))
    a@value <- 1
    grad <- grad(expr)
    expect_equal(grad[[a@id]], grad(obj@.args[[1]])[[a@id]] + 0)
    
    # Optimize over a
    fix_prob <- Problem(obj, list(x + a >= c(5,8), x == 0))
    result <- solve(fix_prob)
    dual_val <- fix_prob@constraints[1]@dual_variable@value
    expr <- partial_optimize(prob, opt_vars = list(a))
    x@value <- c(0,0)
    grad <- grad(expr)
    expect_equal(as.matrix(grad[[x@id]]), dual_val)
    
    # Optimize over x and a
    expr <- partial_optimize(prob, opt_vars = list(x, a))
    grad <- grad(expr)
    expect_equal(grad, list())
  }
})

test_that("Test grad for affine atoms", {
  expr <- -a
  a@value <- 2
  expect_equal(grad(expr)[[a@id]], -1)
  
  expr <- 2*a
  a@value <- 2
  expect_equal(grad(expr)[[a@id]], 2)
  
  expr <- a/2
  a@value <- 2
  expect_equal(grad(expr)[[a@id]], 0.5)
  
  expr <- -x
  x@value <- c(3,4)
  val <- matrix(0, nrow = 2, ncol = 2) - diag(c(1,1))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- -A
  A@value <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 4, ncol = 4) - diag(rep(1,4))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
  
  expr <- A[1,2]
  A@value <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 4, ncol = 1)
  val[3] <- 1
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
  
  z <- Variable(3)
  expr <- VStack(x, z)
  x@value <- c(1,2)
  z@value <- c(1,2,3)
  val <- matrix(0, nrow = 2, ncol = 5)
  val[,1:3] <- diag(rep(1,2))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  val <- matrix(0, nrow = 3, ncol = 5)
  val[,3:ncol(val)] <- diag(rep(1,3))
  expect_equal(as.matrix(grad(expr)[[z@id]]), val)
})
