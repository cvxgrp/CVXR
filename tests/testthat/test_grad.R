a <- Variable(name = "a")

x <- Variable(2, name = "x")
y <- Variable(2, name = "y")

A <- Variable(2, 2, name = "A")
B <- Variable(2, 2, name = "B")
C <- Variable(3, 2, name = "C")

test_that("Test gradient for AffineProd", {
  expr <- AffineProd(C, A)
  value(C) <- rbind(c(1,-2), c(3,4), c(-1,-3))
  value(A) <- rbind(c(3,2), c(-5,1))
  
  expect_equal(as.matrix(grad(expr)[[C@id]]), rbind(c(3,0,0,2,0,0), c(0,3,0,0,2,0), c(0,0,3,0,0,2),
                                                    c(-5,0,0,1,0,0), c(0,-5,0,0,1,0), c(0,0,-5,0,0,1)), tolereance = TOL)
  expect_equal(as.matrix(grad(expr)[[A@id]]), rbind(c(1,3,-1,0,0,0), c(-2,4,-3,0,0,0), c(0,0,0,1,3,-1), 
                                                   c(0,0,0,-2,4,-3)), tolerance = TOL)
})

test_that("Test gradient for Pnorm", {
  expr <- Pnorm(x, 1)
  value(x) <- c(-1,0)
  expect_equal(as.matrix(grad(expr)[[x@id]]), as.matrix(c(-1,0)), tolerance = TOL)
  
  value(x) <- c(0,10)
  expect_equal(as.matrix(grad(expr)[[x@id]]), as.matrix(c(0,1)), tolerance = TOL)
  
  expr <- Pnorm(x, 2)
  value(x) <- c(-3,4)
  expect_equal(as.matrix(grad(expr)[[x@id]]), as.matrix(c(-3.0/5),4.0/5), tolerance = TOL)
  
  expr <- Pnorm(x, 0.5)
  value(x) <- c(-1,2)
  expect_equal(grad(expr)[[x@id]], NA)
  
  expr <- Pnorm(x, 0.5)
  value(x) <- c(0,0)
  expect_equal(grad(expr)[[x@id]], NA)
  
  expr <- Pnorm(x, 2)
  value(x) <- c(0,0)
  expect_equal(as.matrix(grad(expr)[[x@id]]), as.matrix(c(0,0)))
  
  expr <- Pnorm(A, 2)
  value(A) <- rbind(c(2,-2), c(2,2))
  expect_equal(as.matrix(grad(expr)[[A@id]]), matrix(c(0.5,0.5,-0.5,0.5)))
  
  expr <- Pnorm(A, 2, axis = 0)
  value(A) <- rbind(c(3,-3), c(4,4))
  expect_equal(as.matrix(grad(expr)[[A@id]]), rbind(c(0.6,0), c(0.8,0), c(0,-0.6), c(0,0.8)))
  
  expr <- Pnorm(A, 2, axis = 1)
  value(A) <- rbind(c(3,-4), c(4,3))
  expect_equal(as.matrix(grad(expr)[[A@id]]), rbind(c(0.6,0), c(0,0.8), c(-0.8,0), c(0,0.6)))
  
  expr <- Pnorm(A, 0.5)
  value(A) <- rbind(c(3,-4), c(4,3))
  expect_equal(grad(expr)[[A@id]], NA)
})

test_that("Test gradient for LogSumExp", {
  expr <- LogSumExp(x)
  value(x) <- c(0,1)
  e <- exp(1)
  expect_equal(as.matrix(grad(expr)[[x@id]]), c(1.0/(1+e), e/(1+e)))
  
  expr <- LogSumExp(A)
  value(A) <- rbind(c(0,1), c(-1,0))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(1.0/(2+e+1.0/e), 1.0/e/(2+e+1.0/e), e/(2+e+1.0/e), 1.0/(2+e+1.0/e)))
  
  expr <- LogSumExp(A, axis = 0)
  value(A) <- rbind(c(0,1), c(-1,0))
  expect_equal(as.matrix(grad(expr)[[A@id]]), cbind(c(1.0/(1+1.0/e), 1.0/e/(1+1.0/e), 0, 0), c(0, 0, e/(1+e), 1.0/(1+e))))
})

test_that("Test gradient for GeoMean", {
  # expr <- GeoMean(x)
  # value(x) <- c(1,2)
  # expect_equal(as.matrix(grad(expr)[[x@id]]), c(sqrt(2)/2, 1.0/2/sqrt(2)))
  
  # value(x) <- c(0,2)
  # expect_equal(grad(expr)[[x@id]], NA)
  
  # expr <- GeoMean(x, c(1,0))
  # value(x) <- c(1,2)
  # expect_equal(as.matrix(grad(expr)[[x@id]]), c(1,0))
  
  # No exception for single weight
  # value(x) <- c(-1,2)
  # expect_equal(grad(expr)[[x@id]], NA)
})

test_that("Test gradient for LambdaMax", {
  expr <- LambdaMax(A)
  value(A) <- rbind(c(2,0), c(0,1))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(1,0,0,0))
  
  value(A) <- rbind(c(1,0), c(0,2))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(0,0,0,1))
  
  value(A) <- rbind(c(1,0), c(0,1))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(0,0,0,1))
})

test_that("Test gradient for MatrixFrac", {
  expr <- MatrixFrac(A, B)
  value(A) <- diag(rep(1, 2))
  value(B) <- diag(rep(1, 2))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(2,0,0,2))
  expect_equal(as.matrix(grad(expr)[[B@id]]), c(-1,0,0,-1))
  
  value(B) <- matrix(0, nrow = 2, ncol = 2)
  expect_equal(grad(expr)[[A@id]], NA)
  expect_equal(grad(expr)[[B@id]], NA)
  
  expr <- MatrixFrac(x, A)
  value(x) <- c(2,3)
  value(A) <- diag(rep(1,2))
  expect_equal(as.matrix(grad(expr)[[x@id]]), c(4,6))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(-4,-6,-6,-9))
})

test_that("Test gradient for NormNuc", {
  expr <- NormNuc(A)
  value(A) <- rbind(c(10,4), c(4,30))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(1,0,0,1))
})

test_that("Test gradient for LogDet", {
  expr <- LogDet(A)
  value(A) <- 2*diag(rep(1,2))
  expect_equal(as.matrix(grad(expr)[[A@id]]), 1.0/2*diag(rep(1,2)))
  
  mat <- rbind(c(1,2), c(3,5))
  value(A) <- t(mat) %*% mat
  val <- t(solve(value(A)))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
  
  value(A) <- matrix(0, nrow = 2, ncol = 2)
  expect_equal(grad(expr)[[A@id]], NA)
  
  value(A) <- -rbind(c(1,2), c(3,4))
  expect_equal(grad(expr)[[A@id]], NA)
})

test_that("Test gradient for QuadOverLin", {
  expr <- QuadOverLin(x, a)
  value(x) <- c(1,2)
  value(a) <- 2
  expect_equal(as.matrix(grad(expr)[[x@id]]), c(1,2))
  expect_equal(grad(expr)[[a@id]], -1.25)
  
  value(a) <- 0
  expect_equal(grad(expr)[[x@id]], NA)
  expect_equal(grad(expr)[[a@id]], NA)
  
  expr <- QuadOverLin(A, a)
  value(A) <- diag(rep(1,2))
  value(a) <- 2
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(1,0,0,1))
  expect_equal(grad(expr)[[a@id]], -0.5)
  
  expr <- QuadOverLin(x, a) + QuadOverLin(y, a)
  value(x) <- c(1,2)
  value(a) <- 2
  value(y) <- c(1,2)
  value(a) <- 2
  expect_equal(as.matrix(grad(expr)[[x@id]]), c(1,2))
  expect_equal(as.matrix(grad(expr)[[y@id]]), c(1,2))
  expect_equal(grad(expr)[[a@id]], -2.5)
})

test_that("Test gradient for MaxEntries", {
  expr <- MaxEntries(x)
  value(x) <- c(2,1)
  expect_equal(as.matrix(grad(expr)[[x@id]]), c(1,0))
  
  expr <- MaxEntries(A)
  value(A) <- rbind(c(1,2), c(4,3))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(0,1,0,0))
  
  expr <- MaxEntries(A, axis = 0)
  value(A) <- rbind(c(1,2), c(4,3))
  expect_equal(as.matrix(grad(expr)[[A@id]]), rbind(c(0,0), c(1,0), c(0,0), c(0,1)))
  
  expr <- MaxEntries(A, axis = 1)
  value(A) <- rbind(c(1,2), c(4,3))
  expect_equal(as.matrix(grad(expr)[[A@id]]), rbind(c(0,0), c(0,1), c(1,0), c(0,0)))
})

test_that("Test SigmaMax", {
  expr <- SigmaMax(A)
  value(A) <- rbind(c(1,0), c(0,2))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(0,0,0,1))
  
  value(A) <- rbind(c(1,0), c(0,1))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(1,0,0,0))
})

test_that("Test SumLargest", {
  expr <- SumLargest(A, 2)
  
  value(A) <- rbind(c(4,3), c(2,1))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(1,0,1,0))
  
  value(A) <- rbind(c(1,2), c(3,0.5))
  expect_equal(as.matrix(grad(expr)[[A@id]]), c(0,1,1,0))
})

test_that("Test Abs", {
  expr <- Abs(A)
  value(A) <- rbind(c(1,2), c(-1,0))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(1,1,-1,0))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
})

test_that("Test linearize method", {
  # Affine
  expr <- (2*x - 5)[1]
  value(x) <- c(1,2)
  lin_expr <- linearize(expr)
  value(x) <- c(55,22)
  expect_equal(value(lin_expr), value(expr))
  value(x) <- c(-1,-5)
  expect_equal(value(lin_expr), value(expr))
  
  # Convex
  expr <- A^2 + 5
  expect_error(linearize(expr))
  
  value(A) <- rbind(c(1,2), c(3,4))
  lin_expr <- linearize(expr)
  manual <- value(expr) + 2*Reshape(value(Diag(Vec(A))) * Vec(A - value(A)), 2, 2)
  expect_equal(value(lin_expr), value(expr))
  value(A) <- rbind(c(-5,-5), c(8.2,4.4))
  expect_true(all(value(lin_expr) <= value(expr)))
  expect_equal(value(lin_expr), value(manual))
  
  # Concave
  expr <- Log(x)/2
  value(x) <- c(1,2)
  lin_expr <- linearize(expr)
  manual <- value(expr) + value(Diag(0.5*x^-1))*(x - value(x))
  expect_equal(value(lin_expr), value(expr))
  value(x) <- c(3,4.4)
  expect_true(all(value(lin_expr) >= value(expr)))
  expect_equal(value(lin_expr), value(manual))
})

test_that("Test gradient for Log", {
  expr <- Log(a)
  value(a) <- 2
  expect_equal(grad(expr)[[a@id]], 1.0/2)
  
  value(a) <- 3
  expect_equal(grad(expr)[[a@id]], 1.0/3)
  
  value(a) <- -1
  expect_equal(grad(expr)[[a@id]], NA)
  
  expr <- Log(x)
  value(x) <- c(3,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1/3,1/4))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Log(x)
  value(x) <- c(1e-9,4)
  expect_equal(grad(expr)[[x@id]], NA)
  
  expr <- Log(A)
  value(A) <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(1, 1/2, 1/3, 1/4))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
})

test_that("Test domain for Log1p", {
  expr <- Log1p(a)
  value(a) <- 2
  expect_equal(grad(expr)[[a@id]], 1.0/3)
  
  value(a) <- 3
  expect_equal(grad(expr)[[a@id]], 1.0/4)
  
  value(a) <- -1
  expect_equal(grad(expr)[[a@id]], NA)
  
  expr <- Log1p(x)
  value(x) <- c(3,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1/4,1/5))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Log1p(x)
  value(x) <- c(-1e-9-1,4)
  expect_equal(grad(expr)[[x@id]], NA)
  
  expr <- Log1p(A)
  value(A) <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(1/2, 1/3, 1/4, 1/5))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
})

test_that("Test domain for Entr", {
  expr <- Entr(a)
  value(a) <- 2
  expect_equal(grad(expr)[[a@id]], -log(2)-1)
  
  value(a) <- 3
  expect_equal(grad(expr)[[a@id]], -(log(3)+1))
  
  value(a) <- -1
  expect_equal(grad(expr)[[a@id]], NA)
  
  expr <- Entr(x)
  value(x) <- c(3,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(-(log(c(3,4)) + 1))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Entr(x)
  value(x) <- c(-1e-9,4)
  expect_equal(grad(expr)[[x@id]], NA)
  
  expr <- Entr(A)
  value(A) <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(-(log(1:4)+1))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
})

test_that("Test domain for Exp", {
  expr <- Exp(a)
  value(a) <- 2
  expect_equal(grad(expr)[[a@id]], exp(2))
  
  value(a) <- 3
  expect_equal(grad(expr)[[a@id]], exp(3))
  
  value(a) <- -1
  expect_equal(grad(expr)[[a@id]], exp(-1))
  
  expr <- Exp(x)
  value(x) <- c(3,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(exp(c(3,4)))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Exp(x)
  value(x) <- c(-1e-9,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(exp(c(-1e-9,4)))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Exp(A)
  value(A) <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(exp(1:4))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
})

test_that("Test domain for logistic", {
  expr <- Logistic(a)
  value(a) <- 2
  expect_equal(grad(expr)[[a@id]], exp(2)/(1+exp(2)))
  
  value(a) <- 3
  expect_equal(grad(expr)[[a@id]], exp(3)/(1+exp(3)))
  
  value(a) <- -1
  expect_equal(grad(expr)[[a@id]], exp(-1)/(1+exp(-1)))
  
  expr <- Logistic(x)
  value(x) <- c(3,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(exp(c(3,4))/(1+exp(3,4)))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Logistic(x)
  value(x) <- c(-1e-9,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(exp(c(-1e-9,4))/(1+exp(c(-1e-9,4))))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Logistic(A)
  value(A) <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(exp(1:4)/(1+exp(1:4)))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
})

test_that("Test domain for Huber", {
  expr <- Huber(a)
  value(a) <- 2
  expect_equal(grad(expr)[[a@id]], 2)
  
  expr <- Huber(a, M = 2)
  value(a) <- 3
  expect_equal(grad(expr)[[a@id]], 4)
  
  value(a) <- -1
  expect_equal(grad(expr)[[a@id]], -2)
  
  expr <- Huber(x)
  value(x) <- c(3,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(2,2))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Huber(x)
  value(x) <- c(-1e-9,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(0,2))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- Huber(A, M = 3)
  value(A) <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(2,4,6,6))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
})

test_that("Test domain for KLDiv", {
  b <- Variable()
  expr <- KLDiv(a, b)
  value(a) <- 2
  value(b) <- 4
  expect_equal(grad(expr)[[a@id]], log(2/4))
  expect_equal(grad(expr)[[b@id]], 1-(2/4))
  
  value(a) <- 3
  value(b) <- 0
  expect_equal(grad(expr)[[a@id]], NA)
  expect_equal(grad(expr)[[b@id]], NA)
  
  value(a) <- -1
  value(b) <- 2
  expect_equal(grad(expr)[[a@id]], NA)
  expect_equal(grad(expr)[[b@id]], NA)
  
  y <- Variable(2)
  expr <- KLDiv(x, y)
  value(x) <- c(3,4)
  value(y) <- c(5,8)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(log(c(3,4)) - log(c(5,8)))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  val <-matrix(0, nrow = 2, ncol = 2) + diag(c(1-3/5,1-4/8))
  expect_equal(as.matrix(grad(expr)[[y@id]]), val)
  
  expr <- KLDiv(x, y)
  value(x) <- c(-1e-9,4)
  value(y) <- c(1,2)
  expect_equal(grad(expr)[[x@id]], NA)
  expect_equal(grad(expr)[[y@id]], NA)
  
  expr <- KLDiv(A, B)
  value(A) <- rbind(c(1,2), c(3,4))
  value(B) <- rbind(c(5,1), c(3.5,2.3))
  div <- as.vector(value(A) / value(B))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(log(div))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(1-div)
  expect_equal(as.matrix(grad(expr)[[B@id]]), val)
})

test_that("Test domain for MaxElemwise", {
  b <- Variable()
  expr <- MaxElemwise(a, b)
  value(a) <- 2
  value(b) <- 4
  expect_equal(grad(expr)[[a@id]], 0)
  expect_equal(grad(expr)[[b@id]], 1)
  
  value(a) <- 3
  value(b) <- 0
  expect_equal(grad(expr)[[a@id]], 1)
  expect_equal(grad(expr)[[b@id]], 0)
  
  value(a) <- -1
  value(b) <- 2
  expect_equal(grad(expr)[[a@id]], 0)
  expect_equal(grad(expr)[[b@id]], 1)
  
  y <- Variable(2)
  expr <- MaxElemwise(x, y)
  value(x) <- c(3,4)
  value(y) <- c(5,-5)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(0,1))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1,0))
  expect_equal(as.matrix(grad(expr)[[y@id]]), val)
  
  expr <- MaxElemwise(x, y)
  value(x) <- c(-1e-9,4)
  value(y) <- c(1,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(0,1))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1,0))
  expect_equal(as.matrix(grad(expr)[[y@id]]), val)
  
  expr <- MaxElemwise(A, B)
  value(A) <- rbind(c(1,2), c(3,4))
  value(B) <- rbind(c(5,1), c(3,2.3))
  div <- as.vector(value(A) / value(B))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(0,1,1,1))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(1,0,0,0))
  expect_equal(as.matrix(grad(expr)[[B@id]]), val)
})

test_that("Test domain for MinElemwise", {
  b <- Variable()
  expr <- MinElemwise(a, b)
  value(a) <- 2
  value(b) <- 4
  expect_equal(grad(expr)[[a@id]], 1)
  expect_equal(grad(expr)[[b@id]], 0)
  
  value(a) <- 3
  value(b) <- 0
  expect_equal(grad(expr)[[a@id]], 0)
  expect_equal(grad(expr)[[b@id]], 1)
  
  value(a) <- -1
  value(b) <- 2
  expect_equal(grad(expr)[[a@id]], 1)
  expect_equal(grad(expr)[[b@id]], 0)
  
  y <- Variable(2)
  expr <- MinElemwise(x, y)
  value(x) <- c(3,4)
  value(y) <- c(5,-5)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1,0))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(0,1))
  expect_equal(as.matrix(grad(expr)[[y@id]]), val)
  
  expr <- MinElemwise(x, y)
  value(x) <- c(-1e-9,4)
  value(y) <- c(1,4)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(1,1))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  val <- matrix(0, nrow = 2, ncol = 2) + diag(c(0,0))
  expect_equal(as.matrix(grad(expr)[[y@id]]), val)
  
  expr <- MinElemwise(A, B)
  value(A) <- rbind(c(1,2), c(3,4))
  value(B) <- rbind(c(5,1), c(3,2.3))
  div <- as.vector(value(A) / value(B))
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(1,0,1,0))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
  val <- matrix(0, nrow = 4, ncol = 4) + diag(c(0,1,0,1))
  expect_equal(as.matrix(grad(expr)[[B@id]]), val)
})

test_that("Test domain for Power", {
  expr <- Sqrt(a)
  value(a) <- 2
  expect_equal(grad(expr)[[a@id]], 0.5/sqrt(2))
  
  value(a) <- 3
  expect_equal(grad(expr)[[a@id]], 0.5/sqrt(3))
  
  value(a) <- -1
  expect_equal(grad(expr)[[a@id]], NA)
  
  expr <- x^3
  value(x) <- c(3,4)
  expect_equal(as.matrix(grad(expr)[[x@id]]), rbind(c(27,0), c(0,48)))
  
  expr <- x^3
  value(x) <- c(-1e-9,4)
  expect_equal(as.matrix(grad(expr)[[x@id]]), rbind(c(0,0), c(0,48)))
  
  expr <- A^2
  value(A) <- rbind(c(1,-2), c(3,4))
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
    value(a) <- NA
    value(x) <- NA
    grad <- grad(expr)
    expect_equal(grad[[a@id]], NA)
    expect_equal(grad[[x@id]], NA)
    
    # Outside domain
    value(a) <- 1.0
    value(x) <- c(5,5)
    grad <- grad(expr)
    expect_equal(grad[[a@id]], NA)
    expect_equal(grad[[x@id]], NA)
    
    value(a) <- 1
    value(x) <- c(10,10)
    grad <- grad(expr)
    expect_equal(grad[[a@id]], grad(obj@.args[[1]])[[a@id]])
    expect_equal(as.matrix(grad[[x@id]]), rep(0,4))
    
    # Optimize over x
    expr <- partial_optimize(prob, opt_vars = list(x))
    value(a) <- 1
    grad <- grad(expr)
    expect_equal(grad[[a@id]], grad(obj@.args[[1]])[[a@id]] + 0)
    
    # Optimize over a
    fix_prob <- Problem(obj, list(x + a >= c(5,8), x == 0))
    result <- solve(fix_prob)
    dual_val <- fix_prob@constraints[1]@dual_variable@value
    expr <- partial_optimize(prob, opt_vars = list(a))
    value(x) <- c(0,0)
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
  value(a) <- 2
  expect_equal(grad(expr)[[a@id]], -1)
  
  expr <- 2*a
  value(a) <- 2
  expect_equal(grad(expr)[[a@id]], 2)
  
  expr <- a/2
  value(a) <- 2
  expect_equal(grad(expr)[[a@id]], 0.5)
  
  expr <- -x
  value(x) <- c(3,4)
  val <- matrix(0, nrow = 2, ncol = 2) - diag(c(1,1))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  expr <- -A
  value(A) <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 4, ncol = 4) - diag(rep(1,4))
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
  
  expr <- A[1,2]
  value(A) <- rbind(c(1,2), c(3,4))
  val <- matrix(0, nrow = 4, ncol = 1)
  val[3] <- 1
  expect_equal(as.matrix(grad(expr)[[A@id]]), val)
  
  z <- Variable(3)
  expr <- VStack(x, z)
  value(x) <- c(1,2)
  z@value <- c(1,2,3)
  val <- matrix(0, nrow = 2, ncol = 5)
  val[,1:3] <- diag(rep(1,2))
  expect_equal(as.matrix(grad(expr)[[x@id]]), val)
  
  val <- matrix(0, nrow = 3, ncol = 5)
  val[,3:ncol(val)] <- diag(rep(1,3))
  expect_equal(as.matrix(grad(expr)[[z@id]]), val)
})
