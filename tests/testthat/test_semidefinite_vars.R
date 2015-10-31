# X <- Semidef(2)
Y <- Variable(2, 2)
Fmat <- cbind(c(1,0), c(0,-1))

test_that("SDP in objective and constraint", {
  # SDP in objective
  # obj <- Minimize(SumEntries(Square(X - Fmat)))
  # p <- Problem(obj, list())
  
  # SDP in constraint
  # obj <- Minimize(SumEntries(Square(Y - Fmat)))
  # p <- Problem(obj, list(Y == Semidef(2)))
  
  # Index into semidef
  # obj <- Minimize(Square(X[1,1] - 1) + Square(X[2,1] - 2) + Square(X[2,2] - 4))
})