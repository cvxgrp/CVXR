require(Matrix, quietly=TRUE)
require(ECOSolveR, quietly=TRUE)

N <- 5                 # 2N + 1 variables
L <- 1; h <- 2/(N-1)
c <- c(rep(0,N), rep(1,N), 0)

A <- Matrix(0, nrow=5, ncol=2*N+1, sparse=TRUE)
A[1, 2*N+1] <- 1                # x[2*N+1] = 1
A[2, 1] <- 1; A[3, N] <- 1      # x[1] = 0; x[N] = 1
A[4, N+1] <- 1; A[5, 2*N] <- 1  # y[1] = 1; y[N] = 1

b = c(h, 0, 1, 1, 1)

G <- Matrix(0, nrow=3*(N-1), ncol=2*N+1, sparse=TRUE)

for (i in 1:(N-1)) {
    j <- 3*(i-1) + 1
    G[j, 2*N+1] <- -1
    G[j+1, i] <- -1; G[j+1, i+1] <- 1
    G[j+2, N+i] <- -1; G[j+2, N+i+1] <- 1
}

H <- rep(0, 3*(N-1))

quad <- as.integer(rep(3, N-1))

system.time(
  sole <- ECOS_csolve(c, G, H, dims=list(q=quad), A, b)
)
