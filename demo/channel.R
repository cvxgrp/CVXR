# Problem data
n <- 2
m <- 2
P <- rbind(c(0.75, 0.25),    # Channel transition matrix
           c(0.25, 0.75))

# Form problem
x <- Variable(n)   # Probability distribution of input signal x(t)
y <- P %*% x       # Probability distribution of output signal y(t)
c <- apply(P * log2(P), 2, sum)
I <- c %*% x + sum(entr(y))   # Mutual information between x and y
obj <- Maximize(I)
constraints <- list(sum(x) == 1, x >= 0)
prob <- Problem(obj, constraints)

result <- solve(prob)
result$value
result$getValue(x)
